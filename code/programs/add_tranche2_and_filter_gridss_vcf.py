#!/usr/bin/python
# python add_tranche2_and_filter_gridss_vcf.py [-max_BND_length_for_no_split_reads]
# cat 2770-FD02806534-170221-HGM3TALXX-8.gridss.vcf | python add_tranche2_and_filter_gridss_vcf.py -max_BND_length_for_no_split_reads 100 > isks2770_gridss_BND_filtered.vcf

# Input is the vcf of structural variant BND breakend records called by gridss.
# Filter out some of those low quality BND calls.
# Assign a TRANCHE2 value in INFO field. Values are HIGH / INTERMEDIATE / LOW.
# Calculate a BNDVAF (variant allele frequency) for this gridss-called BND record.
# The BNDVAF = (SR + RP + IC + AS) / (REFPAIR + SR + RP + IC + AS)

# This program assigns both TRANCHE and TRANCHE2 to BND records.
# TRANCHE is assigned simply to continue what was done previously, even though the TRANCHE value is not used for filtering by this program.
# TRANCHE2 is assigned and then used for filtering BND records.

# Rules for filtering:
# The version of gridss used does not yet perform multiple test correction or score recalibration and QUAL scores are vastly overestimated for all variants. 
# As a rule of thumb for variants that do not have any flags in the filter field,
# variants with QUAL >= 1000 and have assembles from both sides of the breakpoint (AS > 0 & RAS > 0) are considered of high quality, 
# variant with QUAL >= 500 but can only be assembled from one breakend (AS > 0 | RAS > 0) are considered of intermediate quality, 
# and variants with low QUAL score or lack any supporting assemblies are considered the be of low quality.
# If FILTER contains any of the following SINGLE_ASSEMBLY, ASSEMBLY_ONLY, NO_ASSEMBLY, LOW_QUAL, then don't automatically filter it out. It might be a true positive.
# For these BND calls, if SRQ > 1 (Quality score of split reads supporting breakpoint), then assume it is a true positive.
# If SRQ = 0, then see whether this BND would be expected to have any split-reads. If not, then assume it is a true positive.
# Only if the SINGLE_ASSEMBLY/ASSEMBLY_ONLY/NO_ASSEMBLY/LOW_QUAL BND call fails all that, then assume it is a false positive and filter it out.
# If (abs(POS - ALT) + length_of_insert_in_ALT_field) <= 100, 
# then we don't expect split reads (and thus SINGLE_ASSEMBLY/ASSEMBLY_ONLY/NO_ASSEMBLY/LOW_QUAL will be called a true positive even if SRQ=0).
# The value of 100 is chosen because read length is 100, so a small SV of upto 100 length can be contained within a read without there being any split-reads and thus without there being any high quality split reads. 
# If IHOMPOS=-300,300, then this is a false positive, so filter out this read. 
# Gridss looks at 300 bp in either direction, and if the read sequence is homologous to the reference sequence, albeit inexactly, 
# then this is probably not a breakend afterall, so filter it out.

# Things that are still not filtered very well:
# If there is stuttering in the reads and subsequently two false positive BND reads are called as the breakends of a small INDEL.
# Even though SR=0, these false positives BNDs will make it through the filtering, because the BND_length will be small enough to allow that SR=0 be present.
# There does not seem to be any INFO fields from grids that would allow us to identify that these BNDs are a result of stuttering and are not true BNDs.
# In some cases, the stuttering causes the two BND records to have different IHOMPOS values.
# Thus, when the BNDs are joined to make an INDEL call, they may be identified as a false positive at that stage, 
# due to the lack of concordance of IHOMPOS of the two breakends.

# Fields produced by gridss that are used by this program:
# ##INFO=<ID=AS,Number=1,Type=Integer,Description="Count of assemblies supporting breakpoint">
# ##INFO=<ID=IC,Number=1,Type=Integer,Description="Count of read indels supporting breakpoint">
# ##INFO=<ID=RAS,Number=1,Type=Integer,Description="Count of assemblies supporting breakpoint from remote breakend">
# ##INFO=<ID=REF,Number=1,Type=Integer,Description="Count of reads mapping across this breakend">
# ##INFO=<ID=REFPAIR,Number=1,Type=Integer,Description="Count of reference read pairs spanning this breakpoint supporting the reference allele">
# ##INFO=<ID=RP,Number=1,Type=Integer,Description="Count of read pairs supporting breakpoint">
# ##INFO=<ID=SR,Number=1,Type=Integer,Description="Count of split reads supporting breakpoint">
# According to the above definitions, the following should be true:
# REFPAIR + SR = REF
# SR - num_SR_that_dont_have_read_pair_supporting_this_breakpoint + REFPAIR + IC = RP
# However, sometimes we see REF + 1 = REFPAIR when SR = 0 and we would expect REF = REFPAIR

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2019, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import commands
import argparse
import re

######################################################
def is_integer(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

######################################################
def is_float(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

######################################################
def extract_info_fields( in_info ):

	out_info = {}
        if ((in_info != '') and (in_info != '.')):
		info_fields = in_info.split(";")
		for key_value in info_fields:
			bits = key_value.split("=")
			info_key = bits[0]
			info_value = ""
			if (len(bits) > 1):
				info_value = bits[1]
			out_info[ info_key ] = info_value
        return out_info

######################################################
def append_field( old_field, extra_field, field_delimiter ):

	new_field = old_field
	if ((old_field == "") or (old_field == ".")):
		new_field = extra_field
	else:
		new_field = old_field + field_delimiter + extra_field
	return new_field

######################################################
def determine_tranche_for_one_record( this_qual, this_filter, this_AS, this_RAS ):

	this_tranche = ''
	if (this_filter != '.'):
		this_tranche = ''
	else:
		this_tranche = 'LOW'
		if (this_qual >= float(500)):
			if ((this_AS > 0) or (this_RAS > 0)):
				this_tranche = 'INTERMEDIATE'
		if (this_qual >= float(1000)):
			if ((this_AS > 0) and (this_RAS > 0)):
				this_tranche = 'HIGH'

	return this_tranche

######################################################
def filter_vcf_record_and_calc_new_fields( max_BND_length_for_no_split_reads, this_chrom, this_pos, this_alts, this_qual, this_filter, this_info_fields ):

	keep_this_vcf = False
	new_tranche = ""
	new_tranche2 = ""
	new_BNDVAF = ""

	is_homologous_to_reference = False
	if ("IHOMPOS" in this_info_fields):
		ihompos_1 = 0
		ihompos_2 = 0
		min_gridss_ihompos = -300
		max_gridss_ihompos = 300
		this_ihompos = this_info_fields["IHOMPOS"]
		bits = this_ihompos.split(",")
		if (is_integer(bits[0])):
			ihompos_1 = int(bits[0])
		if (len(bits) > 1):
			if (is_integer(bits[1])):
				ihompos_2 = int(bits[1])
		if ((ihompos_1 <= min_gridss_ihompos) and (ihompos_2 >= max_gridss_ihompos)) :
			is_homologous_to_reference = True
			keep_this_vcf = False

	if (is_homologous_to_reference == False):

		float_qual = 0
		if (is_float(this_qual)):
			float_qual = float(this_qual)

		##INFO=<ID=AS,Number=1,Type=Integer,Description="Count of assemblies supporting breakpoint">
		##INFO=<ID=IC,Number=1,Type=Integer,Description="Count of read indels supporting breakpoint">
		##INFO=<ID=IHOMPOS,Number=2,Type=Integer,Description="Position of inexact homology">
		##INFO=<ID=RAS,Number=1,Type=Integer,Description="Count of assemblies supporting breakpoint from remote breakend">
		##INFO=<ID=REF,Number=1,Type=Integer,Description="Count of reads mapping across this breakend">
		##INFO=<ID=REFPAIR,Number=1,Type=Integer,Description="Count of reference read pairs spanning this breakpoint supporting the reference allele">
		##INFO=<ID=RP,Number=1,Type=Integer,Description="Count of read pairs supporting breakpoint">
		##INFO=<ID=SC,Number=1,Type=String,Description="CIGAR for displaying anchoring alignment of any contributing evidence and microhomologies.">
		##INFO=<ID=SR,Number=1,Type=Integer,Description="Count of split reads supporting breakpoint">
		##INFO=<ID=SRQ,Number=1,Type=Float,Description="Quality score of split reads supporting breakpoint">

		this_AS = 0
		this_IC = 0
		this_RAS = 0
		this_REFPAIR = 0
		this_RP = 0
		this_SR = 0
		this_SRQ = 0
		if ("AS" in this_info_fields):
			if (is_integer(this_info_fields["AS"])):
				this_AS = int(this_info_fields["AS"])
		if ("IC" in this_info_fields):
			if (is_integer(this_info_fields["IC"])):
				this_IC = int(this_info_fields["IC"])
		if ("RAS" in this_info_fields):
			if (is_integer(this_info_fields["RAS"])):
				this_RAS = int(this_info_fields["RAS"])
		if ("REFPAIR" in this_info_fields):
			if (is_integer(this_info_fields["REFPAIR"])):
				this_REFPAIR = int(this_info_fields["REFPAIR"])
		if ("RP" in this_info_fields):
			if (is_integer(this_info_fields["RP"])):
				this_RP = int(this_info_fields["RP"])
		if ("SR" in this_info_fields):
			if (is_integer(this_info_fields["SR"])):
				this_SR = int(this_info_fields["SR"])
		if ("SRQ" in this_info_fields):
			if (is_float(this_info_fields["SRQ"])):
				this_SRQ = float(this_info_fields["SRQ"])

		# For gridss VCF BND records having nothing in filter field, then keep record and assign TRANCHE2 according to AS and RAS

		if (this_filter == "."):
			keep_this_vcf = True

		# If FILTER contains any of the following SINGLE_ASSEMBLY, ASSEMBLY_ONLY, NO_ASSEMBLY, LOW_QUAL, then don't automatically filter it out. It might be a true positive.
		# For these BND calls, if SRQ > 1 (Quality score of split reads supporting breakpoint), then assume it is a true positive.
		# If SRQ = 0, then see whether this BND would be expected to have any split-reads. If not, then assume it is a true positive.
		# Only if the SINGLE_ASSEMBLY/ASSEMBLY_ONLY/NO_ASSEMBLY/LOW_QUAL BND call fails all that, then assume it is a false positive and filter it out.
		# If (abs(POS - ALT) + length_of_insert_in_ALT_field) <= 100, 
		# then we don't expect split reads (and thus SINGLE_ASSEMBLY/ASSEMBLY_ONLY/NO_ASSEMBLY/LOW_QUAL will be called a true positive even if SRQ=0).
		# The value of 100 is chosen because read length is 100, 
		# so a small SV of upto 100 length can be contained within a read without there being any split-reads and thus without there being any high quality split reads. 

	 	else: # this_filter != "."
			this_filter_fields = this_filter.split(";")
			filter_is_in_list_of_possibly_good_filter_values = False
			list_of_possibly_good_filter_values = ( "SINGLE_ASSEMBLY", "ASSEMBLY_ONLY", "NO_ASSEMBLY", "LOW_QUAL" )
			for this_filter_value in this_filter_fields:
				if (this_filter_value in list_of_possibly_good_filter_values):
					filter_is_in_list_of_possibly_good_filter_values = True

			if (filter_is_in_list_of_possibly_good_filter_values):

				if (this_SRQ >= 1):
					keep_this_vcf = True

				else: # this_SRQ == 0, how long is the BND pos, can the BND be contained within a read and that could be why there are no split reads?

					# 1       1002653 gridss0_10117o  T       [1:1004132[T
					# 1       1004137 gridss0_8159o   C       ]1:1004214]C    1088.27
					# 1       1004214 gridss0_8159h   G       G[1:1004137[    1088.27
					# 1       1023510 gridss0_1524o   C       CCCCTCCCCCCGGCCCGCCCCCCCCCCCGCCCC[1:1023705[
					# 1       1023705 gridss0_1524h   C       ]1:1023510]CCCTCCCCCCGGCCCCCCCCCCCCGCCCCCCCCCCCGCCCCC

					this_pos = int(this_pos)
					this_alt_fields = this_alts.split(",")
					for this_alt in this_alt_fields:
						this_alt_chrom = ""
						this_alt_pos = 0
						this_alt_insert_seq_length = 0
						this_alt = this_alt.replace( "]", "[" )
						if (this_alt[0:1] == "["):
							bit = this_alt.replace( "[", ":" )
							bits = bit.split(":")
							this_alt_chrom = bits[1]
							this_alt_pos = int(bits[2])
							this_alt_insert_seq_length = len( bits[3] )
						else: # last character of this_alt is [ or ]
							bit = this_alt.replace( "[", ":" )
							bits = bit.split(":")
							this_alt_insert_seq_length = len( bits[0] )
							this_alt_chrom = bits[1]
							this_alt_pos = int(bits[2])
						if (this_alt_chrom == this_chrom):
							this_BND_length = abs(this_pos - this_alt_pos) + this_alt_insert_seq_length - 1
							if (this_BND_length < max_BND_length_for_no_split_reads):
								keep_this_vcf = True
						# else the BND is split over different chromosomes/contigs, so it is larger than a read, so a high qual BND would generate split reads.

	if (keep_this_vcf):

		# In theory, every BND record can have a TRANCHE value assigned it to, even before filtering.
		# In this program, BND records are filters, so only filtered records will have a TRANCHE assigned to them.
		new_tranche = determine_tranche_for_one_record( float_qual, this_filter, this_AS, this_RAS )

		if ((float_qual >= 1000) and (this_AS > 0) and (this_RAS > 0)):
			new_tranche2 = "HIGH"
		elif ((float_qual >= 500) and ((this_AS > 0) or (this_RAS > 0))):
			new_tranche2 = "INTERMEDIATE"
		else:
			new_tranche2 = "LOW"

		new_BNDVAF = 0
		new_BNDVAF_denominator = this_REFPAIR + this_SR + this_RP + this_IC + this_AS
		if (new_BNDVAF_denominator > 0):
			new_BNDVAF = float(this_SR + this_RP + this_IC + this_AS) / float(new_BNDVAF_denominator)
		new_BNDVAF = round( new_BNDVAF, 2 )

	return keep_this_vcf, new_tranche, new_tranche2, new_BNDVAF

######################################################
def main():

        # what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description="Read gridss structural variant BND VCF file from stdin, filter and add info fields, write to stdout.")
	parser.add_argument('-max_BND_length_for_no_split_reads', action="store", dest="max_BND_length_for_no_split_reads", required=False, help='Max length of BND breakpoint that would not be expected to produce split reads because it can be contained in one read (eg. if read length is 150 then max_BND_length_for_no_split_reads is probably around 100. Default is 100')
	args = parser.parse_args()

	max_BND_length_for_no_split_reads = 100
	if (args.max_BND_length_for_no_split_reads is not None):
		max_BND_length_for_no_split_reads = int(args.max_BND_length_for_no_split_reads)

	# Read in the input VCF file from STDIN

	in_header = True

	for inline in sys.stdin:

		if (in_header == True):
			if (len(inline) >= 1):
				first_char = inline[0:1]
				if (first_char != '#'):
					in_header = False

			if (len(inline) >= 6):
				first_six_chars = inline[0:6]
				if (first_six_chars == '#CHROM'):

					# Before writing out the last VCF header record, which is #CHROM POS...
					# write out the VCF header records for the new INFO fields added by this program.

					sys.stdout.write( '##INFO=<ID=TRANCHE,Number=1,Type=String,Description="Quality category of GRIDSS structural variant calls determined using QUAL,AS,RAS. Values are LOW INTERMEDIATE HIGH">' + "\n" )
					sys.stdout.write( '##INFO=<ID=TRANCHE2,Number=1,Type=String,Description="Quality category of GRIDSS structural variant calls determined using FILTER,SRQ,AS,RAS. Values are LOW INTERMEDIATE HIGH">' + "\n" )
					sys.stdout.write( '##INFO=<ID=BNDVAF,Number=1,Type=Float,Description="VAF of this gridss-called BND calculated as (SR+RP+IC+AS)/(REF+SR+RP+IC+AS)">' + "\n" )
					sys.stdout.write( '##FORMAT=<ID=TRANCHE,Number=1,Type=String,Description="Quality category of GRIDSS structural variant calls determined using QUAL,AS,RAS. Values are LOW INTERMEDIATE HIGH">' + "\n" )
					sys.stdout.write( '##FORMAT=<ID=TRANCHE2,Number=1,Type=String,Description="Quality category of GRIDSS structural variant calls determined using FILTER,SRQ,AS,RAS. Values are LOW INTERMEDIATE HIGH">' + "\n" )
					sys.stdout.write( '##FORMAT=<ID=BNDVAF,Number=1,Type=Float,Description="VAF of this gridss-called BND calculated as (SR+RP+IC+AS)/(REFPAIR+SR+RP+IC+AS)">' + "\n" )

		if (in_header == True): # Write out all the existing VCF header records.
			sys.stdout.write( inline )

		if (in_header == False): # We are processing VCF data records. We are no longer in the header part of the file.

			# Extract filter and gridss info fields.
			
			inline = inline.strip()
			infields = inline.split("\t")
			this_chrom = str(infields[0])
			this_pos = str(infields[1])
			this_id = str(infields[2])
			this_ref = str(infields[3])
			this_alts = str(infields[4])
			this_qual = str(infields[5])
			this_filter = str(infields[6])
			this_info = str(infields[7])

			this_info_fields = extract_info_fields( this_info )

			keep_this_vcf, new_tranche, new_tranche2, new_BNDVAF = filter_vcf_record_and_calc_new_fields( max_BND_length_for_no_split_reads, this_chrom, this_pos, this_alts, this_qual, this_filter, this_info_fields )

			if (keep_this_vcf):

				new_info = ""
				extra_info = "TRANCHE=" + str(new_tranche) + ";TRANCHE2=" + str(new_tranche2) + ";BNDVAF=" + str(new_BNDVAF)
				new_info = append_field( this_info, extra_info, ";" )
				outline = this_chrom + "\t" + this_pos + "\t" + this_id + "\t" + this_ref + "\t" + this_alts + "\t" + this_qual + "\t" + this_filter + "\t" + new_info
				if (len(infields) > 8):
					this_format = str(infields[8])
					extra_format = "TRANCHE:TRANCHE2:BNDVAF"
					new_format = append_field( this_format, extra_format, ":" )
					outline = outline + "\t" + new_format
					extra_sample = str(new_tranche) + ":" + str(new_tranche2) + ":" + str(new_BNDVAF)
					for i in range( 9, len(infields) ):
						this_sample = str(infields[i])
						new_sample = append_field( this_sample, extra_sample, ":" )
						outline = outline + "\t" + new_sample

				outline = outline + "\n"
				sys.stdout.write( outline )

if __name__=='__main__':
    main()


