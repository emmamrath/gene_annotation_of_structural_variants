#!/usr/bin/python
# python merge_a_multi_sample_SV_VCF_for_similar_position_variants.py input_Mobster input_window_in_basepairs input_reference_genome_fastq
# cat AllSamples1953_gridss_INS_DEL_INDEL_DUP_INV_BND.vcf | python merge_a_multi_sample_SV_VCF_for_similar_position_variants.py 100 /nvme/emmrat/reference_genome_for_MGRB/hs37d5x/hs37d5x.fa > AllSamples1953_gridss_INS_DEL_INDEL_DUP_INV_BND_merged100bp.vcf

# This program takes 1 VCF file on input from stdin, and produces 1 VCF file on output to stdout.
# This program assume that all inputs are of the same type of ALT allele and thus subsequence VCF records can be merged if the POS are close enough.
# Thus, make sure that this program receives VCF records of the same ALT allele type.
# Examples are <INS>, <INV>, <INDEL>, <DEL>, <DUP..., <INV...
# The program merges VCF records whose positions are within input_window_in_basepairs of the first variant seen in the group of close variants.
# A chain of variants where the first and last variants are further aparts than input_window_in_basepairs will thus not be merged.
# The final POS of merged variants will be the average POS of the merged variants.
# This means that we have to get the reference sequence's nucleotide at that new POS position.

# Which field to take for QUAL?
# This program has the same behaviour as bcftools merge, that is,
# QUAL receives the highest QUAL value of the merged VCF records.
# Which fields to take for ID and FILTER?
# Let's take the fields from the first VCF record in the merge group (until we find a reason to do otherwise).
# Which fields to take for FORMAT?
# Let's take the field from the first VCF record, and verify that all subsequent records have the same FORMAT.
# If they don't, crash the program.

# Example input:
# 1       83839   gridss0_1647    G       <INDEL> 1189.67 .       AC=1;AF=1.00;AN=1;END=83863;INSSEQ=AAAGAAAGAAAGAAAGAAAGAAAG;SVLEN=0;SVTYPE=INDEL;TRANCHE=HIGH;set=variant1773   GT      ./.     ./.
# 1       668632  gridss0_12662   C       <DUP:TANDEM>    975.94  .       AC=1;AF=1.00;AN=1;END=850205;SVLEN=181573;SVTYPE=DUP;TRANCHE=INTERMEDIATE;set=variant1012       GT      ./.     ./.     ./.
# 1       714060  gridss0_11538   C       <DUP:INS>       508.18  .       AC=1;AF=1.00;AN=1;END=714101;INSSEQ=GGAAGGCG;SVLEN=49;SVTYPE=DUP;TRANCHE=INTERMEDIATE;set=variant1008   GT      ./.     ./.
# 1       721941  gridss0_9781,gridss0_10170      T       <DUP:TANDEM>    719.61  .       AC=2;AF=1.00;AN=2;END=773013;SVLEN=51072;SVTYPE=DUP;set=variant295-variant1244  GT      ./.     ./.     ./.
# 1       755602  gridss0_2281,gridss0_6083       T       <DEL>   552.49  .       AC=2;AF=1.00;AN=2;END=755636;SVLEN=-34;SVTYPE=DEL;TRANCHE=INTERMEDIATE;set=variant234-variant1623       GT      ./.
# 1       11759449        gridss1_11970   C       <INV-DEL>       1116.96 .       AC=1;AF=1.00;AN=1;END=11760277;INVBND1DEL=1;INVBND2DEL=3;SVTYPE=INV;TRANCHE=HIGH;set=variant566 GT:ASQ:ASRP:ASSR:BAQ:BQ
# 1       52348953        gridss5_9118,gridss5_10703,gridss5_8521,gridss5_11355,gridss5_10965     G       <INV-DEL-INS>   1247.78 .       AC=5;AF=1.00;AN=5;INVBND1BDR3INS=GAG;INVBND2BDR3INS=TCA;SVTYPE=
# 1       16407687        gridss1_14409   G       <INV-DEL-OVERLAP>       823.46  .       AC=1;AF=1.00;AN=1;END=16408899;INVBND1DEL=4;INVBND2OVERLAP=138;SVTYPE=INV;TRANCHE=INTERMEDIATE;set=variant1594
# 1       16947758        gridss1_14169   G       <INV-DEL-INS-OVERLAP>   820.32  .       AC=1;AF=1.00;AN=1;END=17278450;INVBND1BDR5INS=CGTCC;INVBND1DEL=1965;INVBND2BDR5INS=AGGAC;INVBND2OVERLAP=391;SVT
# 1       46816730        gridss4_2152    C       <INV-DEL-CIRCULAR>      1714.21 .       AC=1;AF=1.00;AN=1;END=46820664;INVBND2DEL=99;SVTYPE=INV;TRANCHE=HIGH;set=variant1915    GT:ASQ:ASRP:ASSR:BAQ:BQ
# 1       44821990        gridss4_1389    A       <INV-DEL-INS-OVERLAP-CIRCULAR>  915.55  .       AC=1;AF=1.00;AN=1;END=44823691;INVBND1BDR5INS=CTTTTATTTATTTGTTTATTTAAGTCACGGAGTCTTGCCATGTTGCCCAGGCTAGCC

# Which fields to take for INFO?
# Take a TRANCHE=HIGH VCF record's INFO field over a TRANCHE=INTERMEDIATE. 
# If there are no TRANCHE=HIGH VCF records, then take a TRANCHE=INTERMEDIATE record's INFO field over a TRANCHE=LOW.
# Of all the TRANCHE=HIGH (or highest TRANCHE value records present), set INFO.END to the highest value present, and set INFO.SVLEN to the largest value present.

import sys
import os
import commands
import re

DELIMITER_FOR_INFO_MEINFO_FIELD = ','

######################################################
def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False



######################################################
def extract_info_end( info ):

	info_end = ''
	info_fields = info.split(';')
	for key_pair in info_fields:
		this_key = str(key_pair[0])
		if (this_key == 'END'):
			info_end = str(key_pair[1])
        return info_end

######################################################
def decide_if_same_variant( window, chrom1, chrom2, pos1, pos2_array, end1, end2_array ):

	same_variant = False
	if (chrom1 == chrom2):
		pos_is_same = False
		end_is_same = False
		pos1 = int(pos1)
		pos2 = ''
		if (len(pos2_array) > 0):
			pos2 = int(pos2_array[0])
			if ( abs(pos1 - pos2) <= window ):
				pos_is_same = True
				if (pos_is_same):
					if ((end1 == '') or (end1 == '.')):
						end_is_same = True
					else:
						end1 = int(end1)
						end2 = ''
						if (len(end2_array) == 0):
							end_is_same = True
						else:
							if (len(end2_array) > 0):
								end2 = int(end2_array[0])
								if ( abs(end1 - end2) <= window ):
									end_is_same = True
		if (pos_is_same and end_is_same):
			same_variant = True
        return same_variant

######################################################
def choose_qual( qual1, qual2 ):

	qual1 = float(qual1)
	qual2 = float(qual2)
	return_qual = qual1
	if (qual2 > qual1):
		return_qual = qual2
        return return_qual

######################################################
def construct_best_info_field( info_array ):

	best_info_idx = 0
	best_info_tranche = ''
	best_info_tranche2 = ''
	best_info_end = 0
	all_info_end = []
	best_info_svlen = 0
	best_info_fields = []
	points_system = { '':0, 'LOW':1, 'INTERMEDIATE':2, 'HIGH':3 }
	for i in range( 0, len(info_array) ):
		# get fields of this info
		this_info_tranche = ''
		this_info_tranche2 = ''
		this_info_end = 0
		this_info_svlen = 0
		info_fields_string = info_array[i]
		info_fields = info_fields_string.split(';')
		this_info_fields = {}
		for info_entry_pair in info_fields:
			bits = info_entry_pair.split('=')
			info_type = bits[0]
			info_value = ''
			if (len(bits) > 1):
				info_value = bits[1]
			if (info_type == 'END'):
				this_info_end = int(info_value)
				all_info_end.append( this_info_end )
			elif (info_type == 'SVLEN'):
				this_info_svlen = int(info_value)
			elif (info_type == 'TRANCHE'):
				this_info_tranche = info_value
			elif (info_type == 'TRANCHE2'):
				this_info_tranche2 = info_value
			this_info_fields[info_type] = info_value
		# get the best fields of this info
		if (i == 0):
			best_info_tranche = this_info_tranche
			best_info_tranche2 = this_info_tranche2
			best_info_svlen = this_info_svlen
			best_info_fields = this_info_fields
		else:
			best_points = points_system[best_info_tranche2]
			this_points = points_system[this_info_tranche2]
			if (this_points > best_points):
				best_info_tranche = this_info_tranche
				best_info_tranche2 = this_info_tranche2
				best_info_svlen = this_info_svlen
				best_info_fields = this_info_fields
			elif (this_points == best_points):
				if (abs(this_info_svlen) > abs(best_info_svlen)):
					best_info_svlen = this_info_svlen
					best_info_fields['SVLEN'] = this_info_svlen
				if (points_system[this_info_tranche] > points_system[best_info_tranche]):
					best_info_fields['TRANCHE'] = this_info_tranche
				if (points_system[this_info_tranche2] > points_system[best_info_tranche2]):
					best_info_fields['TRANCHE2'] = this_info_tranche2
	if (len(all_info_end) > 0):
		best_info_end = int( float(sum(all_info_end)) / float(len(all_info_end)) )
		best_info_fields['END'] = int(best_info_end)

	# build new INFO record
	new_info_string = ''
	for info_key in best_info_fields:
		if (new_info_string == ''):
			new_info_string = str(info_key) + '=' + str(best_info_fields[info_key])
		else:
			new_info_string = new_info_string + ';' + str(info_key) + '=' + str(best_info_fields[info_key])
        return new_info_string

######################################################
def create_new_format_from_2_different_formats( format1, format2  ):

	new_format = ''
	got_GT = False
	format1_fields = format1.split(':')
	format2_fields = format2.split(':')
	new_format_fields = {}
	for this_field in format1_fields:
		new_format_fields[this_field] = True
	for this_field in format2_fields:
		if this_field not in new_format_fields:
			new_format_fields[this_field] = True
	for this_field in new_format_fields:
		if (this_field == 'GT'):
			got_GT = True
	if (got_GT):
		new_format = 'GT'
	for this_field in new_format_fields:
		if ((this_field == 'GT') and (got_GT)):
			do_nothing = True # We already added this GT genotype field as the first field in the new format
		else:
			if (new_format == ''):
				new_format = str(this_field)
			else:
				new_format = new_format + ':' + str(this_field)
        return new_format

######################################################
def convert_sample_to_new_format( one_sample, old_format, new_format ):

	# It would be nice to assume that all VCF records have the same FORMAT fields in the sample fields
	# because when merging 2 VCF records, they will get the same FORMAT field in the one record, and all samples need to have all those fields.
	# However, it appears that GATK CombineVariants removes some of the FORMAT fields, and its output is input to this program,
	# so we can't assume that all VCF records to be merged have the same FORMAT fields.

	one_sample_reformatted = ''
	if (one_sample == './.'):
		one_sample_reformatted = './.'
	else:
		one_sample_fields = one_sample.split(':')
		old_format_fields = old_format.split(':')
		new_format_fields = new_format.split(':')
		one_sample_reformatted_fields = {}
		for this_field_name in new_format_fields:
			one_sample_reformatted_fields[this_field_name] = '.'
		for i in range( 0, len(old_format_fields) ):
			this_field_name = old_format_fields[i]
			this_field_value = '.'
			if (len(one_sample_fields) > i):
				this_field_value = one_sample_fields[i]
			one_sample_reformatted_fields[this_field_name] = this_field_value
		for this_field_name in new_format_fields:
			this_field_value = one_sample_reformatted_fields[this_field_name]
			if (one_sample_reformatted == ''):
				one_sample_reformatted = str(this_field_value)
			else:
				one_sample_reformatted = one_sample_reformatted + ':' + str(this_field_value)

        return one_sample_reformatted

######################################################
def find_first_separator_position( this_field ):

	separator_position = -1
	separator = ''
	slash_position = this_field.find('/')
	bar_position = this_field.find('|')
	if (slash_position >= 0):
		separator_position = slash_position
		separator = '/'
	if (bar_position >= 0):
		if (bar_position < separator_position):
			separator_position = bar_position
			separator = '|'

	return separator_position, separator

######################################################
def main():

	input_window_in_basepairs = sys.argv[1]
	input_window_in_basepairs = int(input_window_in_basepairs)
	input_reference_genome_fastq_file = sys.argv[2] # needs to have been indexed with bwa index

	# Read in the input VCF file from STDIN

	in_header = True
	previous_SV_inline_chrom = ''
	previous_SV_inline_pos_array = []
	previous_SV_inline_info_end_array = []
	previous_SV_inline_id = ''
	previous_SV_inline_ref = ''
	previous_SV_inline_alt = ''
	previous_SV_inline_qual = ''
	previous_SV_inline_filter = ''
	previous_SV_inline_info = ''
	previous_SV_info_array = []
	previous_SV_inline_format = ''
	previous_SV_inline_multiple_samples_fields = ''
	previous_SV_inline_subsequent_non_SV_inlines_to_write_out = []
	
	for inline in sys.stdin:

		if (in_header == True):
			if (len(inline) >= 1):
				first_char = inline[0:1]
				if (first_char != '#'):
					in_header = False

		if (in_header == True):
			sys.stdout.write( inline )

		else: # We are processing VCF data records. We are no longer in the header part of the file.

			# Is this variant of a type that we are interested in merging?

			inline_stripped = inline.strip()
			infields = inline_stripped.split("\t")
			this_inline_chrom = infields[0]
			this_inline_pos = infields[1]
			this_inline_id = infields[2]
			this_inline_ref = infields[3]
			this_inline_alt = infields[4]
			this_inline_qual = infields[5]
			this_inline_filter = infields[6]
			this_inline_info = infields[7]
			this_inline_info_end = extract_info_end( this_inline_info )
			this_inline_format = infields[8]
			this_inline_multiple_samples_fields = ''
			if (len(infields) >= 10):
				this_inline_multiple_samples_fields = infields[9]
				if (len(infields) >= 11):
					for i in range( 10, len(infields) ):
						this_inline_multiple_samples_fields = this_inline_multiple_samples_fields + "\t" + infields[i]

			variant_is_type_that_can_be_merged_result = True
			if (variant_is_type_that_can_be_merged_result == False):

				# If this variant is not of a type that we are interested in merging,
				# then it will be written out without being merged.
				# It may need to written out after merging of a previous variant prior to it
				# which may not yet be ready and completely merged to write out,
				# so that the output file is still sorted.

				if (previous_SV_inline_chrom == ''):
					sys.stdout.write( inline )
				else:
					previous_SV_inline_subsequent_non_SV_inlines_to_write_out.append( inline )

			else: # This variant is of a type that we will merge if it's POS is close to similar variants.

				similar_position_result = decide_if_same_variant( input_window_in_basepairs, this_inline_chrom, previous_SV_inline_chrom, this_inline_pos, previous_SV_inline_pos_array, this_inline_info_end, previous_SV_inline_info_end_array )

				if (similar_position_result == False):

					# If this variant is not close enough to be merged with previous variant(s)
					# then write out any previous variants, 
					# and write out any non-SV variants being saved to write out after the previous merged variants.
					# Don't write out this variant though,
					# because it might be close to future variants and be merged with them.

					# Write out any previous SV variants being merged

					if (previous_SV_inline_chrom != ''):

						# When writing out a merged variant, get average POS and it's corresponding REF

						previous_SV_inline_avg_pos = previous_SV_inline_pos_array[0]

						if (len(previous_SV_inline_pos_array) > 1):

							# The previous variant is a bunch of variants of similar position.
							# Now that we are writing them out as one merged variant
							# and know all the variants in this merged bunch,
							# get the average POS and corresponding REF

							# Get the average POS

							sum_of_pos = 0
							for this_pos in previous_SV_inline_pos_array:
								sum_of_pos = sum_of_pos + this_pos
							previous_SV_inline_avg_pos = int( float(sum_of_pos) / float(len(previous_SV_inline_pos_array)) )

							# Get the corresponding REF of that POS

							samtools_faidx_position = str(previous_SV_inline_chrom) + ':' + str(previous_SV_inline_avg_pos) + '-' + str(previous_SV_inline_avg_pos)
							samtools_faidx_command = 'samtools faidx ' + input_reference_genome_fastq_file + ' ' + samtools_faidx_position
							command_status, command_output = commands.getstatusoutput( samtools_faidx_command )
							if (command_status != 0):
								raise ValueError("\n\nWas not able to get REF position from reference genome. Thus will not continue processing any more mobile elements.\n")
							command_output_lines = command_output.split("\n")
							# some positions that mobile element mapped to may not be found in reference genome
							# eg. hs37d5:180728
							if (len(command_output_lines) > 1):
								previous_SV_inline_ref = command_output_lines[1].strip()

						previous_SV_inline_avg_info_end = ''
						if (len(previous_SV_inline_info_end_array) >= 1):
							previous_SV_inline_avg_info_end = previous_SV_inline_info_end_array[0]
							if (len(previous_SV_inline_info_end_array) > 1):
								# Get the average INFO.END
								sum_of_info_end = 0
								for this_info_end in previous_SV_inline_info_end_array:
									sum_of_info_end = sum_of_info_end + this_info_end
								previous_SV_inline_avg_info_end = sum_of_info_end / len(previous_SV_inline_info_end_array)
								previous_SV_inline_avg_info_end = int(previous_SV_inline_avg_info_end)

						best_info_field = construct_best_info_field( previous_SV_info_array )

						outline = str(previous_SV_inline_chrom) + "\t" + str(previous_SV_inline_avg_pos) + "\t" + str(previous_SV_inline_id) + "\t" + previous_SV_inline_ref + "\t" + previous_SV_inline_alt + "\t" + str(previous_SV_inline_qual) + "\t" + previous_SV_inline_filter + "\t" + best_info_field + "\t" + previous_SV_inline_format + "\t" + previous_SV_inline_multiple_samples_fields + "\n"
						sys.stdout.write( outline )

						previous_SV_info_array = []

					# Write out any previous non-SV variants being saved to write out after the previous merged variants.

					if ( len(previous_SV_inline_subsequent_non_SV_inlines_to_write_out) > 0 ):
						for previous_SV_inline_subsequent_inline in previous_SV_inline_subsequent_non_SV_inlines_to_write_out:
							sys.stdout.write( previous_SV_inline_subsequent_inline )
						previous_SV_inline_subsequent_non_SV_inlines_to_write_out = []

					# Don't write out this SV variant though.
					# Save it so it can be merged with future SV variants if need be.

					previous_SV_inline_chrom = this_inline_chrom
					previous_SV_inline_pos_array = []
					previous_SV_inline_pos_array.append( int(this_inline_pos) )
					previous_SV_inline_info_end_array = []
					if ((this_inline_info_end != '') and (this_inline_info_end != '.')):
						previous_SV_inline_info_end_array.append( int(this_inline_info_end) )
					previous_SV_inline_id = this_inline_id
					previous_SV_inline_ref = this_inline_ref
					previous_SV_inline_alt = this_inline_alt
					previous_SV_inline_qual = this_inline_qual
					previous_SV_inline_filter = this_inline_filter
					previous_SV_inline_info = this_inline_info
					previous_SV_info_array.append( this_inline_info )
					previous_SV_inline_format = this_inline_format
					previous_SV_inline_multiple_samples_fields = this_inline_multiple_samples_fields

				else: 

					# This ME variant needs to be merged with previous ME variant(s). 
					# There is at least one previous SV variant to merge this one with
					# because we have just seen that this SV variant is within the similarity window
					# of the first of the previous SV variants in the currently saved previous ME variants bunch.

					# Verify that variants to be merged have the same FORMAT field.
					# If they don't, then make them the same by expanding both to have all fields from both.

					if (previous_SV_inline_format != this_inline_format):
						new_inline_format = create_new_format_from_2_different_formats( this_inline_format, previous_SV_inline_format  )
						new_samples_fields = ''
						this_inline_multiple_samples_fields_split = this_inline_multiple_samples_fields.split("\t")
						i = 1
						for one_sample in this_inline_multiple_samples_fields_split:
							one_sample_reformatted = convert_sample_to_new_format( one_sample, this_inline_format, new_inline_format )
							if (new_samples_fields == ''):
								new_samples_fields = one_sample_reformatted
							else:
								new_samples_fields = new_samples_fields + "\t" + one_sample_reformatted
							i = i + 1
						this_inline_multiple_samples_fields = new_samples_fields
						new_samples_fields = ''
						previous_SV_inline_multiple_samples_fields_split = previous_SV_inline_multiple_samples_fields.split("\t")
						for one_sample in previous_SV_inline_multiple_samples_fields_split:
							one_sample_reformatted = convert_sample_to_new_format( one_sample, previous_SV_inline_format, new_inline_format )
							if (new_samples_fields == ''):
								new_samples_fields = one_sample_reformatted
							else:
								new_samples_fields = new_samples_fields + "\t" + one_sample_reformatted
						previous_SV_inline_multiple_samples_fields = new_samples_fields
						this_inline_format = new_inline_format
						previous_SV_inline_format = new_inline_format

					# Choose the highest QUAL from amongst the variants to be merged

					previous_SV_inline_qual = choose_qual( previous_SV_inline_qual, this_inline_qual )

					# Keep a list of all the POS being merged. Will choose the average when we have all the POS to be merged.

					previous_SV_inline_pos_array.append( int(this_inline_pos) )
					if ((this_inline_info_end != '') and (this_inline_info_end != '.')):
						previous_SV_inline_info_end_array.append( int(this_inline_inline_end) )

					# If the new VCF record's ALTs are different to the existing ALTs,
					# then add the new ALTs to the existing ALTS.
					# Then add each individual new sample's formatted fields 
					# as new sample-formatted-field columns
					# such that any Genotype GT field values point to the new ALT value and not to their old ALT value. 

					# Make the new ALT field as a string list of all the ALTs in the the previous ME variant
					# plus any new unique ALT values in this new VCF record.

					this_alt_fields = this_inline_alt.split(',')
					new_alt_fields_string = previous_SV_inline_alt 			# eg. <INS:ME:ALU> or <INS:ME:ALU>,<INS:ME:L1>
					for i in range( 0, len(this_alt_fields) ):
						this_alt_field = this_alt_fields[i]
						does_this_alt_field_exist_in_list = False
						temp_new_alt_fields = new_alt_fields_string.split(',')
						for j in range( 0, len(temp_new_alt_fields) ):
							temp_new_alt_field = temp_new_alt_fields[j]
							if (temp_new_alt_field == this_alt_field):
								does_this_alt_field_exist_in_list = True
						if (does_this_alt_field_exist_in_list == False):
							new_alt_fields_string = new_alt_fields_string + ',' + this_alt_field

					# Keep a list of all the INFO fields. When merging, will pick the best fields to be the INFO fields of merged record.

					previous_SV_info_array.append( this_inline_info )

					# Start processing ALT fields

					new_alt_fields = new_alt_fields_string.split(',')
					map_old_alt_to_new_alt = [0] * (len(this_alt_fields) + 1)
					i = 1
					for this_alt_field in this_alt_fields:
						for j in range( 0, len(new_alt_fields) ):
							if (this_alt_field == new_alt_fields[j]):
								new_i = j + 1
								map_old_alt_to_new_alt[i] = new_i
						i = i + 1

					# Now reformat all the formatted samples in this new VCF variant
					# to point to the ALT variants in the new ALT string.
					# Actually, Genotype GT is the only format field to be reformatted.
					# All other format types are just copied across, they don't point to ALT values.
					# Inside the sample's formatted GT field, each value has to be reformatted to point to correct new ALT value.
					# A new string for the formatted samples needs to be produced.
					# It will be the same as the string of formatted samples of the previous inline, except for the following.
					# Where this new VCF variant contains formatted sample for a sample, 
					# that sample's formatted field must come from this new VCF variant record,
					# not from the previous inline.
					# We assume that the GT field is the first field of FORMAT.
					# We will look at the GT fields of the formatted sample field of this new VCF variant record.
					# If it contains any values other than . then we assume that this new VCF variant record contains data for this sample.
					# Otherwise we assume that it does not, and will take the formatted sample field of the previous inline
					# that we are merging this new VCF variant record with.

					# Crash if the first FORMAT field is not GT

					this_inline_format_fields = this_inline_format.split(':')
					if (this_inline_format_fields[0] != 'GT'):
						raise ValueError("\n\nThe first FORMAT field for variant " + str(this_inline_chrom) + ":" + str(this_inline_pos) + " is not GT for Genotype. We are trying to merge multiple ALT variants into one and we need this genotype field to know which VCF record contains each sample's data. Will not continue processing this file.\n")

					# For each sample, figure out whether its data is in the new VCF variant record and thus needs the GT to be reformatted for the merged ALT field,
					# or is its data in the previous inline field so just copy the entire formatted data across from the previous inline field.

					where_is_sample_data_array = []
					previous_inline_sample_data_array = []
					this_inline_multiple_samples_fields_array = this_inline_multiple_samples_fields.split("\t")
					for sample_idx in range( 0, len(this_inline_multiple_samples_fields_array) ):
						previous_SV_inline_multiple_samples_fields_array = previous_SV_inline_multiple_samples_fields.split("\t")
						previous_inline_this_sample = previous_SV_inline_multiple_samples_fields_array[sample_idx]
						this_inline_this_sample = this_inline_multiple_samples_fields_array[sample_idx]
						previous_inline_this_sample_array = previous_inline_this_sample.split(':')
						this_inline_this_sample_array = this_inline_this_sample.split(':')
						previous_inline_this_sample_GT = previous_inline_this_sample_array[0]
						this_inline_this_sample_GT = this_inline_this_sample_array[0]
						this_inline_this_sample_GT_temp = this_inline_this_sample_GT
						this_inline_this_sample_GT_temp = this_inline_this_sample_GT_temp.replace('|','')
						this_inline_this_sample_GT_temp = this_inline_this_sample_GT_temp.replace('/','')
						this_inline_this_sample_GT_temp = this_inline_this_sample_GT_temp.replace('.','')
						where_is_sample_data = 'this_inline'
						if (this_inline_this_sample_GT_temp == ''):
							where_is_sample_data = 'previous_inline'
						where_is_sample_data_array.append( where_is_sample_data )
						previous_inline_sample_data_array.append( previous_inline_this_sample )

					# Now reconstruct the string of formatted sample fields.
					# The formatted data for each sample either comes from the previous_inline, 
					# or from this new VCF inline after being reformatted for the change in its merged ALT alleles

					new_this_inline_multiple_samples_fields_string = ''
					for sample_idx in range( 0, len(this_inline_multiple_samples_fields_array) ):

						formatted_sample_field_for_this_sample = ''
						if (where_is_sample_data_array[sample_idx] == 'previous_inline'):
							formatted_sample_field_for_this_sample = previous_inline_sample_data_array[sample_idx]

						else: # (where_is_sample_data_array[sample_idx] == 'this_inline')

							this_one_sample_multiple_format_fields = this_inline_multiple_samples_fields_array[sample_idx]
							new_one_sample_multiple_format_fields = ''
							this_inline_format_array = this_inline_format.split(':')
							this_inline_one_sample_multiple_format_fields_array = this_one_sample_multiple_format_fields.split(':')
							format_idx = 0
							for this_one_sample_one_format_field in this_inline_one_sample_multiple_format_fields_array:
								new_one_sample_one_format_field = ''
								this_format_type = this_inline_format_array[format_idx]
								if (this_format_type != 'GT'):
									new_one_sample_one_format_field = this_one_sample_one_format_field
								else: # This is the Genotype field, change values according to new ALT positions
									if (this_one_sample_one_format_field == '.'):
										new_one_sample_one_format_field = this_one_sample_one_format_field
									else:
										separator_position = find_first_separator_position( this_one_sample_one_format_field )
										if (separator_position == -1):
											raise ValueError("\n\nSample formatted genotype field for variant " + str(this_inline_chrom) + ":" + str(this_inline_pos) + " does not have a slash / or bar | and we are trying to merge multiple ALT variants into one and accordingly modify this genotype field. Will not continue processing this file.\n")
										sample_GT_fields = re.split( "\||\/", this_one_sample_one_format_field ) # split on | or /
										new_one_sample_one_format_field = ''
										remaining_GT_field_to_look_for_separators = this_one_sample_one_format_field
										for sample_GT_field in sample_GT_fields:
											new_sample_GT_field = '.'
											separator_position, sample_GT_separator = find_first_separator_position( remaining_GT_field_to_look_for_separators )
											remaining_GT_field_to_look_for_separators = remaining_GT_field_to_look_for_separators[ (separator_position+1) : ]
											if (sample_GT_field != '.'):
												new_sample_GT_field = map_old_alt_to_new_alt[int(sample_GT_field)]
												if (new_sample_GT_field == ''):
													raise ValueError("\n\nSample formatted genotype field for variant " + str(this_inline_chrom) + ":" + str(this_inline_pos) + " has a value that can't be matched to one of its ALT variants, and we are trying to modify it to match the new merged ALTs. Will not continue processing this file.\n")
											new_one_sample_one_format_field = new_one_sample_one_format_field + str(new_sample_GT_field) + sample_GT_separator

								# new_one_sample_multiple_format_fields = new_one_sample_multiple_format_fields + ':' + new_one_sample_one_format_field
								if (new_one_sample_one_format_field == ''):
									new_one_sample_one_format_field = '.'
								if (new_one_sample_multiple_format_fields == ''):
									new_one_sample_multiple_format_fields = new_one_sample_one_format_field
								else:
									new_one_sample_multiple_format_fields = new_one_sample_multiple_format_fields + ':' + new_one_sample_one_format_field
								format_idx = format_idx + 1
							formatted_sample_field_for_this_sample = new_one_sample_multiple_format_fields

						# new_this_inline_multiple_samples_fields_string = new_this_inline_multiple_samples_fields_string + "\t" + new_one_sample_multiple_format_fields
						if (new_this_inline_multiple_samples_fields_string == ''):
							new_this_inline_multiple_samples_fields_string = formatted_sample_field_for_this_sample
						else:
							new_this_inline_multiple_samples_fields_string = new_this_inline_multiple_samples_fields_string + "\t" + formatted_sample_field_for_this_sample

					if (new_this_inline_multiple_samples_fields_string == ''):
						new_this_inline_multiple_samples_fields_string = '.'

					previous_SV_inline_alt = new_alt_fields_string
					previous_SV_inline_multiple_samples_fields = new_this_inline_multiple_samples_fields_string

	# We have read and processed all the input VCF records.
	# If there is a previous SV merge group still hanging around then write it out.
	# If there are non-SV inlines still hanging around after them, then write them out.
	# Actually, there won't be any because we assume that the input is only SVs of a certain same mergable type.

	if (previous_SV_inline_chrom != ''):

		# Write out the previous ME merge group.
		# It has not been written out yet because we were still looking for more ME variants having similar POS
		# when we hit the end of the VCF file.
		# If necessary get the average position of the list of positions, and it's ref.seq. nucleotide

		previous_SV_inline_avg_pos = previous_SV_inline_pos_array[0]
		if (len(previous_SV_inline_pos_array) > 1):

			# The previous variant is a bunch of variants of similar position.
			# Now that we are writing them out as one merged variant
			# and know all the variants in this merged bunch,
			# get the average POS and corresponding REF

			# Get the average POS

			sum_of_pos = 0
			for this_pos in previous_SV_inline_pos_array:
				sum_of_pos = sum_of_pos + this_pos
			previous_SV_inline_avg_pos = int( float(sum_of_pos) / float(len(previous_SV_inline_pos_array)) )

			# Get the corresponding REF of that POS

			samtools_faidx_position = str(previous_SV_inline_chrom) + ':' + str(previous_SV_inline_avg_pos) + '-' + str(previous_SV_inline_avg_pos)
			samtools_faidx_command = 'samtools faidx ' + input_reference_genome_fastq_file + ' ' + samtools_faidx_position
			command_status, command_output = commands.getstatusoutput( samtools_faidx_command )
			if (command_status != 0):
				raise ValueError("\n\nWas not able to get REF position from reference genome. Thus will not continue processing any more mobile elements.\n")
			command_output_lines = command_output.split("\n")
			# some positions that mobile element mapped to may not be found in reference genome
			# eg. hs37d5:180728
			if (len(command_output_lines) > 1):
				previous_SV_inline_ref = command_output_lines[1].strip()

		previous_SV_inline_avg_info_end = ''
		if (len(previous_SV_inline_info_end_array) >= 1):
			previous_SV_inline_avg_info_end = previous_SV_inline_info_end_array[0]
			if (len(previous_SV_inline_info_end_array) > 1):
				# Get the average INFO.END
				sum_of_info_end = 0
				for this_info_end in previous_SV_inline_info_end_array:
					sum_of_info_end = sum_of_info_end + this_info_end
				previous_SV_inline_avg_info_end = sum_of_info_end / len(previous_SV_inline_info_end_array)
				previous_SV_inline_avg_info_end = int(previous_SV_inline_avg_info_end)

		best_info_field = construct_best_info_field( previous_SV_info_array )

		outline = str(previous_SV_inline_chrom) + "\t" + str(previous_SV_inline_avg_pos) + "\t" + str(previous_SV_inline_id) + "\t" + previous_SV_inline_ref + "\t" + previous_SV_inline_alt + "\t" + str(previous_SV_inline_qual) + "\t" + previous_SV_inline_filter + "\t" + best_info_field + "\t" + previous_SV_inline_format + "\t" + previous_SV_inline_multiple_samples_fields + "\n"
		sys.stdout.write( outline )

		# Write out any non-SV variants having positions after the previous SV merge group

		if ( len(previous_SV_inline_subsequent_non_SV_inlines_to_write_out) > 0 ):
			for previous_SV_inline_subsequent_non_SV_inline in previous_SV_inline_subsequent_non_SV_inlines_to_write_out:
				sys.stdout.write( previous_SV_inline_subsequent_non_SV_inline )
			previous_SV_inline_subsequent_non_SV_inlines_to_write_out = []

if __name__=='__main__':
    main()



