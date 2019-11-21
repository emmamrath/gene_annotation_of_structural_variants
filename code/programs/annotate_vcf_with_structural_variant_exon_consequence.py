#!/usr/bin/python
# python annotate_vcf_with_structural_variant_exon_consequence.py -i infile -o outfile -g genes_canonical_transcript_exons -t tmpdir
# python annotate_vcf_with_structural_variant_exon_consequence.py -i test1.vcf -o test1_annotated.vcf -g UCSC_Genes_canonical_exons_20180613.vcf.gz

# Input vcf is 1-based. Start pos in POS. End pos in INFO.END. Annotation field is INFO.SV_CANONICAL_CONSEQUENCE

# This program calls bedtools to find the intersection between the structural variant and a reference file of exons.

# The input file of exon regions to this program should be exons plus 2 bp either side to include the splice sites.
# If the exon splice sites are impacted by DEL, INDEL, INS, INV, then it is probably deleterious.
# However, if the exon splice sites are impacted by a DUP, then it is not deleterious. 
# This program keeps or removes 2 bp either side of the exons where necessary for analysis.
# The input file needs to have the 2 bp added because prior to this program it is used to do an intersect with bedtools to identify SVs overlapping exons+splice-sites.

# Output rules:
# <DUP		==>	if hits exons at start or end of transcript, = probably_benign
#			else if hits exons, = probably_deleterious
# <DEL		==>	if hits exons, = probably_deleterious
# <INDEL	==>	same as <DEL
# <INS		==>	if encompasses exons but doesn't have breakends inside exon, = probably_benign
#			else if has breakend inside exon, = probably deleterious
# <INV		==>	same as <INS

# This program does not compare the replacement sequence with reference sequence, thus some variants reported as deleterious may not be.
# Eg. if an INV at a splice site results in the splice site being conserved, or the INDEL sequence code for the same amino acids.


__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2018, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import datetime
import math
import commands
import argparse
import re
import pybedtools
from datetime import datetime
from uuid import uuid4

######################################################
def is_integer(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

######################################################
def extract_variant_end( info ):

	variant_end = ''

        infields = info.split(';')
        for this_info_key_value in infields:
                bits = this_info_key_value.split('=')
                this_key = bits[0]
                if (this_key == 'END'):
                        variant_end = str(bits[1])

        return variant_end

######################################################
def add_consequence_to_info_fields( info, consequence, info_code ):

	return_info = info
	if ((info != '') and (consequence != '')):
		return_info = info_code + '=' + consequence + ';' + info
	elif (consequence != ''):
		return_info = info_code + '=' + consequence
	return return_info

######################################################
def choose_final_consequence( consequence1, consequence2 ):

	return_consequence = consequence1
	if ((consequence1 != '') and (consequence2 == '')):
		return_consequence = consequence1
	elif ((consequence1 == '') and (consequence2 != '')):
		return_consequence = consequence2
	else:
		if ((consequence1 == 'PROBABLY_DELETERIOUS_DUPLICATION') and (consequence2 == 'PROBABLY_BENIGN_DUPLICATION')):
			return_consequence = 'PROBABLY_BENIGN_DUPLICATION'
		elif ((consequence1 == 'PROBABLY_BENIGN_DUPLICATION') and (consequence2 == 'PROBABLY_DELETERIOUS_DUPLICATION')):
			return_consequence = 'PROBABLY_BENIGN_DUPLICATION'
		if ((consequence1 == 'PROBABLY_DELETERIOUS') or (consequence2 == 'PROBABLY_DELETERIOUS')):
			return_consequence = 'PROBABLY_DELETERIOUS'
		elif ((consequence1 == 'PROBABLY_BENIGN') or (consequence2 == 'PROBABLY_BENIGN')):
			return_consequence = 'PROBABLY_BENIGN'
	return return_consequence

######################################################
def choose_final_consequence_2( consequence ):

	return_consequence = consequence
	if (consequence == 'PROBABLY_DELETERIOUS_DUPLICATION'):
		return_consequence = 'PROBABLY_DELETERIOUS'
	elif (consequence == 'PROBABLY_BENIGN_DUPLICATION'):
		return_consequence = 'PROBABLY_BENIGN'
	return return_consequence

######################################################
def determine_one_consequence( alt, variant_start, variant_end, reference_vcf_line ):

	consequence = ''

	variant_start = int(variant_start)
	variant_end = int(variant_end)

	if (len(reference_vcf_line) > 0):
		infields = reference_vcf_line.split("\t")
		bedtools_ref_chrom = str(infields[0])
		bedtools_ref_start = int(infields[1]) # This is the start of the CDS minus 2 bp for splice-site. For CDS that starts in middle of exon, this extra 2 bp is meaningless.
		bedtools_ref_end = int(infields[2]) # This is the end of the CDS plus 2 bp for splice-site. For CDS that ends in middle of exon, this extra 2 bp is meaningless.
		bedtools_gene_symbol = str(infields[3])
		bedtools_strand = str(infields[4])
		bedtools_ref_this_exon_num = int(infields[5])
		bedtools_ref_num_exons = int(infields[6])
		bedtools_true_exon_start = int(infields[7])
		bedtools_true_exon_end = int(infields[8])

		# These co-ordinates do not include the 2 bp splice-sites
		bedtools_ref_cds_exon_start = bedtools_ref_start + 2
		bedtools_ref_cds_exon_end = bedtools_ref_end - 2

		# These co-ordinates do include the 2 bp splice-sites
		bedtools_ref_cds_exon_start_of_splice_site = bedtools_ref_start
		bedtools_ref_cds_exon_end_of_splice_site = bedtools_ref_end

		# These co-ordinates do include the 2 bp splice-sites, unless this CDS exon starts or ends in the middle of an exon, in which case it doesn't include "splice-sites" inside an exon
		bedtools_ref_cds_exon_start_that_matters = bedtools_ref_cds_exon_start_of_splice_site
		bedtools_ref_cds_exon_end_that_matters = bedtools_ref_cds_exon_end_of_splice_site
		if ((bedtools_true_exon_start - 2) < bedtools_ref_cds_exon_start_of_splice_site):
			bedtools_ref_cds_exon_start_that_matters = bedtools_ref_cds_exon_start
		if ((bedtools_true_exon_end + 2) > bedtools_ref_cds_exon_end_of_splice_site):
			bedtools_ref_cds_exon_end_that_matters = bedtools_ref_cds_exon_end

		svtype = ''
		if (alt == '<DEL>'):
			svtype = 'DEL'
		elif (alt == '<INDEL>'):
			svtype = 'INDEL'
		elif (alt == 'INS'):
			svtype = 'INS'
		else:
			if (len(alt) >= 4):
				if (alt[0:4] == '<DUP'):
					svtype = 'DUP'
				elif (alt[0:4] == '<INV'):
					svtype = 'INV'

		if ((svtype == 'DEL') or (svtype == 'INDEL') or (svtype == 'INV')):

			# If DEL/INDEL/INV overlaps coding exon and/or its 2 bp splice-sites, regardless of whether breakend is inside an exon, then it is probably deleterious.
			consequence = 'PROBABLY_DELETERIOUS'

			# The only case where this SV in the 2 bp splice-sites is not deleterious is when this CDS exon is starting or ending inside the actual exon
			# and thus the splice-sites are not adjacent to this start/end of the CDS exon.
			if ((bedtools_true_exon_start - 2) < bedtools_ref_cds_exon_start_of_splice_site):
				if (variant_end < bedtools_true_exon_start):
					consequence = ''
			if ((bedtools_true_exon_end + 2) > bedtools_ref_cds_exon_end_of_splice_site):
				if (variant_start > bedtools_true_exon_end):
					consequence = ''

		else: # DUP, INS

			# If breakend of DUP falls inside an exon, then deleterious.
			# If breakend of exon falls inside this DUP, then deleterious.
			# If breakend of INS fall inside an exon or splice-site, then deleterious.
			breakend_is_in_exon = False
			if (svtype == 'DUP'):
				if ((bedtools_ref_cds_exon_start <= variant_start) and (variant_start < bedtools_ref_cds_exon_end)):
					breakend_is_in_exon = True
				elif ((bedtools_ref_cds_exon_start <= variant_end) and (variant_end <= bedtools_ref_cds_exon_end)):
					breakend_is_in_exon = True
				elif ((variant_start <= bedtools_ref_cds_exon_start) and (bedtools_ref_cds_exon_start < variant_end)):
					breakend_is_in_exon = True
				elif ((variant_start <= bedtools_ref_cds_exon_end) and (bedtools_ref_cds_exon_end <= variant_end)):
					breakend_is_in_exon = True
			else: # is INS
				# INS has the REF nucleotide followed by the insertion, so the insertion can be on the last nucleotide of the exon and not be deleterious
				if ((bedtools_ref_cds_exon_start_that_matters <= variant_start) and (variant_start < bedtools_ref_cds_exon_end_that_matters)):
					breakend_is_in_exon = True
				elif ((bedtools_ref_cds_exon_start_that_matters <= variant_end) and (variant_end < bedtools_ref_cds_exon_end_that_matters)):
					breakend_is_in_exon = True

			if (breakend_is_in_exon):
				consequence = 'PROBABLY_DELETERIOUS'

			else: # The breakends fall outside of exons, and encompass exon(s)

				if (svtype == 'INS'): # the breakends of insertion fall outside of exons
					consequence = 'PROBABLY_BENIGN'

				elif (svtype == 'DUP'): # the breakends of duplication fall outside of exons

					if ((variant_start < bedtools_ref_cds_exon_start) and (bedtools_ref_cds_exon_end <= variant_end)):
						if ((bedtools_ref_this_exon_num == 1) and (variant_start < bedtools_ref_cds_exon_start)):
							consequence = 'PROBABLY_BENIGN_DUPLICATION'
						elif ((bedtools_ref_this_exon_num == bedtools_ref_num_exons) and (bedtools_ref_cds_exon_end <= variant_end)):
							consequence = 'PROBABLY_BENIGN_DUPLICATION'
					else:
						consequence = 'PROBABLY_DELETERIOUS_DUPLICATION'

	return consequence

######################################################
def parse_out_any_warnings( command_output, command ):

	outlines = []

	all_outlines = command_output.split("\n")
	for outline in all_outlines:
		is_a_warning_line = False
		if (len(outline) >= 7):
			if (outline[0:7] == 'Warning'):
				is_a_warning_line = True
		if (is_a_warning_line == False):
			outlines.append( outline )

	return outlines

######################################################
def read_and_determine_consequence( genefile, alt, chrom, variant_start, variant_end, tmpdir ):

	consequence = ''
	if ((variant_start != '') and (variant_end != '')):
		eventid = datetime.now().strftime('%Y%m-%d%H-%M%S-') + str(uuid4())
		tmpfile1_filename = tmpdir + '/' + 'TEMPFILE1_' + str(chrom) + '_' + str(variant_start) + '_' + str(variant_end) + '_' + str(eventid) + '.txt'
		tmpfile2_filename = tmpdir + '/' + 'TEMPFILE2_' + str(chrom) + '_' + str(variant_start) + '_' + str(variant_end) + '_' + str(eventid) + '.txt'
		bedtools_key = str(chrom) + '\t' + str(variant_start) + '\t' + str(variant_end)
		tmpfile1 = open( tmpfile1_filename, 'w' )
		tmpfile1.write( bedtools_key + '\n' )
		tmpfile1.close()
		a = pybedtools.BedTool( genefile )
		a.intersect( tmpfile1_filename ).saveas( tmpfile2_filename )
		tmpfile2 = open( tmpfile2_filename, 'r' )
		bedtools_command_output_lines = tmpfile2.readlines()
		there_is_bedtools_output = False
		if (len(bedtools_command_output_lines) > 0):
			if (len(bedtools_command_output_lines) > 1):
				there_is_bedtools_output = True
			else:
				if (bedtools_command_output_lines[0] != ''):
					there_is_bedtools_output = True
		if (there_is_bedtools_output):
			for i in range( 0, len(bedtools_command_output_lines) ):
				one_consequence = determine_one_consequence( alt, variant_start, variant_end, bedtools_command_output_lines[i] )
				consequence = choose_final_consequence( consequence, one_consequence )
		consequence = choose_final_consequence_2( consequence )
		os.remove( tmpfile1_filename )
		os.remove( tmpfile2_filename )
	return consequence

######################################################
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in a VCF file of variants, find any annotations for it in a reference VCF, output the annotated variants in the VCF file.')
	parser.add_argument('-i', action="store", dest="infile", required=True, help='Input VCF file')
	parser.add_argument('-o', action="store", dest="outfile", required=True, help='Output annotated VCF file')
	parser.add_argument('-g', action="store", dest="genefile", required=True, help='Reference VCF file of gene exons')
	parser.add_argument('-c', action="store", dest="info_code", required=False, help='VCF.INFO field name for annotating VCF.INFO field')
	parser.add_argument('-t', action="store", dest="tmpdir", required=False, help='Temporary directory to use for temporary files')
	args = parser.parse_args()

	info_code = 'SV_CANONICAL_CONSEQUENCE'
	if (args.info_code is not None):
		info_code = str(args.info_code)
	tmpdir = '.'
	if (args.tmpdir is not None):
		tmpdir = str(args.tmpdir)

	# output the VCF file headers

	genefile = args.genefile
	outfile = open(args.outfile, 'w')

	# read input file and write out each record

	infile = open(args.infile, 'r')
	in_header = True
	for inline in infile:

		inline = inline.strip()
		if (inline != ''):

			if (in_header == True):
				if (len(inline) >= 6):
					if (inline[0:6] == '#CHROM'): # "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
						# This is the last header line. Output a new header line before outputting this last header line.
						outfile = open(args.outfile, 'a')
						outfile.write('##INFO=<ID=' + info_code + ',Number=1,Type=String,Description="Consequence of this structural variant that this a canonical transcript exon">' + "\n")
						outfile.close()
				if (len(inline) >= 1):
					if (inline[0:1] == '#'):
						outfile = open(args.outfile, 'a')
						outfile.write(inline + "\n")
						outfile.close()
					else:
						in_header = False

			if (in_header == False):

				infields = inline.split("\t")
				chrom = str(infields[0])
				pos = str(infields[1])
				vcfid = str(infields[2])
				ref = str(infields[3])
				alt = str(infields[4])
				qual = str(infields[5])
				vcffilter = str(infields[6])
				info = str(infields[7])
				if ( len(infields) >= 9):
					vcfformat = str(infields[8])

				variant_end = extract_variant_end( info )
				consequence = read_and_determine_consequence( genefile, alt, chrom, pos, variant_end, tmpdir )
				if (consequence != ''):
					info = add_consequence_to_info_fields( info, consequence, info_code )

				outline = chrom + "\t" + pos + "\t" + vcfid + "\t" + ref + "\t" + alt + "\t" + qual + "\t" + vcffilter + "\t" + info
				if ( len(infields) >= 9):
					outline = outline + "\t" + vcfformat
				for i in range( 9, len(infields) ):
					outline = outline + "\t" + str(infields[i])
				outline = outline + "\n"

				outfile = open(args.outfile, 'a')
				outfile.write( outline )
				outfile.close()

if __name__=='__main__':
    main()

