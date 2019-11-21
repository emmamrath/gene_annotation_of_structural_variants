#!/usr/bin/python
# cat in_vcf | python extract_vcf_info_annotation.py -info info_field > out_tsv
# cat MEI_183410_B_SV_INS_DEL_INDEL_DUP_INV_BND_merged100bp_annotated_with_genomic_regions_and_1000G.vcf | python extract_vcf_info_annotation.py -info HITEXONS | sort | uniq > MEI_183410_B_HITEXONS.tsv


__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2018, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import argparse
import commands
import re

######################################################
def is_integer(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

######################################################
def extract_info( arg_info, in_info ):

	out_info = ''
	if ((in_info != '') and (in_info != '.')):
		info_fields = in_info.split(";")
		i = 0
		keep_looking = True
		while (keep_looking):
			this_pair = info_fields[i]
			bits = this_pair.split("=")
			this_key = bits[0]
			if (this_key == arg_info):
				out_info = bits[1]
				keep_looking = False
			else:
				i = i + 1
				if (i >= len(info_fields)):
					keep_looking = False
	return out_info

######################################################
def main():

	comparison_window = 100

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description="Read in VCF file from STDIN, extract the info field value, output to STDOUT.")
	parser.add_argument('-info', action="store", dest="outinfo", required=False, help='Tab-delimited VCF.INFO field to be output.')
	args = parser.parse_args()

	arg_info = args.outinfo

	# compare input file

	for inline in sys.stdin:
		inline = inline.strip()
		if (len(inline) >= 1):
			if (inline[0:1] != "#"):
				infields = inline.split("\t")
				in_info = infields[7]
				out_info = extract_info( arg_info, in_info )
				if (out_info != ''):
					sys.stdout.write(out_info + "\n")

if __name__=='__main__':
    main()

