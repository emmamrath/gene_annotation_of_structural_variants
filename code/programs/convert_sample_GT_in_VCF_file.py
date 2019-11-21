#!/usr/bin/python
# python convert_sample_GT_in_VCF_file.py -i input_VCF -o output_VCF [-a all_flag]
# python convert_sample_GT_in_VCF_file.py -i test1.vcf -o test1_out.vcf -a ALL 

# If the sample field has values and yet the first format field is GT and is a dot,
# then change it to 1/. because dot signifies that this sample doesn't have this variant
# whereas the sample-specific information shows that the sample obviously does have this variant.

# If the -a option is present, assume all samples present have the variant,
# so change all . GT to 1/.

# If a sample does have the variant, then we don't want the GT field to be dot 
# because downstream sorting and VCF merging may decide that sample doesn't have the variant when GT field is a dot.

# GRIDSS outputs a GT of dot when reporting the variants of a sample.


__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2017, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import argparse
import commands


######################################################
def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

######################################################
def convert_to_integer(input_string):
	s_int = 0
	if ((input_string != 'NA') and (input_string != '')):
		s_int = int(input_string)
        return s_int

######################################################
def is_valid_ref(input_ref):
	is_valid = False
	valid_chars = 'ACGTNacgtn'
	if (valid_chars.find(input_ref) >= 0):
		is_valid = True
        return is_valid

######################################################
def main():

	parser = argparse.ArgumentParser(description='Read in VCF file, change GT of dot to GT 1/. and then output the VCF')
	parser.add_argument('-i', action="store", dest="in_vcf", required=True, help='Input VCF file')
	parser.add_argument('-o', action="store", dest="out_vcf", required=True, help='Output VCF file')
	parser.add_argument('-a', action="store", dest="all_flag", required=False, help='If this all-flag is present, change dot GT of all samples, otherwise only those that have other sample fields.')
	args = parser.parse_args()

	# Read in and process each VCF file record.
	# Then write out the record.

	out_vcf = open(args.out_vcf, 'w')
	in_vcf = open(args.in_vcf, 'r')
	for inline in in_vcf:
		inline = inline.strip()
		if (inline != ''):
			first_char = inline[0:1]
			if (first_char == '#'): # Just write out header lines, don't process them.
				outline = inline + "\n"
				out_vcf.write( outline )
			else:
				infields = inline.split("\t")
				this_format = infields[8]
				format_fields = this_format.split(':')
				if (format_fields[0] == 'GT'):

					if (len(infields) > 9):
						for i in range( 9, len(infields) ):
							this_sample = infields[i]
							sample_fields = this_sample.split(':')
							change_this_gt_field = False
							if (sample_fields[0] == '.'):
								if (args.all_flag is not None):
									change_this_gt_field = True
								else:
									if (len(sample_fields) > 1):
										change_this_gt_field = True
							if (change_this_gt_field):
								this_sample = '1/.'
								if (len(sample_fields) > 1):
									for j in range( 1, len(sample_fields) ):
										this_sample = this_sample + ':' + str(sample_fields[j])
								infields[i] = this_sample
					inline = infields[0]
					for k in range( 1, len(infields) ):
						inline = inline + "\t" + infields[k]

				outline = inline + "\n"
				out_vcf.write( outline )
	out_vcf.close()

if __name__=='__main__':
    main()


