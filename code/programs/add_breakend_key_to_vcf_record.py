#!/usr/bin/python
# python add_breakend_key_to_vcf_record.py -p position_flag
# cat test1.vcf | python add_breakend_key_to_vcf_record.py -p LEFT > test1_out.vcf
# cat test1.vcf | python add_breakend_key_to_vcf_record.py -p RIGHT > test1_out.vcf
# cat test1.vcf | python add_breakend_key_to_vcf_record.py -p LEFT -nohigh HIGH > test1_out2.vcf
# cat test1.vcf | python add_breakend_key_to_vcf_record.py -p LEFT -nohigh INTERMEDIATE > test1_out3.vcf

# This program reads in a VCF file from stdin, 
# adds a tab-delimited key of chrom pos1 pos2 to the front of each line, 
# then outputs to stdout.

# If position_flag is LEFT, then pos1 pos2 are both VCF.POS +- 100 bp.
# If position_flag is RIGHT, then pos1 pos2 are both VCF.INFO.END +- 100 bp,
# or is the position in the ALT field of a BND record.


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
def highest_tranche( tranche1, tranche2 ):

	highest_tranche = tranche2
	if ((tranche1 == 'HIGH') or (tranche2 == 'HIGH')):
		highest_tranche = 'HIGH'
	elif ((tranche1 == 'INTERMEDIATE') or (tranche2 == 'INTERMEDIATE')):
		highest_tranche = 'INTERMEDIATE'
	elif ((tranche1 == 'LOW') or (tranche2 == 'LOW')):
		highest_tranche = 'LOW'

	return highest_tranche

######################################################
def get_tranche( this_format, infields ):

	this_tranche = ''

	tranche_idx = -1
	format_fields = this_format.split(':')
	for i in range( 0, len(format_fields) ):
		this_format_field = format_fields[i]
		if (this_format_field == 'TRANCHE'):
			tranche_idx = i

	if (tranche_idx >= 0):
		if (len(infields) > 9):
			for sample_upto in range( 9, len(infields) ):
				this_sample = infields[sample_upto]
				if ((this_sample != './.') and (this_sample != '.')):
					this_sample_fields = this_sample.split(':')
					this_sample_tranche = this_sample_fields[tranche_idx]
					this_tranche = highest_tranche( this_sample_tranche, this_tranche )

	return this_tranche

######################################################
def construct_LEFT_key( this_chrom, this_pos, this_alts, this_info, window ):

	this_pos = int(this_pos)
	if (this_pos < 1):
		this_pos = 1
	pos1 = this_pos - 1 - window
	if (pos1 < 0):
		pos1 = 0
	pos2 = this_pos + window
	return_key = str(this_chrom) + "\t" + str(pos1) + "\t" + str(pos2)
	return return_key

######################################################
def construct_RIGHT_key( this_chrom, this_pos, this_alts, this_info, window ):

	# One BND:          [1:937030[G
	# One BND:          C]1:936513]
	# Multiple ALTS:    [1:937030[G,[1:162145669[C,C]1:936513]
	return_key = str(this_chrom) + "\t" + str(this_pos) + "\t" + str(this_pos)
	found_INFO_END = False
	this_info_fields = this_info.split(';')
	for this_info_field in this_info_fields:
		key_and_value = this_info_field.split('=')
		this_key = key_and_value[0]
		if (this_key == 'END'):
			found_INFO_END = True
			this_INFO_END = int(key_and_value[1])
			pos1 = this_INFO_END - 1 - window
			if (pos1 < 0):
				pos1 = 0
			pos2 = this_INFO_END + window
			return_key = str(this_chrom) + "\t" + str(pos1) + "\t" + str(pos2)
	if (found_INFO_END == False):
		# then look at first ALT
		this_alts_fields = this_alts.split(',')
		first_alt = this_alts_fields[0]
		BND_bracket = ']'
		BND_idx = first_alt.find( BND_bracket )
		if (BND_idx == -1):
			BND_bracket = '['
			BND_idx = first_alt.find( BND_bracket )
		if (BND_idx > -1):
			this_alt_field_bits = first_alt.split( BND_bracket )
			for this_bit in this_alt_field_bits:
				colon_idx = this_bit.find( ':' )
				if (colon_idx > -1):
					colon_bits = this_bit.split( ':' )
					this_RIGHT_chrom = str(colon_bits[0])
					this_RIGHT_pos = int(colon_bits[1])
					pos1 = this_RIGHT_pos - 1 - window
					if (pos1 < 0):
						pos1 = 0
					pos2 = this_RIGHT_pos + window
					return_key = str(this_RIGHT_chrom) + "\t" + str(pos1) + "\t" + str(pos2)
	return return_key


######################################################
def main():

	parser = argparse.ArgumentParser(description='Add a bed key to beginning of each line of a VCF file - either POS or INFO.END')
	parser.add_argument('-p', action="store", dest="position", required=True, help='LEFT or RIGHT, to specify whether key to add is POS or INFO.END')
	parser.add_argument('-nohigh', action="store", dest="nohigh", required=False, help='If this flag is present with HIGH, do not include TRANCHE=HIGH variants. If it is INTERMEDIATE, do not include TRANCHE=HIGH or TRANCHE=INTERMEDIATE.')
	parser.add_argument('-window', action="store", dest="window", required=False, help='If this flag is present then it specifies how many basepairs (bp) to add onto either side of the breakend point. Otherwise the default is 100 bp.')
	args = parser.parse_args()

	include_high_in_output = True
	include_intermediate_in_output = True
	if (args.nohigh is not None):
		include_high_in_output = False
		if (args.nohigh == 'INTERMEDIATE'):
			include_intermediate_in_output = False
	window = 100
	if (args.window is not None):
		window = int(args.window)

	# Read in the input VCF file from STDIN

	in_header = True
	for inline in sys.stdin:

		if (in_header == True):
			if (len(inline) >= 1):
				first_char = inline[0:1]
				if (first_char != '#'):
					in_header = False

		if (in_header):
			sys.stdout.write( inline )

		else: # in_header == False # We are processing VCF data records. We are no longer in the header part of the file.

			inline = inline.strip()
			infields = inline.split("\t")
			this_chrom = str(infields[0])
			this_pos = str(infields[1])
			this_alts = str(infields[4])
			this_info = str(infields[7])
			if (len(infields) > 8):
				this_format = str(infields[8])

			include_this_variant_in_output = True
			if (len(infields) > 8):
				if (include_high_in_output == False):
					this_tranche = get_tranche( this_format, infields )
					if (this_tranche == 'HIGH'):
						include_this_variant_in_output = False
					elif ((this_tranche == 'INTERMEDIATE') and (include_intermediate_in_output == False)):
						include_this_variant_in_output = False

			if (include_this_variant_in_output):
				key_to_add = '.' + "\t" + '.' + "\t" + '.'
				if (args.position == 'RIGHT'):
					key_to_add = construct_RIGHT_key( this_chrom, this_pos, this_alts, this_info, window )
				else: # args.position == 'LEFT'
					key_to_add = construct_LEFT_key( this_chrom, this_pos, this_alts, this_info, window )

				outline = key_to_add + "\t" + inline + "\n"
				sys.stdout.write( outline )

if __name__=='__main__':
    main()


