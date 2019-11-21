#!/usr/bin/python
# python images_for_deep_learning_sv02_create_non_white_space_mapping_file.py -i input_bed -o output_mapping_file
# python images_for_deep_learning_sv02_create_non_white_space_mapping_file.py -i non_white_out_all_mgrb_and_isks.txt -o non_white_out_all_mgrb_and_isks_map_350x350.txt

# Sort input bed file by chromosome. It is assumed that positions within chromosome are already sorted. 
# Count number of nucleotides in it.
# Map nucleotide positions to a 350 x 350 array.
# Output the mapping.
# This mapping will be used to convert a VCF file to a 350 x 350 image.
# The entire 3 billion nucleotide genome could have been mapped. However, most of it doesn't contain any variants. 
# This program reads and maps only those genome positions that have variants.
# Even then, there will be more genomic nucleotides to map to pixels than there are available pixels in a 350 x 350 image.
# This program calculates the nucleotide_to_pixel_ratio and maps or compresses multiple bp into one pixel.

# Mapping file is tab-delimited and has no header. First line is nucleotide_to_pixel_ratio, num_rows, and num_cols for pixel image.
# Columns are chrom, start_pos, end_pos, map_row_start, map_col_start, map_row_end, map_col_end
# 14214	350	350
# 1	843216	843248	1	1	1	1
# 1	869460	870342	1	1	1	1
# 1	884041	884110	1	1	1	1
# ...
# Y	22661495	22661520	350	52	350	52
# Y	24417006	24417026	350	52	350	52
# Y	28787561	28802853	350	52	350	53
# MT	1	16569	350	53	350	54

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2019, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import argparse
import math

######################################################
def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

######################################################
def read_input( in_bed ):

	in_chrom = []
	in_start = []
	in_end = []
	infile = open( in_bed, 'r')
	for inline in infile:
		inline = inline.strip()
		if (inline != ''):
			infields = inline.split('\t')
			in_chrom.append( str(infields[0]) )
			in_start.append( int(infields[1]) )
			in_end.append( int(infields[2]) )

	return in_chrom, in_start, in_end

######################################################
def sort_input( in1_chrom, in1_start, in1_end ):

	order_of_chromosomes = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']

	in2_chrom = []
	in2_start = []
	in2_end = []
	for this_chrom in order_of_chromosomes:
		for i in range( 0, len(in1_chrom) ):
			if (in1_chrom[i] == this_chrom):
				in2_chrom.append( in1_chrom[i] )
				in2_start.append( in1_start[i] )
				in2_end.append( in1_end[i] )

	return in2_chrom, in2_start, in2_end

######################################################
def count_nucleotides( in_start, in_end ):

	num_nucleotides = 0
	for i in range( 0, len(in_start) ):
		num_nucleotides = num_nucleotides + in_end[i] - in_start[i] + 1

	return num_nucleotides

######################################################
def map_positions_to_pixels( in_chrom, in_start, in_end, num_rows, num_cols, nucleotide_to_pixel_ratio ):

	map_row_start = [0] * len(in_chrom)
	map_col_start = [0] * len(in_chrom)
	map_row_end = [0] * len(in_chrom)
	map_col_end = [0] * len(in_chrom)

	row_upto = 1
	col_upto = 1
	for i in range( 0, len(in_chrom)):

		chrom = str(in_chrom[i])
		start_pos = int(in_start[i])
		end_pos = int(in_end[i])

		map_row_start[i] = row_upto
		map_col_start[i] = col_upto

		remaining_unmapped_bp = end_pos - start_pos + 1
		remaining_unmapped_pixels = int(remaining_unmapped_bp / nucleotide_to_pixel_ratio)
		if (int(remaining_unmapped_bp % nucleotide_to_pixel_ratio) > 0):
			remaining_unmapped_pixels = remaining_unmapped_pixels + 1
		remaining_cols_in_row = num_cols - col_upto - 1
		if (remaining_unmapped_pixels <= remaining_cols_in_row):
			map_row_end[i] = row_upto
			map_col_end[i] = col_upto + remaining_unmapped_pixels - 1
			col_upto = col_upto + remaining_unmapped_pixels - 1
		else:
			remaining_unmapped_pixels = remaining_unmapped_pixels - remaining_cols_in_row
			additional_rows = int(math.ceil( float(remaining_unmapped_pixels) / float(num_cols) ))
			map_row_end[i] = row_upto + additional_rows
			row_upto = row_upto + additional_rows
			additional_cols = int(remaining_unmapped_pixels % num_cols)
			if (additional_cols == 0):
				col_upto = num_cols
				map_col_end[i] = num_cols
			else:
				col_upto = additional_cols
				map_col_end[i] = additional_cols

	return map_row_start, map_col_start, map_row_end, map_col_end

######################################################
def write_output_map( out_map, nucleotide_to_pixel_ratio, num_rows, num_cols, in_chrom, in_start, in_end, map_row_start, map_col_start, map_row_end, map_col_end ):

	out_map_file = open(out_map, 'w')

	outline = str(nucleotide_to_pixel_ratio) + "\t" + str(num_rows) + "\t" + str(num_cols) + "\n"
	out_map_file.write( outline )

	for i in range( 0, len(in_chrom) ):
		outline = str(in_chrom[i]) + "\t" + str(in_start[i]) + "\t" + str(in_end[i]) + "\t" + str(map_row_start[i]) + "\t" + str(map_col_start[i]) + "\t" + str(map_row_end[i]) + "\t" + str(map_col_end[i]) + "\n"
		out_map_file.write( outline )

	out_map_file.close()

	return

######################################################
def main():

	parser = argparse.ArgumentParser(description='Read in BED file and sort. Map BED file locations to a 350 x 350 array.')
	parser.add_argument('-i', action="store", dest="in_bed", required=True, help='Input BED file')
	parser.add_argument('-o', action="store", dest="out_map", required=True, help='Output mapping file')
	args = parser.parse_args()

	num_rows = 350
	num_cols = 350

	in1_chrom, in1_start, in1_end = read_input( args.in_bed )
	in2_chrom, in2_start, in2_end = sort_input( in1_chrom, in1_start, in1_end )

	num_nucleotides = count_nucleotides( in2_start, in2_end )
	num_pixels = num_rows * num_cols
	nucleotide_to_pixel_ratio = int(num_nucleotides / num_pixels) + 1

	map_row_start, map_col_start, map_row_end, map_col_end = map_positions_to_pixels( in2_chrom, in2_start, in2_end, num_rows, num_cols, nucleotide_to_pixel_ratio )

	write_output_map( args.out_map, nucleotide_to_pixel_ratio, num_rows, num_cols, in2_chrom, in2_start, in2_end, map_row_start, map_col_start, map_row_end, map_col_end )

if __name__=='__main__':
    main()


