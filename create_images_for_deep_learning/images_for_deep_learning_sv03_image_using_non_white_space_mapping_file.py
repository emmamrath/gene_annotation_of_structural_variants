#!/usr/bin/python
# python3 images_for_deep_learning_sv03_image_using_non_white_space_mapping_file.py -i infile.tsv -o outfile.png -m map_file_for_chrom_pos_to_pixel_row_col.tsv
# python3 images_for_deep_learning_sv03_image_using_non_white_space_mapping_file.py -i AAAAA.gridss_annotated_genomic_regions_2_consequences.tsv -o AAAAA_SV_350x350.png -m non_white_out_all_mgrb_and_isks_map_350x350.txt
#
# module load python3/3.5.2
# module load python3/3.5.2-matplotlib
# module unload intel-fc intel-cc
# module load gcc/9.1.0
# pip3 install -v --no-binary :all: --prefix=/g/data3/wq2/users/emr913/software_various/python_libraries/python3.5.2 Pillow
# pip3 install -v --no-binary :all: --prefix=/g/data3/wq2/users/emr913/software_various/python_libraries/python3.5.2 argparse
# export PYTHONPATH=/g/data3/wq2/users/emr913/software_various/python_libraries/python3.5.2/lib/python3.5/site-packages
#
# python3 convert_gridss_SV_to_image.py -i ../vcf_TRANCHE2_annotations1_annotations2/AAAAA.gridss_annotated_genomic_regions_2.tsv -o ../plot_sv_square/AAAAA_sv_square.png

# Colours represent structural variant types: order of priority
# INS = red
# INDEL = darkred
# DEL = black
# INV = blue
# DUP = green
# no sv = white

# Mapping file is tab-delimited and has no header. First line is nucleotide_to_pixel_ratio in chrom:start-pos, num_rows, and num_cols for pixel image.
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
import math
import subprocess
import argparse
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

######################################################
def is_integer(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

######################################################
def init_img( img_width_pixels, img_height_pixels, pix ):

	img_array = np.zeros([img_height_pixels, img_width_pixels, 3], dtype=np.uint8, order='C')
	# 3,000,000,000 bytes = 3 GB memory needed

	for x in range( 0, img_height_pixels ):
		for y in range( 0, img_width_pixels ):
			img_array[x, y] = pix

	return img_array

######################################################
def get_col_numbers_for_col_headings( inline ):

	col = {}
	col['chrom'] = -1
	col['pos'] = -1
	col['end'] = -1
	col['svtype'] = -1
	col['svlen'] = -1
	col['tranche2'] = -1
	col['vaf'] = -1
	infields = inline.split("\t")
	for i in range( 0, len(infields) ):
		this_col = infields[i]
		if (this_col == 'CHROM'):
			col['chrom'] = i
		elif (this_col == 'POS'):
			col['pos'] = i
		elif (this_col == 'END'):
			col['end'] = i
		elif (this_col == 'SVTYPE'):
			col['svtype'] = i
		elif (this_col == 'SVLEN'):
			col['svlen'] = i
		elif (this_col == 'TRANCHE2'):
			col['tranche2'] = i
		elif (this_col == 'VAF'):
			col['vaf'] = i

	return col

######################################################
def extract_SV( inline, col ):

	infields = inline.split("\t")
	svtype = str( infields[ col['svtype'] ] )
	chrom = str( infields[ col['chrom'] ] )
	start = str( infields[ col['pos'] ] )
	end = str( infields[ col['end'] ] )
	svlen = str( infields[ col['svlen'] ] )
	tranche2 = str( infields[ col['tranche2'] ] )
	if (len(infields) > col['vaf']):
		vaf = str( infields[ col['vaf'] ] )

	if (svtype in [ 'INS', 'DEL', 'INDEL', 'DUP', 'INV' ]):
		start = int(start)
		end = int(end)
		if (svlen == ''):
			svlen = end - start + 1
		else:
			svlen = int(svlen)
		vaf = float(vaf)
	else:
		svtype = ''
		chrom = ''
		start = ''
		end = ''
		svlen = ''
		tranche2 = ''
		vaf = ''

	return svtype, chrom, start, end, svlen, tranche2, vaf

######################################################
def choose_pix( svtype, white_pix, black_pix, red_pix, darkred_pix, green_pix, blue_pix ):

	pix = white_pix
	if (svtype == 'INS'):
		pix = red_pix
	elif (svtype == 'DEL'):
		pix = black_pix
	elif (svtype == 'INDEL'):
		pix = darkred_pix
	elif (svtype == 'DUP'):
		pix = green_pix
	elif (svtype == 'INV'):
		pix = blue_pix

	return pix

######################################################
def add_pixels_to_pixel_position( start_y, start_x, pixels_to_add, max_x, max_y ):

	# x is pixel column, y is pixel row
	# This add_pixels_to_pixel_position assume 1 to max_col columns and 1 to max_row rows. It assumes 1-based coordinates.
	x = start_x
	y = start_y
	remaining_pixels_on_this_row = max_x - start_x
	if (pixels_to_add <= remaining_pixels_on_this_row):
		x = start_x + pixels_to_add
		y = start_y
	else:
		remaining_pixels_after_this_row = pixels_to_add - remaining_pixels_on_this_row
		remaining_rows = int( remaining_pixels_after_this_row / max_x )
		y = start_y + remaining_rows
		remaining_pixels_in_target_row = int( remaining_pixels_after_this_row % max_x )
		if (remaining_pixels_in_target_row == 0):
			x = max_x
		else:
			y = y + 1
		x = remaining_pixels_in_target_row

	return x, y

######################################################
def img_position( chrom, pos, map_chrom, map_start_pos, map_end_pos, map_row_start, map_col_start, map_row_end, map_col_end, nucleotide_to_pixel_ratio, img_width_pixels, img_height_pixels ):

	# x is pixel column, y is pixel row
	x = img_width_pixels
	y = img_height_pixels

	found_position = False
	i = 0
	while (found_position == False):
		if (chrom == map_chrom[i]):

			if (map_start_pos[i] <= pos <= map_end_pos[i]):
				# is in here
				found_position = True

				distance_of_pos_from_map_start_pos_in_bp = pos - map_start_pos[i]
				distance_of_pos_from_map_start_pos_in_pixels = int( distance_of_pos_from_map_start_pos_in_bp / nucleotide_to_pixel_ratio )
				x, y = add_pixels_to_pixel_position( map_row_start[i], map_col_start[i], distance_of_pos_from_map_start_pos_in_pixels, img_width_pixels, img_height_pixels )
		i = i + 1
		if ( (found_position == False) and (i >= len(map_chrom)) ):
			print( "Error: didn't find SV position in mapping file:" + str(chrom) + ":" + str(start) + "-" + str(end) )
			found_position = True

	# This img_position function assumes 1 to max_col columns and 1 to max_row rows. It assumes 1-based coordinates.
	# Until these last lines. After working in 1-based, it converts 1-based to 0-based, and returns coordinates in 0-based.
	x = x - 1 # map file is 1 to 350, image is 0 to 349, convert 1-based to 0-based
	y = y - 1 # map file is 1 to 350, image is 0 to 349, convert 1-based to 0-based

	return x, y

######################################################
def add_sv_to_image( img_array, start_x, start_y, end_x, end_y, svtype_pix, img_width_pixels ):

	# x is pixel column, y is pixel row
	for y in range( start_y, (end_y+1) ):
		if ((y == start_y) and (y == end_y)):
			for x in range( start_x, (end_x+1) ):
				img_array[y,x] = svtype_pix
		elif (y == end_y):
			for x in range( 0, (end_x+1) ):
				img_array[y,x] = svtype_pix
		elif (y == start_y):
			for x in range( start_x, img_width_pixels):
				img_array[y,x] = svtype_pix
		else:
			for x in range( 0, img_width_pixels ):
				img_array[y,x] = svtype_pix

	return img_array

######################################################
def read_map_file( in_map_file ):

	img_width_pixels = 0
	img_height_pixels = 0
	map_chrom = []
	map_start_pos = []
	map_end_pos = []
	map_row_start = []
	map_col_start = []
	map_row_end = []
	map_col_end = []
	chroms = {}

	infile = open( in_map_file, 'r')
	is_first_line = True
	for inline in infile:
		inline = inline.strip()
		if (inline != ''):
			infields = inline.split('\t')
			if (is_first_line):
				nucleotide_to_pixel_ratio = int(infields[0])
				img_width_pixels = int(infields[1])
				img_height_pixels = int(infields[2])
				is_first_line = False
			else:
				this_chrom = str(infields[0])
				chroms[this_chrom] = this_chrom
				map_chrom.append( str(infields[0]) )
				map_start_pos.append( int(infields[1]) )
				map_end_pos.append( int(infields[2]) )
				map_row_start.append( int(infields[3]) )
				map_col_start.append( int(infields[4]) )
				map_row_end.append( int(infields[5]) )
				map_col_end.append( int(infields[6]) )

	return nucleotide_to_pixel_ratio, img_width_pixels, img_height_pixels, map_chrom, map_start_pos, map_end_pos, map_row_start, map_col_start, map_row_end, map_col_end, chroms

######################################################
def main():

        # what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Convert a tab-separated file containing structural variants to a square image.')
	parser.add_argument('-i', action="store", dest="infile", required=True, help='Input tab-separated file of SVs')
	parser.add_argument('-o', action="store", dest="outfile", required=True, help='Output PNG image')
	parser.add_argument('-m', action="store", dest="in_map_file", required=True, help='Input file containing the mapping of genomic chrom:start-end to image pixel row-col')
	args = parser.parse_args()

	infile = args.infile
	outfile = args.outfile

	nucleotide_to_pixel_ratio, img_width_pixels, img_height_pixels, map_chrom, map_start_pos, map_end_pos, map_row_start, map_col_start, map_row_end, map_col_end, chroms = read_map_file( args.in_map_file )

	white_pix = np.array([255, 255, 255], dtype=np.uint8) # no sv
	black_pix = np.array([0, 0, 0], dtype=np.uint8) # del
	red_pix = np.array([255, 0, 0], dtype=np.uint8) # ins
	darkred_pix = np.array([127, 0, 0], dtype=np.uint8) # indel
	green_pix = np.array([0, 255, 0], dtype=np.uint8) # dup
	blue_pix = np.array([0, 0, 255], dtype=np.uint8) # inv

	img_array = init_img( img_width_pixels, img_height_pixels, white_pix )

	# Read the input variants tsv file 5 times, one time for each of DUP, INV, DEL, INDEL, INS
	# By doing INS last, that colour will be on top for pixels having multiple SVs.

	list_svtypes = [ 'DUP', 'INV', 'DEL', 'INDEL', 'INS' ]
	list_svtype_pix = [ green_pix, blue_pix, black_pix, darkred_pix, red_pix ]

	for i in range( 0, len(list_svtypes) ):
		this_svtype = list_svtypes[i]
		this_svtype_pix = list_svtype_pix[i]

		infile_vars = open(infile, 'r')
		in_header = True
		for inline in infile_vars:

			inline = inline.strip()
			if (inline != ''):

				if (in_header == True):
					col = get_col_numbers_for_col_headings( inline )
					in_header = False

				else: # (in_header == False):

					svtype, chrom, start, end, svlen, tranche2, vaf = extract_SV( inline, col )
					if ((svtype == this_svtype) and (chrom in chroms) and (tranche2 in [ 'HIGH', 'INTERMEDIATE' ]) and (vaf >= 0.25)):
						start_x, start_y = img_position( chrom, start, map_chrom, map_start_pos, map_end_pos, map_row_start, map_col_start, map_row_end, map_col_end, nucleotide_to_pixel_ratio, img_width_pixels, img_height_pixels )
						end_x, end_y = img_position( chrom, end, map_chrom, map_start_pos, map_end_pos, map_row_start, map_col_start, map_row_end, map_col_end, nucleotide_to_pixel_ratio, img_width_pixels, img_height_pixels )
						img_array = add_sv_to_image( img_array, start_x, start_y, end_x, end_y, this_svtype_pix, img_width_pixels )

	img = Image.fromarray(img_array)
	img.save(outfile)

if __name__=='__main__':
    main()

