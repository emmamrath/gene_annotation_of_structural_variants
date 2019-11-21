#!/usr/bin/python
# python3 convert_gridss_SV_to_image.py -i infile.vcf -o outfile.png -l chrom_lengths.tsv
# python3 convert_gridss_SV_to_image.py -i AAAAA.gridss_annotated_genomic_regions_2_consequences.tsv -o AAAAA_SV_img_INS_DEL_INDEL_DUP_INV.png

# Note, this program produces an image that is too big (~ 53000 x 53000 pixels) to view on some laptops.
# Note, this program has had a fix applied to function add_sv_to_image that has not yet been tested.

# module load python3/3.5.2
# module load python3/3.5.2-matplotlib
# module unload intel-fc intel-cc
# module load gcc/9.1.0
# pip3 install -v --no-binary :all: --prefix=/g/data3/wq2/users/emr913/software_various/python_libraries/python3.5.2 Pillow
# pip3 install -v --no-binary :all: --prefix=/g/data3/wq2/users/emr913/software_various/python_libraries/python3.5.2 argparse
# export PYTHONPATH=/g/data3/wq2/users/emr913/software_various/python_libraries/python3.5.2/lib/python3.5/site-packages
#
# python3 convert_gridss_SV_to_image.py -i ../vcf_TRANCHE2_annotations1_annotations2/AAAAA.gridss_annotated_genomic_regions_2.tsv -o ../plot_sv_square/AAAAA_sv_square.png

# no sv = white
# del = black
# ins = red
# indel = darkred
# dup = green
# inv = blue

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
def add_pixel_bits( bit1, bit2 ):

	bit = bit1
	if (bit1 > bit2):
		bit = int(round(bit1 - bit2)/2) + bit2
	else:
		bit = int(round(bit2 - bit1)/2) + bit1
	bit = bit.astype('uint8')

	return bit

######################################################
def add_pixels( pix1, pix2 ):

	pix = pix1
	r1 = pix1[0]
	g1 = pix1[1]
	b1 = pix1[2]
	r2 = pix2[0]
	g2 = pix2[1]
	b2 = pix2[2]
	r = add_pixel_bits( r1, r2 )
	g = add_pixel_bits( g1, g2 )
	b = add_pixel_bits( b1, b2 )
	pix[0] = r
	pix[1] = g
	pix[2] = b

	return pix

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

	if (svtype in [ 'INS', 'DEL', 'INDEL', 'DUP', 'INV' ]):
		do_nothing = 1
	else:
		svtype = ''
		chrom = ''
		start = ''
		end = ''
		svlen = ''
		tranche2 = ''

	return svtype, chrom, start, end, svlen, tranche2

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
def add_pixel_bits( bit1, bit2 ):

	bit = bit1
	if (bit1 > bit2):
		bit = int(round(bit1 - bit2)/2) + bit2
	else:
		bit = int(round(bit2 - bit1)/2) + bit1
	bit = bit.astype('uint8')

	return bit

######################################################
def merge_pixels( pix1, pix2 ):

	new_pix = pix1
	r1 = pix1[0]
	g1 = pix1[1]
	b1 = pix1[2]
	r2 = pix2[0]
	g2 = pix2[1]
	b2 = pix2[2]
	r = add_pixel_bits( r1, r2 )
	g = add_pixel_bits( g1, g2 )
	b = add_pixel_bits( b1, b2 )
	new_pix[0] = r
	new_pix[1] = g
	new_pix[2] = b

	return new_pix

######################################################
def add_pixel_to_img_pixel( pix1, pix2, white_pix ):

	new_pix = pix1
	#if (pix1 == white_pix):
	if (np.allclose( pix1, white_pix )):
		new_pix = pix2
	else:
		new_pix = merge_pixels( pix1, pix2 )

	return new_pix

######################################################
def img_position( chroms, chrom, pos, img_width_pixels, img_height_pixels ):

	x = img_width_pixels - 1
	y = img_height_pixels - 1

	global_pos = chroms[chrom] + int(pos) - 1

	whole_num = int( global_pos / img_width_pixels )
	remainder = global_pos % img_width_pixels

	x = whole_num
	y = remainder - 1
	if (remainder == 0):
		x = whole_num - 1
		y = img_width_pixels - 1

	return x, y

######################################################
def add_sv_to_image( img_array, img_width_pixels, img_height_pixels, chroms, svtype, chrom, start, end, svlen, svtype_pix, white_pix, black_pix, red_pix, darkred_pix, green_pix, blue_pix ):

	start_x, start_y = img_position( chroms, chrom, start, img_width_pixels, img_height_pixels )
	end_x, end_y = img_position( chroms, chrom, end, img_width_pixels, img_height_pixels )

	# This program had x and y mixed up as row and column. Wrong code has been commented out. Next code below has not been tested.
	#for x in range( start_x, (end_x+1) ):
	#	if ((x == start_x) and (x == end_x)):
	#		for y in range( start_y, (end_y+1) ):
	#			img_array[x,y] = add_pixel_to_img_pixel( img_array[x,y], svtype_pix, white_pix )
	#	elif (x == end_x):
	#		for y in range( 0, (end_y+1) ):
	#			img_array[x,y] = add_pixel_to_img_pixel( img_array[x,y], svtype_pix, white_pix )
	#	elif (x == start_x):
	#		for y in range( start_y, img_width_pixels):
	#			img_array[x,y] = add_pixel_to_img_pixel( img_array[x,y], svtype_pix, white_pix )
	#	else:
	#		for y in range( 0, img_width_pixels ):
	#			img_array[x,y] = add_pixel_to_img_pixel( img_array[x,y], svtype_pix, white_pix )

	# x is column, y is row
	for y in range( start_y, (end_y+1) ):
		if ((y == start_y) and (y == end_y)):
			for x in range( start_x, (end_x+1) ):
				img_array[y,x] = add_pixel_to_img_pixel( img_array[y,x], svtype_pix, white_pix )
		elif (y == end_y):
			for x in range( 0, (end_x+1) ):
				img_array[y,x] = add_pixel_to_img_pixel( img_array[y,x], svtype_pix, white_pix )
		elif (y == start_y):
			for x in range( start_x, img_width_pixels):
				img_array[y,x] = add_pixel_to_img_pixel( img_array[y,x], svtype_pix, white_pix )
		else:
			for x in range( 0, img_width_pixels ):
				img_array[y,x] = add_pixel_to_img_pixel( img_array[y,x], svtype_pix, white_pix )

	return img_array

######################################################
def calc_chrom_positions_in_img( in_lengths_file, img_width_pixels, img_height_pixels ):

	chrs = [ '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT' ]
	l = {}
	l['1'] = 249904550
	l['2'] = 243199373
	l['3'] = 198022430
	l['4'] = 191535534
	l['5'] = 180915260
	l['6'] = 171115067
	l['7'] = 159321559
	l['8'] = 146440111
	l['9'] = 141696573
	l['10'] = 135534747
	l['11'] = 135046619
	l['12'] = 133851895
	l['13'] = 115169878
	l['14'] = 107349540
	l['15'] = 102531392
	l['16'] = 90354753
	l['17'] = 81529607
	l['18'] = 78081510
	l['19'] = 59380841
	l['20'] = 63025520
	l['21'] = 48157577
	l['22'] = 51304566
	l['X'] = 155270560
	l['Y'] = 59373566
	l['MT'] = 16569

	chroms = {}
	pos_upto = 0

	#infile = open(in_lengths_file, 'r')
	#for inline in infile:
	#	inline = inline.strip()
	#	if (inline != ''):
	#		infields = inline.split("\t")
	#		chrom = str(infields[0])
	#		chrom_length = int(infields[1])
	#		chroms[chrom] = pos_upto
	#		pos_upto = pos_upto + chrom_length

	for chrom in chrs:
		chroms[chrom] = pos_upto
		pos_upto = pos_upto + l[chrom]

	return chroms

######################################################
def main():

        # what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Convert a vcf containing structural variants to a square image.')
	parser.add_argument('-i', action="store", dest="infile", required=True, help='Input VCF file')
	parser.add_argument('-o', action="store", dest="outfile", required=True, help='Output PNG image')
	# parser.add_argument('-l', action="store", dest="in_lengths_file", required=False, help='Input file containing the length of each chromosome')
	args = parser.parse_args()

	infile = args.infile
	outfile = args.outfile
	in_lengths_file = 'dummy' # args.in_lengths_file
	img_width_pixels = 55670
	img_height_pixels = 55670

	chroms = calc_chrom_positions_in_img( in_lengths_file, img_width_pixels, img_height_pixels )

	white_pix = np.array([255, 255, 255], dtype=np.uint8) # no sv
	black_pix = np.array([0, 0, 0], dtype=np.uint8) # del
	red_pix = np.array([255, 0, 0], dtype=np.uint8) # ins
	darkred_pix = np.array([127, 0, 0], dtype=np.uint8) # indel
	green_pix = np.array([0, 255, 0], dtype=np.uint8) # dup
	blue_pix = np.array([0, 0, 255], dtype=np.uint8) # inv

	# red = ins/indel/del, green = inv, blue = dup
	# order of priority in showing SV when they overlap on a pixel is: ins, indel, del, inv, dup

	img_array = init_img( img_width_pixels, img_height_pixels, white_pix )

	infile = open(infile, 'r')
	in_header = True
	for inline in infile:

		inline = inline.strip()
		if (inline != ''):

			if (in_header == True):
				col = get_col_numbers_for_col_headings( inline )
				in_header = False

			else: # (in_header == False):

				svtype, chrom, start, end, svlen, tranche2 = extract_SV( inline, col )
				if ((svtype in [ 'INS', 'DEL', 'INDEL', 'DUP', 'INV' ]) and (chrom in chroms) and (tranche2 in [ 'HIGH', 'INTERMEDIATE' ])):
					svtype_pix = choose_pix( svtype, white_pix, black_pix, red_pix, darkred_pix, green_pix, blue_pix )
					img_array = add_sv_to_image( img_array, img_width_pixels, img_height_pixels, chroms, svtype, chrom, start, end, svlen, svtype_pix, white_pix, black_pix, red_pix, darkred_pix, green_pix, blue_pix )

	img = Image.fromarray(img_array)
	img.save(outfile)

if __name__=='__main__':
    main()

