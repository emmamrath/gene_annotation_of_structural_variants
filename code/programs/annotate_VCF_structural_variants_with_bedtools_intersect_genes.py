#!/usr/bin/python
# python annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i intersect_file_of_intersecting_genes -t type_of_intersecting_genes -d description_of_intersecting_genes
# cat BNDtestdata1_highlevel_calls.vcf | python annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i BNDtestdata1_highlevel_calls_intersect_genes.txt -t HITGENES -d 'Intersecting genes' > BNDtestdata1_highlevel_calls_HITGENES.vcf

# This program takes 1 VCF file on input from stdin, and writes out an annotated VCF file to stdout.
# Another input to this program is a bed file from bedtools intersect containing the genes that are hit by the VCF structural variants.
# A VCF record may hit 0, 1 or more genes. There is a record in the bedtools intersect file for each hit.
# The VCF file and intersect file are in the same order, and have the same ID values.
# However, the chromosomal position keys of the 2 files do not have the same values for the same variants.
# This program annotates the VCF.INFO field with type_of_intersecting_genes=gene_id.
# where the type_of_intersecting_genes is given on the input command line and the gene_id is from each line of the intersect_file_of_intersecting_genes.
# If a VCF variant has no entries in the intersect_file_of_intersecting_genes, then it will not get an annotation.
# If a VCF variant has more than one entry in the intersect_file_of_intersecting_gene, then it will get only one field for this type_of_intersecting_genes in the VCF.INFO field
# and the multiple entries will appear separated by a bar | character.

# This program assumes that there is only one ALT per VCF record.

# Here mapping of VCF chromosomal position keys to the bed file chromosomal position keys, which was created by previous program convert_VCF_structural_variants_to_BED_format.py
# VCF.ALT	==>		CHR		START		END
# =======			===		=====		===	
# <DEL>		==>		VCF.CHROM	VCF.POS-1	VCF.INFO.END
# <INS>		==>		VCF.CHROM	VCF.POS-1	VCF.POS
# <INDEL>	==>		VCF.CHROM	VCF.POS-1	VCF.INFO.END
# <DUP:TANDEM>	==>		VCF.CHROM	VCF.INFO.END-1	VCF.INFO.END
# <DUP:INS>	==>		VCF.CHROM	VCF.INFO.END-1	VCF.INFO.END
# <INV>		==>		VCF.CHROM	VCF.POS-1	VCF.INFO.END

# This program does not match up the VCF records to the bed file records by chromosomal positions.
# Instead, this program matches up the VCF records to the bed file records by ID, which are not in alphabetical order.
# However, both file have records in the same ID order.
# The VCF file has 1 record per ID and ID is the 3rd field. 
# The bed file has 0, 1 or multiple records per ID, and ID is the 4th field.

# If the intersect_file_of_intersecting_genes has 0 in the 14th field (which is 0 base-pairs overlapping with a gene), 
# then ignore this intersect_file_of_intersecting_genes record.

# Example input VCF file:
# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  AABAF.sorted.dupmarked.bam
# 1       106761073       gridss10_5537   A       <INDEL> 1363    .       SVTYPE=INDEL;END=106761089;SVLEN=64;TRANCHE=HIGH;INSSEQ=TGCCATGTTTTTGCATTTTTGTGCTTTTTGTTTGTAATTATGCTGTTTAAAATGACACCCAAGCACAGTGCTGAAGTGCT        GT      1/.
# 6       94381988        gridss118_7208  T       <DUP:TANDEM>    1100    .       SVTYPE=DUP;END=94382102;SVLEN=114;TRANCHE=HIGH  GT      1/.
# 6       124259133       gridss121_3513  T       <INS>   1064    .       SVTYPE=INS;END=124259133;SVLEN=35;TRANCHE=HIGH;INSSEQ=GAATCTGATGCACACAGCTTCATGAGAAAACAAGC       GT      1/.
# 7       102213180       gridss137_8268  C       <DUP:INS>       492.93  LOW_QUAL;SINGLE_ASSEMBLY        SVTYPE=DUP;END=115930270;SVLEN=13717100;TRANCHE=LOW;INSSEQ=TTCCCCCCCC   GT      1/.
# 1       157568893       gridss15_7528   A       <DEL>   176.63  LOW_QUAL;SINGLE_ASSEMBLY        SVTYPE=DEL;END=157568929;SVLEN=-36;TRANCHE=LOW  GT      1/.
# 7       10000381        notgridss100_100        G       <INV>   1697    ASSEMBLY_ONLY   SVTYPE=INV;END=10000481;TRANCHE=HIGH    GT      1/.

# Example of input intersecting_genes file:
# 7	97872287	97872288	gridss136_8207	C	<DUP:TANDEM>	598.71	LOW_QUAL;SINGLE_ASSEMBLY	SVTYPE=DUP;END=97872288;SVLEN=196;TRANCHE=INTERMEDIATE	7	97844754	97881563	Q7Z6L1-TECPR1	1
# 7	97872287	97872288	gridss136_8207	C	<DUP:TANDEM>	598.71	LOW_QUAL;SINGLE_ASSEMBLY	SVTYPE=DUP;END=97872288;SVLEN=196;TRANCHE=INTERMEDIATE	7	97855903	97881563	Q7Z6L1-3-TECPR1	1
# 7	115930269	115930270	gridss137_8268	C	<DUP:INS>	492.93	LOW_QUAL;SINGLE_ASSEMBLY	SVTYPE=DUP;END=115930270;SVLEN=13717100;TRANCHE=LOW;INSSEQ=TTCCCCCCCC	7	115928126	115977312	EF070119-	1
# 7	115930269	115930270	gridss137_8268	C	<DUP:INS>	492.93	LOW_QUAL;SINGLE_ASSEMBLY	SVTYPE=DUP;END=115930270;SVLEN=13717100;TRANCHE=LOW;INSSEQ=TTCCCCCCCC	7	115928126	115977312	EF070117-	1
# 7	115930269	115930270	gridss137_8268	C	<DUP:INS>	492.93	LOW_QUAL;SINGLE_ASSEMBLY	SVTYPE=DUP;END=115930270;SVLEN=13717100;TRANCHE=LOW;INSSEQ=TTCCCCCCCC	7	115928126	115977312	EF070122-	1

# Example output file:


__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2017, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import commands
import argparse

######################################################
def is_integer(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

######################################################
def read_intersects_for_this_variant( VCF_ID ):

	global intersect_file_of_intersecting_genes
	global intersect_line
	global intersect_EOF
	global intersect_ID

	genes_that_intersect_this_VCF_record = []
	genes_that_intersect_this_VCF_record_dict = {}
	intersect_fields = intersect_line.split("\t")

	keep_reading_intersect_genes_for_this_variant = True
	while (keep_reading_intersect_genes_for_this_variant):

		# Add this intersect gene to the list of genes for this variant

		this_intersect_gene = str(intersect_fields[12])
		if (this_intersect_gene not in genes_that_intersect_this_VCF_record_dict):
			genes_that_intersect_this_VCF_record.append( this_intersect_gene )
			genes_that_intersect_this_VCF_record_dict[this_intersect_gene] = this_intersect_gene

		# Read the next intersect gene.
		# Is it also a gene for this variant?

		intersect_line = intersect_file_of_intersecting_genes.readline()
		if (intersect_line):
			intersect_fields = intersect_line.split("\t")
			intersect_ID = str(intersect_fields[3])
			if (intersect_ID != VCF_ID):
				keep_reading_intersect_genes_for_this_variant = False
		else:
			intersect_EOF = True
			keep_reading_intersect_genes_for_this_variant = False

	return genes_that_intersect_this_VCF_record

######################################################
def construct_annotated_VCF_record( inline, type_of_intersecting_genes, genes_that_intersect_this_VCF_record ):

	# Insert the intersecting_genes annotation into this VCF record

	infields = inline.split("\t")
	infields[ len(infields)-1 ] = infields[ len(infields)-1 ].strip() # Remove the carriage return on the last field of the VCF record
	list_of_intersection_genes = genes_that_intersect_this_VCF_record[0]
	if (len(genes_that_intersect_this_VCF_record) > 1):
		for i in range( 1, len(genes_that_intersect_this_VCF_record) ):
			list_of_intersection_genes = list_of_intersection_genes + '|' + genes_that_intersect_this_VCF_record[i]
	info_string_for_intersecting_genes = ''
	if (list_of_intersection_genes != ''):
		info_string_for_intersecting_genes = type_of_intersecting_genes + '=' + list_of_intersection_genes
	new_info_field = infields[7]
	if (new_info_field == ''):
		new_info_field = info_string_for_intersecting_genes
	else:
		new_info_field = new_info_field + ';' + info_string_for_intersecting_genes

	outline = infields[0]
	for i in range( 1, len(infields) ):
		if (i == 7):
			outline = outline + "\t" + new_info_field
		else:
			outline = outline + "\t" + str(infields[i])
	outline = outline + "\n"

	return outline

######################################################
def main():

	global intersect_file_of_intersecting_genes
	global intersect_line
	global intersect_EOF
	global intersect_ID

	# What input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in an uncompressed VCF file from stdin, annotate the INFO fields with any intersecting genes, then output annotated VCF file to stdout.')
	parser.add_argument('-i', action="store", dest="intersect_file_of_intersecting_genes", required=True, help='bedtools intersect file containing one line per gene hit, in the same VCF variant order as the input VCF file')
	parser.add_argument('-t', action="store", dest="type_of_intersecting_genes", required=True, help='The INFO fields will be annotated with this key, eg. HITGENES or HITREGULATORIES')
	parser.add_argument('-d', action="store", dest="description_of_intersecting_genes", required=True, help='The description of this INFO field that will appear in the VCF header record for it')
	args = parser.parse_args()
	args.type_of_intersecting_genes = str(args.type_of_intersecting_genes)
	args.description_of_intersecting_genes = str(args.description_of_intersecting_genes)

	# Open the intersect_file_of_intersecting_genes and read its first record

	intersect_file_of_intersecting_genes = open(args.intersect_file_of_intersecting_genes)
	intersect_line = intersect_file_of_intersecting_genes.readline()
	intersect_EOF = False
	intersect_ID = ''
	if (intersect_line):
		intersect_fields = intersect_line.split("\t")
		intersect_ID = str(intersect_fields[3])
	else:
		intersect_EOF = True

	# Read in the input VCF file from STDIN and process each record
	# The intersect_file_of_intersecting_genes file will be read line by line to keep up with the processing of the records in this input VCF file.

	in_header = True
	have_outputted_new_info_headers = False
	previous_header_type = ''
	new_info_headers = []
	new_info_headers.append('##INFO=<ID=' + args.type_of_intersecting_genes + ',Number=A,Type=String,Description="' + args.description_of_intersecting_genes + '">\n')

	for inline in sys.stdin:

		if (in_header == True):
			if (len(inline) >= 1):
				first_char = inline[0:1]
				if (first_char != '#'):
					in_header = False

		if (in_header == True):

			# At the end of existing INFO headers, output new INFO header

			if (len(inline) >= 6):
				if (inline[0:6] == '#CHROM'):
					if (have_outputted_new_info_headers == False):
						for header_line in new_info_headers:
							sys.stdout.write( header_line )
						have_outputted_new_info_headers = True
			split_header = inline.split('=')
			this_header_type = split_header[0]
			if ((previous_header_type == '##INFO') and (this_header_type != '##INFO')):
				if (have_outputted_new_info_headers == False):
					for header_line in new_info_headers:
						sys.stdout.write( header_line )
					have_outputted_new_info_headers = True
			previous_header_type = this_header_type

			# Output header lines as-is

			sys.stdout.write( inline )

		if (in_header == False): # We are processing VCF data records. We are no longer in the header part of the file.

			# Process this VCF record.
			# Collect all the 0, 1 or multiple intersect records for it.
			# Then write out this VCF record, with the intersect genes as annotations in the INFO field.

			infields = inline.split("\t")
			VCF_ID = infields[2]

			if (VCF_ID == intersect_ID):

				# This VCF record has one or more hit-genes in the intersect file.
				# Read all the hit-genes (all the intersect records for this VCF variant),
				# annotate its INFO field with these genes, then write out this annotated VCF record.

				genes_that_intersect_this_VCF_record = read_intersects_for_this_variant( VCF_ID )
				new_VCF_record = construct_annotated_VCF_record( inline, args.type_of_intersecting_genes, genes_that_intersect_this_VCF_record )
				sys.stdout.write( new_VCF_record )

			else:

				# This VCF record doesn't have any hit-genes in the intersect file.
				# So write out this VCF record as-is.

				sys.stdout.write( inline )

	intersect_file_of_intersecting_genes.close()


if __name__=='__main__':
    main()


