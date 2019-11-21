#!/usr/bin/python
# python convert_VCF_structural_variants_to_BED_format.py [-f filter_lengths_less_than_this] [-w window_extension]
# cat BNDtestdata1_highlevel_calls.vcf | python convert_VCF_structural_variants_to_BED_format.py -f 50 > BNDtestdata1_highlevel_calls.bed
# cat BNDtestdata1_highlevel_calls.vcf | python convert_VCF_structural_variants_to_BED_format.py -w 20000 > BNDtestdata1_highlevel_calls_w20000.bed

# This program takes 1 VCF file on input from stdin, and writes out a BED file to stdout containing only certain structural variants (SVs).
# The input VCF is 1-based and the output BED is 0-based.
# It optionally filters out SVs whose VCF.INFO.SVLEN length is less than that specified.
# -f filters out variant less than 50 bp long
# -w extends the window, eg. -w 20000 makes the window be 20000 basepairs on either side of the mobile element borders

# Here are the SV VCF.ALT types that are converted and what is output by this program:
# VCF.ALT	==>		CHR		START				END			...more fields
# =======			===		=====				===	
# <DEL>		==>		VCF.CHROM	VCF.POS-1			VCF.INFO.END		ID	QUAL	REF	ALT	FILTER	INFO	FORMAT	SAMPLE	...rest of VCF file samples
# <INS>		==>		VCF.CHROM	VCF.POS-1			VCF.POS			...
# <INDEL>	==>		VCF.CHROM	VCF.POS-1			VCF.INFO.END		...
# <DUP:TANDEM>	==>		VCF.CHROM	VCF.POS-1			VCF.INFO.END		...
# <DUP:INS>	==>		VCF.CHROM	VCF.POS-1			VCF.INFO.END		...
# <INV>		==>		VCF.CHROM	VCF.POS-1			VCF.INFO.END		...
# <INS:ME:...	==>		VCF.CHROM	VCF.POS+VCF.INFO.MEINFI.START	VCF.POS+VCF.INFO.MEINFI.END

# Example file input:
# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  AABAF.sorted.dupmarked.bam
# 1       106761073       gridss10_5537   A       <INDEL> 1363    .       SVTYPE=INDEL;END=106761089;SVLEN=64;TRANCHE=HIGH;INSSEQ=TGCCATGTTTTTGCATTTTTGTGCTTTTTGTTTGTAATTATGCTGTTTAAAATGACACCCAAGCACAGTGCTGAAGTGCT        GT      1/.
# 6       94381988        gridss118_7208  T       <DUP:TANDEM>    1100    .       SVTYPE=DUP;END=94382102;SVLEN=114;TRANCHE=HIGH  GT      1/.
# 6       124259133       gridss121_3513  T       <INS>   1064    .       SVTYPE=INS;END=124259133;SVLEN=35;TRANCHE=HIGH;INSSEQ=GAATCTGATGCACACAGCTTCATGAGAAAACAAGC       GT      1/.
# 7       102213180       gridss137_8268  C       <DUP:INS>       492.93  LOW_QUAL;SINGLE_ASSEMBLY        SVTYPE=DUP;END=115930270;SVLEN=13717100;TRANCHE=LOW;INSSEQ=TTCCCCCCCC   GT      1/.
# 7       157568893       gridss15_7528   A       <DEL>   176.63  LOW_QUAL;SINGLE_ASSEMBLY        SVTYPE=DEL;END=157568929;SVLEN=-36;TRANCHE=LOW  GT      1/.
# 8       10000381        notgridss100_100        G       <INV>   1697    ASSEMBLY_ONLY   SVTYPE=INV;END=10000481;TRANCHE=HIGH    GT      1/.
# 9       564965  .       A       <INS:ME:ALU>,<INS:ME:L1>        14      .       AC=1;AF=0.500;AN=2;IMPRECISE;MEINFO=ALU,-137,138,.;SVTYPE=INS;TARGETSITEDUPL=unknown;set=variant412     GT:MEI3MB:MEI3PR:MEI

# Example output file:
# 1	106761072	106761089	gridss10_5537	A	<INDEL>	1363	.	SVTYPE=INDEL;END=106761089;SVLEN=64;TRANCHE=HIGH;INSSEQ=TGCCATGTTTTTGCATTTTTGTGCTTTTTGTTTGTAATTATGCTGTTTAAAATGACACCCAAGCACAGTGCTGAAGTGCT
# 6	94381987	94382102	gridss118_7208	T	<DUP:TANDEM>	1100	.	SVTYPE=DUP;END=94382102;SVLEN=114;TRANCHE=HIGH
# 6	124259132	124259133	gridss121_3513	T	<INS>	1064	.	SVTYPE=INS;END=124259133;SVLEN=35;TRANCHE=HIGH;INSSEQ=GAATCTGATGCACACAGCTTCATGAGAAAACAAGC
# 7	102213179	115930270	gridss137_8268	C	<DUP:INS>	492.93	LOW_QUAL;SINGLE_ASSEMBLY	SVTYPE=DUP;END=115930270;SVLEN=13717100;TRANCHE=LOW;INSSEQ=TTCCCCCCCC
# 7	157568892	157568929	gridss15_7528	A	<DEL>	176.63	LOW_QUAL;SINGLE_ASSEMBLY	SVTYPE=DEL;END=157568929;SVLEN=-36;TRANCHE=LOW
# 8	10000380	10000481	notgridss100_100	G	<INV>	1697	ASSEMBLY_ONLY	SVTYPE=INV;END=10000481;TRANCHE=HIGH



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
def format_bed_output( chrom, pos1, pos2, infields ):

	snpid = infields[2]
	ref = infields[3]
	alt = infields[4]
	outline = str(chrom) + "\t" + str(pos1) + "\t" + str(pos2) + "\t" + str(snpid) + "\t" + ref + "\t" + alt
	for i in range( 5, 8 ):
		outline = outline + "\t" + str(infields[i])
	outline = outline + "\n"

	return outline

######################################################
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in an uncompressed VCF file from stdin, convert certain structural variants (SVs) to BED format, then output BED file to stdout.')
	parser.add_argument('-f', action="store", dest="filter_length", required=False, help='filter to retain SVs whose INFO.SVLEN is greater than or equal to this filter length')
	parser.add_argument('-w', action="store", dest="window_extension", required=False, help='extend the window on either side of the variant borders by this amount')
	args = parser.parse_args()
	filter_length = -1
	if (args.filter_length is not None):
		filter_length = int(args.filter_length)
	window_extension = 0
	if (args.window_extension is not None):
		window_extension = int(args.window_extension)

	# Read in the input VCF file from STDIN

	ALTs_to_process = [ '<DEL>', '<INS>', '<INDEL>', '<DUP:TANDEM>', '<DUP:INS>', '<INV>' ]

	in_header = True
	for inline in sys.stdin:
		inline = inline.strip()

		if (in_header == True):
			if (len(inline) >= 1):
				first_char = inline[0:1]
				if (first_char != '#'):
					in_header = False

		if (in_header == False): # We are processing VCF data records. We are no longer in the header part of the file.

			infields = inline.split("\t")
			chrom = str(infields[0])
			pos = int(infields[1])
			snpid = infields[2]
			ref = infields[3]
			alt = infields[4]
			#qual = infields[5]
			#snpfilter = infields[6]
			info = infields[7]

			this_is_MEI_ALT = False
			if (len(alt) >= 8):
				if (alt[0:8] == '<INS:ME:'):
					this_is_MEI_ALT = True

			process_this_ALT = True
			# if ((alt in ALTs_to_process) or this_is_MEI_ALT): # Process and output only those SVs we're interested in
			# 	process_this_ALT = True

			if (process_this_ALT):

				#if (snpid == '.'):
				#	sys.stderr.write( 'This variant will not be properly processed by the annotation pipeline because it does not have a VCF.ID: ' + chrom + ' ' + str(pos) + ' ' + snpid + ' ' + ref + ' ' + alt + "\n" )

				we_know_the_SV_length = False
				this_end = pos + 1
				this_svlen = 0
				this_meinfo_exists = False
				this_meinfo_start = 0
				this_meinfo_end = 0
				info_fields = info.split(';')
				for this_info_field in info_fields:
					bits = this_info_field.split('=')
					if (bits[0] == 'END'):
						this_end = int(bits[1])
					elif (bits[0] == 'SVLEN'):
						this_svlen = abs(int(bits[1]))
						we_know_the_SV_length = True
					elif (bits[0] == 'MEINFO'):
						bits2 = bits[1].split(',')
						meinfo_start = int(bits2[1])
						meinfo_end = int(bits2[2])
						this_svlen = meinfo_end - meinfo_start + 1
						this_meinfo_exists = True
						
				if (alt == '<INV>'): # SVLEN should be there for the other ALTs and according to VCF spec is not there for INV.
					this_svlen = this_end - pos + 1
					we_know_the_SV_length = True

				# Process and output only those SVs of length greater than the minimum.
				# Process and output all SVs if there is no filter specified.
				# Even if there is a filter, all SVTYPE=BND records will be processed because we don't know the SVLEN

				filter_out_this_variant = False
				if (args.filter_length is not None):
					if (we_know_the_SV_length):
						if (abs(this_svlen) < filter_length):
							filter_out_this_variant = True

				if (filter_out_this_variant == False):

					out_pos1 = pos - 1
					out_pos2 = this_end

					if (alt == '<DEL>'):
						out_pos1 = pos - 1
						out_pos2 = this_end
					elif (alt == '<INS>'):
						out_pos1 = pos - 1
						out_pos2 = pos
					elif (alt == '<INDEL>'):
						out_pos1 = pos - 1
						out_pos2 = this_end
					elif (alt == '<DUP:TANDEM>'):
						out_pos1 = pos - 1
						out_pos2 = this_end
					elif (alt == '<DUP:INS>'):
						out_pos1 = pos - 1
						out_pos2 = this_end
					elif (alt == '<INV>'):
						out_pos1 = pos - 1
						out_pos2 = this_end
					elif (this_is_MEI_ALT):
						out_pos1 = pos - 1
						out_pos2 = pos
						# MEINFO should always exist, but GATK merge sometimes removes it
						if (this_meinfo_exists):
							out_pos1 = pos + meinfo_start - 1
							out_pos2 = pos + meinfo_end
						# If there are samples, then get their actual start and end positions, and use the largest values for the annotation of all samples
						if (len(infields) >= 10):
							format_fields_string = infields[8]
							format_fields = format_fields_string.split(':')
							idx1 = -1
							idx2 = -1
							idx3 = -1
							idx4 = -1
							for i in range( 0, len(format_fields) ):
								if (format_fields[i] == 'MEI5MB'): # Mobster border5
									idx1 = i
								elif (format_fields[i] == 'MEI3MB'): # Mobster border3
									idx2 = i
								elif (format_fields[i] == 'MEI5RB'): # Mobelwrap refined border5
									idx3 = i
								elif (format_fields[i] == 'MEI3RB'): # Mobelwrap refined border3
									idx4 = i
							for sample_upto in range( 9, len(infields) ):
								sample_fields_string = infields[sample_upto]
								sample_fields = sample_fields_string.split(':')
								if (len(sample_fields) > 1): # if this sample has more than just the genotype field or just the empty field
									if ((is_integer(sample_fields[idx1])) and (idx1 > -1)):
										out_pos1 = min( out_pos1, int(sample_fields[idx1]) - 1 )
										out_pos2 = max( out_pos2, int(sample_fields[idx1]) )
									if ((is_integer(sample_fields[idx2])) and (idx2 > -1)):
										out_pos1 = min( out_pos1, int(sample_fields[idx2]) - 1 )
										out_pos2 = max( out_pos2, int(sample_fields[idx2]) )
									if ((is_integer(sample_fields[idx3])) and (idx3 > -1)):
										out_pos1 = min( out_pos1, int(sample_fields[idx3]) - 1 )
										out_pos2 = max( out_pos2, int(sample_fields[idx3]) )
									if ((is_integer(sample_fields[idx4])) and (idx4 > -1)):
										out_pos1 = min( out_pos1, int(sample_fields[idx4]) - 1 )
										out_pos2 = max( out_pos2, int(sample_fields[idx4]) )

					if (out_pos1 > out_pos2): # mobile elements can have the end position before the start position, bed file cannot, so swap the two positions for bed file
						temp_pos1 = out_pos1
						out_pos1 = out_pos2 - 1		# pos1 is to be zero-based
						out_pos2 = temp_pos1 + 1	# pos2 is to be 1-based

					out_pos1 = out_pos1 - window_extension
					out_pos2 = out_pos2 + window_extension

					if (out_pos1 < 0):
						out_pos1 = 0
					if (out_pos2 < 1):
						out_pos2 = 1

					outline = format_bed_output( chrom, out_pos1, out_pos2, infields )
					sys.stdout.write( outline )


if __name__=='__main__':
    main()


