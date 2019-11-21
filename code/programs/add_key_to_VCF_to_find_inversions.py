#!/usr/bin/python
# python add_key_to_VCF_to_find_inversions.py -i input_VCF -o output_file
# python add_key_to_VCF_to_find_inversions.py -i AOCS_92_32_tumour_sorted_markedDupl.sv.vcf -o AOCS_92_32_tumour_GRIDSS_with_key.txt

# This program is part of a 3-step pipeline to identify structural variant inversions in VCF output from GRIDSS.
#
# python add_key_to_VCF_to_find_inversions.py -i AOCS_92_32_tumour_sorted_markedDupl.sv.vcf -o AOCS_92_32_tumour_GRIDSS_with_key.txt
# sort -k1,1 -k2,2n -k3,3 -k4,4n AOCS_92_32_tumour_GRIDSS_with_key.txt > AOCS_92_32_tumour_GRIDSS_with_key_sorted.txt
# python identify_inversions_in_VCF_with_key.py -i AOCS_92_32_tumour_GRIDSS_with_key_sorted.txt -o AOCS_92_32_tumour_GRIDSS_inversions_only.vcf -u AOCS_92_32_tumour_GRIDSS_unused_breakend_records.vcf -gap 200

# The sort keys added to the VCF records in this pipeline are
# CHROM, lowest(POS,ALT), INFO.EVENT, POS

# After adding key and sorting, a VCF record looks like:
# #       1       #       1       ##fileformat=VCFv4.2
# #       2       #       2       ##ALT=<ID=INV,Description="Inversion">
# #       3       #       3       ##FILTER=<ID=ASSEMBLY_ONLY,Description="Variant is supported only by assembly evidence.">
# #       4       #       4       ##FILTER=<ID=ASSEMBLY_TOO_FEW_READ,Description="Not enough reads contribute to this assembly as specified by 'assembly.minReads'">
# #       162     #       162     #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  AOCS_92_32_tumour_sorted_markedDupl.bam
# 1       9999    gridss0_17704   9999    1       9999    gridss0_17704o  N       ]1:10472]N      538.49  NO_ASSEMBLY     AS=0;ASQ=0.00;ASRP=0;ASSR=0;BA=0;BAQ=0.00;BQ=
# 1       9999    gridss0_17704   10472   1       10472   gridss0_17704h  G       G[1:9999[       538.49  NO_ASSEMBLY     AS=0;ASQ=0.00;ASRP=0;ASSR=0;BA=2;BAQ=186.82;B
# 1       10116   gridss0_17706   10116   1       10116   gridss0_17706o  A       ]1:10317]A      682.77  SINGLE_ASSEMBLY AS=0;ASQ=0.00;ASRP=22;ASSR=15;BA=1;BAQ=59.92;
# 1       10116   gridss0_17706   10317   1       10317   gridss0_17706h  A       A[1:10116[      682.77  SINGLE_ASSEMBLY AS=1;ASQ=641.21;ASRP=22;ASSR=15;BA=1;BAQ=83.9
# 1       10181   gridss0_5422    10181   1       10181   gridss0_5422o   A       ACTTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCCTA[1:10242[  438.90  LOW_Q
# 1       10181   gridss0_5422    10242   1       10242   gridss0_5422h   A       ]1:10181]CTTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCCTAA  438.90  LOW_Q

# A structural variant inversion appears in GRIDSS output VCF as two pairs of breakends.
# Each pair has the same INFO.EVENT within the pair, but the 2 pairs do not have the same INFO.EVENT.
# Thus GRIDSS outputs an inversion as 2 separate pairs of events.

# <INV>			:	classic inversion (classsic inversion involves the canonical orientation of square brackets in the ALT field)
# <INV-DEL>		:	classic inversion with gap at either end
# <INV-INS>		:	classic inversion with more than one nucleotide in any ALT field
# <INV-OVERLAP>		:	inversion with overlap at either end
# <INV-DEL-INS-OVERLAP>	:	if multiple conditions found, they appear in this order
#
# Thus the full list of possible output ALT fields is:
# <INV>   <INV-DEL>   <INV-INS>   <INV-OVERLAP>   <INV-DEL-INS-OVERLAP>   <INV-DEL-INS>   <INV-DEL-OVERLAP>   <INV-INS-OVERLAP>


# Classic inversion:
# This will be output as <INV>
# 1       177614384       gridss17_6741o  C       C]1:177614422]  1964.05 SINGLE_ASSEMBLY AS=0;ASQ=0.00;ASRP=0;ASSR=52;BA=0;BAQ=0.00;BEID=asm17-33647;BQ=0.00;BSC=0;BSCQ=0.00;BUM=0;BUMQ=0.00;CAS=0;CASQ=0.00;CIPOS=0,5;CIR
# 1       177614422       gridss17_6741h  A       A]1:177614384]  1964.05 SINGLE_ASSEMBLY AS=1;ASQ=1167.83;ASRP=0;ASSR=52;BA=0;BAQ=0.00;BEID=asm17-33647;BQ=425.98;BSC=21;BSCQ=425.98;BUM=0;BUMQ=0.00;CAS=0;CASQ=0.00;CIPOS
# 1       177614385       gridss17_40377o T       [1:177614423[T  2300.11 SINGLE_ASSEMBLY AS=1;ASQ=1323.36;ASRP=0;ASSR=59;BA=0;BAQ=0.00;BEID=asm17-77885;BQ=414.88;BSC=21;BSCQ=414.88;BUM=0;BUMQ=0.00;CAS=0;CASQ=0.00;CIPOS
# 1       177614423       gridss17_40377h C       [1:177614385[C  2278.02 SINGLE_ASSEMBLY AS=0;ASQ=0.00;ASRP=0;ASSR=59;BA=0;BAQ=0.00;BEID=asm17-77885;BQ=0.00;BSC=0;BSCQ=0.00;BUM=0;BUMQ=0.00;CAS=0;CASQ=0.00;CIPOS=-5,0;CI
#
#          +--------------------------------------+
#          |                                      |
# ________ |          _______________________     |
# ________384    +-385_______________________422--+
#                |
#                |                                 ________________
#                +------------------------------423________________


# Inversion with loss of 1 nucleotide: 197756788
# This will be output as <INV-DEL>
# 1       197756787       gridss19_7574o  A       A]1:197757986]  790.77  .       AS=1;ASQ=252.02;ASRP=18;ASSR=14;BA=0;BAQ=0.00;BEID=asm19-35424,asm19-35425;BQ=20.33;BSC=1;BSCQ=20.33;BUM=0;BUMQ=0.00;CAS=0;CASQ=0.00;CIPO
# 1       197757986       gridss19_7574h  A       A]1:197756787]  790.77  .       AS=1;ASQ=307.05;ASRP=18;ASSR=14;BA=0;BAQ=0.00;BEID=asm19-35424,asm19-35425;BQ=103.17;BSC=5;BSCQ=103.17;BUM=0;BUMQ=0.00;CAS=0;CASQ=0.00;CI
# 1       197756789       gridss19_43156o T       [1:197757987[T  685.73  .       AS=1;ASQ=241.37;ASRP=22;ASSR=7;BA=0;BAQ=0.00;BEID=asm19-81069,asm19-81070,asm19-81120;BQ=19.19;BSC=1;BSCQ=19.19;BUM=0;BUMQ=0.00;CAS=0;CAS
# 1       197757987       gridss19_43156h C       [1:197756789[C  685.73  .       AS=2;ASQ=222.18;ASRP=22;ASSR=7;BA=0;BAQ=0.00;BEID=asm19-81069,asm19-81070,asm19-81120;BQ=0.00;BSC=0;BSCQ=0.00;BUM=0;BUMQ=0.00;CAS=0;CASQ=


# Inversion with loss of 1 nucleotide: 27374700
# This will be output as <INV-DEL>
# 21      27374158        gridss281_555o  A       A]21:27374699]  5483.21 .       AS=1;ASQ=1952.81;ASRP=134;ASSR=89;BA=0;BAQ=0.00;BEID=asm281-14050,asm281-14051;BQ=278.15;BSC=14;BSCQ=278.15;BUM=0;BUMQ=0.00;CAS=0;CASQ=0.
# 21      27374699        gridss281_555h  G       G]21:27374158]  5483.21 .       AS=1;ASQ=1855.74;ASRP=134;ASSR=89;BA=0;BAQ=0.00;BEID=asm281-14050,asm281-14051;BQ=200.43;BSC=9;BSCQ=181.08;BUM=1;BUMQ=19.35;CAS=0;CASQ=0.
# 21      27374159        gridss281_10694o        T       [21:27374701[T  5068.22 .       AS=1;ASQ=1726.94;ASRP=127;ASSR=79;BA=0;BAQ=0.00;BEID=asm281-49489,asm281-49504;BQ=201.97;BSC=10;BSCQ=201.97;BUM=0;BUMQ=0.00;CAS=0
# 21      27374701        gridss281_10694h        G       [21:27374159[G  5068.22 .       AS=1;ASQ=1802.40;ASRP=127;ASSR=79;BA=0;BAQ=0.00;BEID=asm281-49489,asm281-49504;BQ=282.87;BSC=13;BSCQ=263.52;BUM=1;BUMQ=19.35;CAS=


# Annotated by GRIDSS similar to an inversion, but may be an insertion where part of inserted sequence just happens to be same an inversion of reference sequence close by
# This will be annotated as <INV-OVERLAP>
# 11      23968142        gridss184_645o  A       A]11:23968183]  1251.39 SINGLE_ASSEMBLY AS=0;ASQ=0.00;ASRP=5;ASSR=31;BA=0;BAQ=0.00;BEID=asm184-6981,asm184-6983;BQ=281.27;BSC=14;BSCQ=281.27;BUM=0;BUMQ=0.00;CAS=0;CASQ=0
# 11      23968183        gridss184_645h  A       A]11:23968142]  1251.39 SINGLE_ASSEMBLY AS=2;ASQ=746.45;ASRP=5;ASSR=31;BA=1;BAQ=812.74;BEID=asm184-6981,asm184-6983;BQ=1054.26;BSC=12;BSCQ=241.52;BUM=0;BUMQ=0.00;CAS=0;C
# 11      23968147        gridss184_13844o        G       [11:23968180[G  2126.38 .       AS=2;ASQ=736.85;ASRP=6;ASSR=66;BA=0;BAQ=0.00;BEID=asm184-35738,asm184-35739,asm184-35751;BQ=160.05;BSC=8;BSCQ=160.05;BUM=0;BUMQ=0
# 11      23968180        gridss184_13844h        T       [11:23968147[T  2126.38 .       AS=1;ASQ=812.74;ASRP=6;ASSR=66;BA=0;BAQ=0.00;BEID=asm184-35738,asm184-35739,asm184-35751;BQ=280.56;BSC=14;BSCQ=280.56;BUM=0;BUMQ=
#
#          +--------------------------------------------+
#          |                                            |
# ________ |          _____________________________     |
# ________142    +-147_____________________________183--+
#                |
#                |                                 ________________
#                +------------------------------180________________


# Annotated by GRIDSS similar to an inversion, but may be an insertion where part of inserted sequence just happens to be same as inversion of reference sequence close by
# The will be annotated as <INV-INS-OVERLAP>
# 10      47023102        gridss172_1697o G       GTGGTGTGTGTGGGTGTGTGGTTGTGTGATGTGTGTATGCAATGTATGTGTGTGCATGTATGTGTGGAGTGTATGATTATGCATGTGTGGTGGT]10:47059590]
# 10      47059590        gridss172_1697h A       AACCACGTGCATCACACAACCACACACACATACACAACACACACATCACATAATCACACACAGTCACACATATACACATCATACAAACAACCCA]10:47023102]
# 10      47023103        gridss172_11565o        G       [10:47059581[CTGTGTGTACTGTATGGGGTGTTAGTATGTGTGTCTGTATGTGATGTGTGTGGTG
# 10      47059581        gridss172_11565h        C       [10:47023103[ACCCCCCACACACCAACCACACATCACAGCACACCACACACACACCACATACACA
#
#          +-------------------GTGTGTGT---------------------------+
#          |                                                      |
# ________ |          _______________________________________     |
# ________102    +-103_______________________________________590--+
#                |
#                |                                 ________________
#                +-----------ACACACAC-----------581________________


# Annotated by GRIDDS as similar to an inversion with breakends the opposite of how they present for inversions usually, with one breakend fragment also completely overlapping the inverted fragment, thus inverted fragment is a duplicate
# This will not be output by this pipeline. It will not be called as an INVERSION.
# 14      46521653        gridss224_12214o        T       [14:46521692[T  1952.65 .       AS=1;ASQ=757.12;ASRP=2;ASSR=61;BA=1;BAQ=738.93;BEID=asm224-55539,asm224-55540;BQ=1005.21;BSC=13;BSCQ=266.27;BUM=0;BUMQ=0.00;CAS=0
# 14      46521692        gridss224_12214h        C       [14:46521653[C  1952.65 .       AS=1;ASQ=676.87;ASRP=2;ASSR=61;BA=0;BAQ=0.00;BEID=asm224-55539,asm224-55540;BQ=158.22;BSC=8;BSCQ=158.22;BUM=0;BUMQ=0.00;CAS=0;CAS
# 14      46521655        gridss224_1761o T       T]14:46521688]  1996.21 .       AS=1;ASQ=738.93;ASRP=14;ASSR=59;BA=0;BAQ=0.00;BEID=asm224-27352,asm224-27353,asm224-27354;BQ=244.85;BSC=12;BSCQ=244.85;BUM=0;BUMQ=0.00;CA
# 14      46521688        gridss224_1761h A       A]14:46521655]  1996.21 .       AS=2;ASQ=763.20;ASRP=14;ASSR=59;BA=1;BAQ=676.87;BEID=asm224-27352,asm224-27353,asm224-27354;BQ=945.99;BSC=13;BSCQ=269.12;BUM=0;BUMQ=0.00;
# 
#                                        ___________________
#            +-------------------------92___________________
#            |
#            |    __
#            +--53__55--------+
#                             |
# ________________________    |
# ________________________88--+



__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2017, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import datetime
import math
import random
import commands
import argparse
import re
# import subprocess
# from multiprocessing import Pool

######################################################
def is_integer(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

######################################################
def is_float(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

######################################################
def extract_alt_pos( this_alt ):

	this_alt_chrom = ''
	this_alt_pos = -1
	this_alt = this_alt.replace( ']', '[' )
	bits = this_alt.split('[')
	for this_bit in bits:
		find_idx = this_bit.find( ':' )
		if (find_idx > -1):
			bits2 = this_bit.split(':')
			this_alt_chrom = str(bits2[0])
			this_alt_pos = int(bits2[1])

	return this_alt_chrom, this_alt_pos

######################################################
def extract_event( this_info ):

	this_event = ' '
	info_fields = this_info.split(';')
	for info_pair in info_fields:
		key_and_value = info_pair.split('=')
		# this_key = key_and_value[0]
		if (key_and_value[0] == 'EVENT'):
			# this_value = key_and_value[1]
			this_event = key_and_value[1]

	return this_event

######################################################
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Add sort keys in front of every VCF records.')
	parser.add_argument('-i', action="store", dest="input_VCF", required=True, help='Input VCF')
	parser.add_argument('-o', action="store", dest="output_file", required=True, help='Output file')
	args = parser.parse_args()

	output_file = open(args.output_file, 'w')

	# Open the input VCF and read its first record

	input_VCF = open(args.input_VCF)
	VCF_EOF = False
	continue_processing = True
	in_header = True
	inline = input_VCF.readline()
	inline = inline.strip()
	if (inline):
		continue_processing = True
	else:
		VCF_EOF = True
		continue_processing = False

	# Process each VCF record

	hdr_count = 0
	if (continue_processing):
		while (continue_processing):
			if (in_header == True):
				if (len(inline) >= 1):
					first_char = inline[0:1]
					if (first_char != '#'):
						in_header = False

			if (in_header):
				hdr_count = hdr_count + 1
				key = '#' + "\t" + str(hdr_count) + "\t" + '#' + "\t" + str(hdr_count)
				outline = key + "\t" + inline + "\n"
				output_file.write( outline )

			else: # (in_header == False)

				# Process this VCF record
				# Only keep it if the POS and ALT are on the same chromosome. If not, then it can't be an inversion.

				infields = inline.split("\t")
				this_chrom = str(infields[0])
				this_pos = int(infields[1])
				this_alt = infields[4]
				this_info = infields[7]

				this_alt_chrom, this_alt_pos = extract_alt_pos( this_alt )
				if (this_chrom == this_alt_chrom):

					lowest_pos = min( this_pos, this_alt_pos )
					this_event = extract_event( this_info )

					key = this_chrom + "\t" + str(lowest_pos) + "\t" + this_event + "\t" + str(this_pos)
					outline = key + "\t" + inline + "\n"
					output_file.write( outline )

			# Now read the next VCF record

			inline = input_VCF.readline()
			inline = inline.strip()
			if (inline):
				continue_processing = True
			else:
				VCF_EOF = True
				continue_processing = False

	output_file.close()

if __name__=='__main__':
    main()


