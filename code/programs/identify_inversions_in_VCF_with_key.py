#!/usr/bin/python
# python identify_inversions_in_VCF_with_key.py -i input_file -o output_file [-gap gap_or_overlap] [-circular look_for_circular_inversions_too] [-u unused_breakend_records]
# python identify_inversions_in_VCF_with_key.py -i AOCS_92_32_tumour_GRIDSS_with_key_sorted.txt -o AOCS_92_32_tumour_GRIDSS_inversions_only.txt



# This program is part of a 3-step pipeline to identify structural variant inversions in VCF output from GRIDSS.
#
# python add_key_to_VCF_to_find_inversions.py -i AOCS_92_32_tumour_sorted_markedDupl.sv.vcf -o AOCS_92_32_tumour_GRIDSS_with_key.txt
# sort -k1,1 -k2,2n -k3,3 -k4,4n AOCS_92_32_tumour_GRIDSS_with_key.txt > AOCS_92_32_tumour_GRIDSS_with_key_sorted.txt
# python identify_inversions_in_VCF_with_key.py -i AOCS_92_32_tumour_GRIDSS_with_key_sorted.txt -o AOCS_92_32_tumour_GRIDSS_inversions_only.vcf -u AOCS_92_32_tumour_GRIDSS_unused_breakend_records.vcf -gap 200



# If circular inversions are to be identified too, then there are more steps to the process, and the circular-inversions should be identified separately and after identifing inversions
# because identifying circular inversions first may cause normal inversions to be misidentified as a circular inversions.
#
# python add_key_to_VCF_to_find_inversions.py -i AOCS_92_32_tumour_sorted_markedDupl.sv.vcf -o AOCS_92_32_tumour_GRIDSS_with_key.txt
# sort -k1,1 -k2,2n -k3,3 -k4,4n AOCS_92_32_tumour_GRIDSS_with_key.txt > AOCS_92_32_tumour_GRIDSS_with_key_sorted.txt
# python identify_inversions_in_VCF_with_key.py -i AOCS_92_32_tumour_GRIDSS_with_key_sorted.txt -o AOCS_92_32_tumour_GRIDSS_inversions_only.vcf -u AOCS_92_32_tumour_GRIDSS_unused_breakend_records.vcf -gap 200
#
# python add_key_to_VCF_to_find_inversions.py -i AOCS_92_32_tumour_GRIDSS_unused_breakend_records.vcf -o AOCS_92_32_tumour_GRIDSS_unused_breakend_records.txt
# sort -k1,1 -k2,2n -k3,3 -k4,4n AOCS_92_32_tumour_GRIDSS_unused_breakend_records.txt > AOCS_92_32_tumour_GRIDSS_unused_breakend_records_with_key_sorted.txt
# python identify_inversions_in_VCF_with_key.py -i AOCS_92_32_tumour_GRIDSS_unused_breakend_records_with_key_sorted.txt -o AOCS_92_32_tumour_GRIDSS_circular_inversions_too.vcf -circular CIRCULAR -gap 300
#
# grep '^#' AOCS_92_32_tumour_sorted_markedDupl.sv.vcf > header_only.vcf
# grep -v '^#' AOCS_92_32_tumour_GRIDSS_inversions_only.vcf > inversions_only.txt
# grep -v '^#' AOCS_92_32_tumour_GRIDSS_circular_inversions_too.vcf > circular_inversions_too.txt
# cat inversions_only.vcf circular_inversions_too.vcf | sort -k1,1 -k2,2n > inversions_sorted.txt
# cat header_only.vcf inversions_sorted.txt > AOCS_92_32_tumour_GRIDSS_inversions_and_circular_inversions.vcf



# When this program produces one output VCF record for a group of input VCF records,
# it assigns a TRANCHE "score" of LOW, INTERMEDIATE, or HIGH corresponding to the confidence in this higher-level structural variant call.



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
# 1       177614384       gridss17_6741   177614384       1       177614384       gridss17_6741o  C       C]1:177614422]  1964.05 SINGLE_ASSEMBLY AS=0;ASQ=0.00;ASRP=0;
# 1       177614384       gridss17_6741   177614422       1       177614422       gridss17_6741h  A       A]1:177614384]  1964.05 SINGLE_ASSEMBLY AS=1;ASQ=1167.83;ASRP
# 1       177614385       gridss17_40377  177614385       1       177614385       gridss17_40377o T       [1:177614423[T  2300.11 SINGLE_ASSEMBLY AS=1;ASQ=1323.36;ASRP
# 1       177614385       gridss17_40377  177614423       1       177614423       gridss17_40377h C       [1:177614385[C  2278.02 SINGLE_ASSEMBLY AS=0;ASQ=0.00;ASRP=0;

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
# This will be output as <INV> ==> 1	177614385	gridss17_40377	T	<INV>	1964.05	SINGLE_ASSEMBLY	SVTYPE=INV;END=177614422;TRANCHE=INTERMEDIATE
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
# This will be output as <INV-DEL> ==> 1	197756789	gridss19_43156	T	<INV-DEL>	790.77	.	SVTYPE=INV;END=197757986;TRANCHE=INTERMEDIATE;INVBND1DEL=1
# 1       197756787       gridss19_7574o  A       A]1:197757986]  790.77  .       AS=1;ASQ=252.02;ASRP=18;ASSR=14;BA=0;BAQ=0.00;BEID=asm19-35424,asm19-35425;BQ=20.33;BSC=1;BSCQ=20.33;BUM=0;BUMQ=0.00;CAS=0;CASQ=0.00;CIPO
# 1       197757986       gridss19_7574h  A       A]1:197756787]  790.77  .       AS=1;ASQ=307.05;ASRP=18;ASSR=14;BA=0;BAQ=0.00;BEID=asm19-35424,asm19-35425;BQ=103.17;BSC=5;BSCQ=103.17;BUM=0;BUMQ=0.00;CAS=0;CASQ=0.00;CI
# 1       197756789       gridss19_43156o T       [1:197757987[T  685.73  .       AS=1;ASQ=241.37;ASRP=22;ASSR=7;BA=0;BAQ=0.00;BEID=asm19-81069,asm19-81070,asm19-81120;BQ=19.19;BSC=1;BSCQ=19.19;BUM=0;BUMQ=0.00;CAS=0;CAS
# 1       197757987       gridss19_43156h C       [1:197756789[C  685.73  .       AS=2;ASQ=222.18;ASRP=22;ASSR=7;BA=0;BAQ=0.00;BEID=asm19-81069,asm19-81070,asm19-81120;BQ=0.00;BSC=0;BSCQ=0.00;BUM=0;BUMQ=0.00;CAS=0;CASQ=


# Inversion with loss of 1 nucleotide: 27374700
# This will be output as <INV-DEL> ==> 21	27374159	gridss281_10694	T	<INV-DEL>	5483.21	.	SVTYPE=INV;END=27374699;TRANCHE=HIGH;INVBND2DEL=1
# 21      27374158        gridss281_555o  A       A]21:27374699]  5483.21 .       AS=1;ASQ=1952.81;ASRP=134;ASSR=89;BA=0;BAQ=0.00;BEID=asm281-14050,asm281-14051;BQ=278.15;BSC=14;BSCQ=278.15;BUM=0;BUMQ=0.00;CAS=0;CASQ=0.
# 21      27374699        gridss281_555h  G       G]21:27374158]  5483.21 .       AS=1;ASQ=1855.74;ASRP=134;ASSR=89;BA=0;BAQ=0.00;BEID=asm281-14050,asm281-14051;BQ=200.43;BSC=9;BSCQ=181.08;BUM=1;BUMQ=19.35;CAS=0;CASQ=0.
# 21      27374159        gridss281_10694o        T       [21:27374701[T  5068.22 .       AS=1;ASQ=1726.94;ASRP=127;ASSR=79;BA=0;BAQ=0.00;BEID=asm281-49489,asm281-49504;BQ=201.97;BSC=10;BSCQ=201.97;BUM=0;BUMQ=0.00;CAS=0
# 21      27374701        gridss281_10694h        G       [21:27374159[G  5068.22 .       AS=1;ASQ=1802.40;ASRP=127;ASSR=79;BA=0;BAQ=0.00;BEID=asm281-49489,asm281-49504;BQ=282.87;BSC=13;BSCQ=263.52;BUM=1;BUMQ=19.35;CAS=


# Annotated by GRIDSS similar to an inversion, but may be an insertion where part of inserted sequence just happens to be same an inversion of reference sequence close by
# This will be annotated as <INV-DEL-OVERLAP> ==> 11	23968147	gridss184_13844	G	<INV-DEL-OVERLAP>	1251.39	SINGLE_ASSEMBLY	SVTYPE=INV;END=23968183;TRANCHE=INTERMEDIATE;INVBND1DEL=4;INVBND2OVERLAP=4
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
# These will be annotated as <INV-INS-OVERLAP> ==> 10	47023103	gridss172_11565	G	<INV-INS-OVERLAP>	5266.49	.	SVTYPE=INV;END=47059590;TRANCHE=HIGH;INVBND1BDR5INS=GTGGTGT...;INVBND2BDR5INS=AACCACGT...;INVBND1BDR3INS=CTGTGTGT...;INVBND2BDR3INS=ACCCCCCACA...;INVBND2OVERLAP=10
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

# Annotated by GRIDDS as an inversion of a circular segment.
# We saw data where a true inversion was reported by gridss as a circular segment due to inaccuracy of the actual breakend positions
# This will be annotated as <INV-CIRCULAR>
# 1       203802400       gridss20_15877  203802400       1       203802400       gridss20_15877o G       [1:203803310[G  226.05  LOW_QUAL;SINGLE_ASSEMBLY        AS=0;ASQ=0.00;ASRP=8;ASSR=0;BA=0;BAQ=0.
# 1       203802400       gridss20_15877  203803310       1       203803310       gridss20_15877h T       [1:203802400[T  226.05  LOW_QUAL;SINGLE_ASSEMBLY        AS=1;ASQ=113.02;ASRP=8;ASSR=0;BA=0;BAQ=
# 1       203802416       gridss20_1048   203802416       1       203802416       gridss20_1048o  C       C]1:203803589]  625.77  .       AS=1;ASQ=287.36;ASRP=13;ASSR=19;BA=0;BAQ=0.00;BEID=asm20-9177,a
# 1       203802416       gridss20_1048   203803589       1       203803589       gridss20_1048h  A       A]1:203802416]  625.77  .       AS=2;ASQ=221.63;ASRP=13;ASSR=19;BA=1;BAQ=154.07;BEID=asm20-9177
#
#         +------------------------------+
#         |                              |
#         |      _____                   |      _____
#         +--2400_____2416--+            +--3310_____3589--+
#                           |                              |
#                           |                              |
#                           +------------------------------+
#
#
# How gridss should have reported this true inversion example:
# 1       203802415       gridss20_1048o  C       C]1:203803589]  625.77  .       AS=1;ASQ=287.36;ASRP=13;ASSR=19;BA=0;BAQ=0.00;BEID=asm20-9177,a
# 1       203803589       gridss20_1048h  A       A]1:203802415]  625.77  .       AS=2;ASQ=221.63;ASRP=13;ASSR=19;BA=1;BAQ=154.07;BEID=asm20-9177
# 1       203802416       gridss20_15877o G       [1:203803590[G  226.05  LOW_QUAL;SINGLE_ASSEMBLY        AS=0;ASQ=0.00;ASRP=8;ASSR=0;BA=0;BAQ=0.
# 1       203803590       gridss20_15877h T       [1:203802416[T  226.05  LOW_QUAL;SINGLE_ASSEMBLY        AS=1;ASQ=113.02;ASRP=8;ASSR=0;BA=0;BAQ=
#
#          +--------------------------------------+
#          |                                      |
# _______  |           _____________________      |
# _______2415    +-2416_____________________3589--+
#                |
#                |                                 ________________
#                +-----------------------------3590________________

# This program does not handle or identify inverted-duplications
# __________>>>>>>>>>><<<<<<<<<<__________



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

	this_alt_pos = -1
	this_alt = this_alt.replace( ']', '[' )
	bits = this_alt.split('[')
	for this_bit in bits:
		find_idx = this_bit.find( ':' )
		if (find_idx > -1):
			bits2 = this_bit.split(':')
			this_alt_chrom = str(bits2[0])
			this_alt_pos = int(bits2[1])

	return this_alt_pos

######################################################
def extract_alt_seq( this_alt ):

	this_alt_seq = ''
	this_alt = this_alt.replace( ']', '[' )
	find_idx = this_alt.find( '[' )
	if (find_idx == 0): # [1:1234144[C
		find_idx = this_alt.find( '[', 1 )
		this_alt_seq = this_alt[ (find_idx+1): ]
	else: # T[1:1362556[
		this_alt_seq = this_alt[ 0: find_idx ]

	return this_alt_seq

######################################################
def extract_alt_bracket( this_alt ):

	this_alt_bracket = ''
	find_idx = this_alt.find( '[' )
	if (find_idx >= 0):
		this_alt_bracket = '['
	else:
		find_idx = this_alt.find( ']' )
		if (find_idx >= 0):
			this_alt_bracket = ']'

	return this_alt_bracket

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
def extract_AS_RAS_from_INFO( this_qual, this_info ):

	if (is_float(this_qual)):
		this_qual = float(this_qual)
	else:
		this_qual = float(0)
	this_AS = 0
	this_RAS = 0
	info_fields = this_info.split(';')
	for this_info_field in info_fields:
		bits = this_info_field.split('=')
		if (bits[0] == 'AS'):
			if (is_integer(str(bits[1]))):
				this_AS = int(bits[1])
		elif (bits[0] == 'RAS'):
			if (is_integer(str(bits[1]))):
				this_RAS = int(bits[1])

	return this_qual, this_AS, this_RAS

######################################################
def determine_tranche_for_one_record( this_qual, this_AS, this_RAS ):

	this_tranche = 'LOW'
	if (this_qual >= float(500)):
		if ((this_AS > 0) or (this_RAS > 0)):
			this_tranche = 'INTERMEDIATE'
	if (this_qual >= float(1000)):
		if ((this_AS > 0) and (this_RAS > 0)):
			this_tranche = 'HIGH'

	return this_tranche

######################################################
def determine_tranche_field( BND_records_fields ):

	# Determine the new TRANCHE field

	this_qual = BND_records_fields[0]['QUAL']
	this_filter = BND_records_fields[0]['FILTER']
	this_info = BND_records_fields[0]['INFO']
	this_qual, this_AS, this_RAS = extract_AS_RAS_from_INFO( this_qual, this_info )
	this_tranche = determine_tranche_for_one_record( this_qual, this_AS, this_RAS )

	highest_qual = this_qual
	highest_filter = this_filter
	highest_info = this_info
	highest_tranche = this_tranche

	for i in range( 1, len(BND_records_fields) ):

		this_qual = BND_records_fields[i]['QUAL']
		this_filter = BND_records_fields[i]['FILTER']
		this_info = BND_records_fields[i]['INFO']
		this_qual, this_AS, this_RAS = extract_AS_RAS_from_INFO( this_qual, this_info )
		this_tranche = determine_tranche_for_one_record( this_qual, this_AS, this_RAS )

		if (this_tranche == 'HIGH'):
			highest_filter = this_filter
			highest_tranche = this_tranche

		if ((this_tranche == 'INTERMEDIATE') and (highest_tranche != 'HIGH')):
			highest_filter = this_filter
			highest_tranche = this_tranche

		if (this_qual > highest_qual):
			highest_qual = this_qual

	return highest_tranche, highest_qual, highest_filter

######################################################
def max_tranche2( tranche2_in1, tranche2_in2, tranche2_in3, tranche2_in4 ):

	tranche2_final = ""
	if ((tranche2_in1 == "HIGH") or (tranche2_in2 == "HIGH") or (tranche2_in3 == "HIGH") or (tranche2_in4 == "HIGH")):
		tranche2_final = "HIGH"
	elif ((tranche2_in1 == "INTERMEDIATE") or (tranche2_in2 == "INTERMEDIATE") or (tranche2_in3 == "INTERMEDIATE") or (tranche2_in4 == "INTERMEDIATE")):
		tranche2_final = "INTERMEDIATE"
	elif ((tranche2_in1 == "LOW") or (tranche2_in2 == "LOW") or (tranche2_in3 == "LOW") or (tranche2_in4 == "LOW")):
		tranche2_final = "LOW"

	return tranche2_final

######################################################
def determine_VAFs( vaf_in1, vaf_in2, vaf_in3, vaf_in4 ):

	vaf_in1 = float(vaf_in1)
	vaf_in2 = float(vaf_in2)
	vaf_in3 = float(vaf_in3)
	vaf_in4 = float(vaf_in4)
	vaf_final = str(max( vaf_in1, vaf_in2, vaf_in3, vaf_in4 ))
	list_of_vafs = str(vaf_in1) + ',' + str(vaf_in2) + ',' + str(vaf_in3) + ',' + str(vaf_in4)

	return vaf_final, list_of_vafs

######################################################
def determine_GT( vaf_in1, vaf_in2, vaf_in3, vaf_in4 ):

	new_GT = '0/1'
	vaf_in1 = float(vaf_in1)
	vaf_in2 = float(vaf_in2)
	vaf_in3 = float(vaf_in3)
	vaf_in4 = float(vaf_in4)
	if ((vaf_in1 > 0.65) or (vaf_in2 > 0.65) or (vaf_in3 > 0.65) or (vaf_in4 > 0.65)):
		new_GT = '1/1'

	return new_GT

######################################################
def init_rec():

	record = {}
	record['inline'] = ''
	record['chrom'] = ''
	record['lowest_pos'] = -1
	record['event'] = ''
	record['pos'] = -1
	record['alt_pos'] = -1
	record['alt_seq'] = ''
	record['alt_bracket'] = ''

	return record

######################################################
def fill_record_from_inline( inline ):

	record = {}
	record['inline'] = inline
	infields = inline.split("\t")
	record['chrom'] = str(infields[4])
	record['lowest_pos'] = int(infields[1])
	record['event'] = str(infields[2])
	record['pos'] = int(infields[5])
	alt_field = infields[8]
	record['alt_pos'] = extract_alt_pos(alt_field)
	record['alt_seq'] = extract_alt_seq(alt_field)
	record['alt_bracket'] = extract_alt_bracket(alt_field)

	return record

######################################################
def same_breakend(record1, record2):

	is_same_breakend = False
	if ((record1['chrom'] == record2['chrom']) and (record1['lowest_pos'] == record2['lowest_pos']) and (record1['event'] == record2['event'])):
		is_same_breakend = True

	return is_same_breakend

######################################################
def pos_is_inbetween_other_2_pos( middle_pos, pos1, pos2 ):

	is_inbetween = False
	if ((pos1 <= middle_pos) and (middle_pos <= pos2)):
		is_inbetween = True

	return is_inbetween

######################################################
def looks_like_INVERSION_batch(gap_or_overlap, look_for_circular_inversions_too, pair1_rec1, pair1_rec2, pair2_rec1, pair2_rec2):

	is_inversion_batch = False
	if ( (pair1_rec1['chrom'] == pair1_rec2['chrom']) and (pair1_rec1['chrom'] == pair2_rec1['chrom']) and (pair1_rec1['chrom'] == pair2_rec2['chrom']) ):
		if ( (same_breakend(pair1_rec1, pair1_rec2)) and (same_breakend(pair2_rec1, pair2_rec2)) ):
			if ( (pair1_rec1['alt_bracket'] == ']') and (pair1_rec2['alt_bracket'] == ']') and (pair2_rec1['alt_bracket'] == '[') and (pair2_rec2['alt_bracket'] == '[') ):
				# brackets are arranged as for inversion arrangement
				if ( pos_is_inbetween_other_2_pos(pair2_rec1['pos'],pair1_rec1['pos'],pair1_rec2['pos']) or pos_is_inbetween_other_2_pos(pair2_rec2['pos'],pair1_rec1['pos'],pair1_rec2['pos']) ): 
					# one or both of pair2 ends are between the pair1 ends

					diff1 = abs(pair1_rec1['pos'] - pair2_rec1['pos'])
					diff2 = abs(pair1_rec2['pos'] - pair2_rec2['pos'])
					if ((diff1 <= gap_or_overlap) and (diff2 <= gap_or_overlap)):
						is_inversion_batch = True

			elif (look_for_circular_inversions_too):
				# brackets are arranged as for circular inversion arrangement
				if ( pos_is_inbetween_other_2_pos(pair2_rec1['pos'],pair1_rec1['pos'],pair1_rec2['pos']) or pos_is_inbetween_other_2_pos(pair2_rec2['pos'],pair1_rec1['pos'],pair1_rec2['pos']) ): 
					# one or both of pair2 ends are between the pair1 ends
					if ( (pair1_rec1['alt_bracket'] == '[') and (pair1_rec2['alt_bracket'] == '[') and (pair2_rec1['alt_bracket'] == ']') and (pair2_rec2['alt_bracket'] == ']') ):

						diff1 = abs(pair1_rec1['pos'] - pair2_rec1['pos'])
						diff2 = abs(pair1_rec2['pos'] - pair2_rec2['pos'])
						if ((diff1 <= gap_or_overlap) and (diff2 <= gap_or_overlap)):
							is_inversion_batch = True

	return is_inversion_batch

######################################################
def assign_inversion_type( pair1_rec1, pair1_rec2, pair2_rec1, pair2_rec2 ):

	inversion_type = 'INV'
	INVBND1BDR5INS = ''
	INVBND2BDR5INS = ''
	INVBND1BDR3INS = ''
	INVBND2BDR3INS = ''
	INVBND1DEL = 0
	INVBND2DEL = 0
	INVBND1OVERLAP = 0
	INVBND2OVERLAP = 0
	pair1_rec1_pos = int(pair1_rec1['pos'])
	pair1_rec2_pos = int(pair1_rec2['pos'])
	pair2_rec1_pos = int(pair2_rec1['pos'])
	pair2_rec2_pos = int(pair2_rec2['pos'])

	add_del = False
	if ( (pair1_rec1_pos+1) < pair2_rec1_pos ):
		add_del = True
		INVBND1DEL = pair2_rec1_pos - pair1_rec1_pos - 1
	if ( (pair1_rec2_pos+1) < pair2_rec2_pos ):
		add_del = True
		INVBND2DEL = pair2_rec2_pos - pair1_rec2_pos - 1
	if (add_del):
		inversion_type = inversion_type + '-DEL'

	pair1_rec1_alt_seq = pair1_rec1['alt_seq']
	pair1_rec2_alt_seq = pair1_rec2['alt_seq']
	pair2_rec1_alt_seq = pair2_rec1['alt_seq']
	pair2_rec2_alt_seq = pair2_rec2['alt_seq']
	if ( (len(pair1_rec1_alt_seq)>1) or (len(pair1_rec2_alt_seq)>1) or (len(pair2_rec1_alt_seq)>1) or (len(pair2_rec2_alt_seq)>1) ):
		inversion_type = inversion_type + '-INS'
		if (len(pair1_rec1_alt_seq) > 1):
			INVBND1BDR5INS = pair1_rec1_alt_seq
		if (len(pair1_rec2_alt_seq) > 1):
			INVBND2BDR5INS = pair1_rec2_alt_seq
		if (len(pair2_rec1_alt_seq) > 1):
			INVBND1BDR3INS = pair2_rec1_alt_seq
		if (len(pair2_rec2_alt_seq) > 1):
			INVBND2BDR3INS = pair2_rec2_alt_seq

	add_overlap = False
	if ( (pair1_rec1_pos+1) > pair2_rec1_pos ):
		add_overlap = True
		INVBND1OVERLAP = pair1_rec1_pos - pair2_rec1_pos + 1
	if ( (pair1_rec2_pos+1) > pair2_rec2_pos ):
		add_overlap = True
		INVBND2OVERLAP = pair1_rec2_pos - pair2_rec2_pos + 1
	if (add_overlap):
		inversion_type = inversion_type + '-OVERLAP'

	pair1_rec1_bracket = pair1_rec1['alt_bracket']
	pair1_rec2_bracket = pair1_rec2['alt_bracket']
	pair2_rec1_bracket = pair2_rec1['alt_bracket']
	pair2_rec2_bracket = pair2_rec2['alt_bracket']
	if ( (pair1_rec1_bracket == '[') and (pair1_rec2_bracket == '[') and (pair2_rec1_bracket == ']') and (pair2_rec2_bracket == ']') ):
		inversion_type = inversion_type + '-CIRCULAR'

	inversion_type = '<' + inversion_type + '>'

	return inversion_type, INVBND1BDR5INS, INVBND2BDR5INS, INVBND1BDR3INS, INVBND2BDR3INS, INVBND1DEL, INVBND2DEL, INVBND1OVERLAP, INVBND2OVERLAP

######################################################
def choose_highest_filter_result( filter1, filter2, filter3, filter4 ):

	return_filter = '.'
	highest_seen = 0 # .=0 LOW=1 MEDIUM=2 HIGH=3
	filters = []
	filters.append( filter1 )
	filters.append( filter2 )
	filters.append( filter3 )
	filters.append( filter4 )
	for this_filter in filters:
		if (highest_seen < 3):
			if (this_filter.find('HIGH_QUAL') > -1):
				highest_seen = 3
				return_filter = this_filter
		if (highest_seen < 2):
			if (this_filter.find('MEDIUM_QUAL') > -1):
				highest_seen = 2
				return_filter = this_filter
		if (highest_seen < 1):
			if (this_filter.find('LOW') > -1):
				highest_seen = 1
				return_filter = this_filter
		if (highest_seen == 0):
			if ( (this_filter != '.') and (return_filter == '.') ):
				return_filter = this_filter

	return return_filter

######################################################
def extract_info_fields( in_fields_string ):

	return_fields = {}
	delim = ';'
	split_fields = in_fields_string.split(delim)
	for key_and_value in split_fields:
		split_key_value = key_and_value.split('=')
		this_key = split_key_value[0]
		if (len(split_key_value) > 1):
			this_value = split_key_value[1]
		else:
			this_value = ''
		return_fields[this_key] = this_value

	return return_fields

######################################################
def choose_highest_field( value1, value2 ):

	return_info_field = value1
	if ( (is_integer(value1)) and (is_integer(value2)) ):
		return_info_field = max( int(value1), int(value2) )
	elif ( (is_float(value1)) and (is_float(value2)) ):
		return_info_field = max( float(value1), float(value2) )

	return return_info_field

######################################################
def choose_highest_info_fields( info1, info2, info3, info4 ):

	return_info_fields = {}
	info1_fields = extract_info_fields( info1 )
	info2_fields = extract_info_fields( info2 )
	info3_fields = extract_info_fields( info3 )
	info4_fields = extract_info_fields( info4 )

	for info_key in info2_fields:
		if info_key not in return_info_fields:
			return_info_fields[info_key] = ''
		if info_key not in info2_fields:
			info2_fields[info_key] = ''
		return_info_fields[info_key] = info2_fields[info_key]
	for info_key in info1_fields:
		if info_key not in return_info_fields:
			return_info_fields[info_key] = ''
		if info_key not in info1_fields:
			info1_fields[info_key] = ''
		return_info_fields[info_key] = choose_highest_field( return_info_fields[info_key], info1_fields[info_key] )
	for info_key in info3_fields:
		if info_key not in return_info_fields:
			return_info_fields[info_key] = ''
		if info_key not in info3_fields:
			info3_fields[info_key] = ''
		return_info_fields[info_key] = choose_highest_field( return_info_fields[info_key], info3_fields[info_key] )
	for info_key in info4_fields:
		if info_key not in return_info_fields:
			return_info_fields[info_key] = ''
		if info_key not in info4_fields:
			info4_fields[info_key] = ''
		return_info_fields[info_key] = choose_highest_field( return_info_fields[info_key], info4_fields[info_key] )

	return_info = ''
	for info_key in return_info_fields:
		key_and_value = str(info_key)
		if (str(return_info_fields[info_key]) != ''):
			key_and_value = str(info_key) + '=' + str(return_info_fields[info_key])
		if (return_info == ''):
			return_info = key_and_value
		else:
			return_info = return_info + ';' + key_and_value
	if (return_info == ''):
		return_info = '.'

	return return_info

######################################################
def extract_sample_fields( in_format_fields, in_sample_fields ):

	return_fields = {}
	delim = ':'

	format_fields = []
	split_fields = in_format_fields.split(delim)
	for this_field in split_fields:
		format_fields.append( this_field )

	sample_fields = []
	split_fields = in_sample_fields.split(delim)
	i = 0
	for this_field in split_fields:
		return_fields[ format_fields[i] ] = this_field
		i = i + 1

	return return_fields

######################################################
def merge_sample_breakends_for_sample_field( inv_format, sample1, sample2, sample3, sample4 ):

	inv_format_fields = inv_format.split(':')

	return_sample_fields = {}
	sample1_fields = extract_sample_fields( inv_format, sample1 )
	sample2_fields = extract_sample_fields( inv_format, sample2 )
	sample3_fields = extract_sample_fields( inv_format, sample3 )
	sample4_fields = extract_sample_fields( inv_format, sample4 )

	for sample_key in sample2_fields:
		return_sample_fields[sample_key] = sample2_fields[sample_key]
	for sample_key in sample1_fields:
		return_sample_fields[sample_key] = choose_highest_field( return_sample_fields[sample_key], sample1_fields[sample_key] )
	for sample_key in sample3_fields:
		return_sample_fields[sample_key] = choose_highest_field( return_sample_fields[sample_key], sample3_fields[sample_key] )
	for sample_key in sample4_fields:
		return_sample_fields[sample_key] = choose_highest_field( return_sample_fields[sample_key], sample4_fields[sample_key] )

	return_sample_fields['GT'] = '1/.'

	return_sample = ''
	for sample_key in inv_format_fields:
		if (return_sample == ''):
			return_sample = str(return_sample_fields[sample_key])
		else:
			return_sample = return_sample + ':' + str(return_sample_fields[sample_key])
	if (return_sample == ''):
		return_sample = '.'

	return return_sample

######################################################
def extract_fields_from_INFO( this_info ):

	this_bndvaf = float(0)
	this_ciend = ','
	this_cipos = ','
	this_cirpos = ','
	this_tranche2 = ''
	this_ihompos = ''
	this_as = ''
	this_asq = ''
	this_asrp = ''
	this_asrr = ''
	this_ba = ''
	this_baq = ''
	this_beid = ''
	this_bq = ''
	this_bsc = ''
	this_bscq = ''
	this_bum = ''
	this_bumq = ''
	this_cas = ''
	this_casq = ''
	this_cq = ''
	this_homlen = ''
	this_homseq = ''
	this_ic = ''
	this_ihompos = ''
	this_iq = ''
	this_ras = ''
	this_rasq = ''
	this_ref_count = ''
	this_refpair = ''
	this_rp = ''
	this_rpq = ''
	this_rsi = ''
	this_si = ''
	this_sr = ''
	this_srq = ''

	info_fields = this_info.split(';')
	for this_info_field in info_fields:
		bits = this_info_field.split('=')
		if (bits[0] == 'BNDVAF'):
			this_bndvaf = float(bits[1])
		elif (bits[0] == 'CIEND'):
			this_ciend = bits[1]
		elif (bits[0] == 'CIPOS'):
			this_cipos = bits[1]
		elif (bits[0] == 'CIRPOS'):
			this_cirpos = bits[1]
		elif (bits[0] == 'TRANCHE2'):
			this_tranche2 = bits[1]
		elif (bits[0] == 'IHOMPOS'):
			this_ihompos = bits[1]
		elif (bits[0] == 'AS'):
			this_as = int(bits[1])
		elif (bits[0] == 'ASQ'):
			this_asq = float(bits[1])
		elif (bits[0] == 'ASRP'):
			this_asrp = int(bits[1])
		elif (bits[0] == 'ASRR'):
			this_asrr = int(bits[1])
		elif (bits[0] == 'BA'):
			this_ba = int(bits[1])
		elif (bits[0] == 'BAQ'):
			this_baq = float(bits[1])
		elif (bits[0] == 'BEID'):
			this_beid = bits[1]
		elif (bits[0] == 'BQ'):
			this_bq = float(bits[1])
		elif (bits[0] == 'BSC'):
			this_bsc = int(bits[1])
		elif (bits[0] == 'BSCQ'):
			this_bscq = float(bits[1])
		elif (bits[0] == 'BUM'):
			this_bum = int(bits[1])
		elif (bits[0] == 'BUMQ'):
			this_bumq = float(bits[1])
		elif (bits[0] == 'CAS'):
			this_cas = int(bits[1])
		elif (bits[0] == 'CASQ'):
			this_casq = float(bits[1])
		elif (bits[0] == 'CQ'):
			this_cq = float(bits[1])
		elif (bits[0] == 'HOMLEN'):
			this_homlen = int(bits[1])
		elif (bits[0] == 'HOMSEQ'):
			this_homseq = bits[1]
		elif (bits[0] == 'IC'):
			this_ic = int(bits[1])
		elif (bits[0] == 'HOMPOS'):
			this_hompos = int(bits[1])
		elif (bits[0] == 'IQ'):
			this_iq = float(bits[1])
		elif (bits[0] == 'RAS'):
			this_ras = int(bits[1])
		elif (bits[0] == 'RASQ'):
			this_rasq = float(bits[1])
		elif (bits[0] == 'REF'):
			this_ref_count = int(bits[1])
		elif (bits[0] == 'REFPAIR'):
			this_refpair = int(bits[1])
		elif (bits[0] == 'RP'):
			this_rp = int(bits[1])
		elif (bits[0] == 'RPQ'):
			this_rpq = float(bits[1])
		elif (bits[0] == 'RSI'):
			this_rsi = int(bits[1])
		elif (bits[0] == 'SI'):
			this_si = int(bits[1])
		elif (bits[0] == 'SR'):
			this_sr = int(bits[1])
		elif (bits[0] == 'SRQ'):
			this_srq = float(bits[1])

	return this_bndvaf, this_ciend, this_cipos, this_cirpos, this_tranche2, this_ihompos, this_as, this_asq, this_asrp, this_asrr, this_ba, this_baq, this_beid, this_bq, this_bsc, this_bscq, this_bum, this_bumq, this_cas, this_casq, this_cq, this_homlen, this_homseq, this_ic, this_ihompos, this_iq, this_ras, this_rasq, this_ref_count, this_refpair, this_rp, this_rpq, this_rsi, this_si, this_sr, this_srq

######################################################
def determine_cipos_ciend( BND_records_fields ):

	# BND_records_fields[0] = pair1_rec1_fields
	# BND_records_fields[1] = pair1_rec2_fields
	# BND_records_fields[2] = pair2_rec1_fields
	# BND_records_fields[3] = pair2_rec2_fields

	this_cipos = BND_records_fields[0]['CIPOS']
	this_ciend = BND_records_fields[2]['CIPOS']

	return this_cipos, this_ciend

######################################################
def parse_BND_record(BND_record):

	this_BND_record = {}
	this_BND_record['CHROM'] = BND_record[4]
	this_BND_record['POS'] = BND_record[5]
	this_BND_record['ID'] = BND_record[6]
	this_BND_record['REF'] = BND_record[7]
	this_BND_record['ALT'] = BND_record[8]
	this_BND_record['QUAL'] = BND_record[9]
	this_BND_record['FILTER'] = BND_record[10]
	this_BND_record['INFO'] = BND_record[11]
	this_BND_record['FORMAT'] = BND_record[12]
	this_bndvaf, this_ciend, this_cipos, this_cirpos, this_tranche2, this_ihompos, this_as, this_asq, this_asrp, this_asrr, this_ba, this_baq, this_beid, this_bq, this_bsc, this_bscq, this_bum, this_bumq, this_cas, this_casq, this_cq, this_homlen, this_homseq, this_ic, this_ihompos, this_iq, this_ras, this_rasq, this_ref_count, this_refpair, this_rp, this_rpq, this_rsi, this_si, this_sr, this_srq = extract_fields_from_INFO( this_BND_record['INFO'] )
	this_BND_record['BNDVAF'] = this_bndvaf
	this_BND_record['CIEND'] = this_ciend
	this_BND_record['CIPOS'] = this_cipos
	this_BND_record['CIRPOS'] = this_cirpos
	this_BND_record['TRANCHE2'] = this_tranche2
	this_BND_record['IHOMPOS'] = this_ihompos
	this_BND_record['AS'] = this_as
	this_BND_record['ASQ'] = this_asq
	this_BND_record['ASRP'] = this_asrp
	this_BND_record['ASRR'] = this_asrr
	this_BND_record['BA'] = this_ba
	this_BND_record['BAQ'] = this_baq
	this_BND_record['BEID'] = this_beid
	this_BND_record['BQ'] = this_bq
	this_BND_record['BSC'] = this_bsc
	this_BND_record['BSCQ'] = this_bscq
	this_BND_record['BUM'] = this_bum
	this_BND_record['BUMQ'] = this_bumq
	this_BND_record['CAS'] = this_cas
	this_BND_record['CASQ'] = this_casq
	this_BND_record['CQ'] = this_cq
	this_BND_record['HOMLEN'] = this_homlen
	this_BND_record['HOMSEQ'] = this_homseq
	this_BND_record['IC'] = this_ic
	this_BND_record['IHOMPOS'] = this_ihompos
	this_BND_record['IQ'] = this_iq
	this_BND_record['RAS'] = this_ras
	this_BND_record['RASQ'] = this_rasq
	this_BND_record['REFCOUNT'] = this_ref_count
	this_BND_record['REFPAIR'] = this_refpair
	this_BND_record['RP'] = this_rp
	this_BND_record['RPQ'] = this_rpq
	this_BND_record['RSI'] = this_rsi
	this_BND_record['SI'] = this_si
	this_BND_record['SR'] = this_sr
	this_BND_record['SRQ'] = this_srq

	return this_BND_record

######################################################
def extract_BND_values( BND_records_fields ):

	new_gridssAS = str(BND_records_fields[0]['AS']) + ',' + str(BND_records_fields[1]['AS']) + ',' + str(BND_records_fields[2]['AS']) + ',' + str(BND_records_fields[3]['AS'])
	new_gridssASQ = str(BND_records_fields[0]['ASQ']) + ',' + str(BND_records_fields[1]['ASQ']) + ',' + str(BND_records_fields[2]['ASQ']) + ',' + str(BND_records_fields[3]['ASQ'])
	new_gridssASRP = str(BND_records_fields[0]['ASRP']) + ',' + str(BND_records_fields[1]['ASRP']) + ',' + str(BND_records_fields[2]['ASRP']) + ',' + str(BND_records_fields[3]['ASRP'])
	new_gridssASRR = str(BND_records_fields[0]['ASRR']) + ',' + str(BND_records_fields[1]['ASRR']) + ',' + str(BND_records_fields[2]['ASRR']) + ',' + str(BND_records_fields[3]['ASRR'])
	new_gridssBA = str(BND_records_fields[0]['BA']) + ',' + str(BND_records_fields[1]['BA']) + ',' + str(BND_records_fields[2]['BA']) + ',' + str(BND_records_fields[3]['BA'])
	new_gridssBAQ = str(BND_records_fields[0]['BAQ']) + ',' + str(BND_records_fields[1]['BAQ']) + ',' + str(BND_records_fields[2]['BAQ']) + ',' + str(BND_records_fields[3]['BAQ'])
	new_gridssBEID = str(BND_records_fields[0]['BEID']) + ',' + str(BND_records_fields[1]['BEID']) + ',' + str(BND_records_fields[2]['BEID']) + ',' + str(BND_records_fields[3]['BEID'])
	new_gridssBQ = str(BND_records_fields[0]['BQ']) + ',' + str(BND_records_fields[1]['BQ']) + ',' + str(BND_records_fields[2]['BQ']) + ',' + str(BND_records_fields[3]['BQ'])
	new_gridssBSC = str(BND_records_fields[0]['BSC']) + ',' + str(BND_records_fields[1]['BSC']) + ',' + str(BND_records_fields[2]['BSC']) + ',' + str(BND_records_fields[3]['BSC'])
	new_gridssBSCQ = str(BND_records_fields[0]['BSCQ']) + ',' + str(BND_records_fields[1]['BSCQ']) + ',' + str(BND_records_fields[2]['BSCQ']) + ',' + str(BND_records_fields[3]['BSCQ'])
	new_gridssBUM = str(BND_records_fields[0]['BUM']) + ',' + str(BND_records_fields[1]['BUM']) + ',' + str(BND_records_fields[2]['BUM']) + ',' + str(BND_records_fields[3]['BUM'])
	new_gridssBUMQ = str(BND_records_fields[0]['BUMQ']) + ',' + str(BND_records_fields[1]['BUMQ']) + ',' + str(BND_records_fields[2]['BUMQ']) + ',' + str(BND_records_fields[3]['BUMQ'])
	new_gridssCAS = str(BND_records_fields[0]['CAS']) + ',' + str(BND_records_fields[1]['CAS']) + ',' + str(BND_records_fields[2]['CAS']) + ',' + str(BND_records_fields[3]['CAS'])
	new_gridssCASQ = str(BND_records_fields[0]['CASQ']) + ',' + str(BND_records_fields[1]['CASQ']) + ',' + str(BND_records_fields[2]['CASQ']) + ',' + str(BND_records_fields[3]['CASQ'])
	new_gridssCIEND = str(BND_records_fields[0]['CIEND']) + ',' + str(BND_records_fields[1]['CIEND']) + ',' + str(BND_records_fields[2]['CIEND']) + ',' + str(BND_records_fields[3]['CIEND'])
	new_gridssCIPOS = str(BND_records_fields[0]['CIPOS']) + ',' + str(BND_records_fields[1]['CIPOS']) + ',' + str(BND_records_fields[2]['CIPOS']) + ',' + str(BND_records_fields[3]['CIPOS'])
	new_gridssCIRPOS = str(BND_records_fields[0]['CIRPOS']) + ',' + str(BND_records_fields[1]['CIRPOS']) + ',' + str(BND_records_fields[2]['CIRPOS']) + ',' + str(BND_records_fields[3]['CIRPOS'])
	new_gridssCQ = str(BND_records_fields[0]['CQ']) + ',' + str(BND_records_fields[1]['CQ']) + ',' + str(BND_records_fields[2]['CQ']) + ',' + str(BND_records_fields[3]['CQ'])
	new_gridssFILTER = str(BND_records_fields[0]['FILTER']) + ',' + str(BND_records_fields[1]['FILTER']) + ',' + str(BND_records_fields[2]['FILTER']) + ',' + str(BND_records_fields[3]['FILTER'])
	new_gridssHOMLEN = str(BND_records_fields[0]['HOMLEN']) + ',' + str(BND_records_fields[1]['HOMLEN']) + ',' + str(BND_records_fields[2]['HOMLEN']) + ',' + str(BND_records_fields[3]['HOMLEN'])
	new_gridssHOMSEQ = str(BND_records_fields[0]['HOMSEQ']) + ',' + str(BND_records_fields[1]['HOMSEQ']) + ',' + str(BND_records_fields[2]['HOMSEQ']) + ',' + str(BND_records_fields[3]['HOMSEQ'])
	new_gridssIC = str(BND_records_fields[0]['IC']) + ',' + str(BND_records_fields[1]['IC']) + ',' + str(BND_records_fields[2]['IC']) + ',' + str(BND_records_fields[3]['IC'])
	new_gridssIHOMPOS = str(BND_records_fields[0]['IHOMPOS']) + ',' + str(BND_records_fields[1]['IHOMPOS']) + ',' + str(BND_records_fields[2]['IHOMPOS']) + ',' + str(BND_records_fields[3]['IHOMPOS'])
	new_gridssIQ = str(BND_records_fields[0]['IQ']) + ',' + str(BND_records_fields[1]['IQ']) + ',' + str(BND_records_fields[2]['IQ']) + ',' + str(BND_records_fields[3]['IQ'])
	new_gridssQUAL = str(BND_records_fields[0]['QUAL']) + ',' + str(BND_records_fields[1]['QUAL']) + ',' + str(BND_records_fields[2]['QUAL']) + ',' + str(BND_records_fields[3]['QUAL'])
	new_gridssRAS = str(BND_records_fields[0]['RAS']) + ',' + str(BND_records_fields[1]['RAS']) + ',' + str(BND_records_fields[2]['RAS']) + ',' + str(BND_records_fields[3]['RAS'])
	new_gridssRASQ = str(BND_records_fields[0]['RASQ']) + ',' + str(BND_records_fields[1]['RASQ']) + ',' + str(BND_records_fields[2]['RASQ']) + ',' + str(BND_records_fields[3]['RASQ'])
	new_gridssREF = str(BND_records_fields[0]['REFCOUNT']) + ',' + str(BND_records_fields[1]['REFCOUNT']) + ',' + str(BND_records_fields[2]['REFCOUNT']) + ',' + str(BND_records_fields[3]['REFCOUNT'])
	new_gridssREFPAIR = str(BND_records_fields[0]['REFPAIR']) + ',' + str(BND_records_fields[1]['REFPAIR']) + ',' + str(BND_records_fields[2]['REFPAIR']) + ',' + str(BND_records_fields[3]['REFPAIR'])
	new_gridssRP = str(BND_records_fields[0]['RP']) + ',' + str(BND_records_fields[1]['RP']) + ',' + str(BND_records_fields[2]['RP']) + ',' + str(BND_records_fields[3]['RP'])
	new_gridssRPQ = str(BND_records_fields[0]['RPQ']) + ',' + str(BND_records_fields[1]['RPQ']) + ',' + str(BND_records_fields[2]['RPQ']) + ',' + str(BND_records_fields[3]['RPQ'])
	new_gridssRSI = str(BND_records_fields[0]['RSI']) + ',' + str(BND_records_fields[1]['RSI']) + ',' + str(BND_records_fields[2]['RSI']) + ',' + str(BND_records_fields[3]['RSI'])
	new_gridssSI = str(BND_records_fields[0]['SI']) + ',' + str(BND_records_fields[1]['SI']) + ',' + str(BND_records_fields[2]['SI']) + ',' + str(BND_records_fields[3]['SI'])
	new_gridssSR = str(BND_records_fields[0]['SR']) + ',' + str(BND_records_fields[1]['SR']) + ',' + str(BND_records_fields[2]['SR']) + ',' + str(BND_records_fields[3]['SR'])
	new_gridssSRQ = str(BND_records_fields[0]['SRQ']) + ',' + str(BND_records_fields[1]['SRQ']) + ',' + str(BND_records_fields[2]['SRQ']) + ',' + str(BND_records_fields[3]['SRQ'])

	return new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ

######################################################
def build_INVERSION_record( pair1_rec1, pair1_rec2, pair2_rec1, pair2_rec2, output_gridss_info_fields, output_gridss_format_fields ):

	pair1_rec1_fields = pair1_rec1['inline'].split("\t")
	pair1_rec2_fields = pair1_rec2['inline'].split("\t")
	pair2_rec1_fields = pair2_rec1['inline'].split("\t")
	pair2_rec2_fields = pair2_rec2['inline'].split("\t")
	BND_records_fields = []
	BND_records_fields.append( parse_BND_record(pair1_rec1_fields) )
	BND_records_fields.append( parse_BND_record(pair1_rec2_fields) )
	BND_records_fields.append( parse_BND_record(pair2_rec1_fields) )
	BND_records_fields.append( parse_BND_record(pair2_rec2_fields) )

	inversion_chrom = str(pair2_rec1_fields[4])
	inversion_start_pos = int(pair2_rec1_fields[5])
	inversion_end_pos = int(pair1_rec2_fields[5])
	if (inversion_end_pos < inversion_start_pos):
		print 'Warning: SVTYPE=INV inversion on chromosome', inversion_chrom, 'has END position', str(inversion_end_pos), 'that is less than start POS', str(inversion_start_pos)

	chrom = inversion_chrom
	pos = str(inversion_start_pos)
	inv_id = str(pair2_rec1_fields[2])
	ref = str(pair2_rec1_fields[7])
	alt, INVBND1BDR5INS, INVBND2BDR5INS, INVBND1BDR3INS, INVBND2BDR3INS, INVBND1DEL, INVBND2DEL, INVBND1OVERLAP, INVBND2OVERLAP = assign_inversion_type( pair1_rec1, pair1_rec2, pair2_rec1, pair2_rec2 )

	#new_qual = max( float(pair1_rec1_fields[9]), float(pair1_rec2_fields[9]), float(pair2_rec1_fields[9]), float(pair2_rec2_fields[9]) )
	#new_filter = choose_highest_filter_result( pair1_rec1_fields[10], pair1_rec2_fields[10], pair2_rec1_fields[10], pair2_rec2_fields[10] )
	#old_info = choose_highest_info_fields( pair1_rec1_fields[11], pair1_rec2_fields[11], pair2_rec1_fields[11], pair2_rec2_fields[11] )

	new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ = extract_BND_values( BND_records_fields )

	new_info_tranche, new_qual, new_filter = determine_tranche_field( BND_records_fields )
	new_cipos, new_ciend = determine_cipos_ciend( BND_records_fields )
	new_info_tranche2 = max_tranche2( BND_records_fields[0]['TRANCHE2'], BND_records_fields[1]['TRANCHE2'], BND_records_fields[2]['TRANCHE2'], BND_records_fields[3]['TRANCHE2'] )
	new_VAF, new_gridssVAF = determine_VAFs( BND_records_fields[0]['BNDVAF'], BND_records_fields[1]['BNDVAF'], BND_records_fields[2]['BNDVAF'], BND_records_fields[3]['BNDVAF'] )
	new_GT = determine_GT( BND_records_fields[0]['BNDVAF'], BND_records_fields[1]['BNDVAF'], BND_records_fields[2]['BNDVAF'], BND_records_fields[3]['BNDVAF'] )

	new_info = 'SVTYPE=INV;END=' + str(inversion_end_pos) + ';TRANCHE=' + str(new_info_tranche) + ';TRANCHE2=' + str(new_info_tranche2) + ';CIPOS=' + str(new_cipos) + ';CIEND=' + str(new_ciend) + ';VAF=' + str(new_VAF) + ';BNDVAFS=' + str(new_gridssVAF)

	if (INVBND1BDR5INS != ''):
		new_info = new_info + ';INVBND1BDR5INS=' + INVBND1BDR5INS
	if (INVBND2BDR5INS != ''):
		new_info = new_info + ';INVBND2BDR5INS=' + INVBND2BDR5INS
	if (INVBND1BDR3INS != ''):
		new_info = new_info + ';INVBND1BDR3INS=' + INVBND1BDR3INS
	if (INVBND2BDR3INS != ''):
		new_info = new_info + ';INVBND2BDR3INS=' + INVBND2BDR3INS
	if (INVBND1DEL > 0):
		new_info = new_info + ';INVBND1DEL=' + str(INVBND1DEL)
	if (INVBND2DEL > 0):
		new_info = new_info + ';INVBND2DEL=' + str(INVBND2DEL)
	if (INVBND1OVERLAP > 0):
		new_info = new_info + ';INVBND1OVERLAP=' + str(INVBND1OVERLAP)
	if (INVBND2OVERLAP > 0):
		new_info = new_info + ';INVBND2OVERLAP=' + str(INVBND2OVERLAP)

	if (output_gridss_info_fields):
		new_info = new_info + ';gridssVAF=' + str(new_gridssVAF) + ';gridssAS=' + str(new_gridssAS) + ';gridssASQ=' + str(new_gridssASQ) + ';gridssASRP=' + str(new_gridssASRP) + ';gridssASRR=' + str(new_gridssASRR)
		new_info = new_info + ';gridssBA=' + str(new_gridssBA) + ';gridssBAQ=' + str(new_gridssBAQ) + ';gridssBEID=' + str(new_gridssBEID) + ';gridssBQ=' + str(new_gridssBQ) + ';gridssBSC=' + str(new_gridssBSC) + ';gridssBSCQ=' + str(new_gridssBSCQ)
		new_info = new_info + ';gridssBUM=' + str(new_gridssBUM) + ';gridssBUMQ=' + str(new_gridssBUMQ) + ';gridssCAS=' + str(new_gridssCAS) + ';gridssCASQ=' + str(new_gridssCASQ)
		new_info = new_info + ';gridssCIEND=' + str(new_gridssCIEND) + ';gridssCIPOS=' + str(new_gridssCIPOS) + ';gridssCIRPOS=' + str(new_gridssCIRPOS) + ';gridssCQ=' + str(new_gridssCQ)
		new_info = new_info + ';gridssFILTER=' + str(new_gridssFILTER) + ';gridssHOMLEN=' + str(new_gridssHOMLEN) + ';gridssHOMSEQ=' + str(new_gridssHOMSEQ)
		new_info = new_info + ';gridssIC=' + str(new_gridssIC) + ';gridssIHOMPOS=' + str(new_gridssIHOMPOS) + ';gridssIQ=' + str(new_gridssIQ)
		new_info = new_info + ';gridssQUAL=' + str(new_gridssQUAL) + ';gridssRAS=' + str(new_gridssRAS) + ';gridssRASQ=' + str(new_gridssRASQ)
		new_info = new_info + ';gridssREF=' + str(new_gridssREF) + ';gridssREFPAIR=' + str(new_gridssREFPAIR) + ';gridssRP=' + str(new_gridssRP) + ';gridssRPQ=' + str(new_gridssRPQ)
		new_info = new_info + ';gridssRSI=' + str(new_gridssRSI) + ';gridssSI=' + str(new_gridssSI) + ';gridssSR=' + str(new_gridssSR) + ';gridssSRQ=' + str(new_gridssSRQ)
	info = new_info

	inv_format = 'GT:POS:END:TRANCHE:TRANCHE2:CIPOS:CIEND:VAF:BNDVAFS:INVBND1BDR5INS:INVBND2BDR5INS:INVBND1BDR3INS:INVBND2BDR3INS:INVBND1DEL:INVBND2DEL:INVBND1OVERLAP:INVBND2OVERLAP'

	if (output_gridss_format_fields):
		inv_format = inv_format + ':gridssVAF:gridssAS:gridssASQ:gridssASRP:gridssASRR:gridssBA:gridssBAQ:gridssBEID:gridssBQ:gridssBSC:gridssBSCQ'
		inv_format = inv_format + ':gridssBUM:gridssBUMQ:gridssCAS:gridssCASQ:gridssCIEND:gridssCIPOS:gridssCIRPOS:gridssCQ'
		inv_format = inv_format + ':gridssFILTER:gridssHOMLEN:gridssHOMSEQ:gridssIC:gridssIHOMPOS:gridssIQ:gridssQUAL:gridssRAS:gridssRASQ'
		inv_format = inv_format + ':gridssREF:gridssREFPAIR:gridssRP:gridssRPQ:gridssRSI:gridssSI:gridssSR:gridssSRQ'

	samples = []
	for i in range( 13, len(pair1_rec2_fields) ):
		#this_sample = merge_sample_breakends_for_sample_field( inv_format, pair1_rec1_fields[i], pair1_rec2_fields[i], pair2_rec1_fields[i], pair2_rec2_fields[i] )
		this_gt = new_GT # not '1/.'
		if (str(new_info_tranche) == ''):
			new_info_tranche = '.'
		if (INVBND1BDR5INS == ''):
			INVBND1BDR5INS = '.'
		if (INVBND2BDR5INS == ''):
			INVBND2BDR5INS = '.'
		if (INVBND1BDR3INS == ''):
			INVBND1BDR3INS = '.'
		if (INVBND2BDR3INS == ''):
			INVBND2BDR3INS = '.'
		if (INVBND1DEL == 0):
			INVBND1DEL = '.'
		if (INVBND2DEL == 0):
			INVBND2DEL = '.'
		if (INVBND1OVERLAP == 0):
			INVBND1OVERLAP = '.'
		if (INVBND2OVERLAP == 0):
			INVBND2OVERLAP = '.'

		this_sample = this_gt + ':' + str(pos) + ':' + str(inversion_end_pos) + ':' + str(new_info_tranche) + ':' + str(new_info_tranche2) + ':' + new_cipos + ':' + new_ciend + ':' + str(new_VAF) + ':' + str(new_gridssVAF)
		this_sample = this_sample + ':' + str(INVBND1BDR5INS) + ':' + str(INVBND2BDR5INS) + ':' + str(INVBND1BDR3INS) + ':' + str(INVBND2BDR3INS) + ':' + str(INVBND1DEL) + ':' + str(INVBND2DEL) + ':' + str(INVBND1OVERLAP) + ':' + str(INVBND2OVERLAP)

	if (output_gridss_format_fields):
		this_sample = this_sample + ':' + str(new_gridssVAF) + ':' + str(new_gridssAS) + ':' + str(new_gridssASQ) + ':' + str(new_gridssASRP) + ':' + str(new_gridssASRR) + ':' + str(new_gridssBA) + ':' + str(new_gridssBAQ) + ':' + str(new_gridssBEID) + ':' + str(new_gridssBQ) + ':' + str(new_gridssBSC) + ':' + str(new_gridssBSCQ)
		this_sample = this_sample + ':' + str(new_gridssBUM) + ':' + str(new_gridssBUMQ) + ':' + str(new_gridssCAS) + ':' + str(new_gridssCASQ) + ':' + str(new_gridssCIEND) + ':' + str(new_gridssCIPOS) + ':' + str(new_gridssCIRPOS) + ':' + str(new_gridssCQ)
		this_sample = this_sample + ':' + str(new_gridssFILTER) + ':' + str(new_gridssHOMLEN) + ':' + str(new_gridssHOMSEQ) + ':' + str(new_gridssIC) + ':' + str(new_gridssIHOMPOS) + ':' + str(new_gridssIQ) + ':' + str(new_gridssQUAL) + ':' + str(new_gridssRAS) + ':' + str(new_gridssRASQ)
		this_sample = this_sample + ':' + str(new_gridssREF) + ':' + str(new_gridssREFPAIR) + ':' + str(new_gridssRP) + ':' + str(new_gridssRPQ) + ':' + str(new_gridssRSI) + ':' + str(new_gridssSI) + ':' + str(new_gridssSR) + ':' + str(new_gridssSRQ)

	samples.append( this_sample )

	return_record = chrom + "\t" + str(pos) + "\t" + inv_id + "\t" + ref + "\t" + alt + "\t" + str(new_qual) + "\t" + str(new_filter) + "\t" + info + "\t" + inv_format
	for this_sample in samples:
		return_record = return_record + "\t" + this_sample

	return_record = return_record

	return return_record

######################################################
def output_info_headers_for_gridss_derived_fields( output_file, output_gridss_info_fields ):

	if (output_gridss_info_fields):
		output_file.write('##INFO=<ID=gridssVAF,Number=.,Type=String,Description="VAF values of the gridss-called BND calls used to call this SV">\n')
		output_file.write('##INFO=<ID=gridssAS,Number=.,Type=String,Description="AS values of the gridss-called BND calls used to call this SV. Count of assemblies supporting breakpoint">\n')
		output_file.write('##INFO=<ID=gridssASQ,Number=.,Type=String,Description="ASQ values of the gridss-called BND calls used to call this SV. Quality score of assemblies supporting breakpoint">\n')
		output_file.write('##INFO=<ID=gridssASRP,Number=.,Type=String,Description="ASRP values of the gridss-called BND calls used to call this SV. Count of read pairs incorporated into any breakpoint assembly">\n')
		output_file.write('##INFO=<ID=gridssASRR,Number=.,Type=String,Description="ASRR values of the gridss-called BND calls used to call this SV. Count of split, soft clipped or indel-containing reads incorporated into any breakpoint assemblies">\n')
		output_file.write('##INFO=<ID=gridssBA,Number=.,Type=String,Description="BA values of the gridss-called BND calls used to call this SV. Count of assemblies supporting just local breakend">\n')
		output_file.write('##INFO=<ID=gridssBAQ,Number=.,Type=String,Description="BAQ values of the gridss-called BND calls used to call this SV. Quality score of assemblies supporting just local breakend">\n')
		output_file.write('##INFO=<ID=gridssBEID,Number=.,Type=String,Description="BEID values of the gridss-called BND calls used to call this SV. Breakend assemblies contributing support to the breakpoint">\n')
		output_file.write('##INFO=<ID=gridssBQ,Number=.,Type=String,Description="BQ values of the gridss-called BND calls used to call this SV. Quality score of breakend evidence">\n')
		output_file.write('##INFO=<ID=gridssBSC,Number=.,Type=String,Description="BSC values of the gridss-called BND calls used to call this SV. Count of soft clips supporting just local breakend">\n')
		output_file.write('##INFO=<ID=gridssBSCQ,Number=.,Type=String,Description="BSCQ values of the gridss-called BND calls used to call this SV. Quality score of soft clips supporting just local breakend">\n')
		output_file.write('##INFO=<ID=gridssBUM,Number=.,Type=String,Description="BUM values of the gridss-called BND calls used to call this SV. Count of read pairs (with one read unmapped) supporting just local breakend">\n')
		output_file.write('##INFO=<ID=gridssBUMQ,Number=.,Type=String,Description="BUMQ values of the gridss-called BND calls used to call this SV. Quality score of read pairs (with one read unmapped) supporting just local breakend">\n')
		output_file.write('##INFO=<ID=gridssCAS,Number=.,Type=String,Description="CAS values of the gridss-called BND calls used to call this SV. Count of complex compound breakpoint assemblies supporting breakpoint from elsewhere">\n')
		output_file.write('##INFO=<ID=gridssCASQ,Number=.,Type=String,Description="CASQ values of the gridss-called BND calls used to call this SV. Quality score of complex compound breakpoint assemblies supporting breakpoint from elsewhere">\n')
		output_file.write('##INFO=<ID=gridssCIEND,Number=.,Type=String,Description="CIEND values of the gridss-called BND calls used to call this SV. Confidence interval around END for imprecise variants">\n')
		output_file.write('##INFO=<ID=gridssCIPOS,Number=.,Type=String,Description="CIPOS values of the gridss-called BND calls used to call this SV. Confidence interval around POS for imprecise variants">\n')
		output_file.write('##INFO=<ID=gridssCIRPOS,Number=.,Type=String,Description="CIRPOS values of the gridss-called BND calls used to call this SV. Confidence interval around remote breakend POS for imprecise variants">\n')
		output_file.write('##INFO=<ID=gridssCQ,Number=.,Type=String,Description="CQ values of the gridss-called BND calls used to call this SV. Breakpoint quality score before evidence reallocation">\n')
		output_file.write('##INFO=<ID=gridssFILTER,Number=.,Type=String,Description="FILTER values of the gridss-called BND calls used to call this SV">\n')
		output_file.write('##INFO=<ID=gridssHOMLEN,Number=.,Type=String,Description="HOMLEN values of the gridss-called BND calls used to call this SV. Length of base pair identical micro-homology at event breakpoints">\n')
		output_file.write('##INFO=<ID=gridssHOMSEQ,Number=.,Type=String,Description="HOMSEQ values of the gridss-called BND calls used to call this SV. Sequence of base pair identical micro-homology at event breakpoints">\n')
		output_file.write('##INFO=<ID=gridssIC,Number=.,Type=String,Description="IC values of the gridss-called BND calls used to call this SV. Count of read indels supporting breakpoint">\n')
		output_file.write('##INFO=<ID=gridssIHOMPOS,Number=.,Type=String,Description="IHOMPOS values of the gridss-called BND calls used to call this SV. Position of inexact homology">\n')
		output_file.write('##INFO=<ID=gridssIQ,Number=.,Type=String,Description="IQ values of the gridss-called BND calls used to call this SV. Quality score of read indels supporting breakpoint">\n')
		output_file.write('##INFO=<ID=gridssQUAL,Number=.,Type=String,Description="QUAL values of the gridss-called BND calls used to call this SV.">\n')
		output_file.write('##INFO=<ID=gridssRAS,Number=.,Type=String,Description="RAS values of the gridss-called BND calls used to call this SV. Count of assemblies supporting breakpoint from remote breakend">\n')
		output_file.write('##INFO=<ID=gridssRASQ,Number=.,Type=String,Description="RASQ values of the gridss-called BND calls used to call this SV. Quality score of assemblies supporting breakpoint from remote breakend">\n')
		output_file.write('##INFO=<ID=gridssREF,Number=.,Type=String,Description="REF values of the gridss-called BND calls used to call this SV. Count of reads mapping across this breakend">\n')
		output_file.write('##INFO=<ID=gridssREFPAIR,Number=.,Type=String,Description="REFPAIR values of the gridss-called BND calls used to call this SV. Count of reference read pairs spanning this breakpoint supporting the reference allele">\n')
		output_file.write('##INFO=<ID=gridssRP,Number=.,Type=String,Description="RP values of the gridss-called BND calls used to call this SV. Count of read pairs supporting breakpoint">\n')
		output_file.write('##INFO=<ID=gridssRPQ,Number=.,Type=String,Description="RPQ values of the gridss-called BND calls used to call this SV. Quality score of read pairs supporting breakpoint">\n')
		output_file.write('##INFO=<ID=gridssRSI,Number=.,Type=String,Description="RSI values of the gridss-called BND calls used to call this SV. Support interval offsets of partner breakend">\n')
		output_file.write('##INFO=<ID=gridssSI,Number=.,Type=String,Description="SI values of the gridss-called BND calls used to call this SV. Support interval offsets from breakend position in which at least one supporting read/read pair/assembly is mapped">\n')
		output_file.write('##INFO=<ID=gridssSR,Number=.,Type=String,Description="SR values of the gridss-called BND calls used to call this SV. Count of split reads supporting breakpoint">\n')
		output_file.write('##INFO=<ID=gridssSRQ,Number=.,Type=String,Description="SRQ values of the gridss-called BND calls used to call this SV. Quality score of split reads supporting breakpoint">\n')

	return

######################################################
def output_format_headers_for_gridss_derived_fields( output_file, output_gridss_format_fields ):

	if (output_gridss_format_fields):
		output_file.write('##FORMAT=<ID=gridssVAF,Number=.,Type=String,Description="VAF values of the gridss-called BND calls used to call this SV">\n')
		output_file.write('##FORMAT=<ID=gridssAS,Number=.,Type=String,Description="AS values of the gridss-called BND calls used to call this SV. Count of assemblies supporting breakpoint">\n')
		output_file.write('##FORMAT=<ID=gridssASQ,Number=.,Type=String,Description="ASQ values of the gridss-called BND calls used to call this SV. Quality score of assemblies supporting breakpoint">\n')
		output_file.write('##FORMAT=<ID=gridssASRP,Number=.,Type=String,Description="ASRP values of the gridss-called BND calls used to call this SV. Count of read pairs incorporated into any breakpoint assembly">\n')
		output_file.write('##FORMAT=<ID=gridssASRR,Number=.,Type=String,Description="ASRR values of the gridss-called BND calls used to call this SV. Count of split, soft clipped or indel-containing reads incorporated into any breakpoint assemblies">\n')
		output_file.write('##FORMAT=<ID=gridssBA,Number=.,Type=String,Description="BA values of the gridss-called BND calls used to call this SV. Count of assemblies supporting just local breakend">\n')
		output_file.write('##FORMAT=<ID=gridssBAQ,Number=.,Type=String,Description="BAQ values of the gridss-called BND calls used to call this SV. Quality score of assemblies supporting just local breakend">\n')
		output_file.write('##FORMAT=<ID=gridssBEID,Number=.,Type=String,Description="BEID values of the gridss-called BND calls used to call this SV. Breakend assemblies contributing support to the breakpoint">\n')
		output_file.write('##FORMAT=<ID=gridssBQ,Number=.,Type=String,Description="BQ values of the gridss-called BND calls used to call this SV. Quality score of breakend evidence">\n')
		output_file.write('##FORMAT=<ID=gridssBSC,Number=.,Type=String,Description="BSC values of the gridss-called BND calls used to call this SV. Count of soft clips supporting just local breakend">\n')
		output_file.write('##FORMAT=<ID=gridssBSCQ,Number=.,Type=String,Description="BSCQ values of the gridss-called BND calls used to call this SV. Quality score of soft clips supporting just local breakend">\n')
		output_file.write('##FORMAT=<ID=gridssBUM,Number=.,Type=String,Description="BUM values of the gridss-called BND calls used to call this SV. Count of read pairs (with one read unmapped) supporting just local breakend">\n')
		output_file.write('##FORMAT=<ID=gridssBUMQ,Number=.,Type=String,Description="BUMQ values of the gridss-called BND calls used to call this SV. Quality score of read pairs (with one read unmapped) supporting just local breakend">\n')
		output_file.write('##FORMAT=<ID=gridssCAS,Number=.,Type=String,Description="CAS values of the gridss-called BND calls used to call this SV. Count of complex compound breakpoint assemblies supporting breakpoint from elsewhere">\n')
		output_file.write('##FORMAT=<ID=gridssCASQ,Number=.,Type=String,Description="CASQ values of the gridss-called BND calls used to call this SV. Quality score of complex compound breakpoint assemblies supporting breakpoint from elsewhere">\n')
		output_file.write('##FORMAT=<ID=gridssCIEND,Number=.,Type=String,Description="CIEND values of the gridss-called BND calls used to call this SV. Confidence interval around POS for imprecise variants">\n')
		output_file.write('##FORMAT=<ID=gridssCIPOS,Number=.,Type=String,Description="CIPOS values of the gridss-called BND calls used to call this SV. Confidence interval around POS for imprecise variants">\n')
		output_file.write('##FORMAT=<ID=gridssCIRPOS,Number=.,Type=String,Description="CIRPOS values of the gridss-called BND calls used to call this SV. Confidence interval around remote breakend POS for imprecise variants">\n')
		output_file.write('##FORMAT=<ID=gridssCQ,Number=.,Type=String,Description="CQ values of the gridss-called BND calls used to call this SV. Breakpoint quality score before evidence reallocation">\n')
		output_file.write('##FORMAT=<ID=gridssFILTER,Number=.,Type=String,Description="FILTER values of the gridss-called BND calls used to call this SV">\n')
		output_file.write('##FORMAT=<ID=gridssHOMLEN,Number=.,Type=String,Description="HOMLEN values of the gridss-called BND calls used to call this SV. Length of base pair identical micro-homology at event breakpoints">\n')
		output_file.write('##FORMAT=<ID=gridssHOMSEQ,Number=.,Type=String,Description="HOMSEQ values of the gridss-called BND calls used to call this SV. Sequence of base pair identical micro-homology at event breakpoints">\n')
		output_file.write('##FORMAT=<ID=gridssIC,Number=.,Type=String,Description="IC values of the gridss-called BND calls used to call this SV. Count of read indels supporting breakpoint">\n')
		output_file.write('##FORMAT=<ID=gridssIHOMPOS,Number=.,Type=String,Description="IHOMPOS values of the gridss-called BND calls used to call this SV. Position of inexact homology">\n')
		output_file.write('##FORMAT=<ID=gridssIQ,Number=.,Type=String,Description="IQ values of the gridss-called BND calls used to call this SV. Quality score of read indels supporting breakpoint">\n')
		output_file.write('##FORMAT=<ID=gridssQUAL,Number=.,Type=String,Description="QUAL values of the gridss-called BND calls used to call this SV.">\n')
		output_file.write('##FORMAT=<ID=gridssRAS,Number=.,Type=String,Description="RAS values of the gridss-called BND calls used to call this SV. Count of assemblies supporting breakpoint from remote breakend">\n')
		output_file.write('##FORMAT=<ID=gridssRASQ,Number=.,Type=String,Description="RASQ values of the gridss-called BND calls used to call this SV. Quality score of assemblies supporting breakpoint from remote breakend">\n')
		output_file.write('##FORMAT=<ID=gridssREF,Number=.,Type=String,Description="REF values of the gridss-called BND calls used to call this SV. Count of reads mapping across this breakend">\n')
		output_file.write('##FORMAT=<ID=gridssREFPAIR,Number=.,Type=String,Description="REFPAIR values of the gridss-called BND calls used to call this SV. Count of reference read pairs spanning this breakpoint supporting the reference allele">\n')
		output_file.write('##FORMAT=<ID=gridssRP,Number=.,Type=String,Description="RP values of the gridss-called BND calls used to call this SV. Count of read pairs supporting breakpoint">\n')
		output_file.write('##FORMAT=<ID=gridssRPQ,Number=.,Type=String,Description="RPQ values of the gridss-called BND calls used to call this SV. Quality score of read pairs supporting breakpoint">\n')
		output_file.write('##FORMAT=<ID=gridssRSI,Number=.,Type=String,Description="RSI values of the gridss-called BND calls used to call this SV. Support interval offsets of partner breakend">\n')
		output_file.write('##FORMAT=<ID=gridssSI,Number=.,Type=String,Description="SI values of the gridss-called BND calls used to call this SV. Support interval offsets from breakend position in which at least one supporting read/read pair/assembly is mapped">\n')
		output_file.write('##FORMAT=<ID=gridssSR,Number=.,Type=String,Description="SR values of the gridss-called BND calls used to call this SV. Count of split reads supporting breakpoint">\n')
		output_file.write('##FORMAT=<ID=gridssSRQ,Number=.,Type=String,Description="SRQ values of the gridss-called BND calls used to call this SV. Quality score of split reads supporting breakpoint">\n')

	return

######################################################
def output_this_inline_record( output_file_ptr, inline_record ):

	# Print out the VCF record. First strip off the leading sort fields.

	infields = inline_record.split("\t")
	outline = infields[4]
	for i in range( 5, len(infields) ):
		outline = outline + "\t" + infields[i]
	outline = outline + "\n"
	output_file_ptr.write( outline )

	return

######################################################
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Identify inversions from sorted GRIDSS VCF output.')
	parser.add_argument('-i', action="store", dest="input_VCF", required=True, help='Input VCF file of breakends')
	parser.add_argument('-o', action="store", dest="output_file", required=True, help='Output VCF file of inversions')
	parser.add_argument('-u', action="store", dest="unused_breakend_records", required=False, help='Output VCF file of breakends that were not used to call inversions')
	parser.add_argument('-gap', action="store", dest="gap_or_overlap", required=False, help='Gap or overlap permitted at either breakpoint, in base pairs, for the event to still be called as an INVERSION')
	parser.add_argument('-circular', action="store", dest="look_for_circular_inversions_too", required=False, help='If present, circular inversions will be called in addition to inversions. These may be actual inversions with incorrect breakend positions that make them look circular.')
	args = parser.parse_args()

	output_gridss_info_fields = True
	output_gridss_format_fields = False

	output_file = open(args.output_file, 'w')
	if (args.unused_breakend_records is not None):
		unused_breakend_records = open(args.unused_breakend_records, 'w')
	gap_or_overlap = 200
	if (args.gap_or_overlap is not None):
		gap_or_overlap = int(args.gap_or_overlap)
	look_for_circular_inversions_too = False
	if (args.look_for_circular_inversions_too is not None):
		look_for_circular_inversions_too = True

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

	# Process each VCF record. Really they need to be processed as a batch of 4 records that identify 1 inversion

	pair1_rec1 = init_rec()
	pair1_rec2 = init_rec()
	pair2_rec1 = init_rec()
	pair2_rec2 = init_rec()

	if (continue_processing):
		while (continue_processing):
			if (in_header == True):
				if (len(inline) >= 1):
					first_char = inline[0:1]
					if (first_char != '#'):
						in_header = False

			if (in_header):

				# Before printing the final header record, print the new header records for the new INFO fields
				infields = inline.split("\t")
				infield5 = infields[4]
				if (len(infield5) >= 6):
					if (infield5[0:6] == '#CHROM'):
						output_file.write( '##ALT=<ID=INV,Description="Inversion">\n' )
						output_file.write( '##ALT=<ID=INV-DEL,Description="Inversion with deletion on one or both sides">\n' )
						output_file.write( '##ALT=<ID=INV-INS,Description="Inversion with insertion on one or both sides">\n' )
						output_file.write( '##ALT=<ID=INV-OVERLAP,Description="Inversion with overlap of uninverted sequence on one or both sides">\n' )
						output_file.write( '##ALT=<ID=INV-DEL-INS-OVERLAP,Description="Inversion with deletion, insertion, and overlap on one or both sides">\n' )
						output_file.write( '##ALT=<ID=INV-DEL-INS,Description="Inversion with deletion and insertion on one or both sides">\n' )
						output_file.write( '##ALT=<ID=INV-DEL-OVERLAP,Description="Inversion with deletion and overlap on one or both sides">\n' )
						output_file.write( '##ALT=<ID=INV-INS-OVERLAP,Description="Inversion with insertion and overlap on one or both sides">\n' )
						if (args.look_for_circular_inversions_too is not None):
							output_file.write( '##ALT=<ID=INV-CIRCULAR,Description="Inversion">\n' )
							output_file.write( '##ALT=<ID=INV-DEL-CIRCULAR,Description="Circular inversion with deletion on one or both sides">\n' )
							output_file.write( '##ALT=<ID=INV-INS-CIRCULAR,Description="Circular inversion with insertion on one or both sides">\n' )
							output_file.write( '##ALT=<ID=INV-OVERLAP-CIRCULAR,Description="Circular inversion with overlap of uninverted sequence on one or both sides">\n' )
							output_file.write( '##ALT=<ID=INV-DEL-INS-OVERLAP-CIRCULAR,Description="Circular inversion with deletion, insertion, and overlap on one or both sides">\n' )
							output_file.write( '##ALT=<ID=INV-DEL-INS-CIRCULAR,Description="Circular inversion with deletion and insertion on one or both sides">\n' )
							output_file.write( '##ALT=<ID=INV-DEL-OVERLAP-CIRCULAR,Description="Circular inversion with deletion and overlap on one or both sides">\n' )
							output_file.write( '##ALT=<ID=INV-INS-OVERLAP-CIRCULAR,Description="Circular inversion with insertion and overlap on one or both sides">\n' )

						# output_file.write('##INFO=<ID=TRANCHE,Number=1,Type=String,Description="Quality category of GRIDSS structural variant calls determined using QUAL,AS,RAS. Values are LOW INTERMEDIATE HIGH">\n') # Assume this header is already present from the input vcf
						# output_file.write('##INFO=<ID=TRANCHE2,Number=1,Type=String,Description="Quality category of GRIDSS structural variant calls determined using FILTER,SRQ,AS,RAS. Values are LOW INTERMEDIATE HIGH">\n') # Assume this header is already present from the input vcf
						# output_file.write( '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">\n' ) # Assume this header is already present from the input vcf
						# output_file.write( '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">\n' ) # Assume this header is already present from the input vcf
						# output_file.write('##INFO=<ID=VAF,Number=1,Type=String,Description="VAF of this SV call, derived from BNDVAF values of BND calls used to call this SV">\n') # Assume this header is already present from the input vcf

						output_file.write( '##INFO=<ID=INVBND1BDR5INS,Number=.,Type=String,Description="Insertion sequence following 5 prime end of breakend at inversion breakend 1 (leftmost breakend of inversion)">\n' )
						output_file.write( '##INFO=<ID=INVBND2BDR5INS,Number=.,Type=String,Description="Insertion sequence following 5 prime end of breakend at inversion breakend 2 (rightmost breakend of inversion)">\n' )
						output_file.write( '##INFO=<ID=INVBND1BDR3INS,Number=.,Type=String,Description="Insertion sequence following 3 prime end of breakend at inversion breakend 1 (leftmost breakend of inversion)">\n' )
						output_file.write( '##INFO=<ID=INVBND2BDR3INS,Number=.,Type=String,Description="Insertion sequence following 3 prime end of breakend at inversion breakend 2 (rightmost breakend of inversion)">\n' )
						output_file.write( '##INFO=<ID=INVBND1DEL,Number=.,Type=String,Description="Deletion at inversion breakend 1 (leftmost breakend of inversion)">\n' )
						output_file.write( '##INFO=<ID=INVBND2DEL,Number=.,Type=String,Description="Deletion at inversion breakend 2 (rightmost breakend of inversion)">\n' )
						output_file.write( '##INFO=<ID=INVBND1OVERLAP,Number=.,Type=String,Description="Overlap at inversion breakend 1 (leftmost breakend of inversion)">\n' )
						output_file.write( '##INFO=<ID=INVBND2OVERLAP,Number=.,Type=String,Description="Overlap at inversion breakend 2 (rightmost breakend of inversion)">\n' )
						output_file.write( '##FORMAT=<ID=POS,Number=1,Type=Integer,Description="Start position of the variant described in this record">\n' )
						output_file.write( '##FORMAT=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n' )

						# output_info_headers_for_gridss_derived_fields( output_file, output_gridss_info_fields ) # not needed, already provided by convert_VCF_BND_records_to_simple_structural_variant_VCF_records.py

						# output_file.write('##FORMAT=<ID=TRANCHE,Number=1,Type=String,Description="Quality category of GRIDSS structural variant calls determined using QUAL,AS,RAS. Values are LOW INTERMEDIATE HIGH">\n') # Assume this header is already present from the input vcf
						# output_file.write('##FORMAT=<ID=TRANCHE2,Number=1,Type=String,Description="Quality category of GRIDSS structural variant calls determined using FILTER,SRQ,AS,RAS. Values are LOW INTERMEDIATE HIGH">\n') # Assume this header is already present from the input vcf
						# output_file.write( '##FORMAT=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">\n' ) # Assume this header is already present from the input vcf
						# output_file.write( '##FORMAT=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">\n' ) # Assume this header is already present from the input vcf
						# output_file.write('##FORMAT=<ID=VAF,Number=1,Type=String,Description="VAF of this SV call, derived from BNDVAF values of BND calls used to call this SV">\n') # Assume this header is already present from the input vcf

						output_file.write( '##FORMAT=<ID=INVBND1BDR5INS,Number=.,Type=String,Description="Insertion sequence following 5 prime end of breakend at inversion breakend 1 (leftmost breakend of inversion)">\n' )
						output_file.write( '##FORMAT=<ID=INVBND2BDR5INS,Number=.,Type=String,Description="Insertion sequence following 5 prime end of breakend at inversion breakend 2 (rightmost breakend of inversion)">\n' )
						output_file.write( '##FORMAT=<ID=INVBND1BDR3INS,Number=.,Type=String,Description="Insertion sequence following 3 prime end of breakend at inversion breakend 1 (leftmost breakend of inversion)">\n' )
						output_file.write( '##FORMAT=<ID=INVBND2BDR3INS,Number=.,Type=String,Description="Insertion sequence following 3 prime end of breakend at inversion breakend 2 (rightmost breakend of inversion)">\n' )
						output_file.write( '##FORMAT=<ID=INVBND1DEL,Number=.,Type=String,Description="Deletion at inversion breakend 1 (leftmost breakend of inversion)">\n' )
						output_file.write( '##FORMAT=<ID=INVBND2DEL,Number=.,Type=String,Description="Deletion at inversion breakend 2 (rightmost breakend of inversion)">\n' )
						output_file.write( '##FORMAT=<ID=INVBND1OVERLAP,Number=.,Type=String,Description="Overlap at inversion breakend 1 (leftmost breakend of inversion)">\n' )
						output_file.write( '##FORMAT=<ID=INVBND2OVERLAP,Number=.,Type=String,Description="Overlap at inversion breakend 2 (rightmost breakend of inversion)">\n' )

						# output_format_headers_for_gridss_derived_fields( output_file, output_gridss_format_fields ) # not needed, already provided by convert_VCF_BND_records_to_simple_structural_variant_VCF_records.py

				# Print out the VCF header records. First strip off the leading sort fields.

				output_this_inline_record( output_file, inline )
				if (args.unused_breakend_records is not None):
					output_this_inline_record( unused_breakend_records, inline )

			else: # (in_header == False)

				# Process this VCF data record.
				# Read in data records until we have 4 of them. 
				# Then see if the 4 records represent an inversion.

				if (pair1_rec1['chrom'] == ''):
					pair1_rec1 = fill_record_from_inline(inline)
				elif (pair1_rec2['chrom'] == ''):
					pair1_rec2 = fill_record_from_inline(inline)
				elif (pair2_rec1['chrom'] == ''):
					pair2_rec1 = fill_record_from_inline(inline)
				elif (pair2_rec2['chrom'] == ''):
					pair2_rec2 = fill_record_from_inline(inline)

				if ( (pair1_rec1['chrom'] != '') and (pair1_rec2['chrom'] != '') and (pair2_rec1['chrom'] != '') and (pair2_rec2['chrom'] != '') ):

					if ( looks_like_INVERSION_batch(gap_or_overlap, look_for_circular_inversions_too, pair1_rec1, pair1_rec2, pair2_rec1, pair2_rec2) ):

						# The current 4 VCF data records represent an INVERSION.
						# Write out the INVERSION record.
						# Then remove the 4 VCF data records so that new records can be read in.

						inversion_record = build_INVERSION_record( pair1_rec1, pair1_rec2, pair2_rec1, pair2_rec2, output_gridss_info_fields, output_gridss_format_fields )
						outline = inversion_record + "\n"
						output_file.write( outline )

						pair1_rec1 = init_rec()
						pair1_rec2 = init_rec()
						pair2_rec1 = init_rec()
						pair2_rec2 = init_rec()

					else:

						# The current 4 VCF data records do not represent an inversion.
						# So remove the one(s) in the first position(s), shuffle the second one(s) to the first position.
						# If we are writing out unused records to a file of unused records, then write them out.
						# New data records will be read into the second position so as to consider them with the first position ones as a possible INVERSION.

						if (same_breakend(pair1_rec1, pair1_rec2) == False):

							if (args.unused_breakend_records is not None):
								output_this_inline_record( unused_breakend_records, pair1_rec1['inline'] )

							pair1_rec1 = pair1_rec2
							pair1_rec2 = pair2_rec1
							pair2_rec1 = pair2_rec2
							pair2_rec2 = init_rec()

						elif ( (same_breakend(pair1_rec1, pair1_rec2)) and (same_breakend(pair1_rec2, pair2_rec1) == False) ):

							if (args.unused_breakend_records is not None):
								output_this_inline_record( unused_breakend_records, pair1_rec1['inline'] )
								output_this_inline_record( unused_breakend_records, pair1_rec2['inline'] )

							pair1_rec1 = pair2_rec1
							pair1_rec2 = pair2_rec2
							pair2_rec1 = init_rec()
							pair2_rec2 = init_rec()

			# Now read the next VCF record

			inline = input_VCF.readline()
			inline = inline.strip()
			if (inline):
				continue_processing = True
			else:
				VCF_EOF = True
				continue_processing = False

	output_file.close()
	if (args.unused_breakend_records is not None):
		unused_breakend_records.close()

if __name__=='__main__':
    main()


