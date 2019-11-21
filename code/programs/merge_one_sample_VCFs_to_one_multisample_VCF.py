#!/usr/bin/python
# python merge_one_sample_VCFs_to_one_multisample_VCF.py sample_list
# cat AllSamples_merged.txt | python merge_one_sample_VCFs_to_one_multisample_VCF.py sample_list.txt > AllSamples_merged.vcf

# This program merges VCF files that have only 1 sample (and no more than one sample)
# to a merge VCF file having all the samples, in order of the sample_list.

# The input to this program is a file that is multiple VCF files, with key preceding each record which is the sample-id, and no VCF headers, sorted by CHROM+POS.
# This program will produce a vcf file, with the #CHROM (and without the other VCF headers)
# that has one line per CHROM+POS+ALT with the multiple samples on the one line.
# Samples will appear in the sample_list order.
# If a sample does not have a variant, then the sample field will be ./.
# For the INFO fields, this program will take the largest or highest or widest or most confident value seen in a given field, for the following fields:
# SVTYPE, END, SVLEN, OOE, LOG2COPYRATIO, TRANCHE, IMPRECISE, MEINFO, TARGETSITEDUPL
# ALL OTHER INFO FIELDS ARE LOST!!!
# This program assumes that the field appearing in the INFO field are all the same for multiple VCF records.
# The highest QUAL value will be output for a given merged VCF record.

# This program assumes that the format fields for samples are all the same.
# Multiple sample records will be merged by this program without changing the format field or sample fields.
# If merged sample records do not have the same format fields, then the sample fields will not all match the one format field of the merged record.

# Example sample list input:
# AAAPA
# AAAPB
# AAAPC

# Example VCF input with sample key:
# AAAPA   1       10027   .       A       <CNV>   50      .       SVTYPE=CNV;END=10447;SVLEN=420;OOE=13;LOG2COPYRATIO=3.66662973467533     GT:END:SVL:CN:L2CR:QUAL 1/.:1
# AAAPB   1       10027   .       A       <CNV>   50      .       SVTYPE=CNV;END=10447;SVLEN=420;OOE=15;LOG2COPYRATIO=3.85886039149222     GT:END:SVL:CN:L2CR:QUAL 1/.:1
# AAAPC   1       10027   .       A       <CNV>   50      .       SVTYPE=CNV;END=10447;SVLEN=420;OOE=13;LOG2COPYRATIO=3.72570121120816     GT:END:SVL:CN:L2CR:QUAL 1/.:1
# AAAPA   1       10448   .       C       <CNV>   0       .       SVTYPE=CNV;END=16884326;SVLEN=16873878;OOE=1.0;LOG2COPYRATIO=-0.00412586019128133        GT:END:SVL:CN
# AAAPB   1       10448   .       C       <CNV>   0       .       SVTYPE=CNV;END=2585035;SVLEN=2574587;OOE=1.0;LOG2COPYRATIO=0.00340339629079111   GT:END:SVL:CN:L2CR:QU
# AAAPC   1       10448   .       C       <CNV>   3       .       SVTYPE=CNV;END=16884326;SVLEN=16873878;OOE=1.0;LOG2COPYRATIO=-0.00205906279490569        GT:END:SVL:CN
# AAAPB   1       2585036 .       C       <CNV>   50      .       SVTYPE=CNV;END=2629608;SVLEN=44572;OOE=3.4;LOG2COPYRATIO=1.74987438660141        GT:END:SVL:CN:L2CR:QU
# AAAPB   1       2629609 .       A       <CNV>   1       .       SVTYPE=CNV;END=16884326;SVLEN=14254717;OOE=0.99;LOG2COPYRATIO=-0.0155558046716885        GT:END:SVL:CN

# Example VCF input with sample key:
# AAAAA   1       829169  gridss0_4477    A       <DEL>   1556.54 .       SVTYPE=DEL;END=829205;SVLEN=-36;TRANCHE=HIGH    GT      1/.
# AAAAB   1       829169  gridss0_4883    A       <DEL>   1462.6  .       SVTYPE=DEL;END=829205;SVLEN=-36;TRANCHE=HIGH    GT      1/.
# AAAAD   1       829171  gridss0_3246    A       <DEL>   863.14  .       SVTYPE=DEL;END=829207;SVLEN=-36;TRANCHE=INTERMEDIATE    GT      1/.
# AAAAB   1       1176029 gridss0_15836   C       <DUP:TANDEM>    541.8   .       SVTYPE=DUP;END=1176071;SVLEN=42;TRANCHE=INTERMEDIATE    GT      1/.
# AAAAA   1       1613427 gridss0_5403    C       <DEL>   795.52  .       SVTYPE=DEL;END=1613462;SVLEN=-35;TRANCHE=INTERMEDIATE   GT      1/.
# AAAAA   1       2324467 gridss0_15184   T       <DUP:INS>       2073.33 .       SVTYPE=DUP;END=2324498;SVLEN=41;TRANCHE=HIGH;INSSEQ=CAGGCTTCAG  GT      1/.
# AAAAB   1       2324467 gridss0_16551   T       <DUP:INS>       1688.01 .       SVTYPE=DUP;END=2324498;SVLEN=41;TRANCHE=HIGH;INSSEQ=CAGGCTTCAG  GT      1/.
# AAAAD   1       2324467 gridss0_13376   T       <DUP:INS>       1367.75 .       SVTYPE=DUP;END=2324498;SVLEN=41;TRANCHE=HIGH;INSSEQ=CAGGCTTCAG  GT      1/.
# AAAAB   1       2533555 gridss0_16659   T       <DUP:INS>       762.25  .       SVTYPE=DUP;END=2533628;SVLEN=79;TRANCHE=INTERMEDIATE;INSSEQ=TTCGTC      GT      1/.
# AAAAB   1       3097279 gridss0_7489    A       <INDEL> 851.88  .       SVTYPE=INDEL;END=3097285;SVLEN=36;TRANCHE=INTERMEDIATE;INSSEQ=GGGGTCTTCCAGACAGAAGAGAACAGCGTTCTCCTGAGGGGT        
# AAAAB   1       3147940 gridss0_16953o  C       ]1:3147978]C    514.96  .       AS=1;ASQ=182.71;ASRP=0;ASSR=20;BA=0;BAQ=0.00;BEID=asm0-173245,asm0-40964;BQ=55.67;BSC=1;BSCQ=15.91;BUM=2
# AAAAB   1       3418582 gridss0_17113   C       <DUP:INS>       2182.32 .       SVTYPE=DUP;END=3418666;SVLEN=114;TRANCHE=HIGH;INSSEQ=GTGGTCTCGTGGCGAGGGCGCCCCCCGCCT     GT      1/.
# AAAAB   2       131886724       gridss38_477    C       <INV-DEL-OVERLAP-CIRCULAR>      1993.44 .       SVTYPE=INV;END
# AAAAB   2       184569193       gridss43_19225  C       <INV-DEL-OVERLAP>       1080.33 .       SVTYPE=INV;END=1845720
# AAAAB   20      18591561        gridss274_6843  A       <INV-DEL-OVERLAP>       842.12  .       SVTYPE=INV;END=1859588
# AAAAA   21      27374159        gridss281_7374  T       <INV-DEL-OVERLAP>       2447.7  .       SVTYPE=INV;END=2737470
# AAAAB   21      27374159        gridss281_7636  T       <INV-DEL-OVERLAP>       3109.95 .       SVTYPE=INV;END=2737470
# AAAAD   21      27374159        gridss281_7586  T       <INV-DEL-OVERLAP>       2803.75 .       SVTYPE=INV;END=2737470
# AAAAD   3       44740979        gridss53_14613  C       <INV>   1257.37 .       SVTYPE=INV;END=44742297;TRANCHE=HIGH
# AAAAB   5       147554239       gridss103_15444 A       <INV-DEL>       3319.56 .       SVTYPE=INV;END=147554615;TRANC
# AAAAA   5       147554240       gridss103_13797 T       <INV-DEL>       2849.55 .       SVTYPE=INV;END=147554615;TRANC
# AAAAD   5       147554240       gridss103_11491 T       <INV-DEL>       1479.56 .       SVTYPE=INV;END=147554615;TRANC

# Example VCF input with sample key:
# D02527433       1       44736470        .       C       <INS:ME:L1>     24      .       SVTYPE=INS;MEINFO=L1,0,3,.;TARGETSITEDUPL=duplication   GT:MEIQUAL:MEI5MB:MEI3MB:MEIDU
# D02527431       1       45963313        .       G       <INS:ME:SVA>    15      .       IMPRECISE;SVTYPE=INS;MEINFO=SVA,-94,573,.;TARGETSITEDUPL=unknown        GT:MEIQUAL:MEI
# D02527432       1       45963394        .       C       <INS:ME:SVA>    18      .       IMPRECISE;SVTYPE=INS;MEINFO=SVA,-20,565,.;TARGETSITEDUPL=unknown        GT:MEIQUAL:MEI
# D02527433       1       45963394        .       C       <INS:ME:SVA>    20      .       IMPRECISE;SVTYPE=INS;MEINFO=SVA,-20,401,.;TARGETSITEDUPL=unknown        GT:MEIQUAL:MEI
# D02527431       1       46169013        .       G       <INS:ME:ALU>    21      .       SVTYPE=INS;MEINFO=ALU,-20,20,.;TARGETSITEDUPL=unknown   GT:MEIQUAL:MEI5MB:MEI3MB:MEIDU
# D02527432       1       46169013        .       G       <INS:ME:ALU>    20      .       SVTYPE=INS;MEINFO=ALU,-20,20,.;TARGETSITEDUPL=unknown   GT:MEIQUAL:MEI5MB:MEI3MB:MEIDU
# D02527433       1       46169013        .       G       <INS:ME:ALU>    21      .       SVTYPE=INS;MEINFO=ALU,-20,20,.;TARGETSITEDUPL=unknown   GT:MEIQUAL:MEI5MB:MEI3MB:MEIDU
# D02527433       1       142915489       .       G       <INS:ME:HERV>   17      .       IMPRECISE;SVTYPE=INS;MEINFO=HERV,-103,633,.;TARGETSITEDUPL=unknown      GT:MEIQUAL:MEI
# D02527431       1       142915491       .       T       <INS:ME:HERV>   14      .       IMPRECISE;SVTYPE=INS;MEINFO=HERV,-105,583,.;TARGETSITEDUPL=unknown      GT:MEIQUAL:MEI
# D02527431       1       149308416       .       G       <INS:ME:L1>     20      .       IMPRECISE;SVTYPE=INS;MEINFO=L1,-2,4,.;TARGETSITEDUPL=unknown    GT:MEIQUAL:MEI5MB:MEI3
# D02527433       1       149308416       .       G       <INS:ME:L1>     20      .       SVTYPE=INS;MEINFO=L1,-6,4,.;TARGETSITEDUPL=unknown      GT:MEIQUAL:MEI5MB:MEI3MB:MEIDU
# D02527432       1       149308432       .       T       <INS:ME:L1>     17      .       SVTYPE=INS;MEINFO=L1,-12,4,.;TARGETSITEDUPL=duplication GT:MEIQUAL:MEI5MB:MEI3MB:MEIDU

# An example of creating the input files to this program:
# cat AAAPA_bicseq_CNV_with_sample_key.txt AAAPB_bicseq_CNV_with_sample_key.txt AAAPC_bicseq_CNV_with_sample_key.txt | sort -k2,2 -k3,3n -k6,6 > AllSamples_merged.txt
# echo AAAPA > sample_list.txt
# echo AAAPB >> sample_list.txt
# echo AAAPC >> sample_list.txt

import sys
import os
import commands
import re

DELIMITER_FOR_INFO_MEINFO_FIELD = ','

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
def extract_info_fields( info_string ):

	global DELIMITER_FOR_INFO_MEINFO_FIELD

	info_fields = {}
	info_string_split = info_string.split(';')
	for info_key_and_value in info_string_split:
		bits = info_key_and_value.split('=')
		info_key = bits[0]
		info_value = ''
		if (len(bits) > 1):
			info_value = bits[1]
		if (info_key == 'SVTYPE'):
			info_fields['SVTYPE'] = info_value
		if (info_key == 'END'):
			info_fields['END'] = info_value
		if (info_key == 'SVLEN'):
			info_fields['SVLEN'] = info_value
		if (info_key == 'OOE'):
			info_fields['OOE'] = info_value
		if (info_key == 'LOG2COPYRATIO'):
			info_fields['LOG2COPYRATIO'] = info_value
		if (info_key == 'MEINFO'):
			info_fields['MEINFO'] = info_value
		if (info_key == 'TARGETSITEDUPL'):
			info_fields['TARGETSITEDUPL'] = info_value
		if (info_key == 'IMPRECISE'):
			info_fields['IMPRECISE'] = info_value

	return info_fields

######################################################
def convert_alt_to_its_simple_type( this_alt ):

	return_alt = str(this_alt)
	if (this_alt.find('[') >= 0):
		return_alt = '<BND>'
	elif (this_alt.find(']') >= 0):
		return_alt = '<BND>'
	elif (this_alt.find('<INS:ME:') >= 0):
		return_alt = '<MEI>'
	
	return return_alt

######################################################
def are_these_two_position_variants_the_same( chrom1, pos1, alt1, chrom2, pos2, alt2 ):

	they_are_the_same = False
	chrom1 = str(chrom1)
	chrom2 = str(chrom2)
	if ((chrom1 == chrom2) and (pos1 == pos2)):
		if (alt1 == alt2):
			they_are_the_same = True
		else:
			alt1_list = alt1.split(',')
			alt2_list = alt2.split(',')
			for i in range( 0, len(alt1_list) ):
				this_alt = alt1_list[i]
				alt1_list[i] = convert_alt_to_its_simple_type(alt1_list[i])
			for i in range( 0, len(alt2_list) ):
				this_alt = alt2_list[int(i)]
				alt2_list[i] = convert_alt_to_its_simple_type(alt2_list[i])
			found_non_match = False
			for this_alt in alt1_list:
				if (this_alt != alt2_list[0]):
					found_non_match = True
			for this_alt in alt2_list:
				if (this_alt != alt1_list[0]):
					found_non_match = True
			they_are_the_same = True
			if (found_non_match):
				they_are_the_same == False

	return they_are_the_same

######################################################
def merge_alts( alt1, alt2 ):

	alt1_list = alt1.split(',')
	alt2_list = alt2.split(',')
	for i in range( 0, len(alt1_list) ):
		this_alt1 = alt1_list[i]
		alt1_list[i] = convert_alt_to_its_simple_type(alt1_list[i])
	for i in range( 0, len(alt2_list) ):
		this_alt = alt2_list[i]
		alt2_list[i] = convert_alt_to_its_simple_type(alt2_list[i])

	return_alt_list = alt1_list
	for i in range( 0, len(alt2_list) ):
		this_alt2 = alt2_list[i]
		found_this_alt2_in_return_list = False
		for j in range( 0, len(return_alt_list) ):
			this_return_alt = return_alt_list[j]
			if (this_alt2 == this_return_alt):
				found_this_alt2_in_return_list = True
		if (found_this_alt2_in_return_list == False):
			return_alt_list.append( this_alt2 )

	return_alt = return_alt_list[0]
	if (len(return_alt_list) > 1):
		for i in range( 1, len(return_alt_list) ):
			return_alt = return_alt + ',' + return_alt_list[i]

	return return_alt

######################################################
def merge_quals( qual1, qual2 ):

	return_qual = qual1
	if (qual1 != '.'):
		if (qual2 != '.'):
			float_qual1 = float(qual1)
			float_qual2 = float(qual2)
			if (float_qual2 > float_qual1):
				return_qual = qual2
	elif (qual2 != '.'):
		return_qual = qual2
	return return_qual

######################################################
def merge_filters( filter1, filter2 ):

	return_filter = filter1
	filter1 = str(filter1)
	filter2 = str(filter2)
	if ((filter1 == 'PASS') or (filter2 == 'PASS')):
		return_filter = 'PASS'
	else:
		if (filter1 == '.'):
			if (filter2 != '.'):
				return_filter = filter2
	return return_filter

######################################################
def merge_info_fields( info_fields1, info_fields2 ):

	return_info_fields = info_fields1
	for info_key in info_fields2:
		info_field1_value = ''
		if (info_key in info_fields1):
			info_field1_value = info_fields1[info_key]
		info_field2_value = info_fields2[info_key]
		if (info_key == 'SVTYPE'):
			if (info_field1_value != info_field2_value):
				raise ValueError("\n\nTwo records to be merged do not have the same INFO.SVTYPE: " + str(info_field1_value) + ' and ' + str(info_field2_value) + ". Will not continue processing.\n")
		elif (info_key == 'END'):
			info_field1_value_numeric = int(info_field1_value)
			info_field2_value_numeric = int(info_field2_value)
			if (info_field2_value_numeric > info_field1_value_numeric):
				return_info_fields['END'] = str(info_field2_value)
		elif (info_key == 'SVLEN'):
			info_field1_value_numeric = abs(int(info_field1_value))
			info_field2_value_numeric = abs(int(info_field2_value))
			if (info_field2_value_numeric > info_field1_value_numeric):
				return_info_fields['SVLEN'] = str(info_field2_value)
		elif (info_key == 'OOE'):
			info_field1_value_numeric = float(info_field1_value)
			info_field2_value_numeric = float(info_field2_value)
			if (info_field2_value_numeric > info_field1_value_numeric):
				return_info_fields['OOE'] = str(info_field2_value)
		elif (info_key == 'LOG2COPYRATIO'):
			info_field1_value_numeric = float(info_field1_value)
			info_field2_value_numeric = float(info_field2_value)
			if (info_field2_value_numeric > info_field1_value_numeric):
				return_info_fields['LOG2COPYRATIO'] = str(info_field2_value)
		elif (info_key == 'MEINFO'):
			return_info_fields['MEINFO'] = choose_best_meinfo( info_field1_value, info_field2_value )
		elif (info_key == 'TARGETSITEDUPL'):
			return_info_fields['TARGETSITEDUPL'] = choose_best_targetsitedupl( info_field1_value, info_field2_value )
	# If any samples are PRECISE, then don't keep the IMPRECISE key
	if 'IMPRECISE' in return_info_fields:
		del return_info_fields['IMPRECISE']
	if (('IMPRECISE' in info_fields1) and ('IMPRECISE' in info_fields2)):
		return_info_fields['IMPRECISE'] = ''

	return return_info_fields

######################################################
def choose_best_meinfo( meinfo1, meinfo2 ): # eg. L1,-179,180,.

	global DELIMITER_FOR_INFO_MEINFO_FIELD

	new_meinfo = [ '', 0, 0, '' ] # MEI-type , start-pos , end-pos , polarity
	new_meinfo_string = ''
	if ((meinfo1 == '') and (meinfo2 == '')):
		new_meinfo_string = ''
	elif (meinfo1 == ''):
		new_meinfo_string = meinfo2
	elif (meinfo2 == ''):
		new_meinfo_string = meinfo1
	else:
		meinfo1_fields = meinfo1.split( DELIMITER_FOR_INFO_MEINFO_FIELD )
		meinfo2_fields = meinfo2.split( DELIMITER_FOR_INFO_MEINFO_FIELD )
		new_meinfo[0] = meinfo1_fields[0]
		new_meinfo[1] = int(meinfo1_fields[1])
		if (int(meinfo2_fields[1]) < int(meinfo1_fields[1])):
			new_meinfo[1] = int(meinfo2_fields[1])
		new_meinfo[2] = int(meinfo1_fields[2])
		if (int(meinfo2_fields[2]) > int(meinfo1_fields[2])):
			new_meinfo[2] = int(meinfo2_fields[2])
		new_meinfo[3] = meinfo1_fields[3]
		if ((meinfo1_fields[3] == '.') and (meinfo2_fields[3] != '.')):
			new_meinfo[3] = meinfo2_fields[3]
		new_meinfo_string = new_meinfo[0] + DELIMITER_FOR_INFO_MEINFO_FIELD + str(new_meinfo[1]) + DELIMITER_FOR_INFO_MEINFO_FIELD + str(new_meinfo[2]) + DELIMITER_FOR_INFO_MEINFO_FIELD + new_meinfo[3]
        return new_meinfo_string

######################################################
def choose_best_targetsitedupl( targetsitedupl1, targetsitedupl2 ): # values are duplication, noTSD, or unknown

	new_targetsitedupl = targetsitedupl1
	if (targetsitedupl1 == 'unknown'):
		new_targetsitedupl = targetsitedupl2
	elif ((targetsitedupl1 == 'noTSD') and (targetsitedupl2 == 'duplication')):
		new_targetsitedupl = 'unknown'
	elif ((targetsitedupl2 == 'noTSD') and (targetsitedupl1 == 'duplication')):
		new_targetsitedupl = 'unknown'
        return new_targetsitedupl

######################################################
def extract_gt_fields( gt_string ):

	gt_fields = []
	gt_separators = []
	rest_of_gt_string = gt_string
	look_for_more_gt_fields = True
	while (look_for_more_gt_fields):

		separator_position = -1
		separator = ''
		slash_position = rest_of_gt_string.find('/')
		bar_position = rest_of_gt_string.find('|')
		if (slash_position >= 0):
			separator_position = slash_position
			separator = '/'
		if (bar_position >= 0):
			if (bar_position < separator_position):
				separator_position = bar_position
				separator = '|'

		if (separator_position == -1):
			look_for_more_gt_fields = False
			gt_fields.append( rest_of_gt_string )
		else:
			this_gt = rest_of_gt_string[ 0:separator_position ]
			gt_fields.append( this_gt )
			gt_separators.append( separator )
			rest_of_gt_string = rest_of_gt_string[ (separator_position+1): ]

	return gt_fields, gt_separators

######################################################
def determine_merged_gt_of_sample( old_gt_fields, old_alt_string, new_alt_string ):

	new_gt_fields = []
	old_alt_fields = old_alt_string.split(',')
	new_alt_fields = new_alt_string.split(',')
	for old_gt_field in old_gt_fields:
		new_gt_field = '.'
		if (old_gt_field != '.'):
			old_gt_idx = int(old_gt_field) - 1
			alt_for_old_gt = old_alt_fields[old_gt_idx]
			new_gt_idx = -1
			for i in range( 0, len(new_alt_fields) ):
				if (alt_for_old_gt == (new_alt_fields[i])):
					new_gt_idx = i
			alt_for_new_gt = new_alt_fields[new_gt_idx]
			new_gt_field = (new_gt_idx + 1)
		new_gt_fields.append( new_gt_field )
	
	return new_gt_fields

######################################################
def write_out_variant_record( this_chrom, this_pos, this_id, this_ref, this_alt, this_qual, this_filter, this_info_fields, this_format, this_sample_list ):

	if (this_chrom != '.'):
		outline = str(this_chrom) + "\t" + str(this_pos) + "\t" + this_id + "\t" + this_ref + "\t" + this_alt + "\t" + this_qual + "\t" + this_filter
		this_info_string = ''
		for info_key in sorted(this_info_fields):
			this_info_pair = str(info_key)
			if (this_info_fields[info_key] != ''):
				this_info_pair = this_info_pair + '=' + str(this_info_fields[info_key])
			if (this_info_string == ''):
				this_info_string = this_info_pair
			else:
				this_info_string = this_info_string + ';' + this_info_pair
		outline = outline + "\t" + this_info_string + "\t" + this_format
		if (len(this_sample_list) > 10):
			for i in range( 10, len(this_sample_list) ):
				this_sample_field = this_sample_list[i]
				outline = outline + "\t" + this_sample_field
		outline = outline + "\n"
		sys.stdout.write( outline )
	return

######################################################
def main():

	input_sample_list = sys.argv[1]

	# Get the list of samples, and output a VCF header for them

	in_sample_list = open(input_sample_list, 'r')
	temp_list_of_samples = in_sample_list.readlines()
	in_sample_list.close()
	list_of_samples = []
	for i in range( 0, len(temp_list_of_samples) ):
		this_sample = temp_list_of_samples[i]
		this_sample = this_sample.strip()
		if (this_sample != ''):
			list_of_samples.append( this_sample )
	list_of_samples_idx = {}
	outline = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"

	for i in range( 0, len(list_of_samples) ):
		this_sample = list_of_samples[i]
		this_sample = this_sample.strip()
		list_of_samples_idx[ this_sample ] = i
		outline = outline + "\t" + this_sample
	outline = outline + "\n"
	sys.stdout.write( outline )

	# Read in the input VCF file from STDIN

	prev_chrom = '.'
	prev_pos = -1
	prev_id = '.'
	prev_ref = '.'
	prev_alt = '.'
	prev_qual = '.'
	prev_filter = '.'
	prev_info_fields = {}
	prev_info_fields['SVTYPE'] = []
	prev_info_fields['END'] = []
	prev_info_fields['SVLEN'] = []
	prev_info_fields['OOE'] = []
	prev_info_fields['LOG2COPYRATIO'] = []
	prev_info_fields['MEINFO'] = []
	prev_info_fields['TARGETSITEDUPL'] = []
	prev_info_fields['IMPRECISE'] = []
	prev_format = '.'
	prev_sample_list = ['./.'] * (len(list_of_samples) + 10) # first 10 positions of array are dummy, are ignored, are there for ease of indices

	in_header = True
	for inline in sys.stdin:

		if (in_header == True):
			if (len(inline) >= 1):
				first_char = inline[0:1]
				if (first_char != '#'):
					in_header = False

		if (in_header == False): # We are processing VCF data records. We are no longer in the header part of the file.

			inline = inline.strip()
			if (inline != ''):

				infields = inline.split("\t")
				this_added_key = str(infields[0])
				this_chrom = str(infields[1])
				this_pos = int(infields[2])
				this_id = infields[3]
				this_ref = infields[4]
				this_alt = infields[5]
				this_qual = str(infields[6])
				this_filter = infields[7]
				this_info_string = infields[8]
				this_info_fields = extract_info_fields( this_info_string )
				this_format = infields[9]
				this_sample_string = infields[10]
				this_sample_list = ['./.'] * (len(list_of_samples) + 10)
				index_of_this_sample_in_sample_list = -1
				index_of_this_sample_in_VCF_file_having_added_key = -1
				for i in range( 0, len(list_of_samples) ):
					if (this_added_key == list_of_samples[i]):
						index_of_this_sample_in_sample_list = i
						index_of_this_sample_in_VCF_file_having_added_key = i + 10
				this_sample_list[index_of_this_sample_in_VCF_file_having_added_key] = this_sample_string

				this_position_variant_same_as_prev = are_these_two_position_variants_the_same( this_chrom, this_pos, this_alt, prev_chrom, prev_pos, prev_alt )

				if (this_position_variant_same_as_prev):

					# Merge this variant's ALT field into the previous merged record's ALT

					prev_alt = merge_alts( this_alt, prev_alt )
					prev_qual = merge_quals( this_qual, prev_qual)
					prev_filter = merge_filters( this_filter, prev_filter )
					prev_info_fields = merge_info_fields( this_info_fields, prev_info_fields )

					# For each sample in this variant, change the GT fields to reflect the merged ALT values

					this_sample_fields = this_sample_string.split(':')
					gt_string = this_sample_fields[0]
					gt_fields, gt_separators = extract_gt_fields( gt_string )
					merged_gt_fields = determine_merged_gt_of_sample( gt_fields, this_alt, prev_alt )
					merged_gt_string = merged_gt_fields[0]
					if (len(merged_gt_fields) > 1):
						for j in range( 1, len(merged_gt_fields) ):
							merged_gt_string = str(merged_gt_string) + gt_separators[j-1] + str(merged_gt_fields[j])
					this_sample_string_converted_to_its_merged_values = merged_gt_string
					if (len(this_sample_fields) > 1):
						for i in range( 1, len(this_sample_fields) ):
							this_sample_string_converted_to_its_merged_values = this_sample_string_converted_to_its_merged_values + ':' + str(this_sample_fields[i])
					this_sample_list[index_of_this_sample_in_VCF_file_having_added_key] = this_sample_string_converted_to_its_merged_values

					# Once each sample has had its GT field changed for correct values of the merged variant record,
					# add this sample to the merged variant record.

					# Some pipelines, such as mobile elements, may end up calling more than one variant at the same position for the same variant.
					# In that case, we'll take only one of them.
					# if (prev_sample_list[index_of_this_sample_in_VCF_file_having_added_key] != './.'):
					# 	raise ValueError("\n\nFor a given variant and sample, there should be only 1 or no records. Sample " + list_of_samples[index_of_this_sample_in_sample_list] + " has more than one record for variant at " + this_chrom + ':' + str(this_pos) + ". Will not continue processing.\n")

					prev_sample_list[index_of_this_sample_in_VCF_file_having_added_key] = this_sample_string

				else: # this_position_variant_same_as_prev == False

					# Write out the previous variant record. 
					# This new variant will become the previous variant record that subsequent variants might merge with.

					write_out_variant_record( prev_chrom, prev_pos, prev_id, prev_ref, prev_alt, prev_qual, prev_filter, prev_info_fields, prev_format, prev_sample_list )

					prev_chrom = this_chrom
					prev_pos = this_pos
					prev_id = this_id
					prev_ref = this_ref
					prev_alt = this_alt
					prev_qual = this_qual
					prev_filter = this_filter
					prev_info_fields = this_info_fields
					prev_format = this_format
					prev_sample_list = this_sample_list

	# Write out the last variant record.

	write_out_variant_record( prev_chrom, prev_pos, prev_id, prev_ref, prev_alt, prev_qual, prev_filter, prev_info_fields, prev_format, prev_sample_list )


if __name__=='__main__':
    main()


