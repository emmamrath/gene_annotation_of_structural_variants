#!/usr/bin/python
# python convert_VCF_BND_records_to_simple_structural_variant_VCF_records.py
# cat input.vcf | python convert_VCF_BND_records_to_simple_structural_variant_VCF_records.py > output.vcf

# This program takes 1 VCF file on input from stdin, and writes out a VCF to stdout.
# It converts low-level BND records to structural variant records.
# The types of higher-level structural variant records that this program produces are:
# deletion, insertion, indel (insertion and deletion), tandem duplication, duplication having an inserted sequence too
# These are output as a new VCF record having ALT respectively of:
# <DEL>, <INS>, <INDEL>, <DUP:TANDEM>, <DUP:INS>
# Low-level BND records that are not used to create the higher-level records are discarded and do not appear in the output.

# This program assumes that the input VCF is sorted by GRIDSS event in the following way:
# grep '^#' mysample_gridss.vcf > mysample_gridss_hdr.vcf
# grep -v '^#' mysample_gridss.vcf | sort -t$'\t' -k 3,3 > mysample_gridss_nohdr_sortedByEvent.txt
# cat mysample_gridss_hdr.vcf mysample_gridss_nohdr_sortedByEvent.txt > mysample_gridss_sortedByEvent.vcf

# After running this program, the output should be sorted in CHROM+POS as per a normal VCF file:
# grep '^#' mysample_INS_DEL_INDEL_DUP_INV_BND.vcf > mysample_gridss_INS_DEL_INDEL_DUP_INV_BND_hdr.txt
# grep -v '^#'mysample_INS_DEL_INDEL_DUP_INV_BND.vcf > mysample_INS_DEL_INDEL_DUP_INV_BND_nohdr.txt
# sort -t$'\t' -k1,1 -k2,2n mysample_INS_DEL_INDEL_DUP_INV_BND_nohdr.txt > mysample_INS_DEL_INDEL_DUP_INV_BND_nohdr_sorted.txt
# cat mysample_INS_DEL_INDEL_DUP_INV_BND_hdr.txt mysample_INS_DEL_INDEL_DUP_INV_BND_nohdr_sorted.txt > mysample_INS_DEL_INDEL_DUP_INV_BND_sorted.vcf

# Although the VCF spec says to put inserted/deleted sequences in the ALT field where possible,
# this program does not do that, becuase this program is designed to identify large structural variants (SV).
# We already have GATK HaplotypeCaller that does a good job of identifying small INDELs.
# This program here is utilising GRIDDS output to identify large SVs and large INDELs.
# This program is designed for intrachromosomal structural variants, not interchromosomal variants.
# Thus, this program does not try to identify SVs resulting from breakends (BND) on different chromosomes.

# This program assumes that multiple VCF records for the same event appear subsequent to each other in the VCF file.
# If this is not the case, then the VCF record needs to be sorted by INFO.EVENT so that it will become the case.
# IF the input is from GRIDSS, it appears that sorting by the VCF.ID field will effectively sort it by VCF.INFO.EVENT
# Eg.
# grep '^#' mySample.gridss.vcf > mySample.gridss_hdr.vcf
# grep -v '^#' mySample.gridss.vcf | sort -t$'\t' -k 3,3 > mySample.gridss_sorted_nohdr.vcf
# cat mySample.gridss_hdr.vcf mySample.gridss_sorted_nohdr.vcf > mySample.gridss_sorted.vcf

# An example of how to run this program, including sorting before and after:
# 
# grep '^#' AOCS_92_32_blood_sorted_markedDupl.sv.vcf > blood_gridss_hdr.vcf
# grep -v '^#' AOCS_92_32_blood_sorted_markedDupl.sv.vcf | sort -t$'\t' -k 3,3 > blood_gridss_sortedByEvent_nohdr.txt
# cat blood_gridss_hdr.vcf blood_gridss_sortedByEvent_nohdr.txt > blood_gridss_sortedByEvent.vcf
# cat blood_gridss_sortedByEvent.vcf | python convert_VCF_BND_records_to_simple_structural_variant_VCF_records.py > blood_gridss_sortedByEvent_INS_DEL_INDEL_DUP.vcf
# 
# grep '^#' blood_gridss_sortedByEvent_INS_DEL_INDEL_DUP.vcf > blood_gridss_INDELetc_hdr.vcf
# grep -v '^#' blood_gridss_sortedByEvent_INS_DEL_INDEL_DUP.vcf > blood_gridss_INDELetc_nohdr.txt
# sort -t$'\t' -k1,1 -k2,2n blood_gridss_INDELetc_nohdr.txt > blood_gridss_INDELetc_nohdr_sorted.txt
# cat blood_gridss_INDELetc_hdr.vcf blood_gridss_INDELetc_nohdr_sorted.txt > blood_gridss_INS_DEL_INDEL_DUP.vcf

# This program assume that each VCF record has one and only one sample.

# When this program produces one output VCF record for a group of input VCF records,
# it assigns a TRANCHE "score" of LOW, INTERMEDIATE, or HIGH corresponding to the confidence in this higher-level structural variant call.
# Other scores, quality and precision fields are not retained.

# When this program produces one output VCF record for a group of input VCF records,
# the output INFO and FORMAT fields fields for counts and quality will contain the highest value from the input records.

# When this program produces one output VCF record for a group of input VCF records,
# at first genotype of the sample (GT) is set to 1/.
# even though the GRIDSS GT appears to be set to a period.
# If any of the BND records have an INFO BNDVAF field, then this is used to set the VAF.
# For INS,        if any BNDVAF > 0.85, then GT is set to 1/1, otherwise GT is set to 0/1.
# For DEL,        if any BNDVAF > 0.85, then GT is set to 1/1, otherwise GT is set to 0/1. (Would be nice to first verify 1/1 with a copy number loss.)
# For INDEL,      if any BNDVAF > 0.85, then GT is set to 1/1, otherwise GT is set to 0/1. (Would be nice to first verify 1/1 with a copy number loss.)
# For DUP:TANDEM, if any BNDVAF > 0.65, then GT is set to 1/1, otherwise GT is set to 0/1. (Would be nice to first verify 1/1 with a copy number gain.)
# For DUP:INS,    if any BNDVAF > 0.65, then GT is set to 1/1, otherwise GT is set to 0/1. (Would be nice to first verify 1/1 with a copy number gain.)
# TRANCHE2 for the INS/DEL/INDEL/DUP/INV call will be set to the highest TRANCHE2 value of the BND records contributing to the call.

# Here are the conversions that this program can do:

# CHROM	POS	ID	REF	ALT			INFO					CHROM	POS	ID	REF	ALT		INFO
# 1	124001	g3	G	G[1:124020[		EVENT=E2;PARID=g3;SVTYPE=BND	} ==>	1	124001	event	G	<DEL>		SVTYPE=DEL;END=124019;SVLEN=-18
# 1	124020	g4	C	]1:124001]C		EVENT=E2;PARID=g4;SVTYPE=BND	}

# CHROM	POS	ID	REF	ALT			INFO					CHROM	POS	ID	REF	ALT		INFO
# 1	1146	g5	G       ]1:1146]ACGGGGGTTCTG   	EVENT=E3;PARID=g5;SVTYPE=BND	} ==>	1	1146	event	G	<INS>		SVTYPE=INS;END=1146;SVLEN=11
# 1	1146	g6	G       GACGGGGGTTCT[1:1146[	EVENT=E3;PARID=g6;SVTYPE=BND	}

# CHROM	POS	ID	REF	ALT			INFO					CHROM	POS	ID	REF	ALT		INFO
# 1	534	g7	G       GGGCTGC[1:533[   	EVENT=E3;PARID=g6;SVTYPE=BND	} ==>	1	533	event	G	<INS>		SVTYPE=INS;END=533;SVLEN=6
# 1	533	g8	G       ]1:534]GGCTGCA		EVENT=E3;PARID=g6;SVTYPE=BND	}

# CHROM	POS	ID	REF	ALT			INFO					CHROM	POS	ID	REF	ALT		INFO
# 1	124001	g3	G	GAAAATTTT[1:124020[	EVENT=E2;PARID=g3;SVTYPE=BND	} ==>	1	124001	event	G	<INDEL>		SVTYPE=INDEL;END=124019;SVLEN=-10
# 1	124020	g4	C	]1:124001]AAAATTTTC	EVENT=E2;PARID=g4;SVTYPE=BND	}

# CHROM	POS	ID	REF	ALT			INFO					CHROM	POS	ID	REF	ALT		INFO
# 1	124001	g3	G	GAAAATTTT[1:124004[	EVENT=E2;PARID=g3;SVTYPE=BND	} ==>	1	124001	event	G	<INDEL>		SVTYPE=INDEL;END=124003;SVLEN=6
# 1	124004	g4	C	]1:124001]AAAATTTTC	EVENT=E2;PARID=g4;SVTYPE=BND	}

# CHROM	POS	ID	REF	ALT			INFO					CHROM	POS	ID	REF	ALT		INFO
# 1	123001	g1	T	]1:123200]T		EVENT=E1;PARID=g1;SVTYPE=BND	} ==>	1	123001	event	T	<DUP:TANDEM>	SVTYPE=DUP;END=123200;SVLEN=200
# 1	123200	g2	A	A[1:123001[		EVENT=E1;PARID=g2;SVTYPE=BND	}

# CHROM	POS	ID	REF	ALT			INFO					CHROM	POS	ID	REF	ALT		INFO
# 1	3001	g9	A       ]1:3200]CCCGGGA	   	EVENT=E3;PARID=g6;SVTYPE=BND	} ==>	1	100	event	A	<DUP:INS>	SVTYPE=DUP;END=3200;SVLEN=206
# 1	3200	g10	T       TCCCGGG[1:3001[		EVENT=E3;PARID=g6;SVTYPE=BND	}



# A DUP (duplication) looks like this.
# When calling a DUP, make sure that the ref.seq. is seen spanning the breakend points of the duplication.
# If the ref.seq. is not present, then it can't be a DUP, it must be some other complex variant.
# These normal duplication events are flagged as DUP:TANDEM or DUP:INS in the ALT field.
#
# ---------------------+----------------------------------------------------------+
#             16376518 | 16376519                                        16386259 +-----+
# ---------------------+----------------------------------------------------------+     |
#                                                                                       |
#                +----------------[may_or_may_not_have_inserted_sequence]---------------+
#                |           
#                |     +----------------------------------------------------------+----------------------------
#                +-----+ 16376319                                        16386259 | 16386260
#                      +----------------------------------------------------------+----------------------------
# 

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2019, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

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
def round_float_to_str( a_float ):
	if (a_float == int(a_float)):
		a_float = int(a_float)
	str_float = str(a_float)
	return str_float

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
		this_tranche2 = BND_records_fields[i]['TRANCHE2']
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
def max_tranche2( tranche2_in1, tranche2_in2 ):

	tranche2_final = ""
	if ((tranche2_in1 == "HIGH") or (tranche2_in2 == "HIGH")):
		tranche2_final = "HIGH"
	elif ((tranche2_in1 == "INTERMEDIATE") or (tranche2_in2 == "INTERMEDIATE")):
		tranche2_final = "INTERMEDIATE"
	elif ((tranche2_in1 == "LOW") or (tranche2_in2 == "LOW")):
		tranche2_final = "LOW"

	return tranche2_final

######################################################
def determine_VAFs( vaf_in1, vaf_in2 ):

	vaf_final = vaf_in1
	if (vaf_in2 > vaf_in1):
		vaf_final = vaf_in2
	vaf_final = str(vaf_final)
	list_of_vafs = str(vaf_in1) + ',' + str(vaf_in2)

	return vaf_final, list_of_vafs

######################################################
def extract_BND_values( index_left, index_right, BND_records_fields ):

	new_gridssAS = str(BND_records_fields[index_left]['AS']) + ',' + str(BND_records_fields[index_right]['AS'])
	new_gridssASQ = str(BND_records_fields[index_left]['ASQ']) + ',' + str(BND_records_fields[index_right]['ASQ'])
	new_gridssASRP = str(BND_records_fields[index_left]['ASRP']) + ',' + str(BND_records_fields[index_right]['ASRP'])
	new_gridssASRR = str(BND_records_fields[index_left]['ASRR']) + ',' + str(BND_records_fields[index_right]['ASRR'])
	new_gridssBA = str(BND_records_fields[index_left]['BA']) + ',' + str(BND_records_fields[index_right]['BA'])
	new_gridssBAQ = str(BND_records_fields[index_left]['BAQ']) + ',' + str(BND_records_fields[index_right]['BAQ'])
	new_gridssBEID = str(BND_records_fields[index_left]['BEID']) + ',' + str(BND_records_fields[index_right]['BEID'])
	new_gridssBQ = str(BND_records_fields[index_left]['BQ']) + ',' + str(BND_records_fields[index_right]['BQ'])
	new_gridssBSC = str(BND_records_fields[index_left]['BSC']) + ',' + str(BND_records_fields[index_right]['BSC'])
	new_gridssBSCQ = str(BND_records_fields[index_left]['BSCQ']) + ',' + str(BND_records_fields[index_right]['BSCQ'])
	new_gridssBUM = str(BND_records_fields[index_left]['BUM']) + ',' + str(BND_records_fields[index_right]['BUM'])
	new_gridssBUMQ = str(BND_records_fields[index_left]['BUMQ']) + ',' + str(BND_records_fields[index_right]['BUMQ'])
	new_gridssCAS = str(BND_records_fields[index_left]['CAS']) + ',' + str(BND_records_fields[index_right]['CAS'])
	new_gridssCASQ = str(BND_records_fields[index_left]['CASQ']) + ',' + str(BND_records_fields[index_right]['CASQ'])
	new_gridssCIEND = str(BND_records_fields[index_left]['CIEND']) + ',' + str(BND_records_fields[index_right]['CIEND'])
	new_gridssCIPOS = str(BND_records_fields[index_left]['CIPOS']) + ',' + str(BND_records_fields[index_right]['CIPOS'])
	new_gridssCIRPOS = str(BND_records_fields[index_left]['CIRPOS']) + ',' + str(BND_records_fields[index_right]['CIRPOS'])
	new_gridssCQ = str(BND_records_fields[index_left]['CQ']) + ',' + str(BND_records_fields[index_right]['CQ'])
	new_gridssFILTER = str(BND_records_fields[index_left]['FILTER']) + ',' + str(BND_records_fields[index_right]['FILTER'])
	new_gridssHOMLEN = str(BND_records_fields[index_left]['HOMLEN']) + ',' + str(BND_records_fields[index_right]['HOMLEN'])
	new_gridssHOMSEQ = str(BND_records_fields[index_left]['HOMSEQ']) + ',' + str(BND_records_fields[index_right]['HOMSEQ'])
	new_gridssIC = str(BND_records_fields[index_left]['IC']) + ',' + str(BND_records_fields[index_right]['IC'])
	new_gridssIHOMPOS = str(BND_records_fields[index_left]['IHOMPOS']) + ',' + str(BND_records_fields[index_right]['IHOMPOS'])
	new_gridssIQ = str(BND_records_fields[index_left]['IQ']) + ',' + str(BND_records_fields[index_right]['IQ'])
	new_gridssQUAL = str(BND_records_fields[index_left]['QUAL']) + ',' + str(BND_records_fields[index_right]['QUAL'])
	new_gridssRAS = str(BND_records_fields[index_left]['RAS']) + ',' + str(BND_records_fields[index_right]['RAS'])
	new_gridssRASQ = str(BND_records_fields[index_left]['RASQ']) + ',' + str(BND_records_fields[index_right]['RASQ'])
	new_gridssREF = str(BND_records_fields[index_left]['REFCOUNT']) + ',' + str(BND_records_fields[index_right]['REFCOUNT'])
	new_gridssREFPAIR = str(BND_records_fields[index_left]['REFPAIR']) + ',' + str(BND_records_fields[index_right]['REFPAIR'])
	new_gridssRP = str(BND_records_fields[index_left]['RP']) + ',' + str(BND_records_fields[index_right]['RP'])
	new_gridssRPQ = str(BND_records_fields[index_left]['RPQ']) + ',' + str(BND_records_fields[index_right]['RPQ'])
	new_gridssRSI = str(BND_records_fields[index_left]['RSI']) + ',' + str(BND_records_fields[index_right]['RSI'])
	new_gridssSI = str(BND_records_fields[index_left]['SI']) + ',' + str(BND_records_fields[index_right]['SI'])
	new_gridssSR = str(BND_records_fields[index_left]['SR']) + ',' + str(BND_records_fields[index_right]['SR'])
	new_gridssSRQ = str(BND_records_fields[index_left]['SRQ']) + ',' + str(BND_records_fields[index_right]['SRQ'])

	return new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ

######################################################
def format_new_fields_for_info_format_sample( output_gridss_info_fields, output_gridss_format_fields, new_pos, new_info_svtype, new_info_end, new_info_svlen, new_info_tranche, new_info_tranche2, new_cipos, new_ciend, new_VAF, new_GT,new_gridssVAF, new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ ):

	new_info = 'SVTYPE=' + new_info_svtype + ';END=' + str(new_info_end) + ';SVLEN=' + str(new_info_svlen) + ';TRANCHE=' + str(new_info_tranche) + ';TRANCHE2=' + str(new_info_tranche2) + ';CIPOS=' + str(new_cipos) + ';CIEND=' + str(new_ciend) + ';VAF=' + str(new_VAF)
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

	new_format = 'GT:POS:END:SVLEN:TRANCHE:TRANCHE2:CIPOS:CIEND:VAF'
	if (output_gridss_format_fields):
		new_format = new_format + ':gridssVAF:gridssAS:gridssASQ:gridssASRP:gridssASRR:gridssBA:gridssBAQ:gridssBEID:gridssBQ:gridssBSC:gridssBSCQ'
		new_format = new_format + ':gridssBUM:gridssBUMQ:gridssCAS:gridssCASQ:gridssCIEND:gridssCIPOS:gridssCIRPOS:gridssCQ'
		new_format = new_format + ':gridssFILTER:gridssHOMLEN:gridssHOMSEQ:gridssIC:gridssIHOMPOS:gridssIQ:gridssQUAL:gridssRAS:gridssRASQ'
		new_format = new_format + ':gridssREF:gridssREFPAIR:gridssRP:gridssRPQ:gridssRSI:gridssSI:gridssSR:gridssSRQ'

	new_sample = new_GT + ':' + str(new_pos) + ':' + str(new_info_end) + ':' + str(new_info_svlen) + ':' + str(new_info_tranche) + ':' + str(new_info_tranche2) + ':' + str(new_cipos) + ':' + str(new_ciend) + str(new_VAF)
	if (output_gridss_format_fields):
		new_sample = new_sample + ':' + str(new_gridssVAF) + ':' + str(new_gridssAS) + ':' + str(new_gridssASQ) + ':' + str(new_gridssASRP) + ':' + str(new_gridssASRR) + ':' + str(new_gridssBA) + ':' + str(new_gridssBAQ) + ':' + str(new_gridssBEID) + ':' + str(new_gridssBQ) + ':' + str(new_gridssBSC) + ':' + str(new_gridssBSCQ)
		new_sample = new_sample + ':' + str(new_gridssBUM) + ':' + str(new_gridssBUMQ) + ':' + str(new_gridssCAS) + ':' + str(new_gridssCASQ) + ':' + str(new_gridssCIEND) + ':' + str(new_gridssCIPOS) + ':' + str(new_gridssCIRPOS) + ':' + str(new_gridssCQ)
		new_sample = new_sample + ':' + str(new_gridssFILTER) + ':' + str(new_gridssHOMLEN) + ':' + str(new_gridssHOMSEQ) + ':' + str(new_gridssIC) + ':' + str(new_gridssIHOMPOS) + ':' + str(new_gridssIQ) + ':' + str(new_gridssQUAL) + ':' + str(new_gridssRAS) + ':' + str(new_gridssRASQ)
		new_sample = new_sample + ':' + str(new_gridssREF) + ':' + str(new_gridssREFPAIR) + ':' + str(new_gridssRP) + ':' + str(new_gridssRPQ) + ':' + str(new_gridssRSI) + ':' + str(new_gridssSI) + ':' + str(new_gridssSR) + ':' + str(new_gridssSRQ)

	return new_info, new_format, new_sample

######################################################
def write_one_summary_VCF_record( BND_records, event_id, output_gridss_info_fields, output_gridss_format_fields ):

	# A duplication has two breakends and at both breakends, there should be both wild-type ref.seq. and split-read.
	# If the ref.seq. is not present at one of the breakends, then this is not a classic duplication.
	# Instead it is probably a duplicated sequence that has been inserted elsewhere.
	# The number of ref.seq. spanning the breakend is REFCOUNT. The number of split-reads supporting the breakend is REFPAIR.
	# This constant is the cut-off point for deciding if ref.seq. has been detected or not, to determine if it is a normal duplication event.
	CONST_refseq_to_splitread_ratio_too_low_for_normal_DUP = 0.1

	count_left_N = 0
	count_left_NNN = 0
	count_N_left = 0
	count_NNN_left = 0
	count_right_N = 0
	count_right_NNN = 0
	count_N_right = 0
	count_NNN_right = 0

	index_left_N = -1
	index_left_NNN = -1
	index_N_left = -1
	index_NNN_left = -1
	index_right_N = -1
	index_right_NNN = -1
	index_N_right = -1
	index_NNN_right = -1
	index2_N_left = -1
	index2_right_N = -1

	re_left_N = re.compile('\][a-zA-Z0-9.]+\:[0-9]+\][ACGTN]$') 		# eg. ]1:123200]T
	re_left_NNN = re.compile('\][a-zA-Z0-9.]+\:[0-9]+\][ACGTN]{2,}$')	# eg. ]1:123200]TA, ]1:123200]TAC, ]1:123200]TACC, etc.
	re_N_left = re.compile('[ACGTN]\][a-zA-Z0-9.]+\:[0-9]+\]$') 		# eg. T]1:123200]
	re_NNN_left = re.compile('[ACGTN]{2,}\][a-zA-Z0-9.]+\:[0-9]+\]$')	# eg. TA]1:123200], TAC]1:123200], TACC]1:123200], etc.
	re_right_N = re.compile('\[[a-zA-Z0-9.]+\:[0-9]+\[[ACGTN]$') 		# eg. [1:123200[T
	re_right_NNN = re.compile('\[[a-zA-Z0-9.]+\:[0-9]+\[[ACGTN]{2,}$')	# eg. [1:123200[TA, [1:123200[TAC, [1:123200[TACC, etc.
	re_N_right = re.compile('[ACGTN]\[[a-zA-Z0-9.]+\:[0-9]+\[$') 		# eg. T[1:123200[
	re_NNN_right = re.compile('[ACGTN]{2,}\[[a-zA-Z0-9.]+\:[0-9]+\[$')	# eg. TA[1:123200[, TAC[1:123200[, TACC[1:123200[, etc.

	BND_records_fields = []

	for i in range( 0, len(BND_records) ):
		inline = BND_records[i]
		infields = inline.split("\t")
		chrom = str(infields[0])
		pos = int(infields[1])
		snpid = str(infields[2])
		ref = str(infields[3])
		alt = str(infields[4])
		qual = str(infields[5])
		vcffilter = str(infields[6])
		info = str(infields[7])
		snpformat = str(infields[8])
		sample = str(infields[9])
		if (len(infields) > 10):
			sys.stderr.write('FOR EVENT ' + event_id + ' THERE IS MORE THAN ONE SAMPLE AND THESE SUBSEQUENT EXAMPLES WILL BE IGNORED' + "\n")
		this_BND_record = {}
		this_BND_record['CHROM'] = chrom
		this_BND_record['POS'] = pos
		this_BND_record['ID'] = snpid
		this_BND_record['REF'] = ref
		this_BND_record['ALT'] = alt
		this_BND_record['QUAL'] = qual
		this_BND_record['FILTER'] = vcffilter
		this_BND_record['INFO'] = info
		this_BND_record['FORMAT'] = snpformat
		this_BND_record['sample'] = sample
		this_bndvaf, this_ciend, this_cipos, this_cirpos, this_tranche2, this_ihompos, this_as, this_asq, this_asrp, this_asrr, this_ba, this_baq, this_beid, this_bq, this_bsc, this_bscq, this_bum, this_bumq, this_cas, this_casq, this_cq, this_homlen, this_homseq, this_ic, this_ihompos, this_iq, this_ras, this_rasq, this_ref_count, this_refpair, this_rp, this_rpq, this_rsi, this_si, this_sr, this_srq = extract_fields_from_INFO( info )
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

		BND_records_fields.append( this_BND_record )
		if (re_left_N.match( alt ) is not None):
			count_left_N = count_left_N + 1
			index_left_N = i
		elif (re_left_NNN.match( alt ) is not None):
			count_left_NNN = count_left_NNN + 1
			index_left_NNN = i
		elif (re_N_left.match( alt ) is not None):
			count_N_left = count_N_left + 1
			if (index_N_left == -1):
				index_N_left = i
			else:
				index2_N_left = i
		elif (re_NNN_left.match( alt ) is not None):
			count_NNN_left = count_NNN_left + 1
			index_NNN_left = i
		elif (re_right_N.match( alt ) is not None):
			count_right_N = count_right_N + 1
			if (index_right_N == -1):
				index_right_N = i
			else:
				index2_right_N = i
		elif (re_right_NNN.match( alt ) is not None):
			count_right_NNN = count_right_NNN + 1
			index_right_NNN = i
		elif (re_N_right.match( alt ) is not None):
			count_N_right = count_N_right + 1
			index_N_right = i
		elif (re_NNN_right.match( alt ) is not None):
			count_NNN_right = count_NNN_right + 1
			index_NNN_right = i
		#else:
		#	sys.stderr.write('FOR EVENT ' + event_id + ' A BND RECORD ALT FIELD DIDN\'T MATCH EXPECTED PATTERNS: ' + inline + "\n")

	number_of_SV_patterns_this_BND_group_matches = 0
	this_group_is_discarded = False

	# Ignore multi-chromosome arrangments for now.
	# Don't process them. However if user has asked to save them, then do save them.

	chromosomes_are_all_the_same = True
	first_chrom = BND_records_fields[0]['CHROM']
	if (len(BND_records_fields) > 1):
		for i in range( 1, len(BND_records_fields) ) :
			next_chrom = this_BND_record['CHROM']
			if (next_chrom != first_chrom):
				chromosomes_are_all_the_same = False

	if (chromosomes_are_all_the_same == False):
		this_group_is_discarded = True
		for this_record in BND_records:
			outline = this_record + "\n"
			sys.stdout.write( outline )

	else: # chromosomes_are_all_the_same == True

		# Identify tandem duplications and deletions

		# DUP:TANDEM	CHROM	POS	ALT
		# 		1	100	]1:200]N
		#		1	200	N[1:100[

		# DEL:		CHROM	POS	ALT
		# 		1	100	N[1:200[
		#		1	200	]1:100]N

		if ((count_left_N == 1) and (count_N_right == 1) and
			(count_left_NNN == 0) and (count_N_left == 0) and (count_NNN_left == 0) and (count_right_N == 0) and (count_right_NNN == 0) and (count_NNN_right == 0)):

			if (BND_records_fields[index_left_N]['CHROM'] == BND_records_fields[index_N_right]['CHROM']):

				if (BND_records_fields[index_left_N]['POS'] < BND_records_fields[index_N_right]['POS']):

					pos_of_left_N = int(BND_records_fields[index_left_N]['POS'])
					bits = BND_records_fields[index_left_N]['ALT'].split(']')
					bits1 = bits[1].split(':')
					pos_of_alt_of_left_N = int(bits1[1])
					pos_of_N_right = int(BND_records_fields[index_N_right]['POS'])
					bits = BND_records_fields[index_N_right]['ALT'].split('[')
					bits1 = bits[1].split(':')
					pos_of_alt_of_N_right = int(bits1[1])

					if (    ((pos_of_left_N == pos_of_alt_of_N_right) and (pos_of_N_right == pos_of_alt_of_left_N)) or
						(pos_of_left_N == pos_of_alt_of_N_right) or (pos_of_N_right == pos_of_alt_of_left_N)    ):

						# This pair of BND records is describing a tandem duplication    or
						# This pair of BND records is probably describing a tandem duplication
						# and one of the POS is not the same as the other record's ALT due to mapping problems caused by low complexity region
						# or due to there being more than one tandem duplication of the same sequence

						new_chrom = BND_records_fields[index_left_N]['CHROM']
						new_pos = BND_records_fields[index_left_N]['POS']
						new_ref = BND_records_fields[index_left_N]['REF']
						new_snpid = event_id
						new_alt = '<DUP:TANDEM>'
						new_info_svtype = 'DUP'
						new_info_end = BND_records_fields[index_N_right]['POS']
						new_info_svlen = new_info_end - new_pos
						new_info_tranche, new_qual, new_filter = determine_tranche_field( BND_records_fields )
						new_cipos = BND_records_fields[index_left_N]['CIPOS']
						new_ciend = BND_records_fields[index_N_right]['CIRPOS']
						new_info_tranche2 = max_tranche2( BND_records_fields[index_left_N]['TRANCHE2'], BND_records_fields[index_N_right]['TRANCHE2'] )
						new_GT = '0/1'
						if ((BND_records_fields[index_left_N]['BNDVAF'] > 0.65) or (BND_records_fields[index_N_right]['BNDVAF'] > 0.65)):
							new_GT = '1/1'
						new_VAF, new_gridssVAF = determine_VAFs( BND_records_fields[index_left_N]['BNDVAF'], BND_records_fields[index_N_right]['BNDVAF'] )
						index_left = index_left_N
						index_right = index_N_right
						new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ = extract_BND_values( index_left, index_right, BND_records_fields )
						new_info, new_format, new_sample = format_new_fields_for_info_format_sample( output_gridss_info_fields, output_gridss_format_fields, new_pos, new_info_svtype, new_info_end, new_info_svlen, new_info_tranche, new_info_tranche2, new_cipos, new_ciend, new_VAF, new_GT, new_gridssVAF, new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ )
						outline = str(new_chrom) + "\t" + str(new_pos) + "\t" + str(event_id) + "\t" + new_ref + "\t" + new_alt + "\t" + round_float_to_str(new_qual) + "\t" + str(new_filter) + "\t" + new_info + "\t" + new_format + "\t" + new_sample + "\n"
						sys.stdout.write( outline )
						number_of_SV_patterns_this_BND_group_matches = number_of_SV_patterns_this_BND_group_matches + 1

					#else: # (pos_of_left_N != pos_of_alt_of_N_right) and (pos_of_N_right != pos_of_alt_of_left_N)
					#	sys.stderr.write('FOR EVENT ' + event_id + ' left_N.POS not= N_right.ALT OR left_N.ALT not= N_right.POS' + "\n")

				elif (BND_records_fields[index_left_N]['POS'] > BND_records_fields[index_N_right]['POS']):

					# This pair of BND records is describing a deletion

					pos_of_left_N = int(BND_records_fields[index_left_N]['POS'])
					bits = BND_records_fields[index_left_N]['ALT'].split(']')
					bits1 = bits[1].split(':')
					pos_of_alt_of_left_N = int(bits1[1])
					pos_of_N_right = int(BND_records_fields[index_N_right]['POS'])
					bits = BND_records_fields[index_N_right]['ALT'].split('[')
					bits1 = bits[1].split(':')
					pos_of_alt_of_N_right = int(bits1[1])

					if ((pos_of_left_N == pos_of_alt_of_N_right) and (pos_of_N_right == pos_of_alt_of_left_N)):

						new_chrom = BND_records_fields[index_left_N]['CHROM']
						new_pos = BND_records_fields[index_N_right]['POS']
						new_ref = BND_records_fields[index_N_right]['REF']
						new_snpid = event_id
						new_alt = '<DEL>'
						new_info_svtype = 'DEL'
						new_info_end = BND_records_fields[index_left_N]['POS'] - 1
						new_info_svlen = new_pos - new_info_end
						new_info_tranche, new_qual, new_filter = determine_tranche_field( BND_records_fields )
						new_cipos = BND_records_fields[index_N_right]['CIPOS']
						new_ciend = BND_records_fields[index_left_N]['CIRPOS']
						new_info_tranche2 = max_tranche2( BND_records_fields[index_left_N]['TRANCHE2'], BND_records_fields[index_N_right]['TRANCHE2'] )
						new_GT = '0/1'
						if ((BND_records_fields[index_left_N]['BNDVAF'] > 0.85) or (BND_records_fields[index_N_right]['BNDVAF'] > 0.85)):
							new_GT = '1/1'
						new_VAF, new_gridssVAF = determine_VAFs( BND_records_fields[index_left_N]['BNDVAF'], BND_records_fields[index_N_right]['BNDVAF'] )
						index_left = index_left_N
						index_right = index_N_right
						new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ = extract_BND_values( index_left, index_right, BND_records_fields )
						new_info, new_format, new_sample = format_new_fields_for_info_format_sample( output_gridss_info_fields, output_gridss_format_fields, new_pos, new_info_svtype, new_info_end, new_info_svlen, new_info_tranche, new_info_tranche2, new_cipos, new_ciend, new_VAF, new_GT, new_gridssVAF, new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ )
						outline = str(new_chrom) + "\t" + str(new_pos) + "\t" + str(event_id) + "\t" + new_ref + "\t" + new_alt + "\t" + round_float_to_str(new_qual) + "\t" + str(new_filter) + "\t" + new_info + "\t" + new_format + "\t" + new_sample + "\n"
						sys.stdout.write( outline )
						number_of_SV_patterns_this_BND_group_matches = number_of_SV_patterns_this_BND_group_matches + 1

					else:
						# sys.stderr.write('FOR EVENT ' + event_id + ' left_N.POS not= N_right.ALT OR left_N.ALT not= N_right.POS' + "\n")
						# Don't output an error. They were inspected and were found to be mappings to low complexity regions, with no apparent structural variants.
						this_group_is_discarded = True
						for this_record in BND_records:
							outline = this_record + "\n"
							sys.stdout.write( outline )
				#else:
				#	sys.stderr.write('FOR EVENT ' + event_id + ' left_N and N_right UNEXPECTEDLY HAVE SAME POS' + "\n")
			#else:
			#	sys.stderr.write('FOR EVENT ' + event_id + ' left_N and N_right DON\'T HAVE THE SAME CHROMOSOME' + "\n")

		# Identify insertions and indels, and also duplications that also have an inserted sequence

		# INS:		CHROM	POS	ALT
		# 		1	100	]1:100]NNN
		#		1	100	NNN[1:100[

		# INS:		CHROM	POS	ALT
		# 		1	100	]1:101]NNN
		#		1	101	NNN[1:100[

		# INDEL:	CHROM	POS	ALT
		# 		1	100	NNN[1:200[
		#		1	200	]1:100]NNN

		# DUP:INS:	CHROM	POS	ALT
		# 		1	100	]1:200]NNN
		#		1	200	NNN[1:100[

		if ((count_left_NNN == 1) and (count_NNN_right == 1) and
			(count_left_N == 0) and (count_N_left == 0) and (count_NNN_left == 0) and (count_right_N == 0) and (count_right_NNN == 0)and (count_N_right == 0)):

			if (BND_records_fields[index_left_NNN]['CHROM'] == BND_records_fields[index_NNN_right]['CHROM']):

				if ( BND_records_fields[index_left_NNN]['POS'] == (BND_records_fields[index_NNN_right]['POS']+1) ):

					pos_of_left_NNN = int(BND_records_fields[index_left_NNN]['POS'])
					bits = BND_records_fields[index_left_NNN]['ALT'].split(']')
					bits1 = bits[1].split(':')
					pos_of_alt_of_left_NNN = int(bits1[1])
					pos_of_NNN_right = int(BND_records_fields[index_NNN_right]['POS'])
					bits = BND_records_fields[index_NNN_right]['ALT'].split('[')
					bits1 = bits[1].split(':')
					pos_of_alt_of_NNN_right = int(bits1[1])

					if ((pos_of_left_NNN == pos_of_alt_of_NNN_right) and (pos_of_NNN_right == pos_of_alt_of_left_NNN)):

						# This pair of BND records is describing an insertion

						new_chrom = BND_records_fields[index_NNN_right]['CHROM']
						new_pos = BND_records_fields[index_NNN_right]['POS']
						new_info_end = new_pos
						new_ref = BND_records_fields[index_NNN_right]['REF']
						new_snpid = event_id
						new_alt = '<INS>'
						new_info_svtype = 'INS'
						bits = BND_records_fields[index_NNN_right]['ALT'].split('[')
						ALT_sequence = bits[0]
						new_info_insseq = ALT_sequence[1:]
						new_info_svlen = len(new_info_insseq)
						bits2 = BND_records_fields[index_left_NNN]['ALT'].split(']')
						ALT2_sequence = bits2[ len(bits2)-1 ]
						#if ( len(ALT_sequence) != len(ALT2_sequence) ):
						#	sys.stderr.write('FOR EVENT ' + event_id + ' left_NNN and NNN_right DON\'T HAVE THE SAME LENGTH OF ALT SEQUENCES (INS WHERE left_NNN.POS = NNN_right.POS + 1)' + "\n")
						new_info_tranche, new_qual, new_filter = determine_tranche_field( BND_records_fields )
						new_cipos = BND_records_fields[index_NNN_right]['CIPOS']
						new_ciend = BND_records_fields[index_NNN_right]['CIRPOS']
						new_info_tranche2 = max_tranche2( BND_records_fields[index_left_NNN]['TRANCHE2'], BND_records_fields[index_NNN_right]['TRANCHE2'] )
						new_GT = '0/1'
						if ((BND_records_fields[index_left_NNN]['BNDVAF'] > 0.85) or (BND_records_fields[index_NNN_right]['BNDVAF'] > 0.85)):
							new_GT = '1/1'
						new_VAF, new_gridssVAF = determine_VAFs( BND_records_fields[index_left_NNN]['BNDVAF'], BND_records_fields[index_NNN_right]['BNDVAF'] )

						index_left = index_left_NNN
						index_right = index_NNN_right
						new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ = extract_BND_values( index_left, index_right, BND_records_fields )
						new_info, new_format, new_sample = format_new_fields_for_info_format_sample( output_gridss_info_fields, output_gridss_format_fields, new_pos, new_info_svtype, new_info_end, new_info_svlen, new_info_tranche, new_info_tranche2, new_cipos, new_ciend, new_VAF, new_GT, new_gridssVAF, new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ )
						new_info = new_info + ';INSSEQ=' + new_info_insseq
						new_format = new_format + ':INSSEQ'
						new_sample = new_sample + ':' + new_info_insseq
						outline = str(new_chrom) + "\t" + str(new_pos) + "\t" + str(event_id) + "\t" + new_ref + "\t" + new_alt + "\t" + round_float_to_str(new_qual) + "\t" + str(new_filter) + "\t" + new_info + "\t" + new_format + "\t" + new_sample + "\n"
						sys.stdout.write( outline )
						number_of_SV_patterns_this_BND_group_matches = number_of_SV_patterns_this_BND_group_matches + 1

					#else:
					#	sys.stderr.write('FOR EVENT ' + event_id + ' left_NNN.POS not= NNN_right.ALT OR left_NNN.ALT not= NNN_right.POS' + "\n")

				elif (BND_records_fields[index_left_NNN]['POS'] == BND_records_fields[index_NNN_right]['POS']):

					pos_of_left_NNN = int(BND_records_fields[index_left_NNN]['POS'])
					bits = BND_records_fields[index_left_NNN]['ALT'].split(']')
					bits1 = bits[1].split(':')
					pos_of_alt_of_left_NNN = int(bits1[1])
					pos_of_NNN_right = int(BND_records_fields[index_NNN_right]['POS'])
					bits = BND_records_fields[index_NNN_right]['ALT'].split('[')
					bits1 = bits[1].split(':')
					pos_of_alt_of_NNN_right = int(bits1[1])

					if ((pos_of_left_NNN == pos_of_alt_of_NNN_right) and (pos_of_NNN_right == pos_of_alt_of_left_NNN)):

						# This pair of BND records is describing an insertion

						new_chrom = BND_records_fields[index_left_NNN]['CHROM']
						new_pos = BND_records_fields[index_left_NNN]['POS']
						new_info_end = new_pos
						new_ref = BND_records_fields[index_left_NNN]['REF']
						new_snpid = event_id
						new_alt = '<INS>'
						new_info_svtype = 'INS'
						bits = BND_records_fields[index_left_NNN]['ALT'].split(']')
						ALT_sequence = bits[ len(bits)-1 ]
						bits2 = BND_records_fields[index_NNN_right]['ALT'].split('[')
						ALT2_sequence = bits2[0]
						#if ( len(ALT_sequence) != len(ALT2_sequence) ):
						#	sys.stderr.write('FOR EVENT ' + event_id + ' left_NNN and NNN_right DON\'T HAVE THE SAME LENGTH OF ALT SEQUENCES (INS WHERE left_NNN.POS <= NNN_right.POS)' + "\n")
						new_info_tranche, new_qual, new_filter = determine_tranche_field( BND_records_fields )
						bits = BND_records_fields[index_NNN_right]['ALT'].split('[')
						new_info_insseq = bits[0]
						new_info_insseq = new_info_insseq[1:] + new_ref
						new_info_svlen = len(new_info_insseq)
						new_cipos = BND_records_fields[index_left_NNN]['CIPOS']
						new_ciend = BND_records_fields[index_left_NNN]['CIRPOS']
						new_info_tranche2 = max_tranche2( BND_records_fields[index_left_NNN]['TRANCHE2'], BND_records_fields[index_NNN_right]['TRANCHE2'] )
						new_GT = '0/1'
						if ((BND_records_fields[index_left_NNN]['BNDVAF'] > 0.85) or (BND_records_fields[index_NNN_right]['BNDVAF'] > 0.85)):
							new_GT = '1/1'
						new_VAF, new_gridssVAF = determine_VAFs( BND_records_fields[index_left_NNN]['BNDVAF'], BND_records_fields[index_NNN_right]['BNDVAF'] )
						index_left = index_left_NNN
						index_right = index_NNN_right
						new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ = extract_BND_values( index_left, index_right, BND_records_fields )
						new_info, new_format, new_sample = format_new_fields_for_info_format_sample( output_gridss_info_fields, output_gridss_format_fields, new_pos, new_info_svtype, new_info_end, new_info_svlen, new_info_tranche, new_info_tranche2, new_cipos, new_ciend, new_VAF, new_GT, new_gridssVAF, new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ )
						new_info = new_info + ';INSSEQ=' + new_info_insseq
						new_format = new_format + ':INSSEQ'
						new_sample = new_sample + ':' + new_info_insseq
						outline = str(new_chrom) + "\t" + str(new_pos) + "\t" + str(event_id) + "\t" + new_ref + "\t" + new_alt + "\t" + round_float_to_str(new_qual) + "\t" + str(new_filter) + "\t" + new_info + "\t" + new_format + "\t" + new_sample + "\n"
						sys.stdout.write( outline )
						number_of_SV_patterns_this_BND_group_matches = number_of_SV_patterns_this_BND_group_matches + 1

					#else: # (pos_of_left_NNN != pos_of_alt_of_NNN_right) or (pos_of_NNN_right != pos_of_alt_of_left_NNN)
					#	sys.stderr.write('FOR EVENT ' + event_id + ' left_NNN.POS = NNN_right.POS BUT left_NNN.POS not= NNN_right.ALT OR left_NNN.ALT not= NNN_right.POS' + "\n")

				elif (BND_records_fields[index_left_NNN]['POS'] < BND_records_fields[index_NNN_right]['POS']):

					pos_of_left_NNN = int(BND_records_fields[index_left_NNN]['POS'])
					bits = BND_records_fields[index_left_NNN]['ALT'].split(']')
					bits1 = bits[1].split(':')
					pos_of_alt_of_left_NNN = int(bits1[1])
					pos_of_NNN_right = int(BND_records_fields[index_NNN_right]['POS'])
					bits = BND_records_fields[index_NNN_right]['ALT'].split('[')
					bits1 = bits[1].split(':')
					pos_of_alt_of_NNN_right = int(bits1[1])

					if ((pos_of_left_NNN == pos_of_alt_of_NNN_right) and (pos_of_NNN_right == pos_of_alt_of_left_NNN)):

						# This pair of BND records is describing a duplication that also has an inserted sequence.

						new_chrom = BND_records_fields[index_left_NNN]['CHROM']
						new_pos = BND_records_fields[index_left_NNN]['POS']
						new_ref = BND_records_fields[index_left_NNN]['REF']
						new_snpid = event_id
						new_alt = '<DUP:INS>'
						new_info_svtype = 'DUP'
						new_info_end = BND_records_fields[index_NNN_right]['POS']
						length_of_insert_in_POS = new_info_end - new_pos
						bits = BND_records_fields[index_left_NNN]['ALT'].split(']')
						ALT_sequence = bits[ len(bits)-1 ]
						length_of_insert_in_ALT = len(ALT_sequence) - 1
						bits2 = BND_records_fields[index_NNN_right]['ALT'].split('[')
						ALT2_sequence = bits2[0]
						#if ( len(ALT_sequence) != len(ALT2_sequence) ):
						#	sys.stderr.write('FOR EVENT ' + event_id + ' left_NNN and NNN_right DON\'T HAVE THE SAME LENGTH OF ALT SEQUENCES (DUP:INS WHERE left_NNN.POS <= NNN_right.POS)' + "\n")
						new_info_svlen = length_of_insert_in_ALT + length_of_insert_in_POS
						new_info_tranche, new_qual, new_filter = determine_tranche_field( BND_records_fields )
						bits = BND_records_fields[index_NNN_right]['ALT'].split('[')
						new_info_insseq = bits[0]
						new_info_insseq = new_info_insseq[1:]
						new_cipos = BND_records_fields[index_left_NNN]['CIPOS']
						new_ciend = BND_records_fields[index_left_NNN]['CIRPOS']
						new_info_tranche2 = max_tranche2( BND_records_fields[index_left_NNN]['TRANCHE2'], BND_records_fields[index_NNN_right]['TRANCHE2'] )
						new_GT = '0/1'
						if ((BND_records_fields[index_left_NNN]['BNDVAF'] > 0.65) or (BND_records_fields[index_NNN_right]['BNDVAF'] > 0.65)):
							new_GT = '1/1'
						new_VAF, new_gridssVAF = determine_VAFs( BND_records_fields[index_left_NNN]['BNDVAF'], BND_records_fields[index_NNN_right]['BNDVAF'] )
						index_left = index_left_NNN
						index_right = index_NNN_right
						new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ = extract_BND_values( index_left, index_right, BND_records_fields )
						new_info, new_format, new_sample = format_new_fields_for_info_format_sample( output_gridss_info_fields, output_gridss_format_fields, new_pos, new_info_svtype, new_info_end, new_info_svlen, new_info_tranche, new_info_tranche2, new_cipos, new_ciend, new_VAF, new_GT, new_gridssVAF, new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ )
						new_info = new_info + ';INSSEQ=' + new_info_insseq
						new_format = new_format + ':INSSEQ'
						new_sample = new_sample + ':' + new_info_insseq
						outline = str(new_chrom) + "\t" + str(new_pos) + "\t" + str(event_id) + "\t" + new_ref + "\t" + new_alt + "\t" + round_float_to_str(new_qual) + "\t" + str(new_filter) + "\t" + new_info + "\t" + new_format + "\t" + new_sample + "\n"
						sys.stdout.write( outline )
						number_of_SV_patterns_this_BND_group_matches = number_of_SV_patterns_this_BND_group_matches + 1

					#else: # not a duplication that also has an inserted sequence

					#	sys.stderr.write('FOR EVENT ' + event_id + ' left_NNN.POS < NNN_right.POS BUT left_NNN.POS not= NNN_right.ALT OR left_NNN.ALT not= NNN_right.POS' + "\n")

				else: # (BND_records_fields[index_left_NNN]['POS'] > BND_records_fields[index_NNN_right]['POS']+1):

					# This pair of BND records is describing an indel, that is, a deletion that has an insertion inserted into it
					# The inserted sequence may be larger or smaller than the deleted sequence

					pos_of_NNN_right = int(BND_records_fields[index_NNN_right]['POS'])
					bits = BND_records_fields[index_NNN_right]['ALT'].split('[')
					bits1 = bits[1].split(':')
					pos_of_alt_of_NNN_right = int(bits1[1])
					pos_of_left_NNN = int(BND_records_fields[index_left_NNN]['POS'])
					bits = BND_records_fields[index_left_NNN]['ALT'].split(']')
					bits1 = bits[1].split(':')
					pos_of_alt_of_left_NNN = int(bits1[1])

					if ((pos_of_left_NNN == pos_of_alt_of_NNN_right) and (pos_of_NNN_right == pos_of_alt_of_left_NNN)):

						new_chrom = BND_records_fields[index_NNN_right]['CHROM']
						new_pos = BND_records_fields[index_NNN_right]['POS']
						new_ref = BND_records_fields[index_NNN_right]['REF']
						new_snpid = event_id
						new_alt = '<INDEL>'
						new_info_svtype = 'INDEL'
						new_info_end = BND_records_fields[index_left_NNN]['POS'] - 1
						length_of_insert_in_POS = new_pos - new_info_end
						bits = BND_records_fields[index_NNN_right]['ALT'].split('[')
						ALT_sequence = bits[0]
						length_of_insert_in_ALT = len(ALT_sequence) - 1
						bits2 = BND_records_fields[index_left_NNN]['ALT'].split(']')
						ALT2_sequence = bits2[ len(bits)-1 ]
						#if ( len(ALT_sequence) != len(ALT2_sequence) ):
						#	sys.stderr.write('FOR EVENT ' + event_id + ' NNN_right and left_NNN DON\'T HAVE THE SAME LENGTH OF ALT SEQUENCES' + "\n")
						new_info_svlen = length_of_insert_in_ALT + length_of_insert_in_POS
						new_info_tranche, new_qual, new_filter = determine_tranche_field( BND_records_fields )
						bits = BND_records_fields[index_NNN_right]['ALT'].split('[')
						new_info_insseq = bits[0]
						new_info_insseq = new_info_insseq[1:]
						new_cipos = BND_records_fields[index_NNN_right]['CIPOS']
						new_ciend = BND_records_fields[index_left_NNN]['CIRPOS']
						new_info_tranche2 = max_tranche2( BND_records_fields[index_left_NNN]['TRANCHE2'], BND_records_fields[index_NNN_right]['TRANCHE2'] )
						new_GT = '0/1'
						if ((BND_records_fields[index_left_NNN]['BNDVAF'] > 0.85) or (BND_records_fields[index_NNN_right]['BNDVAF'] > 0.85)):
							new_GT = '1/1'
						new_VAF, new_gridssVAF = determine_VAFs( BND_records_fields[index_left_NNN]['BNDVAF'], BND_records_fields[index_NNN_right]['BNDVAF'] )
						index_left = index_left_NNN
						index_right = index_NNN_right
						new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ = extract_BND_values( index_left, index_right, BND_records_fields )
						new_info, new_format, new_sample = format_new_fields_for_info_format_sample( output_gridss_info_fields, output_gridss_format_fields, new_pos, new_info_svtype, new_info_end, new_info_svlen, new_info_tranche, new_info_tranche2, new_cipos, new_ciend, new_VAF, new_GT, new_gridssVAF, new_gridssAS, new_gridssASQ, new_gridssASRP, new_gridssASRR, new_gridssBA, new_gridssBAQ, new_gridssBEID, new_gridssBQ, new_gridssBSC, new_gridssBSCQ, new_gridssBUM, new_gridssBUMQ, new_gridssCAS, new_gridssCASQ, new_gridssCIEND, new_gridssCIPOS, new_gridssCIRPOS, new_gridssCQ, new_gridssFILTER, new_gridssHOMLEN, new_gridssHOMSEQ, new_gridssIC, new_gridssIHOMPOS, new_gridssIQ, new_gridssQUAL, new_gridssRAS, new_gridssRASQ, new_gridssREF, new_gridssREFPAIR, new_gridssRP, new_gridssRPQ, new_gridssRSI, new_gridssSI, new_gridssSR, new_gridssSRQ )
						new_info = new_info + ';INSSEQ=' + new_info_insseq
						new_format = new_format + ':INSSEQ'
						new_sample = new_sample + ':' + new_info_insseq
						outline = str(new_chrom) + "\t" + str(new_pos) + "\t" + str(event_id) + "\t" + new_ref + "\t" + new_alt + "\t" + round_float_to_str(new_qual) + "\t" + str(new_filter) + "\t" + new_info + "\t" + new_format + "\t" + new_sample + "\n"
						sys.stdout.write( outline )
						number_of_SV_patterns_this_BND_group_matches = number_of_SV_patterns_this_BND_group_matches + 1

					#else:
					#	sys.stderr.write('FOR EVENT ' + event_id + ' left_NNN.POS not= NNN_right.ALT OR left_NNN.ALT not= NNN_right.POS' + "\n")
			#else:
			#	sys.stderr.write('FOR EVENT ' + event_id + ' left_NNN and NNN_right DON\'T HAVE THE SAME CHROMOSOME' + "\n")

		# Identify various cases that will not be dealt with for now

		# HALF OF INV:	CHROM	POS	ALT
		# 		1	101	N]1:10001]
		#		1	10002	[1:102[N

		# HALF OF INV:	CHROM	POS	ALT
		#		1	102	[1:10002[N
		#		1	10001	N]1:101]

		if ((count_N_left == 1) and (count_right_N == 1) and
			(count_left_N == 0) and (count_left_NNN == 0) and (count_NNN_left == 0) and (count_right_NNN == 0) and (count_N_right == 0) and (count_NNN_right == 0)):
			this_group_is_discarded = True
			for this_record in BND_records:
				outline = this_record + "\n"
				sys.stdout.write( outline )

		# CASE:		CHROM	POS	ALT
		# 		1	100	[1:200[N
		#		1	200	[1:100[N
		if ((count_right_N == 2) and
			(count_left_N == 0) and (count_left_NNN == 0) and (count_N_left == 0) and (count_NNN_left == 0) and (count_right_NNN == 0) and (count_N_right == 0) and (count_NNN_right == 0)):
			this_group_is_discarded = True
			for this_record in BND_records:
				outline = this_record + "\n"
				sys.stdout.write( outline )

		# CASE:		CHROM	POS	ALT
		# 		1	100	[1:200[NNN
		#		1	200	[1:100[NNN
		if ((count_right_NNN == 2) and
			(count_left_N == 0) and (count_left_NNN == 0) and (count_N_left == 0) and (count_NNN_left == 0) and (count_right_N == 0) and (count_N_right == 0) and (count_NNN_right == 0)):
			this_group_is_discarded = True
			for this_record in BND_records:
				outline = this_record + "\n"
				sys.stdout.write( outline )

		# CASE:		CHROM	POS	ALT
		# 		1	100	N]1:200]
		#		1	200	N]1:100]
		if ((count_N_left == 2) and
			(count_left_N == 0) and (count_left_NNN == 0) and (count_NNN_left == 0) and (count_right_N == 0) and (count_right_NNN == 0) and (count_N_right == 0) and (count_NNN_right == 0)):
			this_group_is_discarded = True
			for this_record in BND_records:
				outline = this_record + "\n"
				sys.stdout.write( outline )

		# CASE:		CHROM	POS	ALT
		# 		1	100	NNN]1:200]
		#		1	200	NNN]1:100]
		if ((count_NNN_left == 2) and
			(count_left_N == 0) and (count_left_NNN == 0) and (count_N_left == 0) and (count_right_N == 0) and (count_right_NNN == 0) and (count_N_right == 0) and (count_NNN_right == 0)):
			this_group_is_discarded = True
			for this_record in BND_records:
				outline = this_record + "\n"
				sys.stdout.write( outline )

		if (this_group_is_discarded == False):
			if (number_of_SV_patterns_this_BND_group_matches > 1):
				do_nothing = 1
				#sys.stderr.write('FOR EVENT ' + event_id + ' THERE WAS MORE THAN ONE OUTPUT TYPE' + "\n")
			elif (number_of_SV_patterns_this_BND_group_matches == 0):
				# This group of records has not been output at all, so output them
				for this_record in BND_records:
					outline = this_record + "\n"
					sys.stdout.write( outline )
				#if (len(BND_records) > 1):
				#	sys.stderr.write('FOR EVENT ' + event_id + ' THERE WERE NO OUTPUT TYPES' + "\n")

        return

######################################################
def main():

	output_gridss_info_fields = True
	output_gridss_format_fields = False
	current_event_id = ''
	current_event_BND_records = []

	# Read in the input VCF file from STDIN

	in_header = True
	have_outputted_new_info_headers = False
	previous_header_type = ''
	new_info_headers = []
	new_info_headers.append('##ALT=<ID=DEL,Description="Deletion">\n')
	new_info_headers.append('##ALT=<ID=INS,Description="Insertion">\n')
	new_info_headers.append('##ALT=<ID=INDEL,Description="Insertion and deletion">\n')
	new_info_headers.append('##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">\n')
	new_info_headers.append('##ALT=<ID=DUP:INS,Description="Duplication with an additional inserted sequence">\n')
	# new_info_headers.append('##INFO=<ID=TRANCHE,Number=1,Type=String,Description="Quality category of GRIDSS structural variant calls determined using QUAL,AS,RAS. Values are LOW INTERMEDIATE HIGH">\n') # Assume this header is already present from the input vcf
	# new_info_headers.append('##INFO=<ID=TRANCHE2,Number=1,Type=String,Description="Quality category of GRIDSS structural variant calls determined using FILTER,SRQ,AS,RAS. Values are LOW INTERMEDIATE HIGH">\n') # Assume this header is already present from the input vcf
	new_info_headers.append('##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">\n')
	new_info_headers.append('##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">\n')
	new_info_headers.append('##INFO=<ID=VAF,Number=1,Type=String,Description="VAF of this SV call, derived from BNDVAF values of BND calls used to call this SV">\n')
	if (output_gridss_info_fields):
		new_info_headers.append('##INFO=<ID=gridssVAF,Number=.,Type=String,Description="VAF values of the gridss-called BND calls used to call this SV">\n')
		new_info_headers.append('##INFO=<ID=gridssAS,Number=.,Type=String,Description="AS values of the gridss-called BND calls used to call this SV. Count of assemblies supporting breakpoint">\n')
		new_info_headers.append('##INFO=<ID=gridssASQ,Number=.,Type=String,Description="ASQ values of the gridss-called BND calls used to call this SV. Quality score of assemblies supporting breakpoint">\n')
		new_info_headers.append('##INFO=<ID=gridssASRP,Number=.,Type=String,Description="ASRP values of the gridss-called BND calls used to call this SV. Count of read pairs incorporated into any breakpoint assembly">\n')
		new_info_headers.append('##INFO=<ID=gridssASRR,Number=.,Type=String,Description="ASRR values of the gridss-called BND calls used to call this SV. Count of split, soft clipped or indel-containing reads incorporated into any breakpoint assemblies">\n')
		new_info_headers.append('##INFO=<ID=gridssBA,Number=.,Type=String,Description="BA values of the gridss-called BND calls used to call this SV. Count of assemblies supporting just local breakend">\n')
		new_info_headers.append('##INFO=<ID=gridssBAQ,Number=.,Type=String,Description="BAQ values of the gridss-called BND calls used to call this SV. Quality score of assemblies supporting just local breakend">\n')
		new_info_headers.append('##INFO=<ID=gridssBEID,Number=.,Type=String,Description="BEID values of the gridss-called BND calls used to call this SV. Breakend assemblies contributing support to the breakpoint">\n')
		new_info_headers.append('##INFO=<ID=gridssBQ,Number=.,Type=String,Description="BQ values of the gridss-called BND calls used to call this SV. Quality score of breakend evidence">\n')
		new_info_headers.append('##INFO=<ID=gridssBSC,Number=.,Type=String,Description="BSC values of the gridss-called BND calls used to call this SV. Count of soft clips supporting just local breakend">\n')
		new_info_headers.append('##INFO=<ID=gridssBSCQ,Number=.,Type=String,Description="BSCQ values of the gridss-called BND calls used to call this SV. Quality score of soft clips supporting just local breakend">\n')
		new_info_headers.append('##INFO=<ID=gridssBUM,Number=.,Type=String,Description="BUM values of the gridss-called BND calls used to call this SV. Count of read pairs (with one read unmapped) supporting just local breakend">\n')
		new_info_headers.append('##INFO=<ID=gridssBUMQ,Number=.,Type=String,Description="BUMQ values of the gridss-called BND calls used to call this SV. Quality score of read pairs (with one read unmapped) supporting just local breakend">\n')
		new_info_headers.append('##INFO=<ID=gridssCAS,Number=.,Type=String,Description="CAS values of the gridss-called BND calls used to call this SV. Count of complex compound breakpoint assemblies supporting breakpoint from elsewhere">\n')
		new_info_headers.append('##INFO=<ID=gridssCASQ,Number=.,Type=String,Description="CASQ values of the gridss-called BND calls used to call this SV. Quality score of complex compound breakpoint assemblies supporting breakpoint from elsewhere">\n')
		new_info_headers.append('##INFO=<ID=gridssCIEND,Number=.,Type=String,Description="CIEND values of the gridss-called BND calls used to call this SV. Confidence interval around END for imprecise variants">\n')
		new_info_headers.append('##INFO=<ID=gridssCIPOS,Number=.,Type=String,Description="CIPOS values of the gridss-called BND calls used to call this SV. Confidence interval around POS for imprecise variants">\n')
		new_info_headers.append('##INFO=<ID=gridssCIRPOS,Number=.,Type=String,Description="CIRPOS values of the gridss-called BND calls used to call this SV. Confidence interval around remote breakend POS for imprecise variants">\n')
		new_info_headers.append('##INFO=<ID=gridssCQ,Number=.,Type=String,Description="CQ values of the gridss-called BND calls used to call this SV. Breakpoint quality score before evidence reallocation">\n')
		new_info_headers.append('##INFO=<ID=gridssFILTER,Number=.,Type=String,Description="FILTER values of the gridss-called BND calls used to call this SV">\n')
		new_info_headers.append('##INFO=<ID=gridssHOMLEN,Number=.,Type=String,Description="HOMLEN values of the gridss-called BND calls used to call this SV. Length of base pair identical micro-homology at event breakpoints">\n')
		new_info_headers.append('##INFO=<ID=gridssHOMSEQ,Number=.,Type=String,Description="HOMSEQ values of the gridss-called BND calls used to call this SV. Sequence of base pair identical micro-homology at event breakpoints">\n')
		new_info_headers.append('##INFO=<ID=gridssIC,Number=.,Type=String,Description="IC values of the gridss-called BND calls used to call this SV. Count of read indels supporting breakpoint">\n')
		new_info_headers.append('##INFO=<ID=gridssIHOMPOS,Number=.,Type=String,Description="IHOMPOS values of the gridss-called BND calls used to call this SV. Position of inexact homology">\n')
		new_info_headers.append('##INFO=<ID=gridssIQ,Number=.,Type=String,Description="IQ values of the gridss-called BND calls used to call this SV. Quality score of read indels supporting breakpoint">\n')
		new_info_headers.append('##INFO=<ID=gridssQUAL,Number=.,Type=String,Description="QUAL values of the gridss-called BND calls used to call this SV.">\n')
		new_info_headers.append('##INFO=<ID=gridssRAS,Number=.,Type=String,Description="RAS values of the gridss-called BND calls used to call this SV. Count of assemblies supporting breakpoint from remote breakend">\n')
		new_info_headers.append('##INFO=<ID=gridssRASQ,Number=.,Type=String,Description="RASQ values of the gridss-called BND calls used to call this SV. Quality score of assemblies supporting breakpoint from remote breakend">\n')
		new_info_headers.append('##INFO=<ID=gridssREF,Number=.,Type=String,Description="REF values of the gridss-called BND calls used to call this SV. Count of reads mapping across this breakend">\n')
		new_info_headers.append('##INFO=<ID=gridssREFPAIR,Number=.,Type=String,Description="REFPAIR values of the gridss-called BND calls used to call this SV. Count of reference read pairs spanning this breakpoint supporting the reference allele">\n')
		new_info_headers.append('##INFO=<ID=gridssRP,Number=.,Type=String,Description="RP values of the gridss-called BND calls used to call this SV. Count of read pairs supporting breakpoint">\n')
		new_info_headers.append('##INFO=<ID=gridssRPQ,Number=.,Type=String,Description="RPQ values of the gridss-called BND calls used to call this SV. Quality score of read pairs supporting breakpoint">\n')
		new_info_headers.append('##INFO=<ID=gridssRSI,Number=.,Type=String,Description="RSI values of the gridss-called BND calls used to call this SV. Support interval offsets of partner breakend">\n')
		new_info_headers.append('##INFO=<ID=gridssSI,Number=.,Type=String,Description="SI values of the gridss-called BND calls used to call this SV. Support interval offsets from breakend position in which at least one supporting read/read pair/assembly is mapped">\n')
		new_info_headers.append('##INFO=<ID=gridssSR,Number=.,Type=String,Description="SR values of the gridss-called BND calls used to call this SV. Count of split reads supporting breakpoint">\n')
		new_info_headers.append('##INFO=<ID=gridssSRQ,Number=.,Type=String,Description="SRQ values of the gridss-called BND calls used to call this SV. Quality score of split reads supporting breakpoint">\n')
	new_info_headers.append('##INFO=<ID=INSSEQ,Number=.,Type=String,Description="Insertion sequence of structural variant, not including sequence marked as duplication">\n')
	new_info_headers.append('##FORMAT=<ID=POS,Number=1,Type=Integer,Description="Start position of the variant described in this record">\n')
	new_info_headers.append('##FORMAT=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
	new_info_headers.append('##FORMAT=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
	# new_info_headers.append('##FORMAT=<ID=TRANCHE,Number=1,Type=String,Description="Quality category of GRIDSS structural variant calls determined using QUAL,AS,RAS. Values are LOW INTERMEDIATE HIGH">\n') # Assume this header is already present from the input vcf
	# new_info_headers.append('##FORMAT=<ID=TRANCHE2,Number=1,Type=String,Description="Quality category of GRIDSS structural variant calls determined using FILTER,SRQ,AS,RAS. Values are LOW INTERMEDIATE HIGH">\n') # Assume this header is already present from the input vcf
	new_info_headers.append('##FORMAT=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">\n')
	new_info_headers.append('##FORMAT=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">\n')
	new_info_headers.append('##FORMAT=<ID=VAF,Number=1,Type=String,Description="VAF of this SV call, derived from BNDVAF values of BND calls used to call this SV">\n')
	if (output_gridss_format_fields):
		new_info_headers.append('##FORMAT=<ID=gridssVAF,Number=.,Type=String,Description="VAF values of the gridss-called BND calls used to call this SV">\n')
		new_info_headers.append('##FORMAT=<ID=gridssAS,Number=.,Type=String,Description="AS values of the gridss-called BND calls used to call this SV. Count of assemblies supporting breakpoint">\n')
		new_info_headers.append('##FORMAT=<ID=gridssASQ,Number=.,Type=String,Description="ASQ values of the gridss-called BND calls used to call this SV. Quality score of assemblies supporting breakpoint">\n')
		new_info_headers.append('##FORMAT=<ID=gridssASRP,Number=.,Type=String,Description="ASRP values of the gridss-called BND calls used to call this SV. Count of read pairs incorporated into any breakpoint assembly">\n')
		new_info_headers.append('##FORMAT=<ID=gridssASRR,Number=.,Type=String,Description="ASRR values of the gridss-called BND calls used to call this SV. Count of split, soft clipped or indel-containing reads incorporated into any breakpoint assemblies">\n')
		new_info_headers.append('##FORMAT=<ID=gridssBA,Number=.,Type=String,Description="BA values of the gridss-called BND calls used to call this SV. Count of assemblies supporting just local breakend">\n')
		new_info_headers.append('##FORMAT=<ID=gridssBAQ,Number=.,Type=String,Description="BAQ values of the gridss-called BND calls used to call this SV. Quality score of assemblies supporting just local breakend">\n')
		new_info_headers.append('##FORMAT=<ID=gridssBEID,Number=.,Type=String,Description="BEID values of the gridss-called BND calls used to call this SV. Breakend assemblies contributing support to the breakpoint">\n')
		new_info_headers.append('##FORMAT=<ID=gridssBQ,Number=.,Type=String,Description="BQ values of the gridss-called BND calls used to call this SV. Quality score of breakend evidence">\n')
		new_info_headers.append('##FORMAT=<ID=gridssBSC,Number=.,Type=String,Description="BSC values of the gridss-called BND calls used to call this SV. Count of soft clips supporting just local breakend">\n')
		new_info_headers.append('##FORMAT=<ID=gridssBSCQ,Number=.,Type=String,Description="BSCQ values of the gridss-called BND calls used to call this SV. Quality score of soft clips supporting just local breakend">\n')
		new_info_headers.append('##FORMAT=<ID=gridssBUM,Number=.,Type=String,Description="BUM values of the gridss-called BND calls used to call this SV. Count of read pairs (with one read unmapped) supporting just local breakend">\n')
		new_info_headers.append('##FORMAT=<ID=gridssBUMQ,Number=.,Type=String,Description="BUMQ values of the gridss-called BND calls used to call this SV. Quality score of read pairs (with one read unmapped) supporting just local breakend">\n')
		new_info_headers.append('##FORMAT=<ID=gridssCAS,Number=.,Type=String,Description="CAS values of the gridss-called BND calls used to call this SV. Count of complex compound breakpoint assemblies supporting breakpoint from elsewhere">\n')
		new_info_headers.append('##FORMAT=<ID=gridssCASQ,Number=.,Type=String,Description="CASQ values of the gridss-called BND calls used to call this SV. Quality score of complex compound breakpoint assemblies supporting breakpoint from elsewhere">\n')
		new_info_headers.append('##FORMAT=<ID=gridssCIEND,Number=.,Type=String,Description="CIEND values of the gridss-called BND calls used to call this SV. Confidence interval around POS for imprecise variants">\n')
		new_info_headers.append('##FORMAT=<ID=gridssCIPOS,Number=.,Type=String,Description="CIPOS values of the gridss-called BND calls used to call this SV. Confidence interval around POS for imprecise variants">\n')
		new_info_headers.append('##FORMAT=<ID=gridssCIRPOS,Number=.,Type=String,Description="CIRPOS values of the gridss-called BND calls used to call this SV. Confidence interval around remote breakend POS for imprecise variants">\n')
		new_info_headers.append('##FORMAT=<ID=gridssCQ,Number=.,Type=String,Description="CQ values of the gridss-called BND calls used to call this SV. Breakpoint quality score before evidence reallocation">\n')
		new_info_headers.append('##FORMAT=<ID=gridssFILTER,Number=.,Type=String,Description="FILTER values of the gridss-called BND calls used to call this SV">\n')
		new_info_headers.append('##FORMAT=<ID=gridssHOMLEN,Number=.,Type=String,Description="HOMLEN values of the gridss-called BND calls used to call this SV. Length of base pair identical micro-homology at event breakpoints">\n')
		new_info_headers.append('##FORMAT=<ID=gridssHOMSEQ,Number=.,Type=String,Description="HOMSEQ values of the gridss-called BND calls used to call this SV. Sequence of base pair identical micro-homology at event breakpoints">\n')
		new_info_headers.append('##FORMAT=<ID=gridssIC,Number=.,Type=String,Description="IC values of the gridss-called BND calls used to call this SV. Count of read indels supporting breakpoint">\n')
		new_info_headers.append('##FORMAT=<ID=gridssIHOMPOS,Number=.,Type=String,Description="IHOMPOS values of the gridss-called BND calls used to call this SV. Position of inexact homology">\n')
		new_info_headers.append('##FORMAT=<ID=gridssIQ,Number=.,Type=String,Description="IQ values of the gridss-called BND calls used to call this SV. Quality score of read indels supporting breakpoint">\n')
		new_info_headers.append('##FORMAT=<ID=gridssQUAL,Number=.,Type=String,Description="QUAL values of the gridss-called BND calls used to call this SV.">\n')
		new_info_headers.append('##FORMAT=<ID=gridssRAS,Number=.,Type=String,Description="RAS values of the gridss-called BND calls used to call this SV. Count of assemblies supporting breakpoint from remote breakend">\n')
		new_info_headers.append('##FORMAT=<ID=gridssRASQ,Number=.,Type=String,Description="RASQ values of the gridss-called BND calls used to call this SV. Quality score of assemblies supporting breakpoint from remote breakend">\n')
		new_info_headers.append('##FORMAT=<ID=gridssREF,Number=.,Type=String,Description="REF values of the gridss-called BND calls used to call this SV. Count of reads mapping across this breakend">\n')
		new_info_headers.append('##FORMAT=<ID=gridssREFPAIR,Number=.,Type=String,Description="REFPAIR values of the gridss-called BND calls used to call this SV. Count of reference read pairs spanning this breakpoint supporting the reference allele">\n')
		new_info_headers.append('##FORMAT=<ID=gridssRP,Number=.,Type=String,Description="RP values of the gridss-called BND calls used to call this SV. Count of read pairs supporting breakpoint">\n')
		new_info_headers.append('##FORMAT=<ID=gridssRPQ,Number=.,Type=String,Description="RPQ values of the gridss-called BND calls used to call this SV. Quality score of read pairs supporting breakpoint">\n')
		new_info_headers.append('##FORMAT=<ID=gridssRSI,Number=.,Type=String,Description="RSI values of the gridss-called BND calls used to call this SV. Support interval offsets of partner breakend">\n')
		new_info_headers.append('##FORMAT=<ID=gridssSI,Number=.,Type=String,Description="SI values of the gridss-called BND calls used to call this SV. Support interval offsets from breakend position in which at least one supporting read/read pair/assembly is mapped">\n')
		new_info_headers.append('##FORMAT=<ID=gridssSR,Number=.,Type=String,Description="SR values of the gridss-called BND calls used to call this SV. Count of split reads supporting breakpoint">\n')
		new_info_headers.append('##FORMAT=<ID=gridssSRQ,Number=.,Type=String,Description="SRQ values of the gridss-called BND calls used to call this SV. Quality score of split reads supporting breakpoint">\n')
	new_info_headers.append('##FORMAT=<ID=INSSEQ,Number=.,Type=String,Description="Insertion sequence of structural variant, not including sequence marked as duplication">\n')

	for inline in sys.stdin:
		inline = inline.strip()

		if (in_header == True):
			if (len(inline) >= 1):
				first_char = inline[0:1]
				if (first_char != '#'):
					in_header = False

		if (in_header == True):

			# At the end of existing INFO headers, output new INFO headers

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

			outline = inline + "\n"
			sys.stdout.write( outline )

		else: # We are processing VCF data records. We are no longer in the header part of the file.

			infields = inline.split("\t")
			chrom = str(infields[0])
			pos = infields[1]
			if (is_integer(pos)):
				pos = int(pos)
			snpid = infields[2]
			ref = infields[3]
			alt = infields[4]
			qual = infields[5]
			snpfilter = infields[6]
			info = infields[7]
			snpformat = infields[8]
			if (info.find('SVTYPE=BND') > -1):
				info_fields = info.split(';')
				this_event_id = ''
				for this_info_field in info_fields:
					bits = this_info_field.split('=')
					if (bits[0] == 'EVENT'):
						this_event_id = str(bits[1])
				if (this_event_id != ''):
					if (this_event_id == current_event_id):
						# Continue accumulating BND records for the current event.
						current_event_BND_records.append(inline)
					else:
						# We have seen all the BND records for the current event,
						# so write out the one VCF record for them.
						if (len(current_event_BND_records) > 0):
							write_one_summary_VCF_record( current_event_BND_records, current_event_id, output_gridss_info_fields, output_gridss_format_fields )
						# Then start accumulating BND records for the new event.
						current_event_BND_records = []
						current_event_BND_records.append(inline)
						current_event_id = this_event_id

	# write out the last VCF record if there is a final group of BND records for an event not yet written out
	if (len(current_event_BND_records) > 0):
		write_one_summary_VCF_record( current_event_BND_records, current_event_id, output_gridss_info_fields, output_gridss_format_fields )


if __name__=='__main__':
    main()


