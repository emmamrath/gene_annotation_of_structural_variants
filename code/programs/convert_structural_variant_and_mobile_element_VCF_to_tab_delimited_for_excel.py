#!/usr/bin/python
# python convert_structural_variant_and_mobile_element_VCF_to_tab_delimited_for_excel.py -i input_VCF_file -o output_tsv_file


# This program reads in a VCF file of structural variants from GRIDSS and mobile elements from Mobster
# that has been annotated in the VCF.INFO field with genes and other information.
# This program outputs this information as a tab-separated file, with one output line per input VCF data record.
#
# Multiple alleles are treated as one variant with a long ALT field, eg. <INS:ME:L1>,<INS:ME:ALU>
# This program assumes that an INFO field contains only one entry per VCF data record for a given INFO ID, even if multiple alleles are present
# in which case only the INFO field for the first allele will be output.
# This program assumes that subfields of a given INFO ID are separated in the VCF data record by by comma ,
# Thus, an INFO field having multiple value separated by bar | will end up with all the value and the bar separators in the output field.
# Thus, an INFO field having values for multiple alleles, separated by commas, will have only the first allele's value in the output field.
# This program outputs one INFO column where the INFO field does not have a specific number specified, eg. for Number=R, Number=A, Number=.
# When the INFO field does have a specific number of subfields specified, then this program outputs one INFO column for each INFO subfield.
# The output column heading of each INFO subfield will be the INFO name followed by a number,
# except in the case of MEINFO in which case the output column headings are MEINFO_NAME, MEINFO_START, MEINFO_END, MEINFO_POLARITY
# This program assumes that samples in the VCF record that don't have this particular variant will have ./. in their sample field.
# This program outputs a count of how many samples have each variant. The header field contains a count of the total number of samples as NUM_SAMPLES=1011


# Example input file that has mobile elements that have extra annotation for them:
#
# ##fileformat=VCFv4.2
# ##FILTER=<ID=PASS,Description="All filters passed">
# ##fileDate=20170331
# ##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">
# ##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">
# ##ALT=<ID=INS:ME:HERV,Description="Insertion of HERV element">
# ##ALT=<ID=INS:ME:SVA,Description="Insertion of SVA element">
# ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
# ##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">
# ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
# ##INFO=<ID=TARGETSITEDUPL,Number=1,Type=String,Description="Did a small part of this recognition site get duplicated or deleted during the insertion event. Values:duplication,deletion,noTSD,unknown">
# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# ##FORMAT=<ID=MEI5MB,Number=1,Type=String,Description="Mobster_border_5_one_base_pair_after_MEI_3_prime of Mobster MEI call.">
# ##FORMAT=<ID=MEI3MB,Number=1,Type=String,Description="Mobster_border_3_one_base_pair_before_MEI_5_prime of Mobster MEI call.">
# ##FORMAT=<ID=MEIDU,Number=1,Type=String,Description="Mobster target_site_duplication.">
# ##FORMAT=<ID=MEIPLY,Number=1,Type=String,Description="poly-A or poly-T sequence seen in MEI or near insertion point. Values are A,T (for whether a poly-A or poly-T sequence was seen respectively), indicating gene direction">
# ##FORMAT=<ID=MEI5RB,Number=1,Type=String,Description="refined_border_5_one_base_pair_after_MEI_3_prime of refinement of mobile element insertion (MEI) call from running convert_imprecise_Mobster_predictions_into_precise.py on Mobster MEI call.">
# ##FORMAT=<ID=MEI3RB,Number=1,Type=String,Description="refined_border_3_one_base_pair_before_MEI_5_prime of refinement of mobile element insertion (MEI) call from running convert_imprecise_Mobster_predictions_into_precise.py on Mobster MEI call.">
# ##FORMAT=<ID=MEI5SP,Number=1,Type=String,Description="refined_border_5_MEI_fragment_start_pos of refinement of mobile element insertion (MEI) call from running convert_imprecise_Mobster_predictions_into_precise.py on Mobster MEI call.">
# ##FORMAT=<ID=MEI3SP,Number=1,Type=String,Description="refined_border_3_MEI_fragment_start_pos of refinement of mobile element insertion (MEI) call from running convert_imprecise_Mobster_predictions_into_precise.py on Mobster MEI call.">
# ##FORMAT=<ID=MEI5PR,Number=1,Type=String,Description="refined_border_5_precision of refinement of mobile element insertion (MEI) call from running convert_imprecise_Mobster_predictions_into_precise.py on Mobster MEI call. Values are PRECISE,IMPRECISE">
# ##FORMAT=<ID=MEI3PR,Number=1,Type=String,Description="refined_border_3_precision of refinement of mobile element insertion (MEI) call from running convert_imprecise_Mobster_predictions_into_precise.py on Mobster MEI call. Values are PRECISE,IMPRECISE">
# ##FORMAT=<ID=MEI5SQ,Number=1,Type=String,Description="refined_border_5_MEI_fragment_sequence of refinement of mobile element insertion (MEI) call from running convert_imprecise_Mobster_predictions_into_precise.py on Mobster MEI call.">
# ##FORMAT=<ID=MEI3SQ,Number=1,Type=String,Description="refined_border_3_MEI_fragment_sequence of refinement of mobile element insertion (MEI) call from running convert_imprecise_Mobster_predictions_into_precise.py on Mobster MEI call.">
# ##contig=<ID=1>
# ##contig=<ID=10>
# ##contig=<ID=11>
# ##contig=<ID=12>
# ##contig=<ID=13>
# ##contig=<ID=14>
# ##contig=<ID=GL000215.1>
# ##contig=<ID=GL000207.1>
# ##contig=<ID=GL000210.1>
# ##bcftools_mergeVersion=1.3.1+htslib-1.3.1
# ##bcftools_mergeCommand=merge -m none -o MGRBp1_all_samples_mobile_elements_sameSampleFields.vcf MGRBp1_AAAAL_sample12_mobile_elements_withPolyACounts_sameSampleFields_sorted.vcf.gz MGRBp1_AAAEJ_sample104_mobile_elements_withPolyACounts_sameSampleFields_sorted.vcf.gz MGRBp1_
# ##INFO=<ID=GenomicRegionGene,Number=4,Type=String,Description="Mobile element info of the form TranscriptId|RegionType|GeneName|PfamId. Multiple PfamId separated by & (And sign). Multiple GenomicRegionGene for same ALT allele separated by exclamation mark !.">
# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  AAAAL   AAAEJ   AAAFY   AAAGU   AAAID   AAAIE   AAAIF   AAAIG   AAAIH   AAAII   AAAIJ   AAAIK   AAAIL   AAAIM   AAAIN   AAAIO   AAAIP   AAAIQ   AAAIR   AAAIS   AAAIT   AAAIU   AAAIV   AAAIW   AAAIX   AAA
# 1       10281   .       A       <INS:ME:L1>     16      .       IMPRECISE;SVTYPE=INS;MEINFO=L1,-84,27,.;TARGETSITEDUPL=duplication      GT:MEI5MB:MEI3MB:MEIDU:MEIPLY:MEI5RB:MEI3RB:MEI5SP:MEI3SP:MEI5PR:MEI3PR:MEI5SQ:MEI3SQ   ./.:.:.:.:.:.:.:.:.:.:.:.:.     ./.:.:.:.:.:.:.:.:.
# 1       10450   .       T       <INS:ME:L1>     15      .       IMPRECISE;SVTYPE=INS;MEINFO=L1,-176,710,.;TARGETSITEDUPL=duplication    GT:MEI5MB:MEI3MB:MEIDU:MEIPLY:MEI5RB:MEI3RB:MEI5SP:MEI3SP:MEI5PR:MEI3PR:MEI5SQ:MEI3SQ   ./.:.:.:.:.:.:.:.:.:.:.:.:.     ./.:.:.:.:.:.:.:.:.
# 1       1243933 .       G       <INS:ME:ALU>,<INS:ME:L1>        18      .       IMPRECISE;SVTYPE=INS;MEINFO=ALU,-16,17,.;TARGETSITEDUPL=duplication;GenomicRegionGene=ENST00000472541|retained_intron|ACAP3|!ENST00000379031|protein_coding|PUSL1|PF01416,ENST00000472541|retained_
# 1       1244069 .       G       <INS:ME:L1>,<INS:ME:ALU>        13      .       IMPRECISE;SVTYPE=INS;MEINFO=L1,-156,157,.;TARGETSITEDUPL=unknown        GT:MEI5MB:MEI3MB:MEIDU:MEIPLY:MEI5RB:MEI3RB:MEI5SP:MEI3SP:MEI5PR:MEI3PR:MEI5SQ:MEI3SQ   ./.:.:.:.:.:.:.:.:.:.:.:.:.     ./.
# 1       818025  .       A       <INS:ME:L1>     26      .       IMPRECISE;SVTYPE=INS;MEINFO=L1,-20,20,.;TARGETSITEDUPL=unknown;GenomicRegionGene=ENST00000594233|protein_coding|AL645608.2|     GT:MEI5MB:MEI3MB:MEIDU:MEIPLY:MEI5RB:MEI3RB:MEI5SP:MEI3SP:MEI5PR:MEI3PR:MEI5SQ:MEI3
# 1       16637573        .       A       <INS:ME:ALU>    25      .       IMPRECISE;SVTYPE=INS;MEINFO=ALU,-135,133,.;TARGETSITEDUPL=deletion;GenomicRegionGene=ENST00000375592|protein_coding|FBXO42|PF01344&PF07646&PF00646!ENST00000478089|processed_transcript|FBXO42|PF00646&PF07

# Here are more examples of the types of VCF.INFO fields that may be received on input by this program:
#
# ##INFO=<ID=AS,Number=1,Type=Integer,Description="Count of assemblies supporting breakpoint">
# ##INFO=<ID=ASQ,Number=1,Type=Float,Description="Quality score of assemblies supporting breakpoint">
# ##INFO=<ID=ASRP,Number=1,Type=Integer,Description="Count of read pairs incorporated into any breakpoint assembly">
# ##INFO=<ID=ASSR,Number=1,Type=Integer,Description="Count of split, soft clipped or indel-containing reads incorporated into any breakpoint assemblies">
# ##INFO=<ID=BA,Number=1,Type=Integer,Description="Count of assemblies supporting just local breakend">
# ##INFO=<ID=BAQ,Number=1,Type=Float,Description="Quality score of assemblies supporting just local breakend">
# ##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
# ##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
# ##INFO=<ID=IC,Number=1,Type=Integer,Description="Count of read indels supporting breakpoint">
# ##INFO=<ID=IHOMPOS,Number=2,Type=Integer,Description="Position of inexact homology">
# ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
# ##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
# ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
# ##INFO=<ID=TRANCHE,Number=.,Type=String,Description="Quality category of GRIDSS structural variant calls. Values are LOW INTERMEDIATE HIGH">
# ##INFO=<ID=INSSEQ,Number=.,Type=String,Description="Insertion sequence of structural variant, not including sequence marked as duplication">
# ##INFO=<ID=HITGENES,Number=A,Type=String,Description="Intersecting genes">
# ##INFO=<ID=HITREGULATORY,Number=A,Type=String,Description="Intersecting regulatory regions">
# ##INFO=<ID=HITMOBEL,Number=A,Type=String,Description="Intersecting regulatory regions">
# ##INFO=<ID=HITEXONS,Number=A,Type=String,Description="Intersecting gene exons">
# ##INFO=<ID=HITINTRONS,Number=A,Type=String,Description="Intersecting gene introns">

# Here are some examples of input data:
# left side of data file:
#
# 1       566040  MEI_1_566040    C       <INS:ME:ALU>,<INS:ME:L1>        16      .       AC=1;AF=0.500;AN=2;IMPRECISE;MEINFO=ALU,-807,261,.;SVTYPE=INS;TARGETSITEDUPL=unknown;set=variant352;HITGENES=JA429830--|JA429831--|JB137814--;HITREGULATORY=proximal_39|
# 1       566212  MEI_1_566212    C       <INS:ME:L1>     17      .       AC=1;AF=0.500;AN=2;IMPRECISE;MEINFO=L1,-683,108,.;SVTYPE=INS;TARGETSITEDUPL=unknown;set=variant153;HITGENES=JA429830--|JA429831--|JB137814--;HITREGULATORY=ctcf_14|proximal_40;HITEXONS=
# 1       566264  MEI_1_566264    A       <INS:ME:L1>     20      .       AC=1;AF=0.500;AN=2;IMPRECISE;MEINFO=L1,-528,54,.;SVTYPE=INS;TARGETSITEDUPL=;set=variant110;HITGENES=JA429830--|JA429831--|JB137814--;HITREGULATORY=ctcf_14|proximal_40;HITEXONS=JA429830
# 1       566437  MEI_1_566437    T       <INS:ME:L1>,<INS:ME:ALU>        16      .       AC=1;AF=0.500;AN=2;IMPRECISE;MEINFO=L1,-854,285,.;SVTYPE=INS;TARGETSITEDUPL=unknown;set=variant339;HITGENES=JA429830--|JA429831--|JB137814--;HITREGULATORY=ctcf_14|proxi
# 1       566505  MEI_1_566505    C       <INS:ME:L1>,<INS:ME:ALU>        16      .       AC=1;AF=0.500;AN=2;IMPRECISE;MEINFO=L1,-688,748,.;SVTYPE=INS;TARGETSITEDUPL=unknown;set=variant361;HITGENES=JA429830--|JA429831--|JB137814--;HITREGULATORY=ctcf_14|proxi
# 1       566603  MEI_1_566603    C       <INS:ME:L1>     16      .       AC=1;AF=0.500;AN=2;IMPRECISE;MEINFO=L1,-252,793,.;SVTYPE=INS;TARGETSITEDUPL=duplication;set=variant988;HITREGULATORY=ctcf_14|proximal_40        GT:MEI3MB:MEI3PR:MEI3RB:MEI3SP:MEI3SQ:ME
#
# right side of data file:
# ctcf_14|proximal_40;HITEXONS=JA429830--|JA429831--|JB137814--   GT:MEI3MB:MEI3PR:MEI3RB:MEI3SP:MEI3SQ:MEI5MB:MEI5PR:MEI5RB:MEI5SP:MEI5SQ:MEIDU:MEIPLY   ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     
# JA429830--|JA429831--|JB137814--        GT:MEI3MB:MEI3PR:MEI3RB:MEI3SP:MEI3SQ:MEI5MB:MEI5PR:MEI5RB:MEI5SP:MEI5SQ:MEIDU:MEIPLY   ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     
# --|JA429831--|JB137814--        GT:MEI3MB:MEI3PR:MEI3RB:MEI3SP:MEI3SQ:MEI5MB:MEI5PR:MEI5RB:MEI5SP:MEI5SQ:MEIDU:MEIPLY   ./.     ./.     ./.     ./.     ./.     0/1:566282:PRECISE:566321:566160:ACATCCACAGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTACTAGACCAATGGGACT
# mal_40;HITEXONS=JA429830--|JA429831--|JB137814--        GT:MEI3MB:MEI3PR:MEI3RB:MEI3SP:MEI3SQ:MEI5MB:MEI5PR:MEI5RB:MEI5SP:MEI5SQ:MEIDU:MEIPLY   ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     0/1:566595:PRECISE:566767:566641:CATGGAATGTTCTTC
# mal_40;HITEXONS=JA429830--|JA429831--|JB137814--        GT:MEI3MB:MEI3PR:MEI3RB:MEI3SP:MEI3SQ:MEI5MB:MEI5PR:MEI5RB:MEI5SP:MEI5SQ:MEIDU:MEIPLY   ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     0/1:566697:PRECISE:566296:566135
# I5MB:MEI5PR:MEI5RB:MEI5SP:MEI5SQ:MEIDU:MEIPLY   ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     ./.     

# More input data examples from GRIDSS:
#
# ##INFO=<ID=IHOMPOS,Number=2,Type=Integer,Description="Position of inexact homology">
#
# ;BUMQ=0.00;CAS=0;CASQ=0.00;CQ=220.05;EVENT=gridss0_1994;IC=0;IHOMPOS=9,-1;IQ=0.00;PARID=gridss0_1994h;RAS=0;RASQ=0.00;REF=354;REFPAIR=223;RP=0;RPQ=0.00;SC
# 0;BUMQ=0.00;CAS=0;CASQ=0.00;CQ=220.05;EVENT=gridss0_1994;IC=0;IHOMPOS=44,-44;IQ=0.00;PARID=gridss0_1994o;RAS=2;RASQ=291.63;REF=538;REFPAIR=152;RP=0;RPQ=0.


__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2017, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import argparse
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
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in an uncompressed VCF file containing structural variants (from GRIDSS) and mobile element insertions (from Mobster) that have INFO field annotations, and output this record in tab-separated format.')
	parser.add_argument('-i', action="store", dest="in_VCF", required=True, help='Input VCF file containing structural variants, mobile elements, their annotations, and their samples')
	parser.add_argument('-o', action="store", dest="out_TSV", required=False, help='Output tab-delimited Excel file - one line per VCF record')
	args = parser.parse_args()

	# read each record of input VCF file and process it

	in_VCF = open(args.in_VCF, 'r')
	out_TSV = open(args.out_TSV, 'w')

	in_header = True
	out_TSV_header = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER"
	sample_FORMAT_part_of_out_TSV_header = ""

	# ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
	# ##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
	# ##INFO=<ID=HITGENES,Number=A,Type=String,Description="Intersecting genes">
	# ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
	# ##INFO=<ID=IHOMPOS,Number=2,Type=Integer,Description="Position of inexact homology">
	# ##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">

	INFO_field_order = []
	INFO_field_count_of_subfields = {}
	INFO_field_names_of_subfields = {}

	for inline in in_VCF:

		inline = inline.strip()

		if (in_header == True):

			# process each header line read in

			if (inline[0:1] == '#'):
				in_header = True
				if (len(inline) > 11):

					# Process the INFO header lines.

					if (inline[0:11] == '##INFO=<ID='):

						# Add this INFO field and it's subfields to the Excel header lines

						# ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
						# ##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">
						# ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
						# ##INFO=<ID=TARGETSITEDUPL,Number=1,Type=String,Description="Did a small part of this recognition site get duplicated or deleted during the insertion event. Values:duplication,deletion,noTSD,unknown">
						# ##INFO=<ID=GenomicRegionGene,Number=4,Type=String,Description="Mobile element info of the form TranscriptId|RegionType|GeneName|PfamId. Multiple PfamId separated by & (And sign). Multiple GenomicRegionGene for same ALT allele separated by exclamation mark !.">

						info_regex = re.compile('##INFO=<ID=(.*)>')
						this_INFO_definition = info_regex.search(inline).group(1)

						bits_by_comma = this_INFO_definition.split(',')
						INFO_header = bits_by_comma[0]
						INFO_Number_field = bits_by_comma[1]
						INFO_Number_field = INFO_Number_field.replace('Number=','')
						INFO_field_subfields = []

						if (is_integer(INFO_Number_field)):

							# This INFO field has a fixed number of subfields. There will be multiple output columns, one column per subfield.

							num_output_subfields = int(INFO_Number_field)
							if (num_output_subfields == 0):
								num_output_subfields = 1
							if (num_output_subfields == 1):
								out_TSV_header = out_TSV_header + "\t" + INFO_header
								INFO_field_order.append( INFO_header )
								INFO_field_count_of_subfields[INFO_header] = 1
							else:
								if ((INFO_header == 'MEINFO') and (num_output_subfields == 4)):
									info_subfield_keys = []
									MEINFO_subfields = [ 'MEINFO_NAME', 'MEINFO_START', 'MEINFO_END', 'MEINFO_POLARITY' ]
									for INFO_subfield_name in MEINFO_subfields:
										out_TSV_header = out_TSV_header + "\t" + INFO_subfield_name
										INFO_field_order.append( INFO_subfield_name )
										info_subfield_keys.append( INFO_subfield_name )
									INFO_field_count_of_subfields[INFO_header] = 4
									INFO_field_names_of_subfields[INFO_header] = info_subfield_keys
								else:
									info_subfield_keys = []
									for idx0 in range( 0, num_output_subfields ):
										idx1 = idx0 + 1
										INFO_subfield_name = INFO_header + '_' + str(idx1)
										out_TSV_header = out_TSV_header + "\t" + INFO_subfield_name
										INFO_field_order.append( INFO_subfield_name )
										info_subfield_keys.append( INFO_subfield_name )
									INFO_field_count_of_subfields[INFO_header] = num_output_subfields
									INFO_field_names_of_subfields[INFO_header] = info_subfield_keys
						else:

							# Only one subfield column will be output for this INFO field, regardless of how many each variants has - can vary for each allele

							out_TSV_header = out_TSV_header + "\t" + INFO_header
							INFO_field_order.append( INFO_header )
							INFO_field_count_of_subfields[INFO_header] = 1

					elif (inline[0:6] == '#CHROM'):

						infields = inline.split("\t")
						num_samples = len(infields) - 9
						if (num_samples < 0):
							num_samples = 0
						sample_field_name = 'NUM_SAMPLES(TOTAL=' + str(num_samples) + ')'
						out_TSV_header = out_TSV_header + "\t" + sample_field_name + "\t" + 'NUM_HETEROZYGOUS' + "\t" + 'NUM_HOMOZYGOUS' + "\t" + 'NUM_UNCERTAIN_ZYGOSITY'

			else:

				# when finished reading in VCF header records, write out the Excel header record

				in_header = False

				out_TSV.write(out_TSV_header + "\n")

		if (in_header == False):

			# Process each data line read in.
			# For each VCF record, process the multiple ALTs separately, one at a time.
			# This means process the line once for the first ALT, then reprocess the line a second time for the second ALT, etc.
			# For the output file that has one line per VCF record and not one line per ALT (ie. multiple ALTs at same POS will produce 1 output line),
			# then output only the first of each INFO annotation (for some INFO fields, there is only 1 INFO annotation, for others there is one per ALT).
			# For INFO fields that can have more than one set of entries, separated by !, output only the first one.

			# 1       10281   .       A       <INS:ME:L1>     16      .       IMPRECISE;SVTYPE=INS;MEINFO=L1|-84|,27|,.;TARGETSITEDUPL=duplication      GT:MEI5MB:MEI3MB:MEIDU:MEIPLY:MEI5RB:MEI3RB:MEI5SP:MEI3SP:MEI5PR:MEI3PR:MEI5SQ:MEI3SQ   ./.:.:.:.:.:.:.:.:.:.:.:.:.     ./.:.:.:.:.:.:.:.:.
			# 1       10450   .       T       <INS:ME:L1>     15      .       IMPRECISE;SVTYPE=INS;MEINFO=L1|-176|710|.;TARGETSITEDUPL=duplication    GT:MEI5MB:MEI3MB:MEIDU:MEIPLY:MEI5RB:MEI3RB:MEI5SP:MEI3SP:MEI5PR:MEI3PR:MEI5SQ:MEI3SQ   ./.:.:.:.:.:.:.:.:.:.:.:.:.     ./.:.:.:.:.:.:.:.:.
			# 1       1243933 .       G       <INS:ME:ALU>,<INS:ME:L1>        18      .       IMPRECISE;SVTYPE=INS;MEINFO=ALU|-16|17|.;TARGETSITEDUPL=duplication;GenomicRegionGene=ENST00000472541|retained_intron|ACAP3|!ENST00000379031|protein_coding|PUSL1|PF01416,ENST00000472541|retained_
			# 1       1244069 .       G       <INS:ME:L1>,<INS:ME:ALU>        13      .       IMPRECISE;SVTYPE=INS;MEINFO=L1|-156|157|.;TARGETSITEDUPL=unknown        GT:MEI5MB:MEI3MB:MEIDU:MEIPLY:MEI5RB:MEI3RB:MEI5SP:MEI3SP:MEI5PR:MEI3PR:MEI5SQ:MEI3SQ   ./.:.:.:.:.:.:.:.:.:.:.:.:.     ./.
			# 1       818025  .       A       <INS:ME:L1>     26      .       IMPRECISE;SVTYPE=INS;MEINFO=L1|-20|20|.;TARGETSITEDUPL=unknown;GenomicRegionGene=ENST00000594233|protein_coding|AL645608.2|     GT:MEI5MB:MEI3MB:MEIDU:MEIPLY:MEI5RB:MEI3RB:MEI5SP:MEI3SP:MEI5PR:MEI3PR:MEI5SQ:MEI3
			# 1       16637573        .       A       <INS:ME:ALU>    25      .       IMPRECISE;SVTYPE=INS;MEINFO=ALU|-135|133|.;TARGETSITEDUPL=deletion;GenomicRegionGene=ENST00000375592|protein_coding|FBXO42|PF01344&PF07646&PF00646!ENST00000478089|processed_transcript|FBXO42|PF00646&PF07

			INFO_field_output = [''] * len(INFO_field_order)

			infields = inline.split("\t")
			CHROM = str(infields[0])
			POS = str(infields[1])
			SNV_ID = str(infields[2])
			REF = str(infields[3])
			all_ALTs = str(infields[4])
			QUAL = str(infields[5])
			FILTER = str(infields[6])

			outline = CHROM + "\t" + POS + "\t" + SNV_ID + "\t" + REF + "\t" + all_ALTs + "\t" + QUAL + "\t" + FILTER

			INFO = infields[7]
			INFO_fields_and_subfields = {}
			INFO_fields = INFO.split(';')
			for INFO_field in INFO_fields:
				key_and_value = INFO_field.split('=')
				info_key = key_and_value[0]
				if (len(key_and_value) <= 1):
					info_value = info_key # eg. IMPRECISE key has no value. Value assigned will be the key 'IMPRECISE'
					INFO_fields_and_subfields[info_key] = info_value
				else:
					if (info_key in INFO_field_count_of_subfields): # Ignore this variant's info field if it wasn't in the header
						if (INFO_field_count_of_subfields[info_key] == 1): # eg. SVTYPE has one and only one value. eg. for HITGENES we take only the first value which is for first allele
							info_value = key_and_value[1]
							INFO_fields_and_subfields[info_key] = info_value
						elif (INFO_field_count_of_subfields[info_key] == 0): # for eg. IMPRECISE key that has no value, but should be processed above already, not here
							info_value = info_key # eg. IMPRECISE key has no value. Value assigned will be the key 'IMPRECISE'
							INFO_fields_and_subfields[info_key] = info_value
						else: # eg. MEINFO has 4 subfields. eg. IHOMPOS has 2 subfields
							info_values = key_and_value[1]
							info_subfield_keys = INFO_field_names_of_subfields[info_key]
							info_subfield_values = info_values.split(',')
							for i in range( 0, INFO_field_count_of_subfields[info_key] ):
								info_subfield_key = info_subfield_keys[i]
								info_subfield_value = info_subfield_values[i]
								INFO_fields_and_subfields[info_subfield_key] = info_subfield_value

			# INFO_field_count_of_subfields[INFO_header]

			for i in range( 0, len(INFO_field_order) ): # For every info field in header, write out this variant's value for it
				this_info_key = INFO_field_order[i]
				this_info_value = ''
				if this_info_key in INFO_fields_and_subfields:
					this_info_value = str(INFO_fields_and_subfields[this_info_key])
				outline = outline + "\t" + this_info_value

			count_samples_in_this_variant = 0
			count_heterozygous_in_this_variant = 0
			count_homozygous_in_this_variant = 0
			count_unknownzygous_in_this_variant = 0
			if (len(infields) >= 10):
				for i in range( 9, len(infields) ):
					this_sample_format_field_string = infields[i]
					this_sample_format_fields = this_sample_format_field_string.split(':')
					this_sample = this_sample_format_fields[0].strip()
					if ((this_sample != './.') and (this_sample != '.')):
						count_samples_in_this_variant = count_samples_in_this_variant + 1
						this_sample_split = this_sample.split('/')
						this_sample_simplified = ''
						if (this_sample_split[0] == '.'):
							this_sample_simplified = '.'
						elif (this_sample_split[0] == '0'):
							this_sample_simplified = '0'
						else:
							this_sample_simplified = '1'
						this_sample_simplified = this_sample_simplified + '/'
						# For diploid, there should be 2 alleles.
						# Let's accept only one because we don't want to crash when there are problems in the data
						if (len(this_sample_split) > 1):
							if (this_sample_split[1] == '.'):
								this_sample_simplified = this_sample_simplified + '.'
							elif (this_sample_split[1] == '0'):
								this_sample_simplified = this_sample_simplified + '0'
							else:
								this_sample_simplified = this_sample_simplified + '1'
						if ((this_sample_simplified == '0/1') or (this_sample_simplified == '1/0')):
							count_heterozygous_in_this_variant = count_heterozygous_in_this_variant + 1
						elif (this_sample_simplified == '1/1'):
							count_homozygous_in_this_variant = count_homozygous_in_this_variant + 1
						elif ((this_sample_simplified == '1/.') or (this_sample_simplified == './1')):
							count_unknownzygous_in_this_variant = count_unknownzygous_in_this_variant + 1
						elif (this_sample_simplified == '0/0'):
							do_not_count_this = True
						else:
							print 'ERROR: Unexpected genotype', this_sample, 'encountered for variant', CHROM, POS, SNV_ID, REF, all_ALTs

			outline = outline + "\t" + str(count_samples_in_this_variant) + "\t" + str(count_heterozygous_in_this_variant) + "\t" + str(count_homozygous_in_this_variant) + "\t" + str(count_unknownzygous_in_this_variant)

			out_TSV.write(outline + "\n")

	out_TSV.close()


if __name__=='__main__':
    main()

