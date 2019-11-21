#!/usr/bin/python
# python convert_VCF_INFO_one_sample_to_tab_delimited_and_intact_fields.py -i input_VCF_file -o output_tsv_file

# This program reads in a VCF file and outputs all the existing tab-delimited columns for each VCF record. 
# In addition, for every INFO field defined in the header, this program produces an output column, for every VCF record.
# It assumes there is zero or only one sample. 
# This program outputs one column for the sample GT (assumed to be first sample format field) and one for the entire sample format.
# The output columns are CHROM, POS, ID, REF, ALT, QUAL, FILTER, FORMAT, GT(contains sample GT), <sample_id>(contains entire sample format), INFO columns.
# Thus, no data is lost (apart from the VCF headers) and output could be used to reconstruct the input data.

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2019, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import argparse
import re
import gc

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

	parser = argparse.ArgumentParser(description='Read in an uncompressed VCF file. Output tab-delimited file having one column per input INFO field.')
	parser.add_argument('-i', action="store", dest="in_VCF", required=True, help='Input VCF file containing structural variants, mobile elements, their annotations, and zero or one sample')
	parser.add_argument('-o', action="store", dest="out_TSV", required=False, help='Output tab-delimited file - one line per VCF record')
	args = parser.parse_args()

	# read each record of input VCF file and process it

	in_VCF = open(args.in_VCF, 'r')
	out_TSV = open(args.out_TSV, 'w')

	in_header = True
	has_sample = False
	out_TSV_header = ""
	out_TSV_INFO_header = ""

	INFO_field_order = []

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

						if (INFO_header == 'MEINFO'):
							MEINFO_subfields = [ 'MEINFO_NAME', 'MEINFO_START', 'MEINFO_END', 'MEINFO_POLARITY' ]
							for INFO_subfield_name in MEINFO_subfields:
								if (INFO_subfield_name not in INFO_field_order):
									INFO_field_order.append( INFO_subfield_name )
									if (out_TSV_INFO_header == ""):
										out_TSV_INFO_header = INFO_subfield_name
									else:
										out_TSV_INFO_header = out_TSV_INFO_header + "\t" + INFO_subfield_name
						else:
							if (INFO_header not in INFO_field_order):
								#gc.disable()
								INFO_field_order.append( INFO_header ) # In some python versions, this will take too long when got lots of info fields
								#gc.enable()
								if (out_TSV_INFO_header == ""):
									out_TSV_INFO_header = INFO_header
								else:
									out_TSV_INFO_header = out_TSV_INFO_header + "\t" + INFO_header

					elif (inline[0:6] == '#CHROM'):

						out_TSV_header = inline[1:]
						infields = inline.split("\t")
						if (len(infields) > 8):
							has_sample = True
							out_TSV_header = out_TSV_header + "\t" + "GT"
						out_TSV_header = out_TSV_header + "\t" + out_TSV_INFO_header

			else:

				# when finished reading in VCF header records, write out the tab-delimited file header record

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

			outline = inline

			if (has_sample):
				this_format = str(infields[8])
				this_sample = str(infields[9])
				bits = this_sample.split(':')
				this_gt = str(bits[0])
				outline = outline + "\t" + this_gt

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

					if (info_key == 'MEINFO'): # MEINFO is to be split into 4 output columns
						info_values = str(key_and_value[1])
						bits = info_values.split(',')
						INFO_fields_and_subfields['MEINFO_NAME'] = str(bits[0])
						INFO_fields_and_subfields['MEINFO_START'] = str(bits[1])
						INFO_fields_and_subfields['MEINFO_END'] = str(bits[2])
						INFO_fields_and_subfields['MEINFO_POLARITY'] = str(bits[3])

					elif (info_key in INFO_field_order): # All info fields except MEINFO are output as one column
						info_value = info_key # eg. IMPRECISE key has no value. Value assigned will be the key 'IMPRECISE'
						if (len(key_and_value) > 1):
							info_value = str(key_and_value[1])
						INFO_fields_and_subfields[info_key] = info_value

					# else: Ignore this info field if it wasn't in the header

			for i in range( 0, len(INFO_field_order) ): # For every info field in header, write out this variant's value for it
				this_info_key = INFO_field_order[i]
				this_info_value = ''
				if (this_info_key in INFO_fields_and_subfields):
					this_info_value = str(INFO_fields_and_subfields[this_info_key])
				outline = outline + "\t" + this_info_value

			out_TSV.write(outline + "\n")

	out_TSV.close()


if __name__=='__main__':
    main()

