#!/usr/bin/python
# python add_sample_key_to_VCF_to_merge_multiple_samples.py sample_id
# cat AAAPI_bicseq_CNV.vcf | python add_sample_key_to_VCF_to_merge_multiple_samples.py AAAPI > AAAPI_bicseq_CNV_with_sample_key.txt

# Put the sample id in front of each VCF record so that we can subsequently
# sort the VCF records of multiple samples and then merge into one large multi-sample VCF.
# Which field to take for QUAL?

# cat AAAPA_bicseq_CNV.vcf | python add_sample_key_to_VCF_to_merge_multiple_samples.py AAAPA > AAAPA_bicseq_CNV_with_sample_key.txt
# cat AAAPB_bicseq_CNV.vcf | python add_sample_key_to_VCF_to_merge_multiple_samples.py AAAPB > AAAPB_bicseq_CNV_with_sample_key.txt
# cat AAAPC_bicseq_CNV.vcf | python add_sample_key_to_VCF_to_merge_multiple_samples.py AAAPC > AAAPC_bicseq_CNV_with_sample_key.txt


import sys
import os
import commands

######################################################
def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

######################################################
def main():

	# Read in the input VCF file from STDIN

	sample_id = str(sys.argv[1])
	in_header = True
	
	for inline in sys.stdin:

		if (in_header == True):
			if (len(inline) >= 1):
				first_char = inline[0:1]
				if (first_char != '#'):
					in_header = False

		if (in_header == False): # We are processing VCF data records. We are no longer in the header part of the file.

			# Is this variant of a type that we are interested in merging?

			outline = sample_id + "\t" + inline
			sys.stdout.write( outline )


if __name__=='__main__':
    main()


