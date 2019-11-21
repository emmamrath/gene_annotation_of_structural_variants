#!/usr/bin/python
# python filter_out_large_variants.py sample_id
# cat test_SV.vcf | python filter_out_large_variants.py -f 1000001 > test_SV_small.vcf

# This program removes VCF records whose size is larger than the given size.

# Filter out structural variants whose size is smaller than or equal to the provided filter.
# Use INFO.END to determine the size. If INFO.END is not present, then use abs(INFO.SVLEN).
# If INFO.SVLEN is not present, then use INFO.MEINFO.
# If none of those are present and the record is a BND record, then determine the size from the ALT field.

# If there are multiple BND alleles, then if any of them are larger than the given size, the VCF record is discarded.

# ##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">

# There may be multiple ALT fields for non-BND records, but they will only have one INFO field 
# which is assumed to apply to all the ALT fields.

# There may be multiple ALT fields for BND records and none of the above INFO fields (ie. INFO.END, INFO.SVLEN, INFO.MEINFO).
# If one but not all of multiple ALT fields should be filtered, then keep them all
# because this program is not going to adjust genotypes for multiple alleles.
# If there are multiple ALT BND fields, then filter out the VCF record only if all of them are smaller in size than the filter.


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
def main():

	parser = argparse.ArgumentParser(description='Filter out structural variants larger than given size')
	parser.add_argument('-f', action="store", dest="filter_length", required=True, help='filter to discard SVs whose length (INFO.END,INFO.SVLEN.INFO.MEINFO.END is larger than or equal to this filter length')
	args = parser.parse_args()
	filter_length = int(args.filter_length)

	# Read in the input VCF file from STDIN

	in_header = True
	
	for inline in sys.stdin:

		if (in_header == True):
			if (len(inline) >= 1):
				first_char = inline[0:1]
				if (first_char != '#'):
					in_header = False

		if (in_header):
			sys.stdout.write( inline )

		else: # in_header == False # We are processing VCF data records. We are no longer in the header part of the file.

			inline = inline.strip()
			infields = inline.split("\t")
			this_chrom = str(infields[0])
			this_pos = int(infields[1])
			this_alts = str(infields[4])
			this_info = str(infields[7])

			got_length = False
			this_length = 0

			is_an_insertion = False
			this_alts_fields = this_alts.split(',')
			for this_alt in this_alts_fields:
				if ( this_alt.find('<INS') > -1 ):
					is_an_insertion = True
				elif ( this_alt.find('<INDEL') > -1 ):
					is_an_insertion = True

			if (is_an_insertion):
				info_fields = this_info.split(';')
				for this_field in info_fields:
					bits = this_field.split('=')
					this_key = bits[0]
					if (this_key == 'SVLEN'):
						this_length = abs(int(bits[1]))
						got_length = True

			else:
				info_fields = this_info.split(';')
				for this_field in info_fields:
					bits = this_field.split('=')
					this_key = bits[0]
					if (this_key == 'END'):
						if (got_length == False): # SVLEN takes precedence over END
							this_length = int(bits[1]) - this_pos + 1
							got_length = True
					elif (this_key == 'SVLEN'):
						this_length = abs(int(bits[1]))
						got_length = True
					# MEINFO=L1,-20,20,. gives position of mobile element insertion, doesn't give length

			# T]1:146701554]
			# [1:182274316[TTTTTTTTTTTTT
			# G[1:87502881[
			# AAAA]1:36734642]
			if ((got_length == False) and (is_an_insertion == False)):

				this_alts_fields = this_alts.split(',')
				got_at_least_one_alt_length = False
				smallest_alt_length = 0
				for this_alt in this_alts_fields:

					got_this_alt_length = False
					this_alt_length = 0
					this_alt = this_alt.replace( ']', '[' )
					idx = this_alt.find('[')
					if (idx > -1):
						alt_bits = this_alt.split('[')
						for this_bit in alt_bits:
							if (len(this_bit) > 0):
								idx = this_bit.find(':')
								if (idx > -1):
									bits2 = this_bit.split(':')
									alt_chrom = str(bits2[0])
									alt_pos = int(bits2[1])
									if (this_chrom == alt_chrom):
										this_alt_length = this_alt_length + abs(alt_pos - this_pos + 1)
										got_this_alt_length = True
								else:
									this_alt_length = this_alt_length + len(this_bit) - 1
					if (got_this_alt_length):
						got_at_least_one_alt_length = True
						smallest_alt_length = min(this_alt_length, smallest_alt_length)
				if (got_at_least_one_alt_length):
					got_length = True
					this_length = smallest_alt_length

			if (got_length):
				if (this_length <= filter_length):
					outline = inline + "\n"
					sys.stdout.write( outline )
			else:
				outline = inline + "\n"
				sys.stdout.write( outline )


if __name__=='__main__':
    main()


