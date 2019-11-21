#!/usr/bin/python
# python split_UCSC_color_genes_to_genes_exons_introns.py -l input_genes_file -out_genes output_genes_file -out_exons output_exons_file -out_introns out_introns_file -out_cds output_CDS_file -out_cdsexons output_CDSexons_file -out_utr5 output_utr5_file -out_utr3 output_utr3_file

# This program reads in UCSC genes, and outputs separate files containing genes, exons, and introns.
# Input and output files are in bed format.

# Example input genes:
# #hg19.knownGene.name    hg19.knownGene.chrom    hg19.knownGene.strand   hg19.knownGene.txStart  hg19.knownGene.txEnd    hg19.knownGene.cdsStart hg19.knownGene.cdsEnd   hg19.knownGene.exonCount        hg19.knownGene.exonStarts       hg19.knownGe
# uc001aaa.3      chr1    +       11873   14409   11873   11873   3       11873,12612,13220,      12227,12721,14409,      uc001aaa.3      130     130     210     uc001aaa.3      DDX11L1
# uc010nxr.1      chr1    +       11873   14409   11873   11873   3       11873,12645,13220,      12227,12697,14409,      uc010nxr.1      130     130     210     uc010nxr.1      DDX11L1
# uc010nxq.1      chr1    +       11873   14409   12189   13639   3       11873,12594,13402,      12227,12721,14409,      uc010nxq.1      130     130     210     uc010nxq.1      DDX11L1
# uc009vis.3      chr1    -       14361   16765   14361   14361   4       14361,14969,15795,16606,        14829,15038,15942,16765,        uc009vis.3      130     130     210     uc009vis.3      WASH7P
# uc009vjc.1      chr1    -       16857   17751   16857   16857   2       16857,17232,    17055,17751,    uc009vjc.1      130     130     210     uc009vjc.1      WASH7P
# uc009vjd.2      chr1    -       15795   18061   15795   15795   5       15795,16606,16857,17232,17605,  15947,16765,17055,17368,18061,  uc009vjd.2      130     130     210     uc009vjd.2      WASH7P


__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2017, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

import sys
import os
import subprocess
import argparse
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

	parser = argparse.ArgumentParser(description='Read in UCSC table text file of genes, output bed files for genes, exons and introns.')
	parser.add_argument('-i', action="store", dest="in_UCSC", required=True, help='Input UCSC table text file list of genes')
	parser.add_argument('-out_genes', action="store", dest="out_genes", required=True, help='Output bed format file of genes')
	parser.add_argument('-out_exons', action="store", dest="out_exons", required=True, help='Output bed format file of exons')
	parser.add_argument('-out_introns', action="store", dest="out_introns", required=True, help='Output bed format file of introns')
	parser.add_argument('-out_cds', action="store", dest="out_cds", required=True, help='Output bed format file of CDS')
	parser.add_argument('-out_cdsexons', action="store", dest="out_cdsexons", required=True, help='Output bed format file of CDS exons')
	parser.add_argument('-out_utr5', action="store", dest="out_utr5", required=True, help='Output bed format file of UTR-5')
	parser.add_argument('-out_utr3', action="store", dest="out_utr3", required=True, help='Output bed format file of UTR-3')
	args = parser.parse_args()

	in_UCSC = open(args.in_UCSC, 'r')
	inlines = in_UCSC.readlines()
	in_UCSC.close()

	out_genes = open(args.out_genes, 'w')
	out_exons = open(args.out_exons, 'w')
	out_introns = open(args.out_introns, 'w')
	out_cds = open(args.out_cds, 'w')
	out_cdsexons = open(args.out_cdsexons, 'w')
	out_utr5 = open(args.out_utr5, 'w')
	out_utr3 = open(args.out_utr3, 'w')

	for i in range( 1, len(inlines) ): # skip the first line, it is the header line
		inline = inlines[i]
		stripped_inline = inline.strip()
		if (stripped_inline != ''):
			infields = inline.split("\t")
			infields[ len(infields)-1 ] = infields[ len(infields)-1 ].strip() # remove the carriage return on last field
			input_field_name_1 = str(infields[0])
			input_field_chrom = str(infields[1])
			input_field_strand = str(infields[2])
			input_field_txStart = str(infields[3])
			input_field_txEnd = str(infields[4])
			input_field_cdsStart = str(infields[5])
			input_field_cdsEnd = str(infields[6])
			input_field_exonCount = str(infields[7])
			input_field_exonStarts = str(infields[8])
			input_field_exonEnds = str(infields[9])
			input_field_kgID_1 = str(infields[10])
			input_field_color_r = int(infields[11])
			input_field_color_g = int(infields[12])
			input_field_color_b = int(infields[13])
			input_field_kgID_2 = str(infields[14])
			input_field_geneSymbol = str(infields[15])

			outline_gene = input_field_chrom + "\t" + input_field_txStart + "\t" + input_field_txEnd + "\t" + input_field_geneSymbol + "\t" + input_field_strand + "\n"
			out_genes.write( outline_gene )
			outline_cds = input_field_chrom + "\t" + input_field_cdsStart + "\t" + input_field_cdsEnd + "\t" + input_field_geneSymbol + "\t" + input_field_strand + "\n"
			out_cds.write( outline_cds )
			if (int(input_field_txStart) < int(input_field_cdsStart)):
				outline_utr5 = input_field_chrom + "\t" + input_field_txStart + "\t" + input_field_cdsStart + "\t" + input_field_geneSymbol + "\t" + input_field_strand + "\n"
				out_utr5.write( outline_utr5 )
			if (int(input_field_txEnd) > int(input_field_cdsEnd)):
				outline_utr3 = input_field_chrom + "\t" + input_field_cdsEnd + "\t" + input_field_txEnd + "\t" + input_field_geneSymbol + "\t" + input_field_strand + "\n"
				out_utr3.write( outline_utr3 )

			in_exon_starts = input_field_exonStarts.split(",")
			in_exon_ends = input_field_exonEnds.split(",")
			input_field_cdsStart = int(input_field_cdsStart)
			input_field_cdsEnd = int(input_field_cdsEnd)
			for j in range( 0, len(in_exon_starts) ):
				this_exon_start = in_exon_starts[j]
				this_exon_end = in_exon_ends[j]
				if (this_exon_start != ''):
					this_exon_start = int(this_exon_start)
					this_exon_end = int(this_exon_end)
					outline_exons = input_field_chrom + "\t" + str(this_exon_start) + "\t" + str(this_exon_end) + "\t" + input_field_geneSymbol + "\t" + input_field_strand + "\n"
					out_exons.write( outline_exons )
					this_exon_falls_inside_cds = False
					this_cdsexon_start = this_exon_start
					this_cdsexon_end = this_exon_end
					if ((input_field_cdsStart <= this_exon_start) and (input_field_cdsEnd >= this_exon_end)):
						this_exon_falls_inside_cds = True
					else:
						if ((this_exon_start <= input_field_cdsStart) and (input_field_cdsStart <= this_exon_end)):
							this_exons_falls_inside_cds = True
							this_cdsexon_start = input_field_cdsStart
						if ((this_exon_end >= input_field_cdsEnd) and (input_field_cdsEnd >= this_exon_start)):
							this_exon_falls_inside_cds = True
							this_cdsexon_end = input_field_cdsEnd
					if (this_exon_falls_inside_cds):
						outline_cdsexons = input_field_chrom + "\t" + str(this_cdsexon_start) + "\t" + str(this_cdsexon_end) + "\t" + input_field_geneSymbol + "\t" + input_field_strand + "\n"
						out_cdsexons.write( outline_cdsexons )

			first_exon_starts = str(in_exon_starts[0])
			last_exon_ends = str(in_exon_ends[ len(in_exon_ends)-1 ])
			if (last_exon_ends == ''):
				last_exon_ends = str(in_exon_ends[ len(in_exon_ends)-2 ])
			if ((input_field_txStart != first_exon_starts) or (input_field_txEnd != last_exon_ends)):
				print 'Unexpected mismatch of transcript and exon starts and ends for', input_field_txStart, first_exon_starts, input_field_txEnd, last_exon_ends, 'inline', inline

			for j in range( 1, len(in_exon_starts) ):
				if (in_exon_starts[j] != ''):
					this_intron_start = str( int(in_exon_ends[j-1]) + 1 )
					this_intron_end = str( int(in_exon_starts[j]) - 1 )
					if (int(this_intron_start) < int(this_intron_end)):
						outline_introns = input_field_chrom + "\t" + this_intron_start + "\t" + this_intron_end + "\t" + input_field_geneSymbol + "\t" + input_field_strand + "\n"
						out_introns.write( outline_introns )

	out_genes.close()
	out_exons.close()
	out_introns.close()
	out_cds.close()
	out_cdsexons.close()
	out_utr5.close()
	out_utr3.close()

if __name__=='__main__':
    main()


