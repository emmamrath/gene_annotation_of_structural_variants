#!/usr/bin/python
# python split_UCSC_refseq_genes_to_genes_exons_introns.py -i input_genes_file -out_genes output_genes_file -out_exons output_exons_file -out_introns out_introns_file -out_cds output_CDS_file -out_cdsexons output_CDSexons_file -out_utr5 output_utr5_file -out_utr3 output_utr3_file

# This program reads in UCSC genes, and outputs separate files containing genes, exons, and introns.
# Input and output files are in bed format.

# Example input genes:
# hg19_ncbiRefSeq_bin     refseq_transcript       hg19_ncbiRefSeq_chrom   hg19_ncbiRefSeq_strand  hg19_ncbiRefSeq_txStart hg19_ncbiRefSeq_txEnd   hg19_ncbiRefSeq_cdsStart        hg19_ncbiRefSeq_cdsEnd
# 0       NM_032291       chr1    +       66999638        67216822        67000041        67208778        25      66999638,67091529,67098752,67101626,67105459,67108492,67109226,67126195,67133212,671366
# 1       NM_001080397    chr1    +       8378144 8404227 8378168 8404073 9       8378144,8384365,8385357,8385877,8390268,8395496,8397875,8399552,8403806,        8378246,8384786,8385450,8386102,8390996
# 1       NM_001145277    chr1    +       16767166        16786585        16767256        16785491        7       16767166,16770126,16774364,16774554,16775587,16778332,16785336, 16767348,16770227,16774
# 1       NM_013943       chr1    +       25071759        25170815        25072044        25167428        6       25071759,25124232,25140584,25153500,25166350,25167263,  25072116,25124342,25140710,2515
# 1       NM_052998       chr1    +       33546713        33586132        33547850        33585783        12      33546713,33546988,33547201,33547778,33549554,33557650,33558882,33560148,33562307,335636
# 1       NM_032785       chr1    -       48998526        50489626        48999844        50489468        14      48998526,49000561,49005313,49052675,49056504,49100164,49119008,49128823,49332862,495112
# 2       NM_003243       chr1    -       92145899        92351836        92149295        92327088        17      92145899,92161228,92163645,92174219,92177799,92181792,92182124,92184868,92185449,921875
# 2       NM_001918       chr1    -       100652477       100715409       100661810       100715376       11      100652477,100671785,100672000,100676249,100680372,100681538,100684181,100696288,1007009
# 3       NM_021222       chr1    +       150980866       151008189       150981108       151006710       8       150980866,150990287,150990942,150997086,150997990,150999708,151001261,151006281,

# The header of the refseq genes input file:
# hg19_ncbiRefSeq_bin	refseq_transcript	hg19_ncbiRefSeq_chrom	hg19_ncbiRefSeq_strand	hg19_ncbiRefSeq_txStart	hg19_ncbiRefSeq_txEnd	hg19_ncbiRefSeq_cdsStart	hg19_ncbiRefSeq_cdsEndnum_refseq_exons	hg19_ncbiRefSeq_exonStarts	hg19_ncbiRefSeq_exonEnds	hg19_ncbiRefSeq_score	refseq_gene	hg19_ncbiRefSeq_cdsStartStat	hg19_ncbiRefSeq_cdsEndStat	hg19_ncbiRefSeq_exonFrames	hg19_ncbiRefSeqHgmd_name2	hg19_ncbiRefSeqLink_id	hg19_ncbiRefSeqLink_status	hg19_ncbiRefSeqLink_name	hg19_ncbiRefSeqLink_product	hg19_ncbiRefSeqLink_mrnaAcc	hg19_ncbiRefSeqLink_protAcc	hg19_ncbiRefSeqLink_locusLinkId	hg19_ncbiRefSeqLink_omimId	hg19_ncbiRefSeqLink_hgnc	hg19_ncbiRefSeqLink_genbank	hg19_ncbiRefSeqLink_pseudo	hg19_ncbiRefSeqLink_gbkey	hg19_ncbiRefSeqLink_source	hg19_ncbiRefSeqLink_gene_biotype	hg19_ncbiRefSeqLink_gene_synonym	hg19_ncbiRefSeqLink_ncrna_class	hg19_ncbiRefSeqLink_note	hg19_ncbiRefSeqLink_description	hg19_ncbiRefSeqLink_externalId

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2019, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

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
			in_chrom = str(infields[2])
			in_strand = str(infields[3])
			in_txStart = str(infields[4])
			in_txEnd = str(infields[5])
			in_cdsStart = str(infields[6])
			in_cdsEnd = str(infields[7])
			in_exonCount = str(infields[8])
			in_exonStarts = str(infields[9])
			in_exonEnds = str(infields[10])
			in_geneSymbol = str(infields[12])

			outline_gene = in_chrom + "\t" + in_txStart + "\t" + in_txEnd + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
			out_genes.write( outline_gene )
			outline_cds = in_chrom + "\t" + in_cdsStart + "\t" + in_cdsEnd + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
			out_cds.write( outline_cds )
			if (int(in_txStart) < int(in_cdsStart)):
				outline_utr5 = in_chrom + "\t" + in_txStart + "\t" + in_cdsStart + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
				out_utr5.write( outline_utr5 )
			if (int(in_txEnd) > int(in_cdsEnd)):
				outline_utr3 = in_chrom + "\t" + in_cdsEnd + "\t" + in_txEnd + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
				out_utr3.write( outline_utr3 )

			in_exon_starts = in_exonStarts.split(",")
			in_exon_ends = in_exonEnds.split(",")
			in_cdsStart = int(in_cdsStart)
			in_cdsEnd = int(in_cdsEnd)
			for j in range( 0, len(in_exon_starts) ):
				this_exon_start = in_exon_starts[j]
				this_exon_end = in_exon_ends[j]
				if (this_exon_start != ''):
					this_exon_start = int(this_exon_start)
					this_exon_end = int(this_exon_end)
					outline_exons = in_chrom + "\t" + str(this_exon_start) + "\t" + str(this_exon_end) + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
					out_exons.write( outline_exons )
					this_exon_falls_inside_cds = False
					this_cdsexon_start = this_exon_start
					this_cdsexon_end = this_exon_end
					if ((in_cdsStart <= this_exon_start) and (in_cdsEnd >= this_exon_end)):
						this_exon_falls_inside_cds = True
					else:
						if ((this_exon_start <= in_cdsStart) and (in_cdsStart <= this_exon_end)):
							this_exons_falls_inside_cds = True
							this_cdsexon_start = in_cdsStart
						if ((this_exon_end >= in_cdsEnd) and (in_cdsEnd >= this_exon_start)):
							this_exon_falls_inside_cds = True
							this_cdsexon_end = in_cdsEnd
					if (this_exon_falls_inside_cds):
						outline_cdsexons = in_chrom + "\t" + str(this_cdsexon_start) + "\t" + str(this_cdsexon_end) + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
						out_cdsexons.write( outline_cdsexons )

			first_exon_starts = str(in_exon_starts[0])
			last_exon_ends = str(in_exon_ends[ len(in_exon_ends)-1 ])
			if (last_exon_ends == ''):
				last_exon_ends = str(in_exon_ends[ len(in_exon_ends)-2 ])
			if ((in_txStart != first_exon_starts) or (in_txEnd != last_exon_ends)):
				print 'Unexpected mismatch of transcript and exon starts and ends for', in_txStart, first_exon_starts, in_txEnd, last_exon_ends, 'inline', inline

			for j in range( 1, len(in_exon_starts) ):
				if (in_exon_starts[j] != ''):
					this_intron_start = str( int(in_exon_ends[j-1]) + 1 )
					this_intron_end = str( int(in_exon_starts[j]) - 1 )
					if (int(this_intron_start) < int(this_intron_end)):
						outline_introns = in_chrom + "\t" + this_intron_start + "\t" + this_intron_end + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
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


