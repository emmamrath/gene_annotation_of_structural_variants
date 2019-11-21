# Rscript choose_canonical_exons_from_UCSC_data_for_SV_consequences.R

# head -n 1 UCSC_tables_GRCh37_UCSC_genes_20190716.txt
# #hg19.knownGene.name	hg19.knownGene.chrom	hg19.knownGene.strand	hg19.knownGene.txStart	hg19.knownGene.txEnd	hg19.knownGene.cdsStart	hg19.knownGene.cdsEnd	hg19.knownGene.exonCount	hg19.knownGene.exonStarts	hg19.knownGene.exonEnds	hg19.knownGene.proteinID	hg19.knownGene.alignID	hg19.kgXref.kgID	hg19.kgXref.mRNA	hg19.kgXref.spID	hg19.kgXref.spDisplayID	hg19.kgXref.geneSymbol	hg19.kgXref.refseq	hg19.kgXref.protAcc	hg19.kgXref.description	hg19.kgXref.rfamAcc	hg19.kgXref.tRnaName	hg19.knownCanonical.chrom	hg19.knownCanonical.chromStart	hg19.knownCanonical.chromEnd	hg19.knownCanonical.clusterId	hg19.knownCanonical.transcript	hg19.knownCanonical.protein

# head -n 1 UCSC_tables_GRCh37_UCSC_genes_20190716.txt
# #hg19.knownGene.name	hg19.knownGene.chrom	hg19.knownGene.strand	hg19.knownGene.txStart	hg19.knownGene.txEnd	hg19.knownGene.cdsStart	hg19.knownGene.cdsEnd	hg19.knownGene.exonCount	hg19.knownGene.exonStarts	hg19.knownGene.exonEnds	hg19.knownGene.proteinID	hg19.knownGene.alignID	hg19.kgXref.kgID	hg19.kgXref.mRNA	hg19.kgXref.spID	hg19.kgXref.spDisplayID	hg19.kgXref.geneSymbol	hg19.kgXref.refseq	hg19.kgXref.protAcc	hg19.kgXref.description	hg19.kgXref.rfamAcc	hg19.kgXref.tRnaName	hg19.knownCanonical.chrom	hg19.knownCanonical.chromStart	hg19.knownCanonical.chromEnd	hg19.knownCanonical.clusterId	hg19.knownCanonical.transcript	hg19.knownCanonical.protein

options(width=240)
options(stringsAsFactors = FALSE)
library(data.table)
library(sqldf)
library(reshape2)

# args = commandArgs(trailingOnly=TRUE) # for production
# args=c( 'UCSC_GRCh37_RefSeq_genes_20190708.txt', 'UCSC_tables_UCSC_genes_20190716.txt' )
# refseq_file = args[1]
# ucsc_file = args[2]

refseq_file = 'UCSC_tables_GRCh37_RefSeq_genes_20190716.txt'
ucsc_file = 'UCSC_tables_GRCh37_UCSC_genes_20190716.txt'

output_refseq_transcripts_file = 'UCSC_tables_GRCh37_RefSeq_genes_20190716_canonical_transcripts.txt'

refseq = read.table( refseq_file, sep="\t", header=TRUE, fill=TRUE, quote='"', comment.char = "" )
ucsc = read.table( ucsc_file, sep="\t", header=TRUE, fill=TRUE, quote='"', comment.char = "" )

colnames(refseq) = gsub( "\\.","_", colnames(refseq) )
colnames(refseq) = gsub( "^X_","", colnames(refseq) )
colnames(ucsc) = gsub( "\\.","_", colnames(ucsc) )
colnames(ucsc) = gsub( "^X_","", colnames(ucsc) )

names(refseq)[names(refseq)=='hg19_ncbiRefSeq_name2'] = 'refseq_gene'
names(refseq)[names(refseq)=='hg19_ncbiRefSeq_name'] = 'refseq_transcript'
names(refseq)[names(refseq)=='hg19_ncbiRefSeq_exonCount'] = 'num_refseq_exons'
names(ucsc)[names(ucsc)=='hg19_kgXref_geneSymbol'] = 'ucsc_gene'
names(ucsc)[names(ucsc)=='hg19_knownCanonical_transcript'] = 'ucsc_canonical'
names(ucsc)[names(ucsc)=='hg19_kgXref_refseq'] = 'ucsc_key_for_refseq_transcript'

refseq$refseq_transcript = ifelse( is.na(refseq$refseq_transcript), "", refseq$refseq_transcript )
refseq$refseq_gene = ifelse( is.na(refseq$refseq_gene), "", refseq$refseq_gene )
ucsc$ucsc_gene = ifelse( is.na(ucsc$ucsc_gene), "", ucsc$ucsc_gene )
ucsc$ucsc_canonical = ifelse( is.na(ucsc$ucsc_canonical), "", ucsc$ucsc_canonical )
ucsc$ucsc_canonical = ifelse( (ucsc$ucsc_canonical=="n/a"), "", ucsc$ucsc_canonical )
ucsc$ucsc_key_for_refseq_transcript = ifelse( is.na(ucsc$ucsc_key_for_refseq_transcript), "", ucsc$ucsc_key_for_refseq_transcript )
ucsc$ucsc_key_for_refseq_transcript = ifelse( (ucsc$ucsc_key_for_refseq_transcript=="n/a"), "", ucsc$ucsc_key_for_refseq_transcript )

refseq$refseq_transcript = gsub( "\\..*", "", refseq$refseq_transcript )
refseq$num_refseq_exons = as.numeric(refseq$num_refseq_exons)

list_of_refseq_genes = unique(refseq$refseq_gene)
length(list_of_refseq_genes) # 27115 refseq gene names

refseq_transcripts_counts = sqldf("select refseq_gene, count(*) from refseq group by refseq_gene")
colnames(refseq_transcripts_counts) = c( "refseq_gene", "num_transcripts" )

list_of_refseq_genes_having_one_refseq_transcript = unique(refseq_transcripts_counts[ (refseq_transcripts_counts$num_transcripts==1), c("refseq_gene") ])
list_of_refseq_genes_having_more_than_one_refseq_transcript = unique(refseq_transcripts_counts[ (refseq_transcripts_counts$num_transcripts>1), c("refseq_gene") ])
list_of_refseq_transcripts_having_more_than_one_refseq_transcript = unique( refseq[ (refseq$refseq_gene %in% list_of_refseq_genes_having_more_than_one_refseq_transcript), c("refseq_transcript") ] )

length(list_of_refseq_genes_having_one_refseq_transcript) # 15424 refseq genes having only 1 transcript in refseq so we will use that transcript
refseq_transcripts_for_genes_having_one_refseq_transcript = unique(refseq[ (refseq$refseq_gene %in% list_of_refseq_genes_having_one_refseq_transcript), ])
list_of_refseq_transcripts_for_genes_having_one_refseq_transcript = unique(refseq_transcripts_for_genes_having_one_refseq_transcript$refseq_transcript)
nrow(refseq_transcripts_for_genes_having_one_refseq_transcript) # 15424 refseq transcripts of refseq genes having only 1 transcript in refseq so we will use that transcript
length(list_of_refseq_transcripts_for_genes_having_one_refseq_transcript) ##### <===== 15424 refseq transcripts of refseq genes having only 1 transcript in refseq so we will use that transcript

length(list_of_refseq_genes_having_more_than_one_refseq_transcript) # 11691 refseq genes having more than 1 transcript in refseq so we will go to UCSC table to get the canonical transcript
length(list_of_refseq_transcripts_having_more_than_one_refseq_transcript) # 42941 refseq transcripts belonging to the 11691 refseq genes that have more than 1 transcript in refseq, so we will go to UCSC table to get the canonical transcript

list_of_ucsc_canonical_transcripts_for_refseq_genes_having_more_than_one_refseq_transcript = unique( ucsc[ ((ucsc$ucsc_key_for_refseq_transcript %in% list_of_refseq_transcripts_having_more_than_one_refseq_transcript) & (ucsc$ucsc_canonical != "")), c("ucsc_key_for_refseq_transcript") ] )
list_refseq_genes_having_more_than_one_reseq_transcript_and_have_one_or_more_ucsc_canonical_transcripts = unique( refseq[ ((refseq$refseq_transcript %in% list_of_ucsc_canonical_transcripts_for_refseq_genes_having_more_than_one_refseq_transcript) & !(refseq$refseq_gene %in% list_of_refseq_genes_having_one_refseq_transcript)), c("refseq_gene") ] )

length(list_of_ucsc_canonical_transcripts_for_refseq_genes_having_more_than_one_refseq_transcript) ##### <===== 10896 refseq transcripts that are canonical in USCS table, and we will use these refseq transcripts. There are more transcripts than genes, thus some genes must have more than one canonical transcript.
length(list_refseq_genes_having_more_than_one_reseq_transcript_and_have_one_or_more_ucsc_canonical_transcripts) # 10810 refseq genes for which we will use the transcripts identified as canonical in UCSC table

diff1 = list_of_refseq_genes[ !(list_of_refseq_genes %in% list_of_refseq_genes_having_one_refseq_transcript) ]
length(diff1) # 11691 refseq genes have more than one refseq transcript
diff2 = diff1[ !(diff1 %in% list_refseq_genes_having_more_than_one_reseq_transcript_and_have_one_or_more_ucsc_canonical_transcripts) ]
length(diff2) # 881 refseq genes have more than one refseq transcript and no UCSC canonical transcript when we compare refseq_transcript to UCSC_refseq_transcript_key
list_refseq_genes_that_dont_have_canonical_ucsc_transcript = diff2
length(list_refseq_genes_that_dont_have_canonical_ucsc_transcript) # 881 refseq genes have more than one refseq transcript and no UCSC canonical transcript when we compare refseq_transcript to UCSC_refseq_transcript_key
list_ucsc_transcripts_of_refseq_genes_having_ucsc_canonical_transcript_when_match_by_refseq_gene_name = unique( ucsc[ ((ucsc$ucsc_gene %in% list_refseq_genes_that_dont_have_canonical_ucsc_transcript) & (ucsc$ucsc_canonical != "")), c("ucsc_key_for_refseq_transcript") ] )
length(list_ucsc_transcripts_of_refseq_genes_having_ucsc_canonical_transcript_when_match_by_refseq_gene_name) # 88 of the 881 refseq genes do have a UCSC canonical transcript when we compare by gene_name not refseq_transcript. Let's see if the canonical transcripts exist in refseq.
list_refseq_transcripts_of_refseq_genes_having_ucsc_canonical_transcript_when_match_by_refseq_gene_name = unique( refseq[ ((refseq$refseq_transcript %in% list_ucsc_transcripts_of_refseq_genes_having_ucsc_canonical_transcript_when_match_by_refseq_gene_name) & !(refseq$refseq_gene %in% list_of_refseq_genes_having_one_refseq_transcript) & !(refseq$refseq_gene %in% list_of_ucsc_canonical_transcripts_for_refseq_genes_having_more_than_one_refseq_transcript)), c("refseq_transcript") ] )
length(list_refseq_transcripts_of_refseq_genes_having_ucsc_canonical_transcript_when_match_by_refseq_gene_name) # 13 of those 88 canonical transcripts actually exist in refseq, so we can use only 37 of them.
list_refseq_genes_having_ucsc_canonical_transcript_when_match_by_refseq_gene_name_and_transcript_exists_in_refseq = unique( refseq[ ((refseq$refseq_transcript %in% list_refseq_transcripts_of_refseq_genes_having_ucsc_canonical_transcript_when_match_by_refseq_gene_name) & !(refseq$refseq_gene %in% list_of_refseq_genes_having_one_refseq_transcript) & !(refseq$refseq_gene %in% list_refseq_genes_having_more_than_one_reseq_transcript_and_have_one_or_more_ucsc_canonical_transcripts)), c("refseq_gene") ] )
length(list_refseq_genes_having_ucsc_canonical_transcript_when_match_by_refseq_gene_name_and_transcript_exists_in_refseq) # 0 which means that those 37 transcripts belong to genes that we already know what transcript to take. So we still have 881 refseq genes for which we need to choose one of the multiple refseq transcripts, because they don't have UCSC canonical transcripts (or they do, but the UCSC transcript matches to another gene in refseq)

refseq_having_no_ucsc_canonical = unique( refseq[ (refseq$refseq_gene %in% list_refseq_genes_that_dont_have_canonical_ucsc_transcript), ] )
refseq_highest_num_exons = unique( sqldf("select refseq_gene, max(num_refseq_exons) from refseq_having_no_ucsc_canonical group by refseq_gene") )
colnames(refseq_highest_num_exons) = c( "refseq_gene", "num_refseq_exons" )
nrow(refseq_highest_num_exons) # 881
refseq_transcripts_having_highest_num_exons = merge( x=refseq_having_no_ucsc_canonical, y=refseq_highest_num_exons, by=c( "refseq_gene", "num_refseq_exons" ), all.x=FALSE, all.y=FALSE )
nrow(refseq_transcripts_having_highest_num_exons) # 1663. So the 881 genes have 1663 transcripts corresponding to the most number of exons for that gene. There are multiple longest transcripts for some genes.
refseq_transcripts_having_highest_num_exons_counts = sqldf("select refseq_gene, num_refseq_exons, count(*) from refseq_transcripts_having_highest_num_exons group by refseq_gene, num_refseq_exons")
colnames(refseq_transcripts_having_highest_num_exons_counts) = c( "refseq_gene", "num_refseq_exons", "num_transcripts_having_highest_num_exons" )
length(( refseq_transcripts_having_highest_num_exons_counts[ (refseq_transcripts_having_highest_num_exons_counts$num_transcripts_having_highest_num_exons==1), c("refseq_gene") ] )) # 480 of the 881 genes have only one longest transcript (in terms of number of exons)
length(( refseq_transcripts_having_highest_num_exons_counts[ (refseq_transcripts_having_highest_num_exons_counts$num_transcripts_having_highest_num_exons>1), c("refseq_gene") ] )) # 401 of the 881 genes have multiple longest transcripts (in terms of number of exons)

# Regardless of whether there are one or multiple longest transcripts, we will take all the multiple longest transcripts for these refseq genes that have multiple transcripts and don't have a UCSC canonical transcript.
refseq_transcripts_having_highest_num_exons = merge( x=refseq, y=refseq_highest_num_exons, by=c( "refseq_gene", "num_refseq_exons" ), all.x=FALSE, all.y=FALSE )
nrow(refseq_transcripts_having_highest_num_exons) # 1663
list_longest_refseq_transcripts_for_genes_having_more_than_one_refseq_transcript_and_no_canonical_ucsc_transcripts = unique(refseq_transcripts_having_highest_num_exons$refseq_transcript)
length(list_longest_refseq_transcripts_for_genes_having_more_than_one_refseq_transcript_and_no_canonical_ucsc_transcripts) ##### <===== 1367 refseq transcripts of genes having more than one refseq transcript and no ucsc canonical transcripts. These are the longest transcripts (highest number of exons) for these genes.

list_all_canonical_transcripts_for_all_genes = c( list_of_refseq_transcripts_for_genes_having_one_refseq_transcript, list_of_ucsc_canonical_transcripts_for_refseq_genes_having_more_than_one_refseq_transcript, list_longest_refseq_transcripts_for_genes_having_more_than_one_refseq_transcript_and_no_canonical_ucsc_transcripts )
length(list_all_canonical_transcripts_for_all_genes) # 27687
list_refseq_genes_for_all_canonical_transcripts_for_all_genes = unique(refseq[ (refseq$refseq_transcript %in% list_all_canonical_transcripts_for_all_genes), c("refseq_gene") ])
length(list_refseq_genes_for_all_canonical_transcripts_for_all_genes) # 27115. These 27687 canonical refseq transcripts that we will use cover all 27115 refseq genes.

# A quick check that all refseq genes are covered by these refseq transcripts.
length(unique( refseq[ (refseq$refseq_transcript %in% list_of_refseq_transcripts_for_genes_having_one_refseq_transcript), c("refseq_gene") ] )) # 15424
length(unique( refseq[ (refseq$refseq_transcript %in% list_of_ucsc_canonical_transcripts_for_refseq_genes_having_more_than_one_refseq_transcript), c("refseq_gene") ] )) # 10810
length(unique( refseq[ (refseq$refseq_transcript %in% list_longest_refseq_transcripts_for_genes_having_more_than_one_refseq_transcript_and_no_canonical_ucsc_transcripts), c("refseq_gene") ] )) # 881

# Now list the refseq exons.
all_canonical_refseq_transcripts = unique( refseq[ (refseq$refseq_transcript %in% list_all_canonical_transcripts_for_all_genes), ] )
nrow(all_canonical_refseq_transcripts) # 29157. Some transcripts must have multiple rows because there are only 27687 transcripts but they correspond to 29157 rows.

write.table( all_canonical_refseq_transcripts, file=output_refseq_transcripts_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE )



