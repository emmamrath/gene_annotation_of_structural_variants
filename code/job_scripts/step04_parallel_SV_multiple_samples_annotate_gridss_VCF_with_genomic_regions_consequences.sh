#!/bin/bash

INDIR=$1 # "../vcf_INS_DEL_INDEL_DUP_INV_BND_annotated"
OUTDIR=$2 # "../vcf_INS_DEL_INDEL_DUP_INV_BND_annotated_consequences"

mkdir -p "${OUTDIR}"

#module load python/2.7.6

ls -1 ../vcf_INS_DEL_INDEL_DUP_INV_BND_annotated/*.vcf | parallel -j10 bash ./step04_parallel_SV_one_sample_annotate_gridss_VCF_with_genomic_regions_consequences.sh {} "${OUTDIR}"

echo 'Finished'

