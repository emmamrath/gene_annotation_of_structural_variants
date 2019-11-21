#!/bin/bash

infile=$1 # "../vcf_INS_DEL_INDEL_DUP_INV_BND_annotated/AAAAA_SV.vcf"
OUTDIR=$2 # "../vcf_INS_DEL_INDEL_DUP_INV_BND_annotated_consequences"

#module load python/2.7.6

infile_basename=$(basename "${infile}")
infile_basename_prefix=${infile_basename::-4}
IFS=$'_' read -r -a array <<< "${infile_basename}"
sampleid="${array[0]}"
outfile="${OUTDIR}"/"${infile_basename_prefix}"_consequences.vcf

if [ -e "${outfile}" ]; then
	echo 'sampleid' $sampleid 'already has output' $outfile
else
	echo 'call step04_SV_one_sample_annotate_gridss_VCF_with_genomic_regions_consequences.sh' "${sampleid}" "${infile}" "${outfile}"
	./step04_parallel_SV_one_sample_annotate_gridss_VCF_with_genomic_regions_consequences.sh "${sampleid}" "${infile}" "${outfile}"
fi
echo ''

