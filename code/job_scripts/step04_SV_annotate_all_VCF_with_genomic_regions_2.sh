#!/bin/bash
set -euo pipefail

unique_id=$1
INDIR=$2
OUTDIR=$3

for infile in "${INDIR}"/*.vcf; do

	infile_basename=$(basename "${infile}")
	IFS=$'.' read -r -a array <<< "${infile_basename}"
	sampleid_and_more="${array[0]}"
	IFS=$'_' read -r -a array <<< "${sampleid_and_more}"
	sampleid="${array[0]}"
	outfile="${OUTDIR}"/"${sampleid}".gridss_annotated_genomic_regions_2.vcf
	unique_id_2="${unique_id}"_"${sampleid}"

	#if [ -e "${outfile}" ]; then
	#	do_nothing=1
	#else
		echo 'call step04_SV_annotate_one_VCF_with_genomic_regions_2.sh' "${unique_id_2}" "${infile}" "${outfile}"
		./step04_SV_annotate_one_VCF_with_genomic_regions_2.sh "${unique_id_2}" "${infile}" "${outfile}"
	#fi
done

