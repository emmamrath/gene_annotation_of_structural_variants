#!/bin/bash
set -euo pipefail

unique_id=$1
INDIR='../vcf_TRANCHE2_INS_DEL_INDEL_DUP_INV_BND'
OUTDIR='../vcf_TRANCHE2_INS_DEL_INDEL_DUP_INV_BND_annotated'

for infile in "${INDIR}"/*.vcf; do

	infile_basename=$(basename "${infile}")
	IFS=$'.' read -r -a array <<< "${infile_basename}"
	sampleid="${array[0]}"
	outfile="${OUTDIR}"/"${sampleid}".gridss_BNDVAF_INS_DEL_INDEL_DUP_INV_BND_annotated.vcf
	unique_id_2="${unique_id}"_"${sampleid}"

	if [ -e "${outfile}" ]; then
		do_nothing=1
	else
		echo 'call step04_SV_annotate_one_VCF_with_genomic_regions.sh' "${unique_id_2}" "${infile}" "${outfile}"
		./step04_SV_annotate_one_VCF_with_genomic_regions.sh "${unique_id_2}" "${infile}" "${outfile}"
	fi
done

