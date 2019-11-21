#!/bin/bash

INDIR=$1
OUTDIR=$2

# ./step04_convert_VCF_INFO_to_tab_delimited_and_intact_fields.sh ../vcf_TRANCHE2_INS_DEL_INDEL_DUP_INV_BND ../vcf_TRANCHE2_INS_DEL_INDEL_DUP_INV_BND

SVCODEDIR='/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation/programs'

for infile in "${INDIR}"/*.vcf; do

	infile_basename=$(basename "${infile}")
	IFS=$'.' read -r -a array <<< "${infile_basename}"
	sampleid="${array[0]}"
	outfile_basename="${infile_basename/%.vcf/.tsv}"
	outfile="${OUTDIR}"/"${outfile_basename}"

	echo 'python' "${SVCODEDIR}"'/convert_VCF_INFO_one_sample_to_tab_delimited_and_intact_fields.py -i '"${infile}"' -o '"${outfile}"
	python "${SVCODEDIR}"/convert_VCF_INFO_one_sample_to_tab_delimited_and_intact_fields.py -i "${infile}" -o "${outfile}"
done

