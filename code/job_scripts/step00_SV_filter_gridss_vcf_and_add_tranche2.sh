#!/bin/bash

INDIR="../gridss_results"
OUTDIR="../gridss_results_TRANCHE2"
sv_pgms_dir='/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation/programs'

for infile in "${INDIR}"/*.gridss.vcf; do

	IFS='/' read -ra array1 <<< "$infile"
	filename="${array1[-1]}"
	IFS='.' read -ra array2 <<< "$filename"
	sampleid="${array2[0]}"
	outfile="${OUTDIR}"/"${sampleid}".gridss_BND_filtered.vcf

	if [ -e "${infile}" ]; then
		if [ -e "${outfile}" ]; then
			do_nothing=1
		else
			echo 'cat '"${infile}"' | python '"${sv_pgms_dir}"'/add_tranche2_and_filter_gridss_vcf.py -max_BND_length_for_no_split_reads 100 > '"${outfile}"
			cat "${infile}" | python "${sv_pgms_dir}"/add_tranche2_and_filter_gridss_vcf.py -max_BND_length_for_no_split_reads 100 > "${outfile}"
		fi
	else
		echo 'file not found for' "${infile}"
	fi

done

