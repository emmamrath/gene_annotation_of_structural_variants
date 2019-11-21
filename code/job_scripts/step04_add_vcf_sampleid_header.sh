#!/bin/bash
set -euo pipefail

unique_id=$1
INDIR=$2
OUTDIR=$3

TMPDIR='/nvme/emmrat/tmp_SV'

for infile in "${INDIR}"/*.vcf; do

	infile_basename=$(basename "${infile}")
	infile_basename_prefix=${infile_basename::-4}
	IFS=$'_' read -r -a array <<< "${infile_basename}"
	sampleid="${array[0]}"
	outfile="${OUTDIR}"/"${infile_basename_prefix}"_sampleid.vcf

	echo 'Processing' $infile
	tmp_hdr="${TMPDIR}"/"${unique_id}"_"${sampleid}"_hdr.txt
	echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"$sampleid > "${tmp_hdr}"
	cat "${tmp_hdr}" "${infile}" > "${outfile}"
	rm "${tmp_hdr}"
done

