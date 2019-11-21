#!/bin/bash
set -euo pipefail

INDIR="/g/data3/wq2/users/emr913/Emma_mobile_element_investigations/results/ISKS_results"
OUTDIR=$INDIR

REF="/g/data3/wq2/users/emr913/Emma_mobile_element_investigations/reference_data/GRCh37_Swergold_1990_Speek_2001_Line1_promoters_GenesHitInFirstCodingExon.bed"

for infile in "${INDIR}"/*_DEL_INDEL.bed; do
  filename=$(basename "${infile}")
  IFS='_' read -r -a array <<< "$filename"
  sampleid="${array[0]}"
  outfile="${OUTDIR}"/"${sampleid}"_DEL_INDEL_intersects_Swergold_1990_Speek_2001_Line1_promoters.bed
  bedtools intersect -a "${infile}" -b "${REF}" -wo -F 0.90 | awk 'BEGIN {FS="\t";OFS="\t"} {print $0, $9-$8}' > "${outfile}"
done


