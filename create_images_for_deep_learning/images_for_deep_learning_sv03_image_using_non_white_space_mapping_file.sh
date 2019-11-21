#!/bin/bash
set -euo pipefail

module load python3/3.5.2
module load python3/3.5.2-matplotlib
module unload intel-fc intel-cc
module load gcc/9.1.0

export PYTHONPATH=/g/data3/wq2/users/emr913/software_various/python_libraries/python3.5.2/lib/python3.5/site-packages

QUEUE=normal
INDIR="../vcf_TRANCHE2_annotations1_annotations2"
OUTDIR="../images_for_deep_learning"
map_file="/g/data3/zc9/users/emr913/ISKS_2018jan/hs37d5x/gridss_SV/images_for_deep_learning/non_white_out_all_mgrb_and_isks_map_350x350.txt"

cohort="ISKS"

MAXQUEUED=100
num_queued=$((MAXQUEUED+1))

COUNTER=1
for infile in "${INDIR}"/*.tsv; do

  sampleid=$(basename "${infile}")
  sampleid="${sampleid%%.*}"
  outplot="${OUTDIR}"/"${sampleid}"_SV_350x350.png

  echo 'python3 images_for_deep_learning_sv03_image_using_non_white_space_mapping_file.py -i' "${infile}" '-o' "${outplot}" '-m' "${map_file}"
  python3 images_for_deep_learning_sv03_image_using_non_white_space_mapping_file.py -i "${infile}" -o "${outplot}" -m "${map_file}"

done

