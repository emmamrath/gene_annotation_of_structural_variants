#!/bin/bash

unique_id=$1
INDIR="../gridss_results"
OUTDIR='../vcf_TRANCHE2_INS_DEL_INDEL_DUP_INV_BND'

for infile in "${INDIR}"/*.vcf; do

    # sampleid=${infile:0:5}
    infile2=$(basename "${infile}")
    IFS='.' read -ra infile_array <<< "$infile2"
    sampleid="${infile_array[0]}"
    outfile="${OUTDIR}"/"${sampleid}".gridss_BNDVAF_INS_DEL_INDEL_DUP_INV_BND.vcf
    unique_id_2="${unique_id}"_"${sampleid}"

    if [ -e "${outfile}" ]; then
        do_nothing=1
    else
        echo 'call step01_SV_one_sample_call_INS_DEL_INDEL_DUP_INV_from_gridss_output.sh for' "${unique_id_2}"
        ./step01_SV_one_sample_call_INS_DEL_INDEL_DUP_INV_from_gridss_output.sh "${unique_id_2}" "${infile}" "${outfile}"
    fi
done

echo 'Finished'

