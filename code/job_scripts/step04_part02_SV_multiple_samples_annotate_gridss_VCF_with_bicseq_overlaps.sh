#!/bin/bash

cohortid=$1
INDIR_SV=$2 # "../vcf_TRANCHE2_INS_DEL_INDEL_DUP_INV_BND_annotated"
OUTDIR=$3 # "../vcf_TRANCHE2_INS_DEL_INDEL_DUP_INV_BND_annotated_bicseq"
INDIR_CNV=$4 # "/nvme/emmrat/MGRB_2017/bicseq_CNV/VCF_files_per_sample_filtered_pvalue"

mkdir -p "${OUTDIR}"

outfile_cohort_overlap_stats_tsv="${OUTDIR}"/"${cohortid}"_overlap_stats_for_gridss_and_bicseq.tsv

#module load python/2.7.6

##### For each sample, annotate the gridss SV file with bicseq CNV overlaps from bicseq CNV file

for infile_sv in "${INDIR_SV}"/*.vcf; do

    infile2=$(basename "${infile_sv}")
    IFS='.' read -ra infile_array <<< "$infile2"
    sampleid="${infile_array[0]}"

    infile_cnv="${INDIR_CNV}"/"${sampleid}".bicseq_filtered.vcf

    outfile="${OUTDIR}"/"${sampleid}"_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.vcf
    outfile_tsv="${OUTDIR}"/"${sampleid}"_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv
    outfile_overlap_stats_tsv="${OUTDIR}"/"${sampleid}"_overlap_stats_for_gridss_and_bicseq.tsv

    if [ -e "${infile_cnv}" ]; then # if the bicseq CNV file exists, then incorporate it into the gridss SV file

        echo 'sampleid' $sampleid 'outfile' $outfile
        #if [ -e "${outfile}" ]; then
        #    : # outfile has already been created
        #else

            echo 'call step04_part02_SV_one_sample_annotate_gridss_VCF_with_bicseq_overlaps.sh' "${cohortid}" "${sampleid}" "${infile_sv}" "${infile_cnv}" "${outfile}" "${outfile_tsv}"
            ./step04_part02_SV_one_sample_annotate_gridss_VCF_with_bicseq_overlaps.sh "${cohortid}" "${sampleid}" "${infile_sv}" "${infile_cnv}" "${outfile}" "${outfile_tsv}" "${outfile_overlap_stats_tsv}"
        #fi

    else
        echo 'sampleid' $sampleid 'has no bicseq file'
    fi
    echo ''
done

# Concatenate stats files for all samples into one file for whole cohort

rm -f "${outfile_cohort_overlap_stats_tsv}"
cat "${OUTDIR}"/*_overlap_stats_for_gridss_and_bicseq.tsv | grep '^sample' | head -n 1 > temp_stats_hdr.tsv
cat "${OUTDIR}"/*_overlap_stats_for_gridss_and_bicseq.tsv | grep -v '^sample'| sort > temp_stats_body.tsv
cat temp_stats_hdr.tsv temp_stats_body.tsv > "${outfile_cohort_overlap_stats_tsv}"
rm temp_stats_hdr.tsv temp_stats_body.tsv

echo 'Finished'
