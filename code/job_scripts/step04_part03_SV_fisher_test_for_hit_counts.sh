#!/bin/bash

cohort1="MGRB"
indir1="/nvme/emmrat/MGRB_2017/gridss_SV_2017sep/vcf_INS_DEL_INDEL_DUP_INV_BND_annotated_bicseq"
cohort2="ISKS"
indir2="/nvme/emmrat/ISKS_2018jan/hs37d5x/gridss_SV/vcf_INS_DEL_INDEL_DUP_INV_BND_annotated_bicseq"
cohort_out="MGRBvsISKS"
outdir="/nvme/emmrat/MGRB_2017/gridss_SV_2017sep/analysis_MGRB_vs_ISKS_for_INS_DEL_INDEL_DUP_INV_BND_annotated_bicseq"
sv_type="DELINDEL"

SVCODEDIR='/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation/programs'
TMPDIR='/nvme/emmrat/tmp_SV'

num_mgrb="$(ls -1 ${indir1}/*_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | wc -l)"
num_isks="$(ls -1 ${indir2}/*_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | wc -l)"
param_sv_type='DEL/INDEL'
param_sv_hit_type='HITEXONS'

########## Look at gridss calls that hit genes, and that have a matching bicseq call

echo 'Look at gridss calls that hit genes, and that have a matching bicseq call'

hit_type="HITEXONS"
bicseq_type="matchBicseq"
unique_id="${cohort_out}_${sv_type}_${hit_type}_${bicseq_type}"
outfile="${outdir}/${unique_id}_fisher_tests.tsv"

echo '   create input files'

# subset the MGRB data to look at
cat "${indir1}"/ZAAAA_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | head -n 1 | cut -d$'\t' -f1-2,11,13-16,32 > "${TMPDIR}/${unique_id}"_temp1_hits_hdr.txt
cat "${indir1}"/*_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | grep -P 'INTERMEDIATE|HIGH' | grep -P 'DEL|INDEL' | awk '($16 != 0) ' | cut -d$'\t' -f1-2,11,13-16,32 | awk '($8 != "")' > "${TMPDIR}/${unique_id}"_temp1_hits_body.txt # includes non-bicseq matches of gridss
cat "${TMPDIR}/${unique_id}"_temp1_hits_hdr.txt "${TMPDIR}/${unique_id}"_temp1_hits_body.txt > "${TMPDIR}/${unique_id}"_temp1_hits_input.txt

# subset the ISKS data to look at
cat "${indir2}"/isks1002_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | head -n 1 | cut -d$'\t' -f1-2,11,13-16,32 > "${TMPDIR}/${unique_id}"_temp2_hits_hdr.txt
cat "${indir2}"/*_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | grep -P 'INTERMEDIATE|HIGH' | grep -P 'DEL|INDEL' | awk '($16 != 0) ' | cut -d$'\t' -f1-2,11,13-16,32 | awk '($8 != "")' > "${TMPDIR}/${unique_id}"_temp2_hits_body.txt # includes non-bicseq matches of gridss
cat "${TMPDIR}/${unique_id}"_temp2_hits_hdr.txt "${TMPDIR}/${unique_id}"_temp2_hits_body.txt > "${TMPDIR}/${unique_id}"_temp2_hits_input.txt

echo '   call Rscript SV_fisher_test_for_hit_counts.R'

# run fisher test for each gene
Rscript "${SVCODEDIR}"/SV_fisher_test_for_hit_counts.R "${TMPDIR}/${unique_id}"_temp1_hits_input.txt "${TMPDIR}/${unique_id}"_temp2_hits_input.txt "${outfile}" "${num_mgrb}" "${num_isks}" "${param_sv_type}" "${param_sv_hit_type}"

rm -f "${TMPDIR}/${unique_id}"_*

########## Look at gridss calls that hit genes, regardless of whether they have a matching bicseq call

echo 'Look at gridss calls that hit genes, regardless of whether they have a matching bicseq call'

hit_type="HITEXONS"
bicseq_type="matchOrNotBicseq"
unique_id="${cohort_out}_${sv_type}_${hit_type}_${bicseq_type}"
outfile="${outdir}/${unique_id}_fisher_tests.tsv"

echo '   create input files'

# subset the MGRB data to look at
cat "${indir1}"/ZAAAA_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | head -n 1 | cut -d$'\t' -f1-2,11,13-16,32 > "${TMPDIR}/${unique_id}"_temp1_hits_hdr.txt
cat "${indir1}"/*_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | grep -P 'INTERMEDIATE|HIGH' | grep -P 'DEL|INDEL' | cut -d$'\t' -f1-2,11,13-16,32 | awk '($8 != "")' > "${TMPDIR}/${unique_id}"_temp1_hits_body.txt # includes non-bicseq matches of gridss
cat "${TMPDIR}/${unique_id}"_temp1_hits_hdr.txt "${TMPDIR}/${unique_id}"_temp1_hits_body.txt > "${TMPDIR}/${unique_id}"_temp1_hits_input.txt

# subset the ISKS data to look at
cat "${indir2}"/isks1002_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | head -n 1 | cut -d$'\t' -f1-2,11,13-16,32 > "${TMPDIR}/${unique_id}"_temp2_hits_hdr.txt
cat "${indir2}"/*_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | grep -P 'INTERMEDIATE|HIGH' | grep -P 'DEL|INDEL' | cut -d$'\t' -f1-2,11,13-16,32 | awk '($8 != "")' > "${TMPDIR}/${unique_id}"_temp2_hits_body.txt # includes non-bicseq matches of gridss
cat "${TMPDIR}/${unique_id}"_temp2_hits_hdr.txt "${TMPDIR}/${unique_id}"_temp2_hits_body.txt > "${TMPDIR}/${unique_id}"_temp2_hits_input.txt

echo '   call Rscript SV_fisher_test_for_hit_counts.R'

# run fisher test for each gene
Rscript "${SVCODEDIR}"/SV_fisher_test_for_hit_counts.R "${TMPDIR}/${unique_id}"_temp1_hits_input.txt "${TMPDIR}/${unique_id}"_temp2_hits_input.txt "${outfile}" "${num_mgrb}" "${num_isks}" "${param_sv_type}" "${param_sv_hit_type}"

rm -f "${TMPDIR}/${unique_id}"_*

########## Look at gridss calls that have a matching bicseq call, regardless of whether they hit genes

echo 'Look at gridss calls that have a matching bicseq call, regardless of whether they hit genes'

hit_type="hitGenesOrNot"
bicseq_type="matchBicseq"
unique_id="${cohort_out}_${sv_type}_${hit_type}_${bicseq_type}"
outfile="${outdir}/${unique_id}_fisher_tests.tsv"

echo '   create input files'

# subset the MGRB data to look at
cat "${indir1}"/ZAAAA_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | head -n 1 | cut -d$'\t' -f1-2,11,13-16,32 > "${TMPDIR}/${unique_id}"_temp1_hits_hdr.txt
cat "${indir1}"/*_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | grep -P 'INTERMEDIATE|HIGH' | grep -P 'DEL|INDEL' | awk '($16 != 0) ' | cut -d$'\t' -f1-2,11,13-16,32 > "${TMPDIR}/${unique_id}"_temp1_hits_body.txt # includes non-bicseq matches of gridss
cat "${TMPDIR}/${unique_id}"_temp1_hits_hdr.txt "${TMPDIR}/${unique_id}"_temp1_hits_body.txt > "${TMPDIR}/${unique_id}"_temp1_hits_input.txt

# subset the ISKS data to look at
cat "${indir2}"/isks1002_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | head -n 1 | cut -d$'\t' -f1-2,11,13-16,32 > "${TMPDIR}/${unique_id}"_temp2_hits_hdr.txt
cat "${indir2}"/*_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | grep -P 'INTERMEDIATE|HIGH' | grep -P 'DEL|INDEL' | awk '($16 != 0) ' | cut -d$'\t' -f1-2,11,13-16,32 > "${TMPDIR}/${unique_id}"_temp2_hits_body.txt # includes non-bicseq matches of gridss
cat "${TMPDIR}/${unique_id}"_temp2_hits_hdr.txt "${TMPDIR}/${unique_id}"_temp2_hits_body.txt > "${TMPDIR}/${unique_id}"_temp2_hits_input.txt

echo '   call Rscript SV_fisher_test_for_hit_counts.R'

# run fisher test for each gene
#Rscript "${SVCODEDIR}"/SV_fisher_test_for_hit_counts.R "${TMPDIR}/${unique_id}"_temp1_hits_input.txt "${TMPDIR}/${unique_id}"_temp2_hits_input.txt "${outfile}" "${num_mgrb}" "${num_isks}" "${param_sv_type}" "${param_sv_hit_type}"

#rm -f "${TMPDIR}/${unique_id}"_*

########## Look at gridss calls, regardless of whether they have a matching bicseq call, regardless of whether they hit genes

echo 'Look at gridss calls, regardless of whether they have a matching bicseq call, regardless of whether they hit genes'

hit_type="hitGenesOrNot"
bicseq_type="matchOrNotBicseq"
unique_id="${cohort_out}_${sv_type}_${hit_type}_${bicseq_type}"
outfile="${outdir}/${unique_id}_fisher_tests.tsv"

echo '   create input files'

# subset the MGRB data to look at
cat "${indir1}"/ZAAAA_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | head -n 1 | head -n 1 | cut -d$'\t' -f1-2,11,13-16,32 > "${TMPDIR}/${unique_id}"_temp1_hits_hdr.txt
cat "${indir1}"/*_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | grep -P 'INTERMEDIATE|HIGH' | grep -P 'DEL|INDEL' | cut -d$'\t' -f1-2,11,13-16,32 > "${TMPDIR}/${unique_id}"_temp1_hits_body.txt # includes non-bicseq matches of gridss
cat "${TMPDIR}/${unique_id}"_temp1_hits_hdr.txt "${TMPDIR}/${unique_id}"_temp1_hits_body.txt > "${TMPDIR}/${unique_id}"_temp1_hits_input.txt

# subset the ISKS data to look at
cat "${indir2}"/isks1002_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | head -n 1 | cut -d$'\t' -f1-2,11,13-16,32 > "${TMPDIR}/${unique_id}"_temp2_hits_hdr.txt
cat "${indir2}"/*_SV_INS_DEL_DUP_INV_annotated_bicseqOverlaps.tsv | grep -P 'INTERMEDIATE|HIGH' | grep -P 'DEL|INDEL' | cut -d$'\t' -f1-2,11,13-16,32 > "${TMPDIR}/${unique_id}"_temp2_hits_body.txt # includes non-bicseq matches of gridss
cat "${TMPDIR}/${unique_id}"_temp2_hits_hdr.txt "${TMPDIR}/${unique_id}"_temp2_hits_body.txt > "${TMPDIR}/${unique_id}"_temp2_hits_input.txt

echo '   call Rscript SV_fisher_test_for_hit_counts.R'

# run fisher test for each gene
#Rscript "${SVCODEDIR}"/SV_fisher_test_for_hit_counts.R "${TMPDIR}/${unique_id}"_temp1_hits_input.txt "${TMPDIR}/${unique_id}"_temp2_hits_input.txt "${outfile}" "${num_mgrb}" "${num_isks}" "${param_sv_type}" "${param_sv_hit_type}"

#rm -f "${TMPDIR}/${unique_id}"_*

echo 'Finished'

