#!/bin/bash

# input variables
cohort_id=$1
SVtype=$2
TMPDIR='/nvme/emmrat/tmp_SV'
OUTDIR='.'

cat "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr1_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr10_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr11_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr12_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr13_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr14_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr15_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr16_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr17_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr18_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr19_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr2_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr20_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr21_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr22_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr3_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr4_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr5_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr6_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr7_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr8_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr9_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chrX_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chrY_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chrM_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_nonNumChr_nohdr.vcf > "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_nohdr.vcf

grep '^#' "${TMPDIR}"/"${cohort_id}"_gridss_AllSamples_SV_"${SVtype}"_merged100bp_chr1.vcf > "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr1_hdr.vcf

cat "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_chr1_hdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_nohdr.vcf > "${TMPDIR}"/"${cohort_id}"_gridss_AllSamples_SV_"${SVtype}"_merged100bp.vcf

echo 'Finished'


