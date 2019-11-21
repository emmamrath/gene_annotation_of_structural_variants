#!/bin/bash

# input variables
cohort_id=$1
SVtype=$2
chr="nonNumChr"

merged_VCF_file_100bp="${cohort_id}"_SV_merged100bp.vcf
TMPDIR='/nvme/emmrat/tmp_SV'
OUTDIR='.'
SVannotation_software_directory='/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation'
Mobster_software_directory='/nvme/emmrat/software_for_mobile_elements/mobster_for_mobile_elements'
reference_genome='/nvme/emmrat/reference_genome_for_MGRB/hs37d5x/hs37d5x.fa'

grep '^#' "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_"${SVtype}"_only.vcf > "${TMPDIR}"/"${cohort_id}"_gridss_AllSamples_"${SVtype}"_only_"${chr}"_hdr.vcf
grep -v -P '^1' "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_"${SVtype}"_only.vcf | grep -v -P '^2' | grep -v -P '^3' | grep -v -P '^4' | grep -v -P '^5' | grep -v -P '^6' | grep -v -P '^7' | grep -v -P '^8' | grep -v -P '^9' | grep -v -P '^10' | grep -v -P '^11' | grep -v -P '^12' | grep -v -P '^13' | grep -v -P '^14' | grep -v -P '^15' | grep -v -P '^16' | grep -v -P '^17' | grep -v -P '^18' | grep -v -P '^19' | grep -v -P '^20' | grep -v -P '^21' | grep -v -P '^22' | grep -v -P '^X' | grep -v -P '^Y' | grep -v -P '^M' > "${TMPDIR}"/"${cohort_id}"_gridss_AllSamples_"${SVtype}"_only_"${chr}"_body.vcf

echo 'For ' $SVtype ' records, run merge_a_multi_sample_SV_VCF_for_similar_position_variants.py'
cat "${TMPDIR}"/"${cohort_id}"_gridss_AllSamples_"${SVtype}"_only_"${chr}"_hdr.vcf "${TMPDIR}"/"${cohort_id}"_gridss_AllSamples_"${SVtype}"_only_"${chr}"_body.vcf \
	| python "${SVannotation_software_directory}"/programs/merge_a_multi_sample_SV_VCF_for_similar_position_variants.py \
	100 "${reference_genome}" > "${TMPDIR}"/"${cohort_id}"_gridss_AllSamples_SV_"${SVtype}"_merged100bp_"${chr}".vcf

echo 'Remove header to prepare for merge'
grep -v '^#' "${TMPDIR}"/"${cohort_id}"_gridss_AllSamples_SV_"${SVtype}"_merged100bp_"${chr}".vcf > "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_"${chr}"_nohdr.vcf

echo 'Finished'


