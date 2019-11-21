#!/bin/bash

# input variables
cohort_id=$1
SVtype=$2

TMPDIR='/nvme/emmrat/tmp_SV'
OUTDIR='.'
SVannotation_software_directory='/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation'
Mobster_software_directory='/nvme/emmrat/software_for_mobile_elements/mobster_for_mobile_elements'
reference_genome='/nvme/emmrat/reference_genome_for_MGRB/hs37d5x/hs37d5x.fa'

echo 'For ' $SVtype ' records, run merge_a_multi_sample_SV_VCF_for_similar_position_variants.py'
cat "${cohort_id}"_gridss_AllSamples_"${SVtype}"_only.vcf \
	| python "${SVannotation_software_directory}"/programs/merge_a_multi_sample_SV_VCF_for_similar_position_variants.py \
	100 "${reference_genome}" > "${TMPDIR}"/"${cohort_id}"_gridss_AllSamples_SV_"${SVtype}"_merged100bp.vcf

echo 'Remove header to prepare for merge'
grep -v '^#' "${TMPDIR}"/"${cohort_id}"_gridss_AllSamples_SV_"${SVtype}"_merged100bp.vcf > "${TMPDIR}"/"${cohort_id}"_SV_"${SVtype}"_merged100bp_nohdr.vcf

echo 'Finished'


