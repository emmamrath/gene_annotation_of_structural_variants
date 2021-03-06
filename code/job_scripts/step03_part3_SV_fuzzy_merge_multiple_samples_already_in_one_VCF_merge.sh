#!/bin/bash

# input variables
cohort_id=$1
merged_VCF_file_100bp="${cohort_id}"_SV_merged100bp.vcf
TMPDIR='/nvme/emmrat/tmp_SV'
OUTDIR='.'
SVannotation_software_directory='/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation'
Mobster_software_directory='/nvme/emmrat/software_for_mobile_elements/mobster_for_mobile_elements'
reference_genome='/nvme/emmrat/reference_genome_for_MGRB/hs37d5x/hs37d5x.fa'


echo 'Sort using linux sort'
# Sort the VCF data so that POS is sorted numerically within chromosome, so that white space after short POS are not affecting sort order.
# Non-numeric chromosomes will end up at the beginning of the file.
# Make sure sort key is -k2,2n or -k2n, not -nk2 or -nk2,2 because the latter causes global sorting of -n and will mix up X and Y variants
cat "${TMPDIR}"/"${cohort_id}"_SV_INS_merged100bp_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_DEL_merged100bp_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_INDEL_merged100bp_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_DUP_merged100bp_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_INV_merged100bp_nohdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_BND_merged100bp_nohdr.vcf | \
	sort -k1,1 -k2,2n > "${TMPDIR}"/"${cohort_id}"_SV_all_SVTYPE_merged100bp_sorted_nohdr.vcf

# Put the VCF header back onto the now sorted VCF data
grep '^#' "${TMPDIR}"/"${cohort_id}"_gridss_AllSamples_SV_INS_merged100bp.vcf > "${TMPDIR}"/"${cohort_id}"_SV_INS_merged100bp_hdr.vcf

cat "${TMPDIR}"/"${cohort_id}"_SV_INS_merged100bp_hdr.vcf "${TMPDIR}"/"${cohort_id}"_SV_all_SVTYPE_merged100bp_sorted_nohdr.vcf \
        > "${TMPDIR}"/"${cohort_id}"_SV_all_SVTYPE_merged100bp_sorted.vcf

# Now sort the VCF file by the chromosome sort order of the reference genome.
# Non-numeric chromosomes will appear not in alphabetical order, and instead will appear at the end as X, Y, MT, viruses, etc.
echo 'Sort using chromosome sort order'
python "${Mobster_software_directory}"/extra_programs_for_Mobster/sort_VCF_by_order_found_in_input_list.py \
        -i "${TMPDIR}"/"${cohort_id}"_SV_all_SVTYPE_merged100bp_sorted.vcf \
        -l "${Mobster_software_directory}"/extra_data_for_Mobster/hs37d5x_sort_order.txt \
        -o "${OUTDIR}"/"${merged_VCF_file_100bp}"

echo 'Finished'


