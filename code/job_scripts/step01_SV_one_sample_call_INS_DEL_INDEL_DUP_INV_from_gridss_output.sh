#!/bin/bash
# set -euo pipefail 

# input variables
TMPDIR="/nvme/emmrat/tmp_sv"
SAMPLEID=$1
INFILE=$2
OUTFILE=$3
gridss_vcf_file="${INFILE}"
output_SV_vcf_file="${OUTFILE}"
SVannotation_software_directory='/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation'
GAP=3000
cohort_id=$SAMPLEID

##################################################
# Filter out low quality SV calls and add TRANCHE2 quality value
##################################################

cat "${gridss_vcf_file}" | python "${SVannotation_software_directory}"/programs/add_tranche2_and_filter_gridss_vcf.py -max_BND_length_for_no_split_reads 100 \
	> "${TMPDIR}/${SAMPLEID}_gridss_BND_filtered.vcf"

##################################################
# Convert GRIDSS BND records into <DEL>, <INS>, <INDEL>, <DUP:TANDEM>, <DUP:INS>
##################################################
#echo 'For' ${SAMPLEID} 'convert GRIDSS BND records into <DEL>, <INS>, <INDEL>, <DUP:TANDEM>, <DUP:INS>'

# sort input by GRIDSS EVENT
grep '^#' "${TMPDIR}/${SAMPLEID}_gridss_BND_filtered.vcf" > "${TMPDIR}/${SAMPLEID}_gridss_hdr.vcf"
grep -v '^#' "${TMPDIR}/${SAMPLEID}_gridss_BND_filtered.vcf" | sort -t$'\t' -k 3,3 > "${TMPDIR}/${SAMPLEID}_gridss_nohdr_sortedByEvent.txt"
cat "${TMPDIR}/${SAMPLEID}_gridss_hdr.vcf" "${TMPDIR}/${SAMPLEID}_gridss_nohdr_sortedByEvent.txt" > "${TMPDIR}/${SAMPLEID}_gridss_sortedByEvent.vcf"

# run program to identify <DEL>, <INS>, <INDEL>, <DUP:TANDEM>, <DUP:INS>
cat "${TMPDIR}/${SAMPLEID}_gridss_sortedByEvent.vcf" \
	| python "${SVannotation_software_directory}"/programs/convert_VCF_BND_records_to_simple_structural_variant_VCF_records.py \
	> "${TMPDIR}/${SAMPLEID}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP_includes_LowComplexity.vcf"   

# Remove SV variants having one or both breakends fall in a region of low-complexity or simple-repeats.
# They are probably mapped incorrectly and so the SV is false.
# These regions are the low-complexity or simple-repeat regions defined in UCSC Repeats file, not the entire Repeats file of mobile elements and other repeats.
grep '^#' "${TMPDIR}/${SAMPLEID}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP_includes_LowComplexity.vcf" > "${TMPDIR}/${cohort_id}_VCF_hdr.vcf"
grep -v '^#' "${TMPDIR}/${SAMPLEID}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP_includes_LowComplexity.vcf" | sort -k1,1 -k2,2n > "${TMPDIR}/${cohort_id}_VCF_body.vcf"

cat "${TMPDIR}/${cohort_id}_VCF_body.vcf" | \
	python "${SVannotation_software_directory}"/programs/add_breakend_key_to_vcf_record.py -p LEFT > "${TMPDIR}/${cohort_id}_add_key_LEFT.txt"
bedtools intersect -a "${TMPDIR}/${cohort_id}_add_key_LEFT.txt" \
	-b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_Repeats_20171007_LowComplexity_SimpleRepeats.bed" -v | cut -d$'\t' -f4- > "${TMPDIR}/${cohort_id}_add_key_LEFT_removed.vcf"
cat "${TMPDIR}/${cohort_id}_add_key_LEFT_removed.vcf" | \
	python "${SVannotation_software_directory}"/programs/add_breakend_key_to_vcf_record.py -p RIGHT > "${TMPDIR}/${cohort_id}_add_key_RIGHT.txt"
bedtools intersect -a "${TMPDIR}/${cohort_id}_add_key_RIGHT.txt" \
	-b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_Repeats_20171007_LowComplexity_SimpleRepeats.bed" -v | cut -d$'\t' -f4- > "${TMPDIR}/${cohort_id}_removed_low_complexity_body.vcf"
cat "${TMPDIR}/${cohort_id}_VCF_hdr.vcf" "${TMPDIR}/${cohort_id}_removed_low_complexity_body.vcf" > "${TMPDIR}/${cohort_id}_removed_low_complexity.vcf"
rm "${TMPDIR}/${cohort_id}_VCF_body.vcf"
rm "${TMPDIR}/${cohort_id}_add_key_LEFT.txt" "${TMPDIR}/${cohort_id}_add_key_LEFT_removed.vcf"
rm "${TMPDIR}/${cohort_id}_add_key_RIGHT.txt" "${TMPDIR}/${cohort_id}_removed_low_complexity_body.vcf"

cat "${TMPDIR}/${cohort_id}_VCF_hdr.vcf" "${TMPDIR}/${cohort_id}_removed_low_complexity.vcf" > "${TMPDIR}/${SAMPLEID}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP.vcf"

# sort output by CHROM+POS 
grep -v '^#' "${TMPDIR}/${SAMPLEID}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP.vcf" > "${TMPDIR}/${SAMPLEID}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP_nohdr.txt"
sort -t$'\t' -k1,1 -k2,2n "${TMPDIR}/${SAMPLEID}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP_nohdr.txt" > "${TMPDIR}/${SAMPLEID}_gridss_BND_INS_DEL_INDEL_DUP_nohdr_sorted.txt"

# split output into the new INS_DEL_INDEL_DUP records and the remaining old BND records not used to call the INS_DEL_INDEL_DUP
grep 'SVTYPE=BND' "${TMPDIR}/${SAMPLEID}_gridss_BND_INS_DEL_INDEL_DUP_nohdr_sorted.txt" > "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP_nohdr_sorted.txt"
grep -v 'SVTYPE=BND' "${TMPDIR}/${SAMPLEID}_gridss_BND_INS_DEL_INDEL_DUP_nohdr_sorted.txt" > "${TMPDIR}/${SAMPLEID}_gridss_INS_DEL_INDEL_DUP_nohdr_sorted.txt"
cat "${TMPDIR}/${SAMPLEID}_gridss_hdr.vcf" "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP_nohdr_sorted.txt" > "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP.vcf"

# identify new VCF headers added by INS_DEL_INDEL_DUP program
grep '^#' "${TMPDIR}/${SAMPLEID}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP.vcf" > "${TMPDIR}/${SAMPLEID}_gridss_BND_INS_DEL_INDEL_DUP_hdr.vcf"
sort "${TMPDIR}/${SAMPLEID}_gridss_hdr.vcf" | uniq > "${TMPDIR}/${SAMPLEID}_gridss_hdr_sorted.txt"
sort "${TMPDIR}/${SAMPLEID}_gridss_BND_INS_DEL_INDEL_DUP_hdr.vcf" | uniq > "${TMPDIR}/${SAMPLEID}_gridss_BND_INS_DEL_INDEL_DUP_hdr_sorted.txt"
comm -23 "${TMPDIR}/${SAMPLEID}_gridss_BND_INS_DEL_INDEL_DUP_hdr_sorted.txt" "${TMPDIR}/${SAMPLEID}_gridss_hdr_sorted.txt" > "${TMPDIR}/${SAMPLEID}_gridss_INS_DEL_INDEL_DUP_new_hdrs_only.txt"

##################################################
# Convert GRIDSS BND records into <INV> and other inversion records
##################################################
#echo 'For' ${SAMPLEID} 'convert GRIDSS BND records into <INV> and other inversion records'

# sort input by special key needed to put BND records of one inversion next to each other, in preparation for identifying inversions
python "${SVannotation_software_directory}"/programs/add_key_to_VCF_to_find_inversions.py \
	-i "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP.vcf" \
	-o "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP_with_key.txt"
sort -k1,1 -k2,2n -k3,3 -k4,4n "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP_with_key.txt" > "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP_with_key_sorted.txt"

# run program to identify inversions, but don't yet identify circular-inversions
python "${SVannotation_software_directory}"/programs/identify_inversions_in_VCF_with_key.py \
	-i "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP_with_key_sorted.txt" \
	-o "${TMPDIR}/${SAMPLEID}_gridss_inversions_only.vcf" \
	-u "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_inversions.vcf" \
	-gap "${GAP}"

# sort input by special key needed to put BND records of one inversion next to each other, in preparation for identifying circular-inversions
python "${SVannotation_software_directory}"/programs/add_key_to_VCF_to_find_inversions.py \
	-i "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_inversions.vcf" \
	-o "${TMPDIR}/${SAMPLEID}_gridss_unused_breakend_records_with_key.txt"
sort -k1,1 -k2,2n -k3,3 -k4,4n "${TMPDIR}/${SAMPLEID}_gridss_unused_breakend_records_with_key.txt" > "${TMPDIR}/${SAMPLEID}_gridss_unused_breakend_records_with_key_sorted.txt"
# run program to identify circular-inversions
python "${SVannotation_software_directory}"/programs/identify_inversions_in_VCF_with_key.py \
	-i "${TMPDIR}/${SAMPLEID}_gridss_unused_breakend_records_with_key_sorted.txt" \
	-o "${TMPDIR}/${SAMPLEID}_gridss_circular_inversions_too.vcf" \
	-circular CIRCULAR \
	-u "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_circular_inversions.vcf" \
	-gap "${GAP}"

# identify new VCF headers added by inversion program
# the circular_inversions output has same VCF headers as inversions output plus a few more VCF headers
grep '^#' "${TMPDIR}/${SAMPLEID}_gridss_circular_inversions_too.vcf" > "${TMPDIR}/${SAMPLEID}_gridss_circular_inversions_too_hdr.vcf"
sort "${TMPDIR}/${SAMPLEID}_gridss_circular_inversions_too_hdr.vcf" > "${TMPDIR}/${SAMPLEID}_gridss_circular_inversions_too_hdr_sorted.txt"
comm -23 "${TMPDIR}/${SAMPLEID}_gridss_circular_inversions_too_hdr_sorted.txt" "${TMPDIR}/${SAMPLEID}_gridss_hdr_sorted.txt" > "${TMPDIR}/${SAMPLEID}_gridss_inversions_new_hdrs_only.txt"

# join INS_DEL_INDEL_DUP, inversions and circular-inversions output, and also include remaining unsed BND records, sort it all by CHROM+POS, add headers ready for input to genomic-region-annotation
grep -v '^#' "${TMPDIR}/${SAMPLEID}_gridss_inversions_only.vcf" > "${TMPDIR}/${SAMPLEID}_gridss_inversions_only_nohdr.txt"
grep -v '^#' "${TMPDIR}/${SAMPLEID}_gridss_circular_inversions_too.vcf" > "${TMPDIR}/${SAMPLEID}_gridss_circular_inversions_too_nohdr.txt"
grep -v '^#' "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_circular_inversions.vcf" > "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_circular_inversions_nohdr.txt"
cat "${TMPDIR}/${SAMPLEID}_gridss_INS_DEL_INDEL_DUP_nohdr_sorted.txt" "${TMPDIR}/${SAMPLEID}_gridss_inversions_only_nohdr.txt" "${TMPDIR}/${SAMPLEID}_gridss_circular_inversions_too_nohdr.txt" "${TMPDIR}/${SAMPLEID}_gridss_BND_not_used_for_circular_inversions_nohdr.txt" \
	| sort -k1,1 -k2,2n > "${TMPDIR}/${SAMPLEID}_gridss_INS_DEL_INDEL_DUP_INV_BND_nohdr_sorted.txt"
grep '^#CHROM' "${TMPDIR}/${SAMPLEID}_gridss_hdr.vcf" > "${TMPDIR}/${SAMPLEID}_gridss_hdr_CHROM_only.txt"
grep -v '^#CHROM' "${TMPDIR}/${SAMPLEID}_gridss_hdr.vcf" > "${TMPDIR}/${SAMPLEID}_gridss_hdr_no_CHROM.txt"
cat "${TMPDIR}/${SAMPLEID}_gridss_hdr_no_CHROM.txt" "${TMPDIR}/${SAMPLEID}_gridss_INS_DEL_INDEL_DUP_new_hdrs_only.txt" "${TMPDIR}/${SAMPLEID}_gridss_inversions_new_hdrs_only.txt" "${TMPDIR}/${SAMPLEID}_gridss_hdr_CHROM_only.txt" \
	"${TMPDIR}/${SAMPLEID}_gridss_INS_DEL_INDEL_DUP_INV_BND_nohdr_sorted.txt" > "${TMPDIR}/${SAMPLEID}_gridss_INS_DEL_INDEL_DUP_INV_BND_got_GT_of_dot.vcf"

python "${SVannotation_software_directory}"/programs/convert_sample_GT_in_VCF_file.py -i "${TMPDIR}/${SAMPLEID}_gridss_INS_DEL_INDEL_DUP_INV_BND_got_GT_of_dot.vcf" -o "${output_SV_vcf_file}" -a ALL

