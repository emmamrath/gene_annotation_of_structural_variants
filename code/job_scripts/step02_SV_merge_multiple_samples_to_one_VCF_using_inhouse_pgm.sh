#!/bin/bash

# input variables
cohort_id=$1
INDIR=$2
OUTDIR=$3
TMPDIR='/nvme/emmrat/tmp_SV2'
REFERENCE='/nvme/emmrat/reference_genomes/hs37d5x/hs37d5x.fa'
sv_pgms_dir='/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation/programs'
window_size=100

mkdir -p "${OUTDIR}"
mkdir -p "${TMPDIR}"

echo 'Add key to SV files for all samples, and concatenate ready for sorting'
:>"${TMPDIR}"/"${cohort_id}"_sample_list.txt
:>"${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated.txt
one_sample_vcf_file=''

#for infile in "${INDIR}"/*_gridss_INS_DEL_INDEL_DUP_INV_BND.vcf; do
for infile in "${INDIR}"/*.vcf; do

    infile_basename=$(basename "${infile}")
    #sampleid=${infile:0:5}
    IFS='.' read -ra infile_array <<< "$infile_basename"
    sampleid="${infile_array[0]}"
    #sampleid=${infile/_gridss_INS_DEL_INDEL_DUP_INV_BND\.vcf/}

    echo '    ' "${sampleid}"
    cat "${infile}" | python "${sv_pgms_dir}"/add_sample_key_to_VCF_to_merge_multiple_samples.py "${sampleid}" > "${TMPDIR}"/"${cohort_id}"_"${sampleid}"_gridss_SV_with_sample_key.txt

    echo "${sampleid}" >> "${TMPDIR}"/"${cohort_id}"_sample_list.txt

    cat "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated.txt "${TMPDIR}"/"${cohort_id}"_"${sampleid}"_gridss_SV_with_sample_key.txt > "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_new.txt
    mv "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_new.txt "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated.txt

    one_sample_vcf_file="${infile}"
done
grep '^#' "${one_sample_vcf_file}" | grep -v '^#CHROM' > "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf

echo 'Sort the concatenated SV file ready to merge into one VCF'
cat "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated.txt | sort -k2,2 -k3,3n -k6,6 > "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted.txt

echo 'Split the concatenated SV file into INS, DEL, INDEL, DUP, INV and BND files before merging them into one VCF per SV type'
grep 'SVTYPE=INS' "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted.txt > "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_INS_only.txt
grep 'SVTYPE=DEL' "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted.txt > "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_DEL_only.txt
grep 'SVTYPE=INDEL' "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted.txt > "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_INDEL_only.txt
grep 'SVTYPE=DUP' "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted.txt > "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_DUP_only.txt
grep 'SVTYPE=INV' "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted.txt > "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_INV_only.txt
grep 'SVTYPE=BND' "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted.txt > "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_BND_only.txt

echo 'Merge the SV files so that all samples are in one VCF file for INS only'
cat "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_INS_only.txt | python "${sv_pgms_dir}"/merge_one_sample_VCFs_to_one_multisample_VCF.py "${TMPDIR}"/"${cohort_id}"_sample_list.txt \
    python "${sv_pgms_dir}"/merge_a_multi_sample_SV_VCF_for_similar_position_variants.py "${window_size}" "${REFERENCE}" > "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_INS_only.vcf
echo 'Merge the SV files so that all samples are in one VCF file for DEL only'
cat "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_DEL_only.txt | python "${sv_pgms_dir}"/merge_one_sample_VCFs_to_one_multisample_VCF.py "${TMPDIR}"/"${cohort_id}"_sample_list.txt \
    python "${sv_pgms_dir}"/merge_a_multi_sample_SV_VCF_for_similar_position_variants.py "${window_size}" "${REFERENCE}" > "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_DEL_only.vcf
echo 'Merge the SV files so that all samples are in one VCF file for INDEL only'
cat "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_INDEL_only.txt | python "${sv_pgms_dir}"/merge_one_sample_VCFs_to_one_multisample_VCF.py "${TMPDIR}"/"${cohort_id}"_sample_list.txt \
    python "${sv_pgms_dir}"/merge_a_multi_sample_SV_VCF_for_similar_position_variants.py "${window_size}" "${REFERENCE}" > "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_INDEL_only.vcf
echo 'Merge the SV files so that all samples are in one VCF file for DUP only'
cat "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_DUP_only.txt | python "${sv_pgms_dir}"/merge_one_sample_VCFs_to_one_multisample_VCF.py "${TMPDIR}"/"${cohort_id}"_sample_list.txt \
    python "${sv_pgms_dir}"/merge_a_multi_sample_SV_VCF_for_similar_position_variants.py "${window_size}" "${REFERENCE}" > "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_DUP_only.vcf
echo 'Merge the SV files so that all samples are in one VCF file for INV only'
cat "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_INV_only.txt | python "${sv_pgms_dir}"/merge_one_sample_VCFs_to_one_multisample_VCF.py "${TMPDIR}"/"${cohort_id}"_sample_list.txt \
    python "${sv_pgms_dir}"/merge_a_multi_sample_SV_VCF_for_similar_position_variants.py "${window_size}" "${REFERENCE}" > "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_INV_only.vcf
echo 'Merge the SV files so that all samples are in one VCF file for BND only'
cat "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_BND_only.txt | python "${sv_pgms_dir}"/merge_one_sample_VCFs_to_one_multisample_VCF.py "${TMPDIR}"/"${cohort_id}"_sample_list.txt \
    python "${sv_pgms_dir}"/merge_a_multi_sample_SV_VCF_for_similar_position_variants.py "${window_size}" "${REFERENCE}" > "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_BND_only.vcf

echo 'Put the header back on the merged SV VCF files, one file each for INS, DEL, INDEL, DUP, INV and BND'
cat "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_INS_only.vcf > "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_INS_only_allbp.vcf
cat "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_DEL_only.vcf > "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_DEL_only_allbp.vcf
cat "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_INDEL_only.vcf > "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_INDEL_only_allbp.vcf
cat "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_DUP_only.vcf > "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_DUP_only_allbp.vcf
cat "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_INV_only.vcf > "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_INV_only_allbp.vcf
cat "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_BND_only.vcf > "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_BND_only_allbp.vcf

echo 'filter out variants smaller than 50 bp in length with filter_out_small_variants.py -f 49'
cat "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_INS_only_allbp.vcf | python "${sv_pgms_dir}"/filter_out_small_variants.py -f 49 > "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_INS_only.vcf
cat "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_DEL_only_allbp.vcf | python "${sv_pgms_dir}"/filter_out_small_variants.py -f 49 > "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_DEL_only.vcf
cat "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_INDEL_only_allbp.vcf | python "${sv_pgms_dir}"/filter_out_small_variants.py -f 49 > "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_INDEL_only.vcf
cat "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_DUP_only_allbp.vcf | python "${sv_pgms_dir}"/filter_out_small_variants.py -f 49 > "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_DUP_only.vcf
cat "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_INV_only_allbp.vcf | python "${sv_pgms_dir}"/filter_out_small_variants.py -f 49 > "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_INV_only.vcf
cat "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_BND_only_allbp.vcf | python "${sv_pgms_dir}"/filter_out_small_variants.py -f 49 > "${OUTDIR}"/"${cohort_id}"_gridss_AllSamples_BND_only.vcf

rm "${TMPDIR}"/"${cohort_id}"_sample_list.txt
rm "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated.txt
rm "${TMPDIR}"/"${cohort_id}"_"${sampleid}"_gridss_SV_with_sample_key.txt
rm "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf
rm "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted.txt
rm "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_INS_only.txt
rm "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_DEL_only.txt
rm "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_INDEL_only.txt
rm "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_DUP_only.txt
rm "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_INV_only.txt
rm "${TMPDIR}"/"${cohort_id}"_gridss_SV_concatenated_sorted_BND_only.txt
rm "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_INS_only.vcf
rm "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_DEL_only.vcf
rm "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_INDEL_only.vcf
rm "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_DUP_only.vcf
rm "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_INV_only.vcf
rm "${TMPDIR}"/"${cohort_id}"_temp_hdr.vcf "${TMPDIR}"/"${cohort_id}"_temp_SV_nohdr_BND_only.vcf

echo 'Finished'


