#!/bin/bash
set -euo pipefail

cohortid=$1
sampleid=$2
input_SV_file=$3
input_CNV_file=$4
final_output_VCF_file=$5
final_output_TSV_file=$6
final_output_TSV_overlap_stats_file=$7

uniqueid="${cohortid}"_"${sampleid}"
SVCODEDIR='/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation/programs'
TMPDIR='/nvme/emmrat/tmp_SV'

# Convert VCF files to TSV files for input to R
grep -v 'SVTYPE=BND' "${input_SV_file}" > "${TMPDIR}"/"${uniqueid}".SV_no_BND.vcf
python "${SVCODEDIR}"/convert_VCF_INFO_one_sample_to_tab_delimited_and_intact_fields.py -i "${TMPDIR}"/"${uniqueid}".SV_no_BND.vcf -o "${TMPDIR}"/"${uniqueid}".SV_no_BND.tsv
python "${SVCODEDIR}"/convert_VCF_INFO_one_sample_to_tab_delimited_and_intact_fields.py -i "${input_CNV_file}" -o "${TMPDIR}"/"${uniqueid}".CNV.tsv

# Call R to merge SV and BICSEQ, output VCF fields in TSV format, and output final TSV file.
# Rscript merge_bicseq_cnv_into_annotated_gridss_sv.R ..
echo 'Rscript' "${SVCODEDIR}"/merge_bicseq_cnv_into_gridss_sv.R "${TMPDIR}"/"${uniqueid}".SV_no_BND.tsv "${TMPDIR}"/"${uniqueid}".CNV.tsv "${TMPDIR}"/"${uniqueid}".SV_with_CNV.vcf "${final_output_TSV_file}" "${final_output_TSV_overlap_stats_file}"
Rscript "${SVCODEDIR}"/merge_bicseq_cnv_into_gridss_sv.R "${TMPDIR}"/"${uniqueid}".SV_no_BND.tsv "${TMPDIR}"/"${uniqueid}".CNV.tsv "${TMPDIR}"/"${uniqueid}".SV_with_CNV.vcf "${final_output_TSV_file}" "${final_output_TSV_overlap_stats_file}"

# Prepare files for concatenating the results into one VCF file
grep '^##' "${input_SV_file}" > "${TMPDIR}"/"${uniqueid}".hdr1.txt
grep '^#CHROM' "${input_SV_file}" > "${TMPDIR}"/"${uniqueid}".hdr2.txt
:>"${TMPDIR}"/"${uniqueid}".hdr_new.txt
echo '##INFO=<ID=BICSEQ_POS,Number=1,Type=Integer,Description="The VCF.POS of overlapping bicseq copy number variant">' >> "${TMPDIR}"/"${uniqueid}".hdr_new.txt
echo '##INFO=<ID=BICSEQ_END,Number=1,Type=Integer,Description="The VCF.INFO.END of overlapping bicseq copy number variant">' >> "${TMPDIR}"/"${uniqueid}".hdr_new.txt
echo '##INFO=<ID=BICSEQ_OOE,Number=1,Type=Integer,Description="The VCF.INFO.OOE of overlapping bicseq copy number variant">' >> "${TMPDIR}"/"${uniqueid}".hdr_new.txt
echo '##INFO=<ID=BICSEQ_SVLEN,Number=1,Type=Integer,Description="The VCF.INFO.SVLEN of overlapping bicseq copy number variant">' >> "${TMPDIR}"/"${uniqueid}".hdr_new.txt
echo '##INFO=<ID=BICSEQ_LOG2COPYRATIO,Number=1,Type=Integer,Description="The VCF.INFO.LOG2COPYRATIO of overlapping bicseq copy number variant">' >> "${TMPDIR}"/"${uniqueid}".hdr_new.txt
echo '##INFO=<ID=BICSEQ_PVALUE,Number=1,Type=Integer,Description="The VCF.INFO.PVALUE of overlapping bicseq copy number variant">' >> "${TMPDIR}"/"${uniqueid}".hdr_new.txt
echo '##INFO=<ID=BICSEQ_OBSERVED,Number=1,Type=Integer,Description="The VCF.INFO.OBSERVED of overlapping bicseq copy number variant">' >> "${TMPDIR}"/"${uniqueid}".hdr_new.txt
echo '##INFO=<ID=BICSEQ_EXPECTED,Number=1,Type=Integer,Description="The VCF.INFO.EXPECTED of overlapping bicseq copy number variant">' >> "${TMPDIR}"/"${uniqueid}".hdr_new.txt
grep 'SVTYPE=BND' "${input_SV_file}" > "${TMPDIR}"/"${uniqueid}".SV_only_BND.vcf

# Concatenate the results into one VCF file
cat "${TMPDIR}"/"${uniqueid}".SV_with_CNV.vcf "${TMPDIR}"/"${uniqueid}".SV_only_BND.vcf | sort -t$'\t' -k1,1 -k2,2n > "${TMPDIR}"/"${uniqueid}".sorted_results.txt
cat "${TMPDIR}"/"${uniqueid}".hdr1.txt "${TMPDIR}"/"${uniqueid}".hdr_new.txt "${TMPDIR}"/"${uniqueid}".hdr2.txt "${TMPDIR}"/"${uniqueid}".sorted_results.txt > "${final_output_VCF_file}"

# Clean up temp files
rm "${TMPDIR}"/"${uniqueid}".hdr1.txt
rm "${TMPDIR}"/"${uniqueid}".hdr2.txt
rm "${TMPDIR}"/"${uniqueid}".hdr_new.txt
rm "${TMPDIR}"/"${uniqueid}".SV_no_BND.vcf
rm "${TMPDIR}"/"${uniqueid}".SV_no_BND.tsv
rm "${TMPDIR}"/"${uniqueid}".CNV.tsv
rm "${TMPDIR}"/"${uniqueid}".SV_only_BND.vcf
rm "${TMPDIR}"/"${uniqueid}".SV_with_CNV.vcf
rm "${TMPDIR}"/"${uniqueid}".sorted_results.txt


