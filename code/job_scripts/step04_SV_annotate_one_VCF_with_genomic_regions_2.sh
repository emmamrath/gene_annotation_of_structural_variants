#!/bin/bash
set -euo pipefail

unique_id=$1
input_VCF_file=$2
output_VCF_file=$3

TMPDIR='/nvme/emmrat/tmp_SV4'
SVannotation_software_directory='/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation'

input_VCF_file_with_ID="${TMPDIR}"/"${unique_id}"_SV_annotate_genomic_regions_2_with_ID.vcf
input_VCF_file_with_ID_sorted_by_ID="${TMPDIR}"/"${unique_id}"_SV_annotate_genomic_regions_2_with_ID_sorted_by_ID.vcf
annotated_VCF_file="${TMPDIR}"/"${unique_id}"_SV_annotate_genomic_regions_2.vcf
annotated_VCF_file_with_consequences="${TMPDIR}"/"${unique_id}"_SV_annotate_genomic_regions_2_and_consequences.vcf

##################################################
# Annotate the VCF with gene regions
##################################################
echo 'Annotate the variants with genomic regions'

echo 'Remove SV variants having one or both breakends fall in a region of low-complexity or simple-repeats. Making sure the VCF file is sorted. Annotation by bedtools intersect relies on this.'
# They are probably mapped incorrectly and so the SV is false.
grep '^#' "${input_VCF_file}" > "${TMPDIR}/${unique_id}_VCF_hdr.vcf"
grep -v '^#' "${input_VCF_file}" | sort -k1,1 -k2,2n > "${TMPDIR}/${unique_id}_VCF_body.vcf"

cat "${TMPDIR}/${unique_id}_VCF_body.vcf" | \
	python "${SVannotation_software_directory}"/programs/add_breakend_key_to_vcf_record.py -p LEFT > "${TMPDIR}/${unique_id}_add_key_LEFT.txt"
bedtools intersect -a "${TMPDIR}/${unique_id}_add_key_LEFT.txt" \
	-b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_Repeats_20171007_LowComplexity_SimpleRepeats.bed" -v | cut -d$'\t' -f4- > "${TMPDIR}/${unique_id}_add_key_LEFT_removed.vcf"
cat "${TMPDIR}/${unique_id}_add_key_LEFT_removed.vcf" | \
	python "${SVannotation_software_directory}"/programs/add_breakend_key_to_vcf_record.py -p RIGHT > "${TMPDIR}/${unique_id}_add_key_RIGHT.txt"
bedtools intersect -a "${TMPDIR}/${unique_id}_add_key_RIGHT.txt" \
	-b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_Repeats_20171007_LowComplexity_SimpleRepeats.bed" -v | cut -d$'\t' -f4- > "${TMPDIR}/${unique_id}_removed_low_complexity_body.vcf"
cat "${TMPDIR}/${unique_id}_VCF_hdr.vcf" "${TMPDIR}/${unique_id}_removed_low_complexity_body.vcf" > "${TMPDIR}/${unique_id}_removed_low_complexity.vcf"
rm "${TMPDIR}/${unique_id}_VCF_body.vcf"
rm "${TMPDIR}/${unique_id}_add_key_LEFT.txt" "${TMPDIR}/${unique_id}_add_key_LEFT_removed.vcf"
rm "${TMPDIR}/${unique_id}_add_key_RIGHT.txt" "${TMPDIR}/${unique_id}_removed_low_complexity_body.vcf"

echo 'Making sure the VCF file is sorted. Annotation by bedtools intersect relies on this.'
grep -v '^#' "${TMPDIR}/${unique_id}_removed_low_complexity.vcf" > "${TMPDIR}/${unique_id}_VCF_body2.vcf"
sort -k1,1 -k2,2n "${TMPDIR}/${unique_id}_VCF_body2.vcf" > "${TMPDIR}/${unique_id}_VCF_body_sorted.vcf"
cat "${TMPDIR}/${unique_id}_VCF_hdr.vcf" "${TMPDIR}/${unique_id}_VCF_body_sorted.vcf" > "${TMPDIR}/${unique_id}_VCF_sorted.vcf"
rm "${TMPDIR}/${unique_id}_VCF_body2.vcf"

echo 'generate_VCF_IDs_for_VCF_variants.py'
cat "${TMPDIR}/${unique_id}_VCF_sorted.vcf" | \
	python "${SVannotation_software_directory}"/programs/generate_VCF_IDs_for_VCF_variants.py -p SV_ \
	> "${input_VCF_file_with_ID}"

echo 'For annotation purposes only, remove SV variants larger than 1,000,000 bp. They will be included in the final file but will not be annotated.'
cat "${input_VCF_file_with_ID}" | \
	python "${SVannotation_software_directory}"/programs/filter_out_large_variants.py -f 1000001 \
	> "${TMPDIR}/"${unique_id}"_SV_merged100bp_with_ID_without_large_variants.vcf"

# Convert higher-level VCF variants (but not remaing BND records) to BED format of intervals so that easy BEDTOOLS intersection can be done
echo 'convert_VCF_structural_variants_to_BED_format.py'
cat "${TMPDIR}/"${unique_id}"_SV_merged100bp_with_ID_without_large_variants.vcf" | \
	python "${SVannotation_software_directory}"/programs/convert_VCF_structural_variants_to_BED_format.py > "${TMPDIR}/${unique_id}_VCF.bed"

# Sort by ID so that VCF and the gene bed resulting from bedtools intersect will be in the same order.
# bedtools intersect relies on the same order, and they are ordered and matched by chrom+pos by bedtools intersect.
# More than one variant may have the same chrom+pos (because INS/DEL/etc. are kept in separate variant records).
# The program convert_VCF_structural_variants_to_BED_format.py will gives these different IDs.
# The annotate_VCF_structural_variants_with_bedtools_intersect_genes.py run after bedtools intersect also relies on having the same order in the 2 files.
# Thus, order by ID, so that 2 variants having the same chrom+pos will have the same sort order.
echo 'sort VCF file by ID'
grep '^#' "${input_VCF_file_with_ID}" > "${TMPDIR}/${unique_id}_VCF_with_ID_hdr.vcf"
grep -v '^#' "${input_VCF_file_with_ID}" | sort -k3,3 > "${TMPDIR}/${unique_id}_VCF_with_ID_body.vcf"
cat "${TMPDIR}/${unique_id}_VCF_with_ID_hdr.vcf" "${TMPDIR}/${unique_id}_VCF_with_ID_body.vcf" > "${input_VCF_file_with_ID_sorted_by_ID}"

########## The various annotations ##########

# Run bedtools to identify genes hits by VCF's BED intervals, which compares by chromosomal position (comparing gene interval to VCF's BED interval, output VCF.ID for matching intervals)
echo 'bedtools intersect for UPSTREAM1TO20KB of UCSC Genes'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_UPSTREAM1TO20KB_RefSeq_20170928.txt" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_UPSTREAM1TO20KB.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_UPSTREAM1TO20KB.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_UPSTREAM1TO20KB_sorted.txt"
cat "${input_VCF_file_with_ID_sorted_by_ID}" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_UPSTREAM1TO20KB_sorted.txt" -t UPSTREAM1TO20KB -d 'Intersecting gene' > "${TMPDIR}/${unique_id}_after_UPSTREAM1TO20KB.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_UPSTREAM1TO20KB.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_UPSTREAM1TO20KB_sorted.txt"

echo 'bedtools intersect for UPSTREAM20KBTO50KB of UCSC Genes'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_UPSTREAM20KBTO50KB_RefSeq_20170928.txt" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_UPSTREAM20KBTO50KB.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_UPSTREAM20KBTO50KB.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_UPSTREAM20KBTO50KB_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_UPSTREAM1TO20KB.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_UPSTREAM20KBTO50KB_sorted.txt" -t UPSTREAM20KBTO50KB -d 'Intersecting gene exon' > "${TMPDIR}/${unique_id}_after_UPSTREAM20KBTO50KB.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_UPSTREAM20KBTO50KB.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_UPSTREAM20KBTO50KB_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_UPSTREAM1TO20KB.txt"

echo 'bedtools intersect for UPSTREAM50KBTO100KB of UCSC Genes'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_UPSTREAM50KBTO100KB_RefSeq_20170928.txt" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_UPSTREAM50KBTO100KB.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_UPSTREAM50KBTO100KB.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_UPSTREAM50KBTO100KB_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_UPSTREAM20KBTO50KB.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_UPSTREAM50KBTO100KB_sorted.txt" -t UPSTREAM50KBTO100KB -d 'Intersecting gene intron' > "${TMPDIR}/${unique_id}_after_UPSTREAM50KBTO100KB.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_UPSTREAM50KBTO100KB.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_UPSTREAM50KBTO100KB_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_UPSTREAM20KBTO50KB.txt"

echo 'bedtools intersect for DOWNSTREAM1TO20KB of UCSC Genes'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_DOWNSTREAM1TO20KB_RefSeq_20170928.txt" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_DOWNSTREAM1TO20KB.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_DOWNSTREAM1TO20KB.txt" > "${TMPDIR}/${unique_id}_VCF_DOWNSTREAM1TO20KB_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_UPSTREAM50KBTO100KB.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_DOWNSTREAM1TO20KB_sorted.txt" -t DOWNSTREAM1TO20KB -d 'Intersecting CDS region' > "${TMPDIR}/${unique_id}_after_DOWNSTREAM1TO20KB.txt"
rm "${TMPDIR}/${unique_id}_VCF_DOWNSTREAM1TO20KB.txt"
rm "${TMPDIR}/${unique_id}_VCF_DOWNSTREAM1TO20KB_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_UPSTREAM50KBTO100KB.txt"

echo 'bedtools intersect for DOWNSTREAM20KBTO50KB of UCSC Genes'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_DOWNSTREAM20KBTO50KB_RefSeq_20170928.txt" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_DOWNSTREAM20KBTO50KB.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_DOWNSTREAM20KBTO50KB.txt" > "${TMPDIR}/${unique_id}_VCF_DOWNSTREAM20KBTO50KB_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_DOWNSTREAM1TO20KB.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_DOWNSTREAM20KBTO50KB_sorted.txt" -t DOWNSTREAM20KBTO50KB -d 'Intersecting CDS exons' > "${TMPDIR}/${unique_id}_after_DOWNSTREAM20KBTO50KB.txt"
rm "${TMPDIR}/${unique_id}_VCF_DOWNSTREAM20KBTO50KB.txt"
rm "${TMPDIR}/${unique_id}_VCF_DOWNSTREAM20KBTO50KB_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_DOWNSTREAM1TO20KB.txt"

echo 'bedtools intersect for DOWNSTREAM50KBTO100KB of UCSC Genes'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_DOWNSTREAM50KBTO100KB_RefSeq_20170928.txt" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_DOWNSTREAM20KBTO50KB.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5_sorted.txt" -t DOWNSTREAM50KBTO100KB -d "Intersecting UTR-5' region" > "${TMPDIR}/${unique_id}_after_DOWNSTREAM50KBTO100KB.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_DOWNSTREAM20KBTO50KB.txt"

echo 'bedtools intersect for SURROUND3000BEFORE0AFTER of UCSC Genes'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_SURROUND3000BEFORE0AFTER_RefSeq_20170928.txt" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_SURROUND3000BEFORE0AFTER.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_SURROUND3000BEFORE0AFTER.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_SURROUND3000BEFORE0AFTER_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_DOWNSTREAM50KBTO100KB.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_SURROUND3000BEFORE0AFTER_sorted.txt" -t SURROUND3000BEFORE0AFTER -d "Intersecting UTR-3' region" > "${TMPDIR}/${unique_id}_after_SURROUND3000BEFORE0AFTER.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_SURROUND3000BEFORE0AFTER.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_SURROUND3000BEFORE0AFTER_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_DOWNSTREAM50KBTO100KB.txt"

echo 'sort by chrom + pos (instead of by ID)'
grep '^#' "${TMPDIR}/${unique_id}_after_SURROUND3000BEFORE0AFTER.txt" > "${TMPDIR}/${unique_id}_annotated_hdr.txt"
grep -v '^#' "${TMPDIR}/${unique_id}_after_SURROUND3000BEFORE0AFTER.txt" | sort -k1,1 -k2,2n > "${TMPDIR}/${unique_id}_annotated_and_sorted_body.txt"
cat "${TMPDIR}/${unique_id}_annotated_hdr.txt" "${TMPDIR}/${unique_id}_annotated_and_sorted_body.txt" > "${annotated_VCF_file}"

##### Add probable consequences of the structural variants

#echo 'add probable consequences of the structural variants'
#python "${SVannotation_software_directory}"/programs/annotate_vcf_with_structural_variant_exon_consequence.py -i "${annotated_VCF_file}" -o "${annotated_VCF_file_with_consequences}" -g "${SVannotation_software_directory}/reference_data/UCSC_Genes_canonical_exons_20180613.bed"

#####

cp "${annotated_VCF_file}" "${output_VCF_file}"
rm "${annotated_VCF_file}"

echo 'Finished!'


