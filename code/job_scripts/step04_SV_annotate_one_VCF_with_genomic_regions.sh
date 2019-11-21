#!/bin/bash
set -euo pipefail

unique_id=$1
input_VCF_file=$2
output_VCF_file=$3

TMPDIR='/nvme/emmrat/tmp_SV2'
SVannotation_software_directory='/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation'

input_VCF_file_with_ID="${TMPDIR}"/"${unique_id}"_SV_merged100bp_with_ID.vcf
input_VCF_file_with_ID_sorted_by_ID="${TMPDIR}"/"${unique_id}"_SV_merged100bp_with_ID_sorted_by_ID.vcf
annotated_VCF_file="${TMPDIR}"/"${unique_id}"_SV_INS_DEL_INDEL_DUP_INV_BND_merged100bp_annotated_with_genomic_regions_and_1000G.vcf
annotated_TSV_file="${TMPDIR}"/"${unique_id}"_SV_INS_DEL_INDEL_DUP_INV_BND_merged100bp_annotated_with_genomic_regions_and_1000G.tsv
annotated_VCF_file_with_consequences="${TMPDIR}"/"${unique_id}"_SV_INS_DEL_INDEL_DUP_INV_BND_merged100bp_annotated_with_genomic_regions_and_1000G_and_consequences.vcf

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
echo 'bedtools intersect for UCSC GeneAndGenePredictions'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_genes_RefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_sorted.txt"
cat "${input_VCF_file_with_ID_sorted_by_ID}" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_sorted.txt" -t HITGENES -d 'Intersecting gene' > "${TMPDIR}/${unique_id}_after_HITGENES.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_sorted.txt"

echo 'bedtools intersect for UCSC GeneAndGenePredictions exons'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_exons_RefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_exons.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_exons.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_exons_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_HITGENES.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_exons_sorted.txt" -t HITEXONS -d 'Intersecting gene exon' > "${TMPDIR}/${unique_id}_after_HITEXONS.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_exons.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_exons_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_HITGENES.txt"

echo 'bedtools intersect for UCSC GeneAndGenePredictions introns'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_introns_RefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_introns.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_introns.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_introns_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_HITEXONS.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_introns_sorted.txt" -t HITINTRONS -d 'Intersecting gene intron' > "${TMPDIR}/${unique_id}_after_HITINTRONS.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_introns.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_introns_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_HITEXONS.txt"

echo 'bedtools intersect for UCSC GeneAndGenePredictions CDS'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_CDS_RefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_CDS.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_CDS.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_CDS_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_HITINTRONS.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_CDS_sorted.txt" -t HITCDS -d 'Intersecting CDS region' > "${TMPDIR}/${unique_id}_after_HITCDS.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_CDS.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_CDS_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_HITINTRONS.txt"

echo 'bedtools intersect for UCSC GeneAndGenePredictions CDS exons'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_CDSexons_RefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_CDS_exons.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_CDS_exons.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_CDS_exons_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_HITCDS.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_CDS_exons_sorted.txt" -t HITCDSEXONS -d 'Intersecting CDS exons' > "${TMPDIR}/${unique_id}_after_HITCDSEXONS.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_CDS_exons.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_CDS_exons_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_HITCDS.txt"

echo 'bedtools intersect for UCSC GenesAndGenePredictions UTR5'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_UTR5_RefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_HITCDSEXONS.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5_sorted.txt" -t HITUTR5 -d "Intersecting UTR-5' region" > "${TMPDIR}/${unique_id}_after_HITUTR5.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_HITCDSEXONS.txt"

echo 'bedtools intersect for UCSC GenesAndGenePredictions UTR3'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_UTR3_RefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR3.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR3.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR3_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_HITUTR5.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR3_sorted.txt" -t HITUTR3 -d "Intersecting UTR-3' region" > "${TMPDIR}/${unique_id}_after_HITUTR3.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR3.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR3_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_HITUTR5.txt"

echo 'bedtools intersect for UCSC TranscriptionStartSite + 1500bp downstream'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_TranscriptionStartSite_Plus1500bp_RefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_TSS_Plus1500bp.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_TSS_Plus1500bp.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_TSS_Plus1500bp_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_HITUTR3.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_TSS_Plus1500bp_sorted.txt" -t HITTSS -d 'Intersecting transcription start site and 1500 bp downstream' > "${TMPDIR}/${unique_id}_after_HITTSS.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_TSS_Plus1500bp.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_TSS_Plus1500bp_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_HITUTR3.txt"

echo 'bedtools intersect for UCSC GeneAndGenePredictions flank 1000bp'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_flank1000bp_RefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank1000bp.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank1000bp.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank1000bp_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_HITTSS.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank1000bp_sorted.txt" -t FLANK1000GENES -d 'Intersecting gene flank of 1000 bp (on either of but not including the gene)' > "${TMPDIR}/${unique_id}_after_FLANK1000GENES.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank1000bp.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank1000bp_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_HITTSS.txt"

echo 'bedtools intersect for UCSC GeneAndGenePredictions flank 3000bp'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_flank3000bp_RefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank3000bp.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank3000bp.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank3000bp_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_FLANK1000GENES.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank3000bp_sorted.txt" -t FLANK3000GENES -d 'Intersecting gene flank of 3000 bp (on either of but not including the gene)' > "${TMPDIR}/${unique_id}_after_FLANK3000GENES.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank3000bp.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank3000bp_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_FLANK1000GENES.txt"

echo 'bedtools intersect for UCSC GeneAndGenePredictions flank 10000bp'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_flank10000bp_RefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank10000bp.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank10000bp.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank10000bp_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_FLANK3000GENES.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank10000bp_sorted.txt" -t FLANK10000GENES -d 'Intersecting gene flank of 10000 bp (on either of but not including the gene)' > "${TMPDIR}/${unique_id}_after_FLANK10000GENES.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank10000bp.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank10000bp_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_FLANK3000GENES.txt"

echo 'bedtools intersect for UCSC EnsembleRegulatoryBuild'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_EnsembleRegulatoryBuild_20170526.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_regulatoryregions.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_regulatoryregions.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_regulatoryregions_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_FLANK10000GENES.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_regulatoryregions_sorted.txt" -t HITREGULATORY -d 'Intersecting regulatory region' > "${TMPDIR}/${unique_id}_after_HITREGULATORY.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_regulatoryregions.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_regulatoryregions_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_FLANK10000GENES.txt"

echo 'bedtools intersect for Mobster ExistingMobileElementRegions'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/Mobster_GRCh37_ExistingMobileElementRegions_20170526.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_mobileelements.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_mobileelements.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_mobileelements_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_HITREGULATORY.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_mobileelements_sorted.txt" -t HITMOBEL -d 'Intersecting existing mobile element region' > "${TMPDIR}/${unique_id}_after_HITMOBEL.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_mobileelements.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_mobileelements_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_HITREGULATORY.txt"

echo 'bedtools intersect for UCSC Regulation CpGIslands'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_Regulation_CpGIslands_20170614.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_CpGIslands.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_CpGIslands.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_CpGIslands_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_HITMOBEL.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_CpGIslands_sorted.txt" -t HITCPG -d 'Intersecting CpG island' > "${TMPDIR}/${unique_id}_after_HITCPG.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_CpGIslands.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_CpGIslands_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_HITMOBEL.txt"

echo 'bedtools intersect for subtelomeric repeat region'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_SREboundaries_20170927_Subtelomeric_Repeat_Boundaries.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_subtelomeric.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_subtelomeric.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_subtelomeric_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_HITCPG.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_subtelomeric_sorted.txt" -t HITSRE -d 'Intersecting subtelomeric repeat region' > "${TMPDIR}/${unique_id}_after_HITSRE.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_subtelomeric.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_subtelomeric_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_HITCPG.txt"

echo 'bedtools intersect for Non-RefSeq GeneAndGenePredictions'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_genes_NonRefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_non_refseq.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_non_refseq.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_non_refseq_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_HITSRE.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_non_refseq_sorted.txt" -t NONREFSEQGENES -d 'Intersecting non-RefSeq gene' > "${TMPDIR}/${unique_id}_after_NONREFSEQGENES.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_non_refseq.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_non_refseq_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_HITSRE.txt"

echo 'bedtools intersect for Non-RefSeq GeneAndGenePredictions Exons'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_exons_NonRefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_exons_non_refseq.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_exons_non_refseq.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_exons_non_refseq_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_NONREFSEQGENES.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_exons_non_refseq_sorted.txt" -t NONREFSEQEXONS -d 'Intersecting non-RefSeq gene exon' > "${TMPDIR}/${unique_id}_after_NONREFSEQEXONS.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_exons_non_refseq.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_exons_non_refseq_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_NONREFSEQGENES.txt"

echo 'bedtools intersect for Non-RefSeq GeneAndGenePredictions Introns'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_introns_NonRefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_introns_non_refseq.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_introns_non_refseq.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_introns_non_refseq_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_NONREFSEQEXONS.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_introns_non_refseq_sorted.txt" -t NONREFSEQINTRONS -d 'Intersecting non-RefSeq gene intron' > "${TMPDIR}/${unique_id}_after_NONREFSEQINTRONS.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_introns_non_refseq.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_introns_non_refseq_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_NONREFSEQEXONS.txt"

echo 'bedtools intersect for Non-RefSeq GeneAndGenePredictions CDS'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_CDS_NonRefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_CDS_non_refseq.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_CDS_non_refseq.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_CDS_non_refseq_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_NONREFSEQINTRONS.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_CDS_non_refseq_sorted.txt" -t NONREFSEQCDS -d 'Intersecting non-RefSeq CDS region' > "${TMPDIR}/${unique_id}_after_NONREFSEQCDS.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_CDS_non_refseq.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_CDS_non_refseq_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_NONREFSEQINTRONS.txt"

echo 'bedtools intersect for Non-RefSeq GeneAndGenePredictions CDS exons'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_CDSexons_NonRefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_CDS_exons_non_refseq.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_CDS_exons_non_refseq.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_CDS_exons_non_refseq_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_NONREFSEQCDS.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_CDS_exons_non_refseq_sorted.txt" -t NONREFSEQCDSEXONS -d 'Intersecting non-RefSeq CDS exons' > "${TMPDIR}/${unique_id}_after_NONREFSEQCDSEXONS.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_CDS_exons_non_refseq.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_CDS_exons_non_refseq_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_NONREFSEQCDS.txt"

echo 'bedtools intersect for Non-RefSeq GenesAndGenePredictions UTR5'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_UTR5_NonRefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5_non_refseq.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5_non_refseq.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5_non_refseq_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_NONREFSEQCDSEXONS.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5_non_refseq_sorted.txt" -t NONREFSEQUTR5 -d "Intersecting non-RefSeq UTR-5' region" > "${TMPDIR}/${unique_id}_after_NONREFSEQUTR5.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5_non_refseq.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR5_non_refseq_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_NONREFSEQCDSEXONS.txt"

echo 'bedtools intersect for Non-RefSeq GenesAndGenePredictions UTR3'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_UTR3_NonRefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR3_non_refseq.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR3_non_refseq.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR3_non_refseq_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_NONREFSEQUTR5.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR3_non_refseq_sorted.txt" -t NONREFSEQUTR3 -d "Intersecting non-RefSeq UTR-3' region" > "${TMPDIR}/${unique_id}_after_NONREFSEQUTR3.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR3_non_refseq.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_UTR3_non_refseq_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_NONREFSEQUTR5.txt"

echo 'bedtools intersect for Non-RefSeq TranscriptionStartSite + 1500bp downstream'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_TranscriptionStartSite_Plus1500bp_NonRefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_TSS_Plus1500bp_non_refseq.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_TSS_Plus1500bp_non_refseq.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_TSS_Plus1500bp_non_refseq_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_NONREFSEQUTR3.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_TSS_Plus1500bp_non_refseq_sorted.txt" -t NONREFSEQTSS -d 'Intersecting non-RefSeq transcription start site and 1500 bp downstream' > "${TMPDIR}/${unique_id}_after_NONREFSEQTSS.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_TSS_Plus1500bp_non_refseq.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_TSS_Plus1500bp_non_refseq_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_NONREFSEQUTR3.txt"

echo 'bedtools intersect for Non-RefSeq GeneAndGenePredictions flank 1000bp'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_flank1000bp_NonRefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank1000bp_non_refseq.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank1000bp_non_refseq.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank1000bp_non_refseq_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_NONREFSEQTSS.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank1000bp_non_refseq_sorted.txt" -t FLANK1000NONREFSEQ -d 'Intersecting non-RefSeq gene flank of 1000 bp (on either of but not including the gene)' > "${TMPDIR}/${unique_id}_after_FLANK1000NONREFSEQ.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank1000bp_non_refseq.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank1000bp_non_refseq_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_NONREFSEQTSS.txt"

echo 'bedtools intersect for Non-RefSeq GeneAndGenePredictions flank 3000bp'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_flank3000bp_NonRefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank3000bp_non_refseq.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank3000bp_non_refseq.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank3000bp_non_refseq_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_FLANK1000NONREFSEQ.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank3000bp_non_refseq_sorted.txt" -t FLANK3000NONREFSEQ -d 'Intersecting non-RefSeq gene flank of 3000 bp (on either of but not including the gene)' > "${TMPDIR}/${unique_id}_after_FLANK3000NONREFSEQ.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank3000bp_non_refseq.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank3000bp_non_refseq_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_FLANK1000NONREFSEQ.txt"

echo 'bedtools intersect for Non-RefSeq GeneAndGenePredictions flank 10000bp'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/UCSC_GRCh37_GenesAndGenePredictions_flank10000bp_NonRefSeq_20170928.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank10000bp_non_refseq.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank10000bp_non_refseq.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank10000bp_non_refseq_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_FLANK3000NONREFSEQ.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank10000bp_non_refseq_sorted.txt" -t FLANK10000NONREFSEQ -d 'Intersecting non-RefSeq gene flank of 10000 bp (on either of but not including the gene)' > "${TMPDIR}/${unique_id}_after_FLANK10000NONREFSEQ.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank10000bp_non_refseq.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_genes_flank10000bp_non_refseq_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_FLANK3000NONREFSEQ.txt"

echo 'bedtools intersect for 1000G variants SV regions'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_region.germline.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_1000G_SVregions.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_1000G_SVregions.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_1000G_SVregions_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_FLANK10000NONREFSEQ.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_1000G_SVregions_sorted.txt" -t SVREGION1000G -d '1000 Genomes Project Phase 3 SV analysis, Sudmant et al. 2015, germline variant regions (CNV,INS,INV).' > "${TMPDIR}/${unique_id}_after_SVREGION1000G.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_1000G_SVregions.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_1000G_SVregions_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_FLANK10000NONREFSEQ.txt"

echo 'bedtools intersect for 1000G variants SV calls'
bedtools intersect -a "${TMPDIR}/${unique_id}_VCF.bed" -b "${SVannotation_software_directory}/reference_data/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.bed" -wao | awk '$15 != "0"' > "${TMPDIR}/${unique_id}_VCF_intersect_1000G_SVcalls.txt"
sort -k4,4 "${TMPDIR}/${unique_id}_VCF_intersect_1000G_SVcalls.txt" > "${TMPDIR}/${unique_id}_VCF_intersect_1000G_SVcalls_sorted.txt"
cat "${TMPDIR}/${unique_id}_after_SVREGION1000G.txt" \
	| python "${SVannotation_software_directory}"/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py -i "${TMPDIR}/${unique_id}_VCF_intersect_1000G_SVcalls_sorted.txt" -t SVCALL1000G -d '1000 Genomes Project Phase 3 SV analysis, Sudmant et al. 2015, germline variant calls (CNV,DEL,DUP,INS,INV).' > "${TMPDIR}/${unique_id}_after_SVCALL1000G.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_1000G_SVcalls.txt"
rm "${TMPDIR}/${unique_id}_VCF_intersect_1000G_SVcalls_sorted.txt"
rm "${TMPDIR}/${unique_id}_after_SVREGION1000G.txt"

echo 'sort by chrom + pos (instead of by ID)'
grep '^#' "${TMPDIR}/${unique_id}_after_SVCALL1000G.txt" > "${TMPDIR}/${unique_id}_annotated_hdr.txt"
grep -v '^#' "${TMPDIR}/${unique_id}_after_SVCALL1000G.txt" | sort -k1,1 -k2,2n > "${TMPDIR}/${unique_id}_annotated_and_sorted_body.txt"
cat "${TMPDIR}/${unique_id}_annotated_hdr.txt" "${TMPDIR}/${unique_id}_annotated_and_sorted_body.txt" > "${annotated_VCF_file}"

##### Add probable consequences of the structural variants

#echo 'add probable consequences of the structural variants'
#python "${SVannotation_software_directory}"/programs/annotate_vcf_with_structural_variant_exon_consequence.py -i "${annotated_VCF_file}" -o "${annotated_VCF_file_with_consequences}" -g "${SVannotation_software_directory}/reference_data/UCSC_Genes_canonical_exons_20180613.bed"

##### Convert the VCF file to tab-delimited text file

#python "${SVannotation_software_directory}"/programs/convert_structural_variant_and_mobile_element_VCF_to_tab_delimited_for_excel.py -i "${annotated_VCF_file}" -o "${TMPDIR}/${unique_id}_annotated.tsv"
#cut -d$'\t' -f1-7,29,51- "${TMPDIR}/${unique_id}_annotated.tsv" > "${annotated_TSV_file}"

#####

cp "${annotated_VCF_file}" "${output_VCF_file}"

echo 'Finished!'


