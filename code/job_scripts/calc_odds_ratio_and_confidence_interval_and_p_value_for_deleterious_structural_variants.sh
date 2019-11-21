#!/bin/bash

cohort_name_1=$1 # ISKS
cohort_name_2=$2 # MGRB
cohort_total_1=$3 # 1121
cohort_total_2=$4 # 2845
cohort_data_1=$5 # /nvme/emmrat/ISKS_2018jan/hs37d5x/gridss_SV/tsv_cohort_summary_results/ISKS_SV_output_one_line_per_gene_for_sample_variants_in_cohort.tsv
cohort_data_2=$6 # /nvme/emmrat/MGRB_2017/gridss_SV_2017sep/tsv_cohort_summary_results/MGRB_SV_output_one_line_per_gene_for_sample_variants_in_cohort.tsv
outfile=$7 # /nvme/emmrat/ISKS_2018jan/hs37d5x/gridss_SV/tsv_cohort_summary_results/ISKS_MGRB_odds_ratio_and_confidence_intervals_and_p_value.tsv

svdir="/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation/programs"
tmpdir="/nvme/emmrat/tmp_SV"

acmg_name="ACMG"
acmg_file="/nvme/emmrat/software_various/git_for_helpful_programs_and_scripts/git_for_helpful_programs_and_scripts/helpful_data_files/list_of_genes_ACMG_gene_names.txt"
acmg_extended_name="ACMG_extended"
acmg_extended_file="/nvme/emmrat/software_various/git_for_helpful_programs_and_scripts/git_for_helpful_programs_and_scripts/helpful_data_files/list_of_genes_CGC_GCR_gene_names.txt"
telo_name="Telomere_associated"
telo_file="/nvme/emmrat/software_various/git_for_helpful_programs_and_scripts/git_for_helpful_programs_and_scripts/helpful_data_files/list_of_genes_POT1etal_gene_names.txt"
aml_mds_name="AML_MDS"
aml_mds_file="${tmpdir}"/"${cohort}"_aml_mds_file.txt
echo -e "ANKRD26\nASXL1\nBRAF\nCALR\nCEBPA\nDDX41\nDNMT3A\nETV6\nFLT3\nGATA2\nIDH1\nIDH2\nJAK2\nKIT\nMPL\nNPM1\nRUNX1\nSAMD9L\nSRP72\nTERC\nTERT\nTET2\nTP53\nU2AF1\nZRSR2" > $aml_mds_file

echo 'Rscript calc_odds_ratio_and_confidence_interval_and_p_value_for_deleterious_structural_variants.R' $cohort_name_1 $cohort_name_2 $cohort_total_1 $cohort_total_2 $cohort_data_1 $cohort_data_2 $outfile $acmg_name $acmg_file $acmg_extended_name $acmg_extended_file $telo_name $telo_file $aml_mds_name $aml_mds_file
Rscript "${svdir}"/calc_odds_ratio_and_confidence_interval_and_p_value_for_deleterious_structural_variants.R "${cohort_name_1}" "${cohort_name_2}" "${cohort_total_1}" "${cohort_total_2}" "${cohort_data_1}" "${cohort_data_2}" "${outfile}" "${acmg_name}" "${acmg_file}" "${acmg_extended_name}" "${acmg_extended_file}" "${telo_name}" "${telo_file}" "${aml_mds_name}" "${aml_mds_file}"



