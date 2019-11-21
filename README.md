Gene annotation of structural variants
======================================

The first input to this pipeline is a vcf file genomic breakends (BND records), for which genomic structural variants (SVs) are called (insertions, deletions, indels, duplications, inversions). It also merges multiple samples into one vcf, matching SVs within 100bp of each other. This software was written specifically for structural variant BND records produced by the [GRIDSS](https://genome.cshlp.org/content/27/12/2050) software.  

The second input to this pipeline is a vcf file containing genomic DNA structural variants (INS, DEL, INDEL, DUP, INV; output from above), which are annotated for their relationship to genes. Relationships include HITGENES, HITEXONS, HITINTRONS, HITUTR5, HITUTR3, and more.  

Gene locations were taken from [UCSC Tables](https://genome.ucsc.edu/cgi-bin/hgTables) and updated versions may be substituted.  

This pipeline was created and used for the [MGRB project](https://sgc.garvan.org.au/initiatives). The [MGRB paper](https://www.biorxiv.org/content/10.1101/473348v1) can be considered the citation for this pipeline.  

Please note that this repository contains the component programs and scripts used for the MGRB project, and does not represent a ready-to-install-and-run pipeline.  

### Pipeline components

#### Pipeline scripts to call INS/DEL/INDEL/DUP/INV
* code/step00_SV_filter_gridss_vcf_and_add_tranche2.sh
* code/step01_SV_multiple_samples_call_INS_DEL_INDEL_DUP_INV_from_gridss_output.sh
* code/step01_SV_one_sample_call_INS_DEL_INDEL_DUP_INV_from_gridss_output.sh
* code/step02_SV_merge_multiple_samples_to_one_VCF_using_inhouse_pgm.sh
* code/step03_SV_fuzzy_merge_multiple_samples_already_in_one_VCF.sh
* code/step03_part1_SV_fuzzy_merge_by_chromosome_for_one_SVtype.sh
* code/step03_part1_SV_fuzzy_merge_by_non_number_chromosome_for_one_SVtype.sh
* code/step03_part2_SV_fuzzy_merge_for_one_SVtype.sh
* code/step03_part2_SV_merge_the_fuzzy_merge_the_chromosomes_for_one_SVtype.sh
* code/step03_part3_SV_fuzzy_merge_multiple_samples_already_in_one_VCF_merge.sh

#### Pipeline scripts to annotate with genes
* code/job_scripts/step04_add_vcf_sampleid_header.sh
* code/job_scripts/step04_convert_VCF_INFO_to_tab_delimited_and_intact_fields.sh
* code/job_scripts/step04_MEI_annotate_merged_VCF_with_genomic_regions_for_all_cohort.sh
* code/job_scripts/step04_MEI_annotate_merged_VCF_with_genomic_regions_for_one_sample.sh
* code/job_scripts/step04_MEI_compare_cohort_genes_lists.sh
* code/job_scripts/step04_parallel_SV_multiple_samples_annotate_gridss_VCF_with_genomic_regions_consequences.sh
* code/job_scripts/step04_parallel_SV_one_sample_annotate_gridss_VCF_with_genomic_regions_consequences.sh
* code/job_scripts/step04_part02_SV_multiple_samples_annotate_gridss_VCF_with_bicseq_overlaps.sh
* code/job_scripts/step04_part02_SV_one_sample_annotate_gridss_VCF_with_bicseq_overlaps.sh
* code/job_scripts/step04_part03_SV_fisher_test_for_hit_counts.sh
* code/job_scripts/step04_SV_annotate_all_VCF_with_genomic_regions_2.sh
* code/job_scripts/step04_SV_annotate_all_VCF_with_genomic_regions.sh
* code/job_scripts/step04_SV_annotate_one_VCF_with_genomic_regions_2.sh
* code/job_scripts/step04_SV_annotate_one_VCF_with_genomic_regions.sh
* code/job_scripts/calc_odds_ratio_and_confidence_interval_and_p_value_for_deleterious_structural_variants.sh
* code/job_scripts/convert_SV_VCF_to_tab_delimited_spreadsheet_with_selected_INFO_fields.sh
* code/job_scripts/intersect_DEL_INDEL_to_Swergold_1990_Speek_2001_Line1_promoters.sh
* code/job_scripts/queue_jobs_for_consequences_submit_only_100_entry_point.sh
* code/job_scripts/queue_jobs_for_consequences_submit_only_100.sh
* code/job_scripts/queue_jobs_for_inhouse_structural_variant_annotation.sh

#### Pipeline programs to call INS/DEL/INDEL/DUP/INV and annotate with genes
* code/programs/add_breakend_key_to_vcf_record.py
* code/programs/add_key_to_VCF_to_find_inversions.py
* code/programs/add_sample_key_to_VCF_to_merge_multiple_samples.py
* code/programs/add_tranche2_and_filter_gridss_vcf.py
* code/programs/annotate_VCF_structural_variants_with_bedtools_intersect_genes.py
* code/programs/annotate_vcf_with_structural_variant_exon_consequence.py
* code/programs/convert_sample_GT_in_VCF_file.py
* code/programs/convert_structural_variant_and_mobile_element_VCF_to_tab_delimited_for_excel.py
* code/programs/convert_VCF_BND_records_to_simple_structural_variant_VCF_records.py
* code/programs/convert_VCF_INFO_one_sample_to_tab_delimited_and_intact_fields.py
* code/programs/convert_VCF_structural_variants_to_BED_format.py
* code/programs/extract_vcf_info_annotation.py
* code/programs/filter_out_large_variants.py
* code/programs/filter_out_small_variants.py
* code/programs/generate_VCF_IDs_for_VCF_variants.py
* code/programs/identify_inversions_in_VCF_with_key.py
* code/programs/merge_a_multi_sample_SV_VCF_for_similar_position_variants.py
* code/programs/merge_one_sample_VCFs_to_one_multisample_VCF.py
* code/programs/sort_VCF_by_order_found_in_input_list.py
* code/programs/split_UCSC_color_genes_to_genes_exons_introns.py
* code/programs/split_UCSC_refseq_genes_to_genes_exons_introns.py
* code/programs/choose_canonical_exons_from_UCSC_data_for_SV_consequences.R
* code/programs/add_tranche_to_vcf.awk

#### Pipeline components to prepare structural variants as images for deep learning
* create_images_for_deep_learning/images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh
* create_images_for_deep_learning/images_for_deep_learning_sv02_create_non_white_space_mapping_file.py
* create_images_for_deep_learning/images_for_deep_learning_sv02_create_non_white_space_mapping_file.sh
* create_images_for_deep_learning/images_for_deep_learning_sv03_image_using_non_white_space_mapping_file.py
* create_images_for_deep_learning/images_for_deep_learning_sv03_image_using_non_white_space_mapping_file.sh
* create_images_for_deep_learning/convert_gridss_SV_to_image.py
* create_images_for_deep_learning/plot_sv_square.pbs

### References

* MGRB [(Medical Genome Reference Bank)](https://sgc.garvan.org.au/initiatives) data is a whole-genome data resource of 4000 healthy elderly individuals ([manuscript](https://www.biorxiv.org/content/10.1101/473348v1).  

* GRIDSS software
GRIDSS: sensitive and specific genomic rearrangement detection using positional de Bruijn graph assembly.  
Cameron DL, Schröder J, Penington JS, Do H, Molania R, Dobrovic A, Speed TP, Papenfuss AT.  
Genome Res. 2017 Dec;27(12):2050-2060. doi: 10.1101/gr.222109.117. Epub 2017 Nov 2.  
PMID: [29097403](https://www.ncbi.nlm.nih.gov/pubmed/?term=29097403) PMCID: PMC5741059 DOI: [10.1101/gr.222109.117](https://genome.cshlp.org/content/27/12/2050)  




