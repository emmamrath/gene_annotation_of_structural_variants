#!/bin/bash

cohort1=$1
indir1=$2
cohort2=$3
indir2=$4
outname=$5
outdir=$6
tmpdir=$6

#cohort1='MGRB'
#indir1='/g/data3/wq2/users/emr913/MGRB_mobile_elements_Mobster_mobiome_2017jul/annotated_mobelwrap_results'
#cohort2='ISKS'
#indir2='/g/data3/zc9/users/emr913/ISKS_2018jan/hs37d5x/mobster_MEI/annotated_mobelwrap_results_samples_need_renaming'
#outname='MGRB_vs_ISKS'
#outdir='/g/data3/zc9/users/emr913/ISKS_2018jan/hs37d5x/mobster_MEI/comparing_cohort_mobelwrap_results'
#tmpdir='/g/data3/wq2/users/emr913/temp'

cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_HITGENES.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_HITGENES_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_HITEXONS.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_HITEXONS_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_HITINTRONS.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_HITINTRONS_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_HITCDSEXONS.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_HITCDSEXONS_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_HITUTR5.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_HITUTR5_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_HITUTR3.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_HITUTR3_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_HITTSS.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_HITTSS_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_FLANK1000GENES.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_FLANK1000GENES_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_FLANK3000GENES.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_FLANK3000GENES_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_FLANK10000GENES.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_FLANK10000GENES_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_HITREGULATORY.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_HITREGULATORY_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_HITMOBEL.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_HITMOBEL_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_HITCPG.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_HITCPG_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_HITSRE.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_HITSRE_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_NONREFSEQGENES.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQGENES_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_NONREFSEQEXONS.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQEXONS_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_NONREFSEQINTRONS.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQINTRONS_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_NONREFSEQCDS.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQCDS_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_NONREFSEQCDSEXONS.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQCDSEXONS_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_NONREFSEQUTR5.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQUTR5_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_NONREFSEQUTR3.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQUTR3_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_NONREFSEQTSS.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQTSS_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_FLANK1000NONREFSEQ.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_FLANK1000NONREFSEQ_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_FLANK1000NONREFSEQ.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_FLANK3000NONREFSEQ_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_FLANK10000NONREFSEQ.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_FLANK10000NONREFSEQ_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_SVREGION1000G.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_SVREGION1000G_col1.txt
cut -d$'\t' -f1 "${indir1}"/"${cohort1}"_MEI_SVCALL1000G.txt | sort | uniq > "${TMPDIR}"/"${cohort1}"_MEI_SVCALL1000G_col1.txt

cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_HITGENES.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_HITGENES_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_HITEXONS.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_HITEXONS_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_HITINTRONS.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_HITINTRONS_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_HITCDSEXONS.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_HITCDSEXONS_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_HITUTR5.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_HITUTR5_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_HITUTR3.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_HITUTR3_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_HITTSS.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_HITTSS_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_FLANK1000GENES.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_FLANK1000GENES_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_FLANK3000GENES.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_FLANK3000GENES_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_FLANK10000GENES.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_FLANK10000GENES_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_HITREGULATORY.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_HITREGULATORY_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_HITMOBEL.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_HITMOBEL_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_HITCPG.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_HITCPG_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_HITSRE.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_HITSRE_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_NONREFSEQGENES.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQGENES_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_NONREFSEQEXONS.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQEXONS_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_NONREFSEQINTRONS.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQINTRONS_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_NONREFSEQCDS.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQCDS_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_NONREFSEQCDSEXONS.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQCDSEXONS_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_NONREFSEQUTR5.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQUTR5_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_NONREFSEQUTR3.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQUTR3_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_NONREFSEQTSS.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQTSS_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_FLANK1000NONREFSEQ.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_FLANK1000NONREFSEQ_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_FLANK3000NONREFSEQ.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_FLANK3000NONREFSEQ_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_FLANK10000NONREFSEQ.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_FLANK10000NONREFSEQ_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_SVREGION1000G.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_SVREGION1000G_col1.txt
cut -d$'\t' -f1 "${indir2}"/"${cohort2}"_MEI_SVCALL1000G.txt | sort | uniq > "${TMPDIR}"/"${cohort2}"_MEI_SVCALL1000G_col1.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_HITGENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITGENES_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_HITGENES.txt > "${outdir}"/"${outname}"_MEI_HITGENES_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_HITGENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITGENES_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_HITGENES.txt > "${outdir}"/"${outname}"_MEI_HITGENES_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_HITGENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITGENES_col1.txt > "${outdir}"/"${outname}"_MEI_HITGENES_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_HITEXONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITEXONS_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_HITEXONS.txt > "${outdir}"/"${outname}"_MEI_HITEXONS_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_HITEXONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITEXONS_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_HITEXONS.txt > "${outdir}"/"${outname}"_MEI_HITEXONS_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_HITEXONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITEXONS_col1.txt > "${outdir}"/"${outname}"_MEI_HITEXONS_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_HITINTRONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITINTRONS_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_HITINTRONS.txt > "${outdir}"/"${outname}"_MEI_HITINTRONS_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_HITINTRONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITINTRONS_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_HITINTRONS.txt > "${outdir}"/"${outname}"_MEI_HITINTRONS_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_HITINTRONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITINTRONS_col1.txt > "${outdir}"/"${outname}"_MEI_HITINTRONS_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_HITCDSEXONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITCDSEXONS_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_HITCDSEXONS.txt > "${outdir}"/"${outname}"_MEI_HITCDSEXONS_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_HITCDSEXONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITCDSEXONS_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_HITCDSEXONS.txt > "${outdir}"/"${outname}"_MEI_HITCDSEXONS_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_HITCDSEXONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITCDSEXONS_col1.txt > "${outdir}"/"${outname}"_MEI_HITCDSEXONS_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_HITUTR5_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITUTR5_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_HITUTR5.txt > "${outdir}"/"${outname}"_MEI_HITUTR5_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_HITUTR5_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITUTR5_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_HITUTR5.txt > "${outdir}"/"${outname}"_MEI_HITUTR5_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_HITUTR5_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITUTR5_col1.txt > "${outdir}"/"${outname}"_MEI_HITUTR5_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_HITUTR3_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITUTR3_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_HITUTR3.txt > "${outdir}"/"${outname}"_MEI_HITUTR3_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_HITUTR3_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITUTR3_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_HITUTR3.txt > "${outdir}"/"${outname}"_MEI_HITUTR3_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_HITUTR3_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITUTR3_col1.txt > "${outdir}"/"${outname}"_MEI_HITUTR3_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_HITTSS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITTSS_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_HITTSS.txt > "${outdir}"/"${outname}"_MEI_HITTSS_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_HITTSS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITTSS_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_HITTSS.txt > "${outdir}"/"${outname}"_MEI_HITTSS_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_HITTSS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITTSS_col1.txt > "${outdir}"/"${outname}"_MEI_HITTSS_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_FLANK1000GENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK1000GENES_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_HITFLANK1000GENES.txt > "${outdir}"/"${outname}"_MEI_FLANK1000GENES_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_FLANK1000GENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK1000GENES_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_FLANK1000GENES.txt > "${outdir}"/"${outname}"_MEI_FLANK1000GENES_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_FLANK1000GENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK1000GENES_col1.txt > "${outdir}"/"${outname}"_MEI_FLANK1000GENES_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_FLANK3000GENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK3000GENES_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_FLANK3000GENES.txt > "${outdir}"/"${outname}"_MEI_FLANK3000GENES_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_FLANK3000GENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK3000GENES_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_FLANK3000GENES.txt > "${outdir}"/"${outname}"_MEI_FLANK3000GENES_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_FLANK3000GENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK3000GENES_col1.txt > "${outdir}"/"${outname}"_MEI_FLANK3000GENES_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_FLANK10000GENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK10000GENES_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_FLANK10000GENES.txt > "${outdir}"/"${outname}"_MEI_FLANK10000GENES_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_FLANK10000GENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK10000GENES_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_FLANK10000GENES.txt > "${outdir}"/"${outname}"_MEI_FLANK10000GENES_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_FLANK10000GENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK10000GENES_col1.txt > "${outdir}"/"${outname}"_MEI_FLANK10000GENES_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_HITREGULATORY_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITREGULATORY_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_HITREGULATORY.txt > "${outdir}"/"${outname}"_MEI_HITREGULATORY_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_HITREGULATORY_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITREGULATORY_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_HITREGULATORY.txt > "${outdir}"/"${outname}"_MEI_HITREGULATORY_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_HITREGULATORY_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITREGULATORY_col1.txt > "${outdir}"/"${outname}"_MEI_HITREGULATORY_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_HITMOBEL_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITMOBEL_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_HITMOBEL.txt > "${outdir}"/"${outname}"_MEI_HITMOBEL_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_HITMOBEL_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITMOBEL_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_HITMOBEL.txt > "${outdir}"/"${outname}"_MEI_HITMOBEL_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_HITMOBEL_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITMOBEL_col1.txt > "${outdir}"/"${outname}"_MEI_HITMOBEL_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_HITCPG_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITCPG_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_HITCPG.txt > "${outdir}"/"${outname}"_MEI_HITCPG_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_HITCPG_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITCPG_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_HITCPG.txt > "${outdir}"/"${outname}"_MEI_HITCPG_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_HITCPG_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITCPG_col1.txt > "${outdir}"/"${outname}"_MEI_HITCPG_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_HITSRE_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITSRE_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_HITSRE.txt > "${outdir}"/"${outname}"_MEI_HITSRE_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_HITSRE_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITSRE_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_HITSRE.txt > "${outdir}"/"${outname}"_MEI_HITSRE_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_HITSRE_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_HITSRE_col1.txt > "${outdir}"/"${outname}"_MEI_HITSRE_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQGENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQGENES_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_NONREFSEQGENES.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQGENES_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQGENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQGENES_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_NONREFSEQGENES.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQGENES_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQGENES_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQGENES_col1.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQGENES_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQEXONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQEXONS_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_NONREFSEQEXONS.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQEXONS_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQEXONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQEXONS_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_NONREFSEQEXONS.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQEXONS_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQEXONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQEXONS_col1.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQEXONS_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQINTRONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQINTRONS_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_NONREFSEQINTRONS.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQINTRONS_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQINTRONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQINTRONS_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_NONREFSEQINTRONS.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQINTRONS_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQINTRONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQINTRONS_col1.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQINTRONS_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQCDS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQCDS_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_NONREFSEQCDS.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQCDS_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQCDS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQCDS_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_NONREFSEQCDS.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQCDS_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQCDS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQCDS_col1.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQCDS_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQCDSEXONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQCDSEXONS_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_NONREFSEQCDSEXONS.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQCDSEXONS_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQCDSEXONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQCDSEXONS_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_NONREFSEQCDSEXONS.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQCDSEXONS_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQCDSEXONS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQCDSEXONS_col1.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQCDSEXONS_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQUTR5_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQUTR5_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_NONREFSEQUTR5.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQUTR5_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQUTR5_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQUTR5_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_NONREFSEQUTR5.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQUTR5_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQUTR5_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQUTR5_col1.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQUTR5_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQUTR3_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQUTR3_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_NONREFSEQUTR3.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQUTR3_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQUTR3_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQUTR3_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_NONREFSEQUTR3.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQUTR3_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQUTR3_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQUTR3_col1.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQUTR3_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQTSS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQTSS_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_NONREFSEQTSS.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQTSS_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQTSS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQTSS_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_NONREFSEQTSS.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQTSS_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_NONREFSEQTSS_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_NONREFSEQTSS_col1.txt > "${outdir}"/"${outname}"_MEI_NONREFSEQTSS_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_FLANK1000NONREFSEQ_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK1000NONREFSEQ_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_FLANK1000NONREFSEQ.txt > "${outdir}"/"${outname}"_MEI_FLANK1000NONREFSEQ_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_FLANK1000NONREFSEQ_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK1000NONREFSEQ_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_FLANK1000NONREFSEQ.txt > "${outdir}"/"${outname}"_MEI_FLANK1000NONREFSEQ_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_FLANK1000NONREFSEQ_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK1000NONREFSEQ_col1.txt > "${outdir}"/"${outname}"_MEI_FLANK1000NONREFSEQ_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_FLANK3000NONREFSEQ_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK3000NONREFSEQ_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_FLANK3000NONREFSEQ.txt > "${outdir}"/"${outname}"_MEI_FLANK3000NONREFSEQ_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_FLANK3000NONREFSEQ_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK3000NONREFSEQ_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_FLANK3000NONREFSEQ.txt > "${outdir}"/"${outname}"_MEI_FLANK3000NONREFSEQ_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_FLANK3000NONREFSEQ_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK3000NONREFSEQ_col1.txt > "${outdir}"/"${outname}"_MEI_FLANK3000NONREFSEQ_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_FLANK10000NONREFSEQ_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK10000NONREFSEQ_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_FLANK10000NONREFSEQ.txt > "${outdir}"/"${outname}"_MEI_FLANK10000NONREFSEQ_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_FLANK10000NONREFSEQ_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK10000NONREFSEQ_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_FLANK10000NONREFSEQ.txt > "${outdir}"/"${outname}"_MEI_FLANK10000NONREFSEQ_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_FLANK10000NONREFSEQ_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_FLANK10000NONREFSEQ_col1.txt > "${outdir}"/"${outname}"_MEI_FLANK10000NONREFSEQ_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_SVREGION1000G_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_SVREGION1000G_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_SVREGION1000G.txt > "${outdir}"/"${outname}"_MEI_SVREGION1000G_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_SVREGION1000G_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_SVREGION1000G_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_SVREGION1000G.txt > "${outdir}"/"${outname}"_MEI_SVREGION1000G_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_SVREGION1000G_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_SVREGION1000G_col1.txt > "${outdir}"/"${outname}"_MEI_SVREGION1000G_in_both_cohorts.txt

comm -23 "${TMPDIR}"/"${cohort1}"_MEI_SVCALL1000G_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_SVCALL1000G_col1.txt | join -1 1 -2 1 - "${indir1}"/"${cohort1}"_MEI_SVCALL1000G.txt > "${outdir}"/"${outname}"_MEI_SVCALL1000G_in_"${cohort1}"_only.txt
comm -13 "${TMPDIR}"/"${cohort1}"_MEI_SVCALL1000G_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_SVCALL1000G_col1.txt | join -1 1 -2 1 - "${indir2}"/"${cohort2}"_MEI_SVCALL1000G.txt > "${outdir}"/"${outname}"_MEI_SVCALL1000G_in_"${cohort2}"_only.txt
comm -12 "${TMPDIR}"/"${cohort1}"_MEI_SVCALL1000G_col1.txt "${TMPDIR}"/"${cohort2}"_MEI_SVCALL1000G_col1.txt > "${outdir}"/"${outname}"_MEI_SVCALL1000G_in_both_cohorts.txt


