#!/bin/bash
set -euo pipefail

COHORT=$1
INDIR=$2
OUTDIR=$3
OUT_PREFIX=$4

TMPDIR="/nvme/emmrat/tmp_for_spreadsheets"
SVDIR="/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation/programs"

# Create a dummy VCF header of the INFO fields that we want to output as columns.
temphdr="${TMPDIR}"/"${COHORT}"_temp_INFO_hdr.txt
:>"${temphdr}"

echo '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">' >> "${temphdr}"
echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' >> "${temphdr}"
echo '##INFO=<ID=HITGENES,Number=A,Type=String,Description="Intersecting gene">' >> "${temphdr}"
echo '##INFO=<ID=HITEXONS,Number=A,Type=String,Description="Intersecting gene exon">' >> "${temphdr}"
echo '##INFO=<ID=HITINTRONS,Number=A,Type=String,Description="Intersecting gene intron">' >> "${temphdr}"
echo '##INFO=<ID=HITCDS,Number=A,Type=String,Description="Intersecting CDS region">' >> "${temphdr}"
echo '##INFO=<ID=HITCDSEXONS,Number=A,Type=String,Description="Intersecting CDS exons">' >> "${temphdr}"
echo '##INFO=<ID=HITUTR5,Number=A,Type=String,Description="Intersecting UTR-5 region">' >> "${temphdr}"
echo '##INFO=<ID=HITUTR3,Number=A,Type=String,Description="Intersecting UTR-3 region">' >> "${temphdr}"
echo '##INFO=<ID=HITTSS,Number=A,Type=String,Description="Intersecting transcription start site and 1500 bp downstream">' >> "${temphdr}"
echo '##INFO=<ID=FLANK1000GENES,Number=A,Type=String,Description="Intersecting gene flank of 1000 bp (on either of but not including the gene)">' >> "${temphdr}"
echo '##INFO=<ID=FLANK3000GENES,Number=A,Type=String,Description="Intersecting gene flank of 3000 bp (on either of but not including the gene)">' >> "${temphdr}"
echo '##INFO=<ID=FLANK10000GENES,Number=A,Type=String,Description="Intersecting gene flank of 10000 bp (on either of but not including the gene)">' >> "${temphdr}"
echo '##INFO=<ID=HITREGULATORY,Number=A,Type=String,Description="Intersecting regulatory region">' >> "${temphdr}"
echo '##INFO=<ID=HITMOBEL,Number=A,Type=String,Description="Intersecting existing mobile element region">' >> "${temphdr}"
echo '##INFO=<ID=HITCPG,Number=A,Type=String,Description="Intersecting CpG island">' >> "${temphdr}"
echo '##INFO=<ID=HITSRE,Number=A,Type=String,Description="Intersecting subtelomeric repeat region">' >> "${temphdr}"
echo '##INFO=<ID=NONREFSEQGENES,Number=A,Type=String,Description="Intersecting non-RefSeq gene">' >> "${temphdr}"
echo '##INFO=<ID=NONREFSEQEXONS,Number=A,Type=String,Description="Intersecting non-RefSeq gene exon">' >> "${temphdr}"
echo '##INFO=<ID=NONREFSEQINTRONS,Number=A,Type=String,Description="Intersecting non-RefSeq gene intron">' >> "${temphdr}"
echo '##INFO=<ID=NONREFSEQCDS,Number=A,Type=String,Description="Intersecting non-RefSeq CDS region">' >> "${temphdr}"
echo '##INFO=<ID=NONREFSEQCDSEXONS,Number=A,Type=String,Description="Intersecting non-RefSeq CDS exons">' >> "${temphdr}"
echo '##INFO=<ID=NONREFSEQUTR5,Number=A,Type=String,Description="Intersecting non-RefSeq UTR-5 region">' >> "${temphdr}"
echo '##INFO=<ID=NONREFSEQUTR3,Number=A,Type=String,Description="Intersecting non-RefSeq UTR-3 region">' >> "${temphdr}"
echo '##INFO=<ID=NONREFSEQTSS,Number=A,Type=String,Description="Intersecting non-RefSeq transcription start site and 1500 bp downstream">' >> "${temphdr}"
echo '##INFO=<ID=FLANK1000NONREFSEQ,Number=A,Type=String,Description="Intersecting non-RefSeq gene flank of 1000 bp (on either of but not including the gene)">' >> "${temphdr}"
echo '##INFO=<ID=FLANK3000NONREFSEQ,Number=A,Type=String,Description="Intersecting non-RefSeq gene flank of 3000 bp (on either of but not including the gene)">' >> "${temphdr}"
echo '##INFO=<ID=FLANK10000NONREFSEQ,Number=A,Type=String,Description="Intersecting non-RefSeq gene flank of 10000 bp (on either of but not including the gene)">' >> "${temphdr}"
echo '##INFO=<ID=SVREGION1000G,Number=A,Type=String,Description="1000 Genomes Project Phase 3 SV analysis, Sudmant et al. 2015, germline variant regions (CNV,INS,INV).">' >> "${temphdr}"
echo '##INFO=<ID=SVCALL1000G,Number=A,Type=String,Description="1000 Genomes Project Phase 3 SV analysis, Sudmant et al. 2015, germline variant calls (CNV,DEL,DUP,INS,INV).">' >> "${temphdr}"
echo '##INFO=<ID=TRANCHE,Number=.,Type=String,Description="Quality category of GRIDSS structural variant calls. Values are LOW INTERMEDIATE HIGH">' >> "${temphdr}"
echo '##INFO=<ID=SV_CANONICAL_CONSEQUENCE,Number=.,Type=String,Description="The theoretical consequence to the canoncial transcript of an exon-overlapping DUP/DEL/INDEL. Values are PROBABLY_DELETERIOUS, PROBABLY_BENIGN">' >> "${temphdr}"
echo '##INFO=<ID=NUM_GENES_SPANNED,Number=1,Type=Integer,Description="Number of HITEXONS spanned by this variant">' >> "${temphdr}"
echo '##INFO=<ID=VARIANT_DEPTH_OVER_NEIGHBOUR_DEPTHS,Number=1,Type=Float,Description="samtools depth of variant region divided by samtools depth of regions on either side, for DUP, DEL, or INDEL.">' >> "${temphdr}"
echo '##INFO=<ID=VARIANT_DEPTH_OVER_LEFT_NEIGHBOUR_DEPTH,Number=1,Type=Float,Description="samtools depth of variant region divided by samtools depth of preceding region, for DUP, DEL, or INDEL.">' >> "${temphdr}"
echo '##INFO=<ID=VARIANT_DEPTH_OVER_RIGHT_NEIGHBOUR_DEPTH,Number=1,Type=Float,Description="samtools depth of variant region divided by samtools depth of following region, for DUP, DEL, or INDEL.">' >> "${temphdr}"
echo '##INFO=<ID=VAF_ASSUMED_FROM_ALT_AND_DEPTHS,Number=1,Type=Float,Description="Assumed VAF of 0.5 or 1 for high depth DUP. Assumed VAF of 0.5 or 1 for low depth DEL/INDEL.">' >> "${temphdr}"

for infile in "${INDIR}"/*.vcf; do

	infile_basename=$(basename "${infile}")
	IFS=$'_' read -r -a array <<< "${infile_basename}"
	sampleid="${array[0]}"
	infile2="${TMPDIR}"/"${sampleid}"_SV_with_INFO_hdrs.vcf
	outfile="${OUTDIR}"/"${sampleid}"_"${OUT_PREFIX}"_SV_columns.tsv
	tmp_outfile1="${TMPDIR}"/"${sampleid}"_"${OUT_PREFIX}"_SV_columns_temp1.tsv
        tmp_outfile2="${TMPDIR}"/"${sampleid}"_"${OUT_PREFIX}"_SV_columns_temp2.tsv
        tmp_outfile3="${TMPDIR}"/"${sampleid}"_"${OUT_PREFIX}"_SV_columns_temp3.tsv

	cat "${temphdr}" "${infile}" > "${infile2}"

	echo 'python convert_structural_variant_and_mobile_element_VCF_to_tab_delimited_for_excel.py -i '"${infile2}"' -o '"${tmp_outfile1}"
	python "${SVDIR}"/convert_structural_variant_and_mobile_element_VCF_to_tab_delimited_for_excel.py -i "${infile2}" -o "${tmp_outfile1}"

	grep '^CHROM' "${tmp_outfile1}" | awk 'BEGIN {FS="\t";OFS="\t"} {print "COHORT", "SAMPLE", $0}' > "${tmp_outfile2}"
	grep -v '^CHROM' "${tmp_outfile1}" | awk -v awkcohort="$COHORT" -v awksampleid="$sampleid" 'BEGIN {FS="\t";OFS="\t"} {print awkcohort, awksampleid, $0}' > "${tmp_outfile3}"
	cat "${tmp_outfile2}" "${tmp_outfile3}" > "${outfile}"
done

