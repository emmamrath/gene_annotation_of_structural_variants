#!/bin/bash
set -euo pipefail

cohort_id=$1
INDIR=$2
OUTDIR=$3
TMPDIR=$4

#cohort_id="MGRB"
#INDIR="/g/data3/wq2/users/emr913/MGRB_mobile_elements_Mobster_mobiome_2017jul"
#OUTDIR='/g/data3/wq2/users/emr913/MGRB_mobile_elements_Mobster_mobiome_2017jul/annotated_mobelwrap_results'
#TMPDIR='/g/data3/wq2/users/emr913/temp'

SVannotation_software_directory='/nvme/emmrat/software_various/inhouse_structural_variant_annotation/structural_variant_annotation'

:>"${TMPDIR}"/"${cohort_id}"_MEI_HITGENES.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_HITEXONS.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_HITINTRONS.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_HITCDSEXONS.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_HITUTR5.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_HITUTR3.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_HITTSS.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_FLANK1000GENES.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_FLANK3000GENES.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_FLANK10000GENES.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_HITREGULATORY.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_HITMOBEL.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_HITCPG.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_HITSRE.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQGENES.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQEXONS.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQINTRONS.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQCDS.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQCDSEXONS.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQUTR5.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQUTR3.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQTSS.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_FLANK1000NONREFSEQ.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_FLANK3000NONREFSEQ.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_FLANK10000NONREFSEQ.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_SVREGION1000G.txt
:>"${TMPDIR}"/"${cohort_id}"_MEI_SVCALL1000G.txt

for infile in "${INDIR}"/MEI_*.vcf; do

	filename=$(basename "${infile}")
	IFS='_' read -r -a array <<< "$filename"
	sample_id="${array[1]}"

	echo 'step04_MEI_annotate_merged_VCF_with_genomic_regions_for_one_sample.sh' $sample_id
	input_VCF_file="${infile}"
	annotated_VCF_file="${OUTDIR}"/MEI_"${sample_id}"_annotated.vcf
	annotated_TSV_file="${OUTDIR}"/MEI_"${sample_id}"_annotated.tsv

        if [ -e "${annotated_TSV_file}" ]; then
                :
        else

		./step04_MEI_annotate_merged_VCF_with_genomic_regions_for_one_sample.sh "${sample_id}" "${input_VCF_file}" "${annotated_VCF_file}" "${annotated_TSV_file}" "${TMPDIR}"

		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info HITGENES | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_HITGENES.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info HITEXONS | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_HITEXONS.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info HITINTRONS | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_HITINTRONS.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info HITCDSEXONS | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_HITCDSEXONS.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info HITUTR5 | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_HITUTR5.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info HITUTR3 | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_HITUTR3.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info HITTSS | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_HITTSS.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info FLANK1000GENES | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_FLANK1000GENES.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info FLANK3000GENES | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_FLANK3000GENES.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info FLANK10000GENES | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_FLANK10000GENES.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info HITREGULATORY | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_HITREGULATORY.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info HITMOBEL | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_HITMOBEL.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info HITCPG | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_HITCPG.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info HITSRE | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_HITSRE.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info NONREFSEQGENES | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQGENES.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info NONREFSEQEXONS | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQEXONS.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info NONREFSEQINTRONS | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQINTRONS.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info NONREFSEQCDS | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQCDS.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info NONREFSEQCDSEXONS | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQCDSEXONS.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info NONREFSEQUTR5 | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQUTR5.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info NONREFSEQUTR3 | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQUTR3.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info NONREFSEQTSS | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQTSS.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info FLANK1000NONREFSEQ | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_FLANK1000NONREFSEQ.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info FLANK3000NONREFSEQ | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_FLANK3000NONREFSEQ.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info FLANK10000NONREFSEQ | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_FLANK10000NONREFSEQ.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info SVREGION1000G | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_SVREGION1000G.txt
		cat "${annotated_VCF_file}" | python "${SVannotation_software_directory}"/programs/extract_vcf_info_annotation.py -info SVCALL1000G | sort | uniq >> "${TMPDIR}"/"${cohort_id}"_MEI_SVCALL1000G.txt

	fi
done

echo 'get unique info field values for this cohort' $cohort_id
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_HITGENES.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_HITGENES.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_HITEXONS.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_HITEXONS.txt 
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_HITINTRONS.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_HITINTRONS.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_HITCDSEXONS.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_HITCDSEXONS.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_HITUTR5.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_HITUTR5.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_HITUTR3.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_HITUTR3.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_HITTSS.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_HITTSS.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_FLANK1000GENES.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_FLANK1000GENES.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_FLANK3000GENES.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_FLANK3000GENES.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_FLANK10000GENES.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_FLANK10000GENES.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_HITREGULATORY.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_HITREGULATORY.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_HITMOBEL.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_HITMOBEL.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_HITCPG.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_HITCPG.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_HITSRE.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_HITSRE.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQGENES.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_NONREFSEQGENES.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQEXONS.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_NONREFSEQEXONS.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQINTRONS.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_NONREFSEQINTRONS.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQCDS.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_NONREFSEQCDS.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQCDSEXONS.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_NONREFSEQCDSEXONS.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQUTR5.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_NONREFSEQUTR5.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQUTR3.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_NONREFSEQUTR3.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_NONREFSEQTSS.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_NONREFSEQTSS.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_FLANK1000NONREFSEQ.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_FLANK1000NONREFSEQ.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_FLANK3000NONREFSEQ.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_FLANK3000NONREFSEQ.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_FLANK10000NONREFSEQ.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_FLANK10000NONREFSEQ.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_SVREGION1000G.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_SVREGION1000G.txt
sed -e 's/|/\n/g' "${TMPDIR}"/"${cohort_id}"_MEI_SVCALL1000G.txt | sort | uniq -c | tr -s ' ' | awk '{FS=" ";OFS="\t"} {print $2, $1}' > "${OUTDIR}"/"${cohort_id}"_MEI_SVCALL1000G.txt

echo 'Finished the cohort'


