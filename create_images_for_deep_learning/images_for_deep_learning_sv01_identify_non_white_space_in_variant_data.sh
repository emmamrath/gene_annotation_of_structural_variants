#!/bin/bash

# Non-white space are genomic coordinates that contain structural variants.
# Input: a file containing a list of tab-delimited sample files, each file contains structural variants for a sample.
# Input: a file containing already-identified non-white space coordinates.
# Processing: For each sample file, use bedtools to take the union of the non-white file and sample file to produce the new non-white file to be used.
# Output: a file containing non-white space coordinates.

in_list_file=$1
in_non_white=$2
out_non_white=$3

module load bedtools/2.28.0

tmpdir="./tmp"
mkdir -p "${tmpdir}"
outfile_basename=$(basename $out_non_white)
outfile_basename=${outfile_basename/%????/}
tmp_output_one="${tmpdir}"/"${outfile_basename}"_non_white_output_one.bed
tmp_output_all="${tmpdir}"/"${outfile_basename}"_non_white_output_all.bed
tmp_output_all_plus_one="${tmpdir}"/"${outfile_basename}"_non_white_output_all_plus_one.bed
touch "${in_non_white}" # create file if not already created
cp "${in_non_white}" "${tmp_output_all}"

while read -r infile; do

	# SVTYPE=INS/DEL/INDEL/DUP have SVLEN field. SVLEN is -ve for DEL. SVLEN can be +ve or -ve for INDEL.
	# SVTYPE=INV does not have SVLEN field, can calculate from from END-START.
	# Choose only TRANCHE2=INTERMEDIATE/HIGH and VAF >= 0.25

	# MGRB: CHROM=3, POS=4, END=30, SVLEN=53, SVTYPE=54, TRANCHE2=91, VAF=134
	# ISKS vcf_TRANCHE2_annotations1_annotations2/*.tsv: CHROM=1, POS=2, END=30, SVLEN=51, SVTYPE=52, TRANCHE2=89, VAF=128
	cat "${infile}" | cut -d$'\t' -f1,2,30,51,52,89,128 | grep -v '^CHROM' | grep -P '\tHIGH\t|\tINTERMEDIATE\t' | awk 'BEGIN {FS="\t";OFS="\t"} {if ($7 >= 0.25) print $0}' | cut -d$'\t' -f1-3 > "${tmp_output_one}"

	cat "${tmp_output_one}" "${tmp_output_all}" | bedtools sort -i - | bedtools merge -i - > "${tmp_output_all_plus_one}"
	cp "${tmp_output_all_plus_one}" "${tmp_output_all}"

done < "$in_list_file"

cp "${tmp_output_all}" "${out_non_white}"

##### step 01 - ISKS non-white space
#cd /g/data3/zc9/users/emr913/ISKS_2018jan/hs37d5x/gridss_SV/code

#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 100 > list_images_for_deep_learning_tsv_0001.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 200 | tail -n 100 > list_images_for_deep_learning_tsv_0002.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 300 | tail -n 100 > list_images_for_deep_learning_tsv_0003.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 400 | tail -n 100 > list_images_for_deep_learning_tsv_0004.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 500 | tail -n 100 > list_images_for_deep_learning_tsv_0005.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 600 | tail -n 100 > list_images_for_deep_learning_tsv_0006.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 700 | tail -n 100 > list_images_for_deep_learning_tsv_0007.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 800 | tail -n 100 > list_images_for_deep_learning_tsv_0008.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 900 | tail -n 100 > list_images_for_deep_learning_tsv_0009.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 1000 | tail -n 100 > list_images_for_deep_learning_tsv_0010.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | tail -n 121 > list_images_for_deep_learning_tsv_0011.txt

#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0001.txt ../images_for_deep_learning/non_white_in_0001.txt ../images_for_deep_learning/non_white_out_0001.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0002.txt ../images_for_deep_learning/non_white_in_0002.txt ../images_for_deep_learning/non_white_out_0002.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0003.txt ../images_for_deep_learning/non_white_in_0003.txt ../images_for_deep_learning/non_white_out_0003.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0004.txt ../images_for_deep_learning/non_white_in_0004.txt ../images_for_deep_learning/non_white_out_0004.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0005.txt ../images_for_deep_learning/non_white_in_0005.txt ../images_for_deep_learning/non_white_out_0005.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0006.txt ../images_for_deep_learning/non_white_in_0006.txt ../images_for_deep_learning/non_white_out_0006.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0007.txt ../images_for_deep_learning/non_white_in_0007.txt ../images_for_deep_learning/non_white_out_0007.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0008.txt ../images_for_deep_learning/non_white_in_0008.txt ../images_for_deep_learning/non_white_out_0008.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0009.txt ../images_for_deep_learning/non_white_in_0009.txt ../images_for_deep_learning/non_white_out_0009.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0010.txt ../images_for_deep_learning/non_white_in_0010.txt ../images_for_deep_learning/non_white_out_0010.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0011.txt ../images_for_deep_learning/non_white_in_0011.txt ../images_for_deep_learning/non_white_out_0011.txt

#cat ../images_for_deep_learning/non_white_out_00??.txt | awk 'BEGIN {FS="\t";OFS="\t"} {if ($2==$3) print $1,$2,$2+1; else print $0}' | bedtools sort -i - | bedtools merge -i - > ../images_for_deep_learning/non_white_out_all_isks.txt

##### step 01 - MGRB non-white space
#cd /g/data3/wq2/users/emr913/MGRB_gridss_SV_2017sep/code # 2978

#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 500 > list_images_for_deep_learning_tsv_0001.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 1000 | tail -n 500 > list_images_for_deep_learning_tsv_0002.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 1500 | tail -n 500 > list_images_for_deep_learning_tsv_0003.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 2000 | tail -n 500 > list_images_for_deep_learning_tsv_0004.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | head -n 2500 | tail -n 500 > list_images_for_deep_learning_tsv_0005.txt
#ls -1 ../vcf_TRANCHE2_annotations1_annotations2/*.tsv | tail -n 478 > list_images_for_deep_learning_tsv_0006.txt

#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0001.txt ../images_for_deep_learning/non_white_in_0001.txt ../images_for_deep_learning/non_white_out_0001.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0002.txt ../images_for_deep_learning/non_white_in_0002.txt ../images_for_deep_learning/non_white_out_0002.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0003.txt ../images_for_deep_learning/non_white_in_0003.txt ../images_for_deep_learning/non_white_out_0003.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0004.txt ../images_for_deep_learning/non_white_in_0004.txt ../images_for_deep_learning/non_white_out_0004.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0005.txt ../images_for_deep_learning/non_white_in_0005.txt ../images_for_deep_learning/non_white_out_0005.txt
#./images_for_deep_learning_sv01_identify_non_white_space_in_variant_data.sh list_images_for_deep_learning_tsv_0006.txt ../images_for_deep_learning/non_white_in_0006.txt ../images_for_deep_learning/non_white_out_0006.txt

#cat ../images_for_deep_learning/non_white_out_000?.txt | awk 'BEGIN {FS="\t";OFS="\t"} {if ($2==$3) print $1,$2,$2+1; else print $0}' | bedtools sort -i - | bedtools merge -i - > ../images_for_deep_learning/non_white_out_all_mgrb.txt

#module load bedtools/2.28.0
#cat /g/data3/wq2/users/emr913/MGRB_gridss_SV_2017sep/images_for_deep_learning/non_white_out_all_mgrb.txt /g/data3/zc9/users/emr913/ISKS_2018jan/hs37d5x/gridss_SV/images_for_deep_learning/non_white_out_all_isks.txt | bedtools sort -i - | bedtools merge -i - > /g/data3/zc9/users/emr913/ISKS_2018jan/hs37d5x/gridss_SV/images_for_deep_learning/non_white_out_all_mgrb_and_isks.txt

##### step 02 - bp to pixel mapping file

#cd /g/data3/zc9/users/emr913/ISKS_2018jan/hs37d5x/gridss_SV/code
#module load python3/3.5.2
#python3 images_for_deep_learning_sv02_create_non_white_space_mapping_file.py -i ../images_for_deep_learning/non_white_out_all_mgrb_and_isks.txt -o ../images_for_deep_learning/non_white_out_all_mgrb_and_isks_map_350x350.txt

##### step 03 - create images

#cd /g/data3/zc9/users/emr913/ISKS_2018jan/hs37d5x/gridss_SV/code
#./images_for_deep_learning_sv03_image_using_non_white_space_mapping_file.sh

#cd /g/data3/wq2/users/emr913/MGRB_gridss_SV_2017sep/code
#./images_for_deep_learning_sv03_image_using_non_white_space_mapping_file.sh


