#!/bin/bash
set -euo pipefail

#INDIR="/g/data3/wq2/results/phase2/hs37d5x/bwa_IR_BQSR_gridss-1.4.1"
INDIR="/g/data3/wq2/users/emr913/MGRB_structural_variants_inhouse_annotations_2017jul"
OUTDIR="/g/data3/wq2/users/emr913/MGRB_structural_variants_inhouse_annotations_2017jul"

# AAAAA.gridss.vcf.gz
# AAAAA.gridss.vcf.gz.done
# AAAAB.gridss.vcf.gz
# AAAAB.gridss.vcf.gz.done
# AAAAD.gridss.vcf.gz
# AAAAD.gridss.vcf.gz.done
# AAAAF.gridss.vcf.gz
# AAAAF.gridss.vcf.gz.done
# AAAAG.gridss.vcf.gz
# AAAAG.gridss.vcf.gz.done

for infile in "${INDIR}"/*.gridss.vcf.gz; do

    vcf_gz_file="${infile}"
    #vcf_file="${vcf_gz_file%.gz}"
    sampleid=$(basename "${infile}")
    sampleid="${sampleid%%.*}"

    queue_file="${OUTDIR}/${sampleid}.SVannotation.queued"
    lock_file="${OUTDIR}/${sampleid}.SVannotation.lock"
    done_file="${OUTDIR}/${sampleid}.SVannotation.done"
    term_file="${OUTDIR}/${sampleid}.SVannotation.term"
    log_file="${OUTDIR}/${sampleid}.SVannotation.log"

    if [ -e "${queue_file}" ]; then
        echo "${sampleid} already queued"
    elif [ -e "${lock_file}" ]; then
        echo "${sampleid} already running"
    elif [ -e "${done_file}" ]; then
        echo "${sampleid} already done"
    elif [ -e "${term_file}" ]; then
        echo "${sampleid} was terminated"
    else
        qsub -z -v SAMPLEID="${sampleid}",VCF_GZ_FILE="${vcf_gz_file}",OUTDIR="${OUTDIR}" -N SV${sampleid} inhouse_structural_variant_annotation.pbs
        touch "${queue_file}"
        echo "${sampleid} queued"
    fi
done

