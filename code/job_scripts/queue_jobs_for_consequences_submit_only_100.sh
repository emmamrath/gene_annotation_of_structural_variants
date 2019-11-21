#!/bin/bash
set -euo pipefail

QUEUE=normal
INDIR="../vcf_TRANCHE2_annotations1_annotations2"
OUTDIR="../vcf_TRANCHE2_annotations1_annotations2_tranche2HighMedium_consequences"

cohort="MGRB"

MAXQUEUED=100
num_queued=$((MAXQUEUED+1))

COUNTER=1
for infile in "${INDIR}"/*.vcf; do

  # Only submit 100 jobs at a time
  if [ ${num_queued} -gt ${MAXQUEUED} ]; then
    set +e
    num_queued=$(qstat -u `whoami` | grep "${QUEUE}" | awk '($10 == "Q"){total+=1} END{print total+0}')
    set -e
  fi
  if [ ${num_queued} -gt ${MAXQUEUED} ]; then
    echo -e "\e[31mTerminating: queue full\e[39m"
    break
  fi

  sampleid=$(basename "${infile}")
  sampleid="${sampleid%%.*}"

  queue_file="${OUTDIR}/${sampleid}.svcons.queued"
  lock_file="${OUTDIR}/${sampleid}.svcons.lock"
  done_file="${OUTDIR}/${sampleid}.svcons.done"
  term_file="${OUTDIR}/${sampleid}.svcons.term"
  log_file="${OUTDIR}/${sampleid}.svcons.log"

  if [ -e "${queue_file}" ]; then
      echo "${sampleid} already queued"
  elif [ -e "${lock_file}" ]; then
      echo "${sampleid} already running"
  elif [ -e "${done_file}" ]; then
      echo "${sampleid} already done"
  elif [ -e "${term_file}" ]; then
      echo "${sampleid} was terminated"
  else

    qsub -z -v COHORT="${cohort}",SAMPLEID="${sampleid}",INFILE="${infile}",OUTDIR="${OUTDIR}" -N SVc"${sampleid}" step04_SV_one_sample_annotate_gridss_VCF_with_genomic_regions_consequences.pbs
    touch "${queue_file}"
    echo "${sampleid} queued"
    echo 'COUNTER =' ${COUNTER}
    COUNTER=$((COUNTER +1))
    num_queued=$((num_queued+1))
    sleep 5

  fi
done

