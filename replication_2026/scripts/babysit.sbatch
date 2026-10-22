#!/usr/bin/env bash
#SBATCH --job-name=babysit
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=48:00:00
#SBATCH --output=/work/nfp8/EBV/logs/babysit_%j.log
#
# Unattended manager for the EBI-bound jobs. Two duties:
#
# 1. REBALANCE EBI CONCURRENCY. EBI refuses connections above ~35 concurrent streams. MAGE's
#    downloader holds 20 while alive. The first trio-CN array lost 646 of 698 tasks to
#    "Connection refused" because it ran at throttle 10 against that. Budget is therefore
#    allocated in priority order: trio CN (the scientific blocker) > low-coverage CN > GEUVADIS.
#
# 2. RESUBMIT INCOMPLETE ARRAYS. Every array here is idempotent (per-sample .done markers or
#    output files), so resubmitting only retries what is missing. A task that dies on a refused
#    connection is simply retried later when bandwidth frees up.
#
# Jobs are resolved BY NAME each pass, so resubmission does not invalidate the manager.
set -uo pipefail
B=/work/nfp8/EBV
jid()   { squeue -u "$USER" -h -n "$1" -o "%A" 2>/dev/null | head -1; }
nleft() { local j; j=$(jid "$1"); [[ -n $j ]] && squeue -h -j "$j" -t PD,R -o "%i" 2>/dev/null | wc -l || echo 0; }
setthr(){ local j; j=$(jid "$1"); [[ -n $j ]] && scontrol update JobId="$j" ArrayTaskThrottle="$2" >/dev/null 2>&1; }

trio_done() { comm -12 <(awk -F'\t' '{print $1}' "$B/metadata/cram_manifest_698.tsv" | sort) \
                       <(ls "$B/cn/.done" 2>/dev/null | sort) | wc -l; }
low_done()  { ls "$B"/lowcov/*.tsv 2>/dev/null | wc -l; }
geuv_done() { ls /work/nfp8/GEUVADIS/fastq/*/.done 2>/dev/null | wc -l; }

echo "$(date -Is) babysitter started"
while true; do
  mage=$([[ -n $(jid mage_fastq) ]] && echo 1 || echo 0)
  ont=$([[ -n $(jid ebv_ont) ]] && echo 1 || echo 0)

  budget=30
  (( mage )) && budget=$((budget-20))
  (( ont ))  && budget=$((budget-6))
  (( budget < 2 )) && budget=2

  td=$(trio_done); ld=$(low_done); gd=$(geuv_done)

  # priority 1: trios -- but HOLD while MAGE is alive.
  # The 698 CRAMs exist ONLY on ENA (ftp.sra.ebi.ac.uk): not on the AWS mirror, and the
  # 1000genomes FTP carries only index files for this collection. MAGE's downloader saturates
  # that same host, and the trio job is refused even at throttle 4 -- it lost 646/698 tasks the
  # first time and kept failing on resubmission. So do not submit it at all until MAGE is done;
  # burning array tasks on refused connections achieves nothing and buries the logs.
  if (( td < 698 )) && (( mage == 0 )); then
      if [[ -z $(jid ebv_698) ]]; then
          t=$(( budget > 6 ? 6 : budget ))
          sbatch --array=1-698%${t} --job-name=ebv_698 \
            --export=ALL,MANIFEST=$B/metadata/cram_manifest_698.tsv \
            "$B/scripts/ebv_cn_slice.sbatch" >/dev/null 2>&1
          echo "$(date -Is) resubmitted trio array at throttle $t ($td/698 done)"
      else
          setthr ebv_698 $(( budget > 10 ? 10 : (budget<2?2:budget) ))
      fi
      budget=$(( budget - 6 )); (( budget < 0 )) && budget=0
  elif (( td < 698 )); then
      echo "$(date -Is) trio held: MAGE still running, ENA would refuse ($td/698 done)"
  fi

  # priority 2: low-coverage CN
  if (( ld < 2499 )); then
      if [[ -z $(jid ebv_low) ]]; then
          sbatch --array=1-2499%4 --job-name=ebv_low \
            --export=ALL,MANIFEST=$B/metadata/lowcov_full.tsv,TAG=lowcov \
            "$B/scripts/ebv_cn_alt.sbatch" >/dev/null 2>&1
          echo "$(date -Is) resubmitted low-coverage array ($ld/2499 done)"
      elif (( budget > 0 )); then
          setthr ebv_low $(( budget > 8 ? 8 : budget )); budget=$((budget-4))
      fi
  fi

  # priority 3: GEUVADIS download gets whatever remains
  (( budget < 1 )) && budget=1
  setthr geuv_fastq $(( budget > 14 ? 14 : budget ))

  # LOCAL arrays (no EBI bandwidth): just keep them topped up. These are the ones that can
  # silently under-run -- the BCR evidence array exited COMPLETED on ~280 samples because it ran
  # before TRUST4 had produced their reports, and nothing requeued it afterwards.
  bd=$(ls "$B"/bcr_evidence/*.tsv 2>/dev/null | wc -l)
  if (( bd < 731 )) && [[ -z $(jid bcr_evid) ]] \
     && (( $(ls "$B"/bcr/*/*_report.tsv 2>/dev/null | wc -l) > bd )); then
      sbatch "$B/scripts/bcr_evidence.sbatch" >/dev/null 2>&1
      echo "$(date -Is) resubmitted BCR evidence array ($bd/731 done)"
  fi
  # WGS BCR extended to all 2,504 (AWS-sourced, so no EBI contention). Idempotent on the
  # per-sample TRUST4 report, so resubmission only fills gaps.
  ba=$(ls "$B"/bcr/*/*_report.tsv 2>/dev/null | wc -l)
  if (( ba < 2504 )) && [[ -z $(jid bcr_all) ]]; then
      sbatch --array=1-2504%12 --job-name=bcr_all \
        --export=ALL,MANIFEST=$B/metadata/bcr_manifest_all.tsv \
        "$B/scripts/bcr_wgs.sbatch" >/dev/null 2>&1
      echo "$(date -Is) resubmitted full WGS-BCR array ($ba/2504 done)"
  fi

  gb=$(ls "$B"/gcbias/*.tsv 2>/dev/null | wc -l)
  if (( gb < 2502 )) && [[ -z $(jid ebv_gcb) ]]; then
      sbatch "$B/scripts/viral_gcbias.sbatch" >/dev/null 2>&1
      echo "$(date -Is) resubmitted viral GC-bias array ($gb/2502 done)"
  fi

  echo "$(date -Is) mage=$mage ont=$ont | trio $td/698 | low $ld/2499 | geuv $gd/666 | geuv_thr=$(( budget > 14 ? 14 : budget ))"
  (( td >= 698 && ld >= 2499 && gd >= 666 )) && { echo "all complete; exiting"; break; }
  sleep 600
done
