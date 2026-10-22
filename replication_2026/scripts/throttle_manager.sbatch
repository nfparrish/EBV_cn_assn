#!/usr/bin/env bash
#SBATCH --job-name=throttle_mgr
#SBATCH --account=parrishlab
#SBATCH --partition=common
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=24:00:00
#SBATCH --output=/work/nfp8/EBV/logs/throttle_mgr_%j.log
#
# Hands EBI bandwidth to the next job in priority order as earlier ones drain, so the
# rebalancing does not depend on anyone being awake.
#
# Priority: the 698-sample trio CN job first (it is the scientific blocker and was being refused
# by EBI), then the phase 3 low-coverage CN job, then the GEUVADIS download.
# EBI starts refusing connections somewhere above ~35 concurrent streams, so the total is kept
# near 30. MAGE's downloader holds 20 of those until it finishes.
set -uo pipefail
TRIO=59077948; LOW=59062063; GEUV=59083173; MAGE=59026676; ONT=59067075

running() { squeue -h -j "$1" -o "%A" 2>/dev/null | grep -qc . ; }
nleft()   { squeue -h -j "$1" -t PD,R -o "%i" 2>/dev/null | wc -l ; }
setthr()  { scontrol update JobId="$1" ArrayTaskThrottle="$2" 2>&1 | grep -v "already finished" || true; }

echo "$(date -Is) throttle manager started"
while true; do
  mage_on=$(running "$MAGE" && echo 1 || echo 0)
  ont_on=$(running "$ONT"  && echo 1 || echo 0)
  trio_on=$(running "$TRIO" && echo 1 || echo 0)
  low_on=$(running "$LOW"  && echo 1 || echo 0)

  # budget of ~30 concurrent EBI streams; MAGE takes 20 while alive, ONT ~6
  budget=30
  [[ $mage_on == 1 ]] && budget=$((budget-20))
  [[ $ont_on  == 1 ]] && budget=$((budget-6))
  (( budget < 2 )) && budget=2

  if [[ $trio_on == 1 ]]; then
      t=$(( budget > 12 ? 12 : budget )); setthr "$TRIO" "$t"; budget=$((budget-t))
  fi
  if [[ $low_on == 1 && $budget -gt 0 ]]; then
      t=$(( budget > 8 ? 8 : budget )); setthr "$LOW" "$t"; budget=$((budget-t))
  fi
  if [[ $budget -gt 0 ]]; then
      t=$(( budget > 14 ? 14 : budget )); setthr "$GEUV" "$t"
  fi

  echo "$(date -Is) mage=$mage_on ont=$ont_on trio=$trio_on($(nleft $TRIO) left) low=$low_on($(nleft $LOW)) geuv=$(nleft $GEUV) -> geuv throttle target $((budget>14?14:budget))"
  running "$GEUV" || { echo "GEUVADIS finished; exiting"; break; }
  sleep 600
done
