#!/bin/bash
# One condor job: run ePIC_multikt2.C on one list, inside the container.
# Usage: run_job.sh <truth|reco> <sys> <list> <output.root>
# The output appears at its final path only if the macro returned 0.
set -euo pipefail

MODE=$1; SYS=$2; LIST=$3; OUT=$4
HERE=$(dirname "$(realpath -e "$0")")
source "$HERE/config.sh"
export CTH_DIR=${CTH_DIR:-$HERE/TH}

case "$MODE" in
  truth) USE_RECO=false ;;
  reco)  USE_RECO=true ;;
  *) echo "mode must be truth or reco" >&2; exit 1 ;;
esac

WORK=${_CONDOR_SCRATCH_DIR:-$(mktemp -d)}
TMP=$WORK/out.root
ARGS="$TRIG_LOW,$TRIG_HIGH,$ASSO_LOW,$ASSO_HIGH,\"$LIST\",\"$TMP\",\"$SYS\",$DO_SUBCUT,$CUTDECAY,$USE_RECO,$PID_MODE"

# root -q exits 0 even when the macro returns 1, so read the printed return value
singularity exec -B "$BINDS" "$IMAGE" \
  root -l -b -q "$HERE/ePIC_multikt2.C+($ARGS)" 2>&1 | tee "$WORK/root.log"
grep -q '^(int) 0$' "$WORK/root.log"

mkdir -p "$(dirname "$OUT")"
cp "$TMP" "$OUT.part"
mv "$OUT.part" "$OUT"
