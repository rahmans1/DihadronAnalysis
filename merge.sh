#!/bin/bash
# hadd the job outputs of each sample into runs/<run>/merged/<mode>/<sample>.root.
# Usage: merge.sh <run name> <truth|reco>     (run check.py first)
set -euo pipefail
HERE=$(dirname "$(realpath -e "$0")")
source "$HERE/config.sh"
RUNDIR=$HERE/runs/$1; MODE=$2
mkdir -p "$RUNDIR/merged/$MODE"
for d in "$RUNDIR/output/$MODE"/*/; do
  s=$(basename "$d")
  singularity exec -B "$BINDS" "$IMAGE" hadd -f "$RUNDIR/merged/$MODE/$s.root" "$d"out_*.root
done
