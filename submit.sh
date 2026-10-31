#!/bin/bash
# Submit one condor job per list. Skips jobs whose output already exists or
# that are still in the queue, so running it again submits only what is missing.
# Usage: submit.sh <run name> <truth|reco> [max jobs]
set -euo pipefail

RUN=$1; MODE=$2; MAX=${3:-0}
HERE=$(dirname "$(realpath -e "$0")")
source "$HERE/config.sh"
RUNDIR=$HERE/runs/$RUN
[ "$HERE/ePIC_multikt2_C.so" -nt "$HERE/ePIC_multikt2.C" ] || { echo "macro not compiled: run ./build.sh" >&2; exit 1; }
# LFS pointer files are about 130 bytes, the real theory files about 55 MB
for s in ep eAu eAg; do
  [ "$(stat -Lc %s "$HERE/TH/Cth_$s.root" 2>/dev/null || echo 0)" -gt 1000000 ] ||
    { echo "TH/Cth_$s.root is missing or an LFS pointer: run git lfs pull in eic-shell" >&2; exit 1; }
done
LOGDIR=$RUNDIR/logs/$MODE
mkdir -p "$LOGDIR"

STAMP=$(date +%Y-%m-%dT%H%M%S)
ITEMS=$RUNDIR/$MODE-$STAMP.items
QUEUED=$RUNDIR/$MODE-$STAMP.queued
condor_q "$USER" -af Args Arguments > "$QUEUED"
python3 - "$RUNDIR" "$MODE" "$MAX" "$QUEUED" > "$ITEMS" <<'PY'
import csv, os, sys
rundir, mode, limit, queued = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
queued = set(open(queued).read().split())
n = 0
# by job number, so a limited batch covers every sample
for j in sorted(csv.DictReader(open(f"{rundir}/jobs.csv")), key=lambda j: j["job"]):
    out = f"{rundir}/output/{mode}/{j['sample']}/out_{j['job']}.root"
    if os.path.exists(out) or out in queued:
        continue
    if limit and n >= limit:
        break
    print(j["sys"], j["list"], out, f"{j['sample']}_{j['job']}")
    n += 1
PY
[ -s "$ITEMS" ] || { echo "nothing to submit"; exit 0; }

SUBMIT=$RUNDIR/$MODE-$STAMP.submit
cat > "$SUBMIT" <<EOF
universe            = vanilla
executable          = $HERE/run_job.sh
transfer_executable = False
arguments           = $MODE \$(sys) \$(list) \$(out)

error  = $LOGDIR/\$(tag).\$(Cluster).err
output = $LOGDIR/\$(tag).\$(Cluster).out
log    = $LOGDIR/$STAMP.log

request_cpus   = 1
request_memory = $REQUEST_MEMORY

queue sys,list,out,tag from $ITEMS
EOF

condor_submit "$SUBMIT"
