#!/usr/bin/env python3
"""Compare outputs with jobs.csv. Run inside eic-shell.

    python3 check.py <run name> <truth|reco>

Each output must have read exactly the events Rucio lists for its input files.
Prints the jobs to redo; delete a bad output and run submit.sh again.
"""

import csv
import os
import sys

import ROOT

run, mode = sys.argv[1:3]
rundir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs", run)
done = missing = bad = 0
expected = read = 0
for j in csv.DictReader(open(f"{rundir}/jobs.csv")):
    out = f"{rundir}/output/{mode}/{j['sample']}/out_{j['job']}.root"
    expected += int(j["events"])
    if not os.path.exists(out):
        print("missing", j["sample"], j["job"])
        missing += 1
        continue
    f = ROOT.TFile.Open(out)
    h = f.Get("hCounters") if f and not f.IsZombie() else None
    n = int(h.GetBinContent(4)) if h else -1
    f and f.Close()
    if n != int(j["events"]):
        print("bad", out, "read", n, "expected", j["events"])
        bad += 1
        continue
    done += 1
    read += n
print(f"ok {done}, missing {missing}, bad {bad}; events {read} of {expected}")
sys.exit(1 if missing or bad else 0)
