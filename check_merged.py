#!/usr/bin/env python3
"""Compare merged files with their job outputs. Run inside eic-shell.

    python3 check_merged.py <run name> <truth|reco>

Every histogram in runs/<run>/merged/<mode>/<sample>.root must match the sum of
out_*.root in entries and sum of weights. Prints hCounters totals and pairs per sample.
"""
import glob
import os
import sys

import ROOT

ROOT.gErrorIgnoreLevel = ROOT.kError
run, mode = sys.argv[1:3]
rundir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs", run)


def walk(d, prefix=""):
    for k in d.GetListOfKeys():
        o = k.ReadObj()
        name = prefix + k.GetName()
        if o.InheritsFrom("TDirectory"):
            yield from walk(o, name + "/")
        elif o.InheritsFrom("TH1"):
            yield name, o.GetEntries(), o.GetSumOfWeights()


fail = 0
tot_ev = tot_pairs = 0
for merged in sorted(glob.glob(f"{rundir}/merged/{mode}/*.root")):
    s = os.path.basename(merged)[:-5]
    jobs = sorted(glob.glob(f"{rundir}/output/{mode}/{s}/out_*.root"))
    want = {}
    for j in jobs:
        f = ROOT.TFile.Open(j)
        for n, e, i in walk(f):
            a = want.setdefault(n, [0.0, 0.0])
            a[0] += e
            a[1] += i
        f.Close()
    f = ROOT.TFile.Open(merged)
    got = {n: (e, i) for n, e, i in walk(f)}
    h = f.Get("hCounters")
    ev = int(h.GetBinContent(4))
    pairs = sum(e for n, (e, i) in got.items() if n.endswith("/deltaPhi"))
    decs = {n: int(e) for n, (e, i) in got.items() if n.endswith("/deltaPhi")}
    f.Close()
    bad = [n for n in want if n not in got
           or abs(got[n][0] - want[n][0]) > 1e-6 * max(1, want[n][0])
           or abs(got[n][1] - want[n][1]) > 1e-6 * max(1, abs(want[n][1]))]
    extra = [n for n in got if n not in want]
    fail += len(bad) + len(extra)
    tot_ev += ev
    tot_pairs += pairs
    print(f"{s}: jobs {len(jobs)}, hists {len(got)}, mismatched {len(bad)}, extra {len(extra)}, "
          f"events {ev}, pairs {int(pairs)} {decs}")
    for n in (bad + extra)[:5]:
        print("   ", n, want.get(n), got.get(n))
print(f"{mode}: total events {tot_ev}, pairs {int(tot_pairs)}, problems {fail}")
sys.exit(1 if fail else 0)
