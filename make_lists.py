#!/usr/bin/env python3
"""Build per-job input lists from Rucio. Run inside eic-shell.

    python3 make_lists.py <run name> [--events-per-job 100000] [--files-per-job 100]
                          [--pattern PATTERN ...] [--exclude FILE]

epic_unique_files.py picks one copy of each file across tries. Files are then
grouped per sample with limits on both events and files per job. Standard files
hold about 200 events and event service files about 10,000, and on the workers
each file costs 40-47 s however few events it has.

Files named in excluded_files.txt are left out, matched by file name, and that file is
copied into the run. Writes runs/<run name>/lists/<sample>/list_NNNN.txt and
runs/<run name>/jobs.csv.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import defaultdict

from rucio.client import Client

HERE = os.path.dirname(os.path.abspath(__file__))
UNIQUE = os.path.join(HERE, "epic_unique_files.py")
EXCLUDE = os.path.join(HERE, "excluded_files.txt")
BASE = "RECO/26.07.3/epic_craterlake/SIDIS/pythia6-lhapdf/1.1.0"
# Samples at 99.5% of target or above on 2026-10-30
PATTERNS = [
    f"{BASE}/ep_noradcor/9x130/x_1e-3to1e-2",
    f"{BASE}/ep_noradcor/9x130/x_1e-4to1e-3",
    f"{BASE}/eAg_EPPS16_noradcor/9x115/x_1e-4to1e-3",
    f"{BASE}/eAu_EPPS16_noradcor/9x100/x_1e-2to1e-1",
    f"{BASE}/eAu_EPPS16_noradcor/9x100/x_1e-4to1e-3",
]
SYSTEMS = {"ep_noradcor": "ep", "eAu_EPPS16_noradcor": "eAu", "eAg_EPPS16_noradcor": "eAg"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--events-per-job", type=int, default=100000)
    ap.add_argument("--files-per-job", type=int, default=100)
    ap.add_argument("--pattern", nargs="+", default=PATTERNS)
    ap.add_argument("--exclude", default=EXCLUDE, help="file of URLs to leave out")
    args = ap.parse_args()

    rundir = os.path.join(HERE, "runs", args.run)
    if os.path.exists(rundir):
        sys.exit(f"{rundir} already exists; use a new run name")
    os.makedirs(rundir)
    urls = os.path.join(rundir, "files.txt")
    manifest = os.path.join(rundir, "files.csv")
    subprocess.run([sys.executable, UNIQUE, *args.pattern, "-o", urls, "--manifest", manifest],
                   check=True)

    with open(args.exclude) as f:
        excluded = {line.strip().rsplit("/", 1)[-1] for line in f
                    if line.strip() and not line.startswith("#")}
    shutil.copy(args.exclude, os.path.join(rundir, "excluded_files.txt"))

    with open(manifest) as f:
        found = list(csv.DictReader(f))
    rows = [r for r in found if r["file"].rsplit("/", 1)[-1] not in excluded]
    print(f"left out {len(found) - len(rows)} excluded files")
    client = Client()
    events = {}
    for dataset in sorted({r["dataset"] for r in rows}):
        for f in client.list_files("epic", dataset):
            events[(dataset, f["name"].rsplit("/", 1)[-1])] = f["events"]

    samples = defaultdict(list)
    for r in rows:
        n = events.get((r["dataset"], r["file"].rsplit("/", 1)[-1]))
        if not n:
            sys.exit(f"no event count in Rucio for {r['file']}")
        samples[r["sample"]].append((r["url"], n))

    jobs = []
    for sample, files in sorted(samples.items()):
        species, beam, xbin = sample.split("/")[-3:]
        name = f"{SYSTEMS[species]}_{beam}_{xbin}"
        listdir = os.path.join(rundir, "lists", name)
        os.makedirs(listdir)
        chunks, chunk, total = [], [], 0
        for url, n in files:
            if chunk and (total + n > args.events_per_job or len(chunk) >= args.files_per_job):
                chunks.append((chunk, total))
                chunk, total = [], 0
            chunk.append(url)
            total += n
        chunks.append((chunk, total))
        for i, (chunk, total) in enumerate(chunks):
            path = os.path.join(listdir, f"list_{i:04d}.txt")
            with open(path, "w") as f:
                f.write("\n".join(chunk) + "\n")
            jobs.append([name, f"{i:04d}", SYSTEMS[species], path, len(chunk), total])

    with open(os.path.join(rundir, "jobs.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "job", "sys", "list", "files", "events"])
        w.writerows(jobs)
    print(f"{len(rows)} files, {sum(j[5] for j in jobs)} events, {len(jobs)} jobs in {rundir}")


if __name__ == "__main__":
    main()
