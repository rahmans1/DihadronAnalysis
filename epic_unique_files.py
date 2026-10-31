#!/usr/bin/env python3
"""List the unique, streamable files of an ePIC sample across all production tries.

A campaign task's first attempt writes to

    RECO/<version>/<config>/<sample path>

and every full or residual rerun (.tryN) writes to

    RECO/<version>/<config>/tryN/<sample path>

or, for a background-mixed sample, below its background tag:

    RECO/<version>/<config>/<bkg tag>/tryN/<sample path>

The file name is set by the work unit (<stem>.<chunk>.eicrecon.edm4eic.root),
so the same name in two tries is the same work unit produced twice. This tool
finds the base dataset and every tryN dataset for each requested sample, keeps
one copy of each file name (from the latest try that has an available disk
replica), drops files with no available disk replica (ghosts, tape-only), and
prints xrootd URLs ready for ROOT/uproot.

Run it inside eic-shell, whose rucio client is preconfigured with the
read-only `eicread` account:

    ./eic-shell -- python3 epic_unique_files.py \\
        RECO/26.07.3/epic_craterlake/DIS/BeAGLE1.03.02-3.1/eAu/9x100/q2_1to10 \\
        -o files.txt --manifest files.csv

Any path component may contain shell-style wildcards, e.g.
'RECO/26.07.*/epic_craterlake/DIS/BeAGLE1.03.02-3.1/eAu/9x100/*'. When a
sample exists in several matching campaigns, only the latest campaign is
listed (26.07.3 over 26.07.1); tries are deduplicated within that campaign.
"""

import argparse
import csv
import fnmatch
import re
import sys
from collections import defaultdict

from rucio.client import Client

SCOPE = "epic"
# The production disk storages, in preference order. Other RSEs (tape,
# logs, the closed ASGC-XRD, non-production sites) are used only when
# named with --rse.
DEFAULT_RSES = ["BNL-XRD", "EIC-XRD"]
CHUNK = 1000
TRY_RE = re.compile(r"^try(\d+)$")
# <stem>.<chunk>[.p<attempt mark>].eicrecon.edm4eic.root: a diverted
# registration carries the attempt's mark after the chunk and is another
# copy of the same work unit (pcs.commands.RECO_LFN_TAIL_RE).
DIVERTED_RE = re.compile(r"^(.*\.\d{4})\.p\w+(\.eicrecon\.edm4eic\.root)$")


def work_unit(base):
    """The work unit a file name delivers, with any diverted mark removed."""
    return DIVERTED_RE.sub(r"\1\2", base)


def parse_sample(arg):
    """Split 'RECO/<version>/<config>/<rest>' into (version, config, rest) patterns."""
    name = arg.split(":", 1)[1] if arg.startswith(SCOPE + ":") else arg
    parts = name.strip("/").split("/")
    if len(parts) < 4 or parts[0] != "RECO":
        sys.exit(f"error: expected RECO/<version>/<config>/<sample>, got {arg!r}")
    if any(TRY_RE.match(p) for p in parts[3:]):
        sys.exit(f"error: give the sample without its tryN component: {arg!r}")
    return parts[1], parts[2], "/".join(parts[3:])


def version_key(version):
    """Order campaign versions numerically: 26.07.10 > 26.07.3 > 26.07."""
    return [int(p) if p.isdigit() else -1 for p in re.split(r"[.-]", version)]


def find_datasets(client, version, config, rest):
    """Return {(config, sample): (version, older_versions, [(try, dataset), ...])}.

    When a sample exists in several campaigns matching `version`, only the
    latest campaign is kept; its tries are sorted oldest first.
    """
    by_version = defaultdict(lambda: defaultdict(list))
    names = set()
    # A rerun's tryN segment follows the config, or follows a background
    # tag prefix for background-mixed samples (TAG_PREFIX=<bkg>/tryN), so
    # it is looked for at every depth of the sample path.
    segs = rest.split("/")
    patterns = [f"/RECO/{version}/{config}/{rest}"] + [
        f"/RECO/{version}/{config}/" + "/".join(segs[:k] + ["try*"] + segs[k:])
        for k in range(len(segs))]
    for pattern in patterns:
        names.update(client.list_dids(SCOPE, {"name": pattern}, did_type="dataset"))
    for name in names:
        # list_dids wildcards may cross '/', so recheck each component
        parts = name.strip("/").split("/")
        ver, cfg, tail = parts[1], parts[2], parts[3:]
        at = next((k for k, p in enumerate(tail) if TRY_RE.match(p)), None)
        try_nr = 1 if at is None else int(TRY_RE.match(tail[at]).group(1))
        sample = "/".join(tail if at is None else tail[:at] + tail[at + 1:])
        if (fnmatch.fnmatchcase(ver, version) and fnmatch.fnmatchcase(cfg, config)
                and fnmatch.fnmatchcase(sample, rest)
                and sample.count("/") == rest.count("/")):
            by_version[(cfg, sample)][ver].append((try_nr, name))
    found = {}
    for key, versions in by_version.items():
        ordered = sorted(versions, key=version_key)
        found[key] = (ordered[-1], ordered[:-1], sorted(versions[ordered[-1]]))
    return found


def available_replicas(client, names, rses):
    """Map file name -> (rse, root url) for the most preferred RSE holding it."""
    best = {}
    rank = {rse: i for i, rse in enumerate(rses)}
    for i in range(0, len(names), CHUNK):
        dids = [{"scope": SCOPE, "name": n} for n in names[i:i + CHUNK]]
        for rep in client.list_replicas(dids, schemes=["root"],
                                        rse_expression="|".join(rses)):
            choices = [(rank[rse], rse, pfns[0])
                       for rse, pfns in rep.get("rses", {}).items()
                       if rse in rank and pfns]
            if choices:
                _, rse, url = min(choices)
                best[rep["name"]] = (rse, url)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples", nargs="+",
                    help="RECO/<version>/<config>/<sample path>, wildcards allowed")
    ap.add_argument("--rse", action="append", dest="rses",
                    help=f"disk RSE to read from, in preference order (default: {' '.join(DEFAULT_RSES)})")
    ap.add_argument("--prefer", choices=["latest", "first"], default="latest",
                    help="which try supplies a file produced by more than one try (default: latest)")
    ap.add_argument("-n", "--max-files", type=int,
                    help="keep at most this many files per sample")
    ap.add_argument("-o", "--output", help="write URLs here instead of stdout")
    ap.add_argument("--manifest", help="also write a CSV with try, dataset, size, checksum, RSE")
    args = ap.parse_args()
    rses = args.rses or DEFAULT_RSES

    client = Client()
    rows = []
    for arg in args.samples:
        samples = find_datasets(client, *parse_sample(arg))
        if not samples:
            print(f"warning: no datasets match {arg}", file=sys.stderr)
        for (config, sample), (version, older, datasets) in sorted(samples.items()):
            prefix = f"/RECO/{version}/{config}"
            # file name -> [(try, dataset, bytes, adler32)], one entry per try
            copies = defaultdict(list)
            for try_nr, ds in datasets:
                for f in client.list_files(SCOPE, ds):
                    base = work_unit(f["name"].rsplit("/", 1)[-1])
                    copies[base].append((try_nr, ds, f["name"], f["bytes"], f["adler32"]))

            all_names = [c[2] for cs in copies.values() for c in cs]
            replicas = available_replicas(client, all_names, rses)

            picked, dup, missing = [], 0, 0
            for base in sorted(copies):
                cands = [c for c in copies[base] if c[2] in replicas]
                dup += len(copies[base]) > 1
                if not cands:
                    missing += 1
                    continue
                c = (min if args.prefer == "first" else max)(cands)
                rse, url = replicas[c[2]]
                picked.append(dict(sample=sample, campaign=version,
                                   file=c[2].rsplit("/", 1)[-1],
                                   try_nr=c[0], dataset=c[1],
                                   bytes=c[3], adler32=c[4], rse=rse, url=url))
            if args.max_files:
                picked = picked[:args.max_files]
            rows.extend(picked)

            print(f"\n{prefix}/{sample}", file=sys.stderr)
            if older:
                print(f"  skipped older campaigns: {' '.join(older)}", file=sys.stderr)
            for try_nr, ds in datasets:
                n = sum(1 for cs in copies.values() for c in cs if c[1] == ds)
                print(f"  try{try_nr:<3} {n:>8} files  {SCOPE}:{ds}", file=sys.stderr)
            per_rse = defaultdict(int)
            for r in picked:
                per_rse[r["rse"]] += 1
            print(f"  unique names {len(copies):>7}   in >1 try {dup:>6}   "
                  f"no disk replica {missing:>6}", file=sys.stderr)
            print(f"  selected     {len(picked):>7}   "
                  f"{sum(r['bytes'] for r in picked) / 1e9:.1f} GB   "
                  + "  ".join(f"{k}:{v}" for k, v in sorted(per_rse.items())), file=sys.stderr)

    out = open(args.output, "w") if args.output else sys.stdout
    for r in rows:
        print(r["url"], file=out)
    if args.output:
        out.close()
    if args.manifest:
        with open(args.manifest, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["url"])
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
