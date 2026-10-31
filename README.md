# ePIC multi-trigger dihadron analysis on BNL condor

Runs `ePIC_multikt2.C` on the `pythia6-lhapdf` RECO files in Rucio, one condor
job per input list, inside the `eic_xl:26.07.3-stable` container. Files are read
over XRootD, mostly from BNL-XRD.

Clone onto `/gpfs02`, which the workers can read. Everything else follows from where
the clone is, so no paths need editing. Clone and commit from inside eic-shell, because
`git-lfs` is in the container but not on the submit host. `TH/`, the merged files and
`files.csv` go through LFS. A clone made without LFS gets pointer files instead; run
`git lfs pull` in eic-shell to fetch them.

## Files

| File | Purpose |
|---|---|
| `ePIC_multikt2.C` | The analysis macro |
| `TH/` | Theory correlation functions `Cth_ep`, `Cth_eAu`, `Cth_eAg`, used as saturation weights, stored with LFS |
| `config.sh` | Cuts, options, container image, bind paths, memory request |
| `build.sh` | Compiles the macro; run after any macro change |
| `make_lists.py` | Builds job lists from Rucio, leaving out `excluded_files.txt` (run in eic-shell) |
| `excluded_files.txt` | Files every run must leave out |
| `epic_unique_files.py` | Picks one copy of each file across tries; called by `make_lists.py` |
| `submit.sh` | Submits one job per list, skipping finished ones |
| `run_job.sh` | What each job runs |
| `check.py` | Compares events read with `jobs.csv`, which holds the Rucio counts (run in eic-shell) |
| `merge.sh` | `hadd` per sample |
| `check_merged.py` | Compares merged histograms with the sum of job outputs (run in eic-shell) |
| `plots/make_decade1_plots.py` | First-decade C(Δφ) plots and fits from the merged files, written to `plots/<run>/` (run in eic-shell) |
| `workflow.dot` | Workflow diagram source; `workflow.pdf`, `.png` and `.svg` are the rendered versions |
| `runs/<name>/` | Lists, logs, outputs and merged files of one run |

## Run

```bash
./build.sh                                  # outside eic-shell
# in eic-shell:
python3 make_lists.py r2                    # the five near-complete samples
# outside eic-shell:
./submit.sh r2 truth
./submit.sh r2 reco
condor_q
# in eic-shell:
python3 check.py r2 truth
python3 check.py r2 reco
# outside eic-shell:
./merge.sh r2 truth
./merge.sh r2 reco
# in eic-shell:
python3 check_merged.py r2 truth
python3 check_merged.py r2 reco
python3 plots/make_decade1_plots.py r2
```

`submit.sh` skips finished and queued jobs, so running it again submits only what is
missing. `check.py` lists missing and bad outputs. Delete a bad output and run
`submit.sh` again to redo only those jobs. Use a new run name to pick up newly finished files.

Plots read the merged files, so rerun `merge.sh` and `check_merged.py` after any
output changes. `make_decade1_plots.py` expects the three first-decade samples listed
at its top.

## Settings

From Gabriel's own commands: trigger 2-50 GeV, associate 1-50 GeV, `doSubCut` and
`cutdecay` off, `pidMode=0`. `sys` comes from the sample. Change them in `config.sh`.

## Samples

`make_lists.py` defaults to the five samples at 99.5% of target or above on
2026-10-30: `p4491`, `p4494` (ep), `p4493` (eAg), `p4492`, `p4496` (eAu). Use
`--pattern` for others.

### eAg exclusion

Run `r1` leaves out 25 one-event files from `p4493`, 25 events in total. They are
listed in `runs/r1/excluded_files.txt`; the original lists and `jobs.csv` are in
`runs/r1/backup_before_exclude/`.

Normal `p4493` files store `InclusiveKinematicseSigma` differently from these files,
and a job crashes with exit 129 when it reads from one layout into the other. ep and
eAu are not affected. The macro does not read this collection.

`r1` had them removed by hand. `make_lists.py` now leaves out everything in
`excluded_files.txt`, matched by file name, and copies that file into the run. The
other eAg samples have not been checked for the same mix of layouts. Details are in
`RUN_STATUS.md`.

## Sizes

From test jobs on the workers: about 45 s per file plus 3 minutes of start-up. The
default 100 files per job takes about 1.3 hours. 2 GB of memory is enough for analysis
jobs.

## Macro changes from Gabriel's version

- Theory files are read from `$CTH_DIR`, which `run_job.sh` sets to `TH/` in the clone.
- It returns 1 if the event loop stops before reading every entry.
