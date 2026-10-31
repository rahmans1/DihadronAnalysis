# Run r1 status: pick up from here

Updated 2026-10-31, about 11:00 EDT. Steps and file roles are in `README.md`.
Check the queue before trusting anything below.

## Where it stands

Run `r1`: `ePIC_multikt2.C` on five `pythia6-lhapdf` 1.1.0 RECO samples from Rucio.
197,267 files, 39,186,659 events, 1,975 jobs per mode.

The first pass finished on 2026-10-31 at about 03:10 EDT: 1,953 of 1,975 jobs per
mode succeeded. The 22 failed jobs per mode were all eAg. They were rerun and all
exited 0 by 07:54 EDT. `check.py` prints `missing 0, bad 0` for both modes, with
39,186,659 of 39,186,659 events read.

| Cluster | Mode | Jobs | State |
|---|---|---:|---|
| 320273 | truth | 22 | rerun, done |
| 320274 | reco | 22 | rerun, done |
| 320267, 320269 | truth | 1,975 | done (22 exit 129) |
| 320268, 320270 | reco | 1,975 | done (22 exit 129) |


## The eAg failures

Exit 129 was a segfault in `podio::ROOTReader::getCollectionBuffers`. It was not a
storage problem: all files tested were reachable over XRootD, with the right size and
event counts.

- Normal eAg files store `InclusiveKinematicseSigma` with 1,305 branches in `events`.
  The sample's 25 one-event files use `InclusiveKinematicseSigma_objIdx`, 1,304
  branches. A chain that goes from the normal layout to the other one crashes.
- Reproduced in `eic_xl:26.07.3-stable` with neighbouring files from list 0080.
- ep and eAu are unaffected: their sampled normal files and all their one-event files
  use the 1,304 layout. Neighbouring-file chains from all four samples passed.
- Why normal eAg files differ is unknown.
- The macro does not read this collection, so finished jobs are unaffected.
- `pc6799` and `pc6800` are also eAg and were not checked. Check their layouts before
  running on them.

The fix removed the 25 files from the 22 lists and corrected `jobs.csv`, losing 25
events. The removed files are in `runs/r1/excluded_files.txt`; the original lists and
`jobs.csv` are in `runs/r1/backup_before_exclude/`. New runs leave them out
automatically through `excluded_files.txt`.

## Samples

| Physics tag | Output folder in `runs/r1/output/<mode>/` | Files | Events | Jobs |
|---|---|---:|---:|---:|
| `p4491` | `ep_9x130_x_1e-3to1e-2` | 39,899 | 7,936,923 | 399 |
| `p4494` | `ep_9x130_x_1e-4to1e-3` | 37,888 | 7,529,282 | 379 |
| `p4493` | `eAg_9x115_x_1e-4to1e-3` | 25,000 | 4,975,000 | 251 |
| `p4492` | `eAu_9x100_x_1e-2to1e-1` | 73,157 | 14,502,407 | 732 |
| `p4496` | `eAu_9x100_x_1e-4to1e-3` | 21,323 | 4,243,047 | 214 |

Only these near-complete samples are in the run, as Sakib asked. `p4493` gained one
more file in Rucio after the lists were made (199 events), so a later run could add it.

## Settings

From Gabriel's own commands: trigger 2-50 GeV, associate 1-50 GeV, `doSubCut` and
`cutdecay` off, `pidMode=0`, both truth and reco. Container
`eic_xl:26.07.3-stable`, 1 core, 2 GB.

## Next steps

1. Give the merged files in `runs/r1/merged/<mode>/<sample>.root` to Gabriel.

Merging finished on 2026-10-31 at about 10:57 EDT. `check.py` was rerun first and
passed for both modes. `check_merged.py` then found all 86 histograms in each of the
10 merged files equal to the sum of their job outputs, with the event and pair totals
below.

## Results so far

Final r1 pairs. Totals are 1,478 truth and 1,031 reco.

| Sample | Events read | Truth pairs | Reco pairs |
|---|---:|---:|---:|
| ep 9x130, x 1e-4 to 1e-3 | 7,529,282 | 328 | 198 |
| ep 9x130, x 1e-3 to 1e-2 | 7,936,923 | 336 | 273 |
| eAu 9x100, x 1e-4 to 1e-3 | 4,243,047 | 228 | 147 |
| eAu 9x100, x 1e-2 to 1e-1 | 14,502,407 | 269 | 216 |
| eAg 9x115, x 1e-4 to 1e-3 | 4,975,000 | 317 | 197 |

- Sakib's 283-1667 pairs are per-decade counts from the 50M EVGEN study, before the
  |z1-z2| cut (`../trim_hepmc_context.md`). After that cut and scaling to r1 events,
  truth pairs agree: for example eAu x 1e-4 to 1e-3 expects about 227 and has 228.
- The intended plot is taken from Gabriel's `CDeltaPhi2_ePIC.C`: C(Δφ) for e+p, e+A
  nPDF only and e+A nPDF+saturation, per x decade. Gabriel has not confirmed this for
  these files.
- Only the first decade has all three systems in r1. The other two wait on `pc6798`,
  `pc6799`, `pc6796` and `pc6800`.
- Preliminary first-decade plots are in `plots/r1/`, made from the merged files by
  `plots/make_decade1_plots.py`. Reading the job outputs instead gave identical
  histograms and fit values. They use MC statistics with no lumi
  projection, 20 bins, fit widths limited to 0.2-0.6, and e+p at 9x130 in both panels.
  The ePIC style file did not load.
- Ask Gabriel about the 1 fb⁻¹ scaling, and whether reco should be on the plot.

## Test results

- All eight test jobs (`runs/smoke_std`, `runs/smoke_es`, clusters 320262-320265)
  passed `check.py`.
- Standard files took 40-47 s each, plus about 3 minutes of start-up per job. The
  default 100 files per job should take about 1.3 hours.

## Not yet tested

- First-pass jobs may have stalled on reads: one failed job used 17 minutes of CPU in
  2.9 hours of wall time. Keep an eye on job time and `epicxrd1` load.

## Things to watch

- The macro's output has `hCounters` at the top. `nEvents` is a `TObjString`, so use
  `hCounters` for totals after `hadd`.
- `root -q` exits 0 even when the macro returns 1, so `run_job.sh` reads the printed
  return value.
- In zsh, `C="1 2"; condor_q $C` passes one argument and prints nothing, which fooled
  two watchers into reporting an empty queue. List the clusters out in full.
- `TH/` and the outputs sit inside the `hepmc3ascii2root` git repo, untracked and not
  ignored. Do not commit the ROOT files.
