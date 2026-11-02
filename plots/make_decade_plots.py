#!/usr/bin/env python3
"""C(dphi) plots and fits for all three x decades, combining runs r1 and r2.

    python3 plots/make_decade_plots.py [truth] [reco]     (default truth; run in eic-shell)

Each sample comes from runs/<run>/merged/<mode>/ if merged, otherwise from the sum of
its job outputs so far. A sample with jobs still missing is marked on the plot.
C(dphi) is per trigger, so missing jobs cost statistics but do not change the scale.
Writes decade<n>_<mode>.{png,pdf}, prelim_<mode>.root and fit_<mode>.json to plots/r1_r2/.
"""
import csv, glob, json, math, os, sys
import ROOT

ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kWarning; ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)

MODES = sys.argv[1:] or ['truth']
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = f'{HERE}/../runs'
OUT = f'{HERE}/r1_r2'
os.makedirs(OUT, exist_ok=True); os.chdir(OUT)

REB = 3
Q = ['QCDC', 'LODIS', 'ResQQ', 'ResGQ']; G = ['PGF', 'ResQG', 'ResGG']
NAMES = ['deltaPhi'] + [f'deltaPhi{p}' for p in Q + G] + [f'deltaPhi{p}Sat' for p in G]
BEAM = {'ep': 'ep_9x130', 'eAu': 'eAu_9x100', 'eAg': 'eAg_9x115'}
ENERGY = {'eAu': '9#times100', 'eAg': '9#times115'}
# decade number: (directory, sample suffix, legend header, run of each system)
DECADES = {
    1: ('x_1e-4_1e-3', '1e-4to1e-3', '10^{-4} #leq #it{x} < 10^{-3}', dict(ep='r1', eAu='r1', eAg='r1')),
    2: ('x_1e-3_1e-2', '1e-3to1e-2', '10^{-3} #leq #it{x} < 10^{-2}', dict(ep='r1', eAu='r2', eAg='r2')),
    3: ('x_1e-2_1e-1', '1e-2to1e-1', '10^{-2} #leq #it{x} < 10^{-1}', dict(ep='r2', eAu='r1', eAg='r2')),
}


def job_count(run, sample):
    with open(f'{RUNS}/{run}/jobs.csv') as f:
        return sum(r['sample'] == sample for r in csv.DictReader(f))


def read_file(path, dname, H, tag):
    t = ROOT.TFile.Open(path)
    if not t or t.IsZombie():
        raise RuntimeError(f'cannot read {path}')
    d = t.Get(dname)
    nt = d.Get('TrigPt').GetEntries()
    events = t.Get('hCounters').GetBinContent(4)
    for n in NAMES:
        h = d.Get(n)
        if n in H:
            H[n].Add(h)
        else:
            H[n] = h.Clone(f'{tag}_{n}'); H[n].SetDirectory(0)
    t.Close()
    return nt, events


def load(mode, dec, sysname):
    dname, suffix, _, runs = DECADES[dec]
    run = runs[sysname]
    sample = f'{BEAM[sysname]}_x_{suffix}'
    total = job_count(run, sample)
    merged = f'{RUNS}/{run}/merged/{mode}/{sample}.root'
    paths = [merged] if os.path.exists(merged) else \
        sorted(glob.glob(f'{RUNS}/{run}/output/{mode}/{sample}/out_*.root'))
    done = total if paths == [merged] else len(paths)
    if not paths:
        raise RuntimeError(f'no input for {sample}')
    H = {}; nt = events = 0
    tag = f'{mode}_{dec}_{sysname}'
    for p in paths:
        a, b = read_file(p, dname, H, tag)
        nt += a; events += b
    print(f'{tag}: {sample} from {run}, jobs {done}/{total}, events {int(events)}, '
          f'triggers {int(nt)}, pairs {int(H["deltaPhi"].GetEntries())}')
    return dict(H=H, nt=nt, run=run, done=done, total=total, events=events)


def fit(h, name):
    f = ROOT.TF1(name, 'gaus(0)+gaus(3)+[6]', -math.pi/2, 3*math.pi/2)
    m = h.GetMaximum() or 1
    f.SetParameters(m/3, 0, 0.4, m, math.pi, 0.4, 0.1)
    f.FixParameter(1, 0); f.FixParameter(4, math.pi)
    f.SetParLimits(0, 0, 10*m+1); f.SetParLimits(3, 0, 10*m+1)
    f.SetParLimits(2, 0.2, 0.6); f.SetParLimits(5, 0.2, 0.6); f.SetParLimits(6, 0, m+1)
    r = h.Fit(f, 'WLRNQS' if h.GetSumOfWeights() != h.GetEntries() else 'LRNQS')
    return f, r


def saturated(H, name):
    sat = H[f'deltaPhi{G[0]}Sat'].Clone(name)
    for p in G[1:]:
        sat.Add(H[f'deltaPhi{p}Sat'])
    for p in Q:
        sat.Add(H[f'deltaPhi{p}'])
    return sat


keep = []


def panel(mode, dec, ep, a, s, res, out):
    tag = f'{mode}_{dec}_{s}'
    nuc = s[1:]
    series = [('e+p', ep['H']['deltaPhi'], ep['nt'], ROOT.kBlack, 20, 1),
              (f'e+{nuc}, nPDF only', a['H']['deltaPhi'], a['nt'], ROOT.kMagenta+1, 24, 2),
              (f'e+{nuc}, nPDF+"Saturation"', saturated(a['H'], f'{tag}_sat'), a['nt'], ROOT.kRed+1, 20, 1)]
    ymax = 0; info = []
    for i, (lab, h, nt, c, ms, ls) in enumerate(series):
        raw = h.Clone(f'{tag}_{i}_r'); raw.Rebin(REB)
        f, r = fit(raw, f'f_{tag}_{i}')
        C = h.Clone(f'{tag}_{i}_C'); C.Rebin(REB); C.Scale(1./nt/C.GetBinWidth(1))
        k = 1./nt/C.GetBinWidth(1)
        F = ROOT.TF1(f'F_{tag}_{i}', '[7]*(gaus(0)+gaus(3)+[6])', -math.pi/2, 3*math.pi/2)
        for j in range(7):
            F.SetParameter(j, f.GetParameter(j))
        F.SetParameter(7, k)
        C.SetMarkerColor(c); C.SetLineColor(c); C.SetMarkerStyle(ms); C.SetMarkerSize(.8)
        F.SetLineColor(c); F.SetLineStyle(ls); F.SetLineWidth(2); F.SetNpx(500)
        ymax = max(ymax, max(C.GetBinContent(b)+C.GetBinError(b) for b in range(1, C.GetNbinsX()+1)))
        keep.extend([C, F, raw, f, h]); info.append((lab, C, F))
        y = f.GetParameter(3)*abs(f.GetParameter(5))*math.sqrt(2*math.pi)*k
        ey = y*math.hypot(f.GetParError(3)/max(f.GetParameter(3), 1e-9), f.GetParError(5)/f.GetParameter(5))
        res[f'decade{dec}_{s}_{lab}'] = dict(
            pairs=h.GetEntries(), sum=round(h.GetSumOfWeights(), 1), ntrig=nt,
            near_A=f.GetParameter(0)*k, near_s=f.GetParameter(2), away_A=f.GetParameter(3)*k,
            away_s=f.GetParameter(5), base=f.GetParameter(6)*k, away_yield=y, away_yield_err=ey,
            status=int(r.Status()))
        out.cd(); C.Write()
    return info, ymax


for m in MODES:
    res = {}
    out = ROOT.TFile(f'prelim_{m}.root', 'RECREATE')
    for dec, (_, _, header, _) in DECADES.items():
        data = {s: load(m, dec, s) for s in BEAM}
        c = ROOT.TCanvas(f'c{m}{dec}', '', 1400, 560); c.Divide(2, 1, 0.001, 0.001)
        for ip, s in enumerate(('eAu', 'eAg')):
            c.cd(ip+1)
            ROOT.gPad.SetLeftMargin(.14); ROOT.gPad.SetRightMargin(.03)
            ROOT.gPad.SetTopMargin(.05); ROOT.gPad.SetBottomMargin(.12)
            info, ymax = panel(m, dec, data['ep'], data[s], s, res, out)
            fr = info[0][1]; fr.GetXaxis().SetTitle('#Delta#phi'); fr.GetYaxis().SetTitle('C(#Delta#phi)')
            fr.SetMinimum(0); fr.SetMaximum(ymax*1.9); fr.GetYaxis().SetTitleOffset(1.45); fr.Draw('E1')
            for lab, C, F in info:
                C.Draw('E1 SAME'); F.Draw('L SAME')
            lg = ROOT.TLegend(.60, .73, .97, .93); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextSize(.034)
            lg.SetHeader(header)
            for lab, C, F in info:
                lg.AddEntry(C, lab, 'lp')
            lg.Draw(); keep.append(lg)
            t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(62); t.SetTextSize(.05); t.DrawLatex(.17, .87, 'ePIC Internal')
            t.SetTextFont(42); t.SetTextSize(.034)
            ep_run, a_run = data['ep']['run'], data[s]['run']
            lines = [f'{m.capitalize()}, 26.07.3 (preliminary)',
                     f'e+p 9#times130 ({ep_run}), e+{s[1:]} {ENERGY[s]} GeV ({a_run})',
                     'MC stats only, no lumi projection',
                     '1 #leq Q^{2} < 9 GeV^{2}']
            short = [f'{n} {d["done"]}/{d["total"]}' for n, d in (('e+p', data['ep']), (f'e+{s[1:]}', data[s]))
                     if d['done'] < d['total']]
            if short:
                lines.append('Partial: ' + ', '.join(short) + ' jobs')
            for k, l in enumerate(lines):
                t.DrawLatex(.17, .81-.05*k, l)
        c.SaveAs(f'decade{dec}_{m}.png'); c.SaveAs(f'decade{dec}_{m}.pdf')
    out.Close()
    json.dump(res, open(f'fit_{m}.json', 'w'), indent=1)
    for k, v in res.items():
        print(k, f"pairs {int(v['pairs'])} trig {int(v['ntrig'])} near_s {v['near_s']:.3f} away_s {v['away_s']:.3f} "
                 f"yield {v['away_yield']:.5f} +- {v['away_yield_err']:.5f} status {v['status']}")
