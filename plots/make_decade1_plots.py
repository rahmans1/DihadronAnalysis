import ROOT,math,json,os,sys
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel=ROOT.kWarning; ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0);ROOT.gStyle.SetOptTitle(0)
# Usage: make_decade1_plots.py [run name]. Reads runs/<run>/merged/, writes plots/<run>/.
RUN=sys.argv[1] if len(sys.argv)>1 else 'r1'
HERE=os.path.dirname(os.path.abspath(__file__))
R=f'{HERE}/../runs/{RUN}/merged'
os.makedirs(f'{HERE}/{RUN}',exist_ok=True); os.chdir(f'{HERE}/{RUN}')
D='x_1e-4_1e-3'; REB=3
SMP={'ep':'ep_9x130_x_1e-4to1e-3','eAu':'eAu_9x100_x_1e-4to1e-3','eAg':'eAg_9x115_x_1e-4to1e-3'}
Q=['QCDC','LODIS','ResQQ','ResGQ']; G=['PGF','ResQG','ResGG']
def load(m,s):
    H={}
    t=ROOT.TFile.Open(f'{R}/{m}/{SMP[s]}.root'); d=t.Get(D)
    nt=d.Get('TrigPt').GetEntries()
    for n in ['deltaPhi']+[f'deltaPhi{p}' for p in Q+G]+[f'deltaPhi{p}Sat' for p in G]:
        H[n]=d.Get(n).Clone(f'{m}_{s}_{n}'); H[n].SetDirectory(0)
    t.Close()
    return H,nt
def fit(h,name):
    f=ROOT.TF1(name,'gaus(0)+gaus(3)+[6]',-math.pi/2,3*math.pi/2)
    m=h.GetMaximum() or 1
    f.SetParameters(m/3,0,0.4,m,math.pi,0.4,0.1)
    f.FixParameter(1,0); f.FixParameter(4,math.pi)
    f.SetParLimits(0,0,10*m+1);f.SetParLimits(3,0,10*m+1);f.SetParLimits(2,0.2,0.6);f.SetParLimits(5,0.2,0.6);f.SetParLimits(6,0,m+1)
    r=h.Fit(f,'WLRNQS' if h.GetSumOfWeights()!=h.GetEntries() else 'LRNQS')
    return f,r
data={};res={}
for m in('truth','reco'):
    for s in SMP:
        H,nt=load(m,s); h=H['deltaPhi'].Clone(f'{m}_{s}_np'); 
        data[m,s]=dict(np=h,nt=nt)
        if s!='ep':
            sat=H[f'deltaPhi{G[0]}Sat'].Clone(f'{m}_{s}_sat')
            for p in G[1:]: sat.Add(H[f'deltaPhi{p}Sat'])
            for p in Q: sat.Add(H[f'deltaPhi{p}'])
            data[m,s]['sat']=sat
out=ROOT.TFile('prelim_decade1.root','RECREATE')
def norm(h,nt):
    c=h.Clone(h.GetName()+'_C'); c.Rebin(REB); c.Scale(1./nt/c.GetBinWidth(1)); return c
keep=[]
def panel(m,s):
    ep=data[m,'ep'];a=data[m,s]
    series=[('e+p',ep['np'],ep['nt'],ROOT.kBlack,20,1),(f'e+{s[1:]}, nPDF only',a['np'],a['nt'],ROOT.kMagenta+1,24,2),(f'e+{s[1:]}, nPDF+"Saturation"',a['sat'],a['nt'],ROOT.kRed+1,20,1)]
    ymax=0;info=[]
    for i,(lab,h,nt,c,ms,ls) in enumerate(series):
        raw=h.Clone(h.GetName()+'_r'); raw.Rebin(REB)
        f,r=fit(raw,f'f_{m}_{s}_{i}')
        C=norm(h,nt); k=1./nt/C.GetBinWidth(1)
        F=ROOT.TF1(f'F_{m}_{s}_{i}','[7]*(gaus(0)+gaus(3)+[6])',-math.pi/2,3*math.pi/2)
        for j in range(7): F.SetParameter(j,f.GetParameter(j))
        F.SetParameter(7,k)
        C.SetMarkerColor(c);C.SetLineColor(c);C.SetMarkerStyle(ms);C.SetMarkerSize(.8)
        F.SetLineColor(c);F.SetLineStyle(ls);F.SetLineWidth(2);F.SetNpx(500)
        ymax=max(ymax,max(C.GetBinContent(b)+C.GetBinError(b) for b in range(1,C.GetNbinsX()+1)))
        keep.extend([C,F,raw,f]);info.append((lab,C,F))
        # away-side yield per trigger
        y=f.GetParameter(3)*abs(f.GetParameter(5))*math.sqrt(2*math.pi)*k
        ey=y*math.hypot(f.GetParError(3)/max(f.GetParameter(3),1e-9),f.GetParError(5)/f.GetParameter(5))
        res[f'{m}_{s}_{lab}']=dict(pairs=h.GetEntries(),sum=round(h.GetSumOfWeights(),1),ntrig=nt,near_A=f.GetParameter(0)*k,near_s=f.GetParameter(2),away_A=f.GetParameter(3)*k,away_s=f.GetParameter(5),base=f.GetParameter(6)*k,away_yield=y,away_yield_err=ey,status=int(r.Status()))
        C.Write();
    return info,ymax
for m in('truth','reco'):
    c=ROOT.TCanvas(f'c{m}','',1400,560); c.Divide(2,1,0.001,0.001)
    for ip,s in enumerate(('eAu','eAg')):
        c.cd(ip+1); ROOT.gPad.SetLeftMargin(.13);ROOT.gPad.SetRightMargin(.03);ROOT.gPad.SetTopMargin(.05);ROOT.gPad.SetBottomMargin(.12)
        info,ymax=panel(m,s)
        fr=info[0][1]; fr.SetTitle('');fr.GetXaxis().SetTitle('#Delta#phi');fr.GetYaxis().SetTitle('C(#Delta#phi)')
        fr.SetMinimum(0);fr.SetMaximum(ymax*1.9);fr.GetYaxis().SetTitleOffset(1.25);fr.Draw('E1')
        for lab,C,F in info: C.Draw('E1 SAME');F.Draw('L SAME')
        lg=ROOT.TLegend(.55,.71,.97,.93);lg.SetBorderSize(0);lg.SetFillStyle(0);lg.SetTextSize(.036)
        lg.SetHeader('10^{-4} #leq #it{x} < 10^{-3}')
        for lab,C,F in info: lg.AddEntry(C,lab,'lp')
        lg.Draw();keep.append(lg)
        t=ROOT.TLatex();t.SetNDC();t.SetTextFont(62);t.SetTextSize(.05);t.DrawLatex(.17,.87,'ePIC Internal')
        t.SetTextFont(42);t.SetTextSize(.034)
        e={'eAu':'9#times100','eAg':'9#times115'}[s]
        for k,l in enumerate([f'{m.capitalize()}, 26.07.3, {RUN} (preliminary)',f'e+p 9#times130, e+{s[1:]} {e} GeV',f'MC stats only, no lumi projection','1 #leq Q^{2} < 9 GeV^{2}']):
            t.DrawLatex(.17,.81-.05*k,l)
    c.SaveAs(f'decade1_{m}.png');c.SaveAs(f'decade1_{m}.pdf')
out.Close()
json.dump(res,open('fit.json','w'),indent=1)
for k,v in res.items(): print(k,{a:(round(b,4) if isinstance(b,float) else b) for a,b in v.items()})
