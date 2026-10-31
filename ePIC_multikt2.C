//
// ePIC_multikt2.C -- clean version of ePIC_multikt1.C
//
// Multiple-trigger dihadron Delta-phi correlation on ePIC podio/EDM4eic files. The CGC saturation
// weight Cth_eA/Cth_ep (5D interpolation of ../TH/Cth_*.root) is applied to the gluon-initiated
// channels (PGF, ResQG, ResGG).
//
// Differences w.r.t. ePIC_multikt1.C:
//  - every analysis histogram/counter is booked once per x_Bj decade and filled according to the
//    event x, so all decades come out of a single pass over the events
//  - only phase-space QA plus the histograms and counters needed in post-processing are kept
//  - saturation correction 2 (G1/G2 from qT and Qs) is removed; the 5D Cth weight is kept
//
// Output layout:
//  /                     inputFileList, cuts, config, nEvents (TObjString)
//                        hCounters (sum over decades), Q2_x, Q2_xg (all events, before cuts),
//                        Q2vx, x_hist, y_hist (accepted events)
//  /x_1e-4_1e-3/ ...     one directory per decade, identical content:
//     hCounters            Ntrig, NtrigEvt, Nevt, nEventsRead (raw counts, summed by hadd)
//     procTot              Pythia process id of the accepted events
//     deltaPhi             all pairs, unweighted
//     deltaPhi<Proc>       per process, unweighted (Proc = PGF QCDC LODIS ResQQ ResQG ResGQ ResGG)
//     deltaPhi<Proc>Sat    PGF, ResQG, ResGG weighted with Cth_eA/Cth_ep (= unweighted for ep)
//     TrigPt, TrigZ, TrigPtVsEta, AssoPt, AssoZ, AssoPtVsEta, zCorre   trigger/associate phase space
//
// pT is measured in the gamma*-p CM frame (photon along +z); eta is the lab pseudorapidity.
#include <TSystem.h>
#include <TROOT.h>
#include <TString.h>
#include <TObjString.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THnSparse.h>
#include <TLorentzVector.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>

//==========================================================================================
// ePIC input layer (podio / EDM4hep / EDM4eic) -- replaces eic-smear EventPythia/ParticleMC
//
// The analysis loop below still uses "event->..." and "particle->..." exactly as before.
// erhic::EventPythia / erhic::ParticleMC are now small stand-in classes filled by EpicReader
// from one podio Frame ("events" category) per entry.
//
//  TRUTH mode (useReco=false): hadrons = generator-level final state (generatorStatus==1) of the
//        head-on generator record (MCParticlesHeadOnFrameNoBeamFX), i.e. the same frame the old
//        Pythia-only analysis used. DIS kinematics = Pythia values stored in the Frame parameters.
//  RECO  mode (useReco=true):  hadrons = reconstructed charged particles (lab frame, with crossing
//        angle), scattered electron = ScatteredElectronsEMinusPz, DIS kinematics = electron
//        method from the reconstructed electron + true beams. Pythia-side information (process,
//        tgtparton, xtgtparton, ...) always comes from the generator.
//==========================================================================================
R__LOAD_LIBRARY(libpodio)
R__LOAD_LIBRARY(libpodioRootIO)
R__LOAD_LIBRARY(libedm4hep)
R__LOAD_LIBRARY(libedm4eic)

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <TVector2.h>
#include <TVector3.h>
using namespace std;

#include "podio/Frame.h"
#include "podio/ROOTReader.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4eic/ReconstructedParticleCollection.h"
#include "edm4eic/MCRecoParticleAssociationCollection.h"

namespace erhic
{
  //stand-in for erhic::ParticleMC (same member names as used by the analysis)
  struct ParticleMC
  {
    int    I    = 0;    //index in the record (0-based; replaces the eic-smear line number)
    int    id   = 0;    //PDG code
    int    KS   = 0;    //1 = final-state candidate for the analysis
    int    gs   = 0;    //generatorStatus (MC record only)
    int    orig = -1;   //index of the (first) mother in the record, -1 if none
    double px=0, py=0, pz=0, E=0, m=0;
    double eta=0;
    double z=0;         //z = (P.p)/(P.q)
    double ptVsGamma=0; //pt relative to the virtual photon direction (lab frame)
    TLorentzVector Get4Vector() const { return TLorentzVector(px,py,pz,E); }
  };

  //stand-in for erhic::EventPythia
  class EventPythia
  {
  public:
    //Pythia-side information (podio Frame parameters)
    int    process   = 0;
    int    tgtparton = 0;
    double xtgtparton= 0, EBrems = 0, SigRadCor = 0;
    //DIS kinematics (truth mode: generator values; reco mode: electron method)
    double QSquared=0, x=0, y=0, WSquared=0, nu=0;

    std::vector<ParticleMC> tracks; //record used by the analysis loop
    std::vector<ParticleMC> ho;     //generator record, head-on frame (used by the saturation weights)
    int  iBeamE=-1, iBeamP=-1, iScatE=-1;
    ParticleMC recoE;
    bool hasRecoE = false;
    bool useReco  = false;

    void Reset()
    {
      tracks.clear(); ho.clear();
      iBeamE=iBeamP=iScatE=-1; hasRecoE=false;
      process=tgtparton=0; xtgtparton=EBrems=SigRadCor=0.;
      QSquared=x=y=WSquared=nu=0.;
    }
    int GetNTracks() const { return (int)tracks.size(); }
    int GetProcess() const { return process; }
    ParticleMC* GetTrack(int i)   { return &tracks.at(i); }
    ParticleMC* GetHOTrack(int i) { return &ho.at(i); }
    ParticleMC* GetBeamElectron() { return iBeamE>=0 ? &tracks[iBeamE] : NULL; }
    ParticleMC* GetBeamProton()   { return iBeamP>=0 ? &tracks[iBeamP] : NULL; }
    ParticleMC* GetScatteredElectron()
    {
      if(useReco) return hasRecoE ? &recoE : NULL;
      return iScatE>=0 ? &tracks[iScatE] : NULL;
    }
  };
}

struct EpicConfig
{
  bool        useReco   = false;
  //PID for reco hadrons: 0 = PDG of the associated generator particle (perfect PID)
  //                      1 = PDG stored in the reco particle itself
  int         pidMode   = 0;
  std::string recoColl  = "ReconstructedChargedParticles";
  std::string assocColl = "ReconstructedChargedParticleAssociations"; //must match recoColl
  std::string scatEColl = "ScatteredElectronsEMinusPz";
  std::string mcLabColl = "MCParticles";                        //lab frame (crossing angle)
  std::string mcHOColl  = "MCParticlesHeadOnFrameNoBeamFX";     //head-on frame, no beam effects
  std::string truthColl = "MCParticlesHeadOnFrameNoBeamFX";     //record used in TRUTH mode
};

class EpicReader
{
public:
  EpicReader(const std::vector<std::string>& files, const EpicConfig& c) : cfg(c)
  {
    reader.openFiles(files);
    nEntries = reader.getEntries("events");
  }
  long GetEntries() const { return nEntries; }
  bool Load(long entry, erhic::EventPythia* ev);

private:
  podio::ROOTReader reader;
  EpicConfig cfg;
  long nEntries = 0;
  bool checked  = false;
  bool hasHO    = true;
  bool hasScatE = true;
  bool hasAssoc = true;
  int  nWarn    = 0;

  static double ParD(const podio::Frame& f, const std::string& key, double def=0.)
  {
    auto v = f.getParameter<std::string>(key);   //Pythia info is stored as strings
    return v ? atof(v->c_str()) : def;
  }

  static void FillRecord(const edm4hep::MCParticleCollection& mc, int nGen, bool markFinal,
                         std::vector<erhic::ParticleMC>& out)
  {
    out.clear();
    out.reserve(nGen+64);
    for(int i=0; i<nGen && i<(int)mc.size(); ++i)
    {
      auto p = mc[i];
      erhic::ParticleMC t;
      t.I  = i;
      t.id = p.getPDG();
      t.gs = p.getGeneratorStatus();
      auto mom = p.getMomentum();
      t.px = mom.x; t.py = mom.y; t.pz = mom.z;
      t.m  = p.getMass();
      t.E  = sqrt(t.px*t.px + t.py*t.py + t.pz*t.pz + t.m*t.m);
      for(const auto& par : p.getParents()) { t.orig = par.getObjectID().index; break; }
      t.KS = (markFinal && t.gs==1) ? 1 : 0;
      t.eta = TLorentzVector(t.px,t.py,t.pz,t.E).Eta();
      out.push_back(t);
    }
  }

  static int FindBeam(const std::vector<erhic::ParticleMC>& r, int pdg)
  {
    for(const auto& p : r) if(p.id==pdg && p.gs==4) return p.I;
    return -1;
  }

  //final-state electron whose mother is the scattered-electron line (pdg 11, status 21)
  static int FindScatE(const std::vector<erhic::ParticleMC>& r)
  {
    for(const auto& p : r)
    {
      if(p.id!=11 || p.gs!=1 || p.orig<0) continue;
      const auto& m = r[p.orig];
      if(m.id==11 && m.gs==21) return p.I;
    }
    return -1;
  }

  static void SetKin(erhic::ParticleMC& p, const TLorentzVector& P, const TLorentzVector& q, bool haveQ)
  {
    TLorentzVector v(p.px,p.py,p.pz,p.E);
    p.eta = v.Eta();
    p.z = 0.; p.ptVsGamma = 0.;
    if(!haveQ) return;
    double Pq = P*q;
    if(Pq!=0.) p.z = (P*v)/Pq;
    TVector3 qh = q.Vect().Unit();
    TVector3 pv = v.Vect();
    p.ptVsGamma = (pv - qh*(pv.Dot(qh))).Mag();
  }
};

bool EpicReader::Load(long entry, erhic::EventPythia* ev)
{
  if(entry<0 || entry>=nEntries) return false;

  podio::Frame frame(reader.readEntry("events", entry));
  ev->Reset();
  ev->useReco = cfg.useReco;

  //---- one-time checks on the collections
  if(!checked)
  {
    checked = true;
    auto names = frame.getAvailableCollections();
    auto has = [&](const std::string& n){ return std::find(names.begin(),names.end(),n)!=names.end(); };
    if(!has(cfg.mcLabColl)) { cout<<"EpicReader ERROR: collection "<<cfg.mcLabColl<<" not found"<<endl; return false; }
    hasHO = has(cfg.mcHOColl);
    if(!hasHO) cout<<"EpicReader WARNING: "<<cfg.mcHOColl<<" not found, using "<<cfg.mcLabColl<<" (crossing angle!)"<<endl;
    if(!cfg.useReco && cfg.truthColl!=cfg.mcLabColl && !has(cfg.truthColl))
    { cout<<"EpicReader ERROR: "<<cfg.truthColl<<" not found"<<endl; return false; }
    if(cfg.useReco)
    {
      if(!has(cfg.recoColl))  { cout<<"EpicReader ERROR: "<<cfg.recoColl<<" not found"<<endl; return false; }
      hasAssoc = has(cfg.assocColl);
      if(!hasAssoc && cfg.pidMode==0) { cout<<"EpicReader ERROR: "<<cfg.assocColl<<" not found"<<endl; return false; }
      if(!hasAssoc) cout<<"EpicReader WARNING: no associations -> no mother info, cutdecay will reject all reco hadrons"<<endl;
      hasScatE = has(cfg.scatEColl);
      if(!hasScatE) cout<<"EpicReader ERROR: "<<cfg.scatEColl<<" not found, no reco electron!"<<endl;
    }
  }

  //---- Pythia information (Frame parameters)
  ev->process    = (int)ParD(frame,"process");
  ev->tgtparton  = (int)ParD(frame,"tgtparton");
  ev->xtgtparton = ParD(frame,"xtgtparton");
  ev->EBrems     = ParD(frame,"EBrems");
  ev->SigRadCor  = ParD(frame,"SigRadCor");

  //---- MC records
  const auto& mcLab = frame.get<edm4hep::MCParticleCollection>(cfg.mcLabColl);
  int nGen = (int)ParD(frame,"nTracks",0);     //number of generator particles (no Geant4 secondaries)
  if(nGen<=0 || nGen>(int)mcLab.size()) nGen = (int)mcLab.size();

  const edm4hep::MCParticleCollection& mcHO = hasHO ? frame.get<edm4hep::MCParticleCollection>(cfg.mcHOColl) : mcLab;
  FillRecord(mcHO, nGen, false, ev->ho);

  if(cfg.useReco) FillRecord(mcLab, nGen, false, ev->tracks);
  else
  {
    const auto& mcT = (cfg.truthColl==cfg.mcLabColl) ? mcLab : frame.get<edm4hep::MCParticleCollection>(cfg.truthColl);
    FillRecord(mcT, nGen, true, ev->tracks);
  }

  //sanity check (first events): lab and head-on records must be identical in structure
  if(nWarn<3 && hasHO && ev->ho.size()!=ev->tracks.size()) { cout<<"EpicReader WARNING: record sizes differ"<<endl; nWarn++; }
  for(size_t i=0; nWarn<3 && i<ev->ho.size() && i<ev->tracks.size(); ++i)
    if(ev->ho[i].id!=ev->tracks[i].id) { cout<<"EpicReader WARNING: MC record mismatch at index "<<i<<endl; nWarn++; break; }

  //---- beams and (true) scattered electron
  ev->iBeamE = FindBeam(ev->tracks, 11);
  ev->iBeamP = FindBeam(ev->tracks, 2212);
  ev->iScatE = FindScatE(ev->tracks);
  if(ev->iBeamE<0 || ev->iBeamP<0) { if(nWarn<10) { cout<<"EpicReader WARNING: beam particles not found, entry "<<entry<<endl; nWarn++; } return true; }

  const erhic::ParticleMC& be = ev->tracks[ev->iBeamE];
  const erhic::ParticleMC& bp = ev->tracks[ev->iBeamP];
  TLorentzVector k(be.px,be.py,be.pz,be.E);
  TLorentzVector P(bp.px,bp.py,bp.pz,bp.E);

  //true photon (for z and pt of the MC record)
  bool haveQT = (ev->iScatE>=0);
  TLorentzVector qT;
  if(haveQT)
  {
    const erhic::ParticleMC& se = ev->tracks[ev->iScatE];
    qT = k - TLorentzVector(se.px,se.py,se.pz,se.E);
  }
  for(auto& t : ev->tracks) SetKin(t, P, qT, haveQT);
  for(auto& t : ev->ho)     { TLorentzVector v(t.px,t.py,t.pz,t.E); t.eta = v.Eta(); }

  //---- DIS kinematics
  if(!cfg.useReco)
  {
    ev->QSquared = ParD(frame,"QSquared");
    ev->x        = ParD(frame,"x");
    ev->y        = ParD(frame,"y");
    ev->WSquared = ParD(frame,"WSquared");
    ev->nu       = ParD(frame,"nu");
    return true;
  }

  //---- RECO: scattered electron, electron-method kinematics
  ev->QSquared = -1.;
  if(hasScatE)
  {
    const auto& se = frame.get<edm4eic::ReconstructedParticleCollection>(cfg.scatEColl);
    if(se.size()>0)
    {
      auto e = se[0];
      auto mom = e.getMomentum();
      erhic::ParticleMC& r = ev->recoE;
      r = erhic::ParticleMC();
      r.I = -2; r.id = 11; r.KS = 0;
      r.px = mom.x; r.py = mom.y; r.pz = mom.z; r.E = e.getEnergy(); r.m = e.getMass();
      ev->hasRecoE = true;
    }
  }

  bool haveQR = false;
  TLorentzVector qR;
  if(ev->hasRecoE)
  {
    const erhic::ParticleMC& r = ev->recoE;
    qR = k - TLorentzVector(r.px,r.py,r.pz,r.E);
    double Pq = P*qR;
    if(Pq>0.)
    {
      haveQR = true;
      ev->QSquared = -qR.M2();
      ev->y        = Pq/(P*k);
      ev->x        = ev->QSquared/(2.*Pq);
      ev->WSquared = (P+qR).M2();
      ev->nu       = Pq/P.M();
    }
    SetKin(ev->recoE, P, qR, haveQR);
  }

  //---- RECO: charged hadron candidates appended to the record
  const auto& rc = frame.get<edm4eic::ReconstructedParticleCollection>(cfg.recoColl);
  std::vector<int> simIdx(rc.size(), -1);
  if(hasAssoc)
  {
    const auto& as = frame.get<edm4eic::MCRecoParticleAssociationCollection>(cfg.assocColl);
    std::vector<double> wbest(rc.size(), -1.);
    for(const auto& a : as)
    {
      int ir = a.getRec().getObjectID().index;
      if(ir<0 || ir>=(int)rc.size()) continue;
      if(a.getWeight()>wbest[ir]) { wbest[ir]=a.getWeight(); simIdx[ir]=a.getSim().getObjectID().index; }
    }
  }

  //walk from a Geant4 secondary up to its generator ancestor (index < nGen)
  auto walkUp = [&](int idx)->int
  {
    for(int s=0; s<30 && idx>=nGen; ++s)
    {
      if(idx<0 || idx>=(int)mcLab.size()) return -1;
      int par=-1;
      auto mp = mcLab[idx];
      for(const auto& q : mp.getParents()) { par=q.getObjectID().index; break; }
      idx = par;
    }
    return (idx>=0 && idx<nGen) ? idx : -1;
  };

  for(size_t ir=0; ir<rc.size(); ++ir)
  {
    auto r = rc[ir];
    auto mom = r.getMomentum();
    erhic::ParticleMC t;
    t.px = mom.x; t.py = mom.y; t.pz = mom.z; t.E = r.getEnergy(); t.m = r.getMass();

    //do not use the scattered electron as a hadron candidate
    if(ev->hasRecoE)
    {
      const erhic::ParticleMC& e = ev->recoE;
      double d = sqrt(pow(t.px-e.px,2)+pow(t.py-e.py,2)+pow(t.pz-e.pz,2));
      if(d < 1e-3*std::max(1.,e.E)) continue;
    }

    int gi = hasAssoc ? walkUp(simIdx[ir]) : -1;
    if(gi>=0) t.orig = ev->tracks[gi].orig;  //mother of the matched generator particle (for cutdecay)
    if(cfg.pidMode==0)
    {
      if(gi<0) continue;                     //no generator match -> no truth PID
      t.id = ev->tracks[gi].id;
    }
    else t.id = r.getPDG();

    t.I  = (int)ev->tracks.size();
    t.KS = 1;
    SetKin(t, P, qR, haveQR);
    ev->tracks.push_back(t);
  }

  return true;
}

//==========================================================================================

double Interpolate5D(THnSparseD* h, double Q2, double x, double z1, double z2, double phi)
{
  if (!h) return 0.0;
  if (h->GetNdimensions() != 5) return 0.0;

  double coords[5] = {Q2, x, z1, z2, phi};

  int    binLow[5];
  int    binHigh[5];
  double frac[5];

  // ------------------------------
  // Determine interpolation points
  // ------------------------------
  for (int d = 0; d < 5; ++d)
  {
    TAxis* ax = h->GetAxis(d);

    if (coords[d] < ax->GetXmin() || coords[d] >= ax->GetXmax())
      return 0.0;

    int b = ax->FindBin(coords[d]);

    // ==========================================================
    // Special treatment for x-axis (axis 1)
    // Use low edges as theory grid points
    // ==========================================================
    if (d == 1)
    {
      // Last grid point: cannot interpolate upward
      if (b >= ax->GetNbins())
      {
        binLow[d]  = ax->GetNbins();
        binHigh[d] = ax->GetNbins();
        frac[d]    = 0.0;
        continue;
      }

      double x1 = ax->GetBinLowEdge(b);
      double x2 = ax->GetBinLowEdge(b + 1);

      binLow[d]  = b;
      binHigh[d] = b + 1;

      frac[d] = (coords[d] - x1) / (x2 - x1);

      // numerical safety
      if (frac[d] < 0.0) frac[d] = 0.0;
      if (frac[d] > 1.0) frac[d] = 1.0;

      continue;
    }

    // ==========================================================
    // ROOT-style center interpolation for all other axes
    // ==========================================================
    double center = ax->GetBinCenter(b);

    int b2;
    if (coords[d] >= center)
      b2 = b + 1;
    else
      b2 = b - 1;

    if (b2 < 1 || b2 > ax->GetNbins())
    {
      binLow[d]  = b;
      binHigh[d] = b;
      frac[d]    = 0.0;
      continue;
    }

    double center2 = ax->GetBinCenter(b2);

    if (std::abs(center2 - center) < 1e-12)
    {
      binLow[d]  = b;
      binHigh[d] = b;
      frac[d]    = 0.0;
      continue;
    }

    binLow[d]  = b;
    binHigh[d] = b2;

    frac[d] = (coords[d] - center) / (center2 - center);
  }

  // ------------------------------
  // 5D multilinear interpolation
  // ------------------------------
  double value = 0.0;

  for (int i0 = 0; i0 <= 1; ++i0)
  for (int i1 = 0; i1 <= 1; ++i1)
  for (int i2 = 0; i2 <= 1; ++i2)
  for (int i3 = 0; i3 <= 1; ++i3)
  for (int i4 = 0; i4 <= 1; ++i4)
  {
    int bins[5] =
    {
      (i0 ? binHigh[0] : binLow[0]),
      (i1 ? binHigh[1] : binLow[1]),
      (i2 ? binHigh[2] : binLow[2]),
      (i3 ? binHigh[3] : binLow[3]),
      (i4 ? binHigh[4] : binLow[4])
    };

    double weight =
      (i0 ? frac[0] : 1.0 - frac[0]) *
      (i1 ? frac[1] : 1.0 - frac[1]) *
      (i2 ? frac[2] : 1.0 - frac[2]) *
      (i3 ? frac[3] : 1.0 - frac[3]) *
      (i4 ? frac[4] : 1.0 - frac[4]);

    value += weight * h->GetBinContent(bins);
  }

  return value;
}

//==========================================================================================
//Saturation correction: ratio of the theory correlation functions eA/ep at the pair kinematics.
//Interpolate5D returns 0 outside the Cth grid -> weight 1 there.
double SaturationCorrection(double Q2, double x, double z1, double z2, double phi, THnSparseD* Cth_eA, THnSparseD* Cth_ep)
{
  double C_p = Interpolate5D(Cth_ep,Q2,x,z1,z2,phi);
  if(C_p==0.) return 1.;
  return Interpolate5D(Cth_eA,Q2,x,z1,z2,phi)/C_p;
}

//Pythia6 process classes used for the per-process Delta-phi histograms
enum { kPGF, kQCDC, kLODIS, kResQQ, kResQG, kResGQ, kResGG, nProc };
const char* procName[nProc]  = {"PGF","QCDC","LODIS","ResQQ","ResQG","ResGQ","ResGG"};
const bool  procGluon[nProc] = {true, false, false,  false,  true,   false,  true};  //get the saturation weight

//returns -1 for processes not in any class (diffractive, ...)
int ProcessClass(int process, int tgtparton)
{
  if(process>132)  return kPGF;
  if(process>130)  return kQCDC;
  if(process==99)  return kLODIS;
  if(process>=90)  return -1;
  if(process<20)   return kResQQ;
  if(process==28)  return (tgtparton==21) ? kResQG : kResGQ;
  return kResGG;
}

//direct daughter of a Pythia string (id 92), i.e. not a resonance-decay product
bool FromString(erhic::EventPythia* event, const erhic::ParticleMC* p)
{
  int o = p->orig;
  return o>=0 && o<event->GetNTracks() && event->GetTrack(o)->id==92;
}

std::vector<double> LogBins(int n, double lo, double hi)
{
  std::vector<double> e(n+1);
  for(int i=0; i<=n; i++) e[i] = lo*pow(hi/lo, double(i)/n);
  return e;
}

//raw counters, summed by hadd; needed to normalise after merging
void WriteCounters(long Ntrig, long NtrigEvt, long Nevt, long nRead)
{
  TH1D h("hCounters","Ntrig, NtrigEvt, Nevt, nEventsRead;;counts",4,0,4);
  h.SetDirectory(0);
  const char* label[4] = {"Ntrig","NtrigEvt","Nevt","nEventsRead"};
  double      value[4] = {double(Ntrig),double(NtrigEvt),double(Nevt),double(nRead)};
  for(int i=0; i<4; i++)
  {
    h.GetXaxis()->SetBinLabel(i+1,label[i]);
    h.SetBinContent(i+1,value[i]);
  }
  h.Write();
}

//histograms and counters of one x_Bj decade
//in memory the names carry a "_d<i>" suffix (unique per decade); it is stripped when writing
struct DecadeHists
{
  TString dir;
  double  xLow=0, xHigh=0;
  long    Nevt=0, Ntrig=0, NtrigEvt=0;

  TH1D *procTot, *deltaPhi;
  TH1D *deltaPhiProc[nProc];
  TH1D *deltaPhiProcSat[nProc];   //NULL for the quark-initiated channels
  TH1D *TrigPt, *TrigZ, *AssoPt, *AssoZ;
  TH2D *TrigPtVsEta, *AssoPtVsEta, *zCorre;

  void Book(int i, double lo, double hi)
  {
    const double pi = TMath::Pi();
    xLow = lo; xHigh = hi;
    dir  = Form("x_%.0e_%.0e",lo,hi);
    dir.ReplaceAll("e-0","e-").ReplaceAll("e+00","");   //x_1e-4_1e-3 ... x_1e-1_1
    TString s   = Form("_d%d",i);
    TString tag = Form(" (%.0e #leq x < %.0e)",lo,hi);

    procTot  = new TH1D("procTot"+s, "process"+tag+";id;events",140,0,140);
    deltaPhi = new TH1D("deltaPhi"+s,"#Delta#phi"+tag+";#Delta#phi;pairs",60,-pi/2,3*pi/2);
    for(int p=0; p<nProc; p++)
    {
      deltaPhiProc[p] = new TH1D(Form("deltaPhi%s",procName[p])+s, Form("#Delta#phi %s",procName[p])+tag+";#Delta#phi;pairs",60,-pi/2,3*pi/2);
      deltaPhiProcSat[p] = !procGluon[p] ? NULL :
        new TH1D(Form("deltaPhi%sSat",procName[p])+s, Form("#Delta#phi %s, saturation weighted",procName[p])+tag+";#Delta#phi;pairs",60,-pi/2,3*pi/2);
    }

    TrigPt      = new TH1D("TrigPt"+s,     "trigger"+tag+";p_{T} [GeV];triggers",100,0,20);
    TrigZ       = new TH1D("TrigZ"+s,      "trigger"+tag+";z;triggers",100,0,1);
    TrigPtVsEta = new TH2D("TrigPtVsEta"+s,"trigger"+tag+";#eta_{lab};p_{T} [GeV]",100,-5,5,100,0,20);
    AssoPt      = new TH1D("AssoPt"+s,     "associate"+tag+";p_{T} [GeV];pairs",100,0,20);
    AssoZ       = new TH1D("AssoZ"+s,      "associate"+tag+";z;pairs",100,0,1);
    AssoPtVsEta = new TH2D("AssoPtVsEta"+s,"associate"+tag+";#eta_{lab};p_{T} [GeV]",100,-5,5,100,0,20);
    zCorre      = new TH2D("zCorre"+s,     "pairs"+tag+";z^{trig};z^{asso}",100,0,1,100,0,1);
  }

  static void WriteH(TH1* h)
  {
    TString n = h->GetName();
    n.Remove(n.Last('_'));
    h->Write(n);
  }

  void Write(long nRead)
  {
    WriteCounters(Ntrig,NtrigEvt,Nevt,nRead);
    WriteH(procTot);
    WriteH(deltaPhi);
    for(int p=0; p<nProc; p++)
    {
      WriteH(deltaPhiProc[p]);
      if(deltaPhiProcSat[p]) WriteH(deltaPhiProcSat[p]);
    }
    WriteH(TrigPt); WriteH(TrigZ); WriteH(TrigPtVsEta);
    WriteH(AssoPt); WriteH(AssoZ); WriteH(AssoPtVsEta);
    WriteH(zCorre);
  }
};

//=================================================================================================
//TrigCutLow/High, AssoCutLow/High: pT window of trigger / associate (gamma*-p frame)
//filelist: list of podio files (local or root://); outFile: output root file
//sys: ep, eAu, eAg -> selects ../TH/Cth_<sys>.root for the saturation weight
//doSubCut: keep only gluon-initiated events (tgtparton==21, no diffraction)
//cutdecay: keep only hadrons coming directly from the string (reject resonance decays)
//useReco, pidMode, recoPartColl: see the EpicReader description at the top
int ePIC_multikt2(double TrigCutLow, double TrigCutHigh, double AssoCutLow, double AssoCutHigh,
                  TString filelist, TString outFile, TString sys="ep",
                  bool doSubCut=false, bool cutdecay=false,
                  bool useReco=false, int pidMode=0,
                  TString recoPartColl="ReconstructedChargedParticles")
{
  const double pi = TMath::Pi();

  //---- cuts
  const double Q2Min = 1., Q2Max = 9.;   //Cth Q2 grid
  const double AcceptCut = 4.0;          //|eta_lab| of the scattered electron and of the hadrons
  //Cth is integrated over z1+z2<1, |z1-z2|<=0.1, which implies 0.1 <= z_trig <= 0.525
  const double zTrigLow = 0.1, zTrigHigh = 0.525;
  const double zSumMax  = 1.,  zDiffMax  = 0.1;
  const double zMin     = 0.3;           //on trigger and associate: low z dilutes the saturation signal

  //x_Bj decades, each with its own set of histograms
  const int    nDec = 4;
  const double decEdges[nDec+1] = {1e-4,1e-3,1e-2,1e-1,1.};

  //Cth grid binning (QA of the accepted events)
  double Q2_bins[10] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};
  double x_bins[37]  = {
    1E-4,2E-4,3E-4,4E-4,5E-4,6E-4,7E-4,8E-4,9E-4,
    1E-3,2E-3,3E-3,4E-3,5E-3,6E-3,7E-3,8E-3,9E-3,
    1E-2,2E-2,3E-2,4E-2,5E-2,6E-2,7E-2,8E-2,9E-2,
    1E-1,2E-1,3E-1,4E-1,5E-1,6E-1,7E-1,8E-1,9E-1,
    1.0
  };

  //---- theory correlation functions for the saturation weight
  //run_job.sh sets CTH_DIR to TH/ in this directory, a copy of /eic/u/gabrielg/Saturation_vES/TH/Cth_*.root
  const char* cthEnv = gSystem->Getenv("CTH_DIR");
  const TString thDir = cthEnv ? cthEnv : "TH";
  TFile* Fth_eA = TFile::Open(thDir+"/Cth_"+sys+".root");
  TFile* Fth_ep = TFile::Open(thDir+"/Cth_ep.root");
  THnSparseD* Cth_eA = Fth_eA ? (THnSparseD*)Fth_eA->Get("Cth_"+sys) : NULL;
  THnSparseD* Cth_ep = Fth_ep ? (THnSparseD*)Fth_ep->Get("Cth_ep")   : NULL;
  if(!Cth_eA || !Cth_ep)
  {
    cout<<"cannot load Cth_"<<sys<<" / Cth_ep from "<<thDir<<endl;
    return 1;
  }

  //---- input files
  ifstream fileNames(filelist);
  if(!(fileNames))
  {
    cout<<"cannot open file list!"<<endl;
    return 1;
  }

  char rootfile[1024];
  TObjString inputFileList;
  std::vector<std::string> inputFiles;

  while(fileNames.getline(rootfile, 1024))
  {
    if(strlen(rootfile)==0 || rootfile[0]=='#') continue;

    if( !TString(rootfile).BeginsWith("root://") )
    {
      TFile *fTmp = TFile::Open(rootfile);
      if( !fTmp || fTmp->IsZombie() || !(fTmp->GetNkeys()) )
      {
        cout<<"read in root file open error: "<<rootfile<<endl;
        return 1;
      }
      fTmp->Close();
      delete fTmp;
    }

    inputFileList.String()+=TString(rootfile);
    inputFileList.String()+="\n";
    inputFiles.push_back(rootfile);
    cout<<"input file:"<<rootfile<<endl;
  }

  if(inputFiles.empty())
  {
    cout<<"empty name in filelist!"<<endl;
    return 1;
  }
  gROOT->cd();   //book the histograms in memory, not in the Cth/input files

  //podio reader configuration
  EpicConfig epicCfg;
  epicCfg.useReco  = useReco;
  epicCfg.pidMode  = pidMode;
  epicCfg.recoColl = recoPartColl.Data();
  if(useReco) cout<<"RECO mode: hadrons from "<<epicCfg.recoColl<<", pidMode="<<pidMode<<endl;
  else        cout<<"TRUTH mode: generator record "<<epicCfg.truthColl<<endl;

  EpicReader reader(inputFiles, epicCfg);
  long sampleSize = reader.GetEntries();
  if( !(sampleSize>0) )
  {
    cout<<"no events found in the input files"<<endl;
    return 1;
  }
  cout<<"Total event sample size: "<<sampleSize<<endl;
  gROOT->cd();

  TObjString cuts;
  cuts.String().Form("Pt trig in (%.2f,%.2f), Pt asso in (%.2f,%.2f), Q2 in [%.1f,%.1f], |eta|<%.1f, "
                     "z trig in [%.3f,%.3f], z>%.2f, z1+z2<%.1f, |z1-z2|<=%.2f",
                     TrigCutLow,TrigCutHigh,AssoCutLow,AssoCutHigh,Q2Min,Q2Max,AcceptCut,
                     zTrigLow,zTrigHigh,zMin,zSumMax,zDiffMax);
  TObjString config;
  config.String().Form("sys=%s useReco=%d pidMode=%d recoColl=%s doSubCut=%d cutdecay=%d",
                       sys.Data(),useReco,pidMode,recoPartColl.Data(),doSubCut,cutdecay);
  TObjString evtNumTot;
  evtNumTot.String().Form("%ld",sampleSize);
  cout<<cuts.String()<<endl<<config.String()<<endl;
  cout<<"output file will be:"<<outFile<<endl;

  //-------------book histograms---------------------------------
  TH1::SetDefaultSumw2();

  //all events, before cuts
  std::vector<double> xLog  = LogBins(100,1e-5,1.);
  std::vector<double> Q2Log = LogBins(100,0.5,100.);
  TH2D* hQ2_x  = new TH2D("Q2_x", "all events;x;Q^{2} [GeV^{2}]",    100,xLog.data(),100,Q2Log.data());
  TH2D* hQ2_xg = new TH2D("Q2_xg","all events;x_{g};Q^{2} [GeV^{2}]",100,xLog.data(),100,Q2Log.data());

  //accepted events (Q2, scattered electron, x in [1e-4,1))
  TH2D* Q2vx   = new TH2D("Q2vx",  "accepted events;x;Q^{2} [GeV^{2}]",36,x_bins,9,Q2_bins);
  TH1D* x_hist = new TH1D("x_hist","accepted events;x;events",36,x_bins);
  TH1D* y_hist = new TH1D("y_hist","accepted events;y;events",100,0,1);

  DecadeHists dec[nDec];
  for(int i=0; i<nDec; i++) dec[i].Book(i,decEdges[i],decEdges[i+1]);

  //-------------event loop---------------------------------------
  //hadron candidate: z and 4-vector in the gamma*-p CM frame
  struct Hadron { double z, etaLab; TLorentzVector v; };
  std::vector<Hadron> had;
  std::vector<int>    trig;

  erhic::EventPythia* event = new erhic::EventPythia;
  long entry = 0;

  TStopwatch timer;
  timer.Start();

  while( reader.Load(entry, event) )
  {
    entry++;
    if( entry%100000 == 0 ) cout<<entry<<" done"<<endl;

    const double Q2      = event->QSquared;
    const double x       = event->x;
    const int    process = event->GetProcess();

    hQ2_x ->Fill(x,Q2);
    hQ2_xg->Fill(event->xtgtparton,Q2);

    //---- event selection
    if(Q2<Q2Min || Q2>Q2Max) continue;

    //scattered electron: MC (truth mode) or ScatteredElectronsEMinusPz (reco mode)
    erhic::ParticleMC* electron_out = event->GetScatteredElectron();
    erhic::ParticleMC* electron_in  = event->GetBeamElectron();
    erhic::ParticleMC* proton       = event->GetBeamProton();
    if(electron_out==NULL || electron_in==NULL || proton==NULL) continue;
    if(fabs(electron_out->eta)>AcceptCut) continue;

    int d = -1;
    for(int i=0; i<nDec; i++) if(x>=decEdges[i] && x<decEdges[i+1]) { d=i; break; }
    if(d<0) continue;
    DecadeHists& D = dec[d];

    x_hist->Fill(x);
    Q2vx->Fill(x,Q2);
    y_hist->Fill(event->y);
    D.procTot->Fill(process);

    if(doSubCut && (event->tgtparton!=21 || (process<100 && process>90))) continue;
    D.Nevt++;

    const int pc = ProcessClass(process, event->tgtparton);

    //virtual photon; boost to the gamma*-p CM frame and rotate the photon onto +z
    TLorentzVector vP(proton->px,proton->py,proton->pz,proton->E);
    TLorentzVector q(electron_in->px-electron_out->px, electron_in->py-electron_out->py,
                     electron_in->pz-electron_out->pz, electron_in->E -electron_out->E);
    TVector3 b = -(vP+q).BoostVector();
    TLorentzVector qCM = q;
    qCM.Boost(b);

    //---- charged pi/K/p candidates
    had.clear();
    for(int j=0; j<event->GetNTracks(); j++)
    {
      erhic::ParticleMC* p = event->GetTrack(j);
      if(p->KS!=1) continue;
      int aid = abs(p->id);
      if(aid!=211 && aid!=321 && aid!=2212) continue;
      if(fabs(p->eta)>AcceptCut) continue;
      if(cutdecay && !FromString(event,p)) continue;

      Hadron h;
      h.v.SetPxPyPzE(p->px,p->py,p->pz,p->E);
      h.z      = (h.v*vP)/(vP*q);   //lab frame
      h.etaLab = p->eta;
      h.v.Boost(b);
      h.v.RotateZ(-qCM.Phi());
      h.v.RotateY(-qCM.Theta());
      had.push_back(h);
    }

    //---- triggers, in decreasing pT
    trig.clear();
    for(int i=0; i<(int)had.size(); i++)
    {
      const Hadron& h = had[i];
      double pt = h.v.Pt();
      if(pt<TrigCutLow || pt>TrigCutHigh) continue;
      if(h.z<zTrigLow || h.z>zTrigHigh || h.z<zMin) continue;
      trig.push_back(i);
    }
    std::sort(trig.begin(), trig.end(), [&](int a, int c){ return had[a].v.Pt() > had[c].v.Pt(); });

    //multiple-trigger algorithm: every trigger is paired with all associates of lower pT
    int    eventTrigNum = 0;
    double lastPt = -1;
    for(int it : trig)
    {
      const Hadron& t = had[it];
      const double ptT = t.v.Pt();
      if(eventTrigNum>0 && ptT > lastPt-0.00001) continue;   //same pT as the previous trigger
      lastPt = ptT;
      eventTrigNum++;
      D.Ntrig++;

      D.TrigPt->Fill(ptT);
      D.TrigZ->Fill(t.z);
      D.TrigPtVsEta->Fill(t.etaLab, ptT);

      for(int ia=0; ia<(int)had.size(); ia++)
      {
        if(ia==it) continue;
        const Hadron& a = had[ia];
        const double ptA = a.v.Pt();
        if(ptA<AssoCutLow || ptA>AssoCutHigh || ptA>ptT) continue;
        if(t.z + a.z >= zSumMax) continue;
        if(fabs(t.z - a.z) > zDiffMax) continue;
        if(a.z<zMin) continue;

        double deltaphi = t.v.DeltaPhi(a.v);
        if(deltaphi < -pi/2) deltaphi += 2*pi;

        D.deltaPhi->Fill(deltaphi);
        if(pc>=0)
        {
          D.deltaPhiProc[pc]->Fill(deltaphi);
          if(procGluon[pc]) D.deltaPhiProcSat[pc]->Fill(deltaphi, SaturationCorrection(Q2,x,t.z,a.z,deltaphi,Cth_eA,Cth_ep));
        }

        D.AssoPt->Fill(ptA);
        D.AssoZ->Fill(a.z);
        D.AssoPtVsEta->Fill(a.etaLab, ptA);
        D.zCorre->Fill(t.z, a.z);
      }
    }
    if(eventTrigNum>0) D.NtrigEvt++;

  }//event loop

  timer.Stop();
  if(entry!=sampleSize)   //Load() also returns false on a read problem
  {
    cout<<"event loop stopped early: "<<entry<<" of "<<sampleSize<<" entries read"<<endl;
    return 1;
  }
  cout<<"total run time "<<timer.RealTime()<<" sec, "<<timer.RealTime()/float(entry)<<" sec per event"<<endl;

  //-------------summary and output-------------------------------
  long Nevt=0, Ntrig=0, NtrigEvt=0;
  for(int i=0; i<nDec; i++)
  {
    const DecadeHists& D = dec[i];
    Nevt += D.Nevt; Ntrig += D.Ntrig; NtrigEvt += D.NtrigEvt;
    cout<<D.dir<<": accepted events="<<D.Nevt<<" triggers="<<D.Ntrig<<" trigger events="<<D.NtrigEvt
        <<" pairs="<<D.deltaPhi->GetEntries()<<endl;
  }
  cout<<"total: events read="<<entry<<" accepted="<<Nevt<<" triggers="<<Ntrig<<" trigger events="<<NtrigEvt<<endl;

  TFile fout(outFile,"RECREATE");
  inputFileList.Write("inputFileList");
  cuts.Write("cuts");
  config.Write("config");
  evtNumTot.Write("nEvents");
  WriteCounters(Ntrig,NtrigEvt,Nevt,sampleSize);
  hQ2_x->Write();
  hQ2_xg->Write();
  Q2vx->Write();
  x_hist->Write();
  y_hist->Write();

  for(int i=0; i<nDec; i++)
  {
    fout.mkdir(dec[i].dir)->cd();
    dec[i].Write(sampleSize);
  }

  fout.Close();
  return 0;
}
