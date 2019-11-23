#include "headers.h"

const NuTruthVar nutrinoEnergy([](const caf::SRNeutrinoProxy *nu) {
  return nu->E;
});

const NuTruthVar X([](const caf::SRNeutrinoProxy *nu) {
  return nu->x;
});

const NuTruthVar Y([](const caf::SRNeutrinoProxy *nu) {
  return nu->y;
});

const NuTruthVar pdgValue([](const caf::SRNeutrinoProxy *nu) {
  return nu->pdg;
});

const Var kMuonEnergy([](const caf::SRProxy *sr) {
  if (sr->mc.nu.size() == 0)
    return -5.0;

  double MuonEnergy = -5;

  for (unsigned int i = 0; i < sr->mc.nu[0].prim.size(); i++)
  {

    if (sr->mc.nu[0].prim[i].pdg == 13)
      MuonEnergy = double(sr->mc.nu[0].prim[i].p.E);
  }

  return MuonEnergy;
});

const Var kPionEnergy([](const caf::SRProxy *sr) {
  if (sr->mc.nu.size() == 0)
    return -5.0;
  double PionEnergy;

  for (unsigned int i = 0; i < sr->mc.nu[0].prim.size(); i++)
  {

    if (sr->mc.nu[0].prim[i].pdg == 211)
      PionEnergy = double(sr->mc.nu[0].prim[i].p.E);
  }

  return PionEnergy;
});


const Var kt([](const caf::SRProxy *sr) {
  if (sr->mc.nu.size() == 0)
    return -5.0;

  TLorentzVector Ppi(0,0,0,0);
  TLorentzVector Pmu(0,0,0,0);
  TLorentzVector Pnu(0,0,0,0);
  TLorentzVector Pf(0,0,0,0);

  Pnu.SetPxPyPzE(sr->mc.nu[0].p.px,sr->mc.nu[0].p.py,sr->mc.nu[0].p.pz,sr->mc.nu[0].p.E);

  bool hasMuon=false;
  bool hasPion=false;

  for (unsigned int i = 0; i < sr->mc.nu[0].prim.size(); i++)
  {
    if (sr->mc.nu[0].prim[i].pdg == 13){
      hasMuon=true;
      Pmu.SetPxPyPzE(sr->mc.nu[0].prim[i].p.px,sr->mc.nu[0].prim[i].p.py,sr->mc.nu[0].prim[i].p.pz,sr->mc.nu[0].prim[i].p.E);
    }
    
    if (sr->mc.nu[0].prim[i].pdg == 211){
      hasPion=true;
      Ppi.SetPxPyPzE(sr->mc.nu[0].prim[i].p.px,sr->mc.nu[0].prim[i].p.py,sr->mc.nu[0].prim[i].p.pz,sr->mc.nu[0].prim[i].p.E);
    }
  }

  if(!hasPion) return -5.0;
  if(!hasMuon) return -5.0;

  Pf=Pnu-Pmu-Ppi;

  return TMath::Abs(Pf*Pf);
});

const Var kPionAngle([](const caf::SRProxy *sr) {
  if (sr->mc.nu.size() == 0)
    return -5.0;

  double CosThetaPi=0;

  TVector3 dir(0,0,0);
  TVector3 beamdir = NuMIBeamDirection(caf::kNEARDET);

  bool hasPion = false;
  for (unsigned int i = 0; i < sr->mc.nu[0].prim.size(); i++)
  {

    if (sr->mc.nu[0].prim[i].pdg == 211)
    {
      hasPion = true;
      dir = sr->mc.nu[0].prim[i].p.Vect();
//      beamdir = NuMIBeamDirection(caf::kNEARDET);
      CosThetaPi = double(dir.Unit().Dot(beamdir));
    }
    

  }

    if(!hasPion) return -5.0;

  return (TMath::ACos(CosThetaPi))/(TMath::Pi());
});

const Var kMuonAngle([](const caf::SRProxy *sr) {
  if (sr->mc.nu.size() == 0)
    return -5.0;

  double CosThetaMu=0;

  TVector3 dir(0,0,0);
  TVector3 beamdir=NuMIBeamDirection(caf::kNEARDET);
  bool hasMuon=false;
  for (unsigned int i = 0; i < sr->mc.nu[0].prim.size(); i++)
  {

    if (sr->mc.nu[0].prim[i].pdg == 13)
    {
      hasMuon=true;
      dir = sr->mc.nu[0].prim[i].p.Vect();
      //beamdir = NuMIBeamDirection(caf::kNEARDET);
      CosThetaMu = double(dir.Unit().Dot(beamdir));
    }
  }

  if(!hasMuon) return -5.0;
  return (TMath::ACos(CosThetaMu))/(TMath::Pi());
});

const Var kZeta([](const caf::SRProxy *sr) {
  if (sr->mc.nu.size() == 0)
    return -5.0;
  double Zeta=0;
  TVector3 dir(0,0,0);
  TVector3 beamdir=NuMIBeamDirection(caf::kNEARDET);
  bool hasPion=false;
  for (unsigned int i = 0; i < sr->mc.nu[0].prim.size(); i++)
  {

    if (sr->mc.nu[0].prim[i].pdg == 211)
    {
      hasPion=true;
      dir = sr->mc.nu[0].prim[i].p.Vect();
      //beamdir = NuMIBeamDirection(caf::kNEARDET);
      Zeta = double((sr->mc.nu[0].prim[i].p.E) * (1 - dir.Unit().Dot(beamdir)));
    }
  }
  if(!hasPion) return -5.0;
  //	std::cout<<"The value is: "<<dir.Unit().Dot(beamdir)<<std::endl;
  return Zeta;
});

const Var NumberOfTracks([](const caf::SRProxy *sr){
  if(sr->vtx.nelastic==0)
    return -5.0;
  return double(sr->trk.kalman.ntracks);
});

const Var NumberOfProngs([](const caf::SRProxy *sr){
  if(sr->vtx.nelastic==0)
    return -5.0;

  return double(sr->vtx.elastic[0].fuzzyk.npng);
});


