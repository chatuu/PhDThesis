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

const Var calE([](const caf::SRProxy *sr){
if(sr->vtx.nelastic==0)
  return -5.0;

return double(sr->slc.calE);
});

int GetMuonProngId(const caf::SRProxy *sr){
  int prongNum=0;
  float maxVal=-1.0;

  if(sr->vtx.nelastic==0 || sr->vtx.elastic[0].fuzzyk.npng==0)
    return -5;

  for(size_t i=0;i<sr->vtx.elastic[0].fuzzyk.npng;++i){
    if(sr->vtx.elastic[0].fuzzyk.png[i].cvnpart.muonid>maxVal){
      maxVal=sr->vtx.elastic[0].fuzzyk.png[i].cvnpart.muonid;
      prongNum=i;
    }
    else if(sr->vtx.elastic[0].fuzzyk.png[i].len>=500){
      maxVal=1.0;
      prongNum=i;
    }
  }
  return prongNum;
}

const Var muonCalE([](const caf::SRProxy *sr){
  int prongNum=GetMuonProngId(sr);
  if(prongNum<0)
    return -5.0;
  else  
    return double(sr->vtx.elastic[0].fuzzyk.png[prongNum].calE);
});

const Var pionCalE([](const caf::SRProxy *sr){
  int muonNum=GetMuonProngId(sr);
  if(muonNum<0)
    return -5.0;
  else{
    int prongNum=0;
    int maxVal=-1;
    if(sr->vtx.nelastic==0 || sr->vtx.elastic[0].fuzzyk.npng==0)
      return -5.0;

    for(size_t i=0;i<sr->vtx.elastic[0].fuzzyk.npng;++i){
      if(sr->vtx.elastic[0].fuzzyk.png[i].cvnpart.pionid>maxVal &&
        i != muonNum){
        maxVal=sr->vtx.elastic[0].fuzzyk.png[i].cvnpart.pionid;
        prongNum=i;
      }
    }
    return double(sr->vtx.elastic[0].fuzzyk.png[prongNum].calE);
  }
  //if(sr->vtx.elastic[0].fuzzyk.png[prongNum].truth.pdg == 211)
  //else if (! ())
});

const Var vertexAct=calE-muonCalE-pionCalE;



