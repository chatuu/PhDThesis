#include "headers.h"
#include <stdio.h>
//#include "cuts.h"

// Truth Variables
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

  return (MuonEnergy - 0.105658);
});

const Var kPionEnergy([](const caf::SRProxy *sr) {
  if (sr->mc.nu.size() == 0)
    return -5.0;
  double PionEnergy = -10.0;

  for (unsigned int i = 0; i < sr->mc.nu[0].prim.size(); i++)
  {

    if (sr->mc.nu[0].prim[i].pdg == 211)
      PionEnergy = double(sr->mc.nu[0].prim[i].p.E);
  }

  return (PionEnergy - 0.13957);
});

const Var kt([](const caf::SRProxy *sr) {
  if (sr->mc.nu.size() == 0)
    return -5.0;

  TLorentzVector Ppi(0, 0, 0, 0);
  TLorentzVector Pmu(0, 0, 0, 0);
  TLorentzVector Pnu(0, 0, 0, 0);
  TLorentzVector Pf(0, 0, 0, 0);

  Pnu.SetPxPyPzE(sr->mc.nu[0].p.px, sr->mc.nu[0].p.py, sr->mc.nu[0].p.pz, sr->mc.nu[0].p.E);

  bool hasMuon = false;
  bool hasPion = false;

  for (unsigned int i = 0; i < sr->mc.nu[0].prim.size(); i++)
  {
    if (sr->mc.nu[0].prim[i].pdg == 13)
    {
      hasMuon = true;
      Pmu.SetPxPyPzE(sr->mc.nu[0].prim[i].p.px, sr->mc.nu[0].prim[i].p.py, sr->mc.nu[0].prim[i].p.pz, sr->mc.nu[0].prim[i].p.E);
    }

    if (sr->mc.nu[0].prim[i].pdg == 211)
    {
      hasPion = true;
      Ppi.SetPxPyPzE(sr->mc.nu[0].prim[i].p.px, sr->mc.nu[0].prim[i].p.py, sr->mc.nu[0].prim[i].p.pz, sr->mc.nu[0].prim[i].p.E);
    }
  }

  if (!hasPion)
    return -5.0;
  if (!hasMuon)
    return -5.0;

  Pf = Pnu - Pmu - Ppi;

  return TMath::Abs(Pf * Pf);
});

const Var kPionAngle([](const caf::SRProxy *sr) {
  if (sr->mc.nu.size() == 0)
    return -5.0;

  double CosThetaPi = 0;

  TVector3 dir(0, 0, 0);
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

  if (!hasPion)
    return -5.0;

  return (TMath::ACos(CosThetaPi)) / (TMath::Pi());
});

const Var kMuonAngle([](const caf::SRProxy *sr) {
  if (sr->mc.nu.size() == 0)
    return -5.0;

  double CosThetaMu = 0;

  TVector3 dir(0, 0, 0);
  TVector3 beamdir = NuMIBeamDirection(caf::kNEARDET);
  bool hasMuon = false;
  for (unsigned int i = 0; i < sr->mc.nu[0].prim.size(); i++)
  {

    if (sr->mc.nu[0].prim[i].pdg == 13)
    {
      hasMuon = true;
      dir = sr->mc.nu[0].prim[i].p.Vect();
      //beamdir = NuMIBeamDirection(caf::kNEARDET);
      CosThetaMu = double(dir.Unit().Dot(beamdir));
    }
  }

  if (!hasMuon)
    return -5.0;
  return (TMath::ACos(CosThetaMu)) / (TMath::Pi());
});

const Var kZeta([](const caf::SRProxy *sr) {
  if (sr->mc.nu.size() == 0)
    return -5.0;
  double Zeta = 0;
  TVector3 dir(0, 0, 0);
  TVector3 beamdir = NuMIBeamDirection(caf::kNEARDET);
  bool hasPion = false;
  for (unsigned int i = 0; i < sr->mc.nu[0].prim.size(); i++)
  {

    if (sr->mc.nu[0].prim[i].pdg == 211)
    {
      hasPion = true;
      dir = sr->mc.nu[0].prim[i].p.Vect();
      //beamdir = NuMIBeamDirection(caf::kNEARDET);
      Zeta = double((sr->mc.nu[0].prim[i].p.E) * (1 - dir.Unit().Dot(beamdir)));
    }
  }
  if (!hasPion)
    return -5.0;
  //	std::cout<<"The value is: "<<dir.Unit().Dot(beamdir)<<std::endl;
  return Zeta;
});

const Var NumberOfTracks([](const caf::SRProxy *sr) {
  if (sr->vtx.elastic.IsValid != true)
    return -5.0;
  return double(sr->trk.kalman.ntracks);
});

const Var NumberOfProngs([](const caf::SRProxy *sr) {
  if (sr->vtx.elastic.IsValid != true)
    return -5.0;

  //std::cout << sr->vtx.elastic.fuzzyk.npng << "\n";
  return double(sr->vtx.elastic.fuzzyk.npng);
});

const Var calE([](const caf::SRProxy *sr) {
  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng != 2 || sr->vtx.elastic.fuzzyk.npng2d != 0)
    return -10.0;

  return double(sr->slc.calE);
});

/****************************************** Muon Track Number ********************************************************************/

// unsigned int GetMuonTrackId(const caf::SRProxy *sr)
// {
//   int prongNum = 0;
//   float maxVal = -1.0;

//   if (sr->trk.kalman.tracks. != true)
//     return 10000;

//   for (size_t i = 0; i < sr->trk.kalman.ntracks; ++i)
//   {
//     if (sr->vtx.elastic.fuzzyk.png[i].cvnpart.muonid > maxVal)
//     {
//       maxVal = sr->vtx.elastic.fuzzyk.png[i].cvnpart.muonid;
//       prongNum = i;
//     }
//     else if (sr->vtx.elastic.fuzzyk.png[i].len >= 500)
//     {
//       maxVal = 1.0;
//       prongNum = i;
//     }
//   }
//   return prongNum;
// }

/*********************************************************************************************************************************/
/****************************************** Muon Prong Number ********************************************************************/
/*
unsigned int GetMuonProngId(const caf::SRProxy *sr)
{
  int prongNum = 0;
  float maxVal = -1.0;

  if (sr->vtx.elastic.IsValid != true)
    return 10000;

  for (size_t i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
  {
    if (sr->vtx.elastic.fuzzyk.png[i].cvnpart.muonid > maxVal)
    {
      maxVal = sr->vtx.elastic.fuzzyk.png[i].cvnpart.muonid;
      prongNum = i;
    }
    else if (sr->vtx.elastic.fuzzyk.png[i].len >= 500)
    {
      maxVal = 1.0;
      prongNum = i;
    }
  }
  return prongNum;
}
*/
/*********************************************************************************************************************************/

/****************************************** Pion Prong Number ********************************************************************/
/*
unsigned int GetPionProngId(const caf::SRProxy *sr)
{
  unsigned int muonNum = GetMuonProngId(sr);
  if (muonNum == 10000 || sr->vtx.elastic.IsValid != true)
    return 10000;
  else
  {
    int prongNum = 0;
    float maxVal = -1;

    for (unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
    {
      if (sr->vtx.elastic.fuzzyk.png[i].cvnpart.pionid > maxVal &&
          i != muonNum)
      {
        maxVal = sr->vtx.elastic.fuzzyk.png[i].cvnpart.pionid;
        prongNum = i;
      }
    }
    return prongNum;
  }
}
*/
/*********************************************************************************************************************************/

/****************************************** Reconstructed |t| ********************************************************************/

const Var kRecot([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {
    double Emuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.energy);
    double Pxmuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.x);
    double Pymuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.y);
    double Pzmuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.z);

    double Epion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.energy);
    double Pxpion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.x);
    double Pypion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.y);
    double Pzpion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.z);

    double PtX = Pxmuon + Pxpion;
    double PtY = Pymuon + Pypion;

    double t = (((Emuon - Pzmuon) + (Epion - Pzpion)) * ((Emuon - Pzmuon) + (Epion - Pzpion)) + ((PtX * PtX) + (PtY * PtY)));
    //std::cout << "\n reconstructed values: " << " MuonE: " << Emuon << " muonPx: "  << Pxmuon << " muonPy: "  << Pymuon << " muonPz: " << Pzmuon << " PionE: " << Epion << " pionPx: " << Pxpion << " pionPy: " << Pypion << " pionPz: " << Pzpion << "\n";
    return t;
  }
});

/****************************************** Reconstructed |t| without using bpf ********************************************************************/

const Var kRecotNoBPF([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);
  double Mpion = 0.139570;
  double Mmuon = 0.105658;
  TVector3 muonDir(0, 0, 0);
  TVector3 pionDir(0, 0, 0);
  TVector3 beamdir = NuMIBeamDirection(caf::kNEARDET);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {
    double Emuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].calE);
    double Pmuon = sqrt(pow(Emuon, 2) - pow(Mmuon, 2));

    muonDir.SetXYZ(sr->vtx.elastic.fuzzyk.png[muonNum].dir.x,
                   sr->vtx.elastic.fuzzyk.png[muonNum].dir.y,
                   sr->vtx.elastic.fuzzyk.png[muonNum].dir.z);

    double CosThetaMuon = double(muonDir.Unit().Dot(beamdir));
    double SinThetaMuon = sqrt((1 - pow(CosThetaMuon, 2)));

    double Plmuon = Pmuon * CosThetaMuon;
    double Ptmuon = Pmuon * SinThetaMuon;

    double Epion = double(sr->vtx.elastic.fuzzyk.png[pionNum].calE);
    double Ppion = sqrt(pow(Epion, 2) - pow(Mpion, 2));

    pionDir.SetXYZ(sr->vtx.elastic.fuzzyk.png[pionNum].dir.x,
                   sr->vtx.elastic.fuzzyk.png[pionNum].dir.y,
                   sr->vtx.elastic.fuzzyk.png[pionNum].dir.z);

    double CosThetapion = double(pionDir.Unit().Dot(beamdir));
    double SinThetapion = sqrt((1 - pow(CosThetapion, 2)));

    double Plpion = Ppion * CosThetapion;
    double Ptpion = Ppion * SinThetapion;

    double t = (pow(((Emuon - Plmuon) + (Epion - Plpion)), 2) + pow((Ptmuon + Ptpion), 2));
    std::cout << "\n reconstructed values: "
              << " MuonE: " << Emuon
              << " Pmuon: " << Pmuon
              << " CosThetaMuon: " << CosThetaMuon
              << " SinThetaMuon: " << SinThetaMuon
              << " muonPl: " << Plmuon
              << " muonPt: " << Ptmuon
              << " PionE: " << Epion
              << " Ppion: " << Ppion
              << " CosThetaPion: " << CosThetapion
              << " SinThetaPion: " << SinThetapion
              << " pionPl: " << Plpion
              << " pionPt: " << Ptpion << "\n";
    return t;
  }
});

/*********************************************************************************************************************************/
/****************************************** Reconstructed |t| with E^2 = P^2 +m^2 ********************************************************************/
// const Var kRecot2([](const caf::SRProxy){
//   unsigned int muonNum = GetMuonProngId(sr);
//   unsigned int pionNum = GetPionProngId(sr);

//   if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
//     return -100.0;

//   else{
//     double Emuon =

//   }

// });

/****************************************** transverse momentum balance ********************************************************************/

const Var kRecoPt([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {
    //double Emuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.energy);
    double Pxmuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.x);
    double Pymuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.y);
    //double Pzmuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.z);

    //double Epion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.energy);
    double Pxpion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.x);
    double Pypion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.y);
    //double Pzpion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.z);

    double PtX = Pxmuon + Pxpion;
    double PtY = Pymuon + Pypion;

    if (isnan(Pxmuon) || isnan(Pymuon) || isnan(Pxpion) || isnan(Pypion))
      return -100.0;
    else
    {
      double t = ((PtX * PtX) + (PtY * PtY));
      //printf("\n Pxmuon: %f Pymuon: %f Pxpion: %f Pypion: %f", Pxmuon, Pymuon, Pxpion, Pypion);
      //std::cout << "|t| value: " << t << "\n";
      return t;
    }
  }
});

/*********************************************************************************************************************************/
/****************************************** Reconstructed Muon Energy ********************************************************************/

const Var kRecoEMuon([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  //unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true ||
      sr->vtx.elastic.fuzzyk.npng == 0 ||
      sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.IsValid == false)
    return -100.0;

  else
  {
    double Emuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.energy);
    return Emuon;
  }
});

/****************************************** Reconstructed Muon Kinetic Energy ********************************************************************/

const Var kRecoKEMuon([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  //unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true ||
      sr->vtx.elastic.fuzzyk.npng == 0 ||
      sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.IsValid == false)
    return -100.0;

  else
  {
    double Emuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.energy);
    return (Emuon - 0.105658);
  }
});

/****************************************** Reconstructed Pion Energy ********************************************************************/

const Var kRecoEPion([](const caf::SRProxy *sr) {
  //unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true ||
      sr->vtx.elastic.fuzzyk.npng == 0 ||
      sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.IsValid == false)
    return -100.0;

  else
  {
    double Epion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.muon.energy);
    return Epion;
  }
});

/****************************************** Reconstructed Pion Kinetic Energy ********************************************************************/

const Var kRecoKEPion([](const caf::SRProxy *sr) {
  //unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true ||
      sr->vtx.elastic.fuzzyk.npng == 0 ||
      sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.IsValid == false)
    return -100.0;

  else
  {
    double Epion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.muon.energy);
    return (Epion - 0.139570);
  }
});

/****************************************** Pion dE/dX Logliklihood **************************************************************/

// const Var pionLLDEDX([](const caf::SRProxy *sr) { //This is empty on prod 5
//   unsigned int prongNum = GetPionProngId(sr);
//   if (prongNum == 10000 || sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || sr->vtx.elastic.fuzzyk.png[prongNum].bpf.pion.IsValid != true)
//     return -10.0;
//   else
//     std::cout << double(sr->vtx.elastic.fuzzyk.png[prongNum].bpf.pion.dEdXLL) << "\n";
//   //return double(sr->vtx.elastic.fuzzyk.png[prongNum].bpf.pion.dEdXLL);

// });

/*********************************************************************************************************************************/
/****************************************** Muon dE/dX Logliklihood **************************************************************/

// const Var muonLLDEDX([](const caf::SRProxy *sr) { //This is empty on prod 5
//   unsigned int prongNum = GetMuonProngId(sr);
//   if (prongNum == 10000 || sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || sr->vtx.elastic.fuzzyk.png[prongNum].bpf.muon.IsValid != true)
//     return -10.0;
//   else
//     return double(sr->vtx.elastic.fuzzyk.png[prongNum].bpf.muon.dEdXLL);

// });

/*********************************************************************************************************************************/

/****************************************** Kalman dE/dX Logliklihood **************************************************************/
// const Var kalmandEdXllh([](const caf::SRProxy *sr) {
//   return double(sr->trk.kalman.tracks.dedxllh);
// });
/*********************************************************************************************************************************/
const Var muonCalE([](const caf::SRProxy *sr) {
  unsigned int prongNum = GetMuonProngId(sr);
  if (prongNum == 10000 || sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng != 2 || sr->vtx.elastic.fuzzyk.npng2d != 0)
    return -10.0;
  else
    return double(sr->vtx.elastic.fuzzyk.png[prongNum].calE);
});

const Var pionCalE([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  if (muonNum == 10000 || sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng != 2 || sr->vtx.elastic.fuzzyk.npng2d != 0)
    return -10.0;
  else
  {
    int prongNum = 0;
    float maxVal = -1;
    // if(sr->vtx.nelastic==0 || sr->vtx.elastic.fuzzyk.npng==0)
    //   return -5.0;

    for (unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
    {
      if (sr->vtx.elastic.fuzzyk.png[i].cvnpart.pionid > maxVal &&
          i != muonNum)
      {
        maxVal = sr->vtx.elastic.fuzzyk.png[i].cvnpart.pionid;
        prongNum = i;
      }
    }
    //if (sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
    return double(sr->vtx.elastic.fuzzyk.png[prongNum].calE);

    //else
    //return -10.0;
  }
  //if(sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
  //else if (! ())
});

const Var vertexAct = calE - muonCalE - pionCalE;

const Var track2D([](const caf::SRProxy *sr) {
  return double(sr->trk.kalman.ntracks2d);
});

const Var prongs2D([](const caf::SRProxy *sr) {
  if (sr->vtx.elastic.IsValid == true)
    return double(sr->vtx.elastic.fuzzyk.npng2d);
  else
    return -5.0;
});

const Var Topology([](const caf::SRProxy *sr) {
  size_t npng = sr->vtx.elastic.fuzzyk.npng;
  size_t npng2d = sr->vtx.elastic.fuzzyk.npng2d;

  if (npng == 0 && npng2d == 0)
    return 1.0;
  else if (npng == 0 && npng2d == 1)
    return 2.0;
  else if (npng == 0 && npng2d == 2)
    return 3.0;
  else if (npng == 0 && npng2d == 3)
    return 4.0;
  else if (npng == 0 && npng2d == 4)
    return 5.0;
  else if (npng == 0 && npng2d >= 5)
    return 6.0;
  else if (npng == 1 && npng2d == 0)
    return 7.0;
  else if (npng == 1 && npng2d == 1)
    return 8.0;
  else if (npng == 1 && npng2d == 2)
    return 9.0;
  else if (npng == 1 && npng2d == 3)
    return 10.0;
  else if (npng == 1 && npng2d == 4)
    return 11.0;
  else if (npng == 1 && npng2d >= 5)
    return 12.0;
  else if (npng == 2 && npng2d == 0)
    return 13.0;
  else if (npng == 2 && npng2d == 1)
    return 14.0;
  else if (npng == 2 && npng2d == 2)
    return 15.0;
  else if (npng == 2 && npng2d == 3)
    return 16.0;
  else if (npng == 2 && npng2d == 4)
    return 17.0;
  else if (npng == 2 && npng2d >= 5)
    return 18.0;
  else if (npng == 3 && npng2d == 0)
    return 19.0;
  else if (npng == 3 && npng2d == 1)
    return 20.0;
  else if (npng == 3 && npng2d == 2)
    return 21.0;
  else if (npng == 3 && npng2d == 3)
    return 22.0;
  else if (npng == 3 && npng2d == 4)
    return 23.0;
  else if (npng == 3 && npng2d >= 5)
    return 24.0;
  else if (npng == 4 && npng2d == 0)
    return 25.0;
  else if (npng == 4 && npng2d == 1)
    return 26.0;
  else if (npng == 4 && npng2d == 2)
    return 27.0;
  else if (npng == 4 && npng2d == 3)
    return 28.0;
  else if (npng == 4 && npng2d == 4)
    return 29.0;
  else if (npng == 4 && npng2d >= 5)
    return 30.0;
  else if (npng >= 5)
    return 31.0;

  return -10.0;
});

const Var mostEnergetic([](const caf::SRProxy *sr) {
  if (sr->mc.nu.size() == 0)
    return -5.0;
  double ParticleEnergy = -10.0;

  for (unsigned int i = 0; i < sr->mc.nu[0].prim.size(); i++)
  {

    if (sr->mc.nu[0].prim[i].pdg != 13 && sr->mc.nu[0].prim[i].p.E > ParticleEnergy)
      ParticleEnergy = double(sr->mc.nu[0].prim[i].p.E);
  }

  return ParticleEnergy;
});

const Var vtxX([](const caf::SRProxy *sr) {
  if (sr->vtx.elastic.IsValid == true)
    return double(sr->vtx.elastic.vtx.x);
  else
    return double(-500.0);
});

const Var vtxY([](const caf::SRProxy *sr) {
  if (sr->vtx.elastic.IsValid == true)
    return double(sr->vtx.elastic.vtx.y);
  else
    return double(-500.0);
});

const Var vtxZ([](const caf::SRProxy *sr) {
  if (sr->vtx.elastic.IsValid == true)
    return double(sr->vtx.elastic.vtx.z);
  else
    return double(-500.0);
});

const Var ProngPionAngleSq([](const caf::SRProxy *sr) {
  double CosThetaPi = 0;
  unsigned int muonNum = GetMuonProngId(sr);

  TVector3 dir(0, 0, 0);
  TVector3 beamdir = NuMIBeamDirection(caf::kNEARDET);
  bool hasPion = false;

  if (muonNum == 10000 || sr->vtx.elastic.IsValid == true || sr->vtx.elastic.fuzzyk.npng != 2 || sr->vtx.elastic.fuzzyk.npng2d != 0)
    return -10.0;
  else
  {
    int prongNum = 0;
    float maxVal = -1;
    // if(sr->vtx.nelastic==0 || sr->vtx.elastic.fuzzyk.npng==0)
    //   return -5.0;

    for (unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
    {
      if (sr->vtx.elastic.fuzzyk.png[i].cvnpart.pionid > maxVal &&
          i != muonNum)
      {
        maxVal = sr->vtx.elastic.fuzzyk.png[i].cvnpart.pionid;
        prongNum = i;
      }
    }
    if (sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
    {
      dir.SetXYZ(sr->vtx.elastic.fuzzyk.png[prongNum].dir.x, sr->vtx.elastic.fuzzyk.png[prongNum].dir.y, sr->vtx.elastic.fuzzyk.png[prongNum].dir.z);
      CosThetaPi = double(dir.Unit().Dot(beamdir));
      //return double(sr->vtx.elastic.fuzzyk.png[prongNum].calE);
      //std::cout<<(TMath::ACos(CosThetaPi)) * (TMath::ACos(CosThetaPi))<<"\n";
      return (TMath::ACos(CosThetaPi)) * (TMath::ACos(CosThetaPi));
    }

    else
      return -10.0;
    //if(sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
    //else if (! ())
  }
});

const Var EthetaSq = pionCalE * ProngPionAngleSq;

const Var CosPionAngle([](const caf::SRProxy *sr) {
  double CosThetaPi = 0;
  unsigned int muonNum = GetMuonProngId(sr);

  TVector3 dir(0, 0, 0);
  TVector3 beamdir = NuMIBeamDirection(caf::kNEARDET);
  bool hasPion = false;

  if (muonNum == 10000 || sr->vtx.elastic.IsValid != true)
    return -10.0;
  else
  {
    int prongNum = 0;
    float maxVal = -1;
    // if(sr->vtx.nelastic==0 || sr->vtx.elastic.fuzzyk.npng==0)
    //   return -5.0;

    for (unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
    {
      if (sr->vtx.elastic.fuzzyk.png[i].cvnpart.pionid > maxVal &&
          i != muonNum)
      {
        maxVal = sr->vtx.elastic.fuzzyk.png[i].cvnpart.pionid;
        prongNum = i;
      }
    }
    if (sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
    {
      dir.SetXYZ(sr->vtx.elastic.fuzzyk.png[prongNum].dir.x, sr->vtx.elastic.fuzzyk.png[prongNum].dir.y, sr->vtx.elastic.fuzzyk.png[prongNum].dir.z);
      CosThetaPi = double(dir.Unit().Dot(beamdir));
      //return double(sr->vtx.elastic.fuzzyk.png[prongNum].calE);
      //std::cout<<(TMath::ACos(CosThetaPi)) * (TMath::ACos(CosThetaPi))<<"\n";
      return double(1 - CosThetaPi);
    }

    else
      return -10.0;
    //if(sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
    //else if (! ())
  }
});

//const Var kRecoPiE([](const caf::SRProxy *sr){

//});

const Var RecoEeta = pionCalE * CosPionAngle;

const Var pionProngId([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  if (muonNum == 10000 || sr->vtx.elastic.IsValid != true)
    return -10.0;
  else
  {
    int prongNum = 0;
    float maxVal = -1;
    // if(sr->vtx.nelastic==0 || sr->vtx.elastic.fuzzyk.npng==0)
    //   return -5.0;

    for (unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
    {
      if (sr->vtx.elastic.fuzzyk.png[i].cvnpart.pionid > maxVal &&
          i != muonNum)
      {
        maxVal = sr->vtx.elastic.fuzzyk.png[i].cvnpart.pionid;
        //std::cout<<maxVal<<"\n";
        prongNum = i;
      }
    }
    //if (sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
    return double(maxVal);

    //else
    //return -10.0;
  }
  //if(sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
  //else if (! ())
});

const Var muonProngId([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  if (muonNum == 10000 || sr->vtx.elastic.IsValid != true)
    return -10.0;
  else
  {
    int prongNum = 0;
    float maxVal = -1;
    // if(sr->vtx.nelastic==0 || sr->vtx.elastic.fuzzyk.npng==0)
    //   return -5.0;

    for (unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
    {
      if (sr->vtx.elastic.fuzzyk.png[i].cvnpart.muonid > maxVal)
      {
        maxVal = sr->vtx.elastic.fuzzyk.png[i].cvnpart.muonid;
        //std::cout<<maxVal<<"\n";
        prongNum = i;
      }
    }
    //if (sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
    return double(maxVal);

    //else
    //return -10.0;
  }
  //if(sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
  //else if (! ())
});

const Var protonProngId([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  if (muonNum == 10000 || sr->vtx.elastic.IsValid != true)
    return -10.0;
  else
  {
    int prongNum = 0;
    float maxVal = -1;
    // if(sr->vtx.nelastic==0 || sr->vtx.elastic.fuzzyk.npng==0)
    //   return -5.0;

    for (unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
    {
      if (sr->vtx.elastic.fuzzyk.png[i].cvnpart.protonid > maxVal &&
          i != muonNum)
      {
        maxVal = sr->vtx.elastic.fuzzyk.png[i].cvnpart.protonid;
        prongNum = i;
      }
    }
    //if (sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
    return double(maxVal);

    //else
    //return -10.0;
  }
  //if(sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
  //else if (! ())
});

// const Var protonProngId([](const caf::SRProxy *sr) {
//   unsigned int muonNum = GetMuonProngId(sr);
//   if (muonNum == 10000 || sr->vtx.elastic.IsValid != true)
//     return -10.0;
//   else
//   {
//     int prongNum = 0;
//     float maxVal = -1;
//     // if(sr->vtx.nelastic==0 || sr->vtx.elastic.fuzzyk.npng==0)
//     //   return -5.0;

//     for (unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
//     {
//       if (sr->vtx.elastic.fuzzyk.png[i].cvnpart.protonid > maxVal &&
//           i != muonNum)
//       {
//         maxVal = sr->vtx.elastic.fuzzyk.png[i].cvnpart.protonid;
//         prongNum = i;
//       }
//     }
//     //if (sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
//     return double(maxVal);

//     //else
//     //return -10.0;
//   }
//   //if(sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
//   //else if (! ())
// });
