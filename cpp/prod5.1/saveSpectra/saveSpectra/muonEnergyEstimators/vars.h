#include "headers.h"

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
/****************************************** muon Prong Length ********************************************************************/
const Var kMuonProngLen([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	  sr->vtx.elastic.fuzzyk.npng == 0 || 
	  muonNum == pionNum)
    return -100.0;

  else
  {

    double muonLen = double(sr->vtx.elastic.fuzzyk.png[muonNum].len);

	return muonLen;
  }
});
/****************************************** muon Prong Length Based on Single Particle CVN ********************************************************************/
const Var kMuonProngLenNew([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngIdNew(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	  sr->vtx.elastic.fuzzyk.npng == 0 || 
	  muonNum == pionNum)
    return -100.0;

  else
  {

	if(sr->vtx.elastic.fuzzyk.png[muonNum].truth.pdg != 13 )
		return -100.0;
	else
	{

    double muonLen = double(sr->vtx.elastic.fuzzyk.png[muonNum].len);

	return muonLen;
	}
  }
});

/****************************************** pion Prong Length ********************************************************************/
const Var kPionProngLen([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	  sr->vtx.elastic.fuzzyk.npng == 0 || 
	  muonNum == pionNum)
    return -100.0;

  else
  {

    double pionLen = double(sr->vtx.elastic.fuzzyk.png[pionNum].len);

	return pionLen;
  }
});
/****************************************** muon Prong True KE ********************************************************************/
const Var kMuonProngTrueKE([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	  sr->vtx.elastic.fuzzyk.npng == 0 || 
	  muonNum == pionNum)
    return -100.0;

  else
  {

    double muonKE = double(sr->vtx.elastic.fuzzyk.png[muonNum].truth.p.E);

	return (muonKE - 0.105658);
  }
});
/****************************************** muon Prong True KE based on Sigle Particle CVN ******************************************/
const Var kMuonProngTrueKENew([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngIdNew(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	    sr->vtx.elastic.fuzzyk.npng == 0   || 
      sr->vtx.elastic.fuzzyk.png[muonNum].truth.pdg != 13 ||
	  muonNum == pionNum)
    return -100.0;
  else
  {
    double muonKE = double(sr->vtx.elastic.fuzzyk.png[muonNum].truth.p.E);
	  return (muonKE - 0.105658);
  }
});
/****************************************** Estimated muon Prong True KE based on Sigle Particle CVN ******************************/
const Var kMuonProngEstKE([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngIdNew(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	    sr->vtx.elastic.fuzzyk.npng == 0   || 
      sr->vtx.elastic.fuzzyk.png[muonNum].truth.pdg != 13 ||
	  muonNum == pionNum)
    return -100.0;
  else
  {
    double muonLen          = double(sr->vtx.elastic.fuzzyk.png[muonNum].len);
    double calculatedMuonKE = (0.00206123*muonLen + 0.029106);
	  return calculatedMuonKE;
  }
});

/****************************************** Calculated muon Pl using muon energy estimator ******************************/
const Var kMuonPlCalculated([](const caf::SRProxy *sr) {
	unsigned int muonNum = GetMuonProngIdNew(sr);
  	unsigned int pionNum = GetPionProngId(sr);

  	if (sr->vtx.elastic.IsValid != true    || 
	  	sr->vtx.elastic.fuzzyk.npng == 0   || 
      	sr->vtx.elastic.fuzzyk.png[muonNum].truth.pdg != 13 ||
	  	muonNum == pionNum)
    	return -100.0;

  	else
  	{	
		double muonLen = double(sr->vtx.elastic.fuzzyk.png[muonNum].len);
	    double muonKE  = (0.00206123*muonLen + 0.029106);

		double Pmuon   = sqrt( muonKE*(muonKE+ ( 2*0.105658 )) );
  		TVector3 muonDir(0, 0, 0);
  		TVector3 beamdir = NuMIBeamDirection(caf::kNEARDET);
  		muonDir.SetXYZ(sr->vtx.elastic.fuzzyk.png[muonNum].dir.x,
                   	   sr->vtx.elastic.fuzzyk.png[muonNum].dir.y,
                   	   sr->vtx.elastic.fuzzyk.png[muonNum].dir.z);

    	double CosThetaMuon = double(muonDir.Unit().Dot(beamdir));
    	//double SinThetaMuon = sqrt((1 - pow(CosThetaMuon, 2)));

    	double PlMuon = Pmuon * CosThetaMuon;
    	//double PtMuon = Pmuon * SinThetaMuon;

  		return PlMuon;
  	}
});
/****************************************** Calculated muon Pt using muon energy estimator ******************************/
const Var kMuonPtCalculated([](const caf::SRProxy *sr) {
	unsigned int muonNum = GetMuonProngIdNew(sr);
  	unsigned int pionNum = GetPionProngId(sr);

  	if (sr->vtx.elastic.IsValid != true    || 
	  	sr->vtx.elastic.fuzzyk.npng == 0   || 
      	sr->vtx.elastic.fuzzyk.png[muonNum].truth.pdg != 13 ||
	  	muonNum == pionNum)
    	return -100.0;

  	else
  	{	
		double muonLen = double(sr->vtx.elastic.fuzzyk.png[muonNum].len);
	    double muonKE  = (0.00206123*muonLen + 0.029106);

		double Pmuon   = sqrt( muonKE*(muonKE+ ( 2*0.105658 )) );
  		TVector3 muonDir(0, 0, 0);
  		TVector3 beamdir = NuMIBeamDirection(caf::kNEARDET);
  		muonDir.SetXYZ(sr->vtx.elastic.fuzzyk.png[muonNum].dir.x,
                   	   sr->vtx.elastic.fuzzyk.png[muonNum].dir.y,
                   	   sr->vtx.elastic.fuzzyk.png[muonNum].dir.z);

    	double CosThetaMuon = double(muonDir.Unit().Dot(beamdir));
    	double SinThetaMuon = sqrt((1 - pow(CosThetaMuon, 2)));

    	//double PlMuon = Pmuon * CosThetaMuon;
    	double PtMuon = Pmuon * SinThetaMuon;

  		return PtMuon;
  	}
});
/****************************************** True muon Pl  *****************************************************************/
const Var kTrueMuonPl([](const caf::SRProxy *sr) {
	unsigned int muonNum = GetMuonProngIdNew(sr);
  	unsigned int pionNum = GetPionProngId(sr);

  	if (sr->vtx.elastic.IsValid != true    || 
	  	sr->vtx.elastic.fuzzyk.npng == 0   || 
      	sr->vtx.elastic.fuzzyk.png[muonNum].truth.pdg != 13 ||
	  	muonNum == pionNum)
    	return -100.0;

  	else
  	{	
  		return double(sr->vtx.elastic.fuzzyk.png[muonNum].truth.p.pz);
  	}
});
/****************************************** True muon Pt *****************************************************************/
const Var kTrueMuonPt([](const caf::SRProxy *sr) {
	unsigned int muonNum = GetMuonProngIdNew(sr);
  	unsigned int pionNum = GetPionProngId(sr);

  	if (sr->vtx.elastic.IsValid != true    || 
	  	sr->vtx.elastic.fuzzyk.npng == 0   || 
      	sr->vtx.elastic.fuzzyk.png[muonNum].truth.pdg != 13 ||
	  	muonNum == pionNum)
    	return -100.0;

  	else
  	{	
		double muonPx = double(sr->vtx.elastic.fuzzyk.png[muonNum].truth.p.px);
	    double muonPy = double(sr->vtx.elastic.fuzzyk.png[muonNum].truth.p.py);

		double PtMuon     = sqrt( pow(muonPx, 2) + pow(muonPy, 2) );
  		return PtMuon;
  	}
});
/***************************** New Pl Resolution  ************************************************************************/
const Var kMuonPlResolution([](const caf::SRProxy *sr) {
	unsigned int muonNum = GetMuonProngIdNew(sr);
  	unsigned int pionNum = GetPionProngId(sr);

  	if (sr->vtx.elastic.IsValid != true    || 
	  	sr->vtx.elastic.fuzzyk.npng == 0   || 
      	sr->vtx.elastic.fuzzyk.png[muonNum].truth.pdg != 13 ||
	  	muonNum == pionNum)
    	return -100.0;

  	else
  	{	
		double muonLen = double(sr->vtx.elastic.fuzzyk.png[muonNum].len);
	    double muonKE  = (0.00206123*muonLen + 0.029106);

		double Pmuon   = sqrt( muonKE*(muonKE+ ( 2*0.105658 )) );
  		TVector3 muonDir(0, 0, 0);
  		TVector3 beamdir = NuMIBeamDirection(caf::kNEARDET);
  		muonDir.SetXYZ(sr->vtx.elastic.fuzzyk.png[muonNum].dir.x,
                   	   sr->vtx.elastic.fuzzyk.png[muonNum].dir.y,
                   	   sr->vtx.elastic.fuzzyk.png[muonNum].dir.z);

    	double CosThetaMuon = double(muonDir.Unit().Dot(beamdir));
    	//double SinThetaMuon = sqrt((1 - pow(CosThetaMuon, 2)));

    	double PlMuon = Pmuon * CosThetaMuon;
		double truePlMuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].truth.p.pz);
    	//double PtMuon = Pmuon * SinThetaMuon;
		if (truePlMuon == 0)
			return 0.0;
		else
			{
			double res = (PlMuon - truePlMuon)/truePlMuon;
  			return res;
			}  	
	}
});
/************************** New Pt Resolution ***********************************************************************/
const Var kMuonPtResolution([](const caf::SRProxy *sr) {
	unsigned int muonNum = GetMuonProngIdNew(sr);
  	unsigned int pionNum = GetPionProngId(sr);

  	if (sr->vtx.elastic.IsValid != true    || 
	  	sr->vtx.elastic.fuzzyk.npng == 0   || 
      	sr->vtx.elastic.fuzzyk.png[muonNum].truth.pdg != 13 ||
	  	muonNum == pionNum)
    	return -100.0;

  	else
  	{	
		double muonLen = double(sr->vtx.elastic.fuzzyk.png[muonNum].len);
	    double muonKE  = (0.00206123*muonLen + 0.029106);

		double Pmuon   = sqrt( muonKE*(muonKE+ ( 2*0.105658 )) );
  		TVector3 muonDir(0, 0, 0);
  		TVector3 beamdir = NuMIBeamDirection(caf::kNEARDET);
  		muonDir.SetXYZ(sr->vtx.elastic.fuzzyk.png[muonNum].dir.x,
                   	   sr->vtx.elastic.fuzzyk.png[muonNum].dir.y,
                   	   sr->vtx.elastic.fuzzyk.png[muonNum].dir.z);

    	double CosThetaMuon = double(muonDir.Unit().Dot(beamdir));
    	double SinThetaMuon = sqrt((1 - pow(CosThetaMuon, 2)));

    	//double PlMuon = Pmuon * CosThetaMuon;
    	double PtMuon = Pmuon * SinThetaMuon;
		double muonPx = double(sr->vtx.elastic.fuzzyk.png[muonNum].truth.p.px);
	    double muonPy = double(sr->vtx.elastic.fuzzyk.png[muonNum].truth.p.py);

		double truePtMuon     = sqrt( pow(muonPx, 2) + pow(muonPy, 2) );
		
		if (truePtMuon == 0)
			return 0.0;
		else
		{
			double res = ( PtMuon - truePtMuon )/truePtMuon;
			return res; 
		}
  	}
});

/****************************************** pion Prong Length ********************************************************************/
const Var kPionProngTrueKE([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	  sr->vtx.elastic.fuzzyk.npng == 0 || 
	  muonNum == pionNum)
    return -100.0;

  else
  {

    double pionKE = double(sr->vtx.elastic.fuzzyk.png[pionNum].truth.p.E);

	return ( pionKE - 0.13957 );
  }
});
/****************************************** Interaction Type ********************************************************************/
const Var InteractionType([](const caf::SRProxy *sr) -> int {
  if (sr->mc.nnu == 0)
    return -100.0;
  else
    return (sr->mc.nu[0].mode);
});



/****************************************** muon Pt Resolution ********************************************************************/

const Var kMuonPtRes([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {

    double RecoMuonPx = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.x);
    double TrueMuonPx = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.truth.p.px);
	double RecoMuonPy = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.y);
    double TrueMuonPy = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.truth.p.py);

	double TrueMuonPtSq = pow(TrueMuonPx,2) + pow(TrueMuonPy,2);
	double TrueMuonPt = sqrt(TrueMuonPtSq);

	double RecoMuonPtSq = pow(RecoMuonPx,2) + pow(RecoMuonPy,2);
	double RecoMuonPt = sqrt(RecoMuonPtSq);

	double res = (RecoMuonPt - TrueMuonPt)/TrueMuonPt;
	
	return res;
  }
});

/****************************************** pion Pt Resolution ********************************************************************/

const Var kPionPtRes([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {

    double RecoPionPx = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.x);
    double TruePionPx = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.truth.p.px);
	double RecoPionPy = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.y);
    double TruePionPy = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.truth.p.py);

	double TruePionPtSq = pow(TruePionPx,2) + pow(TruePionPy,2);
	double TruePionPt = sqrt(TruePionPtSq);

	double RecoPionPtSq = pow(RecoPionPx,2) + pow(RecoPionPy,2);
	double RecoPionPt = sqrt(RecoPionPtSq);

	double res = (RecoPionPt - TruePionPt)/TruePionPt;
	
	return res;
  }
});
/****************************************** Muon Angle Resolution ********************************************************************/

const Var kMuonAngleRes([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	  sr->vtx.elastic.fuzzyk.npng == 0 || 
	  muonNum == pionNum               || 
      sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {

    double TrueMuonPx = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.truth.p.px);
    double TrueMuonPy = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.truth.p.py);
    double TrueMuonPz = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.truth.p.pz);

    double RecoMuonPx = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.x);
	double RecoMuonPy = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.y);
	double RecoMuonPz = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.z);

	TVector3 tMuMomentum(TrueMuonPx,
						 TrueMuonPy, 
						 TrueMuonPz); 

  	TVector3 beamdir = NuMIBeamDirection(caf::kNEARDET);

	double trueCosTheta = double( tMuMomentum.Unit().Dot(beamdir) );

	TVector3 rMuMomentum(RecoMuonPx,
						 RecoMuonPy, 
						 RecoMuonPz); 

	double recoCosTheta = double( rMuMomentum.Unit().Dot(beamdir) );

	double res = (recoCosTheta - trueCosTheta)/trueCosTheta;
	
	return res;
  }
});

/****************************************** pion Angle Resolution ********************************************************************/

const Var kPionAngleRes([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	  sr->vtx.elastic.fuzzyk.npng == 0 || 
	  muonNum == pionNum               || 
      sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {

    double TruePionPx = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.truth.p.px);
    double TruePionPy = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.truth.p.py);
    double TruePionPz = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.truth.p.pz);

    double RecoPionPx = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.x);
	double RecoPionPy = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.y);
	double RecoPionPz = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.z);

	TVector3 tPiMomentum(TruePionPx,
						 TruePionPy, 
						 TruePionPz); 

  	TVector3 beamdir = NuMIBeamDirection(caf::kNEARDET);

	double trueCosTheta = double( tPiMomentum.Unit().Dot(beamdir) );

	TVector3 rPiMomentum(RecoPionPx,
						 RecoPionPy, 
						 RecoPionPz); 

	double recoCosTheta = double( rPiMomentum.Unit().Dot(beamdir) );

	double res = (recoCosTheta - trueCosTheta)/trueCosTheta;
	
	return res;
  }
});



/****************************************** muon Px Resolution ********************************************************************/

const Var kMuonPxRes([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {

    double RecoMuonPx = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.x);
    double TrueMuonPx = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.truth.p.px);
	double res = (RecoMuonPx - TrueMuonPx)/TrueMuonPx;
	
	return res;
  }
});

/****************************************** muon Py Resolution ********************************************************************/

const Var kMuonPyRes([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {

    double RecoMuonPy = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.y);
    double TrueMuonPy = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.truth.p.py);
	double res = (RecoMuonPy - TrueMuonPy)/TrueMuonPy;
	
	return res;
  }
});

/****************************************** muon Pz Resolution ********************************************************************/

const Var kMuonPzRes([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {

    double RecoMuonPz = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.momentum.z);
    double TrueMuonPz = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.truth.p.pz);
	double res = (RecoMuonPz - TrueMuonPz)/TrueMuonPz;
	
	return res;
  }
});

/****************************************** pion Px Resolution ********************************************************************/

const Var kPionPxRes([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {

    double RecoPionPx = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.x);
    double TruePionPx = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.truth.p.px);
	double res = (RecoPionPx - TruePionPx)/TruePionPx;
	
	return res;
  }
});

/****************************************** pion Py Resolution ********************************************************************/

const Var kPionPyRes([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {

    double RecoPionPy = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.y);
    double TruePionPy = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.truth.p.py);
	double res = (RecoPionPy - TruePionPy)/TruePionPy;
	
	return res;
  }
});

/****************************************** pion Pz Resolution ********************************************************************/

const Var kPionPzRes([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {

    double RecoPionPz = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.momentum.z);
    double TruePionPz = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.truth.p.pz);
	double res = (RecoPionPz - TruePionPz)/TruePionPz;
	
	return res;
  }
});





/****************************************** Muon E Resolution using BPF ********************************************************************/
const Var kMuonERes([](const caf::SRProxy *sr){

  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {

    double RecoMuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.energy);
    double TrueMuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.truth.p.E);
	double res = (RecoMuon - TrueMuon)/TrueMuon;
	
	return res;
  }
 
});
/****************************************** Muon E Resolution using Muon energy Estimator ********************************************************************/
const Var kMuonEResNew([](const caf::SRProxy *sr){

  unsigned int muonNum = GetMuonProngIdNew(sr);
  unsigned int pionNum = GetPionProngIdNew(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	    sr->vtx.elastic.fuzzyk.npng == 0 ||
      sr->vtx.elastic.fuzzyk.png[muonNum].truth.pdg != 13 ||
      sr->vtx.elastic.fuzzyk.png[muonNum].len <40 || 
      sr->vtx.elastic.fuzzyk.png[muonNum].len >1170 || 
      muonNum == pionNum)
    return -100.0;
  else
  {
    double muonLen          = double(sr->vtx.elastic.fuzzyk.png[muonNum].len);
    double trueMuonKE       = double((sr->vtx.elastic.fuzzyk.png[muonNum].truth.p.E - 0.105658)); 
    double calculatedMuonKE = (0.00206123*muonLen + 0.029106);

      if ( trueMuonKE == 0 )
        return 0.0;
      else
        {
          double res = (calculatedMuonKE - trueMuonKE)/trueMuonKE;
          return res;
        }
  }
});

/****************************************** Muon E Resolution using Muon energy Estimator ********************************************************************/
const Var kMuonEResNewNoLenLimit([](const caf::SRProxy *sr){

  unsigned int muonNum = GetMuonProngIdNew(sr);
  unsigned int pionNum = GetPionProngIdNew(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	    sr->vtx.elastic.fuzzyk.npng == 0 ||
      sr->vtx.elastic.fuzzyk.png[muonNum].truth.pdg != 13 ||
      muonNum == pionNum)
    return -100.0;
  else
  {
    double muonLen          = double(sr->vtx.elastic.fuzzyk.png[muonNum].len);
    double trueMuonKE       = double((sr->vtx.elastic.fuzzyk.png[muonNum].truth.p.E - 0.105658)); 
    double calculatedMuonKE = (0.00206123*muonLen + 0.029106);

      if ( trueMuonKE == 0 )
        return 0.0;
      else
        {
          double res = (calculatedMuonKE - trueMuonKE)/trueMuonKE;
          return res;
        }
  }
});
/****************************************** Muon E Resolution using BPF ********************************************************************/
const Var kMuonEResBPFNew([](const caf::SRProxy *sr){

  unsigned int muonNum = GetMuonProngIdNew(sr);
  unsigned int pionNum = GetPionProngIdNew(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	    sr->vtx.elastic.fuzzyk.npng == 0 ||
      sr->vtx.elastic.fuzzyk.png[muonNum].truth.pdg != 13 ||
      muonNum == pionNum)
    return -100.0;
  else
  {
    double muonKEBPF        = double((sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.energy - 0.105658));
    double trueMuonKE       = double((sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.truth.p.E - 0.105658)); 
    //double calculatedMuonKE = (0.00206123*muonLen + 0.029106);

      if ( trueMuonKE == 0 )
        return 0.0;
      else
        {
          double res = (muonKEBPF - trueMuonKE)/trueMuonKE;
          return res;
        }
  }
});

/****************************************** Pion E Resolution using BPF ********************************************************************/
const Var kPionERes([](const caf::SRProxy *sr){

  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return -100.0;

  else
  {

    double RecoPion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.energy);
    double TruePion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.truth.p.E);
	double res = (RecoPion - TruePion)/TruePion;
	
	return res;
  }
 
});






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
              << " MuonE: "  << Emuon
              << " Pmuon: "  << Pmuon 
              << " CosThetaMuon: " << CosThetaMuon
              << " SinThetaMuon: " << SinThetaMuon
              << " muonPl: " << Plmuon
              << " muonPt: " << Ptmuon
              << " PionE: "  << Epion
              << " Ppion: "  << Ppion
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

/****************************************** BPF True Muon Energy ********************************************************************/

const Var kTrueBPFMuonE([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  //unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true ||
      sr->vtx.elastic.fuzzyk.npng == 0 ||
      sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.IsValid == false)
    return -100.0;

  else
  {
    double Emuon = double(sr->vtx.elastic.fuzzyk.png[muonNum].bpf.muon.truth.p.E);
    return Emuon;
  }
});

/****************************************** BPF True Pion Energy ********************************************************************/

const Var kTrueBPFPionE([](const caf::SRProxy *sr) {
  //unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true ||
      sr->vtx.elastic.fuzzyk.npng == 0 ||
      sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.IsValid == false)
    return -100.0;

  else
  {
    double Epion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.truth.p.E);
    return Epion;
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
    double Epion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.energy);
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
    double Epion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.energy);
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
/************************************************ New Calorimetric Energy for Prongs ******************************************************************/
const Var muonCalENew([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngIdNew(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	  sr->vtx.elastic.fuzzyk.npng == 0 || 
	  muonNum == pionNum)
    return -100.0;

  else
  {

	if(sr->vtx.elastic.fuzzyk.png[muonNum].truth.pdg != 13 )
		return -100.0;
	else
	{
		 double muonCalE = double(sr->vtx.elastic.fuzzyk.png[muonNum].calE);
		return muonCalE;
	}
  }
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
  return double(sr->vtx.elastic.fuzzyk.npng2d);
});
/********************************** Newly created Vertex Energy Variables **************************************/
const Var sliceVertexEVol10([](const caf::SRProxy *sr){
  return double(sr->vtx.elastic.slicevertexenergyvolume10);
});
const Var sliceVertexEVol20([](const caf::SRProxy *sr){
  return double(sr->vtx.elastic.slicevertexenergyvolume20);
});
const Var sliceVertexEVol30([](const caf::SRProxy *sr){
  return double(sr->vtx.elastic.slicevertexenergyvolume30);
});
const Var sliceVertexEVol40([](const caf::SRProxy *sr){
  return double(sr->vtx.elastic.slicevertexenergyvolume40);
});

const Var prong3DVertexEVol10([](const caf::SRProxy *sr){
  return double(sr->vtx.elastic.prong3dvertexenergyvolume10);
});
const Var prong3DVertexEVol20([](const caf::SRProxy *sr){
  return double(sr->vtx.elastic.prong3dvertexenergyvolume20);
});
const Var prong3DVertexEVol30([](const caf::SRProxy *sr){
  return double(sr->vtx.elastic.prong3dvertexenergyvolume30);
});
const Var prong3DVertexEVol40([](const caf::SRProxy *sr){
  return double(sr->vtx.elastic.prong3dvertexenergyvolume40);
});

const Var prong2DVertexEVol10([](const caf::SRProxy *sr){
  return double(sr->vtx.elastic.prong2dvertexenergyvolume10);
});
const Var prong2DVertexEVol20([](const caf::SRProxy *sr){
  return double(sr->vtx.elastic.prong2dvertexenergyvolume20);
});
const Var prong2DVertexEVol30([](const caf::SRProxy *sr){
  return double(sr->vtx.elastic.prong2dvertexenergyvolume30);
});
const Var prong2DVertexEVol40([](const caf::SRProxy *sr){
  return double(sr->vtx.elastic.prong2dvertexenergyvolume40);
});

const Var vertexE10([](const caf::SRProxy *sr){
  double slice   = sr->vtx.elastic.slicevertexenergyvolume10;
  double prong3D = sr->vtx.elastic.prong3dvertexenergyvolume10;
  double prong2D = sr->vtx.elastic.prong2dvertexenergyvolume10;

  //if((slice - prong3D - prong2D) >= 0)
    return (slice - prong3D - prong2D);
});

const Var SelectedVertexE10([](const caf::SRProxy *sr){
  double slice   = sr->vtx.elastic.slicevertexenergyvolume10;
  double prong3D = sr->vtx.elastic.prong3dvertexenergyvolume10;
  double prong2D = sr->vtx.elastic.prong2dvertexenergyvolume10;

  //if((slice - prong3D - prong2D) >= 0)
    if ( slice   != 0 && prong3D  != 0 )
      return (slice - prong3D - prong2D);
});

const Var vertexE20([](const caf::SRProxy *sr){
  double slice   = sr->vtx.elastic.slicevertexenergyvolume20;
  double prong3D = sr->vtx.elastic.prong3dvertexenergyvolume20;
  double prong2D = sr->vtx.elastic.prong2dvertexenergyvolume20;

  //if((slice - prong3D - prong2D) >= 0)
    return (slice - prong3D - prong2D);
});

const Var vertexE30([](const caf::SRProxy *sr){
  double slice   = sr->vtx.elastic.slicevertexenergyvolume30;
  double prong3D = sr->vtx.elastic.prong3dvertexenergyvolume30;
  double prong2D = sr->vtx.elastic.prong2dvertexenergyvolume30;

  //if((slice - prong3D - prong2D) >= 0)
    return (slice - prong3D - prong2D);
});

const Var vertexE40([](const caf::SRProxy *sr){
  double slice   = sr->vtx.elastic.slicevertexenergyvolume40;
  double prong3D = sr->vtx.elastic.prong3dvertexenergyvolume40;
  double prong2D = sr->vtx.elastic.prong2dvertexenergyvolume40;

  //if((slice - prong3D - prong2D) >= 0)
    return (slice - prong3D - prong2D);
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
/*
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
      if (sr->vtx.elastic.fuzzyk.png[i].spprongcvnpart5label.pionid > maxVal &&
          i != muonNum)
      {
        //maxVal = sr->vtx.elastic.fuzzyk.png[i].cvnpart.pionid;
        maxVal = sr->vtx.elastic.fuzzyk.png[i].spprongcvnpart5label.pionid;
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
*/
/********************************** Muon Prong ID  ************************************************************/
const Var muonProngId([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  if (muonNum == 10000 || sr->vtx.elastic.IsValid != true)
    return -10.0;
  else
    return double(sr->vtx.elastic.fuzzyk.png[muonNum].cvnpart.muonid);
});
/********************************** Muon Prong ID using Single Particle CVN ************************************/
const Var muonProngIdNew([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngIdNew(sr);
  if (muonNum == 10000 || sr->vtx.elastic.IsValid != true)
    return -10.0;
  else
    return double(sr->vtx.elastic.fuzzyk.png[muonNum].spprongcvnpart5label.muonid);
});
/********************************** Pion Prong ID  ************************************************************/
const Var pionProngId([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);
  if (muonNum == 10000 || sr->vtx.elastic.IsValid != true || muonNum == pionNum)
    return -10.0;
  else
    return double(sr->vtx.elastic.fuzzyk.png[pionNum].cvnpart.pionid);
});
/********************************** Pion Prong ID using Single Particle CVN ************************************/
const Var pionProngIdNew([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngIdNew(sr);
  unsigned int pionNum = GetPionProngIdNew(sr);
  if (muonNum == 10000 || sr->vtx.elastic.IsValid != true || muonNum == pionNum)
    return -10.0;
  else
    return double(sr->vtx.elastic.fuzzyk.png[pionNum].spprongcvnpart5label.pionid);
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
/********************************** muonID defined in NuMu CC inclusive analysis  ************************************/
//const Var muonIDNuMuCCInc = ana::muonid_classifier::kMuonID;


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
