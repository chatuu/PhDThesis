#include "headers.h"

/********************************** Get Muon Prong ID *********************************************************/

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

/********************************** Get Muon Prong ID based on Single Particle CVN *****************************/

unsigned int GetMuonProngIdNew(const caf::SRProxy *sr)
{
  int prongNum = 0;
  float maxVal = -1.0;

  if (sr->vtx.elastic.IsValid != true)
    return 10000;

  for (size_t i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
  {
    if (sr->vtx.elastic.fuzzyk.png[i].spprongcvnpart5label.muonid > maxVal)
    {
      maxVal = sr->vtx.elastic.fuzzyk.png[i].spprongcvnpart5label.muonid;
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

/********************************** Get Pion Prong ID *********************************************************/

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
/********************************** Get Pion Prong ID using Single Particle CVN *******************************/

unsigned int GetPionProngIdNew(const caf::SRProxy *sr)
{
  unsigned int muonNum = GetMuonProngIdNew(sr);
  if (muonNum == 10000 || sr->vtx.elastic.IsValid != true)
    return 10000;
  else
  {
    int prongNum = 0;
    float maxVal = -1;

    for (unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
    {
      if (sr->vtx.elastic.fuzzyk.png[i].spprongcvnpart5label.pionid > maxVal &&
          i != muonNum)
      {
        maxVal = sr->vtx.elastic.fuzzyk.png[i].spprongcvnpart5label.pionid;
        prongNum = i;
      }
    }
    return prongNum;
  }
}




/********************************** Prong based Muon ID cut using CVN ******************************************/
const Cut muonIDCutCVNBased([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true  || 
	  sr->vtx.elastic.fuzzyk.npng == 0 || 
	  muonNum == pionNum )
    return false;

  else
  {
    double muonID = double(sr->vtx.elastic.fuzzyk.png[muonNum].cvnpart.muonid);
	//std::cout << muonID << std::endl;
    if ( muonID > 0.775 )
      return true;
    else
      return false;
  }
});



/********************************** CC Coherent Event Selection ************************************************/
const Cut kMyNuMuCCCoh([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  return  (sr->mc.nu[0].iscc &&
          (sr->mc.nu[0].pdg == 14) &&
          (sr->mc.nu[0].pdgorig == 14) &&
          (sr->mc.nu[0].mode == caf::kCoh));
});
/********************************** NuMuCCInteractions Selections ************************************************/
const Cut NuMuCCInts([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  return  (sr->mc.nu[0].iscc &&
          (sr->mc.nu[0].pdg == 14) &&
          (sr->mc.nu[0].pdgorig == 14));
});


const Cut NoCuts([](const caf::SRProxy *sr) {
  return true;
});

const NuTruthCut kTrueMyNuMuCCCoh([](const caf::SRNeutrinoProxy *sr) {
  //if(sr->nnu==0) return false;
  return (sr->iscc &&
          (sr->pdg == 14) &&
          (sr->pdgorig == 14) &&
          (sr->mode == caf::kCoh));
});

const NuTruthCut kMyTrueFiducialST([](const caf::SRNeutrinoProxy *sr) {
  const TVector3 vtxmin(-130, -130, 100);
  const TVector3 vtxmax(140, 140, 1000);

  return (sr->vtx.X() < vtxmax.X() &&
          sr->vtx.X() > vtxmin.X() &&
          sr->vtx.Y() > vtxmin.Y() &&
          sr->vtx.Y() < vtxmax.Y() &&
          sr->vtx.Z() > vtxmin.Z() &&
          sr->vtx.Z() < vtxmax.Z());
});

const Cut kNuMuCCCohSig = kMyNuMuCCCoh && CutFromNuTruthCut(kMyTrueFiducialST);
const NuTruthCut kTrueNuMuCCCohSig = kTrueMyNuMuCCCoh && kMyTrueFiducialST;

const Cut kNCBkgd([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  return (!(sr->mc.nu[0].iscc)); // Is a NC interaction
});

const Cut kTwoProng([](const caf::SRProxy *sr) {
  if (sr->vtx.elastic.fuzzyk.npng == 2 && sr->vtx.elastic.fuzzyk.npng2d == 0)
    return true;
  else
    return false;
});

// Pion Energy cut for Signal Definition
const Cut kPionEnergyCut([](const caf::SRProxy *sr){
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);
  if (sr->vtx.elastic.IsValid != true  || 
      sr->vtx.elastic.fuzzyk.npng == 0 ||
      muonNum == pionNum               || 
      sr->vtx.elastic.fuzzyk.npng < 2)
      return false;
  else
  {
    double Epion = double(sr->vtx.elastic.fuzzyk.png[pionNum].bpf.pion.energy);
    if (Epion >= 0.3)
      return true;
    else
      return false;
  }
});

// To Get 1-3D Prong and 1-2D Prong
const Cut kOne3DProngOne2DProng([](const caf::SRProxy *sr) {
  if (sr->vtx.elastic.fuzzyk.npng == 1 && sr->vtx.elastic.fuzzyk.npng2d == 1)
    return true;
  else
    return false;
});
const Cut kTwo3DProngOne2DProng([](const caf::SRProxy *sr) {
  if (sr->vtx.elastic.fuzzyk.npng == 2 && sr->vtx.elastic.fuzzyk.npng2d == 1)
    return true;
  else
    return false;
});

//********************************** Replicating NuMuCC inclusive cuts ****************************************

//*************************************************************************************

const Cut kNumuMyQuality([](const caf::SRProxy *sr) {
  return (sr->trk.kalman.ntracks > 0 &&
          sr->slc.nhit > 20 &&
          sr->slc.ncontplanes > 4);
});

const Cut kIsFiducial([](const caf::SRProxy *sr) {
  const TVector3 vtxmin(-130, -130, 100);
  const TVector3 vtxmax(140, 140, 1000);

  return (sr->vtx.elastic.vtx.x < vtxmax.X() &&
          sr->vtx.elastic.vtx.x > vtxmin.X() &&
          sr->vtx.elastic.vtx.y > vtxmin.Y() &&
          sr->vtx.elastic.vtx.y < vtxmax.Y() &&
          sr->vtx.elastic.vtx.z > vtxmin.Z() &&
          sr->vtx.elastic.vtx.z < vtxmax.Z());
});

const Cut kNumuTightContainND([](const caf::SRProxy *sr) {
  if (sr->vtx.elastic.IsValid == false || sr->trk.kalman.ntracks < 1)
    return false;
  //int ibesttrk = kBestTrack(sr);

  // only primary muon track present in muon catcher
  else
  {
    for (int i = 0; i < int(sr->trk.kalman.ntracks); ++i)
    {
      //if( i == ibesttrk ) continue;
      if (sr->trk.kalman.tracks[i].start.Z() > 1275 ||
          sr->trk.kalman.tracks[i].stop.Z() > 1275)
        return false;
    }
    //printf("\nnshwlid: %i", int(sr->vtx.elastic.fuzzyk.nshwlid));
    //if (sr->vtx.elastic.fuzzyk.nshwlid == 0)
    //std::cout << "ZERO" << std::endl;
    // reconstructed showers all contained
    for (unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.nshwlid; ++i)
    {
      TVector3 start = sr->vtx.elastic.fuzzyk.png[i].shwlid.start;
      TVector3 stop = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop;
      if ((std::min(start.X(), stop.X()) < -180.0) ||
          (std::max(start.X(), stop.X()) > 180.0) ||
          (std::min(start.Y(), stop.Y()) < -180.0) ||
          (std::max(start.Y(), stop.Y()) > 180.0) ||
          (std::min(start.Z(), stop.Z()) < 20.0) ||
          (std::max(start.Z(), stop.Z()) > 1525.0))
        return false;
      else
        return true;
    }
  }
});

const Cut kMuonIDCut([](const caf::SRProxy *sr) {
  int prongNumVal = 0;
  int prongNumLen = 0;
  float maxVal = -1.0;
  float maxLen = -1.0;

  if (sr->vtx.elastic.IsValid == false)
    return false;

  for (size_t i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
  {
    if (sr->vtx.elastic.fuzzyk.png[i].cvnpart.muonid > maxVal)
    {
      maxVal = sr->vtx.elastic.fuzzyk.png[i].cvnpart.muonid;
      prongNumVal = i;
    }
  }
  for (size_t i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
  {
    if (sr->vtx.elastic.fuzzyk.png[i].len > maxLen)
    {
      maxLen = sr->vtx.elastic.fuzzyk.png[i].len;
      prongNumLen = i;
    }
  }

  if (prongNumVal == prongNumLen)
    return true;
  else
    return false;
});

//********************************** NuMu CC MuonID Funtion *********************************************************
const Var kBestMuonIDIndex([](const caf::SRProxy *sr) -> int
{
	int bestidx = sr->trk.kalman.idxmuonid;
	if (bestidx < 0 || bestidx >= (int) sr->trk.kalman.ntracks)
		return -99;
	return bestidx;
});

//********************************** NuMu CC MuonID Cut *************************************************************
const Cut NuMuCCMuonIDCut([](const caf::SRProxy *sr) {
	int bestidx = kBestMuonIDIndex(sr);
	if (bestidx < 0 || bestidx >= (int) sr->trk.kalman.ntracks)
		return false;
	double muonid = sr->trk.kalman.tracks[bestidx].muonid;
	if (muonid > 0.24)
		return true;
	else
		return false;
});


//********************************** |t| Cut ************************************************************************

const Cut ktCut([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return false;

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
    //std::cout << "|t| value: " << t << "\n";
    //return t;
    if (t < 0.1 || t == 0.1)
      return true;
    else
      return false;
  }
});

//********************************** |t| Cut Faliure Check ************************************************************************

const Cut ktCutFaliure([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  unsigned int pionNum = GetPionProngId(sr);

  if (sr->vtx.elastic.IsValid != true || sr->vtx.elastic.fuzzyk.npng == 0 || muonNum == pionNum || sr->vtx.elastic.fuzzyk.npng < 2)
    return false;

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

    //double PtX = Pxmuon + Pxpion;
    //double PtY = Pymuon + Pypion;

    //double t = (((Emuon - Pzmuon) + (Epion - Pzpion)) * ((Emuon - Pzmuon) + (Epion - Pzpion)) + ((PtX * PtX) + (PtY * PtY)));
    //std::cout << "|t| value: " << t << "\n";
    //return t;
    if (isnan(Emuon)  ||
        isnan(Pxmuon) ||
        isnan(Pymuon) ||
        isnan(Pzmuon) ||
        isnan(Epion)  ||
        isnan(Pxpion) ||
        isnan(Pypion) ||
        isnan(Pzpion))
      return true;
    else
      return false;
  }
});

//********************************** Old Pion ID Cut *******************************************************************

const Cut kPionIdOldCut([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  if (muonNum == 10000 || sr->vtx.elastic.IsValid != true)
    return false;
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

    //return double(maxVal);

    if (maxVal > 0.9 || maxVal == 0.9)
      return true;
    else
      return false;
    //else
    //return -10.0;
  }
  //if(sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
  //else if (! ())
});
//********************************** New Pion ID Cut *******************************************************************

const Cut kPionIdCut([](const caf::SRProxy *sr) {
  unsigned int muonNum = GetMuonProngId(sr);
  if (muonNum == 10000 || sr->vtx.elastic.IsValid != true)
    return false;
  else
  {
    int prongNum = 0;
    float maxVal = -1;
    
    for (unsigned int i = 0; i < sr->vtx.elastic.fuzzyk.npng; ++i)
    {
      if (sr->vtx.elastic.fuzzyk.png[i].spprongcvnpart5label.pionid > maxVal &&
          i != muonNum)
      {
        maxVal = sr->vtx.elastic.fuzzyk.png[i].spprongcvnpart5label.pionid;
        prongNum = i;
      }
    }
    if (maxVal > 0.7)
      return true;
    else
      return false;
  }
});

//********************************** Replicating NuMuCC inclusive cuts ****************************************

//********************************** Vertex Energy Cut *******************************************************************

const Cut kVertexECut([](const caf::SRProxy *sr) {
  if (sr->vtx.elastic.IsValid != true)
    return false;
  else
  {
    double slice     = sr->vtx.elastic.slicevertexenergyvolume10;
    double prong3D   = sr->vtx.elastic.prong3dvertexenergyvolume10;
    double prong2D   = sr->vtx.elastic.prong2dvertexenergyvolume10;
    double vertexE10 = (slice - prong3D - prong2D);
    if (vertexE10 > -0.075)
      return true;
    else
      return false;
  }
});

//********************************** Vertex Energy Faliure Check Cut *******************************************************************

const Cut kVertexECheckCut([](const caf::SRProxy *sr) {
  if (sr->vtx.elastic.IsValid != true)
    return false;
  else
  {
    double slice     = sr->vtx.elastic.slicevertexenergyvolume10;
    double prong3D   = sr->vtx.elastic.prong3dvertexenergyvolume10;
    double prong2D   = sr->vtx.elastic.prong2dvertexenergyvolume10;
    //double vertexE10 = (slice - prong3D - prong2D);
    if (isnan(slice)   ||
        isnan(prong3D) ||
        isnan(prong2D))
      return false;
    else
      return true;
  }
});


// const Cut isPion([](const caf::SRProxy *sr){
//   int muonNum=GetMuonProngId(sr);
//   if(muonNum<0)
//     return -5.0;
//   else{
//     int prongNum=0;
//     int maxVal=-1;
//     if(sr->vtx.nelastic==0 || sr->vtx.elastic.fuzzyk.npng==0)
//       return -5.0;

//     for(size_t i=0;i<sr->vtx.elastic.fuzzyk.npng;++i){
//       if(sr->vtx.elastic.fuzzyk.png[i].cvnpart.pionid>maxVal){
//         maxVal=sr->vtx.elastic.fuzzyk.png[i].cvnpart.pionid;
//         prongNum=i;
//       }
//     }
//     return (sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
//     //return double(sr->vtx.elastic.fuzzyk.png[prongNum].calE);
//   }
//   //if(sr->vtx.elastic.fuzzyk.png[prongNum].truth.pdg == 211)
//   //else if (! ())
// });




const Cut kQE([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  else
    return (sr->mc.nu[0].mode == caf::kQE);
});

const Cut kCOH([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  else
    return (sr->mc.nu[0].mode == caf::kCoh);
});

const Cut kMEC([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  else
    return (sr->mc.nu[0].mode == caf::kMEC);
});

const Cut kRES([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  else
    return (sr->mc.nu[0].mode == caf::kRes);
});

const Cut kDIS([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  else
    return (sr->mc.nu[0].mode == caf::kDIS);
});

const Cut kCTarget([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  if (sr->mc.nu[0].tgtZ == 6)
    return true;
  else
    return false;
});
const Cut kHTarget([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  if (sr->mc.nu[0].tgtZ == 1)
    return true;
  else
    return false;
});
const Cut kClTarget([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  if (sr->mc.nu[0].tgtZ == 17)
    return true;
  else
    return false;
});
const Cut kTiTarget([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  if (sr->mc.nu[0].tgtZ == 22)
    return true;
  else
    return false;
});

const Cut kOTarget([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  if (sr->mc.nu[0].tgtZ == 8)
    return true;
  else
    return false;
});

const Cut kNC([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  else
    if (sr->mc.nu[0].iscc == false)
      return true;
    else
      return false;
 });

