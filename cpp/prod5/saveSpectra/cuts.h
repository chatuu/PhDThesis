#include "headers.h"

const Cut kMyNuMuCCCoh([](const caf::SRProxy *sr) {
  if (sr->mc.nnu == 0)
    return false;
  return (sr->mc.nu[0].iscc &&
          (sr->mc.nu[0].pdg == 14) &&
          (sr->mc.nu[0].pdgorig == 14) &&
          (sr->mc.nu[0].mode == caf::kCoh));
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

//********************************** Replicating NuMuCC inclusive cuts ****************************************

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
