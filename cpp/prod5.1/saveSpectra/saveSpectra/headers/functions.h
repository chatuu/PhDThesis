#include "headers.h"
#ifndef FUNCTIONS_H
#define FUNCTIONS_H
/********************************** Get NuMu FourMomenta ******************************************************/
TLorentzVector GetNeutrinoFourMomenta(const caf::SRProxy *sr)
{
  if (sr->mc.nnu != 0)
  {
    TLorentzVector nu;
    nu.SetPxPyPzE( sr->mc.nu[0].p.px,
                   sr->mc.nu[0].p.py,
                   sr->mc.nu[0].p.pz,
                   sr->mc.nu[0].p.E );

    return nu;
  }
}
/********************************** Get Muon FourMomenta ******************************************************/
TLorentzVector GetMuonFourMomenta(const caf::SRProxy *sr)
{
  if (sr->mc.nnu != 0)
  {
  
    TLorentzVector muon;

    for (unsigned int i = 0; i < sr->mc.nu[0].prim.size(); ++i)
    {
      if (sr->mc.nu[0].prim[i].pdg == 13)
      {
        muon.SetPxPyPzE( sr->mc.nu[0].prim[i].p.px,
                         sr->mc.nu[0].prim[i].p.py,
                         sr->mc.nu[0].prim[i].p.pz,
                         sr->mc.nu[0].prim[i].p.E );
      }
    }
    return muon;
  }
}
/********************************** Get Pion FourMomenta ******************************************************/
TLorentzVector GetPionFourMomenta(const caf::SRProxy *sr)
{
  if (sr->mc.nnu != 0)
  {
    TLorentzVector pion;

    for (unsigned int i = 0; i < sr->mc.nu[0].prim.size(); ++i)
    {
      if (sr->mc.nu[0].prim[i].pdg == 211)
      {
        pion.SetPxPyPzE( sr->mc.nu[0].prim[i].p.px,
                         sr->mc.nu[0].prim[i].p.py,
                         sr->mc.nu[0].prim[i].p.pz,
                         sr->mc.nu[0].prim[i].p.E );
      }
    }
    return pion;
  }
}
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
  int prongNum = -1;
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
      return prongNum;
    }
  }
  if (prongNum != -1)
    return prongNum;
  else
    return 10000;
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
    int prongNum = -1;
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
      if (prongNum != -1)
        return prongNum;
      else
        return 10000;
  }
}
/******************************************* New Muon Energy Estimator ***************************************************/
double MuonEnergyEstimator(float muonLen){
  double muonKE = (0.00206646*muonLen + 0.0201737);
  return muonKE;
} 

/******************************************* New Pion Energy Estimator ***************************************************/
/*
****************************************
Minimizer is Linear / Migrad

*/
double PionEnergyEstimator(float calE){
  //float p0                        =     0.249689; //  +/-   0.115564    
  //float p1                        =     -1.13424; //  +/-   3.2228      
  //float p2                        =     -2.50427; //  +/-   30.0966     
  //float p3                        =      62.1546; //  +/-   131.021     
  //float p4                        =     -217.197; //  +/-   310.894     
  //float p5                        =      361.806; //  +/-   433.586     
  //float p6                        =      -335.23; //  +/-   364.923     
  //float p7                        =      176.795; //  +/-   182.094     
  //float p8                        =      -49.726; //  +/-   49.5573     
  //float p9                        =      5.80004; //  +/-   5.66286  

  //float p0                        =     0.316376; //   +/-   0.039659    
  //float p1                        =      -4.1791; //   +/-   1.9122      
  //float p2                        =      33.1168; //   +/-   31.5049     
  //float p3                        =     -114.635; //   +/-   242.257     
  //float p4                        =      239.423; //   +/-   1014.62     
  //float p5                        =      -319.64; //   +/-   2496.1      
  //float p6                        =      286.998; //   +/-   3704.82     
  //float p7                        =       -177.3; //   +/-   3260.21     
  //float p8                        =      69.3147; //   +/-   1565.08     
  //float p9                        =      -12.443; //   +/-   315.578    

//  Minimizer is Linear / Migrad
//Chi2                      =     0.132106
//NDf                       =           95
//float p0                        =     0.141479; //   +/-   0.0542481   
//float p1                        =    -0.516404; //   +/-   0.513912    
//float p2                        =       4.7073; //   +/-   1.5677      
//float p3                        =     -5.52809; //   +/-   1.8926      
//float p4                        =      2.18037; //   +/-   0.784279  

//Minimizer is Linear / Migrad
//Chi2                      =    0.0347367
//NDf                       =           55
float p0                        =    0.0416097;  // +/-   0.0830038   
float p1                        =     0.831222;  // +/-   1.08518     
float p2                        =    -0.283639;  // +/-   4.75972     
float p3                        =      1.03192;  // +/-   8.48411     
float p4                        =     -0.45108;  // +/-   5.27915


	
	return (
    p0             + 
		p1*calE        + 
		p2*pow(calE,2) + 
		p3*pow(calE,3) +
		p4*pow(calE,4) );
		//p5*pow(calE,5) +
		//p6*pow(calE,6) +
		//p7*pow(calE,7) +
		//p8*pow(calE,8) +
		//p9*pow(calE,9) );
}


#endif
