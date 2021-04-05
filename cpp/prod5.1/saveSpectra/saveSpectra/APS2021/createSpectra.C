#include "headers.h"
#include "functions.h"
#include "structs.h"
#include "cuts.h"
#include "vars.h"

using namespace ana;

void createSpectra()
{
    const std::string fname = "prod_flatsumdecaf_development_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_ndphysics_contain_v1";

    Var kRun    = SIMPLEVAR(hdr.run);
    Var kSubrun = SIMPLEVAR(hdr.subrun);
    Var kCycle  = SIMPLEVAR(hdr.cycle);
    Var kEvt    = SIMPLEVAR(hdr.evt);
    Var kSlice  = SIMPLEVAR(hdr.subevt);
    Var kNuE    = VarFromNuTruthVar(nutrinoEnergy);

    std::vector<const Var *> Newvars = {&kRun, &kSubrun, &kCycle, &kEvt, &kSlice, &kNuE, &vtxX, &vtxY, &vtxZ, &prong3DVertexEVol10, &vertexE10};

    MakeTextListFile(fname, {kNuMuCCCohSig && kNumuMyQuality && kIsFiducial && kNumuTightContainND && NuMuCCMuonIDCut && kTwoProng && kPionIdCut}, {"CCCoherent.txt"}, Newvars, &kStandardSpillCuts);
    MakeTextListFile(fname, {kNuMuCCCohBack && kRES && kNumuMyQuality && kIsFiducial && kNumuTightContainND && NuMuCCMuonIDCut && kTwoProng && kPionIdCut}, {"CCRES.txt"}, Newvars, &kStandardSpillCuts);
    MakeTextListFile(fname, {kNuMuCCCohBack && kDIS && kNumuMyQuality && kIsFiducial && kNumuTightContainND && NuMuCCMuonIDCut && kTwoProng && kPionIdCut}, {"CCDIS.txt"}, Newvars, &kStandardSpillCuts);

}
    
