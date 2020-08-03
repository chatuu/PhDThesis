#include "headers.h"
#include "structs.h"
#include "cuts.h"
#include "vars.h"
#include <stdio.h>

using namespace ana;

void createSpectra()
{

    std::vector<Cuts> selectionCuts;
    selectionCuts.push_back({"slicing", kNoCut});
    selectionCuts.push_back({"dqcuts", kNumuMyQuality});
    selectionCuts.push_back({"fiducial", kNumuMyQuality && kIsFiducial});
    selectionCuts.push_back({"containment", kNumuMyQuality && kIsFiducial && kNumuTightContainND});
    selectionCuts.push_back({"muonID", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut});
    selectionCuts.push_back({"TwoProng", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut && kTwoProng});

    std::vector<Cuts> interactionCuts;
    interactionCuts.push_back({"RecoCohSig", kNuMuCCCohSig});
    interactionCuts.push_back({"RecoCohBkgd", !kNuMuCCCohSig});

    std::vector<Vars> vars;
    vars.push_back({"NeutrinoEnergy", {"Neutrino Energy", Binning::Simple(100, 0, 10), VarFromNuTruthVar(nutrinoEnergy)}});

    std::vector<Spectra> spectra;

    SpectrumLoader lNDMC("defname: prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");
    //SpectrumLoader lNDMC("def_snapshot prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");

    lNDMC.SetSpillCut(kStandardSpillCuts);

    //1-D Spectra
    for (auto intit = interactionCuts.begin(); intit != interactionCuts.end(); ++intit)
        for (auto varit = vars.begin(); varit != vars.end(); ++varit)
            for (auto cutit = selectionCuts.begin(); cutit != selectionCuts.end(); ++cutit)
                spectra.push_back({intit->name, varit->name, cutit->name, new Spectrum(lNDMC, varit->axis, cutit->cut && intit->cut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2020)});

    lNDMC.Go();
    const std::string fname = "defname: prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1";
    //const std::string fname = "def_snapshot prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1";

    Var kRun = SIMPLEVAR(hdr.run);
    Var kSubrun = SIMPLEVAR(hdr.subrun);
    Var kCycle = SIMPLEVAR(hdr.cycle);
    Var kEvt = SIMPLEVAR(hdr.evt);
    Var kSlice = SIMPLEVAR(hdr.subevt);
    Var kNuE = VarFromNuTruthVar(nutrinoEnergy);

    std::vector<const Var *> Newvars = {&kRun, &kSubrun, &kCycle, &kEvt, &kSlice, &kNuE, &vtxX, &vtxY, &vtxZ, &Topology};

    MakeTextListFile(fname, {kNuMuCCCohSig && kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut && kTwoProng && ktCutFaliure}, {"SignalEventListAfterTwo3DProngCut.txt"}, Newvars, &kStandardSpillCuts);

    TFile *f = new TFile("Aug_02_2020.root", "recreate");

    std::cout << "\n"
              << std::endl;

    printf("  %-18s %-16s %-22s %-5s \n", "Interaction Type", "Variable Name", "Cut Name", "POT");

    for (auto it = spectra.begin(); it != spectra.end(); ++it)
    {
        it->spectrum->SaveTo(f->mkdir(Form("Interaction_%s_Variable_%s_Cut_%s", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str())));
        printf("  %-18s %-16s %-22s %.6e \n", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str(), it->spectrum->POT());
    }

    f->Close();
    std::cout << "\n"
              << std::endl;
}
