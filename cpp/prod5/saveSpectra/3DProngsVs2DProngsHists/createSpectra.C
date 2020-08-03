#include "headers.h"
#include "structs.h"
#include "cuts.h"
#include "vars.h"
#include <stdio.h>

using namespace ana;

void createSpectra()
{

    std::vector<Cuts> selectionCuts;
    selectionCuts.push_back({"muonID", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut});

    std::vector<Cuts> interactionCuts;
    interactionCuts.push_back({"RecoCohSig", kNuMuCCCohSig});
    interactionCuts.push_back({"RecoCohBkgd", !kNuMuCCCohSig});

    std::vector<Vars2D> vars2D;
    vars2D.push_back({"Number of 3D Prongs and 2D Prongs", {"Number of 3D Prongs", Binning::Simple(10, 0, 10), NumberOfProngs}, {"Number of 2D Prongs", Binning::Simple(12, 0, 12), prongs2D}});

    std::vector<Spectra> spectra2D;

    //SpectrumLoader lNDMC("defname: prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");
    SpectrumLoader lNDMC("def_snapshot prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");

    lNDMC.SetSpillCut(kStandardSpillCuts);

    for (auto intit = interactionCuts.begin(); intit != interactionCuts.end(); ++intit)
        for (auto varit = vars2D.begin(); varit != vars2D.end(); ++varit)
            for (auto cutit = selectionCuts.begin(); cutit != selectionCuts.end(); ++cutit)
                spectra2D.push_back({intit->name, varit->name, cutit->name, new Spectrum(lNDMC, varit->Xaxis, varit->Yaxis, cutit->cut && intit->cut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2020)});

    lNDMC.Go();

    TFile *f = new TFile("Aug-02-2020-3DProngsVs2DProngsProd5.root", "recreate");

    std::cout << "\n"
              << std::endl;

    printf("  %-18s %-16s %-22s %-5s \n", "Interaction Type", "Variable Name", "Cut Name", "POT");

    for (auto it = spectra2D.begin(); it != spectra2D.end(); ++it)
    {
        auto htemp = it->spectrum->ToTH2(it->spectrum->POT(), kPOT);
        htemp->SetName(Form("Interaction_%s_Variable_%s_Cut_%s", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str()));
        htemp->Write();
        it->spectrum->SaveTo(f->mkdir(Form("%s_%s_%s", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str())));
        printf("  %-18s %-16s %-22s %.6e \n", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str(), it->spectrum->POT());
    }

    f->Close();
    std::cout << "\n"
              << std::endl;
    //Spectrum2D->toTH2(Spectrum2D->POT(), kPOT)
}
