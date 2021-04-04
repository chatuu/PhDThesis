#include "headers.h"
#include "structs.h"

#include "cuts.h"
#include "vars.h"

#include <stdio.h>
//#include "switches.cxx"

using namespace ana;

void createSpectra()
{
    std::vector<Cuts> selectionCuts;
    selectionCuts.push_back({"VertexECut",        kNumuMyQuality && kIsFiducial && kNumuTightContainND && NuMuCCMuonIDCut && kTwoProng && kPionIdCut && kVertexECut});
    
    std::vector<Cuts> interactionCuts;
    interactionCuts.push_back({"RecoCohSig",    kNuMuCCCohSig});

    std::vector<Vars> vars;
    vars.push_back({"Est Pion K.E", {"Estimated Pion K.E (GeV)", Binning::Simple(3000, 0, 3), kPionKEEstNew}});
    vars.push_back({"True Pion K.E", {"True Pion K.E (GeV)", Binning::Simple(3000, 0, 3), kTruePionKENew}});
    vars.push_back({"NewPionEstRes", {"(Reco-True)/True", Binning::Simple(4000, -1, 1), kPionKEResNew}});
    vars.push_back({"NewPionBPFRes", {"(Reco-True)/True", Binning::Simple(4000, -1, 1), kPionKEResBPFNew}});
    vars.push_back({"Est muon K.E", {"Estimated muon K.E (GeV)", Binning::Simple(3000, 0, 3), kMuonProngEstKE}});
    vars.push_back({"True muon K.E", {"True muon K.E (GeV)", Binning::Simple(3000, 0, 3), kMuonProngTrueKENew}});
    vars.push_back({"NewmuonEstRes", {"(Reco-True)/True", Binning::Simple(4000, -1, 1), kMuonEResNew}});
    vars.push_back({"NewmuonBPFRes", {"(Reco-True)/True", Binning::Simple(4000, -1, 1), kMuonEResBPFNew}});
     
    // std::vector<Vars2D> vars2D;
    // vars2D.push_back({"Muon Prong Length (cm) Vs True Muon Kinetic Energy (GeV)", {"Muon Prong Length (cm)", Binning::Simple(1600, 0, 1600), kMuonProngLenNew}, {"True Muon K.E (GeV)", Binning::Simple(3000, 0, 3), kMuonProngTrueKENew}});
    // vars2D.push_back({"Pion CalE (GeV) Vs True Pion Kinetic Energy (GeV)", {"Pion Prong calE (GeV)", Binning::Simple(2000, 0, 2), kPionCalENew }, {"True Pion K.E (GeV)", Binning::Simple(2000, 0, 2), kTruePionKENew}});
    
    // std::vector<Spectra> spectra2D;
    std::vector<Spectra> spectra;

    SpectrumLoader lNDMC("prod_flatsumdecaf_development_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_ndphysics_contain_v1");

    lNDMC.SetSpillCut(kStandardSpillCuts);

    //1-D Spectra
    for (auto intit = interactionCuts.begin(); intit != interactionCuts.end(); ++intit)
        for (auto varit = vars.begin(); varit != vars.end(); ++varit)
            for (auto cutit = selectionCuts.begin(); cutit != selectionCuts.end(); ++cutit)
                spectra.push_back({intit->name, varit->name, cutit->name, new Spectrum(lNDMC, varit->axis, cutit->cut && intit->cut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2020)});


 
    //2-D Spectra
    // for (auto intit = interactionCuts.begin(); intit != interactionCuts.end(); ++intit)
    //     for (auto varit = vars2D.begin(); varit != vars2D.end(); ++varit)
    //         for (auto cutit = selectionCuts.begin(); cutit != selectionCuts.end(); ++cutit)
    //             spectra2D.push_back({intit->name, varit->name, cutit->name, new Spectrum(lNDMC, varit->Xaxis, varit->Yaxis, cutit->cut && intit->cut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2020)});

    lNDMC.Go();

    TFile *f = new TFile("Mar_28_2021_FinalizingEnergyEstimators_1D.root", "recreate");

    std::cout << "\n"
              << std::endl;

    printf("  %-18s %-16s %-22s %-5s \n", "Interaction Type", "Variable Name", "Cut Name", "POT");

    // for (auto it = spectra2D.begin(); it != spectra2D.end(); ++it)
    // {
    //     auto htemp = it->spectrum->ToTH2(it->spectrum->POT(), kPOT);
    //     htemp->SetName(Form("Interaction_%s_Variable_%s_Cut_%s", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str()));
    //     htemp->Write();
    //     it->spectrum->SaveTo(f,Form("%s_%s_%s", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str()));
    //     printf("  %-18s %-16s %-22s %.6e \n", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str(), it->spectrum->POT());
    // }
    for (auto it = spectra.begin(); it != spectra.end(); ++it)
    {
        it->spectrum->SaveTo(f,Form("Interaction_%s_Variable_%s_Cut_%s", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str()));
        printf("  %-18s %-16s %-22s %.6e \n", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str(), it->spectrum->POT());
    }

     f->Close();
}
    
