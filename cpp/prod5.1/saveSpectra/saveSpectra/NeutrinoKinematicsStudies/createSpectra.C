#include "../headers/headers.h"
#include "../headers/functions.h"
#include "../headers/structs.h"
#include "../headers/cuts.h"
#include "../headers/vars.h"

using namespace ana;

void createSpectra()
{
    std::vector<Cuts> selectionCuts;
    selectionCuts.push_back({"NewPionIDCut",      kNumuMyQuality && kIsFiducial && kNumuTightContainND && NuMuCCMuonIDCut && kTwoProng && kPionIdCut});
    selectionCuts.push_back({"VertexECut",        kNumuMyQuality && kIsFiducial && kNumuTightContainND && NuMuCCMuonIDCut && kTwoProng && kPionIdCut && kVertexECut});
    
    std::vector<Cuts> interactionCuts;
    interactionCuts.push_back({"RecoCohSig",    kNuMuCCCohSig});
    interactionCuts.push_back({"RecoCohBkgd",  !kNuMuCCCohSig});

    std::vector<Vars> vars;
    vars.push_back({"NeutrinoMass",	 {"M_{#nu_{#mu}} (GeV)", Binning::Simple(10000, -0.005, 0.005), kNuMuMass}});
    vars.push_back({"Truet", {"True |t| (GeV)", Binning::Simple(1000, 0, 1), kTruet}});
    //vars.push_back({"RecoTruet", {"Estimated True |t| (GeV)", Binning::Simple(1000, 0, 1), kRecoTruet}});
    vars.push_back({"RecoTruet", {"Estimated True |t| (GeV)", Binning::Simple(1000, 0, 1), kRecoTruet2}});
    
    //vars.push_back({"Estimated |t|", {"Reconstructed |t| (GeV)", Binning::Simple(1000, 0, 1), kEstt}});
    vars.push_back({"Estimated |t|", {"Reconstructed |t| (GeV)", Binning::Simple(1000, 0, 1), kEstt2}});

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

    TFile *f = new TFile("Apr_01_2021_NeutrinoMass.root", "recreate");

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
    
