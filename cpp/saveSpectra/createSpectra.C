#include "headers.h"
#include "structs.h"
#include "cuts.h"
#include "vars.h"

using namespace ana;

void createSpectra()
{
    std::vector<Cuts> selectionCuts;
    selectionCuts.push_back({"slicing", kNoCut});
    selectionCuts.push_back({"dqcuts", kNumuMyQuality});
    selectionCuts.push_back({"fiducial", kNumuMyQuality && kIsFiducial});
    selectionCuts.push_back({"containment", kNumuMyQuality && kIsFiducial && kNumuTightContainND});
    selectionCuts.push_back({"muonID", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut});

    std::vector<Cuts> interactionCuts;
    interactionCuts.push_back({"slicing",kNoCut});
    interactionCuts.push_back({"RecoCohSig", kNuMuCCCohSig});
    interactionCuts.push_back({"RecoCohBkgd", !kNuMuCCCohSig});
    
    std::vector<TruthCuts> truthInteractionCuts;
    truthInteractionCuts.push_back({"TruthCohSig", kTrueNuMuCCCohSig});
    truthInteractionCuts.push_back({"TruthCohBkgd", !kTrueNuMuCCCohSig});

    std::vector<Vars> vars;
    //vars.push_back({"NeutrinoEnergy", {"Neutrino Energy", Binning::Simple(100, 0, 10), VarFromNuTruthVar(nutrinoEnergy)}});
    vars.push_back({"NumberOfTracks",{"Number of Tracks",Binning::Simple(10,0,10),NumberOfTracks}});
    vars.push_back({"NumberOfProngs",{"Number of Prongs",Binning::Simple(10,0,10),NumberOfProngs}});
    vars.push_back({"calE",{"calE",Binning::Simple(140,0,14),calE}});
    vars.push_back({"muonCalE",{"muonCalE",Binning::Simple(140,0,14),muonCalE}});
    vars.push_back({"pionCalE",{"pionCalE",Binning::Simple(140,0,14),pionCalE}});
    vars.push_back({"vertexAct",{"vertexAct",Binning::Simple(200,-10,10),vertexAct}});

    std::vector<TruthVars> trueVars;
    trueVars.push_back({"NeutrinoEnergy", {"Neutrino Energy", Binning::Simple(100, 0, 10),nutrinoEnergy}});

    std::vector<Spectra> spectra;

    SpectrumLoader lNDMC("def_snapshot prod_caf_R17-11-14-prod4reco.CVNprong-respin.b_nd_genie_nonswap_fhc_nova_v08_full_v1");
    //SpectrumLoader lNDMC("defname: prod_caf_R17-11-14-prod4reco.CVNprong-respin.b_nd_genie_nonswap_fhc_nova_v08_full_v1 with limit 10");
    //SpectrumLoader lNDMC("defname: prod_caf_R17-11-14-prod4reco.CVNprong-respin.b_nd_genie_nonswap_fhc_nova_v08_full_v1");

    lNDMC.SetSpillCut(kStandardSpillCuts);

    for(auto intit=interactionCuts.begin();intit!=interactionCuts.end();++intit)
        for (auto varit = vars.begin(); varit != vars.end(); ++varit)
            for (auto cutit = selectionCuts.begin(); cutit != selectionCuts.end(); ++cutit)
                spectra.push_back({intit->name,varit->name, cutit->name, new Spectrum(lNDMC, varit->axis, cutit->cut && intit->cut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2018)});

    // for(auto intit=truthInteractionCuts.begin();intit!=truthInteractionCuts.end();++intit)
    //     for (auto varit = trueVars.begin(); varit != trueVars.end(); ++varit)
    //         spectra.push_back({intit->name,varit->name, "Truth", new Spectrum(lNDMC, varit->axis, intit->cut, kNoShift, kPPFXFluxCVWgtST * kXSecCVWgt2018_smallerDISScale_NT)});

    lNDMC.Go();

    TFile *f = new TFile("TestNewFileTrackVsProng.root", "recreate");
    
    std::cout<<"\n"<<std::endl;

    printf("  %-18s %-16s %-22s %-5s \n","Interaction Type","Variable Name","Cut Name","POT");

    for (auto it = spectra.begin(); it != spectra.end(); ++it){
        it->spectrum->SaveTo(f->mkdir(Form("Interaction_%s_Variable_%s_Cut_%s",it->intName.c_str() ,it->varName.c_str(), it->cutName.c_str())));
        printf("  %-18s %-16s %-22s %.6e \n",it->intName.c_str(),it->varName.c_str(),it->cutName.c_str(),it->spectrum->POT());
    }

    f->Close();
    std::cout<<"\n"<<std::endl;
}
