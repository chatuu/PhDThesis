#include "headers.h"
#include "structs.h"
#include "cuts.h"
#include "vars.h"
#include <stdio.h>

using namespace ana;

void createSpectra()
{

    std::vector<Cuts> selectionCuts;
    //selectionCuts.push_back({"NoCut", NoCuts});
    //selectionCuts.push_back({"slicing", kNoCut});
    //selectionCuts.push_back({"dqcuts", kNumuMyQuality});
    //selectionCuts.push_back({"fiducial", kNumuMyQuality && kIsFiducial});
    selectionCuts.push_back({"containment", kNumuMyQuality && kIsFiducial && kNumuTightContainND});
    selectionCuts.push_back({"muonID", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut});
    selectionCuts.push_back({"TwoProng", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut && kTwoProng});

    //selectionCuts.push_back({"PionIdCut", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut && kTwoProng && ktCut && kPionIdCut});
    //selectionCuts.push_back({"PionIdCutNotCut", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut && kTwoProng && kPionIdCut});
    selectionCuts.push_back({"PionIdCutNotCut", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut && kTwoProng && kPionIdCut});
    //selectionCuts.push_back({"tCut", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut && kTwoProng && kPionIdCut && ktCut});
    //selectionCuts.push_back({"NewtCut", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut && kTwoProng && kPionIdCut && kRecotNoBPF});

    //selectionCuts.push_back({"tFaliure", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut && kTwoProng && kPionIdCut && ktCutFaliure});

    std::vector<Cuts> interactionCuts;
    //interactionCuts.push_back({"slicing",kNoCut});
    interactionCuts.push_back({"RecoCohSig", kNuMuCCCohSig});
    interactionCuts.push_back({"RecoCohBkgd", !kNuMuCCCohSig});
    //interactionCuts.push_back({"RecoCohBkgdQEOnly", !kNuMuCCCohSig && kQE});
    //interactionCuts.push_back({"RecoCohBkgdMECOnly", !kNuMuCCCohSig && kMEC});
    //interactionCuts.push_back({"RecoCohBkgdRESOnly", !kNuMuCCCohSig && kRES});
    //interactionCuts.push_back({"RecoCohBkgdDISOnly", !kNuMuCCCohSig && kDIS});

    //std::vector<TruthCuts> truthInteractionCuts;
    //truthInteractionCuts.push_back({"TruthCohSig", kTrueNuMuCCCohSig});
    //truthInteractionCuts.push_back({"TruthCohBkgd", !kTrueNuMuCCCohSig});

    std::vector<Vars> vars;
    //vars.push_back({"NeutrinoEnergy", {"Neutrino Energy", Binning::Simple(100, 0, 10), VarFromNuTruthVar(nutrinoEnergy)}});
    //vars.push_back({"NumberOfTracks",{"Number of Tracks",Binning::Simple(10,0,10),NumberOfTracks}});
    //vars.push_back({"NumberOfProngs",{"Number of Prongs",Binning::Simple(10,0,10),NumberOfProngs}});
    //vars.push_back({"calE",{"calE",Binning::Simple(140,0,14),calE}});
    //vars.push_back({"muonCalE", {"muonCalE", Binning::Simple(140, 0, 14), muonCalE}});
    //vars.push_back({"pionCalE", {"pionCalE", Binning::Simple(1000, 0, 10), pionCalE}});
    //vars.push_back({"vertexAct", {"vertexAct", Binning::Simple(200, -1, 1), vertexAct}});
    // vars.push_back({"ntracks2D",{"ntracks2D",Binning::Simple(12,0,12),track2D}});
    // vars.push_back({"nprongs2D",{"nprongs2D",Binning::Simple(12,0,12),track2D}});
    //vars.push_back({"Topology", {"Topology", Binning::Simple(32, 0, 32), Topology}});
    //vars.push_back({"tSquared", {"|t|", Binning::Simple(100, 0, 1), kt}});
    //vars.push_back({"EThetaSquared", {"EThetaSq", Binning::Simple(100, 0, 1), EthetaSq}});
    //vars.push_back({"EEta", {"EEta", Binning::Simple(100, 0, 1), RecoEeta}});
    //vars.push_back({"True Pion Kinetic Energy", {"True Pion Kinetic Energy (GeV)", Binning::Simple(1000, 0, 10), kTrueBPFPionE}});
    //vars.push_back({"True Muon Kinetic Energy", {"True Muon Kinetic Energy (GeV)", Binning::Simple(1000, 0, 10), kTrueBPFMuonE}});
    //vars.push_back({"Reconstructed Pion Kinetic Energy", {"Reconstructed Pion Kinetic Energy (GeV)", Binning::Simple(1000, 0, 10), kRecoKEPion}});
    //vars.push_back({"Reconstructed Muon Kinetic Energy", {"Reconstructed Muon Kinetic Energy (GeV)", Binning::Simple(1000, 0, 10), kRecoKEMuon}});
    //vars.push_back({"Reconstructed |t| without BPF", {"Reconstructed |t| (GeV)", Binning::Simple(1000, 0, 10), kRecotNoBPF}});
    //vars.push_back({"PionE Resolution", {"PionE Resolution", Binning::Simple(2000, -10, 10), kPionERes}});
    //vars.push_back({"muonE Resolution", {"muonE Resolution", Binning::Simple(2000, -10, 10), kMuonERes}});
 
    //vars.push_back({"muonPx Resolution", {"muonPx Resolution", Binning::Simple(2000, -10, 10), kMuonPxRes}});
    //vars.push_back({"muonPy Resolution", {"muonPy Resolution", Binning::Simple(2000, -10, 10), kMuonPyRes}});
    //vars.push_back({"muonPz Resolution", {"muonPz Resolution", Binning::Simple(2000, -10, 10), kMuonPzRes}});
    vars.push_back({"muonPt Resolution", {"muonPt Resolution", Binning::Simple(2000, -10, 10), kMuonPtRes}});

    //vars.push_back({"pionPx Resolution", {"pionPx Resolution", Binning::Simple(2000, -10, 10), kPionPxRes}});
    //vars.push_back({"pionPy Resolution", {"pionPy Resolution", Binning::Simple(2000, -10, 10), kPionPyRes}});
    //vars.push_back({"pionPz Resolution", {"pionPz Resolution", Binning::Simple(2000, -10, 10), kPionPzRes}});
    vars.push_back({"pionPt Resolution", {"pionPt Resolution", Binning::Simple(2000, -10, 10), kPionPtRes}});
 
    //vars.push_back({"Reconstructed |t| without BPF", {"Reconstructed |t| (GeV)", Binning::Simple(1000, 0, 10), kRecotNoBPF}});

    //vars.push_back({"Pion Prong ID", {"Pion Prong ID", Binning::Simple(100, 0, 1), pionProngId}});
    //vars.push_back({"Proton Prong ID", {"Proton Prong ID", Binning::Simple(100, 0, 1), protonProngId}});
    //vars.push_back({"vtxX", {"vtx X", Binning::Simple(800, -200, 200), vtxX}});
    //vars.push_back({"vtxY", {"vtx Y", Binning::Simple(800, -200, 200), vtxY}});
    //vars.push_back({"vtxZ", {"vtx Z", Binning::Simple(3400, -100, 1700), vtxZ}});
    //vars.push_back({,{}});
    //std::vector<Vars2D> vars2D;
    //vars2D.push_back({"Number of Tracks and 2D tracks",{"Number of Tracks",Binning::Simple(10,0,10),NumberOfTracks},{"ntracks2D",Binning::Simple(12,0,12),track2D}});
    //vars2D.push_back({"Number of 3D Prongs and 2D Prongs", {"Number of 3D Prongs", Binning::Simple(10, 0, 10), NumberOfProngs}, {"Number of 2D Prongs", Binning::Simple(12, 0, 12), prongs2D}});
    //vars2D.push_back({"Topology Vs Pion Kinetic Energy", {"Pion Energy (GeV)", Binning::Simple(1000, 0, 10), kPionEnergy}, {"Topology", Binning::Simple(32, 0, 32), Topology}});
    //vars2D.push_back({"Topology Vs Particle Kinetic Energy", {"Energy (GeV)", Binning::Simple(1000, 0, 10), mostEnergetic}, {"Topology", Binning::Simple(32, 0, 32), Topology}});
    //vars2D.push_back({"thetaSq Vs TruePion Kinetic Energy", {"Energy (GeV)", Binning::Simple(1000, 0, 10), pionCalE}, {"RadSquared", Binning::Simple(1000, 0, 10), ProngPionAngleSq}});
    //vars.push_back({"dEdX_Log_likelihood", {"Pion dE/dX log likelihood", Binning::Simple(100, -1, 1), pionLLDEDX}});
    //vars.push_back({"dEdX_Log_likelihood_Kalman", {"Pion dE/dX log likelihood", Binning::Simple(100, -3, 3), kalmandEdXllh}});
    //vars.push_back({"dE/dX Log likelihood", {"Muon dE/dX log likelihood", Binning::Simple(100,-1,1), muonLLDEDX}});
    //vars.push_back({"Reconstructed |t|", {"Reconstructed |t| (GeV)", Binning::Simple(1000, 0, 10), kRecot}});
    //vars.push_back({"Pt", {"Pt (GeV)", Binning::Simple(2000, -10, 10), kRecoPt}});

    //std::vector<TruthVars> trueVars;
    //trueVars.push_back({"NeutrinoEnergy", {"Neutrino Energy", Binning::Simple(100, 0, 10), nutrinoEnergy}});
    //trueVars.push_back({"tSquared", {"|t|", Binning::Simple(100, 0, 1), kt}});
    //trueVars.push_back({"EThetaSquared", {"EThetaSq", Binning::Simple(100, 0, 1), TrueEthetaSq}});

    std::vector<Spectra> spectra;
    //std::vector<Spectra> spectra2D;

    SpectrumLoader lNDMC("defname: prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");
    //SpectrumLoader lNDMC("def_snapshot prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");
    //SpectrumLoader lNDMC("defname: prod_caf_R17-11-14-prod4reco.CVNprong-respin.b_nd_genie_nonswap_fhc_nova_v08_full_v1 with limit 10");
    //SpectrumLoader lNDMC("def_snapshot prod_caf_R17-11-14-prod4reco.CVNprong-respin.b_nd_genie_nonswap_fhc_nova_v08_full_v1");

    lNDMC.SetSpillCut(kStandardSpillCuts);
    //1-D Spectra
    for (auto intit = interactionCuts.begin(); intit != interactionCuts.end(); ++intit)
        for (auto varit = vars.begin(); varit != vars.end(); ++varit)
            for (auto cutit = selectionCuts.begin(); cutit != selectionCuts.end(); ++cutit)
                spectra.push_back({intit->name, varit->name, cutit->name, new Spectrum(lNDMC, varit->axis, cutit->cut && intit->cut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2020)});

    // for (auto intit = interactionCuts.begin(); intit != interactionCuts.end(); ++intit)
    //     for (auto varit = vars2D.begin(); varit != vars2D.end(); ++varit)
    //         for (auto cutit = selectionCuts.begin(); cutit != selectionCuts.end(); ++cutit)
    //             spectra2D.push_back({intit->name, varit->name, cutit->name, new Spectrum(lNDMC, varit->Xaxis, varit->Yaxis, cutit->cut && intit->cut, kNoShift, kPPFXFluxCVWgt * kXSecCVWGT2020)});

    //1 - D Truth Variables Spectra
    //    for (auto intit = truthInteractionCuts.begin(); intit != truthInteractionCuts.end(); ++intit)
    //        for (auto varit = trueVars.begin(); varit != trueVars.end(); ++varit)
    //            for (auto varit = trueVars.begin(); varit != trueVars.end(); ++varit)
    //                spectra.push_back({intit->name, varit->name, "_Truth_", new Spectrum(lNDMC, varit->axis, intit->cut, kNoShift, kPPFXFluxCVWgtST * kXSecCVWgt2020_NT)});

    lNDMC.Go();
    // const std::string fname = "defname: prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1";

    // Var kRun = SIMPLEVAR(hdr.run);
    // Var kSubrun = SIMPLEVAR(hdr.subrun);
    // Var kCycle = SIMPLEVAR(hdr.cycle);
    // Var kEvt = SIMPLEVAR(hdr.evt);
    // Var kSlice = SIMPLEVAR(hdr.subevt);
    // Var kNuE = VarFromNuTruthVar(nutrinoEnergy);

    // std::vector<const Var *> Newvars = {&kRun, &kSubrun, &kCycle, &kEvt, &kSlice, &kNuE, &vtxX, &vtxY, &vtxZ, &Topology};

    // MakeTextListFile(fname, {kNuMuCCCohSig && kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut && kTwoProng && ktCutFaliure}, {"SignaleventList.txt"}, Newvars, &kStandardSpillCuts);

    TFile *f = new TFile("Aug_14_2020_MomentumResolutionPlots.root", "recreate");

    std::cout << "\n"
              << std::endl;

    printf("  %-18s %-16s %-22s %-5s \n", "Interaction Type", "Variable Name", "Cut Name", "POT");

    for (auto it = spectra.begin(); it != spectra.end(); ++it)
    {
        it->spectrum->SaveTo(f->mkdir(Form("Interaction_%s_Variable_%s_Cut_%s", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str())));
        printf("  %-18s %-16s %-22s %.6e \n", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str(), it->spectrum->POT());
    }

    // for (auto it = spectra2D.begin(); it != spectra2D.end(); ++it)
    // {
    //     auto htemp = it->spectrum->ToTH2(it->spectrum->POT(), kPOT);
    //     htemp->SetName(Form("Interaction_%s_Variable_%s_Cut_%s", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str()));
    //     htemp->Write();
    //     it->spectrum->SaveTo(f->mkdir(Form("%s_%s_%s", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str())));
    //     printf("  %-18s %-16s %-22s %.6e \n", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str(), it->spectrum->POT());
    // }

    f->Close();
    std::cout << "\n"
              << std::endl;
    //Spectrum2D->toTH2( Spectrum2D->POT(), kPOT )
}
