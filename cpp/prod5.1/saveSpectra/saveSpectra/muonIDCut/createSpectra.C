#include "headers.h"
#include "structs.h"

#include "cuts.h"
#include "vars.h"

#include <stdio.h>
#include "switches.cxx"

using namespace ana;

void createSpectra()
{

    std::vector<Cuts> selectionCuts;
    //selectionCuts.push_back({"NoCut",           NoCuts});
    //selectionCuts.push_back({"slicing",         kNoCut});
    //selectionCuts.push_back({"dqcuts",          kNumuMyQuality});
    //selectionCuts.push_back({"fiducial",        kNumuMyQuality && kIsFiducial});
    //selectionCuts.push_back({"containment",	      kNumuMyQuality && kIsFiducial && kNumuTightContainND});
    //selectionCuts.push_back({"muonID",          kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut});
    //selectionCuts.push_back({"muonIDCVN",         kNumuMyQuality && kIsFiducial && kNumuTightContainND && muonIDCutCVNBased});
    //selectionCuts.push_back({"muonIDNuMuCC",      kNumuMyQuality && kIsFiducial && kNumuTightContainND && NuMuCCMuonIDCut});
    //selectionCuts.push_back({"TwoProng",        kNumuMyQuality && kIsFiducial && kNumuTightContainND && muonIDCutCVNBased && kTwoProng});
    selectionCuts.push_back({"TwoProngNuMuCC",    kNumuMyQuality && kIsFiducial && kNumuTightContainND && NuMuCCMuonIDCut && kTwoProng});
    selectionCuts.push_back({"TwoProngNoMuonID",        kNumuMyQuality && kIsFiducial && kNumuTightContainND  && kTwoProng});
    //selectionCuts.push_back({"OldPionIdCut",       kNumuMyQuality && kIsFiducial && kNumuTightContainND && NuMuCCMuonIDCut && kTwoProng && kPionIdOldCut});
    selectionCuts.push_back({"NewPionIdCut",       kNumuMyQuality && kIsFiducial && kNumuTightContainND && muonIDCutCVNBased && kTwoProng && kPionIdCut});
    //selectionCuts.push_back({"VertexECut",      kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut && kTwoProng && ktCut && kPionIdCut && kVertexECut});
    //selectionCuts.push_back({"VertexECheckCut", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut && kTwoProng && ktCut && kPionIdCut && kVertexECut && kVertexECheckCut});
    //selectionCuts.push_back({"VertexECheckCut", kNumuMyQuality && kIsFiducial && kNumuTightContainND && kMuonIDCut && kTwoProng && kVertexECheckCut});
    
    std::vector<Cuts> interactionCuts;
    //interactionCuts.push_back({"slicing",kNoCut});
    interactionCuts.push_back({"RecoCohSig",    kNuMuCCCohSig});
    interactionCuts.push_back({"RecoCohBkgd",  !kNuMuCCCohSig});
    //interactionCuts.push_back({"RecoCohBkgdQEOnly", !kNuMuCCCohSig && kQE});
    //interactionCuts.push_back({"RecoCohBkgdMECOnly", !kNuMuCCCohSig && kMEC});
    //interactionCuts.push_back({"RecoCohBkgdRESOnly", !kNuMuCCCohSig && kRES});
    //interactionCuts.push_back({"RecoCohBkgdDISOnly", !kNuMuCCCohSig && kDIS});
    //interactionCuts.push_back({"RecoNuMuCCInts", NuMuCCInts});


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
    //vars.push_back({"EThetaSquared", {"EThetaSq", Binning::Simple(100, 0, 1), EthetaSq}});
    //vars.push_back({"tSquared", {"|t|", Binning::Simple(100, 0, 1), kt}});
    vars.push_back({"Truet", {"True |t| (GeV)", Binning::Simple(1000, 0, 1), kTruet}});
    vars.push_back({"RecoTruet", {"Estimated True |t| (GeV)", Binning::Simple(1000, 0, 1), kRecoTruet}});
    //vars.push_back({"EEta", {"EEta", Binning::Simple(100, 0, 1), RecoEeta}});
    //vars.push_back({"True Pion Kinetic Energy", {"True Pion Kinetic Energy (GeV)", Binning::Simple(1000, 0, 10), kTrueBPFPionE}});
    //vars.push_back({"True Muon Kinetic Energy", {"True Muon Kinetic Energy (GeV)", Binning::Simple(1000, 0, 10), kTrueBPFMuonE}});
    //vars.push_back({"Reconstructed Pion Kinetic Energy", {"Reconstructed Pion Kinetic Energy (GeV)", Binning::Simple(1000, 0, 10), kRecoKEPion}});
    //vars.push_back({"Reconstructed Muon Kinetic Energy", {"Reconstructed Muon Kinetic Energy (GeV)", Binning::Simple(1000, 0, 10), kRecoKEMuon}});
    //vars.push_back({"Reconstructed |t| without BPF", {"Reconstructed |t| (GeV)", Binning::Simple(1000, 0, 10), kRecotNoBPF}});
    vars.push_back({"Estimated |t|", {"Reconstructed |t| (GeV)", Binning::Simple(1000, 0, 1), kEstt}});

    vars.push_back({"Pion Prong ID", {"Pion Prong ID", Binning::Simple(100, 0, 1), newPionID}});
    //vars.push_back({"Pion Prong ID", {"Pion Prong ID", Binning::Simple(100, 0, 1), pionProngId}});
    vars.push_back({"VertexE10", {"Vertex Energy 10", Binning::Simple(6000, -0.3, 0.3), vertexE10}});
    vars.push_back({"ProngVertexE10", {"Prong 3D Vertex Energy 10", Binning::Simple(6000, 0, 1), prong3DVertexEVol10}});
    //vars.push_back({"SelectedVertexE10", {"Cleaned Vertex Energy 10", Binning::Simple(4000, -2, 2), SelectedVertexE10}});
    
    //vars.push_back({"Est Pion K.E", {"Estimated Pion K.E (GeV)", Binning::Simple(3000, 0, 3), kPionKEEstNew}});
    //vars.push_back({"True Pion K.E", {"True Pion K.E (GeV)", Binning::Simple(3000, 0, 3), kTruePionKENew}});
    //vars.push_back({"NewPionEstRes", {"(Reco-True)/True", Binning::Simple(4000, -1, 1), kPionKEResNew}});
    //vars.push_back({"NewPionBPFRes", {"(Reco-True)/True", Binning::Simple(4000, -1, 1), kPionKEResBPFNew}});
    //vars.push_back({"True Muon K.E", {"True Muon K.E (GeV)", Binning::Simple(3000, 0, 3), kMuonProngTrueKE}});
    
    //vars.push_back({"True Muon Prong Length", {"Muon Prong Length (cm)", Binning::Simple(1600, 0, 1600), kMuonProngLen}});
    //vars.push_back({"True Pion Prong Length", {"Pion Prong Length (cm)", Binning::Simple(1600, 0, 1600), kPionProngLen}});




    //vars.push_back({"Proton Prong ID", {"Proton Prong ID", Binning::Simple(100, 0, 1), protonProngId}});
    vars.push_back({"muonID NuMuCC Inc", {"muonID NumuCC Inc", Binning::Simple(200, -1, 1), NuMuCCMuonID}});
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
    //vars2D.push_back({"Muon Prong Length (cm) Vs True Muon Kinetic Energy (GeV)", {"Muon Prong Length (cm)", Binning::Simple(1600, 0, 1600), kMuonProngLen}, {"True Muon K.E (GeV)", Binning::Simple(3000, 0, 3), kMuonProngTrueKE}});
    //vars2D.push_back({"Pion Prong Length (cm) Vs True Pion Kinetic Energy (GeV)", {"Pion Prong Length (cm)", Binning::Simple(1600, 0, 1600), kPionProngLen}, {"True Pion K.E (GeV)", Binning::Simple(3000, 0, 3), kPionProngTrueKE}});
    //vars2D.push_back({"Pion CalE (GeV) Vs True Pion Kinetic Energy (GeV)", {"Pion Prong calE (GeV)", Binning::Simple(300, 0, 3), pionCalE }, {"True Pion K.E (GeV)", Binning::Simple(300, 0, 3), kPionProngTrueKE}});
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
    // std::vector<Spectra> spectra2D;

    

    SwitchInfo switch1;
    SwitchInfo switch2;
    SwitchInfo switch3;
    SwitchInfo switch4;
    SwitchInfo switch5;
    SwitchInfo switch6;
    SwitchInfo switch7;
    SwitchInfo switch8;
    SwitchInfo switch9;
    SwitchInfo switch10;
    SwitchInfo switch11;
    SwitchInfo switch12;
    SwitchInfo switch13;
    SwitchInfo switch14;


    switch1.cut         =  "TwoProngNoMuonID";
    switch1.variable    =  "muonID NuMuCC Inc";
    switch1.interaction =  "RecoCohSig";

    switch2.cut         =  "TwoProngNoMuonID";
    switch2.variable    =  "muonID NuMuCC Inc";
    switch2.interaction =  "RecoCohBkgd";

    switch3.cut         =  "TwoProngNuMuCC";
    switch3.variable    =  "Pion Prong ID";
    switch3.interaction =  "RecoCohSig";

    switch4.cut         =  "TwoProngNuMuCC";
    switch4.variable    =  "Pion Prong ID";
    switch4.interaction =  "RecoCohBkgd";

    switch5.cut         =  "NewPionIdCut";
    switch5.variable    =  "VertexE10";
    switch5.interaction =  "RecoCohSig";

    switch6.cut         =  "NewPionIdCut";
    switch6.variable    =  "VertexE10";
    switch6.interaction =  "RecoCohBkgd";

    switch7.cut         =  "NewPionIdCut";
    switch7.variable    =  "ProngVertexE10";
    switch7.interaction =  "RecoCohSig";

    switch8.cut         =  "NewPionIdCut";
    switch8.variable    =  "ProngVertexE10";
    switch8.interaction =  "RecoCohBkgd";

    switch9.cut         =  "NewPionIdCut";
    switch9.variable    =  "Truet";
    switch9.interaction =  "RecoCohSig";

    switch10.cut         =  "NewPionIdCut";
    switch10.variable    =  "Truet";
    switch10.interaction =  "RecoCohBkgd";

    switch11.cut         =  "NewPionIdCut";
    switch11.variable    =  "RecoTruet";
    switch11.interaction =  "RecoCohSig";

    switch12.cut         =  "NewPionIdCut";
    switch12.variable    =  "RecoTruet";
    switch12.interaction =  "RecoCohBkgd";

    switch13.cut         =  "NewPionIdCut";
    switch13.variable    =  "Estimated |t|";
    switch13.interaction =  "RecoCohSig";

    switch14.cut         =  "NewPionIdCut";
    switch14.variable    =  "Estimated |t|";
    switch14.interaction =  "RecoCohBkgd";



    switches switchingList;
    switchingList.InsertSwitch(switch1);
    switchingList.InsertSwitch(switch2);
    switchingList.InsertSwitch(switch3);
    switchingList.InsertSwitch(switch4);
    switchingList.InsertSwitch(switch5);
    switchingList.InsertSwitch(switch6);
    switchingList.InsertSwitch(switch7);
    switchingList.InsertSwitch(switch8);
    switchingList.InsertSwitch(switch9);
    switchingList.InsertSwitch(switch10);
    switchingList.InsertSwitch(switch11);
    switchingList.InsertSwitch(switch12);
    switchingList.InsertSwitch(switch13);
    switchingList.InsertSwitch(switch14);



    //SpectrumLoader lNDMC("defname: prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");
    //SpectrumLoader lNDMC("prod_caf_R20-10-06-miniprod5.1reco.b_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");
    //SpectrumLoader lNDMC("prod_caf_R20-11-25-prod5.1reco.a_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");
    SpectrumLoader lNDMC("prod_flatsumdecaf_development_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_ndphysics_contain_v1");
    //SpectrumLoader lNDMC("def_snapshot prod_caf_R19-11-18-prod5reco.d.h.l_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1");
    //SpectrumLoader lNDMC("defname: prod_caf_R17-11-14-prod4reco.CVNprong-respin.b_nd_genie_nonswap_fhc_nova_v08_full_v1 with limit 10");
    //SpectrumLoader lNDMC("def_snapshot prod_caf_R17-11-14-prod4reco.CVNprong-respin.b_nd_genie_nonswap_fhc_nova_v08_full_v1");

    
    

    lNDMC.SetSpillCut(kStandardSpillCuts);
    //1-D Spectra
    for (auto intit = interactionCuts.begin(); intit != interactionCuts.end(); ++intit)
        for (auto varit = vars.begin(); varit != vars.end(); ++varit)
            for (auto cutit = selectionCuts.begin(); cutit != selectionCuts.end(); ++cutit)
                if (switchingList.checkToSaveSpectrum(intit->name+varit->name+cutit->name))
                {
                    spectra.push_back({intit->name, varit->name, cutit->name, new Spectrum(lNDMC, varit->axis, cutit->cut && intit->cut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2020)});
                    //std::cout<<"Expected: " << intit->name + varit->name + cutit->name << "\n";
                }

    // //2-D Spectra
    // for (auto intit = interactionCuts.begin(); intit != interactionCuts.end(); ++intit)
    //     for (auto varit = vars2D.begin(); varit != vars2D.end(); ++varit)
    //         for (auto cutit = selectionCuts.begin(); cutit != selectionCuts.end(); ++cutit)
    //             spectra2D.push_back({intit->name, varit->name, cutit->name, new Spectrum(lNDMC, varit->Xaxis, varit->Yaxis, cutit->cut && intit->cut, kNoShift, kPPFXFluxCVWgt * kXSecCVWgt2020)});




    // 1 - D Truth Variables Spectra
    //    for (auto intit = truthInteractionCuts.begin(); intit != truthInteractionCuts.end(); ++intit)
    //        for (auto varit = trueVars.begin(); varit != trueVars.end(); ++varit)
    //            for (auto varit = trueVars.begin(); varit != trueVars.end(); ++varit)
    //                spectra.push_back({intit->name, varit->name, "_Truth_", new Spectrum(lNDMC, varit->axis, intit->cut, kNoShift, kPPFXFluxCVWgtST * kXSecCVWgt2020_NT)});

    lNDMC.Go();
    // const std::string fname = "prod_caf_R20-11-25-prod5.1reco.a_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1";

    // Var kRun    = SIMPLEVAR(hdr.run);
    // Var kSubrun = SIMPLEVAR(hdr.subrun);
    // Var kCycle  = SIMPLEVAR(hdr.cycle);
    // Var kEvt    = SIMPLEVAR(hdr.evt);
    // Var kSlice  = SIMPLEVAR(hdr.subevt);
    // Var kNuE    = VarFromNuTruthVar(nutrinoEnergy);

    // std::vector<const Var *> Newvars = {&kRun, &kSubrun, &kCycle, &kEvt, &kSlice, &kNuE, &vtxX, &vtxY, &vtxZ, &muonPID, &pionPID, &oldPionID, &newPionID};

    // MakeTextListFile(fname, {kNuMuCCCohSig && kNumuMyQuality && kIsFiducial && kNumuTightContainND && NuMuCCMuonIDCut && kTwoProng}, {"SignaleventList.txt"}, Newvars, &kStandardSpillCuts);

    TFile *f = new TFile("Mar_24_2021_FinalizingCuts.root", "recreate");

    std::cout << "\n"
              << std::endl;

    printf("  %-18s %-16s %-22s %-5s \n", "Interaction Type", "Variable Name", "Cut Name", "POT");

    for (auto it = spectra.begin(); it != spectra.end(); ++it)
    {
        it->spectrum->SaveTo(f,Form("Interaction_%s_Variable_%s_Cut_%s", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str()));
        printf("%-18s %-16s %-22s %.6e \n", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str(), it->spectrum->POT());
    }

    //for (auto it = spectra2D.begin(); it != spectra2D.end(); ++it)
    // {
    //     auto htemp = it->spectrum->ToTH2(it->spectrum->POT(), kPOT);
    //     htemp->SetName(Form("Interaction_%s_Variable_%s_Cut_%s", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str()));
    //     htemp->Write();
    //     it->spectrum->SaveTo(f,Form("%s_%s_%s", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str()));
    //     printf("  %-18s %-16s %-22s %.6e \n", it->intName.c_str(), it->varName.c_str(), it->cutName.c_str(), it->spectrum->POT());
    // }

     f->Close();
    // std::cout << "\n"
    //           << std::endl;
    //Spectrum2D->toTH2( Spectrum2D->POT(), kPOT )
}
