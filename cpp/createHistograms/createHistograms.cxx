#include "headers.h"
#include "structs.h"
#include "makeHistograms.h"

int main(int argc, char **argv)
{
    std::vector<unsigned int> SighistNums = {15, 6, 7, 8, 9};
    std::vector<histDetails> lSignal;
    lSignal.push_back({1, "Signal with No Cuts"});
    lSignal.push_back({4, "Signal with DQ Cuts"});
    lSignal.push_back({2, "Signal with Fiducial Cuts"});
    lSignal.push_back({8, "Signal with Containment Cuts"});
    lSignal.push_back({6, "Signal with muonID Cuts"});

    std::vector<unsigned int> BackhistNums = {16, 11, 12, 13, 14};
    std::vector<histDetails> lBackground;
    lBackground.push_back({1, "Background with No Cuts"});
    lBackground.push_back({4, "Background with DQ Cuts"});
    lBackground.push_back({2, "Background with Fiducial Cuts"});
    lBackground.push_back({8, "Background with Containment Cuts"});
    lBackground.push_back({6, "Background with muonID Cuts"});

    std::vector<unsigned int> SigVBackhistNums = {16, 11, 12, 13, 14, 15, 6, 7, 8, 9};
    std::vector<histDetails> lColorAndLeg;
    lColorAndLeg.push_back({1, "Background with No Cuts"});
    lColorAndLeg.push_back({4, "Background with DQ Cuts"});
    lColorAndLeg.push_back({2, "Background with Fiducial Cuts"});
    lColorAndLeg.push_back({8, "Background with Containment Cuts"});
    lColorAndLeg.push_back({6, "Background with muonID Cuts"});
    lColorAndLeg.push_back({9, "Signal with No Cuts"});
    lColorAndLeg.push_back({3, "Signal with DQ Cuts"});
    lColorAndLeg.push_back({5, "Signal with Fiducial Cuts"});
    lColorAndLeg.push_back({7, "Signal with Containment Cuts"});
    lColorAndLeg.push_back({28, "Signal with muonID Cuts"});

    makeHistograms nuESignalPlots;
    nuESignalPlots.SetRootFile(argv[1]);
    nuESignalPlots.PrintHistograms();
    nuESignalPlots.SetHistInfo(nuESignalPlots.GetHists(SighistNums),
                               lSignal,
                               "#nu_{#mu} Energy of the Coherent Signal (GeV)",
                               "#nu_{#mu} Energy (GeV)",
                               "./plots/signal.png");

    makeHistograms nuEBackgroundPlots;
    nuEBackgroundPlots.SetRootFile(argv[1]);
    nuEBackgroundPlots.PrintHistograms();
    nuEBackgroundPlots.SetHistInfo(nuEBackgroundPlots.GetHists(BackhistNums),
                                   lBackground,
                                   "#nu_{#mu} Energy of the Coherent Background (GeV)",
                                   "#nu_{#mu} Energy (GeV)",
                                   "./plots/Background.png");

    makeHistograms nuESigVsBack;
    nuESigVsBack.SetRootFile(argv[1]);
    nuESigVsBack.PrintHistograms();

    nuESigVsBack.SetHistInfo(nuESigVsBack.GetHists(SigVBackhistNums),
                             lColorAndLeg,
                             "#nu_{#mu} Energy Signal Vs Background (GeV)",
                             "#nu_{#mu} Energy (GeV)",
                             "./plots/SigVsBack.png");

    nuESigVsBack.PrintNevents();
    //nuESigVsBack.createCutTableFile("./csvFiles/cutTable.csv");

    makeHistograms createCutTable;
    createCutTable.SetRootFile(argv[1]);
    std::vector<unsigned int> histlist = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    std::vector<histDetails> lLeg;
    lLeg.push_back({1, "Background with No Cuts"});
    lLeg.push_back({4, "Background with DQ Cuts"});
    lLeg.push_back({2, "Background with Fiducial Cuts"});
    lLeg.push_back({8, "Background with Containment Cuts"});
    lLeg.push_back({6, "Background with muonID Cuts"});
    lLeg.push_back({9, "Signal with No Cuts"});
    lLeg.push_back({3, "Signal with DQ Cuts"});
    lLeg.push_back({5, "Signal with Fiducial Cuts"});
    lLeg.push_back({7, "Signal with Containment Cuts"});
    lLeg.push_back({28, "Signal with muonID Cuts"});
    lLeg.push_back({7, "Signal with Containment Cuts"});
    lLeg.push_back({28, "Signal with muonID Cuts"});

    createCutTable.SetHistInfo(createCutTable.GetHists(histlist),
                               lLeg,
                               "#nu_{#mu} Energy Signal Vs Background (GeV)",
                               "#nu_{#mu} Energy (GeV)",
                               "./plots/testing.png");

    createCutTable.createCutTableFile("./csvFiles/FinalcutTable.csv");
   
    std::vector<histDetails> sigdetails;
    sigdetails.push_back({1, "DQ Cuts"});
    sigdetails.push_back({2, "Fiducial Cuts"});
    sigdetails.push_back({3, "Containement Cuts"});
    sigdetails.push_back({4, "Muon ID Cuts"});
    createCutTable.makeEfficiencyCurve(SighistNums, sigdetails,"./plots/efficiency.png");

    std::vector<histDetails> PurityLeg;
    PurityLeg.push_back({1, "DQ Cuts"});
    PurityLeg.push_back({2, "Fiducial Cuts"});
    PurityLeg.push_back({3, "Containement Cuts"});
    PurityLeg.push_back({4, "Muon ID Cuts"});

    createCutTable.makePurityCurve(SighistNums, BackhistNums, PurityLeg,"./plots/purity.png");

    return 0;
}
