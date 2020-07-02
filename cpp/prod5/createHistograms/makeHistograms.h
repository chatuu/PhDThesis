#ifndef MAKEHISTOGRAMS_H
#define MAKEHISTOGRAMS_H

#include "headers.h"

class makeHistograms
{
private:
    TFile *f1 = nullptr;
    unsigned int histNum = 0;
    TKey *key = nullptr;
    TKey *subkey = nullptr;
    TDirectory *dir = nullptr;
    TDirectory *newDir = nullptr;
    histogramInfo histInfo;
    std::vector<histogram> hist;
    std::vector<histogram> selectedHists;
    //std::vector<double> nentries;
    std::vector<cutTable> eventTable;

public:
    makeHistograms(/* args */);
    ~makeHistograms();

    std::vector<histogram> GetHists(std::vector<unsigned int> histNums);
    void SetRootFile(std::string fileName);
    void PrintHistograms();
    std::vector<histogram> makeCumulativeCurve(std::vector<histogram> hist);
    void PrintNevents();
    void createCutTableFile(std::vector<unsigned int>, std::string fileName);
    void drawHistograms(histogramInfo histInfo);
    void makeEfficiencyCurve(std::vector<unsigned int> signalHists, std::vector<histDetails> histinfo, std::vector<int> XAxisRange, std::string XLabel, std::string fileName);
    void makePurityCurve(std::vector<unsigned int> signalHists, std::vector<unsigned int> backgroundHist, std::vector<histDetails> histinfo, std::vector<int> XAxisRange, std::string XLabel, std::string fileName);
    void makeCumulativePCurve(std::vector<histogram> signalPlots, std::vector<histogram> backgroundPlots, std::vector<histDetails> histinfo, std::vector<float> XAxisRange, std::string XLabel, std::string fileName);
    void makeCumulativeECurve(std::vector<histogram> signalHists, std::vector<histDetails> histinfo, std::vector<float> XAxisRange, std::string XLabel, std::string fileName);

    void SetHistInfo(std::vector<histogram> hists,
                     std::vector<histDetails> lColorAndLeg,
                     std::string title,
                     std::string axis,
                     std::vector<int> XAxisRange,
                     double scale,
                     std::string saveAs);

    void DivideHists(std::vector<unsigned int> histNums,
                     std::vector<histDetails> histinfo,
                     std::vector<int> XAxisRange,
                     std::string fileName);
};

#endif //MAKEHISTOGRAMS_H
