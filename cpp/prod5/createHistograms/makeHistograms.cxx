#include "headers.h"
#include "structs.h"
#include "makeHistograms.h"

makeHistograms::makeHistograms(/* args */)
{
}

makeHistograms::~makeHistograms()
{
    delete f1;
    delete key;
    delete subkey;

    hist.clear();
    selectedHists.clear();
    eventTable.clear();
}

void makeHistograms::SetRootFile(std::string fileName = "../rootFiles/Cleaned_20191028_backup.root")
{
    f1 = TFile::Open(fileName.c_str(), "READ");
    if (!f1 || f1->IsZombie())
    {
        exit(-1);
    }
    else
    {
        printf("ROOT file %s successfully opened..! \n", fileName.c_str());
    }
}

void makeHistograms::PrintHistograms()
{
    histNum = 0;
    TIter next(f1->GetListOfKeys());
    while ((key = (TKey *)next()))
    {
        dir = (TDirectory *)f1->Get(key->GetName());
        TIter subnext(dir->GetListOfKeys());
        while ((subkey = (TKey *)subnext()))
            if (std::string(subkey->GetClassName()) == "TH1D" && std::string(subkey->GetName()) == "hist")
            {
                printf("%d %s \n", histNum, std::string(key->GetName()).c_str());
                histNum++;
            }
    }
}

std::vector<histogram> makeHistograms::GetHists(std::vector<unsigned int> histNums)
{
    histNum = 0;
    TIter next(f1->GetListOfKeys());
    while ((key = (TKey *)next()))
    {
        dir = (TDirectory *)f1->Get(key->GetName());
        TIter subnext(dir->GetListOfKeys());
        while ((subkey = (TKey *)subnext()))
            if (std::string(subkey->GetClassName()) == "TH1D" && std::string(subkey->GetName()) == "hist")
            {
                TH1D *htemp = (TH1D *)dir->Get(subkey->GetName());
                htemp->SetDirectory(0);
                hist.push_back({histNum, std::string(key->GetName()), htemp});
                histNum++;
            }
    }
    histNum = 0;
    selectedHists.clear();
    for (auto it = selectedHists.begin(); it != selectedHists.end(); ++it)
        std::cout << it->number << it->name << std::endl;
    for (auto it = histNums.begin(); it != histNums.end(); ++it)
    {
        std::cout << *it << std::endl;
        selectedHists.push_back({histNum, hist.at(*it).name, hist.at(*it).hist});
        histNum++;
    }
    return selectedHists;
}

void makeHistograms::PrintNevents()
{
    for (auto it = selectedHists.begin(); it != selectedHists.end(); ++it)
        std::cout << " " << it->number << " " << it->name << " " << it->hist->GetEntries() << std::endl;
}

void makeHistograms::createCutTableFile(std::vector<unsigned int> histNums, std::string fileName = "defaultCutTable.txt")
{
    std::vector<histogram> selectedHistos = GetHists(histNums);
    std::ofstream output;
    output.open(fileName);
    output << "Histogram Number, Histogram Name, Number of Events\n";
    for (auto it = selectedHistos.begin(); it != selectedHistos.end(); ++it)
        output << it->number + 1 << "," << it->name << "," << it->hist->GetEntries() << "\n";

    output.close();
    std::cout << "File that contains cutTable : " << fileName << " is created..!" << std::endl;
}

void makeHistograms::drawHistograms(histogramInfo histInfo)
{
    int numPlots = histInfo.hist.size();
    auto c = new TCanvas("c", "title", 1024, 768);
    auto l = new TLegend(0.3, 0.7 - 0.04 * numPlots, 0.8, 0.9);
    THStack hs("hs", histInfo.title.c_str());
    TH1D *h[histInfo.hist.size()];
    for (auto it = histInfo.hist.begin(); it != histInfo.hist.end(); ++it)
    {
        h[it->number] = (TH1D *)it->hist;
        //h[it->number]->Scale((8.09e+20) / (3.53871e+21)); //Need to make a function for scaling Factor
        h[it->number]->SetLineColor(histInfo.lColorAndLeg.at(it->number).lineColor);
        h[it->number]->GetXaxis()->SetTitle(histInfo.axis.c_str());
        h[it->number]->GetYaxis()->SetNdivisions(1010);
        h[it->number]->GetXaxis()->SetNdivisions(1010);
        h[it->number]->GetXaxis()->SetRangeUser(histInfo.XRange[0], histInfo.XRange[1]);
        h[it->number]->Scale(histInfo.Scale);
        //h[it->number]->GetYaxis()->SetRangeUser(0, 50);

        eventTable.push_back({histInfo.title, h[it->number]->Integral(1, -1)});

        l->AddEntry(h[it->number], histInfo.lColorAndLeg.at(it->number).legName.c_str(), "l");
        h[it->number]->Draw("hist same");
    }
    l->Draw();
    c->SaveAs(histInfo.saveAs.c_str());
    f1->Close();
    delete c;
    delete l;
}

void makeHistograms::DivideHists(std::vector<unsigned int> histNums,
                                 std::vector<histDetails> histinfo,
                                 std::vector<int> XAxisRange,
                                 std::string fileName)
{

    TCanvas *c1 = new TCanvas("c1", "True #pi^{+} K.E curves", 1024, 768);
    auto l = new TLegend(0.6, 0.7, 0.9, 0.4);
    c1->cd();

    TH1D *h1 = nullptr;
    TH1D *h2 = nullptr;
    TH1D *h3 = nullptr;

    std::vector<histogram> hists = GetHists(histNums);

    //	for(auto it=histNums.begin();it!=histNums.end();++it)
    //		std::cout<<*it<<"\n";

    h1 = (TH1D *)hists.at(0).hist->Clone("h1");
    h2 = (TH1D *)hists.at(1).hist->Clone("h2");
    h3 = (TH1D *)h1->Clone("h3");
    h3->Divide(h2);
    //
    h3->GetYaxis()->SetTitle("Efficiency");
    h3->GetXaxis()->SetTitle("True #pi^{+} K.E (GeV)");
    h3->GetXaxis()->SetRangeUser(XAxisRange[0], XAxisRange[1]);
    h3->SetLineColor(histinfo.at(0).lineColor);
    //            //hEff[it->number]->GetYaxis()->SetRangeUser(0, 500);
    //
    h3->Draw("hist same");
    l->AddEntry(h3, histinfo.at(0).legName.c_str(), "l");
    //l->Draw();

    c1->SaveAs(fileName.c_str());
    delete c1;
    delete l;
}

void makeHistograms::makeCumulativeECurve(std::vector<histogram> signalPlots, std::vector<histDetails> histinfo, std::vector<float> XAxisRange, std::string XLabel, std::string fileName)
{
    //std::vector<histogram> signalPlots = GetHists(signalHists);
    std::vector<histogramInfo> efficiencyCurves;
    TH1D *hTrue = nullptr;
    //TH1D *hTrue2 = nullptr;
    TH1D *hCut = nullptr;

    TCanvas *c1 = new TCanvas("c1", "Efficiency curves", 1024, 768);
    int numPlots = (signalPlots.size() - 1);
    TH1D *hEff[numPlots];
    auto l = new TLegend(0.3, 0.5, 0.5, 0.7);
    c1->cd();
    unsigned int count = 0;
    for (auto it = signalPlots.begin(); it != signalPlots.end(); ++it)
    {
        std::cout << "Signal Hists: " << it->number << it->name << std::endl;
        if (count == 0)
        {

            hTrue = (TH1D *)it->hist;
            hTrue->SetDirectory(0);
            //hTrue2 = (TH1D *)it->hist;
            //hTrue2->SetDirectory(0);
            std::cout << "Truth hist was selected as: " << it->name << std::endl;
            count++;
            //hTrue1->Draw("hist same");

            //std::cout<<"N Bins: "<<nBins<<"\n";

            //delete htemp;
        }

        else
        {
            hCut = (TH1D *)it->hist;
            hCut->SetDirectory(0);
            hEff[count] = (TH1D *)it->hist->Clone(TString::Format("h%d", count));
            hEff[count]->SetDirectory(0);
            hEff[count]->Divide(hTrue);
            std::cout << "Testing: " << hEff[count]->Integral() << "\n";
            hEff[count]->GetYaxis()->SetTitle("Efficiency");
            hEff[count]->GetXaxis()->SetTitle(XLabel.c_str());
            hEff[count]->GetXaxis()->SetRangeUser(XAxisRange[0], XAxisRange[1]);
            hEff[count]->GetYaxis()->SetRangeUser(0, 1);
            hEff[count]->SetLineColor(histinfo.at(count).lineColor);
            std::cout << "line color is: " << histinfo.at(count).lineColor << std::endl;
            //hEff[count]->GetYaxis()->SetRangeUser(0, 500);

            // if (count == 1)
            // {
            //     TH1D *h1 = (TH1D *)it->hist->Clone("h1");
            //     //h1->Divide(hTrue1);
            //     //h1->GetYaxis()->SetRangeUser(0, 1);
            //     h1->Draw("hist same");
            // }
            // if (count == 2)
            // {
            //     TH1D *h2 = (TH1D *)it->hist->Clone("h2");
            //     //h2->Divide(hTrue2);
            //     h2->GetYaxis()->SetRangeUser(0, 1);
            //     h2->Draw("hist same");

            //     unsigned int nBins = h2->GetXaxis()->GetNbins();

            //     for (unsigned int i = 1; i <= nBins; ++i)
            //     {
            //         std::cout << h2->GetBinContent(i) << "\n";
            //     }

            hEff[count]->Draw("hist same");
            std::cout << count << '\n';
            l->AddEntry(hEff[count], histinfo.at(count).legName.c_str(), "l");
            l->Draw();
            count++;
        }
    }

    c1->SaveAs(fileName.c_str());
    delete c1;
    delete l;
    //delete h1;
    //delete h2;
}

void makeHistograms::makeEfficiencyCurve(std::vector<unsigned int> signalHists, std::vector<histDetails> histinfo, std::vector<int> XAxisRange, std::string XLabel, std::string fileName)
{
    std::vector<histogram> signalPlots = GetHists(signalHists);
    std::vector<histogramInfo> efficiencyCurves;
    TH1D *hTrue = nullptr;
    TH1D *hCut = nullptr;

    TCanvas *c1 = new TCanvas("c1", "Efficiency curves", 1024, 768);
    int numPlots = (signalPlots.size() - 1);
    TH1D *hEff[numPlots];
    auto l = new TLegend(0.6, 0.8, 0.9, 0.5);
    c1->cd();
    for (auto it = signalPlots.begin(); it != signalPlots.end(); ++it)
    {
        std::cout << "Signal Hists: " << it->number << it->name << std::endl;
        if (it->number == 0)
        {
            hTrue = (TH1D *)it->hist;
            std::cout << "Truth hist was selected" << std::endl;
        }
        else
        {
            hCut = (TH1D *)it->hist;
            hEff[it->number] = (TH1D *)hCut->Clone(TString::Format("h%d", it->number));
            hEff[it->number]->Divide(hTrue);
            std::cout << "Testing: " << hEff[it->number]->Integral() << "\n";
            hEff[it->number]->GetYaxis()->SetTitle("Efficiency");
            hEff[it->number]->GetXaxis()->SetTitle(XLabel.c_str());
            hEff[it->number]->GetXaxis()->SetRangeUser(XAxisRange[0], XAxisRange[1]);
            hEff[it->number]->GetYaxis()->SetRangeUser(0, 1);
            hEff[it->number]->SetLineColor(histinfo.at(it->number).lineColor);
            std::cout << "line color is: " << histinfo.at(it->number).lineColor << std::endl;
            //hEff[it->number]->GetYaxis()->SetRangeUser(0, 500);

            hEff[it->number]->Draw("hist same");
            std::cout << it->number << '\n';
            l->AddEntry(hEff[it->number], histinfo.at(it->number).legName.c_str(), "l");
            l->Draw();
        }
    }
    c1->SaveAs(fileName.c_str());
    delete c1;
    delete l;
}

std::vector<histogram> makeHistograms::makeCumulativeCurve(std::vector<histogram> hist)
{
    std::vector<histogram> hists;
    unsigned int nBins = 0;
    double total = 0;
    double binVal = 0;
    unsigned int count = 0;
    TH1D *h[hist.size()];
    //std::cout<<"N hists "<<hist.size()<<"\n";

    //hist->Scale(0.71801);

    //TStyle *st1 = new TStyle("st1","my style");
    //st1->SetOptStat(0);
    //st1->cd();
    //TCanvas *c = new TCanvas("c", "Cumulative Distributions", 1024, 768);

    for (auto it = hist.begin(); it != hist.end(); ++it)
    {
        //std::cout<<it->number<<" "<<it->name<<" "<<"\n";
        count++;
        //std::cout<<"count is: "<<count<<"\n";
        it->hist->Scale(0.71801);
        nBins = it->hist->GetXaxis()->GetNbins();
        //std::cout<<"N Bins: "<<nBins<<"\n";
        h[count] = new TH1D(TString::Format("hc%d", count), "", nBins, 0, 1);
        h[count]->SetDirectory(0);
        for (unsigned int i = 1; i <= nBins; ++i)
        {
            //std::cout<<"Bin "<<i<<" bin content "<<it->hist->GetBinContent(i)<<"\n";
            //total = it->hist->Integral(i, nBins);
            total = it->hist->Integral(1, i); // For Reco |t|
            h[count]->SetBinContent(i, total);
        }
        hists.push_back({it->number, it->name, h[count]});
        //delete htemp;
    }

    //unsigned int nBins = hist->GetXaxis()->GetNbins();
    //double total = 0;
    //double binVal = 0;
    ////hist->Scale(0.71801);
    //TH1D* histo = new TH1D("histo", "Cumulative Histogram", nBins,0,1);
    //std::cout<< "Number of Bins: " << nBins <<"\n";

    //for(unsigned int i = 1; i <= nBins; ++i){
    //	total = hist->Integral(i,nBins);
    //	binVal = hist->GetBinContent(i);
    //	std::cout<<i<< " Bin Val is: "<<binVal <<" total is: "<<total <<"\n";
    //	histo->SetBinContent(i,total);
    //}
    //c->cd();
    //histo->Draw();
    //c->SaveAs("cumulativeTest.png");
    return hists;
}

void makeHistograms::makeCumulativePCurve(std::vector<histogram> signalPlots, std::vector<histogram> backgroundPlots, std::vector<histDetails> histinfo, std::vector<float> XAxisRange, std::string XLabel, std::string fileName)
{
    //std::vector<histogram> signalPlots = GetHists(signalHists);
    //std::vector<histogram> backgroundPlots = GetHists(backgroundHists);
    std::vector<purityBlock> purityHists;

    int numPlots = (signalPlots.size());
    unsigned int count = 0;
    std::cout << "Num Plots: " << numPlots << "\n";

    for (auto itSig = signalPlots.begin(); itSig != signalPlots.end(); ++itSig)
    {
        purityHists.push_back({signalPlots.at(count), backgroundPlots.at(count)});
        count++;
    }

    count = 0;
    std::reverse(purityHists.begin(), purityHists.end());
    std::reverse(histinfo.begin(), histinfo.end());

    TCanvas *c1 = new TCanvas("c1", "Purity curves", 1200, 768);
    int numSigPlots = (signalPlots.size());
    int numBackPlots = (backgroundPlots.size());
    if (numSigPlots == numBackPlots)
    {
        TH1D *hSignal[numSigPlots];
        TH1D *hBack[numSigPlots];
        TH1D *hNorm[numSigPlots];
        TH1D *hDenom[numSigPlots];

        auto l = new TLegend(0.3, 0.5, 0.5, 0.7);
        c1->cd();
        for (auto it = purityHists.begin(); it != purityHists.end(); ++it)
        {
            hSignal[count] = (TH1D *)it->signal.hist;
            hBack[count] = (TH1D *)it->background.hist;
            hDenom[count] = (TH1D *)hSignal[count]->Clone(TString::Format("hs%d", count));
            TH1D *htemp = (TH1D *)hBack[count]->Clone(TString::Format("hb%d", count));
            hDenom[count]->Add(htemp);
            hNorm[count] = (TH1D *)it->signal.hist;
            TH1D *htemp2 = (TH1D *)hDenom[count]->Clone(TString::Format("hdenom%d", count));
            hNorm[count]->Divide(htemp2);

            //if (it->signal.number != 0 && it->background.number != 0)
            //{
            hNorm[count]->GetYaxis()->SetTitle("Purity");
            hNorm[count]->GetXaxis()->SetRangeUser(XAxisRange[0], XAxisRange[1]);
            hNorm[count]->GetXaxis()->SetTitle(XLabel.c_str());

            hNorm[count]->Draw("hist same");
            hNorm[count]->SetLineColor(histinfo.at(count).lineColor);
            l->AddEntry(hNorm[count], histinfo.at(count).legName.c_str(), "l");
            std::cout << "Now Drawing: " << it->signal.name << "And it's legend: " << histinfo.at(count).legName.c_str() << '\n';
            count++;
            //}
        }
        l->Draw();
        c1->SaveAs(fileName.c_str());
    }
}

void makeHistograms::makePurityCurve(std::vector<unsigned int> signalHists, std::vector<unsigned int> backgroundHists, std::vector<histDetails> histinfo, std::vector<int> XAxisRange, std::string XLabel, std::string fileName)
{
    std::vector<histogram> signalPlots = GetHists(signalHists);
    std::vector<histogram> backgroundPlots = GetHists(backgroundHists);
    std::vector<purityBlock> purityHists;

    int numPlots = (signalPlots.size() - 1);

    for (auto itSig = signalPlots.begin(); itSig != signalPlots.end(); ++itSig)
        purityHists.push_back({signalPlots.at(itSig->number), backgroundPlots.at(itSig->number)});

    std::reverse(purityHists.begin(), purityHists.end());

    TCanvas *c1 = new TCanvas("c1", "Efficiency curves", 1200, 768);
    int numSigPlots = (signalPlots.size() - 1);
    int numBackPlots = (backgroundPlots.size() - 1);
    if (numSigPlots == numBackPlots)
    {
        TH1D *hSignal[numSigPlots];
        TH1D *hBack[numSigPlots];
        TH1D *hNorm[numSigPlots];
        TH1D *hDenom[numSigPlots];

        auto l = new TLegend(0.7, 0.5, 0.9, 0.9);
        c1->cd();
        for (auto it = purityHists.begin(); it != purityHists.end(); ++it)
        {
            hSignal[it->signal.number] = (TH1D *)it->signal.hist;
            hBack[it->background.number] = (TH1D *)it->background.hist;
            hDenom[it->signal.number] = (TH1D *)hSignal[it->signal.number]->Clone(TString::Format("h%d", it->signal.number));
            TH1D *htemp = (TH1D *)hBack[it->background.number]->Clone(TString::Format("hb%d", it->background.number));
            hDenom[it->signal.number]->Add(htemp);
            hNorm[it->signal.number] = (TH1D *)it->signal.hist;
            TH1D *htemp2 = (TH1D *)hDenom[it->signal.number]->Clone(TString::Format("hdenom%d", it->signal.number));
            hNorm[it->signal.number]->Divide(htemp2);

            if (it->signal.number != 0 && it->background.number != 0)
            {
                hNorm[it->signal.number]->GetYaxis()->SetTitle("Purity");
                hNorm[it->signal.number]->GetXaxis()->SetRangeUser(XAxisRange[0], XAxisRange[1]);
                hNorm[it->signal.number]->GetXaxis()->SetTitle(XLabel.c_str());

                hNorm[it->signal.number]->Draw("hist same");
                hNorm[it->signal.number]->SetLineColor(histinfo.at(it->signal.number - 1).lineColor);
                l->AddEntry(hNorm[it->signal.number], histinfo.at((it->signal.number - 1)).legName.c_str(), "l");
                std::cout << "Now Drawing: " << it->signal.name << "And it's legend: " << histinfo.at(it->signal.number - 1).legName.c_str() << '\n';
            }
        }
        l->Draw();
        c1->SaveAs(fileName.c_str());
    }
}

void makeHistograms::SetHistInfo(std::vector<histogram> hists,
                                 std::vector<histDetails> lColorAndLeg,
                                 std::string title,
                                 std::string axis,
                                 std::vector<int> XAxisRange,
                                 double scale,
                                 std::string saveAs)
{
    histInfo.hist = hists;
    histInfo.lColorAndLeg = lColorAndLeg;
    histInfo.title = title;
    histInfo.axis = axis;
    histInfo.XRange = XAxisRange;
    histInfo.Scale = scale;
    histInfo.saveAs = saveAs;

    drawHistograms(histInfo);
}
