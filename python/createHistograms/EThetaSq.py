#! /usr/bin/env python
import ROOT

f = ROOT.TFile("../../cpp/rootFiles/NewFile.root", "READ")
c1 = ROOT.TCanvas("c1", "Signal", 1024, 768)
c2 = ROOT.TCanvas("c2", "Background", 1024, 768)
c3 = ROOT.TCanvas("c3", "SignalVsBackground", 1024, 768)

sigDir = f.Get("Interaction_RecoCohSig_Variable_EThetaSquared_Cut_TwoProng")
sigDir.cd()

hist1 = sigDir.Get("hist")

backDir = f.Get("Interaction_RecoCohBkgd_Variable_EThetaSquared_Cut_TwoProng")
backDir.cd()


hist2 = backDir.Get("hist")

xTitle = "E_{#pi}(#theta_{#pi})^{2}(GeV)"
yTitle = "Number of Events"

hist1.SetXTitle(xTitle)
hist1.SetYTitle(yTitle)
hist1.GetXaxis().SetRangeUser(0, 1)
# hist1.GetXaxis().SetNdivisions(110)
#hist1.GetYaxis().SetRangeUser(1, 10)
# hist1.GetYaxis().SetNdivisions(112)
hist1.SetTitle("|t| (Signal) ")
hist1.SetLineColor(2)

hist2.SetXTitle(xTitle)
hist2.SetYTitle(yTitle)
hist2.GetXaxis().SetRangeUser(0, 1)
# hist2.GetXaxis().SetNdivisions(110)
#hist2.GetYaxis().SetRangeUser(1, 10)
# hist2.GetYaxis().SetNdivisions(110)
hist2.SetTitle("|t| (Background) ")
hist2.SetLineColor(4)

c1.cd()
hist1.Draw("hist")
c1.SaveAs("SignalEThetaSq.png")

c2.cd()
hist2.Draw("hist")
c2.SaveAs("BackgroundEThetaSq.png")

c3.cd()
hist2.Draw("hist")
hist1.Draw("hist same")
legend = ROOT.TLegend(0.4, 0.7, 0.78, 0.9)
legend.AddEntry(hist1, "Signal", "l")
legend.AddEntry(hist2, "Background", "l")
legend.Draw()
c3.SaveAs("SignalVsBackgroundEThetaSq.png")

f.Close()
