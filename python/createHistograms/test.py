#!/usr/bin/env python

import ROOT
c = ROOT.TCanvas("c")
h1 = ROOT.TH1I("h1", "hist", 10, 0, 10)
h1.Fill(0)
h1.Fill(0)
h1.Fill(1)
h1.Fill(2)
h1.Fill(2)
h1.Fill(2)
h1.Fill(3)
h1.Fill(3)
h1.Fill(3)
h1.Fill(3)
h1.Draw()
c.SaveAs('test.png')
print(h1.GetBinContent(2))
print(h1.GetBin(2))
print(h1.GetNbinsX())

h2 = ROOT.TH1I("h2", "hist", 10, 0, 10)


for bin in range(h1.GetNbinsX()):
    print("Bin Number: %d Bin Content: %d" %
          (bin + 1, h1.GetBinContent(bin + 1)))
    h2.SetBinContent(bin, h1.GetBinContent(bin))
c1 = ROOT.TCanvas("c1")
c1.cd()
h2.Draw()
c1.SaveAs('test1.png')

f = ROOT.TFile("../../cpp/rootFiles/TopologyCorrected.root", "READ")
hist1 = f.Get(
    "Interaction_RecoCohBkgd_Variable_Topology Vs Pion Kinetic Energy_Cut_muonID")
print(hist1.GetNbinsX())
print(hist1.GetNbinsY())
for ix in range(hist1.GetNbinsX()):
    print(ix+1)
