#! /usr/bin/env python

import ROOT

ROOT.gStyle.SetOptStat(0)

f = ROOT.TFile("../../cpp/rootFiles/TopologyCorrected.root", "READ")

c1 = ROOT.TCanvas("c1", "Background", 3024, 968)
hist1 = f.Get(
    "Interaction_RecoCohBkgd_Variable_Topology Vs Pion Kinetic Energy_Cut_muonID")



xTitle = "Second most Energetic Particle's K.E (GeV)"
#yTitle = "Topology"

hist1.SetXTitle(xTitle)
#hist1.SetYTitle( yTitle )
hist1.SetTitle(
    "Topology Vs Second most Energetic particle's Kinetic Energy (Background)")

yAxis = hist1.GetYaxis()
yLabelArray = []
yLabelArray.append("0 3D  0 2D")
yLabelArray.append("0 3D  1 2D")
yLabelArray.append("0 3D  2 2D")
yLabelArray.append("0 3D  3 2D")
yLabelArray.append("0 3D  4 2D")
yLabelArray.append("0 3D  >5 2D")
yLabelArray.append("1 3D  0 2D")
yLabelArray.append("1 3D  1 2D")
yLabelArray.append("1 3D  2 2D")
yLabelArray.append("1 3D  3 2D")
yLabelArray.append("1 3D  4 2D")
yLabelArray.append("1 3D  >5 2D")
yLabelArray.append("2 3D  0 2D")
yLabelArray.append("2 3D  1 2D")
yLabelArray.append("2 3D  2 2D")
yLabelArray.append("2 3D  3 2D")
yLabelArray.append("2 3D  4 2D")
yLabelArray.append("2 3D  >5 2D")
yLabelArray.append("3 3D  0 2D")
yLabelArray.append("3 3D  1 2D")
yLabelArray.append("3 3D  2 2D")
yLabelArray.append("3 3D  3 2D")
yLabelArray.append("3 3D  4 2D")
yLabelArray.append("3 3D  >5 2D")
yLabelArray.append("4 3D  0 2D")
yLabelArray.append("4 3D  1 2D")
yLabelArray.append("4 3D  2 2D")
yLabelArray.append("4 3D  3 2D")
yLabelArray.append("4 3D  4 2D")
yLabelArray.append("4 3D  >5 2D")
yLabelArray.append("5 3D  >0 2D")

i = 1
for label in yLabelArray:
    binIndex = yAxis.FindBin(i)
    yAxis.SetBinLabel(binIndex, label)
    i += 1

ROOT.gPad.SetRightMargin(0.15)

hist1.GetYaxis().SetTitleOffset(999)
hist1.GetYaxis().SetRangeUser(7, 32)
hist1.GetXaxis().SetRangeUser(0.15, 2)

c1.cd()
hist1.Draw("COLZ")
c1.SaveAs("TopologyBackgroundZoomIn.png")

print(hist1.GetNbinsX(), hist1.GetNbinsY())
print(hist1.GetXaxis().GetXmin(), hist1.GetXaxis().GetXmax())
print(hist1.GetYaxis().GetXmin(), hist1.GetYaxis().GetXmax())


hist2 = ROOT.TH2D("hist2", "TopologyBackgroundRatio", hist1.GetNbinsX(
), hist1.GetXaxis().GetXmin(), hist1.GetXaxis().GetXmax(), hist1.GetNbinsY(), hist1.GetYaxis().GetXmin(), hist1.GetYaxis().GetXmax())

for ix in range(hist1.GetNbinsX()):
    sum = 0
    for iy in range(hist1.GetNbinsY()):
        sum += hist1.GetBinContent(ix + 1, iy + 1)
    for iy in range(hist1.GetNbinsY()):
        if sum == 0:
            hist2.SetBinContent(ix + 1, iy + 1, 0)
        else:
            hist2.SetBinContent(
                ix + 1, iy + 1, (hist1.GetBinContent(ix + 1, iy + 1)) / sum)

c2 = ROOT.TCanvas("c2", "BackgroundRatio", 3024, 968)
c2.cd()

yAxis = hist2.GetYaxis()
i = 1
for label in yLabelArray:
    binIndex = yAxis.FindBin(i)
    yAxis.SetBinLabel(binIndex, label)
    i += 1

ROOT.gPad.SetRightMargin(0.15)

hist2.GetYaxis().SetTitleOffset(999)
hist2.GetYaxis().SetRangeUser(7, 32)
hist2.GetXaxis().SetRangeUser(0.15, 2)
xTitle = "Second most Energetic Particle's K.E (GeV)"
#yTitle = "Topology"

hist2.SetXTitle(xTitle)

hist2.Draw('COLZ')
c2.SaveAs('TopologyBackgroundRatio.png')


f.Close()
