#! /usr/bin/env python
import ROOT

f = ROOT.TFile( "../../cpp/rootFiles/EThetaSq.root","READ" )
c1 = ROOT.TCanvas("c1", "Signal", 1024, 768)
c2 = ROOT.TCanvas("c2", "Background", 1024, 768)
c3 = ROOT.TCanvas("c3", "Signal Vs Background", 1024, 768)
sigHist = f.Get('Interaction_RecoCohSig_Variable_thetaSq Vs TruePion Kinetic Energy_Cut_muonID')
backHist = f.Get('Interaction_RecoCohBkgd_Variable_thetaSq Vs TruePion Kinetic Energy_Cut_muonID')

xTitle = 'Energy of #pi (GeV)'
yTitle = '#theta^{2} (Rad^{2})'

sigHist.SetXTitle( xTitle )
sigHist.SetYTitle( yTitle )
sigHist.GetXaxis().SetRangeUser(0, 3)
sigHist.GetXaxis().SetNdivisions(110)
#sigHist.GetYaxis().SetRangeUser(0, 1200)
sigHist.GetYaxis().SetNdivisions(112)
sigHist.SetTitle("#theta^{2} Vs Energy of #pi (Signal)")
#sigHist.SetLineColor(2)

c1.cd()
sigHist.Draw('COLZ')
c1.SaveAs('EVsThetaSquare_Signal.png')

backHist.SetXTitle( xTitle )
backHist.SetYTitle( yTitle )
backHist.GetXaxis().SetRangeUser(0, 3)
backHist.GetXaxis().SetNdivisions(110)
#backHist.GetYaxis().SetRangeUser(0, 1200)
backHist.GetYaxis().SetNdivisions(112)
backHist.SetTitle("#theta^{2} Vs Energy of #pi (Background)")
backHist.SetLineColor(4)

c2.cd()
backHist.Draw('COLZ')

c2.SaveAs('EVsThetaSquare_Background.png')

# c3.cd()

# backHist.Draw('hist')
# sigHist.Draw('hist same')
# legend = ROOT.TLegend(0.1, 0.7, 0.48, 0.9)
# legend.AddEntry(sigHist,"Signal","l")
# legend.AddEntry(backHist, "Background", "l")
# legend.Draw()
# c3.SaveAs('HadronicActivity_SignalVsBackground.png')






f.Close()
