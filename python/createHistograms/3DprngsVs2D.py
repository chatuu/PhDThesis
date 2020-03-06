#! /usr/bin/env python
import ROOT

f = ROOT.TFile( "../../cpp/rootFiles/2DVs3DPongs.root","READ" )
c1 = ROOT.TCanvas("c1", "Signal", 1024, 768)
ROOT.gPad.SetRightMargin(0.15)
c2 = ROOT.TCanvas("c2", "Background", 1024, 768)
ROOT.gPad.SetRightMargin(0.15)
hist1 = f.Get( "Interaction_RecoCohSig_Variable_Number of 3D Prongs and 2D Prongs_Cut_muonID" )
hist2 = f.Get( "Interaction_RecoCohBkgd_Variable_Number of 3D Prongs and 2D Prongs_Cut_muonID" )

xTitle = "Number of 3D Prongs"
yTitle = "Number of 2D Prongs"



hist1.SetXTitle( xTitle )
hist1.SetYTitle( yTitle )
hist1.GetXaxis().SetRangeUser(0, 10)
hist1.GetXaxis().SetNdivisions(110)
hist1.GetYaxis().SetRangeUser(0, 10)
hist1.GetYaxis().SetNdivisions(112)
hist1.SetTitle( "2D Prongs Vs 3D Prongs (Signal) " )

hist2.SetXTitle( xTitle )
hist2.SetYTitle( yTitle )
hist2.GetXaxis().SetRangeUser(0, 10)
hist2.GetXaxis().SetNdivisions(110)
hist2.GetYaxis().SetRangeUser(0, 10)
hist2.GetYaxis().SetNdivisions(110)
hist2.SetTitle( "2D Prongs Vs 3D Prongs (Background) " )

c1.cd()
hist1.Draw( "COLZ" )
c1.SaveAs( "Signal2DVs3D.png" )

c2.cd()
hist2.Draw( "COLZ" )
c2.SaveAs( "Background2DVs3D.png" )

f.Close()

