#! /usr/bin/env python
import ROOT

f = ROOT.TFile( "NewVertexActivity.root","READ" )
c1 = ROOT.TCanvas( "c1","Signal",1024,768 )
c2 = ROOT.TCanvas( "c2","Background",1024,768 )

sigDir = f.Get( "Interaction_RecoCohSig_Variable_vertexAct_Cut_muonID" )
sigDir.cd()

hist1 = sigDir.Get( "hist" )

backDir = f.Get( "Interaction_RecoCohBkgd_Variable_vertexAct_Cut_muonID" )
backDir.cd()


hist2 = backDir.Get( "hist" )

xTitle = "Vertex Activity (GeV)"
yTitle = "Number of Events"

hist1.SetXTitle( xTitle )
hist1.SetYTitle( yTitle )
hist1.GetXaxis().SetRangeUser(-1, 1)
#hist1.GetXaxis().SetNdivisions(110)
#hist1.GetYaxis().SetRangeUser(1, 10)
#hist1.GetYaxis().SetNdivisions(112)
hist1.SetTitle( "Vertex activity (Signal) " )

hist2.SetXTitle( xTitle )
hist2.SetYTitle( yTitle )
hist2.GetXaxis().SetRangeUser(-1, 1)
#hist2.GetXaxis().SetNdivisions(110)
#hist2.GetYaxis().SetRangeUser(1, 10)
#hist2.GetYaxis().SetNdivisions(110)
hist2.SetTitle( "Vertex Activity (Background) " )

c1.cd()
hist1.Draw( "hist" )
c1.SaveAs( "SignalVertex.png" )

c2.cd()
hist2.Draw( "hist" )
c2.SaveAs( "BackgroundVertex.png" )

f.Close()

