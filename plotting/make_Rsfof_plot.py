#! /usr/bin/env python

import ROOT
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

fin = ROOT.TFile("../looper/output/test/top.root")

h_mumu = fin.Get("mumu_zmass/h_DMAvgM_binned")
h_emu = fin.Get("emu_zmass/h_DMAvgM_binned")

h_Rsfof = h_mumu.Clone("h_Rsfof")
h_Rsfof.Divide(h_emu)

c = ROOT.TCanvas()
c.SetCanvasSize(700,504)
ROOT.gStyle.SetNumberContours(255)

h_Rsfof.GetZaxis().SetRangeUser(0.25,0.75)
h_Rsfof.Draw("COLZ")

text = ROOT.TLatex()
text.SetTextAlign(22)
text.SetTextFont(62)
text.SetTextSize(0.022)
text.SetTextColor(ROOT.kBlack)
for i in range(h_Rsfof.GetNbinsX()):
    for j in range(h_Rsfof.GetNbinsY()):
        rat = h_Rsfof.GetBinContent(i+1,j+1)
        err  = h_Rsfof.GetBinError(i+1,j+1)

        x = h_Rsfof.GetXaxis().GetBinCenter(i+1)
        y = h_Rsfof.GetYaxis().GetBinCenter(j+1)
        text.DrawLatex(x, y, "{0:.2f}#pm{1:.2f}".format(rat,err))

text.SetTextFont(42)
text.SetTextSize(0.06)
text.DrawLatex(500, 2100, "R_{#mu#mu/e#mu} in top sample")        
        
c.SaveAs("plots/test/other/Rsfof.png")
c.SaveAs("plots/test/other/Rsfof.pdf")




