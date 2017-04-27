#! /usr/bin/env python

import ROOT
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetNumberContours(255)
ROOT.gStyle.SetPalette(ROOT.kPastel)

prefs = ["mumu","emu"]
samples = ["all","dyjetsll","top","ttdl","ttsl","ttV","singletop","wjets","ww_2l2nu"]

for pref in prefs:
    for samp in samples:

        fin = ROOT.TFile("../looper/output/test/{0}.root".format(samp))
        
        h = fin.Get("{0}_zmass/h_DMAvgM_binned".format(pref))

        c = ROOT.TCanvas()
        c.SetCanvasSize(700,504)
        c.SetLogz(1)

        h.Draw("COLZ TEXT")

        c.SaveAs("plots/test/other/yields_{0}_{1}.png".format(pref,samp))
        c.SaveAs("plots/test/other/yields_{0}_{1}.pdf".format(pref,samp))

        fin.Close()



