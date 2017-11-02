#! /usr/bin/env python

import sys
sys.path.append("../MT2Analysis/scripts/CRplotMaker")
import ROOT
ROOT.gROOT.SetBatch(1)

from MT2PlotMaker import *
import MT2PlotDefs as pd

plotdefs = [
    ("M1",True,None,None,2),
    ("M2",True,None,None,2),
    ("MLL",True,None,None,4),
    ("DM",True,None,None,2),
    ("AvgM",True,None,None,2),
    ("met",True,None,None,2),
    ("Mbb",True,None,None,4),
    ("deltaPhiLL",True,None,None,4),
    ("deltaPhibb",True,None,None,4),
    ("deltaPhiMETJet",True,None,None,4),
    ("bidx2",False,None,None,1),
    ("bidx1",False,None,None,1),
    ("ST",True,None,None,4),
]

input_dir = "../looper/output/test2"
output_dir = "plots/test"

pd.lumi = 36.8
pd.lumiUnit = "fb"

exts = ["pdf","png"]

signals = ["T2tt_RPV_700","T2tt_RPV_900","T2tt_RPV_1100"]

MT2PlotMaker(input_dir, ["top","dyjetsll","wjets","ww_2l2nu"], None, "mumu_baseline", plotdefs, output_dir, exts, signals=signals, opts="noSubtitles")
MT2PlotMaker(input_dir, ["top","dyjetsll","wjets","ww_2l2nu"], None, "emu_baseline", plotdefs, output_dir, exts, signals=signals, opts="noSubtitles")
MT2PlotMaker(input_dir, ["top","dyjetsll","wjets","ww_2l2nu"], None, "mumu_zmass", plotdefs, output_dir, exts, signals=signals, opts="noSubtitles")
MT2PlotMaker(input_dir, ["top","dyjetsll","wjets","ww_2l2nu"], None, "emu_zmass", plotdefs, output_dir, exts, signals=signals, opts="noSubtitles")
MT2PlotMaker(input_dir, ["top","dyjetsll","wjets","ww_2l2nu"], None, "mumu_sr", plotdefs, output_dir, exts, signals=signals, opts="noSubtitles")
MT2PlotMaker(input_dir, ["top","dyjetsll","wjets","ww_2l2nu"], None, "emu_sr", plotdefs, output_dir, exts, signals=signals, opts="noSubtitles")
MT2PlotMaker(input_dir, ["top","dyjetsll","wjets","ww_2l2nu"], None, "mumu_sr_AvgM500", plotdefs, output_dir, exts, signals=signals, opts="noSubtitles")
MT2PlotMaker(input_dir, ["top","dyjetsll","wjets","ww_2l2nu"], None, "emu_sr_AvgM500", plotdefs, output_dir, exts, signals=signals, opts="noSubtitles")
MT2PlotMaker(input_dir, ["top","dyjetsll","wjets","ww_2l2nu"], None, "mumu_sr_AvgM500_Btop2", plotdefs, output_dir, exts, signals=signals, opts="noSubtitles")
MT2PlotMaker(input_dir, ["top","dyjetsll","wjets","ww_2l2nu"], None, "emu_sr_AvgM500_Btop2", plotdefs, output_dir, exts, signals=signals, opts="noSubtitles")
MT2PlotMaker(input_dir, ["top","dyjetsll","wjets","ww_2l2nu"], None, "mumu_crdm", plotdefs, output_dir, exts, signals=signals, opts="noSubtitles")
MT2PlotMaker(input_dir, ["top","dyjetsll","wjets","ww_2l2nu"], None, "emu_crdm", plotdefs, output_dir, exts, signals=signals, opts="noSubtitles")

## make 2D plots

# fin = ROOT.TFile(os.path.join(input_dir,"all.root"))

# for dir in ["mumu_baseline","emu_baseline", "mumu_zmass", "emu_zmass"]:
# # for dir in ["mumu_baseline","mumu_zmass"]:

#     c1 = ROOT.TCanvas()
#     c1.SetCanvasSize(700,504)
#     c1.SetLogz(True)
#     ROOT.gStyle.SetNumberContours(255)

#     h_M1M2 = fin.Get("{0}/h_M1M2".format(dir))
#     h_M1M2.Draw("COLZ")

#     c1.SaveAs(os.path.join(output_dir, dir, "{0}_M1M2.pdf".format(dir)))
#     c1.SaveAs(os.path.join(output_dir, dir, "{0}_M1M2.png".format(dir)))

#     c1.Clear()
    
#     h_DMAvgM = fin.Get("{0}/h_DMAvgM".format(dir))
#     h_DMAvgM.Draw("COLZ")
    
#     c1.SaveAs(os.path.join(output_dir, dir, "{0}_DMAvgM.pdf".format(dir)))
#     c1.SaveAs(os.path.join(output_dir, dir, "{0}_DMAvgM.png".format(dir)))
    
# fin.Close()


