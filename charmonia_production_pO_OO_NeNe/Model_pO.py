import sys
import argparse
import yaml
import ROOT
import numpy as np
import array
import pandas as pd

sys.path.append('../utils')
from plot_library import LoadStyle, SetHistStat, SetHistSyst, SetGraStat, SetGraSyst, SetLegend
from fast_analysis_utils import ExtractFromYaml


def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--run", help="Compute cross section", action="store_true")
    args = parser.parse_args()
    

    with open(args.cfgFileName, 'r') as yml_cfg:
        inputCfg = yaml.load(yml_cfg, yaml.FullLoader)

    if args.run:
        compareToModels(inputCfg)

def compareToModels(config):
    """
    function to compare results with theory models
    """
    LoadStyle()
    ROOT.gStyle.SetOptStat(False)

    dfJpsiRaaThuVsPt = pd.read_csv(config["inputs"]["fInJpsiRaaModelVsPt"],sep=r"\s+",engine="python")
    print(dfJpsiRaaThuVsPt.head())
    print(dfJpsiRaaThuVsPt.columns)
    print(dfJpsiRaaThuVsPt.dtypes)

    jpsiRaaThuVsPt = pd.to_numeric(dfJpsiRaaThuVsPt["R_ave"])
    jpsiRaaThuVsPt_Up = pd.to_numeric(dfJpsiRaaThuVsPt["R_min"])
    jpsiRaaThuVsPt_Low = pd.to_numeric(dfJpsiRaaThuVsPt["R_max"])
    ptCenters = pd.to_numeric(dfJpsiRaaThuVsPt["pT/GeV"])
    ptWidths = np.full_like(ptCenters, 1e-10, dtype=float)
    errUp = jpsiRaaThuVsPt_Up - jpsiRaaThuVsPt
    errLow = jpsiRaaThuVsPt - jpsiRaaThuVsPt_Low

    nPoints = len(ptCenters)

    ptWidthsArray = np.full(nPoints, 1e-8)  # same as your tiny width
    errUpLineArray = np.full(nPoints, 1e-6)
    errDownLineArray = np.full(nPoints, 1e-6)

    graJpsiRaaThuVsPt = ROOT.TGraphAsymmErrors(len(ptCenters),np.array(ptCenters),np.array(jpsiRaaThuVsPt),np.array(ptWidths), np.array(ptWidths),np.array(errLow), np.array(errUp))
    graJpsiRaaThuVsPt_line = ROOT.TGraphAsymmErrors(nPoints,np.array(ptCenters, dtype=float),np.array(jpsiRaaThuVsPt, dtype=float),ptWidthsArray,ptWidthsArray,errDownLineArray,errUpLineArray)

    graJpsiRaaThuVsPt.SetFillColorAlpha(ROOT.kRed-7, 0.1)
    graJpsiRaaThuVsPt.SetFillStyle(1001)
    graJpsiRaaThuVsPt.SetLineColor(ROOT.kBlack+4)
    graJpsiRaaThuVsPt.SetLineWidth(2) 
    graJpsiRaaThuVsPt_line.SetLineColor(ROOT.kRed+1)
    graJpsiRaaThuVsPt_line.SetLineWidth(2)
    graJpsiRaaThuVsPt_line.SetLineWidth(3)    

    histGridJpsiRaaVsPt = ROOT.TH2D("histGridJpsiRaaVsPt", ";#it{p}_{T} (GeV/#it{c});#it{R}_{pO}", 100, 0, 10, 100, 0, 1.5)
    canvasJpsiRaaThuVsPt = ROOT.TCanvas("canvasJpsiRaapOVsPt", "", 800, 600)
    histGridJpsiRaaVsPt.Draw()
    graJpsiRaaThuVsPt.Draw("E3 SAME")
    graJpsiRaaThuVsPt_line.Draw("L SAME")

    lineUnityVsPt = ROOT.TLine(0, 1, 10, 1)
    lineUnityVsPt.SetLineColor(ROOT.kGray+1)
    lineUnityVsPt.SetLineWidth(2)
    lineUnityVsPt.SetLineStyle(ROOT.kDashed)
    lineUnityVsPt.Draw()

    latexTitle = ROOT.TLatex()
    latexTitle.SetNDC()
    latexTitle.SetTextSize(0.05)
    latexTitle.SetTextFont(42)
    latexTitle.DrawLatex(0.20, 0.86, "pO, #sqrt{#it{s}} = 9.62 TeV, 2.15<#it{y}<3.65")

    legendRaaVsPt = ROOT.TLegend(0.66, 0.19, 0.98, 0.44, " ", "brNDC")
    SetLegend(legendRaaVsPt)
    legendRaaVsPt.SetTextSize(0.045)
    legendRaaVsPt.AddEntry(graJpsiRaaThuVsPt_line, "FCEL, color octet (F. Arleo et al.)", "L")
    legendRaaVsPt.SetBorderSize(0) 
    legendRaaVsPt.SetFillStyle(0) 
    legendRaaVsPt.Draw("SAME")

    xMinBox = 7.8
    xMaxBox = 8.0
    yRef = 1.0
    systCommonTotData = 0.050744457825461095 

    yMinBox = yRef * (1.0 - systCommonTotData)
    yMaxBox = yRef * (1.0 + systCommonTotData)

    boxCommon_Data = ROOT.TBox(xMinBox, yMinBox, xMaxBox, yMaxBox)
    boxCommon_Data.SetFillColorAlpha(ROOT.kRed+1, 0.4)
    boxCommon_Data.SetLineColor(ROOT.kRed+1)
    boxCommon_Data.SetLineWidth(1)


    canvasJpsiRaaThuVsPt.Update()
    
    input()


if __name__ == '__main__':
    main()