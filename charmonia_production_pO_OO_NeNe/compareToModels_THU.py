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


    dfJpsiRaaVsPt = pd.read_csv(config["inputs"]["fInJpsiRaaVsPt"], sep=' ')
    ptMin = dfJpsiRaaVsPt["x_min"]
    ptMax = dfJpsiRaaVsPt["x_max"]
    ptCenters = (ptMax + ptMin) / 2.
    ptWidths = (ptMax - ptMin) / 2.
    ptWidths_zero = np.zeros(len(ptWidths))
    ptSystWidths = np.repeat(0.2, len(ptWidths))
    jpsiRaaVsPt = dfJpsiRaaVsPt["val"]
    jpsiStatRaaVsPt = dfJpsiRaaVsPt["stat"]
    jpsiSystRaaVsPt = dfJpsiRaaVsPt["syst"]

    graStatJpsiRaaVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt), np.array(ptWidths_zero), np.array(jpsiStatRaaVsPt))
    graSystJpsiRaaVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt), np.array(ptSystWidths), np.array(jpsiSystRaaVsPt))

    SetHistStat(graStatJpsiRaaVsPt, 24, ROOT.kBlack)
    SetHistSyst(graSystJpsiRaaVsPt, 20, ROOT.kRed+1)
    graStatJpsiRaaVsPt.SetLineWidth(1) 
    graSystJpsiRaaVsPt.SetLineWidth(2) 


    dfJpsiRaaThuVsPt = pd.read_csv(config["inputs"]["fInJpsiRaaThuVsPt"],sep=r"\s+",engine="python")
    ptMin = pd.to_numeric(dfJpsiRaaThuVsPt["x_min"])
    ptMax = pd.to_numeric(dfJpsiRaaThuVsPt["x_max"])
    jpsiRaaThuVsPt = dfJpsiRaaThuVsPt["val"]
    jpsiRaaThuVsPt_Up = dfJpsiRaaThuVsPt["val_up"]
    jpsiRaaThuVsPt_Low = dfJpsiRaaThuVsPt["val_low"]
    jpsiRaaThuVsPt_Regeneration = pd.to_numeric(dfJpsiRaaThuVsPt["Raa-regeneration"])	
    jpsiRaaThuVsPt_InitialStage = pd.to_numeric(dfJpsiRaaThuVsPt["Raa-initial"])
    ptCenters = (ptMax + ptMin) / 2.
    ptWidths = (ptMax - ptMin) / 2.
    ptSystWidths = np.repeat(1.5, len(ptWidths))

    errUp = jpsiRaaThuVsPt_Up - jpsiRaaThuVsPt
    errLow = jpsiRaaThuVsPt - jpsiRaaThuVsPt_Low

    nPoints = len(ptCenters)

    ptWidthsArray = np.full(nPoints, 1e-8)  # tiny width
    errUpLineArray = np.full(nPoints, 1e-6)
    errDownLineArray = np.full(nPoints, 1e-6)

    graJpsiRaaThuVsPt = ROOT.TGraphAsymmErrors(len(ptCenters),np.array(ptCenters),np.array(jpsiRaaThuVsPt),np.array(ptWidths), np.array(ptWidths),np.array(errLow), np.array(errUp))
    graJpsiRaaThuVsPt_line = ROOT.TGraphAsymmErrors(nPoints,np.array(ptCenters, dtype=float),np.array(jpsiRaaThuVsPt, dtype=float),ptWidthsArray,ptWidthsArray,errDownLineArray,errUpLineArray)
    graJpsiRaaThuVsPt_Regeneration_line = ROOT.TGraphAsymmErrors(nPoints,np.array(ptCenters, dtype=float),np.array(jpsiRaaThuVsPt_Regeneration, dtype=float),ptWidthsArray,ptWidthsArray,errDownLineArray,errUpLineArray)
    graJpsiRaaThuVsPt_InitialStage_line = ROOT.TGraphAsymmErrors(nPoints,np.array(ptCenters, dtype=float),np.array(jpsiRaaThuVsPt_InitialStage, dtype=float),ptWidthsArray,ptWidthsArray,errDownLineArray,errUpLineArray)

    graJpsiRaaThuVsPt.SetFillColorAlpha(ROOT.kCyan, 0.3)
    graJpsiRaaThuVsPt.SetFillStyle(1001)
    graJpsiRaaThuVsPt.SetLineColor(ROOT.kBlack+4)
    graJpsiRaaThuVsPt.SetLineWidth(2) 
    graJpsiRaaThuVsPt_line.SetLineColor(ROOT.kAzure+2)
    graJpsiRaaThuVsPt_line.SetLineWidth(2)
    graJpsiRaaThuVsPt_line.SetLineWidth(3)
    graJpsiRaaThuVsPt_Regeneration_line.SetLineStyle(9)
    graJpsiRaaThuVsPt_Regeneration_line.SetLineColor(ROOT.kAzure+2)
    graJpsiRaaThuVsPt_Regeneration_line.SetLineWidth(3)
    graJpsiRaaThuVsPt_InitialStage_line.SetLineStyle(2)
    graJpsiRaaThuVsPt_InitialStage_line.SetLineColor(ROOT.kAzure+2)
    graJpsiRaaThuVsPt_InitialStage_line.SetLineWidth(3)
 
    

    histGridJpsiRaaVsPt = ROOT.TH2D("histGridJpsiRaaVsPt", ";#it{p}_{T} (GeV/#it{c});#it{R}_{OO}", 100, 0, 8, 100, 0, 1.5)
    canvasJpsiRaaThuVsPt = ROOT.TCanvas("canvasJpsiRaaThuVsPt", "", 800, 600)
    histGridJpsiRaaVsPt.Draw()
    graJpsiRaaThuVsPt.Draw("E3 SAME")
    graJpsiRaaThuVsPt_line.Draw("L SAME")
    graStatJpsiRaaVsPt.Draw("EP SAME")
    graSystJpsiRaaVsPt.Draw("E2P SAME")
    graJpsiRaaThuVsPt_Regeneration_line.Draw("L SAME")
    graJpsiRaaThuVsPt_InitialStage_line.Draw("L SAME")

    lineUnityVsPt = ROOT.TLine(0, 1, 8, 1)
    lineUnityVsPt.SetLineColor(ROOT.kGray+1)
    lineUnityVsPt.SetLineWidth(2)
    lineUnityVsPt.SetLineStyle(ROOT.kDashed)
    lineUnityVsPt.Draw()

    latexTitle = ROOT.TLatex()
    latexTitle.SetNDC()
    latexTitle.SetTextSize(0.05)
    latexTitle.SetTextFont(42)
    latexTitle2 = ROOT.TLatex()
    latexTitle2.SetNDC()
    latexTitle2.SetTextSize(0.042)
    latexTitle2.SetTextFont(42)

    latexTitle.DrawLatex(0.20, 0.86, "ALICE Preliminary, OO,")
    latexTitle.DrawLatex(0.57, 0.86, "#sqrt{#it{s}} = 5.36 TeV")
    latexTitle.DrawLatex(0.20, 0.78, "J/#psi")
    latexTitle.DrawLatex(0.26, 0.78, "#rightarrow")
    latexTitle.DrawLatex(0.30, 0.78, "#mu^{#plus}#mu^{#minus},")
    latexTitle.DrawLatex(0.38, 0.78, "2.5<")
    latexTitle.DrawLatex(0.45, 0.78, "#it{y}")
    latexTitle.DrawLatex(0.47, 0.78, "<4,")
    latexTitle.DrawLatex(0.53, 0.78, "0#minus100%")
    latexTitle2.DrawLatex(0.63, 0.34, "Zhuang et al.")

    legendRaaVsPt_Data = ROOT.TLegend(0.44, 0.26, 0.60, 0.28, " ", "brNDC")
    SetLegend(legendRaaVsPt_Data)
    legendRaaVsPt_Data.SetTextSize(0.042)
    legendRaaVsPt_Data.AddEntry(graSystJpsiRaaVsPt, "Data", "P")
    legendRaaVsPt_Data.SetBorderSize(0) 
    legendRaaVsPt_Data.SetFillStyle(0) 
    legendRaaVsPt_Data.Draw("SAME") 

    legendRaaVsPt_Model = ROOT.TLegend(0.62, 0.19, 0.999, 0.36, " ", "brNDC")
    SetLegend(legendRaaVsPt_Model)
    legendRaaVsPt_Model.SetTextSize(0.042)
    legendRaaVsPt_Model.AddEntry(graJpsiRaaThuVsPt_line, "Total", "L")
    legendRaaVsPt_Model.AddEntry(graJpsiRaaThuVsPt_InitialStage_line, "Initial", "L")
    legendRaaVsPt_Model.AddEntry(graJpsiRaaThuVsPt_Regeneration_line, "Regeneration", "L")
    legendRaaVsPt_Model.SetBorderSize(0) 
    legendRaaVsPt_Model.SetFillStyle(0) 
    legendRaaVsPt_Model.Draw("SAME") 

    xMinBox = 7.8
    xMaxBox = 8.0
    yRef = 1.0
    systCommonTotData = 0.050744457825461095 

    yMinBox = yRef * (1.0 - systCommonTotData)
    yMaxBox = yRef * (1.0 + systCommonTotData)

    boxCommon_Data = ROOT.TBox(xMinBox, yMinBox, xMaxBox, yMaxBox)
    boxCommon_Data.SetFillColorAlpha(ROOT.kRed+1, 1)
    boxCommon_Data.SetLineColor(ROOT.kRed+1)
    boxCommon_Data.SetLineWidth(1)
    boxCommon_Data.Draw("EP same")


    canvasJpsiRaaThuVsPt.Update()
    canvasJpsiRaaThuVsPt.SaveAs("SQM2026/Jpsi_ROO_vs_Pt_vs_THUModel.pdf")
    canvasJpsiRaaThuVsPt.SaveAs("SQM2026/Jpsi_ROO_vs_Pt_vs_THUModel.png")
    input()




if __name__ == '__main__':
    main()