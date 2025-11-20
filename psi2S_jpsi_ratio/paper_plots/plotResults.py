import sys
import argparse
import yaml
import ROOT
import numpy as np
import array
import pandas as pd
import math

sys.path.append('../../utils')
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
        plotResults(inputCfg)

def plotResults(config):
    LoadStyle()
    ROOT.gStyle.SetOptStat(False)

    BrJpsiToMuMu =  0.05961
    errBrJpsiToMuMu = 0.033e-2
    BrPsi2sToMuMu =  0.008
    errBrPsi2sToMuMu = 0.6e-3

    systBrRatio = math.sqrt((errBrPsi2sToMuMu/BrPsi2sToMuMu)**2 + (errBrJpsiToMuMu/BrJpsiToMuMu)**2)
    print(f'Run 3 BR systematic uncertainty = {systBrRatio}')
    systBrRatioRun2 = 0.11412162 # taken form Table 1 of https://arxiv.org/pdf/1702.00557

    latexTitle = ROOT.TLatex()
    latexTitle.SetNDC()
    latexTitle.SetTextSize(0.05)
    latexTitle.SetTextFont(42)

    fInFwdRatioVsPt = config["inputs"]["fInFwdRatioVsPt"]
    fInFwdRatioVsRap = config["inputs"]["fInFwdRatioVsRap"]
    fInFwdRatioVsPtRun2 = config["inputs"]["fInFwdRatioVsPtRun2"]
    fInFwdRatioVsRapRun2 = config["inputs"]["fInFwdRatioVsRapRun2"]

    # ------------ #
    # Run 2 results
    # ------------ #

    # pt dependence
    ptCentersRun2, ptWidthsRun2, fwdPsi2sOverJpsiVsPtRun2, statFwdPsi2sOverJpsiVsPtRun2, systFwdPsi2sOverJpsiVsPtRun2 = ExtractFromYaml(fInFwdRatioVsPtRun2)
    ptSystWidthsRun2 = np.repeat(0.2, len(ptCentersRun2))
    systFwdPsi2sOverJpsiVsPtRun2 = np.array(fwdPsi2sOverJpsiVsPtRun2) * np.sqrt((np.array(systFwdPsi2sOverJpsiVsPtRun2)/np.array(fwdPsi2sOverJpsiVsPtRun2))**2 - systBrRatioRun2**2)
    
    graStatFwdPsi2sOverJpsiVsPtRun2 = ROOT.TGraphErrors(len(ptCentersRun2), np.array(ptCentersRun2), np.array(fwdPsi2sOverJpsiVsPtRun2), np.array(ptWidthsRun2), np.array(statFwdPsi2sOverJpsiVsPtRun2))
    graSystFwdPsi2sOverJpsiVsPtRun2 = ROOT.TGraphErrors(len(ptCentersRun2), np.array(ptCentersRun2), np.array(fwdPsi2sOverJpsiVsPtRun2), np.array(ptSystWidthsRun2), np.array(systFwdPsi2sOverJpsiVsPtRun2))
    SetGraStat(graStatFwdPsi2sOverJpsiVsPtRun2, 20, ROOT.kGray+2, 1, 1, 0.7)
    SetGraSyst(graSystFwdPsi2sOverJpsiVsPtRun2, 20, ROOT.kGray+2, 1, 1, 0.7)

    # rapidity dependence
    rapCentersRun2, rapWidthsRun2, fwdPsi2sOverJpsiVsRapRun2, statFwdPsi2sOverJpsiVsRapRun2, systFwdPsi2sOverJpsiVsRapRun2 = ExtractFromYaml(fInFwdRatioVsRapRun2)
    rapSystWidthsRun2 = np.repeat(0.05, len(rapCentersRun2))
    systFwdPsi2sOverJpsiVsRapRun2 = np.array(fwdPsi2sOverJpsiVsRapRun2) * np.sqrt((np.array(systFwdPsi2sOverJpsiVsRapRun2)/np.array(fwdPsi2sOverJpsiVsRapRun2))**2 - systBrRatioRun2**2)
    
    graStatFwdPsi2sOverJpsiVsRapRun2 = ROOT.TGraphErrors(len(rapCentersRun2), np.array(rapCentersRun2), np.array(fwdPsi2sOverJpsiVsRapRun2), np.array(rapWidthsRun2), np.array(statFwdPsi2sOverJpsiVsRapRun2))
    graSystFwdPsi2sOverJpsiVsRapRun2 = ROOT.TGraphErrors(len(rapCentersRun2), np.array(rapCentersRun2), np.array(fwdPsi2sOverJpsiVsRapRun2), np.array(rapSystWidthsRun2), np.array(systFwdPsi2sOverJpsiVsRapRun2))
    SetGraStat(graStatFwdPsi2sOverJpsiVsRapRun2, 20, ROOT.kGray+2, 1, 1, 0.7)
    SetGraSyst(graSystFwdPsi2sOverJpsiVsRapRun2, 20, ROOT.kGray+2, 1, 1, 0.7)

    # ------------ #
    # Run 3 results
    # ------------ #

    # pt dependence
    dfFwdPsi2sOverJpsiVsPt = pd.read_csv(fInFwdRatioVsPt, sep=' ')
    ptMin = dfFwdPsi2sOverJpsiVsPt["x_min"]
    ptMax = dfFwdPsi2sOverJpsiVsPt["x_max"]
    ptCenters = (ptMax + ptMin) / 2.
    ptWidths = (ptMax - ptMin) / 2.
    ptSystWidths = np.repeat(0.2, len(ptWidths))
    ptEdges = np.append(ptMin.to_numpy(), ptMax.to_numpy()[len(ptMax)-1])
    fwdPsi2sOverJpsiVsPt = dfFwdPsi2sOverJpsiVsPt["val"]
    statFwdPsi2sOverJpsiVsPt = dfFwdPsi2sOverJpsiVsPt["stat"]
    systFwdPsi2sOverJpsiVsPt = dfFwdPsi2sOverJpsiVsPt["syst"]

    graStatFwdPsi2sOverJpsiVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(fwdPsi2sOverJpsiVsPt), np.array(ptWidths), np.array(statFwdPsi2sOverJpsiVsPt))
    graSystFwdPsi2sOverJpsiVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(fwdPsi2sOverJpsiVsPt), np.array(ptSystWidths), np.array(systFwdPsi2sOverJpsiVsPt))

    SetGraStat(graStatFwdPsi2sOverJpsiVsPt, 24, ROOT.kBlack)
    SetGraSyst(graSystFwdPsi2sOverJpsiVsPt, 20, ROOT.kRed+1, 1, 2)

    # Rapidity dependence
    dfFwdPsi2sOverJpsiVsRap = pd.read_csv(fInFwdRatioVsRap, sep=' ')
    rapMin = dfFwdPsi2sOverJpsiVsRap["x_min"]
    rapMax = dfFwdPsi2sOverJpsiVsRap["x_max"]
    rapCenters = (rapMax + rapMin) / 2.
    rapWidths = (rapMax - rapMin) / 2.
    rapSystWidths = np.repeat(0.05, len(rapWidths))
    rapEdges = np.append(rapMin.to_numpy(), rapMax.to_numpy()[len(rapMax)-1])
    fwdPsi2sOverJpsiVsRap = dfFwdPsi2sOverJpsiVsRap["val"]
    statFwdPsi2sOverJpsiVsRap = dfFwdPsi2sOverJpsiVsRap["stat"]
    systFwdPsi2sOverJpsiVsRap = dfFwdPsi2sOverJpsiVsRap["syst"]

    graStatFwdPsi2sOverJpsiVsRap = ROOT.TGraphErrors(len(rapCenters), np.array(rapCenters), np.array(fwdPsi2sOverJpsiVsRap), np.array(rapWidths), np.array(statFwdPsi2sOverJpsiVsRap))
    graSystFwdPsi2sOverJpsiVsRap = ROOT.TGraphErrors(len(rapCenters), np.array(rapCenters), np.array(fwdPsi2sOverJpsiVsRap), np.array(rapSystWidths), np.array(systFwdPsi2sOverJpsiVsRap))

    SetGraStat(graStatFwdPsi2sOverJpsiVsRap, 24, ROOT.kBlack)
    SetGraSyst(graSystFwdPsi2sOverJpsiVsRap, 20, ROOT.kRed+1, 1, 2)

    fInMidPsi2sOverJpsiVsPt = ROOT.TFile("mid_ratio_vs_pt.root", "READ")
    histStatMidPsi2sOverJpsiVsPt = fInMidPsi2sOverJpsiVsPt.Get("h_ratio_pt_corr")
    histSystMidPsi2sOverJpsiVsPt = fInMidPsi2sOverJpsiVsPt.Get("h_syst_all")

    ptMinMid, ptMaxMid , midPsi2sOverJpsiVsPt, statMidPsi2sOverJpsiVsPt, systMidPsi2sOverJpsiVsPt = [], [], [], [], []
    for iPt in range(histStatMidPsi2sOverJpsiVsPt.GetXaxis().GetNbins()):
        ptMinMid.append(histStatMidPsi2sOverJpsiVsPt.GetBinLowEdge(iPt+1))
        ptMaxMid.append(histStatMidPsi2sOverJpsiVsPt.GetBinLowEdge(iPt+1) + histStatMidPsi2sOverJpsiVsPt.GetBinWidth(iPt+1))
        midPsi2sOverJpsiVsPt.append(histStatMidPsi2sOverJpsiVsPt.GetBinContent(iPt+1))
        statMidPsi2sOverJpsiVsPt.append(histStatMidPsi2sOverJpsiVsPt.GetBinError(iPt+1))
        systMidPsi2sOverJpsiVsPt.append(histSystMidPsi2sOverJpsiVsPt.GetBinContent(iPt+1))

    ptCentersMid = (np.array(ptMaxMid) + np.array(ptMinMid)) / 2.
    ptWidthsMid = (np.array(ptMaxMid) - np.array(ptMinMid)) / 2.
    ptSystWidthsMid = np.repeat(0.2, len(ptWidthsMid))

    graStatMidPsi2sOverJpsiVsPt = ROOT.TGraphErrors(len(ptCentersMid), np.array(ptCentersMid), np.array(midPsi2sOverJpsiVsPt), np.array(ptWidthsMid), np.array(statMidPsi2sOverJpsiVsPt))
    graSystMidPsi2sOverJpsiVsPt = ROOT.TGraphErrors(len(ptCentersMid), np.array(ptCentersMid), np.array(midPsi2sOverJpsiVsPt), np.array(ptSystWidthsMid), np.array(systMidPsi2sOverJpsiVsPt))

    SetGraStat(graStatMidPsi2sOverJpsiVsPt, 24, ROOT.kBlack)
    SetGraSyst(graSystMidPsi2sOverJpsiVsPt, 20, ROOT.kAzure+4, 1, 2)



    fInMidPsi2sOverJpsiVsRap = ROOT.TFile("mid_ratio_vs_rap.root", "READ")
    histStatMidPsi2sOverJpsiVsRap = fInMidPsi2sOverJpsiVsRap.Get("h_ratio_pt_corr")
    histSystMidPsi2sOverJpsiVsRap = fInMidPsi2sOverJpsiVsRap.Get("h_syst_all")

    rapMinMid, rapMaxMid , midPsi2sOverJpsiVsRap, statMidPsi2sOverJpsiVsRap, systMidPsi2sOverJpsiVsRap = [], [], [], [], []
    for iRap in range(histStatMidPsi2sOverJpsiVsRap.GetXaxis().GetNbins()):
        rapMinMid.append(histStatMidPsi2sOverJpsiVsRap.GetBinLowEdge(iRap+1))
        rapMaxMid.append(histStatMidPsi2sOverJpsiVsRap.GetBinLowEdge(iRap+1) + histStatMidPsi2sOverJpsiVsRap.GetBinWidth(iRap+1))
        midPsi2sOverJpsiVsRap.append(histStatMidPsi2sOverJpsiVsRap.GetBinContent(iRap+1))
        statMidPsi2sOverJpsiVsRap.append(histStatMidPsi2sOverJpsiVsRap.GetBinError(iRap+1))
        systMidPsi2sOverJpsiVsRap.append(histSystMidPsi2sOverJpsiVsRap.GetBinContent(iRap+1))

    rapCentersMid = (np.array(rapMaxMid) + np.array(rapMinMid)) / 2.
    rapWidthsMid = (np.array(rapMaxMid) - np.array(rapMinMid)) / 2.
    rapSystWidthsMid = np.repeat(0.05, len(rapWidthsMid))

    graStatMidPsi2sOverJpsiVsRap = ROOT.TGraphErrors(len(rapCentersMid), np.array(rapCentersMid), np.array(midPsi2sOverJpsiVsRap), np.array(rapWidthsMid), np.array(statMidPsi2sOverJpsiVsRap))
    graSystMidPsi2sOverJpsiVsRap = ROOT.TGraphErrors(len(rapCentersMid), np.array(rapCentersMid), np.array(midPsi2sOverJpsiVsRap), np.array(rapSystWidthsMid), np.array(systMidPsi2sOverJpsiVsRap))

    SetGraStat(graStatMidPsi2sOverJpsiVsRap, 24, ROOT.kBlack)
    SetGraSyst(graSystMidPsi2sOverJpsiVsRap, 20, ROOT.kAzure+4, 1, 2)

    # ------------ #
    # Plot results
    # ------------ #    

    # pt dependence vs Run 2
    canvasFwdPsi2sOverJpsiVsPtVsRun2 = ROOT.TCanvas("canvasFwdPsi2sOverJpsiVsPtVsRun2", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    #ROOT.gPad.SetLogy(True)
    histGridFwdPsi2sOverJpsiVsPtVsRun2 = ROOT.TH2D("histGridFwdPsi2sOverJpsiVsPtVsRun2", ";#it{p}_{T} (GeV/#it{c});#sigma_{#psi(2S)} / #sigma_{J/#psi}", 100, 0, 20, 100, 0.0, 0.8)
    histGridFwdPsi2sOverJpsiVsPtVsRun2.Draw()
    graStatFwdPsi2sOverJpsiVsPtRun2.Draw("EP SAME")
    graSystFwdPsi2sOverJpsiVsPtRun2.Draw("E2P SAME")
    graSystFwdPsi2sOverJpsiVsPt.Draw("E2P SAME")
    graStatFwdPsi2sOverJpsiVsPt.Draw("EP SAME")

    legendFwdPsi2sOverJpsiVsPtVsRun2 = ROOT.TLegend(0.20, 0.60, 0.40, 0.80, " ", "brNDC")
    SetLegend(legendFwdPsi2sOverJpsiVsPtVsRun2)
    legendFwdPsi2sOverJpsiVsPtVsRun2.SetTextSize(0.045)
    legendFwdPsi2sOverJpsiVsPtVsRun2.AddEntry(graSystFwdPsi2sOverJpsiVsPtRun2, "#sqrt{#it{s}} = 13 TeV", "FP")
    legendFwdPsi2sOverJpsiVsPtVsRun2.AddEntry(graSystFwdPsi2sOverJpsiVsPt, "#sqrt{#it{s}} = 13.6 TeV", "FP")
    legendFwdPsi2sOverJpsiVsPtVsRun2.Draw("SAME")

    latexTitle.DrawLatex(0.20, 0.87, "ALICE, pp collisions")
    latexTitle.DrawLatex(0.20, 0.80, "J/#psi, #psi(2S) #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    latexTitle.DrawLatex(0.52, 0.20, "#pm 7.5% BR unc. not shown")
    canvasFwdPsi2sOverJpsiVsPtVsRun2.Update()
    canvasFwdPsi2sOverJpsiVsPtVsRun2.SaveAs("FwdPsi2sOverJpsiVsPtVsRun2.pdf")

    # rapidity dependence vs Run 2
    canvasFwdPsi2sOverJpsiVsRapVsRun2 = ROOT.TCanvas("canvasFwdPsi2sOverJpsiVsRapVsRun2", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    #ROOT.gPad.SetLogy(True)
    histGridFwdPsi2sOverJpsiVsRapVsRun2 = ROOT.TH2D("histGridFwdPsi2sOverJpsiVsRapVsRun2", ";#it{y};#sigma_{#psi(2S)} / #sigma_{J/#psi}", 100, 2.5, 4, 100, 0.0, 0.35)
    histGridFwdPsi2sOverJpsiVsRapVsRun2.Draw()
    graStatFwdPsi2sOverJpsiVsRapRun2.Draw("EP SAME")
    graSystFwdPsi2sOverJpsiVsRapRun2.Draw("E2P SAME")
    graSystFwdPsi2sOverJpsiVsRap.Draw("E2P SAME")
    graStatFwdPsi2sOverJpsiVsRap.Draw("EP SAME")

    legendFwdPsi2sOverJpsiVsRapVsRun2 = ROOT.TLegend(0.20, 0.60, 0.40, 0.80, " ", "brNDC")
    SetLegend(legendFwdPsi2sOverJpsiVsRapVsRun2)
    legendFwdPsi2sOverJpsiVsRapVsRun2.SetTextSize(0.045)
    legendFwdPsi2sOverJpsiVsRapVsRun2.AddEntry(graSystFwdPsi2sOverJpsiVsRapRun2, "#sqrt{#it{s}} = 13 TeV", "FP")
    legendFwdPsi2sOverJpsiVsRapVsRun2.AddEntry(graSystFwdPsi2sOverJpsiVsRap, "#sqrt{#it{s}} = 13.6 TeV", "FP")
    legendFwdPsi2sOverJpsiVsRapVsRun2.Draw("SAME")

    latexTitle.DrawLatex(0.20, 0.87, "ALICE, pp collisions")
    latexTitle.DrawLatex(0.20, 0.80, "J/#psi, #psi(2S) #rightarrow #mu^{#plus}#mu^{#minus}, #it{p}_{T} < 20 GeV/#it{c}")
    latexTitle.DrawLatex(0.52, 0.20, "#pm 7.5% BR unc. not shown")
    canvasFwdPsi2sOverJpsiVsRapVsRun2.Update()
    canvasFwdPsi2sOverJpsiVsRapVsRun2.SaveAs("FwdPsi2sOverJpsiVsRapVsRun2.pdf")

    # pt dependence vs midrapidity
    canvasFwdPsi2sOverJpsiVsPtVsMid = ROOT.TCanvas("canvasFwdPsi2sOverJpsiVsPtVsMid", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    #ROOT.gPad.SetLogy(True)
    histGridFwdPsi2sOverJpsiVsPtVsMid = ROOT.TH2D("histGridFwdPsi2sOverJpsiVsPtVsMid", ";#it{y};#sigma_{#psi(2S)} / #sigma_{J/#psi}", 100, 0, 20, 100, 0, 0.8)
    histGridFwdPsi2sOverJpsiVsPtVsMid.Draw()
    graSystMidPsi2sOverJpsiVsPt.Draw("E2P SAME")
    graStatMidPsi2sOverJpsiVsPt.Draw("EP SAME")
    graSystFwdPsi2sOverJpsiVsPt.Draw("E2P SAME")
    graStatFwdPsi2sOverJpsiVsPt.Draw("EP SAME")

    legendFwdPsi2sOverJpsiVsPtVsMid = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendFwdPsi2sOverJpsiVsPtVsMid)
    legendFwdPsi2sOverJpsiVsPtVsMid.SetTextSize(0.045)
    #legendFwdPsi2sOverJpsiVsPtVsMid.AddEntry(graSystFwdPsi2sOverJpsiVsPtRun2, "#sqrt{#it{s}} = 13 TeV", "FP")
    legendFwdPsi2sOverJpsiVsPtVsMid.AddEntry(graSystFwdPsi2sOverJpsiVsPt, "J/#psi, #psi(2S) #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4", "FP")
    legendFwdPsi2sOverJpsiVsPtVsMid.AddEntry(graSystMidPsi2sOverJpsiVsPt, "J/#psi, #psi(2S) #rightarrow e^{#plus}e^{#minus}, |#it{y}| < 0.8", "FP")
    legendFwdPsi2sOverJpsiVsPtVsMid.Draw("SAME")

    latexTitle.DrawLatex(0.20, 0.87, "ALICE, pp #sqrt{#it{s}} = 13.6 TeV")
    #latexTitle.DrawLatex(0.20, 0.80, "J/#psi, #psi(2S) #rightarrow #mu^{#plus}#mu^{#minus}, #it{p}_{T} < 20 GeV/#it{c}")
    latexTitle.DrawLatex(0.52, 0.20, "#pm 7.5% BR unc. not shown")
    canvasFwdPsi2sOverJpsiVsPtVsMid.Update()
    canvasFwdPsi2sOverJpsiVsPtVsMid.SaveAs("FwdPsi2sOverJpsiVsPtVsMid.pdf")

    # rapidity dependence vs midrapidity
    canvasFwdPsi2sOverJpsiVsRapVsMid = ROOT.TCanvas("canvasFwdPsi2sOverJpsiVsRapVsMid", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    #ROOT.gPad.SetLogy(True)
    histGridFwdPsi2sOverJpsiVsRapVsMid = ROOT.TH2D("histGridFwdPsi2sOverJpsiVsRapVsMid", ";#it{y};#sigma_{#psi(2S)} / #sigma_{J/#psi}", 100, -1, 4.5, 100, 0.0, 0.25)
    histGridFwdPsi2sOverJpsiVsRapVsMid.Draw()
    graSystMidPsi2sOverJpsiVsRap.Draw("E2P SAME")
    graStatMidPsi2sOverJpsiVsRap.Draw("EP SAME")
    graSystFwdPsi2sOverJpsiVsRap.Draw("E2P SAME")
    graStatFwdPsi2sOverJpsiVsRap.Draw("EP SAME")

    legendFwdPsi2sOverJpsiVsRapVsRun2 = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendFwdPsi2sOverJpsiVsRapVsRun2)
    legendFwdPsi2sOverJpsiVsRapVsRun2.SetTextSize(0.045)
    #legendFwdPsi2sOverJpsiVsRapVsRun2.AddEntry(graSystFwdPsi2sOverJpsiVsRapRun2, "#sqrt{#it{s}} = 13 TeV", "FP")
    legendFwdPsi2sOverJpsiVsRapVsRun2.AddEntry(graSystFwdPsi2sOverJpsiVsRap, "J/#psi, #psi(2S) #rightarrow #mu^{#plus}#mu^{#minus}, #it{p}_{T} < 20 GeV/#it{c}", "FP")
    legendFwdPsi2sOverJpsiVsRapVsRun2.AddEntry(graSystMidPsi2sOverJpsiVsRap, "J/#psi, #psi(2S) #rightarrow e^{#plus}e^{#minus}, #it{p}_{T} < 16 GeV/#it{c}", "FP")
    legendFwdPsi2sOverJpsiVsRapVsRun2.Draw("SAME")


    latexTitle.DrawLatex(0.20, 0.87, "ALICE, pp #sqrt{#it{s}} = 13.6 TeV")
    #latexTitle.DrawLatex(0.20, 0.80, "J/#psi, #psi(2S) #rightarrow #mu^{#plus}#mu^{#minus}, #it{p}_{T} < 20 GeV/#it{c}")
    latexTitle.DrawLatex(0.52, 0.20, "#pm 7.5% BR unc. not shown")
    canvasFwdPsi2sOverJpsiVsRapVsMid.Update()
    canvasFwdPsi2sOverJpsiVsRapVsMid.SaveAs("FwdPsi2sOverJpsiVsRapVsMid.pdf")

    input()
    exit()

"""
def raa(config):
    LoadStyle()
    ROOT.gStyle.SetOptStat(False)
    BrJpsiToMuMu =  0.05961
    errBrJpsiToMuMu = 0.033e-2
    BrPsi2sToMuMu =  0.008
    errBrPsi2sToMuMu = 0.6e-3

    latexTitle = ROOT.TLatex()
    latexTitle.SetNDC()
    latexTitle.SetTextSize(0.05)
    latexTitle.SetTextFont(42)

    deltaRap = config["inputs"]["deltaRap"]
    Taa = config["inputs"]["Taa"]
    TaaVsCentr = config["inputs"]["TaaVsCentr"]

    print("***** Compute corrected yield *****")
    # Lumi computed with normalization.C  + check_normalization
    #lumi = 0.170285 # pb-1
    fInLumi = ROOT.TFile(config["inputs"]["fInLumi"], "READ")
    histLumi = fInLumi.Get("histLumi")
    lumi = histLumi.GetBinContent(1)
    histNevMinBias = fInLumi.Get("histNevMinBias")
    nevMinBias = histNevMinBias.GetBinContent(1)
    print("N. min bias events = ", nevMinBias)

    histNevMinBiasCentr = fInLumi.Get("histNeventMinBias_Centr")
    nCentrBins = histNevMinBiasCentr.GetNbinsX()
    centrLabels = ["0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-90%", "0-90%", "0-100%"]
    nevCentrDict = {}  # Dictionary for nEvts vs Centr
    for i in range(nCentrBins):
        nEvts = histNevMinBiasCentr.GetBinContent(i+1)

        if i < len(centrLabels):
            label = centrLabels[i]
        else:
            label = histNevMinBiasCentr.GetXaxis().GetBinLabel(i+1)

        nevCentrDict[label] = nEvts
        print(f"N. events min bias {label} = {nEvts}")

    print("***** Extract pp reference [mub] *****")
    fInPPrefVsPt = ROOT.TFile(config["inputs"]["fInPPredVsPt"], "READ")
    fInPPrefVsPt.ls()
    histStatPPrefXsec = fInPPrefVsPt.Get("histStatJpsiXsecInterp")
    histSystPPrefXsec = fInPPrefVsPt.Get("histSystJpsiXsecInterp")
    histStatPPrefXsecVsSqrts = fInPPrefVsPt.Get("histStatJpsiXsecInterpVsSqrts")
    histSystPPrefXsecVsSqrts = fInPPrefVsPt.Get("histSystJpsiXsecInterpVsSqrts")
    
    SetHistStat(histStatPPrefXsec, 20, ROOT.kAzure+4)
    SetHistSyst(histSystPPrefXsec, 20, ROOT.kAzure+4)

    SetHistStat(histStatPPrefXsecVsSqrts, 20, ROOT.kAzure+4)
    SetHistSyst(histSystPPrefXsecVsSqrts, 20, ROOT.kAzure+4)

    canvasPPrefXsecVsSqrts = ROOT.TCanvas("canvasPPrefXsecVsSqrts", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histStatPPrefXsecVsSqrts.Draw("EP")
    histSystPPrefXsecVsSqrts.Draw("E2P SAME")
    canvasPPrefXsecVsSqrts.Update()

    jpsiPPrefXsecVsCentr = histStatPPrefXsecVsSqrts.GetBinContent(1)
    jpsiStatPPrefXsecVsCentr = histStatPPrefXsecVsSqrts.GetBinError(1)
    jpsiSystPPrefXsecVsCentr = histSystPPrefXsecVsSqrts.GetBinError(1)

    canvasPPrefXsec = ROOT.TCanvas("canvasPPrefXsec", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histStatPPrefXsec.Draw("EP")
    histSystPPrefXsec.Draw("E2P SAME")
    canvasPPrefXsec.Update()


    jpsiPPrefXsecVsPt, jpsiStatPPrefXsecVsPt, jpsiSystPPrefXsecVsPt = [], [], []

    for iPt in range(0, histStatPPrefXsec.GetXaxis().GetNbins()):
        jpsiPPrefXsecVsPt.append(histStatPPrefXsec.GetBinContent(iPt+1))
        jpsiStatPPrefXsecVsPt.append(histStatPPrefXsec.GetBinError(iPt+1))
        jpsiSystPPrefXsecVsPt.append(histSystPPrefXsec.GetBinError(iPt+1))

    print("-------- Centrality dependence --------")
    dfJpsiRawYieldVsCentr = pd.read_csv(config["inputs"]["fInRawYieldVsCentr"], sep=' ')
    centrMin = dfJpsiRawYieldVsCentr["x_min"]
    centrMax = dfJpsiRawYieldVsCentr["x_max"]
    centrCenters = (centrMax + centrMin) / 2.
    centrWidths = (centrMax - centrMin) / pt = np.repeat(0.2, len(centrWidths))
    centrEdges = np.append(centrMin.to_numpy(), centrMax.to_numpy()[len(centrMax)-1])
    jpsiRawYieldVsCentr = dfJpsiRawYieldVsCentr["val"] / deltaRap
    jpsiStatRawYieldVsCentr = dfJpsiRawYieldVsCentr["stat"] / deltaRap
    jpsiSystRawYieldVsCentr = dfJpsiRawYieldVsCentr["syst"] / deltaRap

    histStatRawYieldVsCentr = ROOT.TH1D("histStatRawYieldVsCentr", ";Centrality (%);dN/d#it{y}", len(centrEdges)-1, centrEdges)
    histSystRawYieldVsCentr = ROOT.TH1D("histSystRawYieldVsCentr", ";Centrality (%);dN/d#it{y}", len(centrEdges)-1, centrEdges)
    histSystRelRawYieldVsCentr = ROOT.TH1D("histSystRelRawYieldVsCentr", ";Centrality (%);dN/d#it{y}", len(centrEdges)-1, centrEdges)

    SetHistStat(histStatRawYieldVsCentr, 20, ROOT.kAzure+4)
    SetHistSyst(histSystRawYieldVsCentr, 20, ROOT.kAzure+4)

    for iBin in range(0, len(centrMin)):
        histStatRawYieldVsCentr.SetBinContent(iBin+1, jpsiRawYieldVsCentr[iBin])
        histStatRawYieldVsCentr.SetBinError(iBin+1, jpsiStatRawYieldVsCentr[iBin])
        histSystRawYieldVsCentr.SetBinContent(iBin+1, jpsiRawYieldVsCentr[iBin])
        histSystRawYieldVsCentr.SetBinError(iBin+1, jpsiSystRawYieldVsCentr[iBin])
        histSystRelRawYieldVsCentr.SetBinContent(iBin+1, jpsiSystRawYieldVsCentr[iBin] / jpsiRawYieldVsCentr[iBin])

    canvasRawYieldVsCentr = ROOT.TCanvas("canvasRawYieldVsCentr", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histStatRawYieldVsCentr.Draw("EP SAME")
    canvasRawYieldVsCentr.Update()

    print(f'[INFO] Sum centr-diff. spectrum = {histStatRawYieldVsCentr.Integral():.0f}')

    print("***** Extract Axe *****")
    dfJpsiAxeVsCentr = pd.read_csv(config["inputs"]["fInAxeVsCentr"], sep=' ')
    jpsiAxeVsCentr = dfJpsiAxeVsCentr["val"]
    jpsiStatAxeVsCentr = dfJpsiAxeVsCentr["stat"]

    # Axe -> to compute the array 
    nCentrBins_Axe = len(centrMin)
    jpsiAxeVsCentr = np.full(nCentrBins_Axe, jpsiAxeVsCentr)
    jpsiStatAxeVsCentr = np.full(nCentrBins_Axe, jpsiStatAxeVsCentr)

    histStatAxeVsCentr = ROOT.TH1D("histStatAxeVsCentr", ";#it{p}_{T} (GeV/#it{c});A#times#varepsilon", len(centrEdges)-1, centrEdges)

    SetHistStat(histStatAxeVsCentr, 20, ROOT.kAzure+4)

    for iBin in range(0, len(centrMin)):
        histStatAxeVsCentr.SetBinContent(iBin+1, jpsiAxeVsCentr[iBin])
        histStatAxeVsCentr.SetBinError(iBin+1, jpsiStatAxeVsCentr[iBin])

    canvasAxeVsCentr = ROOT.TCanvas("canvasAxeVsCentr", "", 800, 600)
    histStatAxeVsCentr.Draw("EP SAME")
    canvasAxeVsCentr.Update()

    # pp ref -> to compute the array
    jpsiPPrefXsecVsCentr = np.full(nCentrBins_Axe, jpsiPPrefXsecVsCentr)
    jpsiStatPPrefXsecVsCentr = np.full(nCentrBins_Axe, jpsiStatPPrefXsecVsCentr)
    jpsiSystPPrefXsecVsCentr = np.full(nCentrBins_Axe, jpsiSystPPrefXsecVsCentr)

    # nEvts vs centr to use
    selectedLabels = centrLabels[:nCentrBins_Axe]
    nevCentrArray = np.array([nevCentrDict[label] for label in centrLabels[:nCentrBins_Axe]])

    # checking
    print(f"jpsiRawYieldVsCentr: {jpsiRawYieldVsCentr}")
    print(f"jpsiAxeVsCentr: {jpsiAxeVsCentr}")
    for label in selectedLabels:
        print(f"N. events {label} = {nevCentrDict[label]}")
    print(f"TaaVsCentr: {TaaVsCentr}")
    print(f"jpsiPPrefXsecVsCentr: {jpsiPPrefXsecVsCentr}")
    input()

    print("Relative systematic on raw yield extraction")
    jpsiSystRelRawYieldVsCentr = (np.array(jpsiStatRawYieldVsCentr) / np.array(jpsiRawYieldVsCentr))
    print(jpsiSystRelRawYieldVsCentr)

    print("Relative systematic on pp reference")
    jpsi1SystRelPPrefXsecVsCentr = (np.array(jpsiStatPPrefXsecVsCentr) / np.array(jpsiPPrefXsecVsCentr))
    jpsi2SystRelPPrefXsecVsCentr = (np.array(jpsiSystPPrefXsecVsCentr) / np.array(jpsiPPrefXsecVsCentr))
    jpsiSystRelPPrefXsecVsCentr = np.sqrt(jpsi1SystRelPPrefXsecVsCentr**2 + jpsi2SystRelPPrefXsecVsCentr**2)
    print(jpsiSystRelPPrefXsecVsCentr)


    jpsiRaaVsCentr = (jpsiRawYieldVsCentr) / (jpsiAxeVsCentr * BrJpsiToMuMu * nevCentrArray * np.array(TaaVsCentr) * jpsiPPrefXsecVsCentr)
    jpsiStatRaaVsCentr = jpsiRaaVsCentr * (np.array(jpsiStatRawYieldVsCentr) / np.array(jpsiRawYieldVsCentr))
    jpsiSystRaaVsCentr = jpsiRaaVsCentr * np.sqrt(jpsiSystRelRawYieldVsCentr**2 + jpsiSystRelPPrefXsecVsCentr**2)

    graStatRaaVsCentr = ROOT.TGraphErrors(len(centrCenters), np.array(centrCenters), np.array(jpsiRaaVsCentr), np.array(centrWidths), np.array(jpsiStatRaaVsCentr))
    graSystRaaVsCentr = ROOT.TGraphErrors(len(centrCenters), np.array(centrCenters), np.array(jpsiRaaVsCentr), np.array(centrSystWidths), np.array(jpsiSystRaaVsCentr))

    SetGraStat(graStatRaaVsCentr, 20, ROOT.kRed+1)
    SetGraSyst(graSystRaaVsCentr, 20, ROOT.kRed+1)

    lineUnityVsCentr = ROOT.TLine(0., 1., 90., 1.)
    lineUnityVsCentr.SetLineColor(ROOT.kGray+1)
    lineUnityVsCentr.SetLineWidth(2)
    lineUnityVsCentr.SetLineStyle(ROOT.kDashed)

    canvasRaaVsCentrVsRun2 = ROOT.TCanvas("canvasRaaVsCentrVsRun2", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    histGridRaaVsCentr = ROOT.TH2D("histGridRaaVsCentr", ";Centrality (%);#it{R}_{AA}", 100, 0, 90, 100, 0, 2.2)
    histGridRaaVsCentr.Draw()
    lineUnityVsCentr.Draw("SAME")
    #graStatRaaVsCentrRun2.Draw("EP SAME")
    #graSystRaaVsCentrRun2.Draw("E2P SAME")
    graStatRaaVsCentr.Draw("EP SAME")
    graSystRaaVsCentr.Draw("E2P SAME")

    legendRaaVsCentrVsRun2 = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsCentrVsRun2)
    legendRaaVsCentrVsRun2.SetTextSize(0.045)
    legendRaaVsCentrVsRun2.AddEntry(graStatRaaVsCentr, "#sqrt{#it{s}} = 5.36 TeV, OO, 0 < #it{p}_{T} (GeV/#it{c}) < 20 ", "FP")
    #legendRaaVsCentrVsRun2.AddEntry(graStatRaaVsCentrRun2, "#sqrt{#it{s}} = 5.02 TeV, Pb-Pb, 0#minus90%", "FP")
    legendRaaVsCentrVsRun2.Draw("SAME")

    latexTitle.DrawLatex(0.15, 0.87, "ALICE Work in Progress, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    canvasRaaVsCentrVsRun2.Update()
    canvasRaaVsCentrVsRun2.SaveAs("plots/Jpsi_RAA_vs_Centr.pdf")


    print("-------- Pt dependence --------")

    print("***** Extract raw yield *****")
    dfJpsiRawYieldVsPt = pd.read_csv(config["inputs"]["fInRawYieldVsPt"], sep=' ')
    ptMin = dfJpsiRawYieldVsPt["x_min"]
    ptMax = dfJpsiRawYieldVsPt["x_max"]
    ptCenters = (ptMax + ptMin) / 2.
    ptWidths = (ptMax - ptMin) / 2.
    ptSystWidths = np.repeat(0.2, len(ptWidths))
    ptEdges = np.append(ptMin.to_numpy(), ptMax.to_numpy()[len(ptMax)-1])
    jpsiRawYieldVsPt = dfJpsiRawYieldVsPt["val"] / deltaRap
    jpsiStatRawYieldVsPt = dfJpsiRawYieldVsPt["stat"] / deltaRap
    jpsiSystRawYieldVsPt = dfJpsiRawYieldVsPt["syst"] / deltaRap

    histStatRawYieldVsPt = ROOT.TH1D("histStatRawYieldVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdges)-1, ptEdges)
    histSystRawYieldVsPt = ROOT.TH1D("histSystRawYieldVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdges)-1, ptEdges)
    histSystRelRawYieldVsPt = ROOT.TH1D("histSystRelRawYieldVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdges)-1, ptEdges)

    SetHistStat(histStatRawYieldVsPt, 20, ROOT.kAzure+4)
    SetHistSyst(histSystRawYieldVsPt, 20, ROOT.kAzure+4)

    for iBin in range(0, len(ptMin)):
        histStatRawYieldVsPt.SetBinContent(iBin+1, jpsiRawYieldVsPt[iBin])
        histStatRawYieldVsPt.SetBinError(iBin+1, jpsiStatRawYieldVsPt[iBin])
        histSystRawYieldVsPt.SetBinContent(iBin+1, jpsiRawYieldVsPt[iBin])
        histSystRawYieldVsPt.SetBinError(iBin+1, jpsiSystRawYieldVsPt[iBin])
        histSystRelRawYieldVsPt.SetBinContent(iBin+1, jpsiSystRawYieldVsPt[iBin] / jpsiRawYieldVsPt[iBin])

    canvasRawYieldVsPt = ROOT.TCanvas("canvasRawYieldVsPt", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    #histSystRawYieldInt.Draw("E2P SAME")
    histStatRawYieldVsPt.Draw("EP SAME")
    histSystRawYieldVsPt.Draw("E2P SAME")
    canvasRawYieldVsPt.Update()

    print(f'[INFO] Sum pt-diff. spectrum = {histStatRawYieldVsPt.Integral():.0f}')

    print("***** Extract Axe *****")
    dfJpsiAxeVsPt = pd.read_csv(config["inputs"]["fInAxeVsPt"], sep=' ')
    jpsiAxeVsPt = dfJpsiAxeVsPt["val"]
    jpsiStatAxeVsPt = dfJpsiAxeVsPt["stat"]

    histStatAxeVsPt = ROOT.TH1D("histStatAxeVsPt", ";#it{p}_{T} (GeV/#it{c});A#times#varepsilon", len(ptEdges)-1, ptEdges)

    SetHistStat(histStatAxeVsPt, 20, ROOT.kAzure+4)

    for iBin in range(0, len(ptMin)):
        histStatAxeVsPt.SetBinContent(iBin+1, jpsiAxeVsPt[iBin])
        histStatAxeVsPt.SetBinError(iBin+1, jpsiStatAxeVsPt[iBin])

    canvasAxeVsPt = ROOT.TCanvas("canvasAxeVsPt", "", 800, 600)
    histStatAxeVsPt.Draw("EP SAME")
    canvasAxeVsPt.Update()

    systRelLumi = 0.1
    histSystRelLumiVsPt = ROOT.TH1D("histSystRelLumiVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdges)-1, ptEdges)

    for iBin in range(0, len(ptMin)):
        histSystRelLumiVsPt.SetBinContent(iBin+1, systRelLumi)
    print("Syst. Rel. Luminosity = ", systRelLumi, " pb-1")

    systRelBrJpsiToMuMu = errBrJpsiToMuMu / BrJpsiToMuMu
    histSystRelBrJpsiToMuMuVsPt = ROOT.TH1D("histSystRelBrJpsiToMuMuVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdges)-1, ptEdges)

    for iBin in range(0, len(ptMin)):
        histSystRelBrJpsiToMuMuVsPt.SetBinContent(iBin+1, systRelBrJpsiToMuMu)
    print("Syst. Rel. BR Jpsi->mumu = ", systRelBrJpsiToMuMu)


    systRelTrackingEff = 0.01
    histSystRelTrackingEffVsPt = ROOT.TH1D("histSystRelTrackingEffVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdges)-1, ptEdges)

    for iBin in range(0, len(ptMin)):
        histSystRelTrackingEffVsPt.SetBinContent(iBin+1, systRelTrackingEff)
    print("Syst. Rel. Tracking efficiency = ", systRelTrackingEff)

    systRelMatchingEff = 0.03
    histSystRelMatchingEffVsPt = ROOT.TH1D("histSystRelMatchingEffVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdges)-1, ptEdges)

    for iBin in range(0, len(ptMin)):
        histSystRelMatchingEffVsPt.SetBinContent(iBin+1, systRelMatchingEff)
    print("Syst. Rel. Matching efficiency = ", systRelMatchingEff)


    systRelMcRealistic = 0.02
    histSystRelMcRealisticVsPt = ROOT.TH1D("histSystRelMcRealisticVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdges)-1, ptEdges)

    for iBin in range(0, len(ptMin)):
        histSystRelMcRealisticVsPt.SetBinContent(iBin+1, systRelMcRealistic)
    print("Syst. Rel. MC realisticness = ", systRelMcRealistic)


    print("Relative systematic on raw yield extraction")
    jpsiSystRelRawYieldVsPt = (np.array(jpsiStatRawYieldVsPt) / np.array(jpsiRawYieldVsPt))
    print(jpsiSystRelRawYieldVsPt)

    print("Relative systematic on pp reference")
    jpsi1SystRelPPrefXsecVsPt = (np.array(jpsiStatPPrefXsecVsPt) / np.array(jpsiPPrefXsecVsPt))
    jpsi2SystRelPPrefXsecVsPt = (np.array(jpsiSystPPrefXsecVsPt) / np.array(jpsiPPrefXsecVsPt))
    jpsiSystRelPPrefXsecVsPt = np.sqrt(jpsi1SystRelPPrefXsecVsPt**2 + jpsi2SystRelPPrefXsecVsPt**2)
    print(jpsiSystRelPPrefXsecVsPt)

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    # > > > Pt dependence < < < #
    jpsiXsecVsPt = (jpsiRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * lumi * 1e6 * (2 * ptWidths))
    jpsiStatXsecVsPt = (jpsiStatRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * lumi * 1e6 * (2 * ptWidths))
    jpsiSystXsecVsPt = (jpsiSystRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * lumi * 1e6 * (2 * ptWidths))

    nevMinBias_0_90 = nevCentrDict["0-90%"]
    jpsiRaaVsPt = (jpsiRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * nevMinBias_0_90 * (2 * ptWidths) * Taa * np.array(jpsiPPrefXsecVsPt))
    jpsiStatRaaVsPt = jpsiRaaVsPt * (np.array(jpsiStatRawYieldVsPt) / np.array(jpsiRawYieldVsPt))
    jpsiSystRaaVsPt = jpsiRaaVsPt * np.sqrt(jpsiSystRelRawYieldVsPt**2 + jpsiSystRelPPrefXsecVsPt**2)

    # Add all systematics contributions
    #jpsiSystXsecVsPt = jpsiXsecVsPt * np.sqrt((jpsiSystXsecVsPt / jpsiXsecVsPt)**2 + systRelBrJpsiToMuMu**2 + systRelTrackingEff**2 + systRelMatchingEff**2)

    graStatXsecVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiXsecVsPt), np.array(ptWidths), np.array(jpsiStatXsecVsPt))
    graSystXsecVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiXsecVsPt), np.array(ptSystWidths), np.array(jpsiSystXsecVsPt))

    graStatRaaVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt), np.array(ptWidths), np.array(jpsiStatRaaVsPt))
    graSystRaaVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt), np.array(ptSystWidths), np.array(jpsiSystRaaVsPt))

    SetGraStat(graStatXsecVsPt, 20, ROOT.kGreen+2)
    SetGraSyst(graSystXsecVsPt, 20, ROOT.kGreen+2)

    SetGraStat(graStatRaaVsPt, 20, ROOT.kGreen+2)
    SetGraSyst(graSystRaaVsPt, 20, ROOT.kGreen+2)


    canvasXsecVsPt = ROOT.TCanvas("canvasXsecVsPt", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    ROOT.gPad.SetLogy(True)
    histGridXsecVsPt = ROOT.TH2D("histGridXsecVsPt", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}", 100, 0, 12, 100, 0, 300)
    histGridXsecVsPt.Draw()
    graStatXsecVsPt.Draw("EP SAME")
    graSystXsecVsPt.Draw("EP SAME")
    canvasXsecVsPt.Update()

    histStatXsecVsPt = ROOT.TH1D("histStatXsecVsPt", ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{y} d#it{p}_{T} (#mub / GeV/#it{c})", len(ptEdges)-1, ptEdges)
    histSystXsecVsPt = ROOT.TH1D("histSystXsecVsPt", ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{y} d#it{p}_{T} (#mub / GeV/#it{c})", len(ptEdges)-1, ptEdges)

    SetHistStat(histStatXsecVsPt, 20, ROOT.kRed+1)
    SetHistSyst(histSystXsecVsPt, 20, ROOT.kRed+1)

    for iBin in range(0, len(ptMin)):
        histStatXsecVsPt.SetBinContent(iBin+1, jpsiXsecVsPt[iBin])
        histStatXsecVsPt.SetBinError(iBin+1, jpsiStatXsecVsPt[iBin])
        histSystXsecVsPt.SetBinContent(iBin+1, jpsiXsecVsPt[iBin])
        histSystXsecVsPt.SetBinError(iBin+1, jpsiSystXsecVsPt[iBin])

    # --- Canvas setup ---
    canvasRaaVsPt = ROOT.TCanvas("canvasRaaVsPt", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    
    histGridRaaVsPt = ROOT.TH2D("histGridRaaVsPt",";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}",100, 0, 8, 100, 0, 1.4)
    histGridRaaVsPt.Draw()
    
    # --- Unity line ---
    lineUnityVsPt = ROOT.TLine(0., 1., 12., 1.)
    lineUnityVsPt.SetLineColor(ROOT.kGray+1)
    lineUnityVsPt.SetLineWidth(2)
    lineUnityVsPt.SetLineStyle(ROOT.kDashed)
    lineUnityVsPt.Draw("SAME")
    
    # --- Draw graphs ---
    graStatRaaVsPt.Draw("EP SAME")
    graSystRaaVsPt.Draw("E2P SAME")

    legendRaaVsPtVsRun2 = ROOT.TLegend(0.20, 0.73, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsPtVsRun2)
    legendRaaVsPtVsRun2.SetTextSize(0.045)
    legendRaaVsPtVsRun2.AddEntry(graStatRaaVsPt, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus90%", "FP")
    legendRaaVsPtVsRun2.Draw("SAME")

    latexTitle.DrawLatex(0.22, 0.87, "ALICE Work In Progress, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    
    canvasRaaVsPt.Update()
    canvasRaaVsPt.SaveAs("plots/Jpsi_RAA_vs_Pt.pdf")


    # Compare to Run 2 results
    print("***** Extract RAA from Run 2 *****")

    # --- Centrality 0–20% ---
    filePathJpsi_020 = "HEP_Data/jpsi_raa_PbPb_5TeV_0_20_centr.yaml"
    ptCenters_020, ptWidths_020, jpsiRaaVsPt_020, jpsiStatRaaVsPt_020, jpsiSystRaaVsPt_020 = ExtractFromYaml(filePathJpsi_020)
    jpsiPtSystWidths_020 = np.repeat(0.2, len(ptCenters_020))
    
    graStatRaaVsPt_020 = ROOT.TGraphErrors(len(ptCenters_020), np.array(ptCenters_020), np.array(jpsiRaaVsPt_020),
                                    np.array(ptWidths_020), np.array(jpsiStatRaaVsPt_020))
    graSystRaaVsPt_020 = ROOT.TGraphErrors(len(ptCenters_020), np.array(ptCenters_020), np.array(jpsiRaaVsPt_020),
                                    jpsiPtSystWidths_020, np.array(jpsiSystRaaVsPt_020))
    SetGraStat(graStatRaaVsPt_020, 20, ROOT.kRed+1)
    SetGraSyst(graSystRaaVsPt_020, 20, ROOT.kRed+1)
    
    # --- Centrality 20–40% ---
    filePathJpsi_2040 = "HEP_Data/jpsi_raa_PbPb_5TeV_20_40_centr.yaml"
    ptCenters_2040, ptWidths_2040, jpsiRaaVsPt_2040, jpsiStatRaaVsPt_2040, jpsiSystRaaVsPt_2040 = ExtractFromYaml(filePathJpsi_2040)
    jpsiPtSystWidths_2040 = np.repeat(0.2, len(ptCenters_2040))
    
    graStatRaaVsPt_2040 = ROOT.TGraphErrors(len(ptCenters_2040), np.array(ptCenters_2040), np.array(jpsiRaaVsPt_2040),
                                     np.array(ptWidths_2040), np.array(jpsiStatRaaVsPt_2040))
    graSystRaaVsPt_2040 = ROOT.TGraphErrors(len(ptCenters_2040), np.array(ptCenters_2040), np.array(jpsiRaaVsPt_2040),
                                     jpsiPtSystWidths_2040, np.array(jpsiSystRaaVsPt_2040))
    SetGraStat(graStatRaaVsPt_2040, 21, ROOT.kBlue)
    SetGraSyst(graSystRaaVsPt_2040, 21, ROOT.kBlue)
    
    # --- Centrality 40–60% ---
    filePathJpsi_4060 = "HEP_Data/jpsi_raa_PbPb_5TeV_40_90_centr.yaml"
    ptCenters_4060, ptWidths_4060, jpsiRaaVsPt_4060, jpsiStatRaaVsPt_4060, jpsiSystRaaVsPt_4060 = ExtractFromYaml(filePathJpsi_4060)
    jpsiPtSystWidths_4060 = np.repeat(0.2, len(ptCenters_4060))
    
    graStatRaaVsPt_4060 = ROOT.TGraphErrors(len(ptCenters_4060), np.array(ptCenters_4060), np.array(jpsiRaaVsPt_4060),
                                     np.array(ptWidths_4060), np.array(jpsiStatRaaVsPt_4060))
    graSystRaaVsPt_4060 = ROOT.TGraphErrors(len(ptCenters_4060), np.array(ptCenters_4060), np.array(jpsiRaaVsPt_4060),
                                     jpsiPtSystWidths_4060, np.array(jpsiSystRaaVsPt_4060))
    SetGraStat(graStatRaaVsPt_4060, 22, ROOT.kBlack)
    SetGraSyst(graSystRaaVsPt_4060, 22, ROOT.kBlack)
    
    # --- Canvas setup ---
    canvasRaaVsPt_Comparison1 = ROOT.TCanvas("canvasRaaVsPt_Comparison1", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    
    histGridRaaVsPt_Comparison1 = ROOT.TH2D("histGridRaaV_Comparison1",
                                ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}",
                                100, 0, 12, 100, 0, 1.7)
    histGridRaaVsPt_Comparison1.Draw()
    
    # --- Unity line ---
    lineUnityVsPt = ROOT.TLine(0., 1., 12., 1.)
    lineUnityVsPt.SetLineColor(ROOT.kGray+1)
    lineUnityVsPt.SetLineWidth(2)
    lineUnityVsPt.SetLineStyle(ROOT.kDashed)
    lineUnityVsPt.Draw("SAME")
    
    # --- Draw graphs ---
    graSystRaaVsPt_020.Draw("2 SAME")
    graStatRaaVsPt_020.Draw("P SAME")
    graSystRaaVsPt_2040.Draw("2 SAME")
    graStatRaaVsPt_2040.Draw("P SAME")
    graSystRaaVsPt_4060.Draw("2 SAME")
    graStatRaaVsPt_4060.Draw("P SAME")
    graStatRaaVsPt.Draw("EP SAME")
    graSystRaaVsPt.Draw("E2P SAME")

    legendRaaVsPtVsRun2_Comparison1 = ROOT.TLegend(0.20, 0.66, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsPtVsRun2_Comparison1)
    legendRaaVsPtVsRun2_Comparison1.SetTextSize(0.045)
    legendRaaVsPtVsRun2_Comparison1.AddEntry(graStatRaaVsPt, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus90%", "FP")
    legendRaaVsPtVsRun2_Comparison1.AddEntry(graStatRaaVsPt_020, "#sqrt{#it{s}} = 5.02 TeV, Pb-Pb, 0#minus20%", "FP")
    legendRaaVsPtVsRun2_Comparison1.AddEntry(graStatRaaVsPt_2040, "#sqrt{#it{s}} = 5.02 TeV, Pb-Pb, 20#minus60%", "FP")
    legendRaaVsPtVsRun2_Comparison1.AddEntry(graStatRaaVsPt_4060, "#sqrt{#it{s}} = 5.02 TeV, Pb-Pb, 40#minus90%", "FP")
    legendRaaVsPtVsRun2_Comparison1.Draw("SAME")

    latexTitle.DrawLatex(0.16, 0.90, "ALICE Work In Progress, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    canvasRaaVsPt_Comparison1.Update()
    canvasRaaVsPt_Comparison1.SaveAs("plots/Jpsi_RAA_vs_Pt__Comparison_manycentralities.pdf")


    # --- Centrality 0–90% ---
    filePathJpsi_090 = "HEP_Data/jpsi_raa_PbPb_5TeV_0_90_centr.yaml"
    ptCentersRun2_090, ptWidthsRun2_090, jpsiRaaVsPtRun2_090, jpsiStatRaaVsPtRun2_090, jpsiSystRaaVsPtRun2_090 = ExtractFromYaml(filePathJpsi_090)
    jpsiPtSystWidthsRun2_090 = np.repeat(0.2, len(ptCentersRun2_090))

    graStatRaaVsPtRun2_090 = ROOT.TGraphErrors(len(ptCentersRun2_090), np.array(ptCentersRun2_090), np.array(jpsiRaaVsPtRun2_090), np.array(ptWidthsRun2_090), np.array(jpsiStatRaaVsPtRun2_090))
    graSystRaaVsPtRun2_090 = ROOT.TGraphErrors(len(ptCentersRun2_090), np.array(ptCentersRun2_090), np.array(jpsiRaaVsPtRun2_090), jpsiPtSystWidthsRun2_090, np.array(jpsiSystRaaVsPtRun2_090))

    SetGraStat(graStatRaaVsPtRun2_090, 20, ROOT.kGray+2)
    SetGraSyst(graSystRaaVsPtRun2_090, 20, ROOT.kGray+2)

    lineUnityVsPt = ROOT.TLine(0., 1., 12., 1.)
    lineUnityVsPt.SetLineColor(ROOT.kGray+1)
    lineUnityVsPt.SetLineWidth(2)
    lineUnityVsPt.SetLineStyle(ROOT.kDashed)

    canvasRaaVsPtVsRun2_Comparison2 = ROOT.TCanvas("canvasRaaVsPtVsRun2_Comparison2", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    histGridRaaVsPt_Comparison2 = ROOT.TH2D("histGridRaaVsPt_Comaprison2", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}", 100, 0, 12, 100, 0, 1.7)
    histGridRaaVsPt_Comparison2.Draw()
    lineUnityVsPt.Draw("SAME")
    graStatRaaVsPtRun2_090.Draw("EP SAME")
    graSystRaaVsPtRun2_090.Draw("E2P SAME")
    graStatRaaVsPt.Draw("EP SAME")
    graSystRaaVsPt.Draw("E2P SAME")

    legendRaaVsPtVsRun2_Comparison2 = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsPtVsRun2_Comparison2)
    legendRaaVsPtVsRun2_Comparison2.SetTextSize(0.045)
    legendRaaVsPtVsRun2_Comparison2.AddEntry(graStatRaaVsPt, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus90%", "FP")
    legendRaaVsPtVsRun2_Comparison2.AddEntry(graStatRaaVsPtRun2_090, "#sqrt{#it{s}} = 5.02 TeV, Pb-Pb, 0#minus90%", "FP")
    legendRaaVsPtVsRun2_Comparison2.Draw("SAME")

    latexTitle.DrawLatex(0.22, 0.87, "ALICE Work In Progress, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")

    canvasRaaVsPtVsRun2_Comparison2.Update()
    canvasRaaVsPtVsRun2_Comparison2.SaveAs("plots/Jpsi_RAA_vs_Pt__Comparison_0_90_centr.pdf")

    # ***** Comparison with pO *****
    x_vect0 = np.array([0.5, 1.5, 2.5, 3.5, 5.0, 8.0], dtype='float64')
    y_vect1 = np.array([0.7010977121877087, 0.8096814426116999, 0.821249497674029, 0.9733794606858844, 0.8910243018805722, 1.250866084626591], dtype='float64')
    ex_vect2 = np.array([0, 0, 0, 0, 0, 0], dtype='float64')
    ey_vect3 = np.array([0.05321832374435311, 0.05126399014605561, 0.05428298032696315, 0.0715452447917095, 0.06328270207447632, 0.1384129436083637], dtype='float64')

    graStatRaaVsPtComparison_pO = ROOT.TGraphErrors(len(x_vect0),
                                                 x_vect0,
                                                 y_vect1,
                                                 ex_vect2,
                                                 ey_vect3)

    SetGraStat(graStatRaaVsPtComparison_pO, 21, ROOT.kRed+1)

    canvasRaaVsPtComparison_pO = ROOT.TCanvas("canvasRaaVsPtComparison_pO", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)

    histGridRaaVsPtComparison_pO = ROOT.TH2D("histGridRaaVsPtComparison_pO", 
                                          ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}", 
                                          100, 0, 12, 100, 0, 2)
    histGridRaaVsPtComparison_pO.Draw()

    # --- Unity line ---
    lineUnityComparison_pO = ROOT.TLine(0., 1., 12., 1.)
    lineUnityComparison_pO.SetLineColor(ROOT.kGray+1)
    lineUnityComparison_pO.SetLineWidth(2)
    lineUnityComparison_pO.SetLineStyle(ROOT.kDashed)
    lineUnityComparison_pO.Draw("SAME")

    graStatRaaVsPt.Draw("EP SAME")
    graSystRaaVsPt.Draw("E2P SAME")

    graStatRaaVsPtComparison_pO.Draw("EP SAME")

    legendRaaVsPtComparison_pO = ROOT.TLegend(0.20, 0.70, 0.45, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsPtComparison_pO)
    legendRaaVsPtComparison_pO.SetTextSize(0.045)
    legendRaaVsPtComparison_pO.AddEntry(graStatRaaVsPt, "#sqrt{#it{s}} = 5.36 TeV, OO,  2.5 < #it{y} < 4", "FP")
    legendRaaVsPtComparison_pO.AddEntry(graStatRaaVsPtComparison_pO, "#sqrt{#it{s}} = 9.62 TeV, pO,  2.15 < #it{y} < 3.65", "FP")
    legendRaaVsPtComparison_pO.Draw("SAME")

    latexTitle.DrawLatex(0.22, 0.87, "ALICE Work In Progress, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}")

    canvasRaaVsPtComparison_pO.Update()
    canvasRaaVsPtComparison_pO.SaveAs("plots/Jpsi_RAA_vs_Pt_Comparison_with_pO.pdf")

    # ***** Comparison pO and Run 2 results pPb *****
    filePathJpsi_pPb = "HEP_Data/jpsi_raa_PbPb_5TeV_pPb.yaml"
    ptCentersRun2_pPb, ptWidthsRun2_pPb, jpsiRaaVsPtRun2_pPb, jpsiStatRaaVsPtRun2_pPb, jpsiSystRaaVsPtRun2_pPb = ExtractFromYaml(filePathJpsi_pPb)
    jpsiPtSystWidthsRun2_pPb = np.repeat(0.2, len(ptCentersRun2_pPb))

    graStatRaaVsPtRun2_pPb = ROOT.TGraphErrors(len(ptCentersRun2_pPb), np.array(ptCentersRun2_pPb), np.array(jpsiRaaVsPtRun2_pPb), np.array(ptWidthsRun2_pPb), np.array(jpsiStatRaaVsPtRun2_pPb))
    graSystRaaVsPtRun2_pPb = ROOT.TGraphErrors(len(ptCentersRun2_pPb), np.array(ptCentersRun2_pPb), np.array(jpsiRaaVsPtRun2_pPb), jpsiPtSystWidthsRun2_pPb, np.array(jpsiSystRaaVsPtRun2_pPb))

    SetGraStat(graStatRaaVsPtRun2_pPb, 20, ROOT.kGray+2)
    SetGraSyst(graSystRaaVsPtRun2_pPb, 20, ROOT.kGray+2)

    lineUnityVsPt_pPb = ROOT.TLine(0., 1., 12., 1.)
    lineUnityVsPt_pPb.SetLineColor(ROOT.kGray+1)
    lineUnityVsPt_pPb.SetLineWidth(2)
    lineUnityVsPt_pPb.SetLineStyle(ROOT.kDashed)

    canvasRaaVsPtVsRun2_pPb = ROOT.TCanvas("canvasRaaVsPtVsRun2_pPb", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    histGridRaaVsPt_pPb = ROOT.TH2D("histGridRaaVsPt_pPb", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}", 100, 0, 12, 100, 0, 2.0)
    histGridRaaVsPt_pPb.Draw()
    lineUnityVsPt_pPb.Draw("SAME")
    graStatRaaVsPtRun2_pPb.Draw("EP SAME")
    graSystRaaVsPtRun2_pPb.Draw("E2P SAME")
    #graStatRaaVsPt.Draw("EP SAME")
    #graSystRaaVsPt.Draw("E2P SAME")

    legendRaaVsPtVsRun2_pPb = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsPtVsRun2_pPb)
    legendRaaVsPtVsRun2_pPb.SetTextSize(0.045)
    legendRaaVsPtVsRun2_pPb.AddEntry(graStatRaaVsPtComparison_pO, "#sqrt{#it{s}} = 9.62 TeV, pO, 2.15 < #it{y} < 3.65", "FP")
    legendRaaVsPtVsRun2_pPb.AddEntry(graStatRaaVsPtRun2_pPb, "#sqrt{#it{s}} = 8.016 TeV, p-Pb, 2.03  < #it{y} < 3.53", "FP")
    legendRaaVsPtVsRun2_pPb.Draw("SAME")
    graStatRaaVsPtComparison_pO.Draw("EP SAME")

    latexTitle.DrawLatex(0.22, 0.87, "ALICE Work In Progress, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")

    canvasRaaVsPtVsRun2_pPb.Update()
    canvasRaaVsPtVsRun2_pPb.SaveAs("plots/Jpsi_Comparison_pO_pPbRun2.pdf")

    input()

    print(f"***** Write results in {config['outputs']['fOut']} *****")
    fOut = ROOT.TFile(config["outputs"]["fOut"], "RECREATE")
    histStatRawYieldVsPt.Write()
    histSystRawYieldVsPt.Write()
    histStatAxeVsPt.Write()
    histStatXsecVsPt.Write()
    histSystXsecVsPt.Write()

    fOut.mkdir("systematics")
    fOut.cd("systematics")
    histSystRelRawYieldVsPt.Write("syst_raw_yield_vs_pt")
    fOut.Close()
    """


if __name__ == '__main__':
    main()