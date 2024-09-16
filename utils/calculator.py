import matplotlib.pyplot as plt
import array as arr
import numpy as np
import os
import sys
import argparse
import yaml
import ROOT
from os import path
from plot_library import LoadStyle, SetLegend

def cross_section(inputCfg):
    """
    function to compute the cross section of charmonia
    """
    LoadStyle()

    fInTriggers =  ROOT.TFile(inputCfg["lumi"]["fInName"], "READ")
    dir = fInTriggers.Get("table-maker/Zorro")
    runs = [key.GetName() for key in dir.GetListOfKeys()]
    histInspectedTVX = ROOT.TH1D("histInspectedTVX", "", len(runs), 0, len(runs))
    histSelections = ROOT.TH1D("histSelections", "", len(runs), 0, len(runs))
    histAnalysedTriggersOfInterest = ROOT.TH1D("histAnalysedTriggersOfInterest", "", len(runs), 0, len(runs))

    inspectedTVX = 0
    selections = 0
    analysedTriggersOfInterest = 0

    for iRun, run in enumerate(runs):
        # Inspected TVX triggers
        histTmp1 = fInTriggers.Get(f'table-maker/Zorro/{run}/InspectedTVX')
        inspectedTVX += histTmp1.GetBinContent(1)
        histInspectedTVX.SetBinContent(iRun+1, histTmp1.GetBinContent(1))
        histInspectedTVX.GetXaxis().SetBinLabel(iRun+1, f'{run}')

        # Originally selected TOI
        histTmp2 = fInTriggers.Get(f'table-maker/Zorro/{run}/Selections')
        selections += histTmp2.GetBinContent(histTmp2.GetXaxis().FindBin("fDiMuon"))
        histSelections.SetBinContent(iRun+1, histTmp2.GetBinContent(histTmp2.GetXaxis().FindBin("fDiMuon")))
        histSelections.GetXaxis().SetBinLabel(iRun+1, f'{run}')

        # Analyzed selected TOI
        histTmp3 = fInTriggers.Get(f'table-maker/Zorro/{run}/AnalysedTriggersOfInterest')
        analysedTriggersOfInterest += histTmp3.GetBinContent(histTmp3.GetXaxis().FindBin("fDiMuon"))
        histAnalysedTriggersOfInterest.SetBinContent(iRun+1, histTmp3.GetBinContent(histTmp3.GetXaxis().FindBin("fDiMuon")))
        histAnalysedTriggersOfInterest.GetXaxis().SetBinLabel(iRun+1, f'{run}')

    limiInt = (inspectedTVX * (analysedTriggersOfInterest / selections)) / 59.4e6 # mb
    print("Integrated Luminosity = ", limiInt, "nb-1")
    histLumi = ROOT.TH1D("histLumi", "", 1, 0, 1)
    histLumi.GetYaxis().SetTitle("Lumi. (nb^{-1})")
    histLumi.SetBinContent(1, limiInt)

    # Raw yields
    fInRawYields = ROOT.TFile(inputCfg["yields"]["fInName"], "READ")
    histRawYield = fInRawYields.Get(inputCfg["yields"]["histName"])
    histRawYield.SetLineColor(ROOT.kRed+1)

    canvasRawYelds = ROOT.TCanvas("canvasRawYelds", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histRawYield.GetYaxis().SetTitle("d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}")
    histRawYield.SetTitle("")
    histRawYield.SetStats(False)
    histRawYield.Draw("H")
    canvasRawYelds.Update()

    # Acceptance-times-efficiency
    fInAxe = ROOT.TFile(inputCfg["axe"]["fInName"], "READ")
    listGen = fInAxe.Get("analysis-same-event-pairing/output")

    listGenCand = listGen.FindObject(inputCfg["axe"]["histGenName"])
    histGenPtRap = listGenCand.FindObject("Pt_Rapidity")
    histGenPtRap.GetXaxis().SetRangeUser(inputCfg["axe"]["kineMinCuts"][0], inputCfg["axe"]["kineMaxCuts"][0])
    histGenPtRap.GetYaxis().SetRangeUser(inputCfg["axe"]["kineMinCuts"][1], inputCfg["axe"]["kineMaxCuts"][1])
    histGen = histGenPtRap.ProjectionX("histGen")

    listRec = fInAxe.Get("analysis-same-event-pairing/output")
    listRecCand = listRec.FindObject(inputCfg["axe"]["histRecName"])
    histRecMassPtRap = listRecCand.FindObject("Mass_Pt_Rapidity")
    histRecMassPtRap.SetName("histRecMassPtRap")
    histRecMassPtRap.GetAxis(1).SetRangeUser(inputCfg["axe"]["kineMinCuts"][0], inputCfg["axe"]["kineMaxCuts"][0])
    histRecMassPtRap.GetAxis(2).SetRangeUser(inputCfg["axe"]["kineMinCuts"][1], inputCfg["axe"]["kineMaxCuts"][1])
    histRec = histRecMassPtRap.Projection(1, "histRec")

    histGen.Rebin(inputCfg["axe"]["rebin"])
    histRec.Rebin(inputCfg["axe"]["rebin"])

    histAxe = histRec.Clone("histAxe")
    histAxe.Divide(histGen)
    histAxe.SetLineColor(ROOT.kRed+1)

    canvasAxe = ROOT.TCanvas("canvasAxe", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histAxe.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
    histAxe.GetYaxis().SetTitle("A#times#varepsilon")
    histAxe.SetTitle("")
    histAxe.SetStats(False)
    histAxe.Draw("H")
    canvasAxe.Update()

    # Cross sections
    branchingRatio = inputCfg["crossSec"]["branchingRatio"]
    histCrossSec = histRawYield.Clone("histCrossSec")
    histCrossSec.Divide(histAxe)
    histCrossSec.Scale(1. / (branchingRatio * limiInt))
    histCrossSec.SetLineColor(ROOT.kRed+1)

    ptCentrRun2 = inputCfg["crossSec"]["ptCentrRun2"]
    ptWidthRun2 = inputCfg["crossSec"]["ptWidthRun2"]
    crossSecRun2 = inputCfg["crossSec"]["crossSecRun2"]
    statCrossSecRun2 = inputCfg["crossSec"]["statCrossSecRun2"]

    graCrossSecRun2 = ROOT.TGraphErrors(len(ptCentrRun2), np.array(ptCentrRun2), np.array(crossSecRun2), np.array(ptWidthRun2), np.array(statCrossSecRun2))
    graCrossSecRun2.SetMarkerStyle(20)
    graCrossSecRun2.SetMarkerColor(ROOT.kRed+1)
    graCrossSecRun2.SetLineColor(ROOT.kRed+1)

    canvasCrossSec = ROOT.TCanvas("canvasCrossSec", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histCrossSec.GetYaxis().SetTitle("d#sigma/d#it{p}_{T} (GeV/#it{c})^{-1}")
    histCrossSec.SetTitle("")
    histCrossSec.SetStats(False)
    histCrossSec.Draw("H")
    graCrossSecRun2.Draw("EP SAME")

    legend = ROOT.TLegend(0.65, 0.57, 0.82, 0.92, " ", "brNDC")
    SetLegend(legend)
    legend.SetTextSize(0.04)
    legend.AddEntry(histCrossSec, inputCfg["crossSec"]["label"], "L")
    legend.AddEntry(graCrossSecRun2, inputCfg["crossSec"]["labelRun2"], "PL")
    legend.Draw()

    canvasCrossSec.Update()

    fOut = ROOT.TFile(f'{inputCfg["output"]["fOutName"]}.root', "RECREATE")
    histInspectedTVX.Write("histInspectedTVX")
    histSelections.Write("histSelections")
    histAnalysedTriggersOfInterest.Write("histAnalysedTriggersOfInterest")
    histLumi.Write("histLumi")
    histRawYield.Write("histRawYield")
    histAxe.Write("histAxe")
    histCrossSec.Write("histCrossSec")
    fOut.Close()

    input()



### ### ###
def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--cross_section", help="Compute luminosity for skimmed events", action="store_true")
    args = parser.parse_args()

    print('Loading task configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading task configuration: Done!')

    if args.cross_section:
        cross_section(inputCfg)

main()