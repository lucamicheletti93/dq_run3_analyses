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

    exit()

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
    deltaY = inputCfg["crossSec"]["deltaY"]
    histCrossSec = histRawYield.Clone("histCrossSec")
    histCrossSec.Divide(histAxe)
    histCrossSec.Scale(1. / (branchingRatio * limiInt * deltaY))
    histCrossSec.SetLineColor(ROOT.kRed+1)
    histCrossSec.SetLineWidth(3)

    ptCentrRun2 = inputCfg["crossSec"]["ptCentrRun2"]
    ptWidthRun2 = inputCfg["crossSec"]["ptWidthRun2"]
    crossSecRun2 = inputCfg["crossSec"]["crossSecRun2"]
    statCrossSecRun2 = inputCfg["crossSec"]["statCrossSecRun2"]

    graCrossSecRun2 = ROOT.TGraphErrors(len(ptCentrRun2), np.array(ptCentrRun2), np.array(crossSecRun2), np.array(ptWidthRun2), np.array(statCrossSecRun2))
    graCrossSecRun2.SetMarkerStyle(20)
    graCrossSecRun2.SetMarkerColor(ROOT.kBlack)
    graCrossSecRun2.SetLineColor(ROOT.kBlack)

    canvasCrossSec = ROOT.TCanvas("canvasCrossSec", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histCrossSec.GetYaxis().SetTitle("d^{2}#sigma/d#it{p}_{T}dy (nb^{-1} / GeV/#it{c})")
    histCrossSec.SetTitle("")
    histCrossSec.SetStats(False)
    histCrossSec.Draw("H")
    graCrossSecRun2.Draw("EP SAME")

    legend = ROOT.TLegend(0.65, 0.67, 0.82, 0.92, " ", "brNDC")
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

def simple_ratio(inputCfg):
    """
    function to compute the simple ratio of charmonia states
    """
    LoadStyle()

    # Acceptance-times-efficiency
    fInYields = ROOT.TFile(inputCfg["simpleRatio"]["fInNameYields"], "READ")
    histRawYield1 = fInYields.Get(inputCfg["simpleRatio"]["histNameYields"][0])
    histRawYield2 = fInYields.Get(inputCfg["simpleRatio"]["histNameYields"][1])

    fInAxe1 = ROOT.TFile(inputCfg["simpleRatio"]["fInNameAxe1"], "READ")
    histAxe1 = fInAxe1.Get("histAxe")

    fInAxe2 = ROOT.TFile(inputCfg["simpleRatio"]["fInNameAxe2"], "READ")
    histAxe2 = fInAxe2.Get("histAxe")

    histRawYieldRatios = histRawYield1.Clone("histRawYieldRatios")
    histRawYieldRatios.Divide(histRawYield2)

    histAxeRatios = histAxe1.Clone("histAxeRatios")
    histAxeRatios.Divide(histAxe2)

    histCrossSecRatios = histRawYieldRatios.Clone("histCrossSecRatios")
    histCrossSecRatios.Divide(histAxeRatios)
    histCrossSecRatios.Scale(inputCfg["simpleRatio"]["branchingRatios"][1] / inputCfg["simpleRatio"]["branchingRatios"][0])

    canvasCrossSecRatio = ROOT.TCanvas("canvasCrossSecRatio", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histCrossSecRatios.GetYaxis().SetTitle("#sigma^{#psi(2S)} / #sigma^{J/#psi}")
    histCrossSecRatios.SetTitle("")
    histCrossSecRatios.SetStats(False)
    histCrossSecRatios.Draw("H")
    #graCrossSecRun2.Draw("EP SAME")

    legend = ROOT.TLegend(0.65, 0.67, 0.82, 0.92, " ", "brNDC")
    SetLegend(legend)
    legend.SetTextSize(0.04)
    legend.AddEntry(histCrossSecRatios, inputCfg["simpleRatio"]["label"], "L")
    #legend.AddEntry(graCrossSecRun2, inputCfg["crossSec"]["labelRun2"], "PL")
    legend.Draw()

    canvasCrossSecRatio.Update()

    input()

def luminosity(inputCfg):
    """
    function to compute luminosity for triggered data -> outout of check_triggers.C
    """
    LoadStyle()

    sigmaTVX = inputCfg["luminosity"]["sigmaTVX"] # pb
    path = inputCfg["luminosity"]["path"]
    periods = inputCfg["luminosity"]["periods"]
    trigger = inputCfg["luminosity"]["trigger"]

    fInNames = [f'{path}/{period}_{trigger}_trigger_summary.root' for period in periods]

    lumiInt = 0
    for iPeriod, fInName in enumerate(fInNames):
        counterTVX = 0
        counterScalTrig = 0
        counterSelTrig = 0
        fIn = ROOT.TFile(fInName, "READ")
        histCounterTVX = fIn.Get("histCounterTVX")
        histCounterScalTrig = fIn.Get("histCounterScalTrig")
        histCounterSelTrig = fIn.Get("histCounterSelTrig")

        for iRun in range(0, int(histCounterTVX.GetEntries())):
            counterTVX += histCounterTVX.GetBinContent(iRun+1)
            counterScalTrig += histCounterScalTrig.GetBinContent(iRun+1)
            counterSelTrig += histCounterSelTrig.GetBinContent(iRun+1)

        lumi = (counterTVX * (counterScalTrig / counterSelTrig)) / sigmaTVX # pb
        print(f'{periods[iPeriod]} = {lumi} pb-1')
        lumiInt += lumi
    print(f'--> luminosity LHC24 = {lumiInt}')






### ### ###
def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--cross_section", help="Compute the cross section for skimmed events", action="store_true")
    parser.add_argument("--simple_ratio", help="Compute the simple ratio for skimmed events", action="store_true")
    parser.add_argument("--luminosity", help="Compute the luminosity for skimmed events", action="store_true")
    args = parser.parse_args()

    print('Loading task configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading task configuration: Done!')

    if args.cross_section:
        cross_section(inputCfg)
    
    if args.simple_ratio:
        simple_ratio(inputCfg)
    
    if args.luminosity:
        luminosity(inputCfg)

main()