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

def qa(config):
    """
    function to compute the cross section
    """
    LoadStyle()
    ROOT.gStyle.SetOptStat(False)
    BrJpsiToMuMu =  0.05961
    errBrJpsiToMuMu = 0.033e-2
    BrPsi2sToMuMu =  0.008
    errBrPsi2sToMuMu = 0.6e-3

    recoTags = config["inputs"]["recoTags"]
    colors = {
        "pO": ROOT.kAzure,
        "OO": ROOT.kRed,
        "NeNe": ROOT.kGreen+1,
    }
    styles = [24, 20]

    lineJpsiMass = ROOT.TLine(0, 3096, 10, 3096)
    lineJpsiMass.SetLineStyle(ROOT.kDashed)
    lineJpsiMass.SetLineColor(ROOT.kGray+2)

    latexTitle = ROOT.TLatex()
    latexTitle.SetNDC()
    latexTitle.SetTextSize(0.05)
    latexTitle.SetTextFont(42)

    nCollspO = []
    histCentralitypO = []
    for iProd, fInpO in enumerate(config["inputs"]["fInpO"]):
        fInpO = ROOT.TFile(fInpO, "READ")
        histEvents = fInpO.Get("event-selection-task/hColCounterTVX")
        if not histEvents:
            histEvents = fInpO.Get("eventselection-run3/eventselection/hColCounterTVX")
        nCollspO.append(histEvents.GetEntries())

        hlistCentrality = fInpO.Get("table-maker/output")
        histCentralitypO.append((hlistCentrality.FindObject("Event_AfterCuts")).FindObject("CentFT0C"))

    histStatRawYieldpO = []
    histSystRawYieldpO = []
    histStatMeanpO = []
    histStatWidthpO = []
    for iProd, pathFitResults in enumerate(config["inputs"]["pathFitResultspO"]):
        dfRawYield = pd.read_csv(f'{pathFitResults}/systematic_sig_Jpsi.txt', sep=' ')
        ptMin = dfRawYield["x_min"]
        ptMax = dfRawYield["x_max"]
        ptEdges = np.append(ptMin.to_numpy(), ptMax.to_numpy()[len(ptMax)-1])

        jpsiRawYields = dfRawYield["val"]
        jpsiStatRawYields = dfRawYield["stat"]
        jpsiSystRawYields = dfRawYield["syst"]

        dfMean = pd.read_csv(f'{pathFitResults}/systematic_mean_Jpsi.txt', sep=' ')
        jpsiMeans = dfMean["val"]
        jpsiStatMeans = dfMean["stat"]

        dfWidth = pd.read_csv(f'{pathFitResults}/systematic_width_Jpsi.txt', sep=' ')
        jpsiWidths = dfWidth["val"]
        jpsiStatWidths = dfWidth["stat"]

        histStatRawYieldpO.append(ROOT.TH1D(f'histStatRawYield_pO_{recoTags[iProd]}', ";#it{p}_{T} (GeV/#it{c});1/N_{TVX} #times dN/d#it{p}_{T}", len(ptEdges)-1, ptEdges))
        histSystRawYieldpO.append(ROOT.TH1D(f'histSystRawYield_pO_{recoTags[iProd]}', ";#it{p}_{T} (GeV/#it{c});1/N_{TVX} #times dN/d#it{p}_{T}", len(ptEdges)-1, ptEdges))
        histStatMeanpO.append(ROOT.TH1D(f'histStatMean_pO_{recoTags[iProd]}', ";#it{p}_{T} (GeV/#it{c});#mu_{J/#psi} (MeV/c^{2})", len(ptEdges)-1, ptEdges))
        histStatWidthpO.append(ROOT.TH1D(f'histStatWidth_pO_{recoTags[iProd]}', ";#it{p}_{T} (GeV/#it{c});#sigma_{J/#psi} (MeV/c^{2})", len(ptEdges)-1, ptEdges))
        for iPt, (jpsiRawYield,jpsiStatRawYield,jpsiSystRawYield) in enumerate(zip(jpsiRawYields,jpsiStatRawYields,jpsiSystRawYields)):
            histStatRawYieldpO[iProd].SetBinContent(iPt+1, jpsiRawYield)
            histStatRawYieldpO[iProd].SetBinError(iPt+1, jpsiStatRawYield)
            histSystRawYieldpO[iProd].SetBinContent(iPt+1, jpsiRawYield)
            histSystRawYieldpO[iProd].SetBinError(iPt+1, jpsiSystRawYield)

            histStatMeanpO[iProd].SetBinContent(iPt+1, jpsiMeans[iPt])
            histStatMeanpO[iProd].SetBinError(iPt+1, jpsiStatMeans[iPt])

            histStatWidthpO[iProd].SetBinContent(iPt+1, jpsiWidths[iPt])
            histStatWidthpO[iProd].SetBinError(iPt+1, jpsiStatWidths[iPt])

    for iProd in range(0, len(recoTags)):
        SetHistStat(histStatRawYieldpO[iProd], styles[iProd], colors["pO"]+iProd)
        SetHistSyst(histSystRawYieldpO[iProd], styles[iProd], colors["pO"]+iProd)
        SetHistStat(histStatMeanpO[iProd], styles[iProd], colors["pO"]+iProd)
        SetHistStat(histStatWidthpO[iProd], styles[iProd], colors["pO"]+iProd)
        histStatRawYieldpO[iProd].Scale(1. / nCollspO[iProd], "WIDTH")
        histSystRawYieldpO[iProd].Scale(1. / nCollspO[iProd], "WIDTH")
        histStatMeanpO[iProd].Scale(1000.)
        histStatWidthpO[iProd].Scale(1000.)

        SetHistStat(histCentralitypO[iProd], styles[iProd], colors["pO"]+iProd)
        histCentralitypO[iProd].Scale(1. / histCentralitypO[iProd].Integral())
        histCentralitypO[iProd].SetTitle("")
        histCentralitypO[iProd].SetName(f'histCentrality_pO_{recoTags[iProd]}')
        histCentralitypO[iProd].SetTitleSize(0.045,"X")
        histCentralitypO[iProd].SetTitleSize(0.045,"Y")
        histCentralitypO[iProd].SetLabelSize(0.045,"X")
        histCentralitypO[iProd].SetLabelSize(0.045,"Y")
        histCentralitypO[iProd].SetTitleOffset(1.2,"X")
        histCentralitypO[iProd].SetTitleOffset(1.35,"Y")

    canvasRawYieldpO = ROOT.TCanvas("canvasRawYieldpO", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histStatRawYieldpO[0].GetYaxis().SetRangeUser(0.01 * histStatRawYieldpO[0].GetMaximum(), 2 * histStatRawYieldpO[0].GetMaximum())
    histStatRawYieldpO[0].Draw("EP")
    for iProd in range(0, len(recoTags)):
        if iProd == 0:
            histSystRawYieldpO[iProd].Draw("E2P SAME")
        else:
            histStatRawYieldpO[iProd].Draw("EP SAME")
            histSystRawYieldpO[iProd].Draw("E2P SAME")
    canvasRawYieldpO.Update()

    canvasMeanpO = ROOT.TCanvas("canvasMeanpO", "", 800, 600)
    histStatMeanpO[0].GetYaxis().SetRangeUser(3050, 3150)
    for iProd in range(0, len(recoTags)):
        if iProd == 0:
            histStatMeanpO[iProd].Draw("EP")
        else:
            histStatMeanpO[iProd].Draw("EP SAME")
    lineJpsiMass.Draw()
    canvasMeanpO.Update()

    canvasWidthpO = ROOT.TCanvas("canvasWidthpO", "", 800, 600)
    histStatWidthpO[0].GetYaxis().SetRangeUser(0, 150)
    for iProd in range(0, len(recoTags)):
        if iProd == 0:
            histStatWidthpO[iProd].Draw("EP")
        else:
            histStatWidthpO[iProd].Draw("EP SAME")
    canvasWidthpO.Update()

    #########################################################
    nCollsOO = []
    histCentralityOO = []
    for iProd, fInOO in enumerate(config["inputs"]["fInOO"]):
        fInOO = ROOT.TFile(fInOO, "READ")
        histEvents = fInOO.Get("event-selection-task/hColCounterTVX")
        if not histEvents:
            histEvents = fInOO.Get("eventselection-run3/eventselection/hColCounterTVX")
        nCollsOO.append(histEvents.GetEntries())

        hlistCentrality = fInOO.Get("table-maker/output")
        histCentralityOO.append((hlistCentrality.FindObject("Event_AfterCuts")).FindObject("CentFT0C"))

    histStatRawYieldOO = []
    histSystRawYieldOO = []
    histStatMeanOO = []
    histStatWidthOO = []
    for iProd, pathFitResults in enumerate(config["inputs"]["pathFitResultsOO"]):
        dfRawYield = pd.read_csv(f'{pathFitResults}/systematic_sig_Jpsi.txt', sep=' ')
        ptMin = dfRawYield["x_min"]
        ptMax = dfRawYield["x_max"]
        ptEdges = np.append(ptMin.to_numpy(), ptMax.to_numpy()[len(ptMax)-1])

        jpsiRawYields = dfRawYield["val"]
        jpsiStatRawYields = dfRawYield["stat"]
        jpsiSystRawYields = dfRawYield["syst"]

        dfMean = pd.read_csv(f'{pathFitResults}/systematic_mean_Jpsi.txt', sep=' ')
        jpsiMeans = dfMean["val"]
        jpsiStatMeans = dfMean["stat"]

        dfWidth = pd.read_csv(f'{pathFitResults}/systematic_width_Jpsi.txt', sep=' ')
        jpsiWidths = dfWidth["val"]
        jpsiStatWidths = dfWidth["stat"]

        histStatRawYieldOO.append(ROOT.TH1D(f'histStatRawYield_OO_{recoTags[iProd]}', ";#it{p}_{T} (GeV/#it{c});1/N_{TVX} #times dN/d#it{p}_{T}", len(ptEdges)-1, ptEdges))
        histSystRawYieldOO.append(ROOT.TH1D(f'histSystRawYield_OO_{recoTags[iProd]}', ";#it{p}_{T} (GeV/#it{c});1/N_{TVX} #times dN/d#it{p}_{T}", len(ptEdges)-1, ptEdges))
        histStatMeanOO.append(ROOT.TH1D(f'histStatMean_OO_{recoTags[iProd]}', ";#it{p}_{T} (GeV/#it{c});#mu_{J/#psi} (MeV/c^{2})", len(ptEdges)-1, ptEdges))
        histStatWidthOO.append(ROOT.TH1D(f'histStatWidth_OO_{recoTags[iProd]}', ";#it{p}_{T} (GeV/#it{c});#sigma_{J/#psi} (MeV/c^{2})", len(ptEdges)-1, ptEdges))
        for iPt, (jpsiRawYield,jpsiStatRawYield,jpsiSystRawYield) in enumerate(zip(jpsiRawYields,jpsiStatRawYields,jpsiSystRawYields)):
            histStatRawYieldOO[iProd].SetBinContent(iPt+1, jpsiRawYield)
            histStatRawYieldOO[iProd].SetBinError(iPt+1, jpsiStatRawYield)
            histSystRawYieldOO[iProd].SetBinContent(iPt+1, jpsiRawYield)
            histSystRawYieldOO[iProd].SetBinError(iPt+1, jpsiSystRawYield)

            histStatMeanOO[iProd].SetBinContent(iPt+1, jpsiMeans[iPt])
            histStatMeanOO[iProd].SetBinError(iPt+1, jpsiStatMeans[iPt])

            histStatWidthOO[iProd].SetBinContent(iPt+1, jpsiWidths[iPt])
            histStatWidthOO[iProd].SetBinError(iPt+1, jpsiStatWidths[iPt])

    for iProd in range(0, len(recoTags)):
        SetHistStat(histStatRawYieldOO[iProd], styles[iProd], colors["OO"]+iProd)
        SetHistSyst(histSystRawYieldOO[iProd], styles[iProd], colors["OO"]+iProd)
        SetHistStat(histStatMeanOO[iProd], styles[iProd], colors["OO"]+iProd)
        SetHistStat(histStatWidthOO[iProd], styles[iProd], colors["OO"]+iProd)
        histStatRawYieldOO[iProd].Scale(1. / nCollsOO[iProd], "WIDTH")
        histSystRawYieldOO[iProd].Scale(1. / nCollsOO[iProd], "WIDTH")
        histStatMeanOO[iProd].Scale(1000.)
        histStatWidthOO[iProd].Scale(1000.)

        SetHistStat(histCentralityOO[iProd], styles[iProd], colors["OO"]+iProd)
        histCentralityOO[iProd].Scale(2. / histCentralityOO[iProd].Integral())
        histCentralityOO[iProd].SetTitle("")
        histCentralityOO[iProd].SetName(f'histCentrality_OO_{recoTags[iProd]}')
        histCentralityOO[iProd].SetTitleSize(0.045,"X")
        histCentralityOO[iProd].SetTitleSize(0.045,"Y")
        histCentralityOO[iProd].SetLabelSize(0.045,"X")
        histCentralityOO[iProd].SetLabelSize(0.045,"Y")
        histCentralityOO[iProd].SetTitleOffset(1.2,"X")
        histCentralityOO[iProd].SetTitleOffset(1.35,"Y")

    canvasRawYieldOO = ROOT.TCanvas("canvasRawYieldOO", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histStatRawYieldOO[0].GetYaxis().SetRangeUser(0.01 * histStatRawYieldOO[0].GetMaximum(), 2 * histStatRawYieldOO[0].GetMaximum())
    histStatRawYieldOO[0].Draw("EP")
    for iProd in range(0, len(recoTags)):
        if iProd == 0:
            histSystRawYieldOO[iProd].Draw("E2P SAME")
        else:
            histStatRawYieldOO[iProd].Draw("EP SAME")
            histSystRawYieldOO[iProd].Draw("E2P SAME")
    canvasRawYieldOO.Update()

    canvasMeanOO = ROOT.TCanvas("canvasMeanOO", "", 800, 600)
    histStatMeanOO[0].GetYaxis().SetRangeUser(3050, 3150)
    for iProd in range(0, len(recoTags)):
        if iProd == 0:
            histStatMeanOO[iProd].Draw("EP")
        else:
            histStatMeanOO[iProd].Draw("EP SAME")
    lineJpsiMass.Draw()
    canvasMeanOO.Update()

    canvasWidthOO = ROOT.TCanvas("canvasWidthOO", "", 800, 600)
    histStatWidthOO[0].GetYaxis().SetRangeUser(0, 150)
    for iProd in range(0, len(recoTags)):
        if iProd == 0:
            histStatWidthOO[iProd].Draw("EP")
        else:
            histStatWidthOO[iProd].Draw("EP SAME")
    canvasWidthOO.Update()

    #########################################################
    nCollsNeNe = []
    histCentralityNeNe = []
    for iProd, fInNeNe in enumerate(config["inputs"]["fInNeNe"]):
        fInNeNe = ROOT.TFile(fInNeNe, "READ")
        histEvents = fInNeNe.Get("event-selection-task/hColCounterTVX")
        if not histEvents:
            histEvents = fInNeNe.Get("eventselection-run3/eventselection/hColCounterTVX")
        nCollsNeNe.append(histEvents.GetEntries())

        hlistCentrality = fInNeNe.Get("table-maker/output")
        histCentralityNeNe.append((hlistCentrality.FindObject("Event_AfterCuts")).FindObject("CentFT0C"))

    histStatRawYieldNeNe = []
    histSystRawYieldNeNe = []
    histStatMeanNeNe = []
    histStatWidthNeNe = []
    for iProd, pathFitResults in enumerate(config["inputs"]["pathFitResultsNeNe"]):
        dfRawYield = pd.read_csv(f'{pathFitResults}/systematic_sig_Jpsi.txt', sep=' ')
        ptMin = dfRawYield["x_min"]
        ptMax = dfRawYield["x_max"]
        ptEdges = np.append(ptMin.to_numpy(), ptMax.to_numpy()[len(ptMax)-1])

        jpsiRawYields = dfRawYield["val"]
        jpsiStatRawYields = dfRawYield["stat"]
        jpsiSystRawYields = dfRawYield["syst"]

        dfMean = pd.read_csv(f'{pathFitResults}/systematic_mean_Jpsi.txt', sep=' ')
        jpsiMeans = dfMean["val"]
        jpsiStatMeans = dfMean["stat"]

        dfWidth = pd.read_csv(f'{pathFitResults}/systematic_width_Jpsi.txt', sep=' ')
        jpsiWidths = dfWidth["val"]
        jpsiStatWidths = dfWidth["stat"]

        histStatRawYieldNeNe.append(ROOT.TH1D(f'histStatRawYield_NeNe_{recoTags[iProd]}', ";#it{p}_{T} (GeV/#it{c});1/N_{TVX} #times dN/d#it{p}_{T}", len(ptEdges)-1, ptEdges))
        histSystRawYieldNeNe.append(ROOT.TH1D(f'histSystRawYield_NeNe_{recoTags[iProd]}', ";#it{p}_{T} (GeV/#it{c});1/N_{TVX} #times dN/d#it{p}_{T}", len(ptEdges)-1, ptEdges))
        histStatMeanNeNe.append(ROOT.TH1D(f'histStatMean_NeNe_{recoTags[iProd]}', ";#it{p}_{T} (GeV/#it{c});#mu_{J/#psi} (MeV/c^{2})", len(ptEdges)-1, ptEdges))
        histStatWidthNeNe.append(ROOT.TH1D(f'histStatWidth_NeNe_{recoTags[iProd]}', ";#it{p}_{T} (GeV/#it{c});#sigma_{J/#psi} (MeV/c^{2})", len(ptEdges)-1, ptEdges))
        for iPt, (jpsiRawYield,jpsiStatRawYield,jpsiSystRawYield) in enumerate(zip(jpsiRawYields,jpsiStatRawYields,jpsiSystRawYields)):
            histStatRawYieldNeNe[iProd].SetBinContent(iPt+1, jpsiRawYield)
            histStatRawYieldNeNe[iProd].SetBinError(iPt+1, jpsiStatRawYield)
            histSystRawYieldNeNe[iProd].SetBinContent(iPt+1, jpsiRawYield)
            histSystRawYieldNeNe[iProd].SetBinError(iPt+1, jpsiSystRawYield)

            histStatMeanNeNe[iProd].SetBinContent(iPt+1, jpsiMeans[iPt])
            histStatMeanNeNe[iProd].SetBinError(iPt+1, jpsiStatMeans[iPt])

            histStatWidthNeNe[iProd].SetBinContent(iPt+1, jpsiWidths[iPt])
            histStatWidthNeNe[iProd].SetBinError(iPt+1, jpsiStatWidths[iPt])

    for iProd in range(0, len(recoTags)):
        SetHistStat(histStatRawYieldNeNe[iProd], styles[iProd], colors["NeNe"]+iProd)
        SetHistSyst(histSystRawYieldNeNe[iProd], styles[iProd], colors["NeNe"]+iProd)
        SetHistStat(histStatMeanNeNe[iProd], styles[iProd], colors["NeNe"]+iProd)
        SetHistStat(histStatWidthNeNe[iProd], styles[iProd], colors["NeNe"]+iProd)
        histStatRawYieldNeNe[iProd].Scale(1. / nCollsNeNe[iProd], "WIDTH")
        histSystRawYieldNeNe[iProd].Scale(1. / nCollsNeNe[iProd], "WIDTH")
        histStatMeanNeNe[iProd].Scale(1000.)
        histStatWidthNeNe[iProd].Scale(1000.)

        SetHistStat(histCentralityNeNe[iProd], styles[iProd], colors["NeNe"]+iProd)
        histCentralityNeNe[iProd].Scale(4. / histCentralityNeNe[iProd].Integral())
        histCentralityNeNe[iProd].SetTitle("")
        histCentralityNeNe[iProd].SetName(f'histCentrality_NeNe_{recoTags[iProd]}')
        histCentralityNeNe[iProd].SetTitleSize(0.045,"X")
        histCentralityNeNe[iProd].SetTitleSize(0.045,"Y")
        histCentralityNeNe[iProd].SetLabelSize(0.045,"X")
        histCentralityNeNe[iProd].SetLabelSize(0.045,"Y")
        histCentralityNeNe[iProd].SetTitleOffset(1.2,"X")
        histCentralityNeNe[iProd].SetTitleOffset(1.35,"Y")

    canvasRawYieldNeNe = ROOT.TCanvas("canvasRawYieldNeNe", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histStatRawYieldNeNe[0].GetYaxis().SetRangeUser(0.01 * histStatRawYieldNeNe[0].GetMaximum(), 2 * histStatRawYieldNeNe[0].GetMaximum())
    histStatRawYieldNeNe[0].Draw("EP")
    for iProd in range(0, len(recoTags)):
        if iProd == 0:
            histSystRawYieldNeNe[iProd].Draw("E2P SAME")
        else:
            histStatRawYieldNeNe[iProd].Draw("EP SAME")
            histSystRawYieldNeNe[iProd].Draw("E2P SAME")
    canvasRawYieldNeNe.Update()

    canvasMeanNeNe = ROOT.TCanvas("canvasMeanNeNe", "", 800, 600)
    histStatMeanNeNe[0].GetYaxis().SetRangeUser(3050, 3150)
    for iProd in range(0, len(recoTags)):
        if iProd == 0:
            histStatMeanNeNe[iProd].Draw("EP")
        else:
            histStatMeanNeNe[iProd].Draw("EP SAME")
    lineJpsiMass.Draw()
    canvasMeanNeNe.Update()

    canvasWidthNeNe = ROOT.TCanvas("canvasWidthNeNe", "", 800, 600)
    histStatWidthNeNe[0].GetYaxis().SetRangeUser(0, 150)
    for iProd in range(0, len(recoTags)):
        if iProd == 0:
            histStatWidthNeNe[iProd].Draw("EP")
        else:
            histStatWidthNeNe[iProd].Draw("EP SAME")
    canvasWidthNeNe.Update()


    canvasCentralitySummary = ROOT.TCanvas("canvasCentralitySummary", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    legendCentralitySummary = ROOT.TLegend(0.20, 0.70, 0.90, 0.90, " ", "brNDC")
    SetLegend(legendCentralitySummary)
    legendCentralitySummary.SetTextSize(0.045)
    legendCentralitySummary.SetNColumns(3)
    for iProd in range(0, len(recoTags)):
        legendCentralitySummary.AddEntry(histCentralityOO[iProd], f'OO {recoTags[iProd]}', "EP")
        legendCentralitySummary.AddEntry(histCentralitypO[iProd], f'pO {recoTags[iProd]}', "EP")
        legendCentralitySummary.AddEntry(histCentralityNeNe[iProd], f'Ne-Ne {recoTags[iProd]}', "EP")

        if iProd == 0:
            histCentralityOO[iProd].GetYaxis().SetRangeUser(0.1 * histCentralityOO[0].GetMaximum(), 5 * histCentralityOO[0].GetMaximum())
            histCentralityOO[iProd].Draw("EP")
        else:
            histCentralityOO[iProd].Draw("EP SAME")
        histCentralitypO[iProd].Draw("EP SAME")
        histCentralityNeNe[iProd].Draw("EP SAME")
    legendCentralitySummary.Draw("SAME")
    canvasCentralitySummary.Update()

    canvasRawYieldSummary = ROOT.TCanvas("canvasRawYieldSummary", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histStatRawYieldOO[1].GetYaxis().SetRangeUser(0.005 * histStatRawYieldOO[0].GetMaximum(), 5 * histStatRawYieldOO[0].GetMaximum())
    histStatRawYieldOO[1].Draw("EP")
    histSystRawYieldOO[1].Draw("E2P SAME")
    histStatRawYieldpO[1].Draw("EP SAME")
    histSystRawYieldpO[1].Draw("E2P SAME")
    histStatRawYieldNeNe[1].Draw("EP SAME")
    histSystRawYieldNeNe[1].Draw("E2P SAME")

    legendRawYieldSummary = ROOT.TLegend(0.70, 0.60, 0.90, 0.85, " ", "brNDC")
    SetLegend(legendRawYieldSummary)
    legendRawYieldSummary.SetTextSize(0.045)
    legendRawYieldSummary.AddEntry(histSystRawYieldpO[1], "pO", "FP")
    legendRawYieldSummary.AddEntry(histSystRawYieldOO[1], "OO", "FP")
    legendRawYieldSummary.AddEntry(histSystRawYieldNeNe[1], "Ne-Ne", "FP")
    legendRawYieldSummary.Draw("SAME")

    latexTitle.DrawLatex(0.2, 0.87, "ALICE work in progress")
    latexTitle.DrawLatex(0.2, 0.80, "J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4, 0-90%")

    canvasRawYieldSummary.Update()



    canvasMeanSummary = ROOT.TCanvas("canvasMeanSummary", "", 800, 600)
    histStatMeanOO[1].GetYaxis().SetRangeUser(3050, 3150)
    histStatMeanOO[1].Draw("EP")
    histStatMeanpO[1].Draw("EP SAME")
    histStatMeanNeNe[1].Draw("EP SAME")
    lineJpsiMass.Draw()

    legendMeanSummary = ROOT.TLegend(0.70, 0.60, 0.90, 0.85, " ", "brNDC")
    SetLegend(legendMeanSummary)
    legendMeanSummary.SetTextSize(0.045)
    legendMeanSummary.AddEntry(histStatMeanpO[1], "pO", "FP")
    legendMeanSummary.AddEntry(histStatMeanOO[1], "OO", "FP")
    legendMeanSummary.AddEntry(histStatMeanNeNe[1], "Ne-Ne", "FP")
    legendMeanSummary.Draw("SAME")

    latexTitle.DrawLatex(0.2, 0.87, "ALICE work in progress")
    latexTitle.DrawLatex(0.2, 0.80, "J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4, 0-90%")

    canvasMeanSummary.Update()


    canvasWidthSummary = ROOT.TCanvas("canvasWidthSummary", "", 800, 600)
    histStatWidthOO[1].GetYaxis().SetRangeUser(0, 150)
    histStatWidthOO[1].Draw("EP")
    histStatWidthpO[1].Draw("EP SAME")
    histStatWidthNeNe[1].Draw("EP SAME")
    lineJpsiMass.Draw()

    legendWidthSummary = ROOT.TLegend(0.70, 0.30, 0.90, 0.50, " ", "brNDC")
    SetLegend(legendWidthSummary)
    legendWidthSummary.SetTextSize(0.045)
    legendWidthSummary.AddEntry(histStatWidthpO[1], "pO", "FP")
    legendWidthSummary.AddEntry(histStatWidthOO[1], "OO", "FP")
    legendWidthSummary.AddEntry(histStatWidthNeNe[1], "Ne-Ne", "FP")
    legendWidthSummary.Draw("SAME")

    latexTitle.DrawLatex(0.2, 0.87, "ALICE work in progress")
    latexTitle.DrawLatex(0.2, 0.80, "J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4, 0-90%")

    canvasWidthSummary.Update()

    input()

    if args.dOutName:
        canvasCentralitySummary.SaveAs(f'{args.dOutName}/CentralitySummary.pdf')
        canvasRawYieldSummary.SaveAs(f'{args.dOutName}/RawYieldSummary.pdf')
        canvasMeanSummary.SaveAs(f'{args.dOutName}/MeanSummary.pdf')
        canvasWidthSummary.SaveAs(f'{args.dOutName}/WidthSummary.pdf')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--run", help="Compute cross section", action="store_true")
    parser.add_argument("-o", "--output", dest="dOutName", type=str, help="Name of the output directory", required=False)
    args = parser.parse_args()
    

    with open(args.cfgFileName, 'r') as yml_cfg:
        inputCfg = yaml.load(yml_cfg, yaml.FullLoader)

    if args.run:
        qa(inputCfg)