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
        raa(inputCfg)

def raa(config):
    """
    function to compute the RAA (nuclear modification factor)
    """
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
    centrWidths = (centrMax - centrMin) / 2.
    centrSystWidths = np.repeat(0.2, len(centrWidths))
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

    histStatAxeVsCentr = ROOT.TH1D("histStatAxeVsCentr", ";#it{p}_{T} (GeV/#it{c});A#times#varepsilon", len(centrEdges)-1, centrEdges)

    SetHistStat(histStatAxeVsCentr, 20, ROOT.kAzure+4)

    for iBin in range(0, len(centrMin)):
        histStatAxeVsCentr.SetBinContent(iBin+1, jpsiAxeVsCentr[iBin])
        histStatAxeVsCentr.SetBinError(iBin+1, jpsiStatAxeVsCentr[iBin])

    canvasAxeVsCentr = ROOT.TCanvas("canvasAxeVsCentr", "", 800, 600)
    histStatAxeVsCentr.Draw("EP SAME")
    canvasAxeVsCentr.Update()

    print(jpsiRawYieldVsCentr)
    print(jpsiAxeVsCentr)
    print(nevMinBias)
    print(TaaVsCentr)
    print(jpsiPPrefXsecVsCentr)

    jpsiRaaVsCentr = (jpsiRawYieldVsCentr) / (jpsiAxeVsCentr * BrJpsiToMuMu * (nevMinBias/10) * np.array(TaaVsCentr) * jpsiPPrefXsecVsCentr)
    jpsiStatRaaVsCentr = jpsiRaaVsCentr * (np.array(jpsiStatRawYieldVsCentr) / np.array(jpsiRawYieldVsCentr))
    jpsiSystRaaVsCentr = jpsiRaaVsCentr * (np.array(jpsiSystRawYieldVsCentr) / np.array(jpsiRawYieldVsCentr))

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
    histGridRaaVsCentr = ROOT.TH2D("histGridRaaVsCentr", ";Centrality (%);#it{R}_{AA}", 100, 0, 90, 100, 0, 1.7)
    histGridRaaVsCentr.Draw()
    lineUnityVsCentr.Draw("SAME")
    #graStatRaaVsCentrRun2.Draw("EP SAME")
    #graSystRaaVsCentrRun2.Draw("E2P SAME")
    graStatRaaVsCentr.Draw("EP SAME")
    graSystRaaVsCentr.Draw("E2P SAME")

    legendRaaVsCentrVsRun2 = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsCentrVsRun2)
    legendRaaVsCentrVsRun2.SetTextSize(0.045)
    legendRaaVsCentrVsRun2.AddEntry(graStatRaaVsCentr, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus100%", "FP")
    #legendRaaVsCentrVsRun2.AddEntry(graStatRaaVsCentrRun2, "#sqrt{#it{s}} = 5.02 TeV, Pb-Pb, 0#minus90%", "FP")
    legendRaaVsCentrVsRun2.Draw("SAME")

    latexTitle.DrawLatex(0.22, 0.87, "ALICE Preliminary, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    canvasRaaVsCentrVsRun2.Update()




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

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    # List of systematics computed externally
    """
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
    """

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

    jpsiRaaVsPt = (jpsiRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * nevMinBias * (2 * ptWidths) * Taa * np.array(jpsiPPrefXsecVsPt))
    jpsiStatRaaVsPt = jpsiRaaVsPt * (np.array(jpsiStatRawYieldVsPt) / np.array(jpsiRawYieldVsPt))
    jpsiSystRaaVsPt = jpsiRaaVsPt * np.sqrt(jpsiSystRelRawYieldVsPt**2 + jpsiSystRelPPrefXsecVsPt**2)

    # Add all systematics contributions
    #jpsiSystXsecVsPt = jpsiXsecVsPt * np.sqrt((jpsiSystXsecVsPt / jpsiXsecVsPt)**2 + systRelBrJpsiToMuMu**2 + systRelTrackingEff**2 + systRelMatchingEff**2)

    graStatXsecVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiXsecVsPt), np.array(ptWidths), np.array(jpsiStatXsecVsPt))
    graSystXsecVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiXsecVsPt), np.array(ptSystWidths), np.array(jpsiSystXsecVsPt))

    graStatRaaVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt), np.array(ptWidths), np.array(jpsiStatRaaVsPt))
    graSystRaaVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt), np.array(ptSystWidths), np.array(jpsiSystRaaVsPt))

    SetGraStat(graStatXsecVsPt, 20, ROOT.kRed+1)
    SetGraSyst(graSystXsecVsPt, 20, ROOT.kRed+1)

    SetGraStat(graStatRaaVsPt, 20, ROOT.kRed+1)
    SetGraSyst(graSystRaaVsPt, 20, ROOT.kRed+1)


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

    canvasXsecVsPt = ROOT.TCanvas("canvasXsecVsPt", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histStatXsecVsPt.Draw("EP")
    histSystXsecVsPt.Draw("E2P SAME")
    canvasXsecVsPt.Update()

    # Compare to Run 2 results
    print("***** Extract RAA from Run 2 *****")
    filePathJpsi = "HEP_Data/jpsi_raa_PbPb_5TeV.yaml"
    ptCentersRun2, ptWidthsRun2, jpsiRaaVsPtRun2, jpsiStatRaaVsPtRun2, jpsiSystRaaVsPtRun2 = ExtractFromYaml(filePathJpsi)
    jpsiPtSystWidthsRun2 = np.repeat(0.2, len(ptCentersRun2))

    graStatRaaVsPtRun2 = ROOT.TGraphErrors(len(ptCentersRun2), np.array(ptCentersRun2), np.array(jpsiRaaVsPtRun2), np.array(ptWidthsRun2), np.array(jpsiStatRaaVsPtRun2))
    graSystRaaVsPtRun2 = ROOT.TGraphErrors(len(ptCentersRun2), np.array(ptCentersRun2), np.array(jpsiRaaVsPtRun2), jpsiPtSystWidthsRun2, np.array(jpsiSystRaaVsPtRun2))

    SetGraStat(graStatRaaVsPtRun2, 20, ROOT.kGray+2)
    SetGraSyst(graSystRaaVsPtRun2, 20, ROOT.kGray+2)

    lineUnityVsPt = ROOT.TLine(0., 1., 12., 1.)
    lineUnityVsPt.SetLineColor(ROOT.kGray+1)
    lineUnityVsPt.SetLineWidth(2)
    lineUnityVsPt.SetLineStyle(ROOT.kDashed)

    canvasRaaVsPtVsRun2 = ROOT.TCanvas("canvasRaaVsPtVsRun2", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    histGridRaaVsPt = ROOT.TH2D("histGridRaaVsPt", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}", 100, 0, 12, 100, 0, 1.7)
    histGridRaaVsPt.Draw()
    lineUnityVsPt.Draw("SAME")
    graStatRaaVsPtRun2.Draw("EP SAME")
    graSystRaaVsPtRun2.Draw("E2P SAME")
    graStatRaaVsPt.Draw("EP SAME")
    graSystRaaVsPt.Draw("E2P SAME")

    legendRaaVsPtVsRun2 = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsPtVsRun2)
    legendRaaVsPtVsRun2.SetTextSize(0.045)
    legendRaaVsPtVsRun2.AddEntry(graStatRaaVsPt, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus100%", "FP")
    legendRaaVsPtVsRun2.AddEntry(graStatRaaVsPtRun2, "#sqrt{#it{s}} = 5.02 TeV, Pb-Pb, 0#minus90%", "FP")
    legendRaaVsPtVsRun2.Draw("SAME")

    latexTitle.DrawLatex(0.22, 0.87, "ALICE Preliminary, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")

    canvasRaaVsPtVsRun2.Update()

    input()

    print(f"***** Write results in {config["outputs"]["fOut"]} *****")
    fOut = ROOT.TFile(config["outputs"]["fOut"], "RECREATE")
    histStatRawYieldVsPt.Write()
    histSystRawYieldVsPt.Write()
    histStatAxeVsPt.Write()
    histStatXsecVsPt.Write()
    histSystXsecVsPt.Write()

    fOut.mkdir("systematics")
    fOut.cd("systematics")
    histSystRelRawYieldVsPt.Write("syst_raw_yield_vs_pt")
    #histSystRelLumiVsPt.Write("syst_lumi_vs_pt")
    #histSystRelBrJpsiToMuMuVsPt.Write("syst_br_vs_pt")
    #histSystRelTrackingEffVsPt.Write("syst_tracking_eff_vs_pt")
    #histSystRelMatchingEffVsPt.Write("syst_matching_eff_vs_pt")
    #histSystRelMcRealisticVsPt.Write("syst_mc_realisticness_vs_pt")
    fOut.Close()

    canvasRawYieldVsPt.SaveAs("figures/raw_yield/raw_yeild_jpsi.pdf")
    canvasAxeVsPt.SaveAs("figures/axe/axe_jpsi.pdf")
    canvasRaaVsPtVsRun2.SaveAs("figures/raa/pt_differential_raa_jpsi.pdf")


if __name__ == '__main__':
    main()