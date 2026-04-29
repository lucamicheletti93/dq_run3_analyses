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
    #delptapt = [1,1,1,1,1,1,2]
    Taa = config["inputs"]["Taa"]

    print("***** Compute corrected yield *****")
    # Lumi computed with normalization.C  + check_normalization
    #lumi = 0.170285 # pb-1
    fInLumi = ROOT.TFile(config["inputs"]["fInLumi"], "READ")
    histLumi = fInLumi.Get("histLumi")
    lumi = histLumi.GetBinContent(1)
    histNevMinBias = fInLumi.Get("histNevMinBias")
    nevMinBias = histNevMinBias.GetBinContent(1)
    print("N. min bias events = ", nevMinBias)

    # Normalization: N_minbias
    histNevMinBiasCentr = fInLumi.Get("histnEvtsBcSelCentr_AfterCuts") 
    nCentrBins = histNevMinBiasCentr.GetNbinsX()
    centrLabels = ["0-10%","10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%","70-90%", "0-90%", "0-100%"]
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
    histStatPPrefXsec = fInPPrefVsPt.Get("histStatJpsiXsecInterp") # vs pT
    histSystPPrefXsec = fInPPrefVsPt.Get("histSystJpsiXsecInterp") # vs pT
    histStatPPrefXsecVsSqrts = fInPPrefVsPt.Get("histStatJpsiXsecInterpVsSqrts") # pT integrated
    histSystPPrefXsecVsSqrts = fInPPrefVsPt.Get("histSystJpsiXsecInterpVsSqrts") # pT integrated
    
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

    print("-------- Pt dependence --------")

    print("***** Extract raw yield *****")
    dfJpsiRawYieldVsPt = pd.read_csv(config["inputs"]["fInRawYieldVsPt"], sep=' ')
    ptMin = dfJpsiRawYieldVsPt["x_min"]
    ptMax = dfJpsiRawYieldVsPt["x_max"]
    ptCenters = (ptMax + ptMin) / 2.
    ptWidths = (ptMax - ptMin) / 2.
    ptSystWidths = np.repeat(0.2, 1)
    ptEdges = np.array([ptMin, ptMax])

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
    jpsi1SystRelPPrefXsecVsPt = (np.array(jpsiStatPPrefXsecVsCentr) / np.array(jpsiPPrefXsecVsCentr))
    jpsi2SystRelPPrefXsecVsPt = (np.array(jpsiSystPPrefXsecVsCentr) / np.array(jpsiPPrefXsecVsCentr))
    jpsiSystRelPPrefXsecVsPt = np.sqrt(jpsiSystPPrefXsecVsCentr**2 + jpsi2SystRelPPrefXsecVsPt**2)
    print(jpsiSystRelPPrefXsecVsPt)

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    # > > > Pt dependence < < < #
    jpsiXsecVsPt = (jpsiRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * lumi * 1e6 * (2 * ptWidths))
    jpsiStatXsecVsPt = (jpsiStatRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * lumi * 1e6 * (2 * ptWidths))
    jpsiSystXsecVsPt = (jpsiSystRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * lumi * 1e6 * (2 * ptWidths))

    nevMinBias_0_100 = 1.045 * nevCentrDict["0-100%"] 
    nevMinBias_0_100_Eff = 1.045 * nevCentrDict["0-100%"] / 0.935
    jpsiRaaVsPt = (jpsiRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * nevMinBias_0_100_Eff  * Taa * np.array(jpsiPPrefXsecVsCentr))
    jpsiStatRaaVsPt = jpsiRaaVsPt * (np.array(jpsiStatRawYieldVsPt) / np.array(jpsiRawYieldVsPt))
    jpsiSystRaaVsPt = jpsiRaaVsPt * np.sqrt(jpsiSystRelRawYieldVsPt**2 + jpsiSystRelPPrefXsecVsPt**2)

    jpsiRaaVsPt_Method2 = (jpsiRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * nevMinBias_0_100 * 16*16*(1/1151000) * np.array(jpsiPPrefXsecVsCentr))
    jpsiStatRaaVsPt_Method2 = jpsiRaaVsPt * (np.array(jpsiStatRawYieldVsPt) / np.array(jpsiRawYieldVsPt))
    jpsiSystRaaVsPt_Method2 = jpsiRaaVsPt * np.sqrt(jpsiSystRelRawYieldVsPt**2 + jpsiSystRelPPrefXsecVsPt**2)

    # Add all systematics contributions
    #jpsiSystXsecVsPt = jpsiXsecVsPt * np.sqrt((jpsiSystXsecVsPt / jpsiXsecVsPt)**2 + systRelBrJpsiToMuMu**2 + systRelTrackingEff**2 + systRelMatchingEff**2)

    graStatXsecVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiXsecVsPt), np.array(ptWidths), np.array(jpsiStatXsecVsPt))
    graSystXsecVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiXsecVsPt), np.array(ptSystWidths), np.array(jpsiSystXsecVsPt))

    graStatRaaVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt), np.array(ptWidths), np.array(jpsiStatRaaVsPt))
    graSystRaaVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt), np.array(ptSystWidths), np.array(jpsiSystRaaVsPt))

    graStatRaaVsPt_Method2 = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt_Method2), np.array(ptWidths), np.array(jpsiStatRaaVsPt))
    graSystRaaVsPt_Method2 = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt_Method2), np.array(ptSystWidths), np.array(jpsiSystRaaVsPt))

    SetGraStat(graStatXsecVsPt, 24, ROOT.kBlack)
    SetGraSyst(graSystXsecVsPt, 20, ROOT.kRed+1, 1, 2)

    SetGraStat(graStatRaaVsPt, 24, ROOT.kBlack)
    SetGraSyst(graSystRaaVsPt, 20, ROOT.kRed+1, 1, 2)

    SetGraStat(graStatRaaVsPt_Method2, 24, ROOT.kBlack)
    SetGraSyst(graSystRaaVsPt_Method2, 20, ROOT.kGreen+2, 1, 2)

    output = np.column_stack((ptCenters, jpsiRaaVsPt, jpsiStatRaaVsPt, jpsiSystRaaVsPt))

    np.savetxt("jpsi_RAA_vs_pT_OnlyOneBin_0_20_Train_581158_PbPbQuality_Plus_CorrectionPileUp.txt", output, header="pT_center(GeV/c)   RAA   stat_err   syst_err", fmt="%.3f   %.6f   %.6f   %.6f")
    

    canvasXsecVsPt = ROOT.TCanvas("canvasXsecVsPt", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    ROOT.gPad.SetLogy(True)
    histGridXsecVsPt = ROOT.TH2D("histGridXsecVsPt", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}", 100, 0, 20, 100, 0, 300)
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
    
    histGridRaaVsPt = ROOT.TH2D("histGridRaaVsPt",";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}",100, 0, 20, 100, 0, 1.5)
    histGridRaaVsPt.Draw()
    
    # --- Unity line ---
    lineUnityVsPt = ROOT.TLine(0., 1., 20., 1.)
    lineUnityVsPt.SetLineColor(ROOT.kGray+1)
    lineUnityVsPt.SetLineWidth(2)
    lineUnityVsPt.SetLineStyle(ROOT.kDashed)
    lineUnityVsPt.Draw("SAME")

    # -------------------------------
    # Common systematics (global box)
    # -------------------------------
    systmatchingeffJpsi = 0.02   
    systtrackeffJpsi = 0.015
    systBR = 0.005

    systTOO = 0.06
    systEffTVX = 0.06

    systLumi = 0.05

    systCommonTot_Method1 = np.sqrt(systmatchingeffJpsi**2 + systtrackeffJpsi**2 + systBR**2 + systTOO**2 + systEffTVX**2)
    systCommonTot_Method2 = np.sqrt(systmatchingeffJpsi**2 + systtrackeffJpsi**2 + systBR**2 + systLumi**2)

    xMinBox_Method1 = 7.6
    xMaxBox_Method1 = 7.8
    xMinBox_Method2 = 7.8
    xMaxBox_Method2 = 8.0
    yRef = 1.0

    yMinBox_Method1 = yRef * (1.0 - systCommonTot_Method1)
    yMaxBox_Method1 = yRef * (1.0 + systCommonTot_Method1)
    yMinBox_Method2 = yRef * (1.0 - systCommonTot_Method2)
    yMaxBox_Method2 = yRef * (1.0 + systCommonTot_Method2)

    boxCommon_Method1 = ROOT.TBox(xMinBox_Method1, yMinBox_Method1, xMaxBox_Method1, yMaxBox_Method1)
    boxCommon_Method1.SetFillColorAlpha(ROOT.kRed+1, 0.4)
    boxCommon_Method1.SetLineColor(ROOT.kRed+1)
    boxCommon_Method1.SetLineWidth(1)

    boxCommon_Method2 = ROOT.TBox(xMinBox_Method2, yMinBox_Method2, xMaxBox_Method2, yMaxBox_Method2)
    boxCommon_Method2.SetFillColorAlpha(ROOT.kGreen+2, 0.4)
    boxCommon_Method2.SetLineColor(ROOT.kGreen+2)
    boxCommon_Method2.SetLineWidth(1)

    boxCommon_Method1.Draw("same")
    boxCommon_Method2.Draw("EP same")

    # --- Draw graphs ---
    graStatRaaVsPt.Draw("EP SAME")
    graSystRaaVsPt.Draw("E2P SAME")
    #graStatRaaVsPt_Method2.Draw("EP SAME")
    #graSystRaaVsPt_Method2.Draw("E2P SAME")

    legendRaaVsPtVsRun2 = ROOT.TLegend(0.20, 0.70, 0.40, 0.87, " ", "brNDC")
    SetLegend(legendRaaVsPtVsRun2)
    legendRaaVsPtVsRun2.SetTextSize(0.045)
    legendRaaVsPtVsRun2.AddEntry(graSystRaaVsPt, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus100%", "FP")
    #legendRaaVsPtVsRun2.AddEntry(graSystRaaVsPt_Method2, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus100%, 2nd Method", "FP")
    legendRaaVsPtVsRun2.Draw("SAME")

    latexTitle.DrawLatex(0.20, 0.90, "ALICE Work In Progress")
    latexTitle.DrawLatex(0.20, 0.84, "J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    
    canvasRaaVsPt.Update()
    canvasRaaVsPt.SaveAs("plots/Jpsi_RAA_vs_Pt.pdf")
    canvasRaaVsPt.SaveAs("plots/Jpsi_RAA_vs_Pt.png")

    # --- Different between the two Methods ---
    diffRaaVsPt = np.array(jpsiRaaVsPt) - np.array(jpsiRaaVsPt_Method2)

    # For now: the error on the difference "diffRaaVsPt" is the quadratic sum of the statistical uncertainties
    diffStatRaaVsPt = np.sqrt(np.array(jpsiStatRaaVsPt)**2 + np.array(jpsiStatRaaVsPt_Method2)**2)

    graDiffRaaVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), diffRaaVsPt, np.array(ptWidths), diffStatRaaVsPt)

    graDiffRaaVsPt.SetMarkerStyle(21)
    graDiffRaaVsPt.SetMarkerColor(ROOT.kBlue)
    graDiffRaaVsPt.SetLineColor(ROOT.kBlue)

    canvasDiffRaaVsPt = ROOT.TCanvas("canvasDiffRaaVsPt", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)

    histGridDiffRaaVsPt = ROOT.TH2D("histGridDiffRaaVsPt", ";#it{p}_{T} (GeV/#it{c});R_{AA}^{1stMethod} - R_{AA}^{2ndMethod}", 100, 0, 20, 100, -0.5, 0.5)
    histGridDiffRaaVsPt.Draw()

    lineUnity_DiffRaa = ROOT.TLine(0., 0., 8., 0.)
    lineUnity_DiffRaa.SetLineColor(ROOT.kGray+1)
    lineUnity_DiffRaa.SetLineWidth(2)
    lineUnity_DiffRaa.SetLineStyle(ROOT.kDashed)
    lineUnity_DiffRaa.Draw("SAME")

    graDiffRaaVsPt.Draw("EP SAME")
    canvasDiffRaaVsPt.Update()
    canvasDiffRaaVsPt.SaveAs("plots/Jpsi_RAA_vs_Pt_Diff_Two_Methods.pdf")

    # --- Compare with Previous Results: Method 1 Here and Method 1 in Previous Results ---
    
    jpsiRaaVsPt_Train_591157_Sel8 = np.array([0.726947])
    jpsiStatRaaVsPt_Train_591157_Sel8 = np.array([0.005382])
    jpsiSystRaaVsPt_Train_591157_Sel8 = np.array([0.045589])

    jpsiRaaVsPt_Train_601190_Sel8PlusNoSameBunch = np.array([0.741551])
    jpsiStatRaaVsPt_Train_601190_Sel8PlusNoSameBunch = np.array([0.00582])
    jpsiSystRaaVsPt_Train_601190_Sel8PlusNoSameBunch = np.array([0.046545])

    graStatRaaVsPt_Train_591157_Sel8 = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), jpsiRaaVsPt_Train_591157_Sel8, np.array(ptWidths), jpsiStatRaaVsPt_Train_591157_Sel8)
    graSystRaaVsPt_Train_591157_Sel8 = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), jpsiRaaVsPt_Train_591157_Sel8, np.array(ptSystWidths), jpsiSystRaaVsPt_Train_591157_Sel8)
    SetGraStat(graStatRaaVsPt_Train_591157_Sel8, 24, ROOT.kBlack)
    SetGraSyst(graSystRaaVsPt_Train_591157_Sel8, 20, ROOT.kRed+1,1,2)

    graStatRaaVsPt_Train_601190_Sel8PlusNoSameBunch = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), jpsiRaaVsPt_Train_601190_Sel8PlusNoSameBunch, np.array(ptWidths), jpsiStatRaaVsPt_Train_601190_Sel8PlusNoSameBunch)
    graSystRaaVsPt_Train_601190_Sel8PlusNoSameBunch = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), jpsiRaaVsPt_Train_601190_Sel8PlusNoSameBunch, np.array(ptSystWidths), jpsiSystRaaVsPt_Train_601190_Sel8PlusNoSameBunch)
    SetGraStat(graStatRaaVsPt_Train_601190_Sel8PlusNoSameBunch, 24, ROOT.kBlack)
    SetGraSyst(graSystRaaVsPt_Train_601190_Sel8PlusNoSameBunch, 20, ROOT.kGreen+1,1,2)

    canvasCompare = ROOT.TCanvas("canvasCompare", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)

    histGrid = ROOT.TH2D("histGrid",";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}", 100, 0, 20, 100, 0, 1.5)
    histGrid.Draw()

    lineUnityVsPt.Draw("SAME")

    boxCommon_Method1.Draw("same")

    SetGraSyst(graSystRaaVsPt, 20, ROOT.kBlue+1, 1, 2)

    graSystRaaVsPt.Draw("E2P SAME")
    graStatRaaVsPt.Draw("EP SAME")
    graStatRaaVsPt_Train_591157_Sel8.Draw("EP SAME")
    graSystRaaVsPt_Train_591157_Sel8.Draw("E2P SAME")
    graStatRaaVsPt_Train_601190_Sel8PlusNoSameBunch.Draw("EP SAME")
    graSystRaaVsPt_Train_601190_Sel8PlusNoSameBunch.Draw("E2P SAME")

    legend = ROOT.TLegend(0.2, 0.7, 0.5, 0.87, "", "brNDC")
    SetLegend(legend)
    legend.AddEntry(graSystRaaVsPt_Train_591157_Sel8, "eventSel8", "FP")
    legend.AddEntry(graSystRaaVsPt_Train_601190_Sel8PlusNoSameBunch, "eventSel8NoSameBunch", "FP")
    legend.AddEntry(graSystRaaVsPt, "eventStandardSel8PbPbQuality", "FP")

    
    legend.Draw("SAME")

    canvasCompare.Update()
    canvasCompare.SaveAs("plots/Jpsi_RAA_vs_Pt_Method1_Comparison_here_and_previous_results.pdf")

    ratio_Sel8_over_NoSame = np.array(jpsiRaaVsPt_Train_591157_Sel8 / jpsiRaaVsPt_Train_601190_Sel8PlusNoSameBunch,dtype='float64')

    ratio_Sel8_over_PbPbQ = np.array(jpsiRaaVsPt_Train_591157_Sel8 / jpsiRaaVsPt,dtype='float64')

    ratio_NoSame_over_PbPbQ = np.array(jpsiRaaVsPt_Train_601190_Sel8PlusNoSameBunch / jpsiRaaVsPt,dtype='float64')

    zeroErrors = np.zeros(len(ptCenters), dtype='float64')
    graRatio_Sel8_NoSame = ROOT.TGraphErrors(len(ptCenters),np.array(ptCenters),ratio_Sel8_over_NoSame,np.array(ptWidths),zeroErrors)
    graRatio_Sel8_PbPbQ = ROOT.TGraphErrors(len(ptCenters),np.array(ptCenters),ratio_Sel8_over_PbPbQ,np.array(ptWidths),zeroErrors)
    graRatio_NoSame_PbPbQ = ROOT.TGraphErrors(len(ptCenters),np.array(ptCenters),ratio_NoSame_over_PbPbQ,np.array(ptWidths),zeroErrors)

    graRatio_Sel8_NoSame.SetMarkerStyle(20)
    graRatio_Sel8_NoSame.SetMarkerColor(ROOT.kBlack+1)
    graRatio_Sel8_NoSame.SetLineColor(ROOT.kBlack+1)

    graRatio_Sel8_PbPbQ.SetMarkerStyle(21)
    graRatio_Sel8_PbPbQ.SetMarkerColor(ROOT.kBlue+1)
    graRatio_Sel8_PbPbQ.SetLineColor(ROOT.kBlue+1)

    graRatio_NoSame_PbPbQ.SetMarkerStyle(22)
    graRatio_NoSame_PbPbQ.SetMarkerColor(ROOT.kGray+2)
    graRatio_NoSame_PbPbQ.SetLineColor(ROOT.kGray+2)

    canvasRatio = ROOT.TCanvas("canvasRatio", "", 800, 600)

    histRatioGrid = ROOT.TH2D("histRatioGrid",";#it{p}_{T} (GeV/#it{c});Ratio",100, 0, 20,100, 0.9, 1.1)

    histRatioGrid.Draw()

    lineUnityVsPt.Draw("SAME")   # se già definita prima

    graRatio_Sel8_NoSame.Draw("P SAME")
    graRatio_Sel8_PbPbQ.Draw("P SAME")
    graRatio_NoSame_PbPbQ.Draw("P SAME")

    legendRatio = ROOT.TLegend(0.2, 0.75, 0.5, 0.88, "", "brNDC")
    SetLegend(legendRatio)

    legendRatio.AddEntry(graRatio_Sel8_NoSame, "Sel8 / Sel8NoSameBunch", "P")
    legendRatio.AddEntry(graRatio_Sel8_PbPbQ, "Sel8 / PbPbQuality", "P")
    legendRatio.AddEntry(graRatio_NoSame_PbPbQ, "Sel8NoSameBunch / PbPbQuality", "P")

    legendRatio.Draw("SAME")
    lineUnity_DiffRaa.Draw("SAME")

    canvasRatio.Update()
    canvasRatio.SaveAs("plots/Jpsi_RAA_vs_Pt_Method1_Comparison_here_and_previous_results_RATIO.pdf")
    

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


if __name__ == '__main__':
    main()