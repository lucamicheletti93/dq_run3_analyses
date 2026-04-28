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
    # delptapt = [1,1,1,1,1,1,2]
    Taa = config["inputs"]["Taa"]

    print("***** Compute corrected yield *****")
    # Lumi computed with normalization.C  + check_normalization
    fInLumi = ROOT.TFile(config["inputs"]["fInLumi"], "READ")
    histLumi = fInLumi.Get("histLumi")
    lumi = histLumi.GetBinContent(1)
    # histNevMinBias = fInLumi.Get("histNevMinBias")
    # nevMinBias = histNevMinBias.GetBinContent(1)
    # print("N. min bias events = ", nevMinBias)

    # Normalization: N_minbias
    histNevMinBiasCentr = fInLumi.Get("histnEvtsBcSelCentr_AfterCuts_TMaker") # AfterCuts, YES PileUp
    #histNevMinBiasCentr = fInLumi.Get("histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp") # AfterCuts, NO PileUp
    nCentrBins = histNevMinBiasCentr.GetNbinsX()
    centrLabels = ["0-10%","10-20%", "20-30%", "30-40%", "40-50%", "50-60%","60-70","70-90%", "90-100%", "0-90%", "0-100%"]
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

    SetHistStat(histStatPPrefXsec, 20, ROOT.kAzure+4)
    SetHistSyst(histSystPPrefXsec, 20, ROOT.kAzure+4)

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
    ptWidths_zero = np.zeros(len(ptWidths))
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
    jpsiSystRelRawYieldVsPt = (np.array(jpsiSystRawYieldVsPt) / np.array(jpsiRawYieldVsPt))
    print(jpsiSystRelRawYieldVsPt)

    print("Relative systematic on pp reference")
    jpsi1SystRelPPrefXsecVsPt = (np.array(jpsiStatPPrefXsecVsPt) / np.array(jpsiPPrefXsecVsPt))
    jpsi2SystRelPPrefXsecVsPt = (np.array(jpsiSystPPrefXsecVsPt) / np.array(jpsiPPrefXsecVsPt))
    jpsiSystRelPPrefXsecVsPt = np.sqrt(jpsi1SystRelPPrefXsecVsPt**2 + jpsi2SystRelPPrefXsecVsPt**2)
    print(jpsiSystRelPPrefXsecVsPt)

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    # > > > Pt dependence < < < #

    # RAA
    nevMinBias_0_100 = nevCentrDict["0-100%"] # NOT Divided by the TVX efficiency -> SigmaMethod
    nevMinBias_0_100_Eff = nevCentrDict["0-100%"] / 0.935 # Divided by the TVX efficiency -> GlauberMethod
    jpsiRaaVsPt_MethodGlauber = (jpsiRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * nevMinBias_0_100_Eff * (2 * ptWidths) * Taa * np.array(jpsiPPrefXsecVsPt))
    jpsiStatRaaVsPt_MethodGlauber = jpsiRaaVsPt_MethodGlauber * (np.array(jpsiStatRawYieldVsPt) / np.array(jpsiRawYieldVsPt))
    jpsiSystRaaVsPt_MethodGlauber = jpsiRaaVsPt_MethodGlauber * np.sqrt(jpsiSystRelRawYieldVsPt**2 + jpsiSystRelPPrefXsecVsPt**2)
    #print("jpsiSystRelPPrefXsecVsPt**2")
    #print(jpsiSystRelPPrefXsecVsPt**2)
    #print("jpsiSystRelRawYieldVsPt**2")
    #print(jpsiSystRelRawYieldVsPt**2)

    jpsiRaaVsPt_MethodSigma = (jpsiRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * nevMinBias_0_100 * (2 * ptWidths) * 16*16*(1/1151000) * np.array(jpsiPPrefXsecVsPt))
    jpsiStatRaaVsPt_MethodSigma = jpsiRaaVsPt_MethodGlauber * (np.array(jpsiStatRawYieldVsPt) / np.array(jpsiRawYieldVsPt))
    jpsiSystRaaVsPt_MethodSigma = jpsiRaaVsPt_MethodGlauber * np.sqrt(jpsiSystRelRawYieldVsPt**2 + jpsiSystRelPPrefXsecVsPt**2)


    graStatRaaVsPt_MethodGlauber = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt_MethodGlauber), np.array(ptWidths), np.array(jpsiStatRaaVsPt_MethodGlauber))
    graSystRaaVsPt_MethodGlauber = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt_MethodGlauber), np.array(ptSystWidths), np.array(jpsiSystRaaVsPt_MethodGlauber))

    graStatRaaVsPt_MethodSigma = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt_MethodSigma), np.array(ptWidths_zero), np.array(jpsiStatRaaVsPt_MethodSigma))
    graSystRaaVsPt_MethodSigma = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt_MethodSigma), np.array(ptSystWidths), np.array(jpsiSystRaaVsPt_MethodSigma))

    SetGraStat(graStatRaaVsPt_MethodGlauber, 24, ROOT.kBlack)
    SetGraSyst(graSystRaaVsPt_MethodGlauber, 20, ROOT.kBlack, 1, 2)

    SetGraStat(graStatRaaVsPt_MethodSigma, 24, ROOT.kBlack)
    SetGraSyst(graSystRaaVsPt_MethodSigma, 20, ROOT.kRed+1, 1, 2)

    output = np.column_stack((ptCenters, jpsiRaaVsPt_MethodSigma, jpsiStatRaaVsPt_MethodSigma, jpsiSystRaaVsPt_MethodSigma))
    np.savetxt("SQM2026/ROO_Values/Train_706901_Sel8/jpsi_RAA_vs_pT_Train_607901_MethodSigma_Sel8_so_YESPileUpCorr_ppref_UNC_5.02.txt", output, header="pT_center(GeV/c)   RAA   stat_err   syst_err", fmt="%.3f   %.6f   %.6f   %.6f")

    # --- Canvas setup ---
    canvasRaaVsPt = ROOT.TCanvas("canvasRaaVsPt", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    
    histGridRaaVsPt = ROOT.TH2D("histGridRaaVsPt",";#it{p}_{T} (GeV/#it{c});#it{R}_{OO}",100, 0, 8, 100, 0, 1.5)
    histGridRaaVsPt.Draw()
    
    # --- Unity line ---
    lineUnityVsPt = ROOT.TLine(0., 1., 8., 1.)
    lineUnityVsPt.SetLineColor(ROOT.kGray+1)
    lineUnityVsPt.SetLineWidth(2)
    lineUnityVsPt.SetLineStyle(ROOT.kDashed)
    lineUnityVsPt.Draw("SAME")

    # -------------------------------
    # Common systematics (global box)
    # -------------------------------
    systinputshape = 0.025
    systmatchingeffJpsi = 0.02   
    systtrackeffJpsi = 0.025
    systBR = 0.005

    systTOO = 0.06
    systEffTVX = 0.06

    systLumi = 0.03  

    systCommonTot_MethodGlauber = np.sqrt(systinputshape**2+systmatchingeffJpsi**2 + systtrackeffJpsi**2 + systBR**2 + systTOO**2 + systEffTVX**2)
    systCommonTot_MethodSigma = np.sqrt(systinputshape**2+systmatchingeffJpsi**2 + systtrackeffJpsi**2 + systBR**2 + systLumi**2)
    print("systCommonTot_MethodSigma")
    print(systCommonTot_MethodSigma)
    xMinBox_MethodGlauber = 7.6
    xMinBox_MethodGaluber = 7.8
    xMinBox_MethodSigma = 7.8
    xMaxBox_MethodSigma = 8.0
    yRef = 1.0

    yMinBox_MethodGlauber = yRef * (1.0 - systCommonTot_MethodGlauber)
    yMaxBox_MethodGlauber = yRef * (1.0 + systCommonTot_MethodGlauber)
    yMinBox_MethodSigma = yRef * (1.0 - systCommonTot_MethodSigma)
    yMaxBox_MethodSigma = yRef * (1.0 + systCommonTot_MethodSigma)

    boxCommon_MethodGlauber = ROOT.TBox(xMinBox_MethodGlauber, yMinBox_MethodGlauber, xMinBox_MethodGaluber, yMaxBox_MethodGlauber)
    boxCommon_MethodGlauber.SetFillColorAlpha(ROOT.kBlack+1, 1)
    boxCommon_MethodGlauber.SetLineColor(ROOT.kBlack+1)
    boxCommon_MethodGlauber.SetLineWidth(1)

    boxCommon_MethodSigma = ROOT.TBox(xMinBox_MethodSigma, yMinBox_MethodSigma, xMaxBox_MethodSigma, yMaxBox_MethodSigma)
    boxCommon_MethodSigma.SetFillColorAlpha(ROOT.kRed+1, 1)
    boxCommon_MethodSigma.SetLineColor(ROOT.kRed+1)
    boxCommon_MethodSigma.SetLineWidth(1)

    #boxCommon_MethodGlauber.Draw("same")
    boxCommon_MethodSigma.Draw("EP same")

    # --- Draw graphs ---
    #graStatRaaVsPt_MethodGlauber.Draw("EP SAME")
    #graSystRaaVsPt_MethodGlauber.Draw("E2P SAME")
    #graStatRaaVsPt_MethodGlauber.Fit("pol0")
    graStatRaaVsPt_MethodSigma.Draw("EP SAME")
    graSystRaaVsPt_MethodSigma.Draw("E2P SAME")
    #graSystRaaVsPt_MethodSigma.Fit("pol0")
    #f = graSystRaaVsPt_MethodSigma.GetFunction("pol0")
    #f.SetLineColor(ROOT.kBlack)

    legendRaaVsPt = ROOT.TLegend(0.65, 0.25, 0.85, 0.42, " ", "brNDC")
    SetLegend(legendRaaVsPt)
    legendRaaVsPt.SetTextSize(0.045)
    ##legendRaaVsPt.AddEntry(graSystRaaVsPt_MethodGlauber, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus100%, Glauber ", "FP")
    #legendRaaVsPt.AddEntry(graSystRaaVsPt_MethodSigma, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus100%, Sigma", "FP")
    legendRaaVsPt.AddEntry(graStatRaaVsPt_MethodSigma, "Stat. Uncert.", "PL")
    legendRaaVsPt.AddEntry(graSystRaaVsPt_MethodSigma, "Syst. Uncert.", "F")
    legendRaaVsPt.Draw("SAME")

    latexTitle.DrawLatex(0.20, 0.86, "ALICE Preliminary, OO,")
    latexTitle.DrawLatex(0.57, 0.86, "#sqrt{#it{s}} = 5.36 TeV")
    latexTitle.DrawLatex(0.20, 0.78, "J/#psi")
    latexTitle.DrawLatex(0.26, 0.78, "#rightarrow")
    latexTitle.DrawLatex(0.30, 0.78, "#mu^{#plus}#mu^{#minus},")
    latexTitle.DrawLatex(0.38, 0.78, "2.5<")
    latexTitle.DrawLatex(0.45, 0.78, "#it{y}")
    latexTitle.DrawLatex(0.47, 0.78, "<4,")
    latexTitle.DrawLatex(0.53, 0.78, "0#minus100%")
    canvasRaaVsPt.Update()
    canvasRaaVsPt.SaveAs("SQM2026/Jpsi_RAA_vs_Pt.pdf")
    canvasRaaVsPt.SaveAs("SQM2026/Jpsi_RAA_vs_Pt.png")


    ######################################################
    # ***** Comparison with pO *****
    ptCenters_pO = np.array([0.5, 1.5, 2.5, 3.5, 5, 7])
    ptWidths_pO = np.array([0.5,  0.5, 0.5, 0.5, 1, 1])
    ptSystWidths_pO = np.repeat(0.2, len(ptWidths_pO))
    jpsiRaaVsPt_pO = np.array([0.719, 0.779, 0.797, 0.857, 0.864, 0.996])
    jpsiStatRaaVsPt_pO = np.array([0.0311327, 0.02479557, 0.02758417, 0.03219749, 0.03030912, 0.06258864])
    jpsiSystRaaVsPt_pO = np.array([0.04935935, 0.05198267, 0.04952558, 0.05080296, 0.04625856, 0.060258])


    graStatRaaVsPtComparison_pO = ROOT.TGraphErrors(len(ptCenters_pO), ptCenters_pO, jpsiRaaVsPt_pO, ptWidths_pO, jpsiStatRaaVsPt_pO)
    graSystRaaVsPtComparison_pO = ROOT.TGraphErrors(len(ptCenters_pO), ptCenters_pO, jpsiRaaVsPt_pO, np.array(ptSystWidths_pO), jpsiSystRaaVsPt_pO)
    
    SetGraStat(graStatRaaVsPtComparison_pO, 24, ROOT.kBlack)
    SetGraSyst(graSystRaaVsPtComparison_pO, 20, ROOT.kBlue+1,1,2)

    canvasRaaVsPtComparison_pO = ROOT.TCanvas("canvasRaaVsPtComparison_pO", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)

    histGridRaaVsPtComparison_pO = ROOT.TH2D("histGridRaaVsPtComparison_pO", 
                                          ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}", 
                                          100, 0, 8, 100, 0, 1.5)
    histGridRaaVsPtComparison_pO.Draw()

    # --- Unity line ---
    lineUnityComparison_pO = ROOT.TLine(0., 1., 8., 1.)
    lineUnityComparison_pO.SetLineColor(ROOT.kGray+1)
    lineUnityComparison_pO.SetLineWidth(2)
    lineUnityComparison_pO.SetLineStyle(ROOT.kDashed)
    lineUnityComparison_pO.Draw("SAME")

    graStatRaaVsPt_MethodSigma.Draw("EP SAME")
    graSystRaaVsPt_MethodSigma.Draw("E2P SAME")

    graStatRaaVsPtComparison_pO.Draw("EP SAME")
    graSystRaaVsPtComparison_pO.Draw("E2P SAME")

    legendRaaVsPtComparison_pO = ROOT.TLegend(0.20, 0.70, 0.45, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsPtComparison_pO)
    legendRaaVsPtComparison_pO.SetTextSize(0.045)
    legendRaaVsPtComparison_pO.AddEntry(graSystRaaVsPt_MethodSigma, "#sqrt{#it{s}} = 5.36 TeV, OO,  2.5 < #it{y} < 4", "FP")
    legendRaaVsPtComparison_pO.AddEntry(graSystRaaVsPtComparison_pO, "#sqrt{#it{s}} = 9.62 TeV, pO,  2.15 < #it{y} < 3.65", "FP")
    legendRaaVsPtComparison_pO.Draw("SAME")

    latexTitle.DrawLatex(0.22, 0.87, "ALICE Preliminary, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}")

    xMinBox_pO = 0.0
    xMaxBox_pO = 0.2
    xMinBox_OO = 0.2
    xMaxBox_OO = 0.4

    systCommonTot_pO = 0.067
    yMinBox_pO = yRef * (1.0 - systCommonTot_pO)
    yMaxBox_pO = yRef * (1.0 + systCommonTot_pO)

    boxCommon_OO = ROOT.TBox(xMinBox_OO, yMinBox_MethodSigma, xMaxBox_OO,yMaxBox_MethodSigma )
    boxCommon_OO.SetFillColorAlpha(ROOT.kRed+1, 1)
    boxCommon_OO.SetLineColor(ROOT.kRed+1)
    boxCommon_OO.SetLineWidth(1)

    boxCommon_pO = ROOT.TBox(xMinBox_pO, yMinBox_pO, xMaxBox_pO, yMaxBox_pO)
    boxCommon_pO.SetFillColorAlpha(ROOT.kBlue+1, 1)
    boxCommon_pO.SetLineColor(ROOT.kBlue+1)
    boxCommon_pO.SetLineWidth(1)

    boxCommon_OO.Draw("EP same")
    boxCommon_pO.Draw("EP same")

    canvasRaaVsPtComparison_pO.Update()
    #canvasRaaVsPtComparison_pO.SaveAs("SQM2026/Jpsi_RAA_vs_Pt_Comparison_with_pO.pdf")
    #canvasRaaVsPtComparison_pO.SaveAs("SQM2026/Jpsi_RAA_vs_Pt_Comparison_with_pO.png")

    input()


if __name__ == '__main__':
    main()