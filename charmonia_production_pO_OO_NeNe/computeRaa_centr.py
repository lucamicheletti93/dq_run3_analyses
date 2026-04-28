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
    TaaVsCentr = config["inputs"]["TaaVsCentr"]
    Ncoll = config["inputs"]["Ncoll"]
    errNcoll = config["inputs"]["errNcoll"]
    Npart = config["inputs"]["Npart"]
    errNpart = config["inputs"]["errNpart"]

    print("***** Compute corrected yield *****")
    # Lumi computed with normalization.C  + check_normalization
    fInLumi = ROOT.TFile(config["inputs"]["fInLumi"], "READ")

    # Normalization: N_minbias
    #histNevMinBiasCentr = fInLumi.Get("histnEvtsBcSelCentr_AfterCuts_TMaker") #  AfterCuts, YES PileUpCorrection
    #histNevMinBiasCentr = fInLumi.Get("histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp") #  AfterCuts, NO PileUpCorrection
    histNevMinBiasCentr = fInLumi.Get("histSelectedIntegrals_CentrFT0C_AfterCuts_TReader") #  Only Ncoll, NO others (as no PileUp Correction)
    nCentrBins = histNevMinBiasCentr.GetNbinsX()-6
    centrLabels = ["0-10%","10-20%", "20-30%", "30-40%", "40-50%", "50-60%"]
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
    fInPPrefVsCentr = ROOT.TFile(config["inputs"]["fInPPredVsCentr"], "READ")
    fInPPrefVsCentr.ls()
    histStatPPrefXsecVsSqrts = fInPPrefVsCentr.Get("histStatJpsiXsecInterpVsSqrts") # pT integrated
    histSystPPrefXsecVsSqrts = fInPPrefVsCentr.Get("histSystJpsiXsecInterpVsSqrts") # pT integrated
    print(f"histSystPPrefXsecVsSqrts: {histSystPPrefXsecVsSqrts}")

    SetHistStat(histStatPPrefXsecVsSqrts, 20, ROOT.kAzure+4)
    SetHistSyst(histSystPPrefXsecVsSqrts, 20, ROOT.kAzure+4)

    canvasPPrefXsecVsSqrts = ROOT.TCanvas("canvasPPrefXsecVsSqrts", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histStatPPrefXsecVsSqrts.Draw("EP")
    histSystPPrefXsecVsSqrts.Draw("E2P SAME")
    canvasPPrefXsecVsSqrts.Update()

    jpsiPPrefXsecVsCentr = histStatPPrefXsecVsSqrts.GetBinContent(1)
    print(f"jpsiPPrefXsecVsCentr: {jpsiPPrefXsecVsCentr}")
    jpsiStatPPrefXsecVsCentr = histStatPPrefXsecVsSqrts.GetBinError(1)
    jpsiSystPPrefXsecVsCentr = histSystPPrefXsecVsSqrts.GetBinError(1)
    print(f"jpsiSystPPrefXsecVsCentr: {jpsiSystPPrefXsecVsCentr}")


    print("-------- Centrality dependence --------")
    dfJpsiRawYieldVsCentr = pd.read_csv(config["inputs"]["fInRawYieldVsCentr"], sep=' ')
    centrMin = dfJpsiRawYieldVsCentr["x_min"]
    #print("centrMin:")
    #print(centrMin)
    #print("type:", type(centrMin))
    #print("len:", len(centrMin))
    centrMax = dfJpsiRawYieldVsCentr["x_max"]
    centrCenters = (centrMax + centrMin) / 2.
    centrWidths = (centrMax - centrMin) / 2.
    centrSystWidths = np.repeat(1.5, len(centrWidths))
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
    histSystRawYieldVsCentr.Draw("E2P SAME")
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

    histStatAxeVsCentr = ROOT.TH1D("histStatAxeVsCentr", ";Centrality (%);A#times#varepsilon", len(centrEdges)-1, centrEdges)

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

    print("Relative systematic on raw yield extraction")
    jpsiSystRelRawYieldVsCentr = (np.array(jpsiStatRawYieldVsCentr) / np.array(jpsiRawYieldVsCentr))
    print(jpsiSystRelRawYieldVsCentr)

    print("Relative systematic on pp reference")
    jpsi1SystRelPPrefXsecVsCentr = (np.array(jpsiStatPPrefXsecVsCentr) / np.array(jpsiPPrefXsecVsCentr))
    jpsi2SystRelPPrefXsecVsCentr = (np.array(jpsiSystPPrefXsecVsCentr) / np.array(jpsiPPrefXsecVsCentr))
    jpsiSystRelPPrefXsecVsCentr = np.sqrt(jpsi1SystRelPPrefXsecVsCentr**2 + jpsi2SystRelPPrefXsecVsCentr**2)
    print(jpsiSystRelPPrefXsecVsCentr)


    #jpsiRaaVsCentr = (jpsiRawYieldVsCentr) / (jpsiAxeVsCentr * BrJpsiToMuMu * (nevCentrArray/0.935) * np.array(TaaVsCentr) * jpsiPPrefXsecVsCentr)
    jpsiRaaVsCentr = (jpsiRawYieldVsCentr) / (jpsiAxeVsCentr * BrJpsiToMuMu * (nevCentrArray) * np.array(TaaVsCentr) * jpsiPPrefXsecVsCentr)
    jpsiStatRaaVsCentr = jpsiRaaVsCentr * (np.array(jpsiStatRawYieldVsCentr) / np.array(jpsiRawYieldVsCentr))
    jpsiSystRaaVsCentr = jpsiRaaVsCentr * np.sqrt(jpsiSystRelRawYieldVsCentr**2 + jpsiSystRelPPrefXsecVsCentr**2)

    graStatRaaVsCentr = ROOT.TGraphErrors(len(centrCenters), np.array(centrCenters), np.array(jpsiRaaVsCentr), np.array(centrWidths), np.array(jpsiStatRaaVsCentr))
    graSystRaaVsCentr = ROOT.TGraphErrors(len(centrCenters), np.array(centrCenters), np.array(jpsiRaaVsCentr), np.array(centrSystWidths), np.array(jpsiSystRaaVsCentr))

    SetGraStat(graStatRaaVsCentr, 24, ROOT.kBlack)
    SetGraSyst(graSystRaaVsCentr, 20, ROOT.kBlue, 1, 2)

    lineUnityVsCentr = ROOT.TLine(0., 1., 90., 1.)
    lineUnityVsCentr.SetLineColor(ROOT.kGray+1)
    lineUnityVsCentr.SetLineWidth(2)
    lineUnityVsCentr.SetLineStyle(ROOT.kDashed) 

    canvasRaaVsCentrVsRun2 = ROOT.TCanvas("canvasRaaVsCentrVsRun2", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    histGridRaaVsCentr = ROOT.TH2D("histGridRaaVsCentr", ";Centrality (%);#it{R}_{AA}", 100, 0, 90, 100, 0, 1.8)
    histGridRaaVsCentr.Draw()
    lineUnityVsCentr.Draw("SAME")
    graStatRaaVsCentr.Draw("EP SAME")
    graSystRaaVsCentr.Draw("E2P SAME")

    legendRaaVsCentrVsRun2 = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsCentrVsRun2)
    legendRaaVsCentrVsRun2.SetTextSize(0.045)
    legendRaaVsCentrVsRun2.AddEntry(graSystRaaVsCentr, "#sqrt{#it{s}} = 5.36 TeV, OO, 0 < #it{p}_{T} (GeV/#it{c}) < 20 ", "FP")
    legendRaaVsCentrVsRun2.Draw("SAME")

    latexTitle.DrawLatex(0.20, 0.90, "ALICE Work in Progress")
    latexTitle.DrawLatex(0.20, 0.84, "J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    canvasRaaVsCentrVsRun2.Update()
    #canvasRaaVsCentrVsRun2.SaveAs("Study_centr_ROO/Jpsi_ROO_vs_Centr.pdf")
    #canvasRaaVsCentrVsRun2.SaveAs("Study_centr_ROO/Jpsi_ROO_vs_Centr.png")

    print("-------- RAA vs Ncoll --------") 
    centrSystWidths_NColl = np.repeat(0.4, len(centrWidths))
    graStatRaaVsNcoll = ROOT.TGraphErrors(len(centrMin),
                                          np.array(Ncoll, dtype=float),
                                          np.array(jpsiRaaVsCentr, dtype=float),
                                          np.array(errNcoll),
                                          np.array(jpsiStatRaaVsCentr, dtype=float))

    graSystRaaVsNcoll = ROOT.TGraphErrors(len(centrMin),
                                          np.array(Ncoll, dtype=float),
                                          np.array(jpsiRaaVsCentr, dtype=float),
                                          np.array(centrSystWidths_NColl),
                                          np.array(jpsiSystRaaVsCentr, dtype=float))

    SetGraStat(graStatRaaVsNcoll, 24, ROOT.kBlack)
    SetGraSyst(graSystRaaVsNcoll, 20, ROOT.kRed+1, 1, 2)

    lineUnityVsNcoll = ROOT.TLine(0., 1., max(Ncoll)*1.1, 1.)
    lineUnityVsNcoll.SetLineColor(ROOT.kGray+1)
    lineUnityVsNcoll.SetLineWidth(2)
    lineUnityVsNcoll.SetLineStyle(ROOT.kDashed)

    canvasRaaVsNcoll = ROOT.TCanvas("canvasRaaVsNcoll", "", 800, 600)
    histGridRaaVsNcoll = ROOT.TH2D("histGridRaaVsNcoll",";#it{N}_{coll};R_{AA}",100, 0, max(Ncoll)*1.1,100, 0, 1.8)
    histGridRaaVsNcoll.Draw()
    lineUnityVsNcoll.Draw("SAME")
    graSystRaaVsNcoll.Draw("E2P SAME")
    graStatRaaVsNcoll.Draw("EP SAME")

    legendRaaVsNcoll = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsNcoll)
    legendRaaVsNcoll.SetTextSize(0.045)
    legendRaaVsNcoll.AddEntry(graSystRaaVsNcoll, "#sqrt{#it{s}} = 5.36 TeV, OO, 0 < #it{p}_{T} (GeV/#it{c}) < 20 ", "FP")
    legendRaaVsNcoll.Draw("SAME")

    latexTitle.DrawLatex(0.20, 0.90, "ALICE Work In Progress")
    latexTitle.DrawLatex(0.20, 0.84, "J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    #canvasRaaVsNcoll.SaveAs("Study_centr_ROO/Jpsi_ROO_vs_Ncoll.pdf")
    #canvasRaaVsNcoll.SaveAs("Study_centr_ROO/Jpsi_ROO_vs_Ncoll.png")

    # Compare to Run 2 results, pPb
    print("***** Extract RAA from Run 2, pPb and Pbp *****")

    # --- p-Pb: pT 0–20, 2.03 < y < 3.53  ---
    filePathJpsi_020pT_pPb_pos_y = "HEP_Data/jpsi_raa_pPb_8TeV_0_20_pT_pos_y.yaml"
    ptCenters_020pT_pPb_pos_y, ptWidths_020pT_pPb_pos_y, jpsiRaaVsCentr_020pT_pPb_pos_y, jpsiStatRaaVsCentr_020pT_pPb_pos_y, jpsiSystRaaVsCentr_020pT_pPb_pos_y = ExtractFromYaml(filePathJpsi_020pT_pPb_pos_y)
    jpsiPtSystWidths_020pT_pPb_pos_y = np.repeat(0.2, len(ptCenters_020pT_pPb_pos_y))
    
    graStatRaaVsNColl_020pT_pPb_pos_y = ROOT.TGraphErrors(len(ptCenters_020pT_pPb_pos_y), np.array(ptCenters_020pT_pPb_pos_y), np.array(jpsiRaaVsCentr_020pT_pPb_pos_y),
                                    np.array(ptWidths_020pT_pPb_pos_y), np.array(jpsiStatRaaVsCentr_020pT_pPb_pos_y))
    graSystRaaVsNColl_020pT_pPb_pos_y = ROOT.TGraphErrors(len(ptCenters_020pT_pPb_pos_y), np.array(ptCenters_020pT_pPb_pos_y), np.array(jpsiRaaVsCentr_020pT_pPb_pos_y),
                                    jpsiPtSystWidths_020pT_pPb_pos_y, np.array(jpsiSystRaaVsCentr_020pT_pPb_pos_y))
    SetGraStat(graStatRaaVsNColl_020pT_pPb_pos_y, 20, ROOT.kBlack)
    SetGraSyst(graSystRaaVsNColl_020pT_pPb_pos_y, 20, ROOT.kBlack)

    # --- Pb-p: pT 0–20, -4.46 < y < -2.96  ---
    filePathJpsi_020pT_pPb_neg_y = "HEP_Data/jpsi_raa_pPb_8TeV_0_20_pT_neg_y.yaml"
    ptCenters_020pT_pPb_neg_y, ptWidths_020pT_pPb_neg_y, jpsiRaaVsCentr_020pT_pPb_neg_y, jpsiStatRaaVsCentr_020pT_pPb_neg_y, jpsiSystRaaVsCentr_020pT_pPb_neg_y = ExtractFromYaml(filePathJpsi_020pT_pPb_neg_y)
    jpsiPtSystWidths_020pT_pPb_neg_y = np.repeat(0.2, len(ptCenters_020pT_pPb_neg_y))
    
    graStatRaaVsNColl_020pT_pPb_neg_y = ROOT.TGraphErrors(len(ptCenters_020pT_pPb_neg_y), np.array(ptCenters_020pT_pPb_neg_y), np.array(jpsiRaaVsCentr_020pT_pPb_neg_y),
                                    np.array(ptWidths_020pT_pPb_neg_y), np.array(jpsiStatRaaVsCentr_020pT_pPb_neg_y))
    graSystRaaVsNColl_020pT_pPb_neg_y = ROOT.TGraphErrors(len(ptCenters_020pT_pPb_neg_y), np.array(ptCenters_020pT_pPb_neg_y), np.array(jpsiRaaVsCentr_020pT_pPb_neg_y),
                                    jpsiPtSystWidths_020pT_pPb_neg_y, np.array(jpsiSystRaaVsCentr_020pT_pPb_neg_y))
    SetGraStat(graStatRaaVsNColl_020pT_pPb_neg_y, 20, ROOT.kBlue+1)
    SetGraSyst(graSystRaaVsNColl_020pT_pPb_neg_y, 20, ROOT.kBlue+1)
    
    # --- Canvas setup ---
    canvasRaaVsNcoll_Comparison_pPb = ROOT.TCanvas("canvasRaaVsNcoll_Comparison_pPb", "", 800, 600)
    canvasRaaVsNcoll_Comparison_pPb.SetLogx()
    ROOT.gStyle.SetOptStat(False)
    
    histGridRaaVsNColl_Comparison_pPb = ROOT.TH2D("histGridRaaVsNColl_Comparison_pPb",
                                ";#it{N}_{coll};#it{R}_{AA}",
                                100, 1, max(Ncoll)*1.1, 100, 0, 1.8)
    histGridRaaVsNColl_Comparison_pPb.Draw()
    
    # --- Unity line ---
    lineUnityVsNColl_Comparison_pPb = ROOT.TLine(0., 1., max(Ncoll)*1.1, 1.)
    lineUnityVsNColl_Comparison_pPb.SetLineColor(ROOT.kGray+1)
    lineUnityVsNColl_Comparison_pPb.SetLineWidth(2)
    lineUnityVsNColl_Comparison_pPb.SetLineStyle(ROOT.kDashed)
    lineUnityVsNColl_Comparison_pPb.Draw("SAME")
    
    # --- Draw graphs ---
    graSystRaaVsNColl_020pT_pPb_pos_y.Draw("2 SAME")
    graStatRaaVsNColl_020pT_pPb_pos_y.Draw("P SAME")
    graSystRaaVsNColl_020pT_pPb_neg_y.Draw("2 SAME")
    graStatRaaVsNColl_020pT_pPb_neg_y.Draw("P SAME")
    graStatRaaVsNcoll.Draw("EP SAME")
    graSystRaaVsNcoll.Draw("E2P SAME")

    legendRaaVsNColl_Comparison_pPb = ROOT.TLegend(0.20, 0.63, 0.40, 0.89, " ", "brNDC")
    SetLegend(legendRaaVsNColl_Comparison_pPb)
    legendRaaVsNColl_Comparison_pPb.SetTextSize(0.045)
    legendRaaVsNColl_Comparison_pPb.AddEntry(graSystRaaVsNcoll, "#sqrt{#it{s}} = 5.36 TeV, OO, 2.5 < #it{y} < 4", "FP")
    legendRaaVsNColl_Comparison_pPb.AddEntry(graStatRaaVsNColl_020pT_pPb_pos_y, "#sqrt{#it{s}} = 8.16 TeV, p-Pb, 2.03  < #it{y} < 3.53", "FP")
    legendRaaVsNColl_Comparison_pPb.AddEntry(graStatRaaVsNColl_020pT_pPb_neg_y, "#sqrt{#it{s}} = 8.16 TeV, p-Pb, -4.46  < #it{y} < -2.06", "FP")
    legendRaaVsNColl_Comparison_pPb.Draw("SAME")

    latexTitle.DrawLatex(0.20, 0.90, "ALICE Work In Progress")
    latexTitle.DrawLatex(0.20, 0.84, "J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 0 < #it{p}_{T} (GeV/#it{c}) < 20")
    canvasRaaVsNcoll_Comparison_pPb.Update()
    #canvasRaaVsNcoll_Comparison_pPb.SaveAs("Study_centr_ROO/Jpsi_ROO_vs_Ncoll_Comparison_pPb.pdf")
    #canvasRaaVsNcoll_Comparison_pPb.SaveAs("Study_centr_ROO/Jpsi_ROO_vs_Ncoll_Comparison_pPb.png")

    print("-------- RAA vs Npart --------") 
    centrSystWidths_Npart = np.repeat(0.5, len(centrWidths))
    graStatRaaVsNpart = ROOT.TGraphErrors(len(centrMin),
                                          np.array(Npart, dtype=float),
                                          np.array(jpsiRaaVsCentr, dtype=float),
                                          np.array(errNpart),
                                          np.array(jpsiStatRaaVsCentr, dtype=float))

    graSystRaaVsNpart = ROOT.TGraphErrors(len(centrMin),
                                          np.array(Npart, dtype=float),
                                          np.array(jpsiRaaVsCentr, dtype=float),
                                          np.array(centrSystWidths_Npart),
                                          np.array(jpsiSystRaaVsCentr, dtype=float))

    SetGraStat(graStatRaaVsNpart, 24, ROOT.kBlack)
    SetGraSyst(graSystRaaVsNpart, 20, ROOT.kMagenta, 1, 2)

    lineUnityVsNpart = ROOT.TLine(5., 1., max(Npart)*1.1, 1.)
    lineUnityVsNpart.SetLineColor(ROOT.kGray+1)
    lineUnityVsNpart.SetLineWidth(2)
    lineUnityVsNpart.SetLineStyle(ROOT.kDashed)

    canvasRaaVsNpart = ROOT.TCanvas("canvasRaaVsNpart", "", 800, 600)
    histGridRaaVsNpart = ROOT.TH2D("histGridRaaVsNcoll",";#it{N}_{part};R_{AA}",100, 5, max(Npart)*1.1,100, 0, 1.8)
    histGridRaaVsNpart.Draw()
    lineUnityVsNpart.Draw("SAME")
    graSystRaaVsNpart.Draw("E2P SAME")
    graStatRaaVsNpart.Draw("EP SAME")
    graStatRaaVsNpart.Fit("pol0")
    f = graStatRaaVsNpart.GetFunction("pol0")
    f.SetLineColor(ROOT.kBlack)

    legendRaaVsNpart = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsNpart)
    legendRaaVsNpart.SetTextSize(0.045)

    latexTitle.DrawLatex(0.20, 0.84, "J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    legendRaaVsNpart.AddEntry(graSystRaaVsNpart, "#sqrt{#it{s}} = 5.36 TeV, OO, 0 < #it{p}_{T} (GeV/#it{c}) < 20 ", "FP")
    legendRaaVsNpart.Draw("SAME")

    output = np.column_stack((centrWidths, jpsiRaaVsCentr, jpsiStatRaaVsCentr, jpsiSystRaaVsCentr))

    #np.savetxt("Study_centr_ROO/Train_617522/NEW/jpsi_RAA_vs_centr_Train_617522_PbPbQuality_ONLYNcoll_NOPileUpCorrection_EFF1.txt", output, header="pT_center(GeV/c)   RAA   stat_err   syst_err", fmt="%.3f   %.6f   %.6f   %.6f")
    
    # --- Compare with Previous Results: Method 1 Here and Method 1 in Previous Results ---
    '''
    jpsiRaaVsCentr_Train_591157_Sel8 = np.array([0.738633,0.718966,0.714172,0.795502,0.87903,0.968229,0.982391,0.811491])
    jpsiStatRaaVsCentr_Train_591157_Sel8 = np.array([0.015571,0.016523,0.01522,0.01854,0.025909,0.030025,0.03115,0.031349])
    jpsiSystRaaVsCentr_Train_591157_Sel8 = np.array([0.037341,0.036937,0.036173,0.040985,0.047986,0.053673,0.054849,0.048714])
    # 617522 ma normal n_mb 
    jpsiRaaVsCentr_Train_591157_Sel8 = np.array([0.722434,0.691527,0.698491,0.759254,0.835379,0.916794]) qui usual
    jpsiStatRaaVsCentr_Train_591157_Sel8 = np.array([0.013832,0.01369,0.014802,0.016596,0.025796,0.025291])
    jpsiSystRaaVsCentr_Train_591157_Sel8 = np.array([0.035962,0.034598,0.035344,0.038633,0.046247,0.049134])

    jpsiRaaVsCentr_Train_591157_Sel8 = np.array([0.807427,0.772883,0.780666,0.848578,0.933659,1.024652]) # qui eff = 1 e no corr PileUp
    jpsiStatRaaVsCentr_Train_591157_Sel8 = np.array([0.01546,0.0153,0.016543,0.018549,0.028831,0.028266])
    jpsiSystRaaVsCentr_Train_591157_Sel8 = np.array([0.040192,0.038669,0.039502,0.043178,0.051688,0.054915])

    # 617522 Eff = 93.5%, YES PileUp Correction NEW
    jpsiRaaVsCentr_Train_591157_Sel8 = np.array([0.722434,0.691527,0.698491,0.759254,0.835379,0.916794])
    jpsiStatRaaVsCentr_Train_591157_Sel8 = np.array([0.013832,0.01369,0.014802,0.016596,0.025796,0.025291])
    jpsiSystRaaVsCentr_Train_591157_Sel8 = np.array([0.035962,0.034598,0.035344,0.038633,0.046247,0.049134])

    # 617522 Eff = 93.5%, NO PileUp Correction NEW
    jpsiRaaVsCentr_Train_591157_Sel8 = np.array([0.754944,0.722646,0.729923,0.79342,0.872971,0.95805])
    jpsiStatRaaVsCentr_Train_591157_Sel8 = np.array([0.014455,0.014306,0.015468,0.017343,0.026957,0.026429])
    jpsiSystRaaVsCentr_Train_591157_Sel8 = np.array([0.03758,0.036155,0.036934,0.040372,0.048328,0.051345])

    # 617522 Eff = 1, NO PileUp Correction NEW
    jpsiRaaVsCentr_Train_591157_Sel8 = np.array([0.807427,0.772883,0.780666,0.848578,0.933659,1.024652])
    jpsiStatRaaVsCentr_Train_591157_Sel8 = np.array([0.01546,0.0153,0.016543,0.018549,0.028831,0.028266])
    jpsiSystRaaVsCentr_Train_591157_Sel8 = np.array([0.040192,0.038669,0.039502,0.043178,0.051688,0.054915])

    graStatRaaVsCentr_Train_591157_Sel8 = ROOT.TGraphErrors(len(Npart), np.array(Npart), jpsiRaaVsCentr_Train_591157_Sel8, np.array(errNpart), jpsiStatRaaVsCentr_Train_591157_Sel8)
    graSystRaaVsCentr_Train_591157_Sel8 = ROOT.TGraphErrors(len(Npart), np.array(Npart), jpsiRaaVsCentr_Train_591157_Sel8, np.array(centrSystWidths_Npart), jpsiSystRaaVsCentr_Train_591157_Sel8)
    SetGraStat(graStatRaaVsCentr_Train_591157_Sel8, 24, ROOT.kBlack)
    SetGraSyst(graSystRaaVsCentr_Train_591157_Sel8, 20, ROOT.kBlue,1,2)

    canvasCompare__ = ROOT.TCanvas("canvasCompare__", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    histGrid__ = ROOT.TH2D("histGrid",";#it{N}_{part};#it{R}_{AA}", 100, 5,  max(Npart)*1.1, 100, 0, 1.8)
    histGrid__.Draw()
    lineUnityVsNpart.Draw("SAME")
    #boxCommon_Method1.Draw("same")
    graSystRaaVsNpart.Draw("E2P SAME")
    graStatRaaVsNpart.Draw("EP SAME")
    graStatRaaVsCentr_Train_591157_Sel8.Draw("EP SAME")
    graSystRaaVsCentr_Train_591157_Sel8.Draw("E2P SAME")

    legend_ = ROOT.TLegend(0.2, 0.7, 0.5, 0.87, "", "brNDC")
    SetLegend(legend_)
    legend_.AddEntry(graSystRaaVsCentr_Train_591157_Sel8, "PbPbQuality, No PileUp Corr, Eff = 1", "FP")
    legend_.AddEntry(graSystRaaVsNpart, "PbPbQuality, NO PileUp Corr, Eff = 1, (only NcollSel+PU)", "FP")
    
    legend_.Draw("SAME")


    canvasCompare__.Update()
    canvasCompare__.SaveAs("Study_centr_ROO/Jpsi_RAA_vs_Centr_Method1_Comparison_here_and_previous_results_617522_Nmb.pdf")
    '''


    # -------------------------------
    # Common systematics (global box)
    # -------------------------------
    systmatchingeffJpsi_centr = 0.02   
    systtrackeffJpsi_centr = 0.015
    systBR_centr = 0.005
    systTOO_centr = 0.06
    systEffTVX_centr = 0.06

    systCommonTot_Method1_centr = np.sqrt(systmatchingeffJpsi_centr**2 + systtrackeffJpsi_centr**2 + systBR_centr**2 + systTOO_centr**2 + systEffTVX_centr**2)

    xMinBox_centr = 27.0
    xMaxBox_centr = 27.7
    yRef_centr = 1.0

    yMinBox_centr = yRef_centr * (1.0 - systCommonTot_Method1_centr)
    yMaxBox_centr = yRef_centr * (1.0 + systCommonTot_Method1_centr)

    boxCommon_centr = ROOT.TBox(xMinBox_centr, yMinBox_centr, xMaxBox_centr, yMaxBox_centr)
    boxCommon_centr.SetFillColorAlpha(ROOT.kMagenta, 0.4)
    boxCommon_centr.SetLineColor(ROOT.kMagenta)
    boxCommon_centr.SetLineWidth(1)
    boxCommon_centr.Draw("same")

    latexTitle.DrawLatex(0.20, 0.90, "ALICE Work In Progress")
    latexTitle.DrawLatex(0.20, 0.84, "J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    #canvasRaaVsNpart.SaveAs("Study_centr_ROO/Jpsi_ROO_vs_Npart.pdf")
    #canvasRaaVsNpart.SaveAs("Study_centr_ROO/Jpsi_ROO_vs_Npart.png")



    # Compare to Run 2 results, PbPb
    print("***** Extract RAA from Run 2, PbPb *****")

    # --- PbPb: pT 0.3–8 GeV/c ---
    filePathJpsi_0_8pT_PbPb = "HEP_Data/jpsi_raa_PbPb_5TeV_0_8_pT.yaml"
    ptCenters_0_8pT_PbPb, ptWidths_0_8pT_PbPb, jpsiRaaVsNpart_0_8pT_PbPb, jpsiStatRaaVsNpart_0_8pT_PbPb, jpsiSystRaaVsNpart_0_8pT_PbPb = ExtractFromYaml(filePathJpsi_0_8pT_PbPb)
    jpsiPtSystWidths_0_8pT_PbPb = np.repeat(1.2, len(ptCenters_0_8pT_PbPb))
    
    graStatRaaVsNpart_0_8pT_PbPb = ROOT.TGraphErrors(len(ptCenters_0_8pT_PbPb), np.array(ptCenters_0_8pT_PbPb), np.array(jpsiRaaVsNpart_0_8pT_PbPb),
                                    np.array(ptWidths_0_8pT_PbPb), np.array(jpsiStatRaaVsNpart_0_8pT_PbPb))
    graSystRaaVsNpart_0_8pT_PbPb = ROOT.TGraphErrors(len(ptCenters_0_8pT_PbPb), np.array(ptCenters_0_8pT_PbPb), np.array(jpsiRaaVsNpart_0_8pT_PbPb),
                                    jpsiPtSystWidths_0_8pT_PbPb, np.array(jpsiSystRaaVsNpart_0_8pT_PbPb))
    SetGraStat(graStatRaaVsNpart_0_8pT_PbPb, 20, ROOT.kGray+2)
    SetGraSyst(graSystRaaVsNpart_0_8pT_PbPb, 20, ROOT.kGray+2)
    
    # --- Canvas setup ---
    canvasRaaVsNpart_Comparison_PbPb = ROOT.TCanvas("canvasRaaVsNpart_Comparison_PbPb", "", 800, 600)
    canvasRaaVsNpart_Comparison_PbPb.SetLogx()
    ROOT.gStyle.SetOptStat(False)
    
    histGridRaaVsNpart_Comparison_PbPb = ROOT.TH2D("histGridRaaVsNpart_Comparison_PbPb",
                                ";#it{N}_{part};#it{R}_{AA}",
                                100, 0, 400, 100, 0, 1.8)
    histGridRaaVsNpart_Comparison_PbPb.Draw()
    
    # --- Unity line ---
    lineUnityVsNpart_Comparison_PbPb = ROOT.TLine(0., 1., 400, 1.)
    lineUnityVsNpart_Comparison_PbPb.SetLineColor(ROOT.kGray+1)
    lineUnityVsNpart_Comparison_PbPb.SetLineWidth(2)
    lineUnityVsNpart_Comparison_PbPb.SetLineStyle(ROOT.kDashed)
    lineUnityVsNpart_Comparison_PbPb.Draw("SAME")
    
    # --- Draw graphs ---
    graSystRaaVsNpart_0_8pT_PbPb.Draw("2 SAME")
    graStatRaaVsNpart_0_8pT_PbPb.Draw("P SAME")
    graSystRaaVsNpart.Draw("E2P SAME")
    graStatRaaVsNpart.Draw("EP SAME")

    legendRaaVsNpart_Comparison_PbPb = ROOT.TLegend(0.20, 0.65, 0.40, 0.85, " ", "brNDC")
    SetLegend(legendRaaVsNpart_Comparison_PbPb)
    legendRaaVsNpart_Comparison_PbPb.SetTextSize(0.045)
    legendRaaVsNpart_Comparison_PbPb.AddEntry(graSystRaaVsNpart, "#sqrt{#it{s}} = 5.36 TeV, OO, 0 < #it{p}_{T} < 20 GeV/#it{c}", "FP")
    legendRaaVsNpart_Comparison_PbPb.AddEntry(graStatRaaVsNpart_0_8pT_PbPb, "#sqrt{#it{s}} = 5.02 TeV, Pb#minusPb, 0.3 < #it{p}_{T} < 8 GeV/#it{c}", "FP")
    legendRaaVsNpart_Comparison_PbPb.Draw("SAME")

    latexTitle.DrawLatex(0.20, 0.90, "ALICE Work In Progress")
    latexTitle.DrawLatex(0.20, 0.84, "J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    canvasRaaVsNpart_Comparison_PbPb.Update()
    #canvasRaaVsNpart_Comparison_PbPb.SaveAs("Study_centr_ROO/Jpsi_ROO_vs_Npart_Comparison_PbPb.pdf")
    #canvasRaaVsNpart_Comparison_PbPb.SaveAs("Study_centr_ROO/Jpsi_ROO_vs_Npart_Comparison_PbPb.png")


    input()



if __name__ == '__main__':
    main()