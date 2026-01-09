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
    Ncoll = config["inputs"]["Ncoll"]
    errNcoll = config["inputs"]["errNcoll"]
    Npart = config["inputs"]["Npart"]
    errNpart = config["inputs"]["errNpart"]

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
    #histNevMinBiasCentr = fInLumi.Get("histNeventMinBias_Centr") # old
    #histNevMinBiasCentr = fInLumi.Get("histnEvtsBcSelCentr_BeforeCuts") # Method 1 - BeforeCuts
    #histNevMinBiasCentr = fInLumi.Get("histnEvtsBeforeCutsCentr_CollTVXdivCollAll") # Method 2 - BeforeCuts
    histNevMinBiasCentr = fInLumi.Get("histnEvtsBcSelCentr_AfterCuts") # Method 1 - AfterCuts
    #histNevMinBiasCentr = fInLumi.Get("histnEvtsAfterCutsCentr_CollTVXdivCollAll") # Method 2 - AfterCuts
    nCentrBins = histNevMinBiasCentr.GetNbinsX()
    centrLabels = ["0-5%","5-10%","10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%","70-90%", "0-90%", "0-100%"]
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


    jpsiRaaVsCentr = (jpsiRawYieldVsCentr) / (jpsiAxeVsCentr * BrJpsiToMuMu * nevCentrArray * np.array(TaaVsCentr) * jpsiPPrefXsecVsCentr)
    jpsiStatRaaVsCentr = jpsiRaaVsCentr * (np.array(jpsiStatRawYieldVsCentr) / np.array(jpsiRawYieldVsCentr))
    jpsiSystRaaVsCentr = jpsiRaaVsCentr * np.sqrt(jpsiSystRelRawYieldVsCentr**2 + jpsiSystRelPPrefXsecVsCentr**2)

    graStatRaaVsCentr = ROOT.TGraphErrors(len(centrCenters), np.array(centrCenters), np.array(jpsiRaaVsCentr), np.array(centrWidths), np.array(jpsiStatRaaVsCentr))
    graSystRaaVsCentr = ROOT.TGraphErrors(len(centrCenters), np.array(centrCenters), np.array(jpsiRaaVsCentr), np.array(centrSystWidths), np.array(jpsiSystRaaVsCentr))

    SetGraStat(graStatRaaVsCentr, 24, ROOT.kBlack)
    SetGraSyst(graSystRaaVsCentr, 20, ROOT.kRed+1, 1, 2)

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
    canvasRaaVsCentrVsRun2.SaveAs("plots/Jpsi_RAA_vs_Centr.pdf")
    canvasRaaVsCentrVsRun2.SaveAs("plots/Jpsi_RAA_vs_Centr.png")

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
    canvasRaaVsNcoll.SaveAs("plots/Jpsi_RAA_vs_Ncoll.pdf")
    canvasRaaVsNcoll.SaveAs("plots/Jpsi_RAA_vs_Ncoll.png")

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
    canvasRaaVsNcoll_Comparison_pPb.SaveAs("plots/Jpsi_RAA_vs_Ncoll_Comparison_pPb.pdf")
    canvasRaaVsNcoll_Comparison_pPb.SaveAs("plots/Jpsi_RAA_vs_Ncoll_Comparison_pPb.png")

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
    SetGraSyst(graSystRaaVsNpart, 20, ROOT.kRed+1, 1, 2)

    lineUnityVsNpart = ROOT.TLine(0., 1., max(Npart)*1.1, 1.)
    lineUnityVsNpart.SetLineColor(ROOT.kGray+1)
    lineUnityVsNpart.SetLineWidth(2)
    lineUnityVsNpart.SetLineStyle(ROOT.kDashed)

    canvasRaaVsNpart = ROOT.TCanvas("canvasRaaVsNpart", "", 800, 600)
    histGridRaaVsNpart = ROOT.TH2D("histGridRaaVsNcoll",";#it{N}_{part};R_{AA}",100, 0, max(Npart)*1.1,100, 0, 1.8)
    histGridRaaVsNpart.Draw()
    lineUnityVsNpart.Draw("SAME")
    graSystRaaVsNpart.Draw("E2P SAME")
    graStatRaaVsNpart.Draw("EP SAME")

    legendRaaVsNpart = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsNpart)
    legendRaaVsNpart.SetTextSize(0.045)
    legendRaaVsNpart.AddEntry(graSystRaaVsNpart, "#sqrt{#it{s}} = 5.36 TeV, OO, 0 < #it{p}_{T} (GeV/#it{c}) < 20 ", "FP")
    legendRaaVsNpart.Draw("SAME")

    latexTitle.DrawLatex(0.20, 0.90, "ALICE Work In Progress")
    latexTitle.DrawLatex(0.20, 0.84, "J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    canvasRaaVsNpart.SaveAs("plots/Jpsi_RAA_vs_Npart.pdf")
    canvasRaaVsNpart.SaveAs("plots/Jpsi_RAA_vs_Npart.png")

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
    canvasRaaVsNpart_Comparison_PbPb.SaveAs("plots/Jpsi_RAA_vs_Npart_Comparison_PbPb.pdf")
    canvasRaaVsNpart_Comparison_PbPb.SaveAs("plots/Jpsi_RAA_vs_Npart_Comparison_PbPb.png")

    #input()
    #exit()

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

    SetGraStat(graStatXsecVsPt, 24, ROOT.kBlack)
    SetGraSyst(graSystXsecVsPt, 20, ROOT.kRed+1, 1, 2)

    SetGraStat(graStatRaaVsPt, 24, ROOT.kBlack)
    SetGraSyst(graSystRaaVsPt, 20, ROOT.kRed+1, 1, 2)


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
    legendRaaVsPtVsRun2.AddEntry(graSystRaaVsPt, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus90%", "FP")
    legendRaaVsPtVsRun2.Draw("SAME")

    latexTitle.DrawLatex(0.20, 0.90, "ALICE Work In Progress")
    latexTitle.DrawLatex(0.20, 0.84, "J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    
    canvasRaaVsPt.Update()
    canvasRaaVsPt.SaveAs("plots/Jpsi_RAA_vs_Pt.pdf")
    canvasRaaVsPt.SaveAs("plots/Jpsi_RAA_vs_Pt.png")


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
    SetGraStat(graStatRaaVsPt_020, 20, ROOT.kOrange+7)
    SetGraSyst(graSystRaaVsPt_020, 20, ROOT.kOrange+1)
    
    # --- Centrality 20–40% ---
    filePathJpsi_2040 = "HEP_Data/jpsi_raa_PbPb_5TeV_20_40_centr.yaml"
    ptCenters_2040, ptWidths_2040, jpsiRaaVsPt_2040, jpsiStatRaaVsPt_2040, jpsiSystRaaVsPt_2040 = ExtractFromYaml(filePathJpsi_2040)
    jpsiPtSystWidths_2040 = np.repeat(0.2, len(ptCenters_2040))
    
    graStatRaaVsPt_2040 = ROOT.TGraphErrors(len(ptCenters_2040), np.array(ptCenters_2040), np.array(jpsiRaaVsPt_2040),
                                     np.array(ptWidths_2040), np.array(jpsiStatRaaVsPt_2040))
    graSystRaaVsPt_2040 = ROOT.TGraphErrors(len(ptCenters_2040), np.array(ptCenters_2040), np.array(jpsiRaaVsPt_2040),
                                     jpsiPtSystWidths_2040, np.array(jpsiSystRaaVsPt_2040))
    SetGraStat(graStatRaaVsPt_2040, 20, ROOT.kGreen+2)
    SetGraSyst(graSystRaaVsPt_2040, 20, ROOT.kGreen+2)
    
    # --- Centrality 40–60% ---
    filePathJpsi_4060 = "HEP_Data/jpsi_raa_PbPb_5TeV_40_90_centr.yaml"
    ptCenters_4060, ptWidths_4060, jpsiRaaVsPt_4060, jpsiStatRaaVsPt_4060, jpsiSystRaaVsPt_4060 = ExtractFromYaml(filePathJpsi_4060)
    jpsiPtSystWidths_4060 = np.repeat(0.2, len(ptCenters_4060))
    
    graStatRaaVsPt_4060 = ROOT.TGraphErrors(len(ptCenters_4060), np.array(ptCenters_4060), np.array(jpsiRaaVsPt_4060),
                                     np.array(ptWidths_4060), np.array(jpsiStatRaaVsPt_4060))
    graSystRaaVsPt_4060 = ROOT.TGraphErrors(len(ptCenters_4060), np.array(ptCenters_4060), np.array(jpsiRaaVsPt_4060),
                                     jpsiPtSystWidths_4060, np.array(jpsiSystRaaVsPt_4060))
    SetGraStat(graStatRaaVsPt_4060, 20, ROOT.kAzure+4)
    SetGraSyst(graSystRaaVsPt_4060, 20, ROOT.kAzure+4)
    
    # --- Canvas setup ---
    canvasRaaVsPt_Comparison1 = ROOT.TCanvas("canvasRaaVsPt_Comparison1", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    
    histGridRaaVsPt_Comparison1 = ROOT.TH2D("histGridRaaV_Comparison1",
                                ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}",
                                100, 0, 12, 100, 0, 2)
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
    graSystRaaVsPt.Draw("E2P SAME")
    graStatRaaVsPt.Draw("EP SAME")

    legendRaaVsPtVsRun2_Comparison1 = ROOT.TLegend(0.20, 0.60, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsPtVsRun2_Comparison1)
    legendRaaVsPtVsRun2_Comparison1.SetTextSize(0.045)
    legendRaaVsPtVsRun2_Comparison1.AddEntry(graSystRaaVsPt, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus90%", "FP")
    legendRaaVsPtVsRun2_Comparison1.AddEntry(graStatRaaVsPt_020, "#sqrt{#it{s}} = 5.02 TeV, Pb-Pb, 0#minus20%", "FP")
    legendRaaVsPtVsRun2_Comparison1.AddEntry(graStatRaaVsPt_2040, "#sqrt{#it{s}} = 5.02 TeV, Pb-Pb, 20#minus40%", "FP")
    legendRaaVsPtVsRun2_Comparison1.AddEntry(graStatRaaVsPt_4060, "#sqrt{#it{s}} = 5.02 TeV, Pb-Pb, 40#minus90%", "FP")
    legendRaaVsPtVsRun2_Comparison1.Draw("SAME")

    latexTitle.DrawLatex(0.17, 0.90, "ALICE Work In Progress, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")
    canvasRaaVsPt_Comparison1.Update()
    canvasRaaVsPt_Comparison1.SaveAs("plots/Jpsi_RAA_vs_Pt_Comparison_PbPb_3_centralities.pdf")
    canvasRaaVsPt_Comparison1.SaveAs("plots/Jpsi_RAA_vs_Pt_Comparison_PbPb_3_centralities.png")


    # --- Centrality 0–90% ---
    filePathJpsi_090 = "HEP_Data/jpsi_raa_PbPb_5TeV_0_90_centr.yaml"
    ptCentersRun2_090, ptWidthsRun2_090, jpsiRaaVsPtRun2_090, jpsiStatRaaVsPtRun2_090, jpsiSystRaaVsPtRun2_090 = ExtractFromYaml(filePathJpsi_090)
    jpsiPtSystWidthsRun2_090 = np.repeat(0.2, len(ptCentersRun2_090))

    graStatRaaVsPtRun2_090 = ROOT.TGraphErrors(len(ptCentersRun2_090), np.array(ptCentersRun2_090), np.array(jpsiRaaVsPtRun2_090), np.array(ptWidthsRun2_090), np.array(jpsiStatRaaVsPtRun2_090))
    graSystRaaVsPtRun2_090 = ROOT.TGraphErrors(len(ptCentersRun2_090), np.array(ptCentersRun2_090), np.array(jpsiRaaVsPtRun2_090), jpsiPtSystWidthsRun2_090, np.array(jpsiSystRaaVsPtRun2_090))

    SetGraStat(graStatRaaVsPtRun2_090, 20, ROOT.kGray+1)
    SetGraSyst(graSystRaaVsPtRun2_090, 20, ROOT.kGray+1)

    lineUnityVsPt = ROOT.TLine(0., 1., 12., 1.)
    lineUnityVsPt.SetLineColor(ROOT.kGray+1)
    lineUnityVsPt.SetLineWidth(2)
    lineUnityVsPt.SetLineStyle(ROOT.kDashed)

    canvasRaaVsPtVsRun2_Comparison2 = ROOT.TCanvas("canvasRaaVsPtVsRun2_Comparison2", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    histGridRaaVsPt_Comparison2 = ROOT.TH2D("histGridRaaVsPt_Comaprison2", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}", 100, 0, 12, 100, 0, 1.8)
    histGridRaaVsPt_Comparison2.Draw()
    lineUnityVsPt.Draw("SAME")
    graStatRaaVsPtRun2_090.Draw("EP SAME")
    graSystRaaVsPtRun2_090.Draw("E2P SAME")
    graStatRaaVsPt.Draw("EP SAME")
    graSystRaaVsPt.Draw("E2P SAME")

    legendRaaVsPtVsRun2_Comparison2 = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsPtVsRun2_Comparison2)
    legendRaaVsPtVsRun2_Comparison2.SetTextSize(0.045)
    legendRaaVsPtVsRun2_Comparison2.AddEntry(graSystRaaVsPt, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus90%", "FP")
    legendRaaVsPtVsRun2_Comparison2.AddEntry(graStatRaaVsPtRun2_090, "#sqrt{#it{s}} = 5.02 TeV, Pb-Pb, 0#minus90%", "FP")
    legendRaaVsPtVsRun2_Comparison2.Draw("SAME")

    latexTitle.DrawLatex(0.17, 0.87, "ALICE Work In Progress, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")

    canvasRaaVsPtVsRun2_Comparison2.Update()
    canvasRaaVsPtVsRun2_Comparison2.SaveAs("plots/Jpsi_RAA_vs_Pt_Comparison_PbPb_0_90_centr.pdf")
    canvasRaaVsPtVsRun2_Comparison2.SaveAs("plots/Jpsi_RAA_vs_Pt_Comparison_PbPb_0_90_centr.png")

    # ***** Comparison with pO *****
    x_vect0 = np.array([0.5, 1.5, 2.5, 3.5, 5.0, 8.0], dtype='float64')
    y_vect1 = np.array([0.7010977121877087, 0.8096814426116999, 0.821249497674029, 0.9733794606858844, 0.8910243018805722, 1.250866084626591], dtype='float64')
    ex_vect2 = np.array([0, 0, 0, 0, 0, 0], dtype='float64')
    ey_vect3 = np.array([0.05321832374435311, 0.05126399014605561, 0.05428298032696315, 0.0715452447917095, 0.06328270207447632, 0.1384129436083637], dtype='float64')

    graStatRaaVsPtComparison_pO = ROOT.TGraphErrors(len(x_vect0), x_vect0, y_vect1, ex_vect2, ey_vect3)

    SetGraStat(graStatRaaVsPtComparison_pO, 20, ROOT.kAzure+4)

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
    legendRaaVsPtComparison_pO.AddEntry(graSystRaaVsPt, "#sqrt{#it{s}} = 5.36 TeV, OO,  2.5 < #it{y} < 4", "FP")
    legendRaaVsPtComparison_pO.AddEntry(graStatRaaVsPtComparison_pO, "#sqrt{#it{s}} = 9.62 TeV, pO,  2.15 < #it{y} < 3.65", "FP")
    legendRaaVsPtComparison_pO.Draw("SAME")

    latexTitle.DrawLatex(0.22, 0.87, "ALICE Work In Progress, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}")

    canvasRaaVsPtComparison_pO.Update()
    canvasRaaVsPtComparison_pO.SaveAs("plots/Jpsi_RAA_vs_Pt_Comparison_with_pO.pdf")
    canvasRaaVsPtComparison_pO.SaveAs("plots/Jpsi_RAA_vs_Pt_Comparison_with_pO.png")

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
    histGridRaaVsPt_pPb = ROOT.TH2D("histGridRaaVsPt_pPb", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}", 100, 0, 12, 100, 0, 1.8)
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

    latexTitle.DrawLatex(0.22, 0.87, "ALICE Work In Progress, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}")

    canvasRaaVsPtVsRun2_pPb.Update()
    canvasRaaVsPtVsRun2_pPb.SaveAs("plots/Jpsi_Comparison_pO_pPbRun2.pdf")
    canvasRaaVsPtVsRun2_pPb.SaveAs("plots/Jpsi_Comparison_pO_pPbRun2.png")

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