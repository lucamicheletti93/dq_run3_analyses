import os
import sys
import argparse
import yaml
from pathlib import Path
import math
import numpy as np
import ROOT

def load_config(cfgFileName):
    with open(cfgFileName, 'r') as yml_cfg:
        return yaml.load(yml_cfg, yaml.FullLoader)

# ===================== #
#    Pt distribution    #
# ===================== #

def PtJPsiPbPb5TeV_Func():
    """J/psi pT in Pb-Pb"""
    def func_formula(x, p):
        return p[0] * x[0] / (1. + (x[0] / p[1])**p[2])**p[3]
    return ROOT.TF1("PtPsiPbPb5TeV", func_formula, 0, 20, 4)

# ===================== #
# Rapidity distribution #
# ===================== #

def RapJPsiPbPb5TeV_Func():
    """TF1 for rapidity"""
    def func_formula(x, p):
        return p[0] * math.exp(-0.5 * ((x[0] - p[1]) / p[2])**2)
    return ROOT.TF1("RapPsiPbPb5TeV", func_formula, 2.5, 4, 3)

def RapPsiPbPb5TeV_Original():
    """TF1 for J/psi rapidity in Pb-Pb (original)"""
    def func_formula(x, p):
        return p[0] * math.exp(-0.5 * ((x[0] - p[1]) / p[2])**2)
    return ROOT.TF1("RapPsiPbPb5TeV", func_formula, 2.5, 4, 3)

# ========= #
#   Style   #
# ========= #

def SetHistogram(hist, color):
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetMarkerStyle(20)

def SetLegend(legend):
    legend.SetBorderSize(0)
    legend.SetFillColor(10)
    legend.SetFillStyle(1)
    legend.SetLineStyle(0)
    legend.SetLineColor(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.045)

def load_hist_from_txt(filename, hist_name, hist_title, color=ROOT.kBlack):
    """creating an histo from systematic txt file"""
    data = np.loadtxt(filename, skiprows=1)
    xmins, xmaxs, vals, stat, syst = data.T

    bin_edges = np.concatenate((xmins, [xmaxs[-1]]))
    nbins = len(xmins)

    print(f"--- Contents of {filename} ---")
    for i in range(len(vals)):
        print(f"Bin {i}: x=[{xmins[i]}, {xmaxs[i]}], val={vals[i]}, stat={stat[i]}, syst={syst[i]}")

    hist = ROOT.TH1D(hist_name, hist_title, nbins, bin_edges)

    for i in range(nbins):
        content = vals[i]
        error = np.sqrt(stat[i]**2 + syst[i]**2)
        hist.SetBinContent(i+1, content)
        hist.SetBinError(i+1, error)

    SetHistogram(hist, color)

    # Normalised at 1/area
    area = hist.Integral("width")
    if area > 0:
        hist.Scale(1.0 / area, "width")

    return hist

def do_inputShapes(inputCfg):
    print("Start do_inputShapes")
    
    # Load and compile C++ functions
    print("Loading C++ functions...")
    ROOT.gROOT.ProcessLine(".L treeLoop.C+")

    print("C++ functions compiled successfully!")

    pathFileDataPt = inputCfg["inputs"]["pathFileDataPt"]
    fileNameDataPt = inputCfg["inputs"]["fileNameDataPt"]
    pathFileDataRap = inputCfg["inputs"]["pathFileDataRap"]
    fileNameDataRap = inputCfg["inputs"]["fileNameDataRap"]
    pathFileDataCentr = inputCfg["inputs"]["pathFileDataCentr"]
    fileNameDataCentr = inputCfg["inputs"]["fileNameDataCentr"]
    pathToFileAO2D = inputCfg["inputs"]["pathAO2D"]
    beamType = inputCfg["inputs"]["beamType"]
    ptCut = inputCfg["inputs"]["ptCut"]
    ptMin = inputCfg["inputs"]["ptMin"]
    ptMax = inputCfg["inputs"]["ptMax"]
    centrMin = inputCfg["inputs"]["centrMin"]
    centrMax = inputCfg["inputs"]["centrMax"]
    rapMin = inputCfg["inputs"]["rapMin"]
    rapMax = inputCfg["inputs"]["rapMax"]

    print("Searched parameters")
    print(f"ptMin: {ptMin}")
    print(f"ptMax: {ptMax}")
    print(f"centrMin: {centrMin}")
    print(f"centrMax: {centrMax}")
    print(f"rapMin: {rapMin}")
    print(f"rapMax: {rapMax}")

    isPbPb = (beamType == "PbPb")
    database = ROOT.TDatabasePDG.Instance()
    muPdgCode = 13
    jpsiPdgCode = 443
    massMu = database.GetParticle(muPdgCode).Mass()
    massJpsi = database.GetParticle(jpsiPdgCode).Mass()

    ptBins = np.array(ptMin + [ptMax[-1]], dtype='float64')
    rapBins = np.array(rapMin + [rapMax[-1]], dtype='float64')
    centrBins = np.array(centrMin + [centrMax[-1]], dtype='float64')
    nBinsPt = len(ptBins) - 1
    nBinsRap = len(rapBins) - 1
    nBinsCentr = len(centrBins) - 1

    iterColors = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2]

    histPtJpsiData = load_hist_from_txt(f"{pathFileDataPt}/{fileNameDataPt}", "histPtJpsiData", ";p_{T} (GeV/c);Normalized yield", ROOT.kBlack)
    histRapJpsiData = load_hist_from_txt(f"{pathFileDataRap}/{fileNameDataRap}", "histRapJpsiData", ";y;Normalized yield", ROOT.kBlack)
    if isPbPb:
        histCentrJpsiData = load_hist_from_txt(f"{pathFileDataCentr}/{fileNameDataCentr}", "histCentrJpsiData", ";Centrality (%);Normalized yield", ROOT.kBlack)    
    
    nIterations = 2

    histPtJpsiDataCorr = [None] * (nIterations + 1)
    histRapJpsiDataCorr = [None] * (nIterations + 1)

    histPtJpsiGen = [None] * (nIterations + 1)
    histRapJpsiGen = [None] * (nIterations + 1)
    histCentrJpsiGen = [None] * (nIterations + 1)

    histPtJpsiRec = [None] * (nIterations + 1)
    histRapJpsiRec = [None] * (nIterations + 1)
    histCentrJpsiRec = [None] * (nIterations + 1)

    histPtJpsiAxe = [None] * (nIterations + 1)
    histRapJpsiAxe = [None] * (nIterations + 1)
    histCentrJpsiAxe = [None] * (nIterations + 1)

    fitFunctionPtOriginal = PtJPsiPbPb5TeV_Func()
    fitFunctionPtOriginal.SetParameter(0, 2)
    fitFunctionPtOriginal.SetParameter(1, 3.50274)
    fitFunctionPtOriginal.SetParameter(2, 1.93403)
    fitFunctionPtOriginal.SetParameter(3, 3.96363)
    fitFunctionPtOriginal.SetLineColor(ROOT.kGray+2)

    histFromFuncPtOriginal = fitFunctionPtOriginal.GetHistogram()
    SetHistogram(histFromFuncPtOriginal, ROOT.kGray+2)

    fitFunctionPt = [None] * (nIterations + 1)
    histFromFuncPt = [None] * (nIterations + 1)
    histFromFuncPtRatio = [None] * (nIterations + 1)

    fitFunctionRapOriginal = RapPsiPbPb5TeV_Original()
    fitFunctionRapOriginal.SetParameter(0, 8)
    fitFunctionRapOriginal.SetParameter(1, 0)
    fitFunctionRapOriginal.SetParameter(2, 2.12568)
    fitFunctionRapOriginal.SetLineColor(ROOT.kGray+2)

    histFromFuncRapOriginal = fitFunctionRapOriginal.GetHistogram()
    SetHistogram(histFromFuncRapOriginal, ROOT.kGray+2)

    area = histFromFuncRapOriginal.Integral("width")
    histFromFuncRapOriginal.Scale(1.0 / area)

    fitFunctionRap = [None] * (nIterations + 1)
    histFromFuncRap = [None] * (nIterations + 1)
    histFromFuncRapRatio = [None] * (nIterations + 1)

    fIn = ROOT.TFile(f"{pathToFileAO2D}", "READ")

    for iter in range(nIterations):
        print(f"************* Iteration {iter} *************")

        # Initialize the corrected data distribution to be fitted
        histPtJpsiDataCorr[iter] = histPtJpsiData.Clone(f"histPtJpsiDataCorr_iter_{iter}")
        SetHistogram(histPtJpsiDataCorr[iter], iterColors[iter])
        histRapJpsiDataCorr[iter] = histRapJpsiData.Clone(f"histRapJpsiDataCorr_iter_{iter}")
        SetHistogram(histRapJpsiDataCorr[iter], iterColors[iter])

        # Generated histograms
        histPtJpsiGen[iter] = ROOT.TH1D(f"histPtJpsiGen_iter_{iter}", " ; #it{{p}}_{{T}} (GeV/c)", nBinsPt, ptBins)
        histRapJpsiGen[iter] = ROOT.TH1D(f"histRapJpsiGen_iter_{iter}", " ; #it{{y}}", nBinsRap, rapBins)
        histCentrJpsiGen[iter] = ROOT.TH1D(f"histCentrJpsiGen_iter_{iter}", " ; Centr (%)", nBinsCentr, centrBins)
        SetHistogram(histPtJpsiGen[iter], ROOT.kRed + 1)
        SetHistogram(histRapJpsiGen[iter], ROOT.kRed + 1)
        SetHistogram(histCentrJpsiGen[iter], ROOT.kRed + 1)

        # Reconstructed histograms
        histPtJpsiRec[iter] = ROOT.TH1D(f"histPtJpsiRec_iter_{iter}", " ; #it{{p}}_{{T}} (GeV/c)", nBinsPt, ptBins)
        histRapJpsiRec[iter] = ROOT.TH1D(f"histRapJpsiRec_iter_{iter}", " ; #it{{y}}", nBinsRap, rapBins)
        histCentrJpsiRec[iter] = ROOT.TH1D(f"histCentrJpsiRec_iter_{iter}", " ; Centr (%)", nBinsCentr, centrBins)
        SetHistogram(histPtJpsiRec[iter], ROOT.kBlue)
        SetHistogram(histRapJpsiRec[iter], ROOT.kBlue)
        SetHistogram(histCentrJpsiRec[iter], ROOT.kBlue)

        # Axe histograms
        histPtJpsiAxe[iter] = ROOT.TH1D(f"histPtJpsiAxe_iter_{iter}", " ; #it{{p}}_{{T}} (GeV/c)", nBinsPt, ptBins)
        histRapJpsiAxe[iter] = ROOT.TH1D(f"histRapJpsiAxe_iter_{iter}", " ; #it{{y}}", nBinsRap, rapBins)
        histCentrJpsiAxe[iter] = ROOT.TH1D(f"histCentrJpsiAxe_iter_{iter}", " ; Centr (%)", nBinsCentr, centrBins)
        SetHistogram(histPtJpsiAxe[iter], iterColors[iter])
        SetHistogram(histRapJpsiAxe[iter], iterColors[iter])
        SetHistogram(histCentrJpsiAxe[iter], iterColors[iter])

        # Loop over file trees
        for key in fIn.GetListOfKeys():
            dirName = key.GetName()
            if "DF_" not in dirName:
                continue

            print(f"Processing directory: {dirName}")

            # ================ #
            #  Generated Tree  #
            # ================ #
            treeGen = fIn.Get(f"{dirName}/O2rtdilmtreegen")
            
            if iter == 0:
                ROOT.ProcessGeneratedTree(
                    treeGen,
                    histPtJpsiGen[iter],
                    histRapJpsiGen[iter],
                    histCentrJpsiGen[iter],
                    ROOT.nullptr,
                    ROOT.nullptr,
                    iter,
                    isPbPb,
                    massJpsi
                )
            else:
                ROOT.ProcessGeneratedTree(
                    treeGen,
                    histPtJpsiGen[iter],
                    histRapJpsiGen[iter],
                    histCentrJpsiGen[iter],
                    histFromFuncPtRatio[iter-1],
                    histFromFuncRapRatio[iter-1],
                    iter,
                    isPbPb,
                    massJpsi
                )

            # ==================== #
            #  Reconstructed Tree  #
            # ==================== #
            treeRec = fIn.Get(f"{dirName}/O2rtdilmtreerec")

            if iter == 0:
                ROOT.ProcessReconstructedTree(
                    treeRec,
                    histPtJpsiRec[iter],
                    histRapJpsiRec[iter],
                    histCentrJpsiRec[iter],
                    ROOT.nullptr, 
                    ROOT.nullptr,
                    iter,
                    isPbPb,
                    massMu,
                    ptCut
                )
            else:
                ROOT.ProcessReconstructedTree(
                    treeRec,
                    histPtJpsiRec[iter],
                    histRapJpsiRec[iter],
                    histCentrJpsiRec[iter],
                    histFromFuncPtRatio[iter-1],
                    histFromFuncRapRatio[iter-1],
                    iter,
                    isPbPb,
                    massMu,
                    ptCut
                )

        # JPsi Axe histograms
        print(f"Filling histograms - Iteration {iter}")
        print(f"Gen pT entries: {histPtJpsiGen[iter].GetEntries()}")
        print(f"Rec pT entries: {histPtJpsiRec[iter].GetEntries()}")
        print(f"Gen rap entries: {histRapJpsiGen[iter].GetEntries()}")
        print(f"Rec rap entries: {histRapJpsiRec[iter].GetEntries()}")

        histPtJpsiAxe[iter].Divide(histPtJpsiRec[iter], histPtJpsiGen[iter], 1, 1, "B")
        histRapJpsiAxe[iter].Divide(histRapJpsiRec[iter], histRapJpsiGen[iter], 1, 1, "B")
        
        if isPbPb:
            histCentrJpsiAxe[iter].Divide(histCentrJpsiRec[iter], histCentrJpsiGen[iter], 1, 1, "B")

        # JPsi Axe ratio and normalization
        histPtJpsiDataCorr[iter].Divide(histPtJpsiAxe[iter])
        
        # CORREZIONE: Normalizza dopo la divisione
        area_pt = histPtJpsiDataCorr[iter].Integral("width")
        if area_pt > 0:
            histPtJpsiDataCorr[iter].Scale(1.0 / area_pt, "width")

        fitFunctionPt[iter] = PtJPsiPbPb5TeV_Func()
        fitFunctionPt[iter].SetParameters(1, 2.83941, 2.6687, 2.37032)
        fitFunctionPt[iter].SetLineColor(iterColors[iter])
        histPtJpsiDataCorr[iter].Fit(fitFunctionPt[iter], "R0") 

        histFromFuncPt[iter] = fitFunctionPt[iter].GetHistogram()
        histFromFuncPt[iter].SetName(f"histFromFuncPt_iter_{iter}")

        histFromFuncPtRatio[iter] = histFromFuncPt[iter].Clone(f"histFromFuncPtRatio_iter_{iter}")
        if iter == 0:
            histFromFuncPtRatio[iter].Divide(histFromFuncPtOriginal)
        else:
            histFromFuncPtRatio[iter].Divide(histFromFuncPt[iter-1])
    
        histFromFuncPtRatio[iter].SetLineColor(iterColors[iter])

        histRapJpsiDataCorr[iter].Divide(histRapJpsiAxe[iter])
        
        # CORREZIONE: Normalizza dopo la divisione
        area_rap = histRapJpsiDataCorr[iter].Integral("width")
        if area_rap > 0:
            histRapJpsiDataCorr[iter].Scale(1.0 / area_rap, "width")

        fitFunctionRap[iter] = RapJPsiPbPb5TeV_Func()
        fitFunctionRap[iter].SetParameter(0, 8)       
        fitFunctionRap[iter].SetParameter(1, 0)       
        fitFunctionRap[iter].SetParameter(2, 2.12568)
        fitFunctionRap[iter].SetLineColor(iterColors[iter])
        histRapJpsiDataCorr[iter].Fit(fitFunctionRap[iter], "R0") 

        histFromFuncRap[iter] = fitFunctionRap[iter].GetHistogram()
        histFromFuncRap[iter].SetName(f"histFromFuncRap_iter_{iter}")

        histFromFuncRapRatio[iter] = histFromFuncRap[iter].Clone(f"histFromFuncRapRatio_iter_{iter}")
        if (iter == 0):
            histFromFuncRapRatio[iter].Divide(histFromFuncRapOriginal)
        else:
            histFromFuncRapRatio[iter].Divide(histFromFuncRap[iter-1])

        histFromFuncRapRatio[iter].SetLineColor(iterColors[iter])

    lineUnityPt = ROOT.TLine(0, 1, 20, 1)
    lineUnityRap = ROOT.TLine(2.5, 1, 4, 1)

    print(f"Integral data (raw): {histRapJpsiData.Integral()}")
    print(f"Integral data (after Scale): {histRapJpsiData.Integral("width")}")
    funcIntegral = fitFunctionRapOriginal.Integral(2.5, 4); #TF1 numeric area
    print(f"TF1 integral on [2.5,4]: {funcIntegral}")
    print(f"Histogram from TF1 integral: {histFromFuncRapOriginal.Integral() }")

    #-------------------------Summary of the iterative procedure---------------------------#
    canvasSummaryIterativeTuning = ROOT.TCanvas("canvasSummaryIterativeTuning", "", 1800, 1200)
    canvasSummaryIterativeTuning.Divide(3, 2)

    canvasSummaryIterativeTuning.cd(1)
    ROOT.gPad.SetLogy(True)
    histPtJpsiData.SetStats(False)
    histPtJpsiData.SetTitle("Data")
    histPtJpsiData.Draw("EP")

    canvasSummaryIterativeTuning.cd(2)
    ROOT.gPad.SetLogy(True)
    histPtJpsiDataCorr[0].SetStats(False)
    histPtJpsiDataCorr[0].SetTitle("Data / Ax#epsilon")
    #histPtJpsiDataCorr[0].GetYaxis().SetRangeUser(1e-5, 2)

    legendPtIterativeTuning = ROOT.TLegend(0.70, 0.75, 0.90, 0.90)
    SetLegend(legendPtIterativeTuning)
    legendPtIterativeTuning.AddEntry(histFromFuncPtOriginal, "Original", "L")
    for iter in range(nIterations) :
        histPtJpsiDataCorr[iter].Draw("EP SAME")
        histFromFuncPt[iter].Draw("HIST SAME")
        legendPtIterativeTuning.AddEntry(histFromFuncPt[iter], f"Iteration {iter}")
        
    histFromFuncPtOriginal.Draw("HIST SAME")
    legendPtIterativeTuning.Draw("SAME")

    canvasSummaryIterativeTuning.cd(3)
    histFromFuncPtRatio[0].SetStats(False)
    histFromFuncPtRatio[0].SetTitle("Weights")
    for iter in range(nIterations) :
        #histFromFuncPtRatio[iter].GetYaxis().SetRangeUser(-0.2, 2)
        histFromFuncPtRatio[iter].Draw("HIST SAME")

    lineUnityPt.Draw()

    canvasSummaryIterativeTuning.cd(4)
    histRapJpsiData.SetStats(False)
    histRapJpsiData.SetTitle("Data")
    histRapJpsiData.Draw("EP")

    canvasSummaryIterativeTuning.cd(5)
    histRapJpsiDataCorr[0].SetStats(False)
    histRapJpsiDataCorr[0].SetTitle("Data / Ax#epsilon")
    #histRapJpsiDataCorr[0].GetYaxis().SetRangeUser(0.1, 1.2)

    legendRapIterativeTuning = ROOT.TLegend(0.70, 0.75, 0.90, 0.90)
    SetLegend(legendRapIterativeTuning)
    legendRapIterativeTuning.AddEntry(histFromFuncRapOriginal, "Original", "L")
    for iter in range(nIterations):
        histRapJpsiDataCorr[iter].Draw("EP SAME")
        histFromFuncRap[iter].Draw("HIST SAME")
        legendRapIterativeTuning.AddEntry(histFromFuncRap[iter], f"Iteration {iter}")
    
    histFromFuncRapOriginal.Draw("HIST SAME")
    legendRapIterativeTuning.Draw("SAME")

    canvasSummaryIterativeTuning.cd(6)
    histFromFuncRapRatio[0].SetStats(False)
    histFromFuncRapRatio[0].SetTitle("Weights")
    for iter in range(nIterations):
        #histFromFuncRapRatio[iter].GetYaxis().SetRangeUser(-0.2, 2)
        histFromFuncRapRatio[iter].Draw("HIST SAME")

    lineUnityRap.Draw()

    #-------------------------Summary of the checks---------------------------#
    canvasCheck = ROOT.TCanvas("canvasCheck", "", 1200, 600)
    histFromFuncRapOriginal.Draw("HIST")
    histRapJpsiRec[0].Draw("HIST")

    canvasSummaryAxe = ROOT.TCanvas("canvasSummaryAxe", "", 1200, 600)
    canvasSummaryAxe.Divide(2, 1)

    canvasSummaryAxe.cd(1)
    legendPtAxe = ROOT.TLegend(0.60, 0.75, 0.80, 0.90)
    SetLegend(legendPtAxe)
    histPtJpsiAxe[0].SetStats(False)
    for iter in range(nIterations):
        #histPtJpsiAxe[iter].GetYaxis().SetRangeUser(0, 3.5)
        histPtJpsiAxe[iter].Draw("EP SAME")
        legendPtAxe.AddEntry(histPtJpsiAxe[iter], f"Iteration {iter}")

    legendPtAxe.Draw("SAME")

    canvasSummaryAxe.cd(2)
    legendRapAxe = ROOT.TLegend(0.60, 0.75, 0.80, 0.90)
    SetLegend(legendRapAxe)
    histRapJpsiAxe[0].SetStats(False)
    for iter in range(nIterations):
        #histRapJpsiAxe[iter].GetYaxis().SetRangeUser(0, 3.5)
        histRapJpsiAxe[iter].Draw("EP SAME")
        legendRapAxe.AddEntry(histRapJpsiAxe[iter], f"Iteration {iter}")

    legendRapAxe.Draw("SAME")


    canvasSummaryIterativeTuning.SaveAs("SummaryIterativeTuning.pdf")
    canvasSummaryAxe.SaveAs("SummaryAxe.pdf")

    rootFile = ROOT.TFile("all_canvases.root", "RECREATE")
    canvasSummaryIterativeTuning.Write()
    canvasSummaryAxe.Write()
    c.Write()
    rootFile.Close()



def main():
    parser = argparse.ArgumentParser(description='Input shapes code')
    parser.add_argument('cfgFileName', help='YAML configuration file')
    parser.add_argument('--run', action='store_true', help='Run input shapes code')
    args = parser.parse_args()

    config = load_config(args.cfgFileName)
    if args.run:
        do_inputShapes(config)

if __name__ == "__main__":
    main()