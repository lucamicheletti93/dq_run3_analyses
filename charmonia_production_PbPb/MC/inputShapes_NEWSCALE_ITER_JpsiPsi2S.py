import os
import sys
import argparse
import yaml
from pathlib import Path
import math
import statistics
import itertools
import functools
import operator
import shutil
import logging
from collections import defaultdict, deque, namedtuple, OrderedDict
from typing import List, Tuple, Dict, Any, Optional
import glob
import json
import csv
import numpy as np
from array import array
import matplotlib.pyplot as plt
import ROOT


def load_config(cfgFileName):
    with open(cfgFileName, 'r') as yml_cfg:
        return yaml.load(yml_cfg, yaml.FullLoader)

def PtJPsiPbPb5TeV_Func():
    """
    J/psi pT in Pb—Pb
    """
    return ROOT.TF1(
        "PtPsiPbPb5TeV",
        lambda x, p: p[0] * x[0] / (1. + (x[0] / p[1])**p[2])**p[3], 0, 20, 4
    )

def RapJPsiPbPb5TeV_Func():
    """
    TF1 for rap
    """
    return ROOT.TF1(
        "RapPsiPbPb5TeV",
        lambda x, p: p[0] * math.exp(-0.5 * ((x[0] - p[1]) / p[2])**2), 2.5, 4, 3
    )

# ========= #
#   Style   #
# ========= #

def SetHistogram(hist, color):
    hist.SetLineColor(color)
    hist.SetMarkerColor(color)
    hist.SetMarkerStyle(20)
    hist.SetLineWidth(1)

def SetLegend(legend):
    legend.SetBorderSize(0)
    legend.SetFillColor(10)
    legend.SetFillStyle(1)
    legend.SetLineStyle(0)
    legend.SetLineColor(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.045)

def DrawIntegralText(x, y, value_0, value_last, color=1):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.03)
    latex.SetTextColor(color)
    latex.DrawLatex(x, y, f"Integral(0) = {value_0:.2f}, Integral(last) = {value_last:.2f}")

def load_hist_from_txt(filename, hist_name, hist_title, color=ROOT.kBlack):
    """
    creating an histo from systematic txt file
    """
    data = np.loadtxt(filename, skiprows=1) #skip first row
    xmins, xmaxs, vals, stat, syst = data.T

    bin_edges = np.concatenate((xmins, [xmaxs[-1]]))
    nbins = len(xmins)

    hist = ROOT.TH1D(hist_name, hist_title, nbins, bin_edges)

    for i in range(nbins):
        content = vals[i]
        error = np.sqrt(stat[i]**2 + syst[i]**2)
        hist.SetBinContent(i+1, content)
        hist.SetBinError(i+1, error)

    SetHistogram(hist, color)

    return hist

def do_inputShapes(inputCfg, str_particle):
    print("Start do_inputShapes")

    # Load and compile C++ functions
    print("Loading C++ functions...")
    ROOT.gROOT.ProcessLine(".L treeLoop.C+")
    print("C++ functions compiled successfully!")

    pathFileDataPt = inputCfg["inputs"]["pathFileDataPt"]
    fileNameDataPt_Jpsi = inputCfg["inputs"]["fileNameDataPt_Jpsi"]
    fileNameDataPt_Psi2S = inputCfg["inputs"]["fileNameDataPt_Psi2S"]
    pathFileDataRap = inputCfg["inputs"]["pathFileDataRap"]
    fileNameDataRap_Jpsi = inputCfg["inputs"]["fileNameDataRap_Jpsi"]
    fileNameDataRap_Psi2S = inputCfg["inputs"]["fileNameDataRap_Psi2S"]
    pathFileDataCentr = inputCfg["inputs"]["pathFileDataCentr"]
    fileNameDataCentr = inputCfg["inputs"]["fileNameDataCentr"]
    pathToFileAO2D_GEN = inputCfg["inputs"]["pathAO2D_GEN"]
    pathToFileAO2D_REC = inputCfg["inputs"]["pathAO2D_REC"]
    beamType = inputCfg["inputs"]["beamType"]
    ptCut = inputCfg["inputs"]["ptCut"]
    ptMin = inputCfg["inputs"]["ptMin"]
    ptMax = inputCfg["inputs"]["ptMax"]
    centrMin = inputCfg["inputs"]["centrMin"]
    centrMax = inputCfg["inputs"]["centrMax"]
    rapMin = inputCfg["inputs"]["rapMin"]
    rapMax = inputCfg["inputs"]["rapMax"]

    isPbPb = bool("false")
    database = ROOT.TDatabasePDG.Instance()
    muPdgCode = 13
    jpsiPdgCode = 443
    massMu = database.GetParticle(muPdgCode).Mass()
    massJpsi = database.GetParticle(jpsiPdgCode).Mass()

    fChi2pca   = array('f', [-99999.])
    fSVertex   = array('f', [-99999.])
    fEMC1      = array('f', [-99999.])
    fEMC2      = array('f', [-99999.])
    fPt1       = array('f', [-99999.])
    fPt2       = array('f', [-99999.])
    fPhi1      = array('f', [-99999.])
    fPhi2      = array('f', [-99999.])
    fEta1      = array('f', [-99999.])
    fEta2      = array('f', [-99999.])
    fPtMC1     = array('f', [-99999.])
    fPtMC2     = array('f', [-99999.])
    fPhiMC1    = array('f', [-99999.])
    fPhiMC2    = array('f', [-99999.])
    fEtaMC1    = array('f', [-99999.])
    fEtaMC2    = array('f', [-99999.])
    fPt        = array('f', [-99999.])
    fEta       = array('f', [-99999.])
    fPhi       = array('f', [-99999.])
    fMass      = array('f', [-99999.])
    fSelection = array('i', [-99999])
    fSign      = array('i', [-99999])
    fImpactParameter = array('f', [-99999.])
    fCentFT0C  = array('f', [-99999.])
    fMcDecision= array('I', [0])

    ptBins = np.array(ptMin + [ptMax[-1]], dtype='float64')
    rapBins = np.array(rapMin + [rapMax[-1]], dtype='float64')
    centrBins = np.array(centrMin + [centrMax[-1]], dtype='float64')
    nBinsPt = len(ptBins) - 1
    nBinsRap = len(rapBins) - 1
    nBinsCentr = len(centrBins) - 1

    iterColors = [ROOT.kCyan+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+2, ROOT.kBlack]

    if(str_particle == "Jpsi"):
        histPtJpsiData = load_hist_from_txt(f"{pathFileDataPt}/{fileNameDataPt_Jpsi}", "histPtJpsiData", ";p_{T} (GeV/c);Normalized yield", ROOT.kBlack)
        histRapJpsiData = load_hist_from_txt(f"{pathFileDataRap}/{fileNameDataRap_Jpsi}", "histRapJpsiData", ";y;Normalized yield", ROOT.kBlack)
        if beamType == "PbPb":
            histCentrJpsiData = load_hist_from_txt(f"{pathFileDataCentr}/{fileNameDataCentr}", "histCentrJpsiData", ";Centrality (%);Normalized yield", ROOT.kBlack)    
    elif(str_particle == "Psi2S"):
        histPtJpsiData = load_hist_from_txt(f"{pathFileDataPt}/{fileNameDataPt_Psi2S}", "histPtJpsiData", ";p_{T} (GeV/c);Normalized yield", ROOT.kBlack)
        histRapJpsiData = load_hist_from_txt(f"{pathFileDataRap}/{fileNameDataRap_Psi2S}", "histRapJpsiData", ";y;Normalized yield", ROOT.kBlack)
        if beamType == "PbPb":
            histCentrJpsiData = load_hist_from_txt(f"{pathFileDataCentr}/{fileNameDataCentr}", "histCentrJpsiData", ";Centrality (%);Normalized yield", ROOT.kBlack)  

    """histPtJpsiData.Scale(1.0/histPtJpsiData.Integral(), "width")
    histRapJpsiData.Scale(1.0/histRapJpsiData.Integral(), "width")"""
    histPtJpsiData.Scale(1.0, "width")
    histRapJpsiData.Scale(1.0, "width")
    
    nIterations = 5 ##but for now 1! iteration

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

    histPtJpsiGen_AR = [None] * (nIterations + 1)
    histRapJpsiGen_AR = [None] * (nIterations + 1)
    histCentrJpsiGen_AR = [None] * (nIterations + 1)

    histPtJpsiRec_AR = [None] * (nIterations + 1)
    histRapJpsiRec_AR = [None] * (nIterations + 1)
    histCentrJpsiRec_AR = [None] * (nIterations + 1)

    histPtJpsiAxe_AR = [None] * (nIterations + 1)
    histRapJpsiAxe_AR = [None] * (nIterations + 1)
    histCentrJpsiAxe_AR = [None] * (nIterations + 1)

    histPtJpsiAxe_ratio = [None] * (nIterations + 1)
    histRapJpsiAxe_ratio = [None] * (nIterations + 1)
    histCentrJpsiAxe_ratio = [None] * (nIterations + 1)

    fitFunctionPt = [None] * (nIterations + 1)
    histFromFuncPt = [None] * (nIterations + 1)
    histFromFuncPtRatio = [None] * (nIterations + 1)

    fitFunctionRap = [None] * (nIterations + 1)
    histFromFuncRap = [None] * (nIterations + 1)
    histFromFuncRapRatio = [None] * (nIterations + 1)

    #Correction factor to account for different train efficiencies for generated and reconstructed
    fInGen_AR = ROOT.TFile("/Users/mattei/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/MC/PRODUCTIONS_2024_noVtxCut/AnalysisResults_LHC25c3_518183_GEN.root", "READ")
    fInRec_AR = ROOT.TFile("/Users/mattei/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/MC/PRODUCTIONS_2024_noVtxCut/AnalysisResults_LHC25c3_518184_REC.root", "READ")
    hlistEventsMcGen = fInGen_AR.Get("analysis-event-selection/output")
    histEventsMcGen = (hlistEventsMcGen.FindObject("EventsMC")).FindObject("MCVtxX") 
    hlistEventsMcRec = fInRec_AR.Get("analysis-event-selection/output")
    histEventsMcRec = (hlistEventsMcRec.FindObject("EventsMC")).FindObject("MCVtxX")
    corrFactor = (histEventsMcRec.GetEntries()) / (histEventsMcGen.GetEntries())
    print(f'MC evts. Gen = {histEventsMcGen.GetEntries()} ; MC evts. Rec = {histEventsMcRec.GetEntries()} ; corr. Factor = {corrFactor}')

    fIn_REC = ROOT.TFile(f"{pathToFileAO2D_REC}", "READ")
    fIn_GEN = ROOT.TFile(f"{pathToFileAO2D_GEN}", "READ")

    if(str_particle == "Jpsi"):
        fMcDecison_value = 1
    elif(str_particle == "Psi2S"):
        fMcDecison_value = 2
    else:
        print("error particle!!!")

    for iter in range(nIterations):
        print(f"************* Iteration {iter} *************")

        # ----------- Initialize the corrected data distribution to be fitted ---------- #
        histPtJpsiDataCorr[iter] = histPtJpsiData.Clone(f"histPtJpsiDataCorr_iter_{iter}")
        SetHistogram(histPtJpsiDataCorr[iter], iterColors[iter])
        histRapJpsiDataCorr[iter] = histRapJpsiData.Clone(f"histRapJpsiDataCorr_iter_{iter}")
        SetHistogram(histRapJpsiDataCorr[iter], iterColors[iter])

        # ---------------------------- Generated histograms ---------------------------- #
        histPtJpsiGen[iter] = ROOT.TH1D(f"histPtJpsiGen_iter_{iter}", " ; #it{p}_{T} GeV/c", nBinsPt, ptBins)
        histRapJpsiGen[iter] = ROOT.TH1D(f"histRapJpsiGen_iter_{iter}", " ; #it{y}", nBinsRap, rapBins)
        histCentrJpsiGen[iter] = ROOT.TH1D(f"histCentrJpsiGen_iter_{iter}", " ; Centr %", nBinsCentr, centrBins)
        SetHistogram(histPtJpsiGen[iter], ROOT.kCyan+1)
        SetHistogram(histRapJpsiGen[iter], ROOT.kCyan+1)
        SetHistogram(histCentrJpsiGen[iter], ROOT.kCyan+1)

        if iter == 0:
            # ---------------------------- Generated histograms for "original" reference---------------------------- #
            histPtJpsiGen_allbins = ROOT.TH1D(f"histPtJpsiGen_allbins", " ; #it{p}_{T} GeV/c", 100, 0, 15)
            histRapJpsiGen_allbins = ROOT.TH1D(f"histRapJpsiGen_allbins", " ; #it{y}", nBinsRap, rapBins)
            histCentrJpsiGen_allbins = ROOT.TH1D(f"histCentrJpsiGen_allbins", " ; Centr %", nBinsCentr, centrBins)
            SetHistogram(histPtJpsiGen[iter], ROOT.kCyan+1)
            SetHistogram(histRapJpsiGen[iter], ROOT.kCyan+1)
            SetHistogram(histCentrJpsiGen[iter], ROOT.kCyan+1)

        # -------------------------- Reconstructed histograms -------------------------- #
        histPtJpsiRec[iter] = ROOT.TH1D(f"histPtJpsiRec_iter_{iter}", " ; #it{p}_{T} GeV/c", nBinsPt, ptBins)
        histRapJpsiRec[iter] = ROOT.TH1D(f"histRapJpsiRec_iter_{iter}", " ; #it{y}", nBinsRap, rapBins)
        histCentrJpsiRec[iter] = ROOT.TH1D(f"histCentrJpsiRec_iter_{iter}", " ; Centr %", nBinsCentr, centrBins)
        SetHistogram(histPtJpsiRec[iter], ROOT.kBlue)
        SetHistogram(histRapJpsiRec[iter], ROOT.kBlue)
        SetHistogram(histCentrJpsiRec[iter], ROOT.kBlue)

        # ------------------------------- Axe histograms ------------------------------- #
        histPtJpsiAxe[iter] = ROOT.TH1D(f"histPtJpsiAxe_iter_{iter}", " ; #it{p}_{T} GeV/c", nBinsPt, ptBins)
        histRapJpsiAxe[iter] = ROOT.TH1D(f"histRapJpsiAxe_iter_{iter}", " ; #it{y}", nBinsRap, rapBins)
        histCentrJpsiAxe[iter] = ROOT.TH1D(f"histCentrJpsiAxe_iter_{iter}", " ; Centr %", nBinsCentr, centrBins)
        SetHistogram(histPtJpsiAxe[iter], iterColors[iter])
        SetHistogram(histRapJpsiAxe[iter], iterColors[iter])
        SetHistogram(histCentrJpsiAxe[iter], iterColors[iter])

        #----------------------------Loop over file trees---------------------------#
        index = 0
        for key in fIn_GEN.GetListOfKeys():
            dirName = key.GetName()
            if "DF_" not in dirName:
                continue

            print(f"Processing directory: {dirName}")

            # ================ #
            #  Generated Tree  #
            # ================ #
            treeGen = fIn_GEN.Get(f"{dirName}/O2rtdilmtreegen")
            
            if iter == 0:
                ROOT.ProcessGeneratedTree(
                    fMcDecison_value, #fMcDecision of Jpsi
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
                    fMcDecison_value, #fMcDecision of Jpsi
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

        for key in fIn_REC.GetListOfKeys():
            dirName = key.GetName()
            if "DF_" not in dirName:
                continue

            print(f"Processing directory: {dirName}")
            # ==================== #
            #  Reconstructed Tree  #
            # ==================== #
            treeRec = fIn_REC.Get(f"{dirName}/O2rtdilmtreerec")

            if iter == 0:
                ROOT.ProcessReconstructedTree(
                    fMcDecison_value, #fMcDecision of Jpsi
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
                    fMcDecison_value, #fMcDecision of Jpsi
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
        
        histPtJpsiGen[iter].Scale(corrFactor)
        histRapJpsiGen[iter].Scale(corrFactor)
       

        #-----------------------JPsi Axe histograms------------------------#

        """histPtJpsiAxe[iter].Divide(histPtJpsiRec[iter], histPtJpsiGen[iter], 1, 1, "B")
        histRapJpsiAxe[iter].Divide(histRapJpsiRec[iter], histRapJpsiGen[iter], 1, 1, "B")"""

        nBinsPt = histPtJpsiGen[iter].GetNbinsX()
        nBinsRap = histRapJpsiGen[iter].GetNbinsX()

        # pT efficiency
        for bin in range(1, nBinsPt + 1):
            rec = histPtJpsiRec[iter].GetBinContent(bin)
            gen = histPtJpsiGen[iter].GetBinContent(bin)

            eff = 0.0
            err = 0.0

            if gen != 0:
                eff = rec / gen
                err = math.sqrt(eff * (1.0 - eff) / gen)

            histPtJpsiAxe[iter].SetBinContent(bin, eff)
            histPtJpsiAxe[iter].SetBinError(bin, err)

        # Rapidity efficiency
        for bin in range(1, nBinsRap + 1):
            rec = histRapJpsiRec[iter].GetBinContent(bin)
            gen = histRapJpsiGen[iter].GetBinContent(bin)

            eff = 0.0
            err = 0.0

            if gen != 0:
                eff = rec / gen
                err = math.sqrt(eff * (1.0 - eff) / gen)

            histRapJpsiAxe[iter].SetBinContent(bin, eff)
            histRapJpsiAxe[iter].SetBinError(bin, err)
        

        #-------------------------JPsi Axe ratio---------------------------#
        
        print("*************************************")
        print(f'[{iter}], {histPtJpsiDataCorr[iter].GetBinContent(1)}, {histPtJpsiAxe[iter].GetBinContent(1)}')

        histPtJpsiDataCorr[iter].Divide(histPtJpsiAxe[iter])
        print(f'Corrected yield = {histPtJpsiDataCorr[iter].GetBinContent(1)}')
        print("*************************************")

        fitFunctionPt[iter] = PtJPsiPbPb5TeV_Func()
        fitFunctionPt[iter].SetParameters(histPtJpsiData.Integral("width"), 2.83941, 2.6687, 2.37032)
        fitFunctionPt[iter].SetParLimits(1, 1e-3, 10)   # p1 > 0
        fitFunctionPt[iter].SetParLimits(2, 0, 10)      # p2 >= 0
        fitFunctionPt[iter].SetParLimits(3, 0, 10)      # p3 reasonable
        fitFunctionPt[iter].SetLineColor(iterColors[iter])
        histPtJpsiDataCorr[iter].Fit(fitFunctionPt[iter], "R0") 

        histFromFuncPt[iter] = fitFunctionPt[iter].GetHistogram()
        histFromFuncPt[iter].SetName(f"histFromFuncPt_iter_{iter}")
        SetHistogram(histFromFuncPt[iter], iterColors[iter])

        histFromFuncPtRatio[iter] = histFromFuncPt[iter].Clone(f"histFromFuncPtRatio_iter_{iter}")
        #-----------------------Generated "original" distributions-----------------------
        #The iteration 0 is compared to the original, aka the generated at iteration 0
        #Pt
        histPtGen0original = histPtJpsiGen[0].Clone("histPtGen0original")
        histPtGen0original.Scale(histPtJpsiDataCorr[iter].Integral("width")/histPtGen0original.Integral(), "width")
        print("#####################STO SCALANDO!!!")
        index +=1
        SetHistogram(histPtGen0original, ROOT.kGray+2)
        histPtGen0original.SetLineColor(ROOT.kGray+2)
        fitFunctionPtGen0original = PtJPsiPbPb5TeV_Func()
        fitFunctionPtGen0original.SetParameters(histPtJpsiDataCorr[iter].Integral("width"), 4.52544, 1.78863, 3.99586)
        fitFunctionPtGen0original.SetLineColor(ROOT.kGray+2)
        print("---------fitFunctionPtGen0-------------")
        fitFunctionPtGen0original.SetParLimits(1, 1e-3, 10)  # p[1] positive
        fitFunctionPtGen0original.SetParLimits(2, 0, 10)     # p[2] positive
        fitFunctionPtGen0original.SetParLimits(3, 0, 10)     # p[3] reasonable
        histPtGen0original.Fit(fitFunctionPtGen0original, "R0") 
        histFromFuncPtGen0original = fitFunctionPtGen0original.GetHistogram()
        histFromFuncPtGen0original.SetName(f"histFromFuncPtGen0original")
        histFromFuncPtGen0original.SetLineColor(ROOT.kGray+2)
        SetHistogram(histFromFuncPtGen0original, ROOT.kGray+2)
        if (iter == 0):
            histFromFuncPtRatio[iter].Divide(histFromFuncPtGen0original)
        else:
            #histFromFuncPtRatio[iter].Divide(histFromFuncPt[iter-1])
            histFromFuncPtRatio[iter].Divide(histFromFuncPtGen0original)
    
        histFromFuncPtRatio[iter].SetLineColor(iterColors[iter])

        histRapJpsiDataCorr[iter].Divide(histRapJpsiAxe[iter])

        fitFunctionRap[iter] = RapJPsiPbPb5TeV_Func()
        fitFunctionRap[iter].SetParameters(histRapJpsiData.Integral("width"), 8, 2.12568)
        fitFunctionRap[iter].SetLineColor(iterColors[iter])
        histRapJpsiDataCorr[iter].Fit(fitFunctionRap[iter], "R0") 

        histFromFuncRap[iter] = fitFunctionRap[iter].GetHistogram()
        histFromFuncRap[iter].SetName(f"histFromFuncRap_iter_{iter}")
        SetHistogram(histFromFuncRap[iter], iterColors[iter])

        histFromFuncRapRatio[iter] = histFromFuncRap[iter].Clone(f"histFromFuncRapRatio_iter_{iter}")
        #-----------------------Generated "original" distributions-----------------------
        #The iteration 0 is compared to the original, aka the generated at iteration 0
        #Rap
        histRapGen0original = histRapJpsiGen[0].Clone("histRapGen0original")
        SetHistogram(histRapGen0original, ROOT.kGray+2)
        histRapGen0original.Scale(histRapJpsiDataCorr[iter].Integral("width")/histRapGen0original.Integral(), "width")
        fitFunctionRapGen0original = RapJPsiPbPb5TeV_Func()
        fitFunctionRapGen0original.SetParameters(histRapJpsiDataCorr[iter].Integral("width"), 2.83941, 2.6687, 2.37032)
        fitFunctionRapGen0original.SetLineColor(ROOT.kGray+2)
        print("---------fitFunctionRapGen0-------------")
        histRapGen0original.Fit(fitFunctionRapGen0original, "R0") 
        histFromFuncRapGen0original = fitFunctionRapGen0original.GetHistogram()
        histFromFuncRapGen0original.SetName(f"histFromFuncRapGen0original")
        histFromFuncRapGen0original.SetLineColor(ROOT.kGray+2)
        SetHistogram(histFromFuncRapGen0original, ROOT.kGray+2)
        if (iter == 0):
            histFromFuncRapRatio[iter].Divide(histFromFuncRapGen0original)
        else:
            #histFromFuncRapRatio[iter].Divide(histFromFuncRap[iter-1])
            histFromFuncRapRatio[iter].Divide(histFromFuncRapGen0original)

        histFromFuncRapRatio[iter].SetLineColor(iterColors[iter])

    lineUnityPt = ROOT.TLine(0, 1, 20, 1)
    lineUnityRap = ROOT.TLine(2.5, 1, 4, 1)

    """print(f"Integral data (raw): {histRapJpsiData.Integral()}")
    print(f"Integral data (after Scale): {histRapJpsiData.Integral("width")}")"""

    #-------------------------Summary of the iterative procedure---------------------------#
    canvasSummaryIterativeTuning = ROOT.TCanvas("canvasSummaryIterativeTuning", "", 1800, 1200)
    canvasSummaryIterativeTuning.Divide(4, 2)

    #####
    canvasSummaryIterativeTuning.cd(1)
    ROOT.gPad.SetLogy(True)
    histPtJpsiData.SetStats(False)
    histPtJpsiData.SetTitle("Data")
    #histPtJpsiData.Scale(1.0, "width")
    histPtJpsiData.Draw("EP")

    #####
    canvasSummaryIterativeTuning.cd(2)
    legendPtAxe = ROOT.TLegend(0.60, 0.75, 0.80, 0.90)
    SetLegend(legendPtAxe)
    histPtJpsiAxe[0].SetStats(False)
    for iter in range(nIterations):
        #histPtJpsiAxe[iter].GetYaxis().SetRangeUser(0, 3.5)
        histPtJpsiAxe[iter].Draw("EP SAME")
    legendPtAxe.Draw("SAME")

    #####
    canvasSummaryIterativeTuning.cd(3)
    #ROOT.gPad.SetLogy(True)
    histPtJpsiDataCorr[0].SetStats(False)
    histPtJpsiDataCorr[0].SetTitle("Data / Ax#epsilon")

    legendPtIterativeTuning = ROOT.TLegend(0.45, 0.75, 0.89, 0.89)
    SetLegend(legendPtIterativeTuning)
    legendPtIterativeTuning.AddEntry(histFromFuncPtGen0original, "HistGen, Iteration 0", "L")
    for iter in range(nIterations) :
        histPtJpsiDataCorr[iter].Draw("EP SAME")
        histFromFuncPt[iter].Draw("HIST SAME")
        legendPtIterativeTuning.AddEntry(histFromFuncPt[iter], f"Iteration {iter}")
        fitFunctionPt[iter].Draw("SAME")
        
    histPtGen0original.Draw("EP SAME")
    histFromFuncPtGen0original.Draw("HIST SAME")
    fitFunctionPtGen0original.Draw("SAME")
    legendPtIterativeTuning.Draw("SAME")

    #####
    canvasSummaryIterativeTuning.cd(4)
    histFromFuncPtRatio[0].SetStats(False)
    histFromFuncPtRatio[0].SetTitle("Weights")
    for iter in range(nIterations) :
        histFromFuncPtRatio[iter].GetYaxis().SetRangeUser(0.25, 1.25)
        histFromFuncPtRatio[iter].Draw("HIST SAME")


    lineUnityPt.Draw()

    #####
    canvasSummaryIterativeTuning.cd(5)
    histRapJpsiData.SetStats(False)
    histRapJpsiData.SetTitle("Data")
    histRapJpsiData.Draw("EP")

    #####
    canvasSummaryIterativeTuning.cd(6)
    legendRapAxe = ROOT.TLegend(0.60, 0.75, 0.80, 0.90)
    SetLegend(legendRapAxe)
    histRapJpsiAxe[0].SetStats(False)
    for iter in range(nIterations):
        histRapJpsiAxe[iter].Draw("EP SAME")
    legendRapAxe.Draw("SAME")

    #####
    canvasSummaryIterativeTuning.cd(7)
    histRapJpsiDataCorr[0].SetStats(False)
    histRapJpsiDataCorr[0].SetTitle("Data / Ax#epsilon")

    legendRapIterativeTuning = ROOT.TLegend(0.45, 0.75, 0.89, 0.89)
    SetLegend(legendRapIterativeTuning)
    legendRapIterativeTuning.AddEntry(histFromFuncRapGen0original, "HistGen, Iteration 0", "L")
    for iter in range(nIterations):
        histRapJpsiDataCorr[iter].Draw("EP SAME")
        histFromFuncRap[iter].Draw("HIST SAME")
        legendRapIterativeTuning.AddEntry(histFromFuncRap[iter], f"Iteration {iter}")
        fitFunctionRap[iter].Draw("SAME")

    histRapGen0original.Draw("EP SAME")
    histFromFuncRapGen0original.Draw("HIST SAME")
    fitFunctionRapGen0original.Draw("SAME")
    legendRapIterativeTuning.Draw("SAME")
    
    #####
    canvasSummaryIterativeTuning.cd(8)
    histFromFuncRapRatio[0].SetStats(False)
    histFromFuncRapRatio[0].SetTitle("Weights")
    for iter in range(nIterations):
        histFromFuncRapRatio[iter].GetYaxis().SetRangeUser(0, 2)
        histFromFuncRapRatio[iter].Draw("HIST SAME")

    lineUnityRap.Draw()

    canvasSummaryIterativeTuning.Update()

    canvasSummaryGenRec = ROOT.TCanvas("canvasSummaryGenRec", "", 1200, 600)
    canvasSummaryGenRec.Divide(2, 2)

    canvasSummaryGenRec.cd(1)
    ROOT.gPad.SetLogy(True)
    legendPtAxe = ROOT.TLegend(0.60, 0.75, 0.80, 0.90)
    SetLegend(legendPtAxe)
    histPtJpsiGen[0].SetStats(False)
    histPtJpsiGen[0].SetTitle("Generated")
    histPtJpsiGen[0].GetYaxis().SetRangeUser(1e2,1e7)
    for iter in range(nIterations):
        SetHistogram(histPtJpsiGen[iter], iterColors[iter])
        #histPtJpsiGen[iter].Scale(1.0/histPtJpsiGen[iter].Integral(), "width")
        histPtJpsiGen[iter].Scale(1.0, "width")
        histPtJpsiGen[iter].Draw("EP SAME")
        legendPtAxe.AddEntry(histPtJpsiGen[iter], f"Iteration {iter}")
    last_integral_PtGen = histPtJpsiGen[nIterations-1].Integral("width")
    first_integral_PtGen = histPtJpsiGen[0].Integral("width")
    DrawIntegralText(0.15, 0.85, first_integral_PtGen, last_integral_PtGen)
    legendPtAxe.Draw("SAME")

    canvasSummaryGenRec.cd(2)
    ROOT.gPad.SetLogy(True)
    legendPtAxe = ROOT.TLegend(0.60, 0.75, 0.80, 0.90)
    SetLegend(legendPtAxe)
    histPtJpsiRec[0].SetStats(False)
    histPtJpsiRec[0].SetTitle("Reconstructed")
    histPtJpsiRec[0].GetYaxis().SetRangeUser(1e2,1e7)
    for iter in range(nIterations):
        SetHistogram(histPtJpsiRec[iter], iterColors[iter])
        #histPtJpsiRec[iter].Scale(1.0/histPtJpsiRec[iter].Integral(), "width")
        histPtJpsiRec[iter].Scale(1.0, "width")
        histPtJpsiRec[iter].Draw("EP SAME")
        legendPtAxe.AddEntry(histPtJpsiRec[iter], f"Iteration {iter}")
    last_integral_PtRec = histPtJpsiRec[nIterations-1].Integral("width")
    first_integral_PtRec = histPtJpsiRec[0].Integral("width")
    DrawIntegralText(0.15, 0.85, first_integral_PtRec, last_integral_PtRec)
    legendPtAxe.Draw("SAME")

    canvasSummaryGenRec.cd(3)
    ROOT.gPad.SetLogy(True)
    legendRapAxe = ROOT.TLegend(0.60, 0.75, 0.80, 0.90)
    SetLegend(legendRapAxe)
    histRapJpsiGen[0].SetStats(False)
    histRapJpsiGen[0].SetTitle("Generated")
    for iter in range(nIterations):
        SetHistogram(histRapJpsiGen[iter], iterColors[iter])
        #histRapJpsiGen[iter].Scale(1.0/histRapJpsiGen[iter].Integral(), "width")
        histRapJpsiGen[iter].Scale(1.0, "width")
        histRapJpsiGen[iter].Draw("EP SAME")
        legendRapAxe.AddEntry(histRapJpsiGen[iter], f"Iteration {iter}")
    last_integral_RapGen = histRapJpsiGen[nIterations-1].Integral("width")
    first_integral_RapGen = histRapJpsiGen[0].Integral("width")
    DrawIntegralText(0.15, 0.85, first_integral_RapGen, last_integral_RapGen)
    legendRapAxe.Draw("SAME")

    canvasSummaryGenRec.cd(4)
    ROOT.gPad.SetLogy(True)
    legendRapAxe = ROOT.TLegend(0.60, 0.75, 0.80, 0.90)
    SetLegend(legendRapAxe)
    histRapJpsiRec[0].SetStats(False)
    histRapJpsiRec[0].SetTitle("Reconstructed")
    for iter in range(nIterations):
        SetHistogram(histRapJpsiRec[iter], iterColors[iter])
        #histRapJpsiRec[iter].Scale(1.0/histRapJpsiRec[iter].Integral(), "width")
        histRapJpsiRec[iter].Scale(1.0, "width")
        histRapJpsiRec[iter].Draw("EP SAME")
        legendRapAxe.AddEntry(histRapJpsiRec[iter], f"Iteration {iter}")
    last_integral_RapRec = histRapJpsiRec[nIterations-1].Integral("width")
    first_integral_RapRec = histRapJpsiRec[0].Integral("width")
    DrawIntegralText(0.15, 0.85, first_integral_RapRec, last_integral_RapRec)
    legendRapAxe.Draw("SAME")

    canvasSummaryGenRec.Update()

    canvasSummaryGenRec.SaveAs(f"SummaryGenRec_NEWSCALE_ITER_{str_particle}.pdf")
    canvasSummaryGenRec.SaveAs(f"SummaryGenRec_NEWSCALE_ITER_{str_particle}.root")

    canvasOriginal = ROOT.TCanvas("canvasOriginal", "", 1200, 600)
    canvasOriginal.cd()
    histRapGen0original.Draw("EP")
    canvasOriginal.SaveAs(f"CanvasOriginal_NEWSCALE_ITER_{str_particle}.root")
    """
    print("#########################")
    print(histFromFuncPt[0].GetNbinsX())
    print(histFromFuncRap[0].GetNbinsX())
    print("#########Summary#########")
    for iter in range(nIterations):
        print(iter, histPtJpsiData.GetBinContent(1), histPtJpsiAxe[iter].GetBinContent(1), histPtJpsiDataCorr[iter].GetBinContent(1))
    print("#########################")

    print("#########Summary 2#########")
    for iter in range(1, nIterations):
        print(iter, histPtJpsiAxe[iter].GetBinContent(1)/histPtJpsiAxe[iter-1].GetBinContent(1), histPtJpsiDataCorr[iter-1].GetBinContent(1)/histPtJpsiDataCorr[iter].GetBinContent(1))
    print("#########################")

    print("#########Summary 3#########")
    for iter in range(nIterations):
        print("iteration ", iter, "      Axe, pT, bin1 ", histPtJpsiAxe[iter].GetBinContent(1), "     Axe, rap, bin1 ", histRapJpsiAxe[iter].GetBinContent(1))
    print("#########################")

    print("#########Summary 4#########")
    for iter in range(nIterations):
        print("iteration ", iter, "      (Axe_iter - Axe_0)/Axe_0, pT, bin1 ", ROOT.TMath.Abs(histPtJpsiAxe[iter].GetBinContent(1)-histPtJpsiAxe[0].GetBinContent(1))/histPtJpsiAxe[0].GetBinContent(1), "     (Axe_iter - Axe_0)/Axe_0, rap, bin1 ", ROOT.TMath.Abs(histRapJpsiAxe[iter].GetBinContent(1)-histRapJpsiAxe[0].GetBinContent(1))/histRapJpsiAxe[0].GetBinContent(1))
    print("#########################")
    """
    #input()

    canvasSummaryIterativeTuning.SaveAs(f"SummaryIterativeTuning_NEWSCALE_ITER_{str_particle}.pdf")
    canvasSummaryIterativeTuning.SaveAs(f"SummaryIterativeTuning_NEWSCALE_ITER_{str_particle}.root")

    #After running the procedure, the histogrmas needed are: the iteration 0 and the iteration nIterations-1 for both Pt and Rap
    histAxePt_iter0 = histPtJpsiAxe[0].Clone(f"histAxePt_iter0_{str_particle}")
    histAxePt_iter0.SetDirectory(0)
    histAxePt_iterLast = histPtJpsiAxe[nIterations-1].Clone(f"histAxePt_iterLast_{str_particle}")
    histAxePt_iterLast.SetDirectory(0)
    histAxeRap_iter0 = histRapJpsiAxe[0].Clone(f"histAxeRap_iter0_{str_particle}")
    histAxeRap_iter0.SetDirectory(0)
    histAxeRap_iterLast = histRapJpsiAxe[nIterations-1].Clone(f"histAxeRap_iterLast_{str_particle}")
    histAxeRap_iterLast.SetDirectory(0)

    
    #Integrals
    axePt_int_iter0 = histPtJpsiRec[0].Integral("width")/histPtJpsiGen[0].Integral("width")
    axeRap_int_iter0 = histRapJpsiRec[0].Integral("width")/histRapJpsiGen[0].Integral("width")

    axePt_int_iterLast = histPtJpsiRec[nIterations-1].Integral("width")/histPtJpsiGen[nIterations-1].Integral("width")
    axeRap_int_iterLast = histRapJpsiRec[nIterations-1].Integral("width")/histRapJpsiGen[nIterations-1].Integral("width")

    return histAxePt_iter0, histAxePt_iterLast, histAxeRap_iter0, histAxeRap_iterLast, axePt_int_iter0, axeRap_int_iter0, axePt_int_iterLast, axeRap_int_iterLast



def main():
    parser = argparse.ArgumentParser(description='Input shapes code')
    parser.add_argument('cfgFileName', help='YAML configuration file')
    parser.add_argument('--run', action='store_true', help='Run input shapes code')
    args = parser.parse_args()

    config = load_config(args.cfgFileName)
    if args.run:
        histAxePtJpsi_iter0,  histAxePtJpsi_iterLast,  histAxeRapJpsi_iter0,  histAxeRapJpsi_iterLast,  axePtJpsi_int_iter0,  axeRapJpsi_int_iter0,  axePtJpsi_int_iterLast,  axeRapJpsi_int_iterLast  = do_inputShapes(config, "Jpsi")
        histAxePtPsi2S_iter0, histAxePtPsi2S_iterLast, histAxeRapPsi2S_iter0, histAxeRapPsi2S_iterLast, axePtPsi2S_int_iter0, axeRapPsi2S_int_iter0, axePtPsi2S_int_iterLast, axeRapPsi2S_int_iterLast = do_inputShapes(config, "Psi2S")

    legend = ROOT.TLegend(0.14, 0.76, 0.51, 0.89)
    SetLegend(legend)

    #Summarize Jpsi, Psi2S, ratio
    canvasFinalSummary = ROOT.TCanvas("canvasFinalSummary", "", 2800, 1200)
    canvasFinalSummary.Divide(4, 2)

    canvasFinalSummary.cd(1)
    histAxePtJpsi_iter0.SetTitle("J/#psi, Pt")
    histAxePtJpsi_iter0.GetYaxis().SetTitle("A#times#varepsilon")
    SetHistogram(histAxePtJpsi_iter0, ROOT.kCyan+1)
    histAxePtJpsi_iter0.Draw("ep")
    SetHistogram(histAxePtJpsi_iterLast, ROOT.kBlack)
    histAxePtJpsi_iterLast.Draw("epsame")
    legend.AddEntry(histAxePtJpsi_iter0, "iteration 0")
    legend.AddEntry(histAxePtJpsi_iterLast, "iteration last")
    legend.Draw("same")

    canvasFinalSummary.cd(2)
    histAxePtPsi2S_iter0.SetTitle("#psi(2S), Pt")
    histAxePtPsi2S_iter0.GetYaxis().SetTitle("A#times#varepsilon")
    SetHistogram(histAxePtPsi2S_iter0, ROOT.kCyan+1)
    histAxePtPsi2S_iter0.Draw("ep")
    SetHistogram(histAxePtPsi2S_iterLast, ROOT.kBlack)
    histAxePtPsi2S_iterLast.Draw("epsame")

    canvasFinalSummary.cd(3)
    histAxePtRatio_iter0 = histAxePtPsi2S_iter0.Clone("histAxePtRatio_iter0")
    histAxePtRatio_iter0.Divide(histAxePtJpsi_iter0)
    histAxePtRatio_iterLast = histAxePtPsi2S_iterLast.Clone("histAxePtRatio_iterLast")
    histAxePtRatio_iterLast.Divide(histAxePtJpsi_iterLast)
    histAxePtRatio_iter0.SetTitle("#psi(2S) / J/#psi, Pt")
    histAxePtRatio_iter0.GetYaxis().SetTitle("A#times#varepsilon ratio")
    histAxePtRatio_iter0.GetYaxis().SetRangeUser(0.85, 1.15)
    SetHistogram(histAxePtRatio_iter0, ROOT.kCyan+1)
    histAxePtRatio_iter0.Draw("ep")
    SetHistogram(histAxePtRatio_iterLast, ROOT.kBlack)
    histAxePtRatio_iterLast.GetYaxis().SetRangeUser(0.85, 1.15)
    histAxePtRatio_iterLast.Draw("epsame")

    canvasFinalSummary.cd(4)
    histAxePtRatio_iter0overLast = histAxePtRatio_iter0.Clone("histAxePtRatio_iter0overLast")
    for i in range(1, histAxePtRatio_iter0overLast.GetNbinsX() + 1):
        #histAxePtRatio_iter0overLast.SetBinContent(i, ROOT.TMath.Abs(histAxePtRatio_iter0.GetBinContent(i) - histAxePtRatio_iterLast.GetBinContent(i))/ (histAxePtRatio_iter0.GetBinContent(i) + histAxePtRatio_iterLast.GetBinContent(i)))
        histAxePtRatio_iter0overLast.SetBinContent(i, 0.5*ROOT.TMath.Abs(histAxePtRatio_iter0.GetBinContent(i) - histAxePtRatio_iterLast.GetBinContent(i))/ (histAxePtRatio_iterLast.GetBinContent(i)))
        histAxePtRatio_iter0overLast.SetBinError(i, 0)
    histAxePtRatio_iter0overLast.SetTitle("(iteration0-iterationLast) / 2*iterationLast")
    histAxePtRatio_iter0overLast.SetTitleSize(0.1, "t") 
    histAxePtRatio_iter0overLast.GetYaxis().SetRangeUser(0, 0.015)
    histAxePtRatio_iter0overLast.Draw("hp")

    canvasFinalSummary.cd(5)
    histAxeRapJpsi_iter0.SetTitle("J/#psi, Rap")
    SetHistogram(histAxeRapJpsi_iter0, ROOT.kCyan+1)
    histAxeRapJpsi_iter0.Draw("ep")
    SetHistogram(histAxeRapJpsi_iterLast, ROOT.kBlack)
    histAxeRapJpsi_iterLast.Draw("epsame")

    canvasFinalSummary.cd(6)
    histAxeRapPsi2S_iter0.SetTitle("#psi(2S), Rap")
    SetHistogram(histAxeRapPsi2S_iter0, ROOT.kCyan+1)
    histAxeRapPsi2S_iter0.Draw("ep")
    SetHistogram(histAxeRapPsi2S_iterLast, ROOT.kBlack)
    histAxeRapPsi2S_iterLast.Draw("epsame")

    canvasFinalSummary.cd(7)
    histAxeRapRatio_iter0 = histAxeRapPsi2S_iter0.Clone("histAxeRapRatio_iter0")
    histAxeRapRatio_iter0.Divide(histAxeRapJpsi_iter0)
    histAxeRapRatio_iterLast = histAxeRapPsi2S_iterLast.Clone("histAxeRapRatio_iterLast")
    histAxeRapRatio_iterLast.Divide(histAxeRapJpsi_iterLast)
    histAxeRapRatio_iter0.SetTitle("#psi(2S) / J/#psi, Rap")
    histAxeRapRatio_iter0.GetYaxis().SetRangeUser(1.0, 1.2)
    SetHistogram(histAxeRapRatio_iter0, ROOT.kCyan+1)
    histAxeRapRatio_iter0.Draw("ep")
    histAxeRapRatio_iterLast.GetYaxis().SetRangeUser(1.0, 1.2)
    SetHistogram(histAxeRapRatio_iterLast, ROOT.kBlack)
    histAxeRapRatio_iterLast.Draw("epsame")
    canvasFinalSummary.Update()

    canvasFinalSummary.cd(8)
    histAxeRapRatio_iter0overLast = histAxeRapRatio_iter0.Clone("histAxeRapRatio_iter0overLast")
    for i in range(1, histAxeRapRatio_iter0overLast.GetNbinsX() + 1):
        #histAxeRapRatio_iter0overLast.SetBinContent(i, ROOT.TMath.Abs(histAxeRapRatio_iter0.GetBinContent(i) - histAxeRapRatio_iterLast.GetBinContent(i))/ (histAxeRapRatio_iter0.GetBinContent(i) + histAxeRapRatio_iterLast.GetBinContent(i)))
        histAxeRapRatio_iter0overLast.SetBinContent(i, 0.5*ROOT.TMath.Abs(histAxeRapRatio_iter0.GetBinContent(i) - histAxeRapRatio_iterLast.GetBinContent(i))/ (histAxeRapRatio_iterLast.GetBinContent(i)))
        histAxeRapRatio_iter0overLast.SetBinError(i, 0)
        print(i, histAxeRapRatio_iter0overLast.GetBinContent(i))
    histAxeRapRatio_iter0overLast.SetTitle("(iteration0-iterationLast) / 2*iterationLast")
    histAxeRapRatio_iter0overLast.SetTitleSize(0.1, "t") 
    histAxeRapRatio_iter0overLast.GetYaxis().SetRangeUser(0, 0.025)
    histAxeRapRatio_iter0overLast.Draw("hp")

    canvasFinalSummary.SaveAs("FinalSummary.pdf")
    canvasFinalSummary.SaveAs("FinalSummary.root")

    ptMin = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0]
    ptMax = [0.5 ,1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0]

    rapMin = [2.50, 2.75, 3.00, 3.25, 3.50, 3.75]
    rapMax = [2.75, 3.00, 3.25, 3.50, 3.75, 4.00]

    #Print the Axe values for the ratio with given systematic uncertainty
    print("-------------- Ratio vs pT --------------")
    for i in range(1, histAxePtRatio_iter0overLast.GetNbinsX() +1):
        #axe_val = 0.5*(histAxePtRatio_iter0.GetBinContent(i) + histAxePtRatio_iterLast.GetBinContent(i))
        axe_val = histAxePtRatio_iterLast.GetBinContent(i)
        axe_sys = 0.5*ROOT.TMath.Abs(histAxePtRatio_iter0.GetBinContent(i) - histAxePtRatio_iterLast.GetBinContent(i))
        print(ptMin[i-1], ptMax[i-1], axe_val, axe_sys, axe_sys/axe_val)

    print("-------------- Ratio vs rapidity--------------")
    for i in range(1, histAxeRapRatio_iter0overLast.GetNbinsX() +1):
        #axe_val = 0.5*(histAxeRapRatio_iter0.GetBinContent(i) + histAxeRapRatio_iterLast.GetBinContent(i))
        axe_val = histAxeRapRatio_iterLast.GetBinContent(i)
        axe_sys = 0.5*ROOT.TMath.Abs(histAxeRapRatio_iter0.GetBinContent(i) - histAxeRapRatio_iterLast.GetBinContent(i))
        print(rapMin[i-1], rapMax[i-1], axe_val, axe_sys, axe_sys/axe_val)

    #Compute and compare integrated values vs pt and rapidity
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("-------------- Ratio vs Integral--------------")
    axePtRatio_int_iter0  = axePtPsi2S_int_iter0 / axePtJpsi_int_iter0
    axeRapRatio_int_iter0 = axeRapPsi2S_int_iter0 / axeRapJpsi_int_iter0

    axePtRatio_int_iterLast  = axePtPsi2S_int_iterLast / axePtJpsi_int_iterLast
    axeRapRatio_int_iterLast = axeRapPsi2S_int_iterLast / axeRapJpsi_int_iterLast

    print(f"iter 0, axe Jpsi, pt {axePtJpsi_int_iter0} and y {axeRapJpsi_int_iter0}")
    print(f"iter 0, axe Psi2s, pt {axePtPsi2S_int_iter0} and y {axeRapPsi2S_int_iter0}")

    print(f"iter L, axe Jpsi, pt {axePtJpsi_int_iterLast} and y {axeRapJpsi_int_iterLast}")
    print(f"iter L, axe Psi2s, pt {axePtPsi2S_int_iterLast} and y {axeRapPsi2S_int_iterLast}")

    print(f"iter 0, axe Psi2S/Jpsi, pt {axePtRatio_int_iter0} and y {axeRapRatio_int_iter0}")
    print(f"iter L, axe Psi2S/Jpsi, pt {axePtRatio_int_iterLast} and y {axeRapRatio_int_iterLast}")

    print("----- pol0 fits")
    fit_result_JpsiPt_0 = histAxePtJpsi_iter0.Fit("pol0", "S")
    p0JpsiPt_0 = fit_result_JpsiPt_0.Parameter(0)
    fit_result_JpsiRap_0 = histAxeRapJpsi_iter0.Fit("pol0", "S")
    p0JpsiRap_0 = fit_result_JpsiRap_0.Parameter(0)


    fit_result_JpsiPt_L = histAxePtJpsi_iterLast.Fit("pol0", "S")
    p0JpsiPt_L = fit_result_JpsiPt_L.Parameter(0)
    fit_result_JpsiRap_L = histAxeRapJpsi_iterLast.Fit("pol0", "S")
    p0JpsiRap_L = fit_result_JpsiRap_L.Parameter(0)

    fit_result_Psi2SPt_0 = histAxePtPsi2S_iter0.Fit("pol0", "S")
    p0Psi2SPt_0 = fit_result_Psi2SPt_0.Parameter(0)
    fit_result_Psi2SRap_0 = histAxeRapPsi2S_iter0.Fit("pol0", "S")
    p0Psi2SRap_0 = fit_result_Psi2SRap_0.Parameter(0)


    fit_result_Psi2SPt_L = histAxePtPsi2S_iterLast.Fit("pol0", "S")
    p0Psi2SPt_L = fit_result_Psi2SPt_L.Parameter(0)
    fit_result_Psi2SRap_L = histAxeRapPsi2S_iterLast.Fit("pol0", "S")
    p0Psi2SRap_L = fit_result_Psi2SRap_L.Parameter(0)

    print(f"iter 0, axe Jpsi, fit pol0, pt {p0JpsiPt_0} and y {p0JpsiRap_0}")
    print(f"iter L, axe Jpsi, fit pol0, pt {p0JpsiPt_L} and y {p0JpsiRap_L}")

    print(f"iter 0, axe Psi2S, fit pol0, pt {p0Psi2SPt_0} and y {p0Psi2SRap_0}")
    print(f"iter L, axe Psi2S, fit pol0, pt {p0Psi2SPt_L} and y {p0Psi2SRap_L}")

    print(f"iter 0, axe Psi2S/Jpsi, fit pol0, pt {p0Psi2SPt_0/p0JpsiPt_0} and y {p0Psi2SRap_0/p0JpsiRap_0}")
    print(f"iter L, axe Psi2S/Jpsi, fit pol0, pt {p0Psi2SPt_L/p0JpsiPt_L} and y {p0Psi2SRap_L/p0JpsiRap_L}")


    print("-----")
    print(f"iter VAL, axe Psi2S/Jpsi, pt {axePtRatio_int_iterLast} and y {axeRapRatio_int_iterLast}") #considering last iteration as tandard value
    print(f"iter SYS, axe Psi2S/Jpsi, pt {0.5*(axePtRatio_int_iter0-axePtRatio_int_iterLast)} and y {0.5*(axePtRatio_int_iter0-axePtRatio_int_iterLast)}")

    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

    print("###################################")
    print(r"""\begin{table}[h!]
        \begin{center}
            \begin{tabular}{cc}
                \hline
                \hline
                $\pt$ & $(A\times\varepsilon)_{\psip}$ / $(A\times\varepsilon)_{\jpsi}$ \\
                \hline""")
    for i in range(1, histAxePtRatio_iter0overLast.GetNbinsX() +1):
        #axe_val = 0.5*(histAxePtRatio_iter0.GetBinContent(i) + histAxePtRatio_iterLast.GetBinContent(i))
        #axe_stat = 0.5*(histAxePtRatio_iter0.GetBinError(i) + histAxePtRatio_iterLast.GetBinError(i))
        axe_val = histAxePtRatio_iterLast.GetBinContent(i)
        axe_stat = histAxePtRatio_iterLast.GetBinError(i)
        axe_sys = 0.5*ROOT.TMath.Abs(histAxePtRatio_iter0.GetBinContent(i) - histAxePtRatio_iterLast.GetBinContent(i))
        axe_stat_rel = axe_stat/axe_val *100
        axe_sys_rel = axe_sys/axe_val *100
        print(
            f"            {ptMin[i-1]} - {ptMax[i-1]} & {axe_val:.4f} $\pm$ {axe_stat:.4f} ({axe_stat_rel:.2f} \%) $\pm$ {axe_sys:.4f} ({axe_sys_rel:.2f} \%) \\\\"
        )
    print(r"""            \hline
            \end{tabular}
            \caption{Acceptance-efficency corrections as a function of $\pt$.}
            \label{tab:axe_ratio_pt_new}
        \end{center}
    \end{table}""")
    print("###################################")

    print("###################################")
    print(r"""\begin{table}[h!]
        \begin{center}
            \begin{tabular}{cc}
                \hline
                \hline
                $\y$ & $(A\times\varepsilon)_{\psip}$ / $(A\times\varepsilon)_{\jpsi}$ \\
                \hline""")
    for i in range(1, histAxeRapRatio_iter0overLast.GetNbinsX() +1):
        #axe_val = 0.5*(histAxeRapRatio_iter0.GetBinContent(i) + histAxeRapRatio_iterLast.GetBinContent(i))
        #axe_stat = 0.5*(histAxeRapRatio_iter0.GetBinError(i) + histAxeRapRatio_iterLast.GetBinError(i))
        axe_val = histAxeRapRatio_iterLast.GetBinContent(i)
        axe_stat = histAxeRapRatio_iterLast.GetBinError(i)
        axe_sys = 0.5*ROOT.TMath.Abs(histAxeRapRatio_iter0.GetBinContent(i) - histAxeRapRatio_iterLast.GetBinContent(i))
        axe_stat_rel = axe_stat/axe_val *100
        axe_sys_rel = axe_sys/axe_val *100
        print(
            f"            {rapMin[i-1]:.2f} - {rapMax[i-1]:.2f} & {axe_val:.4f} $\pm$ {axe_stat:.4f} ({axe_stat_rel:.2f} \%) $\pm$ {axe_sys:.4f} ({axe_sys_rel:.2f} \%) \\\\"
        )
    print(r"""            \hline
            \end{tabular}
            \caption{Acceptance-efficency corrections as a function of $\y$.}
            \label{tab:axe_ratio_rap_new}
        \end{center}
    \end{table}""")
    print("###################################")
    input()

if __name__ == "__main__":
    main()