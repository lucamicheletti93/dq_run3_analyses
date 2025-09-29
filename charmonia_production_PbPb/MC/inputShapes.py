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
    
kNsel: int = 51

selectionLabels = [
    "kIsBBV0A", "kIsBBV0C", "kIsBBFDA", "kIsBBFDC", "kIsBBT0A", "kIsBBT0C",
    "kNoBGV0A", "kNoBGV0C", "kNoBGFDA", "kNoBGFDC", "kNoBGT0A", "kNoBGT0C",
    "kIsBBZNA", "kIsBBZNC", "kIsBBZAC", "kNoBGZNA", "kNoBGZNC", "kNoV0MOnVsOfPileup",
    "kNoSPDOnVsOfPileup", "kNoV0Casymmetry", "kIsGoodTimeRange", "kNoIncompleteDAQ",
    "kNoTPCLaserWarmUp", "kNoTPCHVdip", "kNoPileupFromSPD", "kNoV0PFPileup", "kNoSPDClsVsTklBG",
    "kNoV0C012vsTklBG", "kNoInconsistentVtx", "kNoPileupInMultBins", "kNoPileupMV",
    "kNoPileupTPC", "kIsTriggerTVX", "kIsINT1", "kNoITSROFrameBorder", "kNoTimeFrameBorder",
    "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "kIsVertexITSTPC", "kIsVertexTOFmatched",
    "kIsVertexTRDmatched", "kNoCollInTimeRangeNarrow", "kNoCollInTimeRangeStrict",
    "kNoCollInTimeRangeStandard", "kNoCollInTimeRangeVzDependent", "kNoCollInRofStrict",
    "kNoCollInRofStandard", "kNoHighMultCollInPrevRof", "kIsGoodITSLayer3",
    "kIsGoodITSLayer0123", "kIsGoodITSLayersAll"
]

kTriggerMask = (
    (1 << 37) |  # kIsGoodZvtxFT0vsPV
    (1 << 32) |  # kIsTriggerTVX
    (1 << 36) |  # kNoSameBunchPileup
    (1 << 34) |  # kNoITSROFrameBorder
    (1 << 35)    # kNoTimeFrameBorder
)

def eventSelection(selection: int) -> bool:
    return (selection & kTriggerMask) == kTriggerMask

# ===================== #
#    Pt distribution    #
# ===================== #

def PtJPsiPbPb5TeV_Func():
    """
    J/psi pT in Pb—Pb
    """
    return ROOT.TF1(
        "PtPsiPbPb5TeV",
        lambda x, p: p[0] * x[0] / (1. + (x[0] / p[1])**p[2])**p[3],
        0, 20,
        4
    )

def PtJPsiPbPb5TeV_tuned(px):
    """
    J/psi pT in Pb—Pb, tuned on 2015 data -> Castillo embedding https://alice.its.cern.ch/jira/browse/ALIROOT-8174?jql=text%20~%20%22LHC19a2%22
    """
    x = px[0]
    p0 = 1.00715e6
    p1 = 3.50274
    p2 = 1.93403
    p3 = 3.96363
    return p0 * x / (1. + (x / p1)**p2)**p3


# ===================== #
# Rapidity distribution #
# ===================== #

def RapJPsiPbPb5TeV_Func():
    """
    TF1 for rap
    """
    return ROOT.TF1(
        "RapPsiPbPb5TeV",
        #lambda x, p: p[0] * x[0] / (1. + (x[0] / p[1])**p[2])**p[3], 2.5, 4, 4
        lambda x, p: p[0] * math.exp(-0.5 * ((x[0] - p[1]) / p[2])**2), 2.5, 4, 3
    )

def RapPsiPbPb5TeV_Original():
    """
    TF1 for J/psi rapidity in Pb—Pb (original)
    """
    return ROOT.TF1(
        "RapPsiPbPb5TeV",
        lambda x, p: p[0] * math.exp(-0.5 * ((x[0] - p[1]) / p[2])**2), 2.5, 4, 3
    )

def YPsiPbPb5TeV(py):
    """
    J/psi rapidity in Pb—Pb, tuned on 2015 data -> Castillo embedding https://alice.its.cern.ch/jira/browse/ALIROOT-8174?jql=text%20~%20%22LHC19a2%22
    Found at the link  https://github.com/AliceO2Group/O2DPG/blob/master/MC/config/PWGDQ/external/generator/GeneratorCocktailPromptCharmoniaToMuonEvtGen_PbPb5TeV.C#L164
    """
    y = py[0]
    p0 = 1.09886e6
    p1 = 0
    p2 = 2.12568
    return p0 * math.exp(-0.5 * ((y - p1) / p2)**2)

def YJPsipp5TeV_pp(py):
    """
    J/psi rapidity in pp at 5.02 TeV, tuned on HEPData INS1935680
    """
    y = py[0]
    p0 = 1
    p1 = 0.0338222
    p2 = 2.96748
    return p0 * math.exp(-0.5 * ((y - p1) / p2)**2)

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
    """
    creating an histo from systematic txt file
    """
    data = np.loadtxt(filename, skiprows=1) #skip first row
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
    #area = hist.Integral()
    area = hist.Integral("width")
    if area > 0:
        hist.Scale(1.0 / area, "width")

    return hist



def do_inputShapes(inputCfg):
    print("Start do_inputShapes")

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

    isPbPb = bool("false")
    database = ROOT.TDatabasePDG.Instance()
    muPdgCode = 13
    jpsiPdgCode = 443
    massMu = database.GetParticle(muPdgCode).Mass()
    massJpsi = database.GetParticle(jpsiPdgCode).Mass()

    """ fMass = fPt = fEta = fTauz = fTauxy = fU2Q2 = fCos2DeltaPhi = fR2EP = fR2SP = fCentFT0C = fImpactParameter = fMcDecision = -99999
    fChi2pca = fSVertex = fEMC1 = fEMC2 = fPt1 = fPt2 = fPhi1 = fPhi2 = fEta1 = fEta2 = fPtMC1 = fPtMC2 = fPhiMC1 = fPhiMC2 = fPhi = fEtaMC1 = fEtaMC2 = -99999
    fSign = fSign1 = fSign2 = -99999 """

    """ fMass = fPt = fEta = fTauz = fTauxy = fU2Q2 = fCos2DeltaPhi = fR2EP = fR2SP = fCentFT0C = fImpactParameter = fMcDecision = np.array([-99999], dtype=np.int32)
    fChi2pca = fSVertex = fEMC1 = fEMC2 = fPt1 = fPt2 = fPhi1 = fPhi2 = fEta1 = fEta2 = fPtMC1 = fPtMC2 = fPhiMC1 = fPhiMC2 = fPhi = fEtaMC1 = fEtaMC2 = np.array([-99999], dtype=np.int32)
    fSign = fSign1 = fSign2 = np.array([-99999], dtype=np.int32) """

    """ fMass = fPt = fEta = fTauz = fTauxy = fU2Q2 = fCos2DeltaPhi = fR2EP = fR2SP = fCentFT0C = fImpactParameter = np.array([-99999.0], dtype=np.float32)
    fChi2pca = fSVertex = fEMC1 = fEMC2 = fPt1 = fPt2 = fPhi1 = fPhi2 = fEta1 = fEta2 = fPtMC1 = fPtMC2 = fPhiMC1 = fPhiMC2 = fPhi = fEtaMC1 = fEtaMC2 = np.array([-99999.0], dtype=np.float32)
    fMcDecision = fSign = fSign1 = fSign2 = np.array([-99999.0], dtype=np.float32) """

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

    iterColors = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2]

    """ filein = ROOT.TFile(f"{pathName}/{fileName}", "READ")
    histPtJpsiData = filein.Get("histJpsiPt")
    SetHistogram(histPtJpsiData, ROOT.kBlack)
    histRapJpsiData = filein.Get("histJpsiRap")
    SetHistogram(histRapJpsiData, ROOT.kBlack)

    histPtJpsiData.Scale(1.0 / histPtJpsiData.Integral(), "WIDTH")
    histRapJpsiData.Scale(1.0 / histRapJpsiData.Integral(), "WIDTH") """

    histPtJpsiData = load_hist_from_txt(f"{pathFileDataPt}/{fileNameDataPt}", "histPtJpsiData", ";p_{T} (GeV/c);Normalized yield", ROOT.kBlack)
    histRapJpsiData = load_hist_from_txt(f"{pathFileDataRap}/{fileNameDataRap}", "histRapJpsiData", ";y;Normalized yield", ROOT.kBlack)
    if beamType == "PbPb":
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
    #fIn = ROOT.TFile(pathToFileAO2D/"AO2D-3.root", "READ")

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
        SetHistogram(histPtJpsiGen[iter], ROOT.kRed + 1)
        SetHistogram(histRapJpsiGen[iter], ROOT.kRed + 1)
        SetHistogram(histCentrJpsiGen[iter], ROOT.kRed + 1)

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
        for key in fIn.GetListOfKeys():
            """ if index >= 2:
                break """
            dirName = key.GetName()
            if "DF_" not in dirName:
                continue

            # ================ #
            #  Generated Tree  #
            # ================ #
            treeGen = fIn.Get(f"{dirName}/O2rtdilmtreegen")

            treeGen.SetBranchAddress("fMcDecision", fMcDecision)
            if beamType == "PbPb":
                treeGen.SetBranchAddress("fImpactParameter", fImpactParameter)
            treeGen.SetBranchAddress("fPtMC1", fPtMC1)
            treeGen.SetBranchAddress("fEtaMC1", fEtaMC1)
            treeGen.SetBranchAddress("fPhiMC1", fPhiMC1)
            treeGen.SetBranchAddress("fPtMC2", fPtMC2)
            treeGen.SetBranchAddress("fEtaMC2", fEtaMC2)
            treeGen.SetBranchAddress("fPhiMC2", fPhiMC2)
            print("dove mi blocco 1")
            treeGen.GetEntry(0)
            print("DEBUG treeGen first entry:")
            #print(f"fMcDecision={fMcDecision[0]}, fPtMC1={fPtMC1[0]}, fPtMC2={fPtMC2[0]}, fEtaMC1={fEtaMC1[0]}, fPhiMC1={fPhiMC1[0]}")

            """ massJpsi = 3.0969  # GeV/c^2
            pt = 2.0
            eta = 1.2
            phi = 0.5

            vec = ROOT.Math.PtEtaPhiMVector(pt, eta, phi, massJpsi)
            print(f"Rapidity: {vec.Rapidity()}, Pt: {vec.Pt()}, Eta: {vec.Eta()}") """

            for iEntry in range(treeGen.GetEntries()):
             #for iEntry in range(5000):
                treeGen.GetEntry(iEntry)
                if fMcDecision[0] < 1:
                    print(f"fMcDecision[0]: {fMcDecision[0]}")
                    continue
                vecJpsiGentest = ROOT.Math.PtEtaPhiMVector(fPtMC1[0], fEtaMC1[0], fPhiMC1[0], massJpsi)
                #print(f"vecJpsiGentest.Rapidity(): {vecJpsiGentest.Rapidity()}")
                    
                if fPtMC2[0] == -999.0 and fPhiMC2[0] == -999.0:
                    vecJpsiGen = ROOT.Math.PtEtaPhiMVector(fPtMC1[0], fEtaMC1[0], fPhiMC1[0], massJpsi)
                    #print(f"fPtMC1[0]={fPtMC1[0]}, fEtaMC1={fEtaMC1[0]}, fPhiMC1={fPhiMC1[0]}, massJpsi={massJpsi}")
                    #print(f"vecJpsiGen.Rapidity(): {vecJpsiGen.Rapidity()}")
                    if abs(vecJpsiGen.Rapidity()) > 4 or abs(vecJpsiGen.Rapidity()) < 2.5:
                        continue
                    if vecJpsiGen.Pt() > 20:
                        continue

                    if iter == 0:
                        histPtJpsiGen[iter].Fill(vecJpsiGen.Pt())
                        histRapJpsiGen[iter].Fill(-vecJpsiGen.Rapidity())
                        #print(f"vecJpsiGen.Pt()={vecJpsiGen.Pt()}, -vecJpsiGen.Rapidity()={-vecJpsiGen.Rapidity()}")
                    else:
                        binPt = histFromFuncPtRatio[iter-1].FindBin(vecJpsiGen.Pt())
                        binRap = histFromFuncRapRatio[iter-1].FindBin(-vecJpsiGen.Rapidity())
                        weightPt = histFromFuncPtRatio[iter-1].GetBinContent(binPt)
                        weightRap = histFromFuncRapRatio[iter-1].GetBinContent(binRap)
                        weightTot = weightPt * weightRap
                        histPtJpsiGen[iter].Fill(vecJpsiGen.Pt(), weightTot)
                        histRapJpsiGen[iter].Fill(-vecJpsiGen.Rapidity(), weightTot)

                    if beamType == "PbPb":
                        if fImpactParameter[0] < 5.625:
                            histCentrJpsiGen[iter].AddBinContent(1)
                        if 5.625 <= fImpactParameter[0] < 8.375:
                            histCentrJpsiGen[iter].AddBinContent(2)
                        if 8.375 <= fImpactParameter[0] < 10.625:
                            histCentrJpsiGen[iter].AddBinContent(3)
                        if 10.625 <= fImpactParameter[0] <= 13.875:
                            histCentrJpsiGen[iter].AddBinContent(4)

            """ for iEntry in range(min(20, treeGen.GetEntries())):
                treeGen.GetEntry(iEntry)
                print(f"Entry {iEntry}: fPtMC1={fPtMC1[0]}, fPtMC2={fPtMC2[0]}, Rapidity1={ROOT.Math.PtEtaPhiMVector(fPtMC1[0], fEtaMC1[0], fPhiMC1[0], massJpsi).Rapidity()}, Rapidity2={ROOT.Math.PtEtaPhiMVector(fPtMC2[0], fEtaMC2[0], fPhiMC2[0], massJpsi).Rapidity()}") """

            print("dove mi blocco 2")
            # ==================== #
            #  Reconstructed Tree  #
            # ==================== #
            treeRec = fIn.Get(f"{dirName}/O2rtdilmtreerec")

            treeRec.SetBranchAddress("fMcDecision", fMcDecision)
            treeRec.SetBranchAddress("fMass", fMass)
            treeRec.SetBranchAddress("fPt", fPt)
            treeRec.SetBranchAddress("fEta", fEta)
            treeRec.SetBranchAddress("fPhi", fPhi)
            treeRec.SetBranchAddress("fCentFT0C", fCentFT0C)
            treeRec.SetBranchAddress("fPtMC1", fPtMC1)
            treeRec.SetBranchAddress("fEtaMC1", fEtaMC1)
            treeRec.SetBranchAddress("fPhiMC1", fPhiMC1)
            treeRec.SetBranchAddress("fPtMC2", fPtMC2)
            treeRec.SetBranchAddress("fEtaMC2", fEtaMC2)
            treeRec.SetBranchAddress("fPhiMC2", fPhiMC2)
            treeRec.SetBranchAddress("fPt1", fPt1)
            treeRec.SetBranchAddress("fEta1", fEta1)
            treeRec.SetBranchAddress("fPhi1", fPhi1)
            treeRec.SetBranchAddress("fPt2", fPt2)
            treeRec.SetBranchAddress("fEta2", fEta2)
            treeRec.SetBranchAddress("fPhi2", fPhi2)

            
            treeRec.GetEntry(0)
            print("DEBUG treeRec first entry:")
            #print(f"fMcDecision={fMcDecision[0]}, fPtMC1={fPtMC1[0]}, fPtMC2={fPtMC2[0]}, fEtaMC1={fEtaMC1[0]}, fPhiMC1={fPhiMC1[0]}")

            print("dove mi blocco 3")
            for iEntry in range(treeRec.GetEntries()):
             #for iEntry in range(5000):
                treeRec.GetEntry(iEntry)

                if (fMcDecision[0] < 1):
                    print("continuo x fMCdec")
                    continue
                vecMuGen1 = ROOT.Math.PtEtaPhiMVector(fPtMC1[0], fEtaMC1[0], fPhiMC1[0], massMu)
                vecMuGen2 = ROOT.Math.PtEtaPhiMVector(fPtMC2[0], fEtaMC2[0], fPhiMC2[0], massMu)
                vecMuRec1 = ROOT.Math.PtEtaPhiMVector(fPt1[0], fEta1[0], fPhi1[0], massMu)
                vecMuRec2 = ROOT.Math.PtEtaPhiMVector(fPt2[0], fEta2[0], fPhi2[0], massMu)

                vecJpsiGen = vecMuGen1 + vecMuGen2
                vecJpsiRec = ROOT.Math.PtEtaPhiMVector(fPt[0], fEta[0], fPhi[0], fMass[0])

                if abs(vecJpsiRec.Rapidity()) > 4 or abs(vecJpsiRec.Rapidity()) < 2.5: 
                    print(f"continuo rec x rap: {vecJpsiRec.Rapidity()}")
                    continue
                if fPt1[0] < ptCut or fPt2[0] < ptCut: 
                    #print(f"continuo rec x ptCut, fPt1: {fPt1}  and fPt2: {fPt2}")
                    continue
                #if fPt[0] > 20 or fMass[0] > 3.4: 
                if fPt[0] > 20: 
                    #print(f"continuo rec x pt dimuon, fPt[0]: {fPt[0]}")
                    continue
                if(iter == 0):
                    histPtJpsiRec[iter].Fill(vecJpsiRec.Pt())
                    histRapJpsiRec[iter].Fill(-vecJpsiRec.Rapidity())
                    #print(f"vecJpsiRec.Pt()={vecJpsiRec.Pt()}, -vecJpsiRec.Rapidity()={-vecJpsiRec.Rapidity()}")
                else:
                    binPt = histFromFuncPtRatio[iter-1].FindBin(vecJpsiGen.Pt())
                    binRap = histFromFuncRapRatio[iter-1].FindBin(-vecJpsiGen.Rapidity())
                    weightPt = histFromFuncPtRatio[iter-1].GetBinContent(binPt)
                    weightRap = histFromFuncRapRatio[iter-1].GetBinContent(binRap)
                    weightTot = weightPt*weightRap
                    histPtJpsiRec[iter].Fill(vecJpsiRec.Pt(), weightTot)
                    #print(f"vecJpsiRec.Pt()={vecJpsiRec.Pt()}, weightTot={weightTot}")
                    histRapJpsiRec[iter].Fill(-vecJpsiRec.Rapidity(), weightTot)
                
                if beamType == "PbPb":
                    histCentrJpsiRec[iter].Fill(fCentFT0C[0])
                
                #print("all check rec ok")

            print("dove mi blocco 4")
            index += 1

        #-----------------------JPsi Axe histograms------------------------#

        histPtJpsiAxe[iter].Divide(histPtJpsiRec[iter], histPtJpsiGen[iter], 1, 1, "B")
        histRapJpsiAxe[iter].Divide(histRapJpsiRec[iter], histRapJpsiGen[iter], 1, 1, "B")
        print(f"Integral Gen Pt: {histPtJpsiGen[iter].Integral()}")
        print(f"Integral Rec Pt: {histPtJpsiRec[iter].Integral()}")

        print("Pt bin edges:")
        for ibin in range(nBinsPt+1):
            print(f"edge {ibin}: {ptBins[ibin]}")

        print("Rapidity bin edges:")
        for ibin in range(nBinsRap+1):
            print(f"edge {ibin}: {rapBins[ibin]}")

        print(f"Bins rec pT: {histPtJpsiRec[iter].GetNbinsX()} \n")
        print(f"Bins gen pT: {histPtJpsiGen[iter].GetNbinsX()} \n")
        print(f"Bins rec rap: {histRapJpsiRec[iter].GetNbinsX()} \n")
        print(f"Bins gen rap: {histRapJpsiGen[iter].GetNbinsX()}")
        for binAxe in range(histPtJpsiAxe[iter].GetNbinsX()):
            print(f"rec pT: {histPtJpsiRec[iter].GetBinContent(binAxe)}")
            print(f"gen pT: {histPtJpsiGen[iter].GetBinContent(binAxe)}")
            print(f"Axe pT: {histPtJpsiAxe[iter].GetBinContent(binAxe)}")
        for binAxe in range(histRapJpsiAxe[iter].GetNbinsX()):
            print(f"rec rap: {histRapJpsiRec[iter].GetBinContent(binAxe)}")
            print(f"gen rap: {histRapJpsiGen[iter].GetBinContent(binAxe)}")
            print(f"Axe rap: {histRapJpsiAxe[iter].GetBinContent(binAxe)}")
        if beamType == "PbPb":
            histCentrJpsiAxe[iter].Divide(histCentrJpsiRec[iter], histCentrJpsiGen[iter], 1, 1, "B")
            print(f"Bins rec centr: {histCentrJpsiRec[iter].GetNbinsX()} \n")
            print(f"Bins gen centr: {histCentrJpsiGen[iter].GetNbinsX()}")

        #-------------------------JPsi Axe ratio---------------------------#

        histPtJpsiDataCorr[iter].Divide(histPtJpsiAxe[iter])

        fitFunctionPt[iter] = PtJPsiPbPb5TeV_Func()
        fitFunctionPt[iter].SetParameters(1, 2.83941, 2.6687, 2.37032)
        fitFunctionPt[iter].SetLineColor(iterColors[iter])
        histPtJpsiDataCorr[iter].Fit(fitFunctionPt[iter], "R0") 

        histFromFuncPt[iter] = fitFunctionPt[iter].GetHistogram()
        histFromFuncPt[iter].SetName(f"histFromFuncPt_iter_{iter}")

        histFromFuncPtRatio[iter] = histFromFuncPt[iter].Clone(f"histFromFuncPtRatio_iter_{iter}")
        if (iter == 0):
            histFromFuncPtRatio[iter].Divide(histFromFuncPtOriginal)
        else:
            histFromFuncPtRatio[iter].Divide(histFromFuncPt[iter-1])
    
        histFromFuncPtRatio[iter].SetLineColor(iterColors[iter])

        histRapJpsiDataCorr[iter].Divide(histRapJpsiAxe[iter])

        fitFunctionRap[iter] = RapJPsiPbPb5TeV_Func()
        fitFunctionRap[iter].SetParameter(0, 8)       
        fitFunctionRap[iter].SetParameter(1, 0)       
        fitFunctionRap[iter].SetParameter(2, 2.12568) #now 3 parameters
        #fitFunctionRap[iter].SetParameters(1, 2.83941, 2.6687, 2.37032) #modified
        #fitFunctionRap[iter].SetParameters(1, 2.83941, 2.6687, 2.37032)
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


    c = ROOT.TCanvas("c","",800,600)
    histPtJpsiRec[0].SetLineColor(ROOT.kBlue)
    histPtJpsiGen[0].SetLineColor(ROOT.kRed)
    histPtJpsiRec[0].Draw("HIST")
    histPtJpsiGen[0].Draw("HIST SAME")
    c.SaveAs("debug_pt.pdf")

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