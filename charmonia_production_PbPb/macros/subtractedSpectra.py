import yaml
import json
import sys
import argparse
from array import array
import os
import math 
from os import path
import numpy as np
import pandas as pd
import ROOT
from ROOT import TCanvas, TH1F, TH2F

def rebin_centr_pT(pt_bins, centr_bins):
    input_file_path = "/Users/saragaretti/cernbox/JPsi_Psi2S_ratio/Train_321538/data/LHC23_pass4_pTcut_1GeV_new/mergedHistograms.root"
    output_file_path = "/Users/saragaretti/cernbox/JPsi_Psi2S_ratio/Train_321538/data/LHC23_pass4_pTcut_1GeV_new/mergedHistograms_sub.root"

    # Open input file
    f_in = ROOT.TFile.Open(input_file_path, "READ")
    if not f_in or f_in.IsZombie():
        print("Error: Unable to open the source file!")
        exit()

    # Create output file
    f_out = ROOT.TFile(output_file_path, "RECREATE")
    index = 0

    hist_Int_SEPM = f_in.Get(f"Mass_Int_SEPM").Clone()
    hist_Int_SEPP = f_in.Get(f"Mass_Int_SEPP").Clone()
    hist_Int_SEMM = f_in.Get(f"Mass_Int_SEMM").Clone()
    hist_Int_MEPM = f_in.Get(f"Mass_Int_MEPM").Clone()
    hist_Int_MEPP = f_in.Get(f"Mass_Int_MEPP").Clone()
    hist_Int_MEMM = f_in.Get(f"Mass_Int_MEMM").Clone()

    hist_Pt_SEPM = f_in.Get(f"Mass_Pt_SEPM").Clone()
    hist_Pt_SEPP = f_in.Get(f"Mass_Pt_SEPP").Clone()
    hist_Pt_SEMM = f_in.Get(f"Mass_Pt_SEMM").Clone()
    hist_Pt_MEPM = f_in.Get(f"Mass_Pt_MEPM").Clone()
    hist_Pt_MEPP = f_in.Get(f"Mass_Pt_MEPP").Clone()
    hist_Pt_MEMM = f_in.Get(f"Mass_Pt_MEMM").Clone()

    hist_Rap_SEPM = f_in.Get(f"Mass_Rap_SEPM").Clone()
    hist_Rap_SEPP = f_in.Get(f"Mass_Rap_SEPP").Clone()
    hist_Rap_SEMM = f_in.Get(f"Mass_Rap_SEMM").Clone()
    hist_Rap_MEPM = f_in.Get(f"Mass_Rap_MEPM").Clone()
    hist_Rap_MEPP = f_in.Get(f"Mass_Rap_MEPP").Clone()
    hist_Rap_MEMM = f_in.Get(f"Mass_Rap_MEMM").Clone()

    # Initialize variables for normalization
    n_MEPM_Int = n_MEPP_Int = n_MEMM_Int = n_SEPP_Int = n_SEMM_Int = n_SEPM_Int = 0
    #sum_SEPM = sum_SEPP = sum_SEMM = sum_MEPM = sum_MEPP = sum_MEMM = 0
    err_SEPP_Int = err_SEMM_Int = err_MEPM_Int = err_MEPP_Int = err_MEMM_Int = 0
    R_bin_Int = err_R_bin_Int = 0
    F_Int = R_Int = F_num_Int = F_denom_Int = 0
    bin_temp_Int = geom_avg_Int = geom_err_Int = 0

    # Calculate F factor
    hist_R_values_Int = ROOT.TH1D("hist_R_values_centr", "Valori di R bin per bin", hist_Int_SEPM.GetNbinsX(), hist_Int_SEPM.GetXaxis().GetXmin(), hist_Int_SEPM.GetXaxis().GetXmax())
    hist_LS_Int = hist_Int_MEPM.Clone(f"Mass_LS_Int")
    hist_LS_norm_Int = hist_Int_MEPM.Clone(f"Mass_LS_Int")
    for nBin in range(1, hist_Int_MEPM.GetNbinsX() + 1):

        n_SEPM_Int = hist_Int_SEPM.GetBinContent(nBin) 
        n_SEPP_Int = hist_Int_SEPP.GetBinContent(nBin)
        n_SEMM_Int = hist_Int_SEMM.GetBinContent(nBin)
        n_MEPM_Int = hist_Int_MEPM.GetBinContent(nBin)
        n_MEPP_Int = hist_Int_MEPP.GetBinContent(nBin)
        n_MEMM_Int = hist_Int_MEMM.GetBinContent(nBin)

        err_SEMM_Int = hist_Int_SEMM.GetBinError(nBin)
        err_SEPP_Int = hist_Int_SEPP.GetBinError(nBin)
        err_MEPM_Int = hist_Int_MEPM.GetBinError(nBin)
        err_MEPP_Int = hist_Int_MEPP.GetBinError(nBin)
        err_MEMM_Int = hist_Int_MEMM.GetBinError(nBin)

        R_bin_Int = n_MEPM_Int / (2 * math.sqrt(n_MEPP_Int * n_MEMM_Int))
        deriv_MEPM_Int = 1 / (2 * math.sqrt(hist_Int_MEPP.GetBinContent(nBin) * hist_Int_MEMM.GetBinContent(nBin)))
        deriv_MEPP_Int = - hist_Int_MEPM.GetBinContent(nBin) / (4 * hist_Int_MEPP.GetBinContent(nBin) ** (3 / 2) * math.sqrt(hist_Int_MEMM.GetBinContent(nBin)))
        deriv_MEMM_Int = - hist_Int_MEPM.GetBinContent(nBin) / (4 * hist_Int_MEMM.GetBinContent(nBin) ** (3 / 2) * math.sqrt(hist_Int_MEPP.GetBinContent(nBin)))
        err_R_bin = math.sqrt((err_MEPM_Int *deriv_MEPM_Int) ** 2 +(err_MEPP_Int * deriv_MEPP_Int) ** 2 +(err_MEMM_Int * deriv_MEMM_Int) ** 2)
        hist_R_values_Int.SetBinContent(nBin, R_bin_Int)
        hist_R_values_Int.SetBinError(nBin, err_R_bin_Int)
        R_Int += R_bin_Int

        F_num_Int += 2 * R_bin_Int * math.sqrt(n_SEPP_Int * n_SEMM_Int)
        F_denom_Int += n_MEPM_Int

        geom_avg_Int = 2 * math.sqrt(n_SEMM_Int * n_SEPP_Int)
        geom_err_Int = math.sqrt((n_SEPP_Int / n_SEMM_Int) * err_SEMM_Int**2 + (n_SEMM_Int / n_SEPP_Int) * err_SEPP_Int**2)
        LS_bin_err_Int = math.sqrt((geom_avg_Int*err_R_bin_Int) ** 2 + (R_bin_Int*geom_err_Int) ** 2)
        hist_LS_Int.SetBinContent(nBin, geom_avg_Int)
        hist_LS_Int.SetBinError(nBin, geom_err_Int)
        bin_temp_Int = R_bin_Int*geom_avg_Int
        hist_LS_norm_Int.SetBinContent(nBin, bin_temp_Int)
        hist_LS_norm_Int.SetBinError(nBin, LS_bin_err_Int)
        #print(f"bin histo LS centr {nBin}: R = {R_bin_Int:.4f}, geom = {geom_avg_Int:.4f}, R*geom_avg = {bin_temp_Int:.4f}")

    F_Int = F_num_Int / F_denom_Int
    print(f"integrated F = {F_Int:.4f}")
    # Normalize MEPM histogram
    hist_Int_MEPM.Scale(F_Int)
    hist_Pt_MEPM.Scale(F_Int)
    hist_Rap_MEPM.Scale(F_Int)
    hist_Int_MEPP.Scale(F_Int)
    hist_Pt_MEPP.Scale(F_Int)
    hist_Rap_MEPP.Scale(F_Int)
    hist_Int_MEMM.Scale(F_Int)
    hist_Pt_MEMM.Scale(F_Int)
    hist_Rap_MEMM.Scale(F_Int)

    # Ratio calculation integrated
    histRatioSumCentrPM = hist_Int_SEPM.Clone(f"ratio_Int_PM")
    histRatioSumCentrPM.Divide(hist_Int_MEPM)
    histRatioSumCentrPP = hist_Int_SEPP.Clone(f"ratio_Int_PP")
    histRatioSumCentrPP.Divide(hist_Int_MEPP)
    histRatioSumCentrMM = hist_Int_SEMM.Clone(f"ratio_Int_MM")
    histRatioSumCentrMM.Divide(hist_Int_MEMM)
    # Ratio calculation pT
    histRatioSumPtPM = hist_Pt_SEPM.Clone(f"ratio_Pt_PM")
    histRatioSumPtPM.Divide(hist_Pt_MEPM)
    histRatioSumPtPP = hist_Pt_SEPP.Clone(f"ratio_Pt_PP")
    histRatioSumPtPP.Divide(hist_Pt_MEPP)
    histRatioSumPtMM = hist_Pt_SEMM.Clone(f"ratio_Pt_MM")
    histRatioSumPtMM.Divide(hist_Pt_MEMM)
    # Ratio calculation rapidity
    histRatioSumRapPM = hist_Rap_SEPM.Clone(f"ratio_Rap_PM")
    histRatioSumRapPM.Divide(hist_Rap_MEPM)
    histRatioSumRapPP = hist_Rap_SEPP.Clone(f"ratio_Rap_PP")
    histRatioSumRapPP.Divide(hist_Rap_MEPP)
    histRatioSumRapMM = hist_Rap_SEMM.Clone(f"ratio_Rap_MM")
    histRatioSumRapMM.Divide(hist_Rap_MEMM)

    #histograms lists initialization
    hist_SEPM_list = []
    hist_SEPP_list = []
    hist_SEMM_list = []
    hist_MEPM_list = []
    hist_MEPP_list = []
    hist_MEMM_list = []
    hist_SEPM_Pt_list = []
    hist_SEPP_Pt_list = []
    hist_SEMM_Pt_list = []
    hist_MEPM_Pt_list = []
    hist_MEPP_Pt_list = []
    hist_MEMM_Pt_list = []
    hist_SEPM_Rap_list = []
    hist_SEPP_Rap_list = []
    hist_SEMM_Rap_list = []
    hist_MEPM_Rap_list = []
    hist_MEPP_Rap_list = []
    hist_MEMM_Rap_list = []

    #ratio histograms lists initialization
    histRatioCentrPM = [None]*10
    histRatioCentrPP = [None]*10
    histRatioCentrMM = [None]*10
    histRatioPtPM = [None]*10
    histRatioPtPP = [None]*10
    histRatioPtMM = [None]*10
    histRatioRapPM = [None]*10
    histRatioRapPP = [None]*10
    histRatioRapMM = [None]*10

    #Racc histos list fot centrality
    hist_R_valuesCentr_list = []
    hist_R_valuesCentr = [None]*10
    lineCentr_list = []
    lineCentr_index = [None]*10
    #Racc histos list fot pT
    hist_R_valuesPt_list = []
    hist_R_valuesPt = [None]*10
    linePt_list = []
    linePt_index = [None]*10

    # Loop over pT bins defined manually
    # Create a canvas to contain all hist_R_values histograms
    canvaspT = ROOT.TCanvas("canvas_R_allPt", "All pT R Values", 1200, 800)
    canvaspT.Divide(3, 2)  # Divide the canvas into a 3x2 grid (6 pT bins in total)
    # Initialize a counter for the pads in the canvas
    pad_index_pT = 1
    for pt_min, pt_max in pt_bins:
        pt_bin = f"{pt_min}_{pt_max}"
        print(f"Processing pT bin: {pt_bin}")

        # Load histograms
        hist_SEPM = f_in.Get(f"Mass_Pt_{pt_bin}_SEPM").Clone()
        hist_SEPP = f_in.Get(f"Mass_Pt_{pt_bin}_SEPP").Clone()
        hist_SEMM = f_in.Get(f"Mass_Pt_{pt_bin}_SEMM").Clone()
        hist_MEPM = f_in.Get(f"Mass_Pt_{pt_bin}_MEPM").Clone()
        hist_MEPP = f_in.Get(f"Mass_Pt_{pt_bin}_MEPP").Clone()
        hist_MEMM = f_in.Get(f"Mass_Pt_{pt_bin}_MEMM").Clone()

        # Initialize variables for normalization
        n_MEPM = n_MEPP = n_MEMM = n_SEPP = n_SEMM = n_SEPM = 0
        #sum_SEPM = sum_SEPP = sum_SEMM = sum_MEPM = sum_MEPP = sum_MEMM = 0
        err_SEPP = err_SEMM = err_MEPM = err_MEPP = err_MEMM = 0
        F_num = F_denom = 0
        R_bin = err_R_bin = 0
        F = R = F_num = F_denom = 0
        bin_temp = geom_avg = geom_err = 0

        # Calculate F factor
        # Calculate R values bin by bin
        hist_R_values = ROOT.TH1D(f"hist_R_values_{pt_bin}", f"R Values for pT bin {pt_bin}",hist_SEPM.GetNbinsX(), hist_SEPM.GetXaxis().GetXmin(), hist_SEPM.GetXaxis().GetXmax())
        hist_LS = hist_MEPM.Clone(f"Mass_LS_{pt_bin}")
        hist_LS_norm = hist_MEPM.Clone(f"Mass_LS_{pt_bin}")
        for nBin in range(1, hist_MEPM.GetNbinsX() + 1):
            n_SEPM = hist_SEPM.GetBinContent(nBin) 
            n_SEPP = hist_SEPP.GetBinContent(nBin)
            n_SEMM = hist_SEMM.GetBinContent(nBin)
            n_MEPM = hist_MEPM.GetBinContent(nBin)
            n_MEPP = hist_MEPP.GetBinContent(nBin)
            n_MEMM = hist_MEMM.GetBinContent(nBin)

            err_SEMM = hist_SEMM.GetBinError(nBin)
            err_SEPP = hist_SEPP.GetBinError(nBin)
            err_MEPM = hist_MEPM.GetBinError(nBin)
            err_MEPP = hist_MEPP.GetBinError(nBin)
            err_MEMM = hist_MEMM.GetBinError(nBin)

            R_bin = n_MEPM / (2 * math.sqrt(n_MEPP * n_MEMM))
            deriv_MEPM = 1 / (2 * math.sqrt(hist_MEPP.GetBinContent(nBin) * hist_MEMM.GetBinContent(nBin)))
            deriv_MEPP = - hist_MEPM.GetBinContent(nBin) / (4 * hist_MEPP.GetBinContent(nBin) ** (3 / 2) * math.sqrt(hist_MEMM.GetBinContent(nBin)))
            deriv_MEMM = - hist_MEPM.GetBinContent(nBin) / (4 * hist_MEMM.GetBinContent(nBin) ** (3 / 2) * math.sqrt(hist_MEPP.GetBinContent(nBin)))
            err_R_bin = math.sqrt((err_MEPM *deriv_MEPM) ** 2 +(err_MEPP * deriv_MEPP) ** 2 +(err_MEMM * deriv_MEMM) ** 2)
            hist_R_values.SetBinContent(nBin, R_bin)
            hist_R_values.SetBinError(nBin, err_R_bin)
            R += R_bin

            F_num += 2 * R_bin * math.sqrt(n_SEPP * n_SEMM)
            F_denom += n_MEPM

            geom_avg = 2 * math.sqrt(n_SEMM * n_SEPP)
            geom_err = math.sqrt((n_SEPP / n_SEMM) * err_SEMM**2 + (n_SEMM / n_SEPP) * err_SEPP**2)
            LS_bin_err = math.sqrt((geom_avg*err_R_bin) ** 2 + (R_bin*geom_err) ** 2)
            hist_LS.SetBinContent(nBin, geom_avg)
            hist_LS.SetBinError(nBin, geom_err)
            bin_temp = R_bin*geom_avg
            hist_LS_norm.SetBinContent(nBin, bin_temp)
            hist_LS_norm.SetBinError(nBin, LS_bin_err)
            #print(f"bin histo LS {nBin}: R = {R_bin:.4f}, geom = {geom_avg:.4f}, R*geom_avg = {bin_temp:.4f}")

        F = F_num / F_denom
        print(f"pT bin {pt_bin}: F = {F:.4f}")

        # Normalize MEPM histogram
        hist_MEPM.Scale(F)

        # Background subtraction ME
        hist_BkgSub_ME = hist_SEPM.Clone(f"Mass_BkgSub_{pt_bin}_ME")
        hist_BkgSub_ME.Add(hist_MEPM, -1)

        # Calculate Like-Sign Average and Background Subtraction LS

        hist_BkgSub_LS = hist_SEPM.Clone(f"Mass_BkgSub_{pt_bin}_LS")
        hist_BkgSub_LS.Add(hist_LS, -1)

        hist_BkgSub_LS_norm = hist_SEPM.Clone(f"Mass_BkgSub_{pt_bin}_LS_norm")
        hist_BkgSub_LS_norm.Add(hist_LS_norm, -1)

        # Draw the histogram on the divided canvas
        hist_R_valuesPt_list.append(hist_R_values)
        hist_R_valuesPt[index] = hist_R_valuesPt_list[index].Clone(f"hist_R_values_{pt_bin}")
        canvaspT.cd(index+1)
        hist_R_valuesPt[index].SetMarkerStyle(20)
        hist_R_valuesPt[index].SetMarkerSize(0.8)
        hist_R_valuesPt[index].SetLineColor(ROOT.kBlue)
        hist_R_valuesPt[index].SetLineWidth(2)
        hist_R_valuesPt[index].GetYaxis().SetRangeUser(0.9,1.1)
        hist_R_valuesPt[index].Draw("EP")

        # Draw a dashed line at y = 1
        linePt = ROOT.TLine(hist_R_valuesPt[index].GetXaxis().GetXmin(), 1, hist_R_valuesPt[index].GetXaxis().GetXmax(), 1)
        linePt.SetLineColor(ROOT.kRed)
        linePt.SetLineStyle(2)
        linePt.SetLineWidth(2)
        linePt_list.append(linePt)
        linePt_index[index] = linePt_list[index].Clone()
        linePt_index[index].Draw("SAME")
        canvaspT.Update()

        pad_index_pT += 1
        print("pad_index_pT: ", pad_index_pT)

        # Save histograms to the output file
        f_out.cd()
        hist_SEPM.Write(f"Mass_Pt_{pt_bin}_SEPM")
        hist_MEPM.Write(f"Mass_Pt_{pt_bin}_MEPM_normalized")
        hist_LS.Write(f"Mass_LS_{pt_bin}")
        hist_BkgSub_ME.Write(f"Mass_BkgSub_{pt_bin}_ME")
        hist_BkgSub_LS.Write(f"Mass_BkgSub_{pt_bin}_LS")
        hist_BkgSub_LS_norm.Write(f"Mass_BkgSub_{pt_bin}_LS_norm")
        hist_R_values.Write()
        index = index + 1

    # Create a canvas to contain all hist_R_values histograms
    canvasCentr = ROOT.TCanvas("canvas_R_allCentr", "All Centr R Values", 1200, 800)
    canvasCentr.Divide(2, 2)  # Divide the canvas into a 2x2 grid (4 centr bins in total)
    # Initialize a counter for the pads in the canvas
    pad_index_centr = 1
    index = 0
    for centr_min, centr_max in centr_bins:
        centr_bins = f"{centr_min}_{centr_max}"
        print(f"Processing pT bin: {centr_bins}")

        # Load histograms
        hist_SEPM = f_in.Get(f"Mass_CentrFT0C_{centr_bins}_SEPM").Clone()
        hist_SEPP = f_in.Get(f"Mass_CentrFT0C_{centr_bins}_SEPP").Clone()
        hist_SEMM = f_in.Get(f"Mass_CentrFT0C_{centr_bins}_SEMM").Clone()
        hist_MEPM = f_in.Get(f"Mass_CentrFT0C_{centr_bins}_MEPM").Clone()
        hist_MEPP = f_in.Get(f"Mass_CentrFT0C_{centr_bins}_MEPP").Clone()
        hist_MEMM = f_in.Get(f"Mass_CentrFT0C_{centr_bins}_MEMM").Clone()
        # Load histograms pT
        hist_Pt_SEPM = f_in.Get(f"Mass_Pt_{centr_bins}_SEPM").Clone()
        hist_Pt_SEPP = f_in.Get(f"Mass_Pt_{centr_bins}_SEPP").Clone()
        hist_Pt_SEMM = f_in.Get(f"Mass_Pt_{centr_bins}_SEMM").Clone()
        hist_Pt_MEPM = f_in.Get(f"Mass_Pt_{centr_bins}_MEPM").Clone()
        hist_Pt_MEPP = f_in.Get(f"Mass_Pt_{centr_bins}_MEPP").Clone()
        hist_Pt_MEMM = f_in.Get(f"Mass_Pt_{centr_bins}_MEMM").Clone()
        # Load histograms rapidity
        hist_Rap_SEPM = f_in.Get(f"Mass_Rap_{centr_bins}_SEPM").Clone()
        hist_Rap_SEPP = f_in.Get(f"Mass_Rap_{centr_bins}_SEPP").Clone()
        hist_Rap_SEMM = f_in.Get(f"Mass_Rap_{centr_bins}_SEMM").Clone()
        hist_Rap_MEPM = f_in.Get(f"Mass_Rap_{centr_bins}_MEPM").Clone()
        hist_Rap_MEPP = f_in.Get(f"Mass_Rap_{centr_bins}_MEPP").Clone()
        hist_Rap_MEMM = f_in.Get(f"Mass_Rap_{centr_bins}_MEMM").Clone()

        # Initialize variables for normalization
        n_MEPM = n_MEPP = n_MEMM = n_SEPP = n_SEMM = n_SEPM = 0
        #sum_SEPM = sum_SEPP = sum_SEMM = sum_MEPM = sum_MEPP = sum_MEMM = 0
        err_SEPP = err_SEMM = err_MEPM = err_MEPP = err_MEMM = 0
        F_num = F_denom = 0
        R_bin = err_R_bin = 0
        F = R = F_num = F_denom = 0
        bin_temp = geom_avg = geom_err = 0

        # Calculate F factor
        #hist_R_values = ROOT.TH1D("hist_R_values_centr", "Valori di R bin per bin", hist_SEPM.GetNbinsX(), hist_SEPM.GetXaxis().GetXmin(), hist_SEPM.GetXaxis().GetXmax())
        hist_R_values = ROOT.TH1D(f"hist_R_values_{centr_bins}", f"R Values for Centr bin {centr_bins}",hist_SEPM.GetNbinsX(), hist_SEPM.GetXaxis().GetXmin(), hist_SEPM.GetXaxis().GetXmax())
        hist_LS = hist_MEPM.Clone(f"Mass_LS_{centr_bins}")
        hist_LS_norm = hist_MEPM.Clone(f"Mass_LS_{centr_bins}")
        for nBin in range(1, hist_MEPM.GetNbinsX() + 1):
            n_SEPM = hist_SEPM.GetBinContent(nBin) 
            n_SEPP = hist_SEPP.GetBinContent(nBin)
            n_SEMM = hist_SEMM.GetBinContent(nBin)
            n_MEPM = hist_MEPM.GetBinContent(nBin)
            n_MEPP = hist_MEPP.GetBinContent(nBin)
            n_MEMM = hist_MEMM.GetBinContent(nBin)

            err_SEMM = hist_SEMM.GetBinError(nBin)
            err_SEPP = hist_SEPP.GetBinError(nBin)
            err_MEPM = hist_MEPM.GetBinError(nBin)
            err_MEPP = hist_MEPP.GetBinError(nBin)
            err_MEMM = hist_MEMM.GetBinError(nBin)

            R_bin = n_MEPM / (2 * math.sqrt(n_MEPP * n_MEMM))
            deriv_MEPM = 1 / (2 * math.sqrt(hist_MEPP.GetBinContent(nBin) * hist_MEMM.GetBinContent(nBin)))
            deriv_MEPP = - hist_MEPM.GetBinContent(nBin) / (4 * hist_MEPP.GetBinContent(nBin) ** (3 / 2) * math.sqrt(hist_MEMM.GetBinContent(nBin)))
            deriv_MEMM = - hist_MEPM.GetBinContent(nBin) / (4 * hist_MEMM.GetBinContent(nBin) ** (3 / 2) * math.sqrt(hist_MEPP.GetBinContent(nBin)))
            err_R_bin = math.sqrt((err_MEPM *deriv_MEPM) ** 2 +(err_MEPP * deriv_MEPP) ** 2 +(err_MEMM * deriv_MEMM) ** 2)
            hist_R_values.SetBinContent(nBin, R_bin)
            hist_R_values.SetBinError(nBin, err_R_bin)
            R += R_bin

            F_num += 2 * R_bin * math.sqrt(n_SEPP * n_SEMM)
            F_denom += n_MEPM

            geom_avg = 2 * math.sqrt(n_SEMM * n_SEPP)
            geom_err = math.sqrt((n_SEPP / n_SEMM) * err_SEMM**2 + (n_SEMM / n_SEPP) * err_SEPP**2)
            LS_bin_err = math.sqrt((geom_avg*err_R_bin) ** 2 + (R_bin*geom_err) ** 2)
            hist_LS.SetBinContent(nBin, geom_avg)
            hist_LS.SetBinError(nBin, geom_err)
            bin_temp = R_bin*geom_avg
            hist_LS_norm.SetBinContent(nBin, bin_temp)
            hist_LS_norm.SetBinError(nBin, LS_bin_err)
            #print(f"bin histo LS centr {nBin}: R = {R_bin:.4f}, geom = {geom_avg:.4f}, R*geom_avg = {bin_temp:.4f}")

        F = F_num / F_denom
        print(f"centr bin {centr_bins}: F = {F:.4f}")

        # Normalize MEPM histogram
        hist_MEPM.Scale(F)
        hist_Pt_MEPM.Scale(F)
        hist_Rap_MEPM.Scale(F)
        hist_MEPP.Scale(F)
        hist_Pt_MEPP.Scale(F)
        hist_Rap_MEPP.Scale(F)
        hist_MEMM.Scale(F)
        hist_Pt_MEMM.Scale(F)
        hist_Rap_MEMM.Scale(F)

        # Background subtraction ME
        hist_BkgSub_ME = hist_SEPM.Clone(f"Mass_BkgSub_{centr_bins}_ME")
        hist_BkgSub_ME.Add(hist_MEPM, -1)

        # Calculate Like-Sign Average and Background Subtraction LS

        hist_BkgSub_LS = hist_SEPM.Clone(f"Mass_BkgSub_{centr_bins}_LS")
        hist_BkgSub_LS.Add(hist_LS, -1)

        hist_BkgSub_LS_norm = hist_SEPM.Clone(f"Mass_BkgSub_{centr_bins}_LS_norm")
        hist_BkgSub_LS_norm.Add(hist_LS_norm, -1)

        # Draw the histogram on the divided canvas
        canvasCentr.cd(index+1)
        hist_R_values.SetMarkerStyle(20)
        hist_R_values.SetMarkerSize(0.8)
        hist_R_values.SetLineColor(ROOT.kBlue)
        hist_R_values.SetLineWidth(2)
        hist_R_valuesCentr_list.append(hist_R_values)
        hist_R_valuesCentr[index] = hist_R_valuesCentr_list[index].Clone(f"hist_R_values_{centr_bins}")
        hist_R_valuesCentr[index].GetYaxis().SetRangeUser(0.9,1.1)
        hist_R_valuesCentr[index].Draw("EP")

        # Draw a dashed line at y = 1
        lineCentr = ROOT.TLine(hist_R_values.GetXaxis().GetXmin(), 1, hist_R_values.GetXaxis().GetXmax(), 1)
        lineCentr.SetLineColor(ROOT.kRed)
        lineCentr.SetLineStyle(2)
        lineCentr.SetLineWidth(2)
        lineCentr_list.append(lineCentr)
        lineCentr_index[index] = lineCentr_list[index].Clone()
        lineCentr_index[index].Draw("SAME")
        canvasCentr.Update()

        pad_index_centr += 1
        print("pad_index_centr: ", pad_index_centr)

        # Ratio calculation
        hist_SEPM_list.append(hist_SEPM)
        hist_MEPM_list.append(hist_MEPM)
        histRatioCentrPM[index] = hist_SEPM_list[index].Clone(f"ratio_PM_{centr_bins}")
        histRatioCentrPM[index].Divide(hist_MEPM_list[index])
        hist_SEPP_list.append(hist_SEPP)
        hist_MEPP_list.append(hist_MEPP)
        histRatioCentrPP[index] = hist_SEPP_list[index].Clone(f"ratio_PP_{centr_bins}")
        histRatioCentrPP[index].Divide(hist_MEPP_list[index])
        hist_SEMM_list.append(hist_SEMM)
        hist_MEMM_list.append(hist_MEMM)
        histRatioCentrMM[index] = hist_SEMM_list[index].Clone(f"ratio_MM_{centr_bins}")
        histRatioCentrMM[index].Divide(hist_MEMM_list[index])
        # Ratio calculation pT
        hist_SEPM_Pt_list.append(hist_Pt_SEPM)
        hist_MEPM_Pt_list.append(hist_Pt_MEPM)
        histRatioPtPM[index] = hist_SEPM_Pt_list[index].Clone(f"ratio_Pt_PM_{centr_bins}")
        histRatioPtPM[index].Divide(hist_MEPM_Pt_list[index])
        hist_SEPP_Pt_list.append(hist_Pt_SEPP)
        hist_MEPP_Pt_list.append(hist_Pt_MEPP)
        histRatioPtPP[index] = hist_SEPP_Pt_list[index].Clone(f"ratio_Pt_PP_{centr_bins}")
        histRatioPtPP[index].Divide(hist_MEPP_Pt_list[index])
        hist_SEMM_Pt_list.append(hist_Pt_SEMM)
        hist_MEMM_Pt_list.append(hist_Pt_MEMM)
        histRatioPtMM[index] = hist_SEMM_Pt_list[index].Clone(f"ratio_Pt_MM_{centr_bins}")
        histRatioPtMM[index].Divide(hist_MEMM_Pt_list[index])
        # Ratio calculation rapidity
        hist_SEPM_Rap_list.append(hist_Rap_SEPM)
        hist_MEPM_Rap_list.append(hist_Rap_MEPM)
        histRatioRapPM[index] = hist_SEPM_Rap_list[index].Clone(f"ratio_Rap_PM_{centr_bins}")
        histRatioRapPM[index].Divide(hist_MEPM_Rap_list[index])
        hist_SEPP_Rap_list.append(hist_Rap_SEPP)
        hist_MEPP_Rap_list.append(hist_Rap_MEPP)
        histRatioRapPP[index] = hist_SEPP_Rap_list[index].Clone(f"ratio_Rap_PP_{centr_bins}")
        histRatioRapPP[index].Divide(hist_MEPP_Rap_list[index])
        hist_SEMM_Rap_list.append(hist_Rap_SEMM)
        hist_MEMM_Rap_list.append(hist_Rap_MEMM)
        histRatioRapMM[index] = hist_SEMM_Rap_list[index].Clone(f"ratio_Rap_MM_{centr_bins}")
        histRatioRapMM[index].Divide(hist_MEMM_Rap_list[index])

        # Save histograms to the output file
        f_out.cd()
        hist_SEPM.Write(f"Mass_CentrFT0C_{centr_bins}_SEPM")
        hist_MEPM.Write(f"Mass_CentrFT0C_{centr_bins}_MEPM_normalized")
        hist_LS.Write(f"Mass_LS_{centr_bins}")
        hist_BkgSub_ME.Write(f"Mass_BkgSub_{centr_bins}_ME")
        hist_BkgSub_LS.Write(f"Mass_BkgSub_{centr_bins}_LS")
        hist_BkgSub_LS_norm.Write(f"Mass_BkgSub_{centr_bins}_LS_norm")
        histRatioCentrPM[index].Write()
        histRatioCentrPP[index].Write()
        histRatioCentrMM[index].Write()
        histRatioPtPM[index].Write()
        histRatioPtPP[index].Write()
        histRatioPtMM[index].Write()
        histRatioRapPM[index].Write()
        histRatioRapPP[index].Write()
        histRatioRapMM[index].Write()
        hist_R_values.Write()

        index = index + 1
        print("index: ", index)

    f_out.cd()
    histRatioSumCentrPM.Write()
    histRatioSumCentrPP.Write()
    histRatioSumCentrMM.Write()
    histRatioSumPtPM.Write()
    histRatioSumPtPP.Write()
    histRatioSumPtMM.Write()
    histRatioSumRapPM.Write()
    histRatioSumRapPP.Write()
    histRatioSumRapMM.Write()
    # Create a canvas e draw the histogram for the Racc
    canvas_Int = ROOT.TCanvas("canvas", "R bin per bin", 800, 600)
    hist_R_values_Int.SetMarkerStyle(20)
    hist_R_values_Int.SetMarkerSize(0.8)
    hist_R_values_Int.SetLineColor(ROOT.kBlue)
    hist_R_values_Int.SetLineWidth(2)
    hist_R_values_Int.Draw("EP")

    # Aggiungi una linea tratteggiata a y = 1
    line = ROOT.TLine(hist_R_values_Int.GetXaxis().GetXmin(), 1, hist_R_values_Int.GetXaxis().GetXmax(), 1)
    line.SetLineColor(ROOT.kRed)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw("same")
    canvas_Int.Update()
    canvas_Int.Write("Racc Int")
    canvasCentr.Write()
    canvaspT.Write()


    # Close files
    f_in.Close()
    f_out.Close()
    print("All histograms processed and saved.")

# Define manual pT bins as a list of tuples
pt_bins = [(0, 2), (2, 4), (4, 5), (4,6), (5, 12), (6,12)]
centr_bins = [(0, 20), (20, 40), (40, 60), (60, 90)]

if __name__ == '__main__':
    rebin_centr_pT(pt_bins, centr_bins)
