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
    parser.add_argument("--do_xsec", help="Compute cross section", action="store_true")
    args = parser.parse_args()
    

    with open(args.cfgFileName, 'r') as yml_cfg:
        inputCfg = yaml.load(yml_cfg, yaml.FullLoader)

    if args.do_xsec:
        xsec(inputCfg)

def xsec(config):
    """
    function to compute the cross section
    """
    LoadStyle()
    BrJpsiToMuMu =  0.05961
    errBrJpsiToMuMu = 0.033e-2
    BrPsi2sToMuMu =  0.008
    errBrPsi2sToMuMu = 0.6e-3

    deltaRap = config["inputs"]["deltaRap"]

    print("***** Extract raw yield *****")
    dfJpsiRawYieldInt = pd.read_csv(config["inputs"]["fInRawYieldInt"], sep=' ')
    ptMinInt = dfJpsiRawYieldInt["x_min"]
    ptMaxInt = dfJpsiRawYieldInt["x_max"]
    ptEdgesInt = np.append(ptMinInt.to_numpy(), ptMaxInt.to_numpy()[0])
    jpsiRawYieldInt = dfJpsiRawYieldInt["val"] / deltaRap
    jpsiStatRawYieldInt = dfJpsiRawYieldInt["stat"] / deltaRap
    jpsiSystRawYieldInt = dfJpsiRawYieldInt["syst"] / deltaRap

    histStatRawYieldInt = ROOT.TH1D("histStatRawYieldInt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdgesInt)-1, ptEdgesInt)
    histSystRawYieldInt = ROOT.TH1D("histSystRawYieldInt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdgesInt)-1, ptEdgesInt)
    histSystRelRawYieldInt = ROOT.TH1D("histSystRelRawYieldInt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdgesInt)-1, ptEdgesInt)

    SetHistStat(histStatRawYieldInt, 20, ROOT.kRed+1)
    SetHistSyst(histSystRawYieldInt, 20, ROOT.kRed+1)

    for iBin in range(0, len(ptMinInt)):
        histStatRawYieldInt.SetBinContent(iBin+1, jpsiRawYieldInt[iBin])
        histStatRawYieldInt.SetBinError(iBin+1, jpsiStatRawYieldInt[iBin])
        histSystRawYieldInt.SetBinContent(iBin+1, jpsiRawYieldInt[iBin])
        histSystRawYieldInt.SetBinError(iBin+1, jpsiSystRawYieldInt[iBin])
        histSystRelRawYieldInt.SetBinContent(iBin+1, jpsiSystRawYieldInt[iBin]/ jpsiRawYieldInt[iBin])

    canvasRawYieldInt = ROOT.TCanvas("canvasRawYieldInt", "", 800, 600)
    histStatRawYieldInt.Draw("EP")
    histSystRawYieldInt.Draw("E2P SAME")
    canvasRawYieldInt.Update()


    dfJpsiRawYieldVsPt = pd.read_csv(config["inputs"]["fInRawYieldVsPt"], sep=' ')
    ptMin = dfJpsiRawYieldVsPt["x_min"]
    ptMax = dfJpsiRawYieldVsPt["x_max"]
    ptCenters = (ptMax + ptMin) / 2.
    ptWidths = (ptMax - ptMin) / 2.
    ptSystWidths = sqrtsSystWidths = np.repeat(0.2, len(ptWidths))
    ptEdges = np.append(ptMin.to_numpy(), ptMax.to_numpy()[len(ptMax)-1])
    jpsiRawYieldVsPt = dfJpsiRawYieldVsPt["val"] / deltaRap
    jpsiStatRawYieldVsPt = dfJpsiRawYieldVsPt["stat"] / deltaRap
    jpsiSystRawYieldVsPt = dfJpsiRawYieldVsPt["syst"] / deltaRap

    histStatRawYieldVsPt = ROOT.TH1D("histStatRawYieldVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdges)-1, ptEdges)
    histSystRawYieldVsPt = ROOT.TH1D("histSystRawYieldVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdges)-1, ptEdges)
    histSystRelRawYieldVsPt = ROOT.TH1D("histSystRelRawYieldVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdges)-1, ptEdges)

    SetHistStat(histStatRawYieldVsPt, 20, ROOT.kRed+1)
    SetHistSyst(histSystRawYieldVsPt, 20, ROOT.kRed+1)

    for iBin in range(0, len(ptMin)):
        histStatRawYieldVsPt.SetBinContent(iBin+1, jpsiRawYieldVsPt[iBin])
        histStatRawYieldVsPt.SetBinError(iBin+1, jpsiStatRawYieldVsPt[iBin])
        histSystRawYieldVsPt.SetBinContent(iBin+1, jpsiRawYieldVsPt[iBin])
        histSystRawYieldVsPt.SetBinError(iBin+1, jpsiSystRawYieldVsPt[iBin])
        histSystRelRawYieldVsPt.SetBinContent(iBin+1, jpsiSystRawYieldVsPt[iBin] / jpsiRawYieldVsPt[iBin])

    canvasRawYieldVsPt = ROOT.TCanvas("canvasRawYieldVsPt", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histStatRawYieldVsPt.Draw("EP")
    histSystRawYieldVsPt.Draw("E2P SAME")
    canvasRawYieldVsPt.Update()

    print("***** Extract Axe *****")
    dfJpsiAxeInt = pd.read_csv(config["inputs"]["fInAxeInt"], sep=' ')
    jpsiAxeInt = dfJpsiAxeInt["val"]
    jpsiStatAxeInt = dfJpsiAxeInt["stat"]

    histStatAxeInt = ROOT.TH1D("histStatAxeInt", ";#it{p}_{T} (GeV/#it{c});A#times#epsilon", len(ptEdgesInt)-1, ptEdgesInt)

    SetHistStat(histStatAxeInt, 20, ROOT.kRed+1)

    for iBin in range(0, len(ptMinInt)):
        histStatAxeInt.SetBinContent(iBin+1, jpsiAxeInt[iBin])
        histStatAxeInt.SetBinError(iBin+1, jpsiStatAxeInt[iBin])

    canvasAxeInt = ROOT.TCanvas("canvasAxeInt", "", 800, 600)
    histStatAxeInt.Draw("EP")
    canvasAxeInt.Update()

    dfJpsiAxeVsPt = pd.read_csv(config["inputs"]["fInAxeVsPt"], sep=' ')
    jpsiAxeVsPt = dfJpsiAxeVsPt["val"]
    jpsiStatAxeVsPt = dfJpsiAxeVsPt["stat"]

    histStatAxeVsPt = ROOT.TH1D("histStatAxeVsPt", ";#it{p}_{T} (GeV/#it{c});A#times#epsilon", len(ptEdges)-1, ptEdges)

    SetHistStat(histStatAxeVsPt, 20, ROOT.kRed+1)

    for iBin in range(0, len(ptMin)):
        histStatAxeVsPt.SetBinContent(iBin+1, jpsiAxeVsPt[iBin])
        histStatAxeVsPt.SetBinError(iBin+1, jpsiStatAxeVsPt[iBin])

    canvasAxeVsPt = ROOT.TCanvas("canvasAxeVsPt", "", 800, 600)
    histStatAxeVsPt.Draw("EP")
    canvasAxeVsPt.Update()

    print("***** Compute cross section *****")
    # Lumi computed with normalization.C  + check_normalization
    #lumi = 0.170285 # pb-1
    fInLumi = ROOT.TFile(config["inputs"]["fInLumi"], "READ")
    histLumi = fInLumi.Get("histLumi")
    lumi = histLumi.GetBinContent(1) # pb-1
    print("luminosity = ", lumi, " pb-1")

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    # List of systematics computed externally
    systRelLumi = 0.1
    histSystRelLumiInt = ROOT.TH1D("histSystRelLumiInt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdgesInt)-1, ptEdgesInt)
    histSystRelLumiVsPt = ROOT.TH1D("histSystRelLumiVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdges)-1, ptEdges)

    for iBin in range(0, len(ptMin)):
        if iBin < 1:
            histSystRelLumiInt.SetBinContent(iBin+1, systRelLumi)
        histSystRelLumiVsPt.SetBinContent(iBin+1, systRelLumi)
    print("Syst. Rel. Luminosity = ", systRelLumi, " pb-1")

    systRelBrJpsiToMuMu = errBrJpsiToMuMu / BrJpsiToMuMu
    histSystRelBrJpsiToMuMuInt = ROOT.TH1D("histSystRelBrJpsiToMuMuInt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdgesInt)-1, ptEdgesInt)
    histSystRelBrJpsiToMuMuVsPt = ROOT.TH1D("histSystRelBrJpsiToMuMuVsPt", ";#it{p}_{T} (GeV/#it{c});dN/d#it{y}", len(ptEdges)-1, ptEdges)

    for iBin in range(0, len(ptMin)):
        if iBin < 1:
            histSystRelBrJpsiToMuMuInt.SetBinContent(iBin+1, systRelBrJpsiToMuMu)
        histSystRelBrJpsiToMuMuVsPt.SetBinContent(iBin+1, systRelBrJpsiToMuMu)
    print("Syst. Rel. BR Jpsi->mumu = ", systRelBrJpsiToMuMu)
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

    # > > > Pt - integrated < < < #
    sqrtsCenters = np.array([13.600])
    sqrtsWidths = np.array([0.])
    sqrtsSystWidths = np.repeat(0.2, len(sqrtsCenters))

    jpsiXsecInt = (jpsiRawYieldInt) / (jpsiAxeInt * BrJpsiToMuMu * lumi * 1e6)
    jpsiStatXsecInt = (jpsiStatRawYieldInt) / (jpsiAxeInt * BrJpsiToMuMu * lumi * 1e6)
    jpsiSystXsecInt = (jpsiSystRawYieldInt) / (jpsiAxeInt * BrJpsiToMuMu * lumi * 1e6)
    # Add all systematics contributions
    jpsiSystXsecInt = jpsiXsecInt * np.sqrt((jpsiSystXsecInt / jpsiXsecInt)**2 + systRelBrJpsiToMuMu**2 + systRelLumi**2)

    graStatXsecInt = ROOT.TGraphErrors(len(sqrtsCenters), sqrtsCenters, np.array(jpsiXsecInt), sqrtsWidths, np.array(jpsiStatXsecInt))
    graSystXsecInt = ROOT.TGraphErrors(len(sqrtsCenters), sqrtsCenters, np.array(jpsiXsecInt), sqrtsSystWidths, np.array(jpsiSystXsecInt))

    print(f'[INFO] dsigma/dy (J/psi) = {jpsiXsecInt[0]:.3f} +/- {jpsiStatXsecInt[0]:.3f} +/- {jpsiSystXsecInt[0]:.3f} mub')

    SetGraStat(graStatXsecInt, 20, ROOT.kRed+1)
    SetGraSyst(graSystXsecInt, 20, ROOT.kRed+1)

    histStatXsecInt = ROOT.TH1D("histStatXsecInt", ";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{y} (#mub)", len(ptEdgesInt)-1, ptEdgesInt)
    histSystXsecInt = ROOT.TH1D("histSystXsecInt", ";#it{p}_{T} (GeV/#it{c});d#sigma/d#it{y} (#mub)", len(ptEdgesInt)-1, ptEdgesInt)

    SetHistStat(histStatXsecInt, 20, ROOT.kRed+1)
    SetHistSyst(histSystXsecInt, 20, ROOT.kRed+1)

    for iBin in range(0, len(ptMinInt)):
        histStatXsecInt.SetBinContent(iBin+1, jpsiXsecInt[iBin])
        histStatXsecInt.SetBinError(iBin+1, jpsiStatXsecInt[iBin])
        histSystXsecInt.SetBinContent(iBin+1, jpsiXsecInt[iBin])
        histSystXsecInt.SetBinError(iBin+1, jpsiSystXsecInt[iBin])

    canvasXsecInt = ROOT.TCanvas("canvasXsecInt", "", 800, 600)
    histStatXsecInt.Draw("EP")
    histSystXsecInt.Draw("E2P SAME")
    canvasXsecInt.Update()


    # Compare to Run 2 results
    filePathJpsiRun2 = "/Users/lucamicheletti/GITHUB/dq_run3_analyses/psi2S_jpsi_ratio/results/HEPData-ins1511865-v1-yaml/Table13.yaml"
    sqrtsCentersRun2, sqrtsWidthsRun2, jpsiXsecVsSqrtsRun2, jpsiStatXsecVsSqrtsRun2, jpsiSystXsecVsSqrtsRun2 = ExtractFromYaml(filePathJpsiRun2)
    sqrtsSystWidthsRun2 = np.repeat(200, len(sqrtsCentersRun2))

    graStatXsecIntRun2 = ROOT.TGraphErrors(len(sqrtsCentersRun2), np.array(sqrtsCentersRun2)/1e3, np.array(jpsiXsecVsSqrtsRun2)/1e3, np.array(sqrtsWidthsRun2)/1e3, np.array(jpsiStatXsecVsSqrtsRun2)/1e3)
    graSystXsecIntRun2 = ROOT.TGraphErrors(len(sqrtsCentersRun2), np.array(sqrtsCentersRun2)/1e3, np.array(jpsiXsecVsSqrtsRun2)/1e3, np.array(sqrtsSystWidthsRun2)/1e3, np.array(jpsiSystXsecVsSqrtsRun2)/1e3)

    SetGraStat(graStatXsecIntRun2, 20, ROOT.kBlack)
    SetGraSyst(graSystXsecIntRun2, 20, ROOT.kBlack)

    canvasXsecIntVsRun2 = ROOT.TCanvas("canvasXsecIntVsRun2", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    histGridXsecIntVsSqrt = ROOT.TH2D("histGridXsecIntVsSqrt", ";#sqrt{#it{s}} (TeV);d#sigma/d#it{y} (#mub)", 100, 0, 15, 100, 0, 10)
    histGridXsecIntVsSqrt.Draw()
    graStatXsecIntRun2.Draw("EP SAME")
    graSystXsecIntRun2.Draw("E2P SAME")
    graStatXsecInt.Draw("EP SAME")
    graSystXsecInt.Draw("E2P SAME")  
    canvasXsecIntVsRun2.Update()

    legendXsecIntVsRun2 = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendXsecIntVsRun2)
    legendXsecIntVsRun2.AddEntry(graSystXsecIntRun2, "EPJC 77 (2017) 392 ", "FP")
    legendXsecIntVsRun2.AddEntry(graStatXsecInt, "#sqrt{#it{s}} = 13.6 TeV", "FP")
    legendXsecIntVsRun2.Draw("SAME")
    canvasXsecIntVsRun2.Update()

    # > > > Pt dependence < < < #
    jpsiXsecVsPt = (jpsiRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * lumi * 1e6 * (2 * ptWidths))
    jpsiStatXsecVsPt = (jpsiStatRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * lumi * 1e6 * (2 * ptWidths))
    jpsiSystXsecVsPt = (jpsiSystRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * lumi * 1e6 * (2 * ptWidths))
    # Add all systematics contributions
    jpsiSystXsecVsPt = jpsiXsecVsPt * np.sqrt((jpsiSystXsecVsPt / jpsiXsecVsPt)**2 + systRelBrJpsiToMuMu**2 + systRelLumi**2)

    graStatXsecVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiXsecVsPt), np.array(ptWidths), np.array(jpsiStatXsecVsPt))
    graSystXsecVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiXsecVsPt), np.array(ptSystWidths), np.array(jpsiSystXsecVsPt))

    SetGraStat(graStatXsecVsPt, 20, ROOT.kRed+1)
    SetGraSyst(graSystXsecVsPt, 20, ROOT.kRed+1)

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
    filePathJpsi = "/Users/lucamicheletti/GITHUB/dq_run3_analyses/psi2S_jpsi_ratio/results/HEPData-77781-v1-yaml/Table1.yaml"
    ptCentersRun2, ptWidthsRun2, jpsiXsecVsPtRun2, jpsiStatXsecVsPtRun2, jpsiSystXsecVsPtRun2 = ExtractFromYaml(filePathJpsi)
    jpsiPtSystWidthsRun2 = np.repeat(0.2, len(ptCentersRun2))

    graStatXsecVsPtRun2 = ROOT.TGraphErrors(len(ptCentersRun2), np.array(ptCentersRun2), np.array(jpsiXsecVsPtRun2)/1e3, np.array(ptWidthsRun2), np.array(jpsiStatXsecVsPtRun2)/1e3)
    graSystXsecVsPtRun2 = ROOT.TGraphErrors(len(ptCentersRun2), np.array(ptCentersRun2), np.array(jpsiXsecVsPtRun2)/1e3, jpsiPtSystWidthsRun2, np.array(jpsiSystXsecVsPtRun2)/1e3)

    SetGraStat(graStatXsecVsPtRun2, 20, ROOT.kBlack)
    SetGraSyst(graSystXsecVsPtRun2, 20, ROOT.kBlack)

    canvasXsecVsPtVsRun2 = ROOT.TCanvas("canvasXsecVsPtVsRun2", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    ROOT.gPad.SetLogy(True)
    histGridXsecVsPt = ROOT.TH2D("histGridXsecVsPt", ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{y} d#it{p}_{T} (#mub / GeV/#it{c})", 100, 0, 20, 100, 1e-4, 5)
    histGridXsecVsPt.Draw()
    graStatXsecVsPtRun2.Draw("EP SAME")
    graSystXsecVsPtRun2.Draw("E2P SAME")
    graStatXsecVsPt.Draw("EP SAME")
    graSystXsecVsPt.Draw("E2P SAME")  

    legendXsecVsPtVsRun2 = ROOT.TLegend(0.20, 0.20, 0.40, 0.45, " ", "brNDC")
    SetLegend(legendXsecVsPtVsRun2)
    legendXsecVsPtVsRun2.AddEntry(graSystXsecVsPtRun2, "#sqrt{#it{s}} = 13 TeV", "FP")
    legendXsecVsPtVsRun2.AddEntry(graStatXsecVsPt, "#sqrt{#it{s}} = 13.6 TeV", "FP")
    legendXsecVsPtVsRun2.Draw("SAME")
    canvasXsecVsPtVsRun2.Update()

    input()

    print(f"***** Write results in {config["outputs"]["fOut"]} *****")
    fOut = ROOT.TFile(config["outputs"]["fOut"], "RECREATE")
    histStatRawYieldInt.Write()
    histSystRawYieldInt.Write()
    histStatRawYieldVsPt.Write()
    histSystRawYieldVsPt.Write()
    histStatAxeInt.Write()
    histStatAxeVsPt.Write()
    histStatXsecInt.Write()
    histSystXsecInt.Write()
    histStatXsecVsPt.Write()
    histSystXsecVsPt.Write()

    fOut.mkdir("systematics")
    fOut.cd("systematics")
    histSystRelRawYieldInt.Write("syst_raw_yield_int")
    histSystRelRawYieldVsPt.Write("syst_raw_yield_vs_pt")
    histSystRelBrJpsiToMuMuInt.Write("syst_br_int")
    histSystRelBrJpsiToMuMuVsPt.Write("syst_br_vs_pt")
    histSystRelLumiInt.Write("syst_lumi_int")
    histSystRelLumiVsPt.Write("syst_lumi_vs_pt")
    fOut.Close()


if __name__ == '__main__':
    main()