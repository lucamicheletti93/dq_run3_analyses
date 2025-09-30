import sys
import math
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
        interpolation(inputCfg)

def interpolation(config):
    """
    function to perform interpolation
    """
    LoadStyle()
    ROOT.gStyle.SetOptStat(False)

    sqrtsInterp = config["inputs"]["sqrtsInterp"]
    rangeSqrtsInterp = config["inputs"]["rangeSqrtsInterp"]

    ptEdgesInterp = config["inputs"]["ptEdgesInterp"]
    ptCentersInterp = config["inputs"]["ptCentersInterp"]
    ptWidthsInterp = config["inputs"]["ptWidthsInterp"]
    ptSystWidthsInterp = config["inputs"]["ptSystWidthsInterp"]

    sqrts = np.array(list(config["inputs"]["HepDataList"].keys()))
    HepDataList = config["inputs"]["HepDataList"]

    # pp@2.76TeV
    fInPathJpsiXsec2TeV = HepDataList[2.76]
    ptCenters2TeV, ptWidths2TeV, jpsiXsec2TeVVsPt, jpsiStatXsec2TeVVsPt, jpsiSystXsec2TeVVsPt = ExtractFromYaml(fInPathJpsiXsec2TeV)

    ptSystWidths2TeV = np.repeat(0.25, len(ptCenters2TeV))
    
    ptEdges2TeV = np.array(ptCenters2TeV) - np.array(ptWidths2TeV)
    ptEdges2TeV = np.append(ptEdges2TeV, 8.0)
    jpsiXsec2TeVVsPtRebin, jpsiStatXsec2TeVVsPtRebin, jpsiSystXsec2TeVVsPtRebin = Rebin(ptEdges2TeV, jpsiXsec2TeVVsPt, jpsiStatXsec2TeVVsPt, jpsiSystXsec2TeVVsPt, ptEdgesInterp)

    graStatJpsiXsec2TeV = ROOT.TGraphErrors(len(ptCenters2TeV), np.array(ptCenters2TeV), np.array(jpsiXsec2TeVVsPt), np.array(ptWidths2TeV), np.array(jpsiStatXsec2TeVVsPt))
    graSystJpsiXsec2TeV = ROOT.TGraphErrors(len(ptCenters2TeV), np.array(ptCenters2TeV), np.array(jpsiXsec2TeVVsPt), np.array(ptSystWidths2TeV), np.array(jpsiSystXsec2TeVVsPt))

    SetGraStat(graStatJpsiXsec2TeV, 20, ROOT.kGreen+1)
    SetGraSyst(graSystJpsiXsec2TeV, 20, ROOT.kGreen+1)

    graStatJpsiXsec2TeVRebin = ROOT.TGraphErrors(len(ptCentersInterp), np.array(ptCentersInterp), np.array(jpsiXsec2TeVVsPtRebin), np.array(ptWidthsInterp), np.array(jpsiStatXsec2TeVVsPtRebin))
    graSystJpsiXsec2TeVRebin = ROOT.TGraphErrors(len(ptCentersInterp), np.array(ptCentersInterp), np.array(jpsiXsec2TeVVsPtRebin), np.array(ptSystWidthsInterp), np.array(jpsiSystXsec2TeVVsPtRebin))

    SetGraStat(graStatJpsiXsec2TeVRebin, 20, ROOT.kGreen+1)
    SetGraSyst(graSystJpsiXsec2TeVRebin, 20, ROOT.kGreen+1)

    # pp@5TeV
    fInPathJpsiXsec5TeV = HepDataList[5.02]
    ptCenters5TeV, ptWidths5TeV, jpsiXsec5TeVVsPt, jpsiStatXsec5TeVVsPt, jpsiSystXsec5TeVVsPt = ExtractFromYaml(fInPathJpsiXsec5TeV)

    jpsiXsec5TeVVsPt = np.array(jpsiXsec5TeVVsPt)/1e3
    jpsiStatXsec5TeVVsPt = np.array(jpsiStatXsec5TeVVsPt)/1e3
    jpsiSystXsec5TeVVsPt = np.array(jpsiSystXsec5TeVVsPt)/1e3

    ptSystWidths5TeV = np.repeat(0.25, len(ptCenters5TeV))

    ptEdges5TeV = np.array(ptCenters5TeV) - np.array(ptWidths5TeV)
    jpsiXsec5TeVVsPtRebin, jpsiStatXsec5TeVVsPtRebin, jpsiSystXsec5TeVVsPtRebin = Rebin(ptEdges5TeV, jpsiXsec5TeVVsPt, jpsiStatXsec5TeVVsPt, jpsiSystXsec5TeVVsPt, ptEdgesInterp)

    graStatJpsiXsec5TeV = ROOT.TGraphErrors(len(ptCenters5TeV), np.array(ptCenters5TeV), np.array(jpsiXsec5TeVVsPt), np.array(ptWidths5TeV), np.array(jpsiStatXsec5TeVVsPt))
    graSystJpsiXsec5TeV = ROOT.TGraphErrors(len(ptCenters5TeV), np.array(ptCenters5TeV), np.array(jpsiXsec5TeVVsPt), np.array(ptSystWidths5TeV), np.array(jpsiSystXsec5TeVVsPt))

    SetGraStat(graStatJpsiXsec5TeV, 20, ROOT.kBlack)
    SetGraSyst(graSystJpsiXsec5TeV, 20, ROOT.kBlack)

    graStatJpsiXsec5TeVRebin = ROOT.TGraphErrors(len(ptCentersInterp), np.array(ptCentersInterp), np.array(jpsiXsec5TeVVsPtRebin), np.array(ptWidthsInterp), np.array(jpsiStatXsec5TeVVsPtRebin))
    graSystJpsiXsec5TeVRebin = ROOT.TGraphErrors(len(ptCentersInterp), np.array(ptCentersInterp), np.array(jpsiXsec5TeVVsPtRebin), np.array(ptSystWidthsInterp), np.array(jpsiSystXsec5TeVVsPtRebin))

    SetGraStat(graStatJpsiXsec5TeVRebin, 20, ROOT.kBlack)
    SetGraSyst(graSystJpsiXsec5TeVRebin, 20, ROOT.kBlack)


    # pp@7TeV
    fInPathJpsiXsec7TeV = HepDataList[7.00]
    ptCenters7TeV, ptWidths7TeV, jpsiXsec7TeVVsPt, jpsiStatXsec7TeVVsPt, jpsiSystXsec7TeVVsPt = ExtractFromYaml(fInPathJpsiXsec7TeV)
    ptSystWidths7TeV = np.repeat(0.25, len(ptCenters7TeV))
    
    ptEdges7TeV = np.array(ptCenters7TeV) - np.array(ptWidths7TeV)
    jpsiXsec7TeVVsPtRebin, jpsiStatXsec7TeVVsPtRebin, jpsiSystXsec7TeVVsPtRebin = Rebin(ptEdges7TeV, jpsiXsec7TeVVsPt, jpsiStatXsec7TeVVsPt, jpsiSystXsec7TeVVsPt, ptEdgesInterp)

    graStatJpsiXsec7TeV = ROOT.TGraphErrors(len(ptCenters7TeV), np.array(ptCenters7TeV), np.array(jpsiXsec7TeVVsPt), np.array(ptWidths7TeV), np.array(jpsiStatXsec7TeVVsPt))
    graSystJpsiXsec7TeV = ROOT.TGraphErrors(len(ptCenters7TeV), np.array(ptCenters7TeV), np.array(jpsiXsec7TeVVsPt), np.array(ptSystWidths7TeV), np.array(jpsiSystXsec7TeVVsPt))

    SetGraStat(graStatJpsiXsec7TeV, 20, ROOT.kAzure+4)
    SetGraSyst(graSystJpsiXsec7TeV, 20, ROOT.kAzure+4)

    graStatJpsiXsec7TeVRebin = ROOT.TGraphErrors(len(ptCentersInterp), np.array(ptCentersInterp), np.array(jpsiXsec7TeVVsPtRebin), np.array(ptWidthsInterp), np.array(jpsiStatXsec7TeVVsPtRebin))
    graSystJpsiXsec7TeVRebin = ROOT.TGraphErrors(len(ptCentersInterp), np.array(ptCentersInterp), np.array(jpsiXsec7TeVVsPtRebin), np.array(ptSystWidthsInterp), np.array(jpsiSystXsec7TeVVsPtRebin))

    SetGraStat(graStatJpsiXsec7TeVRebin, 20, ROOT.kAzure+4)
    SetGraSyst(graSystJpsiXsec7TeVRebin, 20, ROOT.kAzure+4)

    # pp@8TeV
    fInPathJpsiXsec8TeV = HepDataList[8.00]
    ptCenters8TeV, ptWidths8TeV, jpsiXsec8TeVVsPt, jpsiStatXsec8TeVVsPt, jpsiSystXsec8TeVVsPt = ExtractFromYaml(fInPathJpsiXsec8TeV)
    ptSystWidths8TeV = np.repeat(0.25, len(ptCenters8TeV))

    ptEdges8TeV = np.array(ptCenters8TeV) - np.array(ptWidths8TeV)
    jpsiXsec8TeVVsPtRebin, jpsiStatXsec8TeVVsPtRebin, jpsiSystXsec8TeVVsPtRebin = Rebin(ptEdges8TeV, jpsiXsec8TeVVsPt, jpsiStatXsec8TeVVsPt, jpsiSystXsec8TeVVsPt, ptEdgesInterp)

    graStatJpsiXsec8TeV = ROOT.TGraphErrors(len(ptCenters8TeV), np.array(ptCenters8TeV), np.array(jpsiXsec8TeVVsPt), np.array(ptWidths8TeV), np.array(jpsiStatXsec8TeVVsPt))
    graSystJpsiXsec8TeV = ROOT.TGraphErrors(len(ptCenters8TeV), np.array(ptCenters8TeV), np.array(jpsiXsec8TeVVsPt), np.array(ptSystWidths8TeV), np.array(jpsiSystXsec8TeVVsPt))

    SetGraStat(graStatJpsiXsec8TeV, 20, ROOT.kOrange+7)
    SetGraSyst(graSystJpsiXsec8TeV, 20, ROOT.kOrange+7)

    graStatJpsiXsec8TeVRebin = ROOT.TGraphErrors(len(ptCentersInterp), np.array(ptCentersInterp), np.array(jpsiXsec8TeVVsPtRebin), np.array(ptWidthsInterp), np.array(jpsiStatXsec8TeVVsPtRebin))
    graSystJpsiXsec8TeVRebin = ROOT.TGraphErrors(len(ptCentersInterp), np.array(ptCentersInterp), np.array(jpsiXsec8TeVVsPtRebin), np.array(ptSystWidthsInterp), np.array(jpsiSystXsec8TeVVsPtRebin))

    SetGraStat(graStatJpsiXsec8TeVRebin, 20, ROOT.kOrange+7)
    SetGraSyst(graSystJpsiXsec8TeVRebin, 20, ROOT.kOrange+7)

    # pp@13TeV
    fInPathJpsiXsec13TeV = HepDataList[13.00]
    ptCenters13TeV, ptWidths13TeV, jpsiXsec13TeVVsPt, jpsiStatXsec13TeVVsPt, jpsiSystXsec13TeVVsPt = ExtractFromYaml(fInPathJpsiXsec13TeV)
    ptSystWidths13TeV = np.repeat(0.25, len(ptCenters13TeV))

    jpsiXsec13TeVVsPt = np.array(jpsiXsec13TeVVsPt)/1e3
    jpsiStatXsec13TeVVsPt = np.array(jpsiStatXsec13TeVVsPt)/1e3
    jpsiSystXsec13TeVVsPt = np.array(jpsiSystXsec13TeVVsPt)/1e3

    ptEdges13TeV = np.array(ptCenters13TeV) - np.array(ptWidths13TeV)
    jpsiXsec13TeVVsPtRebin, jpsiStatXsec13TeVVsPtRebin, jpsiSystXsec13TeVVsPtRebin = Rebin(ptEdges13TeV, jpsiXsec13TeVVsPt, jpsiStatXsec13TeVVsPt, jpsiSystXsec13TeVVsPt, ptEdgesInterp)

    graStatJpsiXsec13TeV = ROOT.TGraphErrors(len(ptCenters13TeV), np.array(ptCenters13TeV), np.array(jpsiXsec13TeVVsPt), np.array(ptWidths13TeV), np.array(jpsiStatXsec13TeVVsPt))
    graSystJpsiXsec13TeV = ROOT.TGraphErrors(len(ptCenters13TeV), np.array(ptCenters13TeV), np.array(jpsiXsec13TeVVsPt), np.array(ptSystWidths13TeV), np.array(jpsiSystXsec13TeVVsPt))

    SetGraStat(graStatJpsiXsec13TeV, 20, ROOT.kMagenta)
    SetGraSyst(graSystJpsiXsec13TeV, 20, ROOT.kMagenta)

    graStatJpsiXsec13TeVRebin = ROOT.TGraphErrors(len(ptCentersInterp), np.array(ptCentersInterp), np.array(jpsiXsec13TeVVsPtRebin), np.array(ptWidthsInterp), np.array(jpsiStatXsec13TeVVsPtRebin))
    graSystJpsiXsec13TeVRebin = ROOT.TGraphErrors(len(ptCentersInterp), np.array(ptCentersInterp), np.array(jpsiXsec13TeVVsPtRebin), np.array(ptSystWidthsInterp), np.array(jpsiSystXsec13TeVVsPtRebin))

    SetGraStat(graStatJpsiXsec13TeVRebin, 20, ROOT.kMagenta)
    SetGraSyst(graSystJpsiXsec13TeVRebin, 20, ROOT.kMagenta)

    canvasJpsiXsecVsPt = ROOT.TCanvas("canvasJpsiXsecVsPt", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    ROOT.gPad.SetLogy(True)
    histGridJpsiXsecVsPt = ROOT.TH2D("histGridJpsiXsecVsPt", ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T} d#it{y} (#mub/GeV/#it{c})", 100, 0, 10, 100, 0.01, 12)
    histGridJpsiXsecVsPt.Draw()
    #graStatJpsiXsec2TeV.Draw("EP SAME")
    #graSystJpsiXsec2TeV.Draw("E2P SAME")
    graStatJpsiXsec2TeVRebin.Draw("EP SAME")
    graSystJpsiXsec2TeVRebin.Draw("E2P SAME")
    #graStatJpsiXsec5TeV.Draw("EP SAME")
    #graSystJpsiXsec5TeV.Draw("E2P SAME")
    graStatJpsiXsec5TeVRebin.Draw("EP SAME")
    graSystJpsiXsec5TeVRebin.Draw("E2P SAME")
    #graStatJpsiXsec7TeV.Draw("EP SAME")
    #graSystJpsiXsec7TeV.Draw("E2P SAME")
    graStatJpsiXsec7TeVRebin.Draw("EP SAME")
    graSystJpsiXsec7TeVRebin.Draw("E2P SAME")
    #graStatJpsiXsec8TeV.Draw("EP SAME")
    #graSystJpsiXsec8TeV.Draw("E2P SAME")
    graStatJpsiXsec8TeVRebin.Draw("EP SAME")
    graSystJpsiXsec8TeVRebin.Draw("E2P SAME")
    #graStatJpsiXsec13TeV.Draw("EP SAME")
    #graSystJpsiXsec13TeV.Draw("E2P SAME")
    graStatJpsiXsec13TeVRebin.Draw("EP SAME")
    graSystJpsiXsec13TeVRebin.Draw("E2P SAME")
    #latexTitle.DrawLatex(0.2, 0.87, "ALICE Preliminary")
    #latexTitle.DrawLatex(0.2, 0.80, "J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4, #it{p}_{T} < 20 GeV/#it{c}")
    #latexTitle.DrawLatex(0.2, 0.73, f'#it{{L_{{int}}}} = {lumi:4.3f} pb^{{-1}}')

    legendJpsiXsecVsPt = ROOT.TLegend(0.65, 0.60, 0.90, 0.90, " ", "brNDC")
    SetLegend(legendJpsiXsecVsPt)
    legendJpsiXsecVsPt.SetTextSize(0.045)
    legendJpsiXsecVsPt.AddEntry(graSystJpsiXsec2TeVRebin, "#sqrt{#it{s}} = 2.76 TeV", "FP")
    legendJpsiXsecVsPt.AddEntry(graSystJpsiXsec5TeVRebin, "#sqrt{#it{s}} = 5.02 TeV", "FP")
    legendJpsiXsecVsPt.AddEntry(graSystJpsiXsec7TeVRebin, "#sqrt{#it{s}} = 7.00 TeV", "FP")
    legendJpsiXsecVsPt.AddEntry(graSystJpsiXsec8TeVRebin, "#sqrt{#it{s}} = 8.00 TeV", "FP")
    legendJpsiXsecVsPt.AddEntry(graSystJpsiXsec13TeVRebin, "#sqrt{#it{s}} = 13.00 TeV", "FP")
    legendJpsiXsecVsPt.Draw("SAME")

    canvasJpsiXsecVsPt.Update()

    matrixJpsiXsec = np.column_stack((jpsiXsec2TeVVsPtRebin, jpsiXsec5TeVVsPtRebin, jpsiXsec7TeVVsPtRebin, jpsiXsec8TeVVsPtRebin, jpsiXsec13TeVVsPtRebin))
    matrixStatJpsiXsec = np.column_stack((jpsiStatXsec2TeVVsPtRebin, jpsiStatXsec5TeVVsPtRebin, jpsiStatXsec7TeVVsPtRebin, jpsiStatXsec8TeVVsPtRebin, jpsiStatXsec13TeVVsPtRebin))
    matrixSystJpsiXsec = np.column_stack((jpsiSystXsec2TeVVsPtRebin, jpsiSystXsec5TeVVsPtRebin, jpsiSystXsec7TeVVsPtRebin, jpsiSystXsec8TeVVsPtRebin, jpsiSystXsec13TeVVsPtRebin))

    print(matrixJpsiXsec)
    print(matrixStatJpsiXsec)
    print(matrixSystJpsiXsec)


    graErrXsecJpsis = []

    for iPt in range(0, len(ptEdgesInterp)-1):
        xSecs = matrixJpsiXsec[iPt, :]
        statXsecs = matrixStatJpsiXsec[iPt, :]
        systXsecs = matrixSystJpsiXsec[iPt, :]

        errXsecs = np.sqrt(np.array(statXsecs)**2 + np.array(systXsecs)**2)

        graErrXsecJpsis.append(ROOT.TGraphErrors(len(sqrts), np.array(sqrts), np.array(xSecs), np.zeros(len(sqrts)), np.array(errXsecs)))

        SetGraStat(graErrXsecJpsis[iPt], 20, ROOT.kBlack)
        graErrXsecJpsis[iPt].SetMarkerSize(0.5)


    canvasInterpolation = ROOT.TCanvas("canvasInterpolation", "", 2400, 1200)
    canvasInterpolation.Divide(4, 2)

    histGridInterp = ROOT.TH2D("histGridInterp", ";#sqrt{#it{s}} (TeV);d^{2}#sigma/d#it{p}_{T} d#it{y} (#mub/GeV/#it{c})", 100, 0, 15, 100, 0.01, 3)
    
    latexTitle = ROOT.TLatex()
    latexTitle.SetNDC()
    latexTitle.SetTextSize(0.05)
    latexTitle.SetTextFont(42)

    latexText = ROOT.TLatex()
    latexText.SetNDC()
    latexText.SetTextSize(0.035)
    latexText.SetTextFont(42)

    meanValInterp, meanErrInterp, meanSystInterp  = [], [], []

    for iPt in range(0, len(ptEdgesInterp)-1):
        canvasInterpolation.cd(iPt+1)
        ROOT.gPad.SetLogy(True)

        xSecs = matrixJpsiXsec[iPt, :]
        histGridInterp.GetYaxis().SetRangeUser(0.5 * xSecs.min(), 1.5 * xSecs.max())
        histGridInterp.DrawCopy()

        funcPol1 = ROOT.TF1("funcPol1", "pol1", 2, 15, 2)
        funcPol1.SetLineColor(ROOT.kRed+1)
        fitResPol1 = graErrXsecJpsis[iPt].Fit(funcPol1, "RS")
        covMatrixPol1 = fitResPol1.GetCovarianceMatrix()
        valInterpPol1 = funcPol1.Eval(sqrtsInterp)
        errInterpPol1 = funcPol1.IntegralError(sqrtsInterp-rangeSqrtsInterp,sqrtsInterp+rangeSqrtsInterp, fitResPol1.GetParams() , covMatrixPol1.GetMatrixArray())

        funcPowerLow = ROOT.TF1("funcPowerLow", "[0]* x**([1])", 2, 15)
        funcPowerLow.SetLineColor(ROOT.kBlue+1)
        fitResPowerLow = graErrXsecJpsis[iPt].Fit(funcPowerLow, "RS")
        covMatrixPowerLow = fitResPowerLow.GetCovarianceMatrix()
        valInterpPowerLow = funcPowerLow.Eval(sqrtsInterp)
        errInterpPowerLow = funcPowerLow.IntegralError(sqrtsInterp-rangeSqrtsInterp,sqrtsInterp+rangeSqrtsInterp, fitResPowerLow.GetParams() , covMatrixPowerLow.GetMatrixArray())

        funcExpo = ROOT.TF1("funcExpo", "expo", 2, 15, 2)
        funcExpo.SetLineColor(ROOT.kGreen+1)
        fitResExpo = graErrXsecJpsis[iPt].Fit(funcExpo, "RS")
        covMatrixExpo = fitResExpo.GetCovarianceMatrix()
        valInterpExpo = funcExpo.Eval(sqrtsInterp)
        errInterpExpo = funcExpo.IntegralError(sqrtsInterp-rangeSqrtsInterp,sqrtsInterp+rangeSqrtsInterp, fitResExpo.GetParams() , covMatrixExpo.GetMatrixArray()) / (2 * rangeSqrtsInterp)

        latexText.DrawLatex(0.50, 0.25, f'Pol1    = {valInterpPol1:0.3f} #pm {errInterpPol1:0.3f}')
        latexText.DrawLatex(0.50, 0.30, f'Pow.Low = {valInterpPowerLow:0.3f} #pm {errInterpPowerLow:0.3f}')
        latexText.DrawLatex(0.50, 0.35, f'Expo    = {valInterpExpo:0.3f} #pm {errInterpExpo:0.3f}')

        graErrXsecJpsis[iPt].Draw("EP SAME")
        funcPol1.Draw("SAME")
        funcPowerLow.Draw("SAME")
        funcExpo.Draw("SAME")

        latexTitle.DrawLatex(0.25, 0.85, f'{ptEdgesInterp[iPt]} < #it{{p}}_{{T}} < {ptEdgesInterp[iPt+1]} GeV/#it{{c}}')

        vecInterpVals = [valInterpPol1, valInterpPowerLow, valInterpExpo]
        meanValInterp.append((valInterpPol1 + valInterpPowerLow + valInterpExpo) / 3.)
        meanErrInterp.append((errInterpPol1 + errInterpPowerLow + errInterpExpo) / 3.)
        meanSystInterp.append(ComputeStdDev(vecInterpVals))

        latexText.DrawLatex(0.45, 0.40, f'#color[2]]{{Mean    = {meanValInterp[iPt]:0.3f} #pm {meanErrInterp[iPt]:0.3f} #pm {meanSystInterp[iPt]:0.3f}}}')

    canvasInterpolation.cd(8)

    graStatJpsiXsecInterp = ROOT.TGraphErrors(len(ptCentersInterp), np.array(ptCentersInterp), np.array(meanValInterp), np.array(ptWidthsInterp), np.zeros(len(meanErrInterp)))
    graSystJpsiXsecInterp = ROOT.TGraphErrors(len(ptCentersInterp), np.array(ptCentersInterp), np.array(meanValInterp), np.array(ptSystWidthsInterp), np.array(meanSystInterp))

    histStatJpsiXsecInterp = ROOT.TH1D("histStatJpsiXsecInterp", ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T} d#it{y} (#mub/GeV/#it{c})", len(ptEdgesInterp)-1, np.array(ptEdgesInterp))
    histSystJpsiXsecInterp = ROOT.TH1D("histSystJpsiXsecInterp", ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T} d#it{y} (#mub/GeV/#it{c})", len(ptEdgesInterp)-1, np.array(ptEdgesInterp))

    for iPt in range(0, len(ptEdgesInterp)-1):
        histStatJpsiXsecInterp.SetBinContent(iPt+1, meanValInterp[iPt])
        histStatJpsiXsecInterp.SetBinError(iPt+1, meanErrInterp[iPt])

        histSystJpsiXsecInterp.SetBinContent(iPt+1, meanValInterp[iPt])
        histSystJpsiXsecInterp.SetBinError(iPt+1, meanSystInterp[iPt])

    SetGraStat(graStatJpsiXsecInterp, 20, ROOT.kRed+1)
    SetGraSyst(graSystJpsiXsecInterp, 20, ROOT.kRed+1)

    ROOT.gPad.SetLogy(True)
    histGridJpsiXsecInterp = ROOT.TH2D("histGridJpsiXsecInterp", ";#it{p}_{T} (GeV/#it{c});d^{2}#sigma/d#it{p}_{T} d#it{y} (#mub/GeV/#it{c})", 100, 0, 8, 100, 0.01, 5)
    histGridJpsiXsecInterp.Draw()
    graStatJpsiXsec5TeVRebin.Draw("EP SAME")
    graSystJpsiXsec5TeVRebin.Draw("E2P SAME")
    graStatJpsiXsecInterp.Draw("EP SAME")
    graSystJpsiXsecInterp.Draw("E2P SAME")

    legendJpsiXsecInterp = ROOT.TLegend(0.55, 0.70, 0.80, 0.90, " ", "brNDC")
    SetLegend(legendJpsiXsecInterp)
    legendJpsiXsecInterp.SetTextSize(0.045)
    legendJpsiXsecInterp.AddEntry(graSystJpsiXsec5TeVRebin, "#sqrt{#it{s}} = 5.02 TeV", "FP")
    legendJpsiXsecInterp.AddEntry(graSystJpsiXsecInterp, "#sqrt{#it{s}} = 5.36 TeV", "FP")
    legendJpsiXsecInterp.Draw("SAME")

    canvasInterpolation.Update()

    input()

    canvasInterpolation.SaveAs("figures/interpolation/interpolation_vs_pt.pdf")

    fOut = ROOT.TFile(config["outputs"]["fOut"], "RECREATE")
    histStatJpsiXsecInterp.Write()
    histSystJpsiXsecInterp.Write()
    fOut.Close()


    exit()

def Rebin(oldEdges, vals, stats, systs, ptEdgesInterp):
    oldEdges = np.array(oldEdges)
    vals = np.array(vals)
    stats = np.array(stats)
    systs = np.array(systs)
    ptEdgesInterp = np.array(ptEdgesInterp)

    newVals, newStats, newSysts = [], [], []
    
    for i in range(len(ptEdgesInterp)-1):
        newLedge, newHedge = ptEdgesInterp[i], ptEdgesInterp[i+1]
        val = 0
        stat = 0
        syst = 0
        binWidth = 0
        for j in range(len(oldEdges)-1):
            oldLedge, oldHedge = oldEdges[j], oldEdges[j+1]
            if newLedge == oldLedge and newHedge == oldHedge:
                newVals.append(vals[j])
                newStats.append(stats[j])
                newSysts.append(systs[j])
                break
            if (oldLedge >= newLedge) and (oldHedge <= newHedge):
                val += vals[j] * (oldHedge - oldLedge)
                stat += (stats[j]**2) * (oldHedge - oldLedge)
                syst += (systs[j]**2) * (oldHedge - oldLedge)
                binWidth += (oldHedge - oldLedge)
                if oldHedge == newHedge:
                    newVals.append(val / binWidth)
                    newStats.append(math.sqrt(stat) / binWidth)
                    newSysts.append(math.sqrt(syst) / binWidth)
                    val = 0
                    stat = 0
                    syst = 0

    return newVals, newStats, newSysts

def ComputeStdDev(parValArray):
    '''
    Method to evaluate the dispersion of the data around the mean
    '''
    mean = 0
    for parVal in parValArray:
        mean += parVal
    mean = mean / len(parValArray)
    stdDev = 0
    for parVal in parValArray:
        stdDev += (parVal - mean) * (parVal - mean)
    stdDev = math.sqrt(stdDev / len(parValArray))
    return stdDev

if __name__ == '__main__':
    main()