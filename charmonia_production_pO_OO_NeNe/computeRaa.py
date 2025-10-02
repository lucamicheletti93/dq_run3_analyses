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

    print("***** Extract pp reference [mub] *****")

    fInPPrefVsPt = ROOT.TFile(config["inputs"]["fInPPredVsPt"], "READ")
    histStatPPrefXsec = fInPPrefVsPt.Get("histStatJpsiXsecInterp")
    histSystPPrefXsec = fInPPrefVsPt.Get("histSystJpsiXsecInterp")
    
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
    #histSystRawYieldVsPt.Draw("E2P SAME")
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

    print("***** Compute corrected yield *****")
    # Lumi computed with normalization.C  + check_normalization
    #lumi = 0.170285 # pb-1
    fInLumi = ROOT.TFile(config["inputs"]["fInLumi"], "READ")
    histLumi = fInLumi.Get("histLumi")
    lumi = histLumi.GetBinContent(1)
    histNevMinBias = fInLumi.Get("histNevMinBias")
    nevMinBias = histNevMinBias.GetBinContent(1)
    print("N. min bias events = ", nevMinBias)

    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    # List of systematics computed externally
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
    # * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
    # > > > Pt dependence < < < #
    jpsiXsecVsPt = (jpsiRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * lumi * 1e6 * (2 * ptWidths))
    jpsiStatXsecVsPt = (jpsiStatRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * lumi * 1e6 * (2 * ptWidths))
    jpsiSystXsecVsPt = (jpsiSystRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * lumi * 1e6 * (2 * ptWidths))

    jpsiRaaVsPt = (jpsiRawYieldVsPt) / (jpsiAxeVsPt * BrJpsiToMuMu * nevMinBias * (2 * ptWidths) * Taa * np.array(jpsiPPrefXsecVsPt))
    jpsiStatRaaVsPt = jpsiRaaVsPt * (np.array(jpsiStatXsecVsPt) / np.array(jpsiXsecVsPt))
    jpsiSystRaaVsPt = jpsiRaaVsPt * (np.array(jpsiSystXsecVsPt) / np.array(jpsiXsecVsPt))

    # Add all systematics contributions
    jpsiSystXsecVsPt = jpsiXsecVsPt * np.sqrt((jpsiSystXsecVsPt / jpsiXsecVsPt)**2 + systRelBrJpsiToMuMu**2 + systRelTrackingEff**2 + systRelMatchingEff**2)

    graStatXsecVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiXsecVsPt), np.array(ptWidths), np.array(jpsiStatXsecVsPt))
    graSystXsecVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiXsecVsPt), np.array(ptSystWidths), np.array(jpsiSystXsecVsPt))

    graStatRaaVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt), np.array(ptWidths), np.array(jpsiStatRaaVsPt))
    graSystRaaVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaVsPt), np.array(ptSystWidths), np.array(jpsiSystRaaVsPt))

    SetGraStat(graStatXsecVsPt, 20, ROOT.kRed+1)
    SetGraSyst(graSystXsecVsPt, 20, ROOT.kRed+1)

    SetGraStat(graStatRaaVsPt, 20, ROOT.kRed+1)
    SetGraSyst(graSystRaaVsPt, 20, ROOT.kRed+1)


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

    canvasXsecVsPt = ROOT.TCanvas("canvasXsecVsPt", "", 800, 600)
    ROOT.gPad.SetLogy(True)
    histStatXsecVsPt.Draw("EP")
    histSystXsecVsPt.Draw("E2P SAME")
    canvasXsecVsPt.Update()

    # Compare to Run 2 results
    print("***** Extract RAA from Run 2 *****")
    filePathJpsi = "HEP_Data/jpsi_raa_PbPb_5TeV.yaml"
    ptCentersRun2, ptWidthsRun2, jpsiRaaVsPtRun2, jpsiStatRaaVsPtRun2, jpsiSystRaaVsPtRun2 = ExtractFromYaml(filePathJpsi)
    jpsiPtSystWidthsRun2 = np.repeat(0.2, len(ptCentersRun2))

    graStatRaaVsPtRun2 = ROOT.TGraphErrors(len(ptCentersRun2), np.array(ptCentersRun2), np.array(jpsiRaaVsPtRun2), np.array(ptWidthsRun2), np.array(jpsiStatRaaVsPtRun2))
    graSystRaaVsPtRun2 = ROOT.TGraphErrors(len(ptCentersRun2), np.array(ptCentersRun2), np.array(jpsiRaaVsPtRun2), jpsiPtSystWidthsRun2, np.array(jpsiSystRaaVsPtRun2))

    SetGraStat(graStatRaaVsPtRun2, 20, ROOT.kGray+2)
    SetGraSyst(graSystRaaVsPtRun2, 20, ROOT.kGray+2)

    lineUnity = ROOT.TLine(0., 1., 12., 1.)
    lineUnity.SetLineColor(ROOT.kGray+1)
    lineUnity.SetLineWidth(2)
    lineUnity.SetLineStyle(ROOT.kDashed)

    canvasRaaVsPtVsRun2 = ROOT.TCanvas("canvasRaaVsPtVsRun2", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    histGridRaaVsPt = ROOT.TH2D("histGridRaaVsPt", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}", 100, 0, 12, 100, 0, 1.7)
    histGridRaaVsPt.Draw()
    lineUnity.Draw("SAME")
    graStatRaaVsPtRun2.Draw("EP SAME")
    graSystRaaVsPtRun2.Draw("E2P SAME")
    graStatRaaVsPt.Draw("EP SAME")
    graSystRaaVsPt.Draw("E2P SAME")

    legendRaaVsPtVsRun2 = ROOT.TLegend(0.20, 0.70, 0.40, 0.90, " ", "brNDC")
    SetLegend(legendRaaVsPtVsRun2)
    legendRaaVsPtVsRun2.SetTextSize(0.045)
    legendRaaVsPtVsRun2.AddEntry(graStatRaaVsPt, "#sqrt{#it{s}} = 5.36 TeV, OO, 0#minus100%", "FP")
    legendRaaVsPtVsRun2.AddEntry(graStatRaaVsPtRun2, "#sqrt{#it{s}} = 5.02 TeV, Pb-Pb, 0#minus90%", "FP")
    legendRaaVsPtVsRun2.Draw("SAME")

    latexTitle.DrawLatex(0.22, 0.87, "ALICE Preliminary, J/#psi #rightarrow #mu^{#plus}#mu^{#minus}, 2.5 < #it{y} < 4")

    canvasRaaVsPtVsRun2.Update()

    input()

    print(f"***** Write results in {config["outputs"]["fOut"]} *****")
    fOut = ROOT.TFile(config["outputs"]["fOut"], "RECREATE")
    histStatRawYieldVsPt.Write()
    histSystRawYieldVsPt.Write()
    histStatAxeVsPt.Write()
    histStatXsecVsPt.Write()
    histSystXsecVsPt.Write()

    fOut.mkdir("systematics")
    fOut.cd("systematics")
    histSystRelRawYieldVsPt.Write("syst_raw_yield_vs_pt")
    histSystRelLumiVsPt.Write("syst_lumi_vs_pt")
    histSystRelBrJpsiToMuMuVsPt.Write("syst_br_vs_pt")
    histSystRelTrackingEffVsPt.Write("syst_tracking_eff_vs_pt")
    histSystRelMatchingEffVsPt.Write("syst_matching_eff_vs_pt")
    histSystRelMcRealisticVsPt.Write("syst_mc_realisticness_vs_pt")
    fOut.Close()

    canvasRawYieldVsPt.SaveAs("figures/raw_yield/raw_yeild_jpsi.pdf")
    canvasAxeVsPt.SaveAs("figures/axe/axe_jpsi.pdf")
    canvasRaaVsPtVsRun2.SaveAs("figures/raa/pt_differential_raa_jpsi.pdf")


if __name__ == '__main__':
    main()