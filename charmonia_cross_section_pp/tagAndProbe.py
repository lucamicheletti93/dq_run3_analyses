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
        tagAndProbe(inputCfg)

def tagAndProbe(config):
    """
    function to compute tag and probe
    """
    LoadStyle()

    dfDataAllProbes = pd.read_csv("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/tag_and_probes/LHC24_pass1_min_bias/AllProbes/systematic_sig_Jpsi.txt", sep=' ')
    ptMin = dfDataAllProbes["x_min"]
    ptMax = dfDataAllProbes["x_max"]
    ptEdges = np.append(ptMin.to_numpy(), ptMax.to_numpy()[len(ptMax)-1])

    jpsiDataAllProbes = dfDataAllProbes["val"]
    jpsiStatDataAllProbes = dfDataAllProbes["stat"]
    jpsiSystDataAllProbes = dfDataAllProbes["syst"]

    histJpsiStatDataAllProbes = ROOT.TH1D("histJpsiStatDataAllProbes", ";#it{p}_{T}^{Probe} (GeV/#it{c});Passing probes / All probes", len(ptEdges)-1, ptEdges)
    histJpsiSystDataAllProbes = ROOT.TH1D("histJpsiSystDataAllProbes", ";#it{p}_{T}^{Probe} (GeV/#it{c});Passing probes / All probes", len(ptEdges)-1, ptEdges)

    SetHistStat(histJpsiStatDataAllProbes, 20, ROOT.kRed+1)
    SetHistSyst(histJpsiStatDataAllProbes, 20, ROOT.kRed+1)

    dfDataPassingProbes = pd.read_csv("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/tag_and_probes/LHC24_pass1_min_bias/PassingProbes/systematic_sig_Jpsi.txt", sep=' ')
    jpsiDataPassingProbes = dfDataPassingProbes["val"]
    jpsiStatDataPassingProbes = dfDataPassingProbes["stat"]
    jpsiSystDataPassingProbes = dfDataPassingProbes["syst"]

    histJpsiStatDataPassingProbes = ROOT.TH1D("histJpsiStatDataPassingProbes", ";#it{p}_{T}^{Probe} (GeV/#it{c});Passing probes / All probes", len(ptEdges)-1, ptEdges)
    histJpsiSystDataPassingProbes = ROOT.TH1D("histJpsiSystDataPassingProbes", ";#it{p}_{T}^{Probe} (GeV/#it{c});Passing probes / All probes", len(ptEdges)-1, ptEdges)

    SetHistStat(histJpsiStatDataPassingProbes, 20, ROOT.kRed+1)
    SetHistSyst(histJpsiStatDataPassingProbes, 20, ROOT.kRed+1)

    for iBin in range(0, len(ptMin)):
        histJpsiStatDataAllProbes.SetBinContent(iBin+1, jpsiDataAllProbes[iBin])
        histJpsiStatDataAllProbes.SetBinError(iBin+1, jpsiStatDataAllProbes[iBin])
        histJpsiStatDataAllProbes.SetBinContent(iBin+1, jpsiDataAllProbes[iBin])
        histJpsiStatDataAllProbes.SetBinError(iBin+1, jpsiSystDataAllProbes[iBin])
        histJpsiSystDataAllProbes.SetBinContent(iBin+1, jpsiSystDataAllProbes[iBin]/ jpsiDataAllProbes[iBin])

        histJpsiStatDataPassingProbes.SetBinContent(iBin+1, jpsiDataPassingProbes[iBin])
        histJpsiStatDataPassingProbes.SetBinError(iBin+1, jpsiStatDataPassingProbes[iBin])
        histJpsiStatDataPassingProbes.SetBinContent(iBin+1, jpsiDataPassingProbes[iBin])
        histJpsiStatDataPassingProbes.SetBinError(iBin+1, jpsiSystDataPassingProbes[iBin])
        histJpsiSystDataPassingProbes.SetBinContent(iBin+1, jpsiSystDataPassingProbes[iBin]/ jpsiDataPassingProbes[iBin])

    histPassingAllRatioData = histJpsiStatDataPassingProbes.Clone("histPassingAllRatioData")
    histPassingAllRatioData.Divide(histJpsiStatDataAllProbes)



    dfMcAllProbes = pd.read_csv("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/tag_and_probes/LHC25g8/AllProbes/systematic_sig_Jpsi.txt", sep=' ')
    jpsiMcAllProbes = dfMcAllProbes["val"]
    jpsiStatMcAllProbes = dfMcAllProbes["stat"]
    jpsiSystMcAllProbes = dfMcAllProbes["syst"]

    histJpsiStatMcAllProbes = ROOT.TH1D("histJpsiStatMcAllProbes", ";#it{p}_{T}^{Probe} (GeV/#it{c});Passing probes / All probes", len(ptEdges)-1, ptEdges)
    histJpsiSystMcAllProbes = ROOT.TH1D("histJpsiSystMcAllProbes", ";#it{p}_{T}^{Probe} (GeV/#it{c});Passing probes / All probes", len(ptEdges)-1, ptEdges)

    SetHistStat(histJpsiStatMcAllProbes, 20, ROOT.kAzure+4)
    SetHistSyst(histJpsiStatMcAllProbes, 20, ROOT.kAzure+4)

    dfMcPassingProbes = pd.read_csv("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/tag_and_probes/LHC25g8/PassingProbes/systematic_sig_Jpsi.txt", sep=' ')
    jpsiMcPassingProbes = dfMcPassingProbes["val"]
    jpsiStatMcPassingProbes = dfMcPassingProbes["stat"]
    jpsiSystMcPassingProbes = dfMcPassingProbes["syst"]

    histJpsiStatMcPassingProbes = ROOT.TH1D("histJpsiStatMcPassingProbes", ";#it{p}_{T}^{Probe} (GeV/#it{c});Passing probes / All probes", len(ptEdges)-1, ptEdges)
    histJpsiSystMcPassingProbes = ROOT.TH1D("histJpsiSystMcPassingProbes", ";#it{p}_{T}^{Probe} (GeV/#it{c});Passing probes / All probes", len(ptEdges)-1, ptEdges)

    SetHistStat(histJpsiStatMcPassingProbes, 20, ROOT.kAzure+4)
    SetHistSyst(histJpsiStatMcPassingProbes, 20, ROOT.kAzure+4)

    for iBin in range(0, len(ptMin)):
        histJpsiStatMcAllProbes.SetBinContent(iBin+1, jpsiMcAllProbes[iBin])
        histJpsiStatMcAllProbes.SetBinError(iBin+1, jpsiStatMcAllProbes[iBin])
        histJpsiStatMcAllProbes.SetBinContent(iBin+1, jpsiMcAllProbes[iBin])
        histJpsiStatMcAllProbes.SetBinError(iBin+1, jpsiSystMcAllProbes[iBin])
        histJpsiSystMcAllProbes.SetBinContent(iBin+1, jpsiSystMcAllProbes[iBin]/ jpsiMcAllProbes[iBin])

        histJpsiStatMcPassingProbes.SetBinContent(iBin+1, jpsiMcPassingProbes[iBin])
        histJpsiStatMcPassingProbes.SetBinError(iBin+1, jpsiStatMcPassingProbes[iBin])
        histJpsiStatMcPassingProbes.SetBinContent(iBin+1, jpsiMcPassingProbes[iBin])
        histJpsiStatMcPassingProbes.SetBinError(iBin+1, jpsiSystMcPassingProbes[iBin])
        histJpsiSystMcPassingProbes.SetBinContent(iBin+1, jpsiSystMcPassingProbes[iBin]/ jpsiMcPassingProbes[iBin])

    histPassingAllRatioMc = histJpsiStatMcPassingProbes.Clone("histPassingAllRatioMc")
    histPassingAllRatioMc.Divide(histJpsiStatMcAllProbes)

    histRatioDataMc = histPassingAllRatioData.Clone("histRatioDataMc")
    histRatioDataMc.Divide(histPassingAllRatioMc)
    histRatioDataMc.GetXaxis().SetTitle("#it{p}_{T}^{Probe} (GeV/#it{c})")
    histRatioDataMc.GetYaxis().SetTitle("Data / MC")
    SetHistStat(histRatioDataMc, 20, ROOT.kBlack)

    #funcRatioData = ROOT.TF1("funcRatioData", FuncEff, 0, 20, 2)
    #funcRatioData.SetParameters(0.638376, -1.96188)
    #histPassingAllRatioData.Fit(funcRatioData, "R0")
    #funcRatioData.SetLineColor(ROOT.kRed+1)
    #funcRatioData.SetLineStyle(ROOT.kDashed)

    #funcRatioMc = ROOT.TF1("funcRatioMc", FuncEff, 0, 20, 2)
    #funcRatioMc.SetParameters(0.638376, -1.96188)
    #histPassingAllRatioMc.Fit(funcRatioMc, "R0")
    #funcRatioMc.SetLineColor(ROOT.kAzure+4)
    #funcRatioMc.SetLineStyle(ROOT.kDashed)

    funcPol0 = ROOT.TF1("funcPol0", "[0]", 0, 20)
    funcPol0.SetParameter(0, 1)
    histRatioDataMc.Fit(funcPol0, "R")

    padAspectRatio = 0.6
    canvasTagAndProbes = ROOT.TCanvas("canvasTagAndProbes", "", 800, 1000)
    ROOT.gStyle.SetOptFit(True)
    ROOT.gStyle.SetOptStat(False)

    pad1 = ROOT.TPad("pad1", "pad1", 0.05, 1-padAspectRatio, 0.99, 0.99)
    pad1.Draw()
    pad1.cd()
    histPassingAllRatioData.GetYaxis().SetRangeUser(0, 1.2)
    histPassingAllRatioData.Draw("EP")
    histPassingAllRatioMc.Draw("EP SAME")
    #funcRatioData.Draw("SAME")
    #funcRatioMc.Draw("SAME")
    lineUnity = ROOT.TLine(0., 1., 20., 1.)
    lineUnity.SetLineStyle(ROOT.kDashed)
    lineUnity.SetLineColor(ROOT.kGray+1)
    lineUnity.SetLineWidth(2)
    lineUnity.Draw()

    legendTagAndProbes = ROOT.TLegend(0.50, 0.40, 0.80, 0.60, " ", "brNDC")
    SetLegend(legendTagAndProbes)
    legendTagAndProbes.AddEntry(histPassingAllRatioData, "LHC24 pass1 minBias", "EP")
    legendTagAndProbes.AddEntry(histPassingAllRatioMc, "LHC25g8", "EP")
    legendTagAndProbes.Draw("SAME")
    canvasTagAndProbes.Update()

    canvasTagAndProbes.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0.05, 0.05, 0.99, 1-padAspectRatio)
    pad2.SetBottomMargin(0.25)
    pad2.Draw()
    pad2.cd()
    # Update to the size of the small pad
    histRatioDataMc.GetXaxis().SetTitleSize(0.05 * (1./padAspectRatio))
    histRatioDataMc.GetYaxis().SetTitleSize(0.045 * (1./padAspectRatio))
    histRatioDataMc.GetXaxis().SetLabelSize(0.045 * (1./padAspectRatio))
    histRatioDataMc.GetYaxis().SetLabelSize(0.045 * (1./padAspectRatio))
    histRatioDataMc.GetYaxis().SetRangeUser(0.8, 1.2)
    histRatioDataMc.Draw("H")
    lineUnity.Draw()
    funcPol0.Draw("SAME")
    canvasTagAndProbes.Update()

    input()

    canvasTagAndProbes.SaveAs("figures/matching_eff/tag_and_probes.pdf")

def FuncEff(x, par):
    return 1 / (1 + ROOT.TMath.Exp(-par[0]*(x[0] - par[1])))


if __name__ == '__main__':
    main()