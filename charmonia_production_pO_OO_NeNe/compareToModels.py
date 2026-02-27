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
        compareToModels(inputCfg)

def compareToModels(config):
    """
    function to compare results with theory models
    """
    LoadStyle()
    ROOT.gStyle.SetOptStat(False)

    dfJpsiRaaThuVsCentr = pd.read_csv(config["inputs"]["fInJpsiRaaThuVsCentr"], sep=' ')
    centrMin = dfJpsiRaaThuVsCentr["x_min"]
    centrMax = dfJpsiRaaThuVsCentr["x_max"]
    centrCenters = (centrMax + centrMin) / 2.
    centrWidths = (centrMax - centrMin) / 2.
    centrSystWidths = np.repeat(1.5, len(centrWidths))
    centrEdges = np.append(centrMin.to_numpy(), centrMax.to_numpy()[len(centrMax)-1])
    jpsiRaaThuVsCentr = dfJpsiRaaThuVsCentr["val"]

    npartCenters = [25.21, 21.03, 16.88, 13.05, 9.88, 7.39, 5.51, 4.15, 3.17, 2.42]
    npartWidths = [0.69, 0.84, 0.89, 0.75, 0.61, 0.48, 0.36, 0.26, 0.23, 0.16]

    graJpsiRaaThuVsCentr = ROOT.TGraphErrors(len(npartCenters), np.array(npartCenters), np.array(jpsiRaaThuVsCentr), np.array(npartWidths), np.zeros(len(npartCenters)))
    SetHistStat(graJpsiRaaThuVsCentr, 20, ROOT.kAzure+4)

    histGridJpsiRaaVsCentr = ROOT.TH2D("histGridJpsiRaaVsCentr", ";#it{N}_{part};#it{R}_{OO}", 100, 0, 30, 100, 0, 1.2)

    canvasJpsiRaaThuVsCentr = ROOT.TCanvas("canvasJpsiRaaThuVsCentr", "", 800, 600)
    histGridJpsiRaaVsCentr.Draw()
    graJpsiRaaThuVsCentr.Draw("EP SAME")

    lineUnityVsCentr = ROOT.TLine(0, 1, 30, 1)
    lineUnityVsCentr.SetLineColor(ROOT.kGray+1)
    lineUnityVsCentr.SetLineWidth(2)
    lineUnityVsCentr.SetLineStyle(ROOT.kDashed)
    lineUnityVsCentr.Draw()

    legendRaaVsCentr = ROOT.TLegend(0.20, 0.20, 0.40, 0.40, " ", "brNDC")
    SetLegend(legendRaaVsCentr)
    legendRaaVsCentr.SetTextSize(0.045)
    legendRaaVsCentr.AddEntry(graJpsiRaaThuVsCentr, "Transport, THU", "P")
    legendRaaVsCentr.Draw("SAME")

    canvasJpsiRaaThuVsCentr.Update()




    dfJpsiRaaThuVsPt = pd.read_csv(config["inputs"]["fInJpsiRaaThuVsPt"], sep=' ')
    ptMin = dfJpsiRaaThuVsPt["x_min"]
    ptMax = dfJpsiRaaThuVsPt["x_max"]
    ptCenters = (ptMax + ptMin) / 2.
    ptWidths = (ptMax - ptMin) / 2.
    ptSystWidths = np.repeat(1.5, len(ptWidths))
    ptEdges = np.append(ptMin.to_numpy(), ptMax.to_numpy()[len(ptMax)-1])
    jpsiRaaThuVsPt = dfJpsiRaaThuVsPt["val"]

    graJpsiRaaThuVsPt = ROOT.TGraphErrors(len(ptCenters), np.array(ptCenters), np.array(jpsiRaaThuVsPt), np.array(ptWidths), np.zeros(len(ptCenters)))
    SetHistStat(graJpsiRaaThuVsPt, 20, ROOT.kAzure+4)

    histGridJpsiRaaVsPt = ROOT.TH2D("histGridJpsiRaaVsPt", ";#it{p}_{T} (GeV/#it{c});#it{R}_{OO}", 100, 0, 16, 100, 0, 1.2)

    canvasJpsiRaaThuVsPt = ROOT.TCanvas("canvasJpsiRaaThuVsPt", "", 800, 600)
    histGridJpsiRaaVsPt.Draw()
    graJpsiRaaThuVsPt.Draw("EP SAME")

    lineUnityVsPt = ROOT.TLine(0, 1, 16, 1)
    lineUnityVsPt.SetLineColor(ROOT.kGray+1)
    lineUnityVsPt.SetLineWidth(2)
    lineUnityVsPt.SetLineStyle(ROOT.kDashed)
    lineUnityVsPt.Draw()

    legendRaaVsPt = ROOT.TLegend(0.60, 0.20, 0.80, 0.40, " ", "brNDC")
    SetLegend(legendRaaVsPt)
    legendRaaVsPt.SetTextSize(0.045)
    legendRaaVsPt.AddEntry(graJpsiRaaThuVsPt, "Transport, THU", "P")
    legendRaaVsPt.Draw("SAME")


    canvasJpsiRaaThuVsPt.Update()





    

    input()

    """
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
    """


if __name__ == '__main__':
    main()