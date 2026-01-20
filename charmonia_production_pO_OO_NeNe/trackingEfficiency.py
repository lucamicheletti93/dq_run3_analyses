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
    parser.add_argument("--hist", help="Create histos for trk eff calculation", action="store_true")
    parser.add_argument("--run", help="Compute tracking efficiency", action="store_true")
    args = parser.parse_args()
    

    with open(args.cfgFileName, 'r') as yml_cfg:
        inputCfg = yaml.load(yml_cfg, yaml.FullLoader)

    if args.hist:
        createHistograms(inputCfg)

    if args.run:
        trackingEfficiency(inputCfg)

def createHistograms(config):
    """
    function to create the histograms for the trk efficiency
    """
    LoadStyle()

    hits = ["N12", "N10", "N02", "N34", "N30", "N04", "N56", "N50", "N06", "N78", "N70", "N08", "N910", "N90", "N010"] 

    fIn = ROOT.TFile(config["inputs"]["fIn"], "READ")
    histSparse = fIn.Get("task-muon-mch-trk-efficiency/hHitsEtaPtPhi")

    histProjEta = []
    histProjPt = []
    histProjPhi = []
    for iHit, hit in enumerate(hits):
        histTmp = histSparse.Clone(f'histTmp_{hit}')
        histTmp.GetAxis(0).SetRange(histSparse.GetAxis(0).FindBin(hit), histSparse.GetAxis(0).FindBin(hit))
        histProjEta.append(histTmp.Projection(1, f'histProjEta_{hit}'))
        histProjPt.append(histTmp.Projection(2, f'histProjPt_{hit}'))
        histProjPhi.append(histTmp.Projection(3, f'histProjPhi_{hit}'))

        histProjEta[iHit].SetName(f'histProjEta_{hit}')
        histProjPt[iHit].SetName(f'histProjPt_{hit}')
        histProjPhi[iHit].SetName(f'histProjPhi_{hit}')

        histProjEta[iHit].SetTitle(f'#eta distribution {hit}')
        histProjPt[iHit].SetTitle(f'pT distribution {hit}')
        histProjPhi[iHit].SetTitle(f'#phi distribution {hit}')

    # Eta distribution
    canvasEta = ROOT.TCanvas("canvasEta", "", 1000, 1000)
    canvasEta.Divide(4, 4)
    for i in range(0, 15):
        canvasEta.cd(i+1)
        histProjEta[i].Draw("H")
    canvasEta.Update()

    # Pt distribution
    canvasPt = ROOT.TCanvas("canvasPt", "", 1000, 1000)
    canvasPt.Divide(4, 4)
    for i in range(0, 15):
        canvasPt.cd(i+1)
        histProjPt[i].Draw("H")
    canvasPt.Update()

    # Phi distribution
    canvasPhi = ROOT.TCanvas("canvasPhi", "", 1000, 1000)
    canvasPhi.Divide(4, 4)
    for i in range(0, 15):
        canvasPhi.cd(i+1)
        histProjPhi[i].Draw("H")
    canvasPhi.Update()

    input()

    print(f'[INFO] Writing results in {config["outputs"]["fOut"]} ...')
    fOut = ROOT.TFile(config["outputs"]["fOut"], "RECREATE")
    for iHit, hit in enumerate(hits):
        histProjEta[iHit].Write()
        histProjPt[iHit].Write()
        histProjPhi[iHit].Write()
    fOut.Close()
    fIn.Close()

    exit()

def trackingEfficiency(config):
    """
    function to compute the trk efficiency
    """
    LoadStyle()

    hits = ["N12", "N10", "N02", "N34", "N30", "N04", "N56", "N50", "N06", "N78", "N70", "N08", "N910", "N90", "N010"]
    #hits = ["N12", "N10", "N02", "N34", "N30", "N04", "N56", "N50", "N06"]
    vars = ["Eta", "Pt", "Phi"]
    varNames = ["#eta", "#it{p}_{T} (GeV/#it{c})", "#phi"]

    # Group the stations 
    #stations = list(zip(hits[0::3], hits[1::3], hits[2::3], hits[3::6]))
    group1 = [hits[i:i+3] for i in range(0, 9, 3)]
    group2 = hits[9:]
    stations = group1 + [group2]

    fInMc = ROOT.TFile("/Users/lucamicheletti/GITHUB/dq_run3_analyses/charmonia_production_pO_OO_NeNe/mch_trk_eff/Histograms_trk_eff_LHC25i4.root", "READ")
    fInData = ROOT.TFile("/Users/lucamicheletti/GITHUB/dq_run3_analyses/charmonia_production_pO_OO_NeNe/mch_trk_eff/Histograms_trk_eff_LHC25ae_pass2.root", "READ")

    # eff st1 => 1 - [1 - (N12/(N02+N12))]*[1 - (N12/(N10+N12))]
    histMcEff1 = {}  #[len(vars)][len(hits)]
    histMcEff2 = {}  #[len(vars)][len(hits)]
    histMcEffSt = {}  #[len(vars)][len(hits)]

    for iVar, var in enumerate(vars):
        histMcEff1[var] = {}
        histMcEff2[var] = {}
        histMcEffSt[var] = {}
        for iSt, st in enumerate(stations):
            print("Station", st)

            if iSt < 3:
                print(f'histProj{var}_{st[0]}')
                print(f'histProj{var}_{st[1]}')
                print(f'histProj{var}_{st[2]}')

                histNXY = fInMc.Get(f'histProj{var}_{st[0]}')
                histNX0 = fInMc.Get(f'histProj{var}_{st[1]}')
                histN0Y = fInMc.Get(f'histProj{var}_{st[2]}')
                histMcEff1[var][iSt], histMcEff2[var][iSt], histMcEffSt[var][iSt] = funcSt123Eff(histNXY, histNX0, histN0Y)
            else:
                print(f'histProj{var}_{st[0]}')
                print(f'histProj{var}_{st[1]}')
                print(f'histProj{var}_{st[2]}')
                print(f'histProj{var}_{st[3]}')
                print(f'histProj{var}_{st[4]}')
                print(f'histProj{var}_{st[5]}')

                histN78 = fInMc.Get(f'histProj{var}_{st[0]}')
                histN70 = fInMc.Get(f'histProj{var}_{st[1]}')
                histN08 = fInMc.Get(f'histProj{var}_{st[2]}')
                histN910 = fInMc.Get(f'histProj{var}_{st[3]}')
                histN90 = fInMc.Get(f'histProj{var}_{st[4]}')
                histN010 = fInMc.Get(f'histProj{var}_{st[5]}')
                histMcEffSt[var][iSt] = funcSt45Eff(histN78, histN70, histN08, histN910, histN90, histN010)

    histMcMchTrkEff = []
    for iVar, var in enumerate(vars):
        nBins = histMcEffSt[var][iSt].GetNbinsX()
        binEdges = [histMcEffSt[var][iSt].GetXaxis().GetBinLowEdge(iBin+1) for iBin in range(histMcEffSt[var][iSt].GetNbinsX())]
        binEdges.append(histMcEffSt[var][iSt].GetXaxis().GetBinUpEdge(histMcEffSt[var][iSt].GetNbinsX()))
        histMcMchTrkEff.append(ROOT.TH1D(f'histMcMchTrkEff{var}', f';{varNames[iVar]};Efficiency', nBins, np.array(binEdges)))
        for iBin in range(0, nBins):
            effSt1  = histMcEffSt[var][0].GetBinContent(iBin+1)
            effSt2  = histMcEffSt[var][1].GetBinContent(iBin+1)
            effSt3  = histMcEffSt[var][2].GetBinContent(iBin+1)
            effSt45 = histMcEffSt[var][3].GetBinContent(iBin+1)
            effMchTrkEff = effSt1 * effSt2 * effSt3 * effSt45
            histMcMchTrkEff[iVar].SetBinContent(iBin+1, effMchTrkEff)
            histMcMchTrkEff[iVar].SetBinError(iBin+1, 0)


    histDataEff1 = {}  #[len(vars)][len(hits)]
    histDataEff2 = {}  #[len(vars)][len(hits)]
    histDataEffSt = {}  #[len(vars)][len(hits)]

    for iVar, var in enumerate(vars):
        histDataEff1[var] = {}
        histDataEff2[var] = {}
        histDataEffSt[var] = {}
        for iSt, st in enumerate(stations):
            print("Station", st)

            if iSt < 3:
                print(f'histProj{var}_{st[0]}')
                print(f'histProj{var}_{st[1]}')
                print(f'histProj{var}_{st[2]}')

                histNXY = fInData.Get(f'histProj{var}_{st[0]}')
                histNX0 = fInData.Get(f'histProj{var}_{st[1]}')
                histN0Y = fInData.Get(f'histProj{var}_{st[2]}')
                histDataEff1[var][iSt], histDataEff2[var][iSt], histDataEffSt[var][iSt] = funcSt123Eff(histNXY, histNX0, histN0Y)
            else:
                print(f'histProj{var}_{st[0]}')
                print(f'histProj{var}_{st[1]}')
                print(f'histProj{var}_{st[2]}')
                print(f'histProj{var}_{st[3]}')
                print(f'histProj{var}_{st[4]}')
                print(f'histProj{var}_{st[5]}')

                histN78 = fInData.Get(f'histProj{var}_{st[0]}')
                histN70 = fInData.Get(f'histProj{var}_{st[1]}')
                histN08 = fInData.Get(f'histProj{var}_{st[2]}')
                histN910 = fInData.Get(f'histProj{var}_{st[3]}')
                histN90 = fInData.Get(f'histProj{var}_{st[4]}')
                histN010 = fInData.Get(f'histProj{var}_{st[5]}')
                histDataEffSt[var][iSt] = funcSt45Eff(histN78, histN70, histN08, histN910, histN90, histN010)

    histDataMchTrkEff = []
    for iVar, var in enumerate(vars):
        nBins = histDataEffSt[var][iSt].GetNbinsX()
        binEdges = [histDataEffSt[var][iSt].GetXaxis().GetBinLowEdge(iBin+1) for iBin in range(histDataEffSt[var][iSt].GetNbinsX())]
        binEdges.append(histDataEffSt[var][iSt].GetXaxis().GetBinUpEdge(histDataEffSt[var][iSt].GetNbinsX()))
        histDataMchTrkEff.append(ROOT.TH1D(f'histDataMchTrkEff{var}', f';{varNames[iVar]};Efficiency', nBins, np.array(binEdges)))
        for iBin in range(0, nBins):
            effSt1  = histDataEffSt[var][0].GetBinContent(iBin+1)
            effSt2  = histDataEffSt[var][1].GetBinContent(iBin+1)
            effSt3  = histDataEffSt[var][2].GetBinContent(iBin+1)
            effSt45 = histDataEffSt[var][3].GetBinContent(iBin+1)
            effMchTrkEff = effSt1 * effSt2 * effSt3 * effSt45
            histDataMchTrkEff[iVar].SetBinContent(iBin+1, effMchTrkEff)
            histDataMchTrkEff[iVar].SetBinError(iBin+1, 0)

    # Plot results
    lineUnityEta = ROOT.TLine(2.5, 1, 4, 1)
    lineUnityEta.SetLineColor(ROOT.kGray+1)
    lineUnityEta.SetLineStyle(ROOT.kDashed)

    lineUnityPt = ROOT.TLine(0, 1, 65, 1)
    lineUnityPt.SetLineColor(ROOT.kGray+1)
    lineUnityPt.SetLineStyle(ROOT.kDashed)

    lineUnityPhi = ROOT.TLine(-ROOT.TMath.Pi(), 1, ROOT.TMath.Pi(), 1)
    lineUnityPhi.SetLineColor(ROOT.kGray+1)
    lineUnityPhi.SetLineStyle(ROOT.kDashed)

    legendMchTrkEff = ROOT.TLegend(0.20, 0.78, 0.40, 0.93, " ", "brNDC")
    SetLegend(legendMchTrkEff)
    legendMchTrkEff.AddEntry(histDataMchTrkEff[0], "Data", "L")
    legendMchTrkEff.AddEntry(histMcMchTrkEff[0], "MC", "L")

    # Eta dependence
    canvasEtaStEff = ROOT.TCanvas("canvasEtaStEff", "", 1200, 1200)
    canvasEtaStEff.Divide(2, 2)
    for iSt, st in enumerate(stations):
        canvasEtaStEff.cd(iSt+1)
        histMcEffSt["Eta"][iSt].GetYaxis().SetRangeUser(0.98, 1.01)
        SetHistStat(histMcEffSt["Eta"][iSt], 20, ROOT.kAzure+2)
        histMcEffSt["Eta"][iSt].Draw("H")
        SetHistStat(histDataEffSt["Eta"][iSt], 20, ROOT.kRed+1)
        histDataEffSt["Eta"][iSt].Draw("H SAME")
    canvasEtaStEff.Update()

    canvasEtaMchTrkEff = ROOT.TCanvas("canvasEtaMchTrkEff", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    histMcMchTrkEff[0].GetYaxis().SetRangeUser(0.98, 1.01)
    SetHistStat(histMcMchTrkEff[0], 20, ROOT.kAzure+2)
    histMcMchTrkEff[0].Draw("H")
    SetHistStat(histDataMchTrkEff[0], 20, ROOT.kRed+1)
    histDataMchTrkEff[0].Draw("H SAME")
    legendMchTrkEff.Draw("SAME")
    lineUnityEta.Draw()
    canvasEtaMchTrkEff.Update()

    # Pt dependence
    canvasPt = ROOT.TCanvas("canvasPt", "", 1200, 1200)
    canvasPt.Divide(2, 2)
    for iSt, st in enumerate(stations):
        canvasPt.cd(iSt+1)
        histMcEffSt["Pt"][iSt].GetYaxis().SetRangeUser(0.98, 1.01)
        SetHistStat(histMcEffSt["Pt"][iSt], 20, ROOT.kAzure+2)
        histMcEffSt["Pt"][iSt].Draw("H")
        SetHistStat(histDataEffSt["Pt"][iSt], 20, ROOT.kRed+1)
        histDataEffSt["Pt"][iSt].Draw("H SAME")
    canvasPt.Update()

    canvasPtMchTrkEff = ROOT.TCanvas("canvasPtMchTrkEff", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    histMcMchTrkEff[1].GetXaxis().SetRangeUser(0, 65)
    histMcMchTrkEff[1].GetYaxis().SetRangeUser(0.98, 1.01)
    SetHistStat(histMcMchTrkEff[1], 20, ROOT.kAzure+2)
    histMcMchTrkEff[1].Draw("H")
    SetHistStat(histDataMchTrkEff[1], 20, ROOT.kRed+1)
    histDataMchTrkEff[1].Draw("H SAME")
    legendMchTrkEff.Draw("SAME")
    lineUnityPt.Draw()
    canvasPtMchTrkEff.Update()

    # Phi dependence
    canvasPhi = ROOT.TCanvas("canvasPhi", "", 1200, 1200)
    canvasPhi.Divide(2, 2)
    for iSt, st in enumerate(stations):
        canvasPhi.cd(iSt+1)
        histMcEffSt["Phi"][iSt].GetYaxis().SetRangeUser(0.98, 1.01)
        SetHistStat(histMcEffSt["Phi"][iSt], 20, ROOT.kAzure+2)
        histMcEffSt["Phi"][iSt].Draw("H")
        SetHistStat(histDataEffSt["Phi"][iSt], 20, ROOT.kRed+1)
        histDataEffSt["Phi"][iSt].Draw("H SAME")
    canvasPhi.Update()

    canvasPhiMchTrkEff = ROOT.TCanvas("canvasPhiMchTrkEff", "", 800, 600)
    ROOT.gStyle.SetOptStat(False)
    histMcMchTrkEff[2].GetYaxis().SetRangeUser(0.98, 1.01)
    SetHistStat(histMcMchTrkEff[2], 20, ROOT.kAzure+2)
    histMcMchTrkEff[2].Draw("H")
    SetHistStat(histDataMchTrkEff[2], 20, ROOT.kRed+1)
    histDataMchTrkEff[2].Draw("H SAME")
    legendMchTrkEff.Draw("SAME")
    lineUnityPhi.Draw()
    canvasPhiMchTrkEff.Update()

    histCorrMap = []
    for iVar, var in enumerate(vars):
        histCorrMap.append(histDataMchTrkEff[iVar].Clone(f'histCorrMap_{var}'))
        histCorrMap[iVar].Divide(histMcMchTrkEff[iVar])
    
    input()

    print(f'[INFO] Writing output to {config["outputs"]["fOutCorrMap"]} ...')
    fOutCorrMap = ROOT.TFile(config["outputs"]["fOutCorrMap"], "RECREATE")
    for iVar, var in enumerate(vars):
        histCorrMap[iVar].Write()
        histDataMchTrkEff[iVar].Write()
        histMcMchTrkEff[iVar].Write()
    fOutCorrMap.Close()

    canvasEtaMchTrkEff.SaveAs("figures/mch_trk_eff/eff_eta.pdf")
    canvasPtMchTrkEff.SaveAs("figures/mch_trk_eff/eff_pt.pdf")
    canvasPhiMchTrkEff.SaveAs("figures/mch_trk_eff/eff_phi.pdf")

    exit()

def funcSt123Eff(histXY, histX0, hist0Y):
    histEffX = histXY.Clone("histEffX")
    histMcEtaN0YXY = hist0Y.Clone("histMcEtaN0YXY")
    histMcEtaN0YXY.Add(histXY)
    histEffX.Divide(histMcEtaN0YXY)

    histEffY = histXY.Clone("histEffY")
    histMcEtaNX0XY = histX0.Clone("histMcEtaNX0XY")
    histMcEtaNX0XY.Add(histXY)
    histEffY.Divide(histMcEtaNX0XY)

    nBins = histXY.GetNbinsX()
    binEdges = [histXY.GetXaxis().GetBinLowEdge(iBin+1) for iBin in range(histXY.GetNbinsX())]
    binEdges.append(histXY.GetXaxis().GetBinUpEdge(histXY.GetNbinsX()))
    histEffStX = ROOT.TH1D(f'{histXY.GetName}_st', "", nBins, np.array(binEdges))

    for iBin in range(0, nBins):
        effX = histEffX.GetBinContent(iBin+1)
        effY = histEffY.GetBinContent(iBin+1)
        effSt = 1. - ((1. - effX) * (1. - effY))
        #print(f'[{histEffStX.GetXaxis().GetBinCenter(iBin+1)}] effX = {effX} ; effY = {effY} ; eff st1 = {effSt}')
        histEffStX.SetBinContent(iBin+1, effSt)
        histEffStX.SetBinError(iBin+1, 0)

    return histEffX, histEffY, histEffStX

def funcSt45Eff(hist78, hist70, hist08, hist910, hist90, hist010):
    # station 4
    histEff7 = hist78.Clone("histEff7")
    histMcEtaN0878 = hist08.Clone("histMcEtaN0878")
    histMcEtaN0878.Add(hist78)
    histEff7.Divide(histMcEtaN0878)

    histEff8 = hist78.Clone("histEff8")
    histMcEtaN7078 = hist70.Clone("histMcEtaN7078")
    histMcEtaN7078.Add(hist78)
    histEff8.Divide(histMcEtaN7078)

    # station 5
    histEff9 = hist910.Clone("histEff9")
    histMcEtaN010910 = hist010.Clone("histMcEtaN010910")
    histMcEtaN010910.Add(hist910)
    histEff9.Divide(histMcEtaN010910)

    histEff10 = hist910.Clone("histEff10")
    histMcEtaN9010 = hist90.Clone("histMcEtaN9010")
    histMcEtaN9010.Add(hist910)
    histEff10.Divide(histMcEtaN9010)

    nBins = hist78.GetNbinsX()
    binEdges = [hist78.GetXaxis().GetBinLowEdge(iBin+1) for iBin in range(hist78.GetNbinsX())]
    binEdges.append(hist78.GetXaxis().GetBinUpEdge(hist78.GetNbinsX()))
    histEffSt45 = ROOT.TH1D(f'{hist78.GetName}_st', "", nBins, np.array(binEdges))

    # Prod_i{e_i} + sum_i{(1 - e_i)*prod_j{ej}}
    for iBin in range(0, nBins):
        eff7 = histEff7.GetBinContent(iBin+1)
        eff8 = histEff8.GetBinContent(iBin+1)
        eff9 = histEff9.GetBinContent(iBin+1)
        eff10 = histEff10.GetBinContent(iBin+1)

        prod78910 = eff7 * eff8 * eff9 * eff10
        prod8910  = eff8 * eff9 * eff10
        prod7910  = eff7 * eff9 * eff10
        prod7810  = eff7 * eff8 * eff10
        prod789   = eff7 * eff8 * eff9

        sum = 0
        sum += (1 - eff7) * prod8910
        sum += (1 - eff8) * prod7910
        sum += (1 - eff9) * prod7810
        sum += (1 - eff10) * prod789

        effSt45 = prod78910 + sum

        histEffSt45.SetBinContent(iBin+1, effSt45)
        histEffSt45.SetBinError(iBin+1, 0)

    #return histEff7, histEff8, histEff9, histEff10, histEffSt45
    return histEffSt45



if __name__ == '__main__':
    main()