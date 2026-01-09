import sys
import argparse
import yaml
import ROOT
import numpy as np
import array

sys.path.append('../utils')
from plot_library import LoadStyle, SetHistStat, SetHistSyst, SetGraStat, SetGraSyst, SetLegend
from fast_analysis_utils import ExtractFromYaml

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--run", help="Compute Axe", action="store_true")
    parser.add_argument("--scan", help="Scan of Acxe vs Run", action="store_true")
    args = parser.parse_args()
    

    with open(args.cfgFileName, 'r') as yml_cfg:
        inputCfg = yaml.load(yml_cfg, yaml.FullLoader)

    if args.run:
        axe(inputCfg)
    
    if args.scan:
        scan(inputCfg)

def axe(config):
    """
    function to compute the axe
    """

    pathIn = config["inputs"]["path"]
    production = config["inputs"]["production"]
    prefixInGen = config["inputs"]["prefixGen"]
    prefixInRec = config["inputs"]["prefixRec"]
    assocType = config["inputs"]["assocType"]

    fInNameGen = f'{pathIn}/{production}/{prefixInGen}_{assocType}.root'
    fInNameRec = f'{pathIn}/{production}/{prefixInRec}_{assocType}.root'

    pathOut = config["outputs"]["path"]
    prefixOut = config["outputs"]["prefix"]

    corrFactor = 1. # this factor takes into account the running separately of generated and reconstructed
    if (fInNameGen == fInNameRec):
        print(f'[INFO] Opening {fInNameGen}...')
        fInGen = ROOT.TFile(f'{fInNameGen}', "READ")
        fInRec = fInGen
    else:
        print(f'[INFO] Opening {fInNameGen}...')
        print(f'[INFO] Opening {fInNameRec}...')
        fInGen = ROOT.TFile(f'{fInNameGen}', "READ")
        fInRec = ROOT.TFile(f'{fInNameRec}', "READ")

        hlistEventsMcGen = fInGen.Get("analysis-event-selection/output")
        histEventsMcGen = (hlistEventsMcGen.FindObject("EventsMC")).FindObject("MCVtxX")

        hlistEventsMcRec = fInRec.Get("analysis-event-selection/output")
        histEventsMcRec = (hlistEventsMcRec.FindObject("EventsMC")).FindObject("MCVtxX")

        corrFactor = (histEventsMcRec.GetEntries()) / (histEventsMcGen.GetEntries())

        print(f'MC evts. Gen = {histEventsMcGen.GetEntries()} ; MC evts. Rec = {histEventsMcRec.GetEntries()} ; corr. Factor = {corrFactor}')

    # Generated distributions
    hlistGen = fInGen.Get(config["inputs"]["hlistGen"])
    strListGen = config["inputs"]["listGen"]
    strHistGen = config["inputs"]["histGen"]
    sigGen = config["inputs"]["sigGen"]

    minPtBinInt = config["inputs"]["minPtBinInt"]
    maxPtBinInt = config["inputs"]["maxPtBinInt"]
    minRapBinInt = config["inputs"]["minRapBinInt"]
    maxRapBinInt = config["inputs"]["maxRapBinInt"]
    minCentFT0CBinInt = config["inputs"]["minCentFT0CBinInt"]
    maxCentFT0CBinInt = config["inputs"]["maxCentFT0CBinInt"]

    histPtRapCentGen_original = (hlistGen.FindObject(f'{strListGen}_{sigGen}')).FindObject(strHistGen)
    canvasGen = ROOT.TCanvas("canvasGen", "", 800, 600)
    canvasGen.Divide(3, 4)
    canvasGen.cd(1)
    histPtRapCentGen_original.Draw("COLZ")

    # Copying histo 'histPtRapCentGen_original' in a new one 'histPtRapGen' with tigher ranges
    binXmin = histPtRapCentGen_original.GetXaxis().FindBin(minPtBinInt)
    binXmax = histPtRapCentGen_original.GetXaxis().FindBin(maxPtBinInt - 1e-6)  # To avoid to include the next bin

    binYmin = histPtRapCentGen_original.GetYaxis().FindBin(minRapBinInt)
    binYmax = histPtRapCentGen_original.GetYaxis().FindBin(maxRapBinInt - 1e-6)

    binZmin = histPtRapCentGen_original.GetZaxis().FindBin(minCentFT0CBinInt)
    binZmax = histPtRapCentGen_original.GetZaxis().FindBin(maxCentFT0CBinInt - 1e-6)

    nBinsX = binXmax - binXmin + 1
    nBinsY = binYmax - binYmin + 1
    nBinsZ = binZmax - binZmin + 1

    xMin = histPtRapCentGen_original.GetXaxis().GetBinLowEdge(binXmin)
    xMax = histPtRapCentGen_original.GetXaxis().GetBinUpEdge(binXmax)

    yMin = histPtRapCentGen_original.GetYaxis().GetBinLowEdge(binYmin)
    yMax = histPtRapCentGen_original.GetYaxis().GetBinUpEdge(binYmax)

    zMin = histPtRapCentGen_original.GetZaxis().GetBinLowEdge(binZmin)
    zMax = histPtRapCentGen_original.GetZaxis().GetBinUpEdge(binZmax)

    histPtRapGen = ROOT.TH3D("histPtRapGen","histPtRapGen",nBinsX, xMin, xMax,nBinsY, yMin, yMax,nBinsZ, zMin, zMax)

    for ix_new, ix_original in enumerate(range(binXmin, binXmax + 1), start=1):
        for iy_new, iy_original in enumerate(range(binYmin, binYmax + 1), start=1):
            for iz_new, iz_original in enumerate(range(binZmin, binZmax + 1), start=1):

                val = histPtRapCentGen_original.GetBinContent(ix_original, iy_original, iz_original)
                err = histPtRapCentGen_original.GetBinError(ix_original, iy_original, iz_original)

                histPtRapGen.SetBinContent(ix_new, iy_new, iz_new, val)
                histPtRapGen.SetBinError(ix_new, iy_new, iz_new, err)

    canvasGen.cd(2)
    ROOT.gPad.SetLogy(True)
    histPtRapGen.Draw("COLZ")
    canvasGen.Update()

    print(f'GEN --> ORIGINAL Histo: binXmin = {binXmin} ; binXmax = {binXmax} ; binYmin = {binYmin} ; binYmax = {binYmax} ; binZmin = {binZmin} ; binZmax = {binZmax} ')
    print(f'GEN --> NEW Histo: nBinsX = {nBinsX} ; nBinsY = {nBinsY} ; nBinsZ = {nBinsZ}')
    print(f'GEN --> NEW Histo: xMin = {xMin} ; xMax = {xMax}')
    print(f'GEN --> NEW Histo: yMin = {yMin} ; yMax = {yMax}')
    print(f'GEN --> NEW Histo: zMin = {zMin} ; zMax = {zMax}')

    histPtRapCentGenFor1D = histPtRapGen.Clone("histPtRapCentGenFor1D")

    histPtGen = histPtRapCentGenFor1D.ProjectionX(f'{sigGen}_Pt')
    histRapGen = histPtRapCentGenFor1D.ProjectionY(f'{sigGen}_Rap')
    histCentGen = histPtRapCentGenFor1D.ProjectionZ(f'{sigGen}_Z')

    ptBinsRun3 = array.array('d', config["inputs"]["ptBins"])
    rapBinsRun3 = array.array('d', config["inputs"]["rapBins"])
    centBinsRun3 = array.array('d', config["inputs"]["centFT0CBins"])

    histPtRebinGen = histPtGen.Rebin(len(ptBinsRun3) - 1, f'{histPtGen.GetName()}_rebin', ptBinsRun3)
    histRapRebinGen = histRapGen.Rebin(len(rapBinsRun3) - 1, f'{histRapGen.GetName()}_rebin', rapBinsRun3)
    histCentRebinGen = histCentGen.Rebin(len(centBinsRun3) - 1, f'{histCentGen.GetName()}_rebin', centBinsRun3)

    canvasGen.cd(3)
    histPtRapCentGenFor1D.Draw("COLZ")

    canvasGen.cd(4)
    ROOT.gPad.SetLogy(True)
    histPtGen.Draw()

    canvasGen.cd(5)
    ROOT.gPad.SetLogy(True)
    histPtRebinGen.Draw()

    canvasGen.cd(7)
    ROOT.gPad.SetLogy(True)
    histRapGen.Draw()

    canvasGen.cd(8)
    ROOT.gPad.SetLogy(True)
    histRapRebinGen.Draw()

    canvasGen.cd(10)
    ROOT.gPad.SetLogy(True)
    histCentGen.Draw()

    canvasGen.cd(11)
    ROOT.gPad.SetLogy(True)
    histCentRebinGen.Draw()

    canvasGen.Update()
    input()

    # Reconstructed distributions
    hlistRec = fInRec.Get(config["inputs"]["hlistRec"])
    listRec = config["inputs"]["listRec"]
    histRec = config["inputs"]["histRec"]
    cuts = config["inputs"]["cuts"]
    sigRec = config["inputs"]["sigRec"]

    histRecs = []
    for iCut, cut in enumerate(cuts):
        histRecs.append((hlistRec.FindObject(f'{listRec}_{cut}_{sigRec}')).FindObject(f'{histRec}'))

    histPtRapRecs_original = []
    histPtRapCentRecs = []
    histPtRapCentRecsFor1D = []
    histPtRecs = []
    histRapRecs = []
    histCentRecs = []
    for iCut, histRec in enumerate(histRecs):
        histRec.SetName(f'THnSparse_{sigRec}_{cuts[iCut]}') # Just to avoid memory leak

        histPtRapRecs_original.append(histRec.Projection(1, 2, 3, f'{listRec}_{cuts[iCut]}_{sigRec}_PtRapCent'))
        histPtRapRecs_original[iCut].SetName(f'{listRec}_{cuts[iCut]}_{sigRec}_PtRapCent')

        canvasRec = ROOT.TCanvas("canvasRec", "", 800, 600)
        canvasRec.Divide(3, 4)
        canvasRec.cd(1)
        histPtRapRecs_original[iCut].Draw("COLZ")

        # Copying histo 'histPtRapRecs_original' in a new one 'histPtRapCentRecs' with tigher ranges
        binXmin = histPtRapRecs_original[iCut].GetXaxis().FindBin(minPtBinInt)
        binXmax = histPtRapRecs_original[iCut].GetXaxis().FindBin(maxPtBinInt - 1e-6)

        binYmin = histPtRapRecs_original[iCut].GetYaxis().FindBin(minRapBinInt)
        binYmax = histPtRapRecs_original[iCut].GetYaxis().FindBin(maxRapBinInt - 1e-6)

        binZmin = histPtRapRecs_original[iCut].GetZaxis().FindBin(minCentFT0CBinInt)
        binZmax = histPtRapRecs_original[iCut].GetZaxis().FindBin(maxCentFT0CBinInt - 1e-6)

        nBinsX = binXmax - binXmin + 1
        nBinsY = binYmax - binYmin + 1
        nBinsZ = binZmax - binZmin + 1

        xMin = histPtRapRecs_original[iCut].GetXaxis().GetBinLowEdge(binXmin)
        xMax = histPtRapRecs_original[iCut].GetXaxis().GetBinUpEdge(binXmax)

        yMin = histPtRapRecs_original[iCut].GetYaxis().GetBinLowEdge(binYmin)
        yMax = histPtRapRecs_original[iCut].GetYaxis().GetBinUpEdge(binYmax)

        zMin = histPtRapRecs_original[iCut].GetZaxis().GetBinLowEdge(binZmin)
        zMax = histPtRapRecs_original[iCut].GetZaxis().GetBinUpEdge(binZmax)

        histPtRapCentRecs.append(ROOT.TH3D("histPtRapCentRecs","histPtRapCentRecs",nBinsX, xMin, xMax,nBinsY, yMin, yMax,nBinsZ, zMin, zMax))

        for ix_new, ix_original in enumerate(range(binXmin, binXmax + 1), start=1):
            for iy_new, iy_original in enumerate(range(binYmin, binYmax + 1), start=1):
                for iz_new, iz_original in enumerate(range(binZmin, binZmax + 1), start=1):

                    val = histPtRapRecs_original[iCut].GetBinContent(ix_original, iy_original, iz_original)
                    err = histPtRapRecs_original[iCut].GetBinError(ix_original, iy_original, iz_original)

                    histPtRapCentRecs[iCut].SetBinContent(ix_new, iy_new, iz_new, val)
                    histPtRapCentRecs[iCut].SetBinError(ix_new, iy_new, iz_new, err)

        canvasRec.cd(2)
        ROOT.gPad.SetLogy(True)
        histPtRapCentRecs[iCut].Draw("COLZ")
        canvasRec.Update()

        print(f'RECO--> ORIGINAL Histo: binXmin = {binXmin} ; binXmax = {binXmax} ; binYmin = {binYmin} ; binYmax = {binYmax} ; binZmin = {binZmin} ; binZmax = {binZmax} ')
        print(f'RECO--> NEW Histo: nBinsX = {nBinsX} ; nBinsY = {nBinsY} ; nBinsZ = {nBinsZ}')
        print(f'RECO--> NEW Histo: xMin = {xMin} ; xMax = {xMax}')
        print(f'RECO--> NEW Histo: yMin = {yMin} ; yMax = {yMax}')
        print(f'RECO--> NEW Histo: zMin = {zMin} ; zMax = {zMax}')

        histPtRapCentRecsFor1D = histPtRapCentRecs[iCut].Clone("histPtRapCentRecsFor1D")

        histPtRecs.append(histPtRapCentRecsFor1D.ProjectionX(f'{listRec}_{cuts[iCut]}_{sigRec}_Pt'))
        histRapRecs.append(histPtRapCentRecsFor1D.ProjectionY(f'{listRec}_{cuts[iCut]}_{sigRec}_Rap'))
        histCentRecs.append(histPtRapCentRecsFor1D.ProjectionZ(f'{listRec}_{cuts[iCut]}_{sigRec}_CentFT0C'))

    histPtRebinRecs = []
    for iCut, histPtRec in enumerate(histPtRecs):
        histPtRebinRecs.append(histPtRec.Rebin(len(ptBinsRun3) - 1, f'{histPtRec.GetName()}_rebin', ptBinsRun3))

    histRapRebinRecs = []
    for iCut, histRapRec in enumerate(histRapRecs):
        histRapRebinRecs.append(histRapRec.Rebin(len(rapBinsRun3) - 1, f'{histRapRec.GetName()}_rebin', rapBinsRun3))

    histCentRebinRecs = []
    for iCut, histCentRec in enumerate(histCentRecs):
        histCentRebinRecs.append(histCentRec.Rebin(len(centBinsRun3) - 1, f'{histCentRec.GetName()}_rebin', centBinsRun3))

    canvasRec.cd(3)
    histPtRapCentRecsFor1D.Draw("COLZ")

    canvasRec.cd(4)
    ROOT.gPad.SetLogy(True)
    histRapRecs[0].Draw("EP")

    canvasRec.cd(5)
    ROOT.gPad.SetLogy(True)
    histRapRebinRecs[0].Draw("EP")

    canvasRec.cd(7)
    ROOT.gPad.SetLogy(True)
    histPtRecs[0].Draw("EP")

    canvasRec.cd(8)
    ROOT.gPad.SetLogy(True)
    histPtRebinRecs[0].Draw("EP")

    canvasRec.cd(10)
    ROOT.gPad.SetLogy(True)
    histCentRecs[0].Draw("EP")

    canvasRec.cd(11)
    ROOT.gPad.SetLogy(True)
    histCentRebinRecs[0].Draw("EP")

    canvasRec.Update()
    input()


    print("Compute the ratio")
    nBinsX = 25
    binEdgesX = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 12, 14, 16, 18, 20]
    nBinsY = 15
    binEdgesY = [2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4]
    nBinsZ = 9
    binEdgesY = [0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 90.0]

    ######################################
    # Rebinning generated for 2D Axe map #
    ######################################
    '''
    rapRequredBinWidth = 0.1
    minRapGen = histPtRapGen.GetYaxis().GetBinLowEdge(1)
    maxRapGen = histPtRapGen.GetYaxis().GetBinLowEdge(histPtRapGen.GetNbinsY()) + histPtRapGen.GetYaxis().GetBinWidth(histPtRapGen.GetNbinsY())
    print(f'[INFO] Scalig factor for Rapidity = {rapRequredBinWidth * (histPtRapGen.GetNbinsY() / (maxRapGen - minRapGen))}')

    histPtRapMapGenLowPt = histPtRapGen.Clone(f'{histPtRapGen.GetName()}_low_pt')
    histPtRapMapGenLowPt.RebinX(2)
    histPtRapMapGenLowPt.RebinY(10)
    histPtRapMapGenLowPt.GetXaxis().SetRangeUser(0., 20.)
    histPtRapMapGenLowPt.GetYaxis().SetRangeUser(2.5, 4)

    histPtRapMapGenHighPt = histPtRapGen.Clone(f'{histPtRapGen.GetName()}_high_pt')
    histPtRapMapGenHighPt.RebinX(8)
    histPtRapMapGenHighPt.RebinY(10)
    histPtRapMapGenHighPt.GetXaxis().SetRangeUser(0., 20.)
    histPtRapMapGenHighPt.GetYaxis().SetRangeUser(2.5, 4)

    # There could be an offset in rapidity
    rapGenOffset = histPtRapMapGenLowPt.GetYaxis().FindBin(2.55)
    print(f'[INFO] Generated rapidity offset: {rapGenOffset}')

    histPtRapMapGen = ROOT.TH2D(f'{histPtRapGen.GetName()}_map', "", nBinsX, np.array(binEdgesX), nBinsY, np.array(binEdgesY))
    for xBin in range(0, nBinsX):
        for yBin in range(0, nBinsY):
            if xBin < 20:
                if __debug__:
                    debugCtrBinX1 = histPtRapMapGen.GetXaxis().GetBinCenter(xBin+1)
                    debugCtrBinX2 = histPtRapMapGenLowPt.GetXaxis().GetBinCenter(xBin+1)
                    debugCtrBinY1 = histPtRapMapGen.GetYaxis().GetBinCenter(yBin+1)
                    debugCtrBinY2 = histPtRapMapGenLowPt.GetYaxis().GetBinCenter(rapGenOffset+yBin)
                    print(f"pT bins [{debugCtrBinX1:.3f}, {debugCtrBinX2:.3f}] ; Rap bins [{debugCtrBinY1:.3f}, {debugCtrBinY2:.3f}]")
                histPtRapMapGen.SetBinContent(xBin+1, yBin+1, histPtRapMapGenLowPt.GetBinContent(xBin+1, rapGenOffset+yBin))
                histPtRapMapGen.SetBinError(xBin+1, yBin+1, histPtRapMapGenLowPt.GetBinError(xBin+1, rapGenOffset+yBin))
            else:
                if __debug__:
                    debugCtrBinX1 = histPtRapMapGen.GetXaxis().GetBinCenter(xBin+1)
                    debugCtrBinX2 = histPtRapMapGenHighPt.GetXaxis().GetBinCenter(xBin+1-15)
                    debugCtrBinY1 = histPtRapMapGen.GetYaxis().GetBinCenter(yBin+1)
                    debugCtrBinY2 = histPtRapMapGenHighPt.GetYaxis().GetBinCenter(rapGenOffset+yBin)
                    print(f"pT bins [{debugCtrBinX1:.3f}, {debugCtrBinX2:.3f}] ; rap bins [{debugCtrBinY1:.3f}, {debugCtrBinY2:.3f}]")
                histPtRapMapGen.SetBinContent(xBin+1, yBin+1, histPtRapMapGenHighPt.GetBinContent(xBin+1-15, rapGenOffset+yBin))
                histPtRapMapGen.SetBinError(xBin+1, yBin+1, histPtRapMapGenHighPt.GetBinError(xBin+1-15, rapGenOffset+yBin))
    '''
    ##########################################
    # Rebinning reconstructed for 2D Axe map #
    ##########################################
    '''
    histPtRapMapRecsLowPt = []
    histPtRapMapRecsHighPt = []

    for iCut, histPtRapRec in enumerate(histPtRapCentRecs):
        histPtRapMapRecsLowPt.append(histPtRapRec.Clone(f'{histPtRapRec.GetName()}_low_pt'))
        histPtRapMapRecsLowPt[0].RebinX(2)
        histPtRapMapRecsLowPt[0].RebinY(4)

        histPtRapMapRecsHighPt.append(histPtRapRec.Clone(f'{histPtRapRec.GetName()}_high_pt'))
        histPtRapMapRecsHighPt[0].RebinX(8)
        histPtRapMapRecsHighPt[0].RebinY(4)

    # There could be an offset in rapidity
    rapRecOffset = histPtRapMapRecsLowPt[0].GetYaxis().FindBin(2.55)
    print(f'[INFO] Reconstructed rapidity offset: {rapGenOffset}')

    histPtRapMapRecs = []
    for iCut, histPtRapRec in enumerate(histPtRapCentRecs):
        histPtRapMapRecs.append(ROOT.TH2D(f'{histPtRapRec.GetName()}_map', "", nBinsX, np.array(binEdgesX), nBinsY, np.array(binEdgesY)))
        for xBin in range(0, nBinsX):
            for yBin in range(0, nBinsY):
                if xBin < 20:
                    histPtRapMapRecs[iCut].SetBinContent(xBin+1, yBin+1, histPtRapMapRecsLowPt[iCut].GetBinContent(xBin+1, rapRecOffset+yBin))
                    histPtRapMapRecs[iCut].SetBinError(xBin+1, yBin+1, histPtRapMapRecsLowPt[iCut].GetBinError(xBin+1, rapRecOffset+yBin))
                    if __debug__:
                        debugCtrBinX1 = histPtRapMapRecs[iCut].GetXaxis().GetBinCenter(xBin+1)
                        debugCtrBinX2 = histPtRapMapRecsLowPt[iCut].GetXaxis().GetBinCenter(xBin+1)
                        debugCtrBinY1 = histPtRapMapRecs[iCut].GetYaxis().GetBinCenter(yBin+1)
                        debugCtrBinY2 = histPtRapMapRecsLowPt[iCut].GetYaxis().GetBinCenter(rapRecOffset+yBin)
                        print(f"pT bins [{debugCtrBinX1:.3f}, {debugCtrBinX2:.3f}] ; Rap bins [{debugCtrBinY1:.3f}, {debugCtrBinY2:.3f}]")
                else:
                    histPtRapMapRecs[iCut].SetBinContent(xBin+1, yBin+1, histPtRapMapRecsHighPt[iCut].GetBinContent(xBin+1-15, rapRecOffset+yBin))
                    histPtRapMapRecs[iCut].SetBinError(xBin+1, yBin+1, histPtRapMapRecsHighPt[iCut].GetBinError(xBin+1-15, rapRecOffset+yBin))
                    if __debug__:
                        debugCtrBinX1 = histPtRapMapRecs[iCut].GetXaxis().GetBinCenter(xBin+1)
                        debugCtrBinX2 = histPtRapMapRecsHighPt[iCut].GetXaxis().GetBinCenter(xBin+1-15)
                        debugCtrBinY1 = histPtRapMapRecs[iCut].GetYaxis().GetBinCenter(yBin+1)
                        debugCtrBinY2 = histPtRapMapRecsHighPt[iCut].GetYaxis().GetBinCenter(rapRecOffset+yBin)
                        print(f"pT bins [{debugCtrBinX1:.3f}, {debugCtrBinX2:.3f}] ; Rap bins [{debugCtrBinY1:.3f}, {debugCtrBinY2:.3f}]")
    '''
    print(f'[WARNING] scale the generated for job efficiency (corr. factor = {corrFactor})')
    #histPtRapMapGen.Scale(corrFactor)
    histPtRebinGen.Scale(corrFactor)
    histRapRebinGen.Scale(corrFactor)
    histCentRebinGen.Scale(corrFactor)

    '''
    histPtRapAxes = []
    for iCut, histPtRapRec in enumerate(histPtRapCentRecs):
        histPtRapAxes.append(ROOT.TH2D("Axe_map", "", nBinsX, np.array(binEdgesX), nBinsY, np.array(binEdgesY)))
        histPtRapAxes[iCut].Divide(histPtRapMapRecs[iCut], histPtRapMapGen, 1.0, 1.0, "B")

    canvasAxeTest = ROOT.TCanvas("canvasAxeTest", "", 800, 600)
    canvasAxeTest.Divide(2, 2)

    canvasAxeTest.cd(1)
    histPtRapMapGen.Draw("COLZ")

    canvasAxeTest.cd(2)
    histPtRapMapRecs[0].Draw("COLZ")

    canvasAxeTest.cd(3)
    histPtRapAxes[0].Draw("COLZ")

    canvasAxeTest.Update()
    '''

    histPtRebinAxes = []
    for iCut, histPtRebinRec in enumerate(histPtRebinRecs):
        histPtRebinAxes.append(histPtRebinRec.Clone(f'{histPtRebinRec.GetName()}_rebin_Axe'))
        histPtRebinAxes[iCut].Divide(histPtRebinAxes[iCut], histPtRebinGen, 1.0, 1.0, "B")

    histRapRebinAxes = []
    for iCut, histRapRebinRec in enumerate(histRapRebinRecs):
        histRapRebinAxes.append(histRapRebinRec.Clone(f'{histRapRebinRec.GetName()}_rebin_Axe'))
        histRapRebinAxes[iCut].Divide(histRapRebinAxes[iCut], histRapRebinGen, 1.0, 1.0, "B")

    histCentRebinAxes = []
    for iCut, histCentRebinRec in enumerate(histCentRebinRecs):
        histCentRebinAxes.append(histCentRebinRec.Clone(f'{histCentRebinRec.GetName()}_rebin_Axe'))
        histCentRebinAxes[iCut].Divide(histCentRebinAxes[iCut], histCentRebinGen, 1.0, 1.0, "B")

    canvasAxe = ROOT.TCanvas("canvasAxe", "", 800, 600)
    canvasAxe.Divide(2, 2)

    #canvasAxe.cd(1)
    #histPtRapAxes[0].Draw("COLZ")

    canvasAxe.cd(1)
    ROOT.gPad.SetLogy(False)
    histRapRebinAxes[0].Draw()
    
    canvasAxe.cd(2)
    ROOT.gPad.SetLogy(False)
    histPtRebinAxes[0].GetYaxis().SetRangeUser(0.001, 1.2)
    histPtRebinAxes[0].Draw()

    canvasAxe.cd(3)
    ROOT.gPad.SetLogy(False)
    histCentRebinAxes[0].GetYaxis().SetRangeUser(0.001, 1.2)
    histCentRebinAxes[0].Draw()

    canvasAxe.Update()

    nPtBins = len(ptBinsRun3) - 1
    nCentBins = len(centBinsRun3) - 1
    nRapBins = len(rapBinsRun3) - 1
    if (nCentBins <= 1) and (nPtBins > 1):
        print("x_min x_max val stat syst ")
        for iPt in range(len(ptBinsRun3)-1):
            x_min = histPtRebinAxes[0].GetBinLowEdge(iPt+1)
            x_max = histPtRebinAxes[0].GetBinLowEdge(iPt+1) + histPtRebinAxes[0].GetBinWidth(iPt+1)
            val = histPtRebinAxes[0].GetBinContent(iPt+1)
            stat = histPtRebinAxes[0].GetBinError(iPt+1)
            syst = 0.0
            print("{:3.2f} {:3.2f} {:6.5f} {:6.5f} {:6.5f} ".format(x_min, x_max, val, stat, syst))
    elif (nCentBins > 1) and (nPtBins <= 1):
        print("x_min x_max val stat syst ")
        for iCent in range(len(centBinsRun3)-1):
            x_min = histCentRebinAxes[0].GetBinLowEdge(iCent+1)
            x_max = histCentRebinAxes[0].GetBinLowEdge(iCent+1) + histCentRebinAxes[0].GetBinWidth(iCent+1)
            val = histCentRebinAxes[0].GetBinContent(iCent+1)
            stat = histCentRebinAxes[0].GetBinError(iCent+1)
            syst = 0.0
            print("{:3.2f} {:3.2f} {:6.5f} {:6.5f} {:6.5f} ".format(x_min, x_max, val, stat, syst))
    elif (nCentBins <= 1) and (nPtBins <= 1) and (nRapBins <= 1):
        print("x_min x_max val stat syst ")
        for iRap in range(len(rapBinsRun3)-1):
            x_min = histRapRebinAxes[0].GetBinLowEdge(iRap+1)
            x_max = histRapRebinAxes[0].GetBinLowEdge(iRap+1) + histRapRebinAxes[0].GetBinWidth(iRap+1)
            val = histRapRebinAxes[0].GetBinContent(iRap+1)
            stat = histRapRebinAxes[0].GetBinError(iRap+1)
            syst = 0.0
            print("{:3.2f} {:3.2f} {:6.5f} {:6.5f} {:6.5f} ".format(x_min, x_max, val, stat, syst))



    input()

    fOutNameTxt = f'{pathOut}/{production}/{prefixOut}_{assocType}_{sigGen}.txt'
    print(f'Saving in {fOutNameTxt}...')
    fOutTxt = open(f'{fOutNameTxt}', 'w')
    fOutTxt.write("x_min x_max val stat syst \n")
    if (nCentBins <= 1) and (nPtBins > 1):
        for iPt in range(len(ptBinsRun3)-1):
            x_min = histPtRebinAxes[0].GetBinLowEdge(iPt+1)
            x_max = histPtRebinAxes[0].GetBinLowEdge(iPt+1) + histPtRebinAxes[0].GetBinWidth(iPt+1)
            val = histPtRebinAxes[0].GetBinContent(iPt+1)
            stat = histPtRebinAxes[0].GetBinError(iPt+1)
            syst = 0.0
            if (iPt == (len(ptBinsRun3)-2)):
                fOutTxt.write(f"{x_min:3.2f} {x_max:3.2f} {val:6.5f} {stat:6.5f} {syst:6.5f}")
            else:
                fOutTxt.write(f"{x_min:3.2f} {x_max:3.2f} {val:6.5f} {stat:6.5f} {syst:6.5f} \n")
    elif (nCentBins > 1) and (nPtBins <= 1):
        for iCent in range(len(centBinsRun3)-1):
            x_min = histCentRebinAxes[0].GetBinLowEdge(iCent+1)
            x_max = histCentRebinAxes[0].GetBinLowEdge(iCent+1) + histCentRebinAxes[0].GetBinWidth(iCent+1)
            val = histCentRebinAxes[0].GetBinContent(iCent+1)
            stat = histCentRebinAxes[0].GetBinError(iCent+1)
            syst = 0.0
            if (iCent == (len(centBinsRun3)-2)):
                fOutTxt.write(f"{x_min:3.2f} {x_max:3.2f} {val:6.5f} {stat:6.5f} {syst:6.5f}")
            else:
                fOutTxt.write(f"{x_min:3.2f} {x_max:3.2f} {val:6.5f} {stat:6.5f} {syst:6.5f} \n")
    elif (nCentBins <= 1) and (nPtBins <= 1) and (nRapBins <= 1):
        for iRap in range(len(rapBinsRun3)-1):
            x_min = histRapRebinAxes[0].GetBinLowEdge(iRap+1)
            x_max = histRapRebinAxes[0].GetBinLowEdge(iRap+1) + histRapRebinAxes[0].GetBinWidth(iRap+1)
            val = histRapRebinAxes[0].GetBinContent(iRap+1)
            stat = histRapRebinAxes[0].GetBinError(iRap+1)
            syst = 0.0
            if (iRap == (len(rapBinsRun3)-2)):
                fOutTxt.write(f"{x_min:3.2f} {x_max:3.2f} {val:6.5f} {stat:6.5f} {syst:6.5f}")
            else:
                fOutTxt.write(f"{x_min:3.2f} {x_max:3.2f} {val:6.5f} {stat:6.5f} {syst:6.5f} \n")

    fOutNameRoot = f'{pathOut}/{production}/{prefixOut}_{assocType}_{sigGen}.root'
    print(f'Saving in {fOutNameRoot}...')
    fOutRoot = ROOT.TFile(f'{fOutNameRoot}', "RECREATE")
    histPtRapGen.Write("histPtRapGen")
    histPtRebinGen.Write("histPtRebinGen")
    histRapRebinGen.Write("histRapRebinGen")
    histCentRebinGen.Write("histCentRebinGen")
    for iCut, cut in enumerate(cuts):
        histPtRapCentRecs[iCut].Write(f'histPtRapRecs_{cut}')
        histPtRebinRecs[iCut].Write(f'histPtRebinRecs_{cut}')
        histRapRebinRecs[iCut].Write(f'histRapRebinRecs_{cut}')
        histCentRebinRecs[iCut].Write(f'histCentFT0CRebinRecs_{cut}')
        #histPtRapAxes[iCut].Write(f'histPtRapAxes_{cut}')
        histPtRebinAxes[iCut].Write(f'histPtRebinAxes_{cut}')
        histRapRebinAxes[iCut].Write(f'histRapRebinAxes_{cut}')
        histCentRebinAxes[iCut].Write(f'histCentFT0CRebinAxes_{cut}')

    fOutRoot.Close()
    if (fInNameGen == fInNameRec):
        fInGen.Close()
    else:
        fInGen.Close()
        fInRec.Close()
    
    exit()

def scan(config):
    """
    function to perform a Axe scan vs Run
    """
    LoadStyle()
    ROOT.gStyle.SetOptStat(False)

    path = config["scan"]["path"]
    fileName = config["scan"]["fileName"]

    with open(f'{path}/{config["scan"]["runList"]}', 'r') as file:
        runList = file.read().splitlines()

    histAxePerRun = ROOT.TH1D("histAxePerRun", ";Run;A#times#varepsilon", len(runList), 0, len(runList))

    for iRun, run in enumerate(runList):
        print(f'{path}/{run}/{fileName}.root')
        histAxePerRun.GetXaxis().SetBinLabel(iRun+1, run)
        if ROOT.gSystem.AccessPathName(f'{path}/{run}/{fileName}.root'):
            continue
        else:
            fIn = ROOT.TFile(f'{path}/{run}/{fileName}.root', "READ")
            histAxe = fIn.Get("histPtRebinAxes_matchedMchMid")
            histAxePerRun.SetBinContent(iRun+1, histAxe.GetBinContent(1))
            histAxePerRun.SetBinError(iRun+1, histAxe.GetBinError(1))

    canvasAxePerRun = ROOT.TCanvas("canvasAxePerRun", "", 800, 600)
    histAxePerRun.GetYaxis().SetRangeUser(0, 1)
    histAxePerRun.Draw("EP")
    canvasAxePerRun.Update()

    input()
    exit()

if __name__ == '__main__':
    main()