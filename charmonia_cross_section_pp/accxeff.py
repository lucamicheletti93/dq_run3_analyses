import sys
import argparse
import yaml
import ROOT
import numpy as np
import array

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--do_axe", help="Compute Axe", action="store_true")
    args = parser.parse_args()
    

    with open(args.cfgFileName, 'r') as yml_cfg:
        inputCfg = yaml.load(yml_cfg, yaml.FullLoader)

    if args.do_axe:
        axe(inputCfg)

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

    histPtRapGen = (hlistGen.FindObject(f'{strListGen}_{sigGen}')).FindObject(strHistGen)
    histPtRapGenFor1D = histPtRapGen.Clone("histPtRapGenFor1D")

    # For the 2D map the fill (pT, y) range can be used
    histPtRapGen.GetXaxis().SetRangeUser(0.0, 20.0)
    histPtRapGen.GetYaxis().SetRangeUser(2.5, 4.0)

    # For the 1D correction a subrange has to be defined
    histPtRapGenFor1D.GetXaxis().SetRangeUser(minPtBinInt, maxPtBinInt)
    histPtRapGenFor1D.GetYaxis().SetRangeUser(minRapBinInt, maxRapBinInt)

    histPtGen = histPtRapGenFor1D.ProjectionX(f'{sigGen}_Pt')
    histRapGen = histPtRapGenFor1D.ProjectionY(f'{sigGen}_Rap')

    ptBinsRun3 = array.array('d', config["inputs"]["ptBins"])
    rapBinsRun3 = array.array('d', config["inputs"]["rapBins"])

    histPtRebinGen = histPtGen.Rebin(len(ptBinsRun3) - 1, f'{histPtGen.GetName()}_rebin', ptBinsRun3)
    histRapRebinGen = histRapGen.Rebin(len(rapBinsRun3) - 1, f'{histRapGen.GetName()}_rebin', rapBinsRun3)
        
    canvasGen = ROOT.TCanvas("canvasGen", "", 800, 600)
    canvasGen.Divide(2, 2)

    canvasGen.cd(1)
    histPtRapGenFor1D.Draw("COLZ")

    canvasGen.cd(2)
    ROOT.gPad.SetLogy(True)
    histRapRebinGen.Draw()

    canvasGen.cd(3)
    ROOT.gPad.SetLogy(True)
    histPtRebinGen.Draw()

    canvasGen.cd(4)
    histPtRapGen.Draw("COLZ")

    canvasGen.Update()

    # Reconstructed distributions
    hlistRec = fInRec.Get(config["inputs"]["hlistRec"])
    listRec = config["inputs"]["listRec"]
    histRec = config["inputs"]["histRec"]
    cuts = config["inputs"]["cuts"]
    sigRec = config["inputs"]["sigRec"]

    histRecs = []
    for iCut, cut in enumerate(cuts):
        histRecs.append((hlistRec.FindObject(f'{listRec}_{cut}_{sigRec}')).FindObject(f'{histRec}'))
        #histRecs[iCut].GetAxis(1).SetRangeUser(minPtBinInt, maxPtBinInt)
        #histRecs[iCut].GetAxis(2).SetRangeUser(minRapBinInt, maxRapBinInt)

    histPtRapRecs = []
    histPtRapRecsFor1D = []
    histPtRecs = []
    histRapRecs = []
    for iCut, histRec in enumerate(histRecs):
        histRec.SetName(f'THnSparse_{sigRec}_{cuts[iCut]}') # Just to avoid memory leak

        histPtRapRecs.append(histRec.Projection(2, 1, f'{listRec}_{cuts[iCut]}_{sigRec}_PtRap'))
        histPtRapRecs[iCut].SetName(f'{listRec}_{cuts[iCut]}_{sigRec}_PtRap')

        histPtRapRecsFor1D.append(histRec.Projection(2, 1, f'{listRec}_{cuts[iCut]}_{sigRec}_PtRap_for_1D_proj'))
        histPtRapRecsFor1D[iCut].SetName(f'{listRec}_{cuts[iCut]}_{sigRec}_PtRap')

        histPtRapRecs[iCut].GetXaxis().SetRangeUser(0.0, 20.0)
        histPtRapRecs[iCut].GetYaxis().SetRangeUser(2.5, 4.0)

        histPtRapRecsFor1D[iCut].GetXaxis().SetRangeUser(minPtBinInt, maxPtBinInt)
        histPtRapRecsFor1D[iCut].GetYaxis().SetRangeUser(minRapBinInt, maxRapBinInt)

        #histPtRecs.append(histRec.Projection(1, f'{listRec}_{cuts[iCut]}_{sigRec}_Pt'))
        #histRapRecs.append(histRec.Projection(2, f'{listRec}_{cuts[iCut]}_{sigRec}_Rap'))

        histPtRecs.append(histPtRapRecsFor1D[iCut].ProjectionX(f'{listRec}_{cuts[iCut]}_{sigRec}_Pt'))
        histRapRecs.append(histPtRapRecsFor1D[iCut].ProjectionY(f'{listRec}_{cuts[iCut]}_{sigRec}_Rap'))

    histPtRebinRecs = []
    for iCut, histPtRec in enumerate(histPtRecs):
        histPtRebinRecs.append(histPtRec.Rebin(len(ptBinsRun3) - 1, f'{histPtRec.GetName()}_rebin', ptBinsRun3))

    histRapRebinRecs = []
    for iCut, histRapRec in enumerate(histRapRecs):
        histRapRebinRecs.append(histRapRec.Rebin(len(rapBinsRun3) - 1, f'{histRapRec.GetName()}_rebin', rapBinsRun3))


    canvasRec = ROOT.TCanvas("canvasRec", "", 800, 600)
    canvasRec.Divide(2, 2)

    canvasRec.cd(1)
    histPtRapRecsFor1D[0].Draw("COLZ")

    canvasRec.cd(2)
    ROOT.gPad.SetLogy(True)
    histRapRebinRecs[0].Draw("EP")

    canvasRec.cd(3)
    ROOT.gPad.SetLogy(True)
    histPtRebinRecs[0].Draw("EP")

    canvasRec.cd(4)
    histPtRapRecs[0].Draw("COLZ")

    canvasRec.Update()


    print("Compute the ratio")
    nBinsX = 25
    binEdgesX = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 12, 14, 16, 18, 20]
    nBinsY = 15
    binEdgesY = [2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4]

    ######################################
    # Rebinning generated for 2D Axe map #
    ######################################
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

    ##########################################
    # Rebinning reconstructed for 2D Axe map #
    ##########################################
    histPtRapMapRecsLowPt = []
    histPtRapMapRecsHighPt = []

    for iCut, histPtRapRec in enumerate(histPtRapRecs):
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
    for iCut, histPtRapRec in enumerate(histPtRapRecs):
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

    print(f'[WARNING] scale the generated for job efficiency (corr. factor = {corrFactor})')
    histPtRapMapGen.Scale(corrFactor)
    histPtRebinGen.Scale(corrFactor)
    histRapRebinGen.Scale(corrFactor)

    histPtRapAxes = []
    for iCut, histPtRapRec in enumerate(histPtRapRecs):
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

    histPtRebinAxes = []
    for iCut, histPtRebinRec in enumerate(histPtRebinRecs):
        histPtRebinAxes.append(histPtRebinRec.Clone(f'{histPtRebinRec.GetName()}_rebin_Axe'))
        histPtRebinAxes[iCut].Divide(histPtRebinAxes[iCut], histPtRebinGen, 1.0, 1.0, "B")

    histRapRebinAxes = []
    for iCut, histRapRebinRec in enumerate(histRapRebinRecs):
        histRapRebinAxes.append(histRapRebinRec.Clone(f'{histRapRebinRec.GetName()}_rebin_Axe'))
        histRapRebinAxes[iCut].Divide(histRapRebinAxes[iCut], histRapRebinGen, 1.0, 1.0, "B")

    canvasAxe = ROOT.TCanvas("canvasAxe", "", 800, 600)
    canvasAxe.Divide(2, 2)

    canvasAxe.cd(1)
    histPtRapAxes[0].Draw("COLZ")

    canvasAxe.cd(2)
    ROOT.gPad.SetLogy(False)
    histRapRebinAxes[0].Draw()
    
    canvasAxe.cd(3)
    ROOT.gPad.SetLogy(False)
    histPtRebinAxes[0].GetYaxis().SetRangeUser(0.001, 1.2)
    histPtRebinAxes[0].Draw()

    canvasAxe.Update()

    print("x_min x_max val stat syst ")
    for iPt in range(len(ptBinsRun3)-1):
        x_min = histPtRebinAxes[0].GetBinLowEdge(iPt+1)
        x_max = histPtRebinAxes[0].GetBinLowEdge(iPt+1) + histPtRebinAxes[0].GetBinWidth(iPt+1)
        val = histPtRebinAxes[0].GetBinContent(iPt+1)
        stat = histPtRebinAxes[0].GetBinError(iPt+1)
        syst = 0.0
        print("{:3.2f} {:3.2f} {:6.5f} {:6.5f} {:6.5f} ".format(x_min, x_max, val, stat, syst))

    input()

    fOutNameTxt = f'{pathOut}/{production}/{prefixOut}_{assocType}_{sigGen}.txt'
    print(f'Saving in {fOutNameTxt}...')
    fOutTxt = open(f'{fOutNameTxt}', 'w')
    fOutTxt.write("x_min x_max val stat syst \n")
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

    fOutNameRoot = f'{pathOut}/{production}/{prefixOut}_{assocType}_{sigGen}.root'
    print(f'Saving in {fOutNameRoot}...')
    fOutRoot = ROOT.TFile(f'{fOutNameRoot}', "RECREATE")
    histPtRapGen.Write("histPtRapGen")
    histPtRebinGen.Write("histPtRebinGen")
    histRapRebinGen.Write("histRapRebinGen")
    for iCut, cut in enumerate(cuts):
        histPtRapRecs[iCut].Write(f'histPtRapRecs_{cut}')
        histPtRebinRecs[iCut].Write(f'histPtRebinRecs_{cut}')
        histRapRebinRecs[iCut].Write(f'histRapRebinRecs_{cut}')
        histPtRapAxes[iCut].Write(f'histPtRapAxes_{cut}')
        histPtRebinAxes[iCut].Write(f'histPtRebinAxes_{cut}')
        histRapRebinAxes[iCut].Write(f'histRapRebinAxes_{cut}')
    fOutRoot.Close()
    if (fInNameGen == fInNameRec):
        fInGen.Close()
    else:
        fInGen.Close()
        fInRec.Close()
    
    exit()

if __name__ == '__main__':
    main()