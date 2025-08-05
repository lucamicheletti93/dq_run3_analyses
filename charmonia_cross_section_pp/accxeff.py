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
    prefixIn = config["inputs"]["prefix"]
    assocType = config["inputs"]["assocType"]
    fInName = f'{pathIn}/{production}/{prefixIn}_{assocType}.root'

    pathOut = config["outputs"]["path"]
    prefixOut = config["outputs"]["prefix"]

    print(f'Opening {fInName}...')
    fIn = ROOT.TFile(f'{fInName}', "READ")

    # Generated distributions
    hlistGen = fIn.Get(config["inputs"]["hlistGen"])
    strListGen = config["inputs"]["listGen"]
    strHistGen = config["inputs"]["histGen"]
    sigGen = config["inputs"]["sigGen"]

    minPtBinInt = config["inputs"]["minPtBinInt"]
    maxPtBinInt = config["inputs"]["maxPtBinInt"]
    minRapBinInt = config["inputs"]["minRapBinInt"]
    maxRapBinInt = config["inputs"]["maxRapBinInt"]

    histPtRapGen = (hlistGen.FindObject(f'{strListGen}_{sigGen}')).FindObject(strHistGen)

    histPtRapGen.GetXaxis().SetRangeUser(minPtBinInt, maxPtBinInt)
    histPtRapGen.GetYaxis().SetRangeUser(minRapBinInt, maxRapBinInt)
    histPtGen = histPtRapGen.ProjectionX(f'{sigGen}_Pt')
    histRapGen = histPtRapGen.ProjectionY(f'{sigGen}_Rap')

    ptBinsRun3 = array.array('d', config["inputs"]["ptBins"])
    rapBinsRun3 = array.array('d', config["inputs"]["rapBins"])

    histPtRebinGen = histPtGen.Rebin(len(ptBinsRun3) - 1, f'{histPtGen.GetName()}_rebin', ptBinsRun3)
    histRapRebinGen = histRapGen.Rebin(len(rapBinsRun3) - 1, f'{histRapGen.GetName()}_rebin', rapBinsRun3)
        
    canvasGen = ROOT.TCanvas("canvasGen", "", 800, 600)
    canvasGen.Divide(2, 2)

    canvasGen.cd(1)
    histPtRapGen.Draw("COLZ")

    canvasGen.cd(2)
    ROOT.gPad.SetLogy(True)
    histRapRebinGen.Draw()

    canvasGen.cd(3)
    ROOT.gPad.SetLogy(True)
    histPtRebinGen.Draw()

    canvasGen.Update()

    # Reconstructed distributions
    hlistRec = fIn.Get(config["inputs"]["hlistRec"])
    listRec = config["inputs"]["listRec"]
    histRec = config["inputs"]["histRec"]
    cuts = config["inputs"]["cuts"]
    sigRec = config["inputs"]["sigRec"]

    histRecs = []
    for iCut, cut in enumerate(cuts):
        histRecs.append((hlistRec.FindObject(f'{listRec}_{cut}_{sigRec}')).FindObject(f'{histRec}'))
        histRecs[iCut].GetAxis(1).SetRangeUser(minPtBinInt, maxPtBinInt)
        histRecs[iCut].GetAxis(2).SetRangeUser(minRapBinInt, maxRapBinInt)

    histPtRapRecs = []
    histPtRecs = []
    histRapRecs = []
    for iCut, histRec in enumerate(histRecs):
        histRec.SetName(f'THnSparse_{sigRec}_{cuts[iCut]}') # Just to avoid memory leak
        histPtRapRecs.append(histRec.Projection(2, 1, f'{listRec}_{cuts[iCut]}_{sigRec}_PtRap'))
        histPtRapRecs[iCut].GetXaxis().SetRangeUser(minPtBinInt, maxPtBinInt)
        histPtRapRecs[iCut].GetYaxis().SetRangeUser(minRapBinInt, maxRapBinInt)

        histPtRecs.append(histRec.Projection(1, f'{listRec}_{cuts[iCut]}_{sigRec}_Pt'))
        histRapRecs.append(histRec.Projection(2, f'{listRec}_{cuts[iCut]}_{sigRec}_Rap'))


    histPtRebinRecs = []
    for iCut, histPtRec in enumerate(histPtRecs):
        histPtRebinRecs.append(histPtRec.Rebin(len(ptBinsRun3) - 1, f'{histPtRec.GetName()}_rebin', ptBinsRun3))

    histRapRebinRecs = []
    for iCut, histRapRec in enumerate(histRapRecs):
        histRapRebinRecs.append(histRapRec.Rebin(len(rapBinsRun3) - 1, f'{histRapRec.GetName()}_rebin', rapBinsRun3))


    canvasRec = ROOT.TCanvas("canvasRec", "", 800, 600)
    canvasRec.Divide(2, 2)

    canvasRec.cd(1)
    histPtRapRecs[0].Draw("COLZ")

    canvasRec.cd(2)
    ROOT.gPad.SetLogy(True)
    histRapRebinRecs[0].Draw("COLZ")
    

    canvasRec.cd(3)
    ROOT.gPad.SetLogy(True)
    histPtRebinRecs[0].Draw("COLZ")

    canvasRec.Update()


    # Compute the ratio
    histPtRapAxes = []
    for iCut, histPtRapRec in enumerate(histPtRapRecs):
        histPtRapGen.RebinX(5)
        histPtRapRec.RebinX(5)
        histPtRapGen.RebinY(5)
        histPtRapRec.RebinY(2)
        #print(histPtRapRec.GetNbinsX(), " ", histPtRapRec.GetNbinsY())
        #print(histPtRapGen.GetNbinsX(), " ", histPtRapGen.GetNbinsY())

        histPtRapAxes.append(histPtRapRec.Clone(f'{histPtRapRec.GetName()}_rebin_Axe'))
        #histPtRapAxes[iCut].Divide(histPtRapGens[iCut], histPtRapRec, 1.0, 1.0, "B")
        histPtRapAxes[iCut].Divide(histPtRapGen)
        histPtRapAxes[iCut].GetXaxis().SetRangeUser(minPtBinInt, maxPtBinInt)
        histPtRapAxes[iCut].GetYaxis().SetRangeUser(minRapBinInt, maxRapBinInt)

    histPtRebinAxes = []
    for iCut, histPtRebinRec in enumerate(histPtRebinRecs):
        histPtRebinAxes.append(histPtRebinRec.Clone(f'{histPtRebinRec.GetName()}_rebin_Axe'))
        #histPtRebinAxes[iCut].Divide(histPtRebinGens[iCut], histPtRebinRec, 1.0, 1.0, "B")
        histPtRebinAxes[iCut].Divide(histPtRebinGen)

    histRapRebinAxes = []
    for iCut, histRapRebinRec in enumerate(histRapRebinRecs):
        histRapRebinAxes.append(histRapRebinRec.Clone(f'{histRapRebinRec.GetName()}_rebin_Axe'))
        #histRapRebinAxes[iCut].Divide(histRapRebinRec, histRapRebinGens[iCut], 1.0, 1.0, "B")
        histRapRebinAxes[iCut].Divide(histRapRebinGen)

    canvasAxe = ROOT.TCanvas("canvasAxe", "", 800, 600)
    canvasAxe.Divide(2, 2)

    canvasAxe.cd(1)
    ROOT.gPad.SetLogz(True)
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
    for iPt in range(len(ptBinsRun3)):
        x_min = histPtRebinAxes[0].GetBinLowEdge(iPt+1)
        x_max = histPtRebinAxes[0].GetBinLowEdge(iPt+1) + histPtRebinAxes[0].GetBinWidth(iPt+1)
        val = histPtRebinAxes[0].GetBinContent(iPt+1)
        stat = histPtRebinAxes[0].GetBinError(iPt+1)
        syst = 0.0
        print("{:3.2f} {:3.2f} {:6.5f} {:6.5f} {:6.5f} ".format(x_min, x_max, val, stat, syst))

    input()

    fOutName = f'{pathOut}/{prefixOut}_{assocType}_{sigGen}.root'

    print(f'Saving in {fOutName}...')
    fOut = ROOT.TFile(f'{fOutName}', "RECREATE")
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
    fOut.Close()
    fIn.Close()
    
    exit()

if __name__ == '__main__':
    main()