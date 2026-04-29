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

values = np.array([0.641371, 1.18315, 0.93675, 0.54179, 0.262163, 0.138526, 0.0505504])
statErr = np.array([0.0238, 0.0311, 0.0273, 0.0193, 0.0126, 0.0076, 0.0042])
systErr = np.array([0.0357188861417878, 0.0575403176217337, 0.0406919444422357, 0.0236617786787046, 0.0148124985927425, 0.0026628, 0.00255747469391195, 0.00108786176051923])

ptEdges = np.array([0, 1, 2, 3, 4, 5, 6, 8])

nBins = len(values)

histStat = ROOT.TH1D("histStatJpsiXsecInterp", ";p_{T};Cross section", nBins, ptEdges.astype('double'))
histSyst = ROOT.TH1D("histSystJpsiXsecInterp", ";p_{T};Cross section", nBins, ptEdges.astype('double'))

for i in range(nBins):

    histStat.SetBinContent(i+1, values[i])
    histStat.SetBinError(i+1, statErr[i])

    histSyst.SetBinContent(i+1, values[i])
    histSyst.SetBinError(i+1, systErr[i])

fOut = ROOT.TFile("/home/rebecca/cernbox/dq_run3_analyses/charmonia_production_pO_OO_NeNe/data/ppRef/WITH_5.02_UNCERT_new_interp_jpsi_xsec_pp5dot36TeV.root", "RECREATE")

histStat.Write()
histSyst.Write()

fOut.Close()
