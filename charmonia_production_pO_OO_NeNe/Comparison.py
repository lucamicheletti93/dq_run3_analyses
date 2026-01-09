import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import argparse
import yaml
import ROOT
import numpy as np
import array
import pandas as pd


def read_txt(filename):
    """
    Legge file txt con colonne:
    x_min x_max val stat syst
    """
    x_min, x_max, val, stat, syst = np.loadtxt(filename, unpack=True, skiprows=1)

    x = 0.5 * (x_min + x_max)
    xerr = 0.5 * (x_max - x_min)

    return x, xerr, val, stat


# ==== file di input (con path) ====
file1 = "/home/rebecca/cernbox/dq_run3_analyses/charmonia_production_pO_OO_NeNe/data/MC/LHC25i4/axe_pt_dependence_Rap_2.5_4_Cent_0_90_std_association_Jpsi_100PercStat.txt"
file2 = "/home/rebecca/cernbox/dq_run3_analyses/charmonia_production_pO_OO_NeNe/data/MC/LHC25i4/axe_pt_dependence_Rap_2.5_4_NoCentrConsidered_std_association_Jpsi_100PercStat.txt" # axe_pt_dependence_Rap_2.5_4_BeforeCentrImplementationRECREATINGHistos_std_association_Jpsi_10PercStat, axe_pt_dependence_Rap_2.5_4_NoCentrConsidered_std_association_Jpsi_100PercStat

#file1 = "/home/rebecca/cernbox/dq_run3_analyses/charmonia_production_pO_OO_NeNe/data/MC/LHC25i4/axe_cent_dependence_Rap_2.5_4_Pt_0_20_std_association_Jpsi_100PercStat.txt"
#file2 = "/home/rebecca/cernbox/dq_run3_analyses/charmonia_production_pO_OO_NeNe/data/MC/LHC25i4/axe_integrated_Rap_2.5_4_Pt_0_20_NoCentrConsidered_std_association_Jpsi_100PercStat.txt" # axe_integrated_Rap_2.5_4_Pt_0_20_BeforeCentrImplementation_std_association_Jpsi_10PercStat, axe_integrated_Rap_2.5_4_Pt_0_20_NoCentrConsidered_std_association_Jpsi_100PercStat

# ==== lettura ====
x1, xerr1, y1, yerr1 = read_txt(file1)
x2, xerr2, y2, yerr2 = read_txt(file2)

# ==== label pulite ====
label1 = os.path.basename(file1)
label2 = os.path.basename(file2)
label2 = label2[:len(label2)//2]

# ==== plot ====
plt.figure(figsize=(7,5))

plt.errorbar(
    x1, y1, yerr=yerr1, xerr=xerr1,
    fmt='o', capsize=3, label=label1
)

plt.errorbar(
    x2, y2, yerr=yerr2, xerr=xerr2,
    fmt='s', capsize=3, label=label2
)

#plt.errorbar(
#    45, y2, yerr=yerr2, xerr=45,
#    fmt='s', capsize=3, label=label2
#)

plt.xlabel("pT (GeV/c)")
plt.ylabel("val")
plt.legend()
plt.ylim(0.2, 0.5)

plt.tight_layout()

# ==== salva e mostra ====
plt.savefig("/home/rebecca/cernbox/dq_run3_analyses/charmonia_production_pO_OO_NeNe/data/MC/LHC25i4/axe_pt_dependence_COMPARISON_Between_100PercNoCentrConsidered_and_100Perc0_90Centr.pdf") # axe_pt_dependence_COMPARISON_Between_10PercBeforeCentrImplementationRECREATINGHistos_and_100Perc0_90Centr, axe_pt_dependence_COMPARISON_Between_100PercNoCentrConsidered_and_100Perc0_90Centr
#plt.savefig("/home/rebecca/cernbox/dq_run3_analyses/charmonia_production_pO_OO_NeNe/data/MC/LHC25i4/axe_cent_dependence_COMPARISON_Between_100PercNoCentrConsidered_AND_YesCentrConsidered.pdf")  # axe_cent_dependence_COMPARISON_Between_10PercBeforeCentrImplementation_AND_100PercYesCentrConsidered, axe_cent_dependence_COMPARISON_Between_100PercNoCentrConsidered_AND_YesCentrConsidered