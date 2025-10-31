#from re import TEMPLATE
import matplotlib.pyplot as plt
import array as arr
import yaml
import numpy as np
from array import array
import os
import sys
import math
import re
from statistics import mean
import argparse
from os import path

def ExtractFromYaml(file_path):
    #file_path = "/Users/lucamicheletti/GITHUB/dq_run3_analyses/psi2S_jpsi_ratio/results/HEPData-77781-v1-yaml/Table5.yaml"

    with open(file_path, "r") as f:
        hep_data = yaml.safe_load(f)

        independent_var = hep_data["independent_variables"][0]["values"]
        #centers = [(entry["low"] + entry["high"]) / 2 for entry in independent_var]
        #widths = [(entry["high"] - entry["low"]) / 2 for entry in independent_var]
        centers = []
        widths = []
        for entry in independent_var:
            if "low" in entry and "high" in entry:
                centers.append((entry["low"] + entry["high"]) / 2)
                widths.append((entry["high"] - entry["low"]) / 2)
            elif "value" in entry:
                centers.append(entry["value"])
                widths.append(0.0)
            else:
                raise ValueError(f"Independent variable format not recognized: {entry}")

        # Estrarre le variabili dipendenti ($\mathrm{DSIG}/\mathrm{DPT}/\mathrm{DYRAP}$)
        dependent_var = hep_data["dependent_variables"][0]["values"]
        y_values = [entry["value"] for entry in dependent_var]

        # Calcolo degli errori
        stat_errors = [entry["errors"][0]["symerror"] for entry in dependent_var]
        syst_errors = [
            np.sqrt(
                sum(
                    err["symerror"]**2
                    for err in entry["errors"]
                    if "sys" in err["label"]
                )
            )
            for entry in dependent_var
        ]

        # Errore totale (statistico + sistematico)
        total_errors = np.sqrt(np.array(stat_errors)**2 + np.array(syst_errors)**2)

        return centers, widths, y_values, stat_errors, syst_errors