from telnetlib import Telnet
import matplotlib.pyplot as plt
import array as arr
import numpy as np
import os
import sys
import argparse
import yaml
import ROOT
from os import path

def download(inputCfg):
    if os.path.isfile(("{}/run_list.txt").format(inputCfg["input"]["output_dir_name"])):
        os.system(("rm {}/run_list.txt").format(inputCfg["input"]["output_dir_name"]))

    fOut = open(("{}/run_list.txt").format(inputCfg["input"]["output_dir_name"]), 'x')

    output_dir = inputCfg["input"]["output_dir_name"]
    if not os.path.isdir("%s" % (output_dir)):
        os.system("mkdir -p %s" % (output_dir))

    with open(inputCfg["input"]["alien_run_list"], 'r') as file:
        run_list = file.read().splitlines()

    with open(inputCfg["input"]["alien_input_path"], 'r') as file:
        alien_path = file.read().splitlines()

    print("----- Download and save files in %s -----" % (output_dir))
    for iRun, run in enumerate(run_list):
        file_type = inputCfg["input"]["file_type"]

        os.system("mkdir -p %s/%s" % (output_dir, run))

        #print("alien_cp alien://%s/%s file:%s/%s/." % (alien_path[iRun], file_type, output_dir, run))
        #os.system("alien_cp alien://%s/*/%s file:%s/%s/." % (alien_path[iRun], file_type, output_dir, run))
        os.system("alien_cp alien://%s/%s file:%s/%s/." % (alien_path[iRun], file_type, output_dir, run))
        fOut.write("{}\n".format(run))
    fOut.close()
    print("--- How to run the AO2D merger ---")
    print("o2-aod-merger --input input.txt --output AO2D_merged.root --max-size 1000000000 --skip-non-existing-files --skip-parent-files-list")

    exit()

    print("----- Download and save files in %s -----" % (inputCfg["input"]["output_dir_name"]))
    for iRun in range(0, len(inputCfg["input"]["run_list"])):

        file_type = inputCfg["input"]["file_type"]
        run = inputCfg["input"]["run_list"][iRun]
        alien_path = inputCfg["input"]["alien_input_path"][iRun]

        os.system("mkdir -p %s/%s" % (output_dir, run))
        
        #os.system("alien_cp alien://%s/*/%s file:%s/%s/." % (alien_path, file_type, output_dir, run))
        os.system("alien_cp alien://%s/%s file:%s/%s/." % (alien_path, file_type, output_dir, run))
        fOut.write("{}\n".format(run))
    fOut.close()

def fix(inputCfg):
    output_dir = inputCfg["input"]["output_dir_name"]

    with open(inputCfg["input"]["alien_run_list"], 'r') as file:
        run_list = file.read().splitlines()
    
    with open(inputCfg["input"]["alien_input_path"], 'r') as file:
        alien_path = file.read().splitlines()

    print("----- Download and save missing files in %s -----" % (output_dir))
    for iRun, run in enumerate(run_list):
        file_type = inputCfg["input"]["file_type"]

        if not os.path.isfile(("{}/{}/{}").format(output_dir, run, file_type)):
            print("Run ", run, " does not exist!")
            #os.system("alien_cp alien://%s/*/%s file:%s/%s/." % (alien_path[iRun], file_type, output_dir, run))


def merge(inputCfg):
    output_dir = inputCfg["input"]["output_dir_name"]

    with open(inputCfg["input"]["alien_run_list"], 'r') as file:
        run_list = file.read().splitlines()

    fInPath = inputCfg["input"]["file_path"]
    file_type = inputCfg["input"]["file_type"]
    os.system("mkdir -p {}/merged_files".format(fInPath))

    for run in run_list:
        #print("mkdir -p {}/merged_files/{}".format(fInPath, run))
        os.system(f'mkdir -p {fInPath}/merged_files/{run}')
        #print("hadd {}/merged_files/{}/AnalysisResults.root {}/{}/*/AnalysisResults.root".format(fInPath, run, fInPath, run))
        os.system(f'hadd {fInPath}/merged_files/{run}/{file_type} {fInPath}/{run}/*/{file_type}')

def merge_ttree(inputCfg):
    fInPath = inputCfg["input"]["file_path"]
    input_file = os.path.join(fInPath, "AO2D.root")
    output_file = os.path.join(fInPath, "AO2D_merged.root")

    root_file = ROOT.TFile.Open(input_file, "READ")
    if not root_file or root_file.IsZombie():
        print(f"Error: impossible to open file {input_file}")
        return

    df_dirs = [key.GetName() for key in root_file.GetListOfKeys() if key.GetName().startswith("DF")]
    if not df_dirs:
        print("None DF folder DF found in file AO2D.root.")
        return

    print(f"{len(df_dirs)} folders DF* found: {df_dirs}")

    tree_name = "O2rtdimuonall"

    root_file.cd(df_dirs[0])  
    if not root_file.Get(f"{df_dirs[0]}/{tree_name}"):
        print(f"Error: TTree '{tree_name}' doesn't exist in {df_dirs[0]}")
        return

    print(f"TTree '{tree_name}' found")

    temp_files = [] 

    for df in df_dirs:
        output_temp = os.path.join(fInPath, f"{df}_{tree_name}.root")
        temp_files.append(output_temp)

        output_root = ROOT.TFile(output_temp, "RECREATE")
        tree = root_file.Get(f"{df}/{tree_name}")
        if tree:
            output_root.cd()
            tree.CloneTree().Write() 
            output_root.Close()
        else:
            print(f"{df}/{tree_name} Tree not found!")

    hadd_command = f"hadd -f {output_file} {' '.join(temp_files)}"
    print(f"Executing: {hadd_command}")
    os.system(hadd_command)

    for temp in temp_files:
        os.remove(temp)

    print(f"Merging completed. File saved in: {output_file}")

### ### ###
def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--download", help="Download single files", action="store_true")
    parser.add_argument("--fix", help="Download missing files from download", action="store_true")
    parser.add_argument("--merge", help="Do the merging of the downloaded files", action="store_true")
    args = parser.parse_args()

    print('Loading task configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading task configuration: Done!')

    if args.download:
        download(inputCfg)
    if args.fix:
        fix(inputCfg)
    if args.merge:
        merge(inputCfg)

main()