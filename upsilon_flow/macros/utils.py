import sys
import argparse
import yaml
import ROOT
#sys.path.append('../utils')
#from plot_library import LoadStyle, SetGraStat, SetGraSyst, SetLegend

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--do_prefilter", help="Apply selections", action="store_true")
    args = parser.parse_args()
    

    with open(args.cfgFileName, 'r') as yml_cfg:
        inputCfg = yaml.load(yml_cfg, yaml.FullLoader)

    if args.do_prefilter:
        prefilter(inputCfg)

def prefilter(config):
    """
    function to apply pre-filters on D0 tree according to BDT output (see config_fit.yml)
    """
    fIn = ROOT.TFile(config["prefilter"]["data"], "READ")
    tree = fIn.Get(config["prefilter"]["treeIn"])
    rDataFrame = ROOT.RDataFrame(tree).Filter(config["prefilter"]["cuts"])
    fOutName = config["prefilter"]["data"].replace(".root", config["prefilter"]["suffix"] + ".root")
    rDataFrame.Snapshot(config["prefilter"]["treeOut"], fOutName)
    fIn.Close()

if __name__ == '__main__':
    main()