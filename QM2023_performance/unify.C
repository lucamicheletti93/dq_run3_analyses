//This macro is used to unify the partial outputs coming from jPsiStudies.C in a single output file to be fitted
//using the DQ fit library

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMinuit.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMath.h"
#include <fstream>
#include <TStyle.h>
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TTree.h"
#include "TBranch.h"
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <string>
#include "THStack.h"
#include "TLegend.h"

#include <tuple>

using namespace std;

void unify() {

    //Open the output files from jPsiStudies.C
    string out1 = "out_LHC22o.root";
    //string out2 = "out_LHC22others.root";
    string out2 = "out_LHC22others_noLHC22m.root";
    TFile *fOut1 = new TFile(out1.c_str(),"READ");
    TFile *fOut2 = new TFile(out2.c_str(),"READ");

    fOut1->cd(); //Get invariant mass spectra from LHC22o

    TH1F *hChi2_45_SEPM_1_pt_0_2 = (TH1F*)fOut1->Get("massSEPMetaSelectedChi2_45_pt_0_2");
    TH1F *hChi2_45_SEPM_1_pt_2_4 = (TH1F*)fOut1->Get("massSEPMetaSelectedChi2_45_pt_2_4");
    TH1F *hChi2_45_SEPM_1_pt_4_6 = (TH1F*)fOut1->Get("massSEPMetaSelectedChi2_45_pt_4_6");
    TH1F *hChi2_45_SEPM_1_pt_20 = (TH1F*)fOut1->Get("massSEPMetaSelectedChi2_45_pt_20");

    fOut2->cd(); //Get invariant mass spectra from other periods (LHC22 t + LHC22r)

    TH1F *hChi2_45_SEPM_2_pt_0_2 = (TH1F*)fOut2->Get("massSEPMetaSelectedChi2_45_pt_0_2");
    TH1F *hChi2_45_SEPM_2_pt_2_4 = (TH1F*)fOut2->Get("massSEPMetaSelectedChi2_45_pt_2_4");
    TH1F *hChi2_45_SEPM_2_pt_4_6 = (TH1F*)fOut2->Get("massSEPMetaSelectedChi2_45_pt_4_6");
    TH1F *hChi2_45_SEPM_2_pt_20 = (TH1F*)fOut2->Get("massSEPMetaSelectedChi2_45_pt_20");

    //Create histograms to join the data from all periods
    TH1F *hChi2_45_SEPM_1_pt_0_2_tot = new TH1F("TH1F *hChi2_45_SEPM_1_pt_0_2_tot","TH1F *hChi2_45_SEPM_1_pt_0_2_tot",7000,0.,140.);
    TH1F *hChi2_45_SEPM_1_pt_2_4_tot = new TH1F("TH1F *hChi2_45_SEPM_1_pt_2_4_tot","TH1F *hChi2_45_SEPM_1_pt_2_4_tot",7000,0.,140.);
    TH1F *hChi2_45_SEPM_1_pt_4_6_tot = new TH1F("TH1F *hChi2_45_SEPM_1_pt_4_6_tot","TH1F *hChi2_45_SEPM_1_pt_4_6_tot",7000,0.,140.);
    TH1F *hChi2_45_SEPM_1_pt_20_tot = new TH1F("TH1F *hChi2_45_SEPM_1_pt_20_tot","TH1F *hChi2_45_SEPM_1_pt_20_tot",7000,0.,140.);

    //Add data of the two periods
    hChi2_45_SEPM_1_pt_0_2->Add(hChi2_45_SEPM_2_pt_0_2);
    hChi2_45_SEPM_1_pt_2_4->Add(hChi2_45_SEPM_2_pt_2_4);
    hChi2_45_SEPM_1_pt_4_6->Add(hChi2_45_SEPM_2_pt_4_6);
    hChi2_45_SEPM_1_pt_20->Add(hChi2_45_SEPM_2_pt_20);

    //Gloabl output file, without LHC22m period
    TFile *fout = new TFile("outGlobal_no_LHC22m.root","RECREATE");
    fout->cd();

    hChi2_45_SEPM_1_pt_0_2->Write("mass_pt_0_2"); //0 < pt < 2 GeV/c +/- 
    hChi2_45_SEPM_1_pt_2_4->Write("mass_pt_2_4"); //2 < pt < 4 GeV/c +/- 
    hChi2_45_SEPM_1_pt_4_6->Write("mass_pt_4_6"); //4 < pt < 6 GeV/c +/- 
    hChi2_45_SEPM_1_pt_20->Write("mass_pt_20"); //pt integrated

    fout->Close();
}	