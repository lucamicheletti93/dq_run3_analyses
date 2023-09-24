//This macro produces the starting file (.root) for the fit to the dimuon invariant mass spectrum in the J/Psi
//regione using run3 data and matched tracks between MCH and MFT

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
#include "TDirectory.h"
#include "TKey.h"
#include "TROOT.h"

#include <tuple>

void jPsiStudies() {

    //Path where data will be saved on cernbox
    string cernBoxPath = "/eos/user/l/lquaglia/run3Analysis/"; 

    //Name of input files, uncomment only the line concerning the data you are interested in
    //LHC22o -> saved on Luca Micheletti's cernbox folder, reduced tree containing only the data we care about
    //LHC22r -> divided into two separete files: AO2D_partial1 and AO2D_partial2 -> saved on Luca Quaglia's cernbox folder
    //LHC22m/LHC22t -> saved on Luca Quaglia's cernbox folder
    //string inputFiles[4] = {cernBoxPath+"AO2D_LHC22m.root",cernBoxPath+"AO2D_LHC22t.root",cernBoxPath+"AO2D_partial1.root",cernBoxPath+"AO2D_partial2.root"};
    //string inputFiles[3] = {cernBoxPath+"AO2D_LHC22t.root",cernBoxPath+"AO2D_partial1.root",cernBoxPath+"AO2D_partial2.root"}; 
    string inputFiles[1] = {"/eos/user/l/lmichele/QM2023_DQ_performance/LHC22o/528379/001/AO2D_reduced.root"}; 

    //Output file, according to what you want to analyze
    //out_LHC22others -> LHC22m/t/r
    //out_LHC22others_noLHC22m -> LHC22t/r
    //out_LHC22o -> LHC22o
    //TFile *fOut = new TFile("out_LHC22others.root","RECREATE");
    //TFile *fOut = new TFile("out_LHC22others_noLHC22m.root","RECREATE");
    TFile *fOut = new TFile("out_LHC22o.root","RECREATE");

    //Eta of selected tracks, same binning as O2 physics
    TH1F *etaSelected_pt_0_2_SEPM = new TH1F("etaSelected_pt_0_2_SEPM","etaSelected_pt_0_2_SEPM",500,-5.,5.); //0 < pt < 2 GeV/c dimuon pt +/-
    TH1F *etaSelected_pt_0_2_SEPP = new TH1F("etaSelected_pt_0_2_SEPP","etaSelected_pt_0_2_SEPP",500,-5.,5.); //0 < pt < 2 GeV/c dimuon pt +/+
    TH1F *etaSelected_pt_0_2_SEMM = new TH1F("etaSelected_pt_0_2_SEMM","etaSelected_pt_0_2_SEMM",500,-5.,5.); //0 < pt < 2 GeV/c dimuon pt -/-
    
    TH1F *etaSelected_pt_2_4_SEPM = new TH1F("etaSelected_pt_2_4_SEPM","etaSelected_pt_2_4_SEPM",500,-5.,5.); //2 < pt < 4 GeV/c dimuon pt +/-
    TH1F *etaSelected_pt_2_4_SEPP = new TH1F("etaSelected_pt_2_4_SEPP","etaSelected_pt_2_4_SEPP",500,-5.,5.); //2 < pt < 4 GeV/c dimuon pt +/+
    TH1F *etaSelected_pt_2_4_SEMM = new TH1F("etaSelected_pt_2_4_SEMM","etaSelected_pt_2_4_SEMM",500,-5.,5.); //2 < pt < 4 GeV/c dimuon pt -/-

    TH1F *etaSelected_pt_4_6_SEPM = new TH1F("etaSelected_pt_4_6_SEPM","etaSelected_pt_4_6_SEPM",500,-5.,5.); //4 < pt < 6 GeV/c dimuon pt +/-
    TH1F *etaSelected_pt_4_6_SEPP = new TH1F("etaSelected_pt_4_6_SEPP","etaSelected_pt_4_6_SEPP",500,-5.,5.); //4 < pt < 6 GeV/c dimuon pt +/+
    TH1F *etaSelected_pt_4_6_SEMM = new TH1F("etaSelected_pt_4_6_SEMM","etaSelected_pt_4_6_SEMM",500,-5.,5.); //4 < pt < 6 GeV/c dimuon pt -/-

    TH1F *etaSelected_pt_20_SEPM = new TH1F("etaSelected_pt_20_SEPM","etaSelected_pt_20_SEPM",500,-5.,5.); //pt integrated +/-
    TH1F *etaSelected_pt_20_SEPP = new TH1F("etaSelected_pt_20_SEPP","etaSelected_pt_20_SEPP",500,-5.,5.); //pt integrated +/+
    TH1F *etaSelected_pt_20_SEMM = new TH1F("etaSelected_pt_20_SEMM","etaSelected_pt_20_SEMM",500,-5.,5.); //pt integrated -/-


    //Chi2 of matched MCH/MID tracks
    TH1F *chi2MatchMFTMCHetaSelected_pt_0_2_SEPM = new TH1F("chi2MatchMFTMCHetaSelected_pt_0_2_SEPM","chi2MatchMFTMCHetaSelected_pt_0_2_SEPM",1000,0.,1e6); //0 < pt < 2 GeV/c dimuon pt +/-
    TH1F *chi2MatchMFTMCHetaSelected_pt_0_2_SEPP = new TH1F("chi2MatchMFTMCHetaSelected_pt_0_2_SEPP","chi2MatchMFTMCHetaSelected_pt_0_2_SEPP",1000,0.,1e6); //0 < pt < 2 GeV/c dimuon pt +/+
    TH1F *chi2MatchMFTMCHetaSelected_pt_0_2_SEMM = new TH1F("chi2MatchMFTMCHetaSelected_pt_0_2_SEMM","chi2MatchMFTMCHetaSelected_pt_0_2_SEMM",1000,0.,1e6); //0 < pt < 2 GeV/c dimuon pt -/-

    TH1F *chi2MatchMFTMCHetaSelected_pt_2_4_SEPM = new TH1F("chi2MatchMFTMCHetaSelected_pt_2_4_SEPM","chi2MatchMFTMCHetaSelected_pt_2_4_SEPM",1000,0.,1e6); //2 < pt < 4 GeV/c dimuon pt +/-
    TH1F *chi2MatchMFTMCHetaSelected_pt_2_4_SEPP = new TH1F("chi2MatchMFTMCHetaSelected_pt_2_4_SEPP","chi2MatchMFTMCHetaSelected_pt_2_4_SEPP",1000,0.,1e6); //2 < pt < 4 GeV/c dimuon pt +/+
    TH1F *chi2MatchMFTMCHetaSelected_pt_2_4_SEMM = new TH1F("chi2MatchMFTMCHetaSelected_pt_2_4_SEMM","chi2MatchMFTMCHetaSelected_pt_2_4_SEMM",1000,0.,1e6); //2 < pt < 4 GeV/c dimuon pt -/-

    TH1F *chi2MatchMFTMCHetaSelected_pt_4_6_SEPM = new TH1F("chi2MatchMFTMCHetaSelected_pt_4_6_SEPM","chi2MatchMFTMCHetaSelected_pt_4_6_SEPM",1000,0.,1e6); //4 < pt < 6 GeV/c dimuon pt
    TH1F *chi2MatchMFTMCHetaSelected_pt_4_6_SEPP = new TH1F("chi2MatchMFTMCHetaSelected_pt_4_6_SEPP","chi2MatchMFTMCHetaSelected_pt_4_6_SEPP",1000,0.,1e6); //4 < pt < 6 GeV/c dimuon pt
    TH1F *chi2MatchMFTMCHetaSelected_pt_4_6_SEMM = new TH1F("chi2MatchMFTMCHetaSelected_pt_4_6_SEMM","chi2MatchMFTMCHetaSelected_pt_4_6_SEMM",1000,0.,1e6); //4 < pt < 6 GeV/c dimuon pt

    TH1F *chi2MatchMFTMCHetaSelected_pt_20_SEPM = new TH1F("chi2MatchMFTMCHetaSelected_pt_20_SEPM","chi2MatchMFTMCHetaSelected_pt_20_SEPM",1000,0.,1e6); //pt integrated 
    TH1F *chi2MatchMFTMCHetaSelected_pt_20_SEPP = new TH1F("chi2MatchMFTMCHetaSelected_pt_20_SEPP","chi2MatchMFTMCHetaSelected_pt_20_SEPP",1000,0.,1e6); //pt integrated 
    TH1F *chi2MatchMFTMCHetaSelected_pt_20_SEMM = new TH1F("chi2MatchMFTMCHetaSelected_pt_20_SEMM","chi2MatchMFTMCHetaSelected_pt_20_SEMM",1000,0.,1e6); //pt integrated 
    
    //Invariant mass spectrum, no cuts on chi2
    TH1F *hMassSEPM = new TH1F("hMassSEPM","hMassSEPM",7000,0.,140.); // +/- pairs
    TH1F *hMassSEPP = new TH1F("hMassSEPP","hMassSEPP",7000,0.,140.); // +/+ pairs
    TH1F *hMassSEMM = new TH1F("hMassSEMM","hMassSEMM",7000,0.,140.); // -/- pairs

    //Invariant mass spectrum, chi2 < 45, 0 < pt < 2 GeV/c +/- pairs
    TH1F *hMassSEPMchi2_45_pt_0_2 = new TH1F("hMassSEPMchi2_45_pt_0_2","hMassSEPMchi2_45_pt_0_2",7000,0.,140.); 
    hMassSEPMchi2_45_pt_0_2->SetLineColor(kBlack);
    hMassSEPMchi2_45_pt_0_2->SetLineWidth(2);

    //Invariant mass spectrum, chi2 < 45, 2 < pt < 4 GeV/c +/- pairs
    TH1F *hMassSEPMchi2_45_pt_2_4 = new TH1F("hMassSEPMchi2_45_pt_2_4","hMassSEPMchi2_45_pt_2_4",7000,0.,140.); 
    hMassSEPMchi2_45_pt_2_4->SetLineColor(kRed);
    hMassSEPMchi2_45_pt_2_4->SetLineWidth(2);

    //Invariant mass spectrum, chi2 < 45, 4 < pt < 6 GeV/c +/- pairs
    TH1F *hMassSEPMchi2_45_pt_4_6 = new TH1F("hMassSEPMchi2_45_pt_4_6","hMassSEPMchi2_45_pt_4_6",7000,0.,140.); 
    hMassSEPMchi2_45_pt_4_6->SetLineColor(kSpring);
    hMassSEPMchi2_45_pt_4_6->SetLineWidth(2);

    //Invariant mass spectrum, chi2 < 45, pt integrated +/- pairs
    TH1F *hMassSEPMchi2_45_pt_20 = new TH1F("hMassSEPMchi2_45_pt_20","hMassSEPMchi2_45_pt_20",7000,0.,140.); 
    hMassSEPMchi2_45_pt_20->SetLineColor(kBlue);
    hMassSEPMchi2_45_pt_20->SetLineWidth(2);

    //Invariant mass spectrum, chi2 < 45, 0 < pt < 2 GeV/c +/+ pairs
    TH1F *hMassSEPPchi2_45_pt_0_2 = new TH1F("hMassSEPPchi2_45","hMassSEPPchi2_45",7000,0.,140.); 
    hMassSEPPchi2_45_pt_0_2->SetLineColor(kBlack);
    hMassSEPPchi2_45_pt_0_2->SetLineWidth(2);

    //Invariant mass spectrum, chi2 < 45, 2 < pt < 4 GeV/c +/+ pairs
    TH1F *hMassSEPPchi2_45_pt_2_4 = new TH1F("hMassSEPPchi2_45_pt_2_4","hMassSEPPchi2_45_pt_2_4",7000,0.,140.); 
    hMassSEPPchi2_45_pt_2_4->SetLineColor(kRed);
    hMassSEPPchi2_45_pt_2_4->SetLineWidth(2);

    //Invariant mass spectrum, chi2 < 45, 4 < pt < 6 GeV/c +/+ pairs
    TH1F *hMassSEPPchi2_45_pt_4_6 = new TH1F("hMassSEPPchi2_45_pt_4_6","hMassSEPPchi2_45_pt_4_6",7000,0.,140.); 
    hMassSEPPchi2_45_pt_4_6->SetLineColor(kSpring);
    hMassSEPPchi2_45_pt_4_6->SetLineWidth(2);

    //Invariant mass spectrum, chi2 < 45, pt integrated +/+ pairs
    TH1F *hMassSEPPchi2_45_pt_20 = new TH1F("hMassSEPPchi2_45_pt_20","hMassSEPPchi2_45_pt_20",7000,0.,140.); 
    hMassSEPPchi2_45_pt_20->SetLineColor(kBlue);
    hMassSEPPchi2_45_pt_20->SetLineWidth(2);

    //Invariant mass spectrum, chi2 < 45, 0 < pt < 2 GeV/c -/- pairs
    TH1F *hMassSEMMchi2_45_pt_0_2 = new TH1F("hMassSEMMchi2_45","hMassSEMMchi2_45",7000,0.,140.); 
    hMassSEMMchi2_45_pt_0_2->SetLineColor(kBlack);
    hMassSEMMchi2_45_pt_0_2->SetLineWidth(2);

    //Invariant mass spectrum, chi2 < 45, 2 < pt < 4 GeV/c -/- pairs
    TH1F *hMassSEMMchi2_45_pt_2_4 = new TH1F("hMassSEMMchi2_45_pt_2_4","hMassSEMMchi2_45_pt_2_4",7000,0.,140.); 
    hMassSEMMchi2_45_pt_2_4->SetLineColor(kRed);
    hMassSEMMchi2_45_pt_2_4->SetLineWidth(2);

    //Invariant mass spectrum, chi2 < 45, 4 < pt < 6 GeV/c -/- pairs
    TH1F *hMassSEMMchi2_45_pt_4_6 = new TH1F("hMassSEMMchi2_45_pt_4_6","hMassSEMMchi2_45_pt_4_6",7000,0.,140.); 
    hMassSEMMchi2_45_pt_4_6->SetLineColor(kSpring);
    hMassSEMMchi2_45_pt_4_6->SetLineWidth(2);

    //Invariant mass spectrum, chi2 < 45, pt integrated -/- pairs
    TH1F *hMassSEMMchi2_45_pt_20 = new TH1F("hMassSEMMchi2_45_pt_20","hMassSEMMchi2_45_pt_20",7000,0.,140.); 
    hMassSEMMchi2_45_pt_20->SetLineColor(kBlue);
    hMassSEMMchi2_45_pt_20->SetLineWidth(2);

    //chi2 values legend, entries taken from the +/- spectra just to have correct colors
    TLegend *lChi2 = new TLegend(0.57,0.64,0.87,0.84,"","rNDC");
    lChi2->SetBorderSize(0);	
	lChi2->SetFillStyle(0);
    lChi2->AddEntry(hMassSEPMchi2_45_pt_0_2,"#Chi^{2} <= 45 0 < p_{t} < 2 GeV/c","l");
    lChi2->AddEntry(hMassSEPMchi2_45_pt_2_4,"#Chi^{2} <= 45 2 < p_{t} < 4 GeV/c","l"); 
    lChi2->AddEntry(hMassSEPMchi2_45_pt_4_6,"#Chi^{2} <= 45 4 < p_{t} < 6 GeV/c","l");
    lChi2->AddEntry(hMassSEPMchi2_45_pt_20,"#Chi^{2} <= 45 p_{t} < 20 GeV/c","l");

    //THStacks for invariant mass
    THStack *stackMassChi2 = new THStack(); // +/-
    THStack *stackMassChi2SEPP = new THStack(); // +/+
    THStack *stackMassChi2SEMM = new THStack(); // -/-

    int totEntries = 0; //total number of entries in a given period
    int selectedEntries = 0; //number of entries that pass the eta cut

    //Variables to read TTree entries
    float eta1 = 0., eta2 = 0.; //eta of the muon track
    float chi2MatchMCHMFT1 = 0.,chi2MatchMCHMFT2 = 0.; //chi square of track matching between MFT-MCH
    int sign1 = 0., sign2 = 0.; //charge of dimuons
    float mass = 0.; //invariant mass of the dimuon pair
    float pt1 = 0., pt2 = 0.; //pt of single muon
    float pt = 0.; //pt of dimuon pair

    //for (int j = 0; j < 3; j++) {
    for (int j = 0; j < 1; j++) {    
        TFile *f = new TFile(inputFiles[j].c_str(),"READ"); //Open test file
        cout << "File name: " << inputFiles[j] << endl;
        TIter keyList(f->GetListOfKeys());
        TKey *key;

        while ((key = (TKey*)keyList())) {
            TClass *cl = gROOT->GetClass(key->GetClassName());
            if (!cl->InheritsFrom("TDirectoryFile")) 
                continue;
            TDirectoryFile *d = (TDirectoryFile*)key->ReadObj();
            d->cd();
            TTree *t = (TTree*)d->Get("O2rtdimuonall");
            int nEntries = t->GetEntries();
            totEntries += nEntries; //Sum up all the entries
            cout << "Folder name:" << key->GetTitle() << " entires in the tree: " << nEntries << endl;

            //Set branch address to get eta of the tracks
            t->SetBranchAddress("fEta1",&eta1);
            t->SetBranchAddress("fEta2",&eta2);
            t->SetBranchAddress("fChi2MatchMCHMFT1",&chi2MatchMCHMFT1);
            t->SetBranchAddress("fChi2MatchMCHMFT2",&chi2MatchMCHMFT2);
            t->SetBranchAddress("fSign1",&sign1);
            t->SetBranchAddress("fSign2",&sign2);
            t->SetBranchAddress("fMass",&mass);
            t->SetBranchAddress("fPt1",&pt1);
            t->SetBranchAddress("fPt2",&pt2);
            t->SetBranchAddress("fPt",&pt);

            for (int i = 0; i < nEntries; i++) { //Loop through all the tree entires
                t->GetEntry(i); //Get i-th entry

                //MFT eta range, unlike sign dimuon pair, pt of single muon > 0, 0 < pt < 2 GeV/c (dimuon) 
                if ((eta1 >= -3.6 && eta1 <= -2.5) && (eta2 >= -3.6 && eta2 <= -2.5) && ((sign1 == 1 && sign2 == -1) || (sign1 == -1 && sign2 == 1)) && (pt > 0 && pt < 2)) { //Muon quality cuts

                    etaSelected_pt_0_2_SEPM->Fill(eta1);
                    etaSelected_pt_0_2_SEPM->Fill(eta2);

                    hMassSEPM->Fill(mass);

                    //chi2 < 45
                    if (chi2MatchMCHMFT1 <= 45 && chi2MatchMCHMFT2 <= 45) {
                        chi2MatchMFTMCHetaSelected_pt_0_2_SEPM->Fill(chi2MatchMCHMFT1);
                        chi2MatchMFTMCHetaSelected_pt_0_2_SEPM->Fill(chi2MatchMCHMFT2);
                        hMassSEPMchi2_45_pt_0_2->Fill(mass);
                    }
                }

                //MFT eta range, +/+, pt of single muon > 0, 0 < pt < 2 GeV/c (dimuon) 
                else if ((eta1 >= -3.6 && eta1 <= -2.5) && (eta2 >= -3.6 && eta2 <= -2.5) && (sign1 == 1 && sign2 == 1) && (pt > 0 && pt < 2)) { //Muon quality cuts

                    etaSelected_pt_0_2_SEPP->Fill(eta1);
                    etaSelected_pt_0_2_SEPP->Fill(eta2);

                    hMassSEPP->Fill(mass);

                    //chi2 < 45
                    if (chi2MatchMCHMFT1 <= 45 && chi2MatchMCHMFT2 <= 45) {
                        chi2MatchMFTMCHetaSelected_pt_0_2_SEPP->Fill(chi2MatchMCHMFT1);
                        chi2MatchMFTMCHetaSelected_pt_0_2_SEPP->Fill(chi2MatchMCHMFT2);
                        hMassSEPPchi2_45_pt_0_2->Fill(mass);
                    }
                }

                //MFT eta range, -/-, pt of single muon > 0, 0 < pt < 2 GeV/c (dimuon) 
                else if ((eta1 >= -3.6 && eta1 <= -2.5) && (eta2 >= -3.6 && eta2 <= -2.5) && (sign1 == -1 && sign2 == -1) && (pt > 0 && pt < 2)) { //Muon quality cuts

                    etaSelected_pt_0_2_SEMM->Fill(eta1);
                    etaSelected_pt_0_2_SEMM->Fill(eta2);

                    hMassSEPM->Fill(mass);

                    //chi2 < 45
                    if (chi2MatchMCHMFT1 <= 45 && chi2MatchMCHMFT2 <= 45) {
                        chi2MatchMFTMCHetaSelected_pt_0_2_SEMM->Fill(chi2MatchMCHMFT1);
                        chi2MatchMFTMCHetaSelected_pt_0_2_SEMM->Fill(chi2MatchMCHMFT2);
                        hMassSEMMchi2_45_pt_0_2->Fill(mass);
                    }
                }

                //MFT eta range, unlike sign dimuon pair, pt of single muon > 0, 2 < pt < 4 GeV/c (dimuon)
                else if ((eta1 >= -3.6 && eta1 <= -2.5) && (eta2 >= -3.6 && eta2 <= -2.5) && ((sign1 == 1 && sign2 == -1) || (sign1 == -1 && sign2 == 1)) && (pt > 2 && pt < 4)) { //Muon quality cuts
                    
                    etaSelected_pt_2_4_SEPM->Fill(eta1);
                    etaSelected_pt_2_4_SEPM->Fill(eta2);
 
                    //chi2 < 45
                    if (chi2MatchMCHMFT1 <= 45 && chi2MatchMCHMFT2 <= 45) {
                        chi2MatchMFTMCHetaSelected_pt_2_4_SEPM->Fill(chi2MatchMCHMFT1);
                        chi2MatchMFTMCHetaSelected_pt_2_4_SEPM->Fill(chi2MatchMCHMFT2);
                        hMassSEPMchi2_45_pt_2_4->Fill(mass);
                    }
                }

                //MFT eta range, +/+, pt of single muon > 0, 2 < pt < 4 GeV/c (dimuon)
                else if ((eta1 >= -3.6 && eta1 <= -2.5) && (eta2 >= -3.6 && eta2 <= -2.5) && (sign1 == 1 && sign2 == 1) && (pt > 2 && pt < 4)) { //Muon quality cuts
                    
                    etaSelected_pt_2_4_SEPP->Fill(eta1);
                    etaSelected_pt_2_4_SEPP->Fill(eta2);
 
                    //chi2 < 45
                    if (chi2MatchMCHMFT1 <= 45 && chi2MatchMCHMFT2 <= 45) {
                        chi2MatchMFTMCHetaSelected_pt_2_4_SEPP->Fill(chi2MatchMCHMFT1);
                        chi2MatchMFTMCHetaSelected_pt_2_4_SEPP->Fill(chi2MatchMCHMFT2);
                        hMassSEPPchi2_45_pt_2_4->Fill(mass);
                    }
                }

                //MFT eta range, -/-, pt of single muon > 0, 2 < pt < 4 GeV/c (dimuon)
                else if ((eta1 >= -3.6 && eta1 <= -2.5) && (eta2 >= -3.6 && eta2 <= -2.5) && (sign1 == -1 && sign2 == -1) && (pt > 2 && pt < 4)) { //Muon quality cuts
                    
                    etaSelected_pt_2_4_SEMM->Fill(eta1);
                    etaSelected_pt_2_4_SEMM->Fill(eta2);
 
                    //chi2 < 45
                    if (chi2MatchMCHMFT1 <= 45 && chi2MatchMCHMFT2 <= 45) {
                        chi2MatchMFTMCHetaSelected_pt_2_4_SEMM->Fill(chi2MatchMCHMFT1);
                        chi2MatchMFTMCHetaSelected_pt_2_4_SEMM->Fill(chi2MatchMCHMFT2);
                        hMassSEMMchi2_45_pt_2_4->Fill(mass);
                    }
                }

                //MFT eta range, unlike sign dimuon pair, pt of single muon > 0, 4 < pt < 6 GeV/c (dimuon)
                else if ((eta1 >= -3.6 && eta1 <= -2.5) && (eta2 >= -3.6 && eta2 <= -2.5) && ((sign1 == 1 && sign2 == -1) || (sign1 == -1 && sign2 == 1)) && (pt > 4 && pt < 6)) { //Muon quality cuts
                    
                    etaSelected_pt_4_6_SEPM->Fill(eta1);
                    etaSelected_pt_4_6_SEPM->Fill(eta2);
                    
                    //chi2 < 45
                    if (chi2MatchMCHMFT1 <= 45 && chi2MatchMCHMFT2 <= 45) {
                        chi2MatchMFTMCHetaSelected_pt_4_6_SEPM->Fill(chi2MatchMCHMFT1);
                        chi2MatchMFTMCHetaSelected_pt_4_6_SEPM->Fill(chi2MatchMCHMFT2);
                        hMassSEPMchi2_45_pt_4_6->Fill(mass);
                    }
                }

                //MFT eta range, +/+, pt of single muon > 0, 4 < pt < 6 GeV/c (dimuon)
                else if ((eta1 >= -3.6 && eta1 <= -2.5) && (eta2 >= -3.6 && eta2 <= -2.5) && (sign1 == 1 && sign2 == 1) && (pt > 4 && pt < 6)) { //Muon quality cuts
                    
                    etaSelected_pt_4_6_SEPP->Fill(eta1);
                    etaSelected_pt_4_6_SEPP->Fill(eta2);
                    
                    //chi2 < 45
                    if (chi2MatchMCHMFT1 <= 45 && chi2MatchMCHMFT2 <= 45) {
                        chi2MatchMFTMCHetaSelected_pt_4_6_SEPP->Fill(chi2MatchMCHMFT1);
                        chi2MatchMFTMCHetaSelected_pt_4_6_SEPP->Fill(chi2MatchMCHMFT2);
                        hMassSEPPchi2_45_pt_4_6->Fill(mass);
                    }
                }

                //MFT eta range, -/-, pt of single muon > 0, 4 < pt < 6 GeV/c (dimuon)
                else if ((eta1 >= -3.6 && eta1 <= -2.5) && (eta2 >= -3.6 && eta2 <= -2.5) && (sign1 == -1 && sign2 == -1) && (pt > 4 && pt < 6)) { //Muon quality cuts
                    
                    etaSelected_pt_4_6_SEMM->Fill(eta1);
                    etaSelected_pt_4_6_SEMM->Fill(eta2);
                    
                    //chi2 < 45
                    if (chi2MatchMCHMFT1 <= 45 && chi2MatchMCHMFT2 <= 45) {
                        chi2MatchMFTMCHetaSelected_pt_4_6_SEMM->Fill(chi2MatchMCHMFT1);
                        chi2MatchMFTMCHetaSelected_pt_4_6_SEMM->Fill(chi2MatchMCHMFT2);
                        hMassSEMMchi2_45_pt_4_6->Fill(mass);
                    }
                }

                //MFT eta range, unlike sign dimuon pair, pt of single muon > 0, pt integrated (dimuon)
                if ((eta1 >= -3.6 && eta1 <= -2.5) && (eta2 >= -3.6 && eta2 <= -2.5) && ((sign1 == 1 && sign2 == -1) || (sign1 == -1 && sign2 == 1)) && pt < 20) { //Muon quality cuts
                    
                    etaSelected_pt_20_SEPM->Fill(eta1);
                    etaSelected_pt_20_SEPM->Fill(eta2);
                    
                    //chi2 < 45
                    if (chi2MatchMCHMFT1 <= 45 && chi2MatchMCHMFT2 <= 45) {
                        chi2MatchMFTMCHetaSelected_pt_20_SEPM->Fill(chi2MatchMCHMFT1);
                        chi2MatchMFTMCHetaSelected_pt_20_SEPM->Fill(chi2MatchMCHMFT2);
                        hMassSEPMchi2_45_pt_20->Fill(mass);
                    }
                }

                //MFT eta range, +/+, pt of single muon > 0, pt integrated (dimuon)
                else if ((eta1 >= -3.6 && eta1 <= -2.5) && (eta2 >= -3.6 && eta2 <= -2.5) && (sign1 == 1 && sign2 == 1) && pt < 20) { //Muon quality cuts
                    
                    etaSelected_pt_20_SEPP->Fill(eta1);
                    etaSelected_pt_20_SEPP->Fill(eta2);
                    
                    //chi2 < 45
                    if (chi2MatchMCHMFT1 <= 45 && chi2MatchMCHMFT2 <= 45) {
                        chi2MatchMFTMCHetaSelected_pt_20_SEPP->Fill(chi2MatchMCHMFT1);
                        chi2MatchMFTMCHetaSelected_pt_20_SEPP->Fill(chi2MatchMCHMFT2);
                        hMassSEPPchi2_45_pt_20->Fill(mass);
                    }
                }

                //MFT eta range, -/-, pt of single muon > 0, pt integrated (dimuon)
                else if ((eta1 >= -3.6 && eta1 <= -2.5) && (eta2 >= -3.6 && eta2 <= -2.5) && (sign1 == -1 && sign2 == -1) && pt < 20) { //Muon quality cuts
                    
                    etaSelected_pt_20_SEMM->Fill(eta1);
                    etaSelected_pt_20_SEMM->Fill(eta2);
                    
                    //chi2 < 45
                    if (chi2MatchMCHMFT1 <= 45 && chi2MatchMCHMFT2 <= 45) {
                        chi2MatchMFTMCHetaSelected_pt_20_SEMM->Fill(chi2MatchMCHMFT1);
                        chi2MatchMFTMCHetaSelected_pt_20_SEMM->Fill(chi2MatchMCHMFT2);
                        hMassSEMMchi2_45_pt_20->Fill(mass);
                    }
                }
            }
        }
        f->Close();
    }

    //Set draw option for histograms
    // +/-
    hMassSEPMchi2_45_pt_0_2->SetOption("HISTO");
    hMassSEPMchi2_45_pt_2_4->SetOption("HISTO");
    hMassSEPMchi2_45_pt_4_6->SetOption("HISTO");
    hMassSEPMchi2_45_pt_20->SetOption("HISTO");

    // +/+
    hMassSEPPchi2_45_pt_0_2->SetOption("HISTO");
    hMassSEPPchi2_45_pt_2_4->SetOption("HISTO");
    hMassSEPPchi2_45_pt_4_6->SetOption("HISTO");
    hMassSEPPchi2_45_pt_20->SetOption("HISTO");

    // -/-
    hMassSEMMchi2_45_pt_0_2->SetOption("HISTO");
    hMassSEMMchi2_45_pt_2_4->SetOption("HISTO");
    hMassSEMMchi2_45_pt_4_6->SetOption("HISTO");
    hMassSEMMchi2_45_pt_20->SetOption("HISTO");
    
    //Add to THStack
    // +/-
    stackMassChi2->Add(hMassSEPMchi2_45_pt_0_2);
    stackMassChi2->Add(hMassSEPMchi2_45_pt_2_4);
    stackMassChi2->Add(hMassSEPMchi2_45_pt_4_6);
    stackMassChi2->Add(hMassSEPMchi2_45_pt_20);

    // +/+
    stackMassChi2SEPP->Add(hMassSEPPchi2_45_pt_0_2);
    stackMassChi2SEPP->Add(hMassSEPPchi2_45_pt_2_4);
    stackMassChi2SEPP->Add(hMassSEPPchi2_45_pt_4_6);
    stackMassChi2SEPP->Add(hMassSEPPchi2_45_pt_20);

    // -/-
    stackMassChi2SEMM->Add(hMassSEMMchi2_45_pt_0_2);
    stackMassChi2SEMM->Add(hMassSEMMchi2_45_pt_2_4);
    stackMassChi2SEMM->Add(hMassSEMMchi2_45_pt_4_6);
    stackMassChi2SEMM->Add(hMassSEMMchi2_45_pt_20);

    //Draw THStacks of invariant mass
    TCanvas *cMassChi2 = new TCanvas(); // +/-
    cMassChi2->cd();
    stackMassChi2->Draw("nostack");
    lChi2->Draw("SAME");

    TCanvas *cMassChi2SEPP = new TCanvas(); // +/+
    cMassChi2SEPP->cd();
    stackMassChi2SEPP->Draw("nostack");
    lChi2->Draw("SAME");

    TCanvas *cMassChi2SEMM = new TCanvas(); // -/-
    cMassChi2SEMM->cd();
    stackMassChi2SEMM->Draw("nostack");
    lChi2->Draw("SAME");

    //Draw histograms of eta for pt ranges
    TCanvas *cEtaSelected_pt_0_2_SEPM = new TCanvas(); //0 < pt < 2 GeV/c +/-
    cEtaSelected_pt_0_2_SEPM->cd();
    etaSelected_pt_0_2_SEPM->Draw("HISTO");
    TCanvas *cEtaSelected_pt_0_2_SEPP = new TCanvas(); //0 < pt < 2 GeV/c +/+
    cEtaSelected_pt_0_2_SEPP->cd();
    etaSelected_pt_0_2_SEPP->Draw("HISTO");
    TCanvas *cEtaSelected_pt_0_2_SEMM = new TCanvas(); //0 < pt < 2 GeV/c -/-
    cEtaSelected_pt_0_2_SEMM->cd();
    etaSelected_pt_0_2_SEMM->Draw("HISTO");

    TCanvas *cEtaSelected_pt_2_4_SEPM = new TCanvas(); //2 < pt < 4 GeV/c +/-
    cEtaSelected_pt_2_4_SEPM->cd();
    etaSelected_pt_2_4_SEPM->Draw("HISTO");
    TCanvas *cEtaSelected_pt_2_4_SEPP = new TCanvas(); //2 < pt < 4 GeV/c +/+
    cEtaSelected_pt_2_4_SEPP->cd();
    etaSelected_pt_2_4_SEPP->Draw("HISTO");
    TCanvas *cEtaSelected_pt_2_4_SEMM = new TCanvas(); //2 < pt < 4 GeV/c -/-
    cEtaSelected_pt_2_4_SEMM->cd();
    etaSelected_pt_2_4_SEMM->Draw("HISTO");

    TCanvas *cEtaSelected_pt_4_6_SEPM = new TCanvas(); //4 < pt < 6 GeV/c +/-
    cEtaSelected_pt_4_6_SEPM->cd();
    etaSelected_pt_4_6_SEPM->Draw("HISTO");
    TCanvas *cEtaSelected_pt_4_6_SEPP = new TCanvas(); //4 < pt < 6 GeV/c +/+
    cEtaSelected_pt_4_6_SEPP->cd();
    etaSelected_pt_4_6_SEPP->Draw("HISTO");
    TCanvas *cEtaSelected_pt_4_6_SEMM = new TCanvas(); //4 < pt < 6 GeV/c -/-
    cEtaSelected_pt_4_6_SEMM->cd();
    etaSelected_pt_4_6_SEMM->Draw("HISTO");

    TCanvas *cEtaSelected_pt_20_SEPM = new TCanvas(); //pt integrated +/-
    cEtaSelected_pt_20_SEPM->cd();
    etaSelected_pt_20_SEPM->Draw("HISTO");
    TCanvas *cEtaSelected_pt_20_SEPP = new TCanvas(); //pt integrated +/+
    cEtaSelected_pt_20_SEPP->cd();
    etaSelected_pt_20_SEPP->Draw("HISTO");
    TCanvas *cEtaSelected_pt_20_SEMM = new TCanvas(); //pt integrated -/-
    cEtaSelected_pt_20_SEMM->cd();
    etaSelected_pt_20_SEMM->Draw("HISTO");

    //Draw histograms of chi1 for pt ranges
    TCanvas *cChi2MatchMFTMCHetaSelected_pt_0_2_SEPM = new TCanvas(); //0 < pt < 2 GeV/c +/-
    cChi2MatchMFTMCHetaSelected_pt_0_2_SEPM->cd();
    chi2MatchMFTMCHetaSelected_pt_0_2_SEPM->Draw("HISTO");
    TCanvas *cChi2MatchMFTMCHetaSelected_pt_0_2_SEPP = new TCanvas(); //0 < pt < 2 GeV/c +/+
    cChi2MatchMFTMCHetaSelected_pt_0_2_SEPP->cd();
    chi2MatchMFTMCHetaSelected_pt_0_2_SEPP->Draw("HISTO");
    TCanvas *cChi2MatchMFTMCHetaSelected_pt_0_2_SEMM = new TCanvas(); //0 < pt < 2 GeV/c -/-
    cChi2MatchMFTMCHetaSelected_pt_0_2_SEMM->cd();
    chi2MatchMFTMCHetaSelected_pt_0_2_SEMM->Draw("HISTO");

    TCanvas *cChi2MatchMFTMCHetaSelected_pt_2_4_SEPM = new TCanvas(); //2 < pt < 4 GeV/c +/-
    cChi2MatchMFTMCHetaSelected_pt_2_4_SEPM->cd();
    chi2MatchMFTMCHetaSelected_pt_2_4_SEPM->Draw("HISTO");
    TCanvas *cChi2MatchMFTMCHetaSelected_pt_2_4_SEPP = new TCanvas(); //2 < pt < 4 GeV/c +/+
    cChi2MatchMFTMCHetaSelected_pt_2_4_SEPP->cd();
    chi2MatchMFTMCHetaSelected_pt_2_4_SEPP->Draw("HISTO");
    TCanvas *cChi2MatchMFTMCHetaSelected_pt_2_4_SEMM = new TCanvas(); //2 < pt < 4 GeV/c -/-
    cChi2MatchMFTMCHetaSelected_pt_2_4_SEMM->cd();
    chi2MatchMFTMCHetaSelected_pt_2_4_SEMM->Draw("HISTO");
    
    TCanvas *cChi2MatchMFTMCHetaSelected_pt_4_6_SEPM = new TCanvas(); //4 < pt < 6 GeV/c +/-
    cChi2MatchMFTMCHetaSelected_pt_4_6_SEPM->cd();
    chi2MatchMFTMCHetaSelected_pt_4_6_SEPM->Draw("HISTO");
    TCanvas *cChi2MatchMFTMCHetaSelected_pt_4_6_SEPP = new TCanvas(); //4 < pt < 6 GeV/c +/+
    cChi2MatchMFTMCHetaSelected_pt_4_6_SEPP->cd();
    chi2MatchMFTMCHetaSelected_pt_4_6_SEPP->Draw("HISTO");
    TCanvas *cChi2MatchMFTMCHetaSelected_pt_4_6_SEMM = new TCanvas(); //4 < pt < 6 GeV/c -/-
    cChi2MatchMFTMCHetaSelected_pt_4_6_SEMM->cd();
    chi2MatchMFTMCHetaSelected_pt_4_6_SEMM->Draw("HISTO");

    TCanvas *cChi2MatchMFTMCHetaSelected_pt_20_SEPM = new TCanvas(); //pt integrated +/-
    cChi2MatchMFTMCHetaSelected_pt_20_SEPM->cd();
    chi2MatchMFTMCHetaSelected_pt_20_SEPM->Draw("HISTO");
    TCanvas *cChi2MatchMFTMCHetaSelected_pt_20_SEPP = new TCanvas(); //pt integrated +/+
    cChi2MatchMFTMCHetaSelected_pt_20_SEPP->cd();
    chi2MatchMFTMCHetaSelected_pt_20_SEPP->Draw("HISTO");
    TCanvas *cChi2MatchMFTMCHetaSelected_pt_20_SEMM = new TCanvas(); //pt integrated -/-
    cChi2MatchMFTMCHetaSelected_pt_20_SEMM->cd();
    chi2MatchMFTMCHetaSelected_pt_20_SEMM->Draw("HISTO");

    //Save output to root file
    fOut->cd();
    //eta
    cEtaSelected_pt_0_2_SEPM->Write("etaSelectedTracks_pt_0_2_SEPM");
    cEtaSelected_pt_0_2_SEPP->Write("etaSelectedTracks_pt_0_2_SEPP");
    cEtaSelected_pt_0_2_SEMM->Write("etaSelectedTracks_pt_0_2_SEMM");
    cEtaSelected_pt_2_4_SEPM->Write("etaSelectedTracks_pt_2_4_SEPM");
    cEtaSelected_pt_2_4_SEPP->Write("etaSelectedTracks_pt_2_4_SEPP");
    cEtaSelected_pt_2_4_SEMM->Write("etaSelectedTracks_pt_2_4_SEMM");
    cEtaSelected_pt_4_6_SEPM->Write("etaSelectedTracks_pt_4_6_SEPM");
    cEtaSelected_pt_4_6_SEPP->Write("etaSelectedTracks_pt_4_6_SEPP");
    cEtaSelected_pt_4_6_SEMM->Write("etaSelectedTracks_pt_4_6_SEMM");
    cEtaSelected_pt_20_SEPM->Write("etaSelectedTracks_pt_20_SEPM");
    cEtaSelected_pt_20_SEPP->Write("etaSelectedTracks_pt_20_SEPP");
    cEtaSelected_pt_20_SEMM->Write("etaSelectedTracks_pt_20_SEMM");
    //chi2
    cChi2MatchMFTMCHetaSelected_pt_0_2_SEPM->Write("chi2MatchEtaSelected_pt_0_2_SEPM");
    cChi2MatchMFTMCHetaSelected_pt_0_2_SEPP->Write("chi2MatchEtaSelected_pt_0_2_SEPP");
    cChi2MatchMFTMCHetaSelected_pt_0_2_SEMM->Write("chi2MatchEtaSelected_pt_0_2_SEMM");
    cChi2MatchMFTMCHetaSelected_pt_2_4_SEPM->Write("chi2MatchEtaSelected_pt_2_4_SEPM");
    cChi2MatchMFTMCHetaSelected_pt_2_4_SEPP->Write("chi2MatchEtaSelected_pt_2_4_SEPP");
    cChi2MatchMFTMCHetaSelected_pt_2_4_SEMM->Write("chi2MatchEtaSelected_pt_2_4_SEMM");
    cChi2MatchMFTMCHetaSelected_pt_4_6_SEPM->Write("chi2MatchEtaSelected_pt_4_6_SEPM");
    cChi2MatchMFTMCHetaSelected_pt_4_6_SEPP->Write("chi2MatchEtaSelected_pt_4_6_SEPP");
    cChi2MatchMFTMCHetaSelected_pt_4_6_SEMM->Write("chi2MatchEtaSelected_pt_4_6_SEMM");
    cChi2MatchMFTMCHetaSelected_pt_20_SEPM->Write("chi2MatchEtaSelected_pt_20_SEPM");
    cChi2MatchMFTMCHetaSelected_pt_20_SEPP->Write("chi2MatchEtaSelected_pt_20_SEPP");
    cChi2MatchMFTMCHetaSelected_pt_20_SEMM->Write("chi2MatchEtaSelected_pt_20_SEMM");
    //Invariant mass THStacks
    cMassChi2->Write("massChi2"); // +/-
    cMassChi2SEPP->Write("massChi2SEPP"); // +/+ 
    cMassChi2SEMM->Write("massChi2SEMM"); // -/-
    
    //Invariant mass spectra
    // +/-
    hMassSEPMchi2_45_pt_0_2->Write("massSEPMetaSelectedChi2_45_pt_0_2");
    hMassSEPMchi2_45_pt_2_4->Write("massSEPMetaSelectedChi2_45_pt_2_4");
    hMassSEPMchi2_45_pt_4_6->Write("massSEPMetaSelectedChi2_45_pt_4_6");
    hMassSEPMchi2_45_pt_20->Write("massSEPMetaSelectedChi2_45_pt_20");
    // +/+
    hMassSEPPchi2_45_pt_0_2->Write("massSEPPetaSelectedChi2_45_pt_0_2");
    hMassSEPPchi2_45_pt_2_4->Write("massSEPPetaSelectedChi2_45_pt_2_4");
    hMassSEPPchi2_45_pt_4_6->Write("massSEPPetaSelectedChi2_45_pt_4_6");
    hMassSEPPchi2_45_pt_20->Write("massSEPPetaSelectedChi2_45_pt_20");
    // -/-
    hMassSEMMchi2_45_pt_0_2->Write("massSEMMetaSelectedChi2_45_pt_0_2");
    hMassSEMMchi2_45_pt_2_4->Write("massSEMMetaSelectedChi2_45_pt_2_4");
    hMassSEMMchi2_45_pt_4_6->Write("massSEMMetaSelectedChi2_45_pt_4_6");
    hMassSEMMchi2_45_pt_20->Write("massSEMMetaSelectedChi2_45_pt_20");

    fOut->Close();

    cout << "Done!" << endl;
}