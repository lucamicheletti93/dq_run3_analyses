#include "TProfile.h"
#include <fstream>
#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include <vector>
#include <filesystem>
#include <cstdlib>

void merge() {
    std::string pathToFiles = "/Users/saragaretti/cernbox/JPsi_Psi2S_ratio/Train_295284/data/LHC23_pass4";
    std::ifstream fRunList(Form("%s/run_list_analysis.txt", pathToFiles.c_str()));
    /* std::string pathToFiles = "/Users/saragaretti/cernbox/JPsi_Psi2S_ratio/Train_321538/data/LHC23_pass4_pTcut_1GeV_new";
    std::ifstream fRunList(Form("%s/run_list_K.txt", pathToFiles.c_str())); */
    std::ofstream missingRunsFile(Form("%s/missing_runs.txt", pathToFiles.c_str()));
    //string muonCut = "muonLowPt10SigmaPDCA";
    string muonCut = "muonLowPt210SigmaPDCA";
    //string muonCut = "muonLowPt510SigmaPDCA";
    //string muonCut = "matchedMchMid";

    std::vector<int> runList;
    int runNumber;
    while (fRunList >> runNumber) {
        runList.push_back(runNumber);
    }

    const int nRuns = 200;
    const int nCentrBins = 4;
    double minCentrBins[] = {0.0, 20.0, 40.0, 60.0};
    double maxCentrBins[] = {20.0, 40.0, 60.0, 90.0};
    const int nPtBins = 4;
    double minPtBins[] = {0.0, 2.0, 4.0, 5.0};
    double maxPtBins[] = {2.0, 4.0, 5.0, 12.0};
     const int nPtBins2 = 2;
    double minPtBins2[] = {4.0, 6.0};
    double maxPtBins2[] = {6.0, 12.0};

    // Array for centrality integrated histograms
    TH1D* histSumIntSEPM = nullptr;
    TH1D* histSumIntSEPP = nullptr;
    TH1D* histSumIntSEMM = nullptr;
    TH1D* histSumIntMEPM = nullptr;
    TH1D* histSumIntMEPP = nullptr;
    TH1D* histSumIntMEMM = nullptr;

    // Array for pT integrated histograms
    TH1D* histSumIntPtSEPM = nullptr;
    TH1D* histSumIntPtSEPP = nullptr;
    TH1D* histSumIntPtSEMM = nullptr;
    TH1D* histSumIntPtMEPM = nullptr;
    TH1D* histSumIntPtMEPP = nullptr;
    TH1D* histSumIntPtMEMM = nullptr;

    // Array for rapidity integrated histograms
    TH1D* histSumIntRapSEPM = nullptr;
    TH1D* histSumIntRapSEPP = nullptr;
    TH1D* histSumIntRapSEMM = nullptr;
    TH1D* histSumIntRapMEPM = nullptr;
    TH1D* histSumIntRapMEPP = nullptr;
    TH1D* histSumIntRapMEMM = nullptr;

    // Array of histograms for centrality bins
    TH1D* histSumCentrSEPM[nCentrBins] = {nullptr};
    TH1D* histSumCentrSEPP[nCentrBins] = {nullptr};
    TH1D* histSumCentrSEMM[nCentrBins] = {nullptr};
    TH1D* histSumCentrMEPM[nCentrBins] = {nullptr};
    TH1D* histSumCentrMEPP[nCentrBins] = {nullptr};
    TH1D* histSumCentrMEMM[nCentrBins] = {nullptr};

    //pT
    TH1D* histSumCentrPtSEPM[nCentrBins] = {nullptr};
    TH1D* histSumCentrPtSEPP[nCentrBins] = {nullptr};
    TH1D* histSumCentrPtSEMM[nCentrBins] = {nullptr};
    TH1D* histSumCentrPtMEPM[nCentrBins] = {nullptr};
    TH1D* histSumCentrPtMEPP[nCentrBins] = {nullptr};
    TH1D* histSumCentrPtMEMM[nCentrBins] = {nullptr};

    //rapidity
    TH1D* histSumCentrRapSEPM[nCentrBins] = {nullptr};
    TH1D* histSumCentrRapSEPP[nCentrBins] = {nullptr};
    TH1D* histSumCentrRapSEMM[nCentrBins] = {nullptr};
    TH1D* histSumCentrRapMEPM[nCentrBins] = {nullptr};
    TH1D* histSumCentrRapMEPP[nCentrBins] = {nullptr};
    TH1D* histSumCentrRapMEMM[nCentrBins] = {nullptr};

    // Array of histograms for pT bins
    TH1D* histSumPtSEPM[nPtBins] = {nullptr};
    TH1D* histSumPtSEPP[nPtBins] = {nullptr};
    TH1D* histSumPtSEMM[nPtBins] = {nullptr};
    TH1D* histSumPtMEPM[nPtBins] = {nullptr};
    TH1D* histSumPtMEPP[nPtBins] = {nullptr};
    TH1D* histSumPtMEMM[nPtBins] = {nullptr};
    TH1D* histSumPtSEPM2[nPtBins] = {nullptr};
    TH1D* histSumPtSEPP2[nPtBins] = {nullptr};
    TH1D* histSumPtSEMM2[nPtBins] = {nullptr};
    TH1D* histSumPtMEPM2[nPtBins] = {nullptr};
    TH1D* histSumPtMEPP2[nPtBins] = {nullptr};
    TH1D* histSumPtMEMM2[nPtBins] = {nullptr};

    int index = 0;

    std::cout << "Number of runs in the list: " << runList.size() << std::endl;
    if (!fRunList.is_open()) {
        std::cerr << "Error opnening file run_list_analysis.txt" << std::endl;
        return;
    }
    for (const auto& run : runList) {
        cout << "index: " << index << endl;
        TFile* fIn = new TFile(Form("%s/%d/Histograms_%s.root", pathToFiles.c_str(), run, muonCut.c_str()), "READ");
        std::cout << "Doing the merging of run: " << run << std::endl;
        if (fIn->IsZombie()) {
            std::cout << "Run " << run << " -> ISSUE!" << std::endl;
            missingRunsFile << run << "\n";
            continue;
        }

        if (!fIn->GetListOfKeys()->Contains("Mass_Int_SEPM") || !fIn->GetListOfKeys()->Contains("Mass_Int_MEPM")) {
            std::cout << "Run " << run << " -> Missing SEPM or MEPM histogram, skipping this run." << std::endl;
            fIn->Close();
            missingRunsFile << run << "\n";
            continue;
        }

        TH1D* histIntSEPM = (TH1D*) fIn -> Get("Mass_Int_SEPM"); 
        TH1D* histIntSEPP = (TH1D*) fIn -> Get("Mass_Int_SEPP"); 
        TH1D* histIntSEMM = (TH1D*) fIn -> Get("Mass_Int_SEMM"); 
        TH1D* histIntMEPM = (TH1D*) fIn -> Get("Mass_Int_MEPM"); 
        TH1D* histIntMEPP = (TH1D*) fIn -> Get("Mass_Int_MEPP"); 
        TH1D* histIntMEMM = (TH1D*) fIn -> Get("Mass_Int_MEMM");
        //pT
        TH1D* histPtSEPM = (TH1D*) fIn -> Get("Mass_Pt_SEPM"); 
        TH1D* histPtSEPP = (TH1D*) fIn -> Get("Mass_Pt_SEPP"); 
        TH1D* histPtSEMM = (TH1D*) fIn -> Get("Mass_Pt_SEMM"); 
        TH1D* histPtMEPM = (TH1D*) fIn -> Get("Mass_Pt_MEPM"); 
        TH1D* histPtMEPP = (TH1D*) fIn -> Get("Mass_Pt_MEPP"); 
        TH1D* histPtMEMM = (TH1D*) fIn -> Get("Mass_Pt_MEMM");
        //rapidity
        TH1D* histRapSEPM = (TH1D*) fIn -> Get("Mass_Rap_SEPM"); 
        TH1D* histRapSEPP = (TH1D*) fIn -> Get("Mass_Rap_SEPP"); 
        TH1D* histRapSEMM = (TH1D*) fIn -> Get("Mass_Rap_SEMM"); 
        TH1D* histRapMEPM = (TH1D*) fIn -> Get("Mass_Rap_MEPM"); 
        TH1D* histRapMEPP = (TH1D*) fIn -> Get("Mass_Rap_MEPP"); 
        TH1D* histRapMEMM = (TH1D*) fIn -> Get("Mass_Rap_MEMM");

        // Sum of integrated histograms
        if (index == 0) {
            //centrality
            histSumIntSEPM = (TH1D*)histIntSEPM -> Clone("histSumInt_SEPM");
            histSumIntSEPM -> SetDirectory(0);
            histSumIntSEPP = (TH1D*)histIntSEPP -> Clone("histSumInt_SEPP");
            histSumIntSEPP -> SetDirectory(0);
            histSumIntSEMM = (TH1D*)histIntSEMM -> Clone("histSumInt_SEMM");
            histSumIntSEMM -> SetDirectory(0);
            histSumIntMEPM = (TH1D*)histIntMEPM -> Clone("histSumInt_MEPM");
            histSumIntMEPM -> SetDirectory(0);
            histSumIntMEPP = (TH1D*)histIntMEPP -> Clone("histSumInt_MEPP");
            histSumIntMEPP -> SetDirectory(0);
            histSumIntMEMM = (TH1D*)histIntMEMM -> Clone("histSumInt_MEMM");
            histSumIntMEMM -> SetDirectory(0);
            //pT
            histSumIntPtSEPM = (TH1D*)histPtSEPM -> Clone("histSumPt_SEPM");
            histSumIntPtSEPM -> SetDirectory(0);
            histSumIntPtSEPP = (TH1D*)histPtSEPP -> Clone("histSumPt_SEPP");
            histSumIntPtSEPP -> SetDirectory(0);
            histSumIntPtSEMM = (TH1D*)histPtSEMM -> Clone("histSumPt_SEMM");
            histSumIntPtSEMM -> SetDirectory(0);
            histSumIntPtMEPM = (TH1D*)histPtMEPM -> Clone("histSumPt_MEPM");
            histSumIntPtMEPM -> SetDirectory(0);
            histSumIntPtMEPP = (TH1D*)histPtMEPP -> Clone("histSumPt_MEPP");
            histSumIntPtMEPP -> SetDirectory(0);
            histSumIntPtMEMM = (TH1D*)histPtMEMM -> Clone("histSumPt_MEMM");
            histSumIntPtMEMM -> SetDirectory(0);
            //rapidity
            histSumIntRapSEPM = (TH1D*)histRapSEPM -> Clone("histSumRap_SEPM");
            histSumIntRapSEPM -> SetDirectory(0);
            histSumIntRapSEPP = (TH1D*)histRapSEPP -> Clone("histSumRap_SEPP");
            histSumIntRapSEPP -> SetDirectory(0);
            histSumIntRapSEMM = (TH1D*)histRapSEMM -> Clone("histSumRap_SEMM");
            histSumIntRapSEMM -> SetDirectory(0);
            histSumIntRapMEPM = (TH1D*)histRapMEPM -> Clone("histSumRap_MEPM");
            histSumIntRapMEPM -> SetDirectory(0);
            histSumIntRapMEPP = (TH1D*)histRapMEPP -> Clone("histSumRap_MEPP");
            histSumIntRapMEPP -> SetDirectory(0);
            histSumIntRapMEMM = (TH1D*)histRapMEMM -> Clone("histSumRap_MEMM");
            histSumIntRapMEMM -> SetDirectory(0);
        } 
        else {
            //centrality
            histSumIntSEPM -> Add(histIntSEPM);
            histSumIntSEPP -> Add(histIntSEPP);
            histSumIntSEMM -> Add(histIntSEMM);
            histSumIntMEPM -> Add(histIntMEPM);
            histSumIntMEPP -> Add(histIntMEPP);
            histSumIntMEMM -> Add(histIntMEMM);
            //pT
            histSumIntPtSEPM -> Add(histPtSEPM);
            histSumIntPtSEPP -> Add(histPtSEPP);
            histSumIntPtSEMM -> Add(histPtSEMM);
            histSumIntPtMEPM -> Add(histPtMEPM);
            histSumIntPtMEPP -> Add(histPtMEPP);
            histSumIntPtMEMM -> Add(histPtMEMM);
            //rapidity
            histSumIntRapSEPM -> Add(histRapSEPM);
            histSumIntRapSEPP -> Add(histRapSEPP);
            histSumIntRapSEMM -> Add(histRapSEMM);
            histSumIntRapMEPM -> Add(histRapMEPM);
            histSumIntRapMEPP -> Add(histRapMEPP);
            histSumIntRapMEMM -> Add(histRapMEMM);
            //Sumw2
            //
            //
            histSumIntSEPM -> Sumw2();
            histSumIntSEPP -> Sumw2();
            histSumIntSEMM -> Sumw2();
            histSumIntMEPM -> Sumw2();
            histSumIntMEPP -> Sumw2();
            histSumIntMEMM -> Sumw2();
            //pT
            histSumIntPtSEPM -> Sumw2();
            histSumIntPtSEPP -> Sumw2();
            histSumIntPtSEMM -> Sumw2();
            histSumIntPtMEPM -> Sumw2();
            histSumIntPtMEPP -> Sumw2();
            histSumIntPtMEMM -> Sumw2();
            //rapidity
            histSumIntRapSEPM -> Sumw2();
            histSumIntRapSEPP -> Sumw2();
            histSumIntRapSEMM -> Sumw2();
            histSumIntRapMEPM -> Sumw2();
            histSumIntRapMEPP -> Sumw2();
            histSumIntRapMEMM -> Sumw2();
        }


        // Sum of histograms for centrality bins
        for (int iCentr = 0; iCentr < nCentrBins; iCentr++) {

            TH1D* histCentrSEPM = (TH1D*) fIn -> Get(Form("Mass_CentrFT0C_%.0f_%.0f_SEPM", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histCentrSEPP = (TH1D*) fIn -> Get(Form("Mass_CentrFT0C_%.0f_%.0f_SEPP", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histCentrSEMM = (TH1D*) fIn -> Get(Form("Mass_CentrFT0C_%.0f_%.0f_SEMM", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histCentrMEPM = (TH1D*) fIn -> Get(Form("Mass_CentrFT0C_%.0f_%.0f_MEPM", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histCentrMEPP = (TH1D*) fIn -> Get(Form("Mass_CentrFT0C_%.0f_%.0f_MEPP", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histCentrMEMM = (TH1D*) fIn -> Get(Form("Mass_CentrFT0C_%.0f_%.0f_MEMM", minCentrBins[iCentr], maxCentrBins[iCentr]));
            //pT
            TH1D* histPtSEPM = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_SEPM", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histPtSEPP = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_SEPP", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histPtSEMM = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_SEMM", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histPtMEPM = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_MEPM", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histPtMEPP = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_MEPP", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histPtMEMM = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_MEMM", minCentrBins[iCentr], maxCentrBins[iCentr]));
            //rapidity
            TH1D* histRapSEPM = (TH1D*) fIn -> Get(Form("Mass_Rap_%.0f_%.0f_SEPM", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histRapSEPP = (TH1D*) fIn -> Get(Form("Mass_Rap_%.0f_%.0f_SEPP", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histRapSEMM = (TH1D*) fIn -> Get(Form("Mass_Rap_%.0f_%.0f_SEMM", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histRapMEPM = (TH1D*) fIn -> Get(Form("Mass_Rap_%.0f_%.0f_MEPM", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histRapMEPP = (TH1D*) fIn -> Get(Form("Mass_Rap_%.0f_%.0f_MEPP", minCentrBins[iCentr], maxCentrBins[iCentr])); 
            TH1D* histRapMEMM = (TH1D*) fIn -> Get(Form("Mass_Rap_%.0f_%.0f_MEMM", minCentrBins[iCentr], maxCentrBins[iCentr]));

            if (index == 0) {
                histSumCentrSEPM[iCentr] = (TH1D*)histCentrSEPM -> Clone(Form("histSumCentr_SEPM_%.0f", iCentr));
                histSumCentrSEPM[iCentr] -> SetDirectory(0);
                histSumCentrSEPP[iCentr] = (TH1D*)histCentrSEPP -> Clone(Form("histSumCentr_SEPP_%.0f", iCentr));
                histSumCentrSEPP[iCentr] -> SetDirectory(0);
                histSumCentrSEMM[iCentr] = (TH1D*)histCentrSEMM -> Clone(Form("histSumCentr_SEMM_%.0f", iCentr));
                histSumCentrSEMM[iCentr] -> SetDirectory(0);
                histSumCentrMEPM[iCentr] = (TH1D*)histCentrMEPM -> Clone(Form("histSumCentr_MEPM_%.0f", iCentr));
                histSumCentrMEPM[iCentr] -> SetDirectory(0);
                histSumCentrMEPP[iCentr] = (TH1D*)histCentrMEPP -> Clone(Form("histSumCentr_MEPP_%.0f", iCentr));
                histSumCentrMEPP[iCentr] -> SetDirectory(0);
                histSumCentrMEMM[iCentr] = (TH1D*)histCentrMEMM -> Clone(Form("histSumCentr_MEMM_%.0f", iCentr));
                histSumCentrMEMM[iCentr] -> SetDirectory(0);
                //pT
                histSumCentrPtSEPM[iCentr] = (TH1D*)histPtSEPM -> Clone(Form("histSumPt_SEPM_%.0f", iCentr));
                histSumCentrPtSEPM[iCentr] -> SetDirectory(0);
                histSumCentrPtSEPP[iCentr] = (TH1D*)histPtSEPP -> Clone(Form("histSumPt_SEPP_%.0f", iCentr));
                histSumCentrPtSEPP[iCentr] -> SetDirectory(0);
                histSumCentrPtSEMM[iCentr] = (TH1D*)histPtSEMM -> Clone(Form("histSumPt_SEMM_%.0f", iCentr));
                histSumCentrPtSEMM[iCentr] -> SetDirectory(0);
                histSumCentrPtMEPM[iCentr] = (TH1D*)histPtMEPM -> Clone(Form("histSumPt_MEPM_%.0f", iCentr));
                histSumCentrPtMEPM[iCentr] -> SetDirectory(0);
                histSumCentrPtMEPP[iCentr] = (TH1D*)histPtMEPP -> Clone(Form("histSumPt_MEPP_%.0f", iCentr));
                histSumCentrPtMEPP[iCentr] -> SetDirectory(0);
                histSumCentrPtMEMM[iCentr] = (TH1D*)histPtMEMM -> Clone(Form("histSumPt_MEMM_%.0f", iCentr));
                histSumCentrPtMEMM[iCentr] -> SetDirectory(0);
                //rapidity
                histSumCentrRapSEPM[iCentr] = (TH1D*)histRapSEPM -> Clone(Form("histSumRap_SEPM_%.0f", iCentr));
                histSumCentrRapSEPM[iCentr] -> SetDirectory(0);
                histSumCentrRapSEPP[iCentr] = (TH1D*)histRapSEPP -> Clone(Form("histSumRap_SEPP_%.0f", iCentr));
                histSumCentrRapSEPP[iCentr] -> SetDirectory(0);
                histSumCentrRapSEMM[iCentr] = (TH1D*)histRapSEMM -> Clone(Form("histSumRap_SEMM_%.0f", iCentr));
                histSumCentrRapSEMM[iCentr] -> SetDirectory(0);
                histSumCentrRapMEPM[iCentr] = (TH1D*)histRapMEPM -> Clone(Form("histSumRap_MEPM_%.0f", iCentr));
                histSumCentrRapMEPM[iCentr] -> SetDirectory(0);
                histSumCentrRapMEPP[iCentr] = (TH1D*)histRapMEPP -> Clone(Form("histSumRap_MEPP_%.0f", iCentr));
                histSumCentrRapMEPP[iCentr] -> SetDirectory(0);
                histSumCentrRapMEMM[iCentr] = (TH1D*)histRapMEMM -> Clone(Form("histSumRap_MEMM_%.0f", iCentr));
                histSumCentrRapMEMM[iCentr] -> SetDirectory(0);
            } 
            else {
                histSumCentrSEPM[iCentr] -> Add(histCentrSEPM);
                histSumCentrSEPP[iCentr] -> Add(histCentrSEPP);
                histSumCentrSEMM[iCentr] -> Add(histCentrSEMM);
                histSumCentrMEPM[iCentr] -> Add(histCentrMEPM);
                histSumCentrMEPP[iCentr] -> Add(histCentrMEPP);
                histSumCentrMEMM[iCentr] -> Add(histCentrMEMM);
                //pT
                histSumCentrPtSEPM[iCentr] -> Add(histPtSEPM);
                histSumCentrPtSEPP[iCentr] -> Add(histPtSEPP);
                histSumCentrPtSEMM[iCentr] -> Add(histPtSEMM);
                histSumCentrPtMEPM[iCentr] -> Add(histPtMEPM);
                histSumCentrPtMEPP[iCentr] -> Add(histPtMEPP);
                histSumCentrPtMEMM[iCentr] -> Add(histPtMEMM);
                //rapidity
                histSumCentrRapSEPM[iCentr] -> Add(histRapSEPM);
                histSumCentrRapSEPP[iCentr] -> Add(histRapSEPP);
                histSumCentrRapSEMM[iCentr] -> Add(histRapSEMM);
                histSumCentrRapMEPM[iCentr] -> Add(histRapMEPM);
                histSumCentrRapMEPP[iCentr] -> Add(histRapMEPP);
                histSumCentrRapMEMM[iCentr] -> Add(histRapMEMM);
            }
        }

        // Sum of histograms for pT bins
        for (int iPt = 0; iPt < nPtBins; iPt++) {
            cout << "iPt: " << iPt << endl;

            TH1D* histPtSEPM = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_SEPM", minPtBins[iPt], maxPtBins[iPt])); 
            TH1D* histPtSEPP = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_SEPP", minPtBins[iPt], maxPtBins[iPt])); 
            TH1D* histPtSEMM = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_SEMM", minPtBins[iPt], maxPtBins[iPt])); 
            TH1D* histPtMEPM = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_MEPM", minPtBins[iPt], maxPtBins[iPt])); 
            TH1D* histPtMEPP = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_MEPP", minPtBins[iPt], maxPtBins[iPt])); 
            TH1D* histPtMEMM = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_MEMM", minPtBins[iPt], maxPtBins[iPt]));

            if (index == 0) {
                histSumPtSEPM[iPt] = (TH1D*)histPtSEPM -> Clone(Form("histSumPt_SEPM_%.0f", iPt));
                histSumPtSEPM[iPt] -> SetDirectory(0);
                histSumPtSEPP[iPt] = (TH1D*)histPtSEPP -> Clone(Form("histSumPt_SEPP_%.0f", iPt));
                histSumPtSEPP[iPt] -> SetDirectory(0);
                histSumPtSEMM[iPt] = (TH1D*)histPtSEMM -> Clone(Form("histSumPt_SEMM_%.0f", iPt));
                histSumPtSEMM[iPt] -> SetDirectory(0);
                histSumPtMEPM[iPt] = (TH1D*)histPtMEPM -> Clone(Form("histSumPt_MEPM_%.0f", iPt));
                histSumPtMEPM[iPt] -> SetDirectory(0);
                histSumPtMEPP[iPt] = (TH1D*)histPtMEPP -> Clone(Form("histSumPt_MEPP_%.0f", iPt));
                histSumPtMEPP[iPt] -> SetDirectory(0);
                histSumPtMEMM[iPt] = (TH1D*)histPtMEMM -> Clone(Form("histSumPt_MEMM_%.0f", iPt));
                histSumPtMEMM[iPt] -> SetDirectory(0);
            } 
            else {
                histSumPtSEPM[iPt] -> Add(histPtSEPM);
                histSumPtSEPP[iPt] -> Add(histPtSEPP);
                histSumPtSEMM[iPt] -> Add(histPtSEMM);
                histSumPtMEPM[iPt] -> Add(histPtMEPM);
                histSumPtMEPP[iPt] -> Add(histPtMEPP);
                histSumPtMEMM[iPt] -> Add(histPtMEMM);
            } 
        }
        for (int iPt = 0; iPt < nPtBins2; iPt++) {
            cout << "iPt2: " << iPt << endl;

            TH1D* histPtSEPM = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_SEPM", minPtBins[iPt], maxPtBins[iPt])); 
            TH1D* histPtSEPP = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_SEPP", minPtBins[iPt], maxPtBins[iPt])); 
            TH1D* histPtSEMM = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_SEMM", minPtBins[iPt], maxPtBins[iPt])); 
            TH1D* histPtMEPM = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_MEPM", minPtBins[iPt], maxPtBins[iPt])); 
            TH1D* histPtMEPP = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_MEPP", minPtBins[iPt], maxPtBins[iPt])); 
            TH1D* histPtMEMM = (TH1D*) fIn -> Get(Form("Mass_Pt_%.0f_%.0f_MEMM", minPtBins[iPt], maxPtBins[iPt]));

            if (index == 0) {
                histSumPtSEPM2[iPt] = (TH1D*)histPtSEPM -> Clone(Form("histSumPt_SEPM_%.0f", iPt));
                histSumPtSEPM2[iPt] -> SetDirectory(0);
                histSumPtSEPP2[iPt] = (TH1D*)histPtSEPP -> Clone(Form("histSumPt_SEPP_%.0f", iPt));
                histSumPtSEPP2[iPt] -> SetDirectory(0);
                histSumPtSEMM2[iPt] = (TH1D*)histPtSEMM -> Clone(Form("histSumPt_SEMM_%.0f", iPt));
                histSumPtSEMM2[iPt] -> SetDirectory(0);
                histSumPtMEPM2[iPt] = (TH1D*)histPtMEPM -> Clone(Form("histSumPt_MEPM_%.0f", iPt));
                histSumPtMEPM2[iPt] -> SetDirectory(0);
                histSumPtMEPP2[iPt] = (TH1D*)histPtMEPP -> Clone(Form("histSumPt_MEPP_%.0f", iPt));
                histSumPtMEPP2[iPt] -> SetDirectory(0);
                histSumPtMEMM2[iPt] = (TH1D*)histPtMEMM -> Clone(Form("histSumPt_MEMM_%.0f", iPt));
                histSumPtMEMM2[iPt] -> SetDirectory(0);
            } 
            else {
                histSumPtSEPM2[iPt] -> Add(histPtSEPM);
                histSumPtSEPP2[iPt] -> Add(histPtSEPP);
                histSumPtSEMM2[iPt] -> Add(histPtSEMM);
                histSumPtMEPM2[iPt] -> Add(histPtMEPM);
                histSumPtMEPP2[iPt] -> Add(histPtMEPP);
                histSumPtMEMM2[iPt] -> Add(histPtMEMM);
            } 
        }
        // TH1D* histPtSEPMRun = (TH1D*)histPtSEPM -> Clone();
        std::cout << "Run " << run << " -> OK" << std::endl;
        fIn -> Close();
        index ++;
    }

    // Salvataggio in un file di output .root
    TFile* fOut = new TFile(Form("%s/mergedHistograms_Sumw2.root", pathToFiles.c_str()), "RECREATE");
    histSumIntSEPM -> Write("Mass_Int_SEPM");
    histSumIntSEPP -> Write("Mass_Int_SEPP");
    histSumIntSEMM -> Write("Mass_Int_SEMM");
    histSumIntMEPM -> Write("Mass_Int_MEPM");
    histSumIntMEPP -> Write("Mass_Int_MEPP");
    histSumIntMEMM -> Write("Mass_Int_MEMM");

    //pT
    histSumIntPtSEPM -> Write("Mass_Pt_SEPM");
    histSumIntPtSEPP -> Write("Mass_Pt_SEPP");
    histSumIntPtSEMM -> Write("Mass_Pt_SEMM");
    histSumIntPtMEPM -> Write("Mass_Pt_MEPM");
    histSumIntPtMEPP -> Write("Mass_Pt_MEPP");
    histSumIntPtMEMM -> Write("Mass_Pt_MEMM");

    //rapidity
    histSumIntRapSEPM -> Write("Mass_Rap_SEPM");
    histSumIntRapSEPP -> Write("Mass_Rap_SEPP");
    histSumIntRapSEMM -> Write("Mass_Rap_SEMM");
    histSumIntRapMEPM -> Write("Mass_Rap_MEPM");
    histSumIntRapMEPP -> Write("Mass_Rap_MEPP");
    histSumIntRapMEMM -> Write("Mass_Rap_MEMM");

    for (int iCentr = 0; iCentr < nCentrBins; iCentr++) {
        histSumCentrSEPM[iCentr] -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEPM", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrSEPP[iCentr] -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEPP", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrSEMM[iCentr] -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEMM", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrMEPM[iCentr] -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEPM", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrMEPP[iCentr] -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEPP", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrMEMM[iCentr] -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEMM", minCentrBins[iCentr], maxCentrBins[iCentr]));
        //pT
        histSumCentrPtSEPM[iCentr] -> Write(Form("Mass_Pt_%.0f_%.0f_SEPM", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrPtSEPP[iCentr] -> Write(Form("Mass_Pt_%.0f_%.0f_SEPP", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrPtSEMM[iCentr] -> Write(Form("Mass_Pt_%.0f_%.0f_SEMM", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrPtMEPM[iCentr] -> Write(Form("Mass_Pt_%.0f_%.0f_MEPM", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrPtMEPP[iCentr] -> Write(Form("Mass_Pt_%.0f_%.0f_MEPP", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrPtMEMM[iCentr] -> Write(Form("Mass_Pt_%.0f_%.0f_MEMM", minCentrBins[iCentr], maxCentrBins[iCentr]));
        //rapidity
        histSumCentrRapSEPM[iCentr] -> Write(Form("Mass_Rap_%.0f_%.0f_SEPM", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrRapSEPP[iCentr] -> Write(Form("Mass_Rap_%.0f_%.0f_SEPP", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrRapSEMM[iCentr] -> Write(Form("Mass_Rap_%.0f_%.0f_SEMM", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrRapMEPM[iCentr] -> Write(Form("Mass_Rap_%.0f_%.0f_MEPM", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrRapMEPP[iCentr] -> Write(Form("Mass_Rap_%.0f_%.0f_MEPP", minCentrBins[iCentr], maxCentrBins[iCentr]));
        histSumCentrRapMEMM[iCentr] -> Write(Form("Mass_Rap_%.0f_%.0f_MEMM", minCentrBins[iCentr], maxCentrBins[iCentr]));
    }
    for (int iPt = 0; iPt < nPtBins; iPt++) {
        histSumPtSEPM[iPt] -> Write(Form("Mass_Pt_%.0f_%.0f_SEPM", minPtBins[iPt], maxPtBins[iPt]));
        histSumPtSEPP[iPt] -> Write(Form("Mass_Pt_%.0f_%.0f_SEPP", minPtBins[iPt], maxPtBins[iPt]));
        histSumPtSEMM[iPt] -> Write(Form("Mass_Pt_%.0f_%.0f_SEMM", minPtBins[iPt], maxPtBins[iPt]));
        histSumPtMEPM[iPt] -> Write(Form("Mass_Pt_%.0f_%.0f_MEPM", minPtBins[iPt], maxPtBins[iPt]));
        histSumPtMEPP[iPt] -> Write(Form("Mass_Pt_%.0f_%.0f_MEPP", minPtBins[iPt], maxPtBins[iPt]));
        histSumPtMEMM[iPt] -> Write(Form("Mass_Pt_%.0f_%.0f_MEMM", minPtBins[iPt], maxPtBins[iPt]));
    }
    for (int iPt = 0; iPt < nPtBins2; iPt++) {
        histSumPtSEPM2[iPt] -> Write(Form("Mass_Pt_%.0f_%.0f_SEPM", minPtBins2[iPt], maxPtBins2[iPt]));
        histSumPtSEPP2[iPt] -> Write(Form("Mass_Pt_%.0f_%.0f_SEPP", minPtBins2[iPt], maxPtBins2[iPt]));
        histSumPtSEMM2[iPt] -> Write(Form("Mass_Pt_%.0f_%.0f_SEMM", minPtBins2[iPt], maxPtBins2[iPt]));
        histSumPtMEPM2[iPt] -> Write(Form("Mass_Pt_%.0f_%.0f_MEPM", minPtBins2[iPt], maxPtBins2[iPt]));
        histSumPtMEPP2[iPt] -> Write(Form("Mass_Pt_%.0f_%.0f_MEPP", minPtBins2[iPt], maxPtBins2[iPt]));
        histSumPtMEMM2[iPt] -> Write(Form("Mass_Pt_%.0f_%.0f_MEMM", minPtBins2[iPt], maxPtBins2[iPt]));
    }
    //fOut -> ls();
    fOut -> Close();
    std::cout << "Merged histograms saved in mergedHistograms.root" << std::endl;
}
