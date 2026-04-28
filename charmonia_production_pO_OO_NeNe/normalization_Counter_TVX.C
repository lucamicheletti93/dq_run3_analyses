#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TLegend.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <string>
#include "THashList.h"

// xSecTVX cross section
//const double xSecTVX = 50.3e9; // pb-1, pp@5.36TeV, taken from https://alice-notes.web.cern.ch/system/files/notes/analysis/1671/2025-09-04-Analysis_Note_Pi0_OO.pdf, page 14
const double xSecTVX = 1.13e12; // pb-1, OO@5.36TeV, taken from https://alice-notes.web.cern.ch/system/files/notes/analysis/1671/2025-09-04-Analysis_Note_Pi0_OO.pdf, page 14

void RetrieveTriggerInfo(TString , bool , string, double [11]);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_normalization_from_single_file(string year = "2024", string period = "LHC24aq_pass1_small", string triggerMask = "minBias", string assocType = "std_assoc") {
    string fInName = "/home/rebecca/cernbox/O_O/Pass2/Train_706901/AnalysisResults_merged.root"; // LHC25ae_pass2, std assoc
    string fOutName = Form("data/%s/pass2/Train_706901/normalization/%s_%s_%s_trigger_summary_Train_706901_Sel8.root", year.c_str(), period.c_str(), triggerMask.c_str(), assocType.c_str());    

    TFile *fIn = TFile::Open(fInName.c_str());
    if (!fIn || fIn -> IsZombie()) {
        return;
    }

    TH1D *histEvSelCollisionBeforeCuts_TMaker = nullptr;
    TH1D *histEvSelCollisionAfterCuts_TMaker = nullptr;

    TH1D *histEvSelCollisionBeforeCuts_TReader = nullptr;
    TH1D *histEvSelCollisionAfterCuts_TReader = nullptr;

    TH1D *histBcSelCounterTVX = nullptr;
    TH1D *histBcSelCounterTVXafterBCcuts = nullptr;
    TH1D *histEvSelColCounterAll = nullptr;
    TH1D *histEvSelColCounterTVX = nullptr;

    TH1D *histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker = nullptr; // same as 'histEvSelCollisionBeforeCuts_TMaker' but as a function of centr (The bin '0-100%' of 'histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker' should have exactly the same value of 'histEvSelCollisionBeforeCuts_TMaker')
    TH1D *histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker = nullptr;
    TH1D *histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker = nullptr; // same as 'histEvSelCollisionAfterCuts_TMaker' but as a function of centr (The bin '0-100%' of 'histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker' should have exactly the same value of 'histEvSelCollisionAfterCuts_TMaker')
    TH1D *histRatioIntegrals_CentrFT0C_AfterCuts_TMaker = nullptr;

    TH1D *histSelectedIntegrals_CentrFT0C_BeforeCuts_TReader = nullptr; // same as 'histEvSelCollisionBeforeCuts_TReader' but as a function of centr (The bin '0-100%' of 'histSelectedIntegrals_CentrFT0C_BeforeCuts_TReader' should have exactly the same value of 'histEvSelCollisionBeforeCuts_TReader')
    TH1D *histRatioIntegrals_CentrFT0C_BeforeCuts_TReader = nullptr;
    TH1D *histSelectedIntegrals_CentrFT0C_AfterCuts_TReader = nullptr; // same as 'histEvSelCollisionAfterCuts_TReader' but as a function of centr (The bin '0-100%' of 'histSelectedIntegrals_CentrFT0C_AfterCuts_TReader' should have exactly the same value of 'histEvSelCollisionAfterCuts_TReader')
    TH1D *histRatioIntegrals_CentrFT0C_AfterCuts_TReader = nullptr;


    for (auto const& dirKey : *fIn -> GetListOfKeys()) {

        if (TString(dirKey -> GetName()).Contains("eventselection-run3")) {
            TH1D *histCounterTVX = (TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVX");
            TH1D *histCounterTVXafterBCcuts = (TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVXafterBCcuts");
            TH1D *histColCounterAll = (TH1D*) fIn -> Get("eventselection-run3/eventselection/hColCounterAll");
            TH1D *histColCounterTVX = (TH1D*) fIn -> Get("eventselection-run3/eventselection/hColCounterTVX");

            int nRuns = histCounterTVX -> GetXaxis() -> GetNbins();
            vector <string> vecRunList;
            for (int iRun = 0;iRun < nRuns;iRun++) {
                string tmpRunNumber = histCounterTVX -> GetXaxis() -> GetBinLabel(iRun+1);
                if (!(tmpRunNumber.empty())) {
                    vecRunList.push_back(histCounterTVX -> GetXaxis() -> GetBinLabel(iRun+1));
                }
            }

            histBcSelCounterTVX = new TH1D("histBcSelCounterTVX", "", vecRunList.size(), 0, vecRunList.size());
            histBcSelCounterTVXafterBCcuts = new TH1D("histBcSelCounterTVXafterBCcuts", "", vecRunList.size(), 0, vecRunList.size());
            histEvSelColCounterAll = new TH1D("histEvSelColCounterAll", "", vecRunList.size(), 0, vecRunList.size());
            histEvSelColCounterTVX = new TH1D("histEvSelColCounterTVX", "", vecRunList.size(), 0, vecRunList.size());
            

            for (int iRun = 0;iRun < int (vecRunList.size());iRun++) {
                string runNumber = vecRunList.at(iRun).c_str();
                histBcSelCounterTVX -> SetBinContent(iRun+1, histCounterTVX -> GetBinContent(iRun+1));
                histBcSelCounterTVX -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());

                histBcSelCounterTVXafterBCcuts -> SetBinContent(iRun+1, histCounterTVXafterBCcuts -> GetBinContent(iRun+1));
                histBcSelCounterTVXafterBCcuts -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());

                histEvSelColCounterAll -> SetBinContent(iRun+1, histColCounterAll -> GetBinContent(iRun+1));
                histEvSelColCounterAll -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());

                histEvSelColCounterTVX -> SetBinContent(iRun+1, histColCounterTVX -> GetBinContent(iRun+1));
                histEvSelColCounterTVX -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());
            }
        }
               
        if (TString(dirKey -> GetName()).Contains("table-maker")) {
            fIn->cd("table-maker"); 
            THashList* hashList_BeforeCuts_TMaker = (THashList*)gDirectory->Get("output"); 
            TList* list_output_BeforeCuts_TMaker = (TList*)hashList_BeforeCuts_TMaker->FindObject("Event_BeforeCuts");
            TH1D *histCentFTOC_BeforeCuts_TMaker = (TH1D*) list_output_BeforeCuts_TMaker->FindObject("CentFT0C");

            histEvSelCollisionBeforeCuts_TMaker = new TH1D("histEvSelCollisionBeforeCuts_TMaker", "", 1, 0, 1);

            double total_integral_BeforeCuts_TMaker = histCentFTOC_BeforeCuts_TMaker->Integral();
            histEvSelCollisionBeforeCuts_TMaker -> SetBinContent(1, total_integral_BeforeCuts_TMaker);

            //const int centMin[] = { 0,10,20,30,40,50,60,70, 90, 0,  0}; 
            //const int centMax[] = {10,20,30,40,50,60,70,90,100,90,100};
            const int centMin[] = { 0,10,20,30,40,50,60,70,80,90, 0,  0}; 
            const int centMax[] = {10,20,30,40,50,60,70,80,90,100,90,100};

            int nCentBins = sizeof(centMin) / sizeof(int);

            histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker = new TH1D("histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker", "SelectedIntegral of CentrFT0C BeforeCuts TMaker;Centrality interval;Value", nCentBins, 0, nCentBins);
            histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker = new TH1D("histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker", "SelectedIntegral/TotalIntegral of CentrFT0C BeforeCuts TMaker;Centrality interval;Fraction", nCentBins, 0, nCentBins);
            
            for (int i = 0; i < nCentBins; i++) {
                double selected_integral_BeforeCuts_TMaker = histCentFTOC_BeforeCuts_TMaker->Integral(centMin[i]+1, centMax[i]);
                double ratio_BeforeCuts_TMaker = selected_integral_BeforeCuts_TMaker / total_integral_BeforeCuts_TMaker;

                histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker->SetBinContent(i + 1, selected_integral_BeforeCuts_TMaker);
                histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker->SetBinContent(i + 1, ratio_BeforeCuts_TMaker);

                TString label_selected_integral_BeforeCuts_TMaker = Form("%d-%d", centMin[i], centMax[i]);
                TString label_ratio_integral_BeforeCuts_TMaker = Form("%d-%d", centMin[i], centMax[i]);
                histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker->GetXaxis()->SetBinLabel(i + 1, label_selected_integral_BeforeCuts_TMaker);
                histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker->GetXaxis()->SetBinLabel(i + 1, label_ratio_integral_BeforeCuts_TMaker);

            }

            fIn->cd("table-maker"); 
            THashList* hashList_AfterCuts_TMaker = (THashList*)gDirectory->Get("output"); 
            TList* list_output_AfterCuts_TMaker = (TList*)hashList_AfterCuts_TMaker->FindObject("Event_AfterCuts");
            TH1D *histCentFTOC_AfterCuts_TMaker = (TH1D*) list_output_AfterCuts_TMaker->FindObject("CentFT0C");

            histEvSelCollisionAfterCuts_TMaker = new TH1D("histEvSelCollisionAfterCuts_TMaker", "", 1, 0, 1);

            double total_integral_AfterCuts_TMaker = histCentFTOC_AfterCuts_TMaker->Integral();
            histEvSelCollisionAfterCuts_TMaker -> SetBinContent(1, total_integral_AfterCuts_TMaker);

            // stessi centr bin

            histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker = new TH1D("histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker", "SelectedIntegral of CentrFT0C AfterCuts TMaker;Centrality interval;Value", nCentBins, 0, nCentBins);
            histRatioIntegrals_CentrFT0C_AfterCuts_TMaker = new TH1D("histRatioIntegrals_CentrFT0C_AfterCuts_TMaker", "SelectedIntegral/TotalIntegral of CentrFT0C AfterCuts TMaker;Centrality interval;Fraction", nCentBins, 0, nCentBins);
            
            for (int i = 0; i < nCentBins; i++) {
                double selected_integral_AfterCuts_TMaker = histCentFTOC_AfterCuts_TMaker->Integral(centMin[i]+1, centMax[i]);
                double ratio_AfterCuts_TMaker = selected_integral_AfterCuts_TMaker / total_integral_AfterCuts_TMaker;

                histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker->SetBinContent(i + 1, selected_integral_AfterCuts_TMaker);
                histRatioIntegrals_CentrFT0C_AfterCuts_TMaker->SetBinContent(i + 1, ratio_AfterCuts_TMaker);

                TString label_selected_integral_AfterCuts_TMaker = Form("%d-%d", centMin[i], centMax[i]);
                TString label_ratio_integral_AfterCuts_TMaker = Form("%d-%d", centMin[i], centMax[i]);
                histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker->GetXaxis()->SetBinLabel(i + 1, label_selected_integral_AfterCuts_TMaker);
                histRatioIntegrals_CentrFT0C_AfterCuts_TMaker->GetXaxis()->SetBinLabel(i + 1, label_ratio_integral_AfterCuts_TMaker);

            }
        }

        if (TString(dirKey -> GetName()).Contains("analysis-event-selection")) {
            fIn->cd("analysis-event-selection"); 
            THashList* hashList_BeforeCuts_TReader = (THashList*)gDirectory->Get("output"); 
            TList* list_output_BeforeCuts_TReader = (TList*)hashList_BeforeCuts_TReader->FindObject("Event_BeforeCuts");
            TH1D *histCentFTOC_BeforeCuts_TReader = (TH1D*) list_output_BeforeCuts_TReader->FindObject("CentFT0C");

            //const int centMin[] = { 0,10,20,30,40,50,60,70, 90, 0,  0}; 
            //const int centMax[] = {10,20,30,40,50,60,70,90,100,90,100};
            const int centMin[] = { 0,10,20,30,40,50,60,70,80,90, 0,  0}; 
            const int centMax[] = {10,20,30,40,50,60,70,80,90,100,90,100};
                
            int nCentBins = sizeof(centMin) / sizeof(int);

            histEvSelCollisionBeforeCuts_TReader = new TH1D("histEvSelCollisionBeforeCuts_TReader", "", 1, 0, 1);

            double total_integral_BeforeCuts_TReader = histCentFTOC_BeforeCuts_TReader->Integral();
            histEvSelCollisionBeforeCuts_TReader -> SetBinContent(1, total_integral_BeforeCuts_TReader);

            histSelectedIntegrals_CentrFT0C_BeforeCuts_TReader = new TH1D("histSelectedIntegrals_CentrFT0C_BeforeCuts_TReader", "SelectedIntegral of CentrFT0C BeforeCuts TReader;Centrality interval;Value", nCentBins, 0, nCentBins);
            histRatioIntegrals_CentrFT0C_BeforeCuts_TReader = new TH1D("histRatioIntegrals_CentrFT0C_BeforeCuts_TReader", "SelectedIntegral/TotalIntegral of CentrFT0C BeforeCuts TReader;Centrality interval;Fraction", nCentBins, 0, nCentBins);
            
            for (int i = 0; i < nCentBins; i++) {
                double selected_integral_BeforeCuts_TReader = histCentFTOC_BeforeCuts_TReader->Integral(centMin[i]+1, centMax[i]);
                double ratio_BeforeCuts_TReader = selected_integral_BeforeCuts_TReader / total_integral_BeforeCuts_TReader;

                histSelectedIntegrals_CentrFT0C_BeforeCuts_TReader->SetBinContent(i + 1, selected_integral_BeforeCuts_TReader);
                histRatioIntegrals_CentrFT0C_BeforeCuts_TReader->SetBinContent(i + 1, ratio_BeforeCuts_TReader);

                TString label_selected_integral_BeforeCuts_TReader = Form("%d-%d", centMin[i], centMax[i]);
                TString label_ratio_integral_BeforeCuts_TReader = Form("%d-%d", centMin[i], centMax[i]);
                histSelectedIntegrals_CentrFT0C_BeforeCuts_TReader->GetXaxis()->SetBinLabel(i + 1, label_selected_integral_BeforeCuts_TReader);
                histRatioIntegrals_CentrFT0C_BeforeCuts_TReader->GetXaxis()->SetBinLabel(i + 1, label_ratio_integral_BeforeCuts_TReader);

            }



            fIn->cd("analysis-event-selection");
            THashList* hashList_AfterCuts_TReader = (THashList*)gDirectory->Get("output"); 
            TList* list_output_AfterCuts_TReader = (TList*)hashList_AfterCuts_TReader->FindObject("Event_AfterCuts");
            TH1D *histCentFTOC_AfterCuts_TReader = (TH1D*) list_output_AfterCuts_TReader->FindObject("CentFT0C");

            // stessi centbin

            histEvSelCollisionAfterCuts_TReader = new TH1D("histEvSelCollisionAfterCuts_TReader", "", 1, 0, 1);

            double total_integral_AfterCuts_TReader = histCentFTOC_AfterCuts_TReader->Integral();
            histEvSelCollisionAfterCuts_TReader -> SetBinContent(1, total_integral_AfterCuts_TReader);

            histSelectedIntegrals_CentrFT0C_AfterCuts_TReader = new TH1D("histSelectedIntegrals_CentrFT0C_AfterCuts_TReader", "SelectedIntegral of CentrFT0C AfterCuts TReader;Centrality interval;Value", nCentBins, 0, nCentBins);
            histRatioIntegrals_CentrFT0C_AfterCuts_TReader = new TH1D("histRatioIntegrals_CentrFT0C_AfterCuts_TReader", "SelectedIntegral/TotalIntegral of CentrFT0C AfterCuts TReader;Centrality interval;Fraction", nCentBins, 0, nCentBins);
            
            for (int i = 0; i < nCentBins; i++) {
                double selected_integral_AfterCuts_TReader = histCentFTOC_AfterCuts_TReader->Integral(centMin[i]+1, centMax[i]);
                double ratio_integral_AfterCuts_TReader = selected_integral_AfterCuts_TReader / total_integral_AfterCuts_TReader;

                histSelectedIntegrals_CentrFT0C_AfterCuts_TReader->SetBinContent(i + 1, selected_integral_AfterCuts_TReader);
                histRatioIntegrals_CentrFT0C_AfterCuts_TReader->SetBinContent(i + 1, ratio_integral_AfterCuts_TReader);

                TString label_selected_integral_AfterCuts_TReader = Form("%d-%d", centMin[i], centMax[i]);
                TString label_ratio_integral_AfterCuts_TReader = Form("%d-%d", centMin[i], centMax[i]);
                histSelectedIntegrals_CentrFT0C_AfterCuts_TReader->GetXaxis()->SetBinLabel(i + 1, label_selected_integral_AfterCuts_TReader);
                histRatioIntegrals_CentrFT0C_AfterCuts_TReader->GetXaxis()->SetBinLabel(i + 1, label_ratio_integral_AfterCuts_TReader);

            }
        }
    }

    TFile *fOut = new TFile(fOutName.c_str(), "RECREATE");
    histEvSelCollisionBeforeCuts_TMaker -> Write();
    histEvSelCollisionBeforeCuts_TReader -> Write();
    histEvSelCollisionAfterCuts_TMaker -> Write();
    histEvSelCollisionAfterCuts_TReader -> Write();
    histBcSelCounterTVX -> Write();
    histBcSelCounterTVXafterBCcuts -> Write();
    histEvSelColCounterAll -> Write();
    histEvSelColCounterTVX -> Write();
    histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker->Write();
    histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker -> Write();
    histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker->Write();
    histRatioIntegrals_CentrFT0C_AfterCuts_TMaker -> Write();
    histSelectedIntegrals_CentrFT0C_BeforeCuts_TReader->Write();
    histRatioIntegrals_CentrFT0C_BeforeCuts_TReader -> Write();
    histSelectedIntegrals_CentrFT0C_AfterCuts_TReader->Write();
    histRatioIntegrals_CentrFT0C_AfterCuts_TReader -> Write();
    fOut -> Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void luminosity(string year = "2025", string period = "LHC25ae_pass2", string triggerMask = "minBias", string assocType = "std_assoc") {
    TFile *fIn = new TFile(Form("data/%s/pass2/Train_706901/normalization/%s_%s_%s_trigger_summary_Train_706901_Sel8.root", year.c_str(), period.c_str(), triggerMask.c_str(), assocType.c_str()), "READ");    

    TH1D *histEvSelCollisionBeforeCuts_TMaker = (TH1D*) fIn -> Get("histEvSelCollisionBeforeCuts_TMaker");
    TH1D *histEvSelCollisionBeforeCuts_TReader = (TH1D*) fIn -> Get("histEvSelCollisionBeforeCuts_TReader");
    TH1D *histEvSelCollisionAfterCuts_TMaker = (TH1D*) fIn -> Get("histEvSelCollisionAfterCuts_TMaker");
    TH1D *histEvSelCollisionAfterCuts_TReader = (TH1D*) fIn -> Get("histEvSelCollisionAfterCuts_TReader");
    TH1D *histBcSelCounterTVX = (TH1D*) fIn -> Get("histBcSelCounterTVX");
    TH1D *histBcSelCounterTVXafterBCcuts = (TH1D*) fIn -> Get("histBcSelCounterTVXafterBCcuts");
    TH1D *histEvSelColCounterAll = (TH1D*) fIn -> Get("histEvSelColCounterAll");
    TH1D *histEvSelColCounterTVX = (TH1D*) fIn -> Get("histEvSelColCounterTVX");
    TH1D *histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker = (TH1D*) fIn -> Get("histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker"); // same as 'histEvSelCollisionBeforeCuts_TMaker' but as a function of centr (The bin '0-100%' of 'histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker' should have exactly the same value of 'histEvSelCollisionBeforeCuts_TMaker')
    TH1D *histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker = (TH1D*) fIn -> Get("histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker");
    TH1D *histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker = (TH1D*) fIn -> Get("histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker"); // same as 'histEvSelCollisionAfterCuts_TMaker' but as a function of centr (The bin '0-100%' of 'histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker' should have exactly the same value of 'histEvSelCollisionAfterCuts_TMaker')
    TH1D *histRatioIntegrals_CentrFT0C_AfterCuts_TMaker = (TH1D*) fIn -> Get("histRatioIntegrals_CentrFT0C_AfterCuts_TMaker");
    TH1D *histSelectedIntegrals_CentrFT0C_BeforeCuts_TReader = (TH1D*) fIn -> Get("histSelectedIntegrals_CentrFT0C_BeforeCuts_TReader"); // same as 'histEvSelCollisionBeforeCuts_TReader' but as a function of centr (The bin '0-100%' of 'histSelectedIntegrals_CentrFT0C_BeforeCuts_TReader' should have exactly the same value of 'histEvSelCollisionBeforeCuts_TReader')
    TH1D *histRatioIntegrals_CentrFT0C_BeforeCuts_TReader = (TH1D*) fIn -> Get("histRatioIntegrals_CentrFT0C_BeforeCuts_TReader");
    TH1D *histSelectedIntegrals_CentrFT0C_AfterCuts_TReader = (TH1D*) fIn -> Get("histSelectedIntegrals_CentrFT0C_AfterCuts_TReader"); // same as 'histEvSelCollisionAfterCuts_TReader' but as a function of centr (The bin '0-100%' of 'histSelectedIntegrals_CentrFT0C_AfterCuts_TReader' should have exactly the same value of 'histEvSelCollisionAfterCuts_TReader')
    TH1D *histRatioIntegrals_CentrFT0C_AfterCuts_TReader = (TH1D*) fIn -> Get("histRatioIntegrals_CentrFT0C_AfterCuts_TReader");

    double evSelCollisionBeforeCuts_TMaker = histEvSelCollisionBeforeCuts_TMaker -> Integral();
    double evSelCollisionBeforeCuts_TReader = histEvSelCollisionBeforeCuts_TReader -> Integral();
    double evSelCollisionAfterCuts_TMaker = histEvSelCollisionAfterCuts_TMaker -> Integral();
    double evSelCollisionAfterCuts_TReader = histEvSelCollisionAfterCuts_TReader -> Integral();
    double bcCounterTVX = histBcSelCounterTVX -> Integral();
    double bcCounterTVXafterBCcuts = histBcSelCounterTVXafterBCcuts -> Integral();
    double evColCounterAll = histEvSelColCounterAll -> Integral();
    double evColCounterTVX = histEvSelColCounterTVX -> Integral();

    std::cout << "****************************************************" << std::endl;
    std::cout << "evSelCollisionBeforeCuts_TMaker = " << static_cast<long long>(evSelCollisionBeforeCuts_TMaker) << std::endl;
    std::cout << "evSelCollisionBeforeCuts_TReader= " << static_cast<long long>(evSelCollisionBeforeCuts_TReader) << std::endl;
    std::cout << "evSelCollisionAfterCuts_TMaker = " << static_cast<long long>(evSelCollisionAfterCuts_TMaker) << std::endl;
    std::cout << "evSelCollisionAfterCuts_TReader= " << static_cast<long long>(evSelCollisionAfterCuts_TReader) << std::endl;
    std::cout << "bcCounterTVX = " << static_cast<long long>(bcCounterTVX) << std::endl;
    std::cout << "bcCounterTVXafterBCcuts = " << static_cast<long long>(bcCounterTVXafterBCcuts) << std::endl;
    std::cout << "evColCounterAll = " << static_cast<long long>(evColCounterAll) << std::endl;
    std::cout << "evColCounterTVX = " << static_cast<long long>(evColCounterTVX) << std::endl;

    std::cout << "SelectedIntegral of CentrFT0C BeforeCuts TMaker(per bins):" << std::endl;
    int nCentrBins = histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker->GetNbinsX();
    for (int i = 0; i < nCentrBins; i++) {
        TString label = histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker->GetXaxis()->GetBinLabel(i+1);
        double value = histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker->GetBinContent(i+1);
        std::cout << label << " : " << static_cast<long long>(value) << std::endl;
    }
    std::cout << "SelectedIntegral/TotalIntegral of CentrFT0C BeforeCuts TMaker:" << std::endl; 
    for (int i = 0; i < nCentrBins; i++) {
        TString label = histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker->GetXaxis()->GetBinLabel(i+1);
        double value = histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker->GetBinContent(i+1);
        std::cout << label << " : " << value << std::endl;
    }
    std::cout << "SelectedIntegral of CentrFT0C AfterCuts TMaker (per bins):" << std::endl;
    for (int i = 0; i < nCentrBins; i++) {
        TString label = histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker->GetXaxis()->GetBinLabel(i+1);
        double value = histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker->GetBinContent(i+1);
        std::cout << label << " : " << static_cast<long long>(value) << std::endl;
    }
    std::cout << "SelectedIntegral/TotalIntegral of CentrFT0C AfterCuts TMaker:" << std::endl; 
    for (int i = 0; i < nCentrBins; i++) {
        TString label = histRatioIntegrals_CentrFT0C_AfterCuts_TMaker->GetXaxis()->GetBinLabel(i+1);
        double value = histRatioIntegrals_CentrFT0C_AfterCuts_TMaker->GetBinContent(i+1);
        std::cout << label << " : " << value << std::endl;
    }

    std::cout << "SelectedIntegral of CentrFT0C BeforeCuts TReader(per bins):" << std::endl;
    for (int i = 0; i < nCentrBins; i++) {
        TString label = histSelectedIntegrals_CentrFT0C_BeforeCuts_TReader->GetXaxis()->GetBinLabel(i+1);
        double value = histSelectedIntegrals_CentrFT0C_BeforeCuts_TReader->GetBinContent(i+1);
        std::cout << label << " : " << static_cast<long long>(value) << std::endl;
    }
    std::cout << "SelectedIntegral/TotalIntegral of CentrFT0C BeforeCuts TReader:" << std::endl; 
    for (int i = 0; i < nCentrBins; i++) {
        TString label = histRatioIntegrals_CentrFT0C_BeforeCuts_TReader->GetXaxis()->GetBinLabel(i+1);
        double value = histRatioIntegrals_CentrFT0C_BeforeCuts_TReader->GetBinContent(i+1);
        std::cout << label << " : " << value << std::endl;
    }
    std::cout << "SelectedIntegral of CentrFT0C AfterCuts TReader (per bins):" << std::endl;
    for (int i = 0; i < nCentrBins; i++) {
        TString label = histSelectedIntegrals_CentrFT0C_AfterCuts_TReader->GetXaxis()->GetBinLabel(i+1);
        double value = histSelectedIntegrals_CentrFT0C_AfterCuts_TReader->GetBinContent(i+1);
        std::cout << label << " : " << static_cast<long long>(value) << std::endl;
    }
    std::cout << "SelectedIntegral/TotalIntegral of CentrFT0C AfterCuts TReader:" << std::endl; 
    for (int i = 0; i < nCentrBins; i++) {
        TString label = histRatioIntegrals_CentrFT0C_AfterCuts_TReader->GetXaxis()->GetBinLabel(i+1);
        double value = histRatioIntegrals_CentrFT0C_AfterCuts_TReader->GetBinContent(i+1);
        std::cout << label << " : " << value << std::endl;
    }
    std::cout << "****************************************************" << std::endl;

    int nSummaryBins = 11; 
    TH1D *histLumiSummary = new TH1D("histLumiSummary", "", nSummaryBins, 0, nSummaryBins);
    histLumiSummary->GetXaxis()->SetBinLabel(1, "evSelCollisionBeforeCuts_TMaker");
    histLumiSummary->GetXaxis()->SetBinLabel(2, "evSelCollisionBeforeCuts_TReader");
    histLumiSummary->GetXaxis()->SetBinLabel(3, "evSelCollisionAfterCuts_TMaker");
    histLumiSummary->GetXaxis()->SetBinLabel(4, "evSelCollisionAfterCuts_TReader");
    histLumiSummary->GetXaxis()->SetBinLabel(5, "bcCounterTVX");
    histLumiSummary->GetXaxis()->SetBinLabel(6, "bcCounterTVXafterBCcuts");
    histLumiSummary->GetXaxis()->SetBinLabel(7, "bcCounterEfficiency");
    histLumiSummary->GetXaxis()->SetBinLabel(8, "evColCounterAll");
    histLumiSummary->GetXaxis()->SetBinLabel(9, "evColCounterTVX");
    histLumiSummary->GetXaxis()->SetBinLabel(10, "nEvtsBcSel");
    histLumiSummary->GetXaxis()->SetBinLabel(11, "luminosityBcSel");

    int nSummaryBins_nEvtsBcSelCentr = nCentrBins; 

    //------------------------------

    TH1D *hist_CentrFT0C_N_AfterCuts_TReader_Divided_N_BeforeCuts_TMaker = new TH1D("hist_CentrFT0C_N_AfterCuts_TReader_Divided_N_BeforeCuts_TMaker", "hist_CentrFT0C_N_AfterCuts_TReader_Divided_N_BeforeCuts_TMaker", nSummaryBins_nEvtsBcSelCentr, 0, nSummaryBins_nEvtsBcSelCentr);;
    for (int i = 0; i < nCentrBins; i++) {
        TString originalLabel = histSelectedIntegrals_CentrFT0C_BeforeCuts_TMaker->GetXaxis()->GetBinLabel(i+1);
        hist_CentrFT0C_N_AfterCuts_TReader_Divided_N_BeforeCuts_TMaker->GetXaxis()->SetBinLabel(1 + i, originalLabel);    
    }

    // ----------------------------

    TH1D *histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker = new TH1D("histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker", "histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker", nSummaryBins_nEvtsBcSelCentr, 0, nSummaryBins_nEvtsBcSelCentr);
    TH1D *histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp = new TH1D("histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp", "histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp", nSummaryBins_nEvtsBcSelCentr, 0, nSummaryBins_nEvtsBcSelCentr);
    for (int i = 0; i < nCentrBins; i++) {
        TString originalLabel = histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker->GetXaxis()->GetBinLabel(i+1);
        histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker->GetXaxis()->SetBinLabel(1 + i, originalLabel);        
        histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp->GetXaxis()->SetBinLabel(1 + i, originalLabel);  
    }


    std::vector<double> nEvtsBcSelCentr_AfterCuts_TMaker(nCentrBins, -999);
    std::vector<double> nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp(nCentrBins, -999);

    double bcCounterEfficiency = bcCounterTVXafterBCcuts / bcCounterTVX;
    double nEvtsBcSel = -999;
    double luminosityBcSel = -999;

    double pileup = 1.045; // Value for 0-100% - Da capire vs centrality
    double Eff_Train = (evSelCollisionAfterCuts_TMaker/evSelCollisionBeforeCuts_TReader);
    if (triggerMask.find("minBias") != std::string::npos) {
        nEvtsBcSel = pileup * bcCounterTVX * (evSelCollisionAfterCuts_TReader/evColCounterTVX) * Eff_Train; // N vs pT - ALWAYS add the PileUp Correction vs pT!

        for (int i = 0; i < nCentrBins; i++) { // N vs centr

            nEvtsBcSelCentr_AfterCuts_TMaker[i] = pileup * bcCounterTVX * (histSelectedIntegrals_CentrFT0C_AfterCuts_TReader->GetBinContent(i+1)/evColCounterTVX) * Eff_Train; // Used the pileup correction flat in all the centrality bins
            nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp[i] = bcCounterTVX * (histSelectedIntegrals_CentrFT0C_AfterCuts_TReader->GetBinContent(i+1)/evColCounterTVX) * Eff_Train; // Not Used the pileup correction -> To understand the pileup correction vs centrality bins
        
            // --------------------

            double value_TR_Divided_TM = (histSelectedIntegrals_CentrFT0C_AfterCuts_TReader->GetBinContent(i+1)/evColCounterTVX);
            hist_CentrFT0C_N_AfterCuts_TReader_Divided_N_BeforeCuts_TMaker->SetBinContent(i + 1, value_TR_Divided_TM);
        }

        luminosityBcSel = nEvtsBcSel / xSecTVX;

        std::cout << "nEvtsBcSel            = " << nEvtsBcSel << std::endl; // Da inserire ancora efficienza TVX

        std::cout << "nEvtsBcSelCentr_AfterCuts_TMaker = " << std::endl;
        for (size_t i = 0; i < nEvtsBcSelCentr_AfterCuts_TMaker.size(); i++) {
            TString label = histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker->GetXaxis()->GetBinLabel(i+1);
            std::cout << label << " : " << static_cast<long long>(nEvtsBcSelCentr_AfterCuts_TMaker[i]) << std::endl;
        }
        std::cout << "" << std::endl;
        std::cout << "nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp = " << std::endl;
        for (size_t i = 0; i < nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp.size(); i++) {
            TString label = histSelectedIntegrals_CentrFT0C_AfterCuts_TMaker->GetXaxis()->GetBinLabel(i+1);
            std::cout << label << " : " << static_cast<long long>(nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp[i]) << std::endl;
        }
        std::cout << "" << std::endl;
        std::cout << "luminosityBcSel = " << luminosityBcSel << std::endl;
        std::cout << std::endl;
    } else {

    }

    histLumiSummary -> SetBinContent(1, evSelCollisionBeforeCuts_TMaker);
    histLumiSummary -> SetBinContent(2, evSelCollisionBeforeCuts_TReader);
    histLumiSummary -> SetBinContent(3, evSelCollisionAfterCuts_TMaker);
    histLumiSummary -> SetBinContent(4, evSelCollisionAfterCuts_TReader);
    histLumiSummary -> SetBinContent(5, bcCounterTVX);
    histLumiSummary -> SetBinContent(6, bcCounterTVXafterBCcuts);
    histLumiSummary -> SetBinContent(7, bcCounterEfficiency);
    histLumiSummary -> SetBinContent(8, evColCounterAll);
    histLumiSummary -> SetBinContent(9, evColCounterTVX);
    histLumiSummary -> SetBinContent(10, nEvtsBcSel);
    histLumiSummary -> SetBinContent(11, luminosityBcSel);

    for (int i = 0; i < nCentrBins; i++) {
        histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker->SetBinContent(1 + i, nEvtsBcSelCentr_AfterCuts_TMaker[i]);
        histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp->SetBinContent(1 + i, nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp[i]);
    }

    TFile *fOut = new TFile(Form("data/%s/pass2/Train_706901/normalization/%s_%s_%s_luminosity_Train_706901_Sel8.root", year.c_str(), period.c_str(), triggerMask.c_str(), assocType.c_str()), "RECREATE");
    histLumiSummary -> Write(); 
    hist_CentrFT0C_N_AfterCuts_TReader_Divided_N_BeforeCuts_TMaker -> Write();
    histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker -> Write();
    histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp -> Write();
    fOut -> Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_normalization(string fInName = "LHC25ae_pass2_minBias_std_assoc_luminosity_Train_706901_Sel8.root", string fOutName = "luminosity_jpsi_LHC25ae_pass2_minBias_Train_706901_Sel8.root") {
    TFile *fInLumiMinBias2024StdAssoc = TFile::Open(Form("data/2025/pass2/Train_706901/normalization/%s", fInName.c_str()));
    TH1D *histLumiMinBias2024StdAssoc = (TH1D*) fInLumiMinBias2024StdAssoc -> Get("histLumiSummary");

    TH1D *hist_CentrFT0C_N_AfterCuts_TReader_Divided_N_BeforeCuts_TMaker = (TH1D*) fInLumiMinBias2024StdAssoc -> Get("hist_CentrFT0C_N_AfterCuts_TReader_Divided_N_BeforeCuts_TMaker");

    // Corrected one
    TH1D *histLumiMinBias2024StdAssoc_nEvtsBcSelCentr_AfterCuts_TMaker = (TH1D*) fInLumiMinBias2024StdAssoc -> Get("histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker");
    TH1D *histLumiMinBias2024StdAssoc_nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp = (TH1D*) fInLumiMinBias2024StdAssoc -> Get("histLumiSummary_nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp");

    double nEvtsMinBias2024StdAssoc = histLumiMinBias2024StdAssoc -> GetBinContent(histLumiMinBias2024StdAssoc -> GetXaxis() -> FindBin("nEvtsBcSel")); 
    double lumiMinBias2024StdAssoc = histLumiMinBias2024StdAssoc -> GetBinContent(histLumiMinBias2024StdAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));

    //const char* centLabels[] = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-90%", "90-100%","0-90%", "0-100%"};
    const char* centLabels[] = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%","0-90%", "0-100%"};
    int nCentrBins = sizeof(centLabels)/sizeof(char*);

    std::vector<double> nEvtsBcSelCentr_AfterCuts_TMaker(nCentrBins);
    std::vector<double> nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp(nCentrBins);
    for (int i = 0; i < nCentrBins; i++) {
        nEvtsBcSelCentr_AfterCuts_TMaker[i] = histLumiMinBias2024StdAssoc_nEvtsBcSelCentr_AfterCuts_TMaker->GetBinContent(i+1);
        nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp[i] = histLumiMinBias2024StdAssoc_nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp->GetBinContent(i+1);
    }

    TH1D *histLuminosityMinBias2024StdAssoc = new TH1D("histLumi", "; ; Luminosity (pb-1)", 1, 0, 1);
    histLuminosityMinBias2024StdAssoc -> SetBinContent(1, lumiMinBias2024StdAssoc);

    TH1D *histNeventMinBias2024StdAssoc = new TH1D("histNevMinBias", "; ; N. Events min bias", 1, 0, 1);
    histNeventMinBias2024StdAssoc -> SetBinContent(1, nEvtsMinBias2024StdAssoc);

    TH1D *histnEvtsBcSelCentr_AfterCuts_TMaker = new TH1D("histnEvtsBcSelCentr_AfterCuts_TMaker", "; ; histnEvtsBcSelCentr_AfterCuts_TMaker (N. Events min bias for a given centrality After Cuts TMaker)",nCentrBins, 0, nCentrBins);
    TH1D *histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp = new TH1D("histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp", "; ; histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp (N. Events min bias for a given centrality After Cuts TMaker NoPileUp)",nCentrBins, 0, nCentrBins);
    for (int i = 0; i < nCentrBins; i++) {
        histnEvtsBcSelCentr_AfterCuts_TMaker->SetBinContent(i+1, nEvtsBcSelCentr_AfterCuts_TMaker[i]);
        histnEvtsBcSelCentr_AfterCuts_TMaker->GetXaxis()->SetBinLabel(i+1, centLabels[i]);
        histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp->SetBinContent(i+1, nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp[i]);
        histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp->GetXaxis()->SetBinLabel(i+1, centLabels[i]);
    }

    TFile *fOutLumiMinBias2024StdAssoc = new TFile(Form("data/2025/pass2/Train_706901/normalization/%s", fOutName.c_str()), "RECREATE");
    histLuminosityMinBias2024StdAssoc -> Write();
    histNeventMinBias2024StdAssoc -> Write();
    hist_CentrFT0C_N_AfterCuts_TReader_Divided_N_BeforeCuts_TMaker -> Write();
    histnEvtsBcSelCentr_AfterCuts_TMaker -> Write();
    histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp -> Write();
    fOutLumiMinBias2024StdAssoc -> Close();

    std::cout << "---------------------------" << std::endl;
    std::cout << "[Min. Bias, std. assoc.]" << std::endl;
    std::cout << "nEvtsMinBias2024StdAssoc         = " << static_cast<long long>(nEvtsMinBias2024StdAssoc) << std::endl;

    for (int i = 0; i < nCentrBins; i++) {
        std::cout << "N. evts (nEvtsBcSelCentr_AfterCuts_TMaker) " << centLabels[i] << " = " << static_cast<long long>(nEvtsBcSelCentr_AfterCuts_TMaker[i]) << std::endl;
    }
    for (int i = 0; i < nCentrBins; i++) {
        std::cout << "N. evts (nEvtsBcSelCentr_AfterCuts_TMaker NoPileUp) " << centLabels[i] << " = " << static_cast<long long>(nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp[i]) << std::endl;
    }
    std::cout << "Luminosity        = " << lumiMinBias2024StdAssoc << std::endl;
    std::cout << "---------------------------" << std::endl;


}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void RetrieveTriggerInfo(TString dirName = "path/to/file", bool fromAlien = true, string triggerMask = "fDiMuon", double counters[10] = 0) {

    TString fInName;
    double collisionsBeforeFiltering = 0;
    double collisionsBeforeCuts_TMaker = 0;
    double collisionsAfterCuts_TMaker = 0;
    double infoTVX = 0;
    double infoScalTrig = 0;
    double infoSelTrig = 0;
    double selTOI = 0;
    double analysedTriggers = 0;
    double counterTVX = 0;
    double counterTVXafterBCcuts = 0;
    double colCounterTVX = 0;

    if (fromAlien) {
        TGrid::Connect("alien://");
        if (!dirName.Contains("alien://")) {
            dirName = "alien://" + dirName + "/AOD";
        }
    }

    for (int iDir = 1;iDir < 7000;iDir++) {
        fInName.Clear();
        if (iDir < 10) {fInName = dirName + Form("/00%i/AnalysisResults.root", iDir);}
        if (iDir >= 10 && iDir < 100) {fInName = dirName + Form("/0%i/AnalysisResults.root", iDir);}
        if (iDir >= 100) {fInName = dirName + Form("/%i/AnalysisResults.root", iDir);}

        TFile *fIn = TFile::Open(fInName.Data());
        if (!fIn || fIn -> IsZombie()) {
            std::cout << "ERROR! THE FILE DOES NOT EXIST!" << std::endl;
            break;
        } else {
            std::cout << "THE FILE EXISTS!" << std::endl;
        }

        for (auto dirKey : *fIn -> GetListOfKeys()) {
            if (TString(dirKey -> GetName()).Contains("table-maker")) {
                TList *list = (TList*) fIn -> Get("table-maker/Statistics");
                TH2D *histZorroInfo = (TH2D*) list -> FindObject("ZorroInfo");
                TH2D *histZorroSel = (TH2D*) list -> FindObject("ZorroSel");
                TH2D *histEventStats = (TH2D*) list -> FindObject("EventStats");
    
                //std::cout << "inspected TVX: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin("inspected TVX")) << std::endl;
                //std::cout << triggerMask << " scalers: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str()))) << std::endl;
                //std::cout << triggerMask << " fSingleMuLow selections: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str()))) << std::endl;
                //std::cout << triggerMask << " fSingleMuLow scalers (SEL): " << histZorroSel -> GetBinContent(1, histZorroSel -> GetYaxis() -> FindBin(Form("%s", triggerMask.c_str()))) << std::endl;

                infoTVX += histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin("inspected TVX"));
                infoScalTrig += histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str())));
                infoSelTrig += histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str())));
                selTOI += histZorroSel -> GetBinContent(1, histZorroSel -> GetYaxis() -> FindBin(Form("%s", triggerMask.c_str())));
                analysedTriggers += histZorroSel -> GetBinContent(1, histZorroSel -> GetYaxis() -> FindBin(Form("%s AnalysedTriggers", triggerMask.c_str())));

                collisionsBeforeFiltering += histEventStats -> GetBinContent(2, histEventStats -> GetYaxis() -> FindBin("Total"));
                collisionsBeforeCuts_TMaker += histEventStats -> GetBinContent(3, histEventStats -> GetYaxis() -> FindBin("Total"));
                collisionsAfterCuts_TMaker += histEventStats -> GetBinContent(4, histEventStats -> GetYaxis() -> FindBin("Total"));
            }
            if (TString(dirKey -> GetName()).Contains("bc-selection-task")) {
                TH1D *histCounterTVX = (TH1D*) fIn -> Get("bc-selection-task/hCounterTVX");
                TH1D *histCounterTVXafterBCcuts = (TH1D*) fIn -> Get("bc-selection-task/hCounterTVXafterBCcuts");

                counterTVX += histCounterTVX -> GetBinContent(1);
                counterTVXafterBCcuts += histCounterTVXafterBCcuts -> GetBinContent(1);
            }

            if (TString(dirKey -> GetName()).Contains("event-selection-task")) {
                TH1D *histColCounterTVX = (TH1D*) fIn -> Get("event-selection-task/hColCounterTVX");

                colCounterTVX += histColCounterTVX -> GetBinContent(1);
            }

            if (TString(dirKey -> GetName()).Contains("eventselection-run3")) {
                TH1D *histColCounterTVX = (TH1D*) fIn -> Get("eventselection-run3/eventselection/hColCounterTVX");
                TH1D *histCounterTVX = (TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVX");
                TH1D *histCounterTVXafterBCcuts = (TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVXafterBCcuts");

                colCounterTVX += histColCounterTVX -> GetBinContent(1);
                counterTVX += histCounterTVX -> GetBinContent(1);
                counterTVXafterBCcuts += histCounterTVXafterBCcuts -> GetBinContent(1);
            }
        }

        counters[0] = infoTVX;
        counters[1] = infoScalTrig;
        counters[2] = infoSelTrig;
        counters[3] = selTOI;
        counters[4] = analysedTriggers;
        counters[5] = collisionsBeforeFiltering;
        counters[6] = collisionsBeforeCuts_TMaker;
        counters[7] = collisionsAfterCuts_TMaker;
        counters[8] = counterTVX;
        counters[9] = counterTVXafterBCcuts;
        counters[10] = colCounterTVX;

        fIn -> Close();
    }
}