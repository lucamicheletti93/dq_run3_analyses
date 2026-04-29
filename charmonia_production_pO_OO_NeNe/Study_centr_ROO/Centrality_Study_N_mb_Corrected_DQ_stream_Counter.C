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

// per normalization_Correct_DQ_stream_Counter.C

void check_centrality(string fileName1 = "LHC25ae_pass2_minBias_std_assoc_trigger_summary_NTVX_Plus_CorrEvSelFirstCorrected_Plus_PileUp_in_Integrated_Plus_PileUp_and_Not_in_Centr__Train_706903_Centr_0_100_Sel8NoSameBunch.root", string fileName2 = "luminosity_jpsi_LHC25_LHC25ae_pass2_minBias_NTVX_Plus_CorrEvSelCorrected_Plus_PileUp_in_Integrated_Plus_PileUp_and_Not_in_Centr__Train_706903_Centr_0_100_Sel8NoSameBunch.root") {


    TFile *f1 = TFile::Open(Form("data/2025/pass2/Train_617522/%s", fileName1.c_str()));
    if (!f1 || f1->IsZombie()) {
        std::cout << "ERROR opening file" << std::endl;
        return;
    }

    TFile *f2 = TFile::Open(Form("data/2025/pass2/Train_617522/%s", fileName2.c_str()));
    if (!f2 || f2->IsZombie()) {
        std::cout << "ERROR opening file" << std::endl;
        return;
    }


    // ------------------
    // Ratio Distributions
    // ------------------
    TH1D *histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker = (TH1D*) f1->Get("histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker");
    TH1D *histRatioIntegrals_CentrFT0C_AfterCuts_TMaker = (TH1D*) f1->Get("histRatioIntegrals_CentrFT0C_AfterCuts_TMaker");
    TH1D *histRatioIntegrals_CentrFT0C_BeforeCuts_TReader = (TH1D*) f1->Get("histRatioIntegrals_CentrFT0C_BeforeCuts_TReader");
    TH1D *histRatioIntegrals_CentrFT0C_AfterCuts_TReader = (TH1D*) f1->Get("histRatioIntegrals_CentrFT0C_AfterCuts_TReader");

    const char* centLabels[] = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%","0-90%", "0-100%"};
    int nCentrBins = sizeof(centLabels)/sizeof(char*);

    std::vector<double> RatioIntegrals_CentrFT0C_BeforeCuts_TMaker(nCentrBins);
    std::vector<double> RatioIntegrals_CentrFT0C_AfterCuts_TMaker(nCentrBins);
    std::vector<double> RatioIntegrals_CentrFT0C_BeforeCuts_TReader(nCentrBins);
    std::vector<double> RatioIntegrals_CentrFT0C_AfterCuts_TReader(nCentrBins);

    for (int i = 0; i < nCentrBins; i++) {
        RatioIntegrals_CentrFT0C_BeforeCuts_TMaker[i] = histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker->GetBinContent(i+1);
        RatioIntegrals_CentrFT0C_AfterCuts_TMaker[i] = histRatioIntegrals_CentrFT0C_AfterCuts_TMaker->GetBinContent(i+1);
        RatioIntegrals_CentrFT0C_BeforeCuts_TReader[i] = histRatioIntegrals_CentrFT0C_BeforeCuts_TReader->GetBinContent(i+1);
        RatioIntegrals_CentrFT0C_AfterCuts_TReader[i] = histRatioIntegrals_CentrFT0C_AfterCuts_TReader->GetBinContent(i+1);
    }

    double sumAllCentr_RatioIntegrals_BeforeCuts_TMaker = 0.0;
    double sumAllCentr_RatioIntegrals_AfterCuts_TMaker = 0.0;
    double sumAllCentr_RatioIntegrals_BeforeCuts_TReader = 0.0;
    double sumAllCentr_RatioIntegrals_AfterCuts_TReader = 0.0;

    for (int i = 0; i < 10; i++) {   // solo 0-10 ... 90-100
        sumAllCentr_RatioIntegrals_BeforeCuts_TMaker += RatioIntegrals_CentrFT0C_BeforeCuts_TMaker[i];
        sumAllCentr_RatioIntegrals_AfterCuts_TMaker += RatioIntegrals_CentrFT0C_AfterCuts_TMaker[i];
        sumAllCentr_RatioIntegrals_BeforeCuts_TReader += RatioIntegrals_CentrFT0C_BeforeCuts_TReader[i];
        sumAllCentr_RatioIntegrals_AfterCuts_TReader += RatioIntegrals_CentrFT0C_AfterCuts_TReader[i];
    }

    double integrated_RatioIntegrals_0_100_AC_TMaker = histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker->GetBinContent(histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker->GetXaxis()->FindBin("0-100"));
    std::cout << "sumAllCentr_RatioIntegrals_BeforeCuts_TMaker 0-100% = " << sumAllCentr_RatioIntegrals_BeforeCuts_TMaker << ", 0-100 integrated: " << static_cast<long long>(integrated_RatioIntegrals_0_100_AC_TMaker) << std::endl;
    double integrated_RatioIntegrals_0_100_BC_TMaker = histRatioIntegrals_CentrFT0C_AfterCuts_TMaker->GetBinContent(histRatioIntegrals_CentrFT0C_AfterCuts_TMaker->GetXaxis()->FindBin("0-100"));
    std::cout << "sumAllCentr_RatioIntegrals_AfterCuts_TMaker 0-100% = " << sumAllCentr_RatioIntegrals_AfterCuts_TMaker << ", 0-100 integrated: " << static_cast<long long>(integrated_RatioIntegrals_0_100_BC_TMaker) << std::endl;

    double integrated_RatioIntegrals_0_100_AC_TReader = histRatioIntegrals_CentrFT0C_BeforeCuts_TReader->GetBinContent(histRatioIntegrals_CentrFT0C_BeforeCuts_TReader->GetXaxis()->FindBin("0-100"));
    std::cout << "sumAllCentr_RatioIntegrals_BeforeCuts_TReader 0-100% = " << sumAllCentr_RatioIntegrals_BeforeCuts_TReader << ", 0-100 integrated: " << static_cast<long long>(integrated_RatioIntegrals_0_100_AC_TReader) << std::endl;
    double integrated_RatioIntegrals_0_100_BC_TReader = histRatioIntegrals_CentrFT0C_AfterCuts_TReader->GetBinContent(histRatioIntegrals_CentrFT0C_AfterCuts_TReader->GetXaxis()->FindBin("0-100"));
    std::cout << "sumAllCentr_RatioIntegrals_AfterCuts_TReader 0-100% = " << sumAllCentr_RatioIntegrals_AfterCuts_TReader << ", 0-100 integrated: " << static_cast<long long>(integrated_RatioIntegrals_0_100_BC_TReader) << std::endl;

    TH1D* histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker_Until100 = new TH1D("histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker_Until100", histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker->GetTitle(), 10, 0.5, 10.5);
    TH1D* histRatioIntegrals_CentrFT0C_AfterCuts_TMaker_Until100 = new TH1D("histRatioIntegrals_CentrFT0C_AfterCuts_TMaker_Until100", histRatioIntegrals_CentrFT0C_AfterCuts_TMaker->GetTitle(), 10, 0.5, 10.5);
    TH1D* histRatioIntegrals_CentrFT0C_BeforeCuts_TReader_Until100 = new TH1D("histRatioIntegrals_CentrFT0C_BeforeCuts_TReader_Until100", histRatioIntegrals_CentrFT0C_BeforeCuts_TReader->GetTitle(), 10, 0.5, 10.5);
    TH1D* histRatioIntegrals_CentrFT0C_AfterCuts_TReader_Until100 = new TH1D("histRatioIntegrals_CentrFT0C_AfterCuts_TReader_Until100", histRatioIntegrals_CentrFT0C_AfterCuts_TReader->GetTitle(), 10, 0.5, 10.5);

    // Copy the first 10 bin (0-10, 10-20, ... , 90-100)
    for(int i=1; i<=10; i++){
        histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker_Until100->SetBinContent(i, histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker->GetBinContent(i));
        histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker_Until100->SetBinError(i, histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker->GetBinError(i));
        histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker_Until100->GetXaxis()->SetBinLabel(i, histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker->GetXaxis()->GetBinLabel(i));
    
        histRatioIntegrals_CentrFT0C_AfterCuts_TMaker_Until100->SetBinContent(i, histRatioIntegrals_CentrFT0C_AfterCuts_TMaker->GetBinContent(i));
        histRatioIntegrals_CentrFT0C_AfterCuts_TMaker_Until100->SetBinError(i, histRatioIntegrals_CentrFT0C_AfterCuts_TMaker->GetBinError(i));
        histRatioIntegrals_CentrFT0C_AfterCuts_TMaker_Until100->GetXaxis()->SetBinLabel(i, histRatioIntegrals_CentrFT0C_AfterCuts_TMaker->GetXaxis()->GetBinLabel(i));
    

        histRatioIntegrals_CentrFT0C_BeforeCuts_TReader_Until100->SetBinContent(i, histRatioIntegrals_CentrFT0C_BeforeCuts_TReader->GetBinContent(i));
        histRatioIntegrals_CentrFT0C_BeforeCuts_TReader_Until100->SetBinError(i, histRatioIntegrals_CentrFT0C_BeforeCuts_TReader->GetBinError(i));
        histRatioIntegrals_CentrFT0C_BeforeCuts_TReader_Until100->GetXaxis()->SetBinLabel(i, histRatioIntegrals_CentrFT0C_BeforeCuts_TReader->GetXaxis()->GetBinLabel(i));
    
        histRatioIntegrals_CentrFT0C_AfterCuts_TReader_Until100->SetBinContent(i, histRatioIntegrals_CentrFT0C_AfterCuts_TReader->GetBinContent(i));
        histRatioIntegrals_CentrFT0C_AfterCuts_TReader_Until100->SetBinError(i, histRatioIntegrals_CentrFT0C_AfterCuts_TReader->GetBinError(i));
        histRatioIntegrals_CentrFT0C_AfterCuts_TReader_Until100->GetXaxis()->SetBinLabel(i, histRatioIntegrals_CentrFT0C_AfterCuts_TReader->GetXaxis()->GetBinLabel(i));
    
    }
    
    histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker_Until100->SetDirectory(0);
    histRatioIntegrals_CentrFT0C_AfterCuts_TMaker_Until100->SetDirectory(0);
    histRatioIntegrals_CentrFT0C_BeforeCuts_TReader_Until100->SetDirectory(0);
    histRatioIntegrals_CentrFT0C_AfterCuts_TReader_Until100->SetDirectory(0);

    int nbins = histRatioIntegrals_CentrFT0C_BeforeCuts_TReader_Until100->GetNbinsX();
    TLine *line = new TLine(0.5, 0.1, nbins + 0.5, 0.1);
    line->SetLineStyle(2);   
    line->SetLineWidth(2);
    line->SetLineColor(kRed);

    TCanvas *c1 = new TCanvas("c1"," Ratio Centrality",800,600);
    c1->Divide(2,2); 
    c1->cd(1);
    double ymin = 0.09;
    double ymax = 0.13;   

    histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker_Until100->GetYaxis()->SetRangeUser(ymin, ymax);
    histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker_Until100->SetMarkerStyle(20);
    histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker_Until100->SetMarkerSize(1.2);
    histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker_Until100->SetLineColor(kBlack);
    histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker_Until100->SetStats(0);
    histRatioIntegrals_CentrFT0C_BeforeCuts_TMaker_Until100->Draw("HIST");
    line->Draw("SAME");

    c1->cd(2);
    histRatioIntegrals_CentrFT0C_AfterCuts_TMaker_Until100->GetYaxis()->SetRangeUser(ymin, ymax);
    histRatioIntegrals_CentrFT0C_AfterCuts_TMaker_Until100->SetMarkerStyle(20);
    histRatioIntegrals_CentrFT0C_AfterCuts_TMaker_Until100->SetMarkerSize(1.2);
    histRatioIntegrals_CentrFT0C_AfterCuts_TMaker_Until100->SetLineColor(kBlack);
    histRatioIntegrals_CentrFT0C_AfterCuts_TMaker_Until100->SetStats(0);
    histRatioIntegrals_CentrFT0C_AfterCuts_TMaker_Until100->Draw("HIST");
    line->Draw("SAME");

    c1->cd(3); 
    histRatioIntegrals_CentrFT0C_BeforeCuts_TReader_Until100->GetYaxis()->SetRangeUser(ymin, ymax);
    histRatioIntegrals_CentrFT0C_BeforeCuts_TReader_Until100->SetMarkerStyle(20);
    histRatioIntegrals_CentrFT0C_BeforeCuts_TReader_Until100->SetMarkerSize(1.2);
    histRatioIntegrals_CentrFT0C_BeforeCuts_TReader_Until100->SetLineColor(kBlack);
    histRatioIntegrals_CentrFT0C_BeforeCuts_TReader_Until100->SetStats(0);
    histRatioIntegrals_CentrFT0C_BeforeCuts_TReader_Until100->Draw("HIST");
    line->Draw("SAME");

    c1->cd(4);
    histRatioIntegrals_CentrFT0C_AfterCuts_TReader_Until100->GetYaxis()->SetRangeUser(ymin, ymax);
    histRatioIntegrals_CentrFT0C_AfterCuts_TReader_Until100->SetMarkerStyle(20);
    histRatioIntegrals_CentrFT0C_AfterCuts_TReader_Until100->SetMarkerSize(1.2);
    histRatioIntegrals_CentrFT0C_AfterCuts_TReader_Until100->SetLineColor(kBlack);
    histRatioIntegrals_CentrFT0C_AfterCuts_TReader_Until100->SetStats(0);
    histRatioIntegrals_CentrFT0C_AfterCuts_TReader_Until100->Draw("HIST");
    line->Draw("SAME");
    c1->Update();

    // ---------------
    // NevSel Values
    // ---------------
    TH1D *histnEvtsBcSelCentr_AfterCuts_TMaker = (TH1D*) f2->Get("histnEvtsBcSelCentr_AfterCuts_TMaker");
    TH1D *histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp = (TH1D*) f2->Get("histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp");

    std::vector<double> nEvtsBcSelCentr_AfterCuts_TMaker(nCentrBins);
    std::vector<double> nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp(nCentrBins);

    for (int i = 0; i < nCentrBins; i++) {

        nEvtsBcSelCentr_AfterCuts_TMaker[i] = histnEvtsBcSelCentr_AfterCuts_TMaker->GetBinContent(i+1);
        nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp[i] = histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp->GetBinContent(i+1);
    }

    double sumAllCentr_AfterCuts_TMaker_until90 = 0.0;
    double sumAllCentr_AfterCuts_TMaker_NoPileUp_until90 = 0.0;
    double sumAllCentr_AfterCuts_TMaker = 0.0;
    double sumAllCentr_AfterCuts_TMaker_NoPileUp = 0.0;

    for (int i = 0; i < 9; i++) {   // solo 0-10 ... 90-100
        sumAllCentr_AfterCuts_TMaker_until90 += nEvtsBcSelCentr_AfterCuts_TMaker[i];
        sumAllCentr_AfterCuts_TMaker_NoPileUp_until90 += nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp[i];
    }

    for (int i = 0; i < 10; i++) {   // solo 0-10 ... 90-100
        sumAllCentr_AfterCuts_TMaker += nEvtsBcSelCentr_AfterCuts_TMaker[i];
        sumAllCentr_AfterCuts_TMaker_NoPileUp += nEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp[i];
    }

    double integrated_0_90_TMaker = histnEvtsBcSelCentr_AfterCuts_TMaker->GetBinContent(histnEvtsBcSelCentr_AfterCuts_TMaker->GetXaxis()->FindBin("0-90%"));
    double integrated_0_90_TMaker_NoPileUp = histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp->GetBinContent(histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp->GetXaxis()->FindBin("0-90%"));
    std::cout << "sumAllCentr_AfterCuts_TMaker 0-90% = " << static_cast<long long>(sumAllCentr_AfterCuts_TMaker_until90)  << ", 0-90 integrated: " << static_cast<long long>(integrated_0_90_TMaker) << std::endl;
    std::cout << "sumAllCentr_AfterCuts_TMaker_NoPileUp 0-90% = " << static_cast<long long>(sumAllCentr_AfterCuts_TMaker_NoPileUp_until90)  << ", 0-90 integrated: " << static_cast<long long>(integrated_0_90_TMaker_NoPileUp) << std::endl;

    double integrated_0_100_TMaker = histnEvtsBcSelCentr_AfterCuts_TMaker->GetBinContent(histnEvtsBcSelCentr_AfterCuts_TMaker->GetXaxis()->FindBin("0-100%"));
    double integrated_0_100_TMaker_NoPileUp = histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp->GetBinContent(histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp->GetXaxis()->FindBin("0-100%"));
    std::cout << "sumAllCentr_AfterCuts_TMaker 0-100% = " << static_cast<long long>(sumAllCentr_AfterCuts_TMaker)  << ", 0-100 integrated: " << static_cast<long long>(integrated_0_100_TMaker) << std::endl;
    std::cout << "sumAllCentr_AfterCuts_TMaker_NoPileUp 0-100% = " << static_cast<long long>(sumAllCentr_AfterCuts_TMaker_NoPileUp)  << ", 0-100 integrated: " << static_cast<long long>(integrated_0_100_TMaker_NoPileUp) << std::endl;

    TH1D *histnEvtsBcSelCentr_AfterCuts_TMaker_Until100 = new TH1D("histnEvtsBcSelCentr_AfterCuts_TMaker_Until100", histnEvtsBcSelCentr_AfterCuts_TMaker->GetTitle(), 10, 0.5, 10.5);
    TH1D *histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp_Until100 = new TH1D("histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp_Until100", histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp->GetTitle(), 10, 0.5, 10.5);

    //Copy the first 10 bin (0-10, 10-20, ... , 90-100)
    for(int i=1; i<= 10; i++){
        histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->SetBinContent(i, histnEvtsBcSelCentr_AfterCuts_TMaker->GetBinContent(i));
        histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->SetBinError(i, histnEvtsBcSelCentr_AfterCuts_TMaker->GetBinError(i));
        histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->GetXaxis()->SetBinLabel(i, histnEvtsBcSelCentr_AfterCuts_TMaker->GetXaxis()->GetBinLabel(i));

        histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp_Until100->SetBinContent(i, histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp_Until100->GetBinContent(i));
        histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp_Until100->SetBinError(i, histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp_Until100->GetBinError(i));
        histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp_Until100->GetXaxis()->SetBinLabel(i, histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp_Until100->GetXaxis()->GetBinLabel(i));
    }
    
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->SetDirectory(0);
    histnEvtsBcSelCentr_AfterCuts_TMaker_NoPileUp_Until100->SetDirectory(0); 

    TCanvas *c3 = new TCanvas("c3","Centrality",800,600);
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->SetMarkerStyle(20);
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->SetMarkerSize(1.2);
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->SetLineColor(kBlack);
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->SetStats(0);
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->Draw("HIST");

    // Normalization 
    TH1D* histnEvtsBcSelCentr_AfterCuts_TMaker_Until100_norm = (TH1D*)histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->Clone("histnEvtsBcSelCentr_AfterCuts_TMaker_Until100_norm");
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100_norm->SetTitle("Centrality distribution (normalized)");

    double integral = histnEvtsBcSelCentr_AfterCuts_TMaker_Until100_norm->Integral();
    if (integral != 0) {
        histnEvtsBcSelCentr_AfterCuts_TMaker_Until100_norm->Scale(1.0 / integral);
    }

    TCanvas* c_norm = new TCanvas("c_norm","Centrality normalized",800,600);
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100_norm->SetMarkerStyle(20);
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100_norm->SetMarkerSize(1.2);
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100_norm->SetLineColor(kBlue);
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100_norm->SetStats(0);
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100_norm->GetYaxis()->SetTitle("Fraction of events");
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100_norm->Draw("HIST");
    line->Draw("SAME");

    c_norm->Update();





}
