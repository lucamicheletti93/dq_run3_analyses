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

void check_centrality(string fileName2 = "luminosity_jpsi_LHC25_LHC25ae_pass2_minBias_NTVX_Plus_CorrEvSelCorrected_Plus_PileUp_in_Integrated_Plus_PileUp_and_Not_in_Centr__Train_706903_Centr_0_100_Sel8NoSameBunch.root") {

    TFile *f2 = TFile::Open(Form("data/2025/pass2/Train_706901_706902_706903_610867/Train_706902/Base_CounterCollTVX/%s", fileName2.c_str()));
    if (!f2 || f2->IsZombie()) {
        std::cout << "ERROR opening file" << std::endl;
        return;
    }

    const char* centLabels[] = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-90%", "90-100%","0-90%", "0-100%"};
    int nCentrBins = sizeof(centLabels)/sizeof(char*);


    // ---------------
    // NevSel Values
    // ---------------
    TH1D *histnEvtsBcSelCentr_AfterCuts_TMaker = (TH1D*) f2->Get("histnEvtsBcSelCentr_AfterCuts_TMaker");

    std::vector<double> nEvtsBcSelCentr_AfterCuts_TMaker(nCentrBins);

    for (int i = 0; i < nCentrBins; i++) {
        nEvtsBcSelCentr_AfterCuts_TMaker[i] = histnEvtsBcSelCentr_AfterCuts_TMaker->GetBinContent(i+1);
    }

    double sumAllCentr_AfterCuts_TMaker_until90 = 0.0;
    double sumAllCentr_AfterCuts_TMaker = 0.0;

    for (int i = 0; i < 8; i++) {   // solo 0-10 ... 80-90
        sumAllCentr_AfterCuts_TMaker_until90 += nEvtsBcSelCentr_AfterCuts_TMaker[i];
    }

    for (int i = 0; i < 9; i++) {   // solo 0-10 ... 90-100
        sumAllCentr_AfterCuts_TMaker += nEvtsBcSelCentr_AfterCuts_TMaker[i];
    }

    double integrated_0_90_TMaker = histnEvtsBcSelCentr_AfterCuts_TMaker->GetBinContent(histnEvtsBcSelCentr_AfterCuts_TMaker->GetXaxis()->FindBin("0-90%"));
    std::cout << "sumAllCentr_AfterCuts_TMaker 0-90% = " << static_cast<long long>(sumAllCentr_AfterCuts_TMaker_until90)  << ", 0-90 integrated: " << static_cast<long long>(integrated_0_90_TMaker) << std::endl;

    double integrated_0_100_TMaker = histnEvtsBcSelCentr_AfterCuts_TMaker->GetBinContent(histnEvtsBcSelCentr_AfterCuts_TMaker->GetXaxis()->FindBin("0-100%"));
    std::cout << "sumAllCentr_AfterCuts_TMaker 0-100% = " << static_cast<long long>(sumAllCentr_AfterCuts_TMaker)  << ", 0-100 integrated: " << static_cast<long long>(integrated_0_100_TMaker) << std::endl;

    TH1D *histnEvtsBcSelCentr_AfterCuts_TMaker_Until100 = new TH1D("histnEvtsBcSelCentr_AfterCuts_TMaker_Until100", histnEvtsBcSelCentr_AfterCuts_TMaker->GetTitle(), 9, 0.5, 9.5);

    //Copy the first bins(0-10, 10-20, ... , 90-100)
    for(int i=1; i <= 9; i++){
        histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->SetBinContent(i, histnEvtsBcSelCentr_AfterCuts_TMaker->GetBinContent(i));
        histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->SetBinError(i, histnEvtsBcSelCentr_AfterCuts_TMaker->GetBinError(i));
        histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->GetXaxis()->SetBinLabel(i, histnEvtsBcSelCentr_AfterCuts_TMaker->GetXaxis()->GetBinLabel(i));
    }
    
    histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->SetDirectory(0);

    int nbins = histnEvtsBcSelCentr_AfterCuts_TMaker_Until100->GetNbinsX();
    TLine *line = new TLine(0.5, 0.1, nbins + 0.5, 0.1);
    line->SetLineStyle(2);   
    line->SetLineWidth(2);
    line->SetLineColor(kRed);
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
    std::cout << "integral histnEvtsBcSelCentr_AfterCuts_TMaker_Until100_norm = " << static_cast<long long>(integral)  << std::endl;

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
