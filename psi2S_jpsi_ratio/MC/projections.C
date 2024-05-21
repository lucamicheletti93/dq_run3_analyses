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

void LoadStyle();
void SetLegend(TLegend *);

TH1D* projectHistogram(THnSparseD *, double , double , double , double );

void projections() {
    const int nPtBins = 5;
    double minPtBins[] = {0, 1, 2, 3, 5};
    double maxPtBins[] = {1, 2, 3, 5, 10};
    const int nRapBins = 1;
    double minRapBins[] = {2.5};
    double maxRapBins[] = {4};

    /*
    string fInName = "AnalysisResults_dq_efficiency";
    string cuts[] = {
        "PairsMuonSEPM_matchedMchMid_mumuFromJpsi",
        "PairsMuonSEPM_muonLowPt10SigmaPDCA_mumuFromJpsi",
        "PairsMuonSEPM_muonLowPt210SigmaPDCA_mumuFromJpsi",
        "PairsMuonSEPM_muonLowPt510SigmaPDCA_mumuFromJpsi",
        "PairsMuonSEPM_matchedMchMid_mumuFromPsi2S",
        "PairsMuonSEPM_muonLowPt10SigmaPDCA_mumuFromPsi2S",
        "PairsMuonSEPM_muonLowPt210SigmaPDCA_mumuFromPsi2S",
        "PairsMuonSEPM_muonLowPt510SigmaPDCA_mumuFromPsi2S"
    };
    */
    string fInName = "AnalysisResults_table_reader";
    string cuts[] = {
        "PairsMuonSEPM_matchedMchMid",
        "PairsMuonSEPM_muonLowPt10SigmaPDCA",
        "PairsMuonSEPM_muonLowPt210SigmaPDCA",
        "PairsMuonSEPM_muonLowPt510SigmaPDCA"
    };

    string histName = "Mass_Pt_Rapidity";

    TFile *fIn = new TFile(Form("/Users/lucamicheletti/alice/local_train_test_mc/LHC24e5/%s.root", fInName.c_str()), "READ");
    TList *list1 = (TList*) fIn -> Get("analysis-same-event-pairing/output");

    TFile *fOut = new TFile(Form("Histograms_%s.root", fInName.c_str()), "RECREATE");
    for (auto& cut : cuts) {
        TList *list = (TList*) list1 -> FindObject(cut.c_str());
        THnSparseD *histSparse = (THnSparseD*) list -> FindObject(histName.c_str());
        TH1D *histProjInt = (TH1D*) histSparse -> Projection(0, Form("Proj_%s", cut.c_str()));
        histProjInt -> Write(Form("Proj_%s", cut.c_str()));


        for (int iPt = 0;iPt < nPtBins;iPt++) {
            TH1D *histProjPt = (TH1D*) projectHistogram(histSparse, minPtBins[iPt], maxPtBins[iPt], minRapBins[0], maxRapBins[0]);
            histProjPt -> Write(Form("Proj_%s__Pt_%1.0f_%1.0f", cut.c_str(), minPtBins[iPt], maxPtBins[iPt]));
        }
    }
    fOut -> ls();
    fOut -> Close();

}

////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.04);
}

TH1D* projectHistogram(THnSparseD *histSparse, double minPtRange, double maxPtRange, double minRapRange, double maxRapRange) {
    double minPtBin = histSparse -> GetAxis(1) -> FindBin(minPtRange);
    double maxPtBin = histSparse -> GetAxis(1) -> FindBin(maxPtRange - 0.01);
    double minCentrBin = histSparse -> GetAxis(2) -> FindBin(minRapRange);
    double maxCentrBin = histSparse -> GetAxis(2) -> FindBin(maxRapRange - 0.01);

    Printf("minPtBin = %1.0f, maxPtBin = %1.0f, minCentrBin = %1.0f, maxCentrBin = %1.0f", minPtBin, maxPtBin, minCentrBin, maxCentrBin);
    histSparse -> GetAxis(1) -> SetRange(minPtBin, maxPtBin);
    histSparse -> GetAxis(2) -> SetRange(minCentrBin, maxCentrBin);

    TH1D *histProj = (TH1D*) histSparse -> Projection(0, Form("histProj_%1.0f_%1.0f__%1.0f_%1.0f", minPtRange, maxPtRange, minRapRange, maxRapRange));
    return histProj;
}