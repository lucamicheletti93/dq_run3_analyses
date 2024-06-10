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

TH1D* ProjectTHnSparse(THnSparseD *, double , double , double , double );
TH1D* ProjectTH2(TH2D *, double , double );

void projections() {
    // MC settings
    //string dataset = "MC";
    //string productionName = "LHC24e4";
    //string associationType = "std_association";
    // Data settings
    string dataset = "2022";
    string productionName = "LHC22o_pass6_minBias";
    string associationType = "time_association";

    string histName = "Mass_Pt_Rapidity";
    //string histName = "Mass_Pt";

    // pt differential
    const int nPtBins1 = 8;
    double minPtBins1[] = {0, 1, 2, 3, 4, 5, 7, 10};
    double maxPtBins1[] = {1, 2, 3, 4, 5, 7, 10, 20};
    const int nRapBins1 = 1;
    double minRapBins1[] = {2.5};
    double maxRapBins1[] = {4};

    // rapidity differential
    const int nPtBins2 = 1;
    double minPtBins2[] = {0};
    double maxPtBins2[] = {20};
    const int nRapBins2 = 6;
    double minRapBins2[] = {2.50, 2.75, 3.00, 3.25, 3.50, 3.75};
    double maxRapBins2[] = {2.75, 3.00, 3.25, 3.50, 3.75, 4.00};

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

    TFile *fIn = new TFile(Form("/Users/lucamicheletti/cernbox/JPSI/Run3/%s/%s/%s_%s.root", dataset.c_str(), productionName.c_str(), fInName.c_str(), associationType.c_str()), "READ");
    fIn -> ls();
    TList *list1 = (TList*) fIn -> Get("analysis-same-event-pairing/output");

    TFile *fOut = new TFile(Form("/Users/lucamicheletti/cernbox/JPSI/Run3/%s/%s/Histograms_%s_%s.root", dataset.c_str(), productionName.c_str(), fInName.c_str(), associationType.c_str()), "RECREATE");
    for (auto& cut : cuts) {
        TList *list = (TList*) list1 -> FindObject(cut.c_str());

        if (histName == "Mass_Pt_Rapidity") {
            THnSparseD *histSparse = (THnSparseD*) list -> FindObject(histName.c_str());
            TH1D *histProjInt = (TH1D*) histSparse -> Projection(0, Form("Proj_%s", cut.c_str()));
            histProjInt -> Write(Form("Proj_%s", cut.c_str()));


            for (int iPt = 0;iPt < nPtBins1;iPt++) {
                TH1D *histProjPt = (TH1D*) ProjectTHnSparse(histSparse, minPtBins1[iPt], maxPtBins1[iPt], minRapBins1[0], maxRapBins1[0]);
                histProjPt -> Write(Form("Proj_%s__Pt_%1.0f_%1.0f", cut.c_str(), minPtBins1[iPt], maxPtBins1[iPt]));
            }

            for (int iRap = 0;iRap < nRapBins2;iRap++) {
                TH1D *histProjRap = (TH1D*) ProjectTHnSparse(histSparse, minPtBins2[0], maxPtBins2[0], minRapBins2[iRap], maxRapBins2[iRap]);
                histProjRap -> Write(Form("Proj_%s__Rap_%3.2f_%3.2f", cut.c_str(), minRapBins2[iRap], maxRapBins2[iRap]));
            }
        }

        if (histName == "Mass_Pt") {
            TH2D *hist2D = (TH2D*) list -> FindObject(histName.c_str());
            TH1D *histProjInt = (TH1D*) hist2D -> ProjectionX(Form("Proj_%s", cut.c_str()));
            histProjInt -> Write(Form("Proj_%s", cut.c_str()));


            for (int iPt = 0;iPt < nPtBins1;iPt++) {
                TH1D *histProjPt = (TH1D*) ProjectTH2(hist2D, minPtBins1[iPt], maxPtBins1[iPt]);
                histProjPt -> Write(Form("Proj_%s__Pt_%1.0f_%1.0f", cut.c_str(), minPtBins1[iPt], maxPtBins1[iPt]));
            }
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
////////////////////////////////////////////////////////////////////////////////
TH1D* ProjectTHnSparse(THnSparseD *histSparse, double minVar1, double maxVar1, double minVar2, double maxVar2) {
    double minVar1Bin = histSparse -> GetAxis(1) -> FindBin(minVar1);
    double maxVar1Bin = histSparse -> GetAxis(1) -> FindBin(maxVar1 - 0.01);
    double minVar2Bin = histSparse -> GetAxis(2) -> FindBin(minVar2);
    double maxVar2Bin = histSparse -> GetAxis(2) -> FindBin(maxVar2 - 0.01);

    Printf("minVar1Bin = %1.0f, maxVar1Bin = %1.0f, minVar2Bin = %1.0f, maxVar2Bin = %1.0f", minVar1Bin, maxVar1Bin, minVar2Bin, maxVar2Bin);
    histSparse -> GetAxis(1) -> SetRange(minVar1Bin, maxVar1Bin);
    histSparse -> GetAxis(2) -> SetRange(minVar2Bin, maxVar2Bin);

    TH1D *histProj = (TH1D*) histSparse -> Projection(0, Form("histProj_%1.0f_%1.0f__%1.0f_%1.0f", minVar1, maxVar1, minVar2, maxVar2));
    return histProj;
}
////////////////////////////////////////////////////////////////////////////////
TH1D* ProjectTH2(TH2D *hist2D, double minPtRange, double maxPtRange) {
    double minPtBin = hist2D -> GetYaxis() -> FindBin(minPtRange);
    double maxPtBin = hist2D -> GetYaxis() -> FindBin(maxPtRange - 0.01);

    Printf("minPtBin = %1.0f, maxPtBin = %1.0f", minPtBin, maxPtBin);
    hist2D -> GetYaxis() -> SetRange(minPtBin, maxPtBin);

    TH1D *histProj = (TH1D*) hist2D -> ProjectionX(Form("histProj_%1.0f_%1.0f", minPtRange, maxPtRange));
    return histProj;
}