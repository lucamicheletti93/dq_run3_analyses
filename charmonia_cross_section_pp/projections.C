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
    //string productionName = "LHC24e5";
    //string associationType = "time_association";
    // Data settings
    string dataset = "2024"; //
    //string dataset = "2023";
    string productionName = "LHC24aj_pass1_skimmed"; //
    //string productionName = "LHC23_pass4_thinned";
    //string associationType = "time_association";
    string associationType = "time_assoc_fDiMuon"; // fDiMuon,fSingleMuLow,std_assoc_fSingleMuLow,time_assoc_fSingleMuLow,std_assoc_fDiMuon,time_assoc_fDiMuon
    //string pathName = "/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data"; //
    //string pathName = "/Users/lucamicheletti/cernbox/JPSI/Run3";
    string pathName = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_pp_cross_section_run3/data";

    string histName = "Mass_Pt_Rapidity";
    //string histName = "Mass_Pt";

    bool hasMixedEvent = false;

    // pt differential
    //const int nPtBins1 = 8;
    //double minPtBins1[] = {0, 1, 2, 3, 4, 5, 7, 10};
    //double maxPtBins1[] = {1, 2, 3, 4, 5, 7, 10, 20};
    //const int nPtBins1 = 14;
    //double minPtBins1[] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 10, 15};
    //double maxPtBins1[] = {0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 10, 15, 20};

    //const int nPtBins1 = 13;
    //double minPtBins1[] = {0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 15};
    //double maxPtBins1[] = {0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20};
    const int nPtBins1 = 10;
    double minPtBins1[] = {0, 0.5, 1, 2, 3, 4, 5, 6, 7, 10};
    double maxPtBins1[] = {0.5, 1, 2, 3, 4, 5, 6, 7, 10, 20};
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
        "PairsMuonSEPM_matchedMchMid",
        "PairsMuonSEPP_matchedMchMid",
        "PairsMuonSEMM_matchedMchMid",
        "PairsMuonSEPM_matchedMchMid_mumuFromJpsi",
        "PairsMuonSEPM_matchedMchMid_mumuFromPsi2S",
        "PairsMuonSEPM_muonLowPt210SigmaPDCA",
        "PairsMuonSEPP_muonLowPt210SigmaPDCA",
        "PairsMuonSEMM_muonLowPt210SigmaPDCA"
        //"PairsMuonSEPM_muonLowPt10SigmaPDCA_mumuFromPsi2S",
        //"PairsMuonSEPM_muonLowPt210SigmaPDCA_mumuFromPsi2S",
        //"PairsMuonSEPM_muonLowPt510SigmaPDCA_mumuFromPsi2S"
    };
    */

    string fInName = "AnalysisResults";
    string cuts[] = {
        /*"PairsMuonSEPM_matchedMchMid",
        "PairsMuonSEPP_matchedMchMid",
        "PairsMuonSEMM_matchedMchMid",*/
        "PairsMuonSEPM_muonLowPt210SigmaPDCA",
        "PairsMuonSEPP_muonLowPt210SigmaPDCA",
        "PairsMuonSEMM_muonLowPt210SigmaPDCA"
    };

    TFile *fIn = new TFile(Form("%s/%s/%s/%s_%s.root", pathName.c_str(), dataset.c_str(), productionName.c_str(), fInName.c_str(), associationType.c_str()), "READ");
    fIn -> ls();
    TList *list1 = (TList*) fIn -> Get("analysis-same-event-pairing/output");

    TFile *fOut = new TFile(Form("%s/%s/%s/Histograms_%s_%s.root", pathName.c_str(), dataset.c_str(), productionName.c_str(), fInName.c_str(), associationType.c_str()), "RECREATE");
    for (auto& cut : cuts) {
        TList *listSE = (TList*) list1 -> FindObject(cut.c_str());

        if (histName == "Mass_Pt_Rapidity") {
            THnSparseD *histSparse = (THnSparseD*) listSE -> FindObject(histName.c_str());
            histSparse -> GetAxis(1) -> SetRangeUser(0., 20.);
            histSparse -> GetAxis(2) -> SetRangeUser(2.5, 4);

            TH1D *histProjInt = (TH1D*) histSparse -> Projection(0, Form("Proj_%s", cut.c_str()));
            histProjInt -> Write(Form("Proj_%s", cut.c_str()));

            TH1D *histProjPtInt = (TH1D*) histSparse -> Projection(1, Form("Proj_Pt_%s", cut.c_str()));
            histProjPtInt -> Write(Form("Proj_Pt_%s", cut.c_str()));

            TH1D *histProjRapInt = (TH1D*) histSparse -> Projection(2, Form("Proj_Rap_%s", cut.c_str()));
            histProjRapInt -> Write(Form("Proj_Rap_%s", cut.c_str()));


            for (int iPt = 0;iPt < nPtBins1;iPt++) {
                TH1D *histProjPt = (TH1D*) ProjectTHnSparse(histSparse, minPtBins1[iPt], maxPtBins1[iPt], minRapBins1[0], maxRapBins1[0]);
                histProjPt -> Write(Form("Proj_%s__Pt_%2.1f_%2.1f", cut.c_str(), minPtBins1[iPt], maxPtBins1[iPt]));
            }

            for (int iRap = 0;iRap < nRapBins2;iRap++) {
                TH1D *histProjRap = (TH1D*) ProjectTHnSparse(histSparse, minPtBins2[0], maxPtBins2[0], minRapBins2[iRap], maxRapBins2[iRap]);
                histProjRap -> Write(Form("Proj_%s__Rap_%3.2f_%3.2f", cut.c_str(), minRapBins2[iRap], maxRapBins2[iRap]));
            }
        }

        if (histName == "Mass_Pt") {
            TH2D *hist2D = (TH2D*) listSE -> FindObject(histName.c_str());
            TH1D *histProjInt = (TH1D*) hist2D -> ProjectionX(Form("Proj_%s", cut.c_str()));
            histProjInt -> Write(Form("Proj_%s", cut.c_str()));


            for (int iPt = 0;iPt < nPtBins1;iPt++) {
                TH1D *histProjPt = (TH1D*) ProjectTH2(hist2D, minPtBins1[iPt], maxPtBins1[iPt]);
                histProjPt -> Write(Form("Proj_%s__Pt_%1.0f_%1.0f", cut.c_str(), minPtBins1[iPt], maxPtBins1[iPt]));
            }
        }

        if (hasMixedEvent) {
            TList *list2 = (TList*) fIn -> Get("analysis-event-mixing/output");
            auto cutME = cut.replace(9, 2, "ME");
            TList *listME = (TList*) list2 -> FindObject(cutME.c_str());

            if (histName == "Mass_Pt_Rapidity") {
                THnSparseD *histSparseME = (THnSparseD*) listME -> FindObject(histName.c_str());
                histSparseME -> GetAxis(1) -> SetRangeUser(0., 20.);
                histSparseME -> GetAxis(2) -> SetRangeUser(2.5, 4);

                TH1D *histProjIntME = (TH1D*) histSparseME -> Projection(0, Form("Proj_%s", cut.c_str()));
                histProjIntME -> Write(Form("Proj_%s", cut.c_str()));

                TH1D *histProjPtIntME = (TH1D*) histSparseME -> Projection(1, Form("Proj_Pt_%s", cut.c_str()));
                histProjPtIntME -> Write(Form("Proj_Pt_%s", cut.c_str()));

                TH1D *histProjRapIntME = (TH1D*) histSparseME -> Projection(2, Form("Proj_Rap_%s", cut.c_str()));
                histProjRapIntME -> Write(Form("Proj_Rap_%s", cut.c_str()));

                for (int iPt = 0;iPt < nPtBins1;iPt++) {
                    TH1D *histProjPtME = (TH1D*) ProjectTHnSparse(histSparseME, minPtBins1[iPt], maxPtBins1[iPt], minRapBins1[0], maxRapBins1[0]);
                    histProjPtME -> Write(Form("Proj_%s__Pt_%2.1f_%2.1f", cut.c_str(), minPtBins1[iPt], maxPtBins1[iPt]));
                }

                for (int iRap = 0;iRap < nRapBins2;iRap++) {
                    TH1D *histProjRapME = (TH1D*) ProjectTHnSparse(histSparseME, minPtBins2[0], maxPtBins2[0], minRapBins2[iRap], maxRapBins2[iRap]);
                    histProjRapME -> Write(Form("Proj_%s__Rap_%3.2f_%3.2f", cut.c_str(), minRapBins2[iRap], maxRapBins2[iRap]));
                }
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