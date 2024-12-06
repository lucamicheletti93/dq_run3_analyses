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
void SetHistogram(TH1D *, Color_t , double);
TH1D* ComputeRfactor(TH1D *, TH1D *, TH1D *, string);
double ComputeFfactor(TH1D *, TH1D *, TH1D *, TH1D *);
double ComputeScaleFactor(TH1D *, TH1D *, TH1D *, double [2], double[2]);
double ComputeSideBandsNorm(TH1D *, double [2], double[2]);
TH1D* ComputeLsBackground(TH1D *, TH1D *, string);
TH1D* CreateSideBandsHist(TH1D *, double [2], double[2], string);

void event_mixing() {
    //LoadStyle();

    const int nPtBins = 13;
    double minPtBins[] = {0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 15};
    double maxPtBins[] = {0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20};

    double minRanges[13][2] = {{2.1, 4.9},{2.1, 4.9},{2.1, 4.9},{2.1, 4.9},{2.1, 4.9},{1.9, 4.5},{1.9, 4.5},{1.9, 4.5},{1.9, 4.5},{1.9, 4.5},{1.9, 4.5},{1.9, 4.5},{1.9, 4.5}};
    double maxRanges[13][2] = {{2.5, 5.0},{2.5, 5.0},{2.5, 5.0},{2.4, 5.0},{2.4, 5.0},{2.2, 5.0},{2.2, 5.0},{2.2, 5.0},{2.2, 5.0},{2.2, 5.0},{2.2, 5.0},{2.2, 5.0},{2.2, 5.0}};
    double minRange[2];
    double maxRange[2];
    double fFactor = 0;

    // Same Event files [SE]
    //string fInNameSE = "/Users/lucamicheletti/cernbox/JPSI/Run3/MC/LHC24e5/Histograms_AnalysisResults_dq_efficiency_time_association.root";
    //string fInNameSE = "/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data/2022/LHC22_pass6_minimum_bias/Histograms_AnalysisResults_time_association.root";
    //string fInNameSE = "/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data/2023/LHC23_pass4_thinned/Histograms_AnalysisResults_time_association.root";
    //string fInNameSE = "/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data/2024/LHC24_pass1_skimmed/Histograms_AnalysisResults_fDiMuon.root";
    //string fInNameSE = "/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data/2024/LHC24_pass1_skimmed/Histograms_AnalysisResults_fDiMuon.root"; // fDiMuon,fSingleMuLow
    string fInNameSE = "/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data/2024/LHC24al_pass1_skimmed/Histograms_AnalysisResults_std_assoc_fDiMuon.root"; // fDiMuon,fSingleMuLow
    string cutNameSE = "muonLowPt210SigmaPDCA"; // matchedMchMid,muonLowPt210SigmaPDCA

    // Mixed Event files [ME]
    //string fInNameME = "/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data/2024/LHC24_pass1_skimmed/Histograms_AnalysisResults_fSingleMuLow.root";
    string fInNameME = "/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data/2024/LHC24al_pass1_skimmed/Histograms_AnalysisResults_std_assoc_fSingleMuLow.root";
    string cutNameME = "muonLowPt210SigmaPDCA"; // matchedMchMid,muonLowPt210SigmaPDCA

    TFile *fInSE = new TFile(fInNameSE.c_str(), "READ");
    TFile *fInME = new TFile(fInNameME.c_str(), "READ");

    TFile *fOut = new TFile("Histograms_AnalysisResults_fDiMuon_noBkg.root", "RECREATE");

    //ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo//
    // Integrated spectra
    //ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo//
    TH1D *histMassIntDiMuonSEPM = (TH1D*) fInSE -> Get(Form("Proj_PairsMuonSEPM_%s", cutNameSE.c_str())); SetHistogram(histMassIntDiMuonSEPM, kBlack, 2);
    TH1D *histMassIntDiMuonSEPP = (TH1D*) fInSE -> Get(Form("Proj_PairsMuonSEPP_%s", cutNameSE.c_str())); SetHistogram(histMassIntDiMuonSEPP, kBlack, 2);
    TH1D *histMassIntDiMuonSEMM = (TH1D*) fInSE -> Get(Form("Proj_PairsMuonSEMM_%s", cutNameSE.c_str())); SetHistogram(histMassIntDiMuonSEMM, kBlack, 2);

    TH1D *histMassIntDiMuonMEPM = (TH1D*) fInSE -> Get(Form("Proj_PairsMuonMEPM_%s", cutNameSE.c_str())); SetHistogram(histMassIntDiMuonMEPM, kBlack, 2);
    TH1D *histMassIntDiMuonMEPP = (TH1D*) fInSE -> Get(Form("Proj_PairsMuonMEPP_%s", cutNameSE.c_str())); SetHistogram(histMassIntDiMuonMEPP, kBlack, 2);
    TH1D *histMassIntDiMuonMEMM = (TH1D*) fInSE -> Get(Form("Proj_PairsMuonMEMM_%s", cutNameSE.c_str())); SetHistogram(histMassIntDiMuonMEMM, kBlack, 2);

    TH1D *histMassIntSingleMuLowSEPM = (TH1D*) fInME -> Get(Form("Proj_PairsMuonSEPM_%s", cutNameME.c_str())); SetHistogram(histMassIntSingleMuLowSEPM, kBlack, 2);
    TH1D *histMassIntSingleMuLowSEPP = (TH1D*) fInME -> Get(Form("Proj_PairsMuonSEPP_%s", cutNameME.c_str())); SetHistogram(histMassIntSingleMuLowSEPP, kAzure+2, 2);
    TH1D *histMassIntSingleMuLowSEMM = (TH1D*) fInME -> Get(Form("Proj_PairsMuonSEMM_%s", cutNameME.c_str())); SetHistogram(histMassIntSingleMuLowSEMM, kAzure+2, 2);

    TH1D *histMassIntSingleMuLowMEPM = (TH1D*) fInME -> Get(Form("Proj_PairsMuonMEPM_%s", cutNameME.c_str())); SetHistogram(histMassIntSingleMuLowMEPM, kAzure+2, 2);
    TH1D *histMassIntSingleMuLowMEPP = (TH1D*) fInME -> Get(Form("Proj_PairsMuonMEPP_%s", cutNameME.c_str())); SetHistogram(histMassIntSingleMuLowMEPP, kAzure+2, 2);
    TH1D *histMassIntSingleMuLowMEMM = (TH1D*) fInME -> Get(Form("Proj_PairsMuonMEMM_%s", cutNameME.c_str())); SetHistogram(histMassIntSingleMuLowMEMM, kAzure+2, 2);

    TH1D *histMassIntDiMuonSEPPMM = ComputeLsBackground(histMassIntDiMuonSEPP, histMassIntDiMuonSEMM, Form("histMassIntPPMM_%s", cutNameME.c_str())); SetHistogram(histMassIntDiMuonSEPPMM, kBlue+1, 2);
    TH1D *histMassIntDiMuonPPMM = ComputeLsBackground(histMassIntDiMuonSEPP, histMassIntDiMuonSEMM, "histMassIntDiMuonPPMM"); SetHistogram(histMassIntDiMuonPPMM, kBlue+1, 2);
    TH1D *histMassIntSingleMuLowPPMM = ComputeLsBackground(histMassIntSingleMuLowSEPP, histMassIntSingleMuLowSEMM, "histMassIntSingleMuLowPPMM"); SetHistogram(histMassIntSingleMuLowPPMM, kBlue+1, 2);
    TH1D *histMassIntDiMuonMEPPMM = ComputeLsBackground(histMassIntDiMuonMEPP, histMassIntDiMuonMEMM, "histMassIntDiMuonMEPPMM"); SetHistogram(histMassIntDiMuonMEPPMM, kBlue+1, 2);

    //double normIntDiMuonSEPPMM = ComputeScaleFactor(histMassIntDiMuonSEPM, histMassIntDiMuonSEPP, histMassIntDiMuonSEMM, minRanges, maxRanges);
    //double normIntDiMuonMEPPMM = ComputeScaleFactor(histMassIntDiMuonSEPM, histMassIntDiMuonMEPP, histMassIntDiMuonMEMM, minRanges, maxRanges);
    //double normIntSingleMuLowSEPPMM = ComputeScaleFactor(histMassIntSingleMuLowSEPM, histMassIntSingleMuLowSEPP, histMassIntSingleMuLowSEMM, minRanges, maxRanges);
    //double scaleFactorSingleMuLowMEPM = normIntSingleMuLowSEPPMM / normIntDiMuonSEPPMM;
    //double scaleFactorDiMuonMEPM = normIntDiMuonMEPPMM / normIntDiMuonSEPPMM;

    // Normalize using sidebands -> Not working, the shape is different
    //double sbNormDiMuonSEPM = ComputeSideBandsNorm(histMassIntDiMuonSEPM, minRanges, maxRanges);
    //double sbNormDiMuonMEPM = ComputeSideBandsNorm(histMassIntDiMuonMEPM, minRanges, maxRanges);
    //double scaleFactorDiMuonMEPM = sbNormDiMuonSEPM / sbNormDiMuonMEPM; 

    TH1D *histMassIntSingleMuLowRfactor = ComputeRfactor(histMassIntSingleMuLowMEPM, histMassIntSingleMuLowMEPP, histMassIntSingleMuLowMEMM, "histMassIntSingleMuLowRfactor"); SetHistogram(histMassIntSingleMuLowRfactor, kBlue+1, 2);
    TH1D *histMassIntDiMuonRfactor = ComputeRfactor(histMassIntDiMuonMEPM, histMassIntDiMuonMEPP, histMassIntDiMuonMEMM, "histMassIntDiMuonRfactor"); SetHistogram(histMassIntDiMuonRfactor, kRed+1, 2);

    TCanvas *canvasRfactor = new TCanvas("canvasRfactor", "", 800, 600);
    histMassIntSingleMuLowRfactor -> Draw("EP");
    histMassIntDiMuonRfactor -> Draw("EP SAME");
    TLine *lineUnity = new TLine(0, 1, 5, 1);
    lineUnity -> SetLineColor(kGray+1);
    lineUnity -> SetLineStyle(kDashed);
    lineUnity -> Draw("SAME");

    // Compute F-factor for normalization
    fFactor = ComputeFfactor(histMassIntSingleMuLowMEPM, histMassIntDiMuonSEPP, histMassIntDiMuonSEMM, histMassIntSingleMuLowRfactor);
    std::cout << "F-factor = " << fFactor << std::endl;
    histMassIntSingleMuLowMEPM -> Scale(fFactor); // Warning!!!!

    // Compute scale factors
    std::copy(std::begin(minRanges[0]), std::end(minRanges[0]), minRange);
    std::copy(std::begin(maxRanges[0]), std::end(maxRanges[0]), maxRange);

    double normIntDiMuonSEPM = ComputeSideBandsNorm(histMassIntDiMuonSEPM, minRange, maxRange);
    double normIntDiMuonSEPPMM = ComputeSideBandsNorm(histMassIntDiMuonSEPPMM, minRange, maxRange);
    double normIntDiMuonMEPPMM = ComputeSideBandsNorm(histMassIntDiMuonMEPPMM, minRange, maxRange);
    double normIntSingleMuLowMEPM = ComputeSideBandsNorm(histMassIntSingleMuLowMEPM, minRange, maxRange);
    double scaleFactorSingleMuLowMEPM = normIntDiMuonSEPM / normIntSingleMuLowMEPM;
    double scaleFactorDiMuonSEPPMM = normIntDiMuonSEPM / normIntDiMuonSEPPMM;
    double scaleFactorDiMuonMEPPMM = normIntDiMuonSEPM / normIntDiMuonMEPPMM;

    histMassIntSingleMuLowMEPM -> Scale(scaleFactorSingleMuLowMEPM);  // Warning!!!!
    histMassIntDiMuonSEPPMM -> Scale(scaleFactorDiMuonSEPPMM); // Warning!!!!
    histMassIntDiMuonMEPPMM -> Scale(scaleFactorDiMuonMEPPMM); // Warning!!!!

    // Remove combinatorial background
    TH1D *histMassIntDiMuonBkgSubtrSEPM = (TH1D*) histMassIntDiMuonSEPM -> Clone("histMassIntDiMuonBkgSubtrSEPM");
    histMassIntDiMuonBkgSubtrSEPM -> Add(histMassIntSingleMuLowMEPM, -1);
    SetHistogram(histMassIntDiMuonBkgSubtrSEPM, kRed+1, 2);

    TCanvas *canvasIntMassDiMuon = new TCanvas("canvasIntMassDiMuon", "fDiMuon events", 800, 600);
    gStyle -> SetOptStat(false);
    gPad -> SetLogy(true);
    histMassIntDiMuonSEPM -> SetTitle("");
    histMassIntDiMuonSEPM -> GetXaxis() -> SetRangeUser(1.9, 5);
    histMassIntDiMuonSEPM -> GetYaxis() -> SetRangeUser(1, 1e6);
    histMassIntDiMuonSEPM -> Draw("EP");
    histMassIntSingleMuLowMEPM -> Draw("EP SAME");
    histMassIntDiMuonBkgSubtrSEPM -> Draw("EP SAME");
    histMassIntDiMuonPPMM -> Draw("EP SAME");
    histMassIntDiMuonMEPPMM -> Draw("EP SAME");

    TLegend *legendIntMassDiMuon = new TLegend(0.50, 0.70, 0.70, 0.89);
    SetLegend(legendIntMassDiMuon);
    legendIntMassDiMuon -> AddEntry(histMassIntDiMuonSEPM, "SE PM - fDiMuon", "PL");
    legendIntMassDiMuon -> AddEntry(histMassIntSingleMuLowMEPM, "ME PM - fSingleMuLow", "PL");
    legendIntMassDiMuon -> AddEntry(histMassIntDiMuonPPMM, "SE PPMM - fDiMuon", "PL");
    legendIntMassDiMuon -> AddEntry(histMassIntDiMuonBkgSubtrSEPM, "(SE PM) - (ME PM)", "PL");
    legendIntMassDiMuon -> Draw("SAME");

    TCanvas *canvasIntMassSingleMuLow = new TCanvas("canvasIntMassSingleMuLow", "fSingleMuLow events", 800, 600);
    gStyle -> SetOptStat(false);
    gPad -> SetLogy(true);
    histMassIntSingleMuLowSEPM -> SetTitle("");
    histMassIntSingleMuLowSEPM -> GetXaxis() -> SetRangeUser(1.9, 5);
    histMassIntSingleMuLowSEPM -> GetYaxis() -> SetRangeUser(1, 1e6);
    histMassIntSingleMuLowSEPM -> Draw("EP");
    histMassIntSingleMuLowPPMM -> Draw("EP SAME");

    TLegend *legendIntMassSingleMuLow = new TLegend(0.50, 0.70, 0.70, 0.89);
    SetLegend(legendIntMassSingleMuLow);
    legendIntMassSingleMuLow -> AddEntry(histMassIntSingleMuLowSEPM, "SE PM - fSingleMuLow", "PL");
    legendIntMassSingleMuLow -> AddEntry(histMassIntSingleMuLowPPMM, "SE PPMM - fSingleMuLow", "PL");
    legendIntMassSingleMuLow -> Draw("SAME");

    fOut -> cd();
    histMassIntDiMuonSEPM -> Write();
    histMassIntDiMuonBkgSubtrSEPM -> Write();

    //ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo//
    // Pt differential spectrum
    //ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo ooOOoo//
    TCanvas *canvasMassPt = new TCanvas("canvasMassPt", "SEPM signal only", 1800, 1800); canvasMassPt -> Divide(5, 3);
    TCanvas *canvasMassPtBkgs = new TCanvas("canvasMassPtBkgs", "comparisong OS vs LS / ME backgrounds", 1800, 1800); canvasMassPtBkgs -> Divide(5, 3);
    TCanvas *canvasMassPtSideBands = new TCanvas("canvasMassPtSideBands", "comparisong OS vs LS / ME backgrounds", 1800, 1800); canvasMassPtSideBands -> Divide(5, 3);
    //TCanvas *canvasMassPtRatio = new TCanvas("canvasMassPtRatio", "Ratio OS vs LS / ME backgrounds", 1800, 1800); canvasMassPtRatio -> Divide(5, 3);
    //TCanvas *canvasMassPtRfactor = new TCanvas("canvasMassPtRfactor", "R-factor calculation", 1800, 1800); canvasMassPtRfactor -> Divide(5, 3);

    for (int iPt = 0;iPt < nPtBins;iPt++) {
        TH1D *histMassPtDiMuonSEPM = (TH1D*) fInSE -> Get(Form("Proj_PairsMuonSEPM_%s__Pt_%2.1f_%2.1f", cutNameSE.c_str(), minPtBins[iPt], maxPtBins[iPt])); histMassPtDiMuonSEPM -> GetXaxis() -> SetRangeUser(1.9, 5);
        SetHistogram(histMassPtDiMuonSEPM, kBlack, 2);
        TH1D *histMassPtDiMuonSEPP = (TH1D*) fInSE -> Get(Form("Proj_PairsMuonSEPP_%s__Pt_%2.1f_%2.1f", cutNameSE.c_str(), minPtBins[iPt], maxPtBins[iPt])); histMassPtDiMuonSEPP -> GetXaxis() -> SetRangeUser(1.9, 5);
        TH1D *histMassPtDiMuonSEMM = (TH1D*) fInSE -> Get(Form("Proj_PairsMuonSEMM_%s__Pt_%2.1f_%2.1f", cutNameSE.c_str(), minPtBins[iPt], maxPtBins[iPt])); histMassPtDiMuonSEMM -> GetXaxis() -> SetRangeUser(1.9, 5);

        TH1D *histMassPtSingleMuLowSEPM = (TH1D*) fInME -> Get(Form("Proj_PairsMuonSEPM_%s__Pt_%2.1f_%2.1f", cutNameME.c_str(), minPtBins[iPt], maxPtBins[iPt])); histMassPtSingleMuLowSEPM -> GetXaxis() -> SetRangeUser(1.9, 5);
        TH1D *histMassPtSingleMuLowSEPP = (TH1D*) fInME -> Get(Form("Proj_PairsMuonSEPP_%s__Pt_%2.1f_%2.1f", cutNameME.c_str(), minPtBins[iPt], maxPtBins[iPt])); histMassPtSingleMuLowSEPP -> GetXaxis() -> SetRangeUser(1.9, 5);
        TH1D *histMassPtSingleMuLowSEMM = (TH1D*) fInME -> Get(Form("Proj_PairsMuonSEMM_%s__Pt_%2.1f_%2.1f", cutNameME.c_str(), minPtBins[iPt], maxPtBins[iPt])); histMassPtSingleMuLowSEMM -> GetXaxis() -> SetRangeUser(1.9, 5);

        TH1D *histMassPtSingleMuLowMEPM = (TH1D*) fInME -> Get(Form("Proj_PairsMuonMEPM_%s__Pt_%2.1f_%2.1f", cutNameME.c_str(), minPtBins[iPt], maxPtBins[iPt])); histMassPtSingleMuLowMEPM -> GetXaxis() -> SetRangeUser(1.9, 5);
        SetHistogram(histMassPtSingleMuLowMEPM, kAzure+1, 2);
        TH1D *histMassPtSingleMuLowMEPP = (TH1D*) fInME -> Get(Form("Proj_PairsMuonMEPP_%s__Pt_%2.1f_%2.1f", cutNameME.c_str(), minPtBins[iPt], maxPtBins[iPt])); histMassPtSingleMuLowMEPP -> GetXaxis() -> SetRangeUser(1.9, 5);
        TH1D *histMassPtSingleMuLowMEMM = (TH1D*) fInME -> Get(Form("Proj_PairsMuonMEMM_%s__Pt_%2.1f_%2.1f", cutNameME.c_str(), minPtBins[iPt], maxPtBins[iPt])); histMassPtSingleMuLowMEMM -> GetXaxis() -> SetRangeUser(1.9, 5);

        TH1D *histMassPtDiMuonSEPPMM = ComputeLsBackground(histMassPtDiMuonSEPP, histMassPtDiMuonSEMM, Form("histMassPtPPMM_%s__Pt_%2.1f_%2.1f", cutNameME.c_str(), minPtBins[iPt], maxPtBins[iPt]));
        SetHistogram(histMassPtDiMuonSEPPMM, kBlue+1, 2);
        TH1D *histMassPtScaledDiMuonSEPPMM = (TH1D*) histMassPtDiMuonSEPPMM -> Clone(Form("histMassPtScaledPPMM_%s__Pt_%2.1f_%2.1f", cutNameME.c_str(), minPtBins[iPt], maxPtBins[iPt]));

        TH1D *histMassPtRfactor = ComputeRfactor(histMassPtSingleMuLowMEPM, histMassPtSingleMuLowMEPP, histMassPtSingleMuLowMEMM, Form("histMassPtRfactor_%s__Pt_%2.1f_%2.1f", cutNameME.c_str(), minPtBins[iPt], maxPtBins[iPt]));
        
        //canvasMassPtRfactor -> cd(iPt+1);
        //histMassPtRfactor -> Draw("EP");
        //lineUnity -> Draw("SAME");

        // Compute F-factor for normalization
        fFactor = ComputeFfactor(histMassPtSingleMuLowMEPM, histMassPtDiMuonSEPP, histMassPtDiMuonSEMM, histMassPtRfactor);
        histMassPtSingleMuLowMEPM -> Scale(fFactor); // Warning!!!!

        // Compute scale factors
        //double normPtDiMuonSEPPMM = ComputeScaleFactor(histMassPtDiMuonSEPM, histMassPtDiMuonSEPP, histMassPtDiMuonSEMM, minRanges, maxRanges);
        //double normPtSingleMuLowSEPPMM = ComputeScaleFactor(histMassPtSingleMuLowSEPM, histMassPtSingleMuLowSEPP, histMassPtSingleMuLowSEMM, minRanges, maxRanges);
        //double normPtSingleMuLowMEPPMM = ComputeScaleFactor(histMassPtSingleMuLowMEPM, histMassPtSingleMuLowMEPP, histMassPtSingleMuLowMEMM, minRanges, maxRanges); 
        //double scaleFactorSingleMuLowMEPM = normPtSingleMuLowSEPPMM / normPtDiMuonSEPPMM;

        std::copy(std::begin(minRanges[iPt]), std::end(minRanges[iPt]), minRange);
        std::copy(std::begin(maxRanges[iPt]), std::end(maxRanges[iPt]), maxRange);
        std::cout << minRange[0] << " " << minRange[1] << std::endl;
        double normPtDiMuonSEPM = ComputeSideBandsNorm(histMassPtDiMuonSEPM, minRange, maxRange);
        double normPtDiMuonSEPPMM = ComputeSideBandsNorm(histMassPtScaledDiMuonSEPPMM, minRange, maxRange);
        double normPtSingleMuLowMEPM = ComputeSideBandsNorm(histMassPtSingleMuLowMEPM, minRange, maxRange);
        double scaleFactorSingleMuLowMEPM = normPtDiMuonSEPM / normPtSingleMuLowMEPM;
        double scaleFactorDiMuonSEPPMM = normPtDiMuonSEPM / normPtDiMuonSEPPMM;

        histMassPtSingleMuLowMEPM -> Scale(0.5 * scaleFactorSingleMuLowMEPM); // Warning!!!!
        histMassPtScaledDiMuonSEPPMM -> Scale(0.5 * scaleFactorDiMuonSEPPMM); // Warning!!!!

        std::cout << "F-factor = " << fFactor  << "; Scale-factor SingleMuLowMEPM = " << scaleFactorSingleMuLowMEPM << std::endl;

        // Remove combinatorial background
        TH1D *histMassPtDiMuonBkgSubtrSEPM = (TH1D*) histMassPtDiMuonSEPM -> Clone(Form("histMassPtDiMuonBkgSubtrSEPM_%s__Pt_%2.1f_%2.1f", cutNameME.c_str(), minPtBins[iPt], maxPtBins[iPt]));
        histMassPtDiMuonBkgSubtrSEPM -> Add(histMassPtSingleMuLowMEPM, -1);
        SetHistogram(histMassPtDiMuonBkgSubtrSEPM, kRed+1, 2);

        TH1D *histMassPtRatio_DiMuonSEPM_SingleMuLowMEPM = (TH1D*) histMassPtDiMuonSEPM -> Clone(Form("histMassPtRatio_DiMuonSEPM_SingleMuLowMEPM__Pt_%2.1f_%2.1f", minPtBins[iPt], maxPtBins[iPt]));
        histMassPtRatio_DiMuonSEPM_SingleMuLowMEPM -> Divide(histMassPtSingleMuLowMEPM);

        TH1D *histMassPtRatio_DiMuonSEPM_DiMuonPPMM = (TH1D*) histMassPtDiMuonSEPM -> Clone(Form("histMassPtRatio_DiMuonSEPM_DiMuonPPMM__Pt_%2.1f_%2.1f", minPtBins[iPt], maxPtBins[iPt]));
        histMassPtRatio_DiMuonSEPM_DiMuonPPMM -> Divide(histMassPtScaledDiMuonSEPPMM);

        //------------------------------------------------------//
        // Plot mass histo vs pT for DiMuon only
        canvasMassPt -> cd(iPt+1);
        gStyle -> SetOptStat(false);
        gPad -> SetLogy(false);
        histMassPtDiMuonSEPM -> SetTitle(Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}", minPtBins[iPt], maxPtBins[iPt]));
        histMassPtDiMuonSEPM -> Draw("EP");
        histMassPtDiMuonSEPPMM -> Draw("EP SAME");

        TLegend *legendMassPt = new TLegend(0.60, 0.70, 0.80, 0.89);
        SetLegend(legendMassPt);
        legendMassPt -> AddEntry(histMassPtDiMuonSEPM, "OS [fDiMuon]", "PL");
        legendMassPt -> AddEntry(histMassPtDiMuonSEPPMM, "LS [fDiMuon]", "PL");
        legendMassPt -> Draw("SAME");

        //------------------------------------------------------//
        // Plot mass histo vs pT for DiMuon only with backgrounds
        canvasMassPtBkgs -> cd(iPt+1);
        gStyle -> SetOptStat(false);
        gPad -> SetLogy(true);
        histMassPtDiMuonSEPM -> SetTitle(Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}", minPtBins[iPt], maxPtBins[iPt]));
        histMassPtDiMuonSEPM -> Draw("EP");
        histMassPtSingleMuLowMEPM -> Draw("EP SAME");
        //histMassPtDiMuonBkgSubtrSEPM -> Draw("EP SAME");
        //histMassPtScaledDiMuonSEPPMM -> Draw("EP SAME");

        TLegend *legendMassPtBkgs = new TLegend(0.50, 0.70, 0.70, 0.89);
        SetLegend(legendMassPtBkgs);
        legendMassPtBkgs -> AddEntry(histMassPtDiMuonSEPM, "SE PM - fDiMuon", "PL");
        legendMassPtBkgs -> AddEntry(histMassPtSingleMuLowMEPM, "ME PM - fSingleMuLow", "PL");
        //legendMassPtBkgs -> AddEntry(histMassPtScaledDiMuonSEPPMM, "SE PPMM - fDiMuon", "PL");
        legendMassPtBkgs -> Draw("SAME");

        // Plot hist sidebands
        TH1D *histMassPtDiMuonSEPMSideBands = CreateSideBandsHist(histMassPtDiMuonSEPM, minRange, maxRange, Form("histMassPtDiMuonSEPMSideBands_%s__Pt_%2.1f_%2.1f", cutNameME.c_str(), minPtBins[iPt], maxPtBins[iPt])); 
        SetHistogram(histMassPtDiMuonSEPMSideBands, kRed, 2); 

        canvasMassPtSideBands -> cd(iPt+1);
        gStyle -> SetOptStat(false);
        gPad -> SetLogy(true);
        histMassPtDiMuonSEPM -> SetTitle(Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}", minPtBins[iPt], maxPtBins[iPt]));
        histMassPtDiMuonSEPMSideBands -> GetXaxis() -> SetRangeUser(2, 5);
        //histMassPtDiMuonSEPM -> GetYaxis() -> SetRangeUser(0.1, 2. * histMassPtDiMuonSEPM -> GetMaximum());
        //histMassPtDiMuonSEPM -> Draw("EP");
        histMassPtDiMuonSEPMSideBands -> Draw("EP");
        histMassPtSingleMuLowMEPM -> Draw("EP SAME");
        //histMassPtScaledDiMuonSEPPMM -> Draw("EP SAME");

        /*TLegend *legendMassPt = new TLegend(0.50, 0.70, 0.70, 0.89);
        SetLegend(legendMassPt);
        legendMassPt -> AddEntry(histMassPtDiMuonSEPM, "SE PM - fDiMuon", "PL");
        legendMassPt -> AddEntry(histMassPtSingleMuLowMEPM, "ME PM - fSingleMuLow", "PL");
        legendMassPt -> AddEntry(histMassPtDiMuonSEPPMM, "SE PPMM - fDiMuon", "PL");
        legendMassPt -> Draw("SAME");*/

        // Plot the ratio
        /*histMassPtRatio_DiMuonSEPM_SingleMuLowMEPM -> SetLineColor(kAzure+1);
        histMassPtRatio_DiMuonSEPM_SingleMuLowMEPM -> SetMarkerColor(kAzure+1);
        histMassPtRatio_DiMuonSEPM_DiMuonPPMM -> SetLineColor(kBlue+1);
        histMassPtRatio_DiMuonSEPM_DiMuonPPMM -> SetMarkerColor(kBlue+1);

        canvasMassPtRatio -> cd(iPt+1);
        gStyle -> SetOptStat(false);
        gPad -> SetLogy(false);

        histMassPtRatio_DiMuonSEPM_SingleMuLowMEPM -> SetTitle(Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}", minPtBins[iPt], maxPtBins[iPt]));
        histMassPtRatio_DiMuonSEPM_SingleMuLowMEPM -> GetYaxis() -> SetRangeUser(0, 1000);
        histMassPtRatio_DiMuonSEPM_SingleMuLowMEPM -> Draw("EP");
        histMassPtRatio_DiMuonSEPM_DiMuonPPMM -> Draw("EP SAME");

        TLegend *legendMassPtRatio = new TLegend(0.20, 0.70, 0.40, 0.89);
        SetLegend(legendMassPtRatio);
        legendMassPtRatio -> AddEntry(histMassPtRatio_DiMuonSEPM_SingleMuLowMEPM, "SE PM fDiMuon / ME PM fSingleMuLow", "PL");
        legendMassPtRatio -> AddEntry(histMassPtRatio_DiMuonSEPM_DiMuonPPMM, "SE PM fDiMuon / SE PPMM fDiMuon", "PL");
        legendMassPtRatio -> Draw("SAME");*/

        fOut -> cd();
        histMassPtDiMuonBkgSubtrSEPM -> Write();
    }
    canvasMassPt -> SaveAs("figures/histMassPt.pdf");
    canvasMassPtBkgs -> SaveAs("figures/histMassPtBkgs.pdf");
    //fOut -> Close();


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
////////////////////////////////////////////////////////////////////////////////
void SetHistogram(TH1D *hist, Color_t color, double size) {
    hist -> SetLineColor(color);
    hist -> SetMarkerColor(color);
    hist -> SetMarkerSize(size);
}
////////////////////////////////////////////////////////////////////////////////
TH1D* ComputeRfactor(TH1D *histPM, TH1D *histPP, TH1D *histMM, string histName) {
    int nBins = histPM -> GetNbinsX();
    double lowEdge = histPM -> GetBinLowEdge(1);
    double topEdge =  histPM -> GetBinLowEdge(nBins) + histPM -> GetBinWidth(nBins);

    TH1D *histRfactor = new TH1D(histName.c_str(), "", nBins, lowEdge, topEdge);
    for (int iBin = 0;iBin < nBins;iBin++) {
        double nPM = histPM -> GetBinContent(iBin+1);
        double nPP = histPP -> GetBinContent(iBin+1);
        double nMM = histMM -> GetBinContent(iBin+1);
        double rFactor = nPM / (2. * TMath::Sqrt(nPP * nMM));
        //std::cout << rFactor << " -> " << nPM << " " << nPP << " " << nMM << std::endl;
        histRfactor -> SetBinContent(iBin+1, std::isnan(rFactor) || std::isinf(rFactor) ? 0 : rFactor);
        histRfactor -> SetBinError(iBin+1, 0);
    }
    return histRfactor;
}
////////////////////////////////////////////////////////////////////////////////
double ComputeFfactor(TH1D *histPM, TH1D *histPP, TH1D *histMM, TH1D *histRfactor) {
    int nBins = histPM -> GetNbinsX();
    double lowEdge = histPM -> GetBinLowEdge(1);
    double topEdge =  histPM -> GetBinLowEdge(nBins) + histPM -> GetBinWidth(nBins);
    double numFfactor = 0;
    double denFfactor = 0;

    for (int iBin = 0;iBin < nBins;iBin++) {
        double nPM = histPM -> GetBinContent(iBin+1);
        double nPP = histPP -> GetBinContent(iBin+1);
        double nMM = histMM -> GetBinContent(iBin+1);
        double rFactor = histRfactor -> GetBinContent(iBin+1);
        numFfactor += 2. * rFactor * TMath::Sqrt(nPP * nMM);
        denFfactor += nPM;
    }
    return numFfactor / denFfactor;
}
////////////////////////////////////////////////////////////////////////////////
TH1D* ComputeLsBackground(TH1D *histPP, TH1D *histMM, string histName) {
    int nBins = histPP -> GetNbinsX();
    double lowEdge = histPP -> GetBinLowEdge(1);
    double topEdge =  histPP -> GetBinLowEdge(nBins) + histPP -> GetBinWidth(nBins);

    TH1D *histLsBkg = new TH1D(histName.c_str(), "", nBins, lowEdge, topEdge);
    for (int iBin = 0;iBin < nBins;iBin++) {
        double nPP = histPP -> GetBinContent(iBin+1);
        double errPP = histPP -> GetBinError(iBin+1);
        double nMM = histMM -> GetBinContent(iBin+1);
        double errMM = histMM -> GetBinError(iBin+1);

        double lsBkg = 2 * TMath::Sqrt(nPP * nMM);

        histLsBkg -> SetBinContent(iBin+1, lsBkg);
        histLsBkg -> SetBinError(iBin+1, TMath::Sqrt(lsBkg));
    }
    return histLsBkg;
}
////////////////////////////////////////////////////////////////////////////////
double ComputeScaleFactor(TH1D *histMassSEPM, TH1D *histMassSEPP, TH1D *histMassSEMM, double minRanges[], double maxRanges[]) {
    double sbLeftSEPM = histMassSEPM -> Integral(histMassSEPM -> FindBin(minRanges[0]), histMassSEPM -> FindBin(maxRanges[0]));
    double sbRightSEPM = histMassSEPM -> Integral(histMassSEPM -> FindBin(minRanges[1]), histMassSEPM -> FindBin(maxRanges[1]));
    double sbLeftSEPP = histMassSEPP -> Integral(histMassSEPP -> FindBin(minRanges[0]), histMassSEPP -> FindBin(maxRanges[0]));
    double sbRightSEPP = histMassSEPP -> Integral(histMassSEPP -> FindBin(minRanges[1]), histMassSEPP -> FindBin(maxRanges[1]));
    double sbLeftSEMM = histMassSEMM -> Integral(histMassSEMM -> FindBin(minRanges[0]), histMassSEMM -> FindBin(maxRanges[0]));
    double sbRightSEMM = histMassSEMM -> Integral(histMassSEMM -> FindBin(minRanges[1]), histMassSEMM -> FindBin(maxRanges[1]));

    double sbLeftSEPPMM = 2 * TMath::Sqrt(sbLeftSEPP * sbLeftSEMM);
    double sbRightSEPPMM = 2 * TMath::Sqrt(sbRightSEPP * sbRightSEMM);

    double normSEPPMM = (sbLeftSEPM + sbRightSEPM) / (sbLeftSEPPMM + sbRightSEPPMM);
    return normSEPPMM;
}
////////////////////////////////////////////////////////////////////////////////
double ComputeSideBandsNorm(TH1D *hist, double minRanges[], double maxRanges[]) {
    double sbLeft = hist -> Integral(hist -> FindBin(minRanges[0]), hist -> FindBin(maxRanges[0]));
    double sbRight = hist -> Integral(hist -> FindBin(minRanges[1]), hist -> FindBin(maxRanges[1]));
    return sbLeft + sbRight;
}
////////////////////////////////////////////////////////////////////////////////
TH1D *CreateSideBandsHist (TH1D *hist, double minRanges[], double maxRanges[], string histName) {
    int nBins = hist -> GetNbinsX();
    double lowEdge = hist -> GetBinLowEdge(1);
    double topEdge =  hist -> GetBinLowEdge(nBins) + hist -> GetBinWidth(nBins);
    TH1D *histSideBands = new TH1D(histName.c_str(), "", nBins, lowEdge, topEdge);

    for (int iBin = 0;iBin < nBins;iBin++) {
        double binCenter = hist -> GetBinCenter(iBin+1);

        if ((binCenter < maxRanges[0] && binCenter > minRanges[0]) || (binCenter < maxRanges[1] && binCenter > minRanges[1])) {
            histSideBands -> SetBinContent(iBin+1, hist -> GetBinContent(iBin+1));
            histSideBands -> SetBinError(iBin+1, hist -> GetBinError(iBin+1));
        }
    }
    return histSideBands;
}
////////////////////////////////////////////////////////////////////////////////
void LoadStyle(){
    int font = 42;
    gStyle -> SetFrameBorderMode(0);
    gStyle -> SetFrameFillColor(0);
    gStyle -> SetCanvasBorderMode(0);
    gStyle -> SetPadBorderMode(0);
    gStyle -> SetPadColor(10);
    gStyle -> SetCanvasColor(10);
    gStyle -> SetTitleFillColor(10);
    gStyle -> SetTitleBorderSize(1);
    gStyle -> SetStatColor(10);
    gStyle -> SetStatBorderSize(1);
    gStyle -> SetLegendBorderSize(1);
    gStyle -> SetDrawBorder(0);
    gStyle -> SetTextFont(font);
    gStyle -> SetStatFontSize(0.05);
    gStyle -> SetStatX(0.97);
    gStyle -> SetStatY(0.98);
    gStyle -> SetStatH(0.03);
    gStyle -> SetStatW(0.3);
    gStyle -> SetTickLength(0.02,"y");
    gStyle -> SetEndErrorSize(3);
    gStyle -> SetLabelSize(0.05,"xyz");
    gStyle -> SetLabelFont(font,"xyz");
    gStyle -> SetLabelOffset(0.01,"xyz");
    gStyle -> SetTitleFont(font,"xyz");
    gStyle -> SetTitleOffset(0.9,"x");
    gStyle -> SetTitleOffset(1.02,"y");
    gStyle -> SetTitleSize(0.05,"xyz");
    gStyle -> SetMarkerSize(1.3);
    gStyle -> SetOptStat(0);
    gStyle -> SetEndErrorSize(0);
    gStyle -> SetCanvasPreferGL(kTRUE);
    gStyle -> SetHatchesSpacing(0.5);
    gStyle -> SetPadLeftMargin(0.15);
    gStyle -> SetPadBottomMargin(0.15);
    gStyle -> SetPadTopMargin(0.05);
    gStyle -> SetPadRightMargin(0.05);
    gStyle -> SetEndErrorSize(0.0);
    gStyle -> SetTitleSize(0.05,"X");
    gStyle -> SetTitleSize(0.045,"Y");
    gStyle -> SetLabelSize(0.045,"X");
    gStyle -> SetLabelSize(0.045,"Y");
    gStyle -> SetTitleOffset(1.2,"X");
    gStyle -> SetTitleOffset(1.35,"Y");
}