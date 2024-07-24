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
void SetHistogram(TH1D *, int , int , int , double );

struct Row {
    double x_min, x_max, val, stat, syst;
};

void plot_results() {
    LoadStyle();

    //string productionName = "LHC22o_pass6_minBias";
    string productionName = "LHC22o_pass7_skimmed";

    const double BrJpsiToMuMu = 0.05961;
    const double errBrJpsiToMuMu = 0.00033;
    const double BrPsi2sToMuMu = 8.0e-3;
    const double errBrPsi2sToMuMu = 0.6e-3;

    //const int nPtBinsRun3 = 14;
    //double ptBinsRun3[] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 10, 15, 20};
    const int nPtBinsRun3 = 13;
    double ptBinsRun3[] = {0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20};

    const int nPtBinsRun3Prel = 8;
    double ptBinsRun3Prel[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0};

    double corrPsi2sOverJpsiVsPtRun3Prel[] = {0.114, 0.128, 0.131, 0.203, 0.210, 0.239, 0.241, 0.304};
    double statCorrPsi2sOverJpsiVsPtRun3Prel[] = {0.011, 0.009, 0.010, 0.013, 0.017, 0.016, 0.026, 0.037};
    double systCorrPsi2sOverJpsiVsPtRun3Prel[] = {0.025, 0.025, 0.026, 0.033, 0.037, 0.045, 0.060, 0.052};

    TH1D *histStatCorrPsi2sOverJpsiVsPtRun3Prel = new TH1D("histStatCorrPsi2sOverJpsiVsPtRun3Prel", "", nPtBinsRun3Prel, ptBinsRun3Prel);
    TH1D *histSystCorrPsi2sOverJpsiVsPtRun3Prel = new TH1D("histSystCorrPsi2sOverJpsiVsPtRun3Prel", "", nPtBinsRun3Prel, ptBinsRun3Prel);

    for (int iPt = 0;iPt < nPtBinsRun3;iPt++) {
        histStatCorrPsi2sOverJpsiVsPtRun3Prel -> SetBinContent(iPt+1, corrPsi2sOverJpsiVsPtRun3Prel[iPt]);
        histStatCorrPsi2sOverJpsiVsPtRun3Prel -> SetBinError(iPt+1, statCorrPsi2sOverJpsiVsPtRun3Prel[iPt]);
        histSystCorrPsi2sOverJpsiVsPtRun3Prel -> SetBinContent(iPt+1, corrPsi2sOverJpsiVsPtRun3Prel[iPt]);
        histSystCorrPsi2sOverJpsiVsPtRun3Prel -> SetBinError(iPt+1, systCorrPsi2sOverJpsiVsPtRun3Prel[iPt]);
    }

    SetHistogram(histStatCorrPsi2sOverJpsiVsPtRun3Prel, 922, 1, 24, 1);
    SetHistogram(histSystCorrPsi2sOverJpsiVsPtRun3Prel, 922, 1, 24, 1);
    histSystCorrPsi2sOverJpsiVsPtRun3Prel -> SetFillStyle(0);


    double rapBinsRun3[] = {2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00};
    double corrPsi2sOverJpsiVsRapRun3Prel[] = {0.133, 0.161, 0.154, 0.149, 0.136, 0.101};
    double statCorrPsi2sOverJpsiVsRapRun3Prel[] = {0.021, 0.012, 0.008, 0.009, 0.009, 0.017};
    double systCorrPsi2sOverJpsiVsRapRun3Prel[] = {0.028, 0.028, 0.029, 0.026, 0.033, 0.043};

    TH1D *histStatCorrPsi2sOverJpsiVsRapRun3Prel = new TH1D("histStatCorrPsi2sOverJpsiVsRapRun3Prel", "", 6, rapBinsRun3);
    TH1D *histSystCorrPsi2sOverJpsiVsRapRun3Prel = new TH1D("histSystCorrPsi2sOverJpsiVsRapRun3Prel", "", 6, rapBinsRun3);

    for (int iRap = 0;iRap < 6;iRap++) {
        histStatCorrPsi2sOverJpsiVsRapRun3Prel -> SetBinContent(iRap+1, corrPsi2sOverJpsiVsRapRun3Prel[iRap]);
        histStatCorrPsi2sOverJpsiVsRapRun3Prel -> SetBinError(iRap+1, statCorrPsi2sOverJpsiVsRapRun3Prel[iRap]);
        histSystCorrPsi2sOverJpsiVsRapRun3Prel -> SetBinContent(iRap+1, corrPsi2sOverJpsiVsRapRun3Prel[iRap]);
        histSystCorrPsi2sOverJpsiVsRapRun3Prel -> SetBinError(iRap+1, systCorrPsi2sOverJpsiVsRapRun3Prel[iRap]);
    }

    SetHistogram(histStatCorrPsi2sOverJpsiVsRapRun3Prel, 922, 1, 24, 1);
    SetHistogram(histSystCorrPsi2sOverJpsiVsRapRun3Prel, 922, 1, 24, 1);
    histSystCorrPsi2sOverJpsiVsRapRun3Prel -> SetFillStyle(0);

    // Read J/psi signal extraction
    std::ifstream fileJpsiVsPt(Form("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/%s/time_association/pt_dependence/systematic_sig_Jpsi.txt", productionName.c_str()));

    std::string headerJpsiVsPt;
    std::getline(fileJpsiVsPt, headerJpsiVsPt);
    int counterPt = 0;

    TH1D *histStatJpsiVsPtRun3 = new TH1D("histStatJpsiVsPtRun3", "", nPtBinsRun3, ptBinsRun3);
    TH1D *histSystJpsiVsPtRun3 = new TH1D("histSystPsi2sVsPtRun3", "", nPtBinsRun3, ptBinsRun3);

    for (std::string line; std::getline(fileJpsiVsPt, line);) {
        std::istringstream iss(line);
        Row row;
        iss >> row.x_min >> row.x_max >> row.val >> row.stat >> row.syst;
        histStatJpsiVsPtRun3 -> SetBinContent(counterPt+1, row.val);
        histStatJpsiVsPtRun3 -> SetBinError(counterPt+1, row.stat);
        histSystJpsiVsPtRun3 -> SetBinContent(counterPt+1, row.val);
        histSystJpsiVsPtRun3 -> SetBinError(counterPt+1, row.syst);
        counterPt++;
    }

    SetHistogram(histStatJpsiVsPtRun3, 633, 1, 20, 1);
    SetHistogram(histSystJpsiVsPtRun3, 633, 1, 20, 1);
    histSystJpsiVsPtRun3 -> SetFillStyle(0);


    std::ifstream fileJpsiVsRap("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/LHC22o_pass6_minBias/time_association/rap_dependence/systematic_sig_Jpsi.txt");

    std::string headerJpsiVsRap;
    std::getline(fileJpsiVsRap, headerJpsiVsRap);
    int counterRap = 0;

    TH1D *histStatJpsiVsRapRun3 = new TH1D("histStatJpsiVsRapRun3", "", 6, rapBinsRun3);
    TH1D *histSystJpsiVsRapRun3 = new TH1D("histSystPsi2sVsRapRun3", "", 6, rapBinsRun3);

    for (std::string line; std::getline(fileJpsiVsRap, line);) {
        std::istringstream iss(line);
        Row row;
        iss >> row.x_min >> row.x_max >> row.val >> row.stat >> row.syst;
        histStatJpsiVsRapRun3 -> SetBinContent(counterRap+1, row.val);
        histStatJpsiVsRapRun3 -> SetBinError(counterRap+1, row.stat);
        histSystJpsiVsRapRun3 -> SetBinContent(counterRap+1, row.val);
        histSystJpsiVsRapRun3 -> SetBinError(counterRap+1, row.syst);
        counterRap++;
    }

    SetHistogram(histStatJpsiVsRapRun3, 633, 1, 20, 1);
    SetHistogram(histSystJpsiVsRapRun3, 633, 1, 20, 1);
    histSystJpsiVsRapRun3 -> SetFillStyle(0);


    // Read Psi(2S) signal extraction
    std::ifstream filePsi2sVsPt(Form("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/%s/time_association/pt_dependence/systematic_sig_Psi2s.txt", productionName.c_str()));

    std::string headerPsi2sVsPt;
    std::getline(filePsi2sVsPt, headerPsi2sVsPt);
    counterPt = 0;

    TH1D *histStatPsi2sVsPtRun3 = new TH1D("histStatPsi2sVsPtRun3", "", nPtBinsRun3, ptBinsRun3);
    TH1D *histSystPsi2sVsPtRun3 = new TH1D("histSystPsi2sVsPtRun3", "", nPtBinsRun3, ptBinsRun3);

    for (std::string line; std::getline(filePsi2sVsPt, line);) {
        std::istringstream iss(line);
        Row row;
        iss >> row.x_min >> row.x_max >> row.val >> row.stat >> row.syst;
        histStatPsi2sVsPtRun3 -> SetBinContent(counterPt+1, row.val);
        histStatPsi2sVsPtRun3 -> SetBinError(counterPt+1, row.stat);
        histSystPsi2sVsPtRun3 -> SetBinContent(counterPt+1, row.val);
        histSystPsi2sVsPtRun3 -> SetBinError(counterPt+1, row.syst);
        counterPt++;
    }

    SetHistogram(histStatPsi2sVsPtRun3, 633, 1, 20, 1);
    SetHistogram(histSystPsi2sVsPtRun3, 633, 1, 20, 1);
    histSystPsi2sVsPtRun3 -> SetFillStyle(0);


    std::ifstream filePsi2sVsRap("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/LHC22o_pass6_minBias/time_association/rap_dependence/systematic_sig_Psi2s.txt");

    std::string headerPsi2sVsRap;
    std::getline(filePsi2sVsRap, headerPsi2sVsRap);
    counterRap = 0;

    TH1D *histStatPsi2sVsRapRun3 = new TH1D("histStatPsi2sVsRapRun3", "", 6, rapBinsRun3);
    TH1D *histSystPsi2sVsRapRun3 = new TH1D("histSystPsi2sVsRapRun3", "", 6, rapBinsRun3);

    for (std::string line; std::getline(filePsi2sVsRap, line);) {
        std::istringstream iss(line);
        Row row;
        iss >> row.x_min >> row.x_max >> row.val >> row.stat >> row.syst;
        histStatPsi2sVsRapRun3 -> SetBinContent(counterRap+1, row.val);
        histStatPsi2sVsRapRun3 -> SetBinError(counterRap+1, row.stat);
        histSystPsi2sVsRapRun3 -> SetBinContent(counterRap+1, row.val);
        histSystPsi2sVsRapRun3 -> SetBinError(counterRap+1, row.syst);
        counterRap++;
    }

    SetHistogram(histStatPsi2sVsRapRun3, 633, 1, 20, 1);
    SetHistogram(histSystPsi2sVsRapRun3, 633, 1, 20, 1);
    histSystPsi2sVsRapRun3 -> SetFillStyle(0);


    TH1D *histStatPsi2sOverJpsiVsPtRun3 = (TH1D*) histStatPsi2sVsPtRun3 -> Clone("histStatPsi2sOverJpsiVsPtRun3");
    histStatPsi2sOverJpsiVsPtRun3 -> Divide(histStatJpsiVsPtRun3);

    TH1D *histSystPsi2sOverJpsiVsPtRun3 = (TH1D*) histSystPsi2sVsPtRun3 -> Clone("histSystPsi2sOverJpsiVsPtRun3");
    histSystPsi2sOverJpsiVsPtRun3 -> Divide(histSystJpsiVsPtRun3);

    TH1D *histStatPsi2sOverJpsiVsRapRun3 = (TH1D*) histStatPsi2sVsRapRun3 -> Clone("histStatPsi2sOverJpsiVsRapRun3");
    histStatPsi2sOverJpsiVsRapRun3 -> Divide(histStatJpsiVsRapRun3);

    TH1D *histSystPsi2sOverJpsiVsRapRun3 = (TH1D*) histSystPsi2sVsRapRun3 -> Clone("histSystPsi2sOverJpsiVsRapRun3");
    histSystPsi2sOverJpsiVsRapRun3 -> Divide(histSystJpsiVsRapRun3);

    TFile *fInAxeTimeAssoc = new TFile("/Users/lucamicheletti/GITHUB/dq_run3_analyses/psi2S_jpsi_ratio/MC/Axe_LHC24e5_time_association.root");
    TH1D *histAxePsi2sOverJpsiVsPt = (TH1D*) fInAxeTimeAssoc -> Get("histRatioAxePtCut1");
    TH1D *histAxePsi2sOverJpsiVsRap = (TH1D*) fInAxeTimeAssoc -> Get("histRatioAxeRapCut1");

    TH1D *histStatCorrPsi2sOverJpsiVsPt = (TH1D*) histStatPsi2sOverJpsiVsPtRun3 -> Clone("histStatCorrPsi2sOverJpsiVsPt");
    histStatCorrPsi2sOverJpsiVsPt -> Divide(histAxePsi2sOverJpsiVsPt);
    histStatCorrPsi2sOverJpsiVsPt -> Scale(BrJpsiToMuMu / BrPsi2sToMuMu);

    TH1D *histSystCorrPsi2sOverJpsiVsPt = (TH1D*) histSystPsi2sOverJpsiVsPtRun3 -> Clone("histSystCorrPsi2sOverJpsiVsPt");
    histSystCorrPsi2sOverJpsiVsPt -> Divide(histAxePsi2sOverJpsiVsPt);
    histSystCorrPsi2sOverJpsiVsPt -> Scale(BrJpsiToMuMu / BrPsi2sToMuMu);

    TH1D *histStatCorrPsi2sOverJpsiVsRap = (TH1D*) histStatPsi2sOverJpsiVsRapRun3 -> Clone("histStatCorrPsi2sOverJpsiVsRap");
    histStatCorrPsi2sOverJpsiVsRap -> Divide(histAxePsi2sOverJpsiVsRap);
    histStatCorrPsi2sOverJpsiVsRap -> Scale(BrJpsiToMuMu / BrPsi2sToMuMu);

    TH1D *histSystCorrPsi2sOverJpsiVsRap = (TH1D*) histSystPsi2sOverJpsiVsRapRun3 -> Clone("histSystCorrPsi2sOverJpsiVsRap");
    histSystCorrPsi2sOverJpsiVsRap -> Divide(histAxePsi2sOverJpsiVsRap);
    histSystCorrPsi2sOverJpsiVsRap -> Scale(BrJpsiToMuMu / BrPsi2sToMuMu);

    TCanvas *canvasSigJpsiVsPt = new TCanvas("canvasSigJpsiVsPt", "", 800, 600);
    gPad -> SetLogy(1);
    histStatJpsiVsPtRun3 -> SetTitle("");
    histStatJpsiVsPtRun3 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histStatJpsiVsPtRun3 -> GetYaxis() -> SetTitle("#it{N}_{J/#psi}");
    histStatJpsiVsPtRun3 -> Scale(1, "WIDTH");
    histSystJpsiVsPtRun3 -> Scale(1, "WIDTH");
    histStatJpsiVsPtRun3 -> Draw("EP SAME");
    histSystJpsiVsPtRun3 -> Draw("E2P SAME");

    TCanvas *canvasSigJpsiVsRap = new TCanvas("canvasSigJpsiVsRap", "", 800, 600);
    gPad -> SetLogy(1);
    histStatJpsiVsRapRun3 -> SetTitle("");
    histStatJpsiVsRapRun3 -> GetXaxis() -> SetTitle("#it{y}");
    histStatJpsiVsRapRun3 -> GetYaxis() -> SetTitle("#it{N}_{J/#psi}");
    histStatJpsiVsRapRun3 -> Draw("EP SAME");
    histSystJpsiVsRapRun3 -> Draw("E2P SAME");

    TCanvas *canvasSigPsi2sVsPt = new TCanvas("canvasSigPsi2sVsPt", "", 800, 600);
    gPad -> SetLogy(1);
    histStatPsi2sVsPtRun3 -> SetTitle("");
    histStatPsi2sVsPtRun3 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histStatPsi2sVsPtRun3 -> GetYaxis() -> SetTitle("#it{N}_{#psi(2S)}");
    histStatPsi2sVsPtRun3 -> Draw("EP");
    histStatPsi2sVsPtRun3 -> Scale(1, "WIDTH");
    histSystPsi2sVsPtRun3 -> Scale(1, "WIDTH");
    histStatPsi2sVsPtRun3 -> Draw("EP SAME");
    histSystPsi2sVsPtRun3 -> Draw("E2P SAME");

    TCanvas *canvasSigPsi2sVsRap = new TCanvas("canvasSigPsi2sVsRap", "", 800, 600);
    gPad -> SetLogy(1);
    histStatPsi2sVsRapRun3 -> SetTitle("");
    histStatPsi2sVsRapRun3 -> GetXaxis() -> SetTitle("#it{y}");
    histStatPsi2sVsRapRun3 -> GetYaxis() -> SetTitle("#it{N}_{#psi(2S)}");
    histStatPsi2sVsRapRun3 -> Draw("EP");
    histStatPsi2sVsRapRun3 -> Draw("EP SAME");
    histSystPsi2sVsRapRun3 -> Draw("E2P SAME");

    TCanvas *canvasPsi2sOverJpsiVsPt = new TCanvas("canvasPsi2sOverJpsiVsPt", "", 800, 600);
    histStatPsi2sOverJpsiVsPtRun3 -> SetTitle("");
    histStatPsi2sOverJpsiVsPtRun3 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histStatPsi2sOverJpsiVsPtRun3 -> GetYaxis() -> SetTitle("#it{N}_{#psi(2S)} / #it{N}_{J/#psi}");
    histStatPsi2sOverJpsiVsPtRun3 -> GetYaxis() -> SetRangeUser(0, 0.1);
    histStatPsi2sOverJpsiVsPtRun3 -> Draw("EP");
    histSystPsi2sOverJpsiVsPtRun3 -> Draw("E2P SAME");

    TCanvas *canvasPsi2sOverJpsiVsRap = new TCanvas("canvasPsi2sOverJpsiVsRap", "", 800, 600);
    histStatPsi2sOverJpsiVsRapRun3 -> SetTitle("");
    histStatPsi2sOverJpsiVsRapRun3 -> GetXaxis() -> SetTitle("#it{y}");
    histStatPsi2sOverJpsiVsRapRun3 -> GetYaxis() -> SetTitle("#it{N}_{#psi(2S)} / #it{N}_{J/#psi}");
    histStatPsi2sOverJpsiVsRapRun3 -> Draw("EP");
    histSystPsi2sOverJpsiVsRapRun3 -> Draw("E2P SAME");

    TCanvas *canvasAxePsi2sOverJpsiVsPt = new TCanvas("canvasAxePsi2sOverJpsiVsPt", "", 800, 600);
    histAxePsi2sOverJpsiVsPt -> SetTitle("");
    histAxePsi2sOverJpsiVsPt -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histAxePsi2sOverJpsiVsPt -> GetYaxis() -> SetTitle("A#times#varepsilon_{#psi(2S)} / A#times#varepsilon_{J/#psi}");
    histAxePsi2sOverJpsiVsPt -> Draw("EP");

    TCanvas *canvasCorrPsi2sOverJpsiVsPt = new TCanvas("canvasCorrPsi2sOverJpsiVsPt", "", 800, 600);
    histStatCorrPsi2sOverJpsiVsPt -> SetTitle("");
    histStatCorrPsi2sOverJpsiVsPt -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histStatCorrPsi2sOverJpsiVsPt -> GetYaxis() -> SetTitle("d#sigma^{#psi(2S)}/d#it{p}_{T} / d#sigma^{J/#psi}/d#it{p}_{T}");
    histStatCorrPsi2sOverJpsiVsPt -> GetYaxis() -> SetRangeUser(0, 0.5);
    histStatCorrPsi2sOverJpsiVsPt -> Draw("EP");
    histSystCorrPsi2sOverJpsiVsPt -> Draw("E2P SAME");
    histStatCorrPsi2sOverJpsiVsPtRun3Prel -> Draw("EP SAME");
    histSystCorrPsi2sOverJpsiVsPtRun3Prel -> Draw("E2P SAME");

    TCanvas *canvasCorrPsi2sOverJpsiVsRap = new TCanvas("canvasCorrPsi2sOverJpsiVsRap", "", 800, 600);
    histStatCorrPsi2sOverJpsiVsRap -> SetTitle("");
    histStatCorrPsi2sOverJpsiVsRap -> GetXaxis() -> SetTitle("#it{y}");
    histStatCorrPsi2sOverJpsiVsRap -> GetYaxis() -> SetTitle("d#sigma^{#psi(2S)}/d#it{y} / d#sigma^{J/#psi}/d#it{y}");
    histStatCorrPsi2sOverJpsiVsRap -> GetYaxis() -> SetRangeUser(0, 0.5);
    histStatCorrPsi2sOverJpsiVsRap -> Draw("EP");
    histSystCorrPsi2sOverJpsiVsRap -> Draw("E2P SAME");
    histStatCorrPsi2sOverJpsiVsRapRun3Prel -> Draw("EP SAME");
    histSystCorrPsi2sOverJpsiVsRapRun3Prel -> Draw("E2P SAME");


    canvasSigJpsiVsPt -> SaveAs("plots/SigJpsiVsPt.pdf");
    canvasSigJpsiVsRap -> SaveAs("plots/SigJpsiVsRap.pdf");
    canvasSigPsi2sVsPt -> SaveAs("plots/SigPsi2sVsPt.pdf");
    canvasSigPsi2sVsRap -> SaveAs("plots/SigPsi2sVsRap.pdf");
    canvasPsi2sOverJpsiVsPt -> SaveAs("plots/Psi2sOverJpsiVsPt.pdf");
    canvasPsi2sOverJpsiVsRap -> SaveAs("plots/Psi2sOverJpsiVsRap.pdf");
    canvasCorrPsi2sOverJpsiVsPt -> SaveAs("plots/CorrPsi2sOverJpsiVsPt.pdf");
    canvasCorrPsi2sOverJpsiVsRap -> SaveAs("plots/CorrPsi2sOverJpsiVsRap.pdf");
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
void SetHistogram(TH1D *hist, int color, int lineWidth, int markerStyle, double markerSize) {
    hist -> SetLineColor(color);
    hist -> SetLineWidth(lineWidth);
    hist -> SetMarkerColor(color);
    hist -> SetMarkerStyle(markerStyle);
    hist -> SetMarkerSize(markerSize);
}