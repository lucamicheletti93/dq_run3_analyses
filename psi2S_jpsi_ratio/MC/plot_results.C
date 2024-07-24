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

    const double BrJpsiToMuMu = 0.05961;
    const double errBrJpsiToMuMu = 0.00033;
    const double BrPsi2sToMuMu = 8.0e-3;
    const double errBrPsi2sToMuMu = 0.6e-3;

    double ptBinsRun2[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 26.0, 30.0};
    double widthJpsiRun2[] = {0.067, 0.069, 0.068, 0.069, 0.069, 0.070, 0.071, 0.072, 0.074, 0.075, 0.079, 0.080, 0.087, 0.083, 0.098, 0.083, 0.112, 0.121};
    double errWidthJpsiRun2[] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.002, 0.003, 0.005, 0.007, 0.010, 0.011, 0.019, 0.031};

    TH1D *histWidthJpsiDataRun2 = new TH1D("histWidthJpsiDataRun2","", 18, ptBinsRun2);
    for (int iPt = 0;iPt < 18;iPt++) {
        histWidthJpsiDataRun2 -> SetBinContent(iPt+1, widthJpsiRun2[iPt]);
        histWidthJpsiDataRun2 -> SetBinError(iPt+1, errWidthJpsiRun2[iPt]);
    }
    SetHistogram(histWidthJpsiDataRun2, 1, 1, 20, 1);

    double ptBinsRun3[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0};
    double corrPsi2sOverJpsiRun3Prel[] = {0.114, 0.128, 0.131, 0.203, 0.210, 0.239, 0.241, 0.304};
    double statCorrPsi2sOverJpsiRun3Prel[] = {0.011, 0.009, 0.010, 0.013, 0.017, 0.016, 0.026, 0.037};
    double systCorrPsi2sOverJpsiRun3Prel[] = {0.025, 0.025, 0.026, 0.033, 0.037, 0.045, 0.060, 0.052};

    TH1D *histStatCorrPsi2sOverJpsiRun3Prel = new TH1D("histStatCorrPsi2sOverJpsiRun3Prel", "", 8, ptBinsRun3);
    TH1D *histSystCorrPsi2sOverJpsiRun3Prel = new TH1D("histSystCorrPsi2sOverJpsiRun3Prel", "", 8, ptBinsRun3);

    for (int iPt = 0;iPt < 8;iPt++) {
        histStatCorrPsi2sOverJpsiRun3Prel -> SetBinContent(iPt+1, corrPsi2sOverJpsiRun3Prel[iPt]);
        histStatCorrPsi2sOverJpsiRun3Prel -> SetBinError(iPt+1, statCorrPsi2sOverJpsiRun3Prel[iPt]);
        histSystCorrPsi2sOverJpsiRun3Prel -> SetBinContent(iPt+1, corrPsi2sOverJpsiRun3Prel[iPt]);
        histSystCorrPsi2sOverJpsiRun3Prel -> SetBinError(iPt+1, systCorrPsi2sOverJpsiRun3Prel[iPt]);
    }

    SetHistogram(histStatCorrPsi2sOverJpsiRun3Prel, 862, 1, 20, 1);
    SetHistogram(histSystCorrPsi2sOverJpsiRun3Prel, 862, 1, 20, 1);
    histSystCorrPsi2sOverJpsiRun3Prel -> SetFillStyle(0);

    // Read J/psi signal extraction
    std::ifstream fileJpsi("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/LHC22o_pass6_minBias/time_association/pt_dependence/systematic_sig_Jpsi.txt");

    std::string headerJpsi;
    std::getline(fileJpsi, headerJpsi);
    int counterPt = 0;

    TH1D *histStatJpsiRun3 = new TH1D("histStatJpsiRun3", "", 8, ptBinsRun3);
    TH1D *histSystJpsiRun3 = new TH1D("histSystPsi2sRun3", "", 8, ptBinsRun3);

    for (std::string line; std::getline(fileJpsi, line);) {
        std::istringstream iss(line);
        Row row;
        iss >> row.x_min >> row.x_max >> row.val >> row.stat >> row.syst;
        histStatJpsiRun3 -> SetBinContent(counterPt+1, row.val);
        histStatJpsiRun3 -> SetBinError(counterPt+1, row.stat);
        histSystJpsiRun3 -> SetBinContent(counterPt+1, row.val);
        histSystJpsiRun3 -> SetBinError(counterPt+1, row.syst);
        counterPt++;
    }

    SetHistogram(histStatJpsiRun3, 633, 1, 20, 1);
    SetHistogram(histSystJpsiRun3, 633, 1, 20, 1);
    histSystJpsiRun3 -> SetFillStyle(0);

    // Read Psi(2S) signal extraction
    std::ifstream filePsi2s("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/LHC22o_pass6_minBias/time_association/pt_dependence/systematic_sig_Psi2s.txt");

    std::string headerPsi2s;
    std::getline(filePsi2s, headerPsi2s);
    counterPt = 0;

    TH1D *histStatPsi2sRun3 = new TH1D("histStatPsi2sRun3", "", 8, ptBinsRun3);
    TH1D *histSystPsi2sRun3 = new TH1D("histSystPsi2sRun3", "", 8, ptBinsRun3);

    for (std::string line; std::getline(filePsi2s, line);) {
        std::istringstream iss(line);
        Row row;
        iss >> row.x_min >> row.x_max >> row.val >> row.stat >> row.syst;
        histStatPsi2sRun3 -> SetBinContent(counterPt+1, row.val);
        histStatPsi2sRun3 -> SetBinError(counterPt+1, row.stat);
        histSystPsi2sRun3 -> SetBinContent(counterPt+1, row.val);
        histSystPsi2sRun3 -> SetBinError(counterPt+1, row.syst);
        counterPt++;
    }

    SetHistogram(histStatPsi2sRun3, 633, 1, 20, 1);
    SetHistogram(histSystPsi2sRun3, 633, 1, 20, 1);
    histSystPsi2sRun3 -> SetFillStyle(0);

    TFile *fInDataTimeAssocCB2 = new TFile("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/LHC22o_pass6_minBias/time_association/MC_tails/myAnalysis_signal_extraction_CB2.root");
    TH1D *histSigJpsiDataTimeAssocCB2 = (TH1D*) fInDataTimeAssocCB2 -> Get("hist_sig_Jpsi");
    TH1D *histSigPsi2sDataTimeAssocCB2 = (TH1D*) fInDataTimeAssocCB2 -> Get("hist_sig_Psi2s");
    TH1D *histWidthJpsiDataTimeAssocCB2 = (TH1D*) fInDataTimeAssocCB2 -> Get("hist_width_Jpsi");
    TH1D *histMeanJpsiDataTimeAssocCB2 = (TH1D*) fInDataTimeAssocCB2 -> Get("hist_mean_Jpsi");

    SetHistogram(histSigJpsiDataTimeAssocCB2, 1, 1, 24, 1);
    SetHistogram(histSigPsi2sDataTimeAssocCB2, 1, 1, 24, 1);
    SetHistogram(histWidthJpsiDataTimeAssocCB2, 1, 1, 24, 1);
    SetHistogram(histMeanJpsiDataTimeAssocCB2, 1, 1, 24, 1);

    TH1D *histPsi2sOverJpsi = (TH1D*) histSigPsi2sDataTimeAssocCB2 -> Clone("histPsi2sOverJpsi");
    histPsi2sOverJpsi -> Divide(histSigJpsiDataTimeAssocCB2);

    TFile *fInAxeTimeAssoc = new TFile("/Users/lucamicheletti/GITHUB/dq_run3_analyses/psi2S_jpsi_ratio/MC/Axe_LHC24e5_time_association.root");
    TH1D *histAxePsi2sOverJpsi = (TH1D*) fInAxeTimeAssoc -> Get("histRatioAxePtCut1");

    TH1D *histCorrPsi2sOverJpsi = (TH1D*) histPsi2sOverJpsi -> Clone("histCorrPsi2sOverJpsi");
    histCorrPsi2sOverJpsi -> Divide(histAxePsi2sOverJpsi);
    histCorrPsi2sOverJpsi -> Scale(BrJpsiToMuMu / BrPsi2sToMuMu);

    // QA checks
    TFile *fInTimeAssocCB2 = new TFile("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/LHC24e5/time_association/multi_trial/myAnalysis_width_Jpsi_CB2.root");
    TH1D *histWidthJpsiTimeAssocCB2 = (TH1D*) fInTimeAssocCB2 -> Get("hist_width_Jpsi");
    SetHistogram(histWidthJpsiTimeAssocCB2, 633, 1, 20, 0.8);

    TFile *fInStdAssocCB2 = new TFile("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/LHC24e5/standard_association/multi_trial/myAnalysis_width_Jpsi_CB2.root");
    TH1D *histWidthJpsiStdAssocCB2 = (TH1D*) fInStdAssocCB2 -> Get("hist_width_Jpsi");
    SetHistogram(histWidthJpsiStdAssocCB2, 862, 1, 20, 0.8);

    TFile *fInTimeAssocNA60 = new TFile("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/LHC24e5/time_association/multi_trial/myAnalysis_width_Jpsi_NA60.root");
    TH1D *histWidthJpsiTimeAssocNA60 = (TH1D*) fInTimeAssocNA60 -> Get("hist_width_Jpsi");
    SetHistogram(histWidthJpsiTimeAssocNA60, 417, 1, 20, 0.8);

    TFile *fInStdAssocNA60 = new TFile("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/LHC24e5/standard_association/multi_trial/myAnalysis_width_Jpsi_NA60.root");
    TH1D *histWidthJpsiStdAssocNA60 = (TH1D*) fInStdAssocNA60 -> Get("hist_width_Jpsi");
    SetHistogram(histWidthJpsiStdAssocNA60, 880, 1, 20, 0.8);

    TCanvas *canvasWidthJpsi = new TCanvas("canvasWidthJpsi", "", 800, 600);
    histWidthJpsiTimeAssocCB2 -> SetTitle("");
    histWidthJpsiTimeAssocCB2 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histWidthJpsiTimeAssocCB2 -> GetYaxis() -> SetTitle("#sigma_{J/#psi} (GeV/#it{c}^{2})");
    histWidthJpsiTimeAssocCB2 -> GetYaxis() -> SetRangeUser(0, 0.2);
    histWidthJpsiTimeAssocCB2 -> Draw("EP");
    histWidthJpsiStdAssocCB2 -> Draw("EP SAME");
    histWidthJpsiTimeAssocNA60 -> Draw("EP SAME");
    histWidthJpsiStdAssocNA60 -> Draw("EP SAME");
    histWidthJpsiDataTimeAssocCB2 -> Draw("EP SAME");
    histWidthJpsiDataRun2 -> Draw("EP SAME");

    TFile *fOut = new TFile("jpsi_width.root", "RECREATE");
    histWidthJpsiDataTimeAssocCB2 -> Write("histJpsiWidthRun3");
    histWidthJpsiDataRun2 -> Write("histJpsiWidthRun2");

    TLegend *legendWidthJpsi = new TLegend(0.20, 0.69, 0.60, 0.89);
    SetLegend(legendWidthJpsi);
    legendWidthJpsi -> AddEntry(histWidthJpsiTimeAssocCB2, "MC, time association [CB2]", "PL");
    legendWidthJpsi -> AddEntry(histWidthJpsiStdAssocCB2, "MC, std. association [CB2]", "PL");
    legendWidthJpsi -> AddEntry(histWidthJpsiTimeAssocNA60, "MC, time association [NA60]", "PL");
    legendWidthJpsi -> AddEntry(histWidthJpsiStdAssocNA60, "MC, std. association [NA60]", "PL");
    legendWidthJpsi -> AddEntry(histWidthJpsiDataTimeAssocCB2, "Data, time association", "PL");
    legendWidthJpsi -> Draw("SAME");

    canvasWidthJpsi -> SaveAs("plots/width_Jpsi_pt.pdf");


    TCanvas *canvasSigJpsi = new TCanvas("canvasSigJpsi", "", 800, 600);
    histSigJpsiDataTimeAssocCB2 -> SetTitle("");
    histSigJpsiDataTimeAssocCB2 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histSigJpsiDataTimeAssocCB2 -> GetYaxis() -> SetTitle("#it{N}_{J/#psi}");
    histSigJpsiDataTimeAssocCB2 -> Draw("EP");
    histStatJpsiRun3 -> Draw("EP SAME");
    histSystJpsiRun3 -> Draw("E2P SAME");

    TCanvas *canvasSigPsi2s = new TCanvas("canvasSigPsi2s", "", 800, 600);
    histSigPsi2sDataTimeAssocCB2 -> SetTitle("");
    histSigPsi2sDataTimeAssocCB2 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histSigPsi2sDataTimeAssocCB2 -> GetYaxis() -> SetTitle("#it{N}_{#psi(2S)}");
    histSigPsi2sDataTimeAssocCB2 -> Draw("EP");
    histStatPsi2sRun3 -> Draw("EP SAME");
    histSystPsi2sRun3 -> Draw("E2P SAME");

    TCanvas *canvasPsi2sOverJpsi = new TCanvas("canvasPsi2sOverJpsi", "", 800, 600);
    histPsi2sOverJpsi -> SetTitle("");
    histPsi2sOverJpsi -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histPsi2sOverJpsi -> GetYaxis() -> SetTitle("#it{N}_{#psi(2S)} / #it{N}_{J/#psi}");
    histPsi2sOverJpsi -> Draw("EP");

    TCanvas *canvasAxePsi2sOverJpsi = new TCanvas("canvasAxePsi2sOverJpsi", "", 800, 600);
    histAxePsi2sOverJpsi -> SetTitle("");
    histAxePsi2sOverJpsi -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histAxePsi2sOverJpsi -> GetYaxis() -> SetTitle("A#times#varepsilon_{#psi(2S)} / A#times#varepsilon_{J/#psi}");
    histAxePsi2sOverJpsi -> Draw("EP");

    TCanvas *canvasCorrPsi2sOverJpsi = new TCanvas("canvasCorrPsi2sOverJpsi", "", 800, 600);
    histCorrPsi2sOverJpsi -> SetTitle("");
    histCorrPsi2sOverJpsi -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histCorrPsi2sOverJpsi -> GetYaxis() -> SetTitle("d#sigma^{#psi(2S)}/d#it{p}_{T} / d#sigma^{J/#psi}/d#it{p}_{T}");
    histCorrPsi2sOverJpsi -> Draw("EP");
    histStatCorrPsi2sOverJpsiRun3Prel -> Draw("EP SAME");
    histSystCorrPsi2sOverJpsiRun3Prel -> Draw("E2P SAME");









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