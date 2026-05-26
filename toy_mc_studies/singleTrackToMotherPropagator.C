#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <string>

#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TRandom3.h>

#include "Pythia8/Pythia.h"

using namespace Pythia8;

double FuncPt(double *, double *);
double FuncEta(double *, double *);
double FuncAxe(double *, double *);

void LoadStyle();
void SetLegend(TLegend *);
inline void SetHist(auto *hist, Color_t mkrCol = kBlack, int mkrSty = 20, double mkrSize = 1, Color_t lnCol = kBlack, int lnWidth = 1, int fillSty = 0, double alpha = 1) {
    hist -> SetMarkerColorAlpha(mkrCol, alpha);
    hist -> SetMarkerStyle(mkrSty);
    hist -> SetMarkerSize(mkrSize);
    hist -> SetLineColorAlpha(lnCol, alpha);
    hist -> SetLineWidth(lnWidth);
    hist -> SetFillStyle(fillSty);
}

void singleTrackToMotherPropagator(int nEvts = 1e7, int pdgCodeMom = 443, int pdgCodeDau = 13) {
    LoadStyle();

    Pythia8::Pythia pythia;
    pythia.readString("ProcessLevel:all = off");
    pythia.readString("HadronLevel:Decay = on");
    pythia.readString(Form("%d:mayDecay = on", pdgCodeMom));
    pythia.readString(Form("%d:oneChannel = 1 1.0 0 13 -13", pdgCodeMom)); // BR = 100% in mu+ mu-
    pythia.init();

    double massMom = pythia.particleData.m0(pdgCodeMom);
    double massDau = pythia.particleData.m0(pdgCodeDau);

    string momName = pdgCodeMom == 443 ? "J/#psi" : "#psi(2S)";

    std::cout << "*** Get input shapes ***" << std::endl;
    TFile *fInPtMom = new TFile(Form("xsec_%d_pt_pp13TeV.root", pdgCodeMom), "READ");
    TGraphAsymmErrors *graPtMom = (TGraphAsymmErrors*) fInPtMom -> Get(pdgCodeMom == 443 ? "Table 1/Graph1D_y1" : "Table 3/Graph1D_y1");
    TFile *fInEtaMom = new TFile(Form("xsec_%d_y_pp13TeV.root", pdgCodeMom), "READ");
    TGraphAsymmErrors *graEtaMom = (TGraphAsymmErrors*) fInEtaMom -> Get(pdgCodeMom == 443 ? "Table 2/Graph1D_y1" : "Table 4/Graph1D_y1");

    TF1 *funcPtJpsi = new TF1("funcPtJpsi", FuncPt, 0., 30., 4);
    funcPtJpsi -> SetParameter(0, 1000);
    funcPtJpsi -> SetParameter(1, 4.75208);
    funcPtJpsi -> SetParameter(2, 1.69247);
    funcPtJpsi -> SetParameter(3, 4.49224);
    graPtMom -> Fit(funcPtJpsi, "R");

    TF1 *funcEtaJpsi = new TF1("funcEtaJpsi", FuncEta, -5.0, 5.0, 3);
    funcEtaJpsi -> SetParameter(0, 1000);
    funcEtaJpsi -> FixParameter(1, 0);
    funcEtaJpsi -> SetParameter(2, 2.98887);
    graEtaMom -> Fit(funcEtaJpsi, "R");

    TF1 *funcPhiJpsi = new TF1("funcPhiJpsi", "[0]", 0., 2. * TMath::Pi());
    funcPhiJpsi -> SetParameter(0, 1);


    /*TFile *fInPsi2sPt = new TFile("xsec_psi2s_pt_pp13TeV.root", "READ");
    TGraphAsymmErrors *graPtPsi2s = (TGraphAsymmErrors*) fInPsi2sPt -> Get("Table 3/Graph1D_y1");
    TFile *fInPsi2sEta = new TFile("xsec_psi2s_y_pp13TeV.root", "READ");
    TGraphAsymmErrors *graEtaPsi2s = (TGraphAsymmErrors*) fInPsi2sEta -> Get("Table 4/Graph1D_y1");

    TF1 *funcPtPsi2s = new TF1("funcPtPsi2s", FuncPt, 0., 30., 4);
    funcPtPsi2s -> SetParameter(0, 1000);
    funcPtPsi2s -> SetParameter(1, 4.75208);
    funcPtPsi2s -> SetParameter(2, 1.69247);
    funcPtPsi2s -> SetParameter(3, 4.49224);
    graPtPsi2s -> Fit(funcPtPsi2s, "R");

    TF1 *funcEtaPsi2s = new TF1("funcEtaPsi2s", FuncEta, -5.0, 5.0, 3);
    funcEtaPsi2s -> SetParameter(0, 1000);
    funcEtaPsi2s -> FixParameter(1, 0);
    funcEtaPsi2s -> SetParameter(2, 2.98887);
    graEtaPsi2s -> Fit(funcEtaPsi2s, "R");

    TF1 *funcPhiPsi2s = new TF1("funcPhiPsi2s", "[0]", 0., 2. * TMath::Pi());
    funcPhiPsi2s -> SetParameter(0, 1);*/


    TCanvas *canvasInputShapes = new TCanvas("canvasInputShapes", "", 1200, 600);
    canvasInputShapes -> Divide(2, 1);

    canvasInputShapes -> cd(1);
    gPad -> SetLogy(true);
    TH2D *histGridPtJpsi = new TH2D("histGridPtJpsi", ";#it{p}_{T} (GeV/#it{c});Counts", 100, 0, 30, 100, 0.1, 3e3);
    histGridPtJpsi -> Draw();
    graPtMom -> Draw("EP SAME"); 
    funcPtJpsi -> Draw("SAME");
    
    canvasInputShapes -> cd(2);
    gPad -> SetLogy(true);
    TH2D *histGridEtaJpsi = new TH2D("histGridEtaJpsi", ";#it{#eta};Counts", 100, -5, 5, 100, 1e2, 5e4);
    histGridEtaJpsi -> Draw();
    graEtaMom -> Draw("EP SAME"); 
    funcEtaJpsi -> Draw("SAME");

    TF1 *funcAxe = new TF1("funcAxe", FuncAxe, 0., 15., 5);
    funcAxe -> SetParameter(0, 0.416885);
    funcAxe -> SetParameter(1, -0.0609379);
    funcAxe -> SetParameter(2, 0.0222385);
    funcAxe -> SetParameter(3, -0.00171191);
    funcAxe -> SetParameter(4, 4.08582e-05);

    TFile *fInMchTrkEff = new TFile("/Users/lucamicheletti/cernbox/JPSI/Jpsi_pp_cross_section_run3/mch_trk_eff/mch_trk_eff_LHC24_pass1_minBias.root", "READ");
    TH1D *histCorrMapMchTrkEffEta = (TH1D*) fInMchTrkEff -> Get("histCorrMap_Eta");
    TH1D *histCorrMapMchTrkEffPt = (TH1D*) fInMchTrkEff -> Get("histCorrMap_Pt");
    TH1D *histCorrMapMchTrkEffPhi = (TH1D*) fInMchTrkEff -> Get("histCorrMap_Phi");

    double matchEff = 0.986;

    const int nPtBins = 15;
    double pTBinEdges[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 10.0, 12.0, 16.0, 20.0};

    TH2D *histPtMomPtDau = new TH2D("histPtMomPtDau", ";#it{p}^{J/#psi}_{T} (GeV/#it{c});#it{p}^{#mu}_{T} (GeV/#it{c})", 100, 0, 15, 100, 0, 15);
    TH1D *histPtMomGen = new TH1D("histPtMomGen", Form(";#it{p}^{%s}_{T} (GeV/#it{c});Gen", momName.c_str()), nPtBins, pTBinEdges);
    TH1D *histPtMomRec = new TH1D("histPtMomRec", Form(";#it{p}^{%s}_{T} (GeV/#it{c});Rec", momName.c_str()), nPtBins, pTBinEdges);
    TH1D *histPtMomRecMchTrkEff = new TH1D("histPtMomRecMchTrkEff", Form(";#it{p}^{%s}_{T} (GeV/#it{c});Rec", momName.c_str()), nPtBins, pTBinEdges);
    TH1D *histPtMomRecMatchEff = new TH1D("histPtMomRecMatchEff", Form(";#it{p}^{%s}_{T} (GeV/#it{c});Rec", momName.c_str()), nPtBins, pTBinEdges);

    SetHist(histPtMomGen, kBlue, 20, 1, kBlue);
    SetHist(histPtMomRec, kRed, 20, 1, kRed);

    TLorentzVector *vecMom = new TLorentzVector();
    bool passInAccDaup = false;
    bool passInAccDaum = false;

    double weigthTrkEffDaup = 1.;
    double weigthTrkEffDaum = 1.;

    double weightMatchEffDaup = 1.;
    double weightMatchEffDaum = 1.;

    for (int iEvt = 0;iEvt < nEvts;iEvt++) {
        std::cout << "* * * Event : " << iEvt << " * * *" << std::endl; 
        passInAccDaup = false;
        passInAccDaum = false;

        weigthTrkEffDaup = 1.;
        weigthTrkEffDaum = 1.;

        weightMatchEffDaup = 1.;
        weightMatchEffDaum = 1.;

        double ptMom = funcPtJpsi -> GetRandom();
        double etaMom = funcEtaJpsi -> GetRandom();
        double phiMom = funcPhiJpsi -> GetRandom();
        vecMom -> SetPtEtaPhiM(ptMom, etaMom, phiMom, massMom);

        double pxMom  = vecMom -> Px();
        double pyMom  = vecMom -> Py(); 
        double pzMom  = vecMom -> Pz(); 
        double eMom   = vecMom -> Energy();
        double rapMom = vecMom -> Rapidity();

        if (ptMom > 15 || TMath::Abs(rapMom) < 2.5 || TMath::Abs(rapMom) > 4.0) continue;

        int index = pythia.event.append(pdgCodeMom, 1, 0, 0, pxMom, pyMom, pzMom, eMom, massMom);

        if (!pythia.moreDecays()) {
            std::cerr << "Error in the decay" << std::endl;
            return;
        }

        for (int i = 0; i < pythia.event.size(); ++i) {
            const Particle& dau = pythia.event[i];
            if (dau.isFinal()) {
                double etaDau = TMath::Abs(dau.eta());
                double ptDau = dau.pT();
                double phiDau = dau.phi();
                if (dau.id() == 13) {
                    histPtMomPtDau -> Fill(ptMom, ptDau);
                    double cutEta = histCorrMapMchTrkEffEta -> GetBinContent(histCorrMapMchTrkEffEta -> FindBin(etaDau));
                    double cutPt = histCorrMapMchTrkEffPt -> GetBinContent(histCorrMapMchTrkEffPt -> FindBin(ptDau));
                    double cutPhi = histCorrMapMchTrkEffPhi -> GetBinContent(histCorrMapMchTrkEffPhi -> FindBin(phiDau));
                    weigthTrkEffDaup = std::min({cutEta, cutPt, cutPhi});

                    if (ptDau > 0.5 && etaDau > 2.5 && etaDau < 4.0) {passInAccDaup = true;}
                }
                if (dau.id() == -13) {
                    histPtMomPtDau -> Fill(ptMom, ptDau);
                    double cutEta = histCorrMapMchTrkEffEta -> GetBinContent(histCorrMapMchTrkEffEta -> FindBin(etaDau));
                    double cutPt = histCorrMapMchTrkEffPt -> GetBinContent(histCorrMapMchTrkEffPt -> FindBin(ptDau));
                    double cutPhi = histCorrMapMchTrkEffPhi -> GetBinContent(histCorrMapMchTrkEffPhi -> FindBin(phiDau));
                    weigthTrkEffDaum = std::min({cutEta, cutPt, cutPhi});

                    if (ptDau > 0.5 && etaDau > 2.5 && etaDau < 4.0) {passInAccDaum = true;}
                }
            }
        }

        if (passInAccDaup && passInAccDaum) {
            histPtMomGen -> Fill(ptMom);
            histPtMomRec -> Fill(ptMom, funcAxe -> Eval(ptMom));
            histPtMomRecMchTrkEff -> Fill(ptMom, weigthTrkEffDaup * weigthTrkEffDaum);
            histPtMomRecMatchEff -> Fill(ptMom, (1 + (1 - matchEff)) * (1 + (1 - matchEff)));
        }

        pythia.event.reset();
    }


    TFile *fOut = new TFile(Form("ToyMcResults_%d.root", pdgCodeMom), "RECREATE");
    histPtMomGen -> Write();
    histPtMomRec -> Write();
    histPtMomRecMchTrkEff -> Write();
    histPtMomRecMatchEff -> Write();
    fOut -> Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void plotRatio() {
    LoadStyle();
    TFile *fInJpsi = new TFile("ToyMcResults_443.root", "READ");
    TH1D *histPtJpsiGen = (TH1D*) fInJpsi -> Get("histPtMomGen"); histPtJpsiGen -> Sumw2(true);
    TH1D *histPtJpsiRecMchTrkEff = (TH1D*) fInJpsi -> Get("histPtMomRecMchTrkEff"); histPtJpsiRecMchTrkEff -> Sumw2(true);

    TFile *fInPsi2s = new TFile("ToyMcResults_100443.root", "READ");
    TH1D *histPtPsi2sGen = (TH1D*) fInPsi2s -> Get("histPtMomGen"); histPtPsi2sGen -> Sumw2(true);
    TH1D *histPtPsi2sRecMchTrkEff = (TH1D*) fInPsi2s -> Get("histPtMomRecMchTrkEff"); histPtPsi2sRecMchTrkEff -> Sumw2(true);

    SetHist(histPtJpsiGen, kRed, 20, 1, kRed);
    SetHist(histPtJpsiRecMchTrkEff, kRed+1, 24, 1, kRed+1);
    SetHist(histPtPsi2sGen, kAzure, 20, 1, kAzure);
    SetHist(histPtPsi2sRecMchTrkEff, kBlue+1, 24, 1, kBlue+1);

    histPtJpsiGen -> Scale(1, "WIDTH");
    histPtJpsiRecMchTrkEff -> Scale(1, "WIDTH");
    histPtPsi2sGen -> Scale(1, "WIDTH");
    histPtPsi2sRecMchTrkEff -> Scale(1, "WIDTH");

    TCanvas *canvasPtGen = new TCanvas("canvasPtGen", "", 800, 600);
    gPad -> SetLogy(true);
    histPtJpsiGen -> GetXaxis() -> SetRangeUser(0, 16);
    histPtJpsiGen -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histPtJpsiGen -> GetYaxis() -> SetTitle("d#it{N} / d#it{p}_{T} (GeV/#it{c})^{-1}");
    histPtJpsiGen -> Draw("EP");
    histPtPsi2sGen -> Draw("EP SAME");
    histPtJpsiRecMchTrkEff -> Draw("EP SAME");
    histPtPsi2sRecMchTrkEff -> Draw("EP SAME");

    TLegend *legendPtGen = new TLegend(0.56, 0.68, 0.76, 0.93, " ", "brNDC");
    SetLegend(legendPtGen);
    legendPtGen -> AddEntry(histPtJpsiGen, "J/#psi w/o Trk. eff. corr", "EP");
    legendPtGen -> AddEntry(histPtPsi2sGen, "#psi(2S) w/o Trk. eff. corr", "EP");
    legendPtGen -> AddEntry(histPtJpsiRecMchTrkEff, "J/#psi w/ Trk. eff. corr", "EP");
    legendPtGen -> AddEntry(histPtPsi2sRecMchTrkEff, "#psi(2S) w/ Trk. eff. corr", "EP");
    legendPtGen -> Draw("SAME");

    TH1D *histRatioGen = (TH1D*) histPtPsi2sGen -> Clone("histRatioGen");
    histRatioGen -> Divide(histPtJpsiGen);

    TH1D *histRatioRecMchTrkEff = (TH1D*) histPtPsi2sRecMchTrkEff -> Clone("histRatioRecMchTrkEff");
    histRatioRecMchTrkEff -> Divide(histPtJpsiRecMchTrkEff);

    SetHist(histRatioGen, kRed, 20, 1, kRed);
    SetHist(histRatioRecMchTrkEff, kAzure, 24, 1, kAzure);

    TCanvas *canvasRatio = new TCanvas("canvasRatio", "", 800, 600);
    histRatioGen -> GetXaxis() -> SetRangeUser(0, 16);
    histRatioGen -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histRatioGen -> GetYaxis() -> SetTitle("#psi(2S) / J/#psi");
    histRatioGen -> Draw("EP");
    histRatioRecMchTrkEff -> Draw("EP SAME");

    TLegend *legendRatio = new TLegend(0.56, 0.23, 0.76, 0.38, " ", "brNDC");
    SetLegend(legendRatio);
    legendRatio -> AddEntry(histRatioGen, "#psi(2S) / J/#psi w/ Trk. eff. corr", "EP");
    legendRatio -> AddEntry(histRatioRecMchTrkEff, "#psi(2S) / J/#psi w/o Trk. eff. corr", "EP");
    legendRatio -> Draw("SAME");

    TH1D *histRatioOfRatio = (TH1D*) histRatioGen -> Clone("histRatioOfRatio");
    histRatioOfRatio -> Divide(histRatioRecMchTrkEff);
    histRatioOfRatio -> SetMarkerStyle(24);
    histRatioOfRatio -> SetMarkerColor(kBlack);
    histRatioOfRatio -> SetLineColor(kBlack);
    for (int iPt = 0;iPt < 15;iPt++) {
        histRatioOfRatio -> SetBinError(iPt+1, 0);
    }

    TCanvas *canvasRatioOfRatio = new TCanvas("canvasRatioOfRatio", "", 800, 600);
    histRatioOfRatio -> GetXaxis() -> SetRangeUser(0, 16);
    histRatioOfRatio -> GetYaxis() -> SetRangeUser(0.95, 1.05);
    histRatioOfRatio -> GetYaxis() -> SetTitle("(#psi(2S) / J/#psi) / (#psi(2S) / J/#psi with Trk. eff. corr)");
    histRatioOfRatio -> Draw("P");

    TLine *lineUnity = new TLine(0, 1, 16, 1);
    lineUnity -> SetLineWidth(2);
    lineUnity -> SetLineColor(kGray+1);
    lineUnity -> SetLineStyle(kDashed);
    lineUnity -> Draw();

    canvasPtGen -> SaveAs("figures/PtGenDistributions.pdf");
    canvasRatio -> SaveAs("figures/RatioDistributions.pdf");
    canvasRatioOfRatio -> SaveAs("figures/RatioOfRatioDistributions.pdf");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double FuncPt(double *x, double *par) {
    double p0 = par[0];
    double p1 = par[1];
    double p2 = par[2];
    double p3 = par[3];
    return (p0 * x[0]) / TMath::Power(1. + TMath::Power(x[0] / p1, p2), p3);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double FuncEta(double *x, double *par) {
    double p0 = par[0];
    double p1 = par[1];
    double p2 = par[2];
    return p0 * TMath::Exp(-(1. / 2.) * TMath::Power(((x[0] - p1) / p2), 2));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double FuncAxe(double *x, double *par) {
    double xx = x[0];
    return par[0] + par[1] * xx + par[2] * xx * xx + par[3] * xx * xx * xx + par[4] * xx * xx * xx * xx;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.04);
}