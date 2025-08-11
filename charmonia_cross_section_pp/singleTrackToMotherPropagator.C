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
    pythia.readString("443:mayDecay = on");
    pythia.readString("443:oneChannel = 1 1.0 0 13 -13"); // BR = 100% in mu+ mu-
    pythia.init();

    double massMom = pythia.particleData.m0(pdgCodeMom);
    double massDau = pythia.particleData.m0(pdgCodeDau);

    string momName = "J/#psi";

    TF1 *funcPt = new TF1("funcPt", FuncPt, 0., 15., 3);
    funcPt -> SetParameter(0, 2.85);
    funcPt -> SetParameter(1, 2.81);
    funcPt -> SetParameter(2, 2.43);

    TF1 *funcEta = new TF1("funcEta", FuncEta, -5.0, -2.0, 2);
    funcEta -> SetParameter(0, -2.43);
    funcEta -> SetParameter(1, 1.08);

    TF1 *funcPhi = new TF1("funcPhi", "[0]", 0., 2. * TMath::Pi());
    funcPhi -> SetParameter(0, 1);

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

    const int nPtBins = 10;
    double pTBinEdges[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0, 20.0};

    TH2D *histPtMomPtDau = new TH2D("histPtMomPtDau", ";#it{p}^{J/#psi}_{T} (GeV/#it{c});#it{p}^{#mu}_{T} (GeV/#it{c})", 100, 0, 15, 100, 0, 15);
    TH1D *histPtMomGen = new TH1D("histPtMomGen", Form(";#it{p}^{%s}_{T} (GeV/#it{c});Gen", momName.c_str()), nPtBins, pTBinEdges);
    TH1D *histPtMomRec = new TH1D("histPtMomRec", Form(";#it{p}^{%s}_{T} (GeV/#it{c});Rec", momName.c_str()), nPtBins, pTBinEdges);

    TH1D *histPtMomRecWithWeights = new TH1D("histPtMomRecWithWeights", Form(";#it{p}^{%s}_{T} (GeV/#it{c});Rec", momName.c_str()), nPtBins, pTBinEdges);

    TH1D *histPtMomRecMchTrkEff = new TH1D("histPtMomRecMchTrkEff", Form(";#it{p}^{%s}_{T} (GeV/#it{c});Rec", momName.c_str()), nPtBins, pTBinEdges);
    TH1D *histPtMomRecMchTrkEffWithWeights = new TH1D("histPtMomRecMchTrkEffWithWeights", Form(";#it{p}^{%s}_{T} (GeV/#it{c});Rec", momName.c_str()), nPtBins, pTBinEdges);
    TH1D *histPtMomRecMatchEff = new TH1D("histPtMomRecMatchEff", Form(";#it{p}^{%s}_{T} (GeV/#it{c});Rec", momName.c_str()), nPtBins, pTBinEdges);
    TH1D *histPtMomRecMatchEffWithWeights = new TH1D("histPtMomRecMatchEffWithWeights", Form(";#it{p}^{%s}_{T} (GeV/#it{c});Rec", momName.c_str()), nPtBins, pTBinEdges);

    TH1D *histPtMomAxe = new TH1D("histPtMomAxe", Form(";#it{p}^{%s}_{T} (GeV/#it{c});A#times#varepsilon", momName.c_str()), nPtBins, pTBinEdges);
    TH1D *histPtMomAxeWithWeights = new TH1D("histPtMomAxeWithWeights", Form(";#it{p}^{%s}_{T} (GeV/#it{c});A#times#varepsilon", momName.c_str()), nPtBins, pTBinEdges);

    TH1D *histPtMomAxeMchTrkEff = new TH1D("histPtMomAxeMchTrkEff", Form(";#it{p}^{%s}_{T} (GeV/#it{c});A#times#varepsilon", momName.c_str()), nPtBins, pTBinEdges);
    TH1D *histPtMomAxeMchTrkEffWithWeights = new TH1D("histPtMomAxeMchTrkEffWithWeights", Form(";#it{p}^{%s}_{T} (GeV/#it{c});A#times#varepsilon", momName.c_str()), nPtBins, pTBinEdges);

    TH1D *histPtMomAxeMatchEff = new TH1D("histPtMomAxeMatchEff", Form(";#it{p}^{%s}_{T} (GeV/#it{c});A#times#varepsilon", momName.c_str()), nPtBins, pTBinEdges);
    TH1D *histPtMomAxeMatchEffWithWeights = new TH1D("histPtMomAxeMatchEffWithWeights", Form(";#it{p}^{%s}_{T} (GeV/#it{c});A#times#varepsilon", momName.c_str()), nPtBins, pTBinEdges);

    SetHist(histPtMomGen, kBlue, 20, 1, kBlue);
    SetHist(histPtMomRec, kRed, 20, 1, kRed);
    SetHist(histPtMomAxe, kBlack, 20, 1, kBlack);
    SetHist(histPtMomAxeWithWeights, kGray+2, 20, 1, kGray+2);
    SetHist(histPtMomAxeMchTrkEff, kRed+1, 24, 1, kRed+1);
    SetHist(histPtMomAxeMchTrkEffWithWeights, kOrange+7, 24, 1, kOrange+7);
    SetHist(histPtMomAxeMatchEff, kAzure+4, 24, 1, kAzure+4);
    SetHist(histPtMomAxeMatchEffWithWeights, kBlue, 24, 1, kBlue);

    TLorentzVector *vecMom = new TLorentzVector();
    bool passInAccDaup = false;
    bool passInAccDaum = false;

    bool passMchTrkEffDaup = false;
    bool passMchTrkEffDaum = false;

    bool passMatchEffDaup = false;
    bool passMatchEffDaum = false;

    double weigthTrkEffDaup = 1.;
    double weigthTrkEffDaum = 1.;

    double weightMatchEffDaup = 1.;
    double weightMatchEffDaum = 1.;

    for (int iEvt = 0;iEvt < nEvts;iEvt++) {
        std::cout << "* * * Event : " << iEvt << " * * *" << std::endl; 
        passInAccDaup = false;
        passInAccDaum = false;

        passMchTrkEffDaup = false;
        passMchTrkEffDaum = false;

        passMatchEffDaup = false;
        passMatchEffDaum = false;

        weigthTrkEffDaup = 1.;
        weigthTrkEffDaum = 1.;

        weightMatchEffDaup = 1.;
        weightMatchEffDaum = 1.;

        double ptMom = funcPt -> GetRandom();
        double etaMom = funcEta -> GetRandom();
        double phiMom = funcPhi -> GetRandom();
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
                    double rndmDaup = gRandom -> Rndm();
                    //std::cout << etaDau << " " << ptDau << " " << phiDau << std::endl;
                    double cutEta = histCorrMapMchTrkEffEta -> GetBinContent(histCorrMapMchTrkEffEta -> FindBin(etaDau));
                    double cutPt = histCorrMapMchTrkEffPt -> GetBinContent(histCorrMapMchTrkEffPt -> FindBin(ptDau));
                    double cutPhi = histCorrMapMchTrkEffPhi -> GetBinContent(histCorrMapMchTrkEffPhi -> FindBin(phiDau));
                    weigthTrkEffDaup = std::min({cutEta, cutPt, cutPhi});

                    if (ptDau > 0.5 && etaDau > 2.5 && etaDau < 4.0) {
                        passInAccDaup = true;
                        //std::cout << rndmDaup << " - " << cutEta << " -> eta = " << etaDau << std::endl;
                    }

                    if (rndmDaup < cutEta && rndmDaup < cutPt && rndmDaup < cutPhi) {
                        passMchTrkEffDaup = true;
                    }

                    if (rndmDaup < matchEff) {
                        passMatchEffDaup = true;
                    }
                }
                if (dau.id() == -13) {
                    histPtMomPtDau -> Fill(ptMom, ptDau);
                    double rndmDaum = gRandom -> Rndm();
                    //std::cout << etaDau << " " << ptDau << " " << phiDau << std::endl;
                    double cutEta = histCorrMapMchTrkEffEta -> GetBinContent(histCorrMapMchTrkEffEta -> FindBin(etaDau));
                    double cutPt = histCorrMapMchTrkEffPt -> GetBinContent(histCorrMapMchTrkEffPt -> FindBin(ptDau));
                    double cutPhi = histCorrMapMchTrkEffPhi -> GetBinContent(histCorrMapMchTrkEffPhi -> FindBin(phiDau));
                    weigthTrkEffDaum = std::min({cutEta, cutPt, cutPhi});

                    if (ptDau > 0.5 && etaDau > 2.5 && etaDau < 4.0) {
                        passInAccDaum = true;
                        //std::cout << rndmDaum << " - " << cutEta << " -> eta = " << etaDau << std::endl;
                    }
                    
                    if (rndmDaum < cutEta && rndmDaum < cutPt && rndmDaum < cutPhi) {
                        passMchTrkEffDaum = true;
                    }

                    if (rndmDaum < matchEff) {
                        passMatchEffDaum = true;
                    }
                }
            }
        }

        if (passInAccDaup && passInAccDaum) {
            histPtMomGen -> Fill(ptMom);
            histPtMomRecWithWeights -> Fill(ptMom, funcAxe -> Eval(ptMom));
            histPtMomRecMchTrkEffWithWeights -> Fill(ptMom, weigthTrkEffDaup * weigthTrkEffDaum * funcAxe -> Eval(ptMom));
            histPtMomRecMatchEffWithWeights -> Fill(ptMom, (1 + (1 - matchEff)) * (1 + (1 - matchEff)) * funcAxe -> Eval(ptMom));
            double axeCut = funcAxe -> Eval(ptMom);
            if (gRandom -> Rndm() < axeCut) {
                histPtMomRec -> Fill(ptMom);
                if (passMchTrkEffDaup && passMchTrkEffDaum) {
                    histPtMomRecMchTrkEff -> Fill(ptMom);
                }
                if (passMatchEffDaup && passMatchEffDaum) {
                    histPtMomRecMatchEff -> Fill(ptMom);
                }
            }
        }

        pythia.event.reset();
    }

    histPtMomAxe -> Divide(histPtMomRec, histPtMomGen, 1, 1, "B");
    histPtMomAxeWithWeights -> Divide(histPtMomRecWithWeights, histPtMomGen, 1, 1, "B");
    histPtMomAxeMchTrkEff -> Divide(histPtMomRecMchTrkEff, histPtMomGen, 1, 1, "B");
    histPtMomAxeMchTrkEffWithWeights -> Divide(histPtMomRecMchTrkEffWithWeights, histPtMomGen, 1, 1, "B");
    histPtMomAxeMatchEff -> Divide(histPtMomRecMatchEff, histPtMomGen, 1, 1, "B");
    histPtMomAxeMatchEffWithWeights -> Divide(histPtMomRecMatchEffWithWeights, histPtMomGen, 1, 1, "B");

    TH1D *histPtMomRatioAxeMchTrkEff = (TH1D*) histPtMomAxe -> Clone("histPtMomRatioAxeMchTrkEff");
    histPtMomRatioAxeMchTrkEff -> Divide(histPtMomAxeMchTrkEff);
    SetHist(histPtMomRatioAxeMchTrkEff, kRed+1, 24, 1, kRed+1);

    TH1D *histPtMomRatioAxeMchTrkEffWithWeights = (TH1D*) histPtMomAxe -> Clone("histPtMomRatioAxeMchTrkEffWithWeights");
    histPtMomRatioAxeMchTrkEffWithWeights -> Divide(histPtMomAxeMchTrkEffWithWeights);
    SetHist(histPtMomRatioAxeMchTrkEffWithWeights, kOrange+7, 24, 1, kOrange+7);

    TH1D *histPtMomRatioAxeMatchEff = (TH1D*) histPtMomAxe -> Clone("histPtMomRatioAxeMatchEff");
    histPtMomRatioAxeMatchEff -> Divide(histPtMomAxeMatchEff);
    SetHist(histPtMomRatioAxeMatchEff, kAzure+4, 24, 1, kAzure+4);

    TH1D *histPtMomRatioAxeMatchEffWithWeights = (TH1D*) histPtMomAxe -> Clone("histPtMomRatioAxeMatchEffWithWeights");
    histPtMomRatioAxeMatchEffWithWeights -> Divide(histPtMomAxeMatchEffWithWeights);
    SetHist(histPtMomRatioAxeMatchEffWithWeights, kBlue, 24, 1, kBlue);

    TF1 *funcPol0MchTrkEff = new TF1("funcPol0MchTrkEff", "[0]", 0, 20);
    funcPol0MchTrkEff -> SetParameter(0, 1);
    histPtMomRatioAxeMchTrkEff -> Fit(funcPol0MchTrkEff, "R");
    funcPol0MchTrkEff -> SetLineColor(kRed+1);

    TF1 *funcPol0MatchEff = new TF1("funcPol0MatchEff", "[0]", 0, 20);
    funcPol0MatchEff -> SetParameter(0, 1);
    histPtMomRatioAxeMatchEff -> Fit(funcPol0MatchEff, "R");
    funcPol0MatchEff -> SetLineColor(kAzure+4);

    TCanvas *canvasPt = new TCanvas("canvasPt", "", 800, 600);
    histPtMomGen -> Draw("EP");
    histPtMomRec -> Draw("EP SAME");

    double padAspectRatio = 0.6;

    TLatex latexTitle;
    latexTitle.SetNDC();
    latexTitle.SetTextSize(0.06);
    latexTitle.SetTextFont(42);

    TCanvas *canvasAxe = new TCanvas("canvasTagAndProbes", "", 800, 1000);
    gStyle -> SetOptFit(false);
    gStyle -> SetOptStat(false);

    TPad *pad1 = new TPad("pad1", "pad1", 0.05, 1-padAspectRatio, 0.99, 0.99);
    pad1 -> Draw();
    pad1 -> cd();
    histPtMomAxe -> GetYaxis() -> SetRangeUser(0, 1.2);
    histPtMomAxe -> Draw("EP");
    histPtMomAxeWithWeights -> Draw("EP SAME");
    histPtMomAxeMchTrkEff -> Draw("EP SAME");
    histPtMomAxeMchTrkEffWithWeights -> Draw("EP SAME");
    histPtMomAxeMatchEff -> Draw("EP SAME");
    histPtMomAxeMatchEffWithWeights -> Draw("EP SAME");

    TLegend *legendAxe = new TLegend(0.16, 0.70, 0.36, 0.95, " ", "brNDC");
    SetLegend(legendAxe);
    legendAxe -> AddEntry(histPtMomAxe, "Default", "EP");
    legendAxe -> AddEntry(histPtMomAxeMchTrkEff, "Default + MCH tracking", "EP");
    legendAxe -> AddEntry(histPtMomAxeMatchEff, "Default + Matching", "EP");
    legendAxe -> Draw("SAME");

    canvasAxe -> cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.05, 0.05, 0.99, 1-padAspectRatio);
    pad2 -> SetBottomMargin(0.25);
    pad2 -> Draw();
    pad2 -> cd();
    histPtMomRatioAxeMchTrkEff -> GetXaxis() -> SetTitleSize(0.05 * (1./padAspectRatio));
    histPtMomRatioAxeMchTrkEff -> GetYaxis() -> SetTitleSize(0.045 * (1./padAspectRatio));
    histPtMomRatioAxeMchTrkEff -> GetXaxis() -> SetLabelSize(0.045 * (1./padAspectRatio));
    histPtMomRatioAxeMchTrkEff -> GetYaxis() -> SetLabelSize(0.045 * (1./padAspectRatio));
    histPtMomRatioAxeMchTrkEff -> GetYaxis() -> SetRangeUser(0.8, 1.2);
    histPtMomRatioAxeMchTrkEff -> Draw("EP");
    histPtMomRatioAxeMatchEff -> Draw("EP SAME");
    histPtMomRatioAxeMchTrkEffWithWeights -> Draw("EP SAME");
    histPtMomRatioAxeMatchEffWithWeights -> Draw("EP SAME");

    TLine *lineUnity = new TLine(0., 1., 20., 1.);
    lineUnity -> SetLineStyle(kDashed);
    lineUnity -> SetLineColor(kGray+1);
    lineUnity -> SetLineWidth(2);
    lineUnity -> Draw("SAME");

    funcPol0MchTrkEff -> Draw("L SAME");
    funcPol0MatchEff -> Draw("L SAME");

    latexTitle.DrawLatex(0.70, 0.85, Form("#color[633]{p_{0} = %4.3f #pm %4.3f}", funcPol0MchTrkEff -> GetParameter(0), funcPol0MchTrkEff -> GetParError(0)));
    latexTitle.DrawLatex(0.70, 0.78, Form("#color[864]{p_{0} = %4.3f #pm %4.3f}", funcPol0MatchEff -> GetParameter(0), funcPol0MatchEff -> GetParError(0)));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double FuncPt(double *x, double *par) {
  double arg = 1 + TMath::Power(x[0]/par[0], par[1]);
  return x[0] / TMath::Power(arg, par[2]);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double FuncEta(double *x, double *par) {
  double arg = 0.5 * ((x[0] - par[0]) / par[1]) * ((x[0] - par[0]) / par[1]);
  return TMath::Exp(arg);
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