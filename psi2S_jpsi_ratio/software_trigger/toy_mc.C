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
#include "TLorentzVector.h"

void SetLegend(TLegend *);

void toy_mc() {
    TFile *fIn = new TFile("toy_mc_input_shapes.root", "READ");
    TH1D *histPtGen = (TH1D*) fIn -> Get("Pt");
    TH1D *histEtaGen = (TH1D*) fIn -> Get("Eta");
    TH1D *histPhiGen = (TH1D*) fIn -> Get("Phi");

    TH2D *histDimu = new TH2D("histDimu", "", 100, 0, 5, 20, 0, 10);
    TH2D *histDimuPtCut = new TH2D("histDimuPtCut", "", 100, 0, 5, 20, 0, 10);

    const int nEvs = 100000000;
    const double ptCut = 0.7;

    for (int iEv = 0;iEv < nEvs;iEv++) {
        double ptMu1 = histPtGen -> GetRandom();
        double etaMu1 = histEtaGen -> GetRandom();
        double phiMu1 = histPhiGen -> GetRandom();
        ROOT::Math::PtEtaPhiMVector mu1(ptMu1, etaMu1, phiMu1, 0.105);

        double ptMu2 = histPtGen -> GetRandom();
        double etaMu2 = histEtaGen -> GetRandom();
        double phiMu2 = histPhiGen -> GetRandom();
        ROOT::Math::PtEtaPhiMVector mu2(ptMu2, etaMu2, phiMu2, 0.105);

        auto dimuon = mu1 + mu2;
        histDimu -> Fill(dimuon.M(), dimuon.Pt());

        if (ptMu1 > ptCut && ptMu2 > ptCut) {
            histDimuPtCut -> Fill(dimuon.M(), dimuon.Pt());
        }
    }

    TCanvas *canvasDimu = new TCanvas("canvasDimu", "", 800, 600);
    histDimu -> Draw("COLZ");

    TCanvas *canvasDimuPtCut = new TCanvas("canvasDimuPtCut", "", 800, 600);
    histDimuPtCut -> Draw("COLZ");

    TCanvas *canvasMassProj = new TCanvas("canvasMassProj", "", 1800, 1800);
    canvasMassProj -> Divide(4, 3);
    for (int iPt = 0;iPt < 12;iPt++) {
        canvasMassProj -> cd(iPt+1);

        float minPt = histDimu -> GetYaxis() -> GetBinLowEdge(iPt+1);
        float maxPt = histDimu -> GetYaxis() -> GetBinLowEdge(iPt+1) + histDimu -> GetYaxis() -> GetBinWidth(iPt+1);

        TH1D *histMassProj = (TH1D*) histDimu -> ProjectionX(Form("%2.1f < pT < %2.1f GeV/c", minPt, maxPt), iPt+1, iPt+1);
        histMassProj -> SetLineColor(kBlue);
        histMassProj -> SetMarkerColor(kBlue);

        TH1D *histMassProjPtCut = (TH1D*) histDimuPtCut -> ProjectionX(Form("%2.1f < pT < %2.1f GeV/c [pT cut]", minPt, maxPt), iPt+1, iPt+1);
        histMassProjPtCut -> SetLineColor(kRed);
        histMassProjPtCut -> SetMarkerColor(kRed);

        TLegend *legendMassProj = new TLegend(0.20, 0.20, 0.40, 0.50, " ", "brNDC");
        SetLegend(legendMassProj);
        legendMassProj -> SetTextSize(0.05);
        legendMassProj -> AddEntry(histMassProj, "#it{p}_{T}^{#mu} > 0", "PL");
        legendMassProj -> AddEntry(histMassProjPtCut, "#it{p}_{T}^{#mu} > 0.7 GeV/c", "PL");

        gPad -> SetLogy(1);
        histMassProj -> Draw("EP");
        histMassProjPtCut -> Draw("EP SAME");
        legendMassProj -> Draw("SAME");
    }
    canvasMassProj -> SaveAs("toy_mc_result.pdf");
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