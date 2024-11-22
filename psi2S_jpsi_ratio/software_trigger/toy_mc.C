#include "TProfile.h"
#include <fstream>
#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"
#include "TF1.h"
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

void toy_mc(int nEvs = 100000000, double ptCut = 0.7) {
    TFile *fIn = new TFile("toy_mc_input_shapes.root", "READ");
    TH1D *histPtGen = (TH1D*) fIn -> Get("Pt");
    TH1D *histEtaGen = (TH1D*) fIn -> Get("Eta");
    TH1D *histPhiGen = (TH1D*) fIn -> Get("Phi");

    TH2D *histDimu = new TH2D("histDimu", "", 100, 0, 5, 20, 0, 10);
    TH2D *histDimuPtCutRun3 = new TH2D("histDimuPtCutRun3", "", 100, 0, 5, 20, 0, 10);
    TH2D *histDimuPtCutRun2 = new TH2D("histDimuPtCutRun2", "", 100, 0, 5, 20, 0, 10);

    TH1D *histMuPt = new TH1D("histMuPt", "", 100, 0, 10);
    histMuPt -> SetLineColor(kBlue+1);
    histMuPt -> SetMarkerColor(kBlue+1);

    TH1D *histMuPtCutRun3 = new TH1D("histMuPtCutRun3", "", 100, 0, 10);
    histMuPtCutRun3 -> SetLineColor(kRed+1);
    histMuPtCutRun3 -> SetMarkerColor(kRed+1);

    TH1D *histMuPtCutRun2 = new TH1D("histMuPtCutRun2", "", 100, 0, 10);
    histMuPtCutRun2 -> SetLineColor(kGreen+1);
    histMuPtCutRun2 -> SetMarkerColor(kGreen+1);
    

    TF1 *funcPol1 = new TF1("funcPol1", "pol1", 0, 2);
    funcPol1 -> SetParameter(0, 0);
    funcPol1 -> SetParameter(1, 1);

    double ptCutSmeared = ptCut;

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
        histMuPt -> Fill(ptMu1);
        histMuPt -> Fill(ptMu2);

        if (ptMu1 > ptCut && ptMu2 > ptCut) {
            histDimuPtCutRun3 -> Fill(dimuon.M(), dimuon.Pt());
            histMuPtCutRun3 -> Fill(ptMu1);
            histMuPtCutRun3 -> Fill(ptMu2);
        }

        ptCutSmeared = funcPol1 -> GetRandom(0, 1);
        if (ptMu1 > ptCutSmeared && ptMu2 > ptCutSmeared) {
            histDimuPtCutRun2 -> Fill(dimuon.M(), dimuon.Pt());
            histMuPtCutRun2 -> Fill(ptMu1);
            histMuPtCutRun2 -> Fill(ptMu2);
        }
    }

    TCanvas *canvasDimu = new TCanvas("canvasDimu", "", 800, 600);
    histDimu -> Draw("COLZ");

    TCanvas *canvasDimuPtCutRun3 = new TCanvas("canvasDimuPtCutRun3", "", 800, 600);
    histDimuPtCutRun3 -> Draw("COLZ");

    TCanvas *canvasDimuPtCutRun2 = new TCanvas("canvasDimuPtCutRun2", "", 800, 600);
    histDimuPtCutRun2 -> Draw("COLZ");

    TCanvas *canvasMuPt = new TCanvas("canvasMuPt", "", 800, 600);
    gPad -> SetLogy(1);
    histMuPt -> Draw("EP");
    histMuPtCutRun3 -> Draw("EP SAME");
    histMuPtCutRun2 -> Draw("EP SAME");

    TLegend *legendMuPt = new TLegend(0.60, 0.50, 0.80, 0.80, " ", "brNDC");
    SetLegend(legendMuPt);
    legendMuPt -> SetTextSize(0.03);
    legendMuPt -> AddEntry(histMuPt, "#it{p}_{T}^{#mu} > 0", "PL");
    legendMuPt -> AddEntry(histMuPtCutRun3, Form("#it{p}_{T}^{#mu} > %2.1f GeV/c", ptCut), "PL");
    legendMuPt -> AddEntry(histMuPtCutRun2, Form("#it{p}_{T}^{#mu} > %2.1f GeV/c smearing", ptCut), "PL");
    legendMuPt -> Draw("SAME");

    TCanvas *canvasMassProj = new TCanvas("canvasMassProj", "", 1800, 1800);
    canvasMassProj -> Divide(4, 3);

    TCanvas *canvasRatioMassProj = new TCanvas("canvasRatioMassProj", "", 1800, 1800);
    canvasRatioMassProj -> Divide(4, 3);

    for (int iPt = 0;iPt < 12;iPt++) {
        canvasMassProj -> cd(iPt+1);

        float minPt = histDimu -> GetYaxis() -> GetBinLowEdge(iPt+1);
        float maxPt = histDimu -> GetYaxis() -> GetBinLowEdge(iPt+1) + histDimu -> GetYaxis() -> GetBinWidth(iPt+1);

        TH1D *histMassProj = (TH1D*) histDimu -> ProjectionX(Form("%2.1f < pT < %2.1f GeV/c", minPt, maxPt), iPt+1, iPt+1);
        histMassProj -> SetLineColor(kBlue+1);
        histMassProj -> SetMarkerColor(kBlue+1);

        TH1D *histMassProjPtCutRun3 = (TH1D*) histDimuPtCutRun3 -> ProjectionX(Form("%2.1f < pT < %2.1f GeV/c [pT cut]", minPt, maxPt), iPt+1, iPt+1);
        histMassProjPtCutRun3 -> SetLineColor(kRed+1);
        histMassProjPtCutRun3 -> SetMarkerColor(kRed+1);

        TH1D *histMassProjPtCutRun2 = (TH1D*) histDimuPtCutRun2 -> ProjectionX(Form("%2.1f < pT < %2.1f GeV/c [pT cut smearing]", minPt, maxPt), iPt+1, iPt+1);
        histMassProjPtCutRun2 -> SetLineColor(kGreen+1);
        histMassProjPtCutRun2 -> SetMarkerColor(kGreen+1);

        TLegend *legendMassProj = new TLegend(0.20, 0.20, 0.40, 0.50, " ", "brNDC");
        SetLegend(legendMassProj);
        legendMassProj -> SetTextSize(0.05);
        legendMassProj -> AddEntry(histMassProj, "#it{p}_{T}^{#mu} > 0", "PL");
        legendMassProj -> AddEntry(histMassProjPtCutRun3, Form("#it{p}_{T}^{#mu} > %2.1f GeV/c", ptCut), "PL");
        legendMassProj -> AddEntry(histMassProjPtCutRun2, Form("#it{p}_{T}^{#mu} > %2.1f GeV/c smearing", ptCut), "PL");

        gPad -> SetLogy(1);
        histMassProj -> Draw("EP");
        histMassProjPtCutRun3 -> Draw("EP SAME");
        histMassProjPtCutRun2 -> Draw("EP SAME");
        legendMassProj -> Draw("SAME");

        TH1D *histRatioMassProj = (TH1D*) histMassProj -> Clone(Form("Ratio %2.1f < pT < %2.1f GeV/c [pT cut]", minPt, maxPt));
        histRatioMassProj -> Divide(histMassProjPtCutRun3);

        canvasRatioMassProj -> cd(iPt+1);
        histRatioMassProj -> Draw("EP");
    }
    canvasMassProj -> SaveAs("toy_mc_result.pdf");
    canvasRatioMassProj -> SaveAs("ratio_toy_mc_result.pdf");
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