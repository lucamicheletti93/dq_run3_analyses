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
#include <TDatabasePDG.h>

void SetLegend(TLegend *);

void toy_mc(int nEvs = 100000000, double ptCut = 0.7) {
    TDatabasePDG *db = TDatabasePDG::Instance();
    const int pdgMu = 13;
    double massMu = db -> GetParticle(pdgMu) -> Mass();

    const int nPtBins = 13;
    double minPtBins[] = {0.0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 15};
    double maxPtBins[] = {0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20};
    double ptBins[] = {0.0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20};

    TFile *fIn = new TFile("toy_mc_input_shapes.root", "READ");
    TH1D *histPtGen = (TH1D*) fIn -> Get("Pt");
    TH1D *histEtaGen = (TH1D*) fIn -> Get("Eta");
    TH1D *histPhiGen = (TH1D*) fIn -> Get("Phi");

    TH2D *histDimu = new TH2D("histDimu", "", 100, 1.9, 5, nPtBins, ptBins);
    histDimu -> GetXaxis() -> SetTitle("#it{M} GeV/c^{2}");
    histDimu -> GetYaxis() -> SetTitle("#it{p}_{T} GeV/c");

    TH2D *histDimuPtCutRun3 = new TH2D("histDimuPtCutRun3", "", 100, 1.9, 5, nPtBins, ptBins);
    histDimuPtCutRun3 -> GetXaxis() -> SetTitle("#it{M} GeV/c^{2}");
    histDimuPtCutRun3 -> GetYaxis() -> SetTitle("#it{p}_{T}  GeV/c");

    TH2D *histDimuPtCutRun2 = new TH2D("histDimuPtCutRun2", "", 100, 1.9, 5, nPtBins, ptBins);
    histDimuPtCutRun2 -> GetXaxis() -> SetTitle("#it{M} GeV/c^{2}");
    histDimuPtCutRun2 -> GetYaxis() -> SetTitle("#it{p}_{T}  GeV/c");

    TH1D *histMuPt = new TH1D("histMuPt", "", 100, 0, 10);
    histMuPt -> SetLineColor(kBlue+1);
    histMuPt -> SetMarkerColor(kBlue+1);

    TH1D *histMuPtCutRun3 = new TH1D("histMuPtCutRun3", "", 100, 0, 10);
    histMuPtCutRun3 -> SetLineColor(kRed+1);
    histMuPtCutRun3 -> SetMarkerColor(kRed+1);

    TH1D *histMuPtCutRun2 = new TH1D("histMuPtCutRun2", "", 100, 0, 10);
    histMuPtCutRun2 -> SetLineColor(kGreen+1);
    histMuPtCutRun2 -> SetMarkerColor(kGreen+1);
    

    //TF1 *funcPol1 = new TF1("funcPol1", "pol1", 0, 2);
    //funcPol1 -> SetParameter(0, 0);
    //funcPol1 -> SetParameter(1, 1);

    TF1 *funcPol1 = new TF1("funcPol1", "gaus", 0, 2);
    funcPol1 -> SetParameter(0, 1);
    funcPol1 -> SetParameter(1, ptCut);
    funcPol1 -> SetParameter(2, 0.15);

    TH1D *histPtThr = new TH1D("histPtThr", "", 1000, 0, 2);
    histPtThr -> SetLineColor(kBlue+1);
    histPtThr -> SetMarkerColor(kBlue+1);

    TH1D *histIntegralPtThr = new TH1D("histPtThr", "", 1000, 0, 2);
    histIntegralPtThr -> SetLineColor(kBlack);
    histIntegralPtThr -> SetMarkerColor(kBlack);

    double ptCutSmeared = ptCut;

    for (int iEv = 0;iEv < nEvs;iEv++) {
        double ptMu1 = histPtGen -> GetRandom();
        double etaMu1 = histEtaGen -> GetRandom();
        double phiMu1 = histPhiGen -> GetRandom();
        ROOT::Math::PtEtaPhiMVector mu1(ptMu1, etaMu1, phiMu1, massMu);

        double ptMu2 = histPtGen -> GetRandom();
        double etaMu2 = histEtaGen -> GetRandom();
        double phiMu2 = histPhiGen -> GetRandom();
        ROOT::Math::PtEtaPhiMVector mu2(ptMu2, etaMu2, phiMu2, massMu);

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
        histPtThr -> Fill(ptCutSmeared);
        if (ptMu1 > ptCutSmeared && ptMu2 > ptCutSmeared) {
            histDimuPtCutRun2 -> Fill(dimuon.M(), dimuon.Pt());
            histMuPtCutRun2 -> Fill(ptMu1);
            histMuPtCutRun2 -> Fill(ptMu2);
        }
    }

    int sumPreviousBins = 0;
    for (int iBin = 0;iBin < 1000;iBin++) {
        sumPreviousBins += histPtThr -> GetBinContent(iBin+1);
        histIntegralPtThr -> SetBinContent(iBin, sumPreviousBins);
    }

    TCanvas *canvasPtThr = new TCanvas("canvasPtThr", "", 800, 600);
    histPtThr -> Scale(1. / histPtThr -> GetMaximum());
    histIntegralPtThr -> Scale(1. / histIntegralPtThr -> GetBinContent(999));
    histPtThr -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/c)");
    histPtThr -> Draw("EP");
    histIntegralPtThr -> Draw("EP SAME");

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
    canvasMassProj -> Divide(5, 3);

    TCanvas *canvasMassProjLowPt = new TCanvas("canvasMassProjLowPt", "", 1800, 300);
    canvasMassProjLowPt -> Divide(5, 1);

    TCanvas *canvasMassProjvsRun2 = new TCanvas("canvasMassProjvsRun2", "", 1800, 1800);
    canvasMassProjvsRun2 -> Divide(5, 3);

    TCanvas *canvasRatioMassProj = new TCanvas("canvasRatioMassProj", "canvasRatioMassProj", 1800, 1800);
    canvasRatioMassProj -> Divide(5, 3);

    TCanvas *canvasDiffMassProj = new TCanvas("canvasDiffMassProj", "canvasDiffMassProj", 1800, 1800);
    canvasDiffMassProj -> Divide(5, 3);

    TFile *fOut = new TFile("toy_mc_output.root", "RECREATE");
    for (int iPt = 0;iPt < 12;iPt++) {
        canvasMassProj -> cd(iPt+1);
        gStyle -> SetOptStat(false);

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

        gPad -> SetLogy(1);
        histMassProj -> SetTitle(Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}", minPt, maxPt));
        histMassProj -> Draw("EP");
        histMassProjPtCutRun3 -> Draw("EP SAME");
        legendMassProj -> Draw("SAME");

        TH1D *histRatioMassProj = (TH1D*) histMassProj -> Clone(Form("Ratio %2.1f < pT < %2.1f GeV/c [pT cut]", minPt, maxPt));
        histRatioMassProj -> Divide(histMassProjPtCutRun3);

        canvasRatioMassProj -> cd(iPt+1);
        histRatioMassProj -> Draw("EP");

        TH1D *histDiffMassProj = (TH1D*) histMassProj -> Clone(Form("Diff. %2.1f < pT < %2.1f GeV/c [pT cut]", minPt, maxPt));
        histDiffMassProj -> Add(histMassProjPtCutRun3, -1);

        canvasDiffMassProj -> cd(iPt+1);
        histDiffMassProj -> Draw("EP");

        if (iPt < 5) {
            canvasMassProjLowPt -> cd(iPt+1);
            gPad -> SetLogy(1);
            histMassProj -> SetTitle(Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}", minPt, maxPt));
            histMassProj -> Draw("EP");
            histMassProjPtCutRun3 -> Draw("EP SAME");
            legendMassProj -> Draw("SAME");
        }


        canvasMassProjvsRun2 -> cd(iPt+1);
        gStyle -> SetOptStat(false);

        TLegend *legendMassProjVsRun2 = new TLegend(0.20, 0.20, 0.40, 0.50, " ", "brNDC");
        SetLegend(legendMassProjVsRun2);
        legendMassProjVsRun2 -> SetTextSize(0.05);
        legendMassProjVsRun2 -> AddEntry(histMassProj, "#it{p}_{T}^{#mu} > 0", "PL");
        legendMassProjVsRun2 -> AddEntry(histMassProjPtCutRun3, Form("#it{p}_{T}^{#mu} > %2.1f GeV/c", ptCut), "PL");
        legendMassProjVsRun2 -> AddEntry(histMassProjPtCutRun2, Form("#it{p}_{T}^{#mu} > %2.1f GeV/c smearing", ptCut), "PL");

        gPad -> SetLogy(1);
        histMassProj -> SetTitle(Form("%2.1f < #it{p}_{T} < %2.1f GeV/#it{c}", minPt, maxPt));
        histMassProj -> Draw("EP");
        histMassProjPtCutRun3 -> Draw("EP SAME");
        histMassProjPtCutRun2 -> Draw("EP SAME");
        legendMassProjVsRun2 -> Draw("SAME");

        histMassProj -> Write(Form("histMass_Pt_%2.1f_%2.1f", minPt, maxPt));
        histMassProjPtCutRun3 -> Write(Form("histMassPtCut_Pt_%2.1f_%2.1f", minPt, maxPt));
    }
    canvasPtThr -> SaveAs("figures/run2_pt_threshold.pdf");
    canvasMassProj -> SaveAs("figures/toy_mc_result.pdf");
    canvasMassProjLowPt -> SaveAs("figures/toy_mc_result_low_pt.pdf");
    canvasMassProjvsRun2 -> SaveAs("figures/toy_mc_result_vs_Run2.pdf");
    canvasRatioMassProj -> SaveAs("figures/ratio_toy_mc_result.pdf");
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