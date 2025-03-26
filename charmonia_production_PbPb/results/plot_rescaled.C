#include <TFile.h>
#include <TCanvas.h>
#include <TObject.h>
#include <TLegend.h>
#include <TString.h>
#include <RooCurve.h>
#include <TPaveText.h>
#include <iostream>
#include <RooHist.h> 
#include <TPaveText.h>
#include <TPave.h>
#include "TAxis.h"

void LoadStyle();
void SetLegend(TLegend *);

#include <RooHist.h>
#include <TH1D.h>

inline TH1D* RooHistToTH1D(const RooHist* rooHist, const char* histName = "histo") {
    int nPoints = rooHist->GetN();
    double x, y, errLow, errHigh;

    rooHist -> GetPoint(0, x, y);
    double xMin = x;
    rooHist -> GetPoint(nPoints - 1, x, y);
    double xMax = x;
    TH1D* hist = new TH1D(histName, "", nPoints, xMin, xMax);

    for (int i = 0; i < nPoints; i++) {
        rooHist -> GetPoint(i, x, y);
        int bin = hist -> FindBin(x);
        hist -> SetBinContent(bin, y);
        errLow = rooHist -> GetErrorYlow(i);
        errHigh = rooHist -> GetErrorYhigh(i);
        hist -> SetBinError(bin, (errLow + errHigh) / 2.0);
    }

    return hist;
}


void plot_rescaled() {
    LoadStyle();
    // Apri il file ROOT
    TFile *file = TFile::Open("multi_trial_NA60_VWG_1.05_width_2_8_MC_tails_MC_tails.root");
    if (!file || file -> IsZombie()) {
        std::cerr << "Error: impossible to open ROOT file!" << std::endl;
        return;
    }

    // Prendi il Canvas
    TCanvas *canvasFit = (TCanvas*)file -> Get("fit_plot_NA60_NA60_VWG__2.0_5.0_1.05_2_8_data_tails_Mass_BkgSub_Int_ME");
    if (!canvasFit) {
        std::cerr << "Error: TCanvas not found!" << std::endl;
        return;
    }

    RooCurve *psi2S = nullptr;
    RooCurve *jpsi = nullptr;
    RooCurve *bkg = nullptr;
    RooHist *sum = nullptr;
    RooHist *rooHistData = nullptr;

    TIter next(canvasFit->GetListOfPrimitives());
    TObject *obj;
    while ((obj = next())) {
        std::cout << obj->GetName() << std::endl;

        TString rooHistName = obj->GetName();
        if (rooHistName.Contains("h_data")) {
            rooHistData = (RooHist*)obj;
        }
        
        if (obj->InheritsFrom(RooCurve::Class())) {
            RooCurve *curve = (RooCurve*)obj;
            TString name = curve->GetName();
            if (name.Contains("Psi2sPdf")) {
                psi2S = curve;
            } else if (name.Contains("JpsiPdf")) {
                jpsi = curve;
            } else if (name.Contains("BkgPdf")) {
                bkg = curve;
            } else { 
                sum = (RooHist*)obj;
            }
        }
    }

    TH1D *histData = (TH1D*) RooHistToTH1D(rooHistData, "histData");
    histData -> SetMarkerStyle(20);
    histData -> SetMarkerSize(0.8);
    histData -> SetMarkerColor(kBlack);
    histData -> SetLineColor(kBlack);

    if (!psi2S) {
        std::cerr << "Error: Ψ(2S) RooCurve not found!" << std::endl;
        return;
    }
    if (!jpsi) {
        std::cerr << "Error: J/Ψ RooCurve not found!" << std::endl;
        return;
    }
    if (!bkg) {
        std::cerr << "Error: background RooCurve not found!" << std::endl;
        return;
    }
    RooCurve *totalFit_notRescaled = new RooCurve(*bkg);

    for (int i = 0; i < bkg->GetN(); i++) {
        double x_bkg, y_bkg;
        bkg->GetPoint(i, x_bkg, y_bkg);

        double y_jpsi = 0;
        double y_psi2S = 0;

        for (int j = 0; j < jpsi->GetN(); j++) {
            double x_jpsi, y_jpsi_point;
            jpsi->GetPoint(j, x_jpsi, y_jpsi_point);
            if (x_bkg == x_jpsi) {
                y_jpsi = y_jpsi_point;
                break;
            }
        }
        for (int k = 0; k < psi2S->GetN(); k++) {
            double x_psi2S, y_psi2S_point;
            psi2S->GetPoint(k, x_psi2S, y_psi2S_point);
            if (x_bkg == x_psi2S) {
                y_psi2S = y_psi2S_point;
                break;
            }
        }

        double y_total = y_bkg + y_jpsi + y_psi2S;
        totalFit_notRescaled->SetPoint(i, x_bkg, y_total);
    }

    RooCurve *psi2S_scaled = new RooCurve(*psi2S);
    psi2S_scaled->SetName("Psi2S_scaled");
    psi2S_scaled->SetTitle("Riscalata #Psi(2S)");

    for (int i = 0; i < psi2S_scaled->GetN(); i++) {
        double x, y;
        psi2S_scaled->GetPoint(i, x, y);
        psi2S_scaled->SetPoint(i, x, y * 0.0200823 / 0.00922245);
    }

    RooCurve *totalFit_modified = new RooCurve(*bkg); 

    for (int i = 0; i < bkg->GetN(); i++) {
        double x_bkg, y_bkg;
        bkg->GetPoint(i, x_bkg, y_bkg);

        double y_jpsi = 0;
        double y_psi2S = 0;

        for (int j = 0; j < jpsi->GetN(); j++) {
            double x_jpsi, y_jpsi_point;
            jpsi->GetPoint(j, x_jpsi, y_jpsi_point);
            if (x_bkg == x_jpsi) {
                y_jpsi = y_jpsi_point;
                break;
            }
        }

        for (int k = 0; k < psi2S_scaled->GetN(); k++) {
            double x_psi2S, y_psi2S_point;
            psi2S_scaled->GetPoint(k, x_psi2S, y_psi2S_point);
            if (x_bkg == x_psi2S) {
                y_psi2S = y_psi2S_point;
                break;
            }
        }

        double y_total = y_bkg + y_jpsi + y_psi2S;
        totalFit_modified->SetPoint(i, x_bkg, y_total);
    }

    TCanvas *canvasFigure = new TCanvas("canvasFigure", "", 800, 800);
    canvasFigure->SetLeftMargin(0.12);
    canvasFigure->SetRightMargin(0.05);
    canvasFigure->SetTopMargin(0.1);
    canvasFigure->SetBottomMargin(0.15);

    totalFit_modified->GetXaxis()->SetTitleSize(0.045);
    totalFit_modified->GetYaxis()->SetTitleSize(0.045);
    totalFit_modified->GetXaxis()->SetLabelSize(0.035);
    totalFit_modified->GetYaxis()->SetLabelSize(0.035);

    canvasFigure->cd();
    gPad -> SetLogy(true);
    TH2D *histGrid = new TH2D("histGrid", "", 100, 2, 5, 100, 1e2, 1e7);
    histGrid -> GetXaxis() -> SetTitle("#it{M}_{#mu^{+}#mu^{-}} (GeV/#it{c}^{2})");
    histGrid -> GetYaxis() -> SetTitle("Counts");
    histGrid -> Draw();
    histData -> Draw("EP SAME");

    totalFit_modified->SetLineColor(kBlue);
    totalFit_modified->SetLineStyle(kDashed);  
    totalFit_modified->SetLineWidth(2);
    totalFit_modified->Draw("L SAME");

    totalFit_notRescaled->SetLineColor(kRed);  
    totalFit_notRescaled->SetLineStyle(1);     
    totalFit_notRescaled->SetLineWidth(2);

    bkg->SetLineColor(kGray+1);
    bkg->SetLineStyle(kDotted);
    bkg->SetLineWidth(2);
    bkg->Draw("L SAME");

    totalFit_notRescaled->Draw("L SAME");

    jpsi->SetLineColor(kAzure+7);
    jpsi->SetLineStyle(1);
    jpsi->SetLineWidth(2);
    jpsi->Draw("L SAME");

    psi2S->SetLineColor(kGreen-3);
    psi2S->SetLineStyle(1);
    psi2S->SetLineWidth(2);
    psi2S->Draw("L SAME");

    TLegend *legend = new TLegend(0.70, 0.43, 0.85, 0.68);
    SetLegend(legend);
    legend -> AddEntry(histData, "Data", "lep");
    legend -> AddEntry(totalFit_notRescaled, "Fit", "l");
    legend -> AddEntry(jpsi, "J/#psi", "l");
    legend -> AddEntry(psi2S, "#psi(2S)", "l");
    legend -> AddEntry(bkg, "background", "l");
    legend -> AddEntry(totalFit_modified, "pp overlaid", "l");
    legend -> Draw("SAME");

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.040);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);
    latexTitle -> DrawLatex(0.15, 0.84, "ALICE Preliminary, Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.15, 0.78, "J/#psi#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4, 0#minus90\%");
    latexTitle -> DrawLatex(0.15, 0.72, "#it{p}_{T}^{#mu#mu} < 20 GeV/#it{c}, #it{p}_{T}^{#mu} > 0.7 GeV/#it{c}");

    gPad->SetTicks(1, 1);

    canvasFigure->Update();

    canvasFigure -> SaveAs("inv_mass_pp_overlaid.pdf");
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
    legend -> SetTextSize(0.038);
}



