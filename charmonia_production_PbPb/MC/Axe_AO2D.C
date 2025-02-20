#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMinuit.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMath.h"
#include <fstream>
#include "THnSparse.h"
#include <TStyle.h>
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TTree.h"
#include "TBranch.h"
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <string>
#include <iterator>
#include <numeric>
#include <map>
#include "THStack.h"
#include <list>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>
#include <string_view>
#include "Fit/FitResult.h"
#include "TEllipse.h"
#include <tuple>
#include <TLegend.h>
#include "TProfile.h"
#include <TArrayD.h>

using namespace std;
TH1D* ProjectTHnSparse(THnSparseD *, double , double , double , double , double , double );
THnSparseD* CutTHnSparse(THnSparseD *, double , double , double , double , double , double );
TH2D* CutTH2Y(TH2D *, string , double , double);
TH2D* CutTH2X(TH2D *, string , double , double);
void SetLegend(TLegend * );
TCanvas *produceCanvas(bool logX, bool logY) {
	TCanvas *c = new TCanvas();
	c->SetFrameFillColor(0); //Transparent background
    c->SetFrameFillStyle(0);
    c->SetFrameBorderMode(0);
	c->SetGridx(); //Grid on x
	c->SetGridy(); //Grid on y 
	if (logX) c->SetLogx(); //Log on x, if needed
	if (logY) c->SetLogy(); //Log on y, if needed
	return c;
}
TH1F *produceTH1(string name, string title, string xAxisTitle, string yAxisTitle, int nBins, double xMin, double xMax) {
	TH1F *h = new TH1F(name.c_str(), title.c_str(), nBins, xMin, xMax);
	h->GetXaxis()->SetTitle(xAxisTitle.c_str());
	h->GetYaxis()->SetTitle(yAxisTitle.c_str());
	h->GetXaxis()->SetTitleFont(62);
	h->GetXaxis()->SetLabelFont(62);
	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->SetLabelFont(62);
	h->GetXaxis()->CenterTitle(true);
	h->GetYaxis()->CenterTitle(true);
	return h;
}

void Axe_AO2D(){

	double ptBinsRun3[] = {0.0, 2.0, 4.0, 6.0, 12.0};
    double centrBins[] = {0.0, 20.0, 40.0, 60.0, 90.0};
    const int nBinsPt = 4;
    const int nBinsRap = 6;
    double rapBins[] = {2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00};


    //------------------------------------------------- Generated MC events -------------------------------------------------//
    string pathIn = "/Users/saragaretti/dq_fitter_tails/inputShapes/plots_prova";
    string pathOut = "/Users/saragaretti/dq_fitter_tails/inputShapes/plotsAxe_prova";

	TFile *fIn = new TFile(Form("%s/allPlots_reduced.root",pathIn.c_str()),"READ");

	//----------------------------------------------- Generated J/Psi events ------------------------------------------------//
	TH2D *histGenJPsi = (TH2D*) fIn -> Get("histPtRapGen");
	TH2D *histGenJpsiPtRapClone = (TH2D*) CutTH2X(histGenJPsi, "histGenJPsi",0., 20.);
	TH2D *histGenJpsiPtRapClone2 = (TH2D*) CutTH2Y(histGenJPsi, "histGenJPsi", 2.5, 4);
	TH1D *histGenJpsiPt = (TH1D*) histGenJpsiPtRapClone2 -> ProjectionX("histGenJpsiPt");
	TH1D *histGenJpsiRap = (TH1D*) histGenJpsiPtRapClone2 -> ProjectionY("histGenJpsiRap");

	//-------------------------------------------- Reconstructed J/Psi events ----------------------------------------------//

    TH2D *histRecJPsi = (TH2D*) fIn -> Get("histPtRapRec");
	TH2D *histRecJpsiPtRapClone = (TH2D*) CutTH2X(histRecJPsi, "histRecJPsi",0., 20.);
	TH2D *histRecJpsiPtRapClone2 = (TH2D*) CutTH2Y(histRecJPsi, "histRecJPsi", 2.5, 4);
	TH1D *histRecJpsiPt = (TH1D*) histRecJpsiPtRapClone2 -> ProjectionX("histRecJpsiPt");
	TH1D *histRecJpsiRap = (TH1D*) histRecJpsiPtRapClone2 -> ProjectionY("histRecJpsiRap");

    /* histGenJpsiPtRap -> RebinX(2);
    histRecJpsiPtRap -> RebinX(2);

    histGenJpsiPtRap -> RebinY(5);
    histRecJpsiPtRap -> RebinY(2); */

	//------------------------------------------------ J/Psi Axe pT and rap -------------------------------------------------//

    double errRecJpsi, errGenJpsi;

    // Calcolo integrali con errore
    double intRecJpsi = histRecJpsiPtRapClone2->IntegralAndError(1, histRecJpsiPtRapClone2->GetNbinsX(), 1, histRecJpsiPtRapClone2->GetNbinsY(), errRecJpsi);
    double intGenJpsi = histGenJpsiPtRapClone2->IntegralAndError(1, histGenJpsiPtRapClone2->GetNbinsX(), 1, histGenJpsiPtRapClone2->GetNbinsY(), errGenJpsi);

    // Calcolo di Axe con propagazione degli errori
    double AxeJpsi = intRecJpsi / intGenJpsi;
    double errAxeJpsi = AxeJpsi * TMath::Sqrt(TMath::Power(errRecJpsi / intRecJpsi, 2) + TMath::Power(errGenJpsi / intGenJpsi, 2));

    // Stampa dei risultati
    std::cout << "AxeJpsi: " << AxeJpsi << " ± " << errAxeJpsi << std::endl;

    std::ofstream outFile(Form("%s/acceptance_efficiency.txt", pathOut.c_str()));
    if (outFile.is_open()) {
        outFile << AxeJpsi << "  " << errAxeJpsi << std::endl;
        outFile.close();
        std::cout << "Risultato salvato in 'acceptance_efficiency_result.txt'." << std::endl;
    } else {
        std::cerr << "Errore nell'aprire il file per la scrittura." << std::endl;
    }


    TH2D *histAxeJpsiPtRap = (TH2D*) histRecJpsiPtRapClone2 -> Clone("histAxeJpsiPtRap");
    histAxeJpsiPtRap -> Divide(histGenJpsiPtRapClone2);
	cout << "bins: " << histGenJpsiPt -> GetNbinsX() << endl;


	//------------------------------------------------ J/Psi Axe pT cuts -------------------------------------------------//

    TH1D *histGenJpsiPtRebin = (TH1D*) histGenJpsiPt -> Rebin(nBinsPt, "histGenJpsiPtRebin", ptBinsRun3); 
    cout << "bins after Gen: " << histGenJpsiPtRebin -> GetNbinsX() << endl;
    TH1D *histRecJpsiPtRebin = (TH1D*) histRecJpsiPt -> Rebin(nBinsPt, "histRecJpsiPtRebin", ptBinsRun3); 
    cout << "bins after Rec: " << histRecJpsiPtRebin -> GetNbinsX() << endl;

    TH1D *histAxeJpsiPt = new TH1D("histAxeJpsiPt", "", nBinsPt, ptBinsRun3);
    histAxeJpsiPt -> Divide(histRecJpsiPtRebin, histGenJpsiPtRebin, 1, 1, "B");

    for(int bin = 1; bin <= histAxeJpsiPt->GetNbinsX(); bin++) {
        cout << "Bin " << bin << " (pT range " << histAxeJpsiPt->GetXaxis()->GetBinLowEdge(bin) << "-" << histAxeJpsiPt->GetXaxis()->GetBinUpEdge(bin) << " GeV/c): " 
            << histAxeJpsiPt->GetBinContent(bin) << " ± " << histAxeJpsiPt->GetBinError(bin) << endl;
    }
    
    std::ofstream outFileJPsi(Form("%s/Axe_JPsi_Pt.txt", pathOut.c_str()));
    if (outFileJPsi.is_open()) {
        for(int binPt = 0; binPt < nBinsPt; binPt ++){
            outFileJPsi << histAxeJpsiPt->GetBinContent(binPt+1) << "  " << histAxeJpsiPt->GetBinError(binPt+1) << std::endl;
            //cout << "Axe err pT JPsi: " << histAxeJpsiPt->GetBinError(binPt) << endl;
        }
        /* outFileJPsi << AxeJpsi << "  " << errAxeJpsi << std::endl;
        outFileJPsi << AxePsi << "  " << errAxePsi << std::endl; */
        outFileJPsi.close();
        std::cout << "Risultato salvato in 'Axe_JPsi_result.txt'." << std::endl;
    } else {
        std::cerr << "Errore nell'aprire il file per la scrittura." << std::endl;
    }

	//------------------------------------------------- J/Psi Axe y cuts -------------------------------------------------//
    cout << "bins Rap: " << histGenJpsiRap -> GetNbinsX() << endl;
    TH1D *histGenJpsiRapRebin = (TH1D*) histGenJpsiRap -> Rebin(nBinsRap, "histGenJpsiRapRebin", rapBins); 
    TH1D *histRecJpsiRapRebin = (TH1D*) histRecJpsiRap -> Rebin(nBinsRap, "histRecJpsiRapRebin", rapBins);

    cout << "bins after Rap: " << histGenJpsiRapRebin -> GetNbinsX() << endl;

    TH1D *histAxeJpsiRap = new TH1D("histAxeJpsiRap", "", nBinsRap, rapBins);
    histAxeJpsiRap -> Divide(histRecJpsiRapRebin, histGenJpsiRapRebin, 1, 1, "B");

    std::ofstream outFileJPsiRap(Form("%s/Axe_JPsi_Rap.txt", pathOut.c_str()));
    if (outFileJPsiRap.is_open()) {
        for(int binRap = 0; binRap < nBinsRap; binRap ++){
            outFileJPsiRap << histAxeJpsiRap->GetBinContent(binRap+1) << "  " << histAxeJpsiRap->GetBinError(binRap+1) << std::endl;
        }
        /* outFileJPsiRap << AxeJpsi << "  " << errAxeJpsi << std::endl;
        outFileJPsiRap << AxePsi << "  " << errAxePsi << std::endl; */
        outFileJPsiRap.close();
        std::cout << "Risultato salvato in 'Axe_JPsi_result.txt'." << std::endl;
    } else {
        std::cerr << "Errore nell'aprire il file per la scrittura." << std::endl;
    }


	//-------------------------------------------------- Canvas creation --------------------------------------------------//

	//----------------------------------------------- J/Psi canvas creation -----------------------------------------------//

    TCanvas *canvasJpsi = new TCanvas("canvasJpsi", "", 1200, 1200);
    canvasJpsi -> Divide(2, 2);

    canvasJpsi -> cd(1);
    gPad -> SetLogy(1);
    histGenJpsiPtRebin -> SetTitle("");
    histGenJpsiPtRebin -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGenJpsiPtRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histGenJpsiPtRebin -> GetYaxis() -> SetRangeUser(1, 1e6);
    histGenJpsiPtRebin -> Draw("EP");
    histRecJpsiPtRebin -> Draw("EP SAME");

    canvasJpsi -> cd(2);
    gPad -> SetLogy(1);
    histGenJpsiRapRebin -> SetTitle("");
    histGenJpsiRapRebin -> GetXaxis() -> SetTitle("#it{y}");
    histGenJpsiRapRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histGenJpsiRapRebin -> GetYaxis() -> SetRangeUser(1e3, 1e6);
    histGenJpsiRapRebin -> Draw("EP");
    histRecJpsiRapRebin -> Draw("EP SAME");

    canvasJpsi -> cd(3);
    histAxeJpsiPt -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histAxeJpsiPt -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histAxeJpsiPt -> GetYaxis() -> SetRangeUser(0, 2);
    histAxeJpsiPt -> Draw("EP");

    canvasJpsi -> cd(4);
    histAxeJpsiRap -> GetXaxis() -> SetTitle("#it{y}");
    histAxeJpsiRap -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histAxeJpsiRap -> GetYaxis() -> SetRangeUser(0, 1);
    histAxeJpsiRap -> Draw("EP");

	//---------------------------------------------- J/Psi pT canvas creation ----------------------------------------------//

    TCanvas *canvasAxeJpsiPt = new TCanvas("canvasAxeJpsiPt", "", 800, 600);
    histAxeJpsiPt -> SetTitle("");
    histAxeJpsiPt -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histAxeJpsiPt -> GetYaxis() -> SetTitle("A#times#varepsilon_{J/#psi}");
    //histAxeJpsiPt -> GetXaxis() -> SetRangeUser(0, 20.);
    histAxeJpsiPt -> GetYaxis() -> SetRangeUser(0, 2);
	histAxeJpsiPt -> SetMarkerStyle(20);
	histAxeJpsiPt -> SetMarkerColor(kBlue);
    histAxeJpsiPt -> Draw("EP");

	//---------------------------------------------- J/Psi y canvas creation ----------------------------------------------//

    TCanvas *canvasAxeJpsiRap = new TCanvas("canvasAxeJpsiRap", "", 800, 600);
    histAxeJpsiRap -> SetTitle("");
    histAxeJpsiRap -> GetXaxis() -> SetTitle("#it{y}");
    histAxeJpsiRap -> GetYaxis() -> SetTitle("A#times#varepsilon_{J/#psi}");
    histAxeJpsiRap -> GetYaxis() -> SetRangeUser(0, 1.2);
	histAxeJpsiRap -> SetMarkerStyle(20);
	histAxeJpsiRap -> SetMarkerColor(kBlue);
    histAxeJpsiRap -> Draw("EP");

	

    canvasAxeJpsiPt -> SaveAs(Form("%s/Axe_Jpsi_pt.pdf", pathOut.c_str()));
    canvasAxeJpsiRap -> SaveAs(Form("%s/Axe_Jpsi_rap.pdf", pathOut.c_str()));


    TFile *fOut = new TFile(Form("%s/Axe.root", pathOut.c_str()), "RECREATE");
    histAxeJpsiPt -> Write();
    histAxeJpsiRap -> Write();
    canvasAxeJpsiPt -> Write();
    canvasAxeJpsiRap -> Write();
    fOut -> Close();
}
TH1D* ProjectTHnSparse(THnSparseD *histSparse, double minVar1, double maxVar1, double minVar2, double maxVar2, double minVar3, double maxVar3) {
    // Pt bins
    double minVar1Bin = histSparse -> GetAxis(0) -> FindBin(minVar1);
    double maxVar1Bin = histSparse -> GetAxis(0) -> FindBin(maxVar1 - 0.01);
    // Rapidity bins
    double minVar2Bin = histSparse -> GetAxis(1) -> FindBin(minVar2);
    double maxVar2Bin = histSparse -> GetAxis(1) -> FindBin(maxVar2 - 0.01);
    // Centrality bins
    double minVar3Bin = histSparse -> GetAxis(2) -> FindBin(minVar3);
    double maxVar3Bin = histSparse -> GetAxis(2) -> FindBin(maxVar3 - 0.01);

    histSparse -> GetAxis(0) -> SetRange(minVar1Bin, maxVar1Bin);
    histSparse -> GetAxis(1) -> SetRange(minVar2Bin, maxVar2Bin);
    histSparse -> GetAxis(2) -> SetRange(minVar3Bin, maxVar3Bin);

    Printf("%3.2f - %3.2f §§ %3.2f - %3.2f §§ %3.2f - %3.2f \n", minVar1Bin, maxVar1Bin, minVar2Bin, maxVar2Bin, minVar3Bin, maxVar3Bin);

    TH1D *histProj = (TH1D*) histSparse -> Projection(0, Form("histProj_%3.2f_%3.2f__%3.2f_%3.2f__%3.2f_%3.2f", minVar1, maxVar1, minVar2, maxVar2, minVar3, maxVar3));
    return histProj;
}
THnSparseD* CutTHnSparse(THnSparseD *histSparse, double minVar1, double maxVar1, double minVar2, double maxVar2, double minVar3, double maxVar3) {
    // Pt bins
    double minVar1Bin = histSparse -> GetAxis(0) -> FindBin(minVar1);
    double maxVar1Bin = histSparse -> GetAxis(0) -> FindBin(maxVar1 - 0.01);
    // Rapidity bins
    double minVar2Bin = histSparse -> GetAxis(1) -> FindBin(minVar2);
    double maxVar2Bin = histSparse -> GetAxis(1) -> FindBin(maxVar2 - 0.01);
    // Centrality bins
    double minVar3Bin = histSparse -> GetAxis(2) -> FindBin(minVar3);
    double maxVar3Bin = histSparse -> GetAxis(2) -> FindBin(maxVar3 - 0.01);

    histSparse -> GetAxis(0) -> SetRange(minVar1Bin, maxVar1Bin);
    histSparse -> GetAxis(1) -> SetRange(minVar2Bin, maxVar2Bin);
    histSparse -> GetAxis(2) -> SetRange(minVar3Bin, maxVar3Bin);

    Printf("%3.2f - %3.2f §§ %3.2f - %3.2f §§ %3.2f - %3.2f \n", minVar1Bin, maxVar1Bin, minVar2Bin, maxVar2Bin, minVar3Bin, maxVar3Bin);

    return histSparse;
}
TH2D* CutTH2Y(TH2D *hist2D, string title, double minCentr, double maxCentr) {
    double minCentrBin = hist2D -> GetYaxis() -> FindBin(minCentr);
    double maxCentrBin = hist2D -> GetYaxis() -> FindBin(maxCentr - 0.01);

    Printf("minPtBin = %1.0f, maxPtBin = %1.0f", minCentrBin, maxCentrBin);
    hist2D -> GetYaxis() -> SetRange(minCentrBin, maxCentrBin);

    TH2D *histProj = (TH2D*) hist2D -> Clone();
    return histProj;
}
TH2D* CutTH2X(TH2D *hist2D, string title, double minCentr, double maxCentr) {
    double minCentrBin = hist2D -> GetXaxis() -> FindBin(minCentr);
    double maxCentrBin = hist2D -> GetXaxis() -> FindBin(maxCentr - 0.01);

    Printf("minPtBin = %1.0f, maxPtBin = %1.0f", minCentrBin, maxCentrBin);
    hist2D -> GetXaxis() -> SetRange(minCentrBin, maxCentrBin);

    TH2D *histProj = (TH2D*) hist2D -> Clone();
    return histProj;
}
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.04);
}