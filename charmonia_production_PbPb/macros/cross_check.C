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

void cross_check() {
    TFile *fIn = new TFile("/Users/lucamicheletti/Downloads/AnalysisResults_cross_check.root", "READ");
    TList *list1 = (TList*) fIn -> Get("analysis-same-event-pairing/output");

    TList *listSEPM = (TList*) list1 -> FindObject("PairsMuonSEPM_muonLowPt510SigmaPDCA");
    TList *listSEMM = (TList*) list1 -> FindObject("PairsMuonSEMM_muonLowPt510SigmaPDCA");
    TList *listSEPP = (TList*) list1 -> FindObject("PairsMuonSEPP_muonLowPt510SigmaPDCA");

    TH2F* histMEPM = (TH2F*) listSEPM -> FindObject("Mass_CentFT0C");
    TH2F* histMEMM = (TH2F*) listSEMM -> FindObject("Mass_CentFT0C");
    TH2F* histMEPP = (TH2F*) listSEPP -> FindObject("Mass_CentFT0C");

    TH1F *histMEPM_ProjMass = (TH1F*) histMEPM -> ProjectionX("histMEPM_ProjMass");
    TH1F *histMEMM_ProjMass = (TH1F*) histMEMM -> ProjectionX("histMEMM_ProjMass");
    TH1F *histMEPP_ProjMass = (TH1F*) histMEPP -> ProjectionX("histMEPP_ProjMass");

    TH1F *histMEPPMM_ProjMass = new TH1F("histMEPPMM_ProjMass", "", 750, 0., 15.);

    double MM_count, PP_count = 0;
    int nBins = histMEPM_ProjMass -> GetNbinsX();
    for (int iBin = 0;iBin < nBins;iBin++) {
        MM_count = histMEMM_ProjMass -> GetBinContent(iBin+1);
        PP_count = histMEPP_ProjMass -> GetBinContent(iBin+1);
        histMEPPMM_ProjMass -> SetBinContent(iBin+1, 2 * TMath::Sqrt(MM_count * PP_count));
    }

    //histMEPM_ProjMass
    histMEMM_ProjMass -> Scale(1. / histMEMM_ProjMass -> Integral());
    histMEPP_ProjMass -> Scale(1. / histMEPP_ProjMass -> Integral());
    histMEPPMM_ProjMass -> Scale(1. / histMEPPMM_ProjMass -> Integral());

    TCanvas *canvas = new TCanvas("canvas", "", 800, 600);
    //histMEPM_ProjMass -> Draw("EP");
    histMEMM_ProjMass -> Draw("EP SAME");
    histMEPP_ProjMass -> Draw("EP SAME");
    histMEPPMM_ProjMass -> Draw("EP SAME");
}