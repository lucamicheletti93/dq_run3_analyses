#if !defined(CLING) || defined(ROOTCLING)

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <bitset>

#include "TSystemDirectory.h"
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TKey.h"
#include "THashList.h"
#include "TProfile.h"
#include "TTreeReader.h"
#include "TLine.h"
#include "TGaxis.h"
#include <ROOT/RDataFrame.hxx>
#include "TH2F.h"
#include "TMinuit.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "THnSparse.h"
#include <TStyle.h>
#include "TSystem.h"
#include "TAxis.h"
#include "TBranch.h"
#include <vector>
#include <stdlib.h> 
#include <iterator>
#include <numeric>
#include <map>
#include "THStack.h"
#include <list>
#include <algorithm>
#include <cctype>
#include <string_view>
#include "Fit/FitResult.h"
#include "TEllipse.h"
#include <tuple>
#include <TLegend.h>
#include <TArrayD.h>

#endif

const int kNsel = 51;

// Taken from O2Physics/Common/CCDB/EventSelectionParams.cxx
const char* selectionLabels[kNsel] = {
    "kIsBBV0A", "kIsBBV0C", "kIsBBFDA", "kIsBBFDC", "kIsBBT0A", "kIsBBT0C",
    "kNoBGV0A", "kNoBGV0C", "kNoBGFDA", "kNoBGFDC", "kNoBGT0A", "kNoBGT0C",
    "kIsBBZNA", "kIsBBZNC", "kIsBBZAC", "kNoBGZNA", "kNoBGZNC", "kNoV0MOnVsOfPileup",
    "kNoSPDOnVsOfPileup", "kNoV0Casymmetry", "kIsGoodTimeRange", "kNoIncompleteDAQ",
    "kNoTPCLaserWarmUp", "kNoTPCHVdip", "kNoPileupFromSPD", "kNoV0PFPileup", "kNoSPDClsVsTklBG",
    "kNoV0C012vsTklBG", "kNoInconsistentVtx", "kNoPileupInMultBins", "kNoPileupMV",
    "kNoPileupTPC", "kIsTriggerTVX", "kIsINT1", "kNoITSROFrameBorder", "kNoTimeFrameBorder",
    "kNoSameBunchPileup", "kIsGoodZvtxFT0vsPV", "kIsVertexITSTPC", "kIsVertexTOFmatched",
    "kIsVertexTRDmatched", "kNoCollInTimeRangeNarrow", "kNoCollInTimeRangeStrict",
    "kNoCollInTimeRangeStandard", "kNoCollInTimeRangeVzDependent", "kNoCollInRofStrict",
    "kNoCollInRofStandard", "kNoHighMultCollInPrevRof", "kIsGoodITSLayer3",
    "kIsGoodITSLayer0123", "kIsGoodITSLayersAll"
};

constexpr uint64_t kTriggerMask =
    (1ULL << 37) |  // kIsGoodZvtxFT0vsPV
    (1ULL << 32) |  // kIsTriggerTVX
    (1ULL << 36) |  // kNoSameBunchPileup
    (1ULL << 34) |  // kNoITSROFrameBorder
    (1ULL << 35);   // kNoTimeFrameBorder

inline bool eventSelection(uint64_t selection) {
    return (selection & kTriggerMask) == kTriggerMask;
}

TH1D* ProjectTHnSparse(THnSparseD *, double , double , double , double , double , double );
THnSparseD* CutTHnSparse(THnSparseD *, double , double , double , double , double , double );
TH2D* CutTH2Y(TH2D *, string , double , double);
TH2D* CutTH2X(TH2D *, string , double , double);

TF1* PtJPsiPbPb5TeV_Func() {
    return new TF1("PtPsiPbPb5TeV", [](double* x, double* p) {
        Double_t px = x[0];
        Double_t p0 = p[0];
        Double_t p1 = p[1];
        Double_t p2 = p[2];
        Double_t p3 = p[3];

        return p0 * px / TMath::Power(1. + TMath::Power(px / p1, p2), p3);
    }, 0, 12, 4); // Range [0,20] con 4 parametri

}
static Double_t PtJPsiPbPb5TeV_tuned(const Double_t* px, const Double_t*){
    // jpsi pT in PbPb, tuned on data (2015) -> Castillo embedding https://alice.its.cern.ch/jira/browse/ALIROOT-8174?jql=text%20~%20%22LHC19a2%22
    Double_t x = *px;
    Float_t p0, p1, p2, p3;
    p0 = 1.00715e6;
    p1 = 3.50274;
    p2 = 1.93403;
    p3 = 3.96363;
    return p0 * x / TMath::Power(1. + TMath::Power(x / p1, p2), p3);
}

TF1* RapPsiPbPb5TeV_Func() {
    return new TF1("PtPsiPbPb5TeV", [](double* x, double* p) {
        Double_t px = x[0];
        Double_t p0 = p[0];
        Double_t p1 = p[1];
        Double_t p2 = p[2];
        Double_t p3 = p[3];

        return p0 * px / TMath::Power(1. + TMath::Power(px / p1, p2), p3);
    }, 2.5, 4, 4); // Range [0,20] con 4 parametri

}

inline void SetHistogram(TH1D *hist, Color_t color) {
    hist -> SetLineColor(color); 
    hist -> SetMarkerColor(color);
    hist -> SetMarkerStyle(20);
}

inline void SetLegend(TLegend *legend) {
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.045);
}


void iterative_procedure(double ptCut = 0.7) {
    bool isPbPb = false;
    TDatabasePDG *database = TDatabasePDG::Instance();
    int muPdgCode = 13;
    int jpsiPdgCode = 443;
    double massMu = database -> GetParticle(muPdgCode) -> Mass();
    double massJpsi = database -> GetParticle(jpsiPdgCode) -> Mass();

    uint64_t fSelection;
    float fMass, fPt, fEta, fTauz, fTauxy, fU2Q2, fCos2DeltaPhi, fR2EP, fR2SP, fCentFT0C, fImpactParameter = -99999;
    float fChi2pca, fSVertex, fEMC1, fEMC2, fPt1, fPt2, fPhi1, fPhi2, fEta1, fEta2, fPtMC1, fPtMC2, fPhiMC1, fPhiMC2, fPhi, fEtaMC1, fEtaMC2 = -99999;;
    int fSign, fSign1, fSign2 = -99999;
    UInt_t fMcDecision;
    Int_t fIsAmbig1, fIsAmbig2;

    const int nBinsPt = 8;
    const int nBinsRap = 6;
    const int nBinsCentr = 4;

    double ptBins[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0};
    double rapBins[] = {2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00};
    double centrBins[] = {0.0, 20.0, 40.0, 60.0, 90.0};
    

    TH1D *histPtJpsiGen = new TH1D("histPtJpsiGen", " ; #it{p}_{T} (GeV/c)", nBinsPt, ptBins); SetHistogram(histPtJpsiGen, kRed+1);
    TH1D *histRapJpsiGen = new TH1D("histRapJpsiGen", " ; #it{y}", nBinsRap, rapBins); SetHistogram(histRapJpsiGen, kRed+1);
    TH1D *histCentrJpsiGen = new TH1D("histCentrJpsiGen", " ; Centr (%)", nBinsCentr, centrBins); SetHistogram(histCentrJpsiGen, kRed+1);

    TH1D *histPtJpsiRec = new TH1D("histPtJpsiRec", " ; #it{p}_{T} (GeV/c)", nBinsPt, ptBins); SetHistogram(histPtJpsiRec, kBlue);
    TH1D *histRapJpsiRec = new TH1D("histRapJpsiRec", " ; #it{y}", nBinsRap, rapBins); SetHistogram(histRapJpsiRec, kBlue);
    TH1D *histCentrJpsiRec = new TH1D("histCentrJpsiRec", " ; Centr (%)", nBinsCentr, centrBins); SetHistogram(histCentrJpsiRec, kBlue);

    TH1D *histPtJpsiAxe = new TH1D("histPtJpsiAxe", " ; #it{p}_{T} (GeV/c)", nBinsPt, ptBins); SetHistogram(histPtJpsiAxe, kBlack);
    TH1D *histRapJpsiAxe = new TH1D("histRapJpsiAxe", " ; #it{y}", nBinsRap, rapBins); SetHistogram(histRapJpsiAxe, kBlack);
    TH1D *histCentrJpsiAxe = new TH1D("histCentrJpsiAxe", " ; Centr (%)", nBinsCentr, centrBins); SetHistogram(histCentrJpsiAxe, kBlack);

    string pathToFile = "/Users/lucamicheletti/Downloads/AO2D.root";
    TFile *fIn = new TFile(pathToFile.c_str(), "READ");
    TIter next(fIn -> GetListOfKeys()); 
    TKey *key; 
    while ((key = (TKey*) next())) { 
        TString dirName = key -> GetName();
        if (!dirName.Contains("DF_")) {
            continue;
        }

        // Generated Tree
        TTree *treeGen = (TTree*) fIn -> Get(Form("%s/O2rtdilmtreegen", dirName.Data()));
        treeGen -> SetBranchAddress("fMcDecision", &fMcDecision);
        treeGen -> SetBranchAddress("fImpactParameter", &fImpactParameter);
        treeGen -> SetBranchAddress("fPtMC1", &fPtMC1);
        treeGen -> SetBranchAddress("fEtaMC1", &fEtaMC1);
        treeGen -> SetBranchAddress("fPhiMC1", &fPhiMC1);
        treeGen -> SetBranchAddress("fPtMC2", &fPtMC2);
        treeGen -> SetBranchAddress("fEtaMC2", &fEtaMC2);
        treeGen -> SetBranchAddress("fPhiMC2", &fPhiMC2);

        for (int iEntry = 0;iEntry < treeGen -> GetEntries();iEntry++) {
            treeGen -> GetEntry(iEntry);
            if (fMcDecision < 1) continue;

            if (fPtMC2 == -999 && fPhiMC2 == -999 && fPhiMC2 == -999) {
                ROOT::Math::PtEtaPhiMVector vecJpsiGen(fPtMC1, fEtaMC1, fPhiMC1, massJpsi);
                if (TMath::Abs(vecJpsiGen.Rapidity()) > 4 || TMath::Abs(vecJpsiGen.Rapidity()) < 2.5) continue;
                if (vecJpsiGen.Pt() > 30) continue;
                histPtJpsiGen -> Fill(vecJpsiGen.Pt());
                histRapJpsiGen -> Fill(-vecJpsiGen.Rapidity());

                if(fImpactParameter < 5.625) histCentrJpsiGen -> AddBinContent(1);
                if(fImpactParameter >= 5.625 && fImpactParameter < 8.375) histCentrJpsiGen -> AddBinContent(2);
                if(fImpactParameter >= 8.375 && fImpactParameter < 10.625) histCentrJpsiGen -> AddBinContent(3);
                if(fImpactParameter >= 10.625 && fImpactParameter <= 13.875) histCentrJpsiGen -> AddBinContent(4);
            }
        }


        // Reconstructed tree
        TTree *treeRec = (TTree*) fIn -> Get(Form("%s/O2rtdilmtreerec", dirName.Data()));
        treeRec -> SetBranchAddress("fMcDecision", &fMcDecision);
        treeRec -> SetBranchAddress("fMass", &fMass);
        treeRec -> SetBranchAddress("fPt", &fPt);
        treeRec -> SetBranchAddress("fEta", &fEta);
        treeRec -> SetBranchAddress("fPhi", &fPhi);
        treeRec -> SetBranchAddress("fCentFT0C", &fCentFT0C);
        treeRec -> SetBranchAddress("fPtMC1", &fPtMC1);
        treeRec -> SetBranchAddress("fEtaMC1", &fEtaMC1);
        treeRec -> SetBranchAddress("fPhiMC1", &fPhiMC1);
        treeRec -> SetBranchAddress("fPtMC2", &fPtMC2);
        treeRec -> SetBranchAddress("fEtaMC2", &fEtaMC2);
        treeRec -> SetBranchAddress("fPhiMC2", &fPhiMC2);
        treeRec -> SetBranchAddress("fPt1", &fPt1);
        treeRec -> SetBranchAddress("fEta1", &fEta1);
        treeRec -> SetBranchAddress("fPhi1", &fPhi1);
        treeRec -> SetBranchAddress("fPt2", &fPt2);
        treeRec -> SetBranchAddress("fEta2", &fEta2);
        treeRec -> SetBranchAddress("fPhi2", &fPhi2);

        for (int iEntry = 0;iEntry < treeRec -> GetEntries();iEntry++) {
            treeRec -> GetEntry(iEntry);
            if (fMcDecision < 1) continue;
            ROOT::Math::PtEtaPhiMVector vecMuGen1(fPtMC1, fEtaMC1, fPhiMC1, massMu);
            ROOT::Math::PtEtaPhiMVector vecMuGen2(fPtMC2, fEtaMC2, fPhiMC2, massMu);
            ROOT::Math::PtEtaPhiMVector vecMuRec1(fPt1, fEta1, fPhi1, massMu);
            ROOT::Math::PtEtaPhiMVector vecMuRec2(fPt2, fEta2, fPhi2, massMu);

            auto vecJpsiGen = vecMuGen1 + vecMuGen2;
            ROOT::Math::PtEtaPhiMVector vecJpsiRec(fPt, fEta, fPhi, fMass);

            if (TMath::Abs(vecJpsiRec.Rapidity()) > 4 || TMath::Abs(vecJpsiRec.Rapidity()) < 2.5) continue;
            if (fPt1 < ptCut || fPt2 < ptCut) continue;
            if (fPt > 30) continue;

            histPtJpsiRec -> Fill(vecJpsiRec.Pt());
            histRapJpsiRec -> Fill(-vecJpsiRec.Rapidity());
            histCentrJpsiRec -> Fill(fCentFT0C);

        }

    }

    histPtJpsiAxe -> Divide(histPtJpsiRec, histPtJpsiGen, 1, 1, "B");
    histRapJpsiAxe -> Divide(histRapJpsiRec, histRapJpsiGen, 1, 1, "B");
    histCentrJpsiAxe -> Divide(histCentrJpsiRec, histCentrJpsiGen, 1, 1, "B");

    TCanvas *canvasSim = new TCanvas("canvasSim", "", 1800, 1200);
    canvasSim -> Divide(3, 2);
    canvasSim -> cd(1);
    histPtJpsiGen -> Draw("EP");
    histPtJpsiRec -> Draw("EP SAME");
    canvasSim -> cd(2);
    histRapJpsiGen -> Draw("EP");
    histRapJpsiRec -> Draw("EP SAME");
    canvasSim -> cd(3);
    histCentrJpsiGen -> Draw("EP");
    histCentrJpsiRec -> Draw("EP SAME");
    canvasSim -> cd(4);
    histPtJpsiAxe -> Draw("EP");
    canvasSim -> cd(5);
    histRapJpsiAxe -> Draw("EP");
    canvasSim -> cd(6);
    histCentrJpsiAxe -> Draw("EP");


    string pathToFits = "/Users/lucamicheletti/Downloads/Centr_0.0_20.0";
    TH1D *histPtJpsiData = new TH1D("histPtJpsiData", " ; #it{p}_{T} (GeV/c)", nBinsPt, ptBins); SetHistogram(histPtJpsiData, kBlue);
    TH1D *histRapJpsiData = new TH1D("histRapJpsiData", " ; #it{y}", nBinsRap, rapBins); SetHistogram(histRapJpsiData, kBlue);
    for (int iPt = 0;iPt < nBinsPt;iPt++) {
        TFile *fIn = new TFile(Form("%s/Pt_%.1f_%.1f/output__Mass_BkgSub_CentrFT0C_0_20_Pt_%.0f_%.0f_ME_CB2_VWG__2.2_4.5.root",pathToFits.c_str(), ptBins[iPt], ptBins[iPt+1], ptBins[iPt], ptBins[iPt+1]),"READ");
        TH1F *histDataPt = (TH1F*) fIn -> Get(Form("fit_results_CB2_VWG__2.2_4.5_Mass_BkgSub_CentrFT0C_0_20_Pt_%.0f_%.0f_ME",ptBins[iPt], ptBins[iPt+1]));
        histPtJpsiData -> SetBinContent(iPt+1, histDataPt -> GetBinContent(10));
        histPtJpsiData -> SetBinError(iPt+1, histDataPt -> GetBinError(10));
    }

    histPtJpsiData -> Scale(1 / histPtJpsiData -> Integral(), "WIDTH");
    TH1D *histPtJpsiDataCorr = (TH1D*) histPtJpsiData -> Clone("histPtJpsiDataCorr");
    histPtJpsiDataCorr -> Divide(histPtJpsiAxe);

    TF1 *fitFunctionPt = PtJPsiPbPb5TeV_Func();
    fitFunctionPt -> SetParameter(0, 2.37317e+06);
    fitFunctionPt -> SetParameter(1, 2.83941);
    fitFunctionPt -> SetParameter(2, 2.6687);
    fitFunctionPt -> SetParameter(3, 2.37032);
    fitFunctionPt -> SetLineColor(kRed+1);
    histPtJpsiDataCorr->Fit(fitFunctionPt, "R0"); 
    fitFunctionPt->Draw("SAME");

    TH1D *histFromFuncRecPt = (TH1D*) fitFunctionPt -> GetHistogram();
    histFromFuncRecPt -> SetName("histFromFuncRecPt");
    

    TF1 *fitFunctionPtOriginal = PtJPsiPbPb5TeV_Func();
    fitFunctionPtOriginal -> FixParameter(0, 2);
    fitFunctionPtOriginal -> FixParameter(1, 3.50274);
    fitFunctionPtOriginal -> FixParameter(2, 1.93403);
    fitFunctionPtOriginal -> FixParameter(3, 3.96363);
    fitFunctionPtOriginal -> SetLineColor(kBlue+1);
    fitFunctionPtOriginal->Draw("SAME");

    TH1D *histFromFuncGenPt = (TH1D*) fitFunctionPtOriginal -> GetHistogram();
    histFromFuncGenPt -> SetName("histFromFuncGenPt");

    TH1D *histFromFuncRatio = (TH1D*) histFromFuncRecPt -> Clone("histFromFuncRatio");
    histFromFuncRatio -> Divide(histFromFuncGenPt);
    histFromFuncRatio -> SetTitle("Weights ; #it{p}_{T} (GeV/c)");
    histFromFuncRatio -> GetYaxis() -> SetRangeUser(0., 2.);
    SetHistogram(histFromFuncRatio, kBlack);
    
    TCanvas *canvasData = new TCanvas("canvasData", "", 1800, 1200);
    canvasData -> Divide(3, 2);

    canvasData -> cd(1);
    gPad -> SetLogy(true);
    histPtJpsiData -> SetStats(0);
    histPtJpsiData -> SetTitle("Data ; #it{p}_{T} (GeV/c)");
    histPtJpsiData -> Draw("EP");

    canvasData -> cd(2);
    gPad -> SetLogy(true);
    histPtJpsiDataCorr -> SetStats(0);
    histPtJpsiDataCorr -> SetTitle("Data / Axe ; #it{p}_{T} (GeV/c)");
    histPtJpsiDataCorr -> Draw("EP");
    histFromFuncGenPt -> Draw("HIST SAME");
    histFromFuncRecPt -> Draw("HIST SAME");

    TLegend *legendFromFit = new TLegend(0.50, 0.55, 0.70, 0.75);
    SetLegend(legendFromFit);
    legendFromFit -> AddEntry(histPtJpsiDataCorr, "Data / Axe", "PL");
    legendFromFit -> AddEntry(histFromFuncGenPt, "Original Pt shape", "PL");
    legendFromFit -> AddEntry(histFromFuncRecPt, "0 iter Pt shape", "PL");
    legendFromFit -> Draw("SAME");

    canvasData -> cd(3);
    histFromFuncRatio -> Draw("P");
    TLine *lineUnity = new TLine(0, 1, 12, 1);
    lineUnity -> Draw();

}