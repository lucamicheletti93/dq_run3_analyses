#if !defined(CLING) || defined(ROOTCLING)

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <bitset>

#include "TSystemDirectory.h"
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TKey.h"
#include "THashList.h"
#include "TProfile.h"
#include "TTreeReader.h"
#include "TLine.h"
#include "TGaxis.h"
#include <ROOT/RDataFrame.hxx>

#endif

void LoadStyle();
void SetLegend(TLegend *);

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

inline void SetHistogram(TH1D *hist, Color_t color) {
    hist -> SetLineColor(color); 
    hist -> SetMarkerColor(color);
    hist -> SetMarkerStyle(20);
}

void read_mc_minitree(bool isPbPb = false, double ptCut = 0.0) {
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

    TH1D *histMassJpsiGen = new TH1D("histMassJpsiGen", " ; #it{m}_{#mu#mu} (GeV/c^{2})", 240, 2, 5); SetHistogram(histMassJpsiGen, kRed+1);
    TH1D *histMassJpsiGenFromRec = new TH1D("histMassJpsiGenFromRec", " ; #it{m}_{#mu#mu} (GeV/c^{2})", 240, 2, 5); SetHistogram(histMassJpsiGenFromRec, kMagenta);
    TH1D *histMassJpsiRec = new TH1D("histMassJpsiRec", " ; #it{m}_{#mu#mu} (GeV/c^{2})", 240, 2, 5); SetHistogram(histMassJpsiRec, kBlue+1);

    TH1D *histPtJpsiGen = new TH1D("histPtJpsiGen", " ; #it{p}_{T} (GeV/c)", 300, 0, 30); SetHistogram(histPtJpsiGen, kRed+1);
    TH1D *histPtJpsiRec = new TH1D("histPtJpsiRec", " ; #it{p}_{T} (GeV/c)", 300, 0, 30); SetHistogram(histPtJpsiRec, kBlue+1);

    TH1D *histRapJpsiGen = new TH1D("histRapJpsiGen", " ; #it{y}", 150, 2.5, 4); SetHistogram(histRapJpsiGen, kRed+1);
    TH1D *histRapJpsiRec = new TH1D("histRapJpsiRec", " ; #it{y}", 150, 2.5, 4); SetHistogram(histRapJpsiRec, kBlue+1);

    TH1D *histCentrJpsiGen = new TH1D("histCentrJpsiGen", " ; Centr (%)", 20, 0, 20); SetHistogram(histCentrJpsiGen, kRed+1);
    TH1D *histCentrJpsiRec = new TH1D("histCentrJpsiRec", " ; Centr (%)", 10, 0, 100); SetHistogram(histCentrJpsiRec, kBlue+1);

    TH1D *histNormResoPx1 = new TH1D("histNormResoPx1", " ; (#it{p}_{x}^{rec} - #it{p}_{x}^{gen}) / #it{p}_{x}^{gen}", 500, -1, 1); histNormResoPx1 -> SetLineColor(kBlack);
    TH1D *histNormResoPy1 = new TH1D("histNormResoPy1", " ; (#it{p}_{y}^{rec} - #it{p}_{y}^{gen}) / #it{p}_{y}^{gen}", 500, -1, 1); histNormResoPy1 -> SetLineColor(kBlack);
    TH1D *histNormResoPz1 = new TH1D("histNormResoPz1", " ; (#it{p}_{z}^{rec} - #it{p}_{z}^{gen}) / #it{p}_{z}^{gen}", 500, -1, 1); histNormResoPz1 -> SetLineColor(kBlack);
    TH1D *histNormResoPx2 = new TH1D("histNormResoPx2", " ; (#it{p}_{x}^{rec} - #it{p}_{x}^{gen}) / #it{p}_{x}^{gen}", 500, -1, 1); histNormResoPx2 -> SetLineColor(kBlack);
    TH1D *histNormResoPy2 = new TH1D("histNormResoPy2", " ; (#it{p}_{y}^{rec} - #it{p}_{y}^{gen}) / #it{p}_{y}^{gen}", 500, -1, 1); histNormResoPy2 -> SetLineColor(kBlack);
    TH1D *histNormResoPz2 = new TH1D("histNormResoPz2", " ; (#it{p}_{z}^{rec} - #it{p}_{z}^{gen}) / #it{p}_{z}^{gen}", 500, -1, 1); histNormResoPz2 -> SetLineColor(kBlack);

    string pathToFile = "/Users/lucamicheletti/alice/local_train_test_mc/reducedAO2D.root";
    TFile *fIn = new TFile(pathToFile.c_str(), "READ");
    TIter next(fIn -> GetListOfKeys()); 
    TKey *key; 
    while ((key = (TKey*) next())) { 
        TString dirName = key -> GetName();
        if (!dirName.Contains("DF_")) {
            continue;
        }

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
                histMassJpsiGen -> Fill(vecJpsiGen.M());
                histPtJpsiGen -> Fill(vecJpsiGen.Pt());
                histRapJpsiGen -> Fill(-vecJpsiGen.Rapidity());
                histCentrJpsiGen -> Fill(fImpactParameter);
            } else {
                ROOT::Math::PtEtaPhiMVector vecMuGen1(fPtMC1, fEtaMC1, fPhiMC1, massMu);
                ROOT::Math::PtEtaPhiMVector vecMuGen2(fPtMC2, fEtaMC2, fPhiMC2, massMu);
                auto vecJpsiGen = vecMuGen1 + vecMuGen2;
                if (TMath::Abs(vecJpsiGen.Rapidity()) > 4 || TMath::Abs(vecJpsiGen.Rapidity()) < 2.5) continue;
                histMassJpsiGen -> Fill(vecJpsiGen.M());
            }
        }

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

            double pxGen1 = vecMuGen1.Px();
            double pyGen1 = vecMuGen1.Py();
            double pzGen1 = vecMuGen1.Pz();
            double pxGen2 = vecMuGen2.Px();
            double pyGen2 = vecMuGen2.Py();
            double pzGen2 = vecMuGen2.Pz();

            double pxRec1 = vecMuRec1.Px();
            double pyRec1 = vecMuRec1.Py();
            double pzRec1 = vecMuRec1.Pz();
            double pxRec2 = vecMuRec2.Px();
            double pyRec2 = vecMuRec2.Py();
            double pzRec2 = vecMuRec2.Pz();

            histNormResoPx1 -> Fill((pxRec1 - pxGen1) / pxGen1);
            histNormResoPy1 -> Fill((pyRec1 - pyGen1) / pyGen1);
            histNormResoPz1 -> Fill((pzRec1 - pzGen1) / pzGen1);

            histNormResoPx2 -> Fill((pxRec2 - pxGen2) / pxGen2);
            histNormResoPy2 -> Fill((pyRec2 - pyGen2) / pyGen2);
            histNormResoPz2 -> Fill((pzRec2 - pzGen2) / pzGen2);

            histMassJpsiGenFromRec -> Fill(vecJpsiGen.M());
            histMassJpsiRec -> Fill(vecJpsiRec.M());
            histPtJpsiRec -> Fill(vecJpsiRec.Pt());
            histRapJpsiRec -> Fill(-vecJpsiRec.Rapidity());
            histCentrJpsiRec -> Fill(fCentFT0C);
        }
    }
    fIn -> Close();

    TCanvas *canvasNormReso = new TCanvas("canvasNormReso", "", 1800, 1200);
    canvasNormReso -> Divide(3, 2);

    canvasNormReso -> cd(1);
    gPad -> SetLogy(true);
    histNormResoPx1 -> Draw("EP");
    
    canvasNormReso -> cd(2);
    gPad -> SetLogy(true);
    histNormResoPy1 -> Draw("EP");

    canvasNormReso -> cd(3);
    gPad -> SetLogy(true);
    histNormResoPz1 -> Draw("EP");

    canvasNormReso -> cd(4);
    gPad -> SetLogy(true);
    histNormResoPx2 -> Draw("EP");
    
    canvasNormReso -> cd(5);
    gPad -> SetLogy(true);
    histNormResoPy2 -> Draw("EP");

    canvasNormReso -> cd(6);
    gPad -> SetLogy(true);
    histNormResoPz2 -> Draw("EP");

    TCanvas *canvasMass = new TCanvas("canvasMass", "", 800, 600);
    gPad -> SetLogy(true);
    histMassJpsiGen -> Draw("L");
    histMassJpsiGenFromRec -> Draw("EP SAME");
    histMassJpsiRec -> Draw("EP SAME");

    histPtJpsiGen -> Rebin(4);
    histPtJpsiRec -> Rebin(4);

    TH1D *histAxePtJpsi = new TH1D("histAxePtJpsi", " ; #it{p}_{T} (GeV/c)", histPtJpsiGen -> GetNbinsX(), 0, 30);
    histAxePtJpsi -> Divide(histPtJpsiRec, histPtJpsiGen, 1, 1, "B");
    SetHistogram(histAxePtJpsi, kBlack);

    TH1D *histAxeRapJpsi = new TH1D("histAxeRapJpsi", " ; #it{y}", histRapJpsiGen -> GetNbinsX(), 2.5, 4);
    histAxeRapJpsi -> Divide(histRapJpsiRec, histRapJpsiGen, 1, 1, "B");
    SetHistogram(histAxeRapJpsi, kBlack);

    TCanvas *canvasPtRap = new TCanvas("canvasPtRap", "", 1200, 1200);
    canvasPtRap -> Divide(2, 2);

    canvasPtRap -> cd(1);
    gPad -> SetLogy(true);
    histPtJpsiGen -> GetYaxis() -> SetRangeUser(0.1, 1e5);
    histPtJpsiGen -> Draw("EP");
    histPtJpsiRec -> Draw("EP SAME");

    canvasPtRap -> cd(2);
    gPad -> SetLogy(true);
    histRapJpsiGen -> GetYaxis() -> SetRangeUser(0.1, 1e5);
    histRapJpsiGen -> Draw("EP");
    histRapJpsiRec -> Draw("EP SAME");

    canvasPtRap -> cd(3);
    histAxePtJpsi -> GetYaxis() -> SetRangeUser(0, 1.2);
    histAxePtJpsi -> Draw("EP");

    canvasPtRap -> cd(4);
    histAxeRapJpsi -> GetYaxis() -> SetRangeUser(0, 1.2);
    histAxeRapJpsi -> Draw("EP");


    /*TCanvas *canvasCentr = new TCanvas("canvasCentr", "", 1200, 600);
    canvasCentr -> Divide(2, 1);

    canvasCentr -> cd(1);
    gPad -> SetLogy(true);
    histCentrJpsiGen -> Draw("EPL");

    canvasCentr -> cd(2);
    gPad -> SetLogy(true);
    histCentrJpsiRec -> Draw("EP SAME");*/


}