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

void resolution_smearer(double smearing = 0.03) {
    TDatabasePDG *database = TDatabasePDG::Instance();
    int pdgCode = 13;
    double massMu = database -> GetParticle(pdgCode) -> Mass();

    uint64_t fSelection;
    float fMass, fPt, fEta, fTauz, fTauxy, fU2Q2, fCos2DeltaPhi, fR2EP, fR2SP, fCentFT0C = -99999;
    float fChi2pca, fSVertex, fEMC1, fEMC2, fPt1, fPt2, fPhi1, fPhi2, fEta1, fEta2, fPtMC1, fPtMC2, fPhiMC1, fPhiMC2, fPhi, fEtaMC1, fEtaMC2 = -99999;;
    int fSign, fSign1, fSign2 = -99999;
    UInt_t fMcDecision;
    Int_t fIsAmbig1, fIsAmbig2;

    TH1D *histMass = new TH1D("histMass", "", 120, 2, 5); histMass -> SetLineColor(kBlack);
    TH1D *histMassJpsiGen = new TH1D("histMassJpsiGen", "", 120, 2, 5); histMassJpsiGen -> SetLineColor(kRed);
    TH1D *histMassJpsiRec = new TH1D("histMassJpsiRec", "", 120, 2, 5); histMassJpsiRec -> SetLineColor(kBlue);
    TH1D *histMassJpsiSmearRec = new TH1D("histMassJpsiSmearRec", "", 120, 2, 5); histMassJpsiSmearRec -> SetLineColor(kMagenta);

    TH1D *histNormResoPx1 = new TH1D("histNormResoPx1", " ; (#it{p}_{x}^{rec} - #it{p}_{x}^{gen}) / #it{p}_{x}^{gen}", 500, -1, 1); histNormResoPx1 -> SetLineColor(kBlack);
    TH1D *histNormResoPy1 = new TH1D("histNormResoPy1", " ; (#it{p}_{y}^{rec} - #it{p}_{y}^{gen}) / #it{p}_{y}^{gen}", 500, -1, 1); histNormResoPy1 -> SetLineColor(kBlack);
    TH1D *histNormResoPz1 = new TH1D("histNormResoPz1", " ; (#it{p}_{z}^{rec} - #it{p}_{z}^{gen}) / #it{p}_{z}^{gen}", 500, -1, 1); histNormResoPz1 -> SetLineColor(kBlack);
    TH1D *histNormResoPx2 = new TH1D("histNormResoPx2", " ; (#it{p}_{x}^{rec} - #it{p}_{x}^{gen}) / #it{p}_{x}^{gen}", 500, -1, 1); histNormResoPx2 -> SetLineColor(kBlack);
    TH1D *histNormResoPy2 = new TH1D("histNormResoPy2", " ; (#it{p}_{y}^{rec} - #it{p}_{y}^{gen}) / #it{p}_{y}^{gen}", 500, -1, 1); histNormResoPy2 -> SetLineColor(kBlack);
    TH1D *histNormResoPz2 = new TH1D("histNormResoPz2", " ; (#it{p}_{z}^{rec} - #it{p}_{z}^{gen}) / #it{p}_{z}^{gen}", 500, -1, 1); histNormResoPz2 -> SetLineColor(kBlack);

    TH1D *histNormResoSmearPx1 = new TH1D("histNormResoSmearPx1", " ; (#it{p}_{x}^{rec} - #it{p}_{x}^{gen}) / #it{p}_{x}^{gen}", 500, -1, 1); histNormResoSmearPx1 -> SetLineColor(kRed);
    TH1D *histNormResoSmearPy1 = new TH1D("histNormResoSmearPy1", " ; (#it{p}_{y}^{rec} - #it{p}_{y}^{gen}) / #it{p}_{y}^{gen}", 500, -1, 1); histNormResoSmearPy1 -> SetLineColor(kRed);
    TH1D *histNormResoSmearPz1 = new TH1D("histNormResoSmearPz1", " ; (#it{p}_{z}^{rec} - #it{p}_{z}^{gen}) / #it{p}_{z}^{gen}", 500, -1, 1); histNormResoSmearPz1 -> SetLineColor(kRed);
    TH1D *histNormResoSmearPx2 = new TH1D("histNormResoSmearPx2", " ; (#it{p}_{x}^{rec} - #it{p}_{x}^{gen}) / #it{p}_{x}^{gen}", 500, -1, 1); histNormResoSmearPx2 -> SetLineColor(kRed);
    TH1D *histNormResoSmearPy2 = new TH1D("histNormResoSmearPy2", " ; (#it{p}_{y}^{rec} - #it{p}_{y}^{gen}) / #it{p}_{y}^{gen}", 500, -1, 1); histNormResoSmearPy2 -> SetLineColor(kRed);
    TH1D *histNormResoSmearPz2 = new TH1D("histNormResoSmearPz2", " ; (#it{p}_{z}^{rec} - #it{p}_{z}^{gen}) / #it{p}_{z}^{gen}", 500, -1, 1); histNormResoSmearPz2 -> SetLineColor(kRed);

    string pathToFile = "/Users/lucamicheletti/GITHUB/dq_run3_analyses/charmonia_production_PbPb/MC/data/AO2D.root";
    TFile *fIn = new TFile(pathToFile.c_str(), "READ");
    TIter next(fIn -> GetListOfKeys()); 
    TKey *key; 
    while ((key = (TKey*) next())) { 
        TString dirName = key -> GetName();
        if (!dirName.Contains("DF_")) {
            continue;
        }

        TTree *tree = (TTree*) fIn -> Get(Form("%s/O2rtdimuonall", dirName.Data()));
        tree -> SetBranchAddress("fSelection", &fSelection);
        tree -> SetBranchAddress("fMass", &fMass);
        tree -> SetBranchAddress("fPt", &fPt);
        tree -> SetBranchAddress("fEta", &fEta);
        tree -> SetBranchAddress("fPhi", &fPhi);
        tree -> SetBranchAddress("fSign", &fSign);
        tree -> SetBranchAddress("fPt1", &fPt1);
        tree -> SetBranchAddress("fEta1", &fEta1);
        tree -> SetBranchAddress("fPhi1", &fPhi1);
        tree -> SetBranchAddress("fSign1", &fSign1);
        tree -> SetBranchAddress("fPt2", &fPt2);
        tree -> SetBranchAddress("fEta2", &fEta2);
        tree -> SetBranchAddress("fPhi2", &fPhi2);
        tree -> SetBranchAddress("fSign2", &fSign2);
        tree -> SetBranchAddress("fPtMC1", &fPtMC1);
        tree -> SetBranchAddress("fEtaMC1", &fEtaMC1);
        tree -> SetBranchAddress("fPhiMC1", &fPhiMC1);
        tree -> SetBranchAddress("fPtMC2", &fPtMC2);
        tree -> SetBranchAddress("fEtaMC2", &fEtaMC2);
        tree -> SetBranchAddress("fPhiMC2", &fPhiMC2);
        tree -> SetBranchAddress("fMcDecision", &fMcDecision);
        tree -> SetBranchAddress("fIsAmbig1", &fIsAmbig1);
        tree -> SetBranchAddress("fIsAmbig2", &fIsAmbig2);
        tree -> SetBranchAddress("fChi2pca", &fChi2pca);

        for (int iEntry = 0;iEntry < tree -> GetEntries();iEntry++) {
            tree -> GetEntry(iEntry);
            if (!eventSelection(fSelection)) continue;
            if (fMcDecision < 1) continue;
            if (fSign != 0 || TMath::Abs(fEta) < 2.5 || TMath::Abs(fEta) > 4) continue;
            if (TMath::Abs(fEta1) < 2.5 || TMath::Abs(fEta1) > 4) continue;
            if (TMath::Abs(fEta2) < 2.5 || TMath::Abs(fEta2) > 4) continue;
            if (fSign1 == fSign2) continue;
            if (fChi2pca < 0) continue;

            TLorentzVector vecMuGen1;
            TLorentzVector vecMuGen2;
            TLorentzVector vecMuRec1;
            TLorentzVector vecMuRec2;

            vecMuGen1.SetPtEtaPhiM(fPtMC1, fEtaMC1, fPhiMC1, massMu);
            vecMuGen2.SetPtEtaPhiM(fPtMC2, fEtaMC2, fPhiMC2, massMu);
            vecMuRec1.SetPtEtaPhiM(fPt1, fEta1, fPhi1, massMu);
            vecMuRec2.SetPtEtaPhiM(fPt2, fEta2, fPhi2, massMu);

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

            TLorentzVector vecJpsiGen = vecMuGen1 + vecMuGen2;
            TLorentzVector vecJpsiRec = vecMuRec1 + vecMuRec2;

            double smearFactor = gRandom -> Gaus(0, smearing);
            double pxSmearRec1 = pxRec1 + (pxRec1 * smearFactor);
            double pySmearRec1 = pyRec1 + (pyRec1 * smearFactor);
            double pzSmearRec1 = pzRec1 + (pzRec1 * smearFactor);
            double pxSmearRec2 = pxRec2 + (pxRec2 * smearFactor);
            double pySmearRec2 = pyRec2 + (pyRec2 * smearFactor);
            double pzSmearRec2 = pzRec2 + (pzRec2 * smearFactor);

            histNormResoSmearPx1 -> Fill((pxSmearRec1 - pxGen1) / pxGen1);
            histNormResoSmearPy1 -> Fill((pySmearRec1 - pyGen1) / pyGen1);
            histNormResoSmearPz1 -> Fill((pzSmearRec1 - pzGen1) / pzGen1);

            histNormResoSmearPx2 -> Fill((pxSmearRec2 - pxGen2) / pxGen2);
            histNormResoSmearPy2 -> Fill((pySmearRec2 - pyGen2) / pyGen2);
            histNormResoSmearPz2 -> Fill((pzSmearRec2 - pzGen2) / pzGen2);

            ROOT::Math::PxPyPzMVector vecMuSmearRec1(pxSmearRec1, pySmearRec1, pzSmearRec1, massMu);
            ROOT::Math::PxPyPzMVector vecMuSmearRec2(pxSmearRec2, pySmearRec2, pzSmearRec2, massMu);

            auto vecJpsiSmearRec = vecMuSmearRec1 + vecMuSmearRec2;

            histMass -> Fill(fMass);
            histMassJpsiGen -> Fill(vecJpsiGen.M());
            histMassJpsiRec -> Fill(vecJpsiRec.M());
            histMassJpsiSmearRec -> Fill(vecJpsiSmearRec.M());
        }
    }
    fIn -> Close();

    TCanvas *canvasNormReso = new TCanvas("canvasNormReso", "", 1800, 1200);
    canvasNormReso -> Divide(3, 2);

    canvasNormReso -> cd(1);
    gPad -> SetLogy(true);
    histNormResoPx1 -> Draw("EP");
    histNormResoSmearPx1 -> Draw("EP SAME");
    
    canvasNormReso -> cd(2);
    gPad -> SetLogy(true);
    histNormResoPy1 -> Draw("EP");
    histNormResoSmearPy1 -> Draw("EP SAME");

    canvasNormReso -> cd(3);
    gPad -> SetLogy(true);
    histNormResoPz1 -> Draw("EP");
    histNormResoSmearPz1 -> Draw("EP SAME");

    canvasNormReso -> cd(4);
    gPad -> SetLogy(true);
    histNormResoPx2 -> Draw("EP");
    histNormResoSmearPx2 -> Draw("EP SAME");
    
    canvasNormReso -> cd(5);
    gPad -> SetLogy(true);
    histNormResoPy2 -> Draw("EP");
    histNormResoSmearPy2 -> Draw("EP SAME");

    canvasNormReso -> cd(6);
    gPad -> SetLogy(true);
    histNormResoPz2 -> Draw("EP");
    histNormResoSmearPz2 -> Draw("EP SAME");

    TCanvas *canvasMass = new TCanvas("canvasMass", "", 800, 600);
    gPad -> SetLogy(true);
    histMassJpsiGen -> Draw("EP");
    histMass -> Draw("EP SAME");
    histMassJpsiRec -> Draw("EP SAME");
    histMassJpsiSmearRec -> Draw("EP SAME");
}