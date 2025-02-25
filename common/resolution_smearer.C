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
#include <TRandom3.h>
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

inline void SetHistogram(TH1D *hist, Color_t color, int lineWidth = 1, int markerStyle = 0) {
    hist -> SetLineColor(color); 
    hist -> SetLineWidth(lineWidth);
    hist -> SetMarkerColor(color);
    hist -> SetMarkerStyle(markerStyle);
    hist -> SetMarkerSize(0.8);
}

double FuncCB2(double *,double *);

//----------------------------------------------------------------------------------------------------//
void prepare_tree() {
    std::vector<std::pair<double,double>> ptBins = { 
        {0.0, 0.5}, {0.5, 1.0}, {1.0, 1.5}, {1.5, 2.0}, {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}, 
        {5.0, 6.0}, {6.0, 7.0}, {7.0, 8.0}, {8.0, 10.0}, {10.0, 12.0}, {12.0, 15.0}, {15.0, 20.0} 
    };

    TFile *fOut = new TFile("reducedAO2D_splitted.root", "RECREATE");

    string pathToFile = "/Users/lucamicheletti/alice/local_train_test_mc/reducedAO2D_merged.root";
    TFile *fIn = new TFile(pathToFile.c_str(), "READ");
    TIter next(fIn -> GetListOfKeys()); 
    TKey *key; 
    while ((key = (TKey*) next())) { 
        TString dirName = key -> GetName();
        if (!dirName.Contains("DF_")) {
            continue;
        }

        TTree *tree = (TTree*) fIn -> Get(Form("%s/O2rtdilmtreerec", dirName.Data()));

        for (size_t i = 0;i < ptBins.size();i++) {
            double low = ptBins[i].first;
            double high = ptBins[i].second;

            TString selection = Form("fPt > %f && fPt <= %f", low, high);
            std::cout << "Selezione per range " << low << " - " << high << ": " << selection << std::endl;
            TTree* treeSub = tree -> CopyTree(selection.Data());
            if (!treeSub) {
                std::cerr << "Errore nella copia dell'albero per il range " << low << " - " << high << std::endl;
                continue;
            }

            treeSub -> SetName(Form("tree_pt_%2.1f_%2.1f", low, high));
            fOut -> cd();
            treeSub -> Write();
        }
    }
    fOut -> Close();
    fIn -> Close();
}
//----------------------------------------------------------------------------------------------------//
void resolution_smearer(double smearingStep = 0.001, double minSmearing = 0.01, double maxSmearing = 0.03) {
    TDatabasePDG *database = TDatabasePDG::Instance();
    int muPdgCode = 13;
    int jpsiPdgCode = 443;
    double massMu = database -> GetParticle(muPdgCode) -> Mass();
    double massJpsi = database -> GetParticle(jpsiPdgCode) -> Mass();

    std::vector<std::pair<double,double>> vecPtBins = { 
        {0.0, 0.5}, {0.5, 1.0}, {1.0, 1.5}, {1.5, 2.0}, {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}, 
        {5.0, 6.0}, {6.0, 7.0}, {7.0, 8.0}, {8.0, 10.0}, {10.0, 12.0}, {12.0, 15.0}, {15.0, 20.0} 
    };
    const int nPtBins = 14;
    double ptBins[] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0};
    double refWidth[] = {0.076, 0.076, 0.076, 0.079, 0.080, 0.082, 0.084, 0.087, 0.090, 0.095, 0.099, 0.108, 0.112, 0.120};
    TH1D *histRefWidth = new TH1D("histRefWidth",  " ; #it{p}_{T} (GeV/c) ; #sigma (GeV/c)", nPtBins, ptBins); SetHistogram(histRefWidth, kBlack, 1, 20);
    TH1D *histWidth = new TH1D("histWidth",  " ; #it{p}_{T} (GeV/c) ; #sigma (GeV/c)", nPtBins, ptBins); SetHistogram(histWidth, kRed+1, 1, 20);
    TH1D *histWidthSmear = new TH1D("histWidthSmear",  " ; #it{p}_{T} (GeV/c) ; #sigma (GeV/c)", nPtBins, ptBins); SetHistogram(histWidthSmear, kRed+1, 1, 24);

    uint64_t fSelection;
    float fMass, fPt, fEta, fCentFT0C = -99999;
    float fPt1, fPt2, fPhi1, fPhi2, fEta1, fEta2, fPtMC1, fPtMC2, fPhiMC1, fPhiMC2, fPhi, fEtaMC1, fEtaMC2 = -99999;;
    UInt_t fMcDecision;

    string pathToFile = "reducedAO2D_splitted.root";
    TFile *fIn = new TFile(pathToFile.c_str(), "READ");

    TH1D *histMassJpsiRec[nPtBins];
    TH1D *histMassJpsiSmearRec[nPtBins];
    TF1 *funcMassJpsiRec[nPtBins];
    TF1 *funcMassJpsiSmearRec[nPtBins];

    for (size_t iPt = 0; iPt < vecPtBins.size();iPt++) {
        double width = 0;
        double errWidth = 0;
        double widthSmear = 0;
        double errWidthSmear = 0;
        double low = vecPtBins[iPt].first;
        double high = vecPtBins[iPt].second;

        histMassJpsiRec[iPt] = new TH1D(Form("histMassJpsiRec_Pt_%2.1f_%2.1f", low, high), Form("%2.1f < #it{p}_{T} < %2.1f GeV/c ; #it{m}_{#mu#mu} (GeV/c^{2})", low, high), 120, 2, 5); SetHistogram(histMassJpsiRec[iPt], kBlack, 1, 20);
        histMassJpsiSmearRec[iPt] = new TH1D(Form("histMassJpsiSmearRec_Pt_%2.1f_%2.1f", low, high), Form("%2.1f < #it{p}_{T} < %2.1f GeV/c ; #it{m}_{#mu#mu} (GeV/c^{2})", low, high), 120, 2, 5); SetHistogram(histMassJpsiSmearRec[iPt], kBlack, 1, 24);

        TTree *tree = (TTree*) fIn -> Get(Form("tree_pt_%2.1f_%2.1f", low, high));
        tree -> SetBranchAddress("fMcDecision", &fMcDecision);
        tree -> SetBranchAddress("fMass", &fMass);
        tree -> SetBranchAddress("fPt", &fPt);
        tree -> SetBranchAddress("fEta", &fEta);
        tree -> SetBranchAddress("fPhi", &fPhi);
        tree -> SetBranchAddress("fCentFT0C", &fCentFT0C);
        tree -> SetBranchAddress("fPtMC1", &fPtMC1);
        tree -> SetBranchAddress("fEtaMC1", &fEtaMC1);
        tree -> SetBranchAddress("fPhiMC1", &fPhiMC1);
        tree -> SetBranchAddress("fPtMC2", &fPtMC2);
        tree -> SetBranchAddress("fEtaMC2", &fEtaMC2);
        tree -> SetBranchAddress("fPhiMC2", &fPhiMC2);
        tree -> SetBranchAddress("fPt1", &fPt1);
        tree -> SetBranchAddress("fEta1", &fEta1);
        tree -> SetBranchAddress("fPhi1", &fPhi1);
        tree -> SetBranchAddress("fPt2", &fPt2);
        tree -> SetBranchAddress("fEta2", &fEta2);
        tree -> SetBranchAddress("fPhi2", &fPhi2);

        int iteration = 0;
        double smearing = minSmearing;
        histRefWidth -> SetBinContent(iPt+1, refWidth[iPt]);
        histRefWidth -> SetBinError(iPt+1, 0);
        while (TMath::Abs(widthSmear - refWidth[iPt]) > 0.005 || std::isnan(widthSmear) || std::isnan(errWidthSmear)) {
            histMassJpsiRec[iPt] -> Reset();
            histMassJpsiSmearRec[iPt] -> Reset();

            smearing += (iteration * smearingStep);
            for (int iEntry = 0;iEntry < tree -> GetEntries();iEntry++) {
                tree -> GetEntry(iEntry);
                if (fMcDecision < 1) continue;
                ROOT::Math::PtEtaPhiMVector vecMuRec1(fPt1, fEta1, fPhi1, massMu);
                ROOT::Math::PtEtaPhiMVector vecMuRec2(fPt2, fEta2, fPhi2, massMu);

                double pxRec1 = vecMuRec1.Px();
                double pyRec1 = vecMuRec1.Py();
                double pzRec1 = vecMuRec1.Pz();
                double pxRec2 = vecMuRec2.Px();
                double pyRec2 = vecMuRec2.Py();
                double pzRec2 = vecMuRec2.Pz();

                auto vecJpsiRec = vecMuRec1 + vecMuRec2;

                double smearFactor = gRandom -> Gaus(0, smearing);
                double pxSmearRec1 = pxRec1 + (pxRec1 * smearFactor);
                double pySmearRec1 = pyRec1 + (pyRec1 * smearFactor);
                double pzSmearRec1 = pzRec1 + (pzRec1 * smearFactor);
                double pxSmearRec2 = pxRec2 + (pxRec2 * smearFactor);
                double pySmearRec2 = pyRec2 + (pyRec2 * smearFactor);
                double pzSmearRec2 = pzRec2 + (pzRec2 * smearFactor);

                ROOT::Math::PxPyPzMVector vecMuSmearRec1(pxSmearRec1, pySmearRec1, pzSmearRec1, massMu);
                ROOT::Math::PxPyPzMVector vecMuSmearRec2(pxSmearRec2, pySmearRec2, pzSmearRec2, massMu);

                auto vecJpsiSmearRec = vecMuSmearRec1 + vecMuSmearRec2;

                histMassJpsiRec[iPt] -> Fill(vecJpsiRec.M());
                histMassJpsiSmearRec[iPt] -> Fill(vecJpsiSmearRec.M());
            }

            funcMassJpsiRec[iPt] = new TF1(Form("funcMassJpsiRec_Pt_%2.1f_%2.1f", low, high), FuncCB2, 2.5, 3.5, 7);
            funcMassJpsiRec[iPt] -> SetLineColor(kRed+1);
            funcMassJpsiRec[iPt] -> SetLineStyle(kSolid);
            funcMassJpsiRec[iPt] -> SetParameters(100, 3.103, 0.074, 1.232, 2.079, 3.114, 0.527);
            histMassJpsiRec[iPt] -> Fit(funcMassJpsiRec[iPt], "RQ0");

            funcMassJpsiSmearRec[iPt] = new TF1(Form("funcMassJpsiSmearRec_Pt_%2.1f_%2.1f", low, high), FuncCB2, 2.5, 3.5, 7);
            funcMassJpsiSmearRec[iPt] -> SetLineColor(kRed+1);
            funcMassJpsiSmearRec[iPt] -> SetLineStyle(kDashed);
            funcMassJpsiSmearRec[iPt] -> SetParameters(100, 3.103, 0.074, 1.232, 2.079, 3.114, 0.527);
            histMassJpsiSmearRec[iPt] -> Fit(funcMassJpsiSmearRec[iPt], "RQ0");

            width = funcMassJpsiRec[iPt] -> GetParameter(2);
            errWidth = funcMassJpsiRec[iPt] -> GetParError(2);
            widthSmear = funcMassJpsiSmearRec[iPt] -> GetParameter(2);
            errWidthSmear = funcMassJpsiSmearRec[iPt] -> GetParError(2);
            iteration++;

            if (smearing > maxSmearing) break;
        }
        histWidth -> SetBinContent(iPt+1, width);
        histWidth -> SetBinError(iPt+1, errWidth);
        histWidthSmear -> SetBinContent(iPt+1, widthSmear);
        histWidthSmear -> SetBinError(iPt+1, errWidthSmear);
    }

    TCanvas *canvasMassPt = new TCanvas("canvasMassPt", "", 3000, 1800);
    canvasMassPt -> Divide(5, 3);

    for (int iPt = 0;iPt < nPtBins;iPt++) {
        gPad -> SetLogy(true);
        canvasMassPt -> cd(iPt+1);
        histMassJpsiRec[iPt] -> SetStats(0);

        histMassJpsiRec[iPt] -> Draw("EP");
        histMassJpsiSmearRec[iPt] -> Draw("EP SAME");
        funcMassJpsiRec[iPt] -> Draw("L SAME");
        funcMassJpsiSmearRec[iPt] -> Draw("L SAME");

        TLatex latex1;
        latex1.SetNDC();
        latex1.SetTextSize(0.04);
        latex1.SetTextFont(42);
        latex1.DrawLatex(0.60, 0.85, Form("#mu = %.4f #pm %.4f", funcMassJpsiRec[iPt] -> GetParameter(1), funcMassJpsiRec[iPt] -> GetParError(1)));
        latex1.DrawLatex(0.60, 0.80, Form("#sigma = %.4f #pm %.4f", funcMassJpsiRec[iPt] -> GetParameter(2), funcMassJpsiRec[iPt] -> GetParError(2)));

        TLatex latex2;
        latex2.SetNDC();
        latex2.SetTextSize(0.04);
        latex2.SetTextFont(42);
        latex2.SetTextColor(kRed+1);
        latex2.DrawLatex(0.60, 0.75, Form("#mu = %.4f #pm %.4f", funcMassJpsiSmearRec[iPt] -> GetParameter(1), funcMassJpsiSmearRec[iPt] -> GetParError(1)));
        latex2.DrawLatex(0.60, 0.70, Form("#sigma = %.4f #pm %.4f", funcMassJpsiSmearRec[iPt] -> GetParameter(2), funcMassJpsiSmearRec[iPt] -> GetParError(2)));
    }

    TCanvas *canvasWidth = new TCanvas("canvasWidth", "", 800, 600);
    histRefWidth -> GetYaxis() -> SetRangeUser(0, 0.150);
    histRefWidth -> Draw("EP");
    histWidth -> Draw("EP SAME");
    histWidthSmear -> Draw("EP SAME");
}
//----------------------------------------------------------------------------------------------------//
void check_momentum_resolution(double smearing = 0.015) {
    TDatabasePDG *database = TDatabasePDG::Instance();
    int muPdgCode = 13;
    int jpsiPdgCode = 443;
    double massMu = database -> GetParticle(muPdgCode) -> Mass();
    double massJpsi = database -> GetParticle(jpsiPdgCode) -> Mass();

    //------------------------//
    // smearing configuration
    const int nPtBins = 14;
    double ptBins[] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0};
    TH2D *histMassPtJpsiRec = new TH2D("histMassPtJpsiRec", " ; #it{p}_{T} (GeV/c) ; #it{m}_{#mu#mu} (GeV/c^{2})", nPtBins, ptBins, 120, 2, 5);
    TH2D *histMassPtJpsiSmearRec = new TH2D("histMassPtJpsiSmearRec", " ; #it{p}_{T} (GeV/c) ; #it{m}_{#mu#mu} (GeV/c^{2})", nPtBins, ptBins, 120, 2, 5);

    uint64_t fSelection;
    float fMass, fPt, fEta, fCentFT0C = -99999;
    float fPt1, fPt2, fPhi1, fPhi2, fEta1, fEta2, fPtMC1, fPtMC2, fPhiMC1, fPhiMC2, fPhi, fEtaMC1, fEtaMC2 = -99999;;
    UInt_t fMcDecision;

    TH1D *histMassJpsiRec = new TH1D("histMassJpsiRec", " ; #it{m}_{#mu#mu} (GeV/c^{2})", 120, 2, 5); SetHistogram(histMassJpsiRec, kBlack, 1, 20);
    TH1D *histMassJpsiSmearRec = new TH1D("histMassJpsiSmearRec", " ; #it{m}_{#mu#mu} (GeV/c^{2})", 120, 2, 5); SetHistogram(histMassJpsiSmearRec, kBlack, 1, 24);

    TH1D *histNormResoPx1 = new TH1D("histNormResoPx1", " ; (#it{p}_{x}^{rec} - #it{p}_{x}^{gen}) / #it{p}_{x}^{gen}", 500, -1, 1); SetHistogram(histNormResoPx1, kBlack);
    TH1D *histNormResoPy1 = new TH1D("histNormResoPy1", " ; (#it{p}_{y}^{rec} - #it{p}_{y}^{gen}) / #it{p}_{y}^{gen}", 500, -1, 1); SetHistogram(histNormResoPy1, kBlack);
    TH1D *histNormResoPz1 = new TH1D("histNormResoPz1", " ; (#it{p}_{z}^{rec} - #it{p}_{z}^{gen}) / #it{p}_{z}^{gen}", 500, -1, 1); SetHistogram(histNormResoPz1, kBlack);
    TH1D *histNormResoPx2 = new TH1D("histNormResoPx2", " ; (#it{p}_{x}^{rec} - #it{p}_{x}^{gen}) / #it{p}_{x}^{gen}", 500, -1, 1); SetHistogram(histNormResoPx2, kBlack);
    TH1D *histNormResoPy2 = new TH1D("histNormResoPy2", " ; (#it{p}_{y}^{rec} - #it{p}_{y}^{gen}) / #it{p}_{y}^{gen}", 500, -1, 1); SetHistogram(histNormResoPy2, kBlack);
    TH1D *histNormResoPz2 = new TH1D("histNormResoPz2", " ; (#it{p}_{z}^{rec} - #it{p}_{z}^{gen}) / #it{p}_{z}^{gen}", 500, -1, 1); SetHistogram(histNormResoPz2, kBlack);

    TH1D *histNormResoSmearPx1 = new TH1D("histNormResoSmearPx1", " ; (#it{p}_{x}^{rec} - #it{p}_{x}^{gen}) / #it{p}_{x}^{gen}", 500, -1, 1); SetHistogram(histNormResoSmearPx1, kRed+1);
    TH1D *histNormResoSmearPy1 = new TH1D("histNormResoSmearPy1", " ; (#it{p}_{y}^{rec} - #it{p}_{y}^{gen}) / #it{p}_{y}^{gen}", 500, -1, 1); SetHistogram(histNormResoSmearPy1, kRed+1);
    TH1D *histNormResoSmearPz1 = new TH1D("histNormResoSmearPz1", " ; (#it{p}_{z}^{rec} - #it{p}_{z}^{gen}) / #it{p}_{z}^{gen}", 500, -1, 1); SetHistogram(histNormResoSmearPz1, kRed+1);
    TH1D *histNormResoSmearPx2 = new TH1D("histNormResoSmearPx2", " ; (#it{p}_{x}^{rec} - #it{p}_{x}^{gen}) / #it{p}_{x}^{gen}", 500, -1, 1); SetHistogram(histNormResoSmearPx2, kRed+1);
    TH1D *histNormResoSmearPy2 = new TH1D("histNormResoSmearPy2", " ; (#it{p}_{y}^{rec} - #it{p}_{y}^{gen}) / #it{p}_{y}^{gen}", 500, -1, 1); SetHistogram(histNormResoSmearPy2, kRed+1);
    TH1D *histNormResoSmearPz2 = new TH1D("histNormResoSmearPz2", " ; (#it{p}_{z}^{rec} - #it{p}_{z}^{gen}) / #it{p}_{z}^{gen}", 500, -1, 1); SetHistogram(histNormResoSmearPz2, kRed+1);

    string pathToFile = "/Users/lucamicheletti/alice/local_train_test_mc/reducedAO2D.root";
    TFile *fIn = new TFile(pathToFile.c_str(), "READ");
    TIter next(fIn -> GetListOfKeys()); 
    TKey *key; 
    while ((key = (TKey*) next())) { 
        TString dirName = key -> GetName();
        if (!dirName.Contains("DF_")) {
            continue;
        }

        TTree *tree = (TTree*) fIn -> Get(Form("%s/O2rtdilmtreerec", dirName.Data()));
        tree -> SetBranchAddress("fMcDecision", &fMcDecision);
        tree -> SetBranchAddress("fMass", &fMass);
        tree -> SetBranchAddress("fPt", &fPt);
        tree -> SetBranchAddress("fEta", &fEta);
        tree -> SetBranchAddress("fPhi", &fPhi);
        tree -> SetBranchAddress("fCentFT0C", &fCentFT0C);
        tree -> SetBranchAddress("fPtMC1", &fPtMC1);
        tree -> SetBranchAddress("fEtaMC1", &fEtaMC1);
        tree -> SetBranchAddress("fPhiMC1", &fPhiMC1);
        tree -> SetBranchAddress("fPtMC2", &fPtMC2);
        tree -> SetBranchAddress("fEtaMC2", &fEtaMC2);
        tree -> SetBranchAddress("fPhiMC2", &fPhiMC2);
        tree -> SetBranchAddress("fPt1", &fPt1);
        tree -> SetBranchAddress("fEta1", &fEta1);
        tree -> SetBranchAddress("fPhi1", &fPhi1);
        tree -> SetBranchAddress("fPt2", &fPt2);
        tree -> SetBranchAddress("fEta2", &fEta2);
        tree -> SetBranchAddress("fPhi2", &fPhi2);

        for (int iEntry = 0;iEntry < tree -> GetEntries();iEntry++) {
            tree -> GetEntry(iEntry);
            if (fMcDecision < 1) continue;
            ROOT::Math::PtEtaPhiMVector vecMuGen1(fPtMC1, fEtaMC1, fPhiMC1, massMu);
            ROOT::Math::PtEtaPhiMVector vecMuGen2(fPtMC2, fEtaMC2, fPhiMC2, massMu);
            ROOT::Math::PtEtaPhiMVector vecMuRec1(fPt1, fEta1, fPhi1, massMu);
            ROOT::Math::PtEtaPhiMVector vecMuRec2(fPt2, fEta2, fPhi2, massMu);

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

            auto vecJpsiGen = vecMuGen1 + vecMuGen2;
            auto vecJpsiRec = vecMuRec1 + vecMuRec2;

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

            histMassJpsiRec -> Fill(vecJpsiRec.M());
            histMassJpsiSmearRec -> Fill(vecJpsiSmearRec.M());

            histMassPtJpsiRec -> Fill(vecJpsiRec.Pt(), vecJpsiRec.M());
            histMassPtJpsiSmearRec -> Fill(vecJpsiSmearRec.Pt(), vecJpsiSmearRec.M());
        }
    }
    fIn -> Close();

    TCanvas *canvasNormReso = new TCanvas("canvasNormReso", "", 1800, 1200);
    canvasNormReso -> Divide(3, 2);

    canvasNormReso -> cd(1);
    gPad -> SetLogy(true);
    histNormResoPx1 -> Draw("HIST L");
    histNormResoSmearPx1 -> Draw("H SAME");
    
    canvasNormReso -> cd(2);
    gPad -> SetLogy(true);
    histNormResoPy1 -> Draw("HIST L");
    histNormResoSmearPy1 -> Draw("H SAME");

    canvasNormReso -> cd(3);
    gPad -> SetLogy(true);
    histNormResoPz1 -> Draw("HIST L");
    histNormResoSmearPz1 -> Draw("H SAME");

    canvasNormReso -> cd(4);
    gPad -> SetLogy(true);
    histNormResoPx2 -> Draw("HIST L");
    histNormResoSmearPx2 -> Draw("H SAME");
    
    canvasNormReso -> cd(5);
    gPad -> SetLogy(true);
    histNormResoPy2 -> Draw("HIST L");
    histNormResoSmearPy2 -> Draw("H SAME");

    canvasNormReso -> cd(6);
    gPad -> SetLogy(true);
    histNormResoPz2 -> Draw("HIST L");
    histNormResoSmearPz2 -> Draw("H SAME");


    TF1 *funcMassJpsiRec = new TF1("funcMassJpsiRec", FuncCB2, 2.5, 3.5, 7);
    funcMassJpsiRec -> SetLineColor(kRed+1);
    funcMassJpsiRec -> SetLineStyle(kSolid);
    funcMassJpsiRec -> SetParameters(1779, 3.103, 0.074, 1.232, 2.079, 3.114, 0.527);
    histMassJpsiRec -> Fit(funcMassJpsiRec, "RQ0");

    TF1 *funcMassJpsiSmearRec = new TF1("funcMassJpsiSmearRec", FuncCB2, 2.5, 3.5, 7);
    funcMassJpsiSmearRec -> SetLineColor(kRed+1);
    funcMassJpsiSmearRec -> SetLineStyle(kDashed);
    funcMassJpsiSmearRec -> SetParameters(1779, 3.103, 0.074, 1.232, 2.079, 3.114, 0.527);
    histMassJpsiSmearRec -> Fit(funcMassJpsiSmearRec, "RQ0");

    TCanvas *canvasMass = new TCanvas("canvasMass", "", 800, 600);
    gPad -> SetLogy(true);
    histMassJpsiRec -> SetStats(0);
    histMassJpsiRec -> GetXaxis() -> SetRangeUser(2.5, 4);
    histMassJpsiRec -> Draw("EP SAME");
    histMassJpsiSmearRec -> Draw("EP SAME");
    funcMassJpsiRec -> Draw("L SAME");
    funcMassJpsiSmearRec -> Draw("L SAME");

    TLatex latex1;
    latex1.SetNDC();
    latex1.SetTextSize(0.04);
    latex1.SetTextFont(42);
    latex1.DrawLatex(0.65, 0.85, Form("#mu = %.4f #pm %.4f", funcMassJpsiRec -> GetParameter(1), funcMassJpsiRec -> GetParError(1)));
    latex1.DrawLatex(0.65, 0.80, Form("#sigma = %.4f #pm %.4f", funcMassJpsiRec -> GetParameter(2), funcMassJpsiRec -> GetParError(2)));

    TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextSize(0.04);
    latex2.SetTextFont(42);
    latex2.SetTextColor(kRed+1);
    latex2.DrawLatex(0.65, 0.75, Form("#mu = %.4f #pm %.4f", funcMassJpsiSmearRec -> GetParameter(1), funcMassJpsiSmearRec -> GetParError(1)));
    latex2.DrawLatex(0.65, 0.70, Form("#sigma = %.4f #pm %.4f", funcMassJpsiSmearRec -> GetParameter(2), funcMassJpsiSmearRec -> GetParError(2)));

    // Project all histograms
    TH1D *histMassProj[nPtBins];
    for (int iPt = 0;iPt < nPtBins;iPt++) {
        histMassProj[iPt] = histMassPtJpsiRec -> ProjectionY(Form("histMassPtJpsiRec_Pt_%i", iPt), iPt+1, iPt+1);
        SetHistogram(histMassProj[iPt], kBlack, 1, 20);
    }

    TCanvas *canvasMassPt = new TCanvas("canvasMassPt", "", 3000, 1800);
    canvasMassPt -> Divide(5, 3);

    for (int iPt = 0;iPt < nPtBins;iPt++) {
        canvasMassPt -> cd(iPt+1);
        histMassProj[iPt] -> SetStats(0);
        histMassProj[iPt] -> Draw("EP");

        TF1 *funcMassJpsiRec = new TF1(Form("funcMassJpsiRec_Pt_%i", iPt), FuncCB2, 2.5, 3.5, 7);
        funcMassJpsiRec -> SetLineColor(kRed+1);
        funcMassJpsiRec -> SetLineStyle(kSolid);
        funcMassJpsiRec -> SetParameters(1779, 3.103, 0.074, 1.232, 2.079, 3.114, 0.527);
        histMassProj[iPt] -> Fit(funcMassJpsiRec, "RQ0");
        funcMassJpsiRec -> Draw("L SAME");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.SetTextFont(42);
        latex.DrawLatex(0.60, 0.85, Form("#mu = %.4f #pm %.4f", funcMassJpsiRec -> GetParameter(1), funcMassJpsiRec -> GetParError(1)));
        latex.DrawLatex(0.60, 0.80, Form("#sigma = %.4f #pm %.4f", funcMassJpsiRec -> GetParameter(2), funcMassJpsiRec -> GetParError(2)));
    }

    canvasMassPt -> cd(15);
    histMassPtJpsiRec -> Draw("COLZ");
}
///////////////////////////////////////////////////////
double FuncCB2(double *x,double *par)
{
  // Extended Crystal Ball : 7 parameters
  //
  // par[0] = Normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = alpha
  // par[4] = n
  // par[5] = alpha'
  // par[6] = n'

  double t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;

  double absAlpha = fabs((double)par[3]);
  double absAlpha2 = fabs((double)par[5]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core
  {
    return par[0]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail
  {
    double a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    double b = par[4]/absAlpha - absAlpha;
    return par[0]*(a/TMath::Power(b - t, par[4]));
  }

  if (t >= absAlpha2) //right tail
  {

    double c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
    double d = par[6]/absAlpha2 - absAlpha2;
    return par[0]*(c/TMath::Power(d + t, par[6]));
  }

  return 0.;
}