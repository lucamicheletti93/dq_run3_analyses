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
    }, 0, 20, 4); // Range [0,20] con 4 parametri

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
TF1* RapPsiPbPb5TeV_Original() {
    return new TF1("RapPsiPbPb5TeV", [](double* x, double* p) {
        Double_t y = x[0];
        Double_t p0 = p[0];
        Double_t p1 = p[1];
        Double_t p2 = p[2];

        return p0 * TMath::Exp(-0.5 * TMath::Power((y - p1) / p2, 2));
    }, 2.5, 4, 3); // Range [-5,5] con 3 parametri
}

static Double_t YPsiPbPb5TeV(const Double_t* py, const Double_t* /*dummy*/){
    // jpsi y in PbPb, tuned on data (2015) -> Castillo embedding https://alice.its.cern.ch/jira/browse/ALIROOT-8174?jql=text%20~%20%22LHC19a2%22
    // found at the link  https://github.com/AliceO2Group/O2DPG/blob/master/MC/config/PWGDQ/external/generator/GeneratorCocktailPromptCharmoniaToMuonEvtGen_PbPb5TeV.C#L164
    Double_t y = *py;
    Float_t p0, p1, p2;
    p0 = 1.09886e6;
    p1 = 0;
    p2 = 2.12568;
    return p0 * TMath::Exp(-(1. / 2.) * TMath::Power(((y - p1) / p2), 2));
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


void iterative(double ptCut = 0.7) {
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

    const int nBinsPt = 9;
    const int nBinsRap = 6;
    const int nBinsCentr = 4;

    double ptBins[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 20.0};
    double rapBins[] = {2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00};
    double centrBins[] = {0.0, 20.0, 40.0, 60.0, 90.0};
    int iterColors[] = {kRed+1, kBlue+1, kGreen+2};

    TFile *fInJpsiSpectraData = new TFile("jpsi_spectra_data.root", "READ");
    TH1D *histPtJpsiData = (TH1D*) fInJpsiSpectraData -> Get("histPtJpsiData"); SetHistogram(histPtJpsiData, kBlack);
    TH1D *histRapJpsiData = (TH1D*) fInJpsiSpectraData -> Get("histRapJpsiData"); SetHistogram(histRapJpsiData, kBlack);

    histPtJpsiData -> Scale(1 / histPtJpsiData -> Integral(), "WIDTH");
    histRapJpsiData -> Scale(1 / histRapJpsiData -> Integral(), "WIDTH");

    const int nIterations = 3;
    TH1D *histPtJpsiDataCorr[nIterations+1];
    TH1D *histRapJpsiDataCorr[nIterations+1];

    TH1D *histPtJpsiGen[nIterations+1];
    TH1D *histRapJpsiGen[nIterations+1];
    TH1D *histCentrJpsiGen[nIterations+1];

    TH1D *histPtJpsiRec[nIterations+1];
    TH1D *histRapJpsiRec[nIterations+1];
    TH1D *histCentrJpsiRec[nIterations+1];

    TH1D *histPtJpsiAxe[nIterations+1];
    TH1D *histRapJpsiAxe[nIterations+1];
    TH1D *histCentrJpsiAxe[nIterations+1];

    TF1 *fitFunctionPtOriginal = PtJPsiPbPb5TeV_Func();
    fitFunctionPtOriginal -> FixParameter(0, 2);
    fitFunctionPtOriginal -> FixParameter(1, 3.50274);
    fitFunctionPtOriginal -> FixParameter(2, 1.93403);
    fitFunctionPtOriginal -> FixParameter(3, 3.96363);
    fitFunctionPtOriginal -> SetLineColor(kGray+2);
    TH1D *histFromFuncPtOriginal = (TH1D*) fitFunctionPtOriginal -> GetHistogram();
    SetHistogram(histFromFuncPtOriginal, kGray+2);

    TF1 *fitFunctionPt[nIterations+1];
    TH1D *histFromFuncPt[nIterations+1];
    TH1D *histFromFuncPtRatio[nIterations+1];

    TF1 *fitFunctionRapOriginal = RapPsiPbPb5TeV_Original();
    fitFunctionRapOriginal -> FixParameter(0, 8);
    fitFunctionRapOriginal -> FixParameter(1, 0);
    fitFunctionRapOriginal -> FixParameter(2,  2.12568);
    fitFunctionRapOriginal -> SetLineColor(kGray+2);
    TH1D *histFromFuncRapOriginal = (TH1D*) fitFunctionRapOriginal -> GetHistogram();
    SetHistogram(histFromFuncRapOriginal, kGray+2);

    TF1 *fitFunctionRap[nIterations+1];
    TH1D *histFromFuncRap[nIterations+1];
    TH1D *histFromFuncRapRatio[nIterations+1];

    string pathToFile = "AO2D.root";
    TFile *fIn = new TFile(pathToFile.c_str(), "READ");
    for (int iter = 0; iter < nIterations; iter ++) {
        std::cout << Form("************* Iteration %i *************", iter) << std::endl;
        // Initialize the corrected data ditribution to be fitted
        histPtJpsiDataCorr[iter] = (TH1D*) histPtJpsiData -> Clone(Form("histPtJpsiDataCorr_iter_%i", iter));
        SetHistogram(histPtJpsiDataCorr[iter], iterColors[iter]);
        //histPtJpsiDataCorr[iter] -> SetName(Form("histPtJpsiDataCorr_iter_%i", iter));
        histRapJpsiDataCorr[iter] = (TH1D*) histRapJpsiData -> Clone(Form("histRapJpsiDataCorr_iter_%i", iter));
        SetHistogram(histRapJpsiDataCorr[iter], iterColors[iter]);
        //histRapJpsiDataCorr[iter] -> SetName(Form("histRapJpsiDataCorr_iter_%i", iter));

        histPtJpsiGen[iter] = new TH1D(Form("histPtJpsiGen_iter_%i", iter), " ; #it{p}_{T} (GeV/c)", nBinsPt, ptBins); SetHistogram(histPtJpsiGen[iter], kRed+1);
        histRapJpsiGen[iter] = new TH1D(Form("histRapJpsiGen_iter_%i", iter), " ; #it{y}", nBinsRap, rapBins); SetHistogram(histRapJpsiGen[iter], kRed+1);
        histCentrJpsiGen[iter] = new TH1D(Form("histCentrJpsiGen_iter_%i", iter), " ; Centr (%)", nBinsCentr, centrBins); SetHistogram(histCentrJpsiGen[iter], kRed+1);

        histPtJpsiRec[iter] = new TH1D(Form("histPtJpsiRec_iter_%i", iter), " ; #it{p}_{T} (GeV/c)", nBinsPt, ptBins); SetHistogram(histPtJpsiRec[iter], kBlue);
        histRapJpsiRec[iter] = new TH1D(Form("histRapJpsiRec_iter_%i", iter), " ; #it{y}", nBinsRap, rapBins); SetHistogram(histRapJpsiRec[iter], kBlue);
        histCentrJpsiRec[iter] = new TH1D(Form("histCentrJpsiRec_iter_%i", iter), " ; Centr (%)", nBinsCentr, centrBins); SetHistogram(histCentrJpsiRec[iter], kBlue);

        histPtJpsiAxe[iter] = new TH1D(Form("histPtJpsiAxe_iter_%i", iter), " ; #it{p}_{T} (GeV/c)", nBinsPt, ptBins); SetHistogram(histPtJpsiAxe[iter], iterColors[iter]);
        histRapJpsiAxe[iter] = new TH1D(Form("histRapJpsiAxe_iter_%i", iter), " ; #it{y}", nBinsRap, rapBins); SetHistogram(histRapJpsiAxe[iter], iterColors[iter]);
        histCentrJpsiAxe[iter] = new TH1D(Form("histCentrJpsiAxe_iter_%i", iter), " ; Centr (%)", nBinsCentr, centrBins); SetHistogram(histCentrJpsiAxe[iter], iterColors[iter]);

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
                    if (vecJpsiGen.Pt() > 20) continue;
                    if(iter == 0){
                        histPtJpsiGen[iter] -> Fill(vecJpsiGen.Pt());
                        histRapJpsiGen[iter] -> Fill(-vecJpsiGen.Rapidity());
                    } else {
                        //cout << "QUI" << endl;
                        double binPt = histFromFuncPtRatio[iter-1] -> FindBin(vecJpsiGen.Pt());
                        double binRap = histFromFuncRapRatio[iter-1] -> FindBin(-vecJpsiGen.Rapidity());
                        //cout << "binPt : " << binPt << " binRap: " << binRap << endl;
                        double weightPt = histFromFuncPtRatio[iter-1] -> GetBinContent(binPt);
                        double weightRap = histFromFuncRapRatio[iter-1] -> GetBinContent(binRap);
                        double weightTot = weightPt*weightRap;
                        //cout << "weightPt : " << weightPt << " weightRap: " << weightRap << " weightTot: " << weightTot << endl;
                        histPtJpsiGen[iter] -> Fill(vecJpsiGen.Pt(), weightTot);
                        histRapJpsiGen[iter] -> Fill(-vecJpsiGen.Rapidity(), weightTot);
                    }

                    if(fImpactParameter < 5.625) histCentrJpsiGen[iter] -> AddBinContent(1);
                    if(fImpactParameter >= 5.625 && fImpactParameter < 8.375) histCentrJpsiGen[iter] -> AddBinContent(2);
                    if(fImpactParameter >= 8.375 && fImpactParameter < 10.625) histCentrJpsiGen[iter] -> AddBinContent(3);
                    if(fImpactParameter >= 10.625 && fImpactParameter <= 13.875) histCentrJpsiGen[iter] -> AddBinContent(4);
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
                if (fPt > 20) continue;
                if(iter == 0){
                    histPtJpsiRec[iter] -> Fill(vecJpsiRec.Pt());
                    histRapJpsiRec[iter] -> Fill(-vecJpsiRec.Rapidity());
                } else {
                    //cout << "QUI rec" << endl;
                    double binPt = histFromFuncPtRatio[iter-1] -> FindBin(vecJpsiGen.Pt());
                    double binRap = histFromFuncRapRatio[iter-1] -> FindBin(-vecJpsiGen.Rapidity());
                    double weightPt = histFromFuncPtRatio[iter-1] -> GetBinContent(binPt);
                    double weightRap = histFromFuncRapRatio[iter-1] -> GetBinContent(binRap);
                    double weightTot = weightPt*weightRap;
                    histPtJpsiRec[iter] -> Fill(vecJpsiRec.Pt(), weightTot);
                    histRapJpsiRec[iter] -> Fill(-vecJpsiRec.Rapidity(), weightTot);
                }
                histCentrJpsiRec[iter] -> Fill(fCentFT0C);
            }
        }

        histPtJpsiAxe[iter] -> Divide(histPtJpsiRec[iter], histPtJpsiGen[iter], 1, 1, "B");
        histRapJpsiAxe[iter] -> Divide(histRapJpsiRec[iter], histRapJpsiGen[iter], 1, 1, "B");
        histCentrJpsiAxe[iter] -> Divide(histCentrJpsiRec[iter], histCentrJpsiGen[iter], 1, 1, "B");

        // Divide by Axe
        histPtJpsiDataCorr[iter] -> Divide(histPtJpsiAxe[iter]);

        fitFunctionPt[iter] = PtJPsiPbPb5TeV_Func();
        fitFunctionPt[iter] -> SetParameters(1, 2.83941, 2.6687, 2.37032);
        fitFunctionPt[iter] -> SetLineColor(iterColors[iter]);
        histPtJpsiDataCorr[iter] -> Fit(fitFunctionPt[iter], "R0"); 

        histFromFuncPt[iter] = (TH1D*) fitFunctionPt[iter] -> GetHistogram();
        histFromFuncPt[iter] -> SetName(Form("histFromFuncPt_iter_%i", iter));

        histFromFuncPtRatio[iter] = (TH1D*) histFromFuncPt[iter] -> Clone(Form("histFromFuncPtRatio_iter_%i", iter));
        if (iter == 0) {
            histFromFuncPtRatio[iter] -> Divide(histFromFuncPtOriginal);
        } else {
            histFromFuncPtRatio[iter] -> Divide(histFromFuncPt[iter-1]);
        }
        histFromFuncPtRatio[iter] -> SetLineColor(iterColors[iter]);

        histRapJpsiDataCorr[iter] -> Divide(histRapJpsiAxe[iter]);

        fitFunctionRap[iter] = RapPsiPbPb5TeV_Func();
        fitFunctionRap[iter] -> SetParameters(1, 2.83941, 2.6687, 2.37032);
        fitFunctionRap[iter] -> SetLineColor(iterColors[iter]);
        histRapJpsiDataCorr[iter] -> Fit(fitFunctionRap[iter], "R0"); 

        histFromFuncRap[iter] = (TH1D*) fitFunctionRap[iter] -> GetHistogram();
        histFromFuncRap[iter] -> SetName(Form("histFromFuncRap_iter_%i", iter));

        histFromFuncRapRatio[iter] = (TH1D*) histFromFuncRap[iter] -> Clone(Form("histFromFuncRapRatio_iter_%i", iter));
        if (iter == 0) {
            histFromFuncRapRatio[iter] -> Divide(histFromFuncRapOriginal);
        } else {
            histFromFuncRapRatio[iter] -> Divide(histFromFuncRap[iter-1]);
        }
        histFromFuncRapRatio[iter] -> SetLineColor(iterColors[iter]);
    } // end of loop over iteratrions

    TLine *lineUnityPt = new TLine(0, 1, 20, 1);
    TLine *lineUnityRap = new TLine(2.5, 1, 4, 1);
    
    //------------------------------------//
    // Summary of the iterative procedure //
    TCanvas *canvasSummaryIterativeTuning = new TCanvas("canvasSummaryIterativeTuning", "", 1800, 1200);
    canvasSummaryIterativeTuning -> Divide(3, 2);

    canvasSummaryIterativeTuning -> cd(1);
    gPad -> SetLogy(true);
    histPtJpsiData -> SetStats(false);
    histPtJpsiData -> SetTitle("Data");
    histPtJpsiData -> Draw("EP");

    canvasSummaryIterativeTuning -> cd(2);
    gPad -> SetLogy(true);
    histPtJpsiDataCorr[0] -> SetStats(false);
    histPtJpsiDataCorr[0] -> SetTitle("Data / Ax#epsilon");

    TLegend *legendPtIterativeTuning = new TLegend(0.60, 0.55, 0.80, 0.75);
    SetLegend(legendPtIterativeTuning);
    legendPtIterativeTuning -> AddEntry(histFromFuncPtOriginal, "Original", "L");
    for(int iter = 0;iter < nIterations;iter++) {
        histPtJpsiDataCorr[iter] -> Draw("EP SAME");
        histFromFuncPt[iter] -> Draw("HIST SAME");
        legendPtIterativeTuning -> AddEntry(histFromFuncPt[iter], Form("Iteration %i", iter));
    }
    histFromFuncPtOriginal -> Draw("HIST SAME");
    legendPtIterativeTuning -> Draw("SAME");

    canvasSummaryIterativeTuning -> cd(3);
    histFromFuncPtRatio[0] -> SetStats(false);
    histFromFuncPtRatio[0] -> SetTitle("Weights");
    for(int iter = 0;iter < nIterations;iter++) {
        histFromFuncPtRatio[iter] -> GetYaxis() -> SetRangeUser(-1, 3);
        histFromFuncPtRatio[iter] -> Draw("HIST SAME");
    }
    lineUnityPt -> Draw();

    canvasSummaryIterativeTuning -> cd(4);
    histRapJpsiData -> SetStats(false);
    histRapJpsiData -> SetTitle("Data");
    histRapJpsiData -> Draw("EP");

    canvasSummaryIterativeTuning -> cd(5);
    histRapJpsiDataCorr[0] -> SetStats(false);
    histRapJpsiDataCorr[0] -> SetTitle("Data / Ax#epsilon");

    TLegend *legendRapIterativeTuning = new TLegend(0.60, 0.55, 0.80, 0.75);
    SetLegend(legendRapIterativeTuning);
    legendRapIterativeTuning -> AddEntry(histFromFuncRapOriginal, "Original", "L");
    for(int iter = 0;iter < nIterations;iter++) {
        histRapJpsiDataCorr[iter] -> Draw("EP SAME");
        histFromFuncRap[iter] -> Draw("HIST SAME");
        legendRapIterativeTuning -> AddEntry(histFromFuncRap[iter], Form("Iteration %i", iter));
    }
    histFromFuncRapOriginal -> Draw("HIST SAME");
    legendRapIterativeTuning -> Draw("SAME");

    canvasSummaryIterativeTuning -> cd(6);
    histFromFuncRapRatio[0] -> SetStats(false);
    histFromFuncRapRatio[0] -> SetTitle("Weights");
    for(int iter = 0;iter < nIterations;iter++) {
        histFromFuncRapRatio[iter] -> GetYaxis() -> SetRangeUser(-1, 3);
        histFromFuncRapRatio[iter] -> Draw("HIST SAME");
    }
    lineUnityRap -> Draw();

    //------------------------------------//
    // Summary of the iterative procedure //

    TCanvas *canvasSummaryAxe = new TCanvas("canvasSummaryAxe", "", 1200, 600);
    canvasSummaryAxe -> Divide(2, 1);

    canvasSummaryAxe -> cd(1);
    TLegend *legendPtAxe = new TLegend(0.20, 0.55, 0.60, 0.75);
    SetLegend(legendPtAxe);
    histPtJpsiAxe[0] -> SetStats(false);
    for(int iter = 0;iter < nIterations;iter++) {
        histPtJpsiAxe[iter] -> GetYaxis() -> SetRangeUser(0, 1);
        histPtJpsiAxe[iter] -> Draw("EP SAME");
        legendPtAxe -> AddEntry(histPtJpsiAxe[iter], Form("Iteration %i", iter), "L");
    }
    legendPtAxe -> Draw("SAME");

    canvasSummaryAxe -> cd(2);
    TLegend *legendRapAxe = new TLegend(0.60, 0.55, 0.80, 0.75);
    SetLegend(legendRapAxe);
    histRapJpsiAxe[0] -> SetStats(false);
    for(int iter = 0;iter < nIterations;iter++) {
        histRapJpsiAxe[iter] -> GetYaxis() -> SetRangeUser(0, 1);
        histRapJpsiAxe[iter] -> Draw("EP SAME");
        legendRapAxe -> AddEntry(histRapJpsiAxe[iter], Form("Iteration %i", iter), "L");
    }
    legendRapAxe -> Draw("SAME");


    canvasSummaryIterativeTuning -> SaveAs("SummaryIterativeTuning.pdf");
    canvasSummaryAxe -> SaveAs("SummaryAxe.pdf");
    
}