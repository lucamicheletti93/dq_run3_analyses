#if !defined(CLING) || defined(ROOTCLING)

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <bitset>

#include "TSystemDirectory.h"
#include <TLorentzVector.h>
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

void AO2D_analyzer() {
    //LoadStyle();
    uint64_t fSelection;
    float fMass, fPt, fEta, fTauz, fTauxy, fU2Q2, fCos2DeltaPhi, fR2EP, fR2SP, fCentFT0C = -99999;
    int fSign = -99999;

    const int nMassBins = 120;
    const double minMassRange = 2;
    const double maxMassRange = 5;

    TH2F *histMassPt = new TH2F("histMassPt", "", nMassBins, minMassRange, maxMassRange, 100, 0, 5);
    string pathIn = "/Users/saragaretti/cernbox/JPsi_Psi2S_ratio/MC/inputShapes_0.7";
    string pathPlots = "/Users/saragaretti/dq_fitter_tails/inputShapes/plots";

	TFile *fIn = new TFile(Form("%s/AO2D_merged.root",pathIn.c_str()),"READ");

    if (!fIn || fIn->IsZombie()) {
        cout << "Error: impossible to open the file AO2D_merged.root" << endl;
        return;
    }
    TTree* tree = (TTree*)fIn->Get("O2rtdimuonall");

    if (!tree) {
        cout << "Error: TTree 'O2rtdimuonall' not found!" << endl;
        return;
    }

    cout << "TTree 'O2rtdimuonall' found!" << endl;

    float fChi2pca, fSVertex, fEMC1, fEMC2, fPt1, fPt2, fPhi1, fPhi2, fEta1, fEta2, fPtMC1, fPtMC2, fPhiMC1, fPhiMC2, fPhi, fEtaMC1, fEtaMC2 = -99999;;
    UInt_t fMcDecision;
    //Int_t fIsAmbig1, fIsAmbig2;

    //--------------------------------------------------- Get Branches Addresses  ---------------------------------------------------//
    tree->SetBranchAddress("fChi2pca", &fChi2pca);
    tree->SetBranchAddress("fMcDecision", &fMcDecision);
    tree->SetBranchAddress("fSVertex", &fSVertex);
    tree->SetBranchAddress("fMass", &fMass);
    tree -> SetBranchAddress("fSelection", &fSelection);
    tree -> SetBranchAddress("fSign", &fSign);
    tree->SetBranchAddress("fEMC1", &fEMC1);
    tree->SetBranchAddress("fEMC2", &fEMC2);
    tree->SetBranchAddress("fPt1", &fPt1);
    tree->SetBranchAddress("fPt2", &fPt2);
    tree->SetBranchAddress("fPt", &fPt);
    tree->SetBranchAddress("fEta", &fEta);
    tree->SetBranchAddress("fPhi", &fPhi);
    tree->SetBranchAddress("fPhi1", &fPhi1);
    tree->SetBranchAddress("fPhi2", &fPhi2);
    tree->SetBranchAddress("fEta1", &fEta1);
    tree->SetBranchAddress("fEta2", &fEta2);
    tree->SetBranchAddress("fPt1", &fPt1);
    tree->SetBranchAddress("fPt2", &fPt2);
    tree->SetBranchAddress("fPtMC1", &fPtMC1);
    tree->SetBranchAddress("fPtMC2", &fPtMC2);
    tree->SetBranchAddress("fPhiMC1", &fPhiMC1);
    tree->SetBranchAddress("fPhiMC2", &fPhiMC2);
    tree->SetBranchAddress("fEtaMC1", &fEtaMC1);
    tree->SetBranchAddress("fEtaMC2", &fEtaMC2);

    //--------------------------------------------------- Histograms  ---------------------------------------------------//
    TH1F *hMass = new TH1F("hMass", "fMass", 100, 0, 10);
    TH2F *histPtRapRec = new TH2F("histPtRapRec", "p_T vs rapidity", 100, 0, 20, 100, 2.5, 4);
    TH2F *histPtRapGen1 = new TH2F("histPtRapGen1", "p_T vs rapidity", 100, 0, 20, 100, 2.5, 4);
    TH2F *histPtRapGen2 = new TH2F("histPtRapGen2", "p_T vs rapidity", 100, 0, 20, 100, 2.5, 4);
    TH2F *histPtRapDimuonGen = new TH2F("histPtRapDimuonGen", "p_T vs rapidity", 100, 0, 20, 100, 2.5, 4);

    TH1F *hPt1 = new TH1F("hPt1", "Pt Muon 1", 100, 0, 20);
    TH1F *hPt2 = new TH1F("hPt2", "Pt Muon 2", 100, 0, 20);
    TH1F *hPtDimuonRec = new TH1F("hPtDimuonRec", "Pt Dimuon", 100, 0, 20);
    TH1F *hRap1 = new TH1F("hRap1", "Rapidity Muon 1", 100, 2.5,4);
    TH1F *hRap2 = new TH1F("hRap2", "Rapidity Muon 2", 100, 2.5,4);
    TH1F *hEta1 = new TH1F("hEta1", "Eta Muon 1", 100, -4, -2.5);
    TH1F *hEta2 = new TH1F("hEta2", "Eta Muon 2", 100, -4, -2.5);
    TH1F *hRapDimuonRec = new TH1F("hRapDimuonRec", "Rapidity Dimuon Rec", 100, 2.5,4);
    TH1F *hEtaDimuonRec = new TH1F("hEtaDimuonRec", "Eta Dimuon Rec", 100, -4, -2.5);
    TH1F *hPhiDimuonRec = new TH1F("hPhiDimuonRec", "Phi Dimuon Rec", 50, -3, 3);
    TH1F *hRapDimuonGen = new TH1F("hRapDimuonGen", "Rapidity Dimuon Gen", 100, 2.5,4);
    TH1F *hEtaDimuonGen = new TH1F("hEtaDimuonGen", "Eta Dimuon Gen", 100, -4, -2.5);
    TH1F *hPhiDimuonGen = new TH1F("hPhiDimuonGen", "Phi Dimuon Gen", 50, -3, 3);
    TH1F *hRapRec = new TH1F("hRapRec", "Rapidity Rec", 100, 2.5,4);

    TH1F *hRapRatio = new TH1F("hRapRatio", "Rapidity ratio", 100, 2.5,4);
    TH1F *hPtRatio = new TH1F("hPtRatio", "Pt ratio", 100, 2.5,4);

    //--------------------------------------------------- Variables  ---------------------------------------------------//
    float px1, px2, py1, py2, pz1, pz2, rap1, rap2;
    float px1Rec, px2Rec, py1Rec, py2Rec, pz1Rec, pz2Rec;
    float rapRec, pxRec, pyRec, pzRec, fE;

    //--------------------------------------------------- Loop on tree entries  ---------------------------------------------------//

    for (int iEntry = 0;iEntry < tree -> GetEntries();iEntry++) {
        tree -> GetEntry(iEntry);
        //Apply mask for event selection
        if (!eventSelection(fSelection)) continue;
        //Select signal
        if (fMcDecision > 0 && fSign == 0) {
            //Apply cuts on Pt and Eta for Generated
            if (fPtMC1 > 0 && fPtMC1 < 20 && fPtMC2 > 0 && fPtMC2 < 20 && TMath::Abs(fEtaMC1) > 2.5 && TMath::Abs(fEtaMC1) < 4 && TMath::Abs(fEtaMC2) > 2.5 && TMath::Abs(fEtaMC2) < 4) {
                px1 = fPtMC1*TMath::Cos(fPhiMC1);
                px2 = fPtMC2*TMath::Cos(fPhiMC2);
                py1 = fPtMC1*TMath::Sin(fPhiMC1);
                py2 = fPtMC2*TMath::Sin(fPhiMC2);
                pz1 = TMath::Sqrt(pow(fEMC1,2) - pow(fPtMC1,2) - pow(0.1057,2));
                pz2 = TMath::Sqrt(pow(fEMC2,2) - pow(fPtMC2,2) - pow(0.1057,2));
                //cout << "fPtMC1: " << fPtMC1 << ", fPtMC2: " << fPtMC2 << ", fEMC1: " << fEMC1 << ", fEMC2: " << fEMC2 << endl;
                TLorentzVector muon1(px1, py1, pz1, fEMC1);
                TLorentzVector muon2(px2, py2, pz2, fEMC2);
                TLorentzVector dimuonGen = muon1 + muon2;
                hRapDimuonGen->Fill(dimuonGen.Rapidity());
                ROOT::Math::PtEtaPhiMVector muon1Eta(fPtMC1, fEtaMC1, fPhiMC1, 0.1057);
                ROOT::Math::PtEtaPhiMVector muon2Eta(fPtMC2, fEtaMC2, fPhiMC2, 0.1057);
                rap1 = muon1.Rapidity();
                rap2 = muon2.Rapidity();
                histPtRapGen1->Fill(fPtMC1, rap1);
                histPtRapGen2->Fill(fPtMC2, rap2);
                histPtRapDimuonGen->Fill(dimuonGen.Pt(), dimuonGen.Rapidity());
                hPt1->Fill(fPtMC1);
                hPt2->Fill(fPtMC2);
                hRap1->Fill(rap1);
                hRap2->Fill(rap2);
                hEta1->Fill(fEtaMC1);
                hEta2->Fill(fEtaMC2);
                ROOT::Math::PtEtaPhiMVector dimuonGenEta = muon1Eta + muon2Eta;
                hEtaDimuonGen->Fill(dimuonGenEta.Eta());
                hPhiDimuonGen->Fill(dimuonGenEta.Phi());
            }
            //Apply cuts on Pt and Eta for Reconstructed
            if (fPt > 0 && fPt < 20 && TMath::Abs(fEta) > 2.5 && TMath::Abs(fEta) < 4 && fPt1 > 0 && fPt1 < 20 && TMath::Abs(fEta1) > 2.5 && TMath::Abs(fEta1) < 4 && fPt2 > 0 && fPt2 < 20 && TMath::Abs(fEta2) > 2.5 && TMath::Abs(fEta2) < 4) {
                px1Rec = fPt1*TMath::Cos(fPhi1);
                px2Rec = fPt2*TMath::Cos(fPhi2);
                py1Rec = fPt1*TMath::Sin(fPhi1);
                py2Rec = fPt2*TMath::Sin(fPhi2);
                pz1Rec = fPt1 * TMath::SinH(TMath::Abs(fEta1));
                pz2Rec = fPt2 * TMath::SinH(TMath::Abs(fEta2));
                //cout << "pz2Rec: " << pz2Rec << endl;
                float E1Rec = TMath::Sqrt(px1Rec*px1Rec + py1Rec*py1Rec + pz1Rec*pz1Rec + 0.1057*0.1057);
                float E2Rec = TMath::Sqrt(px2Rec*px2Rec + py2Rec*py2Rec + pz2Rec*pz2Rec + 0.1057*0.1057);
                TLorentzVector muon1Rec(px1Rec, py1Rec, pz1Rec, E1Rec);
                TLorentzVector muon2Rec(px2Rec, py2Rec, pz2Rec, E2Rec);
                TLorentzVector dimuonRec = muon1Rec + muon2Rec;
                histMassPt -> Fill(fMass, fPt);
                hMass->Fill(fMass);
                ROOT::Math::PtEtaPhiMVector dimuonRecEta(fPt, fEta, fPhi, 0.2114);
                histPtRapRec->Fill(dimuonRec.Pt(), dimuonRec.Rapidity());
                hPtDimuonRec->Fill(dimuonRec.Pt());
                hEtaDimuonRec->Fill(dimuonRecEta.Eta());
                hRapDimuonRec->Fill(dimuonRec.Rapidity());
                hPhiDimuonRec->Fill(dimuonRec.Phi());
            }
        }
    }
    TCanvas *canvasMassPt = new TCanvas("canvasMassPt", "", 800, 600);
    histMassPt -> Draw("COLZ");
    canvasMassPt -> SaveAs(Form("%s/Mass_Pt.pdf", pathPlots.c_str()));
    /* tree->SetBranchAddress("fIsAmbig1", &fIsAmbig1);
    tree->SetBranchAddress("fIsAmbig2", &fIsAmbig2); */

    /* TObjArray *branches = tree->GetListOfBranches();

    for (int i = 0; i < branches->GetEntries(); i++) {
        TBranch *branch = (TBranch*)branches->At(i);
        cout << "Nome del branch: " << branch->GetName() << endl;

        TLeaf *leaf = branch->GetLeaf(branch->GetName());
        if (leaf) {
            cout << "Tipo del dato del branch " << branch->GetName() << ": " << leaf->GetTypeName() << endl;
        } else {
            cout << "Errore nel recuperare il leaf del branch" << endl;
        }
    }  */

    //--------------------------------------------------- Show plots statistics  ---------------------------------------------------//
    gStyle->SetOptStat(1111);
    hMass->SetStats(true);
    histMassPt->SetStats(true);
    hPt1->SetStats(true);
    hPt2->SetStats(true);
    hRap1->SetStats(true);
    hRap2->SetStats(true);
    hEta1->SetStats(true);
    hEta2->SetStats(true);
    hRapDimuonGen->SetStats(true);
    hEtaDimuonGen->SetStats(true);
    hPhiDimuonGen->SetStats(true);
    histPtRapGen1->SetStats(true);
    histPtRapGen2->SetStats(true);
    histPtRapDimuonGen->SetStats(true);
    histPtRapRec->SetStats(true);
    hPtDimuonRec->SetStats(true);
    hEtaDimuonRec->SetStats(true);
    hRapDimuonRec->SetStats(true);
    hPhiDimuonRec->SetStats(true);

    //--------------------------------------------------- Create Canvases  ---------------------------------------------------//

    TCanvas *c1 = new TCanvas("c1", "fMass", 800, 600);
    hMass->Draw("EP");
    c1->SaveAs(Form("%s/hMass.pdf", pathPlots.c_str()));

    TCanvas *c2 = new TCanvas("c2", "Pt Muon 1", 800, 600);
    hPt1->Draw("EP");
    c2->SaveAs(Form("%s/hPt1.pdf", pathPlots.c_str()));

    TCanvas *c3 = new TCanvas("c3", "Pt Muon 2", 800, 600);
    hPt2->Draw("EP");
    c3->SaveAs(Form("%s/hPt2.pdf", pathPlots.c_str()));

    TCanvas *c4 = new TCanvas("c4", "Pt Dimuon", 800, 600);
    hPtDimuonRec->Draw("EP");
    c4->SaveAs(Form("%s/hPtDimuonRec.pdf", pathPlots.c_str()));

    TCanvas *c5 = new TCanvas("c5", "Rapidity Muon 1", 800, 600);
    hRap1->Draw("EP");
    c5->SaveAs(Form("%s/hRap1.pdf", pathPlots.c_str()));

    TCanvas *c6 = new TCanvas("c6", "Rapidity Muon 2", 800, 600);
    hRap2->Draw("EP");
    c6->SaveAs(Form("%s/hRap2.pdf", pathPlots.c_str()));

    TCanvas *c7 = new TCanvas("c7", "Rapidity Dimuon Gen", 800, 600);
    hRapDimuonGen->Draw("EP");
    c7->SaveAs(Form("%s/hRapDimuonGen.pdf", pathPlots.c_str()));

    TCanvas *c8 = new TCanvas("c8", "PtRapRec", 800, 600);
    histPtRapRec->Draw("COLZ");
    c8->SaveAs(Form("%s/PtRapRec.pdf", pathPlots.c_str()));

    TCanvas *c9 = new TCanvas("c9", "PtRapGen1", 800, 600);
    histPtRapGen1->Draw("COLZ");
    c9->SaveAs(Form("%s/PtRapGen1.pdf", pathPlots.c_str()));

    TCanvas *c10 = new TCanvas("c10", "PtRapGen2", 800, 600);
    histPtRapGen2->Draw("COLZ");
    c10->SaveAs(Form("%s/PtRapGen2.pdf", pathPlots.c_str()));

    TCanvas *c11 = new TCanvas("c11", "Eta Muon 1", 800, 600);
    hEta1->Draw("EP");
    c11->SaveAs(Form("%s/hEta1.pdf", pathPlots.c_str()));

    TCanvas *c12 = new TCanvas("c12", "Eta Muon 2", 800, 600);
    hEta2->Draw("EP");
    c12->SaveAs(Form("%s/hEta2.pdf", pathPlots.c_str()));

    TCanvas *c13 = new TCanvas("c13", "Eta Gen dimuon", 800, 600);
    hEtaDimuonGen ->Draw("EP");
    c13->SaveAs(Form("%s/hEtaDimuonGen.pdf", pathPlots.c_str()));

    TCanvas *c14 = new TCanvas("c14", "Phi Gen dimuon", 800, 600);
    hPhiDimuonGen->Draw("EP");
    c14->SaveAs(Form("%s/hPhiDimuonGen.pdf", pathPlots.c_str()));

    TCanvas *c15 = new TCanvas("c15", "Eta Rec dimuon", 800, 600);
    hEtaDimuonRec->Draw("EP");
    c15->SaveAs(Form("%s/hEtaDimuonRec.pdf", pathPlots.c_str()));

    TCanvas *c16 = new TCanvas("c16", "Rap Rec dimuon", 800, 600);
    hRapDimuonRec->Draw("EP");
    c16->SaveAs(Form("%s/hRapDimuonRec.pdf", pathPlots.c_str()));

    TCanvas *c17 = new TCanvas("c17", "Phi Rec dimuon", 800, 600);
    hPhiDimuonRec->Draw("EP");
    c17->SaveAs(Form("%s/hPhiDimuonRec.pdf", pathPlots.c_str()));

    TCanvas *c18 = new TCanvas("c18", "PtRapGen", 800, 600);
    histPtRapDimuonGen->Draw("COLZ");
    c18->SaveAs(Form("%s/histPtRapDimuonGen.pdf", pathPlots.c_str()));

    //--------------------------------------------------- Save all plots in allPlots.root ---------------------------------------------------//

    TFile *fOut = new TFile(Form("%s/allPlots.root", pathPlots.c_str()), "RECREATE");
    histMassPt->Write("histMassPt");
    hMass->Write("hMass");
    hPt1->Write("hPt1");
    hPt2->Write("hPt2");
    hPtDimuonRec->Write("hPtDimuonRec");
    hRap1->Write("hRap1");
    hRap2->Write("hRap2");
    hEta1->Write("hEta1");
    hEta2->Write("hEta2");
    hRapDimuonGen->Write("hRapDimuonGen");
    hEtaDimuonGen->Write("hEtaDimuonGen");
    hPhiDimuonGen->Write("hPhiDimuonGen");
    histPtRapRec->Write("histPtRapRec");
    histPtRapGen1->Write("histPtRapGen1");
    histPtRapGen2->Write("histPtRapGen2");
    hEtaDimuonRec->Write("hEtaDimuonRec");
    hRapDimuonRec->Write("hRapDimuonRec");
    hPhiDimuonRec->Write("hPhiDimuonRec");
    histPtRapDimuonGen->Write("histPtRapDimuonGen");

    fOut->Close();

    histMassPt->SetOption("EP");
    hEtaDimuonRec->SetOption("EP");
    hPhiDimuonRec->SetOption("EP");
    hEtaDimuonGen->SetOption("EP");
    hPhiDimuonGen->SetOption("EP");

    //--------------------------------------------------- Save essential plots in allPlots_reduced.root ---------------------------------------------------//

    TFile *fOut2 = new TFile(Form("%s/allPlots_reduced.root", pathPlots.c_str()), "RECREATE");
    histMassPt->Write("histMassPt");
    histPtRapRec->Write("histPtRapRec");
    hEtaDimuonRec->Write("hEtaDimuonRec");
    hPhiDimuonRec->Write("hPhiDimuonRec");
    hEtaDimuonGen->Write("hEtaDimuonGen");
    hPhiDimuonGen->Write("hPhiDimuonGen");
    histPtRapDimuonGen->Write("histPtRapDimuonGen");

    fIn->Close();
    fOut2->Close();

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
    legend -> SetTextSize(0.045);
}