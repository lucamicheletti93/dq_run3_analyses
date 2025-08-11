#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TLegend.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <string>

void RetrieveTriggerInfo(TString , bool , string, double [11]);
void LoadStyle();
void SetLegend(TLegend *);
inline void SetHist(auto *hist, Color_t mkrCol = kBlack, int mkrSty = 20, double mkrSize = 1, Color_t lnCol = kBlack, int lnWidth = 1, int fillSty = 0, double alpha = 1) {
    hist -> SetMarkerColorAlpha(mkrCol, alpha);
    hist -> SetMarkerStyle(mkrSty);
    hist -> SetMarkerSize(mkrSize);
    hist -> SetLineColorAlpha(lnCol, alpha);
    hist -> SetLineWidth(lnWidth);
    hist -> SetFillStyle(fillSty);
}

void check_normalization() {
    LoadStyle();

    string fInLumiNameSkim2024StdAssoc = "data/2024/LHC24af_fDiMuon_std_assoc_luminosity.root";
    string fInNameSkim2024StdAssoc = "data/2024/LHC24af_fDiMuon_std_assoc.root";

    string fInLumiNameSkim2023StdAssoc = "data/2023/LHC23zs_fDiMuon_std_assoc_luminosity.root";
    string fInNameSkim2023StdAssoc = "data/2023/LHC23zs_fDiMuon_std_assoc.root";

    string fInLumiNameSkim2024TimeAssoc = "data/2024/LHC24af_fDiMuon_time_assoc_luminosity.root";
    string fInNameSkim2024TimeAssoc = "data/2024/LHC24af_fDiMuon_time_assoc.root";

    string fInLumiNameSkim2023TimeAssoc = "data/2023/LHC23zs_fDiMuon_time_assoc_luminosity.root";
    string fInNameSkim2023TimeAssoc = "data/2023/LHC23zs_fDiMuon_time_assoc.root";

    // Skimmed 2024, std. assoc.
    TFile *fInLumiSkim2024StdAssoc = TFile::Open(fInLumiNameSkim2024StdAssoc.c_str());
    TH1D *hist2024LumiSkim2024StdAssoc = (TH1D*) fInLumiSkim2024StdAssoc -> Get("histLumiSummary");
    double bcSelEffSkim2024StdAssoc = hist2024LumiSkim2024StdAssoc -> GetBinContent(hist2024LumiSkim2024StdAssoc -> GetXaxis() -> FindBin("bcCounterEfficiency"));
    double nEvtsSkim2024StdAssoc = hist2024LumiSkim2024StdAssoc -> GetBinContent(hist2024LumiSkim2024StdAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double lumiSkim2024StdAssoc = hist2024LumiSkim2024StdAssoc -> GetBinContent(hist2024LumiSkim2024StdAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));
    double collsAfterCuts = hist2024LumiSkim2024StdAssoc -> GetBinContent(hist2024LumiSkim2024StdAssoc -> GetXaxis() -> FindBin("collsAfterCuts"));

    TFile *fInSkim2024StdAssoc = TFile::Open(fInNameSkim2024StdAssoc.c_str());
    // Events histograms
    TList *listSkim2024StdAssocEvSel = (TList*) fInSkim2024StdAssoc -> Get("analysis-event-selection/output");
    TList *listSkim2024StdAssocEvBefCuts = (TList*) listSkim2024StdAssocEvSel -> FindObject("Event_BeforeCuts");
    TList *listSkim2024StdAssocEvAftCuts = (TList*) listSkim2024StdAssocEvSel -> FindObject("Event_AfterCuts");
    TH1D *histSkim2024StdAssocEvBefCuts = (TH1D*) listSkim2024StdAssocEvBefCuts -> FindObject("VtxZ");
    TH1D *histSkim2024StdAssocEvAftCuts = (TH1D*) listSkim2024StdAssocEvAftCuts -> FindObject("VtxZ");
    double evSelCorrSkim2024StdAssoc = (double) histSkim2024StdAssocEvBefCuts -> GetEntries() / (double) histSkim2024StdAssocEvAftCuts -> GetEntries();
    double collSelCorrSkim2024StdAssoc = (double) histSkim2024StdAssocEvBefCuts -> GetEntries() / collsAfterCuts;
    // Dimuon histograms
    TList *listSkim2024StdAssoc = (TList*) fInSkim2024StdAssoc -> Get("analysis-same-event-pairing/output");
    TList *listSkim2024StdAssocSE = (TList*) listSkim2024StdAssoc -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA");
    THnSparseD *histSparseSkim2024StdAssoc = (THnSparseD*) listSkim2024StdAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseSkim2024StdAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseSkim2024StdAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntSkim2024StdAssoc = (TH1D*) histSparseSkim2024StdAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");
    double countsSkimStdAssoc = histProjIntSkim2024StdAssoc -> Integral();
    histProjIntSkim2024StdAssoc -> Scale((evSelCorrSkim2024StdAssoc) / (bcSelEffSkim2024StdAssoc * lumiSkim2024StdAssoc));
    SetHist(histProjIntSkim2024StdAssoc, kRed+1, 20, 0.8, kRed+1);

    std::cout << "---------------------------" << std::endl;
    std::cout << "[Skimmed 2024, std. assoc.]" << std::endl;
    std::cout << "eff. ev. sel.    = " << evSelCorrSkim2024StdAssoc << std::endl;
    std::cout << "eff. job         = " << collSelCorrSkim2024StdAssoc << std::endl;
    std::cout << "Luminosity       = " << lumiSkim2024StdAssoc * bcSelEffSkim2024StdAssoc << std::endl;
    std::cout << "N. evts.         = " << nEvtsSkim2024StdAssoc * bcSelEffSkim2024StdAssoc << std::endl;
    std::cout << "Luminosity       = " << lumiSkim2024StdAssoc * bcSelEffSkim2024StdAssoc << std::endl;
    std::cout << "N. evts. corr.   = " << nEvtsSkim2024StdAssoc * bcSelEffSkim2024StdAssoc << std::endl;
    std::cout << "Luminosity corr. = " << lumiSkim2024StdAssoc * bcSelEffSkim2024StdAssoc << std::endl;
    std::cout << "---------------------------" << std::endl;

    // Skimmed 2024, time assoc.
    // WARNING! for the time association you have to take into account the duplication using zorro
    TFile *fInLumiSkim2024TimeAssoc = TFile::Open(fInLumiNameSkim2024TimeAssoc.c_str());
    TH1D *histLumiSkim2024TimeAssoc = (TH1D*) fInLumiSkim2024TimeAssoc -> Get("histLumiSummary");
    double bcSelEffSkim2024TimeAssoc = histLumiSkim2024TimeAssoc -> GetBinContent(histLumiSkim2024TimeAssoc -> GetXaxis() -> FindBin("bcCounterEfficiency"));
    double nEvtsSkim2024TimeAssoc = histLumiSkim2024TimeAssoc -> GetBinContent(histLumiSkim2024TimeAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double nEvtsAfterCutsSkim2024TimeAssoc = histLumiSkim2024TimeAssoc -> GetBinContent(histLumiSkim2024TimeAssoc -> GetXaxis() -> FindBin("nEvtsBcSelAfterCuts"));
    double lumiSkim2024TimeAssoc = histLumiSkim2024TimeAssoc -> GetBinContent(histLumiSkim2024TimeAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));

    TFile *fInSkim2024TimeAssoc = TFile::Open(fInNameSkim2024TimeAssoc.c_str());
    // Events histograms
    TList *listSkim2024TimeAssocEvSel = (TList*) fInSkim2024TimeAssoc -> Get("analysis-event-selection/output");
    TList *listSkim2024TimeAssocEvBefCuts = (TList*) listSkim2024TimeAssocEvSel -> FindObject("Event_BeforeCuts");
    TList *listSkim2024TimeAssocEvAftCuts = (TList*) listSkim2024TimeAssocEvSel -> FindObject("Event_AfterCuts");
    TH1D *histSkim2024TimeAssocEvBefCuts = (TH1D*) listSkim2024TimeAssocEvBefCuts -> FindObject("VtxZ");
    TH1D *histSkim2024TimeAssocEvAftCuts = (TH1D*) listSkim2024TimeAssocEvAftCuts -> FindObject("VtxZ");
    double evSelCorrSkim2024TimeAssoc = (double) histSkim2024TimeAssocEvBefCuts -> GetEntries() / (double) histSkim2024TimeAssocEvAftCuts -> GetEntries();
    // Dimuon histograms
    TList *listSkim2024TimeAssoc = (TList*) fInSkim2024TimeAssoc -> Get("analysis-same-event-pairing/output");
    TList *listSkim2024TimeAssocSE = (TList*) listSkim2024TimeAssoc -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA");
    THnSparseD *histSparseSkim2024TimeAssoc = (THnSparseD*) listSkim2024TimeAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseSkim2024TimeAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseSkim2024TimeAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntSkim2024TimeAssoc = (TH1D*) histSparseSkim2024TimeAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");

    double countsSkim2024TimeAssoc = histProjIntSkim2024TimeAssoc -> Integral();
    double duplCorrZorroFactorSkim2024TimeAssoc = nEvtsAfterCutsSkim2024TimeAssoc / nEvtsSkim2024TimeAssoc;
    double duplCorrDimuFactorSkim2024TimeAssoc = countsSkimStdAssoc / countsSkim2024TimeAssoc;

    histProjIntSkim2024TimeAssoc -> Scale((duplCorrZorroFactorSkim2024TimeAssoc * evSelCorrSkim2024TimeAssoc) / (bcSelEffSkim2024TimeAssoc * lumiSkim2024TimeAssoc));
    //histProjIntSkim2024TimeAssoc -> Scale(1. / (bcSelEffSkim2024TimeAssoc * lumiSkim2024TimeAssoc));
    SetHist(histProjIntSkim2024TimeAssoc, kOrange+7, 20, 0.8, kOrange+7);

    std::cout << "---------------------------" << std::endl;
    std::cout << "[Skimmed 2024, time assoc.]" << std::endl;
    std::cout << "Duplication [zorro] = " << duplCorrZorroFactorSkim2024TimeAssoc << std::endl;
    std::cout << "Duplication [dimu]  = " << duplCorrDimuFactorSkim2024TimeAssoc << std::endl;
    std::cout << "N. evts.            = " << bcSelEffSkim2024TimeAssoc * nEvtsSkim2024TimeAssoc<< std::endl;
    std::cout << "Luminosity          = " << bcSelEffSkim2024TimeAssoc * lumiSkim2024TimeAssoc << std::endl;
    std::cout << "N. evts. corr.      = " << (bcSelEffSkim2024TimeAssoc * nEvtsSkim2024TimeAssoc) / (duplCorrZorroFactorSkim2024TimeAssoc) << std::endl;
    std::cout << "Luminosity corr.    = " << (bcSelEffSkim2024TimeAssoc * lumiSkim2024TimeAssoc) / (duplCorrZorroFactorSkim2024TimeAssoc) << std::endl;
    std::cout << "---------------------------" << std::endl;


    // Skimmed 2023, std. assoc.
    TFile *fInLumiSkim2023StdAssoc = TFile::Open(fInLumiNameSkim2023StdAssoc.c_str());
    TH1D *hist2023LumiSkim2023StdAssoc = (TH1D*) fInLumiSkim2023StdAssoc -> Get("histLumiSummary");
    double bcSelEffSkim2023StdAssoc = hist2023LumiSkim2023StdAssoc -> GetBinContent(hist2023LumiSkim2023StdAssoc -> GetXaxis() -> FindBin("bcCounterEfficiency"));
    double nEvtsSkim2023StdAssoc = hist2023LumiSkim2023StdAssoc -> GetBinContent(hist2023LumiSkim2023StdAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double lumiSkim2023StdAssoc = hist2023LumiSkim2023StdAssoc -> GetBinContent(hist2023LumiSkim2023StdAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));

    TFile *fInSkim2023StdAssoc = TFile::Open(fInNameSkim2023StdAssoc.c_str());
    // Events histograms
    TList *listSkim2023StdAssocEvSel = (TList*) fInSkim2023StdAssoc -> Get("analysis-event-selection/output");
    TList *listSkim2023StdAssocEvBefCuts = (TList*) listSkim2023StdAssocEvSel -> FindObject("Event_BeforeCuts");
    TList *listSkim2023StdAssocEvAftCuts = (TList*) listSkim2023StdAssocEvSel -> FindObject("Event_AfterCuts");
    TH1D *histSkim2023StdAssocEvBefCuts = (TH1D*) listSkim2023StdAssocEvBefCuts -> FindObject("VtxZ");
    TH1D *histSkim2023StdAssocEvAftCuts = (TH1D*) listSkim2023StdAssocEvAftCuts -> FindObject("VtxZ");
    double evSelCorrSkim2023StdAssoc = (double) histSkim2023StdAssocEvBefCuts -> GetEntries() / (double) histSkim2023StdAssocEvAftCuts -> GetEntries();
    // Dimuon histograms
    TList *listSkim2023StdAssoc = (TList*) fInSkim2023StdAssoc -> Get("analysis-same-event-pairing/output");
    TList *listSkim2023StdAssocSE = (TList*) listSkim2023StdAssoc -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA");
    THnSparseD *histSparseSkim2023StdAssoc = (THnSparseD*) listSkim2023StdAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseSkim2023StdAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseSkim2023StdAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntSkim2023StdAssoc = (TH1D*) histSparseSkim2023StdAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");
    histProjIntSkim2023StdAssoc -> Scale((evSelCorrSkim2023StdAssoc) / (bcSelEffSkim2023StdAssoc * lumiSkim2023StdAssoc));
    SetHist(histProjIntSkim2023StdAssoc, kGreen+4, 20, 0.8, kGreen+4);

    std::cout << "---------------------------" << std::endl;
    std::cout << "[Skimmed 2023, std. assoc.]" << std::endl;
    std::cout << "N. evts.         = " << nEvtsSkim2023StdAssoc * bcSelEffSkim2023StdAssoc << std::endl;
    std::cout << "Luminosity       = " << lumiSkim2023StdAssoc * bcSelEffSkim2023StdAssoc << std::endl;
    std::cout << "N. evts. corr.   = " << nEvtsSkim2023StdAssoc * bcSelEffSkim2023StdAssoc << std::endl;
    std::cout << "Luminosity corr. = " << lumiSkim2023StdAssoc * bcSelEffSkim2023StdAssoc << std::endl;
    std::cout << "---------------------------" << std::endl;

    // Skimmed 2023, time assoc.
    // WARNING! for the time association you have to take into account the duplication using zorro
    TFile *fInLumiSkim2023TimeAssoc = TFile::Open(fInLumiNameSkim2023TimeAssoc.c_str());
    TH1D *histLumiSkim2023TimeAssoc = (TH1D*) fInLumiSkim2023TimeAssoc -> Get("histLumiSummary");
    double bcSelEffSkim2023TimeAssoc = histLumiSkim2023TimeAssoc -> GetBinContent(histLumiSkim2023TimeAssoc -> GetXaxis() -> FindBin("bcCounterEfficiency"));
    double nEvtsSkim2023TimeAssoc = histLumiSkim2023TimeAssoc -> GetBinContent(histLumiSkim2023TimeAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double nEvtsAfterCutsSkim2023TimeAssoc = histLumiSkim2023TimeAssoc -> GetBinContent(histLumiSkim2023TimeAssoc -> GetXaxis() -> FindBin("nEvtsBcSelAfterCuts"));
    double lumiSkim2023TimeAssoc = histLumiSkim2023TimeAssoc -> GetBinContent(histLumiSkim2023TimeAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));

    TFile *fInSkim2023TimeAssoc = TFile::Open(fInNameSkim2023TimeAssoc.c_str());
    // Events histograms
    TList *listSkim2023TimeAssocEvSel = (TList*) fInSkim2023TimeAssoc -> Get("analysis-event-selection/output");
    TList *listSkim2023TimeAssocEvBefCuts = (TList*) listSkim2023TimeAssocEvSel -> FindObject("Event_BeforeCuts");
    TList *listSkim2023TimeAssocEvAftCuts = (TList*) listSkim2023TimeAssocEvSel -> FindObject("Event_AfterCuts");
    TH1D *histSkim2023TimeAssocEvBefCuts = (TH1D*) listSkim2023TimeAssocEvBefCuts -> FindObject("VtxZ");
    TH1D *histSkim2023TimeAssocEvAftCuts = (TH1D*) listSkim2023TimeAssocEvAftCuts -> FindObject("VtxZ");
    double evSelCorrSkim2023TimeAssoc = (double) histSkim2023TimeAssocEvBefCuts -> GetEntries() / (double) histSkim2023TimeAssocEvAftCuts -> GetEntries();
    // Dimuon histograms
    TList *listSkim2023TimeAssoc = (TList*) fInSkim2023TimeAssoc -> Get("analysis-same-event-pairing/output");
    TList *listSkim2023TimeAssocSE = (TList*) listSkim2023TimeAssoc -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA");
    THnSparseD *histSparseSkim2023TimeAssoc = (THnSparseD*) listSkim2023TimeAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseSkim2023TimeAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseSkim2023TimeAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntSkim2023TimeAssoc = (TH1D*) histSparseSkim2023TimeAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");

    double countsSkim2023TimeAssoc = histProjIntSkim2023TimeAssoc -> Integral();
    double duplCorrZorroFactorSkim2023TimeAssoc = nEvtsAfterCutsSkim2023TimeAssoc / nEvtsSkim2023TimeAssoc;
    double duplCorrDimuFactorSkim2023TimeAssoc = countsSkimStdAssoc / countsSkim2023TimeAssoc;

    histProjIntSkim2023TimeAssoc -> Scale((duplCorrZorroFactorSkim2023TimeAssoc * evSelCorrSkim2023TimeAssoc) / (bcSelEffSkim2023TimeAssoc * lumiSkim2023TimeAssoc));
    //histProjIntSkim2023TimeAssoc -> Scale(1. / (bcSelEffSkim2023TimeAssoc * lumiSkim2023TimeAssoc));
    SetHist(histProjIntSkim2023TimeAssoc, kGreen+2, 20, 0.8, kGreen+2);

    std::cout << "---------------------------" << std::endl;
    std::cout << "[Skimmed 2023, time assoc.]" << std::endl;
    std::cout << "Duplication [zorro] = " << duplCorrZorroFactorSkim2023TimeAssoc << std::endl;
    std::cout << "Duplication [dimu]  = " << duplCorrDimuFactorSkim2023TimeAssoc << std::endl;
    std::cout << "N. evts.            = " << bcSelEffSkim2023TimeAssoc * nEvtsSkim2023TimeAssoc<< std::endl;
    std::cout << "Luminosity          = " << bcSelEffSkim2023TimeAssoc * lumiSkim2023TimeAssoc << std::endl;
    std::cout << "N. evts. corr.      = " << (bcSelEffSkim2023TimeAssoc * nEvtsSkim2023TimeAssoc) / (duplCorrZorroFactorSkim2023TimeAssoc * evSelCorrSkim2023TimeAssoc) << std::endl;
    std::cout << "Luminosity corr.    = " << (bcSelEffSkim2023TimeAssoc * lumiSkim2023TimeAssoc) / (duplCorrZorroFactorSkim2023TimeAssoc * evSelCorrSkim2023TimeAssoc) << std::endl;
    std::cout << "---------------------------" << std::endl;

    // Min. Bias 2024, std. assoc.
    // WARNING! event selection already computed, bcSelCutEfficiecny not taken into account (already included in the event selection)
    // this correction is necessary in the trigger for inspectedTVX, for minBias it should not be necessary
    //TFile *fInLumiMinBias2024StdAssoc = TFile::Open("data/2024/LHC24_minBias_std_assoc_luminosity.root");
    TFile *fInLumiMinBias2024StdAssoc = TFile::Open("data/2024/LHC24_minBiasEventStandardSel8_std_assoc_luminosity.root");
    TH1D *histLumiMinBias2024StdAssoc = (TH1D*) fInLumiMinBias2024StdAssoc -> Get("histLumiSummary");
    double nEvtsMinBias2024StdAssoc = histLumiMinBias2024StdAssoc -> GetBinContent(histLumiMinBias2024StdAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double lumiMinBias2024StdAssoc = histLumiMinBias2024StdAssoc -> GetBinContent(histLumiMinBias2024StdAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));

    //TFile *fInMinBias2024StdAssoc = TFile::Open("data/2024/LHC24_minBias_std_assoc.root");
    TFile *fInMinBias2024StdAssoc = TFile::Open("data/2024/LHC24_minBiasEventStandardSel8_std_assoc.root");
    TList *listMinBias2024StdAssoc = (TList*) fInMinBias2024StdAssoc -> Get("analysis-same-event-pairing/output");
    TList *listMinBias2024StdAssocSE = (TList*) listMinBias2024StdAssoc -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA");
    THnSparseD *histSparseMinBias2024StdAssoc = (THnSparseD*) listMinBias2024StdAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseMinBias2024StdAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseMinBias2024StdAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntMinBias2024StdAssoc = (TH1D*) histSparseMinBias2024StdAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");
    double countsMinBias2024StdAssoc = histProjIntMinBias2024StdAssoc -> Integral();
    histProjIntMinBias2024StdAssoc -> Scale(1. / (lumiMinBias2024StdAssoc));
    SetHist(histProjIntMinBias2024StdAssoc, kBlue+1, 24, 0.8, kBlue+1);

    // Added for J/psi - D0 associated production
    TH1D *histLuminosityMinBias2024StdAssoc = new TH1D("histLumi", "; ; Luminosity (pb-1)", 1, 0, 1);
    histLuminosityMinBias2024StdAssoc -> SetBinContent(1, lumiMinBias2024StdAssoc);

    TFile *fOutLumiMinBias2024StdAssoc = new TFile("data/luminosity_jpsi_LHC24_minBias.root", "RECREATE");
    histLuminosityMinBias2024StdAssoc -> Write();
    fOutLumiMinBias2024StdAssoc -> Close();

    std::cout << "---------------------------" << std::endl;
    std::cout << "[Min. Bias, std. assoc.]" << std::endl;
    std::cout << "N. evts.         = " << nEvtsMinBias2024StdAssoc << std::endl;
    std::cout << "Luminosity       = " << lumiMinBias2024StdAssoc << std::endl;
    std::cout << "N. evts. corr.   = " << nEvtsMinBias2024StdAssoc << std::endl;
    std::cout << "Luminosity corr. = " << lumiMinBias2024StdAssoc << std::endl;
    std::cout << "---------------------------" << std::endl;

    // Min. Bias 2024, time assoc.
    TFile *fInLumiMinBias2024TimeAssoc = TFile::Open("data/2024/LHC24_minBias_time_assoc_luminosity.root");
    TH1D *histLumiMinBias2024TimeAssoc = (TH1D*) fInLumiMinBias2024TimeAssoc -> Get("histLumiSummary");
    double nEvtsMinBias2024TimeAssoc = histLumiMinBias2024TimeAssoc -> GetBinContent(histLumiMinBias2024TimeAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double lumiMinBias2024TimeAssoc = histLumiMinBias2024TimeAssoc -> GetBinContent(histLumiMinBias2024TimeAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));
    TFile *fInMinBias2024TimeAssoc = TFile::Open("data/2024/LHC24_minBias_time_assoc.root");
    TList *listMinBias2024TimeAssoc = (TList*) fInMinBias2024TimeAssoc -> Get("analysis-same-event-pairing/output");
    TList *listMinBias2024TimeAssocSE = (TList*) listMinBias2024TimeAssoc -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA");
    THnSparseD *histSparseMinBias2024TimeAssoc = (THnSparseD*) listMinBias2024TimeAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseMinBias2024TimeAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseMinBias2024TimeAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntMinBias2024TimeAssoc = (TH1D*) histSparseMinBias2024TimeAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");

    double countsMinBias2024TimeAssoc = histProjIntMinBias2024TimeAssoc -> Integral();
    double duplCorrFactorMinBias2024TimeAssoc = (countsMinBias2024StdAssoc / nEvtsMinBias2024StdAssoc) / (countsMinBias2024TimeAssoc / nEvtsMinBias2024TimeAssoc);
    
    histProjIntMinBias2024TimeAssoc -> Scale(duplCorrFactorMinBias2024TimeAssoc / (lumiMinBias2024TimeAssoc));
    //histProjIntMinBias2024TimeAssoc -> Scale(1. / (lumiMinBias2024TimeAssoc));
    SetHist(histProjIntMinBias2024TimeAssoc, kAzure+2, 24, 0.8, kAzure+2);

    std::cout << "---------------------------" << std::endl;
    std::cout << "[Min. Bias, time assoc.]" << std::endl;
    std::cout << "Duplication      = " << duplCorrFactorMinBias2024TimeAssoc << std::endl;
    std::cout << "N. evts.         = " << nEvtsMinBias2024TimeAssoc << std::endl;
    std::cout << "Luminosity       = " << lumiMinBias2024TimeAssoc << std::endl;
    std::cout << "N. evts. corr.   = " << nEvtsMinBias2024TimeAssoc / duplCorrFactorMinBias2024TimeAssoc << std::endl;
    std::cout << "Luminosity corr. = " << lumiMinBias2024TimeAssoc / duplCorrFactorMinBias2024TimeAssoc << std::endl;
    std::cout << "---------------------------" << std::endl;


    //std::cout << "Min Bias dimuons time assoc. = " << countsMinBias2024TimeAssoc << " ; coll. after cuts = " << collsAfterCutsMinBiasTimeAssoc << std::endl;
    //std::cout << "Min Bias dimuons std. assoc. = " << countsMinBias2024StdAssoc << " ; coll. after cuts = " << collsAfterCutsMinBiasStdAssoc << std::endl;
    //std::cout << "Correction factor from inv. mass = " << (countsMinBias2024TimeAssoc / collsAfterCutsMinBiasTimeAssoc) / (countsMinBias2024StdAssoc / collsAfterCutsMinBiasStdAssoc) << std::endl;
    //std::cout << "Correction factor from zorro     = " << lumiSkim2024TimeAssoc / lumiSkim2024StdAssoc << std::endl;

    TCanvas *canvasCompMassSpectra = new TCanvas("canvasCompMassSpectra", "", 800, 600);
    gPad -> SetLogy(true);
    histProjIntSkim2024StdAssoc -> SetTitle("");
    histProjIntSkim2024StdAssoc -> GetXaxis() -> SetRangeUser(2, 5);
    histProjIntSkim2024StdAssoc -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c})");
    histProjIntSkim2024StdAssoc -> GetYaxis() -> SetRangeUser(1e1, 1e5);
    histProjIntSkim2024StdAssoc -> GetYaxis() -> SetTitle("(1 / #it{L}) * d#it{N}/d#it{M}_{#mu#mu} (GeV/#it{c}^{2} pb^{-1})^{-1}");
    histProjIntSkim2024StdAssoc -> Draw("EP");
    histProjIntSkim2024TimeAssoc -> Draw("EP SAME");
    //histProjIntSkim2023StdAssoc -> Draw("EP SAME");
    //histProjIntSkim2023TimeAssoc -> Draw("EP SAME");
    histProjIntMinBias2024StdAssoc -> Draw("EP SAME");
    histProjIntMinBias2024TimeAssoc -> Draw("EP SAME");

    TLegend *legendCompMassSpectra = new TLegend(0.16, 0.20, 0.36, 0.45, " ", "brNDC");
    SetLegend(legendCompMassSpectra);
    //legendCompMassSpectra -> AddEntry(histProjIntSkim2023StdAssoc, Form("LHC23 skimmed std , #it{L} = %4.3f pb^{-1}", bcSelEffSkim2023StdAssoc * lumiSkim2023StdAssoc), "EP");
    //legendCompMassSpectra -> AddEntry(histProjIntSkim2023TimeAssoc, Form("LHC23 skimmed time, #it{L} = %4.3f pb^{-1}", bcSelEffSkim2023TimeAssoc * lumiSkim2023TimeAssoc), "EP");
    legendCompMassSpectra -> AddEntry(histProjIntSkim2024StdAssoc, Form("LHC24af skimmed std , #it{L} = %4.3f pb^{-1}", bcSelEffSkim2024StdAssoc * lumiSkim2024StdAssoc), "EP");
    legendCompMassSpectra -> AddEntry(histProjIntSkim2024TimeAssoc, Form("LHC24af skimmed time, #it{L} = %4.3f pb^{-1}", bcSelEffSkim2024TimeAssoc * lumiSkim2024TimeAssoc), "EP");
    legendCompMassSpectra -> AddEntry(histProjIntMinBias2024StdAssoc, Form("LHC24 MB std , #it{L} = %4.3f pb^{-1}", lumiMinBias2024StdAssoc), "EP");
    legendCompMassSpectra -> AddEntry(histProjIntMinBias2024TimeAssoc, Form("LHC24 MB time, #it{L} = %4.3f pb^{-1}", lumiMinBias2024TimeAssoc), "EP");
    legendCompMassSpectra -> Draw();

    // Compute the ratio
    TH1D *histRatioIntSkimStdAssoc2024SkimTimeAssoc2024 = (TH1D*) histProjIntSkim2024StdAssoc -> Clone("histRatioIntSkimStdAssoc2024SkimTimeAssoc2024");
    histRatioIntSkimStdAssoc2024SkimTimeAssoc2024 -> Divide(histProjIntSkim2024TimeAssoc);
    SetHist(histRatioIntSkimStdAssoc2024SkimTimeAssoc2024, kRed+1, 20, 0.8, kRed+1);

    TH1D *histRatioIntMinBiasStdAssoc2024MinBiasTimeAssoc2024 = (TH1D*) histProjIntMinBias2024StdAssoc -> Clone("histRatioIntMinBiasStdAssoc2024MinBiasTimeAssoc2024");
    histRatioIntMinBiasStdAssoc2024MinBiasTimeAssoc2024 -> Divide(histProjIntMinBias2024TimeAssoc);
    SetHist(histRatioIntMinBiasStdAssoc2024MinBiasTimeAssoc2024, kBlue+1, 24, 0.8, kBlue+1);

    TH1D *histRatioIntSkimStdAssoc2024MinBiasStdAssoc2024 = (TH1D*) histProjIntSkim2024StdAssoc -> Clone("histRatioIntSkimStdAssoc2024MinBiasStdAssoc2024");
    histRatioIntSkimStdAssoc2024MinBiasStdAssoc2024 -> Divide(histProjIntMinBias2024StdAssoc);
    SetHist(histRatioIntSkimStdAssoc2024MinBiasStdAssoc2024, kRed+1, 24, 0.8, kRed+1);

    TH1D *histRatioIntSkimTimeAssoc2024MinBiasTimeAssoc2024 = (TH1D*) histProjIntSkim2024TimeAssoc -> Clone("histRatioIntSkimTimeAssoc2024MinBiasTimeAssoc2024");
    histRatioIntSkimTimeAssoc2024MinBiasTimeAssoc2024 -> Divide(histProjIntMinBias2024TimeAssoc);
    SetHist(histRatioIntSkimTimeAssoc2024MinBiasTimeAssoc2024, kOrange+7, 24, 0.8, kOrange+7);

    TH1D *histRatioIntSkimStdAssoc2023SkimTimeAssoc2023 = (TH1D*) histProjIntSkim2023StdAssoc -> Clone("histRatioIntSkimStdAssoc2023SkimTimeAssoc2023");
    histRatioIntSkimStdAssoc2023SkimTimeAssoc2023 -> Divide(histProjIntSkim2023TimeAssoc);
    SetHist(histRatioIntSkimStdAssoc2023SkimTimeAssoc2023, kGreen+2, 20, 0.8, kGreen+2);


    TH1D *histRatioIntSkimStdAssoc2023SkimStdAssoc2024 = (TH1D*) histProjIntSkim2023StdAssoc -> Clone("histRatioIntSkimStdAssoc2023SkimStdAssoc2024");
    histRatioIntSkimStdAssoc2023SkimStdAssoc2024 -> Divide(histProjIntSkim2024StdAssoc);
    SetHist(histRatioIntSkimStdAssoc2023SkimStdAssoc2024, kYellow+2, 20, 0.8, kYellow+2);

    TH1D *histRatioIntSkimTimeAssoc2023SkimTimeAssoc2024 = (TH1D*) histProjIntSkim2023TimeAssoc -> Clone("histRatioIntSkimTimeAssoc2023SkimTimeAssoc2024");
    histRatioIntSkimTimeAssoc2023SkimTimeAssoc2024 -> Divide(histProjIntSkim2024TimeAssoc);
    SetHist(histRatioIntSkimTimeAssoc2023SkimTimeAssoc2024, kYellow+4, 20, 0.8, kYellow+4);

    TCanvas *canvasRatioMassSpectra = new TCanvas("canvasRatioMassSpectra", "", 800, 600);
    histRatioIntSkimStdAssoc2024SkimTimeAssoc2024 -> GetXaxis() -> SetRangeUser(2, 5);
    histRatioIntSkimStdAssoc2024SkimTimeAssoc2024 -> GetYaxis() -> SetRangeUser(0.3, 1.7);
    histRatioIntSkimStdAssoc2024SkimTimeAssoc2024 -> GetYaxis() -> SetTitle("Ratio");
    histRatioIntSkimStdAssoc2024SkimTimeAssoc2024 -> Draw("EP");
    histRatioIntSkimStdAssoc2024MinBiasStdAssoc2024 -> Draw("EP SAME");
    histRatioIntSkimTimeAssoc2024MinBiasTimeAssoc2024 -> Draw("EP SAME");
    histRatioIntMinBiasStdAssoc2024MinBiasTimeAssoc2024 -> Draw("EP SAME");

    TLine *lineUnity = new TLine(2, 1, 5, 1);
    lineUnity -> SetLineColor(kGray+2);
    lineUnity -> SetLineStyle(kDashed);
    lineUnity -> Draw();

    TLegend *legendRatioMassSpectra = new TLegend(0.20, 0.20, 0.60, 0.40, " ", "brNDC");
    SetLegend(legendRatioMassSpectra);
    legendRatioMassSpectra -> AddEntry(histRatioIntSkimStdAssoc2024SkimTimeAssoc2024, "skimmed std. / skimmed time", "EP");
    legendRatioMassSpectra -> AddEntry(histRatioIntSkimStdAssoc2024MinBiasStdAssoc2024, "skimmed std. / Min. Bias std.", "EP");
    legendRatioMassSpectra -> AddEntry(histRatioIntSkimTimeAssoc2024MinBiasTimeAssoc2024, "skimmed time / Min. Bias time", "EP");
    legendRatioMassSpectra -> AddEntry(histRatioIntMinBiasStdAssoc2024MinBiasTimeAssoc2024, "Min. Bias std. / Min. Bias time", "EP");
    legendRatioMassSpectra -> Draw();

    TLatex latexTitle;
    latexTitle.SetNDC();
    latexTitle.SetTextSize(0.05);
    latexTitle.SetTextFont(42);

    TCanvas *canvasCompMassSpectraSummary = new TCanvas("canvasCompMassSpectraSummary", "", 1200, 1200);
    canvasCompMassSpectraSummary -> Divide(2, 2);

    canvasCompMassSpectraSummary -> cd(1);
    gPad -> SetLogy(true);
    histProjIntSkim2023StdAssoc -> SetTitle("");
    histProjIntSkim2023StdAssoc -> GetXaxis() -> SetRangeUser(2, 5);
    histProjIntSkim2023StdAssoc -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c})");
    histProjIntSkim2023StdAssoc -> GetYaxis() -> SetRangeUser(1e2, 1e5);
    histProjIntSkim2023StdAssoc -> GetYaxis() -> SetTitle("(1 / #it{L}) * d#it{N}/d#it{M}_{#mu#mu} (GeV/#it{c}^{2} pb^{-1})^{-1}");
    histProjIntSkim2023StdAssoc -> Draw("EP");
    histProjIntSkim2023TimeAssoc -> Draw("EP SAME");
    latexTitle.DrawLatex(0.70, 0.85, "pp 2023 pass4");

    canvasCompMassSpectraSummary -> cd(2);
    gPad -> SetLogy(true);
    histProjIntSkim2024StdAssoc -> SetTitle("");
    histProjIntSkim2024StdAssoc -> GetXaxis() -> SetRangeUser(2, 5);
    histProjIntSkim2024StdAssoc -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c})");
    histProjIntSkim2024StdAssoc -> GetYaxis() -> SetRangeUser(1e2, 1e5);
    histProjIntSkim2024StdAssoc -> GetYaxis() -> SetTitle("(1 / #it{L}) * d#it{N}/d#it{M}_{#mu#mu} (GeV/#it{c}^{2} pb^{-1})^{-1}");
    histProjIntSkim2024StdAssoc -> Draw("EP");
    histProjIntSkim2024TimeAssoc -> Draw("EP SAME");
    latexTitle.DrawLatex(0.70, 0.85, "pp 2024 pass1");
    
    canvasCompMassSpectraSummary -> cd(3);
    histRatioIntSkimStdAssoc2023SkimTimeAssoc2023 -> SetTitle("");
    histRatioIntSkimStdAssoc2023SkimTimeAssoc2023 -> GetXaxis() -> SetRangeUser(2, 5);
    histRatioIntSkimStdAssoc2023SkimTimeAssoc2023 -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c})");
    histRatioIntSkimStdAssoc2023SkimTimeAssoc2023 -> GetYaxis() -> SetRangeUser(0.5, 1.5);
    histRatioIntSkimStdAssoc2023SkimTimeAssoc2023 -> GetYaxis() -> SetTitle("Ratio");
    histRatioIntSkimStdAssoc2023SkimTimeAssoc2023 -> Draw("EP");
    histRatioIntSkimStdAssoc2023SkimStdAssoc2024 -> Draw("EP SAME");
    histRatioIntSkimTimeAssoc2023SkimTimeAssoc2024 -> Draw("EP SAME");
    latexTitle.DrawLatex(0.70, 0.85, "pp 2023 pass4");
    lineUnity -> Draw();


    canvasCompMassSpectraSummary -> cd(4);
    histRatioIntSkimStdAssoc2024SkimTimeAssoc2024 -> SetTitle("");
    histRatioIntSkimStdAssoc2024SkimTimeAssoc2024 -> GetXaxis() -> SetRangeUser(2, 5);
    histRatioIntSkimStdAssoc2024SkimTimeAssoc2024 -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c})");
    histRatioIntSkimStdAssoc2024SkimTimeAssoc2024 -> GetYaxis() -> SetRangeUser(0.5, 1.5);
    histRatioIntSkimStdAssoc2024SkimTimeAssoc2024 -> GetYaxis() -> SetTitle("Ratio");
    histRatioIntSkimStdAssoc2024SkimTimeAssoc2024 -> Draw("EP");
    latexTitle.DrawLatex(0.70, 0.85, "pp 2024 pass1");
    lineUnity -> Draw();

    canvasCompMassSpectra -> SaveAs("figures/CompMassSpectra.pdf");
    canvasRatioMassSpectra -> SaveAs("figures/RatioMassSpectra.pdf");
    canvasCompMassSpectraSummary -> SaveAs("figures/CompMassSpectraSummary.pdf");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.04);
}