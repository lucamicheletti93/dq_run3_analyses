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

void get_normalization_from_grid(string year = "2024", string period = "LHC24", string subPeriod = "None", string triggerMask = "fDiMuon", string assocType = "time_assoc") {
    ifstream fInAlienInputPath(Form("run_lists/%s/%s/%s/alien_input_path.txt", year.c_str(), triggerMask.c_str(), assocType.c_str()));
    if (!fInAlienInputPath.is_open()) {
        std::cout << "Error opening the file!" << std::endl;
        return;
    }

    ifstream fInAlienRunList(Form("run_lists/%s/%s/%s/alien_run_list.txt", year.c_str(), triggerMask.c_str(), assocType.c_str()));
    if (!fInAlienRunList.is_open()) {
        std::cout << "Error opening the file!" << std::endl;
        return;
    }

    string subRunListName;
    if (subPeriod == "None") {
        subRunListName = Form("run_lists/%s/%s/%s/alien_run_list.txt", year.c_str(), triggerMask.c_str(), assocType.c_str());
    } else {
        subRunListName = Form("run_lists/%s/%s/%s.txt", year.c_str(), triggerMask.c_str(), subPeriod.c_str());
    }

    ifstream fInAlienSubRunList(subRunListName);
    if (!fInAlienSubRunList.is_open()) {
        std::cout << "Error opening the file!" << std::endl;
        return;
    }

    double counters[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    std::vector <double> vecInfoCounterTVX;
    std::vector <double> vecInfoCounterScalTrig;
    std::vector <double> vecInfoCounterSelTrig;
    std::vector <double> vecSelCounterSelTOI;
    std::vector <double> vecSelCounterAnalysedTrig;
    std::vector <string> vecRunList;
    std::vector <string> vecSubRunList;

    string alienRunListName;
    while (getline(fInAlienRunList, alienRunListName)) {
        vecRunList.push_back(alienRunListName);
    }

    string alienSubRunListName;
    while (getline(fInAlienSubRunList, alienSubRunListName)) {
        vecSubRunList.push_back(alienSubRunListName);
    }

    TH1D *histZorroInfoCounterTVX = new TH1D("histZorroInfoCounterTVX", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histZorroInfoCounterScalTrig = new TH1D("histZorroInfoCounterScalTrig", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histZorroInfoCounterSelTrig = new TH1D("histZorroInfoCounterSelTrig", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histZorroSelCounterSelTOI = new TH1D("histZorroSelCounterSelTOI", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histZorroSelCounterAnalysedTrig = new TH1D("histZorroSelCounterAnalysedTrig", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histEvSelCollisionBeforeFiltering = new TH1D("histEvSelCollisionBeforeFiltering", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histEvSelCollisionBeforeCuts = new TH1D("histEvSelCollisionBeforeCuts", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histEvSelCollisionAfterCuts = new TH1D("histEvSelCollisionAfterCuts", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histBcSelCounterTVX = new TH1D("histBcSelCounterTVX", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histBcSelCounterTVXafterBCcuts = new TH1D("histBcSelCounterTVXafterBCcuts", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histEvSelColCounterTVX = new TH1D("histEvSelColCounterTVX", "", vecSubRunList.size(), 0, vecSubRunList.size());

    string alienPathName;
    string runNumber;
    int runCounter = 0;
    int subRunCounter = 0;
    while (getline(fInAlienInputPath, alienPathName)) {
        runNumber = vecRunList.at(runCounter).c_str();

        std::cout << "Processing Run " << runNumber << " ..." << std::endl;
        if (std::find(vecSubRunList.begin(), vecSubRunList.end(), runNumber) == vecSubRunList.end()) {
            runCounter++;
            continue;
        }

        RetrieveTriggerInfo(alienPathName.c_str(), true, triggerMask, counters);
        //RetrieveTriggerInfo(Form("%s/AOD", alienPathName.c_str()), true, triggerMask, counters);

        histZorroInfoCounterTVX -> SetBinContent(subRunCounter+1, counters[0]);
        histZorroInfoCounterTVX -> GetXaxis() -> SetBinLabel(subRunCounter+1, runNumber.c_str());

        histZorroInfoCounterScalTrig -> SetBinContent(subRunCounter+1, counters[1]);
        histZorroInfoCounterScalTrig -> GetXaxis() -> SetBinLabel(subRunCounter+1, runNumber.c_str());

        histZorroInfoCounterSelTrig -> SetBinContent(subRunCounter+1, counters[2]);
        histZorroInfoCounterSelTrig -> GetXaxis() -> SetBinLabel(subRunCounter+1, runNumber.c_str());

        histZorroSelCounterSelTOI -> SetBinContent(subRunCounter+1, counters[3]);
        histZorroSelCounterSelTOI -> GetXaxis() -> SetBinLabel(subRunCounter+1, runNumber.c_str());

        histZorroSelCounterAnalysedTrig -> SetBinContent(subRunCounter+1, counters[4]);
        histZorroSelCounterAnalysedTrig -> GetXaxis() -> SetBinLabel(subRunCounter+1, runNumber.c_str());

        histEvSelCollisionBeforeFiltering -> SetBinContent(subRunCounter+1, counters[5]);
        histEvSelCollisionBeforeFiltering -> GetXaxis() -> SetBinLabel(subRunCounter+1, runNumber.c_str());

        histEvSelCollisionBeforeCuts -> SetBinContent(subRunCounter+1, counters[6]);
        histEvSelCollisionBeforeCuts -> GetXaxis() -> SetBinLabel(subRunCounter+1, runNumber.c_str());

        histEvSelCollisionAfterCuts -> SetBinContent(subRunCounter+1, counters[7]);
        histEvSelCollisionAfterCuts -> GetXaxis() -> SetBinLabel(subRunCounter+1, runNumber.c_str());

        histBcSelCounterTVX -> SetBinContent(subRunCounter+1, counters[8]);
        histBcSelCounterTVX -> GetXaxis() -> SetBinLabel(subRunCounter+1, runNumber.c_str());

        histBcSelCounterTVXafterBCcuts -> SetBinContent(subRunCounter+1, counters[9]);
        histBcSelCounterTVXafterBCcuts -> GetXaxis() -> SetBinLabel(subRunCounter+1, runNumber.c_str());

        histEvSelColCounterTVX -> SetBinContent(subRunCounter+1, counters[10]);
        histEvSelColCounterTVX -> GetXaxis() -> SetBinLabel(subRunCounter+1, runNumber.c_str());

        subRunCounter++;
        runCounter++;
    }

    string fOutName;
    if (subPeriod == "None") {
        fOutName = Form("data/%s/%s_%s_%s_trigger_summary.root", year.c_str(), period.c_str(), triggerMask.c_str(), assocType.c_str());
    } else {
        fOutName = Form("data/%s/%s_%s_%s_trigger_summary.root", year.c_str(), subPeriod.c_str(), triggerMask.c_str(), assocType.c_str());
    }

    TFile *fOut = new TFile(fOutName.c_str(), "RECREATE");
    histZorroInfoCounterTVX -> Write();
    histZorroInfoCounterScalTrig -> Write();
    histZorroInfoCounterSelTrig -> Write();
    histZorroSelCounterSelTOI -> Write();
    histZorroSelCounterAnalysedTrig -> Write();
    histEvSelCollisionBeforeFiltering -> Write();
    histEvSelCollisionBeforeCuts -> Write();
    histEvSelCollisionAfterCuts -> Write();
    histBcSelCounterTVX -> Write();
    histBcSelCounterTVXafterBCcuts -> Write();
    histEvSelColCounterTVX -> Write();
    fOut -> Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_normalization_from_local(string year = "2024", string period = "LHC24af", string triggerMask = "fDiMuon", string assocType = "std_assoc") {
    string fInName = Form("data/%s/%s_%s_%s.root", year.c_str(), period.c_str(), triggerMask.c_str(), assocType.c_str());
    string fOutName = Form("data/%s/%s_%s_%s_trigger_summary.root", year.c_str(), period.c_str(), triggerMask.c_str(), assocType.c_str());

    TFile *fIn = TFile::Open(fInName.c_str());
    if (!fIn || fIn -> IsZombie()) {
        return;
    }

    TH1D *histZorroInfoCounterTVX = nullptr;
    TH1D *histZorroInfoCounterScalTrig = nullptr;
    TH1D *histZorroInfoCounterSelTrig = nullptr;
    TH1D *histZorroSelCounterSelTOI = nullptr;
    TH1D *histZorroSelCounterAnalysedTrig = nullptr;
    TH1D *histEvSelCollisionBeforeFiltering = nullptr;
    TH1D *histEvSelCollisionBeforeCuts = nullptr;
    TH1D *histEvSelCollisionAfterCuts = nullptr;
    TH1D *histBcSelCounterTVX = nullptr;
    TH1D *histBcSelCounterTVXafterBCcuts = nullptr;
    TH1D *histEvSelColCounterTVX = nullptr;

    for (auto const& dirKey : *fIn -> GetListOfKeys()) {
        if (TString(dirKey -> GetName()).Contains("table-maker")) {
            TList *list = (TList*) fIn -> Get("table-maker/Statistics");
            TH2D *histZorroInfo = (TH2D*) list -> FindObject("ZorroInfo");
            TH2D *histZorroSel = (TH2D*) list -> FindObject("ZorroSel");
            TH2D *histEventStats = (TH2D*) list -> FindObject("EventStats");

            int nRuns = histZorroSel -> GetXaxis() -> GetNbins();
            vector <string> vecRunList;
            for (int iRun = 0;iRun < nRuns;iRun++) {
                string tmpRunNumber = histZorroInfo -> GetXaxis() -> GetBinLabel(iRun+1);
                if (!(tmpRunNumber.empty())) {
                    vecRunList.push_back(histZorroInfo -> GetXaxis() -> GetBinLabel(iRun+1));
                }
            }

            double collisionsBeforeFiltering = 0;
            double collisionsBeforeCuts = 0;
            double collisionsAfterCuts = 0;
            double infoTVX = 0;
            double infoScalTrig = 0;
            double infoSelTrig = 0;
            double selTOI = 0;
            double analysedTriggers = 0;

            histZorroInfoCounterTVX = new TH1D("histZorroInfoCounterTVX", "", vecRunList.size(), 0, vecRunList.size());
            histZorroInfoCounterScalTrig = new TH1D("histZorroInfoCounterScalTrig", "", vecRunList.size(), 0, vecRunList.size());
            histZorroInfoCounterSelTrig = new TH1D("histZorroInfoCounterSelTrig", "", vecRunList.size(), 0, vecRunList.size());
            histZorroSelCounterSelTOI = new TH1D("histZorroSelCounterSelTOI", "", vecRunList.size(), 0, vecRunList.size());
            histZorroSelCounterAnalysedTrig = new TH1D("histZorroSelCounterAnalysedTrig", "", vecRunList.size(), 0, vecRunList.size());
            histEvSelCollisionBeforeFiltering = new TH1D("histEvSelCollisionBeforeFiltering", "", 1, 0, 1);
            histEvSelCollisionBeforeCuts = new TH1D("histEvSelCollisionBeforeCuts", "", 1, 0, 1);
            histEvSelCollisionAfterCuts = new TH1D("histEvSelCollisionAfterCuts", "", 1, 0, 1);

            for (int iRun = 0;iRun < int (vecRunList.size());iRun++) {
                string runNumber = vecRunList.at(iRun).c_str();
                // Possible difference in the bin numbers among the 2 histograms
                int binZorroInfo = histZorroInfo -> GetXaxis() -> FindBin(runNumber.c_str());
                int binZorroSel = histZorroSel -> GetXaxis() -> FindBin(runNumber.c_str());

                infoTVX = histZorroInfo -> GetBinContent(binZorroInfo, histZorroInfo -> GetYaxis() -> FindBin("inspected TVX"));
                infoScalTrig = histZorroInfo -> GetBinContent(binZorroInfo, histZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str())));
                infoSelTrig = histZorroInfo -> GetBinContent(binZorroInfo, histZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str())));
                selTOI = histZorroSel -> GetBinContent(binZorroSel, histZorroSel -> GetYaxis() -> FindBin(Form("%s", triggerMask.c_str())));
                analysedTriggers = histZorroSel -> GetBinContent(binZorroSel, histZorroSel -> GetYaxis() -> FindBin(Form("%s AnalysedTriggers", triggerMask.c_str())));

                histZorroInfoCounterTVX -> SetBinContent(iRun+1, infoTVX);
                histZorroInfoCounterTVX -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());

                histZorroInfoCounterScalTrig -> SetBinContent(iRun+1, infoScalTrig);
                histZorroInfoCounterScalTrig -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());

                histZorroInfoCounterSelTrig -> SetBinContent(iRun+1, infoSelTrig);
                histZorroInfoCounterSelTrig -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());

                histZorroSelCounterSelTOI -> SetBinContent(iRun+1, selTOI);
                histZorroSelCounterSelTOI -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());

                histZorroSelCounterAnalysedTrig -> SetBinContent(iRun+1, analysedTriggers);
                histZorroSelCounterAnalysedTrig -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());
            }

            collisionsBeforeFiltering = histEventStats -> GetBinContent(2, histEventStats -> GetYaxis() -> FindBin("Total"));
            collisionsBeforeCuts = histEventStats -> GetBinContent(3, histEventStats -> GetYaxis() -> FindBin("Total"));
            collisionsAfterCuts = histEventStats -> GetBinContent(4, histEventStats -> GetYaxis() -> FindBin("Total"));

            histEvSelCollisionBeforeFiltering -> SetBinContent(1, collisionsBeforeFiltering);
            histEvSelCollisionBeforeCuts -> SetBinContent(1, collisionsBeforeCuts);
            histEvSelCollisionAfterCuts -> SetBinContent(1, collisionsAfterCuts);
        }
        if (TString(dirKey -> GetName()).Contains("bc-selection-task")) {
            TH1D *histCounterTVX = (TH1D*) fIn -> Get("bc-selection-task/hCounterTVX");
            TH1D *histCounterTVXafterBCcuts = (TH1D*) fIn -> Get("bc-selection-task/hCounterTVXafterBCcuts");

            int nRuns = histCounterTVX -> GetXaxis() -> GetNbins();
            vector <string> vecRunList;
            for (int iRun = 0;iRun < nRuns;iRun++) {
                string tmpRunNumber = histCounterTVX -> GetXaxis() -> GetBinLabel(iRun+1);
                if (!(tmpRunNumber.empty())) {
                    vecRunList.push_back(histCounterTVX -> GetXaxis() -> GetBinLabel(iRun+1));
                }
            }

            histBcSelCounterTVX = new TH1D("histBcSelCounterTVX", "", vecRunList.size(), 0, vecRunList.size());
            histBcSelCounterTVXafterBCcuts = new TH1D("histBcSelCounterTVXafterBCcuts", "", vecRunList.size(), 0, vecRunList.size());

            for (int iRun = 0;iRun < int (vecRunList.size());iRun++) {
                string runNumber = vecRunList.at(iRun).c_str();
                histBcSelCounterTVX -> SetBinContent(iRun+1, histCounterTVX -> GetBinContent(iRun+1));
                histBcSelCounterTVX -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());

                histBcSelCounterTVXafterBCcuts -> SetBinContent(iRun+1, histCounterTVXafterBCcuts -> GetBinContent(iRun+1));
                histBcSelCounterTVXafterBCcuts -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());
            }
        }

        if (TString(dirKey -> GetName()).Contains("event-selection-task")) {
            TH1D *histColCounterTVX = (TH1D*) fIn -> Get("event-selection-task/hColCounterTVX");

            int nRuns = histColCounterTVX -> GetXaxis() -> GetNbins();
            vector <string> vecRunList;
            for (int iRun = 0;iRun < nRuns;iRun++) {
                string tmpRunNumber = histColCounterTVX -> GetXaxis() -> GetBinLabel(iRun+1);
                if (!(tmpRunNumber.empty())) {
                    vecRunList.push_back(histColCounterTVX -> GetXaxis() -> GetBinLabel(iRun+1));
                }
            }

            histEvSelColCounterTVX = new TH1D("histEvSelColCounterTVX", "", vecRunList.size(), 0, vecRunList.size());

            for (int iRun = 0;iRun < int (vecRunList.size());iRun++) {
                string runNumber = vecRunList.at(iRun).c_str();
                histEvSelColCounterTVX -> SetBinContent(iRun+1, histColCounterTVX -> GetBinContent(iRun+1));
                histEvSelColCounterTVX -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());
            }
        }
    }

    TFile *fOut = new TFile(fOutName.c_str(), "RECREATE");
    histZorroInfoCounterTVX -> Write();
    histZorroInfoCounterScalTrig -> Write();
    histZorroInfoCounterSelTrig -> Write();
    histZorroSelCounterSelTOI -> Write();
    histZorroSelCounterAnalysedTrig -> Write();
    histEvSelCollisionBeforeFiltering -> Write();
    histEvSelCollisionBeforeCuts -> Write();
    histEvSelCollisionAfterCuts -> Write();
    histBcSelCounterTVX -> Write();
    histBcSelCounterTVXafterBCcuts -> Write();
    histEvSelColCounterTVX -> Write();
    fOut -> Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void luminosity(string year = "2024", string period = "LHC24af", string triggerMask = "fDiMuon", string assocType = "std_assoc") {
    //TFile *fIn = new TFile(Form("run_lists/%s/%s_%s_%s_trigger_summary.root", year.c_str(), period.c_str(), triggerMask.c_str(), assocType.c_str()), "READ");
    TFile *fIn = new TFile(Form("data/%s/%s_%s_%s_trigger_summary.root", year.c_str(), period.c_str(), triggerMask.c_str(), assocType.c_str()), "READ");
    TH1D *histZorroInfoCounterTVX = (TH1D*) fIn -> Get("histZorroInfoCounterTVX");
    TH1D *histZorroInfoCounterScalTrig = (TH1D*) fIn -> Get("histZorroInfoCounterScalTrig");
    TH1D *histZorroInfoCounterSelTrig = (TH1D*) fIn -> Get("histZorroInfoCounterSelTrig");
    TH1D *histZorroSelCounterSelTOI = (TH1D*) fIn -> Get("histZorroSelCounterSelTOI");
    TH1D *histZorroSelCounterAnalysedTrig = (TH1D*) fIn -> Get("histZorroSelCounterAnalysedTrig");
    TH1D *histEvSelCollisionBeforeFiltering = (TH1D*) fIn -> Get("histEvSelCollisionBeforeFiltering");
    TH1D *histEvSelCollisionBeforeCuts = (TH1D*) fIn -> Get("histEvSelCollisionBeforeCuts");
    TH1D *histEvSelCollisionAfterCuts = (TH1D*) fIn -> Get("histEvSelCollisionAfterCuts");
    TH1D *histBcSelCounterTVX = (TH1D*) fIn -> Get("histBcSelCounterTVX");
    TH1D *histBcSelCounterTVXafterBCcuts = (TH1D*) fIn -> Get("histBcSelCounterTVXafterBCcuts");
    TH1D *histEvSelColCounterTVX = (TH1D*) fIn -> Get("histEvSelColCounterTVX");

    const int nRuns = histZorroInfoCounterTVX -> GetEntries();
    double inspectedTVX = histZorroInfoCounterTVX -> Integral();
    double selectionsZorroInfo = histZorroInfoCounterSelTrig -> Integral();
    double scalers = histZorroInfoCounterScalTrig -> Integral();
    double selectionsZorroSel = histZorroSelCounterSelTOI -> Integral();
    double selectionsAnalysedTrig = histZorroSelCounterAnalysedTrig -> Integral();
    double collsBeforeFiltering = histEvSelCollisionBeforeFiltering -> Integral();
    double collsBeforeCuts = histEvSelCollisionBeforeCuts -> Integral();
    double collsAfterCuts = histEvSelCollisionAfterCuts -> Integral();
    double bcCounterTVX = histBcSelCounterTVX -> Integral();
    double bcCounterTVXafterBCcuts = histBcSelCounterTVXafterBCcuts -> Integral();
    double evColCounterTVX = histEvSelColCounterTVX -> Integral();

    TH1D *histLumiSummary = new TH1D("histLumiSummary", "", 9, 0, 9);
    histLumiSummary -> GetXaxis() -> SetBinLabel(1, "bcCounterTVX");
    histLumiSummary -> GetXaxis() -> SetBinLabel(2, "bcCounterTVXafterBCcuts");
    histLumiSummary -> GetXaxis() -> SetBinLabel(3, "bcCounterEfficiency");
    histLumiSummary -> GetXaxis() -> SetBinLabel(4, "nEvtsBcSel");
    histLumiSummary -> GetXaxis() -> SetBinLabel(5, "nEvtsBcSelAfterCuts");
    histLumiSummary -> GetXaxis() -> SetBinLabel(6, "collsBeforeCuts");
    histLumiSummary -> GetXaxis() -> SetBinLabel(7, "collsAfterCuts");
    histLumiSummary -> GetXaxis() -> SetBinLabel(8, "luminosityBcSel");
    histLumiSummary -> GetXaxis() -> SetBinLabel(9, "luminosityBcSelAfterCuts");

    std::cout << "****************************************************" << std::endl;
    std::cout << "nRuns = " << nRuns << std::endl;
    std::cout << "inspectedTVX = " << inspectedTVX << std::endl;
    //std::cout << "selectionsZorroInfo = " << selectionsZorroInfo << std::endl;
    std::cout << "scalers = " << scalers << std::endl;
    std::cout << "selectionsAnalysedTrig = " << selectionsAnalysedTrig << std::endl;
    std::cout << "selectionsZorroSel = " << selectionsZorroSel << std::endl;
    //std::cout << "collsBeforeFiltering = " << collsBeforeFiltering << std::endl;
    std::cout << "collsBeforeCuts = " << collsBeforeCuts << std::endl;
    std::cout << "collsAfterCuts = " << collsAfterCuts << std::endl;
    std::cout << "bcCounterTVX = " << bcCounterTVX << std::endl;
    std::cout << "bcCounterTVXafterBCcuts = " << bcCounterTVXafterBCcuts << std::endl;
    std::cout << "evColCounterTVX = " << evColCounterTVX << std::endl;
    std::cout << "****************************************************" << std::endl;

    const double ppCrossSection = 59.4e9; // pb-1
    double bcCounterEfficiency = bcCounterTVXafterBCcuts / bcCounterTVX;
    double nEvtsBcSel, nEvtsBcSelAfterCuts = -999;
    double luminosityBcSel, luminosityBcSelAfterCuts = -999;
    if (triggerMask.find("minBias") != std::string::npos) {
        nEvtsBcSel = (bcCounterTVX * (collsAfterCuts / collsBeforeCuts));
        nEvtsBcSelAfterCuts = (bcCounterTVXafterBCcuts * (collsAfterCuts / collsBeforeCuts));
        luminosityBcSel = nEvtsBcSel / ppCrossSection;
        luminosityBcSelAfterCuts = nEvtsBcSelAfterCuts / ppCrossSection;

        double luminosityBcSelMinBias = (bcCounterTVX * (collsAfterCuts / collsBeforeCuts)) / 59.4e9;
        double luminosityBcSelAfterCutsMinBias = (bcCounterTVXafterBCcuts * (collsAfterCuts / collsBeforeCuts)) / 59.4e9;
        double luminosityEvSelMinBias = (evColCounterTVX * (collsAfterCuts / collsBeforeCuts)) / 59.4e9;
        std::cout << "BC selection efficiency = " << bcCounterTVXafterBCcuts / bcCounterTVX << std::endl;
        std::cout << "N. events            = " << (bcCounterTVX * (collsAfterCuts / collsBeforeCuts)) << std::endl;
        std::cout << "Colls after cuts     = " << collsAfterCuts << std::endl;
        std::cout << "luminosity Min. Bias [bc-selection-task]                  = " << luminosityBcSelMinBias << " pb-1" << std::endl;
        std::cout << "luminosity Min. Bias after BC cuts [bc-selection-task]    = " << luminosityBcSelAfterCutsMinBias << " pb-1" << std::endl;
        std::cout << "luminosity Min. Bias [event-selection-task]               = " << luminosityEvSelMinBias << " pb-1" << std::endl;
    } else {
        double luminosityZorroInfo = (inspectedTVX * (selectionsZorroInfo / scalers)) / 59.4e9;
        double luminosityZorroSel = (inspectedTVX * (selectionsZorroSel / scalers)) / 59.4e9;
        double luminosityAnalysedTrig = (inspectedTVX * (selectionsAnalysedTrig / scalers)) / 59.4e9;
        std::cout << "BC selection efficiency = " << bcCounterTVXafterBCcuts / bcCounterTVX << std::endl;
        std::cout << "N. events TOI           = " << (inspectedTVX * (selectionsZorroSel / scalers)) << std::endl;
        std::cout << "N. events AnalysedTrig  = " << (inspectedTVX * (selectionsAnalysedTrig / scalers)) << std::endl;
        std::cout << "luminosityZorroInfo    = " << luminosityZorroInfo << " pb-1" << std::endl;
        std::cout << "luminosity TOI          = " << luminosityZorroSel << " pb-1" << std::endl;
        std::cout << "luminosity AnalysedTrig = " << luminosityAnalysedTrig << " pb-1" << std::endl;

        nEvtsBcSel = (inspectedTVX * (selectionsAnalysedTrig / scalers));
        nEvtsBcSelAfterCuts = (inspectedTVX * (selectionsZorroSel / scalers));
        luminosityBcSel = nEvtsBcSel / ppCrossSection;
        luminosityBcSelAfterCuts = nEvtsBcSelAfterCuts / ppCrossSection;
    }
    histLumiSummary -> SetBinContent(1, bcCounterTVX);
    histLumiSummary -> SetBinContent(2, bcCounterTVXafterBCcuts);
    histLumiSummary -> SetBinContent(3, bcCounterEfficiency);
    histLumiSummary -> SetBinContent(4, nEvtsBcSel);
    histLumiSummary -> SetBinContent(5, nEvtsBcSelAfterCuts);
    histLumiSummary -> SetBinContent(6, collsBeforeCuts);
    histLumiSummary -> SetBinContent(7, collsAfterCuts);
    histLumiSummary -> SetBinContent(8, luminosityBcSel);
    histLumiSummary -> SetBinContent(9, luminosityBcSelAfterCuts);

    TFile *fOut = new TFile(Form("data/%s/%s_%s_%s_luminosity.root", year.c_str(), period.c_str(), triggerMask.c_str(), assocType.c_str()), "RECREATE");
    histLumiSummary -> Write();
    fOut -> Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_normalization() {
    LoadStyle();

    //string fInLumiNameSkim2024StdAssoc = "data/2024/LHC24af_fDiMuon_std_assoc_luminosity.root";
    //string fInNameSkim2024StdAssoc = "data/2024/LHC24af_fDiMuon_std_assoc.root";

    string fInLumiNameSkim2024StdAssoc = "data/2024/LHC24_CBT_muon_fDiMuon_std_assoc_luminosity.root";
    string fInNameSkim2024StdAssoc = "data/2024/LHC24_CBT_muon_fDiMuon_std_assoc.root";

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
    TFile *fInLumiMinBias2024StdAssoc = TFile::Open("data/2024/LHC24_minBias_std_assoc_luminosity.root");
    TH1D *histLumiMinBias2024StdAssoc = (TH1D*) fInLumiMinBias2024StdAssoc -> Get("histLumiSummary");
    double nEvtsMinBias2024StdAssoc = histLumiMinBias2024StdAssoc -> GetBinContent(histLumiMinBias2024StdAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double lumiMinBias2024StdAssoc = histLumiMinBias2024StdAssoc -> GetBinContent(histLumiMinBias2024StdAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));

    TFile *fInMinBias2024StdAssoc = TFile::Open("data/2024/LHC24_minBias_std_assoc.root");
    TList *listMinBias2024StdAssoc = (TList*) fInMinBias2024StdAssoc -> Get("analysis-same-event-pairing/output");
    TList *listMinBias2024StdAssocSE = (TList*) listMinBias2024StdAssoc -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA");
    THnSparseD *histSparseMinBias2024StdAssoc = (THnSparseD*) listMinBias2024StdAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseMinBias2024StdAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseMinBias2024StdAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntMinBias2024StdAssoc = (TH1D*) histSparseMinBias2024StdAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");
    double countsMinBias2024StdAssoc = histProjIntMinBias2024StdAssoc -> Integral();
    histProjIntMinBias2024StdAssoc -> Scale(1. / (lumiMinBias2024StdAssoc));
    SetHist(histProjIntMinBias2024StdAssoc, kBlue+1, 24, 0.8, kBlue+1);

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
    std::cout << "[Min. Bias, std. assoc.]" << std::endl;
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void RetrieveTriggerInfo(TString dirName = "path/to/file", bool fromAlien = true, string triggerMask = "fDiMuon", double counters[10] = 0) {

    TString fInName;
    double collisionsBeforeFiltering = 0;
    double collisionsBeforeCuts = 0;
    double collisionsAfterCuts = 0;
    double infoTVX = 0;
    double infoScalTrig = 0;
    double infoSelTrig = 0;
    double selTOI = 0;
    double analysedTriggers = 0;
    double counterTVX = 0;
    double counterTVXafterBCcuts = 0;
    double colCounterTVX = 0;

    if (fromAlien) {
        TGrid::Connect("alien://");
        if (!dirName.Contains("alien://")) {
            dirName = "alien://" + dirName;
        }
    }

    for (int iDir = 0;iDir < 1000;iDir++) {
        fInName.Clear();
        if (iDir < 10) {fInName = dirName + Form("/AOD/00%i/AnalysisResults.root", iDir+1);}
        if (iDir >= 10 && iDir < 100) {fInName = dirName + Form("/AOD/0%i/AnalysisResults.root", iDir+1);}

        TFile *fIn = TFile::Open(fInName.Data());
        if (!fIn || fIn -> IsZombie()) {
            std::cout << "ERROR! THE FILE DOES NOT EXIST!" << std::endl;
            break;
        } else {
            std::cout << "THE FILE EXISTS!" << std::endl;
        }

        for (auto dirKey : *fIn -> GetListOfKeys()) {
            if (TString(dirKey -> GetName()).Contains("table-maker")) {
                TList *list = (TList*) fIn -> Get("table-maker/Statistics");
                TH2D *histZorroInfo = (TH2D*) list -> FindObject("ZorroInfo");
                TH2D *histZorroSel = (TH2D*) list -> FindObject("ZorroSel");
                TH2D *histEventStats = (TH2D*) list -> FindObject("EventStats");
    
                //std::cout << "inspected TVX: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin("inspected TVX")) << std::endl;
                //std::cout << triggerMask << " scalers: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str()))) << std::endl;
                //std::cout << triggerMask << " fSingleMuLow selections: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str()))) << std::endl;
                //std::cout << triggerMask << " fSingleMuLow scalers (SEL): " << histZorroSel -> GetBinContent(1, histZorroSel -> GetYaxis() -> FindBin(Form("%s", triggerMask.c_str()))) << std::endl;

                infoTVX += histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin("inspected TVX"));
                infoScalTrig += histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str())));
                infoSelTrig += histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str())));
                selTOI += histZorroSel -> GetBinContent(1, histZorroSel -> GetYaxis() -> FindBin(Form("%s", triggerMask.c_str())));
                analysedTriggers += histZorroSel -> GetBinContent(1, histZorroSel -> GetYaxis() -> FindBin(Form("%s AnalysedTriggers", triggerMask.c_str())));

                collisionsBeforeFiltering += histEventStats -> GetBinContent(2, histEventStats -> GetYaxis() -> FindBin("Total"));
                collisionsBeforeCuts += histEventStats -> GetBinContent(3, histEventStats -> GetYaxis() -> FindBin("Total"));
                collisionsAfterCuts += histEventStats -> GetBinContent(4, histEventStats -> GetYaxis() -> FindBin("Total"));
            }
            if (TString(dirKey -> GetName()).Contains("bc-selection-task")) {
                TH1D *histCounterTVX = (TH1D*) fIn -> Get("bc-selection-task/hCounterTVX");
                TH1D *histCounterTVXafterBCcuts = (TH1D*) fIn -> Get("bc-selection-task/hCounterTVXafterBCcuts");

                counterTVX += histCounterTVX -> GetBinContent(1);
                counterTVXafterBCcuts += histCounterTVXafterBCcuts -> GetBinContent(1);
            }

            if (TString(dirKey -> GetName()).Contains("event-selection-task")) {
                TH1D *histColCounterTVX = (TH1D*) fIn -> Get("event-selection-task/hColCounterTVX");

                colCounterTVX += histColCounterTVX -> GetBinContent(1);
            }
        }

        counters[0] = infoTVX;
        counters[1] = infoScalTrig;
        counters[2] = infoSelTrig;
        counters[3] = selTOI;
        counters[4] = analysedTriggers;
        counters[5] = collisionsBeforeFiltering;
        counters[6] = collisionsBeforeCuts;
        counters[7] = collisionsAfterCuts;
        counters[8] = counterTVX;
        counters[9] = counterTVXafterBCcuts;
        counters[10] = colCounterTVX;

        fIn -> Close();
    }
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