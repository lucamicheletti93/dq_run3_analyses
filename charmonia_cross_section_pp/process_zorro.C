#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TTreeReader.h>
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
        fOutName = Form("run_lists/%s/%s_%s_%s_trigger_summary.root", year.c_str(), period.c_str(), triggerMask.c_str(), assocType.c_str());
    } else {
        fOutName = Form("run_lists/%s/%s_%s_%s_trigger_summary.root", year.c_str(), subPeriod.c_str(), triggerMask.c_str(), assocType.c_str());
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

    TH1D *histZorroInfoCounterTVX;
    TH1D *histZorroInfoCounterScalTrig;
    TH1D *histZorroInfoCounterSelTrig;
    TH1D *histZorroSelCounterSelTOI;
    TH1D *histZorroSelCounterAnalysedTrig;
    TH1D *histEvSelCollisionBeforeFiltering;
    TH1D *histEvSelCollisionBeforeCuts;
    TH1D *histEvSelCollisionAfterCuts;
    TH1D *histBcSelCounterTVX;
    TH1D *histBcSelCounterTVXafterBCcuts;
    TH1D *histEvSelColCounterTVX;

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

    TH1D *histLumiSummary = new TH1D("histLumiSummary", "", 7, 0, 7);
    histLumiSummary -> GetXaxis() -> SetBinLabel(1, "bcCounterTVX");
    histLumiSummary -> GetXaxis() -> SetBinLabel(2, "bcCounterTVXafterBCcuts");
    histLumiSummary -> GetXaxis() -> SetBinLabel(3, "bcCounterEfficiency");
    histLumiSummary -> GetXaxis() -> SetBinLabel(4, "nEvtsBcSel");
    histLumiSummary -> GetXaxis() -> SetBinLabel(5, "nEvtsBcSelAfterCuts");
    histLumiSummary -> GetXaxis() -> SetBinLabel(6, "luminosityBcSel");
    histLumiSummary -> GetXaxis() -> SetBinLabel(7, "luminosityBcSelAfterCuts");

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
    histLumiSummary -> SetBinContent(6, luminosityBcSel);
    histLumiSummary -> SetBinContent(7, luminosityBcSelAfterCuts);

    TFile *fOut = new TFile(Form("data/%s/%s_%s_%s_luminosity.root", year.c_str(), period.c_str(), triggerMask.c_str(), assocType.c_str()), "RECREATE");
    histLumiSummary -> Write();
    fOut -> Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void check_normalization() {
    LoadStyle();

    // Skimmed, std. assoc.
    TFile *fInLumiSkimStdAssoc = TFile::Open("data/2024/LHC24af_fDiMuon_std_assoc_luminosity.root");
    TH1D *histLumiSkimStdAssoc = (TH1D*) fInLumiSkimStdAssoc -> Get("histLumiSummary");
    double bcSelEffSkimStdAssoc = histLumiSkimStdAssoc -> GetBinContent(histLumiSkimStdAssoc -> GetXaxis() -> FindBin("bcCounterEfficiency"));
    double nEvtsSkimStdAssoc = histLumiSkimStdAssoc -> GetBinContent(histLumiSkimStdAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double lumiSkimStdAssoc = histLumiSkimStdAssoc -> GetBinContent(histLumiSkimStdAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));

    TFile *fInSkimStdAssoc = TFile::Open("data/2024/LHC24af_fDiMuon_std_assoc.root");
    TList *listSkimStdAssoc = (TList*) fInSkimStdAssoc -> Get("analysis-same-event-pairing/output");
    TList *listSkimStdAssocSE = (TList*) listSkimStdAssoc -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA");
    THnSparseD *histSparseSkimStdAssoc = (THnSparseD*) listSkimStdAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseSkimStdAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseSkimStdAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntSkimStdAssoc = (TH1D*) histSparseSkimStdAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");
    TH1D *histProjIntSkimStdAssocNormToColl = (TH1D*) histProjIntSkimStdAssoc -> Clone("histProjIntSkimStdAssocNormToColl");
    histProjIntSkimStdAssoc -> Scale(1 / (bcSelEffSkimStdAssoc * lumiSkimStdAssoc));
    SetHist(histProjIntSkimStdAssoc, kRed+1, 20, 0.8, kRed+1);
    histProjIntSkimStdAssocNormToColl -> Scale(1. / (bcSelEffSkimStdAssoc * nEvtsSkimStdAssoc));
    SetHist(histProjIntSkimStdAssocNormToColl, kRed+1, 20, 0.8, kRed+1);

    // Skimmed, time assoc.
    // WARNING! for the time association you have to take into account the duplication using zorro
    TFile *fInLumiSkimTimeAssoc = TFile::Open("data/2024/LHC24af_fDiMuon_time_assoc_luminosity.root");
    TH1D *histLumiSkimTimeAssoc = (TH1D*) fInLumiSkimTimeAssoc -> Get("histLumiSummary");
    double bcSelEffSkimTimeAssoc = histLumiSkimTimeAssoc -> GetBinContent(histLumiSkimTimeAssoc -> GetXaxis() -> FindBin("bcCounterEfficiency"));
    double nEvtsSkimTimeAssoc = histLumiSkimTimeAssoc -> GetBinContent(histLumiSkimTimeAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double nEvtsAfterCutsSkimTimeAssoc = histLumiSkimTimeAssoc -> GetBinContent(histLumiSkimTimeAssoc -> GetXaxis() -> FindBin("nEvtsBcSelAfterCuts"));
    double lumiSkimTimeAssoc = histLumiSkimTimeAssoc -> GetBinContent(histLumiSkimTimeAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));

    TFile *fInSkimTimeAssoc = TFile::Open("data/2024/LHC24af_fDiMuon_time_assoc.root");
    TList *listSkimTimeAssoc = (TList*) fInSkimTimeAssoc -> Get("analysis-same-event-pairing/output");
    TList *listSkimTimeAssocSE = (TList*) listSkimTimeAssoc -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA");
    THnSparseD *histSparseSkimTimeAssoc = (THnSparseD*) listSkimTimeAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseSkimTimeAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseSkimTimeAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntSkimTimeAssoc = (TH1D*) histSparseSkimTimeAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");
    TH1D *histProjIntSkimTimeAssocNormToColl = (TH1D*) histProjIntSkimTimeAssoc -> Clone("histProjIntSkimTimeAssocNormToColl");
    histProjIntSkimTimeAssoc -> Scale((nEvtsAfterCutsSkimTimeAssoc / nEvtsSkimTimeAssoc) / (bcSelEffSkimTimeAssoc * lumiSkimTimeAssoc));
    SetHist(histProjIntSkimTimeAssoc, kRed+1, 24, 0.8, kRed+1);
    histProjIntSkimTimeAssocNormToColl -> Scale((nEvtsAfterCutsSkimTimeAssoc / nEvtsSkimTimeAssoc) / (bcSelEffSkimTimeAssoc * nEvtsSkimTimeAssoc));
    SetHist(histProjIntSkimTimeAssocNormToColl, kRed+1, 24, 0.8, kRed+1);

    // Min. Bias, std. assoc.
    // WARNING! event selection already computed, bcSelCutEfficiecny not taken into account (already included in the event selection)
    // this correction is necessary in the trigger for inspectedTVX, for minBias it should not be necessary
    TFile *fInLumiMinBiasStdAssoc = TFile::Open("data/2024/LHC24_minBias_std_assoc_luminosity.root");
    TH1D *histLumiMinBiasStdAssoc = (TH1D*) fInLumiMinBiasStdAssoc -> Get("histLumiSummary");
    double nEvtsMinBiasStdAssoc = histLumiMinBiasStdAssoc -> GetBinContent(histLumiMinBiasStdAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double lumiMinBiasStdAssoc = histLumiMinBiasStdAssoc -> GetBinContent(histLumiMinBiasStdAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));

    TFile *fInMinBiasStdAssoc = TFile::Open("data/2024/LHC24_minBias_std_assoc.root");
    TList *listMinBiasStdAssoc = (TList*) fInMinBiasStdAssoc -> Get("analysis-same-event-pairing/output");
    TList *listMinBiasStdAssocSE = (TList*) listMinBiasStdAssoc -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA");
    THnSparseD *histSparseMinBiasStdAssoc = (THnSparseD*) listMinBiasStdAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseMinBiasStdAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseMinBiasStdAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntMinBiasStdAssoc = (TH1D*) histSparseMinBiasStdAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");
    TH1D *histProjIntMinBiasStdAssocNormToColl = (TH1D*) histProjIntMinBiasStdAssoc -> Clone("histProjIntMinBiasStdAssocNormToColl");
    double countsMinBiasStdAssoc = histProjIntMinBiasStdAssoc -> Integral();
    histProjIntMinBiasStdAssoc -> Scale(1. / (lumiMinBiasStdAssoc));
    SetHist(histProjIntMinBiasStdAssoc, kBlue+1, 20, 0.8, kBlue+1);
    histProjIntMinBiasStdAssocNormToColl -> Scale(1. / nEvtsMinBiasStdAssoc);
    SetHist(histProjIntMinBiasStdAssocNormToColl, kBlue+1, 20, 0.8, kBlue+1);

    // Min. Bias, std. assoc. + sel8
    TFile *fInLumiMinBiasSel8StdAssoc = TFile::Open("data/2024/LHC24_minBiasSel8_std_assoc_luminosity.root");
    TH1D *histLumiMinBiasSel8StdAssoc = (TH1D*) fInLumiMinBiasSel8StdAssoc -> Get("histLumiSummary");
    double nEvtsMinBiasSel8StdAssoc = histLumiMinBiasSel8StdAssoc -> GetBinContent(histLumiMinBiasSel8StdAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double lumiMinBiasSel8StdAssoc = histLumiMinBiasSel8StdAssoc -> GetBinContent(histLumiMinBiasSel8StdAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));

    TFile *fInMinBiasSel8StdAssoc = TFile::Open("data/2024/LHC24_minBiasSel8_std_assoc.root");
    TList *listMinBiasSel8StdAssoc = (TList*) fInMinBiasSel8StdAssoc -> Get("analysis-same-event-pairing/output");
    TList *listMinBiasSel8StdAssocSE = (TList*) listMinBiasSel8StdAssoc -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA");
    THnSparseD *histSparseMinBiasSel8StdAssoc = (THnSparseD*) listMinBiasSel8StdAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseMinBiasSel8StdAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseMinBiasSel8StdAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntMinBiasSel8StdAssoc = (TH1D*) histSparseMinBiasSel8StdAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");
    TH1D *histProjIntMinBiasSel8StdAssocNormToColl = (TH1D*) histProjIntMinBiasSel8StdAssoc -> Clone("histProjIntMinBiasSel8StdAssocNormToColl");
    double countsMinBiasSel8StdAssoc = histProjIntMinBiasSel8StdAssoc -> Integral();
    histProjIntMinBiasSel8StdAssoc -> Scale(1. / (lumiMinBiasSel8StdAssoc));
    SetHist(histProjIntMinBiasSel8StdAssoc, kAzure+4, 20, 0.8, kAzure+4);
    histProjIntMinBiasSel8StdAssocNormToColl -> Scale(1. / nEvtsMinBiasSel8StdAssoc);
    SetHist(histProjIntMinBiasSel8StdAssocNormToColl, kAzure+4, 20, 0.8, kAzure+4);

    // Min. Bias, time assoc.
    TFile *fInLumiMinBiasTimeAssoc = TFile::Open("data/2024/LHC24_minBias_time_assoc_luminosity.root");
    TH1D *histLumiMinBiasTimeAssoc = (TH1D*) fInLumiMinBiasTimeAssoc -> Get("histLumiSummary");
    double nEvtsMinBiasTimeAssoc = histLumiMinBiasTimeAssoc -> GetBinContent(histLumiMinBiasTimeAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double lumiMinBiasTimeAssoc = histLumiMinBiasTimeAssoc -> GetBinContent(histLumiMinBiasTimeAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));

    TFile *fInMinBiasTimeAssoc = TFile::Open("data/2024/LHC24_minBias_time_assoc.root");
    TList *listMinBiasTimeAssoc = (TList*) fInMinBiasTimeAssoc -> Get("analysis-same-event-pairing/output");
    TList *listMinBiasTimeAssocSE = (TList*) listMinBiasTimeAssoc -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA");
    THnSparseD *histSparseMinBiasTimeAssoc = (THnSparseD*) listMinBiasTimeAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseMinBiasTimeAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseMinBiasTimeAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntMinBiasTimeAssoc = (TH1D*) histSparseMinBiasTimeAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");
    TH1D *histProjIntMinBiasTimeAssocNormToColl = (TH1D*) histProjIntMinBiasTimeAssoc -> Clone("histProjIntMinBiasTimeAssocNormToColl");

    double countsMinBiasTimeAssoc = histProjIntMinBiasTimeAssoc -> Integral();
    double duplCorrFactorMinBiasTimeAssoc = (countsMinBiasStdAssoc / nEvtsMinBiasStdAssoc) / (countsMinBiasTimeAssoc / nEvtsMinBiasTimeAssoc);
    
    histProjIntMinBiasTimeAssoc -> Scale(duplCorrFactorMinBiasTimeAssoc / (lumiMinBiasTimeAssoc));
    SetHist(histProjIntMinBiasTimeAssoc, kBlue+1, 24, 0.8, kBlue+1);
    histProjIntMinBiasTimeAssocNormToColl -> Scale(duplCorrFactorMinBiasTimeAssoc / nEvtsMinBiasTimeAssoc);
    SetHist(histProjIntMinBiasTimeAssocNormToColl, kBlue+1, 24, 0.8, kBlue+1);

    // Min. Bias, time assoc. + sel8
    TFile *fInLumiMinBiasSel8TimeAssoc = TFile::Open("data/2024/LHC24_minBiasSel8_time_assoc_luminosity.root");
    TH1D *histLumiMinBiasSel8TimeAssoc = (TH1D*) fInLumiMinBiasSel8TimeAssoc -> Get("histLumiSummary");
    double nEvtsMinBiasSel8TimeAssoc = histLumiMinBiasSel8TimeAssoc -> GetBinContent(histLumiMinBiasSel8TimeAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double lumiMinBiasSel8TimeAssoc = histLumiMinBiasSel8TimeAssoc -> GetBinContent(histLumiMinBiasSel8TimeAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));

    TFile *fInMinBiasSel8TimeAssoc = TFile::Open("data/2024/LHC24_minBiasSel8_time_assoc.root");
    TList *listMinBiasSel8TimeAssoc = (TList*) fInMinBiasSel8TimeAssoc -> Get("analysis-same-event-pairing/output");
    TList *listMinBiasSel8TimeAssocSE = (TList*) listMinBiasSel8TimeAssoc -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA");
    THnSparseD *histSparseMinBiasSel8TimeAssoc = (THnSparseD*) listMinBiasSel8TimeAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseMinBiasSel8TimeAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseMinBiasSel8TimeAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntMinBiasSel8TimeAssoc = (TH1D*) histSparseMinBiasSel8TimeAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");
    TH1D *histProjIntMinBiasSel8TimeAssocNormToColl = (TH1D*) histProjIntMinBiasSel8TimeAssoc -> Clone("histProjIntMinBiasSel8TimeAssocNormToColl");

    double countsMinBiasSel8TimeAssoc = histProjIntMinBiasSel8TimeAssoc -> Integral();
    double duplCorrFactorMinBiasSel8TimeAssoc = (countsMinBiasSel8StdAssoc / nEvtsMinBiasSel8StdAssoc) / (countsMinBiasSel8TimeAssoc / nEvtsMinBiasSel8TimeAssoc);
    
    histProjIntMinBiasSel8TimeAssoc -> Scale(duplCorrFactorMinBiasSel8TimeAssoc / (lumiMinBiasSel8TimeAssoc));
    SetHist(histProjIntMinBiasSel8TimeAssoc, kAzure+4, 24, 0.8, kAzure+4);
    histProjIntMinBiasSel8TimeAssocNormToColl -> Scale(duplCorrFactorMinBiasSel8TimeAssoc / nEvtsMinBiasSel8TimeAssoc);
    SetHist(histProjIntMinBiasSel8TimeAssocNormToColl, kAzure+4, 24, 0.8, kAzure+4);

    //std::cout << "Min Bias dimuons time assoc. = " << countsMinBiasTimeAssoc << " ; coll. after cuts = " << collsAfterCutsMinBiasTimeAssoc << std::endl;
    //std::cout << "Min Bias dimuons std. assoc. = " << countsMinBiasStdAssoc << " ; coll. after cuts = " << collsAfterCutsMinBiasStdAssoc << std::endl;
    //std::cout << "Correction factor from inv. mass = " << (countsMinBiasTimeAssoc / collsAfterCutsMinBiasTimeAssoc) / (countsMinBiasStdAssoc / collsAfterCutsMinBiasStdAssoc) << std::endl;
    //std::cout << "Correction factor from zorro     = " << lumiSkimTimeAssoc / lumiSkimStdAssoc << std::endl;

    TCanvas *canvasCompMassSpectra = new TCanvas("canvasCompMassSpectra", "", 800, 600);
    gPad -> SetLogy(true);
    histProjIntSkimStdAssoc -> SetTitle("");
    histProjIntSkimStdAssoc -> GetXaxis() -> SetRangeUser(2, 5);
    histProjIntSkimStdAssoc -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c})");
    histProjIntSkimStdAssoc -> GetYaxis() -> SetTitle("Counts");
    histProjIntSkimStdAssoc -> Draw("EP");
    histProjIntSkimTimeAssoc -> Draw("EP SAME");
    histProjIntMinBiasStdAssoc -> Draw("EP SAME");
    histProjIntMinBiasTimeAssoc -> Draw("EP SAME");
    //histProjIntMinBiasSel8StdAssoc -> Draw("EP SAME");
    //histProjIntMinBiasSel8TimeAssoc -> Draw("EP SAME");

    TLegend *legendCompMassSpectra = new TLegend(0.6, 0.57, 0.8, 0.82, " ", "brNDC");
    SetLegend(legendCompMassSpectra);
    legendCompMassSpectra -> AddEntry(histProjIntSkimStdAssoc, "Skimmed, std. assoc", "EP");
    legendCompMassSpectra -> AddEntry(histProjIntSkimTimeAssoc, "Skimmed, time assoc", "EP");
    legendCompMassSpectra -> AddEntry(histProjIntMinBiasStdAssoc, "MB, std. assoc", "EP");
    legendCompMassSpectra -> AddEntry(histProjIntMinBiasTimeAssoc, "MB, time assoc", "EP");
    //legendCompMassSpectra -> AddEntry(histProjIntMinBiasSel8StdAssoc, "MB, std. assoc + sel8", "EP");
    //legendCompMassSpectra -> AddEntry(histProjIntMinBiasSel8TimeAssoc, "MB, time assoc + sel8", "EP");
    legendCompMassSpectra -> Draw();

    // Compute the ratio
    TH1D *histRatioIntSkimStdSkimTimeAssoc = (TH1D*) histProjIntSkimStdAssoc -> Clone("histRatioIntSkimStdSkimTimeAssoc");
    histRatioIntSkimStdSkimTimeAssoc -> Divide(histProjIntSkimTimeAssoc);
    SetHist(histRatioIntSkimStdSkimTimeAssoc, kRed+1, 20, 0.8, kRed+1);

    TH1D *histRatioIntMinBiasStdMinBiasTimeAssoc = (TH1D*) histProjIntMinBiasStdAssoc -> Clone("histRatioIntMinBiasStdMinBiasTimeAssoc");
    histRatioIntMinBiasStdMinBiasTimeAssoc -> Divide(histProjIntMinBiasTimeAssoc);
    SetHist(histRatioIntMinBiasStdMinBiasTimeAssoc, kBlue+1, 20, 0.8, kBlue+1);

    TH1D *histRatioIntSkimStdMinBiasStdAssoc = (TH1D*) histProjIntSkimStdAssoc -> Clone("histRatioIntSkimStdMinBiasStdAssoc");
    histRatioIntSkimStdMinBiasStdAssoc -> Divide(histProjIntMinBiasStdAssoc);
    SetHist(histRatioIntSkimStdMinBiasStdAssoc, kRed+1, 24, 0.8, kRed+1);

    TH1D *histRatioIntSkimTimeMinBiasTimeAssoc = (TH1D*) histProjIntSkimTimeAssoc -> Clone("histRatioIntSkimTimeMinBiasTimeAssoc");
    histRatioIntSkimTimeMinBiasTimeAssoc -> Divide(histProjIntMinBiasTimeAssoc);
    SetHist(histRatioIntSkimTimeMinBiasTimeAssoc, kBlue+1, 24, 0.8, kBlue+1);

    TCanvas *canvasRatioMassSpectra = new TCanvas("canvasRatioMassSpectra", "", 800, 600);
    histRatioIntSkimStdSkimTimeAssoc -> GetXaxis() -> SetRangeUser(2, 5);
    histRatioIntSkimStdSkimTimeAssoc -> GetYaxis() -> SetRangeUser(0.5, 1.5);
    histRatioIntSkimStdSkimTimeAssoc -> GetYaxis() -> SetTitle("Ratio");
    histRatioIntSkimStdSkimTimeAssoc -> Draw("EP");
    histRatioIntSkimStdMinBiasStdAssoc -> Draw("EP SAME");
    histRatioIntSkimTimeMinBiasTimeAssoc -> Draw("EP SAME");
    histRatioIntMinBiasStdMinBiasTimeAssoc -> Draw("EP SAME");

    TLine *lineUnity = new TLine(2, 1, 5, 1);
    lineUnity -> SetLineColor(kGray+2);
    lineUnity -> SetLineStyle(kDashed);
    lineUnity -> Draw();

    TLegend *legendRatioMassSpectra = new TLegend(0.20, 0.20, 0.60, 0.40, " ", "brNDC");
    SetLegend(legendRatioMassSpectra);
    legendRatioMassSpectra -> AddEntry(histRatioIntSkimStdSkimTimeAssoc, "skimmed std. / skimmed time", "EP");
    legendRatioMassSpectra -> AddEntry(histRatioIntSkimStdMinBiasStdAssoc, "skimmed std. / Min. Bias std.", "EP");
    legendRatioMassSpectra -> AddEntry(histRatioIntSkimTimeMinBiasTimeAssoc, "skimmed time / Min. Bias time", "EP");
    legendRatioMassSpectra -> AddEntry(histRatioIntMinBiasStdMinBiasTimeAssoc, "Min. Bias std. / Min. Bias time", "EP");
    legendRatioMassSpectra -> Draw();

    TCanvas *canvasCompMassSpectraNormalized = new TCanvas("canvasCompMassSpectraNormalized", "", 800, 600);
    gPad -> SetLogy(true);
    histProjIntSkimStdAssocNormToColl -> SetTitle("");
    histProjIntSkimStdAssocNormToColl -> GetXaxis() -> SetRangeUser(2, 5);
    histProjIntSkimStdAssocNormToColl -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c})");
    //histProjIntSkimStdAssocNormToColl -> GetYaxis() -> SetRangeUser(1e-9, 1e-5);
    histProjIntSkimStdAssocNormToColl -> GetYaxis() -> SetTitle("Normalized counts");
    histProjIntSkimStdAssocNormToColl -> Draw("EP");
    histProjIntSkimTimeAssocNormToColl -> Draw("EP SAME");
    histProjIntMinBiasStdAssocNormToColl -> Draw("EP SAME");
    histProjIntMinBiasSel8StdAssocNormToColl -> Draw("EP SAME");
    histProjIntMinBiasTimeAssocNormToColl -> Draw("EP SAME");
    histProjIntMinBiasSel8TimeAssocNormToColl -> Draw("EP SAME");

    TLegend *legendCompMassSpectraNormalized = new TLegend(0.6, 0.57, 0.8, 0.82, " ", "brNDC");
    SetLegend(legendCompMassSpectraNormalized);
    legendCompMassSpectraNormalized -> AddEntry(histProjIntSkimStdAssocNormToColl, "Skimmed, std. assoc", "EP");
    legendCompMassSpectraNormalized -> AddEntry(histProjIntSkimTimeAssocNormToColl, "Skimmed, time assoc", "EP");
    legendCompMassSpectraNormalized -> AddEntry(histProjIntMinBiasStdAssocNormToColl, "MB, std. assoc", "EP");
    legendCompMassSpectraNormalized -> AddEntry(histProjIntMinBiasTimeAssocNormToColl, "MB, time assoc", "EP");
    legendCompMassSpectraNormalized -> AddEntry(histProjIntMinBiasSel8StdAssocNormToColl, "MB, std. assoc + sel8", "EP");
    legendCompMassSpectraNormalized -> AddEntry(histProjIntMinBiasSel8TimeAssocNormToColl, "MB, time assoc + sel8", "EP");
    legendCompMassSpectraNormalized -> Draw();

    canvasCompMassSpectra -> SaveAs("figures/CompMassSpectra.pdf");
    canvasRatioMassSpectra -> SaveAs("figures/RatioMassSpectra.pdf");
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
        if (iDir >= 10 && iDir < 100) {fInName = dirName + "/AOD/0" + Form("/AOD/0%i/AnalysisResults.root", iDir+1);}

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