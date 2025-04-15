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

void RetrieveTriggerInfo(TString , bool , string, double [10]);

void get_info_from_zorro(string year = "2024", string period = "LHC24", string triggerMask = "fDiMuon") {
    ifstream fInAlienInputPath(Form("run_lists/%s/%s/alien_input_path.txt", year.c_str(), triggerMask.c_str()));
    if (!fInAlienInputPath.is_open()) {
        std::cout << "Error opening the file!" << std::endl;
        return;
    }

    ifstream fInAlienRunList(Form("run_lists/%s/%s/alien_run_list.txt", year.c_str(), triggerMask.c_str()));
    if (!fInAlienRunList.is_open()) {
        std::cout << "Error opening the file!" << std::endl;
        return;
    }

    double counters[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    std::vector <double> vecInfoCounterTVX;
    std::vector <double> vecInfoCounterScalTrig;
    std::vector <double> vecInfoCounterSelTrig;
    std::vector <double> vecSelCounterSelTOI;
    std::vector <double> vecSelCounterAnalysedTrig;
    std::vector <string> vecRunList;

    string alienRunListName;
    while (getline(fInAlienRunList, alienRunListName)) {
        vecRunList.push_back(alienRunListName);
    }

    TH1D *histZorroInfoCounterTVX = new TH1D("histZorroInfoCounterTVX", "", vecRunList.size(), 0, vecRunList.size());
    TH1D *histZorroInfoCounterScalTrig = new TH1D("histZorroInfoCounterScalTrig", "", vecRunList.size(), 0, vecRunList.size());
    TH1D *histZorroInfoCounterSelTrig = new TH1D("histZorroInfoCounterSelTrig", "", vecRunList.size(), 0, vecRunList.size());
    TH1D *histZorroSelCounterSelTOI = new TH1D("histZorroSelCounterSelTOI", "", vecRunList.size(), 0, vecRunList.size());
    TH1D *histZorroSelCounterAnalysedTrig = new TH1D("histZorroSelCounterAnalysedTrig", "", vecRunList.size(), 0, vecRunList.size());
    TH1D *histEvSelCollisionBeforeFiltering = new TH1D("histEvSelCollisionBeforeFiltering", "", vecRunList.size(), 0, vecRunList.size());
    TH1D *histEvSelCollisionBeforeCuts = new TH1D("histEvSelCollisionBeforeCuts", "", vecRunList.size(), 0, vecRunList.size());
    TH1D *histEvSelCollisionAfterCuts = new TH1D("histEvSelCollisionAfterCuts", "", vecRunList.size(), 0, vecRunList.size());
    TH1D *histBcSelCounterTVX = new TH1D("histBcSelCounterTVX", "", vecRunList.size(), 0, vecRunList.size());
    TH1D *histBcSelCounterTVXafterBCcuts = new TH1D("histBcSelCounterTVXafterBCcuts", "", vecRunList.size(), 0, vecRunList.size());

    string alienPathName;
    string runNumber;
    int runCounter = 0;
    while (getline(fInAlienInputPath, alienPathName)) {
        RetrieveTriggerInfo(alienPathName.c_str(), true, triggerMask, counters);
        //RetrieveTriggerInfo(Form("%s/AOD", alienPathName.c_str()), true, triggerMask, counters);
        
        runNumber = vecRunList.at(runCounter).c_str();

        histZorroInfoCounterTVX -> SetBinContent(runCounter+1, counters[0]);
        histZorroInfoCounterTVX -> GetXaxis() -> SetBinLabel(runCounter+1, runNumber.c_str());

        histZorroInfoCounterScalTrig -> SetBinContent(runCounter+1, counters[1]);
        histZorroInfoCounterScalTrig -> GetXaxis() -> SetBinLabel(runCounter+1, runNumber.c_str());

        histZorroInfoCounterSelTrig -> SetBinContent(runCounter+1, counters[2]);
        histZorroInfoCounterSelTrig -> GetXaxis() -> SetBinLabel(runCounter+1, runNumber.c_str());

        histZorroSelCounterSelTOI -> SetBinContent(runCounter+1, counters[3]);
        histZorroSelCounterSelTOI -> GetXaxis() -> SetBinLabel(runCounter+1, runNumber.c_str());

        histZorroSelCounterAnalysedTrig -> SetBinContent(runCounter+1, counters[4]);
        histZorroSelCounterAnalysedTrig -> GetXaxis() -> SetBinLabel(runCounter+1, runNumber.c_str());

        histEvSelCollisionBeforeFiltering -> SetBinContent(runCounter+1, counters[5]);
        histEvSelCollisionBeforeFiltering -> GetXaxis() -> SetBinLabel(runCounter+1, runNumber.c_str());

        histEvSelCollisionBeforeCuts -> SetBinContent(runCounter+1, counters[6]);
        histEvSelCollisionBeforeCuts -> GetXaxis() -> SetBinLabel(runCounter+1, runNumber.c_str());

        histEvSelCollisionAfterCuts -> SetBinContent(runCounter+1, counters[7]);
        histEvSelCollisionAfterCuts -> GetXaxis() -> SetBinLabel(runCounter+1, runNumber.c_str());

        histBcSelCounterTVX -> SetBinContent(runCounter+1, counters[8]);
        histBcSelCounterTVX -> GetXaxis() -> SetBinLabel(runCounter+1, runNumber.c_str());

        histBcSelCounterTVXafterBCcuts -> SetBinContent(runCounter+1, counters[9]);
        histBcSelCounterTVXafterBCcuts -> GetXaxis() -> SetBinLabel(runCounter+1, runNumber.c_str());

        runCounter++;
    }

    TFile *fOut = new TFile(Form("run_lists/%s/%s_%s_trigger_summary.root", year.c_str(), period.c_str(), triggerMask.c_str()), "RECREATE");
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
    fOut -> Close();
}

////////////////////////////////////////////////////////////////////////////////
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

        auto fIn = TFile::Open(fInName.Data());
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

        fIn -> Close();
    }
}