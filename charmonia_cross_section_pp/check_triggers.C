#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TGrid.h>
#include <iostream>
#include <vector>
#include <filesystem>

void RetrieveTriggerInfo(TString , bool , string, double [5]);

void check_triggers() {
    string year = "2024";
    string period = "LHC24";
    string triggerMask = "fDiMuon"; // fDiMuon,fSingleMuLow

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

    double counters[] = {0, 0, 0, 0, 0};

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

    TH1D *histInfoCounterTVX = new TH1D("histInfoCounterTVX", "", vecInfoCounterTVX.size(), 0, vecInfoCounterTVX.size());
    TH1D *histInfoCounterScalTrig = new TH1D("histInfoCounterScalTrig", "", vecInfoCounterScalTrig.size(), 0, vecInfoCounterScalTrig.size());
    TH1D *histInfoCounterSelTrig = new TH1D("histInfoCounterSelTrig", "", vecInfoCounterSelTrig.size(), 0, vecInfoCounterSelTrig.size());
    TH1D *histSelCounterSelTOI = new TH1D("histSelCounterSelTOI", "", vecInfoCounterSelTrig.size(), 0, vecInfoCounterSelTrig.size());
    TH1D *histSelCounterAnalysedTrig = new TH1D("histSelCounterAnalysedTrig", "", vecSelCounterAnalysedTrig.size(), 0, vecSelCounterAnalysedTrig.size());

    string alienPathName;
    while (getline(fInAlienInputPath, alienPathName)) {
        RetrieveTriggerInfo(alienPathName.c_str(), true, triggerMask, counters);
        //RetrieveTriggerInfo(Form("%s/AOD", alienPathName.c_str()), true, triggerMask, counters);
        vecInfoCounterTVX.push_back(counters[0]);
        vecInfoCounterScalTrig.push_back(counters[1]);
        vecInfoCounterSelTrig.push_back(counters[2]);
        vecSelCounterSelTOI.push_back(counters[3]);
        vecSelCounterAnalysedTrig.push_back(counters[4]);
    }

    for (int iRun = 0;iRun < vecInfoCounterTVX.size();iRun++) {
        histInfoCounterTVX -> SetBinContent(iRun+1, vecInfoCounterTVX.at(iRun));
        histInfoCounterTVX -> GetXaxis() -> SetBinLabel(iRun+1, vecRunList.at(iRun).c_str());

        histInfoCounterScalTrig -> SetBinContent(iRun+1, vecInfoCounterScalTrig.at(iRun));
        histInfoCounterScalTrig -> GetXaxis() -> SetBinLabel(iRun+1, vecRunList.at(iRun).c_str());

        histInfoCounterSelTrig -> SetBinContent(iRun+1, vecInfoCounterSelTrig.at(iRun));
        histInfoCounterSelTrig -> GetXaxis() -> SetBinLabel(iRun+1, vecRunList.at(iRun).c_str());

        histSelCounterSelTOI -> SetBinContent(iRun+1, vecSelCounterSelTOI.at(iRun));
        histSelCounterSelTOI -> GetXaxis() -> SetBinLabel(iRun+1, vecRunList.at(iRun).c_str());

        histSelCounterAnalysedTrig -> SetBinContent(iRun+1, vecSelCounterAnalysedTrig.at(iRun));
        histSelCounterAnalysedTrig -> GetXaxis() -> SetBinLabel(iRun+1, vecRunList.at(iRun).c_str());
    }

    TFile *fOut = new TFile(Form("run_lists/%s/%s_%s_trigger_summary.root", year.c_str(), period.c_str(), triggerMask.c_str()), "RECREATE");
    histInfoCounterTVX -> Write();
    histInfoCounterScalTrig -> Write();
    histInfoCounterSelTrig -> Write();
    histSelCounterSelTOI -> Write();
    histSelCounterAnalysedTrig -> Write();
    fOut -> Close();
}

////////////////////////////////////////////////////////////////////////////////
void RetrieveTriggerInfo(TString dirName = "path/to/file", bool fromAlien = true, string triggerMask = "fDiMuon", double counters[5] = 0) {

    TString fInName;
    double infoTVX = 0;
    double infoScalTrig = 0;
    double infoSelTrig = 0;
    double selTOI = 0;
    double analysedTriggers = 0;

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

        for (auto const& dirKey : *fIn -> GetListOfKeys()) {
            if (TString(dirKey -> GetName()).Contains("table-maker")) {
                TList *list = (TList*) fIn -> Get("table-maker/Statistics");
                TH2D *histZorroInfo = (TH2D*) list -> FindObject("ZorroInfo");
                TH2D *histZorroSel = (TH2D*) list -> FindObject("ZorroSel");
    
                //std::cout << "inspected TVX: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin("inspected TVX")) << std::endl;
                //std::cout << triggerMask << " scalers: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str()))) << std::endl;
                //std::cout << triggerMask << " fSingleMuLow selections: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str()))) << std::endl;
                //std::cout << triggerMask << " fSingleMuLow scalers (SEL): " << histZorroSel -> GetBinContent(1, histZorroSel -> GetYaxis() -> FindBin(Form("%s", triggerMask.c_str()))) << std::endl;
    
                counters[0] = histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin("inspected TVX"));
                counters[1] = histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str())));
                counters[2] = histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str())));
                counters[3] = histZorroSel -> GetBinContent(1, histZorroSel -> GetYaxis() -> FindBin(Form("%s", triggerMask.c_str())));
                counters[4] = histZorroSel -> GetBinContent(1, histZorroSel -> GetYaxis() -> FindBin(Form("%s AnalysedTriggers", triggerMask.c_str())));

                infoTVX += histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin("inspected TVX"));
                infoScalTrig += histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str())));
                infoSelTrig += histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str())));
                selTOI += histZorroSel -> GetBinContent(1, histZorroSel -> GetYaxis() -> FindBin(Form("%s", triggerMask.c_str())));
                analysedTriggers += histZorroSel -> GetBinContent(1, histZorroSel -> GetYaxis() -> FindBin(Form("%s AnalysedTriggers", triggerMask.c_str())));
            }
        }

        counters[0] = infoTVX;
        counters[1] = infoScalTrig;
        counters[2] = infoSelTrig;
        counters[3] = selTOI;
        counters[4] = analysedTriggers;
    }
}