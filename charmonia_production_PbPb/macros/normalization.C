#include <TFile.h>
#include <TH1.h>
#include <TGrid.h>
#include <TEnv.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>
#include <thread>
#include <chrono>

#include <TDirectoryFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TLegend.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TH2.h>
#include <TGridResult.h>
#include <TSystem.h>

void RetrieveTriggerInfo(TString , bool , string, double [11]);

// xSecTVX cross section
//const double xSecTVX = 50.3e9; // pb-1, pp@5.36TeV, taken from https://alice-notes.web.cern.ch/system/files/notes/analysis/1671/2025-09-04-Analysis_Note_Pi0_OO.pdf, page 14
const double xSecTVX = 1.13e12; // pb-1, OO@5.36TeV, taken from https://alice-notes.web.cern.ch/system/files/notes/analysis/1671/2025-09-04-Analysis_Note_Pi0_OO.pdf, page 14

const int kMaxFileAttempts = 1;
const int kWaitSeconds = 5;


bool waitForAlienConnection(int maxRetries = 5, int waitSeconds = 5) {
    for (int i = 0; i < maxRetries; i++) {
        try {
            if (!gGrid) TGrid::Connect("alien://");
            TGridResult *res = gGrid->Command("whoami");
            if (res) { delete res; return true; }
        } catch (...) {}
        std::cerr << "[WARNING] Alien connection not ready, retrying in " << waitSeconds << " s...\n";
        TGrid::Connect("alien://");
        std::this_thread::sleep_for(std::chrono::seconds(waitSeconds));
    }
    std::cerr << "[ERROR] Could not establish Alien connection.\n";
    return false;
}

void get_normalization_from_multiple_files(std::string year = "2023",
                                           std::string period = "LHC23_PbPb_pass4",
                                           std::string subPeriod = "None",
                                           std::string triggerMask = "minBias",
                                           std::string assocType = "std_assoc",
                                           std::string location = "alien") {

    gEnv->SetValue("XrdSecGSISRVNAMES", "");
    gEnv->SetValue("XrdSecGSISRVNAMESCANONICAL", "");
    gEnv->SetValue("XrdSecGSIUSERPROMPT", 0);

    ifstream fInInputPath(Form("%s_input_path_short.txt", location.c_str()));
    if (!fInInputPath.is_open()) {
        std::cout << "Error opening " << Form("%s_input_path_short.txt", location.c_str()) << std::endl;
        return;
    }

    ifstream fInRunList(Form("%s_run_list_short.txt", location.c_str()));
    if (!fInRunList.is_open()) {
        std::cout << "Error opening " << Form("%s_run_list_short.txt", location.c_str()) << std::endl;
        return;
    }

    string subRunListName;
    if (subPeriod == "None") {
        subRunListName = Form("%s_run_list_short.txt", location.c_str());
    } else {
        subRunListName = Form("%s.txt", subPeriod.c_str());
    }

    ifstream fInSubRunList(subRunListName);
    if (!fInSubRunList.is_open()) {
        std::cout << "Error opening the file!" << std::endl;
        return;
    }

    if (!fInInputPath.is_open() || !fInRunList.is_open() || !fInSubRunList.is_open()) {
        std::cerr << "[ERROR] Cannot open input files.\n"; return;
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
    while (getline(fInRunList, alienRunListName)) {
        vecRunList.push_back(alienRunListName);
    }

    string alienSubRunListName;
    while (getline(fInSubRunList, alienSubRunListName)) {
        vecSubRunList.push_back(alienSubRunListName);
    }

    string fOutName;
    if (subPeriod == "None") {
        fOutName = Form("%s_%s_%s_trigger_summary.root", year.c_str(), period.c_str(), triggerMask.c_str(), assocType.c_str());
    } else {
        fOutName = Form("%s_%s_%s_trigger_summary.root", year.c_str(), subPeriod.c_str(), triggerMask.c_str(), assocType.c_str());
    }
    
    bool resuming = std::filesystem::exists(fOutName);
    TFile *fOut = TFile::Open(fOutName.c_str(), resuming ? "UPDATE" : "RECREATE");
    if (!fOut || fOut->IsZombie()) { std::cerr << "[ERROR] Cannot open/create output ROOT file.\n"; return; }

    TH1D *histZorroInfoCounterTVX = (TH1D*)fOut->Get("histZorroInfoCounterTVX");
    if (!histZorroInfoCounterTVX) histZorroInfoCounterTVX = new TH1D("histZorroInfoCounterTVX", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histZorroInfoCounterScalTrig = (TH1D*)fOut->Get("histZorroInfoCounterScalTrig");
    if (!histZorroInfoCounterScalTrig) histZorroInfoCounterScalTrig = new TH1D("histZorroInfoCounterScalTrig", "", vecSubRunList.size(), 0, vecSubRunList.size());
        TH1D *histZorroInfoCounterSelTrig = (TH1D*)fOut->Get("histZorroInfoCounterSelTrig");
    if (!histZorroInfoCounterSelTrig) histZorroInfoCounterSelTrig = new TH1D("histZorroInfoCounterSelTrig", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histZorroSelCounterSelTOI = (TH1D*)fOut->Get("histZorroSelCounterSelTOI");
    if (!histZorroSelCounterSelTOI) histZorroSelCounterSelTOI = new TH1D("histZorroSelCounterSelTOI", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histZorroSelCounterAnalysedTrig = (TH1D*)fOut->Get("histZorroSelCounterAnalysedTrig");
    if (!histZorroSelCounterAnalysedTrig) histZorroSelCounterAnalysedTrig = new TH1D("histZorroSelCounterAnalysedTrig", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histEvSelCollisionBeforeFiltering = (TH1D*)fOut->Get("histEvSelCollisionBeforeFiltering");
    if (!histEvSelCollisionBeforeFiltering) histEvSelCollisionBeforeFiltering = new TH1D("histEvSelCollisionBeforeFiltering", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histEvSelCollisionBeforeCuts = (TH1D*)fOut->Get("histEvSelCollisionBeforeCuts");
    if (!histEvSelCollisionBeforeCuts) histEvSelCollisionBeforeCuts = new TH1D("histEvSelCollisionBeforeCuts", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histEvSelCollisionAfterCuts = (TH1D*)fOut->Get("histEvSelCollisionAfterCuts");
    if (!histEvSelCollisionAfterCuts) histEvSelCollisionAfterCuts = new TH1D("histEvSelCollisionAfterCuts", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histBcSelCounterTVX = (TH1D*)fOut->Get("histBcSelCounterTVX");
    if (!histBcSelCounterTVX) histBcSelCounterTVX = new TH1D("histBcSelCounterTVX", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histBcSelCounterTVXafterBCcuts = (TH1D*)fOut->Get("histBcSelCounterTVXafterBCcuts");
    if (!histBcSelCounterTVXafterBCcuts) histBcSelCounterTVXafterBCcuts = new TH1D("histBcSelCounterTVXafterBCcuts", "", vecSubRunList.size(), 0, vecSubRunList.size());
    TH1D *histEvSelColCounterTVX = (TH1D*)fOut->Get("histEvSelColCounterTVX");
    if (!histEvSelColCounterTVX) histEvSelColCounterTVX = new TH1D("histEvSelColCounterTVX", "", vecSubRunList.size(), 0, vecSubRunList.size());

    auto savePartial = [&](const char *reason){
        std::cerr << "[INFO] Saving ROOT file: " << reason << std::endl;
        fOut->cd();
        histZorroInfoCounterTVX -> Write("",TObject::kOverwrite);
        histZorroInfoCounterScalTrig -> Write("",TObject::kOverwrite);
        histZorroInfoCounterSelTrig -> Write("",TObject::kOverwrite);
        histZorroSelCounterSelTOI -> Write("",TObject::kOverwrite);
        histZorroSelCounterAnalysedTrig -> Write("",TObject::kOverwrite);
        histEvSelCollisionBeforeFiltering -> Write("",TObject::kOverwrite);
        histEvSelCollisionBeforeCuts -> Write("",TObject::kOverwrite);
        histEvSelCollisionAfterCuts -> Write("",TObject::kOverwrite);
        histBcSelCounterTVX -> Write("",TObject::kOverwrite);
        histBcSelCounterTVXafterBCcuts -> Write("",TObject::kOverwrite);
        histEvSelColCounterTVX -> Write("",TObject::kOverwrite);
        fOut->Flush();
    };

    string pathName;
    string runNumber;
    int runCounter = 0;
    int subRunCounter = 0;

    while (getline(fInInputPath, pathName)) {
        runNumber = vecRunList.at(runCounter).c_str();

        std::cout << "Processing Run " << runNumber << " ..." << std::endl;
        if (std::find(vecSubRunList.begin(), vecSubRunList.end(), runNumber) == vecSubRunList.end()) {
            runCounter++;
            continue;
        }

        if(histZorroInfoCounterTVX->GetBinContent(subRunCounter+1)!=0){ runCounter++; subRunCounter++; continue; }

        std::cout << "[INFO] Processing Run " << runNumber << std::endl;

        if(location=="alien" && !waitForAlienConnection()) {
            savePartial("Lost Alien connection, skipping run");
            runCounter++; subRunCounter++; continue;
        }

        double counters[11] = {0};
        RetrieveTriggerInfo(pathName.c_str(), location=="alien", triggerMask, counters);

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

        savePartial(("Checkpoint run "+runNumber).c_str());

        runCounter++; subRunCounter++;
    }

    savePartial("Normal completion");
    fOut->Close();
    std::cout << "[SUCCESS] Completed and saved results in " << fOutName << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_normalization_from_single_file(string year = "2024", string period = "LHC24aq_pass1_small", string triggerMask = "minBias", string assocType = "std_assoc") {
    //string fInName = Form("data/%s/%s_%s_%s.root", year.c_str(), period.c_str(), triggerMask.c_str(), assocType.c_str());
    //string fInName = Form("/Users/lucamicheletti/cernbox/JPSI/Jpsi_pp_ref_cross_section_run3/data/%s/%s/AnalysisResults_%s.root", year.c_str(), period.c_str(), assocType.c_str());
    string fInName = Form("/Users/lucamicheletti/cernbox/JPSI/Jpsi_pO_OO_NeNe_run3/data/%s/%s/AnalysisResults_%s.root", year.c_str(), period.c_str(), assocType.c_str());
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


        if (TString(dirKey -> GetName()).Contains("eventselection-run3")) {
            TH1D *histCounterTVX = (TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVX");
            TH1D *histCounterTVXafterBCcuts = (TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVXafterBCcuts");

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

        if (TString(dirKey -> GetName()).Contains("eventselection-run3")) {
            TH1D *histCounterTVX = (TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVX");
            TH1D *histCounterTVXafterBCcuts = (TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVXafterBCcuts");
            TH1D *histColCounterTVX = (TH1D*) fIn -> Get("eventselection-run3/eventselection/hColCounterTVX");

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
            histEvSelColCounterTVX = new TH1D("histEvSelColCounterTVX", "", vecRunList.size(), 0, vecRunList.size());

            for (int iRun = 0;iRun < int (vecRunList.size());iRun++) {
                string runNumber = vecRunList.at(iRun).c_str();
                histBcSelCounterTVX -> SetBinContent(iRun+1, histCounterTVX -> GetBinContent(iRun+1));
                histBcSelCounterTVX -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());

                histBcSelCounterTVXafterBCcuts -> SetBinContent(iRun+1, histCounterTVXafterBCcuts -> GetBinContent(iRun+1));
                histBcSelCounterTVXafterBCcuts -> GetXaxis() -> SetBinLabel(iRun+1, runNumber.c_str());

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
void luminosity(string year = "2024", string period = "LHC24_ppref_pass1", string triggerMask = "minBias", string assocType = "std_assoc") {
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

    double bcCounterEfficiency = bcCounterTVXafterBCcuts / bcCounterTVX;
    double nEvtsBcSel, nEvtsBcSelAfterCuts = -999;
    double luminosityBcSel, luminosityBcSelAfterCuts = -999;
    if (triggerMask.find("minBias") != std::string::npos) {
        nEvtsBcSel = (bcCounterTVX * (collsAfterCuts / collsBeforeCuts));
        nEvtsBcSelAfterCuts = (bcCounterTVXafterBCcuts * (collsAfterCuts / collsBeforeCuts));
        luminosityBcSel = nEvtsBcSel / xSecTVX;
        luminosityBcSelAfterCuts = nEvtsBcSelAfterCuts / xSecTVX;

        double luminosityBcSelMinBias = (bcCounterTVX * (collsAfterCuts / collsBeforeCuts)) / xSecTVX;
        double luminosityBcSelAfterCutsMinBias = (bcCounterTVXafterBCcuts * (collsAfterCuts / collsBeforeCuts)) / xSecTVX;
        double luminosityEvSelMinBias = (evColCounterTVX * (collsAfterCuts / collsBeforeCuts)) / xSecTVX;
        std::cout << "BC selection efficiency = " << bcCounterTVXafterBCcuts / bcCounterTVX << std::endl;
        std::cout << "N. events            = " << (bcCounterTVX * (collsAfterCuts / collsBeforeCuts)) << std::endl;
        std::cout << "Colls after cuts     = " << collsAfterCuts << std::endl;
        std::cout << "luminosity Min. Bias [bc-selection-task]                  = " << luminosityBcSelMinBias << " pb-1" << std::endl;
        std::cout << "luminosity Min. Bias after BC cuts [bc-selection-task]    = " << luminosityBcSelAfterCutsMinBias << " pb-1" << std::endl;
        std::cout << "luminosity Min. Bias [event-selection-task]               = " << luminosityEvSelMinBias << " pb-1" << std::endl;
    } else {
        double luminosityZorroInfo = (inspectedTVX * (selectionsZorroInfo / scalers)) / xSecTVX;
        double luminosityZorroSel = (inspectedTVX * (selectionsZorroSel / scalers)) / xSecTVX;
        double luminosityAnalysedTrig = (inspectedTVX * (selectionsAnalysedTrig / scalers)) / xSecTVX;
        std::cout << "BC selection efficiency = " << bcCounterTVXafterBCcuts / bcCounterTVX << std::endl;
        std::cout << "N. events TOI           = " << (inspectedTVX * (selectionsZorroSel / scalers)) << std::endl;
        std::cout << "N. events AnalysedTrig  = " << (inspectedTVX * (selectionsAnalysedTrig / scalers)) << std::endl;
        std::cout << "luminosityZorroInfo    = " << luminosityZorroInfo << " pb-1" << std::endl;
        std::cout << "luminosity TOI          = " << luminosityZorroSel << " pb-1" << std::endl;
        std::cout << "luminosity AnalysedTrig = " << luminosityAnalysedTrig << " pb-1" << std::endl;

        nEvtsBcSel = (inspectedTVX * (selectionsAnalysedTrig / scalers));
        nEvtsBcSelAfterCuts = (inspectedTVX * (selectionsZorroSel / scalers));
        luminosityBcSel = nEvtsBcSel / xSecTVX;
        luminosityBcSelAfterCuts = nEvtsBcSelAfterCuts / xSecTVX;
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
void check_normalization(string fInName = "LHC25ae_minBias_std_assoc_luminosity.root", string fOutName = "luminosity_jpsi_LHC25ae_pass2_minBias.root") {
    TFile *fInLumiMinBias2024StdAssoc = TFile::Open(Form("data/2025/pass2/%s", fInName.c_str()));
    TH1D *histLumiMinBias2024StdAssoc = (TH1D*) fInLumiMinBias2024StdAssoc -> Get("histLumiSummary");
    double nEvtsMinBias2024StdAssoc = histLumiMinBias2024StdAssoc -> GetBinContent(histLumiMinBias2024StdAssoc -> GetXaxis() -> FindBin("nEvtsBcSel"));
    double lumiMinBias2024StdAssoc = histLumiMinBias2024StdAssoc -> GetBinContent(histLumiMinBias2024StdAssoc -> GetXaxis() -> FindBin("luminosityBcSel"));

    TFile *fInMinBias2024StdAssoc = TFile::Open("/Users/lucamicheletti/cernbox/JPSI/Jpsi_pp_ref_cross_section_run3/data/2024/LHC24_ppref_pass1/AnalysisResults_std_assoc.root");
    TList *listMinBias2024StdAssoc = (TList*) fInMinBias2024StdAssoc -> Get("analysis-same-event-pairing/output");
    TList *listMinBias2024StdAssocSE = (TList*) listMinBias2024StdAssoc -> FindObject("PairsMuonSEPM_matchedMchMid");
    THnSparseD *histSparseMinBias2024StdAssoc = (THnSparseD*) listMinBias2024StdAssocSE -> FindObject("Mass_Pt_Rapidity");
    histSparseMinBias2024StdAssoc -> GetAxis(1) -> SetRangeUser(0., 20.);
    histSparseMinBias2024StdAssoc -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH1D *histProjIntMinBias2024StdAssoc = (TH1D*) histSparseMinBias2024StdAssoc -> Projection(0, "Proj_Mass_Pt_Rapidity");
    double countsMinBias2024StdAssoc = histProjIntMinBias2024StdAssoc -> Integral();
    histProjIntMinBias2024StdAssoc -> Scale(1. / (lumiMinBias2024StdAssoc));

    // Added for J/psi - D0 associated production
    TH1D *histLuminosityMinBias2024StdAssoc = new TH1D("histLumi", "; ; Luminosity (pb-1)", 1, 0, 1);
    histLuminosityMinBias2024StdAssoc -> SetBinContent(1, lumiMinBias2024StdAssoc);

    TH1D *histNeventMinBias2024StdAssoc = new TH1D("histNevMinBias", "; ; N. Events min bias", 1, 0, 1);
    histNeventMinBias2024StdAssoc -> SetBinContent(1, nEvtsMinBias2024StdAssoc);

    TFile *fOutLumiMinBias2024StdAssoc = new TFile(Form("data/%s", fOutName.c_str()), "RECREATE");
    histLuminosityMinBias2024StdAssoc -> Write();
    histNeventMinBias2024StdAssoc -> Write();
    fOutLumiMinBias2024StdAssoc -> Close();

    std::cout << "---------------------------" << std::endl;
    std::cout << "[Min. Bias, std. assoc.]" << std::endl;
    std::cout << "N. evts.         = " << nEvtsMinBias2024StdAssoc << std::endl;
    std::cout << "Luminosity       = " << lumiMinBias2024StdAssoc << std::endl;
    std::cout << "N. evts. corr.   = " << nEvtsMinBias2024StdAssoc << std::endl;
    std::cout << "Luminosity corr. = " << lumiMinBias2024StdAssoc << std::endl;
    std::cout << "---------------------------" << std::endl;


}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void RetrieveTriggerInfo(TString dirName = "path/to/file", bool fromAlien = true, std::string triggerMask = "fDiMuon", double counters[11] = nullptr) {
    if (!counters) return;

    for (int i = 0; i < 11; i++) counters[i] = 0;

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
        try {
            if (!gGrid || !gGrid->IsConnected()) TGrid::Connect("alien://");
            if (!dirName.Contains("alien://")) dirName = "alien://" + dirName + "/AOD";
        } catch (...) {
            std::cerr << "[ERROR] Cannot connect to Alien, will retry per file.\n";
        }
    }

    for (int iDir = 1; iDir < 7000; iDir++) {
        TString fInName;
        if (iDir < 10) fInName = dirName + Form("/00%i/AnalysisResults.root", iDir);
        else if (iDir < 100) fInName = dirName + Form("/0%i/AnalysisResults.root", iDir);
        else fInName = dirName + Form("/%i/AnalysisResults.root", iDir);

        TFile* fIn = nullptr;
        int attempt = 0;

        while (attempt < kMaxFileAttempts) {
            try {
                if (fromAlien && (!gGrid || !gGrid->IsConnected())) TGrid::Connect("alien://");

                fIn = TFile::Open(fInName.Data());
                if (fIn && !fIn->IsZombie()) break; 
                else {
                    std::cerr << "[WARNING] File " << fInName << " not available, retrying...\n";
                    if (fIn) { fIn->Close(); fIn = nullptr; }
                }
            } catch (...) {
                std::cerr << "[ERROR] Exception opening file " << fInName << ", retrying...\n";
            }
            attempt++;
            if (attempt < kMaxFileAttempts) std::this_thread::sleep_for(std::chrono::seconds(kWaitSeconds));
        }

        if (!fIn || fIn->IsZombie()) {
            std::cerr << "[ERROR] Could not open file " << fInName << " after " << kMaxFileAttempts << " attempts. Skipping this file.\n";
            return; 
        }

        std::cout << "[INFO] File exists: " << fInName << std::endl;


        for (auto dirKey : *fIn->GetListOfKeys()) {
            TString keyName = dirKey->GetName();

            if (TString(dirKey -> GetName()).Contains("table-maker")) {
                TList *list = (TList*) fIn -> Get("table-maker/Statistics");
                TH2D *histZorroInfo = (TH2D*) list -> FindObject("ZorroInfo");
                TH2D *histZorroSel = (TH2D*) list -> FindObject("ZorroSel");
                TH2D *histEventStats = (TH2D*) list -> FindObject("EventStats");
    
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

            if (TString(dirKey -> GetName()).Contains("eventselection-run3")) {
                TH1D *histColCounterTVX = (TH1D*) fIn -> Get("eventselection-run3/eventselection/hColCounterTVX");
                TH1D *histCounterTVX = (TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVX");
                TH1D *histCounterTVXafterBCcuts = (TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVXafterBCcuts");

                colCounterTVX += histColCounterTVX -> GetBinContent(1);
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
        counters[10] = colCounterTVX;

        fIn->Close();
    }
}