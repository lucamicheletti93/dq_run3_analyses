#include <fstream>
#include <iostream>
#include <string>
#include <TSystem.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TProfile.h"

using namespace std;

TH1D* ProjectTHnSparse(THnSparseD *, double , double , double , double , double , double, double , double);
TH1D* ProjectTHnSparsePt(THnSparseD *, double , double , double , double , double , double, double , double);
TH1D* ProjectTHnSparseRap(THnSparseD *, double , double , double , double , double , double, double , double);
TH1D* ProjectTH2(TH2D *, string, double , double);

void projections(){
    //----------------------------------------------------------------------------------------------//
                                        //centrality cuts and muon cuts
    vector <Double_t> centralityCutsMin = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0};
    vector <Double_t> centralityCutsMax = {10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0};
    vector <Double_t> centralityCutsMin2 = {0.0,20.0,40.0,60.0};
    vector <Double_t> centralityCutsMax2 = {20.0,40.0,60.0,90.0};
    int sizeCentrCuts = centralityCutsMin.size();
    int sizeCentrCuts2 = centralityCutsMin2.size();
    vector <Double_t> PtBinsMin = {0.0,2.0,4.0,5.0};
    vector <Double_t> PtBinsMax = {2.0,4.0,5.0,12.0};
    vector <Double_t> PtBinsMin2 = {4.0,6.0};
    vector <Double_t> PtBinsMax2 = {6.0,12.0};
    int sizePtBins = PtBinsMin.size();
    int sizePtBins2 = PtBinsMin2.size();
    //string muonCut = "muonLowPt10SigmaPDCA";
    //string muonCut = "muonLowPt210SigmaPDCA";
    string muonCut = "muonLowPt510SigmaPDCA";
    //string muonCut = "matchedMchMid";

    //string histName = "Mass_Pt_centrFT0C_V2";
    string histName = "Mass_Pt_Cent_Rap";
    string histNamePt = "Mass_Pt";
    string histNameRap = "Mass_Rapidity";
    //string histName = "Mass_CentFT0C";
    //----------------------------------------------------------------------------------------------//
                                        //open input file

    /* string pathToFiles = "/Users/saragaretti/macros/analysis_J_psi/dq_fitter-main_old/analysis";
    TFile* analysisResults = TFile::Open(Form("%s/AnalysisResultsMC.root",pathToFiles.c_str())); */

    //string pathToFiles = "/Users/saragaretti/macros/analysis_J_psi/dq_fitter-main_old/data/LHC23_full";
    //string pathToFiles = "/Users/saragaretti/cernbox/JPsi_Psi2S_ratio/Train_295284/data/LHC23_pass4";
    //string pathToFiles = "/Users/saragaretti/cernbox/JPsi_Psi2S_ratio/Train_314824/data/LHC23_pass4_pTcut_1GeV";
    string pathToFiles ="/Users/saragaretti/cernbox/JPsi_Psi2S_ratio/Train_321538/data/LHC23_pass4_pTcut_1GeV_new";

    ifstream fRunList(Form("%s/run_list_K2.txt", pathToFiles.c_str()));
    int runNumber;
    vector<int> runList;
    runList.reserve(200);
    while(!fRunList.eof()) {
        fRunList >> runNumber;
        runList.push_back(runNumber);
    }

    int runCounter = 0;
    for (auto const& run : runList) {
        TFile *fIn = new TFile(Form("%s/%d/AnalysisResults.root", pathToFiles.c_str(), run), "READ");
        if (fIn -> IsZombie()) {
            std::cout << "Run = " << run << " -> ISSUE!" << std::endl;
            continue;
        }
        std::cout << "Run = " << run << " -> OK!" << std::endl;
        runCounter++;
        std::cout << "Run counter = " << runCounter << std::endl;

        //----------------------------------------------------------------------------------------------//
                                            //create output file

        TFile *fOut = new TFile(Form("%s/%d/Histograms_%s.root", pathToFiles.c_str(), run, muonCut.c_str()), "RECREATE");
        //TFile *fOut = new TFile(Form("%s/Projections_MC_sparse.root", pathToFiles.c_str()), "RECREATE");

        //---------------------------------------------------------------------------------------------//
                                            //Same event pairing

            TList *hListSE = (TList*)fIn->Get("analysis-same-event-pairing/output");
            if (!hListSE) {
                std::cout << "Run = " << run << " -> ISSUE WITH SAME EVENT!" << std::endl;
                continue;
            }
            
            TList *hListSEPM = (TList*)hListSE->FindObject(Form("PairsMuonSEPM_%s", muonCut.c_str()));
            TList *hListSEMM = (TList*)hListSE->FindObject(Form("PairsMuonSEMM_%s", muonCut.c_str()));
            TList *hListSEPP = (TList*)hListSE->FindObject(Form("PairsMuonSEPP_%s", muonCut.c_str()));
            

        //---------------------------------------------------------------------------------------------//
                                                //Event mixing
            
            TList *hListME = (TList*)fIn->Get("analysis-event-mixing/output");
            if (!hListME) {
                std::cout << "Run = " << run << " -> ISSUE WITH EVENT MIXING!" << std::endl;
                continue;
            }
            TList *hListMEPM = (TList*)hListME->FindObject(Form("PairsMuonMEPM_%s", muonCut.c_str()));
            TList *hListMEMM = (TList*)hListME->FindObject(Form("PairsMuonMEMM_%s", muonCut.c_str()));
            TList *hListMEPP = (TList*)hListME->FindObject(Form("PairsMuonMEPP_%s", muonCut.c_str()));


        if (histName == "Mass_Pt_Cent_Rap") {
            //---------------------------------------------------------------------------------------------//
                                            //Same event pairing

            THnSparseD *histSparseSEPM = (THnSparseD*) hListSEPM -> FindObject(histName.c_str());
            if (!hListSEPM) {
                std::cout << "Run = " << run << " -> ISSUE WITH SAME EVENT!" << std::endl;
                continue;
            }
            THnSparseD *histSparseSEPP = (THnSparseD*) hListSEMM -> FindObject(histName.c_str());
            THnSparseD *histSparseSEMM = (THnSparseD*) hListSEPP -> FindObject(histName.c_str());

            TH1D *histSEPMProjInt = (TH1D*) ProjectTHnSparse(histSparseSEPM, 2, 5, 0, 30, 0, 90, 2.5, 4);
            TH1D *histSEPPProjInt = (TH1D*) ProjectTHnSparse(histSparseSEPP, 2, 5, 0, 30, 0, 90, 2.5, 4);
            TH1D *histSEMMProjInt = (TH1D*) ProjectTHnSparse(histSparseSEMM, 2, 5, 0, 30, 0, 90, 2.5, 4);
            //pT
            TH1D *histSEPMProjPt = (TH1D*) ProjectTHnSparsePt(histSparseSEPM, 2, 5, 0, 30, 0, 90, 2.5, 4);
            TH1D *histSEPPProjPt = (TH1D*) ProjectTHnSparsePt(histSparseSEPP, 2, 5, 0, 30, 0, 90, 2.5, 4);
            TH1D *histSEMMProjPt = (TH1D*) ProjectTHnSparsePt(histSparseSEMM, 2, 5, 0, 30, 0, 90, 2.5, 4);
            //rapidity
            TH1D *histSEPMProjRap = (TH1D*) ProjectTHnSparseRap(histSparseSEPM, 2, 5, 0, 30, 0, 90, 2.5, 4);
            TH1D *histSEPPProjRap = (TH1D*) ProjectTHnSparseRap(histSparseSEPP, 2, 5, 0, 30, 0, 90, 2.5, 4);
            TH1D *histSEMMProjRap = (TH1D*) ProjectTHnSparseRap(histSparseSEMM, 2, 5, 0, 30, 0, 90, 2.5, 4);

            //---------------------------------------------------------------------------------------------//
                                                //Event mixing

            THnSparseD *histSparseMEPM = (THnSparseD*) hListMEPM -> FindObject(histName.c_str());
            if (!hListMEPM) {
                std::cout << "Run = " << run << " -> ISSUE WITH SAME EVENT!" << std::endl;
                continue;
            }
            THnSparseD *histSparseMEPP = (THnSparseD*) hListMEMM -> FindObject(histName.c_str());
            THnSparseD *histSparseMEMM = (THnSparseD*) hListMEPP -> FindObject(histName.c_str());
            
            TH1D *histMEPMProjInt = (TH1D*) ProjectTHnSparse(histSparseMEPM, 2, 5, 0, 30, 0, 90, 2.5, 4);
            TH1D *histMEPPProjInt = (TH1D*) ProjectTHnSparse(histSparseMEPP, 2, 5, 0, 30, 0, 90, 2.5, 4);
            TH1D *histMEMMProjInt = (TH1D*) ProjectTHnSparse(histSparseMEMM, 2, 5, 0, 30, 0, 90, 2.5, 4);

            //pT
            TH1D *histMEPMProjPt = (TH1D*) ProjectTHnSparsePt(histSparseMEPM, 2, 5, 0, 30, 0, 90, 2.5, 4);
            TH1D *histMEPPProjPt = (TH1D*) ProjectTHnSparsePt(histSparseMEPP, 2, 5, 0, 30, 0, 90, 2.5, 4);
            TH1D *histMEMMProjPt = (TH1D*) ProjectTHnSparsePt(histSparseMEMM, 2, 5, 0, 30, 0, 90, 2.5, 4);
            //rapidity
            TH1D *histMEPMProjRap = (TH1D*) ProjectTHnSparseRap(histSparseMEPM, 2, 5, 0, 30, 0, 90, 2.5, 4);
            TH1D *histMEPPProjRap = (TH1D*) ProjectTHnSparseRap(histSparseMEPP, 2, 5, 0, 30, 0, 90, 2.5, 4);
            TH1D *histMEMMProjRap = (TH1D*) ProjectTHnSparseRap(histSparseMEMM, 2, 5, 0, 30, 0, 90, 2.5, 4);

            //---------------------------------------------------------------------------------------------//
                                //Centrality cuts on projected histograms

            for (int iCentr = 0;iCentr < sizeCentrCuts;iCentr++) {
                TH1D *histSEPMProjCentr = (TH1D*) ProjectTHnSparse(histSparseSEPM, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histSEPPProjCentr = (TH1D*) ProjectTHnSparse(histSparseSEPP, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histSEMMProjCentr = (TH1D*) ProjectTHnSparse(histSparseSEMM, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histMEPMProjCentr = (TH1D*) ProjectTHnSparse(histSparseMEPM, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histMEPPProjCentr = (TH1D*) ProjectTHnSparse(histSparseMEPP, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histMEMMProjCentr = (TH1D*) ProjectTHnSparse(histSparseMEMM, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);

                histSEPMProjCentr -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEPM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histSEPPProjCentr -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEPP", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histSEMMProjCentr -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEMM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histMEPMProjCentr -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEPM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histMEPPProjCentr -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEPP", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histMEMMProjCentr -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEMM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                //pT
                TH1D *histSEPMProjPtCentr = (TH1D*) ProjectTHnSparsePt(histSparseSEPM, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histSEPPProjPtCentr = (TH1D*) ProjectTHnSparsePt(histSparseSEPP, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histSEMMProjPtCentr = (TH1D*) ProjectTHnSparsePt(histSparseSEMM, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histMEPMProjPtCentr = (TH1D*) ProjectTHnSparsePt(histSparseMEPM, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histMEPPProjPtCentr = (TH1D*) ProjectTHnSparsePt(histSparseMEPP, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histMEMMProjPtCentr = (TH1D*) ProjectTHnSparsePt(histSparseMEMM, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);

                histSEPMProjPtCentr -> Write(Form("Mass_Pt_%.0f_%.0f_SEPM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histSEPPProjPtCentr -> Write(Form("Mass_Pt_%.0f_%.0f_SEPP", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histSEMMProjPtCentr -> Write(Form("Mass_Pt_%.0f_%.0f_SEMM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histMEPMProjPtCentr -> Write(Form("Mass_Pt_%.0f_%.0f_MEPM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histMEPPProjPtCentr -> Write(Form("Mass_Pt_%.0f_%.0f_MEPP", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histMEMMProjPtCentr -> Write(Form("Mass_Pt_%.0f_%.0f_MEMM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                //Rapidity
                TH1D *histSEPMProjRapCentr = (TH1D*) ProjectTHnSparseRap(histSparseSEPM, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histSEPPProjRapCentr = (TH1D*) ProjectTHnSparseRap(histSparseSEPP, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histSEMMProjRapCentr = (TH1D*) ProjectTHnSparseRap(histSparseSEMM, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histMEPMProjRapCentr = (TH1D*) ProjectTHnSparseRap(histSparseMEPM, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histMEPPProjRapCentr = (TH1D*) ProjectTHnSparseRap(histSparseMEPP, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);
                TH1D *histMEMMProjRapCentr = (TH1D*) ProjectTHnSparseRap(histSparseMEMM, 2, 5, 0, 30, centralityCutsMin[iCentr], centralityCutsMax[iCentr], 2.5, 4);

                histSEPMProjRapCentr -> Write(Form("Mass_Rap_%.0f_%.0f_SEPM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histSEPPProjRapCentr -> Write(Form("Mass_Rap_%.0f_%.0f_SEPP", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histSEMMProjRapCentr -> Write(Form("Mass_Rap_%.0f_%.0f_SEMM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histMEPMProjRapCentr -> Write(Form("Mass_Rap_%.0f_%.0f_MEPM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histMEPPProjRapCentr -> Write(Form("Mass_Rap_%.0f_%.0f_MEPP", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histMEMMProjRapCentr -> Write(Form("Mass_Rap_%.0f_%.0f_MEMM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
            }

            for (int iCentr = 0;iCentr < sizeCentrCuts2;iCentr++) {
                TH1D *histSEPMProjCentr2 = (TH1D*) ProjectTHnSparse(histSparseSEPM, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histSEPPProjCentr2 = (TH1D*) ProjectTHnSparse(histSparseSEPP, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histSEMMProjCentr2 = (TH1D*) ProjectTHnSparse(histSparseSEMM, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histMEPMProjCentr2 = (TH1D*) ProjectTHnSparse(histSparseMEPM, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histMEPPProjCentr2 = (TH1D*) ProjectTHnSparse(histSparseMEPP, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histMEMMProjCentr2 = (TH1D*) ProjectTHnSparse(histSparseMEMM, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);

                histSEPMProjCentr2 -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEPM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histSEPPProjCentr2 -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEPP", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histSEMMProjCentr2 -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEMM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histMEPMProjCentr2 -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEPM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histMEPPProjCentr2 -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEPP", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histMEMMProjCentr2 -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEMM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                //pT
                TH1D *histSEPMProjPtCentr2 = (TH1D*) ProjectTHnSparsePt(histSparseSEPM, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histSEPPProjPtCentr2 = (TH1D*) ProjectTHnSparsePt(histSparseSEPP, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histSEMMProjPtCentr2 = (TH1D*) ProjectTHnSparsePt(histSparseSEMM, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histMEPMProjPtCentr2 = (TH1D*) ProjectTHnSparsePt(histSparseMEPM, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histMEPPProjPtCentr2 = (TH1D*) ProjectTHnSparsePt(histSparseMEPP, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histMEMMProjPtCentr2 = (TH1D*) ProjectTHnSparsePt(histSparseMEMM, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);

                histSEPMProjPtCentr2 -> Write(Form("Mass_Pt_%.0f_%.0f_SEPM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histSEPPProjPtCentr2 -> Write(Form("Mass_Pt_%.0f_%.0f_SEPP", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histSEMMProjPtCentr2 -> Write(Form("Mass_Pt_%.0f_%.0f_SEMM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histMEPMProjPtCentr2 -> Write(Form("Mass_Pt_%.0f_%.0f_MEPM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histMEPPProjPtCentr2 -> Write(Form("Mass_Pt_%.0f_%.0f_MEPP", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histMEMMProjPtCentr2 -> Write(Form("Mass_Pt_%.0f_%.0f_MEMM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                //Rapidity
                TH1D *histSEPMProjRapCentr2 = (TH1D*) ProjectTHnSparseRap(histSparseSEPM, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histSEPPProjRapCentr2 = (TH1D*) ProjectTHnSparseRap(histSparseSEPP, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histSEMMProjRapCentr2 = (TH1D*) ProjectTHnSparseRap(histSparseSEMM, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histMEPMProjRapCentr2 = (TH1D*) ProjectTHnSparseRap(histSparseMEPM, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histMEPPProjRapCentr2 = (TH1D*) ProjectTHnSparseRap(histSparseMEPP, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);
                TH1D *histMEMMProjRapCentr2 = (TH1D*) ProjectTHnSparseRap(histSparseMEMM, 2, 5, 0, 30, centralityCutsMin2[iCentr], centralityCutsMax2[iCentr], 2.5, 4);

                histSEPMProjRapCentr2 -> Write(Form("Mass_Rap_%.0f_%.0f_SEPM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histSEPPProjRapCentr2 -> Write(Form("Mass_Rap_%.0f_%.0f_SEPP", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histSEMMProjRapCentr2 -> Write(Form("Mass_Rap_%.0f_%.0f_SEMM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histMEPMProjRapCentr2 -> Write(Form("Mass_Rap_%.0f_%.0f_MEPM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histMEPPProjRapCentr2 -> Write(Form("Mass_Rap_%.0f_%.0f_MEPP", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histMEMMProjRapCentr2 -> Write(Form("Mass_Rap_%.0f_%.0f_MEMM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
            }
            //---------------------------------------------------------------------------------------------//
                                                //Pt cuts on projected histograms
            
            for (int iPt = 0;iPt < sizePtBins;iPt++) {
                TH1D *histSEPMProjPt = (TH1D*) ProjectTHnSparse(histSparseSEPM, 2, 5, PtBinsMin[iPt], PtBinsMax[iPt], 0, 90, 2.5, 4);
                TH1D *histSEPPProjPt = (TH1D*) ProjectTHnSparse(histSparseSEPP, 2, 5, PtBinsMin[iPt], PtBinsMax[iPt], 0, 90, 2.5, 4);
                TH1D *histSEMMProjPt = (TH1D*) ProjectTHnSparse(histSparseSEMM, 2, 5, PtBinsMin[iPt], PtBinsMax[iPt], 0, 90, 2.5, 4);
                TH1D *histMEPMProjPt = (TH1D*) ProjectTHnSparse(histSparseMEPM, 2, 5, PtBinsMin[iPt], PtBinsMax[iPt], 0, 90, 2.5, 4);
                TH1D *histMEPPProjPt = (TH1D*) ProjectTHnSparse(histSparseMEPP, 2, 5, PtBinsMin[iPt], PtBinsMax[iPt], 0, 90, 2.5, 4);
                TH1D *histMEMMProjPt = (TH1D*) ProjectTHnSparse(histSparseMEMM, 2, 5, PtBinsMin[iPt], PtBinsMax[iPt], 0, 90, 2.5, 4);

                histSEPMProjPt -> Write(Form("Mass_Pt_%.0f_%.0f_SEPM", PtBinsMin[iPt], PtBinsMax[iPt]));
                histSEPPProjPt -> Write(Form("Mass_Pt_%.0f_%.0f_SEPP", PtBinsMin[iPt], PtBinsMax[iPt]));
                histSEMMProjPt -> Write(Form("Mass_Pt_%.0f_%.0f_SEMM", PtBinsMin[iPt], PtBinsMax[iPt]));
                histMEPMProjPt -> Write(Form("Mass_Pt_%.0f_%.0f_MEPM", PtBinsMin[iPt], PtBinsMax[iPt]));
                histMEPPProjPt -> Write(Form("Mass_Pt_%.0f_%.0f_MEPP", PtBinsMin[iPt], PtBinsMax[iPt]));
                histMEMMProjPt -> Write(Form("Mass_Pt_%.0f_%.0f_MEMM", PtBinsMin[iPt], PtBinsMax[iPt]));
            }
            for (int iPt = 0;iPt < sizePtBins2;iPt++) {
                TH1D *histSEPMProjPt2 = (TH1D*) ProjectTHnSparse(histSparseSEPM, 2, 5, PtBinsMin2[iPt], PtBinsMax2[iPt], 0, 90, 2.5, 4);
                TH1D *histSEPPProjPt2 = (TH1D*) ProjectTHnSparse(histSparseSEPP, 2, 5, PtBinsMin2[iPt], PtBinsMax2[iPt], 0, 90, 2.5, 4);
                TH1D *histSEMMProjPt2 = (TH1D*) ProjectTHnSparse(histSparseSEMM, 2, 5, PtBinsMin2[iPt], PtBinsMax2[iPt], 0, 90, 2.5, 4);
                TH1D *histMEPMProjPt2 = (TH1D*) ProjectTHnSparse(histSparseMEPM, 2, 5, PtBinsMin2[iPt], PtBinsMax2[iPt], 0, 90, 2.5, 4);
                TH1D *histMEPPProjPt2 = (TH1D*) ProjectTHnSparse(histSparseMEPP, 2, 5, PtBinsMin2[iPt], PtBinsMax2[iPt], 0, 90, 2.5, 4);
                TH1D *histMEMMProjPt2 = (TH1D*) ProjectTHnSparse(histSparseMEMM, 2, 5, PtBinsMin2[iPt], PtBinsMax2[iPt], 0, 90, 2.5, 4);

                histSEPMProjPt2 -> Write(Form("Mass_Pt_%.0f_%.0f_SEPM", PtBinsMin2[iPt], PtBinsMax2[iPt]));
                histSEPPProjPt2 -> Write(Form("Mass_Pt_%.0f_%.0f_SEPP", PtBinsMin2[iPt], PtBinsMax2[iPt]));
                histSEMMProjPt2 -> Write(Form("Mass_Pt_%.0f_%.0f_SEMM", PtBinsMin2[iPt], PtBinsMax2[iPt]));
                histMEPMProjPt2 -> Write(Form("Mass_Pt_%.0f_%.0f_MEPM", PtBinsMin2[iPt], PtBinsMax2[iPt]));
                histMEPPProjPt2 -> Write(Form("Mass_Pt_%.0f_%.0f_MEPP", PtBinsMin2[iPt], PtBinsMax2[iPt]));
                histMEMMProjPt2 -> Write(Form("Mass_Pt_%.0f_%.0f_MEMM", PtBinsMin2[iPt], PtBinsMax2[iPt]));
            }

            histSEPMProjInt -> Write("Mass_Int_SEPM");
            histSEPPProjInt -> Write("Mass_Int_SEPP");
            histSEMMProjInt -> Write("Mass_Int_SEMM");

            histMEPMProjInt -> Write("Mass_Int_MEPM");
            histMEPPProjInt -> Write("Mass_Int_MEPP");
            histMEMMProjInt -> Write("Mass_Int_MEMM");

            //pT
            histSEPMProjPt -> Write("Mass_Pt_SEPM"); 
            histSEPPProjPt -> Write("Mass_Pt_SEPP"); 
            histSEMMProjPt -> Write("Mass_Pt_SEMM");

            histMEPMProjPt -> Write("Mass_Pt_MEPM");
            histMEPPProjPt -> Write("Mass_Pt_MEPP");
            histMEMMProjPt -> Write("Mass_Pt_MEMM");

            //rapidity 
            histSEPMProjRap -> Write("Mass_Rap_SEPM");
            histSEPPProjRap -> Write("Mass_Rap_SEPP");
            histSEMMProjRap -> Write("Mass_Rap_SEMM");

            histMEPMProjRap -> Write("Mass_Rap_MEPM");
            histMEPPProjRap -> Write("Mass_Rap_MEPP");
            histMEMMProjRap -> Write("Mass_Rap_MEMM");

            //fOut -> ls();
            fOut -> Close();
        }
        else{
            //---------------------------------------------------------------------------------------------//
                                            //Same event pairing

            TH2D* histSEPM = (TH2D*) hListSEPM -> FindObject(histName.c_str());
            TH2D* histSEMM = (TH2D*) hListSEMM -> FindObject(histName.c_str());
            TH2D* histSEPP = (TH2D*) hListSEPP -> FindObject(histName.c_str());
            //pT
            TH2D* histPtSEPM = (TH2D*) hListSEPM -> FindObject(histNamePt.c_str());
            TH2D* histPtSEMM = (TH2D*) hListSEMM -> FindObject(histNamePt.c_str());
            TH2D* histPtSEPP = (TH2D*) hListSEPP -> FindObject(histNamePt.c_str());
            //rapidity
            TH2D* histRapSEPM = (TH2D*) hListSEPM -> FindObject(histNameRap.c_str());
            TH2D* histRapSEMM = (TH2D*) hListSEMM -> FindObject(histNameRap.c_str());
            TH2D* histRapSEPP = (TH2D*) hListSEPP -> FindObject(histNameRap.c_str());

            int nBinsSE = histSEPM -> GetXaxis() -> GetNbins();
            TH1D *histSEPMProjInt = (TH1D*) ProjectTH2(histSEPM, "SEPM", 0., 90.);
            TH1D *histSEMMProjInt = (TH1D*) ProjectTH2(histSEMM, "SEMM", 0., 90.);
            TH1D *histSEPPProjInt = (TH1D*) ProjectTH2(histSEPP, "SEPP",  0., 90.);
            //pT
            TH1D *histSEPMProjPt = (TH1D*) ProjectTH2(histPtSEPM, "SEPM", 0., 20.);
            TH1D *histSEMMProjPt = (TH1D*) ProjectTH2(histPtSEMM, "SEMM", 0., 20.);
            TH1D *histSEPPProjPt = (TH1D*) ProjectTH2(histPtSEPP, "SEPP",  0., 20.);
            //rapidity
            TH1D *histSEPMProjRap = (TH1D*) ProjectTH2(histRapSEPM, "SEPM", 2.5, 4.);
            TH1D *histSEMMProjRap = (TH1D*) ProjectTH2(histRapSEMM, "SEMM", 2.5, 4.);
            TH1D *histSEPPProjRap = (TH1D*) ProjectTH2(histRapSEPP, "SEPP",  2.5, 4.);

            //---------------------------------------------------------------------------------------------//
                                                //Event mixing

            TH2D* histMEPM = (TH2D*) hListMEPM -> FindObject(histName.c_str());
            TH2D* histMEMM = (TH2D*) hListMEMM -> FindObject(histName.c_str());
            TH2D* histMEPP = (TH2D*) hListMEPP -> FindObject(histName.c_str());
            //pT
            TH2D* histPtMEPM = (TH2D*) hListMEPM -> FindObject(histNamePt.c_str());
            TH2D* histPtMEMM = (TH2D*) hListMEMM -> FindObject(histNamePt.c_str());
            TH2D* histPtMEPP = (TH2D*) hListMEPP -> FindObject(histNamePt.c_str());
            //rapidity
            TH2D* histRapMEPM = (TH2D*) hListMEPM -> FindObject(histNameRap.c_str());
            TH2D* histRapMEMM = (TH2D*) hListMEMM -> FindObject(histNameRap.c_str());
            TH2D* histRapMEPP = (TH2D*) hListMEPP -> FindObject(histNameRap.c_str());

            int nBinsME = histMEPM -> GetXaxis() -> GetNbins();
            TH1D *histMEPMProjInt = (TH1D*) ProjectTH2(histMEPM, "MEPM", 0., 90.);
            TH1D *histMEMMProjInt = (TH1D*) ProjectTH2(histMEMM, "MEMM", 0., 90.);
            TH1D *histMEPPProjInt = (TH1D*) ProjectTH2(histMEPP, "MEPP", 0., 90.);
            //pT
            TH1D *histMEPMProjPt = (TH1D*) ProjectTH2(histPtMEPM, "MEPM", 0., 20.);
            TH1D *histMEMMProjPt = (TH1D*) ProjectTH2(histPtMEMM, "MEMM", 0., 20.);
            TH1D *histMEPPProjPt = (TH1D*) ProjectTH2(histPtMEPP, "MEPP",  0., 20.);
            //rapidity
            TH1D *histMEPMProjRap = (TH1D*) ProjectTH2(histRapMEPM, "MEPM", 2.5, 4.);
            TH1D *histMEMMProjRap = (TH1D*) ProjectTH2(histRapMEMM, "MEMM", 2.5, 4.);
            TH1D *histMEPPProjRap = (TH1D*) ProjectTH2(histRapMEPP, "MEPP",  2.5, 4.);

            //---------------------------------------------------------------------------------------------//
                                //Centrality cuts on projected histograms

            for (int iCentr = 0;iCentr < sizeCentrCuts;iCentr++) {
                TH1D *histSEPMProjCentr = (TH1D*) ProjectTH2(histSEPM, "SEPM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]);
                TH1D *histSEPPProjCentr = (TH1D*) ProjectTH2(histSEPP, "SEPP", centralityCutsMin[iCentr], centralityCutsMax[iCentr]);
                TH1D *histSEMMProjCentr = (TH1D*) ProjectTH2(histSEMM, "SEMM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]);
                TH1D *histMEPMProjCentr = (TH1D*) ProjectTH2(histMEPM, "MEPM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]);
                TH1D *histMEPPProjCentr = (TH1D*) ProjectTH2(histMEPP, "MEPP", centralityCutsMin[iCentr], centralityCutsMax[iCentr]);
                TH1D *histMEMMProjCentr = (TH1D*) ProjectTH2(histMEMM, "MEMM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]);

                histSEPMProjCentr -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEPM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histSEPPProjCentr -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEPP", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histSEMMProjCentr -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEMM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histMEPMProjCentr -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEPM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histMEPPProjCentr -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEPP", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
                histMEMMProjCentr -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEMM", centralityCutsMin[iCentr], centralityCutsMax[iCentr]));
            }

            for (int iCentr = 0;iCentr < sizeCentrCuts2;iCentr++) {
                TH1D *histSEPMProjCentr2 = (TH1D*) ProjectTH2(histSEPM, "SEPM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]);
                TH1D *histSEPPProjCentr2 = (TH1D*) ProjectTH2(histSEPP, "SEPP", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]);
                TH1D *histSEMMProjCentr2 = (TH1D*) ProjectTH2(histSEMM, "SEMM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]);
                TH1D *histMEPMProjCentr2 = (TH1D*) ProjectTH2(histMEPM, "MEPM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]);
                TH1D *histMEPPProjCentr2 = (TH1D*) ProjectTH2(histMEPP, "MEPP", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]);
                TH1D *histMEMMProjCentr2 = (TH1D*) ProjectTH2(histMEMM, "MEMM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]);

                histSEPMProjCentr2 -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEPM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histSEPPProjCentr2 -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEPP", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histSEMMProjCentr2 -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_SEMM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histMEPMProjCentr2 -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEPM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histMEPPProjCentr2 -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEPP", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
                histMEMMProjCentr2 -> Write(Form("Mass_CentrFT0C_%.0f_%.0f_MEMM", centralityCutsMin2[iCentr], centralityCutsMax2[iCentr]));
            }

            histSEPMProjInt -> Write("Mass_Int_SEPM");
            histSEPPProjInt -> Write("Mass_Int_SEPP");
            histSEMMProjInt -> Write("Mass_Int_SEMM");
            //pT
            histSEPMProjPt -> Write("Mass_Pt_SEPM");
            histSEPPProjPt -> Write("Mass_Pt_SEPP");
            histSEMMProjPt -> Write("Mass_Pt_SEMM");
            //Rapidity
            histSEPMProjRap -> Write("Mass_Rap_SEPM");
            histSEPPProjRap -> Write("Mass_Rap_SEPP");
            histSEMMProjRap -> Write("Mass_Rap_SEMM");

            histMEPMProjInt -> Write("Mass_Int_MEPM");
            histMEPPProjInt -> Write("Mass_Int_MEPP");
            histMEMMProjInt -> Write("Mass_Int_MEMM");
            //pT
            histMEPMProjPt -> Write("Mass_Pt_MEPM");
            histMEPPProjPt -> Write("Mass_Pt_MEPP");
            histMEMMProjPt -> Write("Mass_Pt_MEMM");
            //Rapidity
            histMEPMProjRap -> Write("Mass_Rap_MEPM");
            histMEPPProjRap -> Write("Mass_Rap_MEPP");
            histMEMMProjRap -> Write("Mass_Rap_MEMM");

            //fOut -> ls();
            fOut -> Close();

        }

        fIn -> Close();
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TCanvas *produceCanvas(bool logX, bool logY,TString name) {
	TCanvas *c = new TCanvas(name,name);
	c->SetFrameFillColor(0); //Transparent background
    c->SetFrameFillStyle(0);
    c->SetFrameBorderMode(0);
	c->SetGridx(); //Grid on x
	c->SetGridy(); //Grid on y 
	if (logX) c->SetLogx(); //Log on x, if needed
	if (logY) c->SetLogy(); //Log on y, if needed
	return c;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TLegend *produceLegend(){
    TLegend *legend = new TLegend(0.22,0.81,0.40,0.91,"","rNDC");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(62);
    legend->SetTextSize(0.04);
    return legend;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D *produceTH1D(string name, string title, string xAxisTitle, string yAxisTitle, int nBins, double xMin, double xMax) {
	TH1D *h = new TH1D(name.c_str(), title.c_str(), nBins, xMin, xMax);
	h->GetXaxis()->SetTitle(xAxisTitle.c_str());
	h->GetYaxis()->SetTitle(yAxisTitle.c_str());
	h->GetXaxis()->SetTitleFont(62);
	h->GetXaxis()->SetLabelFont(62);
	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->SetLabelFont(62);
	h->GetXaxis()->CenterTitle(true);
	h->GetYaxis()->CenterTitle(true);
	return h;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TH2D *produceTH2D(string name, string title, string xAxisTitle, string yAxisTitle, int nBinsX, double xMin, double xMax, int nBinsY, double yMin, double yMax) {
	TH2D *h = new TH2D(name.c_str(), title.c_str(), nBinsX, xMin, xMax, nBinsY, yMin, yMax);
	h->GetXaxis()->SetTitle(xAxisTitle.c_str());
	h->GetYaxis()->SetTitle(yAxisTitle.c_str());
	h->GetXaxis()->SetTitleFont(62);
	h->GetXaxis()->SetLabelFont(62);
	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->SetLabelFont(62);
	h->GetXaxis()->CenterTitle(true);
	h->GetYaxis()->CenterTitle(true);
	return h;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D* ProjectTH2(TH2D *hist2D, string title, double minCentr, double maxCentr) {
    double minCentrBin = hist2D -> GetYaxis() -> FindBin(minCentr);
    double maxCentrBin = hist2D -> GetYaxis() -> FindBin(maxCentr - 0.01);

    //Printf("minPtBin = %.0f, maxPtBin = %.0f", minCentrBin, maxCentrBin);
    hist2D -> GetYaxis() -> SetRange(minCentrBin, maxCentrBin);

    TH1D *histProj = (TH1D*) hist2D -> ProjectionX(Form("histProj_%s_%.0f_%.0f", title.c_str(), minCentrBin, maxCentrBin));
    //histProj -> SetTitle(title.c_str());
    return histProj;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D* ProjectTHnSparse(THnSparseD *histSparse, double minVar1, double maxVar1, double minVar2, double maxVar2, double minVar3, double maxVar3, double minVar4, double maxVar4) {
    // Mass bins
    double minVar1Bin = histSparse -> GetAxis(0) -> FindBin(minVar1);
    double maxVar1Bin = histSparse -> GetAxis(0) -> FindBin(maxVar1 - 0.01);
    // Pt bins
    double minVar2Bin = histSparse -> GetAxis(1) -> FindBin(minVar2);
    double maxVar2Bin = histSparse -> GetAxis(1) -> FindBin(maxVar2 - 0.01);
    // Centrality bins
    double minVar3Bin = histSparse -> GetAxis(2) -> FindBin(minVar3);
    double maxVar3Bin = histSparse -> GetAxis(2) -> FindBin(maxVar3 - 0.01);
    // Rapidity bins
    double minVar4Bin = histSparse -> GetAxis(3) -> FindBin(minVar4);
    double maxVar4Bin = histSparse -> GetAxis(3) -> FindBin(maxVar4 - 0.01);

    histSparse -> GetAxis(0) -> SetRange(minVar1Bin, maxVar1Bin);
    histSparse -> GetAxis(1) -> SetRange(minVar2Bin, maxVar2Bin);
    histSparse -> GetAxis(2) -> SetRange(minVar3Bin, maxVar3Bin);
    histSparse -> GetAxis(3) -> SetRange(minVar4Bin, maxVar4Bin);

    //Printf("%3.2f - %3.2f ## %3.2f - %3.2f ## %3.2f - %3.2f \n", minVar1Bin, maxVar1Bin, minVar2Bin, maxVar2Bin, minVar3Bin, maxVar3Bin);

    TH1D *histProj = (TH1D*) histSparse -> Projection(0, Form("histProj_%.0f_%.0f__%.0f_%.0f__%.0f_%.0f", minVar1, maxVar1, minVar2, maxVar2, minVar3, maxVar3));
    return histProj;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D* ProjectTHnSparsePt(THnSparseD *histSparse, double minVar1, double maxVar1, double minVar2, double maxVar2, double minVar3, double maxVar3, double minVar4, double maxVar4) {
    // Mass bins
    double minVar1Bin = histSparse -> GetAxis(0) -> FindBin(minVar1);
    double maxVar1Bin = histSparse -> GetAxis(0) -> FindBin(maxVar1 - 0.01);
    // Pt bins
    double minVar2Bin = histSparse -> GetAxis(1) -> FindBin(minVar2);
    double maxVar2Bin = histSparse -> GetAxis(1) -> FindBin(maxVar2 - 0.01);
    // Centrality bins
    double minVar3Bin = histSparse -> GetAxis(2) -> FindBin(minVar3);
    double maxVar3Bin = histSparse -> GetAxis(2) -> FindBin(maxVar3 - 0.01);
    // Rapidity bins
    double minVar4Bin = histSparse -> GetAxis(3) -> FindBin(minVar4);
    double maxVar4Bin = histSparse -> GetAxis(3) -> FindBin(maxVar4 - 0.01);

    histSparse -> GetAxis(0) -> SetRange(minVar1Bin, maxVar1Bin);
    histSparse -> GetAxis(1) -> SetRange(minVar2Bin, maxVar2Bin);
    histSparse -> GetAxis(2) -> SetRange(minVar3Bin, maxVar3Bin);
    histSparse -> GetAxis(3) -> SetRange(minVar4Bin, maxVar4Bin);

    //Printf("%3.2f - %3.2f ## %3.2f - %3.2f ## %3.2f - %3.2f \n", minVar1Bin, maxVar1Bin, minVar2Bin, maxVar2Bin, minVar3Bin, maxVar3Bin);

    TH1D *histProj = (TH1D*) histSparse -> Projection(1, Form("histProj_%.0f_%.0f__%.0f_%.0f__%.0f_%.0f", minVar1, maxVar1, minVar2, maxVar2, minVar3, maxVar3));
    return histProj;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D* ProjectTHnSparseRap(THnSparseD *histSparse, double minVar1, double maxVar1, double minVar2, double maxVar2, double minVar3, double maxVar3, double minVar4, double maxVar4) {
    // Mass bins
    double minVar1Bin = histSparse -> GetAxis(0) -> FindBin(minVar1);
    double maxVar1Bin = histSparse -> GetAxis(0) -> FindBin(maxVar1 - 0.01);
    // Pt bins
    double minVar2Bin = histSparse -> GetAxis(1) -> FindBin(minVar2);
    double maxVar2Bin = histSparse -> GetAxis(1) -> FindBin(maxVar2 - 0.01);
    // Centrality bins
    double minVar3Bin = histSparse -> GetAxis(2) -> FindBin(minVar3);
    double maxVar3Bin = histSparse -> GetAxis(2) -> FindBin(maxVar3 - 0.01);
    // Rapidity bins
    double minVar4Bin = histSparse -> GetAxis(3) -> FindBin(minVar4);
    double maxVar4Bin = histSparse -> GetAxis(3) -> FindBin(maxVar4 - 0.01);

    histSparse -> GetAxis(0) -> SetRange(minVar1Bin, maxVar1Bin);
    histSparse -> GetAxis(1) -> SetRange(minVar2Bin, maxVar2Bin);
    histSparse -> GetAxis(2) -> SetRange(minVar3Bin, maxVar3Bin);
    histSparse -> GetAxis(3) -> SetRange(minVar4Bin, maxVar4Bin);

    //Printf("%3.2f - %3.2f ## %3.2f - %3.2f ## %3.2f - %3.2f \n", minVar1Bin, maxVar1Bin, minVar2Bin, maxVar2Bin, minVar3Bin, maxVar3Bin);

    TH1D *histProj = (TH1D*) histSparse -> Projection(3, Form("histProj_%.0f_%.0f__%.0f_%.0f__%.0f_%.0f", minVar1, maxVar1, minVar2, maxVar2, minVar3, maxVar3));
    return histProj;
}