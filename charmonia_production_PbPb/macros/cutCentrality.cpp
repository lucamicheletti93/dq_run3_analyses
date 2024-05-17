#include <stdio.h>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "THStack.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include "TNtuple.h"
#include "TH1F.h"
#include "TH3.h"
#include "TH2.h"
#include <TF1.h>
#include "TGraphErrors.h"
#include <TMultiGraph.h>
#include <TLegend.h>
#include "TStyle.h"
#include <vector>
#include <TGraph.h>
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TBrowser.h"
#include "TOrdCollection.h"
#include "TList.h"
#include <TString.h>
#include <algorithm>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include "TPaveText.h"
#include <numeric>

using namespace std;
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
TLegend *produceLegend(){
    TLegend *legend = new TLegend(0.22,0.81,0.40,0.91,"","rNDC");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(62);
    legend->SetTextSize(0.04);
    return legend;
}
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

struct cuts
{
    //vector <Int_t> cuts_10_20_30 = {10,20,30};
    //vector <Int_t> cutsAll = {0,20,40,60,90};
    //vector <Int_t> cutsBins = {0,10,20,30,40,50,60,70,80,90};
};

void cutCentrality(int numberOfCuts){
    vector <Int_t> cuts_10_20_30 = {10,20,30};
    vector <Int_t> cutsAll = {0,20,40,60,90};
    vector <Int_t> cutsBins = {0,10,20,30,40,50,60,70,80,90};
    int size = cutsAll.size();
    int numHistos = size - 1;
    int sizeBins = cutsBins.size();
    int numHistosBins = sizeBins -1;
    TFile* analysisResults = TFile::Open("AnalysisResults.root");
    //Same event pairing
    TDirectory* dir = (TDirectory*)analysisResults->Get("analysis-same-event-pairing");
    //cout << "TDirectory: " << dir << endl;
    TList *hList = (TList*)dir->Get(Form("output"));
    TList *hList2 = (TList*)hList->At(18);
    //cout << hList << endl;
    //cout << "Name: " << hList->At(18)->GetName() << endl;

    TList *hListPP = (TList*)hList->FindObject("PairsMuonSEPP_muonLowPt510SigmaPDCA");
    TList *hListMM = (TList*)hList->FindObject("PairsMuonSEMM_muonLowPt510SigmaPDCA");
    //cout << hList << endl;
    //cout << "Name PP: " << hList->At(19)->GetName() << endl;
    //cout << "Name MM: " << hList->At(20)->GetName() << endl;
    TH1D *histosBins[numHistosBins];
    TCanvas *canvasBins[numHistosBins];
    TCanvas *cTot = produceCanvas(false,false,"All cuts");
    cTot->Divide(3,3);
    TH2D *massVsCentrality = produceTH2D("massVsCentrality","Mass vs Centrality","mass [GeV/cm^{2}]","centrality FT0C (%)",15,0.,15.,100,0.,100.);
    TH2D *massVsCentralityAll = produceTH2D("massVsCentralityAll","Mass vs Centrality All","mass [GeV/cm^{2}]","centrality FT0C (%)",15,0.,15.,100,0.,100.);
    massVsCentrality= (TH2D*)hList2->At(2);
    massVsCentralityAll = (TH2D*)massVsCentrality->Clone();
    massVsCentralityAll->SetDirectory(0);
    TCanvas *canvasNoCuts = produceCanvas(false,false,"massVsCentralityAll");
    massVsCentralityAll->Draw();
    TH2D *massVsCentrality3Cuts[3];
    TH2D *massVsCentralityBinsCuts[numHistosBins];
    TH2D *massVsCentralityAllCuts[numHistos];
    if(numberOfCuts == 3){
        TCanvas *canvas3Cuts[3];
        for (int i_cuts = 0; i_cuts < 3; i_cuts++){
            //massVsCentrality3Cuts[i_cuts] = produceTH2D(Form("massVsCentrality%i",cuts_10_20_30.at(i_cuts)),Form("Mass vs Centrality %i",cuts_10_20_30.at(i_cuts)),"mass [GeV/cm^{2}]","centrality FT0C (%)",15,0.,15.,100,0.,100.);
            massVsCentrality3Cuts[i_cuts] = new TH2D(Form("massVsCentrality%i",cuts_10_20_30.at(i_cuts)),Form("Mass vs Centrality %i",cuts_10_20_30.at(i_cuts)),15,0.,15.,100,0.,100.);
            massVsCentrality3Cuts[i_cuts] = (TH2D*)massVsCentrality->Clone();
            massVsCentrality3Cuts[i_cuts]->SetDirectory(0);
            canvas3Cuts[i_cuts] = produceCanvas(false,false,Form("massVsCentrality%i",cuts_10_20_30.at(i_cuts)));
            massVsCentrality3Cuts[i_cuts]->SetName(Form("mass_%i",cuts_10_20_30.at(i_cuts)));
            //massVsCentrality3Cuts[i_cuts]->SetAxisRange(cuts_10_20_30.at(i_cuts), 70.,"Y");
            massVsCentrality3Cuts[i_cuts]->GetYaxis()->SetRangeUser(cuts_10_20_30.at(i_cuts),70);
            massVsCentrality3Cuts[i_cuts]->Draw();
        }
    }
    else if(numberOfCuts == 10){
        TCanvas *canvasBinsCuts[numHistosBins];
        for (int i_cuts = 0; i_cuts < numHistosBins; i_cuts++){
            //massVsCentralityBinsCuts[i_cuts] = produceTH2D(Form("massVsCentrality%i",cuts_10_20_30.at(i_cuts)),Form("Mass vs Centrality %i",cuts_10_20_30.at(i_cuts)),"mass [GeV/cm^{2}]","centrality FT0C (%)",15,0.,15.,100,0.,100.);
            massVsCentralityBinsCuts[i_cuts] = new TH2D(Form("massVsCentrality%i",cutsBins.at(i_cuts)),Form("Mass vs Centrality %i_%i",cutsBins.at(i_cuts),cutsBins.at(i_cuts+1)),15,0.,15.,100,0.,100.);
            massVsCentralityBinsCuts[i_cuts] = (TH2D*)massVsCentrality->Clone();
            massVsCentralityBinsCuts[i_cuts]->SetDirectory(0);
            canvasBinsCuts[i_cuts] = produceCanvas(false,false,Form("massVsCentrality%i_%i",cutsBins.at(i_cuts),cutsBins.at(i_cuts+1)));
            massVsCentralityBinsCuts[i_cuts]->SetName(Form("mass_%i_%i",cutsBins.at(i_cuts),cutsBins.at(i_cuts+1)));
            //massVsCentralityBinsCuts[i_cuts]->SetAxisRange(cutsBins.at(i_cuts), 70.,"Y");
            massVsCentralityBinsCuts[i_cuts]->GetYaxis()->SetRangeUser(cutsBins.at(i_cuts),cutsBins.at(i_cuts+1));
            massVsCentralityBinsCuts[i_cuts]->Draw();
        }
    }
    else{
        TCanvas *canvasAllCuts[numHistos];
        for (int i_cuts = 0; i_cuts < numHistos; i_cuts++){
            massVsCentralityAllCuts[i_cuts] = produceTH2D(Form("massVsCentrality%i",cutsAll.at(i_cuts)),Form("Mass vs Centrality %i",cutsAll.at(i_cuts)),"mass [GeV/cm^{2}]","centrality FT0C (%)",15,0.,15.,100,0.,100.);
            massVsCentralityAllCuts[i_cuts] = (TH2D*)massVsCentrality->Clone();
            massVsCentralityAllCuts[i_cuts]->SetDirectory(0);
            canvasAllCuts[i_cuts] = produceCanvas(false,false,Form("massVsCentrality_%i_%i",cutsAll.at(i_cuts),cutsAll.at(i_cuts+1)));
            cout << "i_cuts: " << i_cuts << "   cutsAll.at(i_cuts): " << cutsAll.at(i_cuts) << "   cutsAll.at(i_cuts+1): " << cutsAll.at(i_cuts+1) << endl;
            //massVsCentralityAllCuts[i_cuts]->Draw();
            massVsCentralityAllCuts[i_cuts]->GetYaxis()->SetRangeUser(cutsAll.at(i_cuts),cutsAll.at(i_cuts+1));
        }
    }
    //Event mixing
    TDirectory* dirEvMix = (TDirectory*)analysisResults->Get("analysis-event-mixing");
    TList *hListEvMix = (TList*)dirEvMix->Get(Form("output"));
    TList *hListEvMix2 = (TList*)hListEvMix->At(9);
    TH2D *massVsCentralityEvMix = produceTH2D("massVsCentralityEvMix","Mass vs Centrality Event Mixing","mass [GeV/cm^{2}]","centrality FT0C (%)",15,0.,15.,100,0.,100.);
    TH2D *massVsCentralityEvMixClone = produceTH2D("massVsCentralityEvMixClone","Mass vs Centrality Event Mixing Clone","mass [GeV/cm^{2}]","centrality FT0C (%)",15,0.,15.,100,0.,100.);
    massVsCentralityEvMix= (TH2D*)hListEvMix2->At(2);
    massVsCentralityEvMixClone = (TH2D*)massVsCentralityEvMix->Clone();
    massVsCentralityEvMixClone->SetDirectory(0);

    //PP
    TH2D *massVsCentralityPP = produceTH2D("massVsCentralityPP","Mass vs Centrality PP","mass [GeV/cm2]","centrality FT0C (%)",15,0.,15.,100,0.,100.);
    TH2D *massVsCentralityPPClone = produceTH2D("massVsCentralityPPClone","Mass vs Centrality PP Clone","mass [GeV/cm2]","centrality FT0C (%)",15,0.,15.,100,0.,100.);
    massVsCentralityPP= (TH2D*)hListPP->FindObject("Mass_CentFT0C");
    massVsCentralityPPClone = (TH2D*)massVsCentralityPP->Clone();
    massVsCentralityPPClone->SetDirectory(0);
    TCanvas *cProvaPP = produceCanvas(false,false,"cProvaPP_pro");
    TH1D *provaPPMass = new TH1D("provaPPMass","PP projection mass",15,0.,15.);
    provaPPMass = massVsCentralityPPClone->ProjectionX();
    provaPPMass->SetDirectory(0);
    cProvaPP->cd();
    cProvaPP->Draw();
    provaPPMass->Draw();
    //MM
    TH2D *massVsCentralityMM = produceTH2D("massVsCentralityMM","Mass vs Centrality MM","mass [GeV/cm2]","centrality FT0C (%)",15,0.,15.,100,0.,100.);
    TH2D *massVsCentralityMMClone = produceTH2D("massVsCentralityMMClone","Mass vs Centrality MM Clone","mass [GeV/cm2]","centrality FT0C (%)",15,0.,15.,100,0.,100.);
    massVsCentralityMM= (TH2D*)hListMM->FindObject("Mass_CentFT0C");
    massVsCentralityMMClone = (TH2D*)massVsCentralityMM->Clone();
    massVsCentralityMMClone->SetDirectory(0);
    TCanvas *cProvaMM = produceCanvas(false,false,"cProvaMM_pro");
    TH1D *provaMMMass = new TH1D("provaMMMass","MM projection mass",15,0.,15.);
    provaMMMass = massVsCentralityMMClone->ProjectionX();
    provaMMMass->SetDirectory(0);
    cProvaMM->cd();
    cProvaMM->Draw();
    provaMMMass->Draw();

    analysisResults->Close();
    TFile *centralityCuts = TFile::Open("centralityCuts.root", "UPDATE");
    cProvaPP->cd();
    cProvaPP->Write(0,2,0);
    cProvaMM->cd();
    cProvaMM->Write(0,2,0);
    if(numberOfCuts == 3){
        TH1D *mass[3];
        TH1D *centrality[3];
        TCanvas *cMass[3];
        TCanvas *cCentrality[3];
        TLegend *legend3Cuts = produceLegend();
        for (int i_cuts = 0; i_cuts < 3; i_cuts++){
            massVsCentrality3Cuts[i_cuts]->Write(0,2,0);
            cMass[i_cuts] = produceCanvas(false,false,Form("mass_%i",cuts_10_20_30.at(i_cuts)));
            //mass[i_cuts] = massVsCentrality3Cuts[i_cuts]->ProjectionX();
            mass[i_cuts] = massVsCentrality3Cuts[i_cuts]->ProjectionX(Form("mass_%d",i_cuts), cuts_10_20_30.at(i_cuts),70);
            mass[i_cuts]->SetDirectory(0);
            mass[i_cuts]->SetXTitle("mass [GeV/c^{2}]");
            mass[i_cuts]->SetLineColor(kBlue);
            legend3Cuts->AddEntry(mass[i_cuts],Form("(%i-70)%%",cuts_10_20_30.at(i_cuts)),"lep");
            mass[i_cuts]->Draw();
            mass[i_cuts]->SetName(Form("mass_%i",cuts_10_20_30.at(i_cuts)));
            mass[i_cuts]->Write(0,2,0);
            cCentrality[i_cuts] = produceCanvas(false,false,Form("centrality_%i",cuts_10_20_30.at(i_cuts)));
            centrality[i_cuts] = massVsCentrality3Cuts[i_cuts]->ProjectionY();
            centrality[i_cuts]->SetDirectory(0);
            centrality[i_cuts]->SetXTitle("centrality");
            centrality[i_cuts]->Draw();
            centrality[i_cuts]->SetName(Form("centrality_%i",cuts_10_20_30.at(i_cuts)));
            centrality[i_cuts]->Write(0,2,0);
        }
    }
    else if(numberOfCuts == 10){
        TH1D *mass[sizeBins];
        TH1D *centrality[sizeBins];
        TCanvas *cMass[sizeBins];
        TCanvas *cCentrality[sizeBins];
        TLegend *legend3Cuts = produceLegend();
        for (int i_cuts = 0; i_cuts < numHistosBins; i_cuts++){
            massVsCentralityBinsCuts[i_cuts]->Write(0,2,0);
            cMass[i_cuts] = produceCanvas(false,false,Form("mass_%i_%i",cutsBins.at(i_cuts),cutsBins.at(i_cuts+1)));
            //mass[i_cuts] = massVsCentrality3Cuts[i_cuts]->ProjectionX();
            mass[i_cuts] = massVsCentralityBinsCuts[i_cuts]->ProjectionX(Form("mass_%d",i_cuts), cutsBins.at(i_cuts),cutsBins.at(i_cuts+1));
            mass[i_cuts]->SetDirectory(0);
            mass[i_cuts]->SetXTitle("mass [GeV/c^{2}]");
            mass[i_cuts]->SetLineColor(kBlue);
            legend3Cuts->AddEntry(mass[i_cuts],Form("(%i-%i)%%",cutsBins.at(i_cuts),cutsBins.at(i_cuts+1)),"lep");
            mass[i_cuts]->Draw();
            mass[i_cuts]->SetName(Form("mass_%i_%i",cutsBins.at(i_cuts),cutsBins.at(i_cuts+1)));
            mass[i_cuts]->Write(0,2,0);
            cTot->Update();
        }
        for (int i_cuts = 0; i_cuts < numHistosBins; i_cuts++){
            cTot->cd(i_cuts+1);
            mass[i_cuts]->Draw();
            mass[i_cuts]->SetTitle(Form("Mass (%i-%i)%%",cutsBins.at(i_cuts),cutsBins.at(i_cuts+1)));
            mass[i_cuts]->Write(0,2,0);
            cTot->Update();
        }
        cTot->Draw();
        cTot->SetName("All cuts");
        cTot->Write(0,2,0);
    }
    else{
        TH1D *mass[numHistos];
        TH1D *centrality[numHistos];
        TCanvas *cMass[numHistos];
        TCanvas *cCentrality[numHistos];
        TLegend *legendAllCuts = produceLegend();
        for (int i_cuts = 0; i_cuts < numHistos; i_cuts++){
            massVsCentralityAllCuts[i_cuts]->Write(0,2,0);
            cMass[i_cuts] = produceCanvas(false,false,Form("mass_%i_%i",cutsAll.at(i_cuts),cutsAll.at(i_cuts+1)));
            mass[i_cuts] = massVsCentralityAllCuts[i_cuts]->ProjectionX();
            mass[i_cuts]->SetDirectory(0);
            mass[i_cuts]->SetXTitle("mass [GeV/c^{2}]");
            mass[i_cuts]->SetLineColor(kBlue);
            legendAllCuts->AddEntry(mass[i_cuts],Form("(%i-%i)%%",cutsAll.at(i_cuts),cutsAll.at(i_cuts+1)),"lep");
            //mass[i_cuts]->Draw();
            mass[i_cuts]->SetName(Form("mass_%i_%i",cutsAll.at(i_cuts),cutsAll.at(i_cuts+1)));
            mass[i_cuts]->Write(0,2,0);
            cCentrality[i_cuts] = produceCanvas(false,false,Form("centrality_%i_%i",cutsAll.at(i_cuts),cutsAll.at(i_cuts+1)));
            centrality[i_cuts] = massVsCentralityAllCuts[i_cuts]->ProjectionY();
            centrality[i_cuts]->SetDirectory(0);
            centrality[i_cuts]->SetXTitle("centrality");
            //centrality[i_cuts]->Draw();
            centrality[i_cuts]->SetName(Form("centrality_%i_%i",cutsAll.at(i_cuts),cutsAll.at(i_cuts+1)));
            centrality[i_cuts]->Write(0,2,0);
        }
        int index_cut = 10;
        for (Int_t i_cuts = 0; i_cuts < numHistosBins; i_cuts++){
            //cout << "i_cuts: " << i_cuts << endl;
            int index_cut_up = index_cut;
            //cout << "index_cut_up: " << index_cut_up << endl;
            int index_cut_low = index_cut - 10;
            //cout << "index_cut_low: " << index_cut_low << endl;
            if (index_cut_up < 11){
                massVsCentralityAll->GetYaxis()->SetRangeUser(index_cut_low, index_cut_up);
            }
            else {
                double index_cut_low_corrected = index_cut_low + 0.01;
                massVsCentralityAll->GetYaxis()->SetRangeUser(index_cut_low_corrected, index_cut_up);
            }
            histosBins[i_cuts] = massVsCentralityAll->ProjectionX();
            histosBins[i_cuts]->SetDirectory(0);
            canvasBins[i_cuts] = produceCanvas(false,false,Form("c_%i_%i",index_cut_low,index_cut_up));
            canvasBins[i_cuts]->cd();
            //canvasBins[i_cuts]->Draw();
            histosBins[i_cuts]->SetTitle(Form("Mass_(%i-%i)%%",index_cut_low,index_cut_up));
            histosBins[i_cuts]->SetXTitle("mass [GeV/c^{2}]");
            histosBins[i_cuts]->SetLineColor(kBlue+i_cuts);
            //histosBins[i_cuts]->Draw();
            histosBins[i_cuts]->SetName(("Mass_"+to_string(index_cut_low)+"_"+to_string(index_cut_up)).c_str()); 
            histosBins[i_cuts]->Write(0,2,0);
            index_cut = index_cut + 10;
            //cout << "index_cut: " << index_cut << endl;
        }

        //create PP and MM
        TH1D *histoPPMM[4];
        TH1D *histoPPMMFin[4];
        TH2D *histoPPCloneArray[4];
        TH2D *histoMMCloneArray[4];
        TCanvas *cPPMM[4];
        TCanvas *cPPMM2[4];
        //TCanvas *cPPMM2;
        TCanvas *cPPMMProva[4];
        TH1D *histoPPcuts[4];
        TH1D *histoMMcuts[4];
        TCanvas *canvasPPcuts[4];
        TCanvas *canvasPPcutsProva[4];
        TCanvas *canvasMMcuts[4];
        TCanvas *canvasPPMMcuts[4];


        //Event mixing
        TH1D *histoEvMix[4];
        TH1D *histoEvSame[4];
        TH1D *histoEvPPMM[4];
        TH1D *histoSEME[4];
        TH1D *histoEMPMCheck[4];
        TH1D *histoSEPPMM[4];
        TCanvas *cEvPPMM[4];
        TCanvas *cEvMixPPMM[4];
        TCanvas *cEvNorm[4];
        TCanvas *cNorm[4];
        TCanvas *cFinSEME[4];
        TCanvas *cFinSEPPMM[4];
        TCanvas *cPMEMCheck[4];
        //TPad *TPad_Norm[4];
        vector <double> binContentPP;
        vector <double> binContentMM;
        vector <double> binContentPPMM;
        vector <double> binContentPPMMFin;
        vector <int> binSBLowLowSE;
        vector <int> binSBLowUpSE;
        vector <double> integralSBLowSE;
        vector <int> binSBUpLowSE;
        vector <int> binSBUpUpSE;
        vector <double> integralSBUpSE;
        vector <int> binSBLowLowME;
        vector <int> binSBLowUpME;
        vector <double> integralSBLowME;
        vector <int> binSBUpLowME;
        vector <int> binSBUpUpME;
        vector <double> integralSBUpME;
        vector <double> integralSumME;
        vector <double> integralSumSE;
        vector <double> integralSumPPMM;
        vector <double> ratio_integral_SE_ME;
        vector <int> binSBLowLowPPMM;
        vector <int> binSBLowUpPPMM;
        vector <double> integralSBLowPPMM;
        vector <int> binSBUpLowPPMM;
        vector <int> binSBUpUpPPMM;
        vector <double> integralSBUpPPMM;
        vector <double> ratio_integral_SE_PPMM;

        for (int iEvMix = 0; iEvMix < 4; iEvMix++){
            massVsCentralityEvMixClone->GetYaxis()->SetRangeUser(cutsAll.at(iEvMix),cutsAll.at(iEvMix+1));
            histoEvMix[iEvMix] = produceTH1D(Form("histoEvMix_%i_%i",cutsAll.at(iEvMix),cutsAll.at(iEvMix+1)),"Mass event mixing","mass [GeV/c^{2}]","Counts",750,0.,15.);
            histoEvMix[iEvMix] = massVsCentralityEvMixClone->ProjectionX();
            histoEvMix[iEvMix]->SetDirectory(0);
        }
        int cutLowPPMM;
        int cutHighPPMM;
        for (int j = 0; j < 4; j++){
            histoPPCloneArray[j] = produceTH2D(Form("histoPPCloneArray_%i", j),Form("histoPPCloneArray_%i",j),"mass [GeV/cm^{2}]","centrality FT0C (%)",15,0.,15.,100,0.,100.);
            histoPPCloneArray[j] = (TH2D*)massVsCentralityPP->Clone(Form("histoPPCloneArray_%i", j));
            histoPPCloneArray[j]->SetDirectory(0);
            histoMMCloneArray[j] = produceTH2D(Form("histoMMCloneArray_%i", j),Form("histoMMCloneArray_%i",j),"mass [GeV/cm^{2}]","centrality FT0C (%)",15,0.,15.,100,0.,100.);
            histoMMCloneArray[j] = (TH2D*)massVsCentralityMM->Clone(Form("histoMMCloneArray_%i", j));
            histoMMCloneArray[j]->SetDirectory(0);
        }
        int nBinPPMM = 0;
        int index = 0;
        for (int i_PP_MM = 0; i_PP_MM < 4; i_PP_MM++){
            //produce PP
            cutLowPPMM=0;
            cutHighPPMM=0;
            cutLowPPMM = cutsAll.at(i_PP_MM);
            cutHighPPMM = cutsAll.at(i_PP_MM+1);
            //cout << "i_PP_MM: " << i_PP_MM << endl;
            //cout << "cuts at i: " << cutsAll.at(i_PP_MM) << "  cuts at i+1: " << cutsAll.at(i_PP_MM+1) << endl;
            //cout << "cutLowPPMM: " << cutLowPPMM << "  cutHighPPMM: " << cutHighPPMM << endl;
            canvasPPcuts[i_PP_MM] = produceCanvas(false,false,Form("histo_PP_%i_%i",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            canvasPPcuts[i_PP_MM]->cd();
            canvasPPcuts[i_PP_MM]->Draw();
            //histoPPCloneArray[i_PP_MM]->Draw();
            histoPPCloneArray[i_PP_MM]->GetYaxis()->SetRangeUser(cutLowPPMM,cutHighPPMM);
            //histoPPcuts[i_PP_MM] = histoPPCloneArray[i_PP_MM]->ProjectionX(Form("histo_PP_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            //histoPPcuts[i_PP_MM] = histoPPCloneArray[i_PP_MM]->ProjectionX("histo_PP_all");
            //histoPPcuts[i_PP_MM] = histoPPCloneArray[i_PP_MM]->ProjectionX(Form("histo_PP_%i_%i",cutLowPPMM,cutHighPPMM));
            histoPPcuts[i_PP_MM] = histoPPCloneArray[i_PP_MM]->ProjectionX();
            histoPPcuts[i_PP_MM]->SetDirectory(0);
            histoPPcuts[i_PP_MM]->SetName(("output_PP_"+to_string(cutsAll.at(i_PP_MM))+"_"+to_string(cutsAll.at(i_PP_MM+1))).c_str());
            histoPPcuts[i_PP_MM]->Draw();
            histoPPcuts[i_PP_MM]->Write(0,2,0);
            canvasPPcuts[i_PP_MM]->Write(0,2,0);
            //produce MM
            canvasMMcuts[i_PP_MM] = produceCanvas(false,false,Form("histo_MM_%i_%i",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            canvasMMcuts[i_PP_MM]->cd();
            canvasMMcuts[i_PP_MM]->Draw();
            //histoMMCloneArray[i_PP_MM]->Draw();
            histoMMCloneArray[i_PP_MM]->GetYaxis()->SetRangeUser(cutLowPPMM,cutHighPPMM);
            //histoMMcuts[i_PP_MM] = histoMMCloneArray[i_PP_MM]->ProjectionX("histo_MM_all");
            //histoMMcuts[i_PP_MM] = histoMMCloneArray[i_PP_MM]->ProjectionX(Form("histo_MM_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            histoMMcuts[i_PP_MM] = histoMMCloneArray[i_PP_MM]->ProjectionX();
            histoMMcuts[i_PP_MM]->SetDirectory(0);
            histoMMcuts[i_PP_MM]->SetName(("output_MM_"+to_string(cutsAll.at(i_PP_MM))+"_"+to_string(cutsAll.at(i_PP_MM+1))).c_str());
            histoMMcuts[i_PP_MM]->Draw();
            histoMMcuts[i_PP_MM]->Write(0,2,0);
            canvasMMcuts[i_PP_MM]->Write(0,2,0);
            cPPMM[i_PP_MM] = produceCanvas(false,false,Form("Mass_PP_MM_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            cPPMM[i_PP_MM]->cd();
            cPPMM[i_PP_MM]->Draw();
            //produce PPMM
            //histoPPMM[i_PP_MM] = produceTH1D2("histoPPMM",Form("histoPPMM_%i_%i",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)),"mass [GeV/c^{2}]","Counts",750,0.,15.);
            //histoPPMM[i_PP_MM] = (*histoPPcuts[i_PP_MM])*(*histoMMcuts[i_PP_MM]);
            //histoPPMMFin[i_PP_MM] = produceTH1D("histoPPMMFin",Form("Mass_PP_MM_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)),"mass [GeV/cm^{2}]","centrality FT0C (%)",15,0.,15.);
            /*histoPPMM[i_PP_MM] = (TH1D*)histoPPcuts[i_PP_MM]->Clone();
            histoPPMM[i_PP_MM]->Multiply(histoMMcuts[i_PP_MM]);
            histoPPMM[i_PP_MM]->Draw();
            histoPPMM[i_PP_MM]->Write(0,2,0);*/
        }
        /*nBinPPMM = histoPPMM[0]->GetNbinsX();
        cout << "bins histo ppmm: " << histoPPMM[0]->GetNbinsX() << endl;
        histoPPMMFin[0] = produceTH1D("histoPPMMFin","Mass_PP_MM_Prova","mass [GeV/cm^{2}]","Counts",750,0.,770.);
        for (int nBins = 0; nBins < nBinPPMM; nBins ++){
            binContentPPMM.push_back(histoPPMM[0]->GetBinContent(nBins));
            cout << "binContentPPMM: " << binContentPPMM.at(nBins) << endl;
            binContentPPMMFin.push_back(2*sqrt(binContentPPMM.at(nBins)));
            cout << "binContentPPMMFin: " << binContentPPMMFin.at(nBins) << endl;
            histoPPMMFin[0]->SetBinContent(index, binContentPPMMFin.at(nBins));
            index = index +1;
            cout << "index: " << index << endl;
            cout << "nBinPPMM: " << nBinPPMM << endl;
        }*/
        //cout << "num_Histos: " << endl;
        int nBinPMEM;
        for (int i_PP_MM = 0; i_PP_MM < 4; i_PP_MM++){
            index = 1;
            nBinPPMM = 0;
            nBinPMEM = 0;
            //cout << "nBinPPMM: " << nBinPPMM << endl;
            nBinPPMM = histoMMcuts[i_PP_MM]->GetNbinsX();
            //cout << "nBinPPMM:2 " << nBinPPMM << endl;
            histoSEME[i_PP_MM] = produceTH1D(Form("Mass_EM_Sub_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)),Form("Mass_EM_Sub_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)),"mass [GeV/cm^{2}]","Counts",750,0.,15.);
            histoSEPPMM[i_PP_MM] = produceTH1D(Form("Mass_PM_Sub_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)),Form("Mass_PM_Sub_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)),"mass [GeV/cm^{2}]","Counts",750,0.,15.);
            histoPPMMFin[i_PP_MM] = produceTH1D(Form("Mass_PP_MM_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)),Form("Mass_PP_MM_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)),"mass [GeV/cm^{2}]","Counts",750,0.,15.);
            histoEMPMCheck[i_PP_MM] = produceTH1D(Form("Mass_PM_EM_Check_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)),Form("Mass_PM_EM_Check_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)),"mass [GeV/cm^{2}]","Counts",750,0.,15.);
            for (int nBins = 0; nBins < nBinPPMM; nBins ++){
                double MM_count = histoMMcuts[i_PP_MM]->GetBinContent(nBins+1);
                double PP_count = histoPPcuts[i_PP_MM]->GetBinContent(nBins+1);
                histoPPMMFin[i_PP_MM]->SetBinContent(nBins+1, 2 * sqrt(MM_count * PP_count));
                /*binContentMM.push_back(histoMMcuts[i_PP_MM]->GetBinContent(nBins));
                binContentPP.push_back(histoPPcuts[i_PP_MM]->GetBinContent(nBins));
                binContentPPMM.push_back(binContentPP.at(nBins)*binContentMM.at(nBins));
                cout << "binContentPPMM: " << binContentPPMM.at(nBins) << endl;
                binContentPPMMFin.push_back(2*sqrt(binContentPPMM.at(nBins)));
                cout << "binContentPPMMFin: " << binContentPPMMFin.at(nBins) << endl;
                histoPPMMFin[i_PP_MM]->SetBinContent(index, binContentPPMMFin.at(nBins));
                index = index +1;*/
            }
            nBinPMEM = histoPPMMFin[i_PP_MM]->GetNbinsX();
            for (int nBins = 0; nBins < nBinPMEM; nBins ++){
                double PM_count_Check = histoPPMMFin[i_PP_MM]->GetBinContent(nBins+1);
                double EM_count_Check = histoEvMix[i_PP_MM]->GetBinContent(nBins+1);
                if (PM_count_Check > 0){
                    histoEMPMCheck[i_PP_MM]->SetBinContent(nBins+1,EM_count_Check/PM_count_Check);
                }
                else {
                    continue;
                }
            }
            cPMEMCheck[i_PP_MM] = produceCanvas(false,false,Form("Mass_PM_EM_Check_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            cPMEMCheck[i_PP_MM]->cd();
            cPMEMCheck[i_PP_MM]->Draw();
            //histoPPMMFin[i_PP_MM]->SetName(Form("Mass_PP_MM_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            histoEMPMCheck[i_PP_MM]->Draw();
            histoEMPMCheck[i_PP_MM]->Write(0,2,0);
            //cout << "num_Histos: " << i_PP_MM << endl;
            cPPMM2[i_PP_MM] = produceCanvas(false,false,Form("Mass_PP_MM_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            cPPMM2[i_PP_MM]->cd();
            cPPMM2[i_PP_MM]->Draw();
            //histoPPMMFin[i_PP_MM]->SetName(Form("Mass_PP_MM_%d_%d",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            histoPPMMFin[i_PP_MM]->Draw();
            histoPPMMFin[i_PP_MM]->Write(0,2,0);
            //cout << "bins histo ppmm: " << histoPPMMFin[i_PP_MM]->GetNbinsX() << endl;
            binContentMM.clear();
            binContentPP.clear();
            binContentPPMM.clear();
            binContentPPMMFin.clear();
        }
        double tailLowGaus;
        double tailUpGaus;
        for (int i_PP_MM = 0; i_PP_MM < 4; i_PP_MM++){
            histoEvSame[i_PP_MM] = (TH1D*)centralityCuts->Get(Form("mass_%i_%i",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            histoEvSame[i_PP_MM]->SetDirectory(0);
            binSBLowLowSE.push_back(histoEvSame[i_PP_MM]->FindBin(2.));
            //cout << "binSBLowLowSE: " << binSBLowLowSE.at(i_PP_MM) << endl;
            binSBLowUpSE.push_back(histoEvSame[i_PP_MM]->FindBin(2.4));
            //cout << "binSBLowUpSE: " << binSBLowUpSE.at(i_PP_MM) << endl;
            integralSBLowSE.push_back(histoEvSame[i_PP_MM]->Integral(binSBLowLowSE.at(i_PP_MM), binSBLowUpSE.at(i_PP_MM)));
            binSBUpLowSE.push_back(histoEvSame[i_PP_MM]->FindBin(4.5));
            //cout << "binSBUpLowSE: " << binSBUpLowSE.at(i_PP_MM) << endl;
            binSBUpUpSE.push_back(histoEvSame[i_PP_MM]->FindBin(5.));
            //cout << "binSBUpUpSE: " << binSBUpUpSE.at(i_PP_MM) << endl;
            integralSBUpSE.push_back(histoEvSame[i_PP_MM]->Integral(binSBUpLowSE.at(i_PP_MM), binSBUpUpSE.at(i_PP_MM)));
            binSBLowLowME.push_back(histoEvMix[i_PP_MM]->FindBin(2.));
            binSBLowUpME.push_back(histoEvMix[i_PP_MM]->FindBin(2.4));
            integralSBLowME.push_back(histoEvMix[i_PP_MM]->Integral(binSBLowLowME.at(i_PP_MM), binSBLowUpME.at(i_PP_MM)));
            binSBUpLowME.push_back(histoEvMix[i_PP_MM]->FindBin(4.5));
            binSBUpUpME.push_back(histoEvMix[i_PP_MM]->FindBin(5.));
            integralSBUpME.push_back(histoEvMix[i_PP_MM]->Integral(binSBUpLowME.at(i_PP_MM), binSBUpUpME.at(i_PP_MM)));
            integralSumME.push_back(integralSBLowME.at(i_PP_MM)+integralSBUpME.at(i_PP_MM));
            integralSumSE.push_back(integralSBLowSE.at(i_PP_MM)+integralSBUpSE.at(i_PP_MM));
            ratio_integral_SE_ME.push_back(integralSumME.at(i_PP_MM)/integralSumSE.at(i_PP_MM));
            //cout << "Finisco prima fase ratio_integral_SE_ME: " << ratio_integral_SE_ME.at(i_PP_MM) << endl;
            //PPMM
            binSBLowLowPPMM.push_back(histoPPMMFin[i_PP_MM]->FindBin(2.));
            binSBLowUpPPMM.push_back(histoPPMMFin[i_PP_MM]->FindBin(2.4));
            integralSBLowPPMM.push_back(histoPPMMFin[i_PP_MM]->Integral(binSBLowLowPPMM.at(i_PP_MM), binSBLowUpPPMM.at(i_PP_MM)));
            binSBUpLowPPMM.push_back(histoPPMMFin[i_PP_MM]->FindBin(4.5));
            binSBUpUpPPMM.push_back(histoPPMMFin[i_PP_MM]->FindBin(5.));
            integralSBUpPPMM.push_back(histoPPMMFin[i_PP_MM]->Integral(binSBUpLowPPMM.at(i_PP_MM), binSBUpUpSE.at(i_PP_MM)));
            integralSumPPMM.push_back(integralSBLowPPMM.at(i_PP_MM)+integralSBUpPPMM.at(i_PP_MM));
            ratio_integral_SE_PPMM.push_back(integralSumPPMM.at(i_PP_MM)/integralSumSE.at(i_PP_MM));
            //cout << "Finisco seconda fase ratio_integral_SE_PPMM: " << ratio_integral_SE_PPMM.at(i_PP_MM) << endl;
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cEvMixPPMM[i_PP_MM] = produceCanvas(false,false,Form("cEvMixPPMM__%i_%i",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            cEvMixPPMM[i_PP_MM]->cd();
            cEvMixPPMM[i_PP_MM]->Draw();
            histoEvMix[i_PP_MM]->SetTitle(Form("Mass_With_Event_Mixing_%i_%i",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            histoEvMix[i_PP_MM]->SetLineColor(kRed);
            histoEvMix[i_PP_MM]->Draw("HISTO");
            histoEvSame[i_PP_MM]->SetLineColor(kBlue);
            histoEvSame[i_PP_MM]->Draw("SAME");
            histoPPMMFin[i_PP_MM]->SetLineColor(kGreen);
            histoPPMMFin[i_PP_MM]->Draw("SAME");
            cEvMixPPMM[i_PP_MM]->SetName(Form("Mass_With_Event_Mixing_%i_%i",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            //cout << "Finisco terza fase" << endl;
            cEvMixPPMM[i_PP_MM]->Write(0,2,0);
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cEvNorm[i_PP_MM] = produceCanvas(false,false,Form("cEvNorm_%i_%i",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            cEvNorm[i_PP_MM]->cd();
            cEvNorm[i_PP_MM]->Draw();
            histoEvMix[i_PP_MM]->Scale(1./ratio_integral_SE_ME.at(i_PP_MM));
            histoEvMix[i_PP_MM]->Draw();
            histoEvSame[i_PP_MM]->SetLineColor(kBlue);
            histoEvSame[i_PP_MM]->Draw("SAME");
            histoPPMMFin[i_PP_MM]->Scale(1./ratio_integral_SE_PPMM.at(i_PP_MM));
            histoPPMMFin[i_PP_MM]->SetLineColor(kGreen);
            histoPPMMFin[i_PP_MM]->Draw("SAME");
            cEvNorm[i_PP_MM]->SetName(Form("Ev_Norm_%i_%i",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            //cout << "Finisco quarta fase" << endl;
            cEvNorm[i_PP_MM]->Write(0,2,0);
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cFinSEME[i_PP_MM] = produceCanvas(false,false,Form("cFinSEME__%i_%i",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            cFinSEME[i_PP_MM]->cd();
            cFinSEME[i_PP_MM]->Draw();
            int nBinsES = histoEvSame[i_PP_MM] -> GetNbinsX();
            double EM_count, SE_count_1 = 0;
            for (int nBins = 0; nBins < nBinsES; nBins++) {
                SE_count_1 = histoEvSame[i_PP_MM] -> GetBinContent(nBins+1);
                EM_count = histoEvMix[i_PP_MM] -> GetBinContent(nBins+1);
                histoSEME[i_PP_MM] -> SetBinContent(nBins+1, (SE_count_1 - EM_count));
                histoSEME[i_PP_MM]->SetBinError(nBins+1, (SE_count_1 - EM_count)/TMath::Sqrt(nBinsES));
                //histoSEME[i_PP_MM]->SetBinErrorOption(TH1::kPoisson);
                //cout << "EM_count: " << EM_count << "  SE_count_1: " << SE_count_1 << "Sub: " << SE_count_1 - EM_count << endl;
            }
            //histoSEME[i_PP_MM] = histoEvSame[i_PP_MM]->Clone(Form("histoSEME_%i_%i", cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            //histoSEME[i_PP_MM]->Add(histoEvMix[i_PP_MM],-1.);
            histoSEME[i_PP_MM]->SetLineColor(kBlue);
            tailLowGaus = 2.8;
            tailUpGaus = 3.3;
            TF1* Gaus = new TF1("func","gaus",tailLowGaus,tailUpGaus);
            Gaus->SetParameters(0,2);
            Gaus->SetParameters(1,3);
            Gaus->SetParNames("Norm_{1}","#mu_{1}","#sigma_{1}");
            Gaus->SetLineColor(kRed);
            histoSEME[i_PP_MM]->Fit(Gaus, "RM+");
            gStyle->SetOptFit(1);
            histoSEME[i_PP_MM]->Draw("E,HISTO");
            histoSEME[i_PP_MM]->SetName(Form("Mass_Sub_%i_%i_EvMix",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1))); 
            //cout << "Finisco quinta fase" << endl;
            histoSEME[i_PP_MM]->Write(0,2,0);
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            cFinSEPPMM[i_PP_MM] = produceCanvas(false,false,Form("cFinSEPPMM__%i_%i",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            cFinSEPPMM[i_PP_MM]->cd();
            cFinSEPPMM[i_PP_MM]->Draw();
            double PM_count, SE_count = 0;
            for (int nBins = 0; nBins < nBinsES; nBins++) {
                SE_count = histoEvSame[i_PP_MM] -> GetBinContent(nBins+1);
                PM_count = histoPPMMFin[i_PP_MM] -> GetBinContent(nBins+1);
                histoSEPPMM[i_PP_MM] -> SetBinContent(nBins+1, (SE_count - PM_count));
                //cout << "PM_count: " << PM_count << "  SE_count: " << SE_count << "Sub: " << SE_count - PM_count << endl;
            }
            //histoSEPPMM[i_PP_MM] = histoEvSame[i_PP_MM]->Clone(Form("histoSEPPMM_%i_%i", cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1)));
            //histoSEPPMM[i_PP_MM]->Add(histoPPMMFin[i_PP_MM],-1.);
            histoSEPPMM[i_PP_MM]->SetLineColor(kBlue);
            if (i_PP_MM == 0 || i_PP_MM == 1){
                tailLowGaus = 2.9;
                tailUpGaus = 3.2;
            }
            else {
                tailLowGaus = 2.8;
                tailUpGaus = 3.3;
            }
            histoSEPPMM[i_PP_MM]->Fit(Gaus, "RM+");
            gStyle->SetOptFit(1);
            histoSEPPMM[i_PP_MM]->Draw();
            histoSEPPMM[i_PP_MM]->SetName(Form("Mass_Sub_%i_%i_PPMM",cutsAll.at(i_PP_MM),cutsAll.at(i_PP_MM+1))); 
            //cout << "Finisco sesta fase" << endl;
            histoSEPPMM[i_PP_MM]->Write(0,2,0);
            //return;
        }
        vector <double> checkVectorEM;
        vector <double> checkVector;
        for (int iCheck = 0; iCheck < histoSEME[0]->GetNbinsX(); iCheck++){
            checkVectorEM.push_back(histoEvSame[0]->GetBinContent(iCheck+1)-histoEvMix[0]->GetBinContent(iCheck+1));
            //cout << "binEM: " << iCheck+1 << "   checkVectorEM: " << checkVectorEM.at(iCheck) << endl;
            checkVector.push_back(histoEvSame[0]->GetBinContent(iCheck+1)-histoPPMMFin[0]->GetBinContent(iCheck+1));
            //cout << "binEMAdd: " << iCheck+1 << "   checkVectorEMAdd: " << histoSEME[0]->GetBinContent(iCheck+1) << endl;
            //cout << "binPM: " << iCheck+1 << "   checkVectorPM: " << checkVector.at(iCheck) << endl;
            //cout << "binPMAdd: " << iCheck+1 << "   checkVectorPMAdd: " << histoSEPPMM[0]->GetBinContent(iCheck+1) << endl;
        }
    }
    centralityCuts->Close();    

}  


void cutCentrality_20_40_60_90(){
    vector <float> binsMass;
    vector <float> binsCentrality;
    TLegend *legendCanva_60_90 = new TLegend(0.22,0.81,0.40,0.91,"","rNDC");
    legendCanva_60_90->SetBorderSize(0);
    legendCanva_60_90->SetFillStyle(0);
    legendCanva_60_90->SetTextFont(62);
    legendCanva_60_90->SetTextSize(0.04);
    TLegend *legendCanva_40_60 = new TLegend(0.22,0.81,0.40,0.91,"","rNDC");
    legendCanva_40_60->SetBorderSize(0);
    legendCanva_40_60->SetFillStyle(0);
    legendCanva_40_60->SetTextFont(62);
    legendCanva_40_60->SetTextSize(0.04);
    TLegend *legendCanva_20_40 = new TLegend(0.22,0.81,0.40,0.91,"","rNDC");
    legendCanva_20_40->SetBorderSize(0);
    legendCanva_20_40->SetFillStyle(0);
    legendCanva_20_40->SetTextFont(62);
    legendCanva_20_40->SetTextSize(0.04);
    TLegend *legendCanva_0_20 = new TLegend(0.22,0.81,0.40,0.91,"","rNDC");
    legendCanva_0_20->SetBorderSize(0);
    legendCanva_0_20->SetFillStyle(0);
    legendCanva_0_20->SetTextFont(62);
    legendCanva_0_20->SetTextSize(0.04);
    TFile *centralityCuts = TFile::Open("centralityCuts.root","UPDATE");
    TH1D *histo_60_90 = (TH1D*)centralityCuts->Get("mass_60_90");
    TH1D *histo_40_60 = (TH1D*)centralityCuts->Get("mass_40_60");
    TH1D *histo_20_40 = (TH1D*)centralityCuts->Get("mass_20_40");
    TH1D *histo_0_20 = (TH1D*)centralityCuts->Get("mass_0_20");
    TCanvas *canva_60_90 = (TCanvas*)centralityCuts->Get(Form("Ev_Norm_60_90"));
    TList *list_60_90 = canva_60_90->GetListOfPrimitives();
    TH1D *histo_60_90_EM = (TH1D*)list_60_90->FindObject("Mass_CentFT0C_px");
    TH1D *histo_60_90_SE = (TH1D*)list_60_90->FindObject("mass_60_90");
    TH1D *histo_60_90_PPMM = (TH1D*)list_60_90->FindObject("Mass_PP_MM_60_90");
    TCanvas *canva_40_60 = (TCanvas*)centralityCuts->Get("Ev_Norm_40_60");
    TList *list_40_60 = canva_40_60->GetListOfPrimitives();
    TH1D *histo_40_60_EM = (TH1D*)list_40_60->FindObject("Mass_CentFT0C_px");
    TH1D *histo_40_60_SE = (TH1D*)list_40_60->FindObject("mass_40_60");
    TH1D *histo_40_60_PPMM = (TH1D*)list_40_60->FindObject("Mass_PP_MM_40_60");
    TCanvas *canva_20_40 = (TCanvas*)centralityCuts->Get("Ev_Norm_20_40");
    TList *list_20_40 = canva_20_40->GetListOfPrimitives();
    TH1D *histo_20_40_EM = (TH1D*)list_20_40->FindObject("Mass_CentFT0C_px");
    TH1D *histo_20_40_SE = (TH1D*)list_20_40->FindObject("mass_20_40");
    TH1D *histo_20_40_PPMM = (TH1D*)list_20_40->FindObject("Mass_PP_MM_20_40");
    TCanvas *canva_0_20 = (TCanvas*)centralityCuts->Get("Ev_Norm_0_20");
    TList *list_0_20 = canva_0_20->GetListOfPrimitives();
    TH1D *histo_0_20_EM = (TH1D*)list_0_20->FindObject("Mass_CentFT0C_px");
    TH1D *histo_0_20_SE = (TH1D*)list_0_20->FindObject("mass_0_20");
    TH1D *histo_0_20_PPMM = (TH1D*)list_0_20->FindObject("Mass_PP_MM_0_20");
    TCanvas *cSubAll = new TCanvas("cSubAll","SubAll",200,10,700,500);
    cSubAll->Divide(2,2);
    cSubAll->cd(4);
    histo_60_90->SetLineColor(kBlue);
    histo_60_90_EM->Draw("SAME");
    histo_60_90_SE->Draw("SAME");
    histo_60_90_PPMM->Draw("SAME");
    legendCanva_60_90->AddEntry(histo_60_90_EM,"Event mixing","lep");
    legendCanva_60_90->AddEntry(histo_60_90_SE,"Event pairing","lep");
    legendCanva_60_90->AddEntry(histo_60_90_PPMM,"Same sign events","lep");
    legendCanva_60_90->Draw("SAME");
    histo_60_90_EM->SetTitle("Mass_Norm_(60-90)%%");
    histo_60_90_EM->GetXaxis()->SetRangeUser(2.5,4.5);
    histo_60_90_SE->GetXaxis()->SetRangeUser(2.5,4.5);
    histo_60_90_PPMM->GetXaxis()->SetRangeUser(2.5,4.5);
    cSubAll->Update();
    cSubAll->cd(3);
    histo_40_60->SetLineColor(kBlue);
    histo_40_60_EM->Draw("SAME");
    histo_40_60_SE->Draw("SAME");
    histo_40_60_PPMM->Draw("SAME");
    legendCanva_40_60->AddEntry(histo_40_60_EM,"Event mixing","lep");
    legendCanva_40_60->AddEntry(histo_40_60_SE,"Event pairing","lep");
    legendCanva_40_60->AddEntry(histo_40_60_PPMM,"Same sign events","lep");
    legendCanva_40_60->Draw("SAME");
    histo_40_60_EM->SetTitle("Mass_Norm_(40-60)%%");
    histo_40_60_EM->GetXaxis()->SetRangeUser(2.5,4.5);
    histo_40_60_SE->GetXaxis()->SetRangeUser(2.5,4.5);
    histo_40_60_PPMM->GetXaxis()->SetRangeUser(2.5,4.5);
    cSubAll->Update();
    cSubAll->cd(2);
    histo_20_40->SetLineColor(kBlue);
    histo_20_40_EM->Draw("SAME");
    histo_20_40_SE->Draw("SAME");
    histo_20_40_PPMM->Draw("SAME");
    legendCanva_20_40->AddEntry(histo_20_40_EM,"Event mixing","lep");
    legendCanva_20_40->AddEntry(histo_20_40_SE,"Event pairing","lep");
    legendCanva_20_40->AddEntry(histo_20_40_PPMM,"Same sign events","lep");
    legendCanva_20_40->Draw("SAME");
    histo_20_40_EM->SetTitle("Mass_Norm_(20-40)%%");
    histo_20_40_EM->GetXaxis()->SetRangeUser(2.5,4.5);
    histo_20_40_SE->GetXaxis()->SetRangeUser(2.5,4.5);
    histo_20_40_PPMM->GetXaxis()->SetRangeUser(2.5,4.5);
    cSubAll->Update();
    cSubAll->cd(1);
    histo_0_20->SetLineColor(kBlue);
    histo_0_20_EM->Draw("SAME");
    histo_0_20_SE->Draw("SAME");
    histo_0_20_PPMM->Draw("SAME");
    legendCanva_0_20->AddEntry(histo_0_20_EM,"Event mixing","lep");
    legendCanva_0_20->AddEntry(histo_0_20_SE,"Event pairing","lep");
    legendCanva_0_20->AddEntry(histo_0_20_PPMM,"Same sign events","lep");
    legendCanva_0_20->Draw("SAME");
    histo_0_20_EM->SetTitle("Mass_Norm_(0-20)%%");
    histo_0_20_EM->GetXaxis()->SetRangeUser(2.5,4.5);
    histo_0_20_SE->GetXaxis()->SetRangeUser(2.5,4.5);
    histo_0_20_PPMM->GetXaxis()->SetRangeUser(2.5,4.5);
    cSubAll->Update();
    cSubAll->Draw();
    cSubAll->SetName("Sub All");
    cSubAll->Write(0,2,0);
    centralityCuts->Close();
}