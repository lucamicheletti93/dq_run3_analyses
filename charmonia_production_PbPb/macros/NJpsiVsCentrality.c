#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMinuit.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMath.h"
#include <fstream>
#include <TStyle.h>
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TTree.h"
#include "TBranch.h"
#include <iostream>
#include <vector>
#include <stdlib.h> 
#include <string>
#include <iterator>
#include <numeric>
#include <map>
#include "THStack.h"
#include <list>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>
#include <string_view>
#include "Fit/FitResult.h"
#include "TEllipse.h"
#include <tuple>
#include <TLegend.h>

using namespace std;
TCanvas *produceCanvas(bool logX, bool logY) {
	TCanvas *c = new TCanvas();
	c->SetFrameFillColor(0); //Transparent background
    c->SetFrameFillStyle(0);
    c->SetFrameBorderMode(0);
	c->SetGridx(); //Grid on x
	c->SetGridy(); //Grid on y 
	if (logX) c->SetLogx(); //Log on x, if needed
	if (logY) c->SetLogy(); //Log on y, if needed
	return c;
}
TH1F *produceTH1(string name, string title, string xAxisTitle, string yAxisTitle, int nBins, double xMin, double xMax) {
	TH1F *h = new TH1F(name.c_str(), title.c_str(), nBins, xMin, xMax);
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

void NJpsiVsCentrality(){
    vector <float> parValArray;
    vector <float> parErrArray;
    vector <float> parValArrayAll;
    vector <float> parErrArrayAll;
    float errorJPsiTotAll;
    vector <int> run2Data = {108070, 69808, 45007, 24876, 15010, 7895, 4112, 2042, 932};
    vector <float> run2Err = {1801, 1168, 822, 491, 321, 190, 105, 66, 37};
    vector <float> ratioRun2Data;
    vector <float> ratioRun2Err;
    vector <float> ratioRun3Data;
    vector <float> ratioRun3Err;
    float nJPsiRun3Tot = 418376.3631;
    float errJPsiRun3Tot = 8026.1089;
    float nJPsiRun2Tot = 277007;
    float errJPsiRun2Tot = 8026.1089;
    vector <int> cuts = {10,20,30};
    vector <int> cutsAll = {0,10,20,30,40,50,60,70,80,90};
    const Int_t num = 3;
    char cut1[10][100] = {"(10-70)%", "(20-70)%", "(30-70)%"};
    char cut1All[20][100] = {"(0-10)%", "(10-20)%", "(20-30)%","(30-40)%", "(40-50)%", "(50-60)%", "(60-70)%", "(70-80)%", "(80-90)%", "Sum", "Integrated"};
    char cutsRuns[20][100] = {"(0-10)%", "(10-20)%", "(20-30)%","(30-40)%", "(40-50)%", "(50-60)%", "(60-70)%", "(70-80)%", "(80-90)%"};
    const char *cutName1 = "(10-70)";
    const char *cutName2 = "(20-70)";
    const char *cutName3 = "(30-70)";
    TH1F *histo[3];
    TH1F *histo2[3];
    TH1F *histoPlot[3];
    TH1F *histoPlot2[3];
    TH1F *histoPlotAll[9];
    TH1F *histoPlotAll2[9];
    TCanvas *canvas[3];
    TCanvas *canvasPlot[3];
    TCanvas *canvasPlotAll[9];
    TCanvas *c10 = new TCanvas("c10","Mass_for_3_centrality_cuts",200,10,700,500);
    c10->SetFrameFillColor(0); //Transparent background
    c10->SetFrameFillStyle(0);
    c10->SetFrameBorderMode(0);
	c10->SetGridx(); //Grid on x
	c10->SetGridy(); //Grid on y 
    TCanvas *c11 = new TCanvas("c11", "c11", 200,10,700,500);
    TCanvas *prova = produceCanvas(false,false);
    TCanvas *provaAll = produceCanvas(false,false);
    TCanvas *run2Ratio = produceCanvas(false,false);
    TCanvas *run3Ratio = produceCanvas(false,false);
    TCanvas *c12 = new TCanvas("c12","Mass_for_all_cuts",200,10,700,500);
    c12->SetFrameFillColor(0); //Transparent background
    c12->SetFrameFillStyle(0);
    c12->SetFrameBorderMode(0);
	c12->SetGridx(); //Grid on x
	c12->SetGridy(); //Grid on y 
    TH1F *histoTotJPsi = new TH1F("histoTotJPsi", "Number of J/psi for each centrality cut", 3., 15000., 80000.);
    histoTotJPsi->GetXaxis()->SetTitle("Centrality cuts");
	histoTotJPsi->GetYaxis()->SetTitle("J/psi counts");
	histoTotJPsi->GetXaxis()->SetTitleFont(62);
	histoTotJPsi->GetXaxis()->SetLabelFont(62);
	histoTotJPsi->GetYaxis()->SetTitleFont(62);
	histoTotJPsi->GetYaxis()->SetLabelFont(62);
	histoTotJPsi->GetXaxis()->CenterTitle(true);
	histoTotJPsi->GetYaxis()->CenterTitle(true);
    TH1F *histoTotJPsiAllCuts = new TH1F("histoTotJPsiAllCuts", "Number of J/psi for all centrality cuts", 11., 15000., 1000000.);
    histoTotJPsiAllCuts->GetXaxis()->SetTitle("Centrality cuts");
	histoTotJPsiAllCuts->GetYaxis()->SetTitle("J/psi counts");
	histoTotJPsiAllCuts->GetXaxis()->SetTitleFont(62);
	histoTotJPsiAllCuts->GetXaxis()->SetLabelFont(62);
	histoTotJPsiAllCuts->GetYaxis()->SetTitleFont(62);
	histoTotJPsiAllCuts->GetYaxis()->SetLabelFont(62);
	histoTotJPsiAllCuts->GetXaxis()->CenterTitle(true);
	histoTotJPsiAllCuts->GetYaxis()->CenterTitle(true);
    TH1F *histoRatioRun3 = new TH1F("histoRatioRun3", "Ratio of J/psi Run 3", 9., 15000., 1000000.);
    histoRatioRun3->GetXaxis()->SetTitle("Centrality cuts");
	histoRatioRun3->GetYaxis()->SetTitle("(J/psi cut)/(J/psi tot)");
	histoRatioRun3->GetXaxis()->SetTitleFont(62);
	histoRatioRun3->GetXaxis()->SetLabelFont(62);
	histoRatioRun3->GetYaxis()->SetTitleFont(62);
	histoRatioRun3->GetYaxis()->SetLabelFont(62);
	histoRatioRun3->GetXaxis()->CenterTitle(true);
	histoRatioRun3->GetYaxis()->CenterTitle(true);
    TH1F *histoRatioRun2 = new TH1F("histoRatioRun2", "Ratio of J/psi Run 2 and Run 3", 9., 15000., 1000000.);
    histoRatioRun2->GetXaxis()->SetTitle("Centrality cuts");
	histoRatioRun2->GetYaxis()->SetTitle("(J/psi cut)/(J/psi tot)");
	histoRatioRun2->GetXaxis()->SetTitleFont(62);
	histoRatioRun2->GetXaxis()->SetLabelFont(62);
	histoRatioRun2->GetYaxis()->SetTitleFont(62);
	histoRatioRun2->GetYaxis()->SetLabelFont(62);
	histoRatioRun2->GetXaxis()->CenterTitle(true);
	histoRatioRun2->GetYaxis()->CenterTitle(true);
    //gStyle->SetHistFillColor(kBlue);
	//gStyle->SetHistFillStyle(kRed);
	//gStyle->SetHistLineColor(kBlue);
    TLegend *legend = new TLegend(0.22,0.81,0.40,0.91,"","rNDC");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(62);
    legend->SetTextSize(0.04);
    for (Int_t i = 0; i < cuts.size(); i++){
        TFile* analysisResults = TFile::Open(("output_test_"+to_string(cuts.at(i))+"/output__CB2_VWG__2.5_4.5.root").c_str());
        cout << "output_test_"+to_string(cuts.at(i))+"/output__CB2_VWG__2.5_4.5.root" << endl;
        canvas[i] = new TCanvas();
        //histoPlot[i] = produceTH1("prova","prova","Time [ns]","Counts",100.,2.5,4.5);
        histo[i] = (TH1F*)analysisResults->Get("fit_results_CB2_VWG__2.5_4.5");
        histo2[i] = (TH1F*)histo[i]->Clone();
        histo2[i]->SetDirectory(0);
        histo2[i]->SetTitle(("output_test_"+to_string(cuts.at(i))).c_str());
        TCanvas *c= (TCanvas *)analysisResults->Get("fit_plot_CB2_VWG__2.5_4.5");
        TList *l=c->GetListOfPrimitives();
        canvasPlot[i] = new TCanvas();
        histoPlot[i] = (TH1F*)l->At(3);
        histoPlot2[i] = (TH1F*)histoPlot[i]->Clone();
        histoPlot2[i]->SetDirectory(0);
        parValArray.push_back(histo2[i]->GetBinContent(10));
        parErrArray.push_back(histo2[i]->GetBinError(10));
        canvas[i]->cd();
        canvas[i]->Draw();
        histo2[i]->Draw("HISTO");
    }
    for (Int_t i_bins = 0; i_bins < 9; i_bins++){
        TFile* analysisResults = TFile::Open(("output_"+to_string(cutsAll.at(i_bins))+"_"+to_string(cutsAll.at(i_bins+1))+"/output__CB2_VWG__2.5_4.5.root").c_str());
        cout << ("output_"+to_string(cutsAll.at(i_bins))+"_"+to_string(cutsAll.at(i_bins+1))+"/output__CB2_VWG__2.5_4.5.root") << endl;
        canvasPlotAll[i_bins] = new TCanvas();
        histoPlotAll[i_bins] = (TH1F*)analysisResults->Get("fit_results_CB2_VWG__2.5_4.5");
        //TCanvas *cPlot= (TCanvas *)analysisResults->Get("fit_plot_CB2_VWG__2.5_4.5");
        //TList *listPlot=cPlot->GetListOfPrimitives();
        //histoPlotAll[i_bins] = (TH1F*)listPlot->At(3);
        histoPlotAll2[i_bins] = (TH1F*)histoPlotAll[i_bins]->Clone();
        histoPlotAll2[i_bins]->SetDirectory(0);
        histoPlotAll2[i_bins]->SetTitle(("output_"+to_string(cutsAll.at(i_bins))+"_"+to_string(cutsAll.at(i_bins+1))).c_str());
        parValArrayAll.push_back(histoPlotAll2[i_bins]->GetBinContent(10));
        cout << "Numero J/Psi per cut: " << parValArrayAll.at(i_bins) << "cut: " << i_bins  << endl;
        parErrArrayAll.push_back(histoPlotAll2[i_bins]->GetBinError(10));
        canvasPlotAll[i_bins]->cd();
        canvasPlotAll[i_bins]->Draw();
        histoPlotAll2[i_bins]->Draw("HISTO");
        ratioRun3Data.push_back(parValArrayAll.at(i_bins)/nJPsiRun3Tot);
        cout << "Numero J/Psi per cut: " << parValArrayAll.at(i_bins) << "cut: " << i_bins  << endl;
        ratioRun2Data.push_back(run2Data.at(i_bins)/nJPsiRun2Tot);
        ratioRun3Err.push_back(sqrt(pow((1/nJPsiRun3Tot),2)*pow(errJPsiRun3Tot,2)+(pow((parValArrayAll.at(i_bins)/pow(nJPsiRun3Tot,2)),2)*pow(parErrArrayAll.at(i_bins),2))));
        ratioRun2Err.push_back(sqrt(pow((1/nJPsiRun2Tot),2)*pow(errJPsiRun2Tot,2)+(pow((run2Data.at(i_bins)/pow(nJPsiRun2Tot,2)),2)*pow(run2Err.at(i_bins),2))));
    }
    for (int k = 0; k < 9; k++){
        cout << "Errori run 3 " << ratioRun3Err.at(k) << endl;
    }
    for (int k = 0; k < 9; k++){
        cout << "Errori run 2 " << ratioRun2Err.at(k) << endl;
    }
    float mean = 0;
    for (Int_t i_val = 0; i_val < parValArray.size(); i_val++){
        mean += parValArray.at(i_val);
    }
    mean = mean/parValArray.size();
    float stdDev = 0;
    for (Int_t i_val = 0; i_val < parValArray.size(); i_val++){
        stdDev += (parValArray.at(i_val) - mean) * (parValArray.at(i_val) - mean);
    }
    stdDev = sqrt(stdDev / parValArray.size());
    errorJPsiTotAll = stdDev/sqrt(parValArray.size());
    cout << "errorJPsiTotAll: " << errorJPsiTotAll << endl;
    float sumValAll = parValArrayAll.at(0);
    for (Int_t i_val = 0; i_val < (parValArrayAll.size()-1); i_val++){
        sumValAll = sumValAll + parValArrayAll.at(i_val+1);
    }
    float powErrJPsiTotAll = pow(parErrArrayAll.at(0),2);
    for (Int_t i_err = 0; i_err < (parErrArrayAll.size()-1); i_err++){
        powErrJPsiTotAll = powErrJPsiTotAll + pow(parErrArrayAll.at(i_err+1),2);
    }
    //errorJPsiTotAll = sqrt(powErrJPsiTotAll);
    parValArrayAll.push_back(sumValAll);
    parErrArrayAll.push_back(errorJPsiTotAll);
    parValArrayAll.push_back(nJPsiRun3Tot);
    parErrArrayAll.push_back(errJPsiRun3Tot);
    //histoTotJPsi->GetXaxis()->SetBinLabel(0, "(10-70)");
    int index = 1;
    for (Int_t i = 0; i < parValArray.size(); i++){
        //histoTotJPsi->Fill(parValArray.at(i));
        histoTotJPsi->GetXaxis()->SetBinLabel(index, cut1[i]);
        histoTotJPsi->SetBinContent(index, parValArray.at(i));
        histoTotJPsi->SetBinError(index, parErrArray.at(i));
        index = index + 1;
    }
    int indexAll = 1;
    for (Int_t i_bins = 0; i_bins < parValArrayAll.size(); i_bins++){
        //histoTotJPsi->Fill(parValArrayAll.at(i));
        histoTotJPsiAllCuts->GetXaxis()->SetBinLabel(indexAll, cut1All[i_bins]);
        histoTotJPsiAllCuts->SetBinContent(indexAll, parValArrayAll.at(i_bins));
        histoTotJPsiAllCuts->SetBinError(indexAll, parErrArrayAll.at(i_bins));
        indexAll = indexAll + 1;
    }
    int indexRun3 = 1;
    for (Int_t binRun3 = 0; binRun3 < 9; binRun3++){
        //histoTotJPsi->Fill(parValArrayAll.at(i));
        histoRatioRun3->GetXaxis()->SetBinLabel(indexRun3, cutsRuns[binRun3]);
        histoRatioRun3->SetBinContent(indexRun3, ratioRun3Data.at(binRun3));
        histoRatioRun3->SetBinError(indexRun3, ratioRun3Err.at(binRun3));
        indexRun3 = indexRun3 + 1;
    }
    histoRatioRun3->SetLineColor(kBlue);
    int indexRun2 = 1;
    for (Int_t binRun2 = 0; binRun2 < 9; binRun2++){
        //histoTotJPsi->Fill(parValArrayAll.at(i));
        histoRatioRun2->GetXaxis()->SetBinLabel(indexRun2, cutsRuns[binRun2]);
        histoRatioRun2->SetBinContent(indexRun2, ratioRun2Data.at(binRun2));
        histoRatioRun2->SetBinError(indexRun2, ratioRun2Err.at(binRun2));
        indexRun2 = indexRun2 + 1;
    }
    histoRatioRun2->SetLineColor(kRed);
    legend->AddEntry(histoRatioRun2,"Run 2","lep");
    legend->AddEntry(histoRatioRun3,"Run 3","lep");
    TFile *fout = new TFile("Mean_J_psi.root","RECREATE");
    fout->cd();
    prova->cd();
    //histoTotJPsi->SetTitle("histo JPsi tot")
    histoTotJPsi->Draw("E,HISTO");
    prova->Write("Total_J_psi");
    /*for(int i = 0; i < 3; i++){
        canvasPlot[i]->cd();
        cout << "33" << endl;
        canvasPlot[i]->Draw();
        histoPlot2[i]->Draw();
        //canvasPlot[i]->Write(("output_plot_"+to_string(cuts.at(i))).c_str());
    }*/
    provaAll->cd();
    histoTotJPsiAllCuts->Draw("E,HISTO");
    histoTotJPsiAllCuts->Write("Total_J_psi_all_cuts");
    run2Ratio->cd();
    histoRatioRun2->Draw("E,HISTO");
    //histoRatioRun2->Write("Ratio_Run_2_Run_3");
    //run3Ratio->cd();
    histoRatioRun3->Draw("E,HISTO,SAME");
    legend->Draw("SAME");
    run2Ratio->Write("Ratio_Run_2_Run_3");
    //histoRatioRun3->Write("Ratio_Run_3");
    c11->cd();
    c11->Draw();
    histoPlot2[0]->SetMarkerStyle(22);
    histoPlot2[0]->SetMarkerColor(kRed);
    histoPlot2[0]->Draw();
    histoPlot2[1]->SetMarkerStyle(21);
    histoPlot2[1]->SetMarkerColor(kBlue);
    histoPlot2[1]->Draw("SAME");
    histoPlot2[2]->SetMarkerStyle(20);
    histoPlot2[2]->SetMarkerColor(kGreen);
    histoPlot2[2]->Draw("SAME");  
    c11->Write("output_plot");
    c12->cd();
    c12->Draw();
    histoPlotAll2[0]->SetMarkerStyle(22);
    histoPlotAll2[0]->SetMarkerColor(kRed);
    histoPlotAll2[0]->Draw();
    for (Int_t i_histos = 0; i_histos < 9; i_histos++){
        histoPlotAll2[i_histos]->SetMarkerStyle(22+i_histos);
        histoPlotAll2[i_histos]->SetMarkerColor(kRed+i_histos);
        histoPlotAll2[i_histos]->Draw("SAME");
    }
    c12->Write("output_plot_all");
    fout->Close();
    /////////////////////////////////////
    //c10->cd();
    //c10->Draw();
    //legend->AddEntry(histoPlot2[0],"(10-70)%","p");
    //legend->AddEntry(histoPlot2[1],"(20-70)%","p");
    //legend->AddEntry(histoPlot2[2],"(30-70)%","p");
    //hs->Draw();
    //legend->Draw("SAME");
    //c10->Write("Mass_for_3_centrality_cuts");
    //TCanvas *c1 = new TCanvas("histoResults");
    //histo->Draw();
}