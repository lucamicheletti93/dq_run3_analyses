
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"


void LoadStyle();

void Compare_Ophelie(){

 LoadStyle();


  
//----------------------------------------------------------
// ALICE RAA fw-y 2.5<y<4, 2015 0-90%. sqrt(s_NN) = 5.02 TeV (2015)           // OPHELIE
// published fig 2
// https://arxiv.org/pdf/1909.03158 JHEP 2002 (2020) 041
//----------------------------------------------------------

const int n_ALICE_pub_fw090=12;
Double_t pt_ALICE_pub_fw090[n_ALICE_pub_fw090]={0.2,0.475,0.825,1.5,2.5,3.5,4.5,5.5,6.5,7.5,9.,11.};
Double_t ept_ALICE_pub_fw090[n_ALICE_pub_fw090]={0.1,0.175,0.175,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1.,1.};
Double_t esyst_pt_ALICE_pub_fw090[n_ALICE_pub_fw090]={0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2};
Double_t RAA_ALICE_pub_fw090[n_ALICE_pub_fw090]={10.28,1.33,0.98,0.85,0.79,0.79,0.81,0.74,0.73,0.80,0.75,0.82};
Double_t Stat_RAA_ALICE_pub_fw090[n_ALICE_pub_fw090]={0.59,0.10,0.08,0.03,0.03,0.04,0.04,0.05,0.06,0.09,0.10,0.15};
Double_t Syst_RAA_ALICE_pub_fw090[n_ALICE_pub_fw090]={0.78,0.08,0.06,0.06,0.05,0.05,0.05,0.05,0.04,0.05,0.05,0.06};
Double_t Glob_RAA_ALICE_pub_fw = 0.07587;

 TGraphErrors *gRAA_ALICE_pub_fw090 = new TGraphErrors(n_ALICE_pub_fw090,pt_ALICE_pub_fw090,RAA_ALICE_pub_fw090,ept_ALICE_pub_fw090,Stat_RAA_ALICE_pub_fw090);
  gRAA_ALICE_pub_fw090->SetMarkerColor(kBlack);
  gRAA_ALICE_pub_fw090->SetLineColor(kBlack);
  gRAA_ALICE_pub_fw090->SetMarkerStyle(22);
  gRAA_ALICE_pub_fw090->SetMarkerSize(1.5);
  
  TGraphErrors *gRAA_ALICE_pub_fw090_syst = new TGraphErrors(n_ALICE_pub_fw090,pt_ALICE_pub_fw090,RAA_ALICE_pub_fw090,esyst_pt_ALICE_pub_fw090,Syst_RAA_ALICE_pub_fw090);
  gRAA_ALICE_pub_fw090_syst->SetMarkerColor(kBlack);
  gRAA_ALICE_pub_fw090_syst->SetMarkerStyle(22);
  //RAA_ALICE_pub_fw090_syst->SetFillColorAlpha(kGreen+3,0.3);
  gRAA_ALICE_pub_fw090_syst->SetLineColor(kBlack);
  gRAA_ALICE_pub_fw090_syst->SetMarkerSize(1.5);

  TGraphErrors *gRAA_ALICE_pub_fw090_syst2 = new TGraphErrors(n_ALICE_pub_fw090,pt_ALICE_pub_fw090,RAA_ALICE_pub_fw090,esyst_pt_ALICE_pub_fw090,Syst_RAA_ALICE_pub_fw090);
  gRAA_ALICE_pub_fw090_syst2->SetMarkerColor(kBlack);
  gRAA_ALICE_pub_fw090_syst2->SetMarkerStyle(22);
  gRAA_ALICE_pub_fw090_syst2->SetFillStyle(0);
  gRAA_ALICE_pub_fw090_syst2->SetLineColor(kBlack);
  gRAA_ALICE_pub_fw090_syst2->SetMarkerSize(1.5);

  //TBox *boxglobRAA_ALICE_pub_fw090 = new TBox(14.5,1.-Glob_RAA_ALICE_pub_fw,15,1+Glob_RAA_ALICE_pub_fw);     
  //boxglobRAA_ALICE_pub_fw090->SetFillColor(kBlack);

    TBox *boxglobRAA_ALICE_pub_fw090 = new TBox(12.,1.-Glob_RAA_ALICE_pub_fw,12.5,1+Glob_RAA_ALICE_pub_fw);     
  boxglobRAA_ALICE_pub_fw090->SetFillColor(kBlack);


//----------------------------------------------------------
// ALICE ROO fw-y sqrt(s_NN) = 5.36 TeV 2.5<y<4
// Figure to be approved as preliminary
// AN https://alice-notes.web.cern.ch/node/1743
//----------------------------------------------------------

const int n_ALICE_OO=7;
Double_t pt_ALICE_OO[n_ALICE_OO]={0.5, 1.5,2.5,3.5,4.5,5.5,7.};
Double_t ept_ALICE_OO[n_ALICE_OO]={0.5,0.5,0.5,0.5,0.5,0.5,1};
Double_t esyst_pt_ALICE_OO[n_ALICE_OO]={0.2,0.2,0.2,0.2,0.2,0.2,0.2};
Double_t ROO_ALICE_OO[n_ALICE_OO]={0.581728,0.629188,0.641791,0.656929,0.702026,0.698702,0.706779};
Double_t Stat_ROO_ALICE_OO[n_ALICE_OO]={0.014429,0.011421,0.012623,0.015387,0.021201,0.025506,0.024434};
Double_t Syst_ROO_ALICE_OO[n_ALICE_OO]={0.045811,0.041387,0.040849,0.045579,0.067156,0.058217,0.080002}; //0.034578,0.033524,0.032875,0.035657,0.049118,0.049178,0.047188
Double_t GlobError_RO_ALICE=0.050744457825461095; //0.06461423991660044

 TGraphErrors *gROO_ALICE_OO = new TGraphErrors(n_ALICE_OO,pt_ALICE_OO,ROO_ALICE_OO,ept_ALICE_OO,Stat_ROO_ALICE_OO);
  gROO_ALICE_OO->SetMarkerColor(kRed+1);
  gROO_ALICE_OO->SetLineColor(kBlack);
  gROO_ALICE_OO->SetMarkerStyle(20);
  gROO_ALICE_OO->SetMarkerSize(1.5);
  
  TGraphErrors *gROO_ALICE_syst = new TGraphErrors(n_ALICE_OO,pt_ALICE_OO,ROO_ALICE_OO,esyst_pt_ALICE_OO,Syst_ROO_ALICE_OO);
  gROO_ALICE_syst->SetMarkerColor(kRed+1);
  gROO_ALICE_syst->SetMarkerStyle(20);
  //gROO_ALICE_syst->SetFillColorAlpha(kRed+1,0.3);
  gROO_ALICE_syst->SetLineColor(kRed+1);
  gROO_ALICE_syst->SetMarkerSize(1.5);

  TGraphErrors *gROO_ALICE_syst2 = new TGraphErrors(n_ALICE_OO,pt_ALICE_OO,ROO_ALICE_OO,esyst_pt_ALICE_OO,Syst_ROO_ALICE_OO);
  gROO_ALICE_syst2->SetMarkerColor(kRed+1);
  gROO_ALICE_syst2->SetMarkerStyle(20);
  gROO_ALICE_syst2->SetFillStyle(0);
  gROO_ALICE_syst2->SetLineColor(kRed+1);
  gROO_ALICE_syst2->SetMarkerSize(1.5);

  TBox *boxglobROO_ALICE = new TBox(12.5,1.-GlobError_RO_ALICE,13.0,1+GlobError_RO_ALICE);     
  boxglobROO_ALICE->SetFillColor(kRed+1);
  

//--------------------------------------------------------------------------
// Plots RAA fw-y published results 0-90% 2015
//-------------------------------------------------------------------------- 

 TCanvas *c3 = new TCanvas("c3","c3",20,20,800,600);
 c3->SetFillColor(0);
 c3->SetBorderMode(0);
 c3->SetBorderSize(2);
 c3->SetLeftMargin(0.12);
 c3->SetRightMargin(0.03);
 c3->SetTopMargin(0.03);
 c3->SetBottomMargin(0.13);
 c3->SetFrameBorderMode(0);

 TH2D *hpt3 = new TH2D("hpt3","",100,0.,13.,100,0.,2.0);
 hpt3->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#kern[-1.4]{ }#it{c})");
 hpt3->GetYaxis()->SetTitle("#it{R}_{AA}");
 hpt3->Draw();


 gRAA_ALICE_pub_fw090_syst->Draw("E2P");
 gRAA_ALICE_pub_fw090_syst2->Draw("E2P");
 gRAA_ALICE_pub_fw090->Draw("P");

 gROO_ALICE_syst->Draw("E2P");
 gROO_ALICE_syst2->Draw("E2P");
 gROO_ALICE_OO->Draw("P");

 
 boxglobRAA_ALICE_pub_fw090->Draw();
 boxglobROO_ALICE->Draw();

 TLatex *latexTitle = new TLatex();
 latexTitle -> SetTextSize(0.050);
 latexTitle -> SetNDC();
 latexTitle -> SetTextFont(42);
  latexTitle -> DrawLatex(0.15, 0.90, "ALICE Preliminary, J/#kern[-2.1]{ }#psi #kern[-1.9]{ }#rightarrow #kern[-1.9]{ }#mu^{+}#mu^{-}");

 //latexTitle -> DrawLatex(0.15, 0.90, "ALICE Preliminary, J/#kern[-2.1]{ }#psi #rightarrow #mu^{+}#mu^{-}");
 
 TLegend *legpt4 = new TLegend(0.13,0.75,0.45,0.88);
 legpt4->SetBorderSize(0);
 legpt4->SetFillColor(10);
 legpt4->SetFillStyle(1);
 legpt4->SetLineStyle(0);
 legpt4->SetLineColor(0);
 legpt4->SetTextSize(0.04);
 legpt4->AddEntry(gROO_ALICE_OO,"OO #sqrt{#it{s}_{NN}}= 5.36 TeV, 2.5<#kern[-1.5]{ }#it{y}<4, 0-100%","P");
 legpt4->AddEntry(gRAA_ALICE_pub_fw090,"Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}}= 5.02 TeV, 2.5<#kern[-1.5]{ }#it{y}<4, 70-90%, Journal","P");
 legpt4->Draw();
//TLatex *latexTitle2 = new TLatex();
// latexTitle2 -> SetTextSize(0.028);
// latexTitle2 -> SetNDC();
// latexTitle2 -> SetTextFont(42);
// latexTitle2 -> DrawLatex(0.72, 0.77, "Phys. Lett. B 846 (2023) 137467");

//TLatex *latexTitle2 = new TLatex();
//latexTitle2 -> SetTextSize(0.037);
//latexTitle2 -> SetNDC();
//latexTitle2 -> SetTextFont(42);
//latexTitle2 -> DrawLatex(0.14, 0.64, "Phys. Lett. B 846 (2023) 137467");

 

 TLine *lpt4 = new TLine(0.,1.,13.,1.);
 lpt4->SetLineColor(1);
 lpt4->SetLineWidth(2);
 lpt4->SetLineStyle(2);
 lpt4->Draw();
   
   
 gPad->RedrawAxis();


  
}


void LoadStyle(){
  int font = 42;
  gROOT->SetStyle("Plain");
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(0.9,"y");
  gStyle->SetTitleSize(0.05,"x");  
  gStyle->SetTitleSize(0.06,"yz");  
  gStyle->SetMarkerSize(1.3); 
  gStyle->SetPalette(1,0); 
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(10);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetEndErrorSize(0);
}


