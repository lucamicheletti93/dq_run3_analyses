
// SAME as plot_results_OO_pO_PbPb_pPb.C -> only PbPb and OO shown
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

void plot_results_OO_PbPb(){

 LoadStyle();


 Double_t Glob_RAA_ALICE_pub_fw2 = 0.053;
 TBox *boxglobRAA_ALICE_pub_fw = new TBox(13.0,1.-Glob_RAA_ALICE_pub_fw2,13.5,1+Glob_RAA_ALICE_pub_fw2);     
 boxglobRAA_ALICE_pub_fw->SetFillColor(kGreen+3);
  
  //----------------------------------------------------------
  // ALICE RAA fw-y 2.5<y<4, 2015 0-90%. sqrt(s_NN) = 5.02 TeV (2015)
  // published fig 2
  // https://arxiv.org/pdf/1909.03158 JHEP 2002 (2020) 041
  //----------------------------------------------------------

  const int n_ALICE_pub_fw090=11;
  Double_t pt_ALICE_pub_fw090[n_ALICE_pub_fw090]={0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,11};
  Double_t ept_ALICE_pub_fw090[n_ALICE_pub_fw090]={0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1.};
  Double_t esyst_pt_ALICE_pub_fw090[n_ALICE_pub_fw090]={0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2};
  Double_t RAA_ALICE_pub_fw090[n_ALICE_pub_fw090]={0.75,0.72,0.65,0.54,0.53,0.42,0.39,0.38,0.42,0.29,0.41};
  Double_t Stat_RAA_ALICE_pub_fw090[n_ALICE_pub_fw090]={0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02,0.02,0.02};
  Double_t Syst_RAA_ALICE_pub_fw090[n_ALICE_pub_fw090]={0.06,0.05,0.05,0.04,0.05,0.03,0.04,0.05,0.08,0.06,0.09};
  Double_t Glob_RAA_ALICE_pub_fw = 0.039;

  TGraphErrors *gRAA_ALICE_pub_fw090 = new TGraphErrors(n_ALICE_pub_fw090,pt_ALICE_pub_fw090,RAA_ALICE_pub_fw090,ept_ALICE_pub_fw090,Stat_RAA_ALICE_pub_fw090);
  gRAA_ALICE_pub_fw090->SetMarkerColor(kGreen+3);
  gRAA_ALICE_pub_fw090->SetLineColor(kBlack);
  gRAA_ALICE_pub_fw090->SetMarkerStyle(22);
  gRAA_ALICE_pub_fw090->SetMarkerSize(1.5);
  
  TGraphErrors *gRAA_ALICE_pub_fw090_syst = new TGraphErrors(n_ALICE_pub_fw090,pt_ALICE_pub_fw090,RAA_ALICE_pub_fw090,esyst_pt_ALICE_pub_fw090,Syst_RAA_ALICE_pub_fw090);
  gRAA_ALICE_pub_fw090_syst->SetMarkerColor(kGreen+3);
  gRAA_ALICE_pub_fw090_syst->SetMarkerStyle(22);
  gRAA_ALICE_pub_fw090_syst->SetFillColorAlpha(kGreen+3,0.);
  gRAA_ALICE_pub_fw090_syst->SetLineColor(kGreen+3);
  gRAA_ALICE_pub_fw090_syst->SetMarkerSize(1.5);

  TGraphErrors *gRAA_ALICE_pub_fw090_syst2 = new TGraphErrors(n_ALICE_pub_fw090,pt_ALICE_pub_fw090,RAA_ALICE_pub_fw090,esyst_pt_ALICE_pub_fw090,Syst_RAA_ALICE_pub_fw090);
  gRAA_ALICE_pub_fw090_syst2->SetMarkerColor(kGreen+3);
  gRAA_ALICE_pub_fw090_syst2->SetMarkerStyle(22);
  gRAA_ALICE_pub_fw090_syst2->SetFillStyle(0);
  gRAA_ALICE_pub_fw090_syst2->SetLineColor(kGreen+3);
  gRAA_ALICE_pub_fw090_syst2->SetMarkerSize(1.5);

  TBox *boxglobRAA_ALICE_pub_fw090 = new TBox(14.5,1.-Glob_RAA_ALICE_pub_fw,15,1+Glob_RAA_ALICE_pub_fw);     
  boxglobRAA_ALICE_pub_fw090->SetFillColor(kGreen+3);


  //----------------------------------------------------------
  // ALICE RpA fw-y sqrt(s_NN) = 8.16 TeV 2.03<y<3.53
  // published fig 5
  // https://arxiv.org/pdf/1805.04381 JHEP 1807 (2018) 160
  //----------------------------------------------------------

  const int n_ALICE_pA_fw=9;
  Double_t pt_ALICE_pA_fw[n_ALICE_pA_fw]={0.5,1.5,2.5,3.5,4.5,5.5,7.,9.,11};
  Double_t ept_ALICE_pA_fw[n_ALICE_pA_fw]={0.5,0.5,0.5,0.5,0.5,0.5,1,1,1};
  Double_t esyst_pt_ALICE_pA_fw[n_ALICE_pA_fw]={0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2};
  Double_t RpA_ALICE_pA_fw[n_ALICE_pA_fw]={0.57,0.6,0.69,0.77,0.82,0.87,0.90,0.94,1.01};
  Double_t Stat_RpA_ALICE_pA_fw[n_ALICE_pA_fw]={0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.03};
  Double_t Syst_RpA_ALICE_pA_fw[n_ALICE_pA_fw]={0.04,0.05,0.05,0.05,0.06,0.06,0.06,0.06,0.07};
  Double_t GlobError_RpA_ALICE_pA_fw=0.073;

  TGraphErrors *gRpA_ALICE_pA_fw = new TGraphErrors(n_ALICE_pA_fw,pt_ALICE_pA_fw,RpA_ALICE_pA_fw,ept_ALICE_pA_fw,Stat_RpA_ALICE_pA_fw);
  gRpA_ALICE_pA_fw->SetMarkerColor(kGray+2);
  gRpA_ALICE_pA_fw->SetLineColor(kGray+2);
  gRpA_ALICE_pA_fw->SetMarkerStyle(20);
  gRpA_ALICE_pA_fw->SetMarkerSize(1.5);
  
  TGraphErrors *gRpA_ALICE_pA_fw_syst = new TGraphErrors(n_ALICE_pA_fw,pt_ALICE_pA_fw,RpA_ALICE_pA_fw,esyst_pt_ALICE_pA_fw,Syst_RpA_ALICE_pA_fw);
  gRpA_ALICE_pA_fw_syst->SetMarkerColor(kGray+2);
  gRpA_ALICE_pA_fw_syst->SetMarkerStyle(20);
  gRpA_ALICE_pA_fw_syst->SetFillColorAlpha(kGray+1,0.3);
  gRpA_ALICE_pA_fw_syst->SetLineColor(kGray+2);
  gRpA_ALICE_pA_fw_syst->SetMarkerSize(1.5);

  TGraphErrors *gRpA_ALICE_pA_fw_syst2 = new TGraphErrors(n_ALICE_pA_fw,pt_ALICE_pA_fw,RpA_ALICE_pA_fw,esyst_pt_ALICE_pA_fw,Syst_RpA_ALICE_pA_fw);
  gRpA_ALICE_pA_fw_syst2->SetMarkerColor(kGray+2);
  gRpA_ALICE_pA_fw_syst2->SetMarkerStyle(20);
  gRpA_ALICE_pA_fw_syst2->SetFillStyle(0);
  gRpA_ALICE_pA_fw_syst2->SetLineColor(kGray+2);
  gRpA_ALICE_pA_fw_syst2->SetMarkerSize(1.5);

  TBox *boxglobRpA_ALICE_pA_fw = new TBox(12.0,1.-GlobError_RpA_ALICE_pA_fw,12.5,1+GlobError_RpA_ALICE_pA_fw);     
  boxglobRpA_ALICE_pA_fw->SetFillColor(kGray+2);


  //----------------------------------------------------------
  // ALICE RpA fw-y sqrt(s_NN) = 8.16 TeV -4.46<y<-2.96
  // published fig 5 left
  // https://arxiv.org/pdf/1805.04381 JHEP 1807 (2018) 160
  //----------------------------------------------------------

  const int n_ALICE_Ap_bck=9;
  Double_t pt_ALICE_Ap_bck[n_ALICE_Ap_bck]={0.5,1.5,2.5,3.5,4.5,5.5,7.,9.,11};
  Double_t ept_ALICE_Ap_bck[n_ALICE_Ap_bck]={0.5,0.5,0.5,0.5,0.5,0.5,1,1,1};
  Double_t esyst_pt_ALICE_Ap_bck[n_ALICE_Ap_bck]={0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2};
  Double_t RpA_ALICE_Ap_bck[n_ALICE_Ap_bck]={0.82,0.93,1.06,1.13,1.17,1.19,1.14,1.12,1.16};
  Double_t Stat_RpA_ALICE_Ap_bck[n_ALICE_Ap_bck]={0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.04};
  Double_t Syst_RpA_ALICE_Ap_bck[n_ALICE_Ap_bck]={0.07,0.07,0.08,0.08,0.08,0.08,0.07,0.07,0.08};
  Double_t GlobError_RpA_ALICE_Ap_bck=0.073;

  TGraphErrors *gRpA_ALICE_Ap_bck = new TGraphErrors(n_ALICE_Ap_bck,pt_ALICE_Ap_bck,RpA_ALICE_Ap_bck,ept_ALICE_Ap_bck,Stat_RpA_ALICE_Ap_bck);
  gRpA_ALICE_Ap_bck->SetMarkerColor(kGray+3);
  gRpA_ALICE_Ap_bck->SetLineColor(kGray+3);
  gRpA_ALICE_Ap_bck->SetMarkerStyle(20);
  gRpA_ALICE_Ap_bck->SetMarkerSize(1.5);
  
  TGraphErrors *gRpA_ALICE_Ap_bck_syst = new TGraphErrors(n_ALICE_Ap_bck,pt_ALICE_Ap_bck,RpA_ALICE_Ap_bck,esyst_pt_ALICE_Ap_bck,Syst_RpA_ALICE_Ap_bck);
  gRpA_ALICE_Ap_bck_syst->SetMarkerColor(kGray+3);
  gRpA_ALICE_Ap_bck_syst->SetMarkerStyle(20);
  gRpA_ALICE_Ap_bck_syst->SetFillColorAlpha(kGray+3,0.3);
  gRpA_ALICE_Ap_bck_syst->SetLineColor(kGray+3);
  gRpA_ALICE_Ap_bck_syst->SetMarkerSize(1.5);

  TGraphErrors *gRpA_ALICE_Ap_bck_syst2 = new TGraphErrors(n_ALICE_Ap_bck,pt_ALICE_Ap_bck,RpA_ALICE_Ap_bck,esyst_pt_ALICE_Ap_bck,Syst_RpA_ALICE_Ap_bck);
  gRpA_ALICE_Ap_bck_syst2->SetMarkerColor(kGray+3);
  gRpA_ALICE_Ap_bck_syst2->SetMarkerStyle(20);
  gRpA_ALICE_Ap_bck_syst2->SetFillStyle(0);
  gRpA_ALICE_Ap_bck_syst2->SetLineColor(kGray+3);
  gRpA_ALICE_Ap_bck_syst2->SetMarkerSize(1.5);

  TBox *boxglobRpA_ALICE_pA_bck = new TBox(12.0,1.-GlobError_RpA_ALICE_Ap_bck,12.5,1+GlobError_RpA_ALICE_Ap_bck);     
  boxglobRpA_ALICE_pA_bck->SetFillColor(kGray+3);

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
  gROO_ALICE_syst->SetFillColorAlpha(kRed+1,0.);
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

  //----------------------------------------------------------
  // ALICE RpA fw-y 2.15<y<3.65, 2025 0-100%. sqrt(s_NN) = 9.62 TeV
  // preliminary
  //----------------------------------------------------------
  const int n_ALICE_pO=6;
  Double_t pt_ALICE_pO[n_ALICE_pO]={0.5, 1.5, 2.5, 3.5, 5., 7.};
  Double_t ept_ALICE_pO[n_ALICE_pO]={0.5, 0.5, 0.5, 0.5, 1, 1};
  Double_t esyst_pt_ALICE_pO[n_ALICE_pO]={0.2, 0.2, 0.2, 0.2, 0.2, 0.2};
  Double_t RpO_ALICE_pO[n_ALICE_pO]={0.719, 0.779, 0.797, 0.857, 0.864, 0.996};
  Double_t Stat_RpO_ALICE_pO[n_ALICE_pO]={0.0311327, 0.02479557, 0.02758417, 0.03219749, 0.03030912, 0.06258864};
  Double_t Syst_RpO_ALICE_pO[n_ALICE_pO]={0.04935935, 0.05198267, 0.04952558, 0.05080296, 0.04625856, 0.060258};
  Double_t GlobError_RpO_ALICE=0.067;
  TGraphErrors *gRpO_ALICE_pO = new TGraphErrors(n_ALICE_pO,pt_ALICE_pO,RpO_ALICE_pO,ept_ALICE_pO,Stat_RpO_ALICE_pO);
  gRpO_ALICE_pO->SetMarkerColor(kAzure+2);
  gRpO_ALICE_pO->SetLineColor(kBlack);
  gRpO_ALICE_pO->SetMarkerStyle(20);
  gRpO_ALICE_pO->SetMarkerSize(1.5);
  
  TGraphErrors *gRpO_ALICE_syst = new TGraphErrors(n_ALICE_pO,pt_ALICE_pO,RpO_ALICE_pO,esyst_pt_ALICE_pO,Syst_RpO_ALICE_pO);
  gRpO_ALICE_syst->SetMarkerColor(kAzure+2);
  gRpO_ALICE_syst->SetMarkerStyle(20);
  gRpO_ALICE_syst->SetFillColorAlpha(kAzure+2,0.);
  gRpO_ALICE_syst->SetLineColor(kAzure+2);
  gRpO_ALICE_syst->SetLineWidth(2);
  gRpO_ALICE_syst->SetMarkerSize(1.5);
  TGraphErrors *gRpO_ALICE_syst2 = new TGraphErrors(n_ALICE_pO,pt_ALICE_pO,RpO_ALICE_pO,esyst_pt_ALICE_pO,Syst_RpO_ALICE_pO);
  gRpO_ALICE_syst2->SetMarkerColor(kAzure+2);
  gRpO_ALICE_syst2->SetMarkerStyle(20);
  gRpO_ALICE_syst2->SetFillStyle(0);
  gRpO_ALICE_syst2->SetLineColor(kAzure+2);
  gRpO_ALICE_syst2->SetLineWidth(2);
  gRpO_ALICE_syst2->SetMarkerSize(1.5);
  TBox *boxglobRpO_ALICE = new TBox(12.0,1.-GlobError_RpO_ALICE,12.5,1+GlobError_RpO_ALICE);     
  boxglobRpO_ALICE->SetFillColor(kAzure+2);


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

 TH2D *hpt3 = new TH2D("hpt3","",100,0.,13.5,100,0.,1.5);
 //hpt3->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#kern[-1.4]{ }#it{c})");
  hpt3->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
 hpt3->GetYaxis()->SetTitle("#it{R}_{AA}");
 hpt3->Draw();
  
 //gRpA_ALICE_pA_fw_syst->Draw("E2P");
 //gRpA_ALICE_pA_fw_syst2->Draw("E2P");
 //gRpA_ALICE_pA_fw->Draw("P");

 //gRpA_ALICE_Ap_bck_syst->Draw("E2P");
 //gRpA_ALICE_Ap_bck_syst2->Draw("E2P");
 //gRpA_ALICE_Ap_bck->Draw("P");

 gRAA_ALICE_pub_fw090_syst->Draw("E2P");
 gRAA_ALICE_pub_fw090_syst2->Draw("E2P");
 gRAA_ALICE_pub_fw090->Draw("P");

 gROO_ALICE_syst->Draw("E2P");
 gROO_ALICE_syst2->Draw("E2P");
 gROO_ALICE_OO->Draw("P");

  //gRpO_ALICE_syst->Draw("E2P");
  //gRpO_ALICE_syst2->Draw("E2P");
  //gRpO_ALICE_pO->Draw("P");
 
 boxglobRAA_ALICE_pub_fw->Draw();
 //boxglobRpA_ALICE_pA_fw->Draw(); 
 //boxglobRpA_ALICE_pA_bck->Draw();
 //boxglobRpO_ALICE->Draw(); 
 boxglobROO_ALICE->Draw();

 TLatex *latexTitle = new TLatex();
 latexTitle -> SetTextSize(0.050);
 latexTitle -> SetNDC();
 latexTitle -> SetTextFont(42);
 //latexTitle -> DrawLatex(0.15, 0.90, "ALICE Preliminary, J/#kern[-2.1]{ }#psi #kern[-1.9]{ }#rightarrow #kern[-1.9]{ }#mu^{+}#mu^{-}");
  latexTitle -> DrawLatex(0.15, 0.90, "ALICE Preliminary, J/#psi#rightarrow#mu^{+}#mu^{-}");

 TLegend *legpt4 = new TLegend(0.11,0.73,0.45,0.88);
 legpt4->SetBorderSize(0);
 legpt4->SetFillColor(10);
 legpt4->SetFillStyle(1);
 legpt4->SetLineStyle(0);
 legpt4->SetLineColor(0);
 legpt4->SetTextSize(0.039);
 legpt4->AddEntry(gROO_ALICE_OO,"OO #sqrt{#it{s}_{NN}}= 5.36 TeV, 2.5<#it{y}<4, 0#font[122]{-}100%","P");
 //legpt4->AddEntry(gRpO_ALICE_pO,"pO #sqrt{#it{s}_{NN}} = 9.62 TeV, 2.5<#it{y}<4, 0#font[122]{-}100%","P");
 //legpt4->AddEntry(gROO_ALICE_OO,"OO #sqrt{#it{s}_{NN}}= 5.36 TeV, 2.5<#kern[-1.7]{ }#it{y}#kern[-1.1]{ }<4, 0#font[122]{-}100%","P");
 //legpt4->AddEntry(gRpO_ALICE_pO,"pO #sqrt{#it{s}_{NN}} = 9.62 TeV, 2.5<#kern[-1.7]{ }#it{y}#kern[-1.1]{ }<4, 0#font[122]{-}100%","P");
 legpt4->AddEntry(gRAA_ALICE_pub_fw090,"Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}}= 5.02 TeV, 2.5<#it{y}<4, 0#font[122]{-}90%, JHEP 2002 (2020) 041","P");
 //legpt4->AddEntry(gRpA_ALICE_pA_fw,"p-Pb #sqrt{#it{s}_{NN}}= 8.16 TeV, 2.03<#it{y}<3.53, JHEP 1807 (20s18) 160","P");
 //legpt4->AddEntry(gRpA_ALICE_Ap_bck,"p-Pb #sqrt{#it{s}_{NN}}= 8.16 TeV, -4.46<#it{y}<-2.96, JHEP 1807 (2018) 160","P");
 legpt4->Draw();

 TLine *lpt4 = new TLine(0.,1.,13.5,1.);
 lpt4->SetLineColor(1);
 lpt4->SetLineWidth(2);
 lpt4->SetLineStyle(2);
 lpt4->Draw();
   
   
 gPad->RedrawAxis();

 c3 -> SaveAs("SQM2026/comparisonPbPb_Vs_OO.pdf");
  c3 -> SaveAs("SQM2026/comparisonPbPb_Vs_OO.png");


  
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


