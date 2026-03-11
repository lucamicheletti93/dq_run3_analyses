void LoadStyle();
void SetLegend(TLegend *);

void plot_results() {
    LoadStyle();

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

    TGraphErrors *gRAA_ALICE_pub_fw090 = new TGraphErrors(n_ALICE_pub_fw090-3,pt_ALICE_pub_fw090,RAA_ALICE_pub_fw090,ept_ALICE_pub_fw090,Stat_RAA_ALICE_pub_fw090);
    gRAA_ALICE_pub_fw090->SetMarkerColor(kRed+1);
    gRAA_ALICE_pub_fw090->SetLineColor(kRed+1);
    gRAA_ALICE_pub_fw090->SetMarkerStyle(33);
    gRAA_ALICE_pub_fw090->SetMarkerSize(2);
    
    TGraphErrors *gRAA_ALICE_pub_fw090_syst = new TGraphErrors(n_ALICE_pub_fw090-3,pt_ALICE_pub_fw090,RAA_ALICE_pub_fw090,esyst_pt_ALICE_pub_fw090,Syst_RAA_ALICE_pub_fw090);
    gRAA_ALICE_pub_fw090_syst->SetMarkerColor(kRed+1);
    gRAA_ALICE_pub_fw090_syst->SetMarkerStyle(33);
    gRAA_ALICE_pub_fw090_syst->SetFillColorAlpha(kRed+1,0.3);
    gRAA_ALICE_pub_fw090_syst->SetLineColor(kRed+1);
    gRAA_ALICE_pub_fw090_syst->SetMarkerSize(2);

    TGraphErrors *gRAA_ALICE_pub_fw090_syst2 = new TGraphErrors(n_ALICE_pub_fw090-3,pt_ALICE_pub_fw090,RAA_ALICE_pub_fw090,esyst_pt_ALICE_pub_fw090,Syst_RAA_ALICE_pub_fw090);
    gRAA_ALICE_pub_fw090_syst2->SetMarkerColor(kRed+1);
    gRAA_ALICE_pub_fw090_syst2->SetMarkerStyle(33);
    gRAA_ALICE_pub_fw090_syst2->SetFillStyle(0);
    gRAA_ALICE_pub_fw090_syst2->SetLineColor(kRed+1);
    gRAA_ALICE_pub_fw090_syst2->SetMarkerSize(2);

    TBox *boxglobRAA_ALICE_pub_fw090 = new TBox(14.5,1.-Glob_RAA_ALICE_pub_fw,15,1+Glob_RAA_ALICE_pub_fw);     
    boxglobRAA_ALICE_pub_fw090->SetFillColor(kRed+1);

    Double_t Glob_RAA_ALICE_pub_fw2 = 0.053;
    TBox *boxglobRAA_ALICE_pub_fw = new TBox(8.75,1.-Glob_RAA_ALICE_pub_fw2,8.99,1+Glob_RAA_ALICE_pub_fw2);     
    boxglobRAA_ALICE_pub_fw->SetFillColor(kRed+1);

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

    TGraphErrors *gRpA_ALICE_pA_fw = new TGraphErrors(n_ALICE_pA_fw-2,pt_ALICE_pA_fw,RpA_ALICE_pA_fw,ept_ALICE_pA_fw,Stat_RpA_ALICE_pA_fw);
    gRpA_ALICE_pA_fw->SetMarkerColor(kAzure+2);
    gRpA_ALICE_pA_fw->SetLineColor(kAzure+2);
    gRpA_ALICE_pA_fw->SetMarkerStyle(33);
    gRpA_ALICE_pA_fw->SetMarkerSize(2);
    
    TGraphErrors *gRpA_ALICE_pA_fw_syst = new TGraphErrors(n_ALICE_pA_fw-2,pt_ALICE_pA_fw,RpA_ALICE_pA_fw,esyst_pt_ALICE_pA_fw,Syst_RpA_ALICE_pA_fw);
    gRpA_ALICE_pA_fw_syst->SetMarkerColor(kAzure+2);
    gRpA_ALICE_pA_fw_syst->SetMarkerStyle(33);
    gRpA_ALICE_pA_fw_syst->SetFillColorAlpha(kAzure+2,0.3);
    gRpA_ALICE_pA_fw_syst->SetLineColor(kAzure+2);
    gRpA_ALICE_pA_fw_syst->SetMarkerSize(2);

    TGraphErrors *gRpA_ALICE_pA_fw_syst2 = new TGraphErrors(n_ALICE_pA_fw-2,pt_ALICE_pA_fw,RpA_ALICE_pA_fw,esyst_pt_ALICE_pA_fw,Syst_RpA_ALICE_pA_fw);
    gRpA_ALICE_pA_fw_syst2->SetMarkerColor(kAzure+2);
    gRpA_ALICE_pA_fw_syst2->SetMarkerStyle(33);
    gRpA_ALICE_pA_fw_syst2->SetFillStyle(0);
    gRpA_ALICE_pA_fw_syst2->SetLineColor(kAzure+2);
    gRpA_ALICE_pA_fw_syst2->SetMarkerSize(2);

    TBox *boxglobRpA_ALICE_pA_fw = new TBox(8.5,1.-GlobError_RpA_ALICE_pA_fw,8.75,1+GlobError_RpA_ALICE_pA_fw);     
    boxglobRpA_ALICE_pA_fw->SetFillColor(kAzure+2);

    //----------------------------------------------------------
    // ALICE RAA fw-y 2.5<y<4, 2025 0-100%. sqrt(s_NN) = 5.36 TeV
    // preliminary
    //----------------------------------------------------------
    const int n_ALICE_OO=7;
    Double_t pt_ALICE_OO[n_ALICE_OO]={0.5, 1.5,2.5,3.5,4.5,5.5,7.};
    Double_t ept_ALICE_OO[n_ALICE_OO]={0.5,0.5,0.5,0.5,0.5,0.5,1};
    Double_t esyst_pt_ALICE_OO[n_ALICE_OO]={0.2,0.2,0.2,0.2,0.2,0.2,0.2};
    Double_t ROO_ALICE_OO[n_ALICE_OO]={0.581728,0.629188,0.641791,0.656929,0.702026,0.698702,0.706779};
    Double_t Stat_ROO_ALICE_OO[n_ALICE_OO]={0.014429,0.011421,0.012623,0.015387,0.021201,0.025506,0.024434};
    Double_t Syst_ROO_ALICE_OO[n_ALICE_OO]={0.034578,0.033524,0.032875,0.035657,0.049118,0.049178,0.047188};
    Double_t GlobError_ROO_ALICE=0.06461423991660044;

    TGraphErrors *gROO_ALICE_OO = new TGraphErrors(n_ALICE_OO,pt_ALICE_OO,ROO_ALICE_OO,ept_ALICE_OO,Stat_ROO_ALICE_OO);
    gROO_ALICE_OO->SetMarkerColor(kRed+1);
    gROO_ALICE_OO->SetLineColor(kRed+1);
    gROO_ALICE_OO->SetMarkerStyle(20);
    gROO_ALICE_OO->SetMarkerSize(1.5);
    
    TGraphErrors *gROO_ALICE_syst = new TGraphErrors(n_ALICE_OO,pt_ALICE_OO,ROO_ALICE_OO,esyst_pt_ALICE_OO,Syst_ROO_ALICE_OO);
    gROO_ALICE_syst->SetMarkerColor(kRed+1);
    gROO_ALICE_syst->SetMarkerStyle(20);
    gROO_ALICE_syst->SetFillColorAlpha(kRed+1,0.3);
    gROO_ALICE_syst->SetLineColor(kRed+1);
    gROO_ALICE_syst->SetMarkerSize(1.5);

    TGraphErrors *gROO_ALICE_syst2 = new TGraphErrors(n_ALICE_OO,pt_ALICE_OO,ROO_ALICE_OO,esyst_pt_ALICE_OO,Syst_ROO_ALICE_OO);
    gROO_ALICE_syst2->SetMarkerColor(kRed+1);
    gROO_ALICE_syst2->SetMarkerStyle(20);
    gROO_ALICE_syst2->SetFillStyle(0);
    gROO_ALICE_syst2->SetLineColor(kRed+1);
    gROO_ALICE_syst2->SetMarkerSize(1.5);

    TBox *boxglobROO_ALICE = new TBox(8.75,1.-GlobError_ROO_ALICE,8.998,1+GlobError_ROO_ALICE);     
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
    Double_t GlobError_RpO_ALICE=0.12292274;


    TGraphErrors *gRpO_ALICE_pO = new TGraphErrors(n_ALICE_pO,pt_ALICE_pO,RpO_ALICE_pO,ept_ALICE_pO,Stat_RpO_ALICE_pO);
    gRpO_ALICE_pO->SetMarkerColor(kAzure+2);
    gRpO_ALICE_pO->SetLineColor(kAzure+2);
    gRpO_ALICE_pO->SetMarkerStyle(20);
    gRpO_ALICE_pO->SetMarkerSize(1.5);
    
    TGraphErrors *gRpO_ALICE_syst = new TGraphErrors(n_ALICE_pO,pt_ALICE_pO,RpO_ALICE_pO,esyst_pt_ALICE_pO,Syst_RpO_ALICE_pO);
    gRpO_ALICE_syst->SetMarkerColor(kAzure+2);
    gRpO_ALICE_syst->SetMarkerStyle(20);
    gRpO_ALICE_syst->SetFillColorAlpha(kAzure+2,0.3);
    gRpO_ALICE_syst->SetLineColor(kAzure+2);
    gRpO_ALICE_syst->SetMarkerSize(1.5);

    TGraphErrors *gRpO_ALICE_syst2 = new TGraphErrors(n_ALICE_pO,pt_ALICE_pO,RpO_ALICE_pO,esyst_pt_ALICE_pO,Syst_RpO_ALICE_pO);
    gRpO_ALICE_syst2->SetMarkerColor(kAzure+2);
    gRpO_ALICE_syst2->SetMarkerStyle(20);
    gRpO_ALICE_syst2->SetFillStyle(0);
    gRpO_ALICE_syst2->SetLineColor(kAzure+2);
    gRpO_ALICE_syst2->SetMarkerSize(1.5);

    TBox *boxglobRpO_ALICE = new TBox(8.50,1.-GlobError_RpO_ALICE,8.75,1+GlobError_RpO_ALICE);     
    boxglobRpO_ALICE->SetFillColor(kAzure+2);


    /////////////////////////////////
    TCanvas *canvasPbPbVspPbVsOOVspO = new TCanvas("canvasPbPbVspPbVsOOVspO", "", 1800, 600);

    TH2D *histGrid = new TH2D("histGrid", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}, #it{R}_{pA}", 100, 0.001, 8.999, 100, 0, 1.5);
    TLine *lineUnity = new TLine(0, 1, 9, 1);
    lineUnity -> SetLineColor(kGray+1);
    lineUnity -> SetLineStyle(kDashed);

    canvasPbPbVspPbVsOOVspO -> cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0.05, 0.05, 0.5, 0.95); 
    pad1 -> SetBottomMargin(0.15);
    pad1 -> SetTopMargin(0);
    pad1 -> SetLeftMargin(0.15); 
    pad1 -> SetRightMargin(0); 
    pad1 -> Draw(); 
    pad1 -> cd();
    gPad -> SetTicks();
    histGrid -> Draw();
    lineUnity -> Draw();

    gROO_ALICE_syst->Draw("E2P");
    gROO_ALICE_syst2->Draw("E2P");
    gROO_ALICE_OO->Draw("P");

    gRpO_ALICE_syst->Draw("E2P");
    gRpO_ALICE_syst2->Draw("E2P");
    gRpO_ALICE_pO->Draw("P");

    boxglobROO_ALICE->Draw();
    boxglobRpO_ALICE->Draw();

    TLegend *legend1 = new TLegend(0.16,0.75,0.45,0.90);
    legend1->SetBorderSize(0);
    legend1->SetFillColor(10);
    legend1->SetFillStyle(1);
    legend1->SetLineStyle(0);
    legend1->SetLineColor(0);
    legend1->SetTextSize(0.042);
    legend1->AddEntry(gROO_ALICE_OO,"OO #sqrt{#it{s}_{NN}} = 5.36 TeV, 2.5 < #it{y} < 4, 0#font[122]{-}100%","P");
    legend1->AddEntry(gRpO_ALICE_pO,"pO #sqrt{#it{s}_{NN}} = 9.62 TeV, 2.15 < #it{y} < 3.65, 0#font[122]{-}100%","P");
    legend1->Draw();

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.050);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);
    latexTitle -> DrawLatex(0.20, 0.92, "ALICE Preliminary, J/#psi #rightarrow #mu^{+}#mu^{-}");

    canvasPbPbVspPbVsOOVspO -> cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.5, 0.05, 0.95, 0.95); 
    pad2 -> SetBottomMargin(0.15);
    pad2 -> SetTopMargin(0);
    pad2 -> SetLeftMargin(0); 
    pad2 -> SetRightMargin(0.15); 
    pad2 -> Draw(); 
    pad2 -> cd();
    gPad -> SetTicks();
    histGrid -> Draw();
    lineUnity -> Draw();

    gRAA_ALICE_pub_fw090_syst->Draw("E2P");
    gRAA_ALICE_pub_fw090_syst2->Draw("E2P");
    gRAA_ALICE_pub_fw090->Draw("P");

    gRpA_ALICE_pA_fw_syst->Draw("E2P");
    gRpA_ALICE_pA_fw_syst2->Draw("E2P");
    gRpA_ALICE_pA_fw->Draw("P");

    boxglobRAA_ALICE_pub_fw->Draw();
    boxglobRpA_ALICE_pA_fw->Draw(); 

    TLegend *legend2 = new TLegend(0.01,0.75,0.35,0.90);
    legend2->SetBorderSize(0);
    legend2->SetFillColor(10);
    legend2->SetFillStyle(1);
    legend2->SetLineStyle(0);
    legend2->SetLineColor(0);
    legend2->SetTextSize(0.042);
    legend2->AddEntry(gRAA_ALICE_pub_fw090,"Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, 2.5 < #it{y} < 4, 0#font[122]{-}90%, JHEP 2002 (2020) 041","P");
    legend2->AddEntry(gRpA_ALICE_pA_fw,"p#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 8.16 TeV, 2.03 < #it{y} < 3.53, JHEP 1807 (2018) 160","P");
    legend2->Draw();

    canvasPbPbVspPbVsOOVspO -> SaveAs("SQM2026/comparisonPbPbVspPb_Vs_OOVspO.pdf");
}
////////////////////////////////////////////////////////////////////////////////
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
    gStyle -> SetMarkerSize(1);
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
////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.045);
}