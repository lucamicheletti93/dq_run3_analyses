void LoadStyle();
void SetLegend(TLegend *);

inline void EvalXaxis(int nVarBins, double minVarX[], double maxVarX[], double varCentr[], double varWidthStat[], double varWidthSyst[]) {
    for (int iVar = 0;iVar < nVarBins;iVar++) {
        varCentr[iVar] = (maxVarX[iVar] + minVarX[iVar]) / 2.;
        varWidthStat[iVar] = (maxVarX[iVar] - minVarX[iVar]) / 2.;
        varWidthSyst[iVar] = 0.15;
    }
}

void plot_results_all_centr_and_pt() {
    LoadStyle();

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.050);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);

    TLine *lineZero = new TLine(-0.5, 0, 15.5, 0);
    lineZero -> SetLineStyle(kDashed);
    lineZero -> SetLineWidth(2);
    lineZero -> SetLineColor(kGray+1);

    // 5-50%
    const int nPtBins560Run3 = 3;
    double ptMin560Run3[] = {0.0, 3.0, 6.0};
    double ptMax560Run3[] = {3.0, 6.0, 15.0};
    double ptCentr560Run3[nPtBins560Run3], ptWidthStat560Run3[nPtBins560Run3], ptWidthSyst560Run3[nPtBins560Run3];
    EvalXaxis(nPtBins560Run3, ptMin560Run3, ptMax560Run3, ptCentr560Run3, ptWidthStat560Run3, ptWidthSyst560Run3);

    double v2Ups1sVsPt560Run3[] = {-0.009, 0.079, 0.004}; 
    double statV2Ups1sVsPt560Run3[] = {0.047, 0.052, 0.080};
    double systV2Ups1sVsPt560Run3[] = {0.009, 0.012, 0.022};

    TGraphErrors *graStatV2Ups1sVsPt560Run3 = new TGraphErrors(nPtBins560Run3, ptCentr560Run3, v2Ups1sVsPt560Run3, ptWidthStat560Run3, statV2Ups1sVsPt560Run3);
    graStatV2Ups1sVsPt560Run3 -> SetLineColor(kRed+1);
    graStatV2Ups1sVsPt560Run3 -> SetLineWidth(2);
    graStatV2Ups1sVsPt560Run3 -> SetMarkerStyle(20);
    graStatV2Ups1sVsPt560Run3 -> SetMarkerColor(kRed+1);

    TGraphErrors *graSystV2Ups1sVsPt560Run3 = new TGraphErrors(nPtBins560Run3, ptCentr560Run3, v2Ups1sVsPt560Run3, ptWidthSyst560Run3, systV2Ups1sVsPt560Run3);
    graSystV2Ups1sVsPt560Run3 -> SetLineColor(kRed+1);
    graSystV2Ups1sVsPt560Run3 -> SetLineWidth(2);
    graSystV2Ups1sVsPt560Run3 -> SetMarkerStyle(20);
    graSystV2Ups1sVsPt560Run3 -> SetMarkerColor(kRed+1);
    graSystV2Ups1sVsPt560Run3 -> SetFillStyle(0);

    // J/psi v2 vs pT
    // 5-50%
    TCanvas *canvasV2Ups1sVsPtCentr560 = new TCanvas("canvasV2Ups1sVsPtCentr560", "", 800, 600);
    canvasV2Ups1sVsPtCentr560 -> SetTicks(1, 1);

    TH2D *histGridV2Ups1sVsPtCentr560 = new TH2D("histGridV2Ups1sVsPtCentr560", "", 100, -0.5, 15.5, 100, -0.10, 0.3);
    histGridV2Ups1sVsPtCentr560 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    //histGridV2Ups1sVsPtCentr560 -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2Ups1sVsPtCentr560 -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2Ups1sVsPtCentr560 -> Draw();
    lineZero -> Draw();
    graStatV2Ups1sVsPt560Run3 -> Draw("EP SAME");
    graSystV2Ups1sVsPt560Run3 -> Draw("E2 SAME");

    TLegend *legendV2Ups1sVsPtCentr560 = new TLegend(0.65, 0.60, 0.85, 0.80, " ", "brNDC");
    SetLegend(legendV2Ups1sVsPtCentr560);
    legendV2Ups1sVsPtCentr560 -> SetTextSize(0.05);
    legendV2Ups1sVsPtCentr560 -> AddEntry(graStatV2Ups1sVsPt560Run3, "Stat. Uncert.", "PL");
    legendV2Ups1sVsPtCentr560 -> AddEntry(graSystV2Ups1sVsPt560Run3, "Syst. Uncert.", "F");
    legendV2Ups1sVsPtCentr560 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary, Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV");
    latexTitle -> DrawLatex(0.18, 0.77, "#varUpsilon(1S)#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4, 5#minus60\%");

    /////////////////////////////
    // J/psi v2 vs pT vs Run 2 //
    /////////////////////////////
    const int nPtBinsRun2 = 3;
    double ptCentrRun2[] = {1.8, 4.8, 10.8};
    double ptWidthStatRun2[] = {0, 0, 0};
    double ptWidthSystRun2[] = {0.15, 0.15, 0.15};

    double v2Ups1sVsPtCentr560Run2[] = {0.013, -0.010, 0.003};
    double statV2Ups1sVsPtCentr560Run2[] = {0.041, 0.041, 0.059};
    double systV2Ups1sVsPtCentr560Run2[] = {0.016, 0.008, 0.006};

    TGraphErrors *graStatV2Ups1sVsPtCentr560Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2Ups1sVsPtCentr560Run2, ptWidthStatRun2, statV2Ups1sVsPtCentr560Run2);
    graStatV2Ups1sVsPtCentr560Run2 -> SetLineColor(kGray+2);
    graStatV2Ups1sVsPtCentr560Run2 -> SetLineWidth(2);
    graStatV2Ups1sVsPtCentr560Run2 -> SetMarkerStyle(20);
    graStatV2Ups1sVsPtCentr560Run2 -> SetMarkerColor(kGray+2);

    TGraphErrors *graSystV2Ups1sVsPtCentr560Run2 = new TGraphErrors(nPtBinsRun2, ptCentrRun2, v2Ups1sVsPtCentr560Run2, ptWidthSystRun2, systV2Ups1sVsPtCentr560Run2);
    graSystV2Ups1sVsPtCentr560Run2 -> SetLineColor(kGray+2);
    graSystV2Ups1sVsPtCentr560Run2 -> SetLineWidth(2);
    graSystV2Ups1sVsPtCentr560Run2 -> SetMarkerStyle(20);
    graSystV2Ups1sVsPtCentr560Run2 -> SetMarkerColor(kGray+2);
    graSystV2Ups1sVsPtCentr560Run2 -> SetFillStyle(0);


    TCanvas *canvasV2Ups1sVsPtCentr560Run2VsRun3 = new TCanvas("canvasV2Ups1sVsPtCentr560Run2VsRun3", "", 800, 600);
    canvasV2Ups1sVsPtCentr560Run2VsRun3 -> SetTicks(1, 1);

    TH2D *histGridV2Ups1sVsPtCentr560Run2VsRun3 = new TH2D("histGridV2Ups1sVsPtCentr560Run2VsRun3", "", 100, -0.5, 15.5, 100, -0.10, 0.3);
    histGridV2Ups1sVsPtCentr560Run2VsRun3 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    //histGridV2Ups1sVsPtCentr560Run2VsRun3 -> GetYaxis() -> SetTitle("v_{2} {EP, |#Delta#eta > 2.5|}");
    histGridV2Ups1sVsPtCentr560Run2VsRun3 -> GetYaxis() -> SetTitle("v_{2}");
    histGridV2Ups1sVsPtCentr560Run2VsRun3 -> Draw();
    lineZero -> Draw();
    graStatV2Ups1sVsPtCentr560Run2 -> Draw("EP SAME");
    graSystV2Ups1sVsPtCentr560Run2 -> Draw("E2 SAME");
    graStatV2Ups1sVsPt560Run3 -> Draw("EP SAME");
    graSystV2Ups1sVsPt560Run3 -> Draw("E2 SAME");

    TLegend *legendV2Ups1sVsPtCentr560Run2VsRun3 = new TLegend(0.20, 0.70, 0.85, 0.80, " ", "brNDC");
    SetLegend(legendV2Ups1sVsPtCentr560Run2VsRun3);
    legendV2Ups1sVsPtCentr560Run2VsRun3 -> SetNColumns(2);
    legendV2Ups1sVsPtCentr560Run2VsRun3 -> SetTextSize(0.05);
    legendV2Ups1sVsPtCentr560Run2VsRun3 -> AddEntry(graSystV2Ups1sVsPtCentr560Run2, "#sqrt{#it{s}_{NN}} = 5.02 TeV", "PF");
    legendV2Ups1sVsPtCentr560Run2VsRun3 -> AddEntry(graSystV2Ups1sVsPt560Run3, "#sqrt{#it{s}_{NN}} = 5.36 TeV", "PF"); 
    legendV2Ups1sVsPtCentr560Run2VsRun3 -> Draw();

    latexTitle -> DrawLatex(0.18, 0.87, "ALICE Preliminary");
    latexTitle -> DrawLatex(0.18, 0.80, "Pb#minusPb, #varUpsilon(1S)#rightarrow#mu^{+}#mu^{-}, 2.5 < y < 4, 5#minus60\%");

    canvasV2Ups1sVsPtCentr560 -> SaveAs("v2Ups1sVsPtCentr560.pdf");
    canvasV2Ups1sVsPtCentr560Run2VsRun3 -> SaveAs("v2Ups1sVsPtCentr560Run2VsRun3.pdf");
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
////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.04);
}