void LoadStyle();
void SetLegend(TLegend *);

void accxeff() {
    LoadStyle();
    //-----------------------------------------------------------//
    // Run 2 Axe
    double ptBinsRun2[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 26.0, 30.0};
    double axePtRun2[] = {0.14923, 0.13782, 0.12011, 0.12440, 0.15167, 0.19440, 0.24144, 0.28442, 0.32242, 0.37021, 0.41486, 0.45900, 0.46818, 0.50202, 0.49940, 0.52091, 0.51781, 0.44586};
    double errAxePtRun2[] = {0.00089, 0.00054, 0.00031, 0.00034, 0.00047, 0.00068, 0.00099, 0.00140, 0.00193, 0.00208, 0.00349, 0.00549, 0.00804, 0.01212, 0.01678, 0.02332, 0.02437, 0.03951};

    TH1D *histPtAxeRun2 = new TH1D("histPtAxeRun2", "", 18, ptBinsRun2);
    histPtAxeRun2 -> SetLineColor(kBlack);
    histPtAxeRun2 -> SetMarkerStyle(20);
    histPtAxeRun2 -> SetMarkerColor(kBlack);
    for (int iPt = 0;iPt < 18;iPt++) {
        histPtAxeRun2 -> SetBinContent(iPt+1, axePtRun2[iPt]);
        histPtAxeRun2 -> SetBinError(iPt+1, errAxePtRun2[iPt]);
    }

    double yBinsRun2[] = {2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00};
    double axeYRun2[] = {0.04395, 0.15815, 0.23390, 0.24372, 0.18706, 0.07195};
    double errAxeYRun2[] = {0.00023, 0.00042, 0.00051, 0.00054, 0.00052, 0.00036};

    TH1D *histYAxeRun2 = new TH1D("histYAxeRun2", "", 6, yBinsRun2);
    histYAxeRun2 -> SetLineColor(kBlack);
    histYAxeRun2 -> SetMarkerStyle(20);
    histYAxeRun2 -> SetMarkerColor(kBlack);
    for (int iY = 0;iY < 6;iY++) {
        histYAxeRun2 -> SetBinContent(iY+1, axeYRun2[iY]);
        histYAxeRun2 -> SetBinError(iY+1, errAxeYRun2[iY]);
    }
    //-----------------------------------------------------------//

    TFile *fIn = new TFile("/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/MC/central_production/AnalysisResults.root", "READ");

    TList *listGen1 = (TList*) fIn -> Get("analysis-same-event-pairing/output");
    TList *listGen2 = (TList*) listGen1 -> FindObject("MCTruthGen_Jpsi");
    TH1D *histPtGen = (TH1D*)listGen2 -> FindObject("Pt");
    TH1D *histYGen = (TH1D*)listGen2 -> FindObject("Rapidity");

    histPtGen -> SetLineColor(kBlack);
    histYGen -> SetLineColor(kBlack);
    histPtGen -> SetMarkerColor(kBlack);
    histYGen -> SetMarkerColor(kBlack);

    TList *listRec1 = (TList*) fIn -> Get("analysis-same-event-pairing/output");

    TList *listRecCut1 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_matchedMchMid_mumuFromJpsi");
    TList *listRecCut2 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt10SigmaPDCA_mumuFromJpsi");
    TList *listRecCut3 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA_mumuFromJpsi");
    TList *listRecCut4 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt510SigmaPDCA_mumuFromJpsi");

    TH1D *histPtRecCut1 = (TH1D*)listRecCut1 -> FindObject("Pt");
    TH1D *histYRecCut1 = (TH1D*)listRecCut1 -> FindObject("Rapidity");

    TH1D *histPtRecCut2 = (TH1D*)listRecCut2 -> FindObject("Pt");
    TH1D *histYRecCut2 = (TH1D*)listRecCut2 -> FindObject("Rapidity");

    TH1D *histPtRecCut3 = (TH1D*)listRecCut3 -> FindObject("Pt");
    TH1D *histYRecCut3 = (TH1D*)listRecCut3 -> FindObject("Rapidity");

    TH1D *histPtRecCut4 = (TH1D*)listRecCut4 -> FindObject("Pt");
    TH1D *histYRecCut4 = (TH1D*)listRecCut4 -> FindObject("Rapidity");

    histPtRecCut1 -> SetLineColor(kRed+1);
    histYRecCut1 -> SetLineColor(kRed+1);
    histPtRecCut1 -> SetMarkerColor(kRed+1);
    histYRecCut1 -> SetMarkerColor(kRed+1);

    TH1D *histPtGenRebin = (TH1D*) histPtGen -> Rebin(18, "histPtGenRebin", ptBinsRun2); 
    TH1D *histPtRecCut1Rebin = (TH1D*) histPtRecCut1 -> Rebin(18, "histPtRecCut1Rebin", ptBinsRun2); 
    TH1D *histPtRecCut2Rebin = (TH1D*) histPtRecCut2 -> Rebin(18, "histPtRecCut2Rebin", ptBinsRun2); 
    TH1D *histPtRecCut3Rebin = (TH1D*) histPtRecCut3 -> Rebin(18, "histPtRecCut3Rebin", ptBinsRun2); 
    TH1D *histPtRecCut4Rebin = (TH1D*) histPtRecCut4 -> Rebin(18, "histPtRecCut4Rebin", ptBinsRun2); 

    histPtRecCut1Rebin -> SetLineColor(kRed+1);
    histPtRecCut2Rebin -> SetLineColor(kOrange+7);
    histPtRecCut3Rebin -> SetLineColor(kAzure+2);
    histPtRecCut4Rebin -> SetLineColor(kBlue+1);

    histPtRecCut1Rebin -> SetLineWidth(2);
    histPtRecCut2Rebin -> SetLineWidth(2);
    histPtRecCut3Rebin -> SetLineWidth(2);
    histPtRecCut4Rebin -> SetLineWidth(2);

    histPtRecCut1Rebin -> SetMarkerColor(kRed+1);
    histPtRecCut2Rebin -> SetMarkerColor(kOrange+7);
    histPtRecCut3Rebin -> SetMarkerColor(kAzure+2);
    histPtRecCut4Rebin -> SetMarkerColor(kBlue+1);

    TH1D *histPtAxeCut1 = (TH1D*) histPtRecCut1Rebin -> Clone("histPtAxeCut1");
    histPtAxeCut1 -> Divide(histPtGenRebin);
    histPtAxeCut1 -> SetLineColor(kRed+1);
    histPtAxeCut1 -> SetMarkerColor(kRed+1);

    TH1D *histPtAxeCut2 = (TH1D*) histPtRecCut2Rebin -> Clone("histPtAxeCut2");
    histPtAxeCut2 -> Divide(histPtGenRebin);
    histPtAxeCut2 -> SetLineColor(kOrange+7);
    histPtAxeCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histPtAxeCut3 = (TH1D*) histPtRecCut3Rebin -> Clone("histPtAxeCut3");
    histPtAxeCut3 -> Divide(histPtGenRebin);
    histPtAxeCut3 -> SetLineColor(kAzure+2);
    histPtAxeCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histPtAxeCut4 = (TH1D*) histPtRecCut4Rebin -> Clone("histPtAxeCut4");
    histPtAxeCut4 -> Divide(histPtGenRebin);
    histPtAxeCut4 -> SetLineColor(kBlue+1);
    histPtAxeCut4 -> SetMarkerColor(kBlue+1);

    TH1D *histYGenRebin = (TH1D*) histYGen -> Rebin(6, "histYGenRebin", yBinsRun2); 
    TH1D *histYRecCut1Rebin = (TH1D*) histYRecCut1 -> Rebin(6, "histYRecCut1Rebin", yBinsRun2); 
    TH1D *histYRecCut2Rebin = (TH1D*) histYRecCut2 -> Rebin(6, "histYRecCut2Rebin", yBinsRun2); 
    TH1D *histYRecCut3Rebin = (TH1D*) histYRecCut3 -> Rebin(6, "histYRecCut3Rebin", yBinsRun2); 
    TH1D *histYRecCut4Rebin = (TH1D*) histYRecCut4 -> Rebin(6, "histYRecCut4Rebin", yBinsRun2); 

    histYRecCut1Rebin -> SetLineColor(kRed+1);
    histYRecCut2Rebin -> SetLineColor(kOrange+7);
    histYRecCut3Rebin -> SetLineColor(kAzure+2);
    histYRecCut4Rebin -> SetLineColor(kBlue+1);

    histYRecCut1Rebin -> SetLineWidth(2);
    histYRecCut2Rebin -> SetLineWidth(2);
    histYRecCut3Rebin -> SetLineWidth(2);
    histYRecCut4Rebin -> SetLineWidth(2);

    histYRecCut1Rebin -> SetMarkerColor(kRed+1);
    histYRecCut2Rebin -> SetMarkerColor(kOrange+7);
    histYRecCut3Rebin -> SetMarkerColor(kAzure+2);
    histYRecCut4Rebin -> SetMarkerColor(kBlue+1);

    TH1D *histYAxeCut1 = (TH1D*) histYRecCut1Rebin -> Clone("histYAxeCut1");
    histYAxeCut1 -> Divide(histYGenRebin);
    histYAxeCut1 -> SetLineColor(kRed+1);
    histYAxeCut1 -> SetMarkerColor(kRed+1);

    TH1D *histYAxeCut2 = (TH1D*) histYRecCut2Rebin -> Clone("histYAxeCut2");
    histYAxeCut2 -> Divide(histYGenRebin);
    histYAxeCut2 -> SetLineColor(kOrange+7);
    histYAxeCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histYAxeCut3 = (TH1D*) histYRecCut3Rebin -> Clone("histYAxeCut3");
    histYAxeCut3 -> Divide(histYGenRebin);
    histYAxeCut3 -> SetLineColor(kAzure+2);
    histYAxeCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histYAxeCut4 = (TH1D*) histYRecCut4Rebin -> Clone("histYAxeCut4");
    histYAxeCut4 -> Divide(histYGenRebin);
    histYAxeCut4 -> SetLineColor(kBlue+1);
    histYAxeCut4 -> SetMarkerColor(kBlue+1);


    TCanvas *canvasPtGenRec = new TCanvas("canvasPtGenRec", "", 800, 600);
    gPad -> SetLogy(1);
    histPtGenRebin -> SetTitle("");
    histPtGenRebin -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histPtGenRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histPtGenRebin -> GetYaxis() -> SetRangeUser(1, 1e6);
    histPtGenRebin -> Draw("EP");
    histPtRecCut1Rebin -> Draw("EP SAME");
    histPtRecCut2Rebin -> Draw("EP SAME");
    histPtRecCut3Rebin -> Draw("EP SAME");
    histPtRecCut4Rebin -> Draw("EP SAME");

    TCanvas *canvasYGenRec = new TCanvas("canvasYGenRec", "", 800, 600);
    gPad -> SetLogy(1);
    histYGenRebin -> SetTitle("");
    histYGenRebin -> GetXaxis() -> SetTitle("#it{y}");
    histYGenRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histYGenRebin -> GetYaxis() -> SetRangeUser(1e3, 1e6);
    histYGenRebin -> Draw("EP");
    histYRecCut1Rebin -> Draw("EP SAME");
    histYRecCut2Rebin -> Draw("EP SAME");
    histYRecCut3Rebin -> Draw("EP SAME");
    histYRecCut4Rebin -> Draw("EP SAME");

    TCanvas *canvasPtAxe = new TCanvas("canvasPtAxe", "", 800, 600);
    histPtAxeCut1 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histPtAxeCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histPtAxeCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histPtAxeCut1 -> Draw("EP");
    histPtAxeCut2 -> Draw("EP SAME");
    histPtAxeCut3 -> Draw("EP SAME");
    histPtAxeCut4 -> Draw("EP SAME");
    histPtAxeRun2 -> Draw("EP SAME");

    TLegend *legendPtAxe = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendPtAxe);
    legendPtAxe -> AddEntry(histPtAxeRun2, "Run 2", "PL");
    legendPtAxe -> AddEntry(histPtAxeCut1, "Run 3 - matchedMchMid", "PL");
    legendPtAxe -> AddEntry(histPtAxeCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendPtAxe -> AddEntry(histPtAxeCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendPtAxe -> AddEntry(histPtAxeCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendPtAxe -> Draw("SAME");

    TCanvas *canvasYAxe = new TCanvas("canvasYAxe", "", 800, 600);
    histYAxeCut1 -> GetXaxis() -> SetTitle("#it{y}");
    histYAxeCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histYAxeCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histYAxeCut1 -> Draw("EP");
    histYAxeCut2 -> Draw("EP SAME");
    histYAxeCut3 -> Draw("EP SAME");
    histYAxeCut4 -> Draw("EP SAME");
    histYAxeRun2 -> Draw("EP SAME");

    TLegend *legendYAxe = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendYAxe);
    legendYAxe -> AddEntry(histYAxeRun2, "Run 2", "PL");
    legendYAxe -> AddEntry(histYAxeCut1, "Run 3 - matchedMchMid", "PL");
    legendYAxe -> AddEntry(histYAxeCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendYAxe -> AddEntry(histYAxeCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendYAxe -> AddEntry(histYAxeCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendYAxe -> Draw("SAME");


    canvasPtGenRec -> SaveAs("gen_rec_distributions_vs_pt.pdf");
    canvasYGenRec -> SaveAs("gen_rec_distributions_vs_y.pdf");
    canvasPtAxe -> SaveAs("Axe_distributions_vs_pt.pdf");
    canvasYAxe -> SaveAs("Axe_distributions_vs_y.pdf");
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