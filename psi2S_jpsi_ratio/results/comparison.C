void LoadStyle();
void SetLegend(TLegend *);

void comparison(bool normToIntegral = true) {
    LoadStyle();

    // Run 2 results: https://alice-notes.web.cern.ch/system/files/notes/analysis/497/2017-Aug-11-analysis_note-pp13TeV-analysis_note.pdf
    // (0.0, 0.5), sig (12896, 30703), stat (188, 295), syst (333, 752)
    double ptBinsJpsiRun2[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0};
    double sigJpsiPtRun2[] = {43599, 71984, 62544, 48123, 35498, 24109, 15957, 10366, 10457, 4487, 2006, 960};
    double statSigJpsiPtRun2[] = {483, 435, 391, 328, 268, 206, 173, 134, 139, 93, 62, 48};
    double systSigJpsiPtRun2[] = {1097, 1853, 1883, 1474, 970, 578, 337, 232, 239, 104, 56, 23};

    double ptBinsPsi2sRun2[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0};
    double sigPsi2sPtRun2[] = {547, 1922, 1755, 1072, 829, 615, 384, 288, 281, 123, 88, 45};
    double statSigPsi2sPtRun2[] = {116, 162, 137, 105, 81, 63, 48, 40, 41, 29, 19, 16};
    double systSigPsi2sPtRun2[] = {44, 95, 111, 64, 52, 34, 27, 20, 43, 17, 9, 4};

    TH1D *histSigJpsiPtRun2 = new TH1D("histSigJpsiPtRun2", "", 12, ptBinsJpsiRun2);
    histSigJpsiPtRun2 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histSigJpsiPtRun2 -> SetLineWidth(2);
    histSigJpsiPtRun2 -> SetLineColor(kBlack);
    histSigJpsiPtRun2 -> SetMarkerColor(kBlack);

    TH1D *histSigPsi2sPtRun2 = new TH1D("histSigPsi2sPtRun2", "", 12, ptBinsPsi2sRun2);
    histSigPsi2sPtRun2 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histSigPsi2sPtRun2 -> SetLineWidth(2);
    histSigPsi2sPtRun2 -> SetLineColor(kBlack);
    histSigPsi2sPtRun2 -> SetMarkerColor(kBlack);

    for (int i = 0;i < 12;i++) {
        histSigJpsiPtRun2 -> SetBinContent(i+1, sigJpsiPtRun2[i]);
        histSigJpsiPtRun2 -> SetBinError(i+1, TMath::Sqrt(statSigJpsiPtRun2[i]*statSigJpsiPtRun2[i] + systSigJpsiPtRun2[i]*systSigJpsiPtRun2[i]));
        histSigPsi2sPtRun2 -> SetBinContent(i+1, sigPsi2sPtRun2[i]);
        histSigPsi2sPtRun2 -> SetBinError(i+1, TMath::Sqrt(statSigPsi2sPtRun2[i]*statSigPsi2sPtRun2[i] + systSigPsi2sPtRun2[i]*systSigPsi2sPtRun2[i]));
    }

    TH1D *histRatioPsi2sOverJpsiPtRun2 = (TH1D*) histSigPsi2sPtRun2 -> Clone("histRatioPsi2sOverJpsiPtRun2");
    histRatioPsi2sOverJpsiPtRun2 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histRatioPsi2sOverJpsiPtRun2 -> Divide(histSigJpsiPtRun2);
    histRatioPsi2sOverJpsiPtRun2 -> SetLineWidth(2);
    histRatioPsi2sOverJpsiPtRun2 -> SetLineColor(kBlack);
    histRatioPsi2sOverJpsiPtRun2 -> SetMarkerColor(kBlack);

    // Run 3 preliminary results: https://alice-notes.web.cern.ch/system/files/notes/analysis/1478/2023-11-10-Ratio_Psi2s_to_Jpsi___Analysis_Note-2.pdf
    double ptBinsJpsiPrelRun3[] = {0.00, 1.00, 2.00, 3.00, 4.00, 5.00, 7.00, 10.00, 20.00};
    double sigJpsiPtPrelRun3[] = {199646, 387021, 339776, 239877, 155224, 157660, 73202, 27917};
    double statSigJpsiPtPrelRun3[] = {1023, 1531, 1483, 1282, 904, 798, 463, 289};
    double systSigJpsiPtPrelRun3[] = {11430, 24128, 20357, 16037, 9913, 9569, 4331, 1816};

    double ptBinsPsi2sPrelRun3[] = {0.00, 1.00, 2.00, 3.00, 4.00, 5.00, 7.00, 10.00, 20.00};
    double sigPsi2sPtPrelRun3[] = {3341, 7369, 6516, 6779, 4266, 4760, 2231, 1148};
    double statSigPsi2sPtPrelRun3[] = {329, 507, 464, 440, 338, 306, 235, 136};
    double systSigPsi2sPtPrelRun3[] = {887, 1809, 1645, 1522, 998, 1170, 661, 257};

    double ratioPsi2sOverJpsiPtPrelRun3[] = {0.0166, 0.0188, 0.0190, 0.0280, 0.0272, 0.0299, 0.0301, 0.0408};
    double statRatioPsi2sOverJpsiPtPrelRun3[] = {0.0017, 0.0013, 0.0014, 0.0018, 0.0022, 0.0020, 0.0032, 0.0049};
    double systRatioPsi2sOverJpsiPtPrelRun3[] = {0.0036, 0.0036, 0.0037, 0.0045, 0.0047, 0.0056, 0.0075, 0.0069};

    TH1D *histSigJpsiPtPrelRun3 = new TH1D("histSigJpsiPtPrelRun3", "", 8, ptBinsJpsiPrelRun3);
    histSigJpsiPtPrelRun3 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histSigJpsiPtPrelRun3 -> SetLineWidth(2);
    histSigJpsiPtPrelRun3 -> SetLineColor(kAzure+2);
    histSigJpsiPtPrelRun3 -> SetMarkerColor(kAzure+2);

    TH1D *histSigPsi2sPtPrelRun3 = new TH1D("histSigPsi2sPtPrelRun3", "", 8, ptBinsPsi2sPrelRun3);
    histSigPsi2sPtPrelRun3 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histSigPsi2sPtPrelRun3 -> SetLineWidth(2);
    histSigPsi2sPtPrelRun3 -> SetLineColor(kAzure+2);
    histSigPsi2sPtPrelRun3 -> SetMarkerColor(kAzure+2);

    TH1D *histRatioPsi2sOverJpsiPtPrelRun3 = new TH1D("histRatioPsi2sOverJpsiPtPrelRun3", "", 8, ptBinsPsi2sPrelRun3);
    histRatioPsi2sOverJpsiPtPrelRun3 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histRatioPsi2sOverJpsiPtPrelRun3 -> SetLineWidth(2);
    histRatioPsi2sOverJpsiPtPrelRun3 -> SetLineColor(kAzure+2);
    histRatioPsi2sOverJpsiPtPrelRun3 -> SetMarkerColor(kAzure+2);

    for (int i = 0;i < 8;i++) {
        histSigJpsiPtPrelRun3 -> SetBinContent(i+1, sigJpsiPtPrelRun3[i]);
        histSigJpsiPtPrelRun3 -> SetBinError(i+1, TMath::Sqrt(statSigJpsiPtPrelRun3[i]*statSigJpsiPtPrelRun3[i] + systSigJpsiPtPrelRun3[i]*systSigJpsiPtPrelRun3[i]));
        histSigPsi2sPtPrelRun3 -> SetBinContent(i+1, sigPsi2sPtPrelRun3[i]);
        histSigPsi2sPtPrelRun3 -> SetBinError(i+1, TMath::Sqrt(statSigPsi2sPtPrelRun3[i]*statSigPsi2sPtPrelRun3[i] + systSigPsi2sPtPrelRun3[i]*systSigPsi2sPtPrelRun3[i]));
        histRatioPsi2sOverJpsiPtPrelRun3 -> SetBinContent(i+1, ratioPsi2sOverJpsiPtPrelRun3[i]);
        histRatioPsi2sOverJpsiPtPrelRun3 -> SetBinError(i+1, TMath::Sqrt(statRatioPsi2sOverJpsiPtPrelRun3[i]*statRatioPsi2sOverJpsiPtPrelRun3[i] + systRatioPsi2sOverJpsiPtPrelRun3[i]*systRatioPsi2sOverJpsiPtPrelRun3[i]));
    }

    if (normToIntegral) {
        std::cout << "------- Yield histograms normalized to the integral and bin width -------" << std::endl;
        histSigJpsiPtRun2 -> Scale(1. / histSigJpsiPtRun2 -> Integral(), "WIDTH");
        histSigJpsiPtPrelRun3 -> Scale(1. / histSigJpsiPtPrelRun3 -> Integral(), "WIDTH");
        histSigPsi2sPtRun2 -> Scale(1. / histSigPsi2sPtRun2 -> Integral(), "WIDTH");
        histSigPsi2sPtPrelRun3 -> Scale(1. / histSigPsi2sPtPrelRun3 -> Integral(), "WIDTH");
    }

    TLegend *legend = new TLegend(0.60, 0.75, 0.80, 0.85);
    SetLegend(legend);
    legend -> AddEntry(histSigJpsiPtRun2, "Run 2", "PL");
    legend -> AddEntry(histSigJpsiPtPrelRun3, "Run 3 (Preliminary)", "PL");

    TCanvas *canvasSigJpsiPt = new TCanvas("canvasSigJpsiPt", "", 800, 600);
    gPad -> SetLogy(1);
    histSigJpsiPtRun2 -> Draw("EP");
    histSigJpsiPtPrelRun3 -> Draw("EP SAME");
    legend -> Draw("SAME");

    TCanvas *canvasSigPsi2sPt = new TCanvas("canvasSigPsi2sPt", "", 800, 600);
    gPad -> SetLogy(1);
    histSigPsi2sPtRun2 -> Draw("EP");
    histSigPsi2sPtPrelRun3 -> Draw("EP SAME");
    legend -> Draw("SAME");

    TCanvas *canvasRatioPsi2sOverJpsi = new TCanvas("canvasRatioPsi2sOverJpsi", "", 800, 600);
    histRatioPsi2sOverJpsiPtRun2 -> GetYaxis() -> SetRangeUser(0, 0.12);
    histRatioPsi2sOverJpsiPtRun2 -> Draw("EP");
    histRatioPsi2sOverJpsiPtPrelRun3 -> Draw("EP SAME");
    legend -> Draw("SAME");





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
    legend -> SetTextSize(0.045);
}