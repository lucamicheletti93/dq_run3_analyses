void LoadStyle();
void SetLegend(TLegend *, double);

void pileup(string triggerMask = "all") { // all, tvx, sel8
    LoadStyle();

    TFile *fInData = new TFile("AnalysisResults_data.root", "READ");
    TList *hlistStatisticsData = (TList*) fInData -> Get("table-maker/Statistics");
    TH2D *histBcStatsData = (TH2D*) hlistStatisticsData -> FindObject("BcStats");

    int binTriggerMaskData = histBcStatsData -> GetXaxis() -> FindBin(triggerMask.c_str());
    TH1D *histMuData = (TH1D*) histBcStatsData -> ProjectionY("histMuData", binTriggerMaskData, binTriggerMaskData);
    histMuData -> Scale(1. / histMuData -> Integral());
    histMuData -> SetLineColor(kBlack);
    histMuData -> SetLineWidth(2);

    TFile *fInMc = new TFile("AnalysisResults_mc.root", "READ");
    TList *hlistStatisticsMc = (TList*) fInMc -> Get("table-maker/Statistics");
    TH2D *histBcStatsMc = (TH2D*) hlistStatisticsMc -> FindObject("BcStats");

    int binTriggerMaskMc = histBcStatsMc -> GetXaxis() -> FindBin(triggerMask.c_str());
    TH1D *histMuMc = (TH1D*) histBcStatsMc -> ProjectionY("histMuMc", binTriggerMaskMc, binTriggerMaskMc);
    histMuMc -> Scale(1. / histMuMc -> Integral());
    histMuMc -> SetLineColor(kRed+1);
    histMuMc -> SetMarkerStyle(20);
    histMuMc -> SetMarkerSize(0.5);
    histMuMc -> SetMarkerColor(kRed+1);
    histMuMc -> SetLineWidth(1);

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.040);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);

    TCanvas *canvasMu = new TCanvas("canvasMu", "", 800, 600);
    gStyle -> SetOptStat(false);
    histMuData -> SetTitle("");
    histMuData -> GetXaxis() -> SetRangeUser(0, 0.120);
    histMuData -> GetXaxis() -> SetLabelSize(0.045);
    histMuData -> GetXaxis() -> SetTitleSize(0.045);
    histMuData -> GetYaxis() -> SetRangeUser(0, 0.025);
    histMuData -> GetYaxis() -> SetLabelSize(0.045);
    histMuData -> GetYaxis() -> SetTitleSize(0.045);
    histMuData -> Draw("H");
    histMuMc -> Draw("EP SAME");

    TLegend *legendMuData = new TLegend(0.20, 0.85, 0.45, 0.93, " ", "brNDC");
    SetLegend(legendMuData, 0.040);
    legendMuData -> AddEntry(histMuData, "Data (LHC25ae_pass2)", "L");
    legendMuData -> Draw();

    TLegend *legendMuMc = new TLegend(0.20, 0.65, 0.45, 0.73, " ", "brNDC");
    SetLegend(legendMuMc, 0.040);
    legendMuMc -> AddEntry(histMuMc, "MC (LHC25i4)", "L");
    legendMuMc -> Draw();

    double muData = histMuData -> GetMean();
    double corrFactorData = (muData) / (1 - TMath::Exp(-muData));

    double muMc = histMuMc -> GetMean();
    double corrFactorMc = (muMc) / (1 - TMath::Exp(-muMc));

    latexTitle -> DrawLatex(0.22, 0.80, Form("<#mu>_{Data} = %5.4f", muData));
    latexTitle -> DrawLatex(0.22, 0.74, Form("<#it{N}_{coll}^{TVX} | #it{N}_{coll}^{TVX} #geq 1>_{Data} = %5.4f", corrFactorData));
    latexTitle -> DrawLatex(0.22, 0.60, Form("<#mu>_{MC}   = %5.4f", muMc));
    latexTitle -> DrawLatex(0.22, 0.54, Form("<#it{N}_{coll}^{TVX} | #it{N}_{coll}^{TVX} #geq 1>_{MC} = %5.4f", corrFactorMc));


    canvasMu -> SaveAs("muDistribution.pdf");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend, double textSize = 0.040){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(textSize);
}