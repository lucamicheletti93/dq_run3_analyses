void LoadStyle();
void SetLegend(TLegend *);
inline void SetHist(auto *hist, Color_t mkrCol = kBlack, int mkrSty = 20, double mkrSize = 1, Color_t lnCol = kBlack, int lnWidth = 1, int fillSty = 0, double alpha = 1) {
    hist -> SetMarkerColorAlpha(mkrCol, alpha);
    hist -> SetMarkerStyle(mkrSty);
    hist -> SetMarkerSize(mkrSize);
    hist -> SetLineColorAlpha(lnCol, alpha);
    hist -> SetLineWidth(lnWidth);
    hist -> SetFillStyle(fillSty);
    hist -> SetFillColorAlpha(lnCol, alpha);
}

void readSim() {
    LoadStyle();

    ///////////////////////////////////////////////
    double ptCenters[] = {0.250,0.750,1.500,2.500,3.500,4.500,5.500,6.500,8.500,15.000};
    double ptWidths[] = {0.250,0.250,0.500,0.500,0.500,0.500,0.500,0.500,1.500,5.000};
    double jpsiXsecPtRap254[] = {562.776,1577.271,2273.604,2030.779,1263.990,689.372,406.291,223.584,72.302,8.290};
    double jpsiStatXsecPtRap254[] = {28.877,47.418,44.155,39.866,34.850,17.654,12.878,6.588,2.783,0.438};
    double jpsiSystXsecPtRap254[] = {12.444,19.252,6.688,9.108,6.792,8.046,2.923,4.711,0.071,0.090};

    TGraphErrors *graStatJpsiXsecPtRap254 = new TGraphErrors(10, ptCenters, jpsiXsecPtRap254, ptWidths, jpsiStatXsecPtRap254);
    SetHist(graStatJpsiXsecPtRap254, kBlack, 20, 1, kBlack);
    ///////////////////////////////////////////////

    TFile *fIn = new TFile("output/pythia8_onia_kMonash_kSoftQCD.root", "READ");
    TH1D *histNumEvents = (TH1D*) fIn -> Get("histNumEvents");
    TH1D *histXsec = (TH1D*) fIn -> Get("histXsec_0");

    double nEvents = histNumEvents -> GetBinContent(2);
    double xSec = histXsec -> GetBinContent(1) * 1.e6; // nb

    float pTOnia, yOnia, etaOnia, phiOnia, fromB, absPdg;

    TNtuple *ntuple = (TNtuple*) fIn -> Get("tuplePairs");
    ntuple -> SetBranchAddress("pTOnia", &pTOnia);
    ntuple -> SetBranchAddress("yOnia", &yOnia);
    ntuple -> SetBranchAddress("etaOnia", &etaOnia);
    ntuple -> SetBranchAddress("phiOnia", &phiOnia);
    ntuple -> SetBranchAddress("fromB", &fromB);
    ntuple -> SetBranchAddress("absPdg", &absPdg);

    double ptBinEdges[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0};

    TH1D *histXsecJpsiPt = new TH1D("histXsecJpsiPt", ";#it{p}_{T} (GeV/c);d#sigma/d#it{p}_{T} (GeV/c nb^{-1})^{-1}", 16, ptBinEdges);
    TH1D *histXsecJpsiRap = new TH1D("histXsecJpsiRap", ";y;d#sigma/dy (nb)", 20, 0, 5);
    TH1D *histXsecJpsiPtFwdCuts = new TH1D("histXsecJpsiPtFwdCuts", ";#it{p}_{T} (GeV/c);d#sigma/d#it{p}_{T} (GeV/c nb^{-1})^{-1}", 16, ptBinEdges);
    TH1D *histXsecJpsiRapFwdCuts = new TH1D("histXsecJpsiRapFwdCuts", ";y;d#sigma/dy (nb)", 20, 0, 5);
    TH1D *histXsecPsi2sPt = new TH1D("histXsecPsi2sPt", ";#it{p}_{T} (GeV/c);d#sigma/d#it{p}_{T} (GeV/c nb^{-1})^{-1}", 16, ptBinEdges);
    TH1D *histXsecPsi2sRap = new TH1D("histXsecPsi2sRap", ";y;d#sigma/dy (nb)", 20, 0, 5);
    TH1D *histXsecPsi2sPtFwdCuts = new TH1D("histXsecPsi2sPtFwdCuts", ";#it{p}_{T} (GeV/c);d#sigma/d#it{p}_{T} (GeV/c nb^{-1})^{-1}", 16, ptBinEdges);
    TH1D *histXsecPsi2sRapFwdCuts = new TH1D("histXsecPsi2sRapFwdCuts", ";y;d#sigma/dy (nb)", 20, 0, 5);

    SetHist(histXsecJpsiPt, kRed+1, 20, 1, kRed+1, 3004);
    SetHist(histXsecJpsiRap, kRed+1, 20, 1, kRed+1);
    SetHist(histXsecJpsiPtFwdCuts, kRed+1, 24, 1, kRed+1);
    SetHist(histXsecJpsiRapFwdCuts, kRed+1, 24, 1, kRed+1);
    SetHist(histXsecPsi2sPt, kAzure+4, 20, 1, kAzure+4);
    SetHist(histXsecPsi2sRap, kAzure+4, 20, 1, kAzure+4);
    SetHist(histXsecPsi2sPtFwdCuts, kAzure+4, 24, 1, kAzure+4);
    SetHist(histXsecPsi2sRapFwdCuts, kAzure+4, 24, 1, kAzure+4);

    histXsecJpsiPt -> Sumw2(true);
    histXsecJpsiRap -> Sumw2(true);
    histXsecJpsiPtFwdCuts -> Sumw2(true);
    histXsecJpsiRapFwdCuts -> Sumw2(true);
    histXsecPsi2sPt -> Sumw2(true);
    histXsecPsi2sRap -> Sumw2(true);
    histXsecPsi2sPtFwdCuts -> Sumw2(true);
    histXsecPsi2sRapFwdCuts -> Sumw2(true);

    Long64_t nEntries = ntuple -> GetEntries();
    for (Long64_t i = 0;i < nEntries;++i) {
        ntuple -> GetEntry(i);
        if (absPdg == 443) {
            histXsecJpsiPt -> Fill(pTOnia);
            histXsecJpsiRap -> Fill(yOnia);

            if (yOnia > 2.5 && yOnia < 4) {
                histXsecJpsiPtFwdCuts -> Fill(pTOnia);
                histXsecJpsiRapFwdCuts -> Fill(yOnia);
            }
        }

        if (absPdg == 100443) {
            histXsecPsi2sPt -> Fill(pTOnia);
            histXsecPsi2sRap -> Fill(yOnia);

            if (yOnia > 2.5 && yOnia < 4) {
                histXsecPsi2sPtFwdCuts -> Fill(pTOnia);
                histXsecPsi2sRapFwdCuts -> Fill(yOnia);
            }
        }
    }

    TH1D *histPsi2sOverJpsiPt = (TH1D*) histXsecPsi2sPt -> Clone("histPsi2sOverJpsiPt");
    histPsi2sOverJpsiPt -> Divide(histXsecJpsiPt);
    histPsi2sOverJpsiPt -> GetYaxis() -> SetTitle("#psi(2S) / J/#psi");
    SetHist(histPsi2sOverJpsiPt, kBlack, 20, 1, kBlack);

    TH1D *histPsi2sOverJpsiRap = (TH1D*) histXsecPsi2sRap -> Clone("histPsi2sOverJpsiRap");
    histPsi2sOverJpsiRap -> Divide(histXsecJpsiRap);
    histPsi2sOverJpsiRap -> GetYaxis() -> SetTitle("#psi(2S) / J/#psi");
    SetHist(histPsi2sOverJpsiRap, kBlack, 20, 1, kBlack);

    TH1D *histPsi2sOverJpsiPtFwdCuts = (TH1D*) histXsecPsi2sPtFwdCuts -> Clone("histPsi2sOverJpsiPtFwdCuts");
    histPsi2sOverJpsiPtFwdCuts -> Divide(histXsecJpsiPtFwdCuts);
    histPsi2sOverJpsiPtFwdCuts -> GetYaxis() -> SetTitle("#psi(2S) / J/#psi");
    SetHist(histPsi2sOverJpsiPtFwdCuts, kBlack, 24, 1, kBlack);

    TH1D *histPsi2sOverJpsiRapFwdCuts = (TH1D*) histXsecPsi2sRapFwdCuts -> Clone("histPsi2sOverJpsiRapFwdCuts");
    histPsi2sOverJpsiRapFwdCuts -> Divide(histXsecPsi2sRapFwdCuts);
    histPsi2sOverJpsiRapFwdCuts -> GetYaxis() -> SetTitle("#psi(2S) / J/#psi");
    SetHist(histPsi2sOverJpsiRapFwdCuts, kBlack, 24, 1, kBlack);

    histXsecJpsiPt -> Scale(xSec / nEvents, "WIDTH");
    histXsecPsi2sPt -> Scale(xSec / nEvents, "WIDTH");
    histXsecJpsiPtFwdCuts -> Scale(xSec / nEvents, "WIDTH");
    histXsecPsi2sPtFwdCuts -> Scale(xSec / nEvents, "WIDTH");

    histXsecJpsiRap -> Scale(xSec / nEvents);
    histXsecPsi2sRap -> Scale(xSec / nEvents);
    histXsecJpsiRapFwdCuts -> Scale(xSec / nEvents);
    histXsecPsi2sRapFwdCuts -> Scale(xSec / nEvents);

    TCanvas *canvasOniaPt = new TCanvas("canvasOniaPt", "", 800, 600);
    gPad -> SetLogy(true);
    histXsecJpsiPt -> GetYaxis() -> SetRangeUser(1, 2 * histXsecJpsiPt -> GetMaximum());
    histXsecJpsiPt -> Draw("E2P");
    histXsecPsi2sPt -> Draw("E2P SAME");
    histXsecJpsiPtFwdCuts -> Draw("E2P SAME");
    histXsecPsi2sPtFwdCuts -> Draw("E2P SAME");
    graStatJpsiXsecPtRap254 -> Draw("EP SAME");

    TLegend *legendOniaPt = new TLegend(0.65, 0.70, 0.80, 0.93, " ", "brNDC");
    SetLegend(legendOniaPt);
    legendOniaPt -> AddEntry(histXsecJpsiPt, "J/#psi, all", "EP");
    legendOniaPt -> AddEntry(histXsecPsi2sPt, "#psi(2S), all", "EP");
    legendOniaPt -> AddEntry(histXsecJpsiPtFwdCuts, "J/#psi, 2.5 < y < 4", "EP");
    legendOniaPt -> AddEntry(histXsecPsi2sPtFwdCuts, "#psi(2S), 2.5 < y < 4", "EP");
    legendOniaPt -> Draw();

    return;

    TCanvas *canvasOniaRatioPt = new TCanvas("canvasOniaRatioPt", "", 800, 600);
    histPsi2sOverJpsiPt -> GetYaxis() -> SetRangeUser(0, 0.4);
    histPsi2sOverJpsiPt -> Draw("EP");
    histPsi2sOverJpsiPtFwdCuts -> Draw("EP SAME");

    TLegend *legendOniaRatioPt = new TLegend(0.20, 0.70, 0.40, 0.85, " ", "brNDC");
    SetLegend(legendOniaRatioPt);
    legendOniaRatioPt -> AddEntry(histPsi2sOverJpsiPt, "all", "EP");
    legendOniaRatioPt -> AddEntry(histPsi2sOverJpsiPtFwdCuts, "2.5 < y < 4", "EP");
    legendOniaRatioPt -> Draw();

    TCanvas *canvasOniaRap = new TCanvas("canvasOniaRap", "", 800, 600);
    gPad -> SetLogy(true);
    histXsecJpsiRap -> GetYaxis() -> SetRangeUser(0.01, 2 * histXsecJpsiRap -> GetMaximum());
    histXsecJpsiRap -> Draw("EP");
    histXsecPsi2sRap -> Draw("EP SAME");
    histXsecJpsiRapFwdCuts -> Draw("EP SAME");
    histXsecPsi2sRapFwdCuts -> Draw("EP SAME");

    TCanvas *canvasOniaRatioRap = new TCanvas("canvasOniaRatioRap", "", 800, 600);
    histPsi2sOverJpsiRap -> GetYaxis() -> SetRangeUser(0, 0.3);
    histPsi2sOverJpsiRap -> Draw("EP");
    histPsi2sOverJpsiRapFwdCuts -> Draw("EP SAME");

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
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.045);
}