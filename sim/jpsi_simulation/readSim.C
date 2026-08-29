void LoadStyle();
void SetLegend(TLegend *, double);
inline void SetHist(auto *hist, Color_t mkrCol = kBlack, int mkrSty = 20, double mkrSize = 1, Color_t lnCol = kBlack, int lnWidth = 1, Color_t fillCol = kBlack, int fillSty = 0, double alpha = 1.) {
    hist -> SetMarkerColorAlpha(mkrCol, 1);
    hist -> SetMarkerStyle(mkrSty);
    hist -> SetMarkerSize(mkrSize);
    hist -> SetLineColorAlpha(lnCol, 1);
    hist -> SetLineWidth(lnWidth);
    hist -> SetFillStyle(fillSty);
    hist -> SetFillColorAlpha(fillCol, alpha);
}

double PtJPsipp13TeV(double*, double*);
double YJPsipp13TeV(double*, double*);
double PtPsi2spp13TeV(double*, double*);
double YPsi2spp13TeV(double*, double*);

void readSim() {
    LoadStyle();


    float pTOnia, yOnia, etaOnia, phiOnia, fromB, absPdg, nchMid, nchForward;

    TFile *fIn = TFile::Open("pythia8_onia_kMonash_kSoftQCD.root");
    TNtuple *ntuple = (TNtuple*) fIn->Get("tuplePairs");
    ntuple->Print();
    ntuple->SetBranchAddress("pTOnia", &pTOnia);
    ntuple->SetBranchAddress("yOnia", &yOnia);
    ntuple->SetBranchAddress("etaOnia", &etaOnia);
    ntuple->SetBranchAddress("phiOnia", &phiOnia);
    ntuple->SetBranchAddress("fromB", &fromB);
    ntuple->SetBranchAddress("absPdg", &absPdg);
    ntuple->SetBranchAddress("nchMid", &nchMid);
    ntuple->SetBranchAddress("nchForward", &nchForward);

    TH1F *hPtJpsi = new TH1F("hPtJpsi","", 100, 0, 20);
    TH1F *hPtPsi2s = new TH1F("hPtPsi2s","", 100, 0, 20);

    TH1F *hNchMidJpsi = new TH1F("hNchMidJpsi","", 10, 0, 100);
    TH1F *hNchMidPsi2s = new TH1F("hNchMidPsi2s","", 10, 0, 100);
    TH1F *hNchMidUpsilon1s = new TH1F("hNchMidUpsilon1s","", 10, 0, 100);

    TH1F *hNchForwardJpsi = new TH1F("hNchForwardJpsi","", 10, 0, 100);
    TH1F *hNchForwardPsi2s = new TH1F("hNchForwardPsi2s","", 10, 0, 100);

    Long64_t nEntries = ntuple->GetEntries();
    for (int iEntry = 0;iEntry < nEntries;iEntry++) {
        ntuple->GetEntry(iEntry);
        if (pTOnia > 20) continue;
        if (yOnia < 2.5 && yOnia > 4) continue;
        if (absPdg == 443) {
            hPtJpsi->Fill(pTOnia);
            hNchMidJpsi->Fill(nchMid);
            hNchForwardJpsi->Fill(nchForward);
        }
        if (absPdg == 100443) {
            hPtPsi2s->Fill(pTOnia);
            hNchMidPsi2s->Fill(nchMid);
            hNchForwardPsi2s->Fill(nchForward);
        }
        if (absPdg == 553) {
            hNchMidUpsilon1s->Fill(nchForward);
        }
    }

    TH1F *hNchMidRatio = (TH1F*) hNchMidPsi2s->Clone("hNchMidRatio");
    hNchMidRatio->Sumw2();
    hNchMidRatio->Divide(hNchMidJpsi);
    hNchMidRatio->SetLineColor(kRed);

    TH1F *hNchForwardRatio = (TH1F*) hNchForwardPsi2s->Clone("hNchForwardRatio");
    hNchForwardRatio->Sumw2();
    hNchForwardRatio->Divide(hNchForwardJpsi);
    hNchForwardRatio->SetLineColor(kBlue);

    TCanvas *canvasPt = new TCanvas("canvasPt", "", 800, 600);
    hPtJpsi->Draw("EP");
    hPtPsi2s->Draw("EP SAME");

    TCanvas *canvasNchMid = new TCanvas("canvasNchMid", "", 800, 600);
    gPad->SetLogy(true);
    hNchMidJpsi->Draw("EP");
    hNchMidPsi2s->Draw("EP SAME");
    hNchMidUpsilon1s->Draw("EP SAME");

    TCanvas *canvasRatioNchMid = new TCanvas("canvasRatioNchMid", "", 800, 600);
    hNchMidRatio->Draw("EP");
    hNchForwardRatio->Draw("EP SAME");


    
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double PtJPsipp13TeV(double* x, double* pars) {
    // jpsi y in pp at 13 TeV, tuned on data (2015)
    double xx = x[0];
    float p1, p2, p3;
    p1 = 4.75208;
    p2 = 1.69247;
    p3 = 4.49224;
    return pars[0] * xx / TMath::Power(1. + TMath::Power(xx / p1, p2), p3);
}
//-------------------------------------------------------------------------//
double YJPsipp13TeV(const double* x, const double* pars) {
    // jpsi y in pp at 13 TeV, tuned on data (2015)
    double xx = x[0];
    float p1, p2;
    p1 = 0;
    p2 = 2.98887;
    return pars[0] * TMath::Exp(-(1. / 2.) * TMath::Power(((xx - p1) / p2), 2));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double PtPsi2spp13TeV(double* x, double* pars) {
    // jpsi y in pp at 13 TeV, tuned on data (2015)
    double xx = x[0];
    float p1, p2, p3;
    p1 = 4.75208;
    p2 = 1.69247;
    p3 = 4.49224;
    return pars[0] * xx / TMath::Power(1. + TMath::Power(xx / p1, p2), p3);
}
//-------------------------------------------------------------------------//
double YPsi2spp13TeV(const double* x, const double* pars) {
    // jpsi y in pp at 13 TeV, tuned on data (2015)
    double xx = x[0];
    float p1, p2;
    p1 = 0;
    p2 = 2.98887;
    return pars[0] * TMath::Exp(-(1. / 2.) * TMath::Power(((xx - p1) / p2), 2));
}