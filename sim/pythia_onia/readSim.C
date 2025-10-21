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

void readSim(string simPath = "sim_pythia_onia_13.6TeV", string simName = "pythia8_onia_kMonash_kSoftQCD") {
    LoadStyle();

    ///////////////////////////////////////////////
    TFile *fInFonllXsecNonPromptJpsi = new TFile("FONLL_Jpsi.root", "READ");
    TH1D *histFonllXsecNonPromptJpsi = (TH1D*) fInFonllXsecNonPromptJpsi -> Get("histXsec");
    SetHist(histFonllXsecNonPromptJpsi, kGreen+1, 20, 0, kGreen+1, 1, kGreen+1, 3352, 0.7);
    histFonllXsecNonPromptJpsi -> Scale(1 / 1e6, "WIDTH");

    TFile *fInFonllXsecNonPromptPsi2s = new TFile("FONLL_Psi2s.root", "READ");
    TH1D *histFonllXsecNonPromptPsi2s = (TH1D*) fInFonllXsecNonPromptPsi2s -> Get("histXsec");
    SetHist(histFonllXsecNonPromptPsi2s, kGreen+1, 20, 0, kGreen+1, 1, kGreen+1, 3352, 0.7);
    histFonllXsecNonPromptPsi2s -> Scale(1 / 1e6, "WIDTH");

    double ptCenters[] = {0.250,0.750,1.500,2.500,3.500,4.500,5.500,6.500,8.500,15.000};
    double ptWidths[] = {0.250,0.250,0.500,0.500,0.500,0.500,0.500,0.500,1.500,5.000};
    double jpsiXsecPtRap254[] = {562.776,1577.271,2273.604,2030.779,1263.990,689.372,406.291,223.584,72.302,8.290};
    double jpsiStatXsecPtRap254[] = {28.877,47.418,44.155,39.866,34.850,17.654,12.878,6.588,2.783,0.438};
    double jpsiSystXsecPtRap254[] = {12.444,19.252,6.688,9.108,6.792,8.046,2.923,4.711,0.071,0.090};

    double ptCentersRatio[] = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5, 5.5, 6.5, 7.5, 9, 11, 13, 15, 18};
    double ptWidthsRatio[] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1, 2};
    double psi2sToJpsiRatioPtRap254[] = {0.09548, 0.1101, 0.1175, 0.1338, 0.1411, 0.1363, 0.1514, 0.1780, 0.2177, 0.2336, 0.2192, 0.2654, 0.2934, 0.3248, 0.3342, 0.4244};
    double psi2sToJpsiRatioStatPtRap254[] = {0.0072, 0.0045, 0.0041, 0.0041, 0.0046, 0.0047, 0.004, 0.0049, 0.0058, 0.0076, 0.01, 0.0095, 0.015, 0.029, 0.039, 0.049};
    double psi2sToJpsiRatioSystPtRap254[] = {0.0047, 0.0059, 0.0061, 0.0077, 0.0089, 0.01, 0.014, 0.015, 0.0091, 0.012, 0.012, 0.019, 0.031, 0.022, 0.024, 0.054};

    // Include 10% Lumi + 0.5% BR
    for (int iPt = 0;iPt < 10;iPt++) {
        double relSystSigExtr = jpsiSystXsecPtRap254[iPt] / jpsiXsecPtRap254[iPt];
        double relSystLumi = 0.100;
        double relSystBr = 0.005;
        jpsiSystXsecPtRap254[iPt] = jpsiXsecPtRap254[iPt] * TMath::Sqrt(relSystSigExtr*relSystSigExtr + relSystLumi*relSystLumi + relSystBr*relSystBr);
    }

    TGraphErrors *graStatJpsiXsecPtRap254 = new TGraphErrors(10, ptCenters, jpsiXsecPtRap254, ptWidths, jpsiStatXsecPtRap254);
    SetHist(graStatJpsiXsecPtRap254, kBlack, 20, 1, kBlack);
    TGraphErrors *graSystJpsiXsecPtRap254 = new TGraphErrors(10, ptCenters, jpsiXsecPtRap254, ptWidths, jpsiSystXsecPtRap254);
    SetHist(graSystJpsiXsecPtRap254, kBlack, 20, 1, kBlack);

    TGraphErrors *graStatPsi2sOverJpsiPtRap254 = new TGraphErrors(16, ptCentersRatio, psi2sToJpsiRatioPtRap254, ptWidthsRatio, psi2sToJpsiRatioStatPtRap254);
    SetHist(graStatPsi2sOverJpsiPtRap254, kBlack, 20, 1, kBlack);
    TGraphErrors *graSystPsi2sOverJpsiPtRap254 = new TGraphErrors(16, ptCentersRatio, psi2sToJpsiRatioPtRap254, ptWidthsRatio, psi2sToJpsiRatioSystPtRap254);
    SetHist(graSystPsi2sOverJpsiPtRap254, kBlack, 20, 1, kBlack);
    ///////////////////////////////////////////////

    TFile *fIn = new TFile(Form("%s/%s.root", simPath.c_str(), simName.c_str()), "READ");
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

    const int nPtBins = 16;
    double ptBinEdges[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0};
    //const int nPtBins = 10;
    //double ptBinEdges[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0, 20.0};

    TH1D *histXsecJpsiPt = new TH1D("histXsecJpsiPt", ";#it{p}_{T} (GeV/c);d#sigma/d#it{p}_{T} (#mub / (GeV/c))", nPtBins, ptBinEdges);
    TH1D *histXsecJpsiRap = new TH1D("histXsecJpsiRap", ";y;d#sigma/dy (nb)", 20, 0, 5);
    TH1D *histXsecJpsiPtFwdCuts = new TH1D("histXsecJpsiPtFwdCuts", ";#it{p}_{T} (GeV/c);d#sigma/d#it{p}_{T} (#mub / (GeV/c))", nPtBins, ptBinEdges);
    TH1D *histXsecJpsiRapFwdCuts = new TH1D("histXsecJpsiRapFwdCuts", ";y;d#sigma/dy (nb)", 20, 0, 5);
    TH1D *histXsecPromptJpsiPtFwdCuts = new TH1D("histXsecPromptJpsiPtFwdCuts", ";#it{p}_{T} (GeV/c);d#sigma/d#it{p}_{T} (#mub / (GeV/c))", nPtBins, ptBinEdges);
    TH1D *histXsecPromptJpsiRapFwdCuts = new TH1D("histXsecPromptJpsiRapFwdCuts", ";y;d#sigma/dy (nb)", 20, 0, 5);
    TH1D *histXsecNonPromptJpsiPtFwdCuts = new TH1D("histXsecNonPromptJpsiPtFwdCuts", ";#it{p}_{T} (GeV/c);d#sigma/d#it{p}_{T} (#mub / (GeV/c))", nPtBins, ptBinEdges);
    TH1D *histXsecNonPromptJpsiRapFwdCuts = new TH1D("histXsecNonPromptJpsiRapFwdCuts", ";y;d#sigma/dy (nb)", 20, 0, 5);

    TH1D *histXsecPsi2sPt = new TH1D("histXsecPsi2sPt", ";#it{p}_{T} (GeV/c);d#sigma/d#it{p}_{T} (#mub / (GeV/c))", nPtBins, ptBinEdges);
    TH1D *histXsecPsi2sRap = new TH1D("histXsecPsi2sRap", ";y;d#sigma/dy (nb)", 20, 0, 5);
    TH1D *histXsecPsi2sPtFwdCuts = new TH1D("histXsecPsi2sPtFwdCuts", ";#it{p}_{T} (GeV/c);d#sigma/d#it{p}_{T} (#mub / (GeV/c))", nPtBins, ptBinEdges);
    TH1D *histXsecPsi2sRapFwdCuts = new TH1D("histXsecPsi2sRapFwdCuts", ";y;d#sigma/dy (nb)", 20, 0, 5);
    TH1D *histXsecPromptPsi2sPtFwdCuts = new TH1D("histXsecPromptPsi2sPtFwdCuts", ";#it{p}_{T} (GeV/c);d#sigma/d#it{p}_{T} (#mub / (GeV/c))", nPtBins, ptBinEdges);
    TH1D *histXsecPromptPsi2sRapFwdCuts = new TH1D("histXsecPromptPsi2sRapFwdCuts", ";y;d#sigma/dy (nb)", 20, 0, 5);
    TH1D *histXsecNonPromptPsi2sPtFwdCuts = new TH1D("histXsecNonPromptPsi2sPtFwdCuts", ";#it{p}_{T} (GeV/c);d#sigma/d#it{p}_{T} (#mub / (GeV/c))", nPtBins, ptBinEdges);
    TH1D *histXsecNonPromptPsi2sRapFwdCuts = new TH1D("histXsecNonPromptPsi2sRapFwdCuts", ";y;d#sigma/dy (nb)", 20, 0, 5);

    SetHist(histXsecJpsiPt, kRed+1, 24, 1, kRed+1, 1);
    SetHist(histXsecJpsiRap, kRed+1, 24, 1, kRed+1, 1);
    SetHist(histXsecJpsiPtFwdCuts, kRed+1, 20, 0.5, kRed+1, 1, kRed+1, 1001, 0.7);
    SetHist(histXsecJpsiRapFwdCuts, kRed+1, 20, 0.5, kRed+1, 1, kRed+1, 1001, 0.7);
    SetHist(histXsecPromptJpsiPtFwdCuts, kOrange+7, 20, 0.5, kOrange+7, 1, kOrange+7, 1001, 0.7);
    SetHist(histXsecPromptJpsiRapFwdCuts, kOrange+7, 20, 0.5, kOrange+7, 1, kOrange+7, 1001, 0.7);
    SetHist(histXsecNonPromptJpsiPtFwdCuts, kGray+1, 20, 0.5, kGray+1, 1, kGray+1, 1001, 0.7);
    SetHist(histXsecNonPromptJpsiRapFwdCuts, kGray+1, 20, 0.5, kGray+1, 1, kGray+1, 1001, 0.7);

    SetHist(histXsecPsi2sPt, kRed+1, 24, 1, kRed+1, 1);
    SetHist(histXsecPsi2sRap, kRed+1, 24, 1, kRed+1, 1);
    SetHist(histXsecPsi2sPtFwdCuts, kRed+1, 20, 0.5, kRed+1, 1, kRed+1, 1001, 0.7);
    SetHist(histXsecPsi2sRapFwdCuts, kRed+1, 20, 0.5, kRed+1, 1, kRed+1, 1001, 0.7);
    SetHist(histXsecPromptPsi2sPtFwdCuts, kOrange+7, 20, 0.5, kOrange+7, 1, kOrange+7, 1001, 0.7);
    SetHist(histXsecPromptPsi2sRapFwdCuts, kOrange+7, 20, 0.5, kOrange+7, 1, kOrange+7, 1001, 0.7);
    SetHist(histXsecNonPromptPsi2sPtFwdCuts, kGray+1, 20, 0.5, kGray+1, 1, kGray+1, 1001, 0.7);
    SetHist(histXsecNonPromptPsi2sRapFwdCuts, kGray+1, 20, 0.5, kGray+1, 1, kGray+1, 1001, 0.7);

    histXsecJpsiPt -> Sumw2(true);
    histXsecJpsiRap -> Sumw2(true);
    histXsecJpsiPtFwdCuts -> Sumw2(true);
    histXsecJpsiRapFwdCuts -> Sumw2(true);
    histXsecPromptJpsiPtFwdCuts -> Sumw2(true);
    histXsecPromptJpsiRapFwdCuts -> Sumw2(true);
    histXsecNonPromptJpsiPtFwdCuts -> Sumw2(true);
    histXsecNonPromptJpsiRapFwdCuts -> Sumw2(true);

    histXsecPsi2sPt -> Sumw2(true);
    histXsecPsi2sRap -> Sumw2(true);
    histXsecPsi2sPtFwdCuts -> Sumw2(true);
    histXsecPsi2sRapFwdCuts -> Sumw2(true);
    histXsecPromptPsi2sPtFwdCuts -> Sumw2(true);
    histXsecPromptPsi2sRapFwdCuts -> Sumw2(true);
    histXsecNonPromptPsi2sPtFwdCuts -> Sumw2(true);
    histXsecNonPromptPsi2sRapFwdCuts -> Sumw2(true);

    Long64_t nEntries = ntuple -> GetEntries();
    for (Long64_t i = 0;i < nEntries;++i) {
        ntuple -> GetEntry(i);
        if (absPdg == 443) {
            histXsecJpsiPt -> Fill(pTOnia);
            histXsecJpsiRap -> Fill(yOnia);

            if (yOnia > 2.5 && yOnia < 4 && pTOnia < 20) {
                histXsecJpsiPtFwdCuts -> Fill(pTOnia);
                histXsecJpsiRapFwdCuts -> Fill(yOnia);

                if (fromB) {
                    histXsecNonPromptJpsiPtFwdCuts -> Fill(pTOnia);
                    histXsecNonPromptJpsiRapFwdCuts -> Fill(yOnia);
                } else {
                    histXsecPromptJpsiPtFwdCuts -> Fill(pTOnia);
                    histXsecPromptJpsiRapFwdCuts -> Fill(yOnia);
                }
            }
        }

        if (absPdg == 100443) {
            histXsecPsi2sPt -> Fill(pTOnia);
            histXsecPsi2sRap -> Fill(yOnia);

            if (yOnia > 2.5 && yOnia < 4 && pTOnia < 20) {
                histXsecPsi2sPtFwdCuts -> Fill(pTOnia);
                histXsecPsi2sRapFwdCuts -> Fill(yOnia);

                if (fromB) {
                    histXsecNonPromptPsi2sPtFwdCuts -> Fill(pTOnia);
                    histXsecNonPromptPsi2sRapFwdCuts -> Fill(yOnia);
                } else {
                    histXsecPromptPsi2sPtFwdCuts -> Fill(pTOnia);
                    histXsecPromptPsi2sRapFwdCuts -> Fill(yOnia);
                }
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
    SetHist(histPsi2sOverJpsiPtFwdCuts, kRed+1, 20, 0, kRed+1, 1, kRed+1, 1001, 0.7);

    TH1D *histPromptPsi2sOverPromptJpsiPtFwdCuts = (TH1D*) histXsecPromptPsi2sPtFwdCuts -> Clone("histPromptPsi2sOverPromptJpsiPtFwdCuts");
    histPromptPsi2sOverPromptJpsiPtFwdCuts -> Divide(histXsecPromptJpsiPtFwdCuts);
    histPromptPsi2sOverPromptJpsiPtFwdCuts -> GetYaxis() -> SetTitle("Prompt #psi(2S) / Prompt J/#psi");
    SetHist(histPromptPsi2sOverPromptJpsiPtFwdCuts, kOrange+7, 20, 0, kOrange+7, 1, kOrange+7, 1001, 0.5);

    TH1D *histNonPromptPsi2sOverNonPromptJpsiPtFwdCuts = (TH1D*) histXsecNonPromptPsi2sPtFwdCuts -> Clone("histNonPromptPsi2sOverNonPromptJpsiPtFwdCuts");
    histNonPromptPsi2sOverNonPromptJpsiPtFwdCuts -> Divide(histXsecNonPromptJpsiPtFwdCuts);
    histNonPromptPsi2sOverNonPromptJpsiPtFwdCuts -> GetYaxis() -> SetTitle("Non Prompt #psi(2S) / Non Prompt J/#psi");
    SetHist(histNonPromptPsi2sOverNonPromptJpsiPtFwdCuts, kGray+1, 20, 0, kGray+1, 1, kGray+1, 1001, 0.5);

    TH1D *histPsi2sOverJpsiRapFwdCuts = (TH1D*) histXsecPsi2sRapFwdCuts -> Clone("histPsi2sOverJpsiRapFwdCuts");
    histPsi2sOverJpsiRapFwdCuts -> Divide(histXsecPsi2sRapFwdCuts);
    histPsi2sOverJpsiRapFwdCuts -> GetYaxis() -> SetTitle("#psi(2S) / J/#psi");
    SetHist(histPsi2sOverJpsiRapFwdCuts, kBlack, 24, 1, kBlack);

    histXsecJpsiPt -> Scale(xSec / (nEvents * 1e3), "WIDTH");
    histXsecPsi2sPt -> Scale(xSec / (nEvents * 1e3), "WIDTH");
    histXsecJpsiPtFwdCuts -> Scale(xSec / (nEvents * 1e3), "WIDTH");
    histXsecPsi2sPtFwdCuts -> Scale(xSec / (nEvents * 1e3), "WIDTH");
    histXsecPromptJpsiPtFwdCuts -> Scale(xSec / (nEvents * 1e3), "WIDTH");
    histXsecPromptPsi2sPtFwdCuts -> Scale(xSec / (nEvents * 1e3), "WIDTH");
    histXsecNonPromptJpsiPtFwdCuts -> Scale(xSec / (nEvents * 1e3), "WIDTH");
    histXsecNonPromptPsi2sPtFwdCuts -> Scale(xSec / (nEvents * 1e3), "WIDTH");

    histXsecJpsiRap -> Scale(xSec / (nEvents * 1e3));
    histXsecPsi2sRap -> Scale(xSec / (nEvents * 1e3));
    histXsecJpsiRapFwdCuts -> Scale(xSec / (nEvents * 1e3));
    histXsecPsi2sRapFwdCuts -> Scale(xSec / (nEvents * 1e3));
    histXsecPromptJpsiRapFwdCuts -> Scale(xSec / (nEvents * 1e3));
    histXsecPromptPsi2sRapFwdCuts -> Scale(xSec / (nEvents * 1e3));
    histXsecNonPromptJpsiRapFwdCuts -> Scale(xSec / (nEvents * 1e3));
    histXsecNonPromptPsi2sRapFwdCuts -> Scale(xSec / (nEvents * 1e3));

    TF1 *funcPtJpsi = new TF1("funcPtJpsi", PtJPsipp13TeV, 0, 20, 1);
    funcPtJpsi -> SetLineColor(kRed+1);
    funcPtJpsi -> SetLineStyle(kDashed);
    histXsecJpsiPtFwdCuts -> Fit(funcPtJpsi, "RQ0");
    TH1D *histFuncPtJpsi = (TH1D*) funcPtJpsi -> GetHistogram();

    TF1 *funcPtPsi2s = new TF1("funcPtPsi2s", PtPsi2spp13TeV, 0, 20, 1);
    funcPtPsi2s -> SetLineColor(kAzure+4);
    funcPtPsi2s -> SetLineStyle(kDashed);
    histXsecPsi2sPtFwdCuts -> Fit(funcPtPsi2s, "RQ0");
    TH1D *histFuncPtPsi2s = (TH1D*) funcPtPsi2s -> GetHistogram();


    TCanvas *canvasJpsiOniaPt = new TCanvas("canvasJpsiOniaPt", "", 800, 600);
    gPad -> SetLogy(true);
    gStyle -> SetHatchesSpacing(0.2);
    histXsecJpsiPtFwdCuts -> GetYaxis() -> SetRangeUser(1e-4, 10);
    histXsecJpsiPtFwdCuts -> Draw("E2P");
    histXsecPromptJpsiPtFwdCuts -> Draw("E2P SAME");
    histXsecNonPromptJpsiPtFwdCuts -> Draw("E2P SAME");
    histFonllXsecNonPromptJpsi -> Draw("E2P SAME");
    
    TLegend *legendJpsiOniaPt = new TLegend(0.55, 0.70, 0.75, 0.93, " ", "brNDC");
    SetLegend(legendJpsiOniaPt, 0.040);
    legendJpsiOniaPt -> AddEntry(histXsecJpsiPtFwdCuts, "Inclusive J/#psi - PYTHIA", "F");
    legendJpsiOniaPt -> AddEntry(histXsecPromptJpsiPtFwdCuts, "Prompt J/#psi - PYTHIA", "F");
    legendJpsiOniaPt -> AddEntry(histXsecNonPromptJpsiPtFwdCuts, "Non-Prompt J/#psi - PYTHIA", "F");
    legendJpsiOniaPt -> AddEntry(histFonllXsecNonPromptJpsi, "Non-Prompt J/#psi - FONLL", "F");
    legendJpsiOniaPt -> Draw();


    TCanvas *canvasPsi2sOniaPt = new TCanvas("canvasPsi2sOniaPt", "", 800, 600);
    gPad -> SetLogy(true);
    gStyle -> SetHatchesSpacing(0.2);
    histXsecPsi2sPtFwdCuts -> GetYaxis() -> SetRangeUser(1e-4, 1);
    histXsecPsi2sPtFwdCuts -> Draw("E2P");
    histXsecPromptPsi2sPtFwdCuts -> Draw("E2P SAME");
    histXsecNonPromptPsi2sPtFwdCuts -> Draw("E2P SAME");
    histFonllXsecNonPromptPsi2s -> Draw("E2P SAME");

    TLegend *legendPsi2sOniaPt = new TLegend(0.55, 0.70, 0.75, 0.93, " ", "brNDC");
    SetLegend(legendPsi2sOniaPt, 0.04);
    legendPsi2sOniaPt -> AddEntry(histXsecPsi2sPtFwdCuts, "Inclusive #psi(2S) - PYTHIA", "F");
    legendPsi2sOniaPt -> AddEntry(histXsecPromptPsi2sPtFwdCuts, "Prompt #psi(2S) - PYTHIA", "F");
    legendPsi2sOniaPt -> AddEntry(histXsecNonPromptPsi2sPtFwdCuts, "Non-Prompt #psi(2S) - PYTHIA", "F");
    legendPsi2sOniaPt -> AddEntry(histFonllXsecNonPromptPsi2s, "Non-Prompt #psi(2S) - FONLL", "F");
    legendPsi2sOniaPt -> Draw();
    

    /*TCanvas *canvasOniaPtVsData = new TCanvas("canvasOniaPtVsData", "", 800, 600);
    gPad -> SetLogy(true);
    gStyle -> SetHatchesSpacing(0.2);
    histXsecJpsiPtFwdCuts -> GetYaxis() -> SetRangeUser(1, 2 * histXsecJpsiPtFwdCuts -> GetMaximum());
    //histXsecJpsiPt -> Draw("E2P");
    //histXsecPsi2sPt -> Draw("E2P SAME");
    histXsecJpsiPtFwdCuts -> Draw("E2P");
    histXsecPsi2sPtFwdCuts -> Draw("E2P SAME");

    histXsecPromptJpsiPtFwdCuts -> Draw("E2P SAME");
    histXsecPromptPsi2sPtFwdCuts -> Draw("E2P SAME");

    histXsecNonPromptJpsiPtFwdCuts -> Draw("E2P SAME");
    histXsecNonPromptPsi2sPtFwdCuts -> Draw("E2P SAME");

    graStatJpsiXsecPtRap254 -> Draw("EP SAME");
    graSystJpsiXsecPtRap254 -> Draw("E2P SAME");
    funcPtJpsi -> Draw("SAME");
    funcPtPsi2s -> Draw("SAME");

    TLegend *legendOniaPt = new TLegend(0.60, 0.70, 0.80, 0.93, " ", "brNDC");
    SetLegend(legendOniaPt, 0.045);
    //legendOniaPt -> AddEntry(histXsecJpsiPt, "J/#psi, all", "F");
    //legendOniaPt -> AddEntry(histXsecPsi2sPt, "#psi(2S), all", "F");
    legendOniaPt -> AddEntry(histXsecJpsiPtFwdCuts, "J/#psi, 2.5 < y < 4", "F");
    legendOniaPt -> AddEntry(histXsecPsi2sPtFwdCuts, "#psi(2S), 2.5 < y < 4", "F");
    legendOniaPt -> AddEntry(graStatJpsiXsecPtRap254, "Data", "FP");
    legendOniaPt -> Draw();*/


    TH1D *histRatioFuncPt = (TH1D*) histFuncPtPsi2s -> Clone("histRatioFuncPt");
    histRatioFuncPt -> Divide(histFuncPtJpsi);
    histRatioFuncPt -> SetLineColor(kOrange+7);
    histRatioFuncPt -> SetLineStyle(kDashed);

    TCanvas *canvasOniaRatioPt = new TCanvas("canvasOniaRatioPt", "", 800, 600);
    histPsi2sOverJpsiPtFwdCuts -> GetYaxis() -> SetRangeUser(0, 0.7);
    histPsi2sOverJpsiPtFwdCuts -> Draw("E2");
    histPromptPsi2sOverPromptJpsiPtFwdCuts -> Draw("E2 SAME");
    histNonPromptPsi2sOverNonPromptJpsiPtFwdCuts -> Draw("E2 SAME");
    graStatPsi2sOverJpsiPtRap254 -> Draw("EP SAME");
    graSystPsi2sOverJpsiPtRap254 -> Draw("E2P SAME");
    //histRatioFuncPt -> Draw("L SAME");

    TLegend *legendOniaRatioPt = new TLegend(0.20, 0.63, 0.40, 0.93, " ", "brNDC");
    SetLegend(legendOniaRatioPt, 0.045);
    legendOniaRatioPt -> AddEntry(histPsi2sOverJpsiPtFwdCuts, "Inclusive", "F");
    legendOniaRatioPt -> AddEntry(histPromptPsi2sOverPromptJpsiPtFwdCuts, "Prompt", "F");
    legendOniaRatioPt -> AddEntry(histNonPromptPsi2sOverNonPromptJpsiPtFwdCuts, "Non Prompt", "F");
    legendOniaRatioPt -> AddEntry(graStatPsi2sOverJpsiPtRap254, "Data", "FP");
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

    TFile *fOut = new TFile(Form("%s/%s_results.root", simPath.c_str(), simName.c_str()), "RECREATE");
    histXsecJpsiPt -> Write();
    histXsecPsi2sPt -> Write();
    histXsecJpsiPtFwdCuts -> Write();
    histXsecPsi2sPtFwdCuts -> Write();
    histPsi2sOverJpsiPt -> Write();
    histPsi2sOverJpsiPtFwdCuts -> Write();
    histXsecJpsiRap -> Write();
    histXsecPsi2sRap -> Write();
    histXsecJpsiRapFwdCuts -> Write();
    histXsecPsi2sRapFwdCuts -> Write();
    histPsi2sOverJpsiRap -> Write();
    histPsi2sOverJpsiRapFwdCuts -> Write();
    fOut -> Close();

    canvasJpsiOniaPt -> SaveAs("figures/xSecJpsiVsPt.pdf");
    canvasPsi2sOniaPt -> SaveAs("figures/xSecPsi2sVsPt.pdf");
    canvasOniaRatioPt -> SaveAs("figures/ratioVsDataVsPt.pdf");
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