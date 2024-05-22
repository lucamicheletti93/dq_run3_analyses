void LoadStyle();
void SetLegend(TLegend *);

void accxeff() {
    LoadStyle();
    //-----------------------------------------------------------//
    // Run 2 Axe
    double ptBinsJpsiRun2[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 26.0, 30.0};
    double axeJpsiPtRun2[] = {0.14923, 0.13782, 0.12011, 0.12440, 0.15167, 0.19440, 0.24144, 0.28442, 0.32242, 0.37021, 0.41486, 0.45900, 0.46818, 0.50202, 0.49940, 0.52091, 0.51781, 0.44586};
    double errAxeJpsiPtRun2[] = {0.00089, 0.00054, 0.00031, 0.00034, 0.00047, 0.00068, 0.00099, 0.00140, 0.00193, 0.00208, 0.00349, 0.00549, 0.00804, 0.01212, 0.01678, 0.02332, 0.02437, 0.03951};

    TH1D *histAxeJpsiPtRun2 = new TH1D("histAxeJpsiPtRun2", "", 18, ptBinsJpsiRun2);
    histAxeJpsiPtRun2 -> SetLineColor(kBlack);
    histAxeJpsiPtRun2 -> SetMarkerStyle(20);
    histAxeJpsiPtRun2 -> SetMarkerColor(kBlack);
    for (int iPt = 0;iPt < 18;iPt++) {
        histAxeJpsiPtRun2 -> SetBinContent(iPt+1, axeJpsiPtRun2[iPt]);
        histAxeJpsiPtRun2 -> SetBinError(iPt+1, errAxeJpsiPtRun2[iPt]);
    }

    double ptBinsPsi2SRun2[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0};
    double axePsi2SPtRun2[] = {0.19816, 0.17608, 0.16092, 0.16698, 0.18982, 0.22608, 0.26571, 0.30449, 0.34395, 0.39066, 0.43539, 0.45861, 0.47282};
    double errAxePsi2SPtRun2[] = {0.00055, 0.00038, 0.00040, 0.00051, 0.00071, 0.00101, 0.00144, 0.00200, 0.00215, 0.00362, 0.00570, 0.00854, 0.01009};

    TH1D *histAxePsi2SPtRun2 = new TH1D("histAxePsi2SPtRun2", "", 13, ptBinsPsi2SRun2);
    histAxePsi2SPtRun2 -> SetLineColor(kBlack);
    histAxePsi2SPtRun2 -> SetMarkerStyle(20);
    histAxePsi2SPtRun2 -> SetMarkerColor(kBlack);
    for (int iPt = 0;iPt < 13;iPt++) {
        histAxePsi2SPtRun2 -> SetBinContent(iPt+1, axePsi2SPtRun2[iPt]);
        histAxePsi2SPtRun2 -> SetBinError(iPt+1, errAxePsi2SPtRun2[iPt]);
    }

    double rapBinsJpsiRun2[] = {2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00};
    double axeJpsiRapRun2[] = {0.04395, 0.15815, 0.23390, 0.24372, 0.18706, 0.07195};
    double errAxeJpsiRapRun2[] = {0.00023, 0.00042, 0.00051, 0.00054, 0.00052, 0.00036};

    TH1D *histAxeJpsiRapRun2 = new TH1D("histAxeJpsiRapRun2", "", 6, rapBinsJpsiRun2);
    histAxeJpsiRapRun2 -> SetLineColor(kBlack);
    histAxeJpsiRapRun2 -> SetMarkerStyle(20);
    histAxeJpsiRapRun2 -> SetMarkerColor(kBlack);
    for (int iY = 0;iY < 6;iY++) {
        histAxeJpsiRapRun2 -> SetBinContent(iY+1, axeJpsiRapRun2[iY]);
        histAxeJpsiRapRun2 -> SetBinError(iY+1, errAxeJpsiRapRun2[iY]);
    }

    double rapBinsPsi2SRun2[] = {2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00};
    double axePsi2SRapRun2[] = {0.05670, 0.19625, 0.28420, 0.29129, 0.21998, 0.08255};
    double errAxePsi2SRapRun2[] = {0.00027, 0.00048, 0.00057, 0.00060, 0.00058, 0.00040};

    TH1D *histAxePsi2SRapRun2 = new TH1D("histAxePsi2SRapRun2", "", 6, rapBinsPsi2SRun2);
    histAxePsi2SRapRun2 -> SetLineColor(kBlack);
    histAxePsi2SRapRun2 -> SetMarkerStyle(20);
    histAxePsi2SRapRun2 -> SetMarkerColor(kBlack);
    for (int iY = 0;iY < 6;iY++) {
        histAxePsi2SRapRun2 -> SetBinContent(iY+1, axePsi2SRapRun2[iY]);
        histAxePsi2SRapRun2 -> SetBinError(iY+1, errAxePsi2SRapRun2[iY]);
    }

    //-----------------------------------------------------------//
    //TFile *fIn = new TFile("/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/MC/central_production/AnalysisResults.root", "READ");
    //string productionName = "LHC24e5"; // prompt charmonia at forward https://its.cern.ch/jira/browse/O2-4884
    string productionName = "LHC24e4"; // non-prompt charmonia at forward https://its.cern.ch/jira/browse/O2-4885
    TFile *fIn = new TFile(Form("/Users/lucamicheletti/cernbox/JPSI/Run3/MC/%s/AnalysisResults_dq_efficiency.root", productionName.c_str()), "READ");

    TList *listGen1 = (TList*) fIn -> Get("analysis-same-event-pairing/output");
    TList *listGenJpsi2 = (TList*) listGen1 -> FindObject("MCTruthGen_Jpsi");
    TH2D *histGenJpsiPtRap = (TH2D*) listGenJpsi2 -> FindObject("Pt_Rapidity");
    histGenJpsiPtRap -> GetYaxis() -> SetRangeUser(2.5, 4); // WARNING! constrain the generated particles to be in the 2.5 < y < 4 region

    TH1D *histGenJpsiPt = (TH1D*) histGenJpsiPtRap -> ProjectionX("histGenJpsiPt");
    TH1D *histGenJpsiRap = (TH1D*) histGenJpsiPtRap -> ProjectionY("histGenJpsiRap");

    histGenJpsiPt -> SetLineColor(kBlack);
    histGenJpsiRap -> SetLineColor(kBlack);
    histGenJpsiPt -> SetMarkerColor(kBlack);
    histGenJpsiRap -> SetMarkerColor(kBlack);

    TList *listPsi2SGen2 = (TList*) listGen1 -> FindObject("MCTruthGen_Psi2S");
    TH2D *histPsi2SPtYGen = (TH2D*) listPsi2SGen2 -> FindObject("Pt_Rapidity");
    histPsi2SPtYGen -> GetYaxis() -> SetRangeUser(2.5, 4); // WARNING! constrain the generated particles to be in the 2.5 < y < 4 region
    TH1D *histGenPsi2SPt = (TH1D*) histPsi2SPtYGen -> ProjectionX("histGenPsi2SPt");
    TH1D *histGenPsi2SRap = (TH1D*) histPsi2SPtYGen -> ProjectionY("histGenPsi2SRap");

    histGenPsi2SPt -> SetLineColor(kBlack);
    histGenPsi2SRap -> SetLineColor(kBlack);
    histGenPsi2SPt -> SetMarkerColor(kBlack);
    histGenPsi2SRap -> SetMarkerColor(kBlack);

    TList *listRec1 = (TList*) fIn -> Get("analysis-same-event-pairing/output");

    TList *listRecJpsiCut1 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_matchedMchMid_mumuFromJpsi");
    TList *listRecJpsiCut2 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt10SigmaPDCA_mumuFromJpsi");
    TList *listRecJpsiCut3 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA_mumuFromJpsi");
    TList *listRecJpsiCut4 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt510SigmaPDCA_mumuFromJpsi");

    THnSparseD *histRecJpsiMassPtRapCut1 = (THnSparseD*) listRecJpsiCut1 -> FindObject("Mass_Pt_Rapidity");
    histRecJpsiMassPtRapCut1 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecJpsiPtRapCut1 = (TH2D*) histRecJpsiMassPtRapCut1 -> Projection(2, 1, "histRecJpsiPtRapCut1");
    TH1D *histRecJpsiPtCut1 = (TH1D*) histRecJpsiMassPtRapCut1 -> Projection(1, "histRecJpsiPtCut1");
    TH1D *histRecJpsiRapCut1 = (TH1D*) histRecJpsiMassPtRapCut1 -> Projection(2, "histRecJpsiRapCut1");

    THnSparseD *histRecJpsiMassPtRapCut2 = (THnSparseD*) listRecJpsiCut2 -> FindObject("Mass_Pt_Rapidity");
    histRecJpsiMassPtRapCut2 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecJpsiPtRapCut2 = (TH2D*) histRecJpsiMassPtRapCut2 -> Projection(2, 1, "histRecJpsiPtRapCut2");
    TH1D *histRecJpsiPtCut2 = (TH1D*) histRecJpsiMassPtRapCut2 -> Projection(1, "histRecJpsiPtCut2");
    TH1D *histRecJpsiRapCut2 = (TH1D*) histRecJpsiMassPtRapCut2 -> Projection(2, "histRecJpsiRapCut2");

    THnSparseD *histRecJpsiMassPtRapCut3 = (THnSparseD*) listRecJpsiCut3 -> FindObject("Mass_Pt_Rapidity");
    histRecJpsiMassPtRapCut3 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecJpsiPtRapCut3 = (TH2D*) histRecJpsiMassPtRapCut3 -> Projection(2, 1, "histRecJpsiPtRapCut3");
    TH1D *histRecJpsiPtCut3 = (TH1D*) histRecJpsiMassPtRapCut3 -> Projection(1, "histRecJpsiPtCut3");
    TH1D *histRecJpsiRapCut3 = (TH1D*) histRecJpsiMassPtRapCut3 -> Projection(2, "histRecJpsiRapCut3");

    THnSparseD *histRecJpsiMassPtRapCut4 = (THnSparseD*) listRecJpsiCut4 -> FindObject("Mass_Pt_Rapidity");
    histRecJpsiMassPtRapCut4 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecJpsiPtRapCut4 = (TH2D*) histRecJpsiMassPtRapCut4 -> Projection(2, 1, "histRecJpsiPtRapCut4");
    TH1D *histRecJpsiPtCut4 = (TH1D*) histRecJpsiMassPtRapCut4 -> Projection(1, "histRecJpsiPtCut4");
    TH1D *histRecJpsiRapCut4 = (TH1D*) histRecJpsiMassPtRapCut4 -> Projection(2, "histRecJpsiRapCut4");

    TList *listRecPsi2SCut1 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_matchedMchMid_mumuFromPsi2S");
    TList *listRecPsi2SCut2 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt10SigmaPDCA_mumuFromPsi2S");
    TList *listRecPsi2SCut3 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA_mumuFromPsi2S");
    TList *listRecPsi2SCut4 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt510SigmaPDCA_mumuFromPsi2S");

    THnSparseD *histRecPsi2SMassPtRapCut1 = (THnSparseD*) listRecPsi2SCut1 -> FindObject("Mass_Pt_Rapidity");
    histRecPsi2SMassPtRapCut1 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecPsi2SPtRapCut1 = (TH2D*) histRecPsi2SMassPtRapCut1 -> Projection(2, 1, "histRecPsi2SPtRapCut1");
    TH1D *histRecPsi2SPtCut1 = (TH1D*) histRecPsi2SMassPtRapCut1 -> Projection(1, "histRecPsi2SPtCut1");
    TH1D *histRecPsi2SRapCut1 = (TH1D*) histRecPsi2SMassPtRapCut1 -> Projection(2, "histRecPsi2SRapCut1");

    THnSparseD *histRecPsi2SMassPtRapCut2 = (THnSparseD*) listRecPsi2SCut2 -> FindObject("Mass_Pt_Rapidity");
    histRecPsi2SMassPtRapCut2 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecPsi2SPtRapCut2 = (TH2D*) histRecPsi2SMassPtRapCut2 -> Projection(2, 1, "histRecPsi2SPtRapCut2");
    TH1D *histRecPsi2SPtCut2 = (TH1D*) histRecPsi2SMassPtRapCut2 -> Projection(1, "histRecPsi2SPtCut2");
    TH1D *histRecPsi2SRapCut2 = (TH1D*) histRecPsi2SMassPtRapCut2 -> Projection(2, "histRecPsi2SRapCut2");

    THnSparseD *histRecPsi2SMassPtRapCut3 = (THnSparseD*) listRecPsi2SCut3 -> FindObject("Mass_Pt_Rapidity");
    histRecPsi2SMassPtRapCut3 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecPsi2SPtRapCut3 = (TH2D*) histRecPsi2SMassPtRapCut3 -> Projection(2, 1, "histRecPsi2SPtRapCut3");
    TH1D *histRecPsi2SPtCut3 = (TH1D*) histRecPsi2SMassPtRapCut3 -> Projection(1, "histRecPsi2SPtCut3");
    TH1D *histRecPsi2SRapCut3 = (TH1D*) histRecPsi2SMassPtRapCut3 -> Projection(2, "histRecPsi2SRapCut3");

    THnSparseD *histRecPsi2SMassPtRapCut4 = (THnSparseD*) listRecPsi2SCut4 -> FindObject("Mass_Pt_Rapidity");
    histRecPsi2SMassPtRapCut4 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecPsi2SPtRapCut4 = (TH2D*) histRecPsi2SMassPtRapCut4 -> Projection(2, 1, "histRecPsi2SPtRapCut4");
    TH1D *histRecPsi2SPtCut4 = (TH1D*) histRecPsi2SMassPtRapCut4 -> Projection(1, "histRecPsi2SPtCut4");
    TH1D *histRecPsi2SRapCut4 = (TH1D*) histRecPsi2SMassPtRapCut4 -> Projection(2, "histRecPsi2SRapCut4");

    histGenJpsiPtRap -> RebinX(2);
    histRecJpsiPtRapCut1 -> RebinX(2);
    histRecJpsiPtRapCut2 -> RebinX(2);
    histRecJpsiPtRapCut3 -> RebinX(2);
    histRecJpsiPtRapCut4 -> RebinX(2);

    histGenJpsiPtRap -> RebinY(5);
    histRecJpsiPtRapCut1 -> RebinY(2);
    histRecJpsiPtRapCut2 -> RebinY(2);
    histRecJpsiPtRapCut3 -> RebinY(2);
    histRecJpsiPtRapCut4 -> RebinY(2);

    TH2D *histAxeJpsiPtRapCut1 = (TH2D*) histRecJpsiPtRapCut1 -> Clone("histAxeJpsiPtRapCut1");
    histAxeJpsiPtRapCut1 -> Divide(histGenJpsiPtRap);

    TH2D *histAxeJpsiPtRapCut2 = (TH2D*) histRecJpsiPtRapCut2 -> Clone("histAxeJpsiPtRapCut2");
    histAxeJpsiPtRapCut2 -> Divide(histGenJpsiPtRap);

    TH2D *histAxeJpsiPtRapCut3 = (TH2D*) histRecJpsiPtRapCut3 -> Clone("histAxeJpsiPtRapCut3");
    histAxeJpsiPtRapCut3 -> Divide(histGenJpsiPtRap);

    TH2D *histAxeJpsiPtRapCut4 = (TH2D*) histRecJpsiPtRapCut4 -> Clone("histAxeJpsiPtRapCut4");
    histAxeJpsiPtRapCut4 -> Divide(histGenJpsiPtRap);

    TCanvas *canvasJpsi_2D = new TCanvas("canvasJpsi_2D", "", 1200, 1200);
    canvasJpsi_2D -> Divide(2, 2);
    canvasJpsi_2D -> cd(1);
    histAxeJpsiPtRapCut1 -> GetZaxis() -> SetRangeUser(0, 1);
    histAxeJpsiPtRapCut1 -> Draw("COLZ");
    canvasJpsi_2D -> cd(2);
    histAxeJpsiPtRapCut2 -> GetZaxis() -> SetRangeUser(0, 1);
    histAxeJpsiPtRapCut2 -> Draw("COLZ");
    canvasJpsi_2D -> cd(3);
    histAxeJpsiPtRapCut3 -> GetZaxis() -> SetRangeUser(0, 1);
    histAxeJpsiPtRapCut3 -> Draw("COLZ");
    canvasJpsi_2D -> cd(4);
    histAxeJpsiPtRapCut4 -> GetZaxis() -> SetRangeUser(0, 1);
    histAxeJpsiPtRapCut4 -> Draw("COLZ");

    // Rebin the histograms
    TH1D *histGenJpsiPtRebin = (TH1D*) histGenJpsiPt -> Rebin(18, "histGenJpsiPtRebin", ptBinsJpsiRun2); 
    TH1D *histRecJpsiPtCut1Rebin = (TH1D*) histRecJpsiPtCut1 -> Rebin(18, "histRecJpsiPtCut1Rebin", ptBinsJpsiRun2); 
    TH1D *histRecJpsiPtCut2Rebin = (TH1D*) histRecJpsiPtCut2 -> Rebin(18, "histRecJpsiPtCut2Rebin", ptBinsJpsiRun2); 
    TH1D *histRecJpsiPtCut3Rebin = (TH1D*) histRecJpsiPtCut3 -> Rebin(18, "histRecJpsiPtCut3Rebin", ptBinsJpsiRun2); 
    TH1D *histRecJpsiPtCut4Rebin = (TH1D*) histRecJpsiPtCut4 -> Rebin(18, "histRecJpsiPtCut4Rebin", ptBinsJpsiRun2); 

    histRecJpsiPtCut1Rebin -> SetLineColor(kRed+1);
    histRecJpsiPtCut2Rebin -> SetLineColor(kOrange+7);
    histRecJpsiPtCut3Rebin -> SetLineColor(kAzure+2);
    histRecJpsiPtCut4Rebin -> SetLineColor(kBlue+1);

    histRecJpsiPtCut1Rebin -> SetLineWidth(2);
    histRecJpsiPtCut2Rebin -> SetLineWidth(2);
    histRecJpsiPtCut3Rebin -> SetLineWidth(2);
    histRecJpsiPtCut4Rebin -> SetLineWidth(2);

    histRecJpsiPtCut1Rebin -> SetMarkerColor(kRed+1);
    histRecJpsiPtCut2Rebin -> SetMarkerColor(kOrange+7);
    histRecJpsiPtCut3Rebin -> SetMarkerColor(kAzure+2);
    histRecJpsiPtCut4Rebin -> SetMarkerColor(kBlue+1);

    TH1D *histAxeJpsiPtCut1 = new TH1D("histAxeJpsiPtCut1", "", 18, ptBinsJpsiRun2);
    histAxeJpsiPtCut1 -> Divide(histRecJpsiPtCut1Rebin, histGenJpsiPtRebin, 1, 1, "B");
    histAxeJpsiPtCut1 -> SetLineColor(kRed+1);
    histAxeJpsiPtCut1 -> SetMarkerColor(kRed+1);

    TH1D *histAxeJpsiPtCut2 = new TH1D("histAxeJpsiPtCut2", "", 18, ptBinsJpsiRun2);
    histAxeJpsiPtCut2 -> Divide(histRecJpsiPtCut2Rebin, histGenJpsiPtRebin, 1, 1, "B");
    histAxeJpsiPtCut2 -> SetLineColor(kOrange+7);
    histAxeJpsiPtCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histAxeJpsiPtCut3 = new TH1D("histAxeJpsiPtCut3", "", 18, ptBinsJpsiRun2);
    histAxeJpsiPtCut3 -> Divide(histRecJpsiPtCut3Rebin, histGenJpsiPtRebin, 1, 1, "B");
    histAxeJpsiPtCut3 -> SetLineColor(kAzure+2);
    histAxeJpsiPtCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histAxeJpsiPtCut4 = new TH1D("histAxeJpsiPtCut4", "", 18, ptBinsJpsiRun2);
    histAxeJpsiPtCut4 -> Divide(histRecJpsiPtCut4Rebin, histGenJpsiPtRebin, 1, 1, "B");
    histAxeJpsiPtCut4 -> SetLineColor(kBlue+1);
    histAxeJpsiPtCut4 -> SetMarkerColor(kBlue+1);

    TH1D *histGenJpsiRapRebin = (TH1D*) histGenJpsiRap -> Rebin(6, "histGenJpsiRapRebin", rapBinsJpsiRun2); 
    TH1D *histRecJpsiRapCut1Rebin = (TH1D*) histRecJpsiRapCut1 -> Rebin(6, "histRecJpsiRapCut1Rebin", rapBinsJpsiRun2); 
    TH1D *histRecJpsiRapCut2Rebin = (TH1D*) histRecJpsiRapCut2 -> Rebin(6, "histRecJpsiRapCut2Rebin", rapBinsJpsiRun2); 
    TH1D *histRecJpsiRapCut3Rebin = (TH1D*) histRecJpsiRapCut3 -> Rebin(6, "histRecJpsiRapCut3Rebin", rapBinsJpsiRun2); 
    TH1D *histRecJpsiRapCut4Rebin = (TH1D*) histRecJpsiRapCut4 -> Rebin(6, "histRecJpsiRapCut4Rebin", rapBinsJpsiRun2); 

    histRecJpsiRapCut1Rebin -> SetLineColor(kRed+1);
    histRecJpsiRapCut2Rebin -> SetLineColor(kOrange+7);
    histRecJpsiRapCut3Rebin -> SetLineColor(kAzure+2);
    histRecJpsiRapCut4Rebin -> SetLineColor(kBlue+1);

    histRecJpsiRapCut1Rebin -> SetLineWidth(2);
    histRecJpsiRapCut2Rebin -> SetLineWidth(2);
    histRecJpsiRapCut3Rebin -> SetLineWidth(2);
    histRecJpsiRapCut4Rebin -> SetLineWidth(2);

    histRecJpsiRapCut1Rebin -> SetMarkerColor(kRed+1);
    histRecJpsiRapCut2Rebin -> SetMarkerColor(kOrange+7);
    histRecJpsiRapCut3Rebin -> SetMarkerColor(kAzure+2);
    histRecJpsiRapCut4Rebin -> SetMarkerColor(kBlue+1);

    TH1D *histAxeJpsiRapCut1 = new TH1D("histAxeJpsiRapCut1", "", 6, rapBinsJpsiRun2);
    histAxeJpsiRapCut1 -> Divide(histRecJpsiRapCut1Rebin, histGenJpsiRapRebin, 1, 1, "B");
    histAxeJpsiRapCut1 -> SetLineColor(kRed+1);
    histAxeJpsiRapCut1 -> SetMarkerColor(kRed+1);

    TH1D *histAxeJpsiRapCut2 = new TH1D("histAxeJpsiRapCut2", "", 6, rapBinsJpsiRun2);
    histAxeJpsiRapCut2 -> Divide(histRecJpsiRapCut2Rebin, histGenJpsiRapRebin, 1, 1, "B");
    histAxeJpsiRapCut2 -> SetLineColor(kOrange+7);
    histAxeJpsiRapCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histAxeJpsiRapCut3 = new TH1D("histAxeJpsiRapCut3", "", 6, rapBinsJpsiRun2);
    histAxeJpsiRapCut3 -> Divide(histRecJpsiRapCut3Rebin, histGenJpsiRapRebin, 1, 1, "B");
    histAxeJpsiRapCut3 -> SetLineColor(kAzure+2);
    histAxeJpsiRapCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histAxeJpsiRapCut4 = new TH1D("histAxeJpsiRapCut4", "", 6, rapBinsJpsiRun2);
    histAxeJpsiRapCut4 -> Divide(histRecJpsiRapCut4Rebin, histGenJpsiRapRebin, 1, 1, "B");
    histAxeJpsiRapCut4 -> SetLineColor(kBlue+1);
    histAxeJpsiRapCut4 -> SetMarkerColor(kBlue+1);

    TH1D *histGenPsi2SPtRebin = (TH1D*) histGenPsi2SPt -> Rebin(18, "histGenPsi2SPtRebin", ptBinsJpsiRun2); 
    TH1D *histRecPsi2SPtCut1Rebin = (TH1D*) histRecPsi2SPtCut1 -> Rebin(18, "histRecPsi2SPtCut1Rebin", ptBinsJpsiRun2); 
    TH1D *histRecPsi2SPtCut2Rebin = (TH1D*) histRecPsi2SPtCut2 -> Rebin(18, "histRecPsi2SPtCut2Rebin", ptBinsJpsiRun2); 
    TH1D *histRecPsi2SPtCut3Rebin = (TH1D*) histRecPsi2SPtCut3 -> Rebin(18, "histRecPsi2SPtCut3Rebin", ptBinsJpsiRun2); 
    TH1D *histRecPsi2SPtCut4Rebin = (TH1D*) histRecPsi2SPtCut4 -> Rebin(18, "histRecPsi2SPtCut4Rebin", ptBinsJpsiRun2); 

    histRecPsi2SPtCut1Rebin -> SetLineColor(kRed+1);
    histRecPsi2SPtCut2Rebin -> SetLineColor(kOrange+7);
    histRecPsi2SPtCut3Rebin -> SetLineColor(kAzure+2);
    histRecPsi2SPtCut4Rebin -> SetLineColor(kBlue+1);

    histRecPsi2SPtCut1Rebin -> SetLineWidth(2);
    histRecPsi2SPtCut2Rebin -> SetLineWidth(2);
    histRecPsi2SPtCut3Rebin -> SetLineWidth(2);
    histRecPsi2SPtCut4Rebin -> SetLineWidth(2);

    histRecPsi2SPtCut1Rebin -> SetMarkerColor(kRed+1);
    histRecPsi2SPtCut2Rebin -> SetMarkerColor(kOrange+7);
    histRecPsi2SPtCut3Rebin -> SetMarkerColor(kAzure+2);
    histRecPsi2SPtCut4Rebin -> SetMarkerColor(kBlue+1);

    TH1D *histAxePsi2SPtCut1 = new TH1D("histAxePsi2SPtCut1", "", 18, ptBinsJpsiRun2);
    histAxePsi2SPtCut1 -> Divide(histRecPsi2SPtCut1Rebin, histGenPsi2SPtRebin, 1, 1, "B");
    histAxePsi2SPtCut1 -> SetLineColor(kRed+1);
    histAxePsi2SPtCut1 -> SetMarkerColor(kRed+1);

    TH1D *histAxePsi2SPtCut2 = new TH1D("histAxePsi2SPtCut2", "", 18, ptBinsJpsiRun2);
    histAxePsi2SPtCut2 -> Divide(histRecPsi2SPtCut2Rebin, histGenPsi2SPtRebin, 1, 1, "B");
    histAxePsi2SPtCut2 -> SetLineColor(kOrange+7);
    histAxePsi2SPtCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histAxePsi2SPtCut3 = new TH1D("histAxePsi2SPtCut3", "", 18, ptBinsJpsiRun2);
    histAxePsi2SPtCut3 -> Divide(histRecPsi2SPtCut3Rebin, histGenPsi2SPtRebin, 1, 1, "B");
    histAxePsi2SPtCut3 -> SetLineColor(kAzure+2);
    histAxePsi2SPtCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histAxePsi2SPtCut4 = new TH1D("histAxePsi2SPtCut4", "", 18, ptBinsJpsiRun2);
    histAxePsi2SPtCut4 -> Divide(histRecPsi2SPtCut4Rebin, histGenPsi2SPtRebin, 1, 1, "B");
    histAxePsi2SPtCut4 -> SetLineColor(kBlue+1);
    histAxePsi2SPtCut4 -> SetMarkerColor(kBlue+1);

    TH1D *histGenPsi2SRapRebin = (TH1D*) histGenPsi2SRap -> Rebin(6, "histGenPsi2SRapRebin", rapBinsJpsiRun2); 
    TH1D *histRecPsi2SRapCut1Rebin = (TH1D*) histRecPsi2SRapCut1 -> Rebin(6, "histRecPsi2SRapCut1Rebin", rapBinsJpsiRun2); 
    TH1D *histRecPsi2SRapCut2Rebin = (TH1D*) histRecPsi2SRapCut2 -> Rebin(6, "histRecPsi2SRapCut2Rebin", rapBinsJpsiRun2); 
    TH1D *histRecPsi2SRapCut3Rebin = (TH1D*) histRecPsi2SRapCut3 -> Rebin(6, "histRecPsi2SRapCut3Rebin", rapBinsJpsiRun2); 
    TH1D *histRecPsi2SRapCut4Rebin = (TH1D*) histRecPsi2SRapCut4 -> Rebin(6, "histRecPsi2SRapCut4Rebin", rapBinsJpsiRun2); 

    histRecPsi2SRapCut1Rebin -> SetLineColor(kRed+1);
    histRecPsi2SRapCut2Rebin -> SetLineColor(kOrange+7);
    histRecPsi2SRapCut3Rebin -> SetLineColor(kAzure+2);
    histRecPsi2SRapCut4Rebin -> SetLineColor(kBlue+1);

    histRecPsi2SRapCut1Rebin -> SetLineWidth(2);
    histRecPsi2SRapCut2Rebin -> SetLineWidth(2);
    histRecPsi2SRapCut3Rebin -> SetLineWidth(2);
    histRecPsi2SRapCut4Rebin -> SetLineWidth(2);

    histRecPsi2SRapCut1Rebin -> SetMarkerColor(kRed+1);
    histRecPsi2SRapCut2Rebin -> SetMarkerColor(kOrange+7);
    histRecPsi2SRapCut3Rebin -> SetMarkerColor(kAzure+2);
    histRecPsi2SRapCut4Rebin -> SetMarkerColor(kBlue+1);

    TH1D *histAxePsi2SRapCut1 = new TH1D("histAxePsi2SRapCut1", "", 6, rapBinsPsi2SRun2);
    histAxePsi2SRapCut1 -> Divide(histRecPsi2SRapCut1Rebin, histGenPsi2SRapRebin, 1, 1, "B");
    histAxePsi2SRapCut1 -> SetLineColor(kRed+1);
    histAxePsi2SRapCut1 -> SetMarkerColor(kRed+1);

    TH1D *histAxePsi2SRapCut2 = new TH1D("histAxePsi2SRapCut2", "", 6, rapBinsPsi2SRun2);
    histAxePsi2SRapCut2 -> Divide(histRecPsi2SRapCut2Rebin, histGenPsi2SRapRebin, 1, 1, "B");
    histAxePsi2SRapCut2 -> SetLineColor(kOrange+7);
    histAxePsi2SRapCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histAxePsi2SRapCut3 = new TH1D("histAxePsi2SRapCut3", "", 6, rapBinsPsi2SRun2);
    histAxePsi2SRapCut3 -> Divide(histRecPsi2SRapCut3Rebin, histGenPsi2SRapRebin, 1, 1, "B");
    histAxePsi2SRapCut3 -> SetLineColor(kAzure+2);
    histAxePsi2SRapCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histAxePsi2SRapCut4 = new TH1D("histAxePsi2SRapCut4", "", 6, rapBinsPsi2SRun2);
    histAxePsi2SRapCut4 -> Divide(histRecPsi2SRapCut4Rebin, histGenPsi2SRapRebin, 1, 1, "B");
    histAxePsi2SRapCut4 -> SetLineColor(kBlue+1);
    histAxePsi2SRapCut4 -> SetMarkerColor(kBlue+1);


    TCanvas *canvasJpsi = new TCanvas("canvasJpsi", "", 1200, 1200);
    canvasJpsi -> Divide(2, 2);

    canvasJpsi -> cd(1);
    gPad -> SetLogy(1);
    histGenJpsiPtRebin -> SetTitle("");
    histGenJpsiPtRebin -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGenJpsiPtRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histGenJpsiPtRebin -> GetYaxis() -> SetRangeUser(1, 1e6);
    histGenJpsiPtRebin -> Draw("EP");
    histRecJpsiPtCut1Rebin -> Draw("EP SAME");
    histRecJpsiPtCut2Rebin -> Draw("EP SAME");
    histRecJpsiPtCut3Rebin -> Draw("EP SAME");
    histRecJpsiPtCut4Rebin -> Draw("EP SAME");

    canvasJpsi -> cd(2);
    gPad -> SetLogy(1);
    histGenJpsiRapRebin -> SetTitle("");
    histGenJpsiRapRebin -> GetXaxis() -> SetTitle("#it{y}");
    histGenJpsiRapRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histGenJpsiRapRebin -> GetYaxis() -> SetRangeUser(1e3, 1e6);
    histGenJpsiRapRebin -> Draw("EP");
    histRecJpsiRapCut1Rebin -> Draw("EP SAME");
    histRecJpsiRapCut2Rebin -> Draw("EP SAME");
    histRecJpsiRapCut3Rebin -> Draw("EP SAME");
    histRecJpsiRapCut4Rebin -> Draw("EP SAME");

    canvasJpsi -> cd(3);
    histAxeJpsiPtCut1 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histAxeJpsiPtCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histAxeJpsiPtCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histAxeJpsiPtCut1 -> Draw("EP");
    histAxeJpsiPtCut2 -> Draw("EP SAME");
    histAxeJpsiPtCut3 -> Draw("EP SAME");
    histAxeJpsiPtCut4 -> Draw("EP SAME");
    histAxeJpsiPtRun2 -> Draw("EP SAME");

    TLegend *legendJpsiPtAxe = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendJpsiPtAxe);
    legendJpsiPtAxe -> AddEntry(histAxeJpsiPtRun2, "Run 2", "PL");
    legendJpsiPtAxe -> AddEntry(histAxeJpsiPtCut1, "Run 3 - matchedMchMid", "PL");
    legendJpsiPtAxe -> AddEntry(histAxeJpsiPtCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendJpsiPtAxe -> AddEntry(histAxeJpsiPtCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendJpsiPtAxe -> AddEntry(histAxeJpsiPtCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendJpsiPtAxe -> Draw("SAME");

    canvasJpsi -> cd(4);
    histAxeJpsiRapCut1 -> GetXaxis() -> SetTitle("#it{y}");
    histAxeJpsiRapCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histAxeJpsiRapCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histAxeJpsiRapCut1 -> Draw("EP");
    histAxeJpsiRapCut2 -> Draw("EP SAME");
    histAxeJpsiRapCut3 -> Draw("EP SAME");
    histAxeJpsiRapCut4 -> Draw("EP SAME");
    histAxeJpsiRapRun2 -> Draw("EP SAME");

    TLegend *legendJpsiYAxe = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendJpsiYAxe);
    legendJpsiYAxe -> AddEntry(histAxeJpsiRapRun2, "Run 2", "PL");
    legendJpsiYAxe -> AddEntry(histAxeJpsiRapCut1, "Run 3 - matchedMchMid", "PL");
    legendJpsiYAxe -> AddEntry(histAxeJpsiRapCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendJpsiYAxe -> AddEntry(histAxeJpsiRapCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendJpsiYAxe -> AddEntry(histAxeJpsiRapCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendJpsiYAxe -> Draw("SAME");


    TCanvas *canvasPsi2S = new TCanvas("canvasPsi2S", "", 1200, 1200);
    canvasPsi2S -> Divide(2, 2);

    canvasPsi2S -> cd(1);
    gPad -> SetLogy(1);
    histGenPsi2SPtRebin -> SetTitle("");
    histGenPsi2SPtRebin -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGenPsi2SPtRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histGenPsi2SPtRebin -> GetYaxis() -> SetRangeUser(1, 1e6);
    histGenPsi2SPtRebin -> Draw("EP");
    histRecPsi2SPtCut1Rebin -> Draw("EP SAME");
    histRecPsi2SPtCut2Rebin -> Draw("EP SAME");
    histRecPsi2SPtCut3Rebin -> Draw("EP SAME");
    histRecPsi2SPtCut4Rebin -> Draw("EP SAME");

    canvasPsi2S -> cd(2);
    gPad -> SetLogy(1);
    histGenPsi2SRapRebin -> SetTitle("");
    histGenPsi2SRapRebin -> GetXaxis() -> SetTitle("#it{y}");
    histGenPsi2SRapRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histGenPsi2SRapRebin -> GetYaxis() -> SetRangeUser(1e3, 1e6);
    histGenPsi2SRapRebin -> Draw("EP");
    histRecPsi2SRapCut1Rebin -> Draw("EP SAME");
    histRecPsi2SRapCut2Rebin -> Draw("EP SAME");
    histRecPsi2SRapCut3Rebin -> Draw("EP SAME");
    histRecPsi2SRapCut4Rebin -> Draw("EP SAME");

    canvasPsi2S -> cd(3);
    histAxePsi2SPtCut1 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histAxePsi2SPtCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histAxePsi2SPtCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histAxePsi2SPtCut1 -> Draw("EP");
    histAxePsi2SPtCut2 -> Draw("EP SAME");
    histAxePsi2SPtCut3 -> Draw("EP SAME");
    histAxePsi2SPtCut4 -> Draw("EP SAME");
    histAxePsi2SPtRun2 -> Draw("EP SAME");
    

    TLegend *legendPsi2SPtAxe = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendPsi2SPtAxe);
    legendPsi2SPtAxe -> AddEntry(histAxePsi2SPtCut1, "Run 3 - matchedMchMid", "PL");
    legendPsi2SPtAxe -> AddEntry(histAxePsi2SPtCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendPsi2SPtAxe -> AddEntry(histAxePsi2SPtCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendPsi2SPtAxe -> AddEntry(histAxePsi2SPtCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendPsi2SPtAxe -> Draw("SAME");

    canvasPsi2S -> cd(4);
    histAxePsi2SRapCut1 -> GetXaxis() -> SetTitle("#it{y}");
    histAxePsi2SRapCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histAxePsi2SRapCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histAxePsi2SRapCut1 -> Draw("EP");
    histAxePsi2SRapCut2 -> Draw("EP SAME");
    histAxePsi2SRapCut3 -> Draw("EP SAME");
    histAxePsi2SRapCut4 -> Draw("EP SAME");
    histAxePsi2SRapRun2 -> Draw("EP SAME");

    TLegend *legendPsi2SYAxe = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendPsi2SYAxe);
    legendPsi2SYAxe -> AddEntry(histAxePsi2SRapCut1, "Run 3 - matchedMchMid", "PL");
    legendPsi2SYAxe -> AddEntry(histAxePsi2SRapCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendPsi2SYAxe -> AddEntry(histAxePsi2SRapCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendPsi2SYAxe -> AddEntry(histAxePsi2SRapCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendPsi2SYAxe -> Draw("SAME");

    // Axe ratio
    TH1D *histRatioPtAxeCut1 = (TH1D*) histAxePsi2SPtCut1 -> Clone("histRatioPtAxeCut1");
    histRatioPtAxeCut1 -> Divide(histAxeJpsiPtCut1);
    histRatioPtAxeCut1 -> SetLineColor(kRed+1);
    histRatioPtAxeCut1 -> SetMarkerColor(kRed+1);

    TH1D *histRatioPtAxeCut2 = (TH1D*) histAxePsi2SPtCut2 -> Clone("histRatioPtAxeCut2");
    histRatioPtAxeCut2 -> Divide(histAxeJpsiPtCut2);
    histRatioPtAxeCut2 -> SetLineColor(kOrange+7);
    histRatioPtAxeCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histRatioPtAxeCut3 = (TH1D*) histAxePsi2SPtCut3 -> Clone("histRatioPtAxeCut3");
    histRatioPtAxeCut3 -> Divide(histAxeJpsiPtCut3);
    histRatioPtAxeCut3 -> SetLineColor(kAzure+2);
    histRatioPtAxeCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histRatioPtAxeCut4 = (TH1D*) histAxePsi2SPtCut4 -> Clone("histRatioPtAxeCut4");
    histRatioPtAxeCut4 -> Divide(histAxeJpsiPtCut4);
    histRatioPtAxeCut4 -> SetLineColor(kBlue+1);
    histRatioPtAxeCut4 -> SetMarkerColor(kBlue+1);


    TCanvas *canvasRatioPt = new TCanvas("canvasRatioPt", "", 800, 600);
    histRatioPtAxeCut1 -> SetTitle("");
    histRatioPtAxeCut1 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histRatioPtAxeCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon_{#psi(2S)} / A#times#varepsilon_{J/#psi}");
    histRatioPtAxeCut1 -> GetYaxis() -> SetRangeUser(0, 2);
    histRatioPtAxeCut1 -> Draw("EP");
    histRatioPtAxeCut2 -> Draw("EP SAME");
    histRatioPtAxeCut3 -> Draw("EP SAME");
    histRatioPtAxeCut4 -> Draw("EP SAME");



    TH1D *histRatioYAxeCut1 = (TH1D*) histAxePsi2SRapCut1 -> Clone("histRatioYAxeCut1");
    histRatioYAxeCut1 -> Divide(histAxeJpsiRapCut1);
    histRatioYAxeCut1 -> SetLineColor(kRed+1);
    histRatioYAxeCut1 -> SetMarkerColor(kRed+1);

    TH1D *histRatioYAxeCut2 = (TH1D*) histAxePsi2SRapCut2 -> Clone("histRatioYAxeCut2");
    histRatioYAxeCut2 -> Divide(histAxeJpsiRapCut2);
    histRatioYAxeCut2 -> SetLineColor(kOrange+7);
    histRatioYAxeCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histRatioYAxeCut3 = (TH1D*) histAxePsi2SRapCut3 -> Clone("histRatioYAxeCut3");
    histRatioYAxeCut3 -> Divide(histAxeJpsiRapCut3);
    histRatioYAxeCut3 -> SetLineColor(kAzure+2);
    histRatioYAxeCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histRatioYAxeCut4 = (TH1D*) histAxePsi2SRapCut4 -> Clone("histRatioYAxeCut4");
    histRatioYAxeCut4 -> Divide(histAxeJpsiRapCut4);
    histRatioYAxeCut4 -> SetLineColor(kBlue+1);
    histRatioYAxeCut4 -> SetMarkerColor(kBlue+1);


    TCanvas *canvasRatioY = new TCanvas("canvasRatioY", "", 800, 600);
    histRatioYAxeCut1 -> SetTitle("");
    histRatioYAxeCut1 -> GetXaxis() -> SetTitle("#it{y}");
    histRatioYAxeCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon_{#psi(2S)} / A#times#varepsilon_{J/#psi}");
    histRatioYAxeCut1 -> GetYaxis() -> SetRangeUser(0, 2);
    histRatioYAxeCut1 -> Draw("EP");
    histRatioYAxeCut2 -> Draw("EP SAME");
    histRatioYAxeCut3 -> Draw("EP SAME");
    histRatioYAxeCut4 -> Draw("EP SAME");


    TFile *fOut = new TFile(Form("Axe_%s.root", productionName.c_str()), "RECREATE");
    histAxeJpsiPtCut1 -> Write();
    histAxeJpsiPtCut2 -> Write();
    histAxeJpsiPtCut3 -> Write();
    histAxeJpsiPtCut4 -> Write();
    histAxePsi2SPtCut1 -> Write();
    histAxePsi2SPtCut2 -> Write();
    histAxePsi2SPtCut3 -> Write();
    histAxePsi2SPtCut4 -> Write();
    histAxeJpsiRapCut1 -> Write();
    histAxeJpsiRapCut2 -> Write();
    histAxeJpsiRapCut3 -> Write();
    histAxeJpsiRapCut4 -> Write();
    histAxePsi2SRapCut1 -> Write();
    histAxePsi2SRapCut2 -> Write();
    histAxePsi2SRapCut3 -> Write();
    histAxePsi2SRapCut4 -> Write();
    fOut -> Close();
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