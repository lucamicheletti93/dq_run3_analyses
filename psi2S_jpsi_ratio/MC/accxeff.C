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

    //TFile *fIn = new TFile("/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/MC/central_production/AnalysisResults.root", "READ");
    TFile *fIn = new TFile("/Users/lucamicheletti/alice/local_train_test_mc/LHC24e5/AnalysisResults_dq_efficiency.root", "READ");

    TList *listGen1 = (TList*) fIn -> Get("analysis-same-event-pairing/output");
    TList *listJpsiGen2 = (TList*) listGen1 -> FindObject("MCTruthGen_Jpsi");
    TH2D *histJpsiPtYGen = (TH2D*) listJpsiGen2 -> FindObject("Pt_Rapidity");
    TH1D *histJpsiPtGen = (TH1D*) histJpsiPtYGen -> ProjectionX("histJpsiPtGen");
    TH1D *histJpsiYGen = (TH1D*) histJpsiPtYGen -> ProjectionY("histJpsiYGen");

    histJpsiPtGen -> SetLineColor(kBlack);
    histJpsiYGen -> SetLineColor(kBlack);
    histJpsiPtGen -> SetMarkerColor(kBlack);
    histJpsiYGen -> SetMarkerColor(kBlack);

    TList *listPsi2SGen2 = (TList*) listGen1 -> FindObject("MCTruthGen_Psi2S");
    TH2D *histPsi2SPtYGen = (TH2D*) listPsi2SGen2 -> FindObject("Pt_Rapidity");
    TH1D *histPsi2SPtGen = (TH1D*) histPsi2SPtYGen -> ProjectionX("histPsi2SPtGen");
    TH1D *histPsi2SYGen = (TH1D*) histPsi2SPtYGen -> ProjectionY("histPsi2SYGen");

    histPsi2SPtGen -> SetLineColor(kBlack);
    histPsi2SYGen -> SetLineColor(kBlack);
    histPsi2SPtGen -> SetMarkerColor(kBlack);
    histPsi2SYGen -> SetMarkerColor(kBlack);

    TList *listRec1 = (TList*) fIn -> Get("analysis-same-event-pairing/output");

    TList *listJpsiRecCut1 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_matchedMchMid_mumuFromJpsi");
    TList *listJpsiRecCut2 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt10SigmaPDCA_mumuFromJpsi");
    TList *listJpsiRecCut3 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA_mumuFromJpsi");
    TList *listJpsiRecCut4 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt510SigmaPDCA_mumuFromJpsi");

    THnSparseD *histJpsiMassPtYRecCut1 = (THnSparseD*) listJpsiRecCut1 -> FindObject("Mass_Pt_Rapidity");
    TH1D *histJpsiPtRecCut1 = (TH1D*) histJpsiMassPtYRecCut1 -> Projection(1, "histJpsiPtRecCut1");
    TH1D *histJpsiYRecCut1 = (TH1D*) histJpsiMassPtYRecCut1 -> Projection(2, "histJpsiYRecCut1");

    THnSparseD *histJpsiMassPtYRecCut2 = (THnSparseD*) listJpsiRecCut2 -> FindObject("Mass_Pt_Rapidity");
    TH1D *histJpsiPtRecCut2 = (TH1D*) histJpsiMassPtYRecCut2 -> Projection(1, "histJpsiPtRecCut2");
    TH1D *histJpsiYRecCut2 = (TH1D*) histJpsiMassPtYRecCut2 -> Projection(2, "histJpsiYRecCut2");

    THnSparseD *histJpsiMassPtYRecCut3 = (THnSparseD*) listJpsiRecCut3 -> FindObject("Mass_Pt_Rapidity");
    TH1D *histJpsiPtRecCut3 = (TH1D*) histJpsiMassPtYRecCut3 -> Projection(1, "histJpsiPtRecCut3");
    TH1D *histJpsiYRecCut3 = (TH1D*) histJpsiMassPtYRecCut3 -> Projection(2, "histJpsiYRecCut3");

    THnSparseD *histJpsiMassPtYRecCut4 = (THnSparseD*) listJpsiRecCut4 -> FindObject("Mass_Pt_Rapidity");
    TH1D *histJpsiPtRecCut4 = (TH1D*) histJpsiMassPtYRecCut4 -> Projection(1, "histJpsiPtRecCut4");
    TH1D *histJpsiYRecCut4 = (TH1D*) histJpsiMassPtYRecCut4 -> Projection(2, "histJpsiYRecCut4");

    histJpsiPtRecCut1 -> SetLineColor(kRed+1);
    histJpsiYRecCut1 -> SetLineColor(kRed+1);
    histJpsiPtRecCut1 -> SetMarkerColor(kRed+1);
    histJpsiYRecCut1 -> SetMarkerColor(kRed+1);

    TList *listPsi2SRecCut1 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_matchedMchMid_mumuFromPsi2S");
    TList *listPsi2SRecCut2 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt10SigmaPDCA_mumuFromPsi2S");
    TList *listPsi2SRecCut3 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA_mumuFromPsi2S");
    TList *listPsi2SRecCut4 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt510SigmaPDCA_mumuFromPsi2S");

    THnSparseD *histPsi2SMassPtYRecCut1 = (THnSparseD*) listPsi2SRecCut1 -> FindObject("Mass_Pt_Rapidity");
    TH1D *histPsi2SPtRecCut1 = (TH1D*) histPsi2SMassPtYRecCut1 -> Projection(1, "histPsi2SPtRecCut1");
    TH1D *histPsi2SYRecCut1 = (TH1D*) histPsi2SMassPtYRecCut1 -> Projection(2, "histPsi2SYRecCut1");

    THnSparseD *histPsi2SMassPtYRecCut2 = (THnSparseD*) listPsi2SRecCut2 -> FindObject("Mass_Pt_Rapidity");
    TH1D *histPsi2SPtRecCut2 = (TH1D*) histPsi2SMassPtYRecCut2 -> Projection(1, "histPsi2SPtRecCut2");
    TH1D *histPsi2SYRecCut2 = (TH1D*) histPsi2SMassPtYRecCut2 -> Projection(2, "histPsi2SYRecCut2");

    THnSparseD *histPsi2SMassPtYRecCut3 = (THnSparseD*) listPsi2SRecCut3 -> FindObject("Mass_Pt_Rapidity");
    TH1D *histPsi2SPtRecCut3 = (TH1D*) histPsi2SMassPtYRecCut3 -> Projection(1, "histPsi2SPtRecCut3");
    TH1D *histPsi2SYRecCut3 = (TH1D*) histPsi2SMassPtYRecCut3 -> Projection(2, "histPsi2SYRecCut3");

    THnSparseD *histPsi2SMassPtYRecCut4 = (THnSparseD*) listPsi2SRecCut4 -> FindObject("Mass_Pt_Rapidity");
    TH1D *histPsi2SPtRecCut4 = (TH1D*) histPsi2SMassPtYRecCut4 -> Projection(1, "histPsi2SPtRecCut4");
    TH1D *histPsi2SYRecCut4 = (TH1D*) histPsi2SMassPtYRecCut4 -> Projection(2, "histPsi2SYRecCut4");

    histPsi2SPtRecCut1 -> SetLineColor(kRed+1);
    histPsi2SYRecCut1 -> SetLineColor(kRed+1);
    histPsi2SPtRecCut1 -> SetMarkerColor(kRed+1);
    histPsi2SYRecCut1 -> SetMarkerColor(kRed+1);

    // Rebin the histograms
    TH1D *histJpsiPtGenRebin = (TH1D*) histJpsiPtGen -> Rebin(18, "histJpsiPtGenRebin", ptBinsRun2); 
    TH1D *histJpsiPtRecCut1Rebin = (TH1D*) histJpsiPtRecCut1 -> Rebin(18, "histJpsiPtRecCut1Rebin", ptBinsRun2); 
    TH1D *histJpsiPtRecCut2Rebin = (TH1D*) histJpsiPtRecCut2 -> Rebin(18, "histJpsiPtRecCut2Rebin", ptBinsRun2); 
    TH1D *histJpsiPtRecCut3Rebin = (TH1D*) histJpsiPtRecCut3 -> Rebin(18, "histJpsiPtRecCut3Rebin", ptBinsRun2); 
    TH1D *histJpsiPtRecCut4Rebin = (TH1D*) histJpsiPtRecCut4 -> Rebin(18, "histJpsiPtRecCut4Rebin", ptBinsRun2); 

    histJpsiPtRecCut1Rebin -> SetLineColor(kRed+1);
    histJpsiPtRecCut2Rebin -> SetLineColor(kOrange+7);
    histJpsiPtRecCut3Rebin -> SetLineColor(kAzure+2);
    histJpsiPtRecCut4Rebin -> SetLineColor(kBlue+1);

    histJpsiPtRecCut1Rebin -> SetLineWidth(2);
    histJpsiPtRecCut2Rebin -> SetLineWidth(2);
    histJpsiPtRecCut3Rebin -> SetLineWidth(2);
    histJpsiPtRecCut4Rebin -> SetLineWidth(2);

    histJpsiPtRecCut1Rebin -> SetMarkerColor(kRed+1);
    histJpsiPtRecCut2Rebin -> SetMarkerColor(kOrange+7);
    histJpsiPtRecCut3Rebin -> SetMarkerColor(kAzure+2);
    histJpsiPtRecCut4Rebin -> SetMarkerColor(kBlue+1);

    TH1D *histJpsiPtAxeCut1 = (TH1D*) histJpsiPtRecCut1Rebin -> Clone("histJpsiPtAxeCut1");
    histJpsiPtAxeCut1 -> Divide(histJpsiPtGenRebin);
    histJpsiPtAxeCut1 -> SetLineColor(kRed+1);
    histJpsiPtAxeCut1 -> SetMarkerColor(kRed+1);

    TH1D *histJpsiPtAxeCut2 = (TH1D*) histJpsiPtRecCut2Rebin -> Clone("histJpsiPtAxeCut2");
    histJpsiPtAxeCut2 -> Divide(histJpsiPtGenRebin);
    histJpsiPtAxeCut2 -> SetLineColor(kOrange+7);
    histJpsiPtAxeCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histJpsiPtAxeCut3 = (TH1D*) histJpsiPtRecCut3Rebin -> Clone("histJpsiPtAxeCut3");
    histJpsiPtAxeCut3 -> Divide(histJpsiPtGenRebin);
    histJpsiPtAxeCut3 -> SetLineColor(kAzure+2);
    histJpsiPtAxeCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histJpsiPtAxeCut4 = (TH1D*) histJpsiPtRecCut4Rebin -> Clone("histJpsiPtAxeCut4");
    histJpsiPtAxeCut4 -> Divide(histJpsiPtGenRebin);
    histJpsiPtAxeCut4 -> SetLineColor(kBlue+1);
    histJpsiPtAxeCut4 -> SetMarkerColor(kBlue+1);

    TH1D *histJpsiYGenRebin = (TH1D*) histJpsiYGen -> Rebin(6, "histJpsiYGenRebin", yBinsRun2); 
    TH1D *histJpsiYRecCut1Rebin = (TH1D*) histJpsiYRecCut1 -> Rebin(6, "histJpsiYRecCut1Rebin", yBinsRun2); 
    TH1D *histJpsiYRecCut2Rebin = (TH1D*) histJpsiYRecCut2 -> Rebin(6, "histJpsiYRecCut2Rebin", yBinsRun2); 
    TH1D *histJpsiYRecCut3Rebin = (TH1D*) histJpsiYRecCut3 -> Rebin(6, "histJpsiYRecCut3Rebin", yBinsRun2); 
    TH1D *histJpsiYRecCut4Rebin = (TH1D*) histJpsiYRecCut4 -> Rebin(6, "histJpsiYRecCut4Rebin", yBinsRun2); 

    histJpsiYRecCut1Rebin -> SetLineColor(kRed+1);
    histJpsiYRecCut2Rebin -> SetLineColor(kOrange+7);
    histJpsiYRecCut3Rebin -> SetLineColor(kAzure+2);
    histJpsiYRecCut4Rebin -> SetLineColor(kBlue+1);

    histJpsiYRecCut1Rebin -> SetLineWidth(2);
    histJpsiYRecCut2Rebin -> SetLineWidth(2);
    histJpsiYRecCut3Rebin -> SetLineWidth(2);
    histJpsiYRecCut4Rebin -> SetLineWidth(2);

    histJpsiYRecCut1Rebin -> SetMarkerColor(kRed+1);
    histJpsiYRecCut2Rebin -> SetMarkerColor(kOrange+7);
    histJpsiYRecCut3Rebin -> SetMarkerColor(kAzure+2);
    histJpsiYRecCut4Rebin -> SetMarkerColor(kBlue+1);

    TH1D *histJpsiYAxeCut1 = (TH1D*) histJpsiYRecCut1Rebin -> Clone("histJpsiYAxeCut1");
    histJpsiYAxeCut1 -> Divide(histJpsiYGenRebin);
    histJpsiYAxeCut1 -> SetLineColor(kRed+1);
    histJpsiYAxeCut1 -> SetMarkerColor(kRed+1);

    TH1D *histJpsiYAxeCut2 = (TH1D*) histJpsiYRecCut2Rebin -> Clone("histJpsiYAxeCut2");
    histJpsiYAxeCut2 -> Divide(histJpsiYGenRebin);
    histJpsiYAxeCut2 -> SetLineColor(kOrange+7);
    histJpsiYAxeCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histJpsiYAxeCut3 = (TH1D*) histJpsiYRecCut3Rebin -> Clone("histJpsiYAxeCut3");
    histJpsiYAxeCut3 -> Divide(histJpsiYGenRebin);
    histJpsiYAxeCut3 -> SetLineColor(kAzure+2);
    histJpsiYAxeCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histJpsiYAxeCut4 = (TH1D*) histJpsiYRecCut4Rebin -> Clone("histJpsiYAxeCut4");
    histJpsiYAxeCut4 -> Divide(histJpsiYGenRebin);
    histJpsiYAxeCut4 -> SetLineColor(kBlue+1);
    histJpsiYAxeCut4 -> SetMarkerColor(kBlue+1);



    TH1D *histPsi2SPtGenRebin = (TH1D*) histPsi2SPtGen -> Rebin(18, "histPsi2SPtGenRebin", ptBinsRun2); 
    TH1D *histPsi2SPtRecCut1Rebin = (TH1D*) histPsi2SPtRecCut1 -> Rebin(18, "histPsi2SPtRecCut1Rebin", ptBinsRun2); 
    TH1D *histPsi2SPtRecCut2Rebin = (TH1D*) histPsi2SPtRecCut2 -> Rebin(18, "histPsi2SPtRecCut2Rebin", ptBinsRun2); 
    TH1D *histPsi2SPtRecCut3Rebin = (TH1D*) histPsi2SPtRecCut3 -> Rebin(18, "histPsi2SPtRecCut3Rebin", ptBinsRun2); 
    TH1D *histPsi2SPtRecCut4Rebin = (TH1D*) histPsi2SPtRecCut4 -> Rebin(18, "histPsi2SPtRecCut4Rebin", ptBinsRun2); 

    histPsi2SPtRecCut1Rebin -> SetLineColor(kRed+1);
    histPsi2SPtRecCut2Rebin -> SetLineColor(kOrange+7);
    histPsi2SPtRecCut3Rebin -> SetLineColor(kAzure+2);
    histPsi2SPtRecCut4Rebin -> SetLineColor(kBlue+1);

    histPsi2SPtRecCut1Rebin -> SetLineWidth(2);
    histPsi2SPtRecCut2Rebin -> SetLineWidth(2);
    histPsi2SPtRecCut3Rebin -> SetLineWidth(2);
    histPsi2SPtRecCut4Rebin -> SetLineWidth(2);

    histPsi2SPtRecCut1Rebin -> SetMarkerColor(kRed+1);
    histPsi2SPtRecCut2Rebin -> SetMarkerColor(kOrange+7);
    histPsi2SPtRecCut3Rebin -> SetMarkerColor(kAzure+2);
    histPsi2SPtRecCut4Rebin -> SetMarkerColor(kBlue+1);

    TH1D *histPsi2SPtAxeCut1 = (TH1D*) histPsi2SPtRecCut1Rebin -> Clone("histPsi2SPtAxeCut1");
    histPsi2SPtAxeCut1 -> Divide(histPsi2SPtGenRebin);
    histPsi2SPtAxeCut1 -> SetLineColor(kRed+1);
    histPsi2SPtAxeCut1 -> SetMarkerColor(kRed+1);

    TH1D *histPsi2SPtAxeCut2 = (TH1D*) histPsi2SPtRecCut2Rebin -> Clone("histPsi2SPtAxeCut2");
    histPsi2SPtAxeCut2 -> Divide(histPsi2SPtGenRebin);
    histPsi2SPtAxeCut2 -> SetLineColor(kOrange+7);
    histPsi2SPtAxeCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histPsi2SPtAxeCut3 = (TH1D*) histPsi2SPtRecCut3Rebin -> Clone("histPsi2SPtAxeCut3");
    histPsi2SPtAxeCut3 -> Divide(histPsi2SPtGenRebin);
    histPsi2SPtAxeCut3 -> SetLineColor(kAzure+2);
    histPsi2SPtAxeCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histPsi2SPtAxeCut4 = (TH1D*) histPsi2SPtRecCut4Rebin -> Clone("histPsi2SPtAxeCut4");
    histPsi2SPtAxeCut4 -> Divide(histPsi2SPtGenRebin);
    histPsi2SPtAxeCut4 -> SetLineColor(kBlue+1);
    histPsi2SPtAxeCut4 -> SetMarkerColor(kBlue+1);

    TH1D *histPsi2SYGenRebin = (TH1D*) histPsi2SYGen -> Rebin(6, "histPsi2SYGenRebin", yBinsRun2); 
    TH1D *histPsi2SYRecCut1Rebin = (TH1D*) histPsi2SYRecCut1 -> Rebin(6, "histPsi2SYRecCut1Rebin", yBinsRun2); 
    TH1D *histPsi2SYRecCut2Rebin = (TH1D*) histPsi2SYRecCut2 -> Rebin(6, "histPsi2SYRecCut2Rebin", yBinsRun2); 
    TH1D *histPsi2SYRecCut3Rebin = (TH1D*) histPsi2SYRecCut3 -> Rebin(6, "histPsi2SYRecCut3Rebin", yBinsRun2); 
    TH1D *histPsi2SYRecCut4Rebin = (TH1D*) histPsi2SYRecCut4 -> Rebin(6, "histPsi2SYRecCut4Rebin", yBinsRun2); 

    histPsi2SYRecCut1Rebin -> SetLineColor(kRed+1);
    histPsi2SYRecCut2Rebin -> SetLineColor(kOrange+7);
    histPsi2SYRecCut3Rebin -> SetLineColor(kAzure+2);
    histPsi2SYRecCut4Rebin -> SetLineColor(kBlue+1);

    histPsi2SYRecCut1Rebin -> SetLineWidth(2);
    histPsi2SYRecCut2Rebin -> SetLineWidth(2);
    histPsi2SYRecCut3Rebin -> SetLineWidth(2);
    histPsi2SYRecCut4Rebin -> SetLineWidth(2);

    histPsi2SYRecCut1Rebin -> SetMarkerColor(kRed+1);
    histPsi2SYRecCut2Rebin -> SetMarkerColor(kOrange+7);
    histPsi2SYRecCut3Rebin -> SetMarkerColor(kAzure+2);
    histPsi2SYRecCut4Rebin -> SetMarkerColor(kBlue+1);

    TH1D *histPsi2SYAxeCut1 = (TH1D*) histPsi2SYRecCut1Rebin -> Clone("histPsi2SYAxeCut1");
    histPsi2SYAxeCut1 -> Divide(histPsi2SYGenRebin);
    histPsi2SYAxeCut1 -> SetLineColor(kRed+1);
    histPsi2SYAxeCut1 -> SetMarkerColor(kRed+1);

    TH1D *histPsi2SYAxeCut2 = (TH1D*) histPsi2SYRecCut2Rebin -> Clone("histPsi2SYAxeCut2");
    histPsi2SYAxeCut2 -> Divide(histPsi2SYGenRebin);
    histPsi2SYAxeCut2 -> SetLineColor(kOrange+7);
    histPsi2SYAxeCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histPsi2SYAxeCut3 = (TH1D*) histPsi2SYRecCut3Rebin -> Clone("histPsi2SYAxeCut3");
    histPsi2SYAxeCut3 -> Divide(histPsi2SYGenRebin);
    histPsi2SYAxeCut3 -> SetLineColor(kAzure+2);
    histPsi2SYAxeCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histPsi2SYAxeCut4 = (TH1D*) histPsi2SYRecCut4Rebin -> Clone("histPsi2SYAxeCut4");
    histPsi2SYAxeCut4 -> Divide(histPsi2SYGenRebin);
    histPsi2SYAxeCut4 -> SetLineColor(kBlue+1);
    histPsi2SYAxeCut4 -> SetMarkerColor(kBlue+1);


    TCanvas *canvasJpsi = new TCanvas("canvasJpsi", "", 1200, 1200);
    canvasJpsi -> Divide(2, 2);

    canvasJpsi -> cd(1);
    gPad -> SetLogy(1);
    histJpsiPtGenRebin -> SetTitle("");
    histJpsiPtGenRebin -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histJpsiPtGenRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histJpsiPtGenRebin -> GetYaxis() -> SetRangeUser(1, 1e6);
    histJpsiPtGenRebin -> Draw("EP");
    histJpsiPtRecCut1Rebin -> Draw("EP SAME");
    histJpsiPtRecCut2Rebin -> Draw("EP SAME");
    histJpsiPtRecCut3Rebin -> Draw("EP SAME");
    histJpsiPtRecCut4Rebin -> Draw("EP SAME");

    canvasJpsi -> cd(2);
    gPad -> SetLogy(1);
    histJpsiYGenRebin -> SetTitle("");
    histJpsiYGenRebin -> GetXaxis() -> SetTitle("#it{y}");
    histJpsiYGenRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histJpsiYGenRebin -> GetYaxis() -> SetRangeUser(1e3, 1e6);
    histJpsiYGenRebin -> Draw("EP");
    histJpsiYRecCut1Rebin -> Draw("EP SAME");
    histJpsiYRecCut2Rebin -> Draw("EP SAME");
    histJpsiYRecCut3Rebin -> Draw("EP SAME");
    histJpsiYRecCut4Rebin -> Draw("EP SAME");

    canvasJpsi -> cd(3);
    histJpsiPtAxeCut1 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histJpsiPtAxeCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histJpsiPtAxeCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histJpsiPtAxeCut1 -> Draw("EP");
    histJpsiPtAxeCut2 -> Draw("EP SAME");
    histJpsiPtAxeCut3 -> Draw("EP SAME");
    histJpsiPtAxeCut4 -> Draw("EP SAME");
    histPtAxeRun2 -> Draw("EP SAME");

    TLegend *legendJpsiPtAxe = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendJpsiPtAxe);
    legendJpsiPtAxe -> AddEntry(histPtAxeRun2, "Run 2", "PL");
    legendJpsiPtAxe -> AddEntry(histJpsiPtAxeCut1, "Run 3 - matchedMchMid", "PL");
    legendJpsiPtAxe -> AddEntry(histJpsiPtAxeCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendJpsiPtAxe -> AddEntry(histJpsiPtAxeCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendJpsiPtAxe -> AddEntry(histJpsiPtAxeCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendJpsiPtAxe -> Draw("SAME");

    canvasJpsi -> cd(4);
    histJpsiYAxeCut1 -> GetXaxis() -> SetTitle("#it{y}");
    histJpsiYAxeCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histJpsiYAxeCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histJpsiYAxeCut1 -> Draw("EP");
    histJpsiYAxeCut2 -> Draw("EP SAME");
    histJpsiYAxeCut3 -> Draw("EP SAME");
    histJpsiYAxeCut4 -> Draw("EP SAME");
    histYAxeRun2 -> Draw("EP SAME");

    TLegend *legendJpsiYAxe = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendJpsiYAxe);
    legendJpsiYAxe -> AddEntry(histYAxeRun2, "Run 2", "PL");
    legendJpsiYAxe -> AddEntry(histJpsiYAxeCut1, "Run 3 - matchedMchMid", "PL");
    legendJpsiYAxe -> AddEntry(histJpsiYAxeCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendJpsiYAxe -> AddEntry(histJpsiYAxeCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendJpsiYAxe -> AddEntry(histJpsiYAxeCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendJpsiYAxe -> Draw("SAME");


    TCanvas *canvasPsi2S = new TCanvas("canvasPsi2S", "", 1200, 1200);
    canvasPsi2S -> Divide(2, 2);

    canvasPsi2S -> cd(1);
    gPad -> SetLogy(1);
    histPsi2SPtGenRebin -> SetTitle("");
    histPsi2SPtGenRebin -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histPsi2SPtGenRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histPsi2SPtGenRebin -> GetYaxis() -> SetRangeUser(1, 1e6);
    histPsi2SPtGenRebin -> Draw("EP");
    histPsi2SPtRecCut1Rebin -> Draw("EP SAME");
    histPsi2SPtRecCut2Rebin -> Draw("EP SAME");
    histPsi2SPtRecCut3Rebin -> Draw("EP SAME");
    histPsi2SPtRecCut4Rebin -> Draw("EP SAME");

    canvasPsi2S -> cd(2);
    gPad -> SetLogy(1);
    histPsi2SYGenRebin -> SetTitle("");
    histPsi2SYGenRebin -> GetXaxis() -> SetTitle("#it{y}");
    histPsi2SYGenRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histPsi2SYGenRebin -> GetYaxis() -> SetRangeUser(1e3, 1e6);
    histPsi2SYGenRebin -> Draw("EP");
    histPsi2SYRecCut1Rebin -> Draw("EP SAME");
    histPsi2SYRecCut2Rebin -> Draw("EP SAME");
    histPsi2SYRecCut3Rebin -> Draw("EP SAME");
    histPsi2SYRecCut4Rebin -> Draw("EP SAME");

    canvasPsi2S -> cd(3);
    histPsi2SPtAxeCut1 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histPsi2SPtAxeCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histPsi2SPtAxeCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histPsi2SPtAxeCut1 -> Draw("EP");
    histPsi2SPtAxeCut2 -> Draw("EP SAME");
    histPsi2SPtAxeCut3 -> Draw("EP SAME");
    histPsi2SPtAxeCut4 -> Draw("EP SAME");

    TLegend *legendPsi2SPtAxe = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendPsi2SPtAxe);
    legendPsi2SPtAxe -> AddEntry(histPsi2SPtAxeCut1, "Run 3 - matchedMchMid", "PL");
    legendPsi2SPtAxe -> AddEntry(histPsi2SPtAxeCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendPsi2SPtAxe -> AddEntry(histPsi2SPtAxeCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendPsi2SPtAxe -> AddEntry(histPsi2SPtAxeCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendPsi2SPtAxe -> Draw("SAME");

    canvasPsi2S -> cd(4);
    histPsi2SYAxeCut1 -> GetXaxis() -> SetTitle("#it{y}");
    histPsi2SYAxeCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histPsi2SYAxeCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histPsi2SYAxeCut1 -> Draw("EP");
    histPsi2SYAxeCut2 -> Draw("EP SAME");
    histPsi2SYAxeCut3 -> Draw("EP SAME");
    histPsi2SYAxeCut4 -> Draw("EP SAME");

    TLegend *legendPsi2SYAxe = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendPsi2SYAxe);
    legendPsi2SYAxe -> AddEntry(histPsi2SYAxeCut1, "Run 3 - matchedMchMid", "PL");
    legendPsi2SYAxe -> AddEntry(histPsi2SYAxeCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendPsi2SYAxe -> AddEntry(histPsi2SYAxeCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendPsi2SYAxe -> AddEntry(histPsi2SYAxeCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendPsi2SYAxe -> Draw("SAME");

    // Axe ratio
    TH1D *histRatioPtAxeCut1 = (TH1D*) histPsi2SPtAxeCut1 -> Clone("histRatioPtAxeCut1");
    histRatioPtAxeCut1 -> Divide(histJpsiPtAxeCut1);
    histRatioPtAxeCut1 -> SetLineColor(kRed+1);
    histRatioPtAxeCut1 -> SetMarkerColor(kRed+1);

    TH1D *histRatioPtAxeCut2 = (TH1D*) histPsi2SPtAxeCut2 -> Clone("histRatioPtAxeCut2");
    histRatioPtAxeCut2 -> Divide(histJpsiPtAxeCut2);
    histRatioPtAxeCut2 -> SetLineColor(kOrange+7);
    histRatioPtAxeCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histRatioPtAxeCut3 = (TH1D*) histPsi2SPtAxeCut3 -> Clone("histRatioPtAxeCut3");
    histRatioPtAxeCut3 -> Divide(histJpsiPtAxeCut3);
    histRatioPtAxeCut3 -> SetLineColor(kAzure+2);
    histRatioPtAxeCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histRatioPtAxeCut4 = (TH1D*) histPsi2SPtAxeCut4 -> Clone("histRatioPtAxeCut4");
    histRatioPtAxeCut4 -> Divide(histJpsiPtAxeCut4);
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



    TH1D *histRatioYAxeCut1 = (TH1D*) histPsi2SYAxeCut1 -> Clone("histRatioYAxeCut1");
    histRatioYAxeCut1 -> Divide(histJpsiYAxeCut1);
    histRatioYAxeCut1 -> SetLineColor(kRed+1);
    histRatioYAxeCut1 -> SetMarkerColor(kRed+1);

    TH1D *histRatioYAxeCut2 = (TH1D*) histPsi2SYAxeCut2 -> Clone("histRatioYAxeCut2");
    histRatioYAxeCut2 -> Divide(histJpsiYAxeCut2);
    histRatioYAxeCut2 -> SetLineColor(kOrange+7);
    histRatioYAxeCut2 -> SetMarkerColor(kOrange+7);

    TH1D *histRatioYAxeCut3 = (TH1D*) histPsi2SYAxeCut3 -> Clone("histRatioYAxeCut3");
    histRatioYAxeCut3 -> Divide(histJpsiYAxeCut3);
    histRatioYAxeCut3 -> SetLineColor(kAzure+2);
    histRatioYAxeCut3 -> SetMarkerColor(kAzure+2);

    TH1D *histRatioYAxeCut4 = (TH1D*) histPsi2SYAxeCut4 -> Clone("histRatioYAxeCut4");
    histRatioYAxeCut4 -> Divide(histJpsiYAxeCut4);
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



    /*
    TCanvas *canvasPtGenRec = new TCanvas("canvasPtGenRec", "", 800, 600);
    gPad -> SetLogy(1);
    histJpsiPtGenRebin -> SetTitle("");
    histJpsiPtGenRebin -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histJpsiPtGenRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histJpsiPtGenRebin -> GetYaxis() -> SetRangeUser(1, 1e6);
    histJpsiPtGenRebin -> Draw("EP");
    histJpsiPtRecCut1Rebin -> Draw("EP SAME");
    histJpsiPtRecCut2Rebin -> Draw("EP SAME");
    histJpsiPtRecCut3Rebin -> Draw("EP SAME");
    histJpsiPtRecCut4Rebin -> Draw("EP SAME");

    TCanvas *canvasYGenRec = new TCanvas("canvasYGenRec", "", 800, 600);
    gPad -> SetLogy(1);
    histJpsiYGenRebin -> SetTitle("");
    histJpsiYGenRebin -> GetXaxis() -> SetTitle("#it{y}");
    histJpsiYGenRebin -> GetYaxis() -> SetTitle("N_{J/#psi}");
    histJpsiYGenRebin -> GetYaxis() -> SetRangeUser(1e3, 1e6);
    histJpsiYGenRebin -> Draw("EP");
    histJpsiYRecCut1Rebin -> Draw("EP SAME");
    histJpsiYRecCut2Rebin -> Draw("EP SAME");
    histJpsiYRecCut3Rebin -> Draw("EP SAME");
    histJpsiYRecCut4Rebin -> Draw("EP SAME");

    TCanvas *canvasPtAxe = new TCanvas("canvasPtAxe", "", 800, 600);
    histJpsiPtAxeCut1 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histJpsiPtAxeCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histJpsiPtAxeCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histJpsiPtAxeCut1 -> Draw("EP");
    histJpsiPtAxeCut2 -> Draw("EP SAME");
    histJpsiPtAxeCut3 -> Draw("EP SAME");
    histJpsiPtAxeCut4 -> Draw("EP SAME");
    histPtAxeRun2 -> Draw("EP SAME");

    TLegend *legendPtAxe = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendPtAxe);
    legendPtAxe -> AddEntry(histPtAxeRun2, "Run 2", "PL");
    legendPtAxe -> AddEntry(histJpsiPtAxeCut1, "Run 3 - matchedMchMid", "PL");
    legendPtAxe -> AddEntry(histJpsiPtAxeCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendPtAxe -> AddEntry(histJpsiPtAxeCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendPtAxe -> AddEntry(histJpsiPtAxeCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendPtAxe -> Draw("SAME");

    TCanvas *canvasYAxe = new TCanvas("canvasYAxe", "", 800, 600);
    histJpsiYAxeCut1 -> GetXaxis() -> SetTitle("#it{y}");
    histJpsiYAxeCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histJpsiYAxeCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histJpsiYAxeCut1 -> Draw("EP");
    histJpsiYAxeCut2 -> Draw("EP SAME");
    histJpsiYAxeCut3 -> Draw("EP SAME");
    histJpsiYAxeCut4 -> Draw("EP SAME");
    histYAxeRun2 -> Draw("EP SAME");

    TLegend *legendYAxe = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendYAxe);
    legendYAxe -> AddEntry(histYAxeRun2, "Run 2", "PL");
    legendYAxe -> AddEntry(histJpsiYAxeCut1, "Run 3 - matchedMchMid", "PL");
    legendYAxe -> AddEntry(histJpsiYAxeCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendYAxe -> AddEntry(histJpsiYAxeCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendYAxe -> AddEntry(histJpsiYAxeCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendYAxe -> Draw("SAME");
    */


    //canvasPtGenRec -> SaveAs("gen_rec_distributions_vs_pt.pdf");
    //canvasYGenRec -> SaveAs("gen_rec_distributions_vs_y.pdf");
    //canvasPtAxe -> SaveAs("Axe_distributions_vs_pt.pdf");
    //canvasYAxe -> SaveAs("Axe_distributions_vs_y.pdf");
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