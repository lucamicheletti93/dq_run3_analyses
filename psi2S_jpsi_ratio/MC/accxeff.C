void LoadStyle();
void SetLegend(TLegend * );
void SetHistogram(TH1D *, int , int , int , double );

void accxeff() {
    LoadStyle();
    double ptBinsRun3[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0};
    //-----------------------------------------------------------//
    // https://alice-notes.web.cern.ch/system/files/notes/analysis/497/2017-Aug-11-analysis_note-pp13TeV-analysis_note.pdf
    // Run 2 Axe
    double ptBinsJpsiRun2[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 26.0, 30.0};
    double genJpsiPtRun2[] = {3101603976, 7943246726, 20952356646, 17586753168, 11232729979, 6444395948, 3592011595, 1985088607, 1128815038, 1037545983, 382478109, 158022187, 73955203, 32679107, 17059015, 8815662, 8074816, 3038991};
    double recJpsiPtRun2[] = {462849972, 1094777590, 2516558340, 2187868332, 1703719599, 1252770554, 867272809, 564602236, 363957858, 384111030, 158673139, 72532504, 34624481, 16405609, 8519213, 4592152, 4181192, 1354967};
    double axeJpsiPtRun2[] = {0.14923, 0.13782, 0.12011, 0.12440, 0.15167, 0.19440, 0.24144, 0.28442, 0.32242, 0.37021, 0.41486, 0.45900, 0.46818, 0.50202, 0.49940, 0.52091, 0.51781, 0.44586};
    double errAxeJpsiPtRun2[] = {0.00089, 0.00054, 0.00031, 0.00034, 0.00047, 0.00068, 0.00099, 0.00140, 0.00193, 0.00208, 0.00349, 0.00549, 0.00804, 0.01212, 0.01678, 0.02332, 0.02437, 0.03951};

    TH1D *histGenJpsiPtRun2 = new TH1D("histGenJpsiPtRun2", "", 18, ptBinsJpsiRun2);
    TH1D *histRecJpsiPtRun2 = new TH1D("histRecJpsiPtRun2", "", 18, ptBinsJpsiRun2);
    TH1D *histAxeJpsiPtRun2 = new TH1D("histAxeJpsiPtRun2", "", 18, ptBinsJpsiRun2);
    SetHistogram(histAxeJpsiPtRun2, 1, 1, 20, 1);

    for (int iPt = 0;iPt < 18;iPt++) {
        histGenJpsiPtRun2 -> SetBinContent(iPt+1, genJpsiPtRun2[iPt]);
        histGenJpsiPtRun2 -> SetBinError(iPt+1, TMath::Sqrt(genJpsiPtRun2[iPt]));
        histRecJpsiPtRun2 -> SetBinContent(iPt+1, recJpsiPtRun2[iPt]);
        histRecJpsiPtRun2 -> SetBinError(iPt+1, TMath::Sqrt(recJpsiPtRun2[iPt]));
        histAxeJpsiPtRun2 -> SetBinContent(iPt+1, axeJpsiPtRun2[iPt]);
        histAxeJpsiPtRun2 -> SetBinError(iPt+1, errAxeJpsiPtRun2[iPt]);
    }

    double ptBinsPsi2SRun2[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0};
    double axePsi2SPtRun2[] = {0.19816, 0.17608, 0.16092, 0.16698, 0.18982, 0.22608, 0.26571, 0.30449, 0.34395, 0.39066, 0.43539, 0.45861, 0.47282};
    double errAxePsi2SPtRun2[] = {0.00055, 0.00038, 0.00040, 0.00051, 0.00071, 0.00101, 0.00144, 0.00200, 0.00215, 0.00362, 0.00570, 0.00854, 0.01009};

    TH1D *histAxePsi2SPtRun2 = new TH1D("histAxePsi2SPtRun2", "", 13, ptBinsPsi2SRun2);
    SetHistogram(histAxePsi2SPtRun2, 1, 1, 20, 1);
    for (int iPt = 0;iPt < 13;iPt++) {
        histAxePsi2SPtRun2 -> SetBinContent(iPt+1, axePsi2SPtRun2[iPt]);
        histAxePsi2SPtRun2 -> SetBinError(iPt+1, errAxePsi2SPtRun2[iPt]);
    }

    double rapBinsJpsiRun2[] = {2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00};
    double axeJpsiRapRun2[] = {0.04395, 0.15815, 0.23390, 0.24372, 0.18706, 0.07195};
    double errAxeJpsiRapRun2[] = {0.00023, 0.00042, 0.00051, 0.00054, 0.00052, 0.00036};

    TH1D *histAxeJpsiRapRun2 = new TH1D("histAxeJpsiRapRun2", "", 6, rapBinsJpsiRun2);
    SetHistogram(histAxeJpsiRapRun2, 1, 1, 20, 1);
    for (int iY = 0;iY < 6;iY++) {
        histAxeJpsiRapRun2 -> SetBinContent(iY+1, axeJpsiRapRun2[iY]);
        histAxeJpsiRapRun2 -> SetBinError(iY+1, errAxeJpsiRapRun2[iY]);
    }

    double rapBinsPsi2SRun2[] = {2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00};
    double axePsi2SRapRun2[] = {0.05670, 0.19625, 0.28420, 0.29129, 0.21998, 0.08255};
    double errAxePsi2SRapRun2[] = {0.00027, 0.00048, 0.00057, 0.00060, 0.00058, 0.00040};

    TH1D *histAxePsi2SRapRun2 = new TH1D("histAxePsi2SRapRun2", "", 6, rapBinsPsi2SRun2);
    SetHistogram(histAxePsi2SRapRun2, 1, 1, 20, 1);
    for (int iY = 0;iY < 6;iY++) {
        histAxePsi2SRapRun2 -> SetBinContent(iY+1, axePsi2SRapRun2[iY]);
        histAxePsi2SRapRun2 -> SetBinError(iY+1, errAxePsi2SRapRun2[iY]);
    }

    TH1D *histRatioAxeRapRun2 = (TH1D*) histAxePsi2SRapRun2 -> Clone("histRatioAxeRapRun2");
    histRatioAxeRapRun2 -> Divide(histAxeJpsiRapRun2);
    SetHistogram(histRatioAxeRapRun2, 1, 1, 20, 1);

    // Rebin to compute the ratio Axe Psi(2S) / J/psi
    TH1D *histGenJpsiPtRun2Rebin = (TH1D*) histGenJpsiPtRun2 -> Rebin(13, "histGenJpsiPtRun2Rebin", ptBinsPsi2SRun2); 
    TH1D *histRecJpsiPtRun2Rebin = (TH1D*) histRecJpsiPtRun2 -> Rebin(13, "histRecJpsiPtRun2Rebin", ptBinsPsi2SRun2); 

    TH1D *histAxeJpsiPtRun2Rebin = new TH1D("histAxeJpsiPtRun2Rebin", "", 13, ptBinsPsi2SRun2);
    histAxeJpsiPtRun2Rebin -> Divide(histRecJpsiPtRun2Rebin, histGenJpsiPtRun2Rebin, 1, 1, "B");
    SetHistogram(histAxeJpsiPtRun2Rebin, 1, 1, 20, 1);

    TH1D *histRatioAxePtRun2 = (TH1D*) histAxePsi2SPtRun2 -> Clone("histRatioAxePtRun2");
    histRatioAxePtRun2 -> Divide(histAxeJpsiPtRun2Rebin);
    SetHistogram(histRatioAxePtRun2, 1, 1, 20, 1);

    // Run 3 Axe (QM2023 preliminary)
    double ptBinsJpsiRun3Prel[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0};
    double axeJpsiPtRun3Prel[] = {0.262, 0.252, 0.252, 0.268, 0.307, 0.371, 0.464, 0.546};
    double statAxeJpsiPtRun3Prel[] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.002, 0.004};
    double systAxeJpsiPtRun3Prel[] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.002, 0.004};

    double ptBinsPsi2SRun3Prel[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0};
    double axePsi2SPtRun3Prel[] = {0.284, 0.275, 0.271, 0.275, 0.296, 0.345, 0.432, 0.545};
    double statAxePsi2SPtRun3Prel[] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.002, 0.003};
    double systAxePsi2SPtRun3Prel[] = {0.000, 0.000, 0.000, 0.000, 0.001, 0.001, 0.001, 0.004};
    
    double ratioAxePtRun3Prel[] = {1.082, 1.094, 1.076, 1.027, 0.963, 0.930, 0.932, 1.000};
    double statRatioAxePtRun3Prel[] = {0.005, 0.004, 0.004, 0.004, 0.005, 0.004, 0.006, 0.008};
    double systRatioAxePtRun3Prel[] = {0.010, 0.005, 0.002, 0.004, 0.004, 0.005, 0.003, 0.013};

    TH1D *histAxeJpsiPtRun3Prel = new TH1D("histAxeJpsiPtRun3Prel", "", 8, ptBinsJpsiRun3Prel);
    SetHistogram(histAxeJpsiPtRun3Prel, 1, 1, 24, 1);

    TH1D *histAxePsi2SPtRun3Prel = new TH1D("histAxePsi2SPtRun3Prel", "", 8, ptBinsPsi2SRun3Prel);
    SetHistogram(histAxePsi2SPtRun3Prel, 1, 1, 24, 1);

    TH1D *histRatioAxePtRun3Prel = new TH1D("histRatioAxePtRun3Prel", "", 8, ptBinsPsi2SRun3Prel);
    SetHistogram(histRatioAxePtRun3Prel, 1, 1, 24, 1);

    for (int iPt = 0;iPt < 8;iPt++) {
        histAxeJpsiPtRun3Prel -> SetBinContent(iPt+1, axeJpsiPtRun3Prel[iPt]);
        histAxeJpsiPtRun3Prel -> SetBinError(iPt+1, TMath::Sqrt(statAxeJpsiPtRun3Prel[iPt]*statAxeJpsiPtRun3Prel[iPt] + systAxeJpsiPtRun3Prel[iPt]*systAxeJpsiPtRun3Prel[iPt]));

        histAxePsi2SPtRun3Prel -> SetBinContent(iPt+1, axePsi2SPtRun3Prel[iPt]);
        histAxePsi2SPtRun3Prel -> SetBinError(iPt+1, TMath::Sqrt(statAxePsi2SPtRun3Prel[iPt]*statAxePsi2SPtRun3Prel[iPt] + systAxePsi2SPtRun3Prel[iPt]*systAxePsi2SPtRun3Prel[iPt]));

        histRatioAxePtRun3Prel -> SetBinContent(iPt+1, ratioAxePtRun3Prel[iPt]);
        histRatioAxePtRun3Prel -> SetBinError(iPt+1, TMath::Sqrt(statRatioAxePtRun3Prel[iPt]*statRatioAxePtRun3Prel[iPt] + systRatioAxePtRun3Prel[iPt]*systRatioAxePtRun3Prel[iPt]));
    }

    double rapBinsJpsiRun3Prel[] = {2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00};
    double axeJpsiRapRun3Prel[] = {0.088, 0.306, 0.425, 0.425, 0.317, 0.115};
    double statAxeJpsiRapRun3Prel[] = {0.000, 0.001, 0.001, 0.001, 0.001, 0.001};
    double systAxeJpsiRapRun3Prel[] = {0.002, 0.004, 0.003, 0.003, 0.003, 0.002};

    double rapBinsPsi2SRun3Prel[] = {2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00};
    double axePsi2SRapRun3Prel[] = {0.101, 0.328, 0.449, 0.450, 0.335, 0.121};
    double statAxePsi2SRapRun3Prel[] = {0.000, 0.001, 0.001, 0.001, 0.001, 0.001};
    double systAxePsi2SRapRun3Prel[] = {0.002, 0.005, 0.007, 0.006, 0.004, 0.002};
    
    double ratioAxeRapRun3Prel[] = {1.159, 1.073, 1.058, 1.059, 1.058, 1.050};
    double statRatioAxeRapRun3Prel[] = {0.008, 0.004, 0.003, 0.003, 0.004, 0.008};
    double systRatioAxeRapRun3Prel[] = {0.015, 0.011, 0.012, 0.012, 0.011, 0.012};

    TH1D *histAxeJpsiRapRun3Prel = new TH1D("histAxeJpsiRapRun3Prel", "", 6, rapBinsJpsiRun3Prel);
    SetHistogram(histAxeJpsiRapRun3Prel, 1, 1, 24, 1);

    TH1D *histAxePsi2SRapRun3Prel = new TH1D("histAxePsi2SRapRun3Prel", "", 6, rapBinsPsi2SRun3Prel);
    SetHistogram(histAxePsi2SRapRun3Prel, 1, 1, 24, 1);

    TH1D *histRatioAxeRapRun3Prel = new TH1D("histRatioAxeRapRun3Prel", "", 6, rapBinsPsi2SRun3Prel);
    SetHistogram(histRatioAxeRapRun3Prel, 1, 1, 24, 1);

    for (int iRap = 0;iRap < 6;iRap++) {
        histAxeJpsiRapRun3Prel -> SetBinContent(iRap+1, axeJpsiRapRun3Prel[iRap]);
        histAxeJpsiRapRun3Prel -> SetBinError(iRap+1, TMath::Sqrt(statAxeJpsiRapRun3Prel[iRap]*statAxeJpsiRapRun3Prel[iRap] + systAxeJpsiRapRun3Prel[iRap]*systAxeJpsiRapRun3Prel[iRap]));

        histAxePsi2SRapRun3Prel -> SetBinContent(iRap+1, axePsi2SRapRun3Prel[iRap]);
        histAxePsi2SRapRun3Prel -> SetBinError(iRap+1, TMath::Sqrt(statAxePsi2SRapRun3Prel[iRap]*statAxePsi2SRapRun3Prel[iRap] + systAxePsi2SRapRun3Prel[iRap]*systAxePsi2SRapRun3Prel[iRap]));

        histRatioAxeRapRun3Prel -> SetBinContent(iRap+1, ratioAxeRapRun3Prel[iRap]);
        histRatioAxeRapRun3Prel -> SetBinError(iRap+1, TMath::Sqrt(statRatioAxeRapRun3Prel[iRap]*statRatioAxeRapRun3Prel[iRap] + systRatioAxeRapRun3Prel[iRap]*systRatioAxeRapRun3Prel[iRap]));
    }

    //-----------------------------------------------------------//
    //TFile *fIn = new TFile("/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/MC/central_production/AnalysisResults.root", "READ");
    string productionName = "LHC24e5"; // prompt charmonia at forward https://its.cern.ch/jira/browse/O2-4884
    //string productionName = "LHC24e4"; // non-prompt charmonia at forward https://its.cern.ch/jira/browse/O2-4885
    string associationType = "time_association";
    //string associationType = "std_association";
    TFile *fIn = new TFile(Form("/Users/lucamicheletti/cernbox/JPSI/Run3/MC/%s/AnalysisResults_dq_efficiency_%s.root", productionName.c_str(), associationType.c_str()), "READ");

    TList *listGen1 = (TList*) fIn -> Get("analysis-same-event-pairing/output");
    TList *listGenJpsi2 = (TList*) listGen1 -> FindObject("MCTruthGen_Jpsi");
    TH2D *histGenJpsiPtRap = (TH2D*) listGenJpsi2 -> FindObject("Pt_Rapidity");
    histGenJpsiPtRap -> GetXaxis() -> SetRangeUser(0., 20.); // WARNING! constrain the generated particles to be in the 0 < pT < 20 region
    histGenJpsiPtRap -> GetYaxis() -> SetRangeUser(2.5, 4); // WARNING! constrain the generated particles to be in the 2.5 < y < 4 region

    TH1D *histGenJpsiPt = (TH1D*) histGenJpsiPtRap -> ProjectionX("histGenJpsiPt");
    SetHistogram(histGenJpsiPt, 1, 1, 20, 1);
    TH1D *histGenJpsiRap = (TH1D*) histGenJpsiPtRap -> ProjectionY("histGenJpsiRap");
    SetHistogram(histGenJpsiRap, 1, 1, 20, 1);

    TList *listPsi2SGen2 = (TList*) listGen1 -> FindObject("MCTruthGen_Psi2S");
    TH2D *histPsi2SPtYGen = (TH2D*) listPsi2SGen2 -> FindObject("Pt_Rapidity");
    histPsi2SPtYGen -> GetXaxis() -> SetRangeUser(0., 20.); // WARNING! constrain the generated particles to be in the 0 < pT < 20 region
    histPsi2SPtYGen -> GetYaxis() -> SetRangeUser(2.5, 4); // WARNING! constrain the generated particles to be in the 2.5 < y < 4 region

    TH1D *histGenPsi2SPt = (TH1D*) histPsi2SPtYGen -> ProjectionX("histGenPsi2SPt");
    SetHistogram(histGenPsi2SPt, 1, 1, 20, 1);
    TH1D *histGenPsi2SRap = (TH1D*) histPsi2SPtYGen -> ProjectionY("histGenPsi2SRap");
    SetHistogram(histGenPsi2SRap, 1, 1, 20, 1);

    TList *listRec1 = (TList*) fIn -> Get("analysis-same-event-pairing/output");

    TList *listRecJpsiCut1 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_matchedMchMid_mumuFromJpsi");
    TList *listRecJpsiCut2 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt10SigmaPDCA_mumuFromJpsi");
    TList *listRecJpsiCut3 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA_mumuFromJpsi");
    TList *listRecJpsiCut4 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt510SigmaPDCA_mumuFromJpsi");

    THnSparseD *histRecJpsiMassPtRapCut1 = (THnSparseD*) listRecJpsiCut1 -> FindObject("Mass_Pt_Rapidity");
    histRecJpsiMassPtRapCut1 -> GetAxis(1) -> SetRangeUser(0., 20.);
    histRecJpsiMassPtRapCut1 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecJpsiPtRapCut1 = (TH2D*) histRecJpsiMassPtRapCut1 -> Projection(2, 1, "histRecJpsiPtRapCut1");
    TH1D *histRecJpsiPtCut1 = (TH1D*) histRecJpsiMassPtRapCut1 -> Projection(1, "histRecJpsiPtCut1");
    TH1D *histRecJpsiRapCut1 = (TH1D*) histRecJpsiMassPtRapCut1 -> Projection(2, "histRecJpsiRapCut1");

    THnSparseD *histRecJpsiMassPtRapCut2 = (THnSparseD*) listRecJpsiCut2 -> FindObject("Mass_Pt_Rapidity");
    histRecJpsiMassPtRapCut2 -> GetAxis(1) -> SetRangeUser(0., 20.);
    histRecJpsiMassPtRapCut2 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecJpsiPtRapCut2 = (TH2D*) histRecJpsiMassPtRapCut2 -> Projection(2, 1, "histRecJpsiPtRapCut2");
    TH1D *histRecJpsiPtCut2 = (TH1D*) histRecJpsiMassPtRapCut2 -> Projection(1, "histRecJpsiPtCut2");
    TH1D *histRecJpsiRapCut2 = (TH1D*) histRecJpsiMassPtRapCut2 -> Projection(2, "histRecJpsiRapCut2");

    THnSparseD *histRecJpsiMassPtRapCut3 = (THnSparseD*) listRecJpsiCut3 -> FindObject("Mass_Pt_Rapidity");
    histRecJpsiMassPtRapCut3 -> GetAxis(1) -> SetRangeUser(0., 20.);
    histRecJpsiMassPtRapCut3 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecJpsiPtRapCut3 = (TH2D*) histRecJpsiMassPtRapCut3 -> Projection(2, 1, "histRecJpsiPtRapCut3");
    TH1D *histRecJpsiPtCut3 = (TH1D*) histRecJpsiMassPtRapCut3 -> Projection(1, "histRecJpsiPtCut3");
    TH1D *histRecJpsiRapCut3 = (TH1D*) histRecJpsiMassPtRapCut3 -> Projection(2, "histRecJpsiRapCut3");

    THnSparseD *histRecJpsiMassPtRapCut4 = (THnSparseD*) listRecJpsiCut4 -> FindObject("Mass_Pt_Rapidity");
    histRecJpsiMassPtRapCut4 -> GetAxis(1) -> SetRangeUser(0., 20.);
    histRecJpsiMassPtRapCut4 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecJpsiPtRapCut4 = (TH2D*) histRecJpsiMassPtRapCut4 -> Projection(2, 1, "histRecJpsiPtRapCut4");
    TH1D *histRecJpsiPtCut4 = (TH1D*) histRecJpsiMassPtRapCut4 -> Projection(1, "histRecJpsiPtCut4");
    TH1D *histRecJpsiRapCut4 = (TH1D*) histRecJpsiMassPtRapCut4 -> Projection(2, "histRecJpsiRapCut4");

    TList *listRecPsi2SCut1 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_matchedMchMid_mumuFromPsi2S");
    TList *listRecPsi2SCut2 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt10SigmaPDCA_mumuFromPsi2S");
    TList *listRecPsi2SCut3 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt210SigmaPDCA_mumuFromPsi2S");
    TList *listRecPsi2SCut4 = (TList*) listRec1 -> FindObject("PairsMuonSEPM_muonLowPt510SigmaPDCA_mumuFromPsi2S");

    THnSparseD *histRecPsi2SMassPtRapCut1 = (THnSparseD*) listRecPsi2SCut1 -> FindObject("Mass_Pt_Rapidity");
    histRecPsi2SMassPtRapCut1 -> GetAxis(1) -> SetRangeUser(0., 20.);
    histRecPsi2SMassPtRapCut1 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecPsi2SPtRapCut1 = (TH2D*) histRecPsi2SMassPtRapCut1 -> Projection(2, 1, "histRecPsi2SPtRapCut1");
    TH1D *histRecPsi2SPtCut1 = (TH1D*) histRecPsi2SMassPtRapCut1 -> Projection(1, "histRecPsi2SPtCut1");
    TH1D *histRecPsi2SRapCut1 = (TH1D*) histRecPsi2SMassPtRapCut1 -> Projection(2, "histRecPsi2SRapCut1");

    THnSparseD *histRecPsi2SMassPtRapCut2 = (THnSparseD*) listRecPsi2SCut2 -> FindObject("Mass_Pt_Rapidity");
    histRecPsi2SMassPtRapCut2 -> GetAxis(1) -> SetRangeUser(0., 20.);
    histRecPsi2SMassPtRapCut2 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecPsi2SPtRapCut2 = (TH2D*) histRecPsi2SMassPtRapCut2 -> Projection(2, 1, "histRecPsi2SPtRapCut2");
    TH1D *histRecPsi2SPtCut2 = (TH1D*) histRecPsi2SMassPtRapCut2 -> Projection(1, "histRecPsi2SPtCut2");
    TH1D *histRecPsi2SRapCut2 = (TH1D*) histRecPsi2SMassPtRapCut2 -> Projection(2, "histRecPsi2SRapCut2");

    THnSparseD *histRecPsi2SMassPtRapCut3 = (THnSparseD*) listRecPsi2SCut3 -> FindObject("Mass_Pt_Rapidity");
    histRecPsi2SMassPtRapCut3 -> GetAxis(1) -> SetRangeUser(0., 20.);
    histRecPsi2SMassPtRapCut3 -> GetAxis(2) -> SetRangeUser(2.5, 4);
    TH2D *histRecPsi2SPtRapCut3 = (TH2D*) histRecPsi2SMassPtRapCut3 -> Projection(2, 1, "histRecPsi2SPtRapCut3");
    TH1D *histRecPsi2SPtCut3 = (TH1D*) histRecPsi2SMassPtRapCut3 -> Projection(1, "histRecPsi2SPtCut3");
    TH1D *histRecPsi2SRapCut3 = (TH1D*) histRecPsi2SMassPtRapCut3 -> Projection(2, "histRecPsi2SRapCut3");

    THnSparseD *histRecPsi2SMassPtRapCut4 = (THnSparseD*) listRecPsi2SCut4 -> FindObject("Mass_Pt_Rapidity");
    histRecPsi2SMassPtRapCut4 -> GetAxis(1) -> SetRangeUser(0., 20.);
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
    /*
    TH1D *histGenJpsiPtRebin = (TH1D*) histGenJpsiPt -> Rebin(18, "histGenJpsiPtRebin", ptBinsJpsiRun2); 
    TH1D *histRecJpsiPtCut1Rebin = (TH1D*) histRecJpsiPtCut1 -> Rebin(18, "histRecJpsiPtCut1Rebin", ptBinsJpsiRun2); 
    SetHistogram(histRecJpsiPtCut1Rebin, 633, 1, 20, 0.8);
    TH1D *histRecJpsiPtCut2Rebin = (TH1D*) histRecJpsiPtCut2 -> Rebin(18, "histRecJpsiPtCut2Rebin", ptBinsJpsiRun2); 
    SetHistogram(histRecJpsiPtCut2Rebin, 862, 1, 20, 0.8);
    TH1D *histRecJpsiPtCut3Rebin = (TH1D*) histRecJpsiPtCut3 -> Rebin(18, "histRecJpsiPtCut3Rebin", ptBinsJpsiRun2); 
    SetHistogram(histRecJpsiPtCut3Rebin, 417, 1, 20, 0.8);
    TH1D *histRecJpsiPtCut4Rebin = (TH1D*) histRecJpsiPtCut4 -> Rebin(18, "histRecJpsiPtCut4Rebin", ptBinsJpsiRun2); 
    SetHistogram(histRecJpsiPtCut4Rebin, 880, 1, 20, 0.8);

    TH1D *histAxeJpsiPtCut1 = new TH1D("histAxeJpsiPtCut1", "", 18, ptBinsJpsiRun2);
    histAxeJpsiPtCut1 -> Divide(histRecJpsiPtCut1Rebin, histGenJpsiPtRebin, 1, 1, "B");
    SetHistogram(histAxeJpsiPtCut1, 633, 1, 20, 0.8);

    TH1D *histAxeJpsiPtCut2 = new TH1D("histAxeJpsiPtCut2", "", 18, ptBinsJpsiRun2);
    histAxeJpsiPtCut2 -> Divide(histRecJpsiPtCut2Rebin, histGenJpsiPtRebin, 1, 1, "B");
    SetHistogram(histAxeJpsiPtCut2, 862, 1, 20, 0.8);

    TH1D *histAxeJpsiPtCut3 = new TH1D("histAxeJpsiPtCut3", "", 18, ptBinsJpsiRun2);
    histAxeJpsiPtCut3 -> Divide(histRecJpsiPtCut3Rebin, histGenJpsiPtRebin, 1, 1, "B");
    SetHistogram(histAxeJpsiPtCut3, 417, 1, 20, 0.8);

    TH1D *histAxeJpsiPtCut4 = new TH1D("histAxeJpsiPtCut4", "", 18, ptBinsJpsiRun2);
    histAxeJpsiPtCut4 -> Divide(histRecJpsiPtCut4Rebin, histGenJpsiPtRebin, 1, 1, "B");
    SetHistogram(histAxeJpsiPtCut4, 880, 1, 20, 0.8);
    */

    TH1D *histGenJpsiPtRebin = (TH1D*) histGenJpsiPt -> Rebin(8, "histGenJpsiPtRebin", ptBinsRun3); 
    TH1D *histRecJpsiPtCut1Rebin = (TH1D*) histRecJpsiPtCut1 -> Rebin(8, "histRecJpsiPtCut1Rebin", ptBinsRun3); 
    SetHistogram(histRecJpsiPtCut1Rebin, 633, 1, 20, 0.8);
    TH1D *histRecJpsiPtCut2Rebin = (TH1D*) histRecJpsiPtCut2 -> Rebin(8, "histRecJpsiPtCut2Rebin", ptBinsRun3); 
    SetHistogram(histRecJpsiPtCut2Rebin, 862, 1, 20, 0.8);
    TH1D *histRecJpsiPtCut3Rebin = (TH1D*) histRecJpsiPtCut3 -> Rebin(8, "histRecJpsiPtCut3Rebin", ptBinsRun3); 
    SetHistogram(histRecJpsiPtCut3Rebin, 417, 1, 20, 0.8);
    TH1D *histRecJpsiPtCut4Rebin = (TH1D*) histRecJpsiPtCut4 -> Rebin(8, "histRecJpsiPtCut4Rebin", ptBinsRun3); 
    SetHistogram(histRecJpsiPtCut4Rebin, 880, 1, 20, 0.8);

    TH1D *histAxeJpsiPtCut1 = new TH1D("histAxeJpsiPtCut1", "", 8, ptBinsRun3);
    histAxeJpsiPtCut1 -> Divide(histRecJpsiPtCut1Rebin, histGenJpsiPtRebin, 1, 1, "B");
    SetHistogram(histAxeJpsiPtCut1, 633, 1, 20, 0.8);

    TH1D *histAxeJpsiPtCut2 = new TH1D("histAxeJpsiPtCut2", "", 8, ptBinsRun3);
    histAxeJpsiPtCut2 -> Divide(histRecJpsiPtCut2Rebin, histGenJpsiPtRebin, 1, 1, "B");
    SetHistogram(histAxeJpsiPtCut2, 862, 1, 20, 0.8);

    TH1D *histAxeJpsiPtCut3 = new TH1D("histAxeJpsiPtCut3", "", 8, ptBinsRun3);
    histAxeJpsiPtCut3 -> Divide(histRecJpsiPtCut3Rebin, histGenJpsiPtRebin, 1, 1, "B");
    SetHistogram(histAxeJpsiPtCut3, 417, 1, 20, 0.8);

    TH1D *histAxeJpsiPtCut4 = new TH1D("histAxeJpsiPtCut4", "", 8, ptBinsRun3);
    histAxeJpsiPtCut4 -> Divide(histRecJpsiPtCut4Rebin, histGenJpsiPtRebin, 1, 1, "B");
    SetHistogram(histAxeJpsiPtCut4, 880, 1, 20, 0.8);

    TH1D *histGenJpsiRapRebin = (TH1D*) histGenJpsiRap -> Rebin(6, "histGenJpsiRapRebin", rapBinsJpsiRun2); 
    TH1D *histRecJpsiRapCut1Rebin = (TH1D*) histRecJpsiRapCut1 -> Rebin(6, "histRecJpsiRapCut1Rebin", rapBinsJpsiRun2); 
    SetHistogram(histRecJpsiRapCut1Rebin, 633, 1, 20, 0.8);
    TH1D *histRecJpsiRapCut2Rebin = (TH1D*) histRecJpsiRapCut2 -> Rebin(6, "histRecJpsiRapCut2Rebin", rapBinsJpsiRun2); 
    SetHistogram(histRecJpsiRapCut2Rebin, 862, 1, 20, 0.8);
    TH1D *histRecJpsiRapCut3Rebin = (TH1D*) histRecJpsiRapCut3 -> Rebin(6, "histRecJpsiRapCut3Rebin", rapBinsJpsiRun2); 
    SetHistogram(histRecJpsiRapCut3Rebin, 417, 1, 20, 0.8);
    TH1D *histRecJpsiRapCut4Rebin = (TH1D*) histRecJpsiRapCut4 -> Rebin(6, "histRecJpsiRapCut4Rebin", rapBinsJpsiRun2); 
    SetHistogram(histRecJpsiRapCut4Rebin, 880, 1, 20, 0.8);

    TH1D *histAxeJpsiRapCut1 = new TH1D("histAxeJpsiRapCut1", "", 6, rapBinsJpsiRun2);
    histAxeJpsiRapCut1 -> Divide(histRecJpsiRapCut1Rebin, histGenJpsiRapRebin, 1, 1, "B");
    SetHistogram(histAxeJpsiRapCut1, 633, 1, 20, 0.8);

    TH1D *histAxeJpsiRapCut2 = new TH1D("histAxeJpsiRapCut2", "", 6, rapBinsJpsiRun2);
    histAxeJpsiRapCut2 -> Divide(histRecJpsiRapCut2Rebin, histGenJpsiRapRebin, 1, 1, "B");
    SetHistogram(histAxeJpsiRapCut2, 862, 1, 20, 0.8);

    TH1D *histAxeJpsiRapCut3 = new TH1D("histAxeJpsiRapCut3", "", 6, rapBinsJpsiRun2);
    histAxeJpsiRapCut3 -> Divide(histRecJpsiRapCut3Rebin, histGenJpsiRapRebin, 1, 1, "B");
    SetHistogram(histAxeJpsiRapCut3, 417, 1, 20, 0.8);

    TH1D *histAxeJpsiRapCut4 = new TH1D("histAxeJpsiRapCut4", "", 6, rapBinsJpsiRun2);
    histAxeJpsiRapCut4 -> Divide(histRecJpsiRapCut4Rebin, histGenJpsiRapRebin, 1, 1, "B");
    SetHistogram(histAxeJpsiRapCut4, 880, 1, 20, 0.8);

    /*
    TH1D *histGenPsi2SPtRebin = (TH1D*) histGenPsi2SPt -> Rebin(18, "histGenPsi2SPtRebin", ptBinsJpsiRun2); 
    TH1D *histRecPsi2SPtCut1Rebin = (TH1D*) histRecPsi2SPtCut1 -> Rebin(18, "histRecPsi2SPtCut1Rebin", ptBinsJpsiRun2); 
    SetHistogram(histRecPsi2SPtCut1Rebin, 633, 1, 20, 0.8);
    TH1D *histRecPsi2SPtCut2Rebin = (TH1D*) histRecPsi2SPtCut2 -> Rebin(18, "histRecPsi2SPtCut2Rebin", ptBinsJpsiRun2); 
    SetHistogram(histRecPsi2SPtCut2Rebin, 862, 1, 20, 0.8);
    TH1D *histRecPsi2SPtCut3Rebin = (TH1D*) histRecPsi2SPtCut3 -> Rebin(18, "histRecPsi2SPtCut3Rebin", ptBinsJpsiRun2); 
    SetHistogram(histRecPsi2SPtCut3Rebin, 417, 1, 20, 0.8);
    TH1D *histRecPsi2SPtCut4Rebin = (TH1D*) histRecPsi2SPtCut4 -> Rebin(18, "histRecPsi2SPtCut4Rebin", ptBinsJpsiRun2); 
    SetHistogram(histRecPsi2SPtCut4Rebin, 880, 1, 20, 0.8);

    TH1D *histAxePsi2SPtCut1 = new TH1D("histAxePsi2SPtCut1", "", 18, ptBinsJpsiRun2);
    histAxePsi2SPtCut1 -> Divide(histRecPsi2SPtCut1Rebin, histGenPsi2SPtRebin, 1, 1, "B");
    SetHistogram(histAxePsi2SPtCut1, 633, 1, 20, 0.8);

    TH1D *histAxePsi2SPtCut2 = new TH1D("histAxePsi2SPtCut2", "", 18, ptBinsJpsiRun2);
    histAxePsi2SPtCut2 -> Divide(histRecPsi2SPtCut2Rebin, histGenPsi2SPtRebin, 1, 1, "B");
    SetHistogram(histAxePsi2SPtCut2, 862, 1, 20, 0.8);

    TH1D *histAxePsi2SPtCut3 = new TH1D("histAxePsi2SPtCut3", "", 18, ptBinsJpsiRun2);
    histAxePsi2SPtCut3 -> Divide(histRecPsi2SPtCut3Rebin, histGenPsi2SPtRebin, 1, 1, "B");
    SetHistogram(histAxePsi2SPtCut3, 417, 1, 20, 0.8);

    TH1D *histAxePsi2SPtCut4 = new TH1D("histAxePsi2SPtCut4", "", 18, ptBinsJpsiRun2);
    histAxePsi2SPtCut4 -> Divide(histRecPsi2SPtCut4Rebin, histGenPsi2SPtRebin, 1, 1, "B");
    SetHistogram(histAxePsi2SPtCut4, 880, 1, 20, 0.8);
    */

    TH1D *histGenPsi2SPtRebin = (TH1D*) histGenPsi2SPt -> Rebin(8, "histGenPsi2SPtRebin", ptBinsRun3); 
    TH1D *histRecPsi2SPtCut1Rebin = (TH1D*) histRecPsi2SPtCut1 -> Rebin(8, "histRecPsi2SPtCut1Rebin", ptBinsRun3); 
    SetHistogram(histRecPsi2SPtCut1Rebin, 633, 1, 20, 0.8);
    TH1D *histRecPsi2SPtCut2Rebin = (TH1D*) histRecPsi2SPtCut2 -> Rebin(8, "histRecPsi2SPtCut2Rebin", ptBinsRun3); 
    SetHistogram(histRecPsi2SPtCut2Rebin, 862, 1, 20, 0.8);
    TH1D *histRecPsi2SPtCut3Rebin = (TH1D*) histRecPsi2SPtCut3 -> Rebin(8, "histRecPsi2SPtCut3Rebin", ptBinsRun3); 
    SetHistogram(histRecPsi2SPtCut3Rebin, 417, 1, 20, 0.8);
    TH1D *histRecPsi2SPtCut4Rebin = (TH1D*) histRecPsi2SPtCut4 -> Rebin(8, "histRecPsi2SPtCut4Rebin", ptBinsRun3); 
    SetHistogram(histRecPsi2SPtCut4Rebin, 880, 1, 20, 0.8);

    TH1D *histAxePsi2SPtCut1 = new TH1D("histAxePsi2SPtCut1", "", 8, ptBinsRun3);
    histAxePsi2SPtCut1 -> Divide(histRecPsi2SPtCut1Rebin, histGenPsi2SPtRebin, 1, 1, "B");
    SetHistogram(histAxePsi2SPtCut1, 633, 1, 20, 0.8);

    TH1D *histAxePsi2SPtCut2 = new TH1D("histAxePsi2SPtCut2", "", 8, ptBinsRun3);
    histAxePsi2SPtCut2 -> Divide(histRecPsi2SPtCut2Rebin, histGenPsi2SPtRebin, 1, 1, "B");
    SetHistogram(histAxePsi2SPtCut2, 862, 1, 20, 0.8);

    TH1D *histAxePsi2SPtCut3 = new TH1D("histAxePsi2SPtCut3", "", 8, ptBinsRun3);
    histAxePsi2SPtCut3 -> Divide(histRecPsi2SPtCut3Rebin, histGenPsi2SPtRebin, 1, 1, "B");
    SetHistogram(histAxePsi2SPtCut3, 417, 1, 20, 0.8);

    TH1D *histAxePsi2SPtCut4 = new TH1D("histAxePsi2SPtCut4", "", 8, ptBinsRun3);
    histAxePsi2SPtCut4 -> Divide(histRecPsi2SPtCut4Rebin, histGenPsi2SPtRebin, 1, 1, "B");
    SetHistogram(histAxePsi2SPtCut4, 880, 1, 20, 0.8);

    TH1D *histGenPsi2SRapRebin = (TH1D*) histGenPsi2SRap -> Rebin(6, "histGenPsi2SRapRebin", rapBinsJpsiRun2); 
    TH1D *histRecPsi2SRapCut1Rebin = (TH1D*) histRecPsi2SRapCut1 -> Rebin(6, "histRecPsi2SRapCut1Rebin", rapBinsJpsiRun2); 
    SetHistogram(histRecPsi2SRapCut1Rebin, 633, 1, 20, 0.8);
    TH1D *histRecPsi2SRapCut2Rebin = (TH1D*) histRecPsi2SRapCut2 -> Rebin(6, "histRecPsi2SRapCut2Rebin", rapBinsJpsiRun2); 
    SetHistogram(histRecPsi2SRapCut2Rebin, 862, 1, 20, 0.8);
    TH1D *histRecPsi2SRapCut3Rebin = (TH1D*) histRecPsi2SRapCut3 -> Rebin(6, "histRecPsi2SRapCut3Rebin", rapBinsJpsiRun2); 
    SetHistogram(histRecPsi2SRapCut3Rebin, 417, 1, 20, 0.8);
    TH1D *histRecPsi2SRapCut4Rebin = (TH1D*) histRecPsi2SRapCut4 -> Rebin(6, "histRecPsi2SRapCut4Rebin", rapBinsJpsiRun2); 
    SetHistogram(histRecPsi2SRapCut4Rebin, 880, 1, 20, 0.8);

    TH1D *histAxePsi2SRapCut1 = new TH1D("histAxePsi2SRapCut1", "", 6, rapBinsPsi2SRun2);
    histAxePsi2SRapCut1 -> Divide(histRecPsi2SRapCut1Rebin, histGenPsi2SRapRebin, 1, 1, "B");
    SetHistogram(histAxePsi2SRapCut1, 633, 1, 20, 0.8);

    TH1D *histAxePsi2SRapCut2 = new TH1D("histAxePsi2SRapCut2", "", 6, rapBinsPsi2SRun2);
    histAxePsi2SRapCut2 -> Divide(histRecPsi2SRapCut2Rebin, histGenPsi2SRapRebin, 1, 1, "B");
    SetHistogram(histAxePsi2SRapCut2, 862, 1, 20, 0.8);

    TH1D *histAxePsi2SRapCut3 = new TH1D("histAxePsi2SRapCut3", "", 6, rapBinsPsi2SRun2);
    histAxePsi2SRapCut3 -> Divide(histRecPsi2SRapCut3Rebin, histGenPsi2SRapRebin, 1, 1, "B");
    SetHistogram(histAxePsi2SRapCut3, 417, 1, 20, 0.8);

    TH1D *histAxePsi2SRapCut4 = new TH1D("histAxePsi2SRapCut4", "", 6, rapBinsPsi2SRun2);
    histAxePsi2SRapCut4 -> Divide(histRecPsi2SRapCut4Rebin, histGenPsi2SRapRebin, 1, 1, "B");
    SetHistogram(histAxePsi2SRapCut4, 880, 1, 20, 0.8);


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
    histAxeJpsiPtRun3Prel -> Draw("EP SAME");

    canvasJpsi -> cd(4);
    histAxeJpsiRapCut1 -> GetXaxis() -> SetTitle("#it{y}");
    histAxeJpsiRapCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histAxeJpsiRapCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histAxeJpsiRapCut1 -> Draw("EP");
    histAxeJpsiRapCut2 -> Draw("EP SAME");
    histAxeJpsiRapCut3 -> Draw("EP SAME");
    histAxeJpsiRapCut4 -> Draw("EP SAME");
    histAxeJpsiRapRun2 -> Draw("EP SAME");

    TCanvas *canvasAxeJpsiPt = new TCanvas("canvasAxeJpsiPt", "", 800, 600);
    histAxeJpsiPtCut1 -> SetTitle("");
    histAxeJpsiPtCut1 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histAxeJpsiPtCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon_{J/#psi}");
    histAxeJpsiPtCut1 -> GetXaxis() -> SetRangeUser(0, 20.);
    histAxeJpsiPtCut1 -> GetYaxis() -> SetRangeUser(0, 1.2);
    histAxeJpsiPtCut1 -> Draw("EP");
    histAxeJpsiPtCut2 -> Draw("EP SAME");
    histAxeJpsiPtCut3 -> Draw("EP SAME");
    histAxeJpsiPtCut4 -> Draw("EP SAME");
    histAxeJpsiPtRun2 -> Draw("EP SAME");
    histAxeJpsiPtRun3Prel -> Draw("EP SAME");

    TLegend *legendAxeJpsiPt = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendAxeJpsiPt);
    legendAxeJpsiPt -> AddEntry(histAxeJpsiPtRun2, "Run 2", "PL");
    legendAxeJpsiPt -> AddEntry(histAxeJpsiPtRun3Prel, "Run 3 - Preliminary", "PL");
    legendAxeJpsiPt -> AddEntry(histAxeJpsiPtCut1, "Run 3 - matchedMchMid", "PL");
    legendAxeJpsiPt -> AddEntry(histAxeJpsiPtCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendAxeJpsiPt -> AddEntry(histAxeJpsiPtCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendAxeJpsiPt -> AddEntry(histAxeJpsiPtCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendAxeJpsiPt -> Draw("SAME");


    TCanvas *canvasAxeJpsiRap = new TCanvas("canvasAxeJpsiRap", "", 800, 600);
    histAxeJpsiRapCut1 -> SetTitle("");
    histAxeJpsiRapCut1 -> GetXaxis() -> SetTitle("#it{y}");
    histAxeJpsiRapCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon_{J/#psi}");
    histAxeJpsiRapCut1 -> GetYaxis() -> SetRangeUser(0, 1.2);
    histAxeJpsiRapCut1 -> Draw("EP");
    histAxeJpsiRapCut2 -> Draw("EP SAME");
    histAxeJpsiRapCut3 -> Draw("EP SAME");
    histAxeJpsiRapCut4 -> Draw("EP SAME");
    histAxeJpsiRapRun2 -> Draw("EP SAME");
    histAxeJpsiRapRun3Prel -> Draw("EP SAME");

    TLegend *legendAxeJpsiRap = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendAxeJpsiRap);
    legendAxeJpsiRap -> AddEntry(histAxeJpsiRapRun2, "Run 2", "PL");
    legendAxeJpsiRap -> AddEntry(histAxeJpsiRapRun3Prel, "Run 3 - Preliminary", "PL");
    legendAxeJpsiRap -> AddEntry(histAxeJpsiRapCut1, "Run 3 - matchedMchMid", "PL");
    legendAxeJpsiRap -> AddEntry(histAxeJpsiRapCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendAxeJpsiRap -> AddEntry(histAxeJpsiRapCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendAxeJpsiRap -> AddEntry(histAxeJpsiRapCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendAxeJpsiRap -> Draw("SAME");


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
    histAxePsi2SPtRun3Prel -> Draw("EP SAME");

    canvasPsi2S -> cd(4);
    histAxePsi2SRapCut1 -> GetXaxis() -> SetTitle("#it{y}");
    histAxePsi2SRapCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon");
    histAxePsi2SRapCut1 -> GetYaxis() -> SetRangeUser(0, 1);
    histAxePsi2SRapCut1 -> Draw("EP");
    histAxePsi2SRapCut2 -> Draw("EP SAME");
    histAxePsi2SRapCut3 -> Draw("EP SAME");
    histAxePsi2SRapCut4 -> Draw("EP SAME");
    histAxePsi2SRapRun2 -> Draw("EP SAME");


    TCanvas *canvasAxePsi2SPt = new TCanvas("canvasAxePsi2SPt", "", 800, 600);
    histAxePsi2SPtCut1 -> SetTitle("");
    histAxePsi2SPtCut1 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histAxePsi2SPtCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon_{#psi(2S)}");
    histAxePsi2SPtCut1 -> GetXaxis() -> SetRangeUser(0, 20.);
    histAxePsi2SPtCut1 -> GetYaxis() -> SetRangeUser(0, 1.2);
    histAxePsi2SPtCut1 -> Draw("EP");
    histAxePsi2SPtCut2 -> Draw("EP SAME");
    histAxePsi2SPtCut3 -> Draw("EP SAME");
    histAxePsi2SPtCut4 -> Draw("EP SAME");
    histAxePsi2SPtRun2 -> Draw("EP SAME");
    histAxePsi2SPtRun3Prel -> Draw("EP SAME");

    TLegend *legendAxePsi2SPt = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendAxePsi2SPt);
    legendAxePsi2SPt -> AddEntry(histAxePsi2SPtRun2, "Run 2", "PL");
    legendAxePsi2SPt -> AddEntry(histAxePsi2SPtRun3Prel, "Run 3 - Preliminary", "PL");
    legendAxePsi2SPt -> AddEntry(histAxePsi2SPtCut1, "Run 3 - matchedMchMid", "PL");
    legendAxePsi2SPt -> AddEntry(histAxePsi2SPtCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendAxePsi2SPt -> AddEntry(histAxePsi2SPtCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendAxePsi2SPt -> AddEntry(histAxePsi2SPtCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendAxePsi2SPt -> Draw("SAME");


    TCanvas *canvasAxePsi2SRap = new TCanvas("canvasAxePsi2SRap", "", 800, 600);
    histAxePsi2SRapCut1 -> SetTitle("");
    histAxePsi2SRapCut1 -> GetXaxis() -> SetTitle("#it{y}");
    histAxePsi2SRapCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon_{#psi(2S)}");
    histAxePsi2SRapCut1 -> GetYaxis() -> SetRangeUser(0, 1.2);
    histAxePsi2SRapCut1 -> Draw("EP");
    histAxePsi2SRapCut2 -> Draw("EP SAME");
    histAxePsi2SRapCut3 -> Draw("EP SAME");
    histAxePsi2SRapCut4 -> Draw("EP SAME");
    histAxePsi2SRapRun2 -> Draw("EP SAME");
    histAxePsi2SRapRun3Prel -> Draw("EP SAME");

    TLegend *legendAxePsi2SRap = new TLegend(0.20, 0.65, 0.60, 0.89);
    SetLegend(legendAxePsi2SRap);
    legendAxePsi2SRap -> AddEntry(histAxePsi2SRapRun2, "Run 2", "PL");
    legendAxePsi2SRap -> AddEntry(histAxePsi2SRapRun3Prel, "Run 3 - Preliminary", "PL");
    legendAxePsi2SRap -> AddEntry(histAxePsi2SRapCut1, "Run 3 - matchedMchMid", "PL");
    legendAxePsi2SRap -> AddEntry(histAxePsi2SRapCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendAxePsi2SRap -> AddEntry(histAxePsi2SRapCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendAxePsi2SRap -> AddEntry(histAxePsi2SRapCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendAxePsi2SRap -> Draw("SAME");

    // Axe ratio
    TH1D *histRatioAxePtCut1 = (TH1D*) histAxePsi2SPtCut1 -> Clone("histRatioAxePtCut1");
    histRatioAxePtCut1 -> Divide(histAxeJpsiPtCut1);
    SetHistogram(histRatioAxePtCut1, 633, 1, 20, 0.8);

    TH1D *histRatioAxePtCut2 = (TH1D*) histAxePsi2SPtCut2 -> Clone("histRatioAxePtCut2");
    histRatioAxePtCut2 -> Divide(histAxeJpsiPtCut2);
    SetHistogram(histRatioAxePtCut2, 862, 1, 20, 0.8);

    TH1D *histRatioAxePtCut3 = (TH1D*) histAxePsi2SPtCut3 -> Clone("histRatioAxePtCut3");
    histRatioAxePtCut3 -> Divide(histAxeJpsiPtCut3);
    SetHistogram(histRatioAxePtCut3, 417, 1, 20, 0.8);

    TH1D *histRatioAxePtCut4 = (TH1D*) histAxePsi2SPtCut4 -> Clone("histRatioAxePtCut4");
    histRatioAxePtCut4 -> Divide(histAxeJpsiPtCut4);
    SetHistogram(histRatioAxePtCut4, 880, 1, 20, 0.8);


    TCanvas *canvasRatioAxePt = new TCanvas("canvasRatioAxePt", "", 800, 600);
    histRatioAxePtCut1 -> SetTitle("");
    histRatioAxePtCut1 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histRatioAxePtCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon_{#psi(2S)} / A#times#varepsilon_{J/#psi}");
    histRatioAxePtCut1 -> GetXaxis() -> SetRangeUser(0, 20.);
    histRatioAxePtCut1 -> GetYaxis() -> SetRangeUser(0.4, 1.6);
    histRatioAxePtCut1 -> Draw("EP");
    histRatioAxePtCut2 -> Draw("EP SAME");
    histRatioAxePtCut3 -> Draw("EP SAME");
    histRatioAxePtCut4 -> Draw("EP SAME");
    histRatioAxePtRun2 -> Draw("EP SAME");
    histRatioAxePtRun3Prel -> Draw("EP SAME");

    TLegend *legendRatioAxePt = new TLegend(0.40, 0.65, 0.80, 0.89);
    SetLegend(legendRatioAxePt);
    legendRatioAxePt -> AddEntry(histRatioAxePtRun2, "Run 2", "PL");
    legendRatioAxePt -> AddEntry(histRatioAxePtRun3Prel, "Run 3 - Preliminary", "PL");
    legendRatioAxePt -> AddEntry(histRatioAxePtCut1, "Run 3 - matchedMchMid", "PL");
    legendRatioAxePt -> AddEntry(histRatioAxePtCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendRatioAxePt -> AddEntry(histRatioAxePtCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendRatioAxePt -> AddEntry(histRatioAxePtCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendRatioAxePt -> Draw("SAME");

    TH1D *histRatioAxeRapCut1 = (TH1D*) histAxePsi2SRapCut1 -> Clone("histRatioAxeRapCut1");
    histRatioAxeRapCut1 -> Divide(histAxeJpsiRapCut1);
    SetHistogram(histRatioAxeRapCut1, 633, 1, 20, 0.8);

    TH1D *histRatioAxeRapCut2 = (TH1D*) histAxePsi2SRapCut2 -> Clone("histRatioAxeRapCut2");
    histRatioAxeRapCut2 -> Divide(histAxeJpsiRapCut2);
    SetHistogram(histRatioAxeRapCut2, 862, 1, 20, 0.8);

    TH1D *histRatioAxeRapCut3 = (TH1D*) histAxePsi2SRapCut3 -> Clone("histRatioAxeRapCut3");
    histRatioAxeRapCut3 -> Divide(histAxeJpsiRapCut3);
    SetHistogram(histRatioAxeRapCut3, 417, 1, 20, 0.8);

    TH1D *histRatioAxeRapCut4 = (TH1D*) histAxePsi2SRapCut4 -> Clone("histRatioAxeRapCut4");
    histRatioAxeRapCut4 -> Divide(histAxeJpsiRapCut4);
    SetHistogram(histRatioAxeRapCut4, 880, 1, 20, 0.8);

    TCanvas *canvasRatioRap = new TCanvas("canvasRatioRap", "", 800, 600);
    histRatioAxeRapCut1 -> SetTitle("");
    histRatioAxeRapCut1 -> GetXaxis() -> SetTitle("#it{y}");
    histRatioAxeRapCut1 -> GetYaxis() -> SetTitle("A#times#varepsilon_{#psi(2S)} / A#times#varepsilon_{J/#psi}");
    histRatioAxeRapCut1 -> GetYaxis() -> SetRangeUser(0, 2);
    histRatioAxeRapCut1 -> Draw("EP");
    histRatioAxeRapCut2 -> Draw("EP SAME");
    histRatioAxeRapCut3 -> Draw("EP SAME");
    histRatioAxeRapCut4 -> Draw("EP SAME");
    histRatioAxeRapRun2 -> Draw("EP SAME");
    histRatioAxeRapRun3Prel -> Draw("EP SAME");

    TLegend *legendRatioAxeRap = new TLegend(0.40, 0.25, 0.80, 0.49);
    SetLegend(legendRatioAxeRap);
    legendRatioAxeRap -> AddEntry(histRatioAxeRapRun2, "Run 2", "PL");
    legendRatioAxeRap -> AddEntry(histRatioAxeRapRun3Prel, "Run 3 - Preliminary", "PL");
    legendRatioAxeRap -> AddEntry(histRatioAxeRapCut1, "Run 3 - matchedMchMid", "PL");
    legendRatioAxeRap -> AddEntry(histRatioAxeRapCut2, "Run 3 - muonLowPt10SigmaPDCA", "PL");
    legendRatioAxeRap -> AddEntry(histRatioAxeRapCut3, "Run 3 - muonLowPt210SigmaPDCA", "PL");
    legendRatioAxeRap -> AddEntry(histRatioAxeRapCut4, "Run 3 - muonLowPt510SigmaPDCA", "PL");
    legendRatioAxeRap -> Draw("SAME");

    canvasAxeJpsiPt -> SaveAs(Form("plots/Axe_Jpsi_pt_%s_%s.pdf", productionName.c_str(), associationType.c_str()));
    canvasAxeJpsiRap -> SaveAs(Form("plots/Axe_Jpsi_rap_%s_%s.pdf", productionName.c_str(), associationType.c_str()));
    canvasAxePsi2SPt -> SaveAs(Form("plots/Axe_Psi2S_pt_%s_%s.pdf", productionName.c_str(), associationType.c_str()));
    canvasAxePsi2SRap -> SaveAs(Form("plots/Axe_Psi2S_rap_%s_%s.pdf", productionName.c_str(), associationType.c_str()));
    canvasRatioAxePt -> SaveAs(Form("plots/Axe_ratio_pt_%s_%s.pdf", productionName.c_str(), associationType.c_str()));
    canvasRatioRap -> SaveAs(Form("plots/Axe_ratio_rap_%s_%s.pdf", productionName.c_str(), associationType.c_str()));


    TFile *fOut = new TFile(Form("Axe_%s_%s.root", productionName.c_str(), associationType.c_str()), "RECREATE");
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
    histRatioAxePtCut1 -> Write();
    histRatioAxePtCut2 -> Write();
    histRatioAxePtCut3 -> Write();
    histRatioAxePtCut4 -> Write();
    histRatioAxeRapCut1 -> Write();
    histRatioAxeRapCut2 -> Write();
    histRatioAxeRapCut3 -> Write();
    histRatioAxeRapCut4 -> Write();
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
////////////////////////////////////////////////////////////////////////////////
void SetHistogram(TH1D *hist, int color, int lineWidth, int markerStyle, double markerSize) {
    hist -> SetLineColor(color);
    hist -> SetLineWidth(lineWidth);
    hist -> SetMarkerColor(color);
    hist -> SetMarkerStyle(markerStyle);
    hist -> SetMarkerSize(markerSize);
}