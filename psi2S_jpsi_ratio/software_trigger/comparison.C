void LoadStyle();
void SetLegend(TLegend *);
void SetHistogram(TH1D *, int , int , int , double );

void comparison() {
    LoadStyle();

    TFile *fIn_LHC22Pass7 = new TFile("/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data/LHC22_pass7_skimmed/AnalysisResults_LHC22p.root", "READ");
    TH1D *histEvAccCounter_LHC22Pass7 = (TH1D*) fIn_LHC22Pass7 -> Get("event-selection-task/hColCounterAcc");
    TList *listSE_LHC22Pass7 = (TList*) fIn_LHC22Pass7 -> Get("analysis-same-event-pairing/output");
    TList *listSEPM_LHC22Pass7 = (TList*) listSE_LHC22Pass7 -> FindObject("PairsMuonSEPM_matchedMchMid");
    TH2D *histMassPt_LHC22Pass7 = (TH2D*) listSEPM_LHC22Pass7 -> FindObject("Mass_Pt");
    TH1D *histMass_LHC22Pass7 = (TH1D*) histMassPt_LHC22Pass7 -> ProjectionX("histMass_LHC22Pass7");
    TH1D *histPt_LHC22Pass7 = (TH1D*) histMassPt_LHC22Pass7 -> ProjectionY("histPt_LHC22Pass7");
    histMass_LHC22Pass7 -> Scale(1. / histEvAccCounter_LHC22Pass7 -> Integral());
    histPt_LHC22Pass7 -> Scale(1. / histEvAccCounter_LHC22Pass7 -> Integral());

    TFile *fIn_LHC22Pass7_fDiMuon = new TFile("/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data/LHC22_pass7_skimmed/AnalysisResults_LHC22p_fDiMuon.root", "READ");
    TH1D *histEvAccCounter_LHC22Pass7_fDiMuon = (TH1D*) fIn_LHC22Pass7_fDiMuon -> Get("event-selection-task/hColCounterAcc");
    TList *listSE_LHC22Pass7_fDiMuon = (TList*) fIn_LHC22Pass7_fDiMuon -> Get("analysis-same-event-pairing/output");
    TList *listSEPM_LHC22Pass7_fDiMuon = (TList*) listSE_LHC22Pass7_fDiMuon -> FindObject("PairsMuonSEPM_matchedMchMid");
    TH2D *histMassPt_LHC22Pass7_fDiMuon = (TH2D*) listSEPM_LHC22Pass7_fDiMuon -> FindObject("Mass_Pt");
    TH1D *histMass_LHC22Pass7_fDiMuon = (TH1D*) histMassPt_LHC22Pass7_fDiMuon -> ProjectionX("histMass_LHC22Pass7_fDiMuon");
    TH1D *histPt_LHC22Pass7_fDiMuon = (TH1D*) histMassPt_LHC22Pass7_fDiMuon -> ProjectionY("histPt_LHC22Pass7_fDiMuon");
    histMass_LHC22Pass7_fDiMuon -> Scale(1. / histEvAccCounter_LHC22Pass7_fDiMuon -> Integral());
    histPt_LHC22Pass7_fDiMuon -> Scale(1. / histEvAccCounter_LHC22Pass7_fDiMuon -> Integral());

    TFile *fIn_LHC22Pass7_fSingleMuLow = new TFile("/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data/LHC22_pass7_skimmed/AnalysisResults_LHC22p_fSingleMuLow.root", "READ");
    TH1D *histEvAccCounter_LHC22Pass7_fSingleMuLow = (TH1D*) fIn_LHC22Pass7_fSingleMuLow -> Get("event-selection-task/hColCounterAcc");
    TList *listSE_LHC22Pass7_fSingleMuLow = (TList*) fIn_LHC22Pass7_fSingleMuLow -> Get("analysis-same-event-pairing/output");
    TList *listSEPM_LHC22Pass7_fSingleMuLow = (TList*) listSE_LHC22Pass7_fSingleMuLow -> FindObject("PairsMuonSEPM_matchedMchMid");
    TH2D *histMassPt_LHC22Pass7_fSingleMuLow = (TH2D*) listSEPM_LHC22Pass7_fSingleMuLow -> FindObject("Mass_Pt");
    TH1D *histMass_LHC22Pass7_fSingleMuLow = (TH1D*) histMassPt_LHC22Pass7_fSingleMuLow -> ProjectionX("histMass_LHC22Pass7_fSingleMuLow");
    TH1D *histPt_LHC22Pass7_fSingleMuLow = (TH1D*) histMassPt_LHC22Pass7_fSingleMuLow -> ProjectionY("histPt_LHC22Pass7_fSingleMuLow");
    histMass_LHC22Pass7_fSingleMuLow -> Scale(1. / histEvAccCounter_LHC22Pass7_fSingleMuLow -> Integral());
    histPt_LHC22Pass7_fSingleMuLow -> Scale(1. / histEvAccCounter_LHC22Pass7_fSingleMuLow -> Integral());

    //TFile *fInLHC22Pass7 = new TFile("histogram_LHC22_pass7.root", "READ");
    //TH1D *histMass_LHC22Pass7 = (TH1D*) fInLHC22Pass7 -> Get("Mass_Pt_px");
    //TH1D *histPt_LHC22Pass7 = (TH1D*) fInLHC22Pass7 -> Get("Mass_Pt_py");

    //TFile *fInLHC22Pass7_fDiMuon = new TFile("histogram_LHC22_pass7_fDiMuon.root", "READ");
    //TH1D *histMass_LHC22Pass7_fDiMuon = (TH1D*) fInLHC22Pass7_fDiMuon -> Get("Mass_Pt_px");
    //TH1D *histPt_LHC22Pass7_fDiMuon = (TH1D*) fInLHC22Pass7_fDiMuon -> Get("Mass_Pt_py");

    //TFile *fInLHC22Pass7_fSingleMuLow = new TFile("histogram_LHC22_pass7_fSingleMuLow.root", "READ");
    //TH1D *histMass_LHC22Pass7_fSingleMuLow = (TH1D*) fInLHC22Pass7_fSingleMuLow -> Get("Mass_Pt_px");
    //TH1D *histPt_LHC22Pass7_fSingleMuLow = (TH1D*) fInLHC22Pass7_fSingleMuLow -> Get("Mass_Pt_py");

    TFile *fIn_LHC23Pass4 = new TFile("/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data/LHC23_pass4_skimmed/AnalysisResults.root", "READ");
    TH1D *histEvAccCounter_LHC23Pass4 = (TH1D*) fIn_LHC23Pass4 -> Get("event-selection-task/hColCounterAcc");
    TList *listSE_LHC23Pass4 = (TList*) fIn_LHC23Pass4 -> Get("analysis-same-event-pairing/output");
    TList *listSEPM_LHC23Pass4 = (TList*) listSE_LHC23Pass4 -> FindObject("PairsMuonSEPM_matchedMchMid");
    TH2D *histMassPt_LHC23Pass4 = (TH2D*) listSEPM_LHC23Pass4 -> FindObject("Mass_Pt");
    TH1D *histMass_LHC23Pass4 = (TH1D*) histMassPt_LHC23Pass4 -> ProjectionX("histMass_LHC23Pass4");
    TH1D *histPt_LHC23Pass4 = (TH1D*) histMassPt_LHC23Pass4 -> ProjectionY("histPt_LHC23Pass4");
    histMass_LHC23Pass4 -> Scale(1. / histEvAccCounter_LHC23Pass4 -> Integral());
    histPt_LHC23Pass4 -> Scale(1. / histEvAccCounter_LHC23Pass4 -> Integral());

    TFile *fIn_LHC23Pass4_fDiMuon = new TFile("/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data/LHC23_pass4_skimmed/AnalysisResults_fDiMuon.root", "READ");
    TH1D *histEvAccCounter_LHC23Pass4_fDiMuon = (TH1D*) fIn_LHC23Pass4_fDiMuon -> Get("event-selection-task/hColCounterAcc");
    TList *listSE_LHC23Pass4_fDiMuon = (TList*) fIn_LHC23Pass4_fDiMuon -> Get("analysis-same-event-pairing/output");
    TList *listSEPM_LHC23Pass4_fDiMuon = (TList*) listSE_LHC23Pass4_fDiMuon -> FindObject("PairsMuonSEPM_matchedMchMid");
    TH2D *histMassPt_LHC23Pass4_fDiMuon = (TH2D*) listSEPM_LHC23Pass4_fDiMuon -> FindObject("Mass_Pt");
    TH1D *histMass_LHC23Pass4_fDiMuon = (TH1D*) histMassPt_LHC23Pass4_fDiMuon -> ProjectionX("histMass_LHC23Pass4_fDiMuon");
    TH1D *histPt_LHC23Pass4_fDiMuon = (TH1D*) histMassPt_LHC23Pass4_fDiMuon -> ProjectionY("histPt_LHC23Pass4_fDiMuon");
    histMass_LHC23Pass4_fDiMuon -> Scale(1. / histEvAccCounter_LHC23Pass4_fDiMuon -> Integral());
    histPt_LHC23Pass4_fDiMuon -> Scale(1. / histEvAccCounter_LHC23Pass4_fDiMuon -> Integral());

    TFile *fIn_LHC23Pass4_fSingleMuLow = new TFile("/Users/lucamicheletti/cernbox/JPSI/Psi2S_Jpsi_ratio_run3/data/LHC23_pass4_skimmed/AnalysisResults_fSingleMuLow.root", "READ");
    TH1D *histEvAccCounter_LHC23Pass4_fSingleMuLow = (TH1D*) fIn_LHC23Pass4_fSingleMuLow -> Get("event-selection-task/hColCounterAcc");
    TList *listSE_LHC23Pass4_fSingleMuLow = (TList*) fIn_LHC23Pass4_fSingleMuLow -> Get("analysis-same-event-pairing/output");
    TList *listSEPM_LHC23Pass4_fSingleMuLow = (TList*) listSE_LHC23Pass4_fSingleMuLow -> FindObject("PairsMuonSEPM_matchedMchMid");
    TH2D *histMassPt_LHC23Pass4_fSingleMuLow = (TH2D*) listSEPM_LHC23Pass4_fSingleMuLow -> FindObject("Mass_Pt");
    TH1D *histMass_LHC23Pass4_fSingleMuLow = (TH1D*) histMassPt_LHC23Pass4_fSingleMuLow -> ProjectionX("histMass_LHC23Pass4_fSingleMuLow");
    TH1D *histPt_LHC23Pass4_fSingleMuLow = (TH1D*) histMassPt_LHC23Pass4_fSingleMuLow -> ProjectionY("histPt_LHC23Pass4_fSingleMuLow");
    histMass_LHC23Pass4_fSingleMuLow -> Scale(1. / histEvAccCounter_LHC23Pass4_fSingleMuLow -> Integral());
    histPt_LHC23Pass4_fSingleMuLow -> Scale(1. / histEvAccCounter_LHC23Pass4_fSingleMuLow -> Integral());

    //TFile *fIn_LHC23Pass4 = new TFile("histogram_LHC23_pass4.root", "READ");
    //TH1D *histMass_LHC23Pass4 = (TH1D*) fIn_LHC23Pass4 -> Get("Mass_Pt_px");
    //TH1D *histPt_LHC23Pass4 = (TH1D*) fIn_LHC23Pass4 -> Get("Mass_Pt_py");

    //TFile *fIn_LHC23Pass4_fDiMuon = new TFile("histogram_LHC23_pass4_fDiMuon.root", "READ");
    //TH1D *histMass_LHC23Pass4_fDiMuon = (TH1D*) fIn_LHC23Pass4_fDiMuon -> Get("Mass_Pt_px");
    //TH1D *histPt_LHC23Pass4_fDiMuon = (TH1D*) fIn_LHC23Pass4_fDiMuon -> Get("Mass_Pt_py");

    //TFile *fIn_LHC23Pass4_fSingleMuLow = new TFile("histogram_LHC23_pass4_fSingleMuLow.root", "READ");
    //TH1D *histMass_LHC23Pass4_fSingleMuLow = (TH1D*) fIn_LHC23Pass4_fSingleMuLow -> Get("Mass_Pt_px");
    //TH1D *histPt_LHC23Pass4_fSingleMuLow = (TH1D*) fIn_LHC23Pass4_fSingleMuLow -> Get("Mass_Pt_py");

    histMass_LHC22Pass7 -> Rebin(5);
    histMass_LHC22Pass7_fDiMuon -> Rebin(5);
    histMass_LHC22Pass7_fSingleMuLow -> Rebin(5);

    histMass_LHC23Pass4 -> Rebin(5);
    histMass_LHC23Pass4_fDiMuon -> Rebin(5);
    histMass_LHC23Pass4_fSingleMuLow -> Rebin(5);

    SetHistogram(histMass_LHC22Pass7, 1, 1, 24, 0.8);
    SetHistogram(histMass_LHC22Pass7_fDiMuon, 633, 1, 0, 0.8);
    SetHistogram(histMass_LHC22Pass7_fSingleMuLow, 862, 1, 0, 0.8);

    SetHistogram(histPt_LHC22Pass7, 1, 1, 24, 0.8);
    SetHistogram(histPt_LHC22Pass7_fDiMuon, 633, 1, 0, 0.8);
    SetHistogram(histPt_LHC22Pass7_fSingleMuLow, 862, 1, 0, 0.8);

    SetHistogram(histMass_LHC23Pass4, 1, 1, 24, 0.8);
    SetHistogram(histMass_LHC23Pass4_fDiMuon, 633, 1, 0, 0.8);
    SetHistogram(histMass_LHC23Pass4_fSingleMuLow, 862, 1, 0, 0.8);

    SetHistogram(histPt_LHC23Pass4, 1, 1, 24, 0.8);
    SetHistogram(histPt_LHC23Pass4_fDiMuon, 633, 1, 0, 0.8);
    SetHistogram(histPt_LHC23Pass4_fSingleMuLow, 862, 1, 0, 0.8);

    //////////////////////
    // Mass dependence //
    /////////////////////
    TLegend *legendMass = new TLegend(0.60, 0.75, 0.85, 0.89);
    SetLegend(legendMass);
    legendMass -> AddEntry(histMass_LHC23Pass4, "All triggers", "PL");
    legendMass -> AddEntry(histMass_LHC23Pass4_fDiMuon, "#bf{fDiMuon} trigger", "L");
    legendMass -> AddEntry(histMass_LHC23Pass4_fSingleMuLow, "#bf{fSingleMuLow} trigger", "L");

    TCanvas *canvasMassLHC22Pass7 = new TCanvas("canvasMassLHC22Pass7", "", 800, 600);
    gPad -> SetLogy(true);

    histMass_LHC22Pass7 -> SetStats(false);
    histMass_LHC22Pass7_fDiMuon -> SetStats(false);
    histMass_LHC22Pass7_fSingleMuLow -> SetStats(false);

    histMass_LHC22Pass7 -> SetTitle("");
    histMass_LHC22Pass7 -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c}^{2})");
    histMass_LHC22Pass7 -> GetYaxis() -> SetTitle("Counts");
    histMass_LHC22Pass7 -> GetXaxis() -> SetRangeUser(0, 10);
    histMass_LHC22Pass7 -> Draw("EP");
    histMass_LHC22Pass7_fDiMuon -> Draw("H SAME");
    histMass_LHC22Pass7_fSingleMuLow -> Draw("H SAME");
    legendMass -> Draw("SAME");

    TCanvas *canvasMassLHC23Pass4 = new TCanvas("canvasMassLHC23Pass4", "", 800, 600);
    gPad -> SetLogy(true);

    histMass_LHC23Pass4 -> SetStats(false);
    histMass_LHC23Pass4_fDiMuon -> SetStats(false);
    histMass_LHC23Pass4_fSingleMuLow -> SetStats(false);

    histMass_LHC23Pass4 -> SetTitle("");
    histMass_LHC23Pass4 -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c}^{2})");
    histMass_LHC23Pass4 -> GetYaxis() -> SetTitle("Counts");
    histMass_LHC23Pass4 -> GetXaxis() -> SetRangeUser(0, 10);
    histMass_LHC23Pass4 -> Draw("EP");
    histMass_LHC23Pass4_fDiMuon -> Draw("H SAME");
    histMass_LHC23Pass4_fSingleMuLow -> Draw("H SAME");
    legendMass -> Draw("SAME");

    TH1D *histMass_LHC22Pass7Ratio_fDiMuon = (TH1D*) histMass_LHC22Pass7_fDiMuon -> Clone("histMass_LHC22Pass7Ratio_fDiMuon");   
    histMass_LHC22Pass7Ratio_fDiMuon -> Divide(histMass_LHC22Pass7);

    TH1D *histMass_LHC22Pass7Ratio_fSingleMuLow = (TH1D*) histMass_LHC22Pass7_fSingleMuLow -> Clone("histMass_LHC22Pass7Ratio_fSingleMuLow");   
    histMass_LHC22Pass7Ratio_fSingleMuLow -> Divide(histMass_LHC22Pass7);

    TH1D *histMass_LHC23Pass4Ratio_fDiMuon = (TH1D*) histMass_LHC23Pass4_fDiMuon -> Clone("histMass_LHC23Pass4Ratio_fDiMuon");   
    histMass_LHC23Pass4Ratio_fDiMuon -> Divide(histMass_LHC23Pass4);

    TH1D *histMass_LHC23Pass4Ratio_fSingleMuLow = (TH1D*) histMass_LHC23Pass4_fSingleMuLow -> Clone("histMass_LHC23Pass4Ratio_fSingleMuLow");   
    histMass_LHC23Pass4Ratio_fSingleMuLow -> Divide(histMass_LHC23Pass4);

    SetHistogram(histMass_LHC22Pass7Ratio_fDiMuon, 633, 1, 20, 0.8);
    SetHistogram(histMass_LHC22Pass7Ratio_fSingleMuLow, 862, 1, 20, 0.8);
    SetHistogram(histMass_LHC23Pass4Ratio_fDiMuon, 801, 1, 20, 0.8);
    SetHistogram(histMass_LHC23Pass4Ratio_fSingleMuLow, 417, 1, 20, 0.8);

    TCanvas *canvasMassRatio = new TCanvas("canvasMassRatio", "", 800, 600);
    gPad -> SetLogy(true);
    histMass_LHC23Pass4Ratio_fDiMuon -> SetTitle("");
    histMass_LHC23Pass4Ratio_fDiMuon -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c}^{2})");
    histMass_LHC23Pass4Ratio_fDiMuon -> GetYaxis() -> SetTitle("single trigger / all triggers");
    histMass_LHC23Pass4Ratio_fDiMuon -> GetXaxis() -> SetRangeUser(0, 10);
    histMass_LHC23Pass4Ratio_fDiMuon -> GetYaxis() -> SetRangeUser(0.001, 30);
    histMass_LHC23Pass4Ratio_fDiMuon -> Draw("E");
    histMass_LHC23Pass4Ratio_fSingleMuLow -> Draw("E SAME");
    histMass_LHC22Pass7Ratio_fDiMuon -> Draw("E SAME");
    histMass_LHC22Pass7Ratio_fSingleMuLow -> Draw("E SAME");

    TLine *lineUnityMass = new TLine(0, 1, 10, 1);
    lineUnityMass -> SetLineColor(kGray+1);
    lineUnityMass -> SetLineStyle(kDashed);
    lineUnityMass -> SetLineWidth(2);
    lineUnityMass -> Draw("SAME");

    TLegend *legendMassRatio = new TLegend(0.20, 0.70, 0.45, 0.89);
    SetLegend(legendMassRatio);
    legendMassRatio -> AddEntry(histMass_LHC22Pass7Ratio_fDiMuon, "fDiMuon, LHC22 pass7", "PL");
    legendMassRatio -> AddEntry(histMass_LHC22Pass7Ratio_fSingleMuLow, "fSingleMuLow, LHC22 pass7", "PL");
    legendMassRatio -> AddEntry(histMass_LHC23Pass4Ratio_fDiMuon, "fDiMuon, LHC23 pass4", "PL");
    legendMassRatio -> AddEntry(histMass_LHC23Pass4Ratio_fSingleMuLow, "fSingleMuLow, LHC23 pass4", "PL");
    legendMassRatio -> Draw("SAME");

    ///////////////////
    // Pt dependence //
    ///////////////////
    TLegend *legendPt = new TLegend(0.60, 0.75, 0.85, 0.89);
    SetLegend(legendPt);
    legendPt -> AddEntry(histPt_LHC23Pass4, "All triggers", "PL");
    legendPt -> AddEntry(histPt_LHC23Pass4_fDiMuon, "#bf{fDiMuon} trigger", "L");
    legendPt -> AddEntry(histPt_LHC23Pass4_fSingleMuLow, "#bf{fSingleMuLow} trigger", "L");

    TCanvas *canvasPtLHC22Pass7 = new TCanvas("canvasPtLHC22Pass7", "", 800, 600);
    gPad -> SetLogy(true);

    histPt_LHC22Pass7 -> SetStats(false);
    histPt_LHC22Pass7_fDiMuon -> SetStats(false);
    histPt_LHC22Pass7_fSingleMuLow -> SetStats(false);

    histPt_LHC22Pass7 -> SetTitle("");
    histPt_LHC22Pass7 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histPt_LHC22Pass7 -> GetYaxis() -> SetTitle("Counts");
    histPt_LHC22Pass7 -> GetXaxis() -> SetRangeUser(0, 15);
    histPt_LHC22Pass7 -> Draw("EP");
    histPt_LHC22Pass7_fDiMuon -> Draw("H SAME");
    histPt_LHC22Pass7_fSingleMuLow -> Draw("H SAME");
    legendPt -> Draw("SAME");

    TCanvas *canvasPtLHC23Pass4 = new TCanvas("canvasPtLHC23Pass4", "", 800, 600);
    gPad -> SetLogy(true);

    histPt_LHC23Pass4 -> SetStats(false);
    histPt_LHC23Pass4_fDiMuon -> SetStats(false);
    histPt_LHC23Pass4_fSingleMuLow -> SetStats(false);

    histPt_LHC23Pass4 -> SetTitle("");
    histPt_LHC23Pass4 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histPt_LHC23Pass4 -> GetYaxis() -> SetTitle("Counts");
    histPt_LHC23Pass4 -> GetXaxis() -> SetRangeUser(0, 15);
    histPt_LHC23Pass4 -> Draw("EP");
    histPt_LHC23Pass4_fDiMuon -> Draw("H SAME");
    histPt_LHC23Pass4_fSingleMuLow -> Draw("H SAME");
    legendPt -> Draw("SAME");

    TH1D *histPt_LHC22Pass7Ratio_fDiMuon = (TH1D*) histPt_LHC22Pass7_fDiMuon -> Clone("histPt_LHC22Pass7Ratio_fDiMuon");
    histPt_LHC22Pass7Ratio_fDiMuon -> Divide(histPt_LHC22Pass7);

    TH1D *histPt_LHC22Pass7Ratio_fSingleMuLow = (TH1D*) histPt_LHC22Pass7_fSingleMuLow -> Clone("histPt_LHC22Pass7Ratio_fSingleMuLow");
    histPt_LHC22Pass7Ratio_fSingleMuLow -> Divide(histPt_LHC22Pass7);

    TH1D *histPt_LHC23Pass4Ratio_fDiMuon = (TH1D*) histPt_LHC23Pass4_fDiMuon -> Clone("histPt_LHC23Pass4Ratio_fDiMuon");
    histPt_LHC23Pass4Ratio_fDiMuon -> Divide(histPt_LHC23Pass4);

    TH1D *histPt_LHC23Pass4Ratio_fSingleMuLow = (TH1D*) histPt_LHC23Pass4_fSingleMuLow -> Clone("histPt_LHC23Pass4Ratio_fSingleMuLow");
    histPt_LHC23Pass4Ratio_fSingleMuLow -> Divide(histPt_LHC23Pass4);

    SetHistogram(histPt_LHC22Pass7Ratio_fDiMuon, 633, 1, 20, 0.8);
    SetHistogram(histPt_LHC22Pass7Ratio_fSingleMuLow, 862, 1, 20, 0.8);
    SetHistogram(histPt_LHC23Pass4Ratio_fDiMuon, 801, 1, 20, 0.8);
    SetHistogram(histPt_LHC23Pass4Ratio_fSingleMuLow, 417, 1, 20, 0.8);

    TCanvas *canvasPtRatio = new TCanvas("canvasPtRatio", "", 800, 600);
    gPad ->SetLogy(true);
    histPt_LHC23Pass4Ratio_fDiMuon -> SetTitle("");
    histPt_LHC23Pass4Ratio_fDiMuon -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histPt_LHC23Pass4Ratio_fDiMuon -> GetYaxis() -> SetTitle("single trigger / all triggers");
    histPt_LHC23Pass4Ratio_fDiMuon -> GetXaxis() -> SetRangeUser(0, 15);
    histPt_LHC23Pass4Ratio_fDiMuon -> GetYaxis() -> SetRangeUser(0.001, 30);
    histPt_LHC23Pass4Ratio_fDiMuon -> Draw("EP");
    histPt_LHC23Pass4Ratio_fSingleMuLow -> Draw("EP SAME");
    histPt_LHC22Pass7Ratio_fDiMuon -> Draw("EP SAME");
    histPt_LHC22Pass7Ratio_fSingleMuLow -> Draw("EP SAME");

    TLine *lineUnityPt = new TLine(0, 1, 15, 1);
    lineUnityPt -> SetLineColor(kGray+1);
    lineUnityPt -> SetLineStyle(kDashed);
    lineUnityPt -> SetLineWidth(2);
    lineUnityPt -> Draw("SAME");

    TLegend *legendPtRatio = new TLegend(0.20, 0.70, 0.45, 0.89);
    SetLegend(legendPtRatio);
    legendPtRatio -> AddEntry(histPt_LHC22Pass7Ratio_fDiMuon, "fDiMuon, LHC22 pass7", "PL");
    legendPtRatio -> AddEntry(histPt_LHC22Pass7Ratio_fSingleMuLow, "fSingleMuLow, LHC22 pass7", "PL");
    legendPtRatio -> AddEntry(histPt_LHC23Pass4Ratio_fDiMuon, "fDiMuon, LHC23 pass4", "PL");
    legendPtRatio -> AddEntry(histPt_LHC23Pass4Ratio_fSingleMuLow, "fSingleMuLow, LHC23 pass4", "PL");
    legendPtRatio -> Draw("SAME");

    canvasMassLHC22Pass7 -> SaveAs("mass_LHC22_pass7.pdf");
    canvasMassLHC23Pass4 -> SaveAs("mass_LHC23_pass4.pdf");
    canvasPtLHC22Pass7 -> SaveAs("pt_LHC22_pass7.pdf");
    canvasPtLHC23Pass4 -> SaveAs("pt_LHC23_pass4.pdf");
    canvasMassRatio -> SaveAs("ratio_mass.pdf");
    canvasPtRatio -> SaveAs("ratio_pt.pdf");
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