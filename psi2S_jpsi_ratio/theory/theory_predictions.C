void LoadStyle();
void SetLegend(TLegend *);

void theory_predictions() {
    LoadStyle();

    // taken from https://docs.google.com/spreadsheets/d/1gkvbqAOPlNoyyDsgQ7QtirU1ZQ-LIYxzaXsjY1XG0uI/edit?gid=1699652183#gid=1699652183

    /////////////////////////
    // CGC + NRQCD + FONLL //
    /////////////////////////
    const int nPtBinsCgcNrqcdFonll = 6;
    double ptCentrCgcNrqcdFonll[] = {0.5, 1.5, 2.5, 3.5, 4.5, 6};
    double ptMinCgcNrqcdFonll[] = {0.5, 0.5, 0.5, 0.5, 0.5, 1};
    double ptMaxCgcNrqcdFonll[] = {0.5, 0.5, 0.5, 0.5, 0.5, 1};
    double ptBinWidthCgcNrqcdFonll[] = {1, 1, 1, 1, 1, 2};

    // FONLL
    double csJpsiFonll1Pt[] = {1.45E+02, 3.03E+02, 2.80E+02, 1.99E+02, 1.28E+02, 1.32E+02};

    double csJpsiFonll1PtMinS[] = {6.25E+01, 1.23E+02, 1.04E+02, 6.70E+01, 3.88E+01, 3.44E+01};
    double csJpsiFonll1PtMinM[] = {2.61E+01, 5.06E+01, 4.13E+01, 2.52E+01, 1.40E+01, 1.20E+01};
    double csJpsiFonll1PtMinPdf[] = {2.26E+01, 3.99E+01, 2.76E+01, 1.37E+01, 6.40E+00, 4.71E+00};

    double csJpsiFonll1PtMaxS[] = {6.76E+01, 1.35E+02, 1.21E+02, 8.50E+01, 5.41E+01, 5.33E+01};
    double csJpsiFonll1PtMaxM[] = {3.32E+01, 6.29E+01, 4.96E+01, 2.96E+01, 1.61E+01, 1.34E+01};
    double csJpsiFonll1PtMaxPdf[] = {2.26E+01, 3.99E+01, 2.74E+01, 1.37E+01, 6.43E+00, 4.69E+00};

    // CGC + NRQCD
    double csJpsiCgcNrqcdPt[] = {1.01E+03, 2.47E+03, 2.45E+03, 1.71E+03, 1.02E+03, 9.06E+02};

    double csJpsiCgcNrqcdPtMin[] = {3.29E+02, 8.04E+02, 8.00E+02, 5.58E+02, 3.35E+02, 3.00E+02};
    double csJpsiCgcNrqcdPtMax[] = {3.2920E+02, 8.0424E+02, 7.9951E+02, 5.5786E+02, 3.3482E+02, 2.9998E+02};


    // CGC + NRQCD + FONLL
    double csJpsiCgcNrqcFonllPt[nPtBinsCgcNrqcdFonll], csJpsiCgcNrqcFonllPtMin[nPtBinsCgcNrqcdFonll], csJpsiCgcNrqcFonllPtMax[nPtBinsCgcNrqcdFonll];

    for (int iPt = 0;iPt < nPtBinsCgcNrqcdFonll;iPt++) {
        csJpsiCgcNrqcFonllPt[iPt] = (csJpsiFonll1Pt[iPt] + csJpsiCgcNrqcdPt[iPt]) / ptBinWidthCgcNrqcdFonll[iPt];
        csJpsiCgcNrqcFonllPtMin[iPt] = TMath::Sqrt(TMath::Power(csJpsiCgcNrqcdPtMin[iPt], 2.) + TMath::Power(csJpsiFonll1PtMinS[iPt], 2.) + TMath::Power(csJpsiFonll1PtMinM[iPt], 2.) + TMath::Power(csJpsiFonll1PtMinPdf[iPt], 2.)) / ptBinWidthCgcNrqcdFonll[iPt];
        csJpsiCgcNrqcFonllPtMax[iPt] = TMath::Sqrt(TMath::Power(csJpsiCgcNrqcdPtMax[iPt], 2.) + TMath::Power(csJpsiFonll1PtMaxS[iPt], 2.) + TMath::Power(csJpsiFonll1PtMaxM[iPt], 2.) + TMath::Power(csJpsiFonll1PtMaxPdf[iPt], 2.)) / ptBinWidthCgcNrqcdFonll[iPt];
    }

    ///////////////////
    // NRQCD + FONLL //
    ///////////////////
    const int nPtBinsNrqcdFonll = 5;
    double ptCentrNrqcdFonll[] = {3.5, 4.5, 6, 8.5, 15};
    double ptMinNrqcdFonll[] = {0.5, 0.5, 1, 1.5, 5};
    double ptMaxNrqcdFonll[] = {0.5, 0.5, 1, 1.5, 5};
    double ptBinWidthNrqcdFonll[] = {1, 1, 2, 3, 10};

    // FONLL
    double csJpsiFonll2Pt[] = {1.99E+02, 1.28E+02, 1.32E+02, 6.96E+01, 3.54E+01};

    double csJpsiFonll2PtMinS[] = {6.70E+01, 3.88E+01, 3.44E+01, 1.61E+01, 7.24E+00};
    double csJpsiFonll2PtMinM[] = {2.52E+01, 1.40E+01, 1.20E+01, 4.92E+00, 1.79E+00};
    double csJpsiFonll2PtMinPdf[] = {1.37E+01, 6.40E+00, 4.71E+00, 1.72E+00, 6.26E-01};

    double csJpsiFonll2PtMaxS[] = {8.50E+01, 5.41E+01, 5.33E+01, 2.54E+01, 1.06E+01};
    double csJpsiFonll2PtMaxM[] = {2.96E+01, 1.61E+01, 1.34E+01, 5.38E+00, 1.88E+00};
    double csJpsiFonll2PtMaxPdf[] = {1.37E+01, 6.43E+00, 4.69E+00, 1.72E+00, 6.26E-01};

    // NRQCD
    double csJpsiNrqcdPt[] = {1.98E+03, 8.69E+02, 5.86E+02, 2.05E+02, 6.96E+01};

    double csJpsiNrqcdPtMin[] = {3.04E+02, 2.04E+02, 1.54E+02, 7.00E+01, 2.79E+01};
    double csJpsiNrqcdPtMax[] = {3.04E+02, 2.06E+02, 2.49E+02, 1.24E+02, 5.95E+01};

    // NRQCD + FONLL
    double csJpsiNrqcFonllPt[nPtBinsNrqcdFonll], csJpsiNrqcFonllPtMin[nPtBinsNrqcdFonll], csJpsiNrqcFonllPtMax[nPtBinsNrqcdFonll];

    for (int iPt = 0;iPt < nPtBinsNrqcdFonll;iPt++) {
        csJpsiNrqcFonllPt[iPt] = (csJpsiFonll2Pt[iPt] + csJpsiNrqcdPt[iPt]) / ptBinWidthNrqcdFonll[iPt];
        csJpsiNrqcFonllPtMin[iPt] = TMath::Sqrt(TMath::Power(csJpsiNrqcdPtMin[iPt], 2.) + TMath::Power(csJpsiFonll2PtMinS[iPt], 2.) + TMath::Power(csJpsiFonll2PtMinM[iPt], 2.) + TMath::Power(csJpsiFonll2PtMinPdf[iPt], 2.)) / ptBinWidthCgcNrqcdFonll[iPt];
        csJpsiNrqcFonllPtMax[iPt] = TMath::Sqrt(TMath::Power(csJpsiNrqcdPtMax[iPt], 2.) + TMath::Power(csJpsiFonll2PtMaxS[iPt], 2.) + TMath::Power(csJpsiFonll2PtMaxM[iPt], 2.) + TMath::Power(csJpsiFonll2PtMaxPdf[iPt], 2.)) / ptBinWidthCgcNrqcdFonll[iPt];
    }

    ///////////////////
    // ICEM + FONLL //
    //////////////////
    /*int nPtBinsIcemFonll = 0;
    vector <double> ptCentrIcemFonll;
    vector <double> ptMinIcemFonll;
    vector <double> ptMaxIcemFonll;
    vector <double> ptBinWidthIcemFonll;

    vector <double> csJpsiIcemPt;
    vector <double> csJpsiIcemPtMinRS;
    vector <double> csJpsiIcemPtMinMS;
    vector <double> csJpsiIcemPtMaxRS;
    vector <double> csJpsiIcemPtMaxMS;

    double ptCentrTmp, ptWidthTmp, csTmp, csMinRsTmp, csMaxRsTmp, csMinMsTmp, csMaxMsTmp;
    ifstream fJpsiIcem("jpsi_icem.txt");
    while(!fJpsiIcem.eof()) {
        fJpsiIcem >> ptCentrTmp >> ptWidthTmp >> csTmp >> csMinRsTmp >> csMaxRsTmp >> csMinMsTmp >> csMaxMsTmp;
        ptCentrIcemFonll.push_back(ptCentrTmp - 0.25);
        ptMinIcemFonll.push_back(0.25);
        ptMaxIcemFonll.push_back(0.25);
        ptBinWidthIcemFonll.push_back(0.5);
        csJpsiIcemPt.push_back(csTmp);
        csJpsiIcemPtMinRS.push_back(csTmp - csMinRsTmp);
        csJpsiIcemPtMaxRS.push_back(csMaxRsTmp - csTmp);
        csJpsiIcemPtMinMS.push_back(csTmp - csMinMsTmp);
        csJpsiIcemPtMaxMS.push_back(csMaxMsTmp - csTmp);
        nPtBinsIcemFonll++;
    }*/

    const int nPtBinsIcemFonll = 8;
    double ptCentrIcemFonll[] = {0.5, 1.5, 2.5, 3.5, 4.5, 6, 8.5, 15};
    double ptMinIcemFonll[] = {0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5, 5};
    double ptMaxIcemFonll[] = {0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5, 5};
    double ptBinWidthIcemFonll[] = {1, 1, 1, 1, 1, 2, 3, 10};

    // FONLL
    double csJpsiFonll3Pt[] = {1.45E+02, 3.03E+02, 2.80E+02, 1.99E+02, 1.28E+02, 1.32E+02, 6.96E+01, 3.54E+01};

    double csJpsiFonll3PtMinS[] = {6.25E+01, 1.23E+02, 1.04E+02, 6.70E+01, 3.88E+01, 3.44E+01, 1.61E+01, 7.24E+00};
    double csJpsiFonll3PtMinM[] = {2.61E+01, 5.06E+01, 4.13E+01, 2.52E+01, 1.40E+01, 1.20E+01, 4.92E+00, 1.79E+00};
    double csJpsiFonll3PtMinPdf[] = {2.26E+01, 3.99E+01, 2.76E+01, 1.37E+01, 6.40E+00, 4.71E+00, 1.72E+00, 6.26E-01};

    double csJpsiFonll3PtMaxS[] = {6.76E+01, 1.35E+02, 1.21E+02, 8.50E+01, 5.41E+01, 5.33E+01, 2.54E+01, 1.06E+01};
    double csJpsiFonll3PtMaxM[] = {3.32E+01, 6.29E+01, 4.96E+01, 2.96E+01, 1.61E+01, 1.34E+01, 5.38E+00, 1.88E+00};
    double csJpsiFonll3PtMaxPdf[] = {2.26E+01, 3.99E+01, 2.74E+01, 1.37E+01, 6.43E+00, 4.69E+00, 1.72E+00, 6.26E-01};

    // ICEM
    double csJpsiIcemPt[] = {1052.186373, 2551.68636, 2250.286916, 1309.792997, 646.6338198, 479.9947889, 190.9539356, 95.38180569};

    double csJpsiIcemPtMinRS[] = {376.5235632,898.2031678, 770.3382638, 433.8606155, 207.2702479, 147.7996019, 55.45302669, 25.58886587};
    double csJpsiIcemPtMinMS[] = {409.1338159, 927.2768203, 802.7600056, 469.0383312, 230.079063, 164.0514961, 59.83116239, 27.63605107};

    double csJpsiIcemPtMaxRS[] = {808.7830856, 1893.638049, 1576.317214, 858.9261048, 397.9302254, 273.9155742, 98.0602908, 42.80021763};
    double csJpsiIcemPtMaxMS[] = {135.0836063, 282.5791416, 240.1768786, 139.4607135, 66.55449613, 44.48393196, 14.78375091, 6.408173883};

    // ICEM + FONLL
    double csJpsiIcemFonllPt[nPtBinsIcemFonll], csJpsiIcemFonllPtMin[nPtBinsIcemFonll], csJpsiIcemFonllPtMax[nPtBinsIcemFonll];

    for (int iPt = 0;iPt < nPtBinsIcemFonll;iPt++) {
        csJpsiIcemFonllPt[iPt] = (csJpsiFonll3Pt[iPt] + csJpsiIcemPt[iPt]) / ptBinWidthIcemFonll[iPt];
        csJpsiIcemFonllPtMin[iPt] = TMath::Sqrt(TMath::Power(csJpsiIcemPtMinRS[iPt], 2.) + TMath::Power(csJpsiIcemPtMinMS[iPt], 2.) + TMath::Power(csJpsiFonll3PtMinS[iPt], 2.) + TMath::Power(csJpsiFonll3PtMinM[iPt], 2.) + TMath::Power(csJpsiFonll3PtMinPdf[iPt], 2.)) / ptBinWidthIcemFonll[iPt];
        csJpsiIcemFonllPtMax[iPt] = TMath::Sqrt(TMath::Power(csJpsiIcemPtMaxRS[iPt], 2.) + TMath::Power(csJpsiIcemPtMaxMS[iPt], 2.) + TMath::Power(csJpsiFonll3PtMaxS[iPt], 2.) + TMath::Power(csJpsiFonll3PtMaxM[iPt], 2.) + TMath::Power(csJpsiFonll3PtMaxPdf[iPt], 2.)) / ptBinWidthIcemFonll[iPt];
    }

    // Plot results
    TGraphAsymmErrors *graCsJpsiCgcNrqcFonlldPt = new TGraphAsymmErrors(nPtBinsCgcNrqcdFonll, ptCentrCgcNrqcdFonll, csJpsiCgcNrqcFonllPt, ptMinCgcNrqcdFonll, ptMaxCgcNrqcdFonll, csJpsiCgcNrqcFonllPtMin, csJpsiCgcNrqcFonllPtMax);
    graCsJpsiCgcNrqcFonlldPt -> SetFillColorAlpha(kOrange+7, 0.5);

    TGraphAsymmErrors *graCsJpsiNrqcFonlldPt = new TGraphAsymmErrors(nPtBinsNrqcdFonll, ptCentrNrqcdFonll, csJpsiNrqcFonllPt, ptMinNrqcdFonll, ptMaxNrqcdFonll, csJpsiNrqcFonllPtMin, csJpsiNrqcFonllPtMax);
    graCsJpsiNrqcFonlldPt -> SetFillColorAlpha(kAzure+2, 0.4);

    TGraphAsymmErrors *graCsJpsiIcemFonlldPt = new TGraphAsymmErrors(nPtBinsIcemFonll, ptCentrIcemFonll, csJpsiIcemFonllPt, ptMinIcemFonll, ptMaxIcemFonll, csJpsiIcemFonllPtMin, csJpsiIcemFonllPtMax);
    graCsJpsiIcemFonlldPt -> SetFillStyle(3352);
    graCsJpsiIcemFonlldPt -> SetFillColorAlpha(kMagenta, 0.7);

    TCanvas *canvasCsTheorPt = new TCanvas("canvasCsTheorPt", "", 800, 600);
    canvasCsTheorPt -> SetFillColor(0);
    canvasCsTheorPt -> SetBorderMode(0);
    canvasCsTheorPt -> SetBorderSize(0);
    canvasCsTheorPt -> SetTickx(1);
    canvasCsTheorPt -> SetTicky(1);
    canvasCsTheorPt -> SetLeftMargin(0.15);
    canvasCsTheorPt -> SetBottomMargin(0.1518219);
    canvasCsTheorPt -> SetFrameBorderMode(0);
    canvasCsTheorPt -> SetFrameBorderMode(0);
    gPad -> SetLogy(1);

    TH2D *histGridCsTheorPt = new TH2D("histGridCsTheorPt","",100, 0., 20., 100, 1, 1e5);
    histGridCsTheorPt -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGridCsTheorPt -> GetXaxis() -> SetTitleOffset(1.2);
    histGridCsTheorPt -> GetYaxis() -> SetTitle("d^{2}#sigma / d#it{p}_{T} d#it{y} (nb GeV^{-1})");
    histGridCsTheorPt -> GetYaxis() -> SetTitleOffset(1.3);
    histGridCsTheorPt -> Draw();

    graCsJpsiCgcNrqcFonlldPt -> Draw("E2 SAME");
    graCsJpsiNrqcFonlldPt -> Draw("E2 SAME");
    graCsJpsiIcemFonlldPt -> Draw("E2 SAME");

    TLegend *legendCsTheorPt = new TLegend(0.40, 0.65, 0.60, 0.8);
    SetLegend(legendCsTheorPt);
    legendCsTheorPt -> AddEntry(graCsJpsiCgcNrqcFonlldPt, "CGC + NRQCD + FONLL", "F");
    legendCsTheorPt -> AddEntry(graCsJpsiNrqcFonlldPt, "NRQCD + FONLL", "F");
    legendCsTheorPt -> AddEntry(graCsJpsiIcemFonlldPt, "ICEM + FONLL", "F");
    legendCsTheorPt -> Draw("SAME");
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