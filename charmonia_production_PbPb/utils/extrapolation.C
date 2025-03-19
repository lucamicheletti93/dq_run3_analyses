void LoadStyle();
void SetLegend(TLegend *);

inline void FillPointsForCl(int nPoints, double min, double max, TF1 *func, double xFit[], double yFit[], double yErr[]) {
    for (int i = 0; i < nPoints; i++) {
        xFit[i] = min + i * (max - min) / (nPoints - 1);
        yFit[i] = func -> Eval(xFit[i]);
    }
}

void extrapolation(double mySqrts = 5.36, double confLevel = 0.68) {
    LoadStyle();

    double extrSqrts[1] = {5.36};
    vector<double> extrRatios, errExtrRatios;
    const int nSqrts = 4;

    double sqrts[] = {5.02, 7.00, 8.00, 13.00};
    double errSqrts[] = {0.15, 0.15, 0.15, 0.15};
    double ratioVsSqrts[] = {0.147, 0.170, 0.140, 0.146};
    double statRatioVsSqrts[] = {0.009, 0.011, 0.010, 0.005};
    double systRatioVsSqrts[] = {0.01, 0.013, 0.020, 0.009};
    double systBrRatioVsSqrts[] = {0.011, 0.013, 0.011, 0.017};
    double systAllRatioVsSqrts[nSqrts], errAllRatioVsSqrts[nSqrts];

    for (int iSqrts = 0;iSqrts < nSqrts;iSqrts++) {
        systAllRatioVsSqrts[iSqrts] = ratioVsSqrts[iSqrts] * TMath::Sqrt((systRatioVsSqrts[iSqrts]/ratioVsSqrts[iSqrts])*(systRatioVsSqrts[iSqrts]/ratioVsSqrts[iSqrts]) + (systBrRatioVsSqrts[iSqrts]/ratioVsSqrts[iSqrts])*(systBrRatioVsSqrts[iSqrts]/ratioVsSqrts[iSqrts]));
        errAllRatioVsSqrts[iSqrts] = ratioVsSqrts[iSqrts] * TMath::Sqrt((statRatioVsSqrts[iSqrts]/ratioVsSqrts[iSqrts])*(statRatioVsSqrts[iSqrts]/ratioVsSqrts[iSqrts]) + (systRatioVsSqrts[iSqrts]/ratioVsSqrts[iSqrts])*(systRatioVsSqrts[iSqrts]/ratioVsSqrts[iSqrts]));
    }

    double brJpsi = 5.961;
    double errBrJpsi = 0.033;
    double brPsi2s = 8.0e-3;
    double errBrPsi2s = 0.6e-3;
    double errBrRatio = TMath::Sqrt((errBrJpsi/brJpsi)*(errBrJpsi/brJpsi) + (errBrPsi2s/brPsi2s)*(errBrPsi2s/brPsi2s));


    TGraphErrors *graStatRatioVsSqrts = new TGraphErrors(nSqrts, sqrts, ratioVsSqrts, 0, statRatioVsSqrts);
    graStatRatioVsSqrts -> SetLineColor(kBlack);
    graStatRatioVsSqrts -> SetLineWidth(2);
    graStatRatioVsSqrts -> SetMarkerStyle(20);
    graStatRatioVsSqrts -> SetMarkerColor(kBlack);

    TGraphErrors *graSystRatioVsSqrts = new TGraphErrors(nSqrts, sqrts, ratioVsSqrts, errSqrts, systRatioVsSqrts);
    graSystRatioVsSqrts -> SetLineColor(kBlack);
    graSystRatioVsSqrts -> SetLineWidth(2);
    graSystRatioVsSqrts -> SetMarkerStyle(20);
    graSystRatioVsSqrts -> SetMarkerColor(kBlack);
    graSystRatioVsSqrts -> SetFillStyle(0);

    TGraphErrors *graSystAllRatioVsSqrts = new TGraphErrors(nSqrts, sqrts, ratioVsSqrts, errSqrts, systAllRatioVsSqrts);
    graSystAllRatioVsSqrts -> SetLineColor(kBlack);
    graSystAllRatioVsSqrts -> SetLineWidth(2);
    graSystAllRatioVsSqrts -> SetMarkerStyle(20);
    graSystAllRatioVsSqrts -> SetMarkerColor(kBlack);
    graSystAllRatioVsSqrts -> SetFillStyle(0);

    TGraphErrors *graErrAllRatioVsSqrts = new TGraphErrors(nSqrts, sqrts, ratioVsSqrts, errSqrts, errAllRatioVsSqrts);
    graErrAllRatioVsSqrts -> SetLineColor(kBlack);
    graErrAllRatioVsSqrts -> SetLineWidth(2);
    graErrAllRatioVsSqrts -> SetMarkerStyle(20);
    graErrAllRatioVsSqrts -> SetMarkerColor(kBlack);
    graErrAllRatioVsSqrts -> SetFillStyle(0);

    // Pol1
    TF1 *funcPol1 = new TF1("funcPol1", "pol1", 2, 15);
    funcPol1 -> SetLineColor(kRed+1);
    graErrAllRatioVsSqrts -> Fit(funcPol1, "R0");

    const int nPoints = 1000;
    double xFit[nPoints], yFit[nPoints], yErr[nPoints];
    FillPointsForCl(nPoints, 4, 14, funcPol1, xFit, yFit, yErr);
    TGraphErrors *graFuncPol1Cl = new TGraphErrors(nPoints, xFit, yFit, 0, yErr);
    TVirtualFitter *fitterFuncPol1Cl = TVirtualFitter::GetFitter();
    fitterFuncPol1Cl -> GetConfidenceIntervals(graFuncPol1Cl, 0.68);
    graFuncPol1Cl -> SetFillColorAlpha(kRed+1, 0.3);

    int indexPol1 = graFuncPol1Cl -> GetN() - 1;
    double extrRatioPol1[1], errExtrRatioPol1[1];
    extrRatioPol1[0] = funcPol1 -> Eval(mySqrts);
    errExtrRatioPol1[0] = graFuncPol1Cl -> GetErrorY(indexPol1);
    extrRatios.push_back(extrRatioPol1[0]);
    errExtrRatios.push_back(errExtrRatioPol1[0]);

    TGraphErrors *graExtrRatioPol1 = new TGraphErrors(1, extrSqrts, extrRatioPol1, 0, errExtrRatioPol1);
    graExtrRatioPol1 -> SetLineColor(kRed+1);
    graExtrRatioPol1 -> SetLineWidth(2);
    graExtrRatioPol1 -> SetMarkerStyle(20);
    graExtrRatioPol1 -> SetMarkerColor(kRed+1);
    graExtrRatioPol1 -> SetFillStyle(0);

    // Exponential
    TF1 *funcExp = new TF1("funcExp", "[0]*(1-exp(-[1]*x))", 2, 15);
    funcExp -> SetLineColor(kBlue+1);
    funcExp -> SetParameters(0.01, 1);
    graErrAllRatioVsSqrts -> Fit(funcExp, "R0");

    FillPointsForCl(nPoints, 4, 14, funcExp, xFit, yFit, yErr);
    TGraphErrors *graFuncExpCl = new TGraphErrors(nPoints, xFit, yFit, 0, yErr);
    TVirtualFitter *fitterFuncExpCl = TVirtualFitter::GetFitter();
    fitterFuncExpCl -> GetConfidenceIntervals(graFuncExpCl, 0.68);
    graFuncExpCl -> SetFillColorAlpha(kBlue+1, 0.3);

    int indexExp = graFuncExpCl -> GetN() - 1;
    double extrRatioExp[1], errExtrRatioExp[1];
    extrRatioExp[0] = funcExp -> Eval(mySqrts);
    errExtrRatioExp[0] = graFuncExpCl -> GetErrorY(indexExp);
    extrRatios.push_back(extrRatioExp[0]);
    errExtrRatios.push_back(errExtrRatioExp[0]);

    TGraphErrors *graExtrRatioExp = new TGraphErrors(1, extrSqrts, extrRatioExp, 0, errExtrRatioExp);
    graExtrRatioExp -> SetLineColor(kBlue+1);
    graExtrRatioExp -> SetLineWidth(2);
    graExtrRatioExp -> SetMarkerStyle(20);
    graExtrRatioExp -> SetMarkerColor(kBlue+1);
    graExtrRatioExp -> SetFillStyle(0);

    // Power Low
    TF1 *funcPowerLow = new TF1("funcPowerLow", "[0]*x**([1])", 2, 15);
    funcPowerLow -> SetLineColor(kGreen+1);
    graErrAllRatioVsSqrts -> Fit(funcPowerLow, "R0");

    FillPointsForCl(nPoints, 4, 14, funcPowerLow, xFit, yFit, yErr);
    TGraphErrors *graFuncPowerLowCl = new TGraphErrors(nPoints, xFit, yFit, 0, yErr);
    TVirtualFitter *fitterFuncPowerLowCl = TVirtualFitter::GetFitter();
    fitterFuncPowerLowCl -> GetConfidenceIntervals(graFuncPowerLowCl, 0.68);
    graFuncPowerLowCl -> SetFillColorAlpha(kGreen+1, 0.3);

    int indexPowerLow = graFuncPowerLowCl -> GetN() - 1;
    double extrRatioPowerLow[1], errExtrRatioPowerLow[1];
    extrRatioPowerLow[0] = funcPowerLow -> Eval(mySqrts);
    errExtrRatioPowerLow[0] = graFuncPowerLowCl -> GetErrorY(indexPowerLow);
    extrRatios.push_back(extrRatioPowerLow[0]);
    errExtrRatios.push_back(errExtrRatioPowerLow[0]);

    TGraphErrors *graExtrRatioPowerLow = new TGraphErrors(1, extrSqrts, extrRatioPowerLow, 0, errExtrRatioPowerLow);
    graExtrRatioPowerLow -> SetLineColor(kGreen+1);
    graExtrRatioPowerLow -> SetLineWidth(2);
    graExtrRatioPowerLow -> SetMarkerStyle(20);
    graExtrRatioPowerLow -> SetMarkerColor(kGreen+1);
    graExtrRatioPowerLow -> SetFillStyle(0);

    // Extract the average of the extrapolated points
    double trials[3] = {0.5, 1.5, 2.5};
    TGraphErrors *graEstraRatios = new TGraphErrors(3, trials, &(extrRatios[0]), 0, &(errExtrRatios[0]));
    graEstraRatios -> SetLineColor(kBlack);
    graEstraRatios -> SetLineWidth(2);
    graEstraRatios -> SetMarkerStyle(20);
    graEstraRatios -> SetMarkerColor(kBlack);
    graEstraRatios -> SetFillStyle(0);

    TF1 *funcPol0 = new TF1("funcPol0", "pol0", 0, 3);
    graEstraRatios -> Fit(funcPol0, "R0");

    std::cout << "Extrapolated Ratio: " << funcPol0 -> GetParameter(0) << " +/- " << funcPol0 -> GetParError(0) << std::endl;


    TCanvas *canvasRatioVsSqrts = new TCanvas("canvasRatioVsSqrts", "", 1200, 1200);
    canvasRatioVsSqrts -> Divide(2, 2);

    TH2D *histGridVsSqrts = new TH2D("histGridVsSqrts", "; #sqrt{#it{s}} (GeV) ; #sigma_{#psi(2S)} / #sigma_{J/#psi}", 100, 4, 14, 100, 0, 0.3);
    TH2D *histGridVsTrials = new TH2D("histGridVsTrials", "; Trial ; #sigma_{#psi(2S)} / #sigma_{J/#psi}", 3, 0, 3, 100, 0, 0.3);

    canvasRatioVsSqrts -> cd(1);
    histGridVsSqrts -> Draw();
    graStatRatioVsSqrts -> Draw("EP SAME");
    graSystAllRatioVsSqrts -> Draw("E2 SAME");
    funcPol1 -> Draw("SAME");
    graFuncPol1Cl -> Draw("3 SAME");
    graExtrRatioPol1 -> Draw("EP SAME");

    canvasRatioVsSqrts -> cd(2);
    histGridVsSqrts -> Draw();
    graStatRatioVsSqrts -> Draw("EP SAME");
    graSystAllRatioVsSqrts -> Draw("E2 SAME");
    funcExp -> Draw("SAME");
    graFuncExpCl -> Draw("3 SAME");
    graExtrRatioExp -> Draw("EP SAME");

    canvasRatioVsSqrts -> cd(3);
    histGridVsSqrts -> Draw();
    graStatRatioVsSqrts -> Draw("EP SAME");
    graSystAllRatioVsSqrts -> Draw("E2 SAME");
    funcPowerLow -> Draw("SAME");
    graFuncPowerLowCl -> Draw("3 SAME");
    graExtrRatioPowerLow -> Draw("EP SAME");

    canvasRatioVsSqrts -> cd(4);
    histGridVsTrials -> Draw();
    graEstraRatios -> Draw("EP SAME");
    funcPol0 -> Draw("SAME");

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.050);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);
    latexTitle -> DrawLatex(0.25, 0.87, Form("#sigma_{#psi(2S)} / #sigma_{J/#psi} = %f #pm %f", funcPol0 -> GetParameter(0), funcPol0 -> GetParError(0)));

    canvasRatioVsSqrts -> SaveAs("extrapolation.pdf");
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