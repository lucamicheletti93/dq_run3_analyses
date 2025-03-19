void LoadStyle();
void SetLegend(TLegend *);

inline void FillPointsForCl(int nPoints, double min, double max, TF1 *func, double xFit[], double yFit[], double yErr[]) {
    for (int i = 0; i < nPoints; i++) {
        xFit[i] = min + i * (max - min) / (nPoints - 1);
        yFit[i] = func -> Eval(xFit[i]);
    }
}

inline int GetClosestPointIndex(TGraphErrors *graph, double x_target) {
    int nPoints = graph -> GetN();
    double x, y;
    int bestIndex = -1;
    double minDiff = 1e-2;

    for (int i = 0; i < nPoints; i++) {
        graph -> GetPoint(i, x, y);
        double diff = fabs(x - x_target);
        if (diff < minDiff) {
            minDiff = diff;
            bestIndex = i;
        }
    }
    return bestIndex;
}

void extrapolation(double mySqrts = 5.36, double confLevel = 0.68) {
    LoadStyle();

    double extrSqrts[1] = {5.36};

    TH1D *histRatiosErrAll = new TH1D("histRatiosErrAll", "", 3, 0, 3);
    histRatiosErrAll -> SetLineColor(kBlack);
    histRatiosErrAll -> SetLineWidth(2);
    histRatiosErrAll -> SetMarkerStyle(20);
    histRatiosErrAll -> SetMarkerColor(kBlack);
    histRatiosErrAll -> SetFillStyle(0);

    TH1D *histRatiosErrStat = new TH1D("histRatiosErrStat", "", 3, 0, 3);
    histRatiosErrStat -> SetLineColor(kBlack);
    histRatiosErrStat -> SetLineWidth(2);
    histRatiosErrStat -> SetMarkerStyle(24);
    histRatiosErrStat -> SetMarkerColor(kBlack);
    histRatiosErrStat -> SetFillStyle(0);

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

    const int nPoints = 1000;
    double xFit[nPoints], yFit[nPoints], yErr[nPoints];
    // Pol1
    TF1 *funcErrAllPol1 = new TF1("funcErrAllPol1", "pol1", 2, 15);
    funcErrAllPol1 -> SetLineStyle(kSolid);
    funcErrAllPol1 -> SetLineColor(kRed+1);
    graErrAllRatioVsSqrts -> Fit(funcErrAllPol1, "R0");

    FillPointsForCl(nPoints, 4, 14, funcErrAllPol1, xFit, yFit, yErr);
    TGraphErrors *graFuncErrAllPol1Cl = new TGraphErrors(nPoints, xFit, yFit, 0, yErr);
    graFuncErrAllPol1Cl -> SetName("graFuncErrAllPol1Cl");
    TVirtualFitter *fitterFuncErrAllPol1Cl = TVirtualFitter::GetFitter();
    fitterFuncErrAllPol1Cl -> GetConfidenceIntervals(graFuncErrAllPol1Cl, 0.68);
    graFuncErrAllPol1Cl -> SetFillColorAlpha(kRed+1, 0.3);

    TF1 *funcErrStatPol1 = new TF1("funcErrStatPol1", "pol1", 2, 15);
    funcErrStatPol1 -> SetLineStyle(kDashed);
    funcErrStatPol1 -> SetLineColor(kRed+1);
    graStatRatioVsSqrts -> Fit(funcErrStatPol1, "R0");

    FillPointsForCl(nPoints, 4, 14, funcErrStatPol1, xFit, yFit, yErr);
    TGraphErrors *graFuncErrStatPol1Cl = new TGraphErrors(nPoints, xFit, yFit, 0, yErr);
    graFuncErrStatPol1Cl -> SetName("graFuncErrStatPol1Cl");
    TVirtualFitter *fitterFuncErrStatPol1Cl = TVirtualFitter::GetFitter();
    fitterFuncErrStatPol1Cl -> GetConfidenceIntervals(graFuncErrStatPol1Cl, 0.68);
    graFuncErrStatPol1Cl -> SetFillStyle(3352);
    graFuncErrStatPol1Cl -> SetFillColorAlpha(kRed+1, 0.3);

    int indexPol1 = GetClosestPointIndex(graFuncErrAllPol1Cl, mySqrts);
    double extrRatioErrAllPol1[1], errExtrRatioErrAllPol1[1];
    extrRatioErrAllPol1[0] = funcErrAllPol1 -> Eval(mySqrts);
    errExtrRatioErrAllPol1[0] = graFuncErrAllPol1Cl -> GetErrorY(indexPol1);

    double extrRatioErrStatPol1[1], errExtrRatioErrStatPol1[1];
    extrRatioErrStatPol1[0] = funcErrStatPol1 -> Eval(mySqrts);
    errExtrRatioErrStatPol1[0] = graFuncErrStatPol1Cl -> GetErrorY(indexPol1);

    histRatiosErrAll -> SetBinContent(1, extrRatioErrAllPol1[0]);
    histRatiosErrAll -> SetBinError(1, errExtrRatioErrAllPol1[0]);

    histRatiosErrStat -> SetBinContent(1, extrRatioErrStatPol1[0]);
    histRatiosErrStat -> SetBinError(1, errExtrRatioErrStatPol1[0]);

    TGraphErrors *graExtrRatioErrAllPol1 = new TGraphErrors(1, extrSqrts, extrRatioErrAllPol1, 0, errExtrRatioErrAllPol1);
    graExtrRatioErrAllPol1 -> SetLineColor(kRed+1);
    graExtrRatioErrAllPol1 -> SetLineWidth(2);
    graExtrRatioErrAllPol1 -> SetMarkerStyle(20);
    graExtrRatioErrAllPol1 -> SetMarkerColor(kRed+1);
    graExtrRatioErrAllPol1 -> SetFillStyle(0);

    TGraphErrors *graExtrRatioErrStatPol1 = new TGraphErrors(1, extrSqrts, extrRatioErrStatPol1, 0, errExtrRatioErrStatPol1);
    graExtrRatioErrStatPol1 -> SetLineColor(kRed+1);
    graExtrRatioErrStatPol1 -> SetLineWidth(2);
    graExtrRatioErrStatPol1 -> SetMarkerStyle(24);
    graExtrRatioErrStatPol1 -> SetMarkerColor(kRed+1);
    graExtrRatioErrStatPol1 -> SetFillStyle(0);
    
    // Exponential
    TF1 *funcErrAllExp = new TF1("funcErrAllExp", "[0]*(1-exp(-[1]*x))", 2, 15);
    funcErrAllExp -> SetLineColor(kBlue+1);
    funcErrAllExp -> SetParameters(0.01, 1);
    graErrAllRatioVsSqrts -> Fit(funcErrAllExp, "R0");

    FillPointsForCl(nPoints, 4, 14, funcErrAllExp, xFit, yFit, yErr);
    TGraphErrors *graFuncErrAllExpCl = new TGraphErrors(nPoints, xFit, yFit, 0, yErr);
    graFuncErrAllExpCl -> SetName("graFuncErrAllExpCl");
    TVirtualFitter *fitterFuncErrAllExpCl = TVirtualFitter::GetFitter();
    fitterFuncErrAllExpCl -> GetConfidenceIntervals(graFuncErrAllExpCl, 0.68);
    graFuncErrAllExpCl -> SetFillColorAlpha(kBlue+1, 0.3);

    TF1 *funcErrStatExp = new TF1("funcErrStatExp", "pol1", 2, 15);
    funcErrStatExp -> SetLineStyle(kDashed);
    funcErrStatExp -> SetLineColor(kBlue+1);
    graStatRatioVsSqrts -> Fit(funcErrStatExp, "R0");

    FillPointsForCl(nPoints, 4, 14, funcErrStatExp, xFit, yFit, yErr);
    TGraphErrors *graFuncErrStatExpCl = new TGraphErrors(nPoints, xFit, yFit, 0, yErr);
    graFuncErrStatExpCl -> SetName("graFuncErrStatExpCl");
    TVirtualFitter *fitterFuncErrStatExpCl = TVirtualFitter::GetFitter();
    fitterFuncErrStatExpCl -> GetConfidenceIntervals(graFuncErrStatExpCl, 0.68);
    graFuncErrStatExpCl -> SetFillStyle(3352);
    graFuncErrStatExpCl -> SetFillColorAlpha(kBlue+1, 0.3);

    int indexExp = GetClosestPointIndex(graFuncErrAllExpCl, mySqrts);
    double extrRatioErrAllExp[1], errExtrRatioErrAllExp[1];
    extrRatioErrAllExp[0] = funcErrAllExp -> Eval(mySqrts);
    errExtrRatioErrAllExp[0] = graFuncErrAllExpCl -> GetErrorY(indexExp);

    double extrRatioErrStatExp[1], errExtrRatioErrStatExp[1];
    extrRatioErrStatExp[0] = funcErrStatExp -> Eval(mySqrts);
    errExtrRatioErrStatExp[0] = graFuncErrStatExpCl -> GetErrorY(indexExp);

    histRatiosErrAll -> SetBinContent(2, extrRatioErrAllExp[0]);
    histRatiosErrAll -> SetBinError(2, errExtrRatioErrAllExp[0]);

    histRatiosErrStat -> SetBinContent(2, extrRatioErrStatExp[0]);
    histRatiosErrStat -> SetBinError(2, errExtrRatioErrStatExp[0]);

    TGraphErrors *graExtrRatioErrAllExp = new TGraphErrors(1, extrSqrts, extrRatioErrAllExp, 0, errExtrRatioErrAllExp);
    graExtrRatioErrAllExp -> SetLineColor(kBlue+1);
    graExtrRatioErrAllExp -> SetLineWidth(2);
    graExtrRatioErrAllExp -> SetMarkerStyle(20);
    graExtrRatioErrAllExp -> SetMarkerColor(kBlue+1);
    graExtrRatioErrAllExp -> SetFillStyle(0);

    TGraphErrors *graExtrRatioErrStatExp = new TGraphErrors(1, extrSqrts, extrRatioErrStatExp, 0, errExtrRatioErrStatExp);
    graExtrRatioErrStatExp -> SetLineColor(kBlue+1);
    graExtrRatioErrStatExp -> SetLineWidth(2);
    graExtrRatioErrStatExp -> SetMarkerStyle(24);
    graExtrRatioErrStatExp -> SetMarkerColor(kBlue+1);
    graExtrRatioErrStatExp -> SetFillStyle(0);

    // Power Low
    TF1 *funcErrAllPowerLow = new TF1("funcErrAllPowerLow", "[0]*x**([1])", 2, 15);
    funcErrAllPowerLow -> SetLineColor(kGreen+1);
    graErrAllRatioVsSqrts -> Fit(funcErrAllPowerLow, "R0");

    FillPointsForCl(nPoints, 4, 14, funcErrAllPowerLow, xFit, yFit, yErr);
    TGraphErrors *graFuncErrAllPowerLowCl = new TGraphErrors(nPoints, xFit, yFit, 0, yErr);
    graFuncErrAllPowerLowCl -> SetName("graFuncErrAllPowerLowCl");
    TVirtualFitter *fitterFuncErrAllPowerLowCl = TVirtualFitter::GetFitter();
    fitterFuncErrAllPowerLowCl -> GetConfidenceIntervals(graFuncErrAllPowerLowCl, 0.68);
    graFuncErrAllPowerLowCl -> SetFillColorAlpha(kGreen+1, 0.3);

    TF1 *funcErrStatPowerLow = new TF1("funcErrStatPowerLow", "pol1", 2, 15);
    funcErrStatPowerLow -> SetLineStyle(kDashed);
    funcErrStatPowerLow -> SetLineColor(kGreen+1);
    graStatRatioVsSqrts -> Fit(funcErrStatPowerLow, "R0");

    FillPointsForCl(nPoints, 4, 14, funcErrStatPowerLow, xFit, yFit, yErr);
    TGraphErrors *graFuncErrStatPowerLowCl = new TGraphErrors(nPoints, xFit, yFit, 0, yErr);
    graFuncErrStatPowerLowCl -> SetName("graFuncErrStatPowerLowCl");
    TVirtualFitter *fitterFuncErrStatPowerLowCl = TVirtualFitter::GetFitter();
    fitterFuncErrStatPowerLowCl -> GetConfidenceIntervals(graFuncErrStatPowerLowCl, 0.68);
    graFuncErrStatPowerLowCl -> SetFillStyle(3352);
    graFuncErrStatPowerLowCl -> SetFillColorAlpha(kGreen+1, 0.3);

    int indexPowerLow = GetClosestPointIndex(graFuncErrAllPowerLowCl, mySqrts);
    double extrRatioErrAllPowerLow[1], errExtrRatioErrAllPowerLow[1];
    extrRatioErrAllPowerLow[0] = funcErrAllPowerLow -> Eval(mySqrts);
    errExtrRatioErrAllPowerLow[0] = graFuncErrAllPowerLowCl -> GetErrorY(indexPowerLow);

    double extrRatioErrStatPowerLow[1], errExtrRatioErrStatPowerLow[1];
    extrRatioErrStatPowerLow[0] = funcErrStatPowerLow -> Eval(mySqrts);
    errExtrRatioErrStatPowerLow[0] = graFuncErrStatPowerLowCl -> GetErrorY(indexPowerLow);

    histRatiosErrAll -> SetBinContent(3, extrRatioErrAllPowerLow[0]);
    histRatiosErrAll -> SetBinError(3, errExtrRatioErrAllPowerLow[0]);

    histRatiosErrStat -> SetBinContent(3, extrRatioErrStatPowerLow[0]);
    histRatiosErrStat -> SetBinError(3, errExtrRatioErrStatPowerLow[0]);

    TGraphErrors *graExtrRatioErrAllPowerLow = new TGraphErrors(1, extrSqrts, extrRatioErrAllPowerLow, 0, errExtrRatioErrAllPowerLow);
    graExtrRatioErrAllPowerLow -> SetLineColor(kGreen+1);
    graExtrRatioErrAllPowerLow -> SetLineWidth(2);
    graExtrRatioErrAllPowerLow -> SetMarkerStyle(20);
    graExtrRatioErrAllPowerLow -> SetMarkerColor(kGreen+1);
    graExtrRatioErrAllPowerLow -> SetFillStyle(0);

    TGraphErrors *graExtrRatioErrStatPowerLow = new TGraphErrors(1, extrSqrts, extrRatioErrStatPowerLow, 0, errExtrRatioErrStatPowerLow);
    graExtrRatioErrStatPowerLow -> SetLineColor(kGreen+1);
    graExtrRatioErrStatPowerLow -> SetLineWidth(2);
    graExtrRatioErrStatPowerLow -> SetMarkerStyle(24);
    graExtrRatioErrStatPowerLow -> SetMarkerColor(kGreen+1);
    graExtrRatioErrStatPowerLow -> SetFillStyle(0);

    // Extract the average of the extrapolated points
    TF1 *funcPol0ErrAll = new TF1("funcPol0ErrAll", "pol0", 0, 3);
    histRatiosErrAll -> Fit(funcPol0ErrAll, "R0");

    TF1 *funcPol0ErrStat = new TF1("funcPol0ErrStat", "pol0", 0, 3);
    histRatiosErrStat -> Fit(funcPol0ErrStat, "R0");


    TCanvas *canvasRatioVsSqrts = new TCanvas("canvasRatioVsSqrts", "", 1200, 1200);
    canvasRatioVsSqrts -> Divide(2, 2);

    TH2D *histGridVsSqrts = new TH2D("histGridVsSqrts", "; #sqrt{#it{s}} (GeV) ; #sigma_{#psi(2S)} / #sigma_{J/#psi}", 100, 4, 14, 100, 0.1, 0.2);
    TH2D *histGridVsTrials = new TH2D("histGridVsTrials", "; Trial ; #sigma_{#psi(2S)} / #sigma_{J/#psi}", 3, 0, 3, 100, 0.1, 0.2);

    canvasRatioVsSqrts -> cd(1);
    histGridVsSqrts -> Draw();
    graStatRatioVsSqrts -> Draw("EP SAME");
    graSystAllRatioVsSqrts -> Draw("E2 SAME");
    funcErrAllPol1 -> Draw("SAME");
    //funcErrStatPol1 -> Draw("SAME");
    graFuncErrAllPol1Cl -> Draw("3 SAME");
    //graFuncErrStatPol1Cl -> Draw("3 SAME");
    graExtrRatioErrAllPol1 -> Draw("EP SAME");
    //graExtrRatioErrStatPol1 -> Draw("EP SAME");

    canvasRatioVsSqrts -> cd(2);
    histGridVsSqrts -> Draw();
    graStatRatioVsSqrts -> Draw("EP SAME");
    graSystAllRatioVsSqrts -> Draw("E2 SAME");
    funcErrAllExp -> Draw("SAME");
    //funcErrStatExp -> Draw("SAME");
    graFuncErrAllExpCl -> Draw("3 SAME");
    //graFuncErrStatExpCl -> Draw("3 SAME");
    graExtrRatioErrAllExp -> Draw("EP SAME");
    //graExtrRatioErrStatExp -> Draw("EP SAME");

    canvasRatioVsSqrts -> cd(3);
    histGridVsSqrts -> Draw();
    graStatRatioVsSqrts -> Draw("EP SAME");
    graSystAllRatioVsSqrts -> Draw("E2 SAME");
    funcErrAllPowerLow -> Draw("SAME");
    //funcErrStatPowerLow -> Draw("SAME");
    graFuncErrAllPowerLowCl -> Draw("3 SAME");
    //graFuncErrStatPowerLowCl -> Draw("3 SAME");
    graExtrRatioErrAllPowerLow -> Draw("EP SAME");
    //graExtrRatioErrStatPowerLow -> Draw("EP SAME");

    canvasRatioVsSqrts -> cd(4);
    histGridVsTrials -> Draw();
    funcPol0ErrAll -> Draw("SAME");
    funcPol0ErrStat -> Draw("SAME");
    histRatiosErrAll -> Draw("EP SAME");
    //histRatiosErrStat -> Draw("EP SAME");

    double meanErrStat = 0, meanErrAll = 0;
    for (int iTrial = 0;iTrial < 3;iTrial++) {
        meanErrStat += histRatiosErrStat -> GetBinError(iTrial+1);
        meanErrAll += histRatiosErrAll -> GetBinError(iTrial+1);
    }
    meanErrStat = meanErrStat / 3.;
    meanErrAll = meanErrAll / 3.;

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.050);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);
    //latexTitle -> DrawLatex(0.2, 0.87, Form("#sigma_{#psi(2S)}/#sigma_{J/#psi} = %f #pm %f #pm %f", funcPol0ErrAll -> GetParameter(0), funcPol0ErrAll -> GetParError(0), meanErrAll));
    //latexTitle -> DrawLatex(0.2, 0.82, Form("#sigma_{#psi(2S)}/#sigma_{J/#psi} = %f #pm %f #pm %f", funcPol0ErrStat -> GetParameter(0), funcPol0ErrStat -> GetParError(0), meanErrStat));
    latexTitle -> DrawLatex(0.25, 0.87, Form("#sigma_{#psi(2S)}/#sigma_{J/#psi} = %f #pm %f", funcPol0ErrAll -> GetParameter(0), funcPol0ErrAll -> GetParError(0)));

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