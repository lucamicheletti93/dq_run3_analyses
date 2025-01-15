double FuncBkg1(double *, double *);
double FuncBkg2(double *, double *);
double FuncBkg3(double *, double *);
double FuncBkg4(double *, double *);
double FuncPol2(double *, double *);
double FuncGaus(double *, double *);
double FuncVWG(double *, double *);
double FuncWeight1(double *, double *);
double FuncWeight2(double *, double *);

void extract_bkg_shape() {
    TFile *fIn = new TFile("toy_mc_output.root", "READ");
    fIn -> ls();

    TH1D* histMassProj = (TH1D*) fIn -> Get("histMass_Pt_1.5_2.0");
    histMassProj -> SetLineColor(kRed);
    histMassProj -> SetMarkerColor(kRed);
    histMassProj -> Sumw2(true);

    TH1D* histMassProjPtCutRun3 = (TH1D*) fIn -> Get("histMassPtCut_Pt_1.5_2.0");
    histMassProjPtCutRun3 -> SetLineColor(kBlack);
    histMassProjPtCutRun3 -> SetMarkerColor(kBlack);
    histMassProjPtCutRun3 -> Sumw2(true);

    //histMassProj -> Scale(1. / histMassProj -> Integral());
    //histMassProjPtCutRun3 -> Scale(1. / histMassProjPtCutRun3 -> Integral());

    TH1D *histRatioMass = (TH1D*) histMassProj -> Clone("histRatioMass");
    histRatioMass -> Divide(histMassProjPtCutRun3);
    histRatioMass -> SetLineColor(kBlack);
    histRatioMass -> SetMarkerColor(kBlack);

    TF1 *funcWeight1 = new TF1("funcWeight1", FuncWeight1, 1.9, 5., 4);
    funcWeight1 -> SetParameter(0, 1);
    funcWeight1 -> SetParameter(1, 1);
    funcWeight1 -> SetParameter(2, 1);
    funcWeight1 -> SetParameter(3, 1);

    histRatioMass -> Fit(funcWeight1, "R0");
    funcWeight1 -> SetLineColor(kMagenta);

    TF1 *funcWeight2 = new TF1("funcWeight2", FuncWeight2, 1.9, 4.5, 5);
    funcWeight2 -> SetParameter(0, 1);
    funcWeight2 -> SetParameter(1, 1);
    funcWeight2 -> SetParameter(2, 1);
    funcWeight2 -> SetParameter(3, 1);
    funcWeight2 -> SetParameter(4, 1);

    histRatioMass -> Fit(funcWeight2, "R0");
    funcWeight2 -> SetLineColor(kViolet);

    

    TF1 *funcBkg1 = new TF1("funcBkg1", FuncBkg1, 1.9, 4.5, 10);
    funcBkg1 -> SetParameter(0, 3.5);
    funcBkg1 -> SetParameter(1, -28920.8);
    funcBkg1 -> SetParameter(2, 24892.4);
    funcBkg1 -> SetParameter(3, -4600.53);
    funcBkg1 -> SetParameter(4, 34536.4);
    funcBkg1 -> SetParameter(5, -14899.7);
    funcBkg1 -> SetParameter(6, 1628.73);
    funcBkg1 -> SetParameter(7, 53633.4);
    funcBkg1 -> SetParameter(8, -32966.3);
    funcBkg1 -> SetParameter(9, 5538.45);

    histMassProjPtCutRun3 -> Fit(funcBkg1, "R0");
    funcBkg1 -> SetLineColor(kRed);

    TF1 *funcBkg1Pol2a = new TF1("funcBkg1Pol2a", FuncPol2, 1.9, 4.5, 3);
    funcBkg1Pol2a -> FixParameter(0, funcBkg1 -> GetParameter(1));
    funcBkg1Pol2a -> FixParameter(1, funcBkg1 -> GetParameter(2));
    funcBkg1Pol2a -> FixParameter(2, funcBkg1 -> GetParameter(3));
    funcBkg1Pol2a -> SetLineColor(kOrange+7);
    funcBkg1Pol2a -> SetLineStyle(kDotted);

    TF1 *funcBkg1Pol2b = new TF1("funcBkg1Pol2b", FuncPol2, 1.9, 4.5, 3);
    funcBkg1Pol2b -> FixParameter(0, funcBkg1 -> GetParameter(4));
    funcBkg1Pol2b -> FixParameter(1, funcBkg1 -> GetParameter(5));
    funcBkg1Pol2b -> FixParameter(2, funcBkg1 -> GetParameter(6));
    funcBkg1Pol2b -> SetLineColor(kBlue+1);
    funcBkg1Pol2b -> SetLineStyle(kDashed);

    TF1 *funcBkg1Pol2c = new TF1("funcBkg1Pol2c", FuncPol2, 1.9, 4.5, 3);
    funcBkg1Pol2c -> FixParameter(0, funcBkg1 -> GetParameter(7));
    funcBkg1Pol2c -> FixParameter(1, funcBkg1 -> GetParameter(8));
    funcBkg1Pol2c -> FixParameter(2, funcBkg1 -> GetParameter(9));
    funcBkg1Pol2c -> SetLineColor(kGreen+1);
    funcBkg1Pol2c -> SetLineStyle(kDashed);

    /*TF1 *funcBkg2 = new TF1("funcBkg2", FuncBkg2, 1.9, 4.5, 9);
    funcBkg2 -> SetParameter(0, 3.5);
    //funcBkg2 -> SetParLimits(0, 3, 4);
    funcBkg2 -> SetParameter(1, 0.2);
    //funcBkg2 -> SetParLimits(1, 0, 1);
    funcBkg2 -> SetParameter(2, 1000);
    //funcBkg2 -> SetParLimits(2, 0, 1000000);
    funcBkg2 -> SetParameter(3, 7533.85);
    funcBkg2 -> SetParameter(4, -2467.73);
    funcBkg2 -> SetParameter(5, 193.495);
    funcBkg2 -> SetParameter(6, 12558.8);
    funcBkg2 -> SetParameter(7, -6658.36);
    funcBkg2 -> SetParameter(8, 986.314);*/

    TF1 *funcBkg2 = new TF1("funcBkg2", FuncBkg2, 1.9, 4.5, 8);
    funcBkg2 -> SetParameter(0, 3.5);
    funcBkg2 -> SetParameter(1, 0.2);
    funcBkg2 -> SetParameter(2, 7533.85);
    funcBkg2 -> SetParameter(3, -2467.73);
    funcBkg2 -> SetParameter(4, 193.495);
    funcBkg2 -> SetParameter(5, 12558.8);
    funcBkg2 -> SetParameter(6, -6658.36);
    funcBkg2 -> SetParameter(7, 986.314);

    histMassProjPtCutRun3 -> Fit(funcBkg2, "R0");
    funcBkg2 -> SetLineColor(kRed);

    TF1 *funcBkg2Gaus = new TF1("funcBkg2Gaus", FuncGaus, 1.9, 4.5, 3);
    funcBkg2Gaus -> FixParameter(0, funcBkg2 -> GetParameter(0));
    funcBkg2Gaus -> FixParameter(1, funcBkg2 -> GetParameter(1));
    funcBkg2Gaus -> FixParameter(2, funcBkg2 -> GetParameter(2));
    funcBkg2Gaus -> SetLineColor(kOrange+7);
    funcBkg2Gaus -> SetLineStyle(kDotted);

    TF1 *funcBkg2Pol2a = new TF1("funcBkg2Pol2a", FuncPol2, 1.9, 4.5, 3);
    funcBkg2Pol2a -> FixParameter(0, funcBkg2 -> GetParameter(3));
    funcBkg2Pol2a -> FixParameter(1, funcBkg2 -> GetParameter(4));
    funcBkg2Pol2a -> FixParameter(2, funcBkg2 -> GetParameter(5));
    funcBkg2Pol2a -> SetLineColor(kBlue+1);
    funcBkg2Pol2a -> SetLineStyle(kDashed);

    TF1 *funcBkg2Pol2b = new TF1("funcBkg2Pol2b", FuncPol2, 1.9, 4.5, 3);
    funcBkg2Pol2b -> FixParameter(0, funcBkg2 -> GetParameter(6));
    funcBkg2Pol2b -> FixParameter(1, funcBkg2 -> GetParameter(7));
    funcBkg2Pol2b -> FixParameter(2, funcBkg2 -> GetParameter(8));
    funcBkg2Pol2b -> SetLineColor(kGreen+1);
    funcBkg2Pol2b -> SetLineStyle(kDashed);

    TF1 *funcBkg3 = new TF1("funcBkg3", FuncBkg3, 1.9, 5., 7);
    funcBkg3 -> SetParameter(0, funcWeight1 -> GetParameter(0));
    funcBkg3 -> SetParameter(1, funcWeight1 -> GetParameter(1));
    funcBkg3 -> SetParameter(2, funcWeight1 -> GetParameter(2));
    funcBkg3 -> SetParameter(3, funcWeight1 -> GetParameter(3));
    funcBkg3 -> SetParameter(4, 5.16848e+06);
    funcBkg3 -> SetParameter(5, -5.59407);
    funcBkg3 -> SetParameter(6, 2.52076);

    histMassProjPtCutRun3 -> Fit(funcBkg3, "R0");
    funcBkg3 -> SetLineColor(kMagenta);

    TF1 *funcBkg4 = new TF1("funcBkg4", FuncBkg4, 1.9, 4.5, 9);
    funcBkg4 -> FixParameter(0, funcWeight2 -> GetParameter(0));
    funcBkg4 -> FixParameter(1, funcWeight2 -> GetParameter(1));
    funcBkg4 -> FixParameter(2, funcWeight2 -> GetParameter(2));
    funcBkg4 -> FixParameter(3, funcWeight2 -> GetParameter(3));
    funcBkg4 -> FixParameter(4, funcWeight2 -> GetParameter(4));
    funcBkg4 -> SetParameter(4, 1);
    funcBkg4 -> SetParameter(5, 1);
    funcBkg4 -> SetParameter(6, 1);
    funcBkg4 -> SetParameter(7, 1);
    funcBkg4 -> SetParameter(8, 1);

    histMassProjPtCutRun3 -> Fit(funcBkg4, "R0");
    funcBkg4 -> SetLineColor(kViolet);


    TCanvas *canvasComp = new TCanvas("canvasComp", "", 1200, 1200);
    canvasComp -> Divide(2, 2);

    canvasComp -> cd(1);
    gPad -> SetLogy(true);
    histMassProj -> Draw("EP");
    histMassProjPtCutRun3 -> Draw("EP SAME");

    canvasComp -> cd(2);
    gPad -> SetLogy(true);
    histMassProjPtCutRun3 -> Draw("EP");
    //funcBkg1 -> Draw("SAME");
    //funcBkg1Pol2a -> Draw("SAME");
    //funcBkg1Pol2b -> Draw("SAME");
    //funcBkg1Pol2c -> Draw("SAME");
    funcBkg3 -> Draw("SAME");
    //funcBkg4 -> Draw("SAME");

    canvasComp -> cd(3);
    gPad -> SetLogy(true);
    histRatioMass -> Draw("EP");
    funcWeight1 -> Draw("SAME");
    funcWeight2 -> Draw("SAME");

    canvasComp -> cd(4);
    gPad -> SetLogy(true);
    histMassProjPtCutRun3 -> Draw("EP");
    //funcBkg2 -> Draw("SAME");
    //funcBkg2Gaus -> Draw("SAME");
    //funcBkg2Pol2a -> Draw("SAME");
    //funcBkg2Pol2b -> Draw("SAME");
    //funcBkg3 -> Draw("SAME");
    funcBkg4 -> Draw("SAME");



    //TCanvas *canvasCompRatio = new TCanvas("canvasCompRatio", "", 800, 600);
    //histRatioMass -> Draw("EP");
}
///////////////////////////////////////////////////
double FuncBkg1(double *x, double *par) {
    double xx = x[0];
    double alpha = par[0];

    if (xx > alpha - 0.2 && xx < alpha + 0.2) {
        return par[1] + par[2]*xx + par[3]*xx*xx;
    } else {
        if (xx > alpha + 0.2) {
            return par[4] + par[5]*xx + par[6]*xx*xx;
        }
        if (xx < alpha - 0.2) {
            return par[7] + par[8]*xx + par[9]*xx*xx;
        }
    }
    return 0;
}
///////////////////////////////////////////////////
/*double FuncBkg2(double *x, double *par) {
    double xx = x[0];
    double mean = par[0];
    double sigma = par[1];

    if (xx > mean - 1*sigma && xx <= mean + 1*sigma) {
        double arg = - ((xx - mean) * (xx - mean)) / (2 * sigma * sigma);
        return par[2] * TMath::Exp(arg);
    } else {
        if (xx >= mean + 1*sigma) {
            return par[3] + par[4]*xx + par[5]*xx*xx;
        }
        if (xx < mean - 1*sigma) {
            return par[6] + par[7]*xx + par[8]*xx*xx;
        }
    }
    return 0;
}*/
///////////////////////////////////////////////////
/*double FuncBkg2(double *x, double *par) {
    double xx = x[0];
    double mean = par[0];
    double sigma = par[1];
    double arg = - ((xx - mean) * (xx - mean)) / (2 * sigma * sigma);

    if (xx >= mean) {
        return par[3] + par[4]*xx + par[5]*xx*xx + par[2] * TMath::Exp(arg);;
    }
    if (xx < mean) {
        return par[6] + par[7]*xx + par[8]*xx*xx + par[2] * TMath::Exp(arg);;
    }
    return 0;
}*/
///////////////////////////////////////////////////
double FuncBkg2(double *x, double *par) {
    double xx = x[0];
    double mean = par[0];
    double sigma = par[1];

    if (xx > (mean - sigma) && xx < (mean + sigma)) {
        return par[2] + par[3]*xx + par[4]*xx*xx + par[5] + par[6]*xx + par[7]*xx*xx;
    }
    if (xx >= (mean + sigma)) {
        return par[2] + par[3]*xx + par[4]*xx*xx;
    }
    if (xx < (mean - sigma)) {
        return par[5] + par[6]*xx + par[7]*xx*xx;
    }
    return 0;
}
///////////////////////////////////////////////////
double FuncBkg3(double *x, double *par) {
    double xx = x[0];

    double mean1 = par[1];
    double sigma1 = par[2];
    double arg1 = - ((xx - mean1) * (xx - mean1)) / (2 * sigma1 * sigma1);

    double mean2 = par[5];
    double sigma2 = par[6];
    double arg2 = - ((xx - mean2) * (xx - mean2)) / (2 * sigma2 * sigma2);

    return (par[4] * TMath::Exp(arg2)) / (par[3] + par[0] * TMath::Exp(arg1));
}
///////////////////////////////////////////////////
double FuncBkg4(double *x, double *par) {
    double xx = x[0];

    double mean1 = par[1];
    double sigma1 = par[2] + par[3] * ((xx - mean1) / mean1);
    double arg1 = - ((xx - mean1) * (xx - mean1)) / (2 * sigma1 * sigma1);

    double mean2 = par[6];
    double sigma2 = par[7] + par[8] * ((xx - mean2) / mean2);
    double arg2 = - (xx - mean2) * (xx - mean2) / (2. * sigma2 * sigma2);

    return (par[5] * TMath::Exp(arg2)) / (par[4] + par[0] * TMath::Exp(arg1));
}
///////////////////////////////////////////////////
double FuncPol2(double *x, double *par) {
    double xx = x[0];
    return par[0] + par[1]*xx + par[2]*xx*xx;
}
///////////////////////////////////////////////////
double FuncGaus(double *x, double *par) {
    double xx = x[0];
    double mean = par[0];
    double sigma = par[1];

    double arg = - ((xx - mean) * (xx - mean)) / (2 * sigma * sigma);
    return par[2] * TMath::Exp(arg);
}
///////////////////////////////////////////////////
double FuncVWG(double *x, double *par) {
    double xx = x[0];
    double mean = par[0];
    double sigma = par[1] + par[2] * ((xx - mean) / mean);
    double arg = - (xx - mean) * (xx - mean) / (2. * sigma * sigma);

    return par[3] * TMath::Exp(arg);
}
///////////////////////////////////////////////////
double FuncWeight1(double *x, double *par) {
    double xx = x[0];
    double mean = par[1];
    double sigma = par[2];

    double arg = - ((xx - mean) * (xx - mean)) / (2 * sigma * sigma);
    return par[3] + par[0] * TMath::Exp(arg);
}
///////////////////////////////////////////////////
double FuncWeight2(double *x, double *par) {
    double xx = x[0];
    double mean = par[1];
    double sigma = par[2] + par[3] * ((xx - mean) / mean);
    double arg = - (xx - mean) * (xx - mean) / (2. * sigma * sigma);

    return par[4] + par[0] * TMath::Exp(arg);
}
