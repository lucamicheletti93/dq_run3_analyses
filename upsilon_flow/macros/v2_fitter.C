#include "TProfile.h"
#include <fstream>
#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include <vector>
#include "TSystem.h"

Double_t FitFunctionBackgroundPol2(Double_t *x, Double_t *par);
Double_t FitFunctionBackgroundPol3(Double_t *x, Double_t *par);
Double_t alphaCB2VWG(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2CB2VWGPOL2(Double_t *x, Double_t *par);
void SetFitRejectRange(Double_t a, Double_t b);
Double_t BackgroundVWG(Double_t *x, Double_t *par);
Double_t CrystalBallExtended(Double_t *x,Double_t *par);
Double_t FitFunctionSignalCrystalBallExtended(Double_t *x,Double_t *par);
Double_t FitFunctionBackgroundVWG(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2CB2POL2EXPPOL2(Double_t *x, Double_t *par);
Double_t alphaCB2POL2EXP(Double_t*x, Double_t* par);
Double_t alphaCB2POL4EXP(Double_t*x, Double_t* par);
Double_t FitFunctionBackgroundPol2Exp(Double_t *x, Double_t *par);
Double_t FitFunctionBackgroundPol4Exp(Double_t *x, Double_t *par);
Double_t FitFunctionNA60New(Double_t *x,Double_t *par);
Double_t FitFunctionFlowS2NA60NEWPOL2EXPPOL2(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2NA60NEWVWGPOL2(Double_t *x, Double_t *par);
Double_t alphaNA60NEWVWG(Double_t*x, Double_t* par);
Double_t alphaNA60NEWPOL2EXP(Double_t*x, Double_t* par);
Double_t alphaNA60NEWPOL4EXP(Double_t*x, Double_t* par);
Double_t FitFunctionBackgroundPol2Cheb(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2CB2VWGPol3Cheb(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2NA60NEWPOL4EXPPOL2(Double_t *x, Double_t *par);
Double_t FitFunctionFlowS2CB2POL4EXPPOL2(Double_t *x, Double_t *par);
Double_t fitFunctionCB2VWG(Double_t *x, Double_t *par);
Double_t fitFunctionCB2Pol4Exp(Double_t *x, Double_t *par);
Double_t fitFunctionNA60NEWVWG(Double_t *x, Double_t *par);
Double_t fitFunctionNA60NEWPol4Exp(Double_t *x, Double_t *par);

double* DoFlowFit(TH1D *, TH1D *, TProfile *, double , double, double , double , bool , string , TFile *);
std::vector<TH1*> Uncertainities(std::vector<TH1*> histlist, vector<double> parameter, vector<double> parameter_er, double value[10][3]);

Double_t xmin= 2.3;
Double_t xmax= 4.8;

const int nSigFuncs = 2;
string sigFuncs[] = {"CB2", "NA60"};
const int nBkgFuncs = 2;
string bkgFuncs[] = {"VWG", "Pol4Exp"};
const int nBkgV2Funcs = 2;
string bkgV2Funcs[] = {"Pol2", "Cheb"};
const int nTailSets = 2;
string tailSets[] = {"data", "MC"};

string parNames_CB2_VWG[] = {"bkg","aa","bb","cc","sig_Jpsi","mean_Jpsi","width_Jpsi","a_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi"};
string parNames_CB2_Pol4Exp[] = {"bkg","aa","bb","cc","dd","ee","ff","sig_Jpsi","mean_Jpsi","width_Jpsi","a_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi"};
string parNames_NA60_VWG[] = {"bkg","aa","bb","cc","sig_Jpsi","mean_Jpsi","width_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi","f_Jpsi","g_Jpsi","h_Jpsi","a_Jpsi","e_Jpsi"};
string parNames_NA60_Pol4Exp[] = {"bkg","aa","bb","cc","dd","ee","ff","sig_Jpsi","mean_Jpsi","width_Jpsi","b_Jpsi","c_Jpsi","d_Jpsi","f_Jpsi","g_Jpsi","h_Jpsi","a_Jpsi","e_Jpsi"};

double minCentrBins[] = {10};
double maxCentrBins[] = {50};

const int nPtBins = 10;
double minPtBins[] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12};
double maxPtBins[] = {1, 2, 3, 4, 5, 6, 8, 10, 12, 15};
double ptBins[] = {0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15};
//double minPtBins[] = {3};
//double maxPtBins[] = {4};

const int nFitRanges = 3;
double minFitRanges[] = {2.3, 2.4, 2.5};
double maxFitRanges[] = {4.7, 4.6, 4.5};

int nTrials = 0;

double resolution = 1.28004;
string dirPath = "/Users/lucamicheletti/GITHUB/jpsi_flow_run3/macros/systematics";

void v2_fitter() {
  if (!gSystem -> AccessPathName(dirPath.c_str())) {
    std::cout << "The output directory already exists! " << std::endl;
  } else {
    int status = gSystem -> MakeDirectory(dirPath.c_str());
    if (status == 0) {
      std::cout << "Output directory created!" << std::endl;
    } else {
      std::cout << "Error in the creation of the utput directory" << std::endl;
        return;
    }
  }

  // Histograms for plotting the results
  TH1D *histStatJpsiV2 = new TH1D("histStatJpsiV2", "", 10, ptBins);
  histStatJpsiV2 -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
  histStatJpsiV2 -> GetYaxis() -> SetRangeUser(-0.05, 0.4);
  histStatJpsiV2 -> GetYaxis() -> SetTitle("#it{v}_{2}");
  histStatJpsiV2 -> SetMarkerStyle(20);
  histStatJpsiV2 -> SetMarkerColor(kRed+1);
  histStatJpsiV2 -> SetLineColor(kRed+1);

  TH1D *histSystJpsiV2 = new TH1D("histSystJpsiV2", "", 10, ptBins);
  histSystJpsiV2 -> SetMarkerStyle(20);
  histSystJpsiV2 -> SetMarkerColor(kRed+1);
  histSystJpsiV2 -> SetLineColor(kRed+1);
  histSystJpsiV2 -> SetFillStyle(0);

  for (int iPt = 0;iPt < nPtBins;iPt++) {
    vector<double> jpsiV2s;
    vector<double> errJpsiV2s;
    vector<double> chi2Ndfs;
    vector<string> trialSetNames;

    double minPtBin = minPtBins[iPt];
    double maxPtBin = maxPtBins[iPt];
    nTrials = 0;

    TFile *fIn = new TFile("/Users/lucamicheletti/cernbox/JPSI/Jpsi_flow/data/pass3/LHC23_full/Histograms_Fullpass3matchedMchMid_centr_Add_Pt_bin_10_50.root");
    TFile *fOut = new TFile(Form("systematics/fitResults_Pt_%2.1f_%2.1f.root", minPtBin, maxPtBin), "RECREATE");

    for (int iSigFunc = 0;iSigFunc < nSigFuncs;iSigFunc++) { // loop on mass signal funcs
      for (int iBkgFunc = 0;iBkgFunc < nBkgFuncs;iBkgFunc++) { // loop on mass background funcs
        for (int iTailSet = 0;iTailSet < nTailSets;iTailSet++) { // loop on tail parameters
          string sigFunc = sigFuncs[iSigFunc];
          string bkgFunc = bkgFuncs[iBkgFunc];
          string tailSet = tailSets[iTailSet];

          if (sigFunc == "NA60" && tailSet == "data") {
            continue;
          }

          TFile *fInMassFit = new TFile(Form("/Users/lucamicheletti/GITHUB/dq_fitter/analysis/LHC23_pass3_full/centrality_10_50/pt_dependence/Pt_%2.1f_%2.1f/multi_trial_%s_%s_%s_tails.root", minPtBin, maxPtBin, sigFunc.c_str(), bkgFunc.c_str(), tailSet.c_str()));

          for (int iFitRange = 0;iFitRange < nFitRanges;iFitRange++) {
            double minFitRange = minFitRanges[iFitRange];
            double maxFitRange = maxFitRanges[iFitRange];

            TH1D *histMassFitPars = (TH1D*) fInMassFit -> Get(Form("fit_results_%s_%s__%2.1f_%2.1f_histMassSEPM_%1.0f_%1.0f__10_50", sigFunc.c_str(), bkgFunc.c_str(), minFitRange, maxFitRange, minPtBin, maxPtBin));
            TH1D *histMass = (TH1D*) fIn -> Get(Form("histMassSEPM_%1.0f_%1.0f__10_50", minPtBin, maxPtBin));
            TProfile *profV2 = (TProfile*) fIn -> Get(Form("histV2SEPM_%1.0f_%1.0f__10_50", minPtBin, maxPtBin));

            // Plot options
            histMass -> GetXaxis() -> SetRangeUser(minFitRange, maxFitRange);
            histMass -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c^{2}})");
            histMass -> SetMarkerStyle(20);
            histMass -> SetMarkerSize(0.8);
            histMass -> SetMarkerColor(kBlack);
            histMass -> SetLineColor(kBlack);

            profV2 -> SetTitle("");
            profV2 -> GetXaxis() -> SetRangeUser(minFitRange, maxFitRange);
            profV2 -> GetXaxis() -> SetTitle("#it{M}_{#mu#mu} (GeV/#it{c^{2}})");
            profV2 -> SetMarkerStyle(20);
            profV2 -> SetMarkerSize(0.8);
            profV2 -> SetMarkerColor(kBlack);
            profV2 -> SetLineColor(kBlack);

            for (int iBkgV2Func = 0;iBkgV2Func < nBkgV2Funcs;iBkgV2Func++) { // loop on v2 background types
              string bkgV2Func = bkgV2Funcs[iBkgV2Func];
              double *results; 
              if (iSigFunc == 0 && iBkgFunc == 0 && iTailSet == 0 && iFitRange == 0 && iBkgV2Func == 0) {
                results = DoFlowFit(histMassFitPars, histMass, profV2, minPtBin, maxPtBin, minFitRange, maxFitRange, kTRUE, Form("%s_%s_%s_%s_v2bkg", sigFunc.c_str(), bkgFunc.c_str(), tailSet.c_str(), bkgV2Func.c_str()), fOut);
              } else {
                results = DoFlowFit(histMassFitPars, histMass, profV2, minPtBin, maxPtBin, minFitRange, maxFitRange, kFALSE, Form("%s_%s_%s_%s_v2bkg", sigFunc.c_str(), bkgFunc.c_str(), tailSet.c_str(), bkgV2Func.c_str()), fOut);
              }
              

              // WARNING! v2 from the fit is scaled by the RESOLUTION
              jpsiV2s.push_back(results[0] * resolution);
              errJpsiV2s.push_back(results[1] * resolution);
              chi2Ndfs.push_back(results[2]);
              trialSetNames.push_back(Form("%s + %s + %s tails, %2.1f - %2.1f, %s[v2 bkg]", sigFunc.c_str(), bkgFunc.c_str(), tailSet.c_str(), minFitRange, maxFitRange, bkgV2Func.c_str()));

              // Update the number of trials
              nTrials++;
            }
          }
        }
      }
    }

    double sumJpsiV2s = 0; 
    double statErrJpsiV2s = 0;

    TH1D *histSyst = new TH1D("histSyst", "", nTrials, 0, nTrials);
    histSyst -> SetMarkerStyle(20);
    histSyst -> SetMarkerColor(kBlack);
    histSyst -> SetLineColor(kBlack);

    TH1D *histChi2Ndf = new TH1D("histChi2Ndf", "", nTrials, 0, nTrials);
    histChi2Ndf -> SetMarkerStyle(20);
    histChi2Ndf -> SetMarkerColor(kBlack);
    histChi2Ndf -> SetLineColor(kBlack);

    histSyst -> GetYaxis() -> SetRangeUser(0, 0.2);
    for (int iBin = 0;iBin < nTrials;iBin++) {
      //cout << iBin << ") " << jpsiV2s[iBin] << " +/- " << errJpsiV2s[iBin] << " -> " << trialSetNames[iBin].c_str() << endl;
      histSyst -> SetBinContent(iBin+1, jpsiV2s[iBin]);
      histSyst -> SetBinError(iBin+1, errJpsiV2s[iBin]);
      histSyst -> GetXaxis() -> SetBinLabel(iBin+1, trialSetNames[iBin].c_str());

      sumJpsiV2s = sumJpsiV2s + jpsiV2s[iBin];
	    statErrJpsiV2s = statErrJpsiV2s + errJpsiV2s[iBin];

      histChi2Ndf -> SetBinContent(iBin+1, chi2Ndfs[iBin]);
      histChi2Ndf -> GetXaxis() -> SetBinLabel(iBin+1, trialSetNames[iBin].c_str());
    }
    double meanJpsiV2 = sumJpsiV2s / nTrials;
    double meanStatErrJpsiV2 = statErrJpsiV2s / nTrials;

    double sumSyst = 0;
    for (int iBin = 0;iBin < nTrials;iBin++) {
      double dev = (jpsiV2s[iBin] - meanJpsiV2);
	    sumSyst = sumSyst + dev * dev;
    }
    double meanSystErrJpsiV2 = TMath::Sqrt(sumSyst / nTrials);

    histStatJpsiV2 -> SetBinContent(iPt+1, meanJpsiV2);
    histStatJpsiV2 -> SetBinError(iPt+1, meanStatErrJpsiV2);
    histSystJpsiV2 -> SetBinContent(iPt+1, meanJpsiV2);
    histSystJpsiV2 -> SetBinError(iPt+1, meanSystErrJpsiV2);

    double meanStatErrJpsiV2Perc = (meanStatErrJpsiV2 / meanJpsiV2) * 100;
    double meanSystErrJpsiV2Perc = (meanSystErrJpsiV2 / meanJpsiV2) * 100;

    TLine *lineMean = new TLine(0, meanJpsiV2, nTrials, meanJpsiV2);
    lineMean -> SetLineStyle(1);
    lineMean -> SetLineColor(kBlue);
    lineMean -> SetLineWidth(2);

    TLine *lineSystUp = new TLine(0, meanJpsiV2 + meanSystErrJpsiV2, nTrials, meanJpsiV2 + meanSystErrJpsiV2);
    lineSystUp-> SetLineStyle(2);
    lineSystUp-> SetLineColor(kBlue);
    lineSystUp-> SetLineWidth(2);

    TLine *lineSystLw = new TLine(0, meanJpsiV2 - meanSystErrJpsiV2, nTrials, meanJpsiV2 - meanSystErrJpsiV2);
    lineSystLw -> SetLineStyle(2);
    lineSystLw -> SetLineColor(kBlue);
    lineSystLw -> SetLineWidth(2);

    TLine *lineStatUp = new TLine(0, meanJpsiV2 + meanStatErrJpsiV2, nTrials, meanJpsiV2 + meanStatErrJpsiV2);
    lineStatUp-> SetLineStyle(2);
    lineStatUp-> SetLineColor(kGray+1);
    lineStatUp-> SetLineWidth(2);

    TLine *lineStatLw = new TLine(0, meanJpsiV2 - meanStatErrJpsiV2, nTrials, meanJpsiV2 - meanStatErrJpsiV2);
    lineStatLw -> SetLineStyle(2);
    lineStatLw -> SetLineColor(kGray+1);
    lineStatLw -> SetLineWidth(2);

    TCanvas *canvasSyst = new TCanvas("canvasSyst", "", 800, 900);
    canvasSyst -> SetTopMargin(0.05);
    canvasSyst -> SetBottomMargin(0.5);
    histSyst -> SetStats(0);
    histSyst -> GetXaxis() -> LabelsOption("v");
    if (meanJpsiV2 > 0) {
      histSyst -> GetYaxis() -> SetRangeUser(meanJpsiV2 -  (meanJpsiV2 / 2.), meanJpsiV2 + (meanJpsiV2 / 2.));
    } else {
      histSyst -> GetYaxis() -> SetRangeUser(-0.05, 0.05);
    }
    
    histSyst -> Draw();
    lineMean -> Draw();
    lineSystUp-> Draw();
    lineSystLw -> Draw();
    lineStatUp-> Draw();
    lineStatLw -> Draw();

    //TPaveText *display1 = new TPaveText(0.22, 0.80, 0.80, 0.85, "blNDC");
    TPaveText *display1 = new TPaveText(0.20, 0.80, 0.80, 0.85, "blNDC");
    display1 -> SetTextFont(42);
    display1 -> SetTextSize(0.036);
    display1 -> SetTextColor(kBlack);
    display1 -> SetBorderSize(0);
    display1 -> SetFillColor(0);

    //TText *text1 = display1 -> AddText(Form(" v_{2}^{J/#psi} [%1.0f < #it{p}_{T} < %1.0f GeV/#it{c}] = %5.4f #pm %5.4f (%4.3f %%) #pm %5.4f (%4.3f %%)", minPtBin, maxPtBin, meanJpsiV2, meanStatErrJpsiV2, meanStatErrJpsiV2Perc, meanSystErrJpsiV2, meanSystErrJpsiV2Perc));
    TText *text1 = display1 -> AddText(Form(" v_{2}^{J/#psi} [%1.0f < #it{p}_{T} < %1.0f GeV/#it{c}] = %5.4f #pm %5.4f #pm %5.4f", minPtBin, maxPtBin, meanJpsiV2, meanStatErrJpsiV2, meanSystErrJpsiV2));
    display1 -> Draw("same");

    canvasSyst -> Write();
    histChi2Ndf -> Write();
    canvasSyst -> SaveAs(Form("%s/systematic_plot_Pt_%1.0f_%1.0f.pdf", dirPath.c_str(), minPtBin, maxPtBin));

    fOut -> Close();
  }

  TCanvas *canvasJpsiV2 = new TCanvas("canvasJpsiV2", "", 800, 600);
  histStatJpsiV2 -> Draw("EP");
  histSystJpsiV2 -> Draw("E2P SAME");

  cout << "-------------------------" << endl;
  cout << "x_min x_max val stat syst" << endl;
  for (int iPt = 0;iPt < nPtBins;iPt++) {
    Printf("%3.2f %3.2f %6.5f %6.5f %6.5f ", minPtBins[iPt], maxPtBins[iPt], histStatJpsiV2 -> GetBinContent(iPt+1), histStatJpsiV2 -> GetBinError(iPt+1), histSystJpsiV2 -> GetBinError(iPt+1));
  }
}
////////////////////////////////////////////////////////////
double* DoFlowFit(TH1D *histMassFitPars, TH1D *histMass, TProfile *profV2, double minPtBin, double maxPtBin, double minFitRange, double maxFitRange, bool debug, string fitFuncs, TFile *fOut) {
  // Get signal shape from previous fit
  vector<double> fitPars;
  fitPars.clear();
  double* results = new double[3];

  TF1 *funcMassBkg;
  TF1 *funcMassSig;
  TF1 *funcMassSigBkg;
  TF1 *funcFlowSigBkg;

  //--------------------------------------------------//
  //                     CB2+VWG                      //
  //--------------------------------------------------//
  if (fitFuncs.find("CB2_VWG") != std::string::npos) {
    for (int i = 0; i < 11; i++) {
      fitPars.push_back(histMassFitPars -> GetBinContent(histMassFitPars -> GetXaxis() -> FindBin(Form("%s", parNames_CB2_VWG[i].data()))));
    }

    // STARTING THE FIT TO THE INVARIANT MASS
    funcMassSigBkg = new TF1("funcMassSigBkg", fitFunctionCB2VWG, minFitRange, maxFitRange, 11);
    for (int i = 1; i < 4; i++) {
      funcMassSigBkg->FixParameter(i, fitPars[i]);
    }
    for (int i = 5; i < 11; i++) {
      funcMassSigBkg->FixParameter(i, fitPars[i]);
    }

    for (int iter = 0;iter < 100;iter++) {
      if (iter == 0) {
        funcMassSigBkg -> SetParameter(0, 1e6); // Bkg normalization
        funcMassSigBkg -> SetParameter(4, 300000); // Jpsi normalization
      } else {
        funcMassSigBkg -> SetParLimits(0, 0, 1e10);
        funcMassSigBkg -> SetParLimits(4, 0, 1e10);
        funcMassSigBkg -> SetParameter(0, funcMassSigBkg -> GetParameter(0)); // Bkg normalization
        funcMassSigBkg -> SetParameter(4, funcMassSigBkg -> GetParameter(4)); // Jpsi normalization
      }
      TFitResultPtr fitResult = histMass -> Fit(funcMassSigBkg, "RL0S");
      int fitStatus = fitResult;  // Fit status
      int covMatrixStatus = fitResult -> CovMatrixStatus(); // cov. matrix status

      if (fitStatus == 0 && covMatrixStatus == 3) {
        break;
      }
    }

    funcMassBkg = new TF1("funcMassBkg", FitFunctionBackgroundVWG, minFitRange, maxFitRange, 4);
    funcMassBkg -> SetLineColor(kGray+1);
    funcMassBkg -> SetLineStyle(kDashed);
    for (int iPar = 0;iPar < 4;iPar++) {
      funcMassBkg -> FixParameter(iPar, funcMassSigBkg -> GetParameter(iPar));
    }

    funcMassSig = new TF1("funcMassSig", FitFunctionSignalCrystalBallExtended, minFitRange, maxFitRange, 7);
    funcMassSig -> SetLineColor(kAzure+2);
    funcMassSig -> SetLineStyle(kSolid);
    for (int iPar = 0;iPar < 7;iPar++) {
      funcMassSig -> FixParameter(iPar, funcMassSigBkg -> GetParameter(iPar+4));
    }
        
    funcFlowSigBkg  = new TF1("funcFlowSigBkg ", FitFunctionFlowS2CB2VWGPOL2, minFitRange, maxFitRange, 18);
    funcFlowSigBkg -> SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","alJPsi","nlJPsi","auJPsi","nuJPsi");
    funcFlowSigBkg -> SetParName(11, "kPsiP");
    funcFlowSigBkg -> SetParName(12, "v_{2} JPsi");
    funcFlowSigBkg -> SetParName(13, "v_{2} BG0");
    funcFlowSigBkg -> SetParName(14, "v_{2} BG1");
    funcFlowSigBkg -> SetParName(15, "v_{2} BG2");
    funcFlowSigBkg -> SetParName(16, "v_{2} BG3");
    funcFlowSigBkg -> SetParName(17, "type");
    funcFlowSigBkg -> SetLineColor(kBlue);

    for (int iPar = 0;iPar < 11;iPar++) {
      funcFlowSigBkg -> FixParameter(iPar, funcMassSigBkg -> GetParameter(iPar));
    }
    funcFlowSigBkg -> FixParameter(11, 0);
    if (fitFuncs.find("Pol2_v2bkg") != std::string::npos) {
      funcFlowSigBkg -> SetParameter(12, 0.072);
      funcFlowSigBkg -> SetParameter(13, 0.09);
      funcFlowSigBkg -> SetParameter(14, -0.0252894);
      funcFlowSigBkg -> SetParameter(15, 0.0022179);
      funcFlowSigBkg -> FixParameter(16, 999); // Extra parameter for 3 degree cheby
      funcFlowSigBkg -> FixParameter(17, 0); // v2 bkg is a Pol2
    } else {
      funcFlowSigBkg -> SetParameter(12, 0.072);
      funcFlowSigBkg -> SetParameter(13, 0.0323879);
      funcFlowSigBkg -> SetParameter(14, -0.0117933);
      funcFlowSigBkg -> SetParameter(15, 0.000888077);
      funcFlowSigBkg -> SetParameter(16, -0.000522725); // Extra parameter for 3 degree cheby
      funcFlowSigBkg -> FixParameter(17, 1); // v2 bkg is a Cheby
    }
    profV2 -> Fit (funcFlowSigBkg, "RI"); // Get the normalization from the fit

    results[0] = funcFlowSigBkg -> GetParameter(12);
    results[1] = funcFlowSigBkg -> GetParError(12);
    results[2] = (double) funcFlowSigBkg -> GetChisquare() / (double) funcFlowSigBkg -> GetNDF();
  }

  //--------------------------------------------------//
  //                   CB2+Pol4Exp                    //
  //--------------------------------------------------//
  if (fitFuncs.find("CB2_Pol4Exp") != std::string::npos) {
    for (int i = 0; i < 14; i++) {
      fitPars.push_back(histMassFitPars -> GetBinContent(histMassFitPars -> GetXaxis() -> FindBin(Form("%s",parNames_CB2_Pol4Exp[i].data()))));
    }

    // STARTING THE FIT TO THE INVARIANT MASS
    funcMassSigBkg = new TF1("funcMassSigBkg", fitFunctionCB2Pol4Exp, minFitRange, maxFitRange, 14);
    for (int i = 1; i < 7; i++) {
      funcMassSigBkg->FixParameter(i, fitPars[i]);
    }
    for (int i = 8; i < 14; i++) {
      funcMassSigBkg->FixParameter(i, fitPars[i]);
    }

    for (int iter = 0;iter < 100;iter++) {
      if (iter == 0) {
        funcMassSigBkg -> SetParameter(0, 1e5); // Bkg normalization
        funcMassSigBkg -> SetParameter(7, 100000); // Jpsi normalization
      } else {
        funcMassSigBkg -> SetParLimits(0, 0.1, 1e9);
        funcMassSigBkg -> SetParLimits(7, 0.1, 1e9);
        funcMassSigBkg -> SetParameter(0, funcMassSigBkg -> GetParameter(0)); // Bkg normalization
        funcMassSigBkg -> SetParameter(7, funcMassSigBkg -> GetParameter(7)); // Jpsi normalization
      }
      TFitResultPtr fitResult = histMass -> Fit(funcMassSigBkg, "RL0S");
      int fitStatus = fitResult;  // Fit status
      int covMatrixStatus = fitResult -> CovMatrixStatus(); // cov. matrix status

      if (fitStatus == 0 && covMatrixStatus == 3) {
        break;
      }
    }

    funcMassBkg = new TF1("funcMassBkg", FitFunctionBackgroundPol4Exp, minFitRange, maxFitRange, 7);
    funcMassBkg -> SetLineColor(kGray+1);
    funcMassBkg -> SetLineStyle(kDashed);
    for (int iPar = 0;iPar < 7;iPar++) {
      funcMassBkg -> FixParameter(iPar, funcMassSigBkg -> GetParameter(iPar));
    }

    funcMassSig = new TF1("funcMassSig", FitFunctionSignalCrystalBallExtended, minFitRange, maxFitRange, 7);
    funcMassSig -> SetLineColor(kAzure+2);
    funcMassSig -> SetLineStyle(kSolid);
    for (int iPar = 0;iPar < 7;iPar++) {
      funcMassSig -> FixParameter(iPar, funcMassSigBkg -> GetParameter(iPar+7));
    }
        
    funcFlowSigBkg  = new TF1("funcFlowSigBkg ", FitFunctionFlowS2CB2POL4EXPPOL2, minFitRange, maxFitRange, 21);
    funcFlowSigBkg -> SetParNames("pol0","pol1","pol2","pol3","exp1","exp2","exp3","kJPsi","mJPsi","sJPsi","alJPsi");
    funcFlowSigBkg -> SetParName(11,"nlJPsi");
    funcFlowSigBkg -> SetParName(12,"auJPsi");
    funcFlowSigBkg -> SetParName(13,"nuJPsi");
    funcFlowSigBkg -> SetParName(14,"kPsiP");
    funcFlowSigBkg -> SetParName(15, "v_{2} JPsi");
    funcFlowSigBkg -> SetParName(16, "v_{2} BG0");
    funcFlowSigBkg -> SetParName(17, "v_{2} BG1");
    funcFlowSigBkg -> SetParName(18, "v_{2} BG2");
    funcFlowSigBkg -> SetParName(19, "v_{2} BG3");
    funcFlowSigBkg -> SetParName(20, "type");
    funcFlowSigBkg -> SetLineColor(kBlue);

    for (int iPar = 0;iPar < 14;iPar++) {
      funcFlowSigBkg -> FixParameter(iPar, funcMassSigBkg -> GetParameter(iPar));
    }
    if (fitFuncs.find("Pol2_v2bkg") != std::string::npos) {
      funcFlowSigBkg -> FixParameter(14, 0);
      funcFlowSigBkg -> SetParameter(15, 0.072);
      funcFlowSigBkg -> SetParameter(16, 0.09);
      funcFlowSigBkg -> SetParameter(17, -0.0252894);
      funcFlowSigBkg -> SetParameter(18, 0.0022179);
      funcFlowSigBkg -> FixParameter(19, 999); // Extra parameter for 3 degree cheby
      funcFlowSigBkg -> FixParameter(20, 0); // v2 bkg is a Pol2
    } else {
      funcFlowSigBkg -> FixParameter(14, 0);
      funcFlowSigBkg -> SetParameter(15, 0.072);
      funcFlowSigBkg -> SetParameter(16, 0.0323879);
      funcFlowSigBkg -> SetParameter(17, -0.0117933);
      funcFlowSigBkg -> SetParameter(18, 0.000888077);
      funcFlowSigBkg -> SetParameter(19, -0.000522725); // Extra parameter for 3 degree cheby
      funcFlowSigBkg -> FixParameter(20, 1); // v2 bkg is a Cheby
    }
    profV2 -> Fit (funcFlowSigBkg, "RI"); // Get the normalization from the fit

    results[0] = funcFlowSigBkg -> GetParameter(15);
    results[1] = funcFlowSigBkg -> GetParError(15);
    results[2] = (double) funcFlowSigBkg -> GetChisquare() / (double) funcFlowSigBkg -> GetNDF();
  }

  //--------------------------------------------------//
  //                    NA60+VWG                      //
  //--------------------------------------------------//
  if (fitFuncs.find("NA60_VWG") != std::string::npos) {
    for (int i = 0; i < 15; i++) {
      fitPars.push_back(histMassFitPars -> GetBinContent(histMassFitPars -> GetXaxis() -> FindBin(Form("%s",parNames_NA60_VWG[i].data()))));
    }

    // STARTING THE FIT TO THE INVARIANT MASS
    funcMassSigBkg = new TF1("funcMassSigBkg", fitFunctionNA60NEWVWG, minFitRange, maxFitRange, 15);
    for (int i = 1; i < 4; i++) {
      funcMassSigBkg->FixParameter(i, fitPars[i]);
    }
    for (int i = 5; i < 15; i++) {
      funcMassSigBkg->FixParameter(i, fitPars[i]);
    }

    for (int iter = 0;iter < 100;iter++) {
      if (iter == 0) {
        funcMassSigBkg -> SetParameter(0, 1e6); // Bkg normalization
        funcMassSigBkg -> SetParameter(4, 300000); // Jpsi normalization
      } else {
        funcMassSigBkg -> SetParLimits(0, 0, 1e10);
        funcMassSigBkg -> SetParLimits(4, 0, 1e10);
        funcMassSigBkg -> SetParameter(0, funcMassSigBkg -> GetParameter(0)); // Bkg normalization
        funcMassSigBkg -> SetParameter(4, funcMassSigBkg -> GetParameter(7)); // Jpsi normalization
      }
      TFitResultPtr fitResult = histMass -> Fit(funcMassSigBkg, "RL0S");
      int fitStatus = fitResult;  // Fit status
      int covMatrixStatus = fitResult -> CovMatrixStatus(); // cov. matrix status

      if (fitStatus == 0 && covMatrixStatus == 3) {
        break;
      }
    }

    funcMassBkg = new TF1("funcMassBkg", FitFunctionBackgroundVWG, minFitRange, maxFitRange, 4);
    funcMassBkg -> SetLineColor(kGray+1);
    funcMassBkg -> SetLineStyle(kDashed);
    for (int iPar = 0;iPar < 7;iPar++) {
      funcMassBkg -> FixParameter(iPar, funcMassSigBkg -> GetParameter(iPar));
    }

    funcMassSig = new TF1("funcMassSig", FitFunctionNA60New, minFitRange, maxFitRange, 11);
    funcMassSig -> SetLineColor(kAzure+2);
    funcMassSig -> SetLineStyle(kSolid);
    for (int iPar = 0;iPar < 11;iPar++) {
      funcMassSig -> FixParameter(iPar, funcMassSigBkg -> GetParameter(iPar+7));
    }
          
    funcFlowSigBkg  = new TF1("funcFlowSigBkg ", FitFunctionFlowS2NA60NEWVWGPOL2, minFitRange, maxFitRange, 22);
    funcFlowSigBkg -> SetParNames("kVWG","mVWG","sVWG1","sVWG2","kJPsi","mJPsi","sJPsi","p1LJPsi","p2LJPsi","p3LJPsi","p1RJPsi");
    funcFlowSigBkg -> SetParName(11,"p2RJPsi");
    funcFlowSigBkg -> SetParName(12,"p3RJPsi");
    funcFlowSigBkg -> SetParName(13,"aLJPsi");
    funcFlowSigBkg -> SetParName(14,"aRJPsi");
    funcFlowSigBkg -> SetParName(15,"kPsiP");
    funcFlowSigBkg -> SetParName(16,"v_{2} JPsi");
    funcFlowSigBkg -> SetParName(17,"v_{2} BG0");
    funcFlowSigBkg -> SetParName(18,"v_{2} BG1");
    funcFlowSigBkg -> SetParName(19,"v_{2} BG2");
    funcFlowSigBkg -> SetParName(20,"v_{2} BG3");
    funcFlowSigBkg -> SetParName(21, "type");
    funcFlowSigBkg -> SetLineColor(kBlue);

    for (int iPar = 0;iPar < 15;iPar++) {
      funcFlowSigBkg -> FixParameter(iPar, funcMassSigBkg -> GetParameter(iPar));
    }
    if (fitFuncs.find("Pol2_v2bkg") != std::string::npos) {
      funcFlowSigBkg -> FixParameter(15, 0);
      funcFlowSigBkg -> SetParameter(16, 0.072);
      funcFlowSigBkg -> SetParameter(17, 0.09);
      funcFlowSigBkg -> SetParameter(18, -0.0252894);
      funcFlowSigBkg -> SetParameter(19, 0.0022179);
      funcFlowSigBkg -> FixParameter(20, 999); // Extra parameter for 3 degree cheby
      funcFlowSigBkg -> FixParameter(21, 0); // v2 bkg is a Pol2
    } else {
      funcFlowSigBkg -> FixParameter(15, 0);
      funcFlowSigBkg -> SetParameter(16, 0.072);
      funcFlowSigBkg -> SetParameter(17, 0.0323879);
      funcFlowSigBkg -> SetParameter(18, -0.0117933);
      funcFlowSigBkg -> SetParameter(19, 0.000888077);
      funcFlowSigBkg -> SetParameter(20, -0.000522725); // Extra parameter for 3 degree cheby
      funcFlowSigBkg -> FixParameter(21, 1); // v2 bkg is a Cheby
    }
    profV2 -> Fit (funcFlowSigBkg, "RI"); // Get the normalization from the fit

    results[0] = funcFlowSigBkg -> GetParameter(16);
    results[1] = funcFlowSigBkg -> GetParError(16);
    results[2] = (double) funcFlowSigBkg -> GetChisquare() / (double) funcFlowSigBkg -> GetNDF();
  }


  //--------------------------------------------------//
  //                   NA60+Pol4Exp                    //
  //--------------------------------------------------//
  if (fitFuncs.find("NA60_Pol4Exp") != std::string::npos) {
    for (int i = 0; i < 18; i++) {
      fitPars.push_back(histMassFitPars -> GetBinContent(histMassFitPars -> GetXaxis() -> FindBin(Form("%s",parNames_NA60_Pol4Exp[i].data()))));
    }

    // STARTING THE FIT TO THE INVARIANT MASS
    funcMassSigBkg = new TF1("funcMassSigBkg", fitFunctionNA60NEWPol4Exp, minFitRange, maxFitRange, 18);
    for (int i = 1; i < 7; i++) {
      funcMassSigBkg->FixParameter(i, fitPars[i]);
    }
    for (int i = 8; i < 18; i++) {
      funcMassSigBkg->FixParameter(i, fitPars[i]);
    }

    for (int iter = 0;iter < 100;iter++) {
      if (iter == 0) {
        funcMassSigBkg -> SetParameter(0, 1e6); // Bkg normalization
        funcMassSigBkg -> SetParameter(7, 300000); // Jpsi normalization
      } else {
        funcMassSigBkg -> SetParLimits(0, 0, 1e10);
        funcMassSigBkg -> SetParLimits(7, 0, 1e10);
        funcMassSigBkg -> SetParameter(0, funcMassSigBkg -> GetParameter(0)); // Bkg normalization
        funcMassSigBkg -> SetParameter(7, funcMassSigBkg -> GetParameter(7)); // Jpsi normalization
      }
      TFitResultPtr fitResult = histMass -> Fit(funcMassSigBkg, "RL0S");
      int fitStatus = fitResult;  // Fit status
      int covMatrixStatus = fitResult -> CovMatrixStatus(); // cov. matrix status

      if (fitStatus == 0 && covMatrixStatus == 3) {
        break;
      }
    }

    funcMassBkg = new TF1("funcMassBkg", FitFunctionBackgroundPol4Exp, minFitRange, maxFitRange, 7);
    funcMassBkg -> SetLineColor(kGray+1);
    funcMassBkg -> SetLineStyle(kDashed);
    for (int iPar = 0;iPar < 7;iPar++) {
      funcMassBkg -> FixParameter(iPar, funcMassSigBkg -> GetParameter(iPar));
    }

    funcMassSig = new TF1("funcMassSig", FitFunctionNA60New, minFitRange, maxFitRange, 11);
    funcMassSig -> SetLineColor(kAzure+2);
    funcMassSig -> SetLineStyle(kSolid);
    for (int iPar = 0;iPar < 11;iPar++) {
      funcMassSig -> FixParameter(iPar, funcMassSigBkg -> GetParameter(iPar+7));
    }
        
    funcFlowSigBkg  = new TF1("funcFlowSigBkg ", FitFunctionFlowS2NA60NEWPOL4EXPPOL2, minFitRange, maxFitRange, 25);
    funcFlowSigBkg -> SetParNames("pol0","pol1","pol2","pol3","exp1","exp2","exp3","kJPsi","mJPsi","sJPsi","p1LJPsi");
    funcFlowSigBkg -> SetParName(11,"p2LJPsi");
    funcFlowSigBkg -> SetParName(12,"p3LJPsi");
    funcFlowSigBkg -> SetParName(13,"p1RJPsi");   
    funcFlowSigBkg -> SetParName(14,"p2RJPsi");
    funcFlowSigBkg -> SetParName(15,"p3RJPsi");
    funcFlowSigBkg -> SetParName(16,"aLJPsi");
    funcFlowSigBkg -> SetParName(17,"aRJPsi");
    funcFlowSigBkg -> SetParName(18,"kPsiP");
    funcFlowSigBkg -> SetParName(19,"v_{2} JPsi");
    funcFlowSigBkg -> SetParName(20,"v_{2} BG0");
    funcFlowSigBkg -> SetParName(21,"v_{2} BG1");
    funcFlowSigBkg -> SetParName(22,"v_{2} BG2");
    funcFlowSigBkg -> SetParName(23,"v_{2} BG3");
    funcFlowSigBkg -> SetParName(24, "type");
    funcFlowSigBkg -> SetLineColor(kBlue);

    for (int iPar = 0;iPar < 18;iPar++) {
      funcFlowSigBkg -> FixParameter(iPar, funcMassSigBkg -> GetParameter(iPar));
    }
    if (fitFuncs.find("Pol2_v2bkg") != std::string::npos) {
      funcFlowSigBkg -> FixParameter(18, 0);
      funcFlowSigBkg -> SetParameter(19, 0.072);
      funcFlowSigBkg -> SetParameter(20, 0.09);
      funcFlowSigBkg -> SetParameter(21, -0.0252894);
      funcFlowSigBkg -> SetParameter(22, 0.0022179);
      funcFlowSigBkg -> FixParameter(23, 999); // Extra parameter for 3 degree cheby
      funcFlowSigBkg -> FixParameter(24, 0); // v2 bkg is a Pol2
    } else {
      funcFlowSigBkg -> FixParameter(18, 0);
      funcFlowSigBkg -> SetParameter(19, 0.072);
      funcFlowSigBkg -> SetParameter(20, 0.0323879);
      funcFlowSigBkg -> SetParameter(21, -0.0117933);
      funcFlowSigBkg -> SetParameter(22, 0.000888077);
      funcFlowSigBkg -> SetParameter(23, -0.000522725); // Extra parameter for 3 degree cheby
      funcFlowSigBkg -> FixParameter(24, 1); // v2 bkg is a Cheby
    }
    profV2 -> Fit (funcFlowSigBkg, "RI"); // Get the normalization from the fit

    results[0] = funcFlowSigBkg -> GetParameter(19);
    results[1] = funcFlowSigBkg -> GetParError(19);
    results[2] = (double) funcFlowSigBkg -> GetChisquare() / (double) funcFlowSigBkg -> GetNDF();
  }

  TCanvas *canvasFit = new TCanvas("canvasFit", "", 600, 1200);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1, 1);  // Pad superiore
  TPad *pad2 = new TPad("pad2", "pad1", 0, 0, 1, 0.5);  // Pad inferiore

  pad1 -> SetBottomMargin(0);
  pad2 -> SetTopMargin(0);
  pad2 -> SetBottomMargin(0.1);

  canvasFit -> cd();
  pad1 -> Draw();
  pad2 -> Draw();

  pad1 -> cd();
  gPad -> SetLogy(kTRUE);
  histMass -> Draw();
  funcMassBkg -> Draw("SAME");
  funcMassSig -> Draw("SAME");
  funcMassSigBkg -> Draw("SAME");

  pad2 -> cd();
  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(1111);
  profV2 -> Draw();
  funcFlowSigBkg -> Draw("SAME");

  fOut -> cd();
  canvasFit -> Write(Form("v2_fit_%2.1f_%2.1f_%s", minFitRange, maxFitRange, fitFuncs.c_str()));

  if (debug) {
    canvasFit -> SaveAs(Form("%s/v2_fit_%2.1f_%2.1f_%s_Pt_%1.0f_%1.0f.pdf", dirPath.c_str(), minFitRange, maxFitRange, fitFuncs.c_str(), minPtBin, maxPtBin));
  }
  return results;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
std::vector<TH1*> Uncertainities(std::vector<TH1*> histlist, vector<double> parameter, vector<double> parameter_er, double value[nPtBins][3]) {
    std::vector<TH1*> histograms;
    int ncombination = parameter.size() / nPtBins;
    int nbin = 0;
    for (int i = 0; i < parameter.size(); i++) {
        if ((i) % ncombination != 0) {
            continue;
        }
        if ((i) % ncombination == 0) {
            cout<<" HERE IS THE DIVISIONS "<<i<<endl;
        }

        double sum_v2 = 0; 
        double sum_v2er = 0;
        int ibin =0;
        for (int j = i; j < i+ncombination; j++) {
            histlist[nbin]->SetBinContent(ibin+1,parameter[j]);
	        cout<<" v2 parameter inside  "<<i<<" "<< j<<" "<<parameter[j]<<" +/- "<<parameter_er[j]<<endl;
	        //cout<<" v2 parameter inside2  "<<histlist[nbin]->GetBinContent(ibin+1)<<endl;
	        sum_v2 = sum_v2 + parameter[j];
	        sum_v2er = sum_v2er + parameter_er[j];
	        ibin++;
	    }

        histograms.push_back(histlist[nbin]);
        double mean_v2 = sum_v2/ncombination;
        double mean_v2er = sum_v2er/ncombination;
      
        Double_t sum2 =0.0;
        for (int j = i; j < i+ncombination; j++) {
            Double_t dev = (parameter[i] - mean_v2);
	        sum2 = sum2 + dev*dev;	     
	    }

        double_t sys = TMath::Sqrt(sum2/ncombination);
        value[nbin][0] = mean_v2;
        value[nbin][1] = mean_v2er;
        value[nbin][2] = sys;
      
        Double_t per_stat = (mean_v2er/mean_v2)*100;
        Double_t per_sys = (sys/mean_v2)*100;

        cout<<" Average ==== "<<mean_v2<< " +/- " <<mean_v2er<<" ("<<per_stat<<"\%)"<<" +/- "<<sys<<" ("<<per_sys<<"\%)"<<endl;
        //cout<<" Size "<<histlist.size()<<" Bin Content ==== "<<j<<" "<<histlist[0]->GetBinContent(j)<<endl;
        nbin++;
    }
    return histograms;
}


void SetFitRejectRange(Double_t a, Double_t b)
{
  /// Set a range the fit function(s) can ignore
  if ( a <= TMath::Limits<Double_t>::Max() && b <= TMath::Limits<Double_t>::Max() ) TF1::RejectPoint();
}
//------------------------------------------------------------------------------
Double_t FitFunctionBackgroundPol2(Double_t *x, Double_t *par)
{
  int type = par[4];
  // pol2 3 params

  //if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  //if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
  if(type == 0) {
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];  //return pol2 function
  } else {
    //Return chebyshev Pol2 function
    return FitFunctionBackgroundPol2Cheb(x,par); 
  }
}

//------------------------------------------------------------------------------
Double_t FitFunctionBackgroundPol3(Double_t *x, Double_t *par)
{
  // pol2 3 params

  //if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  //if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0];
}

//____________________________________________________________________________
Double_t FitFunctionBackgroundPol4Exp(Double_t *x, Double_t *par)
{
  // pol4 x exp : 6 params
 //if(x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
 //return par[0]*(par[1]+par[2]*x[0]+par[3]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]+par[5]*x[0]*x[0]*x[0]*x[0])*TMath::Exp(par[6]/x[0]);
 //return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0])*TMath::Exp(par[5]/x[0]);
 //return par[0]*(par[1]+par[2]*x[0]+par[3]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]+par[5]*x[0]*x[0]*x[0]*x[0])*TMath::Exp(par[6]*x[0]);
  return par[0]*(par[1]+par[2]*x[0]+par[3]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]+par[5]*x[0]*x[0]*x[0]*x[0])*TMath::Exp(-par[6]*x[0]);
  
}

//____________________________________________________________________________
// Funzione per calcolare il polinomio di Chebyshev di primo tipo di grado n
double ChebyshevT(int n, double x) {
    if (n == 0)
        return 1.0;
    if (n == 1)
        return x;
    return 2 * x * ChebyshevT(n - 1, x) - ChebyshevT(n - 2, x);
}

// Funzione per mappare x da [a, b] a [-1, 1]
double MapToChebyshevRange(double x, double a, double b) {
    return (2 * x - (b + a)) / (b - a);
}

// Wrapper per TF1: combina i termini di Chebyshev con i coefficienti
double FitFunctionBackgroundPol2Cheb(double *x, double *params) {
    double a = 2.3;  // Limite inferiore dell'intervallo
    double b = 4.7;  // Limite superiore dell'intervallo

    double xPrime = MapToChebyshevRange(x[0], a, b);

    // Somma i termini di Chebyshev con i rispettivi coefficienti
    double result = params[0] * ChebyshevT(0, xPrime) +
                    params[1] * ChebyshevT(1, xPrime) +
                    params[2] * ChebyshevT(2, xPrime) +
                    params[3] * ChebyshevT(3, xPrime);

    return result;
}
//____________________________________________________________________________

Double_t alphaCB2VWG(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[4])/(FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionBackgroundVWG(x,par));
}

//____________________________________________________________________________
Double_t FitFunctionSignalCrystalBallExtended(Double_t *x,Double_t *par)
{
  // Extended Crystal Ball : 7 parameters
  //
  // par[0] = Normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = alpha
  // par[4] = n
  // par[5] = alpha'
  // par[6] = n'

  Double_t t = (x[0]-par[1])/par[2];
  if (par[3] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[3]);
  Double_t absAlpha2 = fabs((Double_t)par[5]);

  if (t >= -absAlpha && t < absAlpha2) // gaussian core
  {
    return par[0]*(exp(-0.5*t*t));
  }

  if (t < -absAlpha) //left tail
  {
    Double_t a =  TMath::Power(par[4]/absAlpha,par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b = par[4]/absAlpha - absAlpha;
    return par[0]*(a/TMath::Power(b - t, par[4]));
  }

  if (t >= absAlpha2) //right tail
  {

    Double_t c =  TMath::Power(par[6]/absAlpha2,par[6])*exp(-0.5*absAlpha2*absAlpha2);
    Double_t d = par[6]/absAlpha2 - absAlpha2;
    return par[0]*(c/TMath::Power(d + t, par[6]));
  }

  return 0. ;
}

//-------------------------------------------------------------------------------
Double_t FitFunctionBackgroundVWG(Double_t *x, Double_t *par)
{
  // gaussian variable width : 4 params

  //if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  //if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
  Double_t sigma = par[2]+par[3]*((x[0]-par[1])/par[1]);
  return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.*sigma*sigma));
}


//-------------------------------------------------------------------------------
Double_t FitFunctionFlowS2CB2VWGPOL2(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = 1.0341;

  Double_t par2[11] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[11], //kPsi'
    par[5]+(3.68609-3.096916),
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };

  return alphaCB2VWG(x,par)*par[12] + (1. - alphaCB2VWG(x,par))*FitFunctionBackgroundPol2(x,&par[13]);
}



//-------------------------------------------------------------------------------
Double_t FitFunctionFlowS2CB2VWGPol3Cheb(Double_t *x, Double_t *par)
{
  // Fit function for Jpsi(Psip) mean pt with alphaJpsi and alphaPsiP
  Double_t SPsiPFactor = 1.0341;

  Double_t par2[11] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[11], //kPsi'
    par[5]+(3.68609-3.096916),
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10]
  };


  return alphaCB2VWG(x,par)*par[12] + (1. - alphaCB2VWG(x,par))*FitFunctionBackgroundPol2Cheb(x,&par[13]);
}


//-------------------------------------------------------------------------------
Double_t FitFunctionFlowS2CB2POL2EXPPOL2(Double_t *x, Double_t *par)
{
  Double_t SPsiPFactor = 1.0341;
  return alphaCB2POL2EXP(x,par)*par[12]  + (1. - alphaCB2POL2EXP(x,par))*FitFunctionBackgroundPol2(x,&par[13]);

}

//-------------------------------------------------------------------------------
Double_t FitFunctionFlowS2CB2POL4EXPPOL2(Double_t *x, Double_t *par)
{
  return alphaCB2POL4EXP(x,par)*par[15]  + (1. - alphaCB2POL4EXP(x,par))*FitFunctionBackgroundPol2(x,&par[16]);

}

//-------------------------------------------------------------------------------
Double_t alphaCB2POL2EXP(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[4])/(FitFunctionSignalCrystalBallExtended(x, &par[4]) + FitFunctionBackgroundPol2Exp(x,par));
}

//-------------------------------------------------------------------------------
Double_t alphaCB2POL4EXP(Double_t*x, Double_t* par)
{
  return FitFunctionSignalCrystalBallExtended(x, &par[7])/(FitFunctionSignalCrystalBallExtended(x, &par[7]) + FitFunctionBackgroundPol4Exp(x,par));
}

//____________________________________________________________________________
Double_t FitFunctionBackgroundPol2Exp(Double_t *x, Double_t *par)
{
  // pol2 x exp : 4 params
  //if (x[0] > 2.9 && x[0] < 3.3) TF1::RejectPoint();
  //if (x[0] > 3.5 && x[0] < 3.7) TF1::RejectPoint();
//  return par[0]*(par[1]+par[2]*x[0]+par[3]*x[0]*x[0])*TMath::Exp(par[4]/x[0]);
  return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0])*TMath::Exp(par[3]*x[0]);
}



//____________________________________________________________________________
Double_t FitFunctionNA60New(Double_t *x,Double_t *par)
{
  // New Formulation of NA60 : 11 parameters
  //
  // par[0] = Normalization
  // par[1] = mean
  // par[2] = sigma
  // par[3] = p1Left
  // par[4] = p2Left
  // par[5] = p3Left
  // par[6] = p1Right
  // par[7] = p2Right
  // par[8] = p3Right
  // par[9] = alphaLeft
  // par[10] = alphaRight


  const Double_t t = (x[0]-par[1])/par[2];

  Double_t sigmaRatio(0.);
  if( t < par[9] ) sigmaRatio = ( 1.0 + TMath::Power( par[3]*(par[9]-t), par[4]-par[5]*TMath::Sqrt(par[9] - t) ) );
  else if( t >= par[9] && t < par[10] ) sigmaRatio = 1;
  else if( t >= par[10] ) sigmaRatio = ( 1.0 + TMath::Power( par[6]*(t-par[10]), par[7]-par[8]*TMath::Sqrt(t - par[10]) ) );

  return par[0]*TMath::Exp( -(1/2.)*TMath::Power(t/sigmaRatio,2.));

}

//Mass Fits
//------------------------------------------------------------------------------
Double_t fitFunctionCB2VWG(Double_t *x, Double_t *par)
{
  return FitFunctionBackgroundVWG(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[4]);
}

//------------------------------------------------------------------------------
Double_t fitFunctionNA60NEWVWG(Double_t *x, Double_t *par)
{
  return FitFunctionBackgroundVWG(x, par) + FitFunctionNA60New(x, &par[4]);
}

//------------------------------------------------------------------------------
Double_t fitFunctionCB2Pol4Exp(Double_t *x, Double_t *par)
{
  return FitFunctionBackgroundPol4Exp(x, par) + FitFunctionSignalCrystalBallExtended(x, &par[7]);
}

//------------------------------------------------------------------------------
Double_t fitFunctionNA60NEWPol4Exp(Double_t *x, Double_t *par)
{
  return FitFunctionBackgroundPol4Exp(x, par) + FitFunctionNA60New(x, &par[7]);
}

//------------------------------------------------------------------------------
Double_t alphaNA60NEWVWG(Double_t*x, Double_t* par)
{
  return FitFunctionNA60New(x, &par[4])/(FitFunctionNA60New(x, &par[4]) + FitFunctionBackgroundVWG(x,par));
}

//------------------------------------------------------------------------------
Double_t alphaNA60NEWPOL2EXP(Double_t*x, Double_t* par)
{
  return FitFunctionNA60New(x, &par[4])/(FitFunctionNA60New(x, &par[4]) + FitFunctionBackgroundPol2Exp(x,par));
}

//------------------------------------------------------------------------------
Double_t alphaNA60NEWPOL4EXP(Double_t*x, Double_t* par)
{
  return FitFunctionNA60New(x, &par[7])/(FitFunctionNA60New(x, &par[7]) + FitFunctionBackgroundPol4Exp(x,par));
}

//____________________________________________________________________________
Double_t FitFunctionFlowS2NA60NEWPOL2EXPPOL2(Double_t *x, Double_t *par)
{
  Double_t SPsiPFactor = 1.0341;
  Double_t par2[15] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[15],
    par[5]+(3.68609-3.096916),
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
  };

  return alphaNA60NEWPOL2EXP(x,par)*par[16] + (1. - alphaNA60NEWPOL2EXP(x,par))*FitFunctionBackgroundPol2(x,&par[17]);

}

//____________________________________________________________________________
Double_t FitFunctionFlowS2NA60NEWPOL4EXPPOL2(Double_t *x, Double_t *par)
{
  return alphaNA60NEWPOL4EXP(x,par)*par[19] + (1. - alphaNA60NEWPOL4EXP(x,par))*FitFunctionBackgroundPol2(x,&par[20]);
}


//____________________________________________________________________________
Double_t FitFunctionFlowS2NA60NEWVWGPOL2(Double_t *x, Double_t *par)
{
  Double_t SPsiPFactor = 1.0341;

  Double_t par2[15] = {
    par[0],
    par[1],
    par[2],
    par[3],
    par[15],
    par[5]+(3.68609-3.096916),
    par[6]*SPsiPFactor, // /3.096916*3.68609,
    par[7],
    par[8],
    par[9],
    par[10],
    par[11],
    par[12],
    par[13],
    par[14],
  };

  return alphaNA60NEWVWG(x,par)*par[16] + (1. - alphaNA60NEWVWG(x,par))*FitFunctionBackgroundPol2(x,&par[17]);

}

