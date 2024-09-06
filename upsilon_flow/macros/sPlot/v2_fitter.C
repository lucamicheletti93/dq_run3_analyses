
#include "TRandom.h"
#include "TF1.h"
#include "TTree.h"
#include "TNtuple.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "TChain.h"
#include "RooCBShape.h"
#include "TFile.h"

// use this order for safety on library loading
using namespace RooFit ;
using namespace RooStats ;

// see below for implementation
void AddModel(RooWorkspace*);
void AddData(RooWorkspace*, string, string);
void DoSPlot(RooWorkspace*);
void MakePlotsAndTree(RooWorkspace*, string, string);

void create_toy_sample(float lumi = 1.0);

//____________________________________
void v2_fitter(bool useToy = false) {
  string fInName, treeName;

  if (useToy) {
    cout << "----------- USE TOY DATA -----------" << endl;
    fInName = "/Users/lucamicheletti/GITHUB/dq_run3_analyses/upsilon_flow/macros/sPlot/toy_MC.root";
    treeName = "toyv2";
  } else {
    cout << "----------- USE REAL DATA -----------" << endl;
    fInName = "/Users/lucamicheletti/GITHUB/dq_run3_analyses/upsilon_flow/data/train_261116/AO2D_centrality_cut.root";
    treeName = "DF_2338657052269120/O2rtdileptmtree";
  }

  //Create a new workspace to manage the project
  RooWorkspace* wspace = new RooWorkspace("myWS");

  //add the signal and background modesl to the workspace.
  AddModel(wspace);
  
  //add toy data
  AddData(wspace, fInName, treeName);

  //inspect the workspace
  wspace->Print();

  //do acutal splot produces new dataset with sWeights for every candidate considered
  DoSPlot(wspace);
  
  //Make control plots
  MakePlotsAndTree(wspace, fInName, treeName);

  //cleanup
  delete wspace;
}
///////////////////////////////////////////////////
void create_toy_sample(float lumi/*/nb*/) {
  double minMassBin = 8.5;
  double maxMassBin = 12.0;
  
  TF1 *funcMassUps = new TF1("funcMassUps", CrystalBall, minMassBin, maxMassBin, 5);
  funcMassUps -> SetParameters(200, 9.46, 0.13, 0.82, 2.44);

  TF1 *funcMassBkg = new TF1("funcMassBkg", "[0]*exp([1]*x)", minMassBin, maxMassBin);
  funcMassBkg -> SetParameter(0,10);
  funcMassBkg -> SetParameter(1,-0.2);

  //v2 as function of mass, just powerlaw from Olitraut
  TF1 *funcFlowUps = new TF1("funcFlowUps", v2eventbyevent, 0.0, 1.0, 1);
  funcFlowUps -> SetParameter(0, 5.64); //from https://arxiv.org/pdf/1312.6555.pdf, equation 4
  
  //v2 as function of mass, take SAME powerlaw, but with different index
  TF1 *funcFlowBkg = new TF1("funcFlowBkg", v2eventbyevent, 0.0, 1.0, 1);//
  funcFlowBkg -> SetParameter(0, 5.64); //only event-by-event shape

  TF1* funcSigDeltaPhi = new TF1("funcSigDeltaPhi","1 + 2*[0]*cos(2*x)");
  funcSigDeltaPhi -> SetParameter(0, 0.05);
  funcSigDeltaPhi -> SetRange(0, TMath::Pi());
  
  TF1* funcBkgDeltaPhi = new TF1("funcBkgDeltaPhi", "1 + 2*[0]*cos(2*x)");
  funcBkgDeltaPhi -> SetParameter(0, 0.025);
  funcBkgDeltaPhi -> SetRange(0, TMath::Pi());

  TNtuple *tuple = new TNtuple("toyv2", "toyv2", "fMass:fCos2DeltaPhi:label");

  // Toy MC input variables
  double sigToBkg = 0.1;
  double sigToSigBkg = 1. / (1 + 1.0/ (sigToBkg) );
  double v2Fraction = 0.8;
  int nGenEvs = int(lumi * 1000) / sigToBkg;//rough number of candidates

  // Toy MC settings
  double seed = 0.0;
  double fMass = 0.0;
  double deltaPhi = 0.0;
  double fCos2DeltaPhi = 0.0;
  double label = 0.0;
  
  for (int iGenEv = 0;iGenEv < nGenEvs;iGenEv++) {
    seed = gRandom -> Rndm();
    if (seed < sigToSigBkg) {
      label = 1.0;
      fMass = funcMassUps -> GetRandom();
      deltaPhi = funcSigDeltaPhi -> GetRandom();
      fCos2DeltaPhi = TMath::Cos(2 * deltaPhi);
    } else {
      label = 0.0;
      fMass = funcMassBkg -> GetRandom();
      deltaPhi = funcBkgDeltaPhi -> GetRandom();
      fCos2DeltaPhi = TMath::Cos(2 * deltaPhi);
    }

    double data[] = {fMass, fCos2DeltaPhi, label};
    
    // Convert to float format
    int nEvents = sizeof(data) / sizeof(data[0]);
    float *dataFloat = new float[nEvents];
    for(int iEvent = 0;iEvent < nEvents;iEvent++) {
      dataFloat[iEvent] = (float) data[iEvent] ;
    }
    
    tuple -> Fill( dataFloat );
    delete[] dataFloat;
  }

  TFile* fOut = new TFile("toy_MC.root", "RECREATE");
  fOut -> cd();
  tuple -> SetDirectory(fOut);
  tuple -> Write();
  fOut -> Close();
}
///////////////////////////////////////////////////
void AddModel(RooWorkspace* ws) {
  double minMassBin = 8.5;
  double maxMassBin = 12.0;

  //make signal model
  RooRealVar fMass("fMass", "m_{#mu#mu}", 10.0, minMassBin, maxMassBin, "GeV");
  RooRealVar meanUps("meanUps", "Upsilon Mass", 9.46, 9.1, 9.9, "GeV");
  RooRealVar widthUps("widthUps", "Upsilon resolution", 0.13, 0.0, 1.0 , "GeV");
  RooRealVar alphaUps("alphaUps","CB alphaUps", 0.82, 0.0, 2.0, "");
  RooRealVar nUps("nUps","CB n", 2.44, 1.0, 4.0, "");
  
  // Fix variables in the fit
  alphaUps.setConstant();
  nUps.setConstant();
  
  RooCBShape modelUps("CBmodel", "CBmodel", fMass, meanUps, widthUps, alphaUps, nUps);
 
  //background model
  RooRealVar decConst("decConst", "decay const", -0.2, -1.0, 0.0);
  //todo implement actual function used
  RooExponential modelBkg("modelBkg", "modelBkg", fMass, decConst);
  //decConst.setConstant();

  //yields
  RooRealVar yieldUps("yieldUps", "fitted Upsilon yield", 6000, 10, 100000, "");
  RooRealVar yieldBkg("yieldBkg", "fitted bkg yield", 10000, 10, 1000000, "");
  
  //combined model
  RooAddPdf model("model", "Upsilon+bkg", RooArgList(modelUps, modelBkg), RooArgList(yieldUps, yieldBkg));
  ws -> import(model);
}
///////////////////////////////////////////////////
void AddData(RooWorkspace *ws, string fInName, string treeName){
  double minMassBin = 8.5;
  double maxMassBin = 12.0;

  RooRealVar fMass("fMass", "m_{#mu#mu}", 10, minMassBin, maxMassBin, "GeV");
  RooRealVar fCos2DeltaPhi("fCos2DeltaPhi", "cos(2#Delta#phi)", 0, -1.0, 1.0, "");

  // Get data from file
  TFile *fIn = new TFile(fInName.c_str(), "READ");
  TTree *realdata = (TTree*) fIn -> Get(treeName.c_str());
  RooDataSet data("data", "data", realdata, RooArgSet(fMass, fCos2DeltaPhi));

  //import data into workspace
  ws->import(data, Rename("data"));
}
///////////////////////////////////////////////////
void DoSPlot(RooWorkspace *ws){
  // Get all ingredients
  RooAbsPdf* model = ws -> pdf("model");
  RooRealVar* yieldUps = ws -> var("yieldUps");
  RooRealVar* yieldBkg = ws -> var("yieldBkg");
  RooDataSet* data = (RooDataSet*) ws -> data("data");

  // Fit the model to the data.
  model -> fitTo(*data, Extended(), Minos(true), Strategy(2));

  // The sPlot technique requires that we fix the parameters of the model that are not yields after doing the fit
  RooRealVar* widthUps = ws -> var("widthUps");
  RooRealVar* meanUps = ws -> var("meanUps");
  RooRealVar* alphaUps = ws -> var("alphaUps");
  RooRealVar* nUps = ws -> var("nUps");
  RooRealVar* decConst = ws -> var("decConst");
  widthUps -> setConstant();
  meanUps -> setConstant();
  alphaUps -> setConstant();
  nUps -> setConstant();
  decConst -> setConstant();

  RooMsgService::instance().setSilentMode(true);

  //Now we use the SPlot class to add SWeight to our data set based on our model and our yield variables
  RooStats::SPlot * sData = new RooStats::SPlot("sData","splot", *data, model, RooArgList(*yieldUps, *yieldBkg) );

  //Check Sweight properties
  std::cout << std::endl << "Yield of Upsilon is " << yieldUps -> getVal() << ". From sWeights it is " << sData -> GetYieldFromSWeight("yieldUps") << std::endl;
  std::cout << std::endl << "Yield of bkg is " << yieldBkg -> getVal() << ". From sWeights it is " << sData -> GetYieldFromSWeight("yieldBkg") << std::endl;

  for (int i = 0; i < 10; i ++) {
    std::cout << "Ups Weight "<< sData -> GetSWeight(i, "yieldUps") << " bkg Yield "<< sData -> GetSWeight(i, "yieldBkg") << "Total Weight "<< sData -> GetSumOfEventSWeight(i) << std::endl;
  }

   //import the new data set with Sweight
   ws -> import(*data, Rename("dataWithSWeights"));
}
///////////////////////////////////////////////////
void MakePlotsAndTree(RooWorkspace *ws, string fInName, string treeName){
  // Open input trees
  //TFile *fIn = new TFile(fInName.c_str(), "READ");
  //TTree *realdata = (TTree*) fIn -> Get(treeName.c_str());
  
  // Get PDFs fron WS
  RooAbsPdf *model = ws -> pdf("model");
  RooAbsPdf *modelUps = ws -> pdf("CBmodel");
  RooAbsPdf *modelBkg = ws -> pdf("modelBkg");

  // Get vars from WS
  RooRealVar *fMass = ws->var("fMass");
  RooRealVar *fCos2DeltaPhi = ws->var("fCos2DeltaPhi");
  
  // Get data from WS
  RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");

  TLatex *latexTitle = new TLatex();
  latexTitle -> SetTextSize(0.07);
  latexTitle -> SetNDC();
  latexTitle -> SetTextFont(42);

  TCanvas* canvasMassFit = new TCanvas("canvasMassFit", "canvasMassFit", 1200, 400);
  canvasMassFit->Divide(3);

  canvasMassFit->cd(1);
  RooPlot* frame = fMass -> frame();
  data -> plotOn(frame);
  model -> plotOn(frame);
  model -> plotOn(frame, Components(*modelUps), LineStyle(kDashed), LineColor(kRed));
  model -> plotOn(frame, Components(*modelUps), LineStyle(kDashed), LineColor(kRed));
  model -> plotOn(frame, Components(*modelBkg), LineStyle(kDashed), LineColor(kGreen));
  frame -> SetTitle("Fit to model to discriminating variable");
  frame -> Draw();

  canvasMassFit -> cd(2);
  data -> plotOn(frame);
  model -> plotOn(frame);
  RooHist *hpull = frame -> pullHist();
  hpull -> Draw();

  canvasMassFit -> cd(3);
  RooHist *hresid = frame -> residHist();
  hresid -> Draw();
  
  //now other variable
  TCanvas *canvasSplot = new TCanvas("Splot_control", "Splot_control", 1200, 400);
  canvasSplot -> Divide(2);

  canvasSplot -> cd(1);

  RooDataSet *dataSwUps = new RooDataSet(data -> GetName(), data -> GetTitle(), data, *data -> get(), 0, "yieldUps_sw");
  cout << "dataSwUps->sumEntries() "  << dataSwUps -> sumEntries() << endl;
  RooPlot *frameFlowUps = fCos2DeltaPhi -> frame();
  dataSwUps -> plotOn(frameFlowUps, DataError(RooAbsData::SumW2));

  TH1F *histFlowUps = (TH1F*) dataSwUps -> createHistogram("fCos2DeltaPhi");
  std::cout << "Upsilon v2 from sWeight " << histFlowUps -> GetMean() << " +/- " << histFlowUps -> GetMeanError() << std::endl;

  frameFlowUps -> SetTitle("cos(2#Delta#phi) #Upsilon(1S) distribution");
  frameFlowUps -> Draw();
  //realdata -> Draw("fCos2DeltaPhi", "SAME");

  latexTitle -> DrawLatex(0.30, 0.65, Form("v_{2}^{#varUpsilon(1S)} = %5.4f #pm %5.4f", histFlowUps -> GetMean(), histFlowUps -> GetMeanError()));

  canvasSplot -> cd(2);
  RooDataSet *dataSwBkg = new RooDataSet(data -> GetName(), data->GetTitle(), data, *data -> get(), 0, "yieldBkg_sw");
  RooPlot *frameFlowBkg = fCos2DeltaPhi -> frame();
  dataSwBkg -> plotOn(frameFlowBkg, DataError(RooAbsData::SumW2));
  TH1F *histFlowBkg = (TH1F*) dataSwBkg -> createHistogram("fCos2DeltaPhi");
  std::cout << "Background v2 from sWeight " << histFlowBkg -> GetMean() << " +/- " << histFlowBkg -> GetMeanError() << std::endl;

  frameFlowBkg -> SetTitle("cos(2#Delta#phi) background distribution");
  frameFlowBkg -> Draw();

  //realdata->Draw("fCos2DeltaPhi","label<0.1","SAME");
  //realdata -> Draw("fCos2DeltaPhi", "same");

  latexTitle -> DrawLatex(0.30, 0.65, Form("v_{2}^{Bkg} = %5.4f #pm %5.4f", histFlowBkg -> GetMean(), histFlowBkg -> GetMeanError()));

  //write resulting tree
  TFile *fOut = new TFile("dataWithSweights.root", "RECREATE");
  fOut -> cd();

  //data->setDefaultStorageType(RooAbsData::Tree);
  (data -> GetClonedTree()) -> Write();
  fOut -> Close();
}



