#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TMCParticle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TSystem.h>
#include <TPythia6Decayer.h>
#endif

double FuncPt(double *, double *);
double FuncEta(double *, double *);
double FuncAxe(double *, double *);

void JpsiToyMC(int nEvts = 100000) {
  std::cout << "To Run this code you need AliGenerators locally or LXPLUS" << std::endl;
  std::cout << "on LXPLUS do /cvmfs/alice.cern.ch/bin/alienv enter VO_ALICE@AliGenerators::v20250107-1" << std::endl;

  gSystem -> Load("liblhapdf.so");      // Parton density functions
  gSystem -> Load("libEGPythia6.so");   // TGenerator interface
  gSystem -> Load("libpythia6.so");     // Pythia
  gSystem -> Load("libAliPythia6.so");  // ALICE specific implementations

  TDatabasePDG *database = TDatabasePDG::Instance();
  int pdgCode = 443;
  double massJpsi = database -> GetParticle(pdgCode) -> Mass();

  TPythia6Decayer *pythiaDec=TPythia6Decayer::Instance();
  pythiaDec -> ForceParticleDecay(443,13,2);
  pythiaDec -> Init();
  
  TClonesArray *array = new TClonesArray("TParticle",100);
  TRandom3 *gener = new TRandom3(0);
  TLorentzVector *vecJpsi = new TLorentzVector();
  TLorentzVector *vecMup = new TLorentzVector();
  TLorentzVector *vecMum = new TLorentzVector();

  TH1F *histGenPtJpsi = new TH1F("histGenPtJpsi", " ; #it{p}_{T} (GeV/c)", 20, 0, 10);
  TH1F *histRecPtJpsi = new TH1F("histRecPtJpsi", " ; #it{p}_{T} (GeV/c)", 20, 0, 10);
  TH1F *histRecMidEffPtJpsi = new TH1F("histRecMidEffPtJpsi", " ; #it{p}_{T} (GeV/c)", 20, 0, 10);
  TH2F *histGenPtJpsiPtMup = new TH2F("histGenPtJpsiPtMup", " #it{p}_{T}^{J/#psi} (GeV/c) ; #it{p}_{T}^{#mu^{+}} (GeV/c)", 20, 0, 10, 20, 0, 10);
  TH2F *histGenPtJpsiPtMum = new TH2F("histGenPtJpsiPtMum", " #it{p}_{T}^{J/#psi} (GeV/c) ; #it{p}_{T}^{#mu^{-}} (GeV/c)", 20, 0, 10, 20, 0, 10);

  TF1 *funcPt = new TF1("funcPt", FuncPt, 0., 15., 3);
  funcPt -> SetParameter(0, 2.85);
  funcPt -> SetParameter(1, 2.81);
  funcPt -> SetParameter(2, 2.43);

  TF1 *funcEta = new TF1("funcEta", FuncEta, -4., -2.5, 2);
  funcEta -> SetParameter(0, -2.43);
  funcEta -> SetParameter(1, 1.08);

  TF1 *funcPhi = new TF1("funcPhi", "[0]", 0., 2. * TMath::Pi());
  funcPhi -> SetParameter(0, 1);

  TF1 *funcAxe = new TF1("funcAxe", FuncAxe, 0., 15., 5);
  funcAxe -> SetParameter(0, 0.416885);
  funcAxe -> SetParameter(1, -0.0609379);
  funcAxe -> SetParameter(2, 0.0222385);
  funcAxe -> SetParameter(3, -0.00171191);
  funcAxe -> SetParameter(4, 4.08582e-05);

  bool passMidEffMup = false;
  bool passMidEffMum = false;
  double axeCut = 0;
  for (int iEvt = 0;iEvt < nEvts;iEvt++) {
    passMidEffMup = false;
    passMidEffMum = false;
    axeCut = 0;

    double ptJpsi = funcPt -> GetRandom();
    double etaJpsi = funcEta -> GetRandom();
    double phiJpsi = funcPhi -> GetRandom();

    vecJpsi -> SetPtEtaPhiM(ptJpsi, etaJpsi, phiJpsi, massJpsi);
    if (vecJpsi -> Pt() > 15) continue;
    pythiaDec -> Decay(pdgCode, vecJpsi);
    array -> Clear();
    int nEntries = pythiaDec -> ImportParticles(array);
    TParticle *daughters = (TParticle*) array -> At(0);
    int nDaughters = daughters -> GetNDaughters();
    for(int j = 0; j < nEntries;j++) {
      TParticle *mcPart = (TParticle*) array -> At(j);
      int pdgDaughter = mcPart -> GetPdgCode();
      if (pdgDaughter == 13) {
        double pxMup = mcPart -> Px();
        double pyMup = mcPart -> Py();
        double pzMup = mcPart -> Pz();
        double enMup = mcPart -> Energy();
        vecMup -> SetPxPyPzE(pxMup, pyMup, pzMup, enMup);
        histGenPtJpsiPtMup -> Fill(vecMup -> Pt(), vecJpsi -> Pt());
        if (gRandom -> Rndm() < 0.91) {
          passMidEffMup = true;
        }
      }
      if (pdgDaughter == -13) {
        double pxMum = mcPart -> Px();
        double pyMum = mcPart -> Py();
        double pzMum = mcPart -> Pz();
        double enMum = mcPart -> Energy();
        vecMum -> SetPxPyPzE(pxMum, pyMum, pzMum, enMum);
        histGenPtJpsiPtMum -> Fill(vecMum -> Pt(), vecJpsi -> Pt());
        if (gRandom -> Rndm() < 0.91) {
          passMidEffMum = true;
        }
      }
    }

    histGenPtJpsi -> Fill(vecJpsi -> Pt());
    axeCut = funcAxe -> Eval(vecJpsi -> Pt());
    if (gRandom -> Rndm() < axeCut) {
      histRecPtJpsi -> Fill(vecJpsi -> Pt());
      if (passMidEffMup && passMidEffMum) {
        histRecMidEffPtJpsi -> Fill(vecJpsi -> Pt());
      }
    }
  }

  TH1F *histAxePtJpsi = (TH1F*) histRecPtJpsi -> Clone("histAxePtJpsi");
  histAxePtJpsi -> Divide(histGenPtJpsi);

  TH1F *histAxeMidEffPtJpsi = (TH1F*) histRecMidEffPtJpsi -> Clone("histAxeMidEffPtJpsi");
  histAxeMidEffPtJpsi -> Divide(histGenPtJpsi);

  TFile *fOut = new TFile("jpsi_toy_mc.root", "RECREATE");
  histGenPtJpsi -> Write();
  histRecPtJpsi -> Write();
  histAxePtJpsi -> Write();
  histAxeMidEffPtJpsi -> Write();
  histGenPtJpsiPtMup -> Write();
  histGenPtJpsiPtMum -> Write();
  fOut -> Close();
}
//////////////////////////////
double FuncPt(double *x, double *par) {
  double arg = 1 + TMath::Power(x[0]/par[0], par[1]);
  return x[0] / TMath::Power(arg, par[2]);
}
//////////////////////////////
double FuncEta(double *x, double *par) {
  double arg = 0.5 * ((x[0] - par[0]) / par[1]) * ((x[0] - par[0]) / par[1]);
  return TMath::Exp(arg);
}
//////////////////////////////
double FuncAxe(double *x, double *par) {
  double xx = x[0];
  return par[0] + par[1] * xx + par[2] * xx * xx + par[3] * xx * xx * xx + par[4] * xx * xx * xx * xx;
}