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

#include "TPythia8.h"

#endif

void decayOnia() {

    TFile *fIn = TFile::Open("pythia8_onia_kMonash_kSoftQCD.root");
    TTree *tree = (TTree*) fIn->Get("tuplePairs");

    float pTOnia; 
    float yOnia; 
    float etaOnia; 
    float phiOnia; 

    tree->SetBranchAddress("pTOnia", &pTOnia); 
    tree->SetBranchAddress("yOnia", &yOnia); 
    tree->SetBranchAddress("etaOnia", &etaOnia); 
    tree->SetBranchAddress("phiOnia", &phiOnia);

    TH1F *hPtJpsi = new TH1F("hPtJpsi", ";#it{p}_{T}^{J/#psi}", 100, 0, 10);
    TH1F *hPtRap = new TH1F("hPtRap", ";#it{y}^{J/#psi}", 100, -5, 5);
    TH2F *hPtJpsiPtMu = new TH2F("hPtJpsiPtMu", ";#it{p}_{T}^{J/#psi};#it{p}_{T}^{#mu^{+}}", 100, 0, 10, 100, 0, 10);

    
    Pythia8::Pythia pythia;

    pythia.readString("ProcessLevel:all = off");

    pythia.readString("443:onMode = off");
    pythia.readString("443:onIfMatch = 13 -13");

    pythia.init();

    double pxMuPlus; double pyMuPlus; double pzMuPlus; double eMuPlus; double pxMuMinus; double pyMuMinus; double pzMuMinus; double eMuMinus;

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) { 
        tree->GetEntry(i);
        TLorentzVector *jpsi = new TLorentzVector();
        jpsi->SetPtEtaPhiM(pTOnia, etaOnia, phiOnia, 3.096);

        tree->GetEntry(i);
        pythia.event.reset();
        
        Pythia8::Particle jpsiParticle(
        443,              // id
        1,                // status
        0,                // mother1
        0,                // mother2
        0,                // daughter1
        0,                // daughter2
        0,                // col
        0,                // acol
        jpsi->Px(),
        jpsi->Py(),
        jpsi->Pz(),
        jpsi->E(),
        jpsi->M()
        );

        pythia.event.append(jpsiParticle);

        pythia.moreDecays();

        bool foundMuPlus = false; 
        bool foundMuMinus = false; 
        for (int j = 0; j < pythia.event.size(); ++j) { 
            const auto &part = pythia.event[j]; 
            if (part.id() == -13) { 
                pxMuPlus = part.px(); 
                pyMuPlus = part.py(); 
                pzMuPlus = part.pz(); 
                eMuPlus = part.e(); 
                foundMuPlus = true; 
            }

            if (part.id() == 13) { 
                pxMuMinus = part.px(); 
                pyMuMinus = part.py(); 
                pzMuMinus = part.pz(); 
                eMuMinus = part.e(); 
                foundMuMinus = true; 
            }

            if (foundMuPlus && foundMuMinus) {
                float ptJpsi = TMath::Sqrt(jpsi->Px()*jpsi->Px() + jpsi->Py()*jpsi->Py());
                float ptMu = TMath::Sqrt(pxMuPlus*pxMuPlus + pyMuPlus*pyMuPlus);
                //if (jpsi->Y() > 2.5 && jpsi->Y() < 4) {
                    hPtJpsi->Fill(ptJpsi);
                    hPtRap->Fill(jpsi->Y());
                    hPtJpsiPtMu->Fill(ptJpsi, ptMu);
                //}
            }
        }
    }

    TCanvas *canvas = new TCanvas("canvas", "", 1200, 1200);
    canvas->Divide(2,2);
    canvas->cd(1);
    hPtJpsiPtMu->Draw("COLZ");
    canvas->cd(2);
    hPtJpsi->Draw("H");
    canvas->cd(3);
    hPtRap->Draw("H");

}