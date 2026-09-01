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

enum Particle {
    kJpsi,
    kPsi2s,
    kUpsilon1s,
    kNparticles
};

const int particlePdg[kNparticles] = {
    443,
    100443,
    553
};

const string particleName[kNparticles] = {
    "jpsi",
    "psi2s",
    "upsilon1s"
};

const string particleTitle[kNparticles] = {
    "J/#psi",
    "#psi(2S)",
    "#Upsilon(1S)"
};

// Opening angle: cosTheta = (E^2 - 2M^2 + 4m^2) / (E^2 - 4m^2)
void decayOnia(Particle particle = kJpsi) {
    TDatabasePDG *database = TDatabasePDG::Instance();
    double momMass = database->GetParticle(particlePdg[particle])->Mass();

    TFile *fIn = TFile::Open("pythia8_onia_kMonash_kSoftQCD.root");
    TTree *tree = (TTree*) fIn->Get("tuplePairs");

    float pTOnia; 
    float yOnia; 
    float etaOnia; 
    float phiOnia; 
    float absPdg;

    tree->SetBranchAddress("pTOnia", &pTOnia); 
    tree->SetBranchAddress("yOnia", &yOnia); 
    tree->SetBranchAddress("etaOnia", &etaOnia); 
    tree->SetBranchAddress("phiOnia", &phiOnia);
    tree->SetBranchAddress("absPdg", &absPdg);

    TH1F *hPt = new TH1F("hPt", Form(";#it{p}_{T}^{%s}", particleTitle[particle].c_str()), 100, 0, 10); hPt->SetLineColor(kBlack);
    TH1F *hRap = new TH1F("hRap", Form(";#it{y}^{%s}", particleTitle[particle].c_str()), 100, -5, 5); hRap->SetLineColor(kBlack);
    TH1F *hPtFwd = new TH1F("hPtFwd", Form(";#it{p}_{T}^{%s}", particleTitle[particle].c_str()), 100, 0, 10);
    TH1F *hRapFwd = new TH1F("hRapFwd", Form(";#it{y}^{%s}", particleTitle[particle].c_str()), 100, -5, 5);
    TH1F *hPtMu = new TH1F("hPtMu", ";#it{p}_{T}^{#mu}", 100, 0, 10); hPtMu->SetLineColor(kBlack);
    TH1F *hRapMu = new TH1F("hRapMu", ";#it{y}^{#mu}", 100, -5, 5); hRapMu->SetLineColor(kBlack);
    TH1F *hPtMuFwd = new TH1F("hPtMuFwd", ";#it{p}_{T}^{#mu}", 100, 0, 10);
    TH1F *hRapMuFwd = new TH1F("hRapMuFwd", ";#it{y}^{#mu}", 100, -5, 5);
    TH1F *hOpAngle = new TH1F("hOpAngle", ";cos#it{#theta}", 100, -1, 1); hOpAngle->SetLineColor(kBlack);
    TH1F *hOpAngleFwd = new TH1F("hOpAngleFwd", ";cos#it{#theta}", 100, -1, 1);
    TH2F *hPtMomPtMuFwd = new TH2F("hPtMomPtMuFwd", Form(";#it{p}_{T}^{%s};#it{p}_{T}^{#mu^{+}}", particleTitle[particle].c_str()), 100, 0, 10, 100, 0, 10);
    TH2F *hPtMomRapMom = new TH2F("hPtMomRapMom", Form(";#it{y}^{%s};#it{p}_{T}^{%s}", particleTitle[particle].c_str(), particleTitle[particle].c_str()), 100, -5, 5, 100, 0, 10);

    
    Pythia8::Pythia pythia;

    pythia.readString("ProcessLevel:all = off");

    pythia.readString("443:onMode = off");
    pythia.readString("443:onIfMatch = 13 -13");
    pythia.readString("100443:onMode = off");
    pythia.readString("100443:onIfMatch = 13 -13");
    pythia.readString("223:onMode = off");
    pythia.readString("223:onIfMatch = 13 -13");

    pythia.init();

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) { 
        tree->GetEntry(i);

        TLorentzVector *mom = new TLorentzVector();
        TLorentzVector *daum = new TLorentzVector();
        TLorentzVector *daup = new TLorentzVector();

        mom->SetPtEtaPhiM(pTOnia, etaOnia, phiOnia, momMass);

        tree->GetEntry(i);
        pythia.event.reset();
        
        Pythia8::Particle momParticle(
        absPdg,           // id
        1,                // status
        0,                // mother1
        0,                // mother2
        0,                // daughter1
        0,                // daughter2
        0,                // col
        0,                // acol
        mom->Px(),
        mom->Py(),
        mom->Pz(),
        mom->E(),
        mom->M()
        );

        pythia.event.append(momParticle);

        pythia.moreDecays();

        bool foundMuPlus = false; 
        bool foundMuMinus = false; 
        for (int j = 0; j < pythia.event.size(); ++j) { 
            const auto &part = pythia.event[j]; 
            if (part.id() == -13) { 
                foundMuPlus = true;
                daum->SetPxPyPzE(part.px(), part.py(), part.pz(), part.e());
            }

            if (part.id() == 13) { 
                foundMuMinus = true; 
                daup->SetPxPyPzE(part.px(), part.py(), part.pz(), part.e());
            }

            if (foundMuPlus && foundMuMinus) {
                hOpAngle->Fill(TMath::Cos(daup->Angle(daum->Vect())));
                hPt->Fill(mom->Pt());
                hRap->Fill(mom->Rapidity());
                hPtMu->Fill(daup->Pt());
                hRapMu->Fill(daup->Rapidity());
                hPtMomRapMom->Fill(mom->Rapidity(), mom->Pt());
                
                if (mom->Rapidity() > 2.5 && mom->Rapidity() < 4) {
                    hOpAngleFwd->Fill(TMath::Cos(daup->Angle(daum->Vect())));
                    hPtFwd->Fill(mom->Pt());
                    hRapFwd->Fill(mom->Rapidity());
                    hPtMuFwd->Fill(daup->Pt());
                    hRapMuFwd->Fill(daup->Rapidity());
                    hPtMomPtMuFwd->Fill(mom->Pt(), daup->Pt());
                }
            }
        }
    }

    TCanvas *canvas = new TCanvas("canvas", "", 1200, 1200);
    canvas->Divide(3,3);
    canvas->cd(1);
    hPtMomPtMuFwd->Draw("COLZ");
    canvas->cd(2);
    hPt->Draw("H");
    hPtFwd->Draw("H SAME");
    canvas->cd(3);
    hRap->Draw("H");
    hRapFwd->Draw("H SAME");
    canvas->cd(4);
    hOpAngle->Draw("H");
    hOpAngleFwd->Draw("H SAME");
    canvas->cd(5);
    hPtMu->Draw("H");
    hPtMuFwd->Draw("H SAME");
    canvas->cd(6);
    hRapMu->Draw("H");
    hRapMuFwd->Draw("H SAME");
    canvas->cd(7);
    hPtMomRapMom->Draw("COLZ");


    TFile *fOut = new TFile(Form("sim_output_%s.root", particleName[particle].c_str()), "RECREATE");
    hPtMomPtMuFwd->Write();
    hPt->Write();
    hPtFwd->Write();
    hRap->Write();
    hRapFwd->Write();
    hOpAngle->Write();
    hOpAngleFwd->Write();
    hPtMu->Write();
    hPtMuFwd->Write();
    hRapMu->Write();
    hRapMuFwd->Write();
    fOut->Close();

}