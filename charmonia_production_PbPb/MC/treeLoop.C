#include "TTree.h"
#include "TH1D.h"
#include "Math/Vector4D.h"
#include "TDatabasePDG.h"
#include <iostream>

// Process Generated Tree

void ProcessGeneratedTree(TTree* treeGen, 
                        TH1D* histPtJpsiGen,
                        TH1D* histRapJpsiGen,
                        TH1D* histCentrJpsiGen,
                        TH1D* histFromFuncPtRatio,
                        TH1D* histFromFuncRapRatio,
                        int iter,
                        bool isPbPb,
                        double massJpsi) {
    
    UInt_t fMcDecision;
    Float_t fImpactParameter;
    Float_t fPtMC1, fPtMC2;
    Float_t fEtaMC1, fEtaMC2;
    Float_t fPhiMC1, fPhiMC2;
    
    treeGen->SetBranchAddress("fMcDecision", &fMcDecision);
    if (isPbPb) {
        treeGen->SetBranchAddress("fImpactParameter", &fImpactParameter);
    }
    treeGen->SetBranchAddress("fPtMC1", &fPtMC1);
    treeGen->SetBranchAddress("fEtaMC1", &fEtaMC1);
    treeGen->SetBranchAddress("fPhiMC1", &fPhiMC1);
    treeGen->SetBranchAddress("fPtMC2", &fPtMC2);
    treeGen->SetBranchAddress("fEtaMC2", &fEtaMC2);
    treeGen->SetBranchAddress("fPhiMC2", &fPhiMC2);
    
    Long64_t nEntries = treeGen->GetEntries();
    std::cout << "Processing " << nEntries << " generated entries" << std::endl;
    
    int filled = 0;
    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {
        treeGen->GetEntry(iEntry);
        
        if (iEntry % 100000 == 0) {
            std::cout << "Processing entry " << iEntry << "/" << nEntries << std::endl;
        }
        
        if (fMcDecision < 1) continue;
        
        if (fPtMC2 == -999.0 && fPhiMC2 == -999.0) {
            ROOT::Math::PtEtaPhiMVector vecJpsiGen(fPtMC1, fEtaMC1, fPhiMC1, massJpsi);
            
            double rap = vecJpsiGen.Rapidity();
            if (TMath::Abs(rap) > 4.0 || TMath::Abs(rap) < 2.5) continue;
            if (vecJpsiGen.Pt() > 20.0) continue;
            
            if (iter == 0) {
                histPtJpsiGen->Fill(vecJpsiGen.Pt());
                histRapJpsiGen->Fill(-rap);
                filled++;
            } else {
                int binPt = histFromFuncPtRatio->FindBin(vecJpsiGen.Pt());
                int binRap = histFromFuncRapRatio->FindBin(-rap);
                double weightPt = histFromFuncPtRatio->GetBinContent(binPt);
                double weightRap = histFromFuncRapRatio->GetBinContent(binRap);
                double weightTot = weightPt * weightRap;
                
                histPtJpsiGen->Fill(vecJpsiGen.Pt(), weightTot);
                histRapJpsiGen->Fill(-rap, weightTot);
                filled++;
            }
            
            if (isPbPb) {
                if (fImpactParameter < 5.625) {
                    histCentrJpsiGen->AddBinContent(1);
                } else if (fImpactParameter < 8.375) {
                    histCentrJpsiGen->AddBinContent(2);
                } else if (fImpactParameter < 10.625) {
                    histCentrJpsiGen->AddBinContent(3);
                } else if (fImpactParameter <= 13.875) {
                    histCentrJpsiGen->AddBinContent(4);
                }
            }
        }
    }
    
    std::cout << "Filled " << filled << " entries in generated histograms" << std::endl;
}

// Process Reconstructed Tree
void ProcessReconstructedTree(TTree* treeRec,
                            TH1D* histPtJpsiRec,
                            TH1D* histRapJpsiRec,
                            TH1D* histCentrJpsiRec,
                            TH1D* histFromFuncPtRatio,
                            TH1D* histFromFuncRapRatio,
                            int iter,
                            bool isPbPb,
                            double massMu,
                            double ptCut) {
    
    UInt_t fMcDecision;
    Float_t fMass, fPt, fEta, fPhi;
    Float_t fCentFT0C;
    Float_t fPtMC1, fPtMC2;
    Float_t fEtaMC1, fEtaMC2;
    Float_t fPhiMC1, fPhiMC2;
    Float_t fPt1, fPt2;
    Float_t fEta1, fEta2;
    Float_t fPhi1, fPhi2;
    
    treeRec->SetBranchAddress("fMcDecision", &fMcDecision);
    treeRec->SetBranchAddress("fMass", &fMass);
    treeRec->SetBranchAddress("fPt", &fPt);
    treeRec->SetBranchAddress("fEta", &fEta);
    treeRec->SetBranchAddress("fPhi", &fPhi);
    treeRec->SetBranchAddress("fCentFT0C", &fCentFT0C);
    treeRec->SetBranchAddress("fPtMC1", &fPtMC1);
    treeRec->SetBranchAddress("fEtaMC1", &fEtaMC1);
    treeRec->SetBranchAddress("fPhiMC1", &fPhiMC1);
    treeRec->SetBranchAddress("fPtMC2", &fPtMC2);
    treeRec->SetBranchAddress("fEtaMC2", &fEtaMC2);
    treeRec->SetBranchAddress("fPhiMC2", &fPhiMC2);
    treeRec->SetBranchAddress("fPt1", &fPt1);
    treeRec->SetBranchAddress("fEta1", &fEta1);
    treeRec->SetBranchAddress("fPhi1", &fPhi1);
    treeRec->SetBranchAddress("fPt2", &fPt2);
    treeRec->SetBranchAddress("fEta2", &fEta2);
    treeRec->SetBranchAddress("fPhi2", &fPhi2);
    
    Long64_t nEntries = treeRec->GetEntries();
    std::cout << "Processing " << nEntries << " reconstructed entries" << std::endl;
    
    int filled = 0;
    for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {
        treeRec->GetEntry(iEntry);
        
        if (iEntry % 100000 == 0) {
            std::cout << "Processing entry " << iEntry << "/" << nEntries << std::endl;
        }
        
        if (fMcDecision < 1) continue;
        
        ROOT::Math::PtEtaPhiMVector vecMuGen1(fPtMC1, fEtaMC1, fPhiMC1, massMu);
        ROOT::Math::PtEtaPhiMVector vecMuGen2(fPtMC2, fEtaMC2, fPhiMC2, massMu);
        ROOT::Math::PtEtaPhiMVector vecJpsiGen = vecMuGen1 + vecMuGen2;
        ROOT::Math::PtEtaPhiMVector vecJpsiRec(fPt, fEta, fPhi, fMass);
        
        double rap = vecJpsiRec.Rapidity();
        if (TMath::Abs(rap) > 4.0 || TMath::Abs(rap) < 2.5) continue;
        if (fPt1 < ptCut || fPt2 < ptCut) continue;
        if (fPt > 20.0) continue;
        
        if (iter == 0) {
            histPtJpsiRec->Fill(vecJpsiRec.Pt());
            histRapJpsiRec->Fill(-rap);
            filled++;
        } else {
            int binPt = histFromFuncPtRatio->FindBin(vecJpsiGen.Pt());
            int binRap = histFromFuncRapRatio->FindBin(-vecJpsiGen.Rapidity());
            double weightPt = histFromFuncPtRatio->GetBinContent(binPt);
            double weightRap = histFromFuncRapRatio->GetBinContent(binRap);
            double weightTot = weightPt * weightRap;
            
            histPtJpsiRec->Fill(vecJpsiRec.Pt(), weightTot);
            histRapJpsiRec->Fill(-rap, weightTot);
            filled++;
        }
        
        if (isPbPb) {
            histCentrJpsiRec->Fill(fCentFT0C);
        }
    }
    
    std::cout << "Filled " << filled << " entries in reconstructed histograms" << std::endl;
}
