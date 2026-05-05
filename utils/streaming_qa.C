void SetLegend(TLegend *);
void SetHist(TH1D *hist, Color_t color, int size) {
    hist->SetMarkerStyle(size);
    hist->SetMarkerColor(color);
    hist->SetMarkerSize(0.5);
    hist->SetLineColor(color);
    hist->Scale(1./hist->Integral());
}
TCanvas* DoPlot(TH1D *hist, string canvasName, string varName) {
    int firstBin = hist->FindFirstBinAbove(0.0);
    int lastBin  = hist->FindLastBinAbove(0.0);
    double xRangeMin = hist->GetBinLowEdge(firstBin);
    double xRangeMax = hist->GetBinCenter(lastBin) + hist->GetBinWidth(lastBin);

    TCanvas *canvas = new TCanvas(canvasName.c_str(), "", 800, 600);
    gPad->SetLogy(true);
    SetHist(hist, kBlack, 20);
    hist->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax);
    hist->SetTitle("");
    hist->SetName(canvasName.c_str());
    hist->Draw("EP");

    return canvas;
}
TCanvas* DoRatioPlot(TH1D *histSum, std::vector<TH1D*> hist, string canvasName, string varName) {
    const int nRuns = int(hist.size());
    gStyle->SetPalette(kRainBow);

    int firstBin = histSum->FindFirstBinAbove(0.0);
    int lastBin  = histSum->FindLastBinAbove(0.0);
    double xRangeMin = histSum->GetBinLowEdge(firstBin);
    double xRangeMax = histSum->GetBinCenter(lastBin) + histSum->GetBinWidth(lastBin);

    TCanvas *canvas = new TCanvas(canvasName.c_str(), "", 600, 800);

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);

    pad1->SetBottomMargin(0);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);

    pad1->Draw();
    pad2->Draw();

    //gPad->BuildLegend();
    pad1->cd();
    gPad->SetLogy(true);
    SetHist(histSum, kBlack, 20);
    histSum->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax);
    histSum->SetTitle("");
    histSum->SetName("Sum");
    histSum->Draw("EP");
    for (int iRun = 0;iRun < nRuns;iRun++) {
        hist[iRun]->Scale(1./hist[iRun]->Integral());
        hist[iRun]->Draw("SAME PLC PMC");
    }
    gPad->BuildLegend();

    pad2->cd();
    std::vector<TH1D*> histRatio(nRuns);
    for (int iRun = 0;iRun < nRuns;iRun++) {
        histRatio[iRun] = (TH1D*) hist[iRun]->Clone(Form("histRatio_%d", iRun));

        if (iRun == 0) {
            histRatio[iRun]->SetTitle("");
            histRatio[iRun]->GetYaxis()->SetTitle("Data/MC");
            histRatio[iRun]->GetXaxis()->SetTitle(varName.c_str());
            histRatio[iRun]->GetXaxis()->SetLabelSize(0.08);
            histRatio[iRun]->GetXaxis()->SetTitleSize(0.1);
            histRatio[iRun]->GetYaxis()->SetNdivisions(505);
            histRatio[iRun]->GetYaxis()->SetTitleSize(0.1);
            histRatio[iRun]->GetYaxis()->SetLabelSize(0.08);
        }

        histRatio[iRun]->SetTitle("");
        histRatio[iRun]->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax);
        histRatio[iRun]->Divide(histSum);
        histRatio[iRun]->GetYaxis()->SetRangeUser(0.5, 1.5);
        histRatio[iRun]->Draw("SAME PLC PMC");
    }
    TLine *line = new TLine(xRangeMin, 1, xRangeMax, 1);
    line->SetLineStyle(2);
    line->Draw("SAME");

    canvas->cd();

    return canvas;
}

void streaming_qa(string process = "table-maker") {
    gStyle->SetPalette(kRainBow);
    string pathIn = "/Users/lucamicheletti/cernbox/DQ_streaming/train_666801";
    string muonCutTm = "muonMinimalCuts10SigmaPDCA";
    string muonCutTr = "muonQualityCutsMUONStandalone";

    TH1D *hSumVtxZ, *hSumCentFT0C, *hSumCentFT0M;
    TH2D *hSumPsi2A_CentFT0C, *hSumPsi2B_CentFT0C, *hSumPsi2C_CentFT0C;
    TH1D *hSumPt, *hSumEta, *hSumPhi, *hSumpdca, *hSumRAtAbsorberEnd, *hSumChi2MCHMID, *hSumChi2MCHMFT;

    if (process == "table-maker") {
        std::vector<int> runList;
        int entry;
        std::ifstream fInRunList(Form("%s/run_list.txt", pathIn.c_str()));
        while (fInRunList >> entry) {runList.push_back(entry);}
        const int nRuns = int(runList.size());
        std::cout << nRuns << std::endl;

        int iRun = 0;

        TH2D *hCountersVsRuns = new TH2D("hCountersVsRuns", ";Run;Counter", nRuns, 0, nRuns, 6, 0, 6);
        hCountersVsRuns->GetYaxis()->SetBinLabel(1, "BCs TVX");
        hCountersVsRuns->GetYaxis()->SetBinLabel(2, "BCs TVX after BC cuts");
        hCountersVsRuns->GetYaxis()->SetBinLabel(3, "Collisions All");
        hCountersVsRuns->GetYaxis()->SetBinLabel(4, "Collisions TVX");
        hCountersVsRuns->GetYaxis()->SetBinLabel(5, "Collisions Acc");
        hCountersVsRuns->GetYaxis()->SetBinLabel(6, "Collisions after cuts");

        std::vector<TH1D*> hVtxZ(nRuns);
        std::vector<TH1D*> hCentFT0C(nRuns);
        std::vector<TH1D*> hCentFT0M(nRuns);
        std::vector<TH2D*> hPsi2A_CentFT0C(nRuns);
        std::vector<TH2D*> hPsi2B_CentFT0C(nRuns);
        std::vector<TH2D*> hPsi2C_CentFT0C(nRuns);
        std::vector<TH1D*> hPt(nRuns);
        std::vector<TH1D*> hEta(nRuns);
        std::vector<TH1D*> hPhi(nRuns);
        std::vector<TH1D*> hpdca(nRuns);
        std::vector<TH1D*> hRAtAbsorberEnd(nRuns);
        std::vector<TH1D*> hChi2MCHMID(nRuns);
        std::vector<TH1D*> hChi2MCHMFT(nRuns);
        for (auto run : runList) {
            std::cout << "-------- Processing run " << run << " --------" << std::endl;
            TFile *fIn = new TFile(Form("%s/%d/AnalysisResults.root", pathIn.c_str(), run));

            TH1D *hCounterTVX = (TH1D*) fIn->Get("eventselection-run3/luminosity/hCounterTVX");
            TH1D *hCounterTVXafterBCcuts = (TH1D*) fIn->Get("eventselection-run3/luminosity/hCounterTVXafterBCcuts");
            TH1D *hColCounterAll = (TH1D*) fIn->Get("eventselection-run3/eventselection/hColCounterAll");
            TH1D *hColCounterTVX = (TH1D*) fIn->Get("eventselection-run3/eventselection/hColCounterTVX");
            TH1D *hColCounterAcc = (TH1D*) fIn->Get("eventselection-run3/eventselection/hColCounterAcc");

            TList *hlistTm = (TList*) fIn->Get("table-maker/output");
            TList *listTmEvents = (TList*) hlistTm->FindObject("Event_AfterCuts");
            hVtxZ[iRun] = (TH1D*) listTmEvents->FindObject("VtxZ");
            hCentFT0C[iRun] = (TH1D*) listTmEvents->FindObject("CentFT0C");
            hCentFT0M[iRun] = (TH1D*) listTmEvents->FindObject("CentFT0M");
            hPsi2A_CentFT0C[iRun] = (TH2D*) listTmEvents->FindObject("Psi2A_CentFT0C");
            hPsi2B_CentFT0C[iRun] = (TH2D*) listTmEvents->FindObject("Psi2B_CentFT0C");
            hPsi2C_CentFT0C[iRun] = (TH2D*) listTmEvents->FindObject("Psi2C_CentFT0C");

            hCountersVsRuns->SetBinContent(iRun+1, 1, hCounterTVX->GetEntries());
            hCountersVsRuns->SetBinContent(iRun+1, 2, hCounterTVXafterBCcuts->GetEntries());
            hCountersVsRuns->SetBinContent(iRun+1, 3, hColCounterAll->GetEntries());
            hCountersVsRuns->SetBinContent(iRun+1, 4, hColCounterTVX->GetEntries());
            hCountersVsRuns->SetBinContent(iRun+1, 5, hColCounterAcc->GetEntries());
            hCountersVsRuns->SetBinContent(iRun+1, 6, hVtxZ[iRun]->GetEntries());
            hCountersVsRuns->GetXaxis()->SetBinLabel(iRun+1, Form("%d", run));

            TList *listTmMuons = (TList*) hlistTm->FindObject(Form("Muons_%s", muonCutTm.c_str()));
            hPt[iRun] = (TH1D*) listTmMuons->FindObject("Pt"); hPt[iRun]->Rebin(10);
            hEta[iRun] = (TH1D*) listTmMuons->FindObject("Eta");
            hPhi[iRun] = (TH1D*) listTmMuons->FindObject("Phi");
            hpdca[iRun] = (TH1D*) listTmMuons->FindObject("pdca");
            hRAtAbsorberEnd[iRun] = (TH1D*) listTmMuons->FindObject("RAtAbsorberEnd");
            hChi2MCHMID[iRun] = (TH1D*) listTmMuons->FindObject("Chi2MCHMID");
            hChi2MCHMFT[iRun] = (TH1D*) listTmMuons->FindObject("Chi2MCHMFT");

            if (iRun == 0) {
                hSumVtxZ = (TH1D*) hVtxZ[iRun]->Clone("hSumVtxZ");
                hSumCentFT0C = (TH1D*) hCentFT0C[iRun]->Clone("hSumCentFT0C");
                hSumCentFT0M = (TH1D*) hCentFT0M[iRun]->Clone("hSumCentFT0M");
                hSumPsi2A_CentFT0C = (TH2D*) hPsi2A_CentFT0C[iRun]->Clone("hSumPsi2A_CentFT0C");
                hSumPsi2B_CentFT0C = (TH2D*) hPsi2B_CentFT0C[iRun]->Clone("hSumPsi2B_CentFT0C");
                hSumPsi2C_CentFT0C = (TH2D*) hPsi2C_CentFT0C[iRun]->Clone("hSumPsi2C_CentFT0C");
                hSumPt = (TH1D*) hPt[iRun]->Clone("hSumPt");
                hSumEta = (TH1D*) hEta[iRun]->Clone("hSumEta");
                hSumPhi = (TH1D*) hPhi[iRun]->Clone("hSumPhi");
                hSumpdca = (TH1D*) hpdca[iRun]->Clone("hSumpdca");
                hSumRAtAbsorberEnd = (TH1D*) hRAtAbsorberEnd[iRun]->Clone("hSumRAtAbsorberEnd");
                hSumChi2MCHMID = (TH1D*) hChi2MCHMID[iRun]->Clone("hSumChi2MCHMID");
                hSumChi2MCHMFT = (TH1D*) hChi2MCHMFT[iRun]->Clone("hSumChi2MCHMFT");
            } else {
                hSumVtxZ->Add(hVtxZ[iRun]);
                hSumCentFT0C->Add(hCentFT0C[iRun]);
                hSumCentFT0M->Add(hCentFT0M[iRun]);
                hSumPsi2A_CentFT0C->Add(hPsi2A_CentFT0C[iRun]);
                hSumPsi2B_CentFT0C->Add(hPsi2B_CentFT0C[iRun]);
                hSumPsi2C_CentFT0C->Add(hPsi2C_CentFT0C[iRun]);
                hSumPt->Add(hPt[iRun]);
                hSumEta->Add(hEta[iRun]);
                hSumPhi->Add(hPhi[iRun]);
                hSumpdca->Add(hpdca[iRun]);
                hSumRAtAbsorberEnd->Add(hRAtAbsorberEnd[iRun]);
                hSumChi2MCHMID->Add(hChi2MCHMID[iRun]);
                hSumChi2MCHMFT->Add(hChi2MCHMFT[iRun]);
            }

            hVtxZ[iRun]->SetTitle(Form("%d", run));
            hCentFT0C[iRun]->SetTitle(Form("%d", run));
            hCentFT0M[iRun]->SetTitle(Form("%d", run));
            hPsi2A_CentFT0C[iRun]->SetTitle(Form("%d", run));
            hPsi2B_CentFT0C[iRun]->SetTitle(Form("%d", run));
            hPsi2C_CentFT0C[iRun]->SetTitle(Form("%d", run));
            hPt[iRun]->SetTitle(Form("%d", run));
            hEta[iRun]->SetTitle(Form("%d", run));
            hPhi[iRun]->SetTitle(Form("%d", run));
            hpdca[iRun]->SetTitle(Form("%d", run));
            hRAtAbsorberEnd[iRun]->SetTitle(Form("%d", run));
            hChi2MCHMID[iRun]->SetTitle(Form("%d", run));
            hChi2MCHMFT[iRun]->SetTitle(Form("%d", run));

            iRun++;
        }

        TCanvas *canvasNormalization = new TCanvas("canvasNormalization", "", 800, 600);
        gStyle->SetOptStat(false);
        hCountersVsRuns->Draw("COLZ");

        TCanvas *canvasVtxZ = DoRatioPlot(hSumVtxZ, hVtxZ, "canvasVtxZ", "vtx_{Z} (cm)");
        TCanvas *canvasCentFT0C = DoRatioPlot(hSumCentFT0C, hCentFT0C, "canvasCentFT0C", "Centr FT0C (%)");
        TCanvas *canvasCentFT0M = DoRatioPlot(hSumCentFT0M, hCentFT0M, "canvasCentFT0M", "Centr FT0M (%)");

        // Rebin EP distribution
        hSumPsi2A_CentFT0C->RebinX(5);
        hSumPsi2B_CentFT0C->RebinX(5);
        hSumPsi2C_CentFT0C->RebinX(5);
        for (int iRun = 0;iRun < nRuns;iRun++) {
            hPsi2A_CentFT0C[iRun]->RebinX(5);
            hPsi2B_CentFT0C[iRun]->RebinX(5);
            hPsi2C_CentFT0C[iRun]->RebinX(5);
        }

        TCanvas *canvasPsi2A_CentFT0C = new TCanvas("canvasPsi2A_CentFT0C", "", 1800, 1800);
        canvasPsi2A_CentFT0C->Divide(3, 3);

        for (int iCent = 0;iCent < 9;iCent++) {
            canvasPsi2A_CentFT0C->cd(iCent+1);
            TH1D *hTmpSumPsi2A_CentFT0C = (TH1D*) hSumPsi2A_CentFT0C->ProjectionY(Form("hTmpSumPsi2A_CentFT0C_%d_%d", iRun, iCent), iCent+1, iCent+1);
            SetHist(hTmpSumPsi2A_CentFT0C, kBlack, 20);
            hTmpSumPsi2A_CentFT0C->Draw("EP");
            for (int iRun = 0;iRun < nRuns;iRun++) {
                TH1D *hTmpPsi2A_CentFT0C = (TH1D*) hPsi2A_CentFT0C[iRun]->ProjectionY(Form("hTmpPsi2A_CentFT0C_%d_%d", iRun, iCent), iCent+1, iCent+1);
                hTmpPsi2A_CentFT0C->Scale(1./hTmpPsi2A_CentFT0C->Integral());
                hTmpPsi2A_CentFT0C->Draw("SAME PLC PMC");
            }
        }

        TCanvas *canvasPsi2B_CentFT0C = new TCanvas("canvasPsi2B_CentFT0C", "", 1800, 1800);
        canvasPsi2B_CentFT0C->Divide(3, 3);

        for (int iCent = 0;iCent < 9;iCent++) {
            canvasPsi2B_CentFT0C->cd(iCent+1);
            TH1D *hTmpSumPsi2B_CentFT0C = (TH1D*) hSumPsi2B_CentFT0C->ProjectionY(Form("hTmpSumPsi2B_CentFT0C_%d_%d", iRun, iCent), iCent+1, iCent+1);
            SetHist(hTmpSumPsi2B_CentFT0C, kBlack, 20);
            hTmpSumPsi2B_CentFT0C->Draw("EP");
            for (int iRun = 0;iRun < nRuns;iRun++) {
                TH1D *hTmpPsi2B_CentFT0C = (TH1D*) hPsi2B_CentFT0C[iRun]->ProjectionY(Form("hTmpPsi2B_CentFT0C_%d_%d", iRun, iCent), iCent+1, iCent+1);
                hTmpPsi2B_CentFT0C->Scale(1./hTmpPsi2B_CentFT0C->Integral());
                hTmpPsi2B_CentFT0C->Draw("SAME PLC PMC");
            }
        }

        TCanvas *canvasPsi2C_CentFT0C = new TCanvas("canvasPsi2C_CentFT0C", "", 1800, 1800);
        canvasPsi2C_CentFT0C->Divide(3, 3);

        for (int iCent = 0;iCent < 9;iCent++) {
            canvasPsi2C_CentFT0C->cd(iCent+1);
            TH1D *hTmpSumPsi2C_CentFT0C = (TH1D*) hSumPsi2C_CentFT0C->ProjectionY(Form("hTmpSumPsi2C_CentFT0C_%d_%d", iRun, iCent), iCent+1, iCent+1);
            SetHist(hTmpSumPsi2C_CentFT0C, kBlack, 20);
            hTmpSumPsi2C_CentFT0C->Draw("EP");
            for (int iRun = 0;iRun < nRuns;iRun++) {
                TH1D *hTmpPsi2C_CentFT0C = (TH1D*) hPsi2C_CentFT0C[iRun]->ProjectionY(Form("hTmpPsi2C_CentFT0C_%d_%d", iRun, iCent), iCent+1, iCent+1);
                hTmpPsi2C_CentFT0C->Scale(1./hTmpPsi2C_CentFT0C->Integral());
                hTmpPsi2C_CentFT0C->Draw("SAME PLC PMC");
            }
        }

        TCanvas *canvasPt = DoRatioPlot(hSumPt, hPt, "canvasPt", "#it{p}_{T} (GeV/#it{c})");
        TCanvas *canvasEta = DoRatioPlot(hSumEta, hEta, "canvasEta", "#eta");
        TCanvas *canvasPhi = DoRatioPlot(hSumPhi, hPhi, "canvasPhi", "#varphi");
        TCanvas *canvaspdca = DoRatioPlot(hSumpdca, hpdca, "canvaspdca", "p x DCA");
        TCanvas *canvasRAtAbsorberEnd = DoRatioPlot(hSumRAtAbsorberEnd, hRAtAbsorberEnd, "canvasRAtAbsorberEnd", "R_{abs}");
        TCanvas *canvasChi2MCHMID = DoRatioPlot(hSumChi2MCHMID, hChi2MCHMID, "canvasChi2MCHMID", "#chi^{2}_{MCH-MID}");
        TCanvas *canvasChi2MCHMFT = DoRatioPlot(hSumChi2MCHMFT, hChi2MCHMFT, "canvasChi2MCHMFT", "#chi^{2}_{MFT-MCH}");


        TFile *fOut = new TFile(Form("%s/streamingResults.root", pathIn.c_str()), "RECREATE");
        hCountersVsRuns->Write();
        canvasVtxZ->Write();
        canvasCentFT0C->Write();
        canvasCentFT0M->Write();
        canvasPsi2A_CentFT0C->Write();
        canvasPsi2B_CentFT0C->Write();
        canvasPsi2C_CentFT0C->Write();
        canvasPt->Write();
        canvasEta->Write();
        canvasPhi->Write();
        canvaspdca->Write();
        canvasRAtAbsorberEnd->Write();
        canvasChi2MCHMID->Write();
        canvasChi2MCHMFT->Write();
        fOut->Close();
    } else if (process == "table-reader") {
        TFile *fIn = new TFile(Form("%s/AnalysisResults.root", pathIn.c_str()));

        TList *hlistTrEvents = (TList*) fIn->Get("analysis-event-selection/output");
        TList *listTrEventsBeforeCuts = (TList*) hlistTrEvents->FindObject("Event_BeforeCuts");
        TH1D *hSumBeforeCutsVtxZ = (TH1D*) listTrEventsBeforeCuts->FindObject("VtxZ");

        TList *listTrEvents = (TList*) hlistTrEvents->FindObject("Event_AfterCuts");
        hSumVtxZ = (TH1D*) listTrEvents->FindObject("VtxZ");
        hSumCentFT0C = (TH1D*) listTrEvents->FindObject("CentFT0C");
        hSumCentFT0M = (TH1D*) listTrEvents->FindObject("CentFT0M");
        hSumPsi2A_CentFT0C = (TH2D*) listTrEvents->FindObject("Psi2A_CentFT0C");
        hSumPsi2B_CentFT0C = (TH2D*) listTrEvents->FindObject("Psi2B_CentFT0C");
        hSumPsi2C_CentFT0C = (TH2D*) listTrEvents->FindObject("Psi2C_CentFT0C");

        TList *hlistTrMuons = (TList*) fIn->Get("analysis-muon-selection/output");
        TList *listTrMuons = (TList*) hlistTrMuons->FindObject(Form("TrackMuon_%s", muonCutTr.c_str()));
        hSumPt = (TH1D*) listTrMuons->FindObject("Pt"); hSumPt->Rebin(10);
        hSumEta = (TH1D*) listTrMuons->FindObject("Eta");
        hSumPhi = (TH1D*) listTrMuons->FindObject("Phi");
        hSumpdca = (TH1D*) listTrMuons->FindObject("pdca");
        hSumRAtAbsorberEnd = (TH1D*) listTrMuons->FindObject("RAtAbsorberEnd");
        hSumChi2MCHMID = (TH1D*) listTrMuons->FindObject("Chi2MCHMID");
        hSumChi2MCHMFT = (TH1D*) listTrMuons->FindObject("Chi2MCHMFT");

        TList *hlistTrDimuons = (TList*) fIn->Get("analysis-same-event-pairing/output");
        TList *listTrDimuons = (TList*) hlistTrDimuons->FindObject(Form("PairsMuonSEPM_%s", muonCutTr.c_str()));
        THnSparseD *hSparseDimuon = (THnSparseD*) listTrDimuons -> FindObject("Mass_Pt_Cent_Rap");
        hSparseDimuon -> GetAxis(1) -> SetRangeUser(0., 20.);
        hSparseDimuon -> GetAxis(3) -> SetRangeUser(2.5, 4);
        TH1D *hDimuonMass = (TH1D*) hSparseDimuon -> Projection(0, "hDimuonMass");

        TH1D *hCounters = new TH1D("hCounters", ";Run;Counter", 2, 0, 2);
        hCounters->GetXaxis()->SetBinLabel(1, "Collisions before cuts");
        hCounters->GetXaxis()->SetBinLabel(2, "Collisions after cuts");
        hCounters->SetBinContent(1, hSumBeforeCutsVtxZ->GetEntries());
        hCounters->SetBinContent(2, hSumVtxZ->GetEntries());

        TCanvas *canvasVtxZ = DoPlot(hSumVtxZ, "canvasVtxZ", "vtx_{Z} (cm)");
        TCanvas *canvasCentFT0C = DoPlot(hSumCentFT0C, "canvasCentFT0C", "Centr FT0C (%)");
        TCanvas *canvasCentFT0M = DoPlot(hSumCentFT0M, "canvasCentFT0M", "Centr FT0M (%)");


        // Rebin EP distribution
        hSumPsi2A_CentFT0C->RebinX(5);
        hSumPsi2B_CentFT0C->RebinX(5);
        hSumPsi2C_CentFT0C->RebinX(5);

        TCanvas *canvasPsi2A_CentFT0C = new TCanvas("canvasPsi2A_CentFT0C", "", 1800, 1800);
        canvasPsi2A_CentFT0C->Divide(3, 3);

        for (int iCent = 0;iCent < 9;iCent++) {
            canvasPsi2A_CentFT0C->cd(iCent+1);
            TH1D *hTmpSumPsi2A_CentFT0C = (TH1D*) hSumPsi2A_CentFT0C->ProjectionY(Form("hTmpSumPsi2A_CentFT0C_%d", iCent), iCent+1, iCent+1);
            SetHist(hTmpSumPsi2A_CentFT0C, kBlack, 20);
            hTmpSumPsi2A_CentFT0C->Draw("EP");
        }

        TCanvas *canvasPsi2B_CentFT0C = new TCanvas("canvasPsi2B_CentFT0C", "", 1800, 1800);
        canvasPsi2B_CentFT0C->Divide(3, 3);

        for (int iCent = 0;iCent < 9;iCent++) {
            canvasPsi2B_CentFT0C->cd(iCent+1);
            TH1D *hTmpSumPsi2B_CentFT0C = (TH1D*) hSumPsi2B_CentFT0C->ProjectionY(Form("hTmpSumPsi2B_CentFT0C_%d", iCent), iCent+1, iCent+1);
            SetHist(hTmpSumPsi2B_CentFT0C, kBlack, 20);
            hTmpSumPsi2B_CentFT0C->Draw("EP");
        }

        TCanvas *canvasPsi2C_CentFT0C = new TCanvas("canvasPsi2C_CentFT0C", "", 1800, 1800);
        canvasPsi2C_CentFT0C->Divide(3, 3);

        for (int iCent = 0;iCent < 9;iCent++) {
            canvasPsi2C_CentFT0C->cd(iCent+1);
            TH1D *hTmpSumPsi2C_CentFT0C = (TH1D*) hSumPsi2C_CentFT0C->ProjectionY(Form("hTmpSumPsi2C_CentFT0C_%d", iCent), iCent+1, iCent+1);
            SetHist(hTmpSumPsi2C_CentFT0C, kBlack, 20);
            hTmpSumPsi2C_CentFT0C->Draw("EP");
        }

        TCanvas *canvasPt = DoPlot(hSumPt, "canvasPt", "#it{p}_{T} (GeV/#it{c})");
        TCanvas *canvasEta = DoPlot(hSumEta, "canvasEta", "#eta");
        TCanvas *canvasPhi = DoPlot(hSumPhi, "canvasPhi", "#varphi");
        TCanvas *canvaspdca = DoPlot(hSumpdca, "canvaspdca", "p x DCA");
        TCanvas *canvasRAtAbsorberEnd = DoPlot(hSumRAtAbsorberEnd, "canvasRAtAbsorberEnd", "R_{abs}");
        TCanvas *canvasChi2MCHMID = DoPlot(hSumChi2MCHMID, "canvasChi2MCHMID", "#chi^{2}_{MCH-MID}");

        TFile *fOut = new TFile(Form("%s/QA.root", pathIn.c_str()), "RECREATE");
        hCounters->Write();
        hDimuonMass->Write("hDimuonMass");
        canvasVtxZ->Write();
        canvasCentFT0C->Write();
        canvasCentFT0M->Write();
        canvasPsi2A_CentFT0C->Write();
        canvasPsi2B_CentFT0C->Write();
        canvasPsi2C_CentFT0C->Write();
        canvasPt->Write();
        canvasEta->Write();
        canvasPhi->Write();
        canvaspdca->Write();
        canvasRAtAbsorberEnd->Write();
        canvasChi2MCHMID->Write();
        fOut->Close();

    } else {
        std::cout << "-------- No process selected --------" << std::endl;
        return;
    }
}