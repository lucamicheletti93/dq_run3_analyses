void readAnalysis(string system = "OO") {
    const int nFiles = 9;

    string outDir = "/Users/lucamicheletti/cernbox/presentazioni/PWG_DQ/O2_DQ/2026/heavy_ion_normalization";
    string dirPathData = "normalization_" + system;
    string dirPathMc = "normalization_" + system + "_MC";

    string fInNames[nFiles] = {
        "eventSel8", "eventSel8NoSameBunch", "eventSel8Ft0CCentral", "eventSel8Ft0CSemiCentral", "eventSel8Ft0CAllCentral", "eventLightIonNoPileUp", "eventLightIonQuality", "eventSel8TriggerZNAZNC", "eventSel8TriggerZNAZNCNoPileUp"
    };
    string histNames[nFiles] = {
        "Sel8", "sel8 + NoPileUp", "Sel8 & Cen", "Sel8 & SCen", "Sel8 & (Cen | SCen)", "Sel8 & (Cen | SCen) & NoPileUp", "Sel8 & (Cen | SCen) & NoPileUp & ZvtxFT0vsPV", "Sel8 & (ZNA & ZNC)", "Sel8 & (ZNA & ZNC) & NoPileUp"
    };
    Color_t colors[nFiles] = {kBlack, kOrange+7, kBlue+1, kRed+1, kGreen+1, kMagenta, kAzure+4, kViolet, kPink};

    //**************************************************//

    TH1D *histDataCentFT0Cs[nFiles];
    TH1D *histDataEvSelSummary = new TH1D("histDataEvSelSummary","", 13+nFiles, 0, 13+nFiles);
    histDataEvSelSummary -> GetXaxis() -> SetBinLabel(1, "hColCounterAll");
    histDataEvSelSummary -> GetXaxis() -> SetBinLabel(2, "hColCounterTVX");
    histDataEvSelSummary -> GetXaxis() -> SetBinLabel(3, "hColCounterAcc");
    histDataEvSelSummary -> GetXaxis() -> SetBinLabel(4, "hCounterTVX");
    histDataEvSelSummary -> GetXaxis() -> SetBinLabel(5, "hCounterTVXZDC");
    histDataEvSelSummary -> GetXaxis() -> SetBinLabel(6, "hCounterTVXafterBCcuts");
    histDataEvSelSummary -> GetXaxis() -> SetBinLabel(7, "hCounterTVXZDCafterBCcuts");

    histDataEvSelSummary -> GetXaxis() -> SetBinLabel(8, "BC: Sel8");
    histDataEvSelSummary -> GetXaxis() -> SetBinLabel(9, "BC: Sel8 & Cen");
    histDataEvSelSummary -> GetXaxis() -> SetBinLabel(10, "BC: Sel8 & SCen");
    histDataEvSelSummary -> GetXaxis() -> SetBinLabel(11, "BC: Sel8 & (Cen | SCen)");
    histDataEvSelSummary -> GetXaxis() -> SetBinLabel(12, "BC: Sel8 & (ZNA & ZNC)");

    histDataEvSelSummary -> GetXaxis() -> SetBinLabel(13, "TM: before cuts");
    for (int iFile = 0;iFile < nFiles;iFile++) {
        histDataEvSelSummary -> GetXaxis() -> SetBinLabel(14+iFile, Form("TM: %s", histNames[iFile].c_str()));
    }

    for (int iFile = 0;iFile < nFiles;iFile++) {
        TFile *fIn = new TFile(Form("%s/AnalysisResults_%s.root", dirPathData.c_str(), fInNames[iFile].c_str()), "READ");
        // eventselection-run3
        histDataEvSelSummary -> SetBinContent(1, ((TH1D*) fIn -> Get("eventselection-run3/eventselection/hColCounterAll")) -> GetEntries());
        histDataEvSelSummary -> SetBinContent(2, ((TH1D*) fIn -> Get("eventselection-run3/eventselection/hColCounterTVX")) -> GetEntries());
        histDataEvSelSummary -> SetBinContent(3, ((TH1D*) fIn -> Get("eventselection-run3/eventselection/hColCounterAcc")) -> GetEntries());
        histDataEvSelSummary -> SetBinContent(4, ((TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVX")) -> GetEntries());
        histDataEvSelSummary -> SetBinContent(5, ((TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVXZDC")) -> GetEntries());
        histDataEvSelSummary -> SetBinContent(6, ((TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVXafterBCcuts")) -> GetEntries());
        histDataEvSelSummary -> SetBinContent(7, ((TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVXZDCafterBCcuts")) -> GetEntries());

        // table-maker
        TList *hlistStatistics = (TList*) fIn -> Get("table-maker/Statistics");
        histDataEvSelSummary -> SetBinContent(8, ((TH1D*) hlistStatistics -> FindObject("BcStats")) -> GetBinContent(1));
        histDataEvSelSummary -> SetBinContent(9, ((TH1D*) hlistStatistics -> FindObject("BcStats")) -> GetBinContent(2));
        histDataEvSelSummary -> SetBinContent(10, ((TH1D*) hlistStatistics -> FindObject("BcStats")) -> GetBinContent(3));
        histDataEvSelSummary -> SetBinContent(11, ((TH1D*) hlistStatistics -> FindObject("BcStats")) -> GetBinContent(4));
        histDataEvSelSummary -> SetBinContent(12, ((TH1D*) hlistStatistics -> FindObject("BcStats")) -> GetBinContent(5));

        TList *hlistOutput = (TList*) fIn -> Get("table-maker/output");
        TList *listBeforeCuts = (TList*) hlistOutput -> FindObject("Event_BeforeCuts");
        histDataEvSelSummary -> SetBinContent(13, ((TH1D*) listBeforeCuts -> FindObject("CentFT0C")) -> GetEntries());

        TList *listAfterCuts = (TList*) hlistOutput -> FindObject("Event_AfterCuts");
        histDataEvSelSummary -> SetBinContent(14+iFile, ((TH1D*) listAfterCuts -> FindObject("CentFT0C")) -> GetEntries());

        histDataCentFT0Cs[iFile] = (TH1D*) listAfterCuts -> FindObject("CentFT0C");
        histDataCentFT0Cs[iFile] -> SetMarkerStyle(20);
        histDataCentFT0Cs[iFile] -> SetMarkerSize(0.8);
        histDataCentFT0Cs[iFile] -> SetMarkerColor(colors[iFile]);
        histDataCentFT0Cs[iFile] -> SetLineColor(colors[iFile]);
        histDataCentFT0Cs[iFile] -> SetName(Form("Data, %s", histNames[iFile].c_str()));
        histDataCentFT0Cs[iFile] -> SetTitle("");
        histDataCentFT0Cs[iFile] -> Rebin(5);
        fIn -> Close();
    }


    TH1D *histDataRatios[nFiles];

    for (int iFile = 1;iFile < nFiles;iFile++) {
        histDataRatios[iFile] = (TH1D*) histDataCentFT0Cs[iFile] -> Clone(Form("%s / sel8", histDataCentFT0Cs[iFile]-> GetName()));
        histDataRatios[iFile] -> Divide(histDataCentFT0Cs[0]);
    }


    double yMax = histDataCentFT0Cs[0] -> GetMaximum();
    TCanvas *canvasData = new TCanvas("canvasData", "", 800, 1200);
    canvasData -> Divide(1, 2);

    canvasData -> cd(1);
    gPad -> SetLogy(true);
    gStyle -> SetOptStat(false);
    for (int iFile = 0;iFile < nFiles;iFile++) {
        histDataCentFT0Cs[iFile] -> GetYaxis() -> SetRangeUser(1, yMax + (0.5 * yMax));
        histDataCentFT0Cs[iFile] -> Draw("SAME");
    }

    gPad -> BuildLegend();

    canvasData -> cd(2);
    for (int iFile = 1;iFile < nFiles;iFile++) {
        histDataRatios[iFile] -> GetYaxis() -> SetRangeUser(0, 1.6);
        histDataRatios[iFile] -> Draw("SAME");
    }
    TLine *lineUnity = new TLine(0, 1, 100, 1);
    lineUnity -> SetLineColor(kGray+2);
    lineUnity -> SetLineStyle(kDashed);
    lineUnity -> Draw();

    gPad -> BuildLegend();

    //**************************************************//

    TH1D *histMcCentFT0Cs[nFiles];
    TH1D *histMcEvSelSummary = new TH1D("histMcEvSelSummary","", 13+nFiles, 0, 13+nFiles);
    histMcEvSelSummary -> GetXaxis() -> SetBinLabel(1, "hColCounterAll");
    histMcEvSelSummary -> GetXaxis() -> SetBinLabel(2, "hColCounterTVX");
    histMcEvSelSummary -> GetXaxis() -> SetBinLabel(3, "hColCounterAcc");
    histMcEvSelSummary -> GetXaxis() -> SetBinLabel(4, "hCounterTVX");
    histMcEvSelSummary -> GetXaxis() -> SetBinLabel(5, "hCounterTVXZDC");
    histMcEvSelSummary -> GetXaxis() -> SetBinLabel(6, "hCounterTVXafterBCcuts");
    histMcEvSelSummary -> GetXaxis() -> SetBinLabel(7, "hCounterTVXZDCafterBCcuts");

    histMcEvSelSummary -> GetXaxis() -> SetBinLabel(8, "BC: sel8");
    histMcEvSelSummary -> GetXaxis() -> SetBinLabel(9, "BC: sel8 & Cent");
    histMcEvSelSummary -> GetXaxis() -> SetBinLabel(10, "BC: sel8 & sCent");
    histMcEvSelSummary -> GetXaxis() -> SetBinLabel(11, "BC: sel8 & (Cent | sCent)");
    histMcEvSelSummary -> GetXaxis() -> SetBinLabel(12, "BC: sel8 & (ZNA & ZNC)");

    histMcEvSelSummary -> GetXaxis() -> SetBinLabel(13, "TM: before cuts");
    for (int iFile = 0;iFile < nFiles;iFile++) {
        histMcEvSelSummary -> GetXaxis() -> SetBinLabel(14+iFile, Form("TM: %s", histNames[iFile].c_str()));
    }
    for (int iFile = 0;iFile < nFiles;iFile++) {
        TFile *fIn = new TFile(Form("%s/AnalysisResults_%s.root", dirPathMc.c_str(), fInNames[iFile].c_str()), "READ");
        
        // eventselection-run3
        histMcEvSelSummary -> SetBinContent(1, ((TH1D*) fIn -> Get("eventselection-run3/eventselection/hColCounterAll")) -> GetEntries());
        histMcEvSelSummary -> SetBinContent(2, ((TH1D*) fIn -> Get("eventselection-run3/eventselection/hColCounterTVX")) -> GetEntries());
        histMcEvSelSummary -> SetBinContent(3, ((TH1D*) fIn -> Get("eventselection-run3/eventselection/hColCounterAcc")) -> GetEntries());
        histMcEvSelSummary -> SetBinContent(4, ((TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVX")) -> GetEntries());
        histMcEvSelSummary -> SetBinContent(5, ((TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVXZDC")) -> GetEntries());
        histMcEvSelSummary -> SetBinContent(6, ((TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVXafterBCcuts")) -> GetEntries());
        histMcEvSelSummary -> SetBinContent(7, ((TH1D*) fIn -> Get("eventselection-run3/luminosity/hCounterTVXZDCafterBCcuts")) -> GetEntries());

        // table-maker
        TList *hlistStatistics = (TList*) fIn -> Get("table-maker/Statistics");
        histMcEvSelSummary -> SetBinContent(8, ((TH1D*) hlistStatistics -> FindObject("BcStats")) -> GetBinContent(1));
        histMcEvSelSummary -> SetBinContent(9, ((TH1D*) hlistStatistics -> FindObject("BcStats")) -> GetBinContent(2));
        histMcEvSelSummary -> SetBinContent(10, ((TH1D*) hlistStatistics -> FindObject("BcStats")) -> GetBinContent(3));
        histMcEvSelSummary -> SetBinContent(11, ((TH1D*) hlistStatistics -> FindObject("BcStats")) -> GetBinContent(4));
        histMcEvSelSummary -> SetBinContent(12, ((TH1D*) hlistStatistics -> FindObject("BcStats")) -> GetBinContent(5));

        TList *hlistOutput = (TList*) fIn -> Get("table-maker/output");

        TList *listBeforeCuts = (TList*) hlistOutput -> FindObject("Event_BeforeCuts");
        histMcEvSelSummary -> SetBinContent(13, ((TH1D*) listBeforeCuts -> FindObject("CentFT0C")) -> GetEntries());

        TList *listAfterCuts = (TList*) hlistOutput -> FindObject("Event_AfterCuts");
        histMcEvSelSummary -> SetBinContent(14+iFile, ((TH1D*) listAfterCuts -> FindObject("CentFT0C")) -> GetEntries());

        histMcCentFT0Cs[iFile] = (TH1D*) listAfterCuts -> FindObject("CentFT0C");
        histMcCentFT0Cs[iFile] -> SetMarkerStyle(24);
        histMcCentFT0Cs[iFile] -> SetMarkerSize(0.8);
        histMcCentFT0Cs[iFile] -> SetMarkerColor(colors[iFile]);
        histMcCentFT0Cs[iFile] -> SetLineColor(colors[iFile]);
        histMcCentFT0Cs[iFile] -> SetName(Form("MC, %s", histNames[iFile].c_str()));
        histMcCentFT0Cs[iFile] -> SetTitle("");
        histMcCentFT0Cs[iFile] -> GetXaxis() -> SetRangeUser(0, 100);
        histMcCentFT0Cs[iFile] -> Rebin(5);
        fIn -> Close();
    }


    TH1D *histMcRatios[nFiles];
    for (int iFile = 1;iFile < nFiles;iFile++) {
        histMcRatios[iFile] = (TH1D*) histMcCentFT0Cs[iFile] -> Clone(Form("%s / sel8", histMcCentFT0Cs[iFile]-> GetName()));
        histMcRatios[iFile] -> Divide(histMcCentFT0Cs[0]);
    }

    yMax = histMcCentFT0Cs[0] -> GetMaximum();
    TCanvas *canvasMc = new TCanvas("canvasMc", "", 800, 1200);
    canvasMc -> Divide(1, 2);

    canvasMc -> cd(1);
    gPad -> SetLogy(true);
    gStyle -> SetOptStat(false);
    for (int iFile = 0;iFile < nFiles;iFile++) {
        histMcCentFT0Cs[iFile] -> GetYaxis() -> SetRangeUser(1, yMax + (0.5 * yMax));
        histMcCentFT0Cs[iFile] -> Draw("SAME");
    }

    gPad -> BuildLegend();

    canvasMc -> cd(2);
    for (int iFile = 1;iFile < nFiles;iFile++) {
        histMcRatios[iFile] -> GetYaxis() -> SetRangeUser(0, 1.6);
        histMcRatios[iFile] -> Draw("SAME");
    }
    lineUnity -> Draw();

    gPad -> BuildLegend();


    TCanvas *canvasCompDataMc = new TCanvas("canvasCompDataMc", "", 1200, 1200);
    canvasCompDataMc -> Divide(3, 3);

    for (int iFile = 1;iFile < nFiles;iFile++) {
        canvasCompDataMc -> cd(iFile); 
        histDataRatios[iFile] -> SetTitle(Form("%s / sel8", histNames[iFile].c_str()));
        histDataRatios[iFile] -> Draw("SAME");
        histMcRatios[iFile] -> Draw("SAME");
        lineUnity -> Draw();
    }

    TCanvas *canvasDataEvSelSummary = new TCanvas("canvasDataEvSelSummary", "", 800, 800);
    gPad -> SetBottomMargin(0.5);
    gPad -> SetRightMargin(0.1);
    histDataEvSelSummary -> LabelsOption("v", "X");
    histDataEvSelSummary -> Draw("EP");

    // Ratios
    const int nRatios = 5;
    string numNames[] = {"TM: Sel8", "TM: Sel8 & Cen", "TM: Sel8 & SCen", "TM: Sel8 & (Cen | SCen)", "TM: Sel8 & (ZNA & ZNC)"};
    string denNames[] = {"BC: Sel8", "BC: Sel8 & Cen", "BC: Sel8 & SCen", "BC: Sel8 & (Cen | SCen)", "BC: Sel8 & (ZNA & ZNC)"};

    TH1D *histDataEvSelRatioSummary = new TH1D("histDataEvSelRatioSummary", "", nRatios, 0, nRatios);
    histDataEvSelRatioSummary -> SetMarkerStyle(20);
    histDataEvSelRatioSummary -> SetMarkerColor(kBlack);
    histDataEvSelRatioSummary -> SetLineColor(kBlack);
    for (int iBin = 0;iBin < nRatios;iBin++) {
        double num = histDataEvSelSummary -> GetBinContent(histDataEvSelSummary -> GetXaxis() -> FindBin(numNames[iBin].c_str()));
        double den = histDataEvSelSummary -> GetBinContent(histDataEvSelSummary -> GetXaxis() -> FindBin(denNames[iBin].c_str()));
        double errRelNum = 1 / TMath::Sqrt(num);
        double errRelDen = 1 / TMath::Sqrt(den);

        histDataEvSelRatioSummary -> GetXaxis() -> SetBinLabel(iBin+1, Form("%s / %s", numNames[iBin].c_str(), denNames[iBin].c_str()));
        histDataEvSelRatioSummary -> SetBinContent(iBin+1, num / den);
        histDataEvSelRatioSummary -> SetBinError(iBin+1, (num / den) * TMath::Sqrt(errRelNum*errRelNum + errRelDen*errRelDen));
        std::cout << Form("%s / %s = ", numNames[iBin].c_str(), denNames[iBin].c_str()) << num / den << std::endl;
    }

    TCanvas *canvasDataEvSelRatioSummary = new TCanvas("canvasDataEvSelRatioSummary", "", 800, 800);
    gPad -> SetBottomMargin(0.6);
    gPad -> SetRightMargin(0.1);
    histDataEvSelRatioSummary -> GetYaxis() -> SetRangeUser(0, 1.2);
    histDataEvSelRatioSummary -> LabelsOption("v", "X");
    histDataEvSelRatioSummary -> Draw("EP");

    TLine *lineRatioUnity = new TLine(0, 1, nRatios, 1);
    lineRatioUnity -> SetLineColor(kGray+2);
    lineRatioUnity -> SetLineStyle(kDashed);
    lineRatioUnity -> Draw();

    // Plotting the impact of the pileup
    TCanvas *canvasPileUp = new TCanvas("canvasPileUp", "", 800, 600);
    histDataRatios[1] -> Draw("EP");
    histMcRatios[1] -> Draw("EP SAME");
    lineRatioUnity -> Draw("SAME");

    double integral = 0;
    for(int iBin = 0;iBin < 20;iBin++) {
        integral += 1 - histDataRatios[1] -> GetBinContent(iBin+1);
    }
    std::cout << "Impact of the Pile-up: " << integral / 20 << std::endl;



    canvasData -> SaveAs(Form("%s/dataVsCentrality.pdf", outDir.c_str()));
    canvasMc -> SaveAs(Form("%s/mcVsCentrality.pdf", outDir.c_str()));
    canvasCompDataMc -> SaveAs(Form("%s/dataVsMcVsCentrality.pdf", outDir.c_str()));
    canvasDataEvSelRatioSummary -> SaveAs(Form("%s/evSelRatioSummary.pdf", outDir.c_str()));



}