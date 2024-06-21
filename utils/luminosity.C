void luminosity() {
    const int nRuns = 13;
    int runList[] = {528531, 528461, 528292, 527899, 527895, 527871, 527850, 527240, 527109, 527057, 
                     527041, 526964, 526641};
    string pathToFiles = "/Users/lucamicheletti/cernbox/JPSI/Run3/2022/LHC22o_pass6_minBias/dq_streaming/merged_files";

    const double csTVX = 59400; // picobarn
    double selEff = 0;
    double nBcTVX = 0;

    int index = 0;
    double counterColAll = 0;
    double counterColAcc = 0;
    double counterTVX = 0;
    double counterLumiTVX = 0;
    
    double counterLumiDQ = 0;

    TH1F *histColCounterAll = new TH1F("histColCounterAll", "", nRuns, 0, nRuns);
    TH1F *histColCounterAcc = new TH1F("histColCounterAcc", "", nRuns, 0, nRuns);
    TH1F *histCounterTVX = new TH1F("histCounterTVX", "", nRuns, 0, nRuns);
    TH1F *histCounterLumiTVX = new TH1F("histCounterLumiTVX", "", nRuns, 0, nRuns);
    histCounterLumiTVX -> SetLineColor(kBlack);
    TH1F *histCounterLumiDQ = new TH1F("histCounterLumiDQ", "", nRuns, 0, nRuns);
    histCounterLumiDQ -> SetLineColor(kBlue);

    for (auto const& run : runList) {
        histColCounterAll -> GetXaxis() -> SetBinLabel(index+1, Form("%i", run));
        histColCounterAcc -> GetXaxis() -> SetBinLabel(index+1, Form("%i", run));
        histCounterTVX -> GetXaxis() -> SetBinLabel(index+1, Form("%i", run));
        histCounterLumiTVX -> GetXaxis() -> SetBinLabel(index+1, Form("%i", run));
        histCounterLumiDQ -> GetXaxis() -> SetBinLabel(index+1, Form("%i", run));

        TFile *fIn = new TFile(Form("%s/%i/AnalysisResults.root", pathToFiles.c_str(), run));
        TH1F *histTmpColCounterAll = (TH1F*) fIn -> Get("event-selection-task/hColCounterAll");
        TH1F *histTmpColCounterAcc = (TH1F*) fIn -> Get("event-selection-task/hColCounterAcc");
        TH1F *histTmpCounterTVX = (TH1F*) fIn -> Get("bc-selection-task/hCounterTVX");
        TH1F *histTmpCounterLumiTVX = (TH1F*) fIn -> Get("bc-selection-task/hLumiTVX");

        histColCounterAll -> SetBinContent(index+1, histTmpColCounterAll -> GetBinContent(1));
        histColCounterAcc -> SetBinContent(index+1, histTmpColCounterAcc -> GetBinContent(1));
        histCounterTVX -> SetBinContent(index+1, histTmpCounterTVX -> GetBinContent(1));
        histCounterLumiTVX -> SetBinContent(index+1, histTmpCounterLumiTVX -> GetBinContent(1));

        selEff = histTmpColCounterAcc -> GetBinContent(1) / histTmpColCounterAll -> GetBinContent(1);
        nBcTVX = histTmpCounterTVX -> GetBinContent(1);
        histCounterLumiDQ-> SetBinContent(index+1, (nBcTVX / csTVX) * selEff);


        counterColAll += histTmpColCounterAll -> GetBinContent(1);
        counterColAcc += histTmpColCounterAcc -> GetBinContent(1);
        counterTVX += histTmpCounterTVX -> GetBinContent(1);
        counterLumiTVX += histCounterLumiTVX -> GetBinContent(1);
        counterLumiDQ += (nBcTVX / csTVX) * selEff;

        std::cout << run << " -> accepted collisions = " << histTmpColCounterAcc -> GetBinContent(1) << " , LUMI TVX = " << histTmpCounterLumiTVX -> GetBinContent(1) << std::endl;
        index++;
    }

    TLatex *latexTitle = new TLatex();
    latexTitle -> SetTextSize(0.065);
    latexTitle -> SetNDC();
    latexTitle -> SetTextFont(42);

    TCanvas *canvasCounterTVX = new TCanvas("canvasCounterTVX", "", 1000, 600);
    gStyle -> SetOptStat(0);
    gPad -> SetLogy();
    histCounterTVX -> Draw("H");
    latexTitle -> DrawLatex(0.20, 0.80, Form("counter TVX = %1.0f", counterTVX));

    TCanvas *canvasCounterLumi = new TCanvas("canvasCounterLumi", "", 1000, 600);
    gStyle -> SetOptStat(0);
    gPad -> SetLogy();
    histCounterLumiTVX -> Draw("H");
    histCounterLumiDQ -> Draw("H SAME");
    latexTitle -> DrawLatex(0.20, 0.80, Form("Luminosity TVX = %3.2f pb^{-1}", counterLumiTVX * 1e-6));
    latexTitle -> DrawLatex(0.20, 0.72, Form("#color[4]{Luminosity DQ = %3.2f pb^{-1}}", counterLumiDQ * 1e-6));
}