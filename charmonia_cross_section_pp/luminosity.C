void luminosity() {
    //string fInName = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_pp_cross_section_run3/data/2024/LHC24af_pass1_skimmed/AnalysisResults_time_assoc_fDiMuon.root";
    //string fInName = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_pp_cross_section_run3/data/2024/LHC24af_pass1_skimmed/AnalysisResults_time_assoc_fDiMuon_fromHF.root";
    string fInName = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_pp_cross_section_run3/data/2024/LHC24_pass1_min_bias/AnalysisResults_time_assoc.root";
    string triggerMask = "fDiMuon";
    bool isSkimmed = false;

    auto fIn = TFile::Open(fInName.c_str());
    if (!fIn || fIn -> IsZombie()) {
        return;
    }

    if (isSkimmed) {
        TH1F *histCounterTVX = (TH1F*) fIn -> Get("bc-selection-task/hCounterTVX");
        TH1F *histCounterTVXafterBCcuts = (TH1F*) fIn -> Get("bc-selection-task/hCounterTVXafterBCcuts");

        double counterTVX = histCounterTVX -> GetEntries();
        double counterTVXafterBCcuts = histCounterTVXafterBCcuts -> GetEntries();
        double effBCcuts = counterTVXafterBCcuts / counterTVX;
        double cunterVtxBeforeCuts;
        double cunterVtxAfterCuts;

        std::cout << "BC cuts efficiency = " << effBCcuts << std::endl;

        for (auto const& dirKey : *fIn -> GetListOfKeys()) {
            if (TString(dirKey -> GetName()).Contains("table-maker")) {
                TList *list1 = (TList*) fIn -> Get("table-maker/Statistics");
                TH2D *histZorroInfo = (TH2D*) list1 -> FindObject("ZorroInfo");
                TH2D *histZorroSel = (TH2D*) list1 -> FindObject("ZorroSel");

                int nRuns = histZorroInfo -> GetXaxis() -> GetNbins();
                double inspectedTVX = 0;
                double scalers = 0;
                double selectionsZorroInfo = 0;
                double selectionsZorroSel = 0;

                for (int iRun = 0;iRun < nRuns;iRun++) {
                    inspectedTVX += histZorroInfo -> GetBinContent(iRun+1, histZorroInfo -> GetYaxis() -> FindBin("inspected TVX"));
                    scalers += histZorroInfo -> GetBinContent(iRun+1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str())));
                    selectionsZorroInfo += histZorroInfo -> GetBinContent(iRun+1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str())));
                    selectionsZorroSel += histZorroSel -> GetBinContent(iRun+1, histZorroSel -> GetYaxis() -> FindBin(Form("%s", triggerMask.c_str())));

                    //std::cout << histZorroSel -> GetBinContent(iRun+1, histZorroSel -> GetYaxis() -> FindBin(Form("%s", triggerMask.c_str()))) << std::endl;
                }
                //std::cout << inspectedTVX << " " << scalers << " " << selectionsZorroInfo << " " << selectionsZorroSel << std::endl;

                TList *list2 = (TList*) fIn -> Get("table-maker/output");

                TList *listBeforeCuts = (TList*) list2 -> FindObject("Event_BeforeCuts");
                TH1F *histVtxBeforeCuts = (TH1F*) listBeforeCuts -> FindObject("VtxZ");

                TList *listAfterCuts = (TList*) list2 -> FindObject("Event_AfterCuts");
                TH1F *histVtxAfterCuts = (TH1F*) listAfterCuts -> FindObject("VtxZ");

                cunterVtxBeforeCuts = histVtxBeforeCuts -> GetEntries();
                cunterVtxAfterCuts = histVtxAfterCuts -> GetEntries();

                double luminosityZorroInfo = (inspectedTVX * (selectionsZorroInfo / scalers)) / 59.4e6;
                double luminosityZorroSel = (inspectedTVX * (selectionsZorroSel / scalers)) / 59.4e6;
                double luminosityZorroInfoCollEff = (inspectedTVX * (selectionsZorroInfo / scalers) * (cunterVtxAfterCuts / cunterVtxBeforeCuts)) / 59.4e6;
                double luminosityZorroSelCollEff = (inspectedTVX * (selectionsZorroSel / scalers) * (cunterVtxAfterCuts / cunterVtxBeforeCuts)) / 59.4e6;
                std::cout << "Number of MB events [ZorroInfo] = " << (inspectedTVX * (selectionsZorroInfo / scalers)) << std::endl;
                std::cout << "Number of MB events [ZorroSel] = " << (inspectedTVX * (selectionsZorroSel / scalers)) << std::endl;
                std::cout << "Number of MB events [ZorroInfo] * collision efficiency = " << (inspectedTVX * (selectionsZorroInfo / scalers)) * (cunterVtxAfterCuts / cunterVtxBeforeCuts) << std::endl;
                std::cout << "Number of MB events [ZorroSel] * collision efficiency = " << (inspectedTVX * (selectionsZorroSel / scalers)) * (cunterVtxAfterCuts / cunterVtxBeforeCuts) << std::endl;
                std::cout << "Integrated Luminosity [ZorroInfo] = " << luminosityZorroInfo << " nb-1" << std::endl;
                std::cout << "Integrated Luminosity [ZorroSel]  = " << luminosityZorroSel << " nb-1" << std::endl;
                std::cout << "Integrated Luminosity [ZorroInfo] * collision efficiency = " << luminosityZorroInfoCollEff << " nb-1" << std::endl;
                std::cout << "Integrated Luminosity [ZorroSel] * collision efficiency  = " << luminosityZorroSelCollEff << " nb-1" << std::endl;


                //std::cout << "inspected TVX: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin("inspected TVX")) << std::endl;
                //std::cout << triggerMask << " scalers: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str()))) << std::endl;
                //std::cout << triggerMask << " fSingleMuLow selections: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str()))) << std::endl;
                //std::cout << triggerMask << " fSingleMuLow scalers (SEL): " << histZorroSel -> GetBinContent(1, histZorroSel -> GetYaxis() -> FindBin(Form("%s", triggerMask.c_str()))) << std::endl;

                //counters[0] = histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin("inspected TVX"));
                //counters[1] = histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str())));
                //counters[2] = histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str())));
            }
        }
    } else {
        TH1F *histCounterTVX = (TH1F*) fIn -> Get("bc-selection-task/hCounterTVX");
        TH1F *histCounterTVXafterBCcuts = (TH1F*) fIn -> Get("bc-selection-task/hCounterTVXafterBCcuts");
        double counterTVX = histCounterTVX -> GetEntries();
        double counterTVXafterBCcuts = histCounterTVXafterBCcuts -> GetEntries();
        double cunterVtxBeforeCuts;
        double cunterVtxAfterCuts;

        for (auto const& dirKey : *fIn -> GetListOfKeys()) {
            if (TString(dirKey -> GetName()).Contains("table-maker")) {
                TList *list = (TList*) fIn -> Get("table-maker/output");

                TList *listBeforeCuts = (TList*) list -> FindObject("Event_BeforeCuts");
                TH1F *histVtxBeforeCuts = (TH1F*) listBeforeCuts -> FindObject("VtxZ");

                TList *listAfterCuts = (TList*) list -> FindObject("Event_AfterCuts");
                TH1F *histVtxAfterCuts = (TH1F*) listAfterCuts -> FindObject("VtxZ");

                cunterVtxBeforeCuts = histVtxBeforeCuts -> GetEntries();
                cunterVtxAfterCuts = histVtxAfterCuts -> GetEntries();
            }
        }

        double luminosityTVX = (counterTVX * (cunterVtxAfterCuts / cunterVtxBeforeCuts)) / 59.4e6;
        double luminosityTVXafterBCcuts = (counterTVXafterBCcuts * (cunterVtxAfterCuts / cunterVtxBeforeCuts)) / 59.4e6;
        std::cout << "Number of MB events [TVX] = " << (counterTVX * (cunterVtxAfterCuts / cunterVtxBeforeCuts)) << std::endl;
        std::cout << "Number of MB events [TVXafterBCcuts] = " << (counterTVXafterBCcuts * (cunterVtxAfterCuts / cunterVtxBeforeCuts)) << std::endl;
        std::cout << "Integrated Luminosity [TVX] = " << luminosityTVX << " nb-1" << std::endl;
        std::cout << "Integrated Luminosity [TVXafterBCcuts]  = " << luminosityTVXafterBCcuts << " nb-1" << std::endl;
    }
}