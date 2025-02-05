void luminosity() {
    string fInName = "/Users/lucamicheletti/cernbox/JPSI/Jpsi_pp_cross_section_run3/data/2024/LHC24aj_pass1_skimmed/AnalysisResults_time_assoc_fDiMuon.root";
    string triggerMask = "fDiMuon";

    auto fIn = TFile::Open(fInName.c_str());
    if (!fIn || fIn -> IsZombie()) {
        return;
    }

    // CORREZIONE PER IL BC!!!!
    for (auto const& dirKey : *fIn -> GetListOfKeys()) {
        if (TString(dirKey -> GetName()).Contains("table-maker")) {
            TList *list = (TList*) fIn -> Get("table-maker/Statistics");
            TH2D *histZorroInfo = (TH2D*) list -> FindObject("ZorroInfo");
            TH2D *histZorroSel = (TH2D*) list -> FindObject("ZorroSel");

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
            std::cout << inspectedTVX << " " << scalers << " " << selectionsZorroInfo << " " << selectionsZorroSel << std::endl;

            double luminosityZorroInfo = (inspectedTVX * (selectionsZorroInfo / scalers)) / 59.4e6;
            double luminosityZorroSel = (inspectedTVX * (selectionsZorroSel / scalers)) / 59.4e6;
            std::cout << "Integrated Luminosity [ZorroInfo] = " << luminosityZorroInfo << " nb-1" << std::endl;
            std::cout << "Integrated Luminosity [ZorroSel]  = " << luminosityZorroSel << " nb-1" << std::endl;


            //std::cout << "inspected TVX: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin("inspected TVX")) << std::endl;
            //std::cout << triggerMask << " scalers: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str()))) << std::endl;
            //std::cout << triggerMask << " fSingleMuLow selections: " << histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str()))) << std::endl;
            //std::cout << triggerMask << " fSingleMuLow scalers (SEL): " << histZorroSel -> GetBinContent(1, histZorroSel -> GetYaxis() -> FindBin(Form("%s", triggerMask.c_str()))) << std::endl;

            //counters[0] = histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin("inspected TVX"));
            //counters[1] = histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str())));
            //counters[2] = histZorroInfo -> GetBinContent(1, histZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str())));
        }
    }
}