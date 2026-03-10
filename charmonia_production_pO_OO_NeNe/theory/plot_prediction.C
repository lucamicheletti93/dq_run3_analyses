void LoadStyle();
void SetLegend(TLegend *);

void plot_prediction(TString model = "THU") {
    LoadStyle();
    double xMin, xMax, yCentral, yMin, yMax;

    std::ifstream fIn(Form("%s/forward/data_Raa_jpsi_0100_y254.dat",model.Data()));

    std::vector<double> xCentrals, exMins, exMaxs, exZeros, yCentrals, eyMins, eyMaxs, eyZeros;
    while (fIn >> xMin >> xMax >> yCentral >> yMin >> yMax) {
        double xCentral = (xMax + xMin) / 2.;

        xCentrals.push_back(xCentral);
        exMins.push_back(xCentral - xMin);
        exMaxs.push_back(xMax - xCentral);
        exZeros.push_back(0);
        yCentrals.push_back(yCentral);
        eyMins.push_back(yCentral - yMin);
        eyMaxs.push_back(yMax - yCentral);
        eyZeros.push_back(0);
    }

    TGraphAsymmErrors *graTheorErrBand = new TGraphAsymmErrors(xCentrals.size(), &(xCentrals[0]), &(yCentrals[0]), &(exMins[0]), &(exMaxs[0]), &(eyMins[0]), &(eyMaxs[0]));
    graTheorErrBand -> SetFillStyle(1001);
    graTheorErrBand -> SetFillColorAlpha(kAzure+4, 0.3);

    TGraphAsymmErrors *graTheorCentrVal = new TGraphAsymmErrors(xCentrals.size(), &(xCentrals[0]), &(yCentrals[0]), &(exZeros[0]), &(exZeros[0]), &(exZeros[0]), &(eyZeros[0]));
    graTheorCentrVal -> SetLineColor(kAzure+4);
    graTheorCentrVal -> SetLineWidth(2);

    TCanvas *canvasTheor = new TCanvas("canvasTheor", "", 800, 600);
    TH2D *histGrid = new TH2D("histGrid", ";#it{p}_{T} (GeV/#it{c});#it{R}_{AA}", 100, 0, 16, 100, 0, 1.2);
    histGrid -> Draw();
    graTheorErrBand -> Draw("E3 SAME");
    graTheorCentrVal -> Draw("L SAME");


    TFile *fOut = new TFile(Form("predictions_%s.root", model.Data()), "RECREATE");
    canvasTheor -> Write("canvasTheor");
    graTheorErrBand -> Write("graTheorErrBand");
    graTheorCentrVal -> Write("graTheorCentrVal");
    fOut -> Close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LoadStyle(){
    int font = 42;
    gStyle -> SetFrameBorderMode(0);
    gStyle -> SetFrameFillColor(0);
    gStyle -> SetCanvasBorderMode(0);
    gStyle -> SetPadBorderMode(0);
    gStyle -> SetPadColor(10);
    gStyle -> SetCanvasColor(10);
    gStyle -> SetTitleFillColor(10);
    gStyle -> SetTitleBorderSize(1);
    gStyle -> SetStatColor(10);
    gStyle -> SetStatBorderSize(1);
    gStyle -> SetLegendBorderSize(1);
    gStyle -> SetDrawBorder(0);
    gStyle -> SetTextFont(font);
    gStyle -> SetStatFontSize(0.05);
    gStyle -> SetStatX(0.97);
    gStyle -> SetStatY(0.98);
    gStyle -> SetStatH(0.03);
    gStyle -> SetStatW(0.3);
    gStyle -> SetTickLength(0.02,"y");
    gStyle -> SetEndErrorSize(3);
    gStyle -> SetLabelSize(0.05,"xyz");
    gStyle -> SetLabelFont(font,"xyz");
    gStyle -> SetLabelOffset(0.01,"xyz");
    gStyle -> SetTitleFont(font,"xyz");
    gStyle -> SetTitleOffset(0.9,"x");
    gStyle -> SetTitleOffset(1.02,"y");
    gStyle -> SetTitleSize(0.05,"xyz");
    gStyle -> SetMarkerSize(1.3);
    gStyle -> SetOptStat(0);
    gStyle -> SetEndErrorSize(0);
    gStyle -> SetCanvasPreferGL(kTRUE);
    gStyle -> SetHatchesSpacing(0.5);
    gStyle -> SetPadLeftMargin(0.15);
    gStyle -> SetPadBottomMargin(0.15);
    gStyle -> SetPadTopMargin(0.05);
    gStyle -> SetPadRightMargin(0.05);
    gStyle -> SetEndErrorSize(0.0);
    gStyle -> SetTitleSize(0.05,"X");
    gStyle -> SetTitleSize(0.045,"Y");
    gStyle -> SetLabelSize(0.045,"X");
    gStyle -> SetLabelSize(0.045,"Y");
    gStyle -> SetTitleOffset(1.2,"X");
    gStyle -> SetTitleOffset(1.35,"Y");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.04);
}