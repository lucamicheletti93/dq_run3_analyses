double FuncAlign(double *, double *);

void alignment_fit() {
    TFile *fIn = new TFile("/Users/lucamicheletti/Downloads/Allineamento_0_60_0_25.root", "READ");
    TGraphErrors *graMassVsDeltaP = (TGraphErrors*) fIn -> Get("Graph");

    TF1 *funcAlign = new TF1("funcAlign", FuncAlign, -90., 90., 2);
    funcAlign -> SetParameter(0, 9.39);
    funcAlign -> SetParameter(1, 0.01);
    graMassVsDeltaP -> Fit(funcAlign, "R");

    TH2D *histGrid = new TH2D("histGrid", " ; #Delta p (GeV/c) ; #it{m} (GeV/c^{2}) ", 100, -90., 90., 100, 9.25, 9.6);
    TCanvas *canvasFit = new TCanvas("canvasFit","", 800, 600);
    histGrid -> Draw();
    graMassVsDeltaP -> Draw("EP SAME");
    
}
///////////////////////////////////////////////////////
double FuncAlign(double *x, double *par) {
    double m0 = par[0];
    double delta = par[1];
    return m0 / (1 + (delta * x[0]));
}