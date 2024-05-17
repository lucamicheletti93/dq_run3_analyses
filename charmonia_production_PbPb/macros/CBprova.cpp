#include "TLorentzVector.h"
//#include "MyEvent.h"
#include <stdio.h>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "THStack.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include "TNtuple.h"
#include "TH1F.h"
#include "TH3.h"
#include "TH2.h"
#include <TF1.h>
#include "TGraphErrors.h"
#include <TMultiGraph.h>
#include <TLegend.h>
#include "TStyle.h"
#include <vector>
#include <TGraph.h>
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TBrowser.h"
#include "TOrdCollection.h"
#include "TList.h"
#include <TString.h>
#include <algorithm>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include "TPaveText.h"
#include <numeric>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TRandom.h"
#include "RooGenericPdf.h"
#include "RooTFnBinding.h"
#include "RooCrystalBall.h"
#include "RooFitResult.h"
using namespace RooFit;



Double_t DoubleSidedCB2(double x, double mu, double width, double a1, double p1, double a2, double p2)
{
  double u   = (x-mu)/width;
  double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(1);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}


double DoubleSidedCB(double* x, double *par)
{
  return(par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]));
}

int DCB(){

TH1F *hmass = new TH1F("hmass","hmass",15,0.,15.);
TFile *f = TFile::Open("/Users/saragaretti/macros/analysis_J_psi/dq_fitter-main/centralityCuts.root", "UPDATE");
TH1D* hist = (TH1D*) f -> Get("Mass_Sub_0_20_EvMix");


TF1 *dcb = new TF1("dcb","DoubleSidedCB",2.9,3.096,0.96e-01);
double par[7];
par[0]=0.3;
par[1]=3.096;
//par[2]=0.96e-01;
par[2]=0.96e-01;
par[3]=1.061;
//par[3]=1.061;
par[4]=23.725;
//par[4]=23.725;
par[5]=3.183;
//par[5]=3.183;
par[6]=21.144;
dcb->SetParameters(&par[0]);
TCanvas *C = new TCanvas("C","C",700,700);
hist->Draw();
gStyle->SetOptStat(0);
hist->SetTitle(";m(KJ/#psi) [GeV];");
hist->Fit("dcb");

/*TLegend *LEG = new TLegend(0.6,0.8,0.89,0.89);
LEG->AddEntry(dcb,"Double Sided Crystal Ball");
LEG->SetBorderSize(0);
LEG->Draw("sames");*/


//Roofit
RooRealVar bMass("m(KJ/#psi)","m(J/#psi) [GeV]",2.6,3.6);
RooDataHist histo("histo","histo",bMass,Import(*hist));

//DCB parameters
RooRealVar mu("mu","mu",par[1],3.096,3.3);
RooRealVar width("width","width",par[2],0.66e-01,0.01);
RooRealVar a1("a1","a1",par[3],0.861,20.);
RooRealVar p1("p1","p1",par[4],2.525,27);
RooRealVar a2("a2","a2",par[5],0.183,20.);
RooRealVar p2("p2","p2",par[6],2.744,27);

RooCrystalBall dcbPdf("dcbPdf","DoubleSidedCB",bMass,mu,width,a1,p1,a2,p2);


auto result = dcbPdf.fitTo(histo, RooFit::Save(true));
//result->Print();
cout << result << endl;

auto pl = bMass.frame();

histo.plotOn(pl);
dcbPdf.plotOn(pl);
pl->Draw(); 

return 0;
}
 
/*void CBprova()
{
   
   
   CB2Pdf::CB2Pdf(const char *name, const char *title,
                        RooAbsReal& _x,
                        RooAbsReal& _A,
                        RooAbsReal& _B,
                        RooAbsReal& _C,
                        RooAbsReal& _D,
                        RooAbsReal& _E,
                        RooAbsReal& _F) :
   RooAbsPdf(name,title),
   x("x","x",this,_x),
   A("A","A",this,_A),
   B("B","B",this,_B),
   C("C","C",this,_C),
   D("D","D",this,_D),
   E("E","E",this,_E),
   F("F","F",this,_F)
 {
 }
   // ----------------------------------------------------
   // G e n e r i c   i n t e r p r e t e d   p . d . f .
   // ====================================================
 
   // Declare observable x
   RooRealVar x("x", "x", -20, 20);
   RooRealVar A("A", "A", 2., 3.096);
   RooRealVar B("B", "B", 0., 0.96e-01);
   RooRealVar C("C", "C", 0, 1.061);
   RooRealVar D("D","D",20.,23.725);
   RooRealVar E("E","E",2.,3.183);
   RooRealVar F("F","F",20.,21.144);
 
   // C o n s t r u c t   g e n e r i c   p d f   f r o m   i n t e r p r e t e d   e x p r e s s i o n
   // -------------------------------------------------------------------------------------------------
 
   // To construct a proper pdf, the formula expression is explicitly normalized internally by dividing
   // it by a numeric integral of the expression over x in the range [-20,20]
   //
   //RooRealVar alpha("alpha", "alpha", 5, 0.1, 10);
   RooGenericPdf genpdf("genpdf", "genpdf", "(1+0.1*abs(x)+sin(sqrt(abs(x*alpha+0.1))))", RooArgSet(x, A,B,C,D,E,F));
 
   // S a m p l e ,   f i t   a n d   p l o t   g e n e r i c   p d f
   // ---------------------------------------------------------------
 
   // Generate a toy dataset from the interpreted pdf
   std::unique_ptr<RooDataSet> data{genpdf.generate(x, 10000)};
 
   // Fit the interpreted pdf to the generated data
   genpdf.fitTo(*data, PrintLevel(-1));
 
   // Make a plot of the data and the pdf overlaid
   RooPlot *xframe = x.frame(Title("Interpreted expression pdf"));
   data->plotOn(xframe);
   genpdf.plotOn(xframe);
 
   // -----------------------------------------------------------------------------------------------------------
   // S t a n d a r d   p . d . f   a d j u s t   w i t h   i n t e r p r e t e d   h e l p e r   f u n c t i o n
   // ==========================================================================================================
   // Make a gauss(x,sqrt(mean2),sigma) from a standard RooGaussian
 
   // C o n s t r u c t   s t a n d a r d   p d f  w i t h   f o r m u l a   r e p l a c i n g   p a r a m e t e r
   // ------------------------------------------------------------------------------------------------------------
 
   // Construct parameter mean2 and sigma
   RooRealVar mean2("mean2", "mean^2", 10, 0, 200);
   RooRealVar sigma("sigma", "sigma", 3, 0.1, 10);
 
   // Construct interpreted function mean = sqrt(mean^2)
   RooFormulaVar mean("mean", "mean", "sqrt(mean2)", mean2);
 
   // Construct a gaussian g2(x,sqrt(mean2),sigma) ;
   RooGaussian g2("g2", "h2", x, mean, sigma);
 
   // G e n e r a t e   t o y   d a t a
   // ---------------------------------
 
   // Construct a separate gaussian g1(x,10,3) to generate a toy Gaussian dataset with mean 10 and width 3
   RooGaussian g1("g1", "g1", x, 10.0, 3.0);
   std::unique_ptr<RooDataSet> data2{g1.generate(x, 1000)};
 
   // F i t   a n d   p l o t   t a i l o r e d   s t a n d a r d   p d f
   // -------------------------------------------------------------------
 
   // Fit g2 to data from g1
   std::unique_ptr<RooFitResult> fitResult{g2.fitTo(*data2, Save(), PrintLevel(-1))};
   fitResult->Print();
 
   // Plot data on frame and overlay projection of g2
   RooPlot *xframe2 = x.frame(Title("Tailored Gaussian pdf"));
   data2->plotOn(xframe2);
   g2.plotOn(xframe2);
 
   // Draw all frames on a canvas
   TCanvas *c = new TCanvas("rf103_interprfuncs", "rf103_interprfuncs", 800, 400);
   c->Divide(2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   xframe->GetYaxis()->SetTitleOffset(1.4);
   xframe->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   xframe2->GetYaxis()->SetTitleOffset(1.4);
   xframe2->Draw();
   

   // S e t u p   m o d e l
   // ---------------------
 
   // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
   RooRealVar x("x", "x", -10, 10);
   RooRealVar A("A", "A of gaussian", 1, -10, 10);
   RooRealVar B("B", "B of gaussian", 1, 0.1, 10);
   RooRealVar C("C", "C of gaussian", 1, 0.1, 10)
    RooRealVar D("D", "D of gaussian", 1, -10, 10);
   RooRealVar E("E", "E of gaussian", 1, 0.1, 10);
    RooRealVar F("F", "F of gaussian", 1, -10, 10);
   
 
   // Build gaussian pdf in terms of x,mean and sigma
   RooGaussian gauss("gauss", "gaussian PDF", x, A, B, C, D, E, F);
 
   // Construct plot frame in 'x'
   RooPlot *xframe = x.frame(Title("Gaussian pdf."));
 
   // P l o t   m o d e l   a n d   c h a n g e   p a r a m e t e r   v a l u e s
   // ---------------------------------------------------------------------------
 
   // Plot gauss in frame (i.e. in x)
   gauss.plotOn(xframe);
 
   // Change the value of sigma to 3
   A.setVal(3.096);
   B.setVal(0.96e-01);
   C.setVal(1.061);
   D.setVal(23.725);
   E.setVal(3.183);
   F.setVal(21.144);
 
   // Plot gauss in frame (i.e. in x) and draw frame on canvas
   gauss.plotOn(xframe, LineColor(kRed));
 
   // G e n e r a t e   e v e n t s
   // -----------------------------
 
   // Generate a dataset of 1000 events in x from gauss
   std::unique_ptr<RooDataSet> data{gauss.generate(x, 10000)};
 
   // Make a second plot frame in x and draw both the
   // data and the pdf in the frame
   RooPlot *xframe2 = x.frame(Title("Gaussian pdf with data"));
   data->plotOn(xframe2);
   gauss.plotOn(xframe2);
 
   // F i t   m o d e l   t o   d a t a
   // -----------------------------
 
   // Fit pdf to data
   gauss.fitTo(*data, PrintLevel(-1));
 
   // Print values of mean and sigma (that now reflect fitted values and errors)
   A.Print();
   B.Print();
   C.Print();
   D.Print();
   E.Print();
   F.Print();
 
   // Draw all frames on a canvas
   TCanvas *c = new TCanvas("rf101_basics", "rf101_basics", 800, 400); c->Divide(2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   xframe->GetYaxis()->SetTitleOffset(1.6);
   xframe->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   xframe2->GetYaxis()->SetTitleOffset(1.6);
   xframe2->Draw();
}*/