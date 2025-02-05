void v2_different_freeze_out_surface() {
    Double_t pt_upsilon[41] = {0.,  0.25, 0.5, 0.75, 1.,  1.25, 1.5, 1.75, 2.,  2.25, 2.5, 2.75, 3.,  3.25,
                               3.5, 3.75, 4.,  4.25, 4.5, 4.75, 5.,  5.25, 5.5, 5.75, 6.,  6.25, 6.5, 6.75,
                               7.,  7.25, 7.5, 7.75, 8.,  8.25, 8.5, 8.75, 9.,  9.25, 9.5, 9.75, 10.};

    Double_t v2_upsilon_fos_A[41] = {
        0.,         0.0000903374, 0.000352131, 0.000759484, 0.00127459, 0.00185423, 0.00245631, 0.00304507, 0.00359415,
        0.00408731, 0.00451741,   0.00488433,  0.0051926,   0.00544918, 0.00566188, 0.00583824, 0.00598504, 0.00610807,
        0.0062123,  0.00630219,   0.00638231,  0.00645832,  0.00653848, 0.00663589, 0.00677175, 0.00697955, 0.00731017,
        0.00783738, 0.00866286,   0.00991978,  0.0117735,   0.0144189,  0.0180731,  0.0229649,  0.0293209,  0.0373505,
        0.0472315,  0.0590972,    0.0730276,   0.0890433,   0.107105};

    Double_t v2_upsilon_fos_B[41] = {0.,           0.0000392938, 0.000150571, 0.000315504, 0.000508022,  0.000699556,
                                     0.000864068,  0.000981709,  0.00104051,  0.00103626,  0.000971071,  0.00085134,
                                     0.000685734,  0.000483539,  0.000253547, 3.45811e-6,  -0.000260324, -0.0005326,
                                     -0.000809088, -0.00108595,  -0.00135905, -0.00162289, -0.00186888,  -0.00208285,
                                     -0.00224158,  -0.00230827,  -0.00222724, -0.00191835, -0.00127207,  -0.000146488,
                                     0.00163283,   0.00426895,   0.00798561,  0.0130154,   0.0195856,    0.0279028,
                                     0.0381387,    0.0504177,    0.0648091,   0.0813219,   0.0999042};

    Double_t v2_upsilon_fos_C[41] = {0.,         0.00014908, 0.000585317, 0.00127753, 0.0021794, 0.00323649, 0.00439367,
                                     0.00560144, 0.00682013, 0.00802161,  0.00918899, 0.0103147, 0.0113981,  0.0124427,
                                     0.0134546,  0.0144404,  0.0154067,   0.0163594,  0.0173038, 0.0182446,  0.0191862,
                                     0.020134,   0.0210945,  0.0220782,   0.0231013,  0.0241893, 0.0253817,  0.0267364,
                                     0.0283351,  0.0302868,  0.0327306,   0.0358341,  0.0397897, 0.0448059,  0.0510961,
                                     0.0588648,  0.0682941,  0.0795296,   0.092671,   0.107764,  0.124796};

    Double_t pt_proton[41] = {0.,  0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,  1.1, 1.2, 1.3,
                              1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.,  2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
                              2.8, 2.9, 3.,  3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.};

    Double_t v2_proton_fos_A[41] = {0.,       0.000593191, 0.0024318, 0.00565663, 0.010407,  0.0167434, 0.0246149,
                                    0.033872, 0.0443049,   0.0556832, 0.0677858,  0.0804169, 0.0934128, 0.106641,
                                    0.119998, 0.133403,    0.146794,  0.160125,   0.173361,  0.186475,  0.199449,
                                    0.212266, 0.224918,    0.237396,  0.249694,   0.26181,   0.27374,   0.285483,
                                    0.297038, 0.308405,    0.319584,  0.330576,   0.341382,  0.352002,  0.362438,
                                    0.37269,  0.382762,    0.392654,  0.402367,   0.411904,  0.421267};

    Double_t v2_proton_fos_B[41] = {0.,        0.000499601, 0.00209458, 0.00501368, 0.00948583, 0.0156269, 0.0234001,
                                    0.0326431, 0.0431217,   0.0545821,  0.0667844,  0.0795203,  0.0926182, 0.105942,
                                    0.119385,  0.132867,    0.146326,   0.159716,   0.173003,   0.186162,  0.199174,
                                    0.212025,  0.224706,    0.23721,    0.24953,    0.261665,   0.273612,  0.28537,
                                    0.296938,  0.308317,    0.319506,   0.330506,   0.34132,    0.351946,  0.362388,
                                    0.372647,  0.382723,    0.392618,   0.402336,   0.411876,   0.421241};

    Double_t v2_proton_fos_C[41] = {0.,        0.00073592, 0.00295663, 0.00668858, 0.0119421, 0.0186819, 0.0268127,
                                    0.0361852, 0.0466152,  0.0579062,  0.0698688,  0.0823323, 0.0951509, 0.108204,
                                    0.121394,  0.134644,   0.147895,   0.1611,     0.174224,  0.187238,  0.200123,
                                    0.212863,  0.225445,   0.237864,   0.250109,   0.262178,  0.274067,  0.285774,
                                    0.297297,  0.308636,   0.31979,    0.33076,    0.341546,  0.352148,  0.362568,
                                    0.372808,  0.382867,   0.392748,   0.402451,   0.41198,   0.421334};

    TGraphErrors* g_ALICE = new TGraphErrors(3);
    g_ALICE->SetPoint(0, 1.89, 0.013); g_ALICE->SetPointError(0, 0, 0.0440114);
    g_ALICE->SetPoint(1, 4.42, -0.01); g_ALICE->SetPointError(1, 0, 0.0417732);
    g_ALICE->SetPoint(2, 8.86, 0.003); g_ALICE->SetPointError(2, 0, 0.0593043);

    TGraph *g_v2_upsilon_fos_A = new TGraph(41, pt_upsilon, v2_upsilon_fos_A);
    TGraph *g_v2_upsilon_fos_B = new TGraph(41, pt_upsilon, v2_upsilon_fos_B);
    TGraph *g_v2_upsilon_fos_C = new TGraph(41, pt_upsilon, v2_upsilon_fos_C);

    TGraph *g_v2_proton_fos_A = new TGraph(41, pt_proton, v2_proton_fos_A);
    TGraph *g_v2_proton_fos_B = new TGraph(41, pt_proton, v2_proton_fos_B);
    TGraph *g_v2_proton_fos_C = new TGraph(41, pt_proton, v2_proton_fos_C);

    gStyle->SetOptStat(kFALSE);

    gStyle->SetLabelOffset(0.005, "x"); // 0.005 = root default
    gStyle->SetLabelOffset(0.005, "y"); // 0.005 = root default
	gStyle->SetLabelSize(.05, "XY");
    gStyle->SetTitleXSize(0.06);        // 0.04  = root default
    gStyle->SetTitleYSize(0.06);        // 0.04  = root default
    gStyle->SetTitleXOffset(1.0);
    gStyle->SetTitleYOffset(1.1);

    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.06); // 0.1 = root default
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetPadBottomMargin(0.15);

    gStyle->SetCanvasColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetOptTitle(0);

    gROOT->ForceStyle();

    TCanvas *c1 = new TCanvas("c1", "c1", 650, 500);

    TH2F *fr = new TH2F("fr", "fr", 1, 0., 10.5, 1, -0.02, 0.12);
    fr->SetXTitle("#it{p}_{T} (GeV/#it{c})");
    fr->SetYTitle("#it{v}_{2}");
    fr->Draw();

    TLine* zero = new TLine(0., 0., 10.5, 0.);
    zero->SetLineStyle(2);
    zero->Draw();

    g_v2_upsilon_fos_A->SetLineColor(kBlue);
    g_v2_upsilon_fos_A->SetLineWidth(3);
    g_v2_upsilon_fos_A->Draw("l");

    g_v2_upsilon_fos_B->SetLineColor(kRed);
    g_v2_upsilon_fos_B->SetLineWidth(3);
    g_v2_upsilon_fos_B->SetLineStyle(2);
    g_v2_upsilon_fos_B->Draw("l");

    g_v2_upsilon_fos_C->SetLineColor(kGreen+3);
    g_v2_upsilon_fos_C->SetLineWidth(3);
    g_v2_upsilon_fos_C->SetLineStyle(3);
    g_v2_upsilon_fos_C->Draw("l");

    g_ALICE->SetMarkerStyle(20);
    g_ALICE->Draw("p");

    auto legend_upsilon = new TLegend(0.35, 0.58, 0.8, 0.88);
    // legend_upsilon->SetHeader("#Upsilon ");
    legend_upsilon->AddEntry(g_ALICE, "#Upsilon, Pb-Pb, 5.02 TeV (ALICE)", "p");
    legend_upsilon->AddEntry(g_v2_upsilon_fos_A, "#tau_{f}(#hat{#it{r}}) = const. = #tau_{f0}");
    legend_upsilon->AddEntry(g_v2_upsilon_fos_B, "#tau_{f}(#hat{#it{r}}) = #tau_{f0} (1 + 0.5 #hat{#it{r}})");
    legend_upsilon->AddEntry(g_v2_upsilon_fos_C, "#tau_{f}(#hat{#it{r}}) = #tau_{f0} (1 - 0.5 #hat{#it{r}})");
    legend_upsilon->SetTextSize(0.04);
    legend_upsilon->SetBorderSize(0);
    legend_upsilon->Draw();

    c1->SaveAs("v2_upsilon_different_freeze_out_surfaces.pdf");

    // TCanvas *c2 = new TCanvas("c2");

    // TH2F *fr2 = new TH2F("fr2", "fr2", 1, 0., 10.5, 1, -0.02, 0.12);
    // fr2->SetXTitle("p_{T} GeV/c");
    // fr2->SetYTitle("v_{2}");
    // fr2->Draw();

    // g_v2_proton_fos_A->SetLineColor(kBlue);
    // g_v2_proton_fos_A->Draw("l");

    // g_v2_proton_fos_B->SetLineColor(kRed);
    // g_v2_proton_fos_B->Draw("l");

    // g_v2_proton_fos_C->SetLineColor(kGreen);
    // g_v2_proton_fos_C->Draw("l");


}