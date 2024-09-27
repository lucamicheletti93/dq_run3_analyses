void LoadStyle();
void SetLegend(TLegend *);

void theory_predictions() {
    LoadStyle();

    // taken from https://docs.google.com/spreadsheets/d/1gkvbqAOPlNoyyDsgQ7QtirU1ZQ-LIYxzaXsjY1XG0uI/edit?gid=1699652183#gid=1699652183

    /////////////////////////
    // CGC + NRQCD + FONLL //
    /////////////////////////
    const int nPtBinsCgcNrqcdFonll = 6;
    double ptCentrCgcNrqcdFonll[] = {0.5, 1.5, 2.5, 3.5, 4.5, 6};
    double ptMinCgcNrqcdFonll[] = {0.5, 0.5, 0.5, 0.5, 0.5, 1};
    double ptMaxCgcNrqcdFonll[] = {0.5, 0.5, 0.5, 0.5, 0.5, 1};
    double ptBinWidthCgcNrqcdFonll[] = {1, 1, 1, 1, 1, 2};

    /////////////////////////
    // J/psi               //
    /////////////////////////
    // FONLL
    double csJpsiFonll1Pt[] = {1.45E+02, 3.03E+02, 2.80E+02, 1.99E+02, 1.28E+02, 1.32E+02};

    double csJpsiFonll1PtMinS[] = {6.25E+01, 1.23E+02, 1.04E+02, 6.70E+01, 3.88E+01, 3.44E+01};
    double csJpsiFonll1PtMinM[] = {2.61E+01, 5.06E+01, 4.13E+01, 2.52E+01, 1.40E+01, 1.20E+01};
    double csJpsiFonll1PtMinPdf[] = {2.26E+01, 3.99E+01, 2.76E+01, 1.37E+01, 6.40E+00, 4.71E+00};

    double csJpsiFonll1PtMaxS[] = {6.76E+01, 1.35E+02, 1.21E+02, 8.50E+01, 5.41E+01, 5.33E+01};
    double csJpsiFonll1PtMaxM[] = {3.32E+01, 6.29E+01, 4.96E+01, 2.96E+01, 1.61E+01, 1.34E+01};
    double csJpsiFonll1PtMaxPdf[] = {2.26E+01, 3.99E+01, 2.74E+01, 1.37E+01, 6.43E+00, 4.69E+00};

    // CGC + NRQCD
    double csJpsiCgcNrqcdPt[] = {1.01E+03, 2.47E+03, 2.45E+03, 1.71E+03, 1.02E+03, 9.06E+02};

    double csJpsiCgcNrqcdPtMin[] = {3.29E+02, 8.04E+02, 8.00E+02, 5.58E+02, 3.35E+02, 3.00E+02};
    double csJpsiCgcNrqcdPtMax[] = {3.2920E+02, 8.0424E+02, 7.9951E+02, 5.5786E+02, 3.3482E+02, 2.9998E+02};


    // CGC + NRQCD + FONLL
    double csJpsiCgcNrqcFonllPt[nPtBinsCgcNrqcdFonll], csJpsiCgcNrqcFonllPtMin[nPtBinsCgcNrqcdFonll], csJpsiCgcNrqcFonllPtMax[nPtBinsCgcNrqcdFonll];

    for (int iPt = 0;iPt < nPtBinsCgcNrqcdFonll;iPt++) {
        csJpsiCgcNrqcFonllPt[iPt] = (csJpsiFonll1Pt[iPt] + csJpsiCgcNrqcdPt[iPt]) / ptBinWidthCgcNrqcdFonll[iPt];
        csJpsiCgcNrqcFonllPtMin[iPt] = TMath::Sqrt(TMath::Power(csJpsiCgcNrqcdPtMin[iPt], 2.) + TMath::Power(csJpsiFonll1PtMinS[iPt], 2.) + TMath::Power(csJpsiFonll1PtMinM[iPt], 2.) + TMath::Power(csJpsiFonll1PtMinPdf[iPt], 2.)) / ptBinWidthCgcNrqcdFonll[iPt];
        csJpsiCgcNrqcFonllPtMax[iPt] = TMath::Sqrt(TMath::Power(csJpsiCgcNrqcdPtMax[iPt], 2.) + TMath::Power(csJpsiFonll1PtMaxS[iPt], 2.) + TMath::Power(csJpsiFonll1PtMaxM[iPt], 2.) + TMath::Power(csJpsiFonll1PtMaxPdf[iPt], 2.)) / ptBinWidthCgcNrqcdFonll[iPt];
    }

    /////////////////////////
    // Psi(2S)             //
    /////////////////////////
    double csPsi2sFonll1Pt[] = {3.09E+01, 6.85E+01, 6.98E+01, 5.46E+01, 3.81E+01, 4.23E+01};

    double csPsi2sFonll1PtMinS[] = {1.36E+01, 2.85E+01, 2.68E+01, 1.91E+01, 1.21E+01, 1.16E+01};
    double csPsi2sFonll1PtMinM[] = {5.70E+00, 1.17E+01, 1.07E+01, 7.29E+00, 4.41E+00, 4.09E+00};
    double csPsi2sFonll1PtMinPdf[] = {5.07E+00, 9.68E+00, 7.49E+00, 4.14E+00, 2.09E+00, 1.64E+00};

    double csPsi2sFonll1PtMaxS[] = {1.48E+01, 3.07E+01, 3.01E+01, 2.35E+01, 1.63E+01, 1.75E+01};
    double csPsi2sFonll1PtMaxM[] = {7.31E+00, 1.48E+01, 1.30E+01, 8.61E+00, 5.09E+00, 4.59E+00};
    double csPsi2sFonll1PtMaxPdf[] = {5.07E+00, 9.76E+00, 7.49E+00, 4.14E+00, 2.09E+00, 1.63E+00};

    // CGC + NRQCD
    double csPsi2sCgcNrqcdPt[] = {1.27E+02, 3.17E+02, 3.33E+02, 2.51E+02, 1.62E+02, 1.57E+02};

    double csPsi2sCgcNrqcdPtMin[] = {4.95E+01, 1.24E+02, 1.29E+02, 9.64E+01, 6.20E+01, 6.00E+01};
    double csPsi2sCgcNrqcdPtMax[] = {4.9537E+01, 1.2357E+02, 1.2889E+02, 9.6384E+01, 6.2025E+01, 6.0044E+01};


    // CGC + NRQCD + FONLL
    double csPsi2sCgcNrqcFonllPt[nPtBinsCgcNrqcdFonll], csPsi2sCgcNrqcFonllPtMin[nPtBinsCgcNrqcdFonll], csPsi2sCgcNrqcFonllPtMax[nPtBinsCgcNrqcdFonll];

    for (int iPt = 0;iPt < nPtBinsCgcNrqcdFonll;iPt++) {
        csPsi2sCgcNrqcFonllPt[iPt] = (csPsi2sFonll1Pt[iPt] + csPsi2sCgcNrqcdPt[iPt]) / ptBinWidthCgcNrqcdFonll[iPt];
        csPsi2sCgcNrqcFonllPtMin[iPt] = TMath::Sqrt(TMath::Power(csPsi2sCgcNrqcdPtMin[iPt], 2.) + TMath::Power(csPsi2sFonll1PtMinS[iPt], 2.) + TMath::Power(csPsi2sFonll1PtMinM[iPt], 2.) + TMath::Power(csPsi2sFonll1PtMinPdf[iPt], 2.)) / ptBinWidthCgcNrqcdFonll[iPt];
        csPsi2sCgcNrqcFonllPtMax[iPt] = TMath::Sqrt(TMath::Power(csPsi2sCgcNrqcdPtMax[iPt], 2.) + TMath::Power(csPsi2sFonll1PtMaxS[iPt], 2.) + TMath::Power(csPsi2sFonll1PtMaxM[iPt], 2.) + TMath::Power(csPsi2sFonll1PtMaxPdf[iPt], 2.)) / ptBinWidthCgcNrqcdFonll[iPt];
    }

    ///////////////////
    // NRQCD + FONLL //
    ///////////////////
    const int nPtBinsNrqcdFonll = 5;
    double ptCentrNrqcdFonll[] = {3.5, 4.5, 6, 8.5, 15};
    double ptMinNrqcdFonll[] = {0.5, 0.5, 1, 1.5, 5};
    double ptMaxNrqcdFonll[] = {0.5, 0.5, 1, 1.5, 5};
    double ptBinWidthNrqcdFonll[] = {1, 1, 2, 3, 10};

    /////////////////////////
    // J/psi               //
    /////////////////////////
    // FONLL
    double csJpsiFonll2Pt[] = {1.99E+02, 1.28E+02, 1.32E+02, 6.96E+01, 3.54E+01};

    double csJpsiFonll2PtMinS[] = {6.70E+01, 3.88E+01, 3.44E+01, 1.61E+01, 7.24E+00};
    double csJpsiFonll2PtMinM[] = {2.52E+01, 1.40E+01, 1.20E+01, 4.92E+00, 1.79E+00};
    double csJpsiFonll2PtMinPdf[] = {1.37E+01, 6.40E+00, 4.71E+00, 1.72E+00, 6.26E-01};

    double csJpsiFonll2PtMaxS[] = {8.50E+01, 5.41E+01, 5.33E+01, 2.54E+01, 1.06E+01};
    double csJpsiFonll2PtMaxM[] = {2.96E+01, 1.61E+01, 1.34E+01, 5.38E+00, 1.88E+00};
    double csJpsiFonll2PtMaxPdf[] = {1.37E+01, 6.43E+00, 4.69E+00, 1.72E+00, 6.26E-01};

    // NRQCD
    double csJpsiNrqcdPt[] = {1.98E+03, 8.69E+02, 5.86E+02, 2.05E+02, 6.96E+01};

    double csJpsiNrqcdPtMin[] = {3.04E+02, 2.04E+02, 1.54E+02, 7.00E+01, 2.79E+01};
    double csJpsiNrqcdPtMax[] = {3.04E+02, 2.06E+02, 2.49E+02, 1.24E+02, 5.95E+01};

    // NRQCD + FONLL
    double csJpsiNrqcFonllPt[nPtBinsNrqcdFonll], csJpsiNrqcFonllPtMin[nPtBinsNrqcdFonll], csJpsiNrqcFonllPtMax[nPtBinsNrqcdFonll];

    for (int iPt = 0;iPt < nPtBinsNrqcdFonll;iPt++) {
        csJpsiNrqcFonllPt[iPt] = (csJpsiFonll2Pt[iPt] + csJpsiNrqcdPt[iPt]) / ptBinWidthNrqcdFonll[iPt];
        csJpsiNrqcFonllPtMin[iPt] = TMath::Sqrt(TMath::Power(csJpsiNrqcdPtMin[iPt], 2.) + TMath::Power(csJpsiFonll2PtMinS[iPt], 2.) + TMath::Power(csJpsiFonll2PtMinM[iPt], 2.) + TMath::Power(csJpsiFonll2PtMinPdf[iPt], 2.)) / ptBinWidthCgcNrqcdFonll[iPt];
        csJpsiNrqcFonllPtMax[iPt] = TMath::Sqrt(TMath::Power(csJpsiNrqcdPtMax[iPt], 2.) + TMath::Power(csJpsiFonll2PtMaxS[iPt], 2.) + TMath::Power(csJpsiFonll2PtMaxM[iPt], 2.) + TMath::Power(csJpsiFonll2PtMaxPdf[iPt], 2.)) / ptBinWidthCgcNrqcdFonll[iPt];
    }

    /////////////////////////
    // Psi(2S)             //
    /////////////////////////
    // FONLL
    double csPsi2sFonll2Pt[] = {5.46E+01, 3.81E+01, 4.23E+01, 2.42E+01, 1.33E+01};

    double csPsi2sFonll2PtMinS[] = {1.91E+01, 1.21E+01, 1.16E+01, 5.74E+00, 2.80E+00};
    double csPsi2sFonll2PtMinM[] = {7.29E+00, 4.41E+00, 4.09E+00, 1.82E+00, 7.18E-01};
    double csPsi2sFonll2PtMinPdf[] = {4.14E+00, 2.09E+00, 1.64E+00, 6.40E-01, 2.46E-01};

    double csPsi2sFonll2PtMaxS[] = {2.35E+01, 1.63E+01, 1.75E+01, 9.14E+00, 4.15E+00};
    double csPsi2sFonll2PtMaxM[] = {8.61E+00, 5.09E+00, 4.59E+00, 1.99E+00, 7.55E-01};
    double csPsi2sFonll2PtMaxPdf[] = {4.14E+00, 2.09E+00, 1.63E+00, 6.38E-01, 2.46E-01};

    // NRQCD
    double csPsi2sNrqcdPt[] = {4.29E+02, 1.97E+02, 1.55E+02, 6.06E+01, 2.44E+01};

    double csPsi2sNrqcdPtMin[] = {5.18E+01, 2.28E+01, 1.68E+01, 6.88E+00, 2.86E+00};
    double csPsi2sNrqcdPtMax[] = {2.73E+01, 3.89E+01, 2.62E+01, 1.12E+01, 2.34E+00};

    // NRQCD + FONLL
    double csPsi2sNrqcFonllPt[nPtBinsNrqcdFonll], csPsi2sNrqcFonllPtMin[nPtBinsNrqcdFonll], csPsi2sNrqcFonllPtMax[nPtBinsNrqcdFonll];

    for (int iPt = 0;iPt < nPtBinsNrqcdFonll;iPt++) {
        csPsi2sNrqcFonllPt[iPt] = (csPsi2sFonll2Pt[iPt] + csPsi2sNrqcdPt[iPt]) / ptBinWidthNrqcdFonll[iPt];
        csPsi2sNrqcFonllPtMin[iPt] = TMath::Sqrt(TMath::Power(csPsi2sNrqcdPtMin[iPt], 2.) + TMath::Power(csPsi2sFonll2PtMinS[iPt], 2.) + TMath::Power(csPsi2sFonll2PtMinM[iPt], 2.) + TMath::Power(csPsi2sFonll2PtMinPdf[iPt], 2.)) / ptBinWidthCgcNrqcdFonll[iPt];
        csPsi2sNrqcFonllPtMax[iPt] = TMath::Sqrt(TMath::Power(csPsi2sNrqcdPtMax[iPt], 2.) + TMath::Power(csPsi2sFonll2PtMaxS[iPt], 2.) + TMath::Power(csPsi2sFonll2PtMaxM[iPt], 2.) + TMath::Power(csPsi2sFonll2PtMaxPdf[iPt], 2.)) / ptBinWidthCgcNrqcdFonll[iPt];
    }

    ///////////////////
    // ICEM + FONLL //
    //////////////////
    /*int nPtBinsIcemFonll = 0;
    vector <double> ptCentrIcemFonll;
    vector <double> ptMinIcemFonll;
    vector <double> ptMaxIcemFonll;
    vector <double> ptBinWidthIcemFonll;

    vector <double> csJpsiIcemPt;
    vector <double> csJpsiIcemPtMinRS;
    vector <double> csJpsiIcemPtMinMS;
    vector <double> csJpsiIcemPtMaxRS;
    vector <double> csJpsiIcemPtMaxMS;

    double ptCentrTmp, ptWidthTmp, csTmp, csMinRsTmp, csMaxRsTmp, csMinMsTmp, csMaxMsTmp;
    ifstream fJpsiIcem("jpsi_icem.txt");
    while(!fJpsiIcem.eof()) {
        fJpsiIcem >> ptCentrTmp >> ptWidthTmp >> csTmp >> csMinRsTmp >> csMaxRsTmp >> csMinMsTmp >> csMaxMsTmp;
        ptCentrIcemFonll.push_back(ptCentrTmp - 0.25);
        ptMinIcemFonll.push_back(0.25);
        ptMaxIcemFonll.push_back(0.25);
        ptBinWidthIcemFonll.push_back(0.5);
        csJpsiIcemPt.push_back(csTmp);
        csJpsiIcemPtMinRS.push_back(csTmp - csMinRsTmp);
        csJpsiIcemPtMaxRS.push_back(csMaxRsTmp - csTmp);
        csJpsiIcemPtMinMS.push_back(csTmp - csMinMsTmp);
        csJpsiIcemPtMaxMS.push_back(csMaxMsTmp - csTmp);
        nPtBinsIcemFonll++;
    }*/

    const int nPtBinsIcemFonll = 8;
    double ptCentrIcemFonll[] = {0.5, 1.5, 2.5, 3.5, 4.5, 6, 8.5, 15};
    double ptMinIcemFonll[] = {0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5, 5};
    double ptMaxIcemFonll[] = {0.5, 0.5, 0.5, 0.5, 0.5, 1, 1.5, 5};
    double ptBinWidthIcemFonll[] = {1, 1, 1, 1, 1, 2, 3, 10};

    /////////////////////////
    // J/psi               //
    /////////////////////////
    // FONLL
    double csJpsiFonll3Pt[] = {1.45E+02, 3.03E+02, 2.80E+02, 1.99E+02, 1.28E+02, 1.32E+02, 6.96E+01, 3.54E+01};

    double csJpsiFonll3PtMinS[] = {6.25E+01, 1.23E+02, 1.04E+02, 6.70E+01, 3.88E+01, 3.44E+01, 1.61E+01, 7.24E+00};
    double csJpsiFonll3PtMinM[] = {2.61E+01, 5.06E+01, 4.13E+01, 2.52E+01, 1.40E+01, 1.20E+01, 4.92E+00, 1.79E+00};
    double csJpsiFonll3PtMinPdf[] = {2.26E+01, 3.99E+01, 2.76E+01, 1.37E+01, 6.40E+00, 4.71E+00, 1.72E+00, 6.26E-01};

    double csJpsiFonll3PtMaxS[] = {6.76E+01, 1.35E+02, 1.21E+02, 8.50E+01, 5.41E+01, 5.33E+01, 2.54E+01, 1.06E+01};
    double csJpsiFonll3PtMaxM[] = {3.32E+01, 6.29E+01, 4.96E+01, 2.96E+01, 1.61E+01, 1.34E+01, 5.38E+00, 1.88E+00};
    double csJpsiFonll3PtMaxPdf[] = {2.26E+01, 3.99E+01, 2.74E+01, 1.37E+01, 6.43E+00, 4.69E+00, 1.72E+00, 6.26E-01};

    // ICEM
    double csJpsiIcemPt[] = {1052.186373, 2551.68636, 2250.286916, 1309.792997, 646.6338198, 479.9947889, 190.9539356, 95.38180569};

    double csJpsiIcemPtMinRS[] = {376.5235632,898.2031678, 770.3382638, 433.8606155, 207.2702479, 147.7996019, 55.45302669, 25.58886587};
    double csJpsiIcemPtMinMS[] = {409.1338159, 927.2768203, 802.7600056, 469.0383312, 230.079063, 164.0514961, 59.83116239, 27.63605107};

    double csJpsiIcemPtMaxRS[] = {808.7830856, 1893.638049, 1576.317214, 858.9261048, 397.9302254, 273.9155742, 98.0602908, 42.80021763};
    double csJpsiIcemPtMaxMS[] = {135.0836063, 282.5791416, 240.1768786, 139.4607135, 66.55449613, 44.48393196, 14.78375091, 6.408173883};

    // ICEM + FONLL
    double csJpsiIcemFonllPt[nPtBinsIcemFonll], csJpsiIcemFonllPtMin[nPtBinsIcemFonll], csJpsiIcemFonllPtMax[nPtBinsIcemFonll];

    for (int iPt = 0;iPt < nPtBinsIcemFonll;iPt++) {
        csJpsiIcemFonllPt[iPt] = (csJpsiFonll3Pt[iPt] + csJpsiIcemPt[iPt]) / ptBinWidthIcemFonll[iPt];
        csJpsiIcemFonllPtMin[iPt] = TMath::Sqrt(TMath::Power(csJpsiIcemPtMinRS[iPt], 2.) + TMath::Power(csJpsiIcemPtMinMS[iPt], 2.) + TMath::Power(csJpsiFonll3PtMinS[iPt], 2.) + TMath::Power(csJpsiFonll3PtMinM[iPt], 2.) + TMath::Power(csJpsiFonll3PtMinPdf[iPt], 2.)) / ptBinWidthIcemFonll[iPt];
        csJpsiIcemFonllPtMax[iPt] = TMath::Sqrt(TMath::Power(csJpsiIcemPtMaxRS[iPt], 2.) + TMath::Power(csJpsiIcemPtMaxMS[iPt], 2.) + TMath::Power(csJpsiFonll3PtMaxS[iPt], 2.) + TMath::Power(csJpsiFonll3PtMaxM[iPt], 2.) + TMath::Power(csJpsiFonll3PtMaxPdf[iPt], 2.)) / ptBinWidthIcemFonll[iPt];
    }

    /////////////////////////
    // Psi(2S)             //
    /////////////////////////
    // FONLL
    double csPsi2sFonll3Pt[] = {3.09E+01, 6.85E+01, 6.98E+01, 5.46E+01, 3.81E+01, 4.23E+01, 2.42E+01, 1.33E+01};

    double csPsi2sFonll3PtMinS[] = {1.36E+01, 2.85E+01, 2.68E+01, 1.91E+01, 1.21E+01, 1.16E+01, 5.74E+00, 2.80E+00};
    double csPsi2sFonll3PtMinM[] = {5.70E+00, 1.17E+01, 1.07E+01, 7.29E+00, 4.41E+00, 4.09E+00, 1.82E+00, 7.18E-01};
    double csPsi2sFonll3PtMinPdf[] = {5.07E+00, 9.68E+00, 7.49E+00, 4.14E+00, 2.09E+00, 1.64E+00, 6.40E-01, 2.46E-01};

    double csPsi2sFonll3PtMaxS[] = {1.48E+01, 3.07E+01, 3.01E+01, 2.35E+01, 1.63E+01, 1.75E+01, 9.14E+00, 4.15E+00};
    double csPsi2sFonll3PtMaxM[] = {7.31E+00, 1.48E+01, 1.30E+01, 8.61E+00, 5.09E+00, 4.59E+00, 1.99E+00, 7.55E-01};
    double csPsi2sFonll3PtMaxPdf[] = {5.07E+00, 9.76E+00, 7.49E+00, 4.14E+00, 2.09E+00, 1.63E+00, 6.38E-01, 2.46E-01};

    // ICEM
    double csPsi2sIcemPt[] = {147.1495306, 373.0399959, 343.9697861, 210.2373504, 110.301308, 88.29372064, 36.75509322, 17.34049981};

    double csPsi2sIcemPtMinRS[] = {51.82077191, 129.3969023, 116.4395633, 69.0612233, 35.13685404, 27.06648265, 10.6503437, 4.652349127};
    double csPsi2sIcemPtMinMS[] = {36.59749233, 78.01303169, 69.99020333, 45.17401467, 24.72548562, 19.76145483, 7.793657102, 3.457841352};

    double csPsi2sIcemPtMaxRS[] = {109.1685206, 268.5622632, 235.5037329, 135.6592502, 67.07627504, 49.97334259, 18.80430293, 7.783197068};
    double csPsi2sIcemPtMaxMS[] = {10.56748012, 20.52681649, 18.09462856, 11.9877668, 6.614028955, 5.192799981, 1.966125503, 0.8396444151};

    // ICEM + FONLL
    double csPsi2sIcemFonllPt[nPtBinsIcemFonll], csPsi2sIcemFonllPtMin[nPtBinsIcemFonll], csPsi2sIcemFonllPtMax[nPtBinsIcemFonll];

    for (int iPt = 0;iPt < nPtBinsIcemFonll;iPt++) {
        csPsi2sIcemFonllPt[iPt] = (csPsi2sFonll3Pt[iPt] + csPsi2sIcemPt[iPt]) / ptBinWidthIcemFonll[iPt];
        csPsi2sIcemFonllPtMin[iPt] = TMath::Sqrt(TMath::Power(csPsi2sIcemPtMinRS[iPt], 2.) + TMath::Power(csPsi2sIcemPtMinMS[iPt], 2.) + TMath::Power(csPsi2sFonll3PtMinS[iPt], 2.) + TMath::Power(csPsi2sFonll3PtMinM[iPt], 2.) + TMath::Power(csPsi2sFonll3PtMinPdf[iPt], 2.)) / ptBinWidthIcemFonll[iPt];
        csPsi2sIcemFonllPtMax[iPt] = TMath::Sqrt(TMath::Power(csPsi2sIcemPtMaxRS[iPt], 2.) + TMath::Power(csPsi2sIcemPtMaxMS[iPt], 2.) + TMath::Power(csPsi2sFonll3PtMaxS[iPt], 2.) + TMath::Power(csPsi2sFonll3PtMaxM[iPt], 2.) + TMath::Power(csPsi2sFonll3PtMaxPdf[iPt], 2.)) / ptBinWidthIcemFonll[iPt];
    }

    // Plot results
    TGraphAsymmErrors *graCsJpsiCgcNrqcFonlldPt = new TGraphAsymmErrors(nPtBinsCgcNrqcdFonll, ptCentrCgcNrqcdFonll, csJpsiCgcNrqcFonllPt, ptMinCgcNrqcdFonll, ptMaxCgcNrqcdFonll, csJpsiCgcNrqcFonllPtMin, csJpsiCgcNrqcFonllPtMax);
    graCsJpsiCgcNrqcFonlldPt -> SetFillColorAlpha(kOrange+7, 0.5);

    TGraphAsymmErrors *graCsPsi2sCgcNrqcFonlldPt = new TGraphAsymmErrors(nPtBinsCgcNrqcdFonll, ptCentrCgcNrqcdFonll, csPsi2sCgcNrqcFonllPt, ptMinCgcNrqcdFonll, ptMaxCgcNrqcdFonll, csPsi2sCgcNrqcFonllPtMin, csPsi2sCgcNrqcFonllPtMax);
    graCsPsi2sCgcNrqcFonlldPt -> SetFillColorAlpha(kOrange+7, 0.5);

    TGraphAsymmErrors *graCsJpsiNrqcFonlldPt = new TGraphAsymmErrors(nPtBinsNrqcdFonll, ptCentrNrqcdFonll, csJpsiNrqcFonllPt, ptMinNrqcdFonll, ptMaxNrqcdFonll, csJpsiNrqcFonllPtMin, csJpsiNrqcFonllPtMax);
    graCsJpsiNrqcFonlldPt -> SetFillColorAlpha(kAzure+2, 0.4);

    TGraphAsymmErrors *graCsPsi2sNrqcFonlldPt = new TGraphAsymmErrors(nPtBinsNrqcdFonll, ptCentrNrqcdFonll, csPsi2sNrqcFonllPt, ptMinNrqcdFonll, ptMaxNrqcdFonll, csPsi2sNrqcFonllPtMin, csPsi2sNrqcFonllPtMax);
    graCsPsi2sNrqcFonlldPt -> SetFillColorAlpha(kAzure+2, 0.4);

    TGraphAsymmErrors *graCsJpsiIcemFonlldPt = new TGraphAsymmErrors(nPtBinsIcemFonll, ptCentrIcemFonll, csJpsiIcemFonllPt, ptMinIcemFonll, ptMaxIcemFonll, csJpsiIcemFonllPtMin, csJpsiIcemFonllPtMax);
    graCsJpsiIcemFonlldPt -> SetFillStyle(3352);
    graCsJpsiIcemFonlldPt -> SetFillColorAlpha(kMagenta, 0.7);

    TGraphAsymmErrors *graCsPsi2sIcemFonlldPt = new TGraphAsymmErrors(nPtBinsIcemFonll, ptCentrIcemFonll, csPsi2sIcemFonllPt, ptMinIcemFonll, ptMaxIcemFonll, csPsi2sIcemFonllPtMin, csPsi2sIcemFonllPtMax);
    graCsPsi2sIcemFonlldPt -> SetFillStyle(3352);
    graCsPsi2sIcemFonlldPt -> SetFillColorAlpha(kMagenta, 0.7);

    TCanvas *canvasCsJpsiTheorPt = new TCanvas("canvasCsJpsiTheorPt", "", 800, 600);
    canvasCsJpsiTheorPt -> SetFillColor(0);
    canvasCsJpsiTheorPt -> SetBorderMode(0);
    canvasCsJpsiTheorPt -> SetBorderSize(0);
    canvasCsJpsiTheorPt -> SetTickx(1);
    canvasCsJpsiTheorPt -> SetTicky(1);
    canvasCsJpsiTheorPt -> SetLeftMargin(0.15);
    canvasCsJpsiTheorPt -> SetBottomMargin(0.1518219);
    canvasCsJpsiTheorPt -> SetFrameBorderMode(0);
    canvasCsJpsiTheorPt -> SetFrameBorderMode(0);
    gPad -> SetLogy(1);

    TH2D *histGridCsJpsiTheorPt = new TH2D("histGridCsJpsiTheorPt","",100, 0., 20., 100, 1, 1e5);
    histGridCsJpsiTheorPt -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGridCsJpsiTheorPt -> GetXaxis() -> SetTitleOffset(1.2);
    histGridCsJpsiTheorPt -> GetYaxis() -> SetTitle("d^{2}#sigma^{J/#psi} / d#it{p}_{T} d#it{y} (nb GeV^{-1})");
    histGridCsJpsiTheorPt -> GetYaxis() -> SetTitleOffset(1.3);
    histGridCsJpsiTheorPt -> Draw();

    graCsJpsiCgcNrqcFonlldPt -> Draw("E2 SAME");
    graCsJpsiNrqcFonlldPt -> Draw("E2 SAME");
    graCsJpsiIcemFonlldPt -> Draw("E2 SAME");

    TLegend *legendCsJpsiTheorPt = new TLegend(0.40, 0.65, 0.60, 0.80);
    SetLegend(legendCsJpsiTheorPt);
    legendCsJpsiTheorPt -> AddEntry(graCsJpsiCgcNrqcFonlldPt, "CGC + NRQCD + FONLL", "F");
    legendCsJpsiTheorPt -> AddEntry(graCsJpsiNrqcFonlldPt, "NRQCD + FONLL", "F");
    legendCsJpsiTheorPt -> AddEntry(graCsJpsiIcemFonlldPt, "ICEM + FONLL", "F");
    legendCsJpsiTheorPt -> Draw("SAME");


    TCanvas *canvasCsPsi2sTheorPt = new TCanvas("canvasCsPsi2sTheorPt", "", 800, 600);
    canvasCsPsi2sTheorPt -> SetFillColor(0);
    canvasCsPsi2sTheorPt -> SetBorderMode(0);
    canvasCsPsi2sTheorPt -> SetBorderSize(0);
    canvasCsPsi2sTheorPt -> SetTickx(1);
    canvasCsPsi2sTheorPt -> SetTicky(1);
    canvasCsPsi2sTheorPt -> SetLeftMargin(0.15);
    canvasCsPsi2sTheorPt -> SetBottomMargin(0.1518219);
    canvasCsPsi2sTheorPt -> SetFrameBorderMode(0);
    canvasCsPsi2sTheorPt -> SetFrameBorderMode(0);
    gPad -> SetLogy(1);

    TH2D *histGridCsPsi2sTheorPt = new TH2D("histGridCsPsi2sTheorPt","",100, 0., 20., 100, 1, 1e5);
    histGridCsPsi2sTheorPt -> GetXaxis() -> SetTitle("#it{p}_{T} (GeV/#it{c})");
    histGridCsPsi2sTheorPt -> GetXaxis() -> SetTitleOffset(1.2);
    histGridCsPsi2sTheorPt -> GetYaxis() -> SetTitle("d^{2}#sigma^{#psi(2S)} / d#it{p}_{T} d#it{y} (nb GeV^{-1})");
    histGridCsPsi2sTheorPt -> GetYaxis() -> SetTitleOffset(1.3);
    histGridCsPsi2sTheorPt -> Draw();

    graCsPsi2sCgcNrqcFonlldPt -> Draw("E2 SAME");
    graCsPsi2sNrqcFonlldPt -> Draw("E2 SAME");
    graCsPsi2sIcemFonlldPt -> Draw("E2 SAME");

    TLegend *legendCsPsi2sTheorPt = new TLegend(0.40, 0.65, 0.60, 0.80);
    SetLegend(legendCsPsi2sTheorPt);
    legendCsPsi2sTheorPt -> AddEntry(graCsPsi2sCgcNrqcFonlldPt, "CGC + NRQCD + FONLL", "F");
    legendCsPsi2sTheorPt -> AddEntry(graCsPsi2sNrqcFonlldPt, "NRQCD + FONLL", "F");
    legendCsPsi2sTheorPt -> AddEntry(graCsPsi2sIcemFonlldPt, "ICEM + FONLL", "F");
    legendCsPsi2sTheorPt -> Draw("SAME");
}
////////////////////////////////////////////////////////////////////////////////
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
}
////////////////////////////////////////////////////////////////////////////////
void SetLegend(TLegend *legend){
    legend -> SetBorderSize(0);
    legend -> SetFillColor(10);
    legend -> SetFillStyle(1);
    legend -> SetLineStyle(0);
    legend -> SetLineColor(0);
    legend -> SetTextFont(42);
    legend -> SetTextSize(0.045);
}