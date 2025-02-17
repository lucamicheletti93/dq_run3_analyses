import math
import ROOT
import numpy as np
from array import array
import statistics
from ROOT import TCanvas, TGraphMultiErrors, TLegend, TObject, TF1, TFitResult, TMatrixDSym

def multiple_y_errors():

    # Define the number of points
    np_Jpsi = 4
    np_psi = 4
    np_ratio = 4

    # PUBLISED RESULTS: https://www.hepdata.net/record/77781
    x = array('d', [5.02, 7, 8, 13])
    y = np.array([3.740* 1.5, 4.460* 1.5, 5.9866667* 1.5, 7.5333333* 1.5]) 
    exh = array('d', [0.3, 0.3, 0.3, 0.3, 0.3, 0.3])
    eyhstat = np.array([0.053333* 1.5, 0.0266667* 1.5, 0.0266667* 1.5, 0.0266667* 1.5]) 
    eyhsys = np.array([0.18551486* 1.5, 0.41927363* 1.5, 0.54566111* 1.5, 0.49156741* 1.5]) 
    eyhsys_BR = np.array([0.02070458* 1.5, 0.02469049* 1.5, 0.03314209* 1.5, 0.04170441* 1.5])

    #Same for psi(2S)
    x_psi = array('d', [5.02, 7, 8, 13])
    exh_psi = array('d', [0.3, 0.3, 0.3, 0.3])
    y_psi = np.array([0.48* 1.5, 0.753333* 1.5, 0.820* 1.5, 1.1* 1.5]) 
    eyhstat_psi = np.array([0.106666* 1.5, 0.0466667* 1.5, 0.053333* 1.5, 0.04* 1.5]) 
    eyhsys_psi = np.array([0.03813338* 1.5, 0.09316058* 1.5, 0.11306742* 1.5, 0.09947756* 1.5]) 
    eyhsys_psi_BR = np.array([0.05468354* 1.5, 0.08582278* 1.5, 0.09341772* 1.5, 0.12531646* 1.5])

    x_ratio = array('d', [5.02, 7, 8, 13])
    exh_ratio = array('d', [0.3, 0.3, 0.3, 0.3])
    y_ratio_list = []
    eyhstat_ratio_list = []
    eyhsys_ratio_list = []
    eyhsys_BR_list = []
    BR_jpsi = 0.05961*1.5
    BR_psi = 0.0077*1.5
    BR_ratio = BR_jpsi/BR_psi
    BR_jpsi_err = 0.00033
    BR_psi_err = 0.0008
    eyhstat_BR = BR_ratio * np.sqrt((BR_jpsi_err / BR_jpsi)**2 + (BR_psi_err / BR_psi)**2)
    print("BR_ratio: ", BR_ratio, " +/- ", eyhstat_BR)

    # Loop to compute ratio
    for energy in range(len(y_psi)):
        ratio = y_psi[energy] / y[energy]
        print("ratio: ", ratio)
        y_ratio_list.append(ratio)
        print("ratio list: ", y_ratio_list[energy])

        # Computation of statistical error
        eyhstat_ratio_list.append(ratio * np.sqrt((eyhstat[energy] / y[energy])**2 + (eyhstat_psi[energy] / y_psi[energy])**2))

        # Computation of systematic error
        eyhsys_ratio_list.append(ratio * np.sqrt((eyhsys[energy] / y[energy])**2 + (eyhsys_psi[energy] / y_psi[energy])**2))

        # Computation of BR systematic error
        eyhsys_BR_list.append(BR_ratio * np.sqrt((eyhsys_BR[energy] / BR_jpsi)**2 + (eyhsys_psi_BR[energy] / BR_psi)**2))

    y_ratio = array('d', y_ratio_list)
    #print("ratio array: ", y_ratio)
    eyhstat_ratio = array('d', eyhstat_ratio_list)
    #print("ratio stat array: ", eyhstat_ratio)
    eyhsys_ratio = array('d', eyhsys_ratio_list)
    #print("ratio syst array: ", eyhsys_ratio)
    eyhsys_BR = array('d', eyhsys_BR_list)
    #print("BR syst array: ", eyhsys_BR)

    #
    # Calculate combined systematic errors for published results
    combined_eyhsys = np.sqrt(np.add(np.square(eyhsys), np.square(eyhsys_BR)))
    combined_eyhsys_psi = np.sqrt(np.add(np.square(eyhsys_psi), np.square(eyhsys_psi_BR)))

    # Convert y and error arrays back to ROOT array format
    y = array('d', y)
    eyhstat = array('d', eyhstat)
    combined_eyhsys = array('d', combined_eyhsys)

    y_psi = array('d', y_psi)
    eyhstat_psi = array('d', eyhstat_psi)
    combined_eyhsys_psi = array('d', combined_eyhsys_psi)

    # Create the TGraphMultiErrors object
    gme = ROOT.TGraphMultiErrors(np_Jpsi, x, y, exh, exh, eyhstat, eyhstat)
    gme.AddYError(np_Jpsi, combined_eyhsys, combined_eyhsys)

    # Set styles
    gme.SetMarkerStyle(20)
    gme.SetMarkerSize(0.5)
    gme.SetMarkerColor(ROOT.kBlue)
    gme.SetLineWidth(0)
    gme.GetAttLine(0).SetLineColor(ROOT.kBlue)
    gme.GetAttLine(1).SetLineColor(ROOT.kBlue)
    gme.GetAttFill(1).SetFillColorAlpha(ROOT.kBlue, 0.1)

    gme.SetTitle("Cross Section vs Center of Mass Energy")
    gme.GetXaxis().SetTitle("#sqrt{#it{s}} (TeV)")
    gme.GetYaxis().SetTitle("#sigma (#mub)")
    gme.GetXaxis().SetLimits(2, 15)
    gme.GetXaxis().SetRangeUser(2, 15)
    gme.GetYaxis().SetRangeUser(2, 15)


    # Create the TGraphMultiErrors object for psi
    gme_psi = ROOT.TGraphMultiErrors(np_psi, x_psi, y_psi, exh_psi, exh_psi, eyhstat_psi, eyhstat_psi)
    gme_psi.AddYError(np_psi, combined_eyhsys_psi, combined_eyhsys_psi)

    # Set styles
    gme_psi.SetMarkerStyle(20)
    gme_psi.SetMarkerSize(0.5)
    gme_psi.SetMarkerColor(ROOT.kRed)
    gme_psi.SetLineWidth(0)
    gme_psi.GetAttLine(0).SetLineColor(ROOT.kRed)
    gme_psi.GetAttLine(1).SetLineColor(ROOT.kRed)
    gme_psi.GetAttFill(1).SetFillColorAlpha(ROOT.kRed, 0.1)

    gme_psi.SetTitle("Cross Section vs Center of Mass Energy")
    gme_psi.GetXaxis().SetTitle("#sqrt{#it{s}} (TeV)")
    gme_psi.GetYaxis().SetTitle("#sigma (#mub)")
    gme_psi.GetXaxis().SetLimits(2, 15)
    gme_psi.GetXaxis().SetRangeUser(2, 15)
    gme_psi.GetYaxis().SetRangeUser(0.4, 2.3)

    # Create the TGraphMultiErrors object for ratio
    gme_ratio = ROOT.TGraphMultiErrors(np_ratio, x_ratio, y_ratio, exh_ratio, exh_ratio, eyhstat_ratio, eyhstat_ratio)
    gme_ratio.AddYError(np_ratio, eyhsys_ratio, eyhsys_ratio)


    # Set styles
    gme_ratio.SetMarkerStyle(20)
    gme_ratio.SetMarkerSize(0.5)
    gme_ratio.SetMarkerColor(ROOT.kGray+2)
    gme_ratio.SetLineWidth(0)
    gme_ratio.GetAttLine(0).SetLineColor(ROOT.kGray+2)
    gme_ratio.GetAttLine(1).SetLineColor(ROOT.kGray+2)
    gme_ratio.GetAttFill(1).SetFillColorAlpha(ROOT.kGray+2, 0.1)

    gme_ratio.SetTitle("Cross Section Ratio vs Center of Mass Energy")
    gme_ratio.GetXaxis().SetTitle("#sqrt{#it{s}} (TeV)")
    gme_ratio.GetYaxis().SetTitle("#sigma_{#psi(2S)}/#sigma_{J/#psi}")
    gme_ratio.GetXaxis().SetLimits(2, 15)
    gme_ratio.GetXaxis().SetRangeUser(2, 15)

    # Create a canvas
    c1 = ROOT.TCanvas("c1", "Jpsi with multiple y-errors", 600, 500)
    c1.GetFrame().SetBorderSize(12)
    gme.Draw("APS ; Z ; 5 s=0.5")

    # Update the canvas to display the plot
    c1.Update()
    c1.Draw()

    # Add legend
    legend = TLegend(0.13, 0.73, 0.55, 0.89)
    legend.SetHeader("ALICE pp, 2.5 < #it{y} < 4.0", "c")
    legend.AddEntry(gme, " J/#psi, published results")
    legend.SetLineWidth(0)
    legend.Draw("same")

    # Save the plot
    c1.SaveAs("JPsivsSqrt(s).pdf")
    c1.SaveAs("JPsivsSqrt(s).root")

    # Create a canvas for psi(2S)
    c2 = ROOT.TCanvas("c2", "psi(2S) with multiple y-errors", 600, 500)
    c2.GetFrame().SetBorderSize(12)
    gme_psi.GetXaxis().SetRangeUser(2, 15)
    gme_psi.Draw("APS ; Z ; 5 s=0.5")

    # Update the canvas to display the plot
    c2.Update()
    c2.Draw()

    # Add legend
    legend_2 = TLegend(0.13, 0.73, 0.55, 0.89)
    legend_2.SetHeader("ALICE pp, 2.5 < #it{y} < 4.0", "c")
    legend_2.AddEntry(gme_psi, "#psi(2S), published results")
    legend_2.SetLineWidth(0)
    legend_2.Draw("same")

    # Save the plot
    c2.SaveAs("psi(2S)vsSqrt(s).pdf")
    c2.SaveAs("psi(2S)vsSqrt(s).root")


    #FITS
    fPol1 = TF1("fPol1", "pol1", 2, 15)
    fPol1.SetLineColor(ROOT.kMagenta+1)
    fExpo = TF1("fExpo", "[0]*(1 - exp(-[1]*x))", 2, 15)
    fExpo.SetParLimits(0, 0, 50)
    fExpo.SetParLimits(1, 0, 10)
    fExpo.SetLineColor(ROOT.kGreen+1)
    fPowerLaw = TF1("fPowerLaw", "[0]* x**([1])", 2, 15)
    fPowerLaw.SetLineColor(ROOT.kOrange+7)

    # Create a canvas for Ratio
    c3 = ROOT.TCanvas("c2", "psi(2S) with multiple y-errors", 600, 500)
    c3.GetFrame().SetBorderSize(12)
    gme_ratio.GetXaxis().SetRangeUser(2, 15)
    gme_ratio.GetYaxis().SetRangeUser(0.095, 0.23)
    gme_ratio.Draw("APS ; Z ; 5 s=0.5")
    gme_ratio.Fit("fPol1", "0", "", 4, 13.1)
    fPol1.Draw("same")
    gme_ratio.Fit("fExpo", "0", "", 4, 13.1)
    fExpo.Draw("same")
    gme_ratio.Fit("fPowerLaw", "0", "", 4, 13.1)
    fPowerLaw.Draw("same")
    fPol1.Draw("same")

    # Add legend
    legend_fit = TLegend(0.13, 0.73, 0.58, 0.89)
    legend_fit.AddEntry(gme_ratio, " Published results")
    legend_fit.AddEntry(fPol1, "pol1")
    legend_fit.AddEntry(fPowerLaw, "power law")
    legend_fit.AddEntry(fExpo, "expo")
    legend_fit.SetLineWidth(0)
    legend_fit.Draw("same")


    # Update the canvas to display the plot
    c3.Update()
    c3.Draw()

    # Save the plot
    c3.SaveAs("RatiovsSqrt(s).pdf")
    c3.SaveAs("RatiovsSqrt(s).root")

    #-------- Ratio prediction from [5.02 - 13] TeV data --------------------
    print("-------- Ratio prediction from [5.02 - 13] TeV data --------------------")
    err_0 = fPol1.GetParError(0)
    err_1 = fPol1.GetParError(1)
    pol1_err = math.sqrt( err_0*err_0 + 5.36*5.36 * err_1*err_1)
    print("pol1: ", fPol1.Eval(5.36), " +/- ", math.sqrt( err_0*err_0 + 5.36*5.36 * err_1*err_1))
    par0 = fPowerLaw.GetParameter(0)
    par1 = fPowerLaw.GetParameter(1)
    par2 = fPowerLaw.GetParameter(2)
    err_0 = fPowerLaw.GetParError(0)
    err_1 = fPowerLaw.GetParError(1)
    err_2 = fPowerLaw.GetParError(2)
    fPowerLaw_err = math.sqrt( 5.36**(2*par1)*err_0*err_0 + (par0*5.36**(par1)*math.log(5.36))*(par0*5.36**(par1)*math.log(5.36))*err_1*err_1 + err_2*err_2)
    print("power law: ", fPowerLaw.Eval(5.36), " +/- ", math.sqrt( 5.36**(2*par1)*err_0*err_0 + (par0*5.36**(par1)*math.log(5.36))*(par0*5.36**(par1)*math.log(5.36))*err_1*err_1 + err_2*err_2))
    err_0 = fExpo.GetParError(0)
    err_1 = fExpo.GetParError(1)
    par0 = fExpo.GetParameter(0)
    par1 = fExpo.GetParameter(1)
    term1 = (1 - math.exp(-par1 * 5.36))**2 * err_0**2
    term2 = (par0 * par1 * math.exp(-par1 * 5.36))**2 * err_1**2
    expo_err = math.sqrt(term1 + term2)
    print("expo: ", fExpo.Eval(5.36), " +/- ", expo_err)
    average_ratio = (fPol1.Eval(5.36) + fPowerLaw.Eval(5.36) + fExpo.Eval(5.36))/3
    average_ratio_stat_err = math.sqrt(pol1_err*pol1_err + fPowerLaw_err*fPowerLaw_err + expo_err*expo_err)/3
    average_ratio_syst_err = (statistics.stdev([fPol1.Eval(5.36), fPowerLaw.Eval(5.36), fExpo.Eval(5.36)])/math.sqrt(3))/BR_ratio
    print("average ratio: ",  average_ratio, " +/- ", average_ratio_stat_err/BR_ratio,  " +/- ", average_ratio_syst_err)
    print("BR_ratio: ", BR_ratio, " +/- ", eyhstat_BR)


    # Create a canvas
    c1_fit = ROOT.TCanvas("c1_fit", "Jpsi with multiple y-errors", 600, 500)
    c1_fit.GetFrame().SetBorderSize(12)
    gme.GetXaxis().SetLimits(0, 15)
    gme.GetXaxis().SetRangeUser(0, 15)
    gme.GetYaxis().SetLimits(0, 15)
    gme.GetYaxis().SetRangeUser(0, 15)
    gme.Draw("APS ; Z ; 5 s=0.5")
    gme.Fit("fPol1", "0", "", 6.75, 13.1)
    fPol1.Draw("same")
    gme.Fit("fExpo", "0", "", 2.75, 13.1)
    fExpo.Draw("same")
    gme.Fit("fPowerLaw", "0", "", 2.75, 13.1)
    fPowerLaw.Draw("same")
    fPol1.Draw("same")

    # Add legend
    legend_fit = TLegend(0.13, 0.73, 0.58, 0.89)
    legend_fit.AddEntry(gme, " J/#psi")
    legend_fit.AddEntry(fPol1, "pol1")
    legend_fit.AddEntry(fPowerLaw, "power law")
    legend_fit.AddEntry(fExpo, "expo")
    legend_fit.SetLineWidth(0)
    legend_fit.Draw("same")

    # Update the canvas to display the plot
    c1_fit.Update()
    c1_fit.Draw()

    # Save the plot
    c1_fit.SaveAs("JPsivsSqrt(s)_fit.pdf")

    #-------- Jpsi prediction from [2.76 - 13] TeV data --------------------
    print("-------- Jpsi prediction from [2.76 - 13] TeV data --------------------")
    err_0 = fPol1.GetParError(0)
    err_1 = fPol1.GetParError(1)
    print("pol1: ", fPol1.Eval(5.36), " +/- ", math.sqrt( err_0*err_0 + 5.36*5.36 * err_1*err_1))
    par0 = fPowerLaw.GetParameter(0)
    par1 = fPowerLaw.GetParameter(1)
    par2 = fPowerLaw.GetParameter(2)
    err_0 = fPowerLaw.GetParError(0) #
    err_1 = fPowerLaw.GetParError(1)
    err_2 = fPowerLaw.GetParError(2)
    print("power law: ", fPowerLaw.Eval(5.36), " +/- ", math.sqrt( 5.36**(2*par1)*err_0*err_0 + (par0*5.36**(par1)*math.log(5.36))*(par0*5.36**(par1)*math.log(5.36))*err_1*err_1 + err_2*err_2))
    err_0 = fExpo.GetParError(0)
    err_1 = fExpo.GetParError(1)
    err_2 = fExpo.GetParError(2)
    par0 = fExpo.GetParameter(0)
    par1 = fExpo.GetParameter(1)
    print("expo: ", fExpo.Eval(5.36), " +/- ", math.sqrt( (1-math.exp(-par1*5.36))*(1-math.exp(-par1*5.36)) * err_0*err_0 + (par0*par1*math.exp(-par1*5.36))*(par0*par1*math.exp(-par1*5.36))* err_1*err_1 + err_2*err_2))


    # Create a canvas
    c2_fit = ROOT.TCanvas("c2_fit", "Psi(2S) with multiple y-errors", 600, 500)
    c2_fit.GetFrame().SetBorderSize(12)
    gme_psi.GetXaxis().SetLimits(0, 15)
    gme_psi.GetXaxis().SetRangeUser(0, 15)
    gme_psi.GetYaxis().SetLimits(0, 2.25)
    gme_psi.GetYaxis().SetRangeUser(0, 2.25)
    gme_psi.Draw("APS ; Z ; 5 s=0.5")
    gme_psi.Fit("fPol1", "0", "", 6.75, 13.1)
    fPol1.Draw("same")
    fit_results_expo = gme_psi.Fit("fExpo", "0", "", 5, 13.1)
    fExpo.SetParameters(2.74, 7.35e-2)
    fExpo.Draw("same")
    gme_psi.Fit("fPowerLaw", "0", "", 5, 13.1)
    fPowerLaw.Draw("same")
    fPol1.Draw("same")


    # Add legend
    legend_fit_psi = TLegend(0.13, 0.73, 0.58, 0.89)
    legend_fit_psi.AddEntry(gme_psi, " #psi(2S)")
    legend_fit_psi.AddEntry(fPol1, "pol1")
    legend_fit_psi.AddEntry(fPowerLaw, "power law")
    legend_fit_psi.AddEntry(fExpo, "expo")
    legend_fit_psi.SetLineWidth(0)
    legend_fit_psi.Draw("same")

    # Update the canvas to display the plot
    c2_fit.Update()
    c2_fit.Draw()

    # Save the plot
    c2_fit.SaveAs("Psi(2S)vsSqrt(s)_fit.pdf")

    #-------- Psi(2S) prediction from [5.02 - 13] TeV data --------------------
    print("-------- Psi(2S) prediction from [5.02 - 13] TeV data --------------------")
    err_0 = fPol1.GetParError(0)
    err_1 = fPol1.GetParError(1)
    print("pol1: ", fPol1.Eval(5.36), " +/- ", math.sqrt( err_0*err_0 + 5.36*5.36 * err_1*err_1))
    par0 = fPowerLaw.GetParameter(0)
    par1 = fPowerLaw.GetParameter(1)
    par2 = fPowerLaw.GetParameter(2)
    err_0 = fPowerLaw.GetParError(0)
    err_1 = fPowerLaw.GetParError(1)
    err_2 = fPowerLaw.GetParError(2)
    print("power law: ", fPowerLaw.Eval(5.36), " +/- ", math.sqrt( 5.36**(2*par1)*err_0*err_0 + (par0*5.36**(par1)*math.log(5.36))*(par0*5.36**(par1)*math.log(5.36))*err_1*err_1 + err_2*err_2))
    err_0 = fExpo.GetParError(0)
    err_1 = fExpo.GetParError(1)
    err_2 = fExpo.GetParError(2)
    #cov_01 = cov_matrix_expo[0][1]
    par0 = fExpo.GetParameter(0)
    par1 = fExpo.GetParameter(1)
    print("expo: ", fExpo.Eval(5.36), " +/- ", math.sqrt( (1-math.exp(-par1*5.36))*(1-math.exp(-par1*5.36)) * err_0*err_0 + (par0*par1*math.exp(-par1*5.36))*(par0*par1*math.exp(-par1*5.36))* err_1*err_1 +err_2*err_2))


# Execute the function
multiple_y_errors()