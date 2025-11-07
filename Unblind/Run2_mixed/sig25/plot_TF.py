from __future__ import print_function, division
import sys
import os
import csv
import rhalphalib as rl
from rhalphalib import AffineMorphTemplate, MorphHistW2
import numpy as np
import scipy.stats
import pickle
import ROOT
from ROOT import *
rl.util.install_roofit_helpers()
#rl.ParametericSample.PreferRooParametricHist = False
import pandas as pd
import json
import matplotlib.pyplot as plt
import mplhep as hep

def AddCMSLumi(pad, fb, extra):
	cmsText     = "CMS " + extra
	cmsTextFont   = 61  
	lumiTextSize     = 0.45
	lumiTextOffset   = 0.15
	cmsTextSize      = 0.5
	cmsTextOffset    = 0.15
	H = pad.GetWh()
	W = pad.GetWw()
	l = pad.GetLeftMargin()
	t = pad.GetTopMargin()
	r = pad.GetRightMargin()
	b = pad.GetBottomMargin()
	e = 0.025
	pad.cd()
	lumiText = str(fb)+" fb^{-1} (13 TeV)"
	latex = TLatex()
	latex.SetNDC()
	latex.SetTextAngle(0)
	latex.SetTextColor(kBlack)	
	extraTextSize = 0.76*cmsTextSize
	latex.SetTextFont(42)
	latex.SetTextAlign(31) 
	latex.SetTextSize(lumiTextSize*t)	
	latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText)
	pad.cd()
	latex.SetTextFont(cmsTextFont)
	latex.SetTextSize(cmsTextSize*t)
	latex.SetTextAlign(11)
	latex.DrawLatex(0.1265, 0.825, cmsText)
	pad.Update()

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-y', type=str, dest='year',
                        action='store', default='all',
                        help='Run 2 year')
    parser.add_argument('-mc', '--MC', action='store_true', default=False, help='TF_MC?')


    args = parser.parse_args()

    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)

#    plt.style.use([hep.style.ROOT, {"font.size": 24}])
#    plt.switch_backend("agg")

    ROOT.gInterpreter.Declare("Double_t widebins_2016[6] = {200, 240, 275, 330, 700};")
    ROOT.gInterpreter.Declare("Double_t widebins_2017[6] = {220, 256, 300, 360, 700};")
    ROOT.gInterpreter.Declare("Double_t widebins_2018[6] = {120, 160, 200, 255, 700};")

    if str(args.year) == 'all': years = ['2016', '2017', '2018']
    else: years = [str(args.year)]

    ptBins = {}
    ptBins['2016'] = np.array([200, 240, 275, 330, 700])
    ptBins['2017'] = np.array([220, 260, 300, 360, 700])
    ptBins['2018'] = np.array([120, 160, 200, 255, 700])

    rho_deg = {}
    if args.MC: rho_deg['2016'] = 4
    else: rho_deg['2016'] = 3
    rho_deg['2017'] = 4
    rho_deg['2018'] = 4

    pT_deg = {}
    if args.MC: 
        pT_deg['2016'] = 4
        pT_deg['2017'] = 4
        pT_deg['2018'] = 4
    else: 
        pT_deg['2016'] = 2
        pT_deg['2017'] = 2
        pT_deg['2018'] = 2

    LUMI = {}
    LUMI['2016'] = 36.31
    LUMI['2017'] = 41.48
    LUMI['2018'] = 59.82
	
    cmsextra = "Preliminary"

    PFratio = {}
    PFratio['2016'] = 0.03582468151021236
    PFratio['2017'] = 0.03278941795956571
    PFratio['2018'] = 0.011821337309360778
    

    TF = {}
    TF['2016'] = TH2F("TF_2016", "Jet pT vs. Softdrop Mass Transfer Factor", 4, widebins_2016, 40, 0, 200)
    TF['2017'] = TH2F("TF_2017", "Jet pT vs. Softdrop Mass Transfer Factor", 4, widebins_2017, 40, 0, 200)
    TF['2018'] = TH2F("TF_2018", "Jet pT vs. Softdrop Mass Transfer Factor", 4, widebins_2018, 40, 0, 200)

    for y in years:
        TF[y].SetXTitle("Jet pT")
        TF[y].SetYTitle("Softdrop Mass")

        if args.MC: TF[y].SetTitle(y+" TF_MC Plotted in Jet pT vs. Softdrop Mass")
        else: TF[y].SetTitle(y+" TF_Residual Plotted in Jet pT vs. Softdrop Mass")

        pts = ptBins[y]
        msdbins = np.linspace(0, 200, 41)
        ptpts, msdpts = np.meshgrid(pts[:-1] + 0.5 * np.diff(pts), msdbins[:-1] + 0.5 * np.diff(msdbins), indexing='ij')

        ptpts_scaled = (ptpts - pts[0]) / (pts[-1] - pts[0])

        rhopts = 2*np.log(msdpts/ptpts)

        if y == '2018': rhopts_scaled = (rhopts - (-7.3)) / ((-2.0) - (-7.3))
        else: rhopts_scaled = (rhopts - (-7.0)) / ((-2.0) - (-7.0))
        

        validbins = (rhopts_scaled >= 0) & (rhopts_scaled <= 1)

        ptpts = ptpts[validbins]
        msdpts = msdpts[validbins]
        rhopts = rhopts[validbins]
        ptpts_scaled = ptpts_scaled[validbins]
        rhopts_scaled = rhopts_scaled[validbins]

        print(y+" pT degree: "+str(pT_deg[y]))
        print(y+" rho degree: "+str(rho_deg[y]))

        # initial values
        if args.MC: print(y+' Initial fit values read from file fit_TFMC_'+y+'.csv')
        else: print(y+' Initial fit values read from file fit_b_final_'+y+'.csv')
        if args.MC: initial_vals = np.genfromtxt('fit_b_TFMC_'+y+'.csv')
        else: initial_vals = np.genfromtxt('fit_b_final_'+y+'.csv')
        initial_vals = initial_vals.reshape(pT_deg[y]+1, rho_deg[y]+1)
        print(initial_vals)

        tf_dataResidual = rl.BasisPoly('CMS_EXO24027_tf_dataResidual_'+y, (pT_deg[y], rho_deg[y]), ['pt', 'rho'] , init_params=initial_vals, basis='Bernstein', limits=(-100, 100))

        for i in range(1, TF[y].GetNbinsX()+1):
            pTval = TF[y].GetXaxis().GetBinCenter(i)
            pTval_scaled = (pTval - pts[0]) / (pts[-1] - pts[0])
            for j in range(1, TF[y].GetNbinsY()+1):
                massVal = TF[y].GetYaxis().GetBinCenter(j)
                rhoVal = 2*np.log(massVal/pTval)
                if y == '2018': rhoVal_scaled = (rhoVal - (-7.3)) / ((-2.0) - (-7.3))
                else: rhoVal_scaled = (rhoVal - (-7.0)) / ((-2.0) - (-7.0))
                
                if rhoVal_scaled <= 0 or rhoVal_scaled >= 1: TF[y].SetBinContent(i, j, 0) # skip bins outside valid range
                else:
                    if y == '2018'and i==4: TF[y].SetBinContent(i, j, tf_dataResidual(pTval_scaled, rhoVal_scaled, nominal=True)*0.4)
                    else: TF[y].SetBinContent(i, j, tf_dataResidual(pTval_scaled, rhoVal_scaled, nominal=True))
        
        c1 = TCanvas()
        c1.cd()
        TF[y].Draw("LEGO2Z")
#        AddCMSLumi(ROOT.gPad, LUMI[y], cmsextra)

        if args.MC: c1.SaveAs(y+"TF_MC.png")
        else: c1.SaveAs(y+"TF_Residual.png")
        c1.Close()



    comment = '''
    # 2018 TF_res
    #pts = np.array([120, 160, 200, 255, 700])
    pts = np.array([120, 160, 200, 255, 700])
#    ptbins = np.array([120, 160, 200, 255, 700])
#    npt = len(ptbins) - 1
    msdbins = np.linspace(0, 200, 41)

    ptpts, msdpts = np.meshgrid(pts[:-1] + 0.5 * np.diff(pts), msdbins[:-1] + 0.5 * np.diff(msdbins), indexing='ij')


    ptpts_scaled= (ptpts - 120.0) / (700.0 - 120.0)
#    ptpts_scaled = (ptpts) / (2000.)
    rhopts = 2*np.log(msdpts/ptpts)

    rhopts_scaled = (rhopts - (-7.3)) / ((-2.0) - (-7.3))
    validbins = (rhopts_scaled >= 0) & (rhopts_scaled <= 1)

    ptpts = ptpts[validbins]
    msdpts = msdpts[validbins]
    rhopts = rhopts[validbins]
    ptpts_scaled = ptpts_scaled[validbins]
    rhopts_scaled = rhopts_scaled[validbins]

    # initial values
    print('2018 Initial fit values read from file fit_b_final_2018.csv')
    initial_vals_2018 = np.genfromtxt('fit_b_final_2018.csv')
    initial_vals_2018 = initial_vals_2018.reshape(3, 5)
    print(initial_vals_2018)

    tf_dataResidual_2018 = rl.BasisPoly('CMS_EXO24027_tf_dataResidual_2018', (2, 4), ['pt', 'rho'] , init_params=initial_vals_2018, basis='Bernstein', limits=(-100, 100))

    tf_dataResidual_2018_vals = tf_dataResidual_2018(ptpts_scaled, rhopts_scaled, nominal=True)
    '''
    comment = '''
    df = pd.DataFrame([])
    df['msd'] = msdpts.reshape(-1)
    df['pt'] = ptpts.reshape(-1)
    df['TF'] = tf_dataResidual_2018_vals.reshape(-1)
    
    mpt = df.to_numpy()

    fig, ax = plt.subplots()
    h = ax.hist2d(x=df["msd"],y=df["pt"],weights=df["TF"], bins=(msdbins,pts))
    plt.xlabel("$m_{sd}$ [GeV]")
    plt.ylabel("$p_{T}$ [GeV]")
    cb = fig.colorbar(h[3],ax=ax)
    cb.set_label("Ratio")
    fig.savefig("TF_msdpt_2018.png",bbox="tight")
    '''

    # Storing vals in ROOT TH2
    

