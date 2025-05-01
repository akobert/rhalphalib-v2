#Copied from https://github.com/osherson/B2GAnalysis_2020/blob/master/analysis/CombineStep1.py
#And modified for my analysis
import ROOT
from ROOT import *
import numpy
import math
import sys
import array
import os
import numpy as np
#import TEMPPAYLOAD
#from TEMPPAYLOAD import *

def GoodPlotFormat(H, *args): # Handy little script for color/line/fill/point/etc...
	try: H.SetStats(0)
	except: print " ------------ [  No stats box found!  ]"
	if args[0] == 'thickline':
		H.SetLineColor(args[1])
		H.SetLineWidth(2)
	if args[0] == 'thinline':
		H.SetLineColor(args[1])
		H.SetLineWidth(1)
		H.SetLineStyle(args[2])
	if args[0] == 'fill':
		H.SetLineColor(args[1])
		H.SetFillColor(args[2])
		H.SetFillStyle(args[3])
	if args[0] == 'markers':
		H.SetLineColor(args[1])
		H.SetMarkerColor(args[1])
		H.SetMarkerStyle(args[2])
	H.GetXaxis().SetTitleSize(0.04)

	if args[0] == 'dashed':
		H.SetLineColor(args[1])
		H.SetLineWidth(1)
		H.SetLineStyle(7)

def convertAsymGraph(TG, template, name):
	Hist = template.Clone(name)
	for i in range(1,Hist.GetNbinsX()+1):
		Hist.SetBinContent(i,0.)

#	print("Number of points: "+str(TG.GetN()))

	for i in range(TG.GetN()):
	#	Hist.SetBinContent(i+1,TG.GetY()[i]*(TG.GetErrorXlow(i)+TG.GetErrorXhigh(i)))
		Hist.SetBinContent(i+1,TG.GetY()[i])
#		print("Y= "+str(TG.GetY()[i]))
#		print("error xlow + error xhigh = "+str(TG.GetErrorXlow(i)+TG.GetErrorXhigh(i)))
		Hist.SetBinError(i+1, TG.GetErrorY(i))
	return Hist

def convertBinNHist(H, template, name):
	Hist = template.Clone(name)
	for i in range(1,Hist.GetNbinsX()+1):
		Hist.SetBinContent(i,H.GetBinContent(i))
		if H.GetBinContent(i) > 0 and (name == "dataprefitnewprefit" or name == "datafit_bnewfit_b" or name == "datafit_snewfit_s"):
	                Hist.SetBinError(i,math.sqrt(H.GetBinContent(i)))
		else:
			Hist.SetBinError(i,H.GetBinError(i))
	return Hist

def Reroll(H, VARS, ZOOM):
	X = TH1F("x_"+H.GetName(), ";"+VARS[2], len(VARS[1])-1, numpy.array(VARS[1]))
	Y = TH1F("y_"+H.GetName(), ";"+VARS[5], len(VARS[4])-1, numpy.array(VARS[4]))
		
	XY = TH2F("xy_"+H.GetName(), ";"+VARS[2]+";"+VARS[5], len(VARS[1])-1, numpy.array(VARS[1]), len(VARS[4])-1, numpy.array(VARS[4]))
	XYe = TH2F("xye_"+H.GetName(), ";"+VARS[2]+";"+VARS[5], len(VARS[1])-1, numpy.array(VARS[1]), len(VARS[4])-1, numpy.array(VARS[4]))

	X.SetStats(0)
	Y.SetStats(0)
	XY.SetStats(0)

	nxb = X.GetNbinsX()
	nyb = Y.GetNbinsX()

	for i in range(0,(nxb)):
		for j in range(0,(nyb)):
			k = XY.GetBin(i+1,j+1)
			index = H.FindBin(1 + j + i*nyb)
			XY.SetBinContent(k, H.GetBinContent(index))
			XY.SetBinError(k, (H.GetBinError(index)))
			XYe.SetBinContent(k, (H.GetBinError(index)))
	zoom_xl = XY.GetYaxis().FindBin(ZOOM[1][0])
	zoom_xh = XY.GetYaxis().FindBin(ZOOM[1][1])
	zoom_yl = XY.GetXaxis().FindBin(ZOOM[0][0])
	zoom_yh = XY.GetXaxis().FindBin(ZOOM[0][1])
	X = XY.ProjectionX("px_"+H.GetName(),1,nxb,"e")
	Y = XY.ProjectionY("py_"+H.GetName(),1,nyb,"e")
	zX = XY.ProjectionX("zpx_"+H.GetName(),zoom_xl,zoom_xh,"e")
	zY = XY.ProjectionY("zpy_"+H.GetName(),zoom_yl,zoom_yh,"e")
	return [X,Y,XY,XYe,zX,zY]

def MakeNBinsFromMinToMax(N,Min,Max): # helper for making large bin arrays makes N bins between Min and Max (same as you're feed to a TH1F)
	BINS = []
	for i in range(N+1):
		BINS.append(Min+(i*(Max-Min)/N))
	return BINS

def DBBW(H):
	for i in range(1,H.GetNbinsX()+1):
		C = H.GetBinContent(i)
		E = H.GetBinError(i)
		W = H.GetBinWidth(i)
		H.SetBinContent(i, C/W)
		H.SetBinError(i, E/W)	
	return H

def FindAndSetMax(*args):
	if len(args) == 1: args = args[0]
	maximum = 0.0
	for i in args:
		i.SetStats(0)
		t = i.GetMaximum()
		if t > maximum:
			maximum = t
	for j in args:
		j.GetYaxis().SetRangeUser(1,maximum*1.35)#should be 1.35 (below as well)
		j.SetLineWidth(2)
	return maximum*1.35
def AddCMSLumi(pad, fb, extra):
	#cmsText     = "CMS " + extra
	#cmsText     = "Private Work (CMS Data/Simulation)"
	cmsText     = "CMS Preliminary"
	cmsTextFont   = 61  
	lumiTextSize     = 0.45
	lumiTextOffset   = 0.15
	cmsTextSize      = 0.42
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
	#latex.DrawLatex(0.1265, 0.825, cmsText)
	latex.DrawLatex(0.165, 0.91, cmsText)
	pad.Update()

def pTbin(pad, Bin):
	cmsText     = "pT Bin: "+Bin+" GeV"
	cmsTextFont   = 61  
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
	latex = TLatex()
	latex.SetNDC()
	latex.SetTextAngle(0)
	latex.SetTextColor(kBlack)	
	extraTextSize = 0.76*cmsTextSize
	pad.cd()
	latex.SetTextFont(cmsTextFont)
	latex.SetTextSize(cmsTextSize*t)
	latex.SetTextAlign(11)
	latex.DrawLatex(0.28, 0.79, cmsText)
	pad.Update()

def GOF(pad, chi, ndf):
	chiText     = "Chi^2 = "+str(chi)
	ndfText     = "NDF = "+str(ndf)
	gofText     = "Chi^2/NDF = "+str(chi/ndf)
	cmsTextFont   = 61  
	cmsTextSize      = 0.42
	cmsTextOffset    = 0.15
	H = pad.GetWh()
	W = pad.GetWw()
	l = pad.GetLeftMargin()
	t = pad.GetTopMargin()
	r = pad.GetRightMargin()
	b = pad.GetBottomMargin()
	e = 0.025
	pad.cd()
	latex = TLatex()
	latex.SetNDC()
	latex.SetTextAngle(0)
	latex.SetTextColor(kBlack)	
	pad.cd()
	latex.SetTextFont(cmsTextFont)
	latex.SetTextSize(cmsTextSize*t)
	latex.SetTextAlign(11)
	#latex.DrawLatex(0.1265, 0.825, cmsText)
	latex.DrawLatex(0.34, 0.825, chiText)
	latex.DrawLatex(0.34, 0.785, ndfText)
	latex.DrawLatex(0.34, 0.745, gofText)
	pad.Update()

def Chi2(data, bkg, nInd):
	# data is data hist, bkg is total background, nInd is number of independent parameters
	# Returns chi2, NDF, and chi2/NDF
	bins = 0 # Non-zero bins
	chi_square = 0 # Running total
	for i in range(1,data.GetNbinsX()+1):
		if (data.GetBinContent(i) > 0):
			bins += 1
			num = data.GetBinContent(i) - bkg.GetBinContent(i)
			den = data.GetBinError(i)
			chi_square += pow(num/den,2)
	if chi_square == 0:
		print("Chi2 Error")
		return 0
	else:
		print("Chi2 = "+str(chi_square))
		print("NDF = "+str(bins - nInd))
#		print("Chi2/NDF = "+str(chi_square/(bins - nInd)))
		return chi_square, (bins-nInd), bins

def drawDiagnostic(name, ifile1, mass):
	gROOT.SetBatch(True)
	Sigs = [["/home/akobert/CMSSW_11_1_0_pre7/src/RData/2016/NanoTool_UL/M25_UL_nano_merged.root"], "4.148*7639.0/567896.0", "", "sig25_nano", "25 GeV Signal", "7639.0"]
	sig_name = "sig25_nano"
	#ptBins = [200, 230, 255, 290, 360, 1000]
    	ptBins = {} # Dictonary of bins for each year
    	ptBins['2016'] = np.array([200, 225, 245, 275, 325, 700])
    	ptBins['2017'] = np.array([220, 245, 270, 300, 350, 700])
    	ptBins['2018'] = np.array([120, 135, 150, 170, 200, 700])
	mBins = MakeNBinsFromMinToMax(40,0,200)
	VAR = ["jpt", ptBins, "Jet pT (GeV)", "sdm", mBins, "Softdrop Mass (GeV)"]

	Blind = True
	LUMI = {}
	LUMI['2016'] = 3.631
	LUMI['2017'] = 4.148
	LUMI['2018'] = 5.982
	LUMI['Run2'] = 13.761 # 10% of Luminosity
	cmsextra = "Preliminary"

	b_pull_count = TH1F("b_pull_count", "Passing Pull Count", 100, -5, 5)
	s_pull_count = TH1F("s_pull_count", "Sig Passing Pull Count", 100, -5, 5)
	year_pull_count = TH1F("year_pull_count", "Yearly Passing Pull Count", 100, -5, 5)
	b_pull_count_total = TH1F("b_pull_count_total", "Background Only Passing Pull Count All pT Bins", 100, -5, 5)
	s_pull_count_total = TH1F("s_pull_count_total", str(mass)+" GeV Passing Pull Count All pT Bins", 100, -5, 5)

	pre_fail_data_total = TH1F("pre_fail_data_total", "Prefit Failing Data All pT Bins", 40, 0, 200)
	pre_pass_data_total = TH1F("pre_pass_data_total", "Prefit Passing Data All pT Bins", 40, 0, 200)
	b_fail_data_total = TH1F("b_fail_data_total", "Background Only Failing Data All pT Bins", 40, 0, 200)
	b_pass_data_total = TH1F("b_pass_data_total", "Background Only Passing Data All pT Bins", 40, 0, 200)
	s_fail_data_total = TH1F("s_fail_data_total", "Failing Data All pT Bins", 40, 0, 200)
	s_pass_data_total = TH1F("s_pass_data_total", "Passing Data All pT Bins", 40, 0, 200)
	pre_fail_nres_total = TH1F("pre_fail_nres_total", "Prefit Failing Non-Resonant Background All pT Bins", 40, 0, 200)
	pre_pass_nres_total = TH1F("pre_pass_nres_total", "Prefit Passing Non-Resonant Background All pT Bins", 40, 0, 200)
	b_fail_nres_total = TH1F("b_fail_nres_total", "Background Only Failing Non-Resonant Background All pT Bins", 40, 0, 200)
	b_pass_nres_total = TH1F("b_pass_nres_total", "Background Only Passing Non-Resonant Background All pT Bins", 40, 0, 200)
	s_fail_nres_total = TH1F("s_fail_nres_total", str(mass)+" GeV Failing Non-Resonant Background All pT Bins", 40, 0, 200)
	s_pass_nres_total = TH1F("s_pass_nres_total", str(mass)+" GeV Passing Non-Resonant Background All pT Bins", 40, 0, 200)
	pre_fail_ttbar_total = TH1F("pre_fail_ttbar_total", "Prefit Failing TTBar All pT Bins", 40, 0, 200)
	pre_pass_ttbar_total = TH1F("pre_pass_ttbar_total", "Prefit Passing TTBar All pT Bins", 40, 0, 200)
	b_fail_ttbar_total = TH1F("b_fail_ttbar_total", "Background Only Failing TTBar All pT Bins", 40, 0, 200)
	b_pass_ttbar_total = TH1F("b_pass_ttbar_total", "Background Only Passing TTBar All pT Bins", 40, 0, 200)
	s_fail_ttbar_total = TH1F("s_fail_ttbar_total", str(mass)+" GeV Failing TTBar All pT Bins", 40, 0, 200)
	s_pass_ttbar_total = TH1F("s_pass_ttbar_total", str(mass)+" GeV Passing TTBar All pT Bins", 40, 0, 200)
	pre_fail_tbkg_total = TH1F("pre_fail_tbkg_total", "Prefit Failing Total Background All pT Bins", 40, 0, 200)
	pre_pass_tbkg_total = TH1F("pre_pass_tbkg_total", "Prefit Passing Total Background All pT Bins", 40, 0, 200)
	b_fail_tbkg_total = TH1F("b_fail_tbkg_total", "Background Only Failing Total Background All pT Bins", 40, 0, 200)
	b_pass_tbkg_total = TH1F("b_pass_tbkg_total", "Background Only Passing Total Background All pT Bins", 40, 0, 200)
	s_fail_tbkg_total = TH1F("s_fail_tbkg_total", str(mass)+" GeV Failing Total Background All pT Bins", 40, 0, 200)
	s_pass_tbkg_total = TH1F("s_pass_tbkg_total", str(mass)+" GeV Passing Total Background All pT Bins", 40, 0, 200)
	pre_fail_tsig10_total = TH1F("pre_fail_tsig10_total", "10 GeV Prefit Failing Signal All pT Bins", 40, 0, 200)
	pre_pass_tsig10_total = TH1F("pre_pass_tsig10_total", "10 GeV Prefit Passing Signal All pT Bins", 40, 0, 200)
	b_fail_tsig10_total = TH1F("b_fail_tsig10_total", "10 GeV Background Only Failing Signal All pT Bins", 40, 0, 200)
	b_pass_tsig10_total = TH1F("b_pass_tsig10_total", "10 GeV Background Only Passing Signal All pT Bins", 40, 0, 200)
	s_fail_tsig10_total = TH1F("s_fail_tsig10_total", "10 GeV Failing Signal All pT Bins", 40, 0, 200)
	s_pass_tsig10_total = TH1F("s_pass_tsig10_total", "10 GeV Passing Signal All pT Bins", 40, 0, 200)
	pre_fail_tsig25_total = TH1F("pre_fail_tsig25_total", "25 GeV Prefit Failing Signal All pT Bins", 40, 0, 200)
	pre_pass_tsig25_total = TH1F("pre_pass_tsig25_total", "25 GeV Prefit Passing Signal All pT Bins", 40, 0, 200)
	b_fail_tsig25_total = TH1F("b_fail_tsig25_total", "25 GeV Background Only Failing Signal All pT Bins", 40, 0, 200)
	b_pass_tsig25_total = TH1F("b_pass_tsig25_total", "25 GeV Background Only Passing Signal All pT Bins", 40, 0, 200)
	s_fail_tsig25_total = TH1F("s_fail_tsig25_total", "25 GeV Failing Signal All pT Bins", 40, 0, 200)
	s_pass_tsig25_total = TH1F("s_pass_tsig25_total", "25 GeV Passing Signal All pT Bins", 40, 0, 200)
	pre_fail_tsig50_total = TH1F("pre_fail_tsig50_total", "50 GeV Prefit Failing Signal All pT Bins", 40, 0, 200)
	pre_pass_tsig50_total = TH1F("pre_pass_tsig50_total", "50 GeV Prefit Passing Signal All pT Bins", 40, 0, 200)
	b_fail_tsig50_total = TH1F("b_fail_tsig50_total", "50 GeV Background Only Failing Signal All pT Bins", 40, 0, 200)
	b_pass_tsig50_total = TH1F("b_pass_tsig50_total", "50 GeV Background Only Passing Signal All pT Bins", 40, 0, 200)
	s_fail_tsig50_total = TH1F("s_fail_tsig50_total", "50 GeV Failing Signal All pT Bins", 40, 0, 200)
	s_pass_tsig50_total = TH1F("s_pass_tsig50_total", "50 GeV Passing Signal All pT Bins", 40, 0, 200)
	pre_fail_tsig75_total = TH1F("pre_fail_tsig75_total", "75 GeV Prefit Failing Signal All pT Bins", 40, 0, 200)
	pre_pass_tsig75_total = TH1F("pre_pass_tsig75_total", "75 GeV Prefit Passing Signal All pT Bins", 40, 0, 200)
	b_fail_tsig75_total = TH1F("b_fail_tsig75_total", "75 GeV Background Only Failing Signal All pT Bins", 40, 0, 200)
	b_pass_tsig75_total = TH1F("b_pass_tsig75_total", "75 GeV Background Only Passing Signal All pT Bins", 40, 0, 200)
	s_fail_tsig75_total = TH1F("s_fail_tsig75_total", "75 GeV Failing Signal All pT Bins", 40, 0, 200)
	s_pass_tsig75_total = TH1F("s_pass_tsig75_total", "75 GeV Passing Signal All pT Bins", 40, 0, 200)
	pre_fail_wgamma_total = TH1F("pre_fail_wgamma_total", "Prefit Failing WGamma All pT Bins", 40, 0, 200)
	pre_pass_wgamma_total = TH1F("pre_pass_wgamma_total", "Prefit Passing WGamma All pT Bins", 40, 0, 200)
	b_fail_wgamma_total = TH1F("b_fail_wgamma_total", "Background Only Failing WGamma All pT Bins", 40, 0, 200)
	b_pass_wgamma_total = TH1F("b_pass_wgamma_total", "Background Only Passing WGamma All pT Bins", 40, 0, 200)
	s_fail_wgamma_total = TH1F("s_fail_wgamma_total", str(mass)+" GeV Failing WGamma All pT Bins", 40, 0, 200)
	s_pass_wgamma_total = TH1F("s_pass_wgamma_total", str(mass)+" GeV Passing WGamma All pT Bins", 40, 0, 200)
	pre_fail_zgamma_total = TH1F("pre_fail_zgamma_total", "Prefit Failing ZGamma All pT Bins", 40, 0, 200)
	pre_pass_zgamma_total = TH1F("pre_pass_zgamma_total", "Prefit Passing ZGamma All pT Bins", 40, 0, 200)
	b_fail_zgamma_total = TH1F("b_fail_zgamma_total", "Background Only Failing ZGamma All pT Bins", 40, 0, 200)
	b_pass_zgamma_total = TH1F("b_pass_zgamma_total", "Background Only Passing ZGamma All pT Bins", 40, 0, 200)
	s_fail_zgamma_total = TH1F("s_fail_zgamma_total", str(mass)+" GeV Failing ZGamma All pT Bins", 40, 0, 200)
	s_pass_zgamma_total = TH1F("s_pass_zgamma_total", str(mass)+" GeV Passing ZGamma All pT Bins", 40, 0, 200)


	ofile = ROOT.TFile("./fit_b_pull_count.root", "RECREATE")
	ofile2 = ROOT.TFile("./diagnositcs.root", "RECREATE")

	CbnF2_10 = ROOT.TFile.Open("sig10_low/higgsCombineTest.AsymptoticLimits.mH120.root", 'read')
	CbnF2_25 = ROOT.TFile.Open("sig25_low/higgsCombineTest.AsymptoticLimits.mH120.root", 'read')
	CbnF2_50 = ROOT.TFile.Open("sig50_low/higgsCombineTest.AsymptoticLimits.mH120.root", 'read')
	CbnF2_75 = ROOT.TFile.Open("sig75_low/higgsCombineTest.AsymptoticLimits.mH120.root", 'read')
	
	tree10 = CbnF2_10.Get("limit")
	limits10 = []
	for quantile in tree10:
		limits10.append(tree10.limit)
	CL95_10 = limits10[2]  # 95% confidence r value
	print("10 GeV 95% CL r value: "+str(CL95_10))

	tree25 = CbnF2_25.Get("limit")
	limits25 = []
	for quantile in tree25:
		limits25.append(tree25.limit)
	CL95_25 = limits25[2]  # 95% confidence r value
	print("25 GeV 95% CL r value: "+str(CL95_25))

	tree50 = CbnF2_50.Get("limit")
	limits50 = []
	for quantile in tree50:
		limits50.append(tree50.limit)
	CL95_50 = limits50[2]  # 95% confidence r value
	print("50 GeV 95% CL r value: "+str(CL95_50))

	tree75 = CbnF2_75.Get("limit")
	limits75 = []
	for quantile in tree75:
		limits75.append(tree75.limit)
	CL95_75 = limits75[2]  # 95% confidence r value
	print("75 GeV 95% CL r value: "+str(CL95_75))

	years = ['2016', '2017', '2018']
	#years = ['2016', '2017']
	
	# Individual Year Summary
	year_summary = {}
	gof_summary = {}
	samps = {'tbkg', 'ttbar', 'wgamma', 'zgamma', 'nres', 'data', 'sig10', 'sig25', 'sig50', 'sig75'}
	for y in years:
		year_summary[y] = {}
		gof_summary[y] = {}
		gof_summary[y]['bins'] = 0
		gof_summary[y]['chi2'] = 0
		for s in samps:
			year_summary[y][s] = TH1F("s_pass_"+y+"_"+s+"_total", y+" Passing "+s+" All pT Bins", 40, 0, 200) 


	# GOF Running Totals
	b_bins = 0
	s_bins = 0

	b_chi = 0
	s_chi = 0

	

	for y in years:
		for pf in ["pass", "fail"]:
			for a in range(0, 5):
				print(y+" ptbin #"+str(a))
				NAME = "bin_"+y+"_ptbin"+str(a)+pf
				refF = ROOT.TFile("FitHist_"+y+".root")
				refH = refF.Get("Data_"+pf+"_soft")
				CbnF = ROOT.TFile("sig25_low/fitDiagnosticsTest.root") # Get non-signal templates from this one
				CbnF_10 = ROOT.TFile("sig10_low/fitDiagnosticsTest.root")
				CbnF_25 = ROOT.TFile("sig25_low/fitDiagnosticsTest.root")
				CbnF_50 = ROOT.TFile("sig50_low/fitDiagnosticsTest.root")
				CbnF_75 = ROOT.TFile("sig75_low/fitDiagnosticsTest.root")
	
				run_code = 0 # Distinguish between prefit (1), background only (2)
				for P in ["prefit", "fit_b"]: 
	#			for P in ["fit_b", "fit_s"]:
	#			for P in ["prefit", "fit_b"]:
						run_code += 1
	
						cDATA = CbnF.Get("shapes_"+P+"/"+NAME+"/data")
						cNRES = CbnF.Get("shapes_"+P+"/"+NAME+"/NonRes"+y)	#Non-Resonant Background
						cTT = CbnF.Get("shapes_"+P+"/"+NAME+"/TTBar")
					#	cTT = cTT.Clone()
						cTBKG = CbnF.Get("shapes_"+P+"/"+NAME+"/total_background")
						if run_code != 2:
							cSIG10 = CbnF_10.Get("shapes_"+P+"/"+NAME+"/total_signal")
							cSIG25 = CbnF_25.Get("shapes_"+P+"/"+NAME+"/total_signal")
							cSIG50 = CbnF_50.Get("shapes_"+P+"/"+NAME+"/total_signal")
							cSIG75 = CbnF_75.Get("shapes_"+P+"/"+NAME+"/total_signal")
						else:
							cSIG10 = CbnF_10.Get("shapes_prefit/"+NAME+"/total_signal")	# Take pre-fit signal to scale
							cSIG25 = CbnF_25.Get("shapes_prefit/"+NAME+"/total_signal")	# Take pre-fit signal to scale
							cSIG50 = CbnF_50.Get("shapes_prefit/"+NAME+"/total_signal")	# Take pre-fit signal to scale
							cSIG75 = CbnF_75.Get("shapes_prefit/"+NAME+"/total_signal")	# Take pre-fit signal to scale
	                                        cWG = CbnF.Get("shapes_"+P+"/"+NAME+"/WGamma")
	                                        cZG = CbnF.Get("shapes_"+P+"/"+NAME+"/ZGamma")
	
						#Testing
	#					print(cDATA.GetN())
#						print(cTT.GetNbinsX())
#						print(cTBKG.GetNbinsX())
#						print(cSIG.GetNbinsX())
	
						cDATA = convertAsymGraph(cDATA, cTT, "data"+P)
						Hvec = []
	
						if P == "prefit":
							for i in [cDATA, cTT, cTBKG, cSIG25, cSIG50, cSIG75, cWG, cZG, cNRES, cSIG10]:
								i.Scale(5.0)
								Hvec.append(convertBinNHist(i, refH, i.GetName()+"new"+P))
						else:
							for i in [cDATA, cTT, cTBKG, cSIG25, cSIG50, cSIG75, cWG, cZG, cNRES, cSIG10]:
								i.Scale(5.0)
								Hvec.append(convertBinNHist(i, refH, i.GetName()+"new"+P))
						if run_code == 2:
							Hvec[9].Scale(1.0/5.0) # Fixing double scaling of cSIG
							Hvec[3].Scale(1.0/5.0) # Fixing double scaling of cSIG
							Hvec[4].Scale(1.0/5.0) # Fixing double scaling of cSIG
							Hvec[5].Scale(1.0/5.0) # Fixing double scaling of cSIG
							print("Scaling Signal By CL95")
							Hvec[3].Scale(CL95_25)
							Hvec[4].Scale(CL95_50)
							Hvec[5].Scale(CL95_75)
	
						if P == "prefit" and pf == "fail":
							pre_fail_data_total.Add(Hvec[0])
							pre_fail_ttbar_total.Add(Hvec[1])
							pre_fail_tbkg_total.Add(Hvec[2])
							pre_fail_tsig10_total.Add(Hvec[9])
							pre_fail_tsig25_total.Add(Hvec[3])
							pre_fail_tsig50_total.Add(Hvec[4])
							pre_fail_tsig75_total.Add(Hvec[5])
							pre_fail_wgamma_total.Add(Hvec[6])
							pre_fail_zgamma_total.Add(Hvec[7])
							pre_fail_nres_total.Add(Hvec[8])
							
						elif P == "prefit" and pf == "pass":
							pre_pass_data_total.Add(Hvec[0])
							pre_pass_ttbar_total.Add(Hvec[1])
							pre_pass_tbkg_total.Add(Hvec[2])
							pre_pass_tsig10_total.Add(Hvec[9])
							pre_pass_tsig25_total.Add(Hvec[3])
							pre_pass_tsig50_total.Add(Hvec[4])
							pre_pass_tsig75_total.Add(Hvec[5])
							pre_pass_wgamma_total.Add(Hvec[6])
							pre_pass_zgamma_total.Add(Hvec[7])
							pre_pass_nres_total.Add(Hvec[8])
						
						elif P == "fit_b" and pf == "fail":
							b_fail_data_total.Add(Hvec[0])
							b_fail_ttbar_total.Add(Hvec[1])
							b_fail_tbkg_total.Add(Hvec[2])
							b_fail_tsig10_total.Add(Hvec[9])
							b_fail_tsig25_total.Add(Hvec[3])
							b_fail_tsig50_total.Add(Hvec[4])
							b_fail_tsig75_total.Add(Hvec[5])
							b_fail_wgamma_total.Add(Hvec[6])
							b_fail_zgamma_total.Add(Hvec[7])
							b_fail_nres_total.Add(Hvec[8])
							s_fail_data_total.Add(Hvec[0])
							s_fail_ttbar_total.Add(Hvec[1])
							s_fail_tbkg_total.Add(Hvec[2])
							s_fail_tsig10_total.Add(Hvec[9])
							s_fail_tsig25_total.Add(Hvec[3])
							s_fail_tsig50_total.Add(Hvec[4])
							s_fail_tsig75_total.Add(Hvec[5])
							s_fail_wgamma_total.Add(Hvec[6])
							s_fail_zgamma_total.Add(Hvec[7])
							s_fail_nres_total.Add(Hvec[8])
						
						elif P == "fit_b" and pf == "pass":
							b_pass_data_total.Add(Hvec[0])
							b_pass_ttbar_total.Add(Hvec[1])
							b_pass_tbkg_total.Add(Hvec[2])
							b_pass_tsig10_total.Add(Hvec[9])
							b_pass_tsig25_total.Add(Hvec[3])
							b_pass_tsig50_total.Add(Hvec[4])
							b_pass_tsig75_total.Add(Hvec[5])
							b_pass_wgamma_total.Add(Hvec[6])
							b_pass_zgamma_total.Add(Hvec[7])
							b_pass_nres_total.Add(Hvec[8])
							s_pass_data_total.Add(Hvec[0])
							s_pass_ttbar_total.Add(Hvec[1])
							s_pass_tbkg_total.Add(Hvec[2])
							s_pass_tsig10_total.Add(Hvec[9])
							s_pass_tsig25_total.Add(Hvec[3])
							s_pass_tsig50_total.Add(Hvec[4])
							s_pass_tsig75_total.Add(Hvec[5])
							s_pass_wgamma_total.Add(Hvec[6])
							s_pass_zgamma_total.Add(Hvec[7])
							s_pass_nres_total.Add(Hvec[8])
							
						
	
			#			cDATA.SetLineColor(kBlack)
			#			cTT.SetLineColor(kRed)
			#			c1 = ROOT.TCanvas()
			#			c1.cd()
						
			#			cDATA.Draw("hist")
			#			cTT.Draw("same hist")
			#			cTBKG.Draw("same hist")
			#			L1 = ROOT.TLegend(0.48,0.6,0.86,0.86)
			#			L1.SetFillColor(0)
			#			L1.SetLineColor(0)
			#			L1.AddEntry(cDATA, "data", "PE")
			#			L1.AddEntry(cTBKG, "total background", "L")
			#			L1.AddEntry(cTT, "t#bar{t} component", "L")
			#			if P != "fit_b":
			#				cSIG.SetLineColor(kGreen+1)
			#				cSIG.Draw("same hist")
			#				L1.AddEntry(cSIG, sig25_name, "F")
			#			c1.SaveAs("./"+P+"_data_test.png")
			#			c1.Close()
						
	#					Pull2D = Hvec[0][2].Clone("pull"+P)
						
	#					Cd = ROOT.TCanvas()
	#					Cd.cd()
	#					Pull2D.Draw("colz")
	#					Cd.Print("results/"+NAME+"/Data2D_"+P+".png")
						
	#					Pull2D.Add(Hvec[2][2], -1.)
	#					Pull2D.Divide(Hvec[0][3])
	#					Hvec[2][3].Divide(Hvec[0][3])
	#					Pull2D.GetZaxis().SetRangeUser(-3.,3.)
						
						
#						C2 = ROOT.TCanvas()
#						C2.cd()
#						Pull2D.Draw("colz")
#						C2.Print("results/"+NAME+"/Pull2D_"+P+".png")
						
#						C2e = ROOT.TCanvas()
#						C2e.cd()
#						Hvec[3][2].Draw("colz")
#						C2e.Print("results/"+NAME+"/Sig2D_"+P+".png")
						savename = P+"_sdm_"+str(ptBins[y][a])+"-"+str(ptBins[y][a+1])
						if P == "fit_b":
							savename10 = "sig10_sdm_"+str(ptBins[y][a])+"-"+str(ptBins[y][a+1])
							savename25 = "sig25_sdm_"+str(ptBins[y][a])+"-"+str(ptBins[y][a+1])
							savename50 = "sig50_sdm_"+str(ptBins[y][a])+"-"+str(ptBins[y][a+1])
							savename75 = "sig75_sdm_"+str(ptBins[y][a])+"-"+str(ptBins[y][a+1])
						ptbin = str(ptBins[y][a])+"-"+str(ptBins[y][a+1])
						#data = DBBW(Hvec[0])
						#GoodPlotFormat(data, "markers", ROOT.kBlack, 20)
						#bkg = DBBW(Hvec[2])
						#GoodPlotFormat(bkg,"thickline", ROOT.kBlue, 1)
						#GoodPlotFormat(data,"thickline", ROOT.kBlue, 1)
						#ttbar = DBBW(Hvec[1])
						#GoodPlotFormat(ttbar,"thickline", ROOT.kRed, 1)
						#sig25 = DBBW(Hvec[3])
						#GoodPlotFormat(sig25, "fill", ROOT.kGreen+1, 3003)
	                                        #s75 = DBBW(Hvec[6])
       	                                 #GoodPlotFormat(s75,"thickline", ROOT.kViolet, 1)
       	                                 #s100 = DBBW(Hvec[7])
       	                                 #GoodPlotFormat(s100,"thickline", ROOT.kOrange, 1)
						data = Hvec[0]
						GoodPlotFormat(data, "markers", ROOT.kBlack, 20)
						bkg = Hvec[2]
						GoodPlotFormat(bkg,"thickline", ROOT.kBlue, 1)
						#GoodPlotFormat(data,"thickline", ROOT.kBlue, 1)
						ttbar = Hvec[1]
						GoodPlotFormat(ttbar,"fill", ROOT.kRed, ROOT.kRed, 1001)
						sig10 = Hvec[9]
						GoodPlotFormat(sig10, "fill", ROOT.kBlue+1, ROOT.kBlue+1, 3003)
						sig25 = Hvec[3]
						GoodPlotFormat(sig25, "fill", ROOT.kGreen+1, ROOT.kGreen+1, 3003)
						sig50 = Hvec[4]
						GoodPlotFormat(sig50, "fill", ROOT.kGray+1, ROOT.kGray+1, 3003)
						sig75 = Hvec[5]
						GoodPlotFormat(sig75, "fill", ROOT.kYellow+1, ROOT.kYellow+1, 3003)
	                                        WG = Hvec[6]
	                                        GoodPlotFormat(WG,"fill", ROOT.kViolet, ROOT.kViolet, 1001)
	                                        ZG = Hvec[7]
						GoodPlotFormat(ZG,"fill", ROOT.kOrange, ROOT.kOrange, 1001)
						nonRes = Hvec[8]
						GoodPlotFormat(nonRes,"fill", ROOT.kCyan, 10, 1)
	
						cheapline = data.Clone("cheapline")
						cheapline.Add(data,-1.)
						cheapline.GetYaxis().SetTitle("#frac{data - bkg}{#sigma_{data}}")
						cheapline.GetYaxis().SetTitleSize(0.175);
						cheapline.GetYaxis().SetNdivisions(6);
						cheapline.GetYaxis().SetLabelSize(0.145);
						cheapline.GetYaxis().SetTitleOffset(0.225);
						cheapline.GetYaxis().CenterTitle(True)
						cheapline.GetYaxis().SetRangeUser(-5.,5.)
						GoodPlotFormat(cheapline, "thinline", ROOT.kGray, 4)
						bkg.GetYaxis().SetTitle("Events / 5 GeV")
						bkg.GetYaxis().SetTitleOffset(0.5);
						bkg.GetYaxis().SetTitleSize(0.075);
#       	                                 if P == "prefit":
 #      	                                         print("PreFit Test sqrt: "+str(math.sqrt(bkg.GetBinContent(3))))
  #     	                                         print("PreFit Test Error: "+str(bkg.GetBinError(3)))
      # 	                                 if P == "fit_b":
       #	                                         print("fit_b Test sqrt: "+str(math.sqrt(bkg.GetBinContent(3))))
        	#                                        print("fit_b Test Error: "+str(bkg.GetBinError(3)))
	
						FindAndSetMax(data, bkg)
							
						E = []
						EP = []
						for i in range(1,data.GetNbinsX()+1):
							Err = bkg.GetBinError(i)
							
							blX = bkg.GetBinLowEdge(i)
							blY = bkg.GetBinContent(i) - Err
							trX = bkg.GetBinWidth(i) + blX
							trY = bkg.GetBinContent(i) + Err
							tBox = ROOT.TBox(blX,blY,trX,trY)
							if  data.GetBinError(i) > 0:
								ue = Err/data.GetBinError(i)
							else:
								ue = Err/2.7
							ue = min(5.0, ue)
							tPBox = ROOT.TBox(blX, -1*ue, trX, ue)
							tBox.SetFillColor(25)
							tBox.SetFillStyle(3144)
							tPBox.SetFillColor(25)
							tPBox.SetFillStyle(3144)
							tBox.SetLineColor(ROOT.kWhite)
							tPBox.SetLineColor(ROOT.kWhite)
							E.append(tBox)
							EP.append(tPBox)
											
						bkg.GetXaxis().SetLabelSize(0)
						
						pull = data.Clone("pull")
						pull.Add(bkg, -1.)
						GoodPlotFormat(pull, "markers", ROOT.kBlack, 20)
						for i in range(pull.GetNbinsX()):
							if not data.GetBinContent(i+1) == 0:
								pull.SetBinContent(i+1, pull.GetBinContent(i+1)/data.GetBinError(i+1))
								pull.SetBinError(i+1, 1)		
								#Fill Pull Counter
								if pf == "pass" and P == "fit_b":
									b_pull_count.Fill(pull.GetBinContent(i+1))
									b_pull_count_total.Fill(pull.GetBinContent(i+1))
									print(y+" fit_b pTbin #"+str(a)+" bin #"+str(i+1)+" Data content: "+str(data.GetBinContent(i+1)))
									print(y+" fit_b pTbin #"+str(a)+" bin #"+str(i+1)+" Background content: "+str(bkg.GetBinContent(i+1)))
									print(y+" fit_b pTbin #"+str(a)+" bin #"+str(i+1)+" Data Error: "+str(data.GetBinError(i+1)))
									print(y+" fit_b pTbin #"+str(a)+" bin #"+str(i+1)+" Pull content: "+str(pull.GetBinContent(i+1)))
								if pf == "pass" and run_code == 3:
									s_pull_count.Fill(pull.GetBinContent(i+1))
									s_pull_count_total.Fill(pull.GetBinContent(i+1)) # Same Pull Count For fit_b and sig25 (double check this)
									print("sig25 pTbin #"+str(a)+" bin #"+str(i+1)+" Signal value: "+str(sig25.GetBinContent(i+1)))
									print("sig25 pTbin #"+str(a)+" bin #"+str(i+1)+" Data content: "+str(data.GetBinContent(i+1)))
									print("sig25 pTbin #"+str(a)+" bin #"+str(i+1)+" Background content: "+str(bkg.GetBinContent(i+1)))
									print("sig25 pTbin #"+str(a)+" bin #"+str(i+1)+" Data Error: "+str(data.GetBinError(i+1)))
									print("sig25 pTbin #"+str(a)+" bin #"+str(i+1)+" Pull content: "+str(pull.GetBinContent(i+1)))
							else:
								pull.SetBinContent(i+1, 0)
								pull.SetBinError(i+1, 0)
						
						# Stack backgrounds for plotting
						ZG.Add(nonRes)
						WG.Add(ZG)
						ttbar.Add(WG)


						L = ROOT.TLegend(0.65,0.7,0.92,0.88)
						#L = ROOT.TLegend(0.58,0.65,0.895,0.895)
						L.SetFillColor(0)
						L.SetLineColor(0)
						if not Blind: L.AddEntry(data, "Data", "PE")
						L.AddEntry(bkg, "Total Background", "L")
						L.AddEntry(ttbar, "t#bar{t} component", "L")
	                                        L.AddEntry(WG, "W+Gamma", "L")
	                                        L.AddEntry(ZG, "Z+Gamma", "L")
						L.AddEntry(nonRes, "Non-Resonant Background", "L")
						L.AddEntry(E[0], "Background Uncertainty", "F")
						if run_code == 2:  L.AddEntry(sig10, "10 GeV Signal", "F")
						if run_code == 2:  L.AddEntry(sig25, "25 GeV Signal", "F")
						if run_code == 2:  L.AddEntry(sig50, "50 GeV Signal", "F")
						if run_code == 2:  L.AddEntry(sig75, "75 GeV Signal", "F")
						C = ROOT.TCanvas()
						C.cd()
						p12 = ROOT.TPad("pad1", "tall",0,0.165,1,1)
						p22 = ROOT.TPad("pad2", "short",0,0.0,1.0,0.23)
						p22.SetBottomMargin(0.35)
						p12.Draw()
						p22.Draw()
						p12.cd() # top
						ROOT.gPad.SetTicks(1,1)
						if P == "prefit" and pf == "pass":
							bkg.SetTitle(str(y)+" Prefit Passing Softdrop in "+str(ptBins[y][a])+"-"+str(ptBins[y][a+1])+" pT Bin")
						if P == "prefit" and pf == "fail":
							bkg.SetTitle(str(y)+" Prefit Failing Softdrop in "+str(ptBins[y][a])+"-"+str(ptBins[y][a+1])+" pT Bin")
						if P == "fit_b" and pf == "pass":
							bkg.SetTitle(str(y)+" Background Only Passing Softdrop in "+str(ptBins[y][a])+"-"+str(ptBins[y][a+1])+" pT Bin")
						if P == "fit_b" and pf == "fail":
							bkg.SetTitle(str(y)+" Background Only Failing Softdrop in "+str(ptBins[y][a])+"-"+str(ptBins[y][a+1])+" pT Bin")
						bkg.Draw("hist")
						ttbar.Draw("histsame")
	                                        WG.Draw("histsame")
	                                        ZG.Draw("histsame")
						nonRes.Draw("histsame")
						if run_code == 2: sig10.Draw("histsame")
						if run_code == 2: sig25.Draw("histsame")
						if run_code == 2: sig50.Draw("histsame")
						if run_code == 2: sig75.Draw("histsame")
						for e in E: e.Draw("same")
						data.Draw("esame")
						L.Draw("same")
						ROOT.TGaxis.SetMaxDigits(3)
						p12.RedrawAxis()
						AddCMSLumi(ROOT.gPad, LUMI[y], cmsextra)
						pTbin(ROOT.gPad, ptbin)

						# Add Chi2 GOF
						chi2, ndf, bins = Chi2(data, bkg, 9) #9 fit parameters for 2,2 Bernstein
						if P == "fit_b" and pf == "pass":
							b_bins += bins
							b_chi += chi2
						#GOF(ROOT.gPad, chi2, ndf)

						p22.cd() # bottom
						ROOT.gPad.SetTicks(1,1)
						cheapline.GetXaxis().SetTitleSize(0.1925);
						cheapline.GetXaxis().SetLabelSize(0.16);
						cheapline.GetXaxis().SetTitleOffset(0.84);
						cheapline.Draw("hist")
						for e in EP: e.Draw("same")
	
						l1 = TLine(0,2,200,2)
						l1.SetLineColor(kGray)
						l1.SetLineWidth(1)
						l1.SetLineStyle(7)
						l2 = TLine(0,-2,200,-2)
						l2.SetLineColor(kGray)
						l2.SetLineWidth(1)
						l2.SetLineStyle(7)
						pull.Draw("esame")
						l1.Draw("same")
						l2.Draw("same")
						p22.RedrawAxis()
	#					C.Print("results/"+NAME+"/"+savename+"_"+pf+str(a)+"_v5.root")
						C.Print("results/"+NAME+"/"+y+"_"+savename+"_"+pf+str(a)+"_"+str(name)+".png")
						p12.cd()
						gPad.SetLogy()
						p12.Draw()
						C.Print("results/"+NAME+"/"+y+"_"+savename+"_"+pf+str(a)+"_"+str(name)+"_logy.png")
						if pf == "pass" and P == "fit_b":
							b_pull_count.SetTitle("fit_b Passing Pull Count "+str(ptBins[y][a])+"-"+str(ptBins[y][a+1])+" GeV")
							pull_name = P+"_"+pf+"_pull_count_bin"+str(a)
							ofile.WriteObject(b_pull_count, pull_name)
							b_pull_count.Reset()
	
						# Only for fit_b make second plot for signal:
						if P == "fit_b":
							L = ROOT.TLegend(0.65,0.62,0.92,0.88)
							#L = ROOT.TLegend(0.62,0.7,0.92,0.86)
							#L = ROOT.TLegend(0.58,0.65,0.895,0.895)
							L.SetFillColor(0)
							L.SetLineColor(0)
							if not Blind: L.AddEntry(data, "Data", "PE")
							L.AddEntry(bkg, "Total Background", "L")
							L.AddEntry(ttbar, "t#bar{t} component", "L")
		                                        L.AddEntry(WG, "W+Gamma", "L")
		                                        L.AddEntry(ZG, "Z+Gamma", "L")
							L.AddEntry(nonRes, "Non-Resonant Background", "L")
							L.AddEntry(E[0], "Background Uncertainty", "F")
							L.AddEntry(sig10, "10 GeV Signal", "F")
							L.AddEntry(sig25, "25 GeV Signal", "F")
							L.AddEntry(sig50, "50 GeV Signal", "F")
							L.AddEntry(sig75, "75 GeV Signal", "F")
							C = ROOT.TCanvas()
							C.cd()
							p12 = ROOT.TPad("pad1", "tall",0,0.165,1,1)
							p22 = ROOT.TPad("pad2", "short",0,0.0,1.0,0.23)
							p22.SetBottomMargin(0.35)
							p12.Draw()
							p22.Draw()
							p12.cd() # top
							ROOT.gPad.SetTicks(1,1)
							if P == "fit_b" and pf == "pass":
								bkg.SetTitle(str(y)+" Background with Signal Passing Softdrop in "+str(ptBins[y][a])+"-"+str(ptBins[y][a+1])+" pT Bin")
							if P == "fit_b" and pf == "fail":
								bkg.SetTitle(str(y)+" Background with Signal Failing Softdrop in "+str(ptBins[y][a])+"-"+str(ptBins[y][a+1])+" pT Bin")
							bkg.Draw("hist")
							ttbar.Draw("histsame")
							WG.Draw("histsame")
	       		                                ZG.Draw("histsame")
							nonRes.Draw("histsame")
							sig10.Draw("histsame")
							sig25.Draw("histsame")
							sig50.Draw("histsame")
							sig75.Draw("histsame")
							for e in E: e.Draw("same")
							data.Draw("esame")
							L.Draw("same")
							ROOT.TGaxis.SetMaxDigits(3)
							p12.RedrawAxis()
							AddCMSLumi(ROOT.gPad, LUMI[y], cmsextra)
							pTbin(ROOT.gPad, ptbin)

							# Add Chi2 GOF
							chi2, ndf, bins = Chi2(data, bkg, 9) #9 fit parameters for 2,2 Bernstein
							if pf == "pass":
								s_bins += bins
								s_chi += chi2
								gof_summary[y]['bins'] += bins
								gof_summary[y]['chi2'] += chi2

								year_summary[y]['ttbar'].Add(ttbar)
								year_summary[y]['wgamma'].Add(WG)
								year_summary[y]['zgamma'].Add(ZG)
								year_summary[y]['nres'].Add(nonRes)
								year_summary[y]['tbkg'].Add(bkg)
								year_summary[y]['data'].Add(data)
								year_summary[y]['sig10'].Add(sig10)
								year_summary[y]['sig25'].Add(sig25)
								year_summary[y]['sig50'].Add(sig50)
								year_summary[y]['sig75'].Add(sig75)

							#GOF(ROOT.gPad, chi2, ndf)

							p22.cd() # bottom
							ROOT.gPad.SetTicks(1,1)
							cheapline.GetXaxis().SetTitleSize(0.1925);
							cheapline.GetXaxis().SetLabelSize(0.16);
							cheapline.GetXaxis().SetTitleOffset(0.84);
							cheapline.Draw("hist")
							for e in EP: e.Draw("same")
		
							l1 = TLine(0,2,200,2)
							l1.SetLineColor(kGray)
							l1.SetLineWidth(1)
							l1.SetLineStyle(7)
							l2 = TLine(0,-2,200,-2)
							l2.SetLineColor(kGray)
							l2.SetLineWidth(1)
							l2.SetLineStyle(7)
							pull.Draw("esame")
							l1.Draw("same")
							l2.Draw("same")
							p22.RedrawAxis()
							C.Print("results/"+NAME+"/"+y+"_multiSig_"+pf+str(a)+"_"+str(name)+".png")
							p12.cd()
							gPad.SetLogy()
							p12.Draw()
							C.Print("results/"+NAME+"/"+y+"_multiSig_"+pf+str(a)+"_"+str(name)+"_logy.png")
							if pf == "pass":
								s_pull_count.SetTitle("multiSig Passing Pull Count "+str(ptBins[y][a])+"-"+str(ptBins[y][a+1])+" GeV")
								pull_name = "multiSig_"+pf+"_pull_count_bin"+str(a)
								ofile.WriteObject(s_pull_count, pull_name)
								s_pull_count.Reset()
	
	

		# Year Summaries
		GoodPlotFormat(year_summary[y]['data'], "markers", ROOT.kBlack, 20)
		GoodPlotFormat(year_summary[y]['tbkg'],"thickline", ROOT.kBlue, 1)
		GoodPlotFormat(year_summary[y]['ttbar'], "fill", ROOT.kRed, ROOT.kRed, 1001)
        	GoodPlotFormat(year_summary[y]['wgamma'],"fill", ROOT.kViolet, ROOT.kViolet, 1001)
        	GoodPlotFormat(year_summary[y]['zgamma'],"fill", ROOT.kOrange, ROOT.kOrange, 1001)
        	GoodPlotFormat(year_summary[y]['nres'],"fill", ROOT.kCyan, 10, 1)
		GoodPlotFormat(year_summary[y]['sig10'], "fill", ROOT.kBlue+1, ROOT.kBlue+1, 3003)
		GoodPlotFormat(year_summary[y]['sig25'], "fill", ROOT.kGreen+1, ROOT.kGreen+1, 3003)
		GoodPlotFormat(year_summary[y]['sig50'], "fill", ROOT.kGray+1, ROOT.kGray+1, 3003)
		GoodPlotFormat(year_summary[y]['sig75'], "fill", ROOT.kYellow+1, ROOT.kYellow+1, 3003)


		year_pull_count.Reset()
		year_cheapline = year_summary[y]['data'].Clone("year_cheapline")
		year_cheapline.Add(year_summary[y]['data'],-1.)
		year_cheapline.SetTitle("")
		year_cheapline.GetYaxis().SetTitle("#frac{data - bkg}{#sigma_{data}}")
		year_cheapline.GetYaxis().SetTitleSize(0.175);
		year_cheapline.GetYaxis().SetNdivisions(6);
		year_cheapline.GetYaxis().SetLabelSize(0.145);
		year_cheapline.GetYaxis().SetTitleOffset(0.225);
		year_cheapline.GetYaxis().CenterTitle(True)
		year_cheapline.GetYaxis().SetRangeUser(-5.,5.)
		GoodPlotFormat(year_cheapline, "thinline", ROOT.kGray, 4)
		year_summary[y]['tbkg'].GetYaxis().SetTitle("Events / 5 GeV")
		year_summary[y]['tbkg'].GetYaxis().SetTitleOffset(0.5);
		year_summary[y]['tbkg'].GetYaxis().SetTitleSize(0.075);
		FindAndSetMax(year_summary[y]['data'], year_summary[y]['tbkg'])
		
							
		E = []
		EP = []
		for i in range(1,year_summary[y]['data'].GetNbinsX()+1):
			Err = year_summary[y]['tbkg'].GetBinError(i)
					
			blX = year_summary[y]['tbkg'].GetBinLowEdge(i)
			blY = year_summary[y]['tbkg'].GetBinContent(i) - Err
			trX = year_summary[y]['tbkg'].GetBinWidth(i) + blX
			trY = year_summary[y]['tbkg'].GetBinContent(i) + Err
			tBox = ROOT.TBox(blX,blY,trX,trY)
			if  year_summary[y]['data'].GetBinError(i) > 0:
				ue = Err/year_summary[y]['data'].GetBinError(i)
			else:
				ue = Err/2.7
			ue = min(5.0, ue)
			tPBox = ROOT.TBox(blX, -1*ue, trX, ue)
			tBox.SetFillColor(25)
			tBox.SetFillStyle(3144)
			tPBox.SetFillColor(25)
			tPBox.SetFillStyle(3144)
			tBox.SetLineColor(ROOT.kWhite)
			tPBox.SetLineColor(ROOT.kWhite)
			E.append(tBox)
			EP.append(tPBox)
											
		year_summary[y]['tbkg'].GetXaxis().SetLabelSize(0)
				
		pull = year_summary[y]['data'].Clone("pull")
		pull.Add(year_summary[y]['tbkg'], -1.)
		GoodPlotFormat(pull, "markers", ROOT.kBlack, 20)
		for i in range(pull.GetNbinsX()):
			if not year_summary[y]['data'].GetBinContent(i+1) == 0:
				pull.SetBinContent(i+1, pull.GetBinContent(i+1)/year_summary[y]['data'].GetBinError(i+1))
				pull.SetBinError(i+1, 1)		
				#Fill Pull Counter
				year_pull_count.Fill(pull.GetBinContent(i+1))
			else:
				pull.SetBinContent(i+1, 0)
				pull.SetBinError(i+1, 0)
						
		# Stack backgrounds for plotting
#		year_summary[y]['zgamma'].Add(year_summary[y]['nres'])
#		year_summary[y]['wgamma'].Add(year_summary[y]['zgamma'])
#		year_summary[y]['ttbar'].Add(year_summary[y]['wgamma'])
	
		L = ROOT.TLegend(0.62,0.54,0.90,0.88)
		#L = ROOT.TLegend(0.6,0.55,0.9,0.86)
		#L = ROOT.TLegend(0.58,0.65,0.895,0.895)
		L.SetFillColor(0)
		L.SetLineColor(0)
		if not Blind: L.AddEntry(year_summary[y]['data'], "Data", "PE")
		L.AddEntry(year_summary[y]['tbkg'], "Total background", "L")
		L.AddEntry(year_summary[y]['ttbar'], "t#bar{t} component", "F")
		L.AddEntry(year_summary[y]['wgamma'], "W+Gamma", "F")
		L.AddEntry(year_summary[y]['zgamma'], "Z+Gamma", "F")
		L.AddEntry(year_summary[y]['nres'], "Non-Resonant Background", "L")
		L.AddEntry(year_summary[y]['sig10'], "10 GeV Z' Signal", "F")
		L.AddEntry(year_summary[y]['sig25'], "25 GeV Z' Signal", "F")
		L.AddEntry(year_summary[y]['sig50'], "50 GeV Z' Signal", "F")
		L.AddEntry(year_summary[y]['sig75'], "75 GeV Z' Signal", "F")
		L.AddEntry(E[0], "Background Uncertainty", "F")
#		if P != "fit_b":  L.AddEntry(year_summary[y]['sig25'], year_summary[y]['sig25']_name, "F")
		C = ROOT.TCanvas()
		C.cd()
		p12 = ROOT.TPad("pad1", "tall",0,0.165,1,1)
		p22 = ROOT.TPad("pad2", "short",0,0.0,1.0,0.23)
		p22.SetBottomMargin(0.35)
		p12.Draw()
		p22.Draw()
		p12.cd() # top
		ROOT.gPad.SetTicks(1,1)
		year_summary[y]['tbkg'].SetTitle(y+" (10%) Multiple Signal Background Fit Plot With Signal All Bins")
		year_summary[y]['tbkg'].Draw("hist")
		year_summary[y]['ttbar'].Draw("histsame")
	        year_summary[y]['wgamma'].Draw("histsame")
	        year_summary[y]['zgamma'].Draw("histsame")
		year_summary[y]['nres'].Draw("histsame")
		year_summary[y]['sig10'].Draw("histsame")
		year_summary[y]['sig25'].Draw("histsame")
		year_summary[y]['sig50'].Draw("histsame")
		year_summary[y]['sig75'].Draw("histsame")
	#	if P != "fit_b": year_summary[y]['sig25'].Draw("histsame")
		for e in E: e.Draw("same")
		year_summary[y]['data'].Draw("esame")
		L.Draw("same")
		ROOT.TGaxis.SetMaxDigits(3)
		p12.RedrawAxis()
		AddCMSLumi(ROOT.gPad, LUMI[y], cmsextra)
		
		# Add Chi2 GOF
		#chi2, ndf, bins = Chi2(year_summary[y]['data'], year_summary[y]['tbkg'], 9) #9 fit parameters for 2,2 Bernstein
		print("Year: "+y)
		print("Chi2: "+str(gof_summary[y]['chi2']))
		print("Bins: "+str(gof_summary[y]['bins']))
		print("NDF: "+str(gof_summary[y]['bins']-9))
		print("Chi2/NDF: "+str(gof_summary[y]['chi2']/(gof_summary[y]['bins']-9)))
		GOF(ROOT.gPad, gof_summary[y]['chi2'], gof_summary[y]['bins']-9)
	
	#	pTbin(ROOT.gPad, "Full")
		p22.cd() # bottom
		ROOT.gPad.SetTicks(1,1)
		year_cheapline.GetXaxis().SetTitle("Softdrop Mass [GeV]")
		year_cheapline.GetXaxis().SetTitleSize(0.1925);
		year_cheapline.GetXaxis().SetLabelSize(0.16);
		year_cheapline.GetXaxis().SetTitleOffset(0.84);
		year_cheapline.Draw("hist")
		for e in EP: e.Draw("same")
		l1 = TLine(0,2,200,2)
		l1.SetLineColor(kGray)
		l1.SetLineWidth(1)
		l1.SetLineStyle(7)
		l2 = TLine(0,-2,200,-2)
		l2.SetLineColor(kGray)
		l2.SetLineWidth(1)
		l2.SetLineStyle(7)
		pull.Draw("esame")
		l1.Draw("same")
		l2.Draw("same")
		p22.RedrawAxis()
		C.Print("results/multiSig_pass_"+y+"_"+str(name)+".png")
		p12.cd()
		gPad.SetLogy()
		p12.Draw()
		C.Print("results/multiSig_pass_"+y+"_"+str(name)+"_logy.png")
		year_pull_count.SetTitle("Passing Pull Count Total Background")
		year_pull_name = "sig_pas_"+y+"_pull_count_totalbkg"
		ofile.WriteObject(year_pull_count, year_pull_name)	
	
	

	GoodPlotFormat(pre_fail_data_total, "markers", ROOT.kBlack, 20)
	GoodPlotFormat(pre_pass_data_total, "markers", ROOT.kBlack, 20)
	GoodPlotFormat(b_fail_data_total, "markers", ROOT.kBlack, 20)
	GoodPlotFormat(b_pass_data_total, "markers", ROOT.kBlack, 20)
	GoodPlotFormat(s_fail_data_total, "markers", ROOT.kBlack, 20)
	GoodPlotFormat(s_pass_data_total, "markers", ROOT.kBlack, 20)
	GoodPlotFormat(pre_fail_tbkg_total,"thickline", ROOT.kBlue, 1)
	GoodPlotFormat(pre_pass_tbkg_total,"thickline", ROOT.kBlue, 1)
#	GoodPlotFormat(pre_fail_tbkg_total,"fill", ROOT.kBlue, 0, 1001)
#	GoodPlotFormat(pre_pass_tbkg_total,"fill", ROOT.kBlue, 0, 1001)
	GoodPlotFormat(b_fail_tbkg_total,"thickline", ROOT.kBlue, 1)
#	GoodPlotFormat(b_fail_tbkg_total,"fill", ROOT.kBlue, 0, 1001)
#	GoodPlotFormat(b_pass_tbkg_total,"fill", ROOT.kBlue, 0, 1001)
#	GoodPlotFormat(s_fail_tbkg_total,"fill", ROOT.kBlue, 0, 1001)
#	GoodPlotFormat(s_pass_tbkg_total,"fill", ROOT.kBlue, 0, 1001)
	GoodPlotFormat(b_pass_tbkg_total,"thickline", ROOT.kBlue, 1)
	GoodPlotFormat(s_fail_tbkg_total,"thickline", ROOT.kBlue, 1)
	GoodPlotFormat(s_pass_tbkg_total,"thickline", ROOT.kBlue, 1)
#	GoodPlotFormat(pre_fail_ttbar_total,"thickline", ROOT.kRed, 1)
#	GoodPlotFormat(pre_pass_ttbar_total,"thickline", ROOT.kRed, 1)
	GoodPlotFormat(pre_fail_ttbar_total, "fill", ROOT.kRed, ROOT.kRed, 1001)
	GoodPlotFormat(pre_pass_ttbar_total, "fill", ROOT.kRed, ROOT.kRed, 1001)
	GoodPlotFormat(b_fail_ttbar_total, "fill", ROOT.kRed, ROOT.kRed, 1001)
	GoodPlotFormat(b_pass_ttbar_total, "fill", ROOT.kRed, ROOT.kRed, 1001)
	GoodPlotFormat(s_fail_ttbar_total, "fill", ROOT.kRed, ROOT.kRed, 1001)
	GoodPlotFormat(s_pass_ttbar_total, "fill", ROOT.kRed, ROOT.kRed, 1001)
#	GoodPlotFormat(b_pass_ttbar_total,"thickline", ROOT.kRed, 1)
#	GoodPlotFormat(s_fail_ttbar_total,"thickline", ROOT.kRed, 1)
#	GoodPlotFormat(s_pass_ttbar_total,"thickline", ROOT.kRed, 1)
	GoodPlotFormat(pre_fail_tsig10_total, "fill", ROOT.kBlue+1, ROOT.kBlue+1, 3003)
	GoodPlotFormat(pre_pass_tsig10_total, "fill", ROOT.kBlue+1, ROOT.kBlue+1, 3003)
	GoodPlotFormat(b_fail_tsig10_total, "fill", ROOT.kBlue+1, ROOT.kBlue+1, 3003)
	GoodPlotFormat(b_pass_tsig10_total, "fill", ROOT.kBlue+1, ROOT.kBlue+1, 3003)
	GoodPlotFormat(s_fail_tsig10_total, "fill", ROOT.kBlue+1, ROOT.kBlue+1, 3003)
	GoodPlotFormat(s_pass_tsig10_total, "fill", ROOT.kBlue+1, ROOT.kBlue+1, 3003)
	GoodPlotFormat(pre_fail_tsig25_total, "fill", ROOT.kGreen+1, ROOT.kGreen+1, 3003)
	GoodPlotFormat(pre_pass_tsig25_total, "fill", ROOT.kGreen+1, ROOT.kGreen+1, 3003)
	GoodPlotFormat(b_fail_tsig25_total, "fill", ROOT.kGreen+1, ROOT.kGreen+1, 3003)
	GoodPlotFormat(b_pass_tsig25_total, "fill", ROOT.kGreen+1, ROOT.kGreen+1, 3003)
	GoodPlotFormat(s_fail_tsig25_total, "fill", ROOT.kGreen+1, ROOT.kGreen+1, 3003)
	GoodPlotFormat(s_pass_tsig25_total, "fill", ROOT.kGreen+1, ROOT.kGreen+1, 3003)
	GoodPlotFormat(pre_fail_tsig50_total, "fill", ROOT.kGray+1, ROOT.kGray+1, 3003)
	GoodPlotFormat(pre_pass_tsig50_total, "fill", ROOT.kGray+1, ROOT.kGray+1, 3003)
	GoodPlotFormat(b_fail_tsig50_total, "fill", ROOT.kGray+1, ROOT.kGray+1, 3003)
	GoodPlotFormat(b_pass_tsig50_total, "fill", ROOT.kGray+1, ROOT.kGray+1, 3003)
	GoodPlotFormat(s_fail_tsig50_total, "fill", ROOT.kGray+1, ROOT.kGray+1, 3003)
	GoodPlotFormat(s_pass_tsig50_total, "fill", ROOT.kGray+1, ROOT.kGray+1, 3003)
	GoodPlotFormat(pre_fail_tsig75_total, "fill", ROOT.kYellow+1, ROOT.kYellow+1, 3003)
	GoodPlotFormat(pre_pass_tsig75_total, "fill", ROOT.kYellow+1, ROOT.kYellow+1, 3003)
	GoodPlotFormat(b_fail_tsig75_total, "fill", ROOT.kYellow+1, ROOT.kYellow+1, 3003)
	GoodPlotFormat(b_pass_tsig75_total, "fill", ROOT.kYellow+1, ROOT.kYellow+1, 3003)
	GoodPlotFormat(s_fail_tsig75_total, "fill", ROOT.kYellow+1, ROOT.kYellow+1, 3003)
	GoodPlotFormat(s_pass_tsig75_total, "fill", ROOT.kYellow+1, ROOT.kYellow+1, 3003)
#        GoodPlotFormat(pre_fail_wgamma_total,"thickline", ROOT.kViolet, 1)
#        GoodPlotFormat(pre_pass_wgamma_total,"thickline", ROOT.kViolet, 1)
        GoodPlotFormat(pre_fail_wgamma_total,"fill", ROOT.kViolet, ROOT.kViolet, 1001)
        GoodPlotFormat(pre_pass_wgamma_total,"fill", ROOT.kViolet, ROOT.kViolet, 1001)
        GoodPlotFormat(b_fail_wgamma_total,"fill", ROOT.kViolet, ROOT.kViolet, 1001)
        GoodPlotFormat(b_pass_wgamma_total,"fill", ROOT.kViolet, ROOT.kViolet, 1001)
        GoodPlotFormat(s_fail_wgamma_total,"fill", ROOT.kViolet, ROOT.kViolet, 1001)
        GoodPlotFormat(s_pass_wgamma_total,"fill", ROOT.kViolet, ROOT.kViolet, 1001)
#        GoodPlotFormat(b_pass_wgamma_total,"thickline", ROOT.kViolet, 1)
#        GoodPlotFormat(s_fail_wgamma_total,"thickline", ROOT.kViolet, 1)
#        GoodPlotFormat(s_pass_wgamma_total,"thickline", ROOT.kViolet, 1)
      #  GoodPlotFormat(pre_fail_zgamma_total,"thickline", ROOT.kOrange, 1)
       # GoodPlotFormat(pre_pass_zgamma_total,"thickline", ROOT.kOrange, 1)
        GoodPlotFormat(pre_fail_zgamma_total,"fill", ROOT.kOrange, ROOT.kOrange, 1001)
        GoodPlotFormat(pre_pass_zgamma_total,"fill", ROOT.kOrange, ROOT.kOrange, 1001)
        GoodPlotFormat(b_fail_zgamma_total,"fill", ROOT.kOrange, ROOT.kOrange, 1001)
        GoodPlotFormat(b_pass_zgamma_total,"fill", ROOT.kOrange, ROOT.kOrange, 1001)
        GoodPlotFormat(s_fail_zgamma_total,"fill", ROOT.kOrange, ROOT.kOrange, 1001)
        GoodPlotFormat(s_pass_zgamma_total,"fill", ROOT.kOrange, ROOT.kOrange, 1001)
#        GoodPlotFormat(b_pass_zgamma_total,"thickline", ROOT.kOrange, 1)
#        GoodPlotFormat(s_fail_zgamma_total,"thickline", ROOT.kOrange, 1)
#        GoodPlotFormat(s_pass_zgamma_total,"thickline", ROOT.kOrange, 1)
#        GoodPlotFormat(pre_fail_nres_total,"dashed", ROOT.kCyan, 1)
#        GoodPlotFormat(pre_pass_nres_total,"dashed", ROOT.kCyan, 1)
        GoodPlotFormat(pre_fail_nres_total,"fill", ROOT.kCyan, 10, 1)
        GoodPlotFormat(pre_pass_nres_total,"fill", ROOT.kCyan, 10, 1)
        GoodPlotFormat(b_fail_nres_total,"fill", ROOT.kCyan, 10, 1)
        GoodPlotFormat(b_pass_nres_total,"fill", ROOT.kCyan, 10, 1)
        GoodPlotFormat(s_fail_nres_total,"fill", ROOT.kCyan, 10, 1)
        GoodPlotFormat(s_pass_nres_total,"fill", ROOT.kCyan, 10, 1)
#        GoodPlotFormat(b_pass_nres_total,"dashed", ROOT.kCyan, 1)
#        GoodPlotFormat(s_fail_nres_total,"dashed", ROOT.kCyan, 1)
#        GoodPlotFormat(s_pass_nres_total,"dashed", ROOT.kCyan, 1)




	ofile.WriteObject(b_pull_count_total, "fit_b_pass_pull_count_total")
	ofile.WriteObject(s_pull_count_total, "multiSig_pass_pull_count_total")
	b_pull_count.Reset()
	s_pull_count.Reset()
	ofile2.WriteObject(pre_fail_data_total, "pre_fail_data_total")
	ofile2.WriteObject(pre_pass_data_total, "pre_pass_data_total")
	ofile2.WriteObject(b_fail_data_total, "fit_b_fail_data_total")
	ofile2.WriteObject(b_pass_data_total, "fit_b_pass_data_total")
	ofile2.WriteObject(s_fail_data_total, "multiSig_fail_data_total")
	ofile2.WriteObject(s_pass_data_total, "multiSig_pass_data_total")
	ofile2.WriteObject(pre_fail_ttbar_total, "pre_fail_ttbar_total")
	ofile2.WriteObject(pre_pass_ttbar_total, "pre_pass_ttbar_total")
	ofile2.WriteObject(b_fail_ttbar_total, "fit_b_fail_ttbar_total")
	ofile2.WriteObject(b_pass_ttbar_total, "fit_b_pass_ttbar_total")
	ofile2.WriteObject(s_fail_ttbar_total, "multiSig_fail_ttbar_total")
	ofile2.WriteObject(s_pass_ttbar_total, "multiSig_pass_ttbar_total")
	ofile2.WriteObject(pre_fail_tbkg_total, "pre_fail_tbkg_total")
	ofile2.WriteObject(pre_pass_tbkg_total, "pre_pass_tbkg_total")
	ofile2.WriteObject(b_fail_tbkg_total, "fit_b_fail_tbkg_total")
	ofile2.WriteObject(b_pass_tbkg_total, "fit_b_pass_tbkg_total")
	ofile2.WriteObject(s_fail_tbkg_total, "multiSig_fail_tbkg_total")
	ofile2.WriteObject(s_pass_tbkg_total, "multiSig_pass_tbkg_total")
	ofile2.WriteObject(pre_fail_tsig10_total, "pre_fail_tsig10_total")
	ofile2.WriteObject(pre_pass_tsig10_total, "pre_pass_tsig10_total")
	ofile2.WriteObject(b_fail_tsig10_total, "fit_b_fail_tsig10_total")
	ofile2.WriteObject(b_pass_tsig10_total, "fit_b_pass_tsig10_total")
	ofile2.WriteObject(s_fail_tsig10_total, "sig10_fail_tsig10_total")
	ofile2.WriteObject(s_pass_tsig10_total, "sig10_pass_tsig10_total")
	ofile2.WriteObject(pre_fail_tsig25_total, "pre_fail_tsig25_total")
	ofile2.WriteObject(pre_pass_tsig25_total, "pre_pass_tsig25_total")
	ofile2.WriteObject(b_fail_tsig25_total, "fit_b_fail_tsig25_total")
	ofile2.WriteObject(b_pass_tsig25_total, "fit_b_pass_tsig25_total")
	ofile2.WriteObject(s_fail_tsig25_total, "sig25_fail_tsig25_total")
	ofile2.WriteObject(s_pass_tsig25_total, "sig25_pass_tsig25_total")
	ofile2.WriteObject(pre_fail_tsig50_total, "pre_fail_tsig50_total")
	ofile2.WriteObject(pre_pass_tsig50_total, "pre_pass_tsig50_total")
	ofile2.WriteObject(b_fail_tsig50_total, "fit_b_fail_tsig50_total")
	ofile2.WriteObject(b_pass_tsig50_total, "fit_b_pass_tsig50_total")
	ofile2.WriteObject(s_fail_tsig50_total, "sig50_fail_tsig50_total")
	ofile2.WriteObject(s_pass_tsig50_total, "sig50_pass_tsig50_total")
	ofile2.WriteObject(pre_fail_tsig75_total, "pre_fail_tsig75_total")
	ofile2.WriteObject(pre_pass_tsig75_total, "pre_pass_tsig75_total")
	ofile2.WriteObject(b_fail_tsig75_total, "fit_b_fail_tsig75_total")
	ofile2.WriteObject(b_pass_tsig75_total, "fit_b_pass_tsig75_total")
	ofile2.WriteObject(s_fail_tsig75_total, "sig75_fail_tsig75_total")
	ofile2.WriteObject(s_pass_tsig75_total, "sig75_pass_tsig75_total")
	ofile2.WriteObject(pre_fail_wgamma_total, "pre_fail_wgamma_total")
	ofile2.WriteObject(pre_pass_wgamma_total, "pre_pass_wgamma_total")
	ofile2.WriteObject(b_fail_wgamma_total, "fit_b_fail_wgamma_total")
	ofile2.WriteObject(b_pass_wgamma_total, "fit_b_pass_wgamma_total")
	ofile2.WriteObject(s_fail_wgamma_total, "multiSig_fail_wgamma_total")
	ofile2.WriteObject(s_pass_wgamma_total, "multiSig_pass_wgamma_total")
	ofile2.WriteObject(pre_fail_zgamma_total, "pre_fail_zgamma_total")
	ofile2.WriteObject(pre_pass_zgamma_total, "pre_pass_zgamma_total")
	ofile2.WriteObject(b_fail_zgamma_total, "fit_b_fail_zgamma_total")
	ofile2.WriteObject(b_pass_zgamma_total, "fit_b_pass_zgamma_total")
	ofile2.WriteObject(s_fail_zgamma_total, "multiSig_fail_zgamma_total")
	ofile2.WriteObject(s_pass_zgamma_total, "multiSig_pass_zgamma_total")
	ofile2.WriteObject(pre_fail_nres_total, "pre_fail_nres_total")
	ofile2.WriteObject(pre_pass_nres_total, "pre_pass_nres_total")
	ofile2.WriteObject(b_fail_nres_total, "fit_b_fail_nres_total")
	ofile2.WriteObject(b_pass_nres_total, "fit_b_pass_nres_total")
	ofile2.WriteObject(s_fail_nres_total, "multiSig_fail_nres_total")
	ofile2.WriteObject(s_pass_nres_total, "multiSig_pass_nres_total")
	
	C = ROOT.TCanvas()
	C.cd()
	s_fail_nres_total.SetTitle("Failing Non-Resonant Background (Closure Test)")
	s_fail_nres_total.GetXaxis().SetTitle("Softdrop Mass")
	s_fail_nres_total.GetYaxis().SetTitle("Events")
	s_fail_nres_total.Draw("hist")
	C.Print("results/nres_fail_allbins_"+str(name)+".png")
	C.Close()

	C = ROOT.TCanvas()
	C.cd()
	s_pass_nres_total.SetTitle("Passing Non-Resonant Background (Closure Test)")
	s_pass_nres_total.GetXaxis().SetTitle("Softdrop Mass")
	s_pass_nres_total.GetYaxis().SetTitle("Events")
	s_pass_nres_total.Draw("hist")
	C.Print("results/nres_pass_allbins_"+str(name)+".png")
	C.Close()

	C = ROOT.TCanvas()
	C.cd()
	pre_fail_data_total.SetTitle("Failing Pre-Fit Data (Closure Test)")
	pre_fail_data_total.GetXaxis().SetTitle("Softdrop Mass")
	pre_fail_data_total.GetYaxis().SetTitle("Events")
	pre_fail_data_total.Draw("hist")
	C.Print("results/prefit_data_fail_allbins_"+str(name)+".png")
	C.Close()

	C = ROOT.TCanvas()
	C.cd()
	pre_pass_data_total.SetTitle("Passing Pre-Fit Data (Closure Test)")
	pre_pass_data_total.GetXaxis().SetTitle("Softdrop Mass")
	pre_pass_data_total.GetYaxis().SetTitle("Events")
	pre_pass_data_total.Draw("hist")
	C.Print("results/prefit_data_pass_allbins_"+str(name)+".png")
	C.Close()
	
	C = ROOT.TCanvas()
	C.cd()
	s_fail_data_total.SetTitle("Failing Data (Closure Test)")
	s_fail_data_total.GetXaxis().SetTitle("Softdrop Mass")
	s_fail_data_total.GetYaxis().SetTitle("Events")
	s_fail_data_total.Draw("hist")
	C.Print("results/multiSig_data_fail_allbins_"+str(name)+".png")
	C.Close()

	C = ROOT.TCanvas()
	C.cd()
	s_pass_data_total.SetTitle("Passing Data (Closure Test)")
	s_pass_data_total.GetXaxis().SetTitle("Softdrop Mass")
	s_pass_data_total.GetYaxis().SetTitle("Events")
	s_pass_data_total.Draw("hist")
	C.Print("results/multiSig_data_pass_allbins_"+str(name)+".png")
	C.Close()

	b_cheapline = b_pass_data_total.Clone("b_cheapline")
	b_cheapline.Add(b_pass_data_total,-1.)
	b_cheapline.SetTitle("")
	b_cheapline.GetYaxis().SetTitle("#frac{data - bkg}{#sigma_{data}}")
	b_cheapline.GetYaxis().SetTitleSize(0.175);
	b_cheapline.GetYaxis().SetNdivisions(6);
	b_cheapline.GetYaxis().SetLabelSize(0.145);
	b_cheapline.GetYaxis().SetTitleOffset(0.225);
	b_cheapline.GetYaxis().CenterTitle(True)
	b_cheapline.GetYaxis().SetRangeUser(-5.,5.)
	GoodPlotFormat(b_cheapline, "thinline", ROOT.kGray, 4)
	b_pass_tbkg_total.GetYaxis().SetTitle("Events / 5 GeV")
	b_pass_tbkg_total.GetYaxis().SetTitleOffset(0.5);
	b_pass_tbkg_total.GetYaxis().SetTitleSize(0.075);
	FindAndSetMax(b_pass_data_total, b_pass_tbkg_total)
						
	E = []
	EP = []
	for i in range(1,b_pass_data_total.GetNbinsX()+1):
		Err = b_pass_tbkg_total.GetBinError(i)
				
		blX = b_pass_tbkg_total.GetBinLowEdge(i)
		blY = b_pass_tbkg_total.GetBinContent(i) - Err
		trX = b_pass_tbkg_total.GetBinWidth(i) + blX
		trY = b_pass_tbkg_total.GetBinContent(i) + Err
		tBox = ROOT.TBox(blX,blY,trX,trY)
		if  b_pass_data_total.GetBinError(i) > 0:
			ue = Err/b_pass_data_total.GetBinError(i)
		else:
			ue = Err/2.7
		ue = min(5.0, ue)
		tPBox = ROOT.TBox(blX, -1*ue, trX, ue)
		tBox.SetFillColor(25)
		tBox.SetFillStyle(3144)
		tPBox.SetFillColor(25)
		tPBox.SetFillStyle(3144)
		tBox.SetLineColor(ROOT.kWhite)
		tPBox.SetLineColor(ROOT.kWhite)
		E.append(tBox)
		EP.append(tPBox)
										
	b_pass_tbkg_total.GetXaxis().SetLabelSize(0)
			
	pull = b_pass_data_total.Clone("pull")
	pull.Add(b_pass_tbkg_total, -1.)
	GoodPlotFormat(pull, "markers", ROOT.kBlack, 20)
	for i in range(pull.GetNbinsX()):
		if not b_pass_data_total.GetBinContent(i+1) == 0:
			pull.SetBinContent(i+1, pull.GetBinContent(i+1)/b_pass_data_total.GetBinError(i+1))
			pull.SetBinError(i+1, 1)		
			#Fill Pull Counter
			b_pull_count.Fill(pull.GetBinContent(i+1))
		else:
			pull.SetBinContent(i+1, 0)
			pull.SetBinError(i+1, 0)
					
	# Stack backgrounds for plotting
	b_pass_zgamma_total.Add(b_pass_nres_total)
	b_pass_wgamma_total.Add(b_pass_zgamma_total)
	b_pass_ttbar_total.Add(b_pass_wgamma_total)

	L = ROOT.TLegend(0.62,0.54,0.90,0.88)
	#L = ROOT.TLegend(0.62,0.7,0.92,0.86)
	#L = ROOT.TLegend(0.6,0.65,0.9,0.9)
	L.SetFillColor(0)
	L.SetLineColor(0)
	if not Blind: L.AddEntry(b_pass_data_total, "Data", "PE")
	L.AddEntry(b_pass_tbkg_total, "Total Background", "L")
	L.AddEntry(b_pass_ttbar_total, "t#bar{t} component", "F")
        L.AddEntry(b_pass_wgamma_total, "W+Gamma", "F")
        L.AddEntry(b_pass_zgamma_total, "Z+Gamma", "F")
	L.AddEntry(b_pass_nres_total, "Non-Resonant Background", "L")
	L.AddEntry(E[0], "Background Uncertainty", "F")
#	if P != "fit_b":  L.AddEntry(b_pass_tsig25_total, b_pass_tsig25_total_name, "F")
	C = ROOT.TCanvas()
	C.cd()
	p12 = ROOT.TPad("pad1", "tall",0,0.165,1,1)
	p22 = ROOT.TPad("pad2", "short",0,0.0,1.0,0.23)
	p22.SetBottomMargin(0.35)
	p12.Draw()
	p22.Draw()
	p12.cd() # top
	ROOT.gPad.SetTicks(1,1)
	b_pass_tbkg_total.SetTitle("Full Run2 (10%) Background Only Fit Plot All Bins")
	b_pass_tbkg_total.Draw("hist")
	b_pass_ttbar_total.Draw("histsame")
        b_pass_wgamma_total.Draw("histsame")
        b_pass_zgamma_total.Draw("histsame")
	b_pass_nres_total.Draw("histsame")
#	if P != "fit_b": b_pass_tsig25_total.Draw("histsame")
	for e in E: e.Draw("same")
	b_pass_data_total.Draw("esame")
	L.Draw("same")
	ROOT.TGaxis.SetMaxDigits(3)
	p12.RedrawAxis()
	AddCMSLumi(ROOT.gPad, LUMI['Run2'], cmsextra)

	# Add Chi2 GOF
	#chi2, ndf, bins = Chi2(b_pass_data_total, b_pass_tbkg_total, 9) #9 fit parameters for 2,2 Bernstein
	GOF(ROOT.gPad, b_chi, b_bins-27)

#	pTbin(ROOT.gPad, "Full")
	p22.cd() # bottom
	ROOT.gPad.SetTicks(1,1)
	b_cheapline.GetXaxis().SetTitle("Softdrop Mass [GeV]")
	b_cheapline.GetXaxis().SetTitleSize(0.1925);
	b_cheapline.GetXaxis().SetLabelSize(0.16);
	b_cheapline.GetXaxis().SetTitleOffset(0.84);
	b_cheapline.Draw("hist")
	for e in EP: e.Draw("same")
	l1 = TLine(0,2,200,2)
	l1.SetLineColor(kGray)
	l1.SetLineWidth(1)
	l1.SetLineStyle(7)
	l2 = TLine(0,-2,200,-2)
	l2.SetLineColor(kGray)
	l2.SetLineWidth(1)
	l2.SetLineStyle(7)
	pull.GetXaxis().SetTitle("Softdrop Mass [GeV]")
	pull.Draw("esame")
	l1.Draw("same")
	l2.Draw("same")
	p22.RedrawAxis()
	C.Print("results/fit_b_pass_allbins_"+str(name)+".png")
	p12.cd()
	gPad.SetLogy()
	p12.Draw()
	C.Print("results/fit_b_pass_allbins_"+str(name)+"_logy.png")
	b_pull_count.SetTitle("Passing Pull Count Total Background")
	pull_name = "fit_b_pass_pull_count_totalbkg"
	ofile.WriteObject(b_pull_count, pull_name)

	C.Close()
	

	s_cheapline = s_pass_data_total.Clone("s_cheapline")
	s_cheapline.Add(s_pass_data_total,-1.)
	s_cheapline.SetTitle("")
	s_cheapline.GetYaxis().SetTitle("#frac{data - bkg}{#sigma_{data}}")
	s_cheapline.GetYaxis().SetTitleSize(0.175);
	s_cheapline.GetYaxis().SetNdivisions(6);
	s_cheapline.GetYaxis().SetLabelSize(0.145);
	s_cheapline.GetYaxis().SetTitleOffset(0.225);
	s_cheapline.GetYaxis().CenterTitle(True)
	s_cheapline.GetYaxis().SetRangeUser(-5.,5.)
	GoodPlotFormat(s_cheapline, "thinline", ROOT.kGray, 4)
	s_pass_tbkg_total.GetYaxis().SetTitle("Events / 5 GeV")
	s_pass_tbkg_total.GetYaxis().SetTitleOffset(0.5);
	s_pass_tbkg_total.GetYaxis().SetTitleSize(0.075);
	FindAndSetMax(s_pass_data_total, s_pass_tbkg_total)
	
						
	E = []
	EP = []
	for i in range(1,s_pass_data_total.GetNbinsX()+1):
		Err = s_pass_tbkg_total.GetBinError(i)
				
		blX = s_pass_tbkg_total.GetBinLowEdge(i)
		blY = s_pass_tbkg_total.GetBinContent(i) - Err
		trX = s_pass_tbkg_total.GetBinWidth(i) + blX
		trY = s_pass_tbkg_total.GetBinContent(i) + Err
		tBox = ROOT.TBox(blX,blY,trX,trY)
		if  s_pass_data_total.GetBinError(i) > 0:
			ue = Err/s_pass_data_total.GetBinError(i)
		else:
			ue = Err/2.7
		ue = min(5.0, ue)
		tPBox = ROOT.TBox(blX, -1*ue, trX, ue)
		tBox.SetFillColor(25)
		tBox.SetFillStyle(3144)
		tPBox.SetFillColor(25)
		tPBox.SetFillStyle(3144)
		tBox.SetLineColor(ROOT.kWhite)
		tPBox.SetLineColor(ROOT.kWhite)
		E.append(tBox)
		EP.append(tPBox)
										
	s_pass_tbkg_total.GetXaxis().SetLabelSize(0)
			
	pull = s_pass_data_total.Clone("pull")
	pull.Add(s_pass_tbkg_total, -1.)
	GoodPlotFormat(pull, "markers", ROOT.kBlack, 20)
	for i in range(pull.GetNbinsX()):
		if not s_pass_data_total.GetBinContent(i+1) == 0:
			pull.SetBinContent(i+1, pull.GetBinContent(i+1)/s_pass_data_total.GetBinError(i+1))
			pull.SetBinError(i+1, 1)		
			#Fill Pull Counter
			s_pull_count.Fill(pull.GetBinContent(i+1))
		else:
			pull.SetBinContent(i+1, 0)
			pull.SetBinError(i+1, 0)
					
	# Stack backgrounds for plotting
	s_pass_zgamma_total.Add(s_pass_nres_total)
	s_pass_wgamma_total.Add(s_pass_zgamma_total)
	s_pass_ttbar_total.Add(s_pass_wgamma_total)

	L = ROOT.TLegend(0.62,0.54,0.90,0.88)
	#L = ROOT.TLegend(0.6,0.55,0.9,0.86)
	#L = ROOT.TLegend(0.58,0.65,0.895,0.895)
	L.SetFillColor(0)
	L.SetLineColor(0)
	if not Blind: L.AddEntry(s_pass_data_total, "Data", "PE")
	L.AddEntry(s_pass_tbkg_total, "Total background", "L")
	L.AddEntry(s_pass_ttbar_total, "t#bar{t} component", "F")
        L.AddEntry(s_pass_wgamma_total, "W+Gamma", "F")
        L.AddEntry(s_pass_zgamma_total, "Z+Gamma", "F")
	L.AddEntry(s_pass_nres_total, "Non-Resonant Background", "L")
	L.AddEntry(s_pass_tsig10_total, "10 GeV Z' Signal", "F")
	L.AddEntry(s_pass_tsig25_total, "25 GeV Z' Signal", "F")
	L.AddEntry(s_pass_tsig50_total, "50 GeV Z' Signal", "F")
	L.AddEntry(s_pass_tsig75_total, "75 GeV Z' Signal", "F")
	L.AddEntry(E[0], "Background Uncertainty", "F")
#	if P != "fit_b":  L.AddEntry(s_pass_tsig25_total, s_pass_tsig25_total_name, "F")
	C = ROOT.TCanvas()
	C.cd()
	p12 = ROOT.TPad("pad1", "tall",0,0.165,1,1)
	p22 = ROOT.TPad("pad2", "short",0,0.0,1.0,0.23)
	p22.SetBottomMargin(0.35)
	p12.Draw()
	p22.Draw()
	p12.cd() # top
	ROOT.gPad.SetTicks(1,1)
	s_pass_tbkg_total.SetTitle("Run 2 (10%) Multiple Signal Background Fit Plot With Signal All Bins")
	s_pass_tbkg_total.Draw("hist")
	s_pass_ttbar_total.Draw("histsame")
        s_pass_wgamma_total.Draw("histsame")
        s_pass_zgamma_total.Draw("histsame")
	s_pass_nres_total.Draw("histsame")
	s_pass_tsig10_total.Draw("histsame")
	s_pass_tsig25_total.Draw("histsame")
	s_pass_tsig50_total.Draw("histsame")
	s_pass_tsig75_total.Draw("histsame")
#	if P != "fit_b": s_pass_tsig25_total.Draw("histsame")
	for e in E: e.Draw("same")
	s_pass_data_total.Draw("esame")
	L.Draw("same")
	ROOT.TGaxis.SetMaxDigits(3)
	p12.RedrawAxis()
	AddCMSLumi(ROOT.gPad, LUMI['Run2'], cmsextra)
	
	# Add Chi2 GOF
	#chi2, ndf, bins = Chi2(s_pass_data_total, s_pass_tbkg_total, 9) #9 fit parameters for 2,2 Bernstein
	GOF(ROOT.gPad, s_chi, s_bins-27)

#	pTbin(ROOT.gPad, "Full")
	p22.cd() # bottom
	ROOT.gPad.SetTicks(1,1)
	s_cheapline.GetXaxis().SetTitle("Softdrop Mass [GeV]")
	s_cheapline.GetXaxis().SetTitleSize(0.1925);
	s_cheapline.GetXaxis().SetLabelSize(0.16);
	s_cheapline.GetXaxis().SetTitleOffset(0.84);
	s_cheapline.Draw("hist")
	for e in EP: e.Draw("same")
	l1 = TLine(0,2,200,2)
	l1.SetLineColor(kGray)
	l1.SetLineWidth(1)
	l1.SetLineStyle(7)
	l2 = TLine(0,-2,200,-2)
	l2.SetLineColor(kGray)
	l2.SetLineWidth(1)
	l2.SetLineStyle(7)
	pull.Draw("esame")
	l1.Draw("same")
	l2.Draw("same")
	p22.RedrawAxis()
	C.Print("results/multiSig_pass_allbins_"+str(name)+".png")
	p12.cd()
	gPad.SetLogy()
	p12.Draw()
	C.Print("results/multiSig_pass_allbins_"+str(name)+"_logy.png")
	s_pull_count.SetTitle("Passing Pull Count Total Background")
	pull_name = "sig_pass_pull_count_totalbkg"
	ofile.WriteObject(s_pull_count, pull_name)

