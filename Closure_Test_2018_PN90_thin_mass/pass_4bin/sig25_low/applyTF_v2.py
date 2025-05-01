#Copied from https://github.com/osherson/B2GAnalysis_2020/blob/master/analysis/CombineStep1.py
#And modified for my analysis
import ROOT
from ROOT import *
import numpy as np
import math
import sys
import array
import os
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
		H.SetFillColor(args[1])
		H.SetFillStyle(args[2])
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
		j.GetYaxis().SetRangeUser(0,maximum*1.35)#should be 1.35 (below as well)
		j.SetLineWidth(2)
	return maximum*1.35

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
	latex.DrawLatex(0.6, 0.525, cmsText)
	pad.Update()

def drawDiagnostic(name, ifile1):
	gROOT.SetBatch(True)
	gStyle.SetOptStat(0)
	Sigs = [["/home/akobert/CMSSW_11_1_0_pre7/src/RData/NanoTool_UL/M25_UL_nano_merged.root"], "5.982*7639.0/567896.0", "", "sig25_nano", "25 GeV Signal", "7639.0"]
	sig_name = "sig25_nano"
    	ptBins = {} # Dictonary of bins for each year
    	ptBins['2016'] = np.array([200, 240, 275, 330, 700])
    	ptBins['2017'] = np.array([220, 260, 300, 360, 700])
    	ptBins['2018'] = np.array([120, 160, 200, 255, 700])
	mBins = MakeNBinsFromMinToMax(40,0,200)
	VAR = ["jpt", ptBins, "Jet pT (GeV)", "sdm", mBins, "Softdrop Mass (GeV)"]

	Blind = True
	LUMI = {}
	LUMI['2016'] = 3.631
	LUMI['2017'] = 4.148
	LUMI['2018'] = 5.982
	LUMI['Run2'] = 13.761 # 10% of Luminosity
	cmsextra = "Preliminary"

	ofile = ROOT.TFile("./TF_plots.root", "RECREATE")


        ROOT.gInterpreter.Declare("Double_t widebins_2016[6] = {200, 240, 275, 330, 700};")
        ROOT.gInterpreter.Declare("Double_t widebins_2017[6] = {220, 256, 300, 360, 700};")
        ROOT.gInterpreter.Declare("Double_t widebins_2018[6] = {120, 160, 200, 255, 700};")


	c0 = TCanvas()
	c0.cd()
	
	CbnF = ROOT.TFile("fitDiagnosticsTest.root")	

	mpt_TF = {}
	pass_data = {}
	fail_data = {}
	pass_nres = {}
	fail_nres = {}
	pass_ttbar = {}
	fail_ttbar = {}
	pass_wgamma = {}
	fail_wgamma = {}
	pass_zgamma = {}
	fail_zgamma = {}
	pass_tbkg = {}
	fail_tbkg = {}

	years = ['2016', '2017', '2018']
	for y in years:
		if y == '2016':
			mpt_TF[y] = TH2F("mpt_TF_"+y, "pT vs. Softdrop Mass Transfer Factor", 4, widebins_2016, 40, 0, 200)
			pass_data[y] = TH2F("pass_data_"+y, "Passing Data "+y, 4, widebins_2016, 40, 0, 200)
			fail_data[y] = TH2F("fail_data_"+y, "Failing Data "+y, 4, widebins_2016, 40, 0, 200)
			pass_nres[y] = TH2F("pass_nres_"+y, "Passing Non-Resonant Background "+y, 4, widebins_2016, 40, 0, 200)
			fail_nres[y] = TH2F("fail_nres_"+y, "Failing Non-Resonant Background "+y, 4, widebins_2016, 40, 0, 200)
			pass_ttbar[y] = TH2F("pass_ttbar_"+y, "Passing TTBar "+y, 4, widebins_2016, 40, 0, 200)
			fail_ttbar[y] = TH2F("fail_ttbar_"+y, "Failing TTBar "+y, 4, widebins_2016, 40, 0, 200)
			pass_wgamma[y] = TH2F("pass_wgamma_"+y, "Passing WGamma "+y, 4, widebins_2016, 40, 0, 200)
			fail_wgamma[y] = TH2F("fail_wgamma_"+y, "Failing WGamma "+y, 4, widebins_2016, 40, 0, 200)
			pass_zgamma[y] = TH2F("pass_zgamma_"+y, "Passing ZGamma "+y, 4, widebins_2016, 40, 0, 200)
			fail_zgamma[y] = TH2F("fail_zgamma_"+y, "Failing ZGamma "+y, 4, widebins_2016, 40, 0, 200)
			pass_tbkg[y] = TH2F("pass_tbkg_"+y, "Passing Total-Background "+y, 4, widebins_2016, 40, 0, 200)
			fail_tbkg[y] = TH2F("fail_tbkg_"+y, "Failing Total-Background "+y, 4, widebins_2016, 40, 0, 200)
		if y == '2017':
			mpt_TF[y] = TH2F("mpt_TF_"+y, "pT vs. Softdrop Mass Transfer Factor", 4, widebins_2017, 40, 0, 200)
			pass_data[y] = TH2F("pass_data_"+y, "Passing Data "+y, 4, widebins_2017, 40, 0, 200)
			fail_data[y] = TH2F("fail_data_"+y, "Failing Data "+y, 4, widebins_2017, 40, 0, 200)
			pass_nres[y] = TH2F("pass_nres_"+y, "Passing Non-Resonant Background "+y, 4, widebins_2017, 40, 0, 200)
			fail_nres[y] = TH2F("fail_nres_"+y, "Failing Non-Resonant Background "+y, 4, widebins_2017, 40, 0, 200)
			pass_ttbar[y] = TH2F("pass_ttbar_"+y, "Passing TTBar "+y, 4, widebins_2017, 40, 0, 200)
			fail_ttbar[y] = TH2F("fail_ttbar_"+y, "Failing TTBar "+y, 4, widebins_2017, 40, 0, 200)
			pass_wgamma[y] = TH2F("pass_wgamma_"+y, "Passing WGamma "+y, 4, widebins_2017, 40, 0, 200)
			fail_wgamma[y] = TH2F("fail_wgamma_"+y, "Failing WGamma "+y, 4, widebins_2017, 40, 0, 200)
			pass_zgamma[y] = TH2F("pass_zgamma_"+y, "Passing ZGamma "+y, 4, widebins_2017, 40, 0, 200)
			fail_zgamma[y] = TH2F("fail_zgamma_"+y, "Failing ZGamma "+y, 4, widebins_2017, 40, 0, 200)
			pass_tbkg[y] = TH2F("pass_tbkg_"+y, "Passing Total-Background "+y, 4, widebins_2017, 40, 0, 200)
			fail_tbkg[y] = TH2F("fail_tbkg_"+y, "Failing Total-Background "+y, 4, widebins_2017, 40, 0, 200)
		if y == '2018':
			mpt_TF[y] = TH2F("mpt_TF_"+y, "pT vs. Softdrop Mass Transfer Factor", 4, widebins_2018, 40, 0, 200)
			pass_data[y] = TH2F("pass_data_"+y, "Passing Data "+y, 4, widebins_2018, 40, 0, 200)
			fail_data[y] = TH2F("fail_data_"+y, "Failing Data "+y, 4, widebins_2018, 40, 0, 200)
			pass_nres[y] = TH2F("pass_nres_"+y, "Passing Non-Resonant Background "+y, 4, widebins_2018, 40, 0, 200)
			fail_nres[y] = TH2F("fail_nres_"+y, "Failing Non-Resonant Background "+y, 4, widebins_2018, 40, 0, 200)
			pass_ttbar[y] = TH2F("pass_ttbar_"+y, "Passing TTBar "+y, 4, widebins_2018, 40, 0, 200)
			fail_ttbar[y] = TH2F("fail_ttbar_"+y, "Failing TTBar "+y, 4, widebins_2018, 40, 0, 200)
			pass_wgamma[y] = TH2F("pass_wgamma_"+y, "Passing WGamma "+y, 4, widebins_2018, 40, 0, 200)
			fail_wgamma[y] = TH2F("fail_wgamma_"+y, "Failing WGamma "+y, 4, widebins_2018, 40, 0, 200)
			pass_zgamma[y] = TH2F("pass_zgamma_"+y, "Passing ZGamma "+y, 4, widebins_2018, 40, 0, 200)
			fail_zgamma[y] = TH2F("fail_zgamma_"+y, "Failing ZGamma "+y, 4, widebins_2018, 40, 0, 200)
			pass_tbkg[y] = TH2F("pass_tbkg_"+y, "Passing Total-Background "+y, 4, widebins_2018, 40, 0, 200)
			fail_tbkg[y] = TH2F("fail_tbkg_"+y, "Failing Total-Background "+y, 4, widebins_2018, 40, 0, 200)

		mpt_TF[y].SetXTitle("Jet pT")
		mpt_TF[y].SetYTitle("Softdrop Mass")
	
		pass_data[y].SetXTitle("Jet pT")
		pass_data[y].SetYTitle("Softdrop Mass")
		fail_data[y].SetXTitle("Jet pT")
		fail_data[y].SetYTitle("Softdrop Mass")
	
		pass_nres[y].SetXTitle("Jet pT")
		pass_nres[y].SetYTitle("Softdrop Mass")
		fail_nres[y].SetXTitle("Jet pT")
		fail_nres[y].SetYTitle("Softdrop Mass")

		pass_ttbar[y].SetXTitle("Jet pT")
		pass_ttbar[y].SetYTitle("Softdrop Mass")
		fail_ttbar[y].SetXTitle("Jet pT")
		fail_ttbar[y].SetYTitle("Softdrop Mass")

		pass_wgamma[y].SetXTitle("Jet pT")
		pass_wgamma[y].SetYTitle("Softdrop Mass")
		fail_wgamma[y].SetXTitle("Jet pT")
		fail_wgamma[y].SetYTitle("Softdrop Mass")
	
		pass_zgamma[y].SetXTitle("Jet pT")
		pass_zgamma[y].SetYTitle("Softdrop Mass")
		fail_zgamma[y].SetXTitle("Jet pT")
		fail_zgamma[y].SetYTitle("Softdrop Mass")

		pass_tbkg[y].SetXTitle("Jet pT")
		pass_tbkg[y].SetYTitle("Softdrop Mass")
		fail_tbkg[y].SetXTitle("Jet pT")
		fail_tbkg[y].SetYTitle("Softdrop Mass")
	for y in years:
		for pf in ["pass", "fail"]:
			for a in range(0, 4):
				print("ptbin #"+str(a))
				NAME = "bin_"+y+"_ptbin"+str(a)+pf
				refF = ROOT.TFile("../FitHist_"+y+".root")
				refH = refF.Get("Data_"+pf+"_soft")
				for P in ["fit_b"]:
					cDATA = CbnF.Get("shapes_"+P+"/"+NAME+"/data")
					cNRES = CbnF.Get("shapes_"+P+"/"+NAME+"/NonRes"+y)	#Non-Resonant Background
					cTT = CbnF.Get("shapes_"+P+"/"+NAME+"/TTBar")
					cTBKG = CbnF.Get("shapes_"+P+"/"+NAME+"/total_background")
                                        cWG = CbnF.Get("shapes_"+P+"/"+NAME+"/WGamma")
                                        cZG = CbnF.Get("shapes_"+P+"/"+NAME+"/ZGamma")


					cDATA = convertAsymGraph(cDATA, cTT, "data"+P)
					Hvec = []

					for i in [cDATA, cTT, cTBKG, cWG, cZG, cNRES]:
						i.Scale(5.0)
						Hvec.append(convertBinNHist(i, refH, i.GetName()+"new"+P))

					if pf == "pass":
						pfname = "Passing"
					if pf == "fail":
						pfname = "Failing"

					Hvec[0].SetTitle("Data "+pfname+" Softdrop Mass")
					Hvec[0].SetXTitle("Softdrop Mass")
					Hvec[1].SetTitle("TTBar "+pfname+" Softdrop Mass")
					Hvec[1].SetXTitle("Softdrop Mass")
					Hvec[2].SetTitle("Total Background "+pfname+" Softdrop Mass")
					Hvec[2].SetXTitle("Softdrop Mass")
					Hvec[3].SetTitle("W+Gamma "+pfname+" Softdrop Mass")
					Hvec[3].SetXTitle("Softdrop Mass")
					Hvec[4].SetTitle("Z+Gamma "+pfname+" Softdrop Mass")
					Hvec[4].SetXTitle("Softdrop Mass")
					Hvec[5].SetTitle("Non-Resonant "+pfname+" Softdrop Mass")
					Hvec[5].SetXTitle("Softdrop Mass")
					Hvec[0].Draw("hist")
					c0.SaveAs("./plots/"+NAME+"_data_"+y+".png")
					Hvec[1].Draw("hist")
					c0.SaveAs("./plots/"+NAME+"_ttbar_"+y+".png")
					Hvec[2].Draw("hist")
					c0.SaveAs("./plots/"+NAME+"_tbkg_"+y+".png")
					Hvec[3].Draw("hist")
					c0.SaveAs("./plots/"+NAME+"_wgamma_"+y+".png")
					Hvec[4].Draw("hist")
					c0.SaveAs("./plots/"+NAME+"_zgamma_"+y+".png")
					Hvec[5].Draw("hist")
					c0.SaveAs("./plots/"+NAME+"_nres_"+y+".png")
#					savename = P+"_sdm_"+str(ptBins[a])+"-"+str(ptBins[a+1])
#					ptbin = str(ptBins[a])+"-"+str(ptBins[a+1])
#					data = Hvec[0]
#					GoodPlotFormat(data, "markers", ROOT.kBlack, 20)
#					bkg = Hvec[2]
#					GoodPlotFormat(bkg,"thickline", ROOT.kBlue, 1)
#					ttbar = Hvec[1]
#					GoodPlotFormat(ttbar,"thickline", ROOT.kRed, 1)
##					sig = Hvec[3]
#					GoodPlotFormat(sig, "fill", ROOT.kGreen+1, 3003)
    #                                    WG = Hvec[4]
   #                                     GoodPlotFormat(WG,"thickline", ROOT.kViolet, 1)
  #                                      ZG = Hvec[5]
 #                                       GoodPlotFormat(ZG,"thickline", ROOT.kOrange, 1)
#					nonRes = Hvec[6]
#                                        GoodPlotFormat(nonRes,"dashed", ROOT.kCyan, 1)


					if pf == "pass":
						for i in range(1, pass_data[y].GetNbinsY()+1):
							pass_data[y].SetBinContent(a+1, i, Hvec[0].GetBinContent(i))
							if Hvec[1].GetBinContent(i) >= 1.0:
								pass_ttbar[y].SetBinContent(a+1, i, Hvec[1].GetBinContent(i))
							if Hvec[2].GetBinContent(i) >= 1.0:
								pass_tbkg[y].SetBinContent(a+1, i, Hvec[2].GetBinContent(i))
							pass_wgamma[y].SetBinContent(a+1, i, Hvec[3].GetBinContent(i))
							pass_zgamma[y].SetBinContent(a+1, i, Hvec[4].GetBinContent(i))
							if Hvec[5].GetBinContent(i) >= 1.0:
								pass_nres[y].SetBinContent(a+1, i, Hvec[5].GetBinContent(i))
	#						print("xbin: "+str(a+1)+" ybin: "+str(i)+" pass_nres[y] content: "+str(Hvec[5].GetBinContent(i)))
					if pf == "fail":
						for i in range(1, fail_data[y].GetNbinsY()+1):
							fail_data[y].SetBinContent(a+1, i, Hvec[0].GetBinContent(i))
							if Hvec[1].GetBinContent(i) >= 1.0:
								fail_ttbar[y].SetBinContent(a+1, i, Hvec[1].GetBinContent(i))
							if Hvec[2].GetBinContent(i) >= 1.0:
								fail_tbkg[y].SetBinContent(a+1, i, Hvec[2].GetBinContent(i))
							fail_wgamma[y].SetBinContent(a+1, i, Hvec[3].GetBinContent(i))
							fail_zgamma[y].SetBinContent(a+1, i, Hvec[4].GetBinContent(i))
							if Hvec[5].GetBinContent(i) >= 1.0:
								fail_nres[y].SetBinContent(a+1, i, Hvec[5].GetBinContent(i))

					
        	ofile.WriteObject(pass_data[y], "pass_data_"+y)
        	ofile.WriteObject(fail_data[y], "fail_data_"+y)
        	ofile.WriteObject(pass_tbkg[y], "pass_tbkg_"+y)
        	ofile.WriteObject(fail_tbkg[y], "fail_tbkg_"+y)
        	ofile.WriteObject(pass_nres[y], "pass_nres_"+y)
        	ofile.WriteObject(fail_nres[y], "fail_nres_"+y)
        	ofile.WriteObject(pass_ttbar[y], "pass_ttbar_"+y)
        	ofile.WriteObject(fail_ttbar[y], "fail_ttbar_"+y)
        	ofile.WriteObject(pass_wgamma[y], "pass_wgamma_"+y)
        	ofile.WriteObject(fail_wgamma[y], "fail_wgamma_"+y)
        	ofile.WriteObject(pass_zgamma[y], "pass_zgamma_"+y)
		ofile.WriteObject(fail_zgamma[y], "fail_zgamma_"+y)
	
		for i in range(1, mpt_TF[y].GetNbinsX()+1):
			for j in range(1, mpt_TF[y].GetNbinsY()+1):
				if pass_nres[y].GetBinContent(i,j) != 0 and fail_nres[y].GetBinContent(i,j) != 0:
					mpt_TF[y].SetBinContent(i, j, pass_nres[y].GetBinContent(i,j)/fail_nres[y].GetBinContent(i,j))
		ofile.WriteObject(mpt_TF[y], "mpt_TF[y]")
		for b in range(1, mpt_TF[y].GetNbinsX()+1):
			mpt_TF_proj = mpt_TF[y].ProjectionY("mpt_TF_"+str(b-1), b, b)
			mpt_TF_proj.SetTitle(y+" Jet pT Bin #"+str(b-1)+" Transfer Factor")
			mpt_TF_proj.SetXTitle("Softdrop Mass")
			mpt_TF_proj.SetYTitle("Transfer Factor")
			mpt_TF_proj.Draw("hist")
			c0.SaveAs("./plots/ptbin"+str(b-1)+"TF_"+y+".png")
			ofile.WriteObject(mpt_TF_proj, "mpt_TF_proj"+str(b-1))
	
		c0_1 = TCanvas()
		c0_1.cd()
		l0_1 =  ROOT.TLegend(0.58,0.69,0.895,0.89)
		l0_1.SetFillColor(0)
		l0_1.SetLineColor(0)
		for b in range(1, mpt_TF[y].GetNbinsX()+1):
			mpt_TF_proj = mpt_TF[y].ProjectionY("mpt_TF_"+str(b-1), b, b)
			l0_1.AddEntry(mpt_TF_proj, "Jet pT Bin #"+str(b-1)+" Transfer Factor")
			mpt_TF_proj.SetLineColor(b+1)
			mpt_TF_proj.SetTitle(y+" Transfer Factor: All pT Bins")
			mpt_TF_proj.SetXTitle("Softdrop Mass")
			mpt_TF_proj.SetYTitle("Transfer Factor")
			if b == 1:
				mpt_TF_proj.GetYaxis().SetRangeUser(0, 0.085)
				mpt_TF_proj.Draw("hist")
			else:
				mpt_TF_proj.Draw("hist same")
		l0_1.Draw("same")
		c0_1.SaveAs("./plots/allTF_"+y+".png")
		c0_1.Close()

		c1 = TCanvas()
		c1.cd()
		pass_data[y].Draw("COLZ")
		c1.SaveAs("./plots/pass_data_COLZ_"+y+".png")
		pass_data[y].Draw("LEGO2Z")
		c1.SaveAs("./plots/pass_data_LEGO2Z_"+y+".png")
	
		fail_data[y].Draw("COLZ")
		c1.SaveAs("./plots/fail_data_COLZ_"+y+".png")
		fail_data[y].Draw("LEGO2Z")
		c1.SaveAs("./plots/fail_data_LEGO2Z_"+y+".png")
		
		pass_tbkg[y].Draw("COLZ")
		c1.SaveAs("./plots/pass_tbkg_COLZ_"+y+".png")
		pass_tbkg[y].Draw("LEGO2Z")
		c1.SaveAs("./plots/pass_tbkg_LEGO2Z_"+y+".png")
		
		fail_tbkg[y].Draw("COLZ")
		c1.SaveAs("./plots/fail_tbkg_COLZ_"+y+".png")
		fail_tbkg[y].Draw("LEGO2Z")
		c1.SaveAs("./plots/fail_tbkg_LEGO2Z_"+y+".png")
		
		pass_nres[y].Draw("COLZ")
		c1.SaveAs("./plots/pass_nres_COLZ_"+y+".png")
		pass_nres[y].Draw("LEGO2Z")
		c1.SaveAs("./plots/pass_nres_LEGO2Z_"+y+".png")
		
		fail_nres[y].Draw("COLZ")
		c1.SaveAs("./plots/fail_nres_COLZ_"+y+".png")
		fail_nres[y].Draw("LEGO2Z")
		c1.SaveAs("./plots/fail_nres_LEGO2Z_"+y+".png")
		
		pass_ttbar[y].Draw("COLZ")
		c1.SaveAs("./plots/pass_ttbar_COLZ_"+y+".png")
		pass_ttbar[y].Draw("LEGO2Z")
		c1.SaveAs("./plots/pass_ttbar_LEGO2Z_"+y+".png")
		
		fail_ttbar[y].Draw("COLZ")
		c1.SaveAs("./plots/fail_ttbar_COLZ_"+y+".png")
		fail_ttbar[y].Draw("LEGO2Z")
		c1.SaveAs("./plots/fail_ttbar_LEGO2Z_"+y+".png")
		
		pass_wgamma[y].Draw("COLZ")
		c1.SaveAs("./plots/pass_wgamma_COLZ_"+y+".png")
		pass_wgamma[y].Draw("LEGO2Z")
		c1.SaveAs("./plots/pass_wgamma_LEGO2Z_"+y+".png")
		
		fail_wgamma[y].Draw("COLZ")
		c1.SaveAs("./plots/fail_wgamma_COLZ_"+y+".png")
		fail_wgamma[y].Draw("LEGO2Z")
		c1.SaveAs("./plots/fail_wgamma_LEGO2Z_"+y+".png")
		
		pass_zgamma[y].Draw("COLZ")
		c1.SaveAs("./plots/pass_zgamma_COLZ_"+y+".png")
		pass_zgamma[y].Draw("LEGO2Z")
		c1.SaveAs("./plots/pass_zgamma_LEGO2Z_"+y+".png")
		
		fail_zgamma[y].Draw("COLZ")
		c1.SaveAs("./plots/fail_zgamma_COLZ_"+y+".png")
		fail_zgamma[y].Draw("LEGO2Z")
		c1.SaveAs("./plots/fail_zgamma_LEGO2Z_"+y+".png")

		mpt_TF[y].SetTitle(y+" pT vs. Softdrop Mass Transfer Factor") 
		
		mpt_TF[y].Draw("COLZ")
		c1.SaveAs("./plots/mpt_TF_COLZ_"+y+".png")
		mpt_TF[y].Draw("LEGO2Z")
		c1.SaveAs("./plots/mpt_TF_LEGO2Z_"+y+".png")
		mpt_TF[y].Draw("LEGO")
		c1.SaveAs("./plots/mpt_TF_LEGO_"+y+".png")
		mpt_TF[y].Draw("SURF")
		c1.SaveAs("./plots/mpt_TF_SURF2Z_"+y+".png")
		mpt_TF[y].Draw("SURF2Z")
		c1.SaveAs("./plots/mpt_TF_SURF_"+y+".png")
		c1.Close()
