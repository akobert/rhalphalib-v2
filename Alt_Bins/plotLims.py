import ROOT
from ROOT import *
import os
from array import array
import math
from math import *
import sys
import glob
import csv
import ctypes
from ctypes import *
import XRootD
from pyxrootd import client
import numpy as np
from optparse import OptionParser

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
	cmsText     = "CMS Preliminary"
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
	#latex.DrawLatex(0.1265, 0.825, cmsText)
	latex.DrawLatex(0.185, 0.845, cmsText)
	pad.Update()
def setDict(year):
	sample_xsec = {} # cross section used to normalize sample (in combine), currently same as theory
	theory_xsec = {} # real theory cross section to compare (no BR)
	legend_entry = {}
	if year == "Combined" or year == "Combined_ParticleNet" or year == "Combined_ParticleNet_All" or year == "Combined_ParticleNet_Allv2" or year == "Combined_ParticleNet_Allv3" or year == "Combined_ParticleNet_Allv4" or year == "Combined_ParticleNet_Allv5" or year == "Combined_ParticleNet_All_200GeV" or year == "Combined_ParticleNet_All_Retrain" or year == "total_5bin" or year == "total_4bin" or year == "pass_5bin" or year == "pass_4bin":	# Still need to update these event counts
		theory_xsec['M10'] = 25390.0/501525.0
		theory_xsec['M20'] = 15660.0/493822.0
		theory_xsec['M25'] = 13290.0/496157.0
		theory_xsec['M50'] = 8049.0/492039.0
		theory_xsec['M75'] = 5843.0/496487.0
		theory_xsec['M100'] = 4780.0/473228.0
		theory_xsec['M125'] = 3671.0/492822.0
		theory_xsec['M150'] = 3006.0/494347.0

		sample_xsec['M10'] = 25390.0/501525.0
		sample_xsec['M20'] = 15660.0/493822.0
		sample_xsec['M25'] = 13290.0/496157.0
		sample_xsec['M50'] = 8049.0/492039.0
		sample_xsec['M75'] = 5843.0/496487.0
		sample_xsec['M100'] = 4780.0/473228.0
		sample_xsec['M125'] = 3671.0/492822.0
		sample_xsec['M150'] = 3006.0/494347.0

	if year == "2018" or year == "2018_ParticleNet" or year == "2018_ParticleNet_Retrain":
		theory_xsec['M10'] = 25390.0/501525.0
		theory_xsec['M20'] = 15660.0/493822.0
		theory_xsec['M25'] = 13290.0/496157.0
		theory_xsec['M50'] = 8049.0/492039.0
		theory_xsec['M75'] = 5843.0/496487.0
		theory_xsec['M100'] = 4780.0/473228.0
		theory_xsec['M125'] = 3671.0/492822.0
		theory_xsec['M150'] = 3006.0/494347.0

		sample_xsec['M10'] = 25390.0/501525.0
		sample_xsec['M20'] = 15660.0/493822.0
		sample_xsec['M25'] = 13290.0/496157.0
		sample_xsec['M50'] = 8049.0/492039.0
		sample_xsec['M75'] = 5843.0/496487.0
		sample_xsec['M100'] = 4780.0/473228.0
		sample_xsec['M125'] = 3671.0/492822.0
		sample_xsec['M150'] = 3006.0/494347.0


	elif year == "2017" or year == "2017_stable" or year == "2017_ParticleNet":
		theory_xsec['M10'] = 25390.0/499072.0
		theory_xsec['M20'] = 15660.0/487317.0
		theory_xsec['M25'] = 13290.0/495785.0
		theory_xsec['M50'] = 8049.0/496133.0
		theory_xsec['M75'] = 5843.0/496846.0
		theory_xsec['M100'] = 4780.0/490390.0
		theory_xsec['M125'] = 3671.0/496239.0
		theory_xsec['M150'] = 3006.0/496687.0
		
		sample_xsec['M10'] = 25390.0/499072.0
		sample_xsec['M20'] = 15660.0/487317.0
		sample_xsec['M25'] = 13290.0/495785.0
		sample_xsec['M50'] = 8049.0/496133.0
		sample_xsec['M75'] = 5843.0/496846.0
		sample_xsec['M100'] = 4780.0/490390.0
		sample_xsec['M125'] = 3671.0/496239.0
		sample_xsec['M150'] = 3006.0/496687.0

	elif year == "2016" or year == "2016_ParticleNet":	
#		theory_xsec['M10'] = 25390.0/254454.0
#		theory_xsec['M20'] = 15660.0/245574.0
#		theory_xsec['M25'] = 13290.0/251926.0
#		theory_xsec['M50'] = 8049.0/251364.0
#		theory_xsec['M75'] = 5843.0/250935.0
#		theory_xsec['M100'] = 4780.0/253395.0
#		theory_xsec['M125'] = 3671.0/242353.0
#		theory_xsec['M150'] = 3006.0/256233.0
		
#		sample_xsec['M10'] = 25390.0/254454.0
#		sample_xsec['M20'] = 15660.0/245574.0
#		sample_xsec['M25'] = 13290.0/251926.0
#		sample_xsec['M50'] = 8049.0/251364.0
#		sample_xsec['M75'] = 5843.0/250935.0
#		sample_xsec['M100'] = 4780.0/253395.0
#		sample_xsec['M125'] = 3671.0/242353.0
#		sample_xsec['M150'] = 3006.0/256233.0
		
		theory_xsec['M10'] = 25390.0/512198.0
		theory_xsec['M20'] = 15660.0/496832.0
		theory_xsec['M25'] = 13290.0/506084.0
		theory_xsec['M50'] = 8049.0/506043.0
		theory_xsec['M75'] = 5843.0/498106.0
		theory_xsec['M100'] = 4780.0/500455.0
		theory_xsec['M125'] = 3671.0/492694.0
		theory_xsec['M150'] = 3006.0/503417.0
		
		sample_xsec['M10'] = 25390.0/512198.0
		sample_xsec['M20'] = 15660.0/496832.0
		sample_xsec['M25'] = 13290.0/506084.0
		sample_xsec['M50'] = 8049.0/506043.0
		sample_xsec['M75'] = 5843.0/498106.0
		sample_xsec['M100'] = 4780.0/500455.0
		sample_xsec['M125'] = 3671.0/492694.0
		sample_xsec['M150'] = 3006.0/503417.0
		
#	elif year == "2016APV":	
#		theory_xsec['M10'] = 25390.0/257744.0
#		theory_xsec['M20'] = 15660.0/251258.0
#		theory_xsec['M25'] = 13290.0/254158.0
#		theory_xsec['M50'] = 8049.0/254679.0
#		theory_xsec['M75'] = 5843.0/247171.0
#		theory_xsec['M100'] = 4780.0/247060.0
#		theory_xsec['M125'] = 3671.0/250341.0
#		theory_xsec['M150'] = 3006.0/247184.0
		
#		sample_xsec['M10'] = 25390.0/257744.0
#		sample_xsec['M20'] = 15660.0/251258.0
#		sample_xsec['M25'] = 13290.0/254158.0
#		sample_xsec['M50'] = 8049.0/254679.0
#		sample_xsec['M75'] = 5843.0/247171.0
#		sample_xsec['M100'] = 4780.0/247060.0
#		sample_xsec['M125'] = 3671.0/250341.0
#		sample_xsec['M150'] = 3006.0/247184.0

	return theory_xsec, sample_xsec

def plotLimits(year, poly1, masses, files, scale):
	theory_xsec, sample_xsec = setDict(year)
	
	print(theory_xsec)
	N = len(masses)
	print("Number of Mass Points: "+str(N))

#	yellowList, greenList, medianList, obsList = [], [], [], []
	yellowList = ROOT.TGraph(2*N)
	greenList = ROOT.TGraph(2*N)
	medianList = ROOT.TGraph(N)
	obsList = ROOT.TGraph(N)
	
	i = -1
	for mass in masses:
#		yellow = ROOT.TGraph(2*N)    # yellow band
#		green = ROOT.TGraph(2*N)     # green band
#		median = ROOT.TGraph(N)      # median line
#		obs = ROOT.TGraph(N)         # observed

		i += 1
		
		tfile = ROOT.TFile.Open(files["M"+mass], 'read')
		tree = tfile.Get("limit")

		limits = []
		for quantile in tree:
			limits.append(tree.limit)
		print("Mass: "+str(mass))
		print("Theory XSection: "+str(theory_xsec["M"+mass]))

		
		# Coupling Limit
#		yellow.SetPoint(i, float(mass), math.sqrt(limits[4] * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))) # +2 sigma
#		green.SetPoint(i, float(mass), math.sqrt(limits[3] * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))) # +1 sigma
#		median.SetPoint(i, float(mass), math.sqrt(limits[2] * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))) # median
#		green.SetPoint(2*N-1-i, float(mass), math.sqrt(limits[1] * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))) # -1 sigma
#		yellow.SetPoint(2*N-1-i, float(mass), math.sqrt(limits[0] * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))) # -2 sigma
		yellowList.SetPoint(i, float(mass), math.sqrt(limits[4] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))) # +2 sigma
		greenList.SetPoint(i, float(mass), math.sqrt(limits[3] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))) # +1 sigma
		medianList.SetPoint(i, float(mass), math.sqrt(limits[2] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))) # median
		greenList.SetPoint(2*N-1-i, float(mass), math.sqrt(limits[1] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))) # -1 sigma
		yellowList.SetPoint(2*N-1-i, float(mass), math.sqrt(limits[0] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))) # -2 sigma
	
		print("+2 Sigma Limit: "+str(math.sqrt(limits[4] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))))
		print("+1 Sigma Limit: "+str(math.sqrt(limits[3] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))))
		print("Median Limit: "+str(math.sqrt(limits[2] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))))
		print("-1 Sigma Limit: "+str(math.sqrt(limits[1] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))))
		print("-2 Sigma Limit: "+str(math.sqrt(limits[0] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))))
			
		if len(limits)>5:
			obsList.SetPoint(i, float(mass), math.sqrt(limits[5] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))
			print("Observed Limit: "+str(math.sqrt(limits[5] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))))
		#print "observed (expected) @ %s: %s (%s)"%( float(mass), math.sqrt(limits[5]*fac/theory), math.sqrt(limits[2]*sample_xsec["M"+mass]/theory_xsec["M"+mass]))
		
#		yellowList.append(yellow)
#		greenList.append(green)
#		medianList.append(median)
#		obsList.append(obs)
		
	return yellowList, greenList, medianList, obsList
	
	
if __name__ == "__main__":
#	parser = OptionParser()
	year = sys.argv[1]
	poly1 = sys.argv[2] #2 for low, 3 for mid (no added tag), 4 for high, 5 for veryhigh
	poly2 = sys.argv[3]
	order = ""	

	if poly1 == "2":
		order = "_low"
	elif poly1 == "3":
		order = ""
	elif poly1 == "4":
		order = "_high"
	elif poly1 == "5":
		order = "_veryhigh"

	print("Order "+order)
	
#	if year == "2016APV":
#		LUMI = 1.95
	if year == "2016" or year == "2016_ParticleNet":
#		LUMI = 1.681
		LUMI = 3.631
	if year == "2017" or year == "2017_stable" or year == "2017_ParticleNet":
		LUMI = 4.148
	if year == "2018" or year == "2018_ParticleNet" or year == "2018_ParticleNet_Retrain":
		LUMI = 5.982
	if year == "Combined":
		LUMI = 13.761
	if year == "Combined_ParticleNet":
		LUMI = 7.779
	if year == "Combined_ParticleNet_All":
		LUMI = 13.761
	if year == "Combined_ParticleNet_Allv2":
		LUMI = 13.761
	if year == "Combined_ParticleNet_Allv3":
		LUMI = 13.761
	if year == "Combined_ParticleNet_Allv4":
		LUMI = 13.761
	if year == "Combined_ParticleNet_Allv5":
		LUMI = 13.761
	if year == "total_5bin" or year == "total_4bin" or year == "pass_5bin" or year == "pass_4bin":
		LUMI = 13.761
	if year == "Combined_ParticleNet_All_200GeV":
		LUMI = 13.761
	if year == "Combined_ParticleNet_All_Retrain":
		LUMI = 13.761

	# LUMI Scaling
	scaled_lumi = sys.argv[4]
	scale = 1.0

	if scaled_lumi != "0":
		scale = 1.0/math.sqrt(float(scaled_lumi)/LUMI)
		LUMI = float(scaled_lumi)
		print("Scaling LUMI to "+str(LUMI))
	else:
		scale = 1.0

	print("Lumi Scaling: "+str(math.sqrt(scale)))

	masses = ["10", "20", "25", "50", "75", "100", "125", "150"]

	files = {}

	for mass in masses:
		files["M"+mass] = "/home/akobert/CMSSW_11_3_4/src/FitTest/Alt_Bins/"+year+"/sig"+mass+""+order+"/higgsCombineTest.AsymptoticLimits.mH120.root"	
		
	yellowList, greenList, medianList, obsList = plotLimits(year, poly1, masses, files, scale)
	
	W = 800
	H  = 600
	T = 0.08*H
	B = 0.12*H
	L = 0.12*W#*scaleleftmargin
	R = 0.04*W#*scalerightmargin
	c = ROOT.TCanvas("c","c",100,100,W,H)
	c.SetFillColor(0)
	c.SetBorderMode(0)
	c.SetFrameFillStyle(0)
	c.SetFrameBorderMode(0)
	c.SetLeftMargin( L/W )
	c.SetRightMargin( R/W )
	c.SetTopMargin( T/H )
	c.SetBottomMargin( B/H )
	c.SetTickx(0)
	c.SetTicky(0)
	c.SetLogy() # Setting Log Y to match Marc
	c.SetLogx() # Setting Log X to match Marc
	if scaled_lumi != "0":
		c.SetTitle("Z' Coupling Limit ("+year+") Bernstein Polynomial Order "+str(poly1)+", "+str(poly2)+" (LUMI Scaled)")
		#c.SetTitle("Z' Coupling Limit (LUMI Scaled)")
	else:
		c.SetTitle("Z' Coupling Limit ("+year+") Bernstein Polynomial Order "+str(poly1)+", "+str(poly2))
	c.cd()
	frame = c.DrawFrame(1.4,0.001, 4.1, 1.0)
	frame.GetYaxis().CenterTitle()
	frame.GetYaxis().SetTitleSize(0.05)
	frame.GetXaxis().SetTitleSize(0.05)
	frame.GetXaxis().SetLabelSize(0.04)
	frame.GetYaxis().SetLabelSize(0.04)
	frame.GetYaxis().SetTitleOffset(0.9)
	frame.GetXaxis().SetNdivisions(508)
	frame.GetYaxis().CenterTitle(True)

	h_limit = ROOT.TMultiGraph()
	#for yellow in yellowList: h_limit.Add(yellow)
	#for green in greenList: h_limit.Add(green)
	#for median in medianList:  h_limit.Add(median)
	h_limit.Add(yellowList)
	h_limit.Add(greenList)
	h_limit.Add(medianList)
	if scaled_lumi != "0":
		h_limit.SetTitle("Z' Coupling Limit ("+year+") Polynomial Order "+str(poly1)+", "+str(poly2)+" (LUMI Scaled)")
#		h_limit.SetTitle("Z' Coupling Limit (LUMI Scaled)")
	else:
		h_limit.SetTitle("Z' Coupling Limit ("+year+") Polynomial Order "+str(poly1)+", "+str(poly2))
		h_limit.Add(obsList) # Add observed limits for NON-LUMI scaled plots
	h_limit.Draw('a3')
	h_limit.GetXaxis().SetLimits(6,500)
#	h_limit.SetMinimum(options.xsecMin)
#	h_limit.SetMaximum(options.xsecMax)
	h_limit.GetXaxis().SetTitle("Z' mass [GeV]")
	h_limit.GetYaxis().SetTitle("g'_{q}")
#	h_limit.GetYaxis().SetMoreLogLabels()
#	h_limit.GetYaxis().SetNdivisions(10);
#	h_limit.GetYaxis().SetLimits(0.02,1.2)
#	h_limit.GetYaxis().SetLimits(0.02,1.2)
	h_limit.SetMinimum(0.03)
	h_limit.SetMaximum(1.4)
	#h_limit.GetYaxis().SetNoExponent()
#	h_limit.GetXaxis().SetMoreLogLabels()
 #       h_limit.GetXaxis().SetNdivisions(10);
#	h_limit.GetXaxis().SetNoExponent()
	h_limit.GetYaxis().SetTitleOffset(0.9)

	
#	for yellow in yellowList:
#		yellow.SetFillColor(ROOT.kOrange)
#		yellow.SetLineColor(ROOT.kBlack)
#		yellow.SetFillStyle(1001)
#		yellow.SetLineWidth(2)
#		yellow.SetLineStyle(2)
#		yellow.Draw('Fsame')

#	for green in greenList:    
#		green.SetFillColor(ROOT.kGreen+1)
#		green.SetLineColor(ROOT.kBlack)
#		green.SetLineWidth(2)
#		green.SetLineStyle(2)
#		green.SetFillStyle(1001)
#		green.Draw('Fsame')
		
#	for median in medianList:
#		median.SetLineColor(1)
#		median.SetLineWidth(3)
#		median.SetLineStyle(2)
#		median.Draw('Csame')
#   		print(median)
 
#	for obs in obsList:
#		obs.SetMarkerStyle(20)
#		obs.SetLineWidth(3)
#		if options.observed:
#			#obs.Draw('PLsame')
#			obs.Draw('Csame')
		
	yellowList.SetFillColor(ROOT.kOrange)
	yellowList.SetLineColor(ROOT.kBlack)
	yellowList.SetFillStyle(1001)
	yellowList.SetLineWidth(2)
	yellowList.SetLineStyle(2)
	yellowList.Draw('Fsame')

	greenList.SetFillColor(ROOT.kGreen+1)
	greenList.SetLineColor(ROOT.kBlack)
	greenList.SetLineWidth(2)
	greenList.SetLineStyle(2)
	greenList.SetFillStyle(1001)
	greenList.Draw('Fsame')
		
	medianList.SetLineColor(1)
	medianList.SetLineWidth(3)
	medianList.SetLineStyle(2)
#	medianList.SetLineStyle(1)
	if scaled_lumi == "0":
		medianList.Draw('Lsame')
	else:
		print("Lumi Scaling Applied")
		medianList.SetMarkerStyle(20)
		medianList.SetMarkerSize(1.5)
		medianList.Draw('PLsame')
   	print("Expected Limit: ")
   	print(medianList)
 
	if scaled_lumi == "0":
		obsList.SetLineColor(1)
		obsList.SetLineWidth(3)
		obsList.SetLineStyle(1)
		obsList.SetMarkerStyle(20)
		obsList.SetMarkerSize(1.5)
		obsList.Draw('LPsame')
   		print("Observed Limit: ")
   		print(obsList)
		

#	obsList.SetMarkerStyle(20)
#	obsList.SetLineWidth(3)
#	if options.observed:
#		#obsList.Draw('PLsame')
#		obsList.Draw('Csame')

	
	
	AddCMSLumi(ROOT.gPad, LUMI, "Preliminary")
	
	legend = ROOT.TLegend(0.35, 0.65, 0.45, 0.85)
	legend.SetFillStyle(0)
	legend.SetBorderSize(0)
	legend.SetTextSize(0.038)
	legend.SetTextFont(42)
#	if options.observed:
#		legend.AddEntry(obs, "Observed",'l')
	#legend.AddEntry(medianList, "Asymptotic CL_{s} expected",'L')
	legend.AddEntry(medianList, "Asymptotic CL_{s} expected",'CL')
	legend.AddEntry(greenList, "Expected #pm 1 s.d.",'lf')
	legend.AddEntry(yellowList,"Expected #pm 2 s.d.",'lf')
	if scaled_lumi == "0":
		legend.AddEntry(obsList, "Observed Limits", 'CL')
	legend.Draw()
	
	legend.Draw("same")
	
	if scaled_lumi != "0":
		c.SaveAs("/home/akobert/CMSSW_11_3_4/src/FitTest/Alt_Bins/limits/limits_"+year+""+order+"_scaled.png")
	else:
		c.SaveAs("/home/akobert/CMSSW_11_3_4/src/FitTest/Alt_Bins/limits/limits_"+year+""+order+".png")
	
	
	
