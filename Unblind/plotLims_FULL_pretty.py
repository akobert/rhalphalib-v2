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


mq = {"u":2.3e-3, "d":4.8e-3, "c":1.275, "s":0.095, "b":4.18, "t":173.21}


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
	lumiTextSize     = 0.71
	lumiTextOffset   = 0.165
	cmsTextSize      = 0.71
	cmsTextOffset    = 0.165
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
	#latex.DrawLatex(0.185, 0.845, cmsText)
	latex.DrawLatex(0.13, 1-t+lumiTextOffset*t, cmsText)
	pad.Update()

# Taken from summary plot code https://gitlab.cern.ch/cms-exo-mci/exo-dmsummaryplots/-/blob/master/dmplots/gqplot/zprime_equations.py
def gq_Z_constraint(x, k):	
	# Zwidth = 2.4952;
	# ZwidthError = 0.0023*3; # times 3 to give 3 sigma
	# relZwidthUnc = ZwidthError/Zwidth;
	# sin2thetaW = 0.2312;
	# sinthetaW = math.sqrt(sin2thetaW);
	# costhetaW = math.sqrt(1 - sin2thetaW);
	# mW = 80.385;
	# vev = 246.;
	# g = mW*2./vev;
	# Vu = 0.25 - (4. * sin2thetaW / 6.);
	# Vd = -0.25 - (2. * sin2thetaW / 6.);
	# mZ = 91.18;
	# mZp = x[0];

	# # y = gZ
	# ynum = relZwidthUnc * 3 * g * -1. * math.fabs(1-(mZp*mZp/(mZ*mZ))) * (2*Vu*Vu + 3*Vd*Vd + 5/16.)
	# yden = 2*0.01*costhetaW*sinthetaW*(2*Vu+3*Vd);
	# # print ynum,yden,x[0],math.sqrt(ynum/yden)
	# y = math.sqrt(ynum/yden);
	# y *= 1.5;

	mZ = 91.1876
	mZp = x[0]
	ynum = 4. * math.sqrt( 4. * math.pi ) * 1.96 * 1.1e-3 * ( 1-(mZp*mZp/(mZ*mZ)) )
	yden = 1.193 * 0.02
	if ynum < 0: 
		ynum *= -1.
	y = math.sqrt(ynum/yden) / 6.

	return y

def gq_upsilon_constraint(x, k):
    mUps = 9.46
    WidUps = 54e-6
    alpha = 1/130.
    Argus = 5.3e-2
    DimuUps = 2.48e-2
    RUps = Argus/DimuUps
    mZp = x[0]

    y = abs(( 2*4*math.pi*alpha ) * ( math.sqrt( 3.*RUps + 1 ) + 1 ) * (1-(mZp*mZp/(mUps*mUps))))

    return math.sqrt(y)/6.


# Convert Z' width (Gamma(Z' -> qq)) to g_q
def width_to_gq(width, mZp, vtype):
    den = 0.
    for qtype in ["u", "d", "s", "c", "b", "t"]:
        if mZp > 2. * mq[qtype]:
            if vtype == "vector":
                den += 1. / (4. * math.pi) * (1. - 4. * mq[qtype]**2 / mZp**2)**0.5 * (1. + 2. * mq[qtype]**2 / mZp**2)
            elif vtype == "axial":
                den += 1. / (4. * math.pi) * (1. - 4. * mq[qtype]**2 / mZp**2)**1.5
            else:
                print(("[width_to_gq] ERROR : Unknown vtype {}".format(vtype)))
                sys.exit(1)
    return math.sqrt(width / mZp / den)


# Convert Gamma-over-mass (Gamma(Z'->qq)) to g_q
def gom_to_gq(gom, mZp, vtype):
	return width_to_gq(gom * mZp, mZp, vtype)


def setDict(year):
	sample_xsec = {} # cross section used to normalize sample (in combine), currently same as theory
	theory_xsec = {} # real theory cross section to compare (no BR)
	legend_entry = {}
	if year == "Combined" or year == "Combined_ParticleNet" or year == "Combined_ParticleNet_All" or year == "Combined_ParticleNet_Allv2" or year == "Combined_ParticleNet_Allv3" or year == "Combined_ParticleNet_Allv4" or year == "Combined_ParticleNet_Allv5" or year == "Combined_ParticleNet_All_200GeV" or year == "Combined_ParticleNet_All_Retrain" or year == "total_5bin" or year == "total_4bin" or year == "pass_5bin" or year == "pass_4bin" or year == "pass_4bin_V2" or year == "pass_4bin_V3" or year == "pass_4bin_V4" or year == "pass_4bin_V5" or year == "Run2" or year == "Run2_23" or year == "Run2_mixed" or year == "Run2_mixed_noTail":	# Still need to update these event counts
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

		# Interpolated samples, need to actually calculate these at some point
		theory_xsec['M15'] = 25390.0/501525.0
		theory_xsec['M30'] = 25390.0/501525.0
		theory_xsec['M35'] = 25390.0/501525.0
		theory_xsec['M40'] = 25390.0/501525.0
		theory_xsec['M45'] = 25390.0/501525.0
		theory_xsec['M55'] = 25390.0/501525.0
		theory_xsec['M60'] = 25390.0/501525.0
		theory_xsec['M65'] = 25390.0/501525.0
		theory_xsec['M70'] = 25390.0/501525.0
		theory_xsec['M80'] = 25390.0/501525.0
		theory_xsec['M85'] = 25390.0/501525.0
		theory_xsec['M90'] = 25390.0/501525.0
		theory_xsec['M95'] = 25390.0/501525.0
		theory_xsec['M105'] = 25390.0/501525.0
		theory_xsec['M110'] = 25390.0/501525.0
		theory_xsec['M115'] = 25390.0/501525.0
		theory_xsec['M120'] = 25390.0/501525.0
		theory_xsec['M130'] = 25390.0/501525.0
		theory_xsec['M135'] = 25390.0/501525.0
		theory_xsec['M140'] = 25390.0/501525.0
		theory_xsec['M145'] = 25390.0/501525.0
		
		sample_xsec['M15'] = 25390.0/501525.0
		sample_xsec['M30'] = 25390.0/501525.0
		sample_xsec['M35'] = 25390.0/501525.0
		sample_xsec['M40'] = 25390.0/501525.0
		sample_xsec['M45'] = 25390.0/501525.0
		sample_xsec['M55'] = 25390.0/501525.0
		sample_xsec['M60'] = 25390.0/501525.0
		sample_xsec['M65'] = 25390.0/501525.0
		sample_xsec['M70'] = 25390.0/501525.0
		sample_xsec['M80'] = 25390.0/501525.0
		sample_xsec['M85'] = 25390.0/501525.0
		sample_xsec['M90'] = 25390.0/501525.0
		sample_xsec['M95'] = 25390.0/501525.0
		sample_xsec['M105'] = 25390.0/501525.0
		sample_xsec['M110'] = 25390.0/501525.0
		sample_xsec['M115'] = 25390.0/501525.0
		sample_xsec['M120'] = 25390.0/501525.0
		sample_xsec['M130'] = 25390.0/501525.0
		sample_xsec['M135'] = 25390.0/501525.0
		sample_xsec['M140'] = 25390.0/501525.0
		sample_xsec['M145'] = 25390.0/501525.0

	return theory_xsec, sample_xsec

def plotLimits(year, masses, files, scale):
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
		
		print("Sigma spread+1: "+str((math.sqrt(limits[3] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))-(math.sqrt(limits[2] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))))+" Sigma spread-1: "+str((math.sqrt(limits[2] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))-(math.sqrt(limits[1] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))))
		print("Sigma spread1 diff: "+str(((math.sqrt(limits[3] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))-(math.sqrt(limits[2] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))) - ((math.sqrt(limits[2] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))-(math.sqrt(limits[1] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))))))
		print("Sigma spread+2: "+str((math.sqrt(limits[4] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))-(math.sqrt(limits[2] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))))+" Sigma spread-2: "+str((math.sqrt(limits[2] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))-(math.sqrt(limits[0] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))))
		print("Sigma spread2 diff: "+str(((math.sqrt(limits[4] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))-(math.sqrt(limits[2] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))) - ((math.sqrt(limits[2] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))-(math.sqrt(limits[0] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))))))
	
		if len(limits)>5:
			obsList.SetPoint(i, float(mass), math.sqrt(limits[5] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass])))
			print("Observed Limit: "+str(math.sqrt(limits[5] * scale * sample_xsec["M"+mass]/(16.0*theory_xsec["M"+mass]))))
		#print "observed (expected) @ %s: %s (%s)"%( float(mass), math.sqrt(limits[5]*fac/theory), math.sqrt(limits[2]*sample_xsec["M"+mass]/theory_xsec["M"+mass]))
		
#		yellowList.append(yellow)
#		greenList.append(green)
#		medianList.append(median)
#		obsList.append(obs)

        marcObsList = ROOT.TGraph(24)
        marcExpectList = ROOT.TGraph(24)

        marcObsList.SetPoint(0, 10, 0.226046843486)
        marcObsList.SetPoint(1, 15, 0.178831789302)
        marcObsList.SetPoint(2, 20, 0.173823977383)
        marcObsList.SetPoint(3, 25, 0.164140246513)
        marcObsList.SetPoint(4, 30, 0.180024516135)
        marcObsList.SetPoint(5, 35, 0.144003433109)
        marcObsList.SetPoint(6, 40, 0.180287823129)
        marcObsList.SetPoint(7, 45, 0.177674059313)
        marcObsList.SetPoint(8, 50, 0.145894943482)
        marcObsList.SetPoint(9, 55, 0.123247280825)
        marcObsList.SetPoint(10, 60, 0.1275632009)
        marcObsList.SetPoint(11, 65, 0.157364058637)
        marcObsList.SetPoint(12, 70, 0.203703446968)
        marcObsList.SetPoint(13, 75, 0.202804712579)
        marcObsList.SetPoint(14, 80, 0.218142261731)
        marcObsList.SetPoint(15, 85, 0.149504985376)
        marcObsList.SetPoint(16, 90, 0.155719871243)
        marcObsList.SetPoint(17, 95, 0.162523865438)
        marcObsList.SetPoint(18, 100, 0.134095335738)
        marcObsList.SetPoint(19, 105, 0.140555208524)
        marcObsList.SetPoint(20, 110, 0.134774291921)
        marcObsList.SetPoint(21, 115, 0.145405916036)
        marcObsList.SetPoint(22, 120, 0.193564822015)
        marcObsList.SetPoint(23, 125, 0.149224103353)

        marcExpectList.SetPoint(0, 10, 0.22750343404)
        marcExpectList.SetPoint(1, 15, 0.180121004852)
        marcExpectList.SetPoint(2, 20, 0.165687226256)
        marcExpectList.SetPoint(3, 25, 0.139560010489)
        marcExpectList.SetPoint(4, 30, 0.141490392783)
        marcExpectList.SetPoint(5, 35, 0.134814365683)
        marcExpectList.SetPoint(6, 40, 0.144149505241)
        marcExpectList.SetPoint(7, 45, 0.138388856829)
        marcExpectList.SetPoint(8, 50, 0.13661330259)
        marcExpectList.SetPoint(9, 55, 0.144149505241)
        marcExpectList.SetPoint(10, 60, 0.147313912747)
        marcExpectList.SetPoint(11, 65, 0.161037758685)
        marcExpectList.SetPoint(12, 70, 0.166991870228)
        marcExpectList.SetPoint(13, 75, 0.17523545668)
        marcExpectList.SetPoint(14, 80, 0.196265038348)
        marcExpectList.SetPoint(15, 85, 0.204919953884)
        marcExpectList.SetPoint(16, 90, 0.199010137235)
        marcExpectList.SetPoint(17, 95, 0.19515618745)
        marcExpectList.SetPoint(18, 100, 0.172111579603)
        marcExpectList.SetPoint(19, 105, 0.191791173328)
        marcExpectList.SetPoint(20, 110, 0.196265038348)
        marcExpectList.SetPoint(21, 115, 0.182514744458)
        marcExpectList.SetPoint(22, 120, 0.195711398209)
        marcExpectList.SetPoint(23, 125, 0.188653857045)


	return yellowList, greenList, medianList, obsList, marcObsList, marcExpectList
	
	
if __name__ == "__main__":
#	parser = OptionParser()

	# Try style macro
	gROOT.LoadMacro("setTDRStyle_teliko.C")
	setTDRStyle_teliko()	

        marc = False # Draw the Marc limits?

	constraint = True # Draw the constraints?

	GOM = False # Draw the Gamma/Mass constraints?


	year = "Run2"
	scale = 1.0
	
#	if year == "2016APV":
#		LUMI = 1.95
	if year == "2016" or year == "2016_ParticleNet" or year == "2016_mixed":
#		LUMI = 1.681
		LUMI = 36.31
	if year == "2017" or year == "2017_stable" or year == "2017_ParticleNet" or year == "2017_mixed":
		LUMI = 41.48
	if year == "2018" or year == "2018_ParticleNet" or year == "2018_ParticleNet_Retrain" or year == "2018_mixed" or year == "2018_mixed_noBin0" or year == "2018_mixed_noBin1" or year == "2018_mixed_noBin2" or year == "2018_mixed_noBin3":
		LUMI = 59.82
	if year == "Combined" or year == "Run2" or year == "Run2_23" or year == "Run2_mixed" or year == "Run2_mixed_noTail":
		LUMI = 138

	generated_masses = ["10", "20", "25", "50", "75", "100", "125", "150"]
	masses = ["10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", "65", "70", "75", "80", "85", "90", "95", "100", "105", "110", "115", "120", "125", "130", "135", "140", "145", "150"]

	files = {}

	for mass in masses:
		if mass in generated_masses:
			files["M"+mass] = "/home/akobert/CMSSW_11_3_4/src/FitTest/Unblind/Run2_mixed/sig"+mass+"/higgsCombineTest.AsymptoticLimits.mH120.root"	
		else:
			files["M"+mass] = "/home/akobert/CMSSW_11_3_4/src/FitTest/Unblind/Run2_Interpolated/sig"+mass+"/higgsCombineTest.AsymptoticLimits.mH120.root"	
	yellowList, greenList, medianList, obsList, marcObsList, marcExpectList = plotLimits(year, masses, files, scale)
	
	W = 800
	H  = 600
	T = 0.08*H
	B = 0.12*H
	L = 0.12*W#*scaleleftmargin
	R = 0.04*W#*scalerightmargin
	c = ROOT.TCanvas("c","c",100,100,W,H)
#	c.SetFillColor(0)
#	c.SetBorderMode(0)
#	c.SetFrameFillStyle(0)
#	c.SetFrameBorderMode(0)
#	c.SetLeftMargin( L/W )
#	c.SetRightMargin( R/W )
	c.SetTopMargin( T/H )
#	c.SetBottomMargin( B/H )
#	c.SetTickx(0)
#	c.SetTicky(0)
	c.SetLogy() # Setting Log Y to match Marc
#	c.SetLogx() # Setting Log X to match Marc
#	if scaled_lumi != "0":
		#c.SetTitle("Z' Coupling Limit ("+year+") Bernstein Polynomial Order "+str(poly1)+", "+str(poly2)+" (LUMI Scaled)")
		#c.SetTitle("Z' Coupling Limit Bernstein Polynomial Order "+str(poly1)+", "+str(poly2)+" (LUMI Scaled)")
#		c.SetTitle("Z' Coupling Limit (LUMI Scaled)")
#	else:
#		c.SetTitle("Z' Coupling Limit ("+year+") Bernstein Polynomial Order "+str(poly1)+", "+str(poly2))
#	c.SetTitle("Z' Coupling Limit")
	c.cd()
	frame = c.DrawFrame(1.4,0.001, 4.1, 1.0)
	frame.GetYaxis().CenterTitle()
	frame.GetYaxis().SetTitleSize(0.05)
	frame.GetXaxis().SetTitleSize(0.05)
	frame.GetXaxis().SetLabelSize(0.04)
	frame.GetYaxis().SetLabelSize(0.04)
	frame.GetYaxis().SetTitleOffset(0.6)
	frame.GetXaxis().SetNdivisions(508)
	frame.GetYaxis().CenterTitle(True)

	h_limit = ROOT.TMultiGraph()
	#for yellow in yellowList: h_limit.Add(yellow)
	#for green in greenList: h_limit.Add(green)
	#for median in medianList:  h_limit.Add(median)
	h_limit.Add(yellowList)
	h_limit.Add(greenList)
	h_limit.Add(medianList)
#	if scaled_lumi != "0":
		#h_limit.SetTitle("Z' Coupling Limit ("+year+") Polynomial Order "+str(poly1)+", "+str(poly2)+" (LUMI Scaled)")
		#h_limit.SetTitle("Z' Coupling Limit Polynomial Order "+str(poly1)+", "+str(poly2)+" (LUMI Scaled)")
#		h_limit.SetTitle("Z' Coupling Limit (LUMI Scaled)")
#	else:
		#h_limit.SetTitle("Z' Coupling Limit ("+year+") Polynomial Order "+str(poly1)+", "+str(poly2))
		#h_limit.SetTitle("Z' Coupling Limit Polynomial Order "+str(poly1)+", "+str(poly2))
#		h_limit.SetTitle("Z' Coupling Limit")
#		h_limit.Add(obsList) # Add observed limits for NON-LUMI scaled plots

#	h_limit.SetTitle("Z' Coupling Limit")
	h_limit.Add(obsList) # Add observed limits for NON-LUMI scaled plots


	h_limit.Draw('a3')
#	h_limit.GetXaxis().SetLimits(6,500)
	h_limit.GetXaxis().SetLimits(5,155)
#	h_limit.SetMinimum(options.xsecMin)
#	h_limit.SetMaximum(options.xsecMax)
	h_limit.GetXaxis().SetTitle("Z' mass [GeV]")
	h_limit.GetYaxis().SetTitle("g'_{q} Coupling Strength")
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
	h_limit.GetYaxis().SetTitleOffset(0.6)
	h_limit.GetXaxis().SetTitleOffset(1.1)

	
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
	#yellowList.SetLineColor(ROOT.kBlack)
	yellowList.SetLineColor(ROOT.kWhite)
	yellowList.SetFillStyle(1001)
	yellowList.SetLineWidth(1)
	yellowList.SetLineStyle(1)
	yellowList.Draw('Fsame')

	greenList.SetFillColor(ROOT.kGreen+1)
	#greenList.SetLineColor(ROOT.kBlack)
	greenList.SetLineColor(ROOT.kWhite)
	greenList.SetLineWidth(1)
	greenList.SetLineStyle(1)
	greenList.SetFillStyle(1001)
	greenList.Draw('Fsame')
		
	medianList.SetLineColor(1)
	medianList.SetLineWidth(3)
	medianList.SetLineStyle(2)
#	medianList.SetLineStyle(1)
#	if scaled_lumi == "0":
#		medianList.Draw('Lsame')
#	else:
#		print("Lumi Scaling Applied")
#		medianList.SetMarkerStyle(20)
#		medianList.SetMarkerSize(1.5)
#		medianList.Draw('PLsame')
	medianList.Draw('Lsame')
   	print("Expected Limit: ")
   	print(medianList)
 
#	if scaled_lumi == "0":
	obsList.SetLineColor(1)
	obsList.SetLineWidth(3)
	obsList.SetLineStyle(1)
	obsList.SetMarkerStyle(20)
	obsList.SetMarkerSize(1.5)
	obsList.Draw('LPsame')
	print("Observed Limit: ")
	print(obsList)
		

        marcObsList.SetLineColor(kMagenta)
        marcObsList.SetLineWidth(3)
        marcObsList.SetLineStyle(1)
        marcObsList.SetMarkerStyle(20)
        marcObsList.SetMarkerSize(1.5)
	if marc: marcObsList.Draw('Lsame')
        
        marcExpectList.SetLineColor(kMagenta)
        marcExpectList.SetLineWidth(3)
        marcExpectList.SetLineStyle(2)
        marcExpectList.SetMarkerStyle(20)
        marcExpectList.SetMarkerSize(1.5)
	if marc: marcExpectList.Draw('Lsame')

#	obsList.SetMarkerStyle(20)
#	obsList.SetLineWidth(3)
#	if options.observed:
#		#obsList.Draw('PLsame')
#		obsList.Draw('Csame')


	x_range = [5., 155.]
	vtype = "vector"
	if constraint:
		# Z constraint
		tf_Z_constraint = TF1("Z_constraint", gq_Z_constraint, x_range[0], x_range[1], 0)
		tf_Z_constraint.SetNpx(1000)
		tf_Z_constraint.SetLineColor(17)
		ROOT.gStyle.SetLineStyleString(9, "40 20");
		tf_Z_constraint.SetLineStyle(9)
		tf_Z_constraint.SetLineWidth(2)
		tf_Z_constraint.Draw("same")
		z_width_legend_entry="Z width (all #Gamma_{Z'}/M_{Z'})"
	
		# Upsilon constraint
		tf_upsilon_constraint = ROOT.TF1("upsilon_constraint", gq_upsilon_constraint, x_range[0], x_range[1], 0)
		tf_upsilon_constraint.SetNpx(1000)
		tf_upsilon_constraint.SetLineColor(17)
		ROOT.gStyle.SetLineStyleString(11, "20 10");
		tf_upsilon_constraint.SetLineStyle(11)
		tf_upsilon_constraint.SetLineWidth(2)
		tf_upsilon_constraint.Draw("same")
		upsilon_width_legend_entry="#Upsilon width (all #Gamma_{Z'}/M_{Z'})"

	gom_x = None
	logx = False
	logy = True
	if GOM:
		# Lines at fixed Gamma / M
		GoM_tf1s = {}
		GoM_pyfuncs = {}
		GoM_labels = {}

		GoMs = [0.05, 0.1, 0.3, 0.5, 1.0]
		for i, GoM in enumerate(GoMs):
			# def gom_to_gq(gom, mZp, vtype):
			GoM_pyfuncs[GoM] = lambda x, p: gom_to_gq(p[0], x[0], vtype)
			#print("[DEBUG] {GoM_pyfuncs[GoM]([100.], [])}")
			GoM_tf1s[GoM] = ROOT.TF1("tf1_gq_{}".format(GoM), GoM_pyfuncs[GoM], x_range[0], x_range[1], 1) # 
			GoM_tf1s[GoM].SetParameter(0, GoM)
			#ROOT.SetOwnership(GoM_tf1s[GoM], False)
			GoM_tf1s[GoM].SetLineColor(ROOT.kGray+1)
			GoM_tf1s[GoM].SetLineStyle(ROOT.kDashed)
			GoM_tf1s[GoM].SetLineWidth(1)
			GoM_tf1s[GoM].Draw("same")

			# TLatex for Gamma / M
			if gom_x:
				label_x = gom_x
			else: 
				if logx:
					label_xfrac = 0.05
				else:
					label_xfrac = 0.864
				label_x = (x_range[1] - x_range[0]) * label_xfrac + x_range[0]
			if logy:
				label_y = GoM_tf1s[GoM].Eval(label_x) * 0.85
				gom_text = "#Gamma_{{Z'}}#kern[-0.6]{{ }}/#kern[-0.7]{{ }}M_{{Z'}}#kern[-0.7]{{ }}=#kern[-0.7]{{ }}{}%".format(int(GoM * 100))
			else:
				# label_y = GoM_tf1s[GoM].Eval(label_x) - 0.085 # For labels under the line
				label_y = GoM_tf1s[GoM].Eval(label_x) + 0.05 # For labels over the line
				gom_text = "#frac{{#Gamma}}{{M_{{Z'}}}} = {}%".format(int(GoM * 100))
			GoM_labels[GoM] = ROOT.TLatex(label_x, label_y, gom_text)
			if logy:
				GoM_labels[GoM].SetTextSize(0.028)
			else:
				GoM_labels[GoM].SetTextSize(0.027)
			GoM_labels[GoM].SetTextColor(ROOT.kGray+1)
			GoM_labels[GoM].Draw("same")
	
	
	AddCMSLumi(ROOT.gPad, LUMI, "Preliminary")
	
	legend = ROOT.TLegend(0.55, 0.6, 0.8, 0.89)
	legend.SetFillStyle(0)
	legend.SetBorderSize(0)
	legend.SetTextSize(0.038)
	legend.SetTextFont(42)
#	if options.observed:
#		legend.AddEntry(obs, "Observed",'l')
	#legend.AddEntry(medianList, "Asymptotic CL_{s} expected",'L')
	legend.AddEntry(obsList, "Observed Limits", 'PL')
	legend.AddEntry(medianList, "Expected",'L')
	legend.AddEntry(greenList, "68% Expected",'F')
	legend.AddEntry(yellowList,"95% Expected",'F')
        if marc: legend.AddEntry(marcObsList, "2016 Analysis Observed Limits", 'CL')
        if marc: legend.AddEntry(marcExpectList, "2016 Analysis Expected Limits", 'CL')
	if constraint:
		legend.AddEntry(tf_Z_constraint, z_width_legend_entry, "l")
		legend.AddEntry(tf_upsilon_constraint, upsilon_width_legend_entry, "l")

	legend.Draw()
	
	legend.Draw("same")
	
	
	if marc: c.SaveAs("/home/akobert/CMSSW_11_3_4/src/FitTest/Unblind/limits/limits_FULL_marc_pretty.png")
	else: c.SaveAs("/home/akobert/CMSSW_11_3_4/src/FitTest/Unblind/limits/limits_FULL_pretty.png")
	
	c.SetLogx() # Setting Log X to match Marc
	h_limit.GetXaxis().SetLimits(6,500)
	if marc: c.SaveAs("/home/akobert/CMSSW_11_3_4/src/FitTest/Unblind/limits/limits_FULL_marc_log_pretty.png")
	else: c.SaveAs("/home/akobert/CMSSW_11_3_4/src/FitTest/Unblind/limits/limits_FULL_log_pretty.png")
	
	
