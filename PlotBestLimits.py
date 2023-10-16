#Limit Plotting Code From Marc

import ROOT
from ROOT import *
from array import array
​
ROOT.gStyle.SetPadRightMargin(0.08);
ROOT.gStyle.SetPadLeftMargin(0.11);
ROOT.gStyle.SetPadTopMargin(0.10);
ROOT.gStyle.SetPalette(1);
​
def AddCMSLumi(pad, fb, extra):
	cmsText     = "#bf{CMS} " + extra
	cmsTextFont   = 42  
	lumiTextSize     = 0.43
	lumiTextOffset   = 0.15
	cmsTextSize      = 0.47
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
	latex.DrawLatex(0.975,0.95,lumiText)
	pad.cd()
	latex.SetTextFont(cmsTextFont)
	latex.SetTextSize(cmsTextSize*t)
	latex.SetTextAlign(11)
	latex.DrawLatex(0.11,0.95, cmsText)
	pad.Update()
​
def SetGfN(G, N, S):
	G.SetPoint(G.GetN(),0.750,45.9397*N*N)
	G.SetPoint(G.GetN(),0.875,20.1308*N*N)
	G.SetPoint(G.GetN(),1.000,9.59447*N*N)
	G.SetPoint(G.GetN(),1.125,4.88278*N*N)
	G.SetPoint(G.GetN(),1.250,2.61745*N*N)
	G.SetPoint(G.GetN(),1.375,1.46371*N*N)
	G.SetPoint(G.GetN(),1.500,0.847454*N*N)
	G.SetPoint(G.GetN(),1.625,0.505322*N*N)
	G.SetPoint(G.GetN(),1.750,0.309008*N*N)
	G.SetPoint(G.GetN(),1.875,0.192939*N*N)
	G.SetPoint(G.GetN(),2.000,0.122826*N*N)
	G.SetPoint(G.GetN(),2.125,0.0795248*N*N)
	G.SetPoint(G.GetN(),2.250,0.0522742*N*N)
	G.SetPoint(G.GetN(),2.375,0.0348093*N*N)
	G.SetPoint(G.GetN(),2.500,0.0235639*N*N)
	G.SetPoint(G.GetN(),2.625,0.0161926*N*N)
	G.SetPoint(G.GetN(),2.750,0.0109283*N*N)
	G.SetPoint(G.GetN(),2.875,0.00759881*N*N)
	G.SetPoint(G.GetN(),3.000,0.00531361*N*N)
	G.SetLineStyle(S)
	G.SetLineColor(kBlue)
	G.SetLineWidth(1)
G1 = TGraph()
G2 = TGraph()
G8 = TGraph()
G4 = TGraph()
SetGfN(G1, 1, 1)
SetGfN(G2, 2, 2)
SetGfN(G4, 4, 3)
SetGfN(G8, 8, 4)
​
Keep = []
​
def GetLimit(A, X, which): # 5=obs, 2=exp
	F = TFile("FINAL/higgsCombineX"+str(X)+"A"+str(A)+".AsymptoticLimits.mH120.root")
	T = F.Get("limit")
	T.GetEntry(which)
	return T.limit*20.
def FormatObs(A, C, M):
	A.SetMarkerColor(C)
	A.SetMarkerStyle(M)
	A.SetMarkerSize(0.7)
def FormatExp(A, C, M):
	A.SetLineColor(C)
	A.SetLineStyle(2)
def makeAFillGraph(listx,listy1,listy2,linecolor, fillcolor, fillstyle):
	a_m = array('f', []);
	a_g = array('f', []);
	for i in range(len(listx)):
		a_m.append(listx[i]);
		a_g.append(listy1[i]);	
	for i in range(len(listx)-1,-1,-1):
		a_m.append(listx[i]);
		a_g.append(listy2[i]);
	gr = ROOT.TGraph(2*len(listx),a_m,a_g);
	gr.SetLineColor(linecolor)
	gr.SetFillColor(fillcolor)
	gr.SetFillStyle(fillstyle)
	return gr
​
Fake1 = ROOT.TGraph()
Fake1.SetLineColor(1)
Fake1.SetLineStyle(2)
Fake1.SetLineWidth(2)
Fake1.SetFillColor(ROOT.kGreen+1)
Fake2 = ROOT.TGraph()
Fake2.SetLineColor(1)
Fake2.SetLineStyle(2)
Fake2.SetLineWidth(2)
Fake2.SetFillColor(ROOT.kOrange)
​
​
C = TCanvas("", "", 800,800)
C.cd()
​
Leg = []
​
def MakeAPad(C, M,x1,y1,x2,y2):
    C.cd()
    P = TPad(str(M), "", x1,y1,x2,y2)
    P.SetFillStyle(0)
    P.Draw()
    P.cd()
    a = str(M)
    OBS = TGraph()
    EXP = TGraph()
    pMin = 10000000
    pMax = 0
    X = []
    M1 = []
    M2 = []
    P1 = []
    P2 = []
    for x in ['1000', '1250', '1500', '1750', '2000', '2250', '2500', '2750', '3000']:
        X.append(float(x)/1000.)
        O = GetLimit(a, x, 5)#*float(x)/1000.
        E = GetLimit(a, x, 2)#*float(x)/1000.
        m2 = GetLimit(a, x, 0)#*float(x)/1000.
        m1 = GetLimit(a, x, 1)#*float(x)/1000.
        p1 = GetLimit(a, x, 3)#*float(x)/1000.
        p2 = GetLimit(a, x, 4)#*float(x)/1000.
        pMin = min(O,pMin, m2)
        pMax = max(O,pMax, p2)
        OBS.SetPoint(OBS.GetN(), float(x)/1000., O)
        EXP.SetPoint(EXP.GetN(), float(x)/1000., E)
        M2.append(m2)
        M1.append(m1)
        P2.append(p2)
        P1.append(p1)
    FormatObs(OBS, 1, 20)
    FormatExp(EXP, 1, 20)
    Onesig = makeAFillGraph(X,M1,P1,kGreen+1, kGreen+1, 1001)
    Twosig = makeAFillGraph(X,M2,P2,kOrange, kOrange, 1001)
    Keep.append(OBS)
    Keep.append(EXP)
    Keep.append(Onesig)
    Keep.append(Twosig)
​
    T = TH2F("T"+a, "", 5, 1, 3, 20, pMin*0.8, pMax*2.25)
    T.GetYaxis().SetLabelSize(0.11)
    T.GetXaxis().SetLabelSize(0.11)
    T.GetXaxis().SetLabelOffset(0.012)
    T.SetStats(0)
    T.Draw()
    Keep.append(T)
    P.SetLogy()
    Twosig.Draw("Fsames")
    Onesig.Draw("Fsames")
    OBS.Draw("PL")
    EXP.Draw("L")
​
    iP = TPaveText(0.5,0.55,0.89,0.89,"NDC")
    iP.SetFillStyle(0)
    iP.AddText("#it{m_{#phi}} = " + a + " GeV")
    iP.SetTextSize(0.165)
    iP.SetFillColor(0)
    iP.SetBorderSize(0)
    iP.SetTextFont(42)
    Keep.append(iP)
    iP.Draw()
​
    gPad.SetTicks(1,1)
    G1.Draw("L")
    G2.Draw("L")
    G8.Draw("L")
    G4.Draw("L")
    P.Update()
    #C.Update()
​
G = [
    [0,0],
    [1,0],
    [2,0],
    [3,0],
    [4,0],
    [5,0],
    [0,1],
    [1,1],
    [2,1],
    [3,1],
    [4,1],
    [5,1],
    [0,2],
    [1,2],
    [2,2],
    [3,2],
    [4,2],
    [5,2],
]
s = 25
for i in range(16):
    c = G[i][1]
    r = G[i][0]    
    print(c)
    print(r)
    x1 = 0.07 + c*(1./3. * 0.93)
    y1 = 0.94 - (r+1)*(1./6. * 0.88)
    x2 = min(1.0, 0.07 + (c+1)*(1./3. * 0.93)+0.025)
    y2 = 0.94 - (r)*(1./6. * 0.88)
    p = MakeAPad(C, s, x1, y1, x2, y2)
    s+=5
​
C.cd()
L = TLegend(0.71,0.075,0.875,0.275)
L.SetLineColor(0)
L.SetFillColor(0)
L.SetHeader("#bf{95% CL limits:}")
L.AddEntry(Keep[0], "observed", "PL")
L.AddEntry(Fake1, "expected #pm 1#sigma", "FL")
L.AddEntry(Fake2, "expected #pm 2#sigma", "FL")
L.Draw()
L2 = TLegend(0.875,0.06,0.99,0.305)
L2.SetLineColor(0)
L2.SetFillColor(0)
L2.AddEntry(G1, "m_{X}N/f = 1", "L")
L2.AddEntry(G2, "m_{X}N/f = 2", "L")
L2.AddEntry(G4, "m_{X}N/f = 4", "L")
L2.AddEntry(G8, "m_{X}N/f = 8", "L")
L2.Draw()
C.Update()
AddCMSLumi(C, 138, "")
​
x1P = TPaveText(0.09,   0.02,   0.07+1./3.*0.93+0.025,   0.055,"NDC")
x1P.AddText("#it{m}_{X} (TeV)")
x1P.SetTextSize(0.031)
x1P.SetFillColor(0)
x1P.SetBorderSize(0)
x1P.SetTextFont(42)
x1P.Draw()
​
x2P = TPaveText(0.07+1./3.*0.93+0.02,   0.02,   0.07+2./3.*0.93+0.025,   0.055,"NDC")
x2P.AddText("#it{m}_{X} (TeV)")
x2P.SetTextSize(0.031)
x2P.SetFillColor(0)
x2P.SetBorderSize(0)
x2P.SetTextFont(42)
x2P.Draw()
​
x3P = TPaveText(0.07+2./3.*0.93+0.02,   0.02+1./3.*0.88,   1.0,   0.055+1./3.*0.88,"NDC")
x3P.AddText("#it{m}_{X} (TeV)")
x3P.SetTextSize(0.031)
x3P.SetFillColor(0)
x3P.SetBorderSize(0)
x3P.SetTextFont(42)
x3P.Draw()
​
P = TPaveText(0.015,0.,0.06,1.0,"NDC")
T = P.AddText("#sigma[pp #rightarrow X #rightarrow #it{#phi#phi} #rightarrow (b#bar{b})(b#bar{b})] (fb)")
P.SetTextSize(0.0523)
P.SetFillColor(0)
T.SetTextAngle(90)
P.SetBorderSize(0)
P.SetTextFont(42)
P.Draw()
​
C.Update()
C.Print("LimitsFinal.pdf")
C.Print("LimitsFinal.png")
C.Print("LimitsFinal.root")
