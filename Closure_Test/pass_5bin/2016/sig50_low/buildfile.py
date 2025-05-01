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

#from future import division

class Build:
	def __init__(self, name, bfile1, bfile2, bfile3, bfile4, bfile5, bfile6, bfile7, bfile8, dfile1, dfile2, sfile1, sfile2):
		gROOT.SetBatch(True)

		#Output Files
	        ofile = ROOT.TFile(name + ".root", "RECREATE")
	        ofile.cd()


		#background files
		self.f = TFile.Open(dfile1, "READ")
		self.f.ls();

		self.Data_p7_w = self.f.Get("jet_pt_soft_pass_wide10_wide")
		self.Data_f7_w = self.f.Get("jet_pt_soft_fail_wide10_wide")
		self.Data_p1 = self.f.Get("pass_soft")
		self.Data_f1 = self.f.Get("fail_soft")
		

		self.g = TFile(bfile1, "READ")
		self.TTBar_p7_w = self.g.Get("jet_pt_soft_pass_wide10_wide")
		self.TTBar_f7_w = self.g.Get("jet_pt_soft_fail_wide10_wide")
		
		self.TTBar_p1 = self.g.Get("pass_soft")	
		self.TTBar_f1 = self.g.Get("fail_soft")	
		
		self.TTBar_p7_w_puUp = self.g.Get("jet_pt_soft_pass_wide10_wide_puUp")
		self.TTBar_f7_w_puUp = self.g.Get("jet_pt_soft_fail_wide10_wide_puUp")
		self.TTBar_p1_puUp = self.g.Get("pass_soft_puUp")	
		self.TTBar_f1_puUp = self.g.Get("fail_soft_puUp")	
		
		self.TTBar_p7_w_puDown = self.g.Get("jet_pt_soft_pass_wide10_wide_puDown")
		self.TTBar_f7_w_puDown = self.g.Get("jet_pt_soft_fail_wide10_wide_puDown")
		self.TTBar_p1_puDown = self.g.Get("pass_soft_puDown")	
		self.TTBar_f1_puDown = self.g.Get("fail_soft_puDown")	
		
		self.TTBar_p7_w_jerUp = self.g.Get("jet_pt_soft_pass_wide10_wide_jerUp")
		self.TTBar_f7_w_jerUp = self.g.Get("jet_pt_soft_fail_wide10_wide_jerUp")
		self.TTBar_p1_jerUp = self.g.Get("pass_soft_jerUp")	
		self.TTBar_f1_jerUp = self.g.Get("fail_soft_jerUp")	
		
		self.TTBar_p7_w_jerDown = self.g.Get("jet_pt_soft_pass_wide10_wide_jerDown")
		self.TTBar_f7_w_jerDown = self.g.Get("jet_pt_soft_fail_wide10_wide_jerDown")
		self.TTBar_p1_jerDown = self.g.Get("pass_soft_jerDown")	
		self.TTBar_f1_jerDown = self.g.Get("fail_soft_jerDown")	
		
		self.TTBar_p7_w_jesUp = self.g.Get("jet_pt_soft_pass_wide10_wide_jesUp")
		self.TTBar_f7_w_jesUp = self.g.Get("jet_pt_soft_fail_wide10_wide_jesUp")
		self.TTBar_p1_jesUp = self.g.Get("pass_soft_jesUp")	
		self.TTBar_f1_jesUp = self.g.Get("fail_soft_jesUp")	
		
		self.TTBar_p7_w_jesDown = self.g.Get("jet_pt_soft_pass_wide10_wide_jesDown")
		self.TTBar_f7_w_jesDown = self.g.Get("jet_pt_soft_fail_wide10_wide_jesDown")
		self.TTBar_p1_jesDown = self.g.Get("pass_soft_jesDown")	
		self.TTBar_f1_jesDown = self.g.Get("fail_soft_jesDown")	
		
		self.h = TFile(bfile2, "READ")
		self.WGamma_p7_w = self.h.Get("jet_pt_soft_pass_wide10_wide")
		self.WGamma_f7_w = self.h.Get("jet_pt_soft_fail_wide10_wide")
		
		self.WGamma_p1 = self.h.Get("pass_soft")	
		self.WGamma_f1 = self.h.Get("fail_soft")	
		
		self.WGamma_p7_w_puUp = self.h.Get("jet_pt_soft_pass_wide10_wide_puUp")
		self.WGamma_f7_w_puUp = self.h.Get("jet_pt_soft_fail_wide10_wide_puUp")
		self.WGamma_p1_puUp = self.h.Get("pass_soft_puUp")	
		self.WGamma_f1_puUp = self.h.Get("fail_soft_puUp")	
		
		self.WGamma_p7_w_puDown = self.h.Get("jet_pt_soft_pass_wide10_wide_puDown")
		self.WGamma_f7_w_puDown = self.h.Get("jet_pt_soft_fail_wide10_wide_puDown")
		self.WGamma_p1_puDown = self.h.Get("pass_soft_puDown")	
		self.WGamma_f1_puDown = self.h.Get("fail_soft_puDown")	
		
		self.WGamma_p7_w_jerUp = self.h.Get("jet_pt_soft_pass_wide10_wide_jerUp")
		self.WGamma_f7_w_jerUp = self.h.Get("jet_pt_soft_fail_wide10_wide_jerUp")
		self.WGamma_p1_jerUp = self.h.Get("pass_soft_jerUp")	
		self.WGamma_f1_jerUp = self.h.Get("fail_soft_jerUp")	
		
		self.WGamma_p7_w_jerDown = self.h.Get("jet_pt_soft_pass_wide10_wide_jerDown")
		self.WGamma_f7_w_jerDown = self.h.Get("jet_pt_soft_fail_wide10_wide_jerDown")
		self.WGamma_p1_jerDown = self.h.Get("pass_soft_jerDown")	
		self.WGamma_f1_jerDown = self.h.Get("fail_soft_jerDown")	
		
		self.WGamma_p7_w_jesUp = self.h.Get("jet_pt_soft_pass_wide10_wide_jesUp")
		self.WGamma_f7_w_jesUp = self.h.Get("jet_pt_soft_fail_wide10_wide_jesUp")
		self.WGamma_p1_jesUp = self.h.Get("pass_soft_jesUp")	
		self.WGamma_f1_jesUp = self.h.Get("fail_soft_jesUp")	
		
		self.WGamma_p7_w_jesDown = self.h.Get("jet_pt_soft_pass_wide10_wide_jesDown")
		self.WGamma_f7_w_jesDown = self.h.Get("jet_pt_soft_fail_wide10_wide_jesDown")
		self.WGamma_p1_jesDown = self.h.Get("pass_soft_jesDown")	
		self.WGamma_f1_jesDown = self.h.Get("fail_soft_jesDown")	
		
		self.j = TFile(bfile3, "READ")
		self.ZGamma_p7_w = self.j.Get("jet_pt_soft_pass_wide10_wide")
		self.ZGamma_f7_w = self.j.Get("jet_pt_soft_fail_wide10_wide")
		
		self.ZGamma_p1 = self.j.Get("pass_soft")	
		self.ZGamma_f1 = self.j.Get("fail_soft")	
		
		self.ZGamma_p7_w_puUp = self.j.Get("jet_pt_soft_pass_wide10_wide_puUp")
		self.ZGamma_f7_w_puUp = self.j.Get("jet_pt_soft_fail_wide10_wide_puUp")
		self.ZGamma_p1_puUp = self.j.Get("pass_soft_puUp")	
		self.ZGamma_f1_puUp = self.j.Get("fail_soft_puUp")	
		
		self.ZGamma_p7_w_puDown = self.j.Get("jet_pt_soft_pass_wide10_wide_puDown")
		self.ZGamma_f7_w_puDown = self.j.Get("jet_pt_soft_fail_wide10_wide_puDown")
		self.ZGamma_p1_puDown = self.j.Get("pass_soft_puDown")	
		self.ZGamma_f1_puDown = self.j.Get("fail_soft_puDown")	
		
		self.ZGamma_p7_w_jerUp = self.j.Get("jet_pt_soft_pass_wide10_wide_jerUp")
		self.ZGamma_f7_w_jerUp = self.j.Get("jet_pt_soft_fail_wide10_wide_jerUp")
		self.ZGamma_p1_jerUp = self.j.Get("pass_soft_jerUp")	
		self.ZGamma_f1_jerUp = self.j.Get("fail_soft_jerUp")	
		
		self.ZGamma_p7_w_jerDown = self.j.Get("jet_pt_soft_pass_wide10_wide_jerDown")
		self.ZGamma_f7_w_jerDown = self.j.Get("jet_pt_soft_fail_wide10_wide_jerDown")
		self.ZGamma_p1_jerDown = self.j.Get("pass_soft_jerDown")	
		self.ZGamma_f1_jerDown = self.j.Get("fail_soft_jerDown")	
		
		self.ZGamma_p7_w_jesUp = self.j.Get("jet_pt_soft_pass_wide10_wide_jesUp")
		self.ZGamma_f7_w_jesUp = self.j.Get("jet_pt_soft_fail_wide10_wide_jesUp")
		self.ZGamma_p1_jesUp = self.j.Get("pass_soft_jesUp")	
		self.ZGamma_f1_jesUp = self.j.Get("fail_soft_jesUp")	
		
		self.ZGamma_p7_w_jesDown = self.j.Get("jet_pt_soft_pass_wide10_wide_jesDown")
		self.ZGamma_f7_w_jesDown = self.j.Get("jet_pt_soft_fail_wide10_wide_jesDown")
		self.ZGamma_p1_jesDown = self.j.Get("pass_soft_jesDown")	
		self.ZGamma_f1_jesDown = self.j.Get("fail_soft_jesDown")	
		
		self.k = TFile(bfile4, "READ")
		self.GJ_p7_w = self.k.Get("jet_pt_soft_pass_wide10_wide")
		self.GJ_f7_w = self.k.Get("jet_pt_soft_fail_wide10_wide")
		
		self.GJ_p1 = self.k.Get("pass_soft")	
		self.GJ_f1 = self.k.Get("fail_soft")	
		
		self.GJ_p7_w_puUp = self.k.Get("jet_pt_soft_pass_wide10_wide_puUp")
		self.GJ_f7_w_puUp = self.k.Get("jet_pt_soft_fail_wide10_wide_puUp")
		self.GJ_p1_puUp = self.k.Get("pass_soft_puUp")	
		self.GJ_f1_puUp = self.k.Get("fail_soft_puUp")	
		
		self.GJ_p7_w_puDown = self.k.Get("jet_pt_soft_pass_wide10_wide_puDown")
		self.GJ_f7_w_puDown = self.k.Get("jet_pt_soft_fail_wide10_wide_puDown")
		self.GJ_p1_puDown = self.k.Get("pass_soft_puDown")	
		self.GJ_f1_puDown = self.k.Get("fail_soft_puDown")	
		
		self.GJ_p7_w_jerUp = self.k.Get("jet_pt_soft_pass_wide10_wide_jerUp")
		self.GJ_f7_w_jerUp = self.k.Get("jet_pt_soft_fail_wide10_wide_jerUp")
		self.GJ_p1_jerUp = self.k.Get("pass_soft_jerUp")	
		self.GJ_f1_jerUp = self.k.Get("fail_soft_jerUp")	
		
		self.GJ_p7_w_jerDown = self.k.Get("jet_pt_soft_pass_wide10_wide_jerDown")
		self.GJ_f7_w_jerDown = self.k.Get("jet_pt_soft_fail_wide10_wide_jerDown")
		self.GJ_p1_jerDown = self.k.Get("pass_soft_jerDown")	
		self.GJ_f1_jerDown = self.k.Get("fail_soft_jerDown")	
		
		self.GJ_p7_w_jesUp = self.k.Get("jet_pt_soft_pass_wide10_wide_jesUp")
		self.GJ_f7_w_jesUp = self.k.Get("jet_pt_soft_fail_wide10_wide_jesUp")
		self.GJ_p1_jesUp = self.k.Get("pass_soft_jesUp")	
		self.GJ_f1_jesUp = self.k.Get("fail_soft_jesUp")	
		
		self.GJ_p7_w_jesDown = self.k.Get("jet_pt_soft_pass_wide10_wide_jesDown")
		self.GJ_f7_w_jesDown = self.k.Get("jet_pt_soft_fail_wide10_wide_jesDown")
		self.GJ_p1_jesDown = self.k.Get("pass_soft_jesDown")	
		self.GJ_f1_jesDown = self.k.Get("fail_soft_jesDown")	

		self.s1 = TFile(sfile1, "READ")
		self.Sig50_p7_w = self.s1.Get("jet_pt_soft_pass_wide10_wide")
		self.Sig50_f7_w = self.s1.Get("jet_pt_soft_fail_wide10_wide")
		
		self.Sig50_p1 = self.s1.Get("pass_soft")	
		self.Sig50_f1 = self.s1.Get("fail_soft")	
		
		self.Sig50_p7_w_puUp = self.s1.Get("jet_pt_soft_pass_wide10_wide_puUp")
		self.Sig50_f7_w_puUp = self.s1.Get("jet_pt_soft_fail_wide10_wide_puUp")
		self.Sig50_p1_puUp = self.s1.Get("pass_soft_puUp")	
		self.Sig50_f1_puUp = self.s1.Get("fail_soft_puUp")	
		
		self.Sig50_p7_w_puDown = self.s1.Get("jet_pt_soft_pass_wide10_wide_puDown")
		self.Sig50_f7_w_puDown = self.s1.Get("jet_pt_soft_fail_wide10_wide_puDown")
		self.Sig50_p1_puDown = self.s1.Get("pass_soft_puDown")	
		self.Sig50_f1_puDown = self.s1.Get("fail_soft_puDown")	
		
		self.Sig50_p7_w_jerUp = self.s1.Get("jet_pt_soft_pass_wide10_wide_jerUp")
		self.Sig50_f7_w_jerUp = self.s1.Get("jet_pt_soft_fail_wide10_wide_jerUp")
		self.Sig50_p1_jerUp = self.s1.Get("pass_soft_jerUp")	
		self.Sig50_f1_jerUp = self.s1.Get("fail_soft_jerUp")	
		
		self.Sig50_p7_w_jerDown = self.s1.Get("jet_pt_soft_pass_wide10_wide_jerDown")
		self.Sig50_f7_w_jerDown = self.s1.Get("jet_pt_soft_fail_wide10_wide_jerDown")
		self.Sig50_p1_jerDown = self.s1.Get("pass_soft_jerDown")	
		self.Sig50_f1_jerDown = self.s1.Get("fail_soft_jerDown")	
		
		self.Sig50_p7_w_jesUp = self.s1.Get("jet_pt_soft_pass_wide10_wide_jesUp")
		self.Sig50_f7_w_jesUp = self.s1.Get("jet_pt_soft_fail_wide10_wide_jesUp")
		self.Sig50_p1_jesUp = self.s1.Get("pass_soft_jesUp")	
		self.Sig50_f1_jesUp = self.s1.Get("fail_soft_jesUp")	
		
		self.Sig50_p7_w_jesDown = self.s1.Get("jet_pt_soft_pass_wide10_wide_jesDown")
		self.Sig50_f7_w_jesDown = self.s1.Get("jet_pt_soft_fail_wide10_wide_jesDown")
		self.Sig50_p1_jesDown = self.s1.Get("pass_soft_jesDown")	
		self.Sig50_f1_jesDown = self.s1.Get("fail_soft_jesDown")	


		self.l = TFile.Open(dfile2, "READ")
		self.l.ls();

		self.Data_APV_p7_w = self.l.Get("jet_pt_soft_pass_wide10_wide")
		self.Data_APV_f7_w = self.l.Get("jet_pt_soft_fail_wide10_wide")
		self.Data_APV_p1 = self.l.Get("pass_soft")
		self.Data_APV_f1 = self.l.Get("fail_soft")
		

		self.m = TFile(bfile5, "READ")
		self.TTBar_APV_p7_w = self.m.Get("jet_pt_soft_pass_wide10_wide")
		self.TTBar_APV_f7_w = self.m.Get("jet_pt_soft_fail_wide10_wide")
		
		self.TTBar_APV_p1 = self.m.Get("pass_soft")	
		self.TTBar_APV_f1 = self.m.Get("fail_soft")	
		
		self.TTBar_APV_p7_w_puUp = self.m.Get("jet_pt_soft_pass_wide10_wide_puUp")
		self.TTBar_APV_f7_w_puUp = self.m.Get("jet_pt_soft_fail_wide10_wide_puUp")
		self.TTBar_APV_p1_puUp = self.m.Get("pass_soft_puUp")	
		self.TTBar_APV_f1_puUp = self.m.Get("fail_soft_puUp")	
		
		self.TTBar_APV_p7_w_puDown = self.m.Get("jet_pt_soft_pass_wide10_wide_puDown")
		self.TTBar_APV_f7_w_puDown = self.m.Get("jet_pt_soft_fail_wide10_wide_puDown")
		self.TTBar_APV_p1_puDown = self.m.Get("pass_soft_puDown")	
		self.TTBar_APV_f1_puDown = self.m.Get("fail_soft_puDown")	
		
		self.TTBar_APV_p7_w_jerUp = self.m.Get("jet_pt_soft_pass_wide10_wide_jerUp")
		self.TTBar_APV_f7_w_jerUp = self.m.Get("jet_pt_soft_fail_wide10_wide_jerUp")
		self.TTBar_APV_p1_jerUp = self.m.Get("pass_soft_jerUp")	
		self.TTBar_APV_f1_jerUp = self.m.Get("fail_soft_jerUp")	
		
		self.TTBar_APV_p7_w_jerDown = self.m.Get("jet_pt_soft_pass_wide10_wide_jerDown")
		self.TTBar_APV_f7_w_jerDown = self.m.Get("jet_pt_soft_fail_wide10_wide_jerDown")
		self.TTBar_APV_p1_jerDown = self.m.Get("pass_soft_jerDown")	
		self.TTBar_APV_f1_jerDown = self.m.Get("fail_soft_jerDown")	
		
		self.TTBar_APV_p7_w_jesUp = self.m.Get("jet_pt_soft_pass_wide10_wide_jesUp")
		self.TTBar_APV_f7_w_jesUp = self.m.Get("jet_pt_soft_fail_wide10_wide_jesUp")
		self.TTBar_APV_p1_jesUp = self.m.Get("pass_soft_jesUp")	
		self.TTBar_APV_f1_jesUp = self.m.Get("fail_soft_jesUp")	
		
		self.TTBar_APV_p7_w_jesDown = self.m.Get("jet_pt_soft_pass_wide10_wide_jesDown")
		self.TTBar_APV_f7_w_jesDown = self.m.Get("jet_pt_soft_fail_wide10_wide_jesDown")
		self.TTBar_APV_p1_jesDown = self.m.Get("pass_soft_jesDown")	
		self.TTBar_APV_f1_jesDown = self.m.Get("fail_soft_jesDown")	
		
		self.n = TFile(bfile6, "READ")
		self.WGamma_APV_p7_w = self.n.Get("jet_pt_soft_pass_wide10_wide")
		self.WGamma_APV_f7_w = self.n.Get("jet_pt_soft_fail_wide10_wide")
		
		self.WGamma_APV_p1 = self.n.Get("pass_soft")	
		self.WGamma_APV_f1 = self.n.Get("fail_soft")	
		
		self.WGamma_APV_p7_w_puUp = self.n.Get("jet_pt_soft_pass_wide10_wide_puUp")
		self.WGamma_APV_f7_w_puUp = self.n.Get("jet_pt_soft_fail_wide10_wide_puUp")
		self.WGamma_APV_p1_puUp = self.n.Get("pass_soft_puUp")	
		self.WGamma_APV_f1_puUp = self.n.Get("fail_soft_puUp")	
		
		self.WGamma_APV_p7_w_puDown = self.n.Get("jet_pt_soft_pass_wide10_wide_puDown")
		self.WGamma_APV_f7_w_puDown = self.n.Get("jet_pt_soft_fail_wide10_wide_puDown")
		self.WGamma_APV_p1_puDown = self.n.Get("pass_soft_puDown")	
		self.WGamma_APV_f1_puDown = self.n.Get("fail_soft_puDown")	
		
		self.WGamma_APV_p7_w_jerUp = self.n.Get("jet_pt_soft_pass_wide10_wide_jerUp")
		self.WGamma_APV_f7_w_jerUp = self.n.Get("jet_pt_soft_fail_wide10_wide_jerUp")
		self.WGamma_APV_p1_jerUp = self.n.Get("pass_soft_jerUp")	
		self.WGamma_APV_f1_jerUp = self.n.Get("fail_soft_jerUp")	
		
		self.WGamma_APV_p7_w_jerDown = self.n.Get("jet_pt_soft_pass_wide10_wide_jerDown")
		self.WGamma_APV_f7_w_jerDown = self.n.Get("jet_pt_soft_fail_wide10_wide_jerDown")
		self.WGamma_APV_p1_jerDown = self.n.Get("pass_soft_jerDown")	
		self.WGamma_APV_f1_jerDown = self.n.Get("fail_soft_jerDown")	
		
		self.WGamma_APV_p7_w_jesUp = self.n.Get("jet_pt_soft_pass_wide10_wide_jesUp")
		self.WGamma_APV_f7_w_jesUp = self.n.Get("jet_pt_soft_fail_wide10_wide_jesUp")
		self.WGamma_APV_p1_jesUp = self.n.Get("pass_soft_jesUp")	
		self.WGamma_APV_f1_jesUp = self.n.Get("fail_soft_jesUp")	
		
		self.WGamma_APV_p7_w_jesDown = self.n.Get("jet_pt_soft_pass_wide10_wide_jesDown")
		self.WGamma_APV_f7_w_jesDown = self.n.Get("jet_pt_soft_fail_wide10_wide_jesDown")
		self.WGamma_APV_p1_jesDown = self.n.Get("pass_soft_jesDown")	
		self.WGamma_APV_f1_jesDown = self.n.Get("fail_soft_jesDown")	
		
		self.o = TFile(bfile7, "READ")
		self.ZGamma_APV_p7_w = self.o.Get("jet_pt_soft_pass_wide10_wide")
		self.ZGamma_APV_f7_w = self.o.Get("jet_pt_soft_fail_wide10_wide")
		
		self.ZGamma_APV_p1 = self.o.Get("pass_soft")	
		self.ZGamma_APV_f1 = self.o.Get("fail_soft")	
		
		self.ZGamma_APV_p7_w_puUp = self.o.Get("jet_pt_soft_pass_wide10_wide_puUp")
		self.ZGamma_APV_f7_w_puUp = self.o.Get("jet_pt_soft_fail_wide10_wide_puUp")
		self.ZGamma_APV_p1_puUp = self.o.Get("pass_soft_puUp")	
		self.ZGamma_APV_f1_puUp = self.o.Get("fail_soft_puUp")	
		
		self.ZGamma_APV_p7_w_puDown = self.o.Get("jet_pt_soft_pass_wide10_wide_puDown")
		self.ZGamma_APV_f7_w_puDown = self.o.Get("jet_pt_soft_fail_wide10_wide_puDown")
		self.ZGamma_APV_p1_puDown = self.o.Get("pass_soft_puDown")	
		self.ZGamma_APV_f1_puDown = self.o.Get("fail_soft_puDown")	
		
		self.ZGamma_APV_p7_w_jerUp = self.o.Get("jet_pt_soft_pass_wide10_wide_jerUp")
		self.ZGamma_APV_f7_w_jerUp = self.o.Get("jet_pt_soft_fail_wide10_wide_jerUp")
		self.ZGamma_APV_p1_jerUp = self.o.Get("pass_soft_jerUp")	
		self.ZGamma_APV_f1_jerUp = self.o.Get("fail_soft_jerUp")	
		
		self.ZGamma_APV_p7_w_jerDown = self.o.Get("jet_pt_soft_pass_wide10_wide_jerDown")
		self.ZGamma_APV_f7_w_jerDown = self.o.Get("jet_pt_soft_fail_wide10_wide_jerDown")
		self.ZGamma_APV_p1_jerDown = self.o.Get("pass_soft_jerDown")	
		self.ZGamma_APV_f1_jerDown = self.o.Get("fail_soft_jerDown")	
		
		self.ZGamma_APV_p7_w_jesUp = self.o.Get("jet_pt_soft_pass_wide10_wide_jesUp")
		self.ZGamma_APV_f7_w_jesUp = self.o.Get("jet_pt_soft_fail_wide10_wide_jesUp")
		self.ZGamma_APV_p1_jesUp = self.o.Get("pass_soft_jesUp")	
		self.ZGamma_APV_f1_jesUp = self.o.Get("fail_soft_jesUp")	
		
		self.ZGamma_APV_p7_w_jesDown = self.o.Get("jet_pt_soft_pass_wide10_wide_jesDown")
		self.ZGamma_APV_f7_w_jesDown = self.o.Get("jet_pt_soft_fail_wide10_wide_jesDown")
		self.ZGamma_APV_p1_jesDown = self.o.Get("pass_soft_jesDown")	
		self.ZGamma_APV_f1_jesDown = self.o.Get("fail_soft_jesDown")	
		
		self.p = TFile(bfile8, "READ")
		self.GJ_APV_p7_w = self.p.Get("jet_pt_soft_pass_wide10_wide")
		self.GJ_APV_f7_w = self.p.Get("jet_pt_soft_fail_wide10_wide")
		
		self.GJ_APV_p1 = self.p.Get("pass_soft")	
		self.GJ_APV_f1 = self.p.Get("fail_soft")	
		
		self.GJ_APV_p7_w_puUp = self.p.Get("jet_pt_soft_pass_wide10_wide_puUp")
		self.GJ_APV_f7_w_puUp = self.p.Get("jet_pt_soft_fail_wide10_wide_puUp")
		self.GJ_APV_p1_puUp = self.p.Get("pass_soft_puUp")	
		self.GJ_APV_f1_puUp = self.p.Get("fail_soft_puUp")	
		
		self.GJ_APV_p7_w_puDown = self.p.Get("jet_pt_soft_pass_wide10_wide_puDown")
		self.GJ_APV_f7_w_puDown = self.p.Get("jet_pt_soft_fail_wide10_wide_puDown")
		self.GJ_APV_p1_puDown = self.p.Get("pass_soft_puDown")	
		self.GJ_APV_f1_puDown = self.p.Get("fail_soft_puDown")	
		
		self.GJ_APV_p7_w_jerUp = self.p.Get("jet_pt_soft_pass_wide10_wide_jerUp")
		self.GJ_APV_f7_w_jerUp = self.p.Get("jet_pt_soft_fail_wide10_wide_jerUp")
		self.GJ_APV_p1_jerUp = self.p.Get("pass_soft_jerUp")	
		self.GJ_APV_f1_jerUp = self.p.Get("fail_soft_jerUp")	
		
		self.GJ_APV_p7_w_jerDown = self.p.Get("jet_pt_soft_pass_wide10_wide_jerDown")
		self.GJ_APV_f7_w_jerDown = self.p.Get("jet_pt_soft_fail_wide10_wide_jerDown")
		self.GJ_APV_p1_jerDown = self.p.Get("pass_soft_jerDown")	
		self.GJ_APV_f1_jerDown = self.p.Get("fail_soft_jerDown")	
		
		self.GJ_APV_p7_w_jesUp = self.p.Get("jet_pt_soft_pass_wide10_wide_jesUp")
		self.GJ_APV_f7_w_jesUp = self.p.Get("jet_pt_soft_fail_wide10_wide_jesUp")
		self.GJ_APV_p1_jesUp = self.p.Get("pass_soft_jesUp")	
		self.GJ_APV_f1_jesUp = self.p.Get("fail_soft_jesUp")	
		
		self.GJ_APV_p7_w_jesDown = self.p.Get("jet_pt_soft_pass_wide10_wide_jesDown")
		self.GJ_APV_f7_w_jesDown = self.p.Get("jet_pt_soft_fail_wide10_wide_jesDown")
		self.GJ_APV_p1_jesDown = self.p.Get("pass_soft_jesDown")	
		self.GJ_APV_f1_jesDown = self.p.Get("fail_soft_jesDown")	

		self.s2 = TFile(sfile2, "READ")
		self.Sig50_APV_p7_w = self.s2.Get("jet_pt_soft_pass_wide10_wide")
		self.Sig50_APV_f7_w = self.s2.Get("jet_pt_soft_fail_wide10_wide")
		
		self.Sig50_APV_p1 = self.s2.Get("pass_soft")	
		self.Sig50_APV_f1 = self.s2.Get("fail_soft")	
		
		self.Sig50_APV_p7_w_puUp = self.s2.Get("jet_pt_soft_pass_wide10_wide_puUp")
		self.Sig50_APV_f7_w_puUp = self.s2.Get("jet_pt_soft_fail_wide10_wide_puUp")
		self.Sig50_APV_p1_puUp = self.s2.Get("pass_soft_puUp")	
		self.Sig50_APV_f1_puUp = self.s2.Get("fail_soft_puUp")	
		
		self.Sig50_APV_p7_w_puDown = self.s2.Get("jet_pt_soft_pass_wide10_wide_puDown")
		self.Sig50_APV_f7_w_puDown = self.s2.Get("jet_pt_soft_fail_wide10_wide_puDown")
		self.Sig50_APV_p1_puDown = self.s2.Get("pass_soft_puDown")	
		self.Sig50_APV_f1_puDown = self.s2.Get("fail_soft_puDown")	
		
		self.Sig50_APV_p7_w_jerUp = self.s2.Get("jet_pt_soft_pass_wide10_wide_jerUp")
		self.Sig50_APV_f7_w_jerUp = self.s2.Get("jet_pt_soft_fail_wide10_wide_jerUp")
		self.Sig50_APV_p1_jerUp = self.s2.Get("pass_soft_jerUp")	
		self.Sig50_APV_f1_jerUp = self.s2.Get("fail_soft_jerUp")	
		
		self.Sig50_APV_p7_w_jerDown = self.s2.Get("jet_pt_soft_pass_wide10_wide_jerDown")
		self.Sig50_APV_f7_w_jerDown = self.s2.Get("jet_pt_soft_fail_wide10_wide_jerDown")
		self.Sig50_APV_p1_jerDown = self.s2.Get("pass_soft_jerDown")	
		self.Sig50_APV_f1_jerDown = self.s2.Get("fail_soft_jerDown")	
		
		self.Sig50_APV_p7_w_jesUp = self.s2.Get("jet_pt_soft_pass_wide10_wide_jesUp")
		self.Sig50_APV_f7_w_jesUp = self.s2.Get("jet_pt_soft_fail_wide10_wide_jesUp")
		self.Sig50_APV_p1_jesUp = self.s2.Get("pass_soft_jesUp")	
		self.Sig50_APV_f1_jesUp = self.s2.Get("fail_soft_jesUp")	
		
		self.Sig50_APV_p7_w_jesDown = self.s2.Get("jet_pt_soft_pass_wide10_wide_jesDown")
		self.Sig50_APV_f7_w_jesDown = self.s2.Get("jet_pt_soft_fail_wide10_wide_jesDown")
		self.Sig50_APV_p1_jesDown = self.s2.Get("pass_soft_jesDown")	
		self.Sig50_APV_f1_jesDown = self.s2.Get("fail_soft_jesDown")	


		# Combine Non-APV and APV
		self.Data_p7_w.Add(self.Data_APV_p7_w)
		self.Data_f7_w.Add(self.Data_APV_f7_w)
		self.Data_p1.Add(self.Data_APV_p1)
		self.Data_f1.Add(self.Data_APV_f1)

		self.TTBar_p7_w.Add(self.TTBar_APV_p7_w)
		self.TTBar_f7_w.Add(self.TTBar_APV_f7_w)
		self.TTBar_p1.Add(self.TTBar_APV_p1)
		self.TTBar_f1.Add(self.TTBar_APV_f1)
		
		self.TTBar_p7_w_puUp.Add(self.TTBar_APV_p7_w_puUp)
		self.TTBar_f7_w_puUp.Add(self.TTBar_APV_f7_w_puUp)
		self.TTBar_p1_puUp.Add(self.TTBar_APV_p1_puUp)
		self.TTBar_f1_puUp.Add(self.TTBar_APV_f1_puUp)
		
		self.TTBar_p7_w_puDown.Add(self.TTBar_APV_p7_w_puDown)
		self.TTBar_f7_w_puDown.Add(self.TTBar_APV_f7_w_puDown)
		self.TTBar_p1_puDown.Add(self.TTBar_APV_p1_puDown)
		self.TTBar_f1_puDown.Add(self.TTBar_APV_f1_puDown)
		
		self.TTBar_p7_w_jerUp.Add(self.TTBar_APV_p7_w_jerUp)
		self.TTBar_f7_w_jerUp.Add(self.TTBar_APV_f7_w_jerUp)
		self.TTBar_p1_jerUp.Add(self.TTBar_APV_p1_jerUp)
		self.TTBar_f1_jerUp.Add(self.TTBar_APV_f1_jerUp)
		
		self.TTBar_p7_w_jerDown.Add(self.TTBar_APV_p7_w_jerDown)
		self.TTBar_f7_w_jerDown.Add(self.TTBar_APV_f7_w_jerDown)
		self.TTBar_p1_jerDown.Add(self.TTBar_APV_p1_jerDown)
		self.TTBar_f1_jerDown.Add(self.TTBar_APV_f1_jerDown)

		self.TTBar_p7_w_jesUp.Add(self.TTBar_APV_p7_w_jesUp)
		self.TTBar_f7_w_jesUp.Add(self.TTBar_APV_f7_w_jesUp)
		self.TTBar_p1_jesUp.Add(self.TTBar_APV_p1_jesUp)
		self.TTBar_f1_jesUp.Add(self.TTBar_APV_f1_jesUp)
		
		self.TTBar_p7_w_jesDown.Add(self.TTBar_APV_p7_w_jesDown)
		self.TTBar_f7_w_jesDown.Add(self.TTBar_APV_f7_w_jesDown)
		self.TTBar_p1_jesDown.Add(self.TTBar_APV_p1_jesDown)
		self.TTBar_f1_jesDown.Add(self.TTBar_APV_f1_jesDown)

		self.WGamma_p7_w.Add(self.WGamma_APV_p7_w)
		self.WGamma_f7_w.Add(self.WGamma_APV_f7_w)
		self.WGamma_p1.Add(self.WGamma_APV_p1)
		self.WGamma_f1.Add(self.WGamma_APV_f1)
		
		self.WGamma_p7_w_puUp.Add(self.WGamma_APV_p7_w_puUp)
		self.WGamma_f7_w_puUp.Add(self.WGamma_APV_f7_w_puUp)
		self.WGamma_p1_puUp.Add(self.WGamma_APV_p1_puUp)
		self.WGamma_f1_puUp.Add(self.WGamma_APV_f1_puUp)
		
		self.WGamma_p7_w_puDown.Add(self.WGamma_APV_p7_w_puDown)
		self.WGamma_f7_w_puDown.Add(self.WGamma_APV_f7_w_puDown)
		self.WGamma_p1_puDown.Add(self.WGamma_APV_p1_puDown)
		self.WGamma_f1_puDown.Add(self.WGamma_APV_f1_puDown)
		
		self.WGamma_p7_w_jerUp.Add(self.WGamma_APV_p7_w_jerUp)
		self.WGamma_f7_w_jerUp.Add(self.WGamma_APV_f7_w_jerUp)
		self.WGamma_p1_jerUp.Add(self.WGamma_APV_p1_jerUp)
		self.WGamma_f1_jerUp.Add(self.WGamma_APV_f1_jerUp)
		
		self.WGamma_p7_w_jerDown.Add(self.WGamma_APV_p7_w_jerDown)
		self.WGamma_f7_w_jerDown.Add(self.WGamma_APV_f7_w_jerDown)
		self.WGamma_p1_jerDown.Add(self.WGamma_APV_p1_jerDown)
		self.WGamma_f1_jerDown.Add(self.WGamma_APV_f1_jerDown)

		self.WGamma_p7_w_jesUp.Add(self.WGamma_APV_p7_w_jesUp)
		self.WGamma_f7_w_jesUp.Add(self.WGamma_APV_f7_w_jesUp)
		self.WGamma_p1_jesUp.Add(self.WGamma_APV_p1_jesUp)
		self.WGamma_f1_jesUp.Add(self.WGamma_APV_f1_jesUp)
		
		self.WGamma_p7_w_jesDown.Add(self.WGamma_APV_p7_w_jesDown)
		self.WGamma_f7_w_jesDown.Add(self.WGamma_APV_f7_w_jesDown)
		self.WGamma_p1_jesDown.Add(self.WGamma_APV_p1_jesDown)
		self.WGamma_f1_jesDown.Add(self.WGamma_APV_f1_jesDown)

		self.ZGamma_p7_w.Add(self.ZGamma_APV_p7_w)
		self.ZGamma_f7_w.Add(self.ZGamma_APV_f7_w)
		self.ZGamma_p1.Add(self.ZGamma_APV_p1)
		self.ZGamma_f1.Add(self.ZGamma_APV_f1)
		
		self.ZGamma_p7_w_puUp.Add(self.ZGamma_APV_p7_w_puUp)
		self.ZGamma_f7_w_puUp.Add(self.ZGamma_APV_f7_w_puUp)
		self.ZGamma_p1_puUp.Add(self.ZGamma_APV_p1_puUp)
		self.ZGamma_f1_puUp.Add(self.ZGamma_APV_f1_puUp)
		
		self.ZGamma_p7_w_puDown.Add(self.ZGamma_APV_p7_w_puDown)
		self.ZGamma_f7_w_puDown.Add(self.ZGamma_APV_f7_w_puDown)
		self.ZGamma_p1_puDown.Add(self.ZGamma_APV_p1_puDown)
		self.ZGamma_f1_puDown.Add(self.ZGamma_APV_f1_puDown)
		
		self.ZGamma_p7_w_jerUp.Add(self.ZGamma_APV_p7_w_jerUp)
		self.ZGamma_f7_w_jerUp.Add(self.ZGamma_APV_f7_w_jerUp)
		self.ZGamma_p1_jerUp.Add(self.ZGamma_APV_p1_jerUp)
		self.ZGamma_f1_jerUp.Add(self.ZGamma_APV_f1_jerUp)
		
		self.ZGamma_p7_w_jerDown.Add(self.ZGamma_APV_p7_w_jerDown)
		self.ZGamma_f7_w_jerDown.Add(self.ZGamma_APV_f7_w_jerDown)
		self.ZGamma_p1_jerDown.Add(self.ZGamma_APV_p1_jerDown)
		self.ZGamma_f1_jerDown.Add(self.ZGamma_APV_f1_jerDown)

		self.ZGamma_p7_w_jesUp.Add(self.ZGamma_APV_p7_w_jesUp)
		self.ZGamma_f7_w_jesUp.Add(self.ZGamma_APV_f7_w_jesUp)
		self.ZGamma_p1_jesUp.Add(self.ZGamma_APV_p1_jesUp)
		self.ZGamma_f1_jesUp.Add(self.ZGamma_APV_f1_jesUp)
		
		self.ZGamma_p7_w_jesDown.Add(self.ZGamma_APV_p7_w_jesDown)
		self.ZGamma_f7_w_jesDown.Add(self.ZGamma_APV_f7_w_jesDown)
		self.ZGamma_p1_jesDown.Add(self.ZGamma_APV_p1_jesDown)
		self.ZGamma_f1_jesDown.Add(self.ZGamma_APV_f1_jesDown)

		self.GJ_p7_w.Add(self.GJ_APV_p7_w)
		self.GJ_f7_w.Add(self.GJ_APV_f7_w)
		self.GJ_p1.Add(self.GJ_APV_p1)
		self.GJ_f1.Add(self.GJ_APV_f1)
		
		self.GJ_p7_w_puUp.Add(self.GJ_APV_p7_w_puUp)
		self.GJ_f7_w_puUp.Add(self.GJ_APV_f7_w_puUp)
		self.GJ_p1_puUp.Add(self.GJ_APV_p1_puUp)
		self.GJ_f1_puUp.Add(self.GJ_APV_f1_puUp)
		
		self.GJ_p7_w_puDown.Add(self.GJ_APV_p7_w_puDown)
		self.GJ_f7_w_puDown.Add(self.GJ_APV_f7_w_puDown)
		self.GJ_p1_puDown.Add(self.GJ_APV_p1_puDown)
		self.GJ_f1_puDown.Add(self.GJ_APV_f1_puDown)
		
		self.GJ_p7_w_jerUp.Add(self.GJ_APV_p7_w_jerUp)
		self.GJ_f7_w_jerUp.Add(self.GJ_APV_f7_w_jerUp)
		self.GJ_p1_jerUp.Add(self.GJ_APV_p1_jerUp)
		self.GJ_f1_jerUp.Add(self.GJ_APV_f1_jerUp)
		
		self.GJ_p7_w_jerDown.Add(self.GJ_APV_p7_w_jerDown)
		self.GJ_f7_w_jerDown.Add(self.GJ_APV_f7_w_jerDown)
		self.GJ_p1_jerDown.Add(self.GJ_APV_p1_jerDown)
		self.GJ_f1_jerDown.Add(self.GJ_APV_f1_jerDown)

		self.GJ_p7_w_jesUp.Add(self.GJ_APV_p7_w_jesUp)
		self.GJ_f7_w_jesUp.Add(self.GJ_APV_f7_w_jesUp)
		self.GJ_p1_jesUp.Add(self.GJ_APV_p1_jesUp)
		self.GJ_f1_jesUp.Add(self.GJ_APV_f1_jesUp)
		
		self.GJ_p7_w_jesDown.Add(self.GJ_APV_p7_w_jesDown)
		self.GJ_f7_w_jesDown.Add(self.GJ_APV_f7_w_jesDown)
		self.GJ_p1_jesDown.Add(self.GJ_APV_p1_jesDown)
		self.GJ_f1_jesDown.Add(self.GJ_APV_f1_jesDown)

		self.Sig50_p7_w.Add(self.Sig50_APV_p7_w)
		self.Sig50_f7_w.Add(self.Sig50_APV_f7_w)
		self.Sig50_p1.Add(self.Sig50_APV_p1)
		self.Sig50_f1.Add(self.Sig50_APV_f1)
		
		self.Sig50_p7_w_puUp.Add(self.Sig50_APV_p7_w_puUp)
		self.Sig50_f7_w_puUp.Add(self.Sig50_APV_f7_w_puUp)
		self.Sig50_p1_puUp.Add(self.Sig50_APV_p1_puUp)
		self.Sig50_f1_puUp.Add(self.Sig50_APV_f1_puUp)
		
		self.Sig50_p7_w_puDown.Add(self.Sig50_APV_p7_w_puDown)
		self.Sig50_f7_w_puDown.Add(self.Sig50_APV_f7_w_puDown)
		self.Sig50_p1_puDown.Add(self.Sig50_APV_p1_puDown)
		self.Sig50_f1_puDown.Add(self.Sig50_APV_f1_puDown)
		
		self.Sig50_p7_w_jerUp.Add(self.Sig50_APV_p7_w_jerUp)
		self.Sig50_f7_w_jerUp.Add(self.Sig50_APV_f7_w_jerUp)
		self.Sig50_p1_jerUp.Add(self.Sig50_APV_p1_jerUp)
		self.Sig50_f1_jerUp.Add(self.Sig50_APV_f1_jerUp)
		
		self.Sig50_p7_w_jerDown.Add(self.Sig50_APV_p7_w_jerDown)
		self.Sig50_f7_w_jerDown.Add(self.Sig50_APV_f7_w_jerDown)
		self.Sig50_p1_jerDown.Add(self.Sig50_APV_p1_jerDown)
		self.Sig50_f1_jerDown.Add(self.Sig50_APV_f1_jerDown)

		self.Sig50_p7_w_jesUp.Add(self.Sig50_APV_p7_w_jesUp)
		self.Sig50_f7_w_jesUp.Add(self.Sig50_APV_f7_w_jesUp)
		self.Sig50_p1_jesUp.Add(self.Sig50_APV_p1_jesUp)
		self.Sig50_f1_jesUp.Add(self.Sig50_APV_f1_jesUp)
		
		self.Sig50_p7_w_jesDown.Add(self.Sig50_APV_p7_w_jesDown)
		self.Sig50_f7_w_jesDown.Add(self.Sig50_APV_f7_w_jesDown)
		self.Sig50_p1_jesDown.Add(self.Sig50_APV_p1_jesDown)
		self.Sig50_f1_jesDown.Add(self.Sig50_APV_f1_jesDown)



		# Closure Test (Add all MC as "Data", "Data" is already GJets)
		self.Data_p7_w.Add(self.TTBar_p7_w)
		self.Data_p7_w.Add(self.WGamma_p7_w)
		self.Data_p7_w.Add(self.ZGamma_p7_w)

		self.Data_f7_w.Add(self.TTBar_f7_w)
		self.Data_f7_w.Add(self.WGamma_f7_w)
		self.Data_f7_w.Add(self.ZGamma_f7_w)
		
		self.Data_p1.Add(self.TTBar_p1)
		self.Data_p1.Add(self.WGamma_p1)
		self.Data_p1.Add(self.ZGamma_p1)
		
		self.Data_f1.Add(self.TTBar_f1)
		self.Data_f1.Add(self.WGamma_f1)
		self.Data_f1.Add(self.ZGamma_f1)

		# Currently using 10% of data
		self.Data_p7_w.Scale(.1)
		self.Data_f7_w.Scale(.1)
		self.Data_p1.Scale(.1)
		self.Data_f1.Scale(.1)

		self.TTBar_p7_w.Scale(.1)
		self.TTBar_f7_w.Scale(.1)
		self.TTBar_p1.Scale(.1)
		self.TTBar_f1.Scale(.1)

		self.TTBar_p7_w_puUp.Scale(.1)
		self.TTBar_f7_w_puUp.Scale(.1)
		self.TTBar_p1_puUp.Scale(.1)
		self.TTBar_f1_puUp.Scale(.1)

		self.TTBar_p7_w_puDown.Scale(.1)
		self.TTBar_f7_w_puDown.Scale(.1)
		self.TTBar_p1_puDown.Scale(.1)
		self.TTBar_f1_puDown.Scale(.1)

		self.TTBar_p7_w_jerUp.Scale(.1)
		self.TTBar_f7_w_jerUp.Scale(.1)
		self.TTBar_p1_jerUp.Scale(.1)
		self.TTBar_f1_jerUp.Scale(.1)

		self.TTBar_p7_w_jerDown.Scale(.1)
		self.TTBar_f7_w_jerDown.Scale(.1)
		self.TTBar_p1_jerDown.Scale(.1)
		self.TTBar_f1_jerDown.Scale(.1)
		
		self.TTBar_p7_w_jesUp.Scale(.1)
		self.TTBar_f7_w_jesUp.Scale(.1)
		self.TTBar_p1_jesUp.Scale(.1)
		self.TTBar_f1_jesUp.Scale(.1)

		self.TTBar_p7_w_jesDown.Scale(.1)
		self.TTBar_f7_w_jesDown.Scale(.1)
		self.TTBar_p1_jesDown.Scale(.1)
		self.TTBar_f1_jesDown.Scale(.1)


		self.WGamma_p7_w.Scale(.1)
		self.WGamma_f7_w.Scale(.1)
		
		self.WGamma_p7_w_puUp.Scale(.1)
		self.WGamma_f7_w_puUp.Scale(.1)
		self.WGamma_p1_puUp.Scale(.1)
		self.WGamma_f1_puUp.Scale(.1)

		self.WGamma_p7_w_puDown.Scale(.1)
		self.WGamma_f7_w_puDown.Scale(.1)
		self.WGamma_p1_puDown.Scale(.1)
		self.WGamma_f1_puDown.Scale(.1)

		self.WGamma_p7_w_jerUp.Scale(.1)
		self.WGamma_f7_w_jerUp.Scale(.1)
		self.WGamma_p1_jerUp.Scale(.1)
		self.WGamma_f1_jerUp.Scale(.1)

		self.WGamma_p7_w_jerDown.Scale(.1)
		self.WGamma_f7_w_jerDown.Scale(.1)
		self.WGamma_p1_jerDown.Scale(.1)
		self.WGamma_f1_jerDown.Scale(.1)
		
		self.WGamma_p7_w_jesUp.Scale(.1)
		self.WGamma_f7_w_jesUp.Scale(.1)
		self.WGamma_p1_jesUp.Scale(.1)
		self.WGamma_f1_jesUp.Scale(.1)

		self.WGamma_p7_w_jesDown.Scale(.1)
		self.WGamma_f7_w_jesDown.Scale(.1)
		self.WGamma_p1_jesDown.Scale(.1)
		self.WGamma_f1_jesDown.Scale(.1)

		self.ZGamma_p7_w.Scale(.1)
		self.ZGamma_f7_w.Scale(.1)
		
		self.ZGamma_p7_w_puUp.Scale(.1)
		self.ZGamma_f7_w_puUp.Scale(.1)
		self.ZGamma_p1_puUp.Scale(.1)
		self.ZGamma_f1_puUp.Scale(.1)

		self.ZGamma_p7_w_puDown.Scale(.1)
		self.ZGamma_f7_w_puDown.Scale(.1)
		self.ZGamma_p1_puDown.Scale(.1)
		self.ZGamma_f1_puDown.Scale(.1)

		self.ZGamma_p7_w_jerUp.Scale(.1)
		self.ZGamma_f7_w_jerUp.Scale(.1)
		self.ZGamma_p1_jerUp.Scale(.1)
		self.ZGamma_f1_jerUp.Scale(.1)

		self.ZGamma_p7_w_jerDown.Scale(.1)
		self.ZGamma_f7_w_jerDown.Scale(.1)
		self.ZGamma_p1_jerDown.Scale(.1)
		self.ZGamma_f1_jerDown.Scale(.1)
		
		self.ZGamma_p7_w_jesUp.Scale(.1)
		self.ZGamma_f7_w_jesUp.Scale(.1)
		self.ZGamma_p1_jesUp.Scale(.1)
		self.ZGamma_f1_jesUp.Scale(.1)

		self.ZGamma_p7_w_jesDown.Scale(.1)
		self.ZGamma_f7_w_jesDown.Scale(.1)
		self.ZGamma_p1_jesDown.Scale(.1)
		self.ZGamma_f1_jesDown.Scale(.1)

		self.GJ_p7_w.Scale(.1)
		self.GJ_f7_w.Scale(.1)
		
		self.GJ_p7_w_puUp.Scale(.1)
		self.GJ_f7_w_puUp.Scale(.1)
		self.GJ_p1_puUp.Scale(.1)
		self.GJ_f1_puUp.Scale(.1)

		self.GJ_p7_w_puDown.Scale(.1)
		self.GJ_f7_w_puDown.Scale(.1)
		self.GJ_p1_puDown.Scale(.1)
		self.GJ_f1_puDown.Scale(.1)

		self.GJ_p7_w_jerUp.Scale(.1)
		self.GJ_f7_w_jerUp.Scale(.1)
		self.GJ_p1_jerUp.Scale(.1)
		self.GJ_f1_jerUp.Scale(.1)

		self.GJ_p7_w_jerDown.Scale(.1)
		self.GJ_f7_w_jerDown.Scale(.1)
		self.GJ_p1_jerDown.Scale(.1)
		self.GJ_f1_jerDown.Scale(.1)
		
		self.GJ_p7_w_jesUp.Scale(.1)
		self.GJ_f7_w_jesUp.Scale(.1)
		self.GJ_p1_jesUp.Scale(.1)
		self.GJ_f1_jesUp.Scale(.1)

		self.GJ_p7_w_jesDown.Scale(.1)
		self.GJ_f7_w_jesDown.Scale(.1)
		self.GJ_p1_jesDown.Scale(.1)
		self.GJ_f1_jesDown.Scale(.1)

		self.Sig50_p7_w.Scale(.1)
		self.Sig50_f7_w.Scale(.1)

		self.Sig50_p7_w_puUp.Scale(.1)
		self.Sig50_f7_w_puUp.Scale(.1)
		self.Sig50_p1_puUp.Scale(.1)
		self.Sig50_f1_puUp.Scale(.1)

		self.Sig50_p7_w_puDown.Scale(.1)
		self.Sig50_f7_w_puDown.Scale(.1)
		self.Sig50_p1_puDown.Scale(.1)
		self.Sig50_f1_puDown.Scale(.1)

		self.Sig50_p7_w_jerUp.Scale(.1)
		self.Sig50_f7_w_jerUp.Scale(.1)
		self.Sig50_p1_jerUp.Scale(.1)
		self.Sig50_f1_jerUp.Scale(.1)

		self.Sig50_p7_w_jerDown.Scale(.1)
		self.Sig50_f7_w_jerDown.Scale(.1)
		self.Sig50_p1_jerDown.Scale(.1)
		self.Sig50_f1_jerDown.Scale(.1)
		
		self.Sig50_p7_w_jesUp.Scale(.1)
		self.Sig50_f7_w_jesUp.Scale(.1)
		self.Sig50_p1_jesUp.Scale(.1)
		self.Sig50_f1_jesUp.Scale(.1)

		self.Sig50_p7_w_jesDown.Scale(.1)
		self.Sig50_f7_w_jesDown.Scale(.1)
		self.Sig50_p1_jesDown.Scale(.1)
		self.Sig50_f1_jesDown.Scale(.1)
		


	
		
		ofile.WriteObject(self.Data_p7_w, "Data_pass_jet_pt_soft_wide10")
		ofile.WriteObject(self.Data_f7_w, "Data_fail_jet_pt_soft_wide10")
		ofile.WriteObject(self.Data_p1, "Data_pass_soft")
		ofile.WriteObject(self.Data_f1, "Data_fail_soft")
		
		
		ofile.WriteObject(self.TTBar_p7_w, "TTBar_pass_jet_pt_soft_wide10")
		ofile.WriteObject(self.TTBar_f7_w, "TTBar_fail_jet_pt_soft_wide10")
		ofile.WriteObject(self.TTBar_p1, "TTBar_pass_soft")
		ofile.WriteObject(self.TTBar_f1, "TTBar_fail_soft")
		
		ofile.WriteObject(self.TTBar_p7_w_puUp, "TTBar_pass_jet_pt_soft_wide10_puUp")
		ofile.WriteObject(self.TTBar_f7_w_puUp, "TTBar_fail_jet_pt_soft_wide10_puUp")
		ofile.WriteObject(self.TTBar_p1_puUp, "TTBar_pass_soft_puUp")
		ofile.WriteObject(self.TTBar_f1_puUp, "TTBar_fail_soft_puUp")

		ofile.WriteObject(self.TTBar_p7_w_puDown, "TTBar_pass_jet_pt_soft_wide10_puDown")
		ofile.WriteObject(self.TTBar_f7_w_puDown, "TTBar_fail_jet_pt_soft_wide10_puDown")
		ofile.WriteObject(self.TTBar_p1_puDown, "TTBar_pass_soft_puDown")
		ofile.WriteObject(self.TTBar_f1_puDown, "TTBar_fail_soft_puDown")
		
		ofile.WriteObject(self.TTBar_p7_w_jerUp, "TTBar_pass_jet_pt_soft_wide10_jerUp")
		ofile.WriteObject(self.TTBar_f7_w_jerUp, "TTBar_fail_jet_pt_soft_wide10_jerUp")
		ofile.WriteObject(self.TTBar_p1_jerUp, "TTBar_pass_soft_jerUp")
		ofile.WriteObject(self.TTBar_f1_jerUp, "TTBar_fail_soft_jerUp")

		ofile.WriteObject(self.TTBar_p7_w_jerDown, "TTBar_pass_jet_pt_soft_wide10_jerDown")
		ofile.WriteObject(self.TTBar_f7_w_jerDown, "TTBar_fail_jet_pt_soft_wide10_jerDown")
		ofile.WriteObject(self.TTBar_p1_jerDown, "TTBar_pass_soft_jerDown")
		ofile.WriteObject(self.TTBar_f1_jerDown, "TTBar_fail_soft_jerDown")
		
		ofile.WriteObject(self.TTBar_p7_w_jesUp, "TTBar_pass_jet_pt_soft_wide10_jesUp")
		ofile.WriteObject(self.TTBar_f7_w_jesUp, "TTBar_fail_jet_pt_soft_wide10_jesUp")
		ofile.WriteObject(self.TTBar_p1_jesUp, "TTBar_pass_soft_jesUp")
		ofile.WriteObject(self.TTBar_f1_jesUp, "TTBar_fail_soft_jesUp")

		ofile.WriteObject(self.TTBar_p7_w_jesDown, "TTBar_pass_jet_pt_soft_wide10_jesDown")
		ofile.WriteObject(self.TTBar_f7_w_jesDown, "TTBar_fail_jet_pt_soft_wide10_jesDown")
		ofile.WriteObject(self.TTBar_p1_jesDown, "TTBar_pass_soft_jesDown")
		ofile.WriteObject(self.TTBar_f1_jesDown, "TTBar_fail_soft_jesDown")


		ofile.WriteObject(self.WGamma_p7_w, "WGamma_pass_jet_pt_soft_wide10")
		ofile.WriteObject(self.WGamma_f7_w, "WGamma_fail_jet_pt_soft_wide10")
		ofile.WriteObject(self.WGamma_p1, "WGamma_pass_soft")
		ofile.WriteObject(self.WGamma_f1, "WGamma_fail_soft")
		
		ofile.WriteObject(self.WGamma_p7_w_puUp, "WGamma_pass_jet_pt_soft_wide10_puUp")
		ofile.WriteObject(self.WGamma_f7_w_puUp, "WGamma_fail_jet_pt_soft_wide10_puUp")
		ofile.WriteObject(self.WGamma_p1_puUp, "WGamma_pass_soft_puUp")
		ofile.WriteObject(self.WGamma_f1_puUp, "WGamma_fail_soft_puUp")

		ofile.WriteObject(self.WGamma_p7_w_puDown, "WGamma_pass_jet_pt_soft_wide10_puDown")
		ofile.WriteObject(self.WGamma_f7_w_puDown, "WGamma_fail_jet_pt_soft_wide10_puDown")
		ofile.WriteObject(self.WGamma_p1_puDown, "WGamma_pass_soft_puDown")
		ofile.WriteObject(self.WGamma_f1_puDown, "WGamma_fail_soft_puDown")
		
		ofile.WriteObject(self.WGamma_p7_w_jerUp, "WGamma_pass_jet_pt_soft_wide10_jerUp")
		ofile.WriteObject(self.WGamma_f7_w_jerUp, "WGamma_fail_jet_pt_soft_wide10_jerUp")
		ofile.WriteObject(self.WGamma_p1_jerUp, "WGamma_pass_soft_jerUp")
		ofile.WriteObject(self.WGamma_f1_jerUp, "WGamma_fail_soft_jerUp")

		ofile.WriteObject(self.WGamma_p7_w_jerDown, "WGamma_pass_jet_pt_soft_wide10_jerDown")
		ofile.WriteObject(self.WGamma_f7_w_jerDown, "WGamma_fail_jet_pt_soft_wide10_jerDown")
		ofile.WriteObject(self.WGamma_p1_jerDown, "WGamma_pass_soft_jerDown")
		ofile.WriteObject(self.WGamma_f1_jerDown, "WGamma_fail_soft_jerDown")
		
		ofile.WriteObject(self.WGamma_p7_w_jesUp, "WGamma_pass_jet_pt_soft_wide10_jesUp")
		ofile.WriteObject(self.WGamma_f7_w_jesUp, "WGamma_fail_jet_pt_soft_wide10_jesUp")
		ofile.WriteObject(self.WGamma_p1_jesUp, "WGamma_pass_soft_jesUp")
		ofile.WriteObject(self.WGamma_f1_jesUp, "WGamma_fail_soft_jesUp")

		ofile.WriteObject(self.WGamma_p7_w_jesDown, "WGamma_pass_jet_pt_soft_wide10_jesDown")
		ofile.WriteObject(self.WGamma_f7_w_jesDown, "WGamma_fail_jet_pt_soft_wide10_jesDown")
		ofile.WriteObject(self.WGamma_p1_jesDown, "WGamma_pass_soft_jesDown")
		ofile.WriteObject(self.WGamma_f1_jesDown, "WGamma_fail_soft_jesDown")



		ofile.WriteObject(self.ZGamma_p7_w, "ZGamma_pass_jet_pt_soft_wide10")
		ofile.WriteObject(self.ZGamma_f7_w, "ZGamma_fail_jet_pt_soft_wide10")
		ofile.WriteObject(self.ZGamma_p1, "ZGamma_pass_soft")
		ofile.WriteObject(self.ZGamma_f1, "ZGamma_fail_soft")
		
		ofile.WriteObject(self.ZGamma_p7_w_puUp, "ZGamma_pass_jet_pt_soft_wide10_puUp")
		ofile.WriteObject(self.ZGamma_f7_w_puUp, "ZGamma_fail_jet_pt_soft_wide10_puUp")
		ofile.WriteObject(self.ZGamma_p1_puUp, "ZGamma_pass_soft_puUp")
		ofile.WriteObject(self.ZGamma_f1_puUp, "ZGamma_fail_soft_puUp")

		ofile.WriteObject(self.ZGamma_p7_w_puDown, "ZGamma_pass_jet_pt_soft_wide10_puDown")
		ofile.WriteObject(self.ZGamma_f7_w_puDown, "ZGamma_fail_jet_pt_soft_wide10_puDown")
		ofile.WriteObject(self.ZGamma_p1_puDown, "ZGamma_pass_soft_puDown")
		ofile.WriteObject(self.ZGamma_f1_puDown, "ZGamma_fail_soft_puDown")
		
		ofile.WriteObject(self.ZGamma_p7_w_jerUp, "ZGamma_pass_jet_pt_soft_wide10_jerUp")
		ofile.WriteObject(self.ZGamma_f7_w_jerUp, "ZGamma_fail_jet_pt_soft_wide10_jerUp")
		ofile.WriteObject(self.ZGamma_p1_jerUp, "ZGamma_pass_soft_jerUp")
		ofile.WriteObject(self.ZGamma_f1_jerUp, "ZGamma_fail_soft_jerUp")

		ofile.WriteObject(self.ZGamma_p7_w_jerDown, "ZGamma_pass_jet_pt_soft_wide10_jerDown")
		ofile.WriteObject(self.ZGamma_f7_w_jerDown, "ZGamma_fail_jet_pt_soft_wide10_jerDown")
		ofile.WriteObject(self.ZGamma_p1_jerDown, "ZGamma_pass_soft_jerDown")
		ofile.WriteObject(self.ZGamma_f1_jerDown, "ZGamma_fail_soft_jerDown")
		
		ofile.WriteObject(self.ZGamma_p7_w_jesUp, "ZGamma_pass_jet_pt_soft_wide10_jesUp")
		ofile.WriteObject(self.ZGamma_f7_w_jesUp, "ZGamma_fail_jet_pt_soft_wide10_jesUp")
		ofile.WriteObject(self.ZGamma_p1_jesUp, "ZGamma_pass_soft_jesUp")
		ofile.WriteObject(self.ZGamma_f1_jesUp, "ZGamma_fail_soft_jesUp")

		ofile.WriteObject(self.ZGamma_p7_w_jesDown, "ZGamma_pass_jet_pt_soft_wide10_jesDown")
		ofile.WriteObject(self.ZGamma_f7_w_jesDown, "ZGamma_fail_jet_pt_soft_wide10_jesDown")
		ofile.WriteObject(self.ZGamma_p1_jesDown, "ZGamma_pass_soft_jesDown")
		ofile.WriteObject(self.ZGamma_f1_jesDown, "ZGamma_fail_soft_jesDown")



		ofile.WriteObject(self.GJ_p7_w, "GJ_pass_jet_pt_soft_wide10")
		ofile.WriteObject(self.GJ_f7_w, "GJ_fail_jet_pt_soft_wide10")
		ofile.WriteObject(self.GJ_p1, "GJ_pass_soft")
		ofile.WriteObject(self.GJ_f1, "GJ_fail_soft")
		
		ofile.WriteObject(self.GJ_p7_w_puUp, "GJ_pass_jet_pt_soft_wide10_puUp")
		ofile.WriteObject(self.GJ_f7_w_puUp, "GJ_fail_jet_pt_soft_wide10_puUp")
		ofile.WriteObject(self.GJ_p1_puUp, "GJ_pass_soft_puUp")
		ofile.WriteObject(self.GJ_f1_puUp, "GJ_fail_soft_puUp")

		ofile.WriteObject(self.GJ_p7_w_puDown, "GJ_pass_jet_pt_soft_wide10_puDown")
		ofile.WriteObject(self.GJ_f7_w_puDown, "GJ_fail_jet_pt_soft_wide10_puDown")
		ofile.WriteObject(self.GJ_p1_puDown, "GJ_pass_soft_puDown")
		ofile.WriteObject(self.GJ_f1_puDown, "GJ_fail_soft_puDown")
		
		ofile.WriteObject(self.GJ_p7_w_jerUp, "GJ_pass_jet_pt_soft_wide10_jerUp")
		ofile.WriteObject(self.GJ_f7_w_jerUp, "GJ_fail_jet_pt_soft_wide10_jerUp")
		ofile.WriteObject(self.GJ_p1_jerUp, "GJ_pass_soft_jerUp")
		ofile.WriteObject(self.GJ_f1_jerUp, "GJ_fail_soft_jerUp")

		ofile.WriteObject(self.GJ_p7_w_jerDown, "GJ_pass_jet_pt_soft_wide10_jerDown")
		ofile.WriteObject(self.GJ_f7_w_jerDown, "GJ_fail_jet_pt_soft_wide10_jerDown")
		ofile.WriteObject(self.GJ_p1_jerDown, "GJ_pass_soft_jerDown")
		ofile.WriteObject(self.GJ_f1_jerDown, "GJ_fail_soft_jerDown")
		
		ofile.WriteObject(self.GJ_p7_w_jesUp, "GJ_pass_jet_pt_soft_wide10_jesUp")
		ofile.WriteObject(self.GJ_f7_w_jesUp, "GJ_fail_jet_pt_soft_wide10_jesUp")
		ofile.WriteObject(self.GJ_p1_jesUp, "GJ_pass_soft_jesUp")
		ofile.WriteObject(self.GJ_f1_jesUp, "GJ_fail_soft_jesUp")

		ofile.WriteObject(self.GJ_p7_w_jesDown, "GJ_pass_jet_pt_soft_wide10_jesDown")
		ofile.WriteObject(self.GJ_f7_w_jesDown, "GJ_fail_jet_pt_soft_wide10_jesDown")
		ofile.WriteObject(self.GJ_p1_jesDown, "GJ_pass_soft_jesDown")
		ofile.WriteObject(self.GJ_f1_jesDown, "GJ_fail_soft_jesDown")


		
		ofile.WriteObject(self.Sig50_p7_w, "Sig50_pass_jet_pt_soft_wide10")
		ofile.WriteObject(self.Sig50_f7_w, "Sig50_fail_jet_pt_soft_wide10")
		ofile.WriteObject(self.Sig50_p1, "Sig50_pass_soft")
		ofile.WriteObject(self.Sig50_f1, "Sig50_fail_soft")
		
		ofile.WriteObject(self.Sig50_p7_w_puUp, "Sig50_pass_jet_pt_soft_wide10_puUp")
		ofile.WriteObject(self.Sig50_f7_w_puUp, "Sig50_fail_jet_pt_soft_wide10_puUp")
		ofile.WriteObject(self.Sig50_p1_puUp, "Sig50_pass_soft_puUp")
		ofile.WriteObject(self.Sig50_f1_puUp, "Sig50_fail_soft_puUp")

		ofile.WriteObject(self.Sig50_p7_w_puDown, "Sig50_pass_jet_pt_soft_wide10_puDown")
		ofile.WriteObject(self.Sig50_f7_w_puDown, "Sig50_fail_jet_pt_soft_wide10_puDown")
		ofile.WriteObject(self.Sig50_p1_puDown, "Sig50_pass_soft_puDown")
		ofile.WriteObject(self.Sig50_f1_puDown, "Sig50_fail_soft_puDown")
		
		ofile.WriteObject(self.Sig50_p7_w_jerUp, "Sig50_pass_jet_pt_soft_wide10_jerUp")
		ofile.WriteObject(self.Sig50_f7_w_jerUp, "Sig50_fail_jet_pt_soft_wide10_jerUp")
		ofile.WriteObject(self.Sig50_p1_jerUp, "Sig50_pass_soft_jerUp")
		ofile.WriteObject(self.Sig50_f1_jerUp, "Sig50_fail_soft_jerUp")

		ofile.WriteObject(self.Sig50_p7_w_jerDown, "Sig50_pass_jet_pt_soft_wide10_jerDown")
		ofile.WriteObject(self.Sig50_f7_w_jerDown, "Sig50_fail_jet_pt_soft_wide10_jerDown")
		ofile.WriteObject(self.Sig50_p1_jerDown, "Sig50_pass_soft_jerDown")
		ofile.WriteObject(self.Sig50_f1_jerDown, "Sig50_fail_soft_jerDown")
		
		ofile.WriteObject(self.Sig50_p7_w_jesUp, "Sig50_pass_jet_pt_soft_wide10_jesUp")
		ofile.WriteObject(self.Sig50_f7_w_jesUp, "Sig50_fail_jet_pt_soft_wide10_jesUp")
		ofile.WriteObject(self.Sig50_p1_jesUp, "Sig50_pass_soft_jesUp")
		ofile.WriteObject(self.Sig50_f1_jesUp, "Sig50_fail_soft_jesUp")

		ofile.WriteObject(self.Sig50_p7_w_jesDown, "Sig50_pass_jet_pt_soft_wide10_jesDown")
		ofile.WriteObject(self.Sig50_f7_w_jesDown, "Sig50_fail_jet_pt_soft_wide10_jesDown")
		ofile.WriteObject(self.Sig50_p1_jesDown, "Sig50_pass_soft_jesDown")
		ofile.WriteObject(self.Sig50_f1_jesDown, "Sig50_fail_soft_jesDown")
		

                ofile.Write()
