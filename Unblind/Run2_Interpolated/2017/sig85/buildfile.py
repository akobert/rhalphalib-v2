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
	def __init__(self, name, bfile1, bfile2, bfile3, bfile4, dfile1, sfile1):
		gROOT.SetBatch(True)

		#Output Files
	        ofile = ROOT.TFile(name + ".root", "RECREATE")
	        ofile.cd()


		#background files
		self.f = TFile.Open(dfile1, "READ")
		self.f.ls();

		self.Data_p7_w = self.f.Get("jet_pt_soft_pass_wide11_wide")
		self.Data_f7_w = self.f.Get("jet_pt_soft_fail_wide11_wide")
		self.Data_p1 = self.f.Get("pass_soft")
		self.Data_f1 = self.f.Get("fail_soft")
		

		self.g = TFile(bfile1, "READ")
		self.TTBar_p7_w = self.g.Get("jet_pt_soft_pass_wide11_wide")
		self.TTBar_f7_w = self.g.Get("jet_pt_soft_fail_wide11_wide")
		
		self.TTBar_p1 = self.g.Get("pass_soft")	
		self.TTBar_f1 = self.g.Get("fail_soft")	
		
		self.TTBar_p7_w_puUp = self.g.Get("jet_pt_soft_pass_wide11_wide_puUp")
		self.TTBar_f7_w_puUp = self.g.Get("jet_pt_soft_fail_wide11_wide_puUp")
		self.TTBar_p1_puUp = self.g.Get("pass_soft_puUp")	
		self.TTBar_f1_puUp = self.g.Get("fail_soft_puUp")	
		
		self.TTBar_p7_w_puDown = self.g.Get("jet_pt_soft_pass_wide11_wide_puDown")
		self.TTBar_f7_w_puDown = self.g.Get("jet_pt_soft_fail_wide11_wide_puDown")
		self.TTBar_p1_puDown = self.g.Get("pass_soft_puDown")	
		self.TTBar_f1_puDown = self.g.Get("fail_soft_puDown")	
		
		self.TTBar_p7_w_jerUp = self.g.Get("jet_pt_soft_pass_wide11_wide_jerUp")
		self.TTBar_f7_w_jerUp = self.g.Get("jet_pt_soft_fail_wide11_wide_jerUp")
		self.TTBar_p1_jerUp = self.g.Get("pass_soft_jerUp")	
		self.TTBar_f1_jerUp = self.g.Get("fail_soft_jerUp")	
		
		self.TTBar_p7_w_jerDown = self.g.Get("jet_pt_soft_pass_wide11_wide_jerDown")
		self.TTBar_f7_w_jerDown = self.g.Get("jet_pt_soft_fail_wide11_wide_jerDown")
		self.TTBar_p1_jerDown = self.g.Get("pass_soft_jerDown")	
		self.TTBar_f1_jerDown = self.g.Get("fail_soft_jerDown")	
		
		self.TTBar_p7_w_jesUp = self.g.Get("jet_pt_soft_pass_wide11_wide_jesUp")
		self.TTBar_f7_w_jesUp = self.g.Get("jet_pt_soft_fail_wide11_wide_jesUp")
		self.TTBar_p1_jesUp = self.g.Get("pass_soft_jesUp")	
		self.TTBar_f1_jesUp = self.g.Get("fail_soft_jesUp")	
		
		self.TTBar_p7_w_jesDown = self.g.Get("jet_pt_soft_pass_wide11_wide_jesDown")
		self.TTBar_f7_w_jesDown = self.g.Get("jet_pt_soft_fail_wide11_wide_jesDown")
		self.TTBar_p1_jesDown = self.g.Get("pass_soft_jesDown")	
		self.TTBar_f1_jesDown = self.g.Get("fail_soft_jesDown")	
		
		self.h = TFile(bfile2, "READ")
		self.WGamma_p7_w = self.h.Get("jet_pt_soft_pass_wide11_wide")
		self.WGamma_f7_w = self.h.Get("jet_pt_soft_fail_wide11_wide")
		
		self.WGamma_p1 = self.h.Get("pass_soft")	
		self.WGamma_f1 = self.h.Get("fail_soft")	
		
		self.WGamma_p7_w_puUp = self.h.Get("jet_pt_soft_pass_wide11_wide_puUp")
		self.WGamma_f7_w_puUp = self.h.Get("jet_pt_soft_fail_wide11_wide_puUp")
		self.WGamma_p1_puUp = self.h.Get("pass_soft_puUp")	
		self.WGamma_f1_puUp = self.h.Get("fail_soft_puUp")	
		
		self.WGamma_p7_w_puDown = self.h.Get("jet_pt_soft_pass_wide11_wide_puDown")
		self.WGamma_f7_w_puDown = self.h.Get("jet_pt_soft_fail_wide11_wide_puDown")
		self.WGamma_p1_puDown = self.h.Get("pass_soft_puDown")	
		self.WGamma_f1_puDown = self.h.Get("fail_soft_puDown")	
		
		self.WGamma_p7_w_jerUp = self.h.Get("jet_pt_soft_pass_wide11_wide_jerUp")
		self.WGamma_f7_w_jerUp = self.h.Get("jet_pt_soft_fail_wide11_wide_jerUp")
		self.WGamma_p1_jerUp = self.h.Get("pass_soft_jerUp")	
		self.WGamma_f1_jerUp = self.h.Get("fail_soft_jerUp")	
		
		self.WGamma_p7_w_jerDown = self.h.Get("jet_pt_soft_pass_wide11_wide_jerDown")
		self.WGamma_f7_w_jerDown = self.h.Get("jet_pt_soft_fail_wide11_wide_jerDown")
		self.WGamma_p1_jerDown = self.h.Get("pass_soft_jerDown")	
		self.WGamma_f1_jerDown = self.h.Get("fail_soft_jerDown")	
		
		self.WGamma_p7_w_jesUp = self.h.Get("jet_pt_soft_pass_wide11_wide_jesUp")
		self.WGamma_f7_w_jesUp = self.h.Get("jet_pt_soft_fail_wide11_wide_jesUp")
		self.WGamma_p1_jesUp = self.h.Get("pass_soft_jesUp")	
		self.WGamma_f1_jesUp = self.h.Get("fail_soft_jesUp")	
		
		self.WGamma_p7_w_jesDown = self.h.Get("jet_pt_soft_pass_wide11_wide_jesDown")
		self.WGamma_f7_w_jesDown = self.h.Get("jet_pt_soft_fail_wide11_wide_jesDown")
		self.WGamma_p1_jesDown = self.h.Get("pass_soft_jesDown")	
		self.WGamma_f1_jesDown = self.h.Get("fail_soft_jesDown")	
		
		self.j = TFile(bfile3, "READ")
		self.ZGamma_p7_w = self.j.Get("jet_pt_soft_pass_wide11_wide")
		self.ZGamma_f7_w = self.j.Get("jet_pt_soft_fail_wide11_wide")
		
		self.ZGamma_p1 = self.j.Get("pass_soft")	
		self.ZGamma_f1 = self.j.Get("fail_soft")	
		
		self.ZGamma_p7_w_puUp = self.j.Get("jet_pt_soft_pass_wide11_wide_puUp")
		self.ZGamma_f7_w_puUp = self.j.Get("jet_pt_soft_fail_wide11_wide_puUp")
		self.ZGamma_p1_puUp = self.j.Get("pass_soft_puUp")	
		self.ZGamma_f1_puUp = self.j.Get("fail_soft_puUp")	
		
		self.ZGamma_p7_w_puDown = self.j.Get("jet_pt_soft_pass_wide11_wide_puDown")
		self.ZGamma_f7_w_puDown = self.j.Get("jet_pt_soft_fail_wide11_wide_puDown")
		self.ZGamma_p1_puDown = self.j.Get("pass_soft_puDown")	
		self.ZGamma_f1_puDown = self.j.Get("fail_soft_puDown")	
		
		self.ZGamma_p7_w_jerUp = self.j.Get("jet_pt_soft_pass_wide11_wide_jerUp")
		self.ZGamma_f7_w_jerUp = self.j.Get("jet_pt_soft_fail_wide11_wide_jerUp")
		self.ZGamma_p1_jerUp = self.j.Get("pass_soft_jerUp")	
		self.ZGamma_f1_jerUp = self.j.Get("fail_soft_jerUp")	
		
		self.ZGamma_p7_w_jerDown = self.j.Get("jet_pt_soft_pass_wide11_wide_jerDown")
		self.ZGamma_f7_w_jerDown = self.j.Get("jet_pt_soft_fail_wide11_wide_jerDown")
		self.ZGamma_p1_jerDown = self.j.Get("pass_soft_jerDown")	
		self.ZGamma_f1_jerDown = self.j.Get("fail_soft_jerDown")	
		
		self.ZGamma_p7_w_jesUp = self.j.Get("jet_pt_soft_pass_wide11_wide_jesUp")
		self.ZGamma_f7_w_jesUp = self.j.Get("jet_pt_soft_fail_wide11_wide_jesUp")
		self.ZGamma_p1_jesUp = self.j.Get("pass_soft_jesUp")	
		self.ZGamma_f1_jesUp = self.j.Get("fail_soft_jesUp")	
		
		self.ZGamma_p7_w_jesDown = self.j.Get("jet_pt_soft_pass_wide11_wide_jesDown")
		self.ZGamma_f7_w_jesDown = self.j.Get("jet_pt_soft_fail_wide11_wide_jesDown")
		self.ZGamma_p1_jesDown = self.j.Get("pass_soft_jesDown")	
		self.ZGamma_f1_jesDown = self.j.Get("fail_soft_jesDown")	
		
		self.k = TFile(bfile4, "READ")
		self.GJ_p7_w = self.k.Get("jet_pt_soft_pass_wide11_wide")
		self.GJ_f7_w = self.k.Get("jet_pt_soft_fail_wide11_wide")
		
		self.GJ_p1 = self.k.Get("pass_soft")	
		self.GJ_f1 = self.k.Get("fail_soft")	
		
		self.GJ_p7_w_puUp = self.k.Get("jet_pt_soft_pass_wide11_wide_puUp")
		self.GJ_f7_w_puUp = self.k.Get("jet_pt_soft_fail_wide11_wide_puUp")
		self.GJ_p1_puUp = self.k.Get("pass_soft_puUp")	
		self.GJ_f1_puUp = self.k.Get("fail_soft_puUp")	
		
		self.GJ_p7_w_puDown = self.k.Get("jet_pt_soft_pass_wide11_wide_puDown")
		self.GJ_f7_w_puDown = self.k.Get("jet_pt_soft_fail_wide11_wide_puDown")
		self.GJ_p1_puDown = self.k.Get("pass_soft_puDown")	
		self.GJ_f1_puDown = self.k.Get("fail_soft_puDown")	
		
		self.GJ_p7_w_jerUp = self.k.Get("jet_pt_soft_pass_wide11_wide_jerUp")
		self.GJ_f7_w_jerUp = self.k.Get("jet_pt_soft_fail_wide11_wide_jerUp")
		self.GJ_p1_jerUp = self.k.Get("pass_soft_jerUp")	
		self.GJ_f1_jerUp = self.k.Get("fail_soft_jerUp")	
		
		self.GJ_p7_w_jerDown = self.k.Get("jet_pt_soft_pass_wide11_wide_jerDown")
		self.GJ_f7_w_jerDown = self.k.Get("jet_pt_soft_fail_wide11_wide_jerDown")
		self.GJ_p1_jerDown = self.k.Get("pass_soft_jerDown")	
		self.GJ_f1_jerDown = self.k.Get("fail_soft_jerDown")	
		
		self.GJ_p7_w_jesUp = self.k.Get("jet_pt_soft_pass_wide11_wide_jesUp")
		self.GJ_f7_w_jesUp = self.k.Get("jet_pt_soft_fail_wide11_wide_jesUp")
		self.GJ_p1_jesUp = self.k.Get("pass_soft_jesUp")	
		self.GJ_f1_jesUp = self.k.Get("fail_soft_jesUp")	
		
		self.GJ_p7_w_jesDown = self.k.Get("jet_pt_soft_pass_wide11_wide_jesDown")
		self.GJ_f7_w_jesDown = self.k.Get("jet_pt_soft_fail_wide11_wide_jesDown")
		self.GJ_p1_jesDown = self.k.Get("pass_soft_jesDown")	
		self.GJ_f1_jesDown = self.k.Get("fail_soft_jesDown")	

		self.s1 = TFile(sfile1, "READ")
		self.Sig85_p7_w = self.s1.Get("jet_pt_soft_pass_wide11_wide")
		self.Sig85_f7_w = self.s1.Get("jet_pt_soft_fail_wide11_wide")
		
		self.Sig85_p1 = self.s1.Get("pass_soft")	
		self.Sig85_f1 = self.s1.Get("fail_soft")	
		
		self.Sig85_p7_w_puUp = self.s1.Get("jet_pt_soft_pass_wide11_wide_puUp")
		self.Sig85_f7_w_puUp = self.s1.Get("jet_pt_soft_fail_wide11_wide_puUp")
		self.Sig85_p1_puUp = self.s1.Get("pass_soft_puUp")	
		self.Sig85_f1_puUp = self.s1.Get("fail_soft_puUp")	
		
		self.Sig85_p7_w_puDown = self.s1.Get("jet_pt_soft_pass_wide11_wide_puDown")
		self.Sig85_f7_w_puDown = self.s1.Get("jet_pt_soft_fail_wide11_wide_puDown")
		self.Sig85_p1_puDown = self.s1.Get("pass_soft_puDown")	
		self.Sig85_f1_puDown = self.s1.Get("fail_soft_puDown")	
		
		self.Sig85_p7_w_jerUp = self.s1.Get("jet_pt_soft_pass_wide11_wide_jerUp")
		self.Sig85_f7_w_jerUp = self.s1.Get("jet_pt_soft_fail_wide11_wide_jerUp")
		self.Sig85_p1_jerUp = self.s1.Get("pass_soft_jerUp")	
		self.Sig85_f1_jerUp = self.s1.Get("fail_soft_jerUp")	
		
		self.Sig85_p7_w_jerDown = self.s1.Get("jet_pt_soft_pass_wide11_wide_jerDown")
		self.Sig85_f7_w_jerDown = self.s1.Get("jet_pt_soft_fail_wide11_wide_jerDown")
		self.Sig85_p1_jerDown = self.s1.Get("pass_soft_jerDown")	
		self.Sig85_f1_jerDown = self.s1.Get("fail_soft_jerDown")	
		
		self.Sig85_p7_w_jesUp = self.s1.Get("jet_pt_soft_pass_wide11_wide_jesUp")
		self.Sig85_f7_w_jesUp = self.s1.Get("jet_pt_soft_fail_wide11_wide_jesUp")
		self.Sig85_p1_jesUp = self.s1.Get("pass_soft_jesUp")	
		self.Sig85_f1_jesUp = self.s1.Get("fail_soft_jesUp")	
		
		self.Sig85_p7_w_jesDown = self.s1.Get("jet_pt_soft_pass_wide11_wide_jesDown")
		self.Sig85_f7_w_jesDown = self.s1.Get("jet_pt_soft_fail_wide11_wide_jesDown")
		self.Sig85_p1_jesDown = self.s1.Get("pass_soft_jesDown")	
		self.Sig85_f1_jesDown = self.s1.Get("fail_soft_jesDown")	

		
		ofile.WriteObject(self.Data_p7_w, "Data_pass_jet_pt_soft_wide11")
		ofile.WriteObject(self.Data_f7_w, "Data_fail_jet_pt_soft_wide11")
		ofile.WriteObject(self.Data_p1, "Data_pass_soft")
		ofile.WriteObject(self.Data_f1, "Data_fail_soft")
		
		
		ofile.WriteObject(self.TTBar_p7_w, "TTBar_pass_jet_pt_soft_wide11")
		ofile.WriteObject(self.TTBar_f7_w, "TTBar_fail_jet_pt_soft_wide11")
		ofile.WriteObject(self.TTBar_p1, "TTBar_pass_soft")
		ofile.WriteObject(self.TTBar_f1, "TTBar_fail_soft")
		
		ofile.WriteObject(self.TTBar_p7_w_puUp, "TTBar_pass_jet_pt_soft_wide11_puUp")
		ofile.WriteObject(self.TTBar_f7_w_puUp, "TTBar_fail_jet_pt_soft_wide11_puUp")
		ofile.WriteObject(self.TTBar_p1_puUp, "TTBar_pass_soft_puUp")
		ofile.WriteObject(self.TTBar_f1_puUp, "TTBar_fail_soft_puUp")

		ofile.WriteObject(self.TTBar_p7_w_puDown, "TTBar_pass_jet_pt_soft_wide11_puDown")
		ofile.WriteObject(self.TTBar_f7_w_puDown, "TTBar_fail_jet_pt_soft_wide11_puDown")
		ofile.WriteObject(self.TTBar_p1_puDown, "TTBar_pass_soft_puDown")
		ofile.WriteObject(self.TTBar_f1_puDown, "TTBar_fail_soft_puDown")
		
		ofile.WriteObject(self.TTBar_p7_w_jerUp, "TTBar_pass_jet_pt_soft_wide11_jerUp")
		ofile.WriteObject(self.TTBar_f7_w_jerUp, "TTBar_fail_jet_pt_soft_wide11_jerUp")
		ofile.WriteObject(self.TTBar_p1_jerUp, "TTBar_pass_soft_jerUp")
		ofile.WriteObject(self.TTBar_f1_jerUp, "TTBar_fail_soft_jerUp")

		ofile.WriteObject(self.TTBar_p7_w_jerDown, "TTBar_pass_jet_pt_soft_wide11_jerDown")
		ofile.WriteObject(self.TTBar_f7_w_jerDown, "TTBar_fail_jet_pt_soft_wide11_jerDown")
		ofile.WriteObject(self.TTBar_p1_jerDown, "TTBar_pass_soft_jerDown")
		ofile.WriteObject(self.TTBar_f1_jerDown, "TTBar_fail_soft_jerDown")
		
		ofile.WriteObject(self.TTBar_p7_w_jesUp, "TTBar_pass_jet_pt_soft_wide11_jesUp")
		ofile.WriteObject(self.TTBar_f7_w_jesUp, "TTBar_fail_jet_pt_soft_wide11_jesUp")
		ofile.WriteObject(self.TTBar_p1_jesUp, "TTBar_pass_soft_jesUp")
		ofile.WriteObject(self.TTBar_f1_jesUp, "TTBar_fail_soft_jesUp")

		ofile.WriteObject(self.TTBar_p7_w_jesDown, "TTBar_pass_jet_pt_soft_wide11_jesDown")
		ofile.WriteObject(self.TTBar_f7_w_jesDown, "TTBar_fail_jet_pt_soft_wide11_jesDown")
		ofile.WriteObject(self.TTBar_p1_jesDown, "TTBar_pass_soft_jesDown")
		ofile.WriteObject(self.TTBar_f1_jesDown, "TTBar_fail_soft_jesDown")


		ofile.WriteObject(self.WGamma_p7_w, "WGamma_pass_jet_pt_soft_wide11")
		ofile.WriteObject(self.WGamma_f7_w, "WGamma_fail_jet_pt_soft_wide11")
		ofile.WriteObject(self.WGamma_p1, "WGamma_pass_soft")
		ofile.WriteObject(self.WGamma_f1, "WGamma_fail_soft")
		
		ofile.WriteObject(self.WGamma_p7_w_puUp, "WGamma_pass_jet_pt_soft_wide11_puUp")
		ofile.WriteObject(self.WGamma_f7_w_puUp, "WGamma_fail_jet_pt_soft_wide11_puUp")
		ofile.WriteObject(self.WGamma_p1_puUp, "WGamma_pass_soft_puUp")
		ofile.WriteObject(self.WGamma_f1_puUp, "WGamma_fail_soft_puUp")

		ofile.WriteObject(self.WGamma_p7_w_puDown, "WGamma_pass_jet_pt_soft_wide11_puDown")
		ofile.WriteObject(self.WGamma_f7_w_puDown, "WGamma_fail_jet_pt_soft_wide11_puDown")
		ofile.WriteObject(self.WGamma_p1_puDown, "WGamma_pass_soft_puDown")
		ofile.WriteObject(self.WGamma_f1_puDown, "WGamma_fail_soft_puDown")
		
		ofile.WriteObject(self.WGamma_p7_w_jerUp, "WGamma_pass_jet_pt_soft_wide11_jerUp")
		ofile.WriteObject(self.WGamma_f7_w_jerUp, "WGamma_fail_jet_pt_soft_wide11_jerUp")
		ofile.WriteObject(self.WGamma_p1_jerUp, "WGamma_pass_soft_jerUp")
		ofile.WriteObject(self.WGamma_f1_jerUp, "WGamma_fail_soft_jerUp")

		ofile.WriteObject(self.WGamma_p7_w_jerDown, "WGamma_pass_jet_pt_soft_wide11_jerDown")
		ofile.WriteObject(self.WGamma_f7_w_jerDown, "WGamma_fail_jet_pt_soft_wide11_jerDown")
		ofile.WriteObject(self.WGamma_p1_jerDown, "WGamma_pass_soft_jerDown")
		ofile.WriteObject(self.WGamma_f1_jerDown, "WGamma_fail_soft_jerDown")
		
		ofile.WriteObject(self.WGamma_p7_w_jesUp, "WGamma_pass_jet_pt_soft_wide11_jesUp")
		ofile.WriteObject(self.WGamma_f7_w_jesUp, "WGamma_fail_jet_pt_soft_wide11_jesUp")
		ofile.WriteObject(self.WGamma_p1_jesUp, "WGamma_pass_soft_jesUp")
		ofile.WriteObject(self.WGamma_f1_jesUp, "WGamma_fail_soft_jesUp")

		ofile.WriteObject(self.WGamma_p7_w_jesDown, "WGamma_pass_jet_pt_soft_wide11_jesDown")
		ofile.WriteObject(self.WGamma_f7_w_jesDown, "WGamma_fail_jet_pt_soft_wide11_jesDown")
		ofile.WriteObject(self.WGamma_p1_jesDown, "WGamma_pass_soft_jesDown")
		ofile.WriteObject(self.WGamma_f1_jesDown, "WGamma_fail_soft_jesDown")



		ofile.WriteObject(self.ZGamma_p7_w, "ZGamma_pass_jet_pt_soft_wide11")
		ofile.WriteObject(self.ZGamma_f7_w, "ZGamma_fail_jet_pt_soft_wide11")
		ofile.WriteObject(self.ZGamma_p1, "ZGamma_pass_soft")
		ofile.WriteObject(self.ZGamma_f1, "ZGamma_fail_soft")
		
		ofile.WriteObject(self.ZGamma_p7_w_puUp, "ZGamma_pass_jet_pt_soft_wide11_puUp")
		ofile.WriteObject(self.ZGamma_f7_w_puUp, "ZGamma_fail_jet_pt_soft_wide11_puUp")
		ofile.WriteObject(self.ZGamma_p1_puUp, "ZGamma_pass_soft_puUp")
		ofile.WriteObject(self.ZGamma_f1_puUp, "ZGamma_fail_soft_puUp")

		ofile.WriteObject(self.ZGamma_p7_w_puDown, "ZGamma_pass_jet_pt_soft_wide11_puDown")
		ofile.WriteObject(self.ZGamma_f7_w_puDown, "ZGamma_fail_jet_pt_soft_wide11_puDown")
		ofile.WriteObject(self.ZGamma_p1_puDown, "ZGamma_pass_soft_puDown")
		ofile.WriteObject(self.ZGamma_f1_puDown, "ZGamma_fail_soft_puDown")
		
		ofile.WriteObject(self.ZGamma_p7_w_jerUp, "ZGamma_pass_jet_pt_soft_wide11_jerUp")
		ofile.WriteObject(self.ZGamma_f7_w_jerUp, "ZGamma_fail_jet_pt_soft_wide11_jerUp")
		ofile.WriteObject(self.ZGamma_p1_jerUp, "ZGamma_pass_soft_jerUp")
		ofile.WriteObject(self.ZGamma_f1_jerUp, "ZGamma_fail_soft_jerUp")

		ofile.WriteObject(self.ZGamma_p7_w_jerDown, "ZGamma_pass_jet_pt_soft_wide11_jerDown")
		ofile.WriteObject(self.ZGamma_f7_w_jerDown, "ZGamma_fail_jet_pt_soft_wide11_jerDown")
		ofile.WriteObject(self.ZGamma_p1_jerDown, "ZGamma_pass_soft_jerDown")
		ofile.WriteObject(self.ZGamma_f1_jerDown, "ZGamma_fail_soft_jerDown")
		
		ofile.WriteObject(self.ZGamma_p7_w_jesUp, "ZGamma_pass_jet_pt_soft_wide11_jesUp")
		ofile.WriteObject(self.ZGamma_f7_w_jesUp, "ZGamma_fail_jet_pt_soft_wide11_jesUp")
		ofile.WriteObject(self.ZGamma_p1_jesUp, "ZGamma_pass_soft_jesUp")
		ofile.WriteObject(self.ZGamma_f1_jesUp, "ZGamma_fail_soft_jesUp")

		ofile.WriteObject(self.ZGamma_p7_w_jesDown, "ZGamma_pass_jet_pt_soft_wide11_jesDown")
		ofile.WriteObject(self.ZGamma_f7_w_jesDown, "ZGamma_fail_jet_pt_soft_wide11_jesDown")
		ofile.WriteObject(self.ZGamma_p1_jesDown, "ZGamma_pass_soft_jesDown")
		ofile.WriteObject(self.ZGamma_f1_jesDown, "ZGamma_fail_soft_jesDown")



		ofile.WriteObject(self.GJ_p7_w, "GJ_pass_jet_pt_soft_wide11")
		ofile.WriteObject(self.GJ_f7_w, "GJ_fail_jet_pt_soft_wide11")
		ofile.WriteObject(self.GJ_p1, "GJ_pass_soft")
		ofile.WriteObject(self.GJ_f1, "GJ_fail_soft")
		
		ofile.WriteObject(self.GJ_p7_w_puUp, "GJ_pass_jet_pt_soft_wide11_puUp")
		ofile.WriteObject(self.GJ_f7_w_puUp, "GJ_fail_jet_pt_soft_wide11_puUp")
		ofile.WriteObject(self.GJ_p1_puUp, "GJ_pass_soft_puUp")
		ofile.WriteObject(self.GJ_f1_puUp, "GJ_fail_soft_puUp")

		ofile.WriteObject(self.GJ_p7_w_puDown, "GJ_pass_jet_pt_soft_wide11_puDown")
		ofile.WriteObject(self.GJ_f7_w_puDown, "GJ_fail_jet_pt_soft_wide11_puDown")
		ofile.WriteObject(self.GJ_p1_puDown, "GJ_pass_soft_puDown")
		ofile.WriteObject(self.GJ_f1_puDown, "GJ_fail_soft_puDown")
		
		ofile.WriteObject(self.GJ_p7_w_jerUp, "GJ_pass_jet_pt_soft_wide11_jerUp")
		ofile.WriteObject(self.GJ_f7_w_jerUp, "GJ_fail_jet_pt_soft_wide11_jerUp")
		ofile.WriteObject(self.GJ_p1_jerUp, "GJ_pass_soft_jerUp")
		ofile.WriteObject(self.GJ_f1_jerUp, "GJ_fail_soft_jerUp")

		ofile.WriteObject(self.GJ_p7_w_jerDown, "GJ_pass_jet_pt_soft_wide11_jerDown")
		ofile.WriteObject(self.GJ_f7_w_jerDown, "GJ_fail_jet_pt_soft_wide11_jerDown")
		ofile.WriteObject(self.GJ_p1_jerDown, "GJ_pass_soft_jerDown")
		ofile.WriteObject(self.GJ_f1_jerDown, "GJ_fail_soft_jerDown")
		
		ofile.WriteObject(self.GJ_p7_w_jesUp, "GJ_pass_jet_pt_soft_wide11_jesUp")
		ofile.WriteObject(self.GJ_f7_w_jesUp, "GJ_fail_jet_pt_soft_wide11_jesUp")
		ofile.WriteObject(self.GJ_p1_jesUp, "GJ_pass_soft_jesUp")
		ofile.WriteObject(self.GJ_f1_jesUp, "GJ_fail_soft_jesUp")

		ofile.WriteObject(self.GJ_p7_w_jesDown, "GJ_pass_jet_pt_soft_wide11_jesDown")
		ofile.WriteObject(self.GJ_f7_w_jesDown, "GJ_fail_jet_pt_soft_wide11_jesDown")
		ofile.WriteObject(self.GJ_p1_jesDown, "GJ_pass_soft_jesDown")
		ofile.WriteObject(self.GJ_f1_jesDown, "GJ_fail_soft_jesDown")


		
		ofile.WriteObject(self.Sig85_p7_w, "Sig85_pass_jet_pt_soft_wide11")
		ofile.WriteObject(self.Sig85_f7_w, "Sig85_fail_jet_pt_soft_wide11")
		ofile.WriteObject(self.Sig85_p1, "Sig85_pass_soft")
		ofile.WriteObject(self.Sig85_f1, "Sig85_fail_soft")
		
		ofile.WriteObject(self.Sig85_p7_w_puUp, "Sig85_pass_jet_pt_soft_wide11_puUp")
		ofile.WriteObject(self.Sig85_f7_w_puUp, "Sig85_fail_jet_pt_soft_wide11_puUp")
		ofile.WriteObject(self.Sig85_p1_puUp, "Sig85_pass_soft_puUp")
		ofile.WriteObject(self.Sig85_f1_puUp, "Sig85_fail_soft_puUp")

		ofile.WriteObject(self.Sig85_p7_w_puDown, "Sig85_pass_jet_pt_soft_wide11_puDown")
		ofile.WriteObject(self.Sig85_f7_w_puDown, "Sig85_fail_jet_pt_soft_wide11_puDown")
		ofile.WriteObject(self.Sig85_p1_puDown, "Sig85_pass_soft_puDown")
		ofile.WriteObject(self.Sig85_f1_puDown, "Sig85_fail_soft_puDown")
		
		ofile.WriteObject(self.Sig85_p7_w_jerUp, "Sig85_pass_jet_pt_soft_wide11_jerUp")
		ofile.WriteObject(self.Sig85_f7_w_jerUp, "Sig85_fail_jet_pt_soft_wide11_jerUp")
		ofile.WriteObject(self.Sig85_p1_jerUp, "Sig85_pass_soft_jerUp")
		ofile.WriteObject(self.Sig85_f1_jerUp, "Sig85_fail_soft_jerUp")

		ofile.WriteObject(self.Sig85_p7_w_jerDown, "Sig85_pass_jet_pt_soft_wide11_jerDown")
		ofile.WriteObject(self.Sig85_f7_w_jerDown, "Sig85_fail_jet_pt_soft_wide11_jerDown")
		ofile.WriteObject(self.Sig85_p1_jerDown, "Sig85_pass_soft_jerDown")
		ofile.WriteObject(self.Sig85_f1_jerDown, "Sig85_fail_soft_jerDown")
		
		ofile.WriteObject(self.Sig85_p7_w_jesUp, "Sig85_pass_jet_pt_soft_wide11_jesUp")
		ofile.WriteObject(self.Sig85_f7_w_jesUp, "Sig85_fail_jet_pt_soft_wide11_jesUp")
		ofile.WriteObject(self.Sig85_p1_jesUp, "Sig85_pass_soft_jesUp")
		ofile.WriteObject(self.Sig85_f1_jesUp, "Sig85_fail_soft_jesUp")

		ofile.WriteObject(self.Sig85_p7_w_jesDown, "Sig85_pass_jet_pt_soft_wide11_jesDown")
		ofile.WriteObject(self.Sig85_f7_w_jesDown, "Sig85_fail_jet_pt_soft_wide11_jesDown")
		ofile.WriteObject(self.Sig85_p1_jesDown, "Sig85_pass_soft_jesDown")
		ofile.WriteObject(self.Sig85_f1_jesDown, "Sig85_fail_soft_jesDown")
		

                ofile.Write()
