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

from buildfile import *

#from future import division

if __name__ == "__main__":
	print("Starting Run")
	bfile1 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2017/NanoTool_UL_ParticleNet/TTBar_UL_nano_2017_merged.root"
	bfile2 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2017/NanoTool_UL_ParticleNet/WGamma_UL_nano_2017_merged.root"
	bfile3 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2017/NanoTool_UL_ParticleNet/ZGamma_UL_nano_2017_merged.root"
	bfile4 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2017/NanoTool_UL_ParticleNet/GJ_UL_2017.root"
	
	dfile1 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2017/NanoTool_UL_ParticleNet/Data_UL_2017.root" 
	
	

	#Signal Files
	sfile1 = "/users/h2/akobert/CMSSW_11_3_4/src/Interpolation/Output/SigFile_2017_105GeV.root"
	
	name = "FitHist"
	Unblind = Build(name, bfile1, bfile2, bfile3, bfile4, dfile1, sfile1)
	print("FitHist Finished")
