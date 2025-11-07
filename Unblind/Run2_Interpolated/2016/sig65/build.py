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
	bfile1 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2016/NanoTool_UL_ParticleNet/TTBar_UL_nano_2016_merged.root"
	bfile2 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2016/NanoTool_UL_ParticleNet/WGamma_UL_nano_2016_merged.root"
	bfile3 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2016/NanoTool_UL_ParticleNet/ZGamma_UL_nano_2016_merged.root"
	bfile4 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2016/NanoTool_UL_ParticleNet/GJ_UL_2016.root"

	bfile5 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2016APV/NanoTool_UL_ParticleNet/TTBar_UL_nano_2016_merged.root"
	bfile6 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2016APV/NanoTool_UL_ParticleNet/WGamma_UL_nano_2016_merged.root"
	bfile7 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2016APV/NanoTool_UL_ParticleNet/ZGamma_UL_nano_2016_merged.root"
	bfile8 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2016APV/NanoTool_UL_ParticleNet/GJ_UL_2016.root"

	
	dfile1 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2016/NanoTool_UL_ParticleNet/Data_UL_2016.root" 
	
	dfile2 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2016APV/NanoTool_UL_ParticleNet/Data_UL_2016.root" 
	

	#Signal Files
	sfile1 = "/users/h2/akobert/CMSSW_11_3_4/src/Interpolation/Output/SigFile_2016_65GeV.root"

	#sfile2 = "/home/akobert/CMSSW_11_1_0_pre7/src/Unblind/2016APV/NanoTool_UL_ParticleNet/M10_UL_nano_2016_merged.root"
	
	name = "FitHist"
	Unblind = Build(name, bfile1, bfile2, bfile3, bfile4, bfile5, bfile6, bfile7, bfile8, dfile1, dfile2, sfile1)
	print("FitHist Finished")
