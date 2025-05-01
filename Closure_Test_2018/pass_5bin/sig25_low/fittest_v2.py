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

from applyTF_v2 import *

#from future import division

if __name__ == "__main__":
	print("Starting Run")
	ifile1 = "./fitDiagnosticsTest.root"
	name = "Closure_Test_2018_total_5bin_sig25_low"
	RData = drawDiagnostic(name, ifile1)
	print("draw Diagnositcs Finished")
