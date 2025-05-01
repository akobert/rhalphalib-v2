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

from drawDiagnostics_new import *

#from future import division

if __name__ == "__main__":
	print("Starting Run")
	ifile1 = "sig25_low/fitDiagnosticsTest.root"
	name = "NEW_multiSig_low_2018_ParticleNet"
	mass = 25
	RData = drawDiagnostic(name, ifile1, mass)
	print("draw Diagnositcs Finished")
