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

from drawMorph import *

#from future import division

if __name__ == "__main__":
	print("Starting Run")
	ifile1 = "./morph_templates.root"
	name = "sig115_low_2017"
	RData = drawMorph(name, ifile1)
	print("Draw Morphing Finished")
