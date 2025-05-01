import ROOT
from ROOT import *
import numpy
import math

def main():
	myFunc = TF1("myFunc", "[0]*(sin([1]*sqrt(x))/x**2)+[2]", 5, 200)
	refF = ROOT.TFile.Open("diagnositcs.root", 'read')

	for i in range(0,4):
		myFunc.SetParameters(2.7,4.6,1)
#		print(myFunc(12.5))

		hist = refF.Get("2018_data_bkg_ratio_pTbin"+str(i))
	
		hist.Fit(myFunc,"","",5,200)

		parameter_0 = myFunc.GetParameter(0)
		error_0 = myFunc.GetParError(0)
		parameter_1 = myFunc.GetParameter(1)
		error_1 = myFunc.GetParError(1)
		parameter_2 = myFunc.GetParameter(2)
		error_2 = myFunc.GetParError(2)
#		parameter_3 = myFunc.GetParameter(3)
#		error_3 = myFunc.GetParError(3)
		print(f"Parameter 0: {parameter_0} +/- {error_0}")
		print(f"Parameter 1: {parameter_1} +/- {error_1}")
		print(f"Parameter 2: {parameter_2} +/- {error_2}")
#		print(f"Parameter 3: {parameter_3} +/- {error_3}")
	
		c1 = TCanvas()
		c1.cd()
		hist.Draw()
		c1.Print("./fit"+str(i)+".png")
		c1.Close()
	

if __name__ == "__main__":
	main()
