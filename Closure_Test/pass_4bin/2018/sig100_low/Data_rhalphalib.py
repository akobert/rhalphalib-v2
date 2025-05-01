from __future__ import print_function, division
import sys
import os
import csv
import rhalphalib as rl
from rhalphalib import AffineMorphTemplate, MorphHistW2
import numpy as np
import scipy.stats
import pickle
import ROOT
from ROOT import *
rl.util.install_roofit_helpers()
#rl.ParametericSample.PreferRooParametricHist = False
import pandas as pd
import json

eps=0.001
do_systematics = True # Quick enable or disable non-lumi systematics

#class AffineMorphTemplate(object):
#    def __init__(self, hist):
#        '''                                                                                                   
#        hist: a numpy-histogram-like tuple of (sumw, edges)                                                       
#        '''
#        from scipy.interpolate import interp1d
#
#        self.sumw = hist[0]
#        self.edges = hist[1]
#        self.norm = self.sumw.sum()
#        self.mean = (self.sumw*(self.edges[:-1] + self.edges[1:])/2).sum() / self.norm
#        self.cdf = interp1d(x=self.edges,
#                            y=np.r_[0, np.cumsum(self.sumw / self.norm)],
#                            kind='linear',
#                            assume_sorted=True,
#                            bounds_error=False,
#                            fill_value=(0, 1),
#                           )
#
#    def get(self, shift=0., scale=1.):
#        '''                                                                                                             
#        Return a shifted and scaled histogram                                                                                      
#        i.e. new edges = edges * scale + shift                                                                                      
#        '''
#        scaled_edges = (self.edges - shift) / scale
#        return np.diff(self.cdf(scaled_edges)) * self.norm, self.edges

def badtemp_ma(hvalues, mask=None):
    # Need minimum size & more than 1 non-zero bins           
    tot = np.sum(hvalues[mask])
    
    count_nonzeros = np.sum(hvalues[mask] > 0)
    if (tot < eps) or (count_nonzeros < 2):
        return True
    else:
        return False

def syst_variation(numerator,denominator):
    """
    Get systematic variation relative to nominal (denominator)
    """
    var = np.divide(numerator,denominator)
    var[np.where(numerator==0)] = 1
    var[np.where(denominator==0)] = 1
    return var

def smass(sName):
    if sName in ['ggF','VBF','WH','ZH','ttH']:
        _mass = 125.
    elif sName in ['Wjets','EWKW','ttbar','singlet','VV', 'WGamma', 'TTBar']:
        _mass = 80.379
    elif sName in ['Zjets','Zjetsbb','EWKZ','EWKZbb', 'ZGamma']:
        _mass = 91.
    elif "Sig" in sName:
	_mass = float(sName[3:])
    else:
        raise ValueError("What is {}".format(sName))
    return _mass

def get_template(sName, passed, ptbin, obs, syst, muon=False):
    """
    Read msd template from root file
    """
    f = ROOT.TFile.Open('FitHist.root')

    name = 'fail'
    if passed:
        name = 'pass'

    name = sName+'_'+name+'_'+"jet_pt_soft_wide15"

    if syst != 'nominal':
	name += '_'+syst

    h = f.Get(name)
#    print("Template: "+name)
    sumw = []
    sumw2 = []
    #for i in range(1, h.GetNbinsX()+1):
    #for i in range(1, 11):
    for i in range(1, 41):
        sumw += [h.GetBinContent(ptbin, i)]
        sumw2 += [h.GetBinError(ptbin, i)*h.GetBinError(ptbin, i)]

    integral = np.array(sumw).sum()
#    print("For template "+name+" integral is "+str(integral))
    
    return (np.array(sumw), obs.binning, obs.name, np.array(sumw2))

def get_template2(sName, passed, ptbin, obs, syst, muon=False):
    """
    Read msd template from root file
    """
    f = ROOT.TFile.Open('FitHist.root')

    name = 'fail'
    if passed:
        name = 'pass'

    name = sName+'_'+name+'_'+"jet_pt_soft_wide15"

    h = f.Get(name)
#    print("Template: "+name)
    sumw = []
    sumw2 = []
    #for i in range(1, 41):
    for i in range(1, 41):
        sumw += [h.GetBinContent(ptbin, i)]
        sumw2 += [h.GetBinError(ptbin, i)*h.GetBinError(ptbin, i)]

    integral = np.array(sumw).sum()
 #   print("For template "+name+" integral is "+str(integral))
    
    return (np.array(sumw), obs.binning, obs.name, np.array(sumw2))

def shape_to_num(var, nom, clip=1.5):
    nom_rate = np.sum(nom)
    var_rate = np.sum(var)

    if abs(var_rate/nom_rate) > clip:
        var_rate = clip*nom_rate

    if var_rate < 0:
        var_rate = 0

    return var_rate/nom_rate

def num_to_hist(num, title):
    hist = TH1F("hist_conv", title, 40, 0, 200) 

    for i in range(1, num.size+1):
	hist.SetBinContent(i, num[i-1])

    return hist

def passfailSF(isPass, sName, ptbin, obs, mask, SF=1, SF_unc=0.1, muon=False):
    """
    Return (SF, SF_unc) for a pass/fail scale factor.
    """
    if isPass:
        return SF, 1. + SF_unc / SF
    else:
        _pass = get_template(sName, 1, ptbin, obs=obs, syst='nominal', muon=muon)
        _pass_rate = np.sum(_pass[0] * mask)

        _fail = get_template(sName, 0, ptbin, obs=obs, syst='nominal', muon=muon)
        _fail_rate = np.sum(_fail[0] * mask)

        if _fail_rate > 0:
            _sf = 1 + (1 - SF) * _pass_rate / _fail_rate
            _sfunc = 1. - SF_unc * (_pass_rate / _fail_rate)
            return _sf, _sfunc
        else:
            return 1, 1

def plot_mctf(tf_MCtempl, msdbins):
    """
    Plot the MC pass / fail TF as function of (pt,rho) and (pt,msd)
    """
    import matplotlib.pyplot as plt
    
    ofile = ROOT.TFile("TF_Data.root", "RECREATE")
    ofile.cd()

    # arrays for plotting pt vs msd                    
    #pts = np.linspace(0,2000,41)
    pts = np.array([0, 120, 160, 200, 255, 700 ,900, 1200, 1500, 2000])
    #ptpts, msdpts = np.meshgrid(pts[:-1] + 0.5 * np.diff(pts), msdbins[:-1] + 0.5 * np.diff(msdbins), indexing='ij')
    ptpts, msdpts = np.meshgrid(pts[:-1] + 0.5 * np.diff(pts), msdbins[:-1] + 0.5 * np.diff(msdbins), indexing='ij')
    ptpts_scaled = (ptpts) / (2000.)
    
    rhopts = 2*np.log(msdpts/ptpts)

    rhopts_scaled = (rhopts - (-7.0)) / ((-2.0) - (-7.0))
    validbins = (rhopts_scaled >= 0) & (rhopts_scaled <= 1)

    ptpts = ptpts[validbins]
    msdpts = msdpts[validbins]
    ptpts_scaled = ptpts_scaled[validbins]
    rhopts_scaled = rhopts_scaled[validbins]

    tf_MCtempl_vals = tf_MCtempl(ptpts_scaled, rhopts_scaled, nominal=True)
    df = pd.DataFrame([])
    df['msd'] = msdpts.reshape(-1)
    df['pt'] = ptpts.reshape(-1)
    df['MCTF'] = tf_MCtempl_vals.reshape(-1)
    
    mpt = df.to_numpy()
    
    ROOT.gInterpreter.Declare("Double_t widebins15[10] = {0, 120, 160, 200, 255, 700 ,900, 1200, 1500, 2000};")
    mpt_TF = TH2F("mpt_TF", "Mass vs. pT Transfer Factor", 9, widebins15, 40, 0, 2000)
    for i in range(0, mpt.shape[0]):
        mpt_TF.Fill(mpt[i][0], mpt[i][1], mpt[i][2])
#	print(mpt[i][2])
    
    ofile.WriteObject(mpt_TF, "mpt_TF")

    fig, ax = plt.subplots()
    h = ax.hist2d(x=df["msd"],y=df["pt"],weights=df["MCTF"], bins=(msdbins,pts))
    plt.xlabel("$m_{sd}$ [GeV]")
    plt.ylabel("$p_{T}$ [GeV]")
    cb = fig.colorbar(h[3],ax=ax)
    cb.set_label("Ratio")
    fig.savefig("MCTF_msdpt_2018_sig10v3.png",bbox="tight")
    plt.clf()

    # arrays for plotting pt vs rho                                          
    rhos = np.linspace(-7.0,-2.0,41)
    #ptpts, rhopts = np.meshgrid(pts[:-1] + 0.5*np.diff(pts), rhos[:-1] + 0.5 * np.diff(rhos), indexing='ij')
    ptpts, rhopts = np.meshgrid(pts[:-1] + 0.5*np.diff(pts), rhos[:-1] + 0.5 * np.diff(rhos), indexing='ij')
    ptpts_scaled = (ptpts) / (2000)
    rhopts_scaled = (rhopts - (-7.0)) / ((-2.0) - (-7.0))
    validbins = (rhopts_scaled >= 0) & (rhopts_scaled <= 1)

    ptpts = ptpts[validbins]
    rhopts = rhopts[validbins]
    ptpts_scaled = ptpts_scaled[validbins]
    rhopts_scaled = rhopts_scaled[validbins]

    tf_MCtempl_vals = tf_MCtempl(ptpts_scaled, rhopts_scaled, nominal=True)

    df = pd.DataFrame([])
    df['rho'] = rhopts.reshape(-1)
    df['pt'] = ptpts.reshape(-1)
    df['MCTF'] = tf_MCtempl_vals.reshape(-1)

    rpt = df.to_numpy()

    rpt_TF = TH2F("rpt_TF", "Rho vs. pT Transfer Factor", 40, -9, 0, 40, 0, 2000)
    for i in range(0, rpt.shape[0]):
        rpt_TF.Fill(rpt[i][0], rpt[i][1], rpt[i][2])

    
    ofile.WriteObject(rpt_TF, "rpt_TF")
    
    
    
    fig, ax = plt.subplots()
    h = ax.hist2d(x=df["rho"],y=df["pt"],weights=df["MCTF"],bins=(rhos,pts))
    plt.xlabel("rho")
    plt.ylabel("$p_{T}$ [GeV]")
    cb = fig.colorbar(h[3],ax=ax)
    cb.set_label("Ratio")
    fig.savefig("MCTF_rhopt_2018_sig10v3.png",bbox="tight")

    return

def test_rhalphabet(tmpdir,
                    throwPoisson = True,
                    fast=0):
    """ 
    Create the data cards!
    """
    with open('/home/akobert/CMSSW_11_3_4/src/FitTest/sf.json') as f:
        SF = json.load(f)
    with open('/home/akobert/CMSSW_11_3_4/src/FitTest/lumi.json') as f:
        lumi = json.load(f)

    # Xsection floating params
    tqq_xs = rl.NuisanceParameter('CMS_EXO24027_cross_section_tt', 'lnN', 1., 0, 10)
    WG_xs = rl.NuisanceParameter('CMS_EXO24027_cross_section_WG', 'lnN', 1., 0, 10)
    ZG_xs = rl.NuisanceParameter('CMS_EXO24027_cross_section_ZG', 'lnN', 1., 0, 10)

    # Systematics
    sys_lumi = rl.NuisanceParameter('lumi_13TeV_correlated', 'lnN')
    sys_phoID = rl.NuisanceParameter('CMS_EXO24027_eff_g_ID', 'lnN')
    sys_trig = rl.NuisanceParameter('CMS_EXO24027_eff_g_trigger', 'lnN') #Trigger Systematic

    sys_dict = {}
    
    sys_dict['jes'] = rl.NuisanceParameter('CMS_scale_j_{}'.format(year[:4]), 'lnN')
    sys_dict['jer'] = rl.NuisanceParameter('CMS_res_j_{}'.format(year[:4]), 'lnN')
    
    sys_dict['pu'] = rl.NuisanceParameter('CMS_pileup_{}'.format(year[:4]), 'lnN')

    sys_ParticleNet_eff = rl.NuisanceParameter('CMS_EXO24027_eff_j_PN_{}'.format(year[:4]), 'lnN')

    sys_smear = rl.NuisanceParameter('CMS_EXO24027_jmr_j_{}'.format(year[:4]), 'shape')
    sys_scale = rl.NuisanceParameter('CMS_EXO24027_jms_j_{}'.format(year[:4]), 'shape')

    # Systematic Tags
    exp_systs = ['pu', 'jes', 'jer']

    # define bins    
    ptbins = np.array([120, 160, 200, 255, 700])
    npt = len(ptbins) - 1
    msdbins = np.linspace(0, 200, 41)
    msd = rl.Observable('msd', msdbins)

    # here we derive these all at once with 2D array
    ptpts, msdpts = np.meshgrid(ptbins[:-1] + 0.5 * np.diff(ptbins), msdbins[:-1] + 0.5 * np.diff(msdbins), indexing='ij')
    print(ptpts)
    print(msdpts)
    rhopts = 2*np.log(msdpts/ptpts)
    print(rhopts)
    ptscaled = (ptpts - 120.0) / (700.0 - 120.0)
    rhoscaled = (rhopts - (-7.0)) / ((-2.0) - (-7.0))
    validbins = (rhoscaled >= 0) & (rhoscaled <= 1)
    rhoscaled[~validbins] = 1  # we will mask these out later

    # Build MC pass+fail model and fit to polynomial
    GJmodel = rl.Model('GJmodel')
    GJpass, GJfail = 0., 0.
    pass_count, fail_count = 0., 0.
    for ptbin in range(npt):
        failCh = rl.Channel('ptbin%d%s' % (ptbin, 'fail'))
        passCh = rl.Channel('ptbin%d%s' % (ptbin, 'pass'))
        GJmodel.addChannel(failCh)
        GJmodel.addChannel(passCh)

        # GJ templates from file
        failTempl = get_template('GJ', 0, ptbin+2, obs=msd, syst='nominal') 
        passTempl = get_template('GJ', 1, ptbin+2, obs=msd, syst='nominal') 
        

        failCh.setObservation(failTempl, read_sumw2=True)
        passCh.setObservation(passTempl, read_sumw2=True)


        GJfail += sum([val for val in failCh.getObservation()[0]])
        fail_count += sum([val for val in failCh.getObservation()[0]])

	print(GJfail)
	print(fail_count)
        GJpass += sum([val for val in passCh.getObservation()[0]])
        pass_count += sum([val for val in passCh.getObservation()[0]])
	print(GJpass)
	print(pass_count)

    GJeff = GJpass / GJfail
    print('Inclusive P/F from Monte Carlo = ' + str(GJeff))

    # initial values
    print('Initial fit values read from file initial_vals.csv')
    initial_vals = np.genfromtxt('initial_vals.csv')
    initial_vals = initial_vals.reshape(3, 3)
    print(initial_vals)

    tf_MCtempl = rl.BasisPoly('CMS_EXO24027_tf_MCtempl', (2, 2), ['pt', 'rho'], basis='Bernstein', init_params=initial_vals, limits=(-100, 100))
    tf_MCtempl_params = GJeff * tf_MCtempl(ptscaled, rhoscaled)
    for ptbin in range(npt):
        failCh = GJmodel['ptbin%dfail' % ptbin]
        passCh = GJmodel['ptbin%dpass' % ptbin]
        failObs = failCh.getObservation()
        passObs = passCh.getObservation()


	print(np.size(failObs))

	if np.any(failObs < 0):
	    for i in range(np.size(failObs[0])):
		if failObs[0][i] < 0:
		   failObs[0][i] = 0
		   failObs[1][i] = 0
	if np.any(passObs < 0):
	    for i in range(np.size(passObs[0])):
		if passObs[0][i] < 0:
		   passObs[0][i] = 0
		   passObs[1][i] = 0


        GJparams = np.array([rl.IndependentParameter('GJparam_ptbin%d_msdbin%d' % (ptbin, i), 0) for i in range(msd.nbins)])
        sigmascale = 10.
        scaledparams = failObs * (1 + sigmascale/np.maximum(1., np.sqrt(failObs)))**GJparams

        fail_GJ = rl.ParametericSample('ptbin%dfail_GJ' % ptbin, rl.Sample.BACKGROUND, msd, scaledparams[0])
        failCh.addSample(fail_GJ)
        pass_GJ = rl.TransferFactorSample('ptbin%dpass_GJ' % ptbin, rl.Sample.BACKGROUND, tf_MCtempl_params[ptbin, :], fail_GJ)
        passCh.addSample(pass_GJ)

        failCh.mask = validbins[ptbin]
        passCh.mask = validbins[ptbin]

    GJfit_ws = ROOT.RooWorkspace('GJfit_ws')

    simpdf, obs = GJmodel.renderRoofit(GJfit_ws)
    GJfit = simpdf.fitTo(obs,
                          ROOT.RooFit.Extended(True),
                          ROOT.RooFit.SumW2Error(True),
                          ROOT.RooFit.Strategy(2),
        		  ROOT.RooFit.Offset(True),
	                  ROOT.RooFit.Save(),
                          ROOT.RooFit.Minimizer('Minuit2', 'migrad'),
                          ROOT.RooFit.PrintLevel(1),
                          )
    GJfit_ws.add(GJfit)
    GJfit_ws.writeToFile(os.path.join(str(tmpdir), 'testModel_GJfit_2018.root'))

	

    # Set parameters to fitted values  
    allparams = dict(zip(GJfit.nameArray(), GJfit.valueArray()))
    for i, p in enumerate(tf_MCtempl.parameters.reshape(-1)):
        p.value = allparams[p.name]
        print(p.name,p.value)

    if GJfit.status() != 0:
        raise RuntimeError('Could not fit GJ')
    
     #Plot the MC P/F transfer factor
    msdbins = np.linspace(0, 200, 41) 
#    plot_mctf(tf_MCtempl,msdbins)
   
    #New Stuff
    # Redefine bins
    ptbins = np.array([120, 160, 200, 255, 700])
    npt = len(ptbins) - 1
    msdbins = np.linspace(0, 200, 41)
    msd = rl.Observable('msd', msdbins)

    # here we derive these all at once with 2D array
    ptpts, msdpts = np.meshgrid(ptbins[:-1] + 0.5 * np.diff(ptbins), msdbins[:-1] + 0.5 * np.diff(msdbins), indexing='ij')
    rhopts = 2*np.log(msdpts/ptpts)
    ptscaled = (ptpts - 120.0) / (700.0 - 120.0)
    rhoscaled = (rhopts - (-7.0)) / ((-2.0) - (-7.0))
    validbins = (rhoscaled >= 0) & (rhoscaled <= 1)
    rhoscaled[~validbins] = 1  # we will mask these out later

    param_names = [p.name for p in tf_MCtempl.parameters.reshape(-1)]
    decoVector = rl.DecorrelatedNuisanceVector.fromRooFitResult(tf_MCtempl.name + '_deco_2018_', GJfit, param_names)
    tf_MCtempl.parameters = decoVector.correlated_params.reshape(tf_MCtempl.parameters.shape)
    tf_MCtempl_params_final = tf_MCtempl(ptscaled, rhoscaled)
    tf_dataResidual = rl.BasisPoly('CMS_EXO24027_tf_dataResidual_2018', (2, 2), ['pt', 'rho'] ,basis='Bernstein', limits=(-100, 100))
    tf_dataResidual_params = tf_dataResidual(ptscaled, rhoscaled)
    tf_params = GJeff * tf_MCtempl_params_final * tf_dataResidual_params

    # build actual fit model now
    model = rl.Model('testModel_2018')

    ofile_morph = ROOT.TFile("morph_templates.root", "RECREATE")
    ofile_morph.cd()
    
    # exclude GJ from MC samps
    samps = ['TTBar','Sig100', 'WGamma', 'ZGamma']
    sigs = ['Sig100']
    for ptbin in range(npt):
        for region in ['pass', 'fail']:

            # drop bins outside rho validity                                                            
            mask = validbins[ptbin]

            ch = rl.Channel('ptbin%d%s' % (ptbin, region))
            model.addChannel(ch)

            isPass = region == 'pass'
            templates = {}

            for sName in samps:
		print("Sample "+region+" "+sName)
                templates[sName] = get_template2(sName, isPass, ptbin+2, obs=msd, syst='nominal')
                nominal = templates[sName][0]

                if(badtemp_ma(nominal)):
                    print("Sample {} is too small, skipping".format(sName))
                    continue

                # expectations
                templ = templates[sName]
                stype = rl.Sample.SIGNAL if sName in sigs else rl.Sample.BACKGROUND

		# Shape Systematics
		MORPHNOMINAL = True
                def smorph(templ):      
                    if templ is None:
                    	return None                  
                            
                    if MORPHNOMINAL and sName not in ['GJ']:
                        return MorphHistW2(templ).get(shift=SF[year]['shift_SF']/smass('Wjets') * smass(sName),
                                                      smear=SF[year]['smear_SF']
                                                  )
                    else:
                        return templ
                templ = smorph(templ)

                sample = rl.TemplateSample(ch.name + '_' + sName, stype, templ)
		print(sName+" "+region+" Expectation: "+str(sample.getExpectation(nominal=True).sum()))

		# Lumi Systematic
		sample.setParamEffect(sys_lumi, lumi['Run2']['UL'])
		sample.setParamEffect(sys_phoID, 1.05) #5% Photon ID Systematic (All Years)
		sample.setParamEffect(sys_trig, 1.05) #5% Trigger Systematic (All Years)

		# Other Systematics
		if do_systematics:
		    sample.autoMCStats(lnN=True)	# Sets MC Statistical Uncertainties
		    for sys in exp_systs:
			syst_up = get_template(sName, isPass, ptbin+2, obs=msd, syst=sys+'Up')[0]
			syst_do = get_template(sName, isPass, ptbin+2, obs=msd, syst=sys+'Down')[0]

			eff_up = shape_to_num(syst_up, nominal)
			eff_do = shape_to_num(syst_do, nominal)

			sample.setParamEffect(sys_dict[sys], eff_up, eff_do)

		# Scale and Smear
		mtempl = AffineMorphTemplate(templ)
                if sName not in ['QCD']:
                    # shift
                    realshift = SF[year]['shift_SF_ERR']/smass('Wjets') * smass(sName)
		    print("realshift: "+str(realshift))
                    _up = mtempl.get(shift=realshift)
                    _down = mtempl.get(shift=-realshift)
                    if badtemp_ma(_up[0]) or badtemp_ma(_down[0]):
                        print("Skipping sample {}, scale systematic would be empty".format(sName))
                    else:
                        sample.setParamEffect(sys_scale, _up, _down, scale=1)
			ofile_morph.WriteObject(num_to_hist(templates[sName][0], sName+" Nominal ptBin "+str(ptbin)),sName+"_nom_"+region+str(ptbin))
			ofile_morph.WriteObject(num_to_hist(templates[sName][0]*_up[0], sName+" ScaleUp ptBin "+str(ptbin)),sName+"_scaleUp_"+region+str(ptbin))
			ofile_morph.WriteObject(num_to_hist(templates[sName][0]*_down[0], sName+" ScaleDown ptBin "+str(ptbin)),sName+"_scaleDown_"+region+str(ptbin))
			ofile_morph.WriteObject(num_to_hist(_up[0], sName+" ScaleUp Weights ptBin "+str(ptbin)),sName+"_scaleUp_weights_"+region+str(ptbin))
			ofile_morph.WriteObject(num_to_hist(_down[0], sName+" ScaleDown Weights ptBin "+str(ptbin)),sName+"_scaleDown_weights_"+region+str(ptbin))

                    # smear
                    _up = mtempl.get(smear=1 + SF[year]['smear_SF_ERR'])
                    _down = mtempl.get(smear=1 - SF[year]['smear_SF_ERR'])
                    if badtemp_ma(_up[0]) or badtemp_ma(_down[0]):
                        print("Skipping sample {}, scale systematic would be empty".format(sName))
                    else:
                        sample.setParamEffect(sys_smear, _up, _down)
			ofile_morph.WriteObject(num_to_hist(templates[sName][0]*_up[0], sName+" SmearUp ptbin "+str(ptbin)),sName+"_smearUp_"+region+str(ptbin))
			ofile_morph.WriteObject(num_to_hist(templates[sName][0]*_down[0], sName+" SmearDown ptbin "+str(ptbin)),sName+"_smearDown_"+region+str(ptbin))
			ofile_morph.WriteObject(num_to_hist(_up[0], sName+" SmearUp Weights ptBin "+str(ptbin)),sName+"_smearUp_weights_"+region+str(ptbin))
			ofile_morph.WriteObject(num_to_hist(_down[0], sName+" SmearDown Weights ptBin "+str(ptbin)),sName+"_smearDown_weights_"+region+str(ptbin))

		# ParticleNet SF
		sample.scale(SF[year]['ParticleNet_SF'])
		if do_systematics:
                    effect = 1.0 + SF[year]['ParticleNet_SF_ERR'] / SF[year]['ParticleNet_SF']
                    sample.setParamEffect(sys_ParticleNet_eff,effect)

                ch.addSample(sample)

            data_obs = get_template2('Data', isPass, ptbin+2, obs=msd, syst='nominal')
            ch.setObservation(data_obs, read_sumw2=True)
    

            # drop bins outside rho validity
            mask = validbins[ptbin]

            ch.mask = mask

    for ptbin in range(npt):
        failCh = model['ptbin%dfail' % ptbin]
        passCh = model['ptbin%dpass' % ptbin]
        Dataparams = np.array([rl.IndependentParameter('CMS_EXO24027_Dataparam_ptbin%d_msdbin%d_2018' % (ptbin, i), 0) for i in range(msd.nbins)])
        initial_Data = failCh.getObservation()[0].astype(float)  # was integer, and numpy complained about subtracting float from it

        for sample in failCh:
	    if sample.sampletype == rl.Sample.BACKGROUND:#Subtract out backgrounds from fail
            	initial_Data -= sample.getExpectation(nominal=True)
		print("Subtract Sample: "+str(sample.name))


        if np.any(initial_Data < 0.):
#            raise ValueError('initial_Data negative for some bins..', initial_Data)
	    for i in range(np.size(initial_Data)):
		if initial_Data[i] < 0:
		   initial_Data[i] = 0
	    

        sigmascale = 10  # to scale the deviation from initial
        scaledparams = initial_Data * (1 + sigmascale/np.maximum(1., np.sqrt(initial_Data)))**Dataparams
        fail_NonRes = rl.ParametericSample('ptbin%dfail_NonRes2018' % ptbin, rl.Sample.BACKGROUND, msd, scaledparams)
        failCh.addSample(fail_NonRes)
        pass_NonRes = rl.TransferFactorSample('ptbin%dpass_NonRes2018' % ptbin, rl.Sample.BACKGROUND, tf_params[ptbin, :], fail_NonRes)
        passCh.addSample(pass_NonRes)

        WGpass = passCh['WGamma']
        WGfail = failCh['WGamma']
        WGpass.setParamEffect(WG_xs, 1.1)
        WGfail.setParamEffect(WG_xs, 1.1)

        ZGpass = passCh['ZGamma']
        ZGfail = failCh['ZGamma']
        ZGpass.setParamEffect(ZG_xs, 1.1)
        ZGfail.setParamEffect(ZG_xs, 1.1)

        tqqpass = passCh['TTBar']
        tqqfail = failCh['TTBar']
        tqqpass.setParamEffect(tqq_xs, 1.1)
        tqqfail.setParamEffect(tqq_xs, 1.1)

    with open(os.path.join(str(tmpdir), 'testModel_2018.pkl'), 'wb') as fout:
        pickle.dump(model, fout)

    model.renderCombine(os.path.join(str(tmpdir), 'testModel_2018'))
    
 

if __name__ == '__main__':

    year = sys.argv[1]
    print("Year: "+str(year))
    if not os.path.exists('output'):
        os.mkdir('output')

    test_rhalphabet('output',year)
    
    print("Rhalphalib Done")
