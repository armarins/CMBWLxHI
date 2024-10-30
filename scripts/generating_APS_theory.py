#sys.path.insert(1, os.getcwd())
import os, sys, time
from copy import deepcopy as dcopy
import cross_functions_theory as cxft
import numpy as np
#import pandas as pd
import camb
from   camb import model, initialpower
import json
import argparse
import warnings
warnings.filterwarnings("ignore")	
#####################################################################################################
#####################################################################################################

#####################################################################################################
# Check the python version and import configparser
#####################################################################################################
if sys.version_info[0]==2:
	import ConfigParser
	config = ConfigParser.RawConfigParser()
elif sys.version_info[0]==3:
	import configparser
	config = configparser.ConfigParser()

#####################################################################################################
###################################################################
# This part is for extracting information from parameters.ini file
###################################################################
timei       = time.time()
INI         = "generating_APS_theory.ini"
name_params = os.path.join(os.getcwd(),INI)
config.read(name_params)
#General
verbose = config.getboolean("General","verbose")
#Cosmology
ombh2   = config.getfloat("Cosmology","ombh2")
omch2   = config.getfloat("Cosmology","omch2")
H0      = config.getfloat("Cosmology","H0")
ns      = config.getfloat("Cosmology","ns")
#CAMB
zmin_integration_camb   = config.getfloat(  "CAMB","zmin_integration_camb")
zmax_integration_camb   = config.getfloat(  "CAMB","zmax_integration_camb")
zsteps_integration_camb = config.getint(    "CAMB","zsteps_integration_camb")
kmax                    = config.getfloat(  "CAMB","kmax")
minkh                   = config.getfloat(  "CAMB","minkh")
maxkh                   = config.getfloat(  "CAMB","maxkh")
npoints                 = config.getint(    "CAMB","npoints")
Pk_nonlinear            = config.getboolean("CAMB","Pk_nonlinear")
#survey
#binning_used       = config.get("survey","binning_used")
freq_min           = config.getfloat("survey","freq_min")
freq_max           = config.getfloat("survey","freq_max")
nchannels          = config.getint(  "survey","nchannels")
#window_function_HI = config.get(     "survey","window_function_HI")
#APS
field_1             = config.get(     "APS","field_1")
field_2             = config.get(     "APS","field_2")
l_min               = config.getint(  "APS","l_min")
l_max               = config.getint(  "APS","l_max")
generate_cl_field_1 = config.getboolean("APS","generate_cl_field_1")
generate_cl_field_2 = config.getboolean("APS","generate_cl_field_2")
generate_cl_cross   = config.getboolean("APS","generate_cl_cross")
pathout             = config.get("APS","pathout")
prefix              = config.get("APS","prefix" )
suffix              = config.get("APS","suffix" )
#21cm
omegaHI_model = config.get("21cm","OmegaHI_model")
biasHI_model  = config.get("21cm","biasHI_model")
#Figures
save_clplot_field1  = config.get("Figures","save_clplot_field1")
save_clplot_field2  = config.get("Figures","save_clplot_field2")
save_clplot_cross   = config.get("Figures","save_clplot_cross")
save_cls_together   = config.get("Figures","save_cls_together")

###############################################################################
# You can modify any options in the parameters.ini file by the command terminal
###############################################################################
parser = argparse.ArgumentParser(description='Modify by the command terminal parameters in {} file'.format(INI))

#General
parser.add_argument('--verbose'  , action = 'store', dest = 'verbose'  , default = verbose  , help = '')

#Cosmology
parser.add_argument('--ombh2'  , action = 'store', dest = 'ombh2'  , default = ombh2  , help = '')
parser.add_argument('--omch2'  , action = 'store', dest = 'omch2'  , default = omch2  , help = '')
parser.add_argument('--H0'     , action = 'store', dest = 'H0'     , default = H0     , help = '')
parser.add_argument('--ns'     , action = 'store', dest = 'ns'     , default = ns     , help = '')

#CAMB
parser.add_argument('--zmin_integration_camb'  , action = 'store', dest = 'zmin_integration_camb'  , default = zmin_integration_camb  , help = '')
parser.add_argument('--zmax_integration_camb'  , action = 'store', dest = 'zmax_integration_camb'  , default = zmax_integration_camb  , help = '')
parser.add_argument('--zsteps_integration_camb', action = 'store', dest = 'zsteps_integration_camb', default = zsteps_integration_camb, help = '')
parser.add_argument('--kmax'        , action = 'store', dest = 'kmax'        , default = kmax        , help = '')
parser.add_argument('--minkh'       , action = 'store', dest = 'minkh'       , default = minkh       , help = '')
parser.add_argument('--maxkh'       , action = 'store', dest = 'maxkh'       , default = maxkh       , help = '')
parser.add_argument('--npoints'     , action = 'store', dest = 'npoints'     , default = npoints     , help = '')
parser.add_argument('--Pk_nonlinear', action = 'store', dest = 'Pk_nonlinear', default = Pk_nonlinear, help = '')

#Survey
parser.add_argument('--freq_min' , action = 'store', dest = 'freq_min' , default = freq_min , help = '')
parser.add_argument('--freq_max' , action = 'store', dest = 'freq_max' , default = freq_max , help = '')
parser.add_argument('--nchannels', action = 'store', dest = 'nchannels', default = nchannels, help = '')

#APS
parser.add_argument('--field_1' , action = 'store', dest = 'field_1' , default = field_1 , help = '')
parser.add_argument('--field_2' , action = 'store', dest = 'field_2' , default = field_2 , help = '')
parser.add_argument('--l_min'   , action = 'store', dest = 'l_min'   , default = l_min   , help = '')
parser.add_argument('--l_max'   , action = 'store', dest = 'l_max'   , default = l_max   , help = '')
parser.add_argument('--generate_cl_field_1', action = 'store', dest = 'generate_cl_field_1', default = generate_cl_field_1, help = '')
parser.add_argument('--generate_cl_field_2', action = 'store', dest = 'generate_cl_field_2', default = generate_cl_field_2, help = '')
parser.add_argument('--generate_cl_cross'  , action = 'store', dest = 'generate_cl_cross'  , default = generate_cl_cross  , help = '')
parser.add_argument('--pathout', action = 'store', dest = 'pathout', default = pathout, help = '')
parser.add_argument('--prefix' , action = 'store', dest = 'prefix' , default = prefix , help = '')
parser.add_argument('--suffix' , action = 'store', dest = 'suffix' , default = suffix , help = '')

#21cm
parser.add_argument('--omegaHI_model', action = 'store', dest = 'omegaHI_model', default = omegaHI_model, help = '')
parser.add_argument('--biasHI_model' , action = 'store', dest = 'biasHI_model' , default = biasHI_model , help = '')

#Figrues
parser.add_argument('--save_clplot_field1', action = 'store', dest = 'save_clplot_field1', default = save_clplot_field1, help = '')
parser.add_argument('--save_clplot_field2', action = 'store', dest = 'save_clplot_field2', default = save_clplot_field2, help = '')
parser.add_argument('--save_clplot_cross' , action = 'store', dest = 'save_clplot_cross' , default = save_clplot_cross , help = '')
parser.add_argument('--save_cls_together' , action = 'store', dest = 'save_cls_together' , default = save_cls_together , help = '')

arguments = parser.parse_args()
###############################################################################
# Variables
###############################################################################
verbose   = bool(arguments.verbose)

H0        = float(arguments.H0)
ns        = float(arguments.ns)
ombh2     = float(arguments.ombh2)
omch2     = float(arguments.omch2)

zmin_integration_camb   = float(arguments.zmin_integration_camb)
zmax_integration_camb   = float(arguments.zmax_integration_camb)
zsteps_integration_camb = int(arguments.zsteps_integration_camb)
kmax         = float(arguments.kmax)
minkh        = float(arguments.minkh)
maxkh        = float(arguments.maxkh)
npoints      = int(arguments.npoints)
Pk_nonlinear = bool(arguments.Pk_nonlinear)

##

freq_min  = float(arguments.freq_min)
freq_max  = float(arguments.freq_max)
nchannels = int(arguments.nchannels)

field_1  = str(arguments.field_1)
field_2  = str(arguments.field_2)
l_min    = float(arguments.l_min)
l_max    = float(arguments.l_max)
generate_cl_field_1 = bool(arguments.generate_cl_field_1)
generate_cl_field_2 = bool(arguments.generate_cl_field_2)
generate_cl_cross   = bool(arguments.generate_cl_cross)
pathout = str(arguments.pathout)
prefix  = str(arguments.prefix)
suffix  = str(arguments.suffix)

omegaHI_model = str(arguments.omegaHI_model)
biasHI_model  = str(arguments.biasHI_model)


save_clplot_field1 = bool(arguments.save_clplot_field1)
save_clplot_field2 = bool(arguments.save_clplot_field2)
save_clplot_cross  = bool(arguments.save_clplot_cross )
save_cls_together  = bool(arguments.save_cls_together )
############################################################################################
#####################################################################################################
if omegaHI_model=="":
    omegaHI_model=None
if biasHI_model=="":
    biasHI_model=None
if prefix=="":
    prefix=None
if suffix=="":
    suffix=None


#print(suffix, type(suffix), len(suffix))
#sys.exit(0)
#####################################################################################################################################################################
#####################################################################################################################################################################
params_general = {'verbose':verbose}

params = { 'H0'     : H0,
           'ombh2'  : ombh2,
           'omch2'  : omch2,
           'ns'     : ns,
           'z'      : np.logspace(np.log10(zmin_integration_camb),np.log10(zmax_integration_camb), zsteps_integration_camb),
           'kmax'   : kmax,
           'minkh'  : minkh,
           'maxkh'  : maxkh,
           'npoints': npoints,
           'Pk_nonlinear' :Pk_nonlinear,
           'omegaHI_model':omegaHI_model,
           'biasHI_model' :biasHI_model}          


params_survey = {"freq_min":freq_min, "freq_max":freq_max, "nchannels":nchannels}
params_APS    = {'field_1'            :field_1, 
                 'field_2'            :field_2, 
                 'l_min'              :l_min, 
                 'l_max'              :l_max, 
                 'generate_cl_field_1':generate_cl_field_1, 
                 'generate_cl_field_2':generate_cl_field_2, 
                 'generate_cl_cross'  :generate_cl_cross, 
                 'pathout'            :pathout, 
                 'prefix'             :prefix,
                 'suffix'             :suffix,
                 'save_clplot_field1' :save_clplot_field1,
                 'save_clplot_field2' :save_clplot_field2,                 
                 'save_clplot_cross'  :save_clplot_cross,
                 'save_cls_together'  :save_cls_together                                  
                }

del H0, ombh2, omch2, ns, kmax, minkh, maxkh, npoints, Pk_nonlinear, omegaHI_model, biasHI_model
del freq_min, freq_max, field_1, field_2, nchannels, l_min, l_max, generate_cl_field_1, generate_cl_field_2, generate_cl_cross, pathout, prefix, suffix, save_clplot_field1, save_clplot_field2, save_clplot_cross, save_cls_together
#####################################################################################################################################################################
#####################################################################################################################################################################
pars = camb.CAMBparams()
pars.set_cosmology(H0=params['H0'], ombh2=params['ombh2'], omch2=params['omch2'])
pars.InitPower.set_params(ns=params['ns'])
pars.set_matter_power(redshifts=params['z'], kmax=params['kmax'])
results   = camb.get_results(pars)
pk_interp = camb.get_matter_power_interpolator(pars, 
                                               nonlinear=params['Pk_nonlinear'], hubble_units=False, k_hunit=False, 
                                               kmax=params['kmax'], zmax=params['z'][-1])
#####################################################################################################
#####################################################################################################
c        = 299792458/1e3
nu_HI    = 1420.405751768 #MHz
nu_vec   = np.linspace(params_survey['freq_min'], params_survey['freq_max'], params_survey['nchannels'] + 1)
z_survey = np.sort(np.flip((nu_HI/nu_vec)-1))

l      = np.arange(params_APS['l_min'], params_APS['l_max']+1, 1)	
#--> binning_l in delta_l

fields          = np.array([params_APS['field_1'].lower(), params_APS['field_2'].lower()])
generate_fields = np.array([params_APS['generate_cl_field_1'], params_APS['generate_cl_field_2']])
#####################################################################################################
#####################################################################################################
if params_general['verbose']:
    print('\n------------------------------------------')
    print('Starting to produce Angular Power Spectrum')
    print('------------------------------------------\n')

plot_together_f1 = 0
plot_together_f2 = 0
plot_together_fx = 0

#if 'cmbwl' in [params_APS['field_1'].lower(), params_APS['field_2'].lower()]:
ind0 = np.where('cmbwl' == fields)[0]
if len(np.where(generate_fields[ind0]==True)[0]):
    if params_general['verbose']: print('Generating CMB-convergence Weak Lensing Ang Power Spectrum...')
    from scipy import integrate
    for i, il in enumerate(l):
        chi_dic   = cxft.chi_vec(results)
        dCkk_dchi = lambda x: pk_interp.P(chi_dic['z_chi'](x),(il)/x)*(cxft.kernel_cmb_chi(camb_params=pars, camb_results=results, chi=x,norm=True)/x)**2
        ckk = integrate.quad(dCkk_dchi, 0, chi_dic['chi_star'])[0]
        if not i:Ckk=ckk
        else: Ckk = np.hstack((Ckk, ckk))
    if params_general['verbose']: print('Generated.\nSaving...')
    import creating_dir
    creating_dir.verification_dir(params_APS['pathout'].split('/')[-1], "/".join(( params_APS['pathout'].split('/')[:-1] )))
    #creating_dir.verification_dir("auxiliary_data",params_APS['pathout'])
    pathname = cxft.savedata(l_=l, Cl_=np.asmatrix(Ckk), filename='CMBWL_cl', prefix=params_APS['prefix'], suffix=params_APS['suffix'], path=params_APS['pathout'], header="[1] l, [2] CMB weak lensing Cl")
    if params_general['verbose']: print('Saved at {}'.format(pathname))
    if params_APS['save_clplot_field1']:
        if params_general['verbose']: 
            print('Saving plot...')
        clf1_dict = cxft.cls_from_matrix2dict(np.vstack((l,Ckk)) )
        pathname  = cxft.save_plots(params_=params_APS, pathname_=pathname, clf_dict=clf1_dict, field='cmbwl')
        if params_general['verbose']: 
            print('Saved at {}'.format(pathname))
    if not params_APS['save_cls_together']: 
        del clf1_dict
    else: 
        plot_together_f1 = 1
    del chi_dic, dCkk_dchi, ckk, pathname, Ckk 
    print('\n')

ind0 = np.where('hi' == fields)[0]
if len(np.where(generate_fields[ind0]==True)[0]):
    if params_general['verbose']: print('Generating HI Ang Power Spectrum...')
    from scipy import integrate
    for i in np.arange(z_survey.size-1):
        zmin= z_survey[i]
        zmax= z_survey[i+1]
        for j, jl in enumerate(l):
            Whi     = lambda zi: cxft.kernel_hi(camb_params=pars, camb_results=results, z=zi, zmin = zmin, zmax=zmax)
            dChi_dz = lambda zi: (results.hubble_parameter(zi)/c)*pk_interp.P(zi,(jl+0.5)/results.comoving_radial_distance(zi))*(Whi(zi)**2/results.comoving_radial_distance(zi)**2)
            c_hi = integrate.quad(dChi_dz, zmin, zmax)[0]
            if not j:C_hi=c_hi
            else: C_hi = np.hstack((C_hi, c_hi))   
        if not i: Chi_bins = C_hi
        else:     Chi_bins = np.vstack(( Chi_bins, C_hi ))
    if params_general['verbose']: print('Generated.\nSaving...')
    import creating_dir
    creating_dir.verification_dir(params_APS['pathout'].split('/')[-1], "/".join(( params_APS['pathout'].split('/')[:-1] )))
    #creating_dir.verification_dir("auxiliary_data",params_APS['pathout'])
    pathname = cxft.savedata(l_=l, Cl_=np.asmatrix(Chi_bins), filename='HI_cl', prefix=params_APS['prefix'], suffix=params_APS['suffix'], path=params_APS['pathout'], header="[1] l, [>2] HI-HI Cl [multipoles, redshift/frequency bins]")    
    if params_general['verbose']: print('Saved at {}'.format(pathname))
    if params_APS['save_clplot_field2']:
        if params_general['verbose']: 
            print('Saving plot...')
        clf2_dict = cxft.cls_from_matrix2dict(np.vstack((l,Chi_bins)) )
        pathname = cxft.save_plots(params_=params_APS, pathname_=pathname, clf_dict=clf2_dict, field='hi')
        if params_general['verbose']: 
            print('Saved at {}'.format(pathname))
    if not params_APS['save_cls_together']: 
        del clf2_dict
    else: 
        plot_together_f2 = 1
    del c_hi, Whi, dChi_dz, pathname, Chi_bins   
    print('\n') 


#Until now, the cross is only for CMBWL x HI
#An ideia is save the kernels and then to use to have the cross-APS
if params_APS['generate_cl_cross']:
    if params_general['verbose']: print('Generating CMBWLxHI Ang Power Spectrum...')
    from scipy import integrate
    for i in np.arange(z_survey.size-1):
        zmin= z_survey[i]
        zmax= z_survey[i+1]
        for j, jl in enumerate(l):
            Whi    = lambda zi: cxft.kernel_hi(camb_params=pars, camb_results=results, z=zi, zmin = zmin, zmax=zmax)
            Wkappa = lambda zi: cxft.kernel_cmb_z(camb_params=pars, camb_results=results, z=zi)
            dCross_dz = lambda zi: (results.hubble_parameter(zi)/c)*pk_interp.P(zi,(jl+0.5)/results.comoving_radial_distance(zi))*(Whi(zi)*Wkappa(zi)/results.comoving_radial_distance(zi)**2)
            cross = integrate.quad(dCross_dz, zmin, zmax)[0]
            if not j:Cross=cross
            else: Cross = np.hstack((Cross, cross))   
        if not i: Cross_bins = Cross
        else:     Cross_bins = np.vstack(( Cross_bins, Cross ))
    if params_general['verbose']: print('Generated.\nSaving...')
    import creating_dir
    creating_dir.verification_dir(params_APS['pathout'].split('/')[-1], "/".join(( params_APS['pathout'].split('/')[:-1] )))
    #creating_dir.verification_dir("auxiliary_data",params_APS['pathout'])
    pathname = cxft.savedata(l_=l, Cl_=np.asmatrix(Cross_bins), filename='CMBWLxHI_cl', prefix=params_APS['prefix'], suffix=params_APS['suffix'], path=params_APS['pathout'], header="[1] l, [>2] CMBWL-HI Cl [multipoles, redshift/frequency bins]")    
    if params_general['verbose']: print('Saved at {}'.format(pathname))
    if params_APS['save_clplot_cross']:
        if params_general['verbose']: 
            print('Saving plot...')
        clcx_dict = cxft.cls_from_matrix2dict(np.vstack((l,Cross_bins)) )
        pathname  = cxft.save_plots(params_=params_APS, pathname_=pathname, clf_dict=clcx_dict, field='cross')
        if params_general['verbose']: 
            print('Saved at {}'.format(pathname))
    if not params_APS['save_cls_together']: 
        del clcx_dict
    else: 
        plot_together_fx = 1        
    del zmin, zmax, dCross_dz, cross, Cross_bins,Cross, Wkappa, Whi  
    print('\n')

if params_APS['save_cls_together'] and plot_together_f1 and plot_together_f2 and plot_together_fx:
    import matplotlib 
    import matplotlib.pyplot as plt
    from matplotlib import cm
    ####################################
    font = {'weight' : 'bold','size'   : 22}
    matplotlib.rc('font', **font)
    plt.rc('font',   size=209)  #set defaults so that the plots are readable
    plt.rc('axes',   titlesize=20)
    plt.rc('axes',   labelsize=20)
    plt.rc('xtick',  labelsize=20)
    plt.rc('ytick',  labelsize=20)
    plt.rc('legend', fontsize =20)
    plt.rc('figure', titlesize=20)
    plt.rc('text',   usetex=True)    
    fig  = plt.figure()
    grid = plt.GridSpec(1,1,top=1.,right=1.5,wspace=0.25)
    ax   = plt.subplot(grid[0,0])	
    fig      = plt.figure()
    grid     = plt.GridSpec(1,1,top=1.,right=1.5,wspace=0.25)
    ax = plt.subplot(grid[0,0])
    A = 1#4/l**4#(l*(l+1))**2/(2*np.pi)
    plt.loglog(clcx_dict['l'],A*clcx_dict['bin1'],linestyle='solid', linewidth=3, color='black', label='bin 1')
    for i,iname in enumerate(clcx_dict.keys()):
        if iname!='l':plt.loglog(clcx_dict['l'],A*clcx_dict[iname],linestyle='solid', color=cm.Purples(i/len(clcx_dict.keys())))
    plt.loglog(clf1_dict['l'], A*clf1_dict['bin1'],linestyle='solid', linewidth=3, color='black', label='bin 1')
    for i,iname in enumerate(clf1_dict.keys()):
        if iname!='l':plt.loglog(clf1_dict['l'], A*clf1_dict[iname],linestyle='solid', color=cm.Blues(i/len(clf1_dict.keys())))
    plt.loglog(clf2_dict['l'], A*clf2_dict['bin1'],linestyle='solid', linewidth=3, color='navy', label='bin 1')
    for i,iname in enumerate(clf2_dict.keys()):
        if iname!='l':plt.loglog(clf2_dict['l'], A*clf2_dict[iname],linestyle='dashed', color=cm.Reds(i/len(clf2_dict.keys())))
    plot_type = ['log','log']
    plt.ylabel(r'$C_{\ell}$')
    plt.xlabel(r'$\ell$')
    if int(params_APS['l_min']): plt.xlim([params_APS['l_min'],params_APS['l_max']])    
    else: plt.xlim([1,params_APS['l_max']]) 
    _lstr_ = np.asarray(pathname.split('/'))
    pathname = '/'.join(( '/'.join( (_lstr_[:-2]) ), 'images', _lstr_[-1].replace('.txt', '.png') ))   
    plt.savefig(pathname, dpi=100, bbox_inches='tight')	
	
if params_general['verbose']:
    print('------------------------------------------')
    print('END')
    print('------------------------------------------\n')

