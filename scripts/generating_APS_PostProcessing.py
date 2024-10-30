#sys.path.insert(1, os.getcwd())
import os, sys, time
from copy import deepcopy as dcopy
import cross_functions_theory as cxft
import cross_functions_simulations as cxfs
import handling_data as hdata
import numpy as np
import healpy as hp
import pandas as pd
import json
import argparse
import warnings
warnings.filterwarnings("ignore")
#####################################################################################################
if 1:
    print('\n===================================================================')
    print('===================================================================')
    print('Starting the code...')
    print('===================================================================')    
    print('===================================================================\n')
#####################################################################################################
timei       = time.time()
#####################################################################################################
# Check the python version and import configparser
#####################################################################################################
if sys.version_info[0]==2:
	import ConfigParser
	config = ConfigParser.RawConfigParser()
elif sys.version_info[0]==3:
	import configparser
	config = configparser.ConfigParser()

#PATH2SCRIPTfiles
#PATH=os.getcwd()
PATH='/data/AMARINS/CMBWLxHI-CODES/scripts'
###################################################################
# This part is for extracting information from parameters.ini file
###################################################################
timei       = time.time()
INI         = "generating_APS_PostProcessing.ini"
name_params = os.path.join(PATH,INI)
config.read(name_params)
#######
#General
verbose = config.getboolean("General","verbose")
project = config.get(       "General","project")
prefix  = config.get(       "General","prefix")
#######
#dataset
nrealizations = config.getint(    "dataset","nrealizations")
lmax          = config.getint(    "dataset","lmax")
mask_used     = config.getboolean("dataset","mask_used" )
ns            = config.getint(    "dataset","ns" )

#######
#binning
#use_namaster = config.getboolean("binning","use_namaster" )
#binning_l    = config.getboolean("binning","binning_l" )
#del_l        = config.getboolean("binning","del_l" )

#######
#PATH
#theory
filepath_field1 = config.get("PATH","filepath_field1")
filepath_field2 = config.get("PATH","filepath_field2")
filepath_cross  = config.get("PATH","filepath_cross" )
#sims
dirpath_sims      = config.get("PATH","dirpath_sims" )
#estimated
dirpath_estimated = config.get("PATH","dirpath_estimated" )
#foreground
dirpath_foregrounds  = config.get("PATH","dirpath_foregrounds" )
filename_foregrounds = config.get("PATH","filename_foregrounds" )
#mask
dirpath_mask         = config.get("PATH","dirpath_mask" )
filename_mask        = config.get("PATH","filename_mask" )
#output_dir
dirpath_output       = config.get("PATH","dirpath_output" )

#######
# FLAGs 
#theory
cl_field1_th = config.getboolean("theory","cl_field1_th" )
cl_field2_th = config.getboolean("theory","cl_field2_th" )
#simulation
cl_field1_sim = config.getboolean("simulated","cl_field1_sim" )
cl_field2_sim = config.getboolean("simulated","cl_field2_sim" )
cl_fg_sim     = config.getboolean("simulated","cl_fg_sim" )
#estimated
cl_field1_rec = config.getboolean("estimated","cl_field1_rec" )
cl_fg_rec     = config.getboolean("estimated","cl_fg_rec" )
filter_fg     = config.getboolean("estimated","filter_fg" )
#leakage
cl_field1_lkg = config.getboolean("leakage","cl_field1_lkg" )
cl_fg_lkg     = config.getboolean("leakage","cl_fg_lkg" )
#cross
cl_cx_field1_fg_sim_sim = config.getboolean("cross","cl_cx_field1_fg_sim_sim" )
cl_cx_field1_fg_rec_rec = config.getboolean("cross","cl_cx_field1_fg_rec_rec" )
cl_cx_field1_fg_lkg_sim = config.getboolean("cross","cl_cx_field1_fg_lkg_sim" )
cl_cx_field1_fg_sim_lkg = config.getboolean("cross","cl_cx_field1_fg_sim_lkg" )
cl_cx_field1_fg_lkg_lkg = config.getboolean("cross","cl_cx_field1_fg_lkg_lkg" )

cl_cx_field2_fg_sim_sim = config.getboolean("cross","cl_cx_field2_fg_sim_sim" )
cl_cx_field2_fg_sim_lkg = config.getboolean("cross","cl_cx_field2_fg_sim_lkg" )

cl_cx_field1_field2_sim_sim = config.getboolean("cross","cl_cx_field1_field2_sim_sim")
cl_cx_field1_field2_rec_sim = config.getboolean("cross","cl_cx_field1_field2_rec_sim")
cl_cx_field1_field2_lkg_sim = config.getboolean("cross","cl_cx_field1_field2_lkg_sim")
#suppression factors
alpha_field1 = config.getboolean("cross","alpha_field1")
alpha_fg     = config.getboolean("cross","alpha_fg")
alpha_cx     = config.getboolean("cross","alpha_cx")

###############################################################################
# You can modify any options in the parameters.ini file by the command terminal
###############################################################################
parser = argparse.ArgumentParser(description='Modify by the command terminal parameters in {} file'.format(INI))
######
#General
parser.add_argument('--verbose', action = 'store', dest = 'verbose', default = verbose, help = '')
parser.add_argument('--project', action = 'store', dest = 'project', default = project, help = '')
parser.add_argument('--prefix' , action = 'store', dest = 'prefix' , default = prefix , help = '')
######
#dataset
parser.add_argument('--nrealizations', action = 'store', dest = 'nrealizations', default = nrealizations, help = '')
parser.add_argument('--lmax'         , action = 'store', dest = 'lmax'         , default = lmax         , help = '')
parser.add_argument('--mask_used'    , action = 'store', dest = 'mask_used'    , default = mask_used    , help = '')
parser.add_argument('--ns'           , action = 'store', dest = 'ns'           , default = ns           , help = '')
######
#PATH
#theory
parser.add_argument('--filepath_field1'     , action = 'store', dest = 'filepath_field1'     , default = filepath_field1     , help = '')
parser.add_argument('--filepath_field2'     , action = 'store', dest = 'filepath_field2'     , default = filepath_field2     , help = '')
parser.add_argument('--filepath_cross'      , action = 'store', dest = 'filepath_cross'      , default = filepath_cross      , help = '')
#sims
parser.add_argument('--dirpath_sims'        , action = 'store', dest = 'dirpath_sims'        , default = dirpath_sims        , help = '')
#estimated
parser.add_argument('--dirpath_estimated'   , action = 'store', dest = 'dirpath_estimated'   , default = dirpath_estimated   , help = '')
#foregrounds
parser.add_argument('--dirpath_foregrounds' , action = 'store', dest = 'dirpath_foregrounds' , default = dirpath_foregrounds , help = '')
parser.add_argument('--filename_foregrounds', action = 'store', dest = 'filename_foregrounds', default = filename_foregrounds, help = '')
#mask
parser.add_argument('--dirpath_mask'        , action = 'store', dest = 'dirpath_mask'        , default = dirpath_mask        , help = '')
parser.add_argument('--filename_mask'       , action = 'store', dest = 'filename_mask'       , default = filename_mask       , help = '')
#output_dir
parser.add_argument('--dirpath_output'      , action = 'store', dest = 'dirpath_output'      , default = dirpath_output      , help = '')

######
#FLAGS
#theory
parser.add_argument('--cl_field1_th' , action = 'store', dest = 'cl_field1_th' , default = cl_field1_th , help = '')
parser.add_argument('--cl_field2_th' , action = 'store', dest = 'cl_field2_th' , default = cl_field2_th , help = '')
#sims
parser.add_argument('--cl_field1_sim', action = 'store', dest = 'cl_field1_sim', default = cl_field1_sim, help = '')
parser.add_argument('--cl_field2_sim', action = 'store', dest = 'cl_field2_sim', default = cl_field2_sim, help = '')
parser.add_argument('--cl_fg_sim'    , action = 'store', dest = 'cl_fg_sim'    , default = cl_fg_sim    , help = '')
#estimated
parser.add_argument('--cl_field1_rec', action = 'store', dest = 'cl_field1_rec', default = cl_field1_rec, help = '')
parser.add_argument('--cl_fg_rec', action = 'store', dest = 'cl_fg_rec', default = cl_fg_rec, help = '')
parser.add_argument('--filter_fg', action = 'store', dest = 'filter_fg', default = filter_fg, help = '')
#leakage
parser.add_argument('--cl_field1_lkg', action = 'store', dest = 'cl_field1_lkg', default = cl_field1_lkg, help = '')
parser.add_argument('--cl_fg_lkg', action = 'store', dest = 'cl_fg_lkg', default = cl_fg_lkg, help = '')
#cross
parser.add_argument('--cl_cx_field1_fg_sim_sim', action = 'store', dest = 'cl_cx_field1_fg_sim_sim', default = cl_cx_field1_fg_sim_sim, help = '')
parser.add_argument('--cl_cx_field1_fg_rec_rec', action = 'store', dest = 'cl_cx_field1_fg_rec_rec', default = cl_cx_field1_fg_rec_rec, help = '')
parser.add_argument('--cl_cx_field1_fg_lkg_sim', action = 'store', dest = 'cl_cx_field1_fg_lkg_sim', default = cl_cx_field1_fg_lkg_sim, help = '')
parser.add_argument('--cl_cx_field1_fg_sim_lkg', action = 'store', dest = 'cl_cx_field1_fg_sim_lkg', default = cl_cx_field1_fg_sim_lkg, help = '')
parser.add_argument('--cl_cx_field1_fg_lkg_lkg', action = 'store', dest = 'cl_cx_field1_fg_lkg_lkg', default = cl_cx_field1_fg_lkg_lkg, help = '')
parser.add_argument('--cl_cx_field2_fg_sim_sim', action = 'store', dest = 'cl_cx_field2_fg_sim_sim', default = cl_cx_field2_fg_sim_sim, help = '')
parser.add_argument('--cl_cx_field2_fg_sim_lkg', action = 'store', dest = 'cl_cx_field2_fg_sim_lkg', default = cl_cx_field2_fg_sim_lkg, help = '')
parser.add_argument('--cl_cx_field1_field2_sim_sim', action = 'store', dest = 'cl_cx_field1_field2_sim_sim', default = cl_cx_field1_field2_sim_sim, help = '')
parser.add_argument('--cl_cx_field1_field2_rec_sim', action = 'store', dest = 'cl_cx_field1_field2_rec_sim', default = cl_cx_field1_field2_rec_sim, help = '')
parser.add_argument('--cl_cx_field1_field2_lkg_sim', action = 'store', dest = 'cl_cx_field1_field2_lkg_sim', default = cl_cx_field1_field2_lkg_sim, help = '')
#suppression factors 
parser.add_argument('--alpha_field1', action = 'store', dest = 'alpha_field1', default = alpha_field1, help = '')
parser.add_argument('--alpha_fg'    , action = 'store', dest = 'alpha_fg'    , default = alpha_fg, help = '')
parser.add_argument('--alpha_cx'    , action = 'store', dest = 'alpha_cx'    , default = alpha_cx, help = '')

arguments = parser.parse_args()
###############################################################################
# Variables
###############################################################################
#General
verbose = bool(arguments.verbose)
project = str( arguments.project)
prefix  = str( arguments.prefix)
#########
#dataset
nrealizations = int( arguments.nrealizations)
lmax          = int( arguments.lmax)
mask_used     = bool(arguments.mask_used)
ns            = int( arguments.ns)
#########
#PATH
#theory
filepath_field1 = str(arguments.filepath_field1)
filepath_field2 = str(arguments.filepath_field2)
filepath_cross  = str(arguments.filepath_cross)
#sims
dirpath_sims = str(arguments.dirpath_sims)
#estimated
dirpath_estimated = str(arguments.dirpath_estimated)
#foregrounds
dirpath_foregrounds  = str(arguments.dirpath_foregrounds)
filename_foregrounds = str(arguments.filename_foregrounds)
#mask
dirpath_mask   = str(arguments.dirpath_mask)
filename_mask  = str(arguments.filename_mask)
#output_dir
dirpath_output = str(arguments.dirpath_output)

######
#FLAGS
#theory
cl_field1_th = bool(arguments.cl_field1_th)
cl_field2_th = bool(arguments.cl_field2_th)
#sims
cl_field1_sim = bool(arguments.cl_field1_sim)
cl_field2_sim = bool(arguments.cl_field2_sim)
cl_fg_sim     = bool(arguments.cl_fg_sim)
#estimated
cl_field1_rec = bool(arguments.cl_field1_rec)
cl_fg_rec     = bool(arguments.cl_fg_rec)
filter_fg     = bool(arguments.filter_fg)
#leakage
cl_field1_lkg = bool(arguments.cl_field1_lkg)
cl_fg_lkg     = bool(arguments.cl_fg_lkg)
#cross
cl_cx_field1_fg_sim_sim = bool(arguments.cl_cx_field1_fg_sim_sim)
cl_cx_field1_fg_rec_rec = bool(arguments.cl_cx_field1_fg_rec_rec)
cl_cx_field1_fg_lkg_sim = bool(arguments.cl_cx_field1_fg_lkg_sim)
cl_cx_field1_fg_sim_lkg = bool(arguments.cl_cx_field1_fg_sim_lkg)
cl_cx_field1_fg_lkg_lkg = bool(arguments.cl_cx_field1_fg_lkg_lkg)

cl_cx_field2_fg_sim_sim = bool(arguments.cl_cx_field2_fg_sim_sim)
cl_cx_field2_fg_sim_lkg = bool(arguments.cl_cx_field2_fg_sim_lkg)

cl_cx_field1_field2_sim_sim = bool(arguments.cl_cx_field1_field2_sim_sim)
cl_cx_field1_field2_rec_sim = bool(arguments.cl_cx_field1_field2_rec_sim)
cl_cx_field1_field2_lkg_sim = bool(arguments.cl_cx_field1_field2_lkg_sim)



dirpath_estimated = os.path.join(dirpath_estimated, "ns{}".format(ns))
#cxfs.verification_dir(dirname="", path=dirpath_output, clear=False, verbose=False)
dirpath_output    = os.path.join(dirpath_output, "ns{}".format(ns))
#####################################################################################################################################################################
#####################################################################################################################################################################
params_general = {'project':project, 'verbose':verbose,'prefix':prefix, 'nreals':nrealizations, 'lmax': lmax, 'mask_used': mask_used, 
                  'ns':ns}#, 'use_namaster':use_namaster, 'binning_l':binning_l, 'del_l':del_l}

params_path = {'dirpath_field1':filepath_field1,     'dirpath_field2':filepath_field2, 'dirpath_cross':filepath_cross,
               'dirpath_sims'  :dirpath_sims,        'dirpath_rec':dirpath_estimated, 
               'dirpath_fg'    :dirpath_foregrounds, 'filename_fg':filename_foregrounds, 
               'dirpath_mask'  :dirpath_mask,        'filename_mask':filename_mask,
               'dirpath_output':dirpath_output}

params_flags = {'cl_field1_th' :cl_field1_th,  'cl_field2_th' :cl_field2_th, 
                'cl_field1_sim':cl_field1_sim, 'cl_field2_sim':cl_field2_sim, 'cl_fg_sim': cl_fg_sim,
                'cl_field1_rec':cl_field1_rec, 'cl_fg_rec': cl_fg_rec,        'filter_fg':filter_fg,
                'cl_field1_lkg':cl_field1_lkg, 'cl_fg_lkg': cl_fg_lkg,
                'cl_cx_field1_fg_sim_sim':cl_cx_field1_fg_sim_sim, 'cl_cx_field1_fg_rec_rec':cl_cx_field1_fg_rec_rec,
                'cl_cx_field1_fg_lkg_sim':cl_cx_field1_fg_lkg_sim, 'cl_cx_field1_fg_sim_lkg':cl_cx_field1_fg_sim_lkg,
                'cl_cx_field1_fg_lkg_lkg':cl_cx_field1_fg_lkg_lkg,
                'cl_cx_field2_fg_sim_sim':cl_cx_field2_fg_sim_sim,         'cl_cx_field2_fg_sim_lkg':cl_cx_field2_fg_sim_lkg,
                'cl_cx_field1_field2_sim_sim':cl_cx_field1_field2_sim_sim, 'cl_cx_field1_field2_rec_sim':cl_cx_field1_field2_rec_sim,
                'cl_cx_field1_field2_lkg_sim':cl_cx_field1_field2_lkg_sim,
                'alpha_field1':alpha_field1, 'alpha_fg':alpha_fg, 'alpha_cx':alpha_cx
               }
if params_general['prefix']=="":params_general['prefix']=None
del project, prefix, nrealizations, lmax, mask_used, ns#, use_namaster, binning_l, del_l
del filepath_field1, filepath_field2, filepath_cross, dirpath_sims, dirpath_estimated
del dirpath_foregrounds, dirpath_mask, dirpath_output, filename_foregrounds, filename_mask
del cl_field1_th, cl_field2_th
del cl_field1_sim, cl_field2_sim, cl_fg_sim
del cl_field1_rec, cl_fg_rec, filter_fg
del cl_field1_lkg, cl_fg_lkg
del cl_cx_field1_fg_sim_sim,     cl_cx_field1_fg_rec_rec,     cl_cx_field1_fg_lkg_sim
del cl_cx_field1_fg_sim_lkg,     cl_cx_field1_fg_lkg_lkg
del cl_cx_field2_fg_sim_sim,     cl_cx_field2_fg_sim_lkg
del cl_cx_field1_field2_sim_sim, cl_cx_field1_field2_rec_sim, cl_cx_field1_field2_lkg_sim
del alpha_field1, alpha_fg, alpha_cx

#if params_general['verbose']:
#    print('\n--Variables to be used:--')
#    print(params_general)    
#    print(params_flags)    
#    print(params_path)
#    print()
if params_general['verbose']:
    #print('\n------------------------------------------------------------------')
    print('Verifying the output directory...')
    #print('------------------------------------------------------------------\n')
cxfs.verification_dir(dirname="", path=params_path['dirpath_output'], clear=False, verbose=False)
_clear_=False

if params_general['verbose']:
    #print('\n------------------------------------------------------------------')
    print('Loading the APS theoretically simulated...')
    #print('------------------------------------------------------------------\n')
clf1 = np.loadtxt(params_path['dirpath_field1']).T[1:,:]
clf2 = np.loadtxt(params_path['dirpath_field2']).T[1:,:]
clcx = np.loadtxt(params_path['dirpath_cross' ]).T[1:,:] 


for jflag,clf in zip(['cl_field1_th', 'cl_field2_th'], [clf1,clf2]):
        params_path['dirpath_output_theory'] = cxfs.verification_dir(dirname='theory', path=params_path['dirpath_output'], clear=False, verbose=False)
        cxfs.savedata_from_dict(Cl_=np.vstack(( np.arange(clf.shape[1]), clf )),  
                                filename=jflag, path=params_path['dirpath_output_theory'], 
                                prefix=params_general['prefix'], clear=False)
        del jflag, clf


if params_general['verbose']:
    #print('\n------------------------------------------------------------------')
    print('Loading the Maps of Foregrounds...')
    #print('------------------------------------------------------------------\n')
MFG  = hdata.getmap(dirpath_=params_path['dirpath_fg'], filename_=params_path['filename_fg'], healpix_readingformat=0, hdu=1)
if params_general['mask_used']:
    MASK = hdata.getmap(dirpath_=params_path['dirpath_mask'], filename_=params_path['filename_mask'], healpix_readingformat=0, hdu=1)
else:
    MASK = np.ones_like(MFG[0])
MFG = MFG*MASK
params_vars = {'map_fg_sim':MFG}
del MASK,MFG

#### storading info
params_general['nside']=hp.get_nside(params_vars['map_fg_sim'][1])
params_general['npix' ]=hp.get_map_size(params_vars['map_fg_sim'][1])
params_general['nch'  ]=int(clf1.shape[0])
del clf1, clf2  #jflag.replace('_field_','hi')
#### loading the simulation names to be used
if params_general['verbose']:
    #print('\n------------------------------------------------------------------')
    print('Loading the simulation names to be used...')
    #print('------------------------------------------------------------------\n')
sim_names = np.array([])
for jname in os.listdir(params_path['dirpath_rec' ]):
    if 'sim' in jname:
        jnum = int(jname.split('sim')[1])
        sim_names = np.hstack(( sim_names , jnum ))
sim_names = np.array([ 'sim{}'.format(jsim) for jsim in np.asarray(np.sort(sim_names), dtype=int) ])
params_path['sim_names'] = sim_names[:params_general['nreals']]
del sim_names, jnum,jname

###############################
### PART1
if params_flags['cl_field1_lkg'] or params_flags['cl_cx_field1_fg_sim_sim'] or params_flags['cl_cx_field1_fg_sim_lkg'] \
   or params_flags['cl_cx_field1_field2_sim_sim']:
    timej  = time.time()
    if params_general['verbose']:
        print('\n------------------------------------------------------------------')
        print('Loading alm coeficients of Field1 simulated...')    
    for j,jsim in enumerate(params_path['sim_names']):
        alm_hi_sim_imag = np.loadtxt(os.path.join(params_path['dirpath_sims' ],jsim, 'alm_hi_sim_imag.txt'))[:,2:].T
        alm_hi_sim_real = np.loadtxt(os.path.join(params_path['dirpath_sims' ],jsim, 'alm_hi_sim_real.txt'))[:,2:].T    
        for jj in range(params_general['nch']):
            if not jj:
                alm_hi_sim = np.array([complex(jreal,jcompl) for (jreal,jcompl) in zip(alm_hi_sim_imag[jj], alm_hi_sim_imag[jj])])
                alm_hi_sim = np.ascontiguousarray(alm_hi_sim)
            else:
                alm_hi_sim = np.vstack(( alm_hi_sim, np.ascontiguousarray(np.array([complex(jreal,jcompl) for (jreal,jcompl) in zip(alm_hi_sim_imag[jj], alm_hi_sim_imag[jj])])) ))
        del alm_hi_sim_imag, alm_hi_sim_real
        alm_hi_sim_nreals = alm_hi_sim.flatten() if not j else np.vstack(( alm_hi_sim_nreals, alm_hi_sim.flatten() ))
        params_vars['alm_field1_sim_matrix'] = alm_hi_sim_nreals
    del alm_hi_sim,jsim,alm_hi_sim_nreals
    if params_general['verbose']: print('Loaded {0:d} sims of field1: {1:.4f} seg'.format(params_general['nreals'], time.time()-timej)) 
    del timej


if params_flags['cl_cx_field2_fg_sim_sim'] or params_flags['cl_cx_field2_fg_sim_lkg'] or params_flags['cl_cx_field1_field2_rec_sim'] \
   or params_flags['cl_cx_field1_field2_lkg_sim']:
    timej  = time.time()
    if params_general['verbose']:
        print('\n------------------------------------------------------------------')
        print('Loading alm coeficients of Field2 simulated...')       
    for j,jsim in enumerate(params_path['sim_names']):
        alm_kp_sim_imag = np.loadtxt(os.path.join(params_path['dirpath_sims' ],jsim, 'alm_kappa_sim_imag.txt')).T[2]
        alm_kp_sim_real = np.loadtxt(os.path.join(params_path['dirpath_sims' ],jsim, 'alm_kappa_sim_real.txt')).T[2]
        alm_kp_sim = np.array([complex(jreal,jcompl) for (jreal,jcompl) in zip(alm_kp_sim_real, alm_kp_sim_imag)])
        alm_kp_sim = np.ascontiguousarray(alm_kp_sim)
        del alm_kp_sim_imag, alm_kp_sim_real
        alm_kp_sim_nreals = alm_kp_sim if not j else np.vstack(( alm_kp_sim_nreals, alm_kp_sim ))
        params_vars['alm_field2_sim_matrix'] = alm_kp_sim_nreals
    del alm_kp_sim,jsim
    if params_general['verbose']: print('Loaded {0:d} sims of field2: {1:.4f} seg'.format(params_general['nreals'], time.time()-timej)) 
    del timej

if params_flags['cl_cx_field1_fg_sim_sim'] or params_flags['cl_cx_field1_fg_lkg_sim'] or \
   params_flags['cl_cx_field2_fg_sim_sim']:
    timej  = time.time()
    if params_general['verbose']:
        print('\n------------------------------------------------------------------')
        print('Loading alm coeficients of FG simulated...')      
    alm_fg_sim  = np.vstack([ hp.map2alm(params_vars['map_fg_sim'][jch],pol=False) for jch in range(params_general['nch']) ])
    params_vars['alm_fg_sim_matrix'] = alm_fg_sim.flatten()
    del alm_fg_sim
    if params_general['verbose']: print('Loaded FG sim: {0:.4f} seg'.format(time.time()-timej)) 
    del timej

###############################
### PART2
if params_flags['cl_field1_sim']:
    timej  = time.time()
    if params_general['verbose']: 
        print('\n------------------------------------------------------------------')
        print('Saving {} ...'.format('cl_field1_sim')) 
    for j,jsim in enumerate(params_path['sim_names']):
        cl_hi_jsim = np.loadtxt(os.path.join(params_path['dirpath_sims' ],jsim, 'cl_hi_sim.txt')).T #[l, CL]
        params_path['dirpath_output_simulations'] = cxfs.verification_dir(dirname='simulations', path=params_path['dirpath_output'], clear=False, verbose=False)
        params_path['dirpath_output_jsim']        = cxfs.verification_dir(dirname=jsim, path=params_path['dirpath_output_simulations'], clear=False, verbose=False)
        cxfs.savedata_from_dict(Cl_=cl_hi_jsim, filename='cl_field1_sim', path=params_path['dirpath_output_jsim'], 
                                prefix=params_general['prefix'], clear=False)
        params_vars['cl_field1_sim_matrix'] = cl_hi_jsim[1:,:].flatten() if not j else np.vstack(( params_vars['cl_field1_sim_matrix'], cl_hi_jsim[1:,:].flatten() ))            
        if params_general['verbose']: print('Saved ({}) at: {}'.format(jsim, params_path['dirpath_output_jsim'])) 
        del j,jsim,cl_hi_jsim,params_path['dirpath_output_jsim']
    params_vars['cl_field1_sim_mean']  = np.average(params_vars['cl_field1_sim_matrix'],axis=0).reshape(params_general['nch'],-1)
    params_vars['cl_field1_sim_mean']  = np.vstack(( np.arange(params_vars['cl_field1_sim_mean'].shape[1]), params_vars['cl_field1_sim_mean'] ))
    params_path['dirpath_output_mean'] = cxfs.verification_dir(dirname='mean', path=params_path['dirpath_output_simulations'], clear=False, verbose=False)
    cxfs.savedata_from_dict(Cl_=params_vars['cl_field1_sim_mean'], filename='cl_field1_sim_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False)   
    if params_general['verbose']: print('Saved (mean) at: {}'.format(params_path['dirpath_output_mean'])) 
    if not params_flags['alpha_field1']:   del params_vars['cl_field1_sim_matrix']    
    if params_general['verbose']: print('Processing time: {0:.4f} seg'.format(time.time()-timej)) 
    del params_path['dirpath_output_mean'],params_vars['cl_field1_sim_mean'], timej


if params_flags['cl_field2_sim']:
    timej  = time.time()
    if params_general['verbose']: 
        print('\n------------------------------------------------------------------')
        print('Saving {} ...'.format('cl_field2_sim'))     
    for j,jsim in enumerate(params_path['sim_names']):
        cl_kp_jsim = np.loadtxt(os.path.join(params_path['dirpath_sims' ],jsim, 'cl_kappa_sim.txt')).T #[l, CL]
        params_path['dirpath_output_simulations'] = cxfs.verification_dir(dirname='simulations', path=params_path['dirpath_output'], clear=False, verbose=False)
        params_path['dirpath_output_jsim']        = cxfs.verification_dir(dirname=jsim, path=params_path['dirpath_output_simulations'], clear=False, verbose=False)
        cxfs.savedata_from_dict(Cl_=cl_kp_jsim, filename='cl_field2_sim', path=params_path['dirpath_output_jsim'], 
                                 prefix=params_general['prefix'], clear=False)
        params_vars['cl_field2_sim_matrix'] = cl_kp_jsim[1:,:] if not j else np.vstack(( params_vars['cl_field2_sim_matrix'], cl_kp_jsim[1:,:] ))            
        if params_general['verbose']: print('Saved ({}) at: {}'.format(jsim, params_path['dirpath_output_jsim'])) 
        del j,jsim,cl_kp_jsim,params_path['dirpath_output_jsim']            
    params_vars['cl_field2_sim_mean']  = np.average(params_vars['cl_field2_sim_matrix'],axis=0)
    params_vars['cl_field2_sim_mean']  = np.vstack(( np.arange(params_vars['cl_field2_sim_mean'].size), params_vars['cl_field2_sim_mean'] ))
    params_path['dirpath_output_mean'] = cxfs.verification_dir(dirname='mean', path=params_path['dirpath_output_simulations'], clear=False, verbose=False)    
    cxfs.savedata_from_dict(Cl_=params_vars['cl_field2_sim_mean'], filename='cl_field2_sim_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False)   
    if params_general['verbose']: print('Saved (mean) at: {}'.format(params_path['dirpath_output_mean'])) 
    if params_general['verbose']: print('Processing time: {0:.4f} seg'.format(time.time()-timej))         
    del params_path['dirpath_output_mean'], params_vars['cl_field2_sim_matrix'], params_vars['cl_field2_sim_mean'], timej

    
if params_flags['cl_cx_field1_field2_sim_sim']:
    timej  = time.time()    
    if params_general['verbose']: 
        print('\n------------------------------------------------------------------')
        print('Saving {} ...'.format('cl_cx_field1_field2_sim_sim'))       
    for j,jsim in enumerate(params_path['sim_names']):
        cl_cx_jsim = np.loadtxt(os.path.join(params_path['dirpath_sims' ],jsim, 'cl_cross_sim.txt')).T #[l, CL]
        params_path['dirpath_output_simulations'] = cxfs.verification_dir(dirname='simulations', path=params_path['dirpath_output'], clear=False, verbose=False)
        params_path['dirpath_output_jsim']        = cxfs.verification_dir(dirname=jsim, path=params_path['dirpath_output_simulations'], clear=False, verbose=False)
        cxfs.savedata_from_dict(Cl_=cl_cx_jsim, filename='cl_cx_field1_field2_sim_sim', path=params_path['dirpath_output_jsim'], 
                                prefix=params_general['prefix'], clear=False)
        params_vars['cl_cx_field1_field2_sim_sim_matrix'] = cl_cx_jsim[1:,:].flatten() if not j else np.vstack(( params_vars['cl_cx_field1_field2_sim_sim_matrix'], cl_cx_jsim[1:,:].flatten() ))            
        if params_general['verbose']: print('Saved ({}) at: {}'.format(jsim, params_path['dirpath_output_jsim']))         
        del j,jsim,cl_cx_jsim,params_path['dirpath_output_jsim']            
    params_vars['cl_cx_field1_field2_sim_sim_mean'] = np.average(params_vars['cl_cx_field1_field2_sim_sim_matrix'],axis=0).reshape(params_general['nch'],-1)
    params_vars['cl_cx_field1_field2_sim_sim_mean'] = np.vstack(( np.arange(params_vars['cl_cx_field1_field2_sim_sim_mean'].shape[1]), params_vars['cl_cx_field1_field2_sim_sim_mean'] ))
    params_path['dirpath_output_mean']              = cxfs.verification_dir(dirname='mean', path=params_path['dirpath_output_simulations'], clear=False, verbose=False)
    cxfs.savedata_from_dict(Cl_=params_vars['cl_cx_field1_field2_sim_sim_mean'], filename='cl_cx_field1_field2_sim_sim_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False)           
    #if params_flags['alpha_cx']: params_vars['cl_cx_sim_matrix'] = cl_cx_sim_matrix
    if params_general['verbose']: print('Saved (mean) at: {}'.format(params_path['dirpath_output_mean']))     
    if not params_flags['alpha_cx']: params_vars['cl_cx_field1_field2_sim_sim_matrix']
    if params_general['verbose']: print('Processing time: {0:.4f} seg'.format(time.time()-timej))         
    del params_path['dirpath_output_mean'],params_vars['cl_cx_field1_field2_sim_sim_mean'], timej


if params_flags['cl_fg_sim']:
    timej  = time.time()    
    if params_general['verbose']: 
        print('\n------------------------------------------------------------------')
        print('Saving {} ...'.format('cl_fg_sim'))      
    params_vars['cl_fg_sim_matrix'] = np.array([ hp.anafast(params_vars['map_fg_sim'][jch],pol=False) for jch in range(params_general['nch']) ]) 
    params_vars['cl_fg_sim_matrix'] = np.vstack(( np.arange(params_vars['cl_fg_sim_matrix'].shape[1]),params_vars['cl_fg_sim_matrix'] ))
    params_path['dirpath_output_simulations'] = cxfs.verification_dir(dirname='simulations', path=params_path['dirpath_output'], clear=False, verbose=False)
    params_path['dirpath_output_fg'] = cxfs.verification_dir(dirname='foregrounds', path=params_path['dirpath_output_simulations'], clear=False, verbose=False)
    cxfs.savedata_from_dict(Cl_=params_vars['cl_fg_sim_matrix'], filename='cl_fg_sim', path=params_path['dirpath_output_fg'], 
                             prefix=params_general['prefix'], clear=False)
    if params_general['verbose']: print('Saved at: {}'.format(params_path['dirpath_output_fg']))         
    if not params_flags['alpha_fg']: params_vars['cl_fg_sim_matrix']
    if params_general['verbose']: print('Processing time: {0:.4f} seg'.format(time.time()-timej))            
    del params_path['dirpath_output_fg'],timej

###############################
### PART3
if params_flags['cl_cx_field1_fg_sim_sim']:
    timej  = time.time()   
    if params_general['verbose']: 
        print('\n------------------------------------------------------------------')
        print('Saving {} ...'.format('cl_cx_field1_fg_sim_sim'))        
    for j,jsim in enumerate(params_path['sim_names']):
        _alm_f1_j_ = params_vars['alm_field1_sim_matrix'][j].reshape(params_general['nch'],-1)
        _alm_fg_   = params_vars['alm_fg_sim_matrix'    ].reshape(params_general['nch'],-1)
        _cx_ = np.array([ hp.alm2cl(alms1=_alm_f1_j_[jch], alms2=_alm_fg_[jch]) for jch in range(params_general['nch']) ]) 
        params_path['dirpath_output_cross'] = cxfs.verification_dir(dirname='cross', path=params_path['dirpath_output'], clear=False, verbose=False)
        params_path['dirpath_output_jcross'] = cxfs.verification_dir(dirname=jsim, path=params_path['dirpath_output_cross'], clear=False, verbose=False)
        cxfs.savedata_from_dict(Cl_=np.vstack(( np.arange(_cx_.shape[1]), _cx_ )), filename='cl_cx_field1_fg_sim_sim', path=params_path['dirpath_output_jcross'], 
                                prefix=params_general['prefix'], clear=False)        
        params_vars['cl_cx_field1_fg_sim_sim_matrix'] = _cx_.flatten() if not j else np.vstack(( params_vars['cl_cx_field1_fg_sim_sim_matrix'],_cx_.flatten() ))
        if params_general['verbose']: print('Saved ({}) at: {}'.format(jsim, params_path['dirpath_output_jcross'])) 
        del _alm_f1_j_, _alm_fg_,jsim,j,_cx_
    params_vars['cl_cx_field1_fg_sim_sim_mean'] = np.average(params_vars['cl_cx_field1_fg_sim_sim_matrix'],axis=0).reshape(params_general['nch'],-1)
    params_vars['cl_cx_field1_fg_sim_sim_mean'] = np.vstack(( np.arange(params_vars['cl_cx_field1_fg_sim_sim_mean'].shape[1]), params_vars['cl_cx_field1_fg_sim_sim_mean'] ))
    params_path['dirpath_output_mean'] = cxfs.verification_dir(dirname='mean', path=params_path['dirpath_output_cross'], clear=False, verbose=False)
    cxfs.savedata_from_dict(Cl_=params_vars['cl_cx_field1_fg_sim_sim_mean'], filename='cl_cx_field1_fg_sim_sim_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False)                   
    if params_general['verbose']: print('Saved (mean) at: {}'.format(params_path['dirpath_output_mean'])) 
    if params_general['verbose']: print('Processing time: {0:.4f} seg'.format(time.time()-timej))       
    del params_vars['cl_cx_field1_fg_sim_sim_matrix'],params_vars['cl_cx_field1_fg_sim_sim_mean']


if params_flags['cl_cx_field2_fg_sim_sim']:
    timej  = time.time()   
    if params_general['verbose']: 
        print('\n------------------------------------------------------------------')
        print('Saving {} ...'.format('cl_cx_field1_fg_sim_sim'))     
    for j,jsim in enumerate(params_path['sim_names']):
        _alm_f2_j_ = params_vars['alm_field2_sim_matrix'][j]
        _alm_fg_   = params_vars['alm_fg_sim_matrix'    ].reshape(params_general['nch'],-1)
        _cx_ = np.array([ hp.alm2cl(alms1=_alm_f2_j_, alms2=_alm_fg_[jch]) for jch in range(params_general['nch']) ]) 
        params_path['dirpath_output_cross'] = cxfs.verification_dir(dirname='cross', path=params_path['dirpath_output'], clear=False, verbose=False)
        params_path['dirpath_output_jcross'] = cxfs.verification_dir(dirname=jsim, path=params_path['dirpath_output_cross'], clear=False, verbose=False)
        cxfs.savedata_from_dict(Cl_=np.vstack(( np.arange(_cx_.shape[1]), _cx_ )), filename='cl_cx_field2_fg_sim_sim', path=params_path['dirpath_output_jcross'], 
                                prefix=params_general['prefix'], clear=False)        
        params_vars['cl_cx_field2_fg_sim_sim_matrix'] = _cx_.flatten() if not j else np.vstack(( params_vars['cl_cx_field2_fg_sim_sim_matrix'],_cx_.flatten() ))
        if params_general['verbose']: print('Saved ({}) at: {}'.format(jsim, params_path['dirpath_output_jcross'])) 
        del _alm_f2_j_, _alm_fg_,jsim,j,_cx_
    params_vars['cl_cx_field2_fg_sim_sim_mean'] = np.average(params_vars['cl_cx_field2_fg_sim_sim_matrix'],axis=0).reshape(params_general['nch'],-1)
    params_vars['cl_cx_field2_fg_sim_sim_mean'] = np.vstack(( np.arange(params_vars['cl_cx_field2_fg_sim_sim_mean'].shape[1]), params_vars['cl_cx_field2_fg_sim_sim_mean'] ))
    params_path['dirpath_output_mean'] = cxfs.verification_dir(dirname='mean', path=params_path['dirpath_output_cross'], clear=False, verbose=False)
    cxfs.savedata_from_dict(Cl_=params_vars['cl_cx_field2_fg_sim_sim_mean'], filename='cl_cx_field2_fg_sim_sim_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False)                   
    if params_general['verbose']: print('Saved (mean) at: {}'.format(params_path['dirpath_output_mean'])) 
    if params_general['verbose']: print('Processing time: {0:.4f} seg'.format(time.time()-timej))           
    del params_vars['cl_cx_field2_fg_sim_sim_matrix'],params_vars['cl_cx_field2_fg_sim_sim_mean']


###############################
### PART4 --- harder part
print('\n------------------------------------------------------------------')
timeii = time.time()
if params_flags['cl_field1_rec'] or params_flags['cl_field1_lkg'] \
    or params_flags['filter_fg'] \
    or params_flags['cl_fg_lkg'] or params_flags['cl_cx_field1_fg_rec_rec'] \
    or params_flags['cl_cx_field1_fg_lkg_sim'] or params_flags['cl_cx_field1_fg_lkg_lkg']\
    or params_flags['cl_cx_field1_field2_rec_sim'] or params_flags['cl_cx_field1_field2_lkg_sim']:
    
    for j,jsim in enumerate(params_path['sim_names']):   
        time_vec = np.array([0])
        params_path['dirpath_jrec'] = os.path.join(params_path['dirpath_rec'],jsim)
        for jname in os.listdir(params_path['dirpath_jrec']):
            if 'extHI'     in jname: jhifilename = jname
            if 'extFG'     in jname: jfgfilename = jname
            if 'mixmatrix' in jname: jmixname    = jname
        _mp_ = hdata.getmap(dirpath_=params_path['dirpath_jrec'], filename_=jhifilename, healpix_readingformat=0, hdu=1)
        if params_flags['cl_field1_rec'] or params_flags['cl_cx_field1_fg_rec_rec'] \
           or params_flags['cl_cx_field1_field2_rec_sim']:
            _alm_hi_rec_ = hp.map2alm(maps=_mp_,pol=False)
            if params_flags['cl_field1_rec']:  
                timej  = time.time()
                _cx_  = np.vstack([hp.alm2cl(_alm_hi_rec_[jch]) for jch in range(params_general['nch'])])
                params_path['dirpath_output_rec']  = cxfs.verification_dir(dirname='estimated', path=params_path['dirpath_output'], clear=False, verbose=False)
                params_path['dirpath_output_jrec'] = cxfs.verification_dir(dirname=jsim   , path=params_path['dirpath_output_rec'], clear=False, verbose=False)                    
                cxfs.savedata_from_dict(Cl_=np.vstack(( np.arange(_cx_.shape[1]), _cx_ )), filename='cl_field1_rec', path=params_path['dirpath_output_jrec'], 
                                            prefix=params_general['prefix'], clear=False)  
                if params_general['verbose']: print('Saved ({}-{}) at              : {}'.format(jsim, 'cl_field1_rec', params_path['dirpath_output_jrec'])) #timer?
                if not j: 
                    params_vars['cl_field1_rec_matrix'] = dcopy(_cx_.flatten())
                else:
                    params_vars['cl_field1_rec_matrix'] = np.vstack(( params_vars['cl_field1_rec_matrix'], dcopy(_cx_.flatten()) ))     
                time_vec = np.hstack((time_vec,time.time()-timej))
                del _cx_,timej
            if params_flags['cl_cx_field1_fg_rec_rec']:  
                timej  = time.time()
                _mfg_ = hdata.getmap(dirpath_=params_path['dirpath_jrec'], filename_=jfgfilename, healpix_readingformat=0, hdu=1)
                _alm_fg_rec_ = hp.map2alm(maps=_mfg_,pol=False)
                _cx_  = np.vstack([hp.alm2cl(alms1=_alm_hi_rec_[jch],alms2=_alm_fg_rec_[jch]) for jch in range(params_general['nch'])])
                params_path['dirpath_output_cross']  = cxfs.verification_dir(dirname='cross',    path=params_path['dirpath_output'], clear=False, verbose=False)
                params_path['dirpath_output_jcross'] = cxfs.verification_dir(dirname=jsim, path=params_path['dirpath_output_cross'], clear=False, verbose=False)                    
                cxfs.savedata_from_dict(Cl_=np.vstack(( np.arange(_cx_.shape[1]), _cx_ )), filename='cl_cx_field1_fg_rec_rec', path=params_path['dirpath_output_jcross'], 
                                            prefix=params_general['prefix'], clear=False)   
                if params_general['verbose']: print('Saved ({}-{}) at    : {}'.format(jsim, 'cl_cx_field1_fg_rec_rec', params_path['dirpath_output_jcross'])) #timer?
                if not j: 
                    params_vars['cl_cx_field1_fg_rec_rec_matrix'] = dcopy(_cx_.flatten())
                else:
                    params_vars['cl_cx_field1_fg_rec_rec_matrix'] = np.vstack(( params_vars['cl_cx_field1_fg_rec_rec_matrix'],dcopy(_cx_.flatten()) ))#_CX_3 = np.vstack(( _CX_3, dcopy(_cx_) ))                
                time_vec = np.hstack((time_vec,time.time()-timej))
                del _cx_,_mfg_,_alm_fg_rec_, timej
            if params_flags['cl_cx_field1_field2_rec_sim']:
                timej  = time.time()
                _cx_  = np.vstack([hp.alm2cl(alms1=_alm_hi_rec_[jch], alms2=params_vars['alm_field2_sim_matrix'][j]) for jch in range(params_general['nch'])])
                params_path['dirpath_output_cross']  = cxfs.verification_dir(dirname='cross', path=params_path['dirpath_output']   , clear=False, verbose=False)
                params_path['dirpath_output_jcross'] = cxfs.verification_dir(dirname=jsim, path=params_path['dirpath_output_cross'], clear=False, verbose=False)                    
                cxfs.savedata_from_dict(Cl_=np.vstack(( np.arange(_cx_.shape[1]), _cx_ )), filename='cl_cx_field1_field2_rec_sim', path=params_path['dirpath_output_jcross'], 
                                            prefix=params_general['prefix'], clear=False)                      
                if params_general['verbose']: print('Saved ({}-{}) at: {}'.format(jsim, 'cl_cx_field1_field2_rec_sim', params_path['dirpath_output_jcross'])) #timer?
                if not j: 
                    params_vars['cl_cx_field1_field2_rec_sim_matrix'] = dcopy(_cx_.flatten())
                else:
                    params_vars['cl_cx_field1_field2_rec_sim_matrix'] = np.vstack(( params_vars['cl_cx_field1_field2_rec_sim_matrix'],dcopy(_cx_.flatten()) ))#_CX_3 = np.vstack(( _CX_3, dcopy(_cx_) ))                
                time_vec = np.hstack((time_vec,time.time()-timej))
                del _cx_,timej
            del _alm_hi_rec_

        if params_flags['cl_field1_lkg'] or params_flags['cl_fg_lkg']\
           or params_flags['filter_fg'] or params_flags['cl_cx_field1_fg_lkg_sim'] \
           or params_flags['cl_cx_field1_fg_lkg_lkg'] or params_flags['cl_cx_field1_field2_lkg_sim']:
            _A_  = hdata.getmap(dirpath_=params_path['dirpath_jrec'], filename_=jmixname, healpix_readingformat=0 , hdu=1)
            _w_   = np.dot( np.linalg.inv(np.dot(_A_.T,_A_)), _A_.T)
            _w_ = np.dot(_A_,_w_)     
            hi_lkg_jsim = np.dot(_w_, hp.alm2map(alms=params_vars['alm_field1_sim_matrix'][j].reshape(params_general['nch'],-1), 
                                                 nside=params_general['nside'],pol=False))
            alm_hi_lkg_jsim = hp.map2alm(maps=hi_lkg_jsim, pol=False)
            del _A_, hi_lkg_jsim
            if params_flags['filter_fg']:
                timej  = time.time()
                params_path['dirpath_output_rec']  = cxfs.verification_dir(dirname='estimated', path=params_path['dirpath_output'], clear=False, verbose=False)
                params_path['dirpath_output_jrec'] = cxfs.verification_dir(dirname=jsim, path=params_path['dirpath_output_rec']   , clear=False, verbose=False)
                cxfs.save_matrix(matrix_=_w_,  filename="filter_fg", path=params_path['dirpath_output_jrec'], prefix=params_general['prefix'], clear=False)         
                if params_general['verbose']: print('Saved ({}-{}) at                  : {}'.format(jsim, 'filter_fg', params_path['dirpath_output_jrec'])) #timer?
                params_vars['filter_fg_matrix']    =  dcopy(_w_.flatten()) if not j else np.vstack(( params_vars['filter_fg_matrix'], dcopy(_w_.flatten()) ))
                time_vec = np.hstack((time_vec,time.time()-timej))
                del timej
            if params_flags['cl_field1_lkg']:
                timej  = time.time()                
                _cx_ = np.vstack([ hp.alm2cl(alm_hi_lkg_jsim[jch]) for jch in range(params_general['nch']) ])
                params_path['dirpath_output_lkg']  = cxfs.verification_dir(dirname='leakage', path=params_path['dirpath_output'], clear=False, verbose=False)
                params_path['dirpath_output_jlkg'] = cxfs.verification_dir(dirname=jsim, path=params_path['dirpath_output_lkg'] , clear=False, verbose=False)
                cxfs.savedata_from_dict(Cl_=np.vstack(( np.arange(_cx_.shape[1]), _cx_ )), filename='cl_field1_lkg', path=params_path['dirpath_output_jlkg'], 
                                        prefix=params_general['prefix'], clear=False)        
                if params_general['verbose']: print('Saved ({}-{}) at              : {}'.format(jsim, 'cl_field1_lkg', params_path['dirpath_output_jlkg'])) #timer?      
                params_vars['cl_field1_lkg_matrix'] = dcopy(_cx_.flatten()) if not j else np.vstack(( params_vars['cl_field1_lkg_matrix'],dcopy(_cx_.flatten()) ))
                time_vec = np.hstack((time_vec,time.time()-timej))
                del _cx_,timej     
            if params_flags['cl_cx_field1_field2_lkg_sim']:
                timej  = time.time()             
                _cx_  = np.vstack([ hp.alm2cl(alms1=alm_hi_lkg_jsim[jch], alms2=params_vars['alm_field2_sim_matrix'][j])  for jch in range(params_general['nch'])])
                params_path['dirpath_output_cross']  = cxfs.verification_dir(dirname='cross', path=params_path['dirpath_output']   , clear=False, verbose=False)
                params_path['dirpath_output_jcross'] = cxfs.verification_dir(dirname=jsim, path=params_path['dirpath_output_cross'], clear=False, verbose=False)                    
                cxfs.savedata_from_dict(Cl_=np.vstack(( np.arange(_cx_.shape[1]), _cx_ )), filename='cl_cx_field1_field2_lkg_sim', path=params_path['dirpath_output_jcross'], 
                                        prefix=params_general['prefix'], clear=False)        
                if params_general['verbose']: print('Saved ({}-{}) at: {}'.format(jsim, 'cl_cx_field1_field2_lkg_sim', params_path['dirpath_output_jcross'])) #timer?
                if not j: 
                    params_vars['cl_cx_field1_field2_lkg_sim_matrix'] = dcopy(_cx_.flatten())
                else:
                    params_vars['cl_cx_field1_field2_lkg_sim_matrix'] = np.vstack(( params_vars['cl_cx_field1_field2_lkg_sim_matrix'],dcopy(_cx_.flatten()) ))               
                time_vec = np.hstack((time_vec,time.time()-timej))
                del _cx_, timej                                
            if params_flags['cl_cx_field1_fg_lkg_sim']:
                timej  = time.time()                        
                _cx_  = np.vstack([ hp.alm2cl(alms1=alm_hi_lkg_jsim[jch], alms2=params_vars['alm_fg_sim_matrix'].reshape(params_general['nch'],-1)[jch])  for jch in range(params_general['nch'])])                
                params_path['dirpath_output_cross']  = cxfs.verification_dir(dirname='cross', path=params_path['dirpath_output']   , clear=False, verbose=False)
                params_path['dirpath_output_jcross'] = cxfs.verification_dir(dirname=jsim, path=params_path['dirpath_output_cross'], clear=False, verbose=False)
                cxfs.savedata_from_dict(Cl_=np.vstack(( np.arange(_cx_.shape[1]), _cx_ )), filename='cl_cx_field1_fg_lkg_sim', path=params_path['dirpath_output_jcross'], 
                                        prefix=params_general['prefix'], clear=False)        
                if not j: 
                    params_vars['cl_cx_field1_fg_lkg_sim_matrix'] = dcopy(_cx_.flatten())
                else:
                    params_vars['cl_cx_field1_fg_lkg_sim_matrix'] = np.vstack(( params_vars['cl_cx_field1_fg_lkg_sim_matrix'], dcopy(_cx_.flatten()) ))                
                if params_general['verbose']: print('Saved ({}-{}) at    : {}'.format(jsim, 'cl_cx_field1_fg_lkg_sim', params_path['dirpath_output_jcross'])) #timer?                                    
                time_vec = np.hstack((time_vec,time.time()-timej))                
                del _cx_,timej
            if params_flags['cl_cx_field1_fg_lkg_lkg'] or params_flags['cl_fg_lkg']:
                fg_lkg_j     = np.dot(1-_w_, hp.alm2map(alms=params_vars['alm_fg_sim_matrix'].reshape(params_general['nch'],-1), nside=params_general['nside'],pol=False))
                alm_fg_lkg_j = hp.map2alm(maps=fg_lkg_j, pol=False)   
                del fg_lkg_j
                if params_flags['cl_fg_lkg']: 
                    timej  = time.time()                        
                    _cx_  = np.vstack([hp.alm2cl(alm_fg_lkg_j[jch]) for jch in range(params_general['nch'])])
                    params_path['dirpath_output_lkg']  = cxfs.verification_dir(dirname='leakage', path=params_path['dirpath_output'], clear=False, verbose=False)
                    params_path['dirpath_output_jlkg'] = cxfs.verification_dir(dirname=jsim, path=params_path['dirpath_output_lkg'] , clear=False, verbose=False)                  
                    cxfs.savedata_from_dict(Cl_=np.vstack(( np.arange(_cx_.shape[1]), _cx_ )), filename='cl_fg_lkg', path=params_path['dirpath_output_jlkg'], 
                                            prefix=params_general['prefix'], clear=False)     
                    params_vars['cl_fg_lkg_matrix'] = dcopy(_cx_.flatten()) if not j else np.vstack(( params_vars['cl_fg_lkg_matrix'], dcopy(_cx_.flatten()) ))
                    if params_general['verbose']: print('Saved ({}-{}) at                  : {}'.format(jsim, 'cl_fg_lkg', params_path['dirpath_output_jlkg'])) #timer?                
                    time_vec = np.hstack((time_vec,time.time()-timej))                
                    del _cx_,timej                      
                if params_flags['cl_cx_field1_fg_lkg_lkg']:                 
                    timej  = time.time()                          
                    _cx_  = np.vstack([ hp.alm2cl(alms1=alm_hi_lkg_jsim[jch], alms2=alm_fg_lkg_j[jch])  for jch in range(params_general['nch'])])
                    params_path['dirpath_output_cross']  = cxfs.verification_dir(dirname='cross', path=params_path['dirpath_output'],    clear=False, verbose=False)
                    params_path['dirpath_output_jcross'] = cxfs.verification_dir(dirname=jsim, path=params_path['dirpath_output_cross'], clear=False, verbose=False)
                    cxfs.savedata_from_dict(Cl_=np.vstack(( np.arange(_cx_.shape[1]), _cx_ )), filename='cl_cx_field1_fg_lkg_lkg', path=params_path['dirpath_output_jcross'], 
                                            prefix=params_general['prefix'], clear=False)        
                    if not j: 
                        params_vars['cl_cx_field1_fg_lkg_lkg_matrix'] = dcopy(_cx_.flatten())
                    else:
                        params_vars['cl_cx_field1_fg_lkg_lkg_matrix'] = np.vstack(( params_vars['cl_cx_field1_fg_lkg_lkg_matrix'], dcopy(_cx_.flatten()) ))               
                    if params_general['verbose']: print('Saved ({}-{}) at    : {}'.format(jsim, 'cl_cx_field1_fg_lkg_lkg', params_path['dirpath_output_jcross'])) #timer?                                    
                    time_vec = np.hstack((time_vec,time.time()-timej))                                        
                    del _cx_                                     
                del alm_fg_lkg_j
        time_matrix = np.around(time_vec,decimals=4) if not j else np.vstack((time_matrix, np.around(time_vec,decimals=4)))
        del jhifilename,jfgfilename,jname,jsim,j,jmixname,alm_hi_lkg_jsim,_w_
if params_general['verbose']:print(np.around(np.average(time_matrix,axis=0),decimals=3)[1:])               
if params_general['verbose']:print('Loop time: {0:.4f} seg'.format(time.time()-timeii)) 
del timeii

###############################
### PART5 
print('\n------------------------------------------------------------------')
timeii = time.time()
if params_flags['filter_fg']:
    params_vars['filter_fg_matrix']    = np.average( params_vars['filter_fg_matrix'],axis=0).reshape(params_general['nch'],params_general['nch'])
    params_path['dirpath_output_rec']  = cxfs.verification_dir(dirname='estimated', path=params_path['dirpath_output'], clear=False, verbose=False)
    params_path['dirpath_output_mean'] = cxfs.verification_dir(dirname='mean', path=params_path['dirpath_output_rec'], clear=False, verbose=False)
    if params_general['verbose']: print('Saved ({}-{}) at                  : {}'.format('mean', 'filter_fg', params_path['dirpath_output_mean'])) #timer?                                    
    cxfs.save_matrix(matrix_= params_vars['filter_fg_matrix'], filename='filter_fg_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False) 
    del params_vars['filter_fg_matrix']
                  
if params_flags['cl_field1_rec']:
    params_vars['cl_field1_rec_mean'] = np.vstack(( np.arange(3*params_general['nside']), 
                                                                np.average(params_vars['cl_field1_rec_matrix'], axis=0).reshape(params_general['nch'],-1) 
                                                 ))
    params_path['dirpath_output_rec']   = cxfs.verification_dir(dirname='estimated', path=params_path['dirpath_output'], clear=False, verbose=False)
    params_path['dirpath_output_mean']  = cxfs.verification_dir(dirname='mean',  path=params_path['dirpath_output_rec'], clear=False, verbose=False)
    cxfs.savedata_from_dict(Cl_=params_vars['cl_field1_rec_mean'], filename='cl_field1_rec_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False) 
    if params_general['verbose']: print('Saved ({}-{}) at              : {}'.format('mean', 'cl_field1_rec', params_path['dirpath_output_mean'])) #timer?                                        
    if not params_flags['alpha_field1']: del params_vars['cl_field1_rec_matrix']
    del params_vars['cl_field1_rec_mean']

if params_flags['cl_cx_field1_fg_rec_rec']:
    params_vars['cl_cx_field1_fg_rec_rec_mean'] = np.vstack(( np.arange(3*params_general['nside']), 
                                                              np.average(params_vars['cl_cx_field1_fg_rec_rec_matrix'], axis=0).reshape(params_general['nch'],-1) 
                                                    ))
    params_path['dirpath_output_cross'] = cxfs.verification_dir(dirname='cross', path=params_path['dirpath_output'], clear=False, verbose=False)
    params_path['dirpath_output_mean']  = cxfs.verification_dir(dirname='mean',  path=params_path['dirpath_output_cross'], clear=False, verbose=False)
    cxfs.savedata_from_dict(Cl_=params_vars['cl_cx_field1_fg_rec_rec_mean'], filename='cl_cx_field1_fg_rec_rec_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False) 
    if params_general['verbose']: print('Saved ({}-{}) at    : {}'.format('mean', 'cl_cx_field1_fg_rec_rec', params_path['dirpath_output_mean'])) #timer?                                            
    del params_vars['cl_cx_field1_fg_rec_rec_mean'], params_vars['cl_cx_field1_fg_rec_rec_matrix']    

if params_flags['cl_cx_field1_field2_rec_sim']:
    params_vars['cl_cx_field1_field2_rec_sim_mean'] = np.vstack(( np.arange(int(params_vars['cl_cx_field1_field2_rec_sim_matrix'].shape[1]/params_general['nch'])), 
                                                                  np.average(params_vars['cl_cx_field1_field2_rec_sim_matrix'], axis=0).reshape(params_general['nch'],-1) 
                                                    ))
    params_path['dirpath_output_cross'] = cxfs.verification_dir(dirname='cross', path=params_path['dirpath_output'], clear=False, verbose=False)
    params_path['dirpath_output_mean' ] = cxfs.verification_dir(dirname='mean',  path=params_path['dirpath_output_cross'], clear=False, verbose=False)
    cxfs.savedata_from_dict(Cl_=params_vars['cl_cx_field1_field2_rec_sim_mean'], filename='cl_cx_field1_field2_rec_sim_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False)     
    if params_general['verbose']: print('Saved ({}-{}) at: {}'.format('mean', 'cl_cx_field1_field2_rec_sim', params_path['dirpath_output_mean'])) #timer?                                            
    if not params_flags['alpha_cx']: del params_vars['cl_cx_field1_field2_rec_sim_matrix']
    del params_vars['cl_cx_field1_field2_rec_sim_mean']

if params_flags['cl_field1_lkg']:
    params_vars['cl_field1_lkg_mean'] = np.vstack(( np.arange(int(params_vars['cl_field1_lkg_matrix'].shape[1]/params_general['nch'])), 
                                                    np.average(params_vars['cl_field1_lkg_matrix'], axis=0).reshape(params_general['nch'],-1)  ))
    params_path['dirpath_output_lkg']  = cxfs.verification_dir(dirname='leakage',   path=params_path['dirpath_output'], clear=False, verbose=False)
    params_path['dirpath_output_mean'] = cxfs.verification_dir(dirname='mean',  path=params_path['dirpath_output_lkg'], clear=False, verbose=False)
    cxfs.savedata_from_dict(Cl_=params_vars['cl_field1_lkg_mean'], filename='cl_field1_lkg_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False)      
    if params_general['verbose']: print('Saved ({}-{}) at              : {}'.format('mean', 'cl_field1_lkg', params_path['dirpath_output_mean'])) #timer?                                                
    del params_vars['cl_field1_lkg_matrix'],params_vars['cl_field1_lkg_mean']

if params_flags['cl_cx_field1_field2_lkg_sim']:
    params_vars['cl_cx_field1_field2_lkg_sim_mean'] = np.vstack(( np.arange(int(params_vars['cl_cx_field1_field2_lkg_sim_matrix'].shape[1]/params_general['nch'])), 
                                                                  np.average(params_vars['cl_cx_field1_field2_lkg_sim_matrix'], axis=0).reshape(params_general['nch'],-1)  ))
    params_path['dirpath_output_cross']= cxfs.verification_dir(dirname='cross', path=params_path['dirpath_output']      , clear=False, verbose=False)
    params_path['dirpath_output_mean'] = cxfs.verification_dir(dirname='mean' , path=params_path['dirpath_output_cross'], clear=False, verbose=False)   
    cxfs.savedata_from_dict(Cl_=params_vars['cl_cx_field1_field2_lkg_sim_mean'], filename='cl_cx_field1_field2_lkg_sim_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False)       
    if params_general['verbose']: print('Saved ({}-{}) at: {}'.format('mean', 'cl_cx_field1_field2_lkg_sim', params_path['dirpath_output_mean'])) #timer?                                                    
    del params_vars['cl_cx_field1_field2_lkg_sim_matrix'], params_vars['cl_cx_field1_field2_lkg_sim_mean']

if params_flags['cl_cx_field1_fg_lkg_sim']:
    params_vars['cl_cx_field1_fg_lkg_sim_mean'] = np.vstack(( np.arange(int(params_vars['cl_cx_field1_fg_lkg_sim_matrix'].shape[1]/params_general['nch'])), 
                                                              np.average(params_vars['cl_cx_field1_fg_lkg_sim_matrix'], axis=0).reshape(params_general['nch'],-1)  ))
    params_path['dirpath_output_cross']= cxfs.verification_dir(dirname='cross', path=params_path['dirpath_output'],       clear=False, verbose=False)
    params_path['dirpath_output_mean'] = cxfs.verification_dir(dirname='mean' , path=params_path['dirpath_output_cross'], clear=False, verbose=False)    
    cxfs.savedata_from_dict(Cl_=params_vars['cl_cx_field1_fg_lkg_sim_mean'], filename='cl_cx_field1_fg_lkg_sim_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False)      
    if params_general['verbose']: print('Saved ({}-{}) at    : {}'.format('mean', 'cl_cx_field1_fg_lkg_sim', params_path['dirpath_output_mean'])) #timer?                                                        
    del params_vars['cl_cx_field1_fg_lkg_sim_matrix'], params_vars['cl_cx_field1_fg_lkg_sim_mean']

if params_flags['cl_fg_lkg']: 
    params_vars['cl_fg_lkg_mean'] = np.vstack(( np.arange(int(params_vars['cl_fg_lkg_matrix'].shape[1]/params_general['nch'])),  
                                                np.average(params_vars['cl_fg_lkg_matrix'], axis=0).reshape(params_general['nch'],-1)  ))
    params_path['dirpath_output_lkg']  = cxfs.verification_dir(dirname='leakage',   path=params_path['dirpath_output'], clear=False, verbose=False)
    params_path['dirpath_output_mean'] = cxfs.verification_dir(dirname='mean',  path=params_path['dirpath_output_lkg'], clear=False, verbose=False)
    cxfs.savedata_from_dict(Cl_=params_vars['cl_fg_lkg_mean'], filename='cl_fg_lkg_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False)     
    if params_general['verbose']: print('Saved ({}-{}) at                  : {}'.format('mean', 'cl_fg_lkg', params_path['dirpath_output_mean'])) #timer?                                                            
    del params_vars['cl_fg_lkg_matrix'], params_vars['cl_fg_lkg_mean']


if params_flags['cl_cx_field1_fg_lkg_lkg']:  
    params_vars['cl_cx_field1_fg_lkg_lkg_mean'] = np.vstack(( np.arange(int(params_vars['cl_cx_field1_fg_lkg_lkg_matrix'].shape[1]/params_general['nch'])),  
                                                              np.average(params_vars['cl_cx_field1_fg_lkg_lkg_matrix'], axis=0).reshape(params_general['nch'],-1)  ))
    params_path['dirpath_output_cross'] = cxfs.verification_dir(dirname='cross', path=params_path['dirpath_output']      , clear=False, verbose=False)
    params_path['dirpath_output_mean']  = cxfs.verification_dir(dirname='mean',  path=params_path['dirpath_output_cross'], clear=False, verbose=False)
    cxfs.savedata_from_dict(Cl_= params_vars['cl_cx_field1_fg_lkg_lkg_mean'], filename='cl_cx_field1_fg_lkg_lkg_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False)       
    if params_general['verbose']: print('Saved ({}-{}) at    : {}'.format('mean', 'cl_cx_field1_fg_lkg_lkg', params_path['dirpath_output_mean'])) #timer?                                                                
    del params_vars['cl_cx_field1_fg_lkg_lkg_matrix'], params_vars['cl_cx_field1_fg_lkg_lkg_mean']
if params_general['verbose']:print('Loop time: {0:.4f} seg'.format(time.time()-timeii)) 
del timeii
###############################
### PART5
print('\n------------------------------------------------------------------')
timej  = time.time()
if params_flags['cl_fg_rec']:
    if params_general['verbose']: 
        print('\n------------------------------------------------------------------')
        print('Saving {} ...'.format('cl_fg_rec'))        
    params_path['dirpath_output_rec']   = cxfs.verification_dir(dirname='estimated', path=params_path['dirpath_output'], clear=False, verbose=False)
    for j,jsim in enumerate(params_path['sim_names']):
        params_path['dirpath_jrec'] = os.path.join(params_path['dirpath_rec'],jsim)
        for jname in os.listdir(params_path['dirpath_jrec']):
            if 'extFG'     in jname: jfgfilename = jname
        _mfg_ = hdata.getmap(dirpath_=params_path['dirpath_jrec'], filename_=jfgfilename, healpix_readingformat=0, hdu=1)
        _cx_  = np.vstack([ hp.anafast(map1=_mfg_[jch], pol=False)  for jch in range(params_general['nch'])]) 
        #params_path['dirpath_output_rec']  = cxfs.verification_dir(dirname='estimated', path=params_path['dirpath_output'], clear=False, verbose=False)
        params_path['dirpath_output_jrec'] = cxfs.verification_dir(dirname=jsim   , path=params_path['dirpath_output_rec'], clear=False, verbose=False)  
        cxfs.savedata_from_dict(Cl_=np.vstack(( np.arange(_cx_.shape[1]), _cx_ )), 
                                filename='cl_fg_rec', path=params_path['dirpath_output_jrec'], 
                                            prefix=params_general['prefix'], clear=False)        
        if not j: _CX_ = dcopy(_cx_.flatten())
        else:     _CX_ = np.vstack(( _CX_, dcopy(_cx_.flatten()) ))
        if params_general['verbose']: print('Saved ({}) at: {}'.format(jsim, params_path['dirpath_output_jrec'])) 
        del _cx_, _mfg_, jname, jsim, j
    params_vars['cl_fg_rec_matrix'] = _CX_
    _CX_ = np.vstack(( np.arange(_CX_.shape[1]), _CX_ ))
    params_vars['cl_fg_rec_mean'] = np.vstack(( np.arange( params_vars['cl_fg_rec_matrix'].shape[1]), params_vars['cl_fg_rec_matrix'] ))
    params_path['dirpath_output_rec']  = cxfs.verification_dir(dirname='estimated', path=params_path['dirpath_output'], clear=False, verbose=False)
    params_path['dirpath_output_mean'] = cxfs.verification_dir(dirname='mean'     , path=params_path['dirpath_output_rec'], clear=False, verbose=False)
    cxfs.savedata_from_dict(Cl_=_CX_, filename='cl_fg_rec_mean', path=params_path['dirpath_output_mean'], 
                            prefix=params_general['prefix'], clear=False)   
    if params_general['verbose']: print('Saved (mean) at: {}'.format(params_path['dirpath_output_mean'])) 
    if params_general['verbose']: print('Processing time: {0:.4f} seg'.format(time.time()-timej))               
    del params_vars['cl_fg_rec_mean'], _CX_, params_path['dirpath_output_mean'], params_path['dirpath_output_rec']

del params_vars['map_fg_sim'], params_vars['alm_field1_sim_matrix'], params_vars['alm_field2_sim_matrix'],params_vars['alm_fg_sim_matrix']
###############################
### PART6
if params_flags['alpha_field1']:
    if params_general['verbose']: 
        print('\n------------------------------------------------------------------')
        print('Saving {} ...'.format('alpha_field1'))     
    alpha = np.average(params_vars[ 'cl_field1_sim_matrix'][1:,:]/params_vars[ 'cl_field1_rec_matrix'][1:,:],axis=0).reshape(params_general['nch'],-1) 
    params_path['dirpath_output_alpha'] = cxfs.verification_dir(dirname='alpha', path=params_path['dirpath_output'], clear=False, verbose=False)
    cxfs.save_matrix(matrix_= alpha, filename='alpha_field1', path=params_path['dirpath_output_alpha'], 
                     header= "<E(Field1)/Field1>  [each column is a differet channel/frequency]",
                            prefix=params_general['prefix'], clear=False) 
    if params_general['verbose']: print('Saved at: {}'.format(params_path['dirpath_output_alpha']))     

if params_flags['alpha_fg']:
    if params_general['verbose']: 
        print('\n------------------------------------------------------------------')
        print('Saving {} ...'.format('alpha_fg'))       
    alpha = np.average(params_vars['cl_fg_rec_matrix']/params_vars['cl_fg_sim_matrix'][1:,:].flatten()[np.newaxis,:],axis=0).reshape(params_general['nch'],-1) 
    params_path['dirpath_output_alpha'] = cxfs.verification_dir(dirname='alpha', path=params_path['dirpath_output'], clear=False, verbose=False)
    cxfs.save_matrix(matrix_= alpha, filename='alpha_fg', path=params_path['dirpath_output_alpha'], 
                     header= "<E(FG)/FG>  [each column is a differet channel/frequency]",
                            prefix=params_general['prefix'], clear=False) 
    if params_general['verbose']: print('Saved at: {}'.format(params_path['dirpath_output_alpha'])) 

if params_flags['alpha_cx']:
    if params_general['verbose']: 
        print('\n------------------------------------------------------------------')
        print('Saving {} ...'.format('alpha_cx'))          
    alpha = np.average(params_vars[ 'cl_cx_field1_field2_sim_sim_matrix'][1:,:]/params_vars[ 'cl_cx_field1_field2_rec_sim_matrix'][1:,:],axis=0).reshape(params_general['nch'],-1) 
    params_path['dirpath_output_alpha'] = cxfs.verification_dir(dirname='alpha', path=params_path['dirpath_output'], clear=False, verbose=False)
    cxfs.save_matrix(matrix_= alpha, filename='alpha_cx', path=params_path['dirpath_output_alpha'], 
                     header= "<E(Field1)xField2/Field1xField2>  [each column is a differet channel/frequency]",
                     prefix=params_general['prefix'], clear=False) 
    if params_general['verbose']: print('Saved at: {}'.format(params_path['dirpath_output_alpha']))     
    

###############################
### OVER                    
print('\n------------------------------------------------------------------')
if params_general['verbose']: print('Total processing time: {0:.2f} seg'.format(time.time()-timei)) 
if 1:
    print('===================================================================')    
    print('END')
    print('===================================================================\n')

