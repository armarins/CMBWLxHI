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
#sys.path.insert(1, '/data/AMARINS/CMBWLxHI-CODES/scripts')
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
INI         = "generating_APS_leakage_estimation.ini"
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
ns            = config.getint(    "dataset","ns" )
nu_min_correlated = config.getint("dataset","nu_min_correlated" )
nu_max_correlated = config.getint("dataset","nu_max_correlated" )
#######
#PATH
#theory
filepath_cross    = config.get("PATH","filepath_cross" )
#sims
dirpath_sims      = config.get("PATH","dirpath_sims" )
#estimated
dirpath_estimated = config.get("PATH","dirpath_estimated" )
#postprocessing
dirpath_postprocessing = config.get("PATH","dirpath_postprocessing" )
#output_dir
dirpath_output    = config.get("PATH","dirpath_output" )

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
parser.add_argument('--ns'           , action = 'store', dest = 'ns'           , default = ns           , help = '')
parser.add_argument('--nu_min_correlated', action = 'store', dest = 'nu_min_correlated', default = nu_min_correlated, help = '')
parser.add_argument('--nu_max_correlated', action = 'store', dest = 'nu_max_correlated', default = nu_max_correlated, help = '')
######
#PATH
#theory
parser.add_argument('--filepath_cross'      , action = 'store', dest = 'filepath_cross'      , default = filepath_cross      , help = '')
#sims
parser.add_argument('--dirpath_sims'        , action = 'store', dest = 'dirpath_sims'        , default = dirpath_sims        , help = '')
#estimated
parser.add_argument('--dirpath_estimated'   , action = 'store', dest = 'dirpath_estimated'   , default = dirpath_estimated   , help = '')
#postprocessing
parser.add_argument('--dirpath_postprocessing' , action = 'store', dest = 'dirpath_postprocessing' , default = dirpath_postprocessing , help = '')
#output_dir
parser.add_argument('--dirpath_output'      , action = 'store', dest = 'dirpath_output'      , default = dirpath_output      , help = '')


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
ns            = int( arguments.ns)
nu_min_correlated = int( arguments.nu_min_correlated )
nu_max_correlated = int( arguments.nu_max_correlated )
#########
#PATH
#theory
filepath_cross  = str(arguments.filepath_cross)
#sims
dirpath_sims = str(arguments.dirpath_sims)
#estimated
dirpath_estimated      = str(arguments.dirpath_estimated)
#postprocessing
dirpath_postprocessing = str(arguments.dirpath_postprocessing)
#output_dir
dirpath_output         = str(arguments.dirpath_output)

dirpath_estimated      = os.path.join(dirpath_estimated     , "ns{}".format(ns))
dirpath_output         = os.path.join(dirpath_output        , "ns{}".format(ns))
dirpath_postprocessing = os.path.join(dirpath_postprocessing, "ns{}".format(ns))
#####################################################################################################################################################################
#####################################################################################################################################################################
params_general = {'project':project, 'verbose':verbose,'prefix':prefix, 
                  'nreals':nrealizations, 'ns':ns, 
                  'nu_min_correlated':nu_min_correlated, 'nu_max_correlated':nu_max_correlated}

params_path = {'filepath_cross':filepath_cross,
               'dirpath_sims'  :dirpath_sims,    'dirpath_estimated':dirpath_estimated, 
               'dirpath_output':dirpath_output,  'dirpath_postprocessing':dirpath_postprocessing, }


if params_general['prefix']=="":params_general['prefix']=None
del project, prefix, nrealizations, ns, verbose
del filepath_cross, dirpath_sims, dirpath_estimated
del dirpath_output, dirpath_postprocessing


if params_general['verbose']:
    print('\n--Variables to be used:--')
    print(params_general)    
    print(params_path)

#########################################
params_vars = {'cl_cx':np.loadtxt(params_path['filepath_cross'] ).T[1:,:] }
params_general['nch'] = int(params_vars['cl_cx'].shape[0])

###############################
sim_names = []
#for j,jsim in enumerate(np.sort( np.asarray( os.listdir(params_path['dirpath_sims']) ) )):
for j,jsim in enumerate(np.sort( np.asarray( os.listdir(params_path['dirpath_estimated']) ) )):
    if 'sim' in jsim:
        sim_names.append(jsim)
params_vars['sim_names'] = np.asarray(sim_names)
del sim_names,j,jsim

for j,jsim in enumerate(params_vars['sim_names'][:params_general['nreals']+1]):
    alm_kp_sim_imag = np.loadtxt(os.path.join(params_path['dirpath_sims'],jsim, 'alm_kappa_sim_imag.txt')).T[2]
    alm_kp_sim_real = np.loadtxt(os.path.join(params_path['dirpath_sims'],jsim, 'alm_kappa_sim_real.txt')).T[2] 
    _cx_            = np.loadtxt(os.path.join(params_path['dirpath_sims'],jsim, 'cl_cross_sim.txt')).T[1:,:]
    alm_kp_sim = np.array([complex(jreal,jcompl) for (jreal,jcompl) in zip(alm_kp_sim_real, alm_kp_sim_imag)])
    alm_kp_sim = np.ascontiguousarray(alm_kp_sim)
    del alm_kp_sim_imag, alm_kp_sim_real
    alm_field2_sim_matrix = alm_kp_sim if not j else np.vstack(( alm_field2_sim_matrix, alm_kp_sim ))
    if not j: 
        cl_cx_field1_field2_sim_sim_matrix = dcopy(_cx_.flatten())
    else:
        cl_cx_field1_field2_sim_sim_matrix = np.vstack(( cl_cx_field1_field2_sim_sim_matrix, dcopy(_cx_.flatten()) ))     
params_vars['cl_cx_f1f2_sim_sim_matrix'] = cl_cx_field1_field2_sim_sim_matrix
params_vars['alm_f2_sim_matrix']         = alm_field2_sim_matrix
del alm_kp_sim,jsim,cl_cx_field1_field2_sim_sim_matrix, alm_field2_sim_matrix



for j,jsim in enumerate(params_vars['sim_names'][:params_general['nreals']+1]):
    dirpath_jrec = os.path.join(params_path['dirpath_estimated'],jsim)
    for jname in os.listdir(dirpath_jrec):
        if 'extHI'     in jname: jhifilename = jname
        if 'mixmatrix' in jname: jmixname    = jname
    _mp_         = hdata.getmap(dirpath_=dirpath_jrec, filename_=jhifilename, healpix_readingformat=0, hdu=1)
    _alm_hi_rec_ = hp.map2alm(maps=_mp_, pol=False)
    _cx_         = np.vstack([hp.alm2cl(alms1=_alm_hi_rec_[jch], alms2=params_vars['alm_f2_sim_matrix'][j]) for jch in range(params_general['nch'])])

    _A_ = hdata.getmap(dirpath_=dirpath_jrec, filename_=jmixname, healpix_readingformat=0 , hdu=1)
    _w_ = np.dot( np.linalg.inv(np.dot(_A_.T,_A_)), _A_.T)
    _w_ = np.dot(_A_,_w_) 
    
    if not j: 
        cl_cx_field1_field2_rec_sim_matrix = dcopy(_cx_.flatten())
        filter_matrix = dcopy(_w_.flatten())
    else:
        cl_cx_field1_field2_rec_sim_matrix = np.vstack(( cl_cx_field1_field2_rec_sim_matrix, dcopy(_cx_.flatten()) ))    
        filter_matrix = np.vstack(( filter_matrix, dcopy(_w_.flatten()) ))    
params_vars['cl_cx_f1f2_rec_sim_matrix'] = cl_cx_field1_field2_rec_sim_matrix
params_vars['filter_matrix'] = filter_matrix
del _w_, _A_, _cx_, _mp_
del _alm_hi_rec_,cl_cx_field1_field2_rec_sim_matrix,filter_matrix,
del dirpath_jrec,jname,jsim,j,jhifilename,jmixname

params_vars['cl_cx_f1f2_sim_sim_mean'] = np.average(params_vars['cl_cx_f1f2_sim_sim_matrix'][1:,:],axis=0).reshape(params_general['nch'],-1)
#params_vars['cl_cx_f1f2_rec_sim_L0'] = params_vars['cl_cx_f1f2_rec_sim_matrix'][0].reshape(params_general['nch'],-1)

cl_cx_lkg_th  = np.zeros_like(params_vars['cl_cx_f1f2_sim_sim_mean'])
cl_cx_lkg_sim = np.zeros_like(params_vars['cl_cx_f1f2_sim_sim_mean'])
params_vars['filter_matrix_L0'] = params_vars['filter_matrix'][0].reshape(params_general['nch'],params_general['nch'])  


nu_min,nu_max=params_general['nu_min_correlated'],params_general['nu_max_correlated']
params_vars['cl_cx_th'] = np.zeros_like(params_vars['cl_cx'])
params_vars['cl_cx_th'][nu_min:nu_max+1,:] = params_vars['cl_cx'][nu_min:nu_max+1,:]

for j in np.arange(0,params_general['nch'],1):
    cl_cx_lkg_th[j]      = np.dot(params_vars['filter_matrix_L0'] [j,:], params_vars['cl_cx_th'][:,:]) 
    cl_cx_lkg_sim[j]     = np.dot(params_vars['filter_matrix_L0'] [j,:], params_vars['cl_cx_f1f2_sim_sim_mean'][:,:])

print('\n------------------------------------------------------------------')
cxfs.savedata_from_dict(Cl_= cl_cx_lkg_th, filename='leakage_theory', path=params_path['dirpath_output'], 
                        prefix=params_general['prefix'], clear=False)       
if params_general['verbose']: print('Saved ({}) at    : {}'.format('leakage_theory', params_path['dirpath_output'])) #timer?        

cxfs.savedata_from_dict(Cl_= cl_cx_lkg_sim, filename='leakage_sim', path=params_path['dirpath_output'], 
                        prefix=params_general['prefix'], clear=False)       
if params_general['verbose']: print('Saved ({}) at    : {}'.format('leakage_theory', params_path['dirpath_output'])) #timer?        

###############################
### OVER                    
print('\n------------------------------------------------------------------')
if params_general['verbose']: print('Total processing time: {0:.2f} seg'.format(time.time()-timei)) 
if 1:
    print('===================================================================')    
    print('END')
    print('===================================================================\n')






