#sys.path.insert(1, os.getcwd())
import os, sys, time
from copy import deepcopy as dcopy
import cross_functions_theory as cxft
import cross_functions_simulations as cxfs
#import cross_functions_simulations_2 as cxfs2
import numpy as np
import healpy as hp
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
PATH        = '/data/AMARINS/CMBWLxHI-CODES/scripts'
INI         = "generating_APS_simulations.ini"
name_params = os.path.join(PATH,INI)
config.read(name_params)
#General
verbose = config.getboolean("General","verbose")
project = config.get(       "General","project")
#PATHs
filepath_field1       = config.get("PATH","filepath_field1")
filepath_field2       = config.get("PATH","filepath_field2")
filepath_cross        = config.get("PATH","filepath_cross" )
filepath_foregrounds  = config.get("PATH","filepath_foregrounds" )
filepath_mask         = config.get("PATH","filepath_mask" )

pathout = config.get("PATH","pathout")
prefix  = config.get("PATH","prefix" )
suffix  = config.get("PATH","suffix" )
#Contaminants
add_foregrounds = config.getboolean("Contaminants","add_foregrounds")
apply_mask      = config.getboolean("Contaminants","apply_mask")
#Simulations
nrealizations = config.getint("Simulations","nrealizations")
seed0         = config.getint("Simulations","seed0")
limited_correlated_channels = config.getboolean("Simulations","limited_correlated_channels")
channel_min_corr            = config.getfloat("Simulations","channel_min_corr")
channel_max_corr            = config.getfloat("Simulations","channel_max_corr")
amplification               = config.getfloat("Simulations","amplification")
channel_tax                 = config.getint(  "Simulations","channel_tax")

###############################################################################
# You can modify any options in the parameters.ini file by the command terminal
###############################################################################
parser = argparse.ArgumentParser(description='Modify by the command terminal parameters in {} file'.format(INI))

#PATHs
parser.add_argument('--verbose'  , action = 'store', dest = 'verbose'  , default = verbose  , help = '')
parser.add_argument('--project'  , action = 'store', dest = 'project'  , default = project  , help = '')

#Cosmology
parser.add_argument('--filepath_field1'     , action = 'store', dest = 'filepath_field1'     , default = filepath_field1     , help = '')
parser.add_argument('--filepath_field2'     , action = 'store', dest = 'filepath_field2'     , default = filepath_field2     , help = '')
parser.add_argument('--filepath_cross'      , action = 'store', dest = 'filepath_cross'      , default = filepath_cross      , help = '')
parser.add_argument('--filepath_foregrounds', action = 'store', dest = 'filepath_foregrounds', default = filepath_foregrounds, help = '')
parser.add_argument('--filepath_mask'       , action = 'store', dest = 'filepath_mask'       , default = filepath_mask       , help = '')
parser.add_argument('--pathout', action = 'store', dest = 'pathout', default = pathout, help = '')
parser.add_argument('--prefix' , action = 'store', dest = 'prefix' , default = prefix , help = '')
parser.add_argument('--suffix' , action = 'store', dest = 'suffix' , default = suffix , help = '')

parser.add_argument('--add_foregrounds', action = 'store', dest = 'add_foregrounds', default = add_foregrounds, help = '')
parser.add_argument('--apply_mask'     , action = 'store', dest = 'apply_mask'     , default = apply_mask     , help = '')

parser.add_argument('--nrealizations', action = 'store', dest = 'nrealizations', default = nrealizations, help = '')
parser.add_argument('--seed0'        , action = 'store', dest = 'seed0'        , default = seed0        , help = '')

parser.add_argument('--limited_correlated_channels', action = 'store', dest = 'limited_correlated_channels', default = limited_correlated_channels , help = '')
parser.add_argument('--channel_min_corr'              , action = 'store', dest = 'channel_min_corr'              , default = channel_min_corr               , help = '')
parser.add_argument('--channel_max_corr'              , action = 'store', dest = 'channel_max_corr'              , default = channel_max_corr               , help = '')
parser.add_argument('--amplification'              , action = 'store', dest = 'amplification'              , default = amplification               , help = '')
parser.add_argument('--channel_tax'                , action = 'store', dest = 'channel_tax'                , default = channel_tax                , help = '')

arguments = parser.parse_args()
###############################################################################
# Variables
###############################################################################
verbose   = bool(arguments.verbose)
project   = str(arguments.project)

filepath_field1      = str(arguments.filepath_field1)
filepath_field2      = str(arguments.filepath_field2)
filepath_cross       = str(arguments.filepath_cross)
filepath_foregrounds = str(arguments.filepath_foregrounds)
filepath_mask        = str(arguments.filepath_mask)

pathout = str(arguments.pathout)
prefix  = str(arguments.prefix)
suffix  = str(arguments.suffix)

add_foregrounds = bool(arguments.add_foregrounds)
apply_mask      = bool(arguments.apply_mask)

nrealizations   = int(arguments.nrealizations)
seed0           = int(arguments.seed0)

limited_correlated_channels = bool(arguments.limited_correlated_channels)

if limited_correlated_channels:
    channel_min_corr = int(arguments.channel_min_corr)
    channel_max_corr = int(arguments.channel_max_corr)
else:
    channel_min_corr = None
    channel_max_corr = None

amplification = float(arguments.amplification)
channel_tax   = int(arguments.channel_tax)
############################################################################################
#####################################################################################################
if prefix=="":
    prefix=None
if suffix=="":
    suffix=None
#####################################################################################################################################################################
#####################################################################################################################################################################
params_general = {'verbose':verbose}


params_path = {'filepath_field1'     : filepath_field1,
               'filepath_field2'     : filepath_field2,
               'filepath_cross'      : filepath_cross,
               'filepath_foregrounds': filepath_foregrounds,
               'filepath_mask'       : filepath_mask,
               'pathout'             : pathout, 
               'prefix'              : prefix,
               'suffix'              : suffix
              }
params_sims = {'add_foregrounds':add_foregrounds,
               'apply_mask'     :apply_mask,
               'nrealizations'  :nrealizations,
               'seed0'          :seed0,
               'limited_correlated_channels':limited_correlated_channels,
               'channel_min_corr':channel_min_corr,
               'channel_max_corr':channel_max_corr,
               'amplification':amplification,
               'channel_tax'  :channel_tax
               }

del verbose, filepath_field1, filepath_field2, filepath_cross, filepath_foregrounds, filepath_mask, pathout, prefix, suffix
del add_foregrounds, apply_mask, nrealizations, seed0, channel_min_corr, channel_max_corr, limited_correlated_channels

#####################################################################################################
#####################################################################################################
#if params_general['verbose']:
if 1:
    print('\n------------------------------------------------------------------')
    print('Starting to produce Angular Power Spectrum theoretical simulations')
    print('------------------------------------------------------------------\n')


if params_general['verbose']:
    print('Loading data from {}'.format(params_path['filepath_field1']))
clf1 = np.loadtxt(params_path['filepath_field1']).T[1:,:]
if params_general['verbose']:
    print('Loading data from {}'.format(params_path['filepath_field2']))
clf2 = np.loadtxt(params_path['filepath_field2']).T[1:,:]
if params_general['verbose']:
    print('Loading data from {}'.format(params_path['filepath_cross']))
clcx = np.loadtxt(params_path['filepath_cross']).T[1:,:]

#clf1_dict = cxft.cls_from_matrix2dict(clf1)
#clf2_dict = cxft.cls_from_matrix2dict(clf2)
#clcx_dict = cxft.cls_from_matrix2dict(clcx)
#####################################################################################################
#####################################################################################################
if params_general['verbose']:
    print('Generating simulations...')
          
#_dict_ = cxfs.dictionary_cross_simulations_quantities_from_dict(clf1_dict_=clf1_dict, clf2_dict_=clf2_dict, clcx_dict_=clcx_dict, 
         #                                                       seed_hi=params_sims['seed0'], seed_k=int(params_sims['seed0']-1000), nsims=params_sims['nrealizations'],
#                                                                show_time=False)                      

params_sims['seed_hi'  ] = params_sims['seed0']
params_sims['seed_kunc'] = params_sims['seed0']-1000
_dict_ = cxfs.dictionary_cross_simulations_quantities_from_matrix(clf1_=clf1, clf2_=clf2, clcx_=clcx, 
                                                                  seed_hi=params_sims['seed_hi'], seed_k=params_sims['seed_kunc'], nsims=params_sims['nrealizations'],
                                                                  channel_min_corr=params_sims['channel_min_corr'],channel_max_corr=params_sims['channel_max_corr'], 
                                                                  fact=params_sims['amplification'], tax=params_sims['channel_tax'], beta=None,
                                                                  show_time=False)                                  
#################################
# LOOP        ###################
#################################
for i,_sim_ in enumerate(_dict_.keys()):   
    params_path['pathout_dir'] = os.path.join( params_path['pathout'], '{ID}/{s}'.format(ID=project, s=_sim_) )
    for j,_var_ in enumerate(_dict_[_sim_].keys()):    
        if not i+j: 
            clear=True
        else:
            clear=False        
        #######
        #### CL
        if ('cl' in  _var_) or ('cross' in  _var_):
            header   =  "[1] l, [>2] Cl [multipoles, redshift/frequency bins] || {ID} simulation {s} || seed-hi {seedhi} ||  seed-kappa-uncorr {seedk} ".format(ID=project, 
                                                                                                                                                                s=_sim_, 
                                                                                                                                                                seedhi=_dict_[_sim_]['seed_hi'],
                                                                                                                                                                seedk =_dict_[_sim_]['seed_kunc'])            
            #if len(_dict_[_sim_][_var_].shape)==1: cl_ = _dict_[_sim_][_var_][np.newaxis,:]
            #else :                                 cl_ = _dict_[_sim_][_var_]
            cl_ = np.vstack(( _dict_[_sim_]['l'] , _dict_[_sim_][_var_] ))
            pathname = cxfs.savedata_from_dict(Cl_=cl_,  filename=_var_, 
                                               path=params_path['pathout_dir'], prefix=params_path['prefix'], suffix=params_path['suffix'], 
                                               header=header, clear=clear, verbose=params_general['verbose'])
            if params_general['verbose']:
                print('{}. saving at: {}'.format(_sim_, pathname) )	
        ########
        #### ALM            
#        elif _var_.split('_')[0].lower()=='alm':
        elif 'alm' in  _var_:
            _lmax_ = _dict_[_sim_]['l'].max()
            # BO <----------------------
            if ('kappa_sim' in _var_) or ('kappa_uncorr' in _var_):
                l_ = np.array([ hp.Alm.getlm(i=il,lmax=_lmax_)[0] for il in range(_dict_[_sim_][_var_].shape[0]) ])
                m_ = np.array([ hp.Alm.getlm(i=il,lmax=_lmax_)[1] for il in range(_dict_[_sim_][_var_].shape[0]) ])
                ##########################################################################
                alm_matrix_isim_comp = np.vstack((l_,m_))
                alm_matrix_isim_real = np.vstack((l_,m_))
                alm_matrix_isim_imag = np.vstack((l_,m_))                
                alm_matrix_isim_comp = np.vstack(( alm_matrix_isim_comp, _dict_[_sim_][_var_]      ))
                alm_matrix_isim_real = np.vstack(( alm_matrix_isim_real, _dict_[_sim_][_var_].real ))
                alm_matrix_isim_imag = np.vstack(( alm_matrix_isim_imag, _dict_[_sim_][_var_].imag ))    
            else:
                l_ = np.array([ hp.Alm.getlm(i=il,lmax=_lmax_)[0] for il in range(_dict_[_sim_][_var_].shape[1]) ])
                m_ = np.array([ hp.Alm.getlm(i=il,lmax=_lmax_)[1] for il in range(_dict_[_sim_][_var_].shape[1]) ])
                ##########################################################################
                alm_matrix_isim_comp = np.vstack((l_,m_))
                alm_matrix_isim_real = np.vstack((l_,m_))
                alm_matrix_isim_imag = np.vstack((l_,m_))                
                for ibin in range(_dict_[_sim_][_var_].shape[0]):
                    alm_matrix_isim_comp = np.vstack(( alm_matrix_isim_comp, _dict_[_sim_][_var_][ibin]      ))
                    alm_matrix_isim_real = np.vstack(( alm_matrix_isim_real, _dict_[_sim_][_var_][ibin].real ))
                    alm_matrix_isim_imag = np.vstack(( alm_matrix_isim_imag, _dict_[_sim_][_var_][ibin].imag ))
        ##########################################################################
            header   =  "[1] l, [2] m, [>3] alm [multipoles, redshift/frequency bins] || {ID} simulation {s}".format(ID=project, s=_sim_)
            for i, (ialm_matrix_isim, ialm_data) in enumerate(zip([alm_matrix_isim_comp, alm_matrix_isim_real, alm_matrix_isim_imag], 
                                                                  ['complex','real', 'imag'])):
                pathname = cxfs.savedata_from_dict(Cl_=ialm_matrix_isim,  filename='_'.join(( _var_, ialm_data )), 
                                                      path=params_path['pathout_dir'],prefix=params_path['prefix'], suffix=params_path['suffix'], 
                                                      header=header, alm_data=ialm_data, clear=False, verbose=params_general['verbose'])        
                if params_general['verbose']:
                    print('{}. saving at: {}'.format(_sim_, pathname) ) 
            del alm_matrix_isim_comp, alm_matrix_isim_real, alm_matrix_isim_imag, pathname        
        ########    
        ####
        else:
           pass

pathname = '/'.join(( pathname.split('/')[:-2]  ))  
for _seed_name_ in ['seed_kunc', 'seed_hi']:
    pf = os.path.join(pathname,_seed_name_)
    if params_general['verbose']:
        print('{}. saving at: {}'.format(_seed_name_, pf) )
    np.savetxt(pf, np.array([params_sims[_seed_name_]]))#, fmt='%d')
    
print('\n------------------------------------------')
print('END')
print('------------------------------------------\n')
