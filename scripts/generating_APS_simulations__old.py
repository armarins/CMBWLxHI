菜鸟

`#sys.path.insert(1, os.getcwd())
import os, sys, time
from copy import deepcopy as dcopy
import cross_functions_theory as cxft
import cross_functions_simulations as cxfs
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
INI         = "generating_APS_simulations.ini"
name_params = os.path.join(os.getcwd(),INI)
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

############################################################################################
#####################################################################################################
if prefix=="":
    prefix=None
if suffix=="":
    suffix=None

#print(generate_cl_field_1, type(generate_cl_field_1))
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
               'seed0'          :seed0
               }

del verbose, filepath_field1, filepath_field2, filepath_cross, filepath_foregrounds, filepath_mask, pathout, prefix, suffix
del add_foregrounds, apply_mask, nrealizations, seed0
#####################################################################################################################################################################
#####################################################################################################################################################################


#####################################################################################################
#####################################################################################################
#if params_general['verbose']:
if 1:
    print('\n------------------------------------------------------------------')
    print('Starting to produce Angular Power Spectrum theoretical simulations')
    print('------------------------------------------------------------------\n')


if params_general['verbose']:
    print('Loading data from {}'.format(params_path['filepath_field1']))
clf1 = np.loadtxt(params_path['filepath_field1']).T
if params_general['verbose']:
    print('Loading data from {}'.format(params_path['filepath_field2']))
clf2 = np.loadtxt(params_path['filepath_field2']).T
if params_general['verbose']:
    print('Loading data from {}'.format(params_path['filepath_cross']))
clcx = np.loadtxt(params_path['filepath_cross']).T

clf1_dict = cxft.cls_from_matrix2dict(clf1)
clf2_dict = cxft.cls_from_matrix2dict(clf2)
clcx_dict = cxft.cls_from_matrix2dict(clcx)


#################################
#####################################################################################################
#####################################################################################################
#jseeds1  = np.arange(params_sims['seed0']                                   , params_sims['seed0'] +   params_sims['nrealizations']  , 1, dtype=np.int16) 
#jseeds2  = np.arange(params_sims['seed0'] + params_sims['nrealizations'] + 1, params_sims['seed0'] + 2*params_sims['nrealizations']+1, 1, dtype=np.int16)
params_sims['seed_kunc'] = np.arange(params_sims['seed0']                                   , params_sims['seed0'] +   params_sims['nrealizations']  , 1, dtype=np.int16) 
params_sims['seed_hi'  ] = np.arange(params_sims['seed0'] + params_sims['nrealizations'] + 1, params_sims['seed0'] + 2*params_sims['nrealizations']+1, 1, dtype=np.int16)

if params_general['verbose']:
    print('Generating simulations...')
_dict_= cxfs.get_dict_simulations(nreals=params_sims['nrealizations'], first_seed=params_sims['seed0'], jseed_kunc=params_sims['seed_kunc'] , jseed_hi=params_sims['seed_hi'],
                                  Chi_dict=clf1_dict, Ckk_vec=clf2_dict['bin1'], Cross_dict=clcx_dict, 
                                  l_          =clf1_dict['l'], 
                                  nbins_      =len(clcx_dict.keys())-1,
                                  fact        =1,
                                  show_time   =False, 
                                  sort_bin_sim=False)

# LOOP        ###################
#################################
for i,_sim_ in enumerate(_dict_.keys()):   
    params_path['pathout_dir'] = os.path.join( params_path['pathout'], '{ID}/{s}'.format(ID=project, s=_sim_) )
    for j,_var_ in enumerate(_dict_['sim0']['bin1'].keys()):    
        if not i+j: 
            clear=True
        else:
            clear=False        
        #######
        #### CL
#        if _var_.split('_')[0].lower() in ['cl', 'cross']:
        if ('cl' in  _var_) or ('cross' in  _var_):
            cl_matrix_isim = _dict_[_sim_]['bin1']['l']
            for k,_bin_ in enumerate(_dict_['sim0'].keys()): #only to know the bin names
                cl_matrix_isim = np.vstack(( cl_matrix_isim, _dict_[_sim_][_bin_][_var_] ))
            header   =  "[1] l, [>2] Cl [multipoles, redshift/frequency bins] || {ID} simulation {s}".format(ID=project, s=_sim_)
            header   =  "[1] l, [>2] Cl [multipoles, redshift/frequency bins] || {ID} simulation {s} || seed-hi {seedhi} ||  seed-kappa-uncorr {seedk} ".format(ID=project, s=_sim_, 
                                                                                                                                            seedhi=_dict_[_sim_]['bin1']['seed-hi'],
                                                                                                                                            seedk =_dict_[_sim_]['bin1']['seed-kunc'])            
            pathname = cxfs.savedata_from_dict(Cl_=cl_matrix_isim,  filename=_var_, 
                                               path=params_path['pathout_dir'], prefix=params_path['prefix'], suffix=params_path['suffix'], 
                                               header=header, clear=clear, verbose=params_general['verbose'])
            if params_general['verbose']:
                print('{}. saving at: {}'.format(_sim_, pathname) )	
            del cl_matrix_isim
        ########
        #### ALM            
#        elif _var_.split('_')[0].lower()=='alm':
        elif 'alm' in  _var_:
            _bin_  = 'bin1' #it`s only fof extracting the
            _lmax_ = _dict_[_sim_][_bin_]['l'].max()
            l_ = np.array([ hp.Alm.getlm(i=il,lmax=_lmax_)[0] for il in range(_dict_[_sim_][_bin_][_var_].shape[0]) ])
            m_ = np.array([ hp.Alm.getlm(i=il,lmax=_lmax_)[1] for il in range(_dict_[_sim_][_bin_][_var_].shape[0]) ])
        ##########################################################################
            alm_matrix_isim_comp = np.vstack((l_,m_))
            alm_matrix_isim_real = np.vstack((l_,m_))
            alm_matrix_isim_imag = np.vstack((l_,m_))
        
            for i,ibin in enumerate(_dict_[_sim_].keys()):
                alm_matrix_isim_comp = np.vstack(( alm_matrix_isim_comp, _dict_[_sim_][ibin][_var_]      ))
                alm_matrix_isim_real = np.vstack(( alm_matrix_isim_real, _dict_[_sim_][ibin][_var_].real ))
                alm_matrix_isim_imag = np.vstack(( alm_matrix_isim_imag, _dict_[_sim_][ibin][_var_].imag ))
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
#p = '/home/user/Programmation/cross/data/simulations/bingo0/sim0/cl_kappa_sim.txt'
pathname = '/'.join(( pathname.split('/')[:-2]  ))
#print(os.path.join(p,'seed_kunc.txt'))
#print(os.path.join(p,'seed_hi.txt'))          
for _seed_name_ in ['seed_kunc', 'seed_hi']:
    pf = os.path.join(pathname,_seed_name_)
    if params_general['verbose']:
        print('{}. saving at: {}'.format(_seed_name_, pf) )
    np.savetxt(pf, params_sims[_seed_name_], fmt='%d')

if 1:
    print('\n------------------------------------------')
    print('END')
    print('------------------------------------------\n')
