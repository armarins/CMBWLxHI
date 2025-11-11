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
INI         = "generating_APS_FGremoval.ini"
name_params = os.path.join(PATH,INI)
config.read(name_params)
#General
verbose          = config.getboolean("General","verbose")
project          = config.get(       "General","project")
#load_simulations = config.getboolean("General","load_simulations")
#Simulations
nrealizations   = config.getint(    "Simulations","nrealizations")
lmax            = config.getint(    "Simulations","lmax")
apply_mask      = config.getboolean("Simulations","apply_mask" )

#Algorithms
method     = config.get(   "Algorithm","method")
wtransform = config.get(   "Algorithm","wtransform")
ns         = config.getint("Algorithm","ns")

#PATH
dirpath_chisel       = config.get("PATH","dirpath_chisel")
dirpath_sims         = config.get("PATH","dirpath_sims" )
dirpath_foregrounds  = config.get("PATH","dirpath_foregrounds" )
dirpath_mask         = config.get("PATH","dirpath_mask" )
filename_foregrounds = config.get("PATH","filename_foregrounds" )
filename_mask        = config.get("PATH","filename_mask" )
#
filepath_field1       = config.get("PATH","filepath_field1")
filepath_field2       = config.get("PATH","filepath_field2")
filepath_cross        = config.get("PATH","filepath_cross" )

#FITS SAVING
dirpath_savedata = config.get("FITS","dirpath_savedata" )
filename_newdir  = config.get("FITS","filename_newdir" )

#FITS
nside  = config.getint("FITS","nside")
numin  = config.getint("FITS","numin")
numax  = config.getint("FITS","numax")
nbands = config.getint("FITS","nbands")
unit   = config.get(   "FITS","unit")
stokes = config.get(   "FITS","stokes")


###############################################################################
# You can modify any options in the parameters.ini file by the command terminal
###############################################################################
parser = argparse.ArgumentParser(description='Modify by the command terminal parameters in {} file'.format(INI))
#General
parser.add_argument('--verbose'         , action = 'store', dest = 'verbose'          , default = verbose          , help = '')
parser.add_argument('--project'         , action = 'store', dest = 'project'          , default = project          , help = '')
#parser.add_argument('--load_simulations', action = 'store', dest = 'load_simulations' , default = load_simulations , help = '')

#Simulations
parser.add_argument('--nrealizations', action = 'store', dest = 'nrealizations', default = nrealizations, help = '')
parser.add_argument('--lmax'         , action = 'store', dest = 'lmax'         , default = lmax         , help = '')
parser.add_argument('--apply_mask'   , action = 'store', dest = 'apply_mask'   , default = apply_mask   , help = '')

#Algorithms
parser.add_argument('--method'    , action = 'store', dest = 'method'    , default = method    , help = '')
parser.add_argument('--wtransform', action = 'store', dest = 'wtransform', default = wtransform, help = '')
parser.add_argument('--ns'        , action = 'store', dest = 'ns'        , default = ns        , help = '')

#PATH
parser.add_argument('--dirpath_chisel'      , action = 'store', dest = 'dirpath_chisel'      , default = dirpath_chisel      , help = '')
parser.add_argument('--dirpath_sims'        , action = 'store', dest = 'dirpath_sims'        , default = dirpath_sims        , help = '')
parser.add_argument('--dirpath_foregrounds' , action = 'store', dest = 'dirpath_foregrounds' , default = dirpath_foregrounds , help = '')
parser.add_argument('--dirpath_mask'        , action = 'store', dest = 'dirpath_mask'        , default = dirpath_mask        , help = '')
parser.add_argument('--filename_foregrounds', action = 'store', dest = 'filename_foregrounds', default = filename_foregrounds, help = '')
parser.add_argument('--filename_mask'       , action = 'store', dest = 'filename_mask'       , default = filename_mask       , help = '')

parser.add_argument('--filepath_field1'     , action = 'store', dest = 'filepath_field1'     , default = filepath_field1     , help = '')
parser.add_argument('--filepath_field2'     , action = 'store', dest = 'filepath_field2'     , default = filepath_field2     , help = '')
parser.add_argument('--filepath_cross'      , action = 'store', dest = 'filepath_cross'      , default = filepath_cross      , help = '')

#FITS
parser.add_argument('--dirpath_savedata', action = 'store', dest = 'dirpath_savedata', default = dirpath_savedata, help = '')
parser.add_argument('--filename_newdir' , action = 'store', dest = 'filename_newdir' , default = filename_newdir , help = '')
 
parser.add_argument('--nside'     , action = 'store', dest = 'nside'     , default = nside     , help = '')
parser.add_argument('--numin'     , action = 'store', dest = 'numin'     , default = numin     , help = '')
parser.add_argument('--numax'     , action = 'store', dest = 'numax'     , default = numax     , help = '')
parser.add_argument('--nbands'    , action = 'store', dest = 'nbands'    , default = nbands    , help = '')
parser.add_argument('--unit'      , action = 'store', dest = 'unit'      , default = unit      , help = '')
parser.add_argument('--stokes'    , action = 'store', dest = 'stokes'    , default = stokes    , help = '')

arguments = parser.parse_args()
###############################################################################
# Variables
###############################################################################
#General
verbose          = bool(arguments.verbose)
project          = str( arguments.project)
#load_simulations = bool(arguments.load_simulations)
#Simulations
nrealizations  = int( arguments.nrealizations)
lmax           = int( arguments.lmax)
apply_mask     = bool(arguments.apply_mask)
#Algorithms
method     = str(arguments.method)
wtransform = str(arguments.wtransform)
ns         = int( arguments.ns)
#PATH
dirpath_chisel       = str(arguments.dirpath_chisel)
dirpath_sims         = str(arguments.dirpath_sims)
dirpath_foregrounds  = str(arguments.dirpath_foregrounds)
dirpath_mask         = str(arguments.dirpath_mask)
filename_foregrounds = str(arguments.filename_foregrounds)
filename_mask        = str(arguments.filename_mask)

filepath_field1      = str(arguments.filepath_field1)
filepath_field2      = str(arguments.filepath_field2)
filepath_cross       = str(arguments.filepath_cross)

#FITS
dirpath_savedata = str(arguments.dirpath_savedata)
filename_newdir  = str(arguments.filename_newdir)

nside      = int(arguments.nside)
numin      = float(arguments.numin)
numax      = float(arguments.numax)
nbands     = int(arguments.nbands)
unit       = str(arguments.unit)
stokes     = str(arguments.stokes)


dirpath_savedata = os.path.join(dirpath_savedata, project)
#####################################################################################################################################################################
#####################################################################################################################################################################
params_general = pd.Series({'verbose':verbose,
                            'project':project
                           })
                            #'load_simulations':load_simulations})

params_sims = pd.Series({'nreals'    :nrealizations,
                         'lmax'      :lmax,
                         'apply_mask':apply_mask
                       })
               
params_chisel = pd.Series({'method'    :method,
                           'wtransform':wtransform,
                           'ns'        :ns
                         })            

params_path = pd.Series({'dirpath_chisel'      : dirpath_chisel,
                         'dirpath_sims'        : dirpath_sims,
                         'dirpath_foregrounds' : dirpath_foregrounds,
                         'dirpath_mask'        : dirpath_mask,
                         'dirpath_savedata'    : dirpath_savedata,                         
                         'filename_foregrounds': filename_foregrounds,                
                         'filename_mask'       : filename_mask,
                         'filename_newdir'     : filename_newdir,
                         'filepath_field1'     : filepath_field1,
                         'filepath_field2'     : filepath_field2,
                         'filepath_cross'      : filepath_cross,                         
                        })

params_fits = pd.Series({'nside'     :nside,
                         'numin'     :numin,
                         'numax'     :numax,
                         'nbands'    :nbands,
                         'unit'      :unit,
                         'stokes'    :stokes
                       })

#del verbose, project#, load_simulations
del nrealizations, lmax, apply_mask
del method, wtransform, ns
del dirpath_chisel, dirpath_sims, dirpath_foregrounds, dirpath_mask, filename_foregrounds, filename_mask, dirpath_savedata,filename_newdir,filepath_field1,filepath_field2,filepath_cross
del nside, numin, numax, nbands, unit, stokes

#####################################################################################################################################################################
#####################################################################################################################################################################
if params_general['verbose']:
    print('\n--Variables to be used:--')
    print(params_general)    
    print(params_sims)    
    print(params_path)
    print(params_fits)
    print()

sys.path.insert(1, params_path['dirpath_chisel'])    
import Extension4BINGO   as cs
import load_standard_params as loadparams
#import pyMRS as pymrs
#import statcosmo as statc
params_maps    = loadparams.load_params_maps()
params_CS      = loadparams.load_params_CS()
params_WT      = loadparams.load_params_WT()
params_path_cs = loadparams.load_params_path('amarins')
####
#modifying chisel params
del params_maps['iseed']    
params_CS['method']     = dcopy(params_chisel['method'])
params_CS['wtransform'] = dcopy(params_chisel['wtransform'])
params_CS['ns']         = dcopy(params_chisel['ns'])
####
if params_general['verbose']:
    print('--CHISEL infos:--')
    print(params_maps)
    print(params_CS)
    print(params_WT)

#####################################################################################################################################################################
#####################################################################################################################################################################
#loading Foregrodun (FG) maps and the mask
Mfg = hdata.getmap(dirpath_=params_path.dirpath_foregrounds, filename_=params_path.filename_foregrounds, healpix_readingformat=False, hdu=1)
if (params_sims.apply_mask) and (params_path.dirpath_mask!="") and (params_path.filename_mask!=""):
    MASK = hdata.getmap(dirpath_=params_path.dirpath_mask, filename_=params_path.filename_mask, healpix_readingformat=False, hdu=1)
    params_chisel['coverage_flag'] = 'masked'
    if params_path['filename_newdir'] =="":
        params_path['filename_newdir'] = 'masked'    
else:
    MASK = np.ones_like(Mfg[0,:])
    params_chisel['coverage_flag'] = 'fullsky'
    if params_path['filename_newdir'] =="":
        params_path['filename_newdir'] = 'fullsky'


hdata.file_verification( params_path['dirpath_savedata'], '', params_path['filename_newdir'] )
params_path['dirpath_savedata'] = os.path.join(params_path['dirpath_savedata'], params_path['filename_newdir'] )


hdata.file_verification(params_path['dirpath_savedata'],'','ns{}'.format(params_CS['ns']))
params_path['dirpath_savedata_ns'] = os.path.join(params_path['dirpath_savedata'],'ns{}'.format(params_CS['ns']))
#####################################################################################################################################################################
#####################################################################################################################################################################
output_info = {'field':'field', 'nside':params_fits.nside, 
               'frequency'   :{'min':params_fits.numin, 'max':params_fits.numax,'nbands':params_fits.nbands,  'unit':params_fits.unit}, 
               'stokes'      :params_fits['stokes'],
               'coverage'    :params_chisel['coverage_flag'], 
               'pathout'     :params_path['dirpath_savedata_ns'],
               'apply_beam'  :False, 
               'namePL'      :False, 
               'nameFG'      :params_path.filename_foregrounds, 
               'apply_mask'  :params_sims.apply_mask, 
               'nameM'       :params_path.filename_mask,               
               'pathdir_N'   :False
              }

output_info['filename'] = hdata.new_formattingnames(F=output_info['field'], 
                                                    S=output_info['stokes'], 
                                                    NSIDE_=output_info['nside'], 
                                                    freq_min=int(output_info['frequency']['min']), 
                                                    freq_max=int(output_info['frequency']['max']),
                                                    freq_unit=output_info['frequency']['unit'],
                                                    Num_=output_info['frequency']['nbands'], 
                                                    C=output_info['coverage'], 
                                                    Bres=None,Bmodel=None,A=None, R='')
#####################################################################################################################################################################
#####################################################################################################################################################################                                                    
if 1:
    print('\n------------------------------------------------------------------')
    print('Loading the APS theoretically simulated')
    print('------------------------------------------------------------------\n')

clf1 = np.loadtxt(params_path['filepath_field1']).T
clf2 = np.loadtxt(params_path['filepath_field2']).T
clcx = np.loadtxt(params_path['filepath_cross'] ).T
    
#clf1_dict = cxft.cls_from_matrix2dict(clf1)
#clf2_dict = cxft.cls_from_matrix2dict(clf2)
#clcx_dict = cxft.cls_from_matrix2dict(clcx)   
#nch = int(len(clf1_dict.keys())-1)
#del clf1, clf2, clcx        

load_sims=1
if load_sims:
    itime         = time.time()
    varnames      = ['alm_hi_sim']#, 'cl_hi_sim', 'cl_kappa_sim', 'cl_cross_sim','cross_correlation_coef', 'alm_kappa_sim']
    dict_all      = cxfs.get_file_simulated_data(params_path['dirpath_sims'], varnames, nsims_to_use=params_sims.nreals)
    dict_all_mean = cxfs.dict_averages_from_loaded_data(dict_all_=dict_all, varnames_=varnames)
    nsims         = int(dict_all_mean[varnames[0]]['bin1']['matrix'].shape[0]    )
    if params_general['verbose']: print('Loading time for {0:d} sims: {1:.2f} seg'.format(nsims, time.time()-itime)) 
    del nsims
#print(nsims)
params_sims['realizations'] = np.asarray(list(dict_all.keys()))#[:nreals]
if params_general['verbose']: print("Using the sims: ",params_sims['realizations'])
params_sims['realizations']
#####################################################################################################################################################################
#####################################################################################################################################################################                                                    
if 1:
    print('\n------------------------------------------------------------------')
    print('Starting the Foreground Removal Process')
    print('------------------------------------------------------------------\n')

for j, jsim in enumerate(params_sims['realizations']):
    print('Job {} -- Sim: {}'.format(j+1, jsim))
    lmax        = dict_all[jsim]['alm_hi_sim'][0].real.astype(int).max()
    mmax        = dict_all[jsim]['alm_hi_sim'][1].real.astype(int).max()
    alm_hi_jsim = dict_all[jsim]['alm_hi_sim'][2:]
    for num in range(alm_hi_jsim.shape[0]):
        alm_hi_num = np.ascontiguousarray(alm_hi_jsim[num,:])
        if not num:  Mhi = hp.alm2map( alm_hi_num,                  params_fits.nside, pol=False)    
        else:        Mhi = np.vstack(( Mhi, hp.alm2map( alm_hi_num, params_fits.nside, pol=False) ))       
    ###
    params_maps["getdata"] = "observed"
    #subdirs = cs.checkdir(params_path_cs.pathout, subdirs=["21cm","foregrounds","mixmatrix"],restart=True)
    params_cs, params_wt = cs.load(params_CS,params_WT)
    X = dcopy(Mhi + Mfg)*MASK
    X = cs.adaptation_maps(X, params_maps, params_path)
    Xr = cs.maps2CSmaps(    X, params_wt  , params_cs)
    ######
    jseed = jsim.split('sim')[1]
    output_info['pathout_dir'] = os.path.join(output_info['pathout'],str(jsim))
    #####
    for i, (ifield_name, ifield_data) in enumerate(zip(['extHI', 'extFG','mixmatrix'],
                                                       [Xr['reconstruction']['21cm'], Xr['reconstruction']['foregrounds'], Xr["mixmatrix"] ]
                                                      )
                                                   ):
        output_info['field']        = ifield_name
        output_info['irealization'] = cxft.get_Lformat(L_=jseed)    
        jfilename = dcopy(output_info['filename'].replace('L',output_info['irealization']))
        jfilename = jfilename.replace('field', output_info['field'])
        output_info['pathout_field'] = os.path.join(output_info['pathout_dir'], jfilename)    
        hdata.file_verification(output_info['pathout'],output_info['filename'],str(jsim))
        cxft.savefits(data_=ifield_data, params_=output_info, save_mixmatrix=True, S_mixmatrix=Xr['mixmatrix'])
    
    print()


if 1:
    print('\n------------------------------------------')
    print('END')
    print('------------------------------------------\n')
