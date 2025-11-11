#sys.path.insert(1, os.getcwd())
import os, sys, time
from copy import deepcopy     as dcopy
import cross_functions_theory as cxft
import cross_functions_simulations as cxfs
import handling_data   as hdata
import noise_functions as fnoise
import numpy  as np
import healpy as hp
import pandas as pd
import json
import argparse
import warnings
warnings.filterwarnings("ignore")
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

###################################################################
# This part is for extracting information from parameters.ini file
###################################################################
timei       = time.time()
#PATH='/data/AMARINS/LSSxHI-CODES/scripts'
PATH        = os.path.dirname(os.path.abspath(__file__)) #Assuming it is in the same directory that this py
INI         = "generating_APS_TransferFunction.ini"
name_params = os.path.join(PATH,INI)
#
if not os.path.exists(name_params):
        print(f"Error: Configuration file not found at {name_params}")
        sys.exit(1)
#
config.read(name_params)
#General
verbose          = config.getboolean("General","verbose")
project          = config.get(       "General","project")
threshold        = config.get(       "General","threshold")
#load_simulations = config.getboolean("General","load_simulations")
#Simulations
nrealizations   = config.getint(    "Simulations","nrealizations")
nside           = config.getint(    "Simulations","nside")
lmin            = config.getint(    "Simulations","lmin")
lmax            = config.getint(    "Simulations","lmax")
apply_mask      = config.getboolean("Simulations","apply_mask" )
#apply_namaster  = config.getboolean("Simulations","apply_namaster" )
#namaster_delell = config.getint(    "Simulations","namaster_delell")
include_wn      = config.getboolean("Simulations","include_wn" )
include_beam    = config.getboolean("Simulations","include_beam" )
wn_sigma_pix    = config.getfloat(  "Simulations","wn_sigma_pix" )
beam_fwhm       = config.getfloat(  "Simulations","beam_fwhm" )
#
numin           = config.getfloat(  "Simulations","numin"  )
numax           = config.getfloat(  "Simulations","numax"  )
nbands          = config.getint(    "Simulations","nbands" )
stokes          = config.get(       "Simulations","stokes" )
unit            = config.get(       "Simulations","unit"   )

#TF
tf_hixhi      = config.getboolean("TransferFunction","tf_hixhi"     )
tf_gxhi       = config.getboolean("TransferFunction","tf_gxhi"      )
cross_g_bin   = config.getint(    "TransferFunction","cross_g_bin"  )
auto_hi_bins  = config.get(       "TransferFunction","auto_hi_bins" )
auto_hi_bins  = np.array([float(item.strip()) for item in auto_hi_bins.split(',')])
cross_hi_bins = config.get(       "TransferFunction","cross_hi_bins")
cross_hi_bins = np.array([float(item.strip()) for item in cross_hi_bins.split(',')])
#
#save_hi_plus_fg                  = config.getboolean("TransferFunction","save_hi_plus_fg"                 )
#save_hi_plus_fg_plus_wn          = config.getboolean("TransferFunction","save_hi_plus_fg_plus_wn"         )
#save_hi_plus_fg_smoothed         = config.getboolean("TransferFunction","save_hi_plus_fg_smoothed"        )
#save_hi_plus_fg_plus_wn_smoothed = config.getboolean("TransferFunction","save_hi_plus_fg_plus_wn_smoothed")
save_simulated   = config.getboolean("TransferFunction","save_simulated"   )
save_theoretical = config.getboolean("TransferFunction","save_theoretical" )
save_hdf5        = config.getboolean("TransferFunction","save_hdf5"        )
save_txt         = config.getboolean("TransferFunction","save_txt"         )


#Algorithms
method     = config.get(   "Algorithm","method")
wtransform = config.get(   "Algorithm","wtransform")
ns         = config.getint("Algorithm","ns")

#Inputs
dirpath_chisel       = config.get("Inputs","dirpath_chisel")
dirpath_sims         = config.get("Inputs","dirpath_sims" )
dirpath_foregrounds  = config.get("Inputs","dirpath_foregrounds" )
dirpath_mask         = config.get("Inputs","dirpath_mask" )
filename_foregrounds = config.get("Inputs","filename_foregrounds" )
filename_mask        = config.get("Inputs","filename_mask" )
#
filepath_field1      = config.get("Inputs","filepath_field1")
filepath_field2      = config.get("Inputs","filepath_field2")
filepath_cross       = config.get("Inputs","filepath_cross" )
#
#Outputs
pathout              = config.get(    "Outputs","pathout"         )
###############################################################################
# You can modify any options in the parameters.ini file by the command terminal
###############################################################################
parser = argparse.ArgumentParser(description='Modify by the command terminal parameters in {} file'.format(INI))
#General
parser.add_argument('--verbose'         , type=hdata.str2bool, dest = 'verbose'  , default = verbose   , help = 'Enable verbose output (True/False).')
parser.add_argument('--project'         , action = 'store'   , dest = 'project'  , default = project   , help = 'Project name to be used as directory name.')
parser.add_argument('--threshold'       , action = 'store'   , dest = 'threshold', default = threshold , help = 'Threshold value for Cls.')
#parser.add_argument('--load_simulations', action = 'store', dest = 'load_simulations' , default = load_simulations , help = '')

#Simulations
parser.add_argument('--nrealizations'   , action = 'store'    , dest = 'nrealizations', default = nrealizations,  help = 'Number of realizations to estimate the transfer function.')
parser.add_argument('--nside'           , action = 'store'    , dest = 'nside'        , default = nside        , help = 'Nside for healpy/Healpix maps.')
parser.add_argument('--lmin'            , action = 'store'    , dest = 'lmin'         , default = lmin         , help = 'Minimum multipole.')
parser.add_argument('--lmax'            , action = 'store'    , dest = 'lmax'         , default = lmax         , help = 'Maximum multipole.')
parser.add_argument('--apply_mask'      , type=hdata.str2bool , dest = 'apply_mask'   , default = apply_mask   , help = 'Apply mask (True/False).')
parser.add_argument('--include_wn'      , type=hdata.str2bool , dest = 'include_wn'   , default = include_wn   , help = 'Include white noise (True/False).')
parser.add_argument('--include_beam'    , type=hdata.str2bool , dest = 'include_beam' , default = include_beam , help = 'Include beam effect (True/False).')
parser.add_argument('--wn_sigma_pix'    , action = 'store'    , dest = 'wn_sigma_pix' , default = wn_sigma_pix , help = 'White noise sigma per pixel.')
parser.add_argument('--beam_fwhm'       , action = 'store'    , dest = 'beam_fwhm'    , default = beam_fwhm    , help = 'Beam FWHM in arcmin.')
parser.add_argument('--numin'           , action = 'store'    , dest = 'numin'        , default = numin        , help = 'Minimum frequency.')
parser.add_argument('--numax'           , action = 'store'    , dest = 'numax'        , default = numax        , help = 'Maximum frequency.')
parser.add_argument('--nbands'          , action = 'store'    , dest = 'nbands'       , default = nbands       , help = 'Number of frequency bands.')
parser.add_argument('--unit'            , action = 'store'    , dest = 'unit'         , default = unit         , help = 'Map units.')
parser.add_argument('--stokes'          , action = 'store'    , dest = 'stokes'       , default = stokes       , help = 'Stokes parameter.')
#parser.add_argument('--apply_mask'     , action = 'store', dest = 'apply_mask'     , default = apply_mask     , help = '')
#parser.add_argument('--include_wn'     , action = 'store', dest = 'include_wn'     , default = include_wn     , help = '')
#parser.add_argument('--include_beam'   , action = 'store', dest = 'include_beam'   , default = include_beam   , help = '')
#parser.add_argument('--apply_namaster' , action = 'store', dest = 'apply_namaster' , default = apply_namaster , help = '')
#parser.add_argument('--namaster_delell', action = 'store', dest = 'namaster_delell', default = namaster_delell, help = '')
#TransferFunction
parser.add_argument('--tf_hixhi'        , type=hdata.str2bool, dest = 'tf_hixhi'        , default=tf_hixhi        , help = 'Calculate HIxHI transfer function (True/False).')
parser.add_argument('--tf_gxhi'         , type=hdata.str2bool, dest = 'tf_gxhi'         , default=tf_gxhi         , help = 'Calculate GxHI transfer function (True/False).' )
parser.add_argument('--cross_g_bin'     , action = 'store'   , dest = 'cross_g_bin'     , default = cross_g_bin   , help = 'Galaxy bin for cross-correlation.')
parser.add_argument('--auto_hi_bins'    , action = 'store'   , dest = 'auto_hi_bins'    , default = auto_hi_bins  , help = 'HI bins for auto-correlation.')
parser.add_argument('--cross_hi_bins'   , action = 'store'   , dest = 'cross_hi_bins'   , default = cross_hi_bins , help = 'HI bins for cross-correlation.')
parser.add_argument('--save_simulated'  , type=hdata.str2bool, dest = 'save_simulated'  , default=save_simulated  , help = 'Save simulated data (True/False).')
parser.add_argument('--save_theoretical', type=hdata.str2bool, dest = 'save_theoretical', default=save_theoretical, help = 'Save theoretical data (True/False).')
parser.add_argument('--save_hdf5'       , type=hdata.str2bool, dest = 'save_hdf5'       , default=save_hdf5       , help = 'Save output as HDF5 (True/False).')
parser.add_argument('--save_txt'        , type=hdata.str2bool, dest = 'save_txt'        , default=save_txt        , help = 'Save output as TXT (True/False).')
#parser.add_argument('--tf_hixhi'        , action = 'store', dest = 'tf_hixhi'        , default = tf_hixhi        , help = '')
#parser.add_argument('--tf_gxhi'         , action = 'store', dest = 'tf_gxhi'         , default = tf_gxhi         , help = '')
#parser.add_argument('--save_hi_plus_fg'                 , action = 'store', dest = 'save_hi_plus_fg'                 , default = save_hi_plus_fg                 , help = '')
#parser.add_argument('--save_hi_plus_fg_smoothed'        , action = 'store', dest = 'save_hi_plus_fg_smoothed'        , default = save_hi_plus_fg_smoothed        , help = '')
#parser.add_argument('--save_hi_plus_fg_plus_wn'         , action = 'store', dest = 'save_hi_plus_fg_plus_wn'         , default = save_hi_plus_fg_plus_wn         , help = '')
#parser.add_argument('--save_hi_plus_fg_plus_wn_smoothed', action = 'store', dest = 'save_hi_plus_fg_plus_wn_smoothed', default = save_hi_plus_fg_plus_wn_smoothed, help = '')
#parser.add_argument('--save_simulated'  , action = 'store', dest = 'save_simulated'  , default = save_simulated  , help = '')
#parser.add_argument('--save_theoretical', action = 'store', dest = 'save_theoretical', default = save_theoretical, help = '')
#parser.add_argument('--save_hdf5' , action = 'store', dest = 'save_hdf5', default = save_hdf5 , help = '')
#parser.add_argument('--save_txt'  , action = 'store', dest = 'save_txt' , default = save_txt  , help = '')

#Algorithms
parser.add_argument('--method'    , action = 'store', dest = 'method'    , default = method    , help = 'Blind Foreground Removal Algorithm.')
parser.add_argument('--wtransform', action = 'store', dest = 'wtransform', default = wtransform, help = 'Wavelet transform type.')
parser.add_argument('--ns'        , action = 'store', dest = 'ns'        , default = ns        , help = 'Number of independent templates.')

#IN/OUTPUT
parser.add_argument('--dirpath_chisel'      , action = 'store', dest = 'dirpath_chisel'      , default = dirpath_chisel      , help = '')
parser.add_argument('--dirpath_sims'        , action = 'store', dest = 'dirpath_sims'        , default = dirpath_sims        , help = '')
parser.add_argument('--dirpath_foregrounds' , action = 'store', dest = 'dirpath_foregrounds' , default = dirpath_foregrounds , help = '')
parser.add_argument('--dirpath_mapython generating_APS_TransferFunction.py --verbose yes --cross_g_bin 1 --nrealizations 3sk'        , action = 'store', dest = 'dirpath_mask'        , default = dirpath_mask        , help = '')
parser.add_argument('--filename_foregrounds', action = 'store', dest = 'filename_foregrounds', default = filename_foregrounds, help = '')
parser.add_argument('--filename_mask'       , action = 'store', dest = 'filename_mask'       , default = filename_mask       , help = '')
parser.add_argument('--filepath_field1'     , action = 'store', dest = 'filepath_field1'     , default = filepath_field1     , help = '')
parser.add_argument('--filepath_field2'     , action = 'store', dest = 'filepath_field2'     , default = filepath_field2     , help = '')
parser.add_argument('--filepath_cross'      , action = 'store', dest = 'filepath_cross'      , default = filepath_cross      , help = '')
parser.add_argument('--pathout'             , action = 'store', dest = 'pathout'             , default = pathout             , help = '')

arguments = parser.parse_args()
###############################################################################
# Variables
###############################################################################
#General
verbose   = bool(arguments.verbose   )
project   = str( arguments.project   )
threshold = float(arguments.threshold)
#load_simulations = bool(arguments.load_simulations)

#Simulations
nrealizations   = int( arguments.nrealizations )
nside           = int( arguments.nside         )
lmin            = int( arguments.lmin          )
lmax            = int( arguments.lmax          )
apply_mask      = bool(arguments.apply_mask    )
#apply_namaster  = bool(arguments.apply_namaster)
#namaster_delell = int(arguments.namaster_delell)
include_wn      = bool(arguments.include_wn    )
include_beam    = bool(arguments.include_beam  )
wn_sigma_pix    = float(arguments.wn_sigma_pix )
beam_fwhm       = float(arguments.beam_fwhm    )
#
numin            = float(arguments.numin)
numax            = float(arguments.numax)
nbands           = int(arguments.nbands )
stokes           = str(arguments.stokes )
unit             = str(arguments.unit   )

#TF
tf_hixhi         = bool(arguments.tf_hixhi)
tf_gxhi          = bool(arguments.tf_gxhi )
#auto_g_bin      = np.asarray(arguments.cross_g_bin  , dtype=np.int32)
cross_g_bin      = np.asarray(arguments.cross_g_bin  , dtype=np.int32)
auto_hi_bins     = np.asarray(arguments.auto_hi_bins , dtype=np.int32)
cross_hi_bins    = np.asarray(arguments.cross_hi_bins, dtype=np.int32)
#save_hi_plus_fg                   = bool(arguments.save_hi_plus_fg )
#save_hi_plus_fg_smoothed          = bool(arguments.save_hi_plus_fg_smoothed )
#save_hi_plus_fg_plus_wn           = bool(arguments.save_hi_plus_fg_plus_wn )
#save_hi_plus_fg_plus_wn_smoothed  = bool(arguments.save_hi_plus_fg_plus_wn_smoothed )
save_simulated   = bool(arguments.save_simulated )
save_theoretical = bool(arguments.save_theoretical )
save_hdf5        = bool(arguments.save_hdf5 )
save_txt         = bool(arguments.save_txt )

#Algorithms
method     = str(arguments.method    )
wtransform = str(arguments.wtransform)
ns         = int( arguments.ns       )

#IN/OUTPUT
dirpath_chisel       = str(arguments.dirpath_chisel)
dirpath_sims         = str(arguments.dirpath_sims)
dirpath_foregrounds  = str(arguments.dirpath_foregrounds)
dirpath_mask         = str(arguments.dirpath_mask)
filename_foregrounds = str(arguments.filename_foregrounds)
filename_mask        = str(arguments.filename_mask)
filepath_field1      = str(arguments.filepath_field1)
filepath_field2      = str(arguments.filepath_field2)
filepath_cross       = str(arguments.filepath_cross)
pathout              = str(arguments.pathout )


dirpath_savedata = os.path.join(pathout, project)
#####################################################################################################################################################################
if filename_mask=="":
    filename_mask=None
lell  = np.arange(int(3*nside))
fell  = lell*(lell+1)/2/np.pi           
if 0:#apply_namaster:   #<-- not implemented yet
    import pymaster as nmt
    b_nmt = nmt.NmtBin.from_nside_linear(nside, nlb=namaster_delell)     
    leff  = b_nmt.get_effective_ells()
    feff  = leff*(leff+1)/2/np.pi    
else:   #<-- not implemented yet
    b_nmt = None
    leff  = dcopy(lell)
    feff  = dcopy(fell)
#mask_nmt = (leff>=lmin)*(leff<=lmax)    #<-- not implemented yet
#mask_ell = (lell>=lmin)*(lell<=lmax)    #<-- not implemented yet    
#####################################################################################################################################################################
params_general = pd.Series({'verbose'  : verbose,
                            'project'  : project,
                            'threshold': threshold,
                           })

params_sims = pd.Series({'nreals'         : nrealizations,
                         'nside'          : nside,                                          
                         'lmin'           : lmin,                         
                         'lmax'           : lmax,
                         'apply_mask'     : apply_mask,
                         #'apply_namaster' : apply_namaster,
                         #'namaster_delell': namaster_delell,
                         'include_wn'     : include_wn,
                         'include_beam'   : include_beam,
                         'wn_sigma_pix'   : wn_sigma_pix,
                         'beam_fwhm'      : beam_fwhm,        
                         'numin'          : numin,
                         'numax'          : numax,
                         'nbands'         : nbands,
                         'stokes'         : stokes,
                         'unit'           : unit,
                        })
               
params_tf = pd.Series({'tf_hixhi'        : tf_hixhi,
                       'tf_gxhi'         : tf_gxhi,
                       'auto_g_bin'      : dcopy( cross_g_bin ),
                       'cross_g_bin'     : cross_g_bin ,
                       'auto_hi_bins'    : auto_hi_bins,
                       'cross_hi_bins'   : cross_hi_bins,
                       'save_simulated'  : save_simulated,
                       'save_theoretical': save_theoretical,
                       'save_hdf5'       : save_hdf5,
                       'save_txt'        : save_txt,
                       #'save_hi_plus_fg'                 :save_hi_plus_fg,
                       #'save_hi_plus_fg_smoothed'        :save_hi_plus_fg_smoothed,
                       #'save_hi_plus_fg_plus_wn'         :save_hi_plus_fg_plus_wn,
                       #'save_hi_plus_fg_plus_wn_smoothed':save_hi_plus_fg_plus_wn_smoothed,
                      })   

params_chisel = pd.Series({'method'    :method,
                           'wtransform':wtransform,
                           'ns'        :ns})            

params_inout = pd.Series({'dirpath_chisel'      : dirpath_chisel,
                          'dirpath_sims'        : dirpath_sims,
                          'dirpath_foregrounds' : dirpath_foregrounds,
                          'dirpath_mask'        : dirpath_mask,
                          'dirpath_savedata'    : dirpath_savedata,                         
                          'filename_foregrounds': filename_foregrounds,                
                          'filename_mask'       : filename_mask,
                          'filepath_field1'     : filepath_field1,
                          'filepath_field2'     : filepath_field2,
                          'filepath_cross'      : filepath_cross,
                          'pathout'             : pathout
                       })


#params_sims['lmin'] = 0
params_sims['l'] = np.arange(params_sims['lmin'], params_sims['lmax']+1, 1)	

try:
    if params_tf['auto_hi_bins']<0:
        params_tf['auto_hi_bins'] = np.arange(params_sims['nbands'], dtype=np.int32)
except:
    pass
### rechosen to have 1<alpha<3
#bins_f1f2 = {
            #0: np.array([2, 3]),    
            #1: np.array([4, 5, 6, 7, 8, 9]),
            #2: np.array([10,11,12,13,14]),
            #1: np.array([7, 8, 9]),
            #2: np.array([11,12,13,14]),
            #3: np.array([15,16,17,18]),
            #4: np.array([19,20,21]),
            #5: np.array([22,23,24]),
            #6: np.array([25,26]),#,27]),
            #7: np.array([28,29])    
#            }    
try:    
    if params_tf['cross_hi_bins']<0:
    #    params_inout['output']['cross_hi_bins'] = np.arange(params_fields['hi']['nbands'], dtype=np.int32)
        if params_tf['auto_g_bin']==0:
            params_tf['cross_hi_bins'] = np.array([0,1,2,3])
        elif params_tf['auto_g_bin']==1:
            params_tf['cross_hi_bins'] = np.array([4,5,6,7,8,9])
        elif params_tf['auto_g_bin']==2:
            params_tf['cross_hi_bins'] = np.array([10,11,12,13,14])
        elif params_tf['auto_g_bin']==3:
            params_tf['cross_hi_bins'] = np.array([15,16,17,18])
        elif params_tf['auto_g_bin']==4:
            params_tf['cross_hi_bins'] = np.array([18,19,20,21])
        elif params_tf['auto_g_bin']==5:
            params_tf['cross_hi_bins'] = np.array([22,23,24])
        elif params_tf['cross_g_bin']==6:
            params_tf['cross_hi_bins'] = np.array([25,26,27])        
        elif params_tf['auto_g_bin']==7:
            params_tf['cross_hi_bins'] = np.array([28,29])       
        else:
            params_tf['cross_hi_bins'] = np.array([28,29])             
except:
    pass
params_tf['cross_g_bin'] = dcopy( np.asarray( params_tf['cross_g_bin'] ) )
params_tf['auto_g_bin' ] = dcopy( np.asarray( params_tf['auto_g_bin' ] ) )  
########################################################################################################################################
del verbose, project, threshold#, load_simulations
del nrealizations, nside, lmin, lmax, apply_mask,  include_wn, include_beam, wn_sigma_pix, beam_fwhm
del numin, numax, nbands, unit, stokes #,apply_namaster, namaster_delell
del tf_hixhi, tf_gxhi, cross_g_bin, auto_hi_bins, cross_hi_bins
#del save_hi_plus_fg, save_hi_plus_fg_smoothed, save_hi_plus_fg_plus_wn, save_hi_plus_fg_plus_wn_smoothed, 
del save_hdf5, save_txt, save_simulated, save_theoretical
del method, wtransform, ns
del dirpath_chisel, dirpath_sims, dirpath_foregrounds, dirpath_mask, filename_foregrounds, filename_mask, pathout, 
del filepath_field1, filepath_field2, filepath_cross, dirpath_savedata
#####################################################################################################################################################################
#####################################################################################################################################################################
try:
    params_general['project' ] = params_general['project' ].format(b=params_tf['auto_g_bin'])    
except:
    pass
try:
    params_inout['dirpath_savedata' ] = params_inout['dirpath_savedata' ].format(b=params_tf['auto_g_bin'])    
except:
    pass         
try:
    params_inout['dirpath_sims' ] = params_inout['dirpath_sims' ].format(b=params_tf['auto_g_bin'])    
except:
    pass     
try:    
    params_inout['filepath_field1'] = params_inout['filepath_field1'].format(b=params_tf['auto_g_bin'])    
except:
    pass
try:
    params_inout['filepath_field2'] = params_inout['filepath_field2'].format(b=params_tf['auto_g_bin'])    
except:
    pass   
try:    
    params_inout['filepath_cross' ] = params_inout['filepath_cross' ].format(b=params_tf['auto_g_bin'])    
except:
    pass       
#####################################################################################################################################################################    
#####################################################################################################################################################################
if 1:
    print('\n===================================================================')
    print('===================================================================')
    print('Starting the code...')
    print('===================================================================')    
    print('===================================================================\n')
#####################################################################################################
if params_general['verbose']:
    print('\n--Variables to be used:--')
    print(params_general)    
    print(params_sims   )     
    print(params_tf     )
    print(params_chisel )
    print(params_inout  )
    #print(params_fits)
    print()

sys.path.insert(1, params_inout['dirpath_chisel'])    
import Extension4BINGO   as cs
import load_standard_params as loadparams
#import pyMRS as pymrs
#import statcosmo as statc
params_maps    = loadparams.load_params_maps()
params_CS      = loadparams.load_params_CS()
params_WT      = loadparams.load_params_WT()
#params_path_cs = loadparams.load_params_path('amarins')
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
    #print(params_path_cs)
#####################################################################################################################################################################
#####################################################################################################################################################################
#loading Foregrodun (FG) maps and the mask
Mfg = hdata.getmap(dirpath_=params_inout.dirpath_foregrounds, filename_=params_inout.filename_foregrounds, healpix_readingformat=False, hdu=1)
if (params_sims.apply_mask) and (params_inout.dirpath_mask!="") and (params_inout.filename_mask!=""):
    MASK = hdata.getmap(dirpath_=params_inout.dirpath_mask, filename_=params_inout.filename_mask, healpix_readingformat=False, hdu=1)
    params_chisel['coverage_flag'] = 'masked'
    if params_inout['pathout'] =="":
        params_inout['pathout'] = 'masked'    
else:
    MASK = np.ones_like(Mfg[0,:])
    params_chisel['coverage_flag'] = 'fullsky'
    if params_inout['pathout'] =="":
        params_inout['pathout'] = 'fullsky'

#hdata.file_verification( params_inout['dirpath_savedata'], '', '')#, params_inout['pathout'] )
#params_inout['dirpath_savedata'] = os.path.join(params_inout['dirpath_savedata'], params_inout['pathout'] )
hdata.file_verification(params_inout['dirpath_savedata'],'','ns{}'.format(params_CS['ns']))
params_inout['dirpath_savedata_ns'] = os.path.join(params_inout['dirpath_savedata'],'ns{}'.format(params_CS['ns']))   
#print(params_general['project' ])
#print(params_inout['dirpath_sims' ])
#print(params_inout['filepath_field1' ])
#print(params_inout['filepath_field2' ])
#print(params_inout['filepath_cross' ])
#print(params_inout['dirpath_savedata' ])
#print(params_inout['dirpath_savedata_ns' ])
#####################################################################################################################################################################
#####################################################################################################################################################################
output_info = {'field':'field', 'nside':params_sims.nside, 
               'frequency'   :{'min':params_sims.numin, 'max':params_sims.numax,'nbands':params_sims.nbands,  'unit':params_sims.unit}, 
               'stokes'      :params_sims['stokes'],
               'coverage'    :None,#params_chisel['coverage_flag'], 
               'pathout'     :params_inout['dirpath_savedata_ns'],
               'apply_beam'  :False, 
               'namePL'      :False, 
               'nameFG'      :params_inout.filename_foregrounds, 
               'apply_mask'  :params_sims.apply_mask, 
               'nameM'       :params_inout.filename_mask,               
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
if params_general['verbose']:
    print('\n------------------------------------------------------------------')
    print('Loading the APS theoretically simulated')
    print('------------------------------------------------------------------\n')

clf1 = np.loadtxt(params_inout['filepath_field1']).T
clf2 = np.loadtxt(params_inout['filepath_field2']).T
clcx = np.loadtxt(params_inout['filepath_cross'] ).T 

load_sims=1
if load_sims:
    itime         = time.time() #['alm_hi_sim', 'cl_hi_sim', 'cl_kappa_sim', 'cl_cross_sim', 'alm_kappa_sim']
    varnames      = ['alm_hi_sim','alm_kappa_sim' , 'cl_hi_sim', 'cl_kappa_sim', 'cl_cross_sim']#,'cross_correlation_coef']
    dict_all      = cxfs.get_file_simulated_data(path_=params_inout['dirpath_sims'], varnames=varnames, nsims_to_use=params_sims.nreals)
    dict_all_mean = cxfs.dict_averages_from_loaded_data(dict_all_=dict_all, varnames_=varnames)
    nsims         = int(dict_all_mean[varnames[0]]['bin1']['matrix'].shape[0]    )
    if params_general['verbose']: print('Loading time for {0:d} sims: {1:.2f} seg'.format(nsims, time.time()-itime)) 
    del nsims
    
params_sims['realizations'] = np.asarray(list(dict_all.keys()))
if params_general['verbose']: print("Using sims: ",params_sims['realizations'])
#####################################################################################################################################################################
#####################################################################################################################################################################
N=params_sims['nbands']
extend_names = [ 'smoothed', 'smoothless',
                 'noisy'   , 'noiseless' ,
               ]
filename = ""
if params_general['verbose']: 
    print('\n------------------------------------------------------------------')
    print('Starting the Foreground Removal Process')
    print('And the Transfer Function estimation   ')
    print('------------------------------------------------------------------\n')
for j, jsim in enumerate(params_sims['realizations'][:]):
    if params_general['verbose']: print('Job {} -- Sim: {}'.format(j+1, jsim))
    alm_hi_jsim = dict_all[jsim]['alm_hi_sim'][2:]
    for num in range(alm_hi_jsim.shape[0]):
        alm_hi_num = np.ascontiguousarray(alm_hi_jsim[num,:])
        if not num:  Mhi = hp.alm2map( alm_hi_num,                  params_sims.nside, pol=False )    
        else:        Mhi = np.vstack(( Mhi, hp.alm2map( alm_hi_num, params_sims.nside, pol=False ) ))       
    ###
    params_maps["getdata"] = "observed"
    if params_sims['include_beam']:
        if not j: filename="{}{}".format(extend_names[0],filename)
        Xsky = np.vstack([ hp.smoothing(iMfg+iMhi, fwhm=params_sims['beam_fwhm'], pol=0) for (iMfg,iMhi) in zip(Mfg,Mhi) ])
        Xhi  = np.vstack([ hp.smoothing(m,         fwhm=params_sims['beam_fwhm'], pol=0) for m in Mhi ])
        Xfg  = np.vstack([ hp.smoothing(m,         fwhm=params_sims['beam_fwhm'], pol=0) for m in Mfg ])
    else:        
        if not j: filename="{}{}".format(extend_names[1],filename)
        Xsky = Mfg+Mhi
        Xhi  = Mhi
        Xfg  = Mfg
    if params_sims['include_wn']:
        if not j: filename="{}_{}".format(extend_names[2],filename)
        WN    = np.vstack([ np.random.normal(loc=0,scale=params_sims['wn_sigma_pix'], size=12*params_sims['nside']**2) for i in range(N) ])        
    else: 
        if not j: filename="{}_{}".format(extend_names[3],filename)
        WN = np.zeros_like(Mhi)        
    params_cs, params_wt = cs.load(params_CS, params_WT)
    ###
    X  = dcopy( Xsky + WN )*MASK
    X  = cs.adaptation_maps(X, params_maps, params_inout)
    Xr = cs.maps2CSmaps(    X, params_wt  , params_cs )# , False)
    ###
    Ae  = Xr['mixmatrix']
    Wfg = np.dot( Ae, np.dot( np.linalg.inv(np.dot(Ae.T,Ae)) ,Ae.T))
    Xfg = np.dot(Wfg,Xfg)
    ##
    alm_k_sim  = dict_all[jsim]['alm_kappa_sim'][2:,:]
    alm_wn_sim = 0  if not params_sims['include_wn'] else np.vstack([ hp.map2alm(maps=m, pol=False) for m in WN*MASK ])
    alm_hi_sim = np.vstack([ hp.map2alm(maps=m, pol=False) for m in Xhi ])#convolved
    alm_hi_est = np.vstack([ hp.map2alm(maps=m, pol=False) for m in Xr['reconstruction']['21cm'] ])
    #alm_fg_est = np.vstack([ hp.map2alm(maps=m, pol=False) for m in Xfg ])#convolved
    ###
    #Alm_k_sim  = alm_k_sim            if not j else np.vstack(( Alm_k_sim , alm_k_sim            ))
    #Alm_hi_sim = alm_hi_sim.flatten() if not j else np.vstack(( Alm_hi_sim, alm_hi_sim.flatten() ))
    #Alm_hi_est = alm_hi_est.flatten() if not j else np.vstack(( Alm_hi_est, alm_hi_est.flatten() ))
    #Alm_fg_est = alm_fg_est.flatten() if not j else np.vstack(( Alm_fg_est, alm_fg_est.flatten() ))
    ###
    ch_sim = np.vstack([ hp.alm2cl(alms1=i_alm) for i_alm in alm_hi_sim+alm_wn_sim ])
    ch_est = np.vstack([ hp.alm2cl(alms1=i_alm) for i_alm in alm_hi_est ])
    cx_sim = np.vstack([ hp.alm2cl(alms1=alm_k_sim, alms2=i_alm) for i_alm in alm_hi_sim+alm_wn_sim ])
    cx_est = np.vstack([ hp.alm2cl(alms1=alm_k_sim, alms2=i_alm) for i_alm in alm_hi_est ])

    if params_general['threshold']>0:
        idx = np.where((ch_sim<params_general['threshold'])*(ch_sim>-params_general['threshold']))
        ch_sim[idx]=0
        idx = np.where((cx_sim<params_general['threshold'])*(cx_sim>-params_general['threshold']))
        cx_sim[idx]=0
        #idx = np.where((ch_est<params_general['threshold'])*(ch_est>-params_general['threshold']))
        #ch_est[idx]=0        
        #idx = np.where((cx_est<params_general['threshold'])*(cx_est>-params_general['threshold']))
        #cx_est[idx]=0    
    #cx_gfg_est = np.vstack([ hp.alm2cl(alms1=alm_k_sim, alms2=i_alm) for i_alm in alm_fg_est ])
    ###
    if params_tf.tf_hixhi: alpha_ch      = ch_sim/ch_est
    if params_tf.tf_gxhi:  alpha_cx      = cx_sim/cx_est
    if params_tf.tf_hixhi: alpha_ch_flat = alpha_ch.flatten() if not j else np.vstack([ alpha_ch_flat, alpha_ch.flatten() ])
    if params_tf.tf_gxhi:  alpha_cx_flat = alpha_cx.flatten() if not j else np.vstack([ alpha_cx_flat, alpha_cx.flatten() ])
    #alpha_gfg_cx = cx_sim/cx_gfg_est
    #alpha_gfg_cx_flat = alpha_gfg_cx.flatten() if not j else np.vstack([ alpha_gfg_cx_flat, alpha_gfg_cx.flatten() ])
del X, Xr, Xfg, Xsky, Xhi, alm_hi_num
del alpha_ch, alpha_cx, ch_sim, ch_est, cx_sim, cx_est, alm_hi_est, alm_hi_sim, alm_wn_sim, alm_k_sim
#
if params_general['verbose']:  print( 'Parameters loaded.' )        

#####################################################################################################################################################################                  
if params_general['verbose']: 
    print('\n------------------------------------------------------------------')
    print('Saving')
    print('------------------------------------------------------------------\n')

#####################
### AUTO
if params_tf.tf_hixhi and params_tf.save_theoretical:     
    prefix = 'alpha_ch_thr'
    fsky  = MASK[MASK>0].size/MASK.size
    Opix  = hp.nside2pixarea(params_sims['nside'])
    Nl_th = fsky*(params_sims['wn_sigma_pix']**2)*Opix
    Nl_th = Nl_th*np.ones(3*params_sims['nside'])
    blf1  = 1 if not params_sims['include_beam'] else fnoise.bl_function(params_sims['beam_fwhm'], input_unit="radians", from_real_space=False)[:-1]
    Cl_th = clf1[1:,:]
    Cl_th = Cl_th*(blf1**2)
    Cl_th = Cl_th if not params_sims['include_wn'] else Cl_th+Nl_th
    for i,nuj in enumerate(params_tf['auto_hi_bins']):
        alpha_thr = 1/( (1-2*Wfg[nuj,nuj]) + np.dot(Wfg**2,Cl_th)/Cl_th) 
        alpha_thr = alpha_thr[nuj,:]
        alpha_ch_thr = alpha_thr if not i else np.vstack([ alpha_ch_thr, alpha_thr ])
    del alpha_thr
    filename_cl = os.path.join( params_inout['dirpath_savedata_ns'], 
                                      '{}_{}'.format(prefix, filename) )
    if params_tf.save_hdf5:
        import h5py
        filename_cl_h5 = "{}.h5".format(filename_cl)
        if params_general['verbose']: print('saving at... {}'.format(filename_cl_h5))
        with h5py.File(filename_cl_h5, "w") as f:
            f.create_dataset("alpha"    , data=alpha_ch_thr             )
            f.create_dataset("mixmatrix", data=Ae                      )
            f.create_dataset("Wfg"      , data=Wfg                      )
            f.create_dataset("hi_bins"  , data=params_tf['auto_hi_bins'])
        del filename_cl_h5 
    if params_tf.save_txt:
        header = "{}{}".format("hi_bins: ", params_tf['auto_hi_bins'])
        filename_cl_txt = "{}.txt".format(filename_cl)
        np.savetxt(filename_cl_txt, alpha_ch_thr.T, fmt=["%e"]*alpha_ch_thr.shape[0], delimiter=" ", header=header)
        print('saving at... {}'.format(filename_cl_txt))                  
        
if params_tf.tf_hixhi and params_tf.save_simulated: 
    prefix = 'alpha_ch_sim'
    alpha_ch_sim = np.average(alpha_ch_flat,axis=0).reshape(N,-1)
    alpha_ch_sim = alpha_ch_sim[params_tf['auto_hi_bins'],:]
    filename_cl  = os.path.join( params_inout['dirpath_savedata_ns'], 
                                '{}_{}'.format(prefix, filename) )
    if params_tf.save_hdf5:
        import h5py
        filename_cl_h5 = "{}.h5".format(filename_cl)
        if params_general['verbose']: print('saving at... {}'.format(filename_cl_h5))
        with h5py.File(filename_cl_h5, "w") as f:
            f.create_dataset("alpha"    , data=alpha_ch_sim             )
            f.create_dataset("mixmatrix", data=Ae                      )
            f.create_dataset("Wfg"      , data=Wfg                      )
            f.create_dataset("hi_bins"  , data=params_tf['auto_hi_bins'])
        del filename_cl_h5    
    if params_tf.save_txt:
        header = "{}{}".format("hi_bins: ", params_tf['auto_hi_bins'])
        filename_cl_txt = "{}.txt".format(filename_cl)
        np.savetxt(filename_cl_txt, alpha_ch_sim.T, fmt=["%e"]*alpha_ch_sim.shape[0], delimiter=" ", header=header)
        print('saving at... {}'.format(filename_cl_txt))          
        
#####################
### CROSS
if params_tf.tf_gxhi and params_tf.save_theoretical:     
    prefix = 'alpha_cx_thr'
    Cl_th  = clcx[1:,:] 
    blf1   = 1 if not params_sims['include_beam'] else fnoise.bl_function(params_sims['beam_fwhm'], input_unit="radians", from_real_space=False)[:-1]
    Cl_th  = Cl_th*blf1
    alpha_cx_thr = 1/(1 - np.dot(Wfg,Cl_th)/(Cl_th))
    alpha_cx_thr = alpha_cx_thr[params_tf['cross_hi_bins'],:]
    del blf1
    filename_cl = os.path.join( params_inout['dirpath_savedata_ns'], 
                                      '{}_{}'.format(prefix, filename) )    
    if params_tf.save_hdf5:
        import h5py
        filename_cl_h5 = "{}.h5".format(filename_cl)
        if params_general['verbose']: print('saving at... {}'.format(filename_cl_h5))
        with h5py.File(filename_cl_h5, "w") as f: 
            f.create_dataset("alpha"    , data=alpha_cx_thr              )
            f.create_dataset("mixmatrix", data=Ae                      )
            f.create_dataset("Wfg"      , data=Wfg                       )
            f.create_dataset("g_bin"    , data=params_tf['cross_g_bin']  )
            f.create_dataset("hi_bins"  , data=params_tf['cross_hi_bins'])
        del filename_cl_h5      
    if params_tf.save_txt:
        header = "{}{}{}{}{}".format("g_bin: ",params_tf['cross_g_bin'], " | ", "hi_bins: ", params_tf['cross_hi_bins'])
        filename_cl_txt = "{}.txt".format(filename_cl)
        np.savetxt(filename_cl_txt, alpha_cx_thr.T, fmt=["%e"]*alpha_cx_thr.shape[0], delimiter=" ", header=header)
        print('saving at... {}'.format(filename_cl_txt))         

if params_tf.tf_gxhi and params_tf.save_simulated:     
    prefix = 'alpha_cx_sim'
    alpha_cx_sim = np.average(alpha_cx_flat,axis=0).reshape(N,-1)
    alpha_cx_sim = alpha_cx_sim[params_tf['cross_hi_bins'],:]
    filename_cl = os.path.join( params_inout['dirpath_savedata_ns'], 
                                      '{}_{}'.format(prefix, filename) )
    if params_tf.save_hdf5:
        import h5py
        filename_cl_h5 = "{}.h5".format(filename_cl)
        if params_general['verbose']: print('saving at... {}'.format(filename_cl_h5))
        with h5py.File(filename_cl_h5, "w") as f:
            f.create_dataset("alpha"    , data=alpha_cx_sim              )
            f.create_dataset("mixmatrix", data=Ae                      )
            f.create_dataset("Wfg"      , data=Wfg                       )
            f.create_dataset("g_bin"    , data=params_tf['cross_g_bin']  )
            f.create_dataset("hi_bins"  , data=params_tf['cross_hi_bins'])
        del filename_cl_h5

    if params_tf.save_txt:
        header = "{}{}{}{}{}".format("g_bin: ",params_tf['cross_g_bin'], " | ", "hi_bins: ", params_tf['cross_hi_bins'])
        filename_cl_txt = "{}.txt".format(filename_cl)
        np.savetxt(filename_cl_txt, alpha_cx_sim.T, fmt=["%e"]*alpha_cx_sim.shape[0], delimiter=" ", header=header)
        print('saving at... {}'.format(filename_cl_txt))                 

#####################################################################################################################################################################                  
if 1:
    print('\n------------------------------------------')
    print('END')
    print('------------------------------------------\n')
 