#sys.path.insert(1, os.getcwd())
import os, sys, time
from copy import deepcopy as dcopy
import cross_functions_theory as cxft
import cross_functions_simulations as cxfs
import cross_functions_theory_4forecast as c4ft
import handling_data  as hdata
import numpy as np
import healpy as hp
import pymaster as nmt
import pandas as pd
import camb
from   camb import model, initialpower
import json
import argparse
import warnings
warnings.filterwarnings("ignore")	
sys.path.insert(1, '/data/AMARINS/chisel/scripts')    
import Extension4BINGO   as cs
import load_standard_params as loadparams
#import pyMRS as pymrs
#import statcosmo as statc
#####################################################################################################
#####################################################################################################
# CONSTANTS
C_LIGHT   = 299792458/1e3 #km/s
nu_HI     = 1420.405751768 #MHz
#THRESHOLD = 1e-20
STERADIAN = np.radians(1/60)**2 #arcmin to steradian
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
PATH        = "/data/AMARINS/LSSxHI-CODES/scripts"
INI         = "generating_data_and_noise_4forecast.ini"
name_params = os.path.join(PATH,INI)
config.read(name_params)
############
#General
verbose   = config.getboolean("General","verbose"  )
threshold = config.getfloat(  "General","threshold")
############
#Simulations
lmin            = config.getint(     "Simulations","lmin"            )
lmax            = config.getint(     "Simulations","lmax"            )
apply_mask      = config.getboolean( "Simulations","apply_mask"      )
apply_tf        = config.getboolean( "Simulations","apply_tf"        )
#use_namaster    = config.getboolean( "Simulations","use_namaster"    )
namaster_delell = config.getint(     "Simulations","namaster_delell" )
ell_nmt_hi      = config.getboolean( "Simulations","ell_nmt_hi"      )
ell_nmt_g       = config.getboolean( "Simulations","ell_nmt_g"       )
ell_nmt_gxhi    = config.getboolean( "Simulations","ell_nmt_gxhi"    )

nside           = config.getint(     "Simulations","nside"           )
seed_hi         = config.getint(     "Simulations","seed_hi"         )
seed_k          = config.getint(     "Simulations","seed_k"          )
seed_tax        = config.getint(     "Simulations","seed_tax"        )
############
#HI
numin       = config.getint(   "HI" , "numin"       )
numax       = config.getint(   "HI" , "numax"       )
nbands      = config.getint(   "HI" , "nbands"      )
fwhm_arcmin = config.getint(   "HI" , "fwhm_arcmin" )
sigma_pix   = config.getfloat( "HI" , "sigma_pix"   )
############
#Galaxy
survey_name   = config.get( "Galaxy" , "survey_name"   )
data_release  = config.get( "Galaxy" , "data_release"  )
type_of_field = config.get( "Galaxy" , "type_of_field" )
path_to_file  = config.get( "Galaxy" , "path_to_file"  )
ntomos        = config.get( "Galaxy" , "ntomos"        )
sigmaz        = config.get( "Galaxy" , "sigmaz"        )
############
#Algorithm
method     = config.get( "Algorithm" , "method"     )
wtransform = config.get( "Algorithm" , "wtransform" )
ns         = config.get( "Algorithm" , "ns"         )
############
#outputs
auto_hi_bins   = config.get( "Outputs" , "auto_hi_bins" ) 
auto_hi_bins   = np.array([float(item.strip()) for item in auto_hi_bins.split(',')])
cross_hi_bins  = config.get( "Outputs" , "cross_hi_bins" ) 
cross_hi_bins  = np.array([float(item.strip()) for item in cross_hi_bins.split(',')])
auto_g_bin     = config.getint("Outputs" , "auto_g_bin"  )
cross_g_bin    = config.getint("Outputs" , "cross_g_bin" )
pathout        = config.get( "Outputs" , "pathout"       )
filename_hi    = config.get( "Outputs" , "filename_hi"   )
filename_g     = config.get( "Outputs" , "filename_g"    )
filename_gxhi  = config.get( "Outputs" , "filename_gxhi" )
#
dirpath_tf        = config.get( "Outputs" , "dirpath_tf"   )
filename_tf_hi    = config.get( "Outputs" , "filename_tf_hi"   )
filename_tf_gxhi  = config.get( "Outputs" , "filename_tf_gxhi" )
#
save_hi   = config.getboolean( "Outputs" , "save_hi"   )
save_g    = config.getboolean( "Outputs" , "save_g"    )
save_gxhi = config.getboolean( "Outputs" , "save_gxhi" )

save_theory    = config.getboolean( "Outputs" , "save_theory"    )
save_simulated = config.getboolean( "Outputs" , "save_simulated" )
save_estimated = config.getboolean( "Outputs" , "save_estimated" )

save_txt       = config.getboolean( "Outputs" , "save_txt" )
save_npz       = config.getboolean( "Outputs" , "save_npz" )
save_hdf5      = config.getboolean( "Outputs" , "save_hdf5")
############
#PATH
dirpath_chisel       = config.get( "PATH" , "dirpath_chisel"       )
dirpath_foregrounds  = config.get( "PATH" , "dirpath_foregrounds"  )
dirpath_mask         = config.get( "PATH" , "dirpath_mask"         )
filename_foregrounds = config.get( "PATH" , "filename_foregrounds" )
filename_mask        = config.get( "PATH" , "filename_mask"        )

filepath_hi   = config.get( "PATH" , "filepath_hi"   )
filepath_g    = config.get( "PATH" , "filepath_g"    )
filepath_gxhi = config.get( "PATH" , "filepath_gxhi" )
###############################################################################
# You can modify any options in the parameters.ini file by the command terminal
###############################################################################
parser = argparse.ArgumentParser(description='Modify by the command terminal parameters in {} file'.format(INI))

#General
parser.add_argument('--verbose'   , action = 'store', dest = 'verbose'  , default = verbose  , help = '')
parser.add_argument('--threshold' , action = 'store', dest = 'threshold', default = threshold, help = '')

#Simulations
parser.add_argument('--lmin'            , action = 'store', dest = 'lmin'           , default = lmin           , help = '')
parser.add_argument('--lmax'            , action = 'store', dest = 'lmax'           , default = lmax           , help = '')
parser.add_argument('--apply_mask'      , action = 'store', dest = 'apply_mask'     , default = apply_mask     , help = '')
parser.add_argument('--apply_tf'        , action = 'store', dest = 'apply_tf'       , default = apply_tf       , help = '')
#parser.add_argument('--use_namaster'    , action = 'store', dest = 'use_namaster'   , default = use_namaster   , help = '')
parser.add_argument('--namaster_delell' , action = 'store', dest = 'namaster_delell', default = namaster_delell, help = '')
parser.add_argument('--ell_nmt_hi'      , action = 'store', dest = 'ell_nmt_hi'     , default = ell_nmt_hi     , help = '')
parser.add_argument('--ell_nmt_g'       , action = 'store', dest = 'ell_nmt_g'      , default = ell_nmt_g      , help = '')
parser.add_argument('--ell_nmt_gxhi'    , action = 'store', dest = 'ell_nmt_gxhi'   , default = ell_nmt_gxhi   , help = '')


parser.add_argument('--nside'           , action = 'store', dest = 'nside'          , default = nside          , help = '')
parser.add_argument('--seed_hi'         , action = 'store', dest = 'seed_hi'        , default = seed_hi        , help = '')
parser.add_argument('--seed_k'          , action = 'store', dest = 'seed_k'         , default = seed_k         , help = '')
parser.add_argument('--seed_tax'        , action = 'store', dest = 'seed_tax'       , default = seed_tax       , help = '')

#HI
parser.add_argument('--numin'      , action = 'store', dest = 'numin'       , default = numin      , help = '')
parser.add_argument('--numax'      , action = 'store', dest = 'numax'       , default = numax      , help = '')
parser.add_argument('--nbands'     , action = 'store', dest = 'nbands'      , default = nbands     , help = '')
parser.add_argument('--fwhm_arcmin', action = 'store', dest = 'fwhm_arcmin' , default = fwhm_arcmin, help = '')
parser.add_argument('--sigma_pix'  , action = 'store', dest = 'sigma_pix'   , default = sigma_pix  , help = '')

#Galaxy
parser.add_argument('--survey_name'  , action = 'store', dest = 'survey_name'  , default = survey_name  , help = '')
parser.add_argument('--data_release' , action = 'store', dest = 'data_release' , default = data_release , help = '')
parser.add_argument('--type_of_field', action = 'store', dest = 'type_of_field', default = type_of_field, help = '')
parser.add_argument('--path_to_file' , action = 'store', dest = 'path_to_file' , default = path_to_file , help = '')
parser.add_argument('--ntomos'       , action = 'store', dest = 'ntomos'       , default = ntomos       , help = '')
parser.add_argument('--sigmaz'       , action = 'store', dest = 'sigmaz'       , default = sigmaz       , help = '')

#Algorithm
parser.add_argument('--method'    , action = 'store', dest = 'method'    , default = method    , help = '')
parser.add_argument('--wtransform', action = 'store', dest = 'wtransform', default = wtransform, help = '')
parser.add_argument('--ns'        , action = 'store', dest = 'ns'        , default = ns        , help = '')

#Outputs
parser.add_argument('--auto_hi_bins' , action = 'store', dest = 'auto_hi_bins' , default = auto_hi_bins , nargs='+', type=float)
parser.add_argument('--cross_hi_bins', action = 'store', dest = 'cross_hi_bins', default = cross_hi_bins, nargs='+', type=float)
parser.add_argument('--auto_g_bin'   , action = 'store', dest = 'auto_g_bin'   , default = auto_g_bin   , help = '')
parser.add_argument('--cross_g_bin'  , action = 'store', dest = 'cross_g_bin'  , default = cross_g_bin  , help = '')
parser.add_argument('--pathout'      , action = 'store', dest = 'pathout'      , default = pathout      , help = '')
parser.add_argument('--filename_hi'  , action = 'store', dest = 'filename_hi'  , default = filename_hi  , help = '')
parser.add_argument('--filename_g'   , action = 'store', dest = 'filename_g'   , default = filename_g   , help = '')
parser.add_argument('--filename_gxhi', action = 'store', dest = 'filename_gxhi', default = filename_gxhi, help = '')
#
parser.add_argument('--dirpath_tf'      , action = 'store', dest = 'dirpath_tf'      , default = dirpath_tf      , help = '')
parser.add_argument('--filename_tf_hi'  , action = 'store', dest = 'filename_tf_hi'  , default = filename_tf_hi  , help = '')
parser.add_argument('--filename_tf_gxhi', action = 'store', dest = 'filename_tf_gxhi', default = filename_tf_gxhi, help = '')
#
parser.add_argument('--save_hi'   , action = 'store', dest = 'save_hi'  , default = save_hi  , help = '')
parser.add_argument('--save_g'    , action = 'store', dest = 'save_g'   , default = save_g   , help = '')
parser.add_argument('--save_gxhi' , action = 'store', dest = 'save_gxhi', default = save_gxhi, help = '')

parser.add_argument('--save_theory'   , action = 'store', dest = 'save_theory'   , default = save_theory    , help = '')
parser.add_argument('--save_simulated', action = 'store', dest = 'save_simulated', default = save_simulated , help = '')
parser.add_argument('--save_estimated', action = 'store', dest = 'save_estimated', default = save_estimated , help = '')

parser.add_argument('--save_txt' , action = 'store', dest = 'save_txt' , default = save_txt , help = '')
parser.add_argument('--save_npz' , action = 'store', dest = 'save_npz' , default = save_npz , help = '')
parser.add_argument('--save_hdf5', action = 'store', dest = 'save_hdf5', default = save_hdf5, help = '')

#PATH
parser.add_argument('--dirpath_chisel'      , action = 'store', dest = 'dirpath_chisel'      , default = dirpath_chisel      , help = '')
parser.add_argument('--dirpath_foregrounds' , action = 'store', dest = 'dirpath_foregrounds' , default = dirpath_foregrounds , help = '')
parser.add_argument('--dirpath_mask'        , action = 'store', dest = 'dirpath_mask'        , default = dirpath_mask        , help = '')
parser.add_argument('--filename_foregrounds', action = 'store', dest = 'filename_foregrounds', default = filename_foregrounds, help = '')
parser.add_argument('--filename_mask'       , action = 'store', dest = 'filename_mask'       , default = filename_mask       , help = '')
parser.add_argument('--filepath_hi'         , action = 'store', dest = 'filepath_hi'         , default = filepath_hi         , help = '')
parser.add_argument('--filepath_g'          , action = 'store', dest = 'filepath_g'          , default = filepath_g          , help = '')
parser.add_argument('--filepath_gxhi'       , action = 'store', dest = 'filepath_gxhi'       , default = filepath_gxhi       , help = '')

arguments = parser.parse_args()
###############################################################################
# Variables
###############################################################################
#GENERAL
verbose   = bool( arguments.verbose)
threshold = float(arguments.threshold)

#Simulations
lmin            = int(arguments.lmin)
lmax            = int(arguments.lmax)
apply_mask      = bool(arguments.apply_mask)
apply_tf        = bool(arguments.apply_tf)
#use_namaster    = bool(arguments.use_namaster)
namaster_delell = int(arguments.namaster_delell)
ell_nmt_hi      = bool(arguments.ell_nmt_hi)
ell_nmt_g       = bool(arguments.ell_nmt_g)
ell_nmt_gxhi    = bool(arguments.ell_nmt_gxhi)

nside           = int(arguments.nside)
seed_hi         = int(arguments.seed_hi)
seed_k          = int(arguments.seed_k)
seed_tax        = int(arguments.seed_tax)

#HI
numin       = int(arguments.numin)
numax       = int(arguments.numax)
nbands      = int(arguments.nbands)
fwhm_arcmin = float(arguments.fwhm_arcmin)
sigma_pix   = float(arguments.sigma_pix)

#Galaxy
survey_name   = str(arguments.survey_name)
data_release  = str(arguments.data_release)
type_of_field = str(arguments.type_of_field)
path_to_file  = str(arguments.path_to_file)
ntomos        = int(arguments.ntomos)
sigmaz        = float(arguments.sigmaz)

#Algorithms
method     = str(arguments.method)
wtransform = str(arguments.wtransform)
ns         = int(arguments.ns)

#Outputs
auto_hi_bins  = np.asarray(arguments.auto_hi_bins, dtype=np.int32)
cross_hi_bins = np.asarray(arguments.cross_hi_bins, dtype=np.int32)
auto_g_bin    = int(arguments.auto_g_bin)#np.asarray(arguments.auto_g_bins, dtype=np.int32)
cross_g_bin   = int(arguments.cross_g_bin)#np.asarray(arguments.cross_g_bins, dtype=np.int32)
pathout       = str(arguments.pathout)
filename_hi   = str(arguments.filename_hi)
filename_g    = str(arguments.filename_g)
filename_gxhi = str(arguments.filename_gxhi)
#
dirpath_tf       = str(arguments.dirpath_tf)
filename_tf_hi   = str(arguments.filename_tf_hi)
filename_tf_gxhi = str(arguments.filename_tf_gxhi)
#
save_hi   = bool(arguments.save_hi  )
save_g    = bool(arguments.save_g   )
savegx_hi = bool(arguments.save_gxhi)

save_theory     = bool(arguments.save_theory)
save_simulated  = bool(arguments.save_simulated)
save_estimated  = bool(arguments.save_estimated)

save_txt  = bool(arguments.save_txt)
save_npz  = bool(arguments.save_npz)
save_hdf5 = bool(arguments.save_hdf5)

#PATH
dirpath_chisel       = str(arguments.dirpath_chisel)
dirpath_foregrounds  = str(arguments.dirpath_foregrounds)
dirpath_mask         = str(arguments.dirpath_mask)
filename_foregrounds = str(arguments.filename_foregrounds)
filename_mask        = str(arguments.filename_mask)
filepath_hi          = str(arguments.filepath_hi)
filepath_g           = str(arguments.filepath_g)
filepath_gxhi        = str(arguments.filepath_gxhi)
#####################################################################################################
#####################################################################################################
if filename_mask=="":
    filename_mask=None

lell  = np.arange(int(3*nside))
fell  = lell*(lell+1)/2/np.pi           
if ell_nmt_hi+ell_nmt_g+ell_nmt_gxhi:
    b_nmt = nmt.NmtBin.from_nside_linear(nside, nlb=namaster_delell)     
    leff  = b_nmt.get_effective_ells()
    feff  = leff*(leff+1)/2/np.pi    
else:
    b_nmt = None
    leff  = dcopy(lell)
    feff  = dcopy(fell)

mask_nmt = (leff>=lmin)*(leff<=lmax)    
mask_ell = (lell>=lmin)*(lell<=lmax)    
#####################################################################################################
params_general = {'verbose'  :verbose,
                  'threshold':threshold}

params_APS = {'lmin'         : lmin, 
              'lmax'         : lmax, 
              'apply_mask'   : apply_mask,  
              'apply_tf'     : apply_tf,               
              'nside'        : nside, 
              'npix'         : 12*nside*nside,
              'seed_hi'      : seed_hi, 
              'seed_k'       : seed_k, 
              'seed_tax'     : seed_tax,
              'leff'         : lell,
              'feff'         : fell,
              'mask_ell'     : mask_ell,
              'namaster':{
                          #'use_namaster' : use_namaster, 
                          'ell_hi'   : ell_nmt_hi,
                          'ell_g'    : ell_nmt_g,
                          'ell_gxhi' : ell_nmt_gxhi,
                          'del_l'    : namaster_delell,                  
                          'b'        : b_nmt,
                          'leff'     : leff,
                          'feff'     : feff,
                          'mask_ell' : mask_nmt
                          }
             }
params_fields = {'hi': {
                        'nside'      : nside,
                        'numin'      : numin,
                        'numax'      : numax,
                        'nbands'     : nbands,
                        'fwhm_arcmin': fwhm_arcmin,
                        'sigma_pix'  : sigma_pix,    
                        'algorithm':{ 
                                     'method'    : method,
                                     'wtransform': wtransform,
                                     'ns'        : ns
                                    }
                       },
                 'galaxy':{
                          'survey_name'  : survey_name,    
                          'data_release' : data_release,    
                          'type_of_field': type_of_field,    
                          'path_to_file' : path_to_file,    
                          'ntomos'       : ntomos,                         
                          'sigmaz'       : sigmaz,                        
                          }
                }

params_inout = {
                'output':{
                        'auto_hi_bins'  : auto_hi_bins,
                        'cross_hi_bins' : cross_hi_bins,
                        'auto_g_bin'    : auto_g_bin,
                        'cross_g_bin'   : cross_g_bin,                    
                        'pathout'       : pathout,
                        'filename_hi'   : filename_hi,
                        'filename_g'    : filename_g,
                        'filename_gxhi' : filename_gxhi,    
                        'dirpath_tf'      : dirpath_tf,
                        'filename_tf_hi'  : filename_tf_hi,
                        'filename_tf_gxhi': filename_tf_gxhi,                      
                        'save_hi'       : save_hi,
                        'save_g'        : save_g,
                        'save_gxhi'     : save_gxhi,                    
                        'save_theory'   : save_theory,
                        'save_simulated': save_simulated,
                        'save_estimated': save_estimated,                    
                        'save_txt'      : save_txt,
                        'save_npz'      : save_npz,                    
                        'save_hdf5'     : save_hdf5,
                        },
                'input':{
                        'dirpath_chisel'      : dirpath_chisel,
                        'dirpath_foregrounds' : dirpath_foregrounds,
                        'dirpath_mask'        : dirpath_mask,
                        'filename_foregrounds': filename_foregrounds,    
                        'filename_mask'       : filename_mask,    
                        'filepath_hi'         : filepath_hi,    
                        'filepath_g'          : filepath_g,    
                        'filepath_gxhi'       : filepath_gxhi,    
                        },    
               }
params_APS['l'] = np.arange(params_APS['lmin'], params_APS['lmax']+1, 1)	
try:
    if params_inout['output']['auto_hi_bins']<0:
        params_inout['output']['auto_hi_bins'] = np.arange(params_fields['hi']['nbands'], dtype=np.int32)
except:
    pass
try:    
    if params_inout['output']['cross_hi_bins']<0:
    #    params_inout['output']['cross_hi_bins'] = np.arange(params_fields['hi']['nbands'], dtype=np.int32)
        if   params_inout['output']['auto_g_bin']==0:
            params_inout['output']['cross_hi_bins'] = np.array([0,1,2,3])
        elif params_inout['output']['auto_g_bin']==1:
            params_inout['output']['cross_hi_bins'] = np.array([4,5,6,7,8,9])
        elif params_inout['output']['auto_g_bin']==2:
            params_inout['output']['cross_hi_bins'] = np.array([10,11,12,13,14])
        elif params_inout['output']['auto_g_bin']==3:
            params_inout['output']['cross_hi_bins'] = np.array([15,16,17,18])
        elif params_inout['output']['auto_g_bin']==4:
            params_inout['output']['cross_hi_bins'] = np.array([18,19,20,21])
        elif params_inout['output']['auto_g_bin']==5:
            params_inout['output']['cross_hi_bins'] = np.array([22,23,24])
        elif params_inout['output']['auto_g_bin']==6:
            params_inout['output']['cross_hi_bins'] = np.array([25,26,27])        
        elif params_inout['output']['auto_g_bin']==7:
            params_inout['output']['cross_hi_bins'] = np.array([28,29])       
        else:
            params_inout['output']['cross_hi_bins'] = np.array([28,29])             
except:
    pass
params_inout['output']['cross_g_bin'] = np.array( [params_inout['output']['cross_g_bin']] )
params_inout['output']['auto_g_bin' ] = np.array( [params_inout['output']['auto_g_bin']] )   
########################################################################################################################################
del verbose, threshold
del lmin   , lmax , apply_mask, ell_nmt_hi, ell_nmt_g, ell_nmt_gxhi, namaster_delell, nside, seed_hi, seed_k, seed_tax
del numin  , numax, nbands, fwhm_arcmin, sigma_pix, method, wtransform, ns, survey_name, data_release, type_of_field, path_to_file, ntomos
del sigmaz , cross_hi_bins, auto_hi_bins, cross_g_bin, auto_g_bin, pathout, filename_hi, filename_g, filename_gxhi, 
del dirpath_tf, filename_tf_gxhi, filename_tf_hi
del dirpath_chisel, dirpath_foregrounds, dirpath_mask, filename_foregrounds, filename_mask, filepath_hi, filepath_g, 
del filepath_gxhi , b_nmt, leff , feff, save_txt, save_npz, save_hdf5, save_theory, save_simulated, save_estimated
del save_hi, save_g, save_gxhi, mask_nmt, mask_ell
########################################################################################################################################
#### 3D MATTER POWER SPECTRUM from CAMB ################################################################################################
if params_general['verbose']:
    print('\n===========================================')
    print('Starting the code')
    print('===========================================\n')
    #
    print('Loading Angular Spectra...')    
#######
### TF
######
if params_APS['apply_tf']:
    import h5py
    pathd = os.path.join(params_inout['output']['dirpath_tf'],params_inout['output']['filename_tf_gxhi'])
    pathd = pathd.format(b=params_inout['output']['cross_g_bin'][0])  
    #pathd = '{p}/LSST_bin{b}_40_0001/ns4/alpha_ch_sim_noiseless_smoothless.h5'.format(p=path,b=g_bin)        
    ###
    # CX
    print("TF-cross: ", pathd )
    with h5py.File(pathd, "r") as f:
        hdf5_dict = {key: f[key][()] for key in f.keys()}        
    alpha_cx = dcopy( hdf5_dict) 
    del hdf5_dict  
    idx = np.array([ np.where(ibin==alpha_cx['hi_bins'])[0][0] for ibin in params_inout['output']['cross_hi_bins'] ])
    try:
        idx = np.array([ np.where(ibin==alpha_cx['hi_bins'])[0][0] for ibin in params_inout['output']['cross_hi_bins'] ])
    except:
        raise ValueError("Cross-correlation bins from Transfer Function and required for analysis do not match.")
    alphacx_dict = {'alpha'  : np.vstack([ alpha_cx['alpha'][i] for i in idx ]), 
                    'idx'    : idx,
                    'hi_bins': np.array([ alpha_cx['hi_bins'][i] for i in idx ])
                   }
    ###
    # CH
    pathd = os.path.join(params_inout['output']['dirpath_tf'],params_inout['output']['filename_tf_hi'])
    pathd = pathd.format(b=params_inout['output']['cross_g_bin'][0])  
    print("TF-auto: ", pathd )
    with h5py.File(pathd, "r") as f:
        hdf5_dict = {key: f[key][()] for key in f.keys()}        
    alpha_ch = dcopy( hdf5_dict) 
    del hdf5_dict  
    alphach_dict = {'alpha'  :dcopy(alpha_ch['alpha']),
                    'idx'    :dcopy(alpha_ch['hi_bins']),
                    'hi_bins':dcopy(alpha_ch['hi_bins'])}
    del alpha_cx, alpha_ch
#sys.exit(0)
#####
try:    
    params_inout['input']['filepath_hi'  ] = params_inout['input']['filepath_hi'  ].format(b=params_inout['output']['auto_g_bin'][0])    
except:
    pass
try:
    params_inout['input']['filepath_g'   ] = params_inout['input']['filepath_g'   ].format(b=params_inout['output']['auto_g_bin' ][0])    
except:
    pass    
try:
    params_inout['input']['filepath_gxhi'] = params_inout['input']['filepath_gxhi'].format(b=params_inout['output']['cross_g_bin'][0])    
except:
    pass    
print('Field 1: ',params_inout['input']['filepath_hi'  ] )    
clf1_th = np.loadtxt(params_inout['input']['filepath_hi'  ]).T[1:,:]
print('Field 2: ',params_inout['input']['filepath_g'  ] )    
clf2_th = np.loadtxt(params_inout['input']['filepath_g'   ]).T[1:,:]
print('Field 1x2: ',params_inout['input']['filepath_gxhi'  ] )    
clcx_th = np.loadtxt(params_inout['input']['filepath_gxhi']).T[1:,:]    
#################
pars_lsst = cxft.load_survey_parameters(survey_name = params_fields['galaxy']['survey_name'  ], 
                                        data_release= params_fields['galaxy']['data_release' ], 
                                        type_field  = params_fields['galaxy']['type_of_field'],
                                        path_file   = params_fields['galaxy']['path_to_file' ])
pars_lsst['nbins']   = params_fields['galaxy']['ntomos']
pars_lsst['sigma_z'] = params_fields['galaxy']['sigmaz']
params_fields['galaxy']['params_survey'] = dcopy(pars_lsst)
del pars_lsst
#################
### THEORY
Nl_hi = c4ft.noise_hydrogen(  sigma_pix  =params_fields['hi']['sigma_pix'],         ell_vec=None, nside=params_APS['nside'], f_sky=1, smoothed=False)
Nl_g  = c4ft.noise_clustering(pars_survey=params_fields['galaxy']['params_survey'], ell_vec=None, nside=params_APS['nside'])
###
if params_inout['output']['save_hi']:
    fsky=1
    bin_=None
    bnmt = params_APS['namaster' ]['b'    ] if params_APS['namaster']['ell_hi'] else None
    leff = params_APS['namaster' ]['leff' ] if params_APS['namaster']['ell_hi'] else params_APS['leff']
    dell = params_APS['namaster' ]['del_l'] if params_APS['namaster']['ell_hi'] else 1
###
if params_inout['output']['save_simulated'] or params_inout['output']['save_estimated']:
    if params_general['verbose']:
        print('Simulating the Angular Power Spectra...')   
    _dict_ = cxfs.dictionary_cross_simulations_quantities_from_matrix(clf1_  = clf1_th, 
                                                                      clf2_  = clf2_th, 
                                                                      clcx_  = dcopy(clcx_th),
                                                                      seed_hi= params_APS['seed_hi'], 
                                                                      seed_k = params_APS['seed_k'], 
                                                                      tax    = params_APS['seed_tax'],
                                                                      fact   = 1, 
                                                                      nsims  = 1,
                                                                      channel_min_corr=1,
                                                                      channel_max_corr=params_fields['hi']['nbands'], 
                                                                      beta=None, show_time=params_general['verbose'])  
        
        
    clf1_sim = dcopy( _dict_['sim0']['cl_hi_sim'   ] )
    clf2_sim = dcopy( np.expand_dims(_dict_['sim0']['cl_kappa_sim'], axis=0) )
    clcx_sim = dcopy( _dict_['sim0']['cl_cross_sim'] )    

if params_inout['output']['save_estimated']:    
    #################
    if params_general['verbose']:
        print('Loading maps...')   
    Mfg = hdata.getmap(dirpath_ = params_inout['input']['dirpath_foregrounds' ],
                       filename_= params_inout['input']['filename_foregrounds'], 
                       healpix_readingformat=False, hdu=1)   
    
    MASK = hdata.getmap(dirpath_ = params_inout['input']['dirpath_mask' ],
                        filename_= params_inout['input']['filename_mask'],
                        healpix_readingformat=False, hdu=1)       
    #################
    ### Estimated
    if params_general['verbose']:
        print('Calculating the Estimated HI maps')   
    params_maps    = loadparams.load_params_maps()
    params_CS      = loadparams.load_params_CS()
    params_WT      = loadparams.load_params_WT()
    params_path    = loadparams.load_params_path('amarins')
    ###
    params_CS['method']     = params_fields['hi']['algorithm']['method']
    params_CS['wtransform'] = params_fields['hi']['algorithm']['wtransform']
    params_CS['ns']         = params_fields['hi']['algorithm']['ns']
    ###            
    params_chisel = pd.Series(params_CS)            
    params_fits = pd.Series({'nside' :params_fields['hi']['nside'],
                             'numin' :params_fields['hi']['numin'],
                             'numax' :params_fields['hi']['numax'],
                             'nbands':params_fields['hi']['nbands']
                            })
    
    _var_='alm_hi_sim'
    mask = np.ones_like(MASK) if not params_APS['apply_mask'] else MASK
    fsky = mask[mask>0].size/mask.size
    for j, jsim in enumerate(_dict_.keys()):
        alm_hi_jsim = _dict_[jsim][_var_]#[1:]
        for num in range(alm_hi_jsim.shape[0]):
            alm_hi_num = np.ascontiguousarray(alm_hi_jsim[num,:])
            if not num:  Mhi = hp.alm2map( alm_hi_num,                  params_fits.nside, pol=False)    
            else:        Mhi = np.vstack(( Mhi, hp.alm2map( alm_hi_num, params_fits.nside, pol=False) ))       
        #
        params_cs, params_wt = cs.load(params_CS,params_WT)
        X    = dcopy(Mhi + Mfg)*mask
        X    = cs.adaptation_maps(X, params_maps, params_path)
        Xr   = cs.maps2CSmaps(    X, params_wt  , params_cs  )
        #
        XG   = hp.alm2map(np.nan_to_num(_dict_[  'sim0']['alm_kappa_sim']), pol=0, nside=params_fits.nside)
        XG  *= mask
        break
    
    for ibin in range(X.shape[0]):
        alm_hi   = hp.map2alm(Xr['reconstruction']['21cm'][ibin,:])
        alm_k    = hp.map2alm(XG)
        clf1_est = hp.alm2cl(alms1=alm_hi) if not ibin else np.vstack((clf1_est, hp.alm2cl(alms1=alm_hi) ))    
        clf2_est = hp.alm2cl(alms1=alm_k)  if not ibin else np.vstack((clf2_est, hp.alm2cl(alms1=alm_k) ))    
        clcx_est = hp.alm2cl(alms1=alm_k, alms2=alm_hi) if not ibin else np.vstack((clcx_est,  hp.alm2cl(alms1=alm_k, alms2=alm_hi ) ))    
    ###!!!!
    #Pay close attention here. The amplitude correction has only been applied to the specific bins, not all
    if params_APS['apply_tf']:
        for i,(idx,ibin) in enumerate(zip(alphacx_dict['idx'],alphacx_dict['hi_bins'])): 
            clcx_est[ibin]=clcx_est[ibin]*alphacx_dict['alpha'][idx]
        for i,(idx,ibin) in enumerate(zip(alphach_dict['idx'],alphach_dict['hi_bins'])): 
            clf1_est[ibin]=clf1_est[ibin]*alphach_dict['alpha'][idx]
        
    del params_CS, alm_hi_num, alm_hi_jsim, alm_hi, alm_k, X, Xr, XG, MASK, Mfg, Mhi
    del params_maps, params_WT, params_path, params_fits, params_chisel    
    ###!!!!
    
if params_inout['output']['save_hi']:
    if params_inout['output']['save_theory']:
        if params_general['verbose']:
            print('Calculating the HI cosmic noise (theoretical)')           
        CV_f1_th = np.array([ c4ft.cosmic_variance_ibin(clf1_     = clf1_th[i],
                                                        clf2_     = clf1_th[i],
                                                        clcx_     = clf1_th[i],
                                                        nlf1_     = Nl_hi     ,    
                                                        nlf2_     = Nl_hi     ,
                                                        nlcx_     = Nl_hi     ,
                                                        fwhmf1_   = params_fields['hi'    ]['fwhm_arcmin'  ],
                                                        fwhmf2_   = params_fields['hi'    ]['fwhm_arcmin'  ],
                                                        fwhmcx_   = params_fields['hi'    ]['fwhm_arcmin'  ],
                                                        l_eff     = leff, 
                                                        delta_ell = dell, 
                                                        b_nmt     = bnmt,
                                                        fsky      = fsky,
                                                        ) for i in np.arange(clf1_th.shape[0]) ])
        if params_general['verbose']:
            print('Theory: (num of C_ell, num of ell):',CV_f1_th.shape, '| delell: ',dell)
        
    if params_inout['output']['save_simulated']:
        if params_general['verbose']:
            print('Calculating the HI cosmic noise (simulated)')   
            
        CV_f1_sim = np.array([ c4ft.cosmic_variance_ibin(clf1_ = clf1_sim[i],
                                                         clf2_ = clf1_sim[i],
                                                         clcx_ = clf1_sim[i],
                                                         nlf1_ = Nl_hi ,
                                                         nlf2_ = Nl_hi ,
                                                         nlcx_ = Nl_hi,
                                                         fwhmf1_   = params_fields['hi'    ]['fwhm_arcmin'  ],
                                                         fwhmf2_   = params_fields['hi'    ]['fwhm_arcmin'  ],
                                                         fwhmcx_   = params_fields['hi'    ]['fwhm_arcmin'  ],
                                                         l_eff     = leff, 
                                                         delta_ell = dell, 
                                                         b_nmt     = bnmt,
                                                         fsky =fsky,
                                                         ) for i in np.arange(clf1_sim.shape[0]) ])    
        if params_general['verbose']:
            print('Simulated: (num of C_ell, num of ell):',CV_f1_sim.shape, '| delell: ',dell)
    if params_inout['output']['save_estimated']:
        CV_f1_est = np.array([ c4ft.cosmic_variance_ibin(clf1_ = clf1_est[i],
                                                         clf2_ = clf1_est[i],
                                                         clcx_ = clf1_est[i],
                                                         nlf1_ = Nl_hi ,
                                                         nlf2_ = Nl_hi ,
                                                         nlcx_ = Nl_hi,
                                                         fwhmf1_   = params_fields['hi'    ]['fwhm_arcmin'  ],
                                                         fwhmf2_   = params_fields['hi'    ]['fwhm_arcmin'  ],
                                                         fwhmcx_   = params_fields['hi'    ]['fwhm_arcmin'  ],
                                                         l_eff     = leff, 
                                                         delta_ell = dell, 
                                                         b_nmt     = bnmt,
                                                         fsky =fsky,
                                                         ) for i in np.arange(clf1_est.shape[0]) ])
        
        if params_general['verbose']:
            print('Simulated: (num of C_ell, num of ell):',CV_f1_sim.shape, '| delell: ',dell)        
    del bnmt,leff
#####################################################################################
### SAVING FILES
#####################################################################################
if params_general['verbose']:
    print('Preparing to save HIxHI APS...')
if params_inout['output']['save_hi'] and\
  (params_inout['output']['save_theory']+\
   params_inout['output']['save_simulated']+\
   params_inout['output']['save_estimated']):
    if params_general['verbose']:
        print('Saving HIxHI APS')    
    mask_ell = dcopy(params_APS['namaster']['mask_ell']) if params_APS['namaster']['ell_hi'] else dcopy(params_APS['mask_ell'])
    l_eff    = dcopy(params_APS['namaster']['leff'])[mask_ell] if params_APS['namaster']['ell_hi'] else dcopy(params_APS['leff'])[mask_ell]
    if params_inout['output']['save_theory']:    
        for i,ich in enumerate(params_inout['output']['auto_hi_bins']):
            CL_TH     = params_APS['namaster' ]['b'].bin_cell(clf1_th[ich,:]) if params_APS['namaster']['ell_hi'] else clf1_th[ich,:]
            CV_TH     = CV_f1_th[ich]    
            cl_eff    = CL_TH[mask_ell]
            errcl_eff = CV_TH[mask_ell]
            ###
            path_dir = params_inout['output']['pathout']
            file_str = params_inout['output']['filename_hi']#'cl_LOWZxLOWZ_{c}x{c}_delell{d}'
            ### 
            if params_inout['output']['save_hdf5']:
                try:
                    import h5py
                    filename_cl_h5 = os.path.join( path_dir, "".join(( file_str.format(b=bin_,c=ich, d=dell), '_theory.h5')) )
                    if params_general['verbose']: print('saving at... {}'.format(filename_cl_h5))
                    with h5py.File(filename_cl_h5, "w") as f:
                            f.create_dataset("l", data=l_eff)
                            f.create_dataset("cl", data=cl_eff)
                            f.create_dataset("errcl", data=errcl_eff)
                except:
                    pass  
        del CL_TH,CV_TH,cl_eff, errcl_eff, path_dir, file_str
            #if  params_inout['output']['save_txt']:
            #if  params_inout['output']['save_npz']:
    if params_inout['output']['save_simulated']: 
        for i,ich in enumerate(params_inout['output']['auto_hi_bins']):
            CL_SM     = params_APS['namaster' ]['b'].bin_cell(clf1_sim[ich,:]) if params_APS['namaster']['ell_hi'] else clf1_sim[ich,:]
            CV_SM     = CV_f1_sim[ich]    
            cl_eff    = CL_SM[mask_ell]
            errcl_eff = CV_SM[mask_ell]
            ###
            path_dir = params_inout['output']['pathout']
            file_str = params_inout['output']['filename_hi']
            ### 
            if params_inout['output']['save_hdf5']:
                try:
                    import h5py
                    filename_cl_h5 = os.path.join( path_dir, "".join(( file_str.format(b=bin_,c=ich, d=dell), '_simulated.h5')) )
                    if params_general['verbose']: print('saving at... {}'.format(filename_cl_h5))
                    with h5py.File(filename_cl_h5, "w") as f:
                            f.create_dataset("l", data=l_eff)
                            f.create_dataset("cl", data=cl_eff)
                            f.create_dataset("errcl", data=errcl_eff)
                except:
                    pass         
        del CL_SM,CV_SM,cl_eff, errcl_eff, path_dir, file_str
            #if  params_inout['output']['save_txt']:
            #if  params_inout['output']['save_npz']:        
    if params_inout['output']['save_estimated']:
        for i,ich in enumerate(params_inout['output']['auto_hi_bins']):
            CL_ES     = params_APS['namaster' ]['b'].bin_cell(clf1_est[ich,:]) if params_APS['namaster']['ell_hi'] else clf1_est[ich,:]
            CV_ES     = CV_f1_est[ich]    
            cl_eff    = CL_ES[mask_ell]
            errcl_eff = CV_ES[mask_ell]
            ###
            path_dir = params_inout['output']['pathout']
            file_str = params_inout['output']['filename_hi']
            ### 
            if params_inout['output']['save_hdf5']:
                try:
                    import h5py
                    filename_cl_h5 = os.path.join( path_dir, "".join(( file_str.format(b=bin_,c=ich, d=dell), '_estimated.h5')) )
                    if params_general['verbose']: print('saving at... {}'.format(filename_cl_h5))
                    with h5py.File(filename_cl_h5, "w") as f:
                            f.create_dataset("l", data=l_eff)
                            f.create_dataset("cl", data=cl_eff)
                            f.create_dataset("errcl", data=errcl_eff)
                except:
                    pass         
        del CL_ES,CV_ES,cl_eff, errcl_eff, path_dir, file_str  
            #if  params_inout['output']['save_txt']:
            #if  params_inout['output']['save_npz']:
    if params_general['verbose']:
        print('Saved.')
else:
    if params_general['verbose']:
        print('Skipping. Not saved.')
####################################################################################
if params_inout['output']['save_g']:
    fsky=1
    bin_ = params_inout['output']['auto_g_bin'][0] 
    bnmt = params_APS['namaster' ]['b'    ] if params_APS['namaster']['ell_g'] else None
    leff = params_APS['namaster' ]['leff' ] if params_APS['namaster']['ell_g'] else params_APS['leff']
    dell = params_APS['namaster' ]['del_l'] if params_APS['namaster']['ell_g'] else 1   
if params_general['verbose']:
    print('Preparing to save gxg APS...')

if params_inout['output']['save_g']:
    clf2_th  = np.loadtxt(params_inout['input']['filepath_g']).T[1:,:]
    if params_inout['output']['save_theory']:
        if params_general['verbose']:
            print('Calculating the galaxy cosmic noise (theoretical)')  
        CV_f2_th = np.array([ c4ft.cosmic_variance_ibin(clf1_     = clf2_th[0],
                                                        clf2_     = clf2_th[0],
                                                        clcx_     = clf2_th[0],
                                                        nlf1_     = Nl_g,
                                                        nlf2_     = Nl_g,
                                                        nlcx_     = Nl_g,
                                                        fwhmf1_   = params_fields['galaxy']['params_survey']['fwhm'],
                                                        fwhmf2_   = params_fields['galaxy']['params_survey']['fwhm'],
                                                        fwhmcx_   = params_fields['galaxy']['params_survey']['fwhm'],
                                                        l_eff     = leff, 
                                                        delta_ell = dell, 
                                                        b_nmt     = bnmt,
                                                        fsky      = fsky,
                                                       ) for i in np.arange(clf2_th.shape[0]) ])    
        if params_general['verbose']:
            print('Theory: (num of C_ell, num of ell):',CV_f2_th.shape, '| delell: ',dell)
        
    if params_inout['output']['save_simulated']:
        if params_general['verbose']:
            print('Calculating the galaxy cosmic noise (simulated)')   
            
        CV_f2_sim = np.array([ c4ft.cosmic_variance_ibin(clf1_     = clf2_sim[0],
                                                         clf2_     = clf2_sim[0],
                                                         clcx_     = clf2_sim[0],
                                                         nlf1_     = Nl_g ,
                                                         nlf2_     = Nl_g ,
                                                         nlcx_     = Nl_g,
                                                         fwhmf1_   = params_fields['galaxy']['params_survey']['fwhm'],
                                                         fwhmf2_   = params_fields['galaxy']['params_survey']['fwhm'],
                                                         fwhmcx_   = params_fields['galaxy']['params_survey']['fwhm'],
                                                         l_eff     = leff, 
                                                         delta_ell = dell, 
                                                         b_nmt     = bnmt,
                                                         fsky      = fsky,
                                                         ) for i in np.arange(clf2_sim.shape[0]) ])    
        if params_general['verbose']:
            print('Simulated: (num of C_ell, num of ell):',CV_f2_sim.shape, '| delell: ',dell)
    del bnmt,leff       
if params_inout['output']['save_g']  and \
   (params_inout['output']['save_theory'  ]+\
   params_inout['output']['save_simulated']): #do not make sense ESTIMATED here
    if params_general['verbose']:
        print('Saving gxg APS')    
    mask_ell = dcopy(params_APS['namaster']['mask_ell']) if params_APS['namaster']['ell_g'] else dcopy(params_APS['mask_ell'])
    l_eff    = dcopy(params_APS['namaster']['leff'])[mask_ell] if params_APS['namaster']['ell_g'] else dcopy(params_APS['leff'])[mask_ell]
    ###
    if params_inout['output']['save_theory']:   
        for i,ich in enumerate(params_inout['output']['auto_g_bin']):
            CL_TH     = params_APS['namaster' ]['b'].bin_cell(clf2_th[0,:]) if params_APS['namaster']['ell_g'] else clf2_th[0,:]
            CV_TH     = CV_f2_th[0]    
            cl_eff    = CL_TH[mask_ell]
            errcl_eff = CV_TH[mask_ell]
            ###
            path_dir = params_inout['output']['pathout']
            file_str = params_inout['output']['filename_g']#'cl_LSSTxLSST_{b}x{b}_delell{d}'
            ###
            if params_inout['output']['save_hdf5']:
                try:
                    import h5py
                    filename_cl_h5 = os.path.join( path_dir, "".join(( file_str.format(b=bin_,c=None, d=dell), '_theory.h5')) )
                    if params_general['verbose']: print('saving at... {}'.format(filename_cl_h5))
                    with h5py.File(filename_cl_h5, "w") as f:
                        f.create_dataset("l"    , data=l_eff    )
                        f.create_dataset("cl"   , data=cl_eff   )
                        f.create_dataset("errcl", data=errcl_eff)
                except:
                    pass     
        del CL_TH,CV_TH,cl_eff, errcl_eff, path_dir, file_str
            #if  params_inout['output']['save_txt']:
            #if  params_inout['output']['save_npz']:
    ###
    if params_inout['output']['save_simulated']:                         
        for i,ich in enumerate(params_inout['output']['auto_g_bin']):
            CL_SM     = params_APS['namaster' ]['b'].bin_cell(clf2_sim[0,:]) if params_APS['namaster']['ell_g'] else clf2_sim[0,:]
            CV_SM     = CV_f2_sim[0]    
            cl_eff    = CL_SM[mask_ell]
            errcl_eff = CV_SM[mask_ell]
            ###
            path_dir = params_inout['output']['pathout']
            file_str = params_inout['output']['filename_g']
            ### 
            if params_inout['output']['save_hdf5']:
                try:
                    import h5py
                    filename_cl_h5 = os.path.join( path_dir, "".join(( file_str.format(b=bin_,c=ich, d=dell), '_simulated.h5')) )
                    if params_general['verbose']: print('saving at... {}'.format(filename_cl_h5))
                    with h5py.File(filename_cl_h5, "w") as f:
                            f.create_dataset("l"    , data=l_eff    )
                            f.create_dataset("cl"   , data=cl_eff   )
                            f.create_dataset("errcl", data=errcl_eff)
                except:
                    pass         
        del CL_SM,CV_SM,cl_eff, errcl_eff, path_dir, file_str
            #if  params_inout['output']['save_txt']:
            #if  params_inout['output']['save_npz']:              
##########################################################
if params_general['verbose']:
    print('Preparing to save gxHI APS...')
if params_inout['output']['save_gxhi']:
    fsky=1
    bin_ = params_inout['output']['cross_g_bin'][0] 
    bnmt = params_APS['namaster' ]['b'    ] if params_APS['namaster']['ell_gxhi'] else None
    leff = params_APS['namaster' ]['leff' ] if params_APS['namaster']['ell_gxhi'] else params_APS['leff']
    dell = params_APS['namaster' ]['del_l'] if params_APS['namaster']['ell_gxhi'] else 1   

if params_inout['output']['save_gxhi']:
    if params_inout['output']['save_theory']:
        if params_general['verbose']:
            print('Calculating the g x HI cosmic noise (theoretical)')      
        CV_cx_th = np.array([ c4ft.cosmic_variance_ibin(clf1_     = clf1_th[i],
                                                        clf2_     = clf2_th[0],
                                                        clcx_     = clcx_th[i],
                                                        nlf1_     = Nl_hi     ,    
                                                        nlf2_     = Nl_g      ,
                                                        nlcx_     = 0         ,
                                                        fwhmf1_   = params_fields['hi'    ]['fwhm_arcmin'  ],
                                                        fwhmf2_   = params_fields['galaxy']['params_survey']['fwhm'],
                                                        fwhmcx_   = None,
                                                        l_eff     = leff, 
                                                        delta_ell = dell, 
                                                        b_nmt     = bnmt,
                                                        fsky      = fsky,
                                                        ) for i in np.arange(clcx_th.shape[0]) ])
        if params_general['verbose']:
            print('Theory: (num of C_ell, num of ell):',CV_cx_th.shape, '| delell: ',dell)
        
    if params_inout['output']['save_simulated']:
        if params_general['verbose']:
            print('Calculating the g x HI cosmic noise (simulated)')   
        CV_cx_sim = np.array([ c4ft.cosmic_variance_ibin(clf1_     = clf1_sim[i],
                                                         clf2_     = clf2_sim[0],
                                                         clcx_     = clcx_sim[i],
                                                         nlf1_     = Nl_hi      ,
                                                         nlf2_     = Nl_g       ,
                                                         nlcx_     = 0          ,
                                                        fwhmf1_    = params_fields['hi'    ]['fwhm_arcmin'  ],
                                                        fwhmf2_    = params_fields['galaxy']['params_survey']['fwhm'],
                                                        fwhmcx_    = None,
                                                         l_eff     = leff, 
                                                         delta_ell = dell, 
                                                         b_nmt     = bnmt,
                                                         fsky      = fsky,
                                                         ) for i in np.arange(clcx_sim.shape[0]) ])    
        if params_general['verbose']:
            print('Simulated: (num of C_ell, num of ell):',CV_cx_sim.shape, '| delell: ',dell)                
    if params_inout['output']['save_estimated']:
        if params_general['verbose']:
            print('Calculating the g x HI cosmic noise (estimated)')   
        CV_cx_est = np.array([ c4ft.cosmic_variance_ibin(clf1_     = clf1_est[i],
                                                         clf2_     = clf2_sim[0],
                                                         clcx_     = clcx_est[i],
                                                         nlf1_     = Nl_hi      ,
                                                         nlf2_     = Nl_g       ,
                                                         nlcx_     = 0          ,
                                                        fwhmf1_    = params_fields['hi'    ]['fwhm_arcmin'  ],
                                                        fwhmf2_    = params_fields['galaxy']['params_survey']['fwhm'],
                                                        fwhmcx_    = None,
                                                         l_eff     = leff, 
                                                         delta_ell = dell, 
                                                         b_nmt     = bnmt,
                                                         fsky      = fsky,
                                                         ) for i in np.arange(clcx_est.shape[0]) ])            

        if params_general['verbose']:
            print('Estimated: (num of C_ell, num of ell):',CV_cx_est.shape, '| delell: ',dell)        
    del bnmt,leff
       
    if params_general['verbose']:
        print('Saved.')
else:
    if params_general['verbose']:
        print('Skipping. Not saved.')
##############################################
if params_inout['output']['save_gxhi'] and\
  (params_inout['output']['save_theory'   ]+\
   params_inout['output']['save_simulated']+\
   params_inout['output']['save_estimated']):
    ###    
    mask_ell = dcopy(params_APS['namaster']['mask_ell']) if params_APS['namaster']['ell_gxhi'] else dcopy(params_APS['mask_ell'])
    l_eff    = dcopy(params_APS['namaster']['leff'])[mask_ell] if params_APS['namaster']['ell_gxhi'] else dcopy(params_APS['leff'])[mask_ell] 
    ###
    if params_inout['output']['save_theory']:   
        for i,ich in enumerate(params_inout['output']['cross_hi_bins']):
            CL_TH     = params_APS['namaster' ]['b'].bin_cell(clcx_th[ich,:]) if params_APS['namaster']['ell_gxhi'] else clcx_th[ich,:]
            CV_TH     = CV_cx_th[ich]    
            cl_eff    = CL_TH[mask_ell]
            errcl_eff = CV_TH[mask_ell]
            ###
            path_dir = params_inout['output']['pathout']
            file_str = params_inout['output']['filename_gxhi']#'cl_gcLSSTxLOWZ_{b}x{c}_40_0001_delell{d}'
            ###
            if params_inout['output']['save_hdf5']:
                try:
                    import h5py
                    filename_cl_h5 = os.path.join( path_dir, "".join(( file_str.format(b=bin_, c=ich, d=dell), '_theory.h5')) )
                    if params_general['verbose']: print('saving at... {}'.format(filename_cl_h5))
                    with h5py.File(filename_cl_h5, "w") as f:
                        f.create_dataset("l"    , data=l_eff    )
                        f.create_dataset("cl"   , data=cl_eff   )
                        f.create_dataset("errcl", data=errcl_eff)
                    del filename_cl_h5
                except:
                    pass    
        del CL_TH,CV_TH,cl_eff, errcl_eff, path_dir, file_str                        
            #if  params_inout['output']['save_txt']:
            #if  params_inout['output']['save_npz']:
    if params_inout['output']['save_simulated']:   
        for i,ich in enumerate(params_inout['output']['cross_hi_bins']):
            CL_SM     = params_APS['namaster' ]['b'].bin_cell(clcx_sim[ich,:]) if params_APS['namaster']['ell_gxhi'] else clcx_sim[ich,:]
            CV_SM     = CV_cx_sim[ich]    
            cl_eff    = CL_SM[mask_ell]
            errcl_eff = CV_SM[mask_ell]
            ###
            path_dir = params_inout['output']['pathout']
            file_str = params_inout['output']['filename_gxhi']#'cl_gcLSSTxLOWZ_{b}x{c}_40_0001_delell{d}'
            ###
            if params_inout['output']['save_hdf5']:
                try:
                    import h5py
                    filename_cl_h5 = os.path.join( path_dir, "".join(( file_str.format(b=bin_, c=ich, d=dell), '_simulated.h5')) )
                    if params_general['verbose']: print('saving at... {}'.format(filename_cl_h5))
                    with h5py.File(filename_cl_h5, "w") as f:
                        f.create_dataset("l"    , data=l_eff    )
                        f.create_dataset("cl"   , data=cl_eff   )
                        f.create_dataset("errcl", data=errcl_eff)
                    del filename_cl_h5
                except:
                    pass    
        del CL_SM,CV_SM,cl_eff, errcl_eff, path_dir, file_str     

    if params_inout['output']['save_simulated']:   
        for i,ich in enumerate(params_inout['output']['cross_hi_bins']):
            CL_ES     = params_APS['namaster' ]['b'].bin_cell(clcx_est[ich,:]) if params_APS['namaster']['ell_gxhi'] else clcx_est[ich,:]
            CV_ES     = CV_cx_est[ich]    
            cl_eff    = CL_ES[mask_ell]
            errcl_eff = CV_ES[mask_ell]
            ###
            path_dir = params_inout['output']['pathout']
            file_str = params_inout['output']['filename_gxhi']#'cl_gcLSSTxLOWZ_{b}x{c}_40_0001_delell{d}'
            ###
            if params_inout['output']['save_hdf5']:
                try:
                    import h5py
                    filename_cl_h5 = os.path.join( path_dir, "".join(( file_str.format(b=bin_, c=ich, d=dell), '_estimated.h5')) )
                    if params_general['verbose']: print('saving at... {}'.format(filename_cl_h5))
                    with h5py.File(filename_cl_h5, "w") as f:
                        f.create_dataset("l"    , data=l_eff    )
                        f.create_dataset("cl"   , data=cl_eff   )
                        f.create_dataset("errcl", data=errcl_eff)
                    del filename_cl_h5
                except:
                    pass    
        del CL_ES, CV_ES, cl_eff, errcl_eff, path_dir, file_str             
       
    if params_general['verbose']:
        print('Saved.')
else:
    if params_general['verbose']:
        print('Skipping. Not saved.')
#sys.exit(0)   
if params_general['verbose']:
    print('\n===========================================')
    print('Done.')
    print('===========================================\n')