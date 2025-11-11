#sys.path.insert(1, os.getcwd())
import os, sys, time
from copy import deepcopy as dcopy
import cross_functions_theory as cxft
import cross_functions_theory_4forecast as c4ft
import handling_data  as hdata
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
PATH        = "/data/AMARINS/LSSxHI-CODES/scripts"
INI         = "generating_APS_theory.ini"
name_params = os.path.join(PATH,INI)
config.read(name_params)
#General
verbose   = config.getboolean("General","verbose")
threshold = config.getfloat(  "General","threshold")
#Cosmology
ombh2   = config.getfloat("Cosmology","ombh2")
omch2   = config.getfloat("Cosmology","omch2")
H0      = config.getfloat("Cosmology","H0")
ns      = config.getfloat("Cosmology","ns")
As      = config.getfloat("Cosmology","As")
mnu     = config.getfloat("Cosmology","mnu")
w0      = config.getfloat("Cosmology","w0")
wa      = config.getfloat("Cosmology","wa")
DE_model= config.get(     "Cosmology","DE_model")
#CAMB
zmin_integration_camb   = config.getfloat(  "CAMB","zmin_integration_camb"  )
zmax_integration_camb   = config.getfloat(  "CAMB","zmax_integration_camb"  )
zsteps_integration_camb = config.getint(    "CAMB","zsteps_integration_camb")
kmax                    = config.getfloat(  "CAMB","kmax"                   )
minkh                   = config.getfloat(  "CAMB","minkh"                  )
maxkh                   = config.getfloat(  "CAMB","maxkh"                  )
npoints                 = config.getint(    "CAMB","npoints"                )
Pk_nonlinear            = config.getboolean("CAMB","Pk_nonlinear"           )
halofit_version         = config.get(       "CAMB","halofit_version"        )
accuracyBoost           = config.get(       "CAMB","accuracyBoost"          )
lsampleBoost            = config.get(       "CAMB","lsampleBoost"           )
laccuracyBoost          = config.get(       "CAMB","laccuracyBoost"         )
save_Pk                 = config.getboolean("CAMB","save_Pk"                )
save_Pk_at_z            = config.get(       "CAMB","save_Pk_at_z"           )
save_Pk_at_z            = np.array([float(item.strip()) for item in save_Pk_at_z.split(',')])

#load_Pk                = config.getboolean("CAMB","load_Pk"               )
#path_to_Pk             = config.get(       "CAMB","path_to_Pk"            )
#path_to_hubble_z      nprocs = comm.Get_size() = config.get(       "CAMB","path_to_hubble_z"      )
#path_to_comovel_dist_z = config.get(       "CAMB","path_to_comovel_dist_z")
############
#APS
field_1             = config.get(        "APS","field_1"            )
field_2             = config.get(        "APS","field_2"            )
l_min               = config.getint(     "APS","l_min"              )
l_max               = config.getint(     "APS","l_max"              )
generate_cl_field_1 = config.getboolean( "APS","generate_cl_field_1")
generate_cl_field_2 = config.getboolean( "APS","generate_cl_field_2")
generate_cl_cross   = config.getboolean( "APS","generate_cl_cross"  )
pathout             = config.get(        "APS","pathout"            )
prefix              = config.get(        "APS","prefix"             )
suffix              = config.get(        "APS","suffix"             )
############
#SHEAR-survey
path_to_shear_params_file = config.get(        "SHEAR","path_to_shear_params_file")
survey_shear              = config.get(        "SHEAR","survey_shear"             )
dataset_shear             = config.get(        "SHEAR","dataset_shear"            )
binning_shear             = config.getboolean( "SHEAR","binning_shear"            )
parameters_shear_set      = config.get(        "SHEAR","parameters_shear_set"     )
photometric_shear_lkl     = config.get(        "SHEAR","photometric_shear_lkl"    )
bin_to_use_shear          = config.getint(     "SHEAR","bin_to_use_shear"         )
include_ia                = config.getboolean( "SHEAR","include_ia"               )
L_model                   = config.get(        "SHEAR","L_model"                  )
L_value                   = config.getfloat(   "SHEAR","L_value"                  )
path_to_L_file            = config.get(        "SHEAR","path_to_L_file"           )
F_model                   = config.get(        "SHEAR","F_model"                  )
F_value                   = config.getfloat(   "SHEAR","F_value"                  )
path_to_F_file            = config.get(        "SHEAR","path_to_F_file"           )
beta_IA                   = config.getfloat(   "SHEAR","beta_IA"                  )
eta_IA                    = config.getfloat(   "SHEAR","eta_IA"                   )
#z_shear_min     = config.getfloat(   "SHEAR","z_shear_min" )
#z_shear_max     = config.getfloat(   "SHEAR","z_shear_max" )
#Changing the yaml file
#None values will not change the standard value
nbins_shear               = config.getint(     "SHEAR","nbins_shear"              )
sigmaz_shear              = config.getfloat(   "SHEAR","sigmaz_shear"             )
#zbins_edges_shear         = config.get(        "SHEAR","zbins_edges_shear"        )

############
#COUNT-survey
path_to_count_params_file = config.get(        "COUNT","path_to_count_params_file")
survey_count              = config.get(        "COUNT","survey_count"             )
dataset_count             = config.get(        "COUNT","dataset_count"            )
mode_count                = config.get(        "COUNT","mode_count"               )
bias_model                = config.get(        "COUNT","bias_model"               )
bias_value                = config.getfloat(   "COUNT","bias_value"               )
path_to_bias_file         = config.get(        "COUNT","path_to_bias_file"        )
binning_count             = config.getboolean( "COUNT","binning_count"            )
parameters_count_set      = config.get(        "COUNT","parameters_count_set"     )
photometric_count_lkl     = config.get(        "COUNT","photometric_count_lkl"    )
bin_to_use_count          = config.getint(     "COUNT","bin_to_use_count"         )
include_rsd               = config.getboolean( "COUNT","include_rsd"              )
include_mb                = config.getboolean( "COUNT","include_mb"               )
mu_rsd                    = config.getfloat(   "COUNT","mu_rsd"                   )
slope_s_model             = config.get(        "COUNT","slope_s_model"            )
slope_s_value             = config.getfloat(   "COUNT","slope_s_value"            )
path_to_slope_s_file      = config.get(        "COUNT","path_to_slope_s_file"     )
#Changing the yaml file
#None values will not change the standard value
nbins_count               = config.getint(     "COUNT","nbins_count"              )
sigmaz_count              = config.getfloat(   "COUNT","sigmaz_count"             )
#zbins_edges_count         = config.get(        "COUNT","zbins_edges_count"        )  
############
#HI-survey
#binning_used       = config.get("survey","binning_used")
freq_min        = config.getfloat("HI","freq_min")
freq_max        = config.getfloat("HI","freq_max")
nchannels       = config.getint(  "HI","nchannels")
#window_function_HI = config.get(     "HI","window_function_HI")
#21cm
omegaHI_model = config.get(       "21cm", "OmegaHI_model")
biasHI_model  = config.get(       "21cm", "biasHI_model" )
fit_omegaHI   = config.getboolean("21cm", "fit_omegaHI"  )
fit_biasHI    = config.getboolean("21cm", "fit_biasHI"   )
fit_omegaHI   = False
fit_biasHI    = False
omegaHI       = config.get("21cm","omegaHI")
biasHI        = config.get("21cm","biasHI")

###############################################################################
# You can modify any options in the parameters.ini file by the command terminal
###############################################################################
parser = argparse.ArgumentParser(description='Modify by the command terminal parameters in {} file'.format(INI))

#General
parser.add_argument('--verbose'   , action = 'store', dest = 'verbose'  , default = verbose  , help = '')
parser.add_argument('--threshold' , action = 'store', dest = 'threshold', default = threshold, help = '')

#Cosmology
parser.add_argument('--ombh2'   , action = 'store', dest = 'ombh2'   , default = ombh2   , help = '')
parser.add_argument('--omch2'   , action = 'store', dest = 'omch2'   , default = omch2   , help = '')
parser.add_argument('--H0'      , action = 'store', dest = 'H0'      , default = H0      , help = '')
parser.add_argument('--ns'      , action = 'store', dest = 'ns'      , default = ns      , help = '')
parser.add_argument('--As'      , action = 'store', dest = 'As'      , default = As      , help = '')
parser.add_argument('--mnu'     , action = 'store', dest = 'mnu'     , default = mnu     , help = '')
parser.add_argument('--w0'      , action = 'store', dest = 'w0'      , default = w0      , help = '')
parser.add_argument('--wa'      , action = 'store', dest = 'wa'      , default = wa      , help = '')
parser.add_argument('--DE_model', action = 'store', dest = 'DE_model', default = DE_model, help = '')

#CAMB
parser.add_argument('--zmin_integration_camb'  , action = 'store', dest = 'zmin_integration_camb'  , default = zmin_integration_camb  , help = '')
parser.add_argument('--zmax_integration_camb'  , action = 'store', dest = 'zmax_integration_camb'  , default = zmax_integration_camb  , help = '')
parser.add_argument('--zsteps_integration_camb', action = 'store', dest = 'zsteps_integration_camb', default = zsteps_integration_camb, help = '')
parser.add_argument('--kmax'           , action = 'store', dest = 'kmax'           , default = kmax           , help = '')
parser.add_argument('--minkh'          , action = 'store', dest = 'minkh'          , default = minkh          , help = '')
parser.add_argument('--maxkh'          , action = 'store', dest = 'maxkh'          , default = maxkh          , help = '')
parser.add_argument('--npoints'        , action = 'store', dest = 'npoints'        , default = npoints        , help = '')
parser.add_argument('--Pk_nonlinear'   , action = 'store', dest = 'Pk_nonlinear'   , default = Pk_nonlinear   , help = '')
parser.add_argument('--halofit_version', action = 'store', dest = 'halofit_version', default = halofit_version, help = '')
parser.add_argument('--accuracyBoost'  , action = 'store', dest = 'accuracyBoost'  , default = accuracyBoost  , help = '')
parser.add_argument('--lsampleBoost'   , action = 'store', dest = 'lsampleBoost'   , default = lsampleBoost   , help = '')
parser.add_argument('--laccuracyBoost' , action = 'store', dest = 'laccuracyBoost' , default = laccuracyBoost , help = '')
parser.add_argument('--save_Pk'        , action = 'store', dest = 'save_Pk'        , default = save_Pk        , help = '')
#parser.add_argument('--save_Pk_at_z'   , action = 'store', dest = 'save_Pk_at_z'   , default = save_Pk_at_z   , help = '')
parser.add_argument('--save_Pk_at_z'   , action = 'store', dest = 'save_Pk_at_z'   , default = save_Pk_at_z   , nargs='+', type=float)#,  help='List of numbers (space separated)')

#parser.add_argument('--load_Pk'               , action = 'store', dest = 'load_Pk'               , default = load_Pk               , help = '')
#parser.add_argument('--path_to_Pk'            , action = 'store', dest = 'path_to_Pk'            , default = path_to_Pk            , help = '')
#parser.add_argument('--path_to_hubble_z'      , action = 'store', dest = 'path_to_hubble_z'      , default = path_to_hubble_z      , help = '')
#parser.add_argument('--path_to_comovel_dist_z', action = 'store', dest = 'path_to_comovel_dist_z', default = path_to_comovel_dist_z, help = '')

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

#SHEARS
parser.add_argument('--path_to_shear_params_file', action = 'store', dest = 'path_to_shear_params_file', default = path_to_shear_params_file, help = '')
parser.add_argument('--survey_shear'             , action = 'store', dest = 'survey_shear'             , default = survey_shear             , help = '')
parser.add_argument('--dataset_shear'            , action = 'store', dest = 'dataset_shear'            , default = dataset_shear            , help = '')
parser.add_argument('--binning_shear'            , action = 'store', dest = 'binning_shear'            , default = binning_shear            , help = '')
parser.add_argument('--parameters_shear_set'     , action = 'store', dest = 'parameters_shear_set'     , default = parameters_shear_set     , help = '')
parser.add_argument('--photometric_shear_lkl'    , action = 'store', dest = 'photometric_shear_lkl'    , default = photometric_shear_lkl    , help = '')
parser.add_argument('--bin_to_use_shear'         , action = 'store', dest = 'bin_to_use_shear'         , default = bin_to_use_shear         , help = '')
parser.add_argument('--include_ia'               , action = 'store', dest = 'include_ia'               , default = include_ia               , help = '')
parser.add_argument('--L_model'                  , action = 'store', dest = 'L_model'                  , default = L_model                  , help = '')
parser.add_argument('--L_value'                  , action = 'store', dest = 'L_value'                  , default = L_value                  , help = '')
parser.add_argument('--F_model'                  , action = 'store', dest = 'F_model'                  , default = F_model                  , help = '')
parser.add_argument('--F_value'                  , action = 'store', dest = 'F_value'                  , default = F_value                  , help = '')
parser.add_argument('--path_to_L_file'           , action = 'store', dest = 'path_to_L_file'           , default = path_to_L_file           , help = '')
parser.add_argument('--path_to_F_file'           , action = 'store', dest = 'path_to_F_file'           , default = path_to_F_file           , help = '')
parser.add_argument('--beta_IA'                  , action = 'store', dest = 'beta_IA'                  , default = beta_IA                  , help = '')
parser.add_argument('--eta_IA'                   , action = 'store', dest = 'eta_IA'                   , default = eta_IA                   , help = '')
#parser.add_argument('--z_shear_min'     , action = 'store', dest = 'z_shear_min'     , default = z_shear_min    , help = '')
#parser.add_argument('--z_shear_max'     , action = 'store', dest = 'z_shear_max'     , default = z_shear_max    , help = '')
parser.add_argument('--nbins_shear'              , action = 'store', dest = 'nbins_shear'              , default = nbins_shear              , help = '')
parser.add_argument('--sigmaz_shear'             , action = 'store', dest = 'sigmaz_shear'             , default = sigmaz_shear             , help = '')
#parser.add_argument('--zbins_edges_shear'        , action = 'store', dest = 'zbins_edges_shear'        , default = zbins_edges_shear        , help = '')

#COUNT
parser.add_argument('--path_to_count_params_file', action = 'store', dest = 'path_to_count_params_file', default = path_to_count_params_file, help = '')
parser.add_argument('--survey_count'             , action = 'store', dest = 'survey_count'             , default = survey_count             , help = '')
parser.add_argument('--dataset_count'            , action = 'store', dest = 'dataset_count'            , default = dataset_count            , help = '')
parser.add_argument('--mode_count'               , action = 'store', dest = 'mode_count'               , default = mode_count               , help = '')
parser.add_argument('--bias_model'               , action = 'store', dest = 'bias_model'               , default = bias_model               , help = '')
parser.add_argument('--bias_value'               , action = 'store', dest = 'bias_value'               , default = bias_value               , help = '')
parser.add_argument('--path_to_bias_file'        , action = 'store', dest = 'path_to_bias_file'        , default = path_to_bias_file        , help = '')
parser.add_argument('--binning_count'            , action = 'store', dest = 'binning_count'            , default = binning_count            , help = '')
parser.add_argument('--parameters_count_set'     , action = 'store', dest = 'parameters_count_set'     , default = parameters_count_set     , help = '')
parser.add_argument('--photometric_count_lkl'    , action = 'store', dest = 'photometric_count_lkl'    , default = photometric_count_lkl    , help = '')
parser.add_argument('--bin_to_use_count'         , action = 'store', dest = 'bin_to_use_count'         , default = bin_to_use_count         , help = '')
parser.add_argument('--include_rsd'              , action = 'store', dest = 'include_rsd'              , default = include_rsd              , help = '')
parser.add_argument('--include_mb'               , action = 'store', dest = 'include_mb'               , default = include_mb               , help = '')
parser.add_argument('--mu_rsd'                   , action = 'store', dest = 'mu_rsd'                   , default = mu_rsd                   , help = '')
parser.add_argument('--slope_s_model'            , action = 'store', dest = 'slope_s_model'            , default = slope_s_model            , help = '')
parser.add_argument('--slope_s_value'            , action = 'store', dest = 'slope_s_value'            , default = slope_s_value            , help = '')
parser.add_argument('--path_to_slope_s_file'     , action = 'store', dest = 'path_to_slope_s_file'     , default = path_to_slope_s_file     , help = '')
parser.add_argument('--nbins_count'              , action = 'store', dest = 'nbins_count'              , default = nbins_count              , help = '')
parser.add_argument('--sigmaz_count'             , action = 'store', dest = 'sigmaz_count'             , default = sigmaz_count             , help = '')
#parser.add_argument('--zbins_edges_count'        , action = 'store', dest = 'zbins_edges_count'        , default = zbins_edges_count        , help = '')
#HI
parser.add_argument('--freq_min'          , action = 'store', dest = 'freq_min'          , default = freq_min          , help = '')
parser.add_argument('--freq_max'          , action = 'store', dest = 'freq_max'          , default = freq_max          , help = '')
parser.add_argument('--nchannels'         , action = 'store', dest = 'nchannels'         , default = nchannels         , help = '')
#parser.add_argument('--window_function_HI', action = 'store', dest = 'window_function_HI', default = window_function_HI, help = '')

#21cm
parser.add_argument('--omegaHI_model', action = 'store', dest = 'omegaHI_model', default = omegaHI_model, help = '')
parser.add_argument('--biasHI_model' , action = 'store', dest = 'biasHI_model' , default = biasHI_model , help = '')
parser.add_argument('--fit_omegaHI'  , action = 'store', dest = 'fit_omegaHI'  , default = fit_omegaHI  , help = '')
parser.add_argument('--fit_biasHI'   , action = 'store', dest = 'fit_biasHI'   , default = fit_biasHI   , help = '')
parser.add_argument('--omegaHI'      , action = 'store', dest = 'omegaHI'      , default = omegaHI      , help = '')
parser.add_argument('--biasHI'       , action = 'store', dest = 'biasHI'       , default = biasHI       , help = '')

arguments = parser.parse_args()
###############################################################################
# Variables
###############################################################################
#GENERAL
verbose   = bool( arguments.verbose)
threshold = float(arguments.threshold)
#COSMOLOGY
H0       = float(arguments.H0)
ombh2    = float(arguments.ombh2)
omch2    = float(arguments.omch2)
ns       = float(arguments.ns)
As       = float(arguments.As)
mnu      = float(arguments.mnu)
w0       = float(arguments.w0)
wa       = float(arguments.wa)
DE_model = str(  arguments.DE_model)
#CAMB
zmin_integration_camb   = float(arguments.zmin_integration_camb)
zmax_integration_camb   = float(arguments.zmax_integration_camb)
zsteps_integration_camb = int(arguments.zsteps_integration_camb)
kmax                    = float(arguments.kmax)
minkh                   = float(arguments.minkh)
maxkh                   = float(arguments.maxkh)
npoints                 = int(arguments.npoints)
Pk_nonlinear            = bool(arguments.Pk_nonlinear)
halofit_version         = str(arguments.halofit_version)
accuracyBoost           = float(arguments.accuracyBoost)
lsampleBoost            = float(arguments.lsampleBoost)
laccuracyBoost          = float(arguments.laccuracyBoost)
save_Pk                 = bool(arguments.save_Pk)
save_Pk_at_z            = np.asarray(arguments.save_Pk_at_z)

#load_Pk                 = bool(arguments.load_Pk)
#path_to_Pk              = str(arguments.path_to_Pk)
#path_to_hubble_z        = str(arguments.path_to_hubble_z)
#path_to_comovel_dist_z  = str(arguments.path_to_comovel_dist_z)

#APS
field_1             = str(arguments.field_1)
field_2             = str(arguments.field_2)
l_min               = float(arguments.l_min)
l_max               = float(arguments.l_max)
generate_cl_field_1 = bool(arguments.generate_cl_field_1)
generate_cl_field_2 = bool(arguments.generate_cl_field_2)
generate_cl_cross   = bool(arguments.generate_cl_cross)
pathout             = str(arguments.pathout)
prefix              = str(arguments.prefix)
suffix              = str(arguments.suffix)

#SHEAR
path_to_shear_params_file = str(arguments.path_to_shear_params_file)
survey_shear              = str(arguments.survey_shear)
dataset_shear             = str(arguments.dataset_shear)
binning_shear             = bool(arguments.binning_shear)
parameters_shear_set      = str(arguments.parameters_shear_set)
photometric_shear_lkl     = str(arguments.photometric_shear_lkl)
bin_to_use_shear          = int(arguments.bin_to_use_shear)
include_ia                = bool(arguments.include_ia)
L_model                   = str(arguments.L_model)
L_value                   = float(arguments.L_value)
path_to_L_file            = str(arguments.path_to_L_file)
F_model                   = str(arguments.F_model)
F_value                   = float(arguments.F_value)
path_to_F_file            = str(arguments.path_to_F_file)
beta_IA                   = float(arguments.beta_IA)
eta_IA                    = float(arguments.eta_IA)
#z_shear_min     = float(arguments.z_shear_min)
#z_shear_max     = float(arguments.z_shear_max)
nbins_shear               = int(arguments.nbins_shear)
sigmaz_shear              = float(arguments.sigmaz_shear)
#zbins_edges_shear

#COUNT
path_to_count_params_file = str(arguments.path_to_count_params_file)
survey_count              = str(arguments.survey_count)
dataset_count             = str(arguments.dataset_count)
mode_count                = str(arguments.mode_count)
bias_model                = str(arguments.bias_model)
bias_value                = float(arguments.bias_value)
path_to_bias_file         = str(arguments.path_to_bias_file)
binning_count             = bool(arguments.binning_count)
parameters_count_set      = str(arguments.parameters_count_set)
photometric_count_lkl     = str(arguments.photometric_count_lkl)
bin_to_use_count          = int(arguments.bin_to_use_count)
include_rsd               = bool(arguments.include_rsd)
include_mb                = bool(arguments.include_mb)
mu_rsd                    = float(arguments.mu_rsd)
slope_s_model             = str(arguments.slope_s_model)
slope_s_value             = float(arguments.slope_s_value)
path_to_slope_s_file      = str(arguments.path_to_slope_s_file)
nbins_count               = int(arguments.nbins_count)
sigmaz_count              = float(arguments.sigmaz_count)
#zbins_edges_count =

#HI
freq_min  = float(arguments.freq_min)
freq_max  = float(arguments.freq_max)
nchannels = int(arguments.nchannels)

#21cm
omegaHI_model = str(arguments.omegaHI_model)
biasHI_model  = str(arguments.biasHI_model)
fit_omegaHI   = bool(arguments.fit_omegaHI)
fit_biasHI    = bool(arguments.fit_biasHI)
omegaHI       = str(arguments.omegaHI)
biasHI        = str(arguments.biasHI)
############################################################################################
#####################################################################################################
if omegaHI_model=="":
    omegaHI_model=None
if biasHI_model=="":
    biasHI_model=None
    
if omegaHI=="":
    omegaHI=None
if biasHI=="":
    biasHI=None
    
#if path_to_Pk=="":
#    path_to_Pk=None
#if path_to_hubble_z=="":
#    path_to_hubble_z=None    
#if path_to_comovel_dist_z=="":
#    path_to_comovel_dist_z=None        
    
if prefix=="":
    prefix=None
if suffix=="":
    suffix=None
    
if path_to_bias_file=="":
    path_to_bias_file=None
if path_to_slope_s_file=="":
    path_to_slope_s_file=None    
if slope_s_model=="":
    slope_s_model=None    
#    
if nbins_shear<=0:
    nbins_shear=None
else:
    nbins_shear=int(nbins_shear)    
if sigmaz_shear<=0:
    sigmaz_shear=None  
else:
    sigmaz_shear=float(sigmaz_shear)       
    
if nbins_count<=0:
    nbins_count=None
else:
    nbins_count=int(nbins_count)
if sigmaz_count<=0:
    sigmaz_count=None  
else:
    sigmaz_count=float(sigmaz_count)    
#print(suffix, type(suffix), len(suffix))
#####################################################################################################################################################################
#####################################################################################################################################################################
params_general = {'verbose'  :verbose,
                  'threshold':threshold}
params = { 'H0'     : H0,
           'ombh2'  : ombh2,
           'omch2'  : omch2,
           'ns'     : ns,
           'As'     : As,
           'mnu'    : mnu,
           'w0'     : w0,
           'wa'     : wa,
           'DE_model': DE_model,
           'z'      : np.logspace(np.log10(zmin_integration_camb),np.log10(zmax_integration_camb), zsteps_integration_camb),
           'kmax'   : kmax,
           'minkh'  : minkh,
           'maxkh'  : maxkh,
           'npoints': npoints,
           'accuracyBoost'  : accuracyBoost,
           'lsampleBoost'   : lsampleBoost,
           'laccuracyBoost' : laccuracyBoost,
           'Pk_nonlinear'   : Pk_nonlinear,
           'halofit_version': halofit_version,
           'omegaHI_model'  : omegaHI_model,
           'biasHI_model'   : biasHI_model,
           'save_Pk'        : save_Pk,
           'save_Pk_at_z'   : save_Pk_at_z,
           #'load_Pk'               : load_Pk,
           #'path_to_Pk'            : path_to_Pk,
           #'path_to_hubble_z'      : path_to_hubble_z,
           #'path_to_comovel_dist_z': path_to_comovel_dist_z,
        }          

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
                }

params_cmbwl_survey = {'field':'cmbwl'}

params_shear_survey = {'field':'shear', 'path_to_params_file':path_to_shear_params_file,
                       "survey":survey_shear, "dataset":dataset_shear, 'binning':binning_shear,
                       "parameters_set":parameters_shear_set, "photometric_lkl":photometric_shear_lkl, "bin_to_use":bin_to_use_shear,
                       "include_ia":include_ia,
                       'L_model':L_model, 'L_value':L_value, 'path_to_L_file':path_to_L_file,
                       'F_model':F_model, 'F_value':F_value, 'path_to_F_file':path_to_F_file,
                       'beta_IA':beta_IA, 'eta_IA':eta_IA,
                       #"z_min":z_shear_min, "z_max":z_shear_max, 
                       #'nbins':nbins_shear, 'sigmaz':sigmaz_shear,                        
                      }

params_count_survey = {'field':'density','path_to_params_file':path_to_count_params_file,
                       "survey":survey_count, "dataset":dataset_count, 'mode_count':mode_count,
                       'bias_model':bias_model, 'bias_value':bias_value,  'path_to_bias_file':path_to_bias_file,
                       'binning':binning_count,
                       "parameters_set":parameters_count_set, "photometric_lkl":photometric_count_lkl, "bin_to_use":bin_to_use_count,
                       "include_rsd":include_rsd, 'include_mb':include_mb,
                       'mu_rsd':mu_rsd, 
                       'slope_s_model':slope_s_model, 'slope_s_value':slope_s_value, 'path_to_slope_s_file':path_to_slope_s_file,
                       #'nbins':nbins_count, 'sigmaz':sigmaz_count, 
                      }

params_hi_survey    = {'field':'hi',
                       "freq_min"     :freq_min,          "freq_max":freq_max,    "nchannels":   nchannels, 
                       'omegaHI_model':omegaHI_model, 'biasHI_model':biasHI_model, 
                       "fit_omegaHI"  :fit_omegaHI,     'fit_biasHI':fit_biasHI,
                       'omegaHI'      :omegaHI,             'biasHI':biasHI}

params_APS['l'] = np.arange(params_APS['l_min'], params_APS['l_max']+1, 1)	
#binning=True
########################################################################################################################################
del verbose, threshold
del H0, ombh2, omch2, As, mnu, w0, wa, DE_model, ns, kmax, minkh, maxkh, npoints, Pk_nonlinear, halofit_version, save_Pk, save_Pk_at_z
del accuracyBoost, lsampleBoost, laccuracyBoost
del field_1, field_2, l_min, l_max, generate_cl_field_1, generate_cl_field_2, generate_cl_cross, pathout, prefix, suffix, 
del path_to_shear_params_file, survey_shear, dataset_shear, binning_shear, parameters_shear_set, photometric_shear_lkl, bin_to_use_shear, 
del include_ia, L_model, L_value, path_to_L_file, F_model, F_value, path_to_F_file, beta_IA, eta_IA, # nbins_shear,sigmaz_shear
del path_to_count_params_file, survey_count, dataset_count, binning_count, bias_model, bias_value, path_to_bias_file, parameters_count_set, 
del photometric_count_lkl, mode_count, #nbins_count,sigmaz_count
del include_rsd, include_mb, mu_rsd, slope_s_model, slope_s_value, path_to_slope_s_file
del freq_min, freq_max, nchannels, omegaHI_model, biasHI_model, fit_biasHI, fit_omegaHI, omegaHI, biasHI
########################################################################################################################################
#### 3D MATTER POWER SPECTRUM from CAMB ################################################################################################
if params_general['verbose']:
    print('\n===========================================')
    print('Starting the code')
    print('===========================================\n')
    #
    print('Calculating the matter Pk...')

pars = camb.set_params(halofit_version=params['halofit_version'],
                       w=params['w0'], wa=params['wa'], 
                       dark_energy_model=params['DE_model'],
                       H0=params['H0'], ombh2=params['ombh2'], 
                       omch2=params['omch2'], mnu=params['mnu'],
                       #WantCls=0,
                       #k_per_logint=0,                       
                      )#, HMCode_A_baryon=3.13, HMCode_eta_baryon=0.603, HMCode_logT_AGN=7.8)    

pars.set_accuracy(AccuracyBoost=params['accuracyBoost'], lSampleBoost=params['lsampleBoost'], lAccuracyBoost=params['laccuracyBoost'])
pars.InitPower.set_params(ns=params['ns'], As=params['As'])
pars.set_matter_power(redshifts=params['z'], kmax=params['kmax'], silent=True)
results   = camb.get_results(pars)
pk_interp = camb.get_matter_power_interpolator(pars, 
                                               nonlinear=params['Pk_nonlinear'], 
                                               hubble_units=False, k_hunit=False, 
                                               kmax=params['kmax'], 
                                               zmax=params['z'][-1])

results.z_input = params['z']
#####################################################################################################
####FIELD1
if params_general['verbose']:
    print('Calculating the LSS kernels for:')
    
if params_APS['field_1'].lower() == 'hi':
    params_APS['field_1_params'] = dcopy(params_hi_survey)
    params_APS['field_1_params'].update(hdata.nu_bins_vector(numin_ =params_APS['field_1_params']['freq_min'], 
                                        numax_ =params_APS['field_1_params']['freq_max'], 
                                        nbands_=params_APS['field_1_params']['nchannels'])
                                       )
    params_APS['field_1_params']['z'] = dcopy(np.flip(params_APS['field_1_params']['z'])) #<------- building F1 wth z-increasing
    #W_f1_func  = lambda x: cxft.kernel_hi(camb_params=pars, camb_results=results, 
    #                                 z=x, zrange=params_APS['field_1_params']['z'], 
    #                                 config_dict=params_APS['field_1_params'] , replace_inf2nan=True)
    W_f1_func = c4ft.kernel_hi_vec(
                #zout=zv, 
                camb_params=pars, 
                camb_results=results,
                zrange=params_APS['field_1_params']['z'], 
                replace_inf2nan=True, 
                use_input=False,
                config_dict=params_APS["field_1_params"], 
                params_interp=None, 
                get_interp_function=1
            ) 
    
    params_APS['field_1_params']['header'] = "[1] l, [>2] {}-{} Cl [multipoles, redshift/frequency bins]".format(params_APS['field_1'].upper(), params_APS['field_1'].upper())
    del params_APS['field_1_params']['min'], params_APS['field_1_params']['max']    
    if params_general['verbose']: print('Field 1: HI')
    del params_hi_survey
    #
elif params_APS['field_1'].lower() in ['cmbwl', 'wlcmb', 'convergence', 'cmb', 'kcmb', 'cmbk']:    
    params_APS['field_1_params'] = dcopy(params_cmbwl_survey)    
    params_APS['field_1_params']['z_min'] = 0
    params_APS['field_1_params']['z_max'] = cxft.chi_vec(results)['z_star']
    W_f1_func  = cxft.kernel_cmb_interp(camb_params=pars, camb_results=results)
    params_APS['field_1_params']['header'] = "[1] l, [2] {f1}-{f1} Cl [multipoles, redshift/frequency bin] --survey: CMB".format(f1=params_APS['field_1'].upper())
    if params_general['verbose']: print('Field 1: CMB WL')    
    del params_cmbwl_survey
    #
elif params_APS['field_1'].lower() in ['lensing','lens', 'shear']:       
    params_APS['field_1_params'] = dcopy(params_shear_survey)  
    params_APS['field_1_params'].update(cxft.load_survey_parameters(survey_name = params_APS['field_1_params']['survey'], 
                                        data_release= params_APS['field_1_params']['dataset'],
                                        type_field='lensing',
                                        path_file=params_APS['field_1_params']['path_to_params_file'])
                                       )    
    if type(None)!=type(nbins_shear):
        params_APS['field_1_params']['nbins']=dcopy(nbins_shear)  
    if type(None)!=type(sigmaz_shear):
        if 'LSST' in params_APS['field_1_params']['survey'].upper():
            params_APS['field_1_params']['sigma_z']=dcopy(sigmaz_shear)
        elif 'EUCLID' in params_APS['field_1_params']['survey'].upper():
            params_APS['field_1_params']['sigma_b']=dcopy(sigmaz_shear)
            params_APS['field_1_params']['sigma_0']=dcopy(sigmaz_shear)
        else:
            pass             
    W_f1_i = cxft.kernel_galaxy_lensing_i_dict(ibin=params_APS['field_1_params']['bin_to_use'],
                                              camb_params=pars, camb_results=results, binning=params_shear_survey['binning'], 
                                              pars=params_APS['field_1_params'], verbose=0)
    W_f1_func = lambda x: W_f1_func['W_interp'](x) #if (x>=params_APS['field_1_params']['z_min'])*(x<=params_APS['field_1_params']['z_max']) else 0    
    params_APS['field_1_params']['header'] = "[1] l, [>2] {f1}-{f1} Cl [multipoles, redshift/frequency bin] --survey: {s} | --dataset: {d} | bin: {b}".format(f1=params_APS['field_1'].upper(),
                                                                                                                                                               s=params_APS['field_1_params']['survey'].upper(), 
                                                                                                                                                               d=params_APS['field_1_params']['dataset'].upper(),
                                                                                                                                                               b=params_APS['field_1_params']['bin_to_use']
                                                                                                                                                               )    
    if params_general['verbose']: print('Field 1: SHEAR')        
    del params_shear_survey

elif params_APS['field_1'].lower() in ['galaxy', 'density', 'dens', 'source', 'overdensity', 'count', 'counting']:           
    params_APS['field_1_params'] = dcopy(params_count_survey)  
    params_APS['field_1_params'].update(cxft.load_survey_parameters(survey_name = params_APS['field_1_params']['survey'], 
                                            data_release= params_APS['field_1_params']['dataset'], 
                                            type_field='density',
                                            path_file=params_APS['field_1_params']['path_to_params_file'])
                                       )       
    
    #print(params_APS['field_1_params']['include_mb' ])
    #print(params_APS['field_1_params']['include_rsd'])
    if type(None)!=type(nbins_count):
        params_APS['field_1_params']['nbins']=dcopy(nbins_count) 
    if type(None)!=type(sigmaz_count):
        if 'LSST' in params_APS['field_1_params']['survey'].upper():
            params_APS['field_1_params']['sigma_z']=dcopy(sigmaz_count)
        elif 'EUCLID' in params_APS['field_1_params']['survey'].upper():
            params_APS['field_1_params']['sigma_b']=dcopy(sigmaz_count)
            params_APS['field_1_params']['sigma_0']=dcopy(sigmaz_count)
        else:
            pass        
    W_f1_i = cxft.kernel_galaxy_clustering_i_dict(camb_params=pars, camb_results=results, pars=params_APS['field_1_params'],
                                                  ibin = params_APS['field_1_params']['bin_to_use'],
                                                  include_rsd=params_APS['field_1_params']['include_rsd'], 
                                                  include_mb = params_APS['field_1_params']['include_mb'] , growth_rate_f=None, 
                                                  bias_model = params_APS['field_1_params']['bias_model'], 
                                                  bias_value = params_APS['field_1_params']['bias_value'],
                                                  path_to_bias_file=params_APS['field_1_params']['path_to_bias_file'],
                                                  mb_model   = params_APS['field_1_params']['slope_s_model'], 
                                                  mb_alpha_value = params_APS['field_1_params']['slope_s_value'], 
                                                  path_to_mb_file=params_APS['field_1_params']['path_to_slope_s_file'], 
                                                  mu_rsd=params_APS['field_1_params']['mu_rsd'], Np=500)        
    W_f1_func = lambda x: W_f1_i['W_interp'](x) #if (x>=params_APS['field_1_params']['z_min'])*(x<=params_APS['field_1_params']['z_max']) else 0    
    params_APS['field_1_params']['header'] = "[1] l, [>2] {f1}-{f1} Cl [multipoles, redshift/frequency bin] --survey: {s} | --dataset: {d} | bin: {b}".format(f1=params_APS['field_1'].upper(),
                                                                                                                                                               s=params_APS['field_1_params']['survey'].upper(), 
                                                                                                                                                               d=params_APS['field_1_params']['dataset'].upper(),
                                                                                                                                                               b=params_APS['field_1_params']['bin_to_use']
                                                                                                                                                               )    
    
    if params_general['verbose']: print('Field 1: COUNT')
    del params_count_survey
else:
    raise NameError('No Field named {}'.format( params_APS['field_1']))
####FIELD2    
if params_APS['field_2'].lower() == 'hi':
    params_APS['field_2_params'] = dcopy(params_hi_survey)
    params_APS['field_2_params'].update(hdata.nu_bins_vector(numin_ =params_APS['field_2_params']['freq_min'], 
                                        numax_ =params_APS['field_2_params']['freq_max'], 
                                        nbands_=params_APS['field_2_params']['nchannels'])
                                       )
    params_APS['field_2_params']['z'] = dcopy(np.flip(params_APS['field_2_params']['z']))    
    #W_f2_func  = lambda x: cxft.kernel_hi(camb_params=pars, camb_results=results, 
    #                                 z=x, zrange=params_APS['field_2_params']['z'], 
    #                                 config_dict=params_APS['field_2_params'] , replace_inf2nan=True)    
    W_f2_func = c4ft.kernel_hi_vec(
                #zout=zv, 
                camb_params=pars, 
                camb_results=results,
                zrange=params_APS['field_2_params']['z'], 
                replace_inf2nan=True, 
                use_input=False,
                config_dict=params_APS["field_2_params"], 
                params_interp=None, 
                get_interp_function=1
            ) 
    
    params_APS['field_2_params']['header'] = "[1] l, [>2] {}-{} Cl [multipoles, redshift/frequency bins]".format(params_APS['field_2'].upper(), params_APS['field_2'].upper())    
    del params_APS['field_2_params']['min'], params_APS['field_2_params']['max']    
    if params_general['verbose']: print('Field 2: HI')    
    del params_hi_survey
    #
elif params_APS['field_2'].lower() in ['cmbwl', 'wlcmb', 'convergence', 'cmb', 'kcmb', 'cmbk']:    
    params_APS['field_2_params'] = dcopy(params_cmbwl_survey)    
    params_APS['field_2_params']['z_min'] = 0
    params_APS['field_2_params']['z_max'] = cxft.chi_vec(results)['z_star']    
    #params_APS['field_2_params']['quad_bound']=[0.1,100]
    #params_APS['field_2_params']['np_ini']=500000
    W_f2_func  = cxft.kernel_cmb_interp(camb_params=pars, camb_results=results, interp=True) #lambda function
    params_APS['field_2_params']['header'] = "[1] l, [2] {f2}-{f2} Cl [multipoles, redshift/frequency bin] --survey: CMB".format(f2=params_APS['field_2'].upper())    
    if params_general['verbose']: print('Field 2: CMB WL')        
    del params_cmbwl_survey
    #
elif params_APS['field_2'].lower() in ['lensing','lens', 'shear']:       
    params_APS['field_2_params'] = dcopy(params_shear_survey)  
    params_APS['field_2_params'].update(cxft.load_survey_parameters(survey_name = params_APS['field_2_params']['survey'], 
                                        data_release= params_APS['field_2_params']['dataset'],
                                        type_field='lensing',
                                        path_file=params_APS['field_2_params']['path_to_params_file'])
                                       )   
    ## NEED TO CHECK
    if type(None)!=type(nbins_shear):
        params_APS['field_2_params']['nbins']=dcopy(nbins_shear)    
    if type(None)!=type(sigmaz_shear):
        if 'LSST' in params_APS['field_2_params']['survey'].upper():
            params_APS['field_2_params']['sigma_z']=dcopy(sigmaz_shear)
        elif 'EUCLID' in params_APS['field_2_params']['survey'].upper():
            params_APS['field_2_params']['sigma_b']=dcopy(sigmaz_shear)
            params_APS['field_2_params']['sigma_0']=dcopy(sigmaz_shear)
        else:
            pass        
        
    W_f2_j = cxft.kernel_galaxy_lensing_i_dict(ibin=params_APS['field_2_params']['bin_to_use'],
                                              camb_params=pars, camb_results=results, binning=params_shear_survey['binning'], 
                                              pars=params_APS['field_2_params'], verbose=0)
    W_f2_func = lambda x: W_f2_j['W_interp'](x) #if (x>=params_APS['field_2_params']['z_min'])*(x<=params_APS['field_2_params']['z_max']) else 0    
    params_APS['field_2_params']['header'] = "[1] l, [>2] {f1}-{f1} Cl [multipoles, redshift/frequency bin] --survey: {s} | --dataset: {d} | bin: {b}".format(f1=params_APS['field_2'].upper(),
                                                                                                                                                               s=params_APS['field_2_params']['survey'].upper(), 
                                                                                                                                                               d=params_APS['field_2_params']['dataset'].upper(),
                                                                                                                                                               b=params_APS['field_2_params']['bin_to_use']
                                                                                                                                                               )    
    #params_APS['field_2_params']['quad_bound']=None#[0.1,100]
    #params_APS['field_2_params']['np_ini']=500    
    if params_general['verbose']: print('Field 2: SHEAR')        
    del params_shear_survey

elif params_APS['field_2'].lower() in ['galaxy', 'density', 'dens', 'source', 'overdensity', 'count', 'counting']:           
    params_APS['field_2_params'] = dcopy(params_count_survey)  
    params_APS['field_2_params'].update(cxft.load_survey_parameters(survey_name  = params_APS['field_2_params']['survey'], 
                                                                    data_release = params_APS['field_2_params']['dataset'], 
                                                                    type_field   = 'density',
                                                                    path_file    = params_APS['field_2_params']['path_to_params_file'])
                                       )       
    if type(None)!=type(nbins_count):
        params_APS['field_2_params']['nbins']=dcopy(nbins_count)
    if type(None)!=type(sigmaz_count):
        if 'LSST' in params_APS['field_2_params']['survey'].upper():
            params_APS['field_2_params']['sigma_z']=dcopy(sigmaz_count)
        elif 'EUCLID' in params_APS['field_2_params']['survey'].upper():
            params_APS['field_2_params']['sigma_b']=dcopy(sigmaz_count)
            params_APS['field_2_params']['sigma_0']=dcopy(sigmaz_count)
        else:
            pass 

    W_f2_func = c4ft.kernel_galaxy_clustering_i_vec(camb_params         = pars, 
                                                    camb_results        = results, 
                                                    growth_rate_f       = None, 
                                                    #zout                = zv,
                                                    Np                  = 30000, 
                                                    pars                = params_APS['field_2_params'],
                                                    ibin                = params_APS['field_2_params']['bin_to_use'   ],
                                                    #include_mb          = params_APS['field_2_params']['include_mb'   ], 
                                                    #bias_model          = params_APS['field_2_params']['bias_model'   ], 
                                                    #bias_value          = params_APS['field_2_params']['bias_value'   ],
                                                    #mb_model            = params_APS['field_2_params']['slope_s_model'], 
                                                    #mb_alpha_value      = params_APS['field_2_params']['slope_s_value'], 
                                                    #path_to_mb_file     = params_APS['field_2_params']['path_to_slope_s_file'], 
                                                    #include_rsd         = params_APS['field_2_params']['include_rsd'  ], 
                                                    #mu_rsd              = params_APS['field_2_params']['mu_rsd'       ], 
                                                    #path_to_bias_file   = params_APS['field_2_params']['path_to_bias_file'],
                                                    get_interp_function = 1,
                                                    get_dict_format     = False,      
                                                    verbose             = False)
         
    params_APS['field_2_params']['header'] = "[1] l, [>2] {f1}-{f1} Cl [multipoles, redshift/frequency bin] --survey: {s} | --dataset: {d} | bin: {b}".format(f1=params_APS['field_2'].upper(),
                                                                                                                                                               s=params_APS['field_2_params']['survey'].upper(), 
                                                                                                                                                               d=params_APS['field_2_params']['dataset'].upper(),
                                                                                                                                                               b=params_APS['field_2_params']['bin_to_use']
                                                                                                                                                               ) 
    #params_APS['field_2_params']['quad_bound']=None#[0.1,100]
    #params_APS['field_2_params']['np_ini']=500        
    if params_general['verbose']: print('Field 2: COUNT')
    del params_count_survey
else:
    raise NameError('No Field named {}'.format( params_APS['field_2']))   
try:    
    params_APS['prefix']=params_APS['prefix'].format(b=params_APS['field_2_params']['bin_to_use'   ])    
except:
    pass
#####################################################################################################
#####################################################################################################
if params_general['verbose']:
    print('\n------------------------------------------')
    print('Starting to produce Angular Power Spectrum')
    print('------------------------------------------\n')

if params['save_Pk']:
    if params_general['verbose']: print('Generating Matter Power Spectrum (Pk)...')
    ksave = np.logspace(np.log10(1e-5), np.log10(params['kmax']), 200)
    zsave = params['save_Pk_at_z']
    zstring = np.array2string(zsave, precision=4, separator=', ')
    if params_general['verbose']: print('Matter Pk at z: {}\n'.format(zstring))
    #np.array2string(zsave, precision=4, separator=', ')
    Pksave = []
    for i,iz in enumerate(zsave):
        Pksave.append(pk_interp.P(iz,ksave))
    Pksave = np.asarray( Pksave )

    _ = cxft.savedata_pk(k_ = ksave, 
                         Pk_= np.asmatrix(Pksave), 
                   filename = 'matterPk', 
                     prefix = params_APS['prefix'], 
                     suffix = params_APS['suffix'], 
                     path   = params_APS['pathout'],                                  
                     header =  "[1] k [1/Mpc], [>2] Pk (k,z) [Mpc3], z:{} ".format(zstring)
                        )        
    del zstring, zsave, ksave, Pksave

if params_APS['generate_cl_field_1']:
    if params_general['verbose']: print('Generating {} Ang Power Spectrum...'.format(params_APS['field_1'].upper()))
    # Mandatory f1=HI -- the same estructure can be use to gamma
    Cl_f1 = np.array([ cxft.angular_power_spectrum_ij(camb_params =pars, 
                                                      camb_results=results, 
                                                      Pk_camb_interp=pk_interp,
                                                      params_aps=params_APS,
                                                      W_f1_i=W_f1_func, 
                                                      W_f2_j=W_f1_func, 
                                                      l=params_APS['l'], 
                                                      z_min=params_APS['field_1_params']['z'][i], 
                                                      z_max=params_APS['field_1_params']['z'][i+1],
                                                      kind_integrator='quad',
                                                      quad_limit=100,
                                                      #quad_bound=None,
                                                      verbose=False)['cl'] 
                       for i in range(params_APS['field_1_params']['z'].size-1) ])   
    
    ind = np.where((Cl_f1<params_general['threshold'])*(Cl_f1>-params_general['threshold']))
    Cl_f1[ind]=0    
    ###
    if params_general['verbose']: print('Generated.\nSaving...')
    import creating_dir
    creating_dir.verification_dir(params_APS['pathout'].split('/')[-1], "/".join(( params_APS['pathout'].split('/')[:-1] )))
    #creating_dir.verification_dir("auxiliary_data",params_APS['pathout'])
    pathname = cxft.savedata(l_ = params_APS['l'], 
                             Cl_= np.asmatrix(Cl_f1), filename='{}_cl'.format(params_APS['field_1'].upper()), 
                             prefix = params_APS['prefix'] , suffix=params_APS['suffix'], 
                             path   = params_APS['pathout'], header=params_APS['field_1_params']['header'])        
    if params_general['verbose']: print('Saved at {}'.format(pathname))
    print('\n')
else:
    if params_general['verbose']: print('skipping {} Ang Power Spectrum'.format(params_APS['field_1'].upper()))
    print('\n')
    
if params_APS['generate_cl_field_2']:
    if params_general['verbose']: print('Generating {} Ang Power Spectrum...'.format(params_APS['field_2'].upper()))

    if params_APS['field_2'].lower() in ['galaxy', 'density', 'dens', 'source', 'overdensity', 'count', 'counting']:
        zrange   = np.linspace(params_APS["field_2_params"]['z_min'], 
                               params_APS["field_2_params"]['z_max'], 
                                                                     1001)
        dNdz_vec = np.array([ cxft.nz_spec(pars=params_APS["field_2_params"])(ix) for ix in zrange ])
        zgbins   = cxft.compute_equal_number_bounds_lsst(redshift_range=zrange, 
                                                  redshift_distribution=dNdz_vec, 
                                                                 n_bins=params_APS["field_2_params"]['nbins'])
        gbin   = params_APS["field_2_params"]["bin_to_use"]
        zgbins = zgbins[gbin:gbin+2]
        quad_limit=5000
        np_ini = 500
        #print(20*'01
        del dNdz_vec,zrange,gbin                
    elif params_APS['field_2'].lower() in ['cmbwl', 'wlcmb', 'convergence', 'cmb', 'kcmb', 'cmbk']:   
        zgbins = [0.1,100]
        quad_limit=5000
        np_ini = 500000
    else:
        zgbins=None
        quad_limit=100
        np_ini=500
    Cl_f2 = cxft.angular_power_spectrum_ij(camb_params =pars, 
                                           camb_results=results, 
                                           Pk_camb_interp=pk_interp,
                                           params_aps=params_APS,
                                           W_f1_i=W_f2_func, 
                                           W_f2_j=W_f2_func, 
                                           l=params_APS['l'], 
                                           z_min=params_APS['field_2_params']['z_min'], 
                                           z_max=params_APS['field_2_params']['z_max'],
                                           kind_integrator='quad',
                                           quad_limit=quad_limit,
                                           quad_bound=zgbins,
                                           #np_ini=np_ini,
                                           verbose=False)['cl'] 
    ###
    ind = np.where((Cl_f2<params_general['threshold'])*(Cl_f2>-params_general['threshold']))
    Cl_f2[ind]=0    
    ###
    if params_general['verbose']: print('Generated.\nSaving...')
    import creating_dir
    creating_dir.verification_dir(params_APS['pathout'].split('/')[-1], "/".join(( params_APS['pathout'].split('/')[:-1] )))
    pathname = cxft.savedata(l_ = params_APS['l'], 
                             Cl_= np.asmatrix(Cl_f2), filename='{}_cl'.format(params_APS['field_2'].upper()), 
                             prefix = params_APS['prefix'] , suffix=params_APS['suffix'], 
                             path   = params_APS['pathout'], header=params_APS['field_2_params']['header'])  
    if params_general['verbose']: print('Saved at {}'.format(pathname))
    print('\n') 
else:
    if params_general['verbose']: print('skipping {} Ang Power Spectrum'.format(params_APS['field_2'].upper()))
if params_APS['generate_cl_cross']:
    if params_general['verbose']: print('Generating {}x{} Ang Power Spectrum...'.format(params_APS['field_1'].upper(), params_APS['field_2'].upper()))
    for i in range(params_APS['field_1_params']['z'].size-1):
        zxmin = np.amax([params_APS['field_1_params']['z'][i:i+2].min(),params_APS['field_2_params']['z_min']])
        zxmax = np.amin([params_APS['field_1_params']['z'][i:i+2].max(),params_APS['field_2_params']['z_max']])    
        if zxmin<=zxmax:
            if params_APS['field_2'].lower() in ['galaxy', 'density', 'dens', 'source', 'overdensity', 'count', 'counting']:
                zrange   = np.linspace(params_APS["field_2_params"]['z_min'], params_APS["field_2_params"]['z_max'], 1001)
                dNdz_vec = np.array([ cxft.nz_spec(pars=params_APS["field_2_params"])(ix) for ix in zrange ])
                zgbins   = cxft.compute_equal_number_bounds_lsst(redshift_range=zrange, 
                                                           redshift_distribution=dNdz_vec, 
                                                           n_bins=params_APS["field_2_params"]['nbins'])
                gbin     = params_APS["field_2_params"]["bin_to_use"]
                zgbins=zgbins[gbin:gbin+2]
                quad_limit=5000
                np_ini = 500
                del dNdz_vec,zrange,gbin 
            elif params_APS['field_2'].lower() in ['cmbwl', 'wlcmb', 'convergence', 'cmb', 'kcmb', 'cmbk']:   
                zgbins = None#[0.1,100]
                quad_limit=100
                np_ini = 500000
            else:
                zgbins=None
                quad_limit=5000
                np_ini = 500
            
            cx = cxft.angular_power_spectrum_ij(camb_params =pars, 
                                                camb_results=results, 
                                                Pk_camb_interp=pk_interp,
                                                params_aps=params_APS,
                                                W_f1_i=W_f1_func, 
                                                W_f2_j=W_f2_func, 
                                                l=params_APS['l'], 
                                                z_min=zxmin, 
                                                z_max=zxmax,
                                                kind_integrator='quad',
                                                quad_limit=quad_limit,
                                                quad_bound=zgbins,
                                                #np_ini=np_ini,
                                                verbose=False)['cl'] 
        else:
            cx = 0*params_APS['l']
        Cl_f1f2=dcopy(cx) if not i else np.vstack(( Cl_f1f2,dcopy(cx) ))
        del cx        
    ind = np.where((Cl_f1f2<params_general['threshold'])*(Cl_f1f2>-params_general['threshold']))
    Cl_f1f2[ind]=0 
    ###
    if params_general['verbose']: print('Generated.\nSaving...')
    import creating_dir
    creating_dir.verification_dir(params_APS['pathout'].split('/')[-1], "/".join(( params_APS['pathout'].split('/')[:-1] )))
    #creating_dir.verification_dir("auxiliary_data",params_APS['pathout'])
    pathname = cxft.savedata(l_ = params_APS['l'], 
                             Cl_= np.asmatrix(Cl_f1f2), filename='{}x{}_cl'.format(params_APS['field_1'].upper(), params_APS['field_2'].upper()), 
                             prefix = params_APS['prefix'] , suffix=params_APS['suffix'], 
                             path   = params_APS['pathout'], 
                             header =  "[1] l, [>2] {}-{} Cl [multipoles, redshift/frequency bins]".format(params_APS['field_1'].upper(), params_APS['field_2'].upper())
                            )      
    if params_general['verbose']: print('Saved at {}'.format(pathname))
    print('\n')
else:
    if params_general['verbose']:  print('skipping {}x{} Ang Power Spectrum...'.format(params_APS['field_1'].upper(), params_APS['field_2'].upper()))
if params_general['verbose']:
    print('------------------------------------------')
    print('END: {0:.2f} seg'.format(time.time()-timei))
    print('------------------------------------------\n')