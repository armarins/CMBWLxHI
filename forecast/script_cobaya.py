#############################################
#MPI
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
#############################################
#COBAYA
from cobaya.log import LoggedError
from cobaya.run import run
success = False
#############################################
#BASIC SYSTEM PACKAGES
import sys,time, os, camb, h5py, yaml
from copy import deepcopy as dcopy
#############################################
#FURTHER BASIC PACKAGES
from scipy import integrate
from scipy.interpolate import interp1d
import numpy as np
import healpy as hp
import pymaster as nmt
#############################################
#PATH TO SCRIPTS FOLDER
path = '/data2/LSSxHI-CODES/scripts/'
sys.path.insert(1, path)
import cross_functions_theory      as cxft
import cross_functions_theory_4forecast as c4ft
import handling_data               as hdata

class CXLikelihood:
    """
    Attributes:
        l_data     (array): observed data multipole
        cl_data    (array): observed data spectrum
        errcl_data (array): observed data spectrum error

    Methods:
        logp(ombh2, tf): Description of method2 with parameters.
    """    
    def __init__(self, kh_data, cl_data, errcl_data):#, params_cosmo):
        self.l_data     = l_data
        self.cl_data    = cl_data
        self.errcl_data = errcl_data
        self.inv_cov = np.diag(1.0/self.errcl_data**2) # Assuming diagonal covariance

    def logp(self, H0, ombh2, tf):
        """
        Calculates the log-likelihood.
        """
        current_params = {"H0": H0, "ombh2": ombh2, "tf": tf}
        params_cosmo_fit = dcopy(params_cosmo)
        params_cosmo_fit['ombh2'] = dcopy(  current_params['ombh2'] )
        params_cosmo_fit['tf'   ] = dcopy(  current_params['tf']    )
        params_cosmo_fit['H0'   ] = dcopy(  current_params['H0']    )
        try:
            pars_camb   = get_camb_results(params_cosmo_fit)
            l_th, cl_th = get_cx_theory(pars_camb)
            cl_eff = params_APS['namaster']['lbinner'](cl_th)
            cl_eff = cl_eff[WSEL]
            if np.isnan(cl_eff).any():
                return -np.inf
            diff = self.cl_data*params_cosmo_fit['tf'] - cl_eff
            chi2 = diff @ self.inv_cov @ diff
        except Exception as e:
            print(f"Error in CAMB or interpolation for params {pars_camb}: {e}")
            return -np.inf # Return negative infinity log-likelihood for invalid points
        del pars_camb
        return -0.5 * chi2
    
def get_camb_results(pars_camb):
    """
    Calculate PK from CAMB
    
    Spectrum and modes in: 
        [k]  = 1/Mpc
        [Pk] = 1/Mpc^3
        
    return :dict: with k,pk,z, and camb results instance
    """
    pars_camb['As']=2e-9
    cpars = camb.set_params(halofit_version=pars_camb['halofit_version'])
    cpars.set_cosmology(H0=pars_camb['H0'], 
                        ombh2=pars_camb['ombh2'], 
                        omch2=pars_camb['omch2'],
                        mnu=0.06)
    cpars.InitPower.set_params(As=pars_camb['As'], 
                               ns=pars_camb['ns'])
    cpars.set_matter_power(redshifts=pars_camb['z'], 
                           kmax=pars_camb['kmax'], silent=1)    
    cresults = camb.get_results(cpars)
    pk_interp = camb.get_matter_power_interpolator(cpars, 
                                                   nonlinear=pars_camb['Pk_nonlinear'], hubble_units=False, k_hunit=False, 
                                                   kmax=pars_camb['kmax'], zmax=pars_camb['z'][-1])
    #kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=10.0, npoints=200)
    kc,zc,pkc = cresults.get_nonlinear_matter_power_spectrum(have_power_spectra=1, hubble_units=False, k_hunit=False, params=cpars)
    zstar=cresults.get_derived_params()['zstar']
    return {'k':kc,'z':zc,'pkc':pkc, 
            'pk': pk_interp.P,
            'comovel_dist_z':cresults.comoving_radial_distance,
            'hubble_z':      cresults.hubble_parameter,
            'fsigma8':       cresults.get_fsigma8(), 
            'sigma8':        cresults.get_sigma8(), 
            'zstar':         zstar,
            'chistar':       cresults.comoving_radial_distance(zstar)
           }
    #return {'k':kc,'z':zc,'pk':pkc, 'results':cresults, 'params':cpars}    
def get_cx_theory(params=None,):
    l_           = params_APS['l']
    hubble_z     = params["hubble_z"]
    Pk_interp    = params["pk"]
    comovel_dist = params["comovel_dist_z"]
    W_f1_vec = c4ft.kernel_hi_vec(zout=zv, zrange=zv,replace_inf2nan=1, use_input=1, config_dict=params_APS['field_1_params'] ,params_interp=params )
    W_f2_vec = c4ft.kernel_galaxy_clustering_i_vec(zout= zv, pars = params_APS['field_2_params'], ibin=params_APS['field_2_params']['bin_to_use'], verbose = False  )  
    for ii, iil in enumerate(l_):
        dCdz =  np.hstack([ (hubble_z(x)/c_light)*Pk_interp(x,(iil+0.5)/comovel_dist(x))/(comovel_dist(x)**2) for x in zv ])[0]
        dCdz *= W_f1_vec*W_f2_vec
        C_ = integrate.simpson(y=dCdz, x=zv)  if not ii else np.hstack(( C_, integrate.simpson(y=dCdz, x=zv)  ))
        #Cz = np.sum(dCdz[1:]*np.diff(zv)) if not ii else np.hstack(( Cz, np.sum(dCdz[1:]*np.diff(zv)) ))
        #Cz = np.trapezoid(dCdz, zv)       if not ii else np.hstack(( Cz, np.trapezoid(dCdz, zv)       ))
    return l_,C_

#############################################
## DEFINITIONS
c_light = 299792458/1e3 #km/s
nu_HI   = 1420.405751768 #MHz
THRESHOLD = 1e-20
STERADIAN = np.radians(1/60)**2 #arcmin to steradian
PATH_LSST   = "/data2/LSSxHI-CODES/scripts/parameters/lsst_desc_parameters.yaml"
PATH_EUCLID = "/data2/LSSxHI-CODES/scripts/parameters/euclid_parameters.yaml"
filename_cl_h5 = 'input/cl_gcLSSTxLOWZ_4x19_40_0001_delell20.h5'
file_loadyaml  = 'input/generating_APS_theory.yaml'
filename_info  = 'input/generating_theory-MH.yaml'
#############################################
## BINS
NUM_Z = 100
HBIN  = 19
## NAMASTER SETTINGS 
DEL_L = 20
NSIDE = 256
NPIX  = 12*NSIDE*NSIDE
LMIN  = 30
LMAX  = 400
##
b    = nmt.NmtBin.from_nside_linear(NSIDE, nlb=DEL_L)
leff = b.get_effective_ells()
feff = leff*(leff+1)/2/np.pi
WSEL = (leff>LMIN)*(leff<LMAX)

params_loaded = c4ft.load_yaml_fileformat(input_filename=file_loadyaml)

params_general = dcopy(params_loaded['general'])
params_APS     = dcopy(params_loaded['APS'])
params_cosmo   = dcopy(params_loaded['cosmology'])
params_cosmo_fit = dcopy(params_cosmo)
del params_loaded

params_APS['del_l']=DEL_L
params_APS['nside']=NSIDE
params_APS['npix' ]=NPIX
params_APS['namaster']={'b':b,'leff':leff,'feff':feff, 'lbinner':b.bin_cell}

#LMIN,LMAX=30,400
with h5py.File(filename_cl_h5, "r") as f:
    l_data     = f["l"][:]
    cl_data    = f["cl"][:]
    errcl_data = f["errcl"][:]
####
zxmin = np.amax([params_APS['field_1_params']['z'][HBIN:HBIN+2].min(),params_APS['field_2_params']['z_min']])
zxmax = np.amin([params_APS['field_1_params']['z'][HBIN:HBIN+2].max(),params_APS['field_2_params']['z_max']])  
####
zv    = np.logspace( np.log10(zxmin), np.log10(zxmax), NUM_Z)

# ---------------------- Create Likelihood Instance ----------------------
likelihood = CXLikelihood(l_data, cl_data, errcl_data)
info_loaded = c4ft.load_yaml_fileformat(input_filename=filename_info, verbose=1)
info_loaded['likelihood']['matter_power_spectrum']['external']=likelihood.logp

print("Starting Cobaya run...")
try:
    upd_info, mcmc = run(info_loaded)
    success = True
except LoggedError as err:
    pass

# Did it work? (e.g. did not get stuck)
success = all(comm.allgather(success))

if not success and rank == 0:
    print("Sampling failed!")
print("Cobaya run finished.")



