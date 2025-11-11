"""
External Cobaya likelihood:
Cross-correlation between HI-intensity-mapping and galaxy clustering.
"""

# ---- standard imports -------------------------------------------------------
import os, sys, camb, h5py
from copy import deepcopy as dcopy
from scipy import integrate
import pymaster as nmt
import numpy as np
from typing import List, Optional
from scipy.interpolate import interp1d
from cobaya.likelihood import Likelihood

path = '/data/AMARINS/LSSxHI-CODES/scripts/'
if path not in sys.path: # Avoid adding the same path multiple times
    sys.path.insert(1, path)
import cross_functions_theory as cxft
import cross_functions_theory_4forecast as c4ft
import handling_data as hdata

# ---- constants --------------------------------------------------------------
C_LIGHT = 299792458 / 1e3  # km s⁻¹
STERADIAN = np.radians(1 / 60) ** 2
EPS=1e-50

##############################################################################
## NEW COSMOLOGICAL PARAMS -- from CAMB
# - (1) Introduce in the YAML file as either fixed or variable. 
#       Params with only 'value' argument are assumed fixed
# - (2) Set params in the get_requirements()
# - (3) If they are variable to sample, set them as logp()'s. 
#       The parameter must be used to update the p_fit and then 
#       the camb.cosmology
# - (4) make sure you are setting _get_camb_results() correctly
#
## NEW DATASET -- created
# - (1) set the name of h5 data in the YAML file, likelihoof. It is necessary write a new variable with its: filename, general APS information about the cross correlation, HI bin to cross correlated.
# - (2) Initialize the parameters in the CXLilkelihood class: f.ex.
#     data_2_h5: str
#     aps_2_yaml: str
#     hi_2_index: int = 19
# - (3) They will be automatically saved in the self

# PS.: Each valor must be spaced as only one null character after ":"
# PS.: Numpy variables are not accepted
# PS.: variable infinity must be: .inf, +.inf, or -.inf
# PS.: Strings without comma
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
class CXLikelihood(Likelihood):
    """Cobaya external likelihood implementing *logp*.
    Uses Cobaya's central CAMB instance and processes two datasets.
    """

    # Parameters from the likelihood's own YAML block
    integrator: str = 'faster'
    include_f1xf1_cl: bool = False
    omegahi_model: str
    nside: int = 256    
    num_z: int = 100
    lmin: int = 30
    lmax: int = 400
    del_l: int = 20      
    
    path_yaml: str
    
    bins_f1f1: Optional[List[float]] = None    
    f1_field: str = 'HI'
    data_f1_h5_1: str 

    # ---------------------------------------------------------------------
    def initialize(self):
        """
        Initialize likelihood
        
        Loading general information
        Set general variables
        
        """  
        p_loaded          = c4ft.load_yaml_fileformat(self.path_yaml, verbose=0)
        self.params_aps   = dcopy(p_loaded["APS"])
        self.params_cosmo = dcopy(p_loaded["cosmology"]) 
        del p_loaded        
            
        if self.include_f1xf1_cl:
            # ---- Load observables for auto-HIxHI ----
            #data_hi_h5_path_1  = os.path.expanduser(self.data_hi_h5_1)                 
            with h5py.File(self.data_f1_h5_1, "r") as f:
                self.lf1_data_1    = f["l"    ][:]
                self.cf1_data_1    = f["cl"   ][:]
                self.errcf1_data_1 = f["errcl"][:]
            self.inv_covcf1_1 = np.diag(1.0 / self.errcf1_data_1**2)
    
        # ---- NaMaster binning for dataset 1 ----
        self.b    = nmt.NmtBin.from_nside_linear(self.nside, nlb=self.del_l)
        self.leff = self.b.get_effective_ells()
        self.WSEL = (self.leff > self.lmin) & (self.leff < self.lmax)
        self.params_aps["namaster"] = {
            "b": self.b, "leff": self.leff, "lbinner": self.b.bin_cell,
        }
        self.num_f1f1_bins = len(self.bins_f1f1)
        self.params_aps["field_1_params"]['bins_to_use'] = {'f1f1': {'rank':np.arange(self.num_f1f1_bins),'bin':np.asarray(self.bins_f1f1)} }       
        
    # ------------------------------------------------------------------
    def get_requirements(self):
        """Exact cosmological outputs needed from CAMB."""

        z_full = self.params_aps["field_1_params"]["z"]
        z_bg   = np.unique(np.concatenate(([0.0], z_full )))   # adds 0 #<-- se add valores
        self.z_full = z_full
        self.z_bg   = z_bg

        k_max_PK = 100.0  # h/Mpc

        return {
            "Pk_interpolator": {"z": z_full.tolist(),  
                                "k_max": k_max_PK,
                                "nonlinear": True,
                                "hubble_units":False,
                                "k_hunit":False},
            "Hubble": {"z": z_bg.tolist()},
            "comoving_radial_distance": {"z": z_bg.tolist()},
            "sigma8": None,
            "w": None,
            "wa":None,
            "ombh2":None,
            "omch2":None,
            "H0":None,
            "ns":None,
            "As":None,
            "mnu":None,
            #"ax_1": None,# "ax_2": None, #"ax_3": None,
            #"ah_1": None,# "ah_2": None, #"ah_3": None,            
            "omegahi": None
        }
    # ---------------------------------------------------------------------
    def logp(self, 
             #ax_1,# ax_2, #ax_3, 
             #ah_1,# ah_2, #ah_3,                
             #omegahi,
            **params_values #it is a dict for the required CAMB quantities and sampled ones
            ):
        """Calculate the log-likelihood using products from Cobaya's CAMB."""

        ombh2_val = self.provider.get_param("ombh2")
        omch2_val = self.provider.get_param("omch2")
        H0_val = self.provider.get_param("H0")
        w_val  = self.provider.get_param("w")     
        omegahi = params_values['omegahi']
        om_val = (ombh2_val+omch2_val) / (H0_val/100)**2
        try:
            Pk_interp_func  = self.provider.get_Pk_interpolator() 
            hubble_z_func   = self.provider.get_Hubble 
            comov_dist_func = self.provider.get_comoving_radial_distance
            sigma8_val      = self.provider.get_param("sigma8")
        except Exception as e:
            self.log.error(f"[LOGP] Failed to get products from provider: {e}")
            return -np.inf    
        #########################################
        print_test_cosmo=False
        if print_test_cosmo:
            z_val = 0.1;k_val = 0.2  # [1/Mpc]
            power_at_point = Pk_interp_func.P(z=z_val, k=k_val)
            print()        
            print("="*60)
            print("Cosmological Parameters & Values for this step:")
            for name, value in params_values.items():
                print(f"  - {name:<10}: {value}")
            power_at_point = Pk_interp_func.P(z=z_val, k=k_val)
            print(f"P(k={k_val}, z={z_val}) = {power_at_point:.4f} (Mpc/h)^3")        
            #print(Pk_interp_func(z_val, k_val))
            print("="*60)
            print()
        #########################################                    
        hubble_z_func  = interp1d(self.z_bg, hubble_z_func(self.z_bg) , kind='cubic')
        self.provider.get_Hubble = dcopy( hubble_z_func )
        comov_dist_func = interp1d(self.z_bg, comov_dist_func(self.z_bg) , kind='cubic')
        self.provider.get_comoving_radial_distance= dcopy( comov_dist_func )        

        # Prepare a dictionary of cosmological functions/parameters for your kernels
        # This needs to be compatible with what `kernel_hi_vec` expects for `params_interp`.
        # If the kernels explicitly need H0, w, or ombh2, pass them.
        # CAMB products (Pk, H(z), chi(z)) are already functions of the full cosmology.
        cosmo_functions_for_kernels = {
            "hubble_z"      : hubble_z_func,
            "comovel_dist_z": comov_dist_func,
            "sigma8"        : sigma8_val,
             "H0"           : H0_val,
             "w"            : w_val,
             "ombh2"        : ombh2_val
                                     }
        try:
            total_chi2=0
            if self.include_f1xf1_cl:
                # --- Calculate chi2 for dataset 1 ---                
                self.params_aps["field_1_params"]['omegaHI'] = omegahi   
                if self.f1_field.lower() in ['ch', 'chi', 'hixhi','hi']:
                    hi_index = self.params_aps['field_1_params']['bins_to_use']['f1f1']['bin'][0]#<--loop over bins
                    self.params_aps['field_1_params']['bin_to_use'] = hi_index
                    zf1  = self.params_aps["field_1_params"]["z"]
                    zmin = zf1[hi_index : hi_index + 2].min()
                    zmax = zf1[hi_index : hi_index + 2].max()
                    self.zv = np.logspace(np.log10(zmin), np.log10(zmax), self.num_z)                    
                    
                    _, chi_theory_1 = self._get_cl_hi_theory(
                                        Pk_interp_func=Pk_interp_func.P, 
                                        hubble_z_func=hubble_z_func, 
                                        comov_dist_func=comov_dist_func, 
                                        cosmo_for_kernels=cosmo_functions_for_kernels, 
                                        params_aps_i=self.params_aps, 
                                        zv_i=self.zv        )
                    chi_eff_1 = self.params_aps["namaster"]["lbinner"](chi_theory_1)[self.WSEL]
                else:
                    print('ERROR')
                    sys.exit(0)
                diff_1 = self.cf1_data_1 - chi_eff_1 
                chi2_1 = diff_1 @ self.inv_covcf1_1 @ diff_1          
                total_chi2 += chi2_1
            return -0.5 * total_chi2

        except Exception as exc:
            print(f"-"*60)
            for name, value in params_values.items():
                print(f"  - {name:<10}: {value}")            
            self.log.warning(
                f"[LOGP] Calculation failed at:"
                f"Error: {exc}"
                            )
            print(f"-"*60)
            import traceback 
            self.log.error(traceback.format_exc()) 
            return -np.inf        

    # ---------------------------------------------------------------------
    def _get_cl_hi_theory(self, 
                          Pk_interp_func, 
                          hubble_z_func, 
                          comov_dist_func, 
                          cosmo_for_kernels, 
                          params_aps_i, 
                          zv_i):#ask for g_index_i to be used in the kernel
        """
        Compute theoretical HixHi-Cl values using products from Cobaya's CAMB
        for a single dataset configuration.
        """
        l_grid = params_aps_i["l"]

        zv_i = np.linspace( zv_i.min(), zv_i.max(), 500 )
        dz   = np.diff(zv_i)[0]
        zv_i = np.linspace( zv_i.min()+2*dz, zv_i.max()-2*dz, 500 )
        W_f1 = c4ft.kernel_hi_vec(
                camb_params=None,
                camb_results=None,
                zout=zv_i,
                zrange=params_aps_i['field_1_params']['z'], 
                replace_inf2nan=True, 
                use_input=True,
                config_dict=params_aps_i["field_1_params"],  #<--- edit Omega HI
                params_interp=cosmo_for_kernels, 
                get_interp_function=0
            )
        
        Cl = []
        for ell in l_grid:
            chi_at_zv = comov_dist_func(zv_i)
            k_values = (ell+0.5)/chi_at_zv 
            pk_values_on_grid = Pk_interp_func(zv_i, k_values, grid=False)
            dCdz = ( (hubble_z_func(zv_i)/C_LIGHT)*pk_values_on_grid/(chi_at_zv**2) )
            dCdz *= W_f1*W_f1
            Cl.append(integrate.simpson(dCdz, zv_i))
        return l_grid, np.asarray(Cl)        