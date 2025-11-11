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
    include_gxhi_cl: bool = True
    include_hixhi_cl: bool = False
    include_gxg_cl: bool = False
    omegahi_model: str
    
    data_cx_h5_1: str
    aps_cx_yaml_1: str
    data_hi_h5_1: str
    g_index_1: int = 4
    hi_index_1: int = 19
    num_z_1: int = 100
    lmin_1: int = 30
    lmax_1: int = 400
    del_l_1: int = 20
    
    #data_cx_h5_2: str
    #aps_cx_yaml_2: str
    #data_hi_h5_2: str    
    #g_index_2: int = 4    
    #hi_index_2: int = 20 
    #num_z_2: int = 100
    #lmin_2: int = 30
    #lmax_2: int = 400
    #del_l_2: int = 20    

    #data_cx_h5_3: str
    #aps_cx_yaml_3: str
    #data_hi_h5_3: str    
    #g_index_3: int = 4    
    #hi_index_3: int = 21 
    #num_z_3: int = 100
    #lmin_3: int = 30
    #lmax_3: int = 400
    #del_l_3: int = 20        
    
    nside: int = 256

    # ---------------------------------------------------------------------
    def initialize(self):
        """
        Initialize likelihood: load data for two datasets, 
        and fixed settings from their respective YAML files.
        """  
        if self.include_gxhi_cl:
            # ---- Paths to input files ----
            data_cx_h5_path_1  = os.path.expanduser(self.data_cx_h5_1)
            aps_cx_yaml_path_1 = os.path.expanduser(self.aps_cx_yaml_1)
            #
            #data_cx_h5_path_2  = os.path.expanduser(self.data_cx_h5_2)
            #aps_cx_yaml_path_2 = os.path.expanduser(self.aps_cx_yaml_2)               
            #
            #data_cx_h5_path_3  = os.path.expanduser(self.data_cx_h5_3)
            #aps_cx_yaml_path_3 = os.path.expanduser(self.aps_cx_yaml_3)                           
            # ---- Load observables for dataset 1 ----
            with h5py.File(data_cx_h5_path_1, "r") as f:
                self.lx_data_1    = f["l"][:]
                self.cx_data_1    = f["cl"][:]
                self.errcx_data_1 = f["errcl"][:]
            self.inv_covcx_1 = np.diag(1.0 / self.errcx_data_1**2)
            # ---- Load observables for dataset 2 ----
            #with h5py.File(data_cx_h5_path_2, "r") as f:
            #    self.lx_data_2    = f["l"][:]
            #    self.cx_data_2    = f["cl"][:]
            #    self.errcx_data_2 = f["errcl"][:]
            #self.inv_covcx_2 = np.diag(1.0 / self.errcx_data_2**2)
            # ---- Load observables for dataset 3 ----
            #with h5py.File(data_cx_h5_path_3, "r") as f:
            #    self.lx_data_3    = f["l"][:]
            #    self.cx_data_3    = f["cl"][:]
            #    self.errcx_data_3 = f["errcl"][:]
            #self.inv_covcx_3 = np.diag(1.0 / self.errcx_data_3**2)
            
            # ---- Load auxiliary configuration for dataset 1 ----
            p_1_loaded         = c4ft.load_yaml_fileformat(aps_cx_yaml_path_1, verbose=0)
            self.params_APS_1  = dcopy(p_1_loaded["APS"])
            self.params_APS_1["field_1_params"]['omegaHI_model'] = self.omegahi_model       
            self.params_APS_1["field_1_params"]['fit_omegaHI']   = True
            self.params_APS_1["field_1_params"]["bin_to_use"]=self.hi_index_1           
            self.params_APS_1["field_2_params"]["bin_to_use"]=self.g_index_1        
            #self.params_APS_1["field_1_params"]['omegaHI']       =           
            # self.params_cosmo_1 = dcopy(p_1_loaded["cosmology"]) # CAMB params now global
            del p_1_loaded
    
            # ---- Load auxiliary configuration for dataset 2 ----
            #p_2_loaded         = c4ft.load_yaml_fileformat(aps_cx_yaml_path_2, verbose=0)
            #self.params_APS_2  = dcopy(p_2_loaded["APS"])
            #self.params_APS_2["field_1_params"]['omegaHI_model'] = self.omegahi_model       
            #self.params_APS_2["field_1_params"]['fit_omegaHI']   = True 
            #self.params_APS_2["field_1_params"]["bin_to_use"]=self.hi_index_2                   
            #self.params_APS_2["field_2_params"]["bin_to_use"]=self.g_index_2
            #del p_2_loaded
            # ---- Load auxiliary configuration for dataset 3 ----
            
            #p_3_loaded         = c4ft.load_yaml_fileformat(aps_cx_yaml_path_3, verbose=0)
            #self.params_APS_3  = dcopy(p_3_loaded["APS"])
            #self.params_APS_3["field_1_params"]['omegaHI_model'] = self.omegahi_model       
            #self.params_APS_3["field_1_params"]['fit_omegaHI']   = True 
            #self.params_APS_3["field_1_params"]["bin_to_use"]=self.hi_index_3
            #self.params_APS_3["field_2_params"]["bin_to_use"]=self.g_index_3
            #del p_3_loaded
            
        if self.include_hixhi_cl:
            data_hi_h5_path_1  = os.path.expanduser(self.data_hi_h5_1)
            data_hi_h5_path_2  = os.path.expanduser(self.data_hi_h5_2)
            #data_hi_h5_path_3  = os.path.expanduser(self.data_hi_h5_3)            
            # ---- Load observables for dataset 1 ----
            with h5py.File(data_hi_h5_path_1, "r") as f:
                self.lhi_data_1    = f["l"][:]
                self.chi_data_1    = f["cl"][:]
                self.errchi_data_1 = f["errcl"][:]
            self.inv_covchi_1 = np.diag(1.0 / self.errchi_data_1**2)
            # ---- Load observables for dataset 2 ----
            #with h5py.File(data_hi_h5_path_2, "r") as f:
            #    self.lhi_data_2    = f["l"][:]
            #    self.chi_data_2    = f["cl"][:]
            #    self.errchi_data_2 = f["errcl"][:]
            #self.inv_covchi_2 = np.diag(1.0 / self.errchi_data_2**2)
            # ---- Load observables for dataset 3 ----
            #with h5py.File(data_hi_h5_path_3, "r") as f:
            #    self.lhi_data_3    = f["l"][:]
            #    self.chi_data_3    = f["cl"][:]
            #    self.errchi_data_3 = f["errcl"][:]
            #self.inv_covchi_3 = np.diag(1.0 / self.errchi_data_3**2)       
    
        # ---- NaMaster binning for dataset 1 ----
        self.b_1    = nmt.NmtBin.from_nside_linear(self.nside, nlb=self.del_l_1)
        self.leff_1 = self.b_1.get_effective_ells()
        self.WSEL_1 = (self.leff_1 > self.lmin_1) & (self.leff_1 < self.lmax_1)
        self.params_APS_1["namaster"] = {
            "b": self.b_1, "leff": self.leff_1, "lbinner": self.b_1.bin_cell,
        }

        # ---- NaMaster binning for dataset 2 ----
        #self.b_2    = nmt.NmtBin.from_nside_linear(self.nside, nlb=self.del_l_2)
        #self.leff_2 = self.b_2.get_effective_ells()
        #self.WSEL_2 = (self.leff_2 > self.lmin_2) & (self.leff_2 < self.lmax_2)
        #self.params_APS_2["namaster"] = {
        #    "b": self.b_2, "leff": self.leff_2, "lbinner": self.b_2.bin_cell,
        #}

        # ---- NaMaster binning for dataset 2 ----
        #self.b_3    = nmt.NmtBin.from_nside_linear(self.nside, nlb=self.del_l_3)
        #self.leff_3 = self.b_3.get_effective_ells()
        #self.WSEL_3 = (self.leff_3 > self.lmin_3) & (self.leff_3 < self.lmax_3)
        #self.params_APS_3["namaster"] = {
        #    "b": self.b_3, "leff": self.leff_3, "lbinner": self.b_3.bin_cell,
        #}

        # ---- Redshift grid for dataset 1 integration ----
        #zf: HI IM redshift edges per bin
        #zmin: the minimum z for the hi_index bin
        #zmax: the maximum z for the hi_index bin
        #zv: logspaced z distribution with num_z points between zmin and zmax
        zf1_1  = self.params_APS_1["field_1_params"]["z"]
        zmin_1 = max(zf1_1[self.hi_index_1 : self.hi_index_1 + 2].min(),
                     self.params_APS_1["field_2_params"]["z_min"])
        zmax_1 = min(zf1_1[self.hi_index_1 : self.hi_index_1 + 2].max(),
                     self.params_APS_1["field_2_params"]["z_max"])
        if zmin_1 <= 0: raise ValueError("zmin_1 for logspace must be positive.")
        if zmax_1 <= zmin_1: raise ValueError("zmax_1 must be > zmin_1 for dataset 1.")
        self.zv_1 = np.logspace(np.log10(zmin_1), np.log10(zmax_1), self.num_z_1)
        #self.log.info(f"{self.name} initialized for dataset 1. Integration z: [{zmin_1:.4f}, {zmax_1:.4f}]")
        self.log.info(f"Dataset 1 initialized for dataset 1. Integration z: [{zmin_1:.4f}, {zmax_1:.4f}]")

        # ---- Redshift grid for dataset 2 integration ----
        #zf1_2  = self.params_APS_2["field_1_params"]["z"]
        #zmin_2 = max(zf1_2[self.hi_index_2 : self.hi_index_2 + 2].min(),
        #             self.params_APS_2["field_2_params"]["z_min"])
        #zmax_2 = min(zf1_2[self.hi_index_2 : self.hi_index_2 + 2].max(),
        #             self.params_APS_2["field_2_params"]["z_max"])
        #if zmin_2 <= 0: raise ValueError("zmin_2 for logspace must be positive.")
        #if zmax_2 <= zmin_2: raise ValueError("zmax_2 must be > zmin_2 for dataset 2.")
        #self.zv_2 = np.logspace(np.log10(zmin_2), np.log10(zmax_2), self.num_z_2)
        ##self.log.info(f"{self.name} initialized for dataset 2. Integration z: [{zmin_2:.4f}, {zmax_2:.4f}]")
        #self.log.info(f"Dataset 2 initialized for dataset 2. Integration z: [{zmin_2:.4f}, {zmax_2:.4f}]")

        # ---- Redshift grid for dataset 3 integration ----
        #zf1_3  = self.params_APS_3["field_1_params"]["z"]
        #zmin_3 = max(zf1_3[self.hi_index_3 : self.hi_index_3 + 2].min(),
        #             self.params_APS_3["field_2_params"]["z_min"])
        #zmax_3 = min(zf1_3[self.hi_index_3 : self.hi_index_3 + 2].max(),
        #             self.params_APS_3["field_2_params"]["z_max"])
        #if zmin_3 <= 0: raise ValueError("zmin_3 for logspace must be positive.")
        #if zmax_3 <= zmin_2: raise ValueError("zmax_3 must be > zmin_3 for dataset 3.")
        #self.zv_3 = np.logspace(np.log10(zmin_3), np.log10(zmax_3), self.num_z_3)
        #self.log.info(f"Dataset 3 initialized for dataset 3. Integration z: [{zmin_3:.4f}, {zmax_3:.4f}]")        
        #    
    # ------------------------------------------------------------------
    def get_requirements(self):
        """Exact cosmological outputs needed from CAMB."""

        # Combine both datasets’ redshifts *and* z=0 for background quantities
        z_full = self.params_APS_1["field_1_params"]["z"]
        #z_full = np.sort(np.unique( np.concatenate( [ 
        #    self.params_APS_1["field_1_params"]["z"], self.params_APS_2["field_1_params"]["z"] 
        #                                          ] )  
        #                          ) 
        #                )
        z_bg   = np.unique(np.concatenate(([0.0], z_full )))   # adds 0 #<-- se add valores
        self.z_full = z_full
        self.z_bg   = z_bg

        k_max_needed = 100.0  # h Mpc⁻¹

        return {
            "Pk_interpolator": {"z": z_full.tolist(),  
                                "k_max": k_max_needed,
                                "nonlinear": True},
            "Hubble": {"z": z_bg.tolist()},
            "comoving_radial_distance": {"z": z_bg.tolist()},
            "sigma8": None,
            "w": None,
            "ax_1": None,# "ax_2": None, #"ax_3": None,
            "ah_1": None,# "ah_2": None, #"ah_3": None,            
            "omegahi": None
        }
    # ---------------------------------------------------------------------
    def logp(self, 
             ax_1,# ax_2, #ax_3, 
             ah_1,# ah_2, #ah_3,              
             omegahi,
             **kwargs):
        """Calculate the log-likelihood using products from Cobaya's CAMB."""

        # Retrieve current cosmology from the provider
        # Fixed parameter (its value comes from the YAML)
        ombh2_val = self.provider.get_param("ombh2")
        omch2_val = self.provider.get_param("omch2")
        # Sampled parameters
        H0_val = self.provider.get_param("H0")
        w_val  = self.provider.get_param("w")       
        om_val = (ombh2_val+omch2_val) / (0.01*H0_val)**2
        try:
            Pk_interp_func  = self.provider.get_Pk_interpolator() 
            hubble_z_func   = self.provider.get_Hubble 
            comov_dist_func = self.provider.get_comoving_radial_distance
            sigma8_val      = self.provider.get_param("sigma8")
        except Exception as e:
            self.log.error(f"[LOGP] Failed to get products from provider: {e}")
            return -np.inf   
        hubble_z_func  = interp1d(self.z_bg, hubble_z_func(self.z_bg) , kind='cubic')
        self.provider.get_Hubble = dcopy( hubble_z_func )
        #Pk_interp_func = interp1d(self.z_bg, Pk_interp_func(self.z_bg) , kind='cubic')
        #self.provider.get_Hubble = dcopy( Pk_interp_func )
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
            if self.include_gxhi_cl:
                # --- Calculate chi2 for dataset 1 --- 
                self.params_APS_1["field_1_params"]['omegaHI'] = omegahi     
                _, cx_theory_1 = self._get_cx_theory(
                    Pk_interp_func, hubble_z_func, comov_dist_func, 
                    cosmo_functions_for_kernels, self.params_APS_1, self.zv_1
                )
                cx_eff_1 = self.params_APS_1["namaster"]["lbinner"](cx_theory_1)[self.WSEL_1]
                diff_1   = self.cx_data_1 * ax_1 - cx_eff_1 
                chi2_1   = diff_1 @ self.inv_covcx_1 @ diff_1
                
                # --- Calculate chi2 for dataset 2 ---
                #self.params_APS_2["field_1_params"]['omegaHI'] = omegahi     
                #_, cl_theory_2 = self._get_cx_theory(
                #    Pk_interp_func, hubble_z_func, comov_dist_func, 
                #    cosmo_functions_for_kernels, self.params_APS_2, self.zv_2
                #    )
                #cl_eff_2 = self.params_APS_2["namaster"]["lbinner"](cl_theory_2)[self.WSEL_2]
                #diff_2 = self.cx_data_2 * ax_2 - cl_eff_2
                #chi2_2 = diff_2 @ self.inv_covcx_2 @ diff_2
    
                # --- Calculate chi2 for dataset 3 ---
                #self.params_APS_3["field_1_params"]['omegaHI'] = omegahi     
                #    _, cl_theory_3 = self._get_cx_theory(
                #        Pk_interp_func, hubble_z_func, comov_dist_func, 
                #        cosmo_functions_for_kernels, self.params_APS_3, self.zv_3
                #    )
                #    cl_eff_3 = self.params_APS_3["namaster"]["lbinner"](cl_theory_3)[self.WSEL_3]
                #    diff_3 = self.cx_data_3 * ax_3 - cl_eff_3
                #    chi2_3 = diff_3 @ self.inv_covcx_3 @ diff_3
                
                total_chi2 += chi2_1 #+ chi2_2 #+ chi2_3      
            if self.include_hixhi_cl:
                # --- Calculate chi2 for dataset 1 ---                
                self.params_APS_1["field_1_params"]['omegaHI'] = omegahi     
                _, chi_theory_1 = self._get_cl_hi_theory(
                    Pk_interp_func, hubble_z_func, comov_dist_func, 
                    cosmo_functions_for_kernels, self.params_APS_1, self.zv_1
                )
                chi_eff_1 = self.params_APS_1["namaster"]["lbinner"](chi_theory_1)[self.WSEL_1]
                diff_1 = self.chi_data_1 * ah_1 - chi_eff_1 
                chi2_1 = diff_1 @ self.inv_covchi_1 @ diff_1          

                # --- Calculate chi2 for dataset 2 ---
                #self.params_APS_2["field_1_params"]['omegaHI'] = omegahi     
                #_, chi_theory_2 = self._get_cl_hi_theory(
                #    Pk_interp_func, hubble_z_func, comov_dist_func, 
                #    cosmo_functions_for_kernels, self.params_APS_2, self.zv_2
                #)
                #chi_eff_2 = self.params_APS_2["namaster"]["lbinner"](chi_theory_2)[self.WSEL_2]
                #diff_2 = self.chi_data_2 * ah_2 - chi_eff_2
                #chi2_2 = diff_2 @ self.inv_covchi_2 @ diff_2 
                
                # --- Calculate chi2 for dataset 3 ---
                #self.params_APS_3["field_1_params"]['omegaHI'] = omegahi     
                #_, chi_theory_3 = self._get_cl_hi_theory(
                #    Pk_interp_func, hubble_z_func, comov_dist_func, 
                #    cosmo_functions_for_kernels, self.params_APS_3, self.zv_3
                #)
                #chi_eff_3 = self.params_APS_3["namaster"]["lbinner"](chi_theory_3)[self.WSEL_3]
                #diff_3 = self.chi_data_3 * ah_3 - chi_eff_3
                #chi2_3 = diff_3 @ self.inv_covchi_3 @ diff_3                 

                total_chi2 += chi2_1 #+ chi2_2 #+ chi2_3  
                
            if self.include_gxg_cl:
                #self.params_APS_1["field_1_params"]['omegaHI'] = omegahi     
                #_, cg_theory_1 = self._get_cl_g_theory(
                #    Pk_interp_func, hubble_z_func, comov_dist_func, 
                #    cosmo_functions_for_kernels, self.params_APS_1, self.zv_1
                #)
                #cg_eff_1 = self.params_APS_1["namaster"]["lbinner"](cg_theory_1)[self.WSEL_1]
                #diff_1 = self.cg_data_1 * ag_1 - cx_eff_1 
                #chi2_1 = diff_1 @ self.inv_covcg_1 @ diff_1          
                total_chi2 += 0 #chi2_1 + chi2_2 #+ chi2_3     
            return -0.5 * total_chi2

        except Exception as exc:
            self.log.warning(
                f"[LOGP] Calculation failed at H0={H0_val:.2f}, w={w_val:.3f}, "
                f"omegaHI={omegahi:.5f}, omegam={om_val:.3f}, ombh2(fixed)={ombh2_val:.4f}"
                #f"ax1={ax_1:.3f}, ax2={ax_2:.3f}, ah1={ah_1:.3f}, ah2={ah_2:.3f}. Error: {exc}"
                f"ax1={ax_1:.3f}, ah1={ah_1:.3f}. Error: {exc}"
            )
            import traceback 
            self.log.error(traceback.format_exc()) 
            return -np.inf        

    # ---------------------------------------------------------------------
    def _get_cx_theory(self, Pk_interp_func, hubble_z_func, comov_dist_func, 
                       cosmo_for_kernels, params_aps_i, zv_i):#ask for g_index_i to be used in the kernel
        """
        Compute theoretical C_l values using products from Cobaya's CAMB
        for a single dataset configuration.
        """
        l_grid = params_aps_i["l"]
        #params_aps_i["field_1_params"]['omegaHI_model']        
        #params_aps_i["field_1_params"]['fit_omegaHI']
        #params_aps_i["field_1_params"]['omegaHI']

        # Ensure kernel functions are compatible with cosmo_for_kernels
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
        
        W_f2 = c4ft.kernel_galaxy_clustering_i_vec(
            zout        = zv_i,
            Np          = 5000, 
            pars        = params_aps_i["field_2_params"],
            ibin        = params_aps_i["field_2_params"]["bin_to_use"],
            bias_model  = 'SQRT', 
            bias_value  = 1  , 
            include_rsd = False ,
            include_mb  = False , #<--- it is not well implemented for forecasting yet
            cut_to_bin        = False , 
            binning           = True  ,# use the full galaxy distribution or make the distribution in bins
            mu_rsd            = 0     , #cosine between los and k versor
            growth_rate_f     = None  , 
            path_to_bias_file = None  ,#if bias values are provided by file
            mb_model          = False ,# magnification bias model
            mb_alpha_value    = 0     ,# alpha parameter of magnification bias 
            path_to_mb_file   = None  ,# listed values of mb
            use_input         = 1 ,
            params_interp     = cosmo_for_kernels , # only if use rsd and/or mb (magnification bias)
            get_interp_function = 0,      
            verbose=False,     
            )   

        Cl = []
        for ell in l_grid:
            chi_at_zv = comov_dist_func(zv_i)
            k_values = (ell + 0.5) / (chi_at_zv + EPS) 
            pk_values_on_grid = Pk_interp_func(zv_i, k_values, grid=False, warn=False)
            dCdz  = ( (hubble_z_func(zv_i) / C_LIGHT)*pk_values_on_grid/(chi_at_zv**2 + EPS) )
            dCdz *= W_f1 * W_f2
            Cl.append(integrate.simpson(dCdz, zv_i))
        return l_grid, np.asarray(Cl)
        
    def _get_cl_hi_theory(self, Pk_interp_func, hubble_z_func, comov_dist_func, 
                         cosmo_for_kernels, params_aps_i, zv_i):#ask for g_index_i to be used in the kernel
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
            k_values = (ell + 0.5) / (chi_at_zv + EPS) 
            pk_values_on_grid = Pk_interp_func(zv_i, k_values, grid=False, warn=False)
            dCdz = ( (hubble_z_func(zv_i) / C_LIGHT)*pk_values_on_grid/(chi_at_zv**2 + EPS) )
            dCdz *= W_f1 * W_f1
            Cl.append(integrate.simpson(dCdz, zv_i))
        return l_grid, np.asarray(Cl)        

    def _get_cl_g_theory(self, Pk_interp_func, hubble_z_func, comov_dist_func, 
                         cosmo_for_kernels, params_aps_i, zv_i):#ask for g_index_i to be used in the kernel
        """
        Compute theoretical gxg-Cl values using products from Cobaya's CAMB
        for a single dataset configuration.
        """
        l_grid = params_aps_i["l"]
        
        W_f2 = c4ft.kernel_galaxy_clustering_i_vec(
            zout        = zv_i,
            Np          = 5000, 
            pars        = params_aps_i["field_2_params"],
            ibin        = params_aps_i["field_2_params"]["bin_to_use"],
            bias_model  = 'SQRT', 
            bias_value  = 1  , 
            include_rsd = False ,
            include_mb  = False , #<--- it is not well implemented for forecasting yet
            cut_to_bin        = False , 
            binning           = True  ,# use the full galaxy distribution or make the distribution in bins
            mu_rsd            = 0     , #cosine between los and k versor
            growth_rate_f     = None  , 
            path_to_bias_file = None  ,#if bias values are provided by file
            mb_model          = False ,# magnification bias model
            mb_alpha_value    = 0     ,# alpha parameter of magnification bias 
            path_to_mb_file   = None  ,# listed values of mb
            use_input         = 1 ,
            params_interp     = cosmo_for_kernels , # only if use rsd and/or mb (magnification bias)
            get_interp_function = 0,      
            verbose=False,     
            )  

        Cl = []
        for ell in l_grid:
            chi_at_zv = comov_dist_func(zv_i)
            k_values = (ell + 0.5) / (chi_at_zv + EPS) 
            
            pk_values_on_grid = Pk_interp_func(zv_i, k_values, grid=False, warn=False)

            dCdz = ( (hubble_z_func(zv_i) / C_LIGHT)*pk_values_on_grid/(chi_at_zv**2 + EPS) 
            )
            dCdz *= W_f2 * W_f2
            
            Cl.append(integrate.simpson(dCdz, zv_i))

        return l_grid, np.asarray(Cl)               
