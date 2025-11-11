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
    data_h5_1: str
    aps_yaml_1: str
    hi_index_1: int = 19
    num_z_1: int = 100
    lmin_1: int = 30
    lmax_1: int = 400
    del_l_1: int = 20
    
    data_h5_2: str
    aps_yaml_2: str
    hi_index_2: int = 20 # Default for dataset 2, can be overridden in YAML
    num_z_2: int = 100
    lmin_2: int = 30
    lmax_2: int = 400
    del_l_2: int = 20
    
    nside: int = 256

    # ---------------------------------------------------------------------
    def initialize(self):
        """
        Initialize likelihood: load data for two datasets, 
        and fixed settings from their respective YAML files.
        """
        # ---- Paths to input files ----
        data_h5_path_1 = os.path.expanduser(self.data_h5_1)
        aps_yaml_path_1 = os.path.expanduser(self.aps_yaml_1)
        data_h5_path_2 = os.path.expanduser(self.data_h5_2)
        aps_yaml_path_2 = os.path.expanduser(self.aps_yaml_2)
        
        # ---- Load observables for dataset 1 ----
        with h5py.File(data_h5_path_1, "r") as f:
            self.l_data_1     = f["l"][:]
            self.cl_data_1    = f["cl"][:]
            self.errcl_data_1 = f["errcl"][:]
        self.inv_cov_1 = np.diag(1.0 / self.errcl_data_1**2)

        # ---- Load observables for dataset 2 ----
        with h5py.File(data_h5_path_2, "r") as f:
            self.l_data_2     = f["l"][:]
            self.cl_data_2    = f["cl"][:]
            self.errcl_data_2 = f["errcl"][:]
        self.inv_cov_2 = np.diag(1.0 / self.errcl_data_2**2)

        # ---- Load auxiliary configuration for dataset 1 ----
        p_1_loaded         = c4ft.load_yaml_fileformat(aps_yaml_path_1, verbose=0)
        self.params_APS_1  = dcopy(p_1_loaded["APS"])
        # self.params_cosmo_1 = dcopy(p_1_loaded["cosmology"]) # CAMB params now global
        del p_1_loaded

        # ---- Load auxiliary configuration for dataset 2 ----
        p_2_loaded         = c4ft.load_yaml_fileformat(aps_yaml_path_2, verbose=0)
        self.params_APS_2  = dcopy(p_2_loaded["APS"])
        # self.params_cosmo_2 = dcopy(p_2_loaded["cosmology"]) # CAMB params now global
        del p_2_loaded

        # ---- NaMaster binning for dataset 1 ----
        self.b_1    = nmt.NmtBin.from_nside_linear(self.nside, nlb=self.del_l_1)
        self.leff_1 = self.b_1.get_effective_ells()
        self.WSEL_1 = (self.leff_1 > self.lmin_1) & (self.leff_1 < self.lmax_1)
        self.params_APS_1["namaster"] = {
            "b": self.b_1, "leff": self.leff_1, "lbinner": self.b_1.bin_cell,
        }

        # ---- NaMaster binning for dataset 2 ----
        self.b_2    = nmt.NmtBin.from_nside_linear(self.nside, nlb=self.del_l_2)
        self.leff_2 = self.b_2.get_effective_ells()
        self.WSEL_2 = (self.leff_2 > self.lmin_2) & (self.leff_2 < self.lmax_2)
        self.params_APS_2["namaster"] = {
            "b": self.b_2, "leff": self.leff_2, "lbinner": self.b_2.bin_cell,
        }

        # ---- Redshift grid for dataset 1 integration ----
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
        zf1_2  = self.params_APS_2["field_1_params"]["z"]
        zmin_2 = max(zf1_2[self.hi_index_2 : self.hi_index_2 + 2].min(),
                     self.params_APS_2["field_2_params"]["z_min"])
        zmax_2 = min(zf1_2[self.hi_index_2 : self.hi_index_2 + 2].max(),
                     self.params_APS_2["field_2_params"]["z_max"])
        if zmin_2 <= 0: raise ValueError("zmin_2 for logspace must be positive.")
        if zmax_2 <= zmin_2: raise ValueError("zmax_2 must be > zmin_2 for dataset 2.")
        self.zv_2 = np.logspace(np.log10(zmin_2), np.log10(zmax_2), self.num_z_2)
        #self.log.info(f"{self.name} initialized for dataset 2. Integration z: [{zmin_2:.4f}, {zmax_2:.4f}]")
        self.log.info(f"Dataset 2 initialized for dataset 2. Integration z: [{zmin_2:.4f}, {zmax_2:.4f}]")
    # ------------------------------------------------------------------
    def get_requirements(self):
        """Exact cosmological outputs needed from CAMB."""

        # Combine both datasets’ redshifts *and* z=0 for background quantities
        z_full = np.unique(np.concatenate([self.zv_1, self.zv_2]))
        z_bg   = np.unique(np.concatenate(([0.0], z_full)))   # adds 0

        k_max_needed = 50.0  # h Mpc⁻¹

        return {
            "Pk_interpolator": {
                "z": z_full.tolist(),    # P(k,z) need not include 0
                "k_max": k_max_needed,
                "nonlinear": True,
                #"var_type": "delta_m",
            },
            "Hubble": {"z": z_bg.tolist()},
            "comoving_radial_distance": {"z": z_bg.tolist()},
            "sigma8": None,
            "tf_1": None,
            "tf_2": None,
        }
    # ---------------------------------------------------------------------
    def logp(self, tf_1, tf_2, **kwargs):
        """Calculate the log-likelihood using products from Cobaya's CAMB."""

        # Retrieve current cosmology from the provider if needed
        H0   = self.provider.get_param("H0")      # optional
        ombh2 = self.provider.get_param("ombh2")  # optional
        try:
            # Pk_interpolator is requested with nonlinear=True in get_requirements
            Pk_interp_func = self.provider.get_Pk_interpolator() 
            hubble_z_func = self.provider.get_Hubble 
            comov_dist_func = self.provider.get_comoving_radial_distance
            sigma8_val = self.provider.get_param("sigma8")
        except Exception as e:
            #self.log.error(f"[{self.name}] Failed to get products from provider: {e}")
            self.log.error(f"[LOGP] Failed to get products from provider: {e}")
            return -np.inf

        # Prepare a dictionary of cosmological functions/parameters for your kernels
        # This needs to be compatible with what `kernel_hi_vec` expects for `params_interp`.
        cosmo_functions_for_kernels = {
            "hubble_z": hubble_z_func,
            "comovel_dist_z": comov_dist_func,
            "sigma8": sigma8_val, # Pass current sigma8 value
            # Add any other derived quantities or functions your kernels might need.
            # E.g., if kernels need H0 or ombh2 explicitly:
            # "H0": H0, "ombh2": ombh2 
        }

        try:
            # --- Calculate chi2 for dataset 1 ---
            _, cl_theory_1 = self._get_cx_theory(
                Pk_interp_func, hubble_z_func, comov_dist_func, 
                cosmo_functions_for_kernels, self.params_APS_1, self.zv_1
            )
            cl_eff_1 = self.params_APS_1["namaster"]["lbinner"](cl_theory_1)[self.WSEL_1]
            diff_1 = self.cl_data_1 * tf_1 - cl_eff_1 
            chi2_1 = diff_1 @ self.inv_cov_1 @ diff_1
            
            # --- Calculate chi2 for dataset 2 ---
            _, cl_theory_2 = self._get_cx_theory(
                Pk_interp_func, hubble_z_func, comov_dist_func, 
                cosmo_functions_for_kernels, self.params_APS_2, self.zv_2
            )
            cl_eff_2 = self.params_APS_2["namaster"]["lbinner"](cl_theory_2)[self.WSEL_2]
            diff_2 = self.cl_data_2 * tf_2 - cl_eff_2
            chi2_2 = diff_2 @ self.inv_cov_2 @ diff_2
            
            total_chi2 = chi2_1 + chi2_2
            return -0.5 * total_chi2

        except Exception as exc:
            self.log.warning(
                #f"[{self.name}] Calculation failed at H0={H0}, ombh2={ombh2}, "
                f"[LOGP] Calculation failed at H0={H0}, ombh2={ombh2}, "
                f"tf1={tf_1}, tf2={tf_2}. Error: {exc}"
            )
            # For debugging, uncomment to see full traceback:
            # import traceback
            # self.log.error(traceback.format_exc())
            return -np.inf

    # ---------------------------------------------------------------------
    def _get_cx_theory(self, Pk_interp_func, hubble_z_func, comov_dist_func, 
                       cosmo_for_kernels, params_aps_i, zv_i):
        """
        Compute theoretical C_l values using products from Cobaya's CAMB
        for a single dataset configuration.
        """
        l_grid = params_aps_i["l"]

        # Ensure kernel functions are compatible with cosmo_for_kernels
        W_f1 = c4ft.kernel_hi_vec(
            zout=zv_i, zrange=zv_i, 
            replace_inf2nan=True, use_input=True,
            config_dict=params_aps_i["field_1_params"], 
            params_interp=cosmo_for_kernels,
        )
        W_f2 = c4ft.kernel_galaxy_clustering_i_vec(
            zout=zv_i,
            pars=params_aps_i["field_2_params"],
            ibin=params_aps_i["field_2_params"]["bin_to_use"],
            verbose=False, 
        )

        Cl = []
        for ell in l_grid:
            chi_at_zv = comov_dist_func(zv_i)
            k_values = (ell + 0.5) / (chi_at_zv + 1e-30) 
            
            pk_values_on_grid = Pk_interp_func(zv_i, k_values, grid=False, warn=False)

            dCdz = (
                (hubble_z_func(zv_i) / C_LIGHT) *
                pk_values_on_grid /
                (chi_at_zv**2 + 1e-30) 
            )
            dCdz *= W_f1 * W_f2
            
            Cl.append(integrate.simpson(dCdz, zv_i))

        return l_grid, np.asarray(Cl)