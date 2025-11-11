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

path = '/data2/LSSxHI-CODES/scripts/'
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
#       The parameter must be used to update the params_cosmo and then 
#       the camb.cosmology
# - (4) make sure you are setting _get_camb_results() correctly
# - (5) update "except Exception as exc" in logp()
#
## NEW DATASET -- created
# - (1) set the name of h5 data in the YAML file, likelihoof. It is necessary write a new variable with its: filename, general APS information about the cross correlation, HI bin to cross correlated.
# - (2) Initialize the parameters in the CXLilkelihood class: f.ex.
#     data_2_h5: str
#     aps_2_yaml: str
#     hi_2_index: int = 19
# - (3) load the data with appropriate names and save them in the self.

# PS.: Each valor must be spaced as only one null character after ":"
# PS.: Numpy variables are not accepted
# PS.: variable infinity must be: .inf, +.inf, or -.inf
# PS.: Strings without comma
# -----------------------------------------------------------------------------
class CXLikelihood(Likelihood):
    """Cobaya external likelihood implementing *logp*.

    All heavy I/O is performed once in :py:meth:`initialize`.
    The sampled parameters (*H0*, *ombh2*, *tf*) are received by
    :py:meth:`logp` each time the sampler proposes a point.
    """

    # Declare parameters that will be read from the YAML file.
    # These will be available as attributes of the instance (e.g., self.data_h5)
    # after the base class __init__ is called.
    # Default values specified here will be used if the option is not in the YAML.
    ## -- params from LIKELIHOOD key (YAML file) --
    # new addition must be implemented here, updated in the initialize(), and check where they are evoked in logp()
    data_1_h5: str
    aps_1_yaml: str
    hi_1_index: int = 19
    num_z: int = 100
    lmin: int = 30
    lmax: int = 400
    del_l: int = 20
    nside: int = 256

    # ---------------------------------------------------------------------
    # 0. Heavy initialisation (runs once before sampling) ------------------
    # ---------------------------------------------------------------------
    def initialize(self):
        """
        This method is called once after __init__ when the likelihood is loaded.
        Use it for heavy setup, loading data, etc.
        Parameters defined above (data_h5, aps_yaml, etc.) are now
        available as attributes of `self`.
        """
        # ---- Paths to input files ----
        # os.path.expanduser is good practice for paths from YAML.
        # Cobaya might do some path resolution, but being explicit is safe.
        data_1_h5_path = os.path.expanduser(self.data_1_h5)
        aps_1_yaml_path = os.path.expanduser(self.aps_1_yaml)

        # ---- load observables -------------------------------------------
        with h5py.File(data_1_h5_path, "r") as f:
            self.l_data     = f["l"][:]
            self.cl_data    = f["cl"][:]
            self.errcl_data = f["errcl"][:]
        self.inv_cov = np.diag(1.0 / self.errcl_data**2)

        # ---- load auxiliary configuration -------------------------------
        p_loaded         = c4ft.load_yaml_fileformat(aps_1_yaml_path, verbose=0)
        self.params_APS   = dcopy(p_loaded["APS"])
        self.params_cosmo = dcopy(p_loaded["cosmology"])
        del p_loaded

        # ---- NaMaster binning -------------------------------------------
        self.b    = nmt.NmtBin.from_nside_linear(self.nside, nlb=self.del_l)
        self.leff = self.b.get_effective_ells()
        self.WSEL = (self.leff > self.lmin) & (self.leff < self.lmax)
        self.params_APS["namaster"] = {
            "b": self.b,
            "leff": self.leff,
            "lbinner": self.b.bin_cell,
        }

        # ---- red‑shift grid ---------------------------------------------
        zf1  = self.params_APS["field_1_params"]["z"]
        zmin = max(zf1[self.hi_1_index : self.hi_1_index + 2].min(),
                   self.params_APS["field_2_params"]["z_min"])
        zmax = min(zf1[self.hi_1_index : self.hi_1_index + 2].max(),
                   self.params_APS["field_2_params"]["z_max"])
        self.zv = np.logspace(np.log10(zmin), np.log10(zmax), self.num_z)

        self.log.info(f"CXLikelihood initialized with lmin={self.lmin}, lmax={self.lmax}")

    # ---------------------------------------------------------------------
    # 1. Inform Cobaya of any external requirements -----------------------
    # ---------------------------------------------------------------------
    def get_requirements(self):
        """
        Specifies requirements from other components.
        For sampled parameters used directly in logp (H0, ombh2, tf),
        their presence in the logp signature should ideally be enough.
        However, explicitly listing them here can resolve detection issues
        if Cobaya isn't picking them up automatically from the logp signature.
        """
        return {
                "H0": None, "ombh2": None, "tf_1": None,
                "omch2": None, "ns": None, "As": None, "mnu": None,
               }

    # ---------------------------------------------------------------------
    # 2. Log‑likelihood ----------------------------------------------------
    # ---------------------------------------------------------------------
    def logp(self, H0, ombh2, tf_1, **kwargs):
        """
        Return **log P(data | parameters)** for the sampler.
        """
        # The _derived dictionary can be accessed via kwargs.get("_derived")
        # if you need to store derived parameters. For example:
        # derived_params_container = kwargs.get("_derived", {})
        params_cosmo = dcopy(self.params_cosmo)
        params_cosmo.update(
            H0=H0,
            ombh2=ombh2,
            omch2=kwargs["omch2"],
            ns=kwargs["ns"],
            As=kwargs["As"],
            mnu=kwargs["mnu"],
        )        

        try:
            pars_camb    = self._get_camb_results(params_cosmo)
            _, cl_theory = self._get_cx_theory(pars_camb)
            
            cl_eff       = self.params_APS["namaster"]["lbinner"](cl_theory)[self.WSEL]
            diff         = self.cl_data * tf_1 - cl_eff 
            chi2         = diff @ self.inv_cov @ diff
            
            return -0.5 * chi2

        except Exception as exc:
            self.log.warning(
                "[CXLikelihood] Calculation failed at H0=%g, ombh2=%g, tf=%g. Error: %s",
                H0, ombh2, tf_1, exc
            )
            return -np.inf

    # ---------------------------------------------------------------------
    # 3. Helpers -----------------------------------------------------------
    # ---------------------------------------------------------------------
    def _get_camb_results(self, pars): # pars is params_cosmo from logp
        """Internal method to run CAMB and get results."""
        # (H0, ombh2, omch2, As, ns, mnu, halofit_version, z, kmax, Pk_nonlinear)
        #current_As = pars.get( "As" , 2e-9)   # Default if not in pars (from self.params_cosmo)
        #current_mnu = pars.get("mnu", 0.06) # ''

        cpars = camb.set_params(halofit_version=pars["halofit_version"])
        cpars.set_cosmology(
            H0=pars["H0"], 
            ombh2=pars["ombh2"], 
            omch2=pars["omch2"], 
            mnu=pars["mnu"]
        )
        cpars.InitPower.set_params(As=pars["As"], ns=pars["ns"])
        cpars.set_matter_power(
            redshifts=pars["z"],
            kmax=pars["kmax"], 
            silent=True
        )

        cres = camb.get_results(cpars)
        pk_interp = camb.get_matter_power_interpolator(
            cpars, 
            nonlinear=pars["Pk_nonlinear"], 
            hubble_units=False, 
            k_hunit=False,
            kmax=pars["kmax"], 
            zmax=pars["z"][-1], 
        )
        
        kc,zc,pkc = cres.get_nonlinear_matter_power_spectrum(have_power_spectra=0, 
                                                             hubble_units=False, 
                                                             k_hunit=False, 
                                                             params=cpars)
                
        zstar = cres.get_derived_params()["zstar"]
        return {
            'k':kc,'z':zc,'pkc':pkc, 
            "pk": pk_interp.P, 
            "z_camb": pars["z"],
            "comovel_dist_z": cres.comoving_radial_distance, 
            "hubble_z": cres.hubble_parameter,
            "fsigma8": cres.get_fsigma8(),
            "sigma8": cres.get_sigma8(), 
            "zstar": zstar,
            "chistar": cres.comoving_radial_distance(zstar),
        }

    def _get_cx_theory(self, params_camb): # receives the dictionary from _get_camb_results
        """Internal method to compute theoretical C_l values."""
        l_grid     = self.params_APS["l"]      
        hubble_z   = params_camb["hubble_z"]   
        Pk_interp  = params_camb["pk"]         
        comov_dist = params_camb["comovel_dist_z"]

        W_f1 = c4ft.kernel_hi_vec(
            zout=self.zv, zrange=self.zv, 
            replace_inf2nan=True, use_input=True,
            config_dict=self.params_APS["field_1_params"], 
            params_interp=params_camb, 
        )
        W_f2 = c4ft.kernel_galaxy_clustering_i_vec(
            zout=self.zv,
            pars=self.params_APS["field_2_params"],
            ibin=self.params_APS["field_2_params"]["bin_to_use"],
            verbose=False, # Keep verbose off during sampling
        )

        Cl = []
        for ell in l_grid:
            chi_zv = comov_dist(self.zv)
            k_values = (ell + 0.5) / chi_zv 
            pk_values_on_grid = Pk_interp(self.zv, k_values, grid=False)
            dCdz = ( (hubble_z(self.zv)/C_LIGHT)*pk_values_on_grid/(chi_zv**2) )
            dCdz *= W_f1 * W_f2
            Cl.append(integrate.simpson(dCdz, self.zv))

        return l_grid, np.asarray(Cl)
