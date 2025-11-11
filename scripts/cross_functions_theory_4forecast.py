#############################################
#PATH TO SCRIPTS FOLDER
import cross_functions_theory      as cxft
import handling_data               as hdata
import noise_functions             as fnoise
from copy import deepcopy as dcopy
import numpy as np
import time, camb
from scipy.interpolate import interp1d
from scipy import integrate
C_LIGHT = 299792458/1e3 #km/s
nu_HI   = 1420.405751768 #MHz
#####################################################################################################
#### HI FUNCTIONS
def z_hi_survey(hi_survey='LOWZ'):
    if hi_survey.upper()=='LOWZ':
        nu_vec = np.linspace(980,1260,30+1) 
    elif hi_survey.upper()=='HIGHZ':
        nu_vec = np.linspace(350,1050,70+1)
    elif hi_survey.upper()=='MIDZ-2':
        nu_vec = np.linspace(470,1420,100+1)    
    elif hi_survey.upper()=='ABDUS':
        nu_vec = np.linspace(450,1420,100+1)            
    else:
        raise NameError('No survey named {}'.format(hi_survey))
    return np.flip((nu_HI/nu_vec)-1)

def hi_brightness_temperature(camb_params=None, camb_results=None, z=None, model = None, fit_OHI=False, OHI=None, use_input=False, params_interp=None ):
    if use_input:
        fact = (1+z)**2/(params_interp['hubble_z'](z)/params_interp['hubble_z'](0))
        h_lower = params_interp['hubble_z'](0)/100
    else:
        fact = (1+z)**2/(camb_results.hubble_parameter(z)/camb_params.H0)
        h_lower = camb_params.h                    
        #try:#if camb_* are CAMB instance;
        #    fact = (1+z)**2/(camb_results.hubble_parameter(z)/camb_params.H0)
        #    h_lower = camb_params.h            
        #except: #otherwise, theyre dictionaries
        #    fact = (1+z)**2/(camb_results['hubble_z'](z)/camb_params['H0'])
        #    h_lower = camb_params['H0']/100       
    if not fit_OHI or type(OHI)==type(None):     
        model = None
        oHI   = omegaHI_biasHI(model=model,z=z)['omegaHI']*h_lower
    else: 
        oHI = OHI*h_lower
    fact = fact*188.8
    return fact*oHI #mK


def omegaHI_biasHI(model=None,z=None, OHI=None, bHI=None):
    """
    ----------------------------------------------    
    model: Omega_HI,bias_HI model to be used.
           - Jiajun-hod (Zhang et al 2022)
           - Jiajun-ham (Zhang et al 2022)
           - padmanabhan (Padmanabhan 2015)
           - cunnington (Cunnington 2019)
           - bull2015 (P. Bull 2015)
           - crighton (I Carucci 2020)
           - 'fixed','cte','const' or 'constant': must provide OHI and/or bHI
    z: numpy.array of redshift values to be used
    ----------------------------------------------
    """
    if model in ["jiajun","jiajun-hod"]: #valido para z < 1.#Zhang et al 2022
        return {"omegaHI":np.array([2.7e-4 + 1e-4*z - 8e-5*z**2]),
                "biasHI" :1*np.ones_like(z)
               }
    if model=="jiajun-ham": #valido para z < 1.#Zhang et al 2022    
        return {"omegaHI":np.array([2.5e-4 - 4e-5*z - 7e-5*z**2]),
                "biasHI" :1*np.ones_like(z)
               }                
    elif model=="padmanabhan":#padmanabhan 2015
        zHI = np.array(  [0.000, 0.250, 0.500, 0.750, 1.000, 1.250, 1.500, 1.750, 2.000, 2.250, 2.500, 2.750, 3.000, 3.250, 3.400])
        oHI = np.array([3.344, 3.443, 4.523, 4.648, 4.710, 4.804, 4.766, 4.804, 4.936, 5.008, 4.750, 5.471, 5.541, 5.756, 5.971])
        bHI = np.array([0.703, 0.972, 1.026, 0.935, 1.005, 1.005, 1.049, 1.099, 1.101, 1.160, 1.261, 1.409, 1.329, 1.498, 1.802])        
        #return {"z":z, "omegaHI":oHI, "bias":bHI}
        return {"z":z, "omegaHI":np.interp(z, zHI, oHI), "bias":np.interp(z, zHI, bHI)}
    elif model=="cunnington":#cunnington2019
        return {"omegaHI":0.00048 + 0.00039*z - 0.000065*z**2}
    elif model=="irfan":#Irfan2021 baseado em bull2015 #o mesmo que IRFAN2021
        omegaHI_0 = 4.86e-4*np.ones_like(z) #Irfan2021
        biasHI_0 = 0.677105*np.ones_like(z) #Irfan2021
        return {"omegaHI": np.array([ (omegaHI_0/4.86)*(4.8304 + 3.8856*z - 0.65119*z**2) ]) ,
                "biasHI" : np.array([ (biasHI_0/0.677105)*(0.66655 + 0.17765*z + 0.050223*z**2) ])}
    elif model=="crighton":#carucci2020
        omegaHI_0 = 4e-4*np.ones_like(z) #Irfan2021
        biasHI_0  = 0.3*np.ones_like(z) #Irfan2021        
        return {"omegaHI": np.array([ omegaHI_0*(1+z)**(0.6) ]),
                "biasHI" : np.array([biasHI_0*(1+z)+0.6 ])}
    elif model in ['fixed','cte','const','constant']:
        if OHI is None: OHI= 4.86e-4*np.ones_like(z) 
        else:
            if (type(OHI)==int) or (type(OHI)==float): 
                OHI = OHI*np.ones_like(z) 
            elif (type(OHI) == list) and (np.array(OHI).size==1):
                OHI = OHI*np.ones_like(z) 
            elif (type(OHI) == list) and (np.array(OHI).size==z.size):
                pass
            else:
                OHI= 4.86e-4*np.ones_like(z) 
        if bHI is None: 1*np.ones_like(z) 
        else:
            if (type(bHI)==int) or (type(bHI)==float): 
                bHI = bHI*np.ones_like(z) 
            elif (type(bHI) == list) and (np.array(bHI).size==1):
                bHI = bHI*np.ones_like(z) 
            elif (type(bHI) == list) and (np.array(bHI).size==z.size):
                pass
            else:
                bHI= 1*np.ones_like(z)  
        return {"omegaHI": OHI,
                "biasHI" : bHI}             
    else:
        return {"omegaHI": 4.86e-4*np.ones_like(z),
                "biasHI" : 1*np.ones_like(z)}             
def kernel_hi(camb_params=None, camb_results=None, z=None, zrange=None, config_dict=None, replace_inf2nan=True, use_input=False, params_interp=None, decimals=6 ):
    #model=None, fit_bHI=False, bHI=None, fit_OHI=False, OHI=None):
    """
    HI-KERNEL for one z value    
    zrange: array of z bins of the survey (not the whole z analyzed)
    z: value to be analyzed (in/out zrange)

    if you choice to determine (fit) HI bias (bHI) and Omega HI (OHI), 
    you must to specify True for the fit_bHI (fit_OHI) parameter and to set the parameter value
    otherwise, the code will assume the model provided.    
    """
    omegaHI_model = config_dict['omegaHI_model']
    biasHI_model  = config_dict['biasHI_model']
    fit_bHI       = config_dict['fit_biasHI']
    fit_OHI       = config_dict['fit_omegaHI']
    bHI           = config_dict['biasHI']
    OHI           = config_dict['omegaHI']
    #if you choice to determine (fit) HI bias (bHI) and Omega HI (OHI), 
    #you must to specify True for the fit_bHI (fit_OHI) parameter and to set the parameter value
    #otherwise, the code will assume the model provided.
    z  = np.around(z, decimals=decimals)
    w1 = cxft.window_2(value=z, arr=zrange)#selection function
    if not w1: return 0
    THI = hi_brightness_temperature(camb_params=camb_params, camb_results=camb_results, z=z, 
                                    model=omegaHI_model, fit_OHI=fit_OHI, OHI=OHI, 
                                    use_input=use_input, params_interp=params_interp )
    if not fit_bHI or type(bHI)==type(None): 
        bHI = omegaHI_biasHI(model=biasHI_model,z=z)['biasHI']
    else:
        pass
    print('THI: | ',THI)
    if replace_inf2nan and (bHI*THI*w1 == np.inf): return np.nan
    else: return bHI*THI*w1 

def kernel_hi_vec(camb_params=None, camb_results=None, 
                  zout=None, zrange=None, 
                  config_dict=None, use_input=False, 
                  params_interp=None, decimals=6, 
                  replace_inf2nan=True, get_interp_function=False):

    """
    HI-KERNEL vetorized (or its interpolation)
    
    Args:
        zrange              :numpy array  : Edge redshifts of the HI IM survey  (standard value: None)
        zout                :numpy array  : redshifts to be provide the HI-kernel from  (standard value: None)    
        camb_params         :CAMB instance: CAMB parameters instance (standard value: None)
        camb_results        :CAMB instance: CAMB results instance    (standard value: None)    

        config_dict         :dict: provide the HI information such as omega_HI model, bias_HI model, ...  (standard value: None)
        params_interp       :dict: ------    (standard value: None)
        use_input           :bool: inform weather it would be used camb instabce or function from the dict params_interp  (standard value: None)
        decimals            :int:  decimal level to be calculate the redshift within the window function  (standard value: 6)
        get_interp_function :bool: to obtain the interpolation instead of   (standard value: False)
        replace_inf2nan     :bool: to replace infinity values to numpy NaN   (standard value: True)
        
    Returns:
           Window function:
               numpy array , if get_interp_function=False
               scipy interp, if get_interp_function=True
        
    Ps.:          
        If you choice to determine (fit) HI bias (bHI) and Omega HI (OHI), 
        you must to specify True for the fit_bHI (fit_OHI) parameter and to set the parameter value
        otherwise, the code will define the standard model.   
    """
    from scipy.interpolate import interp1d
    omegaHI_model = config_dict['omegaHI_model']
    biasHI_model  = config_dict['biasHI_model']
    fit_bHI       = config_dict['fit_biasHI']
    fit_OHI       = config_dict['fit_omegaHI']
    bHI           = config_dict['biasHI']
    OHI           = config_dict['omegaHI']
    
    zrange_ext = []
    for i in range(zrange[:-1].size):
        zrange_ext.append( np.linspace(zrange[i], zrange[i+1], 100) ) ## |delta_zext|<1e-4 for zrange in [0.12,0.45]
    zrange_ext = np.array(  zrange_ext ).flatten()    
    zrange_ext = np.unique( zrange_ext )
    
    w1  = np.array([ cxft.window_2(value=np.around(iz,decimals), arr=zrange) for iz in zrange_ext ]) #selection function
    THI = hi_brightness_temperature(camb_params=camb_params, 
                                         camb_results=camb_results, z=zrange_ext, 
                                         model=omegaHI_model, fit_OHI=fit_OHI, OHI=OHI, 
                                         use_input=use_input, params_interp=params_interp )
    if not fit_bHI or type(bHI)==type(None): 
        bHI = omegaHI_biasHI(model=biasHI_model,z=zrange_ext)['biasHI']
    else:
        pass
    arr = bHI*THI*w1    
    if replace_inf2nan:
        #arr = bHI*THI*w1
        arr_interp = interp1d(x=zrange_ext, y=arr, kind='cubic', fill_value=np.nan)
        arr        = arr_interp(zrange)
        #arr[np.isinf(arr)] = np.nan 
        #return arr
    else: 
        arr_interp = interp1d(x=zrange_ext, y=arr, kind='cubic')
        arr        = arr_interp(zrange)
        
    if not get_interp_function:
        if zout is None:
            return arr#bHI*THI*w1     
        else:
            return arr_interp(zout)#np.interp(zout, zrange, arr)
    else:
        
        return arr_interp#interp1d(x=zrange, y=arr, kind='cubic', fill_value=np.nan)
        

def kernel_galaxy_clustering_i_dict(camb_params       = None  , 
                                    camb_results      = None  , 
                                    pars              = None  , 
                                    cut_to_bin        = False , #binning=True,
                                    ibin              = 0     , 
                                    Np                = 500   , 
                                    binning           = True  ,
                                    include_rsd       = False , 
                                    include_mb        = False , 
                                    mu_rsd            = 0     , 
                                    growth_rate_f     = None  , 
                                    bias_model        = 'SQRT', 
                                    bias_value        = 1     , 
                                    path_to_bias_file = None  ,
                                    mb_model          = False , 
                                    mb_alpha_value    = 0     , 
                                    path_to_mb_file   = None  ,
                                    use_input         = False , 
                                    params_interp     = None  ,
                                    verbose           = False  ):
    #c_light = 299792458/1e3
    #EPS=1e-50
    from scipy.interpolate import interp1d
    if verbose: timei = time.time()     
    #
    delz = np.abs(pars['z_max']-pars['z_min'])/Np
    z_range   = np.linspace(pars['z_min']-delz, pars['z_max']+delz,Np+2)
    bias_gc   = cxft.get_bias_model(bias_model=bias_model, 
                                    bias_value=bias_value, 
                                    pars=pars, 
                                    path_to_file=path_to_bias_file, 
                                    z_range=z_range)
    dNdz_vec  = cxft.dNdz_func(pars=pars, 
                               binning=binning,
                               Np=10*Np,
                               type_field=pars['field'],
                               verbose=False)
    dNdz_intp = interp1d(dNdz_vec['z'], dNdz_vec['dNdz'][ibin,:], kind='cubic')
    kernel_gc = np.interp(z_range[1:-1], dNdz_vec['z'], dNdz_vec['dNdz'][ibin,:])
    #kernel_gc = np.array([ dNdz_intp(zi) for zi in z_range ])
    ######################################################################
    if include_rsd:
        if use_input:
            z_input = params_interp['z']
            fsigma8 = params_interp['fsigma8'](z_input)
            sigma8  = params_interp[ 'sigma8'](z_input)
        else:
            z_input = camb_results.z_input
            fsigma8 = camb_results.get_fsigma8()
            sigma8  = camb_results.get_sigma8()
        f_vec = interp1d(z_input, fsigma8/get_sigma8, kind='cubic', fill_value='extrapolate')(z_range[1:-1])
        kernel_gc+=(f_vec*mu_rsd**2)*kernel_gc
    if 0:#include_mb:
        MB = cxft.magnification_bias(camb_params=camb_params, camb_results=camb_results, pars=pars,
                                slope_s_model=mb_model, z_range=z_range,
                                binning=binning, ibin=ibin, Np=Np, verbose=False)
        kernel_gc+=MB
    ######################################################################        
    kernel_gc_intp = interp1d(z_range[1:-1], kernel_gc, kind='cubic')#, fill_value='extrapolate')
    if verbose: print('Total processing time: {0:.2f} seg'.format(time.time()-timei))             
    return {'W':kernel_gc, 'z':z_range[1:-1], 'W_interp':kernel_gc_intp, 'bin':ibin}      

def kernel_galaxy_clustering_i_vec( camb_params       = None  , 
                                    camb_results      = None  , 
                                    pars              = None  ,
                                    zout              = None  , 
                                    cut_to_bin        = False , #binning=True,
                                    ibin              = 0     , 
                                    Np                = 5000  , 
                                    binning           = True  ,
                                    include_rsd       = False , 
                                    include_mb        = False , 
                                    mu_rsd            = 0     , 
                                    growth_rate_f     = None  , 
                                    bias_model        = 'SQRT', 
                                    bias_value        = 1     , 
                                    path_to_bias_file = None  ,
                                    mb_model          = False , 
                                    mb_alpha_value    = 0     , 
                                    path_to_mb_file   = None  ,
                                    use_input         = False , 
                                    params_interp     = None  ,
                                    get_interp_function = False,
                                    kind_interpolator = 'cubic',
                                    extrap_interpolator = True,
                                    get_dict_format   = False,
                                    verbose           = False ):
    """
    Calculates theoretical Galaxy Clustering Kernel

    Args:
        camb_results (dict): Dictionary containing CAMB outputs like
                                    'pk' (Pk interpolator P(z,k)),
                                    'hubble_z' (H(z) function),
                                    'comovel_dist_z' (chi(z) function), etc.
                                    Matches output of old get_camb_results.
        camb_params (dict): Dictionary with 
        pars (dict): Dictionary with ...
        cut_to_bin (bool): 
        ibin (int): Integer defining the Galaxy redshift bin
        Np (int): Number of points to
        binning (bool): divide or not the survey redshift coverage
        get_interp_function (bool): to obtain the interpolation instead of   (standard value: False)        
        ....
        zout (array): Redshift array for integration.

    Returns:
        array: W kernel
    """    
    #c_light = 299792458/1e3
    from scipy.interpolate import interp1d
    if verbose: timei = time.time()     
    #
    delz      = np.abs(pars['z_max']-pars['z_min'])/Np
    z_range   = np.linspace(pars['z_min']-delz, pars['z_max']+delz,Np+2)
    bias_gc   = cxft.get_bias_model(bias_model=bias_model, 
                                    bias_value=bias_value, 
                                    pars=pars, 
                                    path_to_file=path_to_bias_file, 
                                    z_range=z_range)
    dNdz_vec  = cxft.dNdz_func(pars=pars, 
                               binning=binning,
                               Np=10*Np,
                               type_field=pars['field'],
                               verbose=False)    
    kernel_gc = np.interp(z_range[1:-1], dNdz_vec['z'], dNdz_vec['dNdz'][ibin,:])
    kernel_gc*=bias_gc[1:-1]
    if include_rsd:
        if use_input:
            z_input = params_interp['z']
            fsigma8 = params_interp['fsigma8'](z_input)
            sigma8  = params_interp[ 'sigma8'](z_input)
        else:
            z_input = camb_results.z_input
            fsigma8 = camb_results.get_fsigma8()
            sigma8  = camb_results.get_sigma8()
        f_vec = interp1d(z_input, fsigma8/get_sigma8, kind='cubic', fill_value='extrapolate')(z_range[1:-1])
        #f_vec = interp1d(z_input, fsigma8/get_sigma8, kind='cubic', fill_value='extrapolate')(z_range)
        #dNdz_vec+=(f_vec*mu_rsd**2)*dNdz_vec
        kernel_gc+=(f_vec*mu_rsd**2)*kernel_gc
    if 0:#include_mb:
        MB = cxft.magnification_bias(camb_params=camb_params, 
                                     camb_results=camb_results, 
                                     pars=pars,
                                     slope_s_model=mb_model, 
                                     z_range=z_range,
                                     binning=binning, ibin=ibin, Np=Np, 
                                     verbose=False)
        MB = interp1d(z_range, MB, kind='cubic', fill_value='extrapolate')(z_range[1:-1])
        #dNdz_vec+=MB
        kernel_gc+=MB
    if not get_interp_function:
        if type(zout)==type(None):# or np.around(zout,decimals=6)==np.around(z_range[1:-1],decimals=6): 
            if verbose: print('Total processing time: {0:.2f} seg'.format(time.time()-timei))             
            if get_dict_format:
                return {'z':z_range[1:-1], 'W':kernel_gc}
            else:
                return kernel_gc
        else:
            if verbose: print('Total processing time: {0:.2f} seg'.format(time.time()-timei))             
            if get_dict_format:
                return {'z':z_range[1:-1], 'W':np.interp(x=zout, xp=z_range[1:-1], fp= kernel_gc) }
            else:
                return np.interp(x=zout, xp=z_range[1:-1], fp= kernel_gc) 
    else:
        if verbose: 
            print('Total processing time: {0:.2f} seg'.format(time.time()-timei))        
        if extrap_interpolator:
            if get_dict_format:
                return {'z':z_range[1:-1], 'W':interp1d(z_range[1:-1], kernel_gc, kind=kind_interpolator,  fill_value="extrapolate")}
            else:
                return interp1d(z_range[1:-1], kernel_gc, kind=kind_interpolator,  fill_value="extrapolate")
        else:
            if get_dict_format:
                return {'z':z_range[1:-1], 'W':interp1d(z_range[1:-1], kernel_gc, kind=kind_interpolator,  fill_value=np.nan)}
            else:
                return interp1d(z_range[1:-1], kernel_gc, kind=kind_interpolator,  fill_value=np.nan)       
        
def get_cx_theory(zv, params_camb, params_aps, kind_integ='simpson', quad_limit=100): 
    """
    Calculates theoretical C(l) using pre-calculated CAMB results.

    Args:
        params_camb_results (dict): Dictionary containing CAMB outputs like
                                    'pk' (Pk interpolator P(z,k)),
                                    'hubble_z' (H(z) function),
                                    'comovel_dist_z' (chi(z) function), etc.
                                    Matches output of old get_camb_results.
        params_APS (dict): Dictionary with settings for l-values, fields, etc.
        zv (array): Redshift array for integration.

    Returns:
        tuple: (l_values, C_l_theory)
    """    
    l_grid     = params_aps["l"]      
    hubble_z   = params_camb["hubble_z"]   
    Pk_interp  = params_camb["pk"]         
    comov_dist = params_camb["comovel_dist_z"]
    if kind_integ.lower() in ['quad', 'accurated', 'accuracy']:
        W_f1 =  kernel_hi_vec(
                #zout=zv, 
                camb_params = params_camb["camb_params"], 
                camb_results= params_camb["camb_results"],
                zrange      = params_aps['field_1_params']['z'], 
                config_dict = params_aps["field_1_params"], 
                replace_inf2nan= True, 
                use_input      = False,
                params_interp  = None,             
                get_interp_function=1
            )  
        W_f2 =  kernel_galaxy_clustering_i_vec(
                #zout=zv,
                Np  =5000, 
                pars=params_aps["field_2_params"],
                ibin=params_aps["field_2_params"]["bin_to_use"],
                get_interp_function=1 ,       
                verbose=False, # Keep verbose off during sampling
            )        
        zv = np.linspace( zv.min(), zv.max(), 500 )
        dz = np.diff(zv)[0]
        zv = np.linspace( zv.min()+2*dz, zv.max()-2*dz, 500 )   
        #print(W_f1(zv))
        #print(W_f2(zv))        
        for ii, iil in enumerate(l_grid):
            dCdz = lambda x: (hubble_z(x)/C_LIGHT)*\
                   Pk_interp(x,(iil+0.5)/comov_dist(x))*\
                   (W_f1(x)*W_f2(x))/(comov_dist(x)**2)
            cz   = integrate.quad(dCdz, zv.min(), zv.max(), 
                                  limit=quad_limit,
                                 # points=[zv.min(),zv.max()]
                                 )[0]   
            if not ii: Cz = dcopy(cz)
            else:      Cz = np.hstack(( Cz, dcopy(cz) ))
            del cz, dCdz        
        return l_grid, Cz
        
    elif kind_integ.lower() in ['simpson','fastest','fast']:
        zv = np.linspace( zv.min(), zv.max(), 500 )
        dz = np.diff(zv)[0]
        zv = np.linspace( zv.min()+2*dz, zv.max()-2*dz, 500 )
        W_f1 =  kernel_hi_vec(
                camb_params=None,#params_camb["camb_params"], 
                camb_results=None,#params_camb["camb_results"],
                zout=zv,
                zrange=params_aps['field_1_params']['z'], 
                replace_inf2nan=True, 
                use_input=True,
                config_dict=params_aps["field_1_params"], 
                params_interp=params_camb, 
                get_interp_function=0
            )
        W_f2 =  kernel_galaxy_clustering_i_vec(
                zout=zv,
                Np  =5000, 
                pars=params_aps["field_2_params"],
                ibin=params_aps["field_2_params"]["bin_to_use"],
                get_interp_function=0,        
                verbose=False, 
            )        
        Cl = np.array([])
        for ell in l_grid:
            chi_zv   = comov_dist(zv)
            k_values = (ell+0.5)/chi_zv 
            pk_values_on_grid = Pk_interp(zv, k_values, grid=False)
            dCdz  = ( (hubble_z(zv)/C_LIGHT)*pk_values_on_grid/(chi_zv**2) )
            dCdz *= W_f1*W_f2
            Cl = np.hstack(( Cl, integrate.simpson(y=dCdz, x=zv) ))            
        return l_grid, np.asarray(Cl)
    else:
        return None
    
def get_chi_theory(zv, params_camb, params_aps, kind_integ='simpson', quad_limit=5000): 
    l_grid     = params_aps["l"]      
    hubble_z   = params_camb["hubble_z"]   
    Pk_interp  = params_camb["pk"]         
    comov_dist = params_camb["comovel_dist_z"]
    
    if kind_integ.lower() in ['quad', 'accurated', 'accuracy']:
        W_f1 = c4ft.kernel_hi_vec(
                #zout=zv, 
                camb_params=params_camb["camb_params"], 
                camb_results=params_camb["camb_results"],
                zrange=params_aps['field_1_params']['z'], 
                replace_inf2nan=True, 
                use_input=False,
                config_dict=params_aps["field_1_params"], 
                params_interp=None, 
                get_interp_function=1
            )       

        zv = np.linspace( zv.min(), zv.max(), 500 )
        dz = np.diff(zv)[0]
        zv = np.linspace( zv.min()+2*dz, zv.max()-2*dz, 500 )    
        for ii, iil in enumerate(l_grid):
            dCdz = lambda x: (hubble_z(x)/C_LIGHT)*\
                   Pk_interp(x,(iil+0.5)/comov_dist(x))*\
                   (W_f1(x)*W_f1(x))/(comov_dist(x)**2)
            cz   = integrate.quad(dCdz, zv.min(), zv.max(), limit=quad_limit)[0]   
            if not ii: Cz = dcopy(cz)
            else:      Cz = np.hstack(( Cz, dcopy(cz) ))
            del cz, dCdz        
        return l_grid, Cz
        
    elif kind_integ.lower() in ['simpson','fastest','fast']:
        zv = np.linspace( zv.min(), zv.max(), 501 )
        dz = np.diff(zv)[0]
        zv = np.linspace( zv.min()+2*dz, zv.max()-2*dz, 501 )
        W_f1 = c4ft.kernel_hi_vec(
                camb_params=params_camb["camb_params"], 
                camb_results=params_camb["camb_results"],
                zout=zv,
                zrange=params_aps['field_1_params']['z'], 
                replace_inf2nan=True, 
                use_input=True,
                config_dict=params_aps["field_1_params"], 
                params_interp=params_camb, 
                get_interp_function=0
            )    
        Cl = np.array([])
        for ell in l_grid:
            chi_zv   = comov_dist(zv)
            k_values = (ell+0.5)/chi_zv 
            pk_values_on_grid = Pk_interp(zv, k_values, grid=False)
            dCdz  = ( (hubble_z(zv)/C_LIGHT)*pk_values_on_grid/(chi_zv**2) )
            dCdz *= W_f1*W_f1
            Cl = np.hstack(( Cl, integrate.simpson(y=dCdz, x=zv) ))
        return l_grid, np.asarray(Cl)
    else:
        return None    
      
def get_cg_theory(zv, params_camb, params_aps, kind_integ='simpson', quad_limit=5000):
    
    """Compute theoretical gxg-Cl values using products from Cobaya's CAMB
        for a single dataset configuration.
    """
    l_grid     = params_aps["l"]      
    hubble_z   = params_camb["hubble_z"]   
    Pk_interp  = params_camb["pk"]         
    comov_dist = params_camb["comovel_dist_z"]
    ##
    if kind_integ.lower() in ['quad', 'accurated', 'accuracy']:
        W_f2 =  kernel_galaxy_clustering_i_vec(
                Np  =5000, 
                pars=params_aps["field_2_params"],
                ibin=params_aps["field_2_params"]["bin_to_use"],
                get_interp_function=1 ,       
                verbose=False,
            )        

        zrange = np.linspace(params_aps["field_2_params"]['z_min'], 
                             params_aps["field_2_params"]['z_max'], 5000)
        dNdz_vec = np.array([ cxft.nz_spec(pars=params_aps["field_2_params"])(ix) for ix in zrange ])
        zgbins = cxft.compute_equal_number_bounds_lsst(
                                  redshift_range=zrange, 
                                  redshift_distribution=dNdz_vec, 
                                  n_bins=params_aps["field_2_params"]['nbins']
                                                      )

        zgbins = zgbins[params_aps["field_2_params"]["bin_to_use"]:params_aps["field_2_params"]["bin_to_use"]+2]
        del zrange, dNdz_vec
        zv = np.linspace( zv.min(), zv.max(), 500 )
        dz = np.diff(zv)[0]
        zv = np.linspace( zv.min()+2*dz, zv.max()-2*dz, 500 )
        Cz = np.array([])
        for ii, iil in enumerate(l_grid):
            dCdz = lambda x: (hubble_z(x)/C_LIGHT)*\
                   Pk_interp(x,(iil+0.5)/comov_dist(x))*\
                   (W_f2(x)*W_f2(x))/(comov_dist(x)**2)
            cz   = integrate.quad(dCdz, zv.min(), zv.max(),  
                                  limit=quad_limit,
                                  points=zgbins)[0]   
            Cz = np.hstack(( Cz, dcopy(cz) ))
            del cz, dCdz
        del zgbins       
        return l_grid, Cz   
    elif kind_integ.lower() in ['simpson','fastest','fast']:    
        gbin     = params_aps["field_2_params"]["bin_to_use"]
        dNdz_vec = np.array([ cxft.nz_spec(pars=params_APS["field_2_params"])(ix) for ix in zv ])
        zgbins   = cxft.compute_equal_number_bounds_lsst(redshift_range=zv, 
                                                         redshift_distribution=dNdz_vec, 
                                                         n_bins=params_APS["field_2_params"]['nbins'])
        zgbins=np.array( zgbins[gbin:gbin+2] )
        ##
        dz = np.diff(zv)[0]
        zv = np.linspace( zv.min()+2*dz, zv.max()-2*dz, 500 )
        zv1  = np.linspace( zv.min()    , zgbins.min(),  100 )
        zv12 = np.linspace( zgbins.min(), zgbins.max(), 1000 )
        zv2  = np.linspace( zgbins.max(), zv.max()    ,  100 )
        zv   = np.unique(np.concatenate((zv1, zv12, zv2))) 
        del zv1, zv2, zv12, zgbins, dNdz_vec
        W_f2 =  kernel_galaxy_clustering_i_vec(
                    zout=zv,
                    Np  =10000, 
                    pars=params_aps["field_2_params"],
                    ibin=gbin,          
                    get_interp_function=0 ,       
                    verbose=False
                                                     )        
        Cl = []
        for ell in l_grid:
            chi_zv   = comov_dist(zv)
            k_values = (ell+0.5)/chi_zv 
            pk_values_on_grid = Pk_interp(zv, k_values, grid=False)
            dCdz  = ( (hubble_z(zv)/C_LIGHT)*pk_values_on_grid/(chi_zv**2) )
            dCdz *= W_f2*W_f2
            Cl.append(integrate.simpson(dCdz, zv))        
        return l_grid, np.asarray(Cl)               

def get_camb_results(params): # pars is params_cosmo from logp
    cpars = camb.set_params(halofit_version=params['halofit_version'],
                            w=params['w0'],wa=params['wa'], dark_energy_model=params['DE_model'],
                            WantCls=0,
                            #AccuracyBoost=0.5, lSampleBoost=0.5,
                            #k_per_logint=0,
                           )#, HMCode_A_baryon=3.13, HMCode_eta_baryon=0.603, HMCode_logT_AGN=7.8)
    cpars.set_accuracy(AccuracyBoost=0.5, lSampleBoost=0.5)
    cpars.set_cosmology(
            H0=params["H0"], 
            ombh2=params["ombh2"], 
            omch2=params["omch2"], 
            mnu=params["mnu"],
            #k_per_logint=0,
                       )
    cpars.InitPower.set_params(As=params["As"], ns=params["ns"])
    cpars.set_matter_power(redshifts=params["z"],kmax=params["kmax"],silent=True)

    cres = camb.get_results(cpars)
    pk_interp = camb.get_matter_power_interpolator(
            cpars, 
            nonlinear=params["Pk_nonlinear"], 
            hubble_units=False, 
            k_hunit=False,
            kmax=params["kmax"], 
            zmax=params["z"][-1], 
        )      
                
    zstar = cres.get_derived_params()["zstar"]
    return {
            #'k':kc,'z':zc,'pkc':pkc, 
            "pk_interp"   : pk_interp, 
            "pk"          : pk_interp.P, 
            "z_camb"      : params["z"],
            "comovel_dist_z": cres.comoving_radial_distance, 
            "hubble_z"    : cres.hubble_parameter,
            "fsigma8"     : cres.get_fsigma8(),
            "sigma8"      : cres.get_sigma8(), 
            "zstar"       : zstar,
            "chistar"     : cres.comoving_radial_distance(zstar),
            "camb_params" : cpars,
            "camb_results": cres        
           }
    
def load_params_interp(filename=None):
    try:
        with np.load(filename) as data:
            k   = data['k']
            z   = data['z']
            pk  = data['pk']
            chi = data['comovel_dist']
            Hz  = data['hubble_z']
            fsigma8  = data['fsigma8']
            sigma8   = data['sigma8']
            zstar   = data['zstar']
            chistar = data['chistar']
        from scipy.interpolate import RectBivariateSpline,interp1d
        return {'pk'            : RectBivariateSpline(z, k, pk),
                'comovel_dist_z': interp1d(x=z, y=chi    , kind='cubic',fill_value='extrapolate'),
                'hubble_z'      : interp1d(x=z, y=Hz     , kind='cubic',fill_value='extrapolate'),
                'fsigma8'       : interp1d(x=z, y=fsigma8, kind='cubic',fill_value='extrapolate'),
                'sigma8'        : interp1d(x=z, y=sigma8 , kind='cubic',fill_value='extrapolate'),
                'zstar':zstar, 'chistar':chistar, 'z':z}
        
    except:
        pass
    try:
        with h5py.File(filename, "r") as data:
            k   = data["k"][:]
            pk  = data["pk"][:]
            z   = data["z"][:]
            chi = data["comovel_dist"][:]
            Hz  = data['hubble_z'][:]
            fsigma8  = data['fsigma8'][:]
            sigma8   = data['sigma8'][:]            
            zstar   = data["zstar"]#[:]
            chistar = data["chistar"]#[:]
        from scipy.interpolate import RectBivariateSpline,interp1d
        spline = UnivariateSpline(x, y, s=0, ext='extrapolate')
        return {'pk'            : RectBivariateSpline(z, k, pk),
                'comovel_dist_z': interp1d(x=z, y=chi, kind='cubic',fill_value='extrapolate'),
                'hubble_z'      : interp1d(x=z, y=Hz , kind='cubic',fill_value='extrapolate'),
                'fsigma8'       : interp1d(x=z, y=fsigma8, kind='cubic',fill_value='extrapolate'),
                'sigma8'        : interp1d(x=z, y=sigma8 , kind='cubic',fill_value='extrapolate'),                
                'zstar':zstar, 'chistar':chistar, 'z':z}
    except: 
        pass
    return {'pk'            : None,
            'comovel_dist_z': None,
            'hubble_z'      : None,
            'fsigma8'       : None, 
            'sigma8'        : None, 
            'zstar'  : None, 
            'chistar': None,
            'z'      : None}
                    
    
def load_yaml_fileformat(input_filename=None, verbose=False):
    import yaml
    try:
        with open(input_filename, 'r') as infile:
            # Use safe_load to avoid arbitrary code execution from malicious YAML files
            loaded_data = yaml.safe_load(infile)
        if None: print(f"Successfully loaded data from {input_filename}")

    except FileNotFoundError:
        if None: print(f"Error: File not found at {input_filename}")
        # Decide how to handle the error, e.g., exit or return None
        loaded_data = None
    except yaml.YAMLError as e:
        if None: print(f"An error occurred while parsing YAML: {e}")
        loaded_data = None
    except Exception as e:
        if None: print(f"An unexpected error occurred during loading: {e}")
        loaded_data = None


def convert_numpy_to_list(data):
    """Recursively converts numpy arrays and numpy number types in a dictionary to native Python types."""
    if isinstance(data, dict):
        return {k: convert_numpy_to_list(v) for k, v in data.items()}
    elif isinstance(data, list):
        return [convert_numpy_to_list(item) for item in data]
    elif isinstance(data, np.ndarray):
        # Convert numpy array to list
        return data.tolist()
    elif isinstance(data, (np.intc, np.intp, np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64)):
        # Convert numpy integer types to Python int
        return int(data)
    elif isinstance(data, (np.float16, np.float32, np.float64)):
        # Convert numpy float types to Python float
        return float(data)
    elif isinstance(data, (np.complex64, np.complex128)):
        # Convert numpy complex types to Python complex
        return complex(data)
    elif isinstance(data, (np.bool_)):
        # Convert numpy bool types to Python bool
        return bool(data)
    elif isinstance(data, (np.void)):
         # Handle numpy void types (like structured arrays) - might need custom logic
        return None # Or convert to dict/list as appropriate
    else:
        # Return item as is if it's not a numpy type we need to convert
        return data
    
def write_yaml_fileformat(input_filename=None, data2save=None, verbose=False):
    '''
    FROM DICTIONARY WITH STR,INT,FLOAT,NUMPY to YAML FILE
    '''
    import yaml
    clean_data = convert_numpy_to_list(data2save)
    with open(input_filename, 'w') as file:
        yaml.dump(clean_data, file, default_flow_style=0)
    
    
def convert_lists_to_numpy(data, verbose=0):
    """
    Converts specific lists within the loaded dictionary back to NumPy arrays.
    Modifies the dictionary in place.
    """
    if data is None:
        if verbose: print("Input data is None, cannot perform conversion.")
        return None

    # Define the paths to the keys where lists should be converted to NumPy arrays
    # Each path is a tuple of keys representing the nested structure
    numpy_keys_paths = [
        ('cosmology', 'z'),
        ('APS', 'l'),
        ('APS', 'field_1_params', 'nu'),
        ('APS', 'field_1_params', 'z'),
    ]

    for key_path in numpy_keys_paths:
        try:
            # Navigate through the dictionary using the keys in the path
            current_level = data
            for key in key_path[:-1]: # Iterate up to the second-to-last key
                if isinstance(current_level, dict) and key in current_level:
                    current_level = current_level[key]
                else:
                    # Key path doesn't fully exist
                    raise KeyError(f"Intermediate key '{key}' not found in path {key_path}")

            last_key = key_path[-1] # The final key in the path

            # Check if the final key exists and its value is a list
            if isinstance(current_level, dict) and last_key in current_level:
                if isinstance(current_level[last_key], list):
                    # Convert the list to a numpy array
                    original_list = current_level[last_key]
                    current_level[last_key] = np.array(original_list)
                    if verbose: print(f"Converted list at path {key_path} to NumPy array.")
                elif isinstance(current_level[last_key], np.ndarray):
                    # Already a numpy array (e.g., if run twice), skip.
                     if verbose: print(f"Item at path {key_path} is already a NumPy array. Skipping.")
                else:
                    # It exists but isn't a list
                    if verbose: print(f"Warning: Item at path {key_path} is not a list (Type: {type(current_level[last_key])}). Skipping conversion.")
            else:
                 if verbose: print(f"Warning: Final key '{last_key}' not found at path {key_path}. Skipping conversion.")

        except KeyError as e:
            if verbose: print(f"Warning: {e}. Skipping conversion for path {key_path}.")
        except Exception as e:
            if verbose: print(f"An unexpected error occurred during conversion for path {key_path}: {e}")

    return data

def load_yaml_fileformat(input_filename=None, verbose=False):
    import yaml
    try:
        with open(input_filename, 'r') as infile:
            # Use safe_load to avoid arbitrary code execution from malicious YAML files
            loaded_data = yaml.safe_load(infile)
        if verbose: print(f"Successfully loaded data from {input_filename}")

    except FileNotFoundError:
        if verbose: print(f"Error: File not found at {input_filename}")
        # Decide how to handle the error, e.g., exit or return None
        loaded_data = None
    except yaml.YAMLError as e:
        if verbose: print(f"An error occurred while parsing YAML: {e}")
        loaded_data = None
    except Exception as e:
        if verbose: print(f"An unexpected error occurred during loading: {e}")
        loaded_data = None


    # --- Convert specific lists back to NumPy arrays ---
    if loaded_data is not None:
        # If you want to keep the original loaded_data untouched, make a deep copy first:
        # data_to_convert = copy.deepcopy(loaded_data)
        # data_restored = convert_lists_to_numpy(data_to_convert)

        # Otherwise, modify the loaded_data dictionary in place:
        data_restored = convert_lists_to_numpy(loaded_data,verbose)

        # --- Verification (Optional) ---
        if verbose: print("\n--- Verification ---")
        if verbose: 
            if data_restored:
                try:
                    print(f"Type of data_restored['cosmology']['z']: {type(data_restored['cosmology']['z'])}")
                    # print(f"First 5 elements: {data_restored['cosmology']['z'][:5]}") # Show some data
                    print(f"Type of data_restored['APS']['l']: {type(data_restored['APS']['l'])}")
                    # print(f"First 5 elements: {data_restored['APS']['l'][:5]}") # Show some data
                    print(f"Type of data_restored['APS']['field_1_params']['nu']: {type(data_restored['APS']['field_1_params']['nu'])}")
                    # print(f"First 5 elements: {data_restored['APS']['field_1_params']['nu'][:5]}") # Show some data
                    print(f"Type of data_restored['APS']['field_1_params']['z']: {type(data_restored['APS']['field_1_params']['z'])}")
                    # print(f"First 5 elements: {data_restored['APS']['field_1_params']['z'][:5]}") # Show some data

                    # Check a non-numpy value
                    print(f"Type of data_restored['general']['verbose']: {type(data_restored['general']['verbose'])}")
                    print(f"Value of data_restored['general']['verbose']: {data_restored['general']['verbose']}")

                except KeyError as e:
                    print(f"Verification failed: Could not access expected key {e}")
                except Exception as e:
                    print(f"An error occurred during verification: {e}")

        # Now 'data_restored' holds the dictionary with NumPy arrays reconstituted
        # You can use 'data_restored' like the original dictionary variable
        return data_restored

    else:
        if verbose: print("\nSkipping conversion and verification as data loading failed.")
        return None
    
    
#################################################################################################
### GENERATING DATA & NOISE for forecasting functions ###########################################

def noise_clustering(pars_survey=None, ell_vec=None, nside=256):  
    STERADIAN = np.radians(1/60)**2 #arcmin to steradian
    if type(ell_vec)==type(None): ell_vec = np.arange(3*nside, dtype=np.int32)
    return (4*np.pi*pars_survey['f_sky'])*(STERADIAN/(pars_survey['n_gal']/pars_survey['nbins']))*np.ones_like(ell_vec)

def noise_hydrogen(sigma_pix=60 * 1e-3, ell_vec=None, f_sky=1, smoothed=False, nside=256):
    import healpy as hp
    if type(ell_vec)==type(None): ell_vec = np.arange(3*nside)
    if smoothed: 
        return f_sky*hp.nside2pixarea(nside)*(sigma_pix**2)*np.ones_like(ell_vec)
    else: 
        wn_th = np.random.normal(loc=0,scale=sigma_pix, size=12*nside*nside)
        return hp.anafast(wn_th, pol=0)*f_sky   

def cosmic_variance_ibin(clf1_=None, 
                         clf2_=None, 
                         clcx_=None,
                         nlf1_=None, 
                         nlf2_=None, 
                         nlcx_=None,
                         fwhmf1_  = None,
                         fwhmf2_  = None,
                         fwhmcx_  = None,
                         l_eff    = None, 
                         fsky     = None,                          
                         delta_ell= None,
                         b_nmt    = None,
                         ):
    
    blf1 = fnoise.bl_function(fwhmf1_/60, input_unit="degree", from_real_space=False)[:-1]
    blf2 = fnoise.bl_function(fwhmf2_/60, input_unit="degree", from_real_space=False)[:-1]    
    blcx = 1 if type(fwhmcx_)==type(None) \
             else fnoise.bl_function(fwhmcx_/60,\
                                     input_unit="degree",\
                                     from_real_space=False)[:-1]   
    nlcx_ = np.zeros_like(nlf1_) if type(nlcx_)==type(None) else nlcx_

    if type(b_nmt)==type(None):
        cx_ = clcx_ + nlcx_/(blcx**2)
        c1_ = clf1_ + nlf1_/(blf1**2)
        c2_ = clf2_ + nlf2_/(blf2**2)
    else:
        cx_ = b_nmt.bin_cell(clcx_ + nlcx_/(blcx**2))
        c1_ = b_nmt.bin_cell(clf1_ + nlf1_/(blf1**2))
        c2_ = b_nmt.bin_cell(clf2_ + nlf2_/(blf2**2))    
    err_rec = np.sqrt( cx_**2   + c1_*c2_  )
    return err_rec/np.sqrt( (2*l_eff+1)*fsky*delta_ell )   


def cosmic_variance(Clf1_  =None, 
                    Clf2_  =None, 
                    Clcx_  =None,
                    Nlf1_  =None, 
                    Nlf2_  =None, 
                    Nlcx_  =None,  
                    fwhmf1_=None, 
                    fwhmf2_=None, 
                    fwhmcx_=None,                    
                    l_eff  =None, 
                    fsky   =None, 
                    delta_ell=None,
                    bin_range=None,    
                    b_nmt    =None,                    
                    ):#arcmin
    """
    Cross-Cosmic Variance APS
    tem problemas.... prefira cosmic_variance_ibin ate estar arrumado
    """
    if type(bin_range)==list or type(bin_range)==np.ndarray:
        pass
    else:
        bin_range=np.arange(Clf1_.shape[0])
    return np.array([ cosmic_variance_ibin(clf1_     = Clf1_[ibin,:], 
                                           clf2_     = Clf2_[0,:],
                                           clcx_     = Clcx_[ibin,:],
                                           nlf1_     = Nlf1_,
                                           nlf2_     = Nlf2_,
                                           nlcx_     = Nlcx_,                                           
                                           fwhmf1_   = fwhmf1_,
                                           fwhmf2_   = fwhmf2_,
                                           fwhmcx_   = None,
                                           l_eff     = l_eff,
                                           fsky      = fsky,
                                           delta_ell = delta_ell,
                                           b_nmt     = b_nmt,
                                           ) for ibin in bin_range ])


    
def cosmic_variance_ibin_OLD(clf1_=None, 
                             clf2_=None, 
                             clcx_=None,
                             nlf1_=None, 
                             nlf2_=None,
                             l_eff=None, 
                             fsky =None, 
                             delta_ell=None,
                             b_nmt   =None,
                             fwhmf1_ =None, # in arcmin 
                             fwhmf2_ =None):
    """
    DEPRECATED
    """
    blf1 = fnoise.bl_function(fwhmf1_/60, input_unit="degree", from_real_space=False)[:-1]
    blf2 = fnoise.bl_function(fwhmf2_/60, input_unit="degree", from_real_space=False)[:-1]    
    
    cx_ = b_nmt.bin_cell(clcx_)
    c1_ = b_nmt.bin_cell(clf1_ + nlf1_/(blf1**2))
    c2_ = b_nmt.bin_cell(clf2_ + nlf2_/(blf2**2))    
    err_rec = np.sqrt( cx_**2   + c1_*c2_  )
    return err_rec/np.sqrt( (2*l_eff+1)*fsky*delta_ell )   


def cosmic_variance_OLD(Clf1_=None, 
                        Clf2_=None, 
                        Clcx_=None,
                        Nlf1_=None, 
                        Nlf2_=None, 
                        bin_range=None,
                        l_eff=None, 
                        fsky =None, 
                        delta_ell=None,
                        b_nmt=None,
                        fwhmf1_ =None, 
                        fwhmf2_ =None):
    """
    DEPRECATED
    """
    if type(bin_range)==list or type(bin_range)==np.ndarray:
        pass
    else:
        bin_range=np.arange(Clf1_.shape[0])
    return np.array([ cosmic_variance_ibin(clf1_=Clf1_[ibin,:], clf2_=Clf2_[0,:], clcx_=Clcx_[ibin,:],
                                           nlf1_=Nlf1_, nlf2_=Nlf2_,
                                           l_eff=l_eff, fsky =fsky, delta_ell=delta_ell,
                                           b_nmt=b_nmt,
                                           fwhmf1_ =fwhmf1_,fwhmf2_ =fwhmf2_) for ibin in bin_range ])





def logLKL(cl_data=None,
           errcl_data=None,
           params=None,
           params_aps=None,
           params_camb=None, 
           type_data=None,
           **kwargs):
    """
        Return **log P(data | parameters)** for the sampler.
        
        type_data   :str :     'cx': gxHI-Cl; 'ch': HIxHI-Cl; 'cg': gxg-Cl. If the value is not right, it will be 'cx'.
        params      :dict:     cosmological parameters and z-functions
        params_aps  :dict:     angular Power Spectrum parameters
        theta       :(n)array: params to be sampled 
    """
    if type_data.lower() in ['ch', 'chi', 'hixhi','hi']:
        n_   = params_aps['num_ch_data']
        zf1  = params_aps["field_1_params"]["z"]
        try:
            cl_eff = np.array([])
            for j,_c_ in enumerate(cl_data.reshape(-1,n_).T):        
                hi_index = params_aps['field_1_params']['bins_to_use']['f1f1']['bin'][j]
                params_aps['field_1_params']['bin_to_use'] = hi_index
                zmin = zf1[hi_index : hi_index + 2].min()
                zmax = zf1[hi_index : hi_index + 2].max()
                zv   = np.logspace(np.log10(zmin), np.log10(zmax), NUM_Z)                
                _, cl_th = get_chi_theory(zv, params_camb, params_aps)
                cl_eff = np.hstack(( cl_eff, params_aps["namaster"]["lbinner"](cl_th)[WSEL] ))
                del cl_th,zv,hi_index,zmin,zmax
            chi2 = ( (cl_data-cl_eff)/errcl_data )**2
            chi2 = np.sum(chi2)
            return -0.5 * chi2
        except Exception as exc:
            print("[CHLikelihood] Calculation failed. Error: %s", exc)
            return -np.inf   
    elif type_data.lower() in ['cx', 'cross', 'gxhi','hixg']:
        n_   = params_aps['num_cx_data']
        zf1  = params_aps["field_1_params"]["z"]
        try:
            cl_eff = np.array([])
            for j,_c_ in enumerate(cl_data.reshape(-1,n_).T):        
                hi_index = params_aps['field_1_params']['bins_to_use']['f1f2']['bin'][j]
                g_index  = params_aps['field_2_params']['bins_to_use']['f1f2']['bin'][j]
                params_aps['field_1_params']['bin_to_use'] = hi_index
                params_aps['field_2_params']['bin_to_use'] = g_index                
                zmin = max(zf1[hi_index : hi_index + 2].min(), params_aps["field_2_params"]["z_min"])
                zmax = min(zf1[hi_index : hi_index + 2].max(), params_aps["field_2_params"]["z_max"])  
                zv   = np.logspace(np.log10(zmin), np.log10(zmax), NUM_Z)                
                _, cl_th = get_cx_theory(zv, params_camb, params_aps,  kind_integ='quad', quad_limit=5000)
                cl_eff = np.hstack(( cl_eff, params_aps["namaster"]["lbinner"](cl_th)[WSEL] ))
                del cl_th,zv,hi_index,zmin,zmax
            chi2 = ( (cl_data-cl_eff)/errcl_data )**2
            chi2 = np.sum(chi2)
            return -0.5 * chi2
        except Exception as exc:
            print("[CXLikelihood] Calculation failed. Error: %s", exc)
            return -np.inf           
    elif type_data.lower() in ['cg', 'cgal','gxg','gal','galaxy','galaxies']:
        n_   = params_aps['num_cg_data']
        try:
            cl_eff = np.array([])
            for j,_c_ in enumerate(cl_data.reshape(-1,n_).T):        
                g_index  = params_aps['field_2_params']['bins_to_use']['f2f2']['bin'][j]
                params_aps['field_2_params']['bin_to_use'] = g_index
                zmin = params_aps["field_2_params"]["z_min"]
                zmax = params_aps["field_2_params"]["z_max"]
                zv   = np.logspace(np.log10(zmin), np.log10(zmax), NUM_Z)                
                _, cl_th = get_cg_theory(zv, params_camb, params_aps,  kind_integ='quad', quad_limit=5000)
                cl_eff = np.hstack(( cl_eff, params_aps["namaster"]["lbinner"](cl_th)[WSEL] ))
                del cl_th,zv,g_index,zmin,zmax
            chi2 = ( (cl_data-cl_eff)/errcl_data )**2
            chi2 = np.sum(chi2)
            return -0.5 * chi2
        except Exception as exc:
            print("[CGLikelihood] Calculation failed. Error: %s", exc)
            return -np.inf   
    else:
        return -np.inf 