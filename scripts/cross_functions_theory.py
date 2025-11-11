import time
from copy import deepcopy as dcopy
import numpy as np
from copy import deepcopy as dcopy
c_light = 299792458/1e3 #km/s
nu_HI   = 1420.405751768 #MHz
STERADIAN = np.radians(1/60)**2 #arcmin to steradian
#####################################################################################################
#####################################################################################################
#### HI FUNCTIONS
####
def THI_factor_constants():
    """
    ----------------------------------------------    
    PHYSICAL/COSMOLOGICAL CTEs to provide the 
    HI-APS factor
    ----------------------------------------------
    """      
    G_chbar  = 6.70883e-39 #(GeV/c2)-2
    G_kg     = 6.67430e-11 #N m2/kg2 #m3/kg/s2
    hP_kg    = 6.62607015e-34 #m2 kg s
    kB_kg    = 1.380649e-23 #m2 kg s-2 K-1
    c        = 299792458 #m s-1
    mP       = 1.67262192e-24 #g
    mHI_u    = 1.00784 #u
    u_g      = 1.66054e-24 #g
    A10      = 2.85e-15 #s-1
    H0       = 67.4 #km/s/Mpc.
    Mpc_m    = 3.0857e22 # m
    factor   = 9/(256*np.pi**2)
    nu_MHz   = 1420.405751768 #MHz
    N        = 1 #kg m/s2

    kB     = kB_kg*1e3
    hP     = hP_kg*1e3
    G      = G_kg/1e3
    mHI    = mHI_u*u_g
    H0_ms  = H0*(1e3/Mpc_m)
    nu     = nu_MHz*1e6
    nu2    = nu**2
    c3     = c**3
    h      = H0/100. #admensional 
    fac   = factor*(hP*c3*A10*H0_ms)/(G*kB*nu2*mHI)
    fac_h = factor*(hP*c3*A10)/(G*kB*nu2*mHI)*(100*(1e3/Mpc_m))
    return fac_h

def hi_brightness_temperature(camb_params=None, camb_results=None, z=None, model = None, fit_OHI=False, OHI=None, use_input=False):
    try: #if camb_params and camb_results are directly from camb
        fact = (1+z)**2/(camb_results.hubble_parameter(z)/camb_params.H0)
        if not fit_OHI or type(OHI)==type(None):     
            model = None
            oHI   = omegaHI_biasHI(model=model,z=z)['omegaHI']*camb_params.h
        else: 
            oHI = OHI**camb_params.h
        fact = fact*188.8#(THI_factor_constants()/1e-3)
        return fact*oHI #mK        
    except: #camb_params and camb_results are provide as values and function
        fact = (1+z)**2/(camb_results['hubble_z'](z)/camb_params['H0'])
        if not fit_OHI or type(OHI)==type(None):     
            model = None
            oHI   = omegaHI_biasHI(model=model,z=z)['omegaHI']*(camb_params['H0']/100)
        else: 
            oHI = OHI*(camb_params['H0']/100)
        fact = fact*188.8#(THI_factor_constants()/1e-3)
        return fact*oHI #mK

def omegaHI_biasHI(model=None,z=None):
    """
    ----------------------------------------------    
    DEPRECATED
    ----------------------------------------------
    """    
    if model=="jiajun": #valido para z < 1.#Zhang et al 2022
        return {"omegaHI":{"HOD":2.7e-4 + 1e-4*z - 8e-5*z**2,
                           "HAM":2.5e-4 - 4e-5*z - 7e-5*z**2}
               }
    elif model=="padmanabhan":#padmanabhan 2015
        z = np.array(  [0.000, 0.250, 0.500, 0.750, 1.000, 1.250, 1.500, 1.750, 2.000, 2.250, 2.500, 2.750, 3.000, 3.250, 3.400])
        oHI = np.array([3.344, 3.443, 4.523, 4.648, 4.710, 4.804, 4.766, 4.804, 4.936, 5.008, 4.750, 5.471, 5.541, 5.756, 5.971])
        bHI = np.array([0.703, 0.972, 1.026, 0.935, 1.005, 1.005, 1.049, 1.099, 1.101, 1.160, 1.261, 1.409, 1.329, 1.498, 1.802])        
        return {"z":z, "omegaHI":oHI, "bias":bHI}
    elif model=="cunnington":#cunnington2019
        return {"omegaHI":0.00048 + 0.00039*z - 0.000065*z**2}
    elif model=="irfan":#Irfan2021 baseado em bull2015 #o mesmo que IRFAN2021
        omegaHI_0 = 4.86e-4 #Irfan2021
        biasHI_0 = 0.677105 #Irfan2021
        return {"omegaHI": (omegaHI_0/4.86)*(4.8304 + 3.8856*z - 0.65119*z**2) ,
                "biasHI" : (biasHI_0/0.677105)*(0.66655 + 0.17765*z + 0.050223*z**2)}
    elif model=="crighton":#carucci2020
        omegaHI_0 = 4e-4 #Irfan2021
        biasHI_0  = 0.3 #Irfan2021        
        return {"omegaHI": omegaHI_0*(1+z)**(0.6) ,
                "biasHI" : biasHI_0*(1+z)+0.6}
#    if type(model)==type(None):
#        return {"omegaHI": 4.86e-4,
#                "biasHI" : 1}          
    else:
        return {"omegaHI": 4.86e-4,
                "biasHI" : 1}       
    
def omegaHI_biasHI_test(model=None,z=None, OHI=None, bHI=None):
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
#####################################################################################################
#####################################################################################################
#### KERNELS
####
def window_1(z=None, zmin=None, zmax=None):
    if (z>=zmin)*(z<=zmax): 
        return 1/(zmax-zmin)
    else: 
        return 0

def window_2(value=None, arr=None):
    """
    arr: array containing the binned redshift
    value: a specific redshift to be evaluated
    """
    if (value>=arr.min())*(value<=arr.max()):
        vleft =np.sort(arr[np.where(arr<=value)[0]])[-1]#first number to left
        vright=np.sort(arr[np.where(arr>=value)[0]])[0 ]#first number to right
        return window_1(z=value, zmin=np.minimum(vleft,vright), zmax=np.maximum(vleft,vright))
    else:
        return 0
    
def chi_vec(camb_results=None): #radial comovel distance
    from scipy import interpolate
    zv         = np.linspace(0,1200)
    rv         = camb_results.comoving_radial_distance(zv)
    z_chi      = interpolate.interp1d(rv, zv, kind='linear') #function
    z_star     = camb_results.get_derived_params()['zstar']
    chi_star   = camb_results.comoving_radial_distance(z_star)
    #zs         = results.redshift_at_comoving_radial_distance(np.linspace(0,chi_star,nz))
    z_chi_camb = lambda y: results.redshift_at_comoving_radial_distance(y)    
    return {'zvec':zv, 'chivec':rv, 'z_chi': z_chi, 'z_chi_camb':z_chi_camb, 'z_star':z_star, 'chi_star':chi_star}

def kernel_hi_dict(z_survey=None, camb_params=None, camb_results=None,config_dict=None, replace_inf2nan=False):
    #2set higher z-range, with 0
    Khi_all= np.array([])
    Z_all  = np.array([])
    for i in np.arange(z_survey.size-1):
        zmin = z_survey[i]
        zmax = z_survey[i+1]
        z_range = np.linspace(zmin,zmax,100)
        for j,zj in enumerate(z_range):
            if not j: 
                khi = kernel_hi(camb_params=camb_params, camb_results=camb_results, z=zj, zrange=z_survey, config_dict=config_dict, replace_inf2nan=replace_inf2nan)
                #khi = kernel_hi_vec(camb_params=camb_params, camb_results=camb_results, z=zj, zmin = zmin, zmax=zmax)
            else:
                khi = np.hstack(( khi,kernel_hi(camb_params=camb_params, camb_results=camb_results, z=zj, zrange=z_survey, config_dict=config_dict, replace_inf2nan=replace_inf2nan) ))
                #khi = np.hstack((khi,kernel_hi_vec(camb_params=camb_params, camb_results=camb_results, z=zj, zmin = zmin, zmax=zmax)))
        Khi_all = np.hstack((Khi_all,khi))
        Z_all   = np.hstack((Z_all,z_range))
        if not i:
            Khi_dict = {str(i+1):khi}
            Khi_vec  = khi
            Z_dict   = {str(i+1):z_range}
            Z_vec    = z_range
        else:
            Khi_dict[str(i+1)] = khi 
            Khi_vec = np.vstack(( Khi_vec,khi ))
            Z_dict[str(i+1)]   = z_range 
            Z_vec  = np.vstack(( Z_vec,z_range ))    
    #dNdz_intp = interp1d(dNdz_vec['z'], dNdz_vec['dNdz'][j,:], kind='cubic')
    return {'kernel_vec' :Khi_vec , 'z_vec' :Z_vec,
            'kernel_dict':Khi_dict, 'z_dict':Z_dict,
            'kernel_all' :Khi_all , 'z_all' :Z_all}


def kernel_hi(camb_params=None, camb_results=None, z=None, zrange=None, config_dict=None, replace_inf2nan=False):
    #model=None, fit_bHI=False, bHI=None, fit_OHI=False, OHI=None):
    #zrange: array of z bins of the survey (not the whole z analyzed)
    #z: value to be analyzed (in/out zrange)
    omegaHI_model = config_dict['omegaHI_model']
    biasHI_model  = config_dict['biasHI_model']
    fit_bHI = config_dict['fit_biasHI']
    fit_OHI = config_dict['fit_omegaHI']
    bHI = config_dict['biasHI']
    OHI = config_dict['omegaHI']
    #if you choice to determine (fit) HI bias (bHI) and Omega HI (OHI), 
    #you must to specify True for the fit_bHI (fit_OHI) parameter and to set the parameter value
    #otherwise, the code will assume the model provided.
    w1  = window_2(value=z, arr=zrange)#selection function
    if not w1: return 0
    THI = hi_brightness_temperature(camb_params=camb_params, camb_results=camb_results, z=z, model=omegaHI_model, fit_OHI=fit_OHI, OHI=OHI)
    if not fit_bHI or type(bHI)==type(None): 
        bHI = omegaHI_biasHI(model=biasHI_model,z=z)['biasHI']
    else:
        pass
    if replace_inf2nan and (bHI*THI*w1 == np.inf): return np.nan
    else: return bHI*THI*w1 

def kernel_cmb_z(z_=None, camb_params=None, camb_results=None,verbose=False):
    if verbose: timei = time.time()
    #wl_pars = survey_parameters(survey_name =survey_name, data_release=data_release)
    #zmin    = 0
    #zmax    = 20#camb_results.get_derived_params()['zstar']  
    #zstar   = camb_results.get_derived_params()['zstar']  
    #chi_zprime = lambda x: results.comoving_radial_distance(x)
    fact     = 1.5*camb_params.omegam*camb_params.H0**2/c_light/camb_results.hubble_parameter(z_)
    chi_z    = camb_results.comoving_radial_distance(z_)
    chi_star = chi_vec(camb_results)['chi_star']
    if verbose: print('Total processing time: {0:.2f} seg'.format(time.time()-timei))          
    return fact*(1+z_)*chi_z*(chi_star-chi_z)/chi_star  

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

def kernel_galaxy_lensing_i_dict(camb_params=None, camb_results=None, pars=None,  binning=True, ibin=False, verbose=False, Np=500):#,cut_to_bin=False, ):
    from scipy.interpolate import interp1d
    from scipy import integrate
    if verbose: timei = time.time()    
    #dNdz_vec = dNdz_func(zbins=zbins, pars=lens_paparsrs, verbose=False)
    dNdz_vec = dNdz_func(pars=pars , binning=binning, verbose=False)
    #
    chi_z = lambda x: camb_results.comoving_radial_distance(x)    
    fact  = lambda y: 1.5*camb_params.omegam*camb_params.H0**2/c_light/camb_results.hubble_parameter(y)*(1+y)*chi_z(y)
    ####
    #for j in range(dNdz_vec['dNdz'].shape[0]):
    dNdz_intp = interp1d(dNdz_vec['z'], dNdz_vec['dNdz'][ibin,:], kind='cubic')
    ####
    zi_range = np.linspace(pars['z_min'], pars['z_max'],Np)
    gx = np.array([ fact(zi)*integrate.quad(lambda x: dNdz_intp(x)*( (chi_z(x)-chi_z(zi))/chi_z(x) ), zi, pars['z_max'])[0] for zi in zi_range ])
    gx_intp = interp1d(zi_range, gx, kind='cubic')#, fill_value="extrapolate")
    ####
    if verbose: print('Total processing time: {0:.2f} seg'.format(time.time()-timei))             
    return {'W':gx, 'z':zi_range, 'W_interp':gx_intp, 'bin':ibin}

def kernel_galaxy_lensing_dict(camb_params=None, camb_results=None, pars=None, binning=True, cut_to_bin=False,Np=500, verbose=False):
    c_light = 299792458/1e3
    from scipy.interpolate import interp1d
    from scipy import integrate
    if verbose: timei = time.time()        
    #dNdz_vec = dNdz_func(pars=pars , binning=binning, verbose=False)
    dNdz_vec = dNdz_func(pars=pars, binning=binning,  verbose=False, type_field=pars['field'])
    #
    chi_z = lambda x: camb_results.comoving_radial_distance(x)    
    fact  = lambda y: 1.5*camb_params.omegam*camb_params.H0**2/c_light/camb_results.hubble_parameter(y)*(1+y)*chi_z(y)
    ####
    for j in range(dNdz_vec['dNdz'].shape[0]):
        dNdz_intp = interp1d(dNdz_vec['z'], dNdz_vec['dNdz'][j,:], kind='cubic')
        ####
        zi_range = np.linspace(pars['z_min'], pars['z_max'],Np)
        #if cut_to_bin: zi_range = np.linspace(zbins[j],zbins[j+1],300)
        ####
        igx = np.array([ fact(zi)*integrate.quad(lambda x: dNdz_intp(x)*( (chi_z(x)-chi_z(zi))/chi_z(x) ), zi, pars['z_max'])[0] for zi in zi_range ])
        if not j:
            gx = dcopy(igx)
            zv_range = dcopy(zi_range)
        else:
            gx = np.vstack(( gx,dcopy(igx) ))
            zv_range = np.vstack(( zv_range, dcopy(zi_range) ))
    if verbose: print('Total processing time: {0:.2f} seg'.format(time.time()-timei))             
    return {'W':gx, 'z':zv_range}  

def kernel_cmb_interp(camb_params=None, camb_results=None, interp=False):
    from scipy.interpolate import interp1d
    if not interp: return lambda zi:  kernel_cmb_z(z_=zi, camb_params=camb_params, camb_results=camb_results)
    else:
        zp = np.linspace(0,1150,1000000)
        kcmb = np.array([ kernel_cmb_z(z_=zi, camb_params=camb_params, camb_results=camb_results) for zi in zp ])
        return interp1d(zp, kcmb, kind='cubic')
        
        

def angular_power_spectrum_ij_old(camb_params=None, 
                              camb_results=None, 
                              Pk_camb_interp=None,
                              W_f1_i=None, 
                              W_f2_j=None, 
                              l=None, 
                              z_min=None, 
                              z_max=None, 
                              verbose=False):
    from scipy import integrate
    #''''
    #This function is deprecated and will be removed soon
    #''''
    #       W_f1_i: kernel of field 1 for i-bin
    #       W_f2_j: kernel of field 2 for j-bin    
    #[z_min,z_max]: integration limits
    zvec  = np.linspace( z_min, z_max, 500 )
    dz    = np.diff(zvec)[0]
    zvec  = np.linspace( zvec.min()+2*dz, zvec.max()-2*dz, 500 )        
    z_min = zvec.min()
    z_max = zvec.max()          
    timei=time.time()
    for ii, iil in enumerate(l):
        dCdz = lambda x: (camb_results.hubble_parameter(x)/c_light)*\
               Pk_camb_interp.P(x,(iil+0.5)/camb_results.comoving_radial_distance(x))*\
               (W_f1_i(x)*W_f2_j(x))/(camb_results.comoving_radial_distance(x)**2)
        cz   = integrate.quad(dCdz, z_min, z_max)[0]   
        if not ii: Cz = dcopy(cz)
        else:      Cz = np.hstack(( Cz, dcopy(cz) ))
        del cz, dCdz
    if verbose: print('Total processing time: {0:.2f} seg'.format(time.time()-timei))   
    return {'l':l,'cl':Cz}


    
def angular_power_spectrum_ij(camb_params =None, 
                              camb_results=None, 
                              Pk_camb_interp=None,
                              params_aps=None,
                              W_f1_i=None, 
                              W_f2_j=None, 
                              l=None, 
                              z_min=None, 
                              z_max=None, 
                              kind_integrator='quad',
                              quad_limit=100,
                              quad_bound=None,
                              verbose=False):
    '''
    Args:
       camb_params    :dict      :
       camb_results   :dict      :
       Pk_camb_interp :intp      : interpolation of camb matter PK
       W_f1_i         :intp/array: kernel of field 1 for i-bin
       W_f2_j         :intp/array: kernel of field 2 for j-bin
       l              :array     : multipole range
       z_min          :int/float : lower z bound
       z_max          :int/float : upper z bound
       kind_integrator:str       : sort of integration (quad or simpson)
       quad_limit     :int       : if you select quad, it will set the number of sub intervals       
       quad_bound     :2-array   : A sequence of break points in the bounded integration interval
                                   where local difficulties of the integrand may occur
       verbose        :bool      : 
    Returns:
           Dictionary: 
               l  :array: multipole
               cl :array: angular power spectrum
    Ps.: 
    '''
    import scipy, numpy
    from scipy import integrate
    from scipy.interpolate import interp1d

    if kind_integrator.lower() in ['quad', 'accurated', 'accuracy']:
        '''
          W_f1(2)_i(j) must be a scipy.interpolate._interpolate.interp1d 
          The 'quad' integration 
        '''              
        zvec  = np.linspace( z_min, z_max, 500 )
        dz    = np.diff(zvec)[0]
        zvec  = np.linspace( zvec.min()+2*dz, zvec.max()-2*dz, 500 )        
        z_min = zvec.min(); z_max = zvec.max()          
        if type(W_f1_i) != scipy.interpolate._interpolate.interp1d:
            W_f1_i = interp1d(x=zvec, y=W_f1_i, kind='cubic', fill_value=np.nan)
        if type(W_f2_j) != scipy.interpolate._interpolate.interp1d:
            W_f2_j = interp1d(x=zvec, y=W_f2_j, kind='cubic', fill_value=np.nan)            
        timei = time.time()
        for ii, iil in enumerate(l):
            dCdz = lambda x: (camb_results.hubble_parameter(x)/c_light)*\
                   Pk_camb_interp.P(x,(iil+0.5)/camb_results.comoving_radial_distance(x))*\
                   (W_f1_i(x)*W_f2_j(x))/(camb_results.comoving_radial_distance(x)**2)
            cz   = integrate.quad(dCdz, z_min, z_max, 
                                  limit=quad_limit, 
                                  points=quad_bound)[0]   
            if not ii: Cz = dcopy(cz)
            else:      Cz = np.hstack(( Cz, dcopy(cz) ))
            del cz, dCdz
        if verbose: print('Total processing time: {0:.2f} seg'.format(time.time()-timei))   
        return {'l':l,'cl':Cz}
        
    elif kind_integrator.lower() in ['simpson','fastest','fast']:
        return None
    else:
        return None        
##############################
def kernel_galaxy_clustering_i_dict(camb_params       = None  , 
                                    camb_results      = None  , 
                                    pars              = None  , 
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
                                    verbose           = False  ):
    #c_light = 299792458/1e3
    #EPS=1e-50
    from scipy.interpolate import interp1d
    if verbose: timei = time.time()     
    #
    delz = np.abs(pars['z_max']-pars['z_min'])/Np
    z_range   = np.linspace(pars['z_min']-delz, pars['z_max']+delz,Np+2)
    bias_gc   = get_bias_model(bias_model=bias_model, 
                               bias_value=bias_value, 
                                     pars=pars, 
                             path_to_file=path_to_bias_file, 
                                  z_range=z_range)
    dNdz_vec  = dNdz_func(pars=pars, 
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
        MB = magnification_bias(camb_params=camb_params, camb_results=camb_results, pars=pars,
                                slope_s_model=mb_model, z_range=z_range,
                                binning=binning, ibin=ibin, Np=Np, verbose=False)
        kernel_gc+=MB
    ######################################################################        
    kernel_gc_intp = interp1d(z_range[1:-1], kernel_gc, kind='cubic')#, fill_value='extrapolate')
    if verbose: print('Total processing time: {0:.2f} seg'.format(time.time()-timei))             
    return {'W':kernel_gc, 'z':z_range[1:-1], 'W_interp':kernel_gc_intp, 'bin':ibin}    


def kernel_galaxy_clustering_dict(camb_params=None  , camb_results=None, pars=None, cut_to_bin=False, #binning=True,
                                  include_rsd=False , include_mb=False , mu_rsd=0 , growth_rate_f=None, 
                                  bias_model ='SQRT', bias_value=1     , path_to_bias_file=None,
                                  mb_model   =False , mb_alpha_value=0 , path_to_mb_file=None  ,Np=500,
                                  verbose=False):
    '''
    RSD:Redshift Spatial Distortion
    MB:Magnification Bias
    mu_rsd: <k,n>=cos(theta)
    growth_rate_f: add in the future the possibility to read a file [z,f(z)]
    pars: parameters of Galaxy survey
    bias_model: constant, LSST, EUCLID, READ [to read a file]
    '''
    c_light = 299792458/1e3
    binning = True
    from scipy.interpolate import interp1d
    if verbose: timei = time.time()     
    #
    z_range  = np.linspace(pars['z_min'], pars['z_max'], Np)
    bias_gc  = get_bias_model(bias_model=bias_model, bias_value=bias_value, pars=pars, path_to_file=path_to_bias_file, z_range=z_range)
    dNdz_vec = dNdz_func(pars=pars, binning=binning,  verbose=False, type_field=pars['field'])
    #chi_z = lambda x: camb_results.comoving_radial_distance(x)    
    #fact  = lambda y: camb_results.hubble_parameter(y)/c_light
    for j in range(dNdz_vec['dNdz'].shape[0]):
        dNdz_intp = interp1d(dNdz_vec['z'], dNdz_vec['dNdz'][j,:], kind='cubic')
        #igc       = np.array([ dNdz_intp(zi)*camb_results.hubble_parameter(zi)/c_light for zi in z_range ])
        igc       = np.array([ dNdz_intp(zi)for zi in z_range ])
        igc*=bias_gc
        if not j:
            kernel_gc= dcopy(  igc   )
            zv_range = dcopy(z_range )
        else:
            kernel_gc= np.vstack(( kernel_gc, dcopy(igc)     ))
            zv_range = np.vstack((  zv_range, dcopy(z_range) ))        
    if include_rsd:
        f_vec = interp1d(camb_results.z_input, camb_results.get_fsigma8()/camb_results.get_sigma8(), kind='cubic')(z_range)
        kernel_gc+=(f_vec*mu_rsd**2)*kernel_gc
    if include_mb:
        MB = magnification_bias(camb_params=camb_params, camb_results=camb_results, pars=pars,
                                slope_s_model=mb_model, z_range=z_range,
                                binning=binning, ibin=None, Np=Np, verbose=False)
        kernel_gc+=MB
        pass    
    if verbose: print('Total processing time: {0:.2f} seg'.format(time.time()-timei))             
    return {'W':kernel_gc, 'z':zv_range}  
#####################################################################################################
#####################################################################################################
#### take a matrix and return a dictionary with each z-bin/nu-channel in a key-dictionary
####
def cls_from_matrix2dict(cl_matrix_=None, binning_name_start_as_1=True): #format: [cl_matrix_] = (nu+1), nl. This `1` is due to the l array
    nbins_,_ = cl_matrix_.shape
    nbins_   = nbins_ - 1
    cl_dict_ = {'l': cl_matrix_[0].astype(np.int16)}
    if binning_name_start_as_1:
        for i in range(nbins_):
            cl_dict_["bin{}".format(str(i+1))] = cl_matrix_[i+1]
    else:
        for i in range(nbins_):
            cl_dict_["bin{}".format(str(i))] = cl_matrix_[i+1]            
    return cl_dict_

def get_galaxy_shotnoise(pars=None, lmin=0,lmax=int(3*256)):
    '''
     Calculate the galaxy shot-noise APS from a specific 

     pars (dict): contains galaxy survey information. Mainly, minimum/maximum survey redshift
    '''
    from scipy import integrate
    from scipy.interpolate import CubicSpline    
    zmin, zmax, Np = pars['z_min'], pars['z_max'], 1000
    zrange = np.linspace(zmin, zmax, Np)
    if pars['survey'].upper() in ['LSST']:
        dNdz_vec = np.array([ nz_spec(pars=pars)(ix) for ix in zrange ])
        bins = compute_equal_number_bounds_lsst(zrange, dNdz_vec, pars["nbins"])
        ngbar=[]
        for zimin,zimax in zip(bins[:-1], bins[1:]):
            idNdz = np.array([ integrate.quad(lambda izp: nz_convolved(zp=izp,pars=pars, galaxy_normalization=True)(iz), zimin, zimax)[0] \
                              for iz in zrange ])
            idNdz_intp = CubicSpline(zrange, idNdz)    
            ngbar.append( integrate.quad(lambda iz : idNdz_intp(iz), zmin, zmax)[0] )
        ngbar = np.asarray(ngbar)
        lvec  = np.arange(lmin,lmax,1)
        clvec = np.vstack([ (4*np.pi*pars['f_sky'])*( STERADIAN/(ingbar) )*np.ones_like(lvec) for ingbar in ngbar ])
        return {'l':lvec, 'cl':clvec}
    elif pars['survey'].upper() in ['EUCLID']:
        if int(pars['zbin_edges'].size-1)==int(pars['nbins']):
            bins = pars['zbin_edges']  
        else:
            try:
                zimin=pars['zbin_edges'].min()
                zimax=pars['zbin_edges'].max()
                bins = np.around( np.linspace(zimin, zimax, pars['nbins']), decimals=6)                    
                pars['zbin_edges'] = bins
                del zimin,zimax
            except:
                bins = np.around( np.linspace(z_min, z_max, pars['nbins']), decimals=6)                    
                pars['zbin_edges'] = bins        
        ngbar=[]
        for zimin,zimax in zip(bins[:-1], bins[1:]):
            idNdz = np.array([ integrate.quad(lambda izp: nz_convolved(zp=izp,pars=pars, galaxy_normalization=True)(iz), zimin, zimax)[0] \
                              for iz in zrange ])
            idNdz_intp = CubicSpline(zrange, idNdz)    
            ngbar.append( integrate.quad(lambda iz : idNdz_intp(iz), zmin, zmax)[0] )
        ngbar = np.asarray(ngbar)
        lvec  = np.arange(lmin,lmax,1)
        clvec = np.vstack([ (4*np.pi*pars['f_sky'])*( STERADIAN/(ingbar) )*np.ones_like(lvec) for ingbar in ngbar ])
        return {'l':lvec, 'cl':clvec}                
        
    
#####################################################################################################
#####################################################################################################
#### SAVING
####
def savedata_pk(k_, Pk_, filename=None, path=None, prefix=None, suffix=None, header =''):
    import os
    nu,_ = Pk_.shape
    if not (prefix==None): filename = "_".join((  prefix  , filename  ))
    if not (suffix==None): filename = "_".join((  filename, suffix    )) 
    filename = ".".join((filename, "txt"))
    pathname = os.path.join(path,filename)
    Pk_=np.hstack(( np.asmatrix(k_).T, Pk_.T ))
    np.savetxt(pathname, Pk_, fmt=['%e']+["%e"]*nu, delimiter=" ", header=header)
    return pathname

def savedata(l_, Cl_, filename=None, path=None, prefix=None, suffix=None, header =''):
    import os
    nu,_ = Cl_.shape
    if not (prefix==None): filename = "_".join((  prefix  , filename  ))
    if not (suffix==None): filename = "_".join((  filename, suffix    )) 
    filename = ".".join((filename, "txt"))
    pathname = os.path.join(path,filename)
    Cl_=np.hstack(( np.asmatrix(l_).T, Cl_.T ))
    np.savetxt(pathname, Cl_, fmt=['%d']+["%e"]*nu, delimiter=" ", header=header)
    return pathname
    
    
def save_plots(params_=None, pathname_=None, clf_dict=None, field=None):
    import matplotlib 
    import matplotlib.pyplot as plt
    from matplotlib import cm
    ####################################
    font = {'weight' : 'bold','size'   : 22}
    matplotlib.rc('font', **font)
    plt.rc('font',   size=209)  #set defaults so that the plots are readable
    plt.rc('axes',   titlesize=20)
    plt.rc('axes',   labelsize=20)
    plt.rc('xtick',  labelsize=20)
    plt.rc('ytick',  labelsize=20)
    plt.rc('legend', fontsize =20)
    plt.rc('figure', titlesize=20)
    plt.rc('text',   usetex=True)    
    fig  = plt.figure()
    grid = plt.GridSpec(1,1,top=1.,right=1.5,wspace=0.25)
    ax   = plt.subplot(grid[0,0])
    ####################################
    if field.lower()=='hi':
        plt.loglog(clf_dict['l'],clf_dict['bin1'],linestyle='solid', linewidth=4, color='grey', label='bin1')
        for i,iname in enumerate(clf_dict.keys()):
            if iname!='l':plt.loglog(clf_dict['l'],clf_dict[iname],linestyle='solid', color=cm.Reds(i/len(clf_dict.keys())))
    elif field.lower()=='cmbwl':
        plt.loglog(clf_dict['l'],clf_dict['bin1'],linestyle='solid', linewidth=4, color='grey', label='bin1')
        for i,iname in enumerate(clf_dict.keys()):
            if iname!='l':plt.loglog(clf_dict['l'],clf_dict[iname],linestyle='solid', color=cm.Purples(i/len(clf_dict.keys())))
    elif field.lower()=='shear':
        plt.loglog(clf_dict['l'],clf_dict['bin1'],linestyle='solid', linewidth=4, color='grey', label='bin1')
        for i,iname in enumerate(clf_dict.keys()):
            if iname!='l':plt.loglog(clf_dict['l'],clf_dict[iname],linestyle='solid', color=cm.Purples(i/len(clf_dict.keys())))            
    elif field.lower()=='cross':
        plt.loglog(clf_dict['l'],clf_dict['bin1'],linestyle='solid', linewidth=4, color='grey', label='bin1')
        for i,iname in enumerate(clf_dict.keys()):
            if iname!='l':plt.loglog(clf_dict['l'],clf_dict[iname],linestyle='solid', color=cm.Blues(i/len(clf_dict.keys())))	        
    else: raise Exception('NO field like that.')		    
    plot_type = ['log','log']
    plt.ylabel(r'$C_{\ell}^{\small \textrm{HI} -\small \textrm{HI}}$')
    plt.xlabel(r'$\ell$')
    if int(params_['l_min']): plt.xlim([params_['l_min'],params_['l_max']])    
    else: plt.xlim([1,params_['l_max']]) 
    _lstr_ = np.asarray(pathname_.split('/'))
    pathname_ = '/'.join(( '/'.join( (_lstr_[:-3]) ), 'images', _lstr_[-1].replace('.txt', '.png') ))   
    #pathname_ = '.'.join(( pathname_.split('.txt')[0], 'png' ))
    plt.savefig(pathname_, dpi=100, bbox_inches='tight')
    return pathname_
    
    
def get_Lformat(L_=None):
    return 'L{zeros}{L}'.format(zeros=(4-len(str(L_)))*'0', L=L_)    
    
def savefits(data_=None, params_=None, verbose=True, save_mixmatrix=False, S_mixmatrix=None):#,pathout=None):
    import handling_data   as hdata
    import astropy.io.fits as fits
    if verbose: print('Saving data at {}'.format(params_['pathout_field']))
    nustep = (params_['frequency']['max'] - params_['frequency']['min'])/params_['frequency']['nbands']
    data0  = np.arange(params_['frequency']['min'], params_['frequency']['max'] + nustep, nustep)
    hdr0   = hdata.creating_primary_FITSheader(params_)
    hdu0   = fits.PrimaryHDU(header=hdr0, data=data0)
    hdr1   = fits.Header()
    hdr1   = hdata.creating_observed_FITSheader(params_)
    data1  = data_
    hdu1   = fits.ImageHDU(  header=hdr1, data=data1)
    if not save_mixmatrix:
        hdul = fits.HDUList([hdu0, hdu1])  
    else: 
        hdr2 = hdata.creating_mixmatrix_FITSheader() 
        hdu2 = fits.ImageHDU( header=hdr2, data=S_mixmatrix)
        hdul = fits.HDUList([ hdu0, hdu1, hdu2])  
        
    hdul.writeto(params_['pathout_field'], overwrite=True)     
    

def get_zeff(numin =980, numax =1260, nbands= 30, n_curves=3):
    import handling_data  as hdata
    nu_HI   = 1420.405751768 #MHz
    nu_vec    = hdata.nu_bins_vector(numin_ =numin, numax_ =numax, nbands_= nbands)['nu']
    nuvec_eff = np.array([ 0.5*(nu_vec[i]+nu_vec[i+1]) for i in range(nu_vec.size-1) ])
    ibins     = np.linspace(0, nbands-1, n_curves, dtype=np.int16)
    return {'zeff': np.around(nu_HI/np.array([nuvec_eff[ibins]]) -1, decimals=2)[0],
            'bins':ibins, 
            'nu_eff':nuvec_eff[ibins]}
##############################################################################################################
##############################################################################################################
def prob_euclid_vec(z_=None,zp_=None,pars=None):
    fout  = pars['f_out']#0.1
    cb    = pars['c_b']#1.0
    co    = pars['c_0']#1.0
    zb    = pars['z_b']#0.0
    zo    = pars['z_0']#0.1
    sb    = pars['sigma_b']#0.05
    so    = pars['sigma_0']#0.05
    #prob_ph1 = ((1-fout)/(np.sqrt(2*np.pi)*sigmab*(1+z_)))*(-0.5*((z_-cb*zp_-zb)/(sigmab*(1+z_)))**2)
    #prob_ph2 = ((0+fout)/(np.sqrt(2*np.pi)*sigmab*(1+z_)))*(-0.5*((z_-c0*zp_-z0)/(sigma0*(1+z_)))**2)
    Nz = 1#0.5*special.erfc(-z/(np.sqrt(2)*sb))    
    sb = sb*(1+z_)*Nz
    so = so*(1+z_)*Nz
    prob_ph1 = ((1-fout)/(np.sqrt(2)*sb))*np.exp( -0.5*( (z_-(cb*zp_+zb))/sb )**2 )
    prob_ph2 = ((0+fout)/(np.sqrt(2)*so))*np.exp( -0.5*( (z_-(co*zp_+zo))/so )**2 )     
    return prob_ph1+prob_ph2

def prob_lsst_vec(z_=None, zp_=None, pars=None):
    from scipy import special
    sigma_z = (1+z_)*pars['sigma_z']
    NZ_lsst = 0.5*special.erfc(-z_/(np.sqrt(2)*sigma_z))
    return np.exp(-0.5*((z_-zp_)/sigma_z)**2)/(NZ_lsst*np.sqrt(2*np.pi)*sigma_z)

def prob_euclid_func(zp_=None,pars=None):
    fout  = pars['f_out']#0.1
    cb    = pars['c_b']#1.0
    co    = pars['c_0']#1.0
    zb    = pars['z_b']#0.0
    zo    = pars['z_0']#0.1
    sb    = pars['sigma_b']#0.05
    so    = pars['sigma_0']#0.05
    Nz    = 1#0.5*special.erfc(-z/(np.sqrt(2)*sb))    
    #sb = lambda x: sb*(1+x)*Nz
    #so = lambda x: so*(1+x)*Nz    
    prob_ph1 = lambda x: ((1-fout)/(np.sqrt(2)*sb*(1+x)*Nz))*np.exp( -0.5*( (x-(cb*zp_+zb))/sb )**2 )
    prob_ph2 = lambda x: ((0+fout)/(np.sqrt(2)*so*(1+x)*Nz))*np.exp( -0.5*( (x-(co*zp_+zo))/so )**2 ) 
    return lambda x: prob_ph1(x) + prob_ph2(x)

def prob_lsst_func(zp_=None, pars=None):
    from scipy import special
    sigma_z = lambda x: (1+x)*pars['sigma_z']
    NZ_lsst = lambda x: 0.5*special.erfc(- x/(np.sqrt(2)*sigma_z(x)))
    return lambda x: np.exp(-0.5*((x-zp_)/sigma_z(x))**2)/(NZ_lsst(x)*np.sqrt(2*np.pi)*sigma_z(x))
################
def nz_spec(pars=None):
    #Smail distribution for faint galaxies
    if pars['survey'].upper()=='EUCLID':    
        #z0_ = pars['z_m']/np.sqrt(2)
        z0_ = pars['z_0']        
    elif pars['survey'].upper()=='LSST':   
        z0_ = pars['z_0']        
    else:
        z0_=0.1
    return lambda x: (x/z0_)**pars['beta']*np.exp(-(x/z0_)**pars['alpha'])

def nz_convolved(zp=None, pars=None, galaxy_normalization=False):
    if not galaxy_normalization:
        if pars['survey'].upper()=='LSST': 
            return lambda x: prob_lsst_func(  zp_=zp, pars=pars)(x)*nz_spec(pars=pars)(x)   
        elif pars['survey'].upper()=='EUCLID': 
            return lambda x: prob_euclid_func(zp_=zp, pars=pars)(x)*nz_spec(pars=pars)(x)       
        else:
            raise NameError
    else:
        from scipy import integrate
        zmin, zmax = pars['z_min'], pars['z_max']
        ngal = pars['n_gal']#/STERADIAN        
        A = ngal/integrate.quad(nz_spec(pars), zmin, zmax)[0]
        if pars['survey'].upper()=='LSST': 
            return lambda x: A*prob_lsst_func(  zp_=zp, pars=pars)(x)*nz_spec(pars=pars)(x)   
        elif pars['survey'].upper()=='EUCLID': 
            return lambda x: A*prob_euclid_func(zp_=zp, pars=pars)(x)*nz_spec(pars=pars)(x)       
        else:
            raise NameError        
                
################        
def dNdz_func(pars=None, binning=False, verbose=False, 
              type_field='lensing', Np=500,
              limit=2000, epsabs=1e-20, epsrel=1e-20):
    from scipy.interpolate import CubicSpline
    from scipy import integrate
    if verbose: timei = time.time()
    #zrange = np.linspace(pars['z_min'], pars['z_max'], 1000)
    z_min,z_max = pars['z_min'], pars['z_max']
    zrange = np.linspace(z_min, z_max, Np)
    if type_field.lower() in ['lensing','lens', 'shear']:
        if binning:
            pars['zbin_edges'] = np.linspace( pars['zbin_min'], pars['zbin_max'], pars['nbins']+1 )
            nch   = pars['nbins']
            zbins = pars['zbin_edges']  
            for i in range(nch):
                zi_min, zi_max = zbins[i], zbins[i+1]
                idNdz      = np.array([ integrate.quad(lambda izp: nz_convolved(zp=izp,pars=pars)(iz), zi_min, zi_max)[0] for iz in zrange ])
                idNdz_intp = CubicSpline(zrange, idNdz)    
                idNdz      = idNdz/integrate.quad(lambda iz : idNdz_intp(iz), z_min, z_max, limit=limit, epsabs=epsabs, epsrel=epsrel)[0]
                if not i: dNdz_vec = dcopy(idNdz)
                else: dNdz_vec = np.vstack(( dNdz_vec, dcopy(idNdz) )) 
            if verbose: print('Total processing time: {0:.2f} seg'.format(time.time()-timei)) 
        else:
            dNdz_vec = np.array([ nz_spec(pars=pars)(iz) for iz in zrange ])        
            dNdz_intp = CubicSpline(zrange, dNdz_vec)            
            dNdz_vec=dNdz_vec/integrate.quad( lambda iz: dNdz_intp(iz),z_min,z_max)[0] 
        return {'dNdz':dNdz_vec, 'z':zrange}
    elif type_field.lower() in ['galaxy', 'density', 'dens', 'source', 'overdensity', 'count', 'counting']:    
        dNdz_vec = np.array([ nz_spec(pars=pars)(ix) for ix in zrange ])
        if binning:
            if pars['survey'].upper()=='LSST': 
                dNdz_dict = source_bins_lsst(redshift_range=zrange, redshift_distribution=dNdz_vec, dens_pars=pars)     
                for i,ikey in enumerate(dNdz_dict.keys()):
                    if not i:
                        dNdz_vec = dcopy(dNdz_dict[ikey])
                    else:
                        dNdz_vec = np.vstack(( dNdz_vec,dcopy(dNdz_dict[ikey]) ))                
            elif pars['survey'].upper()=='EUCLID':                 
                nch   = pars['nbins']
                if int(pars['zbin_edges'].size-1)==int(pars['nbins']):
                    zbins = pars['zbin_edges']  
                else:
                    try:
                        zimin=pars['zbin_edges'].min()
                        zimax=pars['zbin_edges'].max()
                        bins = np.around( np.linspace(zimin, zimax, pars['nbins']), decimals=6)                    
                        pars['zbin_edges'] = bins
                        del zimin,zimax
                    except:
                        bins = np.around( np.linspace(z_min, z_max, pars['nbins']), decimals=6)                    
                        pars['zbin_edges'] = bins
                for i in range(nch):
                    zi_min, zi_max = zbins[i], zbins[i+1]
                    idNdz      = np.array([ integrate.quad(lambda izp: nz_convolved(zp=izp,pars=pars)(iz), zi_min, zi_max)[0] for iz in zrange ])
                    idNdz_intp = CubicSpline(zrange, idNdz)    
                    idNdz      = idNdz/integrate.quad(lambda iz : idNdz_intp(iz), z_min, z_max)[0]
                    if not i: dNdz_vec = dcopy(idNdz)
                    else: dNdz_vec = np.vstack(( dNdz_vec, dcopy(idNdz) )) 
                if verbose: print('Total processing time: {0:.2f} seg'.format(time.time()-timei)) 
            else:
                raise NameError(pars['survey']) 
        else:
            dNdz_vec = np.array([ nz_spec(pars=pars)(iz) for iz in zrange ])        
            dNdz_intp = CubicSpline(zrange, dNdz_vec)            
            dNdz_vec=dNdz_vec/integrate.quad( lambda iz: dNdz_intp(iz),z_min,z_max)[0] 
            
        return {'dNdz':dNdz_vec, 'z':zrange}
    else:
        raise NameError('There are no {} field implemented'.format(type_field.lower()))     

##############################################################################################################
##############################################################################################################
## LOAD GALAXY PARAMETERS
def load_survey_parameters(survey_name = 'LSST', data_release='Y10', type_field='density', mode = 'photometric',path_file=None):
    import yaml
    if survey_name.upper()=='LSST':
        path_file = '/data/AMARINS/LSSxHI-CODES/scripts/parameters/lsst_desc_parameters.yaml' \
                   if type(path_file)==type(None) else path_file
        load_pars = yaml.safe_load(open(path_file))
        if type_field.lower() in ['galaxy', 'density', 'dens', 'source', 'overdensity', 'count', 'counting']:
            dt = load_pars['source_sample'][mode.lower()][data_release.upper()] #data['source_sample'][data_release]
            try:
                dt['zbin_edges'] = np.asarray(dt['zbin_edges'])
            except:
                pass
            return dt
        if type_field.lower() in ['lensing','lens', 'shear']: 
            dt= load_pars['lens_sample'  ][data_release.upper()]   
            dt['zbin_edges'] = np.asarray(dt['zbin_edges'])   
            return dt        
    if survey_name.upper()=='EUCLID':
        path_file = '/data/AMARINS/LSSxHI-CODES/scripts/parameters/euclid_parameters.yaml' \
                   if type(path_file)==type(None) else path_file
        load_pars = yaml.safe_load(open(path_file))
        if type_field.lower() in ['galaxy', 'density', 'dens', 'source', 'overdensity', 'count', 'counting']:   
            dt = load_pars['source_sample'][mode.lower()]['Y10'] #data['source_sample'][data_release]
            try:
                dt['zbin_edges'] = np.asarray(dt['zbin_edges'])
            except:
                pass                    
            return dt
        if type_field.lower() in ['lensing','lens', 'shear']: 
            dt = load_pars['lens_sample'][data_release.upper()]  
            dt['zbin_edges']   = np.asarray(dt['zbin_edges'])
            return dt
##############################################################################################################
##############################################################################################################
# LOAD LSST GALAXY DIST
#
# Adapted from CCLX https://github.com/LSSTDESC/CCLX/blob/master/LSST%20SRD%20Redshift%20Distributions%20and%20Binning.ipynb
#
#####
def compute_equal_number_bounds_lsst(redshift_range, redshift_distribution, n_bins):
    from scipy.integrate import simpson, cumulative_trapezoid
    from scipy.special import erf    
    # Calculate the cumulative distribution
    cumulative_distribution = cumulative_trapezoid(redshift_distribution, redshift_range, initial=0)
    total_galaxies = cumulative_distribution[-1]

    # Find the bin edges
    bin_edges = []
    for i in range(1, n_bins):
        fraction = i / n_bins * total_galaxies
        # Find the redshift value where the cumulative distribution crosses this fraction
        bin_edge = np.interp(fraction, cumulative_distribution, redshift_range)
        bin_edges.append(bin_edge)
    return [redshift_range[0]] + bin_edges + [redshift_range[-1]]

def true_redshift_distribution_lsst(upper_edge, lower_edge, variance, bias, redshift_range, redshift_distribution):
    from scipy.special import erf     
    # Calculate the scatter
    scatter = variance * (1 + redshift_range)
    # Calculate the upper and lower limits of the integral
    lower_limit = (upper_edge - redshift_range + bias) / np.sqrt(2) / scatter
    upper_limit = (lower_edge - redshift_range + bias) / np.sqrt(2) / scatter

    # Calculate the true redshift distribution
    true_redshift_distribution = 0.5 * np.array(redshift_distribution) * (erf(upper_limit) - erf(lower_limit))

    return true_redshift_distribution

def source_bins_lsst(redshift_range, redshift_distribution, dens_pars, normalised=True):
    from scipy.integrate import simpson, cumulative_trapezoid
    from scipy.special import erf        
    # Get the bin edges as redshift values directly
    bins = compute_equal_number_bounds_lsst(redshift_range,redshift_distribution, dens_pars["nbins"])

    # Get the bias and variance values for each bin
    source_z_bias_list     = np.repeat(dens_pars["z_bias"] , dens_pars["nbins"])
    source_z_variance_list = np.repeat(dens_pars["sigma_z"], dens_pars["nbins"])

    # Create a dictionary of the redshift distributions for each bin
    source_redshift_distribution_dict = {}
    # Loop over the bins: each bin is defined by the upper and lower edge of the bin
    for index, (x1, x2) in enumerate(zip(bins[:-1], bins[1:])):
        z_bias     = source_z_bias_list[index]
        z_variance = source_z_variance_list[index]
        source_redshift_distribution_dict[index] = true_redshift_distribution_lsst(upper_edge=x1, lower_edge=x2, 
                                                                                   variance=z_variance, bias=z_bias, 
                                                                                   redshift_range=redshift_range, 
                                                                                   redshift_distribution=redshift_distribution)

    # Normalise the distributions
    if normalised:
        norm_factor = []
        for key in sorted(source_redshift_distribution_dict.keys()):
            norm_factor.append(simpson(source_redshift_distribution_dict[key], x=redshift_range))
            source_redshift_distribution_dict[key] /= norm_factor[-1]

    # Create a combined dictionary
    combined_data = {'redshift_range': redshift_range,'bins': source_redshift_distribution_dict}
    return source_redshift_distribution_dict

##############################################################################################################
##############################################################################################################
def get_bias_model(bias_model=None, bias_value=1, pars=None, path_to_file=None, z_range=None):
    #[1]
    if bias_model.upper() in ['LSST','EUCLID', 'SQRT']:
        return np.sqrt(1+z_range)
    #[2]    
    elif bias_model.upper() in ['F', 'FLAGSHIP', "FLAG", 'TUTUSAUS', 'TUTUSAUS2020', "FLAG1"]:
        #Euclid: The importance of galaxy clustering and weak lensing cross-correlations within the photometric Euclid survey - I Tutusaus et al (2020)
        coeffs = [1.0, 2.5, 2.8, 1.6] #[A,B,C,D]
        return coeffs[0] + coeffs[1]/(1+ np.exp(-(z_range-coeffs[3])*coeffs[2]))
    #[3] 
    elif bias_model.upper() in ['F2022', 'FLAGSHIP2022', "FLAG2022", 'F2', 'FLAGSHIP2', "FLAG2"]:
        #Euclid preparation: XIX. Impact of magnification on photometric galaxy clustering - F Lepori et al (2022)
        from numpy.polynomial import Polynomial as nPoly
        coeffs = [0.5125, 1.377, 0.222, -0.249]
        return nPoly(coeffs)(z_range)
    #[4]    
    elif bias_model.upper() in ['READ','FILE']:
        try:
            from scipy.interpolate import interp1d
            zi, bzi = np.loadtxt(path_to_file)            
            return interp1d(zi, bzi, kind='cubic')(z_range)
        except:
            return np.ones_like(z_range)
    #[5]            
    elif  bias_model.upper() in ['CTE','CONST', 'CONSTANT']:
        try:
            return bias_value*np.ones_like(z_range)
        except:
            return np.ones_like(z_range)       
    #[6]        
    else:
        return np.ones_like(z_range)      
                                
def get_slope_s_model(slope_s_model=None, slope_s_value=1, pars=None, path_to_file=None, z_range=None):
    if slope_s_model.upper() in ['F2022', 'FLAGSHIP2022', "FLAG2022", 'F2', 'FLAGSHIP2', "FLAG2",
                                 "LEPORI","LEPORI2022","EUCLID2022","EUCLIDXIX","EXIX","ECXIX"]:
        #Euclid preparation: XIX. Impact of magnification on photometric galaxy clustering - F Lepori et al (2022)
        from numpy.polynomial import Polynomial as nPoly
        coeffs = [0.0842,0.0532,0.298,-0.0113]
        return nPoly(coeffs)(z_range)  
    elif slope_s_model.upper() in ["DESHPANDE", "DESHPANDE2020", "DESH", "DESH2020"]:
        #Euclid: The reduced shear approximation and magnification bias for Stage IV cosmic shear experiments - AC Deshpande et al (2020)
        from numpy.polynomial import Polynomial as nPoly
        coeffs = [ 0.1547362 ,  0.17613424,  0.04163679,  0.35136411, -0.39927333, 0.19266138, -0.03429751]
        return nPoly(coeffs)(z_range)  
    elif slope_s_model.upper() in ['READ','FILE']:
        try:
            from scipy.interpolate import interp1d
            zi, si = np.loadtxt(path_to_file).T            
            return interp1d(zi, si, kind='cubic', fill_value='extrapolate')(z_range)
        except:
            return np.ones_like(z_range)
    elif  slope_s_model.upper() in ['CTE','CONST', 'CONSTANT']:
        try:
            return slope_s_value*np.ones_like(z_range)
        except:
            return np.ones_like(z_range)            
    else:
        return np.ones_like(z_range)                               

    
def magnification_bias(camb_params=None, camb_results=None, pars=None, z_range=None, 
                       slope_s_model=None, slope_s_value=1,  path_to_slope_file=None, 
                       binning=True, ibin=None, Np=500, verbose=False):
    slope_s = get_slope_s_model(slope_s_model=slope_s_model, pars=pars, path_to_file=path_to_slope_file, z_range=z_range)
    if type(ibin)==type(None):
        Wshear = kernel_galaxy_lensing_dict(camb_params=camb_params, Np=Np,
                                            camb_results=camb_results, 
                                            pars=pars,  binning=binning)['W']
    else: 
        Wshear = kernel_galaxy_lensing_i_dict(camb_params=camb_params, camb_results=camb_results, 
                                              Np=Np,pars=pars, ibin=ibin, binning=binning)['W']        
    Wshear*=(5*slope_s-2)
    return Wshear    

def get_charact_luminosity_admensional(L_model=None, L_value=1, pars=None, path_to_file=None, z_range=None):
    #<L>(z)/Lstar(z): Admensional Characteristic Luminosity
    #Euclid: The importance of galaxy clustering and weak lensing cross-correlations within the photometric Euclid survey
    #Eq 15
    #Euclid preparation - VII. Forecast validation for Euclid cosmological probes - Lesgourgues et al 2020
    #Eq 39
    #Euclid preparation: 6x2 pt analysis of Euclid's spectroscopic and photometric data sets - L Paganin et al 2024
    #Eq 13
    if L_model.upper() in ['LESGOURGUES', 'LESGOURGUES2025', 'BLANCHARD', 'BLANCHARD2020', 
                            'EUCLID', 'EUCLID2020','EVII', "ECVII", "EUCLIDVII", "EUCLID2020",
                            'ENLA']:
        from numpy.polynomial import Polynomial as nPoly
        coeffs = [ 0.15682055, -0.95500852,  5.5882817 , -9.60251292,  7.93685053,-2.99044314,  0.4147671 ]
        return nPoly(coeffs)(z_range)      
    elif L_model.upper() in ['READ','FILE']:
        #path_to_file = '../auxiliary/auxiliary/scaledmeanlum-E2Sa.dat'
        try:
            from scipy.interpolate import interp1d
            zc,Lc=np.loadtxt(path_to_file).T
            return interp1d(zc, Lc, kind='cubic', fill_value='extrapolate')(z_range)
        except:
            return np.ones_like(z_range)   
    elif L_model.upper() in ['CTE','CONST', 'CONSTANT']:
        try:
            return L_value*np.ones_like(z_range)
        except:
            return np.ones_like(z_range)           
    else:
        return np.ones_like(z_range)   
    
def get_F_IA_function(F_model=None, F_value=1, pars=None, path_to_file=None, z_range=None,#L_model=None, L_value=1, # A_IA=None, C_IA=None, 
                      beta_IA=None, eta_IA=None, path_to_IA_file=None):
    if F_model.upper() in ['LESGOURGUES',  'LESGOURGUES2025',  'BLANCHARD',  'BLANCHARD2020', 
                           'EUCLID', 'EUCLID2020','EVII', "ECVII", "EUCLIDVII", "EUCLID2020",
                           'ENLA']:
        Lc = get_charact_luminosity_admensional(L_model='EUCLID', pars=pars, z_range=z_range)
        #Lc = get_characteristic_luminosity_admensional(L_model='READ', pars=pars, z_range=z_range, path_to_file='../auxiliary/auxiliary/scaledmeanlum-E2Sa.dat')
        if type(path_to_IA_file)==type(None):
            #C_IA    = 0.0134 if type(   C_IA)==type(None) else C_IA
            #A_IA    = 1.72   if type(   A_IA)==type(None) else A_IA
            eta_IA  = -0.41 if type( eta_IA)==type(None) else eta_IA
            beta_IA = 2.17  if type(beta_IA)==type(None) else beta_IA
        else:
            import yaml
            try:
                F_pars = yaml.safe_load(open(path_to_IA_file))
            except:
                path='../auxiliary/IA_parameters.yaml'
                F_pars = yaml.safe_load(open(path))
            #C_IA    = F_pars[pars['survey']]['eNLA'][   'C_IA']            
            #A_IA    = F_pars[pars['survey']]['eNLA'][   'A_IA']            
            eta_IA  = F_pars[pars['survey']]['eNLA'][ 'eta_IA']            
            beta_IA = F_pars[pars['survey']]['eNLA']['beta_IA']      
        return ( (1+z_range)**eta_IA )*( Lc**beta_IA )
    elif F_model.upper()=='NLA':             
        return np.ones_like(z_range)
    elif F_model.upper() in ['READ','FILE']:
        try:
            from scipy.interpolate import interp1d
            zc,Fia=np.loadtxt(path_to_file).T
            return interp1d(zc, Fia, kind='cubic', fill_value='extrapolate')(z_range)
        except:
            return np.ones_like(z_range)   
    elif F_model.upper() in ['CTE','CONST', 'CONSTANT']:
        try:
            return F_value*np.ones_like(z_range)
        except:
            return np.ones_like(z_range)           
    elif F_model.upper() in ['FIT', 'FITTING']: 
        from numpy.polynomial import Polynomial as nPoly
        coeffs = [-1.44658943e-03, -2.30189882e-03,  4.81011268e+00, -4.11111700e+01,
                   1.41315751e+02, -2.49177926e+02,  2.51029131e+02, -1.50252007e+02,
                   5.29125726e+01, -1.01403713e+01,  8.17115026e-01] #Polynomial O(10)
        return nPoly(coeffs)(z_range)         
    else:
        return np.ones_like(z_range)       