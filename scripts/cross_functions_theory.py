import numpy as np

#####################################################################################################
#####################################################################################################
#### HI FUNCTIONS
####
def THI_factor_constants():
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

def hi_brightness_temperature(camb_params=None, camb_results=None, z=None, model = None):
    fact = (1+z)**2/(camb_results.hubble_parameter(z)/camb_params.H0)
    oHI  = omegaHI_biasHI(model=model,z=z)['omegaHI']*camb_params.h
    fact = fact*188.8#(THI_factor_constants()/1e-3)
    return fact*oHI #mK

def omegaHI_biasHI(model="jiajun",z=None):
    if model=="jiajun": #valido para z < 1.#Zhang et al 2022
        return {"omegaHI":{"HOD":2.7e-4 + 1e-4*z - 8e-5*z**2,
                           "HAM":2.5e-4 - 4e-5*z - 7e-5*z**2}
               }
    if model=="padmanabhan":#padmanabhan 2015
        z   = np.array([0.000, 0.250, 0.500, 0.750, 1.000, 1.250, 1.500, 1.750, 2.000, 2.250, 2.500, 2.750, 3.000, 3.250, 3.400])
        oHI = np.array([3.344, 3.443, 4.523, 4.648, 4.710, 4.804, 4.766, 4.804, 4.936, 5.008, 4.750, 5.471, 5.541, 5.756, 5.971])
        bHI = np.array([0.703, 0.972, 1.026, 0.935, 1.005, 1.005, 1.049, 1.099, 1.101, 1.160, 1.261, 1.409, 1.329, 1.498, 1.802])        
        return {"z":z, "omegaHI":oHI, "bias":bHI}
    if model=="cunnington":#cunnington2019
        return {"omegaHI":0.00048 + 0.00039*z - 0.000065*z**2, 'biasHI': 1}
    if model=="irfan":#Irfan2021 baseado em bull2015 #o mesmo que IRFAN2021
        omegaHI_0 = 4.86e-4 #Irfan2021
        biasHI_0  = 0.677105 #Irfan2021
        return {"omegaHI": (omegaHI_0/4.86)*(4.8304 + 3.8856*z - 0.65119*z**2) ,
                "biasHI" : (biasHI_0/0.677105)*(0.66655 + 0.17765*z + 0.050223*z**2)}
    if model=="crighton":#carucci2020
        omegaHI_0 = 4e-4 #Irfan2021
        biasHI_0  = 0.3 #Irfan2021        
        return {"omegaHI": omegaHI_0*(1+z)**(0.6) ,
                "biasHI" : biasHI_0*(1+z)+0.6}
    else:
        return {"omegaHI": 4.86e-4,
                "biasHI" : 1}        
#####################################################################################################
#####################################################################################################
#### KERNELS
####
def window_1(z=None, zmin=None, zmax=None):
    if (z>=zmin)*(z<=zmax): return 1/(zmax-zmin)
    else: return 0
    
def chi_vec(camb_results=None): #radial comovel distance
    from scipy import interpolate
    zv         = np.linspace(0,1200)
    rv         = camb_results.comoving_radial_distance(zv)
    z_chi      = interpolate.interp1d(rv, zv, kind='linear') #function z(chi)
    z_star     = camb_results.get_derived_params()['zstar']
    chi_star   = camb_results.comoving_radial_distance(z_star)
    #zs         = results.redshift_at_comoving_radial_distance(np.linspace(0,chi_star,nz))
    z_chi_camb = lambda y: results.redshift_at_comoving_radial_distance(y)    
    return {'zvec':zv, 'chivec':rv, 'z_chi': z_chi, 'z_chi_camb':z_chi_camb, 'z_star':z_star, 'chi_star':chi_star}

def kernel_cmb_z(camb_params=None, camb_results=None, z=None, norm=False):
    c        = 299792458/1e3
    fact     = 1.5*camb_params.omegam*camb_params.H0**2/c/camb_results.hubble_parameter(z)
    chi_star = chi_vec(camb_results)['chi_star']
    chi_z    = camb_results.comoving_radial_distance(z)
    chi_z    = fact*(1+z)*chi_z*(chi_star-chi_z)/chi_star
    if norm: A = 1#chi_z.max()
    else:    A = 1
    return chi_z/A

def kernel_cmb_chi(camb_params=None, camb_results=None, chi=None, norm=False):
    c        = 299792458/1e3
    fact     = 1.5*camb_params.omegam*camb_params.H0**2/(c**2)
    chi_dict = chi_vec(camb_results)
    chi_star = chi_dict['chi_star']
    chi_star = fact*(1+chi_dict['z_chi'](chi))*chi*(chi_star-chi)/chi_star
    if norm: A = 1#chi_star.max()
    else:    A = 1
    return chi_star/A    

def kernel_hi(camb_params=None, camb_results=None, z=None, zmin = None, zmax=None, model=None):
    w1  = window_1(z=z, zmin=zmin, zmax=zmax)
    if not w1: return 0
    THI = hi_brightness_temperature(camb_params=camb_params, camb_results=camb_results, z=z, model=model)
    bHI = omegaHI_biasHI(model=model,z=z)['biasHI']
    return bHI*THI*w1
    
def kernel_fields(params=None, camb_params=None, camb_results=None, z=None, type=None):#
    return 0

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

#####################################################################################################
#####################################################################################################
#### SAVING
####
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
    elif field.lower()=='cross':
        plt.loglog(clf_dict['l'],clf_dict['bin1'],linestyle='solid', linewidth=4, color='grey', label='bin1')
        for i,iname in enumerate(clf_dict.keys()):
            if iname!='l':plt.loglog(clf_dict['l'],clf_dict[iname],linestyle='solid', color=cm.Blues(i/len(clf_dict.keys())))	        
    else: raise Exception('NO field like that.')		    
    plot_type = ['log','log']
    plt.ylabel(r'$C_{\ell}^{\small \textrm{HI} -\small \textrm{HI}}$')
    plt.xlabel('$\ell$')
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
    nu_vec    = hdata.nu_bins_vector(numin_ =numin, numax_ =numax, nbands_= nbands)['nu']
    nuvec_eff = np.array([ 0.5*(nu_vec[i]+nu_vec[i+1]) for i in range(nu_vec.size-1) ])
    ibins     = np.linspace(0, nbands-1, n_curves, dtype=np.int16)
    return {'zeff': np.around(1420.405751768/np.array([nuvec_eff[ibins]]) -1, decimals=2)[0],
            'bins':ibins, 
            'nu_eff':nuvec_eff[ibins]}