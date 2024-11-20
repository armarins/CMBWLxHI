import numpy as np
def noise_parameters(nside=256, nu_min_MHz=980, nu_max_MHz=1260, nch=30, 
                     nbeams=28, Tsys=70, Osur=5324, Obeam=0.35, tsur_year=1, 
                     K=1.414  , fsky=0.13, dcycle=1, verbose=False, unit_mK=True):
    from copy import deepcopy as dcopy
    #nside  = 256
    #npix   = 12*nside**2
    #nu_min = 980  #Mhz
    #nu_max = 1260 #Mhz
    #nch    = 30
    #nbeams = 28       #number of beams# here it is the number of feed horns
    #Tsys   = 70       #system temperatyre in K
    #Osur   = 5324     #Full survey area           # sqr deg
    #Obeam  = 0.35     #Telescope beam solid angle # sqr deg
    #tsur   = 1        #Mission duration # yrs
    #K      = 2**(1/2) #two circ polarization contribution: (I-V) and (I+V)
    #fsky   = 0.13
    #dcycle = 0.9
    
    npix   = 12*nside**2
    tsur   = dcopy(tsur_year)
    nu_min = dcopy(nu_min_MHz)
    nu_max = dcopy(nu_max_MHz)

    t_unit  = 24*60*60
    tsur   *= 365*t_unit #yr to sec
    tsur   *= dcycle
    nu_max *= 1e6
    nu_min *= 1e6

    bandwidth   = (nu_max-nu_min)/nch
    N           = Osur/Obeam #number of pixels in the map
    tpix        = (tsur/N)*nbeams
    
    Snoise = K*Tsys/np.sqrt(tpix*bandwidth)
    Spix   = Snoise*np.sqrt(fsky*npix/nbeams/tsur) 
    if verbose:
        #print("tpix {:.2f} sec".format(tpix))
        print("tpix: {:.2f} hour/pix".format(tpix/3600),"\n")
        #print("sigmaN: {:.8f} K".format(Snoise))
        print("sigmaN: {:.2f} mK".format(Snoise*1e6),"\n")
        #print("sigma pix: {:.8f} K".format(Spix))
        print("sigma pix: {:.2f} mK/pix".format(Spix*1e6),"\n")
    if unit_mK: A=1e6
    else:       A=1
    return {'tpix':  tpix/3600,
            'sigmaN': Snoise*A,
            'sigma pix': Spix*A,
            'npix':npix
            }


def theoretical_white_noise_Cl(nside=256, nu_min_MHz=980, nu_max_MHz=1260, nch=30, nbeams=28, Tsys=70, Osur=5324, Obeam=0.35, tsur_year=1, 
                     K=1.414  , fsky=0.13, dcycle=1, lmin=1,lmax=400, del_l=1, verbose=False, unit_mK=True):
    npars = noise_parameters(nside, nu_min_MHz, nu_max_MHz, nch, nbeams, Tsys, Osur, Obeam, tsur_year, K, fsky, dcycle, verbose, unit_mK)
    return 4*np.pi*((npars['sigma pix']**2)/npars['npix'])*np.ones_like(np.arange(lmin, lmax+del_l, del_l))

def healpy_white_noise_Cl(nside=256, nu_min_MHz=980, nu_max_MHz=1260, nch=30, 
                     nbeams=28, Tsys=70, Osur=5324, Obeam=0.35, tsur_year=1, 
                     K=1.414  , fsky=0.13, dcycle=1, verbose=False, unit_mK=True):
    return None


def Tsys_SKA_MID_band1(nu_eff=None):#[nu]=MHz. [Tsys]=K
    Tcmb = 2.73#K
    Tspl = 3.00#K
    Tgal = 25*np.power(408/nu_eff,2.75) #K
    Trx  = 15 + 30*np.power(nu_eff/1e3 - 0.75, 2) #K
    Tsys = Tcmb+Tspl+Tgal+Trx
    return Tsys


def bl_function(fwhm=None, lmax=None, theta_=None, input_unit="radian", from_real_space=False): 
    if input_unit=="radian":
        pass
    elif input_unit=="degree":
        fwhm = np.radians(fwhm)
        try:
            theta_ = np.radians(theta_)
        except:
            pass
    elif input_unit=="arcmin":
        fwhm = np.radians(fwhm/60)
        try:
            theta_ = np.radians(theta_/60)
        except:
            pass
    else:
        pass   
    if from_real_space:
        bl  = model.beam_function(type_,fwhm=fwhm, theta_=theta_)
        bl = hp.beam2bl(bl,theta_,lmax)    
        bl = bl/bl.max()
    else:
        if type(lmax)==type(None):lmax=768
        l  = np.arange(lmax+1)
        bl = np.exp(-0.5*l*(l+1)*fwhm**2/8/np.sqrt(2))
    return bl


###################################################################################
########## Leakage
def leakage_theorectical(ins_=None,
                        pathdir_lkg_gal='/data/AMARINS/CMBWLxHI-DATA/leakage/fullsky/lowz_nch30_980_1260_nch25_30',
                        pathdir_lkg_cmb='/data/AMARINS/CMBWLxHI-DATA/leakage/fullsky/lowz_nch30_980_1260',):
    # LOWZ_gal: '/data/AMARINS/CMBWLxHI-DATA/leakage/fullsky/lowz_nch30_980_1260_nch25_30/{}/'.format(ins_)
    # LOWZ_cmb: '/data/AMARINS/CMBWLxHI-DATA/leakage/fullsky/lowz_nch30_980_1260/{}/'.format(ins_): 
    # HIGHZ gal: '/data/AMARINS/CMBWLxHI-DATA/leakage/fullsky/highz_nch70_350_1050_nch65_70/{}/'.format(ins_)
    # HIGHZ cmb: '/data/AMARINS/CMBWLxHI-DATA/leakage/fullsky/highz_nch70_350_1050/{}/'.format(ins_)
    return {'gal': np.loadtxt(os.path.join(pathdir_lkg_gal,ins,'leakage_theory.txt')).T,
            'cmb': np.loadtxt(os.path.join(pathdir_lkg_cmb,ins,'leakage_theory.txt')).T}
###################################################################################
########## Signal-to-Noise

def SN_per_ell_per_nu(CLs_dict_=None, ins_=None, clcx_=None,leakage=False, leakage_type='gal',
                      ibin=None, cx_label= 'cl_cx_f1_f2_rec_sim', f1_label='cl_f1_rec', f2_label='cl_f2_sim', hi_regime='lowz',
                      noise_f1=0, noise_f2=0, beam_f1=0, beam_f2=0, fsky=1, delell=1, elleff=None, verbose=False):
    timej=time.time()    
    if leakage:
        lkg = leakage_theorectical(ins_=ins_, hi_regime=hi_regime)
        #lkg = leakage_theorectical(CLs_dict_=CLs_dict_, ins_=ins_, clcx_=clcx_, bins_=bins_, calculate=calculate)
        lkg = lkg[leakage_type][ibin]
    else:
        lkg=0
    try:
        noise_f1=noise_f1[ibin]
    except:
        noise_f1=0
    try:
        noise_f2=noise_f2[ibin]
    except:
        noise_f2=0       
    cx_rec_part = b.bin_cell(CLs_dict_[ins_][cx_label][ibin] + lkg)
    c1_rec_part = b.bin_cell(CLs_dict_[ins_][f1_label][ibin] + noise_f1/(beam_f1**2))
    c2_sim_part = b.bin_cell(CLs_dict_[ins_][f2_label][0]    + noise_f2/(beam_f2**2))
    err_rec     = np.sqrt( cx_rec_part**2   + c1_rec_part*c2_sim_part  )
    err_rec     = err_rec/np.sqrt( (2*elleff+1)*fsky*delell )
    if verbose: print('Processing time: {0:.4f} seg'.format(time.time()-timej))     
    return cx_rec_part/err_rec


def cumulative_SN(CLs_dict_=None, clcx_=None, bins_=None, ins_=None, sim_names_=None, 
                  leakage=False, leakage_type='gal', hi_regime='lowz', type_field='gal',
                  cx_label= 'cl_cx_f1_f2_rec_sim', f1_label='cl_f1_rec', f2_label='cl_f2_sim',
                  noise_f1=0,noise_f2=0, beam_f1=0,beam_f2=0, fsky=1, delell=1,elleff=None,verbose=False):
    timej=time.time()    
    for j,jsim in enumerate(sim_names_):
        cls_dict = {ins:cl_dict_function(L0_dir=jsim, ins_=ins_,type_field=type_field, hi_regime=hi_regime, verbose=False)}
        for jj,ibin in enumerate(bins_):
            sn = SN_per_ell_per_nu(CLs_dict_=CLs_dict_, clcx_=clcx_, ibin=ibin, ins_=ins_,
                                          leakage  = leakage , calculate = calculate, leakage_type=leakage_type,
                                          cx_label = cx_label, f1_label  = f1_label , f2_label=f2_label, hi_regime=hi_regime,
                                          noise_f1 = noise_f1, noise_f2  = noise_f2 , 
                                          beam_f1  = beam_f1 , beam_f2   = beam_f2  , 
                                          fsky=fsky, delell=delell, elleff=elleff, 
                                          verbose=verbose)
            #sn = np.sqrt(sn**2)
            sn_m = np.sum(sn) if not jj else np.hstack(( sn_m, np.sum(sn) ))
        SN = sn_m if not j else np.hstack(( SN, np.sum(sn_m) ))
    SN = np.average( SN , axis=0 ).reshape(-1,len(bins_)) 
        
    if verbose: print('Processing time: {0:.4f} seg'.format(time.time()-timej))     
    return SN


def SN_per_ell_per_nu_mean(ins_=None, clcx_=None, bins_=None,
                           leakage=False, leakage_type='gal', hi_regime='lowz', type_field = 'gal', sim_names_=None,
                           cx_label= 'cl_cx_f1_f2_rec_sim', f1_label='cl_f1_rec', f2_label='cl_f2_sim', 
                           noise_f1=0, noise_f2=0, beam_f1=0, beam_f2=0, fsky=1, delell=1, elleff=None, verbose=False):
    for j,jbin in enumerate(bins_):    
        for jj,jsim in enumerate(sim_names_):
            cls_dict_ = {ins_:cl_dict_function(L0_dir=jsim, ins_=ins_,
                                              type_field=type_field, 
                                              hi_regime=hi_regime, 
                                              verbose=verbose)}
            sn_ = SN_per_ell_per_nu(CLs_dict_=cls_dict_, ins_=ins_, 
                                    clcx_=clcx_, hi_regime=hi_regime,
                                    leakage_type=lkg_field, leakage=leakage, ibin=jbin, 
                                    cx_label= cx_label, f1_label=f1_label, f2_label=f2_label,
                                    noise_f1=noise_f1, beam_f1=beam_f1, 
                                    noise_f2=noise_f2, beam_f2=beam_f2, 
                                    delell=delell, elleff=elleff, fsky=fsky,
                                    verbose=verbose)
            sn_vec = sn_ if not jj else np.vstack(( sn_vec, sn_ ))
            del sn_, cls_dict_
        sn_vec = np.average( sn_vec , axis=0 )
        SN_ = sn_vec if not j else np.vstack(( SN_, sn_vec ))
        del sn_vec
    return SN_


def SN_per_nu_mean(ins_=None, clcx_=None, bins_=None, NS=['ns3','ns4','ns5'],
                   leakage=False, leakage_type='gal', hi_regime='lowz', type_field = 'gal', sim_names_=None,
                   cx_label= 'cl_cx_f1_f2_rec_sim', f1_label='cl_f1_rec', f2_label='cl_f2_sim', 
                   noise_f1=0, noise_f2=0, beam_f1=0, beam_f2=0, fsky=1, delell=1, elleff=None, verbose=False):
    for j,jns in enumerate(NS):
        sn_ = SN_per_ell_per_nu_mean(ins_=jns, clcx_=clcx_, bins_=bins_,
                                  leakage=leakage, leakage_type=leakage_type, hi_regime=hi_regime, type_field =type_field, sim_names_=sim_names_,
                                  cx_label=cx_label, f1_label=f1_label, f2_label=f2_label, 
                                  noise_f1=noise_f1, noise_f2=noise_f2, beam_f1=beam_f1, beam_f2=beam_f2, fsky=fsky, delell=delell, elleff=elleff, verbose=verbose)
        sn_ = np.sqrt(sn_**2)
        sn_ = np.sum(sn_,axis=1)
        if not j:
            SN = {jns: sn_}
        else:
            SN[jns] = sn_
        del sn_
    return SN
    
    #for j,jbin in enumerate(bins_):   