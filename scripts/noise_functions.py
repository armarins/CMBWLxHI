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
    