import numpy as np
import os, sys
import healpy as hp
import astropy.io.fits as fits
from copy import deepcopy as dcopy
#HImap-r0023-f1z15.fits
def nu_bins_vector(numin_=None, numax_=None, nbands_=None):
    nustep_ = (numax_-numin_)/nbands_
    nu_     = np.around(np.arange(numin_, numax_ + nustep_, nustep_),decimals=2)
    nbands  = nu_.size-1
    return {"min":numin_,"max":numax_,"bandwidth":nustep_, "nbands":nbands, "nu":nu_}

def get_info_from_directory(dirname_= None, kind="flask"):
    if kind=="flask":  #formatting: ../flask_maps_[NSIDE_NUMBER]
        nside  = int(dirname_.split("_")[-1])
        return {"NSIDE":nside}
    else: raise NameError
    
def get_info_from_filename(dirname_= None, filename_=None, numin_=980,numax_=1260,nbands_=30, kind="flask"):
    if kind=="flask":  #formatting: ../flask_maps_[NSIDE_NUMBER]/HImap-r[REALISATION]-f[FIELD_NUMBER]z[REDSHIFT_BIN].fits
        nudata = nu_bins_vector(numin_=numin_, numax_=numax_, nbands_=nbands_)
        nside  = int(dirname_.split("_")[-1])
        _,realization, zbin = filename_.split("-")
        realization = realization[1:]
        zbin  = int(zbin.split(".fits")[0].split("z")[1])
        nubin = np.flip(np.flip(nudata['nu'])[zbin-1:zbin+1])
        #zbin, nubin        
        return {"NSIDE":nside,
                "realization":realization,
                "frequency":{"min":nudata['min'],"max":nudata['max'], "bandwidth":nudata['bandwidth'], 'nbands':nudata['nbands'],"nu":nudata['nu']},
                "bin":{"numin":nubin[0],"numax":nubin[1],"z":zbin},
                'dirname':dirname_,"filename":filename_
               }
    else: raise NameError
    
def getmap(dirpath_=None, filename_=None, healpix_readingformat=True, hdu=0):
    if healpix_readingformat:
        maps = hp.read_map(os.path.join(dirpath_,filename_))
    else:
        import astropy.io.fits as fits
        with fits.open(os.path.join(dirpath_,filename_)) as h:
            maps = h[hdu].data
    return maps


def get_allfilenames(dirpath_=None):
    names = np.array([], dtype=str)
    for i,iname in enumerate(os.listdir(dirpath_)):
        if "fits"==iname.split(".")[1]:
            names = np.hstack((names,iname))
    return names        

def get_filenames(dirpath_=None):
    filenames = get_allfilenames(dirpath_=dirpath_)
    names = np.array([], dtype=str)
    for i,iname in enumerate(os.listdir(dirpath_)):
        if "fits"==iname.split(".")[1]:
            names = np.hstack((names,iname))
    return names     

def get_realization_id_from_filename(filename_=None, new_format=False):
    if new_format:
        return filename_.split("_L")[1].split(".fits")[0]
    else:
        return filename_.split('-')[1].split('r')[1]
def get_all_realizations_id(filenames_=None, new_format=False):
    IDs = np.array([])
    for i,ifilename in enumerate(filenames_):
        _id_ =get_realization_id_from_filename(filename_=ifilename,  new_format=new_format)
        IDs = np.append(IDs,_id_)
    return np.unique(IDs)    

def get_filenames4realizations(filenames_=None,ids_=None, new_format=False): #ids_: all realization name index. Fex: ['0001', '0002', ..., '0100']
    for i,_id_ in enumerate(ids_):
        if not i:
            idic = {_id_:np.array([])}
        else:
            idic[_id_] = np.array([])
    for i,ifilename in enumerate(filenames_):
        _id_ = get_realization_id_from_filename(filename_=ifilename,  new_format=new_format)
        idic[_id_] = np.hstack((idic[_id_],ifilename))
    for _id_ in idic.keys():
        idic[_id_] = np.sort(idic[_id_])
    return idic

def building_cubemaps(pathdir_=None, dirname_=None, healpix_readingformat=True):
    from copy import deepcopy as dcopy
    filenames = get_allfilenames(dirpath_=pathdir_)
    _ids_     = get_all_realizations_id(filenames)
    idic      = get_filenames4realizations(ids_=_ids_)
    i = np.random.choice(_ids_)
    HI = building_cubemaps_from_realization(pathdir_=pathdir, filenames_=idic[i],)    

def building_cubemaps_from_realization(pathdir_=None, filenames_=None, healpix_readingformat=True):#this's for one id and not for all of them
    from copy import deepcopy as dcopy
    dirname_ = pathdir_.split("/")[-1]
    ginfo    = get_info_from_filename(dirname_= dirname_, filename_=filenames_[0])
    HImaps   = np.ones((ginfo['frequency']['nbands'],12*(ginfo['NSIDE'])**2))
    for i,ifilename in enumerate(filenames_):
        ifilename = os.path.join(pathdir_, ifilename)
        numname   = int(ifilename.split(".fits")[0].split('z')[1])-1 #names in z between 1-30, but cube between 0-29 numeration
        imap      = getmap(dirpath_=pathdir_, filename_=ifilename, healpix_readingformat=healpix_readingformat)
        HImaps[numname,:] = dcopy(imap)   
    return HImaps
        
    
def new_formattingnames(F=None, S=None, NSIDE_=None, 
                        freq_min=None, freq_max=None,freq_unit="mhz",
                        Num_=None, C=None, Bres=None, Bmodel=None,
                        A=None, R=None):
    #STANDARD values:       #OBS: Read add pdf attached
    #F='HI'
    #S='I'
    #NSIDE_=512
    #freq_min=980
    #freq_max=1260
    #freq_unit="mhz"
    #Num_=30
    #C='partial'
    #Bres=40
    #Bmodel='G'
    #A=None
    #R='0001'    
    if Bres==None or Bmodel==None:
        if A==None:
            return "{}_{}_{}_{}{}{}{}_{}bins_{}_L{}.fits".format(F,S,NSIDE_,freq_min, freq_unit, freq_max, freq_unit,
                                                                       Num_, C, R)         
        else:
            return "{}_{}_{}_{}{}{}{}_{}bins_{}_{}_L{}.fits".format(F,S,NSIDE_,freq_min, freq_unit, freq_max, freq_unit,
                                                                       Num_, C,A, R)         
    else:    
        if A==None:
            return "{}_{}_{}_{}{}{}{}_{}bins_{}_{}arcmin{}_beam_L{}.fits".format(F,S,NSIDE_,freq_min, freq_unit, freq_max, freq_unit,
                                                                       Num_, C, Bres, Bmodel, R)         
        else:
            return "{}_{}_{}_{}{}{}{}_{}bins_{}_{}arcmin_{}beam_{}_L{}.fits".format(F,S,NSIDE_,freq_min, freq_unit, freq_max, freq_unit,
                                                                       Num_, C, Bres, Bmodel, A, R) 

def get_beam_model(acron):
    if acron.upper() == "ZN":
        return "Zernike"
    elif acron.upper() == "G":
        return "Gaussian"
    elif acron.upper() == "JC":
        return "Spherical Bessel"
    elif acron.upper() == "CO":
        return "Cosine"
    elif acron.upper() == "LA":
        return "Laguerre"
    elif acron.upper() == "HM":
        return "Hermite"
    elif acron.upper() == "AI":
        return "Airy"
    else:
        return None


def get_beam_model_inv(acron):
    if acron.lower() == "zernike":
        return "ZN"
    elif acron.lower() == "gaussian":
        return 'G'
    elif acron.lower() == "spherical bessel" or acron.lower() == "sphericalbessel":
        return 'JC'
    elif acron.lower() == "cosine":
        return 'CO'        
    elif acron.lower() == "laguerre":
        return 'LA'
    elif acron.lower() == "hermite":
        return 'HM'
    elif acron.lower() == "airy":
        return 'AI'
    else:
        return None

        
#HI_I_512_980mhz1260mhz_30bins_partial_40arcminZNbeam_L0010.fits
def get_info_from_output_filename(filename_=None):
    #HI_I_512_980mhz1260mhz_30bins_partial_40arcminZNbeam_L0010.fits
    iparts = filename_.split("_")
    field_info  = iparts[0]
    stokes_info = iparts[1]
    nside_info  = int(iparts[2])
    freq_info   = np.asarray(iparts[3].split("mhz")[:-1],dtype=np.float64)
    nbins_info  = int(iparts[4].split("bins")[0])
    cover_info  = iparts[5]
    if len(iparts[6].split('beam'))>1:
        beam_info   = iparts[6]
        fwhm_info   = float(beam_info.split("arcmin")[0])
        type_of_beam = get_beam_model(beam_info.split("arcmin")[1].split("beam")[0])
    else:
        type_of_beam   = None
        fwhm_info      = None
    realiz_info = iparts[-1].split("L")[1]
    return {"field": field_info, "stokes":stokes_info, "NSIDE":nside_info, "freq":freq_info,
            "nbins":nbins_info, "coverage":cover_info, "beam":{"fwhm":fwhm_info, "model":type_of_beam}}

def creating_multimap_FITSheader(filename_=None, freq_unit="Mhz", stokes_unit="mk", fwhm_unit="arcmin",realization_=None):
    import datetime
    ginfo = get_info_from_output_filename(filename_)
    hdr = fits.Header()
    hdr['COMMENT'] = "THIS IS A NEW NAME FORMAT FOR THE BINGO PROJECT FITS FILES."
    hdr['COMMENT'] = "THESE FILES ARE USED FOR COMPONENT SEPARATION CODE CALLED --CHISEL--"
    hdr['COMMENT'] = "https://github.com/multinverse/CHISEL---signals-from-LIM"
    hdr['COMMENT'] = "THATS FORMAT ASSUMES FLASK MAPS AS INPUT"
    hdr['COMMENT'] = "CREATED IN: {}".format(str(datetime.datetime.now()).split(".")[0])
    for i,ikey in enumerate(ginfo.keys()):
        if ikey=="field":
            hdr['field'] = ginfo['field']
            #hdr.comments['field'] = "{}".format("...")
        if ikey=="stokes":
            hdr["stokes"] = ginfo["stokes"]
            hdr.comments["stokes"] = "{}".format(stokes_unit)       
        if ikey=="coverage":
            hdr["coverage"] = ginfo["coverage"]
            hdr.comments["coverage"] = "{full- or partial-sky}"
        if ikey=="freq":
            if len(ginfo['freq'])==2:
                hdr["_".join((ikey,"min"))]=np.amin(ginfo['freq'])
                hdr["_".join((ikey,"max"))]=np.amax(ginfo['freq'])
                hdr.comments["_".join((ikey,"min"))]="in {}".format(freq_unit)
                hdr.comments["_".join((ikey,"max"))]="in {}".format(freq_unit)
        elif ikey=="beam":
            if ginfo[ikey]==None or ginfo[ikey]['fwhm']==None or ginfo[ikey]['model']==None:
                hdr["beam"]= 'NONE'
            else:
                hdr['fwhm'] = ginfo['beam']['fwhm']
                hdr['beam'] = ginfo['beam']['model']
                hdr.comments['fwhm']="{}".format(fwhm_unit)
        else:
            hdr[ikey]=ginfo[ikey]
    hdr['REALIZATION']=realization_
    return hdr 
    
def creating_mixmatrix_FITSheader(): 
    hdr = fits.Header()
    hdr['COMMENT']="SECONDARY HDU CONTAINS THE A-MIXMATRIX"
    return hdr
    
def creating_primary_FITSheader(output_info=None):
    hdr = fits.Header()
    hdr['COMMENT']="PRIMARY HDU CONTAINS THE REDSHIFT VECTOR"
    hdr['FREQ_MIN'] = output_info['frequency']['min']
    hdr['FREQ_MAX'] = output_info['frequency']['max']
    hdr['NBANDS'] = output_info['frequency']['nbands']
    return hdr
def return_new_FITSfilename(ginfo=None, output_info=None, add_info=None):
    #ginfo['NSIDE']
    #ginfo['realization']
    #ginfo['frequency']['min']
    #ginfo['frequency']['max']
    #ginfo['nbands']
    #output_info['field']
    #output_info['covering']
    #output_info['beam']['model']
    #output_info['beam']['FWHM']
    #output_info['stokes']    
    output_info['beam_resolution'] = None if not output_info['apply_beam'] else str(int(output_info['fwhm']*60))
    output_info['beam_model']      = None if not output_info['apply_beam'] else output_info['beam']['model']    
    return new_formattingnames(F=output_info['field'], S=output_info['stokes'], 
                        NSIDE_  =ginfo['NSIDE'], 
                        freq_min=ginfo['frequency']['min'],   freq_max=ginfo['frequency']['max'],freq_unit="mhz",
                        Num_    =ginfo['frequency']['nbands'], C=output_info['coverage'], 
                        Bres    =output_info['beam_resolution'], Bmodel=output_info['beam_model'],
                        R       =ginfo['realization'],
                        A       =add_info)
        
        
def file_verification(path,filename,directory):
    if type(directory)==int or type(directory)==float:
        directory = str(directory)
    elif type(directory)==str:
        pass
    else:
        #raise(sys.exc_info()[0])
        print(directory)
    if os.path.isdir(path):
        path = os.path.join(path,directory) 
        if os.path.isdir(path):
            path = os.path.join(path,filename)
            if os.path.isfile(path):
                os.remove(path)        
        else:
            os.mkdir(path)
    else:
        os.mkdir(path)
        os.mkdir(os.path.join(path,directory))

        
def save_FITS_newformat(pathdir_=None, output_info=None):
    from copy import deepcopy as dcopy
    if output_info['field'].upper()=='HI':
        idir      = pathdir_.split('/')[-1]
        filenames = get_allfilenames(dirpath_=pathdir_) #all filenames in the directory
        _ids_     = get_all_realizations_id(filenames) #all realization names
        idic      = get_filenames4realizations(filenames_=filenames, ids_=_ids_) #all filenames for a specific realization
        
        if output_info['output_dir']==None:
            file_verification(os.getcwd(),"","outputs")
            pathout = os.path.join(os.getcwd(),"outputs")
        else:
            pathout = output_info['output_dir']
        if output_info['clear_output_dir']:
            cmd = "rm -rf {}/*".format(pathout)
            os.popen(cmd).readlines()                      
                    
        ###########
        if not (type(output_info['realization'])==type(None) and output_info['realization']=='all'):
            if type(output_info['realization'])==str:
                output_info['realization'] = np.array([output_info['realization']])
           # elif type(output_info['realization'])==np.int and output_info['realization']<=len(_ids_):
            #    reals = _ids_[:int(output_info['realization'])]
            if (type(output_info['realization'])==list or type(output_info['realization'])==np.ndarray):
                reals = output_info['realization']
            else:
                reals = _ids_
            reals = np.sort(reals)
            print("Realizations to be saved: {}\n\n".format(reals))
            for i, ireal in enumerate(reals):
                print("Realization {}".format(ireal))
                ifile     = np.random.choice(idic[ireal]) #Choose a specific realization name just to test
                ginfo     = get_info_from_filename(dirname_= idir, filename_=ifile, 
                                                       numin_ =output_info['frequency']['min'], 
                                                       numax_ =output_info['frequency']['max'], 
                                                       nbands_=output_info['frequency']['nbands']) #info from filenames
                #new_name  = return_new_FITSfilename(ginfo, output_info, add_info=None) #new FITS name for that map cube
                try:
                    new_name = return_new_FITSfilename(ginfo=ginfo, output_info=output_info, add_info=output_info['add_info'])
                except:
                    new_name = return_new_FITSfilename(ginfo=ginfo, output_info=output_info, add_info=None)
                print("Filename: {}".format(new_name))
                vec_hdu0  = ginfo['frequency']['nu']
                hdr_hdu0  = creating_primary_FITSheader(output_info)
                vec_hdu1  = building_cubemaps_from_realization(pathdir_=pathdir_, filenames_=idic[ireal]) #extracting all maps from a specific realization and for building a cube
                hdr_hdu1  = creating_multimap_FITSheader(new_name, 
                                                               freq_unit   = output_info["freq_unit"], 
                                                               stokes_unit = output_info["stokes_unit"], 
                                                               fwhm_unit   = output_info["fwhm_unit"],
                                                               realization_= ireal) #Header to be used 
                try:
                    if int(output_info['savedata_hdu'])==0:            
                        hdu0 = fits.PrimaryHDU(header=hdr_hdu1,data=vec_hdu1)
                        hdul = fits.HDUList([hdu0])
                    elif int(output_info['savedata_hdu'])==1:            
                        hdu0 = fits.PrimaryHDU(header=hdr_hdu0,data=vec_hdu0)
                        hdu1 = fits.ImageHDU(  header=hdr_hdu1,data=vec_hdu1, name="MULTIMAPS")
                        hdul = fits.HDUList([hdu0,hdu1])                
                    else:
                        raise Expection("Savedata HDU must be 0 or 1. Current: {}".format(output_info['savedata_hdu']))
                except:
                    hdu0 = fits.PrimaryHDU(header=hdr_hdu0,data=vec_hdu0)
                    hdu1 = fits.ImageHDU(  header=hdr_hdu1,data=vec_hdu1, name="MULTIMAPS")
                    hdul = fits.HDUList([hdu0,hdu1])                            
                hdul.writeto(os.path.join(output_info['output_dir'],new_name), overwrite=True)                                   
                new_name = os.path.join(pathout, new_name)
                print("Saving in {}...".format(new_name))
                hdul.writeto(new_name,overwrite=True)                
                print("Saved.\n")
                #print('not',hdul.info(),new_name)
                del new_name, hdu0, hdu1, hdul, vec_hdu0,vec_hdu1, hdr_hdu0,hdr_hdu1, ifile, ginfo            
        else:
            raise Exception

    elif  output_info['field'].upper()=='FG':
        print('Checking output directory...')
        if output_info['output_dir']==None:
            file_verification(os.getcwd(),"","outputs")
            pathout = os.path.join(os.getcwd(),"outputs")
        else:
            pathout = output_info['output_dir']
        if output_info['clear_output_dir']:
            cmd = "rm -rf {}/*".format(pathout)
            os.popen(cmd).readlines() 
        print('Output directory ready.')
        ginfo = {'NSIDE' :        output_info['NSIDE'],
              'frequency':{'min': output_info['frequency']['min'],
                           'max': output_info['frequency']['max'],
                        'nbands': output_info['frequency']['nbands']},
            'realization':output_info['realization']}
        
        vec = nu_bins_vector(numin_=ginfo['frequency']['min'], 
                                   numax_=ginfo['frequency']['max'], 
                                   nbands_=ginfo['frequency']['nbands'])
        print('Building HDU data and headers...')
        ginfo['frequency']['nu']= vec['nu']
        try:
            new_name = return_new_FITSfilename(ginfo=ginfo, output_info=output_info, add_info=output_info['add_info'])
        except:
            new_name = return_new_FITSfilename(ginfo=ginfo, output_info=output_info, add_info=None)
        vec_hdu0 = ginfo['frequency']['nu']
        hdr_hdu0 = creating_primary_FITSheader(output_info)
        vec_hdu1 = FG_multiresolution(output_info, verbose=False)
        vec_hdu1 = vec_hdu1[output_info['stokes']]
        hdr_hdu1 = creating_multimap_FITSheader(new_name, 
                                                freq_unit  =output_info["freq_unit"], 
                                                stokes_unit=output_info["stokes_unit"], 
                                                fwhm_unit  =output_info["fwhm_unit"])     
        print('HDUs ready.')
        hdu0 = fits.PrimaryHDU(header=hdr_hdu0,data=vec_hdu0)
        hdu1 = fits.ImageHDU(  header=hdr_hdu1,data=vec_hdu1, name="MULTIMAPS")
        hdul = fits.HDUList([hdu0,hdu1])
        print('Saving FITS file...')
        hdul.writeto(os.path.join(output_info['output_dir'],new_name), overwrite=True)
        print('FITS file saved.')
    
    elif  output_info['field'].upper()=='PL':
        print('Checking output directory...')
        if output_info['output_dir']==None:
            file_verification(os.getcwd(),"","outputs")
            pathout = os.path.join(os.getcwd(),"outputs")
        else:
            pathout = output_info['output_dir']
        if output_info['clear_output_dir']:
            cmd = "rm -rf {}/*".format(pathout)
            os.popen(cmd).readlines() 
        print('Output directory ready.')
        ginfo = {'NSIDE' :        output_info['NSIDE'],
              'frequency':{'min': output_info['frequency']['min'],
                           'max': output_info['frequency']['max'],
                        'nbands': output_info['frequency']['nbands']},
            'realization':output_info['realization']}
        
        vec = nu_bins_vector(numin_=ginfo['frequency']['min'], 
                                   numax_=ginfo['frequency']['max'], 
                                   nbands_=ginfo['frequency']['nbands'])
        print('Building HDU data and headers...')
        ginfo['frequency']['nu']= vec['nu']  
        if output_info['nameFG_Q'] == None:
            FG_Q = 0
        else:
            FG_Q = dcopy(getmap(output_info['pathFG_Q'], output_info['nameFG_Q'], False, 1))
        if output_info['nameFG_U'] == None:
            FG_U = 0
        else:
            FG_U = dcopy(getmap(output_info['pathFG_U'], output_info['nameFG_U'], False, 1))
        if output_info['nameFG_V'] == None:
            FG_V = 0
        else:
            FG_V = dcopy(getmap(output_info['pathFG_V'], output_info['nameFG_V'], False, 1))
        
        new_name = return_new_FITSfilename(ginfo=ginfo, output_info=output_info, add_info=output_info["add_info"])
        vec_hdu0 = ginfo['frequency']['nu']
        hdr_hdu0 = creating_primary_FITSheader(output_info)
        FGpl     = polarization_leakage(mQ=FG_Q, epsQ=output_info['eQ'],
                                        mU=FG_U, epsU=output_info['eQ'], 
                                        mV=FG_V, epsV=output_info['eQ'])
        vec_hdu1  = FGpl
        hdr_hdu1  = creating_multimap_FITSheader(new_name, 
                                                       freq_unit  =output_info["freq_unit"], 
                                                       stokes_unit=output_info["stokes_unit"], 
                                                       fwhm_unit  =output_info["fwhm_unit"])
        print('HDUs ready.')
        try:
            if int(output_info['savedata_hdu'])==0:            
                hdu0 = fits.PrimaryHDU(header=hdr_hdu1,data=vec_hdu1)
                hdu1 = fits.ImageHDU(  header=hdr_hdu0,data=vec_hdu0, name="REDSHIFTS")
                hdul = fits.HDUList([hdu0,hdu1])
            elif int(output_info['savedata_hdu'])==1:            
                hdu0 = fits.PrimaryHDU(header=hdr_hdu0,data=vec_hdu0)
                hdu1 = fits.ImageHDU(  header=hdr_hdu1,data=vec_hdu1, name="MULTIMAPS")
                hdul = fits.HDUList([hdu0,hdu1])                
            else:
                raise Expection("Savedata HDU must be 0 or 1. Current: {}".format(output_info['savedata_hdu']))
        except:
            hdu0 = fits.PrimaryHDU(header=hdr_hdu0,data=vec_hdu0)
            hdu1 = fits.ImageHDU(  header=hdr_hdu1,data=vec_hdu1, name="MULTIMAPS")
            hdul = fits.HDUList([hdu0,hdu1])                            
        hdul.writeto(os.path.join(output_info['output_dir'],new_name), overwrite=True)                            
        print('Saving FITS file...')
        hdul.writeto(os.path.join(output_info['output_dir'],new_name), overwrite=True)       
        print('FITS file saved.')
    else:
        raise Exception('This field is not acceptable')
        
    print("Completed.")

    
    ###########################################
    ###########################################
    ###
    ##
    #    mask
    ##
    ###
    ###########################################
    ###########################################
    
    
def getforegrounds(dirpath=None):
    from copy import deepcopy as dcopy
    namesFG = os.listdir(dirpath)
    for i,iname in enumerate(namesFG):
        mapsiFG = getmap(dirpath,iname,False)
        iFG = iname.split("_")[1]        
        if i==0:
            maps_iFG   = {iFG:mapsiFG}
            maps_iFG["total"] = dcopy(mapsiFG)
        else:
            maps_iFG[iFG] = mapsiFG
            maps_iFG["total"] += dcopy(mapsiFG)  
    return maps_iFG
    
def radec_mask(ra_all=None, dec_all=None, ra_min=None, ra_max=None, dec_min=None, dec_max=None):
    pix_dec = np.where((dec_all<dec_max)*(dec_all>dec_min))[0]
    if (ra_min<0) and (ra_max<0):
        ramax  = np.amax((ra_min, ra_max))+360
        ramin  = np.amin((ra_min, ra_max))+360
        pix_ra = np.where((ra_all<=ramax)*(ra_all>=ramin))[0]    
    elif (ra_min*ra_max)<0:    
        ramax = np.amax((ra_min, ra_max))
        ramin = np.amin((ra_min, ra_max))+360
        pix_ra1 = np.where((ra_all<ramax))[0]
        pix_ra2 = np.where((ra_all<=360)*(ra_all>=ramin))[0]
        pix_ra  = np.union1d(pix_ra1, pix_ra2)
    else:
        ramax = np.amax((ra_min, ra_max))
        ramin = np.amin((ra_min, ra_max))
        pix_ra = np.where((ra_all<=ramax)*(ra_all>=ramin))[0] 
    return np.intersect1d(pix_dec, pix_ra)    

def building_mask(output_info=None):
    from copy import deepcopy as dcopy
    nside      = output_info['nside']
    npix       = 12*nside**2 
    ra,dec     = hp.pix2ang(nside, lonlat=True,ipix=np.arange(npix))
    #pixs       = np.where((dec<output_info['dec']['max'])*(dec>output_info['dec']['min']))[0]
    pixs       = radec_mask(ra_all=ra, dec_all=dec, ra_min=output_info['ra']['min'], ra_max=output_info['ra']['max'], dec_min=output_info['dec']['min'], dec_max=output_info['dec']['max'])
    mask       = np.zeros(npix)
    mask[pixs] = 1
    # second mask -- galactic mask
    if output_info['foreground_cut']:
        #FG         = getmaps_foregrounds(dir_foregrounds=dir_foregrounds)
        #fg         = dcopy(FG['total'][0]) #more intense
        FG = getmap(dirpath_ =output_info['pathdir_FG'], 
                          filename_=output_info['nameFG'], 
                          healpix_readingformat=False, hdu=1)
        fg       = dcopy(FG[0]) #more intense
        npix     = hp.get_map_size(fg)
        n2remove = int(npix*output_info['perc_foreground_cut'])
        sort_fg  = np.flip(np.sort(fg))[:n2remove]
        idex     = np.array([])
        for isort in sort_fg:
            ide  = np.where(fg==isort)[0]
            idex = np.hstack((idex,ide))
            idex = np.int64(idex)
        idex = np.int64(idex)
        idex = idex[:n2remove]
        # combining first and second masks and apodizing the boundary
        mask       = dcopy(mask)
        mask[idex] = 0
    if output_info['apodization']:
        import pymaster as nmt
        mask    = nmt.mask_apodization(mask, output_info['apod_scale'], apotype=output_info['apod_type'])
    return mask

def return_new_FITSfilename_mask(output_info=None):
    #'type':'mask'       #the formatting will follow the mask new formatting name scheme
    #'add_info': None    #None or string name: If None, it'll add in the end nothing. If Str name, the name will be added before '.fits' term
    #'nside': 256        #healpix map resolution  
    #'dec':{'min': -25.48, 'max': -10.17} #declination region does not masked 
    #'galactic_cut':True #If there is some galactic region cutted
    #'foreground_dirpath': "/media/BINGODATA1/ComponentSeparation/MAPS/PAPER/PSM_Components"  #the path to the directory where there are the foreground emission maps 
    #'perc_galactic_cut': 0.2, #Percentage of Foregorund intensity emission to be cutted # It'll be got the first channel because it is higher in intensity medium value
    #'apodization':True,       #True or False: True to apply NaMaster apodization
    #'apod_scale':5,           #Apodization degrees  on the boundary to be applied by NaMaster
    #'apod_type':'C2',         #Type of apodization
    #output_info['apod_scale_unit'] = 'deg'
    if output_info['type'].lower()=='mask':
        output_info['type'] = output_info['type'].lower()
        filename  = "_".join((output_info['type'], str(output_info['nside']) ))
        if output_info['apodization']:
            apod_part = "".join((str(output_info['apod_scale']), 'deg', output_info['apod_type'], 'apod'))
            filename = "_".join((filename, apod_part))
        if output_info['foreground_cut']:
            gal_part = "".join(( str(int(100*output_info['perc_foreground_cut'])),'fgcut' ))
            filename = "_".join((filename, gal_part))
        if (output_info['add_info']!=None)*(type(output_info['add_info'])==str):    
            filename = "_".join((filename, output_info['add_info']))
        return ".".join((filename,'fits'))
    else:
        raise ValueError('output_file type does not accept')
        
def creating_primary_FITSheader_mask():
    hdr = fits.Header()
    hdr['COMMENT']="THIS FILE HAS A SKY MASK."
    return hdr

def creating_mask_FITSheader(output_info=None):
    import datetime
    hdr = fits.Header()
    hdr['COMMENT'] = "THIS IS A NEW NAME FORMAT FOR THE BINGO PROJECT FITS FILES."
    hdr['COMMENT'] = "ALL THESE FILES ARE USED FOR COMPONENT SEPARATION CODE CALLED --CHISEL--"
    hdr['COMMENT'] = "https://github.com/multinverse/CHISEL---signals-from-LIM"
    hdr['COMMENT'] = "THIS IS THE NEW FORMAT FOR MASKS."
    hdr['COMMENT'] = "CREATED IN: {}".format(str(datetime.datetime.now()).split(".")[0])
    for i,ikey in enumerate(output_info.keys()):
        ikey = ikey.lower()
        if ikey=="project":
            hdr[ikey] = output_info[ikey]
            #hdr.comments['field'] = "{}".format("...")
        if ikey=="nside":
            hdr[ikey] = output_info[ikey]
            #hdr.comments["stokes"] = "{}".format(stokes_unit)       
        if ikey=="dec":
            hdr["_".join((ikey,"min"))]=output_info[ikey]['min']
            hdr["_".join((ikey,"max"))]=output_info[ikey]['max']
            hdr.comments["_".join((ikey,"min"))]="in {}".format('deg')
            hdr.comments["_".join((ikey,"max"))]="in {}".format('deg')
        if ikey=="foreground_cut":
            hdr['fgcut'] = output_info[ikey]
            hdr.comments['fgcut']  = "remove part of foreground emission"
            if output_info[ikey]:
                hdr['perc_cut'] = 100*output_info['perc_foreground_cut']
                hdr.comments['perc_cut'] = "% higher foreground removed "
        if ikey=="apodization":
            hdr['apod'] = output_info[ikey]
            hdr.comments['apod'] = 'apply apodization'
            if output_info[ikey]:
                hdr['ap_scale'] = output_info['apod_scale']
                hdr.comments['ap_scale'] = 'apodization scale in deg'
                hdr['ap_type'] = output_info['apod_type']
                hdr.comments['ap_type'] = 'apodization type from NaMASTER'           
    return hdr

def save_FITS_newformat_mask(output_info=None):
    mname    = return_new_FITSfilename_mask(output_info)
    hdr_hdu0 = creating_primary_FITSheader_mask()
    vec_hdu1 = building_mask(output_info) #mask
    hdr_hdu1 = creating_mask_FITSheader(output_info) #Header to be used 

    hdu0 = fits.PrimaryHDU(header=hdr_hdu0)
    hdu1 = fits.ImageHDU(header=hdr_hdu1,data=vec_hdu1, name="MASK")
    hdul = fits.HDUList([hdu0,hdu1])
    if output_info['output_dir']==None:
        filename = os.path.join(os.getcwd()              , mname)
    else:
        filename = os.path.join(output_info['output_dir'], mname)
    if output_info['output_dir']==None:
        file_verification(os.getcwd(),"","outputs")
        pathout = os.path.join(os.getcwd(),"outputs")
    else:
        pathout = output_info['output_dir']   
    print('Filename {}'.format(mname))
    print('Saving in {}'.format(pathout))
    hdul.writeto(filename,overwrite=True)
    print('Saved.\n')


    
    ###########################################
    ###########################################
    ###
    ##
    #    Noise
    ##
    ###
    ###########################################
    ###########################################

def building_realization_strnum(num=100):
    reals = np.array([])
    for i in np.arange(num)+1:
        str_i = str(i)
        len_i = len(str_i)
        num_0 = 4-len_i
        reals = np.append(reals, "".join((num_0*'0',str_i)))
    return reals

def masking_map_by_hitmap(map_=None, hitmap_=None):
    from copy import deepcopy as dcopy
    omap = dcopy(map_)
    if len(map_)>1:
        hitmap_ = hp.ud_grade(hitmap_, hp.get_nside(omap[0]))
        nmap = np.zeros_like(omap)        
        for i,imap in enumerate(omap):
            imap[hitmap_==0]=0
            nmap[i] = imap
        return nmap
    else:
        hitmap_ = hp.ud_grade(hitmap_, hp.get_nside(omap))
        omap[hitmap_==0]=0
        return omap
    
def masking_map_by_pixels(map_=None, pixels_=None):
    from copy import deepcopy as dcopy
    omap = dcopy(map_)
    if len(map_)>1:
        nmap = np.zeros_like(omap)        
        for i,imap in enumerate(omap):
            imap[pixels]=0
            nmap[i] = imap
        return nmap
    else:
        omap[pixels]=0
        return omap    
        
def masking_map_by_declination(map_=None, dec_min=None, dec_max=None):
    npix    = hp.get_map_size(map_)
    nside   = hp.get_nside(map_)
    _,dec   = hp.pix2ang(nside, lonlat=True,ipix=np.arange(npix))
    pixels_ = np.where((dec<dec_max)*(dec>dec_min))[0]  
    return masking_map_by_pixels(map_, pixels_)

def masking_map(map_=None, hitmap_=None, pixels_=None,  dec_min=None, dec_max=None):
    try:
        return masking_map_by_hitmap(map_, hitmap_)
    except: pass
    try:
        return masking_map_by_pixels(map_, pixels_)
    except: pass
    try:
        return masking_map_by_declination(map_, dec_min, dec_max)
    except:
        raise Exception



def save_FITS_newformat_wn(output_info=None, ginfo=None, sigmaN=None, hitmap=None):
    from copy import deepcopy as dcopy
    if type(output_info['realization'])==list or type(output_info['realization'])==np.ndarray:
        pass
    elif type(output_info['realization'])==int:
        output_info['realization'] = dcopy(building_realization_strnum(output_info['realization']))
    else:
        raise Expection      
    nside = output_info['nside']
    nch = output_info['frequency']['nbands']
    
    vec = nu_bins_vector(numin_=ginfo['frequency']['min'], 
                                numax_=ginfo['frequency']['max'], 
                                nbands_=ginfo['frequency']['nbands'])
    print('Building HDU data and headers...')
    ginfo['frequency']['nu']= vec['nu'] 
    
    for i,ireal in enumerate(output_info['realization']):
        ginfo['realization']=ireal
        print("Realization {}".format(ireal))
        np.random.seed(int(ireal))
        wn       = dcopy(np.random.normal(scale = sigmaN, size=(nch,12*nside**2)))# * sigmaN.unit
        try:
            wn   = masking_map(map_=wn, hitmap_=hitmap)
        except:
            pass
        new_name = return_new_FITSfilename(ginfo, output_info)
        vec_hdu0 = ginfo['frequency']['nu']
        hdr_hdu0 = creating_primary_FITSheader(output_info)
        vec_hdu1 = wn
        hdr_hdu1 = creating_multimap_FITSheader(new_name, 
                                                      freq_unit  =output_info["freq_unit"], 
                                                      stokes_unit=output_info["stokes_unit"], 
                                                      fwhm_unit  =output_info["fwhm_unit"]) #Header to be used 
        try:
            if int(output_info['savedata_hdu'])==0:            
                hdu0 = fits.PrimaryHDU(header=hdr_hdu1,data=vec_hdu1)
                hdul = fits.HDUList([hdu0])
            elif int(output_info['savedata_hdu'])==1: 
                hdu0 = fits.PrimaryHDU(header=hdr_hdu0,data=vec_hdu0)
                hdu1 = fits.ImageHDU(  header=hdr_hdu1,data=vec_hdu1, name="MULTIMAPS")
                hdul = fits.HDUList([hdu0,hdu1])                
            else:
                raise Expection("Savedata HDU must be 0 or 1. Current: {}".format(output_info['savedata_hdu']))
        except:
            hdu0 = fits.PrimaryHDU(header=hdr_hdu0,data=vec_hdu0)
            hdu1 = fits.ImageHDU(  header=hdr_hdu1,data=vec_hdu1, name="MULTIMAPS")
            hdul = fits.HDUList([hdu0,hdu1])                            
        
        if output_info['output_dir']==None:
            filename = os.path.join(os.getcwd()              , new_name)
        else:
            filename = os.path.join(output_info['output_dir'], new_name)
        if output_info['output_dir']==None:
            file_verification(os.getcwd(),"","outputs")
            pathout = os.path.join(os.getcwd(),"outputs")
        else:
            pathout = output_info['output_dir']    
        print("Filename: {}".format(new_name))
        print("Saving in {},,,".format(filename))
        hdul.writeto(filename,overwrite=True)
        print("Saved.\n")
        
def get_noise_level(sigma_info=None,verbose=False):
    t_unit = 24*60*60
    tsur   = sigma_info['tsur']
    tsur   *= 365*t_unit #yr --> sec
    tsur   *= sigma_info['dcycle']
    sigma_info['tsur']     = tsur
    sigma_info['nu_max']   *= 1e6 #MHz --> Hz
    sigma_info['nu_min']   *= 1e6    
    sigma_info['bandwidth'] = (sigma_info['nu_max']-sigma_info['nu_min'])/sigma_info['nch']
    
    tpix = tpix_func(tsur = sigma_info['tsur'] , Osur  = sigma_info['Osur'] , 
                     Obeam= sigma_info['Obeam'], nbeams= sigma_info['nbeams'],
                     verbose=verbose)
    sigma_info['tpix']=tpix
    #Snoise = sigma_noise(tpix= sigma_info['tpix'], bandwidth= sigma_info['bandwidth'], 
    #                     Tsys= sigma_info['Tsys'], K= sigma_info['K'],
    #                     verbose=verbose)
    #sigma_info['Snoise']=Snoise
    #sigma_info['sigmaN'] = sigma_info['Spix']*1e6
    Spix = sigma_pix_func(tpix_=sigma_info['tpix'], bandwidth_=sigma_info['bandwidth'] , K_=sigma_info['K'], 
                          Tsys_=sigma_info['Tsys'], verbose=verbose)#unit=Kelvin
    if sigma_info['output_unit']=='mK':
        Spix*=1e3
    elif sigma_info['output_unit']=='muK':
        Spix*=1e6
    else:
        Spix*=1e0    
    sigma_info['sigma_pix']=Spix
    sigma_info['sigmaN']=Spix
    return sigma_info #noise/pix   

def tpix_func(tsur=None, Osur=None, Obeam=None, nbeams=None, verbose=False):
    tpix_ = (tsur/(Osur/Obeam))*nbeams
    if verbose: print("tpix: {:.2f} hour/pix".format(tpix_/3600),"\n")    
    return tpix_
#def sigma_noise(tpix=None, bandwidth=None, Tsys=None, K=None, verbose=False):
    #Snoise_ = K*Tsys/np.sqrt(tpix*bandwidth)
    #if verbose: print("sigmaN: {:.2f} mK".format(Snoise_*1e6),"\n")
#    return Snoise_
def sigma_pix_func(tpix_=None, bandwidth_=None, K_=None, Tsys_=None, output_unit='K', verbose=False): 
    #[tpix]=s
    #[bandwidth]=Hz    
    #[Tsys]=K
    spix_ = K_*Tsys_/np.sqrt(tpix_*bandwidth_)
    if output_unit=='mK':
        spix_*=1e3
    elif output_unit=='muK':
        spix_*=1e6
    else:
        spix_*=1e0
    if verbose: print("sigma pix: {:.6f} mK/pix".format(spix_),"\n")
    return spix_

def Nl_hi_theorectical(sigma_pix=None, nside_=None, output_unit='K'):#[sigma_pix]=K
    if output_unit=='K':
        Nl_= (4/Npix)*np.pi*(Spix*1e0)**2#muK2
    elif output_unit=='mK':
        Nl_= (4/Npix)*np.pi*(Spix*1e3)**2#muK2
    elif output_unit=='muK':
        Nl_= (4/Npix)*np.pi*(Spix*1e6)**2#muK2        
    else:
        Nl_= (4/Npix)*np.pi*(Spix*1e0)**2#muK2
    return Nl_
#########################
#### FG
#########################

def pyPSM_multiresolution(nu_min=None,nu_max=None,nbands=None,sky_model=None, verbose=True):
    import pysm3
    import pysm3.units as u
    stepf = (nu_max-nu_min)/nbands
    freq = np.around(np.arange(nu_min,nu_max+stepf,stepf),decimals=4) * u.MHz  #<---- look at this unit. Pay attetion on it
    for i in range(nbands):
        if verbose: print(i+1,freq[i:i+2])
        iFG = sky_model.get_emission(freq[i:i+2]).value #muK
        iFG = hp.Rotator(coord=['G','C']).rotate_map_pixel(iFG)
        #iFG = hp.Rotator(coord=['G','C']).rotate_map_alms(m=iFG,use_pixel_weights=False)  <---it's accurater
        qFG = iFG[1]/1e3  #muK --> mK
        uFG = iFG[2]/1e3
        iFG = iFG[0]/1e3
        if not i:
            IFG = iFG
            QFG = qFG
            UFG = uFG
        else:
            IFG = np.vstack((IFG,iFG))
            QFG = np.vstack((QFG,qFG))
            UFG = np.vstack((UFG,uFG))
    return {"I":IFG,"Q":QFG,"U":UFG}     

def FG_multiresolution(output_info=None, verbose=None):
    print('Starting to build FG dataset...')
    from copy import deepcopy as dcopy
    import pysm3
    import pysm3.units as u    
    sky = pysm3.Sky(nside=output_info["NSIDE"], preset_strings=output_info["psm_models"])
    fg  = pyPSM_multiresolution(nu_min=output_info["frequency"]["min"],
                                nu_max=output_info["frequency"]["max"],
                                nbands=output_info["frequency"]["nbands"],
                                sky_model=sky, verbose=verbose)
    fg = dcopy(fg) 
    if output_info["include_rps"]:
        ps = dcopy(getmap(output_info["pathrps"], output_info["filename_rps"], False))
        fg["I"]+=ps
#    if output_info["apply_mask"]:
#        pass
    print('FG dataset ready.')
    return fg

#########################
#### PL
#########################
def polarization_leakage(mQ=None, mU=None, mV=None, epsQ=None, epsU=None, epsV=None):
    return epsQ*mQ+epsU*mU+epsV*mV
def Li_Yang_Gao_coeffs(): #高丽阳2022
    eQ = np.array([0.0, 0.5, 1.0, 2.0, 0.0, 0.5, 1.0, 2.0])/100
    eU = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0])/100
    eV = 10*eQ #idk how ill use it yet
    return {"Q":eQ,"U":eU ,"V":eV}


#########################
#### UD GRADE HI maps
#########################

def ud_grade_HImaps(output_info=None):
    from copy import deepcopy as dcopy
    filenames = get_allfilenames(dirpath_=output_info['pathdir_HI'])
    filenames = np.sort(filenames)
    
    if output_info['realization']=='all':
        pass
    elif type(output_info['realization'])==list:
        ids_       = get_all_realizations_id(filenames, new_format=True)
        nfilenames = np.array([])
        for i,ireal in enumerate(output_info['realization']):
            ind = np.where(ids_==ireal)[0]
            nfilenames = np.append(nfilenames, filenames[ind])
        filenames = dcopy(nfilenames)
    elif type(output_info['realization'])==int:  
        filenames = dcopy(filenames[:output_info['realization']])
    else:
        raise Exception
    print('Checking output directory...')
    if output_info['output_dir']==None:
        file_verification(os.getcwd(),"","outputs")
        pathout = os.path.join(os.getcwd(),"outputs")
    else:
        pathout = output_info['output_dir']
    if output_info['clear_output_dir']:
        cmd = "rm -rf {}/*".format(pathout)
        os.popen(cmd).readlines()     
    print('Output directory ready.\n')    
    print('Changing dimensions...')
    for i,ifilename in enumerate(filenames):
        print('File {} as input:  {}'.format(i+1,ifilename))
        omap = getmap(dirpath_ =output_info['pathdir_HI'], 
                            filename_=ifilename, 
                            healpix_readingformat=False, hdu=1)
        vec_hdu1  = hp.ud_grade(omap, output_info['new_nside'])
        old_nside = np.int16(ifilename.split("_")[2]) #hp.get_nside(imap)
        new_name  = ifilename.replace("_{:d}_".format(old_nside),
                                     "_{:d}_".format(output_info['new_nside']))   
        vec_hdu0 = getmap(dirpath_ =output_info['pathdir_HI'], 
                                filename_=ifilename, 
                                healpix_readingformat=False, hdu=0)
        with fits.open(os.path.join(output_info['pathdir_HI'],ifilename)) as h:
            maps = dcopy(h)
        hdr_hdu0 = maps[0].header
        hdr_hdu1 = maps[1].header
        hdr_hdu1['NSIDE'] = output_info['new_nside']
        try:
            if int(output_info['savedata_hdu'])==0:            
                hdu0 = fits.PrimaryHDU(header=hdr_hdu1,data=vec_hdu1)
                hdu1 = fits.ImageHDU(  header=hdr_hdu0,data=vec_hdu0, name="REDSHIFTS")
                hdul = fits.HDUList([hdu0,hdu1])
            elif int(output_info['savedata_hdu'])==1:            
                hdu0 = fits.PrimaryHDU(header=hdr_hdu0,data=vec_hdu0)
                hdu1 = fits.ImageHDU(  header=hdr_hdu1,data=vec_hdu1, name="MULTIMAPS")
                hdul = fits.HDUList([hdu0,hdu1])                
            else:
                raise Expection("Savedata HDU must be 0 or 1. Current: {}".format(output_info['savedata_hdu']))
        except:
            hdu0 = fits.PrimaryHDU(header=hdr_hdu0,data=vec_hdu0)
            hdu1 = fits.ImageHDU(  header=hdr_hdu1,data=vec_hdu1, name="MULTIMAPS")
            hdul = fits.HDUList([hdu0,hdu1])                            
        hdul.writeto(os.path.join(output_info['output_dir'],new_name), overwrite=True)
        print('File {} as output: {}\n'.format(i+1,new_name))
    print("Completed.") 




#########################
#### OBSERVED
#########################

def beam_observed_maps(output_info=None, ifilename=None, pl_position=0):
    import beam_modelling         as model  
    HImap = getmap(dirpath_ =output_info['pathdir_HI'], filename_=ifilename, 
                             healpix_readingformat=False, hdu=1) 
    imap = dcopy(HImap) 
    #PL before beam convolution
    if pl_position==0:
        for i, ifield in enumerate(['FG', 'PL']):
            if type(output_info['name{}'.format(ifield)])!=type(None):
                jmap =getmap(dirpath_ = output_info['pathdir_{}'.format(ifield)], 
                                    filename_= output_info['name{}'.format(ifield)], 
                                   healpix_readingformat=False, hdu=1)
                imap +=jmap 
                if not i:
                    fields = {'{} maps'.format(ifield): jmap}
                else:
                    fields['{} maps'.format(ifield)]=jmap
            else:
                if not i:
                    fields = {'{} maps'.format(ifield): 0*imap}
                else:
                    fields['{} maps'.format(ifield)]=0*imap       
        #Loop2: channels
        omap = np.zeros((imap.shape[0], hp.nside2npix(output_info['new_nside'])))
        bl   = np.zeros((imap.shape[0], int(output_info['lmax']+1)))
        for i, inu in enumerate(output_info['nu'][:-1]):
            output_info['fwhm'] = dcopy(model.fwhm_modelling(nu=inu, type_= output_info["type_fwhm"], D=output_info["main_dish_diameter"], in_degree=True, if_fwhm_fixed=output_info["if_fwhm_fixed"]))
            bl[i,:]             = model.bl_function(type_=output_info['type_beam'], fwhm=output_info['fwhm'] , lmax=output_info['lmax'], theta_=output_info['theta_range'], input_unit="degree")
            omap[i]             = dcopy(hp.ud_grade(imap[i], output_info['new_nside']))
        imap = dcopy(omap)
        del omap
        
        if output_info['apply_beam']:
            imap_mod = np.zeros_like(imap)
            for i,imap_ in enumerate(imap):
                alm_imap     = hp.map2alm(imap_,lmax=output_info['lmax'])
                alm_imap_mod = hp.almxfl(alm_imap, bl[i])
                imap_mod[i]  = hp.alm2map(alm_imap_mod,nside=output_info['new_nside'],lmax=output_info['lmax'], pol=False)
        else:
            imap_mod = np.zeros_like(imap)       
        #del HImap, jmap, alm_imap, alm_imap_mod, output_info['fwhm'], bl, output_info['theta_range'], output_info['nu'], output_info['l']
        return {"unconvolved map": imap,
                "convolved map"  : imap_mod,
                "HI maps": HImap,
                "FG maps": fields['FG maps'],           
                "PL maps": fields['PL maps'],}
    #PL after beam convolution
    if pl_position==1:
        if type(output_info['nameFG'])!=type(None):
            jmap = getmap(dirpath_ = output_info['pathdir_FG'], 
                          filename_= output_info['nameFG'], 
                          healpix_readingformat=False, hdu=1)
            imap+=jmap
            fields = {'FG maps':jmap}   
        #Loop2: channels
        omap = np.zeros((imap.shape[0], hp.nside2npix(output_info['new_nside'])))
        bl   = np.zeros((imap.shape[0], int(output_info['lmax']+1)))
        for i, inu in enumerate(output_info['nu'][:-1]):
            output_info['fwhm'] = dcopy(model.fwhm_modelling(nu=inu, type_= output_info["type_fwhm"], D=output_info["main_dish_diameter"], in_degree=True, if_fwhm_fixed=output_info["if_fwhm_fixed"]))
            bl[i,:]             = model.bl_function(type_=output_info['type_beam'], fwhm=output_info['fwhm'] , lmax=output_info['lmax'], theta_=output_info['theta_range'], input_unit="degree")
            omap[i]             = dcopy(hp.ud_grade(imap[i], output_info['new_nside']))
        imap = dcopy(omap)
        del omap
        
        if output_info['apply_beam']:
            imap_mod = np.zeros_like(imap)            
            for i,imap_ in enumerate(imap):
                alm_imap     = hp.map2alm(imap_,lmax=output_info['lmax'])
                alm_imap_mod = hp.almxfl(alm_imap, bl[i])
                imap_mod[i]  = hp.alm2map(alm_imap_mod,nside=output_info['new_nside'],lmax=output_info['lmax'], pol=False)                
        else:
            imap_mod = np.zeros_like(imap)

        if type(output_info['namePL'])!=type(None):
            jmap = getmap(dirpath_ = output_info['pathdir_PL'], 
                          filename_= output_info['namePL'], 
                          healpix_readingformat=False, hdu=1)
            imap+=hp.ud_grade(jmap, output_info['new_nside'])
            fields['PL maps'] = jmap
            if output_info['apply_beam']:
                imap_mod += hp.ud_grade(jmap, output_info['new_nside'])            
        #del HImap, jmap, alm_imap, alm_imap_mod, output_info['fwhm'], bl, output_info['theta_range'], output_info['nu'], output_info['l']
        return {"unconvolved map": imap,
                "convolved map"  : imap_mod,
                "HI maps": HImap,
                "FG maps": fields['FG maps'],           
                "PL maps": fields['PL maps'],}


def sky_observation(output_info=None, ifilename=None, pl_position=0):
    beamed_map   = beam_observed_maps(output_info, ifilename, pl_position)
    if output_info["apply_beam"]:
        beamed_map['observed map']=dcopy(beamed_map['convolved map'])
    else:
        beamed_map['observed map']=dcopy(beamed_map['unconvolved map'])
    if type(output_info['pathdir_N'])!=type(None) and type(output_info['pathdir_N'])==str:
        _id_                 = ifilename.split("_L")[1].split(".fits")[0]
        filenamesN_template  = get_allfilenames(dirpath_=output_info['pathdir_N'])[0]
        part                 = filenamesN_template.split("_L")[1]
        Nfilename            = filenamesN_template.replace(part, _id_+".fits")
        output_info['nameN'] = Nfilename        
        Nmap = getmap(dirpath_ =output_info['pathdir_N'], 
                      filename_=Nfilename, 
                      healpix_readingformat=False, hdu=1)
        Nmap = hp.ud_grade(Nmap, output_info['new_nside'])
        #imap_mod+=Nmap
        beamed_map['observed map']+=Nmap
        beamed_map['N maps'] = Nmap
    else:
        beamed_map['N maps'] = 0*beamed_map['observed map'][0]        
    if output_info['apply_mask']:
        output_info['coverage'] = 'partial'
        mask = getmap(dirpath_ = output_info['pathdir_M'], 
                            filename_= output_info['nameM'], 
                            healpix_readingformat=False, hdu=1)
        beamed_map['observed map']*=mask
        beamed_map['mask'] = mask
    else:
        beamed_map['mask'] = 0*beamed_map['observed map'][0]
    return beamed_map


#####
#TEST
#####
def creating_observed_FITSheader(output_info=None, freq_unit="Mhz", stokes_unit="mk", fwhm_unit="arcmin"):
    import datetime
    hdr = fits.Header()
    hdr['COMMENT']  = "THIS IS A NEW NAME FORMAT FOR THE BINGO PROJECT FITS FILES."
    hdr['COMMENT']  = "THESE FILES ARE USED FOR COMPONENT SEPARATION CODE CALLED --CHISEL--"
    hdr['COMMENT']  = "https://github.com/multinverse/CHISEL---signals-from-LIM"
    hdr['COMMENT']  = "ITS FORMAT ASSUME FLASK MAPS AS INPUT"
    hdr['COMMENT']  = "CREATED IN: {}".format(str(datetime.datetime.now()).split(".")[0])    
    hdr['field']    = output_info['field']
    hdr["stokes"]   = output_info["stokes"]
    hdr.comments["stokes"]   = "{}".format(stokes_unit) 
    hdr["coverage"] = output_info["coverage"]
    hdr.comments["coverage"] = "{full- or partial-sky}"
    hdr["FREQ_MIN"]= output_info['frequency']['min']
    hdr["FREQ_MAX"]= output_info['frequency']['max']
    hdr.comments["FREQ_MIN"]="in {}".format(freq_unit)
    hdr.comments["FREQ_MAX"]="in {}".format(freq_unit)
    if output_info['apply_beam']:
        hdr['beam'] = output_info['type_beam']
        hdr['fwhm'] = output_info['type_fwhm']
        hdr.comments['fwhm'] = fwhm_unit
        hdr['D'] = output_info['main_dish_diameter']
        hdr.comments["D"] = 'Main dish diameter'
    else:
        hdr["beam"] = 'NONE'
    if type(output_info["namePL"])==str:
        hdr['PL']   = output_info["namePL"]
        hdr.comments['PL'] = 'Polarization Leakage'
    else:
        hdr['PL'] = 'NONE'
    if type(output_info["nameFG"])==str:
        hdr['FG'] = output_info["nameFG"]
        hdr.comments['FG'] = 'Foreground'
    else:
        hdr['FG']   = 'NONE'   
    if output_info['apply_mask']:
        hdr['MASK'] = output_info['nameM']
    else:
        hdr['MASK'] = 'NONE'
    if type(output_info['pathdir_N'])==str:    
        hdr['NOISE'] = output_info['nameN']
    else:
        hdr['NOISE'] = 'NONE'        
    hdr['REALIZATION']= output_info['irealization'] 
    return hdr 


def save_FITS_newformat_observed(output_info=None):   
    from copy import deepcopy as dcopy
    print('Checking output directory...')
    if output_info['output_dir']==None:
        file_verification(os.getcwd(),"","outputs")
        pathout = os.path.join(os.getcwd(),"outputs")
    else:
        pathout = output_info['output_dir']
        #print(pathout.split('/'),"/".join(pathout.split('/')[:-1]), pathout.split('/')[-1])
        file_verification("/".join(pathout.split('/')[:-1]),"",pathout.split('/')[-1])
    if output_info['clear_output_dir']:
        cmd = "rm -rf {}/*".format(pathout)
        os.popen(cmd).readlines()     
    print('Output directory ready.\n')    
    print('Changing dimensions...')  
    filenames = np.sort(get_allfilenames(dirpath_=output_info['pathdir_HI']))    
    if (type(output_info['realization'])==list or type(output_info['realization'])==np.ndarray):
        filenames = np.array(['_'.join((filenames[0].split('_L')[0],'L' +iname+'.fits')) for iname in output_info['realization']])  
    elif type(output_info['realization'])==str and output_info['realization']!='all':
        filenames = np.array(['_'.join((filenames[0].split('_L')[0],'L' +output_info['realization']+'.fits'))])
    elif type(output_info['realization'])==int and output_info['realization']<=len(filenames):
        filenames = filenames[:int(output_info['realization'])]
    elif output_info['realization']=='all':
        pass #filenames=filenames
    else:
        pass #filenames=filenames
    _ids_     = get_all_realizations_id(filenames_=filenames, new_format=True)
    #print("Realizations to be saved: {}".format(filenames))
    print("Realizations to be saved: {}\n\n".format(_ids_))
    #sys.exit(0)
    with fits.open(os.path.join(output_info['pathdir_HI'],filenames[0])) as h:
            output_info['nu'] = dcopy(h[0].data)
    #if there is no new nside just keep the same getting from the first HI name
    if not type(output_info['new_nside'])==int:
        output_info['new_nside'] = int(filenames[0].split("_")[2])
    #LMAX
    if (type(output_info['lmax']) == float or type(output_info['lmax']) == int):
        output_info['lmax'] = np.asarray(output_info['lmax'])
    else:
        output_info['lmax']  = 3*output_info['new_nside']
    #THETA
    if (type(output_info['theta_range'])==np.ndarray or type(output_info['theta_range'])==range or type(output_info['theta_range'])==list):
        output_info['theta_range'] = np.asarray(output_info['theta_range'])
    else:
        output_info['theta_range'] = np.arange(0,10,0.01)
    #
    output_info['l']     = np.arange(output_info['lmax']+1)

    pl_position=output_info['pl_position']
    ###############
    for i,ifilename in enumerate(filenames):
        print("Building maps of L{}...".format(_ids_[i]))
        imap = sky_observation(output_info, ifilename, pl_position) 
        output_info['stokes'] = ifilename.split('_')[1]
        output_info['field' ] = dcopy(output_info['first_name'])
        with fits.open(os.path.join(output_info['pathdir_HI'],ifilename)) as h:
            hdu = dcopy(h)

        output_info['frequency'] = {'min': hdu[1].header['FREQ_MIN'],
                                    'max': hdu[1].header['FREQ_MAX'],
                                 'nbands': hdu[0].header['NBANDS']}
        if output_info['apply_mask']:
            output_info['coverage'] = 'partial'
        elif len(ifilename.split('partial'))>1:
            output_info['coverage'] = 'partial'
        else:
            output_info['coverage'] = 'full'
        output_info['irealization']  = ifilename.split('_L')[1].split('.fits')[0]
        print("Saving the file...")
        hdr_hdu0 = hdu[0].header
        vec_hdu0 = hdu[0].data
        hdr_hdu1 = creating_observed_FITSheader(output_info)
        vec_hdu1 = imap['observed map']
        output_info['beam_resolution'] = None if not output_info['apply_beam'] else str(int(output_info['fwhm']*60))
        output_info['beam_model']      = None if not output_info['apply_beam'] else output_info['type_fwhm']+get_beam_model_inv(output_info['type_beam'])
            
        new_name = new_formattingnames(F=output_info['field'], 
                              S=output_info['stokes'],
                              NSIDE_=output_info['new_nside'], 
                              freq_min=str(int(output_info['frequency']['min'])), freq_max=str(int(output_info['frequency']['max'])),freq_unit="mhz",
                              Num_=str(int(output_info['frequency']['nbands'])), 
                              C=output_info['coverage'], 
                              Bres=output_info['beam_resolution'], Bmodel=output_info['beam_resolution'],
                              A=output_info['add_info'], R=output_info['irealization'])
        try:
            if int(output_info['savedata_hdu'])==0:            
                hdu0 = fits.PrimaryHDU(header=hdr_hdu1,data=vec_hdu1)
                hdu1 = fits.ImageHDU(  header=hdr_hdu0,data=vec_hdu0, name="REDSHIFTS")
                hdul = fits.HDUList([hdu0,hdu1])
            elif int(output_info['savedata_hdu'])==1:            
                hdu0 = fits.PrimaryHDU(header=hdr_hdu0,data=vec_hdu0)
                hdu1 = fits.ImageHDU(  header=hdr_hdu1,data=vec_hdu1, name="MULTIMAPS")
                hdul = fits.HDUList([hdu0,hdu1])                
            else:
                raise Expection("Savedata HDU must be 0 or 1. Current: {}".format(output_info['savedata_hdu']))
        except:
            hdu0 = fits.PrimaryHDU(header=hdr_hdu0,data=vec_hdu0)
            hdu1 = fits.ImageHDU(  header=hdr_hdu1,data=vec_hdu1, name="MULTIMAPS")
            hdul = fits.HDUList([hdu0,hdu1])                            
        hdul.writeto(os.path.join(output_info['output_dir'],new_name), overwrite=True)                    
        new_name = os.path.join(output_info['output_dir'], new_name)
        print("Saving in {}".format(new_name))
    print('Complete.')
