import time, os, sys
import numpy as np
from copy import deepcopy as dcopy
#sys.path.insert(1, '/home/amarins/Programmation/cross/scripts')
import cross_functions_theory      as cxft
import cross_functions_simulations as cxfs
import handling_data as hdata
import healpy as hp
#sys.path.insert(1, '/home/amarins/Programmation/chisel/scripts')    
#import Extension4BINGO   as cs
#####################################################################################################
#####################################################################################################
#### CREATING CONVERGENCE-field SIMULATION
def suppression_fact(dict_rfields=None, map_mask_=None, lmax=None,
                     field_sflag='HI' , field_rflag='extHI', ns_flag='ns2', coverage_flag='fullsky'):
    
    ii=1 #assuming L0=first realzation
    _nsims_,_npix_ = dict_rfields['HI']['maps']['matrix'].shape
    #if map_mask_==None: map_mask_=np.ones_like(_npix_)
    if field_sflag=='HI':
        for i in range(_nsims_):
            if not i: 
                _scl_ = hp.anafast(dict_rfields[field_sflag]['maps']['matrix'][i]*map_mask_                      ,pol=False,lmax=lmax)
                _rcl_ = hp.anafast(dict_rfields[coverage_flag][ns_flag][field_rflag]['maps']['matrix'][i]*map_mask_,pol=False,lmax=lmax)
            else:
                _scl_ = np.vstack(( _scl_,hp.anafast(dict_rfields[field_sflag]['maps']['matrix'][i]*map_mask_, pol=False,lmax=lmax) ))
                _rcl_ = np.vstack(( _rcl_,hp.anafast(dict_rfields[coverage_flag][ns_flag][field_rflag]['maps']['matrix'][i]*map_mask_,pol=False,lmax=lmax) ))
        return np.mean(np.sqrt(_scl_/_rcl_),axis=0)    
    elif field_sflag=='FG':
        _scl_ = hp.anafast(dict_rfields[field_sflag]['maps']['mean']*map_mask_,pol=False,lmax=lmax)
        for i in range(_nsims_):
            if not i: 
                _rcl_ = hp.anafast(dict_rfields[coverage_flag][ns_flag][field_rflag]['maps']['matrix'][i]*map_mask_,pol=False,lmax=lmax)
            else:
                _rcl_ = np.vstack(( _rcl_,hp.anafast(dict_rfields[coverage_flag][ns_flag][field_rflag]['maps']['matrix'][i]*map_mask_,pol=False,lmax=lmax) ))
        return np.mean(np.sqrt(_scl_/_rcl_),axis=0)    
        
        
        
def debiasing_correction(dict_rfields=None, dict_sfields=None,
                         nside=256,channel=None,lmax=None, del_l=1,use_dl=False, 
                         ns_flag='ns2', coverage_flag='fullsky', 
                         map_fg_theo=None, map_mask=None, 
                         apply_namaster=True, show_time=True):
    if show_time:itime = time.time()
    ##########################    
    sim_flag   = 'sim0'
    ii         = 1
    mask = dcopy(map_mask)
    fsky = np.mean(mask*mask)
    ind0 = np.where(mask<=0)[0]
        
    #########################
    #Theory--extracted from the simulation
    alm_hi_theo = dict_sfields[sim_flag]['alm_hi_sim'][int(channel+2)]
    alm_hi_theo = np.ascontiguousarray(alm_hi_theo)#[0=l,1=m,...]
    m_hi_theo   = hp.alm2map(alm_hi_theo, lmax=lmax, pol=False , nside=nside)*mask
    cl_hi_theo  = hp.anafast(m_hi_theo  , lmax=lmax, pol=False)    
    #Average over recoveries
    m_hi_rmean  = dict_rfields[coverage_flag][ns_flag]['extHI']['maps']['mean']*mask
    cl_hi_rmean = hp.anafast(m_hi_rmean, lmax=lmax, pol=False)
    #SIM0 biased
    m_hi_bias  = dict_rfields[coverage_flag][ns_flag]['extHI']['maps']['matrix'][int(ii-1)]*mask
    cl_hi_bias = hp.anafast(m_hi_bias, lmax=lmax, pol=False)
    #Debiased
    supp       = suppression_fact(dict_rfields=dict_rfields, map_mask_=mask, lmax=lmax,
                                  field_sflag='HI' , field_rflag='extHI', 
                                  ns_flag=ns_flag, coverage_flag=coverage_flag)
    #####
    alm_R      = dict_rfields[coverage_flag][ns_flag]['extHI']['alms']['matrix'][0] #realizacao assumida Li=0
    alm_hi_rec = hp.almxfl(alm_R, supp)#, mmax=lmax)
    m_hi_rec   = hp.alm2map(alm_hi_rec, pol=False,  lmax=lmax, nside=nside,)*mask#, mmax=lmax)
    cl_hi_rec  = hp.anafast(m_hi_rec, lmax=lmax, pol=False)    
    #########################
    #Theory
    m_fg_theo   = dcopy(map_fg_theo)*mask
    cl_fg_theo  = hp.anafast(m_fg_theo, lmax=lmax, pol=False)        
    #Average over recoveries
    m_fg_rmean  = dict_rfields[coverage_flag][ns_flag]['extFG']['maps']['mean']*mask
    cl_fg_rmean = hp.anafast(m_fg_rmean, lmax=lmax, pol=False)
    #SIM0 biased
    m_fg_bias   = dict_rfields[coverage_flag][ns_flag]['extFG']['maps']['matrix'][int(ii-1)]*mask
    cl_fg_bias  = hp.anafast(m_fg_bias, lmax=lmax, pol=False)
    #Debiased
    suppfg      = suppression_fact(dict_rfields=dict_rfields, map_mask_=mask, lmax=lmax,
                                   field_sflag='FG' , field_rflag='extFG', 
                                   ns_flag=ns_flag, coverage_flag=coverage_flag)
    #####
    alm_fg     = dict_rfields[coverage_flag][ns_flag]['extFG']['alms']['matrix'][0] #realizacao assumida Li=0
    alm_fg_rec = hp.almxfl(alm_fg, suppfg)#, mmax=lmax)
    m_fg_rec   = hp.alm2map(alm_fg_rec, pol=False, nside=nside)*mask
    cl_fg_rec  = hp.anafast(m_fg_rec, lmax=lmax, pol=False)
    ########################
    m_sky_theo = dcopy(m_hi_theo + m_fg_theo)    
    ####  

##########
    if apply_namaster:
        import pymaster as nmt
        bpw  = nmt.NmtBin.from_nside_linear(nside, nlb=del_l)    
        leff = bpw.get_effective_ells()
        cl_hi_rec   = bpw.bin_cell( cl_hi_rec)
        cl_hi_rmean = bpw.bin_cell( cl_hi_rmean)
        cl_hi_bias  = bpw.bin_cell( cl_hi_bias)
        cl_hi_theo  = bpw.bin_cell( cl_hi_theo)
        
        cl_fg_rec   = bpw.bin_cell( cl_fg_rec)
        cl_fg_rmean = bpw.bin_cell( cl_fg_rmean)
        cl_fg_bias  = bpw.bin_cell( cl_fg_bias)
        cl_fg_theo  = bpw.bin_cell( cl_fg_theo)

        #supp   = bpw.bin_cell( supp)
        #suppfg = bpw.bin_cell( suppfg)
    elif (apply_namaster==False) and (int(del_l)>1):
        include_leftover=False
        leff, cl_hi_rec   = cxfs.clsbinned(cl_hi_rec,   del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)
        cl_hi_rmean = cxfs.clsbinned(cl_hi_rmean, del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        cl_hi_bias  = cxfs.clsbinned(cl_hi_bias,  del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        cl_hi_theo  = cxfs.clsbinned(cl_hi_theo,  del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        
        cl_fg_rec   = cxfs.clsbinned(cl_fg_rec,   del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)
        cl_fg_rmean = cxfs.clsbinned(cl_fg_rmean, del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        cl_fg_bias  = cxfs.clsbinned(cl_fg_bias,  del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        cl_fg_theo  = cxfs.clsbinned(cl_fg_theo,  del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]        
        
        #supp   = cxfs.clsbinned(supp,   del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]        
        #suppfg = cxfs.clsbinned(suppfg, del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]        
    else:
        leff = np.arange(cl_fg_rec.size)
    fact = leff*(leff+1)/2/np.pi if use_dl else 1

    var_hi_rec   = (cl_hi_rec   + 0)/(np.sqrt(2*leff+1)*fsky*del_l)
    var_hi_rmean = (cl_hi_rmean + 0)/(np.sqrt(2*leff+1)*fsky*del_l)
    var_hi_bias  = (cl_hi_bias  + 0)/(np.sqrt(2*leff+1)*fsky*del_l)
    var_hi_theo  = (cl_hi_theo  + 0)/(np.sqrt(2*leff+1)*fsky*del_l)

    var_fg_rec   = (cl_fg_rec   + 0)/(np.sqrt(2*leff+1)*fsky*del_l)
    var_fg_rmean = (cl_fg_rmean + 0)/(np.sqrt(2*leff+1)*fsky*del_l)
    var_fg_bias  = (cl_fg_bias  + 0)/(np.sqrt(2*leff+1)*fsky*del_l)
    var_fg_theo  = (cl_fg_theo  + 0)/(np.sqrt(2*leff+1)*fsky*del_l)    
###########    
    m_sky_theo[ind0] = hp.UNSEEN
    #
    m_hi_theo[ ind0] = hp.UNSEEN
    m_hi_rmean[ind0] = hp.UNSEEN    
    m_hi_bias[ ind0] = hp.UNSEEN    
    m_hi_rec[  ind0] = hp.UNSEEN   
    #
    m_fg_theo[ ind0] = hp.UNSEEN
    m_fg_rmean[ind0] = hp.UNSEEN    
    m_fg_bias[ ind0] = hp.UNSEEN    
    m_fg_rec[  ind0] = hp.UNSEEN    

    if show_time:print('Time {0:.4f} seg\n'.format(time.time()-itime)) 
    return {'l'       :leff,
            'HI':{
            "debiased":{'map':m_hi_rec,  'cl':fact*cl_hi_rec,  'var':fact*var_hi_rec  },
            'rmean'   :{'map':m_hi_rmean,'cl':fact*cl_hi_rmean,'var':fact*var_hi_rmean},
            "biased"  :{'map':m_hi_bias, 'cl':fact*cl_hi_bias, 'var':fact*var_hi_bias },
            'theory'  :{'map':m_hi_theo, 'cl':fact*cl_hi_theo, 'var':fact*var_hi_theo },
            "suppresion":supp
                  },
            'FG':{
            "debiased":{'map':m_fg_rec,  'cl':fact*cl_fg_rec,  'var':fact*var_fg_rec  },
            'rmean'   :{'map':m_fg_rmean,'cl':fact*cl_fg_rmean,'var':fact*var_fg_rmean},
            "biased"  :{'map':m_fg_bias, 'cl':fact*cl_fg_bias, 'var':fact*var_fg_bias },
            'theory'  :{'map':m_fg_theo, 'cl':fact*cl_fg_theo, 'var':fact*var_fg_theo },
            "suppresion":suppfg
                  },
            'SKY' :{'map':m_sky_theo},
            'MASK':{'map':mask, 'fsky':fsky}
           }

#####################################################################################################       
#####################################################################################################       
### Loading t recovered data        
###
def reconstruction_data(path=None, dirnames=None, Ns=None, Sims=None, lmax=None, channel=None, show_time=False):
    #path     = '/data1/cross/FGremoval/bingo_nch30_980_1260'
    #dirnames = ['fullsky', 'mask3']
    #Ns       = ['ns2', 'ns3', 'ns4', 'ns5']
    #Sims     = np.array(['sim{}'.format(i) for i in np.arange(0,nsims,1)])
    #lmax     = 599
    #ich      = 19
    for i,imask in enumerate(dirnames):    
        path1 = os.path.join(path,imask)
        itime = time.time()
        #for j,jns in enumerate(np.sort(np.asarray(os.listdir(path1)))[:3]):
        for j,jns in enumerate(Ns):
            path2 = os.path.join(path1,jns)
            jtime = time.time()
            #for k, ksim in enumerate(np.sort(np.asarray(os.listdir(path2)))):
            for k, ksim in enumerate(Sims):
                path3 = os.path.join(path2, ksim) #simulation
                for l,lfile in enumerate(os.listdir(path3)):
                    lfield = lfile.split('_')[0]
                    if not k:
                        if lfield in 'extHI':
                            mdata  = hdata.getmap(dirpath_=path3, filename_=lfile, healpix_readingformat=False, hdu=1)[channel]
                            hi_map = dcopy(mdata) 
                            hi_alm = hp.map2alm(maps=mdata, pol=False, lmax=lmax)                
                        elif lfield in 'extFG':
                            mdata  = hdata.getmap(dirpath_=path3, filename_=lfile, healpix_readingformat=False, hdu=1)[channel]
                            fg_map = dcopy(mdata) 
                            fg_alm = hp.map2alm(maps=mdata, pol=False, lmax=lmax)   
                        else:
                            print(0)
                    else:
                        if lfield in 'extHI':
                            mdata  = hdata.getmap(dirpath_=path3, filename_=lfile, healpix_readingformat=False, hdu=1)[channel]
                            hi_map = np.vstack(( hi_map,dcopy(mdata) ))
                            hi_alm = np.vstack(( hi_alm,hp.map2alm(maps=mdata, pol=False, lmax=lmax) ))
                        elif lfield in 'extFG':
                            mdata  = hdata.getmap(dirpath_=path3, filename_=lfile, healpix_readingformat=False, hdu=1)[channel]
                            fg_map = np.vstack(( fg_map,dcopy(mdata) ))
                            fg_alm = np.vstack(( fg_alm,hp.map2alm(maps=mdata, pol=False, lmax=lmax) ))
                        else:
                            if lfield in 'mixmatrix':
                                pass
                            else:
                                print(0)  
            if not j:            
                fieldmasked_dict = {str(jns):
                                        {'extHI':
                                                {'maps':{'mean':np.mean(hi_map, axis=0), "matrix": hi_map},
                                                 'alms':{'mean':np.mean(hi_alm, axis=0), "matrix": hi_alm},
                                                },
                                         'extFG':
                                                {'maps':{'mean':np.mean(fg_map, axis=0), "matrix": fg_map},
                                                 'alms':{'mean':np.mean(fg_alm, axis=0), "matrix": fg_alm},
                                                }
                                         }
                                       }
            else:
                fieldmasked_dict[str(jns)] = {'extHI':
                                                     {'maps':{'mean':np.mean(hi_map, axis=0), "matrix": hi_map},
                                                      'alms':{'mean':np.mean(hi_alm, axis=0), "matrix": hi_alm},
                                                },
                                               'extFG':
                                                     {'maps':{'mean':np.mean(fg_map, axis=0), "matrix": fg_map},
                                                      'alms':{'mean':np.mean(fg_alm, axis=0), "matrix": fg_alm},
                                                }
                                                 }
        
            if show_time: print(' time {0:.4f} seg'.format(time.time()-jtime)) 
            del hi_map, hi_alm, fg_map, fg_alm, mdata
        if not i:
            field_dict = {imask:fieldmasked_dict}
        else:
            field_dict[imask] = fieldmasked_dict
        del fieldmasked_dict
        if show_time: print('Time {0:.4f} seg\n'.format(time.time()-itime)) 
    return field_dict
        
#Loading simulations from theory (diff seeds)
#This part does not consider mask effect yet
def include_theoretical_data(dict_rfields=None, dict_sfields=None,
                             nside=None,channel=None,lmax=None,
                             map_fg_theo=None,
                             show_time=False):
    if show_time:itime = time.time()
    #dict_rfields:  contains all recovery data from all sims, ns, coverage, and the results (mean, matriz) for FG and HI using alm_hi_sim/alm_fg 
    #              (ex: field_dict. structure: coverage/ns/extHI(FG)/maps|alms|cls/mean|matrix & also HI|FG/maps|alms|cls/mean|matrix)
    #dict_sfields: contains all simulated data for the 'varnames' variable
    #              (ex: dict_all. structure: sims/varnames)
    ##########################  
    for j, jsim in enumerate(dict_sfields.keys()):
        jsim   = str(jsim)
        jseed  = jsim.split('sim')[1]    
        #l_     = dict_all[jsim]['alm_hi_sim'][0].real
        #m_     = dict_all[jsim]['alm_hi_sim'][1]
        alm_hi = dict_sfields[jsim]['alm_hi_sim'][channel+2]
        alm_hi = np.ascontiguousarray(alm_hi)
        map_hi = hp.alm2map(alm_hi, nside=nside, pol=False)
        if not j:
            ALM = dcopy(alm_hi)
            MAP = dcopy(map_hi)
        else:
            ALM = np.vstack((ALM,dcopy(alm_hi)))
            MAP = np.vstack((MAP,dcopy(map_hi)))
            
    dict_rfields['HI'] = {'alms':{'matrix':ALM}, 
                          'maps':{'matrix':MAP}
                         } 
    
    dict_rfields['FG'] = {'alms':{'mean':hp.map2alm(map_fg_theo[channel], pol=False, lmax=lmax)}, 
                          'maps':{'mean':map_fg_theo[channel]}
                          } 
    if show_time: print('Time {0:.4f} seg\n'.format(time.time()-itime))         
                                                  
#####################################################################################################       
#####################################################################################################       

def cl_cross_debiased(dict_rfields=None,   dict_sfields=None,
                      nside=256,     channel=None,lmax=None,
                      ns_flag='ns2', coverage_flag='fullsky',
                      use_dl=False, del_l=1, apply_namaster=False,
                      map_mask=None, show_time=True):
    if show_time:itime = time.time()
    #dict_rfields:  contains all recovery data from all sims, ns, coverage, and the results (mean, matriz) for FG and HI using alm_hi_sim/alm_fg 
    #              (ex: field_dict. structure: coverage/ns/extHI(FG)/maps|alms|cls/mean|matrix & also HI|FG/maps|alms|cls/mean|matrix)
    #dict_sfields: contains all simulated data for the 'varnames' variable
    #              (ex: dict_all. structure: sims/varnames)
    ##########################    
    sim_flag = 'sim0'
    ii       = 1 #sim position on the matrix
    #####
    l_   = dict_sfields[sim_flag]['alm_hi_sim'][0,:]
    if type(lmax)==type(None): lmax=int(l_.real.max())
    mask = dcopy(map_mask)
    ind0 = np.where(mask<=0)[0]
    fsky = (mask.size-ind0.size)/mask.size
    #########################
    #Simulated
    alm_hi_sim = np.ascontiguousarray(dict_sfields[sim_flag]['alm_hi_sim'][int(channel+2)])#[0=l,1=m,...]
    m_hi_sim   = hp.alm2map(alm_hi_sim, lmax=lmax, pol=False , nside=nside)*mask
    ####  
    #Average over recoveries
    m_hi_rmean = dict_rfields[coverage_flag][ns_flag]['extHI']['maps']['mean']*mask
    ####
    #SIM0 biased
    m_hi_bias  = dict_rfields[coverage_flag][ns_flag]['extHI']['maps']['matrix'][int(ii-1)]*mask
    ####
    #Debiased
    supp       = suppression_fact(dict_rfields=dict_rfields, map_mask_=mask, lmax=lmax,
                                  field_sflag='HI' , field_rflag='extHI', 
                                  ns_flag=ns_flag, coverage_flag=coverage_flag)
    #####
    alm_R      = dict_rfields[coverage_flag][ns_flag]['extHI']['alms']['matrix'][0] #realizacao assumida Li=0
    alm_hi_rec = hp.almxfl(alm_R, supp, mmax=lmax)
    m_hi_rec   = hp.alm2map(alm_hi_rec, pol=False,  lmax=lmax, nside=nside,)*mask#, mmax=lmax)
    #cl_hi_rec  = hp.anafast(m_hi_rec, lmax=lmax, pol=False)    
    #########################
    #########################
    #Theory--extracted from the simulation
    alm_k_sim = np.ascontiguousarray(dict_sfields[sim_flag]['alm_kappa_sim'][0+2])#[0=l,1=m,...]
    m_k_sim   = hp.alm2map(alm_k_sim, lmax=lmax, pol=False , nside=nside)*mask
##########
    if apply_namaster:
        import pymaster as nmt
        bpw  = nmt.NmtBin.from_nside_linear(nside, nlb=del_l)    
        leff = bpw.get_effective_ells()
        cl_k_sim   = bpw.bin_cell( hp.anafast(map1=m_k_sim  ,lmax=lmax, pol=False) )
        
        cl_hi_rec   = bpw.bin_cell( hp.anafast(map1=m_hi_rec   ,lmax=lmax, pol=False) )
        cl_hi_rmean = bpw.bin_cell( hp.anafast(map1=m_hi_rmean ,lmax=lmax, pol=False) )
        cl_hi_bias  = bpw.bin_cell( hp.anafast(map1=m_hi_bias  ,lmax=lmax, pol=False) )
        cl_hi_sim   = bpw.bin_cell( hp.anafast(map1=m_hi_sim   ,lmax=lmax, pol=False) )
        
        cl_cx_sim_rec   = bpw.bin_cell( hp.anafast(map1=m_hi_rec   , map2=m_k_sim  ,lmax=lmax, pol=False))
        cl_cx_sim_rmean = bpw.bin_cell( hp.anafast(map1=m_hi_rmean , map2=m_k_sim  ,lmax=lmax, pol=False))
        cl_cx_sim_bias  = bpw.bin_cell( hp.anafast(map1=m_hi_bias  , map2=m_k_sim  ,lmax=lmax, pol=False))
        cl_cx_sim_sim   = bpw.bin_cell( hp.anafast(map1=m_hi_sim   , map2=m_k_sim  ,lmax=lmax, pol=False))
    elif (apply_namaster==False) and (int(del_l)>1):
        include_leftover = True
        cl_k_sim   = cxfs.clsbinned(hp.anafast(map1=m_k_sim  ,lmax=lmax, pol=False), del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        
        leff, cl_hi_rec   = cxfs.clsbinned(hp.anafast(map1=m_hi_rec   ,lmax=lmax, pol=False),   del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)
        cl_hi_rmean = cxfs.clsbinned(hp.anafast(map1=m_hi_rmean ,lmax=lmax, pol=False), del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        cl_hi_bias  = cxfs.clsbinned(hp.anafast(map1=m_hi_bias  ,lmax=lmax, pol=False),  del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        cl_hi_sim   = cxfs.clsbinned(hp.anafast(map1=m_hi_sim  ,lmax=lmax, pol=False),  del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        
        cl_cx_sim_rec   = cxfs.clsbinned(cls=hp.anafast(map1=m_hi_rec   , map2=m_k_sim  ,lmax=lmax, pol=False), del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        cl_cx_sim_rmean = cxfs.clsbinned(cls=hp.anafast(map1=m_hi_rmean , map2=m_k_sim  ,lmax=lmax, pol=False), del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        cl_cx_sim_bias  = cxfs.clsbinned(cls=hp.anafast(map1=m_hi_bias  , map2=m_k_sim  ,lmax=lmax, pol=False), del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        cl_cx_sim_sim   = cxfs.clsbinned(cls=hp.anafast(map1=m_hi_sim   , map2=m_k_sim  ,lmax=lmax, pol=False), del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
    else:
        leff = np.arange(cl_fg_rec.size)
    fact = leff*(leff+1)/2/np.pi if use_dl else 1
########### 

    var_cx_sim_rec   = np.sqrt((cl_cx_sim_rec**2   + cl_k_sim*cl_hi_rec  )/((2*leff + 1)*fsky*del_l))
    var_cx_sim_rmean = np.sqrt((cl_cx_sim_rmean**2 + cl_k_sim*cl_hi_rmean)/((2*leff + 1)*fsky*del_l))
    var_cx_sim_bias  = np.sqrt((cl_cx_sim_bias**2  + cl_k_sim*cl_hi_bias )/((2*leff + 1)*fsky*del_l))
    var_cx_sim_sim   = np.sqrt((cl_cx_sim_sim**2   + cl_k_sim*cl_hi_sim  )/((2*leff + 1)*fsky*del_l))
    
    #########################
    return {'cls' :{'sim_rec'  : fact*cl_cx_sim_rec,
                    'sim_rmean': fact*cl_cx_sim_rmean,
                    'sim_bias' : fact*cl_cx_sim_bias, 
                    'sim_sim'  : fact*cl_cx_sim_sim},
            'vars':{'sim_rec'  : fact*var_cx_sim_rec,
                    'sim_rmean': fact*var_cx_sim_rmean,
                    'sim_bias' : fact*var_cx_sim_bias, 
                    'sim_sim'  : fact*var_cx_sim_sim},
            'l':leff}
