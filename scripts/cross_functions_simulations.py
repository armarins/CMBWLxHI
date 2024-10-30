import numpy as np
from copy import deepcopy as dcopy
import time
import healpy as hp
#####################################################################################################
#####################################################################################################
#### CREATING CONVERGENCE-field SIMULATION
def crosscorr_coeff(cl_A=None, cl_B=None, cl_AB=None):
    return cl_AB/np.sqrt(cl_A*cl_B)    
    
def cross_simulations_correlated_from_dict(clf1_dict_=None, clf2_dict_=None, clcx_dict_=None, seed_k=None, seed_hi=None, show_time=False):
    import cmath
    tax=1000 #to distinguish the alm_hi
    if show_time: itime = time.time() 
    l             = clf1_dict_['l']
    alm_hi_dict   = {'l': dcopy(l)}
    alm_k_c_dict  = {'l': dcopy(l)}
    alpha_l_ratio = 0
    vec_seeds_hi = np.array([])

    sim_names = dcopy(clf1_dict_)
    del sim_names['l']
    sim_names = np.array(dcopy(list(sim_names.keys())))

    for j,jbin in enumerate(sim_names):
        jseed_hi         = seed_hi + j*tax
        vec_seeds_hi     = np.hstack((vec_seeds_hi, jseed_hi))
        np.random.seed(jseed_hi)
        alm_hi_jbin        = hp.synalm(clf1_dict_[jbin])
        alm_k_cj           = hp.almxfl( alm=alm_hi_jbin, fl=(clcx_dict_[jbin]/clf1_dict_[jbin]) )
        alm_hi_dict[jbin]  = dcopy(alm_hi_jbin)
        alm_k_c_dict[jbin] = dcopy(alm_k_cj)
        alpha_l_ratio      += clcx_dict_[jbin]*clcx_dict_[jbin]/clf1_dict_[jbin]
    
    np.random.seed(seed_k)
    alm_k_uncorr= hp.synalm(clf2_dict_['bin1']-alpha_l_ratio)

    sim_names = dcopy(clf1_dict_)
    del sim_names['l']
    sim_names = np.array(dcopy(list(sim_names.keys())))
    
    alm_k_sim = dcopy(alm_k_uncorr)
    for j,jbin in enumerate(sim_names):
        alm_k_sim += alm_k_c_dict[jbin]
        if not j: 
            alm_hi_sim   = dcopy(alm_hi_dict[jbin])
            alm_k_corr   = dcopy(alm_k_c_dict[jbin])
            cl_hi_sim    = dcopy(hp.alm2cl(alm_hi_dict[jbin]))#, lmax=int(l.max())))
            cl_hi_theory = dcopy(clf1_dict_[jbin])
        else:     
            alm_hi_sim   = np.vstack(( alm_hi_sim  , dcopy(alm_hi_dict[jbin])  ))
            alm_k_corr   = np.vstack((   alm_k_corr, dcopy(alm_k_c_dict[jbin]) ))
            cl_hi_sim    = np.vstack(( cl_hi_sim   , dcopy(hp.alm2cl(alm_hi_dict[jbin])) ))#, lmax=int(l.max()))) ))
            cl_hi_theory = np.vstack(( cl_hi_theory, dcopy(clf1_dict_[jbin])    ))

    sim_names = dcopy(clcx_dict_)
    del sim_names['l']
    sim_names = np.array(dcopy(list(sim_names.keys())))

    
    for j,jbin in enumerate(sim_names):
        if not j:
            cl_cx_sim    = hp.alm2cl(alm_k_sim, alm_hi_sim[j])#, lmax=int(l.max()) )
            cl_cx_theory = dcopy(clcx_dict_[jbin])
        else:
            cl_cx_sim    = np.vstack(( cl_cx_sim   , hp.alm2cl(alm_k_sim, alm_hi_sim[j] ) )) #, lmax=int(l.max()) ) ))
            cl_cx_theory = np.vstack(( cl_cx_theory, dcopy(clcx_dict_[jbin]) ))
    
    cl_k_sim  = hp.alm2cl(alm_k_sim)#,  lmax=int(l.max()) )
    if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))  
    return {'alm_kappa_sim': alm_k_sim , 'alm_hi_sim': alm_hi_sim, 
            'alm_kappa_corr':alm_k_corr, 'alm_kappa_uncorr':alm_k_uncorr,
            'cl_kappa_sim':  cl_k_sim  , 'cl_kappa_theory': clf2_dict_['bin1'], 
            'cl_hi_sim'   :  cl_hi_sim , 'cl_hi_theory'   : cl_hi_theory, 
            'cl_cross_sim':  cl_cx_sim , 'cl_cross_theory': cl_cx_theory, 
            'l':l, 'seed_hi':np.int32(vec_seeds_hi), 'seed_kunc':np.int32(seed_k)
            }
def create_HI_maps_from_dicts(dicts_=None, nside=None):
    for num in range(dicts_['alm_hi_sim'].shape[0]):
        alm_hi_num = np.ascontiguousarray(dicts_['alm_hi_sim'][num,:])
        if not num:  Mhi_ = hp.alm2map( alm_hi_num, nside=nside, pol=False)
        else:        Mhi_ = np.vstack(( Mhi_, hp.alm2map( alm_hi_num, nside=nside, pol=False) ))
    return Mhi_

def create_HI_maps_from_cls(cl_matrix=None, nside=None, seed_hi=None, tax=None):
    nch = cl_matrix.shape[0]
    for jbin in range(nch):
        jseed_hi    = seed_hi + jbin*tax
        np.random.seed(jseed_hi)    
        alm_hi_num = hp.synalm(cl_matrix[jbin,:],verbose=False)
        Mhi_ = hp.alm2map( alm_hi_num, nside=nside, pol=False) if not jbin else np.vstack(( Mhi_, hp.alm2map( alm_hi_num, nside=nside, pol=False) ))
    return Mhi_
    
def cross_simulations_correlated_from_matrix(clf1_=None, clf2_=None, clcx_=None, seed_k=None, seed_hi=None, show_time=False, beta = None, fact=1, tax=300):
    ##############################################
    #clf1: matrix of Cl(HI)        - [nbins, nells]
    #clf2: matrix of Cl(kappa)     - [1, nells]
    #clcx: matrix of Cl(kappa, HI) - [nbins, nells]
    #fact: to amplify the cross artificially
    #tax : to distinguish the alm_hi
    ##############################################
    if show_time: itime = time.time() 
    nch, l        = clf1_.shape[0], np.arange(clf1_.shape[1])
    if type(beta) != type(None):
        clcx_ = beta*np.sqrt(clf1_*clf2_)        
    vec_seeds_hi = np.array([])
    alpha_l_ratio = 0

    for jbin in range(nch):
        jseed_hi     = seed_hi + jbin*tax
        vec_seeds_hi = np.hstack((vec_seeds_hi, jseed_hi))
        np.random.seed(jseed_hi)
        alm_hi_jbin = hp.synalm(clf1_[jbin,:])
        alm_k_cj    = hp.almxfl( alm=alm_hi_jbin, fl=(clcx_[jbin,:]/clf1_[jbin,:]) )
        if not jbin:
            alm_hi_vec  = dcopy(alm_hi_jbin)
            alm_k_c_vec = dcopy(alm_k_cj)
        else:
            alm_hi_vec  = np.vstack(( alm_hi_vec , dcopy(alm_hi_jbin) ))
            alm_k_c_vec = np.vstack(( alm_k_c_vec, dcopy(alm_k_cj) ))
        alpha_l_ratio  += clcx_[jbin,:]*clcx_[jbin,:]/clf1_[jbin,:]
    
    np.random.seed(seed_k)
    alm_k_uncorr = hp.synalm(clf2_[0,:]-fact*alpha_l_ratio)
    
    alm_k_sim   = dcopy(alm_k_uncorr)
    if nch==1:
        alm_k_sim += fact*alm_k_c_vec
        alm_k_sim  = alm_k_sim[np.newaxis,:]
        alm_hi_sim = dcopy(alm_hi_vec)
        alm_hi_sim = alm_hi_sim[np.newaxis,:]
        alm_hi_vec = alm_hi_vec[np.newaxis,:]
        cl_hi_sim  = hp.alm2cl(alm_hi_sim[0])[np.newaxis,:]
        cl_cx_sim  = hp.alm2cl(alm_k_sim[0], alm_hi_sim[0])[np.newaxis,:]
        #cl_hi_sim  = dcopy(cl_cx_sim[np.newaxis,:])
        #cl_cx_sim  = dcopy(cl_cx_sim[np.newaxis,:])

    elif nch>1:
        for jbin in range(nch):
            alm_k_sim += fact*alm_k_c_vec[jbin]
            alm_hi_sim = dcopy(alm_hi_vec[jbin]) if not jbin else np.vstack(( alm_hi_sim, dcopy(alm_hi_vec[jbin]) ))
        for jbin in range(nch):
            cl_hi_sim  = hp.alm2cl(alm_hi_sim[jbin])            if not jbin else np.vstack(( cl_hi_sim, hp.alm2cl(           alm_hi_sim[jbin]) ))  
            cl_cx_sim  = hp.alm2cl(alm_k_sim, alm_hi_sim[jbin]) if not jbin else np.vstack(( cl_cx_sim, hp.alm2cl(alm_k_sim, alm_hi_sim[jbin]) ))  
    else:
        raise NameError
    
    cl_k_sim   = hp.alm2cl(alm_k_sim)
    if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))  
    return {'alm_hi_sim':       alm_hi_vec,            
            'alm_kappa_sim':    alm_k_sim,   
            'alm_kappa_uncorr': alm_k_uncorr,   
            'alm_kappa_corr':   alm_k_c_vec,   
            'alm_kappa_corr_sum': np.sum(alm_k_c_vec,axis=0)[np.newaxis,:],
            'cl_kappa_sim': hp.alm2cl(alm_k_sim),
            'cl_hi_sim'   : cl_hi_sim,
            'cl_cross_sim': cl_cx_sim ,
            'cl_hi_theory'   : clf1_ ,
            'cl_kappa_theory': clf2_[0,:] ,
            'cl_cross_theory': clcx_,
            'l':l, 'seed_hi':np.int32(vec_seeds_hi), 'seed_kunc':np.int32(seed_k)
            } 
    
def dictionary_cross_simulations_quantities_from_dict(clf1_dict_=None, clf2_dict_=None, clcx_dict_=None, seed_k=None, seed_hi=None, nsims=None, show_time=False):
    if show_time: itime = time.time() 
    seeds_k  = np.arange(seed_k , seed_k + nsims +1, 1)
    seeds_hi = np.arange(seed_hi, seed_hi+ nsims +1, 1)
    for j, (jsim_k, jsim_h) in enumerate(zip(seeds_k, seeds_hi)):
        dicts_ = cxfs.cross_simulations_correlated_from_dict(clf1_dict_, clf2_dict_, clcx_dict_, jsim_k, jsim_h, False)
        #_dict = cxfs.simulated_kappa_var(clf1_=clf1_4kappa, clf2_=clf2_4kappa, clcx_=clcx_4kappa, seed_k=8001, seed_hi=int(9000+0), show_time=0,fact=1, tax=300)                                
        if not j:
            dict_sims_ = {'sim{}'.format(j): dicts_}
        else:
            dict_sims_['sim{}'.format(j)] = dicts_
    if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))            
    return dict_sims_
        
def dictionary_cross_simulations_quantities_from_matrix(clf1_=None, clf2_=None, clcx_=None, seed_k=None, seed_hi=None, 
                                                        channel_min_corr=None,channel_max_corr=None, nsims=None, fact=1, 
                                                        beta=None, tax=300, show_time=False):
    if show_time: itime = time.time() 
    clf1_4kappa = dcopy(clf1_)
    clf2_4kappa = dcopy(clf2_)
    clcx_4kappa = dcopy(clcx_)
    #
    if type(channel_min_corr)!=type(None) and type(channel_max_corr)!=type(None):
        clcx_4kappa[  :channel_min_corr-1,:] = 0*clcx_[  :channel_min_corr-1,:]    
        clcx_4kappa[channel_max_corr:,:] = 0*clcx_[channel_max_corr:,:] 
    ####
    #seeds_k  = np.arange(seed_k , seed_k + nsims +1, 1)
    seeds_hi = np.arange(seed_hi, seed_hi+ nsims +1, 1)
    #for j, (jsim_k, jsim_h) in enumerate(zip(seeds_k, seeds_hi)):
    for j, jsim_h in enumerate(seeds_hi):
        dicts_ = cross_simulations_correlated_from_matrix(clf1_=clf1_4kappa, clf2_=clf2_4kappa, clcx_=clcx_4kappa, 
                                                              seed_k=seed_k, seed_hi=jsim_h, fact=fact, beta=beta, tax=tax, show_time=show_time)                                
        if not j:
            dict_sims_ = {'sim{}'.format(j): dicts_}
        else:
            dict_sims_['sim{}'.format(j)] = dicts_
    if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))            
    return dict_sims_



def loop_cls_kappa(jseed=None, del_l=1, nside=256, show_time=False):
    import pymaster as nmt
    if show_time: itime = time.time() 
    dict_vars = cross_simulations_correlated(clf1_dict_=clf1_dict, clf2_dict_=clf2_dict, clcx_dict_=clcx_dict,
                                             seed_k=jseed, seed_hi=jseed+1, show_time=False)
    b = nmt.NmtBin.from_nside_linear(nside, nlb=del_l)

    for i,ialm in enumerate(dict_vars['alm_kappa_corr']):
        icl = b.bin_cell( hp.alm2cl(ialm) )
        if not i: Cl_k_corr= dcopy(icl)        
        else:     Cl_k_corr= np.vstack(( Cl_k_corr, icl ))
        del icl
    
    Cl_k_corr_comps = dcopy(Cl_k_corr) 
    Cl_k_corr_sum   = np.sum(Cl_k_corr, axis=0) 
    Cl_k_uncorr     = b.bin_cell( hp.alm2cl(dict_vars['alm_kappa_uncorr']) )
    Cl_k_theory     = b.bin_cell( dcopy(dict_vars['cl_kappa_theory']) )
    Cl_k_sim        = b.bin_cell( dcopy(dict_vars['cl_kappa_sim']) )
    l               = b.get_effective_ells()

    if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))      
    return {'corr_comps':Cl_k_corr_comps,
            'corr_sum'  :Cl_k_corr_sum,
            'uncorr'    :Cl_k_uncorr,
            'theory'    :Cl_k_theory,
            'sim'       :Cl_k_sim,
            'l'         :l}

def loop_cls_cxKalms(jseed=None, del_l=1, lmin=None, lmax=None, include_leftover=False, nside=256, show_time=False):
    import pymaster as nmt
    if show_time: itime = time.time() 
    dict_vars = cross_simulations_correlated(clf1_dict_=clf1_dict, clf2_dict_=clf2_dict, clcx_dict_=clcx_dict,
                                             seed_k=jseed, seed_hi=jseed+1, show_time=False)
    ##
    b    = nmt.NmtBin.from_nside_linear(nside, nlb=del_l)
    leff = b.get_effective_ells()
    #cl_binned = b.bin_cell(cl_true)
    for i,ialm in enumerate(dict_vars['alm_kappa_corr']):
        cl_cxu = hp.alm2cl(alms1=dict_vars['alm_kappa_uncorr'], alms2=ialm)
        #cl_cxu = cxfs.clsbinned(cls=cl_cxu, del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        cl_cxu = b.bin_cell(cl_cxu)
        if not i: Cl_cxu_comps= cl_cxu
        else:     Cl_cxu_comps= np.vstack(( Cl_cxu_comps, cl_cxu ))

    Cl_cxu_sum = hp.alm2cl(alms1=dict_vars['alm_kappa_uncorr'],alms2=np.sum(dict_vars['alm_kappa_corr'], axis=0))
    Cl_cxu_sum = b.bin_cell(Cl_cxu_sum)
    if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))      
    return {'comps':Cl_cxu_comps,
            'sum'  :Cl_cxu_sum,
            'l'    :b.get_effective_ells()}


def create_HI_maps_from_dicts(dicts_=None, nside=None):
    for num in range(dicts_['alm_hi_sim'].shape[0]):
        alm_hi_num = np.ascontiguousarray(dicts_['alm_hi_sim'][num,:])
        if not num:  Mhi_ = hp.alm2map( alm_hi_num, nside=nside, pol=False)
        else:        Mhi_ = np.vstack(( Mhi_, hp.alm2map( alm_hi_num, nside=nside, pol=False) ))
    return Mhi_

def create_HI_maps_from_cls(cl_matrix=None, nside=None, seed_hi=None, tax=None):
    nch = cl_matrix.shape[0]
    for jbin in range(nch):
        jseed_hi    = seed_hi + jbin*tax
        np.random.seed(jseed_hi)    
        alm_hi_num = hp.synalm(cl_matrix[jbin,:])#,verbose=False)
        Mhi_ = hp.alm2map( alm_hi_num, nside=nside, pol=False) if not jbin else np.vstack(( Mhi_, hp.alm2map( alm_hi_num, nside=nside, pol=False) ))
    return Mhi_


#####################################################################################################
#####################################################################################################
#### STATs
def stat_simulations(dict_sims_=None, lmin=None, lmax=None, nside=256, del_l=1, ich=None, use_dl=True, include_leftover=False, show_time=False):
    if show_time: itime = time.time()
    N = len(dict_sims_.keys())
    if int(del_l)==1:            
        cl_k_mean  = np.zeros_like(dict_sims_['sim0']['cl_kappa_sim'])
        cl_hi_mean = np.zeros_like(dict_sims_['sim0']['cl_hi_sim'   ][ich])
        cl_cx_mean = np.zeros_like(dict_sims_['sim0']['cl_cross_sim'][ich])
        #####
        #l_   = dict_sims_['sim0']['l']
        l   = np.arange(dict_sims_['sim0']['l'].size)
        fact = fact=l*(l+1)/2/np.pi if use_dl else 1        
        #####
        for j, jsim in enumerate(dict_sims_.keys()):
            cl_k_mean  += fact*dict_sims_[jsim]['cl_kappa_sim']/N
            cl_hi_mean += fact*dict_sims_[jsim]['cl_hi_sim'   ][ich]/N
            cl_cx_mean += fact*dict_sims_[jsim]['cl_cross_sim'][ich]/N        
        
        cl_k_std  = np.zeros_like(cl_k_mean )
        cl_hi_std = np.zeros_like(cl_hi_mean)
        cl_cx_std = np.zeros_like(cl_cx_mean)
        for j, jsim in enumerate(dict_sims_.keys()):
            cl_k_std  += ( fact*dict_sims_[jsim]['cl_kappa_sim']      - cl_k_mean )**2/N
            cl_hi_std += ( fact*dict_sims_[jsim]['cl_hi_sim'   ][ich] - cl_hi_mean)**2/N
            cl_cx_std += ( fact*dict_sims_[jsim]['cl_cross_sim'][ich] - cl_cx_mean)**2/N   
            
    elif int(del_l)>=1:
        import pymaster as nmt
        b = nmt.NmtBin.from_nside_linear(nside, nlb=del_l)
        cl_k_mean  = b.bin_cell( dict_sims_['sim0']['cl_kappa_sim']      )
        cl_hi_mean = b.bin_cell( dict_sims_['sim0']['cl_hi_sim'   ][ich] )
        cl_cx_mean = b.bin_cell( dict_sims_['sim0']['cl_cross_sim'][ich] )

        cl_k_mean  = np.zeros_like(cl_k_mean)
        cl_hi_mean = np.zeros_like(cl_hi_mean)
        cl_cx_mean = np.zeros_like(cl_cx_mean)      
        #####
        l    = b.get_effective_ells()
        fact = fact=l*(l+1)/2/np.pi if use_dl else 1
        #####
        for j, jsim in enumerate(dict_sims_.keys()):
            cl_k_mean  += fact*b.bin_cell( dict_sims_[jsim]['cl_kappa_sim']      )/N
            cl_hi_mean += fact*b.bin_cell( dict_sims_[jsim]['cl_hi_sim'   ][ich] )/N
            cl_cx_mean += fact*b.bin_cell( dict_sims_[jsim]['cl_cross_sim'][ich] )/N
        
        cl_k_std  = np.zeros_like(cl_k_mean )
        cl_hi_std = np.zeros_like(cl_hi_mean)
        cl_cx_std = np.zeros_like(cl_cx_mean)
        for j, jsim in enumerate(dict_sims_.keys()):
            cl_k_std  += ( fact*b.bin_cell( dict_sims_[jsim]['cl_kappa_sim']      ) - cl_k_mean )**2/N
            cl_hi_std += ( fact*b.bin_cell( dict_sims_[jsim]['cl_hi_sim'   ][ich] ) - cl_hi_mean)**2/N
            cl_cx_std += ( fact*b.bin_cell( dict_sims_[jsim]['cl_cross_sim'][ich] ) - cl_cx_mean)**2/N
    else:
        raise Exception
    
    cl_k_std  = np.sqrt(cl_k_std )    
    cl_hi_std = np.sqrt(cl_hi_std)    
    cl_cx_std = np.sqrt(cl_cx_std)    
    if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))                
    return {'cl_kappa_sim': {'mean': cl_k_mean ,
                             'std' : cl_k_std
                            },
            'cl_hi_sim'   : {'mean': cl_hi_mean ,
                             'std' : cl_hi_std
                            },
            'cl_cross_sim': {'mean': cl_cx_mean ,
                             'std' : cl_cx_std
                            },
            'l':l
           }

#####################################################################################################
#### kappa components   

def loop_cls_kappa(jseed=None, del_l=1, nside=256, show_time=False):
    import pymaster as nmt
    if show_time: itime = time.time() 
    dict_vars = cross_simulations_correlated(clf1_dict_=clf1_dict, clf2_dict_=clf2_dict, clcx_dict_=clcx_dict,
                                             seed_k=jseed, seed_hi=jseed+1, show_time=False)
    b = nmt.NmtBin.from_nside_linear(nside, nlb=del_l)

    for i,ialm in enumerate(dict_vars['alm_kappa_corr']):
        icl = b.bin_cell( hp.alm2cl(ialm) )
        if not i: Cl_k_corr= dcopy(icl)        
        else:     Cl_k_corr= np.vstack(( Cl_k_corr, icl ))
        del icl
    
    Cl_k_corr_comps = dcopy(Cl_k_corr) 
    Cl_k_corr_sum   = np.sum(Cl_k_corr, axis=0) 
    Cl_k_uncorr     = b.bin_cell( hp.alm2cl(dict_vars['alm_kappa_uncorr']) )
    Cl_k_theory     = b.bin_cell( dcopy(dict_vars['cl_kappa_theory']) )
    Cl_k_sim        = b.bin_cell( dcopy(dict_vars['cl_kappa_sim']) )
    l               = b.get_effective_ells()

    if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))      
    return {'corr_comps':Cl_k_corr_comps,
            'corr_sum'  :Cl_k_corr_sum,
            'uncorr'    :Cl_k_uncorr,
            'theory'    :Cl_k_theory,
            'sim'       :Cl_k_sim,
            'l'         :l}
#####################################################################################################
#####################################################################################################

def get_alm_simulated_file(filename_=None):
    with open(filename_, 'r') as reader:
        almfile = dcopy(reader.readlines())
    #almfile_header = np.asarray(almfile)[0]
    almfile = np.asarray(almfile)[1:]
    for i,f in enumerate(almfile[:]):
        alms_file = f.split('\n')[0].split(' ')
        if not i: 
            l_f   = int(alms_file[0])
            m_f   = int(alms_file[1])
            alm_f = np.array([complex(ialm) for ialm in alms_file[2:]])
        else:
            l_f   = np.vstack((   l_f, int(alms_file[0]) ))
            m_f   = np.vstack((   m_f, int(alms_file[1]) ))
            alm_f = np.vstack(( alm_f, np.array([complex(ialm) for ialm in alms_file[2:]]) )) 
    l_f   = np.asarray(l_f).T
    m_f   = np.asarray(m_f).T
    alm_f = np.asarray(alm_f).T
    
    return {'l':l_f, 'm':m_f, 'alm':alm_f}
    
def get_dict_simulations(nreals=100, first_seed=9000, Chi_dict=None, Ckk_vec=None, Cross_dict=None, l_=None,
                         jseed_kunc=None, jseed_hi=None, nbins_=30, show_time=False, fact=1, sort_bin_sim=True):
    if list(Chi_dict.keys())[1]=='bin1':## to know the first bin name assumed: bin1 or bin0
        bin_start_1 = True
    else: 
        bin_start_1 = False
    ###
    seeds_HIweighted  = np.arange(first_seed              , first_seed + nreals       , 1, dtype=np.int16) #only build the seed vectors
    seeds_Kuncorr     = np.arange(first_seed + nreals + 10, first_seed + 2*nreals + 10, 1, dtype=np.int16) #idem
    ###
    if bin_start_1: 
        bins_vec=np.arange(1, nbins_+1, 1)
    else:
        bins_vec=np.arange(0, nbins_+0, 1)
    if sort_bin_sim:
        if show_time: itime = time.time()
        for b, bin_ in enumerate( bins_vec.astype(str) ):
            for j, (jseed_hi,jseed_kunc) in enumerate( zip(seeds_HIweighted, seeds_Kuncorr) ):
                dict_isim = cmbwl_x_hi_sims(Chi_dict=Chi_dict, Ckk_vec=Ckk_vec, Cross_dict=Cross_dict, bin_='bin{}'.format(bin_), l_=l_,  jseed_kunc=jseed_kunc, jseed_hi=jseed_hi, fact=fact)
                if not j: 
                    dict_sims = {'sim{0}'.format(j):dict_isim}
                else:
                    dict_sims['sim{0}'.format(j)] = dict_isim
            if not b: dict_sim_bins = {'bin{0}'.format(bin_):dict_sims}
            else:     dict_sim_bins['bin{0}'.format(bin_)] = dict_sims
            del dict_sims, dict_isim
        if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))
    ###
    else:    
        if show_time: itime = time.time()    
        for j, (jseed_hi,jseed_kunc) in enumerate( zip(seeds_HIweighted, seeds_Kuncorr) ):
            for b, bin_ in enumerate( bins_vec.astype(str) ):
                dict_isim_bbin = cmbwl_x_hi_sims(Chi_dict=Chi_dict, Ckk_vec=Ckk_vec, Cross_dict=Cross_dict, bin_='bin{}'.format(bin_), l_=l_,  jseed_kunc=jseed_kunc, jseed_hi=jseed_hi, fact=fact)
                if not b:
                    dict_isims_bins = {'bin{0}'.format(bin_):dict_isim_bbin}
                else:
                    dict_isims_bins['bin{0}'.format(bin_)] = dict_isim_bbin
            if not j: 
                dict_sim_bins = {'sim{0}'.format(j):dict_isims_bins}
            else:
                dict_sim_bins['sim{0}'.format(j)] = dict_isims_bins    
            del dict_isims_bins, dict_isim_bbin
        if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))    
    return dict_sim_bins
    
def get_file_simulated_data(path_=None, varnames = ['alm_hi_sim', 'cl_hi_sim', 'cl_kappa_sim', 'cl_cross_sim', 'alm_kappa_sim'], show_time=False, nsims_to_use=100):
    import os
    if show_time: itime = time.time() 
    sims = np.sort(np.array([isim for isim in list(os.listdir(path_)) if 'sim' in isim])) #to avoid other files/directories but 'sim'
    if (sims.size>nsims_to_use)*(type(nsims_to_use)==int or type(nsims_to_use)==float):
        nsims_to_use = int(nsims_to_use)
        sims = sims[:nsims_to_use]
    #################################
    for i, isim in enumerate(sims):
        if 'sim' in isim:	
            pathdir  = '{p}/{s}'.format(p=path_, s=isim)
            for j,jname in enumerate(varnames):
                if 'cl' in jname:
                    filename = os.path.join(pathdir, '.'.join((jname,'txt')) )
                    data     = np.loadtxt(filename).T
#                    l        = data[0]
#                    data     = data[1:]                    
                elif 'cross' in jname:
                    filename = os.path.join(pathdir, '.'.join((jname,'txt')) )
                    data     = np.loadtxt(filename).T                    
#                    l        = data[0]
#                    data     = data[1:]                                        
                elif 'alm' in jname:
                    filename = os.path.join(pathdir, '.'.join((jname,'txt')) )
                    #data = get_alm_simulated_file(filename_=filename) #{'l','m','alm'}                                            
                    alm_real    = np.loadtxt( os.path.join( pathdir, '.'.join(( '_'.join(( jname, 'real' )), 'txt')) ) ).T
                    l,m         = alm_real[0], alm_real[1]
                    alm_real    = alm_real[2:]
                    alm_imag    = np.loadtxt( os.path.join( pathdir, '.'.join(( '_'.join(( jname, 'imag' )), 'txt')) ) ).T[2:]
                    data        = alm_real + alm_imag*1.0j
                    data        = np.vstack(( l, m,  data))      
                else:
                    raise Exception('{} file no exist'.format(jname))
                if not i + j: 
                    dict_all_ = {isim: {jname: data}}
                else:
                    if not j: 
                        dict_all_[isim] =  {jname: data}
                    else:
                        dict_all_[isim][jname] = data  
    if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))  
    return dict_all_	

def get_stat_from_sims(dict_sims=None, varname='cl_cross_sim', sim_name_standard='sim0'):
    #dict_sims_ = dcopy(dict_sims)
    sim_names = np.asarray( list(dict_sims.keys()) )
    nch_      = int( dict_sims['sim0']['cl_hi_sim'].shape[0]-1 )
    if not sim_name_standard in sim_names:
        sim_name_standard = sim_names[0]
    for j, jsim in enumerate(sim_names):
        if jsim==sim_name_standard: cl_sim = dict_sims[jsim][varname][1:,:]
        _cl_ = dict_sims[jsim][varname][1:,:].flatten() if not j else np.vstack(( _cl_,dict_sims[jsim][varname][1:,:].flatten() ))

    leff_ = b_namaster.get_effective_ells() if type(b_namaster)!=type(None) else np.arange(cl_sim.shape[1])
    feff_ = leff_*(leff_+1)/2/np.pi if use_dl else 1

    return {'matrix':feff_*b_namaster.bin_cell(cl_sim), 
            'mean':  feff_*b_namaster.bin_cell( np.average(_cl_, axis=0).reshape(nch_,-1)), 
            'std':np.std(_cl_, axis=0).reshape(nch_,-1)}

def get_stat_from_sims_ibin(dict_sims=None, varname='cl_cross_sim', sim_name_standard='sim0', b_namaster=None, use_dl=False, ibin=1):
    #dict_sims_ = dcopy(dict_sims)
    sim_names = np.asarray( list(dict_sims.keys()) )
    nch_      = int( dict_sims['sim0']['cl_hi_sim'].shape[0]-1 )
    if not sim_name_standard in sim_names:
        sim_name_standard = sim_names[0]
    for j, jsim in enumerate(sim_names):
        if jsim==sim_name_standard: cl_sim = dict_sims[jsim][varname][1:,:]
        _cl_ = dict_sims[jsim][varname][1:,:].flatten() if not j else np.vstack(( _cl_,dict_sims[jsim][varname][1:,:].flatten() ))

    leff_ = b_namaster.get_effective_ells() if type(b_namaster)!=type(None) else np.arange(cl_sim.shape[1])
    feff_ = leff_*(leff_+1)/2/np.pi if use_dl else 1

    return {'matrix':feff_*b_namaster.bin_cell(cl_sim[ibin,:]), 
            'mean':  feff_*b_namaster.bin_cell( np.average(_cl_, axis=0).reshape(nch_,-1)[ibin,:]), 
            'std':   feff_*b_namaster.bin_cell( np.std(_cl_, axis=0).reshape(nch_,-1)[ibin,:] )
           }


def get_data_plot_residuals(ns=3, ich=47, nu_min=45, nu_max=50, seed_k=8001, seed_hi=9000, tax=300, clf1_=None, clf2_=None, clcx_=None, nside=256, fg_maps=None):
    # It is required the sklearn python package
    from sklearn.decomposition import FastICA
    dict_ = dictionary_cross_simulations_quantities_from_matrix(clf1_=clf1_, clf2_=clf2_, clcx_=clcx_, seed_k=seed_k, seed_hi=seed_hi, 
                                                                         channel_min_corr=nu_min, channel_max_corr=nu_max,
                                                                         nsims=1, fact=1, tax=tax, show_time=False)
        
    dicts = dict_['sim0']
    del dict_
    ########
    alm_kappa_sim = dcopy(dicts['alm_kappa_sim'][:])
    Mhi = create_HI_maps_from_cls(cl_matrix=clf1_, nside=nside, seed_hi=seed_hi, tax=tax)
    Mfg = dcopy(fg_maps)
    X   = dcopy(Mhi + Mfg)#*MASK
    ########
    ica = FastICA(n_components = int(ns), whiten = 'arbitrary-variance', algorithm='parallel', fun='logcosh', 
                 fun_args     = None, max_iter=200, tol=0.0001, w_init=None, whiten_solver='svd', 
                 random_state = None)
    S   = ica.fit_transform(X.T) 
    A   = ica.mixing_
    W   = np.dot( np.linalg.inv(np.dot(A.T,A)), A.T)
    Wfg = np.dot(A,W)
    Xfg = np.dot(Wfg, X)
    Xhi = X - Xfg
    Xfg_lkg = Mfg-np.dot(Wfg,Mfg)
    Xhi_lkg =   0+np.dot(Wfg,Mhi)
    ######################
    for jch in range(Mhi.shape[0]):
        Cl_hi_sim_n = hp.anafast(map1=Mhi[jch,:]    , pol=False) if not jch else np.vstack(( Cl_hi_sim_n, hp.anafast(map1=Mhi[jch,:])     ))
        Cl_hi_rec_n = hp.anafast(map1=Xhi[jch,:]    , pol=False) if not jch else np.vstack(( Cl_hi_rec_n, hp.anafast(map1=Xhi[jch,:])     ))
        Cl_hi_lkg_n = hp.anafast(map1=Xhi_lkg[jch,:], pol=False) if not jch else np.vstack(( Cl_hi_lkg_n, hp.anafast(map1=Xhi_lkg[jch,:]) ))
        Cl_fg_sim_n = hp.anafast(map1=Mfg[jch,:]    , pol=False) if not jch else np.vstack(( Cl_fg_sim_n, hp.anafast(map1=Mfg[jch,:])     ))
        Cl_fg_rec_n = hp.anafast(map1=Xfg[jch,:]    , pol=False) if not jch else np.vstack(( Cl_fg_rec_n, hp.anafast(map1=Xfg[jch,:])     ))
        Cl_fg_lkg_n = hp.anafast(map1=Xfg_lkg[jch,:], pol=False) if not jch else np.vstack(( Cl_fg_lkg_n, hp.anafast(map1=Xfg_lkg[jch,:]) ))
    
    return {'HI': {'sim': Cl_hi_sim_n, 'rec': Cl_hi_rec_n, 'lkg':Cl_hi_lkg_n},
            'FG': {'sim': Cl_fg_sim_n, 'rec': Cl_fg_rec_n, 'lkg':Cl_fg_lkg_n} }

#####################################################################################################
#####################################################################################################
#### Changing formating
####
def from_dict2sim_ells(dict_all_ = None, bin_ = None, var_name_=None): #get a vector storage in different keys bin/sim --> return the same in a matrix2D  [nsims, nells] for the respective bin
    #provide bin_ as used in the dictionary: 'bin1', 'bin2', ...
    for i, isim in enumerate(dict_all_[bin_].keys()):
        if not i: var_matrix_ibin = dict_all_[bin_][isim][var_name_]
        else: var_matrix_ibin = np.vstack(( var_matrix_ibin, dict_all_[bin_][isim][var_name_]))
    return var_matrix_ibin

def from_dict2bin_ells(dict_all_ = None, sim_ = None, var_name_=None, sorted_bin_sim=False): #get a vector storage in different keys bin/sim --> return the same in a matrix2D  [nbins, nells] for the respective bin
    #provide sim_ as used in the dictionary: 'sim0', 'sim1', ...
    for i, ibin in enumerate(dict_all_.keys()):
        if not sorted_bin_sim:
            d_ = dict_all_[sim_][ibin][var_name_]        
        else:
            d_ = dict_all_[ibin][sim_][var_name_]
        if not i: var_matrix_isim = d_
        else: var_matrix_isim = np.vstack(( var_matrix_isim, d_ ))
    return var_matrix_isim
    
#####################################################################################################
#####################################################################################################
#### AVERAGES over simulations/realizations
####
def average_over_realization_dict(dict_all_=None, var_name=None): #It will construct a dict with diff bins with matrices of realiza'coes
    for j, jbin in enumerate(dict_all_.keys()):
#        for i in range(dict_all_[jbin]['sim0']['nreals']):
        for i in range( len(dict_all_[jbin].keys()) ):
            if not i: cl_matrix_ibin = dict_all_[jbin]['sim{}'.format(i)][var_name]
            else: cl_matrix_ibin = np.vstack((cl_matrix_ibin, dict_all_[jbin]['sim{}'.format(i)][var_name]))
        if not j: 
            avg_over_reals_Cl   = {jbin: {'mean':np.mean(cl_matrix_ibin,axis=0), 'std':np.std(cl_matrix_ibin,axis=0),  'matrix':cl_matrix_ibin}}
        else: 
            avg_over_reals_Cl[jbin] =    {'mean':np.mean(cl_matrix_ibin,axis=0), 'std':np.std(cl_matrix_ibin,axis=0),  'matrix':cl_matrix_ibin}
        del cl_matrix_ibin
    return avg_over_reals_Cl
    
def dict_averages(dict_all_=None, cl_names_ = ['cl_kappa_sim', 'cl_kappa_uncorrelated', 'cl_cross_sim', 'cl_hi_sim']): 
    for i,icl_name in enumerate(cl_names_):
        if not i: avg_dict = {icl_name: dcopy(average_over_realization_dict(dict_all_=dict_all_, var_name=icl_name))}
        else:     avg_dict[icl_name] =  dcopy(average_over_realization_dict(dict_all_=dict_all_, var_name=icl_name))
    return avg_dict

def get_mean_var_from_loaded_dict(dict_=None, var_=None, bin_=None):
    #bin_ = int(int(bin_[-1]))
    bin_ = int(bin_.split('bin')[-1])
    if ('cl' in var_) or ('cross' in var_):
        for i,isim in enumerate(dict_.keys()):
            if not i:
                l_   = dict_[isim][var_][0,:].astype(int)
                cl_  = dict_[isim][var_][bin_,:]
            else:
                cl_  = np.vstack((cl_,dict_[isim][var_][bin_,:]))
        return {'l':l_, 'mean':np.mean(cl_,axis=0), 'std': np.std(cl_,axis=0), 'matrix':cl_}
    elif 'alm' in var_:
        for i,isim in enumerate(dict_.keys()):
            if not i:
                l_   = dict_[isim][var_][0,:].real.astype(int)
                m_   = dict_[isim][var_][1,:].real.astype(int)
                cl_  = dict_[isim][var_][int(bin_+1),:]
            else:
                cl_  = np.vstack((cl_,dict_[isim][var_][bin_,:]))                
        return {'l':l_, 'm':m_, 'mean':np.mean(cl_,axis=0), 'std': np.std(cl_,axis=0), 'matrix':cl_}
    else:
        raise Exception()

def dict_averages_from_loaded_data(dict_all_=None, varnames_=['cl_hi_sim', 'cl_kappa_sim', 'cl_cross_sim', 'cl_kappa_uncorrelated','cross_correlation_coef']):
    for i,var_ in enumerate(varnames_):
        nch,nl = dict_all_['sim0'][var_].shape
        if ('cl' in var_) or ('cross' in var_):
            nch=nch-1 #first column is l
        elif 'alm' in var_:
            nch=nch-2 #first column is l, second is m            
        else:
            raise Exception
        for j in range(nch):
            bin_ = 'bin{}'.format(j+1)
            if not j: 
                d_mv ={bin_: get_mean_var_from_loaded_dict(dict_=dict_all_, var_=var_, bin_=bin_)}      
            else:
                d_mv[bin_]= get_mean_var_from_loaded_dict(dict_=dict_all_, var_=var_, bin_=bin_)
        if not i:
            dict_all_mean_ = {var_: d_mv}
        else:
            dict_all_mean_[var_] = d_mv
    return dict_all_mean_
    
    
#####################################################################################################
#####################################################################################################
####Binning cls
####   
    
def clsbinned(cls=None,del_l=1,lmin=None, lmax=None, include_leftover=False):
    lvec = np.arange(cls.size)
    lmin = lvec.min() if type(lmin)==type(None) else max(lmin, lvec.min()) 
    lmax = lvec.max() if type(lmax)==type(None) else min(lmax, lvec.max()) 
    lnew = np.arange(lmin, lmax, del_l)
    if include_leftover and lnew.max()<lmax: lnew = np.hstack((lnew, lmax))
    for i in range(lnew.size-1):
        inds = np.where((lvec>=lnew[i])*(lvec<lnew[i+1]))[0]
        if not i: 
            l_weight  = (lnew[i]+lnew[i+1])/2
            cl_weight = np.sum( (2*lvec[inds]+1)*cls[inds])/np.sum(2*lvec[inds]+1)
            #weight    = np.sum(lvec[inds],axis=0)
        else:     
            l_weight  = np.hstack(( l_weight, (lnew[i]+lnew[i+1])/2 ))
            cl_weight = np.hstack(( cl_weight, np.sum( (2*lvec[inds]+1)*cls[inds])/np.sum(2*lvec[inds]+1) ))
    return l_weight, cl_weight
    
#####################################################################################################
#####################################################################################################
#### DIRECTORY VERIFICATION
####
def verification_dir(dirname=None, path=None, clear=True, verbose=False):
    import os, shutil
    pathout = os.path.join(path,dirname)
    if not os.path.isdir(path):
        os.mkdir(path)
    else:
        pass
    if not os.path.isdir(pathout):
        #if verbose:
        #    print("Directory "+ dirname + " does not exist.")
        os.makedirs(pathout)
        #if verbose:
        #    print("Directory "+ dirname + " created.")
    else:
        #if verbose:
        #    print("Directory "+ dirname + " exists.")
        if clear: 
            shutil.rmtree(pathout)
            #if verbose:
            #    print("Directory "+ dirname + " cleaned.")
            os.makedirs(pathout)
        else:
            pass
    return pathout    
    
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
    
    

def savedata_from_dict(Cl_=None,  filename=None, path=None, prefix=None, suffix=None, header= "[1] l, [>2] Cl [multipoles, redshift/frequency bins]", 
                       clear=True, alm_data = None, verbose=False): #alm_data =None, 'real', or 'imag'
    import os
    nu,npix = Cl_.shape #[l, CL]
    if not (prefix==None): filename = "_".join((  prefix  , filename  ))
    if not (suffix==None): filename = "_".join((  filename, suffix    ))     
    pathname = verification_dir(dirname = path.split('/')[-1], 
                                path    = '/'.join(( path.split('/')[:-1] )), 
                                clear   = clear)
    pathname = os.path.join(pathname,filename)
    pathname = ".".join(( pathname, 'txt' ))
    if   alm_data==None:
        np.savetxt(pathname, Cl_.T, fmt=['%d']+["%e"]*int(nu-1), delimiter=" ", header=header)
    elif alm_data.lower() in ['complex']: 
        np.savetxt(pathname, Cl_.T, fmt='%s', delimiter=" ", header=header)
    elif alm_data.lower() in ['real']: 
        np.savetxt(pathname, Cl_.T, fmt=['%d']*2 + ["%e"]*int(nu-2), delimiter=" ", header=header)
    elif alm_data.lower() in ['imag', 'imaginary']: 
        np.savetxt(pathname, Cl_.T, fmt=['%d']*2 + ["%e"]*int(nu-2), delimiter=" ", header=header)     
    else:
        raise Exception ('Name no valid: {}'.format(alm_data.lower()))
    return pathname

def save_matrix(matrix_=None,  filename=None, path=None, prefix=None, suffix=None, header= "[1>] WFG  [redshift/frequency bins, redshift/frequency bins]", clear=True, verbose=False): #alm_data =None, 'real', or 'imag'
    import os
    nu = matrix_.shape[0] #[l, CL]
    if not (prefix==None): filename = "_".join((  prefix  , filename  ))
    if not (suffix==None): filename = "_".join((  filename, suffix    ))     
    pathname = verification_dir(dirname = path.split('/')[-1], 
                                path    = '/'.join(( path.split('/')[:-1] )), 
                                clear   = clear)
    pathname = os.path.join(pathname,filename)
    pathname = ".".join(( pathname, 'txt' ))
    np.savetxt(pathname, matrix_.T, fmt=["%e"]*int(nu), delimiter=" ", header=header)
    return pathname



######################################################################################################
### CEMITERY  ########################################################################################
######################################################################################################
def cmbwl_x_hi_sims_old(Chi_dict=None, Ckk_vec=None, Cross_dict=None, bin_=None, l_=None,  jseed_kunc=None, jseed_hi=None,fact=1):
    import healpy as hp
    #####################################
    Chi        = Chi_dict[bin_]
    Cross      = fact*Cross_dict[bin_]
    Ckk_uncorr = Ckk_vec - Cross*Cross/Chi
    ###########################################
    # [1] Cl(k-unc) --> (s1)alm(k-unc): s1=seed1
    # [2] Cl(hi)    --> (s2)alm(hi)   : s2=seed2    
    # [3] Cl(N)     --> (s3)alm(N)    : s3=seed3    
    # [4] alm(k-sim) = (s1)alm(hi-weigthed) + (s2)alm(k-unc) +(s3)alm(N)
    ##########################################
    # [1] ---- 
    np.random.seed(jseed_kunc)
    auncorr_lm  = hp.synalm(Ckk_uncorr, lmax= int(l_[-1]))
    cl_k_uncorr = dcopy(Ckk_uncorr)#hp.alm2cl(auncorr_lm)
    # [2] ---- 
    np.random.seed(jseed_hi)
    ahi_lm       = hp.synalm(Chi , lmax= int(l_[-1]))
    _phase_      = np.angle(ahi_lm , deg=False)
    _amplitude_  = np.abs(  hp.synalm(Cross*Cross/Chi, lmax= int(l_[-1])))
    import cmath
    ahi_weighted_lm = np.array([cmath.rect(ir, iphi) for ir,iphi in zip(_amplitude_,_phase_)])
    # [4] ---
    akk_lm = ahi_weighted_lm  + auncorr_lm# + akknoise_lm
    
    #cl_k_uncorr = hp.alm2cl(auncorr_lm)
    cl_k_sim        = hp.alm2cl(akk_lm       , lmax_out=l_.max())#[l_.min():l_.max()+1]
    cl_cross_sim    = hp.alm2cl(akk_lm,ahi_lm, lmax_out=l_.max())#[l_.min():l_.maxnchannels()+1]
    cl_hi_sim       = hp.alm2cl(ahi_lm       , lmax_out=l_.max())#[l_.min():l_.max()+1]
    
    return {'cl_kappa_sim' :cl_k_sim, 'cl_kappa_uncorrelated':cl_k_uncorr, 'cl_cross_sim':cl_cross_sim,'cl_hi_sim':cl_hi_sim,
            'alm_kappa_sim':akk_lm,  'alm_kappa_uncorrelated':auncorr_lm , 'alm_hi_sim': ahi_lm, 'alm_hi_weighted': ahi_weighted_lm,
            'cross_correlation_coef': crosscorr_coeff(cl_k_sim, cl_hi_sim, cl_cross_sim),
            'bin':bin_, 'seek_kunc':jseed_kunc, 'seed_hi':jseed_hi, 'l':l_}#np.arange(cl_k_sim.size)}
######################
######################
def cross_simulations_correlated_old(clf1_dict_=None, clf2_dict_=None, clcx_dict_=None, seed_k=None, seed_hi=None, show_time=False):
    import cmath
    tax=1000 #to distinguish the alm_hi
    if show_time: itime = time.time() 
    l             = clf1_dict_['l']
    alm_hi_dict   = {'l': dcopy(l)}
    alm_k_u_dict  = {'l': dcopy(l)}
    alpha_l_ratio = 0
    vec_seeds_hi = np.array([])

    sim_names = dcopy(clf1_dict_)
    del sim_names['l']
    sim_names = np.array(dcopy(list(sim_names.keys())))

    for j,jbin in enumerate(sim_names):
        jseed_hi        = seed_hi + j*tax
        vec_seeds_hi    = np.hstack((vec_seeds_hi, jseed_hi))
        np.random.seed(jseed_hi)
        alm_hi_jbin      = hp.synalm(clf1_dict_[jbin], lmax= int(l.max()))
        _phase_jbin_     = np.angle(alm_hi_jbin , deg=False)
        cl_ratio         = clcx_dict_[jbin]*clcx_dict_[jbin]/clf1_dict_[jbin]
        _amplitude_jbin_ = np.abs(  hp.synalm( cls=cl_ratio, lmax=int(l.max()) )  )
        ############
        alm_k_uj           = np.array([cmath.rect(ir, iphi) for ir,iphi in zip(_amplitude_jbin_,_phase_jbin_)])        
        alm_hi_dict[jbin]  = dcopy(alm_hi_jbin)
        alm_k_u_dict[jbin] = dcopy(alm_k_uj)
        alpha_l_ratio      += cl_ratio
    
    np.random.seed(seed_k)
    alm_k_corr  = hp.synalm(clf2_dict_['bin1']-alpha_l_ratio, lmax=int(l.max()) )

    sim_names = dcopy(clf1_dict_)
    del sim_names['l']
    sim_names = np.array(dcopy(list(sim_names.keys())))
    
    alm_k_sim = dcopy(alm_k_corr)
    for j,jbin in enumerate(sim_names):
        alm_k_sim += alm_k_u_dict[jbin]
        if not j: 
            alm_hi_sim    = dcopy(alm_hi_dict[jbin])
            alm_k_uncorr  = dcopy(alm_k_u_dict[jbin])
            cl_hi_sim     = dcopy(hp.alm2cl(alm_hi_dict[jbin], lmax=int(l.max())))
            cl_hi_theory  = dcopy(clf1_dict_[jbin])
        else:     
            alm_hi_sim   = np.vstack(( alm_hi_sim  , dcopy(alm_hi_dict[jbin])  ))
            alm_k_uncorr = np.vstack(( alm_k_uncorr, dcopy(alm_k_u_dict[jbin]) ))
            cl_hi_sim    = np.vstack(( cl_hi_sim   , dcopy(hp.alm2cl(alm_hi_dict[jbin], lmax=int(l.max()))) ))
            cl_hi_theory = np.vstack(( cl_hi_theory, dcopy(clf1_dict_[jbin])    ))

    sim_names = dcopy(clcx_dict_)
    del sim_names['l']
    sim_names = np.array(dcopy(list(sim_names.keys())))

    
    for j,jbin in enumerate(sim_names):
        if not j:
            cl_cx_sim    = hp.alm2cl(alm_k_sim, alm_hi_sim[j], lmax=int(l.max()) )
            cl_cx_theory = dcopy(clcx_dict_[jbin])
        else:
            cl_cx_sim    = np.vstack(( cl_cx_sim   , hp.alm2cl(alm_k_sim, alm_hi_sim[j], lmax=int(l.max()) ) ))
            cl_cx_theory = np.vstack(( cl_cx_theory, dcopy(clcx_dict_[jbin]) ))
    
    cl_k_sim  = hp.alm2cl(alm_k_sim,  lmax=int(l.max()) )
    if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))  
    return {'alm_kappa_sim': alm_k_sim , 'alm_hi_sim': alm_hi_sim, 
            'alm_kappa_uncorr':alm_k_uncorr, 'alm_kappa_corr':alm_k_corr,
            'cl_kappa_sim':  cl_k_sim  , 'cl_kappa_theory': clf2_dict_['bin1'], 
            'cl_hi_sim'   :  cl_hi_sim , 'cl_hi_theory'   : cl_hi_theory, 
            'cl_cross_sim':  cl_cx_sim , 'cl_cross_theory': cl_cx_theory, 
            'l':l, 'seed-hi':np.int32(vec_seeds_hi), 'seed-kunc':np.int32(seed_k)
            }
######################
######################
def stat_simulations_old1(dict_sims_=None,show_time=False):
    if show_time: itime = time.time()
    N = len(dict_sims_.keys())
    cl_k_mean  = np.zeros_like(dict_sims_['sim0']['cl_kappa_sim'])
    cl_hi_mean = np.zeros_like(dict_sims_['sim0']['cl_hi_sim'])
    cl_cx_mean = np.zeros_like(dict_sims_['sim0']['cl_cross_sim'])
    #####
    for j, jsim in enumerate(dict_sims_.keys()):
        cl_k_mean  += dict_sims_[jsim]['cl_kappa_sim']/N
        cl_hi_mean += dict_sims_[jsim]['cl_hi_sim']/N
        cl_cx_mean += dict_sims_[jsim]['cl_cross_sim']/N
    
    cl_k_std  = np.zeros_like(cl_k_mean )
    cl_hi_std = np.zeros_like(cl_hi_mean)
    cl_cx_std = np.zeros_like(cl_cx_mean)
    for j, jsim in enumerate(dict_sims_.keys()):
        cl_k_std  += (dict_sims_[jsim]['cl_kappa_sim'] - cl_k_mean )**2/N
        cl_hi_std += (dict_sims_[jsim]['cl_hi_sim'   ] - cl_hi_mean)**2/N
        cl_cx_std += (dict_sims_[jsim]['cl_cross_sim'] - cl_cx_mean)**2/N
    cl_k_std  = np.sqrt(cl_k_std )    
    cl_hi_std = np.sqrt(cl_hi_std)    
    cl_cx_std = np.sqrt(cl_cx_std)    
    if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))                
    return {'cl_kappa_sim': {'mean': cl_k_mean ,
                             'std' : cl_k_std
                            },
            'cl_hi_sim'   : {'mean': cl_hi_mean ,
                             'std' : cl_hi_std
                            },
            'cl_cross_sim': {'mean': cl_cx_mean ,
                             'std' : cl_cx_std
                            }
           }
######################
######################
def stat_simulations_old2(dict_sims_=None, lmin=None, lmax=None, del_l=None, ich=None, use_dl=True, include_leftover=False, show_time=False):
    if show_time: itime = time.time()
    N = len(dict_sims_.keys())
    if int(del_l)==1:            
        cl_k_mean  = np.zeros_like(dict_sims_['sim0']['cl_kappa_sim'])
        cl_hi_mean = np.zeros_like(dict_sims_['sim0']['cl_hi_sim'   ][ich])
        cl_cx_mean = np.zeros_like(dict_sims_['sim0']['cl_cross_sim'][ich])
        #####
        #l_   = dict_sims_['sim0']['l']
        l_   = np.arange(dict_sims_['sim0']['l'].size)
        fact = fact=l_*(l_+1)/2/np.pi if use_dl else 1        
        #####
        for j, jsim in enumerate(dict_sims_.keys()):
            cl_k_mean  += fact*dict_sims_[jsim]['cl_kappa_sim']/N
            cl_hi_mean += fact*dict_sims_[jsim]['cl_hi_sim'   ][ich]/N
            cl_cx_mean += fact*dict_sims_[jsim]['cl_cross_sim'][ich]/N        
        
        cl_k_std  = np.zeros_like(cl_k_mean )
        cl_hi_std = np.zeros_like(cl_hi_mean)
        cl_cx_std = np.zeros_like(cl_cx_mean)
        for j, jsim in enumerate(dict_sims_.keys()):
            cl_k_std  += ( fact*dict_sims_[jsim]['cl_kappa_sim']      - cl_k_mean )**2/N
            cl_hi_std += ( fact*dict_sims_[jsim]['cl_hi_sim'   ][ich] - cl_hi_mean)**2/N
            cl_cx_std += ( fact*dict_sims_[jsim]['cl_cross_sim'][ich] - cl_cx_mean)**2/N   
            
    elif int(del_l)>=1:
        l_ ,cl_k_mean = cxfs.clsbinned(cls=dict_sims_['sim0']['cl_kappa_sim']     , del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)
        cl_hi_mean    = cxfs.clsbinned(cls=dict_sims_['sim0']['cl_hi_sim'   ][ich], del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]
        cl_cx_mean    = cxfs.clsbinned(cls=dict_sims_['sim0']['cl_cross_sim'][ich], del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]

        cl_k_mean  = np.zeros_like(cl_k_mean)
        cl_hi_mean = np.zeros_like(cl_hi_mean)
        cl_cx_mean = np.zeros_like(cl_cx_mean)      
        #####
        fact = fact=l_*(l_+1)/2/np.pi if use_dl else 1
        #####
        for j, jsim in enumerate(dict_sims_.keys()):
            cl_k_mean  += fact*cxfs.clsbinned(cls=dict_sims_[jsim]['cl_kappa_sim']     , del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]/N
            cl_hi_mean += fact*cxfs.clsbinned(cls=dict_sims_[jsim]['cl_hi_sim'   ][ich], del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]/N
            cl_cx_mean += fact*cxfs.clsbinned(cls=dict_sims_[jsim]['cl_cross_sim'][ich], del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1]/N
        
        cl_k_std  = np.zeros_like(cl_k_mean )
        cl_hi_std = np.zeros_like(cl_hi_mean)
        cl_cx_std = np.zeros_like(cl_cx_mean)
        for j, jsim in enumerate(dict_sims_.keys()):
            cl_k_std  += ( fact*cxfs.clsbinned(cls=dict_sims_[jsim]['cl_kappa_sim']     , del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1] - cl_k_mean )**2/N
            cl_hi_std += ( fact*cxfs.clsbinned(cls=dict_sims_[jsim]['cl_hi_sim'   ][ich], del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1] - cl_hi_mean)**2/N
            cl_cx_std += ( fact*cxfs.clsbinned(cls=dict_sims_[jsim]['cl_cross_sim'][ich], del_l=del_l, lmin=lmin, lmax=lmax, include_leftover=include_leftover)[1] - cl_cx_mean)**2/N
    else:
        raise Exception
    
    cl_k_std  = np.sqrt(cl_k_std )    
    cl_hi_std = np.sqrt(cl_hi_std)    
    cl_cx_std = np.sqrt(cl_cx_std)    
    if show_time: print('time: {0:.4f} seg'.format(time.time()-itime))                
    return {'cl_kappa_sim': {'mean': cl_k_mean ,
                             'std' : cl_k_std
                            },
            'cl_hi_sim'   : {'mean': cl_hi_mean ,
                             'std' : cl_hi_std
                            },
            'cl_cross_sim': {'mean': cl_cx_mean ,
                             'std' : cl_cx_std
                            },
            'l':l_
           }
