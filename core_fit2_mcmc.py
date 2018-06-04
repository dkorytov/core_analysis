import numpy as np
from core_fit2_util import * 
import emcee
import pickle


def lnprior(theta):
    dis_len, merg_len, infall_m = theta
    if(-3.0 < dis_len < 1.0  and 0.0 < merg_len < 1.0 and 10.0 < infall_m < 15):
        return 0.0
    else:
        return -np.inf

my_zmr_sdss= None
my_zmr_valid = None
my_clstrs = None

def set_my(a,b,c):
    global my_zmr_sdss,my_zmr_valid,my_clstrs
    my_zmr_sdss= a
    my_zmr_valid = b
    my_clstrs = c



def lnlike(theta,zmr_sdss,zmr_valid,clstrs):
    dis_len, merg_len, infall_m = theta
    zmr_core = zmr_from_clusters(dis_len,merg_len,clstrs,zmr_valid,infall_m)
    cost = calc_gal_density_cost2(zmr_core,zmr_sdss)
    return -cost

def lnprob(theta,zmr_sdss,zmr_valid,clstrs):
    lnp = lnprior(theta)
    if not np.isfinite(lnp):
        return -np.inf
    return lnp + lnlike(theta,zmr_sdss,zmr_valid,clstrs)

def lnlike2(theta):
    dis_len, merg_len, infall_m = theta
    zmr_core = zmr_from_clusters(dis_len,merg_len,my_clstrs,my_zmr_valid,infall_m)
    cost = calc_gal_density_cost2(zmr_core,my_zmr_sdss)
    return -cost

def lnprob2(theta):
    lnp = lnprior(theta)
    if not np.isfinite(lnp):
        return -np.inf
    return lnp + lnlike2(theta)



def run_mcmc(theta_0,nwalkers,niter,zmr_sdss,zmr_valid,clstrs):
    dis_len, merg_len, infall_m = theta_0
    ndim = len(theta_0)
    pos = [theta_0 + 1e-2*np.random.randn(ndim) for i in range(nwalkers)]
    set_my(zmr_sdss,zmr_valid,clstrs)
    sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob2,threads=24)
    #args=(zmr_sdss,zmr_valid,clstrs)
    #s1 = pickle.dumps(zmr_sdss)

    #s2 = pickle.dumps(zmr_valid)
    #s3 = pickle.dumps(clstrs)
    sampler.run_mcmc(pos,niter)
    x = range(niter)
    plt.figure()
    plt.ylabel('dis_len')
    for i in range(nwalkers):
        plt.plot(x,sampler.chain[i,:,0],'b')
    plt.figure()
    plt.ylabel('meg_len')
    for i in range(nwalkers):
        plt.plot(x,sampler.chain[i,:,1],'b')
    plt.figure()
    plt.ylabel('infall_m')
    for i in range(nwalkers):
        plt.plot(x,sampler.chain[i,:,2],'b')
    return sampler
