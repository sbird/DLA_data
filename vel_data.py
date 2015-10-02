# -*- coding: utf-8 -*-
"""Plot the metal line data from Prochaska"""

import numpy as np
import matplotlib.pyplot as plt
import leastsq as ls
import os.path as path

datadir = path.dirname(__file__)

def load_data(zrange=None):
    """Load the data tables from Neeleman 2013"""
    data2 = np.loadtxt(path.join(datadir,"apj469315t2_mrt_mod.txt"))
    redshift = data2[:,0]
    met = data2[:,3]
    vel_data = data2[:,5]
    if zrange != None:
        ind = np.where((redshift < zrange[0])*(redshift > zrange[1]))
        return (redshift[ind], met[ind], vel_data[ind])
    return (redshift, met, vel_data)

def _bootstrap_sample(v_table, vel_data, error):
    """Generate a Monte Carlo error sample of the differential distribution."""
    # Generate some Monte Carlo samples where each element is perturbed by
    # a Gaussian, sigma given by error.
    index = np.random.random_integers(0, np.size(vel_data)-1, np.size(vel_data))
    if error > 0.:
        errors = vel_data + np.random.normal(0,error,size=np.size(vel_data))
        bootstrap = errors[index]
    else:
        bootstrap = vel_data[index]
    nn = np.histogram(bootstrap,v_table)[0]
    return nn

def pdf_with_error(vel_data, v_table, lognorm=True, poisson=False, cumulative=False):
    """Plot a pdf of vel_data in v_table bins, with Poisson error bars"""
    center = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
    nn = np.histogram(vel_data,v_table)[0]
    #Normalise so that integral in log space is unity.
    #This emulates np.histogram(np.log10(met), np.log10(bin),density=True)
    if not cumulative:
        if lognorm:
            norm = np.array([(-np.log10(v_table[i])+np.log10(v_table[i+1])) for i in range(np.size(v_table)-1)])
        else:
            norm = np.array([(-v_table[i]+v_table[i+1]) for i in range(np.size(v_table)-1)])
        norm *= np.size(vel_data)
        vels = nn / norm
        if poisson:
            #Use poisson errors
            verr = (np.sqrt(nn))/norm
            #These v_table will have error bars that go to 0
            ind = np.where(nn == 1)
            verr[ind] *= (0.98)
            return (center, vels, verr)
    else:
        vels = np.cumsum(nn)
    return (center, vels)

def plot_prochaska_2008_data(zrange = None, nv_table=11):
    """Plot a velocity width histogram from Neeleman 2013"""
    (_, _, vel_data) = load_data(zrange)
    v_table=np.logspace(np.log10(10),np.log10(np.max(vel_data)+10),nv_table)
    #Bins should not be smaller than error on v90, which is 10.
    (center, vels) = pdf_with_error(vel_data, v_table)
#     plt.semilogx(center, vels,'o', color="black")
    plt.errorbar(center,vels,xerr=[center-v_table[:-1],v_table[1:]-center],fmt='o', color="black")
    plt.xlim(10, 1000)
    return (center, vels)

def plot_cum_vw_data(zrange=None):
    """Plot a cumulative velocity width distribution with error"""
    (_, _, data) = load_data(zrange)
    v_table = 10**np.arange(1, np.log10(np.max(data)+10), 0.1)
    (center, cvels) = pdf_with_error(data, v_table, cumulative=True)
    plt.semilogx(center, cvels, ls="-",color="black", label="N13")
    return (center, cvels)

def plot_prochaska_2008_correlation(zrange = None, color="black"):
    """Plot the observed correlation between velocity widths and metallicity from Neeleman 2013"""
    (_, met, vel_data) = load_data(zrange)
    plt.loglog(vel_data, 10**met,'o',color=color)
    vel = np.log10(vel_data)
    (intercept, slope, _) = ls.leastsq(vel,met)
#     print "obs fit: ",intercept, slope, np.sqrt(var)
#     print "obs correlation: ",ls.pearson(vel, met,intercept, slope)
#     print "obs kstest: ",ls.kstest(vel, met,intercept, slope)
    xx = np.logspace(np.min(vel), np.max(vel),15)
    plt.loglog(xx, 10**intercept*xx**slope, color=color)

def plot_extra_stat_hist(stat=False,zrange=None, nv_table=11):
    """Plot a histogram of the mean-median statistic"""
    data2 = np.loadtxt(path.join(datadir,"apj469315t2_mrt_mod.txt"))
    redshift = data2[:,0] #np.concatenate([data[:,0], data2[:,0]])
    fmm = data2[:,6+int(stat)] #np.concatenate([data[:,1],data2[:,3]])
    if zrange != None:
        ind = np.where((redshift < zrange[0])*(redshift > zrange[1]))
        fmm = fmm[ind]

    v_table=np.linspace(0,1,nv_table)
    (center, vels) = pdf_with_error(fmm, v_table, lognorm=False, cumulative=False)
#     plt.fill_between(center, lower, upper,color="black", alpha=0.3)
    plt.errorbar(center,vels,xerr=[center-v_table[:-1],v_table[1:]-center],fmt='o', color="black")
#     plt.plot(center, vels, 'o', color="black")

def plot_cum_stat_data(stat=False, zrange=None):
    """Plot a cumulative velocity width distribution with error"""
    data2 = np.loadtxt(path.join(datadir,"apj469315t2_mrt_mod.txt"))
    redshift = data2[:,0]
    fmm = data2[:,6+int(stat)]
    if zrange != None:
        ind = np.where((redshift < zrange[0])*(redshift > zrange[1]))
        fmm = fmm[ind]

    v_table = np.arange(0, 1, 0.02)
    (vbin, cfmm) = pdf_with_error(fmm, v_table, lognorm=False, cumulative=True)
    plt.plot(vbin, cfmm, ls="-",color="black", label="N13")
    return (vbin, cfmm)

def load_metal_data(zrange=None):
    """Load the data tables from Rafelski 2012"""
    data2 = np.loadtxt(path.join(datadir,"apj433829t2_mrt.txt"), usecols=(1,11))
    data3 = np.loadtxt(path.join(datadir,"apj433829t3_mrt.txt"),usecols=(1,11))
    redshift = np.concatenate((data2[:,0], data3[:,0]))
    met = np.concatenate((data2[:,1], data3[:,1]))
    if zrange != None:
        ind = np.where((redshift < zrange[0])*(redshift > zrange[1]))
        redshift = redshift[ind]
        met = met[ind]

    return (redshift, met)

def plot_alpha_metal_data(zrange=(3.5,2.5),nv_table=6):
    """
       Plot the metallicities of DLAs from the catalogues found in Rafelski 2012
       redshift  = (4,3) will show only quasars between z=4 and z=3
    """
    (_, met) = load_metal_data(zrange)
    v_table=np.linspace(np.min(met),np.max(met),nv_table)
    #Bins should not be smaller than error on metallicity, which is about 0.2.
    if (- np.min(met) + np.max(met))/ nv_table < 0.2:
        raise ValueError("Requested bins smaller than 0.2, error on metallicity")

    (center, vels, verr) = pdf_with_error(met, v_table, lognorm=False, poisson=True)

    plt.errorbar(center,vels,xerr=[center-v_table[:-1],v_table[1:]-center],yerr=verr,fmt='o', color="black")

    return (center, vels, verr)

def plot_lls_metal_data(nv_table=7):
    """
       Plot the metallicities of DLAs from the catalogues found in Lehner 2013
       at z=0-1
    """
    data = np.loadtxt(path.join(datadir,"lls_metallicity.txt"), usecols=(1,3))
    #Only LLS
    ind = np.where(data[:,0] < 20.3)
    met = data[ind,1]
    v_table=np.linspace(np.min(met),np.max(met),nv_table)
    #Bins should not be smaller than error on metallicity, which is about 0.2.
    if (- np.min(met) + np.max(met))/ nv_table < 0.2:
        raise ValueError("Requested bins smaller than 0.2, error on metallicity")

    (center, vels, verr) = pdf_with_error(met, v_table, lognorm=False, poisson=True)

    plt.errorbar(center,vels,xerr=[center-v_table[:-1],v_table[1:]-center],yerr=verr,fmt='o', color="black")
    return (center, vels, verr)

def plot_si1526_eqw(zrange = None, nv_table = 11):
    """Plot a histogram of the SiII 1526 equivalent width distribution"""
    data2 = np.loadtxt(path.join(datadir,"apj469315t2_mrt_mod.txt"))
    redshift = data2[:,0]
    eqw = data2[:,8]
    ind = np.where(eqw > 0)
    eqw = eqw[ind]
    redshift = redshift[ind]
    if zrange != None:
        ind = np.where((redshift < zrange[0])*(redshift > zrange[1]))
        eqw = eqw[ind]

    v_table=np.linspace(np.log10(np.min(eqw)),np.log10(np.max(eqw)),nv_table)
    (center, bins) = pdf_with_error(np.log10(eqw), v_table, lognorm=False)

    plt.errorbar(center,bins,xerr=[center-v_table[:-1],v_table[1:]-center],fmt='o', color="black")
#     plt.plot(center, bins,"-", color="black")
#     plt.fill_between(center, lower, upper, color="black", alpha=0.3)
    plt.xlim(-1.5, 0.5)
    return (center, bins)

