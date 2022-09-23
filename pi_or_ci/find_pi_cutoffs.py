import numpy as np
from pi_or_ci.utils import name_to_vars,vars_to_name,log

def add_transitions_if_missing(unique_temps,ions,ret_rhos,PI_or_CI):
    for i,t in enumerate(unique_temps):
        for j,ion in enumerate(ions):
            if i==0:
                continue
            if PI_or_CI[i,j] > 0 and PI_or_CI[i-1,j] < 0:
                ret_rhos[i,j] = 1e-13
            if PI_or_CI[i,j] > 0 and PI_or_CI[i-1,j] > 0:
                ret_rhos[i,j] = np.nan
    return ret_rhos      

def find_PI_cutoff(data,threshold = .05,loud = False,min_stable_num=5,negligible = 1e-8,force_convergence = False):
    name,ions,temps,densities,Us,ion_fractions = data
    atom,redshift,radfield = name_to_vars(name)
    unique_temps = np.unique(temps)
    ntemps = len(unique_temps)
    ret_ts = np.zeros(ntemps)
    ret_rhos = np.zeros((ntemps,len(ions)))
    tiny = np.zeros((ntemps,len(ions)),dtype = bool)
    PI_or_CI = np.zeros((ntemps,len(ions))) #-1 = PI, 0 = transition, +1 = CI
    for i,t in enumerate(unique_temps):
        log("T=10**%.2f,z=%.2f"%(np.log10(t),redshift),loud)
        ret_ts[i] = t
        d = densities[temps==t]
        f_i = ion_fractions[temps==t]
        for j,ion in enumerate(ions):
            log_frac = np.log10(f_i[:,j])
            if np.amax(log_frac)<np.log10(negligible):
                tiny[i,j] = True
                log('(negligible) ',loud,end='')
            last_10_average = np.mean(log_frac[-min_stable_num:])
            diffs = np.diff(log_frac)
            relative_stable_max = 10.**last_10_average*(1+threshold)
            relative_stable_min = 10.**last_10_average*(1-threshold)
            relative_stable = np.logical_and(10.**log_frac>relative_stable_min,
                                             10.**log_frac<relative_stable_max)
            if ion.split()[1]=='I':
                ret_rhos[i,j] = np.nan
                nans_or_negligible = np.logical_or(np.isinf(log_frac),log_frac<np.log10(negligible))
                if np.all(relative_stable[-min_stable_num:]) or np.all(nans_or_negligible[-min_stable_num:]):
                    PI_or_CI[i,j] = 1
                    if i>0 and PI_or_CI[i-1,j]==-1 and not tiny[i,j]:
                        log('no transition found, but %s switched from PI to CI at t=%.2e'%(ion,t),loud)
                else:
                    PI_or_CI[i,j] = -1
                    if i>0 and PI_or_CI[i-1,j]==1 and not tiny[i,j] and False:
                        PI_or_CI[:i,j] = -1
                        log('transitioned %s back to PI from CI at t=%.2e, assuming error and overwriting'%(ion,t),loud)
                ietype = 'PIE' if PI_or_CI[i,j] == -1 else 'CIE'
                log('%s is not ionized, but it appears like in %s'%(ion,ietype),loud)
            else:
                maximum = None if not np.any(diffs<-.01) else np.argmax(diffs<0)
                if maximum is None:
                    #no maximum means CI
                    #check if previous was PI
                    if i>0 and PI_or_CI[i-1,j] == -1 and not tiny[i,j]:
                        log('no transition found, but %s switched from PI to CI at t=%.2e'%(ion,t),loud)
                    log('%s all CI'%ion,loud)
                    ret_rhos[i,j] = np.nan
                    PI_or_CI[i,j] = 1
                elif maximum is not None and np.all(relative_stable[-min_stable_num:]):
                    #starts w/PI, goes CI
                    if log_frac[maximum]<=last_10_average+np.log10(2):
                        #PI contribution never bigger than CI
                        if i>0 and PI_or_CI[i-1,j]<0 and not tiny[i,j]:
                            log('no transition found, but %s switched from PI to CI at t=%.2e'%(ion,t),loud)
                        ret_rhos[i,j] = np.nan
                        PI_or_CI[i,j] = 1
                        log('%s has small PI spike at max %.2e, but still CI dominated'%(ion,d[np.argmax(relative_stable)]),loud)
                    else:
                        #PI contribution bigger than CI for some range
                        ret_rhos[i,j] = 10.**np.interp(-(last_10_average+np.log10(2)),-log_frac[maximum:],np.log10(d[maximum:]))
                        PI_or_CI[i,j] = 0
                        log('%s starts PI, goes CI at %.2e'%(ion,ret_rhos[i,j]),loud)
                else:
                    #all PI
                    if i>0 and PI_or_CI[i-1,j] == 1 and not tiny[i,j] and False:
                        ret_rhos[i,j] = 1e-12
                        PI_or_CI[:i,j] = -1
                        log('transitioned %s back to PI from CI at t=%.2e, assuming error and overwriting'%(ion,t),loud)
                    log('%s all PI'%ion,loud)
                    ret_rhos[i,j] = np.nan
                    PI_or_CI[i,j] = -1 
            if ion.split()[1]=='II':
                if PI_or_CI[i,j] == -1 or PI_or_CI[i,j-1] == -1:
                    log('If one of %s and %s is PI, they must both be'%(ion[:-1],ion),loud)
                    if (not tiny[i,j] and not tiny[i,j-1]):
                        PI_or_CI[i,j-1] = -1
                        PI_or_CI[i,j] = -1
                    else:
                        PI_or_CI[i,j-1] = -1
                        PI_or_CI[i,j] = -1
    PI_or_CI[tiny]=1.75
    ret_rhos = add_transitions_if_missing(unique_temps,ions,ret_rhos,PI_or_CI)
    return name,ions,ret_ts,ret_rhos,PI_or_CI      