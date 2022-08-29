from scipy.special import exp1 as e1
from numpy import exp
from pi_or_ci.utils import vars_to_name
import numpy as np

def get_coronal_approx(atom,e = None,tlims = 'default',tspacing = 'default'):
    energies = [0]
    for ion in e:
        myatom,ionization = ion.split(' ')
        if myatom == atom and ionization !='I':
            energies.append(e[ion])
    boltzmann = 8.617333262e-5 #eV/K
    e_H = 13.6 #eV
    B = 1e6
    eh = 13
    energy_diffs = np.diff(energies)

    def ratio(i,T):
        kT = (boltzmann*T)
        return B*kT**2*eh**2/energy_diffs[i]**4*exp(-2*energy_diffs[i]/kT)/e1(energy_diffs[i]/kT)
    def term(j,T):
        product = 1
        for i in range(j):
            product*=ratio(i,T)
        return product

    if tlims == 'default':
        tlims = (2,8.01)
    if tspacing == 'default':
        tspacing = .01
    corona_temps = 10.**np.arange(tlims[0],tlims[1],tspacing)
    all_terms = np.zeros((len(energies),len(corona_temps)))
    for j in range(0,len(energies)):
        for i,t in enumerate(corona_temps):
            all_terms[j,i] = term(j,t)

    corona_cie = all_terms
    corona_cie[np.isnan(corona_cie)] = 0
    return energies,corona_temps,corona_cie
    
def get_cloudy_cie(atom,redshift=0.0,radfield = 'HM12',d=None):
    name = vars_to_name(atom,redshift,radfield)
    name,ions,temps,densities,Us,ion_fractions = d[name]
    cloudy_temps = np.unique(temps)
    cloudy_cie = np.zeros((len(ions),len(cloudy_temps)))
    for rho_index in range(-1,0):
        for k,t in enumerate(cloudy_temps):
            right_ts = ion_fractions[temps == t]
            right_rhos = densities[temps == t]
            rho = right_rhos[rho_index]
            for j,ion in enumerate(ions):
                cloudy_cie[j,k] = right_ts[rho_index,j]
    return cloudy_temps,cloudy_cie
