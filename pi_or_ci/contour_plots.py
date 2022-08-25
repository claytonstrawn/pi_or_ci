# preprocessing step. Actual data generated by Santi from CLOUDY tables

import os
import roman
import numpy as np
import matplotlib.pyplot as plt
from pi_or_ci.utils import log,get_ions,cycle_colors,name_to_vars,vars_to_name
from pi_or_ci.cutoff_plots import get_cutoffs,plot_cutoffs



def get_contour_values(name,mechanism,ions,data=None,cutoffs=None,slope = 1.66):
    name,all_ions,temps,densities,Us,ion_fractions = data[name]
    name,all_ions,ret_ts,ret_rhos,PI_or_CI = cutoffs[name]
    atom,redshift,radfield = name_to_vars(name)
    to_return = []
    for i,ion in enumerate(ions):
        plot_rhos,plot_ts = get_cutoffs(ion,redshift,radfield=radfield,data=data,cutoffs=cutoffs)
        ion_index = all_ions.index(ion)
        unique_temps = np.unique(temps)
        unique_densities = np.unique(densities)
        low_temp = unique_temps[0]
        low_dens = unique_densities[0]
        low_t_ion_fractions = ion_fractions[temps==low_temp,ion_index][:]
        low_rho_ion_fractions = ion_fractions[:,ion_index][densities == low_dens]
        two_dim_ion_fracs = np.zeros((len(unique_temps),len(unique_densities)))
        two_dim_ion_fracs_both = np.zeros((len(unique_temps),len(unique_densities)))
        two_dim_ion_fracs_PI = np.zeros((len(unique_temps),len(unique_densities)))
        two_dim_ion_fracs_CI = np.zeros((len(unique_temps),len(unique_densities)))
        for j,t in enumerate(unique_temps):
            two_dim_ion_fracs_both[j] = ion_fractions[temps==t,ion_index][:]
            for k,rho in enumerate(unique_densities):
                equivalent_density_at_low = rho/((t/low_temp)**(1/slope))
                if equivalent_density_at_low<low_dens:
                    equivalent_temperature_at_low = t/((rho/low_dens)**(slope))
                    frac = np.interp(equivalent_temperature_at_low,unique_temps,low_rho_ion_fractions)
                else:
                    frac = np.interp(equivalent_density_at_low,unique_densities,low_t_ion_fractions)
                two_dim_ion_fracs_PI[j,k] = frac
            two_dim_ion_fracs_CI[j] = ion_fractions[temps==t,ion_index][-1]
        if mechanism == 'both':
            two_dim_ion_fracs = two_dim_ion_fracs_both
        elif mechanism == 'PI':
            two_dim_ion_fracs = two_dim_ion_fracs_PI
        elif mechanism == 'CI':
            two_dim_ion_fracs = two_dim_ion_fracs_CI
        elif mechanism == 'split':
            for j,t in enumerate(unique_temps):
                for k,rho in enumerate(unique_densities):
                    pixel_is_PI = (np.interp(rho,plot_rhos,plot_ts)>t)
                    if pixel_is_PI:
                        two_dim_ion_fracs[j,k] = two_dim_ion_fracs_PI[j,k]
                    else:
                        two_dim_ion_fracs[j,k] = two_dim_ion_fracs_CI[j,k]
        else:
            assert False,'available mechanisms are "both","PI","CI", and "split"'
        to_return.append(two_dim_ion_fracs)
    return to_return
    
def plot_contours(ions,redshift,radfield = 'HM12',ax=None, show_cutoffs= False,
                  levels = 'default',data=None,cutoffs=None,mechanism = 'both',
                  colors = None,cutoff_linewidth = None,fig=None,show_values = False):
    if isinstance(ions,str):
        ions = [ions]
    atom = ions[0].split(' ')[0]
    for ion in ions:
        assert atom==ion.split(' ')[0]
    try:
        name = vars_to_name(atom,redshift,radfield)
        name,all_ions,temps,densities,Us,ion_fractions = data[name]
        name,all_ions,ret_ts,ret_rhos,PI_or_CI = cutoffs[name]
    except KeyError:
        print('no file %s found.'%name)
    if ax==None:
        fig,ax = plt.subplots(1,1)
    if show_cutoffs:
        if isinstance(show_cutoffs,str):
            line_colors = [show_cutoffs]
        else:
            line_colors = None
        plot_cutoffs(atom,redshift,radfield,ax=ax,ions=ions,log='values',data = data,\
                     cutoffs = cutoffs,colors = line_colors,linewidth = cutoff_linewidth)
    if levels == 'default':
        levels = [0.01,0.1]
    to_plot = get_contour_values(name,mechanism,ions,data=data,cutoffs=cutoffs)
    if colors is None:
        colors = cycle_colors
    elif isinstance(colors,str):
        colors = [colors]*len(all_ions)
    for i,ion in enumerate(ions):
        ion_index = all_ions.index(ion)
        unique_temps = np.unique(temps)
        unique_densities = np.unique(densities)
        if show_values:
            if len(ions) == 1:
                C_to_plot = np.log10(to_plot[0])
            elif len(ions) == 2:
                C_to_plot = np.log10(to_plot[0]/to_plot[1])
                if i==0: print(f'log {ions[0]}/{ions[1]}')
            too_small=10
            C_to_plot[np.isnan(C_to_plot)] = -too_small
            C_to_plot[np.isinf(C_to_plot)] = -too_small
            C_to_plot[C_to_plot<-too_small] = -too_small
            C_to_plot[C_to_plot>too_small] = too_small
            im = ax.pcolormesh(np.log10(unique_densities),np.log10(unique_temps),C_to_plot,cmap='seismic')
            if i==0:
                fig.colorbar(im, ax=ax)
        ax.contour(np.log10(unique_densities),np.log10(unique_temps),to_plot[i],\
                   levels = np.array(levels),colors=colors[ion_index],linestyles = 'solid')
    ax.set_title('Contours for %s at z=%.2f'%(ions,redshift))
    return ax

def frac_ratio(ion,redshift,radfield = 'HM12',name_to_compare_to = 'split',ax=None, show_cutoffs= False,\
               data=None,cutoffs=None,log=True):
    atom = ion.split(' ')[0]
    name1 = vars_to_name(atom,redshift,radfield)
    ion_fractions1 = get_contour_values(name1,'both',[ion],data=data,cutoffs=cutoffs)
    tocompare1 = ion_fractions1[0]
    _,_,temps1,densities1,_,_ = data[name1]
    x = np.log10(np.unique(densities1))
    y = np.log10(np.unique(temps1))
    ion_fractions2 = get_contour_values(name1,name_to_compare_to,[ion],data=data,cutoffs=cutoffs)
    tocompare2 = ion_fractions2[0]
    if ax is None:
        fig,ax = plt.subplots(1,1)
    else:
        fig = plt.gcf()

    threshold = 1e-6
    negligible_1 = tocompare1<threshold
    negligible_2 = tocompare2<threshold
    
    if log:
        C = np.log10(tocompare2/tocompare1)
    else:
        C = tocompare2-tocompare1
    
    nan_vals = np.logical_or(np.isnan(C),np.isinf(C))
    nan_replace = -3.0
    C[nan_vals] = nan_replace
    C[np.logical_and(negligible_1,negligible_2)] = nan_replace
    from matplotlib import colors
    if log:
        divnorm=colors.TwoSlopeNorm(vmin=-3.0, vcenter=0., vmax=3.0)
    else:
        divnorm=colors.TwoSlopeNorm(vmin=-0.5, vcenter=0., vmax=0.5)
    im = ax.pcolormesh(x,y,C,cmap='seismic', norm=divnorm)
    fig.colorbar(im, ax=ax)
    if show_cutoffs:
        if show_cutoffs == 'black':
            colors = ['black']
        else:
            colors = None
        plot_cutoffs(atom,redshift,ax=ax,ions=[ion],log='values',data = data,\
                     cutoffs = cutoffs,colors = colors,linewidth = 7)
    ax.set_title('%s ion ratio [%s vs %s]'%(ion,name1,name_to_compare_to))
    return ax