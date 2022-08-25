from pi_or_ci.utils import log,get_ions,cycle_colors,name_to_vars,vars_to_name
import numpy as np


def plot_fracs_at_temp(atom,redshift,temp,radfield = 'HM12',ax=None,ions=None,show_mechanism = False,data=None,\
                       cutoffs=None,legend_params = None,bold_last = False,threshold = .05):
    name = vars_to_name(atom,redshift,radfield)
    try:
        name,all_ions,temps,densities,Us,ion_fractions = data[name]
        name,all_ions,ret_ts,ret_rhos,PI_or_CI = cutoffs[name]
    except KeyError:
        print('no file %s found.'%name)
    if legend_params is None:
        legend_params = {'loc':'upper right','bbox_to_anchor':(1.23,1)}
    if ax==None:
        fig,ax = plt.subplots(1,1)
    closest_T = temps[np.argmin(np.abs(temps-temp))]
    densities_at_temp = densities[temps == closest_T]
    fractions_at_temp = ion_fractions[temps == closest_T]
    if ions is None:
        ions = all_ions
    for i,ion in enumerate(all_ions):
        if ion not in ions:
            continue
        if show_mechanism:
            mechanism = PI_or_CI[ret_ts==closest_T,i]
            if mechanism == -1:
                linestyle = '--'
            elif mechanism == 0:
                linestyle = '-'
                plot_rho = ret_rhos[ret_ts==closest_T,i]
                ax.loglog([plot_rho,plot_rho],[1e-10,1e10],color=cycle_colors[i],linestyle=':')
                star_y = np.interp(plot_rho,densities_at_temp,fractions_at_temp[:,i])
                ax.loglog([plot_rho],[star_y],color=cycle_colors[i],marker = '*',markersize = 10)
            elif mechanism == 1:
                linestyle = '-.'
            else:
                linestyle = ""
        else:
            linestyle = '-'
        ax.loglog(densities_at_temp,fractions_at_temp[:,i],label = ion,linestyle = linestyle,color=cycle_colors[i])
        if bold_last:
            if np.any(fractions_at_temp[-bold_last:,i]<=0):
                continue
            log_last_10_average = np.average(np.log10(fractions_at_temp[-bold_last:,i]))
            relative_stable_top = np.zeros(bold_last)+10**log_last_10_average*(1+threshold)
            relative_stable_bot = np.zeros(bold_last)+10**log_last_10_average*(1-threshold)
            ax.loglog(densities_at_temp[-bold_last:],fractions_at_temp[-bold_last:,i],
                      linestyle = linestyle,color=cycle_colors[i],linewidth = 2,marker = '.')
            ax.loglog([densities_at_temp[-bold_last//2]],[10**log_last_10_average],'ko')
            for line in [relative_stable_top,relative_stable_bot]:
                ax.loglog(densities_at_temp[-bold_last:],line,
                      linestyle = '--',color=cycle_colors[i],linewidth = 1)
    ax.set_ylim(1e-5,1.5)
    ax.legend(**legend_params)
    atom,redshift,radfield = name_to_vars(name)
    ax.set_title('%s ion fractions at z=%.2f, T=10**%1.2f K'%(atom,redshift,np.log10(closest_T)))
    ax.set_xlabel('density (cm-3)')
    ax.set_ylabel('ion fraction')
    return ax