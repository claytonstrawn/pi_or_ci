from pi_or_ci.utils import vars_to_name,cycle_colors
import numpy as np

def get_cutoffs(ion,redshift,radfield = 'HM12',data=None,cutoffs=None,small_offset = .05):
    atom = ion.split(' ')[0]
    name = vars_to_name(atom,redshift,radfield)
    try:
        name,all_ions,ret_ts,ret_rhos,PI_or_CI = cutoffs[name]
    except KeyError:
        print('no file %s found.'%name)
    far_away_lims = [1e15,1e-15]
    i = all_ions.index(ion)
    if np.all(np.isnan(ret_rhos[:,i])):
        plot_rhos,plot_ts = [],[]
    else:
        current_rhos = ret_rhos[:,i][~np.isnan(ret_rhos[:,i])]
        current_ts = ret_ts[~np.isnan(ret_rhos[:,i])]
        min_t = min(current_ts)
        max_t = max(current_ts)
        plot_ts = np.concatenate([np.array([min_t]),current_ts,np.array([max_t])])
        plot_rhos = np.concatenate([np.array([far_away_lims[0]]),current_rhos,np.array([far_away_lims[1]])])
    return plot_rhos,plot_ts

def cutoffs_for_ion_at_redshift(ion,redshift,radfield = 'HM12',cutoffs = None):
    if cutoffs is None:
        assert False, "loading pregenerated cutoffs is not yet supported, "+\
                    "please create some with 'get_all_data_and_cutoffs' and pass as 'cutoffs'"
    atom = ion.split(' ')[0]
    name = vars_to_name(atom,redshift,radfield)
    plot_rhos,plot_ts = get_cutoffs(ion,redshift,radfield,cutoffs=cutoffs)
    return plot_rhos,plot_ts
    

def format_cutoffs(ions,redshifts,data,cutoffs,write_to = None,overwrite = False):
    to_return = ''
    for ion in ions:
        to_return+=ion+'\n'
        for redshift in redshifts:
            try:
                rho,t = get_cutoffs(ion,redshift,data=data,cutoffs=cutoffs)
            except:
                continue
            rho_str_to_print = '%1.2f, rho: ['%redshift
            rho_str_to_print=rho_str_to_print+' 'if len(rho) == 0 else rho_str_to_print
            for v in rho:
                rho_str_to_print+='%e'%v
                rho_str_to_print+=' '
            rho_str_to_print = rho_str_to_print[:-1]+']\n'

            t_str_to_print = '%1.2f,   t: ['%redshift
            t_str_to_print=t_str_to_print+' 'if len(rho) == 0 else t_str_to_print
            for v in t:
                t_str_to_print+='%.0f.'%v
                t_str_to_print+=' '
            t_str_to_print = t_str_to_print[:-1]+']\n'

            to_return+=rho_str_to_print
            to_return+=t_str_to_print
    if write_to is None:
        return to_return
    else:
        fmt = 'w' if overwrite else 'a'
        with open(write_to,fmt) as f:
            f.write(to_return)
        message = 'wrote table to %s.'%write_to
        if overwrite:
            message += ' (overwritten)'
        return message

def plot_cutoffs(atom,redshift,radfield='HM12',ax=None,ions=None,log='ticks',voffset = 0,\
                 marker = '.',linestyle = '-',data=None,cutoffs=None,\
                 markersize = 3,markercolor = 'k',colors = None,first_offset = 0.0,**kwargs):
    try:
        name = vars_to_name(atom,redshift,radfield)
        _,all_ions,_,_,_,_ = data[name]
    except KeyError:
        print('no file %s found.'%name)
    if ions is None:
        ions = all_ions
    if colors is None:
        colors = cycle_colors
    elif isinstance(colors,list):
        repeat = len(all_ions)//len(colors)+1
        colors = colors*repeat
    if ax is None:
        fig,ax = plt.subplots(1,1)
    for i,ion in enumerate(all_ions):
        if ion not in ions:
            continue
        density_lims = np.array([1e-8,3e2])
        temp_lims = np.array([1e2,1e8])
        plot_rhos,plot_ts = get_cutoffs(ion,redshift,radfield=radfield,data=data,cutoffs=cutoffs)
        if log == 'ticks':
            func = ax.loglog
        elif log == 'values':
            func = ax.plot
            plot_ts = np.log10(plot_ts)+voffset
            plot_rhos = np.log10(plot_rhos)
            density_lims = np.log10(density_lims)
            temp_lims = np.log10(temp_lims)
        if ion.split()[1] == 'I':
            if log == 'ticks':
                plot_ts=np.array(plot_ts)*(1-first_offset)
            elif log == 'values':
                plot_ts=np.array(plot_ts)-first_offset
        if linestyle != '':
            func(plot_rhos,plot_ts,linestyle,label=ion,color = colors[i],**kwargs)
        color = markercolor if markercolor is not None else colors[i]
        func(plot_rhos,plot_ts,marker,color=color,markersize=markersize,**kwargs)
    ax.legend(loc = 'lower left')
    ax.set_xlim(density_lims[0],density_lims[1])
    ax.set_ylim(temp_lims[0],temp_lims[1])
    ax.set_title('%s PI-CI cutoffs at z=%.2f'%(atom,redshift))
    ax.set_xlabel('density (cm-3)')
    ax.set_ylabel('Temperature (K)')
    return ax

def connect_cutoffs(atom,redshifts,radfield='HM12',ax=None,voffset = .2,\
                 data=None,cutoffs=None,max_num_ions = 3):
    try:
        name = vars_to_name(atom,redshifts[0],radfield)
        _,all_ions,_,_,_,_ = data[name]
    except KeyError:
        print('no file %s found.'%name)
    skip = len(all_ions)//max_num_ions
    disp_ions = all_ions[::skip]
    for i,ion in enumerate(disp_ions):
        ts = get_cutoffs(ion,redshifts[0],data=data,cutoffs=cutoffs)[1]
        ts = np.log10(ts)
        list_of_lines_xs = {}
        list_of_lines_ys = {}
        which_redshifts = {}
        for j,redshift in enumerate(redshifts):
            plot_rhos,plot_ts = get_cutoffs(ion,redshift,data=data,cutoffs=cutoffs)
            plot_rhos = np.log10(plot_rhos)
            plot_ts = np.log10(plot_ts)
            plot_ts[0]+=voffset/100
            plot_ts[-1]-=voffset/100
            for k,t in enumerate(ts):
                if j==0:
                    list_of_lines_xs[k] = []
                    list_of_lines_ys[k] = []
                    which_redshifts[k] = []
                try:
                    index_of_t = np.where(plot_ts==t)[0][0]
                except IndexError:
                    continue
                which_redshifts[k]+=[j]
                rho = plot_rhos[index_of_t]
                list_of_lines_xs[k]+=[rho]
                list_of_lines_ys[k]+=[t+voffset*j]
        for k,t in enumerate(ts):
            my_color = cycle_colors[all_ions.index(disp_ions[i])]
            ax.plot(list_of_lines_xs[k],list_of_lines_ys[k],color = my_color,linestyle = '-')
            for l in range(len(list_of_lines_xs[k])):
                symbols = ['o','x','*','^','v','>','<','s','1','2','3','4','P','+','d','D']
                marker = symbols[which_redshifts[k][l]]
                ax.plot([list_of_lines_xs[k][l]],[list_of_lines_ys[k][l]],marker,color = my_color)
    return ax