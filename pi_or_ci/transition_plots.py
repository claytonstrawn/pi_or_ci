from pi_or_ci.file_loader import make_list_of_names
from pi_or_ci.utils import vars_to_name
import numpy as np

def plot_transitions(atom,redshift=0.0,radfield = 'HM12',ax=None,data=None,cutoffs=None,cmap = 'rainbow',maxval = 4):
    try:
        name = vars_to_name(atom,redshift,radfield)
        name,ions,ret_ts,ret_rhos,PI_or_CI = cutoffs[name]
    except KeyError:
        assert False,'no file %s found.'%name
    if ax==None:
        fig,ax = plt.subplots(1,1)
    xs = np.arange(len(ions))
    to_plot = np.array(PI_or_CI)
    if not maxval:
        maxval = np.amax(to_plot)
    minval = np.amin(to_plot)
    scaled_for_imshow = (to_plot-minval)/(maxval-minval) 
    try:
        ax.pcolormesh(xs,np.log10(np.array(ret_ts)),scaled_for_imshow,cmap = cmap,vmin = 0,vmax = 1)
    except ValueError as e:
        print(xs,np.log10(np.array(ret_ts)),np.array(PI_or_CI))
        raise e
              
    ax.set_xticks(xs)
    xticklabels = []
    for i,ion in enumerate(ions):
        #56 is sort of arbitrary, just messed with numbers until 
        #the right ones got newlines
        maxlength = 56/len(ions)
        if len(ion)>=maxlength:
            xticklabels.append(ion.replace(' ','\n'))
        else:
            xticklabels.append(ion)
    ax.set_xlim(xs[0]-0.5,xs[-1]+0.5)
    ax.set_xticklabels(xticklabels)
    for x in xs:
        ax.plot([x+0.5,x+0.5],[2,8],color = 'white')
    ax.set_title('%s ionization mechanisms at z=%.2f'%(atom,redshift))
    ax.set_ylabel('Temperature (K)')
    return scaled_for_imshow

def plot_transition_differences(atom,data,cutoffs,comparisons,scaling_of_std = 'default',ax=None):
    #comparisons is a tuple of dictionaries:
    #({'redshifts':[0,1,2,3,4],'radfields':'HM12'},{'redshifts':[0],'radfields':['HM12','a0.28','a0.835','a1.945','a2.5']})
    n = len(comparisons)+1
    first = None
    comparison_names = []
    for comparison in comparisons:
        redshifts = comparison['redshifts']
        radfields = comparison['radfields']
        names = make_list_of_names(atom,redshifts,radfields)
        if first is None:
            first = names[0]
        else:
            assert names[0] == first
        comparison_names.append(names)
    try:
        name,all_ions,ret_ts,ret_rhos,PI_or_CI = cutoffs[first]
    except KeyError:
        print('no file %s found.'%name)
    to_plot_final = np.repeat(PI_or_CI,n,axis=1)*0
    for i,line in enumerate(to_plot_final):
        for j,cell in enumerate(line):
            to_plot_final[i,j] = PI_or_CI[i,j//n]
    xs = np.arange(1,len(all_ions)+1)
    ys = np.arange(2,8.1)
    data_range = 10
    data_offset = ys[0]
    xs_repeats = np.repeat(xs,len(comparisons))
    if scaling_of_std == 'default':
        scaling_default_dict = {1:(3,),2:(4,3)}
        scaling_of_std = scaling_default_dict[len(comparisons)]
    assert len(scaling_of_std) == len(comparisons)
    # this number 3 scales the standard deviation numbers to make it stand out more on the colormap
    # has to be >= 3 to not overlap with negligible CI (=2)
    for k,comparison in enumerate(comparisons):
        names = comparison_names[k]
        start_arrays = []
        for name in names:
            try:
                name,all_ions,ret_ts,ret_rhos,PI_or_CI = cutoffs[name]
            except KeyError:
                print('no file %s found.'%name)
            start_arrays+=[PI_or_CI]
        to_plot_diffs = (np.std(np.array(start_arrays),axis=0)!=0).astype(int)*scaling_of_std[k]
        for i,line in enumerate(to_plot_final):
            for j,cell in enumerate(line):
                if j%n == (k+1) and to_plot_diffs[i,j//n]!=0:
                    to_plot_final[i,j] = to_plot_diffs[i,j//n]
    scaled_for_imshow = (to_plot_final-np.amin(to_plot_final))/(np.amax(to_plot_final)-np.amin(to_plot_final)) 
    ax.imshow(scaled_for_imshow,aspect = 'auto',interpolation = 'none',cmap = 'rainbow')
    ax.set_xticks(xs*n-2.0)
    xticklabels = []
    for i,ion in enumerate(all_ions):
        #56 is sort of arbitrary, just messed with numbers until 
        #the right ones got newlines
        maxlength = 56/len(all_ions)
        if len(ion)>=maxlength:
            xticklabels.append(ion.replace(' ','\n'))
        else:
            xticklabels.append(ion)
    ax.set_xticklabels(xticklabels)
    ax.set_xlim(xs[0]-1.5,xs[-1]*n-.5)
    ax.set_yticks((ys-data_offset)*data_range)
    ax.set_yticklabels(ys)
    ax.set_ylim((ys[0]-data_offset)*data_range,(ys[-1]-data_offset)*data_range)
    for x in (xs*n-1.5+1.0)[:-1]:
        ax.plot([x,x],[(ys[0]-data_offset)*data_range,(ys[-1]-data_offset)*data_range],color = 'white')
    for y in ys[1:]:
        ax.plot([xs[0]-1.5,xs[-1]*n-.5],[(y-data_offset)*data_range,(y-data_offset)*data_range],color = 'white')
    return scaled_for_imshow 