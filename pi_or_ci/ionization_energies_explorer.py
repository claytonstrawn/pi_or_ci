import numpy as np
import roman
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import linregress
from pi_or_ci.utils import add_custom_legend,cycle_colors
from pi_or_ci.utils import all_marker_names as symbols
from pi_or_ci.cutoff_plots import cutoffs_for_ion_at_redshift

#source = https://en.wikipedia.org/wiki/Molar_ionization_energies_of_the_elements
first10 = """
3	Li	lithium	520.2	7298.1	11,815.0
4	Be	beryllium	899.5	1757.1	14,848.7	21,006.6
5	B	boron	800.6	2427.1	3659.7	25,025.8	32,826.7
6	C	carbon	1086.5	2352.6	4620.5	6222.7	37,831	47,277.0
7	N	nitrogen	1402.3	2856	4578.1	7475.0	9444.9	53,266.6	64,360
8	O	oxygen	1313.9	3388.3	5300.5	7469.2	10,989.5	13,326.5	71,330	84,078.0
9	F	fluorine	1681.0	3374.2	6050.4	8407.7	11,022.7	15,164.1	17,868	92,038.1	106,434.3
10	Ne	neon	2080.7	3952.3	6122	9371	12,177	15,238.90	19,999.0	23,069.5	115,379.5	131,432
11	Na	sodium	495.8	4562	6910.3	9543	13,354	16,613	20,117	25,496	28,932	141,362
12	Mg	magnesium	737.7	1450.7	7732.7	10,542.5	13,630	18,020	21,711	25,661	31,653	35,458
13	Al	aluminium	577.5	1816.7	2744.8	11,577	14,842	18,379	23,326	27,465	31,853	38,473
14	Si	silicon	786.5	1577.1	3231.6	4355.5	16,091	19,805	23,780	29,287	33,878	38,726
"""
second10 = """
11	Na	sodium	159,076
12	Mg	magnesium	169,988	189,368
13	Al	aluminium	42,647	201,266	222,316
14	Si	silicon	45,962	50,502	235,196	257,923
"""

def get_energy_dict():
    convert_kJpermol_to_eV = 0.0103642688 #https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
    atoms = []
    energy_dict = {}
    energy_dict_cumulative = {}
    number_dict = {}
    for line in first10.split('\n'):
        if line in ['']:
            continue
        s = line.split('\t')
        number = s[0]
        atom = s[1]
        energies = s[3:]
        energies = [float(x.replace(',',''))*convert_kJpermol_to_eV for x in energies]
        energy_dict[atom] = energies
        number_dict[atom] = int(number)
        atoms+=[atom]
    for line in second10.split('\n'):
        if line == '':
            continue
        s = line.split('\t')
        atom = s[1]
        energies = s[3:]
        energies = [float(x.replace(',',''))*convert_kJpermol_to_eV for x in energies]
        energy_dict[atom] = energy_dict[atom]+energies
    for atom in atoms:
        energy_dict[atom] = np.array(energy_dict[atom])
        energy_dict_cumulative[atom] = np.cumsum(energy_dict[atom])
    ions = []
    for atom in atoms:
        max_ionization = number_dict[atom]+1
        for v in range(1,max_ionization+1):
            ions += [atom+' '+roman.toRoman(v)]
    ionization_energy_dict = {}
    for ion in ions:
        atom,ionization = ion.split(' ')
        ionization = roman.fromRoman(ionization)
        if ionization == 1:
            ionization_energy_dict[ion] = np.nan
        else:
            ionization_energy_dict[ion] = energy_dict_cumulative[atom][ionization-2]
    return atoms,ions,ionization_energy_dict

def get_cutoffs_dict(ions,redshift,t_used,cutoffs = None):
    if t_used == 'median':
        f = np.median
    elif t_used == 'min':
        f = np.amin
    elif t_used == 'mean':
        f = np.average
    else:
        assert False,'t_used = %s is not understood'%t_used
    PI_cutoff_temps_dict = {}
    for ion in ions:
        try:
            PI_cutoff_temps_dict[ion] = f(cutoffs_for_ion_at_redshift(ion,redshift,cutoffs=cutoffs)[1])
        except:
            PI_cutoff_temps_dict[ion] = np.nan
    return PI_cutoff_temps_dict

def get_energy_cutoffs_dicts(redshift=2.0,t_used = 'min',cutoffs = None):
    atoms,ions,ionization_energy_dict = get_energy_dict()
    return atoms,ions,ionization_energy_dict,get_cutoffs_dict(ions,redshift,t_used,cutoffs = cutoffs)

def view_energies(atoms,ions,ionization_energy_dict=None,cutoffs_dict = None,\
                  funky_symbols = False,ax = None,legend = True,grid = False):
    if ax is None:
        _,ax = plt.subplots(1)
    used_atoms = []
    roman_symbols = []
    for ion in ions:
        atom,ionization_r = ion.split(' ')
        i = atoms.index(atom)
        roman_symbols+=[] if ionization_r in roman_symbols else [ionization_r]
        used_atoms+=[] if atom in used_atoms else [atom]
        ionization = roman.fromRoman(ionization_r)
        if ionization == 1:
            continue
        color = cycle_colors[ionization-1]
        if funky_symbols:
            symbol = symbols[i]
        else:
            symbol = 'o'
        if ionization == 1:
            x = i-first_offset
        ax.semilogy(i,ionization_energy_dict[ion],symbol,color = color,markeredgecolor='black',markersize = 10)
    labels_colors_dict1 = {}
    for j,atom in enumerate(used_atoms):
        labels_colors_dict1[atom] = symbols[j]
    labels_colors_dict2 = {}
    for j,ionization_r in enumerate(roman_symbols):
        labels_colors_dict2[ionization_r] = cycle_colors[j]
    if legend:
        add_custom_legend(ax,labels_colors_dict1,keep_old_legend = False,bbox = (1.18,1.02),loc = 'upper right')
        if funky_symbols:
            add_custom_legend(ax,labels_colors_dict2,keep_old_legend = True,bbox = (1.35,1.02),loc = 'upper right')
    if grid:
        ax.grid(which = 'both')
    ax.set_xticks(range(len(atoms)))
    ax.set_xticklabels(atoms)
    ax.tick_params(axis='both',labelsize=15)
    ax.set_ylabel('ionization energy (eV)',size = 15)
    
def view_cutoffs(atoms,ions,ionization_energy_dict,PI_cutoff_temps_dict,units = 'K',legend = True,\
                 funky_symbols = False,ax = None,fig=None,grid = False,first_offset = 0.1):
    if ax is None:
        fig,ax = plt.subplots(1)
    used_atoms = []
    roman_symbols = []
    for ion in ions:
        atom,ionization = ion.split(' ')
        i = atoms.index(atom)
        roman_symbols+=[] if ionization in roman_symbols else [ionization]
        used_atoms+=[] if atom in used_atoms else [atom]
        ionization = roman.fromRoman(ionization)
        boltzmann = 8.617333262e-5 #eV/K
        if units == 'K':
            convert = 1
        elif units == 'eV':
            convert = 1.5*boltzmann
        elif units == 'relative':
            convert = 1.5*boltzmann/ionization_energy_dict[ion]
            convert2 = 1.5*boltzmann
        color = cycle_colors[ionization-1]
        if funky_symbols:
            symbol = symbols[i]
        else:
            symbol = 'o'
        overlaps_left = ['Ne VIII','Al XI']
        overlaps_right = ['Ne IX','Al XII']
        if ionization == 1 or ion in overlaps_left:
            x = i-first_offset
        elif ionization == 2 or ion in overlaps_right:
            x = i+first_offset
        else:
            x = i
        ax.semilogy(x,PI_cutoff_temps_dict[ion]*convert,symbol,color = color,markeredgecolor='black',markersize = 10)
    labels_colors_dict1 = {}
    for j,atom in enumerate(used_atoms):
        labels_colors_dict1[atom] = symbols[j]
    labels_colors_dict2 = {}
    for j,ionization_r in enumerate(roman_symbols):
        labels_colors_dict2[ionization_r] = cycle_colors[j]
    if legend:
        add_custom_legend(ax,labels_colors_dict1,keep_old_legend = False,bbox = (1.18,1.02),loc = 'upper right')
        if funky_symbols:
            add_custom_legend(ax,labels_colors_dict2,keep_old_legend = True,bbox = (1.35,1.02),loc = 'upper right')
    elif not legend and units == 'K':
        ylims = (6e3,7e6)
        ylims2 = (ylims[0]*1.5*boltzmann,ylims[1]*1.5*boltzmann)
        ax2 = ax.twinx()
        ax2.semilogy([],[],'o',color = cycle_colors[i])
        ax2.set_ylim(ylims2[0],ylims2[1])
        ax2.set_ylabel('PI temp threshold (eV)',size = 15)
        ax2.tick_params(axis='both',labelsize=15)
        ax.set_ylim(ylims[0],ylims[1])
    if grid:
        ax.grid(which = 'both')
    ax.set_xticks(range(len(atoms)))
    ax.set_xticklabels(atoms)
    ax.set_ylabel('PI temp threshold (%s)'%units,size = 15)
    ax.tick_params(axis='both',labelsize=15)
    return ax
        
def plot_cutoffs_by_energies(atoms,ions,ionization_energy_dict,PI_cutoff_temps_dict,legend = True,\
                             units = 'K',log = 'default',funky_symbols = True,ax = None,\
                             xlims_for_linreg = (2e1,8e3),loud = False,linecolor = 'black',\
                            lineedgecolor = None):
    boltzmann = 8.617333262e-5 #eV/K
    used_atoms = []
    roman_symbols = []
    xs,ys = [],[]
    if ax is None:
        _,ax = plt.subplots(1)
    if log == 'default':
        if units in ['K','eV']:
            func = ax.loglog
        elif units in ['relative']:
            func = ax.semilogx
    for ion in ions:
        if np.isnan(PI_cutoff_temps_dict[ion]):
            continue
        if units == 'K':
            convert = 1
        elif units == 'eV':
            convert = 1.5*boltzmann
        elif units == 'relative':
            convert = 1.5*boltzmann/ionization_energy_dict[ion]
        xs+=[ionization_energy_dict[ion]]
        ys+=[PI_cutoff_temps_dict[ion]*convert]
        atom,ionization_r = ion.split(' ')
        used_atoms+=[] if atom in used_atoms else [atom]
        roman_symbols+=[]if ionization_r in roman_symbols else [ionization_r]
        i = used_atoms.index(atom)
        ionization = roman.fromRoman(ionization_r)
        color = cycle_colors[ionization-1]
        if funky_symbols:
            symbol = symbols[i]
        else:
            symbol = 'o'
        if ionization == 1:
            continue
        func(ionization_energy_dict[ion],PI_cutoff_temps_dict[ion]*convert,
                   symbol,color = color,markersize = 10,markeredgecolor='black')
        if loud:
            print(f'{ion}, {ionization_energy_dict[ion]}, {PI_cutoff_temps_dict[ion]*convert}')
    xlims = (4e0,1e4)
    ax.set_xlim(xlims[0],xlims[1])
    if units == 'K':
        ylims = (6e3,7e6)
        textloc = (5e0,2e6)
        ylabel = 'PI cutoff temperature (K)'
    elif units == 'eV':
        ylims = (1e-1,7e2)
        textloc = (1.5e1,1e2)
        ylabel = 'PI cutoff temperature (eV)'
    elif units == 'relative':
        ylims = (-.01,.335)
        textloc = None
        ylabel = 'cutoff KE/ionization energy'
        ax.set_yticks([0,.05,.1,.15,.2,.25,.3])
    ax.set_ylim(ylims[0],ylims[1])
    ax.grid(which='both')
    ax.set_xlabel('ionization energy (eV)',size = 15)
    ax.set_ylabel(ylabel,size = 15)
    labels_colors_dict1 = {}
    for j,atom in enumerate(used_atoms):
        labels_colors_dict1[atom] = symbols[j]
    labels_colors_dict2 = {}
    for j,ionization_r in enumerate(roman_symbols):
        labels_colors_dict2[ionization_r] = cycle_colors[j]
    if legend:
        add_custom_legend(ax,labels_colors_dict1,keep_old_legend = False,bbox = (1.18,1.02),loc = 'upper right')
        if funky_symbols:
            add_custom_legend(ax,labels_colors_dict2,keep_old_legend = True,bbox = (1.35,1.02),loc = 'upper right')
    elif not legend and units == 'K':
        ylims2 = (ylims[0]*1.5*boltzmann,ylims[1]*1.5*boltzmann)
        ax2 = ax.twinx()
        ax2.semilogy([],[],'o',color = cycle_colors[i])
        ax2.set_ylim(ylims2[0],ylims2[1])
        ax2.set_ylabel('PI temp threshold (eV)',size = 15)
        ax2.tick_params(axis='both',labelsize=15)
        ax.set_ylim(ylims[0],ylims[1])
    ax.tick_params(axis='both',labelsize=15)
    if units in ['K','eV']:
        xs = np.array(xs)
        ys = np.array(ys)
        logxs = np.log10(np.array(xs[np.logical_and(xs>xlims_for_linreg[0],xs<xlims_for_linreg[1])]))
        logys = np.log10(np.array(ys[np.logical_and(xs>xlims_for_linreg[0],xs<xlims_for_linreg[1])]))
        line = linregress(logxs, logys)
        slope,intercept,_,_,_ = line
        sortxs = np.array(sorted(logxs))
        linexs = 10**sortxs
        lineys = 10**(sortxs*slope+intercept)
        t = ax.text(textloc[0],textloc[1],'y = %.2f x ^ %.3f'%(10**intercept,slope),size= 20)
        t.set_bbox(dict(facecolor='white'))
        import matplotlib.patheffects as pe
        func(linexs,lineys,color=linecolor, lw=2, path_effects=[pe.Stroke(linewidth=5, foreground=lineedgecolor), pe.Normal()])
    return ax

