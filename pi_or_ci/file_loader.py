import os
import roman
import numpy as np
from pi_or_ci.find_pi_cutoffs import find_PI_cutoff
from pi_or_ci.utils import name_to_vars,vars_to_name

def get_data_from_file(filename):
    f = open(filename)
    atom,redshift,radfield = name_to_vars(filename)
    lines = f.readlines()
    ions = lines[0].split()[4:]
    lines = lines[1:]
    numions = len(ions)
    for i in range(numions):
        ions[i] = '%s %s'%(atom,roman.toRoman(i+1))
    temps = np.zeros(len(lines))
    densities = np.zeros(len(lines))
    Us = np.zeros(len(lines))
    ion_fractions = np.zeros((len(lines),numions))
    for i,line in enumerate(lines):
        items = line.split()
        temps[i] = float(items[0])
        densities[i] = float(items[1])
        Us[i] = float(items[2])
        for j in range(numions):
            ion_fractions[i,j] = items[j+3]
    f.close()
    return vars_to_name(atom,redshift,radfield),ions,temps,densities,Us,ion_fractions 

def get_all_data_and_cutoffs(location = 'default',loud=False,\
                             whichatoms = 'all',whichzs = 'all',whichradfields = 'all',**kwargs):
    if location == 'default':
        _ROOT = os.path.abspath(os.path.dirname(__file__))
        location = os.path.join(_ROOT,'cloudy_data')
    names = []
    if whichatoms == 'all':
        whichatoms = ['Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si']
    elif isinstance(whichatoms,str):
        whichatoms = [whichatoms]
    if whichzs == 'all':
        whichzs = [0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0]
    elif whichzs == 'ints': 
        whichzs = [0.0,1.0,2.0,3.0,4.0]
    elif isinstance(whichzs,(int,float)):
        whichzs = [whichzs]
    if whichradfields == 'all':
        whichradfields = ['HM12','NONE','a0.28']
    elif isinstance(whichradfields,str):
        whichradfields = [whichradfields]
    for a in whichatoms:
        for z in whichzs:
            for r in whichradfields:
                names.append(vars_to_name(a,z,r,append = '.dat'))
    filenames = []
    for name in names:
        to_search_for = os.path.join(location,name)
        if os.path.exists(to_search_for):
            filenames.append(to_search_for)
        else:
            if loud:
                print('file %s not found'%to_search_for)
    data = {}
    cutoffs = {}
    for i,f in enumerate(sorted(filenames)):
        name = f.split('/')[-1].replace('.dat','')
        data[name] = get_data_from_file(f)
        cutoffs[name] = find_PI_cutoff(data[name],loud=loud,**kwargs)
    return data,cutoffs



def make_list_of_names(atom,redshifts,radfields):
    if isinstance(redshifts,(int,float)):
        redshifts = [redshifts]
    if isinstance(radfields,(str)):
        radfields = [radfields]
    names = []
    for i,redshift in enumerate(redshifts):
        for j,radfield in enumerate(radfields):
            names.append(vars_to_name(atom,redshift,radfield))

    return names