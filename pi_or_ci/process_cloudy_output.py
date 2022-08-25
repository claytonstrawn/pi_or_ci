import numpy as np
import os
import roman

def weird_inputs(strings):
    for s in strings:
        try:
            fs = float(s)
        except:
            return True
        if float(s)<0 or float(s)>1:
            return True
    return False 

def process_line(i,line,temp,dens,ioniz,count):
    atom_names = {'Li':'Lithium','Be':'Beryllium','B':'Boron',
                  'C':'Carbon','N':'Nitrogen','O':'Oxygen',
                  'F':'Fluorine','Ne':'Neon','Na':'Sodium',
                  'Mg':'Magnesium','Al':'Aluminium','Si':'Silicon'}
    atom_names_r = {'Lithium':'Li','Beryllium':'Be','Boron':'B',
                  'Carbon':'C','Nitrogen':'N','Oxygen':'O',
                  'Fluorine':'F','Neon':'Ne','Sodium':'Na',
                  'Magnesium':'Mg','Aluminium':'Al','Silicon':'Si'}
    atom_numbers = {'Li':3,'Be':4,'B':5,
                  'C':6,'N':7,'O':8,
                  'F':9,'Ne':10,'Na':11,
                  'Mg':12,'Al':13,'Si':14}
    atoms = ['Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si']
    c = line[0:22]
    c1 = line[0:13]
    c2 = line[0:15]
    if not c:
        print("End of file")
        for f_write in outputfiles:
            f_write.close()
        assert False
    elif c == '           U(1.0----):':
        ioniz=float(line[22:33])
        return None,None,temp,dens,ioniz,count
    elif c1 == ' ####  1  Te:':
        count=count+1
        temp=float(line[13:22])
        dens=float(line[28:37])
        return None,None,temp,dens,ioniz,count
    else:
        try:
            atom_fullname = c2.split()[0]
            atom = atom_names_r[atom_fullname]
        except:
            return None,None,temp,dens,ioniz,count
    atom_num = atom_numbers[atom]
    if atom in ['Si','Al']:
        try:
            offset = int(c2.split()[1])
        except:
            return None,None,temp,dens,ioniz,count
        strings = line.split()[2:atom_num+3]
    else:
        offset = 0
        strings = line.split()[1:atom_num+2]
    if weird_inputs(strings):
        return None,None,temp,dens,ioniz,count
    values = [0]*(atom_num+1)
    num_to_iterate = (atom_num+1) if (atom not in ['Si','Al']) else (12+1)
    if count%3 == 0:
        for j in range(num_to_iterate):
            v = float(strings[j])
            values[j+offset] = v
    return atom,values,temp,dens,ioniz,count

def cloudyoutput_to_ionfracs(inputfilename,radfield='HM12',redshift = 0.0,
                             inputpathstart = '/Users/claytonstrawn/research_personal/cloudy_scripts/',
                             outputpathstart = '/Users/claytonstrawn/CLOUDY_data/'):
    atom_names = {'Li':'Lithium','Be':'Berylliu','B':'Boron',
                  'C':'Carbon','N':'Nitrogen','O':'Oxygen',
                  'F':'Fluorine','Ne':'Neon','Na':'Sodium',
                  'Mg':'Magnesium','Al':'Aluminium','Si':'Silicon'}
    atom_names_r = {'Lithium':'Li','Berylliu':'Be','Boron':'B',
                  'Carbon':'C','Nitrogen':'N','Oxygen':'O',
                  'Fluorine':'F','Neon':'Ne','Sodium':'Na',
                  'Magnesium':'Mg','Aluminium':'Al','Silicon':'Si'}
    atom_numbers = {'Li':3,'Be':4,'B':5,
                  'C':6,'N':7,'O':8,
                  'F':9,'Ne':10,'Na':11,
                  'Mg':12,'Al':13,'Si':14}
    atoms = ['Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si']
    f_read_name = inputfilename
    f_write_names = []
    for atom in atoms:
        f_write_names.append(f'{atom}_z{redshift:.02f}_{radfield}.dat')
    outputfiles = []
    for j,atom in enumerate(atoms):
        file = open(outputpathstart+f_write_names[j],'w')
        str_to_write = '# T[k]  Density[cm^-3]  U(1)[] '
        for k in range(atom_numbers[atom]+1):
            if k == 0:
                roman_num = '0'
            else:
                roman_num = roman.toRoman(k)
            str_to_write+=f'{atom}{roman_num}/{atom} '
        str_to_write=str_to_write[:-1]+'\n'
        file.write(str_to_write)
        outputfiles.append(file)
    with open(inputpathstart+f_read_name,'rt') as f_read:
        length_of_file = 27829794 #this is a guess according to my usual parameters, but it does change
        i=0
        count=0
        temp = None
        dens = None
        ioniz = None
        for i,line in enumerate(f_read):
            if i%5000000 == 0:
                print(f'{i}/{length_of_file} ({i/(length_of_file/100):.2f}%)...')

            atom,values,temp,dens,ioniz,count = process_line(i,line,temp,dens,ioniz,count)
            if values is None:
                continue
            if count%3 == 0:
                if atom == 'Li':
                    atom2,values2,_,_,_,_ = process_line(i,'Beryllium'+line.split('Berylliu')[1],temp,dens,ioniz,count)
                    f_write2 = outputfiles[atoms.index(atom2)]
                    str_to_write2 = f"{temp: .03e} {dens: .03e} {ioniz: .03e}"
                    for v in values2:
                        str_to_write2+=f' {v: .05e}'
                    str_to_write2+='\n'
                    f_write2.write(str_to_write2)
                f_write = outputfiles[atoms.index(atom)]
                str_to_write = f"{temp: .03e} {dens: .03e} {ioniz: .03e}"
                for v in values:
                    str_to_write+=f' {v: .05e}'
                str_to_write+='\n'
                f_write.write(str_to_write)
        for file in outputfiles:
            file.close()
        print(f'{i}/{length_of_file} ({i/(length_of_file/100):.2f}%)...')

