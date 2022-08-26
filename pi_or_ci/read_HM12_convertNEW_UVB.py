# Santi Roca-Fabrega 07/20/2022 - Instituto de Astronomia de la UNAM in Ensenada
# Program to read the Cloudy HM12 table and apply the correction proposed by Haislmaier+20
# The output is a new Cloudy table with HM12+correction
# The new table need to be moved to the "data" folder inside Cloudy (same folder as the original HM12 table "hm12_galaxy.ascii") with this same name "hm12_galaxy.ascii"
# REMEMBER to save the old "hm12_galaxy.ascii" before copying the new one!!!

import numpy as np
from matplotlib import pyplot as plt

##########################
### Initial parameters ###
##########################
# See details in the "functions section of this program"

Path_to_data="/Users/claytonstrawn/research_personal/"

#range used here: -0.28 , -0.835, -1.945, -2.5  
beta= 1.0 # default 1.0
import sys

alpha = -float(sys.argv[1]) # default -1.41
E0= 1.0 # defaul 1.0
m= -9.158 # default -9.158
b= -193.15 # default -193.15

favorite_redshift_index = 1 # default 8 (z~0.5)
if len(sys.argv)>2:
    print_plot_vals = bool(sys.argv[2])
else:
    print_plot_vals = False

##########################
# Read Cloudy HM12 tabke #
##########################

f = open(Path_to_data+'hm12_galaxy.ascii', 'r') # 'r' = read
lines = f.readlines()
i=0
for line in lines:
    if line.startswith("#"):
        i=1+i
    else:
        istart=i
        iredshifts=i+4
        ilambdas=i+5
        nredshifts=int(lines[iredshifts])
        nlambdas=int(lines[ilambdas])
        break
full_array=[]
for line in lines[istart+10:]:
    full_array.extend(line.split())
redshifts_list=full_array[0:nredshifts]
lambdas_list=full_array[nredshifts:nredshifts+nlambdas]


##########################
# Haislmaier+20 function #
##########################

# From eq. 5 in https://arxiv.org/pdf/2011.05335.pdf
# The m*E+b has been obtained from fitting a line in Figure 9 (upper panel) to the HM12 (solid black line) between log10(E[eV])=1.13 (=1Ryd), and log10(E[eV])=2.4 (=18.4Ryd)
# the result of this fit is m=-9.158 and b=-193.15
# we assume that E0 is 1Ryd (13.6eV, that is the ionization energy for the Hydrogen)

def f_E(E,func_flux,beta=1.0,alpha=-1.41,E0=1.0):
    if float(func_flux)==0.:
        corr=-999
    else:
        if E <= E0:
            corr= np.log10(float(func_flux))
        else:
            corr=np.log10(float(func_flux))+((alpha+1.41)*np.log10(float(E*13.606)))-(alpha+1.41)
    if  np.isnan(corr):
        return -999
    else:
        return corr
def corr_HM12(lamb,func_flux,beta=1.0,alpha=-1.41,E0=1.0,m=-9.158,b=-193.15):
    E=907.65/float(lamb)
    if f_E(E,func_flux,beta,alpha,E0) == -999:
        return 0.0
    else:
        if 1.0 < E < 18.4:
            correction_HM12=10**(beta*f_E(E,func_flux,beta,alpha,E0)+(1.0-beta)*(m*E+b))
        else:
            correction_HM12=10**(f_E(E,func_flux,beta,alpha,E0))
        return correction_HM12
for i in range(1,nredshifts):
    xs = [np.log10(907.65/float(lambdas_list[i])*13.606) for i in range(0,len(lambdas_list))]
    ys = [np.log10(float(full_array[nredshifts+i*nlambdas+k])) for k in range(0,len(lambdas_list))]
    if print_plot_vals:
        print('z =',redshifts_list[i])
        print('old hm12 xs:')
        print(xs)
        print('old hm12 ys:')
        print(ys)
    if i==1:
        color = 'k'
    else:
        color = None
    plt.plot(xs,ys,color = color,label=f'z={redshifts_list[i]}')
plt.xlim(0,3)
plt.ylim(-25.5,-18.5)
plt.legend()
plt.show()
##########################
###### Main Program ######
##########################

f1 = open(Path_to_data+'hm12_galaxy_new.ascii', 'w')
for line in lines[:istart+88]:
        f1.writelines(line)
for i in range(1,nredshifts+1):
    HM12_new=[corr_HM12(lambdas_list[k],full_array[nredshifts+i*nlambdas+k],beta,alpha,E0,m,b) for k in range(0,len(lambdas_list))]
#####
# plot result for z=0.5 (i=8)  and compare it with the original HM12
    if i == favorite_redshift_index:
        xs = [np.log10(907.65/float(lambdas_list[n])*13.606) for n in range(0,len(lambdas_list))]
        ys = [np.log10(corr_HM12(lambdas_list[k],full_array[nredshifts+favorite_redshift_index*nlambdas+k],beta,alpha,E0,m,b)) for k in range(0,len(lambdas_list))]
        if print_plot_vals:
            print('new hm12 xs:')
            print(xs)
            print('new hm12 ys:')
            print(ys)
        plt.plot(xs,ys,'k--',label='$\\alpha_{UV}$=%2.2f' % alpha + ', $\\beta_{HEII}$=%2.2f' % beta)
        plt.xlabel('log10(E[eV])')
        plt.ylabel('log10(J$\\nu$[ergs$s^{-1}$cm$^{-2}$Hz$^{-1}$sr$^{-1}$]')
        plt.title('z=0.0')
        plt.legend()
        plt.savefig('HM_modif_alpha%2.2f' % alpha + '_beta%2.2f.png' % beta)
        plt.show()
#####

#####
# Write the output file
    line='  '+str()+'  '+str()+''
    residual=len(lambdas_list)%10
    for l in range(0,int(len(lambdas_list)/10)):
        #line = '  '+'  '.join(HM12_new[10*l:10*(l+1)]+'\n'
        line='  '+str('%.3e' % HM12_new[10*l])+'  '+str('%.3e' % HM12_new[10*l+1])+'  '+\
                str('%.3e' % HM12_new[10*l+2])+'  '+str('%.3e' % HM12_new[10*l+3])+'  '+\
                str('%.3e' % HM12_new[10*l+4])+'  '+str('%.3e' % HM12_new[10*l+5])+'  '+\
                str('%.3e' % HM12_new[10*l+6])+'  '+str('%.3e' % HM12_new[10*l+7])+'  '+\
                str('%.3e' % HM12_new[10*l+8])+'  '+str('%.3e' % HM12_new[10*l+9])+'\n'
        f1.writelines(line)
    line='  '
    for l in range(0,residual-1):
        line=line+str('%.3e' % HM12_new[len(lambdas_list)-residual+l])+'  '
    line=line+str('%.3e' % HM12_new[len(lambdas_list)-1])+'\n'
    f1.writelines(line)
f1.close()
