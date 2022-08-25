import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import roman
import matplotlib.pyplot as plt
cycle_colors = (list(plt.rcParams['axes.prop_cycle'].by_key()['color'])+ ["red", "indigo",'blue','gold','magenta'])*4
colors1 = mcolors.BASE_COLORS
colors2 = mcolors.TABLEAU_COLORS
colors3 = mcolors.CSS4_COLORS
colors4 = cycle_colors

all_color_names = []
for name, color in colors1.items():
    all_color_names+=[name]
for name, color in colors2.items():
    all_color_names+=[name]
for name, color in colors3.items():
    all_color_names+=[name]
for name in colors4:
    all_color_names+=[name]
all_marker_names = ['h','o','X','*','^','v','>','<','s','P','H','D','1','2','3','4']

all_linestyle_names = ['','-','--',':','-.']
def iscolor(arg):
    if isinstance(arg,str) and arg in all_color_names:
        return True
    elif isinstance(arg,tuple) and len(arg) == 3:
        for c in arg:
            if not isinstance(c,(int,float)):
                return False
        return True

def add_custom_legend(ax,labels_colors_dict,loc = None,bbox = None,keep_old_legend = True,ncol = 1):
    custom_lines = []
    labels_list = []
    for i,label in enumerate(labels_colors_dict):
        color,marker,linestyle = None,None,None
        labels_list.append(label)
        args = labels_colors_dict[label]
        if isinstance(args,str) or iscolor(args):
            args = [args]
        if args is None:
            args = [cycle_colors[i]]
        for arg in args:
            if iscolor(arg) and color is None:
                color = arg
            elif arg in all_linestyle_names and linestyle is None:
                linestyle = arg
            elif arg in all_marker_names and marker is None:
                marker = arg
            else:
                assert False, "Argument %s not recognized"%(arg,)
        color = 'k' if color is None else color
        linestyle = '' if (color == 'k' and linestyle is None) else linestyle
        custom_lines.append(Line2D([], [], marker = marker,color = color,linestyle = linestyle))
    if keep_old_legend:
        old_legend = ax.get_legend()
        ax.legend(custom_lines,labels_list,bbox_to_anchor = bbox,loc = loc,framealpha = 1.0,ncol=ncol)
        ax.add_artist(old_legend)
    else:
        old_legend = ax.get_legend()
        if old_legend is not None:
            old_legend.remove()
        ax.legend(custom_lines,labels_list,bbox_to_anchor = bbox,loc = loc,framealpha = 1.0,ncol=ncol)
        
def name_to_vars(name):
    if '/' in name:
        name = name.split('/')[-1]
    if name[-4:]=='.dat':
        name = name[:-4]
    splt = name.split('_')
    atom = splt[0]
    redshift = float(splt[1].replace('z',''))
    radfield = splt[2]
    return atom,redshift,radfield

def vars_to_name(atom,redshift,radfield,append = ''):
    return '%s_z%.2f_%s%s'%(atom,redshift,radfield,append)

        
def log(message,loud,end='\n'):
    if loud:
        print(message,end=end)
        
def get_ions(atom):
    number = {'Li':3,'Be':4,'B':5,
              'C':6,'N':7,'O':8,
              'F':9,'Ne':10,'Na':11,
              'Mg':12,'Al':13,'Si':14}[atom]
    ions = []
    for i in range(number+1):
        ions.append("%s %s"%(atom,roman.toRoman(i+1)))
    return ions