from pi_or_ci.find_pi_cutoffs import find_PI_cutoff
from pi_or_ci.file_loader import get_data_from_file,get_all_data_and_cutoffs
from pi_or_ci.utils import vars_to_name

def pi_or_ci(ion,temperature,density = 316,redshift = 0.0,background = 'HM12',cutoffs = None):
    atom = ion.split(' ')[0]
    name = vars_to_name(atom,redshift,radfield)
    if cutoffs is None:
        _,cutoffs = get_all_data_and_cutoffs(whichatoms = atom,whichzs = redshift,whichradfields = background)
    _,ions,ret_ts,ret_rhos,_ = cutoffs[name]
    ion_index = ions.index(ion)
    cutoff_ts,cutoff_rhos = ret_ts[ion_index],ret_rhos[ion_index]
    eq_temp = np.interp(density,cutoff_rhos,cutoff_ts)
    if temperature>=eq_temp:
        return 'CI'
    else:
        return 'PI'