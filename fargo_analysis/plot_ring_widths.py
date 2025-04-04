import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
from scipy.signal import peak_widths
from scipy.signal import find_peaks
from numpy.polynomial import Polynomial
import argparse



# =============== Customisable parameters ===============

sigma_width = 3
r_hill_width = 3
output = 600

# ================== Fitting functions ==================

def gaussian(x, a, x0, b): 
    return a*np.exp(-(x-x0)**2/(2*b**2)) 

    # y = A*np.exp(-1*B*x**2) 
    # return y 


def find_ring_peak(radii, Rp, profile):
    '''
    Written by Josh Adams
    '''
    bounded_cells = np.where((radii>Rp*rlowerbound) & (radii<Rp*rupperbound))
    r_bounded = radii[bounded_cells]
    profile_bounded = profile[bounded_cells]

    # Find the normalized surface density maximum value
    profilemax = np.max(profile_bounded)
    # Find the radius where that occurs
    i_profilemax = np.where(profile_bounded == profilemax)[0][0]
    r_peak = r_bounded[i_profilemax]
    # Find the non-normalized surface density value where *that* occurs
    profilemax = (profile_bounded)[np.where(r_bounded == r_peak)]
    i_profilemax += bounded_cells[0][0]   # get index in full array, not bounded array
    return r_peak, i_profilemax, profilemax  # convenient order for plotting


def find_ring_troughs(radii, dust_mass_normed, i_peak, bins=np.arange(-5,0)):
    grad_dmass = np.gradient(dust_mass_normed, radii)
        # gradients left of peak, in order of increasing distance from peak
    grad_left = grad_dmass[i_peak-1::-1]
        # gradients right of peak, in order of increasing distance from peak
    grad_right = grad_dmass[i_peak+1:]

    # check for change in gradient sign left of peak
    grad_li = 0
    while grad_left[grad_li]*grad_left[grad_li-1] >= 0:
        grad_li += 1
    left_trough_i = i_peak - grad_li

    # check for change in gradient sign right of peak
    grad_ri = 0
    while grad_right[grad_ri]*grad_right[grad_ri-1] >= 0:
        grad_ri += 1
    right_trough_i = i_peak + grad_ri

    return radii[left_trough_i], left_trough_i, dust_mass_normed[left_trough_i], radii[right_trough_i], right_trough_i, dust_mass_normed[right_trough_i]


# ============== Read in data from models ==============

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fit dust rings', prefix_chars='-')

    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/astro/phrkvg/simulations/dusty_fargo_models/"],help="working directory containing simulations")
    parser.add_argument('-sim', metavar='sim', type=str, nargs=1, default=["10Me"] ,help="simulation directory containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images" ,help="directory to save plots to")
    parser.add_argument('-o', metavar='outputs',default=[50], type=int, nargs="*" ,help="outputs to plot")
    parser.add_argument('-plots', metavar='plots',default=["dmass"], type=str, nargs="*" ,help="plots to produce")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-style', metavar='style', type=str, nargs="*", default=["publication"] ,help="style sheet to apply to plots")

    args = parser.parse_args()
    outputs = args.o
    plots = args.plots
    wd = args.wd[0]
    sim = args.sim[0]
    simdir = f"{wd}/{sim}/"
    plot_window = args.plot_window
    plots_savedir = args.savedir
    style = args.style

    if plot_window:
        matplotlib.use('TkAgg')

    params_file = f'{simdir}/variables.par'
    params_dict = {}

    # Load model params into dict
    param_lines = open(params_file).readlines()
    for line in param_lines:
        if line.split():
            param_label, param_value = line.split()[0:2]
            params_dict.update([(param_label, param_value)])

    nphi = int(params_dict['NX'])
    nrad = int(params_dict['NY'])
    hr0 = float(params_dict['ASPECTRATIO'])      # aspect ratio at R=1AU
    ndust = int(params_dict['NDUST'])
    # alpha = float(params_dict['ALPHA'])
    spacing = str(params_dict['SPACING'])
    max_stokes = float(params_dict['STOKES'])
    stokes = np.logspace(np.log10(max_stokes),np.log10(max_stokes*10**(-ndust+1)),ndust)

    dt_orbits = int(float(params_dict['DT'])/(2*np.pi))   # 2pi = 1 orbit = 1 yr
    ninterm = float(params_dict['NINTERM'])               # number of dts between outputs
    dt_outputs = dt_orbits*ninterm                        # time between outputs
    
    planet_data = np.unique(np.loadtxt(f"{simdir}/planet0.dat"), axis=0)
    xp, yp = planet_data[output][:,1], planet_data[output][:,2]
    rp = ((xp**2) + (yp**2))**0.5

    r_cells = np.loadtxt(f'{simdir}/domain_y.dat')[3:-3]             # ignore ghost cells
    phi_cells = np.loadtxt(f'{simdir}/domain_x.dat')
    if spacing == "Linear":
        radii = np.array([(r_cells[n]+r_cells[n+1])/2 for n in range(len(r_cells)-1)])
        delta_r = radii[1]-radii[0]
    else:     # Log grid
        radii = np.array([np.exp((np.log(r_cells[n])+np.log(r_cells[n+1]))/2) for n in range(len(r_cells)-1)])
        delta_log_r = np.log(radii[1]) - np.log(radii[0])
        delta_r = radii*delta_log_r
    phis = np.array([(phi_cells[n]+phi_cells[n+1])/2 for n in range(len(phi_cells)-1)])

    # Get gas and dust sigma
    # sigma_gas = np.zeros((len(outputs), nrad, nphi))
    sigma_dust = np.zeros((ndust, nrad, nphi))

    for n in np.arange(ndust):
        dust_file = f"dustdens{n}_{int(t)}.dat"
        sigma_dust[n] = np.fromfile(simdir+dust_file).reshape(nrad,nphi)/(1.125e-7)
    # sigma_gas_azimsum = np.sum(sigma_gas, axis=2)                  # sum over all phi   
    # sigma_gas_1D = sigma_gas_azimsum/nphi                         # dimensions: (noutputs, nrad) 

    sigma_dust_1D = np.sum(sigma_dust, axis=2)/nphi                         # dimensions: (noutputs, ndust, nrad)   
    

    # ================ Calculate ring width ================





    # =========== Plot ring width vs planet mass ===========


    print(f"-------------------\nPlotting outputs {outputs} for {sim}\n=============")

    if plot_window:
        plt.show()

