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

# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

# ================== Fitting functions ==================

def gaussian(x, a, x0, b): 
    return a*np.exp(-0.5*((x-x0)/b)**2) 
    # return a*np.exp(-0.5*((x-x0)/(b+c(x-x0)))**2)   # skewed gaussian

def find_ring_peak(radii, sigma):
    peak_i = np.argmax(sigma)
    return peak_i


def find_ring_troughs(radii, dust_mass_normed, i_peak, bins=np.arange(-5,0)):
    # TODO This is not bounded correctly!!! Fix it!!
    grad_dmass = np.gradient(dust_mass_normed, radii)
        # gradients left of peak, in order of increasing distance from peak
    grad_left = grad_dmass[i_peak-1::-1]
    print(grad_left)
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

    return left_trough_i, right_trough_i

# def find_ring_edges(radii, sigma, i_peak):
#     # Find ring edges based on 10% drop from peak sigma
#     # Allows for asymmetric/non-Gaussian shape
    
#     sigma_peak = sigma[i_peak]
#     # sigmas left of peak, in order of increasing distance from peak
#     sigma_left = sigma[i_peak-1::-1]
#     # sigmas right of peak, in order of increasing distance from peak
#     sigma_right = sigma[i_peak+1:]

#     # Identify where sigma first drops below 0.1*sigma_peak
#     left_edge_i = 0
#     while sigma_left[left_edge_i] >= 0.5*sigma_peak:
#         left_edge_i += 1
#     left_edge_i = i_peak - left_edge_i

#     right_edge_i = 0
#     while sigma_right[right_edge_i] >= 0.5*sigma_peak:
#         right_edge_i += 1
#     right_edge_i += i_peak
    
#     return left_edge_i, right_edge_i


# ============== Read in data from models ==============

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fit dust rings', prefix_chars='-')

    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/astro/phrkvg/simulations/dusty_fargo_models/"],help="working directory containing simulations")
    parser.add_argument('-sims', metavar='sim', type=str, nargs="*", default=["10Me"] ,help="simulation directory containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images" ,help="directory to save plots to")
    parser.add_argument('-o', metavar='outputs',default=[600], type=int, nargs="*" ,help="outputs to plot")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-style', metavar='style', type=str, nargs="*", default=["publication"] ,help="style sheet to apply to plots")

    args = parser.parse_args()
    output = args.o[0]
    wd = args.wd[0]
    sims = args.sims
    plot_window = args.plot_window
    plots_savedir = args.savedir
    style = args.style

    if plot_window:
        matplotlib.use('TkAgg')

    fig, ax = plt.subplots(figsize=(8,5))
    planet_masses = np.zeros((len(sims)))
    ring_widths = np.zeros((len(sims),3))

    # Iterate through models and calculate ring width for each St
    for s, sim in enumerate(sims):
        simdir = f"{wd}/{sim}/"
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
        alpha = float(params_dict['ALPHA'])
        spacing = str(params_dict['SPACING'])
        max_stokes = float(params_dict['STOKES'])
        stokes = np.logspace(np.log10(max_stokes),np.log10(max_stokes*10**(-ndust+1)),ndust)

        # Calculate timing parameters
        dt_orbits = int(float(params_dict['DT'])/(2*np.pi))   # 2pi = 1 orbit = 1 yr
        ninterm = float(params_dict['NINTERM'])               # number of dts between outputs
        dt_outputs = dt_orbits*ninterm                        # time between outputs
        
        # Calculate planet location and Hill radius
        planet_data = np.unique(np.loadtxt(f"{simdir}/planet0.dat"), axis=0)
        xp, yp = planet_data[output,1], planet_data[output,2]
        rp = ((xp**2) + (yp**2))**0.5
        mp = planet_data[output,7]
        mp = round(mp/(3.0027e-6),0)     # convert to earth masses
        planet_masses[s] = mp
        r_hill = rp*(mp/3)**(1/3)

        # Calculate radial values (cell centres)
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
        sigma_dust0 = np.zeros((ndust, nrad, nphi))

        for n in np.arange(ndust):
            dust_file = f"dustdens{n}_{output}.dat"
            sigma_dust[n] = np.fromfile(simdir+dust_file).reshape(nrad,nphi)/(1.125e-7)
            dust_file0 = f"dustdens{n}_0.dat"
            sigma_dust0[n] = np.fromfile(simdir+dust_file0).reshape(nrad,nphi)/(1.125e-7)
        # sigma_gas_azimsum = np.sum(sigma_gas, axis=2)                  # sum over all phi   
        # sigma_gas_1D = sigma_gas_azimsum/nphi                         # dimensions: (noutputs, nrad) 

        sigma_dust_1D = np.sum(sigma_dust, axis=2)/nphi                         # dimensions: (noutputs, ndust, nrad)   
        sigma_dust0_1D = np.sum(sigma_dust0, axis=2)/nphi                         # dimensions: (noutputs, ndust, nrad)   
        

        # ================ Calculate ring width ================

        for i,st in enumerate(stokes):
            # print("---------- Stokes = ", round(st,3))
            sigma_dust_st = sigma_dust_1D[i]
            # 1) Select data only within region close to planet
            # innerbound = rp + 0.11
            # outerbound = rp + (r_hill*r_hill_width) -0.1 # search for peak from rp to outerbound
            innerbound = 1.1
            outerbound = 1.3
            innerbound_i = np.argmin(np.abs(radii-innerbound))      # index (radial cell number) of lower bound of peak search
            outerbound_i = np.argmin(np.abs(radii-outerbound))      # index (radial cell number) of upper bound of peak search

            sigma_bound = sigma_dust_st[innerbound_i:outerbound_i]
            radii_bound = radii[innerbound_i:outerbound_i]

            peak_i_bound = find_ring_peak(radii_bound,sigma_bound)
            peak_i = peak_i_bound + innerbound_i

            # 2) Fit Gaussian to data to remove small variations
            params, covar = curve_fit(gaussian, radii_bound, sigma_bound/np.max(sigma_bound))
            afit, x0fit, bfit = params
            radii_arr = np.linspace(np.min(radii_bound), np.max(radii_bound),100)
            gaussian_fit = gaussian(radii_arr, afit, x0fit, bfit)
            ring_width = np.abs(4*bfit)     # ring_width = 4sigma
            ring_widths[s,i] = ring_width

            # plt.clf()
            # plt.plot(radii,sigma_dust_st/np.max(sigma_bound), c='k')
            # plt.scatter(radii,sigma_dust_st/np.max(sigma_bound), c='k', marker='x')
            # plt.plot(radii_arr, gaussian_fit, c='r')
            # plt.axvline(innerbound, c='k', linestyle='dashed')
            # plt.axvline(outerbound, c='k', linestyle='dashed')
            # plt.xlim(1,1.5)
            # plt.savefig(f"ring{s}_{i}.png")

    # 3) Plot ring width vs planet mass
    for i,st in enumerate(stokes):
        alpha_st = round(alpha/st, 4)
        ax.scatter(planet_masses, ring_widths[:,i])
        ax.plot(planet_masses, ring_widths[:,i], label=f"$\\alpha/St = {alpha_st}$")
    
    ax.set_xlabel("R")
    ax.set_ylabel("Ring width")
    ax.set_title(f"H/R = {hr0}")
    ax.legend()
    fig.savefig("./ring_widths.png")



    print(f"-------------------\nPlotting output {output} for {sim}\n=============")

    if plot_window:
        plt.show()

