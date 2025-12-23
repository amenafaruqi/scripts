import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
from scipy.interpolate import interp1d

from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D
import argparse
import re

plt.style.use('default')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 12})
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

import warnings
warnings.filterwarnings("ignore")

nstokes = 5
cm_st = plt.get_cmap('viridis')
cm_m = plt.get_cmap('plasma')
plain_clr = "k"
colour_cycler_st = cm_st(np.linspace(0, 1.01, 5))
colour_cycler_m = cm_m(np.linspace(0, 1.1, 7))

# ================== Fitting functions ==================

def gaussian(x, a, x0, b): 
    return a*np.exp(-0.5*((x-x0)/b)**2) 
    # return a*np.exp(-0.5*((x-x0)/(b+c(x-x0)))**2)   # skewed gaussian

def lorentzian(x, a, x0, gamma): 
    return (a*gamma/(2*np.pi))/((x-x0)**2 + (gamma/2)**2)

def find_ring_peak(sigma):
    peak_i = np.argmax(sigma)
    return peak_i

def find_ring_troughs(radii, sigma, i_peak):
    grad_sigma = np.gradient(sigma, radii)
    # gradients left of peak, in order of increasing distance from peak
    grad_left = grad_sigma[i_peak-1::-1]
    # gradients right of peak, in order of increasing distance from peak
    grad_right = grad_sigma[i_peak+1:]

    # check for change in gradient sign left of peak
    grad_li = 0
    while (grad_li < len(grad_left)-1) and (grad_left[grad_li]*grad_left[grad_li+1] >= 0):
        grad_li += 1
    left_trough_i = i_peak - grad_li

    # check for change in gradient sign right of peak
    grad_ri = 0
    while (grad_ri < len(grad_right)-1) and (grad_right[grad_ri]*grad_right[grad_ri+1] >= 0):
        grad_ri += 1
    right_trough_i = i_peak + grad_ri

    return left_trough_i,  right_trough_i


def calculate_Miso(hr, alpha, st, scaling="B18"):
    dlogPdlogR = f - sigmaslope - 2     # taken from eq. 9 from Bitsch et al. 2018, accounting for their s being -ve. 
    f_fit = ((hr/0.05)**3)*(0.34*(np.log10(0.001)/np.log10(alpha))**4 + 0.66)*(1-((dlogPdlogR+2.5)/6))
    alpha_st = alpha/st
    if scaling == "B18":
        M_iso = 25*f_fit
    elif scaling == "L14":
        M_iso = 20*f_fit
    Pi_crit = alpha_st/2
    Lambda = 0.00476/f_fit
    M_iso += Pi_crit/Lambda                  # PIM considering diffusion
    return M_iso


# ============ Plotting Functions ==============

def calculate_ring_widths():
    for i,st in enumerate(stokes):
        print("---------- Stokes = ", round(st,3))
        sigma_dust_st = sigma_dust_1D[i]

        # 1) Select data only within region close to planet
        innerbound = rp + (4*r_hill)
        outerbound = rp + (12*r_hill)  # search for peak from rp to outerbound
        # innerbound = 1.1
        # outerbound = 1.4
        innerbound_i = np.argmin(np.abs(radii-innerbound))      # index (radial cell number) of lower bound of peak search
        outerbound_i = np.argmin(np.abs(radii-outerbound))      # index (radial cell number) of upper bound of peak search

        sigma_bound = sigma_dust_st[innerbound_i:outerbound_i]
        radii_bound = radii[innerbound_i:outerbound_i]

        peak_i_bound = find_ring_peak(sigma_bound)
        peak_i = peak_i_bound + innerbound_i
        peak_r = radii_bound[peak_i_bound]
        init_guess = [1,peak_r,0.04]
        fit_lims = ((0.99,peak_r*0.999,-0.15),(1.01,peak_r*1.001,0.15))

        # 2) Fit Gaussian to data to remove small variations
        try:
            params, covar = curve_fit(gaussian, radii_bound, sigma_bound/np.max(sigma_bound), p0=init_guess, bounds=fit_lims)
            afit, x0fit, bfit = params
            radii_arr = np.linspace(np.min(radii_bound), np.max(radii_bound),100)
            gaussian_fit = gaussian(radii_arr, afit, x0fit, bfit)
            ring_width = np.abs(4*bfit)     # ring_width = 4sigma
        except RuntimeError:
            ring_width = np.nan
        ring_widths[s,i] = ring_width
        
        fig0, ax0 = plt.subplots(figsize=(7,5))
        ax0.cla()
        ax0.plot(radii,sigma_dust_st/np.max(sigma_bound), c=plain_clr)
        ax0.scatter(radii,sigma_dust_st/np.max(sigma_bound), c=plain_clr, marker='x')
        ax0.plot(radii_arr, gaussian_fit, c='r')
        ax0.axvline(innerbound, c=plain_clr, linestyle='dashed')
        ax0.axvline(outerbound, c=plain_clr, linestyle='dashed')
        ax0.set_xlim(1,1.5)
        ax0.set_ylim(0,np.max(sigma_bound/np.max(sigma_bound))*1.05)
        fig0.savefig(f"{plots_savedir}/rings_{mp}Me_{hr0}_{i}.png")


def fit_rings():
    for i,st in enumerate(stokes):
        print("---------- Stokes = ", round(st,3))
        sigma_dust_st = sigma_dust_1D[i]

        # 1) Select data only within region close to planet
        innerbound = rp + (4*r_hill)
        outerbound = rp + (15*r_hill)  # search for peak from rp to outerbound
        innerbound_i = np.argmin(np.abs(radii-innerbound))      # index (radial cell number) of lower bound of peak search
        outerbound_i = np.argmin(np.abs(radii-outerbound))      # index (radial cell number) of upper bound of peak search

        sigma_bound = sigma_dust_st[innerbound_i:outerbound_i]
        radii_bound = radii[innerbound_i:outerbound_i]

        # Find ring peak
        peak_i_bound = find_ring_peak(sigma_bound)
        peak_i = peak_i_bound + innerbound_i
        peak_r = radii_bound[peak_i_bound]

        # 2) Split data into left and right of peak and reflect to get a symmetric distribution
        left_sigma = sigma_dust_st[innerbound_i:peak_i+1]
        left_sigma_symm = np.concatenate((left_sigma, left_sigma[1::-1]))
        left_radii = radii[innerbound_i:innerbound_i+len(left_sigma_symm)]

        right_sigma = sigma_dust_st[peak_i:outerbound_i]
        right_sigma_symm = np.concatenate((right_sigma, right_sigma[1::-1]))
        right_radii = radii[outerbound_i-len(right_sigma_symm):outerbound_i]

        # 3) Fit Gaussian to data left of peak
        init_guess = [1,peak_r,0.01]
        fit_lims = ((0.99,peak_r*0.999,-0.1),(1.01,peak_r*1.001,0.1))

        try:
            params, covar = curve_fit(gaussian, left_radii, left_sigma_symm/np.max(left_sigma_symm), p0=init_guess, bounds=fit_lims)
            afit, x0fit, bfit = params
            radii_arr1 = np.linspace(np.min(radii_bound), peak_r,100)
            gaussian_fit = gaussian(radii_arr1, afit, x0fit, bfit)
            ring_width = 2*bfit*(2*np.log(2))**0.5     # ring_width = FWHM = 2b*sqrt(2ln2)
        except RuntimeError:
            ring_width = np.nan
        ring_widths[s,i,0] = ring_width

        # 4) Fit Lorentzian to data right of peak
        init_guess = [0.15,peak_r,0.01]
        fit_lims = ((0.1,peak_r*0.999,-0.5),(0.2,peak_r*1.001,0.5))

        try:
            params, covar = curve_fit(lorentzian, right_radii, right_sigma_symm/np.max(right_sigma_symm), p0=init_guess, bounds=fit_lims)
            afit, x0fit, gammafit = params
            radii_arr2 = np.linspace(peak_r, np.max(radii_bound),100)
            lorentzian_fit = lorentzian(radii_arr2, afit, x0fit, gammafit)
            ring_width = gammafit     # ring_width = FWHM = gamma
        except RuntimeError:
            ring_width = np.nan
        ring_widths[s,i,1] = ring_width

        # # 4.5) Fit Gaussian to data right of peak
        # init_guess = [1,peak_r,0.01]
        # fit_lims = ((0.99,peak_r*0.999,-0.15),(1.01,peak_r*1.001,0.15))

        # try:
        #     params, covar = curve_fit(gaussian, right_radii, right_sigma_symm/np.max(right_sigma_symm), p0=init_guess, bounds=fit_lims)
        #     afit, x0fit, bfit = params
        #     radii_arr2 = np.linspace(peak_r, np.max(radii_bound),100)
        #     gaussian_fit2 = gaussian(radii_arr2, afit, x0fit, bfit)
        #     ring_width = 2*bfit*(2*np.log(2))**0.5     # ring_width = FWHM = 2b*sqrt(2ln2)
        # except RuntimeError:
        #     ring_width = np.nan
        # ring_widths[s,i,0] = ring_width


        fig0, ax0 = plt.subplots(figsize=(7,5))
        ax0.cla()
        ax0.plot(radii,sigma_dust_st/np.max(sigma_bound), c=plain_clr, alpha=0.5)
        ax0.scatter(radii,sigma_dust_st/np.max(sigma_bound), c=plain_clr, marker='x', alpha=0.5)
        ax0.plot(radii_arr1, gaussian_fit, c='r', linestyle="dashed")
        ax0.plot(radii_arr2, lorentzian_fit, c='b', linestyle="dashed")
        # ax0.plot(radii_arr2, gaussian_fit2, c='b', linestyle="dashed")
        ax0.axvline(innerbound, c=plain_clr, linestyle='dashed')
        ax0.axvline(outerbound, c=plain_clr, linestyle='dashed')
        ax0.set_xlim(1,1.5)
        ax0.set_ylim(0,np.max(sigma_bound/np.max(sigma_bound))*1.05)
        fig0.savefig(f"{plots_savedir}/rings_{mp}Me_{hr0}_{i}.png")


def calculate_ring_masses():
    for i,st in enumerate(stokes):
        print("---------- Stokes = ", round(st,3))
        sigma_dust_st = sigma_dust_1D[i]

        # 1) Select data only within region close to planet
        innerbound = rp + (4*r_hill)
        outerbound = rp + (10*r_hill)  # search for peak from rp to outerbound
        innerbound_i = np.argmin(np.abs(radii-innerbound))      # index (radial cell number) of lower bound of peak search
        outerbound_i = np.argmin(np.abs(radii-outerbound))      # index (radial cell number) of upper bound of peak search

        sigma_bound = sigma_dust_st[innerbound_i:outerbound_i]
        radii_bound = radii[innerbound_i:outerbound_i]

        # 2) Identify ring peaks and troughs in gas
        peak_i_bound = find_ring_peak(sigma_bound)
        peak_i = peak_i_bound + innerbound_i

        l_trough_bound_i, r_trough_bound_i = find_ring_troughs(radii_bound, sigma_bound, peak_i_bound)
        l_trough_i = l_trough_bound_i + innerbound_i
        r_trough_i = r_trough_bound_i + innerbound_i

        # 3) Sum up mass between left and right trough
        # ring_mass = np.sum(dust_mass[i,l_trough_i:r_trough_i])
        ring_mass = np.sum(dust_mass[i,int(inner_edge_i[s,i]):int(outer_edge_i[s,i])])
        ring_masses[s,i] = ring_mass


def calculate_dust_gas_ratios():
    dust_mass_tot = np.sum(dust_mass, axis=0)  # sum over all St
    dgrs.append(dust_mass_tot/gas_mass)
        


def calculate_ring_edges():
    grad_mins=[2e-7,1e-7,2e-7,2e-7, 2e-7]
    for i,st in enumerate(stokes):
        print("---------- Stokes = ", round(st,3))
        sigma_dust_st = sigma_dust_1D[i]

        # 1) Select data only within region close to planet
        innerbound = rp + (3*r_hill)
        outerbound = rp + (20*r_hill)  # search for peak from rp to outerbound
        innerbound_i = np.argmin(np.abs(radii-innerbound))      # index (radial cell number) of lower bound of peak search
        outerbound_i = np.argmin(np.abs(radii-outerbound))      # index (radial cell number) of upper bound of peak search

        sigma_bound = sigma_dust_st[innerbound_i:outerbound_i]
        radii_bound = radii[innerbound_i:outerbound_i]

        # 2) Identify ring peak
        peak_i_bound = find_ring_peak(sigma_bound)
        # peak_i = peak_i_bound + innerbound_i
        peak_r = radii_bound[peak_i_bound]

        # 3) Locate ring edges in either direction based on sigma drop
        sigma_peak = sigma_bound[peak_i_bound]
        grad_bound = np.gradient(sigma_bound)
        # if s==0 or i > 2 :
        #     grad_bound_min = np.min(np.abs(grad_bound)) * 25
        # else:
        #     grad_bound_min = np.min(np.abs(grad_bound)) * 100
        sigma_bound = sigma_bound/sigma_peak
        left_edge_i = peak_i_bound - 5
        right_edge_i = peak_i_bound + 5
        threshold = 0.1
        # Left edge first
        # Iterate while sigma > threshold*sigma_peak AND gradient > than some value
        # while (sigma_bound[left_edge_i] > threshold) and (left_edge_i > 1) and (np.abs(grad_bound[left_edge_i]) > 1e-5):    # 1e-23 or 1e-7 for h=0.05
        #     left_edge_i -= 1
        # while (sigma_bound[right_edge_i] > threshold) and (right_edge_i < outerbound_i-innerbound_i-1) and (np.abs(grad_bound[right_edge_i]) > 1-5):  #2e-7
        #     right_edge_i += 1


        while (sigma_bound[left_edge_i] > threshold) and (left_edge_i > 1) and (np.all(np.abs(grad_bound[left_edge_i:left_edge_i+2])) > 0.5):    # 1e-23 or 1e-7 for h=0.05
            left_edge_i -= 1
        while (sigma_bound[right_edge_i] > threshold) and (right_edge_i < outerbound_i-innerbound_i-1) and (np.all(np.abs(grad_bound[right_edge_i-2:right_edge_i])) > 0.99):  #2e-7
            right_edge_i += 1


        # while (left_edge_i > 1) and (np.all(np.abs(grad_bound[left_edge_i:left_edge_i+3]) > 0.1)):    # 1e-23 or 1e-7 for h=0.05
        #     left_edge_i -= 1
        # while (right_edge_i < outerbound_i-innerbound_i-1) and (np.all(np.abs(grad_bound[right_edge_i-3:right_edge_i])) > 0.2):  #2e-7
        #     right_edge_i += 1


        # Define edges by FWHM
        # threshold = 0.3
        # while (sigma_bound[left_edge_i] > threshold) and (left_edge_i > 1):    # 1e-23 or 1e-7 for h=0.05
        #     left_edge_i -= 1
        # while (sigma_bound[right_edge_i] > threshold) and (right_edge_i < outerbound_i-innerbound_i-1):  #2e-7
        #     right_edge_i += 1

        # left_edge = np.interp(sigma_peak*threshold, sigma_bound[left_edge_i-1:left_edge_i+2], radii_bound[left_edge_i-1:left_edge_i+2])
        # right_edge = np.interp(sigma_peak*threshold, sigma_bound[right_edge_i-1:right_edge_i+2], radii_bound[right_edge_i-1:right_edge_i+2])
        
        # Adjustments for h = 0.05
        # # if (i==0) and (s==5):    # only rmass plots
        # #     right_edge_i -= 3
        # #     left_edge_i += 3
        # if (i==2) and (s==4):
        #     right_edge_i -= 3
        #     left_edge_i -= 2
        # if (i==2) and (s==2):
        #     right_edge_i += 1
        # if (i==3) and (s==0):
        #     right_edge_i -= 7
        # if (i==3) and (s==1):
        #     right_edge_i -= 17
        # if (i==3) and (s==2):
        #     right_edge_i -= 5
        # if (i==3) and (s==4):
        #     right_edge_i -= 8
        # if (i==4) and (s==0):
        #     left_edge_i += 5
        # if (i==4) and (s==1):
        #     left_edge_i -= 5
        #     right_edge_i -= 10
        # if (i==4) and (s==2):
        #     right_edge_i -= 20
        #     left_edge_i -= 5
        # if (i==4) and (s==3):
        #     right_edge_i -= 28
        # if (i==4) and (s==4):
        #     right_edge_i -= 30
        #     left_edge_i -= 3
        # if (i==4) and (s==5):
        #     right_edge_i -= 40


        # # Adjustments for h = 0.06    
        # if (i==2) and (s==0):
        #     left_edge_i -= 2
        #     right_edge_i += 4
        # if (i==2) and (s==1):
        #     left_edge_i -= 2
        #     right_edge_i += 2
        # if (i==2) and (s==2):
        #     right_edge_i += 1
        # if (i==2) and (s==3):
        #     left_edge_i -= 3
        # if (i==3) and (s==0):
        #     right_edge_i -= 7
        # if (i==3) and (s==1):
        #     right_edge_i -= 20
        # if (i==3) and (s==2):
        #     right_edge_i -= 5
        # if (i==4) and (s==0):
        #     left_edge_i += 3
        # if (i==4) and (s==1):
        #     right_edge_i -= 10
        # if (i==4) and (s==2):
        #     right_edge_i -= 25
        # if (i==4) and (s==3):
        #     right_edge_i -= 35
        # if (i==4) and (s==4):
        #     right_edge_i -= 40
        # if (i==4) and (s==5):
        #     right_edge_i -= 52

        # Adustments for h = 0.07
        if (i==1) and (s==0):
            right_edge_i += 3
        if (i==1) and (s==2):
            right_edge_i += 2
        if (i==2) and (s==0):
            left_edge_i -= 4
            right_edge_i += 3
        if (i==2) and (s==1):
            right_edge_i += 4
        if (i==3) and (s==0):
            right_edge_i -= 10
        if (i==3) and (s==1):
            right_edge_i -= 1
            left_edge_i -= 1
        if (i==3) and (s==2):
            right_edge_i -= 4
        if (i==3) and (s==3):
            right_edge_i -= 7
        if (i==4) and (s==1):
            right_edge_i -= 6
        if (i==4) and (s==2):
            right_edge_i -= 18
        if (i==4) and (s==3):
            right_edge_i -= 30
        if (i==4) and (s==4):
            right_edge_i -= 35
        if (i==4) and (s==5):
            right_edge_i -= 45


        # Plots for debugging
        ax_edge[s,i].plot(radii, sigma_dust_st/sigma_peak, color="k")
        if s == 0:
            if i == 2:
                ax_edge[s,i].set_title(f"{int(planet_masses[s])} M$_\oplus$ \n St={round(st,3)}")
            else:
                ax_edge[s,i].set_title(f"St={round(st,3)}")
        else:
            if i == 2:
                ax_edge[s,i].set_title(f"{int(planet_masses[s])} M$_\oplus$")

        ax_edge[s,i].scatter(radii_bound[left_edge_i], sigma_bound[left_edge_i], color='r')
        ax_edge[s,i].scatter(radii_bound[right_edge_i], sigma_bound[right_edge_i], color='r')
        ax_edge[s,i].scatter(radii_bound[peak_i_bound], sigma_bound[peak_i_bound], color='b')
        
        # sigma_data[s,i] = sigma_dust_st/sigma_peak
        # edge_data[s,i] = [
        #     [radii_bound[left_edge_i], sigma_bound[left_edge_i]],
        #     [radii_bound[right_edge_i], sigma_bound[right_edge_i]],
        #     [[peak_i_bound], sigma_bound[peak_i_bound]]
        #     ]
        
        if i==0:
            ax_edge[s,i].set_ylabel("$\\Sigma_{dust}/\\Sigma_{peak}$")
        if s == 5:
            ax_edge[s,i].set_xlabel("Radius ($r_{p}$)")
        ax_edge[s,i].set_xlim(np.min(radii_bound)*0.95, np.max(radii_bound)*0.9)
        ax_edge[s,i].set_ylim(-0.1, 1.2)
        # 3) Populate array with edge locations
        # inner_edges[s,i] = left_edge + innerbound
        # outer_edges[s,i] = right_edge + innerbound
        inner_edges[s,i] = radii_bound[left_edge_i] - peak_r
        outer_edges[s,i] = radii_bound[right_edge_i] - peak_r
        inner_edge_i[s,i] = left_edge_i + innerbound_i    # index of inner edge radial cell in full radii array
        outer_edge_i[s,i] = right_edge_i + innerbound_i
    fig_edge.tight_layout()
    fig_edge.savefig(f"{plots_savedir}/redge_locs_{hr0}.png", dpi=200)
    
def calculate_pressure_gradients():
    # 1) Select data only within region close to planet
    innerbound = rp + (4*r_hill)
    outerbound = rp + (12*r_hill)  # search for peak from rp to outerbound
    innerbound_i = np.argmin(np.abs(radii-innerbound))      # index (radial cell number) of lower bound of peak search
    outerbound_i = np.argmin(np.abs(radii-outerbound))      # index (radial cell number) of upper bound of peak search

    sigma_bound = sigma_gas_1D[innerbound_i:outerbound_i]
    radii_bound = radii[innerbound_i:outerbound_i]

    # 2) Identify ring peaks and troughs in gas
    peak_i_bound = find_ring_peak(sigma_bound)
    peak_i = peak_i_bound + innerbound_i

    l_trough_bound_i, r_trough_bound_i = find_ring_troughs(radii_bound, sigma_bound, peak_i_bound)
    l_trough_i = l_trough_bound_i + innerbound_i
    r_trough_i = r_trough_bound_i + innerbound_i

    # Plot gas rings to check pressure gradient calculations
    fig0, ax0 = plt.subplots(figsize=(7,5))
    ax0.cla()
    ax0.plot(radii,sigma_gas_1D/np.max(sigma_bound), c='k')
    # ax0.scatter(radii,sigma_gas_1D/np.max(sigma_bound), c='k', marker='x')
    ax0.scatter(radii[l_trough_i],(sigma_gas_1D/np.max(sigma_bound))[l_trough_i], c='r', marker='o')
    ax0.scatter(radii[r_trough_i],(sigma_gas_1D/np.max(sigma_bound))[r_trough_i], c='r', marker='o')
    ax0.scatter(radii[peak_i],(sigma_gas_1D/np.max(sigma_bound))[peak_i], c='y', marker='o')
    ax0.axvline(innerbound, c='k', linestyle='dashed')
    ax0.axvline(outerbound, c='k', linestyle='dashed')
    ax0.set_xlim(1,1.6)
    ax0.set_ylim(0,np.max(sigma_bound/np.max(sigma_bound))*1.2)

    # 3) Calculate dP/dr
    hr_bound = hr0*(radii_bound**f)
    pressure = 4*(np.pi**2)*hr_bound*sigma_bound*(radii_bound**(-2))
    dpdr = np.abs(np.gradient(pressure, radii_bound))
    ax0.plot(radii_bound, dpdr/np.max(dpdr), c='r', label="$| \partial P/ \partial r |$")
    ax0.plot(radii_bound, pressure/np.max(pressure), c='cyan', label="P")
    # ax0.scatter(radii_bound, dpdr/np.max(dpdr), c='r')
    ax0.legend()

    # Find steepest point interior/exterior to ring peak
    dpdr_int = dpdr[l_trough_bound_i+3:peak_i_bound-1]
    dpdr_ext = dpdr[peak_i_bound+1:r_trough_bound_i-3]

    max_dpdr_int = np.median(dpdr_int)
    max_dpdr_ext = np.median(dpdr_ext)
    ax0.plot(radii_bound[l_trough_bound_i+3:peak_i_bound-1], dpdr_int/np.max(dpdr), c='b')
    ax0.plot(radii_bound[peak_i_bound+1:r_trough_bound_i-3], dpdr_ext/np.max(dpdr), c='g')
    fig0.savefig(f"{plots_savedir}/gasrings_{mp}Me_{hr0}_.png")

    dpdrs_ext[s] = max_dpdr_ext
    dpdrs_in[s] = max_dpdr_int


def calculate_dlogpdlogr():
    h = hr0*radii**(f+1)
    sigma = sigma_gas_1D
    pressure = h*sigma*(radii**(-3))
    lnp = np.log(pressure)
    lnr = np.log(radii)
    dlogpdlogr = np.gradient(lnp, lnr)
    dlogpdlogrs.append(dlogpdlogr)


def calculate_pressure():
    hr = hr0*(radii**f)
    pressure = 4*(np.pi**2)*hr*sigma_gas_1D*(radii**(-2))
    pressure_arr.append(pressure)


def calculate_dust_flux():
    v_r = v_gas[1,:,:]                     # dimensions: (nrad, nphi)
    R, PHI = np.meshgrid(radii, phis, indexing="ij")   # dimensions: (nrad, nphi)
    gas_flux.append(2*np.pi*v_r*sigma_gas*R*2.1e10)

    for n in np.arange(ndust):
        v_r = v_dust[n,1,:,:]                     # dimensions: (nrad, nphi)
        dust_flux.append(2*np.pi*v_r*sigma_dust[n]*R*2.1e10)


def calculate_ring_peaks():
    for i,st in enumerate(stokes):
        print("---------- Stokes = ", round(st,3))
        sigma_dust_st = sigma_dust_1D[i]

        # 1) Select data only within region close to planet
        innerbound = rp + (4*r_hill)
        outerbound = rp + (12*r_hill)  # search for peak from rp to outerbound
        # innerbound = 1.1
        # outerbound = 1.4
        innerbound_i = np.argmin(np.abs(radii-innerbound))      # index (radial cell number) of lower bound of peak search
        outerbound_i = np.argmin(np.abs(radii-outerbound))      # index (radial cell number) of upper bound of peak search

        sigma_bound = sigma_dust_st[innerbound_i:outerbound_i]
        radii_bound = radii[innerbound_i:outerbound_i]

        peak_i_bound = find_ring_peak(sigma_bound)
        peak_r = radii_bound[peak_i_bound]
        r_peaks[s,i] = peak_r


def calculate_vdrift():
    v_drifts_s = []
    for i,st in enumerate(stokes):
        print("---------- Stokes = ", round(st,3))
        hr = hr0*(radii**f)
        pressure = 4*(np.pi**2)*hr*sigma_gas_1D*(radii**(-2))
        dpdr = np.gradient(pressure, radii)
        omega = radii**(-3/2)
        v_drift = -dpdr*(hr*radii/(omega*sigma_gas_1D))/(st+(1/st))
        v_drifts_s.append(v_drift)

    v_drifts.append(v_drifts_s)


def calculate_ppf():
    delta = 1e-5
    # Calculate mass avged St and Sigma dust
    mass_avg_sum_St = np.sum([stokes[i]*W_St[i] for i in range(len(stokes))])
    St_avg = mass_avg_sum_St/len(stokes)
    print(St_avg)
    # print("St_avg = ", St_avg)
    Sigma_avg = np.sum([sigma_dust_1D[i,:]*W_St[i] for i in range(len(stokes))], axis=0) # divide by ndust then multiply by ndust
    print(radii[list(Sigma_avg[:400]).index(np.min(Sigma_avg[:400]))])
    # Calculate stability criterion for dust clump
    h = hr0*(radii**(1+f))    # h as a function of r
    Qp = ((delta/St_avg)**0.5)*h/(np.pi*(radii**3)*Sigma_avg)
    # print("Qp = ", Qp)
    ppf = 1/(1 + np.exp(10*(Qp-0.75)))  # dimensions = nrad
    print(np.max(ppf[:400]))
    ppfs.append(ppf)


# ============== Read in data from models ==============

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fit dust rings', prefix_chars='-')

    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/amena/scratch/simulations/dusty_fargo/h0.05"],help="working directory containing simulations")
    parser.add_argument('-sims', metavar='sim', type=str, nargs="*", default=["10Me"] ,help="simulation directory containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="/home/amena/scratch/images/" ,help="directory to save plots to")
    parser.add_argument('-o', metavar='outputs',default=[600], type=int, nargs="*" ,help="outputs to plot")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-style', metavar='style', type=str, nargs="*", default=["publication"] ,help="style sheet to apply to plots")
    parser.add_argument('-plots', metavar='plots',default=["rwidth"], type=str, nargs="*" ,help="plots to produce")

    args = parser.parse_args()
    output = args.o[0]
    wd = args.wd[0]
    sims = args.sims
    plot_window = args.plot_window
    plots_savedir = args.savedir
    style = args.style
    plots = args.plots    # opts: rwidth, dpdr, dflux

    if plot_window:
        mpl.use('TkAgg') 
    else:
        for sty in style:
            plt.style.use([f"../styles/{sty}.mplstyle"])
            if "darkbg" in sty:
                plots_savedir = plots_savedir+"/darkbg/"
                plain_clr = "w"

    # ------------------------------------------------------
    
    # Initialise axes and arrays for data to plot
    planet_masses = np.zeros((len(sims)))
    if "rwidth" in plots:
        fig_rw, ax_rw = plt.subplots(figsize=(7,6))
        ring_widths = np.zeros((len(sims),5,2))
        fig_edge,ax_edge = plt.subplots(ncols=5, nrows=6, figsize=(15,15), sharex=True, sharey=True)
        inner_edges = np.zeros((len(sims),5))
        outer_edges = np.zeros((len(sims),5))
        outer_edges = np.zeros((len(sims),5))
        inner_edge_i = np.zeros((len(sims),5))
        outer_edge_i = np.zeros((len(sims),5))

    if "rmass" in plots:
        fig_rm, ax_rm = plt.subplots(figsize=(10,6))
        ring_masses = np.zeros((len(sims),5))

    if "dgr" in plots:
        fig_dgr, ax_dgr = plt.subplots(figsize=(8,6))
        dgrs = []

    if "redge" in plots:
        fig_re, ax_re = plt.subplots(figsize=(12,7))
        fig_edge,ax_edge = plt.subplots(ncols=5, nrows=6, figsize=(15,15), sharex=True)
        inner_edges = np.zeros((len(sims),5))
        outer_edges = np.zeros((len(sims),5))
        inner_edge_i = np.zeros((len(sims),5))
        outer_edge_i = np.zeros((len(sims),5))

    if "dpdr" in plots:
        fig_p, ax_p = plt.subplots(figsize=(6,5))
        fig_plog, ax_plog = plt.subplots(figsize=(6,5))
        fig_press, ax_press = plt.subplots(figsize=(6,5))
        dpdrs_in = np.zeros((len(sims)))
        dpdrs_ext = np.zeros((len(sims)))
        pressure_arr = []
        dlogpdlogrs = []
    
    if "flux" in plots:
        fig_f, ax_f = plt.subplots(figsize=(12,15), ncols=5, nrows=6)
        fig_f1d, ax_f1d = plt.subplots(figsize=(20,5), ncols=6)
        gas_flux = []
        dust_flux = []

    if "rpeak" in plots:
        fig_peak, ax_peak = plt.subplots(figsize=(6,5))
        r_peaks = np.zeros((len(sims),5))

    if "vdrift" in plots:
        fig_v, ax_v = plt.subplots(figsize=(8,5))
        r_peaks = np.zeros((len(sims),5))
        v_drifts = []
    
    if "ppf" in plots:
        fig_pf, ax_pf = plt.subplots(figsize=(12,7))
        ppfs =  []

    # ------------------------------------------------------
    # sigma_data = np.zeros((len(sims),5))
    # edge_data = np.zeros((len(sims),5))

    # Iterate through models and calculate ring width for each St
    for s, sim in enumerate(sims):
        simdir = f"{wd}/{sim}/"
        print(f"Reading variables for {simdir}...")
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
        sigmaslope = float(params_dict['SIGMASLOPE'])
        sigma0 = float(params_dict['SIGMA0'])
        ymin = float(params_dict['YMIN'])
        ymax = float(params_dict['YMAX'])
        f = float(params_dict['FLARINGINDEX'])
        spacing = str(params_dict['SPACING'])
        omegaframe = float(params_dict['OMEGAFRAME'])
        max_stokes = float(params_dict['STOKES'])
        stokes = np.logspace(np.log10(max_stokes),np.log10(max_stokes*10**(-2)),ndust)   # hardcodes St range to be max_stokes - max_stokes*1e-2
        
        # Calculate timing parameters
        dt_orbits = int(float(params_dict['DT'])/(2*np.pi))   # 2pi = 1 orbit = 1 yr
        ninterm = float(params_dict['NINTERM'])               # number of dts between outputs
        dt_outputs = dt_orbits*ninterm                        # time between outputs
        
        # Calculate planet location and Hill radius
        planet_data = np.unique(np.loadtxt(f"{simdir}/planet0.dat"), axis=0)
        xp, yp = planet_data[output,1], planet_data[output,2]
        rp = ((xp**2) + (yp**2))**0.5
        mp = planet_data[output,7]
        r_hill = rp*(mp/3)**(1/3)
        mp = int(round(mp/(3.0027e-6),0))     # convert to earth masses
        planet_masses[s] = mp

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
        gasfile = f"gasdens{output}.dat" 
        sigma_gas = np.fromfile(simdir+gasfile).reshape(nrad,nphi)

        sigma_dust = np.zeros((ndust, nrad, nphi))
        for n in np.arange(ndust):
            dust_file = f"dustdens{n}_{output}.dat"
            sigma_dust[n] = np.fromfile(simdir+dust_file).reshape(nrad,nphi)
            dust_file0 = f"dustdens{n}_0.dat"

        # Average over all phi
        sigma_gas_1D = np.sum(sigma_gas, axis=1)/nphi                        # dimensions: (nrad) 
        sigma_dust_1D = np.sum(sigma_dust, axis=2)/nphi                      # dimensions: (ndust, nrad)   

        # Calculate dust masses in units of disc mass
        m_disc0 = 2*np.pi*sigma0*(ymax-ymin)   # total disc (gas) mass at t=0
        # Calculate MRN weighting coefficient for mass
        W_St = [(st**0.5)/np.sum(stokes**0.5) for st in stokes]

        dust_mass = np.array([2*np.pi*radii*sigma_dust_1D[n,:]*delta_r*W_St[n] for n in range(ndust)])/m_disc0   # ndust x nrad
        gas_mass = 2*np.pi*radii*sigma_gas_1D*delta_r/m_disc0   # nrad


        if "flux" in plots:
            v_dust = np.zeros((ndust, 2, nrad, nphi))   # additional dimension of 2 for x and y velocity
            v_gas = np.zeros((2, nrad, nphi))           # additional dimension of 2 for x and y velocity
            gas_file_x = f"gasvx{output}.dat"
            gas_file_y = f"gasvy{output}.dat"
            v_gas[0] = np.fromfile(simdir+gas_file_x).reshape(nrad,nphi)       # vx
            v_gas[1] = np.fromfile(simdir+gas_file_y).reshape(nrad,nphi)       # vy

            for n in np.arange(ndust):
                dust_file_x = f"dustvx{n}_{output}.dat"
                dust_file_y = f"dustvy{n}_{output}.dat"
                v_dust[n,0] = np.fromfile(simdir+dust_file_x).reshape(nrad,nphi)   # vx (azimuthal v)
                v_dust[n,1] = np.fromfile(simdir+dust_file_y).reshape(nrad,nphi)   # vy (radial v)


        # =========== Compute values to plot and populate arrays ============
        if "rwidth" in plots:
            print(f"Calculating ring widths for {mp} Mearth... \n ~ ~ ~ ~ ~")
            # calculate_ring_widths()
            # fit_rings()
            calculate_ring_edges()
            # np.savetxt(f"{hr0}_1dsigma.txt", [radii, sigma_data, edge_data])
        if "rmass" in plots:
            print(f"Calculating ring masses for {mp} Mearth... \n ~ ~ ~ ~ ~")
            calculate_ring_masses()
        if "dgr" in plots:
            print(f"Calculating dust-gas ratio for {mp} Mearth... \n ~ ~ ~ ~ ~")
            calculate_dust_gas_ratios()
        if "redge" in plots:
            print(f"Calculating ring edge locations for {mp} Mearth... \n ~ ~ ~ ~ ~")
            calculate_ring_edges()
        if "dpdr" in plots:
            print(f"Calculating pressure gradients for {mp} Mearth... \n ~ ~ ~ ~ ~")
            calculate_pressure()
            calculate_pressure_gradients()
            calculate_dlogpdlogr()
        if "flux" in plots:
            print(f"Calculating flux for {mp} Mearth... \n ~ ~ ~ ~ ~")
            calculate_ring_edges()
            calculate_dust_flux()
        if "rpeak" in plots:
            print(f"Calculating ring peak for {mp} Mearth... \n ~ ~ ~ ~ ~")
            calculate_ring_peaks()
        if "vdrift" in plots:
            print(f"Calculating drift velocities for {mp} Mearth... \n ~ ~ ~ ~ ~")
            calculate_vdrift()
            calculate_ring_peaks()
        if "ppf" in plots:
            print(f"Calculating planetesimal formation probability for {mp} Mearth... \n ~ ~ ~ ~ ~")
            calculate_ppf()


    # ================ Generate chosen plots ================
    print(f"-------------------\nPlotting output {output} for {sims}\n=============")

    if "rwidth" in plots:
        # Plot ring width vs planet mass
        print("Plotting ring width against planet mass....")
        for i,st in enumerate(stokes):
            colour = colour_cycler_st[i]
            # alpha_st = round(alpha/st, 4)
            # ax_rw.scatter(planet_masses*(3e-6)/(hr0**3), ring_widths[:,i], color=colour)
            # ax_rw.plot(planet_masses*(3e-6)/(hr0**3), ring_widths[:,i], color=colour, label=f"$\\alpha/St = {alpha_st}$" )
            ax_rw.scatter(planet_masses, outer_edges[:,i]-inner_edges[:,i], color=colour)
            ax_rw.plot(planet_masses, outer_edges[:,i]-inner_edges[:,i], color=colour, label=f"St = {round(st,3)}" )
            M_iso = calculate_Miso(hr0, alpha, st, scaling="B18")
            ax_rw.axvline(M_iso, linestyle='dotted', color=colour)

        # secax = ax_rw.secondary_xaxis('top', functions=(lambda x: x*(3e-6)/(hr0**3), lambda x: x*(hr0**3)/(3e-6)))
        # secax.set_xlabel('Planet mass (M$_{\\rm th}$)')
        # ax_rw.set_ylim(-0.05,0.5)
        ax_rw.set_xlabel("Planet mass (M$_\oplus$)")
        ax_rw.set_ylabel("Ring width ($r_{p}$)")
        ax_rw.set_title(f"$h_{{p}}/r_{{p}}$ = {hr0}", pad=10)
        ax_rw.legend(loc="best")
        fig_rw.tight_layout()
        fig_rw.savefig(f"{plots_savedir}/ring_widths_{hr0}.png", dpi=200)

    if "rmass" in plots:
        # Plot ring mass vs planet mass
        print("Plotting ring mass against planet mass....")
        # ring_masses[-1,0] -= 0.0005
        # print(ring_masses)
        # np.savetxt(f"rmass_{hr0}.txt", (planet_masses, np.sum(ring_masses, axis=1)))
        for i,st in enumerate(stokes):
            colour = colour_cycler_st[i]
            # alpha_st = round(alpha/st, 4)
            ax_rm.scatter(planet_masses, ring_masses[:,i], color=colour)
            ax_rm.plot(planet_masses, ring_masses[:,i], color=colour, label=f"St = {round(st, 3)}" )
            M_iso = calculate_Miso(hr0, alpha, st, scaling="B18")
            ax_rm.axvline(M_iso, linestyle='dotted', color=colour)

        ax_rm.set_xlabel("Planet mass (M$_\oplus$)")
        ax_rm.set_ylabel("Ring mass (M$_{\\rm disc}$)")
        ax_rm.set_title(f"$h_{{p}}/r_{{p}}$ = {hr0}", pad=10)
        # ax_rm.set_title(f"H/R = {hr0}")
        ax_rm.legend(loc="upper right", bbox_to_anchor=(0.98,0.95))
        # Secondary axis showing mass in Mearth
        secaxy = ax_rm.secondary_yaxis('right', functions=(lambda x: x*0.01*333000, lambda x: x/0.01*333000))        
        secaxy.set_ylabel('Ring mass for a 0.01M$_\odot$ disc (M$_\oplus$)')
        # secaxx = ax_rm.secondary_xaxis('top', functions=(lambda x: x*(3e-6)/(hr0**3), lambda x:x*(hr0**3)/(3e-6) ))
        # secaxx.set_xlabel('Planet mass  (M$_{\\rm th}$)')
    
        fig_rm.tight_layout()
        fig_rm.savefig(f"{plots_savedir}/ring_masses_{hr0}.png", dpi=200)

    if "dgr" in plots:
        # Plot dust-gas ratio vs r
        print("Plotting dust-gas ratio for different planet masses....")
        ax_dgr.plot(radii, [1e-2]*len(radii), color="lightgrey", linestyle="dashed")
        ax_dgr.plot(radii, [1]*len(radii), color="lightgrey", linestyle="dotted")
        for s,sim in enumerate(sims):
            label = str(planet_masses[s]) + "$M_\oplus$"
            colour = colour_cycler_m[s]
            ax_dgr.plot(radii, dgrs[s], color=colour, label=label)

        ax_dgr.set_ylabel("Dust-to-gas ratio")
        ax_dgr.set_xlabel("Radius ($r_{p}$)")
        ax_dgr.set_yscale("log")
        ax_dgr.set_ylim(5e-6, 10)
        ax_dgr.set_xlim(0.3,2.2)
        ax_dgr.legend(loc="upper right")
        ax_dgr.set_title(f"$h_{{p}}/r_{{p}}$ = {hr0}", pad=10)

    
        fig_dgr.tight_layout()
        fig_dgr.savefig(f"{plots_savedir}/dust_gas_ratio_{hr0}.png", dpi=200)


    if "dpdr" in plots:
        # Plot pressure against R
        print("Plotting pressure profile.... \n")
        for s, sim in enumerate(sims):
            colour = colour_cycler_m[s]
            mp = re.search(r"(\d+)Me", sim).group(1)
            ax_press.plot(radii, pressure_arr[s], color=colour, label=f"{mp}$M_\oplus$")
        ax_press.set_yscale("log")
        ax_press.set_xlim(0.4,2)
        ax_press.set_ylim(2e-4,4e-2)
        ax_press.set_xlabel("r")
        ax_press.set_ylabel("Gas pressure")
        ax_press.legend()
        fig_press.savefig(f"{plots_savedir}/pressure_profile_{hr0}.png")

        # Plot avg pressure gradient vs planet mass
        print("Plotting dP/dr against planet mass.... \n")
        M_iso_B18_upper = calculate_Miso(hr0, alpha, st=0.002, scaling="L14")
        M_iso_B18_lower = calculate_Miso(hr0, alpha, st=0.2, scaling="L14")
        ax_p.fill_betweenx(x1=M_iso_B18_lower, x2=M_iso_B18_upper, y=np.linspace(0,0.01,10), color='lightgrey',  interpolate=True, alpha=.4, label="$M_{iso}$")

        ax_p.scatter(planet_masses, dpdrs_ext, c='lawngreen')
        ax_p.plot(planet_masses, dpdrs_ext, c='lawngreen', label = "exterior to ring peak")
        ax_p.scatter(planet_masses, dpdrs_in, c='deepskyblue')
        ax_p.plot(planet_masses, dpdrs_in, c='deepskyblue', label="interior to ring peak")
        # ax_p.plot(planet_masses, dpdrs_in/dpdrs_ext, c='orange', label="interior to exterior $\partial P/\partial r$ ratio")
        # ax_p.set_ylim(-0.05,0.7)
        # M_iso_L14 = calculate_Miso(hr0, alpha, st=0.1, scaling="L14")
        # ax_p.axvline(M_iso_B18, linestyle="dotted", color=plain_clr, label="$M_{iso}$ (Bitsch et al. 2018)")
        # ax_p.axvline(M_iso_L14, linestyle="dashed", color=plain_clr, label="$M_{iso}$ (Lambrechts et al. 2014)")
        ax_p.set_xlabel("Planet mass (M$_\oplus$)")
        ax_p.set_ylabel("$ \\rm{max} (|\partial P/\partial r |)$")
        ax_p.set_title(f"$h_{{p}}/r_{{p}}$ = {hr0}", pad=10)
        ax_p.set_ylim(0.85*np.min(dpdrs_in), 1.05*np.max(dpdrs_in))
        ax_p.legend()
        fig_p.tight_layout()
        fig_p.savefig(f"{plots_savedir}/pressure_grad_{hr0}.png", dpi=200)

        print("Plotting dlogP/dlogr for each planet mass....")
        for s, sim in enumerate(sims):
            colour = colour_cycler_m[s]
            mp = re.search(r"(\d+)Me", sim).group(1)
            ax_plog.plot(radii, dlogpdlogrs[s], color=colour, label=f"{mp}M$_\oplus$")
        ax_plog.set_xlabel("r (AU)")
        ax_plog.set_ylabel("$\partial \ln P/ \partial \ln r$")
        ax_plog.set_title(f"H/R = {hr0}")
        ax_plog.set_xlim(0.3,2.0)
        ax_plog.set_ylim(-18,15)
        ax_plog.fill_between(x=radii, y1=0, y2=-20, color='lightgrey',  interpolate=True, alpha=.5)
        ax_plog.axhline(-2.75, linestyle="dashed", color=plain_clr, label="$\\frac{\partial \ln P}{\partial \ln r}\\rvert_{t=0}$")
        ax_plog.legend(loc="upper right")
        fig_plog.tight_layout()
        fig_plog.savefig(f"{plots_savedir}/dlogpdlogr_{hr0}.png", dpi=200)
    
    if "flux" in plots:
        for s, sim in enumerate(sims):
            R, PHI = np.meshgrid(radii, phis, indexing="ij")   # dimensions: (nrad, nphi)
            x = R*np.cos(PHI)
            y = R*np.sin(PHI)
            mp = re.search(r"(\d+)Me", sim).group(1)

            # lims = [[-1e-6,1e-6],[-5e-6,5e-6],[-5e-6,5e-6],[-5e-6,5e-6],[-5e-6,5e-6]]
            # lts = [1e-40,           1e-35,      1e-15,      1e-12,          1e-12]
            for n in np.arange(ndust):
                ax = ax_f[s,n]
                flux = dust_flux[s*ndust:(s+1)*ndust][n] 
                print("------------- St = ", stokes[n])
                print("Min flux = ", np.min(flux))
                print("Max flux = ", np.max(flux))
                print("Median abs flux = ", np.median(np.abs(flux)))

                # 1) Plot 2D flux maps:
                cmin = -1e6
                cmax = 1e6
                im = ax.pcolormesh(x, y, flux, shading="auto",  
                # norm=mpl.colors.SymLogNorm(linthresh=lts[n], linscale=0.0001,vmin=lims[n][0], vmax=lims[n][1]),
                norm=mpl.colors.SymLogNorm(linthresh=1e-12, linscale=1e-15,vmin=cmin, vmax=cmax),
                cmap="seismic", zorder=1)
                ax.scatter([xp], [yp], color='yellow', marker='.', edgecolors='black')
                ax.set_aspect("equal")

                # Plot ring edges:
                # inner_edge = radii[int(inner_edge_i[s,n])]
                # outer_edge = radii[int(outer_edge_i[s,n])]
                # inner = plt.Circle((0,0), inner_edge, color='k', fill=False, linestyle="dashed")
                # outer = plt.Circle((0,0), outer_edge, color='k', fill=False, linestyle="dashed")
                # ax.add_patch(inner)
                # ax.add_patch(outer)
                # ax.set_aspect("equal")
                # ax.set_xlim(-2,2)
                # ax.set_ylim(-2,2)

                if s == 0:
                    ax.set_title(f"St={round(stokes[n],4)}")
                if s != len(sims)-1:
                    ax.set_xticks([])
                if n != 0:
                    ax.set_yticks([])
                if s == len(sims) - 1:
                    ticks = [cmin, 0, cmax]
                    cbar = fig_f.colorbar(im, ax=ax, orientation="horizontal", ticks=ticks)
                    cbar.set_label('$\\Sigma_{d} v_{r} $')
                
                # 2) Plot 1D flux within ring:
                colour = colour_cycler_st[n]
                ring_flux = flux[int(inner_edge_i[s,n]):int(outer_edge_i[s,n])]
                ravg_flux = np.mean(ring_flux,axis=0)    # avg over all radii within ring
                ax_f1d[s].plot(phis, ravg_flux, color=colour, label=f"St={round(stokes[n],4)}")
    
            ax_f1d[s].set_title(f"{mp} $M_\oplus$")
            ax_f1d[s].set_xlabel("$\phi$")
            ax_f1d[s].set_xlim(-np.pi, np.pi)
            ax_f1d[s].set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
            ax_f1d[s].set_xticklabels(["$-\pi$", "$-\pi/2$", "0", "$\pi/2$", "$\pi$"])
            phi_planet = np.arctan(yp/xp)
            ax_f1d[s].axvline(phi_planet, linestyle="dashed", color="k")
        
        ax_f1d[0].set_ylabel("$\Sigma_{d}v_{r}$")
        ax_f1d[0].legend()
        # fig_f.suptitle(f"H/R = {hr0}")
        fig_f1d.suptitle(f"H/R = {hr0}")
        fig_f.tight_layout()
        fig_f1d.tight_layout()
        fig_f.savefig(f"{plots_savedir}/flux_{hr0}.png", dpi=200)
        fig_f1d.savefig(f"{plots_savedir}/ring_flux_{hr0}.png", dpi=200)
    
    if "rpeak" in plots:
        # Plot ring peak vs Hill radius
        print("Plotting ring peak against Hill radius....")
        hill_radii = (planet_masses/(3*333030))**(1/3)
        for i,st in enumerate(stokes):
            colour = colour_cycler_st[i]
            # alpha_st = round(alpha/st, 4)
            ax_peak.scatter(hill_radii, r_peaks[:,i], color=colour, marker="x")
            ax_peak.plot(hill_radii, r_peaks[:,i], color=colour, label=f"St = {round(st,3)}" )
            ax_peak.scatter(hill_radii, inner_edges[:,i]+r_peaks[:,i], color=colour)
            ax_peak.plot(hill_radii, inner_edges[:,i]+r_peaks[:,i], color=colour, linestyle="dashed")
            # M_iso = calculate_Miso(hr0, alpha, st, scaling="B18")
        
        rh_range = np.linspace(np.min(hill_radii)*0.8, np.max(hill_radii)*1.2)
        lodato_gap = 1 + (5.5 * rh_range)
        ax_peak.plot(rh_range, lodato_gap, color=plain_clr, label="Dust gap edge \n (Lodato et al. 2019)")

        ax_peak.set_xlabel("Hill Radius ($r_{p}$)")
        ax_peak.set_ylabel("Location ($r_{p}$)")
        ax_peak.set_title(f"$h_{{p}}/r_{{p}}$ = {hr0}", pad=10)
        np.savetxt(f"{hr0}_peaks.txt", np.array([hill_radii, r_peaks[:,0], r_peaks[:,1], r_peaks[:,2], r_peaks[:,3], r_peaks[:,4], inner_edges[:,0]+r_peaks[:,0],  inner_edges[:,1]+r_peaks[:,1],  inner_edges[:,2]+r_peaks[:,2],  inner_edges[:,3]+r_peaks[:,3],  inner_edges[:,4]+r_peaks[:,4]]))

        ax_peak.set_ylim(np.min(inner_edges+r_peaks)*0.9, np.max(r_peaks)*1.2)
        ax_peak.set_xlim(np.min(hill_radii)-0.002, np.max(hill_radii)+0.0023)
        # for h/r=0.05
#        ax_peak.set_xticks([0.02, 0.025, 0.03, 0.035])
#        ax_peak.set_ylim(1.05,1.4)  
            
        # for h/r=0.06
        # ax_peak.set_ylim(1.07,1.45)

        # for h/r=0.07
        # ax_peak.set_ylim(1.05,1.6)
        # ax_peak.set_xlim(np.min(hill_radii)*0.95, np.max(hill_radii)*1.05)
        # ax_peak.set_ylim(np.min(inner_edges+r_peaks)*0.98, np.max(r_peaks)*1.06)
        # ax_peak.set_title(f"H/R = {hr0}")
        ax_peak.legend(loc="upper left", ncol=2)
        fig_peak.tight_layout()
        fig_peak.savefig(f"{plots_savedir}/ring_peaks_{hr0}.png", dpi=200)
    
    if "vdrift" in plots:
        v_drifts = np.array(v_drifts)
        # ax_v = ax_v.flatten()
        for s, sim in enumerate(sims):
            colour = colour_cycler_m[s]
            mp = re.search(r"(\d+)Me", sim).group(1)
            ax_v.plot(radii, v_drifts[s,0], color=colour, label=f"{mp}$M_\oplus$")
            ax_v.axvline(r_peaks[s,0], linestyle="dashed", color=colour)
        ax_v.set_xlim(0.5,1.5)
        ax_v.set_ylim(np.min(v_drifts[s,0])*0.9,np.max(v_drifts[s,0])*1.1)
        ax_v.set_xlabel("r")
        ax_v.set_ylabel("$v_{drift}$")
        ax_v.fill_between(x=radii, y1=0, y2=-20, color='lightgrey',  interpolate=True, alpha=.5)
        # ax_v[-1].remove()
        ax_v.legend()
        fig_v.tight_layout()
        fig_v.savefig(f"{plots_savedir}/vdrift_{hr0}.png")

    if "ppf" in plots:
        for s, sim in enumerate(sims):
            colour = colour_cycler_m[s]
            mp = re.search(r"(\d+)Me", sim).group(1)
            ax_pf.plot(radii, ppfs[s], color=colour, label=f"{mp}$M_\oplus$")
        
        ax_pf.set_xlabel("Radius ($r_{p}$)")
        ax_pf.set_ylabel("$\mathcal{P}_{pf}$")
        ax_pf.set_xlim(1,1.4)
        ax_pf.set_ylim(0,1.1)
        ax_pf.legend()
        fig_pf.savefig(f"{plots_savedir}/ppf_{hr0}.png")

    if "redge" in plots:
        # Plot ring edge locations vs planet mass
        print("Plotting ring edge locations against planet mass....")
        pms = np.linspace(np.min(planet_masses), np.max(planet_masses), 100)
        hill_radii = (pms/(3*333030))**(1/3)
        for i,st in enumerate(stokes):
            colour = colour_cycler_st[i]
            ax_re.scatter(planet_masses, inner_edges[:,i], color=colour)
            ax_re.plot(planet_masses, inner_edges[:,i], color=colour, label=f"St = {st}$" )
            ax_re.scatter(planet_masses, outer_edges[:,i], color=colour)
            ax_re.plot(planet_masses, outer_edges[:,i], color=colour)
        # ax_re.plot(pms, 1+3.5*hill_radii, color=plain_clr, label="$R_{p} + 3.5R_{Hill}$", linestyle="dotted")
        # ax_re.plot(pms, 1.4-10*(pms/333030)**0.5, color="red", linestyle="dotted")

        # ax_re.set_ylim(-0.25,0.25)
        ax_re.set_ylim(-0.15,0.2)
        ax_re.set_xlim(np.min(planet_masses)-2,np.max(planet_masses)+5)
        ax_re.set_xlabel("Planet mass (M$_\oplus$)")
        ax_re.set_ylabel("Ring edge locations")
        ax_re.set_title(f"H/R = {hr0}")
        ax_re.fill_between(x=np.arange(0,150), y1=0, y2=-20, color='lightgrey',  interpolate=True, alpha=.3)

        ax_re.legend(loc="upper right")
        fig_re.tight_layout()
        fig_re.savefig(f"{plots_savedir}/ring_edges_{hr0}.png", dpi=200)

    if plot_window:
        plt.show()

