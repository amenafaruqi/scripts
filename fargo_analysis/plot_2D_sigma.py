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




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fit dust rings', prefix_chars='-')

    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/amena/scratch/simulations/dusty_fargo/h0.05"],help="working directory containing simulations")
    parser.add_argument('-sims', metavar='sim', type=str, nargs="*", default=["10Me"] ,help="simulation directory containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="/home/amena/scratch/images/" ,help="directory to save plots to")
    parser.add_argument('-o', metavar='outputs',default=[600], type=int, nargs="*" ,help="outputs to plot")
    parser.add_argument('-style', metavar='style', type=str, nargs="*", default=["publication"] ,help="style sheet to apply to plots")

    args = parser.parse_args()
    output = args.o[0]
    wd = args.wd[0]
    sims = args.sims
    plots_savedir = args.savedir
    style = args.style

    for sty in style:
        plt.style.use([f"../styles/{sty}.mplstyle"])
        if "darkbg" in sty:
            plots_savedir = plots_savedir+"/darkbg/"
            plain_clr = "w"

    # ------------------------------------------------------
    
    # Initialise axes and arrays for data to plot
    fig, ax = plt.subplots(nrows=len(sims), ncols=nstokes, figsize=(30,25))
    print(len(sims), nstokes)

    # ------------------------------------------------------

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
        mp = planet_data[output,7]
        mp = int(round(mp/(3.0027e-6),0))     # convert to earth masses

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
        # gasfile = f"gasdens{output}.dat" 
        # sigma_gas = np.fromfile(simdir+gasfile).reshape(nrad,nphi)

        sigma_dust = np.zeros((ndust, nrad, nphi))
        for n in np.arange(ndust):
            dust_file = f"dustdens{n}_{output}.dat"
            sigma_dust[n] = np.fromfile(simdir+dust_file).reshape(nrad,nphi)
            dust_file0 = f"dustdens{n}_0.dat"

        # Plot 2D sigmas
        r,theta=np.meshgrid(radii, phis)
        x = r*np.cos(theta)
        y = r*np.sin(theta)

        vlims = [[-25,-10],[-25,-5],[-15,-2],[-10,-4],[-7,-5]]
        for st_i in np.arange(nstokes):
            sig = np.log10(sigma_dust[st_i].T)
            im = ax[s,st_i].pcolormesh(x,y,sig, shading="auto", vmin=np.min(sig)*1.3, vmax=np.max(sig)*0.9) 
            ax[s,st_i].scatter([xp], [yp], color='r', marker='.')
            # ax[s,st_i].set_xticks([])
            # ax[s,st_i].set_yticks([])
            ax[s,st_i].set_xlim(-1.5,1.5)
            ax[s,st_i].set_ylim(-1.5,1.5)
            ax[s,st_i].set_aspect("equal")
            if s == 0:
                ax[s,st_i].set_title(f"St={round(stokes[st_i],3)}")
            if st_i == 0:
                ax[s,st_i].set_ylabel(f"{mp} $M_\oplus$")

            if s == len(sims)-1:
                divider = make_axes_locatable(ax[s, st_i])
                cax = divider.append_axes("bottom", size="4%", pad=0.2)
                cb = fig.colorbar(im, cax=cax, orientation="horizontal")

    fig.savefig(f"{plots_savedir}/2dsigma_{hr0}.png", dpi=200)


    print(f"-------------------\nPlotting output {output} for {sims}\n=============")

