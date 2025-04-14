import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import pandas as pd
import math
from scipy.optimize import curve_fit
from scipy.signal import peak_widths
from scipy.signal import find_peaks
from numpy.polynomial import Polynomial
import glob, os

# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
# plt.rcParams.update({'font.size': 36})
# plt.style.use('./publication.mplstyle')

Nr = 504
Nphi = 1170

species = [0,1,2]   # make sure to input selected stokes
hs =  [0.05]
alphas = [1e-4]
st = [0.2, 0.02, 0.002]

# sims = [
#     "10Me_B0_h5_a1e-4", "10Me_B1e-1_h5_a1e-4", "10Me_B1_h5_a1e-4", "10Me_B1e1_h5_a1e-4",
#     "20Me_B0_h5_a1e-4", "20Me_B1e-1_h5_a1e-4", "20Me_B1_h5_a1e-4", "20Me_B1e1_h5_a1e-4",
#     "40Me_B0_h5_a1e-4", "40Me_B1e-1_h5_a1e-4", "40Me_B1_h5_a1e-4", "40Me_B1e1_h5_a1e-4",
#     "80Me_B0_h5_a1e-4", "80Me_B1e-1_h5_a1e-4", "80Me_B1_h5_a1e-4", "80Me_B1e1_h5_a1e-4",
#     ]

sims = ["10Me", "20Me", "40Me"]


# simdir = "/home/amena/projects/rrg-pudritz-ab/pudritz_bogert/bogerg1/FARGO3D_new/outputs/fargo_multifluid/"
simdir = "/home/astro/phrkvg/simulations/dusty_fargo_models/"

fig,axes = plt.subplots(figsize=(28, 9*len(sims)),ncols=len(st)+1, nrows=len(sims))

cols = ["Gas", "St=0.2", "St=0.02", "St=0.002"]

# rows = [s.split("Me")[0] + "$M_\oplus$" for s in sims]
# for ri in range(1,len(rows),2):
#     rows[ri] = rows[ri] + " iso"
#     rows[ri-1] = rows[ri-1] + " $beta=10^{-3}$"

# rows = ["iso", "$beta=10^{-3}$", "$beta=10^{-2}$", "$beta=10^{-1}$", "beta=1"]

for ax, col in zip(axes[0], cols):
    ax.set_title(col)
for ax, sim in zip(axes[:,0], sims):
    ax.set_ylabel(sim, rotation=90, size='large')

# fig.tight_layout()
# Extract data (both data_0 and last output, either data_75 or data_30)
for simi,sim in enumerate(sims):
    last_output = 600
    rho = np.fromfile(f"{simdir+sim}/gasdens{last_output}.dat").reshape(Nr,Nphi)

    phi = np.linspace(0,2*np.pi,Nphi)
    r_cells = np.loadtxt(f'{simdir+sim}/domain_y.dat')[3:-3]             # ignore ghost cells
    r = np.array([(r_cells[n]*r_cells[n+1])**0.5 for n in range(len(r_cells)-1)])
    R, Phi = np.meshgrid(r,phi)
    x = R*np.cos(Phi)
    y = R*np.sin(Phi)

    im_g = axes[simi,0].pcolormesh(x, y, np.log10(rho).T, shading='auto', vmin=-0.3, vmax=0.4)
    axes[simi,0].set_aspect('equal')

    planet_data = np.unique(np.loadtxt(f"{simdir+sim}/planet0.dat"), axis=0)
    xp, yp = planet_data[last_output][1], planet_data[last_output][2]
    axes[simi,0].scatter([-xp], [-yp], color='white', marker='.')

    ims = []
    lims = [[-4,-0.5],[-5,-1],[-3,-1.5]]
    for i,s in enumerate(species):
        rho = np.fromfile(f"{simdir+sim}/dustdens{s}_{last_output}.dat").reshape(Nr,Nphi)
        # rho0 = np.fromfile(f"{simdir+sim}/dustdens{s}_0.dat").reshape(Nr,Nphi)

        im = axes[simi,i+1].pcolormesh(x, y, np.log10(rho).T, shading='auto', vmin=lims[i][0], vmax=lims[i][1]) 
        axes[simi,i].set_aspect('equal')
        ims.append(im)
        axes[simi,i].scatter([-xp], [-yp], color='white', marker='.')

        # cbar = fig.colorbar(im, ax=axes[simi,i], orientation='horizontal', aspect=30)
        # cbar.set_label('$log(\\Sigma_{d}) $')

        cbar = fig.colorbar(im, ax=axes[simi,i+1], orientation='horizontal', aspect=10)
        cbar.set_label('$log(\\Sigma_{d}) $')
        # cbar2 = fig.colorbar(im, ax=axes[simi,i], orientation='horizontal', aspect=10)   
        # cbar2.set_label('$log(\\Sigma_{d}) $')
        # cbar3 = fig.colorbar(im, ax=axes[simi,i], orientation='horizontal', aspect=10)
        # cbar3.set_label('$log(\\Sigma_{d}) $')

    cbar_gas = fig.colorbar(im_g, ax=axes[simi,0], orientation='horizontal', aspect=10)
    cbar_gas.set_label('$log(\\Sigma_{g}) $')
    # cbar = fig.colorbar(im, ax=axes[-1,1:], orientation='horizontal', aspect=30)
    # cbar.set_label('$log(\\Sigma_{d}) $')

# fig.suptitle("$\\alpha=10^{-4}, \Beta = 10^{-3}, h=0.05$")
fig.tight_layout()
fig.savefig(f"./2Dsigma.png", dpi=200)

