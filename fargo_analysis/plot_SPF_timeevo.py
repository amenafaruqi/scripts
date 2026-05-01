import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import re

plt.style.use('default')

import warnings
warnings.filterwarnings("ignore")


def calculate_ring_mass(radii, dust_mass, rps, times):
    # Set labels and indexes to populate arrays
    if "stat" in sim:
        if len(rps[0,:]) == 1:
            rowi = 0
            coli = 0
        else:
            rowi = 1
            coli = 0
    else:
        if len(rps[0,:]) == 1:
            rowi = 0
            coli = 1
        else:
            rowi = 1
            coli = 1

    n_size_decades = int(np.log10(maxgsize) - np.log10(mingsize))
    size_decades = np.split(np.arange(ndust), n_size_decades)

    # dust_mass has dimensions of n_outputs x ndust x nrad
    ring_masses = np.zeros((len(times),2,2,n_size_decades))  # timesteps x 2 x 2 x 7 size decades
    
    planetmass = int(re.findall(r"Mp(\d+)_", sim)[0])          # Take inner planet mass only
    for ti, t in enumerate(times):
        r_hill = rps[ti,0]*((planetmass*1e-6)**(1/3))

        # Find radial cells closest to 2 and 20 R_Hill
        i_inner = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[ti,0])+r_hill*2)))
        i_outer = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[ti,0])+r_hill*20)))
        if rowi:     # i.e. if 2-planet model, bound by outer planet location
            rp_p2 = rps[ti,1]      # Outer planet location
            i_p2 = min(range(len(radii)), key=lambda i: abs(radii[i]-((rp_p2))))
            # Use outer planet location as outer bound if it lies within 20 R_Hill (of inner planet)
            i_outer = np.min([i_outer, i_p2])

        dust_mass_ring_t = np.sum(dust_mass[ti,:,i_inner:i_outer], axis=1)       # dimensions = (ndust)

        # Sum within each size decade to go from 70 to 7 ring masses
        for n, size_decade in enumerate(size_decades):
            ring_mass_by_size = np.sum(dust_mass_ring_t[size_decade])      # sum over all grain sizes within size decade to get single value
            ring_masses[ti,rowi, coli, n] = ring_mass_by_size    # ring mass per size decade (for 1 model and timestep)


    # print(ring_masses) 
    return ring_masses


def plot_ring_mass(fig, ax, ring_masses, times):
    # Set labels and indexes to populate arrays
    if "stat" in sim:
        if len(rps[0,:]) == 1:
            label = "S1"
            rowi = 0
            coli = 0
        else:
            label = "S2"
            rowi = 1
            coli = 0
    else:
        if len(rps[0,:]) == 1:
            label = "M1"
            rowi = 0
            coli = 1
        else:
            label = "M2"
            rowi = 1
            coli = 1

    size_labels = ["$10^{-5}-10^{-4}$ cm", "$10^{-4}-10^{-3}$ cm", "$10^{-3}-10^{-2}$ cm", "$10^{-2}-10^{-1}$ cm", "$10^{-1}$-1 cm", "1-10 cm", "10-100 cm"]

    ax0 = ax[rowi, coli]
    ax0.set_prop_cycle(color=colour_cycler)
    # planetmass = int(re.findall(r"Mp(\d+)_", sim)[0])

    for i in np.arange(7):
        color = next(ax0._get_lines.prop_cycler)['color']
        ax0.plot(times, ring_masses[:,rowi, coli,i], color=color, label=size_labels[i])
        ax0.scatter(times, ring_masses[:,rowi, coli,i], color=color)

    ax0.plot(times, np.sum(ring_masses, axis=3)[:,rowi, coli], color="k")
    if (rowi == 0) and (coli == 0):
        ax0.legend()
    
    ax0.set_title(label)

    # if sim == sims[-1]:        
    #     cmap = mpl.cm.viridis_r
    #     bounds = np.arange(-5,3)
    #     norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    #     box = ax.get_position()
    #     ax.set_position([box.x0, box.y0, box.width, box.height*0.9])
    #     cbar_ax = fig.add_axes([box.x0, 0.83, box.width, box.height*0.05])

    #     cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
    #                 cax=cbar_ax, orientation='horizontal',
    #                 label="Dust size (cm)", ticklocation="top")
    #     cbar.set_ticklabels(["$10^{-5}$", "$10^{-4}$", "$10^{-3}$", "$10^{-2}$", "$10^{-1}$", "$10^{0}$", "$10^{1}$", "$10^{2}$"])

    if coli == 0:
        ax0.set_ylabel("Ring dust mass ($M_\oplus$)")
    if rowi == 1:
        ax0.set_xlabel("Time (Myr)")

    fig.tight_layout()



# =================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate 1D plots', prefix_chars='-')
    
    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/astro/phrkvg/simulations/"],help="working directory containing simulations")
    parser.add_argument('-sims', metavar='sims', type=str, nargs="*", default=[] ,help="simulation directories containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images/SPF_plots" ,help="directory to save plots to")
    parser.add_argument('-o', metavar='outputs',default=[10], type=int, nargs="*" ,help="output to plot")
    parser.add_argument('-plots', metavar='plots',default=[], type=str, nargs="*" ,help="plots to produce")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-style', metavar='style', type=str, nargs="*", default=["publication"], help="style sheet to apply to plots")

    args = parser.parse_args()
    o = args.o[-1]                    # Last timestep given
    tevo_o = args.o                   # List of timesteps for time evo plots 
    plots = args.plots
    wd = args.wd[0]
    sims = args.sims
    plot_window = args.plot_window
    plots_savedir = args.savedir
    style = args.style

    cm = plt.get_cmap('viridis_r')
    colour_cycler = cm(np.linspace(0, 1, 7))       # For plotting different mass decades

    if plot_window:
        mpl.use('TkAgg')
    else:
        for s in style:
            plt.style.use([f"../styles/{s}.mplstyle"])

    # =================== Define figures and axes ========================

    if "rmass" in plots:
        fig_m, ax_m = plt.subplots(ncols=2, nrows=2, figsize=(9,8))

    # ================== Read in data at timesteps =======================

    # Array to store ring masses
    # ring_masses = np.zeros((len(tevo_o),2,2,7))           # 4 models x 7 size decades for each timestep
    mps = []

    for s, sim in enumerate(sims):
        print (f"============ Plotting data for {sim} =============")
        params_file = f'{wd+sim}/variables.par'
        params_dict = {}

        planet_masses = re.findall(r"Mp(\d+)_", sim)    # Planet masses, stored as strings in list

        # Load model params into dict
        param_lines = open(params_file).readlines()
        for line in param_lines:
            if line.split():
                param_label, param_value = line.split()[0:2]
                params_dict.update([(param_label, param_value)])

        # Read params from .var file
        nphi = int(params_dict['NX'])
        nrad = int(params_dict['NY'])
        f = float(params_dict['FLARINGINDEX'])
        hr0 = float(params_dict['ASPECTRATIO'])      # aspect ratio at R=1AU
        ndust = int(params_dict['NDUST'])
        alpha = float(params_dict['ALPHA'])
        spacing = str(params_dict['SPACING'])
        mingsize = float(params_dict['MIN_GRAIN_SIZE'])
        maxgsize = float(params_dict['MAX_GRAIN_SIZE'])
        nss_coag = int(params_dict['NUMSUBSTEPS_COAG'])
        rhodust = float(params_dict['RHO_DUST'])
        densfloor = float(params_dict['DENSITY_FLOOR'])
        dt_orbits = int(float(params_dict['DT'])/(2*np.pi))   # 2pi = 1 orbit = 1 yr
        ninterm = float(params_dict['NINTERM'])               # number of dts between outputs
        dt_outputs = dt_orbits*ninterm                        # time between outputs
        t = np.array(tevo_o)*dt_outputs*1e-6                  # time in Myr

        # ----------------------------------------------------

        with open(f"{wd+sim}/planet.cfg", 'r') as pfile:
            num_planets = len(pfile.readlines()) - 5        # ignore 5 lines of headers
        
        rps = np.zeros((len(tevo_o),num_planets))
        xps = np.zeros((len(tevo_o),num_planets))
        yps = np.zeros((len(tevo_o),num_planets))
        for n in range(num_planets):
            planet_data = np.unique(np.loadtxt(f"{wd+sim}/planet{n}.dat"), axis=0)
            xp, yp = planet_data[tevo_o,1], planet_data[tevo_o,2]
            xps[:,n] = xp
            yps[:,n] = yp
            rps[:,n] = ((xp**2) + (yp**2))**0.5


        # FARGO initialises grains with sizes uniformly distributed across ndust bins in logspace
        a = np.logspace(np.log10(mingsize),    
                        np.log10(maxgsize),
                        ndust+1)

        a = (0.5*(a[1:] + a[:-1]))                                   # grain sizes in middles of bins (in cm)

        r_cells = np.loadtxt(f'{wd+sim}/domain_y.dat')[3:-3]             # ignore ghost cells
        phi_cells = np.loadtxt(f'{wd+sim}/domain_x.dat')
        if spacing == "Linear":
            radii = np.array([(r_cells[n]+r_cells[n+1])/2 for n in range(len(r_cells)-1)])
            delta_r = radii[1]-radii[0]
        else:     # Log grid
            radii = np.array([np.exp((np.log(r_cells[n])+np.log(r_cells[n+1]))/2) for n in range(len(r_cells)-1)])
            delta_log_r = np.log(radii[1]) - np.log(radii[0])
            delta_r = radii*delta_log_r
            
        phis = np.array([(phi_cells[n]+phi_cells[n+1])/2 for n in range(len(phi_cells)-1)])

        # Get gas and dust sigma
        # sigma_gas = np.zeros((nrad, nphi))
        # gas_mass = np.zeros((nrad))
        sigma_dust = np.zeros((len(tevo_o), ndust, nrad, nphi))
        dust_mass = np.zeros((len(tevo_o), ndust, nrad))

        # gasfile = f"/gasdens{output}.dat" 
        # sigma_gas = np.fromfile(wd+sim+gasfile).reshape(nrad,nphi)/(1.125e-7)         # convert back to g/cm2
        n_grains = np.zeros((ndust))
        for oi, o in enumerate(tevo_o):
            for n in np.arange(ndust):
                dust_file = f"/dustdens{n}_{o}.dat"
                sigma_dust[oi, n] = np.fromfile(wd+sim+dust_file).reshape(nrad,nphi)/(1.125e-7)

        sigma_dust_azimsum = np.sum(sigma_dust, axis=3)                  # sum over all phi 
        # sigma_gas_azimsum = np.sum(sigma_gas, axis=1)                    # sum over all phi   

        # sigma_gas_1D = sigma_gas_azimsum/nphi                           # dimensions: (nrad) 
        sigma_dust_1D = sigma_dust_azimsum/nphi                         # dimensions: (n_outputs, ndust, nrad)  

        # dust mass for dust of size a as a function of r
        dust_mass = 2*np.pi*radii*sigma_dust_1D*delta_r*333030*(1.125e-7)
        # gas_mass = 2*np.pi*radii*sigma_gas_1D*delta_r*333030*(1.125e-7)      # convert from Msun to Mearth
        dust_mass_tot = np.sum(dust_mass, axis=1)                            # summed over all sizes

        if "rmass" in plots:
            ring_masses = calculate_ring_mass(radii, dust_mass, rps, t)       # assume only 1 planet here
            plot_ring_mass(fig_m, ax_m, ring_masses, t)


    # ======================== Generate Plots ==========================
    print(f"-------------------\nPlotting comparison plots for {sims}\n=============")

    if "rmass" in plots:
        fig_m.savefig(f"{plots_savedir}/SPF_ring_masses.png")


    if plot_window:
        plt.show()
