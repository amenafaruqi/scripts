import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter, NullFormatter
import argparse
import re
plt.style.use('default')


regime_bounds = np.array([
    [0, 0, 10.3, 150],        #  10 Mearth
    [0, 0, 26.0, 150],        #  20 Mearth
    [0, 4.6, 150, 150],       #  40 Mearth
    [0, 14.7, 150, 150],      #  80 Mearth
    # [0, 21.2, 150, 150],    # 100 Mearth
    [0, 150, 150, 150],      # 160 Mearth
    ])

# ====================== Gas Sigma ========================

def overlay_gas_sigmas(fig, ax, radii, sigma_gas_1D, model_num=0):
    print("Plotting gas surface density....")
    if "mig" in sim:
        color = "royalblue"
        label = "Migrating"
    else:
        color = "k"
        label = "Stationary"

    planetmass = re.search(r"Mp(\d+)_", sims[model_num]).group(1)    # get planet mass from file path

    # Indicate regime boundaries on plot
    green_col = plt.get_cmap('summer')(np.linspace(0, 1, 3))[1]
    ax[int(s/2)].fill_betweenx(np.arange(0,250), regime_bounds[int(s/2),2], regime_bounds[int(s/2),3], color="orangered", alpha=0.5)
    ax[int(s/2)].fill_betweenx(np.arange(0,250), regime_bounds[int(s/2),1], regime_bounds[int(s/2),2], color="yellow", alpha=0.4)
    ax[int(s/2)].fill_betweenx(np.arange(0,250), regime_bounds[int(s/2),0], regime_bounds[int(s/2),1], color=green_col, alpha=0.8)

    ax[int(s/2)].plot(radii, sigma_gas_1D, color=color, label=label)
    ax[int(s/2)].set_title(f"{planetmass} $M_\oplus$")

    if planets:
        if "stat" in sim:
            planetcolour = 'k'
        else:
            planetcolour = color
        for rp in rps:
            ax[int(s/2)].axvline(rp, linestyle='dashed', color=planetcolour)
        
    ax[int(s/2)].set_xlabel("Radius (AU)")
    ax[0].set_ylabel("$\Sigma_{gas} (g/cm^{2})$")
    ax[int(s/2)].set_xscale("log")
    ax[int(s/2)].set_yscale("log")
    # ax[int(s/2)].set_ylim(1e-5, 1e2)
    ax[int(s/2)].set_xlim(np.min(radii), np.max(radii))
    ax[int(s/2)].set_ylim(0.1,200)
    ax[0].legend()
    fig.tight_layout()


# ================= Plot dust-gas ratios =====================

def overlay_dust_gas_ratio(fig, ax, radii, dust_mass_tot, gas_mass, model_num=0):
    epsilon = dust_mass_tot/gas_mass
    if "mig" in sim:
        color = "royalblue"
        label = "Migrating"
    else:
        color = "k"
        label = "Stationary"

    planetmass = re.search(r"Mp(\d+)_", sims[model_num]).group(1)    # get planet mass from file path

    # Indicate regime boundaries on plot
    green_col = plt.get_cmap('summer')(np.linspace(0, 1, 3))[1]
    ax[int(s/2)].fill_betweenx(np.arange(0,250), regime_bounds[int(s/2),2], regime_bounds[int(s/2),3], color="orangered", alpha=0.5)
    ax[int(s/2)].fill_betweenx(np.arange(0,250), regime_bounds[int(s/2),1], regime_bounds[int(s/2),2], color="yellow", alpha=0.4)
    ax[int(s/2)].fill_betweenx(np.arange(0,250), regime_bounds[int(s/2),0], regime_bounds[int(s/2),1], color=green_col, alpha=0.8)

    ax[int(s/2)].plot(radii, epsilon, color=color, label=label)
    ax[int(s/2)].set_title(f"{planetmass} $M_\oplus$")

    if planets:
        if "stat" in sim:
            planetcolour = 'k'
        else:
            planetcolour = color
        for rp in rps:
            ax[int(s/2)].axvline(rp, linestyle='dashed', color=planetcolour)
    
    ax[int(s/2)].set_xlabel("Radius (AU)")
    ax[0].set_ylabel("Dust-Gas Ratio")
    ax[int(s/2)].set_xscale("log")
    ax[int(s/2)].set_yscale("log")
    ax[int(s/2)].set_ylim(1e-5, 1e2)
    ax[int(s/2)].set_xlim(np.min(radii), np.max(radii))
    ax[0].legend()
    fig.tight_layout()

# ================= Plot growth timescales =====================

def plot_tgrowth(fig, ax, radii, a, sigma_dust_1D, model_num=0):
    R, A = np.meshgrid(radii, a)
    levels = np.linspace(-4, 20, 7)                   
    print("Plotting growth timescales....")
    ax0 = ax.flatten()[model_num]

    hr = hr0*(radii**0.25)
    h = hr*radii
    St = A*rhodust*np.pi/(sigma_gas_1D*2)
    cs = h*2*np.pi/(radii**0.5)
    tau_growth = rhodust*A*h/(sigma_dust_1D*cs*(alpha*St*3))
    tau_growth = tau_growth*1e-6   # convert to Myr

    con = ax0.contourf(R,A, np.log10(tau_growth), cmap="YlGnBu", levels=levels)
    
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.set_ylabel("a (cm)")
    ax0.set_xlabel("Radius (AU)")
    ax0.set_ylim(np.min(a), np.max(a))

    if planets:
        planetmass = re.search(r"Mp(\d+)_", sims[model_num]).group(1)    # get planet mass from file path
        for rp in rps:
            ax0.axvline(rp, linestyle='dashed', color='black')
        ax0.set_title(f"{planetmass} $M_\oplus$")
        
    if not model_num%plotsizex:
        ax0.set_ylabel("a (cm)")
    else:
        ax0.set_yticks([])
    if model_num < plotsizex and len(sims) > plotsizex:
        ax0.set_xticks([])
    else:
        ax0.set_xlabel("Radius (AU)")

    # fig3.subplots_adjust(right=0.89, hspace=0.35)
    # cbar_ax = fig3.add_axes([0.91, 0.53, 0.02, 0.4])
    # fig3.colorbar(con, cax=cbar_ax, orientation="vertical", label="log$[\\tau_{growth} (Myr)]$")

    fig.tight_layout()
    
    if model_num == len(sims)-1:
        ax_cbar = ax.flatten()[model_num+1]
        ax_cbar.remove()
        fig.subplots_adjust(right=0.89, hspace=0.3)
        cax = fig.add_axes([ax0.get_position().x1+0.09,ax0.get_position().y0,0.02,ax0.get_position().height])
        fig.colorbar(con, cax=cax, orientation="vertical", label="$\\log[\\tau_{growth} (Myr)]$")
    fig.tight_layout()



# =================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate 1D plots', prefix_chars='-')
    
    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/astro/phrkvg/simulations/new_lowres_models/Mdisc0.015/"],help="working directory containing simulations")
    parser.add_argument('-sims', metavar='sims', type=str, nargs="*", default=[] ,help="simulation directories containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images/comparison_plots" ,help="directory to save plots to")
    parser.add_argument('-o', metavar='outputs',default=[30], type=int, nargs=1 ,help="output to plot")
    parser.add_argument('-plots', metavar='plots',default=[], type=str, nargs="*" ,help="plots to produce")
    parser.add_argument('-noplanet', action="store_false")
    parser.add_argument('-nogrog', action="store_false")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-style', metavar='style', type=str, nargs="*", default=["publication"], help="style sheet to apply to plots")

    args = parser.parse_args()
    o = args.o[0]
    plots = args.plots
    wd = args.wd[0]
    sims = args.sims
    planets = args.noplanet
    grog = args.nogrog
    plot_window = args.plot_window
    plots_savedir = args.savedir
    style = args.style

    cm = plt.get_cmap('viridis')
    colour_cycler = cm(np.linspace(0, 1, 5))
    # colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    # linestyles = ['solid', 'dashdot', 'dotted']

    if plot_window:
        matplotlib.use('TkAgg')
    
    if not plot_window:
        for s in style:
            plt.style.use([f"../styles/{s}.mplstyle"])

    # =================== Define figures and axes ========================
    if "gsig" in plots:
        fig_gas_sigma, ax_gas_sigma = plt.subplots(figsize=(15,4), ncols=int(len(sims)/2), sharey=True)
    if "dgr" in plots:
        fig_dgr, ax_dgr = plt.subplots(figsize=(15,4), ncols=int(len(sims)/2), sharey=True)
    if "growth" in plots:   # only specify dcon for grog models
        fig_growth, ax_growth = plt.subplots(figsize=(12,8), nrows=int(len(sims)/2), ncols=2)


    # ================== Read in data at timesteps =======================

    for s, sim in enumerate(sims):
        print (f"============ Plotting data for {sim} =============")
        params_file = f'{wd+sim}/variables.par'
        params_dict = {}

        plotsizex = 3
        plotsizey = 2

        if ("10_" in sim) or ("20_" in sim):
            output = int(o/5)   # need to account for different timestepping in different models
        elif ("100_" in sim):
            outut = int(o*2)
        else:
            output = o

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
        if grog:
            mingsize = float(params_dict['MIN_GRAIN_SIZE'])
            maxgsize = float(params_dict['MAX_GRAIN_SIZE'])
            nss_coag = int(params_dict['NUMSUBSTEPS_COAG'])
            rhodust = float(params_dict['RHO_DUST'])
            densfloor = float(params_dict['DENSITY_FLOOR'])
        else:
            max_stokes = float(params_dict['STOKES'])

        dt_orbits = int(float(params_dict['DT'])/(2*np.pi))   # 2pi = 1 orbit = 1 yr
        ninterm = float(params_dict['NINTERM'])               # number of dts between outputs
        dt_outputs = dt_orbits*ninterm                        # time between outputs
        t = o*dt_outputs*1e-6                                 # time in Myr

        # ----------------------------------------------------

        if planets:
            with open(f"{wd+sim}/planet.cfg", 'r') as pfile:
                num_planets = len(pfile.readlines()) - 5        # ignore 5 lines of headers
            
            rps = np.zeros((num_planets))
            for n in range(num_planets):
                planet_data = np.unique(np.loadtxt(f"{wd+sim}/planet{n}.dat"), axis=0)
                xp, yp = planet_data[output,1], planet_data[output,2]
                # xp, yp = planet_data[np.array(outputs)*5][:,1], planet_data[np.array(outputs)*5][:,2]
                rps[n] = ((xp**2) + (yp**2))**0.5

            planet0_period = (rps[0]**3)**0.5                 # orbital period of planet in yrs
            planet_orbits = t*1e6/planet0_period

        # FARGO initialises grains with sizes uniformly distributed across ndust bins in logspace
        if grog:
            a = np.logspace(np.log10(mingsize),    
                            np.log10(maxgsize),
                            ndust+1)

            a = (0.5*(a[1:] + a[:-1]))                                   # grain sizes in middles of bins (in cm)
        else:
            stokes = np.logspace(np.log10(max_stokes),np.log10(max_stokes*10**(-2)),ndust)

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
        sigma_gas = np.zeros((nrad, nphi))
        gas_mass = np.zeros((nrad))
        sigma_dust = np.zeros((ndust, nrad, nphi))
        dust_mass = np.zeros((ndust, nrad))

        gasfile = f"/gasdens{output}.dat" 
        sigma_gas = np.fromfile(wd+sim+gasfile).reshape(nrad,nphi)/(1.125e-7)         # convert back to g/cm2
        n_grains = np.zeros((ndust))
        for n in np.arange(ndust):
            dust_file = f"/dustdens{n}_{output}.dat"
            sigma_dust[n] = np.fromfile(wd+sim+dust_file).reshape(nrad,nphi)/(1.125e-7)

        sigma_dust_azimsum = np.sum(sigma_dust, axis=2)                  # sum over all phi 
        sigma_gas_azimsum = np.sum(sigma_gas, axis=1)                    # sum over all phi   

        sigma_gas_1D = sigma_gas_azimsum/nphi                           # dimensions: (nrad) 
        sigma_dust_1D = sigma_dust_azimsum/nphi                         # dimensions: (ndust, nrad)  

        # dust mass for dust of size a as a function of r
        dust_mass = [2*np.pi*radii*sigma_dust_1D[n,:]*delta_r*333030*(1.125e-7) for n in range(ndust)]
        gas_mass = 2*np.pi*radii*sigma_gas_1D*delta_r*333030*(1.125e-7)      # convert from Msun to Mearth
        dust_mass_tot = np.sum(dust_mass, axis=0)

        if grog:
            uf = 0.0021                                           # fragmentation velocity 10 m/s in AU/yr
            hr = hr0*(radii**f)                                   # aspect ratio
            b = (uf**2)*radii/(4*(np.pi**2)*alpha*(hr**2))

            cs = hr*(((2e30)*(6.67e-11))/(radii*1.5e11))**0.5     # [m/s]
            p = (sigma_gas_1D*(cs**2)/((2*np.pi)**0.5))*(hr**-1)*((radii*1.5e11)**-1)
            pad = np.empty((1, 1))*np.nan
            gamma = (radii/p)*np.abs(np.append(np.diff(p)/np.diff(radii), pad))
            C = 2/(np.pi*hr)
            a_St1 = (2/np.pi)*(sigma_gas_1D/rhodust)              # plot St=1 line

            # size of largest grains in a fragmentation-dominated distribution
            a_frag = 2*sigma_gas_1D*b/(rhodust*3*np.pi)

            # size of largest grains in a drift-dominated distribution
            a_drift = (2/np.pi)*(np.sum(sigma_dust_1D, axis=0)/(rhodust*gamma*hr**2))    
    
        if "gsig" in plots:
            overlay_gas_sigmas(fig_gas_sigma, ax_gas_sigma, radii, sigma_gas_1D, s)
        if "dgr" in plots:
            overlay_dust_gas_ratio(fig_dgr, ax_dgr, radii, dust_mass_tot, gas_mass, s)
        if "growth" in plots:
            plot_tgrowth(fig_growth, ax_growth, radii, a, sigma_dust_1D, s)
    
    # ======================== Generate Plots ==========================
    print(f"-------------------\nPlotting comparison plots for {sims}\n=============")

    if "gsig" in plots:
        fig_gas_sigma.savefig(f"{plots_savedir}/fullcomp_sigmagas.png")
    if "dgr" in plots:
        fig_dgr.savefig(f"{plots_savedir}/fullcomp_dgr.png")
    if "growth" in plots:
        fig_growth.savefig(f"{plots_savedir}/fullcomp_growth.png")

    if plot_window:
        plt.show()
