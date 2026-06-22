import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter, NullFormatter, LogFormatter
import argparse
import re

plt.style.use('default')

import warnings
warnings.filterwarnings("ignore")

regime_bounds = np.array([
    [0, 0, 10.3, 150],        #  10 Mearth
    [0, 0, 26.0, 150],        #  20 Mearth
    [0, 4.6, 150, 150],       #  40 Mearth
    # [0, 14.7, 150, 150],      #  80 Mearth
    [0, 21.2, 150, 150],    # 100 Mearth
    [0, 150, 150, 150],      # 160 Mearth
    ])

# ====================== Gas Sigma ========================

def overlay_gas_sigmas(fig, ax, radii, sigma_gas_1D, model_num=0):
    print("Plotting gas surface density....")
    if "mig" in sim:
        if colour_bounds:
            color = "royalblue"
        else:
            color = "red"
        label = "Migrating"
    else:
        color = "k"
        label = "Stationary"

    planetmass = re.search(r"Mp(\d+)_", sims[model_num]).group(1)    # get planet mass from file path

    # Indicate regime boundaries on plot
    if colour_bounds:
        green_col = plt.get_cmap('summer')(np.linspace(0, 1, 3))[1]
        ax[int(s/2)].fill_betweenx(np.arange(0,250), regime_bounds[int(s/2),2], regime_bounds[int(s/2),3], color="orangered", alpha=0.5)
        ax[int(s/2)].fill_betweenx(np.arange(0,250), regime_bounds[int(s/2),1], regime_bounds[int(s/2),2], color="yellow", alpha=0.4)
        ax[int(s/2)].fill_betweenx(np.arange(0,250), regime_bounds[int(s/2),0], regime_bounds[int(s/2),1], color=green_col, alpha=0.8)

    ax[int(s/2)].plot(radii, sigma_gas_1D, color=color, label=label)
    ax[int(s/2)].set_title(f"{planetmass} $M_\oplus$")

    if planets:
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
        if colour_bounds:
            color = "royalblue"
        else:
            color = "red"
        label = "Migrating"
    else:
        color = "k"
        label = "Stationary"

    planetmass = re.search(r"Mp(\d+)_", sims[model_num]).group(1)    # get planet mass from file path

    # Indicate regime boundaries on plot
    if colour_bounds:
        green_col = plt.get_cmap('summer')(np.linspace(0, 1, 3))[1]
        ax[int(s/2)].fill_betweenx(np.arange(0,250), regime_bounds[int(s/2),2], regime_bounds[int(s/2),3], color="orangered", alpha=0.5)
        ax[int(s/2)].fill_betweenx(np.arange(0,250), regime_bounds[int(s/2),1], regime_bounds[int(s/2),2], color="yellow", alpha=0.4)
        ax[int(s/2)].fill_betweenx(np.arange(0,250), regime_bounds[int(s/2),0], regime_bounds[int(s/2),1], color=green_col, alpha=0.8)

    ax[int(s/2)].plot(radii, epsilon, color=color, label=label)
    ax[int(s/2)].set_title(f"{planetmass} $M_\oplus$")


    if planets:
        planetcolour = color
        for rp in rps:
            ax[int(s/2)].axvline(rp, linestyle='dashed', color=planetcolour)
    
    ax[int(s/2)].set_xlabel("Radius (AU)")
    ax[0].set_ylabel("Dust-Gas Ratio")
    ax[int(s/2)].set_xscale("log")
    ax[int(s/2)].set_yscale("log")
    ax[int(s/2)].set_ylim(1e-5, 1)
    ax[int(s/2)].set_xlim(6, 140)
    ax[int(s/2)].axhline(1e-2, linestyle="dotted", color="k")
    ax[0].legend()
    fig.tight_layout()

# ================= Plot growth timescales =====================

def plot_tgrowth(fig, ax, radii, a, sigma_dust_1D, model_num=0):
    R, A = np.meshgrid(radii, a)
    levels = np.linspace(-4, 16, 6)                   
    print("Plotting growth timescales....")
    # ax0 = ax.flatten()[model_num]
    if "mig" in sim:
        coli = 1
        miglabel = "$\\bf{Migrating}$ \n"
    else:
        coli = 0
        miglabel = "$\\bf{Stationary}$ \n"

    ax0 = ax[int(s/2),coli]

    hr = hr0*(radii**0.25)
    h = hr*radii
    St = A*rhodust*np.pi/(sigma_gas_1D*2)
    cs = h*2*np.pi/(radii**0.5)
    tau_growth = rhodust*A*h/(sigma_dust_1D*cs*(alpha*St*3))
    tau_growth = tau_growth*1e-6   # convert to Myr

    con = ax0.contourf(R,A, np.log10(tau_growth), cmap="YlGnBu", levels=levels, extend="both")
    
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    if coli == 0:
        ax0.set_ylabel("a (cm)")
    if int(s/2) == int(len(sims)/2)-1:
        ax0.set_xlabel("Radius (AU)")
    ax0.set_ylim(np.min(a), np.max(a))

    if planets:
        planetmass = re.search(r"Mp(\d+)_", sims[model_num]).group(1)    # get planet mass from file path
        for rp in rps:
            ax0.axvline(rp, linestyle='dashed', color='black')
        if int(s/2) == 0:
            ax0.set_title(miglabel + f"{planetmass} $M_\oplus$")
        else:
            ax0.set_title(f"{planetmass} $M_\oplus$")

    if int(s) == len(sims) - 1:
        fig.subplots_adjust(bottom=0.2, hspace=0.05)
        cbar_ax = fig.add_axes([0.13, 0.05, 0.85, 0.02])
        # cax = fig.add_subplot(ax[5, :])
        fig.colorbar(con, cax=cbar_ax, orientation="horizontal", label="log$[\\tau_{growth} (Myr)]$")
        fig.tight_layout(rect=[0, 0.07, 1, 1])


def plot_dust_contours(fig, ax, radii, a, sigma_dust_1D, model_num=0):
    R, A = np.meshgrid(radii, a)
    levels = np.linspace(-7, 0, 8)                  
    print("Plotting dust size contour maps....")

    if "mig" in sim:
        coli = 1
        miglabel = "$\\bf{Migrating}$ \n"
    else:
        coli = 0
        miglabel = "$\\bf{Stationary}$ \n"

    ax0 = ax[int(s/2),coli]

    sigmas = sigma_dust_1D
    con = ax0.contourf(R, A, np.log10(sigmas), cmap="magma", levels=levels)
    # con2 = ax0.contour(R, A, np.log10(sigmas), linewidths=2, linestyles="dotted", levels=[-2])
    ax0.set_facecolor("k")
    ax0.set_ylim(np.min(a), np.max(a))
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    if coli == 0:
        ax0.set_ylabel("a (cm)")
    if int(s/2) == int(len(sims)/2)-1:
        ax0.set_xlabel("Radius (AU)")

    ax0.plot(radii, a_St1, c='white', alpha=0.9, label="St=1", linewidth=1.5)
    ax0.plot(radii, a_drift, c='limegreen', alpha=0.9, label="$a_{{drift}}$", linewidth=1.5)
    ax0.plot(radii, a_frag, c='deepskyblue', alpha=0.9, label="$a_{{frag}}$", linewidth=1.5)
    ax0.grid(which="major", color="lightgrey", linewidth=0.5)

    if planets:
        planetmass = re.search(r"Mp(\d+)_", sims[model_num]).group(1)    # get planet mass from file path
        for rp in rps:
            ax0.axvline(rp, linestyle='dashed', color='white')
        if int(s/2) == 0:
            ax0.set_title(miglabel + f"{planetmass} $M_\oplus$")
        else:
            ax0.set_title(f"{planetmass} $M_\oplus$")
        ax0.set_xlim(5.1, 140)

    if int(s) == 0:
        ax0.legend()
    if int(s) == len(sims) - 1:
        # cbar = fig.colorbar(
        #     con,
        #     ax=ax[int(s/2), :],          # span exactly the same width
        #     orientation="horizontal",
        #     pad=0.15,
        #     fraction=0.05,
        #     label=r"log$[\tau_{growth}\ (\mathrm{Myr})]$"
        # )
        # fig.tight_layout()

        fig.subplots_adjust(bottom=0.2, hspace=0.05)
        cbar_ax = fig.add_axes([0.072, 0.05, 0.92, 0.02])
        # cax = fig.add_subplot(ax[5, :])
        fig.colorbar(con, cax=cbar_ax, orientation="horizontal", label="log$[\Sigma (g/cm^{{2}})]$")
        fig.tight_layout(rect=[0, 0.07, 1, 1])


def plot_2D_sigma(fig, ax, radii, sigma_gas, sigma_dust, model_num=0):
    print("Plotting 2D surface density....")

    if "mig" in sim:
        coli = 1
        miglabel = "$\\bf{Migrating}$ \n"
    else:
        coli = 0
        miglabel = "$\\bf{Stationary}$ \n"

    # Axes for gas and dust sigma
    ax_g = ax[int(s/2),coli*2]
    ax_d = ax[int(s/2),(coli*2)+1]

    R, PHI = np.meshgrid(radii,phis)
    x = R*np.cos(PHI)
    y = R*np.sin(PHI)

    con_g = ax_g.pcolormesh(x, y, np.log10(sigma_gas.T), cmap="afmhot", shading="auto", vmin=-0.4, vmax=1.5)
    ax_g.set_xlim(-np.max(radii), np.max(radii))
    ax_g.set_ylim(-np.max(radii), np.max(radii))
    # ax_g.set_aspect("equal")
    ax_g.set_xlim(-80,80)
    ax_g.set_ylim(-80,80)
    if coli == 0:
        ax_g.set_ylabel("y (AU)")

    sigma_dust_tot = np.sum(sigma_dust, axis=0)
    con_d = ax_d.pcolormesh(x, y, np.log10(sigma_dust_tot.T), cmap="afmhot", shading="auto", vmin=-2.4, vmax=-0.5)
    ax_d.set_xlim(-np.max(radii), np.max(radii))
    ax_d.set_ylim(-np.max(radii), np.max(radii))
    # ax_d.set_aspect("equal")
    ax_d.set_xlim(-80,80)
    ax_d.set_ylim(-80,80)

    # Plot planet locations 
    # if planets:
    #     planetmass = re.search(r"Mp(\d+)_", sims[model_num]).group(1)    # get planet mass from file path
    #     ax_g.scatter([xp], [yp], color='g', marker='.')
    #     ax_d.scatter([xp], [yp], color='g', marker='.')

    fig.tight_layout()

    if int(s/2) == 0:
        title = ax_g.set_title(miglabel + f"{planetmass}$M_\oplus$ \n", loc="center")
        title.set_position((1.3, 1.1))
    else:
        title = ax_g.set_title(f"{planetmass}$M_\oplus$ \n", loc="center")
        title.set_position((1.3, 0.95))


    if int(s) in [len(sims) - 1, len(sims) - 2]:
        fig.subplots_adjust(bottom=0.15, hspace=0.3)
        for i in [0,1]:
            bottom_ax = [ax_g, ax_d][i]
            bottom_ax.set_xlabel("x (AU)")
            left_edge = bottom_ax.get_position().x0
            right_edge = bottom_ax.get_position().x1
            bottom_edge = bottom_ax.get_position().y0
            width = right_edge - left_edge

            cax = fig.add_axes([left_edge, bottom_edge - 0.07, width, 0.02])
            con = [con_g, con_d][i]
            label = ["gas", "dust"][i]
            # formatter = LogFormatter(10, labelOnlyBase=False) 
            fig.colorbar(con, cax=cax, orientation="horizontal", label=f"log$[\Sigma_{{{label}}}\ (g/cm^2)]$")
            # ax[0,0].legend(loc="lower left")


def calculate_ring_mass(radii, dust_mass_tot, rp, s):
    if "mig" in sim:
        coli = 1
    else:
        coli = 0

    planetmass = int(re.search(r"Mp(\d+)_", sim).group(1))    # get planet mass from file path
    r_hill = rp*((planetmass*1e-6)**(1/3))

    # Stat models: sum over all radii from rp*1.1 to rp*1.6
    # if "stat" in sim:
    i_inner = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rp)+r_hill*2)))
    i_outer = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rp)+r_hill*20)))

    # else:
    #     # Mig models: sum over all radii between rp*1.1 and 50 AU
    #     i_inner = min(range(len(radii)), key=lambda i: abs(radii[i]-rp*1.1))
    #     i_outer = min(range(len(radii)), key=lambda i: abs(radii[i]-rp*1.8))

    # print((1.07*rp)+r_hill*3, (1.07*rp)+r_hill*20)
    dust_mass_ring = dust_mass_tot[i_inner:i_outer]
    ring_mass = np.sum(dust_mass_ring)    # in Earth masses
    print(ring_mass)

    # Store value in array
    ring_masses[coli,int(s/2)] = ring_mass


def plot_ring_mass(fig, ax, mps, ring_masses):
    print(mps[::2], ring_masses)
    ax.scatter(mps[::2], ring_masses[0], color="k")
    ax.plot(mps[::2], ring_masses[0], color="k", label="Stationary")
    ax.scatter(mps[::2], ring_masses[1], color="royalblue")
    ax.plot(mps[::2], ring_masses[1], color="royalblue", label="Migrating")

    ax.set_xlabel("Planet mass ($M_\oplus$)")
    ax.set_ylabel("Ring mass ($M_\oplus$)")
    ax.legend()
    
    fig.tight_layout()


def plot_SI_thresholds(fig, ax, radii, a, sigma_gas_1D, sigma_dust_1D, model_num=0):
    R, A = np.meshgrid(radii, a)
    levels = np.linspace(0, 1, 9)
    print("Plotting Z/Z_crit maps....")

    if "mig" in sim:
        coli = 1
        miglabel = "$\\bf{Migrating}$"
    else:
        coli = 0
        miglabel = "$\\bf{Stationary}$"

    ax0 = ax[int(s/2),coli]

    Z = sigma_dust_1D/sigma_gas_1D     # 70x635/635 = 70x635
    St = np.pi*rhodust*A/(sigma_gas_1D*2)    #70x635/635
    Omega_k = radii**(-3/2)            # Keplerian frequency 
    tau_s = St/Omega_k                 # Stopping time
    Z_crit = 10**((0.1*((np.log10(tau_s))**2)) + (0.07*np.log10(tau_s)) - 2.36)          # Eq. 15 from Lim et al. 2025b
    # Z_crit = 0.15*(np.log10(alpha)**2) - 0.24*np.log10(St)*np.log10(alpha) - 1.48*np.log10(St) + 1.18*np.log10(alpha)
    Z0 = Z/np.abs(Z_crit)

    con = ax0.contourf(R, A, Z0, cmap="GnBu", levels=levels, extend="max")
    # con2 = ax0.contour(R, A, np.log10(sigmas), linewidths=2, linestyles="dotted", levels=[-2])
    # ax0.set_facecolor("k")
    ax0.set_ylim(np.min(a), np.max(a))
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    if coli == 0:
        ax0.set_ylabel("a (cm)")
    if int(s/2) == int(len(sims)/2)-1:
        ax0.set_xlabel("Radius (AU)")

    if planets:
        planetmass = re.search(r"Mp(\d+)_", sims[model_num]).group(1)    # get planet mass from file path
        for rp in rps:
            ax0.axvline(rp, linestyle='dashed', color='black')
        if int(s/2) == 0:
            ax0.set_title(miglabel)
        # else:
        if coli == 0:
            ax0.text(41,10, f"{planetmass} $M_\oplus$", fontsize="x-large")
        # r_hill = rp*((int(planetmass)*1e-6)**(1/3))
        ax0.set_xlim(rp-1, rp+20)

    if int(s) == len(sims) - 1:
        # cbar = fig.colorbar(
        #     con,
        #     ax=ax[int(s/2), :],          # span exactly the same width
        #     orientation="horizontal",
        #     pad=0.15,
        #     fraction=0.05,
        #     label=r"log$[\tau_{growth}\ (\mathrm{Myr})]$"
        # )
        # fig.tight_layout()

        fig.subplots_adjust(bottom=0.2, hspace=0.05)
        cbar_ax = fig.add_axes([0.099, 0.05, 0.87, 0.02])
        # cax = fig.add_subplot(ax[5, :])
        fig.colorbar(con, cax=cbar_ax, orientation="horizontal", label="Z/Z$_{crit}$", extend="max")
        fig.tight_layout(rect=[0, 0.07, 1, 1])



# =================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate 1D plots', prefix_chars='-')
    
    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/astro/phrkvg/simulations/highres_models/"],help="working directory containing simulations")
    parser.add_argument('-sims', metavar='sims', type=str, nargs="*", default=[] ,help="simulation directories containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images/comparison_plots" ,help="directory to save plots to")
    parser.add_argument('-o', metavar='outputs',default=[30], type=int, nargs=1 ,help="output to plot")
    parser.add_argument('-plots', metavar='plots',default=[], type=str, nargs="*" ,help="plots to produce")
    parser.add_argument('-noplanet', action="store_false")
    parser.add_argument('-nobounds', action="store_false")
    parser.add_argument('-nogrog', action="store_false")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-style', metavar='style', type=str, nargs="*", default=["publication"], help="style sheet to apply to plots")

    args = parser.parse_args()
    o = args.o[0]
    plots = args.plots
    wd = args.wd[0]
    sims = args.sims     # Put in sequence Mp_stat, Mp_mig, Mp_stat....etc.
    planets = args.noplanet
    colour_bounds = args.nobounds
    grog = args.nogrog
    plot_window = args.plot_window
    plots_savedir = args.savedir
    style = args.style

    # cm = plt.get_cmap('viridis')
    # colour_cycler = cm(np.linspace(0, 1, 5))
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
    if "growth" in plots:   # only specify for grog models
        fig_growth, ax_growth = plt.subplots(figsize=(7,18), nrows=int(len(sims)/2), ncols=2, sharex=True, sharey=True)
    if "dcon" in plots:   # only specify for grog models
        fig_con, ax_con = plt.subplots(figsize=(13,18), nrows=int(len(sims)/2), ncols=2, sharex=True, sharey=True)
    if "2dsig" in plots:   # only specify for grog models
        fig_2d, ax_2d = plt.subplots(figsize=(10,int(len(sims)/2)*3-1), nrows=int(len(sims)/2), ncols=4, sharex=True, sharey=True)
    if "rmass" in plots:
        fig_m, ax_m = plt.subplots(figsize=(6,4))
    if "si" in plots:   # only specify for grog models
        fig_si, ax_si = plt.subplots(figsize=(9,12), nrows=int(len(sims)/2), ncols=2, sharey=True)


    # ================== Read in data at timesteps =======================

    # Array to store ring masses
    ring_masses = np.zeros((2,5))
    mps = []

    for s, sim in enumerate(sims):
        print (f"============ Plotting data for {sim} =============")
        params_file = f'{wd+sim}/variables.par'
        params_dict = {}

        plotsizex = 3
        plotsizey = 2

        # if ("10_" in sim) or ("20_" in sim):
        #     output = int(o/5)   # need to account for different timestepping in different models
        output = o

        planetmass = re.search(r"Mp(\d+)_", sim).group(1)
        mps.append(int(planetmass))

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
            xps = np.zeros((num_planets))
            yps = np.zeros((num_planets))
            for n in range(num_planets):
                planet_data = np.unique(np.loadtxt(f"{wd+sim}/planet{n}.dat"), axis=0)
                if "Mp160_stat" in sim:
                    xp, yp = planet_data[-1,1], planet_data[-1,2]
                else:
                    xp, yp = planet_data[output,1], planet_data[output,2]
                xps[n] = xp
                yps[n] = yp
                # print(xp, yp)
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
        sigma_dust0 = np.zeros((ndust, nrad, nphi))
        dust_mass = np.zeros((ndust, nrad))
        # dust_mass0 = np.zeros((ndust, nrad))

        gasfile = f"/gasdens{output}.dat" 
        sigma_gas = np.fromfile(wd+sim+gasfile).reshape(nrad,nphi)/(1.125e-7)         # convert back to g/cm2
        n_grains = np.zeros((ndust))
        for n in np.arange(ndust):
            dust_file = f"/dustdens{n}_{output}.dat"
            # sigma_dust0[n] = np.fromfile(wd+sim+f"/dustdens{n}_1.dat").reshape(nrad,nphi)/(1.125e-7)
            sigma_dust[n] = np.fromfile(wd+sim+dust_file).reshape(nrad,nphi)/(1.125e-7)

        sigma_dust_azimsum = np.sum(sigma_dust, axis=2)                  # sum over all phi 
        # sigma_dust0_azimsum = np.sum(sigma_dust0, axis=2)                # sum over all phi 
        sigma_gas_azimsum = np.sum(sigma_gas, axis=1)                    # sum over all phi   

        sigma_gas_1D = sigma_gas_azimsum/nphi                           # dimensions: (nrad) 
        sigma_dust_1D = sigma_dust_azimsum/nphi                         # dimensions: (ndust, nrad)  
        # sigma_dust0_1D = sigma_dust0_azimsum/nphi                       # dimensions: (ndust, nrad)  

        # dust mass for dust of size a as a function of r
        dust_mass = [2*np.pi*radii*sigma_dust_1D[n,:]*delta_r*333030*(1.125e-7) for n in range(ndust)]
        # dust_mass0 = [2*np.pi*radii*sigma_dust0_1D[n,:]*delta_r*333030*(1.125e-7) for n in range(ndust)]
        gas_mass = 2*np.pi*radii*sigma_gas_1D*delta_r*333030*(1.125e-7)      # convert from Msun to Mearth
        dust_mass_tot = np.sum(dust_mass, axis=0)
        # dust_mass0_tot = np.sum(dust_mass0, axis=0)

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
        if "dcon" in plots:
            plot_dust_contours(fig_con, ax_con, radii, a, sigma_dust_1D, s)
        if "2dsig" in plots:
            plot_2D_sigma(fig_2d, ax_2d, radii, sigma_gas, sigma_dust, s)
        if "rmass" in plots:
            calculate_ring_mass(radii, dust_mass_tot, rps[0], s)       # assume only 1 planet here
        if "si" in plots:
            plot_SI_thresholds(fig_si, ax_si, radii, a, sigma_gas_1D, sigma_dust_1D, s)


    # ================ Plot ring masses, if selected ===================
    if "rmass" in plots:
        plot_ring_mass(fig_m, ax_m, mps, ring_masses)
        fig_m.savefig(f"{plots_savedir}/ring_masses.png")


    # ======================== Generate Plots ==========================
    print(f"-------------------\nPlotting comparison plots for {sims}\n=============")

    if "gsig" in plots:
        fig_gas_sigma.savefig(f"{plots_savedir}/fullcomp_sigmagas.png")
    if "dgr" in plots:
        fig_dgr.savefig(f"{plots_savedir}/fullcomp_dgr.png")
    if "growth" in plots:
        fig_growth.savefig(f"{plots_savedir}/fullcomp_growth.png")
    if "dcon" in plots:
        fig_con.savefig(f"{plots_savedir}/fullcomp_contours.png")
    if "2dsig" in plots:
        fig_2d.savefig(f"{plots_savedir}/fullcomp_2Dsigma.png")
    if "si" in plots:
        fig_si.savefig(f"{plots_savedir}/fullcomp_SI_thresholds.png")

    if plot_window:
        plt.show()
