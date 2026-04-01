import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import re

plt.style.use('default')

import warnings
warnings.filterwarnings("ignore")

# ====================== Gas Sigma ========================

def overlay_gas_sigmas(fig, ax, radii, sigma_gas_1D):
    print("Plotting gas surface density....")
    if "stat" in sim:
        colour = "k"
        if len(rps) == 1:
            label = "S1"
            ls = "solid"
        else:
            label = "S2"
            ls = "dashdot"
    else:
        colour = "royalblue"
        if len(rps) == 1:
            label = "M1"
            ls = "solid"
        else:
            label = "M2"
            ls = "dashdot"

    # planet_masses = re.findall(r"Mp(\d+)_", sim)    # get planet masses from file path
    # if len(planet_masses) == 2:
    #     planet1_mass = planet_masses[1]             # outer planet mass, if present
    # planet0_mass = planet_masses[0]                 # inner planet mass, assume always present

    ax.plot(radii, sigma_gas_1D, color=colour, label=label, linestyle=ls)
    for rp in rps:
        ax.axvline(rp, linestyle='dashed', color=colour)
        
    ax.set_xlabel("Radius (AU)")
    ax.set_ylabel("$\Sigma_{gas} (g/cm^{2})$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(np.min(radii), np.max(radii))
    ax.legend()
    fig.tight_layout()


# ================= Plot dust-gas ratios =====================

def overlay_dust_gas_ratio(fig, ax, radii, dust_mass_tot, gas_mass):
    epsilon = dust_mass_tot/gas_mass
    if "stat" in sim:
        colour = "k"
        if len(rps) == 1:
            label = "S1"
            ls = "solid"
        else:
            label = "S2"
            ls = "dashdot"
    else:
        colour = "royalblue"
        if len(rps) == 1:
            label = "M1"
            ls = "solid"
        else:
            label = "M2"
            ls = "dashdot"

    # planet_masses = re.findall(r"Mp(\d+)_", sim)    # get planet masses from file path
    # if len(planet_masses) == 2:
    #     planet1_mass = planet_masses[1]             # outer planet mass, if present
    # planet0_mass = planet_masses[0]                 # inner planet mass, assume always present

    ax.plot(radii, epsilon, color=colour, label=label, linestyle=ls)
    for rp in rps:
        ax.axvline(rp, linestyle='dashed', color=colour)
    
    ax.set_xlabel("Radius (AU)")
    ax.set_ylabel("Dust-to-gas ratio")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(np.min(radii), np.max(radii))
    ax.legend()
    fig.tight_layout()

# ================= Plot growth timescales =====================

def plot_tgrowth(fig, ax, radii, a, sigma_dust_1D):
    R, A = np.meshgrid(radii, a)
    levels = np.linspace(-4, 16, 6)                   
    print("Plotting growth timescales....")
    # Set labels and subplots
    if "stat" in sim:
        coli = 0
        if len(rps) == 1:
            label = "S1"
            rowi = 0
        else:
            label = "S2"
            rowi = 1
    else:
        coli = 1
        if len(rps) == 1:
            label = "M1"
            rowi = 0
        else:
            label = "M2"
            rowi = 1
    
    ax0 = ax[rowi,coli]

    hr = hr0*(radii**0.25)
    h = hr*radii
    St = A*rhodust*np.pi/(sigma_gas_1D*2)
    cs = h*2*np.pi/(radii**0.5)
    tau_growth = rhodust*A*h/(sigma_dust_1D*cs*(alpha*St*3))
    tau_growth = tau_growth*1e-6   # convert to Myr

    con = ax0.contourf(R,A, np.log10(tau_growth), cmap="YlGnBu", levels=levels, extend="both")
    
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.set_title(label)

    if coli == 0:
        ax0.set_ylabel("a (cm)")
    if rowi == 1:
        ax0.set_xlabel("Radius (AU)")

    ax0.set_ylim(np.min(a), np.max(a))

    for rp in rps:
            ax0.axvline(rp, linestyle='dashed', color='black')
    if (rowi, coli) == (1,1):
        fig.subplots_adjust(bottom=0.2, hspace=0.05)
        cbar_ax = fig.add_axes([0.13, 0.05, 0.85, 0.02])
        # cax = fig.add_subplot(ax[5, :])
        fig.colorbar(con, cax=cbar_ax, orientation="horizontal", label="log$[\\tau_{growth} (Myr)]$")
        fig.tight_layout(rect=[0, 0.07, 1, 1])


def plot_dust_contours(fig, ax, radii, a, sigma_dust_1D):
    R, A = np.meshgrid(radii, a)
    levels = np.linspace(-7, 1, 9)                  
    print("Plotting dust size contour maps....")

    # Set labels and subplots
    if "stat" in sim:
        coli = 0
        if len(rps) == 1:
            label = "S1"
            rowi = 0
        else:
            label = "S2"
            rowi = 1
    else:
        coli = 1
        if len(rps) == 1:
            label = "M1"
            rowi = 0
        else:
            label = "M2"
            rowi = 1
    
    ax0 = ax[rowi,coli]

    sigmas = sigma_dust_1D
    con = ax0.contourf(R, A, np.log10(sigmas), cmap="Greys", levels=levels)
    ax0.set_ylim(np.min(a), np.max(a))
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    if coli == 0:
        ax0.set_ylabel("a (cm)")
    if rowi == 1:
        ax0.set_xlabel("Radius (AU)")

    ax0.plot(radii, a_St1, c='black', alpha=0.7, label="St=1")
    ax0.plot(radii, a_drift, c='deepskyblue', alpha=0.7, label="$a_{{drift}}$")
    ax0.plot(radii, a_frag, c='red', alpha=0.7, label="$a_{{frag}}$")

    for rp in rps:
        ax0.axvline(rp, linestyle='dashed', color='black')

    ax0.set_title(label)
    fig.tight_layout()

    if (rowi, coli) == (1,1):
        fig.subplots_adjust(bottom=0.3, hspace=0.15)
        cbar_ax = fig.add_axes([0.075, 0.17, 0.91, 0.03])
        fig.colorbar(con, cax=cbar_ax, orientation="horizontal", label="log$[\Sigma_{dust} (g/cm^{{2}})]$")
        # fig.tight_layout(rect=[0, 0.07, 1, 1])


def plot_2D_sigma(fig, ax, radii, sigma_gas, sigma_dust):
    print("Plotting 2D surface density....")

    # Set labels and subplots
    if "stat" in sim:
        if len(rps) == 1:
            label = "S1"
            rowi = 0
        else:
            label = "S2"
            rowi = 1
    else:
        if len(rps) == 1:
            label = "M1"
            rowi = 2
        else:
            label = "M2"
            rowi = 3
    
    # Axes for gas and dust sigma
    ax_g = ax[rowi,0]
    ax_d = ax[rowi,1]

    R, PHI = np.meshgrid(radii,phis)
    x = R*np.cos(PHI)
    y = R*np.sin(PHI)

    # Plot gas surface density
    con_g = ax_g.pcolormesh(x, y, np.log10(sigma_gas.T), cmap="afmhot", shading="auto", vmin=-0.4, vmax=1.5)
    ax_g.set_xlim(-np.max(radii), np.max(radii))
    ax_g.set_ylim(-np.max(radii), np.max(radii))
    # ax_g.set_aspect("equal")
    ax_g.set_xlim(-90,90)
    ax_g.set_ylim(-90,90)
    ax_g.set_ylabel("y (AU)")
    ax_g.set_aspect("equal")

    # Plot dust surface density
    sigma_dust_tot = np.sum(sigma_dust, axis=0)
    con_d = ax_d.pcolormesh(x, y, np.log10(sigma_dust_tot.T), cmap="afmhot", shading="auto", vmin=-2.4, vmax=-0.5)
    ax_d.set_xlim(-np.max(radii), np.max(radii))
    ax_d.set_ylim(-np.max(radii), np.max(radii))
    # ax_d.set_aspect("equal")
    ax_d.set_xlim(-90,90)
    ax_d.set_ylim(-90,90)
    ax_d.set_aspect("equal")

    # Plot planet locations 
    for pi in np.arange(len(rps)):
        ax_g.scatter([xps[pi]], [yps[pi]], color='g', marker='.')
        ax_d.scatter([xps[pi]], [yps[pi]], color='g', marker='.')

    if rowi == 0:
        ax[rowi, 0].set_title("Gas \n" + label, loc="center")
        ax[rowi, 1].set_title("Dust \n" + label, loc="center")

    else:
        ax_g.set_title(label, loc="center")
        ax_d.set_title(label, loc="center")

    fig.tight_layout()

    if rowi == 3:
        # fig.colorbar(con_g, ax=ax[rowi, 0], orientation="horizontal")
        # fig.colorbar(con_d, ax=ax[rowi, 1], orientation="horizontal")
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
            fig.colorbar(con, cax=cax, orientation="horizontal", label=f"log$[\Sigma_{{{label}}}\ (g/cm^2)]$")


# def calculate_ring_mass(radii, dust_mass_tot, rp, s):
#     if "mig" in sim:
#         coli = 1
#     else:
#         coli = 0

#     planetmass = int(re.search(r"Mp(\d+)_", sim).group(1))    # get planet mass from file path
#     r_hill = rp*((planetmass*1e-6)**(1/3))

#     # Stat models: sum over all radii from rp*1.1 to rp*1.6
#     # if "stat" in sim:
#     i_inner = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rp)+r_hill*2)))
#     i_outer = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rp)+r_hill*20)))

#     # else:
#     #     # Mig models: sum over all radii between rp*1.1 and 50 AU
#     #     i_inner = min(range(len(radii)), key=lambda i: abs(radii[i]-rp*1.1))
#     #     i_outer = min(range(len(radii)), key=lambda i: abs(radii[i]-rp*1.8))

#     # print((1.07*rp)+r_hill*3, (1.07*rp)+r_hill*20)
#     dust_mass_ring = dust_mass_tot[i_inner:i_outer]
#     ring_mass = np.sum(dust_mass_ring)    # in Earth masses
#     print(ring_mass)

#     # Store value in array
#     ring_masses[coli,int(s/2)] = ring_mass


# def plot_ring_mass(fig, ax, mps, ring_masses):
#     print(mps[::2], ring_masses)
#     ax.scatter(mps[::2], ring_masses[0], color="k")
#     ax.plot(mps[::2], ring_masses[0], color="k", label="Stationary")
#     ax.scatter(mps[::2], ring_masses[1], color="royalblue")
#     ax.plot(mps[::2], ring_masses[1], color="royalblue", label="Migrating")

#     ax.set_xlabel("Planet mass ($M_\oplus$)")
#     ax.set_ylabel("Ring mass ($M_\oplus$)")
#     ax.legend()
    
#     fig.tight_layout()



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
    o = args.o[-1]                    # Plot single timestep plots at last timestep given
    tevo_o = args.o                   # List of timesteps for time evo plots 
    plots = args.plots
    wd = args.wd[0]
    sims = args.sims
    plot_window = args.plot_window
    plots_savedir = args.savedir
    style = args.style

    cm = plt.get_cmap('viridis')
    colour_cycler = cm(np.linspace(0, 1, 7))       # For plotting different mass decades

    if plot_window:
        matplotlib.use('TkAgg')
    else:
        for s in style:
            plt.style.use([f"../styles/{s}.mplstyle"])

    # =================== Define figures and axes ========================
    if "gsig" in plots:
        fig_gas_sigma, ax_gas_sigma = plt.subplots(figsize=(6,4))
    if "dgr" in plots:
        fig_dgr, ax_dgr = plt.subplots(figsize=(6,4))
    if "growth" in plots:   # only specify for grog models
        fig_growth, ax_growth = plt.subplots(figsize=(7,7), nrows=2, ncols=2, sharex=True, sharey=True)
    if "dcon" in plots:   # only specify for grog models
        fig_con, ax_con = plt.subplots(figsize=(8,8), nrows=2, ncols=2, sharex=True, sharey=True)
    if "2dsig" in plots:   # only specify for grog models
        fig_2d, ax_2d = plt.subplots(figsize=(7,14), nrows=4, ncols=2, sharex=True, sharey=True)
    # if "rmass" in plots:
    #     fig_m, ax_m = plt.subplots(figsize=(6,4))

    # ================== Read in data at timesteps =======================

    # Array to store ring masses
    ring_masses = np.zeros((2,5))
    mps = []

    for s, sim in enumerate(sims):
        print (f"============ Plotting data for {sim} =============")
        params_file = f'{wd+sim}/variables.par'
        params_dict = {}
        output = o

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
        t = o*dt_outputs*1e-6                                 # time in Myr

        # ----------------------------------------------------

        with open(f"{wd+sim}/planet.cfg", 'r') as pfile:
            num_planets = len(pfile.readlines()) - 5        # ignore 5 lines of headers
        
        rps = np.zeros((num_planets))
        xps = np.zeros((num_planets))
        yps = np.zeros((num_planets))
        for n in range(num_planets):
            planet_data = np.unique(np.loadtxt(f"{wd+sim}/planet{n}.dat"), axis=0)
            xp, yp = planet_data[output,1], planet_data[output,2]
            xps[n] = xp
            yps[n] = yp
            rps[n] = ((xp**2) + (yp**2))**0.5


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
            overlay_gas_sigmas(fig_gas_sigma, ax_gas_sigma, radii, sigma_gas_1D)
        if "dgr" in plots:
            overlay_dust_gas_ratio(fig_dgr, ax_dgr, radii, dust_mass_tot, gas_mass)
        if "growth" in plots:
            plot_tgrowth(fig_growth, ax_growth, radii, a, sigma_dust_1D)
        if "dcon" in plots:
            plot_dust_contours(fig_con, ax_con, radii, a, sigma_dust_1D)
        if "2dsig" in plots:
            plot_2D_sigma(fig_2d, ax_2d, radii, sigma_gas, sigma_dust)
        # if "rmass" in plots:
        #     calculate_ring_mass(radii, dust_mass_tot, rps[0])       # assume only 1 planet here

    # ================ Plot ring masses, if selected ===================
    # if "rmass" in plots:
    #     plot_ring_mass(fig_m, ax_m, mps, ring_masses)
    #     fig_m.savefig(f"{plots_savedir}/ring_masses.png")


    # ======================== Generate Plots ==========================
    print(f"-------------------\nPlotting comparison plots for {sims}\n=============")

    if "gsig" in plots:
        fig_gas_sigma.savefig(f"{plots_savedir}/SPF_sigmagas.png")
    if "dgr" in plots:
        fig_dgr.savefig(f"{plots_savedir}/SPF_dgr.png")
    if "growth" in plots:
        fig_growth.savefig(f"{plots_savedir}/SPF_growth.png")
    if "dcon" in plots:
        fig_con.savefig(f"{plots_savedir}/SPF_contours.png")
    if "2dsig" in plots:
        fig_2d.savefig(f"{plots_savedir}/SPF_2Dsigma.png")

    if plot_window:
        plt.show()
