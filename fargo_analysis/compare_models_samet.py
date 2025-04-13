import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter, NullFormatter
import argparse
import re
plt.style.use('default')


# ====================== Gas Sigma ========================

def overlay_gas_sigmas(fig, ax, radii, sigma_gas_1D, model_num=0):
    print("Plotting gas surface density....")

    legend_elements = []
    color = colour_cycler[model_num]
    ax.plot(radii, sigma_gas_1D, color=color)

    for m in range(model_num+1):
        if grog:
            planetmass = re.search(r"Mp(\d+)_", sims[m]).group(1)    # get planet mass from file path
        else:
            planetmass = re.search(r"(\d+)Me", sims[m]).group(1)
        legend_elements.append(Line2D([0], [0], color=colour_cycler[m], label=f"{planetmass} $M_\oplus$"))

    if planets:
        if "stat" in sim:
            planetcolour = 'k'
        else:
            planetcolour = color
        for rp in rps:
            ax.axvline(rp, linestyle='dashed', color=planetcolour)
    
    ax.set_xlabel("R (AU)")
    ax.set_ylabel("$\Sigma_{gas} (g/cm^{2})$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(np.min(radii), np.max(radii))
    ax.legend(handles=legend_elements)
    fig.tight_layout()


# ================= Dust Sigma ===================

def overlay_dust_sigmas(fig, ax, radii, sigma_dust, model_num=0):
    print("Plotting dust distribution by grain size....")
    legend_elements = []
    color = colour_cycler[model_num]

    if grog:
        n_size_decades = int(np.log10(maxgsize) - np.log10(mingsize))   # assumes min and max g size are the same for both models!
        size_decades = np.split(np.arange(ndust), n_size_decades)
        subps = "ABCDEFG"

        for n, size_decade in enumerate(size_decades):
            dust_sigma_binned_tot = np.sum(sigma_dust[size_decades[n],:],axis=0)

            ax[subps[n]].plot(radii, dust_sigma_binned_tot, color=color)

            if planets:
                if "stat" in sim:
                    planetcolour = 'k'
                else:
                    planetcolour = color
                for rp in rps:
                    ax[subps[n]].axvline(rp, linestyle='dashed', color=planetcolour)

            dustsizes = [(10**n) * mingsize, (10**(n+1)) * mingsize]
            dustsizes = [np.format_float_positional(d,3,fractional=False,unique=True) for d in dustsizes]
            ax[subps[n]].set_title(f"{dustsizes[0]}-{dustsizes[1]}cm")
            ax[subps[n]].set_xscale("log")
            ax[subps[n]].set_yscale("log")
            ax[subps[n]].set_xlim(np.min(radii), np.max(radii))
            
        for m in range(model_num+1):
            planetmass = re.search(r"Mp(\d+)_", sims[m]).group(1)    # get planet mass from file path
            legend_elements.append(Line2D([0], [0], color=colour_cycler[m], label=f"{planetmass} $M_\oplus$"))

        ax["A"].set_ylabel("$\Sigma_{dust} (g/cm^{2})$")
        ax["D"].set_ylabel("$\Sigma_{dust} (g/cm^{2})$")
        ax["F"].set_ylabel("$\Sigma_{dust} (g/cm^{2})$")
        ax["A"].set_xlabel("R (AU)")
        ax["B"].set_xlabel("R (AU)")
        ax["C"].set_xlabel("R (AU)")
        ax["F"].set_xlabel("R (AU)")
        ax["G"].set_xlabel("R (AU)")
        # ax["D"].set_xticks([])
        # ax["E"].set_xticks([])
        ax["G"].legend(loc="lower left", handles=legend_elements)
        fig.tight_layout()


    # No coagulation case
    else:
        ax = ax.flatten()
        for n in range(ndust):
            ax[n].plot(radii, sigma_dust[n], color=color)

            if planets:
                for rp in rps:
                    ax[n].axvline(rp, linestyle='dashed', color="k")

            ax[n].set_title(f"St={round(stokes[n],3)}")
            ax[n].set_yscale("log")
            # ax[n].set_xscale("log")
            ax[n].set_xlim(0.7, 1.5)
            
        for m in range(model_num+1):
            planetmass = re.search(r"(\d+)Me", sims[m]).group(1)    # get planet mass from file path
            legend_elements.append(Line2D([0], [0], color=colour_cycler[m], label=f"{planetmass} $M_\oplus$"))


        ax[0].set_ylabel("$\Sigma_{dust}$")
        ax[3].set_ylabel("$\Sigma_{dust}$")
        ax[3].set_xlabel("R")
        ax[4].set_xlabel("R")
        ax[5].set_xlabel("R")
        ax[0].set_xticks([])
        ax[1].set_xticks([])
        ax[2].set_xticks([])
        ax[5].legend(handles=legend_elements)
        fig.tight_layout()

        # Plot total dust sigma
        dust_mass_tot = np.sum(stokes**3)
        mass_weighted_sum = np.sum([sigma_dust[n]*(stokes[n]**3) for n in range(ndust)], axis=0)
        mass_weighted_avg = mass_weighted_sum/dust_mass_tot
        # sigma_dust_tot = np.sum(sigma_dust, axis=0)                  # dimensions: (nrad)  
        ax[5].plot(radii, mass_weighted_avg, color=color)
        ax[5].set_title("All St")
        ax[5].set_yscale("log")
        # ax[5].set_xscale("log")
        ax[5].set_xlim(0.7,1.5)

        if planets:
            for rp in rps:
                ax[5].axvline(rp, linestyle='dashed', color="k")



# =================== Total Dust Mass =====================

def overlay_total_dust_mass(fig, ax, radii, dust_mass_tot, model_num=0):
    legend_elements = []
    color = colour_cycler[model_num]
    ax.plot(radii, dust_mass_tot, color=color)

    for m in range(model_num+1):
        planetmass = re.search(r"Mp(\d+)_", sims[m]).group(1)    # get planet mass from file path
        legend_elements.append(Line2D([0], [0], color=colour_cycler[m], label=f"{planetmass} $M_\oplus$"))

    if planets:
        if "stat" in sim:
            planetcolour = 'k'
        else:
            planetcolour = color
        for rp in rps:
            ax.axvline(rp, linestyle='dashed', color=planetcolour)
    
    ax.set_xlabel("R (AU)")
    ax.set_ylabel("$M_{dust} (M_\oplus)$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(np.min(radii), np.max(radii))
    ax.legend(handles=legend_elements)
    fig.tight_layout()


# ================= Plot final dust contours =====================

def plot_dust_contours(fig, ax, radii, a, sigma_dust_1D, model_num=0):
    R, A = np.meshgrid(radii, a)
    levels = np.linspace(-18, 6, 13) 
                  
    print("Plotting dust size contour maps....")
    ax0 = ax.flatten()[model_num]
    sigmas = sigma_dust_1D
    con = ax0.contourf(R, A, np.log10(sigmas), cmap="Greys", levels=levels)
    ax0.set_ylim(np.min(a), np.max(a))
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.plot(radii, a_St1, c='black', alpha=0.7, label="St=1")
    ax0.plot(radii, a_drift, c='deepskyblue', alpha=0.7, label="$a_{{drift}}$")
    ax0.plot(radii, a_frag, c='red', alpha=0.7, label="$a_{{frag}}$")
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
        ax0.set_xlabel("R (AU)")

    fig.tight_layout()
    
    if model_num == len(sims)-1:
        ax_cbar = ax.flatten()[model_num+1]
        ax_cbar.remove()
        fig.subplots_adjust(right=0.89, hspace=0.3)
        cax = fig.add_axes([ax0.get_position().x1+0.01,ax0.get_position().y0,0.02,ax0.get_position().height])
        fig.colorbar(con, cax=cax, orientation="vertical", label="log$[\Sigma (g/cm^{{2}})]$")
        ax0.legend(loc="upper right")


# ================= Plot dust-gas ratios =====================

def overlay_dust_gas_ratio(fig, ax, radii, dust_mass_tot, gas_mass, model_num=0):
    legend_elements = []
    color = colour_cycler[model_num]
    epsilon = dust_mass_tot/gas_mass
    ax.plot(radii, epsilon, color=color)

    for m in range(model_num+1):
        planetmass = re.search(r"Mp(\d+)_", sims[m]).group(1)    # get planet mass from file path
        legend_elements.append(Line2D([0], [0], color=colour_cycler[m], label=f"{planetmass} $M_\oplus$"))

    if planets:
        if "stat" in sim:
            planetcolour = 'k'
        else:
            planetcolour = color
        for rp in rps:
            ax.axvline(rp, linestyle='dashed', color=planetcolour)
    
    ax.set_xlabel("R (AU)")
    ax.set_ylabel("$M_{dust} (M_\oplus)$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(np.min(radii), np.max(radii))
    ax.legend(handles=legend_elements)
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
    plots_savedir = args.savedir[0]
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
        fig_gas_sigma, ax_gas_sigma = plt.subplots(figsize=(6,5))
    if "dcon" in plots:   # only specify dcon for grog models
        fig_con, ax_con = plt.subplots(figsize=(17,12), nrows=2, ncols=3)
    if "dsig" in plots:
        if grog:   # assume 7 dust size decades
            fig_dust_sigma, ax_dust_sigma = plt.subplot_mosaic("AABBCC;DDDEEE;FFFGGG", figsize=(15,13))
        else:      # 2 cols for 3 St, 3 for 5 St
            fig_dust_sigma, ax_dust_sigma = plt.subplots(figsize=(15,8), ncols=3, nrows=2)
    if "dmass" in plots:
        fig_dust_mass, ax_dust_mass = plt.subplots(figsize=(6,5))
    if "dgr" in plots:
        fig_dgr, ax_dgr = plt.subplots(figsize=(6,5))

    # ================== Read in data at timesteps =======================

    for s, sim in enumerate(sims):
        print (f"============ Plotting data for {sim} =============")
        params_file = f'{wd+sim}/variables.par'
        params_dict = {}

        plotsizex = 3
        plotsizey = 2

        outputs = np.array(o)
        # outputs = outputs.astype(int)

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
                xp, yp = planet_data[o,1], planet_data[o,2]
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

        gasfile = f"/gasdens{o}.dat" 
        sigma_gas = np.fromfile(wd+sim+gasfile).reshape(nrad,nphi)/(1.125e-7)         # convert back to g/cm2
        n_grains = np.zeros((ndust))
        for n in np.arange(ndust):
            dust_file = f"/dustdens{n}_{o}.dat"
            sigma_dust[n] = np.fromfile(wd+sim+dust_file).reshape(nrad,nphi)/(1.125e-7)

        sigma_dust_azimsum = np.sum(sigma_dust, axis=2)                  # sum over all phi 
        sigma_gas_azimsum = np.sum(sigma_gas, axis=1)                    # sum over all phi   

        sigma_gas_1D = sigma_gas_azimsum/nphi                           # dimensions: (nrad) 
        sigma_dust_1D = sigma_dust_azimsum/nphi                         # dimensions: (ndust, nrad)  

        # dust mass for dust of size a as a function of r
        dust_mass = [2*np.pi*radii*sigma_dust_1D[n,:]*delta_r*333030 for n in range(ndust)]
        gas_mass = 2*np.pi*radii*sigma_gas_1D*delta_r*333030      # convert from Msun to Mearth
        dust_mass_tot = np.sum(dust_mass, axis=0)

        if grog:
            uf = 10                                               # fragmentation velocity
            hr = hr0*(radii**f)                                   # aspect ratio
            cs = hr*(((2e30)*(6.67e-11))/(radii*1.5e11))**0.5     # [m/s]
            b = (uf**2/alpha)*(hr**-2)*(radii*1.5e11)/((2e30)*(6.67e-11)) # dimensionless
            p = (sigma_gas_1D*(cs**2)/((2*np.pi)**0.5))*(hr**-1)*((radii*1.5e11)**-1)
            pad = np.empty((1, 1))*np.nan
            gamma = (radii/p)*np.abs(np.append(np.diff(p)/np.diff(radii), pad))
            C = 2/(np.pi*hr)
            a_St1 = (2/np.pi)*(sigma_gas_1D/rhodust)              # plot St=1 line

            # size of largest grains in a fragmentation-dominated distribution
            # a_frag = 100*(2/(3*np.pi))*((uf**2)/(rhodust*1000*alpha))*(hr**-2)*(sigma_gas_1D*10/((2e30)*(6.67e-11)))*(radii*1.5e11)  # from Birnstiel+2012
            a_frag = (sigma_gas_1D/rhodust)*(3-(9-4*(b**2))**0.5)/(np.pi*b)
            # print(np.min(np.sqrt(9-4*(b**2))),np.max(9-4*(b**2)))

            # size of largest grains in a drift-dominated distribution
            # a_drift = 100*(2/(rhodust*1000*np.pi))*(hr**-2)*np.sum(sigma_dust_1D, axis=1)*10*(2/3)   # from Birnstiel+2012 (assume gamma=3/2)
            a_drift = (2/np.pi)*(np.sum(sigma_dust_1D, axis=0)/(rhodust*gamma*hr**2))    
    
        if "gsig" in plots:
            overlay_gas_sigmas(fig_gas_sigma, ax_gas_sigma, radii, sigma_gas_1D, s)
        if "dcon" in plots:
            plot_dust_contours(fig_con, ax_con, radii, a, sigma_dust_1D, s)
        if "dsig" in plots:
            overlay_dust_sigmas(fig_dust_sigma, ax_dust_sigma, radii, sigma_dust_1D, s)
        if "dmass" in plots:
            overlay_total_dust_mass(fig_dust_mass, ax_dust_mass, radii, dust_mass_tot, s)
        if "dgr" in plots:
            overlay_dust_gas_ratio(fig_dgr, ax_dgr, radii, dust_mass_tot, gas_mass, s)
    
    # ======================== Generate Plots ==========================
    print(f"-------------------\nPlotting comparison plots for {sims}\n=============")

    if "gsig" in plots:
        fig_gas_sigma.savefig(f"{plots_savedir}/comparison_sigmagas.png")
    if "dcon" in plots:
        fig_con.savefig(f"{plots_savedir}/comparison_dcontour.png")
    if "dsig" in plots:
        fig_dust_sigma.savefig(f"{plots_savedir}/comparison_sigmadust.png")
    if "dmass" in plots:
        fig_dust_mass.savefig(f"{plots_savedir}/comparison_dustmass.png")
    if "dgr" in plots:
        fig_dgr.savefig(f"{plots_savedir}/comparison_dgr.png")

    if plot_window:
        plt.show()
