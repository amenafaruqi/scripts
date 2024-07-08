import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D
import argparse
plt.style.use('default')

def overlay_gas_sigmas(fig, ax, radii, sigma_gas_1D, model_num=0):
    ax.set_prop_cycle(color=colour_cycler)
    print("Plotting gas surface density....")

    legend_elements = []
    for i, t in enumerate(timesteps):
        color = next(ax._get_lines.prop_cycler)['color']
        ax.plot(radii, sigma_gas_1D[i], label=f"{round(t, 3)} Myr", color=color, linestyle=linestyles[model_num])
        legend_elements.append(Line2D([0], [0], color=color, lw=2, label=f"{round(t, 3)} Myr"))

        if planets:
            for rp in rps[:,i]:
                ax.axvline(rp, linestyle='dashed', color=color)
    for m in range(model_num+1):
        sim = simdirs[m].split('models/')[-1]    # ignore full file path
        legend_elements.append(Line2D([0], [0], color='k', linestyle=linestyles[m], label=f"{sim}"))

    ax.set_xlabel("R (AU)")
    ax.set_ylabel("$\Sigma_{gas} (g/cm^{2})$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(np.min(radii), np.max(radii))
    ax.legend(handles=legend_elements)
    fig.tight_layout()


def overlay_dustgasratio():
    pass

def overlay_dust_mass(fig, ax, radii, dust_mass, model_num=0):
    print("Plotting dust distribution by grain size....")

    n_size_decades = int(np.log10(maxgsize) - np.log10(mingsize))   # assumes these are the same for both models!
    size_decades = np.split(np.arange(ndust), n_size_decades)
    dust_mass_tot_binned = np.zeros((len(outputs), n_size_decades, nrad))
    nrow = 0

    for n, size_decade in enumerate(size_decades):
        for i, t in enumerate(timesteps):
            dust_mass_tot_binned = np.sum(dust_mass[i,size_decade[0]:size_decade[-1],:], axis=0)
            ncol = n%3
            ax[nrow,ncol].set_prop_cycle(color=colour_cycler)
            color = next(ax[nrow,ncol]._get_lines.prop_cycler)['color']
            
            ax[nrow,ncol].plot(radii, np.log10(dust_mass_tot_binned), label=f"{round(t, 3)} Myr", color=color)

            if planets:
                for rp in rps[:,i]:
                    ax[nrow,ncol].axvline(rp, linestyle='dashed', color=color)

        if not n%3:
            ax[nrow,ncol].set_ylabel(f"log[$M_{{dust}} (M_\oplus)$]")
        else:
            ax[nrow,ncol].set_yticks([])
        if n < 3 and n_size_decades > 3:
            ax[nrow,ncol].set_xticks([])
        else:
            ax[nrow,ncol].set_xlabel("R (AU)")

        dustsizes = [(10**n) * mingsize, (10**(n+1)) * mingsize]
        dustsizes = [np.format_float_positional(d,3,fractional=False,unique=True) for d in dustsizes]
        ax[nrow,ncol].set_title(f"{dustsizes[0]}-{dustsizes[1]}cm")
        ax[nrow,ncol].set_xscale("log")
        # ax.set_yscale("log")
        ax[nrow,ncol].set_xlim(min(radii), max(radii))
        if (n+1)%3 == 0:
            nrow += 1
    
    ax[0,0].legend()
    fig.tight_layout()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate 1D plots', prefix_chars='-')

    parser.add_argument('-simdirs', metavar='simdirs', type=str, nargs="*", default=[] ,help="simulation directories containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images" ,help="directory to save plots to")
    parser.add_argument('-o', metavar='outputs',default=[], type=int, nargs="*" ,help="outputs to plot")
    parser.add_argument('-noplanet', action="store_false")
    parser.add_argument('-nogrog', action="store_false")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-style', metavar='style', type=str, nargs="*", default=["publication"], help="style sheet to apply to plots")

    args = parser.parse_args()
    outputs = args.o
    simdirs = [f"/home/astro/phrkvg/simulations/{sim}/" for sim in args.simdirs]
    planets = args.noplanet
    grog = args.nogrog
    plot_window = args.plot_window
    plots_savedir = args.savedir
    style = args.style

    cm = plt.get_cmap('gist_rainbow')
    colour_cycler = [cm(1.*i/5) for i in range(0,len(outputs)+1)]
    linestyles = ['solid', 'dashdot', 'dotted']

    if plot_window:
        matplotlib.use('TkAgg')
    
    if not plot_window:
        for s in style:
            plt.style.use([f"../styles/{s}.mplstyle"])

    # =================== Define figures ========================
    fig_gas_sigma, ax_gas_sigma = plt.subplots(figsize=(6,5))
    fig_dust_mass, ax_dust_mass = plt.subplots(3,3,figsize=(17,16))
    
    # ================== Read in data at timesteps =======================

    for s, simdir in enumerate(simdirs):
        print (f"============ Plotting data for {simdir} =============")
        params_file = f'{simdir}/variables.par'
        params_dict = {}

        plotsizex = 3
        plotsizey = int(len(outputs)/plotsizex)+1

        # Load model params into dict
        param_lines = open(params_file).readlines()
        for line in param_lines:
            if line.split():
                param_label, param_value = line.split()[0:2]
                params_dict.update([(param_label, param_value)])

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
        dt_orbits = int(float(params_dict['DT'])/(2*np.pi))   # 2pi = 1 orbit = 1 yr
        ninterm = float(params_dict['NINTERM'])               # number of dts between outputs
        dt_outputs = dt_orbits*ninterm                        # time between outputs
        timesteps = np.array(outputs)*dt_outputs*1e-6         # time in Myr

        if planets:
            with open(f"{simdir}/planet.cfg", 'r') as pfile:
                num_planets = len(pfile.readlines()) - 5        # ignore 5 lines of headers
            
            rps = np.zeros((num_planets, len(outputs)))
            for n in range(num_planets):
                planet_data = np.unique(np.loadtxt(f"{simdir}/planet{n}.dat"), axis=0)
                xp, yp = planet_data[outputs][:,1], planet_data[outputs][:,2]
                # print(xp,yp)
                rps[n] = ((xp**2) + (yp**2))**0.5

            planet0_period = (rps[0]**3)**0.5                 # orbital period of planet in yrs
            planet_orbits = timesteps*1e6/planet0_period

        # FARGO initialises grains with sizes uniformly distributed across ndust bins in logspace
        if grog:
            a = np.logspace(np.log10(mingsize),    
                            np.log10(maxgsize),
                            ndust+1)

            a = (0.5*(a[1:] + a[:-1]))                                   # grain sizes in middles of bins (in cm)

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
        sigma_gas = np.zeros((len(outputs), nrad, nphi))
        sigma_dust = np.zeros((len(outputs), ndust, nrad, nphi))
        dust_mass = np.zeros((len(outputs), ndust, nrad))
        gas_mass = np.zeros((len(outputs), nrad))

        for i,t in enumerate(outputs):
            gasfile = f"gasdens{int(t)}.dat" 
            sigma_gas[i] = np.fromfile(simdir+gasfile).reshape(nrad,nphi)/(1.125e-7)         # convert back to g/cm2
            n_grains = np.zeros((ndust))
            if grog:
                for n in np.arange(ndust):
                    dust_file = f"dustdens{n}_{int(t)}.dat"
                    sigma_dust[i,n] = np.fromfile(simdir+dust_file).reshape(nrad,nphi)/(1.125e-7)
        sigma_dust_azimsum = np.sum(sigma_dust, axis=3)                  # sum over all phi 
        sigma_gas_azimsum = np.sum(sigma_gas, axis=2)                    # sum over all phi   

        sigma_gas_1D = sigma_gas_azimsum/nphi                           # dimensions: (noutputs, nrad) 
        sigma_dust_1D = sigma_dust_azimsum/nphi                         # dimensions: (noutputs, ndust, nrad)   
        # sigma_dust_sum_1D = avgdustdens_azimsum/nphi                  # dimensions: (noutputs, nrad)
        
        for i,t in enumerate(outputs):
            # dust mass for dust of size a as a function of r
            dust_mass[i,:,:] = [2*np.pi*radii*sigma_dust_1D[i,n,:]*delta_r*333030 for n in range(ndust)]
            gas_mass[i,:] = 2*np.pi*radii*sigma_gas_1D[i,:]*delta_r*333030      # convert from Msun to Mearth

        dust_mass_tot = np.sum(dust_mass, axis=1)

        # Get dust velocities
        # if grog:
        #     v_dust = np.zeros((len(outputs), ndust, 2, nrad, nphi))   # additional dimension of 2 for x and y velocity
        #     for i,t in enumerate(outputs):
        #         # for n in np.arange(ndust):
        #         n=30
        #         dust_file_x = f"dustvx{n}_{int(t)}.dat"
        #         dust_file_y = f"dustvy{n}_{int(t)}.dat"
        #         v_dust[i,n,0] = np.fromfile(simdir+dust_file_x).reshape(nrad,nphi)   # vx
        #         v_dust[i,n,1] = np.fromfile(simdir+dust_file_y).reshape(nrad,nphi)   # vy

        # if grog:
        #     uf = 10                                               # fragmentation velocity
        #     hr = hr0*(radii**f)                                   # aspect ratio
        #     cs = hr*(((2e30)*(6.67e-11))/(radii*1.5e11))**0.5     # [m/s]
        #     b = (uf**2/alpha)*(hr**-2)*(radii*1.5e11)/((2e30)*(6.67e-11)) # dimensionless
        #     p = (sigma_gas_1D*(cs**2)/((2*np.pi)**0.5))*(hr**-1)*((radii*1.5e11)**-1)
        #     pad = np.empty((len(timesteps), 1))*np.nan
        #     gamma = (radii/p)*np.abs(np.append(np.diff(p)/np.diff(radii), pad, axis=1))
        #     C = 2/(np.pi*hr)
        #     a_St1 = (2/np.pi)*(sigma_gas_1D/rhodust)              # plot St=1 line

        #     # size of largest grains in a fragmentation-dominated distribution
        #     # a_frag = 100*(2/(3*np.pi))*((uf**2)/(rhodust*1000*alpha))*(hr**-2)*(sigma_gas_1D*10/((2e30)*(6.67e-11)))*(radii*1.5e11)  # from Birnstiel+2012
        #     a_frag = (sigma_gas_1D/rhodust)*(3-(9-4*(b**2))**0.5)/(np.pi*b)
        #     # print(np.min(np.sqrt(9-4*(b**2))),np.max(9-4*(b**2)))

        #     # size of largest grains in a drift-dominated distribution
        #     # a_drift = 100*(2/(rhodust*1000*np.pi))*(hr**-2)*np.sum(sigma_dust_1D, axis=1)*10*(2/3)   # from Birnstiel+2012 (assume gamma=3/2)
        #     a_drift = (2/np.pi)*(np.sum(sigma_dust_1D, axis=1)/(rhodust*gamma*hr**2))    
    
        overlay_dust_mass(fig_dust_mass, ax_dust_mass, radii, dust_mass, s)
        overlay_gas_sigmas(fig_gas_sigma, ax_gas_sigma, radii, sigma_gas_1D, s)

    # ======================== Generate Plots ==========================
    print(f"-------------------\nPlotting comparison plots for {simdirs}\n=============")

    # fig_dust_mass.savefig(f"{plots_savedir}/comparison_graindist.png")
    fig_gas_sigma.savefig(f"{plots_savedir}/comparison_gassigma.png")

    if plot_window:
        plt.show()
