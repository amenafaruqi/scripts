import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
plt.style.use('default')


# ======================= Misc functions ==============================
def sigma_at_r(r=10, t_index=-1):
    r_in_range = radii[(r*0.8<radii) & (radii<r*1.2)]                # values of r within % of chosen radius
    mini = np.where(radii == np.min(r_in_range))[0][0]
    maxi = np.where(radii == np.max(r_in_range))[0][0]
    sigma_in_r_range = sigma_dust_avg[t_index, :, mini:maxi]         # for chosen timestep, all grain sizes
    sigma_at_r0 = np.median(sigma_in_r_range, axis=1)                # avg sigma at chosen radius for all grain sizes
    return sigma_at_r0, mini, maxi

# ======================= Dust Contour Plots ==============================

def plot_dust_contours():
    fig0 = plt.figure(figsize=(16,12))
    R, A = np.meshgrid(radii, a)
    # levels = np.linspace(-11,1,7)                    # Brauer 2008 levels
    # levels = np.linspace(-7, 2, 10)                  # Birnstiel 2012 levels 
    levels = np.linspace(-18, 6, 13)                   
    print("Plotting dust size contour maps....")

    for i, o in enumerate(outputs):
        ax0 = fig0.add_subplot(plotsizey, plotsizex, i+1)
        sigmas = sigma_dust_avg[i]
        c = ax0.contourf(R, A, np.log10(sigmas), cmap="Greys", levels=levels)
        ax0.set_ylim(np.min(a), np.max(a))
        ax0.set_xscale("log")
        ax0.set_yscale("log")
        ax0.set_title(f"{round(timesteps[i],3)} Myr")
        ax0.plot(radii, a_St1[i], c='black', alpha=0.7, label="St=1")
        ax0.plot(radii, a_drift[i], c='deepskyblue', alpha=0.7, label="$a_{{drift}}$")
        ax0.plot(radii, a_frag[i], c='red', alpha=0.7, label="$a_{{frag}}$")
        if not i%plotsizex:
            ax0.set_ylabel("a (cm)")
        else:
            ax0.set_yticks([])
        if i < plotsizex and len(timesteps) > plotsizex:
            ax0.set_xticks([])
        else:
            ax0.set_xlabel("R (AU)")

    fig0.tight_layout()

    fig0.subplots_adjust(right=0.89, hspace=0.3)
    cbar_ax = fig0.add_axes([0.91, 0.53, 0.02, 0.4])
    fig0.colorbar(c, cax=cbar_ax, orientation="vertical", label="log$[\Sigma (g/cm^{{2}})]$")
    ax0.legend(loc="upper right")
    fig0.savefig(f"{plots_savedir}/{sim}_contour.png")

# ======================= Dust-Gas Ratio =========================

def plot_dustgasratio():
    fig, ax = plt.subplots(figsize=(7,5))
    ax.set_prop_cycle(color=[cm(1.*i/8) for i in range(0,len(timesteps)+1)])
    print("Plotting dust-gas ratio....")

    for i,o in enumerate(timesteps):
        dustgasratio = sigma_dust_tot_avg[i]/sigma_gas_avg[i]
        ax.plot(radii, dustgasratio, label=f"t={round(o, 3)} Myr")
    
    if planets:
        for rp in rps:
            ax.axvline(rp, linestyle='dashed', color='black')

    ax.legend()
    ax.set_xlim(np.min(radii), np.max(radii))
    ax.set_xlabel("R (AU)")
    ax.set_ylabel("dust-gas ratio")
    ax.set_xscale("log")
    ax.set_yscale("log")

    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_dustgasratio.png")

# ================== Gas Sigma Profile ===================

def plot_gas_sigma():
    fig, ax = plt.subplots(figsize=(7,6))
    ax.set_prop_cycle(color=[cm(1.*i/8) for i in range(0,len(timesteps)+1)])
    print("Plotting gas surface density....")

    for i, o in enumerate(timesteps):
        ax.plot(radii, sigma_gas_avg[i], label=f"{round(o,3)} Myr")

    if planets:
        for rp in rps:
            ax.axvline(rp, linestyle='dashed', color='black') 

    ax.set_xlabel("R (AU)")
    ax.set_ylabel("$\Sigma_{gas} (g/cm^{2})$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlim(np.min(radii), np.max(radii))

    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_sigmagas.png")

# ================== Dust Sigma Profile ===================

def plot_dust_sigma():
    fig, ax = plt.subplots(figsize=(7,6))
    ax.set_prop_cycle(color=[cm(1.*i/8) for i in range(0,len(timesteps)+1)])
    print("Plotting dust surface density....")

    for i, o in enumerate(timesteps):
        ax.plot(radii, sigma_dust_tot_avg[i], label=f"{round(o,3)} Myr")

    if planets:
        for rp in rps:
            ax.axvline(rp, linestyle='dashed', color='black') 

    ax.set_xlabel("R (AU)")
    ax.set_ylabel("$\Sigma_{dust} (g/cm^{2})$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlim(np.min(radii), np.max(radii))

    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_sigmadust.png")

# ================== Dust sigma evolution over time (for different a) ===================

def plot_sigma_at_r(rs=[10]):
    for r in rs:
        fig, ax = plt.subplots(1, dpi=150)
        ax.set_prop_cycle(color=[cm(1.*i/8) for i in range(len(timesteps)+1)])
        print(f"Plotting dust surface density at R={r}AU....")

        for t in timesteps:
            c = next(ax._get_lines.prop_cycler)['color']
            sigma_at_rau, mini, maxi  = sigma_at_r(r)
            ax.plot(a, sigma_at_rau, linestyle='solid', label=f"t={round(t,3)} Myr", color=c)
            a_frag_t = np.median(a_frag[i,mini:maxi])
            a_drift_t = np.median(a_drift[i,mini:maxi])
            a1 = np.min(a_frag_t, a_drift_t)
            ax.axvline(a1, linestyle='dashed', color=c)

        ax.set_xlim(np.min(a), np.max(a))
        ax.set_xlabel("a (cm)")
        ax.set_ylabel(f"$\Sigma$ at {r}AU (g/cm$^{{2}}$)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend()
        fig.savefig(f"{plots_savedir}/{sim}_sigmaatr.png")

# ================= Total dust mass evolution over time  ===============

def plot_dust_mass():
    fig, ax = plt.subplots(1, dpi=150)
    ax.set_prop_cycle(color=[cm(1.*i/8) for i in range(0,len(timesteps)+1)])
    print("Plotting total dust mass....")

    delta_r = radii[1]-radii[0]   # works for a linear grid only!!
    M_disc = np.zeros(len(timesteps))

    for i in range(len(timesteps)):
        M_cell = sigma_dust_tot[i]*(radii*1.5e13)*2*np.pi*(delta_r*1.5e13)
        M_disc[i] = np.sum(M_cell)

    ax.plot(timesteps, M_disc)
    ax.set_xlabel("t (Myr)")
    ax.set_ylabel("Total  $M_{{dust}} (g/cm^{{2}})$")
    ax.set_yscale("log")
    ax.set_xscale("log")
    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_Mdust.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate 1D plots', prefix_chars='-')
    
    parser.add_argument('-o', metavar='outputs',default=[], type=int, nargs="*" ,help="outputs to plot")
    parser.add_argument('-noplanet', action="store_false")
    parser.add_argument('-nogrog', action="store_false")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/astro/phrkvg/simulations/planet_growth/test_models"],help="working directory containing simulations")
    parser.add_argument('-sim', metavar='sim', type=str, nargs=1, default=["planettest"] ,help="simulation directory containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images" ,help="directory to save plots to")

    args = parser.parse_args()
    outputs = args.o
    wd = args.wd[0]
    sim = args.sim[0]
    print(f"Plotting outputs {outputs} for {sim}\n =============")
    simdir = f"{wd}/{sim}/"
    planets = args.noplanet
    grog = args.nogrog
    plot_window = args.plot_window
    plots_savedir = args.savedir

    if plot_window:
        matplotlib.use('TkAgg')

    # ================== Read in data at timesteps =======================

    params_file = f'{simdir}/variables.par'
    params_dict = {}

    plotsizex = 3
    plotsizey = int(len(outputs)/plotsizex)+1

    cm = plt.get_cmap('gist_rainbow')

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
    if grog:
        mingsize = float(params_dict['MIN_GRAIN_SIZE'])
        maxgsize = float(params_dict['MAX_GRAIN_SIZE'])
        nss_coag = int(params_dict['NUMSUBSTEPS_COAG'])
        rhodust = float(params_dict['RHO_DUST'])

    densfloor = float(params_dict['DENSITY_FLOOR'])
    dt_orbits = int(float(params_dict['DT'])/(2*np.pi))   # 2pi = 1 orbit = 1 yr
    ninterm = float(params_dict['NINTERM'])               # number of dts between outputs
    dt_outputs = dt_orbits*ninterm                        # time between outputs
    timesteps = np.array(outputs)*dt_outputs*1e-6    # time in Myr

    if planets:
        planets_data = np.genfromtxt(f"{simdir}/planet.cfg").reshape(1,6)
        rps = planets_data[:,1]
        Mps = planets_data[:,2]
        planet_periods = (rps**3)**0.5    # orbital period of planet in yrs
        # planet_orbits = timesteps/planet_periods

    # FARGO initialises grains with sizes uniformly distributed across ndust bins in logspace
    if grog:
        a = np.logspace(np.log10(mingsize),    
                        np.log10(maxgsize),
                        ndust+1)

        a = (0.5*(a[1:] + a[:-1]))                                       # convert cm to um by multiplying by 1e4
    r_cells = np.loadtxt(f'{simdir}/domain_y.dat')[3:-3]             # ignore ghost cells
    phi_cells = np.loadtxt(f'{simdir}/domain_x.dat')[3:-3]
    radii = np.array([(r_cells[n]+r_cells[n+1])/2 for n in range(len(r_cells)-1)])
    phis = np.array([(phi_cells[n]+phi_cells[n+1])/2 for n in range(len(phi_cells)-1)])

    sigma_gas = np.zeros((len(timesteps), nrad, nphi))
    sigma_dust = np.zeros((len(timesteps), ndust, nrad, nphi))

    for i,t in enumerate(outputs):
        gasfile = f"gasdens{int(t)}.dat" 
        sigma_gas[i] = np.fromfile(simdir+gasfile).reshape(nrad,nphi)/(1.125e-7)         # convert back to g/cm2
        if grog:
            for n in np.arange(ndust):
                dust_file = f"dustdens{n}_{int(t)}.dat" 
                sigma_dust[i,n] = np.fromfile(simdir+dust_file).reshape(nrad,nphi)/(1.125e-7)


    sigma_dust_azimsum = np.sum(sigma_dust, axis=3)                  # sum over all phi 
    sigma_gas_azimsum = np.sum(sigma_gas, axis=2)                    # sum over all phi   
    sigma_dust_tot = np.sum(sigma_dust_azimsum, axis=1)              # sum over all dust sizes

    surfdens = sigma_gas[i]
    dens_first_wedge = surfdens[:,0].reshape(nrad,1)
    dens_additional = np.concatenate((surfdens[:,:],dens_first_wedge),axis=1)
    sigma_gas_avg = sigma_gas_azimsum/nphi
    sigma_dust_avg = sigma_dust_azimsum/nphi
    sigma_dust_tot_avg = sigma_dust_tot/nphi

    if grog:
        uf = 10                                               # fragmentation velocity
        hr = hr0*(radii**f)                                   # aspect ratio
        a_St1 = (2/np.pi)*(sigma_gas_avg/rhodust)             # plot St=1 line

        # size of largest grains in a fragmentation-dominated distribution
        a_frag = 100*(2/(3*np.pi))*((uf**2)/(rhodust*1000*alpha))*(hr**-2)*(sigma_gas_avg*10/((2e30)*(6.67e-11)))*(radii*1.5e11)  # from Birnstiel+2012
        # size of largest grains in a drift-dominated distribution
        a_drift = 100*(2/(rhodust*1000*np.pi))*(hr**-2)*sigma_dust_tot*10*(2/3)   # from Birnstiel+2012

    # ======================== Generate Plots ==========================

    if not plot_window:
        plt.style.use(['../styles/publication.mplstyle'])

    plot_gas_sigma()

    if grog:
        plot_dust_contours()
        plot_dust_sigma()
        plot_dustgasratio()
        # plot_sigma_at_r()

    if plot_window:
        plt.show()