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
    sigma_in_r_range = sigma_dust_1D[t_index, :, mini:maxi]         # for chosen timestep, all grain sizes
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
        sigmas = sigma_dust_1D[i]
        con = ax0.contourf(R, A, np.log10(sigmas), cmap="Greys", levels=levels)
        ax0.set_ylim(np.min(a), np.max(a))
        ax0.set_xscale("log")
        ax0.set_yscale("log")
        ax0.plot(radii, a_St1[i], c='black', alpha=0.7, label="St=1")
        ax0.plot(radii, a_drift[i], c='deepskyblue', alpha=0.7, label="$a_{{drift}}$")
        ax0.plot(radii, a_frag[i], c='red', alpha=0.7, label="$a_{{frag}}$")
        
        if not i%plotsizex:
            ax0.set_ylabel("a (cm)")
        else:
            ax0.set_yticks([])
        if i < plotsizex and len(outputs) > plotsizex:
            ax0.set_xticks([])
        else:
            ax0.set_xlabel("R (AU)")
        
        if not p_orbits:
            ax0.set_title(f"{round(timesteps[i],3)} Myr")
        elif planets:
            ax0.set_title(f"{int(round(planet_orbits[i],0))} orbits")

    fig0.tight_layout()

    fig0.subplots_adjust(right=0.89, hspace=0.3)
    cbar_ax = fig0.add_axes([0.91, 0.53, 0.02, 0.4])
    fig0.colorbar(con, cax=cbar_ax, orientation="vertical", label="log$[\Sigma (g/cm^{{2}})]$")
    ax0.legend(loc="upper right")
    fig0.savefig(f"{plots_savedir}/{sim}_contour.png")

# ======================= Dust-Gas Ratio =========================

def plot_dustgasratio():
    fig, ax = plt.subplots(figsize=(7,5))
    ax.set_prop_cycle(color=[cm(1.*i/8) for i in range(0,len(timesteps)+1)])
    print("Plotting dust-gas ratio....")

    for i,t in enumerate(timesteps):
        dustgasratio = sigma_dust_sum_1D[i]/sigma_gas_1D[i]
        if not p_orbits:
            tlabel = f"{round(t, 3)} Myr"
        elif planets:
            tlabel = f"{int(round(planet_orbits[i], 0))} orbits"

        ax.plot(radii, dustgasratio, label=tlabel)
    
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

    for i, t in enumerate(timesteps):
        if not p_orbits:
            tlabel = f"{round(t, 3)} Myr"
        elif planets:
            tlabel = f"{int(round(planet_orbits[i], 0))} orbits"

        ax.plot(radii, sigma_gas_1D[i], label=tlabel)

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

    for i, t in enumerate(timesteps):
        if not p_orbits:
            tlabel = f"{round(t, 3)} Myr"
        elif planets:
            tlabel = f"{int(round(planet_orbits[i], 0))} orbits"

        ax.plot(radii, sigma_dust_sum_1D[i], label=tlabel)

    if planets:
        for rp in rps:
            ax.axvline(rp, linestyle='dashed', color='black') 

    ax.set_xlabel("R (AU)")
    ax.set_ylabel("$ \Sigma_{dust} (g/cm^{2})$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlim(np.min(radii), np.max(radii))

    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_sigmadust.png")

# ================== Dust sigma evolution over time (for different a) ===================

def plot_sigma_at_r(rs=[25]):
    for r in rs:
        fig, ax = plt.subplots(1, dpi=150)
        ax.set_prop_cycle(color=[cm(1.*i/8) for i in range(len(timesteps)+1)])
        print(f"Plotting dust surface density at R={r}AU....")

        for i,t in enumerate(timesteps):
            c = next(ax._get_lines.prop_cycler)['color']
            sigma_at_rau, mini, maxi  = sigma_at_r(r, i)
            ax.plot(a, sigma_at_rau, linestyle='solid', label=f"t={round(t,3)} Myr", color=c)
            # a_frag_t = np.median(a_frag[i,mini:maxi])
            # a_drift_t = np.median(a_drift[i,mini:maxi])
            # a1 = np.min(a_frag_t, a_drift_t)
            # ax.axvline(a1, linestyle='dashed', color=c)

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
        M_cell = sigma_dust_sum_1D[i]*(radii*1.5e13)*2*np.pi*(delta_r*1.5e13)
        M_disc[i] = np.sum(M_cell)

    ax.plot(timesteps, M_disc)
    ax.set_xlabel("t (Myr)")
    ax.set_ylabel("Total $M_{{dust}} (g/cm^{{2}})$")
    ax.set_yscale("log")
    ax.set_xlim(min(timesteps), max(timesteps))
    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_Mdust.png")

# ================= Radial drift timescale  ===============

def plot_t_drift():
    fig, ax = plt.subplots(figsize=(7,6))
    ax.set_prop_cycle(color=[cm(1.*i/8) for i in range(0,len(timesteps)+1)])
    print("Plotting drift timescale....")
    for i,t in  enumerate(timesteps):
        if not p_orbits:
            tlabel = f"{round(t, 3)} Myr"
        elif planets:
            tlabel = f"{int(round(planet_orbits[i], 0))} orbits"
        ax.plot(radii[1:-1], t_drift[i][1:-1], label=tlabel)
    
    ax.legend()
    ax.set_yscale("log")
    ax.set_xlabel("R (AU)")
    ax.set_ylabel("$\\tau_{drift} (Myr)$")

    fig.savefig(f"{plots_savedir}/{sim}_tdrift.png")


# ==========================================================


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate 1D plots', prefix_chars='-')

    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/astro/phrkvg/simulations/planet_growth/test_models"],help="working directory containing simulations")
    parser.add_argument('-sim', metavar='sim', type=str, nargs=1, default=["planettest"] ,help="simulation directory containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images" ,help="directory to save plots to")
    parser.add_argument('-o', metavar='outputs',default=[], type=int, nargs="*" ,help="outputs to plot")
    parser.add_argument('-noplanet', action="store_false")
    parser.add_argument('-nogrog', action="store_false")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-porbits', action="store_true")

    args = parser.parse_args()
    outputs = args.o
    wd = args.wd[0]
    sim = args.sim[0]
    simdir = f"{wd}/{sim}/"
    planets = args.noplanet
    grog = args.nogrog
    plot_window = args.plot_window
    plots_savedir = args.savedir
    p_orbits = args.porbits

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
    timesteps = np.array(outputs)*dt_outputs*1e-6         # time in Myr

    if planets:
        planets_data = np.genfromtxt(f"{simdir}/planet.cfg").reshape(1,6)
        rps = planets_data[:,1]
        Mps = planets_data[:,2]
        planet0_period = (rps[0]**3)**0.5                 # orbital period of planet in yrs
        planet_orbits = timesteps*1e6/planet0_period

    # FARGO initialises grains with sizes uniformly distributed across ndust bins in logspace
    if grog:
        a = np.logspace(np.log10(mingsize),    
                        np.log10(maxgsize),
                        ndust+1)

        a = (0.5*(a[1:] + a[:-1]))                                   # grain sizes in middles of bins (in cm)

    r_cells = np.loadtxt(f'{simdir}/domain_y.dat')[3:-3]             # ignore ghost cells
    phi_cells = np.loadtxt(f'{simdir}/domain_x.dat')[3:-3]
    radii = np.array([(r_cells[n]+r_cells[n+1])/2 for n in range(len(r_cells)-1)])
    phis = np.array([(phi_cells[n]+phi_cells[n+1])/2 for n in range(len(phi_cells)-1)])

    sigma_gas = np.zeros((len(outputs), nrad, nphi))
    sigma_dust = np.zeros((len(outputs), ndust, nrad, nphi))
    avg_dust_dens = np.zeros((len(outputs), nrad, nphi))
    sum_dustvol = np.zeros((len(outputs)))                      # total volume of dust at each timestep

    for i,t in enumerate(outputs):
        gasfile = f"gasdens{int(t)}.dat" 
        sigma_gas[i] = np.fromfile(simdir+gasfile).reshape(nrad,nphi)/(1.125e-7)         # convert back to g/cm2
        n_grains = np.zeros((ndust))
        if grog:
            for n in np.arange(ndust):
                dust_file = f"dustdens{n}_{int(t)}.dat" 
                sigma_dust[i,n] = np.fromfile(simdir+dust_file).reshape(nrad,nphi)/(1.125e-7)
                n_grains[n] = np.sum(sigma_dust[i,n])*3/(a[n]*rhodust)                    # number of dust grains at size a and time t
                avg_dust_dens[i] = avg_dust_dens[i] + sigma_dust[i,n]*(a[n]**3)           # (mass)volume-weighted sum of densities
            sum_dustvol[i] = np.sum(n_grains*(a**3))                                      # numbers of grains of each size * volume of grain
            avg_dust_dens[i] = avg_dust_dens[i]/sum_dustvol[i]                            

    avgdustdens_azimsum = np.sum(avg_dust_dens, axis=2)              # sum over all phi
    sigma_dust_azimsum = np.sum(sigma_dust, axis=3)                  # sum over all phi 
    sigma_gas_azimsum = np.sum(sigma_gas, axis=2)                    # sum over all phi   

    sigma_gas_1D = sigma_gas_azimsum/nphi                           # dimensions: (noutputs, nrad) 
    sigma_dust_1D = sigma_dust_azimsum/nphi                         # dimensions: (noutputs, ndust, nrad)
    sigma_dust_sum_1D = avgdustdens_azimsum/nphi                    # dimensions: (noutputs, nrad)

    if grog:
        uf = 10                                               # fragmentation velocity
        hr = hr0*(radii**f)                                   # aspect ratio
        cs = hr*(((2e30)*(6.67e-11))/(radii*1.5e11))**0.5
        b = (uf**2/alpha)*(hr**-2)*(radii*1.5e11)/((2e30)*(6.67e-11)) # dimensionless
        p = (sigma_gas_1D*(cs**2)/((2*np.pi)**0.5))*(hr**-1)*((radii*1.5e11)**-1)
        pad = np.empty((len(timesteps), 1))
        gamma = (radii/p)*np.abs(np.append(np.diff(p)/np.diff(radii), pad, axis=1))
        C = 2/(np.pi*hr)
        a_vals = np.logspace(np.log10(mingsize),    
                             np.log10(maxgsize),
                             nrad)
        a_St1 = (2/np.pi)*(sigma_gas_1D/rhodust)              # plot St=1 line

        # size of largest grains in a fragmentation-dominated distribution
        a_frag = 100*(2/(3*np.pi))*((uf**2)/(rhodust*1000*alpha))*(hr**-2)*(sigma_gas_1D*10/((2e30)*(6.67e-11)))*(radii*1.5e11)  # from Birnstiel+2012
        # a_frag = (sigma_gas_1D/rhodust)*(3-(9-4*(b**2))**0.5)/(np.pi*b)

        # size of largest grains in a drift-dominated distribution
        a_drift = 100*(2/(rhodust*1000*np.pi))*(hr**-2)*np.sum(sigma_dust_1D, axis=1)*10*(2/3)   # from Birnstiel+2012
        # a_drift = (2/np.pi)*(sigma_dust_sum_1D/rhodust)*(1/gamma)*(hr**-2)
        
        t_drift = ((C/gamma)*(sigma_gas_1D/(rhodust*a_vals))*(radii*1.5e11/cs))*3.17e-14    # drift timescale in Myr
    # ======================== Generate Plots ==========================

    if not plot_window:
        plt.style.use(['../styles/publication.mplstyle'])

    print(f"Plotting outputs {outputs} for {sim}\n =============")

    plot_gas_sigma()

    if grog:
        plot_dust_contours()
        plot_dust_sigma()
        plot_dustgasratio()
        # plot_dust_mass()
        # plot_sigma_at_r([rps[0]+5])
        plot_t_drift()

    if plot_window:
        plt.show()
