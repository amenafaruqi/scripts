import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
plt.style.use('default')

import warnings
warnings.filterwarnings("ignore")

# ======================= Misc functions ==============================
def sigma_at_r(r=10, t_index=-1):
    r_in_range = radii[(r*0.8<radii) & (radii<r*1.2)]                # values of r within % of chosen radius
    mini = np.where(radii == np.min(r_in_range))[0][0]
    maxi = np.where(radii == np.max(r_in_range))[0][0]
    sigma_in_r_range = sigma_dust_1D[t_index, :, mini:maxi]         # for chosen timestep, all grain sizes
    sigma_at_r0 = np.median(sigma_in_r_range, axis=1)                # avg sigma at chosen radius for all grain sizes
    return sigma_at_r0, mini, maxi

def calculate_e(r, v):
    mu = 4*(np.pi**2)              # in code units, otherwise mu=4pi**2 
    h = np.cross(r, v)
    e_vec = (np.cross(v, h)/mu) - (r/np.linalg.norm(r))   # vector e
    e = np.linalg.norm(e_vec)
    return e
    
    # a = np.norm(h)/(mu*(1-(np.norm(e)**2)))   # scalar
    # cos_E = (1 - np.norm(r)/a)/np.norm(e)    # scalar
    # sin_E = np.dot(r, v)/(np.norm(e)*np.norm(r))   # scalar
    # E = np.arctan2(sin_E, cos_E)


# ======================= Dust Contour Plots ==============================

def plot_dust_contours():
    fig0 = plt.figure(figsize=(17,14))
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
        if planets:
            for rp in rps[:,i]:
                ax0.axvline(rp, linestyle='dashed', color='black')
            
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
    ax.set_prop_cycle(color=colour_cycler)
    print("Plotting dust-gas ratio....")

    for i,t in enumerate(timesteps):

        dustgasratio = dust_mass_tot[i]/gas_mass[i]
        if not p_orbits:
            tlabel = f"{round(t, 3)} Myr"
        elif planets:
            tlabel = f"{int(round(planet_orbits[i], 0))} orbits"
        
        color = next(ax._get_lines.prop_cycler)['color']
        ax.plot(radii, dustgasratio, label=tlabel, color=color)
    
        if planets:
            for rp in rps[:,i]:
                ax.axvline(rp, linestyle='dashed', color=color)

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
    fig, ax = plt.subplots(figsize=(6,5))
    ax.set_prop_cycle(color=colour_cycler)
    print("Plotting gas surface density....")

    for i, t in enumerate(timesteps):
        if not p_orbits:
            tlabel = f"{round(t, 3)} Myr"
        elif planets:
            tlabel = f"{int(round(planet_orbits[i], 0))} orbits"

        color = next(ax._get_lines.prop_cycler)['color']
        ax.plot(radii, sigma_gas_1D[i], label=tlabel, color=color)

        if planets:
            for rp in rps[:,i]:
                ax.axvline(rp, linestyle='dashed', color=color)

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
    ax.set_prop_cycle(color=colour_cycler)
    print("Plotting dust surface density....")

    for i, t in enumerate(timesteps):
        if not p_orbits:
            tlabel = f"{round(t, 3)} Myr"
        elif planets:
            tlabel = f"{int(round(planet_orbits[i], 0))} orbits"

        color = next(ax._get_lines.prop_cycler)['color']
        ax.plot(radii, sigma_dust_sum_1D[i], label=tlabel, color=color)

    if planets:
        for rp in rps[:,i]:
            ax.axvline(rp, linestyle='dashed', color=color)

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
        ax.set_prop_cycle(color=colour_cycler)
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
    ax.set_prop_cycle(color=colour_cycler)
    print("Plotting total dust mass....")

    for i, t in enumerate(timesteps):
        if not p_orbits:
            tlabel = f"{round(t, 3)} Myr"
        elif planets:
            tlabel = f"{int(round(planet_orbits[i], 0))} orbits"

        color = next(ax._get_lines.prop_cycler)['color']
        ax.plot(radii, dust_mass_tot[i], label=tlabel, color=color)

        if planets:
            for rp in rps[:,i]:
                ax.axvline(rp, linestyle='dashed', color=color)

    ax.set_xlabel("R (AU)")
    ax.set_ylabel("$M_{{dust}} (M_\oplus)$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlim(min(radii), max(radii))
    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_Mdust.png")


def plot_dust_mass_by_grain():
    fig = plt.figure(figsize=(17,16))
    print("Plotting dust distribution by grain size....")

    n_size_decades = int(np.log10(maxgsize) - np.log10(mingsize))
    size_decades = np.split(np.arange(ndust), n_size_decades)
    dust_mass_tot_binned = np.zeros((len(outputs), n_size_decades, nrad))
    psizex = 3
    psizey = int(n_size_decades/psizex)+1

    for n, size_decade in enumerate(size_decades):
        ax = fig.add_subplot(psizey, psizex, n+1)
        ax.set_prop_cycle(color=colour_cycler)

        for i, t in enumerate(timesteps):
            dust_mass_tot_binned = np.sum(dust_mass[i,size_decade[0]:size_decade[-1],:], axis=0)
            color = next(ax._get_lines.prop_cycler)['color']

            if not p_orbits:
                tlabel = f"{round(t, 3)} Myr"
            elif planets:
                tlabel = f"{int(round(planet_orbits[i], 0))} orbits"
            
            ax.plot(radii, np.log10(dust_mass_tot_binned), label=tlabel, color=color)

            if planets:
                for rp in rps[:,i]:
                    ax.axvline(rp, linestyle='dashed', color=color)

        if not n%psizex:
            ax.set_ylabel(f"log[$M_{{dust}} (M_\oplus)$]")
        else:
            ax.set_yticks([])
        if n < psizex and n_size_decades > psizex:
            ax.set_xticks([])
        else:
            ax.set_xlabel("R (AU)")

        dustsizes = [(10**n) * mingsize, (10**(n+1)) * mingsize]
        dustsizes = [np.format_float_positional(d,3,fractional=False,unique=True) for d in dustsizes]
        ax.set_title(f"{dustsizes[0]}-{dustsizes[1]}cm")
        ax.set_xscale("log")
        # ax.set_yscale("log")
        ax.legend()
        ax.set_xlim(min(radii), max(radii))

    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_graindist.png")


def plot_dust_size_distribution():
    fig, ax = plt.subplots(1, dpi=150)
    ax.set_prop_cycle(color=colour_cycler)
    print("Plotting dust mass distribution....")

    for i, t in enumerate(timesteps):
        if not p_orbits:
            tlabel = f"{round(t, 3)} Myr"
        elif planets:
            tlabel = f"{int(round(planet_orbits[i], 0))} orbits"

        color = next(ax._get_lines.prop_cycler)['color']
        dust_mass_in_bin = np.sum(dust_mass[i], axis=1)
        ax.scatter(a, dust_mass_in_bin, label=tlabel, color=color)
        ax.plot(a, dust_mass_in_bin, label=tlabel, color=color)

    ax.set_xlabel("a (cm)")
    ax.set_ylabel("$M_{{dust}} (g/cm^{{2}})$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlim(min(a), max(a))
    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_dustsizes.png")

# ============== Plot dust eccentricities ==================
def plot_dust_e():
    fig = plt.figure(figsize=(17,12))
    R, PHI = np.meshgrid(radii, phis, indexing="ij")   # dimensions: (nrad, nphi)
    X = R*np.cos(PHI)
    Y = R*np.sin(PHI)
    Z = np.zeros((nrad, nphi))
    pos = np.array((X, Y, Z))    # dimensions: (2, nrad, nphi)
    e = np.zeros((len(outputs),nrad,nphi))
    n = 30
    # for n in np.arange(ndust):
    for i, o in enumerate(outputs):
        ax = fig.add_subplot(plotsizey, plotsizex, i+1)
        v_phi = v_dust[i,n,0,:,:]  # dimensions: (nrad, nphi)
        v_r = v_dust[i,n,1,:,:]
        # print(np.min(v_r), np.max(v_r))
        # print(np.min(v_phi), np.max(v_phi))
        # print(v_r,v_phi)
        v_x = v_r*np.cos(PHI) - v_phi*np.sin(PHI)    # dimensions: (nrad, nphi)
        v_y = v_r*np.sin(PHI) + v_phi*np.cos(PHI)
        # print(np.min(v_x), np.max(v_x))
        # print(np.min(v_y), np.max(v_y))
        v_z = np.zeros((nrad, nphi))
        vel = np.array((v_x,v_y,v_z))/(2*np.pi)     #*0.00021  # dimensions: (3, nrad, nphi)  convert from m/s to AU/yr ???
        # print(vel, "\n =================")
        for nr in np.arange(nrad):
            for nph in np.arange(nphi):
                v_t = vel[:,nr,nph]    # 3 x 1 velocity vector for 1 grain size at 1 timestep in 1 cell
                r_t = pos[:,nr,nph]
                e[i,nr,nph] = abs(calculate_e(r_t,v_t))
        print("ecc: ", np.min(e), np.max(e))
        con = ax.contourf(X, Y, e[i,:,:])
    fig.subplots_adjust(right=0.89, hspace=0.3)
    cbar_ax = fig.add_axes([0.91, 0.53, 0.02, 0.4])
    fig.colorbar(con, cax=cbar_ax, orientation="vertical", label="e")
    ax.legend(loc="upper right")
    fig.savefig(f"{plots_savedir}/{sim}_ecc.png")


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
    parser.add_argument('-style', metavar='style', type=str, nargs="*", default=["publication"] ,help="style sheet to apply to plots")

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
    style = args.style

    if plot_window:
        matplotlib.use('TkAgg')

    # ================== Read in data at timesteps =======================

    params_file = f'{simdir}/variables.par'
    params_dict = {}

    plotsizex = 3
    plotsizey = int(len(outputs)/plotsizex)+1

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

    cm = plt.get_cmap('gist_rainbow')
    colour_cycler = [cm(1.*i/5) for i in range(0,len(timesteps)+1)]

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
    if grog:
        v_dust = np.zeros((len(outputs), ndust, 2, nrad, nphi))   # additional dimension of 2 for x and y velocity
        for i,t in enumerate(outputs):
            # for n in np.arange(ndust):
            n=30
            dust_file_x = f"dustvx{n}_{int(t)}.dat"
            dust_file_y = f"dustvy{n}_{int(t)}.dat"
            v_dust[i,n,0] = np.fromfile(simdir+dust_file_x).reshape(nrad,nphi)   # vx
            v_dust[i,n,1] = np.fromfile(simdir+dust_file_y).reshape(nrad,nphi)   # vy

    if grog:
        uf = 10                                               # fragmentation velocity
        hr = hr0*(radii**f)                                   # aspect ratio
        cs = hr*(((2e30)*(6.67e-11))/(radii*1.5e11))**0.5     # [m/s]
        b = (uf**2/alpha)*(hr**-2)*(radii*1.5e11)/((2e30)*(6.67e-11)) # dimensionless
        p = (sigma_gas_1D*(cs**2)/((2*np.pi)**0.5))*(hr**-1)*((radii*1.5e11)**-1)
        pad = np.empty((len(timesteps), 1))*np.nan
        gamma = (radii/p)*np.abs(np.append(np.diff(p)/np.diff(radii), pad, axis=1))
        C = 2/(np.pi*hr)
        a_St1 = (2/np.pi)*(sigma_gas_1D/rhodust)              # plot St=1 line

        # size of largest grains in a fragmentation-dominated distribution
        # a_frag = 100*(2/(3*np.pi))*((uf**2)/(rhodust*1000*alpha))*(hr**-2)*(sigma_gas_1D*10/((2e30)*(6.67e-11)))*(radii*1.5e11)  # from Birnstiel+2012
        a_frag = (sigma_gas_1D/rhodust)*(3-(9-4*(b**2))**0.5)/(np.pi*b)
        # print(np.min(np.sqrt(9-4*(b**2))),np.max(9-4*(b**2)))

        # size of largest grains in a drift-dominated distribution
        # a_drift = 100*(2/(rhodust*1000*np.pi))*(hr**-2)*np.sum(sigma_dust_1D, axis=1)*10*(2/3)   # from Birnstiel+2012 (assume gamma=3/2)
        a_drift = (2/np.pi)*(np.sum(sigma_dust_1D, axis=1)/(rhodust*gamma*hr**2))    
    
    # ======================== Generate Plots ==========================

    if not plot_window:
        for s in style:
            plt.style.use([f"../styles/{s}.mplstyle"])

    print(f"-------------------\nPlotting outputs {outputs} for {sim}\n=============")

    plot_gas_sigma()
    if grog:
        plot_dust_contours()
        # plot_dust_e()
        # plot_dust_sigma()
        # plot_dustgasratio()
        # plot_dust_size_distribution()
        # plot_dust_mass()
        plot_dust_mass_by_grain()
        # plot_sigma_at_r([rps[0]+5])

    if plot_window:
        plt.show()
