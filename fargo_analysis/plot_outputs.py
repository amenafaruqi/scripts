import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
plt.style.use('default')

import warnings
warnings.filterwarnings("ignore")

# -------------- PLOTTING FUNCTIONS ----------------
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
    # cbar_ax = fig0.add_axes([0.91, 0.53, 0.02, 0.4])
    cax = fig0.add_axes([ax0.get_position().x1+0.01,ax0.get_position().y0,0.02,ax0.get_position().height])
    fig0.colorbar(con, cax=cax, orientation="vertical", label="log$[\Sigma (g/cm^{{2}})]$")
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
    fig = plt.figure(figsize=(17,16))
    print("Plotting dust surface density by grain size....")

    n_size_decades = int(np.log10(maxgsize) - np.log10(mingsize))  # e.g. 7
    size_decades = np.split(np.arange(ndust), n_size_decades)  # e.g. 70/7 = 10
    sigma_dust_binned = np.zeros((len(outputs), n_size_decades, nrad))
    psizex = 3
    psizey = int(n_size_decades/psizex)+1

    # Plot sigma dust for each size decade
    for n, size_decade in enumerate(size_decades):
        ax = fig.add_subplot(psizey, psizex, n+1)
        ax.set_prop_cycle(color=colour_cycler)

        for i, t in enumerate(timesteps):
            sigma_dust_binned = np.sum(sigma_dust_1D[i,size_decade,:], axis=0)
            color = next(ax._get_lines.prop_cycler)['color']

            if not p_orbits:
                tlabel = f"{round(t, 3)} Myr"
            elif planets:
                tlabel = f"{int(round(planet_orbits[i], 0))} orbits"
            
            ax.plot(radii, sigma_dust_binned, label=tlabel, color=color)

            if planets:
                for rp in rps[:,i]:
                    ax.axvline(rp, linestyle='dashed', color=color)
            
        if not n%psizex:
            ax.set_ylabel("$\Sigma_{dust} (g/cm^{2})$")
        if n < psizex and n_size_decades > psizex:
            ax.set_xticks([])
        else:
            ax.set_xlabel("R (AU)")

        dustsizes = [(10**n) * mingsize, (10**(n+1)) * mingsize]
        dustsizes = [np.format_float_positional(d,3,fractional=False,unique=True) for d in dustsizes]
        ax.set_title(f"{dustsizes[0]}-{dustsizes[1]}cm")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(min(radii), max(radii))

    # Plot sigma dust for all grains summed
    ax = fig.add_subplot(psizey, psizex, n+2)
    ax.set_prop_cycle(color=colour_cycler)
    for i, t in enumerate(timesteps):
        if not p_orbits:
            tlabel = f"{round(t, 3)} Myr"
        elif planets:
            tlabel = f"{int(round(planet_orbits[i], 0))} orbits"

        color = next(ax._get_lines.prop_cycler)['color']
        ax.plot(radii, sigma_dust_tot[i], label=tlabel, color=color)

        if planets:
            for rp in rps[:,i]:
                ax.axvline(rp, linestyle='dashed', color=color)

    ax.set_xlabel("R (AU)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_title("All Dust")
    ax.set_xlim(np.min(radii), np.max(radii))

    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_sigmadust.png")


# ================= Total dust mass evolution over time  ===============

def plot_total_dust_mass():
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
    ax.set_ylabel("$M_{{dust}}/M_\oplus$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlim(min(radii), max(radii))
    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_Mdust.png")


# ============ Dust mass evolution in individual bins =============

def plot_dust_mass():
    fig = plt.figure(figsize=(17,16))
    print("Plotting dust mass by grain size....")

    n_size_decades = int(np.log10(maxgsize) - np.log10(mingsize))  # e.g. 7
    size_decades = np.split(np.arange(ndust), n_size_decades)  # e.g. 70/7 = 10
    dust_mass_tot_binned = np.zeros((len(outputs), n_size_decades, nrad))
    psizex = 3
    psizey = int(n_size_decades/psizex)+1

    for n, size_decade in enumerate(size_decades):
        ax = fig.add_subplot(psizey, psizex, n+1)
        ax.set_prop_cycle(color=colour_cycler)

        for i, t in enumerate(timesteps):
            dust_mass_tot_binned = np.sum(dust_mass[i,size_decade,:], axis=0)
            if t == 0:
                dust_mass_tot_binned0 = dust_mass_tot_binned.copy()

            # dust_mass_tot_binned = np.sum(dust_mass[i,size_decade[0]:size_decade[-1],:], axis=0)
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
            ax.set_ylabel(f"log ($M_{{dust}}/M_\oplus$)")
        if n < psizex and n_size_decades > psizex:
            ax.set_xticks([])
        else:
            ax.set_xlabel("R (AU)")

        dustsizes = [(10**n) * mingsize, (10**(n+1)) * mingsize]
        dustsizes = [np.format_float_positional(d,3,fractional=False,unique=True) for d in dustsizes]
        ax.set_title(f"{dustsizes[0]}-{dustsizes[1]}cm")
        ax.set_xscale("log")
        # ax.set_yscale("log")
        ax.set_xlim(min(radii), max(radii))
    
    ax.legend()

    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_dustmass.png")


# ============ Dust mass evolution over time per size decade =============

def plot_dust_mass_per_bin(bin0=-12,bin1=-1):
    print("Plotting dust mass by size bin....")
    fig, ax = plt.subplots(nrows=3, ncols=4, figsize=(17,16))
    size_bins = np.arange(bin0,bin1+1)

    for s,size_bin in enumerate(size_bins):
        ax = ax.flatten()
        ax[s].set_prop_cycle(color=colour_cycler)
        a_bin = a[size_bin]   # physical size of grain in chosen bin

        for i, t in enumerate(timesteps):
            dust_mass_bin = dust_mass[i,size_bin,:]
            if t == 0:
                dust_mass_bin0 = dust_mass_bin.copy()

            color = next(ax[s]._get_lines.prop_cycler)['color']
            tlabel = f"{round(t, 3)} Myr"
            ax[s].plot(radii, np.log10(dust_mass_bin), label=tlabel, color=color)

            if planets:
                for rp in rps[:,i]:
                    ax[s].axvline(rp, linestyle='dashed', color=color)

        if not s%4:
            ax[s].set_ylabel(f"log ($M_{{dust}}/M_\oplus$)")
    
        if s < 4:
            ax[s].set_xticks([])
        else:
            ax[s].set_xlabel("R (AU)")

        ax[s].set_title(f"{round(a_bin,2)} cm")
        ax[s].set_xscale("log")
        ax[s].legend()
        ax[s].set_xlim(min(radii), max(radii))

    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_dustmassperbin.png")


# ============ n(a) vs a  at each timestep ===============

def plot_dust_size_distribution():
    fig, ax = plt.subplots(1, dpi=150)
    ax.set_prop_cycle(color=colour_cycler)
    print("Plotting dust size distribution....")

    for i, t in enumerate(timesteps):
        if not p_orbits:
            tlabel = f"{round(t, 3)} Myr"
        elif planets:
            tlabel = f"{int(round(planet_orbits[i], 0))} orbits"

        color = next(ax._get_lines.prop_cycler)['color']

        # sum dust masses over all radii at a chosen timestep
        m_total = np.sum(dust_mass[i], axis=1)   # dimensions: len(a)

        # mass of a single dust grain of size a
        m_a = (4*np.pi/3)*rhodust*(a**3)

        # number of grains = total mass/mass of one grain
        n_dust = m_total/m_a

        ax.plot(a, n_dust, label=tlabel, color=color)

    # Plot MRN size distribution (Mathis et al. 1977)
    n_mrn = a**(-3.5)
    ax.plot(a, n_mrn, color='k', label="MRN")

    ax.set_xlabel("a (cm)")
    ax.set_ylabel("n(a)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_xlim(min(a), max(a))
    fig.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_dustsizes.png")


# ============== Plot dust eccentricities ==================
def plot_ecc():
    print("Plotting dust and gas eccentricities....")
    # Initialise figures for eccentrities
    fig = plt.figure(figsize=(26,10))
    fig_hist = plt.figure(figsize=(20,10))

    R, PHI = np.meshgrid(radii, phis, indexing="ij")   # dimensions: (nrad, nphi)
    X = R*np.cos(PHI)
    Y = R*np.sin(PHI)
    n_range = np.arange(ndust)
    e_dust = np.zeros((len(outputs),len(n_range), nrad,nphi))
    e_gas = np.zeros((len(outputs),nrad,nphi))

    for i, o in enumerate(outputs):
        ax = fig.add_subplot(len(outputs),1+ndust, ((1+ndust)*(i+1))-ndust)   # 1 row per output, 1 column for gas + 1 column per dust size
        ax_hist = fig_hist.add_subplot(len(outputs),1+ndust, ((1+ndust)*(i+1))-ndust)
        # Apply correction to v_phi when in the guiding-centre ref frame
        # Explained here: https://fargo3d.bitbucket.io/def_setups.html
        v_phi = v_gas[i,0,:,:]+(omegaframe*R)  # dimensions: (nrad, nphi)
        v_r = v_gas[i,1,:,:]    # dimensions: (nrad, nphi)
        
        # All scalars below have dimensions of (nrad, nphi)
        r_scalar = R
        v_scalar = ((v_r**2)+(v_phi**2))**0.5
        h_scalar = R*v_phi

        # Iterate through cells to calculate gas ecc in each cell
        for nr in np.arange(nrad):
            for nph in np.arange(nphi):
                r = np.abs(r_scalar[nr,nph])
                v = np.abs(v_scalar[nr,nph])
                h = np.abs(h_scalar[nr,nph])
                e_scalar = (1+(h**2*((v**2) - (2/r))))**0.5
                e_gas[i,nr,nph] = e_scalar
        
        # Plot histogram of eccentricities across entire grid
        ax_hist.hist(
            e_gas[i].flatten(),
            edgecolor="black",)
        ax_hist.set_title("Gas")
        ax_hist.set_ylabel("Frequency")
        ax_hist.set_xlabel("e")

        # Plot contour of gas eccentricities
        con_g = ax.contourf(X, Y, e_gas[i,:,:])    
        # cbar_ax = fig.add_axes([0.91, 0.53, 0.02, 0.4])
        fig.colorbar(con_g, orientation="vertical", label="e")
        ax.set_title("Gas")

        # Eccentricities for all dust species
        for ni,n in enumerate(n_range):
            # Iterate through dust species
            ax = fig.add_subplot(len(outputs),1+ndust, ((1+ndust)*(i+1))-ndust+ni+1)
            ax_hist = fig_hist.add_subplot(len(outputs),1+ndust, ((1+ndust)*(i+1))-ndust+ni+1)
            v_phi = v_dust[i,n,0,:,:]+(omegaframe*R)   # dimensions: (nrad, nphi)
            v_r = v_dust[i,n,1,:,:]    # dimensions: (nrad, nphi)

            # all scalars below have dimensions of (nrad, nphi)
            r_scalar = R
            v_scalar = ((v_r**2)+(v_phi**2))**0.5
            h_scalar = R*v_phi

            # Dust eccentricity calculation
            for nr in np.arange(nrad):
                for nph in np.arange(nphi):
                    r = np.abs(r_scalar[nr,nph])
                    v = np.abs(v_scalar[nr,nph])
                    h = np.abs(h_scalar[nr,nph])
                    e_scalar = (1+(h**2*((v**2) - (2/r))))**0.5
                    e_dust[i,ni,nr,nph] = e_scalar
            
            ax_hist.hist(
                e_dust[i,ni].flatten(),
                edgecolor="black")
            con_d = ax.contourf(X, Y, e_dust[i,ni,:,:])    
            
            st = st_max/(10**ni)
            ax.set_title(f"St={st}")
            ax_hist.set_title(f"St={st}")
        
            # cbar_ax = fig.add_axes([0.91, 0.53, 0.02, 0.4])
            fig.colorbar(con_d, orientation="vertical", label="e")

    # fig.subplots_adjust(right=0.89, hspace=0.3)
    fig.tight_layout()
    fig_hist.tight_layout()
    fig.savefig(f"{plots_savedir}/{sim}_ecc.png")
    fig_hist.savefig(f"{plots_savedir}/{sim}_ecchist.png")


# ==========================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate 1D plots', prefix_chars='-')

    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/astro/phrkvg/simulations"],help="working directory containing simulations")
    parser.add_argument('-sim', metavar='sim', type=str, nargs=1, default=["planettest"] ,help="simulation directory containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images" ,help="directory to save plots to")
    parser.add_argument('-o', metavar='outputs',default=[], type=int, nargs="*" ,help="outputs to plot")
    parser.add_argument('-plots', metavar='plots',default=["gsig", "dcon", "dmass"], type=str, nargs="*" ,help="plots to produce")
    parser.add_argument('-noplanet', action="store_false")
    parser.add_argument('-nogrog', action="store_false")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-porbits', action="store_true")
    parser.add_argument('-style', metavar='style', type=str, nargs="*", default=["publication"] ,help="style sheet to apply to plots")

    args = parser.parse_args()
    outputs = args.o
    plots = args.plots
    wd = args.wd[0]
    sim = args.sim[0]
    simdir = f"{wd}/{sim}/"
    planets = args.noplanet
    grog = args.nogrog
    plot_window = args.plot_window
    plots_savedir = args.savedir
    p_orbits = args.porbits
    style = args.style

    cm = plt.get_cmap('gist_rainbow')
    colour_cycler = [cm(1.*i/5) for i in range(0,len(outputs)+1)]

    if plot_window:
        matplotlib.use('TkAgg')

    # ================== Read in data at timesteps =======================

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
    omegaframe = float(params_dict['OMEGAFRAME'])
    if grog:
        mingsize = float(params_dict['MIN_GRAIN_SIZE'])
        maxgsize = float(params_dict['MAX_GRAIN_SIZE'])
        nss_coag = int(params_dict['NUMSUBSTEPS_COAG'])
        rhodust = float(params_dict['RHO_DUST'])
    else:
        st_max = float(params_dict['STOKES'])

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
    if spacing == "Log":        
        radii = np.array([np.exp((np.log(r_cells[n])+np.log(r_cells[n+1]))/2) for n in range(len(r_cells)-1)])
        delta_log_r = np.log(radii[1]) - np.log(radii[0])
        delta_r = radii*delta_log_r
    else:     # Lin or default grid
        radii = np.array([(r_cells[n]+r_cells[n+1])/2 for n in range(len(r_cells)-1)])
        delta_r = radii[1]-radii[0]

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
    sigma_dust_tot = np.sum(sigma_dust_1D, axis=1)                     # dimensions: (noutputs, nrad)   
    # sigma_dust_sum_1D = avgdustdens_azimsum/nphi                  # dimensions: (noutputs, nrad)
    
    for i,t in enumerate(outputs):
        # dust mass for dust of size a as a function of r
        dust_mass[i,:,:] = [2*np.pi*radii*sigma_dust_1D[i,n,:]*delta_r*333030 for n in range(ndust)]
        gas_mass[i,:] = 2*np.pi*radii*sigma_gas_1D[i,:]*delta_r*333030      # convert from Msun to Mearth

    dust_mass_tot = np.sum(dust_mass, axis=1)

    # Get dust and gas velocities
    if "ecc" in plots:
        v_dust = np.zeros((len(outputs), ndust, 2, nrad, nphi))   # additional dimension of 2 for x and y velocity
        v_gas = np.zeros((len(outputs), 2, nrad, nphi))           # additional dimension of 2 for x and y velocity
        for i,t in enumerate(outputs):
            for n in np.arange(ndust):
                dust_file_x = f"dustvx{n}_{int(t)}.dat"
                dust_file_y = f"dustvy{n}_{int(t)}.dat"
                gas_file_x = f"gasvx{int(t)}.dat"
                gas_file_y = f"gasvy{int(t)}.dat"
                v_dust[i,n,0] = np.fromfile(simdir+dust_file_x).reshape(nrad,nphi)   # vx (azimuthal v)
                v_dust[i,n,1] = np.fromfile(simdir+dust_file_y).reshape(nrad,nphi)   # vy (radial v)
                v_gas[i,0] = np.fromfile(simdir+gas_file_x).reshape(nrad,nphi)       # vx
                v_gas[i,1] = np.fromfile(simdir+gas_file_y).reshape(nrad,nphi)       # vy

    if "dcon" in plots and grog:
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

    if "gsig" in plots:
        plot_gas_sigma()
    if "ecc" in plots:
        plot_ecc()

    if grog:
        if "dcon" in plots:
            plot_dust_contours()
        if "dgr" in plots:
            plot_dustgasratio()
        if "dsizes" in plots:
            plot_dust_size_distribution()
        if "dmasstot" in plots:
            plot_total_dust_mass()
        if "dmass" in plots:
            plot_dust_mass()
        if "dmassbin" in plots:
            plot_dust_mass_per_bin()
        if "dsig" in plots:
            plot_dust_sigma()

    if plot_window:
        plt.show()
