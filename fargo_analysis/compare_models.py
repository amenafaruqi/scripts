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


# ================= Dust Mass by Grain Size ===================
def overlay_dust_mass(fig, ax, radii, dust_mass, model_num=0):
    print("Plotting dust distribution by grain size....")

    n_size_decades = int(np.log10(maxgsize) - np.log10(mingsize))   # assumes min and max g size are the same for both models!
    size_decades = np.split(np.arange(ndust), n_size_decades)

    dust_mass_tot_binned = np.zeros((len(outputs), n_size_decades, nrad))
    subps = "ABCDEFG"
    legend_elements = []
    for n, size_decade in enumerate(size_decades):
        ax[subps[n]].set_prop_cycle(color=colour_cycler)
        for i, t in enumerate(timesteps):
            dust_mass_tot_binned = np.sum(dust_mass[i,size_decade,:], axis=0)
            if t == 0:
                dust_mass_tot_binned0 = dust_mass_tot_binned.copy()

            color = next(ax[subps[n]]._get_lines.prop_cycler)['color']
            
            ax[subps[n]].plot(radii, (dust_mass_tot_binned/dust_mass_tot_binned0), label=f"{round(t, 3)} Myr", color=color, linestyle=linestyles[model_num])
            if n == len(size_decades)-1:
                legend_elements.append(Line2D([0], [0], color=color, lw=2, label=f"{round(t, 3)} Myr"))

            if planets:
                for rp in rps[:,i]:
                    ax[subps[n]].axvline(rp, linestyle='dashed', color=color)

        dustsizes = [(10**n) * mingsize, (10**(n+1)) * mingsize]
        dustsizes = [np.format_float_positional(d,3,fractional=False,unique=True) for d in dustsizes]
        ax[subps[n]].set_title(f"{dustsizes[0]}-{dustsizes[1]}cm")
        ax[subps[n]].set_xscale("log")
        ax[subps[n]].set_xticks([10,20,40,80])
        ax[subps[n]].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax[subps[n]].xaxis.set_minor_formatter(NullFormatter())
        # ax[subps[n]].ticklabel_format(axis='x',style="plain")
        # ax[subps[n]].set_yscale("log")
        # ax[subps[n]].set_ylim(0, 250)
        ax[subps[n]].set_xlim(10, 80)

    
    for m in range(model_num+1):
        sim = simdirs[m].split('models/')[-1]    # ignore full file path
        if "stat" in sim:
            mlabel = "stationary"
        else:
            mlabel = "migrating"
        legend_elements.append(Line2D([0], [0], color='k', linestyle=linestyles[m], label=f"{mlabel}"))

    ax["A"].set_ylabel(f"$M_{{dust}}/M_{{dust,0}}$")
    ax["D"].set_ylabel(f"$M_{{dust}}/M_{{dust,0}}$")
    ax["F"].set_ylabel(f"$M_{{dust}}/M_{{dust,0}}$")
    ax["A"].set_xlabel("R (AU)")
    ax["B"].set_xlabel("R (AU)")
    ax["C"].set_xlabel("R (AU)")
    ax["F"].set_xlabel("R (AU)")
    ax["G"].set_xlabel("R (AU)")
    # ax["D"].set_xticks([])
    # ax["E"].set_xticks([])
    ax["G"].legend(loc="upper right", handles=legend_elements)
    fig.tight_layout()


# =================== Total Dust Mass =====================
def overlay_total_dust_mass(fig, ax, radii, dust_mass_tot, model_num=0):
    ax.set_prop_cycle(color=colour_cycler)
    print("Plotting total dust mass....")

    legend_elements=[]
    for i, t in enumerate(timesteps):
        color = next(ax._get_lines.prop_cycler)['color']
        ax.plot(radii, dust_mass_tot[i], label=f"{round(t, 3)} Myr", color=color, linestyle=linestyles[model_num])
        legend_elements.append(Line2D([0], [0], color=color, lw=2, label=f"{round(t, 3)} Myr"))

        if planets:
            for rp in rps[:,i]:
                ax.axvline(rp, linestyle='dashed', color=color)

    for m in range(model_num+1):
        sim = simdirs[m].split('models/')[-1]    # ignore full file path
        legend_elements.append(Line2D([0], [0], color='k', linestyle=linestyles[m], label=f"{sim}"))

    ax.set_xlabel("R (AU)")
    ax.set_ylabel("$M_{{dust}}/M_\oplus$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(handles=legend_elements)
    ax.set_xlim(min(radii), max(radii))
    fig.tight_layout()


# ================= Plot final dust mass in 10-100cm bin =====================

def final_large_dust(fig, ax, radii, dust_mass, model_num=0):
    n_size_decades = int(np.log10(maxgsize) - np.log10(mingsize))   # assumes min and max g size are the same for both models!
    size_decades = np.split(np.arange(ndust), n_size_decades)
    dust_mass_largest_bin = np.sum(dust_mass[-1,size_decades[-1],:], axis=0) # sum within bin, dimensions: nrad
    
    c = colours[model_num]
    sim = simdirs[model_num]
    mlabel = "migrating"
    if "stat" in sim:
        mlabel = "stationary"

    x =  re.search("Mp\d+",simdirs[model_num])
    planet_mass =  x.group()[2:]
    ax.plot(radii, dust_mass_largest_bin, color=c, label=f"{planet_mass}$M_\oplus$ {mlabel}")
    for rp in rps[:,-1]:
        print(simdir, rp)
        ax.axvline(rp, linestyle='dashed', color='k')

    if model_num == len(simdirs)-1:
        ax.legend(loc="lower right")
        ax.set_xlabel("R (AU)")
        ax.set_ylabel("$M_{{dust}} (M_\oplus)$")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(min(radii), max(radii))
        fig.tight_layout()

# ================= Plot final dust contours =====================
def plot_final_contours(fig, radii, sigma_dust_1D, model_num=0):
    R, A = np.meshgrid(radii, a)
    levels = np.linspace(-18, 6, 13) 
    plotsizey = int(len(simdirs)/plotsizex)+1
                  
    print("Plotting dust size contour maps....")
    if model_num == 0:
        fig.suptitle(f"t={timesteps[-1]}Myr")
    ax = fig.add_subplot(plotsizey, plotsizex, model_num+1)
    sigmas = sigma_dust_1D[i]
    con = ax.contourf(R, A, np.log10(sigmas), cmap="Greys", levels=levels)
    ax.set_ylim(np.min(a), np.max(a))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.plot(radii, a_St1[i], c='black', alpha=0.7, label="St=1")
    ax.plot(radii, a_drift[i], c='deepskyblue', alpha=0.7, label="$a_{{drift}}$")
    ax.plot(radii, a_frag[i], c='red', alpha=0.7, label="$a_{{frag}}$")
    if planets:
        for rp in rps[:,i]:
            ax.axvline(rp, linestyle='dashed', color='black')
        
    if not model_num%plotsizex:
        ax.set_ylabel("a (cm)")
    else:
        ax.set_yticks([])
    if model_num < plotsizex and len(simdirs) > plotsizex:
        ax.set_xticks([])
    else:
        ax.set_xlabel("R (AU)")
    sim = simdir.split("models/")[-1]
    ax.set_title(sim)

    fig.tight_layout()
    
    if model_num == len(simdirs)-1:
        fig.subplots_adjust(right=0.89, hspace=0.3)
        # cbar_ax = fig0.add_axes([0.91, 0.53, 0.02, 0.4])
        cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
        fig.colorbar(con, cax=cax, orientation="vertical", label="log$[\Sigma (g/cm^{{2}})]$")
        ax.legend(loc="upper right")


# =================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate 1D plots', prefix_chars='-')

    parser.add_argument('-simdirs', metavar='simdirs', type=str, nargs="*", default=[] ,help="simulation directories containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images" ,help="directory to save plots to")
    parser.add_argument('-o', metavar='outputs',default=[], type=int, nargs="*" ,help="outputs to plot")
    parser.add_argument('-plots', metavar='plots',default=[], type=str, nargs="*" ,help="plots to produce")
    parser.add_argument('-noplanet', action="store_false")
    parser.add_argument('-nogrog', action="store_false")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-style', metavar='style', type=str, nargs="*", default=["publication"], help="style sheet to apply to plots")

    args = parser.parse_args()
    o = args.o
    plots = args.plots
    simdirs = args.simdirs
    planets = args.noplanet
    grog = args.nogrog
    plot_window = args.plot_window
    plots_savedir = args.savedir
    style = args.style

    cm = plt.get_cmap('gist_rainbow')
    colour_cycler = [cm(1.*i/5) for i in range(0,len(o)+1)]
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    linestyles = ['solid', 'dashdot', 'dotted']

    if plot_window:
        matplotlib.use('TkAgg')
    
    if not plot_window:
        for s in style:
            plt.style.use([f"../styles/{s}.mplstyle"])

    # =================== Define figures and axes ========================
    if "gs" in plots:
        fig_gas_sigma, ax_gas_sigma = plt.subplots(figsize=(6,5))
    if "dcon" in plots:
        fig_con = plt.figure(figsize=(17,16))
    if "dm" in plots:
        fig_dust_mass, ax_dust_mass = plt.subplot_mosaic("AABBCC;DDDEEE;FFFGGG", figsize=(17,16))
    if "dmt" in plots:
        fig_dm_tot, ax_dm_tot = plt.subplots(figsize=(6,5))
    if "dmf" in plots:
        fig_dm_f, ax_dm_f = plt.subplots(figsize=(6,5))
    if "st" in plots:
        fig_st, ax_st = plt.subplots(nrows=3, ncols=4, figsize=(17,16))

    # ================== Read in data at timesteps =======================

    for s, simdir in enumerate(simdirs):
        print (f"============ Plotting data for {simdir} =============")
        params_file = f'{simdir}/variables.par'
        params_dict = {}

        plotsizex = 3
        plotsizey = int(len(o)/plotsizex)+1

        if (("Mp120" in simdir or "Mp160" in simdir) and "lowres" in simdir) or ("Mp25" in simdir and "test" in simdir):
            outputs = np.array(o)/5
            print("blah1")
            print(outputs)
        elif "Mp160" in simdir and "final" in simdir:
            outputs = np.array(o)*2
            print("blah2")
        else:
            outputs = np.array(o)
        outputs = outputs.astype(int)
        print(outputs)

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
                # xp, yp = planet_data[np.array(outputs)*5][:,1], planet_data[np.array(outputs)*5][:,2]
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
            gasfile = f"/gasdens{int(t)}.dat" 
            sigma_gas[i] = np.fromfile(simdir+gasfile).reshape(nrad,nphi)/(1.125e-7)         # convert back to g/cm2
            n_grains = np.zeros((ndust))
            if grog:
                for n in np.arange(ndust):
                    dust_file = f"/dustdens{n}_{int(t)}.dat"
                    sigma_dust[i,n] = np.fromfile(simdir+dust_file).reshape(nrad,nphi)/(1.125e-7)
            # if s == 1:
            #     for g in range(len(sav_sigma[0,:])):
        sigma_dust_azimsum = np.sum(sigma_dust, axis=3)                  # sum over all phi 
        sigma_gas_azimsum = np.sum(sigma_gas, axis=2)                    # sum over all phi   

        sigma_gas_1D = sigma_gas_azimsum/nphi                           # dimensions: (noutputs, nrad) 
        sigma_dust_1D = sigma_dust_azimsum/nphi                         # dimensions: (noutputs, ndust, nrad)  
        # sigma_dust_sum_1D = avgdustdens_azimsum/nphi                  # dimensions: (noutputs, nrad)
        
        for i,t in enumerate(outputs):
            if "Mp60" in simdir:
                t = t*5
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
    
        if "gs" in plots:
            overlay_gas_sigmas(fig_gas_sigma, ax_gas_sigma, radii, sigma_gas_1D, s)
        if "dcon" in plots:
            plot_final_contours(fig_con, radii, sigma_dust_1D, s)
        if "dm" in plots:
            overlay_dust_mass(fig_dust_mass, ax_dust_mass, radii, dust_mass, s)
        if "dmt" in plots:
            overlay_total_dust_mass(fig_dm_tot, ax_dm_tot, radii, dust_mass_tot, s)
        if "dmf" in plots:
            final_large_dust(fig_dm_f, ax_dm_f, radii, dust_mass, s)
        if "st" in plots:
            plot_st(fig_st, ax_st, radii, a, s)
    
    # ======================== Generate Plots ==========================
    print(f"-------------------\nPlotting comparison plots for {simdirs}\n=============")

    sim = simdirs[0].split("simulations/")[-1].split("_")[0]    # ignore full file path
    print(sim)
    if "gs" in plots:
        fig_gas_sigma.savefig(f"{plots_savedir}/{sim}_comparison_gassigma.png")
    if "dcon" in plots:
        fig_con.savefig(f"{plots_savedir}/comparison_contour.png")
    if "dm" in plots:
        fig_dust_mass.savefig(f"{plots_savedir}/{sim}_comparison_graindist.png")
    if "dmt" in plots:
        fig_dm_tot.savefig(f"{plots_savedir}/{sim}_comparison_Mdust.png")
    if "dmf" in plots:
        fig_dm_f.savefig(f"{plots_savedir}/comparison_largedust.png")
    if "st" in plots:
        fig_st.savefig(f"{plots_savedir}/{sim}_comparison_St.png")

    if plot_window:
        plt.show()
