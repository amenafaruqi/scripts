import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.lines as mlines
import matplotlib.patches as patches
from matplotlib.transforms import Bbox
from mpl_toolkits.axes_grid1.inset_locator import TransformedBbox, BboxPatch, BboxConnector
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter, NullFormatter
import matplotlib.ticker as mticker
import argparse
import re
plt.style.use('default')


def add_x_zoom_box_with_connectors(fig, ax_main, ax_zoom, xlims, edgecolor='k', lw=1.5):
    """
    Draw a zoom box in x only on ax_main, and connect:
    - bottom-left of the box to top-left of ax_zoom
    - bottom-right of the box to top-right of ax_zoom
    """
    # Get y-limits (should be the same on both axes)
    ylims = ax_main.get_ylim()

    # Draw zoom box (vertical rectangle)
    rect = patches.Rectangle(
        (xlims[0], ylims[0]),
        xlims[1] - xlims[0],
        ylims[1] - ylims[0],
        linewidth=lw,
        edgecolor=edgecolor,
        facecolor='none',
        linestyle='dotted',
        transform=ax_main.transData,
        zorder=5
    )
    ax_main.add_patch(rect)

    # --- Get box corners in display coords ---
    bl_disp = ax_main.transData.transform((xlims[0], ylims[0]))  # bottom-left
    br_disp = ax_main.transData.transform((xlims[1], ylims[0]))  # bottom-right

    # Convert to figure coordinates
    bl_fig = fig.transFigure.inverted().transform(bl_disp)
    br_fig = fig.transFigure.inverted().transform(br_disp)

    # --- Get zoomed axes top corners in figure coords ---
    tl_disp = ax_zoom.transAxes.transform((0, 1))  # top-left of zoomed plot
    tr_disp = ax_zoom.transAxes.transform((1, 1))  # top-right

    tl_fig = fig.transFigure.inverted().transform(tl_disp)
    tr_fig = fig.transFigure.inverted().transform(tr_disp)

    # --- Draw connectors from main plot to zoomed plot ---
    connectors = [
        (bl_fig, tl_fig),  # bottom-left to top-left
        (br_fig, tr_fig)   # bottom-right to top-right
    ]

    for (x1, y1), (x2, y2) in connectors:
        line = plt.Line2D(
            [x1, x2], [y1, y2],
            transform=fig.transFigure,
            color=edgecolor,
            linewidth=lw,
            linestyle='dotted',
            zorder=6
        )
        fig.add_artist(line)


# ====================== Gas Sigma ========================

def overlay_gas_sigmas(fig, ax, radii, sigma_gas_1D, model_num=0):
    print("Plotting gas surface density....")

    legend_elements = []
    for i, t in enumerate(timesteps):
        if "stat" in simdir:
            color = "k"
            lb = "Stationary"
        else:
            color = "r"
            lb = "Migrating"
        ax[i].plot(radii, sigma_gas_1D[i], label=lb, color=color)

        if planets:
            for rp in rps[:,i]:
                ax[i].axvline(rp, linestyle='dashed', color=color)
        
        ax[i].set_xlabel("Radius (AU)")
        ax[i].set_title(f"{round(t, 3)} Myr")
        ax[i].set_xscale("log")
        ax[i].set_yscale("log")
        ax[i].set_xlim(np.min(radii), np.max(radii))    

    ax[0].set_ylabel("$\Sigma_{gas} (g/cm^{2})$")
    ax[0].legend()
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

            color = next(ax[subps[n]]._get_lines.prop_cycler)['color']
            ax[subps[n]].plot(radii, (dust_mass_tot_binned), label=f"{round(t, 3)} Myr", color=color, linestyle=linestyles[model_num])
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
            [ax[subps[n]].axvline(rp, linestyle='dashed', color=color) for n in range(len(n_size_decade))]
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

    if "stat" in sim:
        ax.axvline(rp, linestyle='dashed', color="k")

    ax.set_xlabel("R (AU)")
    ax.set_ylabel("$M_{{dust}}/M_\oplus$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(handles=legend_elements)
    ax.set_xlim(min(radii), max(radii))
    fig.tight_layout()


# ================= Plot final dust mass in 10-100cm bin =====================

# def final_large_dust(fig, ax, radii, dust_mass, model_num=0):
#     n_size_decades = int(np.log10(maxgsize) - np.log10(mingsize))   # assumes min and max g size are the same for both models!
#     size_decades = np.split(np.arange(ndust), n_size_decades)
#     dust_mass_largest_bin = np.sum(dust_mass[-1,size_decades[-1],:], axis=0) # sum within bin, dimensions: nrad
    
#     c = colours[model_num]
#     sim = simdirs[model_num]
#     mlabel = "migrating"
#     if "stat" in sim:
#         mlabel = "stationary"

#     x =  re.search("Mp\d+",simdirs[model_num])
#     planet_mass =  x.group()[2:]
#     ax.plot(radii, dust_mass_largest_bin, color=c, label=f"{planet_mass}$M_\oplus$ {mlabel}")
#     for rp in rps[:,-1]:
#         ax.axvline(rp, linestyle='dashed', color='k')

#     if model_num == len(simdirs)-1:
#         ax.legend(loc="lower right")
#         ax.set_xlabel("R (AU)")
#         ax.set_ylabel("$M_{{dust}} (M_\oplus)$")
#         ax.set_xscale("log")
#         ax.set_yscale("log")
#         ax.set_xlim(min(radii), max(radii))
#         fig.tight_layout()


# ================= Plot final dust contours =====================

def plot_dust_contours(fig, ax, radii, sigma_dust_1D, model_num=0):
    R, A = np.meshgrid(radii, a)
    # levels = np.linspace(-18, 6, 13) 
    # levels = np.linspace(-9, 1, 6)  
    levels = np.linspace(-7, 1, 9)  
                  
    print("Plotting dust size contour maps....")
    if "stat" in simdir:
        lb = "Stationary"
    else:
        lb = "Migrating"

    for i, t in enumerate(timesteps):
        # Full and zoomed axes
        ax_full = ax[s*2, i]
        ax_zoom = ax[(s*2)+1, i]

        con = ax_full.contourf(R, A, np.log10(sigma_dust_1D[i]), cmap="Greys", levels=levels)
        con_zoom = ax_zoom.contourf(R, A, np.log10(sigma_dust_1D[i]), cmap="Greys", levels=levels)
        
        # Set axis limits and scales
        ax_zoom.set_ylim(np.min(a), np.max(a))      
        ax_full.set_ylim(np.min(a), np.max(a))        
        ax_full.set_xlim(np.min(radii), np.max(radii))
        ax_full.set_xscale("log")
        ax_full.set_yscale("log")
        ax_zoom.set_xscale("log")
        ax_zoom.set_yscale("log")

        # Plot curves
        ax_full.plot(radii, a_St1[i], c='black', alpha=0.7, label="St=1")
        ax_full.plot(radii, a_drift[i], c='deepskyblue', alpha=0.7, label="$a_{{drift}}$")
        ax_full.plot(radii, a_frag[i], c='red', alpha=0.7, label="$a_{{frag}}$")
        ax_zoom.plot(radii, a_St1[i], c='black', alpha=0.7, label="St=1")
        ax_zoom.plot(radii, a_drift[i], c='deepskyblue', alpha=0.7, label="$a_{{drift}}$")
        ax_zoom.plot(radii, a_frag[i], c='red', alpha=0.7, label="$a_{{frag}}$")

        if planets:
            for rp in rps[:,i]:
                ax_full.axvline(rp, linestyle='dashed', color='k')
                ax_zoom.axvline(rp, linestyle='dashed', color='k')
                ax_zoom.set_xlim(rp-5, rp+20)
                ax_zoom.set_xticks(np.arange(round(rp-5,-1), round(rp+20,-1), 10))
                ax_full.set_xticks([10, 100])
                ax_full.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
                ax_full.xaxis.set_minor_formatter(NullFormatter())
                ax_zoom.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
                ax_zoom.xaxis.set_minor_formatter(NullFormatter())

        # Set labels and titles
        if i == 0:
            ax_full.set_ylabel("a (cm)")
            ax_zoom.set_ylabel("a (cm)")
        else:
            ax_full.yaxis.set_major_formatter(matplotlib.ticker.NullFormatter())  # removes major tick labels
            ax_full.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())  # removes minor tick labels (optional)
            ax_zoom.yaxis.set_major_formatter(matplotlib.ticker.NullFormatter())  # removes major tick labels
            ax_zoom.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())  # removes minor tick labels (optional)


        ax[3,i].set_xlabel("Radius (AU)")

        if s == 0:
            if i == int(len(o)/2):
                title = "$\\bf{Migrating}$" if "mig" in simdir else "$\\bf{Stationary}$"
                ax_full.set_title(f"{title}\n {round(t,3)} Myr")
            else:
                ax_full.set_title(f"{round(t,3)} Myr")
        else:
            if i == int(len(o)/2):
                title = "$\\bf{Migrating}$" if "mig" in simdir else "$\\bf{Stationary}$"
                ax_full.set_title(title, pad=10.0)

        # Annotate max grain size in densest region of ring
        rp = np.array(planet_rs)[s,0,i]    # Look exterior to planet

        xmin, xmax = rp, rp+15
        max_dens = int(np.max(np.log10(sigma_dust_1D[i])))

        level_index = list(con.levels).index(max_dens) - 1
        outer_level = con.levels[level_index - len(con.levels)]
        outer_paths = con.collections[level_index - len(con.levels)].get_paths()
        max_y = -np.inf
        max_y_point = None

        for path in outer_paths:
            vertices = path.vertices  # shape (N, 2), where each row is (x, y)
            x_vals = vertices[:, 0]
            y_vals = vertices[:, 1]

            # Filter for points within the x-range
            mask = (x_vals >= xmin) & (x_vals <= xmax)
            if np.any(mask):
                y_vals_in_range = y_vals[mask]
                x_vals_in_range = x_vals[mask]
                local_max_y = y_vals_in_range.max()
                if local_max_y > max_y:
                    max_y = local_max_y
                    max_index = mask.nonzero()[0][y_vals_in_range.argmax()]
                    max_y_point = vertices[max_index]

        # Annotate if a valid point was found
        if max_y_point is not None:
            ax_zoom.annotate(
                f"{round(max_y_point[1], 1):.1f}cm",
                xy=(max_y_point[0], max_y_point[1]),
                xytext=(4, 10),
                textcoords="offset points",
                color="darkgreen",
                fontsize=11,
                arrowprops=dict(arrowstyle="->", color="darkgreen", lw=0.8),
                ha='left'
            )


    fig.tight_layout()
    
    if s == 1:
        fig.subplots_adjust(right=0.89, hspace=0.3)
        top_ax = ax[0, len(timesteps)-1]
        bottom_ax = ax[3, len(timesteps)-1]
        top = top_ax.get_position().y1
        bottom = bottom_ax.get_position().y0
        height = top - bottom
        right = top_ax.get_position().x1

        cax = fig.add_axes([right + 0.01, bottom, 0.02, height])
        fig.colorbar(con, cax=cax, orientation="vertical", label="log$[\Sigma_{d}\ (g/cm^2)]$")
        ax[0,0].legend(loc="lower left")


        for i, t in enumerate(timesteps):
            ax_full0 = ax[0, i]
            ax_zoom0 = ax[1, i]
            ax_full1 = ax[2, i]
            ax_zoom1 = ax[3, i]
            # print(np.array(planet_rs).shape)
            rp0 = np.array(planet_rs)[0,0,i]    # sim 1
            rp1 = np.array(planet_rs)[1,0,i]    # sim 2
            xlims0 = (rp0 - 5, rp0 + 20)
            xlims1 = (rp1 - 5, rp1 + 20)
            add_x_zoom_box_with_connectors(fig, ax_full0, ax_zoom0, xlims0)
            add_x_zoom_box_with_connectors(fig, ax_full1, ax_zoom1, xlims1)



# ================= Plot dust-gas ratios =====================

def overlay_dust_gas_ratio(fig, ax, radii, dust_mass_tot, gas_mass, model_num=0):
    ax.set_prop_cycle(color=colour_cycler)
    print("Plotting dust-gas ratio....")

    c = colours[model_num]
    sim = simdirs[model_num]
    mlabel = "migrating"
    if "stat" in sim:
        mlabel = "stationary"

    x =  re.search("Mp\d+",simdirs[model_num])
    planet_mass =  x.group()[2:]

    for i,t in enumerate(timesteps):
        dustgasratio = dust_mass_tot[i]/gas_mass[i]
        ax.plot(radii, dustgasratio, label=f"{planet_mass}$M_\oplus$ {mlabel}", color=c)
    
        if planets:
            if mlabel != "stationary":
                for rp in rps[:,i]:
                    ax.axvline(rp, linestyle='dashed', color=c)
    
    if mlabel == "stationary":
        for rp in rps[:,i]:
            ax.axvline(rp, linestyle='dashed', color="k")


    if model_num == len(simdirs)-1:
        ax.legend()
        ax.set_xlim(np.min(radii), np.max(radii))
        ax.set_xlabel("R (AU)")
        ax.set_ylabel("dust-gas ratio")
        ax.set_xscale("log")
        ax.set_yscale("log")

        fig.tight_layout()


# =================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate 1D plots', prefix_chars='-')

    parser.add_argument('-simdirs', metavar='simdirs', type=str, nargs="*", default=[] ,help="simulation directories containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images/comparison_plots" ,help="directory to save plots to")
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
    planet_rs = []

    if plot_window:
        matplotlib.use('TkAgg')
    
    if not plot_window:
        for s in style:
            plt.style.use([f"../styles/{s}.mplstyle"])

    # =================== Define figures and axes ========================
    if "gsig" in plots:
        fig_gas_sigma, ax_gas_sigma = plt.subplots(ncols=len(o), figsize=(int(len(o)*3),4), sharey=True)    # pls do not input > 5 timesteps.
    if "dcon" in plots:
        fig_con, ax_con = plt.subplots(ncols=len(o), nrows=4, figsize=(int(len(o)*3),11))#, sharey=True, sharex=True)
    if "dm" in plots:
        fig_dust_mass, ax_dust_mass = plt.subplot_mosaic("AABBCC;DDDEEE;FFFGGG", figsize=(17,16))
    if "dmt" in plots:
        fig_dm_tot, ax_dm_tot = plt.subplots(figsize=(6,5))
    if "dmf" in plots:
        fig_dm_f, ax_dm_f = plt.subplots(figsize=(6,5))
    if "dgr" in plots:
        fig_dgr, ax_dgr = plt.subplots(ncols=len(o), figsize=(10,5), sharey=True)

    # ================== Read in data at timesteps =======================

    for s, simdir in enumerate(simdirs):
        print (f"============ Plotting data for {simdir} =============")
        params_file = f'{simdir}/variables.par'
        params_dict = {}

        plotsizex = 3
        plotsizey = int(len(o)/plotsizex)+1

        # if ("10_" in simdir) or ("20_" in simdir):
        #     outputs = o/5   # need to account for different timestepping in different models
        # if ("100_" in simdir):
        #     outputs = np.array(o)*2
        # else:
        #     outputs = np.array(o)
        # outputs = outputs.astype(int)
        outputs = np.array(o)
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
            planet_rs.append(rps)

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
        gas_mass = np.zeros((len(outputs), nrad))
        sigma_dust = np.zeros((len(outputs), ndust, nrad, nphi))
        dust_mass = np.zeros((len(outputs), ndust, nrad))

        for i,t in enumerate(outputs):
            gasfile = f"/gasdens{int(t)}.dat" 
            sigma_gas[i] = np.fromfile(simdir+gasfile).reshape(nrad,nphi)/(1.125e-7)         # convert back to g/cm2
            n_grains = np.zeros((ndust))
            if grog:
                for n in np.arange(ndust):
                    dust_file = f"/dustdens{n}_{int(t)}.dat"
                    sigma_dust[i,n] = np.fromfile(simdir+dust_file).reshape(nrad,nphi)/(1.125e-7)

        sigma_dust_azimsum = np.sum(sigma_dust, axis=3)                  # sum over all phi 
        sigma_gas_azimsum = np.sum(sigma_gas, axis=2)                    # sum over all phi   

        sigma_gas_1D = sigma_gas_azimsum/nphi                           # dimensions: (noutputs, nrad) 
        sigma_dust_1D = sigma_dust_azimsum/nphi                         # dimensions: (noutputs, ndust, nrad)  
        # sigma_dust_sum_1D = avgdustdens_azimsum/nphi                  # dimensions: (noutputs, nrad)
        # for i,t in enumerate(outputs):
        #     # if "Mp60" in simdir:
        #     #     t = t*5
        #     # dust mass for dust of size a as a function of r
        #     # dust_mass[i,:,:] = [2*np.pi*radii*sigma_dust_1D[i,n,:]*delta_r*333030 for n in range(ndust)]
        #     gas_mass[i,:] = 2*np.pi*radii*sigma_gas_1D[i,:]*delta_r*333030      # convert from Msun to Mearth


        if "dcon" in plots and grog:
            uf = 0.0021                                           # fragmentation velocity 10 m/s in AU/yr
            hr = hr0*(radii**f)                                   # aspect ratio
            b = (uf**2)*radii/(4*(np.pi**2)*alpha*(hr**2))

            cs = hr*2*np.pi/(radii**0.5)
            p = (sigma_gas_1D*(cs**2)/((2*np.pi)**0.5))*(hr**-1)*(radii**-1)
            pad = np.empty((len(timesteps), 1))*np.nan
            gamma = (radii/p)*np.abs(np.append(np.diff(p)/np.diff(radii), pad, axis=1))
            C = 2/(np.pi*hr)
            a_St1 = (2/np.pi)*(sigma_gas_1D/rhodust)              # plot St=1 line

            # size of largest grains in a fragmentation-dominated distribution
            a_frag = 2*sigma_gas_1D*b/(rhodust*3*np.pi)

            # size of largest grains in a drift-dominated distribution
            a_drift = (2/np.pi)*(np.sum(sigma_dust_1D, axis=1)/(rhodust*gamma*hr**2))    
    
        if "gsig" in plots:
            overlay_gas_sigmas(fig_gas_sigma, ax_gas_sigma, radii, sigma_gas_1D, s)
        if "dcon" in plots:
            plot_dust_contours(fig_con, ax_con, radii, sigma_dust_1D, s)
        if "dm" in plots:
            overlay_dust_mass(fig_dust_mass, ax_dust_mass, radii, dust_mass, s)
        if "dmt" in plots:
            overlay_total_dust_mass(fig_dm_tot, ax_dm_tot, radii, dust_mass_tot, s)
        if "dmf" in plots:
            final_large_dust(fig_dm_f, ax_dm_f, radii, dust_mass, s)
        if "dgr" in plots:
            overlay_dust_gas_ratio(fig_dgr, ax_dgr, radii, dust_mass_tot, gas_mass, s)
    
    # ======================== Generate Plots ==========================
    print(f"-------------------\nPlotting comparison plots for {simdirs}\n=============")

    sim = simdirs[0].split("/")[-1].split("_")[0]    # ignore full file path, just get planet mass
    if "gsig" in plots:
        fig_gas_sigma.savefig(f"{plots_savedir}/{sim}_comparison_gassigma.png")
    if "dcon" in plots:
        fig_con.savefig(f"{plots_savedir}/{sim}_comparison_contour.png")
    if "dm" in plots:
        fig_dust_mass.savefig(f"{plots_savedir}/{sim}_comparison_graindist.png")
    if "dmt" in plots:
        fig_dm_tot.savefig(f"{plots_savedir}/{sim}_comparison_Mdust.png")
    if "dmf" in plots:
        fig_dm_f.savefig(f"{plots_savedir}/comparison_largedust.png")
    if "dgr" in plots:
        fig_dgr.savefig(f"{plots_savedir}/{sim}_comparison_dgr.png")

    if plot_window:
        plt.show()
