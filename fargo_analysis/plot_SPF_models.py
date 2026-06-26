import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import SymLogNorm
import argparse
import re
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import NullFormatter

plt.style.use('default')

import warnings
warnings.filterwarnings("ignore")

# ====================== Gas Sigma ========================

def overlay_gas_sigmas(fig, ax, radii, sigma_gas_1D):
    print("Plotting gas surface density....")
    if "stat" in sim:
        ax0 = ax[0]
        if len(rps) == 1:
            label = "S1"
            colour = "k"
        else:
            label = "S2"
            colour = "r"
    else:
        ax0 = ax[1]
        if len(rps) == 1:
            label = "M1"
            colour = "k"
        else:
            label = "M2"
            colour = "r"

    # planet_masses = re.findall(r"Mp(\d+)_", sim)    # get planet masses from file path
    # if len(planet_masses) == 2:
    #     planet1_mass = planet_masses[1]             # outer planet mass, if present
    # planet0_mass = planet_masses[0]                 # inner planet mass, assume always present

    ax0.plot(radii, sigma_gas_1D, color=colour, label=label)
    for rp in rps:
        ax0.axvline(rp, linestyle='dashed', color=colour)
        
    ax[1].set_xlabel("Radius (AU)")
    ax0.set_ylabel("$\Sigma_{gas} (g/cm^{2})$")
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.set_xlim(np.min(radii), np.max(radii))
    ax0.legend()
    fig.tight_layout()


# ====================== Dust Sigma ========================
def overlay_dust_sigmas(fig, ax, radii, sigma_dust_1D):
    print("Plotting dust surface density....")
    # colours = ["red", "green", "blue"]
    titles = ["$< 10^{-3}$ cm", "$10^{-3} - 10^{-1}$ cm", "> $10^{-1}$ cm"]
    if "stat" in sim:
        coli = 0
        if len(rps) == 1:
            label = "S1"
            colour = "k"
        else:
            label = "S2"
            colour = "r"
    else:
        coli = 1
        if len(rps) == 1:
            label = "M1"
            colour = "k"
        else:
            label = "M2"
            colour = "r"

    sigma_dust_size0 = np.sum(sigma_dust_1D[0:20], axis=0)
    sigma_dust_size1 = np.sum(sigma_dust_1D[20:40], axis=0)
    sigma_dust_size2 = np.sum(sigma_dust_1D[40:70], axis=0)

    sigma_dusts = [sigma_dust_size0, sigma_dust_size1, sigma_dust_size2]

    for rowi in np.arange(3):
        ax[rowi, coli].plot(radii, sigma_dusts[rowi], color=colour, linestyle="dotted")
        ax[rowi, coli].set_xscale("log")
        ax[rowi, coli].set_yscale("log")
        if coli == 0:
            ax[rowi,coli].set_ylabel("$\Sigma_{dust} (g/cm^{2})$")
            if label ==  "S2":
                ax[rowi,coli].text(
                    0.94, 0.94, 
                    s=titles[rowi], 
                    transform=ax[rowi,coli].transAxes,
                    ha="right",
                    va="top",
                    fontsize="large", 
                    backgroundcolor="lightgrey",
                    )


        for rp in rps:
            ax[rowi,coli].axvline(rp, linestyle='dashed', color=colour)
    
        if len(rps)==2:
            ax[rowi,coli].set_xlim(rps[0]-3, rps[1]+45)

        # Plot ring boundaries
        planetmass0 = int(re.findall(r"(?<=Mp)\d+", sim)[0])
        r_hill0 = rps[0]*((planetmass0*1e-6)**(1/3))

        # Find radial cells closest to 2 and 20 R_Hill
        i_inner0 = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[0])+r_hill0*1)))
        i_outer0 = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[0])+r_hill0*20)))

        if len(rps) == 2:     # i.e. if 2-planet model, bound by outer planet location
            rp_p2 = rps[1]      # Outer planet location
            i_p2 = min(range(len(radii)), key=lambda i: abs(radii[i]-((rp_p2))))
            # Use outer planet location as outer bound if it lies within 20 R_Hill (of inner planet)
            i_outer0 = np.min([i_outer0, i_p2])

            # Mark ring edges for outer ring too
            planetmass1 = int(re.findall(r"(?<=Mp)\d+", sim)[1])          # Take inner planet mass only
            r_hill1 = rps[1]*((planetmass1*1e-6)**(1/3))

            # Find radial cells closest to 2 and 20 R_Hill
            i_inner1 = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[1])+r_hill1*1)))
            i_outer1 = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[1])+r_hill1*12)))
        
            # Plot thick line to indicate ring on plot
            ax[rowi,coli].plot(radii[i_inner1:i_outer1], sigma_dusts[rowi][i_inner1:i_outer1], color=colour, linewidth=2.5)
            
        ax[rowi,coli].plot(radii[i_inner0:i_outer0], sigma_dusts[rowi][i_inner0:i_outer0], color=colour, linewidth=2.5, label=label)


    if rowi == 2:    
        ax[rowi,coli].set_xlabel("Radius (AU)")

    ax[2,0].legend()
    ax[2,1].legend()
    ax[2,0].set_ylim(7e-17,100)

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

    # For testing ---------->
    # planetmass0 = int(re.findall(r"(?<=Mp)\d+", sim)[0])          # Take inner planet mass only
    # r_hill0 = rps[0]*((planetmass0*1e-6)**(1/3))

    # # Find radial cells closest to 2 and 20 R_Hill
    # i_inner0 = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[0])+r_hill0*2)))
    # i_outer0 = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[0])+r_hill0*20)))

    # if rowi:     # i.e. if 2-planet model, bound by outer planet location
    #     rp_p2 = rps[1]      # Outer planet location
    #     i_p2 = min(range(len(radii)), key=lambda i: abs(radii[i]-((rp_p2))))
    #     # Use outer planet location as outer bound if it lies within 20 R_Hill (of inner planet)
    #     i_outer0 = np.min([i_outer0, i_p2])

    #     # Mark ring edges for outer ring too
    #     planetmass1 = int(re.findall(r"(?<=Mp)\d+", sim)[1])          # Take inner planet mass only
    #     r_hill1 = rps[1]*((planetmass1*1e-6)**(1/3))

    #     # Find radial cells closest to 2 and 20 R_Hill
    #     i_inner1 = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[1])+r_hill1*2)))
    #     i_outer1 = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[1])+r_hill1*15)))

    #     ax0.axvline(radii[i_inner1], linestyle='dashed', color='red')
    #     ax0.axvline(radii[i_outer1], linestyle='dashed', color='red')

    # ax0.axvline(radii[i_inner0], linestyle='dashed', color='red')
    # ax0.axvline(radii[i_outer0], linestyle='dashed', color='red')
    # <------------

    sigmas = sigma_dust_1D
    con = ax0.contourf(R, A, np.log10(sigmas), cmap="magma", levels=levels)
    ax0.set_ylim(np.min(a), np.max(a))
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    if coli == 0:
        ax0.set_ylabel("a (cm)")
    if rowi == 1:
        ax0.set_xlabel("Radius (AU)")

    ax0.set_facecolor("k")
    ax0.plot(radii, a_St1, c='white', alpha=0.9, label="St=1", linewidth=1.5)
    ax0.plot(radii, a_drift, c='limegreen', alpha=0.9, label="$a_{{drift}}$", linewidth=1.5)
    ax0.plot(radii, a_frag, c='deepskyblue', alpha=0.9, label="$a_{{frag}}$", linewidth=1.5)
    ax0.grid(which="major", color="lightgrey", linewidth=0.5)

    for rp in rps:
        ax0.axvline(rp, linestyle='dashed', color='white', linewidth=1.5)

    ax0.set_title(label)
    # fig.suptitle(f"t={round(t,3)}Myr")
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


def calculate_ring_mass(radii, dust_mass, rps):
    # Set labels and indexes to populate arrays
    if "stat" in sim:
        if len(rps) == 1:
            rowi = 0
            coli = 0
        else:
            rowi = 1
            coli = 0
    else:
        if len(rps) == 1:
            rowi = 0
            coli = 1
        else:
            rowi = 1
            coli = 1
    
    planetmass = int(re.findall(r"Mp(\d+)_", sim)[0])          # Take inner planet mass only
    r_hill = rps[0]*((planetmass*1e-6)**(1/3))

    # Find radial cells closest to 2 and 20 R_Hill
    i_inner = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[0])+r_hill*2)))
    i_outer = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[0])+r_hill*20)))

    if rowi:     # i.e. if 2-planet model, bound by outer planet location
        rp_p2 = rps[1]      # Outer planet location
        i_p2 = min(range(len(radii)), key=lambda i: abs(radii[i]-((rp_p2))))
        # Use outer planet location as outer bound if it lies within 20 R_Hill (of inner planet)
        i_outer = np.min([i_outer, i_p2])

    # Dust mass in ring, broken down by dust bin
    # print(sim)
    # print(radii[i_inner], radii[i_outer])
    # print(dust_mass.shape)
    dust_mass_ring = np.sum(dust_mass[:,i_inner:i_outer], axis=1)       # dimensions = ndust = 70

    n_size_decades = int(np.log10(maxgsize) - np.log10(mingsize))   # assumes min and max g size are the same for both models!
    size_decades = np.split(np.arange(ndust), n_size_decades)

    for n, size_decade in enumerate(size_decades):
        ring_mass_by_size = np.sum(dust_mass_ring[size_decade])
        ring_masses[rowi, coli, n] = ring_mass_by_size

    # print(ring_masses) 
    return ring_masses

def calculate_outer_ring_mass(radii, dust_mass, rps):
    # Set labels and indexes to populate arrays
    if "stat" in sim:
        rowi = 2
        coli = 0
    else:
        rowi = 2
        coli = 1

    n_size_decades = int(np.log10(maxgsize) - np.log10(mingsize))
    size_decades = np.split(np.arange(ndust), n_size_decades)

    # dust_mass has dimensions of n_outputs x ndust x nrad
    # print(re.findall(r"Mp(\d+)", sim))
    planetmass = 50          # Take outer planet mass only
    r_hill = rps[1]*((planetmass*1e-6)**(1/3))

    # Find radial cells closest to 2 and 20 R_Hill
    i_inner = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[1])+r_hill*2)))
    i_outer = min(range(len(radii)), key=lambda i: abs(radii[i]-((1.07*rps[1])+r_hill*15)))
    dust_mass_ring = np.sum(dust_mass[:,i_inner:i_outer], axis=1)       # dimensions = (ndust)

    # Sum within each size decade to go from 70 to 7 ring masses
    for n, size_decade in enumerate(size_decades):
        ring_mass_by_size = np.sum(dust_mass_ring[size_decade])      # sum over all grain sizes within size decade to get single value
        ring_masses[rowi, coli, n] = ring_mass_by_size    # ring mass per size decade (for 1 model and timestep)

    # print(ring_masses) 
    return ring_masses


def plot_ring_mass(fig, ax, ring_masses, outer=False):
    # Set labels and indexes to populate arrays
    if "stat" in sim:
        coli = 0
        if len(rps) == 1:
            label = "S1"
            rowi = 0
        else:
            if not outer:
                label = "S2 \n (inner ring)"
                rowi = 1
            else:
                label = "S2 \n (outer ring)"
                rowi = 2
    else:
        coli = 1
        if len(rps) == 1:
            label = "M1"
            rowi = 0
        else:
            if not outer:
                rowi = 1
                label = "M2 \n (inner ring)"
            else:
                rowi = 2
                label = "M2 \n (outer ring)"


    size_labels = ["$10^{-5}-10^{-4}$ cm", "$10^{-4}-10^{-3}$ cm", "$10^{-3}-10^{-2}$ cm", "$10^{-2}-10^{-1}$ cm", "$10^{-1}$-1 cm", "1-10 cm", "10-100 cm"]
    ring_masses_sim = ring_masses[rowi, coli]    # dimensions = 7

    bottom=0

    total_mass = np.sum(ring_masses_sim)

    # TODO: Check outer ring bounds!!! 
    for i, rmass in enumerate(ring_masses_sim):
        if sim == sims[-1]:        
            p = ax.bar(label, rmass, 0.5, bottom=bottom, color=colour_cycler[i], label=size_labels[i])
        else:
            p = ax.bar(label, rmass, 0.5, bottom=bottom, color=colour_cycler[i])

        # --- Compute percentage ---
        if total_mass > 0:
            percent = (rmass / total_mass) * 100
        else:
            percent = 0

        # --- Add text in the middle of each segment ---
        y_pos = bottom + rmass / 2

        # Only label if the segment is large enough to see
        if percent > 4.5:
            ax.text(label, y_pos, f"{percent:.1f}%",
                    ha='center', va='center', fontsize=7, color='white')

        bottom += rmass
    
    # for i, rmass in enumerate(ring_masses_sim):
    #     if sim == sims[-1]:        
    #         p = ax.bar(label, rmass, 0.4, bottom=bottom, color=colour_cycler[i], label=size_labels[i])
    #     else:
    #         p = ax.bar(label, rmass, 0.4, bottom=bottom, color=colour_cycler[i])
    #     bottom += rmass

    if sim == sims[-1]:        
        cmap = mpl.cm.viridis_r
        bounds = np.arange(-5,3)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        box = ax.get_position()
        # ax.set_position([box.x0, box.y0, box.width, box.height*0.9])
        print(box)
        cbar_ax = fig.add_axes([box.x0, 0.89, box.width, box.height*0.05])

        cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                    cax=cbar_ax, orientation='horizontal',
                    label="Dust size (cm)", ticklocation="top")
        cbar.set_ticklabels(["$10^{-5}$", "$10^{-4}$", "$10^{-3}$", "$10^{-2}$", "$10^{-1}$", "$10^{0}$", "$10^{1}$", "$10^{2}$"])

    ax.set_ylabel("Ring dust mass ($M_\oplus$)")
    ax.set_xlabel("Model")
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0, box.width, box.height*0.96])
    # ax.legend(loc="upper center", ncol=3, bbox_to_anchor=(0.1, 0.83, 0.8, 0.5))

    # fig.tight_layout()



def plot_Macc(fig, ax, radii, a, sigma_dust_1D, v_dust):
    print("Plotting dust accretion rates....")

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

    R, A = np.meshgrid(radii, a)
    levels = np.linspace(-18, -10, 9)    
    # Azimuthally average radial velocities of all dust:              
    v_rad_avg = np.mean(v_dust[:,1,:,:], axis=2)     # Dimensions: ndust x nrad
    
    Macc = -2*np.pi*R*sigma_dust_1D*v_rad_avg*(1.125e-7)
    # Macc = v_rad_avg
    # print("Macc = \n", Macc)
    print("--------------------------")
    print(np.min(Macc), np.max(Macc))
    normed = SymLogNorm(linthresh=1e-13, linscale=0.0001, vmin=np.min(Macc)*0.7, vmax=np.max(Macc)*0.7)
    # normed = SymLogNorm(linthresh=1e-4, linscale=0.1, vmin=np.min(Macc)*0.8, vmax=np.max(Macc)*0.8)

    con = ax0.contourf(R, A, Macc, levels=10, norm=normed, cmap='coolwarm')
    ax0.plot(radii, a_St1, c='white', alpha=0.9, label="St=1", linewidth=1)

    # con = ax0.contourf(R, A, np.log10(Macc), cmap="Blues")
    ax0.set_ylim(np.min(a), 80)
    ax0.set_xlim(5, 140)
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    if coli == 0:
        ax0.set_ylabel("a (cm)")
    if rowi == 1:
        ax0.set_xlabel("Radius (AU)")

    for rp in rps:
        ax0.axvline(rp, linestyle='dashed', color='black')

    ax0.set_title(label)
    fig.tight_layout()

    if (rowi, coli) == (1,1):
        fig.subplots_adjust(bottom=0.1, hspace=0.1)
        cbar_ax = fig.add_axes([0.075, 0.07, 0.89, 0.03])
        fig.colorbar(con, cax=cbar_ax, orientation="horizontal", label="$\dot{M}$")
        # fig.colorbar(con, cax=cbar_ax, orientation="horizontal", label="$v_{rad}$")
        fig.tight_layout(rect=[0, 0.1, 0.98, 1])
        # ax0.legend()


def plot_SI_thresholds(fig, ax, radii, a, sigma_gas_1D, sigma_dust_1D):
    R, A = np.meshgrid(radii, a)
    levels = np.linspace(0, 1, 9)
    print("Plotting Z/Z_crit maps....")

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

    Z = sigma_dust_1D/sigma_gas_1D     # 70x635/635 = 70x635
    St = np.pi*rhodust*A/(sigma_gas_1D*2)    #70x635/635
    Omega_k = radii**(-3/2)            # Keplerian frequency 
    tau_s = St/Omega_k                 # Stopping time
    Z_crit = 10**((0.1*((np.log10(tau_s))**2)) + (0.07*np.log10(tau_s)) - 2.36)          # Eq. 15 from Lim et al. 2025a
    # Z_crit = 0.15*(np.log10(alpha)**2) - 0.24*np.log10(St)*np.log10(alpha) - 1.48*np.log10(St) + 1.18*np.log10(alpha)   # from Lim et al. 2026
    Z0 = Z/np.abs(Z_crit)

    con = ax0.contourf(R, A, Z0, cmap="GnBu", levels=levels, extend="max")
    # con2 = ax0.contour(R, A, np.log10(sigmas), linewidths=2, linestyles="dotted", levels=[-2])
    # ax0.set_facecolor("k")
    ax0.set_ylim(1e-4, 10)
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.set_xlim(6, 100)
    ax0.set_xticks([10, 40, 70, 100])


    if coli == 0:
        ax0.set_ylabel("a (cm)")
        ax0.set_xlim(25, 120)
        ax0.set_xticks([30, 60, 90, 120])
    if rowi == 1:
        ax0.set_xlabel("Radius (AU)")

    ax0.xaxis.set_minor_formatter(NullFormatter())
    ax0.xaxis.set_major_formatter(ScalarFormatter())
    ax0.ticklabel_format(style='plain', axis='x')

    for rp in rps:
        ax0.axvline(rp, linestyle='dashed', color='black')

    ax0.set_title(label)
    fig.tight_layout()

    if (rowi, coli) == (1,1):
        fig.subplots_adjust(bottom=0.3, hspace=0.3)
        cbar_ax = fig.add_axes([0.075, 0.17, 0.91, 0.03])
        fig.colorbar(con, cax=cbar_ax, orientation="horizontal", label="Z/Z$_{crit}$", extend="max")

        # fig.subplots_adjust(bottom=0.2, hspace=0.05)
        # cbar_ax = fig.add_axes([0.071, 0.05, 0.902, 0.02])
        # # cax = fig.add_subplot(ax[5, :])
        # fig.colorbar(con, cax=cbar_ax, orientation="horizontal", label="Z/Z$_{crit}$", extend="max")
        # fig.tight_layout(rect=[0, 0.07, 1, 1])



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
    o = args.o[-1]                    # Last timestep given
    plots = args.plots
    wd = args.wd[0]
    sims = args.sims
    plot_window = args.plot_window
    plots_savedir = args.savedir
    style = args.style

    cm = plt.get_cmap('viridis_r')
    colour_cycler = cm(np.linspace(0, 1, 7))       # For plotting different mass decades

    if plot_window:
        mpl.use('TkAgg')
    else:
        for s in style:
            plt.style.use([f"../styles/{s}.mplstyle"])

    # =================== Define figures and axes ========================
    if "gsig" in plots:
        fig_gas_sigma, ax_gas_sigma = plt.subplots(nrows=2, ncols=1, figsize=(6,9))
    if "dsig" in plots:
        fig_dust_sigma, ax_dust_sigma = plt.subplots(nrows=3, ncols=2, figsize=(9,10))
    if "dgr" in plots:
        fig_dgr, ax_dgr = plt.subplots(figsize=(6,4))
    if "growth" in plots:   # only specify for grog models
        fig_growth, ax_growth = plt.subplots(figsize=(7,7), nrows=2, ncols=2, sharex=True, sharey=True)
    if "dcon" in plots:   # only specify for grog models
        fig_con, ax_con = plt.subplots(figsize=(10,10), nrows=2, ncols=2, sharex=True, sharey=True)
    if "2dsig" in plots:   # only specify for grog models
        fig_2d, ax_2d = plt.subplots(figsize=(7,14), nrows=4, ncols=2, sharex=True, sharey=True)
    if "rmass" in plots:
        fig_m, ax_m = plt.subplots(figsize=(9,6))
    if "macc" in plots:
        fig_acc, ax_acc = plt.subplots(figsize=(8,8), nrows=2, ncols=2, sharex=True, sharey=True)
    if "si" in plots:   # only specify for grog models
        fig_si, ax_si = plt.subplots(figsize=(10,10), nrows=2, ncols=2, sharey=True)


    # ================== Read in data at timesteps =======================

    # Array to store ring masses
    ring_masses = np.zeros((3,2,7))           # 4 models x 7 size decades
    mps = []

    for s, sim in enumerate(sims):
        print (f"============ Plotting data for {sim} =============")
        params_file = f'{wd+sim}/variables.par'
        params_dict = {}
        output = o

        planet_masses = re.findall(r"(?<=Mp)\d+", sim)    # Planet masses, stored as strings in list
        print(planet_masses)

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

        # ====== Get gas and dust sigma from output files ======
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

        # Compute dust mass for dust of size a as a function of r
        dust_mass = 2*np.pi*radii*sigma_dust_1D*delta_r*333030*(1.125e-7)
        gas_mass = 2*np.pi*radii*sigma_gas_1D*delta_r*333030*(1.125e-7)      # convert from Msun to Mearth
        dust_mass_tot = np.sum(dust_mass, axis=0)                            # summed over all sizes

        # ===== Get velocity data, if needed =====
        if "macc" in plots:
            v_dust = np.zeros((ndust, 2, nrad, nphi))                  # additional dimension of 2 for x and y velocity
            # v_gas = np.zeros((len(outputs), 2, nrad, nphi))           # additional dimension of 2 for x and y velocity
            for n in np.arange(ndust):
                dust_file_x = f"/dustvx{n}_{output}.dat"
                dust_file_y = f"/dustvy{n}_{output}.dat"
                # gas_file_x = f"gasvx{int(t)}.dat"
                # gas_file_y = f"gasvy{int(t)}.dat"
                v_dust[n,0] = np.fromfile(wd+sim+dust_file_x).reshape(nrad,nphi)   # vx (azimuthal v)
                v_dust[n,1] = np.fromfile(wd+sim+dust_file_y).reshape(nrad,nphi)   # vy (radial v)
                # v_gas[i,0] = np.fromfile(simdir+gas_file_x).reshape(nrad,nphi)       # vx
                # v_gas[i,1] = np.fromfile(simdir+gas_file_y).reshape(nrad,nphi)       # vy

        # =========================================================================

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
        if "dsig" in plots:
            overlay_dust_sigmas(fig_dust_sigma, ax_dust_sigma, radii, sigma_dust_1D)
        if "dgr" in plots:
            overlay_dust_gas_ratio(fig_dgr, ax_dgr, radii, dust_mass_tot, gas_mass)
        if "growth" in plots:
            plot_tgrowth(fig_growth, ax_growth, radii, a, sigma_dust_1D)
        if "dcon" in plots:
            plot_dust_contours(fig_con, ax_con, radii, a, sigma_dust_1D)
        if "2dsig" in plots:
            plot_2D_sigma(fig_2d, ax_2d, radii, sigma_gas, sigma_dust)
        if "rmass" in plots:
            ring_masses = calculate_ring_mass(radii, dust_mass, rps)       
            plot_ring_mass(fig_m, ax_m, ring_masses)
            if len(planet_masses) > 1:
                ring_masses = calculate_outer_ring_mass(radii, dust_mass, rps)      
                plot_ring_mass(fig_m, ax_m, ring_masses, outer=True)
        if "macc" in plots:
            plot_Macc(fig_acc, ax_acc, radii, a, sigma_dust_1D, v_dust)
        if "si" in plots:
            plot_SI_thresholds(fig_si, ax_si, radii, a, sigma_gas_1D, sigma_dust_1D)


    # ======================== Generate Plots ==========================
    print(f"-------------------\nPlotting comparison plots for {sims}\n=============")

    if "gsig" in plots:
        fig_gas_sigma.savefig(f"{plots_savedir}/SPF_sigmagas.png")
    if "dsig" in plots:
        fig_dust_sigma.savefig(f"{plots_savedir}/SPF_sigmadust.png")
    if "dgr" in plots:
        fig_dgr.savefig(f"{plots_savedir}/SPF_dgr.png")
    if "growth" in plots:
        fig_growth.savefig(f"{plots_savedir}/SPF_growth.png")
    if "dcon" in plots:
        fig_con.savefig(f"{plots_savedir}/SPF_contours.png")
    if "2dsig" in plots:
        fig_2d.savefig(f"{plots_savedir}/SPF_2Dsigma.png")
    if "rmass" in plots:
        fig_m.savefig(f"{plots_savedir}/SPF_ring_masses.png")
    if "macc" in plots:
        fig_acc.savefig(f"{plots_savedir}/SPF_Macc.png")
    if "si" in plots:
        fig_si.savefig(f"{plots_savedir}/SPF_SI_thresholds.png")


    if plot_window:
        plt.show()
