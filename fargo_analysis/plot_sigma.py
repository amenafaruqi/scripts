import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.interpolate as interp
# matplotlib.use('TkAgg')
plt.style.use('default')
# plt.style.use(['../styles/publication.mplstyle'])
plt.style.use(['../styles/presentation.mplstyle', '../styles/darkbg2.mplstyle'])

timesteps = np.array([0, 5e-3, 2.6e-2])  # specify in Myr
print(timesteps)

# ================== Read in data at timesteps =======================
sim =  "planettest_highres"
simdir = f"/home/astro/phrkvg/simulations/planet_growth/{sim}/"
params_file = f'{simdir}/variables.par'
params_dict = {}

plotsizex = 3
plotsizey = int(len(timesteps)/plotsizex)+1

cm = plt.get_cmap('gist_rainbow')

param_lines = open(params_file).readlines()
for line in param_lines:
    if line.split():
        param_label, param_value = line.split()[0:2]
        params_dict.update([(param_label, param_value)])

nphi = int(params_dict['NX'])
nrad = int(params_dict['NY'])
rhodust = float(params_dict['RHO_DUST'])
f = float(params_dict['FLARINGINDEX'])
hr0 = float(params_dict['ASPECTRATIO'])      # aspect ratio at R=1AU
ndust = int(params_dict['NDUST'])
alpha = float(params_dict['ALPHA'])
mingsize = float(params_dict['MIN_GRAIN_SIZE'])
maxgsize = float(params_dict['MAX_GRAIN_SIZE'])
nss_coag = int(params_dict['NUMSUBSTEPS_COAG'])
densfloor = float(params_dict['DENSITY_FLOOR'])
dt_orbits = int(float(params_dict['DT'])/(2*np.pi))   # 2pi = 1 orbit = 1 yr
ninterm = float(params_dict['NINTERM'])               # number of dts between outputs
dt_outputs = dt_orbits*ninterm                        # time between outputs

outputs = timesteps*1e6/dt_outputs
print(outputs)

# FARGO initialises grains with sizes uniformly distributed across ndust bins in logspace
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
    for n in np.arange(ndust):
        dust_file = f"dustdens{n}_{int(t)}.dat" 
        sigma_dust[i,n] = np.fromfile(simdir+dust_file).reshape(nrad,nphi)/(1.125e-7)
    # print(f"t={t}: ", np.min(sigma_dust_azimsum[i]), np.max(sigma_dust_azimsum[i]))


sigma_dust_azimsum = np.sum(sigma_dust, axis=3)                  # sum over all phi 
sigma_gas_azimsum = np.sum(sigma_gas, axis=2)                    # sum over all phi   
sigma_dust_tot = np.sum(sigma_dust_azimsum, axis=1)              # sum over all dust sizes

# surfdens = sigma_gas[i]
# dens_first_wedge = surfdens[:,0].reshape(nrad,1)
# dens_additional = np.concatenate((surfdens[:,:],dens_first_wedge),axis=1)
sigma_gas_avg = sigma_gas_azimsum/nphi
sigma_dust_avg = sigma_dust_azimsum/nphi
sigma_dust_tot_avg = sigma_dust_tot/nphi


# ======================== Data arrays initialised ==============================

def sigma_at_r(r=10, t_index=-1):
    r_in_range = radii[(r*0.8<radii) & (radii<r*1.2)]                # values of r within % of chosen radius
    mini = np.where(radii == np.min(r_in_range))[0][0]
    maxi = np.where(radii == np.max(r_in_range))[0][0]
    sigma_in_r_range = sigma_dust_avg[t_index, :, mini:maxi]     # for chosen timestep, all grain sizes
    sigma_at_r0 = np.median(sigma_in_r_range, axis=1)                # avg sigma at chosen radius for all grain sizes
    return sigma_at_r0, mini, maxi

# ============ Fig 13, Brauer et al. 2008 / Fig 5, Birnstiel et al. 2012 ================

# fig, ax = plt.subplots(1, dpi=150)

# sigma_a_10au  = sigma_at_r(10)
# sigma_a_100au  = sigma_at_r(100)
# ax.plot(a, sigma_a_10au, label="R=10 AU")
# ax.plot(a, sigma_a_100au, label="R=100 AU")
# # ax.set_ylim(1e-12, 1e-2)
# # ax.set_xlim(2e-5, 1e2)
# ax.set_xlim(np.min(a), np.max(a))
# ax.set_xlabel("a (cm)")
# ax.set_ylabel("$\Sigma$ (g/cm$^{2}$)")
# ax.set_xscale("log")
# ax.set_yscale("log")
# ax.set_title('t = 1 Myr')
# ax.legend()
# # ax.set_title("NumSubsteps_Coag="+str(nss_coag))
# fig.savefig(f"./images/{sim}_sigmaatr.png")

# ======================= Fig 1, Birnstiel et al. 2012 ==============================

fig0 = plt.figure(figsize=(16,12))

R, A = np.meshgrid(radii, a)
# levels = np.linspace(-11,1,7)                   # Brauer 2008 levels
# levels = np.linspace(-7, 2, 10)                  # Birnstiel 2012 levels 
levels = np.linspace(-18, 6, 13)                   
a_St1 = (2/np.pi)*(sigma_gas_avg/rhodust)             # plot St=1 line
print(a_St1.shape)
fd = 0.55
ff = 0.37
uf = 10
hr = hr0*(radii**f)
cs = hr*(((2e30)*(6.67e-11))/(radii*1.5e11))**0.5

# size of largest grains in a fragmentation-dominated distribution
# a_frag = 100*ff*(2/(3*np.pi))*((uf**2)/(rhodust*1000*alpha))*(hr**-2)*(sigma_gas*10/((2e30)*(6.67e-11)))*(radii*1.5e11)  # from Birnstiel+2012

b = (uf**2/alpha)*(hr**-2)*(radii*1.5e11)/((2e30)*(6.67e-11)) # dimensionless
a_frag = (sigma_gas_avg/rhodust)*(3-(9-4*(b**2))**0.5)/(np.pi*b)

# size of largest grains in a drift-dominated distribution
# a_drift = 100*(2/(rhodust*1000*np.pi))*(hr**-2)*sigma_dust_tot*10*(2/3)   # from Birnstiel+2012

p = (sigma_gas_avg*(cs**2)/((2*np.pi)**0.5))*(hr**-1)*((radii*1.5e11)**-1)
pad = np.empty((len(timesteps), 1))
gamma = (radii/p)*np.abs(np.append(np.diff(p)/np.diff(radii), pad, axis=1))
a_drift = (2/np.pi)*(sigma_dust_tot/rhodust)*(1/gamma)*(hr**-2)

for i, o in enumerate(outputs):
    ax0 = fig0.add_subplot(plotsizey, plotsizex, i+1)
    sigmas = sigma_dust_avg[i]
    c = ax0.contourf(R, A, np.log10(sigmas), cmap="Greys", levels=levels)
    ax0.set_ylim(1e-4, 1e2)
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.set_title(f"{int(o*dt_outputs)} orbits")
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

# fig0.suptitle(sim)
fig0.tight_layout()

fig0.subplots_adjust(right=0.89, hspace=0.3)
cbar_ax = fig0.add_axes([0.91, 0.53, 0.02, 0.4])
fig0.colorbar(c, cax=cbar_ax, orientation="vertical", label="log$[\Sigma (g/cm^{{2}})]$")
ax0.legend(loc="upper right")

fig0.savefig(f"./images/{sim}_contour.png")

# ======================= Fig 8, Birnstiel et al. 2012 =========================
fig1, ax1 = plt.subplots(figsize=(7,5))
ax1.set_prop_cycle(color=[cm(1.*i/8) for i in range(0,len(timesteps)+1)])

for i,o in enumerate(timesteps):
    dustgasratio = sigma_dust_tot_avg[i]/sigma_gas_avg[i]
    ax1.plot(radii, dustgasratio, label=f"t={o} Myr")
# ax1.axvline(20, linestyle='dashed', color='black')

ax1.legend()
ax1.set_xlim(np.min(radii), np.max(radii))
ax1.set_xlabel("R (AU)")
ax1.set_ylabel("dust-gas ratio")
ax1.set_xscale("log")
ax1.set_yscale("log")

fig1.tight_layout()
fig1.savefig(f"./images/{sim}_dustgasratio.png")


# ================== Gas Sigma Profile ===================

fig2, ax2 = plt.subplots(figsize=(7,6))
ax2.set_prop_cycle(color=[cm(1.*i/8) for i in range(0,len(timesteps)+1)])

for i, o in enumerate(outputs):
    ax2.plot(radii, sigma_gas_avg[i], label=f"{int(o*dt_outputs)} orbits")

    # surfdens = sigma_gas[i]
    # dens_first_wedge = surfdens[:,0].reshape(nrad,1)
    # dens_additional = np.concatenate((surfdens[:,:],dens_first_wedge),axis=1)
    # sigmag_tot = np.sum(dens_additional, axis=1)
    # sigmag_azi = sigmag_tot/nphi
    # ax2.plot(radii, np.log10(sigmag_azi), label=f"{t} Myr")

ax2.axvline(20, linestyle='dashed', color='black')

ax2.set_xlabel("R (AU)")
ax2.set_ylabel("$\Sigma_{{gas}} (g/cm^{{2}})$")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.legend()

ax2.set_xlim(np.min(radii), np.max(radii))
# ax2.set_ylim(np.min(sigma_gas_azimsum), np.max(sigma_gas_azimsum))

fig2.tight_layout()
fig2.savefig(f"./images/{sim}_sigmagas.png")

# ================== Dust sigma evolution over time (for different a) ===================

fig3, ax3 = plt.subplots(figsize=(9,6))
ax3.set_prop_cycle(color=[cm(1.*i/8) for i in range(0,len(timesteps)+1)])

for i, t in enumerate(timesteps):
    c = next(ax3._get_lines.prop_cycler)['color']
    sigma_at_100, mini100, maxi100 = sigma_at_r(20, i)
    ax3.plot(a, sigma_at_100, linestyle='solid', label=f"t={t} Myr", color=c)
    # a_frag100 = np.median(a_frag[i, mini100:maxi100])
    # a_drift100 = np.median(a_drift[i, mini100:maxi100])
    # a1 = np.min(a_frag100, a_drift100)
    # ax3.axvline(a1, linestyle='dashed', color=c)

ax3.set_xlabel("a (cm) ")
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlim(np.min(a), np.max(a))
ax3.set_ylim(densfloor/(1.125e-7),)
ax3.set_ylabel("$\Sigma_{{dust}} (g/cm^{{2}})$ at R=20AU")
ax3.legend()

# ax3.set_title(sim)
fig3.tight_layout()
fig3.savefig(f"./images/{sim}_sigmaatr_timeevo.png")

# ================= Total  dust sigma evolution over time  ===============

fig4, ax4 = plt.subplots(1, dpi=150)
ax4.set_prop_cycle(color=[cm(1.*i/8) for i in range(0,len(timesteps)+1)])

delta_r = radii[1]-radii[0]   # works for a linear grid only!!
M_disc = np.zeros(len(timesteps))

for i in range(len(timesteps)):
    M_cell = sigma_dust_tot[i]*(radii*1.5e13)*2*np.pi*(delta_r*1.5e13)
    M_disc[i] = np.sum(M_cell)

ax4.plot(timesteps, M_disc)
# print(np.min(M_disc), np.max(M_disc), M_disc[-10:])
ax4.set_xlabel("t (Myr)")
ax4.set_ylabel("Total  $M_{{dust}} (g/cm^{{2}})$")
ax4.set_yscale("log")
ax4.set_xscale("log")
# ax4.set_xlim(0,1)
# ax4.set_ylim(3, 40)

ax4.set_title(sim)
fig4.tight_layout()


# ======================= Compare Sigma profiles =========================

# fig2, ax2 = plt.subplots(2)

# r = np.linspace(1,150,500)

# p_birn = 1
# sigma0_birn = 10.4                                  # value read off from fig.3 Birnstiel 2012, Sigma at 60AU
# rc_birn = 60

# p_br = 0.8
# sigma0_br = 20

# p = float(params_dict['SIGMASLOPE'])
# sigma0 = float(params_dict['SIGMA0'])/(1.125e-7)             # convert to g/cm2
# rc = float(params_dict['SIGMACUTOFF'])

# sigma = sigma0*((r)**(-p))*(np.exp(-(r/rc)**(2-p)))                                        # Sigma from simulation
# sigma_birn = sigma0_birn*((r/rc_birn)**(-p_birn))*(np.exp(-(r/rc_birn)**(2-p_birn)))       # Sigma from Birnstiel 2012 (approx.)
# sigma_br = sigma0_br*(r**-p_br)                                                            # Sigma from Brauer 2008

# ax2[0].plot(r, sigma, label="Simulation input")
# ax2[0].plot(r, sigma_birn, label="Birnstiel")
# ax2[0].plot(r, sigma_br, label="Brauer")
# ax2[0].plot(radii, sigma_gas[-1]+sigma_dust_tot[-1], label=f"Simulation $\Sigma$ at t={timesteps[-1]} Myr")

# ax2[0].set_xscale("log")
# ax2[0].set_yscale("log")
# ax2[0].legend()
# ax2[0].set_xlabel('R [AU]')
# ax2[0].set_ylabel('$\Sigma$ [$gcm^{-2}$]')


# # Trialling Sigma profiles based on fig.3 in Birnstiel 2012
# rc_trial = 60               # stated in Birnstiel 2012
# p_trial = 1                 # stated in Birnstiel 2012
# birnstieldata = np.loadtxt("BirnstielData.csv", delimiter=",")
# sigma0_trial = interp.interp1d(birnstieldata[:,0], birnstieldata[:,1])(rc_trial)*(1/rc_trial)**(-p_trial)
# sigma_trial = sigma0_trial*((r)**(-p_trial))*(np.exp(-(r/rc_trial)**(2-p_trial)))

# ax2[1].scatter(birnstieldata[:,0], birnstieldata[:,1], c="r",  label="Birnstiel")
# ax2[1].plot(r, sigma_trial, label=f"p={p_trial}, $\Sigma_{{0}}$={round(sigma0_trial)} $gcm^{{-2}}$, $R_{{c}}$={rc_trial} AU")
# ax2[1].plot(radii, sigma_gas[-1]+sigma_dust_tot[-1], label=f"Simulation $\Sigma$ at t={timesteps[-1]} Myr")
# ax2[1].set_xscale("log")
# ax2[1].set_yscale("log")
# ax2[1].legend()
# ax2[1].set_xlabel('R [AU]')
# ax2[1].set_ylabel('$\Sigma$ [$gcm^{-2}$]')

# =========================================================================

# plt.show()


