import numpy as np
import matplotlib.pyplot as plt
import argparse
import matplotlib as mpl
plt.style.use('default')

import warnings
warnings.filterwarnings("ignore")

# ======================= Compute migration timescales ============================

def calc_t_migration(mass, radius):
        hr = hr0*(radius**f)
        R_h = radius*((mass*1e-6)**(1/3))
        p_crida_Rp = 0.75*(hr*radius)/R_h + 50*alpha*(hr**2)/(mass*3e-6)   # Crida parameter for R=Rp

        p = -1*sigmaslope
        q = -2*f
        gamma = 5/3
        zeta = q - (gamma-1)*p
        Omega = np.sqrt(4*(np.pi**2)/(radius**3))
        Sigmap = sigma0*(float(radius)**-sigmaslope)*np.exp(-radius/Rc)     # assuming sigma profile is not changing!
        if p_crida_Rp > 1:     # type I migration regime if planet is not gap-opening
            Gamma = (-2.5-(1.7*q)+(p/10)+(1.1*(1.5-p))+(7.9*zeta/gamma))*((mass*3e-6/hr)**2)*Sigmap*(radius**4)*(Omega**2)/gamma
            tau = ((radius**2)*Omega*mass*3e-6/(2*Gamma))*1e-6    # convert to Myr
        else:                  # type II migration regime
            h = hr*radius
            nu = alpha*Omega*(h**2)
            B = mass*3e-6/(np.pi*Sigmap*(radius**2))
            tau = ((radius**2)*(1+B)/nu)*1e-6
        return tau


# ======================= Plot Regimes for Gap-Opening ============================

def plot_Mp_regimes():
    print("Plotting mass regimes....")

    mps = np.array([12, 25, 60, 120, 160]) 
    fig0,ax0 = plt.subplots(figsize=(12,8))
    m_range = np.linspace(0,200,50)
    r_range = np.linspace(0,200,201)

    # Plot pebble isolation criterion (eq. 26 from Bitsch+2018)
    hr1d = hr0*(r_range**f)
    dlogPdlogR = f - sigmaslope - 2     # taken from eq. 9 from Bitsch et al. 2018, accounting for their s being -ve. 
    f_fit = ((hr1d/0.05)**3)*(0.34*(np.log10(0.001)/np.log10(alpha))**4 + 0.66)*(1-((dlogPdlogR+2.5)/6))
    alpha_St = 0.01
    M_iso = 25*f_fit
    Pi_crit = alpha_St/2
    Lambda = 0.00476/f_fit
    M_iso += Pi_crit/Lambda                  # PIM considering diffusion
    ax0.fill_between(r_range, M_iso, 200, color='yellow', alpha=0.55)
    ax0.fill_between(r_range, M_iso, 0, color='orangered', alpha=0.73)
    
    # Plot PIM for different alpha/St values
    alpha_St = np.array([0.001, 0.1])
    clrs = ["dotted", "dashdot"]
    for ai,a in enumerate(alpha_St):
        M_iso = 25*f_fit
        Pi_crit = a/2
        Lambda = 0.00476/f_fit
        M_iso += Pi_crit/Lambda                  # PIM considering diffusion
        ax0.plot(r_range, M_iso, linestyle=clrs[ai], label=f"St = {0.001/a}", color="dodgerblue")

    # Plot gap-opening criterion region (Crida et al. 2006)
    Ms, Rs = np.meshgrid(m_range, r_range)
    hr = hr0*(Rs**f)
    R_h = Rs*((Ms*1e-6)**(1/3))
    p_crida = 0.75*(hr*Rs)/R_h + 50*alpha*(hr**2)/(Ms*3e-6)
    ax0.contourf(Rs, Ms, p_crida, levels=[0,1], cmap="summer")

    # Plot planet migration tracks
    Rp0 = 40     # remove this afterwards, for testing only
    dt = 0.1     # timestep to recalculate planet's location

    # Plot radius where dust drift velocity = planet migration velocity
    St = 0.01                            # pebble Stokes number
    q = -2*f
    gamma = 5/3
    zeta = q - (gamma-1)*-sigmaslope
    Omega = np.sqrt(4*(np.pi**2)/(Rp0**3))
    Sigmap = sigma0*(float(Rp0)**-sigmaslope)*np.exp(-Rp0/Rc)
    A = (-2.5-(1.7*q)+(-sigmaslope/10)+(1.1*(1.5+sigmaslope))+(7.9*zeta/gamma))
    Sigmap = sigma0*(r_range**-sigmaslope)*np.exp(-r_range/Rc)

    M_drift_mig = (gamma*St/(2*A))*((hr0*(r_range**f))**4)*dlogPdlogR/(Sigmap*(r_range**2))
    M_drift_mig = np.abs(M_drift_mig)/3e-6    # convert to Earth masses
    ax0.plot(r_range,  M_drift_mig, c='mediumblue', label=f"$v_{{drift}}=v_{{mig}}$ for St={St}")
    ax0.legend(facecolor='grey', framealpha=0.2)

    # Calculate location of inner damping zone
    Rid = (2*(rmin**-1.5)/3)**(-2/3)

    # Calculate migration timescale for each planet mass
    for mp in mps:
        tau = calc_t_migration(mp, Rp0)     # migration timescale at starting position of planet
        Rp_t = Rp0
        R_21 = (2**(-2/3))*Rp_t             # location  of  2:1 resonance
        t = 0
        while R_21 > Rid and t < 0.5:
            t += dt
            Rp_new = Rp_t*np.exp(-dt/tau)
            if Rp_new > rmin:
                Rp_t = Rp_new
            R_21 = (2**(-2/3))*Rp_t              # location  of  2:1 resonance
            tau = calc_t_migration(mp, Rp_t)
        ax0.scatter(Rp_t, mp, marker='<', color='k')   # mark final location of planet

        # Plot migration tracks of planet
        ax0.hlines(y=mp, xmin=Rp_t, xmax=Rp0, linewidth=2, color='k', linestyle='--')
        # Plot planet's initial location
        ax0.scatter(Rp0, mp, marker='o', color='k')
    
    # Plot boundary of inner damping zone
    ax0.vlines(x=Rid, ymin=0, ymax = 200, linewidth=2, color='crimson', linestyle=':')

    ax0.set_xlim(1,70)
    ax0.set_ylim(0,170)
    ax0.set_xlabel("$R_{{p}}$ (AU)")
    ax0.set_ylabel("$M_{{p}} (M_\oplus)$")
    fig0.tight_layout()
    fig0.savefig(f"{plots_savedir}/Mp_regimes.png", dpi=200)


# ========== Plot individual PIM terms for different grain sizes ============

def plot_Miso_vs_R():
    print("Plotting isolation mass by grain size....")

    fig0,ax0 = plt.subplots(figsize=(12,8))
    r_range = np.linspace(0,200,201)

    # term 1
    hr1d = hr0*(r_range**f)
    dlogPdlogR = f - sigmaslope - 2     # taken from eq. 9 from Bitsch et al. 2018, accounting for their s being -ve. 
    f_fit = ((hr1d/0.05)**3)*(0.34*(np.log10(0.001)/np.log10(alpha))**4 + 0.66)*(1-((dlogPdlogR+2.5)/6))
    M_iso = 25*f_fit
    ax0.plot(r_range, M_iso, label="$M_{iso}^{\dag}$")

    # term 2, introduced in eq. 26
    grain_sizes = np.array([0.01,0.1,1])
    sigma_g = sigma0*(r_range**-sigmaslope)/(1.125e-7)
    Lambda = 0.00476/f_fit
    for g in grain_sizes:
        St = g*rhodust*np.pi/(2*sigma_g)

        Pi_crit = alpha/(2*St)
        M_iso_newterm = Pi_crit/Lambda                  # PIM considering diffusion
        ax0.plot(r_range, M_iso_newterm, label=f"$\Pi/\lambda$, a={g}cm")
        c = plt.gca().lines[-1].get_color()
        ax0.plot(r_range, M_iso+M_iso_newterm, label=f"$M_{{iso}}$, a={g}cm", linestyle="dashed", color=c, alpha=0.7)

    ax0.legend()
    ax0.set_xlim(1,180)
    ax0.set_ylim(0,100)
    ax0.set_xlabel("$R_{{p}}$ (AU)")
    ax0.set_ylabel("$M_{iso} (M_\oplus)$")
    fig0.savefig(f"{plots_savedir}/analytics/Miso.png", dpi=200)


# ================ Plot mass regimes from Dipierro + 2017 ====================

def plot_Dipierro_Mp_regimes():  # figure 2 from dipierro 2017 i.e. as function of Mp and a
    print("Plotting Dipierro mass regimes....")

    hr = hr0*(Rp**f)                                   # aspect ratio at R=Rp
    factor = 0.4                                       # can be 0.1 (Rafikov+ 2012), 0.15 (Lin+ 1979), 0.4 (Goldreich+ 1979)
    Mp_th = 3*((hr)**3)                                # Mp in solar masses
    Mp_visc = ((3/(2*factor))**0.5)*(hr**2.5)*alpha**0.5
    Mp_gap = np.max((0.2*Mp_th, Mp_visc))       # convert Mp to Mearth
    Mp_lim = 0.1*Mp_th
    
    St = np.arange(1e-1,1e2)                # chosen arbitrarily based on the condition that "large" grains is St>alpha and St>dgr.
    sig = -(sigmaslope + f + 1.5)
    zeta = ((2*factor)**-1.5)/9
    z = (-sig+(6+3*sig)*alpha/St)/(1+dgr+(dgr/(St**2)))
    Mp_dustonly = (3**1.5)*zeta*((z/St)**1.5)*(hr**3)

    fig0,ax0 = plt.subplots(figsize=(9,6))

    ax0.set_xlim(1e-1,1e2)
    ax0.set_ylim(1e-6,1e-2)
    ax0.plot(St,Mp_dustonly, color="red")
    ax0.axhline(Mp_gap, color="green")
    ax0.axhline(Mp_lim, color="skyblue")
    ax0.fill_between(St, Mp_dustonly, 1e-2, color="red")
    ax0.fill_between(St, Mp_lim, 1e-2, color="skyblue")
    ax0.fill_between(St, Mp_gap, 1e-2, color="green")
    ax0.text(2e-1,Mp_gap*1.1,f"$M_{{p,gap}}$ = {round(Mp_gap/3e-6,1)} $M_\oplus$")
    ax0.text(2e-1,Mp_lim*1.05,f"$M_{{p,lim}}$ = {round(Mp_lim/3e-6,1)} $M_\oplus$")
    ax0.text(50,Mp_lim*0.7, f"$M_{{dustonly}}$")   
    ax0.set_xlabel("St")
    ax0.set_ylabel("$M_{p}/M_\odot$")
    ax0.set_xscale("log")
    ax0.set_yscale("log")

    def Msun_to_Mearth(M):
        return M/3e-6
    
    def Mearth_to_Msun(M):
        return M*3e-6
    
    secax_y = ax0.secondary_yaxis('right', functions=(Msun_to_Mearth, Mearth_to_Msun))
    secax_y.set_ylabel("$M_{p}/M_\oplus$")
    secax_y.set_yscale("log")

    fig0.savefig(f"{plots_savedir}/Dipierro_regimes.png")


# ======================= Plot Drift Timescale ============================

def plot_tdrift(): #as function of R and a
    print("Plotting drift timescale....")

    A,R = np.meshgrid(a,radii)
    hr = hr0*(R**f)                                   # aspect ratio
    cs = hr*(((2e30)*(6.67e-11))/(R*1.5e11))**0.5     # [m/s]
    # p = (sigma_gas_1D*(cs**2)/((2*np.pi)**0.5))*(hr**-1)*((radii*1.5e11)**-1)
    # gamma = (radii/p)*np.abs(np.append(np.diff(p)/np.diff(radii), 0))
    gamma = np.abs(f - sigmaslope - 2)     # taken from eq. 9 from Bitsch et al. 2018, accounting for their s being -ve. 
    # tau_drift = ((2/(np.pi*rhodust))*(np.dot((radii*sigma_gas_azimsum/(gamma*hr*cs)).reshape(nrad,1), a.reshape(1,ndust))*1.5e11))*3.17e-14    # drift timescale in Myr
    fig,ax = plt.subplots(figsize=(11,6))
    levels = np.linspace(-6, 6, 13)                   
    sigma_gas_1D_rep = np.tile(sigma_gas_1D, ndust).reshape(nrad,ndust)  # repeat sigma gas for each dust size for ease of calculation

    tau_drift = (3.1688e-14)*(2/np.pi)*(sigma_gas_1D_rep/(rhodust*A))*(R*1.5e11/cs)/gamma
    print(tau_drift)
    con = ax.contourf(R, A, np.log10(tau_drift), cmap="YlGnBu")#, levels=levels)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("a (cm)")
    ax.set_xlabel("R (AU)")

    fig.subplots_adjust(right=0.89, hspace=0.3)
    cbar_ax = fig.add_axes([0.91, 0.53, 0.02, 0.4])
    fig.colorbar(con, cax=cbar_ax, orientation="vertical", label="log$(\\tau_{drift}/1Myr)$")
    for rp in rps:
        ax.axvline(rp, linestyle='dashed', color="k")    
    
    fig.savefig(f"{plots_savedir}/{sim}_tdrift.png")


# ======================= Plot Migration Timescale ============================

def plot_tmigration():  # as function of R
    p = -1*sigmaslope
    q = -2*f
    gamma = 5/3
    zeta = q - (gamma-1)*p
    Rs = np.arange(10,101,5)                         # in units of AU
    Ms = np.arange(1,301,10)*3e-6                   # in units of Msun
    tauI_arr_R = np.zeros(len(Rs))
    tauI_arr_M = np.zeros(len(Ms))
    tauII_arr_R = np.zeros(len(Rs))
    tauII_arr_M = np.zeros(len(Ms))
    for i,R in enumerate(Rs):
        Omega = np.sqrt(4*(np.pi**2)/(R**3))
        hr = hr0*(R**f)                                   # aspect ratio H/R at Rp
        Sigmap = sigma0*(float(R)**-sigmaslope)*np.exp(-R/Rc)
        Gamma = (-2.5-(1.7*q)+(p/10)+(1.1*(1.5-p))+(7.9*zeta/gamma))*((Mp/hr)**2)*Sigmap*(R**4)*(Omega**2)/gamma
        tauI_arr_R[i] = ((R**2)*Omega*Mp/(2*Gamma))*1e-6    # convert to Myr

        h = hr*R
        nu = alpha*Omega*(h**2)
        B = Mp/(np.pi*Sigmap*(R**2))
        B=0
        tauII_arr_R[i] = ((R**2)*(1+B)/nu)*1e-6

    for i,M in enumerate(Ms):
        Omega = np.sqrt(4*(np.pi**2)/(Rp**3))
        hr = hr0*(Rp**f)                                   # aspect ratio H/R at Rp
        Sigmap = sigma0*(float(Rp)**-sigmaslope)*np.exp(-Rp/Rc)
        Gamma = (-2.5-(1.7*q)+(p/10)+(1.1*(1.5-p))+(7.9*zeta/gamma))*((M/hr)**2)*Sigmap*(Rp**4)*(Omega**2)/gamma
        tauI_arr_M[i] = ((Rp**2)*Omega*M/(2*Gamma))*1e-6

        h = hr*Rp
        nu = alpha*Omega*(h**2)
        B = M/(np.pi*Sigmap*(Rp**2))
        B=0
        print(B,nu, Omega)
        tauII_arr_M[i] = ((Rp**2)*(1+B)/nu)*1e-6
    
    print(tauI_arr_M[::4])
    print(tauII_arr_M[::4])
    print("Plotting migration timescales....")

    fig2, ax2 = plt.subplots(nrows=2, figsize=(7,6))
    ax2[0].plot(Rs,tauI_arr_R, color='blue', label="Type I")
    ax2[0].plot(Rs,tauII_arr_R, color='red', label="Type II")
    ax2[0].set_xlabel("R (AU)")
    ax2[0].set_ylabel("$\\tau$ (Myr)")
    # ax2[0].set_xlim(np.min(np.log10(Rs)), np.max(np.log10(Rs)))
    ax2[0].axvline(Rp, linestyle='dashed', color='black')
    ax2[0].set_yscale("log")
    ax2[0].set_xscale("log")
    ax2[0].legend()

    ax2[1].plot(Ms/(3e-6),tauI_arr_M, color='blue')
    ax2[1].plot(Ms/(3e-6),tauII_arr_M, color='red')
    ax2[1].set_xlabel("$M (M_\oplus)$")
    ax2[1].set_ylabel("$\\tau$ (Myr)")
    # ax2[1].set_xlim(np.min(np.log10(Ms/3e-6)), np.max(np.log10(Ms/3e-6)))
    ax2[1].axvline(Mp/3e-6, linestyle='dashed', color='black')
    ax2[1].set_yscale("log")
    ax2[1].set_xscale("log")

    fig2.tight_layout()
    fig2.savefig(f"{plots_savedir}/{sim}_tmig.png")


# ======================= Plot Growth Timescale ============================

def plot_tgrowth():  #as function of R and a
    R,A = np.meshgrid(radii,a)
    hr = hr0*(radii**0.25)
    h = hr*radii
    cs = hr*(((2e30)*(6.67e-11))/(radii*1.5e11))**0.5     # [m/s]
    St = A*rhodust*np.pi/(sigma_gas_1D*2)
    const = ((4/3)*(2*np.pi/3)**0.5)
    tau_growth = const*A*rhodust*h/(sigma_dust_1D*cs*((alpha*St)**0.5))
    print("Plotting growth timescale....")

    fig3,ax3 = plt.subplots(figsize=(10,6))
    levels = np.linspace(-4, 16, 11)                   
    con = ax3.contourf(R,A, np.log10(tau_growth), cmap="YlGnBu", levels=levels)
    
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_ylabel("a (cm)")
    ax3.set_xlabel("R (AU)")

    fig3.subplots_adjust(right=0.89, hspace=0.35)
    cbar_ax = fig3.add_axes([0.91, 0.53, 0.02, 0.4])
    fig3.colorbar(con, cax=cbar_ax, orientation="vertical", label="log$[\\tau_{growth} (Myr)]$")
    ax3.axvline(Rp, linestyle='dashed', color='black')

    fig3.savefig(f"{plots_savedir}/{sim}_tgrowth.png")



# ==========================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analytic calculations based on model inputs (not outputs)', prefix_chars='-')

    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/astro/phrkvg/simulations/lowres_models"],help="working directory containing simulations")
    parser.add_argument('-sim', metavar='sim', type=str, nargs=1, default=["Mp60_stat"] ,help="simulation directory containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images/analytics" ,help="directory to save plots to")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-plots', metavar='plots',default=["mpreg"], type=str, nargs="*" ,help="plots to produce")
    parser.add_argument('-porbits', action="store_true")
    parser.add_argument('-style', metavar='style', type=str, nargs="*", default=["publication"] ,help="style sheet to apply to plots")

    args = parser.parse_args()
    wd = args.wd[0]
    sim = args.sim[0]
    simdir = f"{wd}/{sim}/"
    plot_window = args.plot_window
    plots = args.plots
    plots_savedir = args.savedir
    p_orbits = args.porbits
    style = args.style

    if plot_window:
        mpl.use('TkAgg') 

    params_file = f'{simdir}/variables.par'
    params_dict = {}

    param_lines = open(params_file).readlines()
    for line in param_lines:
        if line.split():
            param_label, param_value = line.split()[0:2]
            params_dict.update([(param_label, param_value)])

    nphi = int(params_dict['NX'])
    nrad = int(params_dict['NY'])
    rmin = int(params_dict['YMIN'])
    rmax = int(params_dict['YMAX'])
    f = float(params_dict['FLARINGINDEX'])
    hr0 = float(params_dict['ASPECTRATIO'])      # aspect ratio at R=1AU
    alpha = float(params_dict['ALPHA'])
    spacing = str(params_dict['SPACING'])
    ndust = int(params_dict['NDUST'])
    sigma0 = float(params_dict['SIGMA0'])
    Rc = int(params_dict['SIGMACUTOFF'])
    sigmaslope = int(params_dict['SIGMASLOPE'])
    mingsize = float(params_dict['MIN_GRAIN_SIZE'])
    maxgsize = float(params_dict['MAX_GRAIN_SIZE'])
    rhodust = float(params_dict['RHO_DUST'])
    dgr = float(params_dict['DUST_TO_GAS'])

    with open(f"{simdir}/planet.cfg", 'r') as pfile:
        num_planets = len(pfile.readlines()) - 5        # ignore 5 lines of headers
    
    rps = np.zeros((num_planets))
    for n in range(num_planets):
        planet_data = np.unique(np.loadtxt(f"{simdir}/planet{n}.dat"), axis=0)
        xp, yp = planet_data[0][1], planet_data[0][2]
        # print(xp,yp)
        rps[n] = ((xp**2) + (yp**2))**0.5

    a = np.logspace(np.log10(mingsize),    
                    np.log10(maxgsize),
                    ndust+1)

    a = (0.5*(a[1:] + a[:-1]))                                   # grain sizes in middles of bins (in cm)

    r_cells = np.loadtxt(f'{simdir}/domain_y.dat')[3:-3]             # ignore ghost cells
    phi_cells = np.loadtxt(f'{simdir}/domain_x.dat')[3:-3]
    if spacing == "Linear":
        radii = np.array([(r_cells[n]+r_cells[n+1])/2 for n in range(len(r_cells)-1)])
    else:     # Log grid
        radii = np.array([np.exp((np.log(r_cells[n])+np.log(r_cells[n+1]))/2) for n in range(len(r_cells)-1)])
    phis = np.array([(phi_cells[n]+phi_cells[n+1])/2 for n in range(len(phi_cells)-1)])

    gasfile = f"gasdens1.dat" 
    sigma_gas = np.fromfile(simdir+gasfile).reshape(nrad,nphi)/(1.125e-7)         # convert back to g/cm2
    sigma_gas_azimsum = np.sum(sigma_gas, axis=1)                                 # sum over all phi   
    sigma_gas_1D = sigma_gas_azimsum/nphi                                         # dimensions: (nrad) 

    sigma_dust = np.zeros((ndust, nrad, nphi))
    n_grains = np.zeros((ndust))

    for n in np.arange(ndust):
        dust_file = f"dustdens{n}_1.dat" 
        sigma_dust[n] = np.fromfile(simdir+dust_file).reshape(nrad,nphi)/(1.125e-7)
        n_grains[n] = np.sum(sigma_dust[n])*3/(a[n]*rhodust*4)                  # number of dust grains at size a and time t

    sigma_dust_azimsum = np.sum(sigma_dust, axis=2)                  # sum over all phi 
    sigma_dust_1D = sigma_dust_azimsum/nphi                          # dimensions: (ndust, nrad)
    sigma_gas_azimsum = np.sum(sigma_gas, axis=1)                    # sum over all phi   

    if not plot_window:
        for s in style:
            plt.style.use([f"../styles/{s}.mplstyle"])

    if "drift" in plots:
        plot_tdrift()
    if "growth" in plots:
        plot_tgrowth()
    if "mig" in plots:
        plot_tmigration()
    if "mpreg" in plots:
        plot_Mp_regimes()
    if "miso" in plots:
        plot_Miso_vs_R()

    if plot_window:
        plt.show()
