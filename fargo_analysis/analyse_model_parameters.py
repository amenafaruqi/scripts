import numpy as np
import matplotlib.pyplot as plt
import argparse
import matplotlib as mpl
plt.style.use('default')

import warnings
warnings.filterwarnings("ignore")

# ======================= Plot Regimes for Gap-Opening ============================

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
        # print(mass)
        if p_crida_Rp > 1:     # type I migration regime
            # print("Type I")
            Gamma = (-2.5-(1.7*q)+(p/10)+(1.1*(1.5-p))+(7.9*zeta/gamma))*((mass*3e-6/hr)**2)*Sigmap*(radius**4)*(Omega**2)/gamma
            tau = ((radius**2)*Omega*mass*3e-6/(2*Gamma))*1e-6    # convert to Myr
        else:                  # type II migration regime
            h = hr*radius
            nu = alpha*Omega*(h**2)
            B = mass*3e-6/(np.pi*Sigmap*(radius**2))
            tau = ((radius**2)*(1+B)/nu)*1e-6
            # print("Type II")
        return tau


def plot_Mp_regimes():
    mps = np.array([12, 25, 35, 120, 160]) 
    # mps = np.array([12,60,120])
    fig0,ax0 = plt.subplots(figsize=(9,6))
    m_range = np.linspace(0,200,50)
    r_range = np.linspace(0,200,201)

    # Plot pebble isolation criterion (Lambrechts et al. 2014)
    hr1d = hr0*(r_range**f)
    dlogPdlogR = f - sigmaslope - 2     # taken from eq. 9 from Bitsch et al. 2018, accounting for their s being -ve. 
    f_fit = ((hr1d/0.05)**3)*(0.34*(np.log10(0.001)/np.log10(alpha))**4 + 0.66)*(1-((dlogPdlogR+2.5)/6))
    M_iso = 25*f_fit
    alpha_St = 0.01
    Pi_crit = alpha_St/2
    Lambda = 0.00476/f_fit
    M_iso += Pi_crit/Lambda                  # PIM considering diffusion
    ax0.fill_between(r_range, M_iso, 200, color='yellow', alpha=0.4)
    ax0.fill_between(r_range, M_iso, 0, color='coral', alpha=0.8)

    # Plot gap-opening criterion region (Crida et al. 2006)
    Ms, Rs = np.meshgrid(m_range, r_range)
    hr = hr0*(Rs**f)
    R_h = Rs*((Ms*1e-6)**(1/3))
    p_crida = 0.75*(hr*Rs)/R_h + 50*alpha*(hr**2)/(Ms*3e-6)
    ax0.contourf(Rs, Ms, p_crida, levels=[0,1])

    # Plot planet migration tracks
    Rp0 = 40     # remove this afterwards, for testing only
    dt = 0.1     # timestep to recalculate planet's location

    # Plot radius where dust drift velocity = planet migration velocity

    St = 0.01                            # pebble Stokes number
    # u_drift = dlogPdlogR*(hr0**2)*2*np.pi/(St+(1/St))   # in AU/yr
    # print(u_drift)
    a_pebble = 1*6.68e-14                # 1 cm pebble size in AU
    q = -2*f
    gamma = 5/3
    zeta = q - (gamma-1)*-sigmaslope
    Omega = np.sqrt(4*(np.pi**2)/(Rp0**3))
    # Sigmap = sigma0*(float(Rp0)**-sigmaslope)*np.exp(-Rp0/Rc)
    A = (-2.5-(1.7*q)+(-sigmaslope/10)+(1.1*(1.5+sigmaslope))+(7.9*zeta/gamma))
    # M_drift_mig = u_drift*((hr0*(Rp0**f))**2)*gamma/(2*A*Sigmap*Omega*Rp0**3)
    # print(M_drift_mig/3e-6)  # convert to Earth masses
    Sigmap = sigma0*(r_range**-sigmaslope)*np.exp(-r_range/Rc)

    M_drift_mig = (gamma*St/(2*A))*((hr0*(r_range**f))**4)*dlogPdlogR/(Sigmap*(r_range**2))
    # M_drift_mig = (np.pi*(rhodust*1.683e6)*gamma*a_pebble/(4*A))*((hr0*(r_range**f))**4)*dlogPdlogR/((Sigmap*r_range)**2)
    M_drift_mig = np.abs(M_drift_mig)/3e-6    # convert to Earth masses
    ax0.plot(r_range,  M_drift_mig, 'b')

    # Calculate location of inner damping zone
    Rid = (2*(rmin**-1.5)/3)**(-2/3)

    # Calculate migration timescale for each planet mass
    for mp in mps:
        tau = calc_t_migration(mp, Rp0)     # tau at starting position of planet
        Rp_t = Rp0
        R_21 = (2**(-2/3))*Rp_t
        t = 0
        while R_21 > Rid and t < 0.5:
            t += dt
            Rp_new = Rp_t*np.exp(-dt/tau)
            if Rp_new > rmin:
                Rp_t = Rp_new
            R_21 = (2**(-2/3))*Rp_t              # location  of  2:1 resonance
            tau = calc_t_migration(mp, Rp_t)
            # print(mp, " ---- \n")
            # print(Rp_t, R_21)
            ax0.scatter(Rp_t, mp, marker='|', color='k')

        ax0.hlines(y=mp, xmin=Rp_t, xmax=Rp0, linewidth=2, color='k', linestyle='--')
        # ax0.vlines(x=R_21, ymin=0, ymax = 200, linewidth=2, color='r', linestyle=':')
        ax0.scatter(Rp0, mp, marker='x', color='k')
    
    ax0.set_xlim(1,100)
    ax0.set_ylim(0,200)
    ax0.set_xlabel("R (AU)")
    ax0.set_ylabel("$M_{{p}} (M_\oplus)$")
    fig0.savefig(f"{plots_savedir}/Mp_regimes.png", dpi=200)


def plot_Dipierro_Mp_regimes():  #figure 2 from dipierro 2017 i.e. as function of Mp and a
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
    # ax0.plot(St,Mp_dustonly, color="red")
    # ax0.axhline(Mp_gap, color="green")
    # ax0.axhline(Mp_lim, color="skyblue")
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

    fig0.savefig(f"{plots_savedir}/{sim}_Mpregimes.png")


def plot_eta():
    r_range = np.linspace(0,150,151)
    dlogPdlogR = f - sigmaslope - 2
    hr = hr0*(r_range**0.25)
    eta =  -0.5*(hr**2)*dlogPdlogR



# ======================= Plot Drift Timescale ============================

def plot_tdrift(): #as function of R and a
    hr = hr0*(radii**f)                                   # aspect ratio
    cs = hr*(((2e30)*(6.67e-11))/(radii*1.5e11))**0.5     # [m/s]
    p = (sigma_gas_1D*(cs**2)/((2*np.pi)**0.5))*(hr**-1)*((radii*1.5e11)**-1)
    gamma = (radii/p)*np.abs(np.append(np.diff(p)/np.diff(radii), 0))
    tau_drift = ((2/(np.pi*rhodust))*(np.dot((radii*sigma_gas_azimsum/(gamma*hr*cs)).reshape(nrad,1), a.reshape(1,ndust))*1.5e11))*3.17e-14    # drift timescale in Myr

    fig,ax = plt.subplots(figsize=(9,6))
    levels = np.linspace(-6, 4, 11)                   
    print("Plotting drift timescale....")
    A,R = np.meshgrid(a,radii)
    con = ax.contourf(R, A, np.log10(tau_drift), cmap="YlGnBu", levels=levels)
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("a (cm)")
    ax.set_xlabel("R (AU)")

    fig.subplots_adjust(right=0.89, hspace=0.3)
    cbar_ax = fig.add_axes([0.91, 0.53, 0.02, 0.4])
    fig.colorbar(con, cax=cbar_ax, orientation="vertical", label="log$[\\tau_{drift} (Myr)]$")
    ax.axvline(Rp, linestyle='dashed', color='black')

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
    parser = argparse.ArgumentParser(description='Generate 1D plots', prefix_chars='-')

    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/astro/phrkvg/simulations/planet_growth/convergence_tests"],help="working directory containing simulations")
    parser.add_argument('-sim', metavar='sim', type=str, nargs=1, default=["Mp120_420x750_log"] ,help="simulation directory containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images" ,help="directory to save plots to")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-porbits', action="store_true")
    parser.add_argument('-style', metavar='style', type=str, nargs=1, default=["publication"] ,help="style sheet to apply to plots")
    parser.add_argument('-o', metavar='output',default=[1], type=int, nargs=1 ,help="output to plot")

    args = parser.parse_args()
    wd = args.wd[0]
    sim = args.sim[0]
    simdir = f"{wd}/{sim}/"
    plot_window = args.plot_window
    plots_savedir = args.savedir
    p_orbits = args.porbits
    style = args.style
    output = args.o[0]

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

    planets_data = np.genfromtxt(f"{simdir}/planet.cfg").reshape(1,6)
    Rp = planets_data[:,1][0]
    Mp = planets_data[:,2][0]
    planet_period = (Rp**3)**0.5                 # orbital period of planet in yrs

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

    gasfile = f"gasdens{output}.dat" 
    sigma_gas = np.fromfile(simdir+gasfile).reshape(nrad,nphi)/(1.125e-7)         # convert back to g/cm2
    sigma_gas_azimsum = np.sum(sigma_gas, axis=1)                                 # sum over all phi   
    sigma_gas_1D = sigma_gas_azimsum/nphi                                         # dimensions: (nrad) 

    sigma_dust = np.zeros((ndust, nrad, nphi))
    n_grains = np.zeros((ndust))

    for n in np.arange(ndust):
        dust_file = f"dustdens{n}_{output}.dat" 
        sigma_dust[n] = np.fromfile(simdir+dust_file).reshape(nrad,nphi)/(1.125e-7)
        n_grains[n] = np.sum(sigma_dust[n])*3/(a[n]*rhodust*4)                  # number of dust grains at size a and time t

    sigma_dust_azimsum = np.sum(sigma_dust, axis=2)                  # sum over all phi 
    sigma_dust_1D = sigma_dust_azimsum/nphi                          # dimensions: (ndust, nrad)
    sigma_gas_azimsum = np.sum(sigma_gas, axis=1)                    # sum over all phi   

    if not plot_window:
        for s in style:
            plt.style.use([f"../styles/{s}.mplstyle"])

    # plot_tdrift()
    # plot_tgrowth()
    # plot_tmigration()
    plot_Mp_regimes()

    if plot_window:
        plt.show()
