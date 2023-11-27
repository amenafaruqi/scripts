import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
plt.style.use('default')


def plot_Mp_regimes():  #figure 2 from dipierro 2017 i.e. as function of Mp and a
    pass


def plot_tdrift(): #as function of R and a
    hr = hr0*(radii**f)                                   # aspect ratio
    cs = hr*(((2e30)*(6.67e-11))/(radii*1.5e11))**0.5     # [m/s]
    p = (sigma_gas_1D*(cs**2)/((2*np.pi)**0.5))*(hr**-1)*((radii*1.5e11)**-1)
    gamma = (radii/p)*np.abs(np.append(np.diff(p)/np.diff(radii), 0))
    t_drift = ((2/(np.pi*rhodust))*(np.dot((radii*sigma_gas_azimsum/(gamma*hr*cs)).reshape(nrad,1), a.reshape(1,ndust))*1.5e11))*3.17e-14    # drift timescale in Myr

    fig,ax = plt.subplots(figsize=(16,12))
    levels = np.linspace(-6, 4, 11)                   
    print("Plotting drift timescale....")
    A,R = np.meshgrid(a,radii)
    con = ax.contourf(R, A, np.log10(t_drift), cmap="Greys", levels=levels)
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("a (cm)")
    ax.set_xlabel("R (AU)")

    fig.subplots_adjust(right=0.89, hspace=0.3)
    cbar_ax = fig.add_axes([0.91, 0.53, 0.02, 0.4])
    fig.colorbar(con, cax=cbar_ax, orientation="vertical", label="log$[\\tau_{drift} (Myr)]$")

    fig.savefig(f"{plots_savedir}/{sim}_tdrift.png")



def plot_tmigration():  # as function of R
    alpha = -1
    beta = -0.5
    gamma = 5/3
    zeta = beta - (gamma-1)*alpha
    Rs = np.arange(10,101,5)
    Ms = np.arange(10,201,10)*3e-6                   # in units of Msun
    tauI_arr_R = np.zeros(len(Rs))
    tauI_arr_M = np.zeros(len(Ms))
    for i,R in enumerate(Rs):
        Omega = np.sqrt(4*(np.pi**2)/(R**3))
        h = hr0*(R**f)                                   # aspect ratio H/R at Rp
        Sigmap = sigma0*(float(R)**-sigmaslope)*np.exp(-R/Rc)
        Gamma = (-2.5-(1.7*beta)+(alpha/10)+(1.1*(1.5-alpha))+(7.9*zeta/gamma))*((Mp/h)**2)*Sigmap*(R**4)*(Omega**2)/gamma
        tauI_arr_R[i] = (R**2)*Omega*Mp/(2*Gamma)*1e-6    # convert to Myr

    for i,M in enumerate(Ms):
        Omega = np.sqrt(4*(np.pi**2)/(Rp**3))
        h = hr0*(Rp**f)                                   # aspect ratio H/R at Rp
        Sigmap = sigma0*(float(Rp)**-sigmaslope)*np.exp(-Rp/Rc)
        Gamma = (-2.5-(1.7*beta)+(alpha/10)+(1.1*(1.5-alpha))+(7.9*zeta/gamma))*((M/h)**2)*Sigmap*(Rp**4)*(Omega**2)/gamma
        tauI_arr_M[i] = (Rp**2)*Omega*M/(2*Gamma)*1e-6

    fig2, ax2 = plt.subplots(nrows=2, figsize=(7,6))
    ax2[0].plot(Rs,tauI_arr_R)
    ax2[0].set_xlabel("R (AU)")
    ax2[0].set_ylabel("$\\tau_{I}$ (Myr)")
    ax2[0].set_xlim(np.min(Rs), np.max(Rs))
    ax2[0].axvline(Rp, linestyle='dashed', color='black')

    ax2[1].plot(Ms/3e-6,tauI_arr_M)
    ax2[1].set_xlabel("$M (M_\oplus)$")
    ax2[1].set_ylabel("$\\tau_{I}$ (Myr)")
    ax2[1].set_xlim(np.min(Ms/3e-6), np.max(Ms/3e-6))
    ax2[1].axvline(Mp/3e-6, linestyle='dashed', color='black')

    fig2.tight_layout()
    fig2.savefig(f"{plots_savedir}/{sim}_tmig.png")



def plot_tgrowth():  #as function of R and a
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate 1D plots', prefix_chars='-')

    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/astro/phrkvg/simulations/planet_growth/convergence_tests"],help="working directory containing simulations")
    parser.add_argument('-sim', metavar='sim', type=str, nargs=1, default=["Mp120_420x750_log"] ,help="simulation directory containing output files")
    parser.add_argument('-savedir', metavar='savedir', type=str, nargs=1, default="./images" ,help="directory to save plots to")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-porbits', action="store_true")
    parser.add_argument('-style', metavar='style', type=str, nargs=1, default=["publication"] ,help="style sheet to apply to plots")

    args = parser.parse_args()
    wd = args.wd[0]
    sim = args.sim[0]
    simdir = f"{wd}/{sim}/"
    plot_window = args.plot_window
    plots_savedir = args.savedir
    p_orbits = args.porbits
    style = args.style

    if plot_window:
        matplotlib.use('TkAgg') 

    params_file = f'{simdir}/variables.par'
    params_dict = {}

    param_lines = open(params_file).readlines()
    for line in param_lines:
        if line.split():
            param_label, param_value = line.split()[0:2]
            params_dict.update([(param_label, param_value)])

    nphi = int(params_dict['NX'])
    nrad = int(params_dict['NY'])
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

    gasfile = "gasdens0.dat" 
    sigma_gas = np.fromfile(simdir+gasfile).reshape(nrad,nphi)/(1.125e-7)         # convert back to g/cm2
    sigma_gas_azimsum = np.sum(sigma_gas, axis=1)                                 # sum over all phi   
    sigma_gas_1D = sigma_gas_azimsum/nphi                                         # dimensions: (noutputs, nrad) 


    if not plot_window:
        for s in style:
            plt.style.use([f"../styles/{s}.mplstyle"])

    plot_tdrift()
    plot_tmigration()

