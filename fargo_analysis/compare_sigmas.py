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