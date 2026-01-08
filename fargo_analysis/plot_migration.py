import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D
import argparse
import re
plt.style.use('default')
plt.style.use([f"../styles/publication.mplstyle"])

cm = plt.get_cmap('viridis')
colour_cycler = cm(np.linspace(0, 1, 5))

wd = "/home/astro/phrkvg/simulations/highres_models/"
sims = ["Mp10_mig/", "Mp20_mig/", "Mp40_mig/", "Mp100_mig/", "Mp160_mig/"]
typeii_bounds = [-1, -1, 4.6, 21.2, -1]

fig, ax = plt.subplots(figsize=(7,5))

for s, sim in enumerate(sims):
    colour = colour_cycler[s]
    planetmass = re.search(r"Mp(\d+)_", sim).group(1)

    orbit_data = np.loadtxt(wd+sim+"orbit0.dat")
    orbit_data = np.unique(orbit_data, axis=0)[::10]
    times = (orbit_data[:,0]/(2*np.pi))*1e-6
    orbital_radii = orbit_data[:,2]

    ax.plot(times, orbital_radii, color=colour, label=f"{planetmass} $M_\oplus$")
    ax.axhline(typeii_bounds[s], color=colour, linestyle="dashed")


ax.set_xlim(0,0.25)
ax.set_ylim(0,42)
ax.set_xlabel("Time (Myr)")
ax.set_ylabel("Orbital radius (AU)")
ax.legend()
fig.tight_layout()
fig.savefig(f"./images/comparison_plots/migration_tracks.png", dpi=300)



