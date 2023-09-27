import numpy as np
import matplotlib.pyplot as plt
import matplotlib   


outputs = [0,100, 200]
sim="Mp200_vlimBC"
simdir = f"/home/astro/phrkvg/simulations/planet_growth/lowres_models/{sim}"
nrad = 250
nphi = 50

gasvy = np.zeros((len(outputs), nrad, nphi))
radii = np.linspace(5,150,nrad)

for i,t in enumerate(outputs):
    gasvy[i] = np.fromfile(f"{simdir}/gasvy{i}.dat").reshape(nrad,nphi)

# gasvy[gasvy>=0] = np.nan

fig, ax = plt.subplots(1, len(outputs), tight_layout=True)

for i,o in enumerate(outputs):
    print('Innermost annulus: ', gasvy[i,0,:])
    for n in range(nphi):
        ax[i].scatter(radii, gasvy[i,:,n], c='b', alpha=0.2)
    # Min value in innermost annulus i.e. at boundary
    ax[i].set_title(f"output: {o} \n Min: {round(np.nanmin(gasvy[i,0,:]),8)}")
    ax[i].set_xlim(4.75,7)
    # ax[i].axhline(-6e-7, c='k')
    
    # vr = np.nansum(gasvy[i], axis=1)/nphi
    # ax[i].plot(radii, vr)

fig.savefig(f"./images/Vr_{sim}.png")
