import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pylab import*
import argparse
import glob
from scipy.interpolate import RectBivariateSpline
import scipy.ndimage
from matplotlib import rcParams
import matplotlib.ticker as mtick

import warnings
warnings.filterwarnings("ignore")

plt.style.use('default')

#____________________________________PLOTTING FUNCTIONS ____________________________________#
                                                                                                      
# Plot the dust surface mass density at any one time                                                   
def plot_surf_dens(wd,simdir,dustnums,outputnumber,lin_scaling,cbmin,cbmax,plot_gas,beam,plot_planet,zoom):
# Plots a 2D graph of the dust surface density                
    filepaths = []

    params_file = f'{wd}/{simdir}/variables.par'
    params_dict = {}
    param_lines = open(params_file).readlines()
    for line in param_lines:
        if line.split():
            param_label, param_value = line.split()[0:2]
            params_dict.update([(param_label, param_value)])

    nphi = int(params_dict['NX'])
    nrad = int(params_dict['NY'])
    ndust = int(params_dict['NDUST'])
    mingsize = float(params_dict['MIN_GRAIN_SIZE'])
    maxgsize = float(params_dict['MAX_GRAIN_SIZE'])
    sigma0 = float(params_dict['SIGMA0'])
    Rc = float(params_dict['SIGMACUTOFF'])
    dt_orbits = int(float(params_dict['DT'])/(2*np.pi))   # 2pi = 1 orbit = 1 yr
    ninterm = float(params_dict['NINTERM'])               # number of dts between outputs
    dt_outputs = dt_orbits*ninterm                        # time between outputs
    time = round(outputnumber*dt_outputs*1e-6,2)          # time in Myrs
    n_size_decades = np.log10(maxgsize) - np.log10(mingsize)
    dustsizes = mingsize*10**((np.array(dustnums)/ndust)*(n_size_decades))
    dustsizes = [np.format_float_positional(d,3,fractional=False,unique=True) for d in dustsizes]

    phi = np.linspace(0,2*np.pi,nphi+1)
    r_cells = np.loadtxt(f'{wd}/{simdir}/domain_y.dat')[3:-3]             # ignore ghost cells
    radii = np.array([(r_cells[n]+r_cells[n+1])/2 for n in range(len(r_cells)-1)])

    r,theta=np.meshgrid(radii,phi)
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    a = np.logspace(np.log10(mingsize),    
                    np.log10(maxgsize),
                    ndust+1)

    a = (0.5*(a[1:] + a[:-1]))                                   # grain sizes in middles of bins (in cm)
    sum_dustvol = np.sum(a**3)

    if plot_gas in ["yes", "y"]:
        filepaths.append(f'{wd}/{simdir}/gasdens{outputnumber}.dat')

    if dustnums[0] == -1:    # use -1 to plot mass-weighted sum of multiple grain sizes
        if len(dustnums)==1:
            dust_files = glob.glob(f"{wd}/{simdir}/dustdens*_{outputnumber}.dat")
        else:
            dustnum_range = np.arange(dustnums[1], dustnums[-1])
            dust_files = [f"{wd}/{simdir}/dustdens{d}_{outputnumber}.dat" for d in dustnum_range]
            
        weighted_sum_dens = np.zeros((nrad, nphi)) 
        for i,file in enumerate(dust_files):
            surfdens = np.fromfile(file).reshape(nrad,nphi)/(1.125e-7)
            weighted_sum_dens = weighted_sum_dens + surfdens*(a[i]**3)
        avgdustdens = weighted_sum_dens/sum_dustvol
        dens_first_wedge = avgdustdens[:,0].reshape(nrad,1)
        dens_additional = np.concatenate((avgdustdens[:,:],dens_first_wedge),axis=1)
        
        if  np.isnan(cbmin):
            vmin = np.percentile(avgdustdens, 10)
        else:
            vmin = cbmin

        if  np.isnan(cbmax):
            vmax = np.percentile(avgdustdens, 90)
        else: 
            vmax = cbmax

        fig = plt.figure()
        if (lin_scaling not in ["no", "n"]):
            plt.pcolormesh(x,y,dens_additional.T,cmap=cm.Oranges_r,shading="auto",vmin=vmin, vmax=vmax)
            cbar_label = '$\Sigma$ [$g/cm^{2}$] '

        else:
            im = plt.pcolormesh(x,y,np.log10(dens_additional.T),cmap=cm.Oranges_r,shading="auto",vmin=np.log10(vmin), vmax=np.log10(vmax))
            im.set_rasterized(True)
            cbar_label = 'log $\Sigma$ [$g/cm^{2}$]'

        cbar = plt.colorbar()
        cbar.set_label(cbar_label)
        plt.xlabel('x [AU]')
        plt.ylabel('y [AU]')

        if plot_planet not in ["no", "n"]:
            with open(f"{wd}/{simdir}/planet.cfg", 'r') as f:
                num_planets = len(f.readlines()) - 5        # ignore 5 lines of headers

            for n in range(num_planets):
                planet_data = np.unique(np.loadtxt(f"{wd}/{simdir}/planet{n}.dat"), axis=0)
                xp, yp = planet_data[outputnumber][1], planet_data[outputnumber][2]
                if n == 0:
                    rp = (xp**2 + yp**2)**0.5
                    planet0_period = rp**1.5                     # orbital period of first planet in yrs
                plt.scatter([-xp], [-yp], color='g', marker='.')
        
        planet_orbits = int(round(time*1e6/planet0_period,0))

        if zoom:
            plt.xlim(-zoom,zoom)
            plt.ylim(-zoom,zoom)

        if len(dustnums)==1:
            plot_title = f"Mass-averaged dust density at t = {time}Myr = {planet_orbits} orbits"
        else:
            plot_title = f"Mass-averaged dust density ({dustnum_range[0]}-{dustnum_range[-1]}cm) at t = {time}Myr = {planet_orbits} orbits"

        plt.title(plot_title)
        plt.tight_layout()
        plt.savefig(f'./images/dustavg_{simdir}_{outputnumber}.png', dpi=150)

    else:
        for dustnum in dustnums:
            filepaths.append(f'{wd}/{simdir}/dustdens{dustnum}_{outputnumber}.dat')

    for file_i, filepath in enumerate(filepaths):
        print(f"Plotting {filepath.split('/')[-1]}...")
        surfdens = np.fromfile(filepath).reshape(nrad,nphi)/(1.125e-7)
        if  np.isnan(cbmin):
            vmin = np.percentile(surfdens, 10)
        else: 
            vmin = cbmin

        if  np.isnan(cbmax):
            vmax = np.percentile(surfdens, 90)
        else: 
            vmax = cbmax

        dens_first_wedge = surfdens[:,0].reshape(nrad,1)
        dens_additional = np.concatenate((surfdens[:,:],dens_first_wedge),axis=1)

        fig = plt.figure()

        if (lin_scaling not in ["no", "n"]):
            plt.pcolormesh(x,y,dens_additional.T,cmap=cm.Oranges_r,shading="auto",vmin=vmin, vmax=vmax)
#            plt.pcolormesh(x,y,dens_additional.T,cmap=cm.Oranges_r,vmin=cbmin,vmax=cbmax)
            cbar_label = '$\Sigma$ [$g/cm^{2}$] '

        else:
            im = plt.pcolormesh(x,y,np.log10(dens_additional.T),cmap=cm.Oranges_r,shading="auto",vmin=np.log10(vmin), vmax=np.log10(vmax))
            im.set_rasterized(True)
            cbar_label = 'log $\Sigma$ [$g/cm^{2}$]'

        # cbar = colorbar(format=ticker.FuncFormatter(fmt))
        cbar = plt.colorbar()
        cbar.set_label(cbar_label)
        plt.xlabel('x [AU]')
        plt.ylabel('y [AU]')

        if plot_planet not in ["no", "n"]:
            with open(f"{wd}/{simdir}/planet.cfg", 'r') as f:
                num_planets = len(f.readlines()) - 5        # ignore 5 lines of headers

            for n in range(num_planets):
                planet_data = np.unique(np.loadtxt(f"{wd}/{simdir}/planet{n}.dat"), axis=0)
                xp, yp = planet_data[outputnumber][1], planet_data[outputnumber][2]
                if n == 0:
                    rp = (xp**2 + yp**2)**0.5
                    planet0_period = rp**1.5                     # orbital period of first planet in yrs
                plt.scatter([-xp], [-yp], color='g', marker='.')
        
        planet_orbits = int(round(time*1e6/planet0_period,0))

        if zoom:
            plt.xlim(-zoom,zoom)
            plt.ylim(-zoom,zoom)

        plt.tight_layout()

        # r_celledge = [None]*nrad
        # ncelledge = nrad+1

        rcell = [None] * (len(radii)-1)
        for i in range(len(radii)-1):
            rcell[i] = (radii[i]+radii[i+1])/2.

        for i in range(4):
            j = int(nphi/4.*i)
            surfdensi = [None]*nrad
            for k in range(nrad):
                surfdensi[k] = surfdens[k,j]

        if (beam in ["yes", "y"]):
            # Put the image in cartesian co-ordinates
            interp = RectBivariateSpline(np.linspace(0,2*np.pi,nphi+1),rcell,dens_additional.T,kx=1,ky=1)
            npix=2048
            x_im=np.linspace(-3,3,npix)
            y_im=np.linspace(-3,3,npix)
            X_im,Y_im=np.meshgrid(x_im,y_im)
            R_im=np.sqrt(X_im**2+Y_im**2)
            Theta_im=np.arctan2(Y_im,X_im)+np.pi
            image=interp(Theta_im,R_im,grid=False)

            # Convolve with a gaussian and plot
            sigma = 50.0
            convolved=scipy.ndimage.filters.gaussian_filter(image,sigma)
            fig3 = plt.figure()
            ax3 = fig3.add_subplot(111) 
#            plt.pcolormesh(x,y,log10(dens_additional.T),cmap=cm.Oranges_r,vmin=cbmin,vmax=cbmax)
            plt.pcolormesh(-x_im,-y_im,np.log10(convolved),cmap=cm.Oranges_r,vmin=vmin,vmax=vmax)
#            plt.contourf(image, levels=10)
            ax3.set_xlim([-2.5,2.5])
            ax3.set_ylim([-2.5,2.5])
            cbar_label = 'log surface density [code units] (convolved)'
            cbar = plt.colorbar()
            cbar.set_label(cbar_label,color='black')

        if plot_gas in ["yes", "y"] and file_i == 0:
            plt.title(f"Gas density at t = {time}Myr = {planet_orbits} orbits")
            plt.savefig(f'./images/gas_{simdir}_{outputnumber}.png', dpi=150)
        else:
            plt.title(f"Dust density of {str(dustsizes[file_i-1])}cm grains at t = {time}Myr = {planet_orbits} orbits")
            plt.savefig(f'./images/dust{file_i}_{simdir}_{outputnumber}.png', dpi=150)

        # show()

if __name__ == "__main__":
        #___________________________________TAKE COMMAND LINE ARGUMENTS__________________________#     

    parser = argparse.ArgumentParser(description='Plot semi major axis', prefix_chars='-')
    parser.add_argument('-o', metavar='output_number', type=int,default=0, nargs=1 ,help="The output number for a certain orbit")
    parser.add_argument('-d', metavar='dust_number',default=[-1], type=int, nargs="*" ,help="The dust type(s) which is an integer")
    parser.add_argument('-cb', metavar='colourbar scale', type=float, nargs=2, default=[np.nan,np.nan], help="define the maximum and minimum colour bar scale")
    parser.add_argument('-lin', metavar='linear colour scale', type=str, nargs=1, default="no", help="linear plot")
    parser.add_argument('-gas', metavar='gas density', type=str, nargs=1, default="no", help="plot gas density")
    parser.add_argument('-con', metavar='convolved', type=str, nargs=1, default="no", help="convolved with a gaussian beam")
    parser.add_argument('-wd', metavar='wd', type=str, nargs=1, default=["/home/astro/phrkvg/simulations/planet_growth/lowres_models"],help="working directory containing simulations")
    parser.add_argument('-sim', metavar='sim', type=str, nargs=1, default=["planettest"] ,help="simulation directory containing output files")
    parser.add_argument('-planet', metavar='plot planet', type=str, nargs=1, default="no", help="plot planet")
    parser.add_argument('-zoom', metavar='zoomed plot', type=float, nargs=1, default=[0], help="set lims to zoom to")
    parser.add_argument('-plot_window', action="store_true")
    parser.add_argument('-style', metavar='style', type=str, nargs=1, default=["publication"] ,help="style sheet to apply to plots")

    args = parser.parse_args()

    DUSTNUMS = args.d
    print("dust bins = ", DUSTNUMS)

    o = args.o[0]

    cbmin = args.cb[0]
    cbmax = args.cb[1]

    if (args.lin == 0.):
        lin_scale = args.lin
    else:
        lin_scale = args.lin[0]
    
    print("lin_scale = ",lin_scale)

    if (args.gas == 0.):
        gasdens = args.gas
    else:
        gasdens = args.gas[0]

    print("plot_gas = ",gasdens)

    if (args.con == 0.):
        convolve = args.con
    else:
        convolve = args.con[0]

    wd = args.wd[0]
    simdir = args.sim[0]
    plot_planet = args.planet[0]
    zoom = args.zoom[0]
    plot_window = args.plot_window
    style = args.style

    if plot_window:
        matplotlib.use('TkAgg')

        #________________________________END COMMAND LINE ARGS__________________________________#      

    if not plot_window:
        for s in style:
            plt.style.use([f"../styles/{s}.mplstyle"])

    plot_surf_dens(wd,simdir,DUSTNUMS,o,lin_scale,cbmin,cbmax,gasdens,convolve,plot_planet,zoom)

    if plot_window:
        plt.show()
