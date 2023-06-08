from numpy import*
import matplotlib.pyplot as plt
import matplotlib as plt
import matplotlib.ticker as ticker
from pylab import*
import argparse
import sys
import os
import math as mt
from os import path
#from plot_section import create_simulation
#from draggable_colorbar import DraggableColorbar
#import mynormalize
from scipy.interpolate import RectBivariateSpline
import scipy.ndimage
#import visual.graph as vg
from matplotlib import rcParams
import matplotlib.ticker as mtick

plt.style.use('default')
plt.style.use(['../styles/presentation.mplstyle', '../styles/darkbg.mplstyle'])

#____________________________________PLOTTING FUNCTIONS ____________________________________#

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
                                                                                                      
# Plot the dust surface mass density at any one time                                                   
def plot_surf_dens(simdir,dustnums,outputnumber,plottype,lin_scaling,cbmin,cbmax,plot_gas,beam,plot_planet):
# Plots a 2D graph of the dust surface density                
    filepaths = []

    if plot_gas in ["yes", "y"]:
        filepaths.append(f'/home/astro/phrkvg/simulations/planet_growth/{simdir}/gasdens{outputnumber}.dat')

    for dustnum in dustnums:
        filepaths.append(f'/home/astro/phrkvg/simulations/planet_growth/{simdir}/dustdens{dustnum}_{outputnumber}.dat')

    params_file = f'/home/astro/phrkvg/simulations/planet_growth/{simdir}/variables.par'
    params_dict = {}
    param_lines = open(params_file).readlines()
    for line in param_lines:
        if line.split():
            param_label, param_value = line.split()[0:2]
            params_dict.update([(param_label, param_value)])

    nphi = int(params_dict['NX'])
    nrad = int(params_dict['NY'])
    ndust = int(params_dict['NDUST'])

    for file_i, filepath in enumerate(filepaths):
        print(f"Plotting {filepath.split('/')[-1]}...")
        surfdens = fromfile(filepath).reshape(nrad,nphi)/(1.125e-7)
        if  np.isnan(cbmin):
            vmin = np.percentile(surfdens, 10)
        else: 
            vmin = cbmin

        if  np.isnan(cbmax):
            vmax = np.percentile(surfdens, 90)
        else: 
            vmax = cbmax

    #    nx_array = 
    #    ny_array = 
    #    r = 0.5+2.5/ny*ny_array
    #    phi = 2.*pi/nx*nx_array
    #    if (plottype == "contour"):

    #    imshow(log(surfdens),origin='lower',cmap=cm.Oranges_r,aspect='auto',vmin=cbmin,vmax=cbmax)
    #    cbar = colorbar()
    #    cbar.set_label('dust surface density')
    #    show()

        if (plottype == "pcontour"):
            phi = linspace(0,2*pi,nphi+1)
            radii=loadtxt(f'/home/astro/phrkvg/simulations/planet_growth/{simdir}/used_rad.dat')
            r,theta=meshgrid(radii,phi)
            
            dens_first_wedge = surfdens[:,0].reshape(nrad,1)
            # print(dens_first_wedge.shape)
            dens_additional = concatenate((surfdens[:,:],dens_first_wedge),axis=1)

            x = r*cos(theta)
            y = r*sin(theta)

            fig = plt.figure()

    #        plt.rcParams['axes.facecolor'] = 'none'

            if (lin_scaling not in ["no", "n"]):
                plt.pcolormesh(x,y,dens_additional.T,cmap=cm.Oranges_r,shading="auto",vmin=vmin, vmax=vmax)
    #            plt.pcolormesh(x,y,dens_additional.T,cmap=cm.Oranges_r,vmin=cbmin,vmax=cbmax)
                cbar_label = '$\Sigma$ [$g/cm^{2}$] '

            else:
                im = plt.pcolormesh(x,y,log10(dens_additional.T),cmap=cm.Oranges_r,shading="auto",vmin=log10(vmin), vmax=log10(vmax))
                im.set_rasterized(True)
                cbar_label = 'log $\Sigma$ [$g/cm^{2}$]'

            # cbar = colorbar(format=ticker.FuncFormatter(fmt))
            cbar = colorbar()
            cbar.set_label(cbar_label)
            plt.xlabel('x [AU]')
            plt.ylabel('y [AU]')

            if plot_planet not in ["no", "n"]:
                planet_data = loadtxt(f"/home/astro/phrkvg/simulations/planet_growth/{simdir}/planet0.dat")
                xp, yp = planet_data[outputnumber][1], planet_data[outputnumber][2]
                plt.scatter([-xp], [-yp], color='g', marker='.')
            # plt.xlim(-50,50)
            # plt.ylim(-50,50)
            plt.tight_layout()

    #        cbar.ax.yaxis.set_tick_params(color='white')
    #        cbar.ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
            # plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'))#, color='w')
    #        cbar.ax.set_xticklabels(color='white')
    #        cbar.ax.get_ticklabels('white')
    #        spines['top'].set_color('red')
    #        ax1.spines['top'].set_color('red')
    #        ax1.set_xlabel('X-axis')
    #        ax1.xaxis.label.set_color('red')
            

    #        ax1 = plt.gca()
            # ax = fig.add_subplot(111)
            # plt.axis('scaled')
            # ax.set_xlabel('x [code units]')
            # ax.set_ylabel('y [code units]')
            # ax.set_xlim([-1.3,1.3])
            # ax.set_ylim([-1.3,1.3])
    #        ax.tick_params(axis='y', which='both', color='w')
    #        ax.xaxis.label.set_color('white')
    #        ax.yaxis.label.set_color('white')
    #        ax.xticks.set_color('red')
    #        ax.tick_params(axis='x', colors='red')

    #        for spine in ax.axes.spines.values():
    #            spine.set_edgecolor('white')

    #        [t.set_color('white') for t in ax.xaxis.get_ticklines()]
    #        [t.set_color('white') for t in ax.yaxis.get_ticklines()]


    #        [t.set_color('white') for t in ax.xaxis.get_ticklabels()]
    #        [t.set_color('white') for t in ax.yaxis.get_ticklabels()]
    #        [t.set_color('white') for t in cbar.set_ticklabels()]
    #        cbar.ax.label.set_color('white')
    #        [t.set_color('white') for t in cbar.ax.yaxis.set_tick_params()]
    #        cbar.ax.yaxis.set_tick_params(color='white')
    #        cbar.tick(color='white')
    #        cbar.outline.set_edgecolor('white')
    #        cbar.color = 'white'

    #        for tick in cbar.get_ticklabels():
    #            tick.set_color('red')

    # Uncomment the following text, as well as any future references to ax2, if a second plot is needed with the azimuthal averaged sigmad and four slices of the disc
    #        ax2 = fig.add_subplot(212)
    #        plt.figure()
    #        ax2 = plt.gca()
    #        ax2.set_xlabel('r')
    #        ax2.set_ylabel("dust surface density")
    #        sigmad_tot = sum(surfdens,axis=1)
    #        sigmad_azi = sigmad_tot/nphi


            r_celledge = [None]*nrad
            ncelledge = nrad+1
    #        print type(ncelledge), ncelledge

    #        filepath_r = 'link/used_rad.dat'
    #        f = open(filepath_r,'r')
    #        lines = f.readlines()
    #        print type(lines), lines

    #        for line in xrange(ncelledge-1):
    #            data = line.split()
    #            r_celledge.append(float(data[0]))
    #
    #            r_celledge.append(

            rcell = [None] * (len(radii)-1)
            for i in range(len(radii)-1):
    #            rcell[i] = (radii[i]*radii[i+1])**0.5
                rcell[i] = (radii[i]+radii[i+1])/2.

    #        print rcell
    #        print len(rcell), len(sigmad_azi)

    #        ax2.plot(rcell,sigmad_azi)
    
            for i in range(4):
                j = int(nphi/4.*i)
    #            print nx/4.*i, j
                surfdensi = [None]*nrad
                for k in range(nrad):
    #                print "inside k = ",k
    #                print surfdens
    #                surfdens = []
    #                print type(surfdens)                
    #                print surfdens[0,0], type(surfdens[0,0])
    #                sys.exit()
    #                surfdensi = float(surfdens[k,j])
    #                print surfdensi
                    surfdensi[k] = surfdens[k,j]
    #            print len(rcell), len(surfdensi)
    #            ax2.plot(rcell,surfdensi)

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
                print(image)
                print(type(image), len(image), shape(image), shape(x), shape(y), shape(x_im), shape(y_im))
                plt.pcolormesh(-x_im,-y_im,log10(convolved),cmap=cm.Oranges_r,vmin=vmin,vmax=vmax)
    #            plt.contourf(image, levels=10)
                ax3.set_xlim([-2.5,2.5])
                ax3.set_ylim([-2.5,2.5])
                cbar_label = 'log surface density [code units] (convolved)'
                cbar = colorbar()
                cbar.set_label(cbar_label,color='black')


        else:
            plt.figure()
            imshow(log10(surfdens),origin='lower',cmap=cm.Oranges,shading="auto",aspect='auto',vmin=vmin,vmax=vmax)
            cbar = colorbar()
            cbar.set_label('surface density') 


    #    set_axis_bgcolor('red')
        if plot_gas in ["yes", "y"] and file_i == 0:
            plt.savefig(f'./images/gas_{simdir}_{outputnumber}.png', dpi=150)
        else:
            plt.savefig(f'./images/dust{file_i}_{simdir}_{outputnumber}.png', dpi=150)

        # show()

if __name__ == "__main__":
        #___________________________________TAKE COMMAND LINE ARGUMENTS__________________________#     

    parser = argparse.ArgumentParser(description='Plot semi major axis', prefix_chars='-')
    parser.add_argument('-o', metavar='out_put_number', type=int,default=-1, nargs=1 ,help="The output\
 number for a certain orbit")
    parser.add_argument('-d', metavar='dust_number',default=[0], type=int, nargs="*" ,help="The dust type(s) \
which is an integer")
    parser.add_argument('-cb', metavar='colourbar scale', type=float, nargs=2, default=[np.nan,np.nan], help="define the maximum and minimum colour bar scale")
    parser.add_argument('-ptype', metavar='plot type', type=str, nargs=1, default="pcontour", help="type of plot")
    parser.add_argument('-lin', metavar='linear colour scale', type=str, nargs=1, default="no", help="linear plot")
    parser.add_argument('-gas', metavar='gas density', type=str, nargs=1, default="no", help="plot gas density")
    parser.add_argument('-con', metavar='convolved', type=str, nargs=1, default="no", help="convolved with a gaussian beam")
    parser.add_argument('-dir', metavar='dir', type=str, nargs=1, default="planettest" ,help="simulation directory containing output files")
    parser.add_argument('-planet', metavar='plot planet', type=str, nargs=1, default="no", help="plot planet")

    args = parser.parse_args()

    DUSTNUMS = args.d
    print(DUSTNUMS)

    Out = args.o[0]

    cbmin = args.cb[0]
    cbmax = args.cb[1]

    if (args.ptype == 0.):
        plot_type = args.ptype
    else:
        plot_type = args.ptype[0]

    print("plot type = ",plot_type)

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

    simdir = args.dir[0]
    plot_planet = args.planet[0]

    print(simdir)
        #________________________________END COMMAND LINE ARGS__________________________________#      

    # Change matplotlib.rcParams
#    rcParams["text.color"] = 'white'
#    rcParams["text.fontsize"] = 14

    plot_surf_dens(simdir,DUSTNUMS,Out,plot_type,lin_scale,cbmin,cbmax,gasdens,convolve,plot_planet)
