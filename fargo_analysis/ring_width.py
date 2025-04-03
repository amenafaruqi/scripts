#calculates ring position and ring width 

import sys, os
sys.path.append("/home/astro/phrkvg/codes/dusty_fargo3d/python")

from fargo_lib import FARGO_Sim
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
from scipy.signal import peak_widths
from scipy.signal import find_peaks
from numpy.polynomial import Polynomial

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

Mearth = 3e-6
#2.5029077e-6 # from q in 1.2Msol to earth masses
maxearthmass = 40 # what all the radial profiles will be normalized to

# --------------------------------------------------------

def findgapbottom(profile, profileinit, planetloc, usedrads):
    '''
    Find the value of the surface density of the bottom of a gap, and its orbital radius.
    profile is a np array; can be gas or dust surface density profile at some time 
    *** easiest if it's normalized to the initial array, so that natural decrease of profile doesn't interfere
        hence profileinit.
    planetloc is a float (roughly around 1 in code units); of the planet's radial location
    usedrads is a np array; a list outputted by fargo of the cells' radial values, should be same length as profile
    We want to search less than half a code unit exterior to the planet and ID the minimum.
    Returns radial loc, sigma
    '''
    radseg = usedrads[np.where((usedrads>planetloc-0.01) & (usedrads<planetloc+0.4))]

    sigseg = profile[np.where((usedrads>planetloc-0.01) & (usedrads<planetloc+0.4))]
    siginitseg = profileinit[np.where((usedrads>planetloc-0.01) & (usedrads<planetloc+0.4))]

    # Find the normalized surface density minimum value
    gapsignorm = np.min(sigseg/siginitseg)
    # Find the radius where that occurs
    gaprad = radseg[np.where(sigseg/siginitseg == gapsignorm)]
    # Find the non-normalized surface density value where *that* occurs
    gapsig = sigseg[np.where(radseg == gaprad)]

    print('Sigma at gap, gap radius = ', gapsig, gaprad)

    return gaprad, gapsig 


def findringpeak(profile, profileinit, planetloc, usedrads):
    '''
    Same as findgapbottom, but for the ring peak.
    '''
    radseg = usedrads[np.where((usedrads>planetloc-0.01) & (usedrads<planetloc+0.4))]

    sigseg = profile[np.where((usedrads>planetloc-0.01) & (usedrads<planetloc+0.4))]
    siginitseg = profileinit[np.where((usedrads>planetloc-0.01) & (usedrads<planetloc+0.4))]

    # Find the normalized surface density maximum value
    ringsignorm = np.max(sigseg/siginitseg)
    # Find the radius where that occurs
    ringrad = radseg[np.where(sigseg/siginitseg == ringsignorm)]
    # Find the non-normalized surface density value where *that* occurs
    ringsig = sigseg[np.where(radseg == ringrad)]

    #print('Sigma at ring, ring radius = ', ringsig, ringrad)

    return ringrad, ringsig #convenient order for plotting


def gaussian(x, A, mean, width, base):
    '''
    gaussian to be used when finding ring widths - written by Joshua Adams, 2024
    '''
    return A*np.exp(-(x-mean)**2/(2*width**2)) + base
    
def disttocell(x, Nr, Ymin, Ymax):
    '''
    function that converts code unit distance to number of cells from centre - JA 2024
    '''
    B = (Ymax/Ymin)**(1/(Nr-1))
    A = Ymin/B
    return round(np.log10(x/A)/np.log10(B))

def celltodist(x, Nr, Ymin, Ymax):
    A = Ymin
    B = np.log(Ymax/Ymin)/Nr
    return A*np.exp(B*x)

def distcalc(xa, xb, Nr, Ymin, Ymax):
    return abs(celltodist(xa,Nr,Ymin,Ymax) - celltodist(xb,Nr,Ymin,Ymax))

def nearestflats(ydat,ringindex): # finds flat points nearest to ring peak
    zeroes = []
    idx = []
    slope = np.diff(ydat)
    for i in range(len(slope)-1):
        if slope[i]*slope[i+1] <= 0:
            zeroes.append(i)
    for i in zeroes:
        if np.abs(i-ringindex) < 5:
            zeroes.remove(i)
    lowzeroes = [x for x in zeroes if x < ringindex]
    upzeroes = [x for x in zeroes if x > ringindex]
    lowid = (np.abs(lowzeroes-ringindex)).argmin()
    upid = (np.abs(upzeroes-ringindex)).argmin()
    lowid = lowzeroes[lowid]
    upid = upzeroes[upid]
    print(lowid,upid)
    return lowid, upid

# --------------------------------------------------------

label_size = 20

fig = plt.figure(figsize=(18,30))#7,20))

ax1 = fig.add_subplot(521)
ax1.tick_params(direction='in',which='both',labelbottom=True,top=True,bottom=True,right=True,left=True,labelsize=label_size)
div1 = make_axes_locatable(ax1)
ax1.ticklabel_format(axis='y')#, scilimits=(0,1))
ax1.set_ylabel(r"Ring Width (AU)", fontsize=30)
ax1.set_xlabel(r"Planet Mass ($M_{\text{earth}}$)", fontsize=30)
'''
ax2 = fig.add_subplot(522)
ax2.tick_params(direction='in',which='both',labelbottom=True,top=True,bottom=True,right=True,left=True,labelsize=label_size)
div2 = make_axes_locatable(ax2)
ax2.ticklabel_format(axis='y')#, scilimits=(0,1))
ax2.set_xlabel(r"Planet Mass ($M_{\text{earth}}$)", fontsize=30)
'''
ax3 = fig.add_subplot(523)
ax3.tick_params(direction='in',which='both',labelbottom=True,top=True,bottom=True,right=True,left=True,labelsize=label_size)
div3 = make_axes_locatable(ax3)
ax3.ticklabel_format(axis='y')#, scilimits=(0,0))
ax3.set_ylabel(r"Ring Width", fontsize=30)
'''
ax4 = fig.add_subplot(524)
ax4.tick_params(direction='in',which='both',labelbottom=True,top=True,bottom=True,right=True,left=True,labelsize=label_size)
div4 = make_axes_locatable(ax4)
ax4.ticklabel_format(axis='y')#, scilimits=(0,0))

ax5 = fig.add_subplot(525)
ax5.tick_params(direction='in',which='both',top=True,bottom=True,right=True,left=True,labelsize=label_size)
div5 = make_axes_locatable(ax5)
ax5.ticklabel_format(axis='y')#, scilimits=(0,0))
ax5.set_ylabel(r"Ring Width", fontsize=30)
ax5.set_xlabel(r"Planet Mass ($M_{\text{earth}}$)", fontsize=30)

ax6 = fig.add_subplot(526)
ax6.tick_params(direction='in',which='both',top=True,bottom=True,right=True,left=True,labelsize=label_size)
div6 = make_axes_locatable(ax6)
ax6.ticklabel_format(axis='y')#, scilimits=(0,0))

ax7 = fig.add_subplot(527)
ax7.tick_params(direction='in',which='both',labelbottom=True,top=True,bottom=True,right=True,left=True,labelsize=label_size)
div7 = make_axes_locatable(ax7)
ax7.ticklabel_format(axis='y')#, scilimits=(0,0))
ax7.set_ylabel(r"$r_{p} - r_{r} (r_{H})$", fontsize=30)

ax8 = fig.add_subplot(528)
ax8.tick_params(direction='in',which='both',labelbottom=True,top=True,bottom=True,right=True,left=True,labelsize=label_size)
div8 = make_axes_locatable(ax8)
ax8.ticklabel_format(axis='y')#, scilimits=(0,0))

ax9 = fig.add_subplot(529)
ax9.tick_params(direction='in',which='both',top=True,bottom=True,right=True,left=True,labelsize=label_size)
div9 = make_axes_locatable(ax9)
ax9.ticklabel_format(axis='y')#, scilimits=(0,0))
ax9.set_ylabel(r"$r_{p} - r_{r} (r_{H})$", fontsize=30)
ax9.set_xlabel(r"planet mass ($M_{\text{earth}}$)", fontsize=30)

ax10 = fig.add_subplot(5,2,10)
ax10.tick_params(direction='in',which='both',top=True,bottom=True,right=True,left=True,labelsize=label_size)
div10 = make_axes_locatable(ax10)
ax10.ticklabel_format(axis='y')#, scilimits=(0,0))
ax10.set_xlabel(r"planet mass ($M_{\text{earth}}$)", fontsize=30)
'''
# --------------------------------------------------------

def maincode(simnames,alphavalue,ax,xticks_list,stokesindex):
    #savedata = np.zeros((len(simnames), 26)) # columns: Mp, R_hill, R_ring 
    #Rgap, Sgap, Rring, Sring * 6 species + h; rows = number of planet masses simulated for this alpha
    #print(savedata.shape)
    
    peak_planet_distance = [] 
    peak_planet_distance_au = []
    planet_masses = []
    ring_widths = []
    for row,simname in enumerate(simnames): 
        print("We are on: ", simname)
        # try:
        path = '/home/astro/phrkvg/simulations/dusty_fargo_models/' +simname #sys.argv[1]+'/'+simname
        Sim = FARGO_Sim(path)
        # except IndexError:
        #     Sim = FARGO_Sim('outputs/fargo_multifluid')
        # try:
        #     orb_num = int(sys.argv[1])#int(sys.argv[2]) # give the final snapshot in terminal
        # except IndexError:
        #     print("What time?")
        #     orb_num = 0

        planet = pd.read_csv(path+'/planet0.dat', usecols=[0,1,2,3,4,5,6,7], names=['out','x','y','z','vx','vy','vz','m'], sep='\t')
        planet['r'] = np.sqrt(planet['x']**2 + planet['y']**2 + planet['z']**2)
        planet['ny'] = round(planet['y']*(753.0/math.pi)+753.0)
        usedrads = np.loadtxt(path+'/domain_y.dat')
        #print('These are usedrads: ', usedrads)
        
        #savedata[row][0] = planet['m'].iloc[-1]/Mearth
        print("planet mass: ", planet['m'].iloc[-1]/Mearth) 
        #hill radii 
        r_hill = np.cbrt(planet['m'].iloc[-1]/3)
        #print("r_hill: ", r_hill) #Mearth = mass of earth / 1.2 * mass of sun
        # rhill = cbrt (Mearth * mass planet in earth masses /3) 
        
        #find ring position in plot this snap 

        r   = Sim.Rmed
        phi = Sim.phimed 
        maxR = Sim.get_parameter('Ymax')
        minR = Sim.get_parameter('Ymin')
        # frac =  Sim.get_parameter('NINTERM') / ((2*np.pi)/Sim.get_parameter('DT')) # eg. 19/20 orbits
        Ninterm = Sim.get_parameter('NINTERM')
        snap_num = int(orb_num/Ninterm)
        hvalue = Sim.get_parameter('ASPECTRATIO')
        Nr = Sim.get_parameter('NY')
        Nphi = Sim.get_parameter('NX')
        #savedata[row][5] = hvalue
        
        planet['d'] = disttocell(planet['x'],Nr,minR,maxR) # location of the planet
        ringsearch = 3*r_hill
        lowersearch = int(disttocell(planet['x'].iloc[-1] - ringsearch, 1300, 0.2, 3.0)) # the lower bound for the search (cell #)
        uppersearch = int(disttocell(planet['x'].iloc[-1] + 2.3*ringsearch, 1300, 0.2, 3.0)) # the upper bound for the search (cell #)
        print('lowerbound: ', lowersearch, ' upper bound: ', uppersearch)


        ##################
        # Get the Stoke's Number
        St = 1. / np.genfromtxt(path+'/drag_coefficients.dat')

        for i in [0,1,2]: # Josh Adams added this for loop to iterate multiple St
            #trying one stokes at a time 
            stokes = stokesindex#0#2 # og stokes = 0 (1,2,3,4)
	    #for num in [0,1,2,3,4]:
	     #   print("St are: ",St[num])
            print('Current St: ',St[stokes])

            #rho_d  = Sim.load_field('dust{}dens'.format(stokes+1), snap_num).transpose()
            #rho_d0  = Sim.load_field('dust{}dens'.format(stokes+1), 0).transpose()
            
            rho_d = np.fromfile(path+f"/dustdens{i}_600.dat").reshape(Nr,Nphi)
            rho_d0 = np.fromfile(path+f"/dustdens{i}_0.dat").reshape(Nr,Nphi)

            rho_d0 = np.average(rho_d0, axis=1)
            rho_d = np.average(rho_d, axis=1)
                
            rho_d = rho_d / rho_d0
            
            #make a gaussian fit
            xdat = np.linspace(1,Nr,Nr)
            ydat = rho_d
            
            # scipy fitting
            param, dud = curve_fit(lambda t,a,b: a*np.exp(b*t), xdat[150:], ydat[150:],p0=(40, -0.005))
            afit, bfit = param
            
            # flat is flat i guess im not sure man
            sflat = np.subtract(ydat,afit*np.exp(bfit*xdat))
            print(sflat)
            print(lowersearch, uppersearch)
            ring = np.argmax(sflat[lowersearch:uppersearch])+lowersearch
            print('ring peak:', ring)
            # gaussian fit to derivative cuts
            zeroes = nearestflats(sflat,ring)
            try:
                popt, pcov = curve_fit(gaussian,xdat[zeroes[0]:zeroes[1]],sflat[zeroes[0]:zeroes[1]],p0=[max(sflat[zeroes[0]:zeroes[1]]),np.argmax(sflat[zeroes[0]:zeroes[1]])+zeroes[0],5,1])
                Afit, meanfit, widthfit, basefit = popt
                y_fit = gaussian(xdat,Afit,meanfit,widthfit,basefit)

                # NEW WIDTH FINDER
                peaks3, _ = find_peaks(y_fit)
                print(peaks3)
                if int(max(peaks3)) < lowersearch or int(max(peaks3)) > uppersearch:
                    sfwhm = 0.0
                    dist = 0.0
                    print('there is NONE ring >:C')
                else:
                    goodpeak = [max([x for x in peaks3 if x > lowersearch and x < uppersearch])]
                    widthdata = peak_widths(y_fit, goodpeak, rel_height=0.5)
                    sfwhm = widthdata[0]
                    dist = distcalc(widthdata[3],widthdata[2],Nr,0.2,3.0)
                    print('Half-max width of ring: (cells):', sfwhm, ' (au): ',dist)
            
            except RuntimeError:
                print("Couldn't find appropriate Gaussian fit")
                sfwhm = 0.0

	        #find ring position and distance to planet
            ringrad, ringsig = findringpeak(rho_d, rho_d0, planet['r'].iloc[snap_num], usedrads)
            print("ring,planet distance: ") #using function ring position  
            print("in au: ", (ringrad - 1))
            print("in hill radii: ", (ringrad - 1)/r_hill)
            print("")

	        #data for plotting 
            planet_masses.append(int(planet['m'].iloc[-1]/Mearth))
            peak_planet_distance_au.append(ringrad - 1) 
            peak_planet_distance.append((ringrad - 1)/r_hill)
            ring_widths.append(dist) 
            stokesindex += 1

            # this statement takes the code back to the first stokes number for the next mass
            if stokesindex == len(St):
                stokesindex = 0

    #plot ring position vs ring peak density 
    #ax.scatter(ringrad, ringsig, color='k', facecolor='None', zorder=1000)
    
    #next: x axis: Mass, y axis: ring distance to planet 
    #ax.scatter(planet_masses, peak_planet_distance)

    # make scatter of the ring widths
    #ax.scatter(planet_masses, ring_widths)
    
    # divide the ring widths by stokes number
    #print(ring_widths) 
    stokes1 = ring_widths[0::5]
    stokes2 = ring_widths[1::5]
    stokes3 = ring_widths[2::5]
    stokes4 = ring_widths[3::5]
    stokes5 = ring_widths[4::5]
    masses = planet_masses[0::5]
    masses1 = masses[:]
    masses2 = masses[:]
    masses3 = masses[:]
    masses4 = masses[:]
    masses5 = masses[:]
    print(ring_widths)

    # whoever sees this: im tired and not super good at coding - please dont judge me (JA 2024)
    
    # remove zeroes
    for s in range(len(stokes1)-1,-1,-1):
        if stokes1[s] == 0:
            stokes1.pop(s)
            masses1.pop(s)
        if stokes2[s] == 0:
            stokes2.pop(s)
            masses2.pop(s)
        if stokes3[s] == 0:
            stokes3.pop(s)
            masses3.pop(s)
        if stokes4[s] == 0:
            stokes4.pop(s)
            masses4.pop(s)
        if stokes5[s] == 0:
            stokes5.pop(s)
            masses5.pop(s)

    # make fits!
    print(ring_widths)
    if len(masses1) > 0:
        print(stokes1)
        m1, b1 = np.polyfit(masses1,stokes1,deg=1)
        ax.scatter(masses1, stokes1)
        ax.plot(masses1, m1*masses1+b1,label="St = 0.2")
        ax.legend(loc="upper right")
    if len(masses2) > 0:
        print(stokes2)
        m2, b2 = np.polyfit(masses2,stokes2,deg=1)
        ax.scatter(masses2, stokes2)
        ax.plot(masses2, m2*masses2+b2,label="St = 0.06")
        ax.legend(loc="upper right")
    if len(masses3) > 0:
        print(stokes3)
        m3, b3 = np.polyfit(masses3,stokes3,deg=1)
        ax.scatter(masses3, stokes3)
        ax.plot(masses3, m3*masses3+b3,label="St = 0.02")
        ax.legend(loc="upper right")
    if len(masses4) > 0:
        print(stokes4)
        m4, b4 = np.polyfit(masses4,stokes4,deg=1)
        ax.scatter(masses4, stokes4)
        ax.plot(masses4, m4*masses4+b4,label="St = 0.006")
        ax.legend(loc="upper right")
    if len(masses5) > 0:
        print(stokes5)
        m5, b5 = np.polyfit(masses5,stokes5,deg=1)
        ax.scatter(masses5, stokes5)
        ax.plot(masses5, m5*masses5+b5,label="St = 0.002")
        ax.legend(loc="upper right")
    
    #ax.plot(xdat,testplot)
    #ax.plot(planet_masses, ring_widths)
    ax.set_xticks(xticks_list)
    #for i in range(len(peak_planet_distance)): 
        #ax.annotate(str(peak_planet_distance[i])[1:5], xy = (planet_masses[i],peak_planet_distance[i]), xytext = (planet_masses[i],peak_planet_distance[i]), textcoords = 'offset pixels', size = 20) 
        
     
#---------------------------------------
stokes = [0.2,0.063,0.02,0.0063,0.002]
stokesindex = 0 
#----------------------------------------

'''
# h = 0.02
orbnum=1500
alphavalue = 3
simsetname = 'a3-h2'
simnames = ['HLT-a3-1Me-24h-h2/','HLT-a3-2Me-24h-h2/','HLT-a3-3Me-24h-h2/','HLT-a3-5Me-24h-h2/','HLT-a3-7Me-24h-h2/','HLT-a3-10Me-24h-h2/','HLT-a3-13Me-24h-h2/']
xticks_list = [1,2,3,5,7,10,13]
maincode(simnames,alphavalue,ax1,xticks_list,stokesindex)
simsetname = 'a4-h2'
alphavalue = 4
simnames = ['HLT-a4-1Me-3day-h2/','HLT-a4-2Me-3day-h2/','HLT-a4-5Me-24h-h2/','HLT-a4-7Me-24h-h2/','HLT-a4-10Me-24h-h2/']
xticks_list = [1,2,3,5,7,10]
maincode(simnames,alphavalue,ax2,xticks_list,stokesindex)

#h = 0.05
alphavalue = 3
simsetname = 'a3-h5'
simnames = ['HLT-a3-1Me-3day-h5/','HLT-a3-2Me-3day-h5/','HLT-a3-3Me-3day-h5/','HLT-a3-5Me-3day-h5/','HLT-a3-7Me-3day-h5/','HLT-a3-10Me-3day-h5/','HLT-a3-13Me-3day-h5/','HLT-a3-17Me-3day-h5/','HLT-a3-21Me-3day-h5/','HLT-a3-26Me-3day-h5/','HLT-a3-31Me-3day-h5/']#,'HLT-a3-100Me-3day-h5/','HLT-a3-200Me-3day-h5/']
xticks_list = [1,2,3,5,7,10,13,17,21,26,31]
maincode(simnames,alphavalue,ax3,xticks_list,stokesindex)
alphavalue = 4
simsetname = 'a4-h5'
simnames = ['HLT-a4-1Me-3day-h5/','HLT-a4-2Me-3day-h5/','HLT-a4-3Me-3day-h5/','HLT-a4-5Me-3day-h5/','HLT-a4-7Me-3day-h5/','HLT-a4-10Me-3day-h5/','HLT-a4-13Me-3day-h5/']
xticks_list = [1,2,3,5,7,10,13]
maincode(simnames,alphavalue,ax4,xticks_list,stokesindex)

#h = 0.06
alphavalue = 3
simsetname = 'a3-h6'
simnames = ['HLT-a3-1Me-3day-h6/','HLT-a3-2Me-3day-h6/','HLT-a3-3Me-3day-h6/','HLT-a3-5Me-3day-h6/','HLT-a3-7Me-3day-h6-2/','HLT-a3-10Me-3day-h6-2/','HLT-a3-13Me-3day-h6-2/','HLT-a3-17Me-3day-h6-2/','HLT-a3-21Me-3day-h6-2/','HLT-a3-26Me-3day-h6/','HLT-a3-31Me-3day-h6/','HLT-a3-37Me-3day-h6/','HLT-a3-43Me-3day-h6/']#,'HLT-a3-200Me-3day-h6/']
xticks_list = [1,2,3,5,7,10,13,17,21,26,31,37,43]
maincode(simnames,alphavalue,ax5,xticks_list,stokesindex)
alphavalue = 4
simsetname = 'a4-h6'
simnames = ['HLT-a4-1Me-3day-h6/','HLT-a4-2Me-3day-h6/','HLT-a4-3Me-3day-h6/','HLT-a4-5Me-3day-h6/','HLT-a4-7Me-3day-h6/','HLT-a4-10Me-3day-h6/','HLT-a4-13Me-3day-h6/','HLT-a4-17Me-3day-h6/']
xticks_list = [1,2,3,5,7,10,13,17]
maincode(simnames,alphavalue,ax6,xticks_list,stokesindex)

#h = 0.08
alphavalue = 3
simsetname = 'a3-h8'
simnames = ['HLT-a3-1Me-3day-h8/','HLT-a3-2Me-3day-h8/','HLT-a3-3Me-3day-h8/','HLT-a3-5Me-3day-h8/','HLT-a3-7Me-3day-h8/','HLT-a3-10Me-3day-h8/','HLT-a3-13Me-3day-h8/','HLT-a3-17Me-3day-h8/','HLT-a3-21Me-3day-h8/','HLT-a3-26Me-3day-h8/','HLT-a3-31Me-3day-h8/','HLT-a3-37Me-3day-h8/','HLT-a3-43Me-3day-h8/']
xticks_list = [1,2,3,5,7,10,13,17,21,26,31,37,43]
maincode(simnames,alphavalue,ax7,xticks_list,stokesindex)
alphavalue = 4
simsetname = 'a4-h8'
simnames = ['HLT-a4-1Me-3day-h8/','HLT-a4-2Me-3day-h8/','HLT-a4-3Me-3day-h8/','HLT-a4-5Me-3day-h8/','HLT-a4-7Me-3day-h8/','HLT-a4-10Me-3day-h8/','HLT-a4-13Me-3day-h8/','HLT-a4-17Me-3day-h8/','HLT-a4-21Me-3day-h8/','HLT-a4-26Me-3day-h8/','HLT-a4-31Me-3day-h8/','HLT-a4-37Me-3day-h8/','HLT-a4-43Me-3day-h8/']#,'HLT-a4-100Me-3day-h8/']
xticks_list = [1,2,3,5,7,10,13,17,21,26,31,37,43]
maincode(simnames,alphavalue,ax8,xticks_list,stokesindex)
#h = 0.1
alphavalue = 3
simsetname = 'a3-h1'
simnames = ['HLT-a3-1Me-3day-h1/','HLT-a3-2Me-3day-h1/','HLT-a3-3Me-3day-h1/','HLT-a3-5Me-3day-h1/','HLT-a3-7Me-3day-h1/','HLT-a3-10Me-3day-h1/','HLT-a3-13Me-3day-h1/','HLT-a3-17Me-3day-h1/','HLT-a3-21Me-3day-h1/','HLT-a3-26Me-3day-h1/','HLT-a3-31Me-3day-h1/','HLT-a3-37Me-3day-h1/','HLT-a3-43Me-3day-h1/']#,'HLT-a3-200Me-7day-h1/']
xticks_list = [1,2,3,5,7,10,13,17,21,26,31,37,43]
maincode(simnames,alphavalue,ax9,xticks_list,stokesindex)
alphavalue = 4
simsetname = 'a4-h1'
simnames = ['HLT-a4-1Me-3day-h1/','HLT-a4-2Me-3day-h1/','HLT-a4-3Me-3day-h1/','HLT-a4-5Me-3day-h1/','HLT-a4-7Me-3day-h1/','HLT-a4-10Me-3day-h1/','HLT-a4-13Me-3day-h1/','HLT-a4-17Me-3day-h1/','HLT-a4-21Me-3day-h1/','HLT-a4-26Me-3day-h1/','HLT-a4-31Me-3day-h1/','HLT-a4-37Me-3day-h1/','HLT-a4-43Me-3day-h1/']
xticks_list = [1,2,3,5,7,10,13,17,21,26,31,37,43]
maincode(simnames,alphavalue,ax10,xticks_list,stokesindex)
'''

# JOSH PLOTS
# h = 0.02, BETA = 0.001
'''
orbnum=1500
alphavalue = 4
simsetname = 'h2-st5'
simnames = ['10Me-3day-h2-b3-st5_500','15Me-3day-h2-b3-st5_500','25Me-3day-h2-b3-st5_500']
xticks_list = [10,15,20,25,30]
maincode(simnames,alphavalue,ax1,xticks_list,stokesindex)
'''
# h = 0.05, BETA = 0.001

# orbnum=1500
# alphavalue = 4
# simsetname = 'h5-st5'
# simnames = ['10Me-h5-b3-st5_500','15Me-h5-b3-st5_500','25Me-h5-b3-st5_500','30Me-h5-b3-st5_500','40Me-h5-b3-st5_500','50Me-h5-b3-st5_500']
# #xticks_list = []
# xticks_list = [10,15,20,25,30,35,40,45,50,55,60]
# maincode(simnames,alphavalue,ax1,xticks_list,stokesindex)

# # h = 0.08, BETA = 0.001

# orbnum=1500
# alphavalue = 4
# simsetname = 'h8-st5_500'
# simnames = ['10Me-h8-b3-st5_500','15Me-h8-b3-st5_500','25Me-h8-b3-st5_500','30Me-h8-b3-st5_500','40Me-h8-b3-st5_500','50Me-h8-b3-st5_500']
# #xticks_list = [10,15,20,25,30,35,40,45,50,55,60]
# maincode(simnames,alphavalue,ax3,xticks_list,stokesindex)


# ISOTHERMAL SIMS
'''
orbnum=1500
alphavalue = 4
simsetname = 'h5-iso'
simnames = ['10Me-h5-iso','15Me-h5-iso','25Me-h5-iso','30Me-h5-iso','40Me-h5-iso','50Me-h5-iso']
xticks_list = [10,15,20,25,30,45,40,45,50,55,60]
maincode(simnames,alphavalue,ax3,xticks_list,stokesindex)
'''

orbnum=600
alphavalue = 4
# simsetname = 'h5-iso'
simnames = ['10Me', '20Me', '40Me']
xticks_list = [0,10,20,30,40,50,60]
maincode(simnames,alphavalue,ax3,xticks_list,stokesindex)

# -------------------------------------------------------------------------------------------------

filename = './finalfit.png'
plt.show()
#plt.savefig(filename, dpi=400)


sys.exit()


