import numpy as np
from scipy.interpolate import interp2d, interp1d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy import time
from astropy.io import fits
import pandas as pd
from scipy.interpolate import RegularGridInterpolator, interp1d
from datetime import datetime, timedelta


# home_dir = "/storage/astro2/phrkvg"
home_dir = '/home/astro/phrkvg'
wd = "hd98800_30m_transit"
dm = 0.33
wl = 0.5
odms_dir = f'./odms_{wl}um/odms_{wd}_{str(dm)}'
files = np.arange(0,1901,1)
dt = 0.05416
start_file = 3600
start_dt = 0.05416

# ========================================================

# Binary parameters (Zuniga-Fernandez et al. 2021)
t_ref = time.Time(2023, format='decimalyear').mjd   # t_ref = T0_AB 

e_B = 0.805
i_B = -np.radians(66.3)
a_B = 1.01
m_Ba = 0.77
m_Bb = 0.62
omega_Ba = np.radians(104.5)
omega_Bb = omega_Ba - np.pi
Omega_B = np.radians(342.7)+np.pi/2
T0_B = 48707.5
P_B = 314.86
mean_B = (2*np.pi*((t_ref - T0_B)/P_B))%(2*np.pi)

e_A = 0.4808
i_A = -np.radians(135.6)
a_A = 0.86
m_Aa = 0.93
m_Ab = 0.29
omega_Aa = np.radians(68.7)
omega_Ab = omega_Aa - np.pi
Omega_A = np.radians(170.2)+np.pi/2
T0_A = 48742.5
P_A = 264.51
mean_A = (2*np.pi*((t_ref - T0_A)/P_A))%(2*np.pi)

e_AB = 0.46
i_AB = -np.radians(88.1)
a_AB = 51
m_A = 1.22
m_B = 1.4
omega_A = np.radians(65)
omega_B = omega_A - np.pi
Omega_AB = np.radians(184.5)+np.pi/2
T0_AB = 2023
P_AB = 230
mean_AB = 0.0


# Disc parameters (Kennedy et al. 2019)
a_inner = 2.5
a_outer = 4.6
e_disc = 0.03
pos_disc = np.radians(15.6)
inc_disc = np.radians(26)
omega_disc = np.radians(-73)



def merge_ev_files(sink, wd='./hd98800', ev_files=1):
    data = pd.DataFrame()
    for file_number in range(1, ev_files+1):
        data_subset = pd.read_csv(
        f'{wd}/discSink000{sink}N0{file_number}.ev', sep='\s+', skiprows=[0], engine='python', names=[str(i) for i in range(20)]
    )
        data = pd.concat([data, data_subset], ignore_index=True)
    
    return data

data_A = np.array([[0.429, 0.525, 0.954, 0.955, 1.082, 1.251, 1452, 1.658, 1.875],
                   [266, 700, 2080, 2000, 2440, 2000, 2470, 2300, 2190]])
data_B = np.array([[0.429, 0.525, 0.954, 0.955, 1.082, 1.251, 1.658, 2.20],
                   [125, 430,1880, 1890, 2300, 2200, 2500, 1790]])
def interpolate_flux(wavelength, roundto=-1):
    interpA = interp1d(data_A[0], data_A[1], bounds_error=False, fill_value=0)
    interpB = interp1d(data_B[0], data_B[1], bounds_error=False, fill_value=0)
    f_A = np.round(interpA(wavelength),roundto)
    f_B = np.round(interpB(wavelength),roundto)
    
    if wavelength == 0.8:
        fr = 3.9   # flux ratio Aa/Ab at 0.8um
    else:
        fr = 12.3

    f0_Aa = (f_A*fr)/(1+fr)
    f0_Ab = f_A - f0_Aa
    return f_A, f_B, f0_Aa, f0_Ab

def month_year_date(decyear):
    year = int(decyear)
    rem = decyear - year
    base = datetime(year, 1, 1)
    result = base + timedelta(seconds=(base.replace(year=base.year + 1) - base).total_seconds() * rem)
    datemy = datetime.strftime(result, "%m-%Y")
    return datemy


# Note only using discSink files for part 2 of dataset
data1 = merge_ev_files(1, wd=wd, ev_files=1)    # Ab
data3 = merge_ev_files(3, wd=wd, ev_files=1)    # Aa
data2 = merge_ev_files(2, wd=wd, ev_files=1)    # Bb
data4 = merge_ev_files(4, wd=wd, ev_files=1)    # Ba

# times = ((start_file*start_dt)+((files-start_file)*dt))/(2*np.pi)
times = ((files*dt) + (start_file*start_dt))/(2*np.pi)
tau_Aa = np.zeros(len(files))
tau_Ab = np.zeros(len(files))
# No need to store tau at Ba and Bb

Aa_coords = np.empty((len(files), 2))
Ab_coords = np.empty((len(files), 2))
Ba_coords = np.empty((len(files), 2))
Bb_coords = np.empty((len(files), 2))

# Get x,y positions of stars Aa and Ab over time
Aa_pos = np.stack((np.array(data3['1']),np.array(data3['2'])))
Ab_pos = np.stack((np.array(data1['1']),np.array(data1['2'])))
Ba_pos = np.stack((np.array(data4['1']),np.array(data4['2'])))
Bb_pos = np.stack((np.array(data2['1']),np.array(data2['2'])))

# Shift all positions relative to the CoM of B
com_B = (m_Bb*Bb_pos + m_Ba*Ba_pos)/m_B
Aa_pos = Aa_pos - com_B
Ab_pos = Ab_pos - com_B
Ba_pos = Ba_pos - com_B
Bb_pos = Bb_pos - com_B

# Create interpolators for x,y positions of stars at time t
interp_Ba = interp1d(np.array(data4['0']),Ba_pos, bounds_error=False, fill_value=0)
interp_Bb = interp1d(np.array(data2['0']),Bb_pos, bounds_error=False, fill_value=0)
interp_Aa = interp1d(np.array(data3['0']),Aa_pos, bounds_error=False, fill_value=0)
interp_Ab = interp1d(np.array(data1['0']),Ab_pos, bounds_error=False, fill_value=0)

# Get r, phi values for ODM grid
grid_file = 'data_disk/grid.fits.gz'
hdu_list = fits.open(grid_file)
grid_data = hdu_list[0].data

n_phi = len(grid_data[0,:,0,0])
n_rad = len(grid_data[0,0,0,:])
phi = np.linspace(0, np.max(grid_data[2,:,:,:]), n_phi)
r = np.logspace(np.log10(np.min(grid_data[0,:,:,:])), np.log10(np.max(grid_data[0,:,:,:])), n_rad)


for i, file in enumerate(files):
    file_no = str(file).rjust(5, '0')
    code_time = times[i]*2*np.pi
    
    # Get x,y coords at time t using interpolation 
    x_Aa, y_Aa = interp_Aa(code_time)
    x_Ab, y_Ab = interp_Ab(code_time)
    x_Ba, y_Ba = interp_Ba(code_time)
    x_Bb, y_Bb = interp_Bb(code_time)
    
    Aa_coords[i] = (x_Aa, y_Aa)
    Ab_coords[i] = (x_Ab, y_Ab)
    Ba_coords[i] = (x_Ba, y_Ba)
    Bb_coords[i] = (x_Bb, y_Bb)

    
    # Convert x,y to r,phi
    r_Aa = (x_Aa**2 + y_Aa**2)**0.5
    r_Ab = (x_Ab**2 + y_Ab**2)**0.5
    phi_Aa = np.arctan2(-y_Aa,-x_Aa)
    phi_Ab = np.arctan2(-y_Ab,-x_Ab)
    if phi_Aa < 0: phi_Aa = phi_Aa + 2*np.pi 
    if phi_Ab < 0: phi_Ab = phi_Ab + 2*np.pi 
    
    r_Ba = (x_Ba**2 + y_Ba**2)**0.5
    r_Bb = (x_Bb**2 + y_Bb**2)**0.5
    phi_Ba = np.arctan2(-y_Ba,-x_Ba)
    phi_Bb = np.arctan2(-y_Bb,-x_Bb)
    if phi_Ba < 0: phi_Ba = phi_Ba + 2*np.pi 
    if phi_Bb < 0: phi_Bb = phi_Bb + 2*np.pi 
    
    # Create interpolator for r,phi,tau from ODM
    try:
        opt_depth = fits.open(f'{odms_dir}/optical_depth_map_{file_no}.fits.gz')    # optical depth file
        opt_depth_data = opt_depth[0].data    # data in format [direction, phi, z, r]
        tau = opt_depth_data[1,:,0,:]

        fn = RegularGridInterpolator((phi,r), tau, bounds_error=False,  fill_value=0)

        # Get tau at r,phi of Aa and Ab using interpolation
        tau_Aa[i] = fn((phi_Aa, r_Aa))
        tau_Ab[i] = fn((phi_Ab, r_Ab))
        print(file, tau_Aa[i], tau_Ab[i])

    except:
        # FileNotFoundError or ValueError or OSError
        tau_Aa[i] = np.nan
        tau_Ab[i] = np.nan
        print(f'Error, corrupt or missing FITS file: {file}')
        


full_data = np.stack((times, tau_Aa, tau_Ab, Aa_coords[:,0], Aa_coords[:,1], Ab_coords[:,0], Ab_coords[:,1], Ba_coords[:,0], Ba_coords[:,1], Bb_coords[:,0], Bb_coords[:,1]), axis=1)
np.savetxt(f'tau_{wd}_{dm}_{wl}um.txt',full_data) # Data in format [time, tau_Aa, tau_Ab, x_Aa, y_Aa, x_Ab, y_Ab]

