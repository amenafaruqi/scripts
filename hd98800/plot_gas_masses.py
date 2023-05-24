import numpy as np
from scipy.interpolate import interp2d, interp1d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy import time
from astropy.io import fits
import pandas as pd
from scipy.interpolate import RegularGridInterpolator, interp1d
from datetime import datetime, timedelta

import seaborn as sns

sns.set(rc={"figure.dpi":150, 'savefig.dpi':150})

sns.set(font='Arial',
        rc={
 'axes.axisbelow': False,
 'axes.edgecolor': 'black',
 'axes.facecolor': 'None',
 'axes.grid': False,
 'axes.labelcolor': 'black',
 'axes.spines.right': True,
 'axes.spines.top': True,
 'figure.facecolor': 'white',
 'lines.solid_capstyle': 'round',
 'patch.edgecolor': 'w',
 'patch.force_edgecolor': True,
 'text.color': 'black',
 'xtick.bottom': True,
 'xtick.color': 'black',
 'xtick.direction': 'out',
 'xtick.top': False,
 'ytick.color': 'black',
 'ytick.direction': 'out',
 'ytick.left': True,
 'ytick.right': False,
 'figure.dpi':100,
 'savefig.dpi':50,
        })
sns.set_context("notebook", rc={"font.size":13,
                                "axes.titlesize":13,
                                "axes.labelsize":13})
sns.set_palette("colorblind", color_codes=True)


t_ref = time.Time(2023, format='decimalyear').mjd   # t_ref = T0_AB
P_B = 314.86


def fill_disc(inner_r, outer_r, ang, axis, colour='k', alpha=0.5):
    a, b = [inner_r, outer_r], [inner_r*np.cos(inc_disc)*(1-e_disc**2)**0.5,
                                outer_r*np.cos(inc_disc)*(1-e_disc**2)**0.5]

    x = np.outer(a, np.cos(ang))
    y = np.outer(b, np.sin(ang))
    xs = x*np.cos(pos_disc) - y*np.sin(pos_disc)
    ys = x*np.sin(pos_disc) + y*np.cos(pos_disc)
    xs[1,:] = xs[1,::-1]
    ys[1,:] = ys[1,::-1]
    axis.fill(np.ravel(xs), np.ravel(ys), color=colour, edgecolor=None, alpha=alpha)

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



fig, ax = plt.subplots(1, figsize=(10,7))

f_A, f_B, f0_Aa, f0_Ab = interpolate_flux(0.8)

data_alpha_3= np.loadtxt(f'tau_hd98800_3m_transit_0.33_0.8um.txt')
data_alpha_33 = np.loadtxt(f'tau_hd98800_30m_transit_0.33_0.8um.txt')
data_alpha_330 = np.loadtxt(f'tau_hd98800_300m_transit_0.33_0.8um.txt')

times = data_alpha_3[:,0]
tau_Aa_3 = data_alpha_3[:,1]
tau_Ab_3 = data_alpha_3[:,2]
tau_Aa_33 = data_alpha_33[:,1]
tau_Ab_33 = data_alpha_33[:,2]
tau_Aa_330 = data_alpha_330[:,1]
tau_Ab_330 = data_alpha_330[:,2]


f_tot_3 = f0_Aa*np.exp(-tau_Aa_3) + f0_Ab*np.exp(-tau_Ab_3) + f_B
f_tot_33 = f0_Aa*np.exp(-tau_Aa_33) + f0_Ab*np.exp(-tau_Ab_33) + f_B
f_tot_330 = f0_Aa*np.exp(-tau_Aa_330) + f0_Ab*np.exp(-tau_Ab_330) + f_B
t_year = time.Time(t_ref - (40*P_B) + (times)*365.25, format='mjd').decimalyear

ax.plot(t_year, f_tot_3, label='$M_{gas}$=3.3 $M_{{\oplus}}$')
ax.plot(t_year, f_tot_33, label='$M_{gas}$=33 $M_{{\oplus}}$')
ax.plot(t_year, f_tot_330, label='$M_{gas}$=330 $M_{{\oplus}}$')
ax.set_ylabel('Flux (mJy)')
ax.set_xlabel('Year')
ax.set_xlim(np.ceil(min(t_year))+1, (max(t_year)))
# print(np.ceil(min(t_year)))

# f_norm = lambda f: f/(f_A+f_B)
# f_denorm = lambda f: f*(f_A+f_B)


# secax0 = ax.secondary_yaxis(location='right', functions=(f_norm, f_denorm))
# secax0.set_ylabel('Nomalised flux')
plt.legend()

plt.tight_layout()
plt.savefig('gas_masses.png')
