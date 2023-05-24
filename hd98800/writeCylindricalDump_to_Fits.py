import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.interpolate import RegularGridInterpolator, interpn
from tqdm import tqdm

def readData(filename):
    with open(filename) as f:
        y_idx = 0
        z_idx = 0
        for i, line in enumerate(f):
            if i == 21:
                dims = line.split('#')[-1].strip()
                dims = np.array(dims.split(), dtype=np.int)
                density = np.zeros(dims)
                # exit()
            if i == 7:
                limits = line.split('#')[1].strip()
                limits = np.array(limits.split(), dtype=np.float)
            if ('#' not in line) and (i != 21):
                linedata = np.array(line.split(), dtype=np.float)
                density[:, y_idx, z_idx] = linedata
                y_idx += 1
                if y_idx == dims[1]:
                    y_idx = 0
                    z_idx += 1
                # print(y_idx,z_idx)
    return density, limits, dims    


grid_file = 'data_disk/grid.fits.gz'
hdu_list = fits.open(grid_file)
grid_data = hdu_list[0].data
print('Grid shape: ', grid_data.shape)

n_az = len(grid_data[0,:,0,0])
n_z = len(grid_data[0,0,:,0])
n_rad = len(grid_data[0,0,0,:])

files = np.arange(401,551)

for file in files:
    print(f'Generating  density file number: {file}')
    density, limits, dims = readData(f'./hd98800/newdump_00{file}_density_g_cm3_grid.dat')

    r = np.linspace(limits[0], limits[1], dims[0])
    phi = np.linspace(limits[2], limits[3], dims[1])
    z = np.linspace(limits[4], limits[5], dims[2])
    print(r.shape, phi.shape, z.shape)
    print(density.shape)

    R, Phi = np.meshgrid(r, phi)
    R, Phi = R.T, Phi.T

#     fig, ax = plt.subplots()
#     ax.pcolormesh(R*np.cos(Phi), R*np.sin(Phi), np.log10(np.sum(density, axis=2)))
#     plt.savefig('density_map1.png')
    fn = RegularGridInterpolator((r, phi, z), density, bounds_error=False, fill_value=0)

    density_arr = np.zeros((1,n_az,n_z,n_rad))
    pbar = tqdm(total=n_az*n_z*n_rad)

    for ii in range(n_rad):
        for jj in range(n_z):
            for kk in range(n_az):
                rad = grid_data[0,kk,jj,ii]
                z = grid_data[1,kk,jj,ii]
                azimuth = grid_data[2,kk,jj,ii]-np.pi
                density_arr[0, kk, jj, ii] = fn((rad, azimuth, z))
                pbar.update()

    pbar.close()

    hdu = fits.PrimaryHDU(density_arr)
    hdu.writeto(f'density_file_00{file}.fits')


#     fig, ax = plt.subplots()
#     ax.imshow(np.log10(np.sum(density[0,:,:,:], axis=1)), origin='lower', cmap='gist_heat')
#     plt.savefig('density_map2.png')

