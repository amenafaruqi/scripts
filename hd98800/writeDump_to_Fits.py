
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.interpolate import RegularGridInterpolator, interpn


model = 'hd98800_30m_transit'
files = np.arange(799,800)
wd = "."

skipped_files = []

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

for file in files:
    file_no = str(file).rjust(5, '0')
    print(f'Generating  density file number: {file_no}')
    try:
        density, limits, dims = readData(f'{wd}/{model}/newdump_{file_no}_density_g_cm3_grid.dat')    
        # limits = [0.000000E+00,4.908957E+01,-3.141593E+00,3.141593E+00,-1.768642E+01,2.112447E+01]
    except ValueError or FileNotFoundError:
        skipped_files.append(file)
        print('SKIPPING FILE: ',  file)
        continue

    x = np.linspace(limits[0], limits[1], dims[0])
    y = np.linspace(limits[2], limits[3], dims[1])
    z = np.linspace(limits[4], limits[5], dims[2])
    print(x.shape, y.shape, z.shape)
    print(density.shape)

    X, Y = np.meshgrid(x, y)
    X, Y = X.T, Y.T

#     fig, ax = plt.subplots()
#     ax.pcolormesh(X, Y, np.log10(np.sum(density, axis=2)))
#     plt.savefig('density_map1.png')
    fn = RegularGridInterpolator((x, y, z), density, bounds_error=False, fill_value=0) # Scale up here for bigger mcfost grid
    density_arr = np.zeros((1,n_az,n_z,n_rad))

    for ii in range(n_rad):
        for jj in range(n_z):
            rad = grid_data[0,0,jj,ii]
            z = grid_data[1,0,jj,ii]
            azimuth = grid_data[2,:,jj,ii]-np.pi
            density_arr[0, :, jj, ii] = fn((rad*np.cos(azimuth), rad*np.sin(azimuth), z))

    hdu = fits.PrimaryHDU(density_arr)
    hdu.writeto(f'./density_files_{model}/density_file_{file_no}.fits')

print('SKIPPED  FILES: ', skipped_files)

# fig, ax = plt.subplots()
# ax.imshow(np.log10(np.sum(density[0,:,:,:], axis=1)), origin='lower', cmap='gist_heat')
# plt.savefig('density_map2.png')
