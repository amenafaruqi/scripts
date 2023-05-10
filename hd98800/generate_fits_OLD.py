from matplotlib import use
import numpy as np
from astropy.io import fits 
import matplotlib.pyplot as plt
from subprocess import call, check_output
import os
import re

# Constants - edit these before running
sim = 'hd98800_30m_transit'
parafile = 'hd98800.para'     # make sure this is inside wd
dumpfiles = np.arange(350,551,1)     # For efficiency, only run after disc has stabilised
density_files = np.arange(400,401)
wavelength = 0.5   # [microns]
wd = f'./{sim}/'
old_dirs = ['data_th_old', f'data_{wavelength}_old']
fits_dir = 'fits_no_scatt'
density_files_dir = f'density_files_{sim}'
# odms_dir = f'odms_{sim}'
phantom_dt = 0.05416
dt = phantom_dt/(2*np.pi) # [years]
file_prefix = 'newdump_'
limits_file = 'limits.txt'
model_limits = True
gas_mass =  33
dust_masses = [0.33]



def delete_old_data(old_dirs):
    print(f'---------------- DELETING FITS DATA FROM PREVIOUS RUN -----------------')
    for old_dir in old_dirs:
        call(['rm', '-rf', old_dir])

def store_optical_depths(output):
    taus = []
    for star in range(1,5):
        tau = output.decode('unicode_escape').split(f'Optical depth from star # {star} is ')[1].split('\n')[0]
        taus.append(tau)
        print(tau)
    return taus
 
def generate_fits_data(dumpfiles, wavelength, taus_file=None):
    if not dumpfiles[0]:
        delete_old_data([f'./{fits_dir}/*'])

    for i, dumpfile in enumerate(dumpfiles):
        if not i:
            delete_old_data(old_dirs)      # delete previous run FITS files
            file_no = str(dumpfile).rjust(5, '0')
            filename = 'newdump_' + file_no
            output = check_output(['mcfost', wd + parafile, '-phantom', wd + filename])
            update_limits(output)

        delete_old_data(old_dirs)      # delete previous run FITS files
        file_no = str(dumpfile).rjust(5, '0')
        filename = 'newdump_' + file_no
        call(['mcfost', wd + parafile, '-limits', 'limits.txt', '-phantom', wd + filename, '-no_scattering'])
        print(f'------- TEMPERATURE CALCULATION COMPLETE FOR FILE: {filename} ---------')
        output = check_output(['mcfost', wd + parafile, '-limits', 'limits.txt', '-phantom', wd + filename, '-img', str(wavelength), '-no_scattering'])
        print(f'------- SCATTERED LIGHT IMAGE GENERATED FOR FILE:  {filename} ---------')
        update_limits(output)
        if taus_file:
            taus = store_optical_depths(output)
            taus_str = ' '.join([str(t) for t in taus]) + '\n'  
            file = open(taus_file, 'a')
            file.write(taus_str)
            file.close() 

        if i > 0:
            prev_file_no = str(dumpfiles[i-1]).rjust(5, '0')
            print(f'-------------------- MOVING FITS FILE NUMBER {prev_file_no} --------------------')
            call(['cp', f'./data_{wavelength}_old/RT.fits.gz', f'./{fits_dir}/RT_{prev_file_no}.fits.gz'])
    
    print(f'-------------------- MOVING FITS FILE NUMBER {file_no} --------------------')
    call(['cp', f'./data_{wavelength}/RT.fits.gz', f'./{fits_dir}/RT_{file_no}.fits.gz'])


def read_existing_fits(files, directory='./fits', log=True):
    full_data = np.zeros((201,201,len(files)))     # array to store all data in

    for i, file in enumerate(files):
        file_no = str(file).rjust(5, '0')
        data = fits.open(f'{directory}/RT_{file_no}.fits.gz')[0].data[0, :, 0, :,][0,:,]
        if log:
            data = np.log(data, out=np.zeros_like(data), where=(data!=0))
        full_data[:,:,i] = data
        
    return full_data

def update_limits(mcfost_output, use_model_limits=False):
    if use_model_limits:
        start_str = '# Model limits :\n'
        end_str = '\n Using'
    else:
        start_str = '# Farthest particules :\n'
        end_str = '\n Found'

    str_output = mcfost_output.decode('unicode_escape').split(start_str)[1].split(end_str)[0]
    lower_lims = np.array([float(s) for s in str_output.split()[2::4]]) - 3
    upper_lims = np.array([float(s) for s in str_output.split()[3::4]]) + 3
    limits = np.concatenate((lower_lims, upper_lims))
    str_limits = ''.join([str(lim) + '\n' for lim in limits])
    with open(limits_file, 'w') as f:
        f.write(str_limits)
    print(f'UPDATING LIMITS TO {limits}')


def generate_lightcurve(data):
    n_pts = data.shape[2]
    flux = np.zeros(n_pts)
    time = np.arange(0, dt*n_pts, dt)
    for i in range(n_pts):
        pixel_data = data[:,:,i]
        flux[i] = np.sum(pixel_data)

    flux /= np.max(abs(flux))
    
    plt.figure()
    plt.scatter(time, flux)
    plt.xlabel('Time (days)')
    plt.ylabel('Normalised Flux')
    plt.show()
        

def make_gif(data, filename='transit.gif'):
    n_imgs = data.shape[2]
    images = []
    for n in range(n_imgs):
        img_data = data[:,:,n]
        img_name = f'img_{str(n*10).rjust(5, "0")}.png'
        plt.imsave(img_name, img_data)
        images.append(imageio.imread(img_name))
    
    imageio.mimsave(filename, images)
    call('rm *.png', shell=True)


def mod_phantom_dumps(dumpfiles, modfile='moddump_CoM.f90'):
    call(['make', 'moddump', 'MODFILE=' + modfile])
    for i, dumpfile in enumerate(dumpfiles):
        file_no = str(dumpfile).rjust(5, '0')
        call([f'./{wd}phantommoddump', f'./{wd}{file_prefix}{file_no}', f'./{wd}newdump_{file_no}'])


def make_ODMs(files, odms_dir, wavelength=0.8):
    for i, file in enumerate(files):
        if i > 0:
            prev_file_no = str(files[i-1]).rjust(5, '0')
            print(f'-------------------- MOVING FITS FILE NUMBER {prev_file_no} --------------------')
            call(['cp', 'optical_depth_map.fits.gz', f'./{odms_dir}/optical_depth_map_{prev_file_no}.fits.gz'])

        print(f'============== Progress: {i}/{len(files)} ==================') 
        file_no = str(files[i]).rjust(5, '0')
        filename = f'{density_files_dir}/density_file_{file_no}.fits'
        delete_old_data(old_dirs)      # delete previous run FITS files
        call(['mcfost', parafile,  '-density_file', filename,  '-3D', '-od'])
        call(['mcfost', parafile,  '-density_file', filename,  '-3D', '-img', str(wavelength), '-od'])

    print(f'-------------------- MOVING FITS FILE NUMBER 00{file} --------------------')
    call(['cp', 'optical_depth_map.fits.gz', f'./{odms_dir}/optical_depth_map_{file_no}.fits.gz'])


delete_old_data(['data_th_old'])
gm = gas_mass*3.0027e-6

for dust_mass in dust_masses:
    dm = dust_mass*3.0027e-6      # convert units from  Mearth to Msun
    gdr = round(gas_mass/dust_mass, 3)
    odms_dir = f'odms_{str(wavelength)}um/odms_{sim}_{str(dust_mass)}'
    if not os.path.exists(f'./{odms_dir}'):
        call(['mkdir', f'./{odms_dir}'])
    # elif os.listdir(f'./{odms_dir}'):
    #     existing_odms = os.listdir(f'./{odms_dir}')
    #     existing_odms.sort()
        # last_file_no = int(existing_odms[-1].split('.fits')[0][-5:])
        # if last_file_no == 5500:
        #     continue
        # density_files = np.arange(last_file_no+1,5501)
        # print('Skipping to file number: ' + str(last_file_no+1))
        
    print('Gas-dust ratio: ', gdr)
    print('Gas mass: ', gas_mass)
    print('Dust mass: ', dust_mass)

    contents = open(parafile).read()
    edited_contents = re.sub('wall\n.*?\t\t  dust mass',f'wall\n  {dm}    {gdr}\t\t  dust mass', contents, flags=re.DOTALL)
    with open(parafile, 'w') as f:
        f.write(edited_contents)

    make_ODMs(density_files, odms_dir, wavelength=wavelength)

# generate_fits_data(dumpfiles, wavelength)
