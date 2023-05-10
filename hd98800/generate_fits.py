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
dumpfiles = np.arange(350,551,1)     # For efficiency, only run after transients have dissipated
density_files = np.arange(400,401)
wl = 0.8             # [microns]
wd = f'./{sim}/'
old_dirs = ['data_th_old', f'data_{wavelength}_old']
density_files_dir = f'density_files_{sim}'
# odms_dir = f'odms_{sim}'
phantom_dt = 0.05416
dt = phantom_dt/(2*np.pi) # [years]
file_prefix = 'newdump_'
gas_mass =  33
dust_masses = [0.33]


def delete_old_data(old_dirs):
    print(f'---------------- DELETING FITS DATA FROM PREVIOUS RUN -----------------')
    for old_dir in old_dirs:
        call(['rm', '-rf', old_dir])


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
    # elif os.listdir(f'./{odms_dir}'):                   # pick up where last left off
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

    make_ODMs(density_files, odms_dir, wavelength=wl)

