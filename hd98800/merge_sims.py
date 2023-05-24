from subprocess import call
import numpy as np


files = np.arange(1203,1399)
sim0 =  'longsim0'
sim = 'longsim'
dm = 0.33

for i, file in enumerate(files):
    file_no = str(file).rjust(5, '0')
    new_no = str(file-397).rjust(5, '0')
    call(['cp', f'odms_0.8um/odms_{sim0}_{dm}/optical_depth_map_{file_no}.fits.gz', f'odms_0.8um/odms_{sim}_{dm}/optical_depth_map_{new_no}.fits.gz'])

# for i, file in enumerate(files):
#     file_no = str(file).rjust(5, '0')
#     new_no = str(file-3600).rjust(5, '0')
#     call(['cp', f'density_files_{sim0}/density_file_{file_no}.fits', f'density_files_{sim}/density_file_{new_no}.fits'])

