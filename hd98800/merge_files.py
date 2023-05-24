from subprocess import call
import numpy as np


files = np.arange(150)
sim0 =  'hd98800_3m_transit0'
sim = 'hd98800_3m_transit'

for i, file in enumerate(files):
    file_no = str(file).rjust(4, '0')
    new_no = str(file+41).rjust(4, '0')
    call(['cp', f'{sim}/splash_{file_no}.png', f'hd98800_3m_imgs/splash_{new_no}.png'])

# for i, file in enumerate(files):
#     file_no = str(file).rjust(5, '0')
#     new_no = str(file-3600).rjust(5, '0')
#     call(['cp', f'density_files_{sim0}/density_file_{file_no}.fits', f'density_files_{sim}/density_file_{new_no}.fits'])

