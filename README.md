## FARGO Analysis Scripts
#### `2Dsurface_density.py`

Plots the 2D surface density of gas and/or dust in g/cm^2.

Run using 
```
python 2Dsurface_density.py -o [output] -wd [simulations working directory]  -dir [simulation subdirectory] -d [dust size] -cb [colorbar min value, max value] -lin [y/n] -gas [y/n]  -planet [y/n] -zoom [limit]
```

- output: FARGO output number (required)
- wd: directory where simulations subdirectories are stored (required)
- dir: name  of  simulation  subdirectory (required)
- dust size: dust size bin number, default value is -1, which plots the mass averaged sum of all dust grain sizes. 
- colorbar values: minimum and maximum limits on colourbar. By default, the min and  max values are set according to the min and  max density values.
- linear scale: yes or no depending on if  a  log or linear density  scale is desired. Default is no.
- gas: yes  or no depending  on  if a plot of  gas surface density is  desired. Default is no.
- planet: yes or no depending on if plotting the planet's location is  desired. Default is no.
- zoom: distance in AU to zoom  to e.g. 50 will zoom to the inner 100 AU (-50  to +50) of the grid in  the x and y directions. By default, there no  zoom is applied. 


#### `plot_sigma.py`

Plots the 1D surface density of gas and/or dust in g/cm^2 against different physical quantities.

Run using  
```
python plot_sigma.py -o [outputs] -wd [simulations working directory] -sim [simulation subdirectory]  -savedir [directory to save plots to] -noplanet  -nogrog -plot_window -porbits -style [style]
```

- outputs: FARGO output numbers, can  input more than one (required)
- wd: directory where simulations subdirectories are stored (required)
- sim: name  of  simulation  subdirectory (required)
- savedir: directory to save plots to. Default is `./images` (assuming it exists).
- noplanet: use this flag if there is no planet in the simulation.
- nogrog:  use this  flag if plotting a `dusty_fargo` simulation (not `grain_growth`)
- plot_window: use this  flag to open an interactive  plot window rather than saving  the plots as images.
- porbits:  use this flag to express timestepsin planet orbits rather than Myr.
- style: choen  matplotlib style, default is "publication" style.

