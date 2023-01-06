# Speedysedfit 

Speedysedfit is modified from [speedyfit](https://github.com/vosjo/speedyfit) of Joris Vos 

A python package to fit the photometric spectral energy distribution of stars. Uses a Markov chain Monte Carlo approach 
to determine the errors on the derived parameters.

Speedyfit is a command line tool writen in Python 3 that allows you to search the most common online databases for 
photometric observations of your target, and fit theoretical atmosphere models to the obtained photometry. Speedyfit can
deal with both single and binary stars, and allows for the inclusion of constraints from other sources, as for example
the distance or reddening. 

Installing speedysedfit, please look at ./ipynb/install_speedysedfit.ipynb

- $git clone https://github.com/lidihei/speedysedfit.git
- $cd speedysedfit
- $python setup.py intall

#setup configure file download [SED grids](http://www.astro.physik.uni-potsdam.de/~jorisvos/Speedyfit/modelgrids.tar.gz)

- make a directory which you want to store the SED grids and the integrated grids

- $ mkdir /share/speedysedfit  
- $ mkdir /share/speedysedfit/modelgrids
- $ export SPEEDYFIT_MODELS="share/speedyfit/modelgrids/"


- $ cp -r speedysedfit/transmission_curves $SPEEDYFIT_MODELS/../
- $ cp -r speedysedfit/redlaws $SPEEDYFIT_MODELS/../
- $ cp speedysedfit/zeropoints.dat $SPEEDYFIT_MODELS/../

# edit a default grid configure file
- $ cd $SPEEDYFIT_MODELS
- create or edit a file named grid_description.yaml which contains(TLUSTY and COMARCS grids should be created by youself):


kurucz:
    filename: 'kurucz93_z0.0_k2odfnew_sed'
munari:
    filename: 'Munari2005_extended'
tmap:
    filename: 'TMAP2012_sdOB_extended'
blackbody:
    filename: 'blackbody_discint'
tlustyB:
    filename: 'TlustyB_G_v2'
tlustyO:
    filename: 'TlustyO_G_v10'
comarcs:
    filename: 'COMARCS'


# Create [TLUSTY SED grid](creat_tlusty_sed_gride.ipynb) and [COMARCS SED grid](creat_comarcs_sed_gride.ipynb)

If you want calculate and save Integrated TLUSTY and COMARCS grids, or calculate a your own integrated grid, you can look at [calc_integrated_grid.ipynb](calc_integrated_grid.ipynb)
