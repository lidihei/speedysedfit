# Speedysedfit 

Speedysedfit is modified from [speedyfit](https://github.com/vosjo/speedyfit) of Joris Vos 

A python package to fit the photometric spectral energy distribution of stars. Uses a Markov chain Monte Carlo approach 
to determine the errors on the derived parameters.

Speedyfit is a command line tool writen in Python 3 that allows you to search the most common online databases for 
photometric observations of your target, and fit theoretical atmosphere models to the obtained photometry. Speedyfit can
deal with both single and binary stars, and allows for the inclusion of constraints from other sources, as for example
the distance or reddening. 

Installing speedysedfit, please look at ./ipynb/install_speedysedfit.ipynb