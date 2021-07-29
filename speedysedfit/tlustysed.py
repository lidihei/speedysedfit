import numpy as np
import re

re_dbl_fort = re.compile(r'(\d*\.\d+)[dD]([-+]?\d+)')
re_dbl_fort1 = re.compile(r'(\d*\.\d+)[-](\d+)')
re_dbl_fort2 = re.compile(r'(\d*\.\d+)[+](\d+)')

def rsed(fname):
    '''read sed file of Tlusty dowloaded from http://tlusty.oca.eu/Tlusty2002/OS02-SED.html
    parameters:
    --------------------
    fname: [str] e.g. BG28000g475v2.flux
    returns:
    -------------------
    wavelength: [array] in units of angstrom
    flux: [array] in units of erg/cm^2/s/A
    '''
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    nlines = len(lines)
    freq = np.zeros(nlines)
    flux_freq = np.zeros(nlines)
    for _i, line0 in enumerate(lines):
       try:
          line0 =re_dbl_fort.sub(r'\1e\2', line0)
          line = line0.split()
          freq[_i] = float(line[0])
          flux_freq[_i] = float(line[1])
       except:
          line0 = re_dbl_fort1.sub(r'\1e-\2', line0)
          line0 = re_dbl_fort2.sub(r'\1e+\2', line0)
          line = line0.split()
          freq[_i] = float(line[0])
          flux_freq[_i] = float(line[1])
    c = 2.99792458e18 # speed of light A/sec
    lamA = c/freq
    flux  = 4.*np.pi * c/lamA**2 * flux_freq
    ind = np.argsort(lamA)
    wavelength = lamA[ind]
    flux = flux[ind]
    return wavelength, flux
