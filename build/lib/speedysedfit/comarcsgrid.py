import numpy as np
import os

def read_comarcs_spec(fname):
    '''read spec of comarcs
    returns:
    ------------
    wave [array] - wavelength [Å]; 
    flux_norm: [array] F/F(Cont) - continuum normalized flux (normalized to calculation without atomic and molecular line opacities); 
    L: [array] νL(ν) - frequency times specific luminosity [erg/s]; 
    fmean [array] flux devided by its mean value
    flam [array] flux [erg/s/cm2/AA]
    '''
    wave, flux_norm, nuLnu, fmean = np.loadtxt(fname, unpack=True)
    basename = os.path.basename(fname)
    paras = basename.split('_')
    teff = np.float(paras[1][1:])
    logg = np.float(paras[2][1:])/100
    m = np.float(paras[3][1:])/100

    ratio = 5.996233e-28 # 1/(4pi R^2) = (1 cm/s2) /G * 1*Msun/4*pi =5.996233e−28
    flam = nuLnu*10**logg/m *ratio/wave

    return wave, flux_norm, nuLnu, fmean, flam
