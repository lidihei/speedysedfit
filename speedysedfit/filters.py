import glob
import os
from astropy import units
import numpy as np

from astropy.io import ascii, fits
#basedir = os.path.dirname(__file__)
basedir = os.path.join(os.environ.get('SPEEDYFIT_MODELS', None), '..')  # modified by lijiao


def Portalsednv2lambda(Portalsed, filterpara='filter_parameters.fit'):
    '''
    Transform the frequence (nv) and flux_nv of the SED download from http://cdsportal.u-strasbg.fr/ to wave and flux_lambda
    parameters:
    ---------------
    Portalsed: [astropy Tab] the SED downloaded from http://cdsportal.u-strasbg.fr/ 
    filterpara: [str] the *.fit file which can be downloaded from http://vizier.u-strasbg.fr/viz-bin/VizieR-4 
    returns:
    ---------------
    sed_lambda [array] in units of Angstrom
    sed_flux_lambda [array] flux in units of erg/s/cm/cm/AA'
    sed_eflux_lambda [array] flux error in units of erg/s/cm/cm/AA
    '''
    filterpara = fits.getdata(filterpara)
    nsed = len(Portalsed)
    sed = Portalsed
    sed_flux_lambda = np.zeros(nsed)*np.nan
    sed_eflux_lambda = np.zeros(nsed)*np.nan
    sed_lambda = np.zeros(nsed)*np.nan
    for i, sed_filter in enumerate(sed['sed_filter']): 
         system, filteri = sed_filter.split(':')
         ind_filter = (filterpara['system'] == system) & (filterpara['filter'] == filteri )
         dlambda = filterpara[ind_filter]['dlambda']
         dfreq = filterpara[ind_filter]['dfreq']
         sed_flux = sed['sed_flux'][i]
         sed_eflux = sed['sed_eflux'][i]
         if np.sum(ind_filter) >0:
            sed_lambda[i] =  filterpara[ind_filter]['lambda0'][0]
            sed_flux_lambda[i] = sed_flux*dfreq/dlambda
            sed_eflux_lambda[i] = sed_eflux*dfreq/dlambda
    sed_lambda = sed_lambda*units.um
    sed_flux_lambda = sed_flux_lambda*units.Jy*units.GHz/units.um
    sed_eflux_lambda = sed_eflux_lambda*units.Jy*units.GHz/units.um
    return sed_lambda.to('AA').value, sed_flux_lambda.to('erg/s/cm/cm/AA').value, sed_eflux_lambda.to('erg/s/cm/cm/AA').value


def PortalsedJy2mag(Portalsed, filterpara='filter_parameters.fit'):
    '''
    transform the frequence (nv) and flux_nv of the SED download from http://cdsportal.u-strasbg.fr/ to wave and flux_lambda
    parameters:
    ---------------
    Portalsed: [astropy Table] the SED downloaded from http://cdsportal.u-strasbg.fr/ 
    filterpara: [str] the *.fit file which can be downloaded from http://vizier.u-strasbg.fr/viz-bin/VizieR-4 
    returns:
    ---------------
    mag [array] magnitude
    mag_err [array] magnitude error
    '''
    filterpara = fits.getdata(filterpara)
    nsed = len(Portalsed)
    sed =Portalsed
    mag = np.zeros(nsed)*np.nan 
    mag_e = np.zeros(nsed)*np.nan 
    sed_lambda = np.zeros(nsed)*np.nan 
    for i, sed_filter in enumerate(sed['sed_filter']): 
         system, filteri = sed_filter.split(':')
         ind_filter = (filterpara['system'] == system) & (filterpara['filter'] == filteri )
         dlambda = filterpara[ind_filter]['dlambda']
         dfreq = filterpara[ind_filter]['dfreq']
         sed_flux = sed['sed_flux'][i]
         sed_eflux = sed['sed_eflux'][i]
         Fmag0 =  filterpara[ind_filter]['Fmag0']
         if np.sum(ind_filter) >0:
            mag[i] = -2.5*np.log10(sed_flux/Fmag0[0])
         mag_e[i] = 2.5*sed_eflux/sed_flux/np.log(10)
    return mag, mag_e


def is_color(photband):
    """
    Return true if the photometric passband is actually a color.

    @param photband: name of the photometric passband
    @type photband: string
    @return: True or False
    @rtype: bool
    """
    if '-' in photband.split('.')[1]:
        return True
    elif photband.split('.')[1].upper() in ['M1','C1']:
        return True
    else:
        return False


def eff_wave(photband):
    """
    Returns the effective wavelength of the pass band in angstrom

    @param photband: name of the photometric passband
    @type photband: string
    @return: effective wavelength
    @rtype: float
    """

    data = ascii.read(os.path.join(basedir, 'zeropoints.dat'), comment="\s*#")

    s = np.where(data['photband'] == photband)

    return data['eff_wave'][s][0]


def mag2flux(mag, error, photband):
    """
    Converts a magnitude in a given photband to a flux

    Flux is returned in units of erg/s/cm2/AA
    """

    # Todo: Deal with stromgren colours!!!

    data = ascii.read(os.path.join(basedir, 'zeropoints.dat'), comment="\s*#")

    if not isinstance(photband, str):
        flux, err = np.zeros_like(mag), np.zeros_like(mag)
        for i, (m, e, b) in enumerate(zip(mag, error, photband)):

            s = np.where(data['photband'] == b)
            try:
              F0 = data['Flam0'][s][0]
              zpcor = data['zp_corr'][s][0]
            except:
              _file  = os.path.join(basedir, 'zeropoints.dat')
              print(f'the photband {b} does not exsit in the {_file}')

            f = 10**(-(m-zpcor)/2.5)*F0
            fe = np.log(10) * e / 2.5 * f

            flux[i] = f
            err[i] = fe

    else:
        s = np.where(data['photband'] == photband)
        data = data[s]

        F0 = data['Flam0'][0]
        zpcor = data['zp_corr'][0]

        flux = 10**(-(mag-zpcor)/2.5)*F0
        err = np.log(10) * error / 2.5 * flux

        if hasattr(err, '__iter__'):
            err = np.where(err == 0, flux / 10, err)
        elif err == 0:
            err = flux / 10

    return flux, err



def flux2mag(flux, error, photband):
    """
    Converts a flux in a given photband to a magnitude

    Flux has to be provided in units of erg/s/cm2/AA
    """
    data = ascii.read(os.path.join(basedir, 'zeropoints.dat'), comment="\s*#")

    if not isinstance(photband, str):
        mag, err = np.zeros_like(flux), np.zeros_like(flux)
        for i, (f, e, b) in enumerate(zip(flux, error, photband)):

            s = np.where(data['photband'] == b)

            F0 = data['Flam0'][s][0]
            zpcor = data['zp_corr'][s][0]

            m = -2.5*np.log10(f/F0)+zpcor
            me = 2.5 / np.log(10) * e / f

            mag[i] = m
            err[i] = me

    else:
        s = np.where(data['photband'] == photband)
        data = data[s]

        F0 = data['Flam0'][0]
        zpcor = data['zp_corr'][0]

        mag = -2.5*np.log10(flux/F0)+zpcor
        err = 2.5 / np.log(10) * error / flux

    return mag, err

def Jy2mag(jy, jyerr, photband, filterpara=None):
    '''jy has to be provided in units of Jy
    '''
    if filterpara == None:
       filterpara = fits.getdata(os.path.join(basedir, 'filter_parameters_2021.fit'))
    else:
       filterpara = filterpara
    mag, err = np.zeros_like(jy), np.zeros_like(jy)
    for i, (f, e, b) in enumerate(zip(jy, jyerr, photband)):
         system, filteri = b.split('.')
         ind_filter = (filterpara['system'] == system) & (filterpara['filter'] == filteri )
         Fmag0 =  filterpara[ind_filter]['Fmag0']
         mag[i] = -2.5*np.log10(f/Fmag0)
         err[i] = 2.5*e/f/np.log(10)
    return mag, err



def list_response(name='*', wave_range=(-np.inf, +np.inf)):
    """
    List available response curves.

    Specify a glob string C{name} and/or a wavelength range to make a selection
    of all available curves. If nothing is supplied, all curves will be returned.

    :param name: list all curves containing this string
    :type name: str
    :param wave_range: list all curves within this wavelength range (A)
    :type wave_range: (float, float)
    :return: list of curve files
    :rtype: list of str
    """
    # -- collect all curve files
    if not '*' in name:
        name_ = '*' + name + '*'
    else:
        name_ = name

    curve_files = sorted(glob.glob(os.path.join(basedir, 'transmission_curves', name_.upper())))

    # -- select in correct wavelength range
    curve_files = [os.path.basename(curve_file) for curve_file in curve_files if
                   (wave_range[0] <= eff_wave(os.path.basename(curve_file)) <= wave_range[1])]

    return curve_files


def get_response(photband):
    """
    Returns the response/transmission curve of the provided photometric pass band.
    returns wave, transmission
    Wave in AA
    transmission is unitless or more specifically unit independent.
    """

    transmission_file = os.path.join(basedir, 'transmission_curves', photband)

    if not os.path.isfile(transmission_file):
        raise FileNotFoundError('Can not find transmission curve for photband {}: {}'.format(photband,transmission_file))

    table  = ascii.read(transmission_file, names=['wave', 'flux'])

    return table['wave'], table['flux']


def synthetic_flux(wave, flux, photbands):
    """
    Extract flux measurements from a synthetic SED

    """
    energys = np.zeros(len(photbands))

    for i, photband in enumerate(photbands):

        waver, transr = get_response(photband)
        # -- make wavelength range a bit bigger, otherwise F25 from IRAS has only one Kurucz model point in its
        # wavelength range... this is a bit 'ad hoc' but seems to work.
        region = ((waver[0] - 0.4 * waver[0]) <= wave) & (wave <= (2 * waver[-1]))

        # todo: check if the model needs to be reinterpolated in log scale for some filters
        # # -- if we're working in infrared (>4e4A) and the model is not of high enough resolution (100000 points over
        # # wavelength range), interpolate the model in logscale on to a denser grid (in logscale!)
        # filter_info = filters.get_info()
        # if filter_info['eff_wave'][i] >= 4e4 and 1e5 > sum(region) > 1:
        #     print('%10s: Interpolating model to integrate over response curve' % (photband))
        #     wave_ = np.logspace(np.log10(wave[region][0]), np.log10(wave[region][-1]), 1e5)
        #     flux_ = 10 ** np.interp(np.log10(wave_), np.log10(wave[region]), np.log10(flux[region]), )
        # else:
        #     wave_ = wave[region]
        #     flux_ = flux[region]

        wave_ = wave[region]
        flux_ = flux[region]

        if not len(wave_):
            energys[i] = np.nan
            continue

        # -- perhaps the entire response curve falls in between model points (happends with narrowband UV filters), or
        # there's very few model points covering it
        if (np.searchsorted(wave_, waver[-1]) - np.searchsorted(wave_, waver[0])) < 5:
            wave__ = np.sort(np.hstack([wave_, waver]))
            flux_ = np.interp(wave__, wave_, flux_)
            wave_ = wave__

        # -- interpolate response curve onto model grid
        transr = np.interp(wave_, waver, transr, left=0, right=0)

        # -- integrated flux: different for bolometers and CCDs
        # -- WE WORK IN FLAMBDA
        energys[i] = np.trapz(flux_ * transr * wave_, x=wave_) / np.trapz(transr * wave_, x=wave_)

        # todo: support bolometric filters
        # if photband == 'OPEN.BOL':
        #     energys[i] = np.trapz(flux_, x=wave_)
        # elif filter_info['type'][i] == 'BOL':
        #     energys[i] = np.trapz(flux_ * transr, x=wave_) / np.trapz(transr, x=wave_)
        # elif filter_info['type'][i] == 'CCD':
        #     energys[i] = np.trapz(flux_ * transr * wave_, x=wave_) / np.trapz(transr * wave_, x=wave_)

    # -- that's it!
    return energys


if __name__=="__main__":

    print((mag2flux(13.0, 0.01, '2MASS.H')))
