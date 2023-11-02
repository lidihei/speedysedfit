from astropy.table import Table, vstack, hstack
import numpy as np
from astropy import units
from astropy.io import fits
import os
from .filters import Jy2mag, mag2flux, eff_wave
from astropy import units

basedir = os.path.join(os.environ.get('SPEEDYFIT_MODELS', None), '..')

class portalsed():


    def __init__(self, filterpara=None):
       '''
       Transform the SED donwloaded from http://cdsportal.u-strasbg.fr for speedysedfit
       filterpara: [str] the *.fit file which can be downloaded from http://vizier.u-strasbg.fr/viz-bin/VizieR-4 
       '''
       if filterpara is None: filterpara = os.path.join(basedir, 'portal_filter_parameters_2021.fit')
       self.filterpara = filterpara

    def Portalsednv2lambda(self, Portalsed, filterpara=None):
        '''
        Transform the frequence (nv) and flux_nv of the SED download from http://cdsportal.u-strasbg.fr/ to wave and flux_lambda
        parameters:
        ---------------
        Portalsed: [astropy Table] the SED downloaded from http://cdsportal.u-strasbg.fr/ 
        filterpara: [str] the *.fit file which can be downloaded from http://vizier.u-strasbg.fr/viz-bin/VizieR-4 
        returns:
        ---------------
        sed_lambda [array] in units of Angstrom
        sed_flux_lambda [array] flux in units of erg/s/cm/cm/AA'
        sed_eflux_lambda [array] flux error in units of erg/s/cm/cm/AA
        self.dlambda [array] in AA
        self.dfreq [array] in Hz
        '''
        if filterpara is None: filterpara = self.filterpara
        filterpara = fits.getdata(filterpara)
        nsed = len(Portalsed)
        sed = Portalsed
        sed_flux_lambda = np.zeros(nsed)*np.nan
        sed_eflux_lambda = np.zeros(nsed)*np.nan
        sed_lambda = np.zeros(nsed)*np.nan
        sed_dlambda = np.zeros(nsed)*np.nan
        sed_dfreq = np.zeros(nsed)*np.nan
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
                sed_dlambda[i] = dlambda*1.e4
                sed_dfreq[i] = dfreq
        self.dlambda = sed_dlambda
        self.dfreq = sed_dfreq
        sed_lambda = sed_lambda*units.um
        sed_flux_lambda = sed_flux_lambda*units.Jy*units.GHz/units.um
        sed_eflux_lambda = sed_eflux_lambda*units.Jy*units.GHz/units.um
        return sed_lambda.to('AA').value, sed_flux_lambda.to('erg/s/cm/cm/AA').value, sed_eflux_lambda.to('erg/s/cm/cm/AA').value


    def PortalsedJy2mag(self, Portalsed, filterpara=None):
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
        if filterpara is None: filterpara = self.filterpara
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

    def portalsed2phot(self, ra, dec, tab_portsed,
                       radius=3/3600, filename=None,
                       filterpara =None, ra_key = '_RAJ2000', dec_key='_DEJ2000',
                       tabfrom_key = '_tabname',
                       skipphotband=['Cousins.R', 'Cousins.U', 'Cousins.V', 'POSS-II.i', ]):
        '''
        parameters:
        -----------------
        ra: [float]
        dec: [float]
        tab_portsed [Table] the SED downloaded from http://cdsportal.u-strasbg.fr/
        filterpara: [str] the *.fit file which can be downloaded from http://vizier.u-strasbg.fr/viz-bin/VizieR-4
        tabfrom_key: [str] the table name from which the photometry is.
        radius: [float] in degree
        filename: [stri] file out name, if not None write a *.phot file of speedysedfit
        returns:
        tab
        '''
        from PyAstronomy import pyasl
        filterpara = self.filterpara if filterpara is None else filterpara
        filterpara = fits.getdata(filterpara)
        photbands = []
        bands = []
        jys = []
        jyerrs = []
        angdis = pyasl.getAngDist(ra, dec, tab_portsed[ra_key], tab_portsed[dec_key])
        ind = angdis <= radius
        tab = tab_portsed[ind]
        angdis = angdis[ind]*3600
        distance = []
        tabname = []
        waves = []
        lambdas =[]
        dlambdas = []
        dfreqs = []
        tabs = []
        for _i, _b in enumerate(tab['sed_filter']):
            _b = _b.strip(' ')
            if _b in skipphotband: continue
            _sys, _filter = _b.split(':')
            ind_filter = (filterpara['system'] == _sys) & (filterpara['filter'] == _filter )
            dlambda = filterpara[ind_filter]['dlambda'][0]
            dfreq = filterpara[ind_filter]['dfreq'][0]
            sed_lambda =  filterpara[ind_filter]['lambda0'][0]
            self._sys = _sys
            self._filter = _filter
            #_band = f'{_sys}.{_filter}'
            #bands.append(_band)
            if 'PAN-STARRS' in _sys:
                _sys = 'PANSTARRS'
                _filter = _filter.upper()
            elif 'Johnson' in _sys:
                _sys = _sys.upper()
            elif 'GAIA2' in _sys:
                _sys = 'GAIA2'
                if len(_filter)>1: _filter = _filter[0]
                _filter= _filter.upper()
            elif 'GAIA3' in _sys:
                _sys = 'GAIA3E'
                if len(_filter)>1: _filter = _filter[0]
                _filter= _filter.upper()
            elif 'Gaia' in _sys:
                _sys = _sys.upper()
                if len(_filter)>1: _filter = _filter[0]
                _filter= _filter.upper()
            elif 'SDSS' in _sys:
                _filter = _filter.upper()
            elif '2MASS' in _sys:
                _filter = _filter.upper()
                if _filter == 'K': _filter = 'KS'
            _photband = f'{_sys}.{_filter}'
            lambda_eff = eff_wave(_photband)
            if ~np.isfinite(lambda_eff):
                print(f'There is not {_photband} in the zeropoints.dat ; portal filter is {_b}')
                continue
            waves.append(lambda_eff)
            photbands.append(_photband)
            jys.append(tab['sed_flux'][_i])
            jyerrs.append(tab['sed_eflux'][_i])
            tabname.append(tab[tabfrom_key][_i])
            distance.append(angdis[_i])
            dlambdas.append(dlambda)
            dfreqs.append(dfreq)
            lambdas.append(sed_lambda)
            tabs.append(tab[_i])
        mag, magerr = Jy2mag(jys, jyerrs,  bands, filterpara=filterpara)
        flux, err = mag2flux(mag, magerr, np.array(photbands))
        flux = flux*units.erg/units.s/units.A/units.cm**2
        err = err*units.erg/units.s/units.A/units.cm**2
        tab = vstack(tabs)
        dtypes = [('band', '<U20'), ('meas', 'f8'), ('emeas', 'f8'),
                  ('unit', '<U10'), ('distance', 'f8'), ('bibcode', '<U20'),
                  ('flux', 'f8'), ('eflux', 'f8')
                 ]
        unit = ['mag']*len(photbands)
        print(len(mag), len(magerr), len(unit), len(photbands), len(tabname), len(flux), len(distance))
        photometry = Table(data=[photbands, mag, magerr,unit,  distance, tabname, flux, err],
                          names = ['band', 'meas', 'emeas', 'unit', 'distance', 'tabname','flux', 'eflux'],
                          dtype = ['<U20', 'f8', 'f8', '<U10',  'f4', '<U20','f8','f8' ])
        if filename is not None:
            from astropy.io import ascii
            ascii.write(photometry, filename, format='fixed_width', overwrite=True)
        lambdas = np.array(lambdas)*units.um
        dfreqs = np.array(dfreqs)*units.Hz
        dlambdas = np.array(dlambdas)*units.um
        _tab = Table(data=[waves, photbands, mag, magerr, unit,  distance, flux, err, lambdas, dlambdas, dfreqs],
                    names=('wave_eff', 'photbands', 'mag', 'magerr', 'unit', 'distance', 'flux', 'eflux', 'lambda0', 'dlambda', 'freq'))
        tab = hstack([tab, _tab])
        return tab
