import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
import re
from .model import get_grid_file

class readspecgrids():

    def __init__(self):
       '''
       read and dowload the model spectra, such as TLUSTY, ATLAS9 fluxes (Castelli flux grids);
       obtain parameters table from file name
       '''

    @staticmethod
    def get_lam(R, lam_start, lam_end, N=3):
        ''' get log10(lambda) with a proper sample interval
        parameters:
        ---------------
        R: [float] spectral resolution (BFOSC: R ~ 1600)
        lam_start: [float]: the start of wavelength
        lam_end: [float]: the end of wavelength
        N: [int] oversampling, (typical oversampling: N=3 --> 3pixel=FWHM)
           R = lambda/FWHM --> (log10(lambda))' =1/(lambda*ln(10)) --> FWHM = 1/(R*ln(10))
        returns:
        wave: [array] in angstrom
        '''
        deltax = 1/R/N/np.log(10)
        log10lam = np.arange(np.log10(lam_start), np.log10(lam_end), deltax)
        return 10**log10lam

    def creat_regli_for_interpolate(self, keys= ['TEFF', 'LOGG'],grid='tlustyO',  wavelength=None, verbose=True, fout=None):
        ''' create a grid of relig for interpolating flux
        grid [stri]: the grid name in the grid_description.yaml, e.g. 'tlustyO', 'kurucz', 'tmap'
        note: transform the teff into log10(teff); some points might not be interpolated as the neigbor grids points are missed
        returns:
        interpolater

        example:
        flux = interpolater.interpn([np.log10(30618), 4])
        '''
        from regli import Regli
        import joblib
        gridfilename = get_grid_file(integrated=False, grid=grid)
        if fout is None: fout = f'{gridfilename[:-5]}_regli.z'
        hdus = fits.open(gridfilename)
        N = len(hdus)-1
        dim = len(keys)
        parameters = np.zeros((N, dim))
        if wavelength is None:
           wave1 = hdus[1].data['wavelength']
           wave_start, wave_end = np.min(wave1), np.max(wave1)
           wavelength = self.get_lam(400, wave_start, wave_end, N=3)
        fluxs = np.zeros((N, len( wavelength)))
        for _i, hdu in enumerate(hdus[1:]):
            fluxs[_i] = np.interp(wavelength, hdu.data['wavelength'], hdu.data['flux'], left=np.nan, right=np.nan,)
            for _j, key in enumerate(keys):
                parameters[_i, _j] = hdu.header[key]
                if key == 'TEFF': parameters[_i, _j] = np.log10(parameters[_i, _j])
        interpolater = Regli.init_from_flats(parameters, verbose= verbose)
        interpolater.set_wave(wavelength)
        interpolater.set_values(fluxs)
        dumpdic = {'interpolater' : interpolater,
                   'wavelength': wavelength,
                   'code': os.path.abspath(__file__)
                  }
        joblib.dump(dumpdic, fout)
        return interpolater

    def read_castelli_fluxgrids(self, fname):
        '''read the spectrum of Castelli flux grids (https://wwwuser.oats.inaf.it/castelli/grids.html)
        #description see http://kurucz.cfa.harvard.edu/grids.html
        parameters:
        -------------
        fname: [str]
        returns:
        wave: [array] wavelength in A
        nv: [array] frequency nu (s^-1)
        flambda:[array]  (erg/cm^2/s/Hz/ster)
        fnv: [array] erg/cm^2/s/A/ster)
        '''
        ff = open(fname, 'r')
        lines = ff.readlines()
        nn = len(lines) - 3
        def floatstr(stri):
            try: num = np.float32(stri)
            except:
               if '-' in stri:
                   stri= stri.replace('-', 'E-')
               elif '+' in stri:
                   stri= stri.replace('-', 'E+')
               num = np.float(stri)
            return num
        #wv, nv, fnv, fcont, fnorm = np.zeros(5, nn, dtype=np.float32)
        data = np.zeros((nn, 5), dtype=np.float32)
        for i, line in enumerate(lines[2:]):
            if (len(line) >5):
               _line = np.zeros(5, dtype=np.float32)
               line =  line[9:].split()
               _line[0] = np.float32(line[0])
               _line[1] = np.float32(line[1])
               _line[2] = floatstr(line[2])
               _line[3] = floatstr(line[3])
               _line[4] = floatstr(line[4])
               data[i] = _line
        wave, nv, fnv, fcont, fnorm = data.T
        wave = wave*10
        c = 2.9979246e18
        flambda = 4*fnv*c/wave**2
        return wave, nv, flambda, fnv

    def fname2paratab_castelli(self, fname, subgrid='fp00k2c125odfnew'):
        '''
        transform the castelli file name to parameter tab (the castelli grids can be downloaded from https://wwwuser.oats.inaf.it/castelli/grids.html)
        paramters:
        ----------------
        fname: [str] e.g. 'fp00t16000g45k2odfnew.dat'
        returns:
        ----------------
        tab [astropy table], include ('[MH]', 'Teff', 'logg', 'Ki', 'mixing_length', 'fname', 'subgrid')
        '''
        mixing_length = np.asarray(re.findall(r'c\d+', subgrid)[0][1:], float)/100
        if fname[0:2] == 'fp': factor = np.array([0.1, 1., 0.1, 1.])
        elif fname[0:2] == 'fm': factor = np.array([-0.1, 1., 0.1, 1.])
        MH, Teff, logg, Ki = np.asarray(re.findall(r'\d+', fname), float) * factor
        tab = Table(data=[[MH], [Teff], [logg], [Ki], [mixing_length], [fname], [subgrid]], names = ('[MH]', 'Teff', 'logg', 'Ki', 'mixing_length', 'fname', 'subgrid'))
        return tab

    def download_grids_castelli(self, teff, logg,
                   downloaddir = '../fp00k2c125odfnew',
                   prefix = 'https://wwwuser.oats.inaf.it/castelli/grids/',
                   subgrid = 'gridp00k2odfnew',
                   suffix = 'fp00t{}g{}k2odfnew.dat',
                   dir_current =None):
        '''
        Using wget to download castelli grids (https://wwwuser.oats.inaf.it/castelli/grids.html), and store it in downloaddir
        dir_current [str] the current working directory
        '''
        fname = suffix.format(teff, logg)
        print(fname)
        url = f'{prefix}/{subgrid}/{fname}'
        os.chdir(downloaddir)
        os.system(f'wget {url}')
        os.chdir(dir_current)

    def read_tlusty_fluxgrids(self, fname):
        '''read sed file of Tlusty dowloaded from http://tlusty.oca.eu/Tlusty2002/OS02-SED.html
        parameters:
        --------------------
        fname: [str] e.g. BG28000g475v2.flux
        returns:
        -------------------
        wavelength: [array] in units of angstrom
        flux: [array] in units of erg/cm^2/s/A
        '''
        re_dbl_fort = re.compile(r'(\d*\.\d+)[dD]([-+]?\d+)')
        re_dbl_fort1 = re.compile(r'(\d*\.\d+)[-](\d+)')
        re_dbl_fort2 = re.compile(r'(\d*\.\d+)[+](\d+)')
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
