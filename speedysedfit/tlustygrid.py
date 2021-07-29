import numpy as np
import os
from astropy.io import fits
from .tlustysed import rsed
from .broadsed import vgconv

clight = 299792.458


def get_primaryhdu(LOGZ='G', REF='TLUSTY2002', Z = 0.):
    n = np.zeros((1, 2), dtype=np.int32)
    primary_hdu = fits.PrimaryHDU(n)
    header = primary_hdu.header
    header['LOGZ'] = (LOGZ, 'Metallicity wrt Solar')
    header['REF'] =(REF, 'Table reference')
    header['WAVUNIT'] = ('angstrom', 'wavelength units')
    header['FLXUNIT'] = ('erg/s/cm2/A', 'flux units')
    header['Z'] = Z
    return primary_hdu

def get_tablehdu(fname, R =300, wave_start= 9.09e+01, wave_end = 1.60e+06):
    '''
    R: [float] is resolution
    '''
    extname = os.path.basename(fname)
    _extname = extname.replace('B', '')
    TEFF = float(_extname[1:6])
    LOGG = float(_extname[7:10])/100
    Vt = float(_extname[11:13])
    vfhmw = clight/R
    wcor, fcor = rsed(fname)
    w1, fbin = vgconv(wcor, fcor,vfhmw, ppr=1)
    ind =( w1 >= wave_start) & (w1 <=wave_end )
    _w1 = np.array(w1[ind], dtype=np.float32)
    _fbin = np.array(fbin[ind], dtype = np.float32)
    cw = fits.Column(name='wavelength', array=_w1, format='E')
    cf = fits.Column(name='flux', array=_fbin, format='E')
    table_hdu = fits.BinTableHDU.from_columns([cw,cf])
    header = table_hdu.header
    header['EXTNAME'] = (extname, 'name of the extension')
    header['TEFF']    = (TEFF, 'Effective temperature (K)')
    header['LOGG']    = (LOGG, 'Log(g)')
    header['Vt']    = (Vt, 'microturbulence velocity (km/s)')
    table_hdu.header = header
    return table_hdu

if __name__=="__main__":
   import joblib
   from tqdm import tqdm
   filelist = glob('/share/lijiao/TLUSTY_grid/wangluqian/Tlusty/G*.flux')
   dirout = '/share/lijiao/speedyfit/modelgrids/'
   primary_hdu = get_primaryhdu(LOGZ='G', REF='TLUSTY2002', Z = 0.)
   table_hdus = joblib.Parallel(n_jobs=10, )(joblib.delayed(get_tablehdu)(fname, R =300, wave_start= 9.09e+01, wave_end = 1.60e+06) for fname in tqdm(filelist))
   hdul = fits.HDUList([primary_hdu] + table_hdus)
   hdul.writeto(dirout+'TlustyO_G_v10.fits', overwrite=True)
