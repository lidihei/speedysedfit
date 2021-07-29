import numpy as np
import re
clight = 299792.458


def vgconv(xinput,yinput,fwhm, ppr=None):
  """convolution with a Gaussian in log lambda scale
  for a constant resolving power
  form [synple.py](https://github.com/callendeprieto/synple)
  Parameters
  ----------
  xinput: numpy float array
      wavelengths 
  yinput: numpy array of floats
      fluxes
  fwhm: float
      FWHM of the Gaussian (km/s)
  ppr: float, optional
      Points per resolution element to downsample the convolved spectrum
      (default None, to keep the original sampling)

  Returns
  -------
  x: numpy float array
      wavelengths after convolution, will be a subset of xinput when that is equidistant
      in log lambda, otherwise a subset of the resampled version
  y: numpy array of floats
      fluxes after convolution

  """
  #resampling to ln(lambda) if need be
  xx = np.diff(np.log(xinput))
  if np.max(xx) - np.min(xx) > 1.e-7:  #input not equidist in loglambda
    nel = len(xinput)
    minx = np.log(xinput[0])
    maxx = np.log(xinput[-1])
    x = np.linspace(minx,maxx,nel)
    step = x[1] - x[0]
    x = np.exp(x)
    y = np.interp( x, xinput, yinput)
    #y = interp_spl( x, xinput, yinput)
  else:
    x = xinput
    y = yinput
    step = np.log(xinput[1])-np.log(xinput[0])
  fwhm = fwhm/clight # inverse of the resolving power
  sigma=fwhm/2.0/np.sqrt(-2.0*np.log(0.5))
  npoints = 2*int(3*fwhm/2./step)+1
  half = npoints * step /2.
  xx = np.linspace(-half,half,npoints)
  kernel = np.exp(-(xx-np.mean(xx))**2/2./sigma**2)
  kernel = kernel/np.sum(kernel)
  
  y = np.convolve(y,kernel,'valid')
  edge = int(npoints/2)
  x = x[edge:-edge]

  #print(xinput.size,x.size,y.size)

  if ppr != None:
    fac = int(fwhm / step / ppr)
    #print(fwhm,step,ppr,fac)
    subset = np.arange(x.size / fac, dtype=int) * fac 
    x = x[subset]
    y = y[subset]
  return x, y
