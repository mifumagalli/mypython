"""

Collection of procedures that drive the creation of final cubes using
eso recipies. Check out the MUSE reduction manual for details. 

"""
#global import
import numpy as np
import glob 
import subprocess
import time
import os,sys
import copy
from astropy.io import fits, ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import wcs
from astropy.table import Table
from astropy.nddata import NDData
import socket
import matplotlib.pyplot as mp
import multiprocessing as mp
import photutils.aperture as apphot
import scipy.interpolate as sinterp
from photutils import psf

try:
  from bottleneck import nanstd as clipstd
except:
  from numpy import std as clipstd
try:
  from bottleneck import nanmedian as clipmedian 
except:
  from numpy import median as clipmedian

def sigclip_cubex(array, sigma_lower=3, sigma_upper=3, maxiters=1000, minnum=3, cenfunc=None, stdfunc=None, return_bounds=False):
    
    '''
    Function that replicates the sigma clipping scheme used in CubEx
    Written for speed, no checks are done on the consistency of the inputs
    
    Returns clipped array (outliers removed), if return_bounds is True returns 
    also the lower and upper bounds.
    
    '''
    
    npoints = len(array)
    niter   = 0
    stopit  = 0
    
    if cenfunc==None:
       from numpy import nanmedian as cenfunc
    if stdfunc==None:
       from numpy import nanstd as stdfunc   
    
    while(niter<maxiters and npoints>minnum and stopit==0):
       bounds = cenfunc(array)+np.array((-1*sigma_lower,sigma_upper))*stdfunc(array)
       array = array[(array>bounds[0]) & (array<bounds[1])]
       if len(array)<npoints:
          npoints=len(array)
          niter += 1
       else:
          stopit = 1      
    
    if return_bounds:
      return array, bounds[0], bounds[1]     
    else:
      return array
    
def readdata(filelist, ny, nx, minrow, maxrow):
    
    nfiles = len(filelist)
    outstruct = np.zeros((nfiles, maxrow-minrow, nx))
    
    for cnt, ff in enumerate(filelist):
        
        tmphdu = fits.open(ff)
        
        tmpdata = np.zeros((ny, nx))
        minx = tmphdu[0].header['COMIN1']-1
        miny = tmphdu[0].header['COMIN2']-1
        maxx = tmphdu[0].header['COMIN1']-1+tmphdu[0].header['NAXIS1']
        maxy = tmphdu[0].header['COMIN2']-1+tmphdu[0].header['NAXIS2']
        
        tmpdata[miny:maxy,minx:maxx] = tmphdu[0].data
        
        outstruct[cnt,...] = tmpdata[minrow:maxrow,:]
    
    return outstruct

def hawki_combine(filepath, outwcs, memcompress=2, nproc=10, median=False, sigclip=True):
    
    #First read outwcs file and prepare the structure
    outhdu = fits.open(outwcs)
    
    outdata = np.zeros_like(outhdu[0].data)+np.nan
    outstd  = np.zeros_like(outhdu[0].data)+np.nan
    outexp  = np.zeros_like(outhdu[0].data)
    
    ny, nx = np.shape(outdata)
    
    print('Read output WCS file with {}x{} pixels'.format(nx,ny))
    
    #Read file list and find number of files
    filelist = glob.glob(filepath+'/*.fits')
    nfiles = len(filelist)

    print('Found {} files to be combined'.format(nfiles))
    
    nrowbins = memcompress*nproc
    rowbins = np.linspace(0,ny, num=nrowbins+1, dtype=int)

    #loop over processors and start parallel function
    processes = []
    sema = mp.Semaphore(nproc)
    for wrun in range(nrowbins):
        sema.acquire()
        p = mp.Process(target=workercombine, args=(wrun, filelist, ny, nx, rowbins[wrun], rowbins[wrun+1], median, sigclip, sema))
        processes.append(p)
        p.start()
    
    for p in processes:
        p.join()    
    
    #reconstruct variance array 
    for wrun in range(nrowbins):
        thisproc=np.load("hwkcombine_tmpout_batch{}.npz".format(wrun))
        outdata[rowbins[wrun]:rowbins[wrun+1],:]=thisproc['newsci']
        outstd[rowbins[wrun]:rowbins[wrun+1],:]=thisproc['newstd']
        outexp[rowbins[wrun]:rowbins[wrun+1],:]=thisproc['newexp']
        
  
    hduout = fits.PrimaryHDU(outdata, header=outhdu[0].header) 
    hduout.writeto('HAWKI_combine_sci.fits', overwrite=True)

    hduout = fits.PrimaryHDU(outstd, header=outhdu[0].header) 
    hduout.writeto('HAWKI_combine_std.fits', overwrite=True)
                 
    hduout = fits.PrimaryHDU(outexp, header=outhdu[0].header) 
    hduout.writeto('HAWKI_combine_expmap.fits', overwrite=True)

    #clean tmp files
    alltmp=glob.glob("hwkcombine_tmpout_batch*.npz")
    for tmpfile in alltmp:
        os.remove(tmpfile)


    #b /cosma/home/durham/dc-foss1/local/mypython/hawki/hawki_combine.py:68

def workerbootnoise(wrun, filelist, ny, nx, minrow, maxrow, nsamp, median, sigclip, semaphore, sigma_lo=3, sigma_up=3, maxiters=5, minnum=4):
    
    #select evaluator
    if(median): 
         from numpy import median as evalcent
         if wrun ==0:
             print("INFO: Using numpy for statistical functions")
    else:
         from numpy import mean as evalcent
         if wrun ==0:
             print("INFO: Using numpy for statistical functions")
          
    from numpy import std as evalstd
    from numpy import nanmedian as clipmedian
    from numpy import nanstd as clipstd
    
    if not os.path.isfile("bootstrapstd_tmpout_batch{}.npz".format(wrun)):
    
       nrows = maxrow-minrow
   
       tmpdata = readdata(filelist, ny, nx, minrow, maxrow) 
       #make space for output
       newstd=np.zeros((maxrow-minrow,nx))
       
       print('Worker {}. Data read'.format(wrun))

       for yy in range(nrows):

          if yy % int(nrows/5) == 0:
                print('Worker {}. Step {}'.format(wrun, yy))

          for xx in range(nx):
                      
             _dataarr  = tmpdata[:,yy,xx]
             _dataarr  = _dataarr[_dataarr!=0]
             npix=len(_dataarr)

             if(sigclip) and (npix>4):
                _dataarr = sigclip_cubex(_dataarr, sigma_lower=sigma_lo, sigma_upper=sigma_up, maxiters=maxiters, minnum=minnum, cenfunc=clipmedian, stdfunc=clipstd)
                npix=len(_dataarr)
             
             #bootstrap
             if(npix > 1):
                 #run bootstrap
                 rindex=np.random.randint(npix,size=(nsamp,npix))
                 newstd[yy,xx]=evalstd(evalcent(_dataarr[rindex], axis=1))
             else:
                 newstd[yy,xx]=np.nan
   
       #save output
       np.savez("bootstrapstd_tmpout_batch{}.npz".format(wrun),rowstart=minrow,rowend=maxrow,newstd=newstd)
     
    else:
       print('Already done...Skipping.')
     
    sys.stdout.flush()
    semaphore.release()  
       
     
     
def workercombine(wrun, filelist, ny, nx, minrow, maxrow, median, sigclip, semaphore, sigma_lo=3, sigma_up=3, maxiters=5, minnum=4):
    
    #select evaluator
    if(median): 
         from numpy import median as evalcent
         if wrun ==0:
             print("INFO: Using numpy for statistical functions")
    else:
         from numpy import mean as evalcent
         if wrun ==0:
             print("INFO: Using numpy for statistical functions")
          
    from numpy import std as evalstd
    from numpy import nanmedian as clipmedian
    from numpy import nanstd as clipstd
    
    if not os.path.isfile("hwkcombine_tmpout_batch{}.npz".format(wrun)):
    
       nrows = maxrow-minrow
   
       tmpdata = readdata(filelist, ny, nx, minrow, maxrow) 
       #make space for output
       newsci=np.zeros((maxrow-minrow,nx))+np.nan
       newstd=np.zeros((maxrow-minrow,nx))+np.nan
       newexp=np.zeros((maxrow-minrow,nx))+np.nan
       
       print('Worker {}. Data read'.format(wrun))

       for yy in range(nrows):

          if yy % int(nrows/5) == 0:
                print('Worker {}. Step {}'.format(wrun, yy))

          for xx in range(nx):
                      
             _dataarr  = tmpdata[:,yy,xx]
             _dataarr  = _dataarr[_dataarr!=0]
             npix=len(_dataarr)
             
             if npix==0:
                continue
             
             if(sigclip) and (npix>4):
                _dataarr = sigclip_cubex(_dataarr, sigma_lower=sigma_lo, sigma_upper=sigma_up, maxiters=maxiters, minnum=minnum, cenfunc=clipmedian, stdfunc=clipstd)
               
             newsci[yy,xx] = np.mean(_dataarr)
             newexp[yy,xx] = len(_dataarr)
             newstd[yy,xx] = np.std(_dataarr)/np.sqrt(len(_dataarr))
                   
       #save output
       np.savez("hwkcombine_tmpout_batch{}.npz".format(wrun),rowstart=minrow,rowend=maxrow,newstd=newstd, newsci=newsci, newexp=newexp)
       print('Worker {}. DONE'.format(wrun))
     
    else:
       print('Already done...Skipping.')
     
    sys.stdout.flush()
    semaphore.release()  

     
       
def hawki_bootnoise(filepath, outwcs, nsamp=1000, nproc=1, median=False, sigclip=True):
    
    #First read outwcs file and prepare the structure
    outhdu = fits.open(outwcs)
    outstdboot  = np.zeros_like(outhdu[0].data)+np.nan
    
    ny, nx = np.shape(outstdboot)
    
    print('Read output WCS file with {}x{} pixels'.format(nx,ny))
    
    #Read file list and find number of files
    filelist = glob.glob(filepath+'/*.fits')
    nfiles = len(filelist)

    print('Found {} files to be combined'.format(nfiles))
    
    rowbins = np.linspace(0,ny, num=1+nproc, dtype=int)
    nrowbins = nproc
    
    #loop over processors and start parallel function
    processes = []
    sema = mp.Semaphore(nproc)
    for wrun in range(nproc):
        sema.acquire()
        p = mp.Process(target=workerbootnoise, args=(wrun, filelist, ny, nx, rowbins[wrun], rowbins[wrun+1], nsamp, median, sigclip, sema))
        processes.append(p)
        p.start()
    
    for p in processes:
        p.join()    
    
    #reconstruct variance array 
    for wrun in range(nproc):
        thisproc=np.load("bootstrapstd_tmpout_batch{}.npz".format(wrun))
        outstdboot[rowbins[wrun]:rowbins[wrun+1],:]=thisproc['newstd']
    #save to fits file
    hdu=fits.PrimaryHDU(outstdboot)
    hdu.writeto('HAWKI_combine_bootstd.fits',overwrite=True)    

    #clean tmp files
    alltmp=glob.glob("bootstrapstd_tmpout_batch*.npz")
    for tmpfile in alltmp:
        os.remove(tmpfile)
    


def hawki_combinephotom(photcat, scifile, stdfile, magcol='Kmag', ZP=28, addfile=None):
    
    if 'K' in magcol:
       ABoffset= 1.85
    elif 'H' in magcol:
       ABoffset= 1.39
    elif 'J' in magcol:
       ABoffset= 0.91
    else:       
       ABoffset= 0
    
    #Read photom catalogue
    phottab = fits.open(photcat)[1].data
    
    okstars = phottab[magcol]>0
    raphotcat  = phottab['RA'][okstars]
    decphotcat = phottab['Dec'][okstars]
    magphotcat = phottab[magcol][okstars]
    
    scihdu = fits.open(scifile)
    head = scihdu[0].header
    data = scihdu[0].data
    datawcs = wcs.WCS(head)
    
    if ZP is None:
      ZP = head['PHOTZP']
    
    positions = SkyCoord(raphotcat, decphotcat, frame='icrs', unit=(u.deg, u.deg))
    apertures = apphot.SkyCircularAperture(positions, r=2.2 *u.arcsec)
    apertures_pix = apertures.to_pixel(datawcs)
    
    photom = apphot.aperture_photometry(data, apertures_pix)
    maghawki = -2.5*np.log10(photom['aperture_sum'])+ZP
    
    zpoffset = np.median(magphotcat-maghawki)
    fscale = 10**(-0.4*zpoffset)
    
    head['ZPVEGA'] = ZP
    head['ZPAB']   = ZP + ABoffset
    
    data *= fscale
    
    hduout = fits.PrimaryHDU(data, header=head)
    hduout.writeto(scifile, overwrite=True)
    
    stdhdu = fits.open(stdfile)
    
    stddata = stdhdu[0].data
    
    hduout = fits.PrimaryHDU(stddata, header=head)
    hduout.writeto(stdfile, overwrite=True)
    
    #import matplotlib.pyplot as mp
    #mp.scatter(magphotcat, magphotcat-maghawki)
    #mp.show()
    
def hawki_buildepsf(scifile, stars, pscale=0.06, size=69):
    
    oversampling = 0.106/pscale
    
    #read star positions
    strfile = ascii.read(stars)
    
    data = fits.open(scifile)[0].data
    nddata = NDData(data=data)
    
    strdata = psf.extract_stars(nddata, strfile, size=50)
    
    import matplotlib.pyplot as mp
    from astropy.visualization import simple_norm
    nrows = 5
    ncols = 5
    
    fig, ax = mp.subplots(nrows=nrows, ncols=ncols, figsize=(20,20), squeeze=True)
    ax = ax.ravel()
    for i in range(len(strdata)):
        norm = simple_norm(strdata[i], 'log', percent=99.)
        ax[i].imshow(strdata[i], norm=norm, origin='lower', cmap='viridis')
    
    mp.show()        
    
    epsf_builder = psf.EPSFBuilder(oversampling=2, shape=91, maxiters=20, progress_bar=True, smoothing_kernel='quadratic' )
    epsf, fitted_stars = epsf_builder(strdata)
    
    pos_oversamp = (np.arange(91)-45)*(0.106/2)
    pos_requested = (np.arange(size)-np.floor(size/2.))*pscale
    
    finterp = sinterp.interp2d(pos_oversamp, pos_oversamp, epsf.data, kind='linear')
    epsf_interp = finterp(pos_requested, pos_requested)
    epsf_interp[epsf_interp<0] = 0
        
    norm = simple_norm(epsf.data, 'log', percent=99.)
    mp.imshow(epsf.data, norm=norm, origin='lower', cmap='viridis')
    mp.show()
    mp.imshow(epsf_interp, norm=norm, origin='lower', cmap='viridis')
    mp.show()
    
    hduout = fits.PrimaryHDU(epsf_interp/np.nansum(epsf_interp))
    hduout.header['PSCALE' ] = pscale
    hduout.header['OVERSAMP'] = oversampling
    
    hduout.writeto('HAWKI_epsf.fits', overwrite=True)
    
    #b /cosma/home/durham/dc-foss1/local/mypython/hawki/hawki_combine.py:68
