"""

Collection of procedures that drive the creation of final cubes using
eso recipies. Check out the MUSE reduction manual for details. 

"""
#global import
import numpy as np
import glob 
import subprocess
import time
import os 
import copy
from astropy.io import fits
import socket
import matplotlib.pyplot as mp
import multiprocessing as mp

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

def hawki_combine(filepath, outwcs, memcompress=1, method='sigclip'):
    
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
    
    rowbins = np.linspace(0,ny, num=1+memcompress, dtype=int)
    nrowbins = memcompress
    
    for rbin in range(nrowbins):
        
        minrow = rowbins[rbin]
        maxrow = rowbins[rbin+1]
        nrows = maxrow-minrow
        
        tmpdata = readdata(filelist, ny, nx, minrow, maxrow)

        for yy in range(nrows):
           
           #if yy+minrow != 5359:
           #  continue
           if (yy+minrow) % 100 ==0:
              print('Row {}'.format(yy+minrow))
           
           for xx in range(nx):
               
               _dataarr  = tmpdata[:,yy,xx]
               _dataarr  = _dataarr[_dataarr!=0]
               
               if len(_dataarr)==0:
                  continue
               
               if method=='sigclip':
                 _dataclip = sigclip_cubex(_dataarr, sigma_lower=3, sigma_upper=3, cenfunc=clipmedian, stdfunc=clipstd) 

                 outdata[yy+minrow,xx] = np.mean(_dataclip)
                 outexp[yy+minrow,xx] = len(_dataclip)
                 outstd[yy+minrow,xx] = np.std(_dataclip)/np.sqrt(len(_dataclip))
                 
                 #if xx == 5180 or xx ==4892:
                 #  stop=2
    
    hduout = fits.PrimaryHDU(outdata, header=outhdu[0].header) 
    hduout.writeto('HAWKI_combine_sci.fits', overwrite=True)

    hduout = fits.PrimaryHDU(outstd, header=outhdu[0].header) 
    hduout.writeto('HAWKI_combine_std.fits', overwrite=True)
                 
    hduout = fits.PrimaryHDU(outexp, header=outhdu[0].header) 
    hduout.writeto('HAWKI_combine_expmap.fits', overwrite=True)

    #b /cosma/home/durham/dc-foss1/local/mypython/hawki/hawki_combine.py:68

def workerbootnoise(wrun, filelist, ny, nx, minrow, maxrow, nsamp, median, sigclip, semaphore, sigma_low=3, sigma_up=3, maxiters=5, minnum=4):
    
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
    
    if not os.path.isfile("bootstrapstd_tmpout_batch{}.npz".format(wrun)):
    
       nrows = maxrow-minrow
   
       tmpdata = readdata(filelist, ny, nx, minrow, maxrow) 
       #make space for output
       newstd=np.zeros((ny,nx))

       for yy in range(nrows):
          for xx in range(nx):
             
             _dataarr  = tmpdata[:,yy,xx]
             _dataarr  = _dataarr[_dataarr!=0]
             npix=len(fluxpix)

             if(sigclip) and (npix>4):
                _dataarr = sigclip_cubex(_dataarr, sigma_lower=sigma_lo, sigma_upper=sigma_up, maxiters=maxiters, minnum=minnum, cenfunc=clipmedian, stdfunc=clipstd)
                npix=len(_dataarr)
             
             #bootstrap
             if(npix > 1):
                 #run bootstrap
                 rindex=np.random.randint(npix,size=(nsamp,npix))
                 newstd[xx,yy]=evalstd(evalcent(_dataarr[rindex], axis=1))
             else:
                 newstd[xx,yy]=np.nan
   
       #save output
       np.savez("bootstrapstd_tmpout_batch{}.npz".format(wrun),rowstart=minrow,rowend=maxrow,newstd=newstd)

def hawki_bootnoise(filepath, outwcs, nsamp=1000, nproc=1, method='sigclip'):
    
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
        p = mp.Process(target=workerbootnoise, args=(wrun, filelist, ny, nx, rowmin, rowmax, nsamp, median, sigclip, sema))
        processes.append(p)
        p.start()
    
    for p in processes:
        p.join()    
    
    #reconstruct variance array 
    for wrun in range(nproc):
        thisproc=np.load("bootstrapstd_tmpout_batch{}.npz".format(wrun))
        outstdboot[rowbins[wrun]:rowbins[wrun+1],:]=thisproc['newvar']

    #save to fits file
    hdu=fits.PrimaryHDU(outstdboot)
    hdu.writeto('HAWKI_combine_bootstd.fits',overwrite=True)    

    #clean tmp files
    alltmp=glob.glob("bootstrapstd_tmpout_batch*.npz")
    for tmpfile in alltmp:
        os.remove(tmpfile)

    print('INFO: All done at {}'.format(datetime.datetime.now()))
    

    #b /cosma/home/durham/dc-foss1/local/mypython/hawki/hawki_combine.py:68
