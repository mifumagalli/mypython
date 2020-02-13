import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits 
import sep
import datetime
from scipy import ndimage
from scipy import signal
from scipy import stats
from scipy import interpolate as inter
import multiprocessing as mp
import os, sys
import glob

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

def evaluatenoise(wrun,nx,ny,nexp,nsamp,allexposures,allmasks,masks,median,sigclip,semaphore,sigma_lo=3,sigma_up=3,maxiters=5,minnum=4):
    
    """
    Utility function that evaluates boostrap noise
    
    """
    
    if not os.path.isfile("boostrapvar_tmpout_w{}.npz".format(wrun)):
       
       #now load the exposures
       dataexp=np.zeros((nexp,nx,ny))
    
       for exp in range(nexp):
           dataexp[exp]=fits.open(allexposures[exp])[1].data[wrun,:,:]
           #If requested apply mask on the fly
           if(masks):
             tmpmask=fits.open(allmasks[exp])[0].data[:,:]
             dataexp[exp,(tmpmask==0)] = np.nan

       print('PROG: Slice {}. All data loaded, now running job...'.format(wrun))        
           
       #make space for output
       newvar=np.zeros((nx,ny))

       #select evaluator
       if(median): 
           try:
             from bottleneck import nanmedian as evalcent 
             if wrun ==0:
                print("INFO: Using bottleneck for statistical functions")
           except:
             from numpy import median as evalcent
             if wrun ==0:
                print("INFO: Using numpy for statistical functions")
       else:
           try:
             from bottleneck import nanmean as evalcent
             if wrun ==0:
                print("INFO: Using bottleneck for statistical functions")
           except:
             from numpy import mean as evalcent
             if wrun ==0:
                print("INFO: Using numpy for statistical functions")
             
    
       try:
         from bottleneck import nanvar as evalvar
       except:
         from numpy import var as evalvar
    
       if(sigclip):
           try:
             from bottleneck import nanstd as clipstd
           except:
             from numpy import std as clipstd
           try:
             from bottleneck import nanmedian as clipmedian 
           except:
             from numpy import median as clipmedian
        
       #giant loop on pixels - necessary otherwise the memory goes off the roof
       for xx in range(nx):
           for yy in range(ny):
                   
               #ingest flux 
               fluxpix=dataexp[:,xx,yy]
               fluxpix=fluxpix[np.isfinite(fluxpix)]
               npix=len(fluxpix)
               #print(ww,npix)
               if(sigclip) and (npix>4):
                  #dummy, low, high = stats.sigmaclip(fluxpix, low=sigma_lo, high=sigma_up)
                  #fluxpix=fluxpix[(fluxpix>low) & (fluxpix<high)]
                  fluxpix = sigclip_cubex(fluxpix, sigma_lower=sigma_lo, sigma_upper=sigma_up, maxiters=maxiters, minnum=minnum, cenfunc=clipmedian, stdfunc=clipstd)
                  npix=len(fluxpix)
               
               #bootstrap
               if(npix > 1):
                   #run bootstrap
                   rindex=np.random.randint(npix,size=(nsamp,npix))
                   newvar[xx,yy]=evalvar(evalcent(fluxpix[rindex], axis=1))
               else:
                   newvar[xx,yy]=np.nan
                       
       #save output
       np.savez("boostrapvar_tmpout_w{}.npz".format(wrun),wstart=wrun,wend=wrun,newvar=newvar)
    
    else:
       print('PROG: Slice {} already done. Skipping...'.format(wrun))     
    
    sys.stdout.flush()
    semaphore.release()

def bootstrapnoise(cubes,masks=None,nsamp=10000,outvar="bootstrap_variance.fits",nproc=50,median=False,sigclip=False):

    """
    
    Take a list of exposures and estimate variance with boostrap
    at pixel level. Generates an output bootstrap cube.

    cubes -> list of input cubes
    masks -> list of associated masks
    nsamp -> number of samples to draw [ideally 10000]
    outvar -> where to store output
    nproc -> number of proc to run this over 
    median -> if True, switches to median estimator 
    
    """

    print('Start at {}'.format(datetime.datetime.now()))

    #load the exposures and stash them in pointer stack
    allexposures=[]
    allmasks=[]
    nexp=0
    for exp in open(cubes):
        allexposures.append(exp.strip())
        nexp=nexp+1
    print('INFO: Found {} exposures'.format(nexp))

    
    if(masks):
        for exp in open(masks):
            allmasks.append(exp.strip())
    print('INFO: Found {} masks'.format(nexp))
    
    #find format of data and create empty var
    nw,nx,ny=fits.open(allexposures[0])[1].data.shape
    print('INFO: Data format {} {} {}'.format(nw,nx,ny))
    
    #now prepare batches for parallel run
    print('INFO: Running on {} proc '.format(nproc))

    #loop over processors and start parallel function
    processes=[]
    sema = mp.Semaphore(nproc)
    for wrun in range(nw):
        #make sure does not run over index
        sema.acquire() 
        p=mp.Process(target=evaluatenoise,args=(wrun,nx,ny,nexp,nsamp,allexposures,allmasks,masks,median,sigclip,sema))
        processes.append(p)
        p.start()
        
    for p in processes:
        p.join()
    
    #reconstruct variance array 
    allvar=np.zeros((nw,nx,ny), dtype=np.float32)
    for wrun in range(nw):
        thisproc=np.load("boostrapvar_tmpout_w{}.npz".format(wrun))
        allvar[wrun,:,:]=thisproc['newvar']

    #save to fits file
    hdu=fits.PrimaryHDU(allvar)
    hdu.writeto(outvar,overwrite=True)    

    #clean tmp files
    alltmp=glob.glob("boostrapvar_tmpout_w*.npz")
    for tmpfile in alltmp:
        os.remove(tmpfile)

    print('INFO: All done at {}'.format(datetime.datetime.now()))


def applybootnoise(cube, bootcube, outcube, varscale=1.):
    
    """
    
    Get and input cube (.e.g from cubex pipeline) and a bootstrap
    pre-computed variance cube to place as the new variance in the
    input cube. An additional global scaling factor to the variance 
    can be applied.
    
    cube -> the datacube including pipeline variance
    bootcube -> the variance datacube from bootstrap resampling
    outcube -> name of the output cube
    varscale -> multiplicative factor to rescale variance before saving
    
    """
    #open the data
    data=fits.open(cube)
    boot= fits.open(bootcube)
    
    data[2].data = boot[0].data
    
    if varscale != 1:
       data[2].data *= varscale
       
    data.writeto(outcube, overwrite=True)   
    print('All data saved!')
    

def globalscalenoise(cube, outcube, memmap=False):

    """
    Read a cube and compute a global scaling factor to the variance
    such that the variance in the data is consistent with the 
    variance extension. The scaling factor is calculated in regions
    of the spectrum free from skylines.
    
    """

    #Read cube 
    hdu = fits.open(cube, memmap=memmap)

    data = hdu[1].data
    std = np.sqrt(hdu[2].data)

    nz, ny, nx = np.shape(data)
    
    wave = hdu[1].header['CRVAL3']+np.arange(nz)*hdu[1].header['CD3_3']
  
    #compress into image
    image=np.nanmedian(data, axis=0)
    nx,ny=image.shape

    #mask edges
    edges=np.isfinite(image)
    badmask=np.zeros((nx,ny))+1
    badmask[edges]=0.0
    badmask=ndimage.gaussian_filter(badmask,1.5)
    badmask[np.where(badmask > 0)]=1.0

    #mask sources
    bkg = sep.Background(image,mask=badmask)    
    thresh = 1.5 * bkg.globalrms
    segmap = np.zeros((nx,ny))
    objects,segmap=sep.extract(image,thresh,segmentation_map=True,
                               minarea=10,clean=True,mask=badmask)
    badmask[np.where(segmap > 0)]=1.0
  
    tonan = (badmask>0)
    badmask[tonan] = np.nan
    badmask[np.logical_not(tonan)] = 1
  
    mask3d = np.broadcast_to(badmask,(nz,)+badmask.shape)
    fsig = data/std*mask3d

    fsig_1d = np.nanstd(fsig, axis=(1,2))
    
    okwave = ((wave>4700) & (wave<5800)) | ((wave>6600) & (wave<6800))
    
    #Average does not like nans, mask them out
    global_offset = np.nanmedian(fsig_1d[okwave])

    hdu[2].data *= global_offset**2
    hdu[2].header['VARSCALE'] = global_offset**2
    
    hdu.writeto(outcube, overwrite=True)

    

def rescalenoise(cube,rescaleout="rescale_variance.txt",outvar="CUBE_rmsvar.fits",cut=10,smooth=1,block=65,disp=0.07,sthre=1.0,bootstrap=None,expmap=None,expmap_range=[0,0],memmap=False,savechecks=None):
    
    """

    Get an input cube (e.g. from cubex pipeline) and compute rescaling
    factor to match the data rms

    cube -> the datacube including pipeline variance
    rescaleout -> file where rescaling coefficients are stored
    outvar -> where to store rescaled variance 
    cut -> percentile to use in line rejection
    smooth -> s parameter in spline evaluation
    block -> width of wave window for outlier rejections
    disp -> define the dispersion in block above which rejection is applied
    sthre -> absolute threshold below which data are rejected 
    bootstrap -> if set to a bootstrap variance cube perform rescaling using bootstrap resampling
                 rather than imposing rms
    expmap -> image file with number of exposures per pixel
    expmap_range -> array with min/max values of the exposure map to be used. Range is capped to [1,nmax].                
    memmap -> If false read the cube in memory, faster for numerical operations but more memory intensive

    savechecks-> if set to a root name, save check plots in files
    
    """
 
    #open the data
    data=fits.open(cube, memmap=memmap)
    
    if(bootstrap):
        bscube=fits.open(bootstrap, memmap=memmap)
        
    # #check if AO data (not really needed but that's ok)
    # mode=(data[0].header['ESO INS MODE']).strip()
    # if(('WFM-AO-N' in mode)|('WFM-AO-E' in mode)):
    #     aoflag=True
    # else:
    #     aoflag=False
           
    #compress into image
    image=np.nanmedian(data[1].data,axis=0)
    nx,ny=image.shape
  
    #mask edges
    edges=np.isfinite(image)
    badmask=np.zeros((nx,ny))+1
    badmask[edges]=0.0
    badmask=ndimage.gaussian_filter(badmask,1.5)
    badmask[np.where(badmask > 0)]=1.0

        
    #mask sources
    bkg = sep.Background(image,mask=badmask)    
    thresh = 2.0 * bkg.globalrms
    segmap = np.zeros((nx,ny))
    objects,segmap=sep.extract(image,thresh,segmentation_map=True,
                               minarea=10,clean=True,mask=badmask)
    badmask[np.where(segmap > 0)]=1.0
   
    #If expmap is supplied and range is appropriate take only pixels in expmap range
    if (expmap):
       exphdu = fits.open(expmap)
       expdata = exphdu[0].data
       
       nexp_max = np.nanmax(expdata)
       expmap_range = np.clip(expmap_range, 1, nexp_max)
       
       badmask[(expdata < expmap_range[0]) | (expdata > expmap_range[1])] = 1.0       
       
    goodpix=np.where(badmask < 1)
    
    print('Done with source masking!')
    
    #hdu = fits.PrimaryHDU(badmask)
    #hdu.writeto('temp.fits', overwrite=True)
    
    #check shape of distribution
    #ratio=bscube[0].data[10:50,:,:]/data[2].data[10:50,:,:]
    #ww,xx,yy=ratio.shape
    #allrat=[]
    #for ii in range(ww):
    #    allrat.append(ratio[ii][goodpix].flatten())
    #allrat=np.array(allrat)
    #allrat=allrat[np.where(np.isfinite(allrat))]
    #plt.hist(allrat,range=(-1,25),bins=100)
    #plt.axvline(np.mean(allrat))
    #plt.axvline(np.median(allrat),ls=":")
    #plt.show()
    
    #loop over wave
    nw,nx,ny=data[1].data.shape 
    rescale=[]
    varspec=[]
    wave=[]

    for ww in range(nw):

        #get slices
        slicecube=data[1].data[ww]
        slicevar=data[2].data[ww]
        if(bootstrap):
            sliceboot=bscube[0].data[ww]

        #find good values omitting all nan slices
        if(len((np.where(np.isfinite(slicecube)))[0]) > 0):
                
            #utilities
            varspec.append(np.nanmedian(slicevar[goodpix]))
            wave.append(ww)
            
            #compute scaling factor
            if(bootstrap):
                #rescale.append(np.nanmedian(np.sqrt(sliceboot[goodpix]/slicevar[goodpix])))
                rescale.append(np.nanmean(np.sqrt(sliceboot[goodpix]/slicevar[goodpix])))
                #rescale.append(np.nanmean(np.sqrt(sliceboot[goodpix]))/np.nanmean(np.sqrt(slicevar[goodpix])))
            else:
                #use rms value
                normslice=slicecube[goodpix]/np.sqrt(slicevar[goodpix])
                stddata=np.nanstd(normslice)
                rescale.append(stddata)
        
        #checks
        #plt.hist(normslice[np.logical_not(np.isnan(normslice))],bins=100)
        #plt.title(stddata)
        #plt.show()
        #normslice=normslice/stddata
        #stddata2=np.std(normslice)
        #plt.hist(normslice[np.logical_not(np.isnan(normslice))],bins=100)
        #plt.title(stddata2)
        #plt.show()

    
    rescale=np.array(rescale)
    varspec=np.array(varspec)
    wave=np.array(wave)
    
    print('Computed scaling for individual slices!')

    #plt.plot(wave,rescale)
    #plt.plot(wave,varspec)
    #plt.show()
    

    #now do blocks of wave
    endw=block
    starw=0
    
    selectr=np.array([])
    selectw=np.array([])
    
    while endw < max(wave):
        #trigger selection where scatter is big
        disper=np.std(rescale[starw:endw])
        if(disper > disp):
            #handle AO gap
            try:
                keep=np.where(varspec[starw:endw] < np.percentile(varspec[starw:endw],cut))
                selectr=np.append(selectr,rescale[starw:endw][keep])
                selectw=np.append(selectw,wave[starw:endw][keep])
            except:
                pass
        else:
            selectr=np.append(selectr,rescale[starw:endw][:])
            selectw=np.append(selectw,wave[starw:endw][:])

        #move to next block
        starw=endw+1
        endw=starw+block


    #add last chunk
    startw=endw+1
    endw=max(wave)
    #handle AO gap 
    try:
        keep=np.where(varspec[starw:endw] < np.percentile(varspec[starw:endw],cut))
        selectr=np.append(selectr,rescale[starw:endw][keep])
        selectw=np.append(selectw,wave[starw:endw][keep])
    except:
        pass

    print('Flagged lines!')
    

    #filetred version
    pnt=inter.splrep(selectw,selectr,s=smooth)    
    filterscale=inter.splev(wave,pnt,der=0)

    print('Constructed interpolation!')
    
    #do one more rejection of significant outliers
    bestv=inter.interp1d(wave,filterscale)(selectw)
    dist=abs(bestv-selectr)/bestv
    
    #plt.figure()
    #plt.scatter(selectw,bestv,label='Best')
    #plt.scatter(selectw,selectr,label='Data')
    #plt.scatter(selectw,dist,label='Dist')
    #plt.scatter(selectw[np.where(dist < 0.02)],selectr[np.where(dist < 0.02)],label='Sel')
    #plt.scatter(selectw[np.where(dist < 0.02)],dist[np.where(dist < 0.02)],label='Dist Sel')
    #plt.legend()
    #plt.show()
    
    selectw=selectw[np.where((dist < 0.02) & (selectr > sthre))]
    selectr=selectr[np.where((dist < 0.02) & (selectr > sthre))]
    
    #selectw=selectw[np.where((dist < 0.02) & (selectr > 0.88))]
    #selectr=selectr[np.where((dist < 0.02) & (selectr > 0.88))]
    
    print('Performed final rejection!')
     
    #re-filetred version trimming edges
    pnt=inter.splrep(selectw[2:-2],selectr[2:-2],s=smooth)    
    filterscale=inter.splev(np.arange(nw),pnt,der=0)

    print('Constructed final interpolation function!')

    #some plots to show what's going on
    plt.figure()
    plt.ylabel('Std rescale')
    plt.xlabel('Wave')
    plt.scatter(wave,rescale,label='All')    
    plt.scatter(selectw[2:-2],selectr[2:-2],label='Used')
    plt.plot(np.arange(nw),filterscale,color='black',label='Final')
    plt.legend()
    if(savechecks):
        plt.ylim([sthre*0.9,2.0])
        plt.savefig(savechecks+"_scaling.pdf")
    else:
        plt.show()

    #now rescale variance
    newvar=data[2].data
    newrms=np.array([])

    txtout=open(rescaleout,"w+")
    txtout.write("#SliceNumber  VarRescale\n")
    
    for ww in range(nw):
        #First thing write scale txt file
        txtout.write("{} {}\n".format(ww+1,filterscale[ww]**2))
        
        #get slices - handling AO gap, MFossati, unclear why we do need try statement
        #try:
        slicecube=data[1].data[ww]
        slicevar=data[2].data[ww]
        newvar[ww]=newvar[ww]*filterscale[ww]**2
        
        if(bootstrap):
            sliceboot=bscube[0].data[ww]
            newrms=np.append(newrms,np.nanmean(np.sqrt(sliceboot[goodpix]/newvar[ww][goodpix])))
            #newrms=np.append(newrms,np.nanmedian(np.sqrt(sliceboot[goodpix]/newvar[ww][goodpix])))
            #newrms=np.append(newrms,np.nanmean(np.sqrt(sliceboot[goodpix]))/np.nanmean(np.sqrt(newvar[ww][goodpix])))
        else:
            pix=slicecube[goodpix]/np.sqrt(newvar[ww][goodpix])
            newrms=np.append(newrms,np.nanstd(pix))
                
        #except:
        #    pass

    plt.figure()
    plt.ylabel('Rescaled rms')
    plt.xlabel('Wave')
    plt.scatter(np.arange(nw),newrms)    
    plt.legend()
    if(savechecks):
        plt.ylim([0.7,1.2])
        plt.savefig(savechecks+"_rescaled.pdf")
    else:
        plt.show()

    print('Rescaled arrays!')

    #save fits
    data[2].data=newvar
    data.writeto(outvar,overwrite=True)

    print('All data saved!')

    data.close()



