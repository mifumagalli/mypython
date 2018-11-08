import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits 
import sep
from scipy import ndimage
from scipy import signal
from scipy import interpolate as inter
import multiprocessing as mp
import os
import glob
import datetime

def evaluatenoise(iproc,wstart,wend,nx,ny,nexp,nsamp,allexposures,allmasks,masks,median):
    
    """
    Utility function that evaluates boostrap noise
    
    """
    
    #now load the exposures
    dataexp=np.zeros((nexp,wend-wstart+1,nx,ny))
    datamas=np.zeros((nexp,nx,ny))
    
    iexp=0
    for exp in allexposures:
        dataexp[iexp]=fits.open(exp)[1].data[wstart:wend+1,:,:]
        iexp=iexp+1

    iexp=0
    if(masks):
        for exp in allmasks:
            datamas[iexp]=fits.open(exp)[0].data[:,:]
            iexp=iexp+1

    print('Proc {}: All data loaded'.format(iproc))        

    #make space for output
    newvar=np.zeros((wend-wstart+1,nx,ny))


    #select evaluator
    if(median): 
        from numpy import median as evaluator
    else:
        from numpy import mean as evaluator


    #loop over slices (include tail)
    for ww in range(wstart,wend+1):
        print('Proc {}: Working on slice {}/{}'.format(iproc,ww,wend))
       
        #giant loop on pixels - not all pixels are the same length 
        for xx in range(nx):
            for yy in range(ny):
                #ingest mask pixel
                if(masks):
                    maskpix=datamas[:,xx,yy]
                else:
                    maskpix=np.zeros(nexp)+1
                    
                #ingest flux 
                fluxpix=dataexp[:,ww-wstart,xx,yy]
                
                #apply mask
                fluxpix=fluxpix[np.where((maskpix > 0) & (np.isfinite(fluxpix)))]
                npix=len(fluxpix)

                #bootstrap
                if(npix > 1):
                    #bootstrap
                    rindex=np.random.randint(npix,size=(nsamp,npix))
                    newvar[ww-wstart,xx,yy]=np.std(evaluator(fluxpix[rindex],axis=1))**2
                else:
                    newvar[ww-wstart,xx,yy]=np.nan
                    
    #save output
    np.savez("boostrapvar_tmpout_proc{}".format(iproc),wstart=wstart,wend=wend,newvar=newvar)
    print("Proc {}: Done!".format(iproc))


def bootstrapnoise(cubes,masks=None,nsamp=10000,outvar="bootstrap_variance.fits",nproc=50,median=False):

    """
    
    Take a list of exposures and estimate variance with boostrap
    at pixel level

    cubes -> list of input cubes
    masks -> list of associated masks
    nsamp -> number of samples to draw [ideally 500000]
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
    print('Found {} exposures'.format(nexp))

    if(masks):
        for exp in open(masks):
            allmasks.append(exp.strip())
    print('Found {} masks'.format(nexp))
    
    #find format of data and create empty var
    nw,nx,ny=fits.open(allexposures[0])[1].data.shape
    print('Data format {} {} {}'.format(nw,nx,ny))
    
    #now prepare batches for parallel run
    itempbatch=nw/nproc
    print('Running on {} proc with batches of {}'.format(nproc,itempbatch))

    #loop over processors and start parallel function
    wstart=0
    wend=np.round(itempbatch)
    processes=[]
    for iproc in range(nproc):
        #make sure does not run over index 
        wend=np.minimum(wend,nw-1)
        print('Proc {}: Start slice {} End slice {}'.format(iproc,wstart,wend))
        p=mp.Process(target=evaluatenoise,args=(iproc,wstart,wend,nx,ny,nexp,
                                                nsamp,allexposures,allmasks,masks,median))
        processes.append(p)
        p.start()
        
        #update for next loop
        wstart=wend+1
        wend=wstart+itempbatch
        if(iproc == nproc -1):
            wend=nw-1
            
    for p in processes:
        p.join()
    
    #reconstruct variance array 
    allvar=np.zeros((nw,nx,ny))
    for iproc in range(nproc):
        thisproc=np.load("boostrapvar_tmpout_proc{}.npz".format(iproc))
        allvar[thisproc['wstart']:thisproc['wend']+1,:,:]=thisproc['newvar']

    #save to fits file
    hdu=fits.PrimaryHDU(allvar)
    hdu.writeto(outvar,overwrite=True)    

    #clean tmp files
    alltmp=glob.glob("boostrapvar_tmpout_proc*.npz")
    for tmpfile in alltmp:
        os.remove(tmpfile)

    print('All done at {}'.format(datetime.datetime.now()))


def rescalenoise(cube,rescaleout="rescale_variance.txt",outvar="CUBE_rmsvar.fits",cut=10,smooth=1,block=65,disp=0.07,bootstrap=None):
    
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
    bootstrap -> if set to a bootstrap variance cube perform rescaling using bootstrap resampling
                 rather than imposing rms 

    """
 
    #open the data
    data=fits.open(cube)
    
    if(bootstrap):
        bscube=fits.open(bootstrap)
        
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
   
    goodpix=np.where(badmask < 1)
    
    print('Done with source masking!')

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
        
        ##checks
        #plt.hist(normslice,bins=100)
        #plt.title(stddata)
        #plt.show()
        #normslice=normslice/stddata
        #stddata2=np.std(normslice)
        #plt.hist(normslice,bins=100)
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
    
    selectw=selectw[np.where((dist < 0.02) & (selectr > 1))]
    selectr=selectr[np.where((dist < 0.02) & (selectr > 1))]
    
    print('Performed final rejection!')
     
    #re-filetred version trimming edges
    pnt=inter.splrep(selectw[2:-2],selectr[2:-2],s=smooth)    
    filterscale=inter.splev(wave,pnt,der=0)

    print('Constructed final interpolation function!')

    #some plots to show what's going on
    plt.figure()
    plt.ylabel('Std rescale')
    plt.xlabel('Wave')
    plt.scatter(wave,rescale,label='All')    
    plt.scatter(selectw[2:-2],selectr[2:-2],label='Used')
    plt.plot(wave,filterscale,color='black',label='Final')
    plt.legend()


    #now rescale variance
    newvar=data[2].data
    newrms=np.array([])

    txtout=open(rescaleout,"w+")
    txtout.write("SliceNumber  VarRescale\n")
    
    for ww in range(nw):
        #get slices - handling AO gap
        try:
            slicecube=data[1].data[ww]
            slicevar=data[2].data[ww]
            newvar[ww]=newvar[ww]*filterscale[ww]**2
            txtout.write("{} {}\n".format(ww+1,filterscale[ww]**2))
            
            if(bootstrap):
                sliceboot=bscube[0].data[ww]
                newrms=np.append(newrms,np.nanmean(np.sqrt(sliceboot[goodpix]/newvar[ww][goodpix])))
                #newrms=np.append(newrms,np.nanmedian(np.sqrt(sliceboot[goodpix]/newvar[ww][goodpix])))
                #newrms=np.append(newrms,np.nanmean(np.sqrt(sliceboot[goodpix]))/np.nanmean(np.sqrt(newvar[ww][goodpix])))
            else:
                pix=slicecube[goodpix]/np.sqrt(newvar[ww][goodpix])
                newrms=np.append(newrms,np.nanstd(pix))
                
        except:
            pass

    plt.figure()
    plt.ylabel('Rescaled rms')
    plt.xlabel('Wave')
    plt.scatter(wave,newrms)    
    plt.legend()

    print('Rescaled arrays!')

    #save fits
    data[2].data=newvar
    data.writeto(outvar,overwrite=True)

    print('All data saved!')

    data.close()
    plt.show()



