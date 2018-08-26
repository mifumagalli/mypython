import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits 
import sep
from scipy import ndimage
from scipy import signal
from scipy import interpolate as inter
import multiprocessing as mp
import os

def evaluatenoise(ww,dictin):
    
    """
    Utility function that evaluates boostrap noise
    
    """
    
    print(ww,mp.current_process())
    
    
    for i in range(10000):
        a=0

    thisslice=np.zeros((dictin['nx'],dictin['ny']))

    

    return thisslice


def bootstrapnoise(cubes,masks=None,nsamp=10000,outvar="bootstrap_variance.fits",nproc=10):

    """
    
    Take a list of exposures and estimate variance with boostrap
    at pixel level

    cubes -> list of input cubes
    masks -> list of associated masks
    nsamp -> number of samples to draw
    outvar -> where to store output
    nproc -> number of proc to run this over 

    
    """

    #load the exposures and stash them in pointer stack
    allexposures=[]
    allmasks=[]
    nexp=0
    for exp in open(cubes):
        allexposures.append(fits.open(exp.strip()))
        nexp=nexp+1
    print('Opened {} exposures'.format(nexp))

    if(masks):
        for exp in open(masks):
            allmasks.append(fits.open(exp.strip()))
    print('Opened {} masks'.format(nexp))

        
    #find format of data and create empty var
    nw,nx,ny=allexposures[0][1].data.shape
    print('Data format {} {} {}'.format(nw,nx,ny))
    newvar=np.zeros((nw,nx,ny))

    
    #loop over pixels
    for ww in range(nw):
        print('Working on slice {}/{}'.format(ww+1,nw))
        for xx in range(nx):
            for yy in range(ny):
                #initialise arrays
                fluxstack=np.array([])
                npix=0
                #ingest pixel
                for exp in range(nexp):
                    if(masks):
                        maskpix=allmasks[exp][0].data[xx,yy]
                    else:
                        maskpix=1
                    fluxpix=allexposures[exp][1].data[ww,xx,yy]
                    if((maskpix > 0) & (np.isfinite(fluxpix))):
                        fluxstack=np.append(fluxstack,fluxpix)
                        npix=npix+1
                #bootstrap
                if(npix > 0):
                    #bootstrap
                    rindex=np.random.randint(npix,size=nsamp)
                    newvar[ww,xx,yy]=np.std(fluxstack[rindex])**2
                else:
                    newvar[ww,xx,yy]=np.nan
                    
    #save fits
    hdu=fits.PrimaryHDU(newvar)
    hdu.writeto(outvar,overwrite=True)

    
    #nproc=4
    #nw=50
    #pool=mp.Pool(processes=nproc)
    #argsdic={'nx':nx,'ny':ny,'cubes':cubes,'mask':mask}    
    #varslices=[pool.apply(evaluatenoise,args=(ww,argsdic,)) for ww in range(nw)]

        


def rescalenoise(cube,rescaleout="rescale_variance.txt",outvar="rms_rescaled_var.fits",cut=10,smooth=1,block=65,disp=0.07):
    
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

    """
 
    #open the data
    data=fits.open(cube)

    #compress into image
    image=np.median(data[1].data,axis=0)
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
   
    #loop over wave
    nw,nx,ny=data[1].data.shape 
    rescale=[]
    varspec=[]
    wave=[]

    for ww in range(nw):
        #get slices
        slicecube=data[1].data[ww]
        slicevar=data[2].data[ww]
        normslice=slicecube[goodpix]/np.sqrt(slicevar[goodpix])

        #utilities
        varspec.append(np.median(slicevar[goodpix]))
        wave.append(ww)

        #compute scaling factor
        stddata=np.std(normslice)
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
    
    #now do blocks of wave
    endw=block
    starw=0
    
    selectr=np.array([])
    selectw=np.array([])
    
    while endw < max(wave):
        #trigger selection where scatter is big
        disper=np.std(rescale[starw:endw])
        if(disper > disp):
            keep=np.where(varspec[starw:endw] < np.percentile(varspec[starw:endw],cut))
            selectr=np.append(selectr,rescale[starw:endw][keep])
            selectw=np.append(selectw,wave[starw:endw][keep])
        else:
            selectr=np.append(selectr,rescale[starw:endw][:])
            selectw=np.append(selectw,wave[starw:endw][:])

        #move to next block
        starw=endw+1
        endw=starw+block


    #add last chunk
    startw=endw+1
    endw=max(wave)
    keep=np.where(varspec[starw:endw] < np.percentile(varspec[starw:endw],cut))
    selectr=np.append(selectr,rescale[starw:endw][keep])
    selectw=np.append(selectw,wave[starw:endw][keep])
            
    #filetred version
    pnt=inter.splrep(selectw,selectr,s=smooth)    
    filterscale=inter.splev(wave,pnt,der=0)

    #do one more rejection of significant outliers
    bestv=inter.interp1d(wave,filterscale)(selectw)
    dist=abs(bestv-selectr)/bestv
    
    selectr=selectr[np.where(dist < 0.02)]
    selectw=selectw[np.where(dist < 0.02)]
         
    #re-filetred version trimming edges
    pnt=inter.splrep(selectw[2:-2],selectr[2:-2],s=smooth)    
    filterscale=inter.splev(wave,pnt,der=0)

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
        #get slices
        slicecube=data[1].data[ww]
        slicevar=data[2].data[ww]
        newvar[ww]=newvar[ww]*filterscale[ww]**2
        pix=slicecube[goodpix]/np.sqrt(newvar[ww][goodpix])
        newrms=np.append(newrms,np.std(pix))
        txtout.write("{} {}\n".format(ww+1,filterscale[ww]**2))

    plt.figure()
    plt.ylabel('Rescaled rms')
    plt.xlabel('Wave')
    plt.scatter(wave,newrms)    
    plt.legend()

    #save fits
    hdu=fits.PrimaryHDU(newvar)
    hdu.writeto(outvar,overwrite=True)

    data.close()
    plt.show()

    





