def aligntocat(cube,catalogue):

    """ 

    Take a cube and compare the sources in the field to the positions in the reference 
    catalogue. Return offsets in ra/dec that can be used for re-projecting the cube

    This makes use of pyguide to find stars and centroid on them

    """
       
    import PyGuide
    import numpy as np
    from astropy import wcs

    print "Processing cube {} for offsets".format(cube)

    #first project the cube to create a white image with wcs
    img, var, wcsimg = cube2img(cube)

    #load the reference stars
    starcat=np.loadtxt(catalogue,dtype={'names': ('ra','dec'),'formats': ('f4', 'f4')})

    #loop on reference stars and centroid 
    nref=len(starcat['ra'])

    #space for offsets
    raoffcurrent=[]
    decoffcurrent=[]
       
    for sr in range(nref):
        #find first guess
        pix = wcsimg.wcs_world2pix(starcat['ra'][sr],starcat['dec'][sr],0)
        #centroid within 20 pixels 
        ccd=PyGuide.CCDInfo(0,0,1)
        
        #handle weird errors for bas pixels
        usethis=True
        try:
            centroid=PyGuide.centroid(img,None,None,[pix[0],pix[1]],20,ccd)
        except:        
            usethis=False

        #back to ra/dec if centroid is good
        if ((centroid.isOK) & (usethis)): 
            coord = wcsimg.wcs_pix2world(centroid.xyCtr[0],centroid.xyCtr[1],0)
            raoffcurrent.append(coord[0]-starcat['ra'][sr])
            decoffcurrent.append(coord[1]-starcat['dec'][sr])
                
    #stack offsets
    raoff=np.median(np.array(raoffcurrent))
    decoff=np.median(np.array(decoffcurrent))
     
    print "Offsets for cube {} RA: {} Dec: {}".format(cube,raoff*3600.,decoff*3600.)
    print "Error offsets for cube {} RA: {} Dec: {}".format(cube,np.std(np.array(raoffcurrent))*3600.,
                                                            np.std(np.array(decoffcurrent))*3600.)
                             
    return raoff,decoff

def cube2img(cube,write=None,wrange=None,helio=0,filt=None):
    
    """
    Take a cube and make a projection, extracting also WCS information.
    The projection is constructued with a mean, and the var is 
    consistently propagated.  
    
    wrange -> if set to a (minl,maxl) compress only part of the cube
    write -> if set to string, write image in output
    helio -> passes heliocentric correction in km/s [not needed for pipeline v1.2.1 and later]
    filt -> The ID of a filter transmission curve to use for convolution,
            which filter package can understand (e.g. 129 for SDSS r-band) 
            Filt overwrite wrange

    """

    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits 
    from astropy import wcs
    from mypython.filters import filter as fil
    from scipy import interpolate
    from scipy import integrate

    #read the cube
    cubdata,vardata,wcsc,wavec,regions=readcube(cube,helio=helio)
        
    #find delta lambda
    delta_lambda=wavec-np.roll(wavec,1)
    delta_lambda[0]=wavec[1]-wavec[0]

    #########################################
    #compute the desired transmission curves#
    #########################################

    #if no filter or wrange, default to mix max cube
    if not (wrange) and not (filt):
        wrange=[np.min(wavec),np.max(wavec)]

    #implement filtering in wrange
    if(wrange):
        
        #create a good mask escluding NaN
        mask=np.isfinite(cubdata)
        mask=mask.astype(int)
        
        #trim data in desired range
        trim=np.where((wavec < wrange[0]) | (wavec > wrange[1]))
        mask[trim[0],:,:]=0
        
        #set the trasmission T=1 for white/tophat image masking 5 pixel at edges
        trans=np.zeros(len(wavec))+1
        trans[0:5]=0
        trans[-5:]=0

        #set zeropoint 
        #now compute the zero point of the image
        lefffn=interpolate.interp1d(wavec,wavec*trans,bounds_error=False,fill_value=0)
        trans_fnc=interpolate.interp1d(wavec,trans,bounds_error=False,fill_value=0)
        num,err=integrate.quad(lefffn,wrange[0],wrange[1])
        den,err=integrate.quad(trans_fnc,wrange[0],wrange[1])
        lmean=num/den
        print 'Filter mean wavelength ', lmean 
        #Compute the ZP in AB - native image is in 1e-20 erg/s/cm2/A
        ZP=-2.5*np.log10(lmean*lmean/29979245800.*1e-8*1e-20)-48.6
        print 'Filter zeropoint ',  ZP
        
    #implement filter transmission
    if(filt):
        #create a good mask - this overwrites wrange
        mask=np.isfinite(cubdata)
        mask=mask.astype(int)
        
        #load the filter 
        myfilter=fil.Filter(filt)
        myfilter.loadtrans()

        #restrict to sensible range in transmission 
        good=np.where(myfilter.filter['tran'] > 5e-4)
        good=good[0]
        fltn=myfilter.filter['tran'][good]
        flwn=myfilter.filter['wave'][good]

        #interpolate transmission on cube
        trans_fnc=interpolate.interp1d(flwn,fltn,bounds_error=False,fill_value=0)
        trans=trans_fnc(wavec)
        
        #plt.plot(flwn,fltn)
        #plt.scatter(wavec,trans)
        #plt.show()
        #exit()

        #now compute the zero point of the image
        lefffn=interpolate.interp1d(flwn,flwn*fltn)
        num,err=integrate.quad(lefffn,min(flwn),max(flwn))
        den,err=integrate.quad(trans_fnc,min(flwn),max(flwn))
        lmean=num/den
        print 'Filter mean wavelength ', lmean 
        #Compute the ZP AB - image is in 1e-20 erg/s/cm2/A
        ZP=-2.5*np.log10(lmean*lmean/29979245800.*1e-8*1e-20)-48.6
        print 'Filter zeropoint ',  ZP


    ###############################
    #now do the actual combination#
    ###############################

    #fix nans [already masked above]
    cubdata=np.nan_to_num(cubdata)
    vardata=np.nan_to_num(vardata)

    #combine trans with delta lambda
    trans=trans*np.tile(delta_lambda,((cubdata.shape)[2],(cubdata.shape)[1],1))
    trans=np.transpose(trans)

    #do not vectorise - as this requires too much memory
    nx=cubdata.shape[1]
    ny=cubdata.shape[2]
    img=np.zeros((nx,ny))
    var=np.zeros((nx,ny))
    wgt=np.zeros((nx,ny))

    #Do for loop for better memory usage
    for xx in range(nx):
        for yy in range(ny):
            img[xx,yy]=np.sum(cubdata[:,xx,yy]*trans[:,xx,yy]*mask[:,xx,yy])
            var[xx,yy]=np.sum(vardata[:,xx,yy]*mask[:,xx,yy]*trans[:,xx,yy]**2)
            wgt[xx,yy]=np.sum(mask[:,xx,yy]*trans[:,xx,yy])
            if(wgt[xx,yy] <= 0):
                wgt[xx,yy]=1e-20

    #reassemble - this is constrcuting a mean image weighted by the transmission curve
    img=np.nan_to_num(img/wgt)
    var=np.nan_to_num(var/wgt**2)

    #grab 2D wcs
    wcsimg=wcsc.dropaxis(2)
    
    #if write, write
    if(write):
        print 'Writing to ', write
        header=wcsimg.to_header()
        header["ZPAB"]=ZP
        hduhead = fits.PrimaryHDU(img,header=header)
        hduimg  = fits.ImageHDU(img)
        hduvar  = fits.ImageHDU(var)
        hdulist = fits.HDUList([hduhead,hduimg,hduvar])
        hdulist.writeto(write,clobber=True)

    return img, var, wcsimg

def cube2spec(cube,x,y,s,write=None,shape='box',helio=0,mask=None,twod=True,tovac=False,idsource=None):

    """ 
    Extract a 1D spectrum from a cube at position x,y in box or circle of radius s 

    If shape = 'mask', then mask is a boolean mask and pixels within it will be extracted form 
    argument mask. Mask is a datacube [e.g. from cubex]

    idsource -> if > 0, then only pixels in mask with that ID will be extracted

    helio passes an heliocentric correction in km/s [should be 0 with pipeline v1.2.1]

    twod -> also reconstruct a 2D spec

    tovac -> if true, return wavelengths in vacuum 

    write -> output file 

    """
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.io import fits 

    #read the cube
    cubdata,vardata,wcsc,wavec,regi=readcube(cube,helio=helio)
    cubdata=np.nan_to_num(cubdata)

    #if mask extract all True pixels 
    if('mask' in shape):
        if(idsource):
            goodpix=np.nonzero(mask == idsource)
        else:
            goodpix=np.nonzero(mask)
        xpix=goodpix[1]
        ypix=goodpix[2]
    else:
        #If user defined region, grab inner pixels
        #cut region of interest according to shape
        xpix=[]
        ypix=[]
        xside=np.arange(x-s-1,x+s+1,1)
        yside=np.arange(y-s-1,y+s+1,1)
        for xx in xside:
            for yy in yside:
                if('box' in shape):
                    if((abs(xx-x) <= s) & (abs(yy-y) <= s)):
                        xpix.append(xx)
                        ypix.append(yy)
                if('circ' in shape):
                    dist=np.sqrt((xx-x)**2+(yy-y)**2)
                    if(dist <= s):
                        xpix.append(xx)
                        ypix.append(yy)
                        
    #Some checks...
    #cbmed=np.median(cubdata,axis=0)
    #cbmed[xpix,ypix]=100000
    #imgplot=plt.imshow(cbmed,origin='lower')
    #imgplot.set_clim(-5,5)
    #plt.show()

    #now make space for the 1d spectrum 
    spec_flx=np.zeros(len(wavec))
    spec_var=np.zeros(len(wavec))
    spec_med=np.zeros(len(wavec))

    #if want 2d, prepapre space for it
    #This simulates a slit in the x direction 
    #adding up all the flux on the y
    if(twod):
        #find unique pixels (not all x,y) - need to sort for later
        uxpix=np.sort(list(set(xpix)))
        uypix=np.sort(list(set(ypix)))
        npix=len(uxpix)
        nwv=len(wavec)
        twodspec=np.zeros((nwv,npix))
        twoderr=np.zeros((nwv,npix))

    #loop over all wavelength to fill in spectrum
    for ii,ww in enumerate(wavec):
        #get the total spec in the aperture, 
        #summing over all the pixels 
        spec_flx[ii]=np.sum(cubdata[ii,xpix,ypix])
        spec_var[ii]=np.sum(vardata[ii,xpix,ypix])
        spec_med[ii]=np.median(cubdata[ii,xpix,ypix])
        
        #fill in 2D spectrum in a box 
        if(twod):
            #sum only on the full x-extent
            for jj in range(npix):
                #add all the pixels in y
                twodspec[ii,jj]=np.sum(cubdata[ii,uxpix[jj],uypix])
                twoderr[ii,jj]=np.sum(vardata[ii,uxpix[jj],uypix])
   
    #extract the 2D image with a small buffer around
    if(twod):
        twodimg=np.median(cubdata[:,uxpix[0]-5:uxpix[-1]+6,uypix[0]-5:uypix[-1]+6],axis=0)
        #from variance to error
        twoderr=np.sqrt(twoderr)

    #mean in aperture
    #totpix=len(xpix)
    #spec_flx=spec_flx/totpix
    #spec_err=np.sqrt(spec_var/totpix)

    #keep total spectrum and not mean
    totpix=len(xpix)
    spec_err=np.sqrt(spec_var)
    
    #if set, convert to vacuum using airtovac.pro conversion
    if(tovac):
        #save current wave
        wavec=np.array(wavec,dtype=np.float64)
        wave_air=wavec
        
        sigma2 = (1e4/wavec)**2.    
        fact = 1.+5.792105e-2/(238.0185-sigma2)+1.67917e-3/(57.362-sigma2)
        wavec = wavec*fact  

    #tested and working
    #fl=open('test.txt','w') 
    #for rr in range(len(wavec)):
    #    fl.write("{} {}\n".format(wave_air[rr],wavec[rr]))
    #fl.close()

    #if write, write
    if(write):
        prihdr = fits.Header()
        prihdr['NPIX'] = totpix
        hduflx  = fits.PrimaryHDU(spec_flx,header=prihdr) #total in region
        hduerr  = fits.ImageHDU(spec_err) #associated errors
        hduwav  = fits.ImageHDU(wavec)    #wave
        hdumed  = fits.ImageHDU(spec_med) #median spectrum 
        if(twod): #twod 
            hdu2flx  = fits.ImageHDU(twodspec)
            hdu2err  = fits.ImageHDU(twoderr)
            hduimg   = fits.ImageHDU(twodimg)
            hdulist = fits.HDUList([hduflx,hduerr,hduwav,hdumed,hdu2flx,hdu2err,hduimg])
        else:
            hdulist = fits.HDUList([hduflx,hduerr,hduwav,hdumed])
        hdulist.writeto(write,clobber=True)

    return wavec, spec_flx, spec_err, spec_med

def cubestat(cube,region=None,delta=10,mask=None):

    """
    Take the cube and measure the pixel rms in chunks of 10A
    within a specified region 

    region -> False, use the entire cube
              or set to min x,y max x,y or region to be used
    mask -> if set to a segmentation map, exclude soruces and 
            gap 

    delta  -> wavelength window in A

    """

    import numpy as np
    import matplotlib.pyplot as plt
    
    #read the cube
    cubdata,vardata,wcsc,wavec,reg=readcube(cube)
    
    #grab info on grid
    lambdabin=wcsc.pixel_scale_matrix[2,2]*1e10 #in A
    pixbin=wcsc.pixel_scale_matrix[1,1]*3600 #in as

    #find blocks of lambda within the cube 
    nblocks=int(np.floor((np.max(wavec)-np.min(wavec))/delta))
    
    #make empty mask if not provided
    if(mask == None):
        nz,nx,ny=cubdata.shape
        mask=np.zeros((nx,ny))
        
    #carve out box if needed in spatial direction
    if(region):
        cubdata=cubdata[:,region[0]:region[2],region[1]:region[3]]
        mask=mask[region[0]:region[2],region[1]:region[3]]
            
    #init arrays
    rms=np.zeros(nblocks)
    wrms=np.zeros(nblocks)

    #loop over the blocks 
    for ii in range(nblocks):

        #find pixels in this block
        wmin=np.min(wavec)+ii*delta
        wmax=wmin+delta
        wcent=wmin+0.5*delta
        wpix,=np.where((wavec >= wmin) & (wavec < wmax))
        allpix=[]
        #make stack of pixels
        for ww in wpix:
            pixgood=np.where(mask < 1)
            layer=cubdata[ww,:,:]
            allpix.append(layer[pixgood])
            
        rms[ii]=np.std(np.nan_to_num(allpix))
        wrms[ii]=wcent
        
    #normalise units from pixel to as^2 and from pix to A
    rms=rms*1e-20/lambdabin/pixbin**2
    return wrms,rms

def readcube(cube, helio=0):

    """
    Read a cube, expanding wcs and wavelegth

    If setting helio!=0 (in km/s), then helio corrections 
    are applied to the data 

    Note that from v1.2, the ESO pipeline applies correction to 
    data regarding baryocentric heliocentric correction, so helio should be left a 0
  
    MUSE wavelength solution is in air!!

    """

    from astropy.io import fits
    import astropy as aspy
    from astropy.wcs import WCS
    import numpy as np

    #open file 
    cfits=fits.open(cube)
    
    #grab the data
    cubdata=cfits['DATA'].data
    vardata=cfits['STAT'].data

    #compute the helio correction on the fly
    if(helio != 0):
        hel_corr = np.sqrt( (1. + helio/299792.458) / (1. - helio/299792.458) )
        print 'Helio centric correction of {} km/s and lambda {}'.format(helio,hel_corr) 
    else:
        hel_corr=1.0

    #reconstruct wave array
    #wave in air
    sz=cubdata.shape
    delta_lambda=cfits[1].header["CD3_3"]*hel_corr
    zero_lambda=cfits[1].header["CRVAL3"]*hel_corr
    wavec=np.arange(0,sz[0],1)*delta_lambda+zero_lambda
 
    #compute ds9 style regions
    regions=np.arange(0,sz[0],1)+1

    #suck up the wcs 
    wcsc=WCS(header=cfits[1].header)
    #wcscube.printwcs()
    
    #close unit
    cfits.close()

    return cubdata,vardata,wcsc,wavec,regions


def adjust_wcsoffset(data,xpix,ypix,rag,deg):


    """

    Small script that corrects the overall zero-point of the 
    wcs solution. 

    data -> cube or image to which correction should be applied
    xpix -> new reference x pixel, or CRPIX1 [best if close to centre of field]
    ypix -> same but for y, or CRPIX2
    rag  -> ra in deg for reference pixel, or CRVAL1
    deg  -> dec in deg for reference pixel, or CRVAL2 

    """

    from astropy.io import fits

    #open 
    fithdu=fits.open(data,mode='update')
    
    try:

        #save old 
        fithdu[1].header['OLDCRV1']=fithdu[1].header['CRVAL1']
        fithdu[1].header['OLDCRV2']=fithdu[1].header['CRVAL2'] 
        fithdu[1].header['OLDCPX1']=fithdu[1].header['CRPIX1'] 
        fithdu[1].header['OLDCPX2']=fithdu[1].header['CRPIX2'] 

        #write new 
        fithdu[1].header['CRVAL1']=rag 
        fithdu[1].header['CRVAL2']=deg
        fithdu[1].header['CRPIX1']=xpix 
        fithdu[1].header['CRPIX2']=ypix
        
    except:

        #save old
        fithdu[0].header['OLDCRV1']=fithdu[0].header['CRVAL1'] 
        fithdu[0].header['OLDCRV2']=fithdu[0].header['CRVAL2'] 
        fithdu[0].header['OLDCPX1']=fithdu[0].header['CRPIX1'] 
        fithdu[0].header['OLDCPX2']=fithdu[0].header['CRPIX2'] 
        
        #write new 
        fithdu[0].header['CRVAL1']=rag 
        fithdu[0].header['CRVAL2']=deg
        fithdu[0].header['CRPIX1']=xpix 
        fithdu[0].header['CRPIX2']=ypix

    #save 
    fithdu.flush()
    fithdu.close()


def unpack_pixtab(flag):
    
    """

    Take the flag extension of a pixel table and return unpacked information 
    on the origin of the pixels 

    flag - > the pixel origin extension of a pixel table

    """

    xpix = flag >> 24
    ypix = (flag >> 11) & 8191
    ifu  = (flag >> 6) & 31
    islice = flag & 63 

    return ifu, islice


def check_flux_scaling(reference,listexp,maskedges=None,verbose=True,flxlim=150.):

    """

    Perform basic photomety checks to make sure that 
    none of the exposures in a coadd has a very weird calibration error
    relative to the median

    reference -> the refence cube to use as gold standard
    listexp   -> a text file containing the file names of the cubes to check
    maskedges -> if set to a number clips Npixels from the edge
    verbose -> make a bunch of check plots
    flxlim  -> how bright sources to be considered 

    """
    
    from astropy.io import fits 
    import numpy as np 
    import matplotlib.pyplot as plt
    import sep

    #open reference and create image
    refcube=fits.open(reference)
    refimage=np.nanmedian(refcube[1].data,axis=0)
    header=refcube[1].header

    #run source id on reference image - pick brightest only
    bkg = sep.Background(refimage)    
    thresh = 5. * bkg.globalrms
    minarea=20.
    clean=True
    segmap = np.zeros((header["NAXIS2"],header["NAXIS1"]))
    mask=np.zeros((header["NAXIS2"],header["NAXIS1"]))
    
    #create mask if desired 
    if(maskedges):
        mask[0:maskedges,:]=1
        mask[-maskedges:-1,:]=1
        mask[:,0:maskedges]=1
        mask[:,-maskedges:-1]=1
    
    if(verbose):
        plt.imshow(mask)
        plt.show()
        
    
    #extract objects
    objects,segmap=sep.extract(refimage,thresh,segmentation_map=True,
                               minarea=minarea,clean=clean,mask=mask)
    if(verbose):
        plt.imshow(segmap)
        plt.show()
    
    #grab phot on ref imge
    ref_phot=[]
    for ii in range(len(objects)):
        pix=np.where(segmap == ii+1)
        ref_phot.append(np.nansum(refimage[pix]))
    ref_phot=np.array(ref_phot)

    use=np.where(ref_phot > flxlim)

    #get photometry for check images
    fl=open(listexp)
    for ff in fl:
        #create image 
        cubck=fits.open(ff.strip())
        imck=np.nanmedian(cubck[1].data,axis=0)
        #get photometry
        chk_phot=[]
        for ii in range(len(objects)):
            pix=np.where(segmap == ii+1)
            chk_phot.append(np.nansum(imck[pix]))
        chk_phot=np.array(chk_phot)
        if(verbose):
            plt.scatter(ref_phot[use],chk_phot[use]/ref_phot[use])
            plt.show()
        print ("Scaling for {} is {}".format(ff,np.median(chk_phot[use]/ref_phot[use])))
        
    fl.close()

   

    


    
