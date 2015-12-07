

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

def cube2img(cube,write=None,wrange=None,helio=0.0,filt=None):
    
    """
    Take a cube and make a projection, extracting also WCS information.
    The projection is constructued with a mean, and the var is 
    consistently propagated.  
    
    wrange -> if set to a (minl,maxl) compress only part of the cube
    write -> if set to string, write image in output
    helio -> passes heliocentric correction in km/s
    filt -> The ID of a filter transmission curve to use for convolution,
            which filter package can understand (e.g. 129 for SDSS r-band) 
            Filt overwrite wrange

    """


    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits 
    from astropy import wcs
    import mypython.filters as fil
    from scipy import interpolate
    from scipy import integrate

    #read the cube
    cubdata,vardata,wcsc,wavec=readcube(cube,helio=helio)
        
    #find delta lambda
    delta_lambda=wavec-np.roll(wavec,1)
    delta_lambda[0]=wavec[1]-wavec[0]

    #if no filter or wrange, default to mix man cube
    if not (wrange) and not (filt):
        wrange=[np.min(wavec),np.max(wavec)]

    #implement filtering in wrange
    if(wrange):
        #create a good mask
        mask=np.isfinite(cubdata)
        mask=mask.astype(int)
        trim=np.where((wavec < wrange[0]) | (wavec > wrange[1]))
        print 'Trim within range ', wrange
        mask[trim[0],:,:]=0
        
        #set the trasmission T=1 for white image, trim edges
        trans=np.zeros(len(wavec))+1
        trans[0:5]=0
        trans[-5:]=0

        #set zeropoint 
        #now compute the zero point of the image
        lefffn=interpolate.interp1d(wavec,wavec*trans)
        trans_fnc=interpolate.interp1d(wavec,trans,bounds_error=False,fill_value=0)
        num,err=integrate.quad(lefffn,wrange[0],wrange[1])
        den,err=integrate.quad(trans_fnc,wrange[0],wrange[1])
        lmean=num/den
        print 'Filter mean wavelength ', lmean 
        #Compute the ZP - image is in 1e-20 erg/s/cm2/A
        ZP=-2.5*np.log10(lmean*lmean/29979245800.*1e-8*1e-20)-48.6
        print 'Filter zeropoint ',  ZP
        
    #implement filter transmission
    if(filt):
        #create a good mask - this overwrites wrange
        mask=np.isfinite(cubdata)
        mask=mask.astype(int)
        
        #load the filter 
        myfilter=fil.filter.Filter(filt)
        myfilter.loadtrans()

        #restrict to sensible range
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
        #Compute the ZP - image is in 1e-20 erg/s/cm2/A
        ZP=-2.5*np.log10(lmean*lmean/29979245800.*1e-8*1e-20)-48.6
        print 'Filter zeropoint ',  ZP

    #fix nans
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

    #reassembled
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

def cube2spec(cube,x,y,s,write=None,shape='box',helio=0.0,mask=None,twod=True):

    """ 
    Extract a 1D spectrum from a cube at position x,y 
    in box/circle of radius s 

    If shape = 'mask', then mask is a boolean mask and
    pixels within it will be extracted form there

    helio passes an heliocentric correction in km/s

    twod -> reconstruct a 2D spec
    

    """
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.io import fits 

    #read the cube
    cubdata,vardata,wcsc,wavec=readcube(cube,helio=helio)
 
    #if mask extract all True pixels 
    if('mask' in shape):
        goodpix=np.nonzero(mask)
        xpix=goodpix[0]
        ypix=goodpix[1]
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
    spec_opt=np.zeros(len(wavec))
    spec_var=np.zeros(len(wavec))
    spec_varopt=np.zeros(len(wavec))
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
    
    #find the lightweight for optimal extraction 
    weight=np.median(cubdata[:,xpix,ypix],axis=0)
    weight=weight/np.sum(weight)

    for ii,ww in enumerate(wavec):
        #get the total spec in the aperture, 
        #summing over all the pixels 
        spec_flx[ii]=np.sum(cubdata[ii,xpix,ypix])
        spec_opt[ii]=np.sum(cubdata[ii,xpix,ypix]*weight)
        spec_var[ii]=np.sum(vardata[ii,xpix,ypix])
        spec_varopt[ii]=np.sum(vardata[ii,xpix,ypix]*weight**2)
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
    spec_err=np.sqrt(spec_var)
    spec_erropt=np.sqrt(spec_varopt)
    twoderr=np.sqrt(twoderr)

    #plt.plot(wavec, spec_flx)
    #plt.plot(wavec, spec_err)
    #plt.show()

    #if write, write
    if(write):
        hduflxopt  = fits.PrimaryHDU(spec_opt)
        hduerropt  = fits.ImageHDU(spec_erropt)
        hduwav  = fits.ImageHDU(wavec)
        hduflx  = fits.ImageHDU(spec_flx)
        hduerr  = fits.ImageHDU(spec_err)
        if(twod):
            hdu2flx  = fits.ImageHDU(twodspec)
            hdu2err  = fits.ImageHDU(twoderr)
            hduimg   = fits.ImageHDU(twodimg)
            hdulist = fits.HDUList([hduflxopt,hduerropt,hduwav,hduflx,hduerr,hdu2flx,hdu2err,hduimg])
        else:
            hdulist = fits.HDUList([hduflxopt,hduerropt,hduwav,hduflx,hduerr])
        hdulist.writeto(write,clobber=True)

    return wavec, spec_opt, spec_erropt, spec_med



def cubestat(cube,region=False,delta=10,lambdabin=1.25,pixbin=0.2):

    """
    Take the cube and measure the pixel rms in chunks of 10A
    within a specified region 

    region -> False, use the entire cube
              or set to min x,y max x,y or region to be used

    delta  -> wavelength window in A

    lambdabin  -> sampling of wave of final cube  (A)
    pixbin     -> sampling of pixel of final cube (as)

    """

    import numpy as np
    import matplotlib.pyplot as plt

    #read the cube
    cubdata,vardata,wcsc,wavec=readcube(cube)

    #find blocks of lambda
    nblocks=int(np.floor((np.max(wavec)-np.min(wavec))/delta))

    #carve out box if needed
    if(region):
        cubdata=cubdata[:,region[0]:region[2],region[1]:region[3]]

    rms=np.zeros(nblocks)
    wrms=np.zeros(nblocks)

    #loop over the blocks 
    for ii in range(nblocks):

        #find pixels in this block
        wmin=np.min(wavec)+ii*delta
        wmax=wmin+delta
        wcent=wmin+0.5*delta
        wpix=np.where((wavec >= wmin) & (wavec < wmax))

        rms[ii]=np.std(cubdata[wpix,:,:])
        wrms[ii]=wcent

    #normalise units
    rms=rms*1e-20/lambdabin/pixbin**2
   
    return wrms,rms

def readcube(cube, helio=0.0):

    """
    Read a cube, expanding wcs and wavelegth

    If setting helio!=0 (in km/s), then helio corrections 
    are applied to the data 

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
    hel_corr = np.sqrt( (1. + helio/299792.458) / (1. - helio/299792.458) )
    print 'Helio centric correction of {} km/s and lambda {}'.format(helio,hel_corr) 

    #reconstruct wave array
    #wave in air
    sz=cubdata.shape
    delta_lambda=cfits[1].header["CD3_3"]*hel_corr
    zero_lambda=cfits[1].header["CRVAL3"]*hel_corr
    wavec=np.arange(0,sz[0],1)*delta_lambda+zero_lambda
 
    #suck up wcs 
    wcsc=WCS(header=cfits[1].header)
    #wcscube.printwcs()
    
    #close unit
    cfits.close()

    return cubdata,vardata,wcsc,wavec

