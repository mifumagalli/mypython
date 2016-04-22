"""
These are sets of utilities to handle muse sources

"""


def findsources(image,cube,check=False,output='./',spectra=False,helio=0,nsig=2.,
                minarea=10.,regmask=None,clean=True,outspec='Spectra'):

    """      

    Take a detection image (collapse of a cube), or median 
    of an RGB, or whatever you want (but aligned to the cube)
    and run sourceextractor 

   
    Use SEP utilities http://sep.readthedocs.org/en/stable/

    image   -> fits file of image to process
    check   -> if true write a bunch of check mages
    output  -> where to dump the output
    cube    -> the cube used to extract spectra
    spectra -> if True, extract spectra in VACUUM wave!!
    helio   -> pass additional heliocentric correction
    nsig    -> number of skyrms used for source id 
    minarea -> minimum area for extraction 
    regmask -> ds9 region file (image) of regions to be masked before extraction [e.g. edges]
    clean   -> clean souces 
    outspec -> where to store output spectra 

    """

    import sep
    from astropy.io import fits
    import numpy as np
    import os
    from mypython.ifu import muse_utils as utl
    from mypython.fits import pyregmask as msk

    #open image
    img=fits.open(image)
    header=img[0].header
    try:
        #this is ok for narrow band images 
        data=img[1].data
    except:
        #white cubex images
        data=img[0].data
    data=data.byteswap(True).newbyteorder()
    #close fits
    img.close()

    #create bad pixel mask
    if(regmask):
        Mask=msk.PyMask(header["NAXIS1"],header["NAXIS2"],regmask)
        for ii in range(Mask.nreg):
            Mask.fillmask(ii)
            if(ii == 0):
                badmask=Mask.mask
            else:
                badmask+=Mask.mask
            badmask=1.*badmask
    else:
        badmask=np.zeros((header["NAXIS1"],header["NAXIS2"]))

    if(check):
        print('Dumping badmask')
        hdumain  = fits.PrimaryHDU(badmask,header=header)
        hdulist = fits.HDUList([hdumain])
        hdulist.writeto(output+"/badmask.fits",clobber=True)
    
    #check background level, but do not subtract it
    print 'Checking background levels'
    bkg = sep.Background(data,mask=badmask)    
    print 'Residual background level ', bkg.globalback
    print 'Residual background rms ', bkg.globalrms

    if(check):
        print 'Dumping sky...'
        #dump sky properties
        back = bkg.back() 
        rms = bkg.rms()  
        hdumain  = fits.PrimaryHDU(back,header=header)
        hdubk  = fits.ImageHDU(back)
        hdurms  = fits.ImageHDU(rms)
        hdulist = fits.HDUList([hdumain,hdubk,hdurms])
        hdulist.writeto(output+"/skyprop.fits",clobber=True)

    #extracting sources at nsigma
    thresh = nsig * bkg.globalrms
    segmap = np.zeros((header["NAXIS1"],header["NAXIS2"]))
    objects,segmap=sep.extract(data,thresh,segmentation_map=True,
                               minarea=minarea,clean=clean,mask=badmask)
    print "Extracted {} objects... ".format(len(objects))
    
    if(spectra):
        if not os.path.exists(outspec):
            os.makedirs(outspec)

    if((check) | (spectra)):
        #create a detection mask alla cubex
        srcmask=np.zeros((1,data.shape[0],data.shape[1]))
        nbj=1
        print('Generating spectra...')
        #loop over detections
        for obj in objects:
            #init mask
            tmpmask=np.zeros((data.shape[0],data.shape[1]),dtype=np.bool)
            tmpmask3d=np.zeros((1,data.shape[0],data.shape[1]),dtype=np.bool)
            #fill this mask
            sep.mask_ellipse(tmpmask,obj['x'],obj['y'],obj['a'],obj['b'],obj['theta'],r=2)
            tmpmask3d[0,:,:]=tmpmask[:,:]
            srcmask=srcmask+tmpmask3d*nbj
            if(spectra):
                savename="{}/id{}.fits".format(outspec,nbj)
                utl.cube2spec(cube,obj['x'],obj['y'],None,write=savename,
                              shape='mask',helio=helio,mask=tmpmask3d,tovac=True)
            #go to next
            nbj=nbj+1

    if(check):
        print 'Dumping source mask...'
        hdumain  = fits.PrimaryHDU(srcmask,header=header)
        hdubk  = fits.ImageHDU(srcmask)
        hdulist = fits.HDUList([hdumain,hdubk])
        hdulist.writeto(output+"/source.fits",clobber=True)
        
        print 'Dumping segmentation map'
        hdumain  = fits.PrimaryHDU(segmap,header=header)
        hdubk  = fits.ImageHDU(segmap)
        hdulist = fits.HDUList([hdumain,hdubk])
        hdulist.writeto(output+"/segmap.fits",clobber=True)
    

    #write source catalogue
    print 'Writing catalogue..'
    cols = fits.ColDefs(objects)
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(output+'/catalogue.fits',clobber=True)

    print 'All done'
    return objects



def sourcephot(catalogue,image,segmap,detection,instrument='MUSE',dxp=0.,dyp=0.,noise=[False],zpab=False, kn=1.0):

    """ 

    Get a source catalogue from findsources and a fits image with ZP
    and compute magnitudes in that filter 

    catalogue -> source cat from findsources
    image     -> fits image with ZP in header
    segmap    -> fits of segmentation map 
    detection -> the detection image, used to compute Kron radius 

    instrument -> if not MUSE, map positions from detection to image

    dxp,dyp    -> shifts in pixel of image to register MUSE and image astrometry   

    
    noise      -> if set to a noise model, use equation noise[0]*noise[1]*npix**noise[2]
                  to compute the error

    zpab  -> if ZPAB (zeropoint AB) not stored in header, must be supplied

    kn   -> factor to be used in Kron aperture

  
    """  

    from astropy.io import fits
    import numpy as np
    import sep
    import matplotlib.pyplot as plt
    from astropy.table import Table
    from astropy import wcs 

    #open the catalogue/fits 
    cat=fits.open(catalogue)
    img=fits.open(image)
    seg=fits.open(segmap)
    det=fits.open(detection)

    #grab reference wcs
    wref=wcs.WCS(det[0].header)
    psref=wref.wcs.pc[1,1]*3600.
     
    #if not handling MUSE, special cases
    if('MUSE' not in instrument):
        #handle instrument cases
        if('LRIS' in instrument):
            imgdata=img[1].data
            #place holder.. will use noise model
            #afterwards...
            vardata=imgdata*0+1
            vardata=vardata.byteswap(True).newbyteorder()
            #grab wcs image
            wimg=wcs.WCS(img[1].header)
            psimg=wimg.wcs.cd[1,1]*3600.
            #store the ZP 
            if(zpab):
                img[0].header['ZPAB']=zpab
        else:
            print 'Instrument not supported!!'
            exit()
    else:
        imgdata=img[0].data
        vardata=img[1].data
        psimg=psref

    #grab flux and var
    dataflx=np.nan_to_num(imgdata.byteswap(True).newbyteorder())
    datavar=np.nan_to_num(vardata.byteswap(True).newbyteorder())
    #grab detection and seg mask 
    detflx=np.nan_to_num(det[0].data.byteswap(True).newbyteorder())
    segmask=np.nan_to_num(seg[0].data.byteswap(True).newbyteorder())

    #if needed, map the segmap to new image with transformation
    if('MUSE' not in instrument):
        #allocate space for transformed segmentation map
        segmasktrans=np.zeros(dataflx.shape)
        print "Remapping segmentation map to new image..."

        #loop over original segmap and map to trasformed one
        #Just use nearest pixel, and keep only 1 when multiple choices 
        for xx in range(segmask.shape[0]):
            for yy in range(segmask.shape[1]):
                #go to world
                radec=wref.wcs_pix2world([[yy,xx]],0)
                #back to new instrument pixel 
                newxy=wimg.wcs_world2pix(radec,0)
                #apply shift to register WCS
                newxy[0][1]=newxy[0][1]+dyp
                newxy[0][0]=newxy[0][0]+dxp
                segmasktrans[newxy[0][1],newxy[0][0]]=segmask[xx,yy]
                
                #grow buffer as needed by individual instruments
                #This accounts for resampling to finer pixel size
                if('LRIS' in instrument):
                    segmasktrans[newxy[0][1]+1,newxy[0][0]+1]=segmask[xx,yy]
                    segmasktrans[newxy[0][1]-1,newxy[0][0]-1]=segmask[xx,yy]
                    segmasktrans[newxy[0][1]+1,newxy[0][0]-1]=segmask[xx,yy]
                    segmasktrans[newxy[0][1]-1,newxy[0][0]+1]=segmask[xx,yy]
                    segmasktrans[newxy[0][1]+1,newxy[0][0]]=segmask[xx,yy]
                    segmasktrans[newxy[0][1]-1,newxy[0][0]]=segmask[xx,yy]
                    segmasktrans[newxy[0][1],newxy[0][0]-1]=segmask[xx,yy]
                    segmasktrans[newxy[0][1],newxy[0][0]+1]=segmask[xx,yy]
                 
        #dump the transformed segmap for checking 
        hdumain  = fits.PrimaryHDU(segmasktrans,header=img[1].header)
        hdulist = fits.HDUList(hdumain)
        hdulist.writeto("segremap_check.fits",clobber=True)
    else:
        segmasktrans=segmask

    #source to extract
    nsrc=len(cat[1].data)
    phot = Table(names=('ID', 'MAGAP', 'MAGAP_ERR','FXAP', 'FXAP_ERR', 
                        'RAD', 'MAGSEG', 'MAGSEG_ERR', 'FXSEG', 'FXSEG_ERR','ZP'), 
                 dtype=('i4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4'))
    
   
    #create check aperture mask 
    checkaperture=np.zeros(dataflx.shape)
    print 'Computing photometry for objects...'

    for idobj in range(nsrc):
        
        #########
        #Find positions etc and transform as appropriate
        #########
                
        #extract MUSE source paramaters 
        x= cat[1].data['x'][idobj]
        y= cat[1].data['y'][idobj]
        a= cat[1].data['a'][idobj]
        b= cat[1].data['b'][idobj]
        theta= cat[1].data['theta'][idobj]

        #compute kron radius on MUSE detection image 
        #Kron rad in units of a,b
        tmpdata=np.copy(detflx)
        tmpmask=np.copy(segmask)
        kronrad, flg = sep.kron_radius(tmpdata, x, y, a, b, theta, 6.0, mask=tmpmask)

        #now check if size is sensible in units of MUSE data 
        rmin = 2.0  #MUSE pix 
        use_circle = kronrad * np.sqrt(a*b) < rmin

        #now use info to compute photometry and apply 
        #spatial transformation if needed
        if('MUSE' not in instrument):
            #map centre of aperture - +1 reference
            #go to world
            radec=wref.wcs_pix2world([[x,y]],1)
            #back to new instrument pixel 
            newxy=wimg.wcs_world2pix(radec,1)
            #apply shift to register WCS
            xphot=newxy[0][0]+dxp
            yphot=newxy[0][1]+dyp
                      
            #scale radii to new pixel size 
            rminphot=rmin*psref/psimg
            aphot=a*psref/psimg
            bphot=b*psref/psimg
            #Kron radius in units of a,b

        else:
            #for muse, transfer to same units
            xphot=x
            yphot=y
            rminphot=rmin
            aphot=a
            bphot=b     
            
        #####
        #Compute local sky 
        #####
        skyreg=kn*kronrad*np.sqrt(aphot*bphot)+10
        cutskymask=segmasktrans[xphot-skyreg:xphot+skyreg,yphot-skyreg:yphot+skyreg]
        cutskydata=dataflx[xphot-skyreg:xphot+skyreg,yphot-skyreg:yphot+skyreg]
        skymedian=np.nan_to_num(np.median(cutskydata[np.where(cutskymask < 1.0)]))

        #########
        #Now grab the Kron mag computed using detection image
        #########
   
        #mask all other objects to avoid blending   
        tmpdata=np.copy(dataflx)
        #apply local sky subtraction 
        tmpdata=tmpdata-skymedian
        tmpvar=np.copy(datavar)
        tmpmask=np.copy(segmasktrans)
        pixels=np.where(tmpmask == idobj+1) 
        tmpmask[pixels]=0

        if(use_circle):
            
            #circular aperture
            flux_kron, err, flg = sep.sum_circle(tmpdata, xphot, yphot, rminphot,subpix=1,mask=tmpmask)
            fluxvar, err, flg = sep.sum_circle(tmpvar,xphot,yphot,rminphot,subpix=1,mask=tmpmask)
            #store Rused in arcsec
            rused=rminphot*psimg

            #update check aperture
            tmpcheckaper=np.zeros(dataflx.shape,dtype=bool)
            sep.mask_ellipse(tmpcheckaper,xphot,yphot,1.,1.,0.,r=rminphot)
            checkaperture=checkaperture+tmpcheckaper*(idobj+1)

        else:

            #kr mag
            flux_kron, err, flg = sep.sum_ellipse(tmpdata,xphot, yphot, aphot, bphot, theta, kn*kronrad,
                                                  subpix=1,mask=tmpmask)            
            fluxvar, err, flg = sep.sum_ellipse(tmpvar,xphot,yphot, aphot, bphot, theta, kn*kronrad,
                                                subpix=1,mask=tmpmask)
            rused=kn*kronrad*psimg*np.sqrt(aphot*bphot)

            #update check aperture
            tmpcheckaper=np.zeros(dataflx.shape,dtype=bool)
            sep.mask_ellipse(tmpcheckaper,xphot,yphot,aphot,bphot,theta,r=kn*kronrad)
            checkaperture=checkaperture+tmpcheckaper*(idobj+1)

        #compute error for aperture
        if(noise[0]):
            appix=np.where(tmpcheckaper > 0)
            errflux_kron=noise[0]*noise[1]*len(appix[0])**noise[2]
        else:
            errflux_kron=np.sqrt(fluxvar)

        #go to mag 
        try:
            mag_aper=-2.5*np.log10(flux_kron)+img[0].header['ZPAB']
        except:
            mag_aper=99.0
        try:
            errmg_aper=2.5*np.log10(1.+errflux_kron/flux_kron)
        except:
            errmg_aper=99.0
    
        #find out if non detections
        if(errflux_kron >= flux_kron):
            errmg_aper=9
            mag_aper=-2.5*np.log10(2.*errflux_kron)+img[0].header['ZPAB']
        
  
        #######
        #grab the photometry in the segmentation map 
        #####

        #This may not work well for other instruments 
        #if images are not well aligned
        pixels=np.where(segmasktrans == idobj+1) 
        #add flux in pixels
        tmpdata=np.copy(dataflx)
        #apply sky sub
        tmpdata=tmpdata-skymedian
        flux_seg=np.sum(tmpdata[pixels])
        
        #compute noise from model or adding variance 
        if(noise[0]):
            #from model 
            errfx_seg=noise[0]*noise[1]*len(pixels[0])**noise[2]
        else:
            #add variance in pixels to compute error
            errfx_seg=np.sqrt(np.sum(datavar[pixels]))
  
        #go to mag with calibrations 
        try:
            mag_seg=-2.5*np.log10(flux_seg)+img[0].header['ZPAB']
        except:
            mag_seg=99.0
        try:
            errmg_seg=2.5*np.log10(1.+errfx_seg/flux_seg)
        except:
            errmg_seg=99.0
      
        #find out if non detections
        if(errfx_seg >= flux_seg):
            errmg_seg=9
            mag_seg=-2.5*np.log10(2.*errfx_seg)+img[0].header['ZPAB']
        
        #stash by id
        phot.add_row((idobj+1,mag_aper,errmg_aper,flux_kron,errflux_kron,rused,mag_seg,errmg_seg,
                      flux_seg,errfx_seg,img[0].header['ZPAB']))

    #dump the aperture check image 
    hdumain  = fits.PrimaryHDU(checkaperture,header=img[1].header)
    hdulist = fits.HDUList(hdumain)
    hdulist.writeto("aperture_check.fits",clobber=True)

    #close
    cat.close()
    img.close()
    seg.close()
    det.close()

    return phot
