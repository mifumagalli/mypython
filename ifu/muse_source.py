"""
These are sets of utilities to handle muse sources

"""

def findsources(image,cube,varima=None,check=False,output='./',spectra=False,helio=0,nsig=2.,
                minarea=10.,deblend_cont=0.0001,regmask=None,invregmask=False,fitsmask=None,clean=True,
                outspec='Spectra',marz=False,rphot=False, detphot=False, sname='MUSE'):

    """      

    Take a detection image (collapse of a cube), or median 
    of an RGB, or whatever you want (but aligned to the cube)
    and run sourceextractor 

   
    Use SEP utilities http://sep.readthedocs.org/en/stable/

    image   -> fits file of image to process
    cube    -> the cube used to extract spectra
    varima  -> the noise image corresponding to the science image (std), optional
    check   -> if true write a bunch of check mages
    output  -> where to dump the output
    spectra -> if True, extract spectra in VACUUM wave!!
    helio   -> pass additional heliocentric correction
    nsig    -> number of skyrms used for source id 
    minarea -> minimum area for extraction 
    regmask -> ds9 region file (image) of regions to be masked before extraction [e.g. edges]
    invregmask -> if True invert the mask (region defines good area)
    fitsmask -> Fits file with good mask, overrides regmask
    clean   -> clean souces 
    outspec -> where to store output spectra 
    marz    -> write spectra in also marz format (spectra needs to be true). 
               If set to numerical value, this is used as an r-band magnitude limit.
    detphot -> perform aperture phtometry on the detection image and add magnitues to the catalogue            
    rphot   -> perform r-band aperture photometry and add r-band magnitudes to the catalogue
    sname   -> prefix for the source names. Default = MUSE

    """

    import sep
    from astropy.io import fits
    from astropy import wcs
    from astropy import coordinates
    from astropy import units as u
    from astropy import table
    import numpy as np
    import os
    from mypython.ifu import muse_utils as utl
    from mypython.fits import pyregmask as msk
    from shutil import copyfile
    import glob

    #open image
    img=fits.open(image)
    header=img[0].header
    imgwcs = wcs.WCS(header)
    try:
        #this is ok for narrow band images 
        data=img[1].data
    except:
        #white cubex images
        data=img[0].data
        
    data=data.byteswap(True).newbyteorder()
    #grab effective dimension
    nex,ney=data.shape
    #close fits
    img.close()
    
    if(varima):
      var=fits.open(varima)
      try:
          datavar=var[1].data
      except:
          datavar=var[0].data
          
      datavar=datavar.byteswap(True).newbyteorder()
      #grab effective dimension
      stdx,stdy=datavar.shape
      #close fits
      var.close()
      
      if (stdx != nex) or (stdy != ney):
        print("The noise image does not have the same dimensions as the science image")
        return -1

    #create bad pixel mask
    if(fitsmask):
        print("Using FITS image for badmask")
        hdumsk = fits.open(fitsmask)
        try:
          badmask = hdumsk[1].data
        except:
          badmask = hdumsk[0].data
        badmask=badmask.byteswap(True).newbyteorder()
    elif(regmask):
        print("Using region file for badmask")
        Mask=msk.PyMask(ney,nex,regmask,header=img[0].header)
        for ii in range(Mask.nreg):
            Mask.fillmask(ii)
            if(ii == 0):
                badmask=Mask.mask
            else:
                badmask+=Mask.mask
            badmask=1.*badmask
    else:
        badmask=np.zeros((nex,ney))
    
    if (regmask) and (invregmask) and not (fitsmask):
        badmask = 1-badmask
    
    if(check):
        print('Dumping badmask')
        hdumain  = fits.PrimaryHDU(badmask,header=header)
        hdulist = fits.HDUList([hdumain])
        hdulist.writeto(output+"/badmask.fits",overwrite=True)
    

    #check background level, but do not subtract it
    print('Checking background levels')
    bkg = sep.Background(data,mask=badmask)    
    print('Residual background level ', bkg.globalback)
    print('Residual background rms ', bkg.globalrms)

    if(check):
        print('Dumping sky...')
        #dump sky properties
        back = bkg.back() 
        rms = bkg.rms()  
        hdumain  = fits.PrimaryHDU(back,header=header)
        hdubk  = fits.ImageHDU(back)
        hdurms  = fits.ImageHDU(rms)
        hdulist = fits.HDUList([hdumain,hdubk,hdurms])
        hdulist.writeto(output+"/skyprop.fits",overwrite=True)

    if(varima):
        #Use nsigma threshold and a pixel by pixel effective treshold based on variance map
        thresh = nsig
        objects,segmap=sep.extract(data,thresh,var=datavar,segmentation_map=True,
                               minarea=minarea,clean=clean,mask=badmask,deblend_cont=deblend_cont,deblend_nthresh=32)
    else:
        #extracting sources at nsigma, use constant threshold
        thresh = nsig * bkg.globalrms
        objects,segmap=sep.extract(data,thresh,segmentation_map=True,
                               minarea=minarea,clean=clean,mask=badmask,deblend_cont=deblend_cont,deblend_nthresh=32)
    
    
    print("Extracted {} objects... ".format(len(objects)))
    ids  = np.arange(len(objects))+1    
    
    if(spectra):
        if not os.path.exists(outspec):
            os.makedirs(outspec)

    if((check) | (spectra)):
        #create a detection mask a'la cubex
        srcmask=np.zeros((data.shape[0],data.shape[1]))
        print('Generating spectra...')
        #loop over detections
        for nbj in ids:
            obj = objects[nbj-1]
            #init mask
            tmpmask=np.zeros((data.shape[0],data.shape[1]),dtype=np.bool)
            #fill this mask
            sep.mask_ellipse(tmpmask,obj['x'],obj['y'],obj['a'],obj['b'],obj['theta'],r=2)
            #add in global mask
            srcmask=srcmask+tmpmask*nbj
            #verify conflicts, resolve using segmentation map
            if np.nanmax(srcmask)>nbj:
               blended = (srcmask>nbj)
               srcmask[blended] = segmap[blended]
            
        #Now loop again and extract spectra if required
        if(spectra):
           #Verify that the source mask has the same number of objects as the object list
           if not len(np.unique(srcmask[srcmask>0])) == len(objects):
              print("Mismatch between number of objects and number of spectra to extract.")
           for nbj in ids:
              savename="{}/id{}.fits".format(outspec,nbj)
              tmpmask3d=np.zeros((1,data.shape[0],data.shape[1]))
              tmpmask3d[0,:,:]=srcmask[:,:]
              tmpmask3d[tmpmask3d != nbj] = 0
              tmpmask3d[tmpmask3d > 0] = 1
              tmpmask3d = np.array(tmpmask3d, dtype=np.bool)
              utl.cube2spec(cube,None,None,None,write=savename,shape='mask',helio=helio,mask=tmpmask3d,tovac=True)

    if(check):
        print('Dumping source mask...')
        hdumain  = fits.PrimaryHDU(srcmask,header=header)
        hdubk  = fits.ImageHDU(srcmask)
        hdulist = fits.HDUList([hdumain,hdubk])
        hdulist.writeto(output+"/source.fits",overwrite=True)
        
        print('Dumping segmentation map')
        hdumain  = fits.PrimaryHDU(segmap,header=header)
        hdubk  = fits.ImageHDU(segmap)
        hdulist = fits.HDUList([hdumain,hdubk])
        hdulist.writeto(output+"/segmap.fits",overwrite=True)
    
    #Generate source names using coordinates and name prefix
    ra, dec = imgwcs.wcs_pix2world(objects['x'], objects['y'],0)
    coord = coordinates.FK5(ra*u.degree, dec*u.degree)
    rastr  = coord.ra.to_string(u.hour, precision=2, sep='', pad=True)
    decstr = coord.dec.to_string(u.degree, precision=1, sep='', alwayssign=True, pad=True)
    name = [sname+'J{0}{1}'.format(rastr[k], decstr[k]) for k in range(len(rastr))]
    
    
    #Generate a column to be used to flag the sources to be used in the analysis
    #True for all sources at this point
    use_source = np.ones_like(name, dtype=bool)
    
    #write source catalogue
    print('Writing catalogue..')
    tab = table.Table(objects)
    tab.add_column(table.Column(dec),0,name='DEC')
    tab.add_column(table.Column(ra),0,name='RA')
    tab.add_column(table.Column(name),0,name='name')
    tab.add_column(table.Column(ids),0,name='ID')
    tab.add_column(table.Column(use_source),name='use_source')
    tab.write(output+'/catalogue.fits',overwrite=True)
    
    if (detphot):
       #Run source photometry on the extraction image
       whiteimg, whitevar, whitewcsimg = utl.cube2img(cube, write=output+'/Image_white.fits')
       phot_det = sourcephot(output+'/catalogue.fits', output+'/Image_white.fits', output+'/segmap.fits', image, zpab=28.35665)
       phot_det.add_column(table.Column(name),1,name='name')
       tbhdu = fits.open(output+'/catalogue.fits')
       tbhdu.append(fits.BinTableHDU(phot_det))
       tbhdu[-1].header['PHOTBAND'] = 'Detection'
       tbhdu.writeto(output+'/catalogue.fits',overwrite=True)
    
    #rband photometry
    if (rphot):
        rimg, rvar, rwcsimg = utl.cube2img(cube, filt=129, write=output+'/Image_R.fits')
        phot_r = sourcephot(output+'/catalogue.fits', output+'/Image_R.fits', output+'/segmap.fits', image)
        phot_r.add_column(table.Column(name),1,name='name')

        tbhdu = fits.open(output+'/catalogue.fits')
        tbhdu.append(fits.BinTableHDU(phot_r))
        tbhdu[-1].header['PHOTBAND'] = 'SDSS_r'
        tbhdu.writeto(output+'/catalogue.fits',overwrite=True)  
    
    if((marz) & (spectra)):
        #if marz is True but no magnitude limit set, create marz file for whole catalogue
        if marz>10 and (rphot):
            #Requires testing
            hdu = fits.open(output+'/catalogue.fits')
            hdu[1].data['use_source'][hdu[2].data['MAGAP']>marz] = False
            hdu.writeto(output+'/catalogue.fits',overwrite=True)
            
            marz_file(output+'/catalogue.fits', outspec, output, r_lim=marz)
        else:
            marz_file(output+'/catalogue.fits', outspec, output)
        
    print('All done')
    return objects
    

def marz_file(catalogue, specdir, outdir, specfmt='id', r_lim=False):
    
    import glob
    from astropy.io import fits
    from astropy import wcs
    import numpy as np
    from mypython.ifu import muse_utils as utl

    
    #Makes a list of spectra files ** MUST include the spectra from catalogue**
    filelist = glob.glob(specdir+'/{}*'.format(specfmt))
    
    #Need catalogue of spectra objects            
    catalog = fits.open(catalogue)   
    
    #name of MARZ file
    if r_lim:
        marzfile = outdir+'/spectra_marz_r' + str(r_lim) + '.fits'    
    else:
        marzfile = outdir+'/spectra_marz.fits'                                  
    
        
    s0 = fits.open(filelist[0])
    naxis3=len(s0[0].data)

    id         = catalog[1].data['ID']         # Integer array - id of object
    name       = catalog[1].data['name']
    x          = catalog[1].data['x']          # Array of object image x-coords
    y          = catalog[1].data['y']          # Array of object image y-coords
    ra         = catalog[1].data['RA']         # Array of object R.A.s
    dec        = catalog[1].data['DEC']        # Array of object Dec.s
    use_source = catalog[1].data['use_source']
    
    nspectra = len(id[use_source])
    
    id   = id[use_source]
    name = name[use_source]
    x    = x[use_source]
    y    = y[use_source]
    ra   = ra[use_source]
    dec  = dec[use_source]   
    
    print("Generating MARZ file for {} spectra".format(nspectra))
        
    # Initialize flux, variance & sky arrays (sky not mandatory)
    intensity = np.zeros((nspectra,naxis3))
    variance = np.zeros((nspectra,naxis3))
    sky = np.zeros((nspectra,naxis3))
    wavelength = np.zeros((nspectra,naxis3))

    # Fill the flux, variance and sky arrays here.
    # Apparently MARZ does not really like vacuum waves so we turn them back into air
    type = []

    for ii in range(nspectra):
        data = fits.open(specdir+'/{}{}.fits'.format(specfmt,id[ii]))
        type.append('P')
        intensity[ii,:]  = data[0].data
        variance[ii,:]   = data[1].data
        sky[ii,:]        = data[0].data*0.0
        wavelength[ii,:] = utl.vactoair(data[2].data)

    intensity[np.logical_not(np.isfinite(intensity))] = np.nan
    variance[np.logical_not(np.isfinite(variance))] = np.nan
    sky[np.logical_not(np.isfinite(sky))] = 0.0

    # Set-up the fits file
    marz_hdu = fits.HDUList()
    marz_hdu.append(fits.ImageHDU(intensity))
    marz_hdu.append(fits.ImageHDU(variance))
    marz_hdu.append(fits.ImageHDU(sky))
    marz_hdu.append(fits.ImageHDU(wavelength))
    marz_hdu[0].header.set('extname', 'INTENSITY')
    marz_hdu[1].header.set('extname', 'VARIANCE')
    marz_hdu[2].header.set('extname', 'SKY')
    marz_hdu[3].header.set('extname', 'WAVELENGTH')
    
    # Add in source parameters as fits table
    c1 = fits.Column(name='NAME', format='80A', array=name)
    c2 = fits.Column(name='RA', format='D', array=ra)
    c3 = fits.Column(name='DEC', format='D',array=dec)
    c4 = fits.Column(name='X', format='J',array=x)
    c5 = fits.Column(name='Y', format='J', array=y)
    c6 = fits.Column(name='TYPE', format='1A', array=type)
    coldefs = fits.ColDefs([c1, c2, c3, c4, c5,c6])
    marz_hdu.append(fits.BinTableHDU.from_columns(coldefs))
    marz_hdu[4].header.set('extname', 'FIBRES')

    # And write out.
    marz_hdu.writeto(marzfile,overwrite=True)


def sourcephot(catalogue,image,segmap,detection,instrument='MUSE',dxp=0.,dyp=0.,
               noise=[False],zpab=False, kn=2.5, circap=1.0):

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

    kn   -> factor to be used when scaling Kron apertures [sextractor default 2.5]
  
    circap -> radius in arcsec for aperture photmetry to be used when Kron aperture fails 

    """  

    from astropy.io import fits
    import numpy as np
    import sep
    import matplotlib.pyplot as plt
    from astropy.table import Table
    from astropy import wcs 


    #grab root name 
    rname=((image.split('/')[-1]).split('.fits'))[0]
    print('Working on {}'.format(rname))

    #open the catalogue/fits 
    cat=fits.open(catalogue)
    img=fits.open(image)
    seg=fits.open(segmap)
    det=fits.open(detection)

    #grab reference wcs from detection image 
    wref=wcs.WCS(det[0].header)
    psref=wref.pixel_scale_matrix[1,1]*3600.
    print('Reference pixel size {}'.format(psref))


    #if not handling MUSE, special cases for format of data
    if('MUSE' not in instrument):
        #handle instrument cases
        if('LRIS' in instrument):
            #data 
            imgdata=img[1].data
            #place holder for varaince as will use noise model below
            vardata=imgdata*0+1
            vardata=vardata.byteswap(True).newbyteorder()
            #grab wcs image
            wimg=wcs.WCS(img[1].header)
            psimg=wimg.pixel_scale_matrix[1,1]*3600.
            #store the ZP 
            if(zpab):
                img[0].header['ZPAB']=zpab
        else:
            print('Instrument not supported!!')
            exit()
    else:
        #for muse, keep eveything the same
        psimg=psref
        try:
            imgdata=img[1].data
            vardata=img[2].data
        except:
            imgdata=img[0].data
            vardata=img[1].data

    #grab flux and var
    dataflx=np.nan_to_num(imgdata.byteswap(True).newbyteorder())
    datavar=np.nan_to_num(vardata.byteswap(True).newbyteorder())
    #grab detection and seg mask 
    detflx=np.nan_to_num(det[0].data.byteswap(True).newbyteorder())
    #go back to 1d
    if(len(seg[0].data.shape)>2):
        segmask=(np.nan_to_num(seg[0].data.byteswap(True).newbyteorder()))[0,:,:]
    else:
        segmask=(np.nan_to_num(seg[0].data.byteswap(True).newbyteorder()))


    #if needed, map the segmap to new image with transformation
    if('MUSE' not in instrument):
        #allocate space for transformed segmentation map
        segmasktrans=np.zeros(dataflx.shape)
        print("Remapping segmentation map to new image...")

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
        hdulist.writeto("{}_segremap.fits".format(rname),overwrite=True)
    else:
        #no transformation needed
        segmasktrans=segmask

    #source to extract
    nsrc=len(cat[1].data)
    print('Extract photometry for {} sources'.format(nsrc))
    phot = Table(names=('ID', 'MAGAP', 'MAGAP_ERR','FXAP', 'FXAP_ERR', 
                        'RAD', 'MAGSEG', 'MAGSEG_ERR', 'FXSEG', 'FXSEG_ERR','ZP'), 
                 dtype=('i4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4'))
    
   
    #create check aperture mask 
    checkaperture=np.zeros(dataflx.shape)
    print('Computing photometry for objects...')

    #loop over each source
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
        #mask all other sources to avoid overlaps but keep desired one
        pixels=np.where(tmpmask == idobj+1) 
        tmpmask[pixels]=0

        #compute kron radius [pixel of reference image]
        kronrad, flg = sep.kron_radius(tmpdata,x,y,a,b,theta,6.0,mask=tmpmask)

        #plt.imshow(np.log10(tmpdata+1),origin='low')
        #plt.show()
        #exit()

        #now check if size is sensible in units of MUSE data 
        rmin = 2.0  #MUSE pix 
        use_circle = kronrad * np.sqrt(a*b) < rmin
      
        #use circular aperture of 2" in muse pixel unit
        rcircap = circap/psref
        
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
            rminphot=rcircap*psref/psimg
            aphot=a*psref/psimg
            bphot=b*psref/psimg
            #Kron radius in units of a,b

        else:
            #for muse, transfer to same units
            xphot=x
            yphot=y
            rminphot=rcircap
            aphot=a
            bphot=b     
            
        #####
        #Compute local sky 
        #####
        skyreg=kn*kronrad*np.sqrt(aphot*bphot)+15
        if (yphot-skyreg < 0.0): yphot=skyreg
        if (xphot-skyreg < 0.0): xphot=skyreg
        if (yphot+skyreg > segmasktrans.shape[0]-1): yphot=segmasktrans.shape[0]-1-skyreg
        if (xphot+skyreg > segmasktrans.shape[1]-1): xphot=segmasktrans.shape[1]-1-skyreg
        #print(int(yphot-skyreg),int(yphot+skyreg),int(xphot-skyreg),int(xphot+skyreg))
        cutskymask=segmasktrans[int(yphot-skyreg):int(yphot+skyreg),int(xphot-skyreg):int(xphot+skyreg)]
        cutskydata=dataflx[int(yphot-skyreg):int(yphot+skyreg),int(xphot-skyreg):int(xphot+skyreg)]
        skymedian=np.nan_to_num(np.median(cutskydata[np.where(cutskymask < 1.0)]))

        #print skymedian    

        #plt.imshow(cutskymask,origin='low')
        #plt.show()
        #if(idobj > 30):
        #    exit()


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

        #plt.imshow(tmpmask,origin='low')
        #plt.show()
        #exit()

        #circular aperture
        if(use_circle):        
           
            #flux in circular aperture
            flux_kron, err, flg = sep.sum_circle(tmpdata,xphot,yphot,rminphot,mask=tmpmask)
            #propagate variance
            fluxvar, err, flg = sep.sum_circle(tmpvar,xphot,yphot,rminphot,mask=tmpmask)
            #store Rused in arcsec
            rused=rminphot*psimg

            #update check aperture
            tmpcheckaper=np.zeros(dataflx.shape,dtype=bool)
            sep.mask_ellipse(tmpcheckaper,xphot,yphot,1.,1.,0.,r=rminphot)
            checkaperture=checkaperture+tmpcheckaper*(idobj+1)

        #kron apertures 
        else:
            #kron flux 
            flux_kron, err, flg = sep.sum_ellipse(tmpdata,xphot, yphot, aphot, bphot, theta, kn*kronrad,
                                                  mask=tmpmask)            
            #propagate variance 
            fluxvar, err, flg = sep.sum_ellipse(tmpvar,xphot,yphot, aphot, bphot, theta, kn*kronrad,
                                                mask=tmpmask)
            #translate in radius
            rused=kn*kronrad*psimg*np.sqrt(aphot*bphot)

            #update check aperture
            tmpcheckaper=np.zeros(dataflx.shape,dtype=bool)
            sep.mask_ellipse(tmpcheckaper,xphot,yphot,aphot,bphot,theta,r=kn*kronrad)
            checkaperture=checkaperture+tmpcheckaper*(idobj+1)

        #compute error for aperture
        if(noise[0]):
            #use model 
            appix=np.where(tmpcheckaper > 0)
            errflux_kron=noise[0]*noise[1]*len(appix[0])**noise[2]
        else:
            #propagate variance 
            errflux_kron=np.sqrt(fluxvar)

        #go to mag 
        if(flux_kron > 0):
            mag_aper=-2.5*np.log10(flux_kron)+img[0].header['ZPAB']
            errmg_aper=2.5*np.log10(1.+errflux_kron/flux_kron)
        else:
            mag_aper=99.0
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
        if(flux_seg > 0):
            mag_seg=-2.5*np.log10(flux_seg)+img[0].header['ZPAB']
            errmg_seg=2.5*np.log10(1.+errfx_seg/flux_seg)     
        else:
            mag_seg=99.0
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
    hdulist.writeto("{}_aper.fits".format(rname),overwrite=True)

    #close
    cat.close()
    img.close()
    seg.close()
    det.close()

    return phot


def forcephot(image,x,y,rad,skyreg=[10,20],show=False,mask=None):
    
    """

    Compute the photometry on a MUSE reconstructed image within an aperture 
    Treat limits as 2sigma

    x,y        -> position of the aperture in pixels
    radap      -> radius of the aperture in pixels
    skyreg     -> radius of inner/outer sky region in pixel
    mask       -> mask to avoid pixels
    
    show       -> if true, show what's appening


    """
    
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import sep
    import numpy as np

    #open the image
    image=fits.open(image)
    data = image[1].data.byteswap().newbyteorder()
    var = image[2].data.byteswap().newbyteorder()

    #grab the zp
    zp=image[0].header['ZPAB']
    flux, err, flg = sep.sum_circle(data,x,y,rad,var=var,bkgann=skyreg,mask=mask)

    #compute magnitudes
    if(flux < 2*err):
        mag=-2.5*np.log10(2*err)+zp
        errmag=99
    else:
        mag=-2.5*np.log10(flux)+zp
        errmag=2.5*np.log10(1.+err/flux)     

    #if show
    if(show):
        

        #display field
        median=np.median(data)
        stddev=np.std(data)
        plt.imshow(data,clim=[median-0.2*stddev,median+stddev],origin='low')
        

        #now the aperture 
        circ=plt.Circle((x,y),radius=rad, color='red', fill=False)
        ax=plt.gca()
        ax.add_patch(circ)
    
        if(skyreg):
             circ1=plt.Circle((x,y),radius=skyreg[0], color='red', fill=False)
             circ2=plt.Circle((x,y),radius=skyreg[1], color='red', fill=False)
             ax.add_patch(circ1)
             ax.add_patch(circ2)

        

        plt.show()



    return mag, errmag


def sourceimgspec(cubes,tags,objidmask,specube=None,vacuum=True,outpath='imgspecs',\
                     variance=None, scalevar=None):
    
    
    """

    This code takes one (cube) or more cubes together with a 3D segmentation 
    map produced by cubex and generates images and spectra for inspection 
    
    Cube2Im need to be installed and in the path

    cubes      -> list of short cubes (matched to objid dimenion) for image extraction
    tags       -> name tags to be used for *_imf.fits output
    objidmask  -> the 3d segmentation map used for extraction
    specube -> if set, used to extract spectrum (e.g. full cube). If not 
                available, use first cube in the list of spectrum generation
    vacuum  -> if true, convert wave to wacuum 
    outpath -> where to store products

    variance -> if set to a list of variance cubes [1st extension], propagate variance inside the 
                projected 2D images

    scalevar -> if set to cubex variance rescale file (z,scalefactor), 
                apply variance scaling before propagating variance 

    """
    
    from astropy.io import fits 
    import os
    import numpy as np
    import subprocess
    from mypython.ifu import muse_utils as mutl
    from astropy.table import Table

    
    #sanity checks 
    #turn in list what should be list and not string
    if(isinstance(cubes, basestring)):
        cubes=[cubes]
        tags=[tags]
    #check tags match size of cubes
    if(len(cubes) != len(tags)):
        raise ValueError("Cubes and tags do not match in size!")

    #set cube to be used for spectrum extraction
    if(specube):
        cubespec=specube
    else:
        cubespec=cubes[0]

    #stash the variance if needed
    allvar=[]
    if(variance):
        print('Prepare variance data...')
        #check if need to rescale
        if(scalevar):
            #read scaling factor
            vscale=Table.read(scalevar,format='ascii')
            #grab cube dimension and construct rescaled variance
            thisv=fits.open(variance[0])
            nl,nx,ny=thisv[0].data.shape
            scale=np.transpose(np.tile(vscale['rescaling_fact'],(ny,nx,1)))
        else:
            scale=1.0
        #now store with scaling
        for varname in variance:
            thisv=fits.open(varname)
            allvar.append(thisv[0].data*scale)
            thisv.close()
            

    #make output folder
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    #grab the indexes
    segmap=fits.open(objidmask)
    objid=np.max(segmap[0].data)
    
    print('Loop over IDs...')

    for ii in range(objid):
        #first check if index is legit
        inside=np.where(segmap[0].data == ii+1)
        if(len(inside[0]) > 1):
            print('Work on ID {}'.format(ii+1))
            #make outfolder
            if not os.path.exists(outpath+'/id{}'.format(ii+1)):
                os.makedirs(outpath+'/id{}'.format(ii+1))
            currentpath=outpath+'/id{}'.format(ii+1)
            #next generate image
            for cc,thiscube in enumerate(cubes):
                subprocess.call(["Cube2Im","-cube",thiscube,"-out","{}/{}_img.fits".format(currentpath,tags[cc]),"-id","{}".format(ii+1),"-idcube",objidmask,"-nl","-1","-idpad","20","-sbscale",".true."])
                #append variance if needed
                if(variance):
                    #propagate the variance along lambda axis
              
                    #find star-end in each dimension
                    ilstart=np.min(inside[0])
                    ilend=np.max(inside[0])
                    dl=ilend-ilstart+1
                    ixstart=np.min(inside[1])
                    ixend=np.max(inside[1])
                    dx=ixend-ixstart+1
                    iystart=np.min(inside[2])
                    iyend=np.max(inside[2])
                    dy=iyend-iystart+1

                    #now build a mask array
                    thismask=np.zeros((dl,dx,dy))
                    thismask[inside[0]-ilstart,inside[1]-ixstart,inside[2]-iystart]=1.0
                    
                    #build variance slice
                    varslice=np.zeros((dl,dx,dy))
                    varslice[inside[0]-ilstart,inside[1]-ixstart,inside[2]-iystart]=\
                        allvar[cc][inside[0],inside[1],inside[2]]
                    
                    #now with flux 
                    flux=fits.open(thiscube)
                    fluxslice=np.zeros((dl,dx,dy))
                    fluxslice[inside[0]-ilstart,inside[1]-ixstart,inside[2]-iystart]=\
                        flux[0].data[inside[0],inside[1],inside[2]]
                    flux.close()
                    
                    #compute integral flux and SN
                    totflux=np.sum(np.nan_to_num(fluxslice)*thismask)*1e-20*1.25
                    toterr=np.sqrt(np.sum(np.nan_to_num(varslice)*thismask))*1e-20*1.25
                    print('Flux: {}'.format(totflux))
                    print('S/N: {}'.format(totflux/toterr))
                    
                    #open image to update
                    imgupdate=fits.open("{}/{}_img.fits".format(currentpath,tags[cc]),mode='update')
                    imgupdate[0].header['ISOFLUX']=np.nan_to_num(totflux)
                    imgupdate[0].header['ISOERR']=np.nan_to_num(toterr)
                    imgupdate[0].header['S2N']=np.nan_to_num(totflux/toterr)
                    imgupdate.flush()
                    imgupdate.close()
                    
            #finally call spectrum generation
            print('Extracting 1d spectrum from {}'.format(cubespec))
            mutl.cube2spec(cubespec,0,0,0,shape='mask',tovac=vacuum,idsource=ii+1,mask=segmap[0].data,\
                               write='{}/spectrum.fits'.format(currentpath))
        else:
            print('Skip index {} as not in ID map'.format(ii+1))

    segmap.close()
     
         

    




    
    




    

    
