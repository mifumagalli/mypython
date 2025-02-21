"""
These are sets of utilities to handle muse sources

"""


def mocklines(cube,segmap,fluxlimits,badmask=None,output='./',num=500,real=None,wavelimits=None,spatwidth=3.5,wavewidth=2,outprefix='mocks',fill=8., exp=False):

    """

    Inject mock line emission in a cube using a flat distribution between fluxlimits
    and a three dimensional Gaussian with x/y FWHM = spatwidth and lambda FWHM = wavewidth

    cube -> a MUSE cube [filename] to use for mock 
    fluxlimits -> mock sources will be drawn in range [min,max] in units of CUBE [default 1e-20]
    num -> number of mock sources [high numbers give better statistics but can lead to shadowing]
    wavelimits -> [min,max] slices in which mocks are populated in pixel coordinates
    spatwidth -> FWHM in spatial direction, in pixels 
    wavewidth -> FWHM in spectral direction, in slices
    outprefix -> prefix for output  
    fill -> multiple of sigma to evaluate Gaussian. Larger number is more accurate but slower
    exp -> use an exponential profile in x,y. If false, use point sources. If exp is a list the 
           values will be drawn randomly

    """
    
    from astropy.io import fits 
    import os
    import numpy as np
    from astropy import units as u
    from scipy.ndimage import filters
    
    #open the cube 
    cubhdu=fits.open(cube)
    
    try:
        newcube=cubhdu[0].data
        ext=0
    except:
        newcube=cubhdu[1].data
        var=cubhdu[2].data
        ext=1

    #Open segmentation image
    
    seghdu = fits.open(segmap)
    segima = seghdu[0].data
    seghdu.close()
    
    if segima.ndim==2:
       segcube = np.repeat(segima[np.newaxis,...], np.shape(newcube)[0], axis=0)
    else:
       segcube = segima
    
    if np.isscalar(real):
       real = [real]

    try:
      dummy = len(real)
    except:
      real = [real]  
    
    if real[0] is not None:
       
       realcubes = []
       for rr in real:
           try:
             realcubes.append(fits.open(rr)[0].data)
           except:
             pass
               
       nreal = len(realcubes)
       nzreal, nyreal, nxreal = np.shape(realcubes[0])
           
    if(badmask):
      #Open badmask image
      badhdu = fits.open(badmask)
      badmask = badhdu[0].data
      badhdu.close()      
    else:
      badmask = np.zeros_like(segima)
    
    #find ranges
    if not wavelimits:
        #auto select with gap 
        minw=10
        maxw=cubhdu[0].header['NAXIS3']-10
    if wavelimits:
        minw = wavelimits[0]
        maxw = wavelimits[1]
    
    nx, ny, nw = cubhdu[0].header['NAXIS1'], cubhdu[0].header['NAXIS2'], cubhdu[0].header['NAXIS3']
    
    minx=20
    maxx=nx-20
    miny=20
    maxy=ny-20

    #open output to store mocks
    fl=open(output+"{}_{}_catalogue.txt".format(outprefix,os.path.basename(cube).split('.fits')[0]),'w')
        
    #loop on mocks
    print('Injecting mock sources')
    
    mockcube = np.zeros_like(newcube)
    
    if np.isscalar(exp):
      exp = [exp]
      nexp = 1
    else:
      nexp = len(exp)  
        
    ind = 0
    
    while ind<num:
      
      #now draw random distributions  
      mflux    = 10**np.random.uniform(np.log10(fluxlimits[0]),np.log10(fluxlimits[1])) 
      xc       = np.random.uniform(minx,maxx)                                           
      yc       = np.random.uniform(miny,maxy)                                           
      wc       = np.random.uniform(minw,maxw)                                           

      if real[0] is None:                           
                                                                                          
        #compute Gaussian parameters                                                      
        if exp[0] is not False:                                                           
           sigmax=exp[ind%nexp]                                                           
           sigmay=sigmax                                                                  
        else:                                                                             
           sigmax=spatwidth/(2*np.sqrt(2*np.log(2)))                                      
           sigmay=sigmax                                                                  
                                                                                          
        sigmaw=wavewidth/(2*np.sqrt(2*np.log(2)))                                         

        sizex = np.clip(int(np.ceil(fill*sigmax)),3,15)                                   
        sizey = np.clip(int(np.ceil(fill*sigmay)),3,15)                                   
        sizew = np.clip(int(np.ceil(fill*sigmaw)),3,15)                                   
      
      else:
        
        sizex = int((nxreal-3)/2)                                   
        sizey = int((nyreal-3)/2)                                   
        sizew = int((nzreal-3)/2)                                   
      
      if np.ceil(xc)<=sizex or np.ceil(yc)<=sizey:
         continue
      
      #Verify availablity in seg map
      thisseg = np.sum(segcube[int(np.ceil(wc))-sizew:int(np.ceil(wc))+sizew,int(np.ceil(yc))-sizey:int(np.ceil(yc))+sizey,int(np.ceil(xc))-sizex:int(np.ceil(xc))+sizex])
      thiscub = np.sum(newcube[int(np.ceil(wc))-sizew:int(np.ceil(wc))+sizew,int(np.ceil(yc))-sizey:int(np.ceil(yc))+sizey,int(np.ceil(xc))-sizex:int(np.ceil(xc))+sizex])
      thisbad = np.sum(badmask[int(np.ceil(yc))-sizey:int(np.ceil(yc))+sizey,int(np.ceil(xc))-sizex:int(np.ceil(xc))+sizex])

      if thisseg > 0 or thisbad>0 or np.isnan(thiscub):
         continue
      else:
        #Book space in segcube with some space around to avoid overlapping injections
        segcube[int(np.ceil(wc))-2*sizew:int(np.ceil(wc))+2*sizew,int(np.ceil(yc))-2*sizey:int(np.ceil(yc))+2*sizey,int(np.ceil(xc))-2*sizex:int(np.ceil(xc))+2*sizex] = 1 
        ind += 1  

      if real[0] is None:                           

        norm=mflux/(sigmax*sigmay*sigmaw*(2*np.pi)**(3./2.))
      
        #now evaluate Gaussian or exponential at pixel center
        fl.write("{} {} {} {} {} {} {}\n".format(xc+0.5, yc+0.5, wc+0.5,mflux,sigmax, sigmay, sigmaw))

        allx=np.round(np.arange(np.floor(xc-fill*sigmax),xc+fill*sigmax,1))
        ally=np.round(np.arange(np.floor(yc-fill*sigmay),yc+fill*sigmay,1))
        allw=np.round(np.arange(np.floor(wc-fill*sigmaw),wc+fill*sigmaw,1))
        
        
        for xx in allx:
                for yy in ally:
                    for ww in allw:
                         if (xx>0) & (yy>0) & (ww>0) & (xx<nx-1) & (yy<ny-1) & (ww<nw-1):
                            if exp[0] is not False:
                                #For exponential we need to compute the radius explicitly.
                                #The code only allows symmetric exponentials in x and y, sorry about that.
                                rad = np.sqrt((xx-xc)**2+(yy-yc)**2)
                                pix= norm*np.exp(-1.*(( rad /sigmax )+
                                                      (((ww-wc)**2)/(2.*sigmaw**2))))
                            
                            else:
                                pix=norm*np.exp(-1.*((((xx-xc)**2)/(2.*sigmax**2))+
                                                 (((yy-yc)**2)/(2.*sigmay**2))+
                                                 (((ww-wc)**2)/(2.*sigmaw**2))))      
                         
                            #store
                            mockcube[int(ww),int(yy),int(xx)]=mockcube[int(ww),int(yy),int(xx)]+pix

      else:        
        
        realcubeid = ind%nreal
        realhsizew = int((nzreal-1)/2)
        realhsizey = int((nyreal-1)/2) 
        realhsizex = int((nxreal-1)/2)

        fl.write("{} {} {} {} {} {} {}\n".format(int(xc)+0.5,int(yc)+0.5,int(wc)+0.5,mflux,realcubeid,0,0))
        
        mockcube[int(wc)-realhsizew:int(wc)+realhsizew+1,int(yc)-realhsizey:int(yc)+realhsizey+1,int(xc)-realhsizex:int(xc)+realhsizex+1] += realcubes[realcubeid]*mflux                           

    
    if(exp[0] is not False and real[0] is None):
      #The mock exponential profiles need to be convolved with a gaussian 2D filter to simulate PSF effects.
      #go from FWHM to sigma for the Kernel
      kern = (0,spatwidth/(2*np.sqrt(2*np.log(2))),spatwidth/(2*np.sqrt(2*np.log(2))))
      mockcube =  filters.gaussian_filter(mockcube, kern, order=0)    
    
    newcube += mockcube

    #store output
    if(ext==0):
        hdufirst = fits.PrimaryHDU(newcube)
        hdulist = fits.HDUList([hdufirst])
        hdulist[0].header=cubhdu[0].header
    elif(ext==1):
        hdufirst = fits.PrimaryHDU([])
        hdusecond = fits.ImageHDU(newcube)
        hduthird = fits.ImageHDU(var)
        hdulist = fits.HDUList([hdufirst,hdusecond,hduthird])
        hdulist[0].header=cubhdu[0].header
        hdulist[1].header=cubhdu[1].header
        hdulist[2].header=cubhdu[2].header

    write=  output + "{}_{}".format(outprefix,os.path.basename(cube))   
    hdulist.writeto(write,overwrite=True)
    
    fl.close()
              
    return



def mockcont(image,segmap,fluxlimits,badmask=None,output='./',num=100,ZP=-1,spatwidth=3.5,outprefix='cmocks',fill=6.,exp=False, expscale=1.5):

    """

    Inject mock continuum emission in a image using a flat distribution between fluxlimits
    and a three dimensional Gaussian with x/y FWHM = spatwidth and lambda FWHM = wavewidth


    image -> a MUSE image [filename] to use for mock 
    segmap -> Segmentation map of the input image
    fluxlimits -> mock sources will be drawn in range [min,max] in units of 
                  image [default 1e-20] if ZP is -1, in mag units otherwise
    
    badmask -> (Optional) injection will not happen on masked pixels (1)
    ZP -> Zero point for magnitude to flux conversion, if -1 do the mock in flux
    num -> number of mock sources [high numbers give better statistics but can lead to shadowing]
    spatwidth -> spatial FWHM for Gaussian model in pixel
    outprefix -> prefix for temporary output files  
    fill -> multiple of sigma to evaluate Gaussian. Larger number is more accurate but slower
    exp -> use an exponential profile in x,y. If false, use point sources
    expscale -> exponential scale lenght in pixels, this is convolved with the Gaussian FWHM (seeing)

    """
    
    from astropy.io import fits 
    from astropy.convolution import convolve, Gaussian2DKernel
    import numpy as np
    import os
    
    #open the image 
    imahdu=fits.open(image)
    
    try:
        newimage=imahdu[0].data
        ext=0
    except:
        newimage=imahdu[1].data
        ext=1
    
    #Open segmentation image
    
    seghdu = fits.open(segmap)
    segima = seghdu[0].data
    seghdu.close()
    
    if(badmask):
      #Open badmask image
      badhdu = fits.open(badmask)
      badmask = badhdu[0].data
      badhdu.close()      
    else:
      badmask = np.zeros_like(segima)

    nx, ny = imahdu[0].header['NAXIS1'], imahdu[0].header['NAXIS2']
    
    minx=20
    maxx=nx-20
    miny=20
    maxy=ny-20
    
    #now draw random distributions

    #open output to store mocks
    fl=open(output+"{}_{}_catalogue.txt".format(outprefix,os.path.basename(image).split('.fits')[0]),'w')
    
    if(exp):
      #input is scalelength directly
      sigmax=expscale
      sigmay=expscale
    else:
      #go from fwhm to sigma for Gaussian
      sigmax=spatwidth/(2*np.sqrt(2*np.log(2)))
      sigmay=spatwidth/(2*np.sqrt(2*np.log(2)))
    
    ind = 0
   
    #loop on mocks
    print('Injecting mock sources')
    
    mockimage = np.zeros_like(newimage)
    
    while ind<num:
        
        mflux=np.random.uniform(fluxlimits[0],fluxlimits[1])
        if ZP != -1:
          mflux=10**(-0.4*(mflux-ZP))
        
        xc=np.random.uniform(minx,maxx)
        yc=np.random.uniform(miny,maxy)
        
        sizex = int(np.ceil(3*sigmax))
        sizey = int(np.ceil(3*sigmay)) 
        
        #Verify availablity in seg map
        thisseg = np.sum(segima[int(np.ceil(yc))-sizey:int(np.ceil(yc))+sizey,int(np.ceil(xc))-sizex:int(np.ceil(xc))+sizex])
        thisima = np.sum(newimage[int(np.ceil(yc))-sizey:int(np.ceil(yc))+sizey,int(np.ceil(xc))-sizex:int(np.ceil(xc))+sizex])
        thisbad = np.sum(badmask[int(np.ceil(yc))-sizey:int(np.ceil(yc))+sizey,int(np.ceil(xc))-sizex:int(np.ceil(xc))+sizex])
        
        if thisseg > 0 or thisbad>0 or np.isnan(thisima):
           continue
        else:
          segima[int(np.ceil(yc))-sizey:int(np.ceil(yc))+sizey,int(np.ceil(xc))-sizex:int(np.ceil(xc))+sizex] = 1 
          ind += 1  
        
        norm=mflux/(sigmax*sigmay*(2*np.pi))

        #now evaluate model (Gaussian or Exponential) at pixel center 
        allx=np.round(np.arange(np.floor(xc-fill*sigmax),xc+fill*sigmax,1))
        ally=np.round(np.arange(np.floor(yc-fill*sigmay),yc+fill*sigmay,1))

        fl.write("{} {} {} \n".format(xc,yc,mflux))

        for xx in allx:
            for yy in ally:
                     if (xx>0) & (yy>0) & (xx<nx-1) & (yy<ny-1):
                        #evaluate model at pixel
                        if(exp):
                            rad= np.sqrt((xx-xc)**2+(yy-yc)**2)
                            pix=norm*np.exp(-1.*(( rad /sigmax )))
                        else:
                            pix=norm*np.exp(-1.*((((xx-xc)**2)/(2.*sigmax**2))+(((yy-yc)**2)/(2.*sigmay**2))))
                              
                        #store
                        mockimage[int(yy),int(xx)]=mockimage[int(yy),int(xx)]+pix
    
    
    if(exp):
      #The mock exponential profiles need to be convolved with a gaussian 2D filter to simulate PSF effects.
      #go from FWHM to sigma for the Kernel
      kern = Gaussian2DKernel(spatwidth/(2*np.sqrt(2*np.log(2))))
      mockimage = convolve(mockimage, kern)
      
    
    newimage += mockimage
    
    print('Done.. writing!')

    #store output
    if(ext==0):
        hdufirst = fits.PrimaryHDU(newimage)
        hdulist = fits.HDUList([hdufirst])
        hdulist[0].header=imahdu[0].header
    elif(ext==1):
        hdufirst = fits.PrimaryHDU([])
        hdusecond = fits.ImageHDU(newimage)
        hdulist = fits.HDUList([hdufirst,hdusecond])
        hdulist[0].header=imahdu[0].header
        hdulist[1].header=imahdu[1].header

    write='{}_{}'.format(outprefix,os.path.basename(image))
    hdulist.writeto(output+write,overwrite=True)
    
    hdusegout = fits.ImageHDU(segima)
    write='{}_{}'.format(outprefix,os.path.basename(segmap))
    hdusegout.writeto(output+write, overwrite=True)

    fl.close()
              
    return


def run_mockcont(image,  segmap, iters, outdir, outfile,  varima=None, badmask=None, expmap=None, bkgsub=False, magrange=[23,29], SNRdet=3., FWHM_pix=3., EXP_scale=1.3, exp=False, ZP=None, num=80, fill=10., minarea=10.,append=False):
   
   """

    Run several iterations of the mock injection and detection loop to build up a statistical
    sample of simulated objects.

    iters -> Number of injection/detection iterations to run
    outfile -> Name of the output file
    image -> a MUSE image [filename] to use for mock, must have ZPAB keyword in header
    segmap -> Segmentation map of the input image

    varima - > (Optional) variance image for detection thresholding     
    badmask -> (Optional) injection will not happen on masked pixels (1)
    expmap -> (Optional) if supplied, the final catalogue will report the value of the expmap at
              the position of the injected source
    magrange -> mock sources will be drawn in range [min,max] in AB mag units
    SNRdet -> SNR for detection
    FWHM_pix -> spatial FWHM for Gaussian model in pixel
    EXP_scale -> exponential scale lenght in pixels, this is convolved with the Gaussian FWHM (seeing)
    exp -> use an exponential profile in x,y. If false, use point sources
    num -> number of mock sources [high numbers give better statistics but can lead to shadowing]
    fill -> multiple of sigma to evaluate Gaussian. Larger number is more accurate but slower
    minarea -> minimum area for detection
    append -> If true append to the existing outfile instead of generating a new one.

   """


   from mypython.ifu import muse_source as source
   from astropy.io import fits, ascii
   from astropy.table import Table
   import numpy as np
   import os

   hduima = fits.open(image)
   if ZP is None:
     ZP = hduima[0].header['ZPAB']
   
   #Read expmap
   if expmap:
    hduexp = fits.open(expmap)
    expima = hduexp[0].data
   else:
    expima = np.ones_like(hduima[0].data)

   mockxc   = []
   mockyc   = []
   mockflux = []
   mockexp  = []
   mockmag  = []
   sexdet   = []
   sexxc    = []
   sexyc    = []
   sexflux  = []
   sexmag   = []
   
   if not append:
      with open(outdir+outfile, mode='w') as f:
        f.write('#mockxc mockyc mockexp mockflux mockmag sexdet sexxc sexyc sexflux sexmag \n')
   
   #Define more filenames
   cmocks_injcat = outdir+'cmocks_{}_catalogue.txt'.format(os.path.basename(image).split('.fits')[0])
   cmocks_image  = outdir+'cmocks_{}'.format(os.path.basename(image))
   cmocks_segmap = outdir+'cmocks_{}'.format(os.path.basename(segmap))
      
   for repeat in range(iters):
         
         print("Iteration {}".format(repeat))
         
         #Inject sources and read master catalogue
         mockcont(image, segmap, magrange, output=outdir, badmask=badmask, num=num, ZP=ZP, spatwidth=FWHM_pix, fill=fill, exp=exp, expscale=EXP_scale)
         tmpmock = ascii.read(cmocks_injcat, names=['Xpos', 'Ypos', 'Flux'])
         
         #Run sextractor
         source.findsources(cmocks_image, image, nsig=SNRdet, output=outdir, fitsmask=badmask, varima=varima, minarea=minarea, bkgsub=bkgsub)
         sexcat = fits.open(outdir+'catalogue.fits')[1].data
         
         for mockob in range(len(tmpmock)):
           
           thisxc = tmpmock['Xpos'][mockob]
           thisyc = tmpmock['Ypos'][mockob]
           
           mockxc.append(thisxc)
           mockyc.append(thisyc)
           mockflux.append(tmpmock['Flux'][mockob])
           mockmag.append(-2.5*np.log10(tmpmock['Flux'][mockob])+ZP)
           mockexp.append(expima[int(thisyc), int(thisxc)])
           
           distarray = np.sqrt((sexcat.x-thisxc)**2+(sexcat.y-thisyc)**2)
                   
           if np.min(distarray)<1.0:
              indmatch = np.argmin(distarray)
              
              sexdet.append(1)
              sexxc.append(sexcat.x[indmatch])
              sexyc.append(sexcat.y[indmatch])
              sexflux.append(sexcat.flux[indmatch])
              sexmag.append(-2.5*np.log10(sexcat.flux[indmatch])+ZP)   

           else:
              sexdet.append(0)
              sexxc.append(-1)
              sexyc.append(-1)
              sexflux.append(-1)  
              sexmag.append(-1) 
           
         if (repeat%5 ==0) and (repeat>0):
              data = Table([mockxc, mockyc, mockexp, mockflux, mockmag, sexdet, sexxc, sexyc, sexflux, sexmag], names=[\
              'mockxc','mockyc','mockexp','mockflux','mockmag','sexdet','sexxc','sexyc','sexflux','sexmag'])
   
              with open(outdir+outfile, mode='a') as f:
                   data.write(f, format='ascii.no_header')
              
              mockxc   = []
              mockyc   = []
              mockflux = []
              mockexp  = []
              mockmag  = []
              sexdet   = []
              sexxc    = []
              sexyc    = []
              sexflux  = []
              sexmag   = [] 
   
   os.remove(outdir+'catalogue.fits')
   os.remove(cmocks_image)
   os.remove(cmocks_injcat)
   os.remove(cmocks_segmap)


def run_mocklines(cube, varcube, segmap, cubexfile, iters, outdir, outfile, real=None, extended=False, badmask=None, expmap=None, cov_poly=None, tmpdir='tmpdir', FWHM_pix=3.0, SNRdet=5.0, fluxrange=[10,3000], num=1000):

    from mypython.ifu import muse_emitters as emi
    from astropy.io import fits, ascii
    from astropy.table import Table
    import subprocess
    import numpy as np
    import os, shutil

    """

    Run several iterations of the mock injection and detection loop to build up a statistical
    sample of simulated line emitters.

    cube -> a MUSE cube [filename] to use for emitters injection
    varcube -> variance image for detection thresholding
    segmap -> Segmentation map of the input image
    cubexfile -> file with cubex run parameters
    iters -> Number of injection/detection iterations to run
    outdir -> Name of the output directory
    outfile -> Name of the output filename
    
    real -> (Optional) path to a single or a list of small cubes with real sources. Theflux in these
            cubes must be normalized to 1. If set extended, FWHM_pix are disregarded
    badmask -> (Optional) injection will not happen on masked pixels (1)
    expmap -> (Optional) if supplied, the final catalogue will report the value of the expmap at
              the position of the injected source
    extended -> use an exponential profile in x,y with scalelength equal to the value. If false, use point sources.
                If the value is a list of values they will be randomly drawn in the mock run
    fluxrange -> mock sources will be drawn in range [min,max] in 1E-20 erg/s/cm2
    num -> number of mock sources [default=1000, high numbers give better statistics but can lead to shadowing]
    cov_poly -> covariance array. If 1D assumed to be cov vs aperture size and valid at all waves. 
                If 2D, the polynomial fit at the closest wave is used.
    SNRdet -> SNR for detection (default = 5). Covariance is applied on top of this.
    FWHM_pix -> spatial FWHM for Gaussian model in pixel
    
    """
    
    mockxc   = []
    mockyc   = []
    mockwc   = []
    mockflux = []
    mockexp  = []
    mocksig  = []
    sexdet   = []
    sexxc    = []
    sexyc    = []
    sexwc    = []
    sexflux1 = []
    sexflux2 = []
    
    folder = '/'+tmpdir+'/'
    
    if expmap is not None:
      #Read expmap
      hduexp = fits.open(expmap)
      expima = hduexp[0].data

    if not os.path.isdir(outdir+folder):
       os.makedirs(outdir+folder)
            
    for repeat in range(iters):
        print("####################################")
        print("Run " + str(repeat+1) + " of " + str(iters))

        mocklines(cube,segmap,fluxrange,real=real,badmask=badmask,output=outdir+folder,num=num, spatwidth=FWHM_pix, wavewidth=2, outprefix='lmocks', exp=extended)
        
        # run on mock cube 
        cube_mock = outdir+folder+'lmocks_'+os.path.basename(cube)
        cubexcat = outdir+folder+'lmocks_catalogue.txt'
        
        #extract sources
        print("Running CubEx")
        subprocess.call(["CubEx",cubexfile,"-Catalogue",cubexcat,"-cube",cube_mock, "-var", varcube,"-SourceMask", segmap, "-NCheckCubes", "0", "-isn", "{}".format(SNRdet)],stdout=open(os.devnull, 'w'))

        #Load extraction catalog
        f = open(cubexcat, 'r')
        o = open(cubexcat.replace('.txt', '_cln.txt'), 'w')
        for line in f:
           if '**' not in line and 'nan' not in line and 'NaN' not in line:
              o.write(line)
        f.close()
        o.close()
               
        extracted = emi.read_cubex_catalog(cubexcat.replace('.txt', '_cln.txt'))
        
        # read in catalogue of known mock sources
        tmpmock = ascii.read(outdir+folder+"{}_catalogue.txt".format(os.path.basename(cube_mock).split('.fits')[0]))
        
        for mockob in range(len(tmpmock)):
        
          thisxc = tmpmock['col1'][mockob]
          thisyc = tmpmock['col2'][mockob]
          thiswc = tmpmock['col3'][mockob]
          thisflux = tmpmock['col4'][mockob]
          thissigx = tmpmock['col5'][mockob]
          thissigy = tmpmock['col6'][mockob]
          
          mockxc.append(thisxc)
          mockyc.append(thisyc)
          mockwc.append(thiswc)
          mockflux.append(thisflux)
          mocksig.append(thissigx)
          
          if expmap is not None:
            mockexp.append(expima[int(thisyc), int(thisxc)])
          else:
            mockexp.append(1)

          sexx = np.array(extracted['x_fluxw'])
          sexy = np.array(extracted['y_fluxw'])
          sexw = np.array(extracted['z_fluxw'])
          sexflux_aper = np.array(extracted['flux_aper'])
          sexflux_iso  = np.array(extracted['IsoFlux'])
          sexerr_iso   = np.array(extracted['IsoErr'])
          sexsize      = np.sqrt(extracted['area_isoproj'])*0.2
          
          distarray = np.sqrt((sexx-thisxc)**2+(sexy-thisyc)**2+(sexw-thiswc)**2)
          
          tmpindmatch = np.nanargmin(distarray)
            
          if np.nanmin(distarray)<3.0:
             indmatch = tmpindmatch  
          elif np.nanmin(distarray)<9.0:
              if (sexflux_iso[tmpindmatch] > 0.2*thisflux) and (sexflux_iso[tmpindmatch] < 5.0*thisflux):
                indmatch = tmpindmatch
              else:
                indmatch = -1  
          else:
             indmatch = -1
          
          if indmatch>=0:
                 
             #Calculate SN including covariance
             if cov_poly is None:
                 covariance = 1.0
             elif cov_poly.ndim == 1 :
                 size = sexsize[indmatch]
                 covariance = np.polyval(cov_poly,size)
             elif cov_poly.ndim == 2 :
                 size = sexsize[indmatch]
                 try:
                    okind = np.where(extracted['lambda_geow'][indmatch]>cov_poly[:,0])[0][-1]
                 except:
                    okind = 0 
                 covariance = np.polyval(cov_poly[okind,2:],size)
             
             
             sexdet.append((sexflux_iso[indmatch]/sexerr_iso[indmatch])/covariance)
             
             sexxc.append(sexx[indmatch])
             sexyc.append(sexy[indmatch])
             sexwc.append(sexw[indmatch])
             sexflux1.append(sexflux_iso[indmatch])
             sexflux2.append(sexflux_aper[indmatch])

          else:
             sexdet.append(0)
             sexxc.append(-1)
             sexyc.append(-1)
             sexwc.append(-1) 
             sexflux1.append(-1) 
             sexflux2.append(-1) 
                   
        if (repeat%5 ==0):
             
             print('Output mode')
             
             if expmap is not None:
               data = Table([mockxc, mockyc, mockwc, mockexp, mockflux, sexdet, sexxc, sexyc, sexwc, sexflux1])
               header = '#mockxc mockyc mockwc mockexp mockflux sexdet sexxc sexyc sexwc sexfluxiso \n'
             else:
               data = Table([mockxc, mockyc, mockwc, mockflux, sexdet, sexxc, sexyc, sexwc, sexflux1])
               header = '#mockxc mockyc mockwc mockflux sexdet sexxc sexyc sexwc sexfluxiso \n'
             
             if not os.path.isfile(outdir+outfile):
                 with open(outdir+outfile, mode='w') as f:
                   f.write(header)
             
             with open(outdir+outfile, mode='a') as f:
                 data.write(f, format='ascii.no_header')
             
             mockxc   = []
             mockyc   = []
             mockwc   = []
             mockflux = []
             mockexp  = []
             mocksig  = []
             sexdet   = []
             sexxc    = []
             sexyc    = []
             sexwc    = []
             sexflux1 = []
             sexflux2 = []
                 
    shutil.rmtree(outdir+folder)


def run_mocklines_shine(cube, varcube, segmap, iters, outdir, outfile, fftconv=False, real=None, extended=False, badmask=None, expmap=None, cov_poly=None, tmpdir='tmpdir', FWHM_pix=3.0, SNRdet=5.0, fluxrange=[10,3000], num=1000):

    from mypython.ifu import muse_emitters as emi
    from astropy.io import fits, ascii
    from astropy.table import Table
    import subprocess
    import numpy as np
    import os, shutil
    
    from SHINE.shine import SHINE

    """

    Run several iterations of the mock injection and detection loop to build up a statistical
    sample of simulated line emitters. 
    
    IMPORTANT: THIS VERSION USES THE SHINE EXTRACTION CODE.
    Spectral Highlighting and Identification of Emission (SHINE) is publicy available on 
    Zenodo: doi: 10.5281/zenodo.14710518
    
    Please install the code from github.com/matteofox/SHINE.git

    
    cube      -> a MUSE cube [filename] to use for emitters injection
    varcube   -> variance image for detection thresholding
    segmap    -> Segmentation map of the input image
    iters     -> Number of injection/detection iterations to run
    outdir    -> Name of the output directory
    outfile   -> Name of the output filename
    
    fftconv   -> (Optional) if True, use fft convolution rather than the direct algorithm.
    real      -> (Optional) path to a single or a list of small cubes with real sources. Theflux in these
                 cubes must be normalized to 1. If set extended, FWHM_pix are disregarded
    badmask   -> (Optional) injection will not happen on masked pixels (1)
    expmap    -> (Optional) if supplied, the final catalogue will report the value of the expmap at
                 the position of the injected source
    extended  -> use an exponential profile in x,y with scalelength equal to the value. If false, use point sources.
                 If the value is a list of values they will be randomly drawn in the mock run
    fluxrange -> mock sources will be drawn in range [min,max] in 1E-20 erg/s/cm2
    num       -> number of mock sources [default=1000, high numbers give better statistics but can lead to shadowing]
    cov_poly  -> covariance array. If 1D assumed to be cov vs aperture size and valid at all waves. 
                 If 2D, the polynomial fit at the closest wave is used.
    SNRdet    -> SNR for detection (default = 5). Covariance is applied on top of this.
    FWHM_pix  -> spatial FWHM for Gaussian model in pixel
    
    """
    
    mockxc   = []
    mockyc   = []
    mockwc   = []
    mockflux = []
    mockexp  = []
    mocksig  = []
    sexdet   = []
    sexxc    = []
    sexyc    = []
    sexwc    = []
    sexflux1 = []
    
    folder = '/'+tmpdir+'/'
    
    if expmap is not None:
      #Read expmap
      hduexp = fits.open(expmap)
      expima = hduexp[0].data

    if not os.path.isdir(outdir+folder):
       os.makedirs(outdir+folder)
            
    for repeat in range(iters):
        print("####################################")
        print("Run " + str(repeat+1) + " of " + str(iters))

        mocklines(cube,segmap,fluxrange,real=real,badmask=badmask,output=outdir+folder,num=num, spatwidth=FWHM_pix, wavewidth=2, outprefix='lmocks', exp=extended)
        
        # run on mock cube 
        cube_mock = outdir+folder+'lmocks_'+os.path.basename(cube)
        shinecat  = outdir+folder+'lmocks_catalogue.txt'

        
        #-------------------------------------------
        #extract sources
        print("Running SHINE")
        SHINE.runextraction(cube_mock, varcube, mask2d=segmap, mask2dpost=segmap, snthreshold=3, maskspedge=0, extvardata=1, spatsmooth=2, usefftconv=fftconv, connectivity=26, mindz=3, maxdz=50, minvox = 27, minarea=9, outdir=outdir+folder, writelabels=False, writesmdata=False, writesmvar=False, writesmsnr=False, writesubcube=False, writevardata=False)

        extracted = Table.read(outdir + folder + os.path.basename(cube_mock).split('.fits')[0] + '.CATALOGUE_out.fits', format='fits')

        # read in catalogue of known mock sources
        tmpmock = ascii.read(outdir+folder+"{}_catalogue.txt".format(os.path.basename(cube_mock).split('.fits')[0]))
        
        
        for mockob in range(len(tmpmock)):
        
          thisxc = tmpmock['col1'][mockob]
          thisyc = tmpmock['col2'][mockob]
          thiswc = tmpmock['col3'][mockob]
          thisflux = tmpmock['col4'][mockob]
          thissigx = tmpmock['col5'][mockob]
          thissigy = tmpmock['col6'][mockob]
          
          mockxc.append(thisxc)
          mockyc.append(thisyc)
          mockwc.append(thiswc)
          mockflux.append(thisflux)
          mocksig.append(thissigx)
          
          if expmap is not None:
            mockexp.append(expima[int(thisyc), int(thisxc)])
          else:
            mockexp.append(1)

          sexx = np.array(extracted['Xcent'])
          sexy = np.array(extracted['Ycent'])
          sexw = np.array(extracted['Zcent'])
          sexflux_iso  = np.array(extracted['Flux'])
          sexerr_iso   = np.array(extracted['Flux_err'])
          sexsize      = np.sqrt(extracted['Nspat'])*0.2
          
          distarray = np.sqrt((sexx-thisxc)**2+(sexy-thisyc)**2+(sexw-thiswc)**2)
          
          tmpindmatch = np.nanargmin(distarray)
            
          if np.nanmin(distarray)<3.0:
             indmatch = tmpindmatch  
          elif np.nanmin(distarray)<9.0:
              if (sexflux_iso[tmpindmatch] > 0.2*thisflux) and (sexflux_iso[tmpindmatch] < 5.0*thisflux):
                indmatch = tmpindmatch
              else:
                indmatch = -1  
          else:
             indmatch = -1
          
          if indmatch>=0:
                 
             #Calculate SN including covariance
             if cov_poly is None:
                 covariance = 1.0
             elif cov_poly.ndim == 1 :
                 size = sexsize[indmatch]
                 covariance = np.polyval(cov_poly,size)
             elif cov_poly.ndim == 2 :
                 size = sexsize[indmatch]
                 try:
                    okind = np.where(extracted['lambda_geow'][indmatch]>cov_poly[:,0])[0][-1]
                 except:
                    okind = 0 
                 covariance = np.polyval(cov_poly[okind,2:],size)
             
             
             sexdet.append((sexflux_iso[indmatch]/sexerr_iso[indmatch])/covariance)
             
             sexxc.append(sexx[indmatch])
             sexyc.append(sexy[indmatch])
             sexwc.append(sexw[indmatch])
             sexflux1.append(sexflux_iso[indmatch])

          else:
             sexdet.append(0)
             sexxc.append(-1)
             sexyc.append(-1)
             sexwc.append(-1) 
             sexflux1.append(-1) 
    
                   
        if (repeat%5 ==0):
             
             print('Output mode')
             
             if expmap is not None:
               data = Table([mockxc, mockyc, mockwc, mockexp, mockflux, sexdet, sexxc, sexyc, sexwc, sexflux1])
               header = '#mockxc mockyc mockwc mockexp mockflux sexdet sexxc sexyc sexwc sexfluxiso \n'
             else:
               data = Table([mockxc, mockyc, mockwc, mockflux, sexdet, sexxc, sexyc, sexwc, sexflux1])
               header = '#mockxc mockyc mockwc mockflux sexdet sexxc sexyc sexwc sexfluxiso \n'
             
             if not os.path.isfile(outdir+outfile):
                 with open(outdir+outfile, mode='w') as f:
                   f.write(header)
             
             with open(outdir+outfile, mode='a') as f:
                 data.write(f, format='ascii.no_header')
             
             mockxc   = []
             mockyc   = []
             mockwc   = []
             mockflux = []
             mockexp  = []
             mocksig  = []
             sexdet   = []
             sexxc    = []
             sexyc    = []
             sexwc    = []
             sexflux1 = []
                 
    shutil.rmtree(outdir+folder)



def analyse_mockcont(infile, exprange=None, nbins=250, complfrac=0.9):
  
  from astropy.io import fits, ascii
  import numpy as np
  import datetime
  
  """

    Analyse the mock table and return the detection fraction histogram in 
    nbins.

    infile -> Mock table
    exprange -> (Optional) Limit the analysis to the sources injected where
    the expmap is (min,max) exposures
    nbins -> Number of bins for the histogram
    complfrac -> fraction of detected sources that defines the magnitude limit
    
  """
  
  #Format fast_basic + supplying column names is 2x faster but the first line of valid
  #data is lost
  data = ascii.read(infile, format='fast_commented_header')
  
  mockxc   = data['mockxc']
  mockyc   = data['mockyc']   
  mockflux = data['mockflux']
  mockexp  = data['mockexp'] 
  sexdet   = data['sexdet']   
  sexxc    = data['sexxc']     
  sexyc    = data['sexyc']     
  sexflux  = data['sexflux']  
  mockmag  = data['mockmag']
  sexmag   = data['sexmag'] 
  
  detected = (sexdet==1)

  if exprange==None:
     exprange = (np.min(mockexp),np.max(mockexp))
  
  okexp = (mockexp>exprange[0]) & (mockexp<=exprange[1])
 
  bin_edges_mag = np.linspace(np.min(mockmag),np.max(mockmag), nbins)
  bin_cent_mag  = (bin_edges_mag[1:]+bin_edges_mag[:-1])/2.
  
  hist_all_mag, bin_dummy = np.histogram(mockmag[okexp], bins=bin_edges_mag)
  hist_det_mag, bin_dummy = np.histogram(mockmag[(detected) & (okexp)], bins=bin_edges_mag)

  det_frac_mag = 1.*hist_det_mag/hist_all_mag

  #Find first time from the left that the poly goes above 0.9
  maglim = bin_cent_mag[np.argmin(det_frac_mag>complfrac)]
  
  return det_frac_mag, bin_edges_mag, maglim
  

    
