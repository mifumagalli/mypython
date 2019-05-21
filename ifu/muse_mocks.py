"""
These are sets of utilities to handle muse sources

"""


def mocklines(cube,segmap,fluxlimits,badmask=None,output='./',num=500,wavelimits=None,spatwidth=3.5,wavewidth=2,outprefix='mocks',fill=6., exp=False):

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
    exp -> use an exponential profile in x,y. If false, use point sources

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
    
    segcube = np.repeat(segima[np.newaxis,...], np.shape(newcube)[0], axis=0)
    
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
    
    minx=20
    maxx=cubhdu[0].header['NAXIS1']-20
    miny=20
    maxy=cubhdu[0].header['NAXIS2']-20

    #open output to store mocks
    fl=open(output+"{}_{}_catalogue.txt".format(outprefix,os.path.basename(cube).split('.fits')[0]),'w')
    
    #compute Gaussian parameters 
    if exp:
       sigmax=exp
       sigmay=exp
    else:
       sigmax=spatwidth/(2*np.sqrt(2*np.log(2)))
       sigmay=spatwidth/(2*np.sqrt(2*np.log(2)))
    
    sigmaw=wavewidth/(2*np.sqrt(2*np.log(2)))
    
    #loop on mocks
    print('Injecting mock sources')
    
    mockcube = np.zeros_like(newcube)
    
    ind = 0
    
    while ind<num:

      #now draw random distributions
      mflux    = 10**np.random.uniform(np.log10(fluxlimits[0]),np.log10(fluxlimits[1]))
      xc       = np.random.uniform(minx,maxx)
      yc       = np.random.uniform(miny,maxy)
      wc       = np.random.uniform(minw,maxw)

      sizex = int(np.ceil(3*sigmax))
      sizey = int(np.ceil(3*sigmay))
      sizew = int(np.ceil(3*sigmaw)) 

      #Verify availablity in seg map
      thisseg = np.sum(segcube[int(np.ceil(wc))-sizew:int(np.ceil(wc))+sizew,int(np.ceil(yc))-sizey:int(np.ceil(yc))+sizey,int(np.ceil(xc))-sizex:int(np.ceil(xc))+sizex])
      thiscub = np.sum(newcube[int(np.ceil(wc))-sizew:int(np.ceil(wc))+sizew,int(np.ceil(yc))-sizey:int(np.ceil(yc))+sizey,int(np.ceil(xc))-sizex:int(np.ceil(xc))+sizex])
      thisbad = np.sum(badmask[int(np.ceil(yc))-sizey:int(np.ceil(yc))+sizey,int(np.ceil(xc))-sizex:int(np.ceil(xc))+sizex])

      if thisseg > 0 or thisbad>0 or np.isnan(thiscub):
         continue
      else:
        segcube[int(np.ceil(wc))-sizew:int(np.ceil(wc))+sizew,int(np.ceil(yc))-sizey:int(np.ceil(yc))+sizey,int(np.ceil(xc))-sizex:int(np.ceil(xc))+sizex] = 1 
        ind += 1  

      norm=mflux/(sigmax*sigmay*sigmaw*(2*np.pi)**(3./2.))
      
      #now evaluate Gaussian at pixel center [can do better with Gauss integral]
      
      allx=np.round(np.arange(np.floor(xc-fill*sigmax),xc+fill*sigmax,1))
      ally=np.round(np.arange(np.floor(yc-fill*sigmay),yc+fill*sigmay,1))
      allw=np.round(np.arange(np.floor(wc-fill*sigmaw),wc+fill*sigmaw,1))
    
      fl.write("{} {} {} {}\n".format(xc,yc,wc,mflux))
      
      for xx in allx:
              for yy in ally:
        	  for ww in allw:
    
        	      if((xx > minx) & (yy > miny) & (ww > minw) & 
        		 (xx < maxx) & (yy < maxy) & (ww < maxw)):
        		  #evaluate 
    
        		  if exp:
        		      pix= norm*np.exp(-1.*(( abs((xx-xc)) /sigmax )+
        					    ( abs((yy-yc)) /sigmay )+
        					    (((ww-wc)**2)/(2.*sigmaw**2))))
        		  
        		  else:
        		      pix=norm*np.exp(-1.*((((xx-xc)**2)/(2.*sigmax**2))+
        				       (((yy-yc)**2)/(2.*sigmay**2))+
        				       (((ww-wc)**2)/(2.*sigmaw**2))))      
        	       
        		  #store
        		  mockcube[int(ww),int(yy),int(xx)]=mockcube[int(ww),int(yy),int(xx)]+pix
    
    if(exp):
      #The mock exponential profiles need to be convolved with a gaussian 2D filter to simulate PSF effects.
      #go from FWHM to sigma for the Kernel
      kern = kern = (0,spatwidth/(2*np.sqrt(2*np.log(2))),spatwidth/(2*np.sqrt(2*np.log(2))))
      mockcube =  filters.gaussian_filter(mockcube, kern, order=0)    
    
    newcube += mockcube

    
    print('Done.. writing!')
    

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


def mockcont(image,segmap,fluxlimits,badmask=None,num=100,ZP=-1,spatwidth=3.5,outprefix='cmocks',fill=6.,exp=False, expscale=1.5):

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
    
    minx=20
    maxx=imahdu[0].header['NAXIS1']-20
    miny=20
    maxy=imahdu[0].header['NAXIS2']-20
    
    #now draw random distributions

    #open output to store mocks
    fl=open("{}_{}_catalogue.txt".format(outprefix,os.path.basename(image).split('.fits')[0]),'w')
    
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
                    if((xx > minx) & (yy > miny) & (xx < maxx) & (yy < maxy)):
                        #evaluate model at pixel
			if(exp):
			    pix=norm*np.exp(-1.*(( abs((xx-xc)) /sigmax )+( abs((yy-yc)) /sigmay )))
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
    hdulist.writeto(write,overwrite=True)
    
    hdusegout = fits.ImageHDU(segima)
    write='{}_{}'.format(outprefix,os.path.basename(segmap))
    hdusegout.writeto(write, overwrite=True)

    fl.close()
              
    return


def run_mockcont(iters, outfile, image, varima, segmap, badmask=None, expmap=None, magrange=[23,29], SNRdet=3., FWHM_pix=3., EXP_scale=1.3, exp=False, num=80, fill=10., minarea=10.,overwrite=False):
   
   """

    Run several iterations of the mock injection and detection loop to build up a statistical
    sample of simulated objects.

    iters -> Number of injection/detection iterations to run
    outfile -> Name of the output file
    image -> a MUSE image [filename] to use for mock, must have ZPAB keyword in header
    varima - > variance image for detection thresholding
    segmap -> Segmentation map of the input image
    
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
    overwrite -> If true append to the existing outfile instead of generating a new one.

   """


   from mypython.ifu import muse_source as source
   from astropy.io import fits, ascii
   from astropy.table import Table
   import numpy as np
   import os

   hduima = fits.open(image)
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
   
   if not overwrite:
      with open(outfile, mode='w') as f:
        f.write('#mockxc mockyc mockexp mockflux mockmag sexdet sexxc sexyc sexflux sexmag \n')
   
   #Define more filenames
   cmocks_injcat = 'cmocks_{}_catalogue.txt'.format(os.path.basename(image).split('.fits')[0])
   cmocks_image  = 'cmocks_{}'.format(os.path.basename(image))
   cmocks_segmap = 'cmocks_{}'.format(os.path.basename(segmap))
      
   for repeat in range(iters):
   	 
   	 print("Iteration {}".format(repeat))
   	 
   	 #Inject sources and read master catalogue
   	 mockcont(image, segmap, magrange, badmask=badmask, num=num, ZP=ZP, spatwidth=FWHM_pix, fill=fill, exp=exp, expscale=EXP_scale)
   	 tmpmock = ascii.read(cmocks_injcat, names=['Xpos', 'Ypos', 'Flux'])
   	 
   	 #Run sextractor
   	 source.findsources(cmocks_image, image, nsig=SNRdet, fitsmask=badmask, varima=varima, minarea=minarea)
   	 sexcat = fits.open('catalogue.fits')[1].data
   	 
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
   	   
	 if (repeat%10 ==0) and (repeat>0):
	      data = Table([mockxc, mockyc, mockexp, mockflux, mockmag, sexdet, sexxc, sexyc, sexflux, sexmag], names=[\
              'mockxc','mockyc','mockexp','mockflux','mockmag','sexdet','sexxc','sexyc','sexflux','sexmag'])
   
              with open(outfile, mode='a') as f:
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
 
   os.remove('catalogue.fits')
   os.remove(cmocks_image)
   os.remove(cmocks_injcat)
   os.remove(cmocks_segmap)


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
  

    
