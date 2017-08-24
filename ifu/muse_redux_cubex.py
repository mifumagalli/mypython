
def fixandsky_firstpass(cube,pixtab,noclobber,skymask=None):
    
    
    """ 
    Take a cube and pixel table and fix the cube, then skysub and produce white image,
    using CubEx utils

    is skymask is set, mask regions when computing normalisation 

    """

    import os 
    import subprocess
    import numpy as np
    from astropy.io import fits
    from mypython.fits import pyregmask as pmk
    import matplotlib.pyplot as plt
               
    #make some intermediate names
    fixed=cube.split('.fits')[0]+"_fix.fits"
    skysub=cube.split('.fits')[0]+"_skysub.fits"
    white=cube.split('.fits')[0]+"_white.fits"
    sharpmsk=cube.split('.fits')[0]+"_sharpmask.fits"

    #if told to mask sky do it.. otherwise leave image empty
    cb=fits.open(cube)
    nx=cb[1].header['NAXIS1']
    ny=cb[1].header['NAXIS2']
    sharpmask=np.zeros((ny,nx))
    if(skymask):
        #construct the sky region mask
        mysky=pmk.PyMask(nx,ny,"../../"+skymask,header=cb[1].header)
        for ii in range(mysky.nreg):
            mysky.fillmask(ii)
            sharpmask=sharpmask+mysky.mask
    #write mask
    hdu = fits.PrimaryHDU(sharpmask)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(sharpmsk,clobber=True)
    cb.close()

    #now fix the cube
    if ((os.path.isfile(fixed)) & (noclobber)):
        print "Cube {0} already fixed".format(cube)
    else:
        print 'Cubefix ', cube
        subprocess.call(["CubeFix","-cube", cube,"-pixtable", pixtab,"-out", fixed])

    #now run cube skysub
    if ((os.path.isfile(skysub)) & (noclobber)):
        print "Cube {0} already skysub".format(fixed)
    else:
        print 'Sky sub ', fixed
        subprocess.call(["CubeSharp","-cube",fixed,"-out",skysub,"-sourcemask",sharpmsk,"-lcheck",".false."])
                               
    #create a white image
    if ((os.path.isfile(white)) & (noclobber)):
        print "White image for cube {0} already exists".format(skysub)
    else:
        print 'Create white image for ', skysub
        subprocess.call(["Cube2Im","-cube",skysub,"-out",white])
                

def fixandsky_secondpass(cube,pixtab,noclobber,highsn=None,skymask=None):
        
    """ 
 
    Similar to first pass, but operates on cubes that have been realigned and uses masks 
    as appropriate

    If highsn cube is provided, use it to mask the sources and for skysubtraction
 
    """
      
    import os 
    import subprocess
    import numpy as np
    from astropy.io import fits
    from mypython.fits import pyregmask as pmk

    if(highsn):
        #prepare final names
        fixed=cube.split('.fits')[0]+"_fixhsn.fits"
        skysub=cube.split('.fits')[0]+"_skysubhsn.fits"
        white=cube.split('.fits')[0]+"_whitehsn.fits"
        sharpmsk=cube.split('.fits')[0]+"_sharpmasksn.fits"
    else:
        #prepare intermediate names
        fixed=cube.split('.fits')[0]+"_fix2.fits"
        skysub=cube.split('.fits')[0]+"_skysub2.fits"
        white=cube.split('.fits')[0]+"_white2.fits"
        sharpmsk=cube.split('.fits')[0]+"_sharpmask2.fits"

    #assign names for source mask 
    mask_source=cube.split('.fits')[0]+"_white.Objects_Id.fits"
    white_source=cube.split('.fits')[0]+"_white.fits"

    #now fix the cube using masks
    if ((os.path.isfile(fixed)) & (noclobber)):
        print "Cube {0} already fixed".format(cube)
    else:

        print 'Create source mask ', white_source
        #if high cube provide, overwrite white image 
        if(highsn):
            print 'Using high SN cube...'
            subprocess.call(["Cube2Im","-cube",highsn,"-out",white_source])
            subprocess.call(["CubEx",white_source,'-MultiExt','.false.','-SN_Threshold','3','-RescaleVar','.true.'])
        else:
            print 'Using white image from previous loop'
            #create source mask 
            subprocess.call(["CubEx",white_source,'-MultiExt','.false.','-SN_Threshold','5','-RescaleVar','.true.'])
            
        print 'Cubefix ', cube
        subprocess.call(["CubeFix","-cube", cube,"-pixtable", pixtab,"-out", fixed,"-sourcemask",mask_source]) 

        #At this step, check out cubeAdd2Mask if want to fix edges or weird ifus/slices 

    #if told to mask sky do it.. otherwise leave image empty
    cb=fits.open(cube)
    nx=cb[1].header['NAXIS1']
    ny=cb[1].header['NAXIS2']
    sharpmask=np.zeros((ny,nx))
    if(skymask):
        #construct the sky region mask
        mysky=pmk.PyMask(nx,ny,"../../"+skymask,header=cb[1].header)
        for ii in range(mysky.nreg):
            mysky.fillmask(ii)
            sharpmask=sharpmask+mysky.mask
    cb.close()

    #inject src mask in sharpmask
    srcmk=fits.open(mask_source)
    sharpmask=sharpmask+srcmk[0].data

    #write mask
    hdu = fits.PrimaryHDU(sharpmask)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(sharpmsk,clobber=True)

    #now run cube skysub
    if ((os.path.isfile(skysub)) & (noclobber)):
        print "Cube {0} already skysub".format(fixed)
    else:
        print 'Sky sub ', fixed
        if(highsn):
            #now few more options to control sky sub 
            subprocess.call(["CubeSharp","-cube",fixed,"-out",skysub,"-sourcemask",sharpmsk,"-hsncube",highsn,"-lcheck",".false."])
        else:
            subprocess.call(["CubeSharp","-cube",fixed,"-out",skysub,"-sourcemask",sharpmsk])
                               
    #create a white image
    if ((os.path.isfile(white)) & (noclobber)):
        print "White image for cube {0} already exists".format(skysub)
    else:
        print 'Create white image for ', skysub
        subprocess.call(["Cube2Im","-cube",skysub,"-out",white])
                
def cubex_driver(listob,last=False,highsn=None,skymask=None):
    
    """
    Procedures that drives the loops of cubex within each OB folder

    listob -> the list of OBs to process
    last  -> set to True for final pass with high-sn cube 
    highsn -> name of the highsn cube used for masking 
    skymask -> mask this region in source mask before running cubesharp

    """
    
    import os
    import glob
    import subprocess
    import multiprocessing
    import numpy as np

    #grab top dir
    topdir=os.getcwd()

    #now loop over each folder and make the final sky-subtracted cubes
    for ob in listob:
        
        #change dir
        os.chdir(ob+'/Proc/')
        print('Processing {} with cubex '.format(ob))

        #Search how many exposures are there
        scils=glob.glob("OBJECT_RED_0*.fits*")
        nsci=len(scils)
        
        #this is the final pass with highsn cube
        if(last):
            #############################################################
            # do a final loop of cubex fix and sky with proper masking  #
            #############################################################
            print ('Final pass of cubex')
            workers=[]
            for dd in range(nsci):
                #reconstruct the name 
                pixtab="PIXTABLE_REDUCED_LINEWCS_EXP{0:d}.fits".format(dd+1)
                cube="DATACUBE_FINAL_LINEWCS_EXP{0:d}.fits".format(dd+1)
                #now launch the task
                p = multiprocessing.Process(target=fixandsky_secondpass,args=(cube,pixtab,True,highsn,skymask))
                workers.append(p)
                p.start()
   
            #wait for completion of all of them 
            for p in workers:
                if(p.is_alive()):
                    p.join()  
            
        else:
            
            ########################################
            # do a first loop of cubex fix and sky #
            ########################################
            #do it in parallel on exposures
            print ('First pass of cubex')
            workers=[]
            for dd in range(nsci):
                #reconstruct the name 
                pixtab="PIXTABLE_REDUCED_LINEWCS_EXP{0:d}.fits".format(dd+1)
                cube="DATACUBE_FINAL_LINEWCS_EXP{0:d}.fits".format(dd+1)
                #now launch the task
                p = multiprocessing.Process(target=fixandsky_firstpass,args=(cube,pixtab,True,skymask))
                workers.append(p)
                p.start()
   
            #wait for completion of all of them 
            for p in workers:
                if(p.is_alive()):
                    p.join()

                    
            #############################################################
            # do a second loop of cubex fix and sky with proper masking #
            #############################################################
            print ('Second pass of cubex')
            workers=[]
            for dd in range(nsci):
                #reconstruct the name 
                pixtab="PIXTABLE_REDUCED_LINEWCS_EXP{0:d}.fits".format(dd+1)
                cube="DATACUBE_FINAL_LINEWCS_EXP{0:d}.fits".format(dd+1)
                #now launch the task
                p = multiprocessing.Process(target=fixandsky_secondpass,args=(cube,pixtab,True,None,skymask))
                workers.append(p)
                p.start()
   
            #wait for completion of all of them 
            for p in workers:
                if(p.is_alive()):
                    p.join()  
                
        #back to top
        os.chdir(topdir)

def combine_cubes(cubes,masks,regions=True,final=False,halfset=False,halfsetfinal=False):

    """
    Combine a bunch of cubes using masks with CubeCombine
        
    cubes    -> a list of cubes to use in the combine
    masks    -> a list of goodpix masks from the pipeline
    regions  -> if True, code searches for ds9 region files inside path with same name as 
                pipeline mask (.reg), to mask additional area that one wants to clip
    final    -> is True, append final tag to name and prepare median cubes
    
    halfset       -> if set to tag name, append/uses suffix for coadding indepenent halfs 
    halfsetfinal  -> if set to tag name in final loop, append/uses suffix for coadding indepenent halfs 

    """
    import subprocess
    import os
    import numpy as np
    from astropy.io import fits
    from mypython.fits import pyregmask as msk
    
    #define some names for the cubes
    if(final):
        cname="COMBINED_CUBE_FINAL.fits"
        iname="COMBINED_IMAGE_FINAL.fits"
        cmed="COMBINED_CUBE_MED_FINAL.fits"
        imed="COMBINED_IMAGE_MED_FINAL.fits"
        scriptname='runcombine_final.sh'
    elif(halfset):
        cname="COMBINED_CUBE_{}.fits".format(halfset)
        iname="COMBINED_IMAGE_{}.fits".format(halfset)
        cmed="COMBINED_CUBE_MED_{}.fits".format(halfset)
        imed="COMBINED_IMAGE_MED_{}.fits".format(halfset)
        scriptname='runcombine_{}.sh'.format(halfset)
    elif(halfsetfinal):
        cname="COMBINED_CUBE_FINAL_{}.fits".format(halfsetfinal)
        iname="COMBINED_IMAGE_FINAL_{}.fits".format(halfsetfinal)
        cmed="COMBINED_CUBE_MED_FINAL_{}.fits".format(halfsetfinal)
        imed="COMBINED_IMAGE_MED_FINAL_{}.fits".format(halfsetfinal)
        scriptname='runcombine_final_{}.sh'.format(halfsetfinal)
    else:
        cname="COMBINED_CUBE.fits"
        iname="COMBINED_IMAGE.fits"
        cmed="COMBINED_CUBE_MED.fits"
        imed="COMBINED_IMAGE_MED.fits"
        scriptname='runcombine.sh'


    if(os.path.isfile(cname)):
        print ('Cube {} already exists... skip!'.format(cname))
    else:
        print ('Creating combined cube {}'.format(cname))

        if(regions):
            print "Updating the masks"
        
            #loads list
            listmask=np.loadtxt(masks,dtype=np.dtype('a'))

            #redefine new mask 
            mask_new="new_"+masks
            llms=open(mask_new,"w")
            
             #check more than one mask                      Ruari 24/05/17 to deal with data with only 1 OB
            if listmask.shape == ():
                i=0
                cmask=np.array([listmask])[0]
               #create region name
                regname=(cmask.split(".fits")[0])+".reg"
                
                #search if file exist
                if(os.path.isfile(regname)):
                    
                    #update the mask 
                    print ("Updating mask for {}".format(regname))

                    #open fits
                    cfits=fits.open(cmask)
             
                    #init reg mask
                    Mask = msk.PyMask(cfits[0].header["NAXIS1"],cfits[0].header["NAXIS2"],regname)
                    for ii in range(Mask.nreg):
                        Mask.fillmask(ii)
                        if(ii == 0):
                            totmask=Mask.mask
                        else:
                            totmask+=Mask.mask
            
                    #update the mask
                    cfits[0].data=cfits[0].data*1*np.logical_not(totmask)   
                    savename=cmask.split(".fits")[0]+'_wreg.fits'
                    cfits.writeto(savename,clobber=True)
                    llms.write(savename+'\n')
                
                
                else:
                    #keep current mask 
                    llms.write(cmask+'\n')
            else:
                #loop over and update with regions
                for i,cmask in enumerate(listmask):
                
                    #create region name
                    regname=(cmask.split(".fits")[0])+".reg"
                    
                    #search if file exist
                    if(os.path.isfile(regname)):
                        
                        #update the mask 
                        print ("Updating mask for {}".format(regname))

                        #open fits
                        cfits=fits.open(cmask)
                 
                        #init reg mask
                        Mask = msk.PyMask(cfits[0].header["NAXIS1"],cfits[0].header["NAXIS2"],regname)
                        for ii in range(Mask.nreg):
                            Mask.fillmask(ii)
                            if(ii == 0):
                                totmask=Mask.mask
                            else:
                                totmask+=Mask.mask
                
                        #update the mask
                        cfits[0].data=cfits[0].data*1*np.logical_not(totmask)   
                        savename=cmask.split(".fits")[0]+'_wreg.fits'
                        cfits.writeto(savename,clobber=True)
                        llms.write(savename+'\n')
                    
                    
                    else:
                        #keep current mask 
                        llms.write(cmask+'\n')


            #done with new masks
            llms.close()

        else:
            print 'Using original masks'
            mask_new=masks

        #now run combine
        print 'Combine the cube...'    
      
        #make mean cube - write this as script that can be ran indepedently 
        scr=open(scriptname,'w')
        scr.write("export OMP_NUM_THREADS=1\n")
        scr.write("CubeCombine -list "+cubes+" -out "+cname+" -masklist "+mask_new+"\n")
        scr.write("Cube2Im -cube "+cname+" -out "+iname+"\n")
        scr.write("CubeCombine -list "+cubes+" -out "+cmed+" -masklist "+mask_new+" -comb median\n")
        scr.write("Cube2Im -cube "+cmed+" -out "+imed)
        scr.close()
        subprocess.call(["sh",scriptname])

def dataquality(cubeslist,maskslist):

    
    """
    Perform bunch of QA checks to asses the final data quality of reduction
    
    """
    
    import os
    import numpy as np
    from astropy.io import fits
    from mypython.fits import pyregmask as msk
    from mypython.ifu import muse_utils as mutil
    from mypython.ifu import muse_source as msrc
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from astropy.stats import sigma_clipped_stats
    
    try:
        from photutils import CircularAperture, aperture_photometry,\
            data_properties, properties_table, centroids
    except:
        print("To run checks need photutils package")
        return

    print("Perform QA checks...")

    #make QA folder
    if not os.path.exists('QA'):
        os.makedirs('QA')

    #cube names
    cname="COMBINED_CUBE.fits"
    iname="COMBINED_IMAGE.fits"
   
    #first identify bright sources in final white image 
    catsrc=msrc.findsources(iname,cname,check=True,output='QA',nsig=5,minarea=20)

    #make rsdss images
    if not(os.path.isfile('QA/sdssr.fits')):
        mutil.cube2img(cname,write='QA/sdssr.fits',wrange=None,helio=0,filt=129)
    rsdssall=fits.open('QA/sdssr.fits')
    segmask=fits.open('QA/segmap.fits')
    whiteref=fits.open(iname)

    #select round and bright objects 
    shapesrc=catsrc['a']/catsrc['b']
    roundsrc=catsrc[np.where((shapesrc < 1.1) & (catsrc['cflux'] > 50))]
    imgfield=fits.open(iname)
    rms=np.std(imgfield[0].data)

    #perform aperture photometry on rband image - data already skysub
    positions = [roundsrc['x'],roundsrc['y']]
    apertures = CircularAperture(positions, r=10.)
    phot_table = aperture_photometry(rsdssall[1].data, apertures)
    phot_table_white = aperture_photometry(whiteref[0].data, apertures)
    rmag=-2.5*np.log10(phot_table['aperture_sum'])+rsdssall[0].header['ZPAB']
    wmag=-2.5*np.log10(phot_table_white['aperture_sum'])


    #find FWHM on rband image
    fwhm=np.zeros(len(rmag))
    for ii in range(len(rmag)):
        subdata=rsdssall[1].data[roundsrc['y'][ii]-10:roundsrc['y'][ii]+10,\
                                     roundsrc['x'][ii]-10:roundsrc['x'][ii]+10]
        tdfit=centroids.fit_2dgaussian(subdata, error=None, mask=None)
        fwhm[ii]=2.3548*0.5*(tdfit.x_stddev+tdfit.y_stddev)*rsdssall[0].header['PC2_2']*3600.

    #find rms of cube - mask sources and add edge buffer
    maskwbuffer=np.copy(segmask[1].data)
    maskwbuffer[0:30,:]=9999
    maskwbuffer[-31:-1,:]=9999
    maskwbuffer[:,0:30]=9999
    maskwbuffer[:,-31:-1]=9999
    cwrms,crms=mutil.cubestat(cname,mask=maskwbuffer)

    #open diagnostic output
    with PdfPages('QA/QAfile.pdf') as pdf:

        ###########################
        #display field with r mag #
        ###########################

        plt.figure(figsize=(10,10))
        plt.imshow(imgfield[0].data,origin='low',clim=[-0.5*rms,0.5*rms],cmap='gray_r')
        #mark round soruces
        plt.scatter(roundsrc['x'],roundsrc['y'],color='red')
        for ii in range(len(rmag)):
            plt.text(roundsrc['x'][ii],roundsrc['y'][ii]," "+str(rmag[ii]),color='red')
        plt.title('Round sources with SDSS r mag')
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

        ###########################
        #display FWHM             #
        ###########################
        plt.figure(figsize=(10,10))
        plt.scatter(rmag,fwhm,color='red')
        plt.xlabel('Source rmag')
        plt.ylabel('FWHM (arcsec)')
        plt.title('Median FWHM {}'.format(np.median(fwhm)))
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

        ###########################
        #check centroid           #
        ###########################
        plt.figure(figsize=(10,10))
        
        #loop on exposures
        for tmpc in open(cubeslist,'r'):
            thisob=tmpc.split('/')[1]
            thisexp=tmpc.split('_')[3]
            wname='../{}/Proc/DATACUBE_FINAL_LINEWCS_{}_white2.fits'.format(thisob,thisexp)
            wfits=fits.open(wname)
            
            #now loop on sources
            delta_x=np.zeros(len(rmag))
            delta_y=np.zeros(len(rmag))
            for ii in range(len(rmag)):
                subdata=wfits[0].data[roundsrc['y'][ii]-10:roundsrc['y'][ii]+10,\
                                          roundsrc['x'][ii]-10:roundsrc['x'][ii]+10]
                x1,y1=centroids.centroid_2dg(subdata)
                delta_x[ii]=10.5-x1
                delta_y[ii]=10.5-y1
                
            #plot for this subunit
            plt.scatter(delta_x*rsdssall[0].header['PC2_2']*3600.,delta_y*rsdssall[0].header['PC2_2']*3600.)

        plt.xlabel('Delta x (arcsec)')
        plt.ylabel('Delta y (arcsec)')
        plt.title('Check exposure aligment')
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
        
        ###########################
        #check fluxing            #
        ###########################
   
        #make a check on fluxing 
        plt.figure(figsize=(10,10))
                   
        #loop on exposures
        for tmpc in open(cubeslist,'r'):
            thisob=tmpc.split('/')[1]
            thisexp=tmpc.split('_')[3]
            wname='../{}/Proc/DATACUBE_FINAL_LINEWCS_{}_white2.fits'.format(thisob,thisexp)
            wfits=fits.open(wname)
            
            phot_this_white = aperture_photometry(wfits[0].data, apertures)
            wmag_this=-2.5*np.log10(phot_this_white['aperture_sum'])

            #plot for this subunit
            ii=np.argsort(rmag)
            dd=wmag-wmag_this
            plt.plot(rmag[ii],dd[ii],label=thisob+thisexp)

        plt.xlabel('SDSS R mag')
        plt.ylabel('Delta White Mag')
        plt.title('Check exposure photometry')
        plt.legend()
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

        #display rms stats + compute stats over cubes
        plt.figure(figsize=(10,10))
        plt.semilogy(cwrms,crms,label='Coadd')
        
        for tmpc in open(cubeslist,'r'):
            cwrms_this,crms_this=mutil.cubestat(tmpc.strip(),mask=maskwbuffer)
            plt.semilogy(cwrms_this,crms_this,label=tmpc.split('/')[1]+tmpc.split('_')[3])

        plt.xlabel('Wave (A)')
        plt.ylabel('RMS (SB units)')
        plt.legend()
        plt.title('RMS in cubex')
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
        
        

