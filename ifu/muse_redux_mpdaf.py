"""

These are sets of procedures optimised for almost empty fields using mpdaf procedures


"""
from __future__ import print_function
import matplotlib.pyplot as plt
from mypython.fits import pyregmask as msk
    
def coaddcubes(listob,nclip=2.5):

    """

    Loop over each OB and make final combined (mean, median) cubes
    
    listob -> OBs to process
    nclip -> threshold for sigmaclipping

 
    """

    
    import os
    import glob
    import subprocess
    import shutil
    from astropy.io import fits
    import muse_utils as mut 
    import numpy as np
    import matplotlib.pyplot as plt

    
    #names for median/mean output
    cubemed="mpdafcombine/COMBINED_CUBE_MED_FINAL.fits"
    imagemed="mpdafcombine/COMBINED_IMAGE_MED_FINAL.fits"
    cubemean="mpdafcombine/COMBINED_CUBE_FINAL.fits"
    imagemean="mpdafcombine/COMBINED_IMAGE_FINAL.fits"

    if ((not os.path.isfile(cubemed)) | (not os.path.isfile(cubemean))) :
        print('Compute mean/median coadd {} {}'.format(cubemed,cubemean))
          
        #first collect all relevant exposures
        allexposures=[]
        for ob in listob:
            finalexp=glob.glob("{}/Proc/MPDAF/DATACUBE_RESAMPLED_EXP*_zap.fits".format(ob))
            allexposures+=finalexp

        #loop over exposures for  combine
        for i,exp in enumerate(allexposures):
            print('Ingesting data {}'.format(exp))
            #open exposure
            data=fits.open(exp)
            
            #operations to do only once 
            if(i == 0):
                nw,nx,ny=data[1].data.shape
                nexp=len(allexposures)
                datashape=(nexp,nw,nx,ny)
                alldata=np.ma.zeros(datashape)
                alldata.mask=np.full(datashape,False)
                headermain=data[0].header
                headerext=data[1].header
                headervar=data[2].header

            #store data
            alldata.data[i]=data[1].data

            #catch nans in this exposure
            alldata.mask[i]=np.logical_not(np.isfinite(data[1].data))

            data.close()

            #now handle mask for this exposure if available in cubex or mpdaf region
            mpdafmask=exp.split("DATACUBE")[0]+'IMAGE'+exp.split("DATACUBE")[1]
            mpdafreg=mpdafmask.split('fits')[0]+'reg'
            cubexmask=exp.split("_RESAMPLED")
            cubexmask=cubexmask[0]+"_FINAL_RESAMPLED"+cubexmask[1]
            cubexmask=cubexmask.replace('MPDAF','Cubex')
            cubexmask1=cubexmask.replace('_zap.fits','_fixhsn_SliceEdgeMask_wreg.fits')
            cubexmask2=cubexmask.replace('_zap.fits','_fixhsn_SliceEdgeMask.fits')
            
            #handle mpdaf first 
            if os.path.isfile(mpdafreg):
                #now fill using region 
                Mask = msk.PyMask(ny,nx,mpdafreg)
                for ii in range(Mask.nreg):
                    Mask.fillmask(ii)
                    if(ii == 0):
                        slicemask=Mask.mask
                    else:
                        slicemask+=Mask.mask
                slicecube=np.tile(slicemask,(nw,1,1))
                alldata.mask[i]=slicecube
                cubexmask=None
            #else, look at cubex and wor as selector
            elif os.path.isfile(cubexmask1):
                cubexmask=cubexmask1
            elif os.path.isfile(cubexmask2):
                cubexmask=cubexmask2
            else:
                cubexmask=None

            #fill in cubex mask if available 
            if(cubexmask):
                #mask edges
                fitsmask=fits.open(cubexmask)
                slicemask=np.full(fitsmask[0].data.shape,False)
                masked=np.where(fitsmask[0].data < 1)
                slicemask[masked]=True
                #grow cube
                slicecube=np.tile(slicemask,(nw,1,1))
                alldata.mask[i]=slicecube
                fitsmask.close()
                
        #sanitise output
        alldata=np.ma.masked_invalid(alldata)

        #iterative sigmaclipping
        for niter in range(3):
            #now catch significant outliers
            print('Catching outliers with sigma clipping; loop {}'.format(niter+1))
            stdevtmp=alldata.std(axis=0)
            medtmp=np.ma.median(alldata,axis=0)
    
            #loop over exposures
            for i,exp in enumerate(allexposures):
                print('Masking outliers in {}'.format(exp))
                alldata[i]=np.ma.masked_where(abs(alldata[i] - medtmp) >= nclip*stdevtmp,alldata[i])
                            

        #save exposure time map
        print('Compute exposure map')
        expcube=alldata.count(axis=0)
        expmap=np.sum(expcube,axis=0)/nw
      
        #at last perform median combine
        print('Computing median')
        medcube=np.ma.median(alldata,axis=0)
        
        #at last perform mean combine
        print('Computing mean')
        meancube=alldata.mean(axis=0)

        #now ingest variance
        for i,exp in enumerate(allexposures):
            print('Ingesting variance {}'.format(exp))
            #open exposure
            data=fits.open(exp)
            #store data
            alldata.data[i]=data[2].data
            data.close()
            
        print('Computing variance')
        varcube=alldata.sum(axis=0)/expcube/expcube

        #save exposure time map
        hdu1 = fits.PrimaryHDU([])
        hdu2 = fits.ImageHDU(expmap)
        hdu2.header=headerext
        hdulist = fits.HDUList([hdu1,hdu2])
        hdulist.writeto("mpdafcombine/FINAL_COADD_EXPOSUREMAP.fits",overwrite=True)
 

        #save output median
        hdu = fits.PrimaryHDU([])
        hdu1= fits.ImageHDU(medcube.data)
        hdu2= fits.ImageHDU(varcube.data)
        hdulist = fits.HDUList([hdu,hdu1,hdu2])
        hdulist[0].header=headermain
        hdulist[1].header=headerext
        hdulist[2].header=headervar
        hdulist.writeto(cubemed,overwrite=True)
        
        #now make white image
        print ('Creating final white image from median cube')
        white_new=np.zeros((nx,ny))
        for xx in range(nx):
            for yy in range(ny):
                #skip couple of channells for cosmetics 
                white_new[xx,yy]=np.nansum(medcube.data[2:-3,xx,yy])/(nw-4)  
                
        #save projected image 
        hdu1 = fits.PrimaryHDU([])
        hdu2 = fits.ImageHDU(white_new)
        hdu2.header=headerext
        hdulist = fits.HDUList([hdu1,hdu2])
        hdulist.writeto(imagemed,overwrite=True)
        
        #save output mean
        hdu = fits.PrimaryHDU([])
        hdu1= fits.ImageHDU(meancube.data)
        hdu2= fits.ImageHDU(varcube.data)
        hdulist = fits.HDUList([hdu,hdu1,hdu2])
        hdulist[0].header=headermain
        hdulist[1].header=headerext
        hdulist[2].header=headervar
        hdulist.writeto(cubemean,overwrite=True)

        #now make white image
        print ('Creating final white image from mean cube')
        white_new=np.zeros((nx,ny))
        for xx in range(nx):
            for yy in range(ny):
                #skip couple of channells for cosmetics 
                white_new[xx,yy]=np.nansum(meancube.data[2:-3,xx,yy])/(nw-4)  
                
        #save projected image 
        hdu1 = fits.PrimaryHDU([])
        hdu2 = fits.ImageHDU(white_new)
        hdu2.header=headerext
        hdulist = fits.HDUList([hdu1,hdu2])
        hdulist.writeto(imagemean,overwrite=True)
    
    else:
        
        print('Median and mean coadd {} already exists!'.format(cubemed,cubemean))



def zapskysub(listob, skymask=None):
    
    """

    Loop over each OB and performs zap sky subtraction
    

    listob -> OBs to process
 
    """

    import os
    import glob
    import subprocess
    import shutil
    from astropy.io import fits
    import muse_utils as mut 
    import numpy as np
    import sep
    import zap
    from mypython.fits import pyregmask as pmk

    #grab top dir
    topdir=os.getcwd()

    #now loop over each folder and make the final sky-subtracted cubes
    for ob in listob:
        
        #change dir
        os.chdir(ob+'/Proc/MPDAF')

        #use source mask already available 
        srcmask='selfcalib_mask.fits'
        
	if(skymask):
            skyfitsmask = 'skysrc_mask.fits'
	    
	    srchdu = fits.open(srcmask)
	    srcmsk = srchdu[0].data
	    nx=srchdu[0].header['NAXIS1']
            ny=srchdu[0].header['NAXIS2']
	    
	    #construct the sky region mask
            mysky=pmk.PyMask(nx,ny,"../../"+skymask,header=srchdu[0].header)
            for ii in range(mysky.nreg):
                mysky.fillmask(ii)
                srcmask=srcmsk+mysky.mask
	    
	    outhdu = fits.PrimaryHDU(srcmsk)
	    outhdu.writeto(skyfitsmask)
	    
	    srchdu.close()

	else:
	    skyfitsmask=srcmask    
	    	
        #now loop over exposures and apply self calibration
        scils=glob.glob("../Basic/OBJECT_RED_0*.fits*")
        nsci=len(scils)

        print("Processing {} with ZAP".format(ob))
        
        #loop on exposures and apply self calibration
        for exp in range(nsci):
             
            #these are the self-calibrated data
            cubeselfcal="DATACUBE_RESAMPLED_EXP{0:d}_fix.fits".format(exp+1)
            imageselfcal="IMAGE_RESAMPLED_EXP{0:d}_fix.fits".format(exp+1)

            #these are the sky subtracted products
            cubezap="DATACUBE_RESAMPLED_EXP{0:d}_zap.fits".format(exp+1)
            imagezap="IMAGE_RESAMPLED_EXP{0:d}_zap.fits".format(exp+1)

            if not os.path.isfile(cubezap):
                print('Reconstruct cube {} with ZAP'.format(cubezap))
                zap.process(cubeselfcal,outcubefits=cubezap,clean=True,mask=skyfitsmask)

                #create white image from zap cube
                cube=fits.open(cubezap)
                #define geometry 
                nwave=cube[1].header["NAXIS3"]
                nx=cube[1].header["NAXIS1"]
                ny=cube[1].header["NAXIS2"]
   
                print ('Creating final white image from ZAP')
                white_new=np.zeros((ny,nx))
                for xx in range(nx):
                    for yy in range(ny):
                        white_new[yy,xx]=np.nansum(cube[1].data[:,yy,xx])/nwave  
                    
                #save projected image 
                hdu1 = fits.PrimaryHDU([])
                hdu2 = fits.ImageHDU(white_new)
                hdu2.header=cube[1].header
                hdulist = fits.HDUList([hdu1,hdu2])
                hdulist.writeto(imagezap,overwrite=True)

            else:
                print('Cube {} already exists! Skip...'.format(cubezap))
            
        #back to top
        os.chdir(topdir)


def selfcalibrate(listob,deepwhite,refpath='esocombine',nproc=24):
    
    """

    Loop over each OB and performs self-calibration after masking sources as 

    listob -> OBs to process
    deepwhite -> the best white image available to mask sources
    refpath -> where to find a white image to be used as reference wcs grid
    nproc -> number of processors

    """
    import os
    import glob
    import subprocess
    import shutil
    from astropy.io import fits
    import muse_utils as mut 
    import numpy as np
    import sep
    from mpdaf.drs import PixTable

    #grab top dir
    topdir=os.getcwd()

    #now loop over each folder and make the final sky-subtracted cubes
    for ob in listob:
        
        #change dir
        os.chdir(ob+'/Proc/MPDAF')


        print("Processing {} for self-calibration".format(ob))

        #make source mask 
        srcmask='selfcalib_mask.fits'
        if(os.path.isfile(srcmask)):
            print("Source mask already exists!")
        else:
            print("Create source mask for selfcalibration")
  
            #open the ifu mask to create a good mask 
            deepdata=fits.open("../../../"+deepwhite)

            ##define geometry 
            #nx=deepdata[0].header["NAXIS1"]
            #ny=deepdata[0].header["NAXIS2"]

            #now flag the sources (allow cubex/eso images)
            try:
                image=deepdata[0].data.byteswap().newbyteorder()
                datheader=deepdata[0].header
            except:
                image=deepdata[1].data.byteswap().newbyteorder()
                datheader=deepdata[1].header

            bkg=sep.Background(image)
            bkg.subfrom(image)
            obj,segmap=sep.extract(image,5.*bkg.globalrms,minarea=6,segmentation_map=True)
            segmap[np.where(segmap > 0)]=1
            
            #write source mask to disk 
            hdu=fits.PrimaryHDU(segmap,header=datheader)
            hdu.writeto(srcmask,overwrite=True)


        #now loop over exposures and apply self calibration
        scils=glob.glob("../Basic/OBJECT_RED_0*.fits*")
        nsci=len(scils)
        
        #loop on exposures and apply self calibration
        for exp in range(nsci):
            
            pixselfcal="PIXTABLE_REDUCED_RESAMPLED_EXP{0:d}_fix.fits".format(exp+1)
            if not os.path.isfile(pixselfcal):
                print("Apply self-calibration to {}".format(pixselfcal))
                
                #open reduced pixel table
                pix=PixTable("PIXTABLE_REDUCED_RESAMPLED_EXP{0:d}.fits".format(exp+1))
                
                #create mask
                maskpix=pix.mask_column(srcmask)
                maskpix.write("PIXTABLE_REDUCED_RESAMPLED_EXP{0:d}_mask.fits".format(exp+1))
                #selfcalibrate
                autocalib=pix.selfcalibrate(pixmask=maskpix)
              
                #write to disk
                autocalib.write("PIXTABLE_REDUCED_RESAMPLED_EXP{0:d}_autocalib.fits".format(exp+1))
                pix.write(pixselfcal)
                                    
            else:
                print("Self-calibration for {} already done! Skip...".format(pixselfcal))

            cubeselfcal="DATACUBE_RESAMPLED_EXP{0:d}_fix.fits".format(exp+1)
            imageselfcal="IMAGE_RESAMPLED_EXP{0:d}_fix.fits".format(exp+1)
            if not os.path.isfile(cubeselfcal):
                print('Reconstruct cube {}'.format(cubeselfcal))
     
                #now reconstruct cube and white image on right reference frame 
                #handle sof file
                sof_name="../../Script/scipost_mpdaf_self{0:d}.sof".format(exp+1)
                sofedit=open(sof_name,'w')
                sofedit.write('../../../{}/DATACUBE_FINAL.fits OUTPUT_WCS\n'.format(refpath))
                sofedit.write("PIXTABLE_REDUCED_RESAMPLED_EXP{0:d}_fix.fits PIXTABLE_OBJECT\n".format(exp+1))
                sofedit.close()
                
                #now run script
                scr=open("../../Script/make_scipost_mpdaf_self{0:d}.sh".format(exp+1),"w")
                scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
                scr.write('esorex --log-file=scipost_mpdaf_self{0:d}.log muse_scipost_make_cube  ../../Script/scipost_mpdaf_self{0:d}.sof'.format(exp+1))
                scr.close()
                
                #Run pipeline 
                subprocess.call(["sh", "../../Script/make_scipost_mpdaf_self{0:d}.sh".format(exp+1)])    
                subprocess.call(["mv","DATACUBE_FINAL.fits",cubeselfcal])
                subprocess.call(["mv","IMAGE_FOV_0001.fits",imageselfcal])

            else:
                print('Cube {} already exists! Skip...'.format(cubeselfcal))

            
        #back to top
        os.chdir(topdir)


def individual_resample(listob,refpath='./',nproc=24):

    """
    Loop over each OB and re-run scipost using a final coadded cube as
    a reference for WCS. This produces calibrated pixeltable that are all comaptible with the
    astrometric solution of the reference cube.

    If desired, one can then reconstruct a cube by running 
    > esorex  muse_scipost_make_cube makecube.sof
    with "PIXTABLE_REDUCED_RESAMPLED_EXP1.fits PIXTABLE_OBJECT" 
         "../../.././esocombine/DATACUBE_FINAL.fits OUTPUT_WCS"
    in sof file 

    listob -> OBs to process
    refpath -> where reference path is for WCS resampling
    nproc -> numer of processors in parallel runs 

    """
      
    import os
    import glob
    import subprocess
    import shutil
    from astropy.io import fits
    import muse_utils as mut 
    import numpy as np

    #grab top dir
    topdir=os.getcwd()

    #now loop over each folder and make the final sky-subtracted cubes
    for ob in listob:
        
        #change dir
        os.chdir(ob+'/Proc/')

        #make cubex folder
        if not os.path.exists('MPDAF'):
            os.makedirs('MPDAF')

        #change dir
        os.chdir('MPDAF')

        print('Processing {} for resampling on reference cube'.format(ob))
 
        #Search how many exposures are there
        scils=glob.glob("../Basic/OBJECT_RED_0*.fits*")
        nsci=len(scils)
        
        #loop on exposures and reduce frame with sky subtraction 
        for exp in range(nsci):
            
            #define some output names for final cube 
            pname="PIXTABLE_REDUCED_RESAMPLED_EXP{0:d}.fits".format(exp+1)
 
            if not os.path.isfile(pname):
                print("Processing exposure {0:d} to align to reference".format(exp+1))
                
                #copy sof file written for basic reduction
                sof_old=open("../../Script/scipost_{0:d}.sof".format(exp+1))
                sof_name="../../Script/scipost_mpdaf_{0:d}.sof".format(exp+1)
                sofedit=open(sof_name,'w')
                

                #now apply offsets to pixel table
                print ('Apply offsets...')
                pixtablist=[]
                for ll in sof_old:
                   if('STD_' in ll or 'PIXTABLE_OBJECT' in ll):
                        fil,tag=ll.split(' ')
                        sofedit.write("../Basic/"+fil+" "+tag)
                   else:
                        sofedit.write(ll)
		
		#Check existence of ABSOLUTE offset list otherwise fall back onto the relative one
		if os.path.isfile('../../../{}/OFFSET_LIST_ABS.fits'.format(refpath)):
		   sofedit.write('../../../{}/OFFSET_LIST_ABS.fits OFFSET_LIST\n'.format(refpath))
		else:
		   sofedit.write('../../../{}/OFFSET_LIST.fits OFFSET_LIST\n'.format(refpath))   

                #append reference frame to sof file 
                sofedit.write('../../../{}/DATACUBE_FINAL.fits OUTPUT_WCS\n'.format(refpath))
                sofedit.close()
                sof_old.close()
                
                #Write the command file 
                scr=open("../../Script/make_scipost_mpdaf_{0:d}.sh".format(exp+1),"w")
                scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
                
                scr.write('esorex --log-file=scipost_mpdaf_{0:d}.log muse_scipost --skymethod="none" --save=individual ../../Script/scipost_mpdaf_{0:d}.sof'.format(exp+1))
                scr.close()
                
                #Run pipeline 
                subprocess.call(["sh", "../../Script/make_scipost_mpdaf_{0:d}.sh".format(exp+1)])    
                subprocess.call(["mv","PIXTABLE_REDUCED_0001.fits",pname])
            else:
                print("Exposure {0:d} exists.. skip! ".format(exp+1))
     

        #back to top
        os.chdir(topdir)


def old_individual_resample(listob,refpath='./',nproc=24):

    """
    Loop over each OB and re-run scipost using a final coadded cube as
    a reference for WCS. This produces calibrated pixeltable that are all comaptible with the
    astrometric solution of the reference cube.

    If desired, one can then reconstruct a cube by running 
    > esorex  muse_scipost_make_cube makecube.sof
    with "PIXTABLE_REDUCED_RESAMPLED_EXP1.fits PIXTABLE_OBJECT" 
         "../../.././esocombine/DATACUBE_FINAL.fits OUTPUT_WCS"
    in sof file 


    listob -> OBs to process
    refpath -> where reference path is for WCS resampling
    nproc -> numer of processors in parallel runs 

    """
      
    import os
    import glob
    import subprocess
    import shutil
    from astropy.io import fits
    import muse_utils as mut 
    import numpy as np

    #grab top dir
    topdir=os.getcwd()

    #now loop over each folder and make the final sky-subtracted cubes
    for ob in listob:
        
        #change dir
        os.chdir(ob+'/Proc/')

        #make cubex folder
        if not os.path.exists('MPDAF'):
            os.makedirs('MPDAF')

        #change dir
        os.chdir('MPDAF')

        print('Processing {} for resampling on reference cube'.format(ob))
 
        #Search how many exposures are there
        scils=glob.glob("../Basic/OBJECT_RED_0*.fits*")
        nsci=len(scils)
        
        #loop on exposures and reduce frame with sky subtraction 
        for exp in range(nsci):
            
            if not os.path.isfile('OFFSET_LIST_EXP{0:d}.fits'.format(exp+1)):
                print("Compute offsets...")
                            
                #create align file 
                alignsof=open('../../Script/align_toref_{0:d}.sof'.format(exp+1),'w')
                alignsof.write("../../../{}/IMAGE_FOV_0001.fits IMAGE_FOV\n".format(refpath))
                alignsof.write("../Basic/IMAGE_FOV_EXP{0:d}.fits IMAGE_FOV\n".format(exp+1))
                alignsof.close()
                
                #run script align with respect to registered reference cube  
                alignscr=open('../../Script/make_align_toref_{0:d}.sh'.format(exp+1),'w')
                alignscr.write("esorex --log-file=align_toref_{0:d}.log muse_exp_align --threshold=4. ../../Script/align_toref_{0:d}.sof".format(exp+1))
                alignscr.close()
                subprocess.call(["sh","../../Script/make_align_toref_{0:d}.sh".format(exp+1)])    
                
                #copy the offsets 
                alig=fits.open('OFFSET_LIST.fits')
                alig.writeto('OFFSET_LIST_EXP{0:d}.fits'.format(exp+1),clobber=True)

            else:
                print('Offsets exist.. skip')

            #define some output names for final cube 
            pname="PIXTABLE_REDUCED_RESAMPLED_EXP{0:d}.fits".format(exp+1)
 
            if not os.path.isfile(pname):
                print("Processing exposure {0:d} to align to reference".format(exp+1))
                
                #copy sof file written for basic reduction
                sof_old=open("../../Script/scipost_{0:d}.sof".format(exp+1))
                sof_name="../../Script/scipost_mpdaf_{0:d}.sof".format(exp+1)
                sofedit=open(sof_name,'w')
                
                #read the offsets 
                alig=fits.open('OFFSET_LIST_EXP{0:d}.fits'.format(exp+1))
                offsets=alig[1].data[1]

                #now apply offsets to pixel table
                print ('Apply offsets...')
                pixtablist=[]
                for ll in sof_old:
                    if('PIXTABLE_OBJECT' in ll):
                        pixtab=ll.split(' ')[0]
                        pxt=fits.open('../Basic/'+pixtab)
                        pxt[0].header['RA']=pxt[0].header['RA']-offsets[2]
                        pxt[0].header['DEC']=pxt[0].header['DEC']-offsets[3]
                        
                        pxt.writeto("WCS_"+pixtab,clobber=True)
                        pixtablist.append("WCS_"+pixtab)
                        sofedit.write("WCS_"+pixtab+" PIXTABLE_OBJECT\n")
                    elif('STD_' in ll):
                        fil,tag=ll.split(' ')
                        sofedit.write("../Basic/"+fil+" "+tag)
                    else:
                        sofedit.write(ll)

                #append reference frame to sof file 
                sofedit.write('../../../{}/DATACUBE_FINAL.fits OUTPUT_WCS\n'.format(refpath))
                sofedit.close()
                sof_old.close()
            
                
                #Write the command file 
                scr=open("../../Script/make_scipost_mpdaf_{0:d}.sh".format(exp+1),"w")
                scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
                
                scr.write('esorex --log-file=scipost_mpdaf_{0:d}.log muse_scipost --skymethod="none" --save=individual ../../Script/scipost_mpdaf_{0:d}.sof'.format(exp+1))
                scr.close()
                
                #Run pipeline 
                subprocess.call(["sh", "../../Script/make_scipost_mpdaf_{0:d}.sh".format(exp+1)])    
                subprocess.call(["mv","PIXTABLE_REDUCED_0001.fits",pname])
            else:
                print("Exposure {0:d} exists.. skip! ".format(exp+1))
     

        #clean dir for unwanted stuff...
        print ('Clean directory!')
        garbage=glob.glob("WCS_PIXTABLE_OBJECT*")
        for gg in garbage:
            os.remove(gg)

        #back to top
        os.chdir(topdir)
