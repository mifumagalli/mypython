"""

These are sets of procedures optimised for almost empty fields using mpdaf procedures


"""
from __future__ import print_function


def coaddcubes(listob):

    """

    Loop over each OB and make final combined (mean, median) cubes
    
    listob -> OBs to process
 
    """

    
    import os
    import glob
    import subprocess
    import shutil
    from astropy.io import fits
    import muse_utils as mut 
    import numpy as np
    

    #first collect all relevant exposures
    allexposures=[]
    for ob in listob:
        finalexp=glob.glob("{}/Proc/MPDAF/DATACUBE*ZAP*fits".format(ob))
        allexposures.append(allexposures)

        


def zapskysub(listob):
    
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
    

    #grab top dir
    topdir=os.getcwd()

    #now loop over each folder and make the final sky-subtracted cubes
    for ob in listob:
        
        #change dir
        os.chdir(ob+'/Proc/MPDAF')

        #use source mask already available 
        srcmask='selfcalib_mask.fits'

        #now loop over exposures and apply self calibration
        scils=glob.glob("../Basic/OBJECT_RED_0*.fits*")
        nsci=len(scils)
        
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
                zap.process(cubeselfcal,outcubefits=cubezap,clean=True,mask=srcmask)

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

            #now flag the sources
            image=deepdata[0].data.byteswap().newbyteorder()
            bkg=sep.Background(image)
            bkg.subfrom(image)
            obj,segmap=sep.extract(image,5.*bkg.globalrms,minarea=6,segmentation_map=True)
            segmap[np.where(segmap > 0)]=1

            #write source mask to disk 
            hdu=fits.PrimaryHDU(segmap,header=deepdata[0].header)
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
