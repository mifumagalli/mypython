def individual_resample(listob,refpath='./',nproc=24):

    """
    Loop over each OB and re-run scipost using a final coadded cube as
    a reference for WCS. This produces cubes that are all regridded to 
    a common 3D grid with a single interpolation. 


    listob -> OBs to process
    refpath -> where reference path is for WCS resampling
    nproc -> numer of processors in parallel runs 

    """
      
    import os
    import glob
    import subprocess
    import shutil
    from astropy.io import fits
    from . import muse_utils as mut 
    import numpy as np

    #grab top dir
    topdir=os.getcwd()

    #now loop over each folder and make the final sky-subtracted cubes
    for ob in listob:
        
        #change dir
        os.chdir(ob+'/Proc/')

        #make cubex folder
        if not os.path.exists('Cubex'):
            os.makedirs('Cubex')

        #change dir
        os.chdir('Cubex')

        print('Processing {} for resampling on reference cube'.format(ob))
 
        #Search how many exposures are there
        scils=glob.glob("../Basic/OBJECT_RED_0*.fits*")
        nsci=len(scils)
        
        #loop on exposures and reduce frame with sky subtraction 
        for exp in range(nsci):
            

            #define some output names for final cube 
            cname="DATACUBE_FINAL_RESAMPLED_EXP{0:d}.fits".format(exp+1)
            pname="PIXTABLE_REDUCED_RESAMPLED_EXP{0:d}.fits".format(exp+1)
            iname="IMAGE_FOV_RESAMPLED_EXP{0:d}.fits".format(exp+1)
 
            if not os.path.isfile(cname):
                print("Processing exposure {0:d} to align to reference".format(exp+1))
                
                #copy sof file written for basic reduction
                sof_old=open("../../Script/scipost_{0:d}.sof".format(exp+1))
                sof_name="../../Script/scipost_line_{0:d}.sof".format(exp+1)
                sofedit=open(sof_name,'w')
                
                #now apply offsets to pixel table
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
                scr=open("../../Script/make_scipost_line_{0:d}.sh".format(exp+1),"w")
                scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
                
                scr.write('esorex --log-file=scipost_line_{0:d}.log muse_scipost --filter=white  --skymethod="none" --save=cube,individual ../../Script/scipost_line_{0:d}.sof'.format(exp+1))
                scr.close()
                
                #Run pipeline 
                subprocess.call(["sh", "../../Script/make_scipost_line_{0:d}.sh".format(exp+1)])    
                subprocess.call(["mv","DATACUBE_FINAL.fits",cname])
                subprocess.call(["mv","IMAGE_FOV_0001.fits",iname])
                subprocess.call(["mv","PIXTABLE_REDUCED_0001.fits",pname])
            else:
                print("Exposure {0:d} exists.. skip! ".format(exp+1))
     

        #back to top
        os.chdir(topdir)




def old_individual_resample(listob,refpath='./',nproc=24):

    """
    Loop over each OB and re-run scipost using a final coadded cube as
    a reference for WCS. This produces cubes that are all regridded to 
    a common 3D grid with a single interpolation. 


    listob -> OBs to process
    refpath -> where reference path is for WCS resampling
    nproc -> numer of processors in parallel runs 

    """
      
    import os
    import glob
    import subprocess
    import shutil
    from astropy.io import fits
    from . import muse_utils as mut 
    import numpy as np

    #grab top dir
    topdir=os.getcwd()

    #now loop over each folder and make the final sky-subtracted cubes
    for ob in listob:
        
        #change dir
        os.chdir(ob+'/Proc/')

        #make cubex folder
        if not os.path.exists('Cubex'):
            os.makedirs('Cubex')

        #change dir
        os.chdir('Cubex')

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
                alig.writeto('OFFSET_LIST_EXP{0:d}.fits'.format(exp+1),overwrite=True)
                


            else:
                print('Offsets exist.. skip')

            #define some output names for final cube 
            cname="DATACUBE_FINAL_RESAMPLED_EXP{0:d}.fits".format(exp+1)
            pname="PIXTABLE_REDUCED_RESAMPLED_EXP{0:d}.fits".format(exp+1)
            iname="IMAGE_FOV_RESAMPLED_EXP{0:d}.fits".format(exp+1)
 
            if not os.path.isfile(cname):
                print("Processing exposure {0:d} to align to reference".format(exp+1))
                
                #copy sof file written for basic reduction
                sof_old=open("../../Script/scipost_{0:d}.sof".format(exp+1))
                sof_name="../../Script/scipost_line_{0:d}.sof".format(exp+1)
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
                        
                        pxt.writeto("WCS_"+pixtab,overwrite=True)
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
                scr=open("../../Script/make_scipost_line_{0:d}.sh".format(exp+1),"w")
                scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
                
                scr.write('esorex --log-file=scipost_line_{0:d}.log muse_scipost --filter=white  --skymethod="none" --save=cube,individual ../../Script/scipost_line_{0:d}.sof'.format(exp+1))
                scr.close()
                
                #Run pipeline 
                subprocess.call(["sh", "../../Script/make_scipost_line_{0:d}.sh".format(exp+1)])    
                subprocess.call(["mv","DATACUBE_FINAL.fits",cname])
                subprocess.call(["mv","IMAGE_FOV_0001.fits",iname])
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



def fixandsky_firstpass(cube,pixtab,noclobber,skymask=None,version='1.6'):
    
    
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


    #define the step to catch error
    if('1.8' in version):
        errorstep='4'
    else:
        errorstep='3'
    

    #now fix the cube
    if ((os.path.isfile(fixed)) & (noclobber)):
        print("Cube {0} already fixed".format(cube))
    else:
        print('Cubefix ', cube)
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
        hdulist.writeto(sharpmsk,overwrite=True)
        cb.close()

        #run cubefix
        subprocess.call(["CubeFix","-cube", cube,"-pixtable", pixtab,"-out", fixed])
        
        #catch werid cases in which cubefix crashes
        if not os.path.isfile(fixed):
            print('Redo fix with last step only')
            #try again with last step only
            subprocess.call(["CubeFix","-cube", cube,"-pixtable", pixtab,"-out", fixed,"-step",errorstep])


    #now run cube skysub
    if ((os.path.isfile(skysub)) & (noclobber)):
        print("Cube {0} already skysub".format(fixed))
    else:
        print('Sky sub ', fixed)
        subprocess.call(["CubeSharp","-cube",fixed,"-out",skysub,"-sourcemask",sharpmsk,"-lcheck",".false."])
                               
    #create a white image
    if ((os.path.isfile(white)) & (noclobber)):
        print("White image for cube {0} already exists".format(skysub))
    else:
        print('Create white image for ', skysub)
        subprocess.call(["Cube2Im","-cube",skysub,"-out",white])
                

def fixandsky_secondpass(cube,pixtab,noclobber,highsn=None,skymask=None,exthsnmask=None,version='1.6'):
        
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


    #define the step to catch error
    if('1.8' in version):
        errorstep='4'
    else:
        errorstep='3'


    if(highsn):
        #prepare final names
        fixed=cube.split('.fits')[0]+"_fixhsn.fits"
        skysub=cube.split('.fits')[0]+"_skysubhsn.fits"
        white=cube.split('.fits')[0]+"_whitehsn.fits"
        sharpmsk=cube.split('.fits')[0]+"_sharpmaskhsn.fits"
    else:
        #prepare intermediate names
        fixed=cube.split('.fits')[0]+"_fix2.fits"
        skysub=cube.split('.fits')[0]+"_skysub2.fits"
        white=cube.split('.fits')[0]+"_white2.fits"
        sharpmsk=cube.split('.fits')[0]+"_sharpmask2.fits"

    #now fix the cube using masks
    if ((os.path.isfile(fixed)) & (noclobber)):
        print("Cube {0} already fixed".format(cube))
    else:
        print('Create source mask...')
        if(highsn):
            #The final mask will always have the same name and will overwrite existing file
            mask_source=cube.split('.fits')[0]+"_whitedeep.Objects_Id.fits"
            if exthsnmask:
              print('Using external HSN mask...')
              subprocess.call(["cp", exthsnmask, mask_source])
            else:
              print('Using high SN cube...')
              #create source mask from deep exposure 
              white_source=cube.split('.fits')[0]+"_whitedeep.fits"
              subprocess.call(["Cube2Im","-cube",highsn,"-out",white_source])
              subprocess.call(["CubEx",white_source,'-MultiExt','.false.','-ApplyFilter','.true.','-ApplyFilterVar','.true.','-FilterXYRad','1','-SN_Threshold','7','-MinNSpax','5'])
        else:
            print('Using white image from previous loop')
            #create source mask from previous step
            mask_source=cube.split('.fits')[0]+"_white.Objects_Id.fits"
            white_source=cube.split('.fits')[0]+"_white.fits"
            subprocess.call(["CubEx",white_source,'-MultiExt','.false.','-SN_Threshold','4.5','-RescaleVar','.true.'])
            
        print('Cubefix ', cube)
        subprocess.call(["CubeFix","-cube", cube,"-pixtable", pixtab,"-out", fixed,"-sourcemask",mask_source])         
        #catch werid cases in which cubefix crashes
        if not os.path.isfile(fixed):
            #try again with last step only
            print('Redo fix with last step only')
            subprocess.call(["CubeFix","-cube", cube,"-pixtable", pixtab,"-out", fixed,"-sourcemask",mask_source,"-step",errorstep]) 


        #At this step, check out cubeAdd2Mask if want to fix edges or weird ifus/slices 

        #if told to mask a particular sky region do it.. otherwise leave image empty
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
        hdulist.writeto(sharpmsk,overwrite=True)

    #now run cube skysub
    if ((os.path.isfile(skysub)) & (noclobber)):
        print("Cube {0} already skysub".format(fixed))
    else:
        print('Sky sub ', fixed)
        if(highsn):
            #now few more options to control sky sub 
            subprocess.call(["CubeSharp","-cube",fixed,"-out",skysub,"-sourcemask",sharpmsk,"-hsncube",highsn,"-lsig","5"])
        else:
            subprocess.call(["CubeSharp","-cube",fixed,"-out",skysub,"-sourcemask",sharpmsk])
                               
    #create a white image
    if ((os.path.isfile(white)) & (noclobber)):
        print("White image for cube {0} already exists".format(skysub))
    else:
        print('Create white image for ', skysub)
        subprocess.call(["Cube2Im","-cube",skysub,"-out",white])
                
def cubex_driver(listob,last=False,highsn=None,skymask=None,exthsnmask=None,version='1.6'):
    
    """
    Procedures that drives the loops of cubex within each OB folder

    listob -> the list of OBs to process
    last  -> set to True for final pass with high-sn cube 
    highsn -> name of the highsn cube used for masking 
    skymask -> mask this region in source mask before running cubesharp
    exthsnmask -> name of external hsn mask file that overrides the internally computed one
    version -> the version of cubex 
    

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
        os.chdir(ob+'/Proc/Cubex')
        print('Processing {} with cubex '.format(ob))

        #Search how many exposures are there
        scils=glob.glob("PIXTABLE_REDUCED_RESAMPLED_EXP*.fits")
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
                pixtab="PIXTABLE_REDUCED_RESAMPLED_EXP{0:d}.fits".format(dd+1)
                cube="DATACUBE_FINAL_RESAMPLED_EXP{0:d}.fits".format(dd+1)
                #now launch the task
                p = multiprocessing.Process(target=fixandsky_secondpass,args=(cube,pixtab,True,highsn,skymask,exthsnmask,version))
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
                pixtab="PIXTABLE_REDUCED_RESAMPLED_EXP{0:d}.fits".format(dd+1)
                cube="DATACUBE_FINAL_RESAMPLED_EXP{0:d}.fits".format(dd+1)
                #now launch the task
                p = multiprocessing.Process(target=fixandsky_firstpass,args=(cube,pixtab,True,skymask,version))
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
                pixtab="PIXTABLE_REDUCED_RESAMPLED_EXP{0:d}.fits".format(dd+1)
                cube="DATACUBE_FINAL_RESAMPLED_EXP{0:d}.fits".format(dd+1)
                #now launch the task
                p = multiprocessing.Process(target=fixandsky_secondpass,args=(cube,pixtab,True,None,skymask,None,version))
                workers.append(p)
                p.start()
   
            #wait for completion of all of them 
            for p in workers:
                if(p.is_alive()):
                    p.join()  
                
        #back to top
        os.chdir(topdir)

def tweak_edges(rootname):

    """

    Procedure that identifies the edges of the IFU mask and applies them to the SliceEdgeMask

    rootname: the base name to process

    """

    import numpy as np
    from astropy.io import fits
    from scipy import ndimage
    import matplotlib.pyplot as plt
    

    #open IFU image
    ifu=fits.open(rootname+"_IFUSliceMap.fits")

    #make all the same
    ifuimg=ifu[0].data
    ifuimg[np.where(ifuimg > 0)]=1e3
    ifuimg=ndimage.gaussian_filter(ifu[0].data,3)
    
    #find egdes
    sx=ndimage.sobel(ifuimg,axis=0,mode='constant')
    sy=ndimage.sobel(ifuimg,axis=1,mode='constant')
    edges=np.hypot(sx,sy)
    
    #update slice mask
    newmask=np.where(edges > 350)
    slicemask=fits.open(rootname+"_SliceEdgeMask.fits",mode='update')
    slicemask[0].data[newmask]=0
    slicemask.flush()

    #plt.imshow(ifuimg,origin='low')
    #plt.show()
    #plt.imshow(slicemask[0].data,origin='low')
    #plt.imshow(edges,origin='low',alpha=0.5)
    #plt.show()

    #close fits
    ifu.close()
    slicemask.close()
    

def drive_combine(option,listob):

    """

    Utility function that drives the procedure to combine cubes 

    
    option: -> the options to combine
               INTERMEDIATE: intermadiate products following a first pass of cubefix/cubesharp
               HIGHSN: last coadd, following the new high SN iteration 
               INDEPENDENT: split the data in two batches

    listob: -> list of OBs

    """

    import os 
    import glob
           

    #get number of OBs
    nobs=len(listob)

    if('INTERMEDIATE' in option):
        
        print('Coadd intermediate products')
        #dump to disk file lists
        topdir=os.getcwd()
        os.chdir('cubexcombine')

        fl1=open('cubes.lst','w')
        fl2=open('masks.lst','w')

        #loop over OBs
        for oob in range(nobs):
            #count how many science exposures
            nsci=len(glob.glob("../{}/Proc/Cubex/PIXTABLE_REDUCED_RESAMPLED_EXP*.fits".format(listob[oob])))
            #reconstruct names 
            for ll in range(nsci):
                fl1.write('../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_EXP{}_skysub2.fits\n'.format(listob[oob],ll+1))
                fl2.write('../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_EXP{}_fix2_SliceEdgeMask.fits\n'.format(listob[oob],ll+1))
                #tweak IFU edges - this can manually be adjusted using GUI 
                tweak_edges('../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_EXP{}_fix2'.format(listob[oob],ll+1))
                
        fl1.close()
        fl2.close()
        
        #make the temp combine
        combine_cubes("cubes.lst","masks.lst")

        #back to top 
        os.chdir(topdir)


    elif('HIGHSN' in option):
        
        print('Coadd final products')
            
        #dump to disk file lists
        topdir=os.getcwd()
        os.chdir('cubexcombine')

        fl1=open('cubes_final.lst','w')
        fl2=open('masks_final.lst','w')

        #loop over OBs
        for oob in range(nobs):
            #count how many science exposures
            nsci=len(glob.glob("../{}/Proc/Cubex/PIXTABLE_REDUCED_RESAMPLED_EXP*.fits".format(listob[oob])))
            #reconstruct names 
            for ll in range(nsci):
                fl1.write('../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_EXP{}_skysubhsn.fits\n'.format(listob[oob],ll+1))
                fl2.write('../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_EXP{}_fixhsn_SliceEdgeMask.fits\n'.format(listob[oob],ll+1))
                #tweak IFU edges - this can manually be adjusted using GUI 
                tweak_edges('../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_EXP{}_fixhsn'.format(listob[oob],ll+1))
              
        fl1.close()
        fl2.close()
        
        #make the temp combine
        combine_cubes("cubes_final.lst","masks_final.lst",final=True)
        os.chdir(topdir)


    elif('INDEPENDENT' in option):

        print('Perform half coadds')
      
        #dump to disk file lists
        topdir=os.getcwd()
        os.chdir('cubexcombine')

        #now make two independent halves 
        fl1cube=open('cubes_half1.lst','w')
        fl1mask=open('masks_half1.lst','w')
        fl2cube=open('cubes_half2.lst','w')
        fl2mask=open('masks_half2.lst','w')

        #loop over OBs
        counter=0
        for oob in range(nobs):
            #count how many science exposures
            nsci=len(glob.glob("../{}/Proc/Cubex/PIXTABLE_REDUCED_RESAMPLED_EXP*.fits".format(listob[oob])))
            #reconstruct names 
            for ll in range(nsci):
                counter=counter+1
                if(counter % 2 == 0):
                    fl1cube.write('../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_EXP{}_skysubhsn.fits\n'.format(listob[oob],ll+1))
                    fl1mask.write('../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_EXP{}_fixhsn_SliceEdgeMask.fits\n'.format(listob[oob],ll+1))
                    #tweak IFU edges - this can manually be adjusted using GUI 
                    tweak_edges('../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_EXP{}_fixhsn'.format(listob[oob],ll+1))
                else:
                    fl2cube.write('../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_EXP{}_skysubhsn.fits\n'.format(listob[oob],ll+1))
                    fl2mask.write('../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_EXP{}_fixhsn_SliceEdgeMask.fits\n'.format(listob[oob],ll+1))
                    #tweak IFU edges - this can manually be adjusted using GUI 
                    tweak_edges('../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_EXP{}_fixhsn'.format(listob[oob],ll+1))
          
    
        #close files
        fl1cube.close()
        fl1mask.close()
        fl2cube.close()      
        fl2mask.close()

        #now combine
        combine_cubes("cubes_half1.lst","masks_half1.lst",halfsetfinal='half1')
        combine_cubes("cubes_half2.lst","masks_half2.lst",halfsetfinal='half2')
        
 
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
        expmap="COMBINED_EXPMAP_FINAL.fits"
        scriptname='runcombine_final.sh'
    elif(halfset):
        cname="COMBINED_CUBE_{}.fits".format(halfset)
        iname="COMBINED_IMAGE_{}.fits".format(halfset)
        cmed="COMBINED_CUBE_MED_{}.fits".format(halfset)
        imed="COMBINED_IMAGE_MED_{}.fits".format(halfset)
        expmap="COMBINED_EXPMAP_{}.fits".format(halfset)
        scriptname='runcombine_{}.sh'.format(halfset)
    elif(halfsetfinal):
        cname="COMBINED_CUBE_FINAL_{}.fits".format(halfsetfinal)
        iname="COMBINED_IMAGE_FINAL_{}.fits".format(halfsetfinal)
        cmed="COMBINED_CUBE_MED_FINAL_{}.fits".format(halfsetfinal)
        imed="COMBINED_IMAGE_MED_FINAL_{}.fits".format(halfsetfinal)
        expmap="COMBINED_EXPMAP_FINAL_{}.fits".format(halfsetfinal)
        scriptname='runcombine_final_{}.sh'.format(halfsetfinal)
    else:
        cname="COMBINED_CUBE.fits"
        iname="COMBINED_IMAGE.fits"
        cmed="COMBINED_CUBE_MED.fits"
        imed="COMBINED_IMAGE_MED.fits"
        expmap="COMBINED_EXPMAP.fits"
        scriptname='runcombine.sh'


    if(os.path.isfile(cname)):
        print ('Cube {} already exists... skip!'.format(cname))
    else:
        print ('Creating combined cube {}'.format(cname))

        if(regions):
            print("Updating the masks")
        
            #loads list
            listmask=np.loadtxt(masks,dtype=np.dtype('a'))

            #redefine new mask 
            mask_new="new_"+masks
            llms=open(mask_new,"w")
            
            #if scalar, make it 1 element list
            if(listmask.shape == ()):
                listmask=[listmask]
                
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
                    cfits.writeto(savename,overwrite=True)
                    llms.write(savename+'\n')
                
                
                else:
                    #keep current mask 
                    llms.write(cmask+'\n')

            #done with new masks
            llms.close()

        else:
            print('Using original masks')
            mask_new=masks

        #now run combine
        print('Combine the cube...')
      
        #make mean cube - write this as script that can be ran indepedently 
        scr=open(scriptname,'w')
        scr.write("export OMP_NUM_THREADS=1\n")
        scr.write("CubeCombine -list "+cubes+" -out "+cname+" -masklist "+mask_new+" -outexp "+expmap+"\n")
        scr.write("Cube2Im -cube "+cname+" -out "+iname+"\n")
        scr.write("CubeCombine -list "+cubes+" -out "+cmed+" -masklist "+mask_new+" -comb median \n")
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
    from . import muse_utils as mutil
    from . import muse_source as msrc
    import matplotlib
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from astropy.stats import sigma_clipped_stats
    matplotlib.use('Agg')

    try:
        from photutils import CircularAperture, aperture_photometry,\
            data_properties, properties_table, centroids
    except:
        print("To run checks need photutils package")
        return

    print("Perform QA checks...")

    #move in cubex dir 
    topdir=os.getcwd()
    os.chdir('cubexcombine')

    #make QA folder
    if not os.path.exists('QA'):
        os.makedirs('QA')

    #cube names
    cname="COMBINED_CUBE_FINAL.fits"
    iname="COMBINED_IMAGE_FINAL.fits"
   
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
        subdata=rsdssall[1].data[int(roundsrc['y'][ii])-10:int(roundsrc['y'][ii])+10,\
                                     int(roundsrc['x'][ii])-10:int(roundsrc['x'][ii])+10]
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
            wname='../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_{}_whitehsn.fits'.format(thisob,thisexp)
            wfits=fits.open(wname)
            
            #now loop on sources
            delta_x=np.zeros(len(rmag))
            delta_y=np.zeros(len(rmag))
            for ii in range(len(rmag)):
                subdata=wfits[0].data[int(roundsrc['y'][ii])-10:int(roundsrc['y'][ii])+10,\
                                          int(roundsrc['x'][ii])-10:int(roundsrc['x'][ii])+10]
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
            wname='../{}/Proc/Cubex/DATACUBE_FINAL_RESAMPLED_{}_whitehsn.fits'.format(thisob,thisexp)
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
        
    #back to top dir
    os.chdir(topdir)
