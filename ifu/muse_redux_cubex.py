
def fixandsky_firstpass(cube,pixtab,noclobber):
    
    
    """ 

    Take a cube and pixel table and fix the cube, then skysub and produce white image, using CubEx utils

    """

    import os 
    import subprocess

    #make some intermediate names
    fixed=cube.split('.fits')[0]+"_fix.fits"
    skysub=cube.split('.fits')[0]+"_skysub.fits"
    white=cube.split('.fits')[0]+"_white.fits"

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
        subprocess.call(["CubeSharp","-cube",fixed,"-out",skysub,"-lcheck",".false."])
                               
    #create a white image
    if ((os.path.isfile(white)) & (noclobber)):
        print "White image for cube {0} already exists".format(skysub)
    else:
        print 'Create white image for ', skysub
        subprocess.call(["Cube2Im","-cube",skysub,"-out",white])
                

def fixandsky_secondpass(cube,pixtab,noclobber,highsn=None):
        
    """ 
 
    Similar to first pass, but operates on cubes that have been realigned and uses masks 
    as appropriate

    If highsn cube is provided, use it to mask the sources and for skysubtraction
 
    """
      
    import os 
    import subprocess

    if(highsn):
        #make some intermediate names
        fixed=cube.split('.fits')[0]+"_fixhsn.fits"
        skysub=cube.split('.fits')[0]+"_skysubhsn.fits"
        white=cube.split('.fits')[0]+"_whitehsn.fits"
    else:
        #make some intermediate names
        fixed=cube.split('.fits')[0]+"_fix2.fits"
        skysub=cube.split('.fits')[0]+"_skysub2.fits"
        white=cube.split('.fits')[0]+"_white2.fits"
        
    #create a source mask 
    mask_source=cube.split('.fits')[0]+"_white.Objects_Id.fits"
    white_source=cube.split('.fits')[0]+"_white.fits"
    print 'Create source mask ', white_source
    
    #if high cube provide, overwrite white image 
    if(highsn):
        print 'Using high SN cube...'
        subprocess.call(["Cube2Im","-cube",highsn,"-out",white_source])
        subprocess.call(["CubEx-1.5",white_source,'-MultiExt','.false.'])
    else:
        print 'Using white image from previous loop'
        #create source mask 
        subprocess.call(["CubEx-1.5",white_source,'-MultiExt','.false.','-SN_Threshold','5','-RescaleVar','.true.'])
    
    #now fix the cube using masks
    if ((os.path.isfile(fixed)) & (noclobber)):
        print "Cube {0} already fixed".format(cube)
    else:
        print 'Cubefix ', cube
        subprocess.call(["CubeFix","-cube", cube,"-pixtable", pixtab,"-out", fixed,"-sourcemask",mask_source])

    #At this step, check out cubeAdd2Mask if want to fix edges or weird ifus/slices 

    #now run cube skysub
    if ((os.path.isfile(skysub)) & (noclobber)):
        print "Cube {0} already skysub".format(fixed)
    else:
        print 'Sky sub ', fixed
        if(highsn):
            #now few more options to control sky sub 
            subprocess.call(["CubeSharp","-cube",fixed,"-out",skysub,"-sourcemask",mask_source,"-hsncube",highsn])
        else:
            subprocess.call(["CubeSharp","-cube",fixed,"-out",skysub,"-sourcemask",mask_source])
                               
    #create a white image
    if ((os.path.isfile(white)) & (noclobber)):
        print "White image for cube {0} already exists".format(skysub)
    else:
        print 'Create white image for ', skysub
        subprocess.call(["Cube2Im","-cube",skysub,"-out",white])
                
def cubex_driver(listob):
    
    """
    Procedures that drives the loops of cubex within each OB folder
    listob -> the list of OBs to process
    final  -> set to True for final pass with high-sn cube 

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
            p = multiprocessing.Process(target=fixandsky_firstpass,args=(cube,pixtab,True,))
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
            p = multiprocessing.Process(target=fixandsky_secondpass,args=(cube,pixtab,True,))
            workers.append(p)
            p.start()
   
        #wait for completion of all of them 
        for p in workers:
            if(p.is_alive()):
                p.join()  
                
        #back to top
        os.chdir(topdir)

def combine_cubes(cubes,masks,regions=None,final=False):

    """
    Combine a bunch of cubes using masks with CubeCombine
        
    cubes    -> a list of cubes to use in the combine
    masks    -> a list of goodpix masks from the pipeline
    regions  -> if any, additional ds9 regions (image units) used to 
                mask bad regions of the cube

    final    -> is True, append final tag to name 

    """
    import subprocess
    import os
    import numpy as np
    from astropy.io import fits
    from mypython.fits import pyregmask as msk

    if(regions):
        print "Updating the masks"
        
        listmask=np.loadtxt(masks,dtype=np.dtype('a'))
        listreg=np.loadtxt(regions,dtype=np.dtype('a'))
        
        mask_new="new_"+masks
        llms=open(mask_new,"w")


        #loop over and update with regions
        for i,cmask in enumerate(listmask):
            
            #open fits
            cfits=fits.open(cmask)
            
            #init reg mask
            Mask = msk.PyMask(cfits[0].header["NAXIS1"],cfits[0].header["NAXIS2"],listreg[i])
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
        
        #done with new masks
        llms.close()

    else:
        print 'Using original masks'
        mask_new=masks

    #now run combine
    print 'Combine the cube...'    
    my_env = os.environ
    my_env["OMP_NUM_THREADS"] = "1"

    if(final):
        subprocess.call(["CubeCombine","-list",cubes,"-out","COMBINED_CUBE_FINAL.fits","-masklist",mask_new])
        subprocess.call(["Cube2Im","-cube","COMBINED_CUBE_FINAL.fits","-out","COMBINED_IMAGE_FINAL.fits"])
    else:
        subprocess.call(["CubeCombine","-list",cubes,"-out","COMBINED_CUBE.fits","-masklist",mask_new])
        subprocess.call(["Cube2Im","-cube","COMBINED_CUBE.fits","-out","COMBINED_IMAGE.fits"])
