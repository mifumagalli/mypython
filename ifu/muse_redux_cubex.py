
def fixandsky_firstpass(cube,pixtab,noclobber):
    
    
    """ 

    Take a cube and pixel table and fix the cube, then skysub and produce white image, 
    using CubEx utils

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
        subprocess.call(["CubeFix","-cube", cube,"-pixtable", pixtab,"-out", fixed, "-edgemaskpix", '1'])

    #now run cube skysub
    if ((os.path.isfile(skysub)) & (noclobber)):
        print "Cube {0} already skysub".format(fixed)
    else:
        print 'Sky sub ', fixed
        subprocess.call(["CubeSharp","-cube",fixed,"-out",skysub,"-lcheck","False"])
                               
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

    #make some intermediate names
    fixed=cube.split('.fits')[0]+"_fix2.fits"
    skysub=cube.split('.fits')[0]+"_skysub2.fits"
    white=cube.split('.fits')[0]+"_white2.fits"


    #create a source mask 
    white_source=cube.split('.fits')[0]+"_white.fits"
    mask_source=cube.split('.fits')[0]+"_white.Objects_Id.fits"
    print 'Create source mask ', white_source
        
    if(highsn):
        print 'Using high SN cube...'
        subprocess.call(["Cube2Im","-cube",highsn,"-out",white_source])
    else:
        print 'Using white image from previous loop'
    subprocess.call(["CubEx-1.4",white_source,'-MultiExt','.false.'])
    
    #now fix the cube using masks
    if ((os.path.isfile(fixed)) & (noclobber)):
        print "Cube {0} already fixed".format(cube)
    else:
        print 'Cubefix ', cube
        subprocess.call(["CubeFix","-cube", cube,"-pixtable", pixtab,"-out", fixed,"-sourcemask",mask_source,"-edgemaskpix", '1'])

    #now run cube skysub
    if ((os.path.isfile(skysub)) & (noclobber)):
        print "Cube {0} already skysub".format(fixed)
    else:
        print 'Sky sub ', fixed
        if(highsn):
            subprocess.call(["CubeSharp","-cube",fixed,"-out",skysub,"-sourcemask",mask_source,"-hsncube",highsn])
        else:
            subprocess.call(["CubeSharp","-cube",fixed,"-out",skysub,"-sourcemask",mask_source])
                               
    #create a white image
    if ((os.path.isfile(white)) & (noclobber)):
        print "White image for cube {0} already exists".format(skysub)
    else:
        print 'Create white image for ', skysub
        subprocess.call(["Cube2Im","-cube",skysub,"-out",white])
                
    
def combine_cubes(cubes,masks,regions=None):

    """
    Combine a bunch of cubes using masks with CubeCombine
        
    cubes    -> a list of cubes to use in the combine
    masks    -> a list of goodpix masks from the pipeline
    regions  -> if any, additional ds9 regions (image units) used to 
                mask bad regions of the cube

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
    subprocess.call(["CubeCombine","-list",cubes,"-out","COMBINED_CUBE.fits","-masklist",mask_new])
    subprocess.call(["Cube2Im","-cube","COMBINED_CUBE.fits","-out","COMBINED_IMAGE.fits"])

        

    

