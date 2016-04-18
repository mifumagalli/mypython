
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
        #prepare final names
        fixed=cube.split('.fits')[0]+"_fixhsn.fits"
        skysub=cube.split('.fits')[0]+"_skysubhsn.fits"
        white=cube.split('.fits')[0]+"_whitehsn.fits"
    else:
        #prepare intermediate names
        fixed=cube.split('.fits')[0]+"_fix2.fits"
        skysub=cube.split('.fits')[0]+"_skysub2.fits"
        white=cube.split('.fits')[0]+"_white2.fits"
        
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
            subprocess.call(["CubEx-1.5",white_source,'-MultiExt','.false.','-SN_Threshold','3','-RescaleVar','.true.'])
        else:
            print 'Using white image from previous loop'
            #create source mask 
            subprocess.call(["CubEx-1.5",white_source,'-MultiExt','.false.','-SN_Threshold','5','-RescaleVar','.true.'])
            
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
                
def cubex_driver(listob,last=False,highsn=None):
    
    """
    Procedures that drives the loops of cubex within each OB folder

    listob -> the list of OBs to process
    last  -> set to True for final pass with high-sn cube 
    highsn -> name of the highsn cube used for masking 

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
                p = multiprocessing.Process(target=fixandsky_secondpass,args=(cube,pixtab,True,highsn))
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

def combine_cubes(cubes,masks,regions=True,final=False):

    """
    Combine a bunch of cubes using masks with CubeCombine
        
    cubes    -> a list of cubes to use in the combine
    masks    -> a list of goodpix masks from the pipeline
    regions  -> if True, code searches for ds9 region files inside path with same name as 
                pipeline mask (.reg), to mask additional area that one wants to clip
    final    -> is True, append final tag to name and prepare median cubes

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
    else:
        cname="COMBINED_CUBE.fits"
        iname="COMBINED_IMAGE.fits"

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
        if(final):
            scriptname='runcombine_final.sh'
        else:
            scriptname='runcombine.sh'
        scr=open(scriptname,'w')
        scr.write("export OMP_NUM_THREADS=1\n")
        scr.write("CubeCombine -list "+cubes+" -out "+cname+" -masklist "+mask_new+"\n")
        scr.write("Cube2Im -cube "+cname+" -out "+iname+"\n")
        if(final):
            #also create median cube 
            scr.write("CubeCombine -list "+cubes+" -out "+cmed+" -masklist "+mask_new+" -comb median\n")
            scr.write("Cube2Im -cube "+cmed+" -out "+imed)
        scr.close()
        subprocess.call(["sh",scriptname])
