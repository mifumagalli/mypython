"""

These are sets of procedures optimised for almost empty fields using mpdaf procedures


"""
from __future__ import print_function

def selfcalibrate(listob,deepwhite):
    
    """

    Loop over each OB and performs self-calibration after masking sources as 

    listob -> OBs to process
    deepwhite -> the best white image available to mask sources


    """

    



def individual_resample(listob,refpath='./',nproc=24):

    """
    Loop over each OB and re-run scipost using a final coadded cube as
    a reference for WCS. This produces calibrated pixeltable that are all comaptible with the
    astrometric solution of the reference cube.

    If desired, one can then reconstruct a cube by running 
    > esorex  muse_scipost_make_cube makecube.sof
    with "PIXTABLE_REDUCED_RESAMPLED_EXP1.fits PIXTABLE_OBJECT" in sof file 


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
