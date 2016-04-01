"""

These are sets of procedures optimised for almost empty fields 
but with extended line emission  

"""


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
    import muse_utils as mut 
    import numpy as np

    #grab top dir
    topdir=os.getcwd()

    #now loop over each folder and make the final sky-subtracted cubes
    for ob in listob:
        
        #change dir
        os.chdir(ob+'/Proc/')
        print('Processing {} for resampling on reference cube'.format(ob))
 
        #Search how many exposures are there
        scils=glob.glob("OBJECT_RED_0*.fits*")
        nsci=len(scils)

        #loop on exposures and reduce frame with sky subtraction 
        for exp in range(nsci):
            
            if not os.path.isfile('OFFSET_LIST_EXP{0:d}.fits'.format(exp+1)):
                print("Compute offsets...")
                            
                #create align file 
                alignsof=open('../Script/align_toref_{0:d}.sof'.format(exp+1),'w')
                alignsof.write("../../{}/IMAGE_FOV_0001.fits IMAGE_FOV\n".format(refpath))
                alignsof.write("IMAGE_FOV_EXP{0:d}.fits IMAGE_FOV\n".format(exp+1))
                alignsof.close()
                
                #run script align with respect to registered reference cube  
                alignscr=open('../Script/make_align_toref_{0:d}.sh'.format(exp+1),'w')
                alignscr.write("esorex --log-file=align_toref_{0:d}.log muse_exp_align --threshold=4. ../Script/align_toref_{0:d}.sof".format(exp+1))
                alignscr.close()
                subprocess.call(["sh","../Script/make_align_toref_{0:d}.sh".format(exp+1)])    
                
                #copy the offsets 
                alig=fits.open('OFFSET_LIST.fits')
                alig.writeto('OFFSET_LIST_EXP{0:d}.fits'.format(exp+1),clobber=True)

            else:
                print('Offsets exist.. skip')

            #define some output names for final cube 
            cname="DATACUBE_FINAL_LINEWCS_EXP{0:d}.fits".format(exp+1)
            pname="PIXTABLE_REDUCED_LINEWCS_EXP{0:d}.fits".format(exp+1)
            iname="IMAGE_FOV_LINEWCS_EXP{0:d}.fits".format(exp+1)
 
            if not os.path.isfile(cname):
                print("Processing exposure {0:d} to align to reference".format(exp+1))
                
                #copy sof file written for basic reduction
                sof_old=open("../Script/scipost_{0:d}.sof".format(exp+1))
                sof_name="../Script/scipost_line_{0:d}.sof".format(exp+1)
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
                        pxt=fits.open(pixtab)
                        pxt[0].header['RA']=pxt[0].header['RA']-offsets[2]
                        pxt[0].header['DEC']=pxt[0].header['DEC']-offsets[3]
                        
                        pxt.writeto("WCS_"+pixtab,clobber=True)
                        pixtablist.append("WCS_"+pixtab)
                        sofedit.write("WCS_"+pixtab+" PIXTABLE_OBJECT\n")
                    else:
                        sofedit.write(ll)
                        
                #append reference frame to sof file 
                sofedit.write('../../{}/DATACUBE_FINAL.fits OUTPUT_WCS\n'.format(refpath))
                sofedit.close()
                sof_old.close()
                
                #Write the command file 
                scr=open("../Script/make_scipost_line_{0:d}.sh".format(exp+1),"w")
                scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
                
                scr.write('esorex --log-file=scipost_line_{0:d}.log muse_scipost --filter=white  --skymethod="none" --save=cube,individual ../Script/scipost_line_{0:d}.sof'.format(exp+1))
                scr.close()
                
                #Run pipeline 
                subprocess.call(["sh", "../Script/make_scipost_line_{0:d}.sh".format(exp+1)])    
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



def make_ifumasks(listob,refpath='./',nproc=24):

    """
    Loop over each OB and run scipost_make_cube on the reduced pixel tables
    to produce a final IFU masks for the resampled cube 

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
        print('Processing {} for IFU mask'.format(ob))
 
        #Search how many exposures are there
        scils=glob.glob("OBJECT_RED_0*.fits*")
        nsci=len(scils)

        #loop on exposures and reduce frame with sky subtraction 
        for exp in range(nsci):
            
            #define some output names for final cube 
            cname="DATACUBE_IFUMASK_LINEWCS_EXP{0:d}.fits".format(exp+1)
            iname="IMAGE_IFUMASK_LINEWCS_EXP{0:d}.fits".format(exp+1)
 
            if not os.path.isfile(cname):
                print("Processing exposure {0:d} to produce IFU mask".format(exp+1))
                
                #make a sof file 
                sof_name="../Script/scipost_ifumask_{0:d}.sof".format(exp+1)
                sofedit=open(sof_name,'w')
                                        
                #append reduced pixel table and reference frame to sof file 
                origpix='PIXTABLE_REDUCED_LINEWCS_EXP{0:d}.fits'.format(exp+1)
                newpix='IFUMASK_PIXTABLE_LINEWCS_EXP{0:d}.fits'.format(exp+1)

                sofedit.write(newpix+' PIXTABLE_OBJECT\n')
                sofedit.write('../../{}/DATACUBE_FINAL.fits OUTPUT_WCS\n'.format(refpath))
                sofedit.close()
            
                #Write the command file 
                scr=open("../Script/make_scipost_ifumask_{0:d}.sh".format(exp+1),"w")
                scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
                
                scr.write('esorex --log-file=scipost_ifumask_{0:d}.log muse_scipost_make_cube --filter=white ../Script/scipost_ifumask_{0:d}.sof'.format(exp+1))
                scr.close()

                #create the IFU mask 
                #unpack ifu origin
                pxt=fits.open(origpix)
                ifu,islice=mut.unpack_pixtab(pxt[7].data)
                
                #loop over ifu
                for iff in range(24): 
                    #group slices in 4 stacks
                    for i in range(4):
                        imin=i*12+1
                        imax=(i+1)*12 
                        
                        #find pixels and set value to flag 
                        pixinside=np.where((islice >= imin) & (islice <= imax) & (ifu == (iff+1)))
                        pxt[4].data[pixinside] = (iff+1)*100.+i+1
                        
                    print 'Done with IFU:', iff+1

                #save updated
                pxt.writeto(newpix,clobber=True)
                pxt.close()

                #Run pipeline 
                subprocess.call(["sh", "../Script/make_scipost_ifumask_{0:d}.sh".format(exp+1)])    
                subprocess.call(["mv","DATACUBE_FINAL.fits",cname])
                subprocess.call(["mv","IMAGE_FOV_0001.fits",iname])
            else:
                print("IFU mask {0:d} exists.. skip! ".format(exp+1))
     
        #back to top
        os.chdir(topdir)




def make_illcorr(listob):

    """
    Loop over each OB and perform illumination correction 

    listob -> OBs to process

    """
      
    import os
    import glob
    import subprocess
    import shutil
    from astropy.io import fits
    import muse_utils as mut 
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt


    #grab top dir
    topdir=os.getcwd()

    #now loop over each folder and make the final sky-subtracted cubes
    for ob in listob:
        
        #change dir
        os.chdir(ob+'/Proc/')
        print('Processing {} for illumination correction'.format(ob))
 
        #Search how many exposures are there
        scils=glob.glob("OBJECT_RED_0*.fits*")
        nsci=len(scils)

        #loop on exposures and reduce frame with sky subtraction 
        for exp in range(nsci):
            
            #define some names for final cube 
            ifumask_cname="DATACUBE_IFUMASK_LINEWCS_EXP{0:d}.fits".format(exp+1)
            ifumask_iname="IMAGE_IFUMASK_LINEWCS_EXP{0:d}.fits".format(exp+1)
            data_cname="DATACUBE_FINAL_LINEWCS_EXP{0:d}.fits".format(exp+1)
            data_iname="IMAGE_FOV_LINEWCS_EXP{0:d}.fits".format(exp+1)
 


            #open the fov to create a good mask 
            data=fits.open(data_iname)
            ifimask=fits.open(ifumask_iname)
            goodmask=np.nan_to_num(data[1].data)*0.

            
            #grab muse rotator
            rotation=data[0].header["HIERARCH ESO INS DROT POSANG"] 


            #start with central ifu pixels 
            for iff in range(24): 
                #group slices in 4 stacks
                for i in range(4):
                    flagvalue = (iff+1)*100.+i+1
              
                    goodpx=np.where(ifimask[1].data == flagvalue)
                    goodmask[goodpx]=1.0
            
                    

            #plt.imshow(goodmask,origin='lower')
            #plt.show()
                      

        #back to top for next OB 
        os.chdir(topdir)
