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
    import scipy.signal as sgn
    from scipy.stats import sigmaclip

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
            
            #these are the ifu masks
            ifumask_cname="DATACUBE_IFUMASK_LINEWCS_EXP{0:d}.fits".format(exp+1)
            ifumask_iname="IMAGE_IFUMASK_LINEWCS_EXP{0:d}.fits".format(exp+1)
 
            #these are the data cubes
            data_cname="DATACUBE_FINAL_LINEWCS_EXP{0:d}.fits".format(exp+1)
            data_iname="IMAGE_FOV_LINEWCS_EXP{0:d}.fits".format(exp+1)
 
            #open the ifu mask to create a good mask 
            data=fits.open(data_cname)
            ifimask=fits.open(ifumask_iname)
            nwave=data[1].header["NAXIS3"]
            nx=data[1].header["NAXIS1"]
            ny=data[1].header["NAXIS2"]
            goodmask=np.zeros((ny,nx))
            edgemask=np.zeros((ny,nx))+1

            #grab muse rotator for this exposure
            rotation=data[0].header["HIERARCH ESO INS DROT POSANG"] 

            #make space for illumination corrections 
            illcorr=np.zeros((nwave,24,4))
            illcorrfilt=np.zeros((nwave,24,4))
            illcorrnorm=np.zeros((nwave,24,4))            
            pixinstack=np.zeros(nwave)
            
            #some save names
            outcorr="ILLCORR_EXP{0:d}.fits".format(exp+1)
            outcorrnorm="ILLCORRNORM_EXP{0:d}.fits".format(exp+1)

            if not os.path.isfile(outcorr):
                #start with central ifu pixels 
                #this loops over ifus
                for iff in range(24): 
                    print ('Computing correction for IFU {}'.format(iff+1))
                    #next loop over each stack of slices inside ifu
                    for i in range(4):
                        #reconstruct id of pixels in this stack 
                        flagvalue = (iff+1)*100.+i+1
                        #grab the pixels in that family 
                        goodpx=np.where(ifimask[1].data == flagvalue)
                        goodmask[goodpx]=1.0
                    
                        #now reconstruct spectrum
                        pixinstack=pixinstack*0.
                        for ww in range(nwave):
                            #now generate a median spectrum in this stack
                            img=np.nan_to_num(data[1].data[ww,:,:])
                            sclp,lw,up=sigmaclip(img[goodpx],3.,3.)
                            pixinstack[ww]=sclp.mean()
                
                        #stash solution
                        illcorr[:,iff,i]=pixinstack
                #save 
                hdu = fits.PrimaryHDU(illcorr)
                hdulist = fits.HDUList([hdu])
                hdulist.writeto(outcorr,clobber=True)
            else:
                print('Loading pre-computed corrections')
                illcorr=(fits.open(outcorr))[0].data
                
            #compute the actual correction by normalising to the superstack
            superstack=np.median(illcorr,axis=(1,2))
            illcorrnorm=illcorrnorm*0.
            illcorrfilt=illcorrfilt*0.
            for iff in range(24): 
                print ('Compute filtered correction for IFU {}'.format(iff+1))
                for i in range(4):
                    #normalise to superstack 
                    illcorrnorm[:,iff,i]=illcorr[:,iff,i]/superstack
                    #filter the sky lines
                    sky=sgn.medfilt(illcorrnorm[:,iff,i],99)
                    illcorrfilt[:,iff,i]=sgn.savgol_filter(sky,99,1)
                    #plt.plot(illcorrnorm[:,iff,i])
                    #plt.plot(illcorrfilt[:,iff,i])
                    #plt.plot(bb)
                    #plt.show()
                                        
            #save corrections
            hdu1 = fits.PrimaryHDU(illcorrnorm)
            hdu2 = fits.ImageHDU(illcorrfilt)
            hdulist = fits.HDUList([hdu1,hdu2])
            hdulist.writeto(outcorrnorm,clobber=True)

            #now apply corrections accounting for interpolation
            for iff in range(24): 
                print ('Apply correction for IFU {}'.format(iff+1))
                for i in range(3):
                    if(iff < 23):
                        #first, interpolation along ifus
                        x_current_ifu=(iff+1)*100.+i+1
                        x_next_ifu=(iff+2)*100.+i+1
                        #grab relevant pixels 
                        goodpx=np.where((ifimask[1].data >=x_current_ifu) & (ifimask[1].data < x_next_ifu))
                        
                        for ww in range(nwave):
                            y_current=illcorrfilt[ww,iff,i]
                            y_next=illcorrfilt[ww,iff,i]
                            slope=((y_next-y_current)/(x_next_ifu-x_current_ifu))
                            correction=y_current+slope*(ifimask[1].data[goodpx]-x_current_ifu)
                            img=np.nan_to_num(data[1].data[ww,:,:])
                            img[goodpx]=img[goodpx]/correction
                            data[1].data[ww,:,:]=img
                            
                    else:
                        #deal with last 
                        x_current_ifu=(iff+1)*100.+i+1
                        goodpx=np.where((ifimask[1].data ==x_current_ifu))
                       
                        for ww in range(nwave):
                            img=np.nan_to_num(data[1].data[ww,:,:])
                            img[goodpx]=img[goodpx]/illcorrfilt[ww,iff,i]
                            data[1].data[ww,:,:]=img

                        
                
            #edges=np.where((ifimask[1].data > current) & (ifimask[1].data < next))
            #edgemask[edges]=0.0        
            #stick in nans
            #nanpix=np.where(np.isnan(ifimask[1].data))
            #edgemask[nanpix]=0.0     
            #plt.imshow(edgemask,origin='lower')
            #plt.show()

            #save new cubes 
            data.writeto("test.fits",clobber=True)

            exit()

        #back to top for next OB 
        os.chdir(topdir)
