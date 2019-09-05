"""

Collection of procedures that drive the creation of final cubes using
eso recipies. Check out the MUSE reduction manual for details. 

"""
    

def individual_skysub(listob,autocalib=False,nproc=12):

    """
    rerun scipost performing sky subtraction with eso recipies 
    

    listobs -> a list of folders corresponding to individual OBs that needs to be 
               processed enabling skysubtraction

    nproc -> number of processors to be used in parallel mode 

    """
    
    import os
    import glob
    import subprocess
    from astropy.io import fits
    import numpy as np

    #grab top dir
    topdir=os.getcwd()

    #now loop over each folder and make the final sky-subtracted cubes
    for ob in listob:
        
        #change dir
        os.chdir(ob+'/Proc/Basic')
        print('Processing {} for sky subtraction'.format(ob))
        
        #Search how many exposures are there
        scils=glob.glob("OBJECT_RED_0*.fits*")
        nsci=len(scils)

        #Search if there are sky models from emtpy regions
        skyls=glob.glob("SKY_SPECTRUM*.fits")
        nsky=len(skyls)
        
        #Define appropriate string for autocalibration
        if autocalib:
           autocalibstr='deepfield'
        else:
           autocalibstr='none'   

        #collect sky tags and sky times 
        skytags=[]
        skytimes=[]
        for skyframes in skyls:
            skytags.append(skyframes.split('SPECTRUM_')[1].split('.fits')[0])
            thissk=fits.open(skyframes)
            skytimes.append(float(thissk[0].header['MJD-OBS']))
            
        skytimes=np.array(skytimes)

        #handle the case of sky offsets separately
        if(nsky > 0):
            #sky offsets are present
            print('Use sky offsets for sky subtraction')
            #loop on exposures and reduce frame with closest sky model 
            for exp in range(nsci):

                #check if true science and not sky 
                thisframe=scils[exp]
                thisid=thisframe.split('RED_')[1].split('.fits')[0]
                if not (thisid in skytags):
                    #you are science - proceed
                    thissc=fits.open(thisframe)
                    scitimes=float(thissc[0].header['MJD-OBS'])
                    #Find closest sky - works ok for OSOOSOOS.. 
                    #would be better to compute mid point rather than start of exposure
                    matchsky=skytags[np.argmin(abs(skytimes-scitimes))]

                    #update sof file written for basic reduction
                    sof_old=open("../../Script/scipost_{0:d}.sof".format(exp+1),"r")
                    sof_name="../../Script/scipost_{0:d}_sky.sof".format(exp+1)
                    sof=open(sof_name,"w")
                    for ll in sof_old:
                        if('SKY_LINES' in ll):
                            sof.write("SKY_LINES_{}.fits SKY_LINES\n".format(matchsky))
                        else:
                            #copy
                            sof.write(ll)
                    #append other stuff
                    sof.write("SKY_CONTINUUM_{}.fits SKY_CONTINUUM\n".format(matchsky))
                    sof.write("SKY_MASK_{}.fits SKY_MASK\n".format(matchsky))
                    sof.close()

                    #define some output names
                    cname="DATACUBE_FINAL_ESOSKY_EXP{0:d}.fits".format(exp+1)
                    pname="PIXTABLE_REDUCED_ESOSKY_EXP{0:d}.fits".format(exp+1)
                    iname="IMAGE_FOV_ESOSKY_EXP{0:d}.fits".format(exp+1)
                    
                    if not os.path.isfile(cname):
                        print("Processing exposure {0:d}".format(exp+1))

                        #Write the command file 
                        scr=open("../../Script/make_scipost_esosky_{0:d}.sh".format(exp+1),"w")
                        scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 

                        scr.write('esorex --log-file=scipost_esosky_{0:d}.log muse_scipost --filter=white --save=cube,individual --skymethod=subtract-model --autocalib={1} ../../Script/scipost_{0:d}_sky.sof'.format(exp+1,autocalibstr))
                        scr.close()

                        #Run pipeline 
                        subprocess.call(["sh", "../../Script/make_scipost_esosky_{0:d}.sh".format(exp+1)])    
                        subprocess.call(["mv","DATACUBE_FINAL.fits",cname])
                        subprocess.call(["mv","IMAGE_FOV_0001.fits",iname])
                        subprocess.call(["mv","PIXTABLE_REDUCED_0001.fits",pname])

                    else:
                        print("Exposure {0:d} exists.. skip! ".format(exp+1))
                        
            
        else:
            
            #loop on exposures and reduce frame with internal sky subtraction 
            for exp in range(nsci):
                print('Perform internal sky subtraction')
                #use sof file written for basic reduction
                sof_name="../../Script/scipost_{0:d}.sof".format(exp+1)

                #define some output names
                cname="DATACUBE_FINAL_ESOSKY_EXP{0:d}.fits".format(exp+1)
                pname="PIXTABLE_REDUCED_ESOSKY_EXP{0:d}.fits".format(exp+1)
                iname="IMAGE_FOV_ESOSKY_EXP{0:d}.fits".format(exp+1)

                if not os.path.isfile(cname):
                    print("Processing exposure {0:d}".format(exp+1))

                    #Write the command file 
                    scr=open("../../Script/make_scipost_esosky_{0:d}.sh".format(exp+1),"w")
                    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 

                    scr.write('esorex --log-file=scipost_esosky_{0:d}.log muse_scipost --filter=white --save=cube,individual --autocalib={1} ../../Script/scipost_{0:d}.sof'.format(exp+1, autocalibstr))
                    scr.close()

                    #Run pipeline 
                    subprocess.call(["sh", "../../Script/make_scipost_esosky_{0:d}.sh".format(exp+1)])    
                    subprocess.call(["mv","DATACUBE_FINAL.fits",cname])
                    subprocess.call(["mv","IMAGE_FOV_0001.fits",iname])
                    subprocess.call(["mv","PIXTABLE_REDUCED_0001.fits",pname])

                else:
                    print("Exposure {0:d} exists.. skip! ".format(exp+1))
            

        #back to top
        os.chdir(topdir)

def coaddall(listob,nproc=24):

    """
    Grab all the OBs, align them with muse_exp_align, 
    and coadd them with muse_exp_combine 
    
    nproc -> number of proc to use in parallel runs

    """

    import os
    import glob
    import subprocess

    #create align and coadd scripts and sofs
    alignsof=open('align.sof','w')
    alignscr=open('align.sh','w')
    coaddsof=open('combine.sof','w')
    coaddscr=open('combine.sh','w')
 
    for ob in listob:
        
        #grab fov
        fovimg=glob.glob("../{}/Proc/Basic/IMAGE_FOV_ESOSKY_*".format(ob))
        for ff in fovimg:
            alignsof.write("{} IMAGE_FOV\n".format(ff))
    
        #grab pixel tables 
        pixtab=glob.glob("../{}/Proc/Basic/PIXTABLE_REDUCED_ESOSKY_*".format(ob))
        for ff in pixtab:
            coaddsof.write("{} PIXTABLE_REDUCED\n".format(ff))
        
    #final touches
    coaddsof.write("OFFSET_LIST.fits OFFSET_LIST\n")
    coaddsof.close()
    alignsof.close()
    
    #script files 
    alignscr.write("esorex --log-file=exp_align.log muse_exp_align --threshold=4. align.sof")
    alignscr.close()

    coaddscr.write("OMP_NUM_THREADS={}\n".format(nproc))
    coaddscr.write("esorex --log-file=exp_combine.log muse_exp_combine combine.sof")
    coaddscr.close()
    
    #Run pipeline 
    if not os.path.isfile('DATACUBE_FINAL.fits'):
        print("Aliging exposures...")
        subprocess.call(["sh", "align.sh"])    
    else:
        print("Offset already computed!")

    if not os.path.isfile('DATACUBE_FINAL.fits'):
        print("Coadding exposures...")
        subprocess.call(["sh", "combine.sh"])
    else:
        print("Cube already combined!")


