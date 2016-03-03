
class Muse(object):

    """

    This is a class that bundles some of the procedures to handle muse data


    """


    def __init__(self):

        """ Stuff at init .. not much for now """
        
        print "Initialisation of Muse object"
        self.musepip=1.0
        

    def redux_basic(self,path='./',nproc=12):
        
        """ 

        This is a basic engine that performs redux of MUSE data using the eso pipeline in a basic 
        form, that is apply all the basic calibs but stop before sky subtraction and coaddition.
        This can be done in a later step, after some post-processing of the cube for enhanced 
        data quality

        path - the top level folder where data reduction has to be performed.
               This includes a folder Raw where the data have been downloaded using the 
               eso shell script. It assumes that basic calibrations are also included, as 
               provided by eso archive. 
        
        nproc - the number of processors to use during the reduction 

        This code is designed to handle a single OB or groups of OBs that share the same sets of calibrations
 
        """
        
        import muse_redux_basic as rdx
        import os

        print 'Starting reduction...'
        
        #First, make sure the various folders exist as needed 
        if not os.path.exists(path+"Raw"):
            print "Cannot find Raw data..."
            exit()
        if not os.path.exists(path+"Script"):
            os.makedirs(path+"Script")
        if not os.path.exists(path+"Proc"):
            os.makedirs(path+"Proc")

        #parse the xml file(s) 
        xml_info=rdx.parse_xml(path=path,nproc=nproc)
        
        #now start reduction. Enter the proc folder
        currdir=os.getcwd()
        os.chdir(path+'Proc')
        print 'Changing dir to proc...'
        
        #First handle the bias
        if not os.path.isfile("MASTER_BIAS.fits"):
            print 'Creating bias...'
            rdx.make_bias(xml_info,nproc=nproc)
            print 'All done with the bias...'
        else:
            print 'Bias already exist'
            
        #Next handle the dark
        if not os.path.isfile("MASTER_DARK.fits"):
            print 'Creating dark...'
            rdx.make_dark(xml_info,nproc=nproc)
            print 'All done with the dark...'
        else:
            print 'Dark already exist'
            
        #Next handle the flats
        if not os.path.isfile("MASTER_FLAT.fits"):
            print 'Creating flat...'
            rdx.make_flat(xml_info,nproc=nproc)
            print 'All done with flat...'
        else:
            print 'Flat already exist'
  
        #Next handle the arcs
        if not os.path.isfile("WAVECAL_RESIDUALS.fits"):
            print 'Processing the arcs...'
            rdx.make_arcs(xml_info,nproc=nproc)
            print 'All done with arcs...'
        else:
            print 'Arcs already processed'
            
        #Next handle the twilight flat
        if not os.path.isfile("DATACUBE_SKYFLAT.fits"):
            print 'Processing the twiflat...'
            rdx.make_twiflat(xml_info,nproc=nproc)
            print 'All done with twiflat...'
        else:
            print 'Twiflat already processed'

        #Next calibrate standard star
        if not os.path.isfile("STD_RED_0001.fits"):
            print 'Processing the standard star...'
            rdx.make_stdstar(xml_info,nproc=nproc)
            print 'All done with standard star...'
        else:
            print 'Standard star already processed'
                
        #Next generate flux table
        if not os.path.isfile("STD_FLUXES_0001.fits"):
            print 'Processing the flux table...'
            rdx.make_stdflux(xml_info,nproc=nproc)
            print 'All done with flux table...'
        else:
            print 'Flux table already processed'
      
        #Next calibrate objects
        if not os.path.isfile("OBJECT_RED_0001.fits"):
            print 'Processing the objects...'
            rdx.make_objects(xml_info,nproc=nproc)
            print 'All done with objects...'
        else:
            print 'Objects already processed'

        #Finally, process science
        print 'Preparing intermediate data cubes...'
        rdx.make_cubes(xml_info,nproc=nproc)
        
        #Done - back to original directory!
        print 'All done with basic redux...'
        os.chdir(currdir)
        
        return xml_info

    def cubex_process(self,catalogue,xml_info,path='./',refcube=None,highsn=None):
        
        """ 
        
        Take individual cubes and process them with CubEx to tweak the illumination and perform sky subtraction
        Needs a catalogue of sources in the field, (ra, dec) one per line, so as to realign the WCS to 
        absolute references
        
        CubEx is a privite code by S. Cantalupo and cannot be redistributed. Contact the author directly. 

        path - the top level folder where data reduction has to be performed.
               This includes a folder Proc with the processed datacube  
         
        catalogue is catalogue used for ra dec alignement
        xml_info  is xml_info used in reduction
        refcube  a reference cube used to resample WCS
        highsn  a combined, sky subtracted, high signal to noise cube used for the last pass 
                of sky subtraction 

        """

        #now start reduction. Enter the proc folder
        import os
        import glob
        import subprocess
        import muse_redux_cubex as cx 
        import muse_utils as utl
        import muse_redux_basic as rdx
        import multiprocessing
        import numpy as np

        currdir=os.getcwd()
        os.chdir(path+'Proc')
        print 'Changing dir to proc...'

        ########################################
        # do a first loop of cubex fix and sky #
        ########################################

        #Search how many exposures are there
        scils=glob.glob("OBJECT_RED_0*.fits*")
        nsci=len(scils)

        #do it in parallel on exposures
        workers=[]
        for dd in range(nsci):
            #reconstruct the name 
            pixtab="PIXTABLE_REDUCED_EXP{0:d}.fits".format(dd+1)
            cube="DATACUBE_FINAL_EXP{0:d}.fits".format(dd+1)
            #now launch the task
            p = multiprocessing.Process(target=cx.fixandsky_firstpass,args=(cube,pixtab,True,))
            workers.append(p)
            p.start()
   
        #wait for completion of all of them 
        for p in workers:
            if(p.is_alive()):
                p.join()

        
        ########################################
        # compute offsets and re-project cubes #
        ########################################
   

        #recipy provided for this in v1.2
    
        #space
        ra_off=[]
        dec_off=[]
        
        if not os.path.isfile('radecoffsets.txt'):
            print 'Compute offsets to reference...'

            for dd in range(nsci):
                #reconstruct the name 
                cube="DATACUBE_FINAL_EXP{0:d}_skysub.fits".format(dd+1)
                #Compute the offsets
                offra,offdec=utl.aligntocat(cube,currdir+'/'+catalogue)
                ra_off.append(offra)
                dec_off.append(offdec)
                
            #write to disk offsets
            np.savetxt('radecoffsets.txt',(ra_off,dec_off))

        else:
            print 'Offsets file exists.. loading'
            off=open('radecoffsets.txt','r')
            ra_off=[float(ii) for ii in off.readline().strip().split()]
            dec_off=[float(ii) for ii in off.readline().strip().split()]
            off.close()
            
        #now reproject cubes with offsets
        if(refcube):
            print 'Using reference cube for wcs...'
            rdx.make_cubes(xml_info,nproc=12,wcsoff=[ra_off,dec_off],refcube=currdir+'/'+refcube) 
        else:
            print 'Computing reference WCS on data themselves...'
            rdx.make_cubes(xml_info,nproc=12,wcsoff=[ra_off,dec_off])  
            
        ######################################################
        # run a second iter of cubex fix and sky             #
        # with realigned images  mostly to find sources for  #
        # masking                                            #
        ###################################################### 
        
        #do it in parallel on exposures
        workers=[]
        for dd in range(nsci):
            #reconstruct the name 
            pixtab="PIXTABLE_REDUCED_EXP{0:d}_off.fits".format(dd+1)
            cube="DATACUBE_FINAL_EXP{0:d}_off.fits".format(dd+1)
            #now launch the task
            p = multiprocessing.Process(target=cx.fixandsky_firstpass,args=(cube,pixtab,True,))
            workers.append(p)
            p.start()
   
        #wait for completion of all of them 
        for p in workers:
            if(p.is_alive()):
                p.join()

        ##########################################
        # run a third iter of cubex fix and sky  #
        # with better masking on sources         #
        ##########################################
        if(highsn):
            #run the loop
            print 'Do third pass of cubex..'
            #do it in parallel on exposures
            workers=[]
            for dd in range(nsci):
                #reconstruct the name 
                pixtab="PIXTABLE_REDUCED_EXP{0:d}_off.fits".format(dd+1)
                cube="DATACUBE_FINAL_EXP{0:d}_off.fits".format(dd+1)
                #now launch the task
                p = multiprocessing.Process(target=cx.fixandsky_secondpass,args=(cube,pixtab,True,currdir+'/'+highsn))
                workers.append(p)
                p.start()
   
            #wait for completion of all of them 
            for p in workers:
                if(p.is_alive()):
                    p.join()  
        else:
            print 'High SN cube not provided.. skip step 3'
      
        #Done - back to original directory!
        print 'All done with fancy redux... Ready to coadd'
        os.chdir(currdir)
    
#### End of the clean and tested code #####
#### End of the clean and tested code #####
#### End of the clean and tested code #####
#### End of the clean and tested code #####
#### End of the clean and tested code #####

def oldcode():

    """ This is old code that needs to be sorted out """


    #now find the illumination correction if not there yet: this owerwrite the object pixel tables 
    if (illcorr == True):
        for exp in range(nsci):
            if not os.path.isfile("PIXTABLE_REDUCED_EXP{0:d}_ill.txt".format(exp+1)):
                muse_illumcor("DATACUBE_FINAL_EXP{0:d}.fits".format(exp+1),"PIXTABLE_REDUCED_EXP{0:d}.fits".format(exp+1))

                #now apply t he illumination correction: this overwrites object pixels tables
                #read back the corrections
                icor = open("PIXTABLE_REDUCED_EXP{0:d}_ill.txt".format(exp+1),'r')
                corrill=[]
                corrifuid=[]
                for ll in icor:
                    if('#' not in ll):
                        corrill.append(float(ll.split(" ")[0]))
                        corrifuid.append(int(ll.split(" ")[2]))

                #now loop on pixel tables to apply the correction
                for thisifu in range(24):
                    
                    #open
                    hdulist = fits.open("PIXTABLE_OBJECT_{0:04d}-{1:02d}.fits".format(exp+1,thisifu+1),mode='update')
                    
                    #sanity check
                    if(corrifuid[thisifu] != thisifu+1):
                        print "Check order of illumination correction"
                        exit()
                
                    hdulist['data'].data=hdulist['data'].data*corrill[thisifu]
                    hdulist['stat'].data=hdulist['stat'].data*(corrill[thisifu]**2)
                        
                    #write to disk
                    hdulist.flush()
                    hdulist.close()
                    
                #once done, re-run the pipeline again to get the non-sky sub cube but with illumination corrections
                subprocess.call(["sh", "../Script/make_scipost_{0:d}.sh".format(exp+1)])    
                subprocess.call(["mv","DATACUBE_FINAL.fits","DATACUBE_FINAL_EXP{0:d}.fits".format(exp+1)])
                subprocess.call(["mv","IMAGE_FOV_0001.fits","IMAGE_FOV_EXP{0:d}.fits".format(exp+1)])
                subprocess.call(["mv","PIXTABLE_REDUCED_0001.fits","PIXTABLE_REDUCED_EXP{0:d}.fits".format(exp+1)])
                subprocess.call(["mv","PIXTABLE_POSITIONED_0001.fits","PIXTABLE_POSITIONED_EXP{0:d}.fits".format(exp+1)])
    
    #Now create the skysubtracted cube using a mask
    #Going for the pixel fract will mean that the majority of the pixels are 
    #selected from a few IFUs... No good. Could do better.
    for exp in range(nsci):
        
        if not os.path.isfile("PIXTABLE_REDUCEDSKY_EXP{0:d}.fits".format(exp+1)):

            print "Processing exposure {0:d} for sky subtraction".format(exp+1)
            
            ##open the datacube
            #project_cube("DATACUBE_FINAL_EXP{0:d}.fits".format(exp+1),[4800,8500],"fov_exp_{0:d}.fits".format(exp+1))
            ##load the fov projected 
            #hdulist=fits.open("fov_exp_{0:d}.fits".format(exp+1))
            #
            ##construct SN image
            #snimage=np.copy(hdulist[1].data)
            #good=np.where(hdulist[2].data > 0)
            #snimage[good]=(snimage[good]-np.median(snimage[good]))/np.sqrt(hdulist[2].data[good])
            #highsn=np.where(snimage > 5)
            #xpix=highsn[0]
            #ypix=highsn[1]
            #
            #mask=np.copy(hdulist[1].data)*0+1
            # 
            #hdulist.close()
            #
            ##mask +/- 3 pixels around every pixel with SN>value
            #pixbuf=3
            #for pp in range(xpix.size):
            #    mask[xpix[pp]-pixbuf:xpix[pp]+pixbuf,ypix[pp]-pixbuf:ypix[pp]+pixbuf]=0
            #
            #hdu = fits.PrimaryHDU(mask)
            #hdu.writeto("skymask_exp_{0:d}.fits".format(exp+1),clobber=True)

            #Write the sof file 
            sof=open("../Script/scipost_{0:d}_sky.sof".format(exp+1),"w")
            for i in range(len(flname)):
                if(fltype[i] == 'ASTROMETRY_WCS'):
                    sof.write("../Raw/{0}.fits {1}\n".format(flname[i],fltype[i])) 
                if(fltype[i] == 'SKY_LINES'):
                    sof.write("../Raw/{0}.fits {1}\n".format(flname[i],fltype[i])) 
                if(fltype[i] == 'EXTINCT_TABLE'):
                    sof.write("../Raw/{0}.fits {1}\n".format(flname[i],fltype[i])) 
                if(fltype[i] == 'FILTER_LIST'):
                    sof.write("../Raw/{0}.fits {1}\n".format(flname[i],fltype[i])) 
            for ifu in range(24):
                sof.write("PIXTABLE_OBJECT_{0:04d}-{1:02d}.fits PIXTABLE_OBJECT\n".format(exp+1,ifu+1)) 
            #sof.write("skymask_exp_{0:d}.fits SKY_MASK\n".format(exp+1))
            sof.write("STD_RESPONSE_0001.fits STD_RESPONSE\n")
            sof.write("STD_TELLURIC_0001.fits STD_TELLURIC\n")
            sof.write("../Calibs/lsf_profile.fits LSF_PROFILE\n") 
            sof.close()

            #Write the command file 
            scr=open("../Script/make_scipostsky_{0:d}.sh".format(exp+1),"w")
            scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
            scr.write('esorex --log-file=scipost_{0:d}.log muse_scipost --filter=white --save=cube,skymodel,individual --skymodel_fraction=0.15 ../Script/scipost_{0:d}_sky.sof'.format(exp+1))
            scr.close()
    
            #Run pipeline 
            subprocess.call(["sh", "../Script/make_scipostsky_{0:d}.sh".format(exp+1)])    
            subprocess.call(["mv","DATACUBE_FINAL.fits","DATACUBE_FINALSKY_EXP{0:d}.fits".format(exp+1)])
            subprocess.call(["mv","IMAGE_FOV_0001.fits","IMAGE_FOVSKY_EXP{0:d}.fits".format(exp+1)])
            subprocess.call(["mv","PIXTABLE_REDUCED_0001.fits","PIXTABLE_REDUCEDSKY_EXP{0:d}.fits".format(exp+1)])
            subprocess.call(["mv","SKY_CONTINUUM_0001.fits","SKY_CONTINUUM_EXP{0:d}.fits".format(exp+1)])
            subprocess.call(["mv","SKY_LINES_0001.fits","SKY_LINES_EXP{0:d}.fits".format(exp+1)])
            subprocess.call(["mv","SKY_MASK_0001.fits","SKY_MASK_EXP{0:d}.fits".format(exp+1)])
            subprocess.call(["mv","SKY_SPECTRUM_0001.fits","SKY_SPECTRUM_EXP{0:d}.fits".format(exp+1)])
            
    print "Done with single OB reduction!"

    return


def muse_illumcor(cube,pixeltab):

    """ 
    This is a function where the magic of the illumination correction happens 

    """
    import astropy
    from astropy.io import fits
    import numpy as np
    import mypython as mp
    #import matplotlib
    #import matplotlib.pyplot as plt
    import scipy 
    from scipy.optimize import curve_fit
    from scipy import asarray as ar,exp
    from scipy import interpolate
    from scipy.ndimage import filters


    #open the pixel table
    hdulist = fits.open(pixeltab)

    #unpack IFU and slice as in muse_pixtable.c (ESO pipe v1)
    MUSE_ORIGIN_SHIFT_IFU = 6
    ifu=(hdulist['origin'].data >> MUSE_ORIGIN_SHIFT_IFU)  & 0x1f
    slic=hdulist['origin'].data & 0x3f

    #define a set of linelists to be used for illumination correction 
    #linelist = [5198.5977,5577.339,6300.304,6363.776,6533.0527,6923.1528,6948.7515,7341.0840,7821.2354,\
    #                ,7993.343,8399.202,8430.215,8885.914,8919.696]
   
    #This is a set that works ok: if you add more lines, check that are good ones 
    linelist = [5577.339,6300.304,6363.776,7341.0840,7821.2354,7993.343,8399.202,8430.215]
    
    #linelist = [5577.339,6300.304,6363.776,7993.343,8399.202,8430.215,8885.914,8919.696]
    #init the illcorrection - one for each IFU  (could add slice stack)
    illcorr=np.zeros((len(linelist),24))
    illcorr_spline=np.zeros((len(linelist),24))
    masterflux=np.zeros(len(linelist))
    masterflux_spline=np.zeros(len(linelist))

    def gaus(x,a,x0,sigma):
        """ Simple utility function """
        return a*exp(-(x-x0)**2/(2*sigma**2))
    
    for li in range(len(linelist)):
    
        #find min max range for a given line 
        line = linelist[li]
        wstart=line-5
        wend=line+5

        print "Computing illumination correction over each IFU for line ", line
   
        #find relevant pixels
        index=np.nonzero((hdulist['lambda'].data >= wstart) & (hdulist['lambda'].data < wend))
        index=index[0]
    
        #sky subtract the data using edges 
        indxsky=np.nonzero(((hdulist['lambda'].data >= wstart) & (hdulist['lambda'].data < wstart+2.5)) |
                           ((hdulist['lambda'].data >= wend-2.5) & (hdulist['lambda'].data < wend)))
        skyval=np.median(hdulist['data'].data[indxsky[0]])
                
        #copy data in decent format withinthe range of interest and sky subtract
        x=(hdulist['lambda'].data[index])[:,0]
        y=(hdulist['data'].data[index])[:,0]-skyval
        currifu=(ifu[index])

        #find the spline/gaussian representation for the data
        #First create a data structure - deal with only unique data
        wunique,windex=np.unique(x,return_index=True)
        spdata_type=[('wave',float), ('flux', float)]
        spdata_val = np.zeros((wunique.size,), dtype=spdata_type)
        spdata_val['wave']=x[windex]
        spdata_val['flux']=y[windex]

        #now sort in wavelength 
        spdata_val=np.sort(spdata_val, order='wave') 

        #now do some median filtering on 5% of data
        smwin=int(spdata_val['flux'].size*0.05)
        spdata_val['flux']=filters.median_filter(spdata_val['flux'], size=smwin)

        #feed to interpolation function 
        linespl = interpolate.splrep(spdata_val['wave'],spdata_val['flux'])
        #fit the gaussian
        popt,pcov = curve_fit(gaus,spdata_val['wave'],spdata_val['flux'],p0=[1000,line,3])
  
        #compute the gaussian integral
        masterflux[li]=popt[0]*popt[2]*np.sqrt(np.pi)
        #integrate the spline at 4sigma
        nsg=3.5
        masterflux_spline[li]=interpolate.splint(line-nsg*popt[2],line+nsg*popt[2],linespl)

        #This is to check the fit
        #xnew = np.arange(wstart,wend,0.1)
        #yspl=interpolate.splev(xnew, linespl)
        #ygaus=gaus(xnew,*popt)
        #plt.plot(spdata_val['wave'],spdata_val['flux'], linewidth=5)
        #plt.plot(xnew,ygaus, color='black', linewidth=2)
        #plt.plot(xnew,yspl, color='red', linewidth=2)
        #plt.plot([line-nsg*popt[2],line-nsg*popt[2]],[0,np.max(spdata_val['flux'])], linestyle=':', linewidth=2, color='black')
        #plt.plot([line+nsg*popt[2],line+nsg*popt[2]],[0,np.max(spdata_val['flux'])], linestyle=':', linewidth=2, color='black')
        #plt.show()
            
        #loop over each IFU
        for thisifu in range(24):
            
            #find the pixels in each IFU
            ifuid=thisifu+1

            #repeat the procedure with only a subset of pixels
            indxifu=np.nonzero(currifu == ifuid)
            indxifu=indxifu[0]
            xifu=x[indxifu]
            yifu=y[indxifu]
                     

            #find the spline/gaussian representation
            #First create a data structure - deal with only unique data
            wunique,windex=np.unique(xifu,return_index=True)
            spdata_type=[('wave',float), ('flux', float)]
            spdata_val = np.zeros((wunique.size,), dtype=spdata_type)
            spdata_val['wave']=xifu[windex]
            spdata_val['flux']=yifu[windex]
            #now sort in wavelength 
            spdata_val=np.sort(spdata_val, order='wave') 
            #now do some median filtering on 5% of data
            smwin=int(spdata_val['flux'].size*0.05)
            spdata_val['flux']=filters.median_filter(spdata_val['flux'], size=smwin)
            #feed to interpolation
            linesplifu = interpolate.splrep(spdata_val['wave'],spdata_val['flux'])
            #fit the gaussian
            poptifu,pcov = curve_fit(gaus,spdata_val['wave'],spdata_val['flux'],p0=[1000,line,3])
  
            #compute the gaussian integral            
            illcorr[li,thisifu]=poptifu[0]*poptifu[2]*np.sqrt(np.pi)
            #integrate the spline at 4sigma
            illcorr_spline[li,thisifu]=interpolate.splint(line-nsg*poptifu[2],line+nsg*poptifu[2],linesplifu)

            #This is to check the fit
            #xnew = np.arange(wstart,wend,0.1)
            #yspl=interpolate.splev(xnew, linesplifu)
            #ygaus=gaus(xnew,*poptifu)
            #plt.plot(spdata_val['wave'],spdata_val['flux'], linewidth=5)
            #plt.plot(xnew,ygaus, color='black', linewidth=2)
            #plt.plot(xnew,yspl, color='red', linewidth=2)
            #plt.plot([line-nsg*poptifu[2],line-nsg*poptifu[2]],[0,np.max(spdata_val['flux'])], linestyle=':', linewidth=2, color='black')
            #plt.plot([line+nsg*poptifu[2],line+nsg*poptifu[2]],[0,np.max(spdata_val['flux'])], linestyle=':', linewidth=2, color='black')
            #plt.show()
                    
    #Now apply the illumination correction to pixels
    print "Writing illumination correction over each IFU"
    newname=pixeltab.split('.fits')[0]+'_ill.txt'
    icf=open(newname,"w")
    icf.write("# Correction Error IFU Nlines\n")
     
    #loop over each IFU
    for thisifu in range(24):
        
        #define IFU id
        ifuid=thisifu+1
        
        #Check corrections as a function of wavelengths
        #plt.scatter(linelist,masterflux/illcorr[:,thisifu])
        #plt.scatter(linelist,masterflux_spline/illcorr_spline[:,thisifu], color='red')
        #plt.show()
            
        #find the actual correction, rejecting bad fits
        thisill=masterflux/illcorr[:,thisifu]
        googlines=np.nonzero((thisill < 1.1) & (thisill > 0.9))
        correction=np.median(thisill[googlines[0]])
        error=np.std(thisill[googlines[0]])/np.sqrt(len(googlines[0]))
        numl=len(googlines[0])

        thisill_spline=masterflux_spline/illcorr_spline[:,thisifu]
        googlines_spline=np.nonzero((thisill_spline < 1.1) & (thisill_spline > 0.9))
        correction_spline=np.median(thisill_spline[googlines_spline[0]])
        error_spline=np.std(thisill_spline[googlines_spline[0]])/np.sqrt(len(googlines_spline[0]))
        numl_spline=len(googlines_spline[0])

        #go for the spline version
        icf.write("{0:f} {1:f} {2:d} {3:d}\n".format(correction_spline,error_spline,thisifu+1,numl_spline))
        
    #close io
    icf.close()
    hdulist.close()
     
    return

def project_cube(cubefits,filt,write):

    """ 

    This is a function that projects a cube in a given filter 

    cubefits -> the name of a resampled cube
    filter -> the filter used for transmission curve or a wavelength range for top-hat 
    write -> name of output (False will skip writing)

    
    Check out the more basic but vectorized version at muse_utils.cube2img


    """

    import astropy
    from astropy.io import fits 
    import numpy as np
    from scipy.interpolate import interp1d
    from astropy import wcs 

    #prepare the shape of the filter
    if(len(filt) > 1):

        #this is for top-hat
        lambdaf=np.arange(filt[0],filt[1]+1,1.)
        transf=lambdaf-lambdaf+1.
        ledge=filt[0]
        ridge=filt[1]

    else:

        #this is for filter curve
        #Not coded yet...
        
        print "Sorry need to code this up..."
        exit()
        

    #now open the cube
    hdulist = fits.open(cubefits)
    header=hdulist[1].header
    
    #reconstruct wavelength 
    delta_lambda=hdulist[1].header["CD3_3"]
    zero_lambda=hdulist[1].header["CRVAL3"]
    nwav=len(hdulist[1].data[:,0,0])
    #make wavelength array in air 
    wavecube=np.arange(0,nwav,1)*delta_lambda+zero_lambda
    #check for helio correction 
    try: 
        heliocorr=hdulist[1].header["HELIO"]
        obswave=wavecube/heliocorr #undo to observed frame
    except KeyError:
        print "Helio correction not applied"
        obswave=wavecube #this is air no helio
    
    #check if the cube is in range 
    if ((np.max(wavecube) < ridge) | (np.min(wavecube) > ledge)):
        print 'Cube does not overlap fully with filters....'
        
    #if yes, then interpolate filter over cube in range of interest 
    overlap=np.nonzero((wavecube >= ledge) & (wavecube <= ridge))
    overlap=overlap[0]
    fint = interp1d(lambdaf,transf)
    
    #now prepare for compression 
    sz=[hdulist[1].header["NAXIS1"],hdulist[1].header["NAXIS2"]]
    colorimage=np.zeros((sz[1],sz[0]))
    colorwgt=np.zeros((sz[1],sz[0]))
    colorvar=np.zeros((sz[1],sz[0]))

    #kill nans in data/variance
    nans=np.nonzero(np.isnan(hdulist["DATA"].data))
    mask=np.zeros((hdulist["DATA"].data).shape)+1
    hdulist["DATA"].data[nans]=1
    hdulist["STAT"].data[nans]=1
    mask[nans]=0
      
    #find dltal
    deltal=wavecube-np.roll(wavecube,1)
    deltal[0]=wavecube[1]-wavecube[0]
    
    #now fill images 
    for ii in range(overlap.size):
        
        #map wave index 
        windex=overlap[ii]
        colorimage=colorimage+hdulist["DATA"].data[windex,:,:]*fint(wavecube[windex])*deltal[windex]*mask[windex,:,:]
        colorwgt=colorwgt+fint(wavecube[windex])*mask[windex,:,:]*deltal[windex]
        colorvar=colorvar+hdulist["STAT"].data[windex,:,:]*mask[windex,:,:]*(fint(wavecube[windex])*deltal[windex])**2
        
    #compress  
    colorwgt[np.where(colorwgt <= 0)]=1
    outimg=colorimage/colorwgt
    outvar=colorvar/(colorwgt)**2
    
    #extract wcs for later use
    w=wcs.WCS(header)
    #close
    hdulist.close()
    
    if (write != False):

        #write 
        header = w.to_header()
        
        #now append CD to header
        header['CD1_1']=header['PC1_1']
        header['CD2_2']=header['PC2_2']
        header['CD1_2']=0.0
        header['CD2_1']=0.0
      
        hdu = fits.PrimaryHDU(header=header)
        hdu1 = fits.ImageHDU(outimg,header=header)
        hdu2 = fits.ImageHDU(outvar,header=header)
        hdulist = fits.HDUList([hdu,hdu1,hdu2])
        hdulist.writeto(write,clobber=True)

    
    return  {"img":outimg, "var":outvar, "header":header}


def muse_combine(pixtab,cubes,wcslist,lmin=4600,lmax=10000):


    """ 
    Take a list of reduced pixel table and cubes and combine them. 
   
    pixtab   = reduced pixel table for alignment  
    wcslist  = list of ra dec of bright sources in the field for centroid
               that can be used for relative source alignment. 
    cubes    = reduced cubes, which are not directly combined, but they are used to compute alignemnt 
    lmin,lmax = min max wave interval for combining 

    Make sure filter list is in working directory

    """
    
    import astropy 
    from astropy.table import Table
    from astropy.io import fits
    import PyGuide
    import numpy as np
    from astropy import wcs
    from astropy.coordinates import SkyCoord
    import subprocess

    #load the wcs
    wcstab = Table.read(wcslist,format='ascii')

    raoff=[]
    decoff=[]

    #loop over the pixel table to compute the offsets
    for cb in cubes:
        
        print "Processing cube {} for offsets".format(cb)

        #project a cube
        header=''
        fov=project_cube(cb,[4800,8000],False)
        
        #load the wcs
        w = wcs.WCS(fov['header'])
        
        raoffcurrent=[]
        decoffcurrent=[]
        
        #loop over reference stars 
        for sr in range(len(wcstab)):

            #go to pixels 
            pix = w.wcs_world2pix(wcstab['col1'][sr],wcstab['col2'][sr],0,1)
            #centroid
            ccd=PyGuide.CCDInfo(0,0,1)
            centroid=PyGuide.centroid(fov["img"],None,None,[pix[0],pix[1]],20,ccd)
            #back to ra/dec
            if centroid.isOK: 
                coord = w.wcs_pix2world(centroid.xyCtr[0],centroid.xyCtr[1],0,1)
                raoffcurrent.append(coord[0]-wcstab['col1'][sr])
                decoffcurrent.append(coord[1]-wcstab['col2'][sr])
                
        #stack offsets
        thisdra=np.mean(np.array(raoffcurrent))
        thisdde=np.mean(np.array(decoffcurrent))
        raoff.append(thisdra)
        decoff.append(thisdde)

        print "Offsets for this cube RA: {} Dec: {}".format(thisdra,thisdde)
                                                            

    #print offsets
    print "All Delta RA ", raoff
    print "All Delta Dec ", decoff

    rastring = ','.join(str(e) for e in raoff)
    decstring = ','.join(str(e) for e in decoff)

    #now prepare script for stacking and run the pipeline 
    sof=open("combine.sof","w")
    for pxtb in pixtab:
        sof.write(pxtb+' PIXTABLE_REDUCED\n') 
    sof.write('filter_list.fits FILTER_LIST') 
    sof.close()

    #Write the command file 
    scr=open("make_combine.sh","w")
    scr.write('MUSE_XCOMBINE_RA_OFFSETS="'+rastring+'"\n')
    scr.write('MUSE_XCOMBINE_DEC_OFFSETS="'+decstring+'"\n')
    scr.write("esorex --log-file=exp_combine.log muse_exp_combine --lambdamin={0:f} --lambdamax={1:f} combine.sof".format(lmin,lmax))
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "make_combine.sh"])
    
    return



