
class Muse(object):

    """
    This is a class that bundles some of the procedures to handle muse data

    """

    def __init__(self):

        """ Stuff at init .. not much for now """
        
        print "Initialisation of Muse object"
        self.musepip=1.0
        

    def redux_basic(self,path='./',nproc=12,pipecal=False):
        
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

        pipecal - if set to true, static calibrations provided with the pipeline
                  will be used. This pplies to ALL static calibrations

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
        xml_info=rdx.parse_xml(path=path,nproc=nproc,pipecal=pipecal)
        
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
        print('Preparing intermediate data cubes...')
        rdx.make_cubes(xml_info,nproc=nproc)
        
        #In the end, handle sky offsets if present
        print('Checking if sky offsets are present and preparing sky model')
        rdx.make_skymodel(xml_info,nproc=nproc)

        #Done - back to original directory!
        print('All done with basic redux...')
        os.chdir(currdir)
        
        return xml_info

    def cubex_process(self,refpath='./esocombine/',skymask=None):

        """
  
        Take individual cubes and process them with CubEx to tweak the illumination and perform sky subtraction
        CubEx is a privite code by S. Cantalupo and cannot be redistributed. Contact the author directly. 
        
        refpath -> where the reference cubes for wcs resempling are 
        
        skymask -> mask this region before running cubesharp (ds9 region in image units)

        """

        import os
        import glob
        import subprocess
        import muse_redux_line as ex 
        import muse_redux_cubex as cx 
        import multiprocessing
        import numpy as np

        #first, list how many OBs are there
        listob=glob.glob('OB*')
        listob.sort()
        nobs=len(listob)
        print('Process {} OBs'.format(nobs))
    
        #now make space as needed for final products
        if not os.path.exists('cubexcombine'):
            os.makedirs('cubexcombine')

        #rerun pipeline enabling resampling on final ESO cube using modules coded for line_process
        ex.individual_resample(listob,refpath=refpath)
        
        #now do the first two passes of cubex on each OB to prepare a temporary cube
        cx.cubex_driver(listob,skymask=skymask)
       
        #prepare for intermediate combine 
        #dump to disk file lists
        topdir=os.getcwd()
        os.chdir('cubexcombine')

        fl1=open('cubes.lst','w')
        fl2=open('masks.lst','w')

        #loop over OBs
        for oob in range(nobs):
            #count how many science exposures
            nsci=len(glob.glob("../{}/Proc/OBJECT_RED_0*.fits*".format(listob[oob])))
            #reconstruct names 
            for ll in range(nsci):
                fl1.write('../{}/Proc/DATACUBE_FINAL_LINEWCS_EXP{}_skysub2.fits\n'.format(listob[oob],ll+1))
                fl2.write('../{}/Proc/DATACUBE_FINAL_LINEWCS_EXP{}_fix2_SliceEdgeMask.fits\n'.format(listob[oob],ll+1))
        fl1.close()
        fl2.close()
        
        #make the temp combine
        cx.combine_cubes("cubes.lst","masks.lst")

        #now make two independent halves 
        fl1cube=open('cubes_half1.lst','w')
        fl1mask=open('masks_half1.lst','w')
        fl2cube=open('cubes_half2.lst','w')
        fl2mask=open('masks_half2.lst','w')

        #loop over OBs
        counter=0
        for oob in range(nobs):
            #count how many science exposures
            nsci=len(glob.glob("../{}/Proc/OBJECT_RED_0*.fits*".format(listob[oob])))
            #reconstruct names 
            for ll in range(nsci):
                counter=counter+1
                if(counter % 2 == 0):
                    fl1cube.write('../{}/Proc/DATACUBE_FINAL_LINEWCS_EXP{}_skysub2.fits\n'.format(listob[oob],ll+1))
                    fl1mask.write('../{}/Proc/DATACUBE_FINAL_LINEWCS_EXP{}_fix2_SliceEdgeMask.fits\n'.format(listob[oob],ll+1))
                else:
                    fl2cube.write('../{}/Proc/DATACUBE_FINAL_LINEWCS_EXP{}_skysub2.fits\n'.format(listob[oob],ll+1))
                    fl2mask.write('../{}/Proc/DATACUBE_FINAL_LINEWCS_EXP{}_fix2_SliceEdgeMask.fits\n'.format(listob[oob],ll+1))
                    
        #close files
        fl1cube.close()
        fl1mask.close()
        fl2cube.close()      
        fl2mask.close()

        #now combine
        cx.combine_cubes("cubes_half1.lst","masks_half1.lst",halfset='half1')
        cx.combine_cubes("cubes_half2.lst","masks_half2.lst",halfset='half2')
        
        #back to top 
        os.chdir(topdir)

        #now do the final pass of cubex using the tmp combined cube for better masking
        cx.cubex_driver(listob,last=True,highsn='../../cubexcombine/COMBINED_CUBE.fits',skymask=skymask)

        #make the final combined cube
        #dump to disk file lists
        topdir=os.getcwd()
        os.chdir('cubexcombine')

        fl1=open('cubes_final.lst','w')
        fl2=open('masks_final.lst','w')

        #loop over OBs
        for oob in range(nobs):
            #count how many science exposures
            nsci=len(glob.glob("../{}/Proc/OBJECT_RED_0*.fits*".format(listob[oob])))
            #reconstruct names 
            for ll in range(nsci):
                fl1.write('../{}/Proc/DATACUBE_FINAL_LINEWCS_EXP{}_skysubhsn.fits\n'.format(listob[oob],ll+1))
                fl2.write('../{}/Proc/DATACUBE_FINAL_LINEWCS_EXP{}_fixhsn_SliceEdgeMask.fits\n'.format(listob[oob],ll+1))
        fl1.close()
        fl2.close()
        
        #make the temp combine
        cx.combine_cubes("cubes_final.lst","masks_final.lst",final=True)


        #now make two independent halves 
        fl1cube=open('cubes_final_half1.lst','w')
        fl1mask=open('masks_final_half1.lst','w')
        fl2cube=open('cubes_final_half2.lst','w')
        fl2mask=open('masks_final_half2.lst','w')

        #loop over OBs
        counter=0
        for oob in range(nobs):
            #count how many science exposures
            nsci=len(glob.glob("../{}/Proc/OBJECT_RED_0*.fits*".format(listob[oob])))
            #reconstruct names 
            for ll in range(nsci):
                counter=counter+1
                if(counter % 2 == 0):
                    fl1cube.write('../{}/Proc/DATACUBE_FINAL_LINEWCS_EXP{}_skysubhsn.fits\n'.format(listob[oob],ll+1))
                    fl1mask.write('../{}/Proc/DATACUBE_FINAL_LINEWCS_EXP{}_fixhsn_SliceEdgeMask.fits\n'.format(listob[oob],ll+1))
                else:
                    fl2cube.write('../{}/Proc/DATACUBE_FINAL_LINEWCS_EXP{}_skysubhsn.fits\n'.format(listob[oob],ll+1))
                    fl2mask.write('../{}/Proc/DATACUBE_FINAL_LINEWCS_EXP{}_fixhsn_SliceEdgeMask.fits\n'.format(listob[oob],ll+1))
                    
        #close files
        fl1cube.close()
        fl1mask.close()
        fl2cube.close()      
        fl2mask.close()

        #now combine
        cx.combine_cubes("cubes_final_half1.lst","masks_final_half1.lst",halfsetfinal='half1')
        cx.combine_cubes("cubes_final_half2.lst","masks_final_half2.lst",halfsetfinal='half2')

        #now run quality checks on final redux products
        #Typically one uses first pass, so check those
        cx.dataquality("cubes.lst","masks.lst")

        #back to top level
        os.chdir(topdir)
        print 'All done with cubex redux'
        

    def eso_process(self):

        """

        After running the basic reduction, this sequence generates 
        a skysubtracted and combined cube following only eso recepies 
             
        """

        import os
        import glob
        import subprocess
        import muse_redux_eso as ex 
        
        #first, list how many OBs are there
        listob=glob.glob('OB*')
        listob.sort()
        nobs=len(listob)
        print('Process {} OBs'.format(nobs))

        #rerun pipe enabling skysubtraction and 
        #dumping fully reduced pixel table 
        ex.individual_skysub(listob)
        
        #now make space as needed for final products
        if not os.path.exists('esocombine'):
            os.makedirs('esocombine')

        #change dir
        currdir=os.getcwd()
        os.chdir('esocombine')
        print('Changing dir to esocombine...')

        #coadd all obs after aligment 
        ex.coaddall(listob)

        #back to original place
        print('Back to top level...')
        os.chdir(currdir)
        print('All done!')
        

    def line_process(self,skymode='internal',refpath='./esocombine/',skymask=None,lmin=4900,lmax=9000,
                     deepwhite=None):

        """

        Produces final cubes optimised for fields that are relatively 
        empty in continuum sources but that may have very extended 
        emission lines.

        lmin -> the minimum wavelength to consider
        lmax -> the maximum wavelength to consider
        refpath -> where the reference cubes for wcs resempling are 
        skymask -> a skymask to be used for identify good regions for skysubtraction 
                   (expected in image coordinates)
        skymode -> internal: use good pixels (i.e. not containing sources or defined in skymask) to perform 
                   skysubtraction plus run ZAP on it.
        deepwhite -> if set to an image, this is used to mask sources during sky subtraction.
                     otherwise the cube itself is used 
                   
        """

        import os
        import glob
        import subprocess
        import muse_redux_line as ex 

        #first, list how many OBs are there
        listob=glob.glob('OB*')
        listob.sort()
        nobs=len(listob)
        print('Process {} OBs'.format(nobs))
        
    
        #now make space as needed for final products
        if not os.path.exists('linecombine'):
            os.makedirs('linecombine')

        #rerun pipe enabling resampling on final ESO cube
        ex.individual_resample(listob,refpath=refpath)
        
        #next construct the ifu mask for each exposure 
        ex.make_ifumasks(listob,refpath=refpath)

        #compute illumination correction 
        ex.make_illcorr(listob)
        
        #now do background subtraction
        if('internal' in skymode):
            ex.internalskysub(listob,skymask,deepwhite=deepwhite)
        else:
            print ("Sky subtraction mode {} not supported".format(skymode))
        
        #change dir
        currdir=os.getcwd()
        os.chdir('linecombine')
        print('Changing dir to linecombine...')
            
        fl1=open('cubes.lst','w')
        fl2=open('masks.lst','w')

        #loop over OBs
        for oob in range(nobs):
            #count how many science exposures
            nsci=len(glob.glob("../{}/Proc/OBJECT_RED_0*.fits*".format(listob[oob])))
            #reconstruct names 
            for ll in range(nsci):
                fl1.write('../{}/Proc/DATACUBE_FINAL_LINEWCS_EXP{}_zapsky.fits\n'.format(listob[oob],ll+1))
                fl2.write('../{}/Proc/MASK_EXP{}_ILLCORR_edges.fits\n'.format(listob[oob],ll+1))
        fl1.close()
        fl2.close()
        
        #make the temp combine
        ex.combine_cubes("cubes.lst","masks.lst")
        os.chdir(currdir)

        print("All done!")



