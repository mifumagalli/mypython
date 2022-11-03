
class Hawki(object):

    """
    This is a class that bundles some of the procedures to handle hawki data

    """

    def __init__(self):

        """ Stuff at init .. not much for now """
        
        print("Initialisation of HAWKI object")
        self.hawkipip=1.0
        

    def redux_basic(self,path='./',nproc=12,skyalgo='pawsky_mask',astrocat='wise',pipecal=False):
        
        """ 

        This is a basic engine that performs redux of HAWKI data using the eso pipeline in a basic 
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
        
        from . import hawki_redux_basic as rdx
        import os

        print('Starting reduction...')
        
        #First, make sure the various folders exist as needed 
        if not os.path.exists(path+"Raw"):
            print("Cannot find Raw data...")
            exit()
        if not os.path.exists(path+"Script"):
            os.makedirs(path+"Script")
        if not os.path.exists(path+"Proc"):
            os.makedirs(path+"Proc")


        #parse the xml file(s) 
        xml_info=rdx.parse_xml(path=path,nproc=nproc,pipecal=pipecal)
        nobs = xml_info['NOBS']

        #now start reduction. Enter the proc folder
        currdir=os.getcwd()
        if not os.path.exists("./Calib"):
            os.makedirs("./Calib")
        os.chdir("./Calib")
        print('Changing dir to Calib...')
   
        #First handle the darks
        if not os.path.isfile("MASTER_DARK_OBJ.fits"):
            print('Creating dark for obj...')
            rdx.make_dark_obj(xml_info,nproc=nproc)
            print('All done with the dark obj...')
        else:
            print('Dark obj already exists')
            
        if not os.path.isfile("MASTER_DARK_FLT.fits"):
            print('Creating dark for flt...')
            rdx.make_dark_flt(xml_info,nproc=nproc)
            print('All done with the dark flt...')
        else:
            print('Dark flt already exists')
            
        #Next handle the flats
        if not os.path.isfile("MASTER_TWILIGHT_FLAT.fits"):
            print('Creating flat...')
            rdx.make_flat(xml_info,nproc=nproc)
            print('All done with flat...')
        else:
            print('Flat already exists')
       
        os.chdir('../Proc')
        print('Changing dir to proc...')
        
        for ob in range(1,nobs+1):
          
          if not os.path.exists("./OBmin{}".format(ob)):
              os.makedirs("./OBmin{}".format(ob))
          os.chdir("./OBmin{}".format(ob))
          print('Changing dir to OBmin{}...'.format(ob))
         
          #Next calibrate objects
          if not os.path.isfile("exp_1.fits"):
             print('Processing the objects...')
             rdx.make_objects(xml_info, ob, skyalgo=skyalgo,astrocat=astrocat,nproc=nproc)
             rdx.make_fluxcal(xml_info)
             print('All done with objects...')
          else:
            print('Objects already processed')
          
          #This will always be run, since it internally checks for file existence
          print('Second order sky subtraction...')
          rdx.make_skysub(xml_info)
          
          
          os.chdir('../')  
        
        #Done - back to original directory!
        print('All done with frames redux...')
        os.chdir(currdir)
        
        return xml_info

