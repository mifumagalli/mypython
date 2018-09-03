#global import
import numpy as np
import glob 
import subprocess
import time
import os 
from distutils.spawn import find_executable
import socket

def grabtime(namelist):

    namelist=list(namelist)
    timelist=[]
    for name in namelist:
        #grab time sequence
        tmstr=name.split("MUSE.")[1]
        tm=time.strptime(tmstr.split(".")[0],'%Y-%m-%dT%H:%M:%S')
        timelist.append(time.mktime(tm))

    #return time in seconds
    return timelist
        
def parse_xml(path='./',nproc=12,pipecal=False):
    
    import glob
    from astropy.io import fits

    #Parse the xml file to extract all that is needed
    xmllist=glob.glob("{0}/Raw/*.xml".format(path))
    
    if (len(xmllist) < 1):
        print("xml file not found!")
    elif(len(xmllist) > 1):
        print("I found multiple xml files! Merging them...")

    #prepare space
    xml_info={}
    
    #loop over files
    for thisxml in xmllist:
        #open this one
        xml = open(thisxml, 'r')
        #loop over files
        for line in xml:
            #mark valid data
            if("<file category=" in line):
                tmptype=(line.split('"'))[1]
                tmpname=(line.split('"'))[3]
                #if new kew make a new set:
                if tmptype not in xml_info.keys():
                    xml_info[tmptype]=set()
                #next add to the set only if file exists
                fileexist=glob.glob("{0}/Raw/{1}*".format(path,tmpname))
                if(len(fileexist)>0):
                    xml_info[tmptype].add(tmpname)                   
        xml.close()

    print('xml files loaded')

    #now try to pick the best calib files available based on matches of dates etc...
 
    #start by setting the time of the observations
    time_obj=grabtime(xml_info['OBJECT'])
    reference_time=time_obj[0]
 
    #check max time lag between observations
    if(len(xml_info['OBJECT']) > 1):
        delta_time=np.max(abs(np.roll(time_obj,1)-time_obj)/3600.)
        print('Object observations are taken {0} h apart'.format(delta_time))
        if(delta_time > 4):
            print('Large time lag bewteen observations! Check!')

    #now loop over calibrations and grab the best 
    allkey=xml_info.keys()
    for kk in allkey:
        if('OBJECT' not in kk):
            #find list
            currentlist=list(xml_info[kk])                                
                
            #find time offset in seconds
            times=np.array(grabtime(currentlist))
            delta_time=np.abs(times-reference_time)/3600.
            currentlist=np.array(currentlist)
            #now handle by keyword according to calibration plan
            if((kk == 'ARC') or (kk == 'BIAS') or (kk == 'FLAT')):
                #grab daily calibrations 
                recent=np.where(delta_time <= 14.)
                xml_info[kk]=currentlist[recent[0]]
                print 'Found {0} {1} taken within 1 day'.format(len(recent[0]),kk)
            elif((kk == 'SKYFLAT') or (kk == 'DARK')):
                #grab within 20 days
                recent=np.where(delta_time <= 20*12.)
                xml_info[kk]=currentlist[recent[0]]
                print 'Found {0} {1} taken within 20 days'.format(len(recent[0]),kk)
            #This is when you want only the best one
            else:
                #pick closest
                mintm=np.argmin(delta_time)
                xml_info[kk]=[currentlist[mintm]]
                print 'Best {0} taken within {1} days'.format(kk,delta_time[mintm]/24.)

    #set the calibration path relative and suffix
    xml_info["PATHCAL"]='../../Raw/'
    xml_info["SUFFIXCAL"]='.fits'



    #grab some info on executable dir
    hostname=socket.gethostname()
    if('mad' in hostname):
        #find version of pipe loaded
        output=os.popen("esorex --man-page muse_bias").readlines()
        for i in output:
            if('muse_bias -- version' in i):
                pipeversion=i.split(" ")[-1].strip()
        staticalpath='/usr/share/esopipes/datastatic/muse-'+pipeversion+'/'
    elif('zwicky' in hostname):
        esorexpath=find_executable('esorex')
        staticalpath=esorexpath.split('/bin')[0]
        pipeversion=staticalpath.split('/')[-1]
        #handle special case of pipe version 2.1.1
        if('2.1.1-1' in pipeversion):
            pipeversion='muse-2.1.1' 
        staticalpath=staticalpath+'/calib/'+pipeversion+'/cal/'
    elif('cosma' in hostname):
        esorexpath=find_executable('esorex')
        staticalpath=esorexpath.split('/bin')[0]
        pipeversion=staticalpath.split('/')[-1]
        staticalpath='/cosma/local/muse/'+pipeversion+'/calib/muse-'+pipeversion+'/cal/'
    else:
        print('Please specify location of static calibrations for {}'.format(hostname))


    #Here sort out things with static calibrations: GEOMETRY & ASTROMETRY 
    #This is the largest time at one should worry about legacy products 
    legacy_time=time.mktime(time.strptime("14 Feb 16", "%d %b %y"))       

    if(reference_time < legacy_time):
        print('Using legacy static calibrations')
        #pick the right one
        #Static cal are assume dto be top levl dir in reduction folder
        tedge1=time.mktime(time.strptime("01 Dec 14", "%d %b %y"))
        tedge2=time.mktime(time.strptime("15 Apr 15", "%d %b %y"))
        tedge3=time.mktime(time.strptime("09 Sep 15", "%d %b %y"))
        if(reference_time <= tedge1):
            #use commissioning static - you may need even older ones so check 
            geometrystatic='../../staticcal/geometry_table_wfm_comm2b.fits'   
            astrostatic='../../staticcal/astrometry_wcs_wfm_comm2b.fits'
        elif((reference_time > tedge1) & (reference_time <= tedge2)):
            geometrystatic='../../staticcal/geometry_table_wfm_2014-12-01.fits'
            astrostatic='../../staticcal/astrometry_wcs_wfm_2014-12-01.fits'
        elif((reference_time > tedge2) & (reference_time <= tedge3)):
            geometrystatic='../../staticcal/geometry_table_wfm_2015-04-16.fits'
            astrostatic='../../staticcal/astrometry_wcs_wfm_2015-04-16.fits'
        else:
            geometrystatic='../../staticcal/geometry_table_wfm_2015-09-10.fits'
            astrostatic='../../staticcal/astrometry_wcs_wfm_2015-09-10.fits'

    else:
        print('Using pipeline static calibrations for astrometry and geometry')
        astrostatic=staticalpath+'astrometry_wcs_wfm.fits'
        geometrystatic=staticalpath+'geometry_table_wfm.fits'
    
    #update from default
    xml_info['ASTROMETRY_WCS']=[astrostatic]
    xml_info['GEOMETRY_TABLE']=[geometrystatic]


    #if so desired, force the use of calibrations provided in the pipeline package
    if(pipecal):
        print('Switch to static calibrations provided by pipeline!')
        #reset calibration path and suffix to absolute
        xml_info["PATHCAL"]=''
        xml_info["SUFFIXCAL"]=''

        #update
        xml_info['ASTROMETRY_WCS']=[staticalpath+'astrometry_wcs_wfm.fits']
        xml_info['VIGNETTING_MASK']=[staticalpath+'vignetting_mask.fits']
        xml_info['SKY_LINES']=[staticalpath+'sky_lines.fits']
        xml_info['ASTROMETRY_REFERENCE']=[staticalpath+'astrometry_reference.fits']
        xml_info['EXTINCT_TABLE']=[staticalpath+'extinct_table.fits']
        xml_info['GEOMETRY_TABLE']=[staticalpath+'geometry_table_wfm.fits']
        xml_info['FILTER_LIST']=[staticalpath+'filter_list.fits']
        xml_info['STD_FLUX_TABLE']=[staticalpath+'std_flux_table.fits']
        xml_info['BADPIX_TABLE']=[staticalpath+'badpix_table.fits']
        xml_info['LINE_CATALOG']=[staticalpath+'line_catalog.fits']
        #find out mode used for science - up to you not to mix them up!
        #ESO-INS-MODE: WFM-AO-E WFM-AO-N WFM-NOAO-E WFM-NOAO-N
        objheader=fits.open("Raw/"+(list(xml_info['OBJECT'])[0])+".fits.fz")
        mode=(objheader[0].header['ESO INS MODE']).strip()
        if(('WFM-AO-N' in mode)|('WFM-NOAO-N' in mode)):
            xml_info['LSF_PROFILE']=[staticalpath+'lsf_profile_slow_wfm-n.fits']
        elif(('WFM-AO-E' in mode)|('WFM-NOAO-E' in mode)):
            xml_info['LSF_PROFILE']=[staticalpath+'lsf_profile_slow_wfm-e.fits']
        else:
            raise ValueError("Instrument mode {} is not know for calibration selection".format(mode))
        
    #now perform healthy checks for calibrations with AO mode
    objheader=fits.open("Raw/"+(list(xml_info['OBJECT'])[0])+".fits.fz")
    scimode=(objheader[0].header['ESO INS MODE']).strip()
    if('WFM-AO-' in scimode):
        #loop over standards and check
        for std in list(xml_info['STD']):
            objheader=fits.open("Raw/"+std+".fits.fz")
            thismode=(objheader[0].header['ESO INS MODE']).strip()
            if(thismode not in scimode):
                raise ValueError("AO objects require AO standards! Using {} and {}".format(scimode,thismode))
        #now loop over flats
        for flt in list(xml_info['FLAT']):
            objheader=fits.open("Raw/"+flt+".fits.fz")
            thismode=(objheader[0].header['ESO INS MODE']).strip()
            if(thismode not in scimode):
                raise ValueError("AO objects require AO flats! Using {} and {}".format(scimode,thismode))
            
    #now dump cal plan
    print('Writing calibration plan in calplan.txt')
    cl=open('calplan.txt','w')
    for kk in xml_info.keys():
        for ll in xml_info[kk]:
            if(('PATHCAL' not in kk) & ('SUFFIXCAL' not in kk)):
                cl.write('{0} {1}\n'.format(kk,ll))
    cl.close()

    return xml_info  


def make_bias(xml_info,nproc=12):
    
    #grab the bias
    bias_list=xml_info["BIAS"]
    
    #Write the sof file 
    sof=open("../../Script/bias.sof","w")
    for ii in bias_list:
        sof.write("../../Raw/{0}.fits.fz BIAS\n".format(ii)) 
    sof.close()
    
    #Write the command file 
    scr=open("../../Script/make_bias.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=bias.log muse_bias --nifu=-1 --merge ../../Script/bias.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../../Script/make_bias.sh"])
    

def make_dark(xml_info,nproc=12):
    
    #grab the dark
    dark_list=xml_info["DARK"]

    #handle nicels case of no dark
    if(len(dark_list) > 0):

        #Write the sof file 
        sof=open("../../Script/dark.sof","w")
        for ii in dark_list:
            sof.write("../../Raw/{0}.fits.fz DARK\n".format(ii)) 
        sof.write("MASTER_BIAS.fits MASTER_BIAS\n")        
        sof.close()

        #Write the command file 
        scr=open("../../Script/make_dark.sh","w")
        scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
        scr.write("esorex --log-file=dark.log muse_dark --nifu=-1 --merge ../../Script/dark.sof")
        scr.close()

        #Run pipeline 
        subprocess.call(["sh", "../../Script/make_dark.sh"])
    
    else:
        print ("No DARK found... return!")
        return

def make_flat(xml_info,nproc=12):

    #grab the flat
    flat_list=xml_info["FLAT"]
    pix_tab=xml_info["BADPIX_TABLE"][0]
    
    #Write the sof file 
    sof=open("../../Script/flat.sof","w")
    for ii in flat_list:
        sof.write("../../Raw/{0}.fits.fz FLAT\n".format(ii)) 
    sof.write("{}{}{} BADPIX_TABLE\n".format(xml_info["PATHCAL"],pix_tab,xml_info["SUFFIXCAL"])) 
    sof.write("MASTER_BIAS.fits MASTER_BIAS\n")        
    sof.close()
              
    #Write the command file 
    scr=open("../../Script/make_flat.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=flat.log muse_flat --nifu=-1 --merge ../../Script/flat.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../../Script/make_flat.sh"])

def make_arcs(xml_info,nproc=12):
    
    #grab the arcs
    arc_list=xml_info["ARC"]
    line_cat=xml_info["LINE_CATALOG"][0]

    #Write the sof file 
    sof=open("../../Script/wavecal.sof","w")
    for ii in arc_list:
        sof.write("../../Raw/{0}.fits.fz ARC\n".format(ii))
    sof.write("{}{}{} LINE_CATALOG\n".format(xml_info["PATHCAL"],line_cat,xml_info["SUFFIXCAL"])) 
    sof.write("MASTER_BIAS.fits MASTER_BIAS\n") 
    sof.write("TRACE_TABLE.fits TRACE_TABLE\n") 
    sof.close()
    
    #Write the command file 
    scr=open("../../Script/make_wavecal.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=wavecal.log muse_wavecal --nifu=-1 --resample --residuals --merge ../../Script/wavecal.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../../Script/make_wavecal.sh"])


def make_twiflat(xml_info,nproc=12):

    #grab the twiflat
    flat_list=xml_info["SKYFLAT"]
    geom_cat=xml_info["GEOMETRY_TABLE"][0]

    #check time stamp for vignetting
    time_flat=grabtime(flat_list)
    vignetting_cat=xml_info["VIGNETTING_MASK"][0]
  
    #Write the sof file 
    sof=open("../../Script/twilight.sof","w")
    for ii in flat_list:
        sof.write("../../Raw/{0}.fits.fz SKYFLAT\n".format(ii))
    sof.write("{0} GEOMETRY_TABLE\n".format(geom_cat)) 
    sof.write("MASTER_BIAS.fits MASTER_BIAS\n") 
    sof.write("MASTER_FLAT.fits MASTER_FLAT\n") 
    sof.write("TRACE_TABLE.fits TRACE_TABLE\n") 
    sof.write("WAVECAL_TABLE.fits WAVECAL_TABLE\n")

    #add vignetting as appropriate before March 2017
    legacy_time=time.mktime(time.strptime("11 Mar 17", "%d %b %y"))       
    if(time_flat[0] < legacy_time):
        sof.write("../../Raw/{0}.fits VIGNETTING_MASK\n".format(vignetting_cat)) 

    sof.close()

    #Write the command file 
    scr=open("../../Script/make_twilight.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=twilight.log muse_twilight ../../Script/twilight.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../../Script/make_twilight.sh"])



def make_stdstar(xml_info,nproc=12):

    #grab the std
    std_list=xml_info["STD"][0]
    geom_cat=xml_info["GEOMETRY_TABLE"][0]
    pix_tab=xml_info["BADPIX_TABLE"][0]
        
    sof=open("../../Script/object_std.sof","w")
    sof.write("../../Raw/{0}.fits.fz STD\n".format(std_list)) 
    sof.write("{0} GEOMETRY_TABLE\n".format(geom_cat)) 
    sof.write("{}{}{} BADPIX_TABLE\n".format(xml_info["PATHCAL"],pix_tab,xml_info["SUFFIXCAL"])) 
    sof.write("MASTER_BIAS.fits MASTER_BIAS\n") 
    sof.write("MASTER_FLAT.fits MASTER_FLAT\n") 
    sof.write("TRACE_TABLE.fits TRACE_TABLE\n") 
    sof.write("WAVECAL_TABLE.fits WAVECAL_TABLE\n") 
    sof.write("TWILIGHT_CUBE.fits TWILIGHT_CUBE\n") 
    sof.close()

    #Write the command file 
    scr=open("../../Script/make_scibasic_std.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=object_std.log muse_scibasic --nifu=-1 --merge ../../Script/object_std.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../../Script/make_scibasic_std.sh"])

def make_stdflux(xml_info,nproc=12):

    #grab the table
    ext_tab=xml_info["EXTINCT_TABLE"][0]
    flx_tab=xml_info["STD_FLUX_TABLE"][0]
    
    #Write the sof file 
    sof=open("../../Script/std.sof","w")
    sof.write("{}{}{} EXTINCT_TABLE\n".format(xml_info["PATHCAL"],ext_tab,xml_info["SUFFIXCAL"]))
    sof.write("{}{}{} STD_FLUX_TABLE\n".format(xml_info["PATHCAL"],flx_tab,xml_info["SUFFIXCAL"])) 
    for ifu in range(24):
        sof.write("PIXTABLE_STD_0001-{0:02d}.fits PIXTABLE_STD\n".format(ifu+1)) 
    sof.close()

    #Write the command file 
    scr=open("../../Script/make_std.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=std.log muse_standard  --filter=white ../../Script/std.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../../Script/make_std.sh"])


def make_objects(xml_info,nproc=12):

    #grab the objects
    obj_list=list(xml_info["OBJECT"])
    #sort by time - useful is offsets fields available
    obj_list.sort()
    
    ill=xml_info["ILLUM"][0]
    geom_cat=xml_info["GEOMETRY_TABLE"][0]
    pix_tab=xml_info["BADPIX_TABLE"][0]

    #Write the sof file 
    sof=open("../../Script/object.sof","w")
       
    for ii in obj_list:
        sof.write("../../Raw/{0}.fits.fz OBJECT\n".format(ii))
    sof.write("../../Raw/{0}.fits.fz ILLUM\n".format(ill))
    sof.write("{0} GEOMETRY_TABLE\n".format(geom_cat)) 
    sof.write("{}{}{} BADPIX_TABLE\n".format(xml_info["PATHCAL"],pix_tab,xml_info["SUFFIXCAL"])) 
    sof.write("MASTER_BIAS.fits MASTER_BIAS\n") 
    sof.write("MASTER_FLAT.fits MASTER_FLAT\n") 
    sof.write("TRACE_TABLE.fits TRACE_TABLE\n") 
    sof.write("WAVECAL_TABLE.fits WAVECAL_TABLE\n") 
    sof.write("TWILIGHT_CUBE.fits TWILIGHT_CUBE\n") 
    sof.close()

    #Write the command file 
    scr=open("../../Script/make_scibasic.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=object.log muse_scibasic --nifu=-1 --merge ../../Script/object.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../../Script/make_scibasic.sh"])


def make_cubes(xml_info,nproc=12,wcsoff=None,refcube=None,scilist=None):

    import shutil
    from astropy.io import fits

    #Now process the science

    #if scilist not provided, create manaully
    if not (scilist):
        #Search how many exposures are there
        scils=glob.glob("OBJECT_RED_0*.fits*")
        scilist = np.arange(len(scils))+1

    #how many exposures
    nsci=len(scilist)

    #create suffixes as needed
    suffix = ''
    if(refcube):
        print 'Using external WCS structure for cube output'
        suffix += '_pos'
    if(wcsoff):
        suffix += '_off'  
    
    #loop on exposures and reduce frame without sky subtraction 
    for exp in range(nsci):

        #reconstruct names
        cname="DATACUBE_FINAL_EXP{0:d}{1}.fits".format(scilist[exp],suffix)
        pname="PIXTABLE_REDUCED_EXP{0:d}{1}.fits".format(scilist[exp],suffix)
        iname="IMAGE_FOV_EXP{0:d}{1}.fits".format(scilist[exp],suffix)
        
        if not os.path.isfile(cname):

            print "Processing exposure {0:d}".format(exp+1)
            
            #Write the sof file 
            sof=open("../../Script/scipost_{0:d}.sof".format(scilist[exp]),"w")
            sof.write("{0} ASTROMETRY_WCS\n".format(xml_info["ASTROMETRY_WCS"][0])) 
            sof.write("{}{}{} SKY_LINES\n".format(xml_info["PATHCAL"],xml_info["SKY_LINES"][0],xml_info["SUFFIXCAL"])) 
            sof.write("{}{}{} EXTINCT_TABLE\n".format(xml_info["PATHCAL"],xml_info["EXTINCT_TABLE"][0],xml_info["SUFFIXCAL"])) 
            sof.write("{}{}{} FILTER_LIST\n".format(xml_info["PATHCAL"],xml_info["FILTER_LIST"][0],xml_info["SUFFIXCAL"])) 
            sof.write("{}{}{} LSF_PROFILE\n".format(xml_info["PATHCAL"],xml_info["LSF_PROFILE"][0],xml_info["SUFFIXCAL"])) 

            #if using reference set it
            if(refcube):
                sof.write("{0} OUTPUT_WCS\n".format(refcube)) 

            if(wcsoff):
	        sof.write("OFFSET_LIST.fits OFFSET_LIST\n") 
	    
	    for ifu in range(24):
		ifupixtab="PIXTABLE_OBJECT_{0:04d}-{1:02d}.fits".format(scilist[exp],ifu+1)
                #now write the pix tab in sof
                sof.write("{} PIXTABLE_OBJECT\n".format(ifupixtab)) 
                
            #finish writing sof
            sof.write("STD_RESPONSE_0001.fits STD_RESPONSE\n")
            sof.write("STD_TELLURIC_0001.fits STD_TELLURIC\n")
            sof.close()

            #Write the command file 
            scr=open("../../Script/make_scipost_{0:d}.sh".format(scilist[exp]),"w")
            scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 

            scr.write('esorex --log-file=scipost_{0:d}.log muse_scipost --skymethod="none" --filter=white --save=cube,individual ../../Script/scipost_{0:d}.sof'.format(scilist[exp]))
            scr.close()
    
            #Run pipeline 
            subprocess.call(["sh", "../../Script/make_scipost_{0:d}.sh".format(scilist[exp])])    
            subprocess.call(["mv","DATACUBE_FINAL.fits",cname])
            subprocess.call(["mv","IMAGE_FOV_0001.fits",iname])
            subprocess.call(["mv","PIXTABLE_REDUCED_0001.fits",pname])


def make_skymodel(xml_info,nproc=12):


    """

    Check if sky exposures are present - if so, handle them before moving to the 
    next stage 
    
    """

    from astropy.io import fits
    import glob
    import os
    
    #grab object list and build a objs time dictionary 
    allobjs=glob.glob('OBJECT_RED*')
    timedic={}
    for obj in allobjs:
        tagobj=obj.split("RED_")[1].split(".fits")[0]
        objheader=fits.open(obj)
        time=(objheader[0].header['DATE-OBS']).strip()
        timedic[time]=tagobj
    
        
    #now loop and find sky 
    for exp in list(xml_info['OBJECT']):
        objheader=fits.open("../../Raw/"+exp+".fits.fz")
        objorsky=(objheader[0].header['HIERARCH ESO DPR TYPE']).strip()
    
        #process sky exposure if exists
        if('SKY' in objorsky):
            #first identify the ID of the sky exposures
            skytime=(objheader[0].header['DATE-OBS']).strip()
            currentid=timedic[skytime]
          
            #run if needed
            if not os.path.isfile("IMAGE_FOV_{}.fits".format(currentid)):
  
                #Write the sof file 
                sof=open("../../Script/sky_{}.sof".format(currentid),"w")
                for ii in range(24):
                    sof.write("PIXTABLE_OBJECT_{0}-{1:02d}.fits PIXTABLE_SKY\n".format(currentid,ii+1)) 
                sof.write("STD_RESPONSE_0001.fits STD_RESPONSE\n")
                sof.write("STD_TELLURIC_0001.fits STD_TELLURIC\n")
                sof.write("{}{}{} LSF_PROFILE\n".format(xml_info["PATHCAL"],xml_info["LSF_PROFILE"][0],xml_info["SUFFIXCAL"]))
                sof.write("{}{}{} SKY_LINES\n".format(xml_info["PATHCAL"],xml_info["SKY_LINES"][0],xml_info["SUFFIXCAL"]))
                sof.write("{}{}{} EXTINCT_TABLE\n".format(xml_info["PATHCAL"],xml_info["EXTINCT_TABLE"][0],xml_info["SUFFIXCAL"]))
                sof.close()

                #Write the command file 
                scr=open("../../Script/make_sky_{}.sh".format(currentid),"w")
                scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
                scr.write("esorex --log-file=sky_{}.log muse_create_sky ../../Script/sky_{}.sof".format(currentid,currentid))
                scr.close()

                #Run pipeline 
                subprocess.call(["sh", "../../Script/make_sky_{}.sh".format(currentid)])

                #Rename outputs
                subprocess.call(["mv","IMAGE_FOV.fits","IMAGE_FOV_{}.fits".format(currentid)])
                subprocess.call(["mv","SKY_CONTINUUM.fits","SKY_CONTINUUM_{}.fits".format(currentid)])
                subprocess.call(["mv","SKY_LINES.fits","SKY_LINES_{}.fits".format(currentid)])
                subprocess.call(["mv","SKY_MASK.fits","SKY_MASK_{}.fits".format(currentid)])
                subprocess.call(["mv","SKY_SPECTRUM.fits","SKY_SPECTRUM_{}.fits".format(currentid)])

            
    print('All done with sky models!')
    
            
