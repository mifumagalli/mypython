#global import
import numpy as np
import glob 
import subprocess
import time
import os 
from distutils.spawn import find_executable

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
        
def parse_xml(path='./',nproc=12):
    
    import glob

    #Parse the xml file to extract all that is needed
    xmllist=glob.glob("{0}/Raw/*.xml".format(path))
    
    if (len(xmllist) < 1):
        print "xml file not found!"
    elif(len(xmllist) > 1):
        print "I found multiple xml files! Merging them..."

    #prepare space
    xml_info={}
    
    #loop over files
    for thisxml in xmllist:
        #open this one
        xml = open(thisxml, 'r')
        #loop over files
        for line in xml:
            #mark valid data
            if "<file category=" in line:
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

    print 'xml files loaded'

    #now try to pick the best calib files available based on matches of dates etc...
 
    #start by setting the time of the observations
    time_obj=grabtime(xml_info['OBJECT'])
    reference_time=time_obj[0]
 
    #check max time lag between observations
    if(len(xml_info['OBJECT']) > 1):
        delta_time=np.max(abs(np.roll(time_obj,1)-time_obj)/3600.)
        print 'Object observations are taken {0} h apart'.format(delta_time)
        if(delta_time > 4):
            print 'Large time lag bewteen observations! Check!'

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
            #now handle by keyword
            if((kk == 'SKYFLAT') or (kk == 'FLAT') or (kk == 'ARC') or (kk == 'DARK') or (kk == 'BIAS')):
            #This is when you allow for multiple
                #grab within 20 days
                recent=np.where(delta_time <= 20*24.)
                xml_info[kk]=currentlist[recent[0]]
                
                print 'Found {0} {1} taken within 20 days'.format(len(recent[0]),kk)
            #This is when you want only the best one
            else:
                #pick closest
                mintm=np.argmin(delta_time)
                xml_info[kk]=[currentlist[mintm]]
                print 'Best {0} taken within {1} days'.format(kk,delta_time[mintm]/24.)


    #here sort out things with static calibrations: GEOMETRY & ASTROMETRY 

    #This is the largest time at one should worry about legacy products 
    legacy_time=time.mktime(time.strptime("08 Sep 15", "%d %b %y"))       

    if(reference_time < legacy_time):
        print 'Using legacy static calibrations'
        #pick the right one
        tedge1=time.mktime(time.strptime("01 Dec 14", "%d %b %y"))
        tedge2=time.mktime(time.strptime("15 Apr 15", "%d %b %y"))
        if(reference_time <= tedge1):
            #use commissioning static - you may need even older ones so check 
            geometrystatic='staticcal/geometry_table_wfm_comm2b.fits'   
            astrostatic='staticcal/astrometry_wcs_wfm_comm2b.fits'
        elif((reference_time > tedge1) & (reference_time <= tedge2)):
            geometrystatic='staticcal/geometry_table_wfm_2014-12-01.fits'
            astrostatic='staticcal/astrometry_wcs_wfm_2014-12-01.fits'
        else:
            geometrystatic='staticcal/geometry_table_wfm_2015-04-16.fits'
            astrostatic='staticcal/astrometry_wcs_wfm_2015-04-16.fits'
    else:
        print 'Using pipeline static calibrations'
        esorexpath=find_executable('esorex')
        staticalpath=esorexpath.split('/bin')[0]
        pipeversion=staticalpath.split('/')[-1]
        staticalpath=staticalpath+'/calib/'+pipeversion+'/cal/'
        astrostatic=staticalpath+'astrometry_wcs_wfm.fits'
        geometrystatic=staticalpath+'geometry_table_wfm.fits'
        
    #update form default
    xml_info['ASTROMETRY_WCS']=[astrostatic]
    xml_info['GEOMETRY_TABLE']=[geometrystatic]

    #now dump cal plan
    print 'Writing calibration plan in calplan.txt'            
    cl=open('calplan.txt','w')
    for kk in xml_info.keys():
        for ll in xml_info[kk]:
            cl.write('{0} {1}\n'.format(kk,ll))
    cl.close()
    return xml_info  


def make_bias(xml_info,nproc=12):
    
    #grab the bias
    bias_list=xml_info["BIAS"]
    
    #Write the sof file 
    sof=open("../Script/bias.sof","w")
    for ii in bias_list:
        sof.write("../Raw/{0}.fits.fz BIAS\n".format(ii)) 
    sof.close()
    
    #Write the command file 
    scr=open("../Script/make_bias.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=bias.log muse_bias --nifu=-1 --merge ../Script/bias.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../Script/make_bias.sh"])
    

def make_dark(xml_info,nproc=12):
    
    #grab the dark
    dark_list=xml_info["DARK"]

    #Write the sof file 
    sof=open("../Script/dark.sof","w")
    for ii in dark_list:
        sof.write("../Raw/{0}.fits.fz DARK\n".format(ii)) 
    sof.write("MASTER_BIAS.fits MASTER_BIAS\n")        
    sof.close()

    #Write the command file 
    scr=open("../Script/make_dark.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=dark.log muse_dark --nifu=-1 --merge ../Script/dark.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../Script/make_dark.sh"])
   
def make_flat(xml_info,nproc=12):

    #grab the flat
    flat_list=xml_info["FLAT"]
    pix_tab=xml_info["BADPIX_TABLE"][0]
    
    #Write the sof file 
    sof=open("../Script/flat.sof","w")
    for ii in flat_list:
        sof.write("../Raw/{0}.fits.fz FLAT\n".format(ii)) 
    sof.write("../Raw/{0}.fits BADPIX_TABLE\n".format(pix_tab)) 
    sof.write("MASTER_BIAS.fits MASTER_BIAS\n")        
    sof.close()
              
    #Write the command file 
    scr=open("../Script/make_flat.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=flat.log muse_flat --nifu=-1 --merge ../Script/flat.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../Script/make_flat.sh"])

def make_arcs(xml_info,nproc=12):
    
    #grab the arcs
    arc_list=xml_info["ARC"]
    line_cat=xml_info["LINE_CATALOG"][0]

    #Write the sof file 
    sof=open("../Script/wavecal.sof","w")
    for ii in arc_list:
        sof.write("../Raw/{0}.fits.fz ARC\n".format(ii))
    sof.write("../Raw/{0}.fits LINE_CATALOG\n".format(line_cat)) 
    sof.write("MASTER_BIAS.fits MASTER_BIAS\n") 
    sof.write("TRACE_TABLE.fits TRACE_TABLE\n") 
    sof.close()
    
    #Write the command file 
    scr=open("../Script/make_wavecal.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=wavecal.log muse_wavecal --nifu=-1 --resample --residuals --merge ../Script/wavecal.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../Script/make_wavecal.sh"])

def make_twiflat(xml_info,nproc=12):

    #grab the twiflat
    flat_list=xml_info["SKYFLAT"]
    geom_cat=xml_info["GEOMETRY_TABLE"][0]
    
    #Write the sof file 
    sof=open("../Script/twilight.sof","w")
    for ii in flat_list:
        sof.write("../Raw/{0}.fits.fz SKYFLAT\n".format(ii))
    sof.write("../Raw/{0}.fits GEOMETRY_TABLE\n".format(geom_cat)) 
    sof.write("MASTER_BIAS.fits MASTER_BIAS\n") 
    sof.write("MASTER_FLAT.fits MASTER_FLAT\n") 
    sof.write("TRACE_TABLE.fits TRACE_TABLE\n") 
    sof.write("WAVECAL_TABLE.fits WAVECAL_TABLE\n") 
    sof.close()

    #Write the command file 
    scr=open("../Script/make_twilight.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=twilight.log muse_twilight ../Script/twilight.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../Script/make_twilight.sh"])

def make_stdstar(xml_info,nproc=12):

    #grab the std
    std_list=xml_info["STD"][0]
    geom_cat=xml_info["GEOMETRY_TABLE"][0]
    pix_tab=xml_info["BADPIX_TABLE"][0]
        
    sof=open("../Script/object_std.sof","w")
    sof.write("../Raw/{0}.fits.fz STD\n".format(std_list)) 
    sof.write("../Raw/{0}.fits GEOMETRY_TABLE\n".format(geom_cat)) 
    sof.write("../Raw/{0}.fits BADPIX_TABLE\n".format(pix_tab)) 
    sof.write("MASTER_BIAS.fits MASTER_BIAS\n") 
    sof.write("MASTER_FLAT.fits MASTER_FLAT\n") 
    sof.write("TRACE_TABLE.fits TRACE_TABLE\n") 
    sof.write("WAVECAL_TABLE.fits WAVECAL_TABLE\n") 
    sof.write("TWILIGHT_CUBE.fits TWILIGHT_CUBE\n") 
    sof.close()

    #Write the command file 
    scr=open("../Script/make_scibasic_std.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=object_std.log muse_scibasic --nifu=-1 --merge ../Script/object_std.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../Script/make_scibasic_std.sh"])

def make_stdflux(xml_info,nproc=12):

    #grab the table
    ext_tab=xml_info["EXTINCT_TABLE"][0]
    flx_tab=xml_info["STD_FLUX_TABLE"][0]
    
    #Write the sof file 
    sof=open("../Script/std.sof","w")
    sof.write("../Raw/{0}.fits EXTINCT_TABLE\n".format(ext_tab))
    sof.write("../Raw/{0}.fits STD_FLUX_TABLE\n".format(flx_tab)) 
    for ifu in range(24):
        sof.write("PIXTABLE_STD_0001-{0:02d}.fits PIXTABLE_STD\n".format(ifu+1)) 
    sof.close()

    #Write the command file 
    scr=open("../Script/make_std.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=std.log muse_standard  --filter=white ../Script/std.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../Script/make_std.sh"])


def make_objects(xml_info,nproc=12):

    #grab the objects
    obj_list=xml_info["OBJECT"]
    ill=xml_info["ILLUM"][0]
    geom_cat=xml_info["GEOMETRY_TABLE"][0]
    pix_tab=xml_info["BADPIX_TABLE"][0]

    #Write the sof file 
    sof=open("../Script/object.sof","w")
       
    for ii in obj_list:
        sof.write("../Raw/{0}.fits.fz OBJECT\n".format(ii))
    sof.write("../Raw/{0}.fits.fz ILLUM\n".format(ill))
    sof.write("../Raw/{0}.fits GEOMETRY_TABLE\n".format(geom_cat)) 
    sof.write("../Raw/{0}.fits BADPIX_TABLE\n".format(pix_tab)) 
    sof.write("MASTER_BIAS.fits MASTER_BIAS\n") 
    sof.write("MASTER_FLAT.fits MASTER_FLAT\n") 
    sof.write("TRACE_TABLE.fits TRACE_TABLE\n") 
    sof.write("WAVECAL_TABLE.fits WAVECAL_TABLE\n") 
    sof.write("TWILIGHT_CUBE.fits TWILIGHT_CUBE\n") 
    sof.close()

    #Write the command file 
    scr=open("../Script/make_scibasic.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=object.log muse_scibasic --nifu=-1 --merge ../Script/object.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../Script/make_scibasic.sh"])


def make_cubes(xml_info,nproc=12,wcsoff=None,refcube=None):

    import shutil
    from astropy.io import fits

    #Now process the science     
    #Search how many exposures are there
    scils=glob.glob("OBJECT_RED_0*.fits*")
    nsci=len(scils)

    #loop on exposures and reduce frame without sky subtraction 
    for exp in range(nsci):
        
        if(wcsoff):
            cname="DATACUBE_FINAL_EXP{0:d}_off.fits".format(exp+1)
            pname="PIXTABLE_REDUCED_EXP{0:d}_off.fits".format(exp+1)
            iname="IMAGE_FOV_EXP{0:d}_off.fits".format(exp+1)
        else:
            cname="DATACUBE_FINAL_EXP{0:d}.fits".format(exp+1)
            pname="PIXTABLE_REDUCED_EXP{0:d}.fits".format(exp+1)
            iname="IMAGE_FOV_EXP{0:d}.fits".format(exp+1)

        if not os.path.isfile(cname):

            print "Processing exposure {0:d}".format(exp+1)
            
            #Write the sof file 
            sof=open("../Script/scipost_{0:d}.sof".format(exp+1),"w")
            sof.write("../Raw/{0}.fits ASTROMETRY_WCS\n".format(xml_info["ASTROMETRY_WCS"][0])) 
            sof.write("../Raw/{0}.fits SKY_LINES\n".format(xml_info["SKY_LINES"][0])) 
            sof.write("../Raw/{0}.fits EXTINCT_TABLE\n".format(xml_info["EXTINCT_TABLE"][0])) 
            sof.write("../Raw/{0}.fits FILTER_LIST\n".format(xml_info["FILTER_LIST"][0])) 
            sof.write("../Raw/{0}.fits LSF_PROFILE\n".format(xml_info["LSF_PROFILE"][0])) 

            #if using reference set it
            if(refcube):
                sof.write("{0} OUTPUT_WCS\n".format(refcube)) 

            for ifu in range(24):
                if(wcsoff):
                    #check if RA/Dec corrected pix tab exists
                    oldpixtab="PIXTABLE_OBJECT_{0:04d}-{1:02d}.fits".format(exp+1,ifu+1)
                    ifupixtab="PIXTABLE_OBJECT_{0:04d}-{1:02d}_off.fits".format(exp+1,ifu+1)
                    
                    if not os.path.isfile(ifupixtab):
                        print 'Correcting RA/Dec in pix table for ifu ', ifu+1
                        #make a copy
                        shutil.copyfile(oldpixtab,ifupixtab)
                        #update header with RA/Dec
                        pixtabfits=fits.open(ifupixtab, mode='update')
                        pixtabfits[0].header['RA']=pixtabfits[0].header['RA']-wcsoff[0][exp]
                        pixtabfits[0].header['DEC']=pixtabfits[0].header['DEC']-wcsoff[1][exp]
                        pixtabfits.flush()
                        pixtabfits.close()
                    else:
                        print 'Using existing corrected pixel tables for ifu', ifu+1
                else:
                    #handle case of no offset
                    ifupixtab="PIXTABLE_OBJECT_{0:04d}-{1:02d}.fits".format(exp+1,ifu+1)
                #now write the pix tab in sof
                sof.write("{} PIXTABLE_OBJECT\n".format(ifupixtab)) 
                
            #finish writing sof
            sof.write("STD_RESPONSE_0001.fits STD_RESPONSE\n")
            sof.write("STD_TELLURIC_0001.fits STD_TELLURIC\n")
            sof.close()

            #Write the command file 
            scr=open("../Script/make_scipost_{0:d}.sh".format(exp+1),"w")
            scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 

            scr.write('esorex --log-file=scipost_{0:d}.log muse_scipost --skymethod="none" --filter=white --save=cube,individual ../Script/scipost_{0:d}.sof'.format(exp+1))
            scr.close()
    
            #Run pipeline 
            subprocess.call(["sh", "../Script/make_scipost_{0:d}.sh".format(exp+1)])    
            subprocess.call(["mv","DATACUBE_FINAL.fits",cname])
            subprocess.call(["mv","IMAGE_FOV_0001.fits",iname])
            subprocess.call(["mv","PIXTABLE_REDUCED_0001.fits",pname])
