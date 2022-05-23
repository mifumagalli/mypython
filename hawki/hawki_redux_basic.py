#global import
import numpy as np
import glob 
import subprocess
import time
import os 
import copy
from distutils.spawn import find_executable
from astropy.io import fits
import socket
from astropy import stats
from astropy.modeling import models, fitting
import sep
import matplotlib.pyplot as mp

def grabtime(namelist):

    namelist=list(namelist)
    timelist=[]
    for name in namelist:
        #grab time sequence
        tmstr=name.split("HAWKI.")[1]
        tm=time.strptime(tmstr.split(".")[0],'%Y-%m-%dT%H:%M:%S')
        timelist.append(time.mktime(tm))

    #return time in seconds
    return timelist
    
def grabditndit(namelist, path='./Raw/'):
    
    namelist=list(namelist)
    ditlist=[]
    nditlist=[]
    for name in namelist:
      hdu = fits.open(path+name+'.fits')
      head = hdu[0].header
      ditlist.append(head['ESO DET DIT'])
      nditlist.append(head['ESO DET NDIT'])
      hdu.close()
    return np.array(ditlist), np.array(nditlist)   

def grabband(namelist, path='./Raw/'):
    
    namelist=list(namelist)
    bandlist=[]
    for name in namelist:
      hdu = fits.open(path+name+'.fits')
      head = hdu[0].header
      bandlist.append(head['ESO INS FILT1 ID'])
      hdu.close()
    return np.array(bandlist)

def grabob(namelist, path='./Raw/'):
    
    namelist=list(namelist)
    oblist=[]
    for name in namelist:
      hdu = fits.open(path+name+'.fits')
      head = hdu[0].header
      oblist.append(head['ESO OBS ID'])
      hdu.close()
    return np.array(oblist)

def grabtplstart(namelist, path='./Raw/'):
    
    namelist=list(namelist)
    tplstartlist=[]
    for name in namelist:
      hdu = fits.open(path+name+'.fits')
      head = hdu[0].header
      tplstartlist.append(head['ESO TPL START'])
      hdu.close()
    return np.array(tplstartlist)
        
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
        print('Object observations are taken {0:4.3f} h apart'.format(delta_time))
        if(delta_time > 4):
            print('Large time lag bewteen observations! Check!')
    
    #Check DIT NDIT for objects
    dit, ndit = grabditndit(xml_info['OBJECT'])
    band = grabband(xml_info['OBJECT'])
    obname = grabob(xml_info['OBJECT'])
    obname = grabtplstart(xml_info['OBJECT'])    
    if np.any(np.diff(dit)>0) or np.any(np.diff(ndit)>0):
       print('Not all science exposures have been obtained with the same DIT NDIT combination.')
       exit()
    else:
       obj_dit = dit[0]
       obj_ndit = ndit[0]
    if not np.all(band==band[0]):
       print('Not all science exposures have been obtained with the same filter.')
       exit()
    else:
       obj_band = band[0]
       
    #Split objects by obname
    uniqobs = np.unique(obname)
    nobs = len(uniqobs)   
    
    xml_info['NOBS'] = nobs
    objlist = np.array(list(xml_info['OBJECT']))
    for ob in range(nobs):
        xml_info['OBJECT{}'.format(ob+1)] = objlist[obname==uniqobs[ob]]
    
    #Check DIT NDIT for flats
    dit, ndit = grabditndit(xml_info['FLAT_TWILIGHT'])
    band = grabband(xml_info['FLAT_TWILIGHT'])
    
    thisband = band==obj_band
    #b /cosma/home/durham/dc-foss1/local/mypython/hawki/hawki_redux_basic.py:122
    
    fltlist = np.array(list(xml_info['FLAT_TWILIGHT']))
    fltlist = fltlist[thisband]
    dit = dit[thisband]
    ndit = ndit[thisband]
    
    if len(fltlist)==0:
       print('Flats using the same filter as science not found')
       exit()
    else:
       if np.any(np.diff(dit)>0) or np.any(np.diff(ndit)>0):
         print('Not all flat exposures have been obtained with the same DIT NDIT combination.')
         print('Dealing with this situation will be made available in a future release')
         exit()
       else:
         flt_dit = dit[0]
         flt_ndit = ndit[0]
    
    flttime = np.array(grabtime(fltlist))
    delta_time=np.abs(flttime-reference_time)/3600.
    offtime = 1000
    
    recent= (delta_time <= offtime) & (delta_time >0.5)
    
    xml_info['TWILIGHT_FLAT'] = fltlist[recent]
    
    
    #Select DARK for obj and flat
    offtime = 24
    
    dit, ndit = grabditndit(xml_info['DARK'])
    darktime = np.array(grabtime(list(xml_info['DARK'])))
    darklist = np.array(list(xml_info['DARK']))
    
    obj_ditndit = (dit==obj_dit) & (ndit==obj_ndit)
    
    darklist_obj = darklist[obj_ditndit]
    darktime_obj = darktime[obj_ditndit]
    delta_time=np.abs(darktime_obj-reference_time)/3600.
    
    recent= (delta_time <= offtime) & (delta_time >0.5)
    
    xml_info['DARK_OBJ'] = darklist_obj[recent]
    print('Found {0} {1} taken within {2} hours'.format(len(darklist_obj[recent]),'DARKS FOR OBJS', offtime))
    
    offtime = 1e3
    flt_ditndit = (dit==flt_dit) & (ndit==flt_ndit)
    
    darklist_flt = darklist[flt_ditndit]
    darktime_flt = darktime[flt_ditndit]
    delta_time=np.abs(darktime_flt-darktime_flt[0])/3600.
    
    recent= (delta_time <= offtime) 
    
    xml_info['DARK_FLT'] = darklist_flt[recent]
    print('Found {0} {1} taken within {2} hours'.format(len(darklist_flt[recent]),'DARKS FOR FLATS', offtime))
         
    

    #now loop over calibrations and grab the best 
    allkey=xml_info.keys()
    for kk in allkey:
        if('OBJECT' not in kk):
            if(kk == 'DARK_FLT') or (kk == 'DARK_OBJ') or (kk=='NOBS') or (kk=='TWILIGHT_FLAT') or (kk=='FLAT_TWILIGHT'):
                continue               
            #find list
            currentlist=list(xml_info[kk])                                
            #find time offset in seconds
            times=np.array(grabtime(currentlist))
            delta_time=np.abs(times-reference_time)/3600.
            currentlist=np.array(currentlist)
            #now handle by keyword according to calibration plan
            #This is when you want only the best one
            try:
                #pick closest
                mintm=np.argmin(delta_time)
                xml_info[kk]=[currentlist[mintm]]
                print('Best {0} taken within {1} days'.format(kk,delta_time[mintm]/24.))
            except:
                print('No suitable calibration found for {}. Abort.'.format(kk))
                exit()
                
    #set the calibration path relative and suffix
    xml_info["PATHCAL"]='../../Raw/'
    xml_info["SUFFIXCAL"]='.fits'


    staticalpath='dummy'


    #Here sort out things with static calibrations: GEOMETRY & ASTROMETRY 
    #This is the largest time at one should worry about legacy products 
    #print('Using pipeline static calibrations for astrometry and geometry')
    #astrostatic=staticalpath+'astrometry_wcs_wfm.fits'
    #geometrystatic=staticalpath+'geometry_table_wfm.fits'
    
    #update from default
    #xml_info['ASTROMETRY_WCS']=[astrostatic]
    #xml_info['GEOMETRY_TABLE']=[geometrystatic]
                
    #now dump cal plan
    print('Writing calibration plan in calplan.txt')
    cl=open('calplan.txt','w')
    for kk in xml_info.keys():
        if (kk=='NOBS'):
            continue
        for ll in xml_info[kk]:
            if(('PATHCAL' not in kk) & ('SUFFIXCAL' not in kk)):
                cl.write('{0} {1}\n'.format(kk,ll))
    cl.close()

    return xml_info  


def make_dark_obj(xml_info,nproc=12):
    
    #grab the dark
    dark_list=xml_info["DARK_OBJ"]
    
    #Write the sof file 
    sof=open("../Script/dark_obj.sof","w")
    for ii in dark_list:
        sof.write("../Raw/{0}.fits DARK\n".format(ii)) 
    sof.close()
    
    #Write the command file 
    scr=open("../Script/make_dark_obj.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=dark_obj.log hawki_dark_combine ../Script/dark_obj.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../Script/make_dark_obj.sh"])
    subprocess.call(["mv","darkcomb.fits","MASTER_DARK_OBJ.fits"])

def make_dark_flt(xml_info,nproc=12):
    
    #grab the dark
    dark_list=xml_info["DARK_FLT"]
    
    #Write the sof file 
    sof=open("../Script/dark_flt.sof","w")
    for ii in dark_list:
        sof.write("../Raw/{0}.fits DARK\n".format(ii)) 
    sof.close()
    
    #Write the command file 
    scr=open("../Script/make_dark_flt.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=dark_flt.log hawki_dark_combine ../Script/dark_flt.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../Script/make_dark_flt.sh"])
    subprocess.call(["mv","darkcomb.fits","MASTER_DARK_FLT.fits"])
    


def make_flat(xml_info,nproc=12):

    #grab the flat
    flat_list=xml_info["TWILIGHT_FLAT"]
    
    #Write the sof file 
    sof=open("../Script/flat.sof","w")
    for ii in flat_list:
        sof.write("../Raw/{0}.fits FLAT_TWILIGHT\n".format(ii)) 
    sof.write("MASTER_DARK_FLT.fits MASTER_DARK\n")        
    sof.close()
              
    #Write the command file 
    scr=open("../Script/make_flat.sh","w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=flat.log hawki_twilight_flat_combine ../Script/flat.sof")
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../Script/make_flat.sh"])
    subprocess.call(["mv","twilightcomb.fits","MASTER_TWILIGHT_FLAT.fits"])
    subprocess.call(["mv","twilightconf.fits","MASTER_CONF.fits"])
    subprocess.call(["mv","bpmmap.fits","MASTER_BPM.fits"])

def make_objects(xml_info, ob, skyalgo='pawsky_mask', astrocat='wise', nproc=12):

    #grab the objects
    obj_list=list(xml_info["OBJECT{}".format(ob)])
    #sort by time - useful is offsets fields available
    obj_list.sort()
    
    #Write the sof file 
    sof=open("../../Script/object{}.sof".format(ob),"w")
       
    for ii in obj_list:
        sof.write("../../Raw/{0}.fits OBJECT\n".format(ii))
    sof.write("../../Raw/{0}.fits MASTER_READGAIN\n".format(xml_info["MASTER_READGAIN"][0])) 
    sof.write("../../Raw/{0}.fits SCHLEGEL_MAP_NORTH\n".format(xml_info["SCHLEGEL_MAP_NORTH"][0])) 
    sof.write("../../Raw/{0}.fits SCHLEGEL_MAP_SOUTH\n".format(xml_info["SCHLEGEL_MAP_SOUTH"][0])) 
    #sof.write("../../Raw/{0}.fits MASTER_2MASS_CATALOGUE_PHOTOM\n".format(xml_info["MASTER_2MASS_CATALOGUE_PHOTOM"][0])) 
    sof.write("../../Raw/{0}.fits PHOTCAL_TAB\n".format(xml_info["PHOTCAL_TAB"][0])) 
    sof.write("../../Calib/MASTER_DARK_OBJ.fits MASTER_DARK\n") 
    sof.write("../../Calib/MASTER_CONF.fits MASTER_CONF\n") 
    sof.write("../../Calib/MASTER_TWILIGHT_FLAT.fits MASTER_TWILIGHT_FLAT\n") 
    if skyalgo=='pawsky_mask_pre' or skyalgo=='pawsky_minus':
        sof.write("../../../skymask/skymask.fits MASTER_OBJMASK\n") 
    if '.fits' in astrocat:
        sof.write("../../../astrom/{} MASTER_LOCAL_CATALOGUE_ASTROM\n".format(astrocat))
    sof.close()
    
    if '.fits' in astrocat:
       astrocds = 'none'
    else:
       astrocds = astrocat   
    
    #Write the command file 
    scr=open("../../Script/make_stack{}.sh".format(ob),"w")
    scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
    scr.write("esorex --log-file=object.log hawki_science_process --skyalgo=none --cdssearch_astrom={1} --cdssearch_photom=2mass --savemstd=true --savecat=true --magerrcut=0.12 ../../Script/object{2}.sof".format(skyalgo,astrocds,ob))
    scr.close()
    
    #Run pipeline 
    subprocess.call(["sh", "../../Script/make_stack{}.sh".format(ob)])
    
    
def make_fluxcal(xml_info):    
    
    #Now shift ZP for all the objects to a common value of 28
    objlist = glob.glob('exp_?.fits')
    objlist2 = glob.glob('exp_??.fits')
    
    objlist = [*objlist, *objlist2]
    
    for obj in objlist:
        
        hdu = fits.open(obj)
        hduvar = fits.open(obj.replace('exp_', 'exp_var_'))
        
        fluxfac = 10**(-0.4*hdu[1].header['PHOTZP']+11.2)
        
        for ext in range(1,5):
            hdu[ext].header['PHOTZPOR'] = copy.copy(hdu[1].header['PHOTZP'])
            hdu[ext].header['PHOTZP'] = 28
            hdu[ext].header['FLUXFAC'] = fluxfac
            hdu[ext].data *= fluxfac
            hduvar[ext].data *= fluxfac
            
            badpix = hduvar[ext].data==0
            
            hduvar[ext].data = 1./hduvar[ext].data
            hduvar[ext].data[badpix] = 0

        hdu.writeto(obj, overwrite=True)
        hduvar.writeto(obj.replace('exp_', 'exp_var_'), overwrite=True)
        
        os.replace(obj.replace('exp_', 'exp_var_'), obj.replace('.fits', '_var.fits')  )

def make_skysub(xml_info):    
    
    #Subtract sky
    
    objlist = glob.glob('exp_?.fits')
    objlist2 = glob.glob('exp_??.fits')
    
    objlist = [*objlist, *objlist2]
    
    fitter = fitting.LinearLSQFitter()
    
    for obj in objlist:
        
        outfile = obj.replace('.fits', '_bk.fits')
        
        if os.path.isfile(outfile):
           print('Skipping {}. File exists.'.format(obj))
           continue
        else:
           print('Subtracting from {}.'.format(obj))
        
        hdu = fits.open(obj)
        
        for ext in range(1,5):
                        
            data = hdu[ext].data
            swdata = data.byteswap().newbyteorder()
            bkg1 = sep.Background(swdata, bw=128, bh=128)
            obj, segmap = sep.extract(swdata-bkg1, 1.5*bkg1.globalrms, minarea=14, segmentation_map=True, deblend_cont=1.0)
            bkg2 = sep.Background(swdata,bw=128, bh=128, mask=segmap)
            obj, segmap = sep.extract(swdata-bkg2, 1.5*bkg2.globalrms, minarea=14, segmentation_map=True, deblend_cont=1.0)
            bkg3 = sep.Background(swdata,bw=128, bh=128, mask=segmap)
            
            hdu[ext].data -= bkg3
            
            hdu[ext].data[:24,:] = 0
            hdu[ext].data[-24:,:] = 0
            
            hdu[ext].data[:,:24] = 0
            hdu[ext].data[:,-24:] = 0
            
            '''
            filtered  = stats.sigma_clip(data, sigma=2.5, maxiters=10)
            yy,xx = np.mgrid[:hdu[ext].header['NAXIS2'],:hdu[ext].header['NAXIS1']]
            
            p_init = models.Polynomial2D(degree=7)
            p_fit = fitter(p_init, xx, yy, filtered)
            
            hdu[ext].data -= p_fit(xx,yy)
            
            cclip = 5 
            
            
            from astropy.convolution import convolve_fft, Box2DKernel
            smdata = convolve_fft(filtered-np.average(p_fit(xx,yy)), Box2DKernel(20))
            smout = convolve_fft(data, Box2DKernel(20))
            
            fig, ax = mp.subplots(nrows=1, ncols=3, figsize=(15,5))
            ax[0].imshow(smdata, origin='lower', vmin=-1*cclip, vmax=1*cclip)
            ax[1].imshow(p_fit(xx,yy)-np.average(p_fit(xx,yy)), origin='lower', vmin=-1*cclip, vmax=1*cclip)
            ax[2].imshow(smout, origin='lower', vmin=-1*cclip, vmax=1*cclip)
            mp.show()
            '''
            '''
            fig, ax = mp.subplots(nrows=1, ncols=3, figsize=(15,5))
            ax[0].imshow(bkg1, origin='lower')
            ax[1].imshow(bkg2, origin='lower')
            ax[2].imshow(bkg3, origin='lower')
            mp.show()
            
            fig, ax = mp.subplots(nrows=1, ncols=3, figsize=(15,5))
            ax[0].imshow(data, origin='lower')
            ax[1].imshow(data-bkg1, origin='lower', vmin=-10, vmax=10)
            ax[2].imshow(data-bkg3, origin='lower', vmin=-10, vmax=10)
            mp.show()
            '''
                     
            #b /cosma/home/durham/dc-foss1/local/mypython/hawki/hawki_redux_basic.py:412
            

        hdu.writeto(outfile, overwrite=True)
        
        
