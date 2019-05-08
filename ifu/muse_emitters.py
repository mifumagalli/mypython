"""

General code to handle the ID of emitters
See e.g. Lofthouse et al. 2019, Fossati et al. 2019

Depend on proprietary code [cubex]

"""

import subprocess
import os 
import numpy as np
import shutil
import mypython as mp
from mypython.ifu import muse
from mypython.ifu import muse_utils as utl
from mypython.ifu import muse_source as msc
from mypython.spectra import onedutils as spc
from astropy.io import fits, ascii
from astropy.table import Column, vstack, Table, join
import sys    
import matplotlib.pyplot as plt

def preprocess_cubes(cubeslist,zmin,zmax,xpsf=None,ypsf=None,inpath='./',outpath='./'):

    """

    Utility that prepare the cubes for cubex extraction
      > select range
      > substract psf
      > subtract background

    cubeslist -> array of cubes to process
    zmin -> min slice to select 
    zmax -> max slice to select
    x,ypsf -> x,y centroid of psf to subtract [can be None]
    inpath -> path where original cubes reside
    outpath -> path where cubes are written 
    

    """

    for cube in cubeslist:
        
        rootname=cube.split('.fits')[0]
        print('Process ',rootname)
        
        #select cube
        outname="{}/{}{}".format(outpath,rootname,"_trim.fits")
        inname="{}/{}".format(inpath,cube)
        subprocess.call(["CubeSel","-cube",inname,"-out",outname,"-multiext",".false.","-OutVarCube",".true.","-zmin","{}".format(zmin),"-zmax","{}".format(zmax)])

        #psf if desired
        inname=outname
        if(xpsf):
            outname="{}/{}{}".format(outpath,rootname,"_trim_psfsub.fits")
            subprocess.call(["CubePSFSub","-cube",inname,"-out",outname,"-withvar",".false.","-x","{}".format(xpsf),"-y","{}".format(ypsf)])
            subprocess.call(["Cube2Im","-cube",outname])

        #background
        inname=outname
        outname="{}/{}{}".format(outpath,rootname,"_trim_psfsub_bkg.fits")
        subprocess.call(["CubeBKGSub","-cube",inname,"-out",outname,"-bpsize","1 1 40","-bfrad","0 0 3"])
        subprocess.call(["Cube2Im","-cube",outname])
        
        
def run_cubex(cubeslist,varlist,zrange,cubexpar,mask=None,outdir='./'):

    """

    Utility that prepare the cubes for cubex extraction
      > select range
      > substract psf
      > subtract background

    cubeslist -> array of cubes to process [trim, and extract from [0]]
    varlist   -> associated variance cube list
    zrange -> min,max slice to select 
    cubexpar -> parameter file for cubex
    outpath -> path where cubes are written 
    
    """
    
    #make space for output if needed
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
 

    #trim cubes to desired lenght
    for cub in cubeslist:
        subprocess.call(["CubeSel", "-InpFile", cub, "-out", outdir+cub.replace('.fits', '_sel.fits'), "-zmin", "{}".format(zrange[0]), "-zmax", "{}".format(zrange[1])])

    #trim associated variance to desired lenght 
    for var in varlist:
        subprocess.call(["CubeSel", "-InpFile", var, "-out", outdir+var.replace('.VAR.fits', '_sel.VAR.fits'), "-zmin", "{}".format(zrange[0]), "-zmax", "{}".format(zrange[1])])

        
    #mv to new location
    olddir=os.getcwd()
    os.chdir(outdir)
    
    #no run with option of masking 
    if(mask):
        subprocess.call(["CubEx", cubexpar, "-InpFile", cubeslist[0].replace('.fits', '_sel.fits'), "-Catalogue", "MUDF_cubex_{}_{}.cat".format(zrange[0], zrange[1]), "-SourceMask", mask, "-var", varlist[0].replace('.VAR.fits', '_sel.VAR.fits')])
    else:
        subprocess.call(["CubEx", cubexpar, "-InpFile", cubeslist[0].replace('.fits', '_sel.fits'), "-Catalogue", "MUDF_cubex_{}_{}.cat".format(zrange[0], zrange[1]),"-var", varlist[0].replace('.VAR.fits', '_sel.VAR.fits')])
     

    #back to where we were
    os.chdir(olddir)
    

def read_cubex_catalog(catname):

    """

    Utility function that reads the emitter catalogue

    catname -> name of cat to read

    """
    
    fields=['id', 'n_spax_thresh', 'x_geow', 'y_geow', 'z_geow', 'x_fluxw', 'y_fluxw', 'z_fluxw', 'x_min', 'y_min', 'z_min', 'x_max', 'y_max', 'z_max', 'area_isoproj', 'area_z', 'IsoFlux', 'IsoErr', 'flux_minmax', 'err_minmax', 'flux_aper', 'err_aper', 'ra_geow', 'Dec_geow', 'lambda_geow', 'ra_fluxw', 'Dec_fluxw', 'lambda_fluxw', 'assoc_Id']

    try:
        catalog = ascii.read(catname, names=fields)
        
    except:
        #Open file manually and perge bad rows
        f = open(catname, 'r')
        o = open(catname, 'w')
        for line in f:                              
            if (line[0] == '#')or(len(line.split()) == 29): o.write(line)
        o.close()
        f.close()
        catalog = ascii.read(catname, names=fields)

    return catalog
            

def independent_SNR(catalog, covariance, segmap, cube, var, cube_med=None, cube_odd=None, cube_even=None, 
                    var_med=None, var_odd=None, var_even=None,apermap=None,contzcat=None):

    """
    
    Routine that computes effective SN of detections and other metrics in various cubes, 
    also accounting for correction due to correlated noise
  
    catalog       -> source catalog
    covariance    -> correction term to the noise for correlation, once per source
    segmap        -> object segmentation cube generated by extraction 
    cube          -> main analysis cube
    cube_med      -> median cube [optional]
    cube_odd      -> odd cube containing 1/2 exposures [optional]
    cube_even     -> even cube containing 1/2 exposures [optional]
    var           -> main analysis cube variance 
    var_med       -> median cube variance [optional]
    var_odd       -> odd cube variance containing 1/2 exposures [optional]
    var_even      -> even cube variance containing 1/2 exposures [optional]
    apermap       -> 2D aperture map of continuum detected sources [optional] 
    contzcat      -> redshift catalogue written with marz format

    """

    #loop over sources
    for j in range(len(catalog)):

        #find pixels in seg map
	Id = catalog['id'][j]
        okpix = (segmap == Id)

        #compute fraction of pixels in 5x5x5 region
        x=int(catalog['x_geow'][j])
        y=int(catalog['y_geow'][j])
        z=int(catalog['z_geow'][j])
        cutsegmap=segmap[z-2:z+3,y-2:y+3,x-2:x+3]
        cutokpix = (cutsegmap == Id)
        catalog['BoxFraction'][j]=1.*np.sum(cutokpix)/np.sum(okpix)
        
        
        #compute corrected SN
        SNR      = np.nansum(cube[okpix])/np.sqrt(np.nansum(var[okpix]))
        catalog['SNR'][j] = SNR/covariance[j]
        catalog['covfac'][j] = covariance[j]

        if(cube_med is not None):
            try:
                SNR_med  = np.nansum(cube_med[okpix])/np.sqrt(np.nansum(var_med[okpix]))
                catalog['SNR_med'][j] = SNR_med/covariance[j]
            except:
                print('Have you correctly set median cubes/variance?')
                exit()

        if(cube_odd is not None):
            try:
                SNR_odd  = np.nansum(cube_odd[okpix])/np.sqrt(np.nansum(var_odd[okpix]))
                catalog['SNR_odd'][j]  = SNR_odd/covariance[j]
            except:
                print('Have you correctly set odd cubes/variance?')
                exit()

        if(cube_even is not None):
            try:
                SNR_even = np.nansum(cube_even[okpix])/np.sqrt(np.nansum(var_even[okpix]))	
                catalog['SNR_even'][j] = SNR_even/covariance[j]
            except:
                print('Have you correctly set even cubes/variance?')
                exit()
                
        #compute if overlap with continuum source with known redshift 
        if(apermap is not None):
            #find if source exist at position 
            cntid=apermap[y,x]
            if(cntid > 0):
                #now id redshift of this source
                zid=(contzcat['#ID']==int(cntid))
                if(contzcat['QOP'][zid]>2):
                    catalog['OverContinuum'][j] = True
                
    #return updated catalogue            
    return catalog

    
def make_cubex_images(cubelist, fsegmap, id, outdir, outformat, padding=-1):
    
    """
    Utility function that loops over cubes and generate cutouts 

    cubelist  -> list of cubes to process
    id        -> obj id to process
    outdir    -> output directory 
    outformat -> tags that identify the cubes in the list 
    padding   -> padding to be used by cubex cube2im

    """

    for ii in range(len(cubelist)):
        outima = outdir+outformat[ii].format(id)+".fits"
        subprocess.call(["Cube2Im", "-cube", cubelist[ii], "-idcube", fsegmap, "-out", outima,"-id", '{}'.format(id), "-idpad", '{}'.format(padding), "-nl", "-1"])
        

def velocityoffset(wave_air, ztarg, rest_line):
    
    """

    Utility that compute velocity offsets from a reference redshift 

    wave_air  -> air wavelength in A 
    ztarg  -> desired redshift
    rest_line -> vacuum rest frame line in A

    """

    #go to vacuum frame for redshift measurement
    wave_vac=spc.airtovac(wave_air)
    z_vac = (wave_vac/rest_line)-1
        
    vel_offset = (z_vac-ztarg)/(1.+ztarg)*299792.458
    
    return vel_offset


def finalcatalogue(fcube,fcube_var,catname,target_z=None,rest_line=None,
                       cov_poly=None,working_dir='./',fcube_median=None,fcube_odd=None,
                       fcube_even=None,fcube_median_var=None,fcube_odd_var=None,
                       fcube_even_var=None,fcube_orig=None,fsource_img=None,marzred=None,SNcut=[7,5],
                       DeltaEOSNcut=[0.5,0.5],SNEOcut=[3,3],fracnpix=None,derived=True,checkimg=True):
    
    """
    
    Function that processes the cubex catalogue to: 
    1. compute actual SN including covariance term
    2. compute velocity offsets relative to target redshift 
    3. applies cuts on SN, continuum source, and fraction of pixels in segmap
    4. generate inspection images and full spectra
    

    fcube        -> cube used for cubex extraction
    fcube_var    -> var cube for cubex extraction
    catname      -> catalogue produced by cubex 
    target_z     -> [optional] if set and rest_line set, computes velocity offset from target z
                   assuming every detection is at rest_line
    rest_line    -> [optional] see target_z
    cov_poly     -> [optional] third order polinomial model for noise set as array with (N2,N1,N0)
                   If set, covariance is modelled
    working_dir  -> where cubex output is generated
    fcube_median -> [optional] median cube for extraction of QA images
    fcube_odd    -> [optional] odd cube for extraction of QA image and even/odd SN cut
    fcube_even   -> [optional] even cube for extraction of QA image and even/odd SN cut
    fcube_median_var -> [optional] associated variance
    fcube_odd_var -> [optional] associated variance
    fcube_even_var -> [optional] associated variance
    fcube_orig -> [optional] if set to full cube, used for spectral extraction
    source_img -> [optional] if set to aperture image (and marzred set), exclude emitters projected against 
                  continuum source of known redshift
    marzred  -> [optional] redshift file for continuum sources, ala marz (see above)
    SNcut    -> array of SN cuts on main detection, that defines classes
    DeltaEOSNcut -> array of percentage change in even/odd SN allowed (if even/odd cubes provided)
    SNEOcut  -> array of SN cut for even odd exposure
    fracnpix -> [optional] if set, cuts object with less than fracnpix in segmap within 5x5x5 region
    derived  -> if true, compute derived quantities (SN etc..)
    checkimg -> if true generate image cut outs 

    """

    
    #Load cubex catalog
    catalog = read_cubex_catalog(working_dir+catname)

    #create space for bunch of keywords [some may not be used]
    ksig	   = Column(np.zeros(len(catalog), dtype=float), name='SNR')
    ksig_odd    = Column(np.zeros(len(catalog), dtype=float), name='SNR_odd')
    ksig_even   = Column(np.zeros(len(catalog), dtype=float), name='SNR_even')
    ksig_med    = Column(np.zeros(len(catalog), dtype=float), name='SNR_med')
    kcov_fac    = Column(np.zeros(len(catalog), dtype=float), name='covfac')
    kconfidence = Column(np.zeros(len(catalog), dtype=float), name='confidence')
    kveloffset  = Column(np.zeros(len(catalog), dtype=float), name='veloffset')
    kdeltasn    = Column(np.zeros(len(catalog), dtype=float), name='EODeltaSN')
    kfraction   = Column(np.zeros(len(catalog), dtype=float), name='BoxFraction')
    kcontinuum  = Column(np.zeros(len(catalog), dtype=bool), name='OverContinuum')
    
    catalog.add_columns([ksig, ksig_odd, ksig_even, ksig_med, kcov_fac, kconfidence, kveloffset,
                         kdeltasn,kfraction,kcontinuum])
    
    #catalog=catalog[0:100]
    
    #Calculate covariance
    if cov_poly is None:
        covariance = np.zeros(len(catalog), dtype=float)+1.0
    else:
        size = np.sqrt(catalog['area_isoproj'])*0.2
        covariance = np.polyval(cov_poly,size)

    #Open cubes
    cube      = fits.open(working_dir+fcube)[0].data
    cube_var  = fits.open(working_dir+fcube_var)[0].data

    #reconstruct name of cubex segmentation cube and open
    fsegmap = os.path.basename(fcube).split('.fits')[0]+".Objects_Id.fits"
    segmap  = fits.open(working_dir+fsegmap)[0].data
    
    if(fcube_odd):
        try:
            cube_odd          = fits.open(working_dir+fcube_odd)[0].data
            cube_odd_var      = fits.open(working_dir+fcube_odd_var)[0].data
        except:
            print('Have you set the right odd cube/variance?')
            exit()
    else:
        cube_odd=None
        cube_odd_var=None

    if(fcube_even):
        try:
            cube_even         = fits.open(working_dir+fcube_even)[0].data
            cube_even_var     = fits.open(working_dir+fcube_even_var)[0].data
        except:
            print('Have you set the right even cube/variance?')
            exit()
    else:
        cube_even=None
        cube_even_var=None

    if(fcube_median):
        try:
            cube_median       = fits.open(working_dir+fcube_median)[0].data
            cube_median_var   = fits.open(working_dir+fcube_median_var)[0].data
        except:
            print('Have you set the right median cube/variance?')
            exit()
    else:
        cube_median=None
        cube_median_var=None

    #open source image
    if(fsource_img):
        apermap=fits.open(fsource_img)[0].data
        try:
            contzcat=ascii.read(marzred,format='csv',header_start=2)
        except:
            print('Have you set the marz redshift file???')
            exit()
    else:
        apermap=None
        contzcat=None


    #now loop over sources and update the catalogue info 
    if(derived):
        print("Calculating independent SNR and additional metrics")
        catalog = independent_SNR(catalog, covariance, segmap, cube, cube_var, cube_med=cube_median,cube_odd=cube_odd, cube_even=cube_even, var_med=cube_median_var,var_odd=cube_odd_var, var_even=cube_even_var,apermap=apermap,contzcat=contzcat)
        
    #Computing and adding velocity offset
    if(target_z is not None):
        veloff = velocityoffset(catalog['lambda_fluxw'], target_z, rest_line)
        catalog['veloffset'] = veloff

    #Compute EOSN 
    if(fcube_even is not None):
        rel_diff_halves = np.abs(catalog['SNR_even']-catalog['SNR_odd'])/np.minimum(catalog['SNR_even'],catalog['SNR_odd'])
        catalog['EODeltaSN'] = rel_diff_halves

        
    #Write full catalogue with SNR and derived quantities
    print("Writing full catalog to disk")
    catalog.write(working_dir+catname.split(".cat")[0]+"_all_SNR.fits", format="fits", overwrite=True)

    
    #make simplest cut to catalogue to reject unwanted sources
    select=catalog['SNR']>=np.amin(SNcut)
    catalog=catalog[select]
    
    #loop over classes and assign
    print("Assigning classes")
    for iclass,iSN in enumerate(SNcut):

        #case of SN,DeltaSN,fractpix,continnum
        if((fcube_even is not None) & (fracnpix is not None) & (fsource_img is not None)):
            thisclass=((catalog['SNR'] >= iSN) & (catalog['SNR_odd'] >= SNEOcut[iclass]) & (catalog['SNR_even'] >= SNEOcut[iclass]) & (catalog['EODeltaSN'] <= DeltaEOSNcut[iclass]) & (catalog['OverContinuum'] == False) & (catalog['BoxFraction'] >= fracnpix))
            catalog['confidence'][thisclass] = iclass+1

        #case of SN,DeltaSN,fractpix
        elif ((fcube_even is not None) & (fracnpix is not None)):
            thisclass=((catalog['SNR'] > iSN) & (catalog['SNR_odd'] > SNEOcut[iclass]) & (catalog['SNR_even'] > SNEOcut[iclass]) & (catalog['EODeltaSN']<DeltaEOSNcut[iclass]) & (catalog[BoxFraction] > fracnpix))
            catalog['confidence'][thisclass] = iclass+1
        
        #case of SN,DeltaSN
        elif(fcube_even is not None):
            thisclass=((catalog['SNR'] > iSN) & (catalog['SNR_odd'] > SNEOcut[iclass]) & (catalog['SNR_even'] > SNEOcut[iclass]) & (catalog['EODeltaSN']<DeltaEOSNcut[iclass]))
            catalog['confidence'][thisclass] = iclass+1
        
        #remaining cases
        else:
            thisclass=(catalog['SNR'] > iSN) 
            catalog['confidence'][thisclass] = iclass+1
        
        
    #Write full catalogue with SNR and derived quantities
    catalog.write(working_dir+catname.split(".cat")[0]+"_select_SNR.fits", format="fits", overwrite=True)

    #make space for checks
    if not os.path.isdir(working_dir+"/objs"):
        os.mkdir(working_dir+"/objs")

    if(checkimg):
        print("Extracting images of sources")
        #loop over detections
        
        for ii in range(len(catalog)):
            
            #folder for this object
            objid = catalog['id'][ii]
            
            objdir = working_dir+"/objs/id{}/".format(objid)
            if not os.path.isdir(objdir):
                os.mkdir(objdir)

            #make images
            hcubelist=[fsegmap,fcube,fcube_median,fcube_odd,fcube_even]
            namelist=[working_dir+thc for thc in hcubelist]
            taglist=['Image_id{}_det','Image_id{}_mean', 'Image_id{}_median', 'Image_id{}_half1', 'Image_id{}_half2']
            make_cubex_images(namelist, namelist[0], objid, objdir,taglist, padding=-1)
            taglist=['Pstamp_id{}_det','Pstamp_id{}_mean', 'Pstamp_id{}_median', 'PStamp_id{}_half1', 'Pstamp_id{}_half2']
            make_cubex_images(namelist, namelist[0], objid, objdir,taglist, padding=50)

            #Extract spectrum
            if(fcube_orig is not None):
                savename = objdir+"/spectrum.fits".format(objid)
                utl.cube2spec(fcube_orig, 0.0, 0.0, 0.0 , shape='mask', helio=0, mask=segmap, twod=True, tovac=True, write=savename, idsource=objid)





