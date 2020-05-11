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
from datetime import datetime as tm

def preprocess_cubes(cubeslist,zmin,zmax,xpsf=None,ypsf=None,inpath='./',outpath='./'):

    """

    Utility that prepare the cubes for cubex extraction
      > select range
      > substract psf
      > subtract background

    cubeslist -> array of cubes to process
    zmin -> min slice to select 
    zmax -> max slice to select
    x,ypsf -> x,y centroid of psf to subtract [can be None or array]
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
        if(xpsf is not None):
            xpsf=np.asarray(xpsf)
            ypsf=np.asarray(ypsf)
            
            for ii in range(len(xpsf)): 
                outname="{}/{}{}".format(outpath,rootname,"_trim_psfsub.fits")
                subprocess.call(["CubePSFSub","-cube",inname,"-out",outname,"-withvar",".false.","-x","{}".format(xpsf[ii]),"-y","{}".format(ypsf[ii])])
                subprocess.call(["Cube2Im","-cube",outname])
                inname=outname

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
        print 'HERE'
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
        
        x=int(catalog['x_geow'][j])
        y=int(catalog['y_geow'][j])
        z=int(catalog['z_geow'][j])
        
        okpix = (segmap == Id)

        #compute fraction of pixels in 5x5x5 region
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

def independent_SNR_fast(catalog, covariance, segmap, cube, var, cube_med=None, cube_odd=None, cube_even=None, 
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
        
        x=int(catalog['x_geow'][j])
        y=int(catalog['y_geow'][j])
        z=int(catalog['z_geow'][j])
        x1 = int(catalog['x_min'][j])-1
        x2 = int(catalog['x_max'][j])
        y1 = int(catalog['y_min'][j])-1
        y2 = int(catalog['y_max'][j])
        z1 = int(catalog['z_min'][j])-1
        z2 = int(catalog['z_max'][j])
        
        cutsegmap=segmap[z1:z2,y1:y2,x1:x2]
        okpix = (cutsegmap == Id)
        
        #compute fraction of pixels in 5x5x5 region
        cutcutsegmap=cutsegmap[z-2-z1:z+3-z1,y-2-y1:y+3-y1,x-2-x1:x+3-x1]
        cutokpix = (cutcutsegmap == Id)
        catalog['BoxFraction'][j]=1.*np.sum(cutokpix)/np.sum(okpix)
        
        #compute corrected SN
        cutcube = cube[z1:z2,y1:y2,x1:x2]
        SNR      = np.nansum(cube[z1:z2,y1:y2,x1:x2][okpix])/np.sqrt(np.nansum(var[z1:z2,y1:y2,x1:x2][okpix]))
        catalog['SNR'][j] = SNR/covariance[j]
        catalog['covfac'][j] = covariance[j]
        
        if(cube_med is not None):
            try:
                SNR_med  = np.nansum(cube_med[z1:z2,y1:y2,x1:x2][okpix])/np.sqrt(np.nansum(var_med[z1:z2,y1:y2,x1:x2][okpix]))
                catalog['SNR_med'][j] = SNR_med/covariance[j]
            except:
                print('Have you correctly set median cubes/variance?')
                exit()
        
        if(cube_odd is not None):
            try:
                SNR_odd  = np.nansum(cube_odd[z1:z2,y1:y2,x1:x2][okpix])/np.sqrt(np.nansum(var_odd[z1:z2,y1:y2,x1:x2][okpix]))
                catalog['SNR_odd'][j]  = SNR_odd/covariance[j]
            except:
                print('Have you correctly set odd cubes/variance?')
                exit()

        if(cube_even is not None):
            try:
                SNR_even = np.nansum(cube_even[z1:z2,y1:y2,x1:x2][okpix])/np.sqrt(np.nansum(var_even[z1:z2,y1:y2,x1:x2][okpix]))        
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

def make_images_fast(cubelist, segcube, header, catentry, Id, outdir, outnamelist, padding=0):
    
    """
    Utility function that loops over cubes and generate cutouts 

    cubelist  -> list of cubes to process
    id        -> obj id to process
    outdir    -> output directory 
    outformat -> tags that identify the cubes in the list 
    padding   -> padding to be used to make a postage stamp

    """
    
    pixmask  = (segcube == Id)*1.
    mz,my,mx=segcube.shape
    
    zgeo = int(catentry['z_geow'])-1
    xgeo = int(catentry['x_geow'])-1
    ygeo = int(catentry['y_geow'])-1
    
    x1   = int(catentry['x_min'])-1
    x2   = int(catentry['x_max'])
    y1   = int(catentry['y_min'])-1
    y2   = int(catentry['y_max'])
    z1   = int(catentry['z_min'])-1
    z2   = int(catentry['z_max'])
    
    xpad1 = max(x1-padding,0)
    xpad2 = min(x1+padding,mx)
    ypad1 = max(y1-padding,0)
    ypad2 = min(y1+padding,my)
    zpad1 = max(z1-2,0)
    zpad2 = min(z2+2,mz)
    
    
    imahdr = header.copy()
    del imahdr['C*3']
    
    imapadhdr = imahdr.copy()
    imapadhdr['CRPIX1'] -= xpad1
    imapadhdr['CRPIX2'] -= ypad1
    
    cubpadhdr = header.copy()
    cubpadhdr['CRPIX1'] -= xpad1
    cubpadhdr['CRPIX2'] -= ypad1
    cubpadhdr['CRPIX3'] -= zpad1
    
    #Make segmentation image 
    segima = np.nansum(pixmask, axis=0)
    outima1 = outdir+'Image_id{}_det'.format(Id)+".fits"
    hdu = fits.PrimaryHDU(segima, header=imahdr)
    hdu.writeto(outima1, overwrite=True)
    
    if padding>0: 
       outima2 = outdir+'Pstamp_id{}_det'.format(Id)+".fits"
       hdu = fits.PrimaryHDU(segima[ypad1:ypad2,xpad1:xpad2], header=imapadhdr)
       hdu.writeto(outima2, overwrite=True)
    
    #trim segcube and save it in 3D
    segmapshort=pixmask[zpad1:zpad2,ypad1:ypad2,xpad1:xpad2]
    savename = outdir+"/segcube.fits"
    hdu=fits.PrimaryHDU(segmapshort, header=cubpadhdr)
    hdu.writeto(savename,overwrite=True)

    
    for ii in range(len(cubelist)):
        
        thisima = np.nansum(cubelist[ii]*pixmask, axis=0)
        empty = (thisima==0)
        thisima[empty] = cubelist[ii][zgeo,empty]
        
        outima1 = outdir+'Image_id{}'.format(Id)+outnamelist[ii]+".fits"
        hdu = fits.PrimaryHDU(thisima, header=imahdr)
        hdu.writeto(outima1, overwrite=True)
        
        if padding>0: 
           outima2 = outdir+'Pstamp_id{}'.format(Id)+outnamelist[ii]+".fits"
           hdu = fits.PrimaryHDU(thisima[ypad1:ypad2,xpad1:xpad2], header=imapadhdr)
           hdu.writeto(outima2, overwrite=True)
           
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


def finalcatalogue(fcube,fcube_var,catname,target_z=None,rest_line=None,vel_cut=None,
                       cov_poly=None,working_dir='./',output_dir='./',fcube_median=None,fcube_odd=None,
                       fcube_even=None,fcube_median_var=None,fcube_odd_var=None,
                       fcube_even_var=None,fcube_orig=None,fsource_img=None,marzred=None,SNcut=[7,5],
                       DeltaEOSNcut=[0.5,0.5],SNEOcut=[3,3],fracnpix=None,derived=True,checkimg=True,mask=None,startind=0):
    
    """
    
    Function that processes the cubex catalogue to: 
    1. compute actual SN including covariance term
    2. compute velocity offsets relative to target redshift 
    3. applies cuts on SN, continuum source, and fraction of pixels in segmap
    4. generate inspection images and full spectra
    

    fcube        -> cube used for cubex extraction, assumed to be in working_dir
    fcube_var    -> var cube for cubex extraction
    catname      -> catalogue produced by cubex 
    target_z     -> [optional] if set and rest_line set, computes velocity offset from target z
                   assuming every detection is at rest_line
    rest_line    -> [optional] see target_z
    vel_cut      -> [optional] if set, and target_z and rest_line are set as well, it will cut 
                    the catalogue to include only detections within abs(target_z-source_z) <= vel_cut
    cov_poly     -> [optional] third order polinomial model for noise set as array with (N2,N1,N0)
                   If set, covariance is modelled
    working_dir  -> where cubex output is generated
    output_dir   -> where the output of this function is generated
    fcube_median -> [optional] median cube for extraction of QA images
    fcube_odd    -> [optional] odd cube for extraction of QA image and even/odd SN cut
    fcube_even   -> [optional] even cube for extraction of QA image and even/odd SN cut
    fcube_median_var -> [optional] associated variance
    fcube_odd_var -> [optional] associated variance
    fcube_even_var -> [optional] associated variance
    fcube_orig -> [optional] if set to full cube, used for spectral extraction, absolute path can be used.
    source_img -> [optional] if set to aperture image (and marzred set), exclude emitters projected against 
                  continuum source of known redshift
    marzred  -> [optional] redshift file for continuum sources, ala marz (see above)
    SNcut    -> array of SN cuts on main detection, that defines classes. NB classes are assigned 
                in order of which condition they meet first.
    DeltaEOSNcut -> array of percentage change in even/odd SN allowed (if even/odd cubes provided)
    SNEOcut  -> array of SN cut for even odd exposure
    fracnpix -> [optional] if set, cuts object with less than fracnpix in segmap within 5x5x5 region
    derived  -> if true, compute derived quantities (SN etc..)
    checkimg -> if true generate image cut outs 
    mask --> [optional] If a mask is provided (1=bad, 0=ok) the sources with their center on a masked 
             pixel will be rejected
    
    """

    
    #Load cubex catalog
    catalog = read_cubex_catalog(working_dir+catname)

    
    #create space for bunch of keywords [some may not be used]
    ksig           = Column(np.zeros(len(catalog), dtype=float), name='SNR')
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
    elif cov_poly.ndim == 1 :
        size = np.sqrt(catalog['area_isoproj'])*0.2
        covariance = np.polyval(cov_poly,size)
    elif cov_poly.ndim == 2 :
        size = np.sqrt(catalog['area_isoproj'])*0.2
        covariance = np.zeros(len(catalog), dtype=float)
        for ii in range(len(catalog)):
            try:
             okind = np.where(catalog['lambda_geow'][ii]>cov_poly[:,0])[0][-1]
            except:
             okind = 0 
            covariance[ii] = np.polyval(cov_poly[okind,2:],size[ii])
       
    #Open cubes, take header from data extension of average cube
    #We assume all the other products share the same WCS!!!
    print("Reading cubes")
    cubehdu = fits.open(working_dir+fcube)
    try:
       cube      = cubehdu[1].data
       cubehdr   = cubehdu[1].header
    except:
       cube      = cubehdu[0].data
       cubehdr   = cubehdu[0].header
       
    try:
       cube_var  = fits.open(working_dir+fcube_var)[1].data
    except:
       cube_var  = fits.open(working_dir+fcube_var)[0].data
        
    #reconstruct name of cubex segmentation cube and open
    fsegmap = os.path.basename(fcube).split('.fits')[0]+".Objects_Id.fits"
    segmap  = fits.open(working_dir+fsegmap)[0].data
    
    if(fcube_odd):
        try:
            try:
              cube_odd          = fits.open(working_dir+fcube_odd)[1].data
            except:
              cube_odd          = fits.open(working_dir+fcube_odd)[0].data
            try:   
              cube_odd_var      = fits.open(working_dir+fcube_odd_var)[1].data
            except:
              cube_odd_var      = fits.open(working_dir+fcube_odd_var)[0].data
        except:
            print('Have you set the right odd cube/variance?')
            exit()
    else:
        cube_odd=None
        cube_odd_var=None

    if(fcube_even):
        try:
            try:
              cube_even          = fits.open(working_dir+fcube_even)[1].data
            except:
              cube_even          = fits.open(working_dir+fcube_even)[0].data
            try:   
              cube_even_var      = fits.open(working_dir+fcube_even_var)[1].data
            except:
              cube_even_var      = fits.open(working_dir+fcube_even_var)[0].data
        except:
            print('Have you set the right even cube/variance?')
            exit()
    else:
        cube_even=None
        cube_even_var=None

    if(fcube_median):
        try:
            try:
              cube_median       = fits.open(working_dir+fcube_median)[1].data
            except:
              cube_median       = fits.open(working_dir+fcube_median)[0].data
            try:
              cube_median_var   = fits.open(working_dir+fcube_median_var)[1].data
            except:
              cube_median_var   = fits.open(working_dir+fcube_median_var)[0].data
        except:
            print('Have you set the right median cube/variance?')
            exit()
    else:
        cube_median=None
        cube_median_var=None

    #open source image
    if(fsource_img) and (marzred):
        apermap=fits.open(fsource_img)[0].data
        try:
            contzcat=ascii.read(marzred,format='csv',header_start=2)
        except:
            print('Have you set the marz redshift file???')
            exit()
    else:
        apermap=None
        contzcat=None

    #Computing and adding velocity offset
    if(target_z is not None):
        veloff = velocityoffset(catalog['lambda_fluxw'], target_z, rest_line)
        catalog['veloffset'] = veloff
        
        #Trim in velocity if requested
        if(vel_cut is not None):
           select=np.abs(catalog['veloffset'])<=vel_cut
           catalog=catalog[select]
           
    #now loop over sources and update the catalogue info 
    if(derived):
        print("Calculating independent SNR and additional metrics for {} sources".format(len(catalog)))
        catalog = independent_SNR_fast(catalog, covariance, segmap, cube, cube_var, cube_med=cube_median,cube_odd=cube_odd, cube_even=cube_even, var_med=cube_median_var,var_odd=cube_odd_var, var_even=cube_even_var,apermap=apermap,contzcat=contzcat)
        
    #Compute EOSN 
    if(fcube_even is not None):
        rel_diff_halves = np.abs(catalog['SNR_even']-catalog['SNR_odd'])/np.minimum(catalog['SNR_even'],catalog['SNR_odd'])
        catalog['EODeltaSN'] = rel_diff_halves

    
    #make space for checks
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
        
    #Write full catalogue with SNR and derived quantities
    print("Writing full catalog to disk")
    catalog.write(output_dir+catname.split(".cat")[0]+"_all_SNR.fits", format="fits", overwrite=True)

    
    #make simplest cut to catalogue to reject unwanted sources
    select=catalog['SNR']>=np.amin(SNcut)
    catalog=catalog[select]
 
    #If mask is provided apply it
    if (mask):
       hdu = fits.open(mask)
       try:
         msk = hdu[0].data
       except:
         msk = hdu[1].data
       
       masked = msk[np.array(catalog['y_geow'], dtype=int),np.array(catalog['x_geow'], dtype=int)]  
       catalog=catalog[masked==0]
       
    #loop over classes and assign
    print("Assigning classes")
    for iclass,iSN in enumerate(SNcut):

        #case of SN,DeltaSN,fractpix,continnum
        if((fcube_even is not None) & (fracnpix is not None) & (fsource_img is not None)):
            thisclass=((catalog['SNR'] >= iSN) & (catalog['SNR_odd'] >= SNEOcut[iclass]) & (catalog['SNR_even'] >= SNEOcut[iclass]) & (catalog['EODeltaSN'] <= DeltaEOSNcut[iclass]) & (catalog['OverContinuum'] == False) & (catalog['BoxFraction'] >= fracnpix) & (catalog['confidence'] == 0))
            catalog['confidence'][thisclass] = iclass+1

        #case of SN,DeltaSN,fractpix
        elif ((fcube_even is not None) & (fracnpix is not None)):
            thisclass=((catalog['SNR'] > iSN) & (catalog['SNR_odd'] > SNEOcut[iclass]) & (catalog['SNR_even'] > SNEOcut[iclass]) & (catalog['EODeltaSN']<DeltaEOSNcut[iclass]) & (catalog[BoxFraction] > fracnpix) & (catalog['confidence'] == 0))
            catalog['confidence'][thisclass] = iclass+1
        
        #case of SN,DeltaSN
        elif(fcube_even is not None):
            thisclass=((catalog['SNR'] > iSN) & (catalog['SNR_odd'] > SNEOcut[iclass]) & (catalog['SNR_even'] > SNEOcut[iclass]) & (catalog['EODeltaSN']<DeltaEOSNcut[iclass]) & (catalog['confidence'] == 0))
            catalog['confidence'][thisclass] = iclass+1
        
        #remaining cases
        else:
            thisclass=((catalog['SNR'] > iSN) & (catalog['confidence'] == 0)) 
            catalog['confidence'][thisclass] = iclass+1
        
        
    #Write full catalogue with SNR and derived quantities
    catalog.write(output_dir+catname.split(".cat")[0]+"_select_SNR.fits", format="fits", overwrite=True)

    #make space for checks
    if not os.path.isdir(output_dir+"/objs"):
        os.mkdir(output_dir+"/objs")

    if(checkimg):
        print("Extracting images for {} sources".format(len(catalog)))
        #loop over detections
        
        for ii in range(startind,len(catalog)):
            
            #folder for this object
            objid = catalog['id'][ii]
            
            objdir = output_dir+"/objs/id{}/".format(objid)
            if not os.path.isdir(objdir):
                os.mkdir(objdir)

            #make images
            hdulist=[cube, cube_median, cube_odd, cube_even]
            outnamelist=['_mean', '_median', '_half1', '_half2']
            make_images_fast(hdulist, segmap, cubehdr, catalog[ii], objid, objdir, outnamelist, padding=50)
            
            #Uncomment if you want to use cubex to make the images, (SLOW!)
            #hcubelist=[fsegmap,fcube,fcube_median,fcube_odd,fcube_even]
            #namelist=[working_dir+thc for thc in hcubelist]
            #taglist=['Image_id{}_det','Image_id{}_mean', 'Image_id{}_median', 'Image_id{}_half1', 'Image_id{}_half2']
            #make_cubex_images(namelist, namelist[0], objid, objdir,taglist, padding=-1)
            #taglist=['Pstamp_id{}_det','Pstamp_id{}_mean', 'Pstamp_id{}_median', 'Pstamp_id{}_half1', 'Pstamp_id{}_half2']
            #make_cubex_images(namelist, namelist[0], objid, objdir,taglist, padding=50)

            #Extract spectrum
            if(fcube_orig is not None):
                savename = objdir+"/spectrum.fits".format(objid)
                utl.cube2spec(fcube_orig, 0.0, 0.0, 0.0 , shape='mask', helio=0, mask=segmap, twod=True, tovac=True, write=savename, idsource=objid)


def emi_cogphot(fcube, fcube_var, fsegcube, fcatalog, idlist, dz=24, maxrad=15, growthlim=1.025, plots=False):
    
    try:
      cube = fits.open(fcube)[1].data
    except:
      cube = fits.open(fcube)[0].data
    try:
      cubevar = fits.open(fcube_var)[1].data
    except:
      cubevar = fits.open(fcube_var)[0].data
      
    
    segcube = fits.open(fsegcube)[0].data
    catalog = fits.open(fcatalog)[1].data
    
    fluxarr = np.zeros(len(idlist))
    errarr  = np.zeros(len(idlist))
    radarr  = np.zeros(len(idlist))
    
    for ind, eid in enumerate(idlist):
    
      #Find xc, yc, zc in catalog
      emirec = catalog['id'] == eid
      
      if emirec.sum() <1:
         print('ID {} not found in Cubex catalog. Skipping..'.format(eid))
         continue
      
      xc, yc, zc = catalog['x_fluxw'][emirec][0],  catalog['y_fluxw'][emirec][0],  catalog['z_fluxw'][emirec][0]
      ny, nx = np.shape(cube)[1:]

      tmpnbima = np.zeros((ny, nx))
      tmpnbvar = np.zeros((ny, nx))

      for zz in np.arange(np.floor(zc-dz/2.), np.ceil(zc+dz/2.), dtype=int):
          neigh = (segcube[zz,...] != eid) & (segcube[zz,...] > 0)
          
          tmpslice = np.nan_to_num(cube[zz,...])
          tmpslice[neigh] = 0
          tmpnbima += tmpslice
          
          tmpslice = np.nan_to_num(cubevar[zz,...])
          tmpslice[neigh] = 0
          tmpnbvar += tmpslice
      
      #Go in flux 
      tmpnbima *= 1.25
      tmpnbvar *= 1.25
      
      rad = np.arange(1,maxrad+1)
      phot    = np.zeros_like(rad, dtype=float)
      photerr = np.zeros_like(rad, dtype=float)
      growth  = np.zeros_like(rad, dtype=float)
      
      skyaper = CircularAnnulus((xc-1, yc-1), r_in=maxrad, r_out = 1.5*maxrad)
      skymask = skyaper.to_mask(method='center')
      skydata = skymask.multiply(tmpnbima)[skymask.data>0]

      bkg_avg, bkg_med, _  = scl(skydata)       
      
      for ii, r in enumerate(rad):
        
        aper = CircularAperture((xc-1, yc-1), r=r)
        phot[ii]    = (aperture_photometry(tmpnbima, aper))['aperture_sum'][0]-bkg_med*aper.area
        photerr[ii] = np.sqrt((aperture_photometry(tmpnbvar,aper))['aperture_sum'][0])
        if ii==0:
           growth[ii] = 100
        else:
           growth[ii] = phot[ii]/phot[ii-1]    
          
      rlim = np.argmin(growth>growthlim)-1
      
      fluxarr[ind] = phot[rlim]
      errarr[ind] = photerr[rlim]
      radarr[ind] = rad[rlim]
      
      if plots:
        fig, ax = mp.subplots(nrows=1, ncols=1, figsize=(5,5))
        ax.imshow(convolve(tmpnbima, Box2DKernel(5)), vmin=-2, vmax=30, origin='lower')
        #ax.imshow(tmpnbima, vmin=-2, vmax=30, origin='lower')
        circ = mp.Circle((xc,yc), rad[rlim], color='r', fill=False, lw=3)
        ax.add_artist(circ)
        ax.set_xlim(xc-50,xc+50)
        ax.set_ylim(yc-50,yc+50)
        mp.show() 
      
        mp.plot(rad, phot)
        mp.show()
      
    return fluxarr, errarr, radarr    
