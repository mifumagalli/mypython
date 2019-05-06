"""

General code to handle the ID of emitters
See e.g. Lofthouse et al. 2019, Fossati et al. 2019

Depend on proprietary code [cubex]

"""

import subprocess
import os 

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
