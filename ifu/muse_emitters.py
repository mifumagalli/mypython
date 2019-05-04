"""

General code to handle the ID of emitters
See e.g. Lofthouse et al. 2019, Fossati et al. 2019

Depend on proprietary code [cubex]

"""

import subprocess

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
        
        
