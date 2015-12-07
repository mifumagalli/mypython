Some info on muse redux [MF, Dec 2015]
-------------------------------------

> Tested with eso muse pipeline 1.2.1
> Tested with Cubex 1.4

* Step 1

To run the eso pipeline, download data with calibrations from the ESO archive.
Dump all the data in folder ./Raw, including the *xml file.

Run the basic reduction like this:

>>muse=mp.ifu.muse.Muse()
>>xml=muse.redux_basic()                          

This will process the calibrations and the science with the eso pipeline,
producing cubes without sky subtraction for each exposure.

Step 1 is best run is data from a single OB or group of OBs taken within a 
few days apart. The calibrations used are stored in calplan.txt for inspection. 

Typical failures arise if not all calibrations are included or if Raw data are missing

At the end, there should be a IMAGE_FOV and DATACUBE_FINAL for each science 
exposure in Proc. Check if they look ok. 
 
* Step 2

Step 2 uses cubex from S. Cantalupo to perform illumination correction and sky subtraction. Cubex is currently private code, so please contact Cantalupo directly if interested in using it. 

First, prepare a list "radec.txt" of 

rag1 deg1
rag2 deg2

of sources in the filed that can be used to refine the wcs position. 

Next, call cubex with

>>refcat='radec.txt'
>>muse.cubex_process(refcat,xml)

with xml from step 1

This produces a first pass of illumination correction and sky subtraction for each science exposure.

A bunch of data cubes are produced in the process:
 *_fix.fits is the cube with illumination correction
 *_skysub.fits is the suky sub cube
 *_white.fits is a projection of the processed cube

The pipe next proceeds in computing the optimal aligment of the exposures, 
written in radecoffsets.txt. New pixel table *_off.fits are written to disk.
And a cube/image _off.fits is written for each exposures. Now, you can 
check whether the wcs are properly aligned on the white images.

OPTIONAL: At this step, the cube is trimmed to the overlapping region in the 
exposure considered in this batch. If one is handling a large mosaic spanning more than one OB and it is reducing invividual OBs, the following steps should be followed instead.

 . After the first pass of cubex and wcs offsets are found, kill the pipleine.
 . Process all the OBs up to this first step, and combine all the cubes 
   with the eso recipe muse_exp_combine
 . With the reference cube in hand "refcube.fits", restart the pipeline
   with the call

   >>refcat='radec.txt'				
   >>refcube='refcube.fits'			
   >>muse.cubex_process(refcat,xml,refcube=refcube)

   The header of refcube (containing wcs for the entire mosaic), are then 
   used by cubex to ensure that all the data are combined on a uniform 
   grid.

Next, the pipe makes a second pass of cubex to reconstruct the cube with the 
updated wcs. This step is required to avoid resampling the cube more than once.
In fact, this is just a repeat of step 1, but with a single resampling on a 
final WCS reference. 

After this, illumination corrected, sky subtracted, and WCS referenced cubes
are written to disk, one for each exposure. 
   *_off_fix.fits for illumination correction
   *_off_skysub.fits for sky sub cube
   *_off_white.fits for FOV image

OPTIONAL:

   Once all the OBs have been processed up to this point, one may want to 
   combine all the cubes [see below] to produce a high sn cube "highsn.fits" 
   that can be fed back in cubex for optimal sky subtraction and masking.

   To perfom this optional (but recommended) step, just restart the pipeline
   with the call

   >>refcat='radec.txt'				
   >>refcube='refcube.fits'			
   >> highsn='highsn.fits'   						
   >>muse.cubex_process(refcat,xml,refcube=refcube,highsn=highsn)
    
   A second set of cubes "2" are written to disk.
  
* Step 3. At this point, all the data for each science exposure are processed in 
  a final cube. The last step is simply the coaddition of all the cubes.
  With the cubes being written on the same wcs reference, this step is trivial
  as it involves a simple sum of the cubes.

  In cubex, this can be done with the following command

  import mypython 
  from mypython.ifu import muse_redux_cubex as cbs
  cbs.combine_cubes("cubeall.in","maskall.in","regionsall.in")
  
  where cubeall.in is a list of the skysub cubes
        maskall.in is a list of the associated SliceEdgeMask.fits
	regionsall.in [optional] is a list of DS9 regions with boxes for 
	              additional masking of defects, etc..



  
  

 











