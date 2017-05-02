Some info on muse redux [MF, March 2016]
---------------------------------------

> Tested with eso muse pipeline 1.2.1
> Tested with eso muse pipeline 1.6.2
> Tested with Cubex 1.5
> Tested with ZAP v1

* Step 1

To run the eso pipeline, download data with calibrations from the ESO archive.
Dump all the data in folder ./Raw, including the *xml file.

If data require legacy calibrations, make sure there is a folder containing them. 
./staticcal
More info on legacy calibrations here: http://www.eso.org/sci/software/pipelines/#pipelines_table

These are updated regularly with pipeline releases, so check this before running the code!

Run the basic reduction like this:

>>muse=mp.ifu.muse.Muse()
>>xml=muse.redux_basic()                          

This will process the calibrations and the science with the eso pipeline,
producing cubes without sky subtraction for each exposure.

Step 1 is best run is data from a single OB or group of OBs taken within a 
few days apart. The calibrations used are stored in calplan.txt for inspection. 

Typical failures arise if not all calibrations are included or if Raw data are missing.
Also there are problems with old calibrations.

At the end, there should be a IMAGE_FOV and DATACUBE_FINAL for each science 
exposure in Proc. Check if they look ok. 


A. ESO REDUCTION
----------------

This step uses eso recipies to generate final cube. See MUSE manual for details.
Products from this step are also used in subsequent reduction/processing, so it 
is good norm to run this by default

At top level (above OB#) folders, run the following script

muse=mp.ifu.muse.Muse()
muse.eso_process()

This script crawls through the OB# folders, findings science frames.
It performs sky subtraction using ESO recipies, aligns exposures, and coadd them.

Possible points of failure are bad aligment of the exposures before coadding.
That may require manual interaction with the algnment step in muse_exp_align
such as tweaking the threshold. If all attempts fail, run muse_redux_eso.manual_align
for manual aligment.

E.g.

>> from mypython.ifu import muse_redux_eso as eso
>> eso.manual_align('align.sof')


Finally, one can register the absolute calibrations with code in muse_utils.py, 
where a procedure adjust_wcsoffset allows to reset the WCS zero point for cube and image.

E.g.

>> #correct the overall WCS solution 
>> x=165.994
>> y=129.625
>> ra=184.31369
>> dec=22.51928
>> data='esocombine/IMAGE_FOV_0001.fits'
>> mut.adjust_wcsoffset(data,x,y,ra,dec)
>> data='esocombine/DATACUBE_FINAL.fits'
>> mut.adjust_wcsoffset(data,x,y,ra,dec)



B. CUBEX_REDUCTION
------------------

This part of the scrip produces enhanced dataproducts with illumination correction and
sky subtraction following the procedures developed by S. Cantalupo. 
Cubex is currently private code, so please contact Cantalupo directly if interested in using it. 

* Step 2

Step 2 uses cubex from S. Cantalupo to perform illumination correction and sky subtraction. 

Call cubex from top level directory above OB reduction with

#run cubex process
>>muse=mp.ifu.muse.Muse()
>>muse.cubex_process()


First, the pipeline uses the eso reduction and the muse pipeline to align each expsoure to the 
reference, and produces a resampled cube on this frame. It's worth checking is the files 
DATACUBE_FINAL_LINEWCS_EXP* look sesnibly alligned to absolute reference wcs.  

Next, produces a first pass of illumination correction and sky subtraction for each science exposure in the OB#/Proc
subfolders. 

A bunch of data cubes are produced in the process:
 *_fix.fits is the cube with illumination correction
 *_skysub.fits is the suky sub cube
 *_white.fits is a projection of the processed cube

Next, the pipe makes a second pass of cubex to on the cubes with the 
updated wcs and clean source mask form the previous step.

After this, illumination corrected, sky subtracted, and WCS referenced cubes
are written to disk, one for each exposure. 
   *_fix2.fits for illumination correction
   *_skysub2.fits for sky sub cube
   *_white2.fits for FOV image

Once all the OBs have been processed up to this point, a final cube is reconstructed in the folder

cubexcombine/COMBINED_CUBE.fits

If there are obvious large scale artifact that should be mask, this can simply achieved by creating
a ds9 region file inside the OB#?proc folder, with same name as the mask from the pipeline but
extension .reg. The region file should be in ds9 format, with image coordinate.

Once the coadded cube is finally ready, the pipeline proceeds to perform a third iteration of 
cubex with illumination correction and skysubtraction that uses source masking from the high SN cubes.

In the end, a final coadded cube is reconstructed in 

cubexcombine/COMBINED_CUBE_FINAL.fits
cubexcombine/COMBINED_CUBE_MED_FINAL.fits

using both mean and median statistics. 


Optional: one can pass a ds9 region (image coordinates) containing regions to be masked when 
	  computing the sky normalisation in cubesharp. This is useful in case there are 
	  extended sources or extended line emission one wishes to mask

	  The syntax in this case is 
	  
	  >> muse.cubex_process(skymask='cubexcombine/skymask.reg')


C. LINE OPTIMISED REDUCTION
---------------------------

This set of utilities produces a final cube with illumination correction and 
skysubtraction that is optimised for cubes that are relatively empty in terms
of continuum sources, but that may have faint extended emission lines filling 
a good fraction of the field of view.

Sky surbtraction is performed with the public ZAP package 
https://github.com/ktsoto/zap
http://zap.readthedocs.io/en/latest/

At top level (above OB#) folders, run the following script, after ESO reduction 
as described in step A (it can be used also after cubex reduction, as it is independent code).

muse=mp.ifu.muse.Muse()
muse.line_process()

This perform illumination correction at IFU level as a function of coarse bins of wavelength.
Next, residual illumination correction on white images is performed at the stack level.

Finally, a zeroth order sky model is subtracted from the data and residual corrections
are perfomed with ZAP. 

A mask can be passed to avoid regions of extended flux emission.
The mask is passed in the form of a ds9 region file (image coordinate) that contains 
**good** regions to be used in PCA analysis (note the opposite behaviour with cubex, 
where a mask excludes pixels from sky normalisation).
In this case the syntax becomes:

   muse.line_process(skymask='linecombine/pca_region.reg')

An optional deepwhite keyword can be used to pass a deep white image for better masking
of continuum detected sources.

Next, the code comabined exposures into final cubes. Masks with IFU edges are automatically
propagated.

If there are obvious large scale artifacts that should be masked, this can simply achieved by creating a ds9 region file inside the OB#/Proc folder, with same name as the mask from the pipeline but
extension .reg. The region file should be in ds9 format, with image coordinate. In case no masks are found, as a rule, the code also searches for masks from the cubex reduction (above). These masks are not used if new masks are found. 


In the end, a final coadded cube is reconstructed in 

linecombine/COMBINED_CUBE_FINAL.fits
linecombine/COMBINED_CUBE_MED_FINAL.fits

Datacubes are produced with median and mean statistics. 

