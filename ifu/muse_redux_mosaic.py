# muse_redux_mosaic.py
# 
# Created by Guido Consolandi. Modified and improved by Matteo Fossati (June 2017)
#
# This module is a compilation of routines to deal with mosaiced observations with MUSE
# A three step procedure is implemented. First make_expmap generates an exposure map 
# mosaic frame. Then read_expmap produces output_wcs files for the reconstruction of 
# individual exposures via muse_scipost. Lastly reconstruct_mosaic reads in the single
# exposures and combines them on the common grid defined by make_expmap.
# Routines to spatially trim a reconstructed mosaic and to paste a submosaic into a larger
# one are also available.
#


import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy.io import ascii
import string as str
from scipy import asarray as ar, exp
from scipy import ndimage
import glob


def powers_of2_finder(n):
    result = []
    binary = bin(n)[:1:-1]
    for x in range(len(binary)):
        if int(binary[x]):
            result.append(2**x)
    return (np.log2(result)).astype(int)


def make_expmap(RAG, DEG, NAXIS_x, NAXIS_y, imgfv, mosaicpos, outname):

  """
  STEP 1) Generate initial exposure map mosaic frame.
  
  Inputs are the center for the mosaic in RAG and DEG, the size in pixels,
  a MUSE IMAGE_FOV from which the header is taken and an ascii file listing
  the positions and sizes of each exposure. An example is given below.
  Recommended dimensions are 350 pixels for setups aligned N-S or E-W and 490
  pixels for exposures rotated by 45 degrees from these axes.
  The output file is written to outname
  
  #ID  dim[pix]    CRVAL1[deg] CRVAL2[deg]     comments
  1         350       175.80264 20.006542      MYMOSAIC1
  2         350       175.73517 19.97025       MYMOSAIC2
  3         350       175.7857  20.0145        MYMOSAIC3
  4         350       175.76889 20.024283      MYMOSAIC4
  5         350       175.75305 20.036153      MYMOSAIC5
  6         350       175.73672 20.040716      MYMOSAIC6
  7         490       175.73437 19.984127      MYMOSAIC7
  8         490       175.73583 20.00156       MYMOSAIC8
  
  """

  #Set the MOSAIC central coords and axis dimensions, 
  #start from an existent header of an IMAGEFOV and update the WCS
  
  print('Reading header of the file: '+imgfv)
  cubo=fits.open(imgfv)

  header_data = cubo[1].header
  header_data['CRPIX1'] = NAXIS_x/2
  header_data['CRPIX2'] = NAXIS_y/2
  header_data['CRVAL1'] = RAG
  header_data['CRVAL2'] = DEG

  header_data['NAXIS1']= NAXIS_x 
  header_data['NAXIS2']= NAXIS_y 
  header_data['EXTNAME']='BITMAP'

  empty_cube = np.zeros([NAXIS_y, NAXIS_x], dtype=np.float32)
  hdu_prim = fits.PrimaryHDU(header=cubo[0].header)  
  hdu_img  = fits.ImageHDU(empty_cube, header=header_data)           
  hdu_list = fits.HDUList([hdu_prim,hdu_img])        

  ###READ the FILE WITH OB number AND OB coords
  filecoord = np.array(ascii.read(mosaicpos, format='commented_header'))
  wcsmap = wcs.WCS(header_data)

  #go from crvals to pixel centers

  for i in range(np.size(filecoord['ID'])):
      
   coords = np.array([[filecoord['CRVAL1[deg]'][i] ,filecoord['CRVAL2[deg]'][i] ]], np.float_) 
   dim = filecoord['dim[pix]'][i]
   pixel_coo =wcsmap.wcs_world2pix(coords,1)
   xpix =  int(round(pixel_coo[0,0]))
   ypix =  int(round(pixel_coo[0,1]))
   
   ##Fill the bitmap
    
   hdu_list[1].data[ypix-int(dim)/2:ypix+int(dim)/2, xpix-int(dim)/2:xpix+int(dim)/2] +=  2**int(filecoord['ID'][i])
   
  hdu_list.writeto(outname, overwrite=True)
  
  print('Exposure map saved into file: '+outname)



def read_expmap(mosaicID, mosaic_file, nz=3680):    
    
    """
    STEP 2) Return square position of individual mosaic frame
    for OUTPUT_WCS structure.
    
    This function reads a master image (with WCS) which is filled
    with 2**mosaicID for each frame contributing to a given pixel
    The function returns the square position for a given frameID
    which can be fed into make_empty_cube from muse_redux_gc
    """
    
    print("Reading the bit map {0} ".format(mosaic_file))
    
    mastermap = fits.open(mosaic_file)
    
    xmin, ymin, xmax, ymax, mask = get_qbposition(mosaicID, mastermap)
    
    masterwcs = wcs.WCS(mastermap[1].header)

    xc, yc=(xmax+xmin)/2, (ymax+ymin)/2
    nx, ny= (xmax-xmin), ymax-ymin

    coo = np.array([[xc,yc]],np.float_)

    deg_coo =masterwcs.wcs_pix2world(coo,1)

    rag = deg_coo[0,0]
    deg = deg_coo[0,1]
    
    print("Mosaic ID {0:d} is mapped on a square area centered on RAG: {1:f} DEG: {2:f}".format(mosaicID, rag, deg))
            
    square = [rag,deg,nx,ny,nz] 
    
    return square



def reconstruct_mosaic(inname, outname, expmapfile, pix_coords, mode='all', ctype='mean', memmap=True):
    
    """
    STEP 3) Reconstruct the full mosaic structure by reading the individual cubes
    and positioning them at the right place.
    
    inname is the path of the input files, should contain variables to be filled internally
           e.g. '../DIR{0}/Proc/DATACUBE_FINAL_EXP{1}_pos_off_ill_bk_cl.fits'
    
    outname is the name of the output file, overwriteing is active
    
    expmapfile is the path to the exposure map file obtained with make_expmap
    
    pix_coords is a list of six elements the first four are the x1,x2,y1,y2 
    in pixels of the portion of the exposure map that you want to reconstruct.
    The fifth and sixth elements are lambdamin and lambdamax of the reconstruction
    
    mode is a switch to decide if you want a normal MUSE cube with data and variance (all)
    if you want two separate files for data and variance (separate) and if you want only 
    the data (data_only).
    
    ctype can be mean or wmean for weighted mean combine. Other methods will be implemented
    
    """
    
    xi=int(pix_coords[0])
    xf=int(pix_coords[1])
    yi=int(pix_coords[2])
    yf=int(pix_coords[3])
    
    lambdamin = pix_coords[4]
    lambdamax = pix_coords[5]
    
    print("Reading the Exposure mosaic map: {0} ".format(expmapfile))
    
    expmaphdu     = fits.open(expmapfile)
    expmaphdu_wcs = wcs.WCS(expmaphdu[1].header)
    
    crpixx, crpixy= (xi+xf)/2, (yi+yf)/2
    nx, ny= (xf-xi), yf-yi

    coo = np.array([[crpixx,crpixy]],np.float_)

    deg_coo =expmaphdu_wcs.wcs_pix2world(coo,1)

    crvalx = deg_coo[0,0]
    crvaly = deg_coo[0,1]
    
    infos_mos = [crvalx,crvaly,crpixx,crpixy,nx,ny] 
    
    print('Preparing an empty cube for selected area.... ')
    
    QBexpmap                       = make_empty_ima(infos_mos, inname.format(1,1))
    QBmosaic,     lamind1, lamind2 = make_empty_cube(infos_mos, lambdamin, lambdamax, inname.format(1,1))
    QBmosaic_var, lamind1, lamind2 = make_empty_cube(infos_mos, lambdamin, lambdamax, inname.format(1,1))
    
    
    print('Now read the bit map and fill the empty cube with data ')
          
    uniqbits = np.unique(expmaphdu[1].data[yi:yf,xi:xf]).astype(int)
    uniqbits = uniqbits[np.where(uniqbits >0)]
    
    print('The wavelength range {0} - {1} is extracted at pixels = {2} - {3} in the individual cubes'.format(lambdamin,lambdamax, lamind1, lamind2))
    print('The region defined by x_i={0}, x_f={1}, y_i={2}, y_f={3} is composed of {4} different regions of overlap'.format(xi,xf,yi,yf,np.size(uniqbits)))
    print('It has CRPIX1= {0} and CRPIX2 = {1} corresponding to CRVAL1={2:7.4f} and CRVAL2={3:7.4f}'.format(crpixx, crpixy, crvalx, crvaly))
    print(' --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---')
    
    for i in range(len(uniqbits)):                                  
        
        yreg, xreg= np.where(expmaphdu[1].data[yi:yf,xi:xf] == uniqbits[i])
        
        contrib_exp =  powers_of2_finder(uniqbits[i])   
        
        if len(contrib_exp) > 1: 
             print('Region number {0}, coded as {1}, is the combination of the Mosaic QBs:'.format(i+1, uniqbits[i]), contrib_exp)
        else: 
             print('Region number {0}, coded as {1}, is covered only by the Mosaic QB number '.format(i+1, uniqbits[i]), contrib_exp)
        
        #Now find relevant datacubes
        explist = []
        obslist = []
        for k in range(len(contrib_exp)):
          files = glob.glob(inname.format(contrib_exp[k], '*'))
          for j in range(len(files)):
            explist.append(files[j])
            obslist.append(contrib_exp[k])
        
        if ctype=='mean' or ctype=='wmean':
           QBmosaic, QBmosaic_var, QBexpmap  = feedQB(QBmosaic, QBmosaic_var, QBexpmap, explist, obslist, yreg, xreg, pix_coords, lamind1, lamind2,  expmaphdu, uniqbits[i], ctype=ctype, memmap=memmap)
        else:
           'Combine type not understood'
           return   
        
    for wave in range(len(QBmosaic[1].data[:,0,0])):
      nans = (QBmosaic[1].data[wave,:,:] == 0) | (np.isinf(QBmosaic[1].data[wave,:,:]))
      QBmosaic[1].data[wave,nans] = np.nan
      nans = (QBmosaic_var[1].data[wave,:,:] == 0) | (np.isinf(QBmosaic_var[1].data[wave,:,:]))
      QBmosaic_var[1].data[wave,nans] = np.nan
    
    QBexpmap.writeto(outname.replace('.fits', '_expmap.fits') ,overwrite=True)
    
    if mode == 'separate':
        
       #write output       
       outvar= outname.replace('.fits', '_var.fits')       
       QBmosaic.writeto(outname,overwrite=True)
       QBmosaic_var.writeto(outvar,overwrite=True)   
    
    elif mode == 'all':
               
       h0 = QBmosaic[0].header
       hdat = QBmosaic[1].header
       hvar = QBmosaic_var[1].header  
       
       hd_prim = fits.PrimaryHDU(header=h0)
       hd_im =fits.ImageHDU(QBmosaic[1].data)
       hd_im.header=hdat
       hd_var = fits.ImageHDU(QBmosaic_var[1].data)
       hd_var.header = hvar
       
       hdlist =fits.HDUList([hd_prim, hd_im, hd_var]) 
       hdlist.writeto(outname,clobber=True)
       

    elif mode == 'data_only':
               
       QBmosaic.writeto(outname,clobber=True)
       QBexpmap.writeto('exptemp.fits',clobber=True)
     
    else:
     
        print('Mode not understood')

def mosaic_clean(QBinput, QBout, size=1.0):
    
    #open the data
    data=fits.open(QBinput, memmap=False)

    #compress into image
    image=np.nanmedian(data[1].data,axis=0)
    nz,ny,nx=(data[1].data).shape

    #mask edges
    badmask=np.ones((ny,nx))
    badmask[np.isfinite(image)]=0
    badmask=ndimage.gaussian_filter(badmask,size)
    badmask[np.where(badmask > 0)]=1
    okmask = np.copy(badmask)
    okmask[badmask==0] = 1
    okmask[badmask==1] = np.nan
    
    
    for ii in range(nz):
       data[1].data[ii,...] *= okmask
       data[2].data[ii,...] *= okmask
    
    data.writeto(QBout, overwrite=True)
    
    
def mosaic_extract(QBinput, QBout, limits):

    '''
    This function allows to extract sub-QBs from a mother cube
    Highest limit is not included as per Python convention 
    '''
        
    qbin = fits.open(QBinput)
    wcs_map = wcs.WCS(qbin[1].header)
    
    header_new = qbin[1].header
    header_prim =  qbin[0].header
    
    xmin = limits[0]
    xmax = limits[1]
    ymin = limits[2]
    ymax = limits[3]
    
    xdim = xmax - xmin
    ydim = ymax - ymin
    
    crpix1 = xmin+xdim/2
    crpix2 = ymin+ydim/2
    
    coo = np.array([[crpix1,crpix2,1]],np.float_)

    deg_coo =wcs_map.wcs_pix2world(coo,1)

    crvalx = deg_coo[0,0]
    crvaly = deg_coo[0,1]
       
    print('Extracting a sub-cube of {0}x{1} pixels'.format(xdim, ydim))
    
    nlam = header_new['NAXIS3']
    header_new['CRPIX1'] = xdim/2
    header_new['CRPIX2'] = ydim/2
    header_new['CRVAL1'] = crvalx
    header_new['CRVAL2'] = crvaly
    header_prim['RA']  = crvalx
    header_prim['DEC'] = crvaly    
    header_prim['XI'] = xmin
    header_prim['XF'] = xmax    
    header_prim['YI'] = ymin
    header_prim['YF'] = ymax    
    
    empty_cube = np.zeros([nlam, ydim, xdim], dtype=np.float32)
    empty_var = np.zeros([nlam, ydim, xdim], dtype=np.float32)
    
    hdu_pri = fits.PrimaryHDU(header=header_prim)
    hdu_im  = fits.ImageHDU(empty_cube)
    hdu_im.header= header_new
    hdu_var = fits.ImageHDU(empty_var)
    hdu_var.header= header_new
    newQB = fits.HDUList([hdu_pri,hdu_im,hdu_var])
    
    newQB.info()
    
    newQB[1].data[:,:,:]=qbin[1].data[:,ymin:ymax,xmin:xmax]
    newQB[2].data[:,:,:]=qbin[2].data[:,ymin:ymax,xmin:xmax]
    
    
    newQB.writeto(QBout, clobber=True)
    
    qbin.close()


def mosaic_insert(mastercube, childcubes, outputcube):
    
    '''
    Nothing fancy, just pastes subQBs into the mastercube.
    Needs XI, XF, YI, TF in primary header
    '''
    
    #im_gand = np.genfromtxt(child_file,skip_header=1,dtype="S",usecols=0) 
    #xi, xf, yi, yf = np.genfromtxt(child_file,skip_header=1,usecols=(1,2,3,4),unpack=True)

    mother = fits.open(mastercube)
    
    for i in range(np.size(childcubes)):
           
           child = fits.open(childcubes[i])
           
           try:
             xi = child[0].header['XI'] 
             xf = child[0].header['XF'] 
             yi = child[0].header['YI'] 
             yf = child[0].header['YF'] 
          
             print('Pasting {0} at coordinates {1} {2} {3} {4} of master cube'.format(childcubes[i], xi, xf, yi, yf))
           
             mother[1].data[:,yi:yf,xi:xf] = child[1].data
             mother[2].data[:,yi:yf,xi:xf] = child[2].data
           
           except:
             
             print('Coordinate Error. Skipping {0} '.format(childcubes[i]))
           
    mother.writeto(outputcube,clobber=True)


#-------------------------------------------
#           INTERNAL ROUTINES
#-------------------------------------------

def get_qbposition(mosaicID, expmaphdu):
    
    #Find where mosaic XX is contributing data in the master map
    
    expmaphdu_wcs = wcs.WCS(expmaphdu[1].header)

    uniqbits = np.unique(expmaphdu[1].data).astype(int)
    contributing = np.zeros_like(uniqbits)
    
    for i in range(len(uniqbits)):
          
         contrib_exp =  powers_of2_finder(uniqbits[i])  
         if mosaicID in contrib_exp:
             contributing[i] += 1
                    
    selectedbits = uniqbits[(contributing == 1)]

    for kk in range(len(selectedbits)):
            
        if kk == 0: 
          target = (expmaphdu[1].data == selectedbits[kk]) 
        else: 
          target = target | (expmaphdu[1].data == selectedbits[kk])

    yy, xx =np.where(target==True)

    xmin_onmap=np.amin(xx) 
    ymin_onmap=np.amin(yy)
    xmax_onmap=np.amax(xx) 
    ymax_onmap=np.amax(yy)
        
    mask = np.isnan(expmaphdu[1].data[ymin_onmap:ymax_onmap+1,xmin_onmap:xmax_onmap+1])
    
    return xmin_onmap, ymin_onmap, xmax_onmap+1, ymax_onmap+1, mask


def make_empty_cube(infos, lambdamin, lambdamax, musefile):
    
    #Infos is an array of [CRVAL1, CRVAL2, NAXIS1, NAXIS2]
    
    prim_head = fits.getheader(musefile, 0) 
    header_data = fits.getheader(musefile, 1)
    
    RAG=infos[0]
    DEG=infos[1]
    
    allwaves = header_data['CRVAL3'] + header_data['CD3_3'] * np.arange(header_data['NAXIS3'])
    
    okwaves = (allwaves>lambdamin-1.24) & (allwaves<lambdamax+1.24)
    
    lamind1  =  np.argmax(okwaves)
    lamind2  =  len(okwaves)-np.argmax(okwaves[::-1])-1
    
    npix_lam = lamind2-lamind1+1

    header_data['CRPIX1'] = infos[2]
    header_data['CRPIX2'] = infos[3]
    header_data['CRVAL1'] = RAG
    prim_head['RA'] = RAG
    header_data['CRVAL2'] = DEG
    prim_head['DEC'] = DEG
    header_data['CRVAL3'] = lambdamin
    
    #header_data['NAXIS1']= infos[2] 
    #header_data['NAXIS2']= infos[3]
    #header_data['NAXIS3']= npix_lam
    header_data['EXTNAME']='MOSAIC'

    empty_cube = np.zeros([int(npix_lam), int(infos[5]), int(infos[4])], dtype=np.float32)

    hdu_pri = fits.PrimaryHDU(header=prim_head)
    hdu_im  = fits.ImageHDU(empty_cube, header=header_data)
    HDU_2fill = fits.HDUList([hdu_pri,hdu_im])
    
    return HDU_2fill, lamind1, lamind2  


def make_empty_ima(infos, musefile):
    
    prim_head = fits.getheader(musefile, 0) 
    header_data = (fits.getheader(musefile, 1))
    
    RAG=infos[0]
    DEG=infos[1]
    
    header_data['CRPIX1'] = infos[2]
    header_data['CRPIX2'] = infos[3]
    header_data['CRVAL1'] = RAG
    prim_head['RA'] = RAG
    header_data['CRVAL2'] = DEG
    prim_head['DEC'] = DEG
    
    #header_data['NAXIS1']=infos[2] 
    #header_data['NAXIS2']= infos[3]
    header_data['EXTNAME']='MOSAIC'
    header_data.remove('NAXIS3') 
    header_data.remove('CRVAL3')
    header_data.remove('CRPIX3')
    header_data.remove('CD3_3')

    empty_cube = np.zeros([int(infos[5]), int(infos[4])], dtype=np.float32)

    hdu_pri = fits.PrimaryHDU(header=prim_head)
    hdu_im  = fits.ImageHDU(empty_cube, header=header_data)
    HDU_2fill = fits.HDUList([hdu_pri,hdu_im])
    
    return HDU_2fill 

def get_bounds(xiqb, yiqb, xfqb, yfqb, ps, min_x, max_x, min_y, max_y):

    deltax_i = xiqb - ps[0]
    deltax_f = xfqb - ps[1]
    deltay_i = yiqb - ps[2]
    deltay_f = yfqb - ps[3]
    
    if deltax_i > 0 : xi=0
    else: xi=0-deltax_i     
    
    if deltay_i > 0 : yi = 0
    else: yi=0-deltay_i     
    
    if deltax_f > 0 : xf = xfqb - deltax_f -xiqb
    else: xf=xfqb-xiqb+1     
    
    if deltay_f > 0 : yf = yfqb - deltay_f -yiqb
    else: yf=yfqb-yiqb+1   
    
    dxi = min_x - xi
    dxf = max_x - xf
    dyi = min_y - yi
    dyf = max_y - yf
    
    if dxi > 0 : xi = min_x
    if dxf < 0 : xf = max_x+1
    if dyi > 0 : yi = min_y
    if dyf < 0 : yf = max_y+1
    
    return xi, xf, yi, yf

def feedQB(mosdata, mosvar, mosexp, filelist, obslist, yreg, xreg, pix_coords, lamind1, lamind2, expmaphdu, magicnum, ctype='mean', memmap=True):
    
    npix_lam = lamind2-lamind1+1
        
    expformean = mosdata[1].data[0:npix_lam,yreg,xreg]*0
    
    for kk in range(np.size(obslist)):
      
      #print('Adding file to the mosaic: {0}'.format(filelist[kk]))
         
      qb = fits.open(filelist[kk], memmap=memmap)
      
      xiqb, yiqb, xfqb, yfqb, msk  = get_qbposition(obslist[kk], expmaphdu)
            
      indices= np.where(expmaphdu[1].data[yiqb:yfqb,xiqb:xfqb] == magicnum)  
      min_x  = np.amin(indices[1])
      min_y  = np.amin(indices[0])
      max_x  = np.amax(indices[1])
      max_y  = np.amax(indices[0])     
      
      xi,xf,yi,yf = get_bounds(xiqb, yiqb, xfqb, yfqb, pix_coords, min_x, max_x, min_y, max_y)          
      
      msk[yi:yf,xi:xf]=True
      
      good = (expmaphdu[1].data[yiqb:yfqb,xiqb:xfqb] == magicnum) & msk
      
      goody, goodx =  np.where(good == True)

      this_expmap = qb[1].data[lamind1:lamind2+1, goody, goodx]
      nans = np.isnan(this_expmap)
      this_expmap[np.logical_not(nans)] = 1
      this_expmap[nans] = 0
      
      expformean += this_expmap
      
      mosexp[1].data[yreg, xreg] = np.nanmedian(expformean, axis=0) 
      
      if ctype=='mean':
        
        mosdata[1].data[0:npix_lam, yreg, xreg] = np.nansum([mosdata[1].data[0:npix_lam, yreg, xreg], qb[1].data[lamind1:lamind2+1, goody, goodx]], axis=0)
        mosvar[1].data[0:npix_lam, yreg, xreg]  = np.nansum([mosvar[1].data[0:npix_lam, yreg, xreg], qb[2].data[lamind1:lamind2+1, goody, goodx]], axis=0)
      
      elif ctype=='wmean':
        
        data = qb[1].data[lamind1:lamind2+1, goody, goodx]
        weight = 1./(qb[2].data[lamind1:lamind2+1, goody, goodx]) #1/variance
            
        mosdata[1].data[0:npix_lam, yreg, xreg] = np.nansum( [mosdata[1].data[0:npix_lam, yreg, xreg],data*weight], axis=0)
        mosvar[1].data[0:npix_lam, yreg, xreg] = np.nansum( [mosvar[1].data[0:npix_lam, yreg, xreg],weight], axis=0)
      
      qb.close() 
    
    if ctype=='mean':
      mosdata[1].data[:, yreg, xreg] /= expformean
      mosvar[1].data[:, yreg, xreg]  /= expformean**2
    elif ctype=='wmean':
      mosdata[1].data[:, yreg, xreg] /= mosvar[1].data[:, yreg, xreg]
      mosvar[1].data[:, yreg, xreg]  = 1./mosvar[1].data[:, yreg, xreg]

    
    return mosdata, mosvar, mosexp

   
