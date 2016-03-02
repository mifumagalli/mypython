"""
Make finders at telescope

"""

import astropy
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.wcs as wcs
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import mad_std

class Finder():


    def __init__(self,obj,off=None,fits=None):
        
        """

        obj Init the main object position 
            ra,dec either degree or hh,dd
        
        off same for offset star

        fits if set read from file 

        """

        #assume flots are in deg
        if isinstance(obj['ra'],float):
            self.obj=SkyCoord(obj['ra'],obj['dec'],unit='deg')
        else:
            self.obj=SkyCoord(obj['ra'],obj['dec'],unit=(u.hourangle, u.deg))
            
        #reference
        if(off):
            if isinstance(off['ra'],float):
                self.off=SkyCoord(off['ra'],off['dec'],unit='deg')
            else:
                self.off=SkyCoord(off['ra'],off['dec'],unit=(u.hourangle, u.deg))

            try:
                self.offmag=off['mag']
            except:
                self.offmag=0.0

        #deal with image
        if(fits):
            self.loadimg(fits)



    def loadimg(self,image):

        """
        Load a fits from file with name 'image'
        
        """
        
        #load
        data=fits.open(image)
        try:
            self.img=data[1].data
            self.w=wcs.WCS(data[1].header)
        except:
            self.img=data[0].data
            self.w=wcs.WCS(data[0].header)


        #check if East is to the left
        if(self.w.wcs.cd[0,0] > 0):
            print ('Flip image left/right')
            self.fliplr()

        #grab the scale
        self.scale=self.w.wcs.cd[1,1]*3600.



    def printfc(self,outname,annotate=True,cntr=4.0):

        """

        Generate finder

        outname: name of output
        annotate: if True, add info to finder
        cntr: number of sigma sky for display contrast  

        """


        #display image 
        fig = plt.figure(figsize=(8,10),dpi=80)
        ax = fig.add_subplot(111,projection=self.w)
        mean=np.median(self.img)
        stdev=mad_std(self.img) 
        imgplot = ax.imshow(self.img,origin='low',cmap='gray_r',
                            clim=(mean-cntr*stdev,mean+cntr*stdev))

        
        if(annotate):

            #make annotations
            #mark position obj 
            Ax,Ay = wcs.utils.skycoord_to_pixel(self.obj,self.w)
            p=patches.Circle((Ax,Ay),2.0/self.scale,edgecolor='red',
                             linewidth=3,facecolor='none')
            ax.add_patch(p)
  
            txt=self.obj.to_string(style='hmsdms',sep=':',precision=3)
            plt.text(0.5, 0.95,"Object: "+txt, horizontalalignment='center',
                     verticalalignment='center',transform=fig.transFigure,\
                         fontsize=22,color='red')

            #mark position offset 
            if(self.off):
                Bx,By = wcs.utils.skycoord_to_pixel(self.off,self.w)
                p=patches.Circle((Bx,By),2.0/self.scale,edgecolor='blue',
                                 linewidth=3,facecolor='none')
                ax.add_patch(p)
  
                txt=self.off.to_string(style='hmsdms',sep=':',precision=3)
                plt.text(0.5, 0.9,"Offset: "+txt, horizontalalignment='center',
                         verticalalignment='center',transform=fig.transFigure,
                         fontsize=22,color='blue')
                         
                         
                #offsets
                ra_offset = (((self.off.ra - self.obj.ra) * np.cos(self.obj.dec.to('radian'))).to('arcsec')).value
                dec_offset = (((self.off.dec - self.obj.dec)).to('arcsec')).value
                
                posang=self.obj.position_angle(self.off).degree

                strng=r'$\alpha$:{:8.2f}" $\delta$:{:8.2f}" PA:{:6.1f} m:{:4.1f}'.\
                    format(ra_offset,dec_offset,posang,self.offmag)
                plt.text(0.5, 0.84,strng, horizontalalignment='center',
                         verticalalignment='center',transform=fig.transFigure,
                         fontsize=25)
  
            #mark north east 
            plt.text(self.img.shape[0]/2.,self.img.shape[1],'North',
                     horizontalalignment='center',
                     verticalalignment='top',fontsize=25)
            plt.text(0.,self.img.shape[1]/2.,'East',
                     horizontalalignment='left',rotation='vertical',
                     verticalalignment='center',fontsize=25)
  
        #save
        plt.savefig(outname)


        
    def fliplr(self):
        
        
        """
        
        Flip the image left-right to have the North East orientation
        
        
        """
        
        #flip the actual array
        self.img=np.fliplr(self.img)

        #correct the WCS 
        self.neww=wcs.WCS(naxis=2)
        #redefine centre
        self.neww.wcs.crpix = [self.img.shape[0]/2.,self.img.shape[1]/2.]
            
        #reference position at centre
        self.neww.wcs.crval = self.w.wcs_pix2world(self.neww.wcs.crpix[0],self.neww.wcs.crpix[1],0)
            
        #copy over
        self.neww.wcs.cd = self.w.wcs.cd
        self.neww.wcs.ctype = self.w.wcs.ctype
        self.w=self.neww
        
        #flip wcs
        self.w.wcs.cd[0,0]=-1*self.w.wcs.cd[0,0]
        self.w.wcs.cd[1,0]=-1*self.w.wcs.cd[1,0]
            
