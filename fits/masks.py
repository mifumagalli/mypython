class Mask:
    """Class to handle masks in fits files """

    def __init__(self):
        #specify y/x dimen of mask
        dim=[1,1] #NB (y,x)
        data=1
        xdata=1
        ydata=1

    def maskregions(self,regions):
        """For an input regions list [ds9 format], create a mask by selecting regions
           inputreg is a ds9 region file  
        """
        import numpy as np

        #initialize the mask
        self.maskinit()

        #loop over regions 
        for rr in regions:
            #parse them 
            sprr=rr.split("(")
            regtype=sprr[0]            
            regshap=sprr[1].split(")")
            regcomm=regshap[1]
            regshap=regshap[0].split(",")
            #fill regions
            self.fillmask(regtype,regshap)
        

        #import matplotlib.pyplot as plt
        #import matplotlib.image as mpimg
        #plt.imshow(self.data)
        #plt.show()
     

    def maskinit(self):
        """ Initialize a mask with zeros everywhere """

        import numpy as np
       
        print "Initializing mask..."
        
        #data region 
        self.data=np.zeros((self.dim[0],self.dim[1]))
      
        #x/y image
        self.ydata=np.zeros((self.dim[0],self.dim[1]))
        for cc in range(self.dim[1]):
            self.ydata[:,cc]=np.arange(self.dim[0])+1


        self.xdata=np.zeros((self.dim[0],self.dim[1]))
        for cc in range(self.dim[0]):
            self.xdata[cc,:]=np.arange(self.dim[1])+1
        

        #import matplotlib.pyplot as plt
        #import matplotlib.image as mpimg
        #plt.imshow(self.xdata)
        #plt.show()
        
            
    def fillmask(self,regtype,regshap):
        """ Return index of pixels masked for a given shape """
        
        if regtype == 'circle':
           
            #retunr index in circle (xc,yc,rad)
            import numpy as np
            distance=np.sqrt(np.square(self.xdata-float(regshap[0]))+
                             np.square(self.ydata-float(regshap[1])))
            index=np.nonzero(distance < float(regshap[2]))
            self.data[index]=self.data[index]+1


        elif regtype == 'box':
            
            #return index in box (x y width height angle)
            import numpy as np
            from scipy.ndimage import interpolation
        
           

            #construct box with no rotation
            thisbox=np.zeros((self.dim[0],self.dim[1]))
            index=np.nonzero((self.xdata > self.dim[1]*0.5-0.5*float(regshap[2]))
                             &(self.xdata < self.dim[1]*0.5+0.5*float(regshap[2]))
                             &(self.ydata > self.dim[0]*0.5-0.5*float(regshap[3]))
                             &(self.ydata < self.dim[0]*0.5+0.5*float(regshap[3])))
            thisbox[index]=thisbox[index]+1
            thisbox=interpolation.rotate(thisbox,-float(regshap[4]),reshape=False,order=0)
            thisbox=interpolation.shift(thisbox,[float(regshap[1])-self.dim[0]*0.5,
                                                 float(regshap[0])-self.dim[1]*0.5],order=0)
            
            self.data=self.data+thisbox
            
            #import matplotlib.pyplot as plt
            #import matplotlib.image as mpimg
            #plt.imshow(thisbox)
            #plt.show()
            #exit()
    
        else: 
            print "No region "+regtype+" known..."
            exit()


    def write(self,outname):
        """Write the mask output"""

        from astropy.io import fits
        
        #PrimaryHDU object to encapsulate the data:
        hdu = fits.PrimaryHDU(self.data)
        #write
        hdu.writeto(outname,clobber=True)



        
