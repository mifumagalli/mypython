
"""

Make an RGB image


"""

import numpy as np

class RGB:
    

    def __init__(self,b,g,r):
        
        """
        Constructur
        
        """


        #store data
        self.b=b
        self.g=g
        self.r=r


        
    def makergb(self,rgbrange):
        
        """ Map the fits into RGB image """
    
        #create space for rgb
        self.img = np.zeros((self.b.shape[0], self.b.shape[1], 3), dtype=float)
        
        #map each channel and append 
        self.img[:,:,0]=self.logmap(self.r,rgbrange["R"])
        self.img[:,:,1]=self.logmap(self.g,rgbrange["G"])
        self.img[:,:,2]=self.logmap(self.b,rgbrange["B"])
        

    def logmap(self,fits,chanrange):
    
        """ 
        Map a fits into a log scale 
        Follow basic idea of this code:
        http://www.astrobetter.com/wiki/tiki-index.php?page=RGB+Images+with+matplotlib
        """
        #make a copy to avoid overwrite 
        data=np.nan_to_num(np.array(fits, copy=True))
        #interval in log
        factor = 1.*(np.log10(chanrange[1])-np.log10(chanrange[0]))
        
        #find data above/below range
        indx0 = np.where(data < chanrange[0])
        indx1 = np.where((data >= chanrange[0]) & (data <= chanrange[1]))
        indx2 = np.where(data > chanrange[1])
        data[indx0] = 0.0
        data[indx2] = 1.0
    
        #map log strecth in 0-1 range
        try :
            data[indx1] = (np.log10(data[indx1])-np.log10(chanrange[0]))/factor
        except :
            print "Error on log10!"
            
        return data
        
    def show(self):


        """
        Display the RGB image 
        
        """

        import matplotlib.pyplot as plt
        
        ax1=plt.subplot(1,1,1)
        is1=ax1.imshow(self.img,origin='lower')
        plt.show()
        
