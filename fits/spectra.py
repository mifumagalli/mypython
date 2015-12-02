class Spectra:
    """Class to handle spectra """
    def __init__(self):
        self.wave=0
        self.npix=0
        self.flux=0
        self.error=0
        self.redshift=0

    def load(self,file):
        """Handle basic open procedure"""
        
        from astropy.io import fits
        import numpy as np

        #handle fits files 
        if file.endswith(".fits"):
            
            #open 
            specunit = fits.open(file)
            
            #now handle only simple case 
            #header wave, flux ext 0, err ext 1

            self.flux = specunit[0].data
            self.error = specunit[1].data
            self.npix=self.flux.size
            header = specunit[0].header
            self.setwave(header)
            specunit.close()


        else: 
            print "Format not supported"
            exit(1)
            
    def setwave(self,header):
        """ Handle the creation of wave based on header info """
        
        import numpy as np
        #case of linear lambda 
        if(header["CTYPE1"] == "AWAV"):
            self.wave=header["CRVAL1"]+header["CDELT1"]*(np.array(range(self.npix))+1-header["CRPIX1"])
            

    def vplot(self,wlines,vmin=-250,vmax=250):
        """Return velocity flux errors for a list of lines"""
        from astropy import constants as const
        import numpy as np

        #init list 
        vel=[]
        
        for line in wlines:
            tmpvel=np.array((self.wave-line*(1+self.redshift))/(line*(1+self.redshift))*const.c.to('km/s'))
            index=np.nonzero((tmpvel <= vmax) & (tmpvel >= vmin))
            vel.append({"vel":tmpvel[index],"flx":self.flux[index],"err":self.error[index],"wline":line})
            
            

        return vel
    
    
