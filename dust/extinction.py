class Extinction:
    """ Define an extinction class for reddening calculations """
    
    def __init__(self):
        #Coordinate at which ext is computed
        self.coord = -1
        self.ebv = -1
        self.correction = -1

    def getebv(self):
        
        """
        Given an RA/DEC position return EBV in tables
        Written by MF in Durham, Oct 2014. 
        
        """
        #query IRAS service 
        #http://irsa.ipac.caltech.edu/applications/DUST/
        #parse table to get the mean/std of 
        #SandF = Schlafly & Finkbeiner 2011 (ApJ 737, 103) [ext SandF mean,ext SandF std]
        #SFD   = Schlegel et al. 1998 (ApJ 500, 525) [ext SFD mean,ext SFD std]

        if (self.coord == -1):
            print("Error: coordinates must be set")
            exit()


        from astroquery.irsa_dust import IrsaDust
        ebv = IrsaDust.get_query_table(self.coord,section='ebv')
        self.ebv=ebv

        return ebv


    def specext(self):

        """ Compute the wavelength dependence correction factor assuming a MW extinction law

            The output is a table [wave,corr] where wave is wavelegths in Ang and corr
            is correction fraction that can be multiplied to observed spectrum to recover 
            intrinsic spectrum: f_int = f_obs * corr

        """
        import os
        import numpy as np
        from astropy.io import ascii

        #load MW curve 
        path=os.environ['MYPYTHON']
        data=ascii.read(path+'dust/data/MW_EXT_SLUG.dat')
        data.rename_column('col1', 'Wave(A)')
        data.rename_column('col2', 'Kappa(cm^2)')

        if (self.ebv == -1):
            #query for ebv
            self.getebv()
            
        #get the V band for AV
        from mypython.filters import filter
        myfil=filter.Filter()
        myfil.id=95
        myfil.loadtrans()
        
        #MW AV
        AV=self.ebv['ext SandF mean']*3.1
        #MW Kappa_V 
        kappaV=myfil.convolve(data['Kappa(cm^2)'],data['Wave(A)'])
        #find corrections and add it to table
        data['Correction']=np.exp((data['Kappa(cm^2)']/kappaV)*(AV/1.086))
        self.correction=data     
        
        #import matplotlib.pyplot as plt    
        #plt.plot(data['Wave(A)'],data['Correction'])
        #plt.show()

        return data
        
        
    def afilter(self,filterid):
        
        """ Return the extinction in mag in a given filter known to mypython """

        ##query the filter 
        #import mypython as mp
        #myfil=mp.filters.filter.Filter()
        #myfil.id=filterid
        #myfil.loadtrans()
        
        
    def atable(self):
        
        """ Return the extinction in mag in starndard filters """

        #query IRAS service 
        #http://irsa.ipac.caltech.edu/applications/DUST/
        
        if (self.coord == -1):
            print("Error: coordinates must be set")
            exit()

        from astroquery.irsa_dust import IrsaDust
        tab = IrsaDust.get_extinction_table(self.coord)
    
        return tab

    def info(self):
        
        """ Print a bunch of info """

        print("Position: ", self.coord.ra, self.coord.dec)
        print("EBV SDF: ", self.ebv['ext SFD mean'][0], self.ebv['ext SFD std'][0])
        print("EBV S&F: ", self.ebv['ext SandF mean'][0], self.ebv['ext SandF std'][0])

def example():

    """Example program """
    
    from astropy.coordinates import SkyCoord
    
    #set ext class
    ext=Extinction()
   
    #assign coordinates
    coord=SkyCoord('00h58m20s','+45d12m00s')
    ext.coord=coord

    #query for ebv 
    ext.getebv()

    #query for table 
    tab=ext.atable()

    #get a given A
    ext.afilter(154)
    
    #compute wave corrections
    ext.specext()

    #print info
    #ext.info()
 
