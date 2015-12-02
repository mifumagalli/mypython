class Seds: 
    """ Define a class to handle SEDs """

    def __init__(self,name):
        self.name=name
        self.type=""
        self.collection=""
        self.filename=""
        self.wave=1
        self.flux=1
        

    def loadseds(self):
        from astropy.io import ascii
        from astropy.io import fits
        import os 
        import numpy as np

        path=os.environ['MYPYTHON']
        listsed = ascii.read(path+'seds/seds_lists.txt')
        
        #find index
        for index in range(len(listsed["sedname"])):
            if self.name in listsed["sedname"][index]:
                self.type=listsed["type"][index]
                self.collection=listsed["collec"][index]
                self.filename=listsed["sedname"][index]
                
        #now load 
        if self.collection == "calspec":
            hdulist = fits.open(path+'seds/'+self.type+'/'+self.collection+
                                '/'+self.filename)
            hdulist.info()
            tbdata = hdulist[1].data 
            self.wave=np.array(tbdata.field("WAVELENGTH")) #Ang
            self.flux=np.array(tbdata.field("FLUX"))  #erg/s/cm2/A

            ##convert into erg/s/cm2/Hz
            snu=(2.9979e10/self.wave/1e-8) 
            self.flux=self.flux*self.wave*self.wave/2.9979e10*1e-8
            
        else: 
            "No SED known with this name"
            exit()
                
                

                
                                          
                

