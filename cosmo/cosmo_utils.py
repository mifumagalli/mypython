"""

Exapands cosmo utils from astropy.cosmo
Currently works only with FlatLambdaCDM

"""

import numpy as np
from astropy import units as u

class Cosmocal:


    def __init__(self,basic='Planck13'):
        
        """ 
        
        Initilaize basic cosmology from astropy 

        basic -> defines the astropy flavor of FlatLambdaCDM
        See astropy.cosmology manual for more

        """

        if('Planck13' in basic):
            from astropy.cosmology import Planck13 as apcosmo
        elif('WMAP5' in basic):
            from astropy.cosmology import WMAP5 as apcosmo
        elif('WMAP7' in basic):
            from astropy.cosmology import WMAP7 as apcosmo
        elif('WMAP9' in basic):
            from astropy.cosmology import WMAP9 as apcosmo
        else:
            raise ValueError('Cosmology flavor not recognised')
            

        self.cosmo=apcosmo


    def comovingvolume(self,zmax,zmin=False):

        """
        Compute a volume in flat cosmology

        zmax -> maximum redshift 
        zmin -> if set, compute volume between zmin/zmax

        """

        #Compute the comoving distance to redshift 
        disCom = self.cosmo.comoving_distance(zmax)
        #Now volume
        VolCom=((4.*np.pi/3.)*disCom.to('Gpc')**3)

        if(zmin):
            #Find the volume to zmin
            disCom_min = self.cosmo.comoving_distance(zmin)
            VolCom_min=((4.*np.pi/3.)*disCom_min.to('Gpc')**3)
            self.volume=VolCom-VolCom_min
        else:
            self.volume=VolCom

        return self.volume

            

        

