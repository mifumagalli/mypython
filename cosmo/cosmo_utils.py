"""

Exapands cosmo utils from astropy.cosmo
Currently works only with FlatLambdaCDM

"""

import numpy as np
from astropy import units as u
from astropy.table import Table
from scipy.interpolate import interp1d
from scipy.integrate import quad
import matplotlib.pyplot as plt
import os

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

        self.cosmotag=basic
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

    def Dz(self,redshift):

        """

        Compute Dz at a given redshift

        """

        OmegaM=self.cosmo.Om0
        OmegaL=(1-OmegaM)
        z=redshift

        Ez=np.sqrt(OmegaL+OmegaM*(1+z)**3)
        OmegaL_z=OmegaL/Ez**2
        OmegaM_z=1-OmegaL_z
        g0=5/2*OmegaM/(OmegaM**(4./7.)-OmegaL+(1.+0.5*OmegaM)*(1+OmegaL/70.))
        gz=5/2*OmegaM_z/(OmegaM_z**(4./7.)-OmegaL_z+(1+0.5*OmegaM_z)*(1+OmegaL_z/70))
        Dz=gz/(g0*(1+z))

        return Dz
        
    def dmC(self,mass200,redshift):

        """
        
        Compute a concentration parameter using interpolation of results by
        http://adsabs.harvard.edu/abs/2016MNRAS.460.1214L
        
        mass200 -> log10(m200/Msun) 
        redshift -> redshift 

        """

        filepath=os.environ['MYPYTHON']
        
        if(self.cosmotag is 'Planck13'):
            incmz=filepath+'/cosmo/Planck_cMz.dat'
            insigma0=filepath+'/cosmo/Planck_sigma0.dat'            
        elif(self.cosmotag is 'WMAP9'):
            incmz=filepath+'/cosmo/WMAP9_cMz.dat'
            insigma0=filepath+'/cosmo/WMAP9_sigma0.dat'            
        else:
            raise ValueError('Not ready for this cosmology')

        #convert in log (1e10 M/h)
        mass=mass200+np.log10((self.cosmo.H(0)/100.).value)-10
        
        #load data
        #col1: z=0 masses; col2: z=0 concentrations
        #The mass units are log10(M200 [10^10 Msun/h]); c=r200/r_s.
        cmz = Table.read(incmz,format='ascii')

        #mass variance extrapolated to z=0.
        #col1: mass are again [10^10 Msun/h], col2: sigma0(M,z=0) 
        sigma0 = Table.read(insigma0,format='ascii')

        #construct interpolation function sigma0 vs Mass 
        fsigma0=interp1d(sigma0['col1'],sigma0['col2'],bounds_error=False,fill_value='extrapolate')
        #compute nu at z=0 matched to the mass in cz0
        nu0=1.686/fsigma0(cmz['col1'])

        #construct concentration vs nu0 function
        fcnu0=interp1d(nu0,cmz['col2'])

         #get nu at the redshift and mass of interest 
        nuzM=1.686/(fsigma0(mass)*self.Dz(redshift))         
        
        #get c at desired mass
        cMz=fcnu0(nuzM)

        #m=np.arange(-5,5,0.1)
        #plt.scatter(cmz['col1'],cmz['col2'])
        #plt.plot(m,fcnu0(1.686/(fsigma0(m)*self.Dz(0.))))
        #plt.show()
        #exit()
        
        return cMz
    
    def meandensity(self,redshift):

        """
        Compute the mean baryon density at a given redshift 
        following Bryan and Norman 1998

        return values in g/cm^3 and 1/cm^3

        """

        #matter density of the universe today 
        rho_crit_zero=self.cosmo.critical_density(0.)        
        #in baryon
        rho_mean_zero=rho_crit_zero*self.cosmo.Ob0
        rho_mean=rho_mean_zero.value*(1+redshift)**3

        #in number density (25% He, 75% H)
        mp=1.67262e-24 #proton mass in g
        X=0.75
        Y=1-X
        mu=1./(2.*X+3./4.*Y) #mean molecular w.
        nh_mean=rho_mean/(mu*mp)
        
        return rho_mean, nh_mean
    
    def r200(self,mass,redshift):

        """
        Compute R200 (in kpc) for a given redshift and mass (M200)
        
        mass -> log10(M/Msun)
        
        """

        #some constants
        msun=1.99e33 #sun mass in g
        cm2pc=3.086e18 #pc in cm

        #desired overdensity
        Delta=200.

        #get R200 in kpc
        r200cube=3.*((10**mass)*msun)/(4.*np.pi*self.cosmo.critical_density(redshift).value*Delta)
        r200=r200cube**(1./3)/cm2pc/1e3
        
        return r200

    
    def getNFW(self,radius,mass,redshift):

        """
        
        Evaluate a NFW profile given a mass and redshift (density g/cm^3)

        mass -> log10(M200/Msun)
        radius -> radius in kpc for evaluation

        """

        #get R200 in kpc
        r200kpc=self.r200(mass,redshift)
        
        #get concetration
        dmC=self.dmC(mass,redshift)
        
        #characteristic radius
        rskpc=r200kpc/dmC        

        #define x=r/rs
        x=radius/rskpc

        #now get the characteristic density
        fMx=lambda x: 3.*x**2/(x*(1+x)**2)
        fM200,efM200=quad(fMx,1e-8,dmC)
        
        rhos=(200.*self.cosmo.critical_density(redshift)*dmC**3/fM200).value
  
        #compute NFW
        rhonfw=4*rhos/(x*(1+x)**2)
        
        return rhonfw
        
