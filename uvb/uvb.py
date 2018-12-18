import os
import numpy as np
from matplotlib import pyplot as plt
import scipy
from scipy.integrate import quad
from scipy.interpolate import interp1d


class UVB:

    """
    Class to perform operations on UVB models
    
    """

    def __init__(self,redshift,tag="MH2015"):

        """
        
        Init the UVB at a given redshift
        
        """
        
        self.redshift = redshift
        self.tag = tag
        self.load()

    def load(self):

        """

        Utility function to load the desired UVB at a given redshift

        """

        #select how to parse
        if("MH2015" in self.tag):
            uvbfile=os.getenv("MYPYTHON")+"uvb/MH2015.txt"
            self.cubareader(uvbfile)    
        elif("HM2012" in self.tag):
            uvbfile=os.getenv("MYPYTHON")+"uvb/HM2012.txt"
            self.cubareader(uvbfile)    

        
    def cubareader(self,filename):

        #open
        uvbfile=open(filename)

        #select counter
        counter = 0
        #parse
        wavelength=np.array([])
        juvb=np.array([])
        for line in uvbfile:
            if("#" not in line):
                counter=counter +1
                fields=line.split()
                if(counter == 1):
                    #grab the redshifts
                    redshifts=np.array(fields).astype(float)
                    #find inxed of nearest neight in redshift
                    zind=np.argmin(np.abs(redshifts-self.redshift))

                else:
                    wavelength=np.append(wavelength,np.float(fields[0]))
                    juvb=np.append(juvb,np.float(fields[zind+1]))

        uvbfile.close()

        #store
        self.wave=wavelength #AA
        self.juvb=juvb #erg cm−2 s−1 Hz−1 sr−1

    def jintegral(self,minwave,maxwave):

        """
        
        Given two wave in A, performs intergation of J_nu d_nu

        """

        #from wave to freq
        clight=2.99792458e10
        freq=np.flip(clight/(self.wave*1e-8)) # go to Hz
        jnuflip=np.flip(self.juvb)
        minnu=clight/(maxwave*1e-8)
        maxnu=clight/(minwave*1e-8)

        #construct interpolation array 
        jintegral=interp1d(freq,jnuflip)    
        jint=quad(jintegral,minnu,maxnu)
        return jint
        
    def jrate(self,minwave,maxwave):

        """
        
        Given two wave in A, performs intergation of J_nu/h nu d_nu

        """

        #from wave to freq
        clight=2.99792458e10
        hplank=6.6260755e-27
        freq=np.flip(clight/(self.wave*1e-8)) # go to Hz
        jnuflip=np.flip(self.juvb)/(hplank*freq)
        minnu=clight/(maxwave*1e-8)
        maxnu=clight/(minwave*1e-8)

        #construct interpolation array 
        jintegral=interp1d(freq,jnuflip)    
        jint=quad(jintegral,minnu,maxnu)
        return jint      
            
if __name__ == "__main__":

    #load UVB
    red=3.0
    uvb=UVB(red,tag="HM2012")

    #compute integral
    jintegral,jerror=uvb.jintegral(912./20,912.)
    

    #compute rate and SB
    jrate,jerror=uvb.jrate(912./20,912.)
    RateLya= 0.65*jrate
    SBLya=RateLya*6.6260755e-27*2.99792458e10/(1215.67*1e-8)/(1+red)**4/4.25e10 

    SBHa=SBLya/8.7
    
    print(jrate,RateLya,SBLya,SBHa)
    
