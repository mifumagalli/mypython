#Class to evaluate the luminosity function 
import numpy as np

class LumFun:


    def __init__(self,redshift=3.5,whichlf='Buow15'):

        """
        General bookiping of the luminosity function class
        
        """

        self.redshift=redshift
        self.whichlf=whichlf
        
        if('Buow15' in self.whichlf):
            print ('Init Buow15 LF at z={}'.format(self.redshift))
            self.param=self.Bouwens15()
        elif('Cassata11' in self.whichlf):
            print ('Init Cassata11 LF at z={}'.format(self.redshift))
            self.param=self.Cassata11()
        elif('Grove09' in self.whichlf):
            print ('Init Grove09 LF at z={}'.format(self.redshift))
            self.param=self.Grove09()
        elif('Reddy09' in self.whichlf):
            self.param=self.Reddy09()
        else:
            raise TypeError('LF {} not known'.format(self.whichlf))
            
    def Bouwens15(self):

        """
        This evaluates the luminosity function from Bouwens et al. 2015, ApJ, 803, 34
        using the fitting formulae of Sec. 5.1. Genrally valid between z~3-8
        
        """


        alpha=-1.87-0.1*(self.redshift-6)
        phistar=0.47*1e-3*10**(-0.27*(self.redshift-6)) #Mpc^-3
        Mstar=-20.95+0.01*(self.redshift-6) #AB abs mag 
        ltype='Schechter'
        
        return {'alpha':alpha,'type':ltype,'phi':phistar,'Mstar':Mstar}
    
    def Reddy09(self):

        """
        This evaluates the luminosity function from Reddy & Steidel 2009 ApJ 692, 778
        Use table 3 paramaters with switch on redshift
        
        """
        
        if(self.redshift < 2.7):
            #this is valid between z = 1.9 - 2.7
            alpha=-1.73
            phistar=2.75e-3
            Mstar=-20.70 #AB abs mag 
        else:
            #this is valid between z = 2.7 - 3.4
            alpha=-1.73
            phistar=1.71e-3
            Mstar=-20.97 #AB abs mag 

        ltype='Schechter'        
        return {'alpha':alpha,'type':ltype,'phi':phistar,'Mstar':Mstar}
    


    def Cassata11(self):

        """
        This evaluates the Lya luminosity function from Cassata et al. 2011, A&A 525, A143
        Generlly valid between z~2-6. Use corrected IGM values 
        
        """
        
        if(self.redshift < 3.0):
            #use redshift 2-3
            alpha=-1.6
            phistar=7.1e-4
            Mstar=42.70
            ltype='SchLum'
        elif((self.redshift >= 3.0) & (self.redshift < 4.55)):
            #use redshift 3-4.5
            alpha=-1.78
            phistar=4.8e-4
            Mstar= 42.70 
            ltype='SchLum'
        elif((self.redshift >= 4.55) & (self.redshift < 6.6)):
            #use redshift 4.55-6.6
            alpha=-1.69
            phistar=9.2e-4
            Mstar=42.72
            ltype='SchLum'


        return {'alpha':alpha,'type':ltype,'phi':phistar,'Mstar':Mstar}

    def Grove09(self):

        """
        This evaluates the Lya luminosity function from Grove et al. 2009, AA 497, 689-702
        Generally valid around z~3 
        
        """
        alpha=-1.74
        phistar=3.2e-4
        Mstar=43.3
        ltype='SchLum'

        return {'alpha':alpha,'type':ltype,'phi':phistar,'Mstar':Mstar}

    def plot(self):

        """
        Plot the luminosity function 
        """

        import matplotlib.pyplot as plt


        #make the     
        delta=0.1
        if('Schechter' in self.param['type']):
            Marr=np.arange(-24,-15,delta)
        elif('SchLum' in self.param['type']):
            Marr=np.arange(40,44,delta)
                
        lf=self.eval(Marr)

        ax=plt.subplot(111)
        ax.plot(Marr,lf)
        ax.set_yscale("log")
        plt.show()
        
        
    def eval(self,Marr):

        """
        Evaluate the luminosity function at a given Mag/lum value

        """

        lf=None

        if('Schechter' in self.param['type']):
            #in magnitude
            lf=self.param['phi']*(np.log(10.)/2.5)*\
                10**(-0.4*(Marr-self.param['Mstar'])*(self.param['alpha']+1))\
                *np.exp(-10**(-0.4*(Marr-self.param['Mstar'])))
            
        elif('SchLum' in self.param['type']):
            #in luminosity 
            lf=self.param['phi']*np.log(10.)*\
                10**((Marr-self.param['Mstar'])*(self.param['alpha']+1))\
                *np.exp(-10**(Marr-self.param['Mstar']))
            
        else:
            raise TypeError('LF {} not known'.format(self.param['type']))

        return lf

    def integrate(self,bright,faint):

        """
        Compute the integral of the luminosity function between two values
        
        """
        
        import scipy.integrate as integrate
        
        result = integrate.quad(lambda x: self.eval(x),bright,faint)
        
        return result

