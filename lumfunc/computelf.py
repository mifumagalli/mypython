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


    def plot(self):

        """
        Plot the luminosity function 
        """

        import matplotlib.pyplot as plt


        #make the mag array 
        Marr=np.arange(-24,-15,0.1)

        lf=self.eval(Marr)
        
        plt.plot(Marr,np.log10(lf))
        plt.show()
        
        
    def eval(self,Marr):

        """
        Evaluate the luminosity function at a given Mag value

        """

        lf=None

        if('Schechter' in self.param['type']):
            
            lf=self.param['phi']*(np.log(10.)/2.5)*\
            10**(-0.4*(Marr-self.param['Mstar'])*(self.param['alpha']+1))\
            *np.exp(-10**(-0.4*(Marr-self.param['Mstar'])))
            
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

