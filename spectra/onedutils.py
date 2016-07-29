
"""

Collection of 1D untilities for spectra analysis 

"""

def flux_extraction(wave,flux,low,high,error=None,profile='box',par=None):
    
    """
    Performs spectral integration of a spectrum with optimal avarage  
    
    wave  -> wavelenght in any unit WW
    flux  -> flux in any unit, but assumed to be x1/WW
    low   -> lower end in wave of integral
    high  -> upper end in wave of integral 
    error -> 1 sigma error in flux unit
    profile -> profile to use in integration:
    box   -> simple box extraction
    gauss -> gaussian profile. Expect center and sigma as par[0] and par[1].
    
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.integrate as integrate
    import scipy.interpolate as interpol

    #name sure is all in numpy array
    wave=np.array(wave)
    flux=np.array(flux)
    
    #grab length of spectrum
    nwav=len(wave)
    
    #if error not defined, zero it
    if not (error):
        error=np.zeros(nwav)
    else:
        error=np.array(error)

    #check and work on profile 
    if('box' in profile):
        #set weight to unity
        weight=np.zeros(nwav)+1.
    elif('gauss' in profile):

        #sanity check
        if not(par):
            print ('For gauss profile, I am expecting centre and width!')
            return
        else:
            print ('Setting up gaussian with mu={} and sigma={}'.format(par[0],par[1]))
            weight=1./(np.sqrt(2*np.pi)*par[1])*np.exp(-1.*((wave-par[0])**2/(2.*par[1]*par[1])))
    else:
        print ('Profile {} is not recognized!'.format(profile))
        return

    #Check gaussian integrate to 1
    #farg=interpol.interp1d(wave,weight)
    #result = integrate.quad(farg,5000,7000)
    #print result

    #check shape of weight
    #plt.plot(wave,weight)
    #plt.show()

    print ('Still in the works..')
    exit()


def flux_integrator(wave,flux,low,high,error=None,silent=False):
    
    """
    Performs spectral integration of a spectrum with optimal avarage  
    
    wave  -> wavelenght in any unit WW
    flux  -> flux in any unit, but assumed to be x1/WW
    low   -> lower end in wave of integral
    high  -> upper end in wave of integral 
    error -> 1 sigma error in flux unit
    silent -> turn off screen output

    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.integrate as integrate
    import scipy.interpolate as interpol

    #name sure is all in numpy array
    wave=np.array(wave)
    flux=np.array(flux)
    
    #grab length of spectrum
    nwav=len(wave)
    
    #if error not defined, zero it
    if (len(error) > 1):
        error=np.array(error)
    else:
        error=np.zeros(nwav)
   
    #set up the integrator function 
    flxint=interpol.interp1d(wave,flux)
    dlambda=wave-np.roll(wave,1)
    errint=interpol.interp1d(wave,error*error*dlambda)

    #integrate the flux
    result_flux,error = integrate.quad(flxint,low,high)
    #print result
    if not (silent):
        print('Integrating flux between {} and {} gives {}'.format(low,high,result_flux))
    
    #deal with error
    result_error,error = integrate.quad(errint,low,high)
    result_error=np.sqrt(result_error)
    #print result
    if not (silent):
        print('Propagating error between {} and {} gives {}'.format(low,high,result_error))
    
    return [result_flux,result_error]

    #Check gaussian integrate to 1
    #farg=interpol.interp1d(wave,weight)
    #result = integrate.quad(farg,5000,7000)
    #print result

    #check shape of weight
    #plt.plot(wave,weight)
    #plt.show()

 


    

