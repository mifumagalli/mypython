
def flux_integrator(wave,flux,low,high,error=None,profile='box',par=None):
    
    """
    Performs spectral integration of a spectrum 
    
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
    
    #grab length of spectrum
    nwav=len(wave)
    
    #if error not defined, zero it
    if not (error):
        error=np.zeros(nwav)
            
    #check and work on profile 
    if('box' in profile):
        weight=np.zeros(nwav)+1.
    elif('gauss' in profile):

        #sanity check
        if(len(par) < 2):
            print ('For gauss profile, I am expecting centre and width!')
            return
        else:
            print ('Setting up gaussian with mu={} and sigma={}'.format(par[0],par[1]))
        aa=0.
    else:
        print ('Profile {} is not recognized!'.format(profile))
        return
        

