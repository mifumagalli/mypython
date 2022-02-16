
"""

Takes input file written as below and run systems through mcmcm

# name  zabs  ezabs ion logN elogN Flag Sample
J012156.03+144823.8_z2.662 2.66231 0.0001 HI 19.25 0.2 0 P15
J012156.03+144823.8_z2.662 2.66231 0.0001 CI 12.412 0.05 -1 P15
J012156.03+144823.8_z2.662 2.66231 0.0001 CII 14.794 0.05 -2 P15
J012156.03+144823.8_z2.662 2.66231 0.0001 CIV 14.505 0.05 -2 P15
J012156.03+144823.8_z2.662 2.66231 0.0001 OI 14.671 0.05 -2 P15
J012156.03+144823.8_z2.662 2.66231 0.0001 AlII 13.016 0.05 0 P15
J012156.03+144823.8_z2.662 2.66231 0.0001 SiII 14.007 0.05 0 P15
J012156.03+144823.8_z2.662 2.66231 0.0001 SiIV 14.094 0.05 0 P15
J012156.03+144823.8_z2.662 2.66231 0.0001 FeII 13.435 0.05 0 P15
J012156.03+144823.8_z2.662 2.66231 0.0001 NiII 12.877 0.05 -1 P15
J015741.56-010629.6_z2.631 2.63129 0.0001 HI 19.45 0.2 0 P15
J015741.56-010629.6_z2.631 2.63129 0.0001 CI 12.49 0.05 -1 P15
J015741.56-010629.6_z2.631 2.63129 0.0001 CIV 14.658 0.05 0 P15
J015741.56-010629.6_z2.631 2.63129 0.0001 AlII 12.69 0.05 0 P15
J015741.56-010629.6_z2.631 2.63129 0.0001 FeII 13.649 0.05 0 P15
J015741.56-010629.6_z3.385 3.38541 0.0001 HI 18.35 0.75 -3 P15
J015741.56-010629.6_z3.385 3.38541 0.0001 CII 14.597 0.061 -2 P15
J015741.56-010629.6_z3.385 3.38541 0.0001 OI 13.905 0.05 0 P15
J015741.56-010629.6_z3.385 3.38541 0.0001 SiII 13.789 0.05 0 P15
J015741.56-010629.6_z3.385 3.38541 0.0001 SiIV 14.2 0.05 0 P15
J015741.56-010629.6_z3.385 3.38541 0.0001 NiII 13.349 0.05 -1 P15


This is the code that loop over data and runs mcmc chains.
For information/acknowledgment see https://ui.adsabs.harvard.edu/abs/2016MNRAS.455.4100F/abstract 


Written by Michele Fumagalli in Durham, Summer 2015
michele.fumagalli@durham.ac.uk

Copyright (C) 2015 Michele Fumagalli

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


""" 


 
def run_mcmc(argv):


    import argparse
    import mcmc_hd5 as mcmc_ions
    import numpy as np
    from astropy.table import Table
    import os

    #get the call
    parser = argparse.ArgumentParser(description='Running grid on shared memory system')
    parser.add_argument('datain',help='data file containing ions')
    parser.add_argument('grid',help='file containing models')
    parser.add_argument('--outdir',default='output',help='location to store output')
    parser.add_argument('--proc', type=int,default=12,help='number of processors')
    parser.add_argument('--restartdir',help='Folder with output to be used for generating new guess for chains')
    args=parser.parse_args(argv)
 
 
    #read in the data file
    print('Loading the data...')
    data=Table.read(args.datain,format='ascii')

    #find unique filenames
    allsystems=set(data['name'])


    #loop on systems to fit in this batch 
    for sy in allsystems:
     
        #find the ions for the current system
        sets=np.where(data['name'] == sy.strip())[0]
        
        #now loop over ions for this system and initialise the tuples
        observ=[]
        for ii in range(len(sets)):
            if(data['ion'][sets[ii]] == 'HI'): 
                               
                #isolate hydrogen and info
                #print ii, sets[ii], data['ion'][sets[ii]], sy
                obsinfo={'NHI':data['logN'][sets[ii]],'eNHI':data['elogN'][sets[ii]],'hiflag':data['Flag'][sets[ii]],
                         'z':data['zabs'][sets[ii]],'errz':data['ezabs'][sets[ii]],'name':sy.strip()}
            else:
                #append ion value to list 
                observ.append((data['ion'][sets[ii]],data['logN'][sets[ii]],data['elogN'][sets[ii]],data['Flag'][sets[ii]]))
                
                
        #run the mcmc
        print('Ready to run mcmc for {}'.format(sy))
        
        #check if read to restart
        restartchain=None
        if(args.restartdir):
            #check if output exists for this system
            rest_name=args.restartdir+'/'+sy+'_emcee.hd5'
            if(os.path.isfile(rest_name)):
                restartchain=rest_name
            
        #pick optimised values for 12 processors - cosma (proc*NN walkers, proc*YY samples)
        mcmc=mcmc_ions.mcmc_ions(observ,obsinfo,args.grid,nwalkers=(args.proc*80),nsamp=(args.proc*40),
                                 optim=False,effnhi=True,threads=args.proc,outsave=args.outdir,restartchn=restartchain)

    print('All done with this batch')
    
    return

if __name__ == "__main__":
    import sys
    run_mcmc(sys.argv[1:])
    
  











