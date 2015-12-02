def halofind(path="./",output="ascii/",rockfile='rockstar.cfg',rockout="rockstar/"):

    """ 
    Perfom a number of initialization steps to run rockstar halo finder on 
    a ramses simulation. Specifically, conver the dark matter particles from
    binary files to ascii lists, and generate a script that can be used to 
    run rockstar on each snapshot

    Force resolution is computed assuming unigird. 
    Keep track of units throughout.

    """
    
    import yt 
    import mypython as mp
    import os 


    #create rock output folder
    if not os.path.exists(rockout):     
        os.mkdir(rockout)

    #find how many snapshots are there and process them
    for out in os.listdir(path):
        #check if valid output 
        if "output_" in out:
            #grap snap number 
            snap=int(out.split("output_")[1])
            print "Processing snap ", snap
            #process in ascii 
            mp.ramses.ramses2ascii.ramses2ascii(snap,path=path,output=output)

            #create rockstar output 
            rockfilename='output_'+format(snap, '05')+rockfile
            if os.path.isfile(rockfilename):
                print "rockstar config file exist..."      
            else:   
                
                #get sim properties
                snapname=path+'output_'+format(snap, '05')+'/info_'+format(snap, '05')+'.txt'
                sim = yt.load(snapname)  
                ad=sim.all_data()

                #write config file 
                fle = open(rockout+rockfilename, 'w')
                fle.write('#Configuration file for rockstar\n') 
                fle.write('FILE_FORMAT = "ASCII"\n')
                
                #mass particle Msun/h
                mpart=sim.mass_unit.v/2e33*ad['io', 'particle_mass'][0].v*sim.hubble_constant
                fle.write('PARTICLE_MASS = '+str(mpart)+'\n')
                
                #cosmology 
                aexp=1./(1.+sim.current_redshift)
                fle.write('SCALE_NOW = '+str(aexp)+'\n')
                fle.write('h0 = '+str(sim.hubble_constant)+'\n')
                fle.write('Ol = '+str(sim.omega_lambda)+'\n') 
                fle.write('Om = '+str(sim.omega_matter)+'\n')

                #box size 
                boxsize=(sim.domain_width.in_units('Mpccm/h'))[0].v
                fle.write('BOX_SIZE = '+str(boxsize)+'\n')     
                fle.write('FORCE_RES = '+str(boxsize/sim.domain_dimensions[0]/2.)+'\n')   
                
                #utility 
                fle.write('OUTBASE = '+path+'output_'+format(snap, '05')+'/\n')
                fle.write('MIN_HALO_OUTPUT_SIZE = 100\n') 
    
                fle.close()
                exit()
                                

                
