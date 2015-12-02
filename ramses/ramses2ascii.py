def ramses2ascii(snap,path='./',output="ascii/"):
   
    """
    Read an unformatted binary file for ramses particles and 
    convert it into an ascii file with X Y Z VX VY VZ ID
    Use yt for this
    
    snap is the snapshot number
    output is where you want to write output 


    Use units for rockstar:
    position in comoving Mpc/h
    velocity in physical units km/s
    mass in Msun/h
    box size in comoving Mpc/h
    
    """
    
    import yt 
    from astropy.io import ascii
    from astropy.table import Table
    import os 

    #compose snapname 
    snapname=path+'output_'+format(snap, '05')+'/info_'+format(snap, '05')+'.txt'
    
    #open sim 
    sim = yt.load(snapname)
    
    #create folder 
    if not os.path.exists(output):
        os.mkdir(output) 
    
    #do work only if files is not there
    fname=output+'asciipart_'+format(snap, '05')+'.txt'
    if os.path.isfile(fname):
        print "ascii files exist..."
    
    else:
        
        print "Cearing table.."
        ad=sim.all_data()
        
        #conversion in physical cgs units
        #print sim.length_unit
        #print sim.mass_unit
        #print sim.velocity_unit
        #print sim.current_redshift
        #print sim.hubble_constant
        #print "Box size (cMpc/h) ", sim.domain_width.in_units('Mpccm/h')
        #print "Redshift ", sim.current_redshift 
        #print "Mass particle (Msun) ", sim.mass_unit.v/2e33*ad['io', 'particle_mass'][0].v

        #box size in comoving Mpc/h
        boxsize=(sim.domain_width.in_units('Mpccm/h'))[0].v
        
        #conversion to velocity in km/s
        velocity=sim.velocity_unit/1e5
        
        #create datatable
        data=Table([ad['io', 'particle_position_x'][:].v*boxsize,ad['io', 'particle_position_y'][:].v*boxsize,
                    ad['io', 'particle_position_z'][:].v*boxsize,ad['io', 'particle_velocity_x'][:].v*velocity,
                    ad['io', 'particle_velocity_y'][:].v*velocity,ad['io', 'particle_velocity_z'][:].v*velocity,
                    ad['io', 'particle_identifier'][:].v],names=['PX','PY','PZ','VX','VY','VZ','ID'])
        
        #output 
        print "Writing table to "+fname
        ascii.write(data,fname)
        
