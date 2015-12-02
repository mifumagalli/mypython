"""
This program eneable the selection of a zoom-in region 
from a ramses output 

haloid - id of the halo to consider 
nrvir  - zoom-in region in rvir units 
halocat - rockstar halo catalogue
snap    - simulaiton snapshots
icpath - location of initial conditions 
output - rootname for output 
asciifile - ascii particle list 
ascii_ic  - ascii particle list for initial conditions
snap_ic   - simulation snapshot for initial conditions 

"""

def zoominic(haloid,nrvir,halocat='',snap='',icpath='',output='',asciifile='',
             ascii_ic='',snap_ic=''):
    
    import yt
    import mypython as mp 
    import numpy as np
    from astropy.io import ascii
    import cPickle
    import os 
    
    #load sim/halos
    sim=yt.load(snap)

    #load rockstar ascii 
    halos=mp.ramses.parserockstar.parserockstar(halocat)
    nhalos=len(halos['id'])
    
    #find right halo 
    if(halos['id'][haloid] != haloid):
        print "Incorrect halo!"
        exit()
        
    print "Selected halo ID ", halos['id'][haloid], "Check (", haloid, ")"

    #selection distance in rvir (cMpc/h) 
    distance=(nrvir*halos['rvir'][haloid]*1e-3)**2
    iddump=icpath+output+'_idpart.dat'

    #check if particles id already exist 
    if os.path.isfile(iddump):
        print "Particle id file exists..."      

        #load back 
        outid=open(iddump,'rb')
        idgood=cPickle.load(outid)
        outid.close()

    else:   
        
    #find particle ids within region of interest 
    #ascii file is in units of cMpc/h and km/s
    #PX PY PZ VX VY VZ ID
    
        parfile = open(asciifile,'r')
        comment=parfile.readline()
        idgood=[]
        print "Vetting particles... "
        for line in parfile:
            partline=(line.rstrip('\n')).split(" ")
            #if partline == '': break
            dist=(float(partline[0])-halos['x'][haloid])**2+\
                (float(partline[1])-halos['y'][haloid])**2+\
                (float(partline[2])-halos['z'][haloid])**2
            if(dist <= distance):
                idgood.append(long(float(partline[6])))
        parfile.close()
        
        #write output to list 
        outid=open(iddump,'w')
        cPickle.dump(idgood,outid)
        outid.close()
        
    
    #now match to initial conditions 
    #use first output in ramses simulation which replicates the 
    #music ic in more digestible format 
    sim_ic=yt.load(snap_ic)
    ellipsregion=icpath+output+'_ellipse.txt'
    #box size in comoving Mpc/h
    boxsize_ic=(sim_ic.domain_width.in_units('Mpccm/h'))[0].v
        
    #check if ellipsoid region  already exist 
    if os.path.isfile(ellipsregion):
        print "Ellipsoid definition file already exists..."      
    else:   

        #open for wirite
        out = open(ellipsregion,'w')
        #open for read
        read = open(ascii_ic,'r')
        comment=read.readline()
        npartfound=0

        print "Vetting IC particles... "
        for line in read:
            partline=(line.rstrip('\n')).split(" ")
            idpart=long(float(partline[6]))
            
            #if in list keep 
            if idpart in idgood:
                #back to box units 
                out.write(str(float(partline[0])/boxsize_ic)+' ')
                out.write(str(float(partline[1])/boxsize_ic)+' ')
                out.write(str(float(partline[2])/boxsize_ic)+'\n')
                npartfound=npartfound+1
        #close
        out.close()
        read.close()
        
        #parity check 
        print "Written ", npartfound, "particles out of ", len(idgood)
    



    
    

    

