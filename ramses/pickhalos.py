"""
Select a good halo for zoom in 

mass = desired mass of the halo (Msun)
halocat = rockstar halo cat 
snap = ramses snapshop
output = root of output names 
deltamass = tolerance of mass range (fraction) 

"""

def pickhalos(mass=1e12,halocat='halos_0.0.ascii',snap='info_00001.txt',
              deltamass=0.05,output='pickhalo'):
    
    import yt 
    import mypython as mp 
    import numpy as np 
    
    #load the simulation
    sim = yt.load(snap)

    #load rockstar ascii 
    halos=mp.ramses.parserockstar.parserockstar(halocat)
    nhalos=len(halos['id'])
    
    #select halos in correct mass range 

    #masses in Msun 
    masses=np.array(halos['mvir'])/sim.hubble_constant

    #filter 
    mmin=mass-mass*deltamass
    mmax=mass+mass*deltamass
    
    select=np.nonzero((masses >= mmin) & (masses <= mmax))
    select=select[0]
    nselect=select.size
    
    print "Selected ", nselect, " halos with masses [ ", mmin, ",", mmax, " ]"

    #to slice 
    #print halos['mvir'][select]
      
    #open output 
    outfile=open(output+'_haloinfo.txt','w')
    outfile.write('#HaloID,Nnearest,MaxM,NearestID\n')

    #now do some stats on the halos 
    for hhh in range(nselect):

        outfile.write(str(halos['id'][select[hhh]])+'\t')

        #find halos within 4Rvir 
        #distance to halo under exam. in comoving Mpc/h
        radius=(halos['x']-halos['x'][select[hhh]])**2+(halos['y']-halos['y'][select[hhh]])**2+\
            (halos['z']-halos['z'][select[hhh]])**2
        radius=np.sqrt(radius)
        
        #3Rvir in comoving Mpc/h
        fourrvir=4*halos['rvir'][select[hhh]]*1e-3
        
        #find nearest 
        nearest=np.nonzero((radius <= fourrvir) & (radius > 0))
        nearest=nearest[0]
        outfile.write(str(nearest.size)+'\t')
        
        #most massive nearest 
        massnear=np.array(halos['mvir'][nearest])/sim.hubble_constant 
        if nearest.size > 0:
            outfile.write(str(np.max(massnear))+'\t')
        else:
            outfile.write('0.0\t')
            
        #write id nearest 
        outfile.write(','.join([str(x) for x in halos['id'][nearest]])+'\n')  
        
    #close     
    outfile.close()

    #now onto display options 
    plotselecthalos('x',sim,halos,select,output)
    plotselecthalos('y',sim,halos,select,output)
    plotselecthalos('z',sim,halos,select,output)
    

#formatter functions to map grid units to box size
resol=0
boxsize=0
def mjrFormatter(x, pos):
    global resol
    global boxsize
    #The two args are the value and tick position
    off=x/resol*boxsize #phys 
    return "%.1f" % off


def plotselecthalos(axis,sim,halos,select,output):

    """
    Function that plots a projection of the cube with selected halos    

    """

    import matplotlib as mpl
    import yt 
    from yt.visualization.fixed_resolution import FixedResolutionBuffer
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    import numpy as np
    
    
    #define axis and permutation
    ax={'x':0,'y':1,'z':2}
    ax1={'x':'y','y':'z','z':'x'}
    ax2={'x':'z','y':'x','z':'y'}
    
    #make a projection by hand along one axis 
    global resol
    resol=512  
    global boxsize
    boxsize=(sim.domain_width.in_units('Mpccm/h'))[0].v #comoving Mpc/h

    projx=sim.proj(("deposit","all_density"),ax[axis])
    image2d = FixedResolutionBuffer(projx, (0.0, 0.999, 0.0, 0.999),(resol,resol))

    #some display stuff
    fig=plt.figure(1,figsize=(8,8), dpi=300)
    imgplot = plt.imshow(np.log10(image2d["deposit","all_density"].v),origin='lower')
    imgplot.set_clim(-3,-1)

    #cbar=plt.colorbar()
    #cbar.set_label('log Density'+' (g cm$^{-2}$)', rotation=90)

    #loop over the halos and plot them 
    print "Looping over selected ", select.size
    for hh in range(select.size):

        radius=10*halos['rvir'][select[hh]]*1e-3 #cMpc/h

        circle = plt.Circle((halos[ax1[axis]][select[hh]]*resol/boxsize,
                             halos[ax2[axis]][select[hh]]*resol/boxsize), 
                            radius=radius,fc='None', ec='black', linewidth=1)
        plt.gca().add_patch(circle)
        plt.gca().text(halos[ax1[axis]][select[hh]]*resol/boxsize-10,
                       halos[ax2[axis]][select[hh]]*resol/boxsize-10,
                       str(halos['id'][select[hh]]),fontsize=8,color='black')

    #map the axis to box size 
    plt.xlabel(ax1[axis]+' (cMpc h$^{-1}$)')
    plt.ylabel(ax2[axis]+' (cMpc h$^{-1}$)')

    axe1 = plt.gca()
    axe1.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjrFormatter))
    axe1.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjrFormatter))

    plt.savefig(output+"_"+axis+"_proj.ps")
    fig.clf()










    


    



    

    


