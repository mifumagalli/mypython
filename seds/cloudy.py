 
def init_hiireg(param,path='./',pressure=False,extinction=False):
    
    
    """ 
    This is a function that initialise a cloudy model of an HII region
    
    #Units for input paramaters [param] 
    TIME Myr
    QHO  log N/s 
    DENS log cm^-3
    ZSTAR Mass fraction
    ZGAS Mass fraction
    LGFREQ Log of frequency in Hz
    LGFLUX Log of f_nu
    MASS The mass of the stellar pop in the burst - not used directly by cloudy

    The input cloudy file is written to disk at "path" location
    
    """
    
    import numpy as np
    from scipy.constants import c
    from scipy.constants import k as kB
    
    #precompute some relevant quantities
    
    #In the following, treat the Z values as absolute mass fraction
    #Since we are using the Asplung/Grevesse 2010 metallicities, solar 
    #abundance is set to be Z = 0.0134
    #The scaling for gas-phase metallicity thus becomes
    Zsolar=0.0134
    Zgasscale=param["ZGAS"]/Zsolar
   
    if extinction == True:
        
        #The ISM grains are defined as dust-gas ratio in the MW. 
        #In this case, scale arbitrarily
        Grscale=10**param["EXT"]
        
    else:
        
        #The ISM grains are defined as dust-gas ratio in the MW. To keep the same
        #dust to metal ratio, scale consistently
        Grscale=Zgasscale

    #Compute the inner radius and density self-consistently 
    #following the HII region prescription of Kru&Maz 2009
    # Get characteristic radius and time; all quantities
    # defined as in Krumholz & Matnzer (2009)
    age  = param["TIME"]*1e6 #Time in yr
    qH0  = 10**param["QHO"]
    hden = 10**param["DENS"]
        
    alphaB = 2.59e-13       # Case B recombination coefficient
    mu = 2.34e-24           # Mean mass per H nucleus, in g
    eps0 = 2.179e-11        # H ionization potential, in erg
    TII = 1e4               # HII region temp, in K
    psi = 3.2               # Mean photon energy / eps0
    ft = 2.0                # Trapping factor
    phi = 0.73              # Dust absorption fraction
        
    rch = alphaB/(12.0*np.pi*phi) * (eps0/(2.2*kB*1e7*TII))**2 * ft**2 * psi**2 * qH0 / (c*1e2)**2
    tch = (4.0*np.pi * mu*hden*c * rch**4 / (3.0*ft*qH0*psi*eps0))**0.5
        
    # Get xIIgas, xIIrad
    tau = age*365.25*24.*3600./tch
    xIIrad = (2.0*tau**2)**0.25
    xIIgas = (49.0*tau**2/36.0)**(2.0/7.0)
    
    # Get outer radius, inner radius, density
    r = rch*(xIIrad**3.5 + xIIgas**3.5)**(2.0/7.0)
    r0 = r/1e2
    nH = (3.0*qH0 / (4.0*np.pi*alphaB*r**3))**0.5
        
    #assume a stopping radius/outer radius of 20.188 log[cm]
    rstop=20.188
    
    #contruct sed names with self consistent metallicity 
    strz=str(param["ZSTAR"])
    
    #create the model name
    root=param["SEDS"].split("/")
    
    if extinction == True:
        modname_basic='{0}_t{1:.1f}_q{2:.1f}_n{3:.1f}_zs{4:.4f}_zg{5:.4f}_ex{6:.2f}'\
            .format(root[-1],param["TIME"],param["QHO"],param["DENS"],param["ZSTAR"],param["ZGAS"],param["EXT"])
    else:
        modname_basic='{0}_t{1:.1f}_q{2:.1f}_n{3:.1f}_zs{4:.4f}_zg{5:.4f}'\
            .format(root[-1],param["TIME"],param["QHO"],param["DENS"],param["ZSTAR"],param["ZGAS"])
        

    #open for write the parameter files 
    modname=path+modname_basic+".in"
    flin=open(modname,"w")

    #Store info on the model that are not used directly in cloudy
    flin.write("#MassPop {0}\n".format(param["MASS"]))
    flin.write("#RStrom {0:f}\n".format(np.log10(r)))
    if extinction == True:
        flin.write("#Ext {0:f}\n".format(param["EXT"]))

    #set abundances
    flin.write("#This sets the metallicity and grains\n")
    flin.write("abundances GASS10 no grains\n")
    flin.write("metals {0:f}\n".format(Zgasscale))
    flin.write("grains ISM {0:f}\n".format(Grscale))
    #turn off quantuum heating for stability - it affects only thermal continuum 
    flin.write("no grains qheat\n") 

    #set the geometry and density 
    flin.write("#This sets the geometry/density\n")
    flin.write("hden {0:f}\n".format(np.log10(nH)))
    flin.write("sphere static\n")
    flin.write("radius {0:f}\n".format(np.log10(r0))) # log in cm
    if (pressure == True):
        flin.write("constant pressure no abort\n") 

    #set the radiation field
    flin.write("#This sets the radiation\n")
    flin.write("cmb\n")
    flin.write("cosmic rays background\n")

    #This is to use a table
    #flin.write('table star "{0}" age={1}e6 years\n'.format(sedname,param["TIME"]))

    #This is to call the function that write interpolation commands for cloudy 
    interpolate_sb99(param["LGFREQ"],param["LGFLUX"],flin)
    flin.write("Q(H) = {0:f}\n".format(param["QHO"]))
    
    #set the stop/save
    flin.write("#This sets the stop/save\n")
    flin.write("stop temperature 1000\n")
    flin.write("stop efract -2\n")
    flin.write("stop radius {0:f}\n".format(rstop))
    flin.write("iterate to convergence\n")
    flin.write('set save prefix "{0}"\n'.format(modname_basic))

    #Minimum output 
    #NB on units: nu * L(nu)/ (4pi ro^2) [ erg cm-2 s-1 ] relative to the inner radius of the cloud 
    flin.write('save last continuum units Angstrom ".cont"\n')
    flin.write('save last line list emergent absolute column ".lines" "LineList_HII.dat"\n')
    flin.write('save last line list absolute column ".lint" "LineList_HII.dat"\n')
    flin.write('save last grain extinction ".ext"\n')

    #Other optional output 
    #flin.write('save last physical conditions ".phys"\n')
    #flin.write('save last overview ".over"\n')

    flin.close()

    return

def compile_sb99(path='./'):

    """ This is a simple script that compiles sb99 seds """

    import os
    import glob
    import subprocess

    #store current wd 
    cwd=os.getcwd()
    #go to where models are
    os.chdir(path)

    #query for sb99 files 
    sb99list=glob.glob("*.stb99")

    #compile each of them
    for f in sb99list:
        
        #prepare input file 
        clfl=open("cloudy.in","w")
        clfl.write('compile star "'+f+'" stellar')
        clfl.close()

        #run cloudy 
        cloudydir=os.environ['CLOUDY_DIR']
        subprocess.call([cloudydir+"/cloudy.exe","-r","cloudy"])
    
    #get back to current wd
    os.chdir(cwd)
    
    return

def interpolate_sb99(lognu,logfnu,fil):

    """ 
        This is a simple script that interpolates sb99 seds and creates 
        cloudy input commands. Take care of the range that cloudy wants


        lognu  -> the log of freqeuncy in Hz
        logfnu -> the log of flux per unit frequency 
        fil    -> a pointer to a file open for writing 

    """

    import numpy as np

    #make copies     
    spec = 10**np.array(logfnu)
    specclean = np.copy(spec)
    logfreq=np.copy(lognu)
    
    #first point
    specclean[spec == 0.0] = np.amin(spec[spec > 0])*1e-6
    logL_nu = np.log10(specclean)

    #put a floor at -40 when writing to avoid cloudy complaints 

    fil.write("interpolate")
    fil.write(" ({0:f} {1:f})".format(7.51, 
                                        max(np.amin(logL_nu)-6,-40)))
    fil.write(" ({0:f} {1:f})".format(logfreq[-1]-0.01, 
                                        max(np.amin(logL_nu)-6,-40)))
    
    #stuff in between 
    for i in range(len(logL_nu)):
        #if i % 4 == 0:
        fil.write("\ncontinue")
        fil.write(" ({0:f} {1:f})".format(logfreq[-i-1],
                                          max(logL_nu[-i-1],-40)))
    #last point 
    fil.write("\ncontinue ({0:f} {1:f})".
                format(logfreq[0]+0.01, max(np.amin(logL_nu)-6,-40)))

    fil.write(" ({0:f} {1:f})\n".format(22.4,max(np.amin(logL_nu)-6,-40)))

    return

    
def run_grid(nproc=-1,noclobber=1):
    
    """
    
        This script queries the local folder for .in files and run a 
        all the jobs on the available number of processors 
    
        nproc:       useful to set the number of processors to be used
        noclobber:   grids are not overwritten if .out if found 
        
    """
    
    import multiprocessing
    from multiprocessing import Process
    import subprocess
    import threading
    import warnings
    import os 
    import sys
    import glob
    import numpy as np
    
    #get the number of available processors without overshooting
    if (nproc == -1):
        nproc = multiprocessing.cpu_count()
    elif (nproc > multiprocessing.cpu_count()):
        nproc = multiprocessing.cpu_count()
        
    print "Ready to start running on {0:d} processors".format(nproc)
    
    #find the actual input files 
    listin=glob.glob("*.in")
    listout=glob.glob("*.out")

    #find which models you actually want to run 
    if (noclobber == 1):
        listfree=[]
        for f in listin:
            isout=f.split(".in")[0]+".out"
            if isout not in listout:
                listfree.append(isout.split(".out")[0])
                
    else: 
        listfree=listin
        
    #now we now how many models are remaining 
    ntrials=len(listfree) 

    #check if actually need computations
    if (ntrials == 0):
        print "No need to run models!"
    else: 
        print "Running {0:d} models".format(ntrials)

        
    #define the dispaly out 
    def display_out(out, procnum):
        for line in iter(out.readline, ''):
            print("thread {:d}: ".format(procnum+1) + line.split('\n')[0])
            #out.close()
            
    completed_trials = 0
    proc_list = [None] * nproc
    io_threads = [None] * nproc
    err_threads = [None] * nproc
    ON_POSIX = 'posix' in sys.builtin_module_names
    

    #now loop over tasks continously till all have been performed
    while completed_trials < ntrials:
        
        # At each loop, check which processes have completed, so 
        # we now how many processors are idles and ready to accept
        # new jobs
        #
        # See if any running processes have completed; flag those
        # that have, and increment their counters
        proc_avail_list = []
        for i, proc in enumerate(proc_list):
            if proc is None:
                # Process not started, flag as available
                proc_avail_list.append(i)
            elif proc.poll() != None:
                # Process completed, flag as available
                proc_avail_list.append(i)
                
                
        #At this point, we know which processors are ready for more
        #Loop on those and assign a job
        for p in proc_avail_list:

            #While assigning, we may be running out of jobs.
            #If so, stop if all tasks have already been assigned
            if completed_trials >= ntrials:
                break
            
            #Now queue a new job: this is the part that is problem 
            #specific, and generally the only thing that should change 
            #if the code that is being called changes...
            print "Running model {0:d} of {1:d} on processor {2:d}".format(completed_trials+1,ntrials,p)
            cmd = os.environ['CLOUDY_DIR']+'cloudy.exe -r '+listfree[completed_trials]            

            #launch a new process 
            proc_list[p] \
                = subprocess.Popen(cmd, bufsize=0, shell=True,
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE,
                                   close_fds=ON_POSIX)
                
    
            # Increment number of trials assigned
            completed_trials = completed_trials + 1
            
            # Start threads to monitor stdout and stderr from this process
            io_threads[p] \
                = threading.Thread(target=display_out,
                                   args=(proc_list[p].stdout, p))
            io_threads[p].daemon = True
            io_threads[p].start()
            err_threads[p] \
                = threading.Thread(target=display_out,
                                   args=(proc_list[p].stderr, p))
            err_threads[p].daemon = True
            err_threads[p].start()


    # Step 6: wait for final processes to complete before exiting
    running = [True]*len(proc_list)
    while True in running:
        for i, proc in enumerate(proc_list):
            if running[i]:
                if proc is None:
                    # Process never started
                    running[i] = False
                elif proc.poll() != None:
                    running[i] = False         

    print "All done with jobs... exiting nicely!"


def parse(modename,extinction=False):
    
    """ 
    Parse the various outputs given a model name 

    Outmod is a list of dictionaris with info

    0 : warnings and paramaters
    1 : continuum
    2 : lines emergent 
    3 : lines intrinisc  
    4 : extinction 

    """

    import os
    import numpy as np
    
    #init structure of models
    outmod=[]

  
    ##EXT 0: check if out file exists, if so parse it
    if(os.path.isfile(modename+".out")): 
        
        #print 'Parse output for warnings'
        
        file1=open(modename+".out","r")
                
        for line in file1:
            if "Cloudy ends:" in line:
                runinfo = line
            if "[Stop " in line:
                warninfo = line
            if "Dust to gas ratio" in line:
                tmp=line.split("AV(ext):")
                tmp=tmp[1].split("(pnt):")
                avext=float(tmp[0])
                avpnt=float(tmp[1])
        file1.close()

        #now load some of the input paramaters 
        file1=open(modename+".in","r")

        #Init in case of no extinction 
        ext=False
        
        for line in file1:
            if "hden" in line:
                hdens = float((line.split(" "))[1])
            if "radius" in line and "stop" not in line:
                rad = float((line.split(" "))[1])
            if "#MassPop" in line:
                mass = float((line.split(" "))[1])
            if "#RStrom" in line:
                rstr = float((line.split(" "))[1])
            if "#Ext" in line:
                ext = float((line.split(" "))[1])
               
        file1.close()

        try:
            outmod.append({"RUN":runinfo.rstrip('\n'),"WARN":warninfo.rstrip('\n'),"HII":hdens,"R0":rad,"MASS":mass,"RStr":rstr,"Ext":ext,
                           "AV_EXT":avext,"AV_PNT":avpnt})
        except UnboundLocalError:
            print "Error with model... ", modename
            outmod.append(0)


    else: 
        outmod.append(0)

    ##EXT 1: check if continuum file exists, if so parse it
    if(os.path.isfile(modename+".cont")): 
        #print 'Parse continuum'
        
        file1=open(modename+".cont","r")
        
        wave=[]
        incid=[]
        transm=[]
        diffuse=[]
        nettransm=[]
        
        for line in file1:
            if "#" not in line:
                fields=line.split("\t")
                wave.append(float(fields[0]))
                incid.append(float(fields[1]))
                transm.append(float(fields[2]))
                diffuse.append(float(fields[3]))
                nettransm.append(float(fields[4]))

        file1.close()
                                            
        outmod.append({"wave":wave,"incident":incid,"transmitted":transm,"diffuse":diffuse,"net":nettransm})
        
    else: 
        outmod.append(0)


    ##check if line file exists, if so parse it
    if(os.path.isfile(modename+".lines")): 
        #print 'Parse lines'
        
             
        file1=open(modename+".lines","r")
        
        label=[]
        flux=[]
                
        for line in file1:
            fields=line.split("\t")
            
            if (len(fields) > 1):
                
                label.append(fields[0])
                flux.append(np.log10(float(fields[1])))
        
        file1.close()

        outmod.append({"label":label,"flux":flux})

    else: 
        outmod.append(0)
 
    ##check if line intrinisc file exists, if so parse it
    if(os.path.isfile(modename+".lint")): 
        #print 'Parse intrinisc lines'
                     
        file1=open(modename+".lint","r")
        
        label=[]
        flux=[]
        
        for line in file1:
            fields=line.split("\t")
            
            if (len(fields) > 1):
                
                label.append(fields[0])
                flux.append(np.log10(float(fields[1])))
        
        file1.close()

        outmod.append({"label":label,"flux":flux})

    else: 
        outmod.append(0)
 

    ##check if extinction file exists, if so parse it
    if(os.path.isfile(modename+".ext")): 
        #print 'Parse intrinisc lines'
                     
        file1=open(modename+".ext","r")
        
        depth=[]
        avpoint=[]
        avexten=[]

        for line in file1:
            fields=line.split("\t")
            
            if ('#' not in line):
                depth.append(np.log10(float(fields[0])))
                avexten.append(float(fields[1]))
                avpoint.append(float(fields[2]))

        file1.close()

        outmod.append({"depth":depth,"avexten":avexten,"avpoint":avpoint})

    else: 
        outmod.append(0)
 
    return outmod


def compress_grid(namefits,path=['./'],extinction=False): 
    

    """
    Compress a cloudy grid in a fits file with name namefits.
    If path is a list, dive in subdir.
    """

    import glob
    import mypython as mp
    from mypython.seds import cloudy
    from astropy.io import fits
    import numpy as np


    #init the fits file
    hdulist = fits.HDUList()
    nmod=0

    #loop over the paths 
    for dirr in path:
        
        #grb the models
        listfile = glob.glob(dirr+"*.in")
        
        #loop over the files 
        for fil in listfile:
            
            fields=(fil[:-3]).split("_t")
            tag=fields[0]
            
            fields=fields[1].split("_q")
            time = fields[0]
 
            fields=fields[1].split("_n")
            qh0 = fields[0]

            fields=fields[1].split("_zs")
            nH = fields[0]

            if extinction:
                fields=fields[1].split("_zg")
                zs = fields[0]

                fields=fields[1].split("_ex")
                zg = fields[0]
                ex = fields[1]

            else:
                
                fields=fields[1].split("_zg")
                zs = fields[0]
                zg = fields[1]
    

            mod=cloudy.parse(fil[:-3],extinction=extinction)
            

            #create columns with checks for line matches in emergent/intrins 
            linetags=np.array(mod[2]["label"])
            linetags2=np.array(mod[3]["label"])

            if not (linetags==linetags2).all():
                print "Something wrong with line list"
                exit()
            
            emerlum=np.array(mod[2]["flux"])
            intrlum=np.array(mod[3]["flux"])

            #write as table
            tbhdu = fits.BinTableHDU.from_columns(
                [fits.Column(name='LINES', format='20A', array=linetags),
                 fits.Column(name='EMERLUM', format='E', array=emerlum),
                 fits.Column(name='INTRLUM', format='E', array=intrlum)])
            
           
            #append to stack updating header
            tbhdu.header['SED']=tag
            tbhdu.header['TIME']=float(time)
            tbhdu.header['QH0']=float(qh0)
            tbhdu.header['NH']=float(nH)
            tbhdu.header['ZS']=float(zs)
            tbhdu.header['ZG']=float(zg)
            tbhdu.header['NII']=mod[0]["HII"]
            tbhdu.header['R0']=mod[0]["R0"]
            tbhdu.header['RStr']=mod[0]["RStr"]
            tbhdu.header['MASS']=mod[0]["MASS"]
            tbhdu.header['WARN']=mod[0]["WARN"]
            tbhdu.header['AVEXT']=mod[0]["AV_EXT"]
            tbhdu.header['AVPNT']=mod[0]["AV_PNT"]

            if extinction:
                tbhdu.header['EXT']=float(ex)
                
            hdulist.append(tbhdu)
            
            #count models 
            nmod = nmod+1

    #at the end write to disk adding number of models 
    hdulist[0].header['NUMMODEL']=nmod
    hdulist.writeto(namefits,clobber=True)

    return

