"""

Collection of procedures that drive the creation of final cubes using
eso recipies. Check out the MUSE reduction manual for details. 

"""

def individual_skysub(listob,nproc=12):

    """
    rerun scipost performing sky subtraction with eso recipies 
    

    listobs -> a list of folders corresponding to individual OBs that needs to be 
               processed enabling skysubtraction

    nproc -> number of processors to be used in parallel mode 

    """
    
    import os
    import glob
    import subprocess

    #grab top dir
    topdir=os.getcwd()

    #now loop over each folder and make the final sky-subtracted cubes
    for ob in listob:
        
        #change dir
        os.chdir(ob+'/Proc/')
        print('Processing {} for sky subtraction'.format(ob))
        
        #Search how many exposures are there
        scils=glob.glob("OBJECT_RED_0*.fits*")
        nsci=len(scils)

        #loop on exposures and reduce frame with sky subtraction 
        for exp in range(nsci):
            
            #use sof file written for basic reduction
            sof_name="../Script/scipost_{0:d}.sof".format(exp+1)
           
            #define some output names
            cname="DATACUBE_FINAL_ESOSKY_EXP{0:d}.fits".format(exp+1)
            pname="PIXTABLE_REDUCED_ESOSKY_EXP{0:d}.fits".format(exp+1)
            iname="IMAGE_FOV_ESOSKY_EXP{0:d}.fits".format(exp+1)
 
            if not os.path.isfile(cname):
                print("Processing exposure {0:d}".format(exp+1))
                
                #Write the command file 
                scr=open("../Script/make_scipost_esosky_{0:d}.sh".format(exp+1),"w")
                scr.write("OMP_NUM_THREADS={0:d}\n".format(nproc)) 
                
                scr.write('esorex --log-file=scipost_esosky_{0:d}.log muse_scipost --filter=white --save=cube,individual ../Script/scipost_{0:d}.sof'.format(exp+1))
                scr.close()
    
                #Run pipeline 
                subprocess.call(["sh", "../Script/make_scipost_esosky_{0:d}.sh".format(exp+1)])    
                subprocess.call(["mv","DATACUBE_FINAL.fits",cname])
                subprocess.call(["mv","IMAGE_FOV_0001.fits",iname])
                subprocess.call(["mv","PIXTABLE_REDUCED_0001.fits",pname])
                                
            else:
                print("Exposure {0:d} exists.. skip! ".format(exp+1))
            

        #back to top
        os.chdir(topdir)
    

def coaddall(listob,nproc=24):

    """
    Grab all the OBs, align them with muse_exp_align, 
    and coadd them with muse_exp_combine 
    
    nproc -> number of proc to use in parallel runs

    """

    import os
    import glob
    import subprocess

    #create align and coadd scripts and sofs
    alignsof=open('align.sof','w')
    alignscr=open('align.sh','w')
    coaddsof=open('combine.sof','w')
    coaddscr=open('combine.sh','w')
 
    for ob in listob:
        
        #grab fov
        fovimg=glob.glob("../{}/Proc/IMAGE_FOV_ESOSKY_*".format(ob))
        for ff in fovimg:
            alignsof.write("{} IMAGE_FOV\n".format(ff))
    
        #grab pixel tables 
        pixtab=glob.glob("../{}/Proc/PIXTABLE_REDUCED_ESOSKY_*".format(ob))
        for ff in pixtab:
            coaddsof.write("{} PIXTABLE_REDUCED\n".format(ff))
        
    #final touches
    coaddsof.write("OFFSET_LIST.fits OFFSET_LIST\n")
    coaddsof.close()
    alignsof.close()
    
    #script files 
    alignscr.write("esorex --log-file=exp_align.log muse_exp_align --threshold=4. align.sof")
    alignscr.close()

    coaddscr.write("OMP_NUM_THREADS={}\n".format(nproc))
    coaddscr.write("esorex --log-file=exp_combine.log muse_exp_combine combine.sof")
    coaddscr.close()
    
    #Run pipeline 
    if not os.path.isfile('DATACUBE_FINAL.fits'):
        print("Aliging exposures...")
        subprocess.call(["sh", "align.sh"])    
    else:
        print("Offset already computed!")

    if not os.path.isfile('DATACUBE_FINAL.fits'):
        print("Coadding exposures...")
        subprocess.call(["sh", "combine.sh"])
    else:
        print("Cube already combined!")


#global variable
global starra
global stardec
global starid

def manual_align(listimg):

    """

    Take a list of white images and fires up a gui to help the process of image alignement 
    in case the default pipeline fails

    """

    import Tkinter
    import tkFont
    from astropy.io import fits
    import matplotlib 
    matplotlib.use("TkAgg")
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
    from matplotlib.figure import Figure
    from matplotlib.colors import LogNorm
    from astropy import wcs
    import numpy as np
    import PyGuide
    import os
    from shutil import copyfile

    ##define the widget class here
    class align_tk(Tkinter.Tk):
     
        def __init__(self,filename,parent,refra=[],refdec=[],refid=[]):
            self.root=Tkinter.Tk.__init__(self,parent)
            self.parent = parent
            self.filename=filename
            self.refra=refra
            self.refdec=refdec
            self.refid=refid
            #set min and preferred size of main gui
            self.minwinwidth=300
            self.minwinheight=300
            screen_width = self.winfo_screenwidth()
            screen_height = self.winfo_screenheight()
            self.preferwinwidth=int(screen_height*0.8)
            self.preferwinheight=int(screen_height*0.8)
            self.minsize(width=self.minwinwidth, height=self.minwinheight)
            self.geometry("{}x{}".format(self.preferwinwidth,self.preferwinheight))
            #tweak the aspect ratio of the menu and data gui
            self.menuaspect=[1.,0.1]
            self.dataaspect=[1.,1.-0.1]
            self.dpi=80
        

            #Fiddle with font
            default_font = tkFont.nametofont("TkDefaultFont")
            default_font.configure(size=14)


            self.initialize()
            
        def initialize(self):


            #create a menu frame
            self.menuframe=Tkinter.Frame(self,width=int(self.preferwinwidth*self.menuaspect[0]),
                                         height=int(self.preferwinheight*self.menuaspect[1]))
            self.menuframe.grid_propagate(0)
            self.menuframe.grid()
            
            #create a data frame
            self.dataframe=Tkinter.Frame(self,width=int(self.preferwinwidth*self.dataaspect[0]), 
                                         height=int(self.preferwinheight*self.dataaspect[1]))
            self.dataframe.grid_propagate(0)
            self.dataframe.grid()
            
            #update for later use of units 
            self.update()
            
            #now initialise the menu frame
            self.init_menuframe()
            #now initialise the data frame
            self.init_dataframe()


        def init_menuframe(self):

            self.mouse_position=Tkinter.StringVar()
            self.mouse_position.set('Mouse:(None,None)')
            self.mouse_position_w=Tkinter.Label(self.menuframe,textvariable = self.mouse_position)
            self.mouse_position_w.grid(column=0,row=0,columnspan=3)

            self.help=Tkinter.StringVar()
            self.help.set('Mark stars with "c"; save and exit with "q"')
            self.help_w=Tkinter.Label(self.menuframe,textvariable = self.help)
            self.help_w.grid(column=0,row=1,columnspan=3)


        def init_dataframe(self):

            #work out dimensions for twod image
            self.twodimg_width=self.dataframe.winfo_width()
            self.twodimg_height=self.dataframe.winfo_height()
           
            #now open image
            self.fits=fits.open(self.filename)
            self.fitimg=np.nan_to_num(self.fits[1].data)+0.1
            self.wcs=wcs.WCS(self.fits[1].header)

            #init vmin/vmax and variables 
            self.vmin=1.
            self.vmax=10.
            self.centx=[]
            self.centy=[]
            self.nstars=1
            global starra
            global stardec
            global starid
            starra=[]
            stardec=[]
            starid=[]
            self.drawdata()
            

        def drawdata(self,refresh=False):
                
            if(refresh):
            #now destroy all data canvas 
                self.twodimagePlot.get_tk_widget().destroy()
            #refresh 2D image
            self.init_twodimage()

        def init_twodimage(self):

            """ Draw the 2D image """
            
            #create properties for this plot
            self.twodimagePlot_prop={}
        
            #figure staff
            self.twodimagePlot_prop["figure"] = Figure(figsize=(1.*self.twodimg_width/self.dpi,\
                                                                   1.*self.twodimg_height/self.dpi),
                                                       dpi=self.dpi)
            self.twodimagePlot_prop["axis"] = self.twodimagePlot_prop["figure"].add_subplot(111)
            
            #call plotting routine
            self.update_twodimage()
     
            #send it to canvas - connect event 
            self.twodimagePlot = FigureCanvasTkAgg(self.twodimagePlot_prop["figure"],\
                                                       master=self.dataframe)
            self.twodimagePlot.show()
            #need to set tight layout after showing 
            #self.twodimagePlot_prop["figure"].tight_layout()
        
            #enable event on click
            #self.twodimagePlot.mpl_connect("button_press_event", self.pressbutton)
            self.twodimagePlot.mpl_connect("motion_notify_event", self.movemouse)
            self.twodimagePlot.mpl_connect("key_press_event", self.presskey)
      
            self.twodimagePlot.get_tk_widget().grid()

        def update_twodimage(self,update=False):
        
            """
            Code that updates the 2D image
            Update = True, redraw
            
            """
            if(update):
                self.twodimagePlot_prop["axis"].cla()

            self.twodimagePlot_prop["image"] =self.twodimagePlot_prop["axis"].imshow(self.fitimg,origin='lower',norm=LogNorm(vmin=self.vmin,vmax=self.vmax),aspect='auto')
            self.twodimagePlot_prop["image"].set_cmap('gray')
            #self.twodimagePlot_prop["axis"].set_xlabel('Pix')
            #self.twodimagePlot_prop["axis"].set_ylabel('Pix')
            
            if(len(self.centx) > 0):
                self.twodimagePlot_prop["axis"].plot(self.centx,self.centy,'o',color='red')
        
            if(len(self.refid) > 0):
                for ii,ll in enumerate(self.refid):
                    rx,ry=self.wcs.wcs_world2pix(self.refra[ii],self.refdec[ii],1)
                    self.twodimagePlot_prop["axis"].text(rx,ry,"{}".format(ll),color='red')

            xmax,ymax=self.fitimg.shape
            self.twodimagePlot_prop["axis"].set_xlim(0,ymax)
            self.twodimagePlot_prop["axis"].set_ylim(0,xmax)


            #finally draw
            if(update):
                self.twodimagePlot.draw()
    

        def movemouse(self,event):
            
            """ Do stuff when  mouse moves """
            
            try:
                ra,dec=self.wcs.wcs_pix2world(event.xdata,event.ydata,1)
            except:
                ra=0
                dec=0
            self.mouse_position.set('Cursor: xy=({},{}); radec=({},{})'.format(event.xdata,event.ydata,ra,dec))
                 
            #rescale image when mouse 3 down
            if(event.button == 3):
                
                #update vmin vmax and send to update                
                self.vmin=1*(1.-(self.preferwinwidth-event.x)/(1.*self.preferwinwidth))
                self.vmax=10*(1.-(self.preferwinheight-event.y)/(1.*self.preferwinheight))
                if(self.vmin <= 0):
                    self.vmin = 1e-5
                elif(self.vmin > self.vmax):
                    self.vmin=self.vmax
                self.update_twodimage(update=True)
                
        
        def presskey(self,event):
            
            """ Do stuff when press  key  """
            
            #quit on q
            if(event.key == "q"):
                self.OnExit()
            #centroid on c
            elif(event.key == "c"):
                global starra
                global stardec
                global starid
                ccd=PyGuide.CCDInfo(0,0,1)
                centroid=PyGuide.centroid(self.fitimg,None,None,[event.xdata,event.ydata],20,ccd)
                self.centx.append(centroid.xyCtr[0])
                self.centy.append(centroid.xyCtr[1])
                thisra,thisdec=self.wcs.wcs_pix2world(centroid.xyCtr[0],centroid.xyCtr[1],1)
                starra.append(thisra)
                stardec.append(thisdec)
                starid.append(self.nstars)
                self.nstars=self.nstars+1
                self.update_twodimage(update=True)
                

        def OnExit(self):
            """ Quit all on exit """
            self.quit()
            self.destroy()

    #continue with the main programme
    #read image list and loop over aligment
    ra_off=[]
    dec_off=[]
    global starid
    global starra
    global stardec

    starid=[]
    starra=[]
    stardec=[]
    first=True
    
    for ii in open(listimg):
        name=ii.split(" ")[0]
        if(first):
            ##main loop 
            app = align_tk(name,None)
            app.title('Align {}'.format(name))
            app.mainloop()
            ra_off.append(0.0)
            dec_off.append(0.0)
            ra_ref=np.copy(starra)
            dec_ref=np.copy(stardec)
            id_ref=np.copy(starid)
            first=False
        else:
            pass
            ##main loop 
            app = align_tk(name,None,refra=ra_ref,refdec=dec_ref,refid=id_ref)
            app.title('Align {}'.format(name))
            app.mainloop()
            ra_off.append(np.mean(-ra_ref+starra))
            dec_off.append(np.mean(-dec_ref+stardec))
        
        
        #save in fits after making backup copy 
        copyfile('OFFSET_LIST.fits','OFFSET_LIST.fits.bck')
        strcut=fits.open('OFFSET_LIST.fits',mode='update')
       
        for xx in range(len(dec_off)):
            strcut[1].data[xx]['RA_OFFSET']=ra_off[xx]
            strcut[1].data[xx]['DEC_OFFSET']=dec_off[xx]

        strcut.flush()
        strcut.close()

        #print ra_off, dec_off
