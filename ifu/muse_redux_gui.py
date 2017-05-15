#global variable
global starra
global stardec
global starid
global boxlx
global boxly
global boxrx
global boxry

def reduxgui(listimg,mode='align'):

    """

    GUI designed to help with semi-automatic tasks involved in the reduction of 
    muse data.

    listimg --> an ascii file containing the images to perfom operations on
    mode    --> the mode that defines the code behaviour
    
                align -> take a list of white images (generally form align.sof) 
                         and let the user mark stars in interactive mode to compute and
                         write offsets to disk in "OFFSET_LIST.fits"

                maskcubex -> take a list of cubex (generally from cubes.lst) and 
                             enable masking of boxes which are written in format 
                             that can be handled by cubex combine procedure. 

    Interactive behaviour: 

         . 'c': mark star
         . 'q': quit and save
         . RMB: hold and move: stretch image scaling
         . n/m: mark region to mask by selecting left-bottom and top-right corner
         . d  : delete last region

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
    import math

    ##define the widget class here
    class align_tk(Tkinter.Tk):
     
        def __init__(self,filename,parent,mode='align',refra=[],refdec=[],refid=[]):
            self.root=Tkinter.Tk.__init__(self,parent)
            self.parent = parent
            self.filename=filename
            self.refra=refra
            self.refdec=refdec
            self.refid=refid
            self.mode=mode
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
            self.menuaspect=[1.,0.15]
            self.dataaspect=[1.,1.-0.15]
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
            self.mouse_position_w.grid(column=0,row=2,sticky='W',columnspan=4)

            self.helpm=Tkinter.Label(self.menuframe, text='Mark stars with "c"; save and exit with "q"; draw boxes with n/m; d to delete last')
            self.helpm.grid(column=0,row=0,sticky='W',columnspan=4)

            self.modev=Tkinter.Label(self.menuframe,text='Running in mode "{}"'.format(self.mode))
            self.modev.grid(column=0,row=1,sticky='W',columnspan=4)

            
            self.zlab=Tkinter.Label(self.menuframe,text='ZRANGE=')
            self.zlab.grid(column=0,row=3,sticky='W')

            self.vminv=Tkinter.StringVar()
            self.vminv.set('0.1')
            self.vminv_w=Tkinter.Entry(self.menuframe,textvariable = self.vminv)
            self.vminv_w.grid(column=1,row=3,sticky='W')

            self.vmaxv=Tkinter.StringVar()
            self.vmaxv.set('10.0')
            self.vmaxv_w=Tkinter.Entry(self.menuframe,textvariable = self.vmaxv)
            self.vmaxv_w.grid(column=2,row=3,sticky='W')

            self.apply=Tkinter.Button(self.menuframe,text="Apply",command=self.setscale)
            self.apply.grid(column=3,row=3)

            
            self.warning=Tkinter.StringVar()
            self.warning.set('STATUS: All good!')
            self.warning_w=Tkinter.Label(self.menuframe,textvariable = self.warning)
            self.warning_w.grid(column=0,row=4,sticky='W',columnspan=4)

        def init_dataframe(self):

            #work out dimensions for twod image
            self.twodimg_width=self.dataframe.winfo_width()
            self.twodimg_height=self.dataframe.winfo_height()
           
            #now open image
            self.fits=fits.open(self.filename)
            try:
                self.fitimg=np.nan_to_num(self.fits[1].data)+0.1
                self.wcs=wcs.WCS(self.fits[1].header)
            except:
                self.fitimg=np.nan_to_num(self.fits[0].data)+0.1
                self.wcs=wcs.WCS(self.fits[0].header)

            #init vmin/vmax and variables 
            self.vmin=np.median(self.fitimg)-1.*np.std(self.fitimg)
            if(self.vmin <= 0):
                self.vmin=0.05
            self.vmax=np.median(self.fitimg)+1.*np.std(self.fitimg)
            self.vminv.set('{}'.format(self.vmin))
            self.vmaxv.set('{}'.format(self.vmax))

            self.centx=[]
            self.centy=[]
            self.nstars=1
            self.nbox=1
            global starra
            global stardec
            global starid
            global boxlx
            global boxly
            global boxrx
            global boxry
            starra=[]
            stardec=[]
            starid=[]
            boxlx=[]
            boxly=[]
            boxrx=[]
            boxry=[]
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

        def setscale(self,*args):

            #update the vmin/vmax
            self.vmin=float(self.vminv.get())
            self.vmax=float(self.vmaxv.get())

            #trigger update
            self.update_twodimage(update=True)

        def update_twodimage(self,update=False):
        
            """
            Code that updates the 2D image
            Update = True, redraw
            
            """
            if(update):
                self.twodimagePlot_prop["axis"].cla()
                self.vminv.set('{}'.format(self.vmin))
                self.vmaxv.set('{}'.format(self.vmax))

            self.twodimagePlot_prop["image"] =self.twodimagePlot_prop["axis"].imshow(self.fitimg,origin='lower',norm=LogNorm(vmin=self.vmin,vmax=self.vmax),aspect='auto')
            self.twodimagePlot_prop["image"].set_cmap('gray')
            
            if(len(self.centx) > 0):
                self.twodimagePlot_prop["axis"].plot(self.centx,self.centy,'o',color='red')
            
            for ii in range(len(boxrx)):
                self.twodimagePlot_prop["axis"].plot([boxlx[ii],boxrx[ii]],[boxly[ii],boxly[ii]],color='red')
                self.twodimagePlot_prop["axis"].plot([boxrx[ii],boxrx[ii]],[boxly[ii],boxry[ii]],color='red')
                self.twodimagePlot_prop["axis"].plot([boxlx[ii],boxrx[ii]],[boxry[ii],boxry[ii]],color='red')
                self.twodimagePlot_prop["axis"].plot([boxlx[ii],boxlx[ii]],[boxly[ii],boxry[ii]],color='red')
                self.twodimagePlot_prop["axis"].text(boxrx[ii],boxry[ii],\
                                                         "{}".format(ii+1),color='red')
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
                minstarts=np.median(self.fitimg)-1.*np.std(self.fitimg)
                if(minstarts <= 0):
                    minstarts=0.1
                maxstart=np.median(self.fitimg)+1.*np.std(self.fitimg)
                
                self.vmin=minstarts*(1.-(self.preferwinwidth-event.x)/(1.*self.preferwinwidth))
                self.vmax=maxstart*(1.-(self.preferwinheight-event.y)/(1.*self.preferwinheight))
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
            #store bottom left corner on L 
            elif(event.key == "n"):
                global boxlx
                global boxly
                try:
                    check=float(event.xdata)+float(event.ydata)
                    boxlx.append(event.xdata)
                    boxly.append(event.ydata)
                    self.warning.set("STATUS: Now mark top/right corner!")
                except:
                    boxlx.append(event.xdata)
                    boxly.append(event.ydata)
                    self.warning.set("STATUS: Missed data region... try again!")
            elif(event.key == "m"):
                global boxrx
                global boxry
                try:
                    check=float(event.xdata)+float(event.ydata)
                    boxrx.append(event.xdata)
                    boxry.append(event.ydata)
                    self.nbox=self.nbox+1
                    self.update_twodimage(update=True)
                    self.warning.set("STATUS: All good!")
                except:
                    self.warning.set("STATUS: Missed data region... try again!")
            elif(event.key == "d"):
                global boxrx
                global boxry
                global boxlx
                global boxly
                boxlx=boxlx[0:-1]
                boxly=boxly[0:-1]
                boxrx=boxrx[0:-1]
                boxry=boxry[0:-1]
                self.nbox=self.nbox-1
                self.update_twodimage(update=True)
                
        def OnExit(self):
            """ Quit all on exit """
            self.quit()
            self.destroy()


    #continue with the main programme
    #read image list and loop over aligment
    global starid
    global starra
    global stardec
    global boxlx
    global boxly
    global boxrx
    global boxry

    #check modes and run accordingly:
    
    if(mode is 'align'): 
        print('REDUXGUI: Run in align mode')
        ra_off=[]
        dec_off=[]
        starid=[]
        starra=[]
        stardec=[]
        first=True
    
        print('REDUXGUI: Loop over images to align')
        for ii in open(listimg):
            name=ii.split(" ")[0]
            if(first):
                ##main loop for first image to select  
                app = align_tk(name,None,mode=mode)
                app.title('Align {}'.format(name))
                app.mainloop()
                ra_off.append(0.0)
                dec_off.append(0.0)
                ra_ref=np.copy(starra)
                dec_ref=np.copy(stardec)
                id_ref=np.copy(starid)
                first=False
            else:
                ##main loop for other exposures
                app = align_tk(name,None,mode=mode,refra=ra_ref,refdec=dec_ref,refid=id_ref)
                app.title('Align {}'.format(name))
                app.mainloop()
                ra_off.append(np.mean(-ra_ref+starra))
                dec_off.append(np.mean(-dec_ref+stardec))
                
        
        #save in fits after making backup copy 
        print('REDUXGUI: write offsets to OFFSET_LIST.fits')
        copyfile('OFFSET_LIST.fits','OFFSET_LIST.fits.bck')
        strcut=fits.open('OFFSET_LIST.fits',mode='update')
       
        for xx in range(len(dec_off)):
            strcut[1].data[xx]['RA_OFFSET']=ra_off[xx]
            strcut[1].data[xx]['DEC_OFFSET']=dec_off[xx]

        strcut.flush()
        strcut.close()
    elif(mode is 'maskcubex'): 
        print('REDUXGUI: Run in maskcubex mode')
        for ii in open(listimg):
            whiteimg="_".join(ii.split("_")[0:-1])+"_white2.fits"
            region="_".join(ii.split("_")[0:-1])+"_fix2_SliceEdgeMask.reg"

            #run gui
            app = align_tk(whiteimg,None,mode=mode)
            app.title('Mask {}'.format(whiteimg))
            app.mainloop()
        
            #on exit, dump to region file
            if(len(boxrx) > 0):
                rfl=open(region,'w')
                rfl.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
                rfl.write('image\n')
                for bb in range(len(boxrx)):
                    dx=(boxrx[bb]-boxlx[bb])/2.
                    dy=(boxry[bb]-boxly[bb])/2.
                    rfl.write('box({},{},{},{},0)\n'.format(boxlx[bb]+dx,boxly[bb]+dy,2.*dx,2*dy))
                rfl.close()
                    
    else:
        print ('Mode {} not found!'.format(mode))
