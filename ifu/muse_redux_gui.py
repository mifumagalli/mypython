#Redux gui

def reduxgui(listimg,mode='align',refcat='None',cubexsuffix='2'):

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

                finalmask -> let's you select a region on the coadded file 
                             and applies this mask to all the region files. Useful to trim the 
                             edges of all the exposures to look the same. Works as maskcubex.

                maskmpdaf -> add mask option to mpdaf
    
    refcat (optional) --> an ascii file containing a list of RA and DEC of 
                          reference objects (stars) which acts as an absolute 
                          reference for each frame
                 
    cubexsuffix  --> the suffix that specifies the cycle of cubex reduction 

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
    #from matplotlib.colors import LogNorm
    from matplotlib.colors import Normalize
    from astropy.visualization import ZScaleInterval
    from astropy import wcs
    import numpy as np
    import PyGuide
    import os
    from shutil import copyfile
    import math
    import string

    ##define the widget class here
    class guivalues():
       
       def __init__(self):
        
        self.mjd     = []
        self.starra  = []
        self.stardec = []
        self.starid  = []
        self.boxlx   = []
        self.boxly   = []
        self.boxrx   = []
        self.boxry   = []
    
    
    class align_tk(Tkinter.Tk):
     
        def __init__(self,filename,parent,guivalues,mode='align',refra=[],refdec=[],refid=[]):
            
            self.root=Tkinter.Tk.__init__(self,parent)
            self.parent = parent
            self.filename=filename
            self.refra=refra
            self.refdec=refdec
            self.refid=refid
            self.mode=mode
            self.guivalues=guivalues
            
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
            self.menuaspect=[1.,0.24]
            self.dataaspect=[1.,1.-0.24]
            self.dpi=80
        
            #Fiddle with font
            default_font = tkFont.nametofont("TkDefaultFont")
            scalefont = int(screen_height/1080.0*14)     #Ruari 29/08 fixes bug where different resolutions cause the menu to be cut off 
            default_font.configure(size=scalefont)
            
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
                self.mjd=self.fits[0].header['MJD-OBS'] 
            except:
                self.mjd=None
            try:
                self.fitimg=np.nan_to_num(self.fits[1].data)+0.1
                self.wcs=wcs.WCS(self.fits[1].header)
            except:
                self.fitimg=np.nan_to_num(self.fits[0].data)+0.1
                self.wcs=wcs.WCS(self.fits[0].header)

            #init vmin/vmax and variables 
            zscale = ZScaleInterval()
            self.vmin, self.vmax = zscale.get_limits(self.fitimg)
            #self.vmin=np.median(self.fitimg)-1.*np.std(self.fitimg)
            #if(self.vmin <= 0):
            #    self.vmin=0.05
            #self.vmax=np.median(self.fitimg)+1.*np.std(self.fitimg)
            self.vminv.set('{}'.format(self.vmin))
            self.vmaxv.set('{}'.format(self.vmax))

            self.starx=[]
            self.stary=[]
            self.starra=[]
            self.stardec=[]
            self.starid=[]
            self.nstars=0
            self.boxlx=self.guivalues.boxlx
            self.boxly=self.guivalues.boxly
            self.boxrx=self.guivalues.boxrx
            self.boxry=self.guivalues.boxry
            self.nbox=len(self.guivalues.boxry)
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
            #upgraded to zscale
            #zscale = ZScaleInterval()
            #zmin, zmax = scale.get_limits(self.fitimg)
            self.twodimagePlot_prop["image"] =self.twodimagePlot_prop["axis"].imshow(self.fitimg,origin='lower', norm=Normalize(vmin=self.vmin,vmax=self.vmax) ,aspect='auto') 
            #norm=SymLogNorm(linthresh=2.0, linscale=0.2, vmin=self.vmin,vmax=self.vmax)
            #norm=LogNorm(vmin=self.vmin,vmax=self.vmax)
            self.twodimagePlot_prop["image"].set_cmap('gray')
            
            if(len(self.starx) > 0):
                self.twodimagePlot_prop["axis"].scatter(self.starx,self.stary,s=60,color='red')
            
            for ii in range(self.nbox):
                self.twodimagePlot_prop["axis"].plot([self.boxlx[ii],self.boxrx[ii]],[self.boxly[ii],self.boxly[ii]],color='red')
                self.twodimagePlot_prop["axis"].plot([self.boxrx[ii],self.boxrx[ii]],[self.boxly[ii],self.boxry[ii]],color='red')
                self.twodimagePlot_prop["axis"].plot([self.boxlx[ii],self.boxrx[ii]],[self.boxry[ii],self.boxry[ii]],color='red')
                self.twodimagePlot_prop["axis"].plot([self.boxlx[ii],self.boxlx[ii]],[self.boxly[ii],self.boxry[ii]],color='red')
                self.twodimagePlot_prop["axis"].text(self.boxrx[ii],self.boxry[ii],"{}".format(ii+1),color='red')
            
            if(len(self.refid) > 0):
                for ii in range(len(self.refid)):
                    rx,ry=self.wcs.wcs_world2pix(self.refra[ii],self.refdec[ii],1)
                    if(rx>0 and ry>0 and rx<=np.shape(self.fitimg)[1] and ry<=np.shape(self.fitimg)[0]):
                       self.twodimagePlot_prop["axis"].text(rx,ry,"{}".format(self.refid[ii]),color='red', ha='center', va='center')
                       self.twodimagePlot_prop["axis"].scatter(rx,ry,s=200,facecolors='none',edgecolors='red')

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
                x = event.xdata
                y = event.ydata
                self.mouse_position.set('Cursor: xy=({0:7.2f},{1:7.2f}); radec=({2},{3})'.format(x,y,ra,dec))
            except:
                ra=0
                dec=0
                self.mouse_position.set('Cursor: xy=(None,None); radec=(None,None)')
                 
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
        
        def find_id(self,ra,dec):
            
            #Distance squared
            try:
              dist = (self.refra-ra)**2 + (self.refdec-dec)**2
              return self.refid[np.argmin(dist)]
            except:
              return self.nstars
             
                
        def presskey(self,event):
            
            """ Do stuff when press  key  """
            
            #quit on q
            if(event.key == "q"):
                self.OnExit()
            
            #centroid on c
            elif(event.key == "c"):
                ccd=PyGuide.CCDInfo(0,0,1)
                centroid=PyGuide.centroid(self.fitimg,None,None,[event.xdata,event.ydata],20,ccd)
                self.starx.append(centroid.xyCtr[0])
                self.stary.append(centroid.xyCtr[1])
                thisra,thisdec=self.wcs.wcs_pix2world(centroid.xyCtr[0],centroid.xyCtr[1],1)
                self.starra.append(thisra)
                self.stardec.append(thisdec)
                self.starid.append(self.find_id(thisra,thisdec))
                self.nstars=self.nstars+1
                self.update_twodimage(update=True)
            
            #store bottom left corner on L 
            elif(event.key == "n"):
                try:
                    check=float(event.xdata)+float(event.ydata)
                    self.boxlx.append(event.xdata)
                    self.boxly.append(event.ydata)
                    self.warning.set("STATUS: Now mark top/right corner!")
                except:
                    self.boxlx.append(event.xdata)
                    self.boxly.append(event.ydata)
                    self.warning.set("STATUS: Missed data region... try again!")
            elif(event.key == "m"):
                try:
                    check=float(event.xdata)+float(event.ydata)
                    self.boxrx.append(event.xdata)
                    self.boxry.append(event.ydata)
                    self.nbox=self.nbox+1
                    self.update_twodimage(update=True)
                    self.warning.set("STATUS: All good!")
                except:
                    self.warning.set("STATUS: Missed data region... try again!")
            
            elif(event.key == "d"):
                if(self.nbox > 0):
                  self.boxlx=self.boxlx[0:-1]
                  self.boxly=self.boxly[0:-1]
                  self.boxrx=self.boxrx[0:-1]
                  self.boxry=self.boxry[0:-1]
                  self.nbox=self.nbox-1
                if(self.nstars > 0):
                  self.starra=self.starra[0:-1] 
                  self.stardec=self.stardec[0:-1]
                  self.starid=self.starid[0:-1]
                  self.starx=self.starx[0:-1] 
                  self.stary=self.stary[0:-1] 
                  self.nstars=self.nstars-1
                self.update_twodimage(update=True)
                
        def OnExit(self):
            """ Quit all on exit """
            
            #Sort star IDs before return
            ind = np.argsort(np.array(self.starid))
            
            self.guivalues.starra  = np.array(self.starra)[ind]  
            self.guivalues.stardec = np.array(self.stardec)[ind]
            self.guivalues.starid  = np.array(self.starid)[ind] 
            self.guivalues.boxlx = self.boxlx
            self.guivalues.boxly = self.boxly
            self.guivalues.boxrx = self.boxrx
            self.guivalues.boxry = self.boxry
            self.guivalues.mjd   = self.mjd
            
            self.quit()
            self.destroy()


    #continue with the main programme
    #read image list and loop over aligment
    #check modes and run accordingly:
    
    if(mode is 'align'): 
        
        #Check if an external star catalogue exists and is readable
        if refcat != 'None':
          if not os.path.isfile(refcat):
             print('File {0} not found. Abort!'.format(refcat))
             return
          
          print('REDUXGUI: Run in align mode with external catalogue')
          ra_ref = []
          dec_ref = []
          id_ref= []
          first=False
          
          with open(refcat) as f:
             for ii, line in enumerate(f):
                if(line[0] not in ['#',' ','','\n']):
                  temp = string.strip(line).split(None)
                  ra_ref.append(float(temp[0])) 
                  dec_ref.append(float(temp[1])) 
                  id_ref.append(ii)
          
          ra_ref  = np.array(ra_ref)
          dec_ref = np.array(dec_ref)
          id_ref  = np.array(id_ref)
           
        else:
          print('REDUXGUI: Run in align mode with first exp as reference')
          first=True
        
        ra_off=[]
        dec_off=[]
        mjd=[]
        
        print('REDUXGUI: Loop over images to align')
        for ii in open(listimg):
            name=ii.split(" ")[0]
            GUIvalues = guivalues()
            if(first):
                ##main loop for first image to select  
                app = align_tk(name,None,GUIvalues,mode=mode)
                app.title('Align {}'.format(name))
                app.mainloop()
                ra_off.append(0.0)
                dec_off.append(0.0)
                ra_ref=np.copy(GUIvalues.starra)
                dec_ref=np.copy(GUIvalues.stardec)
                id_ref=np.copy(GUIvalues.starid)
                mjd.append(GUIvalues.mjd)
                first=False
            else:
                ##main loop for other exposures
                app = align_tk(name,None,GUIvalues,mode=mode,refra=ra_ref,refdec=dec_ref,refid=id_ref)
                app.title('Align {}'.format(name))
                app.mainloop()
                thisraoff = np.mean(-ra_ref[np.in1d(id_ref,GUIvalues.starid)]+GUIvalues.starra)
                thisdecoff = np.mean(-dec_ref[np.in1d(id_ref,GUIvalues.starid)]+GUIvalues.stardec)
                print('REDUXGUI: offset derived for {0}: RA = {1:6.4f}" DEC = {2:6.4f}"'.format(name,thisraoff*3600,thisdecoff*3600))
                ra_off.append(thisraoff)
                dec_off.append(thisdecoff)
                mjd.append(GUIvalues.mjd)
        
        #save in fits after making backup copy 
        print('REDUXGUI: write offsets to OFFSET_LIST.fits')
        copyfile('OFFSET_LIST.fits','OFFSET_LIST.fits.bck')
        strcut=fits.open('OFFSET_LIST.fits',mode='update')
        
        mjdfile = (strcut[1].data[:]['MJD_OBS'])
         
        for xx in range(len(dec_off)):
            offlistind = np.argmin(abs(mjdfile-mjd[xx]))
            strcut[1].data[offlistind]['RA_OFFSET']=ra_off[xx]
            strcut[1].data[offlistind]['DEC_OFFSET']=dec_off[xx]

        strcut.flush()
        strcut.close()
        
    elif(mode is 'maskcubex'): 
        print('REDUXGUI: Run in maskcubex mode')
	print(listimg)
        for ii in open(listimg):
            #check if running on intermediate or final step
            if('hsn' in ii):
                #override flag to right extension
                cubexsuffix='hsn'
        
            whiteimg="_".join(ii.split("_")[0:-1])+"_white{}.fits".format(cubexsuffix)
            region="_".join(ii.split("_")[0:-1])+"_fix{}_SliceEdgeMask.reg".format(cubexsuffix)
            
            GUIvalues = guivalues()
            
            #if region file exists load it back 
            if os.path.isfile(region):
                print('Loading existing region file {}'.format(region))
                regload=open(region,'r')
                for line in regload:
                    if('box' in line):
                        segment=line.split('(')[1].split(')')[0]
                        xcent,ycent,twodx,twody,zero=segment.split(',')
                        GUIvalues.boxlx.append(float(xcent)-float(twodx)/2.)
                        GUIvalues.boxly.append(float(ycent)-float(twody)/2.)
                        GUIvalues.boxrx.append(float(xcent)+float(twodx)/2.)
                        GUIvalues.boxry.append(float(ycent)+float(twody)/2.)
                regload.close()      
                
            #run gui
            app = align_tk(whiteimg,None,GUIvalues,mode=mode)
            app.title('Mask {}'.format(whiteimg))
            app.mainloop()
        
            #on exit, dump to region file
            if(len(GUIvalues.boxrx) > 0):
                rfl=open(region,'w')
                rfl.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
                rfl.write('image\n')
                for bb in range(len(GUIvalues.boxrx)):
                    dx=(GUIvalues.boxrx[bb]-GUIvalues.boxlx[bb])/2.
                    dy=(GUIvalues.boxry[bb]-GUIvalues.boxly[bb])/2.
                    rfl.write('box({},{},{},{},0)\n'.format(GUIvalues.boxlx[bb]+dx,GUIvalues.boxly[bb]+dy,2.*dx,2*dy))
                rfl.close()
                    
    elif(mode is 'finalmask'): 
        print('REDUXGUI: Run in finalmask mode')                #RUari 25/05/17 added to do final step where the field is finally cropped, this runs off the cubexCOMBINED_IMAGE.fits and appends the regions to all exopsures

        app = align_tk('COMBINED_IMAGE.fits',None,mode='maskcubex')
        app.title('Mask {}'.format('COMBINED_IMAGE.fits'))
        app.mainloop()

        for ii in open(listimg):
            region="_".join(ii.split("_")[0:-1])+"_fix2_SliceEdgeMask.reg"

            #if region file exists then load and used for inital start
            if (os.path.isfile(region)):
                rfl=open(region,'a')
                for bb in range(len(boxrx)):
                    dx=(boxrx[bb]-boxlx[bb])/2.
                    dy=(boxry[bb]-boxly[bb])/2.
                    rfl.write('box({},{},{},{},0)\n'.format(boxlx[bb]+dx,boxly[bb]+dy,2.*dx,2*dy))
                rfl.close()
        
    elif(mode is 'maskmpdaf'): 
        print('REDUXGUI: Run in maskmpdaf mode')

        for ii in open(listimg):
            
            whiteimg=ii.strip()
            region=ii.split("fits")[0]+"reg"
            
            GUIvalues = guivalues()
            
            #if region file exists load it back 
            if os.path.isfile(region):
                print('Loading existing region file {}'.format(region))
                regload=open(region,'r')
                for line in regload:
                    if('box' in line):
                        segment=line.split('(')[1].split(')')[0]
                        xcent,ycent,twodx,twody,zero=segment.split(',')
                        GUIvalues.boxlx.append(float(xcent)-float(twodx)/2.)
                        GUIvalues.boxly.append(float(ycent)-float(twody)/2.)
                        GUIvalues.boxrx.append(float(xcent)+float(twodx)/2.)
                        GUIvalues.boxry.append(float(ycent)+float(twody)/2.)
                regload.close()      
                
            #run gui
            app = align_tk(whiteimg,None,GUIvalues,mode=mode)
            app.title('Mask {}'.format(whiteimg))
            app.mainloop()
        
            #on exit, dump to region file
            if(len(GUIvalues.boxrx) > 0):
                rfl=open(region,'w')
                rfl.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
                rfl.write('image\n')
                for bb in range(len(GUIvalues.boxrx)):
                    dx=(GUIvalues.boxrx[bb]-GUIvalues.boxlx[bb])/2.
                    dy=(GUIvalues.boxry[bb]-GUIvalues.boxly[bb])/2.
                    rfl.write('box({},{},{},{},0)\n'.format(GUIvalues.boxlx[bb]+dx,GUIvalues.boxly[bb]+dy,2.*dx,2*dy))
                rfl.close()
            
        
    else:
        print ('Mode {} not found!'.format(mode))
