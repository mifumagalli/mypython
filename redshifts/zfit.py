"""

Gui to inspect spectra in 1/2D 


"""

import Tkinter
import tkFont
from Tkinter import Tk
import tkFileDialog

from astropy.io import fits
import matplotlib 
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import os
import numpy as np
import scipy
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy import interpolate
from scipy import signal
from astropy.io import fits
from astropy.table import Table



class zfitwin(Tkinter.Tk):
    
    """ The basic class of the widget """

    def __init__(self,parent):
        
        """ My constructor """
        
        self.tk = Tk()

        #set min and preferred size of main gui
        self.minwinwidth=300
        self.minwinheight=300
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        self.preferwinwidth=int(screen_width*0.8)
        self.preferwinheight=int(screen_height*0.8)
        self.minsize(width=self.minwinwidth, height=self.minwinheight)
        self.geometry("{}x{}".format(self.preferwinwidth,self.preferwinheight))
        
        #tweak the aspect ratio of the menu and data gui
        self.menuaspect=[1,0.24*(screen_height/1080.0)]     #Ruari 24/05 fixes bug where different resolutions cause the menu to be cut off 
        self.dataaspect=[1,1-0.24*(screen_height/1080.0)]     #Ruari 24/05 fixes bug where different resolutions cause the menu to be cut off 
        self.dpi=80

        #find exect dir
        self.execdir=__file__.split('zfit.py')[0]
        if(len(self.execdir)==0):
            self.execdir='./'
            
        #Fiddle with font
        default_font = tkFont.nametofont("TkDefaultFont")
        default_font.configure(size=14)

        #init gui frame
        self.initialize()

    def initialize(self):
        """ This init the basic gui """ 
        
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

        #stick the 2D image in a separate window
        self.imgframe=Tkinter.Toplevel(width=600,height=600)
        
        #update for later use of units 
        self.update()

        #now initialise the menu frame
        self.init_menuframe()
        #now initialise the data frame
        self.init_dataframe()

    def init_menuframe(self):

        """ This init the menu specific part """ 
        
        #exit button
        self.menu_exit = Tkinter.Button(self.menuframe,text=u"EXIT",command=self.OnExit)
        self.menu_exit.grid(column=0,row=0)
      
        #save button
        self.menu_save = Tkinter.Button(self.menuframe,text=u"Save",command=self.OnSave)
        self.menu_save.grid(column=0,row=1)

        #choice of spectra
        self.menu_select = Tkinter.Button(self.menuframe,text=u"Open Spectrum",
                                          command=self.SelectFile)
        self.menu_select.grid(column=0,row=2)
        
        #current spectrum
        self.currspec=Tkinter.StringVar()
        self.currspec.set('Spect: Demo')
        self.current=Tkinter.Label(self.menuframe,textvariable = self.currspec)
        self.current.grid(column=0,row=3)

        self.mouse_position=Tkinter.StringVar()
        self.mouse_position.set('Mouse:(None,None)')
        self.mouse_position_w=Tkinter.Label(self.menuframe,textvariable = self.mouse_position)
        self.mouse_position_w.grid(column=0,row=4,columnspan=3)

        #Message window
        self.generic_message=Tkinter.StringVar()
        self.generic_message.set('zfit-> Ready to go!')
        self.generic_message_w=Tkinter.Label(self.menuframe,textvariable = self.generic_message)
        self.generic_message_w.grid(column=5,row=3,columnspan=3)

        #line control stuff
        self.init_linecontrol()

        #templates control stuff 
        self.init_templcontrol()


    def init_dataframe(self):

        """ This init the data specific part """ 

        #Work out the geometry of the different data parts

        #canvas for spectrum ... 
        self.pltspec_width=self.dataframe.winfo_width()
        self.pltspec_height=int(self.dataframe.winfo_height()*0.6)
    
        #canvas for twod spec
        self.twodspc_width=self.dataframe.winfo_width()
        self.twodspc_height=int((self.dataframe.winfo_height()-self.pltspec_height)*0.6)

        #canvas for twod err
        self.twoderr_width=self.dataframe.winfo_width()
        self.twoderr_height=int((self.dataframe.winfo_height()-self.pltspec_height)*0.4)
        
        #work out dimensions for twod image
        self.twodimg_width=self.imgframe.winfo_width()
        self.twodimg_height=self.imgframe.winfo_height()

        #now open with default spectrum and plot
        self.filename=os.path.abspath(self.execdir)+"/test_spectrum.fits"
        self.fits=fits.open(self.filename)
       
        #unpack
        self.fitwav1d=self.fits[2].data
        self.fitspe1d=self.fits[0].data
        self.fitspe1d_original=np.copy(self.fitspe1d)
        self.fiterr1d=self.fits[1].data
        self.fitspe2d=self.fits[4].data
        self.fiterr2d=self.fits[5].data
        self.fitimg=self.fits[6].data
        
        self.drawdata()

        #set tmpfitxcorr to None to avoid error or later init
        self.tmpfitxcorr=None
        #set smoothwindow
        self.smooth=3
        

    def init_linecontrol(self):
        
        """ This controls operation with emission lines """

        #just say what it is 
        linelabel=Tkinter.Label(self.menuframe,text = "Emission lines")
        linelabel.grid(column=1,row=0,columnspan=2)

        #drop down menu to select emission lines
        llab = Tkinter.Label(self.menuframe, text="Select Lines: ")
        llab.grid(column=1,row=1)
        self.linelist = Tkinter.StringVar(self.menuframe)
        self.linelist.set("gal_vac") # default value
        self.lineselect = Tkinter.OptionMenu(self.menuframe, self.linelist,"gal_vac","gal_air","lbg","lls","tell")
        self.lineselect.grid(column=2,row=1)
        #set the linelist in trace state
        self.linelist.trace("w",self.displaylines)

        #line redshift window 
        zlab = Tkinter.Label(self.menuframe, text="z = ")
        zlab.grid(column=1,row=2)
        self.redshiftline = Tkinter.StringVar()
        self.redlinecntr = Tkinter.Entry(self.menuframe,textvariable=self.redshiftline)
        self.redlinecntr.grid(column=2,row=2)
        self.redshiftline.set("0.0000")
        #set the redshift in a trace state
        self.redshiftline.trace("w",self.displaylines)

        #display lines
        self.shwlinstate=Tkinter.IntVar()
        self.lineshow = Tkinter.Checkbutton(self.menuframe, text="Show Lines",
                                            variable=self.shwlinstate,command=self.displaylines)
        self.lineshow.grid(column=1,row=3)
        
        #fit lines
        self.line_fit = Tkinter.Button(self.menuframe,text=u"FitLines",command=self.fitlines)
        self.line_fit.grid(column=2,row=3)
      
    def init_templcontrol(self):

        """ Control the options for template fitting """
        
        #just say what it is 
        templabel=Tkinter.Label(self.menuframe,text = "Templates")
        templabel.grid(column=3,row=0,columnspan=4)

        #drop down menu to select template family 
        llab = Tkinter.Label(self.menuframe, text="Pick template: ")
        llab.grid(column=3,row=1)
        self.tempgroup= Tkinter.StringVar(self.menuframe)
        self.tempgroup.set("Select")
        self.tempselect = Tkinter.OptionMenu(self.menuframe,self.tempgroup,"kinney","lbgs","sdss")
        self.tempselect.grid(column=4,row=1)
        self.tempgroup.trace("w",self.loadtemplate)

        #just say what it is 
        self.currenttemplate=Tkinter.StringVar(self.menuframe)
        self.currenttemplate.set("Current: None")
        self.tempchoice=Tkinter.Label(self.menuframe,textvariable = self.currenttemplate)
        self.tempchoice.grid(column=5,row=1,columnspan=2)

        #D not use trace for template, as these are expensive to compute 
        #template redshift window 
        zlab = Tkinter.Label(self.menuframe, text="z = ")
        zlab.grid(column=3,row=2)
        self.redshifttemp = Tkinter.StringVar()
        self.redtempcntr = Tkinter.Entry(self.menuframe,textvariable=self.redshifttemp)
        self.redtempcntr.grid(column=4,row=2)
        self.redshifttemp.set("0.0000")
         
        #rmag window
        rmg = Tkinter.Label(self.menuframe, text="flux = ")
        rmg.grid(column=3,row=3)
        self.magtemp = Tkinter.StringVar()
        self.magtemcntr = Tkinter.Entry(self.menuframe,textvariable=self.magtemp)
        self.magtemcntr.grid(column=4,row=3)
        self.magtemp.set("1.00")
    
        #display template
        self.shwtempstate=Tkinter.IntVar()
        self.tempshow = Tkinter.Button(self.menuframe,text="Show Template",command=self.displaytemplate)
        self.tempshow.grid(column=3,row=4)
        self.temphide = Tkinter.Button(self.menuframe,text="Hide Template",command=self.hidetemplate)
        self.temphide.grid(column=4,row=4)
       
        #fit template
        self.template_fit = Tkinter.Button(self.menuframe,text=u"FitTemplate",command=self.fittemplate)
        self.template_fit.grid(column=5,row=2)
    
    def OnExit(self):
        """ Quit all on exit """
        self.fits.close()
        self.quit()
        self.destroy()

    def OnSave(self):
        """ Save screen """
        print 'Placeholder'


    def SelectFile(self):
        """ Select and open file as one wishes """
        #select file 
        self.filename=tkFileDialog.askopenfilename(initialdir='./')
        #update name
        self.currspec.set("Spec: "+self.filename.split("/")[-1])
        
        #close old and reopen
        self.fits.close()
        self.fits=fits.open(self.filename)
        
        #unpack
        self.fitwav1d=self.fits[2].data
        self.fitspe1d=self.fits[0].data
        self.fitspe1d_original=np.copy(self.fits[0].data)
        self.fiterr1d=self.fits[1].data
        self.fitspe2d=self.fits[4].data
        self.fiterr2d=self.fits[5].data
        self.fitimg=self.fits[6].data

        #redraw
        self.drawdata(refresh=True)
        
    def drawdata(self,refresh=False):
        
        """
        Once the spectrum is set, populate the data part of the gui
        
        refresh -> True, wipe all canvas before redrawing
  
        """
        
        if(refresh):
            #now destroy all data canvas 
            self.twodimagePlot.get_tk_widget().destroy()
            self.spectrumPlot.get_tk_widget().destroy()
            self.twodspcPlot.get_tk_widget().destroy()
            self.twoderrPlot.get_tk_widget().destroy()

        #refresh 2D image
        self.init_twodimage()

        #refresh the spectrum
        self.init_spectrum()

        #refresh 2D spec
        self.init_twodspec()
        
        #refresh 2D err
        self.init_twoderr()
        

    def init_twodimage(self):

        """ Draw the 2D image """

        #create properties for this plot
        self.twodimagePlot_prop={}
        
        #figure staff
        self.twodimagePlot_prop["figure"] = Figure(figsize=(self.twodimg_width/self.dpi,self.twodimg_height/self.dpi),
                                                   dpi=self.dpi)
        self.twodimagePlot_prop["axis"] = self.twodimagePlot_prop["figure"].add_subplot(111)
 
        #call plotting routine
        self.update_twodimage()
     
        #send it to canvas - connect event 
        self.twodimagePlot = FigureCanvasTkAgg(self.twodimagePlot_prop["figure"],master=self.imgframe)
        self.twodimagePlot.show()
        #need to set tight layout after showing 
        self.twodimagePlot_prop["figure"].tight_layout()
        #enable event on click
        self.twodimagePlot.mpl_connect("button_press_event", self.pressbutton)
        self.twodimagePlot.mpl_connect("key_press_event", self.presskey)

        self.twodimagePlot.get_tk_widget().grid()

    def update_twodimage(self,update=False):
        
        """
        Code that updates the 2D image
        Update = True, redraw
                
        """
        self.twodimagePlot_prop["image"] =self.twodimagePlot_prop["axis"].imshow(self.fitimg,origin='lower',aspect='auto')
        self.twodimagePlot_prop["image"].set_cmap('hot')
        #self.twodimagePlot_prop["axis"].set_xlabel('Pix')
        #self.twodimagePlot_prop["axis"].set_ylabel('Pix')

        
    def init_spectrum(self):

        """ Draw the spectrum """
        
        #create properties for this plot
        self.spectrumPlot_prop={}
        self.spectrumPlot_prop["xmin"]=np.min(np.nan_to_num(self.fitwav1d))
        self.spectrumPlot_prop["xmax"]=np.max(np.nan_to_num(self.fitwav1d))
        self.spectrumPlot_prop["ymin"]=np.min(np.nan_to_num(self.fitspe1d))
        self.spectrumPlot_prop["ymax"]=np.max(np.nan_to_num(self.fitspe1d))

        #figure stuff
        self.spectrumPlot_prop["figure"]= Figure(figsize=(0.99*self.pltspec_width/self.dpi,0.96*self.pltspec_height/self.dpi),
                                                 dpi=self.dpi)
        self.spectrumPlot_prop["axis"]= self.spectrumPlot_prop["figure"].add_subplot(111)

        #call plotting routine
        self.update_spectrum()
        #send it to canvas
        self.spectrumPlot = FigureCanvasTkAgg(self.spectrumPlot_prop["figure"],master=self.dataframe)
        self.spectrumPlot.show()
        #enable event on click
        self.spectrumPlot_prop["figure"].tight_layout()
        self.spectrumPlot.mpl_connect("button_press_event", self.pressbutton)
        self.spectrumPlot.mpl_connect("motion_notify_event", self.movemouse)
        self.spectrumPlot.mpl_connect("key_press_event", self.presskey)
        self.spectrumPlot.get_tk_widget().grid(column=0,row=0)

    def update_spectrum(self,update=False):
        
        """

        Code that updates the spectrum
        
        Update = True, redraw

        """
        if(update):
            self.spectrumPlot_prop["axis"].cla()
            
        #plot main data
        self.spectrumPlot_prop["axis"].step(self.fitwav1d,self.fitspe1d,where='mid')
        self.spectrumPlot_prop["axis"].step(self.fitwav1d,self.fiterr1d,color='red',\
                                            linestyle='--',zorder=1,where='mid')
        self.spectrumPlot_prop["axis"].set_xlim(self.spectrumPlot_prop["xmin"],self.spectrumPlot_prop["xmax"])
        self.spectrumPlot_prop["axis"].set_ylim(self.spectrumPlot_prop["ymin"],self.spectrumPlot_prop["ymax"])
        self.spectrumPlot_prop["axis"].set_xlabel('Wavelength')
        #self.spectrumPlot_prop["axis"].set_ylabel('Flux')

        #if needed, plot lines
        if(self.shwlinstate.get()):
             #set redshift
            try:
                redsh=float(self.redshiftline.get())
            except:
                redsh=0.0
            #loop over lines and draw
            for lw,lnam in self.infoline:
                #find the obs wave
                lwplot=lw*(1+redsh)
                if((lwplot > self.spectrumPlot_prop["xmin"]) & (lwplot < self.spectrumPlot_prop["xmax"])):
                    self.spectrumPlot_prop["axis"].axvline(lwplot, color='grey', linestyle='--')
                    self.spectrumPlot_prop["axis"].text(lwplot,self.spectrumPlot_prop["ymax"],lnam,
                                                        verticalalignment='top',rotation=90,fontsize=12)
  
        #if needed, plot template
        if(self.shwtempstate.get()):
            self.spectrumPlot_prop["axis"].plot(self.fitwav1d,self.templatedata_current,color='black',zorder=3)
                  
        #plot zero line
        self.spectrumPlot_prop["axis"].plot([self.spectrumPlot_prop["xmin"],self.spectrumPlot_prop["xmax"]],
                                            [0,0],color='green',zorder=2,linestyle=':')


        #finally draw
        if(update):
            self.spectrumPlot.draw()
    

    def init_twodspec(self):

        """ Draw the 2D spectrum """

        #create properties for this plot
        self.twodspcPlot_prop={}

        #figure staff
        self.twodspcPlot_prop["figure"]= Figure(figsize=(0.99*self.twodspc_width/self.dpi,0.96*self.twodspc_height/self.dpi),
                                                dpi=self.dpi)
        self.twodspcPlot_prop["axis"] = self.twodspcPlot_prop["figure"].add_subplot(111)

        #call plotting routine
        self.update_twodspec()
    
        #send it to canvas
        self.twodspcPlot = FigureCanvasTkAgg(self.twodspcPlot_prop["figure"],master=self.dataframe)
        self.twodspcPlot.show()
        #enable event on click
        self.twodspcPlot_prop["figure"].tight_layout()
        self.twodspcPlot.mpl_connect("button_press_event", self.pressbutton)
        self.twodspcPlot.mpl_connect("key_press_event", self.presskey)
        self.twodspcPlot.mpl_connect("motion_notify_event", self.movemouse)

        self.twodspcPlot.get_tk_widget().grid(column=0,row=1,sticky='NW')
        
    def wavemap(self,x,pos):
    
        """ Utility to map the pixel in 2D image to wavelegth """
        
        #wavelength mapper 
        index=np.arange(0,len(self.fitwav1d))
        wave=np.interp(x,index,self.fitwav1d)
        'The two args are the value and tick position'
        return "%.1f" % wave
    
    def inv_wavemap(self,x):
    
        """ Utility to map wavelegth to pixel in 2D mage """
        
        #wavelength mapper 
        index=np.arange(0,len(self.fitwav1d))
        pix=np.interp(x,self.fitwav1d,index,left=0,right=len(self.fitwav1d))
        return  pix

    def update_twodspec(self,update=False):
        
        """

        Code that updates the 2D spectrum
        
        Update = True, redraw

        """
        
        if(update):
            self.twodspcPlot_prop["axis"].cla()
    
        self.twodspcPlot_prop["image"]=self.twodspcPlot_prop["axis"].imshow(np.rot90(self.fitspe2d),origin='lower',aspect='auto')
        self.twodspcPlot_prop["image"].set_cmap('hot')
       
        
        #control level
        medianlevel=np.median(np.nan_to_num(self.fitspe2d))
        stdlevel=np.std(np.nan_to_num(self.fitspe2d))
        self.twodspcPlot_prop["image"].set_clim(medianlevel-3.*stdlevel,medianlevel+3*stdlevel)

        #wave mapper
        self.twodspcPlot_prop["axis"].xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(self.wavemap))
        
        #now set X axis as in 1d spectrum
        xpixmin=self.inv_wavemap(self.spectrumPlot_prop["xmin"])
        xpixmax=self.inv_wavemap(self.spectrumPlot_prop["xmax"])
        
        #force minimum maximum
        if(xpixmin == xpixmax):
            xpixmin = xpixmax-1
        if(xpixmax == 0):
            xpixmax = 1

        self.twodspcPlot_prop["axis"].set_xlim(xpixmin,xpixmax)
        self.twodspcPlot_prop["axis"].set_xlabel('Wavelength')

        if(update):
            self.twodspcPlot.draw()
    

    def init_twoderr(self):

        """ Draw the 2D error """

        #create properties for this plot
        self.twoderrPlot_prop={}
        
        #figure staff
        #self.twoderr.grid(column=1,row=2,sticky='NW')
        self.twoderrPlot_prop['figure'] = Figure(figsize=(0.99*self.twoderr_width/self.dpi,0.96*self.twoderr_height/self.dpi),
                                                 dpi=self.dpi)
        self.twoderrPlot_prop['axis'] = self.twoderrPlot_prop['figure'].add_subplot(111)
        
        #call plotting routine
        self.update_twoderr()
    
        #send it to canvas
        self.twoderrPlot = FigureCanvasTkAgg(self.twoderrPlot_prop['figure'],master=self.dataframe)
        self.twoderrPlot.show()
        #enable event on click
        self.twoderrPlot_prop['figure'].tight_layout()
        self.twoderrPlot.mpl_connect("button_press_event", self.pressbutton)
        self.twoderrPlot.mpl_connect("key_press_event", self.presskey)
        self.twoderrPlot.mpl_connect("motion_notify_event", self.movemouse)
        self.twoderrPlot.get_tk_widget().grid(column=0,row=2,sticky='NW')


    def update_twoderr(self,update=False):
        
        """

        Code that updates the 2D error
        
        Update = True, redraw

        """

        if(update):
            self.twoderrPlot_prop["axis"].cla()
        
        self.twoderrPlot_prop['image'] =self.twoderrPlot_prop['axis'].imshow(np.rot90(self.fiterr2d),origin='lower',aspect='auto')
        self.twoderrPlot_prop['image'].set_cmap('hot')
        

        #control level
        medianlevel=np.median(np.nan_to_num(self.fiterr2d))
        stdlevel=np.std(np.nan_to_num(self.fiterr2d))
        self.twoderrPlot_prop["image"].set_clim(medianlevel-3.*stdlevel,medianlevel+3*stdlevel)


        #wave mapper
        self.twoderrPlot_prop["axis"].xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(self.wavemap))
        
        #now set X axis as in 1d spectrum
        xpixmin=self.inv_wavemap(self.spectrumPlot_prop["xmin"])
        xpixmax=self.inv_wavemap(self.spectrumPlot_prop["xmax"])
        
        #force minimum maximum
        if(xpixmin == xpixmax):
            xpixmin = xpixmax-1
        if(xpixmax == 0):
            xpixmax = 1

        self.twoderrPlot_prop["axis"].set_xlim(xpixmin,xpixmax)
        self.twoderrPlot_prop["axis"].set_xlabel('Wavelength')

        if(update):
            self.twoderrPlot.draw()


    def displaylines(self,*args):

        """ Display the line list by refreshing plot in update state """
 
        #first parse the line lists
        linefile=self.execdir+"/lines/"+self.linelist.get()+".lst"
        self.infoline = Table.read(linefile, format='ascii.basic')
        #self.infoline=np.loadtxt(linefile, dtype={'names': ('wave', 'tag'),
        #                                          'formats': ('f4', 'S4')})
        
        #refresh plot
        self.update_spectrum(update=True)
    
    def loadtemplate(self,*args):
        
        """ Load template from disk and preselect some
            useful default                   
        """

        #if so, start dialogue to pick the desired template
        self.picktemplate=tkFileDialog.askopenfilename(initialdir='{}/templates/{}'.format(self.execdir,self.tempgroup.get()))
        #set current template 
        self.currenttemplate.set("Current: "+self.picktemplate.split("/")[-1])
        
        #load current template
        if('sdss' in self.tempgroup.get()):
            #load fits
            fitstemp=fits.open(self.picktemplate)
            #grab flux
            self.templatedata={'flux':fitstemp[0].data[0,:]}
            #cosntruct wave
            waveinx=np.arange(0,len(self.templatedata['flux']),1)
            wavevac=10**(waveinx*1.*fitstemp[0].header['COEFF1']+1.*fitstemp[0].header['COEFF0'])
            ##go to air
            #self.templatedata['wave']= wavevac/(1.0+2.735182e-4+131.4182/wavevac**2+2.76249e8/wavevac**4)
            #remain in vac 
            self.templatedata['wave']= wavevac
        else:
            #load text
            #self.templatedata=np.loadtxt(self.picktemplate, dtype={'names': ('wave', 'flux'),
            #                                                       'formats': ('f10', 'f10')},usecols=(0,1))
            self.templatedata = Table.read(self.picktemplate, format='ascii.basic')
            
        #set sensible pick in redshift and adjust data as needed
        if('lbg' in self.tempgroup.get()):
            self.redshifttemp.set("3.000")
        elif('kinney' in self.tempgroup.get()):
            self.templatedata['flux']=self.templatedata['flux']/1e-14
        elif('sdss' in self.tempgroup.get()):
            self.templatedata['flux']=self.templatedata['flux']*100.
        else:
            self.redshifttemp.set("0.000")
                
    def displaytemplate(self,*args):

        """ Compute and display template """
        self.shwtempstate.set(1)

        #compute template given current values
        self.adapttemplate()
                        
        #refresh plot
        self.update_spectrum(update=True)

    def hidetemplate(self,*args):

        """ Hide template """
        self.shwtempstate.set(0)
        #refresh plot
        self.update_spectrum(update=True)

    def adapttemplate(self):
        
        """ Interpolate a template over the data """

        #redshift factor 
        redhfactor=(1+float(self.redshifttemp.get()))
   
        #now construct interpolation 
        thisw=self.templatedata['wave']*redhfactor
        thisf=self.templatedata['flux']
        intflx = interp1d(thisw,thisf,kind='linear',bounds_error=False,fill_value=0.0)
        #apply normalisation
        self.templatedata_current=intflx(self.fitwav1d)*float(self.magtemp.get())
    

    def fitlines(self):

        """ Fit the line list """

        #loop over lines inside spectrum 
        #lounch a new window
        self.lnfit=Tkinter.Toplevel(self.tk)
        
        #add a display
        fig=Figure(figsize=(self.preferwinwidth/self.dpi,self.preferwinheight/self.dpi),dpi=self.dpi)
                 
        #pick z
        try:
            redsh=float(self.redshiftline.get())
        except:
            redsh=0.0
            
        lines_good_wave_rest=[]
        lines_good_wave_obs=[]
        lines_good_name=[]

        for lw,lnam in self.infoline:
            lwplot=lw*(1+redsh)
            if((lwplot > min(self.fitwav1d)+8) & (lwplot < max(self.fitwav1d)-8)):
                
                #do a boxchart in 6A bin to see if line exists 
                inside=np.where((self.fitwav1d > lwplot-4)& (self.fitwav1d < lwplot+4))
                continuum=np.where(((self.fitwav1d > lwplot-20)& (self.fitwav1d < lwplot-10)) |
                                   ((self.fitwav1d > lwplot+10)& (self.fitwav1d < lwplot+20)))
                clevel=np.median(self.fitspe1d[continuum])
                flux=np.sum((self.fitspe1d[inside]-clevel))
                noise=np.sqrt(np.sum(self.fiterr1d[inside]**2))
                
                #cut in SN
                if(flux/noise > 2):
                    
                    #stash 
                    lines_good_wave_rest.append(lw)
                    lines_good_wave_obs.append(lwplot)
                    lines_good_name.append(lnam)
                    

        #generate a 4x? grid of plots
        nlines=len(lines_good_wave_rest)
        ncol=4
        nraw=int(nlines/ncol)
        if(nlines%ncol > 0):
            nraw=nraw+1
        
        czall=[]

        #loop on good stuff for fits
        for ii in range(nlines):

            #select region to fit
            fitwindow=np.where((self.fitwav1d > lines_good_wave_obs[ii]-10) & (self.fitwav1d < lines_good_wave_obs[ii]+10))
            continuum=np.where(((self.fitwav1d > lines_good_wave_obs[ii]-20)& (self.fitwav1d < lines_good_wave_obs[ii]-10)) |
                               ((self.fitwav1d > lines_good_wave_obs[ii]+10)& (self.fitwav1d < lines_good_wave_obs[ii]+20)))
            clevel=np.median(self.fitspe1d[continuum])
            p0=np.array([10.,1.*float(lines_good_wave_obs[ii]),2.,0.])

            #fit a Gaussian 
            yval=np.nan_to_num(self.fitspe1d[fitwindow]-clevel)
            yerr=np.nan_to_num(self.fiterr1d[fitwindow]*1.)
            xval=np.nan_to_num(self.fitwav1d[fitwindow]*1.)
            popt,pcov=curve_fit(self.gauss,xval,yval,p0=p0, sigma=yerr)
            perr = np.sqrt(np.diag(pcov))
            #eval fit 
            xg=np.arange(min(xval)-2,max(xval)+2,0.2)
            fitg=self.gauss(xg,*popt)

            #grab fits
            czfit=popt[1]/lines_good_wave_rest[ii]-1.
            czfiterr=perr[1]/lines_good_wave_rest[ii]
            czall.append(czfit)

            #display
            ax = fig.add_subplot(nraw,ncol,ii+1)
            ax.plot(xval,yval)
            ax.plot(xval,yerr,color='red',linestyle="--",zorder=1)
            ax.plot(xg,fitg,color='black',linestyle=":")
            ax.set_title("{0}{1} z = {2:.5} +/- {3:.5}".format(lines_good_name[ii],int(lines_good_wave_rest[ii]),czfit,czfiterr))

        #send message to user  and reset redshift 
        bestz=np.median(np.array(czall))
        bestez=np.std(np.array(czall))
        self.generic_message.set(r'zfit-> Best fit is {:6.5f}+/-{:6.5f}'.format(bestz,bestez))
        self.redshiftline.set(bestz)

        #send figure to canvas
        self.linefitplot = FigureCanvasTkAgg(fig,master=self.lnfit)
        self.linefitplot.show()
        #fig.tight_layout()
        self.linefitplot.get_tk_widget().grid()
 
    def fittemplate(self):

        """ Fit the template  """

        #init the template correlation 
        realdata={'wave':self.fitwav1d,'flux':self.fitspe1d,'error':self.fiterr1d}
        
        ##Testing sequence
        #realdata={'wave':self.templatedata['wave']*(1+0.4329),'flux':self.templatedata['flux'], 
        #          'error':self.templatedata['flux']}
        
        print 'Computing correlation... be patient!'
            
        #find the wavelength range covering the min/max extent 
        absmin=np.min([np.min(self.templatedata['wave']),np.min(realdata['wave'])])
        absmax=np.max([np.max(self.templatedata['wave']),np.max(realdata['wave'])])

        #resample in log 
        deltal=5e-4
        lnwave=np.arange(np.log(absmin),np.log(absmax),deltal)
        
        #resample with spline (s controls the smoothing)
        x=np.nan_to_num(self.templatedata['wave'])
        y=np.nan_to_num(self.templatedata['flux'])
        resamp_templ=interpolate.splrep(np.log(x),y,s=0)
        x=np.nan_to_num(realdata['wave'])
        y=np.nan_to_num(realdata['flux'])
        resamp_real=interpolate.splrep(np.log(x),y,s=0)
       
        #put everything on the same array - zero padding the extrapolation
        flux_templ=interpolate.splev(lnwave,resamp_templ,der=0,ext=1)
        flux_real=interpolate.splev(lnwave,resamp_real,der=0,ext=1)

        #masking strong sky lines
        mask=np.where((lnwave > np.log(5569.)) & (lnwave < np.log(5584.)))
        flux_real[mask]=0
        mask=np.where((lnwave > np.log(6292.)) & (lnwave < np.log(6308.)))
        flux_real[mask]=0
        mask=np.where((lnwave > np.log(6356.)) & (lnwave < np.log(6369.)))
        flux_real[mask]=0
        mask=np.where((lnwave > 8.6752) & (lnwave < 8.6860))
        flux_real[mask]=0
        mask=np.where((lnwave > 8.8274) & (lnwave < 8.8525))
        flux_real[mask]=0
        mask=np.where((lnwave > 8.8862) & (lnwave < np.log(12000.)))
        flux_real[mask]=0

        #correlate 
        xcorr=np.correlate(flux_real,flux_templ,mode='full')

        #find the peak in the second half in units of redshift
        indxmax=np.argmax(xcorr)-len(xcorr)/2
        peakz=np.exp(indxmax*deltal)-1
        #print peakz
        
        #find the reshift axis
        indxarr=np.arange(0,len(lnwave),1)
        self.xcorr_redax=np.exp(indxarr*deltal)-1
        self.xcorr_xcorr=xcorr[len(xcorr)/2:]
        self.xcorr_redshift=peakz

        #set the redshift in template window
        self.redshifttemp.set("{}".format(self.xcorr_redshift))

        #trigger display options 
        #lounch a new window
        self.tmlfit=Tkinter.Toplevel(self.tk)
        
        #add xcorr to display
        #create properties for this plot
        self.tmpfitxcorr_prop={}
        self.tmpfitxcorr_prop["xmin"]=np.min(np.nan_to_num(self.xcorr_redax))
        self.tmpfitxcorr_prop["xmax"]=np.max(np.nan_to_num(self.xcorr_redax))
        self.tmpfitxcorr_prop["ymin"]=np.min(np.nan_to_num(self.xcorr_xcorr))
        self.tmpfitxcorr_prop["ymax"]=np.max(np.nan_to_num(self.xcorr_xcorr))


        self.tmpfitxcorr_prop["figure"]=Figure(figsize=(self.preferwinwidth/self.dpi*0.75,self.preferwinheight/self.dpi*0.75),dpi=self.dpi)
        self.tmpfitxcorr_prop["axis"]= self.tmpfitxcorr_prop["figure"].add_subplot(111)

        #call plotting routine
        self.update_xcorrplot()
   
        #send it to canvas
        self.tmpfitxcorr = FigureCanvasTkAgg(self.tmpfitxcorr_prop["figure"],master=self.tmlfit)
        self.tmpfitxcorr.show()
        #enable event on click
        self.tmpfitxcorr_prop["figure"].tight_layout()
        self.tmpfitxcorr.mpl_connect("button_press_event", self.pressbutton)
        self.tmpfitxcorr.mpl_connect("key_press_event", self.presskey)
        self.tmpfitxcorr.get_tk_widget().grid(column=0,row=0)

        
    def update_xcorrplot(self,update=False):

        """ Update plot for xcorrplot """

        if(update):
            self.tmpfitxcorr_prop["axis"].cla()
            
        #plot main data
        self.tmpfitxcorr_prop["axis"].plot(self.xcorr_redax,self.xcorr_xcorr)
        self.tmpfitxcorr_prop["axis"].axvline(self.xcorr_redshift, color='grey', linestyle='--')
        self.tmpfitxcorr_prop["axis"].set_xlim(self.tmpfitxcorr_prop["xmin"],self.tmpfitxcorr_prop["xmax"])
        self.tmpfitxcorr_prop["axis"].set_ylim(self.tmpfitxcorr_prop["ymin"],self.tmpfitxcorr_prop["ymax"])
        self.tmpfitxcorr_prop["axis"].set_xlabel('Redshift')
        self.tmpfitxcorr_prop["axis"].set_ylabel('XCORR')

        #finally draw
        if(update):
            self.tmpfitxcorr.draw()
    

    def movemouse(self,event):
        
        """ Do stuff when  mouse moves """
        if(event.canvas == self.spectrumPlot): 
            self.mouse_position.set('Mouse:({},{})'.format(event.xdata,event.ydata))
        elif(event.canvas == self.twodspcPlot):
            try:
                self.mouse_position.set('Mouse:({},{})'.format(self.wavemap(event.xdata,0.0),event.ydata))
            except:
                self.mouse_position.set('Mouse:(None,None)')
        elif(event.canvas == self.twoderrPlot):
            try:
                self.mouse_position.set('Mouse:({},{})'.format(self.wavemap(event.xdata,0.0),event.ydata))
            except:
                self.mouse_position.set('Mouse:(None,None)')
                
        
    def pressbutton(self,event):
        
        """ Do stuff when data plot is pressed with mouse """

        #this is how to redirect events
        if(event.canvas == self.twoderrPlot): 
            #set focus
            self.twoderrPlot.get_tk_widget().focus_set()
        if(event.canvas == self.twodspcPlot): 
            #set focus
            self.twodspcPlot.get_tk_widget().focus_set()
        if(event.canvas == self.twodimagePlot): 
            #set focus
            self.twodimagePlot.get_tk_widget().focus_set()
        if(event.canvas == self.spectrumPlot): 
            #set focus
            self.spectrumPlot.get_tk_widget().focus_set()
            #for right click, trigger line selector
            if(event.button == 3):
                self.lineselectorwidget(event)
        if(event.canvas == self.tmpfitxcorr): 
            #set focus
            self.tmpfitxcorr.get_tk_widget().focus_set()

    def presskey(self,event):
        
        """ Do stuff when data plot is pressed with key """

        #quit on q
        if(event.key == "q"):
            self.OnExit()

        #keyboard event when focus on spectrum
        if(event.canvas == self.spectrumPlot): 
            self.spetrumPlot_events(event)
            
        #keyboard event when focus on xcorr
        if(event.canvas == self.tmpfitxcorr): 
            self.tmpfitxcorr_events(event)
 
    def tmpfitxcorr_events(self,event):

        """ Handle events of xcorr plot """
        
        

        #set bottom plot
        if(event.key == "b"):
            self.tmpfitxcorr_prop["ymin"]=event.ydata
            self.update_xcorrplot(update=True)
            #set top plot
        if(event.key == "t"):
            self.tmpfitxcorr_prop["ymax"]=event.ydata
            self.update_xcorrplot(update=True)
            #set left plot
        if(event.key == "l"):
            self.tmpfitxcorr_prop["xmin"]=event.xdata
            self.update_xcorrplot(update=True)
         #set right plot
        if(event.key == "r"):
            self.tmpfitxcorr_prop["xmax"]=event.xdata
            self.update_xcorrplot(update=True)
         #zoom in
        if(event.key == "i"):
             #find the current width in x
            currentwidth=self.tmpfitxcorr_prop["xmax"]-self.tmpfitxcorr_prop["xmin"]
             #zoom in by factor of 2
            currentwidth=currentwidth*0.5
             #zoom around selected wave
            self.tmpfitxcorr_prop["xmin"]=event.xdata-currentwidth/2.
            self.tmpfitxcorr_prop["xmax"]=event.xdata+currentwidth/2.
            self.update_xcorrplot(update=True)
         #zoom out
        if(event.key == "o"):
             #find the current width in x
            currentwidth=self.tmpfitxcorr_prop["xmax"]-self.tmpfitxcorr_prop["xmin"]
             #zoom out by factor of 2
            currentwidth=currentwidth*2
             #zoom around selected wave
            self.tmpfitxcorr_prop["xmin"]=event.xdata-currentwidth/2.
            self.tmpfitxcorr_prop["xmax"]=event.xdata+currentwidth/2.
            self.update_xcorrplot(update=True)
    
         #pan left 
        if(event.key == "["):
             #find the current width in x
            currentwidth=self.tmpfitxcorr_prop["xmax"]-self.tmpfitxcorr_prop["xmin"]
             #pan left
            self.tmpfitxcorr_prop["xmin"]=self.tmpfitxcorr_prop["xmin"]-currentwidth/2
            self.tmpfitxcorr_prop["xmax"]=self.tmpfitxcorr_prop["xmax"]-currentwidth/2
            self.update_xcorrplot(update=True)
  
         #pan right 
        if(event.key == "]"):
             #find the current width in x
            currentwidth=self.tmpfitxcorr_prop["xmax"]-self.tmpfitxcorr_prop["xmin"]
             #pan right
            self.tmpfitxcorr_prop["xmin"]=self.tmpfitxcorr_prop["xmin"]+currentwidth/2
            self.tmpfitxcorr_prop["xmax"]=self.tmpfitxcorr_prop["xmax"]+currentwidth/2
            self.update_xcorrplot(update=True)
         #set reset plot
        if(event.key == "W"):
            self.tmpfitxcorr_prop["xmin"]=np.min(np.nan_to_num(self.xcorr_redax))
            self.tmpfitxcorr_prop["xmax"]=np.max(np.nan_to_num(self.xcorr_redax))
            self.tmpfitxcorr_prop["ymin"]=np.min(np.nan_to_num(self.xcorr_xcorr))
            self.tmpfitxcorr_prop["ymax"]=np.max(np.nan_to_num(self.xcorr_xcorr))
            self.update_xcorrplot(update=True)
            
        #mark new redshift
        if(event.key == "z"):
            #update relevent info
            self.xcorr_redshift=event.xdata
            self.redshifttemp.set("{}".format(self.xcorr_redshift))
            #refresh plot
            self.update_xcorrplot(update=True)
            #display template
            self.displaytemplate()


    def spetrumPlot_events(self,event):
    
         """" Handle events of spectrum plot """
         
          #set bottom plot
         if(event.key == "b"):
             self.spectrumPlot_prop["ymin"]=event.ydata
             self.update_spectrum(update=True)
            #set top plot
         if(event.key == "t"):
             self.spectrumPlot_prop["ymax"]=event.ydata
             self.update_spectrum(update=True)
            #set left plot
         if(event.key == "l"):
             self.spectrumPlot_prop["xmin"]=event.xdata
             self.update_spectrum(update=True)
             #update 2d spectra accordingly
             self.update_twodspec(update=True)
             self.update_twoderr(update=True)
         #set right plot
         if(event.key == "r"):
             self.spectrumPlot_prop["xmax"]=event.xdata
             self.update_spectrum(update=True)
             #update 2d spectra accordingly
             self.update_twodspec(update=True)             
             self.update_twoderr(update=True)        
         #zoom in
         if(event.key == "i"):
             #find the current width in x
             currentwidth=self.spectrumPlot_prop["xmax"]-self.spectrumPlot_prop["xmin"]
             #zoom in by factor of 2
             currentwidth=currentwidth*0.5
             #zoom around selected wave
             self.spectrumPlot_prop["xmin"]=event.xdata-currentwidth/2.
             self.spectrumPlot_prop["xmax"]=event.xdata+currentwidth/2.
             self.update_spectrum(update=True)
             #update 2d spectra accordingly
             self.update_twodspec(update=True)
             self.update_twoderr(update=True)              
         #zoom out
         if(event.key == "o"):
             #find the current width in x
             currentwidth=self.spectrumPlot_prop["xmax"]-self.spectrumPlot_prop["xmin"]
             #zoom out by factor of 2
             currentwidth=currentwidth*2
             #zoom around selected wave
             self.spectrumPlot_prop["xmin"]=event.xdata-currentwidth/2.
             self.spectrumPlot_prop["xmax"]=event.xdata+currentwidth/2.
             self.update_spectrum(update=True)
             #update 2d spectra accordingly
             self.update_twodspec(update=True)
             self.update_twoderr(update=True)

         #pan left 
         if(event.key == "["):
             #find the current width in x
             currentwidth=self.spectrumPlot_prop["xmax"]-self.spectrumPlot_prop["xmin"]
             #pan left
             self.spectrumPlot_prop["xmin"]=self.spectrumPlot_prop["xmin"]-currentwidth/2
             self.spectrumPlot_prop["xmax"]=self.spectrumPlot_prop["xmax"]-currentwidth/2
             self.update_spectrum(update=True)
             #update 2d spectra accordingly
             self.update_twodspec(update=True)
             self.update_twoderr(update=True)


         #pan right 
         if(event.key == "]"):
             #find the current width in x
             currentwidth=self.spectrumPlot_prop["xmax"]-self.spectrumPlot_prop["xmin"]
             #pan right
             self.spectrumPlot_prop["xmin"]=self.spectrumPlot_prop["xmin"]+currentwidth/2
             self.spectrumPlot_prop["xmax"]=self.spectrumPlot_prop["xmax"]+currentwidth/2
             self.update_spectrum(update=True)
             #update 2d spectra accordingly
             self.update_twodspec(update=True)
             self.update_twoderr(update=True)
         #set reset plot
         if(event.key == "W"):
             self.spectrumPlot_prop["xmin"]=np.min(np.nan_to_num(self.fitwav1d))
             self.spectrumPlot_prop["xmax"]=np.max(np.nan_to_num(self.fitwav1d))
             self.spectrumPlot_prop["ymin"]=np.min(np.nan_to_num(self.fitspe1d))
             self.spectrumPlot_prop["ymax"]=np.max(np.nan_to_num(self.fitspe1d))
             self.update_spectrum(update=True)
             #update 2d spectra accordingly
             self.update_twodspec(update=True)
             self.update_twoderr(update=True)
        #smooth plot
         if(event.key == "S"):
             self.fitspe1d=signal.medfilt(self.fitspe1d,self.smooth)
             self.smooth=self.smooth+2
             self.update_spectrum(update=True)
            
        #unsmooth smooth 
         if(event.key == "U"):
             self.fitspe1d=self.fitspe1d_original
             self.smooth=3
             self.update_spectrum(update=True)

         
    def lineselectorwidget(self,event):
        """  Control what happens when right-click on 1D spectrum 
        
             - trigger construction of line list selector
        
        """

        #refresh lines as needed
        self.displaylines()

        #lounch a new window
        self.lnsel=Tkinter.Toplevel(self.tk)
     
        #pick z
        try:
            redsh=float(self.redshiftline.get())
        except:
            redsh=0.0
            
        #create line buttons for those visibles
        self.wlineselect = Tkinter.DoubleVar()
        self.wlinepos = event.xdata
        i=0
        for lw,lnam in self.infoline:
            lwplot=lw*(1+redsh)
            Tkinter.Radiobutton(self.lnsel, text=lnam+"{}".format(int(lw)), 
                                variable=self.wlineselect, value=lw,
                                command=self.pickedline).grid(row = i%30, column = i/30, sticky = "NWSE")
            i=i+1

        self.tk.wait_window(self.lnsel)

    def pickedline(self):
        
        """ If one pick a line, find redshift """

        #find the redshift
        redshift=self.wlinepos/self.wlineselect.get()-1

        #set it - auto trigger refresh 
        self.shwlinstate.set(1)        
        self.redshiftline.set("{}".format(redshift))

        #destroy window
        self.lnsel.destroy()
        
    def gauss(self,x, *p):
        """ Gaussian model for line fit """
        
        A, mu, sigma, zero = p
        gg=A*np.exp(-1.*(x-mu)*(x-mu)/(2.*sigma*sigma))+zero
        
        return gg


    


def zfit():

    """ Mains that runs the gui """
    app = zfitwin(None)
    app.title('Fit your redshift!')
    app.mainloop()

if __name__ == "__main__": 
    zfit()

