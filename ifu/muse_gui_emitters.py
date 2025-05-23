"""

GUI to handle visual inspection of emitters from MUSE
Assumes file written with the muse_emitters.py procedures

"""

try:
    import Tkinter as tkinter
    import tkFont as tkfont
    from Tkinter import Tk
    import tkFileDialog as filedialog
    import tkMessageBox
    import tkFileDialog
    
except:
    import tkinter
    from tkinter import font as tkfont
    from tkinter import Tk
    from tkinter import filedialog 
    from tkinter import filedialog as tkFileDialog
    from tkinter import messagebox as tkMessageBox

import argparse
import sys
from astropy.io import fits
from astropy.table import Column
from astropy.table import Table
import numpy as np
import subprocess
import os
import signal
from mypython.ifu import muse_utils as utl
import numpy as np


class TableEdit(tkinter.Frame):
    def __init__(self, parent, rows=10, columns=2):
        # use black background so it "peeks through" to 
        # form grid lines
        tkinter.Frame.__init__(self, parent)
        self._table_labels = []
        self._table_cells = []
        self.clicked=None
        self.rows=rows
        self.columns=columns
        for row in range(self.rows):
            current_row_lab  = []
            current_row_cel  = []
            for column in range(self.columns):
                if column == self.columns-1:
                   width=24
                else:
                   width=8   
                label = tkinter.StringVar()
                cell  = tkinter.Entry(self,textvariable=label,width=width)
                cell.grid(row=row, column=column, sticky="nsew", padx=1, pady=1)
                cell.bind("<Button-1>",self.lastclicked)
                label.set("0.0")
                current_row_lab.append(label)
                current_row_cel.append(cell)
                
            self._table_labels.append(current_row_lab)
            self._table_cells.append(current_row_cel)

        for column in range(self.columns):
            self.grid_columnconfigure(column, weight=1)


    def set(self, row, column, value):
        widget = self._table_labels[row][column]
        widget.set(value)

    def get(self, row, column):
        widget = self._table_labels[row][column]
        return widget.get()

    def lock(self, row, column):
        widget = self._table_cells[row][column]
        widget.configure(state='readonly',readonlybackground='white')

    def unlock(self, row, column):
        widget = self._table_cells[row][column]
        widget.configure(state='normal')
    
    def hasfocus(self):
        for row in range(self.rows):
            for column in range(self.columns):
                if(self.clicked == self._table_cells[row][column]):
                    return [row,column]
        #return widget.cget('') 
        #return widget.select_present()
        
    def lastclicked(self,event):
        self.clicked=event.widget



class TableFix(tkinter.Frame):
    def __init__(self, parent, rows=10, columns=2):
        # use black background so it "peeks through" to 
        # form grid lines
        self.frame = tkinter.Frame.__init__(self, parent)
        self._table_labels = []
        self._table_cells = []
        for row in range(rows):
            current_row  = []
            for column in range(columns):
                if column == columns-1:
                   width=24
                else:
                   width=8   
                cell  = tkinter.Label(self,text='Text',width=width)
                cell.grid(row=row, column=column, sticky="nsew", padx=1, pady=1)
                current_row.append(cell)
                
            self._table_labels.append(current_row)

        for column in range(columns):
            self.grid_columnconfigure(column, weight=1)


    def set(self, row, column, value):
        widget = self._table_labels[row][column]
        widget.configure(text=value)



class Window(tkinter.Tk):

    def __init__(self,parent,startfile=None,white=None):
        
        #store file name
        self.startfile=startfile
        self.white=white

        #start tk
        self.tk=Tk()

        #change getboolean: in Py3 return "??" when using go to page option. With this is ok.
        orig_getbool = self.tk.getboolean
        def safe_getboolean(s):
            try:
                return orig_getbool(s)
            except (ValueError, tkinter.TclError):
                return False
        self.tk.getboolean = safe_getboolean
        
        
        #set min and preferred size of main gui
        self.minwinwidth=400
        self.minwinheight=300
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        self.preferwinwidth=int(screen_width*0.7)
        self.preferwinheight=int(screen_height*0.5)
        self.minsize(width=self.minwinwidth, height=self.minwinheight)
        self.geometry("{}x{}".format(self.preferwinwidth,self.preferwinheight))

        #tweak the aspect ratio of the menu and data gui
        self.menuaspect=[1,0.2]        #Ruari 24/05 fixes bug where different resolutions cause the menu to be cut off 
        self.dataaspect=[1,1-0.2]     #Ruari 24/05 fixes bug where different resolutions cause the menu to be cut off 
        self.dpi=80

        #list for background processes
        self.processes=[]

        #Fiddle with font
        #default_font = tkFont.nametofont("TkDefaultFont")
        #scalefont = int(screen_height/1080.0*14)
        #default_font.configure(size=scalefont)

        #init gui frame
        self.initialize()

    def initialize(self):
        """ This init the basic gui """ 
        
        #create a menu frame
        self.menuframe=tkinter.Frame(self,width=int(self.preferwinwidth*self.menuaspect[0]),
                                     height=int(self.preferwinheight*self.menuaspect[1]))
        self.menuframe.grid_propagate(0)
        self.menuframe.grid()

        #create a data frame
        self.dataframe=tkinter.Frame(self,width=int(self.preferwinwidth*self.dataaspect[0]), 
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
       

        # creating a quite command
        self.quitButton_w = tkinter.Button(self.menuframe, text="Quit",command=self.client_exit).grid(row=0,column=0)
        
        # creating an abort command
        self.abortButton_w = tkinter.Button(self.menuframe, text="Abort",command=self.client_abort).grid(row=1,column=0)


        # creating a save command
        self.saveButton_w = tkinter.Button(self.menuframe, text="Save",command=self.write_current).grid(row=0,column=1)

        #create the control page
        self.nextButton_w = tkinter.Button(self.menuframe, text=">>",command=self.next_page).grid(row=0,column=4)
        self.previousButton_w = tkinter.Button(self.menuframe, text="<<",command=self.previous_page).grid(row=0,column=2)

        self.currentpage=tkinter.StringVar()
        self.currentpage.set(1)
        self.allpages=tkinter.StringVar()
        self.allpages.set(1)
        self.statuspage=tkinter.StringVar()
        self.statuspage.set("Page {}/{}".format(self.currentpage.get(),self.allpages.get()))
        self.currentpage_w=tkinter.Label(self.menuframe,textvariable = self.statuspage)
        self.currentpage_w.grid(column=3,row=0)

        # creating a save command
        self.firstButton_w = tkinter.Button(self.menuframe, text="First",command=self.gotofirst).grid(row=0,column=5)
        self.lastButton_w = tkinter.Button(self.menuframe, text="Last",command=self.gotolast).grid(row=0,column=6)


        #create inspect option
        self.inspectButton_w = tkinter.Button(self.menuframe, text="Inspect Current",command=self.inspect_current).grid(row=1,column=1,columnspan=2)


        #set header table properties
        self.keycol=['id','x_geow','y_geow','lambda_fluxw','SNR','SNR_odd','SNR_even','SNR_med','EODeltaSN','BoxFraction','OverContinuum','relvel','inspect','confidence','redshift','CogFlux','CogErr','CogRad','type','notes']
        self.shortkey=['ID','X','Y','Lambda','SNR','SNRodd','SNReven','SNRmed','EODSN','BoxFrac','Continu','relvel','inspect','confid','redshift','CogFlux','CogErr','CogRad','type','notes']
        self.tabncol=len(self.keycol)

        #create sort by option
        llab = tkinter.Label(self.menuframe, text="Sort by:")
        llab.grid(column=3,row=1)
        self.sortlist = tkinter.StringVar(self.menuframe)
        self.sortlist.set("id") # default value
        self.sortlist_w = tkinter.OptionMenu(self.menuframe, self.sortlist,'id','optimal','SNR','confidence','redshift','type','x_geow','y_geow','lambda_fluxw','notes','relvel','inspect')
        self.sortlist_w.grid(column=4,row=1)
        #set the linelist in trace state
        self.sortlist.trace("w",self.sorttab)

        #create reverse option
        self.reverselist_w = tkinter.Button(self.menuframe,text="Reverse",command=self.reversetab)
        self.reverselist_w.grid(column=5,row=1)
        
        #------------------------------------------
        # added by D.T. 
        
        # create go to page option
        page_lab = tkinter.Label(self.menuframe, text="Go to page:")
        page_lab.grid(row=0, column=7)
        self.page_input = tkinter.StringVar()
        self.page_inputcntr = tkinter.Entry(self.menuframe,textvariable=self.page_input, width=6)
        self.page_inputcntr.grid(row=0, column=8)
        
        self.page_inputcntr.bind('<Return>', self.goto_page_command) 
        #------------------------------------------
        
        
    
        #redshift reference
        zlab = tkinter.Label(self.menuframe, text="z_ref = ")
        zlab.grid(column=9,row=0)
        self.relvel = tkinter.StringVar()
        self.relvelcntr = tkinter.Entry(self.menuframe,textvariable=self.relvel, width=6)
        self.relvelcntr.grid(column=10,row=0)
        self.relvel.set("0.0000")
        #set the redshift in a trace state
        self.relvel.trace("w",self.updaterelvel)

        
        #reference lambda
        llab = tkinter.Label(self.menuframe, text="lam_ref = ")
        llab.grid(column=7,row=1)
        self.rellam = tkinter.StringVar()
        self.rellamcntr = tkinter.Entry(self.menuframe,textvariable=self.rellam, width=6)
        self.rellamcntr.grid(column=8,row=1)
        self.rellam.set("1215.6701")
        #set the redshift in a trace state
        self.rellam.trace("w",self.updaterelvel)

        #start flux analysis GUI by D.T.
        #self.FluxAnalysisCurrentButton = tkinter.Button(self.menuframe, text="GUI Flux Analysis Current ",command=self.GUIFluxAnalysisCurrent).grid(row=1,column=9,columnspan=2)

    def init_dataframe(self):

        #open the file
        self.catdata=Table.read(self.startfile,format='fits')
        #save the time stamp
        self.timeio=os.path.getmtime(self.startfile)
        #check if lock file exists 
        if os.path.exists(self.startfile+'.lock'):
            tkMessageBox.showinfo('I/O Warning','This catalogue is currently locked by another user. Best not to continue.')
            self.client_abort(keeplock=True)
        else:
            #make one
            open(self.startfile+'.lock','a').close()

        if('redshift' not in self.catdata.colnames):
            newcl=Column(np.zeros(len(self.catdata)))
            self.catdata.add_column(newcl,name='redshift')
        if('CogFlux' not in self.catdata.colnames):
            newcl=Column(np.zeros(len(self.catdata)))
            self.catdata.add_column(newcl,name='CogFlux')
        if('CogErr' not in self.catdata.colnames):
            newcl=Column(np.zeros(len(self.catdata)))
            self.catdata.add_column(newcl,name='CogErr')
        if('CogRad' not in self.catdata.colnames):
            newcl=Column(np.zeros(len(self.catdata)))
            self.catdata.add_column(newcl,name='CogRad')
        if('inspect' not in self.catdata.colnames):
            newcl=Column(np.zeros(len(self.catdata)))
            self.catdata.add_column(newcl,name='inspect')
        if('type' not in self.catdata.colnames):
            newcl=Column(np.full(len(self.catdata),'None',dtype="S50"))
            self.catdata.add_column(newcl,name='type')
        if('notes' not in self.catdata.colnames):
            newcl=Column(np.full(len(self.catdata),'None',dtype="S50"))
            self.catdata.add_column(newcl,name='notes')
        if('relvel' not in self.catdata.colnames):
            newcl=Column(np.zeros(len(self.catdata)))
            self.catdata.add_column(newcl,name='relvel')
            self.updaterelvel('firsttime')


        #find number of pages
        self.rowppage=15
        self.pagenum=len(self.catdata) // self.rowppage
        if((len(self.catdata) % self.rowppage) > 0):
            self.pagenum=self.pagenum+1
        
        #update page
        self.allpages.set(self.pagenum)
        self.statuspage.set("Page {}/{}".format(self.currentpage.get(),self.allpages.get()))


        #add table
        self.table_header = TableFix(self.dataframe, 1,self.tabncol)
        self.table_header.pack(side="top", fill="x")
        for jj,txt in enumerate(self.shortkey):
            self.table_header.set(0,jj,txt)

        self.table_data = TableEdit(self.dataframe,self.rowppage,self.tabncol)
        self.table_data.pack(side="top", fill="x")
        
        #write table entries
        self.update_tabdisplay()

    def record_changes(self):
        
        #grab ids
        idred=self.keycol.index('redshift')
        idtype=self.keycol.index('type')
        idid=self.keycol.index('id')
        idconf=self.keycol.index('confidence')
        idnote=self.keycol.index('notes')
        idflux=self.keycol.index('CogFlux')
        idferr=self.keycol.index('CogErr')
        idfrad=self.keycol.index('CogRad')

        #scan them 
        for r in range(self.rowppage):
            myid=self.table_data.get(r,idid)
            myred=self.table_data.get(r,idred)
            mytype=self.table_data.get(r,idtype)
            myconf=self.table_data.get(r,idconf)
            mynote=self.table_data.get(r,idnote)
            myflux=self.table_data.get(r,idflux)
            myferr=self.table_data.get(r,idferr)
            myfrad=self.table_data.get(r,idfrad)
            #store in table by searching for id (allow for blanks)
            try:
                cidx=np.where(self.catdata['id'] == int(myid))
                self.catdata['redshift'][cidx]=myred
                self.catdata['type'][cidx]=mytype
                self.catdata['confidence'][cidx]=myconf
                self.catdata['notes'][cidx]=mynote
                self.catdata['CogFlux'][cidx]=myflux
                self.catdata['CogErr'][cidx]=myferr
                self.catdata['CogRad'][cidx]=myfrad
            except:
                pass

    def update_tabdisplay(self):
        
        #num columns
        ncol=len(self.keycol)
                
        #clear page
        nrow=self.rowppage
        for r in range(nrow):
            for c in range(ncol):
                self.table_data.set(r,c,'')
                
        #update current entries for this page
        self.current_starix=(int(self.currentpage.get())-1)*self.rowppage
        self.current_endix=self.current_starix+self.rowppage
        if(self.current_endix > len(self.catdata)):
            self.current_endix=len(self.catdata)
        nrow=self.current_endix-self.current_starix
        
        for r in range(nrow):
            for c in range(ncol):
                rowindx=self.current_starix+r
                self.table_data.set(r,c,self.catdata[self.keycol[c]][rowindx])

                #lock fields that do not need edits 
                if('redshift' in self.keycol[c]):
                    self.table_data.unlock(r,c)
                elif('type' in self.keycol[c]):
                    self.table_data.unlock(r,c)
                elif('confidence' in self.keycol[c]):
                    self.table_data.unlock(r,c)
                elif('CogFlux' in self.keycol[c]):
                    self.table_data.unlock(r,c)
                elif('CogErr' in self.keycol[c]):
                    self.table_data.unlock(r,c)
                elif('CogRad' in self.keycol[c]):
                    self.table_data.unlock(r,c)
                elif('notes' in self.keycol[c]):
                    self.table_data.unlock(r,c)
                else:
                    self.table_data.lock(r,c)

    def updaterelvel(self,*args):


        #compute relative velocity with respect to reference 
        currentz=float(self.relvel.get())
        currentl=float(self.rellam.get())

        zline=utl.airtovac(self.catdata['lambda_fluxw'])/currentl-1
        ratio=(1+currentz)/(1+zline)
        self.catdata['relvel']=np.abs(299792.458*(ratio**2-1)/(1+ratio**2))
        
        if 'firsttime' in args:
            pass
        else:
            #update entries if table has been drawn already
            self.record_changes()
            self.update_tabdisplay()
            

    

    def inspect_current(self,oldformat=False):
        
        #kill open processes
        for pp in self.processes:
            pp.kill()
        self.processes=[]

        #find the row where focus is
        row=self.table_data.hasfocus()[0]
        #query which ID it is
        idid=self.keycol.index('id')
        idred=self.keycol.index('redshift')
        idlambda=self.keycol.index('lambda_fluxw')
        focusid=self.table_data.get(row,idid)        
        #launch displays
        
        ds9_path = 'ds9'
        
        try:
            #control ds9
            if(oldformat):
                rtname='objs/id{}/Pstamp_id{}'.format(focusid,focusid)
            else:
                rtname  ='objs/id{}/id{}_img.fits'.format(focusid,focusid)
                hdulist = fits.open(rtname)
                nexten  = len(hdulist)
                
            if(self.white is not None):
                if(oldformat):
                    ds9=subprocess.Popen([ds9_path,'-scale','clezscale','-lock','smooth','-lock','frame','wcs',self.white,rtname+'_mean.fits',rtname+'_median.fits',rtname+'_half1.fits',rtname+'_half2.fits',rtname+'_det.fits','-smooth'])
                else:
                    if nexten > 7:
                        ds9=subprocess.Popen([ds9_path, '-scale','zscale','-lock','smooth','-lock','frame','wcs',self.white,rtname+'[8]',rtname+'[9]',rtname+'[10]',rtname+'[11]',rtname+'[7]','-smooth'])
                    else:
                        ds9=subprocess.Popen([ds9_path, '-scale','zscale','-lock','smooth','-lock','frame','wcs',self.white,rtname+'[3]',rtname+'[4]',rtname+'[5]','-smooth'])
                        
            else:
                if(oldformat):
                    ds9=subprocess.Popen([ds9_path, '-scale','zscale','-lock','smooth','-lock','frame','wcs',rtname+'_mean.fits',rtname+'_median.fits',rtname+'_half1.fits',rtname+'_half2.fits',rtname+'_det.fits','-smooth'])
                else:
                    if nexten > 7:
                        ds9=subprocess.Popen([ds9_path, '-scale','zscale','-lock','smooth','-lock','frame','wcs',rtname+'[8]',rtname+'[9]',rtname+'[10]',rtname+'[11]',rtname+'[7]','-smooth'])
                    else:
                        ds9=subprocess.Popen([ds9_path, '-scale','zscale','-lock','smooth','-lock','frame','wcs',rtname+'[3]',rtname+'[4]',rtname+'[5]','-smooth'])
                        
            #collect processes
            self.processes.append(ds9)
            
            #control spectra gui
            tmpz=float(self.table_data.get(row,idred))
            if(tmpz > 0):
                pass
            else:
                currentl=float(self.rellam.get())
                tmpz=float(self.table_data.get(row,idlambda))/currentl-1
                print(tmpz)
            spc=subprocess.Popen(['python','{}/redshifts/zfit.py'.format(os.environ['MYPYTHON']),'-i','objs/id{}/spectrum.fits'.format(focusid),'-z','{}'.format(tmpz)])
            self.processes.append(spc)
          
  
            #view cube
            if(oldformat):
                ds93d=subprocess.Popen([ds9_path,'-3d','objs/id{}/segcube.fits'.format(focusid),'-3d','vp','90','0','-cmap','color'])
            else:
                ds93d=subprocess.Popen([ds9_path,'-3d','objs/id{}/id{}_img.fits[1]'.format(focusid,focusid),'-3d','vp','90','0','-cmap','color'])
            self.processes.append(ds93d)

        except:
            pass



    def gotofirst(self):
        self.record_changes()
        self.currentpage.set(1)
        self.statuspage.set("Page {}/{}".format(self.currentpage.get(),self.allpages.get()))
        self.update_tabdisplay()

    def gotolast(self):
        self.record_changes()
        self.currentpage.set(int(self.allpages.get()))
        self.statuspage.set("Page {}/{}".format(self.currentpage.get(),self.allpages.get()))
        self.update_tabdisplay()
            

    def previous_page(self):
    
        if(int(self.currentpage.get()) > 1):
            #freeze current changes
            self.record_changes()
            #step page and update
            self.currentpage.set(int(self.currentpage.get())-1)
            self.statuspage.set("Page {}/{}".format(self.currentpage.get(),self.allpages.get()))
            self.update_tabdisplay()
        else:
            pass


    def next_page(self):
        
        if(int(self.currentpage.get()) < self.pagenum):
            #freeze current changes
            self.record_changes()
            #step page and update
            self.currentpage.set(int(self.currentpage.get())+1)
            self.statuspage.set("Page {}/{}".format(self.currentpage.get(),self.allpages.get()))
            self.update_tabdisplay()
        else:
            pass
        
    #added by D.T.: this is a fucntion that allows to go to a specific page if the number is provided   
    def goto_page(self, page_number):
        try:
            
            #Convert the input to an integer
            page_number = int(page_number)
            
            if 1 <= page_number <= int(self.allpages.get()):
                
                #freeze current changes
                self.record_changes()
                #step page and update
                self.currentpage.set(page_number)
                self.statuspage.set("Page {}/{}".format(self.currentpage.get(), self.allpages.get()))
                self.update_tabdisplay()
                
            else:
                print("Invalid page number. Please enter a number between 1 and {}.".format(self.allpages.get()))
                pass
        
        except ValueError:
            print("Invalid input. Please enter a valid page number.")
            pass
    
    
    def goto_page_command(self, event):
        #the event is the <Return> I give from the "enter" on the keyboard
        #this function define the command to go to a certain page
        page_number = self.page_input.get()
        self.goto_page(page_number)
     # -- -- -- -- -- -- -- -- -- -- -- --  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
        
            
    def sorttab(self,*args):
    
        #record changes
        self.record_changes()
                
        #sort table
        if('optimal' in self.sortlist.get()):
            self.catdata.sort(['inspect','OverContinuum','SNR'])
            self.catdata.reverse()
        else:
            self.catdata.sort(self.sortlist.get())
        
        #reset page number and update display 
        self.currentpage.set(1)
        self.statuspage.set("Page {}/{}".format(self.currentpage.get(),self.allpages.get()))
        self.update_tabdisplay()
        
    def reversetab(self):
        
        #record changes
        self.record_changes()

        #reverse and update
        self.catdata.reverse()

        self.currentpage.set(1)
        self.statuspage.set("Page {}/{}".format(self.currentpage.get(),self.allpages.get()))
        self.update_tabdisplay()

    def client_exit(self):
        
        #record changes
        self.record_changes()

        #trigger write current on exit
        self.write_current()

        #exit
        self.client_abort()


    def client_abort(self,keeplock=False):
        
        """
        Kill without saving
        """

        #kill shell processes
        for pp in self.processes:
            pp.kill()
        
            
        #clean lock file
        if(keeplock):
            pass
        else:
            os.remove(self.startfile+'.lock')

        #exit for good
        self.destroy()


    def write_current(self):
        
        #record changes
        self.record_changes()

        #sort again table in id -- make copy not to interfere with current sorting 
        tmpcopy=self.catdata.copy()
        tmpcopy.sort('id')
        
        #strip relvel column
        tmpcopy.remove_column('relvel')

        #check time stamp to avoid overwite of catalogue that has changed to disk
        currenttime=os.path.getmtime(self.startfile)
        if(currenttime > self.timeio):
            #possible issues due to io
            tkMessageBox.showinfo('I/O Warning','The catalogue may have changed to disk. Writing conflicting copy to disk. Resolve issue before proceding')
            tmpcopy.write(self.startfile+'_conflicted',format='fits')
            self.client_abort()
        else:
            #write to disk overwriting
            tmpcopy.write(self.startfile,format='fits',overwrite=True)
            #reset clock
            self.timeio=os.path.getmtime(self.startfile)

def start(startfile,white):

    app = Window(None,startfile=startfile,white=white)
    app.title('Emitter classifier')
    app.mainloop()

    
if __name__=="__main__":

    
    parser = argparse.ArgumentParser(description='Process input')
    parser.add_argument('ifile',default=None,help='Input file')
    parser.add_argument('-w','--white',default=None,help='White image')

    args = parser.parse_args()
    start(args.ifile,args.white)





