"""

GUI to handle visual inspection of emitters from MUSE
Assumes file written with the muse_emitters.py procedures

"""

from tkinter import *
import Tkinter
import tkFont
from Tkinter import Tk
import tkFileDialog
import argparse
import sys
from astropy.io import fits
from astropy.table import Column
from astropy.table import Table
import numpy as np
import subprocess
import os
import signal

class TableEdit(Tkinter.Frame):
    def __init__(self, parent, rows=10, columns=2):
        # use black background so it "peeks through" to 
        # form grid lines
        Tkinter.Frame.__init__(self, parent)
        self._table_labels = []
        self._table_cells = []
        self.clicked=None
        self.rows=rows
        self.columns=columns
        for row in range(self.rows):
            current_row_lab  = []
            current_row_cel  = []
            for column in range(self.columns):
                label = Tkinter.StringVar()
                cell  = Tkinter.Entry(self,textvariable=label,width=8)
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
        widget.configure(state=NORMAL)
    
    def hasfocus(self):
        for row in range(self.rows):
            for column in range(self.columns):
                if(self.clicked == self._table_cells[row][column]):
                    return [row,column]
        #return widget.cget('') 
        #return widget.select_present()
        
    def lastclicked(self,event):
        self.clicked=event.widget



class TableFix(Tkinter.Frame):
    def __init__(self, parent, rows=10, columns=2):
        # use black background so it "peeks through" to 
        # form grid lines
        self.frame = Tkinter.Frame.__init__(self, parent)
        self._table_labels = []
        self._table_cells = []
        for row in range(rows):
            current_row  = []
            for column in range(columns):
                cell  = Tkinter.Label(self,text='Text',width=8)
                cell.grid(row=row, column=column, sticky="nsew", padx=1, pady=1)
                current_row.append(cell)
                
            self._table_labels.append(current_row)

        for column in range(columns):
            self.grid_columnconfigure(column, weight=1)


    def set(self, row, column, value):
        widget = self._table_labels[row][column]
        widget.configure(text=value)



class Window(Tkinter.Tk):

    def __init__(self,parent,startfile=None):
        
        #store file name
        self.startfile=startfile

        #start tk
        self.tk=Tk()
        
        #set min and preferred size of main gui
        self.minwinwidth=300
        self.minwinheight=300
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        self.preferwinwidth=int(screen_width*0.5)
        self.preferwinheight=int(screen_height*0.5)
        self.minsize(width=self.minwinwidth, height=self.minwinheight)
        self.geometry("{}x{}".format(self.preferwinwidth,self.preferwinheight))

        #tweak the aspect ratio of the menu and data gui
        self.menuaspect=[1,0.2]     #Ruari 24/05 fixes bug where different resolutions cause the menu to be cut off 
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
       

        # creating a quite command
        self.quitButton_w = Tkinter.Button(self.menuframe, text="Quit",command=self.client_exit).grid(row=0,column=0)
        
        # creating a save command
        self.saveButton_w = Tkinter.Button(self.menuframe, text="Save",command=self.write_current).grid(row=0,column=1)

        #create the control page
        self.nextButton_w = Tkinter.Button(self.menuframe, text=">>",command=self.next_page).grid(row=0,column=4)
        self.previousButton_w = Tkinter.Button(self.menuframe, text="<<",command=self.previous_page).grid(row=0,column=2)

        self.currentpage=Tkinter.StringVar()
        self.currentpage.set(1)
        self.allpages=Tkinter.StringVar()
        self.allpages.set(1)
        self.statuspage=Tkinter.StringVar()
        self.statuspage.set("Page {}/{}".format(self.currentpage.get(),self.allpages.get()))
        self.currentpage_w=Tkinter.Label(self.menuframe,textvariable = self.statuspage)
        self.currentpage_w.grid(column=3,row=0)

        # creating a save command
        self.firstButton_w = Tkinter.Button(self.menuframe, text="First",command=self.gotofirst).grid(row=0,column=4)
        self.lastButton_w = Tkinter.Button(self.menuframe, text="Last",command=self.gotolast).grid(row=0,column=5)


        #create inspect option
        self.inspectButton_w = Tkinter.Button(self.menuframe, text="Inspect Current",command=self.inspect_current).grid(row=1,column=0,columnspan=2)



        #set header table properties
        self.keycol=['id','x_geow','y_geow','lambda_fluxw','SNR','SNR_odd','SNR_even','SNR_med','confidence','EODeltaSN','BoxFraction','OverContinuum','redshift','type']
        self.shortkey=['ID','X','Y','Lambda','SNR','SNRodd','SNReven','SNRmed','Confid','EODSN','BoxFrac','Continu','redshift','type']
        self.tabncol=len(self.keycol)

        #create sort by option
        llab = Tkinter.Label(self.menuframe, text="Sort by:")
        llab.grid(column=2,row=1)
        self.sortlist = Tkinter.StringVar(self.menuframe)
        self.sortlist.set("id") # default value
        self.sortlist_w = Tkinter.OptionMenu(self.menuframe, self.sortlist,'id','SNR','redshift','type')
        self.sortlist_w.grid(column=3,row=1)
        #set the linelist in trace state
        self.sortlist.trace("w",self.sorttab)

        #create reverse option
        self.reverselist_w = Tkinter.Button(self.menuframe,text="Reverse",command=self.reversetab)
        self.reverselist_w.grid(column=4,row=1)
    

    def init_dataframe(self):

        #open the file
        self.catdata=Table.read(self.startfile,format='fits')
        
        if('redshift' not in self.catdata.colnames):
            newcl=Column(np.zeros(len(self.catdata)))
            self.catdata.add_column(newcl,name='redshift')
        if('type' not in self.catdata.colnames):
            newcl=Column(np.full(len(self.catdata),'None',dtype="S25"))
            self.catdata.add_column(newcl,name='type')
            
        #find number of pages
        self.rowppage=15
        self.pagenum=len(self.catdata)/self.rowppage
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

        #scan them 
        for r in range(self.rowppage):
            myid=self.table_data.get(r,idid)
            myred=self.table_data.get(r,idred)
            mytype=self.table_data.get(r,idtype)
            #store in table by searching for id (allow for blanks)
            try:
                cidx=np.where(self.catdata['id'] == int(myid))
                self.catdata['redshift'][cidx]=myred
                self.catdata['type'][cidx]=mytype
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
                else:
                    self.table_data.lock(r,c)

    def inspect_current(self):
        
        #kill open processes
        for pp in self.processes:
            pp.kill()
        self.processes=[]

        #find the row where focus is
        row=self.table_data.hasfocus()[0]
        #query which ID it is
        idid=self.keycol.index('id')
        idred=self.keycol.index('redshift')
        focusid=self.table_data.get(row,idid)        
        #launch displays
        try:
            ds9=subprocess.Popen(['ds9','-scale','zscale','-lock','smooth','-lock','frame','wcs','objs/id{}/Pstamp_id{}_mean.fits'.format(focusid,focusid),'-smooth'])
        
            spc=subprocess.Popen(['python','{}/redshifts/zfit.py'.format(os.environ['MYPYTHON']),'-i','objs/id{}/spectrum.fits'.format(focusid),'-z','{}'.format(self.table_data.get(row,idred))])
            
            #collect processes
            self.processes.append(ds9)
            self.processes.append(spc)

        except:
            pass

#ds9 mean_detection.Objects_Id.fits &

#python $MYPYTHON/redshifts/zfit.py -i id$1/spectrum.fits -z 4.9

    


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
            
    def sorttab(self,*args):
    
        #record changes
        self.record_changes()
                
        #sort table
        self.catdata.sort(self.sortlist.get())
        
        #reset page number and update display 
        self.currentpage.set(1)
        self.update_tabdisplay()
        
    def reversetab(self):
        
        #record changes
        self.record_changes()

        #reverse and update
        self.catdata.reverse()
        self.update_tabdisplay()

    def client_exit(self):
        
        #record changes
        self.record_changes()

        #trigger write current on exit
        self.write_current()

        #kill shell processes
        for pp in self.processes:
            pp.kill()
        

        #exit for good
        exit()

    def write_current(self):
        
        #record changes
        self.record_changes()

        #sort again table in id 
        self.catdata.sort('id')
        #write to disk 
        self.catdata.write(self.startfile,format='fits',overwrite=True)

def start(startfile):

    app = Window(None,startfile=startfile)
    app.title('Emitter classifier')
    app.mainloop()

    
if __name__=="__main__":

    
    parser = argparse.ArgumentParser(description='Process input')
    parser.add_argument('ifile',default=None,help='Input file')
    args = parser.parse_args()
    start(args.ifile)





