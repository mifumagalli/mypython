import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import font as tkfont
from PIL import Image, ImageTk
import numpy as np
import sys
import subprocess
import matplotlib.pyplot as plt
import argparse

from mypython.ifu import muse_emitters as mf



"""
######################################################################################################################################################################

This is the GUI by D.T. that allow you to estimate the flux of an emitter, source by source, using the Curve of Growth (Cog) method (see Fossati et al 2021) and varying all the parameters to obtain the best estimate.

Input:

config.txt: in the same directory of the muse_gui_cog.py (default) or in a separate directory that enters in --config option (NB: must to be named "config").
The configuration file must be a .txt file with (in the following order!!):
    * path of the cube (.fits) to obtain the NB image (psf-bkg subtract);
    * path of the variance cube (.fits) associated to the cube;
    * path of the segmentation map (.fits) with all the detected sources in the cube;
    * path of the catalogue (.fits) of sources with the properties (especially x,y,z);
    * path of the list of ids (.txt) of the sources we want to inspect;
    * path of the output directory (save the flux parameters for each source);

- It uses the functions to calculate the flux (emi_cogphot and nb_cogphot): import from mypython.ifu library in muse_emitters as mf.


Output:

- in the terminal: the flux, flux error and radius of the measured photometry                          
- a .txt file *for each source* with the same results and the parameters used to obtain the results in the output directory; 
                                                                                                       
######################################################################################################################################################################
"""

def parse_args():
    parser = argparse.ArgumentParser(description='FluxGUI Arguments')
    parser.add_argument('--config', type=str, help='Path to the configuration file (.txt)')
    return parser.parse_args()



class FluxGUI:

    #class builder
    def __init__(self, root, args=None):
        self.root = root
        self.root.title("Flux Analysis GUI")
        #self.root.geometry("800x600")
        self.root.resizable(False, False)

        png_icon_path = 'icon_gui.png'

        img = Image.open(png_icon_path)
        img = img.resize((128, 128))  
        self.icon = ImageTk.PhotoImage(img)
        self.root.iconphoto(False, self.icon)
        

        #Check message when you close the window
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

        # Initial parameters
        self.pathCube = ""
        self.pathVar = ""
        self.pathSeg = ""
        self.catalog_path = ""
        self.idlist = []
        self.outdir = ""
        self.num = 0
        self.dz = tk.DoubleVar(value=24)
        self.maxrad = tk.DoubleVar(value=15)
        self.offx = tk.DoubleVar(value=0)
        self.offy = tk.DoubleVar(value=0)
        self.offz = tk.DoubleVar(value=0)
        self.growthlim = tk.DoubleVar(value=1.025)

        # UI Elements
        self.create_widgets()

        # Initial configuration
        self.load_initial_configuration(args)

    
    def load_initial_configuration(self, args):

        if not args.config:  
            #read the configuration file
            config_file = 'config.txt'  # in the same directory of the GUI
        else:
            config_file = f'{args.config}/config.txt'

        try:
            with open(config_file, 'r') as f:
                lines = f.readlines()
                if len(lines) >= 6:
                    self.pathCube = lines[0].strip()
                    self.pathVar = lines[1].strip()
                    self.pathSeg = lines[2].strip()
                    self.catalog_path = lines[3].strip()
                    idlist_file_path = lines[4].strip()
                    self.outdir      = lines[5].strip()

                    # Load the id list from the file
                    try:
                        with open(idlist_file_path, 'r') as id_file:
                            id_lines = id_file.readlines()
                            self.idlist = [int(id.strip()) for id in id_lines if id.strip().isdigit()]

                        # Assign the entries to the GUI (I want to see what I'm using in the GUI itself)
                        self.cube_entry.delete(0, tk.END)
                        self.cube_entry.insert(0, self.pathCube)
                        self.var_entry.delete(0, tk.END)
                        self.var_entry.insert(0, self.pathVar)
                        self.seg_entry.delete(0, tk.END)
                        self.seg_entry.insert(0, self.pathSeg)
                        self.catalog_entry.delete(0, tk.END)
                        self.catalog_entry.insert(0, self.catalog_path)
                        self.id_entry.delete(0, tk.END)
                        self.id_entry.insert(0, ','.join(map(str, self.idlist)))
                        self.outdir_entry.delete(0, tk.END)
                        self.outdir_entry.insert(0, self.outdir)

                        #set the current id
                        if self.idlist:
                            self.current_id_label.config(text=str(self.idlist[self.num]))
                        else:
                            self.current_id_label.config(text="None")
                            
                    except Exception as e:
                        messagebox.showerror("Error", f"Failed to load ID list: {str(e)}")
                        
                else:
                    messagebox.showerror("Error", "Configuration file format is incorrect.")
                    
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load configuration: {str(e)}")

    def create_widgets(self):
        # Font configurations
        label_font = tkfont.Font(family="Helvetica", size=14)
        button_font = tkfont.Font(family="Helvetica", size=12)
        entry_font = tkfont.Font(family="Helvetica", size=14)  # Font per le Entry
        status_font = tkfont.Font(family="Helvetica", size=14, weight="bold")

        header_frame = tk.Frame(self.root)
        header_frame.grid(row=0, column=0, sticky="ew", pady=20)

        #Label for the icon
        logo_label = tk.Label(header_frame, image=self.icon)
        logo_label.grid(row=0, column=0, sticky="w")

        #Label for the title
        title_label = tk.Label(header_frame, text="Flux Analysis GUI", font=("Helvetica", 24, "bold"))
        title_label.grid(row=0, column=1, sticky="w")

        # Labels and Inputs (see the configuration)
        tk.Label(self.root, text="Path to Cube:", font=label_font).grid(row=1, column=0, sticky=tk.W)
        self.cube_entry = tk.Entry(self.root, width=80, font=entry_font)
        self.cube_entry.grid(row=1, column=1)
        #tk.Button(self.root, text="Browse", command=self.load_cube, font=button_font).grid(row=1, column=2)

        tk.Label(self.root, text="Path to Variance:", font=label_font).grid(row=2, column=0, sticky=tk.W)
        self.var_entry = tk.Entry(self.root, width=80, font=entry_font)
        self.var_entry.grid(row=2, column=1)
        #tk.Button(self.root, text="Browse", command=self.load_var, font=button_font).grid(row=2, column=2)

        tk.Label(self.root, text="Path to Segmentation:", font=label_font).grid(row=3, column=0, sticky=tk.W)
        self.seg_entry = tk.Entry(self.root, width=80, font=entry_font)
        self.seg_entry.grid(row=3, column=1)
        #tk.Button(self.root, text="Browse", command=self.load_seg, font=button_font).grid(row=3, column=2)

        tk.Label(self.root, text="Catalog Path:", font=label_font).grid(row=4, column=0, sticky=tk.W)
        self.catalog_entry = tk.Entry(self.root, width=80, font=entry_font)
        self.catalog_entry.grid(row=4, column=1)
        #tk.Button(self.root, text="Browse", command=self.load_catalog, font=button_font).grid(row=4, column=2)

        tk.Label(self.root, text="Outdir Path:", font=label_font).grid(row=5, column=0, sticky=tk.W)
        self.outdir_entry = tk.Entry(self.root, width=80, font=entry_font)
        self.outdir_entry.grid(row=5, column=1)

        tk.Label(self.root, text="ID List:", font=label_font).grid(row=6, column=0, sticky=tk.W)
        self.id_entry = tk.Entry(self.root, width=80, font=entry_font)
        self.id_entry.grid(row=6, column=1)

        # Parameters (see the parameters)
        tk.Label(self.root, text='Go to ID: ', font=label_font).grid(row=7, column=0, sticky=tk.W)
        self.go_to_id_entry = tk.Entry(self.root, font=entry_font)
        self.go_to_id_entry.grid(row=7, column=1)

        # Go to ID Button
        tk.Button(self.root, text="Go", command=self.go_to_id, font=button_font).grid(row=7, column=2)

        # ID Current Display
        tk.Label(self.root, text="Current ID:", font=label_font).grid(row=7, column=3, sticky=tk.W)
        self.current_id_label = tk.Label(self.root, text="None", font=label_font)
        self.current_id_label.grid(row=7, column=4)

        tk.Label(self.root, text="dz:", font=label_font).grid(row=8, column=0, sticky=tk.W)
        tk.Entry(self.root, textvariable=self.dz, font=entry_font).grid(row=8, column=1)

        tk.Label(self.root, text="maxrad:", font=label_font).grid(row=9, column=0, sticky=tk.W)
        tk.Entry(self.root, textvariable=self.maxrad, font=entry_font).grid(row=9, column=1)

        tk.Label(self.root, text="offx:", font=label_font).grid(row=10, column=0, sticky=tk.W)
        tk.Entry(self.root, textvariable=self.offx, font=entry_font).grid(row=10, column=1)

        tk.Label(self.root, text="offy:", font=label_font).grid(row=11, column=0, sticky=tk.W)
        tk.Entry(self.root, textvariable=self.offy, font=entry_font).grid(row=11, column=1)

        tk.Label(self.root, text="offz:", font=label_font).grid(row=12, column=0, sticky=tk.W)
        tk.Entry(self.root, textvariable=self.offz, font=entry_font).grid(row=12, column=1)

        tk.Label(self.root, text="growthlim:", font=label_font).grid(row=13, column=0, sticky=tk.W)
        tk.Entry(self.root, textvariable=self.growthlim, font=entry_font).grid(row=13, column=1)

        # Empty row
        tk.Label(self.root, text="").grid(row=14, column=0, pady=10)

        # Buttons
        tk.Button(self.root, text="Run Analysis", command=self.run_analysis, font=button_font).grid(row=15, column=1, columnspan=1)
        tk.Button(self.root, text="Next ID", command=self.next_id, font=button_font).grid(row=15, column=2, columnspan=1)
        tk.Button(self.root, text="Previous ID", command=self.previous_id, font=button_font).grid(row=15, column=0, columnspan=1)

        # Empty row
        tk.Label(self.root, text="").grid(row=16, column=0, pady=10)

        # Status
        self.status = tk.Label(self.root, text="", fg="red", font=status_font)
        self.status.grid(row=17, column=0, columnspan=3)
        

    def load_cube(self):
        self.pathCube = filedialog.askopenfilename(title="Select Cube File")
        self.cube_entry.delete(0, tk.END)
        self.cube_entry.insert(0, self.pathCube)

    def load_var(self):
        self.pathVar = filedialog.askopenfilename(title="Select Variance File")
        self.var_entry.delete(0, tk.END)
        self.var_entry.insert(0, self.pathVar)

    def load_seg(self):
        self.pathSeg = filedialog.askopenfilename(title="Select Segmentation File")
        self.seg_entry.delete(0, tk.END)
        self.seg_entry.insert(0, self.pathSeg)

    def load_catalog(self):
        self.catalog_path = filedialog.askopenfilename(title="Select Catalog File")
        self.catalog_entry.delete(0, tk.END)
        self.catalog_entry.insert(0, self.catalog_path)

    #defin the main program
    def run_analysis(self):

        #close the plot opened
        plt.close('all')
        
        if not self.pathCube or not self.pathVar or not self.pathSeg or not self.catalog_path:
            self.status.config(text="Error: All paths must be set.", fg='red')
            return

        if not self.idlist:
            self.status.config(text="Error: ID list is empty.", fg='red')
            return

        try:
            num       = self.num
            dz        = self.dz.get()
            maxrad    = self.maxrad.get()
            offx      = self.offx.get()
            offy      = self.offy.get()
            offz      = self.offz.get()
            growthlim = self.growthlim.get()

            subprocess.run('clear', shell=True)
            
            # Call the flux calculation function
            fluxarr, errarr, radarr = mf.emi_cogphot(
                self.pathCube, self.pathVar, self.pathSeg, self.catalog_path, [self.idlist[num]],
                dz=dz, maxrad=maxrad, offx=offx, offy=offy, offz=offz, growthlim=growthlim, writeNB=False, plots=True
            )
            
            # Save the result
            header = f'dz={dz}, maxrad={maxrad}, offx={offx}, offy={offy}, offz={offz}, growthlim={growthlim}, \n flux, flux_err, radius'

            np.savetxt(self.outdir+f'Flux_source_{self.idlist[num]}.txt',
                       np.column_stack((fluxarr, errarr, radarr)),
                       header=header, comments='')

            
            self.status.config(text="Analysis complete. Results saved.", fg='limegreen')
        except Exception as e:
            print(e)
            self.status.config(text=f"Error during analysis: {str(e)}", fg='red')

    def next_id(self):
        if self.idlist:
            self.num = (self.num + 1) % len(self.idlist)
            self.id_entry.delete(0, tk.END)
            self.id_entry.insert(0, ','.join(map(str, self.idlist)))
            
            # Update the label with the current id
            current_id = self.idlist[self.num]
            self.current_id_label.config(text=str(current_id))

            # Set again the default parameters
            self.dz.set(24)
            self.maxrad.set(15)
            self.offx.set(0)
            self.offy.set(0)
            self.growthlim.set(1.025)

            self.status.config(text="Ready for new ID. Click 'Run Analysis' to proceed.", fg='limegreen')
        else:
            self.status.config(text="Error: ID list is empty.", fg='red')
    
    def go_to_id(self):
        try:
            id_to_go = int(self.go_to_id_entry.get().strip())
            if id_to_go in self.idlist:
                self.num = self.idlist.index(id_to_go)
                self.current_id_label.config(text=str(id_to_go))
    
                # Set default parameters
                self.dz.set(24)
                self.maxrad.set(15)
                self.offx.set(0)
                self.offy.set(0)
                self.offz.set(0)
                self.growthlim.set(1.025)
                self.go_to_id_entry.delete(0, tk.END)
                
                self.status.config(text=f"Jumped to ID {id_to_go}. Click 'Run Analysis' to proceed.", fg='limegreen')
            else:
                self.status.config(text=f"ID {id_to_go} not found in the list.", fg='red')
        except ValueError:
            self.status.config(text="Invalid ID. Please enter a valid integer ID.", fg='red')
            
    
    def previous_id(self):
        if not self.idlist:
            self.status.config(text="Error: ID list is empty.", fg='red')
            return

        if self.num > 0:
            self.num -= 1
            self.current_id_label.config(text=str(self.idlist[self.num]))
            
            # Update the label with the current id
            current_id = self.idlist[self.num]
            self.current_id_label.config(text=str(current_id))

            # Set again the default parameters
            self.dz.set(24)
            self.maxrad.set(15)
            self.offx.set(0)
            self.offy.set(0)
            self.offz.set(0)
            self.growthlim.set(1.025)
            self.status.config(text="Ready for new ID. Click 'Run Analysis' to proceed.", fg='limegreen')
        else:
            self.status.config(text="No previous ID.", fg='orange')

    
    def on_closing(self):
        if messagebox.askokcancel("Quit", "Are you sure you want to close the GUI?"):
            self.root.destroy()

if __name__ == "__main__":
    args = parse_args()
    root = tk.Tk()                   #create the main GUI window
    app = FluxGUI(root, args=args)   #create an instance of the GUI class, initializing the elements and the GUI logic 
    root.mainloop()                  #enters the loop of tkinter and wait for commands, update the GUI

