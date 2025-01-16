import tkinter as tk
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from linetools.lists.linelist import LineList        


class SpectrumPlotter:
    def __init__(self, master):
        self.master = master
        self.master.title("Spectrum Plotter")

        # Create the figure and canvas for matplotlib
        self.figure, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.figure, master)
        self.canvas_widget = self.canvas.get_tk_widget()
    
        # Pack the canvas widget
        self.canvas_widget.pack()

        # Add the navigation toolbar
        self.toolbar = NavigationToolbar2Tk(self.canvas, master)
        self.toolbar.update()
        self.toolbar.pack()

        # Create buttons and entries
        self.load_button = tk.Button(master, text="Load Spectrum", command=self.load_spectrum)
        self.redshift_label = tk.Label(master, text="Redshift:")
        textred = tk.StringVar() 
        textred.set("0.0") 
        self.redshift_entry = tk.Entry(master,textvariable = textred)
        self.lines_type_label=tk.Label(master, text="Line type:")
        textline = tk.StringVar() 
        textline.set("Strong") 
        self.lines_type_entry=tk.Entry(master,textvariable = textline)
        self.plot_button = tk.Button(master, text="Plot lines", command=self.plot_lines)

        # Pack the buttons and entries
        self.load_button.pack()
        self.redshift_label.pack()
        self.redshift_entry.pack()
        self.lines_type_label.pack()
        self.lines_type_entry.pack()
        self.plot_button.pack()

        self.spectrum_data = None
        self.lines_data = None

    def load_spectrum(self):
        filepath = filedialog.askopenfilename(filetypes=[("ASCII files", "*.txt *.dat *fits")])
        if "fits" in filepath:
            self.spectrum_data = Table.read(filepath)

            self.wave=self.spectrum_data['OPT_WAVE']
            self.flux=self.spectrum_data['OPT_FLAM']
            self.err=self.spectrum_data['OPT_FLAM_SIG']
            self.plot_spectrum(self.wave,self.flux,self.err)
        else:
            self.spectrum_data = np.loadtxt(filepath)
            self.wave,self.flux,self.err = self.spectrum_data.T 
            self.plot_spectrum(self.wave,self.flux,self.err)


    def plot_lines(self):
        
        self.lines_data = LineList(self.lines_type_entry.get())
        self.plot_spectrum(self.wave,self.flux,self.err, lines=self.lines_data)

            
    def plot_spectrum(self,wavelength, flux, error, lines=None):
                
        if self.spectrum_data is not None:

            # Clear the current plot
            self.ax.clear()

            # Plot the flux and error
            self.ax.plot(wavelength, flux, label='Flux')
            self.ax.fill_between(wavelength, flux - error, flux + error, alpha=0.2, label='Error')

            # Plot vertical lines for redshift
            redshift = float(self.redshift_entry.get() or 0)
            if lines is not None:
                lines=self.lines_data.available_transitions((min(self.wave/(1+redshift)),max(self.wave/(1+redshift)))*u.AA)
                for one_line in lines:
                    shifted_wavelength = one_line['wrest'] * (1 + redshift)
                    self.ax.axvline(x=shifted_wavelength, color='r', linestyle='--')
                    self.ax.text(shifted_wavelength,0,one_line['ion_name'],rotation=90)
                    
            self.ax.set_xlabel('Wavelength')
            self.ax.set_ylabel('Flux')
            self.ax.legend()

            # Redraw the canvas
            self.canvas.draw()

    def load_lines(self, filepath):
        self.lines_data = np.loadtxt(filepath)

# Create the main window
root = tk.Tk()
app = SpectrumPlotter(root)
root.mainloop()
