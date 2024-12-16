"""

Script that handles a reduced MUSE cube and extracts subcubes to be handled by Ai classifiers

INPUT:
- A parameter file with the following structure:

[Paths]
pathcat = ../emitter_sources_358/extracted_SN7/
catalogue = COMBINED_CUBE_FINAL_bootvar_psfsub_bkgsub_select_SNR.fits
pathcube = ../emitter_sources_358/
cube = COMBINED_CUBE_FINAL_bootvar_psfsub_bkgsub.fits

[Wavelength]
minw = 4850
maxw = 9500

"""

from astropy.table import Table
import astropy.units as u
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import configparser
import argparse


parser = argparse.ArgumentParser(description='Process a parameter file.')
parser.add_argument('param_file', type=str, help='Path to the parameter file')
args = parser.parse_args()

# Create a ConfigParser object
config = configparser.ConfigParser()

# Read the parameter file
config.read(args.param_file)

# Access the parameters
pathcat = config['Paths']['pathcat']
catalogue = config['Paths']['catalogue']
pathcube = config['Paths']['pathcube']
cube = config['Paths']['cube']
minw = int(config['Wavelength']['minw'])
maxw = int(config['Wavelength']['maxw'])
augment = bool(config['Data']['augment'])
lineonly = bool(config['Data']['lineonly'])

#pathcat = '../emitter_sources_358/extracted_SN7/'
#catalogue='COMBINED_CUBE_FINAL_bootvar_psfsub_bkgsub_select_SNR.fits'
#pathcube= '../emitter_sources_358/'
#cube='COMBINED_CUBE_FINAL_bootvar_psfsub_bkgsub.fits'
#minw=4850
#maxw=9500

#open fits
datacube=fits.open(pathcube+cube)[0].data
header=fits.open(pathcube+cube)[0].header

#reconstruct wave and trim in defined range
#wave in air
sz=datacube.shape
delta_lambda=header["CD3_3"]
zero_lambda=header["CRVAL3"]
wavec=np.arange(0,sz[0],1)*delta_lambda+zero_lambda


#trim between 4850 and 9500
indexw=np.where((wavec>=minw) & (wavec<=maxw))

#open LAE
laetab=Table.read(pathcat+catalogue)
nlae=len(laetab)
laetab.add_column(np.zeros(nlae,dtype='int64')+np.min(indexw),name='trim_minindex')
laetab.add_column(np.zeros(nlae,dtype='int64')+np.max(indexw),name='trim_maxindex')

for lae in laetab:

    #select in range
    if((lae['lambda_geow']>minw) & (lae['lambda_geow']<maxw)):

        #round
        yc=np.int64(lae['x_geow'])
        xc=np.int64(lae['y_geow'])
        zc=np.int64(lae['z_geow'])
        
        #try:
        #size_pix=25 i.e., 5 arcsec on a side
        subcube=datacube[:,xc-12:xc+13,yc-12:yc+13]

        #find a 41 channel centered on emitter and copy as head part of the data
        zoomcube=subcube[zc-20:zc+21,:,:]
        if lineonly:
            finalcut=zoomcube
        else:
            finalcut=np.concatenate((zoomcube,subcube[indexw[0],:,:]))

        #print(zoomcube.shape)
        #print(subcube[indexw].shape)
        #print(subcube.shape)
        #print(finalcut.shape)
        #exit()
        
        #remove laser gap
        #subcube=np.concatenate((subcubeall[0:925,:,:],subcubeall[1126:-1,:,:]),axis=0)

        #img=np.nansum(subcube,axis=0)
        #plt.imshow(img)
        #plt.show()
        outfile = 'data/subcube_{}.fits'.format(lae['id'])
        hduc = fits.PrimaryHDU(data=finalcut)
        hdut = fits.BinTableHDU(data=Table(lae))
        hdu = fits.HDUList([hduc, hdut])
        hdu.writeto(outfile, overwrite=True)
        

        #augment by rotating if needed
        if((augment) & (('LAE' in lae['type']) | ('lae' in lae['type']))):
            for i in range(1,4):
                #rotate
                finalcut=np.rot90(finalcut,axes=(1,2))
                outfile = 'data/subcube_{}_{}.fits'.format(lae['id'],i)
                hduc = fits.PrimaryHDU(data=finalcut)
                hdut = fits.BinTableHDU(data=Table(lae))
                hdu = fits.HDUList([hduc, hdut])
                hdu.writeto(outfile, overwrite=True)


        #except:
        #    #drop sources close to edges that give size problem
        #    pass
        
            
