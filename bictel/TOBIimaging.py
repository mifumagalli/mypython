import glob
import astropy
from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal as signal 
import os

from photutils.centroids import centroid_com
from photutils.aperture import CircularAperture,CircularAnnulus
from photutils.aperture import aperture_photometry
from photutils.aperture import ApertureStats


class TOBI_img_redux:

    def __init__(self,rawpath,output='./output/',mbiasname='masterbias.fit',mflatname='masterflat.fit'):
        
        #some default imaging numbers
        self.gain=0.6
        self.nfilters=7
        
        #path and alike
        self.rawpath=rawpath
        self.output=output
        
        #self
        self.mbiasname=mbiasname
        self.mflatname=mflatname
        
        if not os.path.exists(self.output):
            os.makedirs(self.output)

            
            
    
    def makebias(self):

        print("*********************")
        print("Processing bias now")
        print("*********************")

        savename=self.mbiasname
        
        biasindex=np.where(self.datatype == "B")[0]

        allbias=[]
        for ii,kk in enumerate(self.alldata[biasindex]):
            bias=fits.open(kk)
            #check if true bias
            if(bias[0].header['EXPTIME']<1):
                allbias.append(bias[0].data)
            

        #make masterbias
        allbias=np.array(allbias)
        self.masterbias=np.median(allbias,0)
        
        #write output
        hdu = fits.PrimaryHDU(self.masterbias)
        hdul = fits.HDUList([hdu])
        hdul.writeto(self.output+'/'+savename)

        print("Masterbias written to ",self.output+'/'+savename)
                
        return
        

    def makeflats(self):

        print("*********************")
        print("Processing flats now")
        print("*********************")

        savename=self.mflatname

        #loop over 7 filter positions
        flatdict=['Ha','Hb','OIII','SII','g','r','i']

        flatindex=np.where(self.datatype == "F")[0]
        self.masterflat={}
        hdulist=[]
        hdulist.append(fits.PrimaryHDU([]))
        
        for findex in range(self.nfilters):
            print("Processing flat position ",findex+1, ": ",flatdict[findex])

            allflat=[]
            
            for ii,kk in enumerate(self.alldata[flatindex]):
                flat=fits.open(kk)

                #check if good flat
                medval=np.median(flat[0].data)
                position=flat[0].header['FILTER']
                if((medval > 10000) & (medval < 40000) & ('{}'.format(findex+1) in position)):
                    allflat.append((flat[0].data-self.masterbias)/flat[0].header['EXPTIME'])
                    
            allflat=np.array(allflat)
            if((allflat.shape)[0] <3):
                print("Skip this position, not enough flats")
                self.masterflat['{}'.format(findex+1)]=None
                thishdu=fits.ImageHDU([])
            else:
                thisflat=np.median(allflat,0)
                thisflat=thisflat/np.median(thisflat)
                self.masterflat['{}'.format(findex+1)]=thisflat
                thishdu=fits.ImageHDU(thisflat)
            hdulist.append(thishdu)
            

        hdul = fits.HDUList(hdulist)
        hdul.writeto(self.output+'/'+savename)

        print("Masterflat written to ",self.output+'/'+savename)
                
        return


    def makescience(self):


        print("*********************")
        print("Processing science now")
        print("*********************")
        
        sciindex=np.where(self.datatype == "S")[0]
    
        for ii,kk in enumerate(self.alldata[sciindex]):

            savename=(kk.split("/")[-1]).split('.fit')[0]+"_redux.fit"

            if not (os.path.isfile(self.output+"/"+savename)):
                        
                
                #open
                science=fits.open(kk)
                findex=((science[0].header['FILTER']).split(' ')[-1])
                texp=science[0].header['EXPTIME']
                
                img=((science[0].data-self.masterbias)/self.masterflat[findex]*self.gain)/texp
                
                
                #write output
                hdu = fits.PrimaryHDU(img)
                hdul = fits.HDUList([hdu])
                hdul.writeto(self.output+'/'+savename)

                print('Done with ', savename)
                
        return
    
    def redux(self,verbose=True):

        """

        Engine that drives the reduction. A tag "flat" and "dark" and "bias" is expected in the file 
        name for this to work. Anything else is assumed to be science 

        """

        #grab data
        self.alldata=np.array(glob.glob(self.rawpath+"/*fit"))

        #tag all as science
        self.datatype=np.full(len(self.alldata), "S")

        #now identify the flats, dark and bias
        for ii,myfile in enumerate(self.alldata):
            if("flat" in myfile):
                self.datatype[ii]="F"
            elif("bias" in myfile):
                self.datatype[ii]="B"
            elif("dark" in myfile):
                self.datatype[ii]="D"

        if(verbose is True):
            print("******************************")
            print("HERE IS THE CURRENT DATA REDUCTION PLAN")
            print("******************************")
            for a, b in zip(self.alldata,self.datatype):
                print(a, b)

        #read the dimension of the data
        mheader=(fits.open(self.alldata[0]))[0].header
        self.nx=int(mheader['NAXIS2'])
        self.ny=int(mheader['NAXIS1'])

        #go to bias processing
        if (os.path.isfile(self.output+"/"+self.mbiasname)):
            print("Masterbias exists, do not overwrite!")
            self.masterbias=fits.open(self.output+"/"+self.mbiasname)[0].data
        else:
            self.makebias()
            
        #now process the dark
        print("TDB... dark will be added in the future")

        #now handle the flats
        if (os.path.isfile(self.output+"/"+self.mflatname)):
            print("Masterflat exists, do not overwrite!")
            restoreflat=fits.open(self.output+"/"+self.mflatname)
            self.masterflat={}
            for ii in range(self.nfilters):
                self.masterflat['{}'.format(ii+1)]=restoreflat[ii+1].data
        else:
            self.makeflats()

        #now handle the science
        self.makescience()
