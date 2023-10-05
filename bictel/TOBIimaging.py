import glob
import astropy
from astropy.io import fits,ascii
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal as signal 
import os

from mypython.filters import filter as flt

from photutils.centroids import centroid_com
from photutils.aperture import CircularAperture,CircularAnnulus
from photutils.aperture import aperture_photometry
from photutils.aperture import ApertureStats

from astropy.stats import sigma_clipped_stats

class TOBI_img_redux:

    def __init__(self,rawpath,output='./output/',mbiasname='masterbias.fit',mflatname='masterflat.fit'):
        
        #some default imaging numbers
        self.gain=0.6
        self.flatdict=['Ha','Hb','OIII','SII','g','r','i']
        self.mypython_filters=[162,168,164,163,165,166,167]
        self.nfilters=len(self.flatdict)
        self.pixsize=0.44
        
        #path and alike
        self.rawpath=rawpath
        self.output=output
        
        #self
        self.mbiasname=mbiasname
        self.mflatname=mflatname
        
        if not os.path.exists(self.output):
            os.makedirs(self.output)
                 
    
    def makebias(self):


        """

        Routine that generate master bias
        
        """

        
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

        """

        Routine that generates coadded flats in each filter
        
        """

        
        print("*********************")
        print("Processing flats now")
        print("*********************")

        savename=self.mflatname

        flatdict=self.flatdict
        
        #loop over 7 filter positions
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

        """
        
        Routine that turns the raw science frame in reduced frames (e/s)

        """
        

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
                hdu = fits.PrimaryHDU(img,header=science[0].header)
                hdul = fits.HDUList([hdu])
                hdul.writeto(self.output+'/'+savename)

                print('Done with ', savename)


        self.allredux=glob.glob(self.output+"/*redux.fit")
                
        print("All done with science")
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



    def zeropoint(self,stdname,listfile=None,guess_x=None,guess_y=None,check=False):


        """
        
        Procedure that computes the zeropoint


        """

        #if list not given, use all redux files


        #load standard
        #col1 wave in A
        #col2 f_lambda 1e-16 erg/s/cm2/A
        std=ascii.read(os.getenv("MYPYTHON")+"/bictel/data/"+stdname)


        allzp=[]
        allfilt=[]
        allsky=[]
        
        if(listfile is None):
            listfile=self.allredux


        if((guess_x is None) | (guess_y is None)):
            guess_x=int(self.nx/2)
            guess_y=int(self.ny/2)

        for ffile in listfile:
            
            print("Processing zeropoint from image ", ffile)

            #data 
            fdata=fits.open(ffile)
            skymean, skymedian, skystd = sigma_clipped_stats(fdata[0].data, sigma=3.0)
            
            #select star
            trimdata=fdata[0].data[guess_x-100:guess_x+100,guess_y-100:guess_y+100]-skymedian
            x1, y1 = centroid_com(trimdata)
            
            #recenter
            trimdata=fdata[0].data[int(y1+guess_x-100)-100:int(y1+guess_x-100)+100,int(x1+guess_y-100)-100:int(x1+guess_y-100)+100]-skymedian
            x1, y1 = centroid_com(trimdata)


            if(check):
                #display
                plt.imshow(trimdata)
                plt.scatter([x1],[y1])
                plt.show()
        

            #now get phot profile
            
            #find bkg in annulus
            annulus_aperture = CircularAnnulus((y1,x1), r_in=70, r_out=90)
            aperstats = ApertureStats(trimdata, annulus_aperture)
            
            radii=[2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,40,50,60]
            phot=[]
            for r in radii:
                aperture = CircularAperture((y1,x1), r=r)
                phot_table = aperture_photometry(trimdata, aperture)
                phot.append(phot_table['aperture_sum']-aperstats.mean*aperture.area)
                

            phot=np.array(phot)
            phot_std=phot[-2][0]


            if(check):
                plt.plot(radii,phot)
                plt.show()
                
            
            #grab filter
            findex=int((fdata[0].header['FILTER']).split(' ')[-1])
            fdata.close()

            #find standrad flux in filter
            filT=flt.Filter(self.mypython_filters[findex-1])
            filT.loadtrans()
            std_fl=filT.convolve(std['col2']*1e-16,std['col1'])
            
            zp=np.log10(phot_std/std_fl)


            print("STD flux [e/s]",phot_std)
            print("STD flux [erg/s/cm2/A]",std_fl)
            print("{} ZP log([e/s]/[erg/s/cm2/A])".format(self.flatdict[findex-1]),zp)


            #go to SB (erg/s/cm2/A/arcsec2)
            skySB=np.log10(skymedian)-zp-2*self.pixsize
            print("Sky SB log([erg/s/cm2/A/arcsec2])",skySB)
            

            allfilt.append(findex)
            allsky.append(skySB)
            allzp.append(zp)
            
        return np.array(allzp),np.array(allsky),np.array(allfilt)
