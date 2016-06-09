import argparse
try:
    import pyregion as pyreg
except:
    print 'Cannot import pyregion'
import numpy as np

class PyMask(object):
    
    """Class to handle masks in fits files """

    def __init__(self, x, y, regfile, header=None):
        
	""" Parse the ds9 region file here. Also save the number 
	of regions and the x y dimensions of the data"""
	if(header):
            self.reg = pyreg.open(regfile).as_imagecoord(header=header)
        else:
            self.reg = pyreg.open(regfile)

	self.filter = self.reg.get_filter()
	self.x = x
	self.y = y
	self.nreg = len(self.reg)
    
    def fillmask(self, ind):
        
	""" Given a region object and an index fill the mask"""
	
	if (self.reg[ind]).coord_format == 'image':
	
	  self.mask = self.filter[ind].mask((self.y, self.x))
          self.maskname = (self.reg[ind]).name
	  self.maskcomm = (self.reg[ind]).comment
	
	else:
	 
	 raise TypeError('The region coordinates are not in image units')
	
    def write(self,outname):
        """Write the mask output"""

        from astropy.io import fits
        
        #PrimaryHDU object to encapsulate the data:
	#Mask is a boolean --> convert on the fly
        hdu = fits.PrimaryHDU(self.mask*1)
        hdu.header['REGTYPE'] = self.maskname
	hdu.header['REGNAME'] = self.maskcomm
	
	#write
        hdu.writeto(outname,clobber=True)


def main():
  
  #MFo If the code is used as a main code then perform the basic operations:
  #Parse ds9 file, extract regions, generate a fits mask for each region and save them
  
  p = argparse.ArgumentParser()
  p.add_argument("-x", action="store", default=500, type=int)
  p.add_argument("-y", action="store", default=500, type=int)
  p.add_argument("Ds9file")

  args = p.parse_args() 
  
  Mask = PyMask(args.x, args.y, args.Ds9file)
  for ii in range(Mask.nreg):
    Mask.fillmask(ii)
    Mask.write('Mask_{0}.fits'.format(ii))


if __name__ == "__main__":
  main()  

        
