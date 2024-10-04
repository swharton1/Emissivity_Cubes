from astropy.io import fits as pyfits
import numpy as np 
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Polygon, Circle

from . import read_batsrus 
from . import get_meridians as gm 

class convert_to_fits():
    '''This class will read in an ASCII BATSRUS file and convert it to a fits file.'''
    
    def __init__(self, filename='px_uni_1000.dat'):
        '''This takes in the filename and reads in the data from the ASCII file.'''
        
        self.filename = filename 
        self.data_path = os.environ.get('BATSRUS_PATH')
        
        self.batsrus = read_batsrus.read_batsrus_cube(filename=self.filename) 
    
    def make_fits(self):
        '''This will make the fits file''' 
        
        #Create filename. 
        self.filename_fits = os.path.join(self.data_path,self.filename[0:-4]+'.fits')
        
        #Create a new HDU object. 
        self.hdu = pyfits.PrimaryHDU()
        
        #Add emissivity data. 
        self.hdu.data = self.batsrus.eta_3d
        self.hdu.header 
        
        #Add the solar wind information. 
        self.hdu.header['bx'] = self.batsrus.bx
        self.hdu.header['by'] = self.batsrus.by
        self.hdu.header['bz'] = self.batsrus.bz
        self.hdu.header['vx'] = self.batsrus.vx
        self.hdu.header['vy'] = self.batsrus.vy
        self.hdu.header['vz'] = self.batsrus.vz
        self.hdu.header['density'] = self.batsrus.density
        self.hdu.header['pdyn'] = self.batsrus.dyn_pressure 
        self.hdu.header['pmag'] = self.batsrus.mag_pressure
        self.hdu.header['temp'] = self.batsrus.temp 
        
        #Add HDUs for x, y and z arrays. 
        self.hdux = pyfits.ImageHDU(data=self.batsrus.x, name='x')
        self.hduy = pyfits.ImageHDU(data=self.batsrus.y, name='y')
        self.hduz = pyfits.ImageHDU(data=self.batsrus.z, name='z')
        
        #Create a HDU list. 
        self.hdul = pyfits.HDUList([self.hdu, self.hdux, self.hduy, self.hduz]) 
        
        #Write the fits file. 
        self.hdul.writeto(self.filename_fits, overwrite=True) 
        print ('Created: {}'.format(self.filename_fits)) 
        
        
        
        
        
        
