from astropy.io import fits as pyfits
import numpy as np 
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Polygon, Circle

from . import read_batsrus 
from SXI_Core import get_meridians as gm 

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
        self.hdu.header['bx'] = (self.batsrus.bx, 'IMF Bx [nT]')
        self.hdu.header['by'] = (self.batsrus.by, 'IMF By [nT]')
        self.hdu.header['bz'] = (self.batsrus.bz, 'IMF Bz [nT]')
        self.hdu.header['vx'] = (self.batsrus.vx, 'Solar Wind Vx [km/s]')
        self.hdu.header['vy'] = (self.batsrus.vy, 'Solar Wind Vy [km/s]')
        self.hdu.header['vz'] = (self.batsrus.vz, 'Solar Wind Vz [km/s]')
        self.hdu.header['density'] = (self.batsrus.density, 'Solar Wind Density [/cm^3]')
        self.hdu.header['pdyn'] = (self.batsrus.dyn_pressure, 'Solar Wind Dynamic Pressure [nPa]') 
        self.hdu.header['pmag'] = (self.batsrus.mag_pressure, 'Solar Wind Magnetic Pressure [nPa]')
        self.hdu.header['temp'] = (self.batsrus.temp, 'Solar Wind Temperature [K]')  
        
        #Add solar wind flux in cm2/s to header. 
        flux = self.batsrus.vx*1e5*self.batsrus.density
        self.hdu.header['flux'] = (flux, 'Solar Wind Proton Flux [cm^2/s]') 
        
        #Add history comment of where the original file came from. 
        self.hdu.header['original'] = self.filename 
        
        #Add HDUs for x, y and z arrays. 
        self.hdux = pyfits.ImageHDU(data=self.batsrus.x, name='x')
        self.hduy = pyfits.ImageHDU(data=self.batsrus.y, name='y')
        self.hduz = pyfits.ImageHDU(data=self.batsrus.z, name='z')
        
        #Create a HDU list. 
        self.hdul = pyfits.HDUList([self.hdu, self.hdux, self.hduy, self.hduz]) 
        
        #Write the fits file. 
        self.hdul.writeto(self.filename_fits, overwrite=True) 
        print ('Created: {}'.format(self.filename_fits)) 
        
        
        
        
        
        
