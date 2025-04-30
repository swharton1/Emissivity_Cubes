#This will convert a FITS cube, which is by default in GSM coordinates, to a cube in GSE coordinates. 

import numpy as np 
import matplotlib.pyplot as plt 
import datetime as dt 
from scipy.interpolate import interpn 
import os
from astropy.io import fits as pyfits
import glob
from SXI_Core import read_fits_cube
from . import transformations as trans 


def convert_all_openggcm_cubes(season='Summer'):
    '''This will list all the OpenGGCM cubes and loop through them to convert them all.'''
    
    #Get filenames. 
    path = '/data/smile/OpenGGCM/' 
    filenames = sorted(glob.glob1(path, f'{season}_*.00_fieldlines.fits'))
    
    
    #Convert each cube, and make the plot for sanity. 
    for f in filenames: 
        print (f)
        convert = convert_cube(filename=f, filetype='OpenGGCM') 
        convert.plot_yz_plane(save=True)#Optional. 
        convert.make_fits() 
        plt.close("all")



class convert_cube():
    '''This will read in a cube, calculate its a new grid of coordinates in GSE instead of GSM, then create a new file. Designed for OpenGGCM cubes, which have dipole information in them. '''
    
    def __init__(self, filename='Summer_00.00_fieldlines.fits', filetype='OpenGGCM'):
        '''This reads in the original ecube in GSM coordinates.''' 
        
        self.filename = filename 
        self.filetype = filetype 
        
        self.plot_path = os.environ.get('PLOT_PATH')
        
        #Read in the cube. Read in the whole cube except the tail. 
        print ('Read in original cube: ', self.filename) 
        self.ecube = read_fits_cube.read_fits_cube(filename=filename, filetype=filetype, xmin=-33) 
        
        #Extract dipole angles. These are not normally extracted. 
        self.dipole = self.ecube.primary_header['DIPOLE']
        self.dipole_y = self.ecube.primary_header['DIPOLE_Y'] 
        self.time = str(self.ecube.primary_header['TIME'])
        
        #Get datetime time. 
        self.dtime = dt.datetime.strptime(self.time, '%Y/%m/%d %H:%M')
        
        #Use same x grid as GSM grid. Use 3D grids. 
        x_1d = self.ecube.x
        y_1d = self.ecube.x
        z_1d = self.ecube.x
        
        #Now make 3D. 
        self.y_gse, self.z_gse, self.x_gse = np.meshgrid(y_1d, z_1d, x_1d) 
        
        shape = self.x_gse.shape 
        
        #Package up x, y and z then rotate around the z axis. 
        rot_angle = -np.deg2rad(self.dipole_y) 
        
        #Package x, y and z into a single vector.    
        positions = np.zeros((3,self.x_gse.size))
        positions[0] = self.x_gse.flatten()
        positions[1] = self.y_gse.flatten()
        positions[2] = self.z_gse.flatten()
 
        #Rotate coordinates to GSM. 
        print ('Rotate to GSM')
        pos_rot = trans.rotate_x(rot_angle, positions.T)        
        
        self.x_gsm = pos_rot[0].reshape(shape)
        self.y_gsm = pos_rot[1].reshape(shape)
        self.z_gsm = pos_rot[2].reshape(shape)
        

        #Now do the interpolation. 
        self.interpolate_to_new_grid() 
        

    def interpolate_to_new_grid(self):
        '''This will interpolate the old array to the new grid.'''
        
        #Set up empty array to fill with interpolated emissivity values. 
        self.peta = np.zeros((self.x_gsm.shape))
        
        #You will need this for interpn. 
        self.points_original = (self.ecube.z, self.ecube.y, self.ecube.x) 
              
        #You now need to setup the smile points in the correct format. 
        new_points = (self.z_gsm, self.y_gsm, self.x_gsm)
        
        #Calculate the interpolated emissivity values. 
        self.peta = interpn(self.points_original, self.ecube.eta_3d, new_points, method='linear') 
    

    def make_fits(self):
        '''This will make the fits file but with the new GSE interpolated cube instead of the original cube. ''' 
        
        #Create filename. 
        self.filename_fits = os.path.join(self.ecube.data_path,self.filename[0:-5]+'_GSE.fits')
        
        #Create a new HDU object. 
        self.hdu = pyfits.PrimaryHDU()
        
        #Add emissivity data. 
        self.hdu.data = self.peta
        self.hdu.header = self.ecube.primary_header

        #Add history comment of where the original file came from. 
        self.hdu.header['original'] = self.filename 
        self.hdu.header['history'] = 'Created using: Emissivity_Cubes/convert_cube_to_gse.py'
        
        #Add HDUs for x, y and z arrays. 
        self.hdux = pyfits.ImageHDU(data=self.ecube.x, name='x')
        self.hduy = pyfits.ImageHDU(data=self.ecube.x, name='y')
        self.hduz = pyfits.ImageHDU(data=self.ecube.x, name='z')
        
        #Create a HDU list. 
        self.hdul = pyfits.HDUList([self.hdu, self.hdux, self.hduy, self.hduz]) 
        
        #Write the fits file. 
        self.hdul.writeto(self.filename_fits, overwrite=True) 
        print ('Created: {}'.format(self.filename_fits)) 
        
              

    def plot_yz_plane(self, vmin=-8, vmax=-4, cmap='hot', levels=1000, save=False):
        '''This plotting is to check it has all been done correctly.''' 
        
        #Get the plane where x is near zero. 
        #Start with original. 
        xmin = abs(self.ecube.x).min() 
        self.i = np.where(self.ecube.x == xmin)[0][0] 
        
        self.xplane = self.ecube.x_3d[:,:,self.i]
        self.yplane = self.ecube.y_3d[:,:,self.i]
        self.zplane = self.ecube.z_3d[:,:,self.i]
        self.eta_plane = self.ecube.eta_3d[:,:,self.i]
        
        leta_plane = np.zeros(self.eta_plane.shape)+vmin
        i = np.where(self.eta_plane != 0)
        leta_plane[i] = np.log10(self.eta_plane[i])
        j = np.where(leta_plane < vmin)
        leta_plane[j] = vmin 
        
        #Now do the interpolated version. Cusps should align with the other side. 
        #Otherwise, you rotated it wrong! 
        yplane_int = self.y_gsm[:,:,self.i]
        zplane_int = self.z_gsm[:,:,self.i]
        eta_plane_int = self.peta[:,:,self.i]
        
        leta_int = np.zeros(eta_plane_int.shape)+vmin
        i = np.where(eta_plane_int != 0)
        leta_int[i] = np.log10(eta_plane_int[i])
        j = np.where(leta_int < vmin)
        leta_int[j] = vmin 
        
        
        # Get contour levels. 
        levels = np.linspace(vmin, vmax, levels+1)
        
        #Make basic figure. 
        fig = plt.figure(figsize=(8,8))
        fig.text(0.5, 0.95, f'{self.filename}', ha='center', fontsize=12)
        fig.subplots_adjust(hspace=0.4, wspace=0.4) 
        lims = (-50,50)
        
         
        ax1 = fig.add_subplot(221)
        ax1.contourf(self.yplane, self.zplane, leta_plane, cmap='hot', levels=levels, vmin=vmin, vmax=vmax)
        
        ax1.set(xlabel='Y GSM', ylabel='Z GSM') 
        ax1.set_aspect('equal')
        ax1.set_title(f'Original GSM')
        ax1.grid()
        ax1.set_xlim(lims)
        ax1.set_ylim(lims) 
        
        ax2 = fig.add_subplot(222)
        ax2.contourf(yplane_int, zplane_int, leta_int, cmap='hot', levels=levels, vmin=vmin, vmax=vmax)
        
        ax2.set(xlabel='Y GSM') 
        ax2.set_aspect('equal')
        ax2.set_xlim(lims)
        ax2.set_ylim(lims)
        ax2.set_title('Interpolated cube: GSM')
        ax2.grid()
        
        #Do the equivalent plots now in GSE coordinates. 
        ax3 = fig.add_subplot(223)
        
        #Package x, y and z into a single vector. 
        shape = self.xplane.shape   
        positions = np.zeros((3,self.xplane.size))
        positions[0] = self.xplane.flatten()
        positions[1] = self.yplane.flatten()
        positions[2] = self.zplane.flatten()
 
        #Rotate coordinates to GSM. 
        print ('Rotate to GSM')
        pos_rot = trans.rotate_x(np.deg2rad(self.dipole_y), positions.T)        
        
        self.xplane_gse = pos_rot[0].reshape(shape)
        self.yplane_gse = pos_rot[1].reshape(shape)
        self.zplane_gse = pos_rot[2].reshape(shape)
        
        ax3.contourf(self.yplane_gse, self.zplane_gse, leta_plane, cmap='hot', levels=levels, vmin=vmin, vmax=vmax)
        
        ax3.set(xlabel='Y GSE', ylabel='Z GSE') 
        ax3.set_aspect('equal')
        ax3.set_title(f'Original GSE')
        ax3.grid()
        ax3.set_xlim(lims)
        ax3.set_ylim(lims)
        
        
        ax4 = fig.add_subplot(224) 
        
        
        #Now do the interpolated version. Cusps should align with the other side. 
        #Otherwise, you rotated it wrong! 
        yplane_int_gse = self.y_gse[:,:,self.i]
        zplane_int_gse = self.z_gse[:,:,self.i]
        
        ax4.contourf(yplane_int_gse, zplane_int_gse, leta_int, cmap='hot', levels=levels, vmin=vmin, vmax=vmax)
        
        ax4.set(xlabel='Y GSE') 
        ax4.set_aspect('equal')
        ax4.set_xlim(ax1.get_xlim())
        ax4.set_ylim(ax1.get_ylim())
        ax4.set_title('Interpolated cube: GSE')
        ax4.grid()
        ax4.set_xlim(lims)
        ax4.set_ylim(lims)
        
        if save: 
            print ('Saving: ', self.plot_path+f'convert_to_gse_{self.filename}.png')
            fig.savefig(self.plot_path+f'convert_to_gse_{self.filename}.png')
