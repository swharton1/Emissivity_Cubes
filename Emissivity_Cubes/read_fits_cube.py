import numpy as np 
from matplotlib.patches import Wedge, Polygon, Circle
import os
import matplotlib.pyplot as plt 
from astropy.io import fits as pyfits 

from SXI_Core import get_meridians as gm 
from SXI_Core import get_earth 
from SXI_Core import calc_pressures 


#THIS IS DIFFERENT TO THE ONE IN SXI_CORE. 

class read_fits_cube():
    '''This class will read in either the PPMLR, OpenGGCM or BATSRUS fits file and add the data to a python object in a common format.''' 
    
    def __init__(self, filename='OG_flow_cusps2_n05_vx-400vy30vz00_bx00by00bz05.fits', filetype='OpenGGCM', xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None): 
        '''This takes in the fits filename and reads it.
        
        Parameters
        ----------
        filename - Name of the FITS Emissivity Cube File. 
        filetype - Type of MHD Simulation. 'PPMLR', 'OpenGGCM', 'BATSRUS'
        
        xmin - Minimum x value to include. 
        xmax - Maximum x value to include. 
        ymin - Minimum y value to include. 
        ymax - Maximum y value to include. 
        zmin - Minimum z value to include. 
        zmax - Maximum z value to include. 
        
        ''' 
        
        #Select path to correct type of emissivity cube. 
        self.filetype = filetype 
        if filetype.upper() == 'PPMLR':
            self.data_path = os.environ.get("PPMLR_PATH")
        elif filetype.upper() == 'OPENGGCM':
            self.data_path = os.environ.get("OPENGGCM_PATH")
        elif filetype.upper() == 'BATSRUS':
            self.data_path = os.environ.get("BATSRUS_PATH") 
        else: 
            raise ValueError("Invalid filetype selected. Choose 'PPMLR', 'OPENGGCM' or 'BATSRUS'. Other MHD filetypes may be added later.")  
        
        self.filename=filename 
        self.plot_path = os.environ.get("PLOT_PATH") 
        
        self.filename_fits = os.path.join(self.data_path, filename) 
        
        # Check file exists. If it does, open it. 
        try:
            # Open the FITS file with this command. 
            with pyfits.open(self.filename_fits) as hdul:
                self.hdul = hdul 
                
                # Add the headers and data to the object. 
                self.primary_header = self.hdul[0].header 
                self.data_header = self.hdul[0].header 

                #Extract header information and add to object. 
                self.bx = np.float32(self.primary_header['BX']) 
                self.by = np.float32(self.primary_header['BY']) 
                self.bz = np.float32(self.primary_header['BZ']) 
                self.vx = np.float32(self.primary_header['VX']) 
                self.vy = np.float32(self.primary_header['VY']) 
                self.vz = np.float32(self.primary_header['VZ']) 
                self.density = np.float32(self.primary_header['DENSITY']) 
                self.temp = np.float32(self.primary_header['TEMP']) 


                # Get dynamic pressure.
                self.dyn_pressure = calc_pressures.calc_dynamic_pressure(self.vx, self.vy, self.vz, self.density)

                # Get magnetic pressure. 
                self.mag_pressure = calc_pressures.calc_magnetic_pressure(self.bx, self.by, self.bz)
                
                #Get the number of indices in the x, y and z directions. 
                self.n = np.array([self.primary_header['NAXIS1'], self.primary_header['NAXIS2'], self.primary_header['NAXIS3']]).astype('int32')

                #Get emissivity data. 
                self.eta_3d = self.hdul[0].data.astype('float32')  
                
                #Get x, y and z data. 
                self.x = self.hdul['x'].data.astype('float32')
                self.y = self.hdul['y'].data.astype('float32')
                self.z = self.hdul['z'].data.astype('float32')
                
                # Use meshgrid to get corresponding 3D x, y and z arrays.
                # Meshgrid associates y-value with 0-axis and x-value with 1-axis.  
                self.y_3d, self.z_3d, self.x_3d = np.meshgrid(self.y, self.z, self.x)

                # Bound the data in x,y,z here. 
                if xmin is not None: 
                    i = np.where(self.x_3d[0][0] >= xmin) 
                    self.x_3d = self.x_3d[:,:,i[0]]
                    self.y_3d = self.y_3d[:,:,i[0]]
                    self.z_3d = self.z_3d[:,:,i[0]]
                    self.eta_3d = self.eta_3d[:,:,i[0]]
                    
                    #Bound 1D x array too. 
                    j = np.where(self.x >= xmin) 
                    self.x = self.x[j]

                if xmax is not None: 
                    i = np.where(self.x_3d[0][0] < xmax) 
                    self.x_3d = self.x_3d[:,:,i[0]]
                    self.y_3d = self.y_3d[:,:,i[0]]
                    self.z_3d = self.z_3d[:,:,i[0]]
                    self.eta_3d = self.eta_3d[:,:,i[0]]
                
                    #Bound 1D x array too. 
                    j = np.where(self.x < xmax) 
                    self.x = self.x[j]
                    
                if ymin is not None: 
                    j = np.where(self.y_3d[0,:,0] >= ymin)
                    self.x_3d = self.x_3d[:,j[0],:]
                    self.y_3d = self.y_3d[:,j[0],:]
                    self.z_3d = self.z_3d[:,j[0],:]
                    self.eta_3d = self.eta_3d[:,j[0],:]

                    #Bound 1D x array too. 
                    j = np.where(self.y >= ymin) 
                    self.y = self.y[j]
                    
                if ymax is not None: 
                    j = np.where(self.y_3d[0,:,0] < ymax)
                    self.x_3d = self.x_3d[:,j[0],:]
                    self.y_3d = self.y_3d[:,j[0],:]
                    self.z_3d = self.z_3d[:,j[0],:]
                    self.eta_3d = self.eta_3d[:,j[0],:]

                    #Bound 1D x array too. 
                    j = np.where(self.y < ymax) 
                    self.y = self.y[j]
                    
                if zmin is not None:
                    k = np.where(self.z_3d[:,0,0] >= zmin)
                    self.x_3d = self.x_3d[k[0],:,:]
                    self.y_3d = self.y_3d[k[0],:,:]
                    self.z_3d = self.z_3d[k[0],:,:]
                    self.eta_3d = self.eta_3d[k[0],:,:]

                    #Bound 1D x array too. 
                    j = np.where(self.z >= zmin) 
                    self.z = self.z[j]
                    
                if zmax is not None:
                    k = np.where(self.z_3d[:,0,0] < zmax)
                    self.x_3d = self.x_3d[k[0],:,:]
                    self.y_3d = self.y_3d[k[0],:,:]
                    self.z_3d = self.z_3d[k[0],:,:]
                    self.eta_3d = self.eta_3d[k[0],:,:]

                    #Bound 1D x array too. 
                    j = np.where(self.z < zmax) 
                    self.z = self.z[j]

        except (FileNotFoundError, IOError):
            print ("Filename not found: {}".format(self.filename))
            
        self.get_subsolar_magnetopauses()

    def __repr__(self):
        return ("Class to read in the PPMLR emissivity cube from the fits file: {}".format(self.filename))



    def get_subsolar_magnetopauses(self):
        '''This will get a few definitions of the subsolar magnetopause.'''
        
        # For the slice with constant y. 
        y_uniq = abs(self.y_3d[0,:,0])
        i_y = np.where(y_uniq == min(y_uniq))[0][0]

        # For the slice with constant z. 
        z_uniq = abs(self.z_3d[:,0,0])
        i_z = np.where(z_uniq == min(z_uniq))[0][0]

        # Get data along sun-earth line. 
        xp = self.x_3d[i_z,i_y]
        yp = self.y_3d[i_z,i_y]
        zp = self.z_3d[i_z,i_y]
        etad = self.eta_3d[i_z,i_y]
        
        # Get max Ix second. 
        ix_index = np.where(etad == etad.max())
        self.maxIx = xp[ix_index][0]
        self.maxIx_eta = etad.max()
        
        # Get max dIx third. 
        # Get difference between etad values. 
        if self.filetype == 'OPENGGCM':
            dIx = -np.array([etad[i+1] - etad[i] for i in range(len(etad) - 1)])
        else:
            dIx = np.array([etad[i+1] - etad[i] for i in range(len(etad) - 1)])
            
        # Get centre point positions for radial direction. 
        xp_cent = xp + (xp[1]-xp[0])/2
        xp_cent = xp_cent[0:-1]

        dix_index = np.where(dIx == dIx.max())
        self.maxdIx = xp[dix_index][0]

        # Get f=0.25 
        dr = xp[ix_index] - xp[dix_index]
        self.f = xp[dix_index] + 0.25*dr[0] 
        
       
    def plot_both_planes(self, cmap="hot", levels=100, vmin=-8, vmax=-4, save=False, savetag=""):
        '''This will plot in the X-Z and X-Y planes side by side. 
        
        Parameters
        ----------
        cmap - matplotlib colourmap. Def = 'hot' 
        levels - number of levels in contour map. Def = 100. 
        vmin - minimum logged value on colourmap. All values below this will be set to this value. Def = -8 
        vmax - maximum logged value on colourmap. 
        save - boolean to save the plot to the PLOT_PATH variable. 
        savetag - string to add additional information to the end of the default file name. 
        
        '''

        #Get PPMLR data. 
        xp_y, zp_y, letad_y, xp_z, yp_z, letad_z = gm.get_log_emissivity(self.x_3d, self.y_3d, self.z_3d, self.eta_3d, vmin=vmin, vmax=vmax) 


        # Create a filename label so you know which file you plotted. 
        file_label = self.filename.split("/")[-1]

         # Get contour levels. 
        levels = np.linspace(vmin, vmax, levels+1)

        # Now you can make the contour plot. 
        fig = plt.figure()
        sw_units = ["cm"+r"$^{-3}$", "km s"+r"$^{-1}$", "nT"]
        label = file_label+"\nn = {:.2f} {}".format(self.density, sw_units[0])+", "+r"$v_x$ = {:.2f} {}".format(self.vx, sw_units[1])+r", $B_z$ = {:.2f} {}".format(self.bz, sw_units[2])
        fig.text(0.5,0.9, label, ha="center")
        fig.subplots_adjust(wspace=0.4, top=0.85)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)

        # Constant y. 
        cont1 = ax1.contourf(xp_y, zp_y, letad_y, cmap='hot', levels=levels, vmin=vmin, vmax=vmax)
        ax1.set_xlabel('X [RE]')
        ax1.set_ylabel('Z [RE]')
        ax1.set_title("XZ Plane")
        ax1.set_aspect("equal")

        # Add sketch of the Earth on top. 
        get_earth.make_earth(ax1, rotation=-90)

         # Constant z. 
        cont2 = ax2.contourf(xp_z, yp_z, letad_z, cmap='hot', levels=levels, vmin=vmin, vmax=vmax)
        ax2.set_xlabel('X [RE]')
        ax2.set_ylabel('Y [RE]')
        ax2.set_title("XY Plane")
        ax2.set_aspect("equal")

        # Add colourbar and earth. 
        cbar = plt.colorbar(cont2, ax=ax2, shrink=0.5)
        cbar.set_label(r"eV cm$^{-3}$ s$^{-1}$")
        level_min = int(np.ceil(cont2.levels.min()))
        level_max = int(np.floor(cont2.levels.max()))
        cticks = np.arange(level_min, level_max+1)
        cbar.set_ticks(cticks)
        cbar.set_ticklabels([r'$10^{'+str(i)+'}$' for i in cticks])

        # Add sketch of the Earth on top. 
        get_earth.make_earth(ax2, rotation=-90)

        # Option to save figure. 
        if save: 
            print ('Saved to: ', self.plot_path+"{}_both_planes{}.png".format(file_label, savetag))
            fig.savefig(self.plot_path+"{}_both_planes{}.png".format(file_label, savetag))




    def make_earth(self, ax, rotation=0):
        '''This will add a little plot of the Earth on top for reference. '''

        # Add white circle first. 
        r=1
        circle = Circle((0,0), r, facecolor='w', edgecolor='navy')
        ax.add_patch(circle)

        # Add nightside. 
        theta2 = np.arange(181)-180+rotation
        xval2 = np.append(r*np.cos(theta2*(np.pi/180)),0)
        yval2 = np.append(r*np.sin(theta2*(np.pi/180)),0)
        verts2 = [[xval2[i],yval2[i]] for i in range(len(xval2))]
        
        polygon2 = Polygon(verts2, closed=True, edgecolor='navy', facecolor='navy', alpha=1) 
        ax.add_patch(polygon2)    
        
        


