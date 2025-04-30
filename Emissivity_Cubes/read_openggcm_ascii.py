import numpy as np 
from matplotlib.patches import Wedge, Polygon, Circle
import os
import matplotlib.pyplot as plt 
from astropy.io import fits as pyfits 

from SXI_Core import get_meridians as gm 
from SXI_Core import get_earth 
from SXI_Core import calc_pressures 



           
class read_openggcm_cube():
    '''This function will read in the OpenGGCM ASCII files to get the emissivity 
    data from model runs. 
    '''

    def __init__(self, filename="OpenGGCM_emissivity_simple_n05_vx-300vy30vz00_bx00by00bz05.txt", xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None):

        self.openggcm_path = os.environ.get("OPENGGCM_PATH")
        print (self.openggcm_path) 
        self.plot_path = os.environ.get("PLOT_PATH")
        
        self.filename=self.openggcm_path+filename

        try: 
            with open(self.filename, 'r') as f: 
                lines = f.readlines()

                header = lines[0:9]
                self.temp = np.float32(header[1][:-1])
                self.density = np.float32(header[2][:-1])
                self.vx = np.float32(header[3][:-1])
                self.vy = np.float32(header[4][:-1])
                self.vz = np.float32(header[5][:-1])
                self.bx = np.float32(header[6][:-1])
                self.by = np.float32(header[7][:-1])
                self.bz = np.float32(header[8][:-1])
                
                # Get dynamic pressure.
                self.dyn_pressure = calc_pressures.calc_dynamic_pressure(self.vx, self.vy, self.vz, self.density)

                # Get magnetic pressure. 
                self.mag_pressure = calc_pressures.calc_magnetic_pressure(self.bx, self.by, self.bz)
                
                
                # Get the number of x, y and z coords.
                self.n = np.array([np.int32(w) for w in lines[9][:-1].split(" ") if w.isdigit()])
                
                # Get the number of lines in the file with either x, y or z data in them. 
                xl, yl, zl = np.int32(np.ceil(self.n/6))
                self.lines = lines[10:]

                # Get the x coords. 
                self.x = []
                self.y = []
                self.z = []
                self.eta = []
                
                for i in range(len(self.lines)):
                    # Extract numbers from string. 
                    if i < xl:
                        [self.x.append(a) for a in self.lines[i][:-1].split(" ") if self.is_float(a)]
                    elif (i >= xl) & (i < xl+yl):
                        [self.y.append(a) for a in self.lines[i][:-1].split(" ") if self.is_float(a)]
                    elif (i >= xl+yl) & (i < xl+yl+zl):
                        [self.z.append(a) for a in self.lines[i][:-1].split(" ") if self.is_float(a)]
                    else:
                        [self.eta.append(a) for a in self.lines[i][:-1].split(" ") if self.is_float(a)]
                    i +=1 
                
                # Convert to numpy arrays. 
                self.x = np.array(self.x, dtype="float32")
                self.y = np.array(self.y, dtype="float32")
                self.z = np.array(self.z, dtype="float32")
                self.eta = np.array(self.eta, dtype="float64")
                
                # Reshape eta. The file has eta as eta[z][y][x]
                self.eta_3d = self.eta.reshape(self.n[::-1])

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

                if xmax is not None: 
                    i = np.where(self.x_3d[0][0] < xmax) 
                    self.x_3d = self.x_3d[:,:,i[0]]
                    self.y_3d = self.y_3d[:,:,i[0]]
                    self.z_3d = self.z_3d[:,:,i[0]]
                    self.eta_3d = self.eta_3d[:,:,i[0]]
                
                if ymin is not None: 
                    j = np.where(self.y_3d[0,:,0] >= ymin)
                    self.x_3d = self.x_3d[:,j[0],:]
                    self.y_3d = self.y_3d[:,j[0],:]
                    self.z_3d = self.z_3d[:,j[0],:]
                    self.eta_3d = self.eta_3d[:,j[0],:]

                if ymax is not None: 
                    j = np.where(self.y_3d[0,:,0] < ymax)
                    self.x_3d = self.x_3d[:,j[0],:]
                    self.y_3d = self.y_3d[:,j[0],:]
                    self.z_3d = self.z_3d[:,j[0],:]
                    self.eta_3d = self.eta_3d[:,j[0],:]

                if zmin is not None:
                    k = np.where(self.z_3d[:,0,0] >= zmin)
                    self.x_3d = self.x_3d[k[0],:,:]
                    self.y_3d = self.y_3d[k[0],:,:]
                    self.z_3d = self.z_3d[k[0],:,:]
                    self.eta_3d = self.eta_3d[k[0],:,:]

                if zmax is not None:
                    k = np.where(self.z_3d[:,0,0] < zmax)
                    self.x_3d = self.x_3d[k[0],:,:]
                    self.y_3d = self.y_3d[k[0],:,:]
                    self.z_3d = self.z_3d[k[0],:,:]
                    self.eta_3d = self.eta_3d[k[0],:,:]

        except (FileNotFoundError, IOError):
             print ("Filename not found: {}".format(self.filename))

    def __repr__(self):
        return ("Custom read_openggcm object for the file: {}".format(self.filename))
    
    def is_float(self, string):
        try:
            float(string)
            return True
        except ValueError:
            return False
  

    def eliminate_cusps(self, rmin=8):
        '''This will set any emissivity inside a certain radius to zero, as a way to eliminate
        some of the emission from the cusp regions.'''

        # Calculate r. 
        self.rmin = rmin
        r_3d = (self.x_3d**2 + self.y_3d**2 + self.z_3d**2)**0.5 
        i = np.where(r_3d < self.rmin)
        self.eta_3d[i] = 0 


    
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
        file_label = self.filename.split("/")[-1][:-4]
        
        figure_label = "OpenGGCM Emissivity Simulation: {} mask".format(file_label.split("_")[2]) 

         # Get contour levels. 
        levels = np.linspace(vmin, vmax, levels+1)

        # Now you can make the contour plot. 
        fig = plt.figure()
        sw_units = ["cm"+r"$^{-3}$", "km s"+r"$^{-1}$", "nT"]
        label = figure_label+"\nn = {:.2f} {}".format(self.density, sw_units[0])+", "+r"$v_x$ = {:.2f} {}".format(self.vx, sw_units[1])+r", $B_z$ = {:.2f} {}".format(self.bz, sw_units[2])
        fig.text(0.5,0.9, label, ha="center")
        fig.subplots_adjust(wspace=0.4, top=0.8)
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
            fig.savefig(self.plot_path+"{}_both_planes{}.png".format(file_label, savetag))


