import numpy as np 
from matplotlib.patches import Wedge, Polygon, Circle
import os
import matplotlib.pyplot as plt 
from astropy.io import fits as pyfits 

from . import get_meridians as gm 

# This function will read the OpenGGCM ASCII files in exactly the same way as the PPMLR files. 

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
                self.dyn_pressure = self.calc_dynamic_pressure()

                # Get magnetic pressure. 
                self.mag_pressure = self.calc_magnetic_pressure()
                
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
        return ("Class to read in the PPMLR emissivity cube from the fits file: {}".format(self.filename))


    def calc_dynamic_pressure(self):
        '''Calculate this as it's a parameter in some models.'''

         # Assume average ion mass = mass of proton. 
        mp = 0.00000000000000000000000000167
        
        # Calculate v in m/s 
        v = (self.vx**2 + self.vy**2 + self.vz**2)**0.5
        v = v*1000 

        # Convert number of particles /cm^3 to /m^3. 
        n = self.density*1000000

        # Calculate dynamic pressure first in Pascals, then nPa. 
        dyn_pressure = 0.5*mp*n*(v**2)
        dyn_pressure = dyn_pressure*1000000000

        return (dyn_pressure)

    def calc_magnetic_pressure(self):
        '''Calculate the magnetic pressure'''

        # Calculate magnitude of B in T. 
        B = (self.bx**2 + self.by**2 + self.bz**2)**0.5
        B = B*0.000000001

        # mu0
        mu0 = 4*np.pi*0.0000001

        # Calculate magnetic pressure in Pa, then nPa. 
        mag_pressure = (B**2)/(2*mu0)
        mag_pressure = mag_pressure*1000000000

        return mag_pressure
        
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

        #Get meridian data for etad. 
        xp_y, yp_y, zp_y, etad_y, xp_z, yp_z, zp_z, etad_z = gm.calculate_meridian_planes(self.x_3d, self.y_3d, self.z_3d, self.eta_3d)

        # Calculate log10 eta values. If eta = 0, set log(eta) = vmin  
        letad_y = np.zeros(etad_y.shape)+vmin
        i = np.where(etad_y != 0)
        letad_y[i] = np.log10(etad_y[i])
        j = np.where(letad_y < vmin)
        letad_y[j] = vmin 


         # Calculate log10 eta values. If eta = 0, set log(eta) = vmin 
        letad_z = np.zeros(etad_z.shape)+vmin
        i = np.where(etad_z != 0)
        letad_z[i] = np.log10(etad_z[i])
        j = np.where(letad_z < vmin)
        letad_z[j] = vmin 


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
        self.make_earth(ax1, rotation=-90)

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
        self.make_earth(ax2, rotation=-90)

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
                self.dyn_pressure = self.calc_dynamic_pressure()

                # Get magnetic pressure. 
                self.mag_pressure = self.calc_magnetic_pressure()
                
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

    def calc_dynamic_pressure(self):
        '''Calculate this as it's a parameter in some models.'''

         # Assume average ion mass = mass of proton. 
        mp = 0.00000000000000000000000000167
        
        # Calculate v in m/s 
        v = (self.vx**2 + self.vy**2 + self.vz**2)**0.5
        v = v*1000 

        # Convert number of particles /cm^3 to /m^3. 
        n = self.density*1000000

        # Calculate dynamic pressure first in Pascals, then nPa. 
        dyn_pressure = 0.5*mp*n*(v**2)
        dyn_pressure = dyn_pressure*1000000000

        return (dyn_pressure)

    def calc_magnetic_pressure(self):
        '''Calculate the magnetic pressure'''

        # Calculate magnitude of B in T. 
        B = (self.bx**2 + self.by**2 + self.bz**2)**0.5
        B = B*0.000000001

        # mu0
        mu0 = 4*np.pi*0.0000001

        # Calculate magnetic pressure in Pa, then nPa. 
        mag_pressure = (B**2)/(2*mu0)
        mag_pressure = mag_pressure*1000000000

        return mag_pressure
    
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

        #Get meridian data for etad. 
        xp_y, yp_y, zp_y, etad_y, xp_z, yp_z, zp_z, etad_z = gm.calculate_meridian_planes(self.x_3d, self.y_3d, self.z_3d, self.eta_3d)

        # Calculate log10 eta values. If eta = 0 or eta < vmin, set log(eta) = vmin  
        letad_y = np.zeros(etad_y.shape)+vmin
        i = np.where(etad_y > 0)
        letad_y[i] = np.log10(etad_y[i])
        j = np.where(letad_y < vmin)
        letad_y[j] = vmin 

         # Calculate log10 eta values. If eta = 0 or eta < vmin, set log(eta) = vmin 
        letad_z = np.zeros(etad_z.shape)+vmin
        i = np.where(etad_z > 0)
        letad_z[i] = np.log10(etad_z[i])
        j = np.where(letad_z < vmin)
        letad_z[j] = vmin 

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
        self.make_earth(ax1, rotation=-90)

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
        self.make_earth(ax2, rotation=-90)

        # Option to save figure. 
        if save: 
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


#    def plot_scatter3D(self, s=2, marker=",", cmap="Greys", alpha=0.5, elev=45, azim=45):
#        '''This will just plot a basic scatter plot of the data. 
#        It could look messy...'''

#        import matplotlib.pyplot as plt 

#        fig = plt.figure()
#        ax = fig.add_subplot(111, projection='3d')
        # ax.contourf(self.x_3d[0], self.y_3d[0], self.eta_3d[0], cmap="Greys")
#        scat = ax.scatter3D(self.x_3d, self.y_3d, self.z_3d, c=np.log(self.eta_3d), s=s, marker=marker,cmap=cmap, alpha=alpha)
#        ax.set_xlabel("x")
#        ax.set_ylabel("y")
#        ax.set_zlabel("z")
        # Prev. y,z,x
#        ax.view_init(elev, azim)
#        cbar = plt.colorbar(scat, ax=ax)
#        cbar.set_label("Log Eta")
        
   
