#This will resample the PPMLR data into the SMILE FOV and work out an image. 

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interpn 

class ppmlr_image():
    '''This class takes in the ppmlr simulation object and the smile fov object and calculates an image through the simulation.'''
    
    def __init__(self, ppmlr, smile):
        '''This takes in the ppmlr and smile objects.'''
        
        self.ppmlr = ppmlr
        self.smile = smile 
        
        # Extract any useful solar wind parameters
        self.temp = ppmlr.temp
        self.density = ppmlr.density
        self.vx = ppmlr.vx
        self.vy = ppmlr.vy
        self.vz = ppmlr.vz
        self.bx = ppmlr.bx
        self.by = ppmlr.by
        self.bz = ppmlr.bz
        self.pdyn = self.calc_dynamic_pressure()
        self.pmag = self.calc_magnetic_pressure()
        
        #Calculate the emissivity along the LOS. 
        self.get_eta_in_fov()
        
        #Calculate the LOS intensity that forms the image. 
        self.calc_model_image()
        
        self.plot_path = os.environ.get("PLOT_PATH") 
    
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
        
    #THIS CODE WILL WORK OUT THE EMISSIVITY ALONG THE LOS FROM THE PPMLR MODEL.
    ###########################################################################
    def get_eta_in_fov(self):
        '''This will be a new method to do it using scipy's linear interpolation function. Should be cleaner and faster.''' 
        
        #Set up empty array to fill with interpolated emissivity values. 
        self.peta = np.zeros((self.smile.xpos.shape))
        
        #You will need this for interpn. 
        self.points_original = (self.ppmlr.z, self.ppmlr.y, self.ppmlr.x) 
        
        #Define edges of cube. 
        self.xmin = self.ppmlr.x.min()
        self.xmax = self.ppmlr.x.max()
        self.ymin = self.ppmlr.y.min()
        self.ymax = self.ppmlr.y.max()
        self.zmin = self.ppmlr.z.min()
        self.zmax = self.ppmlr.z.max()
        
        
        pxn = self.smile.xpos
        pyn = self.smile.ypos
        pzn = self.smile.zpos 
        
        #You should only do this on points inside the cube. 
        inside = np.where((pxn < self.xmax) & (pxn > self.xmin) & (pyn < self.ymax) & (pyn > self.ymin) & (pzn < self.zmax) & (pzn > self.zmin))
        
        #Work out which points are outside the cube. 
        outside = np.where((pxn >= self.xmax) | (pxn <= self.xmin) | (pyn >= self.ymax) | (pyn <= self.ymin) | (pzn >= self.zmax) | (pzn <= self.zmin)) 
            
        #You now need to setup the smile points in the correct format. 
        new_points = (pzn[inside], pyn[inside], pxn[inside]) 
        
        #Calculate the interpolated emissivity values. 
        self.peta[inside] = interpn(self.points_original, self.ppmlr.eta_3d, new_points, method='linear')
        
        #Set the values outside the cube to zero. 
        self.peta[outside] = 0 
    
 
        
    #THIS CODE WILL WORK OUT THE INTENSITY IMAGE
    ######################################################
    
    def trapezium_rule(self, p_spacing, eta_LOS):
        '''This will integrate a function using the trapezium rule. '''

        return (p_spacing/2)*(eta_LOS[0] + eta_LOS[-1] + 2*sum(eta_LOS[1:-1]))
    
    def calc_model_image(self):
        '''This is the function that will actually work out the LOS intensities for the given spacecraft viewing direction.''' 
        
        print ("Calculating LOS intensity...") 
        #Calculate the LOS intensity. 
        self.los_intensity = np.zeros((self.smile.n_pixels, self.smile.m_pixels))
        
        # For each pixel: 
        for i in range(self.smile.n_pixels):
            for j in range(self.smile.m_pixels):
            
                #Added unit conversion factor from ev.RE to kev.cm
                self.los_intensity[i][j] = ((1/(4*np.pi))*self.trapezium_rule(self.smile.p_spacing, self.peta[i][j]))*637100
    
    def calc_model_image_sections(self, r_bounds = [0,8,16,24,32,40]):
        '''This will create a set of images of the emission from different distances, like Andy's code.'''
        
        print ("Calculating LOS intensity for a range of distances...")
        
        self.r_bounds = r_bounds
        
        self.los_intensity_sections = np.zeros((len(self.r_bounds)-1, self.smile.n_pixels, self.smile.m_pixels)) 
        
        #Loop through each section. 
        for r in range(len(self.r_bounds)-1):
            
            # For each pixel: 
            for i in range(self.smile.n_pixels):
                for j in range(self.smile.m_pixels):
                    
                    #Filter to just use emissivities in this range. 
                    idx1 = int(self.r_bounds[r]/self.smile.p_spacing)
                    idx2 = int(self.r_bounds[r+1]/self.smile.p_spacing)
                    
                    #Added unit conversion factor from ev.RE to kev.cm
                    self.los_intensity_sections[r][i][j] = ((1/(4*np.pi))*self.trapezium_rule(self.smile.p_spacing, self.peta[i][j][idx1:idx2]))*637100
    
    
    
    #PLOTTING FUNCTIONS 
    ###################
    
    def plot_image(self, elev=45, azim=45, cmap='hot', vmin=-8, vmax=-4, levels=100, colour_cap=0, los_max=6, image_name=None, ellipse=None, save=False):
        '''This will plot the simulated image it has created. 
        
        Parameters
        ----------
        elev - viewing elevation in degrees for 3d viewing model
        azim - viewing azimuth in degrees for 3d viewing model 
        cmap - colourmap
        vmin - min emissivity on colour bar (logged)
        vmax - max emissivity on colour bar (logged)
        levels - number of levels on contour maps
        colour_cap - order of magnitude above vmin to start plotting emissivity 
        los_max - max los intensity on colourbar
        image_name - will default to standard name if not specified. must be full path
        ellipse - ellipse object if you wish to add orbit ellipse. 
        
        '''
        
        fig = plt.figure(figsize=(8,5))
        fig.subplots_adjust(left=0.05, wspace=0.2, bottom=0.20) 
        ax = fig.add_subplot(121) 
        
        
        # Make pixel arrays for plotting. 
        i_array = np.linspace(0,self.smile.n_pixels, self.smile.n_pixels+1)-0.5
        j_array = np.linspace(0,self.smile.m_pixels, self.smile.m_pixels+1)-0.5
        
        J, I = np.meshgrid(j_array, i_array)
        
        #Convert to degrees. 
        theta_pixels = - (self.smile.theta_fov/2.) + (self.smile.theta_fov/self.smile.n_pixels)*(I+0.5)
        phi_pixels = -(self.smile.phi_fov/2.) + (self.smile.phi_fov/self.smile.m_pixels)*(J+0.5)
        
        #Convert to degrees. Minus sign is so you look away from camera, not towards. 
        theta_pixels = -np.rad2deg(theta_pixels)
        phi_pixels = -np.rad2deg(phi_pixels)
        
        # Get contour levels. 
        levels = np.linspace(vmin, vmax, levels+1)
        
        mesh = ax.pcolormesh(phi_pixels, theta_pixels, self.los_intensity, cmap=cmap, vmin=0, vmax=los_max)
        ax.set_title("LOS integration through PPMLR\n simulation from SMILE")
        ax.set_xlabel('deg')
        ax.set_ylabel('deg')
        ax.set_aspect('equal')
        cbar = plt.colorbar(mesh, ax=ax, shrink=0.8)
        cbar.set_label('SWCX LOS Intensity (keV cm'+r'$^{-2}$ s'+r'$^{-1}$ sr'+r'$^{-1}$)') 
        
        ax2 = fig.add_subplot(122, projection='3d')

        #Add the emissivity along the LOS.     
        los_log = np.zeros(self.peta.shape)+vmin
        i = np.where(self.peta > 0)
        los_log[i] = np.log10(self.peta[i])
        
        #Only plot values above a certain emissivity. 
        bright = np.where(los_log > vmin+colour_cap) 
        
        emit = ax2.scatter(self.smile.xpos[bright], self.smile.ypos[bright], self.smile.zpos[bright], c=los_log[bright], cmap="hot", s=0.05, vmin=vmin, vmax=vmax)
        cbar2 = plt.colorbar(emit, ax=ax2, shrink=0.8)
        
        cticks = np.arange(vmin, vmax)
        cbar2.set_ticks(cticks)
        cbar2.set_ticklabels([r'$10^{'+str(i)+'}$' for i in cticks])
        cbar2.set_label('SWCX Emissivity (eV cm'+r'$^{-3}$ s'+r'$^{-1}$)') 
        #Add FOV boundaries. 
        self.add_fov_boundaries(ax2)
        
        #Add the Earth on. 
        self.add_earth(ax2) 
        
        #Add orbit ellipse if it is provided.
        if ellipse is not None: 
            ax2.plot(ellipse.x3/ellipse.RE, ellipse.y3/ellipse.RE, ellipse.z3/ellipse.RE, 'k')
        
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_zlabel('z')
        ax2.set_xlim(-10,30)
        ax2.set_ylim(-30,30)
        ax2.set_zlim(-30,30)
        ax2.set_title('n = {} cm'.format(self.density)+r'$^{-3}$'+'\nSMILE Coords: ({:.2f},{:.2f},{:.2f})\nAim Point: ({},{},{})'.format(self.smile.smile_loc[0], self.smile.smile_loc[1], self.smile.smile_loc[2], self.smile.target_loc[0], self.smile.target_loc[1], self.smile.target_loc[2]))
        ax2.set_aspect('equal')
        ax2.view_init(elev,azim) 

                    
        #Save the image to a standard name. 
        if save:
            if image_name is None: 
            
        
                fig.savefig(self.plot_path+"PPMLR_image_sim_n_{}_SMILE_{}_{}_{}_Target_{}_{}_{}_nxm_{}_{}.png".format(self.density, self.smile.smile_loc[0], self.smile.smile_loc[1], self.smile.smile_loc[2], self.smile.target_loc[0], self.smile.target_loc[1], self.smile.target_loc[2],self.smile.n_pixels, self.smile.m_pixels)) 
        
            else:
                fig.savefig(image_name) 
    
    def plot_image_sections(self, elev=45, azim=45, cmap='hot', vmin=-8, vmax=-4, levels=100, colour_cap=0, los_max=None, image_name=None, ellipse=None, save=False, savetag=''):
        '''This will plot the series of images created from emission at different distances, plus the total 
        at the bottom. Akin to Andy Read's code. 
        
        Parameters
        ----------
        elev - viewing elevation in degrees for 3d viewing model
        azim - viewing azimuth in degrees for 3d viewing model 
        cmap - colourmap
        vmin - min emissivity on colour bar (logged)
        vmax - max emissivity on colour bar (logged)
        levels - number of levels on contour maps
        colour_cap - order of magnitude above vmin to start plotting emissivity 
        los_max - max los intensity on colourbar
        image_name - will default to standard name if not specified. must be full path
        ellipse - ellipse object if you wish to add orbit ellipse. 
        
        '''
        
        plt.close("all") 
        fig = plt.figure(figsize=(8,8))
        fig.subplots_adjust(left=0.10, wspace=0.3, bottom=0.10) 
        
        n_axes = len(self.r_bounds) 
        rows = n_axes/2
        
        # Make pixel arrays for plotting. 
        i_array = np.linspace(0,self.smile.n_pixels, self.smile.n_pixels+1)-0.5
        j_array = np.linspace(0,self.smile.m_pixels, self.smile.m_pixels+1)-0.5
        
        J, I = np.meshgrid(j_array, i_array)
        
        #Convert to degrees. 
        theta_pixels = - (self.smile.theta_fov/2.) + (self.smile.theta_fov/self.smile.n_pixels)*(I+0.5)
        phi_pixels = -(self.smile.phi_fov/2.) + (self.smile.phi_fov/self.smile.m_pixels)*(J+0.5)
        
        #Convert to degrees. Minus sign is so you look away from camera, not towards. 
        theta_pixels = -np.rad2deg(theta_pixels)
        phi_pixels = -np.rad2deg(phi_pixels)
        
        # Get contour levels. 
        levels = np.linspace(vmin, vmax, levels+1)
    
        print (n_axes)
        for n in range(n_axes):
            #print (n, self.r_bounds[n], n%3, n//3)
            ax = fig.add_subplot(int(n_axes/rows), int(rows), n+1)
            if n == n_axes-1:
                #Do the total image. 
                mesh = ax.pcolormesh(phi_pixels, theta_pixels, self.los_intensity, cmap=cmap, vmin=0, vmax=los_max*3)
                ax.set_title("{} < r < {}".format(0, self.smile.p_max))
                
            else: 
                #Do the sections images. 
                mesh = ax.pcolormesh(phi_pixels, theta_pixels, self.los_intensity_sections[n], cmap=cmap, vmin=0, vmax=los_max)
                ax.set_title("{} < r < {}".format(self.r_bounds[n],self.r_bounds[n+1]))
            
            #Add customisation to the plot. 
            if n//3 == 1: ax.set_xlabel('deg')
            if n%3 == 0 : ax.set_ylabel('deg')
            ax.set_aspect('equal')
            cbar = plt.colorbar(mesh, ax=ax, shrink=0.8)
            if n%3 == int(n_axes/rows): cbar.set_label('SWCX LOS Intensity (keV cm'+r'$^{-2}$ s'+r'$^{-1}$ sr'+r'$^{-1}$)') 
    
        #Add title with information about FOV and simulation. 
        fig.text(0.5, 0.95, 'n = {} cm'.format(self.density)+r'$^{-3}$'+'\nSMILE Coords: ({:.2f},{:.2f},{:.2f}),   Aim Point: ({},{},{})'.format(self.smile.smile_loc[0], self.smile.smile_loc[1], self.smile.smile_loc[2], self.smile.target_loc[0], self.smile.target_loc[1], self.smile.target_loc[2]), ha='center', va='top', fontsize=12)
    
        if save:
            print ('Saved to: ', self.plot_path+'PPMLR_Image_Sections_n_{:.1f}_SMILE_({:.2f},{:.2f},{:.2f})_Target_({:.2f},{:.2f},{:.2f})_{}.png'.format(self.density, self.smile.smile_loc[0], self.smile.smile_loc[1], self.smile.smile_loc[2], self.smile.target_loc[0], self.smile.target_loc[1], self.smile.target_loc[2], savetag))  
            fig.savefig(self.plot_path+'PPMLR_Image_Sections_n_{:.1f}_SMILE_({:.2f},{:.2f},{:.2f})_Target_({:.2f},{:.2f},{:.2f})_{}.png'.format(self.density, self.smile.smile_loc[0], self.smile.smile_loc[1], self.smile.smile_loc[2], self.smile.target_loc[0], self.smile.target_loc[1], self.smile.target_loc[2], savetag))
    
       
    def add_fov_boundaries(self, ax2):
        '''This will add the FOV boundaries in black. '''
        
        #For corner pixels only. 
        ax2.plot(self.smile.xpos[0][0], self.smile.ypos[0][0], self.smile.zpos[0][0], 'k', lw=2)
        ax2.plot(self.smile.xpos[0][-1], self.smile.ypos[0][-1], self.smile.zpos[0][-1], 'k', lw=2)
        ax2.plot(self.smile.xpos[-1][0], self.smile.ypos[-1][0], self.smile.zpos[-1][0], 'k', lw=2)
        ax2.plot(self.smile.xpos[-1][-1], self.smile.ypos[-1][-1], self.smile.zpos[-1][-1], 'k', lw=2)
        
        #Join corners together. 
        ax2.plot([self.smile.xpos[0][0][-1],self.smile.xpos[0][-1][-1]], [self.smile.ypos[0][0][-1],self.smile.ypos[0][-1][-1]], [self.smile.zpos[0][0][-1],self.smile.zpos[0][-1][-1]], 'k')
        ax2.plot([self.smile.xpos[0][-1][-1],self.smile.xpos[-1][-1][-1]], [self.smile.ypos[0][-1][-1],self.smile.ypos[-1][-1][-1]], [self.smile.zpos[0][-1][-1],self.smile.zpos[-1][-1][-1]], 'k')
        ax2.plot([self.smile.xpos[-1][-1][-1],self.smile.xpos[-1][0][-1]], [self.smile.ypos[-1][-1][-1],self.smile.ypos[-1][0][-1]], [self.smile.zpos[-1][-1][-1],self.smile.zpos[-1][0][-1]], 'k')
        ax2.plot([self.smile.xpos[-1][0][-1],self.smile.xpos[0][0][-1]], [self.smile.ypos[-1][0][-1],self.smile.ypos[0][0][-1]], [self.smile.zpos[-1][0][-1],self.smile.zpos[0][0][-1]], 'k')
    

        
        
    def add_earth(self, ax):
        '''This will add a sphere for the Earth. '''
        
        #Create a spherical surface. 
        radius = 1
        u = np.linspace(0, 2*np.pi, 100) 
        v = np.linspace(0, np.pi, 100) 
        x = radius* np.outer(np.cos(u), np.sin(v))
        y = radius* np.outer(np.sin(u), np.sin(v))
        z = radius* np.outer(np.ones(np.size(u)), np.cos(v))

        ax.plot_surface(x, y, z, color='k', lw=0, alpha=1)
        
    def sig_figs(self, x: float, precision: int):
        """
        Rounds a number to number of significant figures
        Parameters:
        - x - the number to be rounded
        - precision (integer) - the number of significant figures
        Returns:
        - float
        """

        x = float(x)
        precision = int(precision)

        return np.round(x, -int(np.floor(np.log10(abs(x)))) + (precision - 1))
        
        
