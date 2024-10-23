#This version includes the constraints on pointing direction but calculates the orientation of the image afterwards instead of working out a tilt angle and applying it in retrospect. 

import numpy as np 
import matplotlib.pyplot as plt
from time import process_time
import os
from matplotlib.patches import Wedge, Polygon, Circle, Arc
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class smile_limb():
    '''This object will use the spacecraft position and limb angle to work out 
    the pointing and target directions, along with everything else.''' 
    
    def __init__(self, theta_fov=27, phi_fov=16, n_pixels=4, m_pixels=2, smile_loc=(0,-10,10), p_spacing=0.5, p_max=80, limb=20.3):
        '''This takes in all the initial parameters 
        
        Parameters
        ----------
        theta_fov - FOV angle (deg) in the theta direction (camera coords)
        phi_fox - FOV angle (deg) in the phi direction (camera coords)
        n_pixels - Number of pixels in the theta direction (camera coords)
        m_pixels - Number of pixels in the phi direction (camera coords)
        smile_loc - vector for the position of smile in magnetospheric xyz coords. 
        p_spacing - space in RE along LOS between points at which to calculate. 
        p_max - maximum distance in RE from spacecraft it will integrate to. 
            If it is set to None, it will be made equal to the length of vector L, 
            the look vector to the aim point.  
        
        '''
        
        self.plot_path = os.environ.get('PLOT_PATH')
        
        self.theta_fov = np.deg2rad(theta_fov)
        self.phi_fov = np.deg2rad(phi_fov)
        self.n_pixels = n_pixels
        self.m_pixels = m_pixels
        self.smile_loc = np.array(smile_loc)
        self.limb = np.deg2rad(limb) 
        self.p_spacing = p_spacing
        self.p_max = p_max 
        
        #Get magnitude of spacecraft vector - its radial position. 
        self.smag = np.sqrt(self.smile_loc[0]**2 + self.smile_loc[1]**2 + self.smile_loc[2]**2)
        
        #Get alpha angle. 
        self.get_alpha_angle()
        
        #Get earth angle. 
        self.r_angle = np.arcsin(1/self.smag) 
        
        #Get combined limb and r angle. 
        #Added correction - 20.3 deg to centre of FOV, not edge. 
        self.limb_c = self.limb + self.r_angle #+ self.phi_fov/2.
        
        #Get look direction vector. 
        self.get_look_direction()
        
        #Get target vector (aim). 
        self.target_loc = self.smile_loc + self.L 
        
        #Get b vector. 
        self.get_b_vector()
        
        #Get the colatitude of the look direction. 
        cos_colat = self.L[2]/self.Lmag 
        self.sxi_theta = np.arccos(cos_colat) 
        
        #Get the longitude of the look direction. 
        self.sxi_phi = np.arctan2(self.L[1], self.L[0])
        
        #Set LoS calculations to only got to the target point if p_max is None. 
        if self.p_max is None:
            self.p_max = self.Lmag
        
        #Add constraint variables. 
        if self.smag > 7.84: 
            self.radial_constraint = True 
        else:
            self.radial_constraint = False
            
        if self.limb_c-self.alpha_angle < np.deg2rad(90-35.83):
            self.solar_constraint = True
        else:
            self.solar_constraint = False 
        
        
            
        #Get unit vectors for the image. 
        self.xi_unit = self.b_unit

        self.yi = np.cross(self.b, self.L) 
        self.yi_unit = self.yi/(self.yi[0]**2 + self.yi[1]**2 + self.yi[2]**2)**0.5 
        
        #Get the arrays of angles for the pixels. 
        self.get_theta_and_phi_all_pixels() 
        
        #Get the unit vectors for each pixel. 
        self.get_pixel_unit_vectors() 
        
        #Get the colatitude and longitude of each pixel. 
        self.get_pixel_colat_longitude()
        
        #Get LOS coordinates. 
        self.get_LOS_coords()
        
     #FUNCTIONS TO GET SOME ANGLES/VECTORS FOR LOOK DIRECTION. 
                
    def get_alpha_angle(self):
        '''This will calculate alpha, the angle between the spacecraft vector and the perpendicular to the x axis. '''
        self.alpha_angle = np.arctan2(self.smile_loc[0],np.sqrt(self.smile_loc[1]**2 + self.smile_loc[2]**2))
        
         
    def get_look_direction(self):
        '''This will get the vector from the spacecraft to the target point on the x axis. '''
        
        lx = np.sqrt(self.smile_loc[1]**2 + self.smile_loc[2]**2)*np.tan(self.limb_c-self.alpha_angle) 
        
        ly = -self.smile_loc[1]
        lz = -self.smile_loc[2]
        
        self.L = np.array([lx, ly, lz]) 
        self.Lmag = np.sqrt(lx**2 + ly**2 + lz**2) 
        self.L_unit = self.L/self.Lmag
    
    def get_b_vector(self):
        '''This is a vector perpendicular to the look direction that 
        points towards the Earth.'''
        
        self.b = self.smile_loc + (self.smag*np.cos(self.limb_c)*self.L)/self.Lmag
        self.bmag = np.sqrt(self.b[0]**2 + self.b[1]**2 + self.b[2]**2) 
        self.b_unit = self.b/self.bmag 

    #CODE FOR GETTING LOOK DIRECTIONS OF EACH PIXEL. 
    
    def get_theta_and_phi_all_pixels(self):
        '''This will calculate theta and phi for all pixels. 
        It uses the method in Jorgensen et al. 
        These angles are in the directions of xi_unit for phi and yi_unit for theta.'''
 
        # Create 2D arrays for i and j. 
        self.J, self.I = np.meshgrid(np.arange(self.m_pixels), np.arange(self.n_pixels))

        # Calculate theta and phi for each pixel. 
        #THETA has changed as it's now relative to the look direction, not a colatitude. 
        self.theta_pixels = - (self.theta_fov/2.) + (self.theta_fov/self.n_pixels)*(self.I+0.5)
        self.phi_pixels = -(self.phi_fov/2.) + (self.phi_fov/self.m_pixels)*(self.J+0.5)
    
    def get_pixel_unit_vectors(self):
        '''This will work out the pixel unit vectors relative to the look direction.
        NEW METHOD.'''
        
        #Get vectors to flat plane of image. 
        self.pixels_x = self.L_unit[0] + np.tan(self.phi_pixels)*self.xi_unit[0] + np.tan(self.theta_pixels)*self.yi_unit[0]
        self.pixels_y = self.L_unit[1] + np.tan(self.phi_pixels)*self.xi_unit[1] + np.tan(self.theta_pixels)*self.yi_unit[1]
        self.pixels_z = self.L_unit[2] + np.tan(self.phi_pixels)*self.xi_unit[2] + np.tan(self.theta_pixels)*self.yi_unit[2]
        
        #Get magnitude. 
        self.pixels_mag = np.sqrt(self.pixels_x**2 + self.pixels_y**2 + self.pixels_z**2) 
        
        #Get unit vectors to unit sphere. 
        self.pixels_x_final = self.pixels_x/self.pixels_mag 
        self.pixels_y_final = self.pixels_y/self.pixels_mag 
        self.pixels_z_final = self.pixels_z/self.pixels_mag 
    
    def get_pixel_colat_longitude(self):
        '''Get the colatitude and longitude of each pixel. May be useful in future.'''

        #Colatitude. 
        self.pixels_colatitude = np.arccos(self.pixels_z_final/self.pixels_mag)  
        
        #Longitude. 
        self.pixels_longitude = np.arctan2(self.pixels_y_final,self.pixels_x_final)  
    
    def get_LOS_coords(self):
        '''This will calculate the coordinates along the LOS for a given 
        coordinate spacing. Same as before. '''

        # Array of points along any given LOS. 
        p = np.arange(0,self.p_max, step=self.p_spacing)+self.p_spacing
       
        self.xpos = np.zeros((self.n_pixels, self.m_pixels, p.size))
        self.ypos = np.zeros((self.n_pixels, self.m_pixels, p.size))
        self.zpos = np.zeros((self.n_pixels, self.m_pixels, p.size))

        # For each pixel: 
        for i in range(self.n_pixels):
            for j in range(self.m_pixels):

                #Positions along LOS. 
                self.xpos[i][j] = self.smile_loc[0] + p*self.pixels_x_final[i][j]
                self.ypos[i][j] = self.smile_loc[1] + p*self.pixels_y_final[i][j]
                self.zpos[i][j] = self.smile_loc[2] + p*self.pixels_z_final[i][j]
                
                
                
    #PLOTTING FUNCTIONS. 
    
        
    def plot_vectors(self, elev=45, azim=45):
        '''This will plot all the vectors to make sense of them.'''
        #plt.close("all")
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d') 
        
        #Add SMILE vector. 
        ax.plot([0, self.smile_loc[0]], [0, self.smile_loc[1]], [0, self.smile_loc[2]], 'b-', label='SMILE') 
        
        #Add Look vector. 
        ax.plot([self.smile_loc[0], self.smile_loc[0]+self.L[0]], [self.smile_loc[1], self.smile_loc[1]+self.L[1]], [self.smile_loc[2], self.smile_loc[2]+self.L[2]], 'r-', label='Look')
        
        #Add Target vector. 
        ax.plot([0, self.target_loc[0]], [0, self.target_loc[1]], [0, self.target_loc[2]], 'g-', label='Aim') 
        
        #Add b vector. 
        ax.plot([0, self.b[0]], [0, self.b[1]], [0, self.b[2]], 'c-', label='b') 
        
        #Add xi unit vector. 
        ax.plot([self.b[0], self.b[0]+self.xi_unit[0]], [self.b[1], self.b[1]+self.xi_unit[1]], [self.b[2], self.b[2]+self.xi_unit[2]], 'k-', label='xi')
        
        #Add yi unit vector. 
        ax.plot([self.b[0], self.b[0]+self.yi_unit[0]], [self.b[1], self.b[1]+self.yi_unit[1]], [self.b[2], self.b[2]+self.yi_unit[2]], c='orange', label='yi')
        
        self.add_earth(ax)
        
        #Sort legend and labels. 
        ax.legend(loc='best')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z') 
        ax.set_aspect('equal')
        ax.view_init(elev=elev, azim=azim)
        
        #Save figure. 
        fig.savefig(self.plot_path+'limb_example.png')

    def plot_LOS_vectors(self, elev=45, azim=45):
        '''This will plot the LOS vectors from the spacecraft position.'''
        plt.close("all")
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        #Add axes. 
        ax.plot([-2,12],[0,0],[0,0],'y')
        ax.plot([0,0],[-10,10],[0,0],'y')
        ax.plot([0,0],[0,0],[-10,10],'y')

        # For each pixel: 
        for i in range(self.n_pixels):
            for j in range(self.m_pixels):
                ax.plot(self.xpos[i][j], self.ypos[i][j], self.zpos[i][j], 'k', lw=0.2)
        
        
        #Add SMILE vector. 
        ax.plot([0, self.smile_loc[0]], [0, self.smile_loc[1]], [0, self.smile_loc[2]], 'b-', label='SMILE') 
        
        #Add Look vector. 
        ax.plot([self.smile_loc[0], self.smile_loc[0]+self.L[0]], [self.smile_loc[1], self.smile_loc[1]+self.L[1]], [self.smile_loc[2], self.smile_loc[2]+self.L[2]], 'r-', label='Look')
        
        #Add Target vector (aim). 
        ax.plot([0, self.target_loc[0]], [0, self.target_loc[1]], [0, self.target_loc[2]], 'g-', label='Target') 
        
        #Add b vector. 
        ax.plot([0, self.b[0]], [0, self.b[1]], [0, self.b[2]], 'c-', label='b') 
        
        #Add the FOV boundaries. 
        self.add_fov_boundaries(ax, lw=2) 
        
        self.add_fov_rectangle(ax, color='gray')
        
        #Add title to show smile location. 
        ax.set_title('SMILE: ({:.2f},{:.2f},{:.2f})\nAim: ({:.2f},{:.2f},{:.2f}) '.format(self.smile_loc[0], self.smile_loc[1], self.smile_loc[2], self.target_loc[0], self.target_loc[1], self.target_loc[2]))
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        self.add_earth(ax)
        ax.set_aspect('equal') 
        
        ax.view_init(elev,azim) 
    
    
    
        
    def add_fov_boundaries(self, ax2, color='k', lw=2):
        '''This will add the FOV boundaries in black/white. '''
        
        #For corner pixels only. 
        ax2.plot(self.xpos[0][0], self.ypos[0][0], self.zpos[0][0], color, lw=lw)
        ax2.plot(self.xpos[0][-1], self.ypos[0][-1], self.zpos[0][-1], color, lw=lw)
        ax2.plot(self.xpos[-1][0], self.ypos[-1][0], self.zpos[-1][0], color, lw=lw)
        ax2.plot(self.xpos[-1][-1], self.ypos[-1][-1], self.zpos[-1][-1], color, lw=lw)
        
        #Join corners together. 
        ax2.plot([self.xpos[0][0][-1],self.xpos[0][-1][-1]], [self.ypos[0][0][-1],self.ypos[0][-1][-1]], [self.zpos[0][0][-1],self.zpos[0][-1][-1]], color, lw=lw)
        ax2.plot([self.xpos[0][-1][-1],self.xpos[-1][-1][-1]], [self.ypos[0][-1][-1],self.ypos[-1][-1][-1]], [self.zpos[0][-1][-1],self.zpos[-1][-1][-1]], color, lw=lw)
        ax2.plot([self.xpos[-1][-1][-1],self.xpos[-1][0][-1]], [self.ypos[-1][-1][-1],self.ypos[-1][0][-1]], [self.zpos[-1][-1][-1],self.zpos[-1][0][-1]], color, lw=lw)
        ax2.plot([self.xpos[-1][0][-1],self.xpos[0][0][-1]], [self.ypos[-1][0][-1],self.ypos[0][0][-1]], [self.zpos[-1][0][-1],self.zpos[0][0][-1]], color, lw=lw)
    
    def add_fov_rectangle(self, ax, color='gray'):
        '''This will hopefully add a rectangle to the end of the FOV to make its shape clearer.'''
        
        v1 = [self.xpos[0][0][-1], self.ypos[0][0][-1], self.zpos[0][0][-1]] 
        v2 = [self.xpos[0][-1][-1], self.ypos[0][-1][-1], self.zpos[0][-1][-1]] 
        v3 = [self.xpos[-1][-1][-1], self.ypos[-1][-1][-1], self.zpos[-1][-1][-1]] 
        v4 = [self.xpos[-1][0][-1], self.ypos[-1][0][-1], self.zpos[-1][0][-1]] 
        
        rects = [[v1, v2, v3, v4, v1]]
        ax.add_collection3d(Poly3DCollection(rects, color=color, alpha=0.5, edgecolor=None))
        
        return 
        
           
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
