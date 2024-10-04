#This version will base the target location and pointing directions on the correct spacecraft constraints around the limb angle. 

import numpy as np 
import matplotlib.pyplot as plt
from time import process_time
import os

class smile_limb():
    '''This object will use the spacecraft position and limb angle to work out the pointing and 
    target directions, along with everything else.''' 
    
    def __init__(self, theta_fov=27, phi_fov=16, n_pixels=4, m_pixels=2, smile_loc=(0,-10,10), p_spacing=0.5, p_max=80, limb=20.5):
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
        
        #Get magnitude of spacecraft vector. 
        self.smag = np.sqrt(self.smile_loc[0]**2 + self.smile_loc[1]**2 + self.smile_loc[2]**2)
        
        #Get alpha angle. 
        self.get_alpha_angle()
        
        #Get earth angle. 
        self.r_angle = np.arcsin(1/self.smag) 
        
        #Get combined limb, r angle and half FOV angle. 
        self.limb_c = self.limb + self.r_angle + self.phi_fov/2.
        
        #Get look direction vector. 
        self.get_look_direction()
        
        #Get target vector. 
        self.target_loc = self.smile_loc + self.L 
        
        #Get b vector. 
        self.get_b_vector()
        
        #Define GSE unit vectors. 
        self.x_unit = np.array([1,0,0])
        self.y_unit = np.array([0,1,0])
        self.z_unit = np.array([0,0,1])
        
        #Get nhat, the unit vector perpendicular to the vertical meridian containing the look vector. 
        self.n = np.cross(self.z_unit, self.L) 
        self.n_unit = self.n/np.sqrt(self.n[0]**2 + self.n[1]**2 + self.n[2]**2) 
        
        #Get the SXI tilt angle. 
        cos_tilt = np.dot(-self.b_unit, self.n_unit)
        self.sxi_tilt = -np.arccos(cos_tilt) 
        
        #Get the colatitude of the look direction. 
        cos_colat = self.L[2]/self.Lmag 
        self.sxi_theta = np.arccos(cos_colat) 
        
        #Get the longitude of the look direction. 
        self.sxi_phi = np.arctan2(self.L[1], self.L[0])
        
        print ('Tilt = ', np.rad2deg(self.sxi_tilt))
        print ('Colat = ', np.rad2deg(self.sxi_theta))
        print ('Long. = ', np.rad2deg(self.sxi_phi)) 
        
        #Set LoS calculations to only got to the target point. 
        self.p_max = self.Lmag
        
        ts = process_time()
        print ("Get theta and phi for each pixel:")
        self.get_theta_and_phi_all_pixels()
        te = process_time()
        print ("Time = {:.1f}s".format(te-ts))

        ts = process_time()
        print ("Get vector for each pixel:")
        self.get_vector_for_all_pixels()
        te = process_time()
        print ("Time = {:.1f}s".format(te-ts))

        ts = process_time()
        print ("Tilt camera: ")
        self.tilt_sxi_camera() 
        te = process_time()
        print ("Time = {:.1f}s".format(te-ts))

        ts = process_time()
        print ("Rotate camera: ")
        self.rotate_camera()
        te = process_time()
        print ("Time = {:.1f}s".format(te-ts))

        ts = process_time()
        print ("Get LOS coordinates: ")
        self.get_LOS_coords()
        te = process_time()
        print ("Time = {:.1f}s".format(te-ts))
        
    def get_alpha_angle(self):
        '''This will calculate alpha, the angle between the spacecraft vector and the perpendicular to the x axis. '''
        #cos_alpha = np.sqrt(self.smile_loc[1]**2 + self.smile_loc[2]**2)/self.smag
        #self.alpha_angle = np.arccos(cos_alpha) 
        #tan_alpha = self.smile_loc[0]/np.sqrt(self.smile_loc[1]**2 + self.smile_loc[2]**2)
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
        '''This is a vector perpendicular to the look direction that points towards the Earth.'''
        
        self.b = self.smile_loc + (self.smag*np.cos(self.limb_c)*self.L)/self.Lmag
        self.bmag = np.sqrt(self.b[0]**2 + self.b[1]**2 + self.b[2]**2) 
        self.b_unit = self.b/self.bmag 

    def plot_vectors(self, elev=45, azim=45):
        '''This will plot all the vectors to make sense of them.'''
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d') 
        
        #Add SMILE vector. 
        ax.plot([0, self.smile_loc[0]], [0, self.smile_loc[1]], [0, self.smile_loc[2]], 'k-', label='SMILE') 
        
        #Add Look vector. 
        ax.plot([self.smile_loc[0], self.smile_loc[0]+self.L[0]], [self.smile_loc[1], self.smile_loc[1]+self.L[1]], [self.smile_loc[2], self.smile_loc[2]+self.L[2]], 'r-', label='Look')
        
        #Add Target vector. 
        ax.plot([0, self.target_loc[0]], [0, self.target_loc[1]], [0, self.target_loc[2]], 'g-', label='Target') 
        
        #Add b vector. 
        ax.plot([0, self.b[0]], [0, self.b[1]], [0, self.b[2]], 'c-', label='b') 
        
        #Add vector n unit. 
        ax.plot([self.b[0], self.b[0]+self.n_unit[0]], [self.b[1], self.b[1]+self.n_unit[1]], [self.b[2], self.b[2]+self.n_unit[2]], 'm-', label='n_unit')
        
        #Sort legend and labels. 
        ax.legend(loc='best')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z') 
        ax.set_aspect('equal')
        ax.view_init(elev=elev, azim=azim)
        
        #Save figure. 
        fig.savefig(self.plot_path+'limb_example.png')

#THIS IS NOW THE SAME CODE AS BEFORE WORKING OUT THE UNIT VECTORS OF EACH PIXEL. 

    def get_theta_and_phi_all_pixels(self):
        '''This will calculate theta and phi for all pixels. 
        It uses the method in Jorgensen et al. '''
 
        # Create 2D arrays for i and j. 
        self.J, self.I = np.meshgrid(np.arange(self.m_pixels), np.arange(self.n_pixels))

        # Calculate theta and phi for each pixel. 
        self.theta_pixels = (np.pi/2.) - (self.theta_fov/2.) + (self.theta_fov/self.n_pixels)*(self.I+0.5)
        self.phi_pixels = -(self.phi_fov/2.) + (self.phi_fov/self.m_pixels)*(self.J+0.5)



    def get_vector_for_all_pixels(self):
        '''This will calculate a unit vector in xyz in camera coords for each pixel using its theta and phi. 
        '''

        # Create x, y and z arrays for the position vector of each pixel. 
        self.pixels_x = np.sin(self.theta_pixels)*np.cos(self.phi_pixels)
        self.pixels_y = np.sin(self.theta_pixels)*np.sin(self.phi_pixels)
        self.pixels_z = np.cos(self.theta_pixels)

    def tilt_sxi_camera(self):
        '''This will apply a camera tilt to the pixels that rotates them around the x-axis from the x-z plane.'''

        # This will rotate about the x axis. 
        self.pixels_x_tilted = self.pixels_x
        self.pixels_y_tilted = self.pixels_y*np.cos(self.sxi_tilt) - self.pixels_z*np.sin(self.sxi_tilt)
        self.pixels_z_tilted = self.pixels_y*np.sin(self.sxi_tilt) + self.pixels_z*np.cos(self.sxi_tilt)
                 

    def rotate_camera(self):
        '''This function will rotate the camera to the correct viewing direction, and rotate all the unit vectors for the pixels. '''   

        # Calculate the rotation angle a from theta. a is the increase in colatitude. 
        a = -(np.pi/2. - self.sxi_theta)
        
        # This will rotate about the y axis. 
        self.pixels_x_roty = self.pixels_x_tilted*np.cos(a) + self.pixels_z_tilted*np.sin(a)
        self.pixels_y_roty = self.pixels_y_tilted
        self.pixels_z_roty = -self.pixels_x_tilted*np.sin(a) + self.pixels_z_tilted*np.cos(a)
            
        # This will rotate about the z axis. 
        self.pixels_x_final = self.pixels_x_roty*np.cos(self.sxi_phi) - self.pixels_y_roty*np.sin(self.sxi_phi)
        self.pixels_y_final = self.pixels_x_roty*np.sin(self.sxi_phi) + self.pixels_y_roty*np.cos(self.sxi_phi)
        self.pixels_z_final = self.pixels_z_roty
              
    def get_LOS_coords(self):
        '''This will calculate the coordinates along the LOS for a given coordinate spacing.'''

        # Array of points along any given LOS. 
        p = np.arange(0,self.p_max, step=self.p_spacing)+self.p_spacing
        # pall = np.zeros((self.n_pixels, self.m_pixels, p.size))
        # pall[:,:] = p

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

    def plot_unit_vectors(self, elev=45, azim=0):
        '''This will plot the position vectors in 3D space.
        elev - sets the viewing elevation. 0 looks along the x-y plane. 
        azim - 0 looks along the x axis. 90 looks along the y axis. 
        
        '''

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        ax.scatter3D(self.pixels_x, self.pixels_y, self.pixels_z, c='k', marker='o')
        ax.scatter3D(self.pixels_x_tilted, self.pixels_y_tilted, self.pixels_z_tilted, c='r', marker='o')
        ax.scatter3D(self.pixels_x_roty, self.pixels_y_roty, self.pixels_z_roty, c='b', marker='o')
        ax.scatter3D(self.pixels_x_final, self.pixels_y_final, self.pixels_z_final, c='g', marker='o')
       
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_xlim(-1.2,1.2)
        ax.set_ylim(-1.2,1.2)
        ax.set_zlim(-1.2,1.2) 

        ax.view_init(elev,azim) 
    
        fig.savefig(self.plot_path+'limb_unit_vectors.png')
         
        
    def plot_LOS_vectors(self, elev=45, azim=45):
        '''This will plot the LOS vectors from the spacecraft position.'''

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
        
        #Add Target vector. 
        ax.plot([0, self.target_loc[0]], [0, self.target_loc[1]], [0, self.target_loc[2]], 'g-', label='Target') 
        
        #Add b vector. 
        ax.plot([0, self.b[0]], [0, self.b[1]], [0, self.b[2]], 'c-', label='b') 
        
        #Add the FOV boundaries. 
        self.add_fov_boundaries(ax, lw=2) 
        
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        self.add_earth(ax)
        ax.set_aspect('equal') 
        
        ax.view_init(elev,azim) 
        
        fig.savefig(self.plot_path+'limb_los_vectors_SMILE_({},{},{})_elev_{}_azim_{}.png'.format(self.smile_loc[0], self.smile_loc[1], self.smile_loc[2], elev, azim))
    
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
