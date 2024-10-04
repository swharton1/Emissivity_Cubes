import numpy as np 
import matplotlib.pyplot as plt
from time import process_time

# THIS IS AN UPDATED VERSION THAT SHOULD BE BETTER. 

class smile_fov():
    '''This object will contain all the informaton the orientation of SXI 
    and its pixels, and provide the vectors for the line of sight of each 
    pixel. '''
    def __init__(self, theta_fov=27, phi_fov=16, n_pixels=2, m_pixels=4, sxi_tilt=0, \
                 sxi_theta=60, sxi_phi=0, smile_loc=(4,-2,10), target_loc=None, \
                 p_spacing=0.5, p_max=80):
        '''This contains all the information needed to calculate everything.
        
        Parameters
        ----------
        theta_fov - FOV angle (deg) in the theta direction (camera coords)
        phi_fox - FOV angle (deg) in the phi direction (camera coords)
        n_pixels - Number of pixels in the theta direction (camera coords)
        m_pixels - Number of pixels in the phi direction (camera coords)
        sxi_tilt - Rotation of the camera (deg) in the anticlockwise direction 
                looking along the x axis (camera coords)
        sxi_theta - theta angle for direction of centre of SXI FOV (magnetospheric coords.). 
        sxi_phi - phi angle for direction of centre of SXI FOV (magnetospheric coords.).
        smile_loc - vector for the position of smile in magnetospheric xyz coords.
        target_loc - vector for the location SXI is looking at in magnetospheric xyz coords. Redundant of theta and phi are specified. 
        p_spacing - space in RE along LOS between points at which to calculate. 
        p_max - maximum distance in RE from spacecraft it will integrate to. 
        '''

        self.theta_fov = np.deg2rad(theta_fov)
        self.phi_fov = np.deg2rad(phi_fov)
        self.n_pixels = n_pixels
        self.m_pixels = m_pixels
        self.sxi_tilt = np.deg2rad(sxi_tilt)
        self.sxi_theta = np.deg2rad(sxi_theta)
        self.sxi_phi = np.deg2rad(sxi_phi)
        self.smile_loc = np.array(smile_loc)
        
        self.p_spacing = p_spacing
        self.p_max = p_max 
        
        # In the scenario where a target location has been given, 
        # calculate sxi_theta and sxi_phi. 
        if target_loc is not None: 
            
            p = target_loc - self.smile_loc
            self.sxi_theta = np.arccos(p[2]/((p[0]**2 + p[1]**2 + p[2]**2)**0.5))
            self.sxi_phi = np.arctan2(p[1], p[0])
        else:
            pass
            #target loc is just used to get sxi_theta and sxi_phi. If target_loc is not set, there is no 
            #guarantee it lies along the x axis. 
        self.target_loc = np.array(target_loc) 
        
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

        # For the last two functions, the program uses the 3D rotation matrices around each axis.
        # These can be found here: https://en.wikipedia.org/wiki/Rotation_matrix 


    def get_theta_and_phi_all_pixels(self):
        '''This will calculate theta and phi for all pixels. 
        It uses the method in Jorgensen et al. '''
 
        # Create 2D arrays for i and j. 
        self.J, self.I = np.meshgrid(np.arange(self.m_pixels), np.arange(self.n_pixels))

        # Calculate theta and phi for each pixel. 
        self.theta_pixels = (np.pi/2.) - (self.theta_fov/2.) + (self.theta_fov/self.n_pixels)*(self.I+0.5)
        self.phi_pixels = -(self.phi_fov/2.) + (self.phi_fov/self.m_pixels)*(self.J+0.5)



    def get_vector_for_all_pixels(self):
        '''This will calculate a unit vector in xyz in camera coords for each pixel 
        using its theta and phi. 
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
        '''This function will rotate the camera to the correct viewing direction, 
        and rotate all the unit vectors for the pixels. '''   

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
        '''This will calculate the coordinates along the LOS for a given 
        coordinate spacing.'''

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
    
        self.fig = fig 
        
    def plot_LOS_vectors(self, elev=45, azim=45):
        '''This will plot the LOS vectors from the spacecraft position.'''

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # For each pixel: 
        for i in range(self.n_pixels):
            for j in range(self.m_pixels):
                ax.plot(self.xpos[i][j], self.ypos[i][j], self.zpos[i][j], 'k')
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        self.add_earth(ax)
        ax.set_aspect('equal') 
        
        ax.view_init(elev,azim) 
        
        self.fig2 = fig
       
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
