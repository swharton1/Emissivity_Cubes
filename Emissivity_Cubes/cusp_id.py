import numpy as np
import matplotlib.pyplot as plt 

import os

class cusps():
    '''This class will contain the methods needed to calculate flowlines for a given point in a simulation.'''
    def __init__(self, openggcm=None):
        '''Takes in initial parameters
        
        Parameters
        ----------
        openggcm - openggcm object from read_openggcm
        
        '''
        
        self.openggcm = openggcm
        
        self.plot_path = os.environ.get('PLOT_PATH')+'masking/'
        
        #Now create an array of ones and zeros that states whether 
        #a flowline is solar wind or not. 
        self.origin = np.zeros((self.openggcm.x_3d.shape))
        self.origin.fill(np.nan)   


    #CUSP IDENTIFICATION USING PRESSURE CONDITIONS. 
    
    
    def get_cusp_mask(self, pthresh=None, r_min=None, hemi='north', figure=True):
        '''This will use the method of Sun et al. (2019) to identify the cusps. Also used by Samsonov et al. (2019a). 
        Parameters
        ----------
        pthresh - limit is fraction of maximum pressure.
        r_min - lowest radial distance for a spherical shell. def = 5RE
        hemi - whether to do the north or south cusp
        figure - boolean to make the full figure showing the details. 
        
        '''
        
        #Set pthresh here based on solar wind speed. 
        if pthresh is None:
            self.pthresh = (0.1/100)*self.openggcm.vx + 1.1 
            print ('pthresh = ', self.pthresh) 
        else:
            self.pthresh = pthresh
        
        #Set r_min here based on solar wind speed. 
        if r_min is None:
            if self.openggcm.vx < 400:
                self.r_min = 5.5
            else:
                self.r_min = 5 
        else: 
            self.r_min = r_min
        
        #Check if you need to convert to spherical. 
        try: 
            print (self.r.shape)
        except AttributeError:
            print ('Getting spherical coords...') 
            self.convert_xyz_to_spherical(self.openggcm.x_3d, self.openggcm.y_3d, self.openggcm.z_3d) 
            
        #First need to define a set of radial shells. 
        dr = 0.5 
        r_sh = np.linspace(self.r_min, self.r_min+6, 13)
        #Colatitude limit 
        if hemi.lower() == 'north':
            colat = np.deg2rad(60)
        elif hemi.lower() == 'south':
            colat = np.deg2rad(120)
        else:
            raise ValueError("Invalid hemisphere chosen: 'north' or 'south'")
        
        self.latitude = np.pi/2 - self.theta
        self.local_time = (self.phi+np.pi)*(12/np.pi)
        
        if figure:
            
            #Make figure to put data on. 
            fig = plt.figure(figsize=(8,10))
            fig.subplots_adjust(hspace=0.5)
        
        #Loop through each shell. 
        for r, rval in enumerate(r_sh):
            print (r, rval)
            
            #Select all values inside shell. 
            if hemi == 'north':
                i = np.where((self.r >= rval) & (self.r < rval+dr) & (self.theta <= colat) & (self.phi >= -np.pi/4) & (self.phi <= np.pi/4))
            else:
                i = np.where((self.r >= rval) & (self.r < rval+dr) & (self.theta >= colat) & (self.phi >= -np.pi/4) & (self.phi <= np.pi/4))
                
            
            #Extract this data. 
            x_sh = self.openggcm.x_3d[i]
            y_sh = self.openggcm.y_3d[i]
            z_sh = self.openggcm.z_3d[i]
            theta_sh = self.theta[i]
            phi_sh = self.phi[i]
            lat_sh = self.latitude[i]
            lt_sh = self.local_time[i]
            or_sh = self.origin[i] 
            pr_sh = self.openggcm.data['p'][i]
            
            lt_sh = np.squeeze(lt_sh)
            lat_sh = np.squeeze(lat_sh) 
            pr_sh = np.squeeze(pr_sh) 
            
            #Add axes for this shell. 
            if figure:
                ax = fig.add_subplot(len(r_sh),4,1+4*r, projection='3d')
                ax.scatter(x_sh, y_sh, z_sh, c=pr_sh)
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_zlabel('z')
                ax.set_xlim(-14,14)
                ax.set_ylim(-14,14)
                if hemi == 'north': ax.set_zlim(0,14)
                else: ax.set_zlim(-14,0)
                ax.set_aspect('equal')
                ax.view_init(elev=30, azim=60)
            
                ax2 = fig.add_subplot(len(r_sh),4,2+4*r) 
                ax2.scatter(lt_sh, np.rad2deg(lat_sh), c=pr_sh)
                if r == len(r_sh)-1: ax2.set_xlabel('Local Time')
                ax2.set_ylabel('latitude (deg)') 
                ax2.set_title('r = {} RE shell'.format(rval), fontsize=10)
            
                #Try to work out an actual contour map. 
            
            
                ax3 = fig.add_subplot(len(r_sh),4,3+4*r)
                tricont = ax3.tricontourf(lt_sh, np.rad2deg(lat_sh), pr_sh)
                if r == len(r_sh)-1: ax3.set_xlabel('Local Time')
                ax3.set_title('r = {} RE shell'.format(rval), fontsize=10)
            
            else:
                tricont = plt.tricontour(lt_sh, np.rad2deg(lat_sh), pr_sh)
                self.tricont = tricont
                
            #Get maximum p value in contour map. 
            psh_max = tricont.get_clim()[1] 
            print ('Max Shell Pressure: ', psh_max)
            j = np.where(pr_sh > self.pthresh*psh_max)
            
            if figure:
                #Plot only value above this pressure threshold. 
                ax4 = fig.add_subplot(len(r_sh),4,4+4*r) 
                tricont2 = ax4.tricontourf(lt_sh[j], np.rad2deg(lat_sh[j]), pr_sh[j], levels=tricont.levels)
                if r == len(r_sh)-1: ax4.set_xlabel('Local Time')
                ax4.set_title('r = {} RE shell'.format(rval), fontsize=10)
                ax4.set_xlim(ax3.get_xlim())
                ax4.set_ylim(ax3.get_ylim())
            
                #self.tricont = tricont
                #self.tricont2 = tricont2
            
            #Set these values to 1. They are inside the cusp. 
            or_sh[j] = 1
            self.origin[i] = or_sh
            
            
        if figure: 
            print ('Making cusp figure:')
            print (self.plot_path+'pressure_contours_{}_{}_{}.png'.format(hemi, self.openggcm.vx, pthresh))
            plt.close("all")
            fig.savefig(self.plot_path+'pressure_contours_{}_{}_{}.png'.format(hemi, self.openggcm.vx, pthresh))
        
            
    def convert_xyz_to_spherical(self, x, y, z):
        '''This will convert the data to spherical coords.'''
        
        self.r = np.sqrt(x**2 + y**2 + z**2)
        
        self.theta = np.arccos(z/self.r)
        
        self.phi = np.arctan2(y,x)
