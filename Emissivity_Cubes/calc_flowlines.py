#This is my new attempt to do the flowlines with some advice from Zac. 
#This is updated not to keep the massive arrays describing all the flowlines, 
#because that crashed the terminal in masking2.py! 

import numpy as np 
from scipy.interpolate import interpn 
import os
import matplotlib.pyplot as plt
from time import process_time

class flowlines():
    '''This will try to be a replacement for the really slow function in masking.py. 
    
    '''
    
    def __init__(self, openggcm=None, step=0.5, n_steps=2, pd=1):
        '''Takes in initial parameters
        
        Parameters
        ----------
        openggcm - openggcm object from read_openggcm
        step - step size to move along the flowline in. Smaller = accurate. 
        
        '''
        
        self.openggcm = openggcm
        self.step = step 
        self.n_steps = n_steps 
        self.pd = pd
        self.plot_path = os.environ.get('PLOT_PATH')+'masking/'
        
        #Get cube limits.  
        self.get_cube_limits() 
        
        #You will need this for interpn. 
        self.points_original = (self.openggcm.data['z'], self.openggcm.data['y'], self.openggcm.data['x']) 
        
        #Setup arrays for first step.  
        self.setup_initial_step(pd=pd) 
        
        #Get all the next steps. 
        self.take_remaining_steps()
        
        #Assign an origin to each point. 
        self.get_origin() 

    def get_cube_limits(self):
        '''This works out a uniform grid to use'''
        
        #Get minima and maxima in x, y and z directions.
        self.xmin = self.openggcm.data['x'].min() 
        self.xmax = self.openggcm.data['x'].max() 
        self.ymin = self.openggcm.data['y'].min() 
        self.ymax = self.openggcm.data['y'].max() 
        self.zmin = self.openggcm.data['z'].min() 
        self.zmax = self.openggcm.data['z'].max() 
        
               
    def setup_initial_step(self, pd):
        '''This will take the initial step. It does not need the uniform cube for any calculations. 
        
        Parameters
        ----------
        pd - Point density. It will only calculate for every pd-th point. A value of 1 will do all values in the cube. Set to a higher number for testing on a whole cube. 
        
        '''
        print ('Take initial step...')    
        
        #Add initial positions to array. 
        self.px_all = np.copy(self.openggcm.x_3d)[::pd,::pd,::pd]
        self.py_all = np.copy(self.openggcm.y_3d)[::pd,::pd,::pd]
        self.pz_all = np.copy(self.openggcm.z_3d)[::pd,::pd,::pd] 
        
        #Extract velocity information. 
        self.vx_all = np.copy(self.openggcm.data['ux'])[::pd,::pd,::pd] 
        self.vy_all = np.copy(self.openggcm.data['uy'])[::pd,::pd,::pd] 
        self.vz_all = np.copy(self.openggcm.data['uz'])[::pd,::pd,::pd] 
        
        #Get velocity magnitude. 
        self.vmag_all = np.sqrt(self.vx_all**2 + self.vy_all**2 + self.vz_all**2) 
        
        #Setup final position arrays. Initialise with original values in case you don't get any further for some points. 
        self.px_final = np.copy(self.px_all) 
        self.py_final = np.copy(self.py_all) 
        self.pz_final = np.copy(self.pz_all)  
        
        #Now create an array of ones and zeros that states whether 
        #a flowline is solar wind or not. 
        self.origin = np.zeros((self.px_final.shape))
 
    def take_remaining_steps(self):
        '''This will loop through the number of steps to take and work out the remaining positions along each field line.''' 
        print ('Take remaining steps...') 
        ts = process_time() 
        #Loop through positions, but don't include initial. 
        for step in range(0,self.n_steps):
            print (step) 
            
            #Calcualate next position from previous position and previous velocity. 
            pxn = self.px_all - self.step*(self.vx_all/self.vmag_all) 
            pyn = self.py_all - self.step*(self.vy_all/self.vmag_all)
            pzn = self.pz_all - self.step*(self.vz_all/self.vmag_all)
            
            #Now work out which new positions are still inside the cube. 
            inside = np.where((pxn < self.xmax) & (pxn > self.xmin) & (pyn < self.ymax) & (pyn > self.ymin) & (pzn < self.zmax) & (pzn > self.zmin))
            
            #Work out which points are outside the cube. 
            outside = np.where((pxn >= self.xmax) | (pxn <= self.xmin) | (pyn >= self.ymax) | (pyn <= self.ymin) | (pzn >= self.zmax) | (pzn <= self.zmin)) 
            
            #Calculate next velocity at these new positions. 
            new_points = (pzn[inside], pyn[inside], pxn[inside]) 
            
            #Only calculate velocities for positions inside the box. 
            vxn = interpn(self.points_original, self.openggcm.data['ux'], new_points, method='linear')
            vyn = interpn(self.points_original, self.openggcm.data['uy'], new_points, method='linear')
            vzn = interpn(self.points_original, self.openggcm.data['uz'], new_points, method='linear')
            
            #REPLACE POSITION ARRAYS WITH NEW POSITIONS, BUT ONLY NEW POINTS. 
            self.px_all[inside] = pxn[inside]
            self.py_all[inside] = pyn[inside]
            self.pz_all[inside] = pzn[inside]
            
            #POSITIONS OUTSIDE BOX GET FILLED WITH NANS. 
            self.px_all[outside] = np.nan
            self.py_all[outside] = np.nan
            self.pz_all[outside] = np.nan
            
            #REPLACE VELOCITY ARRAYS WITH NEW VELOCITIES, BUT ONLY NEW POINTS. 
            self.vx_all[inside] = vxn
            self.vy_all[inside] = vyn
            self.vz_all[inside] = vzn 
            
            #VELOCITIES OUTSIDE BOX GET FILLED WITH NANS. 
            self.vx_all[outside] = np.nan
            self.vy_all[outside] = np.nan
            self.vz_all[outside] = np.nan
            
            #UPDATE VELOCITY MAGNITUDES. 
            self.vmag_all = np.sqrt(self.vx_all**2 + self.vy_all**2 + self.vz_all**2)
        
            #UPDATE FINAL POSITIONS.
            #Just update those inside the box. 
            #The rest can be left with their previous 
            #values as they have already left the box. 
            self.px_final[inside] = pxn[inside]
            self.py_final[inside] = pyn[inside]
            self.pz_final[inside] = pzn[inside] 
            
        #Check time. 
        te = process_time()
        print ('Time to calc. flowlines = {:} s'.format(te-ts))

    def get_origin(self):
        '''This uses the final position arrays to determine whether a point is
        from the solar wind or not.'''
        
        print ('Identifying Origin...') 
        
        sw = np.where((self.px_final > self.xmax-1) | (self.py_final > self.ymax-0.1) | (self.py_final < self.ymin+0.1) | (self.pz_final > self.zmax-0.1) | (self.pz_final < self.zmin+0.1))
        
        self.origin[sw] = 1 
        
    def plot_origins(self, plot_sw=True, plot_ms=True, s=0.2):
        '''This will make a colour coded plot to show whether plasma is from the solar wind or not.
        
        Parameters
        ----------
        plot_sw - Boolean to plot SW plasma
        plot_ms - Boolean to plot MS plasma 
        point_density - reduces the number of points plotted by this factor. 
        
        '''
        
        #Get pd for x,y,z arrays. 
        pd = self.pd 
        
        plt.close("all")
        fig = plt.figure(figsize=(8,4))
        ax1 = fig.add_subplot(131, projection='3d') 
        ax2 = fig.add_subplot(132, projection='3d') 
        ax3 = fig.add_subplot(133, projection='3d') 
        
        #Set all labels. 
        ax1.set(xlabel='x', ylabel='y', zlabel='z', xbound=[self.xmin,self.xmax], ybound=[self.ymin,self.ymax], zbound=[self.zmin,self.zmax]) 
        ax2.set(xlabel='x', ylabel='y', zlabel='z', xbound=[self.xmin,self.xmax], ybound=[self.ymin,self.ymax], zbound=[self.zmin,self.zmax]) 
        ax3.set(xlabel='x', ylabel='y', zlabel='z', xbound=[self.xmin,self.xmax], ybound=[self.ymin,self.ymax], zbound=[self.zmin,self.zmax]) 
        
        #Separate solar wind and magnetospheric scatter. 
        sw = np.where(self.origin == 1)
        ms = np.where(self.origin == 0) 
        
        #Plot data as coloured scatter. 
        if plot_sw:
            ax1.scatter(self.openggcm.x_3d[::pd,::pd,::pd][sw], self.openggcm.y_3d[::pd,::pd,::pd][sw], self.openggcm.z_3d[::pd,::pd,::pd][sw], c='r', marker='o', s=s)
            ax2.scatter(self.openggcm.x_3d[::pd,::pd,::pd][sw], self.openggcm.y_3d[::pd,::pd,::pd][sw], self.openggcm.z_3d[::pd,::pd,::pd][sw], c='r', marker='o', s=s)
            ax3.scatter(self.openggcm.x_3d[::pd,::pd,::pd][sw], self.openggcm.y_3d[::pd,::pd,::pd][sw], self.openggcm.z_3d[::pd,::pd,::pd][sw], c='r', marker='o', s=s)
        
        
        if plot_ms:
            ax1.scatter(self.openggcm.x_3d[::pd,::pd,::pd][ms], self.openggcm.y_3d[::pd,::pd,::pd][ms], self.openggcm.z_3d[::pd,::pd,::pd][ms], c='b', marker='o', s=s)
            ax2.scatter(self.openggcm.x_3d[::pd,::pd,::pd][ms], self.openggcm.y_3d[::pd,::pd,::pd][ms], self.openggcm.z_3d[::pd,::pd,::pd][ms], c='b', marker='o', s=s)
            ax3.scatter(self.openggcm.x_3d[::pd,::pd,::pd][ms], self.openggcm.y_3d[::pd,::pd,::pd][ms], self.openggcm.z_3d[::pd,::pd,::pd][ms], c='b', marker='o', s=s)
        
        ax1.set_aspect('equal')
        ax1.view_init(elev=0, azim=0)  
        ax2.set_aspect('equal')
        ax2.view_init(elev=0, azim=90)  
        ax3.set_aspect('equal')
        ax3.view_init(elev=90, azim=0)  
        
        
#    def get_origin(self):
#        '''This will work out if plasma is from the solar wind or not from its final x position.'''
        
#        print ('Identifying Origin...') 
#        ts = process_time()
#        dimen = self.px_all.shape 
        #Loop through all the start points. 
#        for z in range(dimen[1]):
#            for y in range(dimen[2]):
#                for x in range(dimen[3]): 
#                    xpos = self.px_all[:,z,y,x]
#                    realx = np.isfinite(xpos).size 
#                    lastx = xpos[realx-1] 
                    
#                    if lastx > self.xmax-10:
#                        self.origin[z,y,x] = 1 
#                    else:
#                        self.origin[z,y,x] = 0 
#        te = process_time()
#        print ('Time to assign origin = {:} s'.format(te-ts))  
        
                  
#    def plot_flowlines(self, point_density=10, elev=45, azim=45):
#        '''This will make a plot of all the flowlines.'''
#        plt.close("all")
#        fig = plt.figure(figsize=(8,4))
#        ax1 = fig.add_subplot(131, projection='3d') 
#        ax2 = fig.add_subplot(132, projection='3d') 
#        ax3 = fig.add_subplot(133, projection='3d') 
        
#        ax1.set_xlabel('x')
#        ax1.set_ylabel('y')
#        ax1.set_zlabel('z') 
#        ax2.set_xlabel('x')
#        ax2.set_ylabel('y')
#        ax2.set_zlabel('z') 
#        ax3.set_xlabel('x')
#        ax3.set_ylabel('y')
#        ax3.set_zlabel('z') 
        
#        dimen = self.px_all.shape 
        
        #Loop through all the start points. 
#        for z in range(dimen[1]):
#            for y in range(dimen[2]):
#                for x in range(dimen[3]): 
                    #Only plot these points or it gets silly. 
#                    if (z%point_density == 0) and (y%point_density == 0) and (-x%point_density == 0):
#                        print (z,y,x)
                        #Get x, y and z data for this flowline. 
#                        xpos = self.px_all[:,z,y,x] 
#                        ypos = self.py_all[:,z,y,x]
#                        zpos = self.pz_all[:,z,y,x] 
                        
#                        if self.origin[z,y,x] == 1: 
#                            c = 'r' 
#                        else:
#                            c = 'b' 
                        
#                        ax1.plot(xpos, ypos, zpos, color=c) 
#                        ax2.plot(xpos, ypos, zpos, color=c)  
#                        ax3.plot(xpos, ypos, zpos, color=c)   
                        #ax.scatter(xpos[0], ypos[0], zpos[0], marker='x', c='r')    

#        ax1.set_aspect('equal')
#        ax1.view_init(elev=0, azim=0)  
##        ax2.set_aspect('equal')
#        ax2.view_init(elev=0, azim=90)  
#        ax3.set_aspect('equal')
#        ax3.view_init(elev=90, azim=0)  
        
        
        


