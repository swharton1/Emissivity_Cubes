import numpy as np
from spacepy import pycdf
import os
import json
from matplotlib.patches import Wedge, Polygon, Circle
import matplotlib.pyplot as plt

from . import get_meridians as gm 

class openggcm(): 
    '''This class will contain all the code to read the openggcm MHD files and calculate the emissivity.
    '''
    def __init__(self, filename="OpenGGCM_sim_n05_vx400vy30vz00_bx00by00bz05.cdf", print_info=False, read_metadata=True, xmin=-10, xmax=None, ymin=-40, ymax=40, zmin=-40, zmax=40):
        '''This function will open a cdf file and return the contents as 
        a dictionary. 
        
        Parameters
        ----------
        filename - CDF filename of the OpenGGCM simulation
        print_info - Boolean to print the dictionary information in the CDF file 
        read_metadata - Boolean to read in the metadata from the json file of the same name. 
        xmin - Minimum x coordinate. def = -10
        xmax - Maximum x coordinate. def = None
        ymin - Minimum y coordinate. def = None
        ymax - Maximum y coordinate. def = None
        zmin - Minimum z coordinate. def = None
        zmax - Maximum z coordinate. def = None 
        
        '''

        self.filename = filename 

        #Get paths. 
        self.openggcm_path = os.environ.get("OPENGGCM_PATH") 
        self.plot_path = os.environ.get("PLOT_PATH") 
        
        self.variable_units = {}
        
        print ("Get CDF file...") 
        with pycdf.CDF(os.path.join(self.openggcm_path,filename)) as cdf:

            #This just creates a list of variable names. 
            self.variable_names = list(cdf)
            
            for v in self.variable_names:
                variable = cdf[v]
                if print_info:
                    print(f"\nVariable: {v}")
                for attr_name, attr_value in variable.attrs.items():
                    if attr_name == "units":
                        self.variable_units[v] = attr_value 
                    if print_info:
                        print (f" {attr_name}:{attr_value}")
                
            
                
                
            #Add information to a normal dictionary as 
            #this format is annoying. 
            self.data = {} 
            for k in cdf.keys():
                self.data[k] = np.array(cdf[k][...][0]).squeeze()
        
        #Get the metadata if selected. 
        if read_metadata:        
            self.read_metadata() 
        
        #Get 3D arrays 
        self.get_3d_arrays()
        
        #Apply coordinate limits on simulation. 
        self.apply_limit(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax) 
        
            
    def read_metadata(self, metafile=None):
        '''This will read the json metadata file that goes with the simulation. It should have the same name format as the cdf file. '''
        
        #Get name of file. 
        print ("Get metadata...") 
        if metafile is None: 
            metafile = os.path.join(self.openggcm_path,self.filename.split(".")[0]+".json")
        else:
            metafile = os.path.join(self.openggcm_path,metafile)
        #Check it exists. 
        try: 
            with open(metafile) as f: 
                self.metadata = json.load(f)
                self.input_parameters = self.metadata["inputParameters"] 
                
                #Add the most important SW parameters to self. 
                self.density = float(self.input_parameters["N"])
                self.vx = float(self.input_parameters["Vx"])
                self.vy = float(self.input_parameters["Vy"])
                self.vz = float(self.input_parameters["Vz"])
                self.bx = float(self.input_parameters["Bx"])
                self.by = float(self.input_parameters["By"])
                self.bz = float(self.input_parameters["Bz"]) 
                self.temp = float(self.input_parameters["T"]) 
                self.dipole = float(self.input_parameters["dipole_tilt"]) 
                
        except FileNotFoundError:
            print ("Metadata file not found: {}".format(metafile)) 
    
    def get_3d_arrays(self):
        '''This gets the 3D cartesian arrays for x, y and z and also reshapes 
        all data products on grid system 1. '''
        
        print ('Reshape arrays...')
        #Get X, Y and Z data in 3D. 
        YY, ZZ, XX = np.meshgrid(self.data['y'], self.data['z'], self.data['x'])
        #Add ZZ, YY, and XX to dictionary. 
        
        #The OpenGGCM coordinate system is not GSM, but spun backwards. x and y values need multiplying by -1 to switch them to GSM. 
        # Put -1 because this simulation has x pointing the opposite way to PPMLR!!! 
        self.x_3d = XX*(-1)
        self.y_3d = YY*(-1)
        self.z_3d = ZZ
        self.data['x'] = -self.data['x'] 
        self.data['ux'] = -self.data['ux']
        self.data['bx'] = -self.data['bx']
        self.data['jx'] = -self.data['jx']
        self.data['y'] = -self.data['y'] 
        self.data['uy'] = -self.data['uy']
        self.data['by'] = -self.data['by']
        self.data['jy'] = -self.data['jy']
        
        #Reshape all the other arrays to 3D. 
        self.data['bx'] = self.data['bx'].reshape(self.data['z'].size, self.data['y'].size, self.data['x'].size) 
        self.data['by'] = self.data['by'].reshape(self.data['z'].size, self.data['y'].size, self.data['x'].size)
        self.data['bz'] = self.data['bz'].reshape(self.data['z'].size, self.data['y'].size, self.data['x'].size)
        self.data['ux'] = self.data['ux'].reshape(self.data['z'].size, self.data['y'].size, self.data['x'].size) 
        self.data['uy'] = self.data['uy'].reshape(self.data['z'].size, self.data['y'].size, self.data['x'].size)
        self.data['uz'] = self.data['uz'].reshape(self.data['z'].size, self.data['y'].size, self.data['x'].size)
        self.data['jx'] = self.data['jx'].reshape(self.data['z'].size, self.data['y'].size, self.data['x'].size) 
        self.data['jy'] = self.data['jy'].reshape(self.data['z'].size, self.data['y'].size, self.data['x'].size)
        self.data['jz'] = self.data['jz'].reshape(self.data['z'].size, self.data['y'].size, self.data['x'].size)
        
        #Resistivity (eta), atomic mass density (rho) and pressure (p) 
        self.data['eta'] = self.data['eta'].reshape(self.data['z'].size, self.data['y'].size, self.data['x'].size)
        self.data['rho'] = self.data['rho'].reshape(self.data['z'].size, self.data['y'].size, self.data['x'].size)
        self.data['p'] = self.data['p'].reshape(self.data['z'].size, self.data['y'].size, self.data['x'].size)
                
    def apply_limit(self, xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None):
        '''This will apply a limit to the x, y or z data to exclude some.'''

        print ('Apply coordinate limits...') 
        
        if xmin is not None: 
            i = np.where(self.x_3d[0][0] >= xmin) 
            self.x_3d = self.x_3d[:,:,i[0]]
            self.y_3d = self.y_3d[:,:,i[0]]    
            self.z_3d = self.z_3d[:,:,i[0]]
            self.data['bx'] = self.data['bx'][:,:,i[0]]
            self.data['by'] = self.data['by'][:,:,i[0]]
            self.data['bz'] = self.data['bz'][:,:,i[0]]
            self.data['ux'] = self.data['ux'][:,:,i[0]]
            self.data['uy'] = self.data['uy'][:,:,i[0]]
            self.data['uz'] = self.data['uz'][:,:,i[0]]
            self.data['jx'] = self.data['jx'][:,:,i[0]]
            self.data['jy'] = self.data['jy'][:,:,i[0]]
            self.data['jz'] = self.data['jz'][:,:,i[0]]
            self.data['eta'] = self.data['eta'][:,:,i[0]]
            self.data['rho'] = self.data['rho'][:,:,i[0]]
            self.data['p'] = self.data['p'][:,:,i[0]]
            
            j = np.where(self.data['x'] >= xmin)
            self.data['x'] = self.data['x'][j] 
            

        if xmax is not None: 
            i = np.where(self.x_3d[0][0] < xmax) 
            self.x_3d = self.x_3d[:,:,i[0]]
            self.y_3d = self.y_3d[:,:,i[0]]
            self.z_3d = self.z_3d[:,:,i[0]]
            self.data['bx'] = self.data['bx'][:,:,i[0]]
            self.data['by'] = self.data['by'][:,:,i[0]]
            self.data['bz'] = self.data['bz'][:,:,i[0]]
            self.data['ux'] = self.data['ux'][:,:,i[0]]
            self.data['uy'] = self.data['uy'][:,:,i[0]]
            self.data['uz'] = self.data['uz'][:,:,i[0]]
            self.data['jx'] = self.data['jx'][:,:,i[0]]
            self.data['jy'] = self.data['jy'][:,:,i[0]]
            self.data['jz'] = self.data['jz'][:,:,i[0]]
            self.data['eta'] = self.data['eta'][:,:,i[0]]
            self.data['rho'] = self.data['rho'][:,:,i[0]]
            self.data['p'] = self.data['p'][:,:,i[0]]

            j = np.where(self.data['x'] < xmax)
            self.data['x'] = self.data['x'][j] 
            
        if ymin is not None: 
            j = np.where(self.y_3d[0,:,0] >= ymin)
            self.x_3d = self.x_3d[:,j[0],:]
            self.y_3d = self.y_3d[:,j[0],:]
            self.z_3d = self.z_3d[:,j[0],:]
            self.data['bx'] = self.data['bx'][:,j[0],:]
            self.data['by'] = self.data['by'][:,j[0],:]
            self.data['bz'] = self.data['bz'][:,j[0],:]
            self.data['ux'] = self.data['ux'][:,j[0],:]
            self.data['uy'] = self.data['uy'][:,j[0],:]
            self.data['uz'] = self.data['uz'][:,j[0],:]
            self.data['jx'] = self.data['jx'][:,j[0],:]
            self.data['jy'] = self.data['jy'][:,j[0],:]
            self.data['jz'] = self.data['jz'][:,j[0],:]
            self.data['eta'] = self.data['eta'][:,j[0],:]
            self.data['rho'] = self.data['rho'][:,j[0],:]
            self.data['p'] = self.data['p'][:,j[0],:]

            j = np.where(self.data['y'] >= ymin)
            self.data['y'] = self.data['y'][j] 
            
        if ymax is not None: 
            j = np.where(self.y_3d[0,:,0] < ymax)
            self.x_3d = self.x_3d[:,j[0],:]
            self.y_3d = self.y_3d[:,j[0],:]
            self.z_3d = self.z_3d[:,j[0],:]
            self.data['bx'] = self.data['bx'][:,j[0],:]
            self.data['by'] = self.data['by'][:,j[0],:]
            self.data['bz'] = self.data['bz'][:,j[0],:]
            self.data['ux'] = self.data['ux'][:,j[0],:]
            self.data['uy'] = self.data['uy'][:,j[0],:]
            self.data['uz'] = self.data['uz'][:,j[0],:]
            self.data['jx'] = self.data['jx'][:,j[0],:]
            self.data['jy'] = self.data['jy'][:,j[0],:]
            self.data['jz'] = self.data['jz'][:,j[0],:]
            self.data['eta'] = self.data['eta'][:,j[0],:]
            self.data['rho'] = self.data['rho'][:,j[0],:]
            self.data['p'] = self.data['p'][:,j[0],:]

            j = np.where(self.data['y'] < ymax)
            self.data['y'] = self.data['y'][j] 
            
        if zmin is not None:
            k = np.where(self.z_3d[:,0,0] >= zmin)
            self.x_3d = self.x_3d[k[0],:,:]
            self.y_3d = self.y_3d[k[0],:,:]
            self.z_3d = self.z_3d[k[0],:,:]
            self.data['bx'] = self.data['bx'][k[0],:,:]
            self.data['by'] = self.data['by'][k[0],:,:]
            self.data['bz'] = self.data['bz'][k[0],:,:]
            self.data['ux'] = self.data['ux'][k[0],:,:]
            self.data['uy'] = self.data['uy'][k[0],:,:]
            self.data['uz'] = self.data['uz'][k[0],:,:]
            self.data['jx'] = self.data['jx'][k[0],:,:]
            self.data['jy'] = self.data['jy'][k[0],:,:]
            self.data['jz'] = self.data['jz'][k[0],:,:]
            self.data['eta'] = self.data['eta'][k[0],:,:]
            self.data['rho'] = self.data['rho'][k[0],:,:]
            self.data['p'] = self.data['p'][k[0],:,:]

            j = np.where(self.data['z'] >= zmin)
            self.data['z'] = self.data['z'][j] 
            
        if zmax is not None:
            k = np.where(self.z_3d[:,0,0] < zmax)
            self.x_3d = self.x_3d[k[0],:,:]
            self.y_3d = self.y_3d[k[0],:,:]
            self.z_3d = self.z_3d[k[0],:,:]
            self.data['bx'] = self.data['bx'][k[0],:,:]
            self.data['by'] = self.data['by'][k[0],:,:]
            self.data['bz'] = self.data['bz'][k[0],:,:]
            self.data['ux'] = self.data['ux'][k[0],:,:]
            self.data['uy'] = self.data['uy'][k[0],:,:]
            self.data['uz'] = self.data['uz'][k[0],:,:]
            self.data['jx'] = self.data['jx'][k[0],:,:]
            self.data['jy'] = self.data['jy'][k[0],:,:]
            self.data['jz'] = self.data['jz'][k[0],:,:]
            self.data['eta'] = self.data['eta'][k[0],:,:]
            self.data['rho'] = self.data['rho'][k[0],:,:]
            self.data['p'] = self.data['p'][k[0],:,:]
    
            j = np.where(self.data['z'] < zmax)
            self.data['z'] = self.data['z'][j] 
    
    def plot_both_planes_variable(self, variable="bx", cmap="hot", levels=100, save=False, savetag=""):
        '''This will plot the effective velocity, proton density and neutral hydrogen
         density in the X-Z and X-Y planes side by side. 
        '''

        #Get meridian data for veff. 
        xp_y, yp_y, zp_y, var_y, xp_z, yp_z, zp_z, var_z = gm.calculate_meridian_planes(self.x_3d, self.y_3d, self.z_3d, self.data[variable])

        # Create a filename label so you know which file you plotted. 
        file_label = self.filename.split("/")[-1]
        file_label = file_label.split(".")[0]

        #Create a label so you can see which simulation this is. 
        label = "n = {} cm".format(self.density)+r"$^{-3}$"+"   v = ({}, {}, {}) km/s   B = ({}, {}, {}) nT".format(self.vx, self.vy, self.vz, self.bx, self.by, self.bz) 
        

        # Now you can make the contour plot. 
        fig = plt.figure(figsize=(6,4))
        fig.subplots_adjust(wspace=0.5)

        #Add title label 
        fig.text(0.5, 0.95, label, ha='center', fontsize=10) 
        
        # Effective velocity 
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)

         # Get contour levels. 
        #levels_veff = np.linspace(vmin_veff, vmax_veff, levels+1)
        
        # Constant y. 
        cont1 = ax1.contourf(xp_y, zp_y, var_y, cmap='hot')
        ax1.set_xlabel('X [RE]')
        ax1.set_ylabel('Z [RE]')
        ax1.set_title("{}: XZ Plane".format(variable), fontsize=10)
        ax1.set_aspect("equal")

        # Add sketch of the Earth on top. 
        self.make_earth(ax1, rotation=-90)

         # Constant z. 
        cont2 = ax2.contourf(xp_z, yp_z, var_z, cmap='hot')
        ax2.set_xlabel('X [RE]')
        ax2.set_ylabel('Y [RE]')
        ax2.set_title("{}: XY Plane".format(variable), fontsize=10)
        ax2.set_aspect("equal")

        # Add colourbar and earth. 
        
        cbar = plt.colorbar(cont2, ax=ax2, shrink=0.5)
        cbar.set_label("{} ({})".format(variable, self.variable_units[variable]))
        
        #level_min = np.ceil(cont2.levels.min())
        #level_max = np.floor(cont2.levels.max())
        level_min = cont2.levels.min()
        level_max = cont2.levels.max()
        
        print (level_min, level_max)
        cticks = np.linspace(level_min, level_max,5)
        print(cticks)
        cbar.set_ticks(cticks)
        # cbar.set_ticklabels([r'$10^{'+str(i)+'}$' for i in cticks])

        # Add sketch of the Earth on top. 
        self.make_earth(ax2, rotation=-90)


        # Option to save figure. 
        if save: 
            fig.savefig(self.plot_path+"/{}_{}_{}.png".format(file_label, variable, savetag))

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
    
    
