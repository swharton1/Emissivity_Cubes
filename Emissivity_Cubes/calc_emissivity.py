import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Polygon, Circle
from astropy.io import fits as pyfits 

#from . import masking 
from . import calc_flowlines as flowlines
from . import cusp_id 
from SXI_Core import get_meridians as gm
from SXI_Core import get_earth 
from time import process_time

class emissivity():
    '''This class will take the openggcm object and try to calculate the emissivity. It may contain several methods to do this.''' 
    
    def __init__(self, openggcm):
        '''This takes in the initial parameters.
        
        Parameters
        ----------
        openggcm - openggcm object from read_openggcm. '''
        
        self.openggcm = openggcm 
        
        #Attach common data variables. Leave everything else inside self.openggcm 
        self.x_3d = openggcm.x_3d
        self.y_3d = openggcm.y_3d 
        self.z_3d = openggcm.z_3d 
        self.data = openggcm.data 
        
        self.openggcm_path = os.environ.get("OPENGGCM_PATH") 
        self.plot_path = os.environ.get("PLOT_PATH") 
        
    
    def calc_emissivity(self, alpha=1e-15, masking_method='None', mask_proton=0.2, mask_radius=5, point_density=1, n_steps=4):
        '''This function will calculate the emissivity at all points in space. . 
        
        Parameters
        ----------
        alpha - value of alpha you are using. def = 1e-15. 
        masking_method - 'simple' (def) uses simple conditions 
                        based on plasma conditions or radial position. 
                        - 'flowlines' uses the flowline method in the 
                        masking file. Takes a lot longer. 
                        - 'cusps' only keeps the cusps from the cusp pressure method. 
                        - 'flowlines_cusps' uses the flowline method and also the cusp 
                        identification method. 
                        - 'flow_cusps2' Uses new method to calculate flowlines. 
        mask_proton - def = 0.2. Any coordinate with a proton density less than mask_proton of the max density 
                    will have its emissivity set to zero. 
        mask_radius - def = 5. Any coordinate with a radius less than mask_proton will have its 
                    emissivity set to zero. 
        point_density - How often along each axis to calculate a flowline in the flowline method. 
        n_steps - Number of steps backwards along flowline to calculate in new method. 
        
        '''
        self.alpha = alpha
        self.masking_method = masking_method.lower() 
        self.point_density = point_density 
        
        #First, calculate the average velocity. 
        vsw = (self.data['ux']**2 + self.data['uy']**2 + self.data['uz']**2)**0.5
        # In m/s. 
        self.data['vsw'] = vsw*1000
        
        k = 1.38*10**-23
        mp = 1.672*10**-27 
        
        p_pa = self.data['p']*1e-9
        n_m3 = self.data['rho']*1e6
        
        #Thermal velocity in m/s. 
        self.data['vt'] = ((3*p_pa)/(n_m3*mp))**0.5
        
        #Calculate v effective in m/s. 
        veff = (self.data['vsw']**2 + self.data['vt']**2)**0.5 
        
        #Convert to cm/s. 
        self.data['veff'] = veff*100 

        

        #Next, get the neutral hydrogen density. 
        
        #Get r coordinate. 
        RR = (self.x_3d**2 + self.y_3d**2 + self.z_3d**2)**0.5
        
        #Get density of neutral hydrogen with Hodges' model. (3D array)  
        nH = 25*(10/RR)**3 
        self.data['nH'] = nH

        #Make empty emissivity array. 
        self.emis_3d = np.zeros(self.data['veff'].shape)
        #self.emis_3d.fill(np.nan)
        
        # Add in a simple mask based on proton density. 
        if self.masking_method == 'simple': 
            mask1 = np.where((self.data['rho'] > mask_proton*self.data['rho'].max()) & (RR > mask_radius))
            #Apply the mask. 
            self.emis_3d[mask1] = self.alpha*self.data['veff'][mask1]*self.data['rho'][mask1]*self.data['nH'][mask1] 
        
        elif self.masking_method == 'flowlines': 
            flow = masking.flowlines(openggcm=self.openggcm, step='adaptive')
            ts = process_time()
            flow.calc_all_flowlines(point_density=point_density, max_flow_points=5000, r_inner=mask_radius)
            te = process_time()
            print ('Time to calc. flowlines: {}s'.format(te-ts))
            #if point_density > 1: flow.interpolate_origin()
            
            #Add origin to data. 
            self.data['origin'] = flow.origin
            
            mask1 = np.where(flow.origin == 1) 
            self.mask1 = mask1
            #Apply the mask. 
            self.emis_3d[mask1] = self.alpha*self.data['veff'][mask1]*self.data['rho'][mask1]*self.data['nH'][mask1] 
        
        elif self.masking_method == 'cusps':
            #Get cusps. 
            cusp_n = cusp_id.cusps(openggcm=self.openggcm)
            cusp_n.get_cusp_mask(pthresh=None, r_min=None, hemi='north', figure=False)
            cusp_s = cusp_id.cusps(openggcm=self.openggcm)
            cusp_s.get_cusp_mask(pthresh=None, r_min=None, hemi='south', figure=False)
            
            #Add cusps together. 
            origin = cusp_n.origin 
            isouth = np.where(cusp_s.origin == 1) 
            origin[isouth] = 1
            
            mask1 = np.where(origin == 1) 
            print (mask1)
            #Apply the mask. 
            self.emis_3d[mask1] = self.alpha*self.data['veff'][mask1]*self.data['rho'][mask1]*self.data['nH'][mask1] 
        
        
        elif self.masking_method == 'flow_cusps2':
            #Get flowlines and assign origins based on this. 
            flow = flowlines.flowlines(openggcm=self.openggcm, step=0.1, n_steps=n_steps, pd=1)
            
            #Get cusps. 
            cusp_n = cusp_id.cusps(openggcm=self.openggcm)
            cusp_n.get_cusp_mask(pthresh=None, r_min=None, hemi='north', figure=False)
            cusp_s = cusp_id.cusps(openggcm=self.openggcm)
            cusp_s.get_cusp_mask(pthresh=None, r_min=None, hemi='south', figure=False)
            
            #Put together origin arrays from flowlines and cusps. 
            origin = flow.origin 
            inorth = np.where(cusp_n.origin == 1) 
            isouth = np.where(cusp_s.origin == 1) 
            
            #Add cusp points in to main origin array. 
            origin[inorth] = 1
            origin[isouth] = 1 
            
            #Create the mask. 
            mask1 = np.where(origin == 1) 
            
            #Apply the mask. 
            self.emis_3d[mask1] = self.alpha*self.data['veff'][mask1]*self.data['rho'][mask1]*self.data['nH'][mask1] 
            
        
              
        elif self.masking_method == 'none': 
            #No mask. 
            self.emis_3d = self.alpha*self.data['veff']*self.data['rho']*self.data['nH']
        else:
            raise ValueError("No valid masking method chosen. Options are: 'simple', 'flowlines' or 'none'.") 
        
        
        self.n = (self.data['x'].size, self.data['y'].size, self.data['z'].size)

        # Add veff, rho and nH to the object as well, as you may wish to plot them too. 
        # Store in km/s 
        self.veff_3d = self.data['veff']/1e5
        self.rho_3d = self.data['rho']
        self.nH_3d = self.data['nH']

    def write_fits_file(self):
        '''This will write all the emissivity data to a FITS file in the same format as PPMLR'''
        
        print ("Writing FITS file for emissivity cube...") 
        #Sort out the filename. 
        self.filename_fits = self.openggcm_path+"OG_{}_n{:02d}_vx{:03d}vy{:02d}vz{:02d}_bx{:02d}by{:02d}bz{:02d}.fits".format(self.masking_method, int(self.openggcm.density), int(self.openggcm.vx), int(self.openggcm.vy), int(self.openggcm.vz), int(self.openggcm.bx), int(self.openggcm.by), int(self.openggcm.bz)) 
        
        #Create a new HDU object. 
        self.hdu = pyfits.PrimaryHDU()
        
        #Add emissivity data. 
        self.hdu.data = self.emis_3d
        self.hdu.header 
        
        #Add the solar wind information. 
        self.hdu.header['bx'] = self.openggcm.bx
        self.hdu.header['by'] = self.openggcm.by
        self.hdu.header['bz'] = self.openggcm.bz
        self.hdu.header['vx'] = self.openggcm.vx
        self.hdu.header['vy'] = self.openggcm.vy
        self.hdu.header['vz'] = self.openggcm.vz
        self.hdu.header['density'] = self.openggcm.density
        self.hdu.header['temp'] = self.openggcm.temp 
        
        #Add HDUs for x, y and z arrays. 
        self.hdux = pyfits.ImageHDU(data=self.openggcm.data['x'], name='x')
        self.hduy = pyfits.ImageHDU(data=self.openggcm.data['y'], name='y')
        self.hduz = pyfits.ImageHDU(data=self.openggcm.data['z'], name='z')
        
        #Create a HDU list. 
        self.hdul = pyfits.HDUList([self.hdu, self.hdux, self.hduy, self.hduz]) 
        
        #Write the fits file. 
        self.hdul.writeto(self.filename_fits, overwrite=True) 
        print ('Created: {}'.format(self.filename_fits)) 
        
        
    def write_ascii_file(self):
        '''This will write all the emissivity data to an ascii file in the same format as PPMLR''' 
        
        print ("Writing ASCII file for emissivity cube...") 
        
        #Sort out the filename. 
        filename = self.openggcm_path+"OpenGGCM_emissivity_{}_n{:02d}_vx{:03d}vy{:02d}vz{:02d}_bx{:02d}by{:02d}bz{:02d}.txt".format(self.masking_method, int(self.openggcm.density), int(self.openggcm.vx), int(self.openggcm.vy), int(self.openggcm.vz), int(self.openggcm.bx), int(self.openggcm.by), int(self.openggcm.bz)) 
        
        with open(filename, 'w') as f: 
            #Start writing things in the file. 
            f.write("Emissivity Cube derived from OpenGGCM with {} masking method: alpha = {}\n".format(self.masking_method, self.alpha))
            f.write("{}\n".format(self.openggcm.temp))
            f.write("{}\n".format(self.openggcm.density))
            f.write("{}\n".format(self.openggcm.vx))
            f.write("{}\n".format(self.openggcm.vy))
            f.write("{}\n".format(self.openggcm.vz))
            f.write("{}\n".format(self.openggcm.bx))
            f.write("{}\n".format(self.openggcm.by))
            f.write("{}\n".format(self.openggcm.bz)) 
            
            #Next add on the lengths of the x, y and z arrays. 
            f.write("{}    {}    {}\n".format(self.data['x'].size, self.data['y'].size, self.data['z'].size)) 
    
            #Next, add the x array. 
            for x, xval in enumerate(self.data['x']):
                if x%6 == 5:
                    f.write("{:.5f}\n".format(xval))
                else:
                    if x == self.data['x'].size-1:
                        f.write("{:.5f}\n".format(xval))
                    else:
                        f.write("{:.5f}    ".format(xval))
            #Next, add the y array. 
            for y, yval in enumerate(self.data['y']):
                if y%6 == 5:
                    f.write("{:.5f}\n".format(yval))
                else:
                    if y == self.data['y'].size-1:
                        f.write("{:.5f}\n".format(yval))
                    else:
                        f.write("{:.5f}    ".format(yval))
            #Next, add the z array. 
            for z, zval in enumerate(self.data['z']):
                if z%6 == 5:
                    f.write("{:.5f}\n".format(zval))
                else:
                    if z == self.data['z'].size-1:
                        f.write("{:.5f}\n".format(zval))
                    else:
                        f.write("{:.5f}    ".format(zval))
            
            #Next flatten 3D emissivity array and read into file. 
            eflat = self.emis_3d.flatten() 
            for p, pval in enumerate(eflat):
                f.write("{:.10f}\n".format(pval)) 
                        
                        
    def plot_eta_3d(self, elev=45, azim=45, markersize=1):
        '''Plot eta points in 3D space with eta being the colour.'''
        
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        #Get all the real origin points that are zero. 
        zero = np.where(self.data['origin'] == 0)
        ones = np.where(self.data['origin'] == 1)
        
        #ax.scatter3D(self.openggcm.x_3d[ones], self.openggcm.y_3d[ones], self.openggcm.z_3d[ones], c=np.log10(self.emis_3d[ones]), alpha=1, s=markersize)
        ax.scatter3D(self.openggcm.x_3d[zero], self.openggcm.y_3d[zero], self.openggcm.z_3d[zero], c='k', alpha=1, s=markersize)
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z') 
        ax.set_aspect('equal')
        ax.view_init(elev=elev, azim=azim)
        
        #Save figure. 
        print (self.plot_path+'3d_check.png')
        fig.savefig(self.plot_path+'3d_check.png')                    
                        
            
    def plot_both_planes_eta(self, cmap="hot", levels=100, vmin=-8, vmax=-4, save=False, savetag=""):
        '''This will plot the emissivity in the X-Z and X-Y planes side by side. 
        '''

        #Get meridian data. 
        xp_y, yp_y, zp_y, etad_y, xp_z, yp_z, zp_z, etad_z = gm.calculate_meridian_planes(self.x_3d, self.y_3d, self.z_3d, self.emis_3d)
        
        
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
        file_label = self.openggcm.filename.split("/")[-1]
        file_label = file_label.split(".")[0]

        #Create a label so you can see which simulation this is. 
        label = "n = {} cm".format(self.openggcm.density)+r"$^{-3}$"+"   v = ({}, {}, {}) km/s   B = ({}, {}, {}) nT".format(self.openggcm.vx, self.openggcm.vy, self.openggcm.vz, self.openggcm.bx, self.openggcm.by, self.openggcm.bz) 
        
         # Get contour levels. 
        levels = np.linspace(vmin, vmax, levels+1)

        # Now you can make the contour plot. 
        fig = plt.figure()
        
        #Add title label 
        fig.text(0.5, 0.95, label, ha='center', fontsize=10) 
        
        
        sw_units = [r"$cm^{-3}$", r"$km s^{-1}$", r"nT"]
        # label = file_label+"\n"+r"$n = {} $".format(self.density)+r"$cm^{-3}, $"+r"$v_x = {} $".format(self.vx)+r"$km s^{-1}, $"+r"$B_z = {} nT$".format(self.bz)
        # fig.text(0.5,0.9, label, ha="center")
        fig.subplots_adjust(wspace=0.4)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)

        # Constant y. 
        cont1 = ax1.contourf(xp_y, zp_y, letad_y, cmap=cmap, levels=levels, vmin=vmin, vmax=vmax)
        ax1.set_xlabel('X [RE]')
        ax1.set_ylabel('Z [RE]')
        ax1.set_title("XZ Plane")
        ax1.set_aspect("equal")

        # Add sketch of the Earth on top. 
        get_earth.make_earth(ax1, rotation=-90)

         # Constant z. 
        cont2 = ax2.contourf(xp_z, yp_z, letad_z, cmap=cmap, levels=levels, vmin=vmin, vmax=vmax)
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
            print (self.plot_path+"{}_emissivity_{}{}.png".format(file_label,self.masking_method,savetag))
            fig.savefig(self.plot_path+"{}_emissivity_{}{}.png".format(file_label,self.masking_method,savetag))

    def plot_both_planes_eta_components(self, cmap="hot", levels=100, vmin_veff=0, vmax_veff=2000, vmin_rho=0, vmax_rho=20, vmin_nH=0, vmax_nH=4, save=False, savetag=""):
        '''This will plot the effective velocity, proton density and neutral hydrogen
         density in the X-Z and X-Y planes side by side. These are the key quantities used 
         to calculate the emissivity.  
        '''

        #Get meridian data for veff. 
        xp_y, yp_y, zp_y, veffd_y, xp_z, yp_z, zp_z, veffd_z = gm.calculate_meridian_planes(self.x_3d, self.y_3d, self.z_3d, self.veff_3d)
        
        #Get meridian data for rho. 
        xp_y, yp_y, zp_y, rhod_y, xp_z, yp_z, zp_z, rhod_z = gm.calculate_meridian_planes(self.x_3d, self.y_3d, self.z_3d, self.rho_3d)
        
        #Get meridian data for nH. 
        xp_y, yp_y, zp_y, nHd_y, xp_z, yp_z, zp_z, nHd_z = gm.calculate_meridian_planes(self.x_3d, self.y_3d, self.z_3d, self.nH_3d)
        
        
        # Calculate log10 eta values for nH. If nH = 0, set log(nH) = vmin  
        lnHd_y = np.zeros(nHd_y.shape)+vmin_nH
        i = np.where(nHd_y != 0)
        lnHd_y[i] = np.log10(nHd_y[i])
        j = np.where(lnHd_y < vmin_nH)
        lnHd_y[j] = vmin_nH 


        #  # Calculate log10 eta values. If eta = 0, set log(eta) = vmin 
        lnHd_z = np.zeros(nHd_z.shape)+vmin_nH
        i = np.where(nHd_z != 0)
        lnHd_z[i] = np.log10(nHd_z[i])
        j = np.where(lnHd_z < vmin_nH)
        lnHd_z[j] = vmin_nH 

        # Create a filename label so you know which file you plotted. 
        file_label = self.openggcm.filename.split("/")[-1]
        file_label = file_label.split(".")[0]

        #Create a label so you can see which simulation this is. 
        label = "n = {} cm".format(self.openggcm.density)+r"$^{-3}$"+"   v = ({}, {}, {}) km/s   B = ({}, {}, {}) nT".format(self.openggcm.vx, self.openggcm.vy, self.openggcm.vz, self.openggcm.bx, self.openggcm.by, self.openggcm.bz)

        # Now you can make the contour plot. 
        fig = plt.figure(figsize=(5,8))
        sw_units = [r"$cm^{-3}$", r"$km s^{-1}$", r"nT"]
        # label = file_label+"\n"+r"$n = {} $".format(self.density)+r"$cm^{-3}, $"+r"$v_x = {} $".format(self.vx)+r"$km s^{-1}, $"+r"$B_z = {} nT$".format(self.bz)
        # fig.text(0.5,0.9, label, ha="center")
        fig.subplots_adjust(wspace=0.2, right=0.8, top=0.95)

        #Add title label 
        fig.text(0.5, 0.97, label, ha='center', fontsize=10) 


        # Effective velocity 
        ax1 = fig.add_subplot(321)
        ax2 = fig.add_subplot(322)

         # Get contour levels. 
        levels_veff = np.linspace(vmin_veff, vmax_veff, levels+1)
        
        # Constant y. 
        cont1 = ax1.contourf(xp_y, zp_y, veffd_y, cmap=cmap, levels=levels_veff, vmin=vmin_veff, vmax=vmax_veff)
        #ax1.set_xlabel('X [RE]')
        ax1.set_ylabel('Z [RE]')
        #ax1.set_title("veff: Y = {:.2f}".format(plane_value_y), fontsize=10)
        ax1.set_aspect("equal")

        # Add sketch of the Earth on top. 
        get_earth.make_earth(ax1, rotation=-90)

         # Constant z. 
        cont2 = ax2.contourf(xp_z, yp_z, veffd_z, cmap=cmap, levels=levels_veff, vmin=vmin_veff, vmax=vmax_veff)
        #ax2.set_xlabel('X [RE]')
        ax2.set_ylabel('Y [RE]')
        #ax2.set_title("veff: Z = {:.2f}".format(plane_value_z), fontsize=10)
        ax2.set_aspect("equal")

        # Add colourbar and earth. 
        
        cbar = plt.colorbar(cont2, ax=ax2, shrink=0.5)
        cbar.set_label("Effective Velocity\nkm s"+r"$^{-1}$")
        level_min = int(np.ceil(cont2.levels.min()))
        level_max = int(np.floor(cont2.levels.max()))
        cticks = np.linspace(level_min, level_max,5)
        print (cticks)
        cbar.set_ticks(cticks)
        # cbar.set_ticklabels([r'$10^{'+str(i)+'}$' for i in cticks])

        # Add sketch of the Earth on top. 
        get_earth.make_earth(ax2, rotation=-90)

        # Solar wind density 
        ax3 = fig.add_subplot(323)
        ax4 = fig.add_subplot(324)

         # Get contour levels. 
        levels_rho = np.linspace(vmin_rho, vmax_rho, levels+1)
        
        # Constant y. 
        cont3 = ax3.contourf(xp_y, zp_y, rhod_y, cmap=cmap, levels=levels_rho, vmin=vmin_rho, vmax=vmax_rho)
        #ax3.set_xlabel('X [RE]')
        ax3.set_ylabel('Z [RE]')
        #ax3.set_title("rho: Y = {:.2f}".format(plane_value_y), fontsize=10)
        ax3.set_aspect("equal")

        # Add sketch of the Earth on top. 
        get_earth.make_earth(ax3, rotation=-90)
    
        # Constant z. 
        cont4 = ax4.contourf(xp_z, yp_z, rhod_z, cmap=cmap, levels=levels_rho, vmin=vmin_rho, vmax=vmax_rho)
        #ax4.set_xlabel('X [RE]')
        ax4.set_ylabel('Y [RE]')
        #ax4.set_title("rho: Z = {:.2f}".format(plane_value_z), fontsize=10)
        ax4.set_aspect("equal")

        # Add colourbar and earth. 
        
        cbar = plt.colorbar(cont4, ax=ax4, shrink=0.5)
        cbar.set_label("Proton density\ncm"+r"$^{-3}$")
        level_min = int(np.ceil(cont4.levels.min()))
        level_max = int(np.floor(cont4.levels.max()))
        cticks = np.linspace(level_min, level_max,5)
        cbar.set_ticks(cticks)
        
        # SNeutral hydrogen density
        ax5 = fig.add_subplot(325)
        ax6 = fig.add_subplot(326)

         # Get contour levels. 
        levels_nH = np.linspace(vmin_nH, vmax_nH, levels+1)
        
        # Constant y. 
        cont5 = ax5.contourf(xp_y, zp_y, lnHd_y, cmap=cmap, levels=levels_nH, vmin=vmin_nH, vmax=vmax_nH)
        ax5.set_xlabel('X [RE]')
        ax5.set_ylabel('Z [RE]')
        #ax5.set_title("nH: Y = {:.2f}".format(plane_value_y), fontsize=10)
        ax5.set_aspect("equal")

        # Add sketch of the Earth on top. 
        get_earth.make_earth(ax5, rotation=-90)
        
        # Constant z. 
        cont6 = ax6.contourf(xp_z, yp_z, lnHd_z, cmap=cmap, levels=levels_nH, vmin=vmin_nH, vmax=vmax_nH)
        ax6.set_xlabel('X [RE]')
        ax6.set_ylabel('Y [RE]')
        #ax6.set_title("nH: Z = {:.2f}".format(plane_value_z), fontsize=10)
        ax6.set_aspect("equal")

        # Add colourbar and earth. 
        
        cbar = plt.colorbar(cont6, ax=ax6, shrink=0.5)
        cbar.set_label("H density\ncm"+r"$^{-3}$")
        level_min = int(np.ceil(cont6.levels.min()))
        level_max = int(np.floor(cont6.levels.max()))
        cticks = np.linspace(level_min, level_max,5)
        cbar.set_ticks(cticks)
        cbar.set_ticklabels([r'$10^{'+str(int(i))+'}$' for i in cticks])

        # Add sketch of the Earth on top. 
        get_earth.make_earth(ax6, rotation=-90)
        
        # Option to save figure. 
        if save: 
            fig.savefig(self.plot_path+"{}_em_key_comp{}.png".format(file_label,savetag))
    

