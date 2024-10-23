# Emissivity_Cubes
All the code needed to read and produce emissivity cubes. 

Setup: 
The only setup you should need to do is open Emissivity_Cubes/__init__.py and edit the filepaths to where you keep the data, and where you want any plots to go. 


File Descriptions: 

read_fits_cube.py - This will read in an emissivity cube stored in the FITS format. It will eventually be able to read in cubes created from any type of MHD simulation. Originally based on emissivity cubes derived from OpenGGCM.  

read_openggcm_mhd.py - This reads in the CDF files produced by the OpenGGCM code and the json header files and add all the data to a python object in the same coordinate system as PPMLR. 

calc_emissivity.py - This takes in the object from read_openggcm_mhd.py and calculates the emissivity using standard assumptions. It also applies a mask to hide false emission from the magnetosphere. Several options are available. 

calc_flowlines.py - This takes in the object from read_openggcm_mhd.py and calculates a flowline for each point. It assigns a 1 or 0 to the origin array to determine whether that plasma originated from the solar wind. These values are used to mask false emission from the magnetosphere. 

cusp_id.py - Implementation of the cusp identification algorithm from Sun et al. (2019). Used because the flowlines method can miss the cusps. Used to unmask the cusps after the flowlines method by calc_emissivity.py 

get_meridians.py - Gets the XY and XZ meridian planes out of a 3D datacube for plotting. 

read_ppmlr.py - Code to read in the ASCII PPMLR files. 

ppmlr_fits.py - Code to read in the FITS PPMLR files. Will be made to do the same job as read_fits_cube.py so they are consistent. This file is effectively redundant. 

ppmlr_image.py - Code to produce an image through a PPMLR cube. Will be made to do the same for all cubes. 

smile_fov.py - Object to create the pointing directions and LOS coordinates for SXI. This version is unconstrained and can be pointed anywhere. 

smile_fov_limb.py Object to create the pointing directions and LOS coordinates for SXI, including the constraints to keep a constant angle with the limb of the earth, always pointing at the GSE x axis and having the x axis of the image pointing towards the Earth. UPDATED VERSION HAS BEEN ADDED. 

smile_fov_limb_old.py Previous constrained SMILE object that had an error in the orientation of SXI. DON'T USE THIS ONE. It's here for reference for those that have used the previous code.  

read_batsrus.py - Code to read in emissivity cubes derived from BATSRUS simulations. Based on emissivit cubes created by Andrey Samsonov in the same format as the emissivity cubes derived from PPMLR by Tianran Sun. 

batsrus_fits.py - Code to convert the BATSRUS ASCII files to FITS files in the common format. FITS files can be read in with read_fits_cube.py. 
