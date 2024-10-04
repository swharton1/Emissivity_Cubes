# Emissivity_Cubes
All the code needed to read and produce emissivity cubes. 

File Descriptions: 
read_fits_cube.py - This will read in an emissivity cube stored in the FITS format. It will eventually be able to read in cubes created from any type of MHD simulation. 

read_openggcm_mhd.py - This reads in the CDF files produced by the OpenGGCM code and the json header files and add all the data to a python object in the same coordinate system as PPMLR. 

calc_emissivity.py - This takes in the object from read_openggcm_mhd.py and calculates the emissivity using standard assumptions. It also applies a mask to hide false emission from the magnetosphere. Several options are available. 

calc_flowlines.py - This takes in the object from read_openggcm_mhd.py and calculates a flowline for each point. It assigns a 1 or 0 to the origin array to determine whether that plasma originated from the solar wind. These values are used to mask false emission from the magnetosphere. 

cusp_id.py - Implementation of the cusp identification algorithm from Sun et al. (2019). Used because the flowlines method can miss the cusps. Used to unmask the cusps after the flowlines method by calc_emissivity.py 

get_meridians.py - Gets the XY and XZ meridian planes out of a 3D datacube for plotting. 


