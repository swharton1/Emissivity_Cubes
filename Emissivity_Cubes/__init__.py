#This is the initialisation for Emissivity_Cubes. 

from . import read_openggcm_mhd
from . import calc_emissivity
from . import calc_flowlines
from . import cusp_id 
from . import get_meridians
from . import read_fits_cube
from . import read_ascii_cube
from . import calc_pressures
from . import get_earth 

from . import read_ppmlr
from . import ppmlr_fits
from . import ppmlr_image
from . import smile_fov
from . import smile_fov_limb 

from . import batsrus_fits 
from . import read_batsrus

from . import gse_gsm 
from . import convert_cube_to_gse
from . import transformations

import os 
#Set your paths to the data here and where you want any plots saving too. 
if "PPMLR_PATH" not in os.environ:
    os.environ["PPMLR_PATH"] = "/data/smile/PPMLR/"

if "OPENGGCM_PATH" not in os.environ:
    os.environ["OPENGGCM_PATH"] = "/data/smile/OpenGGCM/" 
    
if "BATSRUS_PATH" not in os.environ:
    os.environ["BATSRUS_PATH"] = "/data/smile/BATSRUS/"
    
if "PLOT_PATH" not in os.environ:
    #os.environ["PLOT_PATH"] = "/scratch/smile/sw682/"
    os.environ["PLOT_PATH"] = "/home/s/sw682/Code/plots/Emissivity_plots/" 
