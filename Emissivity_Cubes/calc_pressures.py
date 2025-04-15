#These functions will calculate either dynamic or magnetic pressures. 

import numpy as np

def calc_dynamic_pressure(vx, vy, vz, density):
    '''Calculate dynamic pressure. Assumes solar wind is just protons.'''

    # Assume average ion mass = mass of proton. 
    mp = 0.00000000000000000000000000167
        
    # Calculate v in m/s 
    v = (vx**2 + vy**2 + vz**2)**0.5
    v = v*1000 

    # Convert number of particles /cm^3 to /m^3. 
    n = density*1000000

    # Calculate dynamic pressure first in Pascals, then nPa. 
    dyn_pressure = 0.5*mp*n*(v**2)
    dyn_pressure = dyn_pressure*1000000000

    return dyn_pressure

def calc_magnetic_pressure(bx, by, bz):
    '''Calculate the magnetic pressure'''

    # Calculate magnitude of B in T. 
    B = (bx**2 + by**2 + bz**2)**0.5
    B = B*0.000000001

    # mu0
    mu0 = 4*np.pi*0.0000001

    # Calculate magnetic pressure in Pa, then nPa. 
    mag_pressure = (B**2)/(2*mu0)
    mag_pressure = mag_pressure*1000000000

    return mag_pressure
