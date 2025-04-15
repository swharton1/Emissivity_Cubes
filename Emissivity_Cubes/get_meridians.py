#This will work out the XZ and XY planes from a 3D dataset. 
import numpy as np 

def get_log_emissivity(x, y, z, eta, vmin, vmax):
    '''This will get the log of the emissivity in XZ and XY planes.

    Parameters
    ----------
    x - 3D grid.
    y - 3D grid. 
    z - 3D grid. 
    eta - 3D grid. 
    vmin - minimum logged emissivity to show. 
    vmax - maximum logged emissivity to show. 

    Returns
    -------
    xp_y - 2D grid of x values for XZ plane. 
    zp_y - 2D grid of z values for XZ plane. 
    leta_y - 2D grid of logged eta values for XZ plane. 
    xp_z - 2D grid of x values for XY plane. 
    yp_z - 2D grid of y values for XY plane. 
    leta_z - 2D grid of logged eta values for XY plane.
 
    '''

    #Get meridian data for etad. 
    xp_y, yp_y, zp_y, eta_y, xp_z, yp_z, zp_z, eta_z = calculate_meridian_planes(x, y, z, eta)

    # Calculate log10 eta values. If eta = 0, set log(eta) = vmin  
    #XY plane. 
    leta_y = np.zeros(eta_y.shape)+vmin
    i = np.where(eta_y != 0)
    leta_y[i] = np.log10(eta_y[i])
    j = np.where(leta_y < vmin)
    leta_y[j] = vmin 

    #XZ plane. 
    leta_z = np.zeros(eta_z.shape)+vmin
    i = np.where(eta_z != 0)
    leta_z[i] = np.log10(eta_z[i])
    j = np.where(leta_z < vmin)
    leta_z[j] = vmin 

    return xp_y, zp_y, leta_y, xp_z, yp_z, leta_z


def calculate_meridian_planes(x_3d, y_3d, z_3d, var_3d):
    '''This will actually work out the XZ and XY plane data properly by taking means between the nearest planes
        
    Parameters
    ----------
    x_3d - 3D array of x positions
    y_3d - 3D array of y positions
    z_3d - 3D array of z positions 
    var_3d - 3D array of data (e.g emissivity)
        
    Returns
    -------
    xp_y - 2D array of x values in XZ plane. 
    yp_y - 2D array of y values in XZ plane. 
    zp_y - 2D array of z values in XZ plane. 
    var_y - 2D array of data values in XZ plane. 
    xp_z - 2D array of x values in XY plane. 
    yp_z - 2D array of y values in XY plane. 
    zp_z - 2D array of z values in XY plane. 
    var_z - 2D array of data values in XY plane. 
        
    ''' 
            
    #XZ plane. 
    #Get lowest positive y value and highest negative y value. 
    i_yl = np.where(y_3d[0,:,0] < 0, y_3d[0,:,0], -np.inf).argmax()
    i_yu = np.where(y_3d[0,:,0] > 0, y_3d[0,:,0], np.inf).argmin()
    
    xp_yl = x_3d[:,i_yl]
    xp_yu = x_3d[:,i_yu]
    xp_y = (xp_yl+xp_yu)/2.
        
    yp_yl = y_3d[:,i_yl]
    yp_yu = y_3d[:,i_yu]
    yp_y = (yp_yl+yp_yu)/2.
        
    zp_yl = z_3d[:,i_yl]
    zp_yu = z_3d[:,i_yu]
    zp_y = (zp_yl+zp_yu)/2.
        
    var_yl = var_3d[:,i_yl]
    var_yu = var_3d[:,i_yu]
    var_y = (var_yl+var_yu)/2.
        
    #XY plane. 
    #Get lowest positive z value and highest negative z value. 
    i_zl = np.where(z_3d[:,0,0] < 0, z_3d[:,0,0], -np.inf).argmax()
    i_zu = np.where(z_3d[:,0,0] > 0, z_3d[:,0,0], np.inf).argmin()
        
    xp_zl = x_3d[i_zl]
    xp_zu = x_3d[i_zu]
    xp_z = (xp_zl+xp_zu)/2. 
        
    yp_zl = y_3d[i_zl]
    yp_zu = y_3d[i_zu]
    yp_z = (yp_zl+yp_zu)/2. 
        
    zp_zl = z_3d[i_zl]
    zp_zu = z_3d[i_zu]
    zp_z = (zp_zl+zp_zu)/2. 
        
    var_zl = var_3d[i_zl]
    var_zu = var_3d[i_zu]
    var_z = (var_zl+var_zu)/2. 
        
    return xp_y, yp_y, zp_y, var_y, xp_z, yp_z, zp_z, var_z 
    
    
def calculate_sunearth_line(x_3d, y_3d, z_3d, var_3d):
    '''This will correctly calculate the Earth-Sun line data along the x axis. 
    
    Parameters
    ----------
    x_3d - 3D array of x positions
    y_3d - 3D array of y positions
    z_3d - 3D array of z positions 
    var_3d - 3D array of data (e.g emissivity)
    
    Returns
    -------
    xp_mean - 1D array of x positions along x axis. 
    yp_mean - 1D array of y positions along x axis. 
    zp mean - 1D array of z positions along x axis. 
    varp_mean - 1D array of var values along x axis. 
    
    '''
    
    #Get lowest positive y value and highest negative y value. 
    i_yl = np.where(y_3d[0,:,0] < 0, y_3d[0,:,0], -np.inf).argmax()
    i_yu = np.where(y_3d[0,:,0] > 0, y_3d[0,:,0], np.inf).argmin()
    
    #Get lowest positive z value and highest negative z value. 
    i_zl = np.where(z_3d[:,0,0] < 0, z_3d[:,0,0], -np.inf).argmax()
    i_zu = np.where(z_3d[:,0,0] > 0, z_3d[:,0,0], np.inf).argmin()
    
    #Get mean x values. 
    xp_1 = x_3d[i_zl,i_yl]
    xp_2 = x_3d[i_zl,i_yu]
    xp_3 = x_3d[i_zu,i_yl]
    xp_4 = x_3d[i_zu,i_yu]
    xp_mean = (xp_1+xp_2+xp_3+xp_4)/4.
    
    #Get mean y values. 
    yp_1 = y_3d[i_zl,i_yl]
    yp_2 = y_3d[i_zl,i_yu]
    yp_3 = y_3d[i_zu,i_yl]
    yp_4 = y_3d[i_zu,i_yu]
    yp_mean = (yp_1+yp_2+yp_3+yp_4)/4.
    
    #Get mean z values. 
    zp_1 = z_3d[i_zl,i_yl]
    zp_2 = z_3d[i_zl,i_yu]
    zp_3 = z_3d[i_zu,i_yl]
    zp_4 = z_3d[i_zu,i_yu]
    zp_mean = (zp_1+zp_2+zp_3+zp_4)/4.
    
    #Get mean z values. 
    varp_1 = var_3d[i_zl,i_yl]
    varp_2 = var_3d[i_zl,i_yu]
    varp_3 = var_3d[i_zu,i_yl]
    varp_4 = var_3d[i_zu,i_yu]
    varp_mean = (varp_1+varp_2+varp_3+varp_4)/4.
    
    return xp_mean, yp_mean, zp_mean, varp_mean 
