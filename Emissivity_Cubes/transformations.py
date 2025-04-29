#This contains a load of functions for doing a range of rotations and translations. 

import numpy as np 

#FUNCTIONS TO ROTATE VECTORS WITH MATRICES. 

def rotate_x(t, vector):
    '''This will rotate a vector around the x-axis by the angle t.
    
    Parameters
    ----------
    t - rotation angle in radians. 
    vector - (x,y,z) vector/matrix of dimension (3, n). Input vector.T if matrix.  
    
    Returns
    -------
    vx - rotated vector/matrix of vectors. ''' 
    
    #Make the rotation matrix. 
    rotx = np.array([[1,0,0],[0,np.cos(t),-np.sin(t)],[0,np.sin(t),np.cos(t)]])
    
    vx = np.matmul(rotx, vector.T)
    
    return vx 

def rotate_y(a, vector):
    '''This will rotate a vector around the y-axis by the angle a.
    
    Parameters
    ----------
    a - rotation angle in radians. 
    vector - (x,y,z) vector/matrix of dimension (3, n). Input vector.T if matrix. 
    
    Returns
    -------
    vy - rotated vector/matrix of vectors. ''' 
    
    #Make the rotation matrix. 
    roty = np.array([[np.cos(a),0,np.sin(a)],[0,1,0],[-np.sin(a),0,np.cos(a)]])
    
    vy = np.matmul(roty, vector.T)
    
    return vy

def rotate_z(phi, vector):
    '''This will rotate a vector around the z-axis by the angle phi.
    
    Parameters
    ----------
    phi - rotation angle in radians. 
    vector - (x,y,z) vector/matrix of dimension (3, n). Input vector.T if matrix. 
    
    Returns
    -------
    vz - rotated vector/matrix of vectors. ''' 
    
    #Make the rotation matrix. 
    rotz = np.array([[np.cos(phi),-np.sin(phi),0],[np.sin(phi),np.cos(phi),0],[0,0,1]])
    
    vz = np.matmul(rotz, vector.T)
    
    return vz

def rotate_3d(alpha, beta, gamma, vector):
    '''This will rotate a vector in 3D using t, a and phi. 
    Angles use the same names as wikipedia. 
    
    https://en.wikipedia.org/wiki/Rotation_matrix 
    
    Parameters
    ----------
    gamma - rotation angle around x axis in radians. (tilt)
    beta - rotation angle around y axis in radian. (-latitude)
    alpha - rotation angle around z axis in radians. (longitude)
    vector - (x,y,z) vector/matrix of dimension (3, n). Input vector.T if matrix. 
    
    Returns
    -------
    v - rotated vector/matrix of vectors. '''   
    
    top = np.array([np.cos(alpha)*np.cos(beta), np.cos(alpha)*np.sin(beta)*np.sin(gamma) - np.sin(alpha)*np.cos(gamma), np.cos(alpha)*np.sin(beta)*np.cos(gamma) + np.sin(alpha)*np.sin(gamma)]) 
    
    mid = np.array([np.sin(alpha)*np.cos(beta), np.sin(alpha)*np.sin(beta)*np.sin(gamma) + np.cos(alpha)*np.cos(gamma), np.sin(alpha)*np.sin(beta)*np.cos(gamma) - np.cos(alpha)*np.sin(gamma)])
    
    bottom = np.array([-np.sin(beta), np.cos(beta)*np.sin(gamma), np.cos(beta)*np.cos(gamma)]) 
    
    rot = np.zeros((3,3))
    rot[0] = top
    rot[1] = mid
    rot[2] = bottom
    
    v = np.matmul(rot, vector.T)
    
    return v
    
#FUNCTIONS TO CONVERT JUPITER COORDINATES TO IMAGE COORDINATES FOR OVERLAYS. 

def convert_cartesian_to_image_coords(fov, vector, sxi_loc):
    '''This will convert a vector in cartesian Jupiter coordinates to Cartesian image coordinates.
    
    Parameters
    ----------
    vector - (x,y,z) vector/matrix of dimensions (3,n) in Jupiter Cartesian coordinates. If just a single vector, input vector.T. 
    sxi_loc - SXI position vector for translation. 
    
    Returns 
    vi - (x,y,z) vector in Image Cartesian coordinates. 
    
    ''' 
    
    #Translate the vector first. 
    vector_trans = vector - sxi_loc
    
    #Now rotate the vector.  
    vix = np.dot(fov.xi, fov.x0)*vector_trans[0] + np.dot(fov.xi, fov.y0)*vector_trans[1] + np.dot(fov.xi, fov.z0)*vector_trans[2] 
    viy = np.dot(fov.yi, fov.x0)*vector_trans[0] + np.dot(fov.yi, fov.y0)*vector_trans[1] + np.dot(fov.yi, fov.z0)*vector_trans[2] 
    viz = np.dot(fov.zi, fov.x0)*vector_trans[0] + np.dot(fov.zi, fov.y0)*vector_trans[1] + np.dot(fov.zi, fov.z0)*vector_trans[2]    
        
    vi = np.array([vix, viy, viz]) 
    
    #print (vi)    
    return vi 

def convert_cartesian_to_image_coords_2(fov, vector, sxi_loc):
    '''This will convert a vector in cartesian Jupiter coordinates to Cartesian image coordinates using matrices. 
    
    Parameters
    ----------
    vector - (x,y,z) vector in Jupiter Cartesian coordinates. Input vector.T if matrix. 
    sxi_loc - SXI position vector for translation. 
    
    Returns 
    vi - (x,y,z) vector in Image Cartesian coordinates. 
    
    ''' 
    
    #Do translation first. 
    vector_trans = vector - sxi_loc
    
    #Create rotation matrix. 
    convert = np.array([[np.dot(fov.xi, fov.x0), np.dot(fov.xi, fov.y0), np.dot(fov.xi, fov.z0)],[np.dot(fov.yi, fov.x0), np.dot(fov.yi, fov.y0), np.dot(fov.yi, fov.z0)], [np.dot(fov.zi, fov.x0), np.dot(fov.zi, fov.y0), np.dot(fov.zi, fov.z0)]]) 
    
    #Now do the matrix multiplication with the vector. 
    
    vi = np.matmul(convert, vector_trans.T) 
    
    #print (vi)
    return vi
    
def get_angular_position(image_vector):
    '''This converts a vector in Cartesian image coordinates to angular positions within the image. 
        
    Parameters
    ----------
    image_vector - (x,y,z) in the image coordinate system. 
        
    Returns
    -------
    alpha - alpha position in image in radians. 
    beta - beta position in image in radians. 
        
    '''
        
    alpha = np.arctan2(-image_vector[1], image_vector[0])
    #beta = np.arctan2(image_vector[2],(image_vector[0]**2 + image_vector[1]**2)**0.5)
    beta = np.arctan2(image_vector[2], image_vector[0])   
    return alpha, beta 
