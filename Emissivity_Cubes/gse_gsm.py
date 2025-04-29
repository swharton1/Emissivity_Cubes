#This contains the code to convert vectors between GSE and GSM. 

from spacepy import coordinates as coord
from spacepy.time import Ticktock 
import numpy as np 

def convert(x, y, z, time, calc_gsm=True):
    '''This function uses spacepy to convert vectors, or lists of vectors, between GSE and GSM. 
    
    Parameters
    ----------
    x - float or array of x positions. 
    y - float or array of y positions. 
    z - float or array of z positions. 
    time - datetime object. 
    calc_gsm - Boolean. If True, goes from GSE to GSM. Else goes from GSM to GSE. 
    
    Returns 
    -------
    x_end - float or array of x positions in other coordinate system. 
    y_end - float or array of x positions in other coordinate system. 
    z_end - float or array of x positions in other coordinate system. 
    
    '''
    
    #First, get direction. 
    start = 'GSE' if calc_gsm else 'GSM'
    end = 'GSM' if calc_gsm else 'GSE' 
           
    print (f'Convert coordinates from {start} to {end} with spacepy...')
        
    #Get size of xgse. 
    size = 1 if type(x) == float else len(x)
        
    coords_start = np.zeros((size,3))
    coords_start[:,0] = x
    coords_start[:,1] = y
    coords_start[:,2] = z
        
    #Needs to take data in a list of xyz points. 
    coord_obj = coord.Coords(coords_start, start, 'car')
    
    #There must be a time stamp for every datapoint. Turn into an array. 
    time_list = [time for i in range(size)]
       
    #Add time information. 
    coord_obj.ticks = Ticktock(time_list, 'UTC') 
       
    #To convert. 
    coords_end = coord_obj.convert(end, 'car') 
        
    x_end = coords_end.x[0]
    y_end = coords_end.y[0]
    z_end = coords_end.z[0] 
        
    return x_end, y_end, z_end 
