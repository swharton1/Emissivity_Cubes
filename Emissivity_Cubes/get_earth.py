#This contains functions to add the Earth to plots. 

import numpy as np 
from matplotlib.patches import Wedge, Polygon, Circle

def make_earth(ax, rotation=0):
    '''This will add a little plot of the Earth on top for reference. This if for 2D   pictures. 
    
    Parameters
    ----------
    ax - axis object to add the earth too. 
    rotation - rotation of the Earth. def = 0. 
    '''

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
    

