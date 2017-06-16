"""
Module for fitting 2d Gaussians to things.
"""

import numpy as np

def moments(data):
    total = np.abs(data).sum()
    Y,X = np.indices(data.shape)
    y = np.argmax((X*np.abs(data)).sum(axis=1)/total)
    x = np.argmax((Y*np.abs(data)).sum(axis=0)/total)
    col = data[int(y),:]
    row = data[:,int(x)]
    #First moment.
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)*col).sum()/np.abs(col).sum())
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)*row).sum()/np.abs(row).sum())
    width = ( width_x + width_y ) / 2.
    height = np.median(data.ravel())
    amplitude = data.max() - height
    mylist = [amplitude,x,y]
    if np.isnan(width_y) or np.isnan(width_x) or np.isnan(height) or np.isnan(amplitude):
        raise ValueError("Somehthing is nan")
    mylist = [height] + mylist
    mylist = mylist + [width_x,width_y]
    return mylist

def twodgaussian(inpars,shape=None):
    inpars_old = inpars
    inpars = list(inpars)
    height = inpars.pop(0)
    height = float(height)
    amplitude,center_y,center_x = inpars.pop(0),inpars.pop(0),inpars.pop(0)
    amplitude = float(amplitude)
    center_x = float(center_x)
    center_y = float(center_y)
    width_x,width_y = inpars.pop(0),inpars.pop(0)
    rcen_x = center_x
    rcen_y = center_y
    if len(inpars) > 0:
        raise ValueError("There are still input parameters:" + str(inpars) + \
                        " and you've input: " + str(inpars_old) + \
                        " circle=%d, rotate=%d, vheight=%d" % (circle,rotate,vheight) )
    def rotgauss(x,y):
        xp = x
        yp = y
        g = amplitude*np.exp( 
            -(((rcen_x-xp)/width_x)**2 + 
            ((rcen_y-yp)/width_y)**2)/2.)
        return g
    if shape is not None:
        return rotgauss(*np.indices(shape))
    else:
        return rotgauss
