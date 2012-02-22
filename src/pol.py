"""
Module for adding polarization information to models.
"""

from aipy import coord,fit,miriad
import numpy as np

#  _   ___     __
# | | | \ \   / /
# | | | |\ \ / /
# | |_| | \ V /
#  \___/   \_/
#

class UV(miriad.UV):
    def read_pol(self):
        """ Reliably read polarization metadata."""
        return miriad.pol2str[self._rdvr('pol','i')]
    def write_pol(self,pol):
        """Reliably write polarization metadata."""
        try: return self._wrvr('pol','i',miriad.str2pol[pol])
        except(KeyError): 
            print pol,"is not a reasonable polarization value!"
            return

#  _   _ _   _ _ _ _           _____                 _   _                 
# | | | | |_(_) (_) |_ _   _  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
# | | | | __| | | | __| | | | | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# | |_| | |_| | | | |_| |_| | |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
#  \___/ \__|_|_|_|\__|\__, | |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#                      |___/        

xy2s_m = np.array([[1.,   0.,  0.,  1.],
                   [1.,   0.,  0., -1.],
                   [0.,   1.,  1.,  0.],
                   [0., -1.j, 1.j,  0.]])

s2xy_m = np.linalg.inv(xy2s_m)

def ParAng(ha,dec,lat):
    """
    For any hour angle, declenation in an image, calculate the paralactic angle at that point. Remember to multiply this by 2 when you're
    doing anything with it...
    """
    tanX = (np.cos(lat)*np.sin(ha))/((np.sin(lat)*np.cos(dec))-(np.cos(lat)*np.sin(dec)*np.cos(ha)))
    up = (np.cos(lat)*np.sin(ha))
    down = (np.sin(lat)*np.cos(dec))-(np.cos(lat)*np.sin(dec)*np.cos(ha))
    return np.arctan2(up,down)

def stokes2xy(V_s):
    """Rotate a Stokes visibility to an XY visibility."""
    if type(V_s) == dict:
        try:
            V_s = np.array([V_s['I'],V_s['Q'],V_s['U'],V_s['V']])
            V_xy_arr = np.dot(s2xy_m,V_s)
            V_xy = {}
            for i,prm in enumerate(('xx','xy','yx','yy')):
                V_xy[prm] = V_xy_arr[i]
            return V_xy
        except(KeyError):
            print 'Label your data array differently!',V_s.keys()
            return None
    else: return np.dot(s2xy_m,V_xy)

def xy2stokes(V_xy):
    """Rotate an XY visibility into a Stokes' visibility."""
    if type(V_xy) == dict:
        try:
            V_xy = np.array([V_xy['xx'],V_xy['xy'],V_xy['yx'],V_xy['yy']])
            V_s_arr = np.dot(xy2s_m,V_xy)
            V_s = {}
            for i,prm in enumerate(('I','Q','U','V')):
                V_s[prm] = V_s_arr[i]
            return V_s
        except(KeyError):
            print 'Label your data array differently!',V_xy.keys()
            return None
    else: return np.dot(xy2s_m,V_xy)

def QU2p(V):
    """If you can't get an absolute polarization calibration, p = \sqrt{Q^2+U^2}/I may be useful. Do that transformation. Make sure input visibility is stored as a dictionary!!!"""
    V = normalizeI(V)
    try: V['p'] = np.sqrt(np.abs(V['Q'])**2 + np.abs(V['U'])**2)
    except(KeyError):
        V = xy2stokes(V)
        V['p'] = np.sqrt(np.abs(V['Q'])**2 + np.abs(V['U'])**2)
    return V

def normalizeI(V):
    """ Divide each visibility by Stokes' I."""
    try: I = V['I']
    except(KeyError):
        V_s = xy2stokes(V)
        I = V_s['I']
    for prm in V:
        V[prm] /= I
    return V

#  ____
# | __ )  ___  __ _ _ ___ ___
# |  _ \ / _ \/ _` | '_  `_  \
# | |_) |  __/ (_| | | | | | |
# |____/ \___|\__,_|_| |_| |_|

#     _          _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  / ___ \| | | | ||  __/ | | | | | | (_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

class Antenna(fit.Antenna):
    def bm_response(self,top,pol='x'):
        """Introduce Stoke' parameters in to the definition of the beam."""
        if pol in ('x','y'):
            return fit.Antenna.bm_response(self,top,pol)
        else:
            assert(pol in ('I','Q','U','V'))
            if pol in ('I','Q'): return np.sqrt(fit.Antenna.bm_response(self,top,pol='x')**2+fit.Antenna.bm_response(self,top,pol='y')**2)
            if pol in ('U','V'): return np.sqrt(2.*fit.Antenna.bm_response(self,top,pol='x')*fit.Antenna.bm_response(self,top,pol='y'))

#     _          _                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class AntennaArray(fit.AntennaArray):
    def get_phs_offset(self,i,j,pol):
        """This assumes you've run apply_cal.py before callihg this function."""
        if pol in ('xx','xy','yx','yy'): return fit.AntennaArray.get_phs_offset(self,i,j,pol)
        if pol in ('I','Q','U','V'): return np.zeros_like(self.get_afreqs()) 
    def passband(self,i,j,pol):
        """This assumes you've run apply_cal.py before calling this function."""
        if pol in ('xx','xy','yx','yy'): return fit.AntennaArray.passband(self,i,j,pol)
        if pol in ('I','Q','U','V'): return np.ones_like(self.get_afreqs()) 
    def bm_response(self,i,j,pol):
        """Introduce Stokes' parameters into the definition of the beam."""
        try: return fit.AntennaArray.bm_response(self,i,j,pol)
        except(AssertionError):
            assert(pol in ('I','Q','U','V'))
            if pol in ('I','Q'): return fit.AntennaArray.bm_response(self,i,j,'xx')+fit.AntennaArray.bm_response(self,i,j,'yy')
            if pol in ('U','V'): return 2.* fit.AntennaArray.bm_response(self,i,j,'xy')
