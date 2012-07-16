"""
Module for adding polarization information to models.
"""

from aipy import coord,fit,miriad
import numpy as n

#  _   ___     __
# | | | \ \   / /
# | | | |\ \ / /
# | |_| | \ V /
#  \___/   \_/
#

def ijp2blp(i,j,pol):
    return miriad.ij2bl(i,j) * 16 + (pol + 9)

def blp2ijp(blp):
    bl,pol = int(blp) / 16, (blp % 16) - 9
    i,j = miriad.bl2ij(bl)
    return i,j,pol

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

p2i = {'x':0,'y':1} #indices for each polarization

#Pauli Spin-Matrices
Sigma = {'t': n.matrix([[1,0],[0,1]]),
         'x': n.matrix([[0,1],[1,0]]),
         'y': n.matrix([[0,-1.j],[1.j,0]]),
         'z': n.matrix([[1,0],[0,-1]])}


def ParAng(ha,dec,lat):
    """
    For any hour angle, declenation in an image, calculate the paralactic angle at that point. Remember to multiply this by 2 when you're
    doing anything with it...
    """
    up = (n.cos(lat)*n.sin(ha))
    down = (n.sin(lat)*n.cos(dec))-(n.cos(lat)*n.sin(dec)*n.cos(ha))
    return n.arctan2(up,down)

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
    def __init__(self,x,y,z,beam,d=0., **kwargs):
        pol.Antenna.__init__(self,x,y,z,beam,**kwargs)
        self.d = d #I may want to update this to be a polynomial or something later (dfm)
    def G_i(self):
        """2x2 gain matrix"""
        amp_i = self.passband()
        phs_i = self.phsoff
        g_ix = amp_i[0]*n.exp(-2.j*n.pi*phs_i[0])
        g_iy = amp_i[1]*n.exp(-2.j*n.pi*phs_i[1])
        return [n.array([[g_ix[i],0.],[0.,g_iy[i]]]) for i in range(len(g_ix))]
    def D_i(self):
        """2x2 rotation matrix for this antenna -- to first order in rot_angle."""
        return n.array([[1.,self.d[1]],[-1.*n.conjugate(self.d[0]),1.]])
    def J_i(self):
        """Compute the Jones' matrix for this antenna."""
        return [n.dot(G_i[i],D_i) for i in range(len(G_i))]


#     _          n                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class AntennaArray(fit.AntennaArray):
    def gen_phs_nocal(self,src,i,j,pol,mfreq=0.150,ionref=None,srcshape=None,resolve_src=False):
        """Do the same thing as aa.gen_phs(), but don't apply delay/offset terms. This gets done in the Jones matrices."""
        if ionref is None:
            try: ionref = src.ionref
            except(AttributeError): pass
        if not ionref is None or resolve_src: u,v,w = self.gen_uvw(i,j,src=src)
        else: w = self.gen_uvw(i,j,src=src, w_only=True)
        if not ionref is None: w += self.refract(u, v, mfreq=mfreq, ionref=ionref)
        if resolve_src:
            if srcshape is None:
                try: res = self.resolve_src(u, v, srcshape=src.srcshape)
                except(AttributeError): res = 1
            else: res = self.resolve_src(u, v, srcshape=srcshape)
        else: res = 1
        phs = res * n.exp(-1j*2*n.pi*(w))
        return phs.squeeze()

