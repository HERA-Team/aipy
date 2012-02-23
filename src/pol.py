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
    def __init__(self, x, y, z, beam, pol='x', num=-1, bp_r=n.array([1]),
             bp_i=n.array([0]), amp=1, pointing=(0.,n.pi/2,0), **kwargs):
        fit.Antenna.__init__(self, x, y, z, beam, **kwargs)
        self.pol = pol
        self.num = num
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
    def get_ant_list(self):
        # XXX this could be done in init
        # could also save a lot of code repetition later if we change how this
        # and getitem below work...
        """Define a consistent numbering system for dual-pol antenna array. 
        Return a dictionary of antenna names a their corresponding indices."""
        try: 
            ants = {}
            for i,ant in enumerate(self): ants[str(ant.num)+str(ant.pol)] = i
            return ants
        except(NameError): return [str(i) for i in self.ants]
    def __getitem__(self, item): 
        # this should follow *args syntax of phs.py
        if type(item) is str:
            return self.ants.__getitem__(self.get_ant_list()[item])
        else:
            return self.ants.__getitem__(item)
    def get_phs_offset(self,i,j,*args):
        if len(args)>0:
            pol = args[0]
        else: 
            return fit.AntennaArray.get_phs_offset(self,i,j)
        """This assumes you've run apply_cal.py before callihg this function."""
        if pol in ('xx','xy','yx','yy'):
            ants = self.get_ant_list()
            try: #if we have pol info, use it
                return self[str(j)+pol[1]].phsoff - self[str(i)+pol[0]].phsoff
            except(KeyError):
                return self[j].phsoff - self[i].phsoff
            except(UnboundLocalError):
                return self[j].phsoff - self[i].phsoff
        elif pol in ('I','Q','U','V'): return np.zeros_like(self.get_afreqs()) 
    def gen_phs(self, src, i, j, pol, mfreq=.150, ionref=None, srcshape=None, 
             resolve_src=False):
        """Return phasing that is multiplied to data to point to src."""
        if ionref is None:
            try: ionref = src.ionref
            except(AttributeError): pass
        if not ionref is None or resolve_src: u,v,w = self.gen_uvw(i,j,src=src)
        else: w = self.gen_uvw(i,j,src=src, w_only=True)
        if not ionref is None: w += self.refract(u, v, mfreq=mfreq, ionref=ionref)
        o = self.get_phs_offset(i,j,pol)
        phs = n.exp(-1j*2*n.pi*(w + o))
        if resolve_src:
            if srcshape is None:
                try: res = self.resolve_src(u, v, srcshape=src.srcshape)
                except(AttributeError): res = 1
            else: res = self.resolve_src(u, v, srcshape=srcshape)
        else: res = 1
        o = self.get_phs_offset(i,j,pol)
        phs = res * n.exp(-1j*2*n.pi*(w + dw + o))
        return phs.squeeze()
    def phs2src(self, data, src, i, j, pol='xx', mfreq=.150, ionref=None, srcshape=None):
        """Apply phasing to zenith-phased data to point to src."""
        return data * self.gen_phs(src, i, j, pol,
            mfreq=mfreq, ionref=ionref, srcshape=srcshape, resolve_src=False)
    def unphs2src(self,data,src, i, j, pol='xx', mfreq=.150, ionref=None, srcshape=None):
        """Remove phasing from src-phased data to point to zenith."""
        return data / self.gen_phs(src, i, j, pol,
            mfreq=mfreq, ionref=ionref, srcshape=srcshape, resolve_src=False)
    def passband(self,i,j,*args):
        if len(args)>0:
            pol = args[0]
            ants = self.get_ant_list()
            if pol in ('xx','xy','yx','yy'):
                return self[ants[str(i)+pol[0]]].passband() * self[ants[str(j)+pol[1]]].passband(conj=True)
            #This assumes you've run apply_cal.py before calling this function for IQUV.
            elif pol in ('I','Q','U','V'): return np.ones_like(self.get_afreqs())
        else: return a.fit.AntennaArray.passband(self, i, j)
    def bm_response(self,i,j,*args):
        if len(args)>0:
            pol = args[0]
        else: 
            return fit.AntennaArray.bm_response(self,i,j)
        """Introduce Stokes' parameters into the definition of the beam."""
        try: return fit.AntennaArray.bm_response(self,i,j,pol)
        except(AssertionError):
            assert(pol in ('I','Q','U','V'))
            if pol in ('I','Q'): return fit.AntennaArray.bm_response(self,i,j,'xx')+fit.AntennaArray.bm_response(self,i,j,'yy')
            if pol in ('U','V'): return 2.* fit.AntennaArray.bm_response(self,i,j,'xy')
    def sim(self, i, j, pol='xx',resolve_src=True):
        """Simulate visibilites for the specified (i,j) baseline and 
        polarization.  sim_cache() must be called at each time step before 
        this will return valid results."""
        assert(pol in ('xx','yy','xy','yx'))
        if self._cache is None:
            raise RuntimeError('sim_cache() must be called before the first sim() call at each time step.')
        elif self._cache == {}:
            return n.zeros_like(self.passband(i,j,pol))
        s_eqs = self._cache['s_eqs']
        u,v,w = self.gen_uvw(i, j, src=s_eqs)
        I_sf = self._cache['jys']
        Gij_sf = self.passband(i,j,pol)
        Bij_sf = self.bm_response(i,j,pol=pol)
        if len(Bij_sf.shape) == 2: Gij_sf = n.reshape(Gij_sf, (1, Gij_sf.size))
        # Get the phase of each src vs. freq, also does resolution effects
        E_sf = n.conjugate(self.gen_phs(s_eqs, i, j, pol, mfreq=self._cache['mfreq'],
            srcshape=self._cache['s_shp'], ionref=self._cache['i_ref'],
            resolve_src=resolve_src))
        try: E_sf.shape = I_sf.shape
        except(AttributeError): pass
        # Combine and sum over sources
        GBIE_sf = Gij_sf * Bij_sf * I_sf * E_sf
        Vij_f = GBIE_sf.sum(axis=0)
        return Vij_f
    def get_params(self, ant_prms={'*':'*'}):
        """Return all fitable parameters in a dictionary."""
        prms = {}
        for k in ant_prms:
            ants = self.get_ant_list()
            if k.startswith('*'): ants = self.get_ant_list()
            else: ants = {k:ants[k]}
            prm_list = ant_prms[k]
            if type(prm_list) is str: prm_list = [prm_list]
            for ant,i in ants.iteritems():
                try: prms[ant] = self.ants[i].get_params(prm_list)
                except(ValueError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        ants = self.get_ant_list()
        for ant,i in ants.iteritems():
            try: self.ants[i].set_params(prms[ant])
            except(KeyError): pass
        self.update()

