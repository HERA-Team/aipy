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

def ParAng(ha,dec,lat):
    """
    For any hour angle, declenation in an image, calculate the paralactic angle at that point. Remember to multiply this by 2 when you're
    doing anything with it...
    """
    up = (np.cos(lat)*np.sin(ha))
    down = (np.sin(lat)*np.cos(dec))-(np.cos(lat)*np.sin(dec)*np.cos(ha))
    return np.arctan2(up,down)

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
    '''This is like fit.Antenna, except now it is expected that phsoff, amp,
    bp_r, and bp_i are dicts of {pol:value, ...}.'''
    def _update_phsoff(self):
        self.phsoff = {}
        for pol in self._phsoff: 
            self.phsoff[pol] = np.polyval(self._phsoff[pol], self.beam.afreqs)
    def _update_gain(self):
        self._gain = {}
        for pol in self.bp_r:
            bp = np.polyval(self.bp_r[pol], self.beam.afreqs) + \
                1j*np.polyval(self.bp_i[pol], self.beam.afreqs)
            self._gain[pol] = self.amp[pol] * bp
    def passband(self, conj=False, pol='x'):
        if conj: return np.conjugate(self._gain[pol])
        else: return self._gain[pol]
    def bm_response(self,top,pol='x'):
        """Introduce Stokes parameters in to the definition of the beam."""
        if pol in 'xy':
            return fit.Antenna.bm_response(self,top,pol)
        else:
            assert(pol in 'IQUV')
            if pol in 'IQ': return np.sqrt(0.5*fit.Antenna.bm_response(self,top,pol='x')**2+0.5*fit.Antenna.bm_response(self,top,pol='y')**2)
            if pol in 'UV': return np.sqrt(fit.Antenna.bm_response(self,top,pol='x')*fit.Antenna.bm_response(self,top,pol='y'))
    def get_params(self, prm_list=['*']):
        """Return all fitable parameters in a dictionary."""
        x,y,z = self.pos
        aprms = {'x':x, 'y':y, 'z':z}
        for p in self._phsoff:
            aprms.update({
                'dly_'+p:self._phsoff[p][-2],
                'off_'+p:self._phsoff[p][-1],
                'phsoff_'+p:self._phsoff[p]})
            aprms['bp_r_'+p] = list(self.bp_r[p])
            aprms['bp_i_'+p] = list(self.bp_i[p])
            aprms['amp_'+p] = self.amp[p]
        aprms.update(self.beam.get_params(prm_list))
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        changed = False
        self.beam.set_params(prms)
        try: self.pos[0], changed = prms['x'], True
        except(KeyError): pass
        try: self.pos[1], changed = prms['y'], True
        except(KeyError): pass
        try: self.pos[2], changed = prms['z'], True
        except(KeyError): pass
        for p in self._phsoff:
            try: self._phsoff[p][-2], changed = prms['dly_'+p], True
            except(KeyError): pass
            try: self._phsoff[p][-1], changed = prms['off_'+p], True
            except(KeyError): pass
            try: self._phsoff[p], changed = prms['phsoff_'+p], True
            except(KeyError): pass
            try: self.bp_r[p], changed = prms['bp_r_'+p], True
            except(KeyError): pass
            try: self.bp_i[p], changed = prms['bp_i'+p], True
            except(KeyError): pass
            try: self.amp[p], changed = prms['amp_'+p], True
            except(KeyError): pass
        if changed: self.update()
        return changed

#     _          n                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class AntennaArray(fit.AntennaArray):
    def set_active_pol(self, pol):
        assert(pol in ('xx','yy','xy','yx','I','Q','U','V')) # Now supports Stokes parameters as well as linear pols
        self.active_pol = pol
    def get_phs_offset(self, i, j):
        pol = self.get_active_pol()
        if pol in ['xx','yy','xy','yx']:
            return self[j].phsoff[pol[-1]] - self[i].phsoff[pol[0]]
        elif pol in ['I','Q','U','V']:
            if pol == 'I':
                return (self[j].phsoff['y'] - self[i].phsoff['y'])+(self[j].phsoff['x'] - self[i].phsoff['x'])
            if pol == 'Q':
                return (self[j].phsoff['y'] - self[i].phsoff['y'])-(self[j].phsoff['x'] - self[i].phsoff['x'])
            if pol == 'U':
                return (self[j].phsoff['y'] - self[i].phsoff['x'])+(self[j].phsoff['x'] - self[i].phsoff['y'])
            if pol == 'V':
                return np.conj((self[j].phsoff['y'] - self[i].phsoff['x'])-(self[j].phsoff['x'] - self[i].phsoff['y']))
    def passband(self, i, j):
        pol = self.get_active_pol()
        return self[j].passband(pol=pol[-1]) * self[i].passband(conj=True, pol=pol[0])
