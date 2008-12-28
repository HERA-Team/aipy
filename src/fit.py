"""
A module for fitting visibilities in a Miriad UV file.

Author: Aaron Parsons
Date: 01/14/2007
Revisions:
    03/13/2007  arp Modified to support other fitting programs in scipy,
                    including ability to fit w/o using a gradient.  Also
                    added functionality for only fitting some parameters
                    but having a standard fit file which still keeps track
                    of all info.  Only computes data relevant to selected
                    parameters.  Added docs.
"""
import ant, sim, numpy as n
from interp import interpolate

#  _   _ _   _ _ _ _           _____                 _   _                 
# | | | | |_(_) (_) |_ _   _  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
# | | | | __| | | | __| | | | | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# | |_| | |_| | | | |_| |_| | |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
#  \___/ \__|_|_|_|\__|\__, | |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#                      |___/        

def flatten_prms(prms, prm_list=[]):
    """Generate a list of parameters suitable for passing to a fitting
    algorithm from the heirarchical parameter dictionary 'prms', along
    with 'key_list' information for reconstructing such a dictionary from
    a list.  'prm_list' is only for recursion."""
    key_list = {}
    keys = prms.keys()
    keys.sort()
    for k in keys:
        if type(prms[k]) == dict:
            prm_list, new_key_list = flatten_prms(prms[k], prm_list)
            key_list[k] = new_key_list
        else:
            try:
                key_list[k] = (len(prm_list), len(prms[k]))
                prm_list += list(prms[k])
            except(TypeError):
                key_list[k] = (len(prm_list), 1)
                prm_list.append(prms[k])
    return prm_list, key_list

def reconstruct_prms(prm_list, key_list):
    """Generate a heirarchical parameter dictionary from a parameter
    list (prm_list) and 'key_list' information from flatten_prms."""
    prms = {}
    for k in key_list:
        v = key_list[k]
        if type(v) == dict: prms[k] = reconstruct_prms(prm_list, v)
        else:
            i, L = v
            if L > 1: prms[k] = prm_list[i:i+L]
            else: prms[k] = prm_list[i]
    return prms

def print_params(prms, indent='', grad=None):
    """Print a nice looking representation of a parameter dictionary."""
    keys = prms.keys()
    keys.sort()
    for k in keys:
        v = prms[k]
        if type(v) == dict:
            print indent, k
            if grad is None: print_params(v, indent + '  ')
            else: print_params(v, indent + '  ', grad[k])
        else:
            print indent, k,
            if grad is None:
                if not type(v) is list:
                    try: v = [list(v)]
                    except(TypeError): v = [v]
                if len(v) == 1: print v[0]
                else:
                    print
                    for i in v: print indent, ' ', i
            else:
                print
                print indent, v, '\t<', grad[k], '>'
                if not type(v) is list:
                    try: v = [list(v)]
                    except(TypeError): v = [v]
                for i in len(v):
                    print indent, ' ', v[i], '\t<', grad[k][i], '>'


#  ____           _ _       _____ _              _ ____            _       
# |  _ \ __ _  __| (_) ___ |  ___(_)_  _____  __| | __ )  ___   __| |_   _ 
# | |_) / _` |/ _` | |/ _ \| |_  | \ \/ / _ \/ _` |  _ \ / _ \ / _` | | | |
# |  _ < (_| | (_| | | (_) |  _| | |>  <  __/ (_| | |_) | (_) | (_| | |_| |
# |_| \_\__,_|\__,_|_|\___/|_|   |_/_/\_\___|\__,_|____/ \___/ \__,_|\__, |
#                                                                    |___/ 

class RadioFixedBody(sim.RadioFixedBody):
    """A class adding parameter fitting to sim.RadioFixedBody"""
    def get_params(self, prm_list=None):
        """Return all fitable parameters in a dictionary."""
        aprms = {
            'str':      list(self._janskies),
            'index':    list(self._index),
            'ra':       float(self._ra),
            'dec':      float(self._dec),
            'angsize':  float(self.angsize),
        }
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        if prms.has_key('str'):
            try:
                len(prms['str'])
                self._janskies = prms['str']
            except: self._janskies = [prms['str']]
        if prms.has_key('index'):
            try:
                len(prms['index'])
                self._index = prms['index']
            except: self._index = [prms['index']]
        try: self._ra = prms['ra']
        except(KeyError): pass
        try: self._dec = prms['dec']
        except(KeyError): pass
        try: self.angsize = prms['angsize']
        except(KeyError): pass

#  ____           _ _      ____                  _       _ 
# |  _ \ __ _  __| (_) ___/ ___| _ __   ___  ___(_) __ _| |
# | |_) / _` |/ _` | |/ _ \___ \| '_ \ / _ \/ __| |/ _` | |
# |  _ < (_| | (_| | | (_) |__) | |_) |  __/ (__| | (_| | |
# |_| \_\__,_|\__,_|_|\___/____/| .__/ \___|\___|_|\__,_|_|
#                               |_|                        

class RadioSpecial(sim.RadioSpecial):
    """A class adding parameter fitting to sim.RadioSpecial"""
    def get_params(self, prm_list=None):
        """Return all fitable parameters in a dictionary."""
        aprms = {
            'str':      list(self._janskies),
            'index':    list(self._index),
            'angsize':  float(self.angsize),
        }
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        if prms.has_key('str'):
            try:
                len(prms['str'])
                self._janskies = prms['str']
            except: self._janskies = [prms['str']]
        if prms.has_key('index'):
            try:
                len(prms['index'])
                self._index = prms['index']
            except: self._index = [prms['index']]
        try: self.angsize = prms['angsize']
        except(KeyError): pass

#  ____            ____      _        _             
# / ___| _ __ ___ / ___|__ _| |_ __ _| | ___   __ _ 
# \___ \| '__/ __| |   / _` | __/ _` | |/ _ \ / _` |
#  ___) | | | (__| |__| (_| | || (_| | | (_) | (_| |
# |____/|_|  \___|\____\__,_|\__\__,_|_|\___/ \__, |
#                                             |___/ 

class SrcCatalog(sim.SrcCatalog):
    """A class for fitting several celestial sources simultaneously."""
    def get_params(self, src_prms={'*':'*'}):
        """Return all fitable parameters in a dictionary."""
        prms = {}
        for k in src_prms:
            if k.startswith('*'): srcs = self.keys()
            else: srcs = [k]
            prm_list = src_prms[k]
            if type(prm_list) is str: prm_list = [prm_list]
            for s in srcs:
                try: prms[s] = self[s].get_params(prm_list)
                except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        for s in prms:
            try: self[s].set_params(prms[s])
            except(KeyError): pass

#  ____
# | __ )  ___  __ _ _ __ ___
# |  _ \ / _ \/ _` | '_ ` _ \
# | |_) |  __/ (_| | | | | | |
# |____/ \___|\__,_|_| |_| |_|

class BeamFlat(sim.BeamFlat):
    def get_params(self, prm_list=None):
        return {}
    def set_params(self, prms):
        pass

class Beam2DGaussian(sim.Beam2DGaussian):
    """A class adding parameter fitting to sim.Beam2DGaussian"""
    def get_params(self, prm_list=None):
        """Return all fitable parameters in a dictionary."""
        aprms = {'bm_xwidth':self.xwidth, 'bm_ywidth':self.ywidth}
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        xwidth, ywidth = None, None
        try: xwidth = prms['bm_xwidth']
        except(KeyError): pass
        try: ywidth = prms['bm_ywidth']
        except(KeyError): pass
        self.update(xwidth, ywidth)

class BeamPolynomial(sim.BeamPolynomial):
    """A class adding parameter fitting to sim.BeamPolynomial"""
    def get_params(self, prm_list=None):
        """Return all fitable parameters in a dictionary."""
        aprms = {'bm_poly':self.poly.flatten()}
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        try: self.update(prms['bm_poly'])
        except(KeyError): pass

class BeamCosSeries(sim.BeamCosSeries):
    """A class adding parameter fitting to sim.BeamCosSeries"""
    def get_params(self, prm_list=None):
        """Return all fitable parameters in a dictionary."""
        aprms = {'bm_poly_cos':self.poly_cos.flatten(),
                 'bm_poly_wid':self.poly_wid.flatten()}
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        try: self.update(poly_cos=prms['bm_poly_cos'])
        except(KeyError): pass
        try: self.update(poly_wid=prms['bm_poly_wid'])
        except(KeyError): pass

class BeamAlm(sim.BeamAlm):
    """A class adding parameter fitting to sim.BeamAlm"""
    def get_params(self, prm_list=[]):
        """Return all fitable parameters in a dictionary."""
        aprms = {}
        for i, a in enumerate(self.alm):
            k = 'alm%d' % i
            data = a.get_data()
            aprms[k] = n.array([data.real, data.imag]).transpose().flatten()
        prms = {}
        for p in prm_list:
            if p.startswith('*'): prms = aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        for p in prms:
            if not p.startswith('alm'): continue
            data = n.array(prms[p])
            p = int(p[3:])
            data.shape = (data.size/2, 2)
            data = data[:,0] + data[:,1] * 1j
            self.update(coeffs={p:data})

#     _          _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  / ___ \| | | | ||  __/ | | | | | | (_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

class Antenna(sim.Antenna):
    """A class adding parameter fitting to sim.Antenna"""
    def get_params(self, prm_list=None):
        """Return all fitable parameters in a dictionary."""
        x,y,z = self.pos
        aprms = {'x':x, 'y':y, 'z':z, 'delay':self.delay}
        aprms['bp_r'] = list(self.bp_r)
        aprms['bp_i'] = list(self.bp_i)
        aprms['amp'] = self.amp
        aprms.update(self.beam.get_params(prm_list))
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        self.beam.set_params(prms)
        try: self.pos[0] = prms['x']
        except(KeyError): pass
        try: self.pos[1] = prms['y']
        except(KeyError): pass
        try: self.pos[2] = prms['z']
        except(KeyError): pass
        try: self.delay = prms['delay']
        except(KeyError): pass
        self.update_gain(bp_r=prms.get('bp_r', None),
            bp_i=prms.get('bp_i', None), amp=prms.get('amp', None))

#     _          _                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class AntennaArray(sim.AntennaArray):
    """A class adding parameter fitting to sim.AntennaArray"""
    def get_params(self, ant_prms={'*':'*'}):
        """Return all fitable parameters in a dictionary."""
        prms = {}
        for k in ant_prms:
            if type(k) is str:
                if k.startswith('*'): ants = range(len(self.ants))
                else: continue
            else: ants = [k]
            prm_list = ant_prms[k]
            if type(prm_list) is str: prm_list = [prm_list]
            for a in ants: prms[a] = self.ants[a].get_params(prm_list)
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        for i, a in enumerate(self.ants):
            try: a.set_params(prms[i])
            except(KeyError): pass
        self.update_antennas(self.ants)

