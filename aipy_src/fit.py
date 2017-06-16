"""
Module for reading and setting parameters in components of an AntennaArray
simulation for purpose of fitting.
"""
import amp, numpy as np

#  _   _ _   _ _ _ _           _____                 _   _                 
# | | | | |_(_) (_) |_ _   _  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
# | | | | __| | | | __| | | | | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# | |_| | |_| | | | |_| |_| | |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
#  \___/ \__|_|_|_|\__|\__, | |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#                      |___/        

def flatten_prms(prms, prm_list=None):
    """Generate list of parameters suitable for passing to fitting
    algorithm from heirarchical parameter dictionary 'prms', along
    with 'key_list' information for reconstructing such a dictionary from
    a list.  'prm_list' is only for recursion."""
    if prm_list is None: prm_list = []
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
    """Generate a heirarchical parameter dictionary from parameter
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
    """Print nice looking representation of a parameter dictionary."""
    keys = prms.keys()
    keys.sort()
    for k in keys:
        v = prms[k]
        if (type(v) is dict and v == {}) or v is None or \
                (type(v) is list and v == []):
            continue
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

class RadioFixedBody(amp.RadioFixedBody):
    """Class representing a source at fixed RA,DEC.  Adds get_params() and
    set_params() to amp.RadioFixedBody."""
    def get_params(self, prm_list=['*']):
        """Return all fitable parameters in a dictionary."""
        aprms = {
            'jys':      float(self._jys),
            'index':    float(self.index),
            'ra':       float(self._ra),
            'dec':      float(self._dec),
            'a1':       float(self.srcshape[0]),
            'a2':       float(self.srcshape[1]),
            'th':       float(self.srcshape[2]),
            'dra':      float(self.ionref[0]),
            'ddec':     float(self.ionref[1]),
            'mfreq':    float(self.mfreq),
        }
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        try: self._jys = prms['jys']
        except(KeyError): pass
        try: self.index = prms['index']
        except(KeyError): pass
        try: self._ra = prms['ra']
        except(KeyError): pass
        try: self._dec = prms['dec']
        except(KeyError): pass
        try: self.srcshape[0] = prms['a1']
        except(KeyError): pass
        try: self.srcshape[1] = prms['a2']
        except(KeyError): pass
        try: self.srcshape[2] = prms['th']
        except(KeyError): pass
        try: self.ionref[0] = prms['dra']
        except(KeyError): pass
        try: self.ionref[1] = prms['ddec']
        except(KeyError): pass

#  ____           _ _      ____                  _       _ 
# |  _ \ __ _  __| (_) ___/ ___| _ __   ___  ___(_) __ _| |
# | |_) / _` |/ _` | |/ _ \___ \| '_ \ / _ \/ __| |/ _` | |
# |  _ < (_| | (_| | | (_) |__) | |_) |  __/ (__| | (_| | |
# |_| \_\__,_|\__,_|_|\___/____/| .__/ \___|\___|_|\__,_|_|
#                               |_|                        

class RadioSpecial(amp.RadioSpecial):
    """Class representing moving sources (Sun,Moon,planets). Adds get_params() 
    and set_params() to amp.RadioSpecial."""
    def get_params(self, prm_list=['*']):
        """Return all fitable parameters in a dictionary."""
        aprms = {
            'jys':      float(self._jys),
            'index':    float(self.index),
            'a1':       float(self.srcshape[0]),
            'a2':       float(self.srcshape[1]),
            'th':       float(self.srcshape[2]),
            'dra':      float(self.ionref[0]),
            'ddec':     float(self.ionref[1]),
        }
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        try: self._jys = prms['jys']
        except(KeyError): pass
        try: self.index = prms['index']
        except(KeyError): pass
        try: self.srcshape[0] = prms['a1']
        except(KeyError): pass
        try: self.srcshape[1] = prms['a2']
        except(KeyError): pass
        try: self.srcshape[2] = prms['th']
        except(KeyError): pass
        try: self.ionref[0] = prms['dra']
        except(KeyError): pass
        try: self.ionref[1] = prms['ddec']
        except(KeyError): pass

#  ____            ____      _        _             
# / ___| _ __ ___ / ___|__ _| |_ __ _| | ___   __ _ 
# \___ \| '__/ __| |   / _` | __/ _` | |/ _ \ / _` |
#  ___) | | | (__| |__| (_| | || (_| | | (_) | (_| |
# |____/|_|  \___|\____\__,_|\__\__,_|_|\___/ \__, |
#                                             |___/ 

class SrcCatalog(amp.SrcCatalog):
    """Class for holding a catalog of celestial sources.  Adds get_params()
    and set_params() to amp.SrcCatalog."""
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

class Beam(amp.Beam):
    """Representation of a flat (gain=1) antenna beam pattern."""
    def get_params(self, prm_list=['*']):
        return {}
    def set_params(self, prms):
        return False

class Beam2DGaussian(amp.Beam2DGaussian):
    """Representation of a 2D Gaussian beam pattern, with default setting for 
    a flat beam."""
    def get_params(self, prm_list=['*']):
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
        changed = False
        try: self.xwidth, changed = prms['bm_xwidth'], True
        except(KeyError): pass
        try: self.ywidth, changed = prms['bm_ywidth'], True
        except(KeyError): pass
        if changed: self.update()
        return changed

class BeamPolynomial(amp.BeamPolynomial):
    """Representation of a gaussian beam model whose width varies with azimuth
    angle and with frequency."""
    def get_params(self, prm_list=['*']):
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
        changed = False
        try: 
            poly_azfreq = prms['bm_poly']
            poly_azfreq.shape = self.poly.shape
            self.poly = poly_azfreq
            changed = True
        except(KeyError): pass
        if changed: self.update()
        return changed

class BeamAlm(amp.BeamAlm):
    """Representation of a beam model where each pointing has a response
    defined as a polynomial in frequency, and the spatial distributions of 
    these coefficients decomposed into spherical harmonics."""
    def get_params(self, prm_list=['*']):
        """Return all fitable parameters in a dictionary."""
        aprms = {}
        for i, a in enumerate(self.alm):
            k = 'alm%d' % i
            data = a.get_data()
            aprms[k] = np.array([data.real, data.imag]).transpose().flatten()
        prms = {}
        for p in prm_list:
            if p.startswith('*'): prms = aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        changed = False
        for p in prms:
            if not p.startswith('alm'): continue
            changed = True
            data = np.array(prms[p])
            c = int(p[3:])
            data.shape = (data.size/2, 2)
            data = data[:,0] + data[:,1] * 1j
            if c < len(self.alm): self.alm[-1-c].set_data(data)
        if changed: self.update()
        return changed

#     _          _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  / ___ \| | | | ||  __/ | | | | | | (_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

class Antenna(amp.Antenna):
    """Representation of physical location and beam pattern of individual 
    antenna in array.  Adds get_params() and set_params() to amp.Antenna."""
    def get_params(self, prm_list=['*']):
        """Return all fitable parameters in a dictionary."""
        x,y,z = self.pos
        aprms = {'x':x, 'y':y, 'z':z, 'dly':self._phsoff[-2], 
            'off':self._phsoff[-1], 'phsoff':self._phsoff}
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
        changed = False
        self.beam.set_params(prms)
        try: self.pos[0], changed = prms['x'], True
        except(KeyError): pass
        try: self.pos[1], changed = prms['y'], True
        except(KeyError): pass
        try: self.pos[2], changed = prms['z'], True
        except(KeyError): pass
        try: self._phsoff[-2], changed = prms['dly'], True
        except(KeyError): pass
        try: self._phsoff[-1], changed = prms['off'], True
        except(KeyError): pass
        try: self._phsoff, changed = prms['phsoff'], True
        except(KeyError): pass
        try: self.bp_r, changed = prms['bp_r'], True
        except(KeyError): pass
        try: self.bp_i, changed = prms['bp_i'], True
        except(KeyError): pass
        try: self.amp, changed = prms['amp'], True
        except(KeyError): pass
        if changed: self.update()
        return changed

#     _          _                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class AntennaArray(amp.AntennaArray):
    """Representation of location and time of observation, and response of
    array of antennas as function of pointing and frequency.  Adds get_params()
    and set_params() to amp.AntennaArray."""
    def get_params(self, ant_prms={'*':'*'}):
        """Return all fitable parameters in a dictionary."""
        prms = {}
        for k in ant_prms:
            if k.startswith('*'): ants = map(str, range(len(self)))
            else: ants = [k]
            prm_list = ant_prms[k]
            if type(prm_list) is str: prm_list = [prm_list]
            for a in ants:
                try: prms[a] = self.ants[int(a)].get_params(prm_list)
                except(ValueError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        changed = False
        for i, a in enumerate(self):
            try: changed |= a.set_params(prms[str(i)])
            except(KeyError): pass
        if changed: self.update()
        return changed
