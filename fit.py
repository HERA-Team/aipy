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
import ant, sim, numpy

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
            print indent, k
            if grad is None:
                if not type(v) is list:
                    try: v = [list(v)]
                    except(TypeError): v = [v]
                for i in v: print indent, ' ', i
            else:
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
    """A class adding parameter fitting to RadioFixedBody"""
    def get_params(self, prm_list=None):
        """Return all fitable parameters in a dictionary."""
        aprms = {
            'str':      self._strength,
            'index':    self._spec_index,
            'ra':       float(self._ra),
            'dec':      float(self._dec),
            'f_c':      self.f_c,
            'size':     self._ang_size,
        }
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        try: self._strength = prms['str']
        except(KeyError): pass
        try: self._spec_index = prms['index']
        except(KeyError): pass
        try: self._ra = prms['ra']
        except(KeyError): pass
        try: self._dec = prms['dec']
        except(KeyError): pass
        try: self.f_c = prms['f_c']
        except(KeyError): pass
        try: self._ang_size = prms['size']
        except(KeyError): pass

#  ____           _ _      ____                  _       _ 
# |  _ \ __ _  __| (_) ___/ ___| _ __   ___  ___(_) __ _| |
# | |_) / _` |/ _` | |/ _ \___ \| '_ \ / _ \/ __| |/ _` | |
# |  _ < (_| | (_| | | (_) |__) | |_) |  __/ (__| | (_| | |
# |_| \_\__,_|\__,_|_|\___/____/| .__/ \___|\___|_|\__,_|_|
#                               |_|                        

class RadioSpecial(sim.RadioSpecial):
    """A class adding parameter fitting to RadioSpecial"""
    def get_params(self, prm_list=None):
        """Return all fitable parameters in a dictionary."""
        aprms = {
            'str':      self._strength,
            'index':    self._spec_index,
            'f_c':      self.f_c,
            'size':     self._ang_size,
        }
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        try: self._strength = prms['str']
        except(KeyError): pass
        try: self._spec_index = prms['index']
        except(KeyError): pass
        try: self.f_c = prms['f_c']
        except(KeyError): pass
        try: self._ang_size = prms['size']
        except(KeyError): pass

#  ____            ____      _        _             
# / ___| _ __ ___ / ___|__ _| |_ __ _| | ___   __ _ 
# \___ \| '__/ __| |   / _` | __/ _` | |/ _ \ / _` |
#  ___) | | | (__| |__| (_| | || (_| | | (_) | (_| |
# |____/|_|  \___|\____\__,_|\__\__,_|_|\___/ \__, |
#                                             |___/ 

class SrcCatalog(ant.SrcCatalog):
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

#     _          _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
#  / ___ \| | | | ||  __/ | | | | | | (_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

class Antenna(sim.Antenna):
    """A class adding parameter fitting to SimAntenna"""
    def get_params(self, prm_list=None):
        """Return all fitable parameters in a dictionary."""
        x,y,z = self.pos
        aprms = {'x':x, 'y':y, 'z':z, 'delay':self.delay, 'offset':self.offset}
        #aprms['gain_poly'] = list(self.gain_poly)
        aprms['spline'] = list(self.spline)
        aprms['amp'] = self.amp
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        try: self.pos[0] = prms['x']
        except(KeyError): pass
        try: self.pos[1] = prms['y']
        except(KeyError): pass
        try: self.pos[2] = prms['z']
        except(KeyError): pass
        try: self.delay = prms['delay']
        except(KeyError): pass
        try: self.offset = prms['offset']
        except(KeyError): pass
        #try: self.update_gain(gain_poly=prms['gain_poly'])
        try: self.update_gain(spline=prms['spline'])
        except(KeyError): pass
        try: self.update_gain(amp=prms['amp'])
        except(KeyError): pass

#     _          _                            _                         
#    / \   _ __ | |_ ___ _ __  _ __   __ _   / \   _ __ _ __ __ _ _   _ 
#   / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` | / _ \ | '__| '__/ _` | | | |
#  / ___ \| | | | ||  __/ | | | | | | (_| |/ ___ \| |  | | | (_| | |_| |
# /_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_/_/   \_\_|  |_|  \__,_|\__, |
#                                                                 |___/ 

class AntennaArray(sim.AntennaArray):
    """A class adding parameter fitting to AntennaArray"""
    def get_params(self, ant_prms={'*':'*'}):
        """Return all fitable parameters in a dictionary."""
        prms = {}
        for k in ant_prms:
            if type(k) is str:
                assert(k.startswith('*'))
                ants = range(len(self.ants))
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

