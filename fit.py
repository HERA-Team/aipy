"""
A module for fitting visibilities in a Miriad UV file.

Author: Aaron Parsons
Date: 01/14/2007
Revisions:
    None
"""

# Copyright (C) 2007 Aaron Parsons
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

import miriad, antennas, numpy, lbfgsb, os, params, pickle
from math import pi

#  _____ _ _   ____           _ _       ____            _       
# |  ___(_) |_|  _ \ __ _  __| (_) ___ | __ )  ___   __| |_   _ 
# | |_  | | __| |_) / _` |/ _` | |/ _ \|  _ \ / _ \ / _` | | | |
# |  _| | | |_|  _ < (_| | (_| | | (_) | |_) | (_) | (_| | |_| |
# |_|   |_|\__|_| \_\__,_|\__,_|_|\___/|____/ \___/ \__,_|\__, |
#                                                         |___/ 

class FitRadioBody(antennas.RadioBody):
    """A class adding parameter fitting to antennas.RadioBody"""
    def get_params(self, prm_list=None):
        """Return all fitable parameters in a dictionary."""
        p = {'strength':self._strength, 'spec_index':self._spec_index}
        if prm_list is None: return p
        else:
            np = {}
            try:
                for i in prm_list: np[i] = p[i]
            except(KeyError): pass
            return np
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        try: strength = prms['strength']
        except(KeyError): strength = self._strength
        try: spec_index = prms['spec_index']
        except(KeyError): spec_index = self._spec_index
        self.update(strength, spec_index)

class FitRadioFixedBody(antennas.RadioFixedBody, FitRadioBody):
    """A class adding parameter fitting to antennas.RadioFixedBody"""
    pass

#  _____ _ _      _          _                         
# |  ___(_) |_   / \   _ __ | |_ ___ _ __  _ __   __ _ 
# | |_  | | __| / _ \ | '_ \| __/ _ \ '_ \| '_ \ / _` |
# |  _| | | |_ / ___ \| | | | ||  __/ | | | | | | (_| |
# |_|   |_|\__/_/   \_\_| |_|\__\___|_| |_|_| |_|\__,_|

class FitAntenna(antennas.Antenna):
    """A class adding parameter fitting to antennas.Antenna"""
    def get_params(self, prm_list=None):
        """Return all fitable parameters in a dictionary."""
        p = {'pos':list(self.pos), 'delay':self.delay, 'offset':self.offset}
        p['gain_poly'] = list(self.gain_poly)
        if prm_list is None: return p
        else:
            np = {}
            try:
                for i in prm_list: np[i] = p[i]
            except(KeyError): pass
            return np
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        try: self.pos = numpy.array(prms['pos'], numpy.float64)
        except(KeyError): pass
        try: self.delay = prms['delay']
        except(KeyError): pass
        try: self.offset = prms['offset']
        except(KeyError): pass
        try: self.update_gain(prms['gain_poly'])
        except(KeyError): pass

# ,---.o|    ,---.     |                        ,---.                    
# |__. .|--- |---|,---.|--- ,---.,---.,---.,---.|---|,---.,---.,---.,   .
# |    ||    |   ||   ||    |---'|   ||   |,---||   ||    |    ,---||   |
# `    ``---'`   '`   '`---'`---'`   '`   '`---^`   '`    `    `---^`---|
#                                                                   `---'
class FitAntennaArray(antennas.SimAntennaArray):
    """A class adding parameter fitting to antennas.SimAntennaArray"""
    def get_params(self, ant_list=None, ant_prm_list=None):
        """Return all fitable parameters in a dictionary."""
        p = {}
        for i, a in enumerate(self.antennas):
            p[i] = a.get_params(ant_prm_list)
        if ant_list is None: return p
        else:
            np = {}
            try:
                for i in ant_list: np[i] = p[i]
            except(KeyError): pass
            return np
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        for i, a in enumerate(self.antennas):
            try: a.set_params(prms[i])
            except(KeyError): pass
        self.update_antennas(self.antennas)

class FitSourceList(antennas.SourceList):
    """A class for fitting several celestial sources simultaneously."""
    def get_params(self, src_list=None, src_prm_list=None):
        """Return all fitable parameters in a dictionary."""
        p = {}
        for n, s in zip(self.names, self.sources):
            p[n] = s.get_params(src_prm_list)
        if src_list is None: return p
        else:
            np = {}
            try:
                for i in src_list: np[i] = p[i]
            except(KeyError): pass
            return np
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        for n, s in zip(self.names, self.sources):
            try: s.set_params(prms[n])
            except(KeyError): pass

#  _____ _ _   ____                               
# |  ___(_) |_|  _ \ __ _ _ __ __ _ _ __ ___  ___ 
# | |_  | | __| |_) / _` | '__/ _` | '_ ` _ \/ __|
# |  _| | | |_|  __/ (_| | | | (_| | | | | | \__ \
# |_|   |_|\__|_|   \__,_|_|  \__,_|_| |_| |_|___/

class FitParams:
    """A structure for managing all of the fitting parameters associated
    with a simulation."""
    def __init__(self, antenna_array, source_list, uvfiles, 
            active_chans=None, fit_file=None):
        """antenna_array:    An AntennaArray stucture.
        source_list:      A SourceList structure.
        uvfiles:          A list of Miriad UV files.
        active_chans:     Channels to be included in simulation.  Default: all.
        fit_file:         File for retrieving and storing fit parameters."""
        if type(uvfiles) == str: uvfiles = [uvfiles]
        self.aa = antenna_array
        self.sl = source_list
        self.uvfiles = uvfiles
        # Initialize some parameters from the first of the UV files
        uv = miriad.UV(uvfiles[0])
        self.sdf = uv['sdf']
        self.sfreq = uv['sfreq']
        self.aa.update_location((uv['latitud'], uv['longitu']))
        del(uv)
        self.select_chans(active_chans)
        self.fit_file = fit_file
        self.best_fit = None
        #self.grad = None
        # Information about which parameters are active
        self.ant_list = None
        self.ant_prm_list = None
        self.src_list = None
        self.src_prm_list = None
        self.active_stokes = None
        self.fromfile()
    def select_chans(self, active_chans):
        """Choose only 'active_chans' for future freq calculations."""
        self.aa.select_chans(active_chans)
        self.sl.select_chans(active_chans)
        self.chans = self.sl.sources[0].chans
        self.freqs = self.chans * self.sdf + self.sfreq
    def tofile(self):
        """Save all current parameters to the fit_file specified in 
        initialization."""
        if self.fit_file is not None:
            prms = self.get_params()
            f = open(self.fit_file, 'w')
            pickle.dump(prms, f)
            f.close()
    def fromfile(self):
        """Retrieve parameters from the fit_file specified in initialization."""
        prms = self.get_params()        # Get current paremeters
        try:
            f = open(self.fit_file)
            new_prms = pickle.load(f)
            f.close()
        except(EOFError, IOError): new_prms = {}
        prms.update(new_prms)           # Update from any parameters listed
        self.set_params(prms)           # in the fit_file.
    def get_params(self, aslist=False, allparams=False):
        """Return a heirarchical dictionary of all parameters being fit."""
        if allparams:
            prms = self.aa.get_params()
            prms.update(self.sl.get_params())
        else:
            prms = self.aa.get_params(self.ant_list, self.ant_prm_list)
            prms.update(self.sl.get_params(self.src_list, self.src_prm_list))
        if aslist:
            prm_list, self.key_list = self.flatten_prms(prms)
            return prm_list
        else: return prms
    def flatten_prms(self, prms, prm_list=[]):
        """Generate a list of parameters suitable for passing to a fitting
        algorithm from the heirarchical parameter dictionary 'prms', along 
        with 'key_list' information for reconstructing such a dictionary from 
        a list.  'prm_list' is only for recursion."""
        key_list = {}
        keys = prms.keys()
        keys.sort()
        for k in keys:
            if type(prms[k]) == dict:
                prm_list, new_key_list = self.flatten_prms(prms[k], prm_list)
                key_list[k] = new_key_list
            else:
                try:
                    key_list[k] = (len(prm_list), len(prms[k]))
                    prm_list += list(prms[k])
                except(TypeError):
                    key_list[k] = (len(prm_list), 1)
                    prm_list.append(prms[k])
        return prm_list, key_list
    def reconstruct_prms(self, prm_list, key_list=None):
        """Generate a heirarchical parameter dictionary from a parameter
        list (prm_list) and 'key_list' information from flatten_prms."""
        if key_list is None: key_list = self.key_list
        prms = {}
        for k in key_list:
            v = key_list[k]
            if type(v) == dict: prms[k] = self.reconstruct_prms(prm_list, v)
            else:
                i, L = v
                if L > 1: prms[k] = prm_list[i:i+L]
                else: prms[k] = prm_list[i]
        return prms
    def get_vars(self):
        """Return a tuple of vars usable by lbfgsb for fitting constraints."""
        prm_list = self.get_params(aslist=True)
        return [lbfgsb.Var(p) for p in prm_list]
    def set_params(self, prms, fromlist=False):
        """Propagate parameter values from 'prms' to the simulators."""
        if fromlist: prms = self.reconstruct_prms(prms)
        self.aa.set_params(prms)
        self.sl.set_params(prms)
    def print_params(self, indent='', prms=None, grad=None):
        """Print a nice looking representation of a parameter dictionary."""
        if prms is None: prms = self.get_params()
        keys = prms.keys()
        keys.sort()
        for k in keys:
            v = prms[k]
            if type(v) == dict:
                print indent, k
                if grad is None: self.print_params(indent + '  ', v)
                else: self.print_params(indent + '  ', v, grad[k])
            else:
                print indent, k
                if grad is None:
                    if not type(v) is list: v = [v]
                    for i in v: print indent, ' ', i
                else:
                    print indent, v, '\t<', grad[k], '>'
                    if not type(v) is list: v = [v]
                    for i in len(v):
                        print indent, ' ', v[i], '\t<', grad[k][i], '>'
    def set_activity(self, ant_list=None, ant_prm_list=None, 
            src_list=None, src_prm_list=None, active_stokes=None):
        """Choose which antennas (ant_list), antenna parameters (ant_prm_list),
        sources (src_list), source parameters (src_prm_list), and stokes
        parameters are to be fit.  Valid antenna parameters are: 'pos',
        'gain_poly', 'offset', 'delay'.  Valid source parameters are:
        'strength', 'spec_index'."""
        self.ant_list = ant_list
        self.ant_prm_list = ant_prm_list
        self.src_list = src_list
        self.src_prm_list = src_prm_list
        self.active_stokes = active_stokes
    def is_active(self, bl, stokes):
        """Check if a given baseline, stokes parameter involves antennas and
        polarizations which are being fit for."""
        i, j = self.aa.gen_ij(bl)
        if not self.ant_list is None and \
                (i not in self.ant_list or j not in self.ant_list): 
            return False
        if not self.active_stokes is None and stokes not in self.active_stokes:
            return False
        return True
    def calc_fit(self, *prms):
        """Calculate the sum-of-squares fit between the simulation model and
        the data in the UV files using the provided parameters.  This function
        is meant to passed to a fitting algorithm which does not need gradient
        information."""
        prms = prms[0]
        self.set_params(prms, fromlist=True)
        self.print_params()
        fit = 0
        grad = numpy.zeros((len(prms),), dtype=numpy.double)
        plot_vec = {'data':[],'data_mask':[],'sim':[],'diff':[]}
        for uvf in self.uvfiles:
            uv = miriad.UV(uvf)
            uv.configure_preamble('time/baseline/pol')
            while True:
                preamble, data = uv.read_data()
                if data.size == 0: break
                t, b, p = preamble
                if not self.is_active(b, p): continue
                i, j = self.aa.gen_ij(b)
                # Skip auto-correlations
                if i == j: continue
                data = data.take(self.chans)
                plot_vec['data'].append(data.data[0])
                plot_vec['data_mask'].append(data.mask[0])
                # Create data
                self.aa.set_jultime(t)
                V_f, V_sf = self.aa.sim_data(self.sl.sources, b, stokes=p,
                        calc_grad=True)
                plot_vec['sim'].append(V_f[0])
                # Find fit (sum square)
                diff = data - V_f
                plot_vec['diff'].append(diff.data[0])
                # Ignore contributions from masked data
                diff = diff.filled(0)
                fit += sum(abs(diff)**2)
            del(uv)
        print '<', fit, '>', '...', 'Previous best:', self.best_fit
        print '----------------------------------------------------'
        if self.best_fit is None: self.best_fit = fit
        elif fit < self.best_fit:
            self.best_fit = fit
            self.tofile()
        # Optional plotting stuff
        if True:
            plot_vec['data'] = numpy.ma.array(plot_vec['data'], 
                mask=plot_vec['data_mask'])
            plot_vec['sim'] = numpy.array(plot_vec['sim'])
            plot_vec['diff'] = numpy.ma.array(plot_vec['diff'], 
                mask=plot_vec['data_mask'])
            import pylab
            pylab.clf()
            pylab.plot(plot_vec['data'].real, 'b')
            pylab.plot(plot_vec['sim'].real, 'r')
            pylab.plot(plot_vec['diff'].real, 'k')
            #pylab.plot(plot_vec['data'].imag, 'g')
            #pylab.plot(plot_vec['sim'].imag, 'k')
            pylab.show()
        return fit
    def calc_fit_grad(self, *prms):
        """Calculate the fit and the gradient of the fit given values of 
        prms."""
        prms = prms[0]
        if hash(str(self.get_params())) == hash(str(prms)): return
        self.set_params(prms, fromlist=True)
        self.print_params()
        fit = 0
        grad = numpy.zeros((len(prms),), dtype=numpy.double)
        plot_vec = {'data':[], 'data_mask':[], 'sim':[], 'diff':[]}
        for uvf in self.uvfiles:
            uv = miriad.UV(uvf)
            uv.configure_preamble('time/baseline/pol')
            while True:
                preamble, data = uv.read_data()
                if data.size == 0: break
                t, b, p = preamble
                if not self.is_active(b, p): continue
                i, j = self.aa.gen_ij(b)
                if i == j: continue
                data = data.take(self.chans)
                plot_vec['data'].append(data.data[0])
                plot_vec['data_mask'].append(data.mask[0])
                # Create data
                self.aa.set_jultime(t)
                V_f, V_sf = self.aa.sim_data(self.sl.sources, b, stokes=p,
                        calc_grad=True)
                plot_vec['sim'].append(V_f[0])
                # Find fit (sum square)
                diff = data - V_f
                plot_vec['diff'].append(diff.data[0])
                # Ignore contributions from masked data
                diff = diff.filled(0)
                fit += sum(abs(diff)**2)
                # Antenna positions
                #oi = self.param_len * (i-1)
                #oj = self.param_len * (j-1)
                #for s, src in enumerate(self.sl.sources):
                #    g_a = 2*pi*self.freqs*1j*V_sf[s]
                #    for n, sn in enumerate(self.aa.src_eq_vec(src)):
                #        g_an = g_a * sn
                #        g_an = -2*sum(diff.real*g_an.real+diff.imag*g_an.imag)
                #        # Position gradient for antenna j
                #        grad[oj+n] += g_an
                #        # Position gradient for antenna i
                #        grad[oi+n] -= g_an
                ## Antenna delay
                #g_tau = 2*pi*self.freqs*1j*V_f
                #g_tau = -2*sum(diff.real*g_tau.real+diff.imag*g_tau.imag)
                #grad[oj+3] += g_tau
                #grad[oi+3] -= g_tau
                ## Antenna phase offset
                #g_off = 2*pi*1j*V_f
                #g_off = -2*sum(diff.real*g_off.real+diff.imag*g_off.imag)
                #grad[oj+4] += g_off
                #grad[oi+4] -= g_off
                ## Antenna gain
                #nfreqs = len(self.freqs)
                #gj = self.aa.antennas[j-1].get_gain(self.chans)
                #Vgj_f = numpy.where(gj != 0, V_f / gj, 0)
                #g_gj_f = -2*(diff.real*Vgj_f.real+diff.imag*Vgj_f.imag)
                #grad[oj+5:oj+5+nfreqs] += numpy.where(gj > 0, g_gj_f, 0)
                #gi = self.aa.antennas[i-1].get_gain(self.chans)
                #Vgi_f = numpy.where(gi != 0, V_f / gi, 0)
                #g_gi_f = -2*(diff.real*Vgi_f.real+diff.imag*Vgi_f.imag)
                #grad[oi+5:oi+5+nfreqs] += numpy.where(gi > 0, g_gi_f, 0)
                ## Gradient for Sources
                #os = self.param_len * len(self.aa.antennas)
                #for s, src in enumerate(self.src_vals):
                #    #VI_f = VI_sf[s] / src._strength
                #    VI_f = V_sf[s] / src._strength
                #    grad[os+2*s] += -2*sum(diff.real*VI_f.real + \
                #        diff.imag*VI_f.imag)
                #    lnV_f = numpy.log(self.freqs / src._meas_freq) * V_sf[s]
                #    grad[os+2*s+1] += -2*sum(diff.real*lnV_f.real + \
                #        diff.imag*lnV_f.imag)
            del(uv)
        #print grad
        print '<', fit, '>'
        #if self.fit_file is not None: self.aa.tofile(self.fit_file)
        plot_vec['data'] = numpy.ma.array(plot_vec['data'], 
            mask=plot_vec['data_mask'])
        plot_vec['sim'] = numpy.array(plot_vec['sim'])
        plot_vec['diff'] = numpy.ma.array(plot_vec['diff'], 
            mask=plot_vec['data_mask'])
        import pylab
        pylab.clf()
        pylab.plot(plot_vec['data'].real, 'b')
        pylab.plot(plot_vec['sim'].real, 'r')
        pylab.plot(plot_vec['diff'].real, 'k')
        #pylab.plot(plot_vec['data'].imag, 'g')
        #pylab.plot(plot_vec['sim'].imag, 'k')
        #pylab.plot(diff.real, 'g')
        pylab.show()
        print '--------------------'
        #self.fit = fit
        #self.grad = grad
        #return fit, list(grad)
        return fit
    #def calc_fit(self, *prms):
    #    self.calc_fit_grad(prms)
    #    return self.fit
    #def calc_grad(self, *prms):
    #    self.calc_fit_grad(prms)
    #    return self.grad

#  _            _   _                     _
# | |_ ___  ___| |_| |__   ___ _ __   ___| |__
# | __/ _ \/ __| __| '_ \ / _ \ '_ \ / __| '_ \
# | ||  __/\__ \ |_| |_) |  __/ | | | (__| | | |
#  \__\___||___/\__|_.__/ \___|_| |_|\___|_| |_|

if __name__ == '__main__':
    import params, sys
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('fit.py [options] *.uv')
    p.set_description(__doc__)
    p.add_option('-f', '--fit_file', dest='fit_file',
        help='Output file where latest fit is stored.')
    opts, args = p.parse_args(sys.argv[1:])

    fp = FitParams(params.antenna_array, params.source_list, args,
        active_chans=(155,), fit_file=opts.fit_file)
        #active_chans=params.active_chans, fit_file=opts.fit_file)

    # Manually override some parameters
    #prms = fp.get_params(aslist=True)
    #prms['Cass A']['strength'] = [1000., 1000., 1120., 13060.]
    #prms['Cygnus A']['strength'] = [-1000., -1000., -1060., 15870.]
    #fp.set_params(prms)

    # Choose which parameters to fit
    fp.set_activity(ant_list=[2, 7], 
        ant_prm_list=['pos', 'offset'],
        src_list=['Cygnus A', 'Cass A'],
        src_prm_list=['strength'],
        active_stokes=[-5])

    # Do the fit
    #prms, fit = lbfgsb.minimize(fp.calc_fit_grad, fp.get_vars(),
    #    pgtol=1e-2, factr=1e8)
    import scipy.optimize
    #print scipy.optimize.fmin_cg(fp.calc_fit, fp.get_params(), 
    #    fprime=fp.calc_grad, gtol=1e6, norm=100)
    scipy.optimize.fmin(fp.calc_fit, fp.get_params(aslist=True))
    print 'Done!'
