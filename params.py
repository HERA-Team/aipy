"""
A module which stores parameters for various observing runs by
the PAPER array.  Can be used as a script to update an existing UV file
with the current antenna positions (without touching the data).

Author: Aaron Parsons
Date: 11/07/2006
Revisions:
    11/10/2006  arp     Added veldop and vsource to miriad params.
                        ver = 0.0.2
    11/12/2006  arp     Used better fit antenna positions from antpos.06sepB
                        ver = 0.0.3
    11/13/2006  arp     Fixed ant,pos ordering in miriad header.
                        ver = 0.0.4
    12/05/2006  arp     Use RadioFixedBodys instead of RadioSources.  Added
                        RadioSun.
    01/10/2007  arp     Incorporated K.Peek's beam code into python.
    02/13/2007  arp     Took sqrt of beam response to make it voltage gain.
    03/05/2007  arp     Added select_chans to speed beam response computation.
"""

# Copyright (C) 2006 Aaron Parsons
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

__version__ = '0.0.4'

import numpy, rfi, antennas, math, fit, os

#def get_receiver_response(filename, nchan, sfreq, sdf):
#    """Retrieve the passband of a receiver from a file, and interpolate
#    to the number of frequency bins in the dataset's parameters."""
#    data = [L for L in open(filename).readlines() if not L[0] in ('!', '#')]
#    data = ' '.join(data)
#    data = numpy.array(map(float, data.split()))
#    data.shape = (len(data)/2, 2)
#    data[:,1] = 10**(data[:,1]/10)  # Convert from dB
#    data[:,0] /= 1e9                     # Convert to GHz
#    frequencies = numpy.arange(nchan) * sdf + sfreq
#    # Resample
#    pband = numpy.zeros((nchan,))
#    for i, f in enumerate(frequencies):
#        if f < data[0,0] or f >= data[-1,0]: pband[i] = 0
#        for j in range(data.shape[0]):
#            if data[j,0] > f:
#                pband[i] = data[j-1,1] + \
#                    (f-data[j-1,0])/(data[j,0]-data[j-1,0]) * \
#                    (data[j,1]-data[j-1,1])
#                break
#    return pband / 3e6

#  _____           ____
# | ____|___  _ __| __ )  ___  __ _ _ __ ___
# |  _| / _ \| '__|  _ \ / _ \/ _` | '_ ` _ \
# | |__| (_) | |  | |_) |  __/ (_| | | | | | |
# |_____\___/|_|  |____/ \___|\__,_|_| |_| |_|

class EorBeam(antennas.Beam):
    def __init__(self, freqs, active_chans=None):
        CsAm = numpy.array([
           [ 4.3197202e+00,-3.6236564e-02, 1.9216586e-04,-3.4582946e-07,], 
           [-6.4876376e-01, 1.0860822e-02,-6.0990572e-05, 1.1385919e-07,],
           [-6.8774619e-01, 9.4422778e-03,-4.3093392e-05, 6.5261223e-08,],
           [-1.1915388e-03, 1.9346063e-05, 2.1812495e-07,-1.3628154e-09,],
           [ 1.9266323e-01,-2.8355537e-03, 1.3838885e-05,-2.2424197e-08,],
           [ 4.0644458e-02,-7.0476657e-04, 4.0026507e-06,-7.4626130e-09,],
        ])
        CsXo = numpy.array([
           [-3.7377001e+02, 5.9933023e+00,-3.2753724e-02, 6.0860359e-05,], 
           [ 1.0590908e+02,-1.7623625e+00, 9.8737879e-03,-1.8423046e-05,],
           [ 7.2094475e+01,-9.2672015e-01, 3.8927716e-03,-5.3087386e-06,],
           [ 2.0078893e+00,-4.5110799e-02, 2.4801570e-04,-3.7067584e-07,],
           [-1.8188548e+01, 2.6325797e-01,-1.2543987e-03, 1.9708660e-06,],
           [-3.7902866e+00, 6.8599528e-02,-4.0342212e-04, 7.7325363e-07,],
        ])
        CsSd = numpy.array([
           [ 1.8399317e+02,-1.6249565e+00, 6.9672259e-03,-9.6826907e-06,], 
           [-2.9244673e+01, 5.7229158e-01,-4.1730996e-03, 9.4038248e-06,],
           [-7.4733475e+01, 9.9345539e-01,-4.3074107e-03, 6.1619539e-06,],
           [ 3.3041833e+00,-3.7977133e-02, 1.7285586e-04,-3.4421417e-07,],
           [ 9.1105618e+00,-1.2763799e-01, 5.8072483e-04,-8.5291801e-07,],
           [ 2.6146875e+00,-4.6834018e-02, 2.7398547e-04,-5.2419876e-07,],
        ])
        mhz_freqs = 1e3 * freqs # GHz -> MHz
        mhz_freqs = numpy.array([
           numpy.ones_like(mhz_freqs), 
           mhz_freqs, 
           mhz_freqs**2, 
           mhz_freqs**3,
        ])
        self.BAm = numpy.dot(CsAm, mhz_freqs)
        self.BXo = numpy.dot(CsXo, mhz_freqs)
        self.BSd = numpy.dot(CsSd, mhz_freqs)
        antennas.Beam.__init__(self, freqs, active_chans)
    def select_chans(self, active_chans):
        if active_chans is None: active_chans = numpy.arange(self.freqs.size)
        antennas.Beam.select_chans(self, active_chans)
        self.BAm_sel = self.BAm.take(active_chans, axis=1)
        self.BXo_sel = self.BXo.take(active_chans, axis=1)
        self.BSd_sel = self.BSd.take(active_chans, axis=1)
    def fresponse(self, zang, az):
        """Return the beam response across the band for input zenith angle
        (zang) and azimuth (az).  Using Katy Peek's original antenna
        modelling code."""
        import antresp
        func = lambda f: antresp.response(zang, az, f, \
            'Cs.amp.dat', 'Cs.xoff.dat', 'Cs.sigma.dat')
        return numpy.array(map(func, 1e3*self.active_freqs))
    def response(self, zang, az, pol=1):
        """Return the beam response across the band for input zenith angle
        (zang) and azimuth (az).  Rotate beam model 90 degrees if pol == 2."""
        #return 1
        if pol == 2: az += math.pi/2
        a = numpy.cos(numpy.array([0, 2*az, 4*az, 6*az, 8*az, 10*az]))
        a[0] = 0.5
        a1 = numpy.dot(a, self.BAm_sel)
        a2 = numpy.dot(a, self.BXo_sel)
        a3 = numpy.dot(a, self.BSd_sel)
        z = (180*zang/math.pi - a2) / a3
        return numpy.sqrt(a1 * numpy.exp(-z**2/2))

#  __  __ _      _           _ ____                               
# |  \/  (_)_ __(_) __ _  __| |  _ \ __ _ _ __ __ _ _ __ ___  ___ 
# | |\/| | | '__| |/ _` |/ _` | |_) / _` | '__/ _` | '_ ` _ \/ __|
# | |  | | | |  | | (_| | (_| |  __/ (_| | | | (_| | | | | | \__ \
# |_|  |_|_|_|  |_|\__,_|\__,_|_|   \__,_|_|  \__,_|_| |_| |_|___/

class MiriadParams(dict):
    """A class for passing information in this parameter file to Miriad.
    This is used in initializing new UV files."""
    def __init__(self):
        dict.__init__(self)
        self['source']  = ('a', 'zenith')
        self['operator']= ('a', 'C2M Python')
        self['version'] = ('a', __version__)
        self['telescop']= ('a', 'CASPER-8')
        # Fix a C-For indexing convention for ant,pos order
        a = numpy.array([i.pos for i in ants], dtype=numpy.float64)
        a = a.transpose()
        self['antpos']  = ('d', a.flatten())
        self['freq']    = ('d', SFREQ)
        self['inttime'] = ('r', INTTIME)
        self['nants']   = ('i', NANTS)
        self['nchan']   = ('i', NCHAN)
        self['nspect']  = ('i', 1)
        self['sfreq']   = ('d', SFREQ)
        self['sdf']     = ('d', BANDWIDTH / NCHAN)
        self['ischan']  = ('i', 1)
        self['nschan']  = ('i', NCHAN)
        self['restfreq']= ('d', SFREQ)
        self['npol']    = ('i', NPOL)
        self['epoch']   = ('r', 2000)
        self['veldop']  = ('r', 0.)
        self['vsource'] = ('r', 0.)
        self['longitu'] = ('d', ant_array.long)
        self['latitud'] = ('d', ant_array.lat)
        self['dec']     = ('d', ant_array.lat)
        self['obsdec']  = ('d', ant_array.lat)
        self['ra']      = ('d') 
        self['obsra']   = ('d') 
        self['lst']     = ('d') 
        self['pol']     = ('i')
        self.baselines = {}
        for i in range(NANTS):
            for j in range(i, NANTS):
                self.baselines[(i+1) << 8 | (j+1)] = tuple(ants[i] - ants[j])
        self.statics = (
            'antpos' , 'freq'   , 'inttime', 'nants'  , 'nchan'  , 'nspect' ,
            'sfreq'  , 'sdf'    , 'ischan' , 'nschan' , 'restfreq', 'npol'   ,
            'epoch'  , 'veldop' , 'vsource', 'longitu', 'latitud', 'dec'    ,)

#   ____                _              _       
#  / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___ 
# | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __|
# | |__| (_) | | | \__ \ || (_| | | | | |_\__ \
#  \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/

# Define constants which are likely to change often
ADC_CLK_RATE = .592         # GHz
ACC_WINDOW = 8              # sync's per accumulation window
CLK_PER_SYNC = 2**27        # clocks per sync
NCHAN = 256
NANTS = 8
NPOL = 2

# Derived contants
SFREQ = ADC_CLK_RATE / 8.
BANDWIDTH = ADC_CLK_RATE / 4.
SDF = BANDWIDTH / NCHAN
INTTIME = CLK_PER_SYNC / (BANDWIDTH * 1e9) * ACC_WINDOW

freqs = numpy.arange(NCHAN) * SDF + SFREQ
active_chans = rfi.range2list((70,182), (194,235))
# A default passband fit from PAPER's measured receiver gain
# (in Receiver_Gain.txt)
#GAIN_POLY = [-8.25e1, 8.34e1, -3.50e1, 7.79e1, -9.71e-1, 6.41e-2, -1.75e-2]
GAIN_POLY = [.054]


#  ___       _ _   _       _ _          _   _             
# |_ _|_ __ (_) |_(_) __ _| (_)______ _| |_(_) ___  _ __  
#  | || '_ \| | __| |/ _` | | |_  / _` | __| |/ _ \| '_ \ 
#  | || | | | | |_| | (_| | | |/ / (_| | |_| | (_) | | | |
# |___|_| |_|_|\__|_|\__,_|_|_/___\__,_|\__|_|\___/|_| |_|

# Current location
#location = ('38:25:59.24', '-79:50:23:41', 806)     # Green Bank
location = ('38:25:59.24', '-79:51:02.1', 806)     # Green Bank

# Beam to use
#beam = antennas.Beam(freqs)
beam = EorBeam(freqs)

# Default passband
#mypath = os.path.dirname(__file__)
#default_pband = get_receiver_response(mypath + '/data/Receiver_Gain.txt',
#    NCHAN, SFREQ, SDF)

# Antenna positions
ants = (
    fit.FitAntenna(  -8.48,     455.28,       9.82, beam, gain_poly=GAIN_POLY ), # 1
    #fit.FitAntenna(     -7.9050,   459.2426,    10.3881, beam, 
    #    delay=-0.9362, offset=-0.3088), # 1 
    fit.FitAntenna(   205.47,     319.53,    -251.71, beam, gain_poly=GAIN_POLY ), # 2
    #fit.FitAntenna(    203.0435,   322.4549,  -256.7504, beam, 
    #    delay=+4.8460, offset=-0.9116), # 2 
    fit.FitAntenna(   187.10,    -352.95,    -232.59, beam, 
        gain_poly=GAIN_POLY, offset=0.), # 3
    #fit.FitAntenna(    189.2293,  -351.9378,  -234.2475, beam, 
    #    delay=-13.6300, offset=+43.5261), # 3 
    fit.FitAntenna(  -262.70,    -219.07,     318.70, beam, gain_poly=GAIN_POLY ), # 4
    #fit.FitAntenna(   -262.2435,  -215.0549,   319.2504, beam, 
    #    delay=+0.0460, offset=-0.0116), # 4 
    fit.FitAntenna(  -293.44,      -7.66,     360.20, beam, gain_poly=GAIN_POLY ), # 5
    #fit.FitAntenna(   -292.9435,    -3.6549,   360.7504, beam, 
    #    delay=+0.0460, offset=-0.0116), # 5 
    fit.FitAntenna(  -286.04,      93.20,     352.23, beam, gain_poly=GAIN_POLY ), # 6
    #fit.FitAntenna(   -285.5435,    97.7549,   352.7504, beam, 
    #    delay=+0.0460, offset=-0.0116), # 6 
    fit.FitAntenna(-182.66,     353.23,     227.56, beam, gain_poly=GAIN_POLY ), # 7
    #fit.FitAntenna(   -184.5802,   359.1503,   229.6533, beam,
    #    delay=+0.9505, offset=-19.7864), # 7
    fit.FitAntenna( -84.27,     434.46,     107.19, beam, 
        gain_poly=GAIN_POLY, offset=1.53), # 8
    #fit.FitAntenna(    -84.0551,   439.3056,   108.2283, beam,
    #    delay=2.0485, offset= 1.6285), # 8
)

# Antenna array
antenna_array = fit.FitAntennaArray(ants, location, active_chans)

# Sources for simulation
sources = {
    'Cygnus A':     fit.FitRadioFixedBody('19:57:44.5', '40:35.0',  freqs,
        strength=[15000.], spec_index=-1.),
    #'0 Cygnus A':     fit.FitRadioFixedBody('19:57:44.5', '40:35.0',  freqs,
    #    strength=10500, spec_index=-1.244),
    'Cass A':       fit.FitRadioFixedBody('23:21:12.0', '58:32.1', freqs,
        strength=[10000.], spec_index=-1.),
    #'2 Cass A':       fit.FitRadioFixedBody('23:21:12.0', '58:32.1', freqs,
    #    strength=12800, spec_index=-0.770),
    #'4 Crab':         antennas.RadioFixedBody('05:31:31.5', '21:59.2', freqs,
    #   strength=1500),
    #'3 Virgo':        antennas.RadioFixedBody('12:28:18.0', '12:40.1', freqs,
    #   strength=1100),
    #'1 Sun':          antennas.RadioSun(freqs, strength=19996.)
}

source_list = fit.FitSourceList(sources, active_chans)

# RFI
bad_bins = rfi.range2list((0,69), (235,255), (183,194))
rfi_freq_flagger = rfi.gen_freq_mfunc(bad_bins)


#  ____            _       _     ___       _             __
# / ___|  ___ _ __(_)_ __ | |_  |_ _|_ __ | |_ ___ _ __ / _| __ _  ___ ___
# \___ \ / __| '__| | '_ \| __|  | || '_ \| __/ _ \ '__| |_ / _` |/ __/ _ \
#  ___) | (__| |  | | |_) | |_   | || | | | ||  __/ |  |  _| (_| | (_|  __/
# |____/ \___|_|  |_| .__/ \__| |___|_| |_|\__\___|_|  |_|  \__,_|\___\___|
#                   |_|

if __name__ == '__main__':
    import sys, os, miriad
    from optparse import OptionParser

    p = OptionParser()
    p.set_usage('params.py [options] *.uv')
    p.set_description(__doc__)
    opts, args = p.parse_args(sys.argv[1:])

    mirparams = MiriadParams()

    # Make a map function which changes uvw coordinates
    def mfunc(p, d, v, c):
        bl = int(p[-1])
        u, v, w = mirparams.baselines[bl] 
        return (u, v, w) + p[-2:], d
        
    for a in args:
        print 'Working on file: ' + a
        uvi = miriad.UV(a)
        uvo = miriad.UV(a + '.new', 'new')
        # Initialize uvo
        for k in uvi.items: uvo.items[k] = uvi.items[k]
        for v in mirparams:
            val = mirparams[v]
            uvo.vars.add_var(v, val[0])
            try: uvo.vars[v] = val[1]
            except: pass
        hist_lines = 'PARAMS: Changed ant/baseline coords to version %s.' \
                % (__version__)
        # Fix uvw coordinates
        miriad.map_uv(uvi, uvo, mfunc, append2history=hist_lines, 
            ignorevars=('coord', 'baseline') + mirparams.statics)
        del(uvi); del(uvo)
        print '    Done.'
