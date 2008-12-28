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
__version__ = '0.0.4'

import numpy, sim, math, fit

#  _____           ____
# | ____|___  _ __| __ )  ___  __ _ _ __ ___
# |  _| / _ \| '__|  _ \ / _ \/ _` | '_ ` _ \
# | |__| (_) | |  | |_) |  __/ (_| | | | | | |
# |_____\___/|_|  |____/ \___|\__,_|_| |_| |_|

class Beam(sim.Beam):
    def __init__(self, freqs, active_chans=None, **kwargs):
        CsAm = numpy.array([        # N. Gugliucci 08/07
            [ 2.3541326  ,-0.0039133512 , 1.6055088e-05,-2.7468911e-08], 
            [-0.46909345 , 0.0084471178 ,-5.1260711e-05, 1.0299793e-07],
            [ 0.32753617 ,-0.0081176326 , 5.8822952e-05,-1.3303273e-07],
            [-0.046844105, 0.00076223627,-3.5474502e-06, 4.4132035e-09],
            [-0.073523813, 0.0018151892 ,-1.3435102e-05, 3.1225928e-08],
            [ 0.047340855,-0.00085097424, 4.9799933e-06,-9.5164123e-09],
        ])
        #CsAm = numpy.array([       # Original K. Peek coeffs
        #   [ 4.3197202e+00,-3.6236564e-02, 1.9216586e-04,-3.4582946e-07,], 
        #   [-6.4876376e-01, 1.0860822e-02,-6.0990572e-05, 1.1385919e-07,],
        #   [-6.8774619e-01, 9.4422778e-03,-4.3093392e-05, 6.5261223e-08,],
        #   [-1.1915388e-03, 1.9346063e-05, 2.1812495e-07,-1.3628154e-09,],
        #   [ 1.9266323e-01,-2.8355537e-03, 1.3838885e-05,-2.2424197e-08,],
        #   [ 4.0644458e-02,-7.0476657e-04, 4.0026507e-06,-7.4626130e-09,],
        #])
        CsXo = numpy.array([        # N. Gugliucci 08/07
           [-121.29224, 1.9851554 ,-0.011876889  , 2.5222526e-05], 
           [ 76.969303,-1.3947796 , 0.0085644354 ,-1.7448153e-05],
           [-36.638691, 0.93699466,-0.0068616164 , 1.5544311e-05],
           [ 10.189859,-0.18212180, 0.00098309486,-1.6152395e-06],
           [ 5.9997050,-0.15737420, 0.0012090764 ,-2.8862905e-06],
           [-5.6561847, 0.10468756,-0.00063126068, 1.2444705e-06],
        ])
        #CsXo = numpy.array([       # Original K. Peek coeffs
        #   [-3.7377001e+02, 5.9933023e+00,-3.2753724e-02, 6.0860359e-05,], 
        #   [ 1.0590908e+02,-1.7623625e+00, 9.8737879e-03,-1.8423046e-05,],
        #   [ 7.2094475e+01,-9.2672015e-01, 3.8927716e-03,-5.3087386e-06,],
        #   [ 2.0078893e+00,-4.5110799e-02, 2.4801570e-04,-3.7067584e-07,],
        #   [-1.8188548e+01, 2.6325797e-01,-1.2543987e-03, 1.9708660e-06,],
        #   [-3.7902866e+00, 6.8599528e-02,-4.0342212e-04, 7.7325363e-07,],
        #])
        CsSd = numpy.array([        # N. Gugliucci 08/07
            [ 143.84525,-1.1088605  , 0.0048397670 ,-7.1054741e-06],
            [-104.00886, 1.9980993  ,-0.013304344  , 2.8955473e-05],
            [ 28.304230,-0.75088201 , 0.0056338561 ,-1.2898564e-05],
            [-8.7149717, 0.16108215 ,-0.00090283393, 1.5386691e-06],
            [-3.4672940, 0.091929321,-0.00071502397, 1.7311496e-06],
            [ 3.4123240,-0.063083812, 0.00038093617,-7.5356570e-07],
        ])
        #CsSd = numpy.array([       # Original K. Peek coeffs
        #   [ 1.8399317e+02,-1.6249565e+00, 6.9672259e-03,-9.6826907e-06,], 
        #   [-2.9244673e+01, 5.7229158e-01,-4.1730996e-03, 9.4038248e-06,],
        #   [-7.4733475e+01, 9.9345539e-01,-4.3074107e-03, 6.1619539e-06,],
        #   [ 3.3041833e+00,-3.7977133e-02, 1.7285586e-04,-3.4421417e-07,],
        #   [ 9.1105618e+00,-1.2763799e-01, 5.8072483e-04,-8.5291801e-07,],
        #   [ 2.6146875e+00,-4.6834018e-02, 2.7398547e-04,-5.2419876e-07,],
        #])
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
        sim.Beam.__init__(self, freqs, active_chans)
    def select_chans(self, active_chans):
        if active_chans is None: active_chans = numpy.arange(self.freqs.size)
        sim.Beam.select_chans(self, active_chans)
        self.BAm_sel = self.BAm.take(active_chans, axis=1)
        self.BXo_sel = self.BXo.take(active_chans, axis=1)
        self.BSd_sel = self.BSd.take(active_chans, axis=1)
    #def fresponse(self, zang, az):
    #    """Return the beam response across the band for input zenith angle
    #    (zang) and azimuth (az).  Using Katy Peek's original antenna
    #    modelling code."""
    #    import antresp
    #    func = lambda f: antresp.response(zang, az, f, \
    #        'Cs.amp.dat', 'Cs.xoff.dat', 'Cs.sigma.dat')
    #    return numpy.array(map(func, 1e3*self.active_freqs))
    def response(self, zang, az, pol=1):
        """Return the beam response across the band for input zenith angle
        (zang) and azimuth (az).  Rotate beam model 90 degrees if pol == 2."""
        if pol == 2: az += math.pi/2
        a = numpy.cos(numpy.array([0, 2*az, 4*az, 6*az, 8*az, 10*az]))
        a[0] = 0.5
        a1 = numpy.dot(a, self.BAm_sel)
        a2 = numpy.dot(a, self.BXo_sel)
        a3 = numpy.dot(a, self.BSd_sel)
        z = (180*zang/math.pi - a2) / a3
        return numpy.sqrt(a1 * numpy.exp(-z**2/2))

##   ____                _              _       
##  / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___ 
## | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __|
## | |__| (_) | | | \__ \ || (_| | | | | |_\__ \
##  \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/
#
## Define constants which are likely to change often
#ADC_CLK_RATE = .592         # GHz
#ACC_WINDOW = 8              # sync's per accumulation window
#CLK_PER_SYNC = 2**27        # clocks per sync
#NCHAN = 256
#NANTS = 8
#NPOL = 2
#
## Derived contants
#SFREQ = ADC_CLK_RATE / 8.
#BANDWIDTH = ADC_CLK_RATE / 4.
#SDF = BANDWIDTH / NCHAN
#INTTIME = CLK_PER_SYNC / (BANDWIDTH * 1e9) * ACC_WINDOW
#
#freqs = numpy.arange(NCHAN) * SDF + SFREQ
##active_chans = rfi.range2list((70,182), (194,235))
#active_chans = (155,)
#
## A default passband fit from PAPER's measured receiver gain
## (in Receiver_Gain.txt)
##GAIN_POLY = [-8.25e1, 8.34e1, -3.50e1, 7.79e1, -9.71e-1, 6.41e-2, -1.75e-2]
#GAIN_POLY = [.054]
#
#
##  ___       _ _   _       _ _          _   _             
## |_ _|_ __ (_) |_(_) __ _| (_)______ _| |_(_) ___  _ __  
##  | || '_ \| | __| |/ _` | | |_  / _` | __| |/ _ \| '_ \ 
##  | || | | | | |_| | (_| | | |/ / (_| | |_| | (_) | | | |
## |___|_| |_|_|\__|_|\__,_|_|_/___\__,_|\__|_|\___/|_| |_|
#
## Current location
##location = ('38:25:59.24', '-79:50:23:41', 806)     # Green Bank
#location = ('38:25:59.24', '-79:51:02.1', 806)     # Green Bank
#
## Beam to use
#beam = sim.Beam(freqs)
##beam = EorBeam(freqs)
#
## Antenna positions
#antennas = (
#    fit.FitAntenna(  -8.48, 455.28,   9.82 ), # 1
#    fit.FitAntenna( 205.47, 319.53,-251.71 ), # 2
#    fit.FitAntenna( 187.10,-352.95,-232.59 ), # 3
#    fit.FitAntenna(-262.70,-219.07, 318.70 ), # 4
#    fit.FitAntenna(-293.44,  -7.66, 360.20 ), # 5
#    fit.FitAntenna(-286.04,  93.20, 352.23 ), # 6
#    fit.FitAntenna(-182.66, 353.23, 227.56 ), # 7
#    #fit.FitAntenna( -84.27, 434.46, 107.19, beam, gain_poly=GAIN_POLY, 
#    fit.FitAntenna( -75.51, 433.83,  97.02, beam, gain_poly=GAIN_POLY, 
#        offset=1.4967), # 8
#)
#
## Sources for simulation
#src_dict = {
#    'Cygnus A':     fit.FitRadioFixedBody('19:57:44.5', '40:35.0',  freqs,
#        strength=[0, 13300.], spec_index=-1.),
#    #'0 Cygnus A':     fit.FitRadioFixedBody('19:57:44.5', '40:35.0',  freqs,
#    #    strength=10500, spec_index=-1.244),
#    'Cass A':       fit.FitRadioFixedBody('23:21:12.0', '58:32.1', freqs,
#        strength=[0, 11700.], spec_index=-1.),
#    #'2 Cass A':       fit.FitRadioFixedBody('23:21:12.0', '58:32.1', freqs,
#    #    strength=12800, spec_index=-0.770),
#    #'4 Crab':         fit.FitRadioFixedBody('05:31:31.5', '21:59.2', freqs,
#    #   strength=1500),
#    #'3 Virgo':        fit.FitRadioFixedBody('12:28:18.0', '12:40.1', freqs,
#    #   strength=1100),
#    #'1 Sun':          fit.FitRadioSun(freqs, strength=19996.)
#}
#
#fit_sim = fit.FitSimulator(antennas, location, src_dict, active_chans=active_chans)
#
## RFI
##bad_bins = rfi.range2list((0,69), (235,255), (183,194))
##rfi_freq_flagger = rfi.gen_freq_mfunc(bad_bins)


#  ____            _       _     ___       _             __
# / ___|  ___ _ __(_)_ __ | |_  |_ _|_ __ | |_ ___ _ __ / _| __ _  ___ ___
# \___ \ / __| '__| | '_ \| __|  | || '_ \| __/ _ \ '__| |_ / _` |/ __/ _ \
#  ___) | (__| |  | | |_) | |_   | || | | | ||  __/ |  |  _| (_| | (_|  __/
# |____/ \___|_|  |_| .__/ \__| |___|_| |_|\__\___|_|  |_|  \__,_|\___\___|
#                   |_|

#if __name__ == '__main__':
#    import sys, os, miriad
#    from optparse import OptionParser
#
#    p = OptionParser()
#    p.set_usage('params.py [options] *.uv')
#    p.set_description(__doc__)
#    opts, args = p.parse_args(sys.argv[1:])
#
#    mirparams = MiriadParams()
#
#    # Make a map function which changes uvw coordinates
#    def mfunc(p, d, v, c):
#        bl = int(p[-1])
#        u, v, w = mirparams.baselines[bl] 
#        return (u, v, w) + p[-2:], d
#        
#    for a in args:
#        print 'Working on file: ' + a
#        uvi = miriad.UV(a)
#        uvo = miriad.UV(a + '.new', 'new')
#        # Initialize uvo
#        for k in uvi.items: uvo.items[k] = uvi.items[k]
#        for v in mirparams:
#            val = mirparams[v]
#            uvo.vars.add_var(v, val[0])
#            try: uvo.vars[v] = val[1]
#            except: pass
#        hist_lines = 'PARAMS: Changed ant/baseline coords to version %s.' \
#                % (__version__)
#        # Fix uvw coordinates
#        miriad.map_uv(uvi, uvo, mfunc, append2history=hist_lines, 
#            ignorevars=('coord', 'baseline') + mirparams.statics)
#        del(uvi); del(uvo)
#        print '    Done.'
