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

import numpy as n, sim, coord

#  ____
# | __ )  ___  __ _ _ __ ___
# |  _ \ / _ \/ _` | '_ ` _ \
# | |_) |  __/ (_| | | | | | |
# |____/ \___|\__,_|_| |_| |_|

class Beam(sim.Beam):
    """A specific beam model for the PAPER experiment.  This model is for
    a single dipole element."""
    def __init__(self, freqs, active_chans=None, **kwargs):
        """The axes of the Cs matrices are polynomials in frequency going
        to the right, and polys in cos(2*az) going down."""
        CsAm = n.array([        # N. Gugliucci 08/07
            [ 2.3541326  ,-0.0039133512 , 1.6055088e-05,-2.7468911e-08], 
            [-0.46909345 , 0.0084471178 ,-5.1260711e-05, 1.0299793e-07],
            [ 0.32753617 ,-0.0081176326 , 5.8822952e-05,-1.3303273e-07],
            [-0.046844105, 0.00076223627,-3.5474502e-06, 4.4132035e-09],
            [-0.073523813, 0.0018151892 ,-1.3435102e-05, 3.1225928e-08],
            [ 0.047340855,-0.00085097424, 4.9799933e-06,-9.5164123e-09],
        ])
        #CsAm = n.array([       # Original K. Peek coeffs
        #   [ 4.3197202e+00,-3.6236564e-02, 1.9216586e-04,-3.4582946e-07,], 
        #   [-6.4876376e-01, 1.0860822e-02,-6.0990572e-05, 1.1385919e-07,],
        #   [-6.8774619e-01, 9.4422778e-03,-4.3093392e-05, 6.5261223e-08,],
        #   [-1.1915388e-03, 1.9346063e-05, 2.1812495e-07,-1.3628154e-09,],
        #   [ 1.9266323e-01,-2.8355537e-03, 1.3838885e-05,-2.2424197e-08,],
        #   [ 4.0644458e-02,-7.0476657e-04, 4.0026507e-06,-7.4626130e-09,],
        #])
        CsXo = n.array([        # N. Gugliucci 08/07
           [-121.29224, 1.9851554 ,-0.011876889  , 2.5222526e-05], 
           [ 76.969303,-1.3947796 , 0.0085644354 ,-1.7448153e-05],
           [-36.638691, 0.93699466,-0.0068616164 , 1.5544311e-05],
           [ 10.189859,-0.18212180, 0.00098309486,-1.6152395e-06],
           [ 5.9997050,-0.15737420, 0.0012090764 ,-2.8862905e-06],
           [-5.6561847, 0.10468756,-0.00063126068, 1.2444705e-06],
        ])
        #CsXo = n.array([       # Original K. Peek coeffs
        #   [-3.7377001e+02, 5.9933023e+00,-3.2753724e-02, 6.0860359e-05,], 
        #   [ 1.0590908e+02,-1.7623625e+00, 9.8737879e-03,-1.8423046e-05,],
        #   [ 7.2094475e+01,-9.2672015e-01, 3.8927716e-03,-5.3087386e-06,],
        #   [ 2.0078893e+00,-4.5110799e-02, 2.4801570e-04,-3.7067584e-07,],
        #   [-1.8188548e+01, 2.6325797e-01,-1.2543987e-03, 1.9708660e-06,],
        #   [-3.7902866e+00, 6.8599528e-02,-4.0342212e-04, 7.7325363e-07,],
        #])
        CsSd = n.array([        # N. Gugliucci 08/07
            [ 143.84525,-1.1088605  , 0.0048397670 ,-7.1054741e-06],
            [-104.00886, 1.9980993  ,-0.013304344  , 2.8955473e-05],
            [ 28.304230,-0.75088201 , 0.0056338561 ,-1.2898564e-05],
            [-8.7149717, 0.16108215 ,-0.00090283393, 1.5386691e-06],
            [-3.4672940, 0.091929321,-0.00071502397, 1.7311496e-06],
            [ 3.4123240,-0.063083812, 0.00038093617,-7.5356570e-07],
        ])
        #CsSd = n.array([       # Original K. Peek coeffs
        #   [ 1.8399317e+02,-1.6249565e+00, 6.9672259e-03,-9.6826907e-06,], 
        #   [-2.9244673e+01, 5.7229158e-01,-4.1730996e-03, 9.4038248e-06,],
        #   [-7.4733475e+01, 9.9345539e-01,-4.3074107e-03, 6.1619539e-06,],
        #   [ 3.3041833e+00,-3.7977133e-02, 1.7285586e-04,-3.4421417e-07,],
        #   [ 9.1105618e+00,-1.2763799e-01, 5.8072483e-04,-8.5291801e-07,],
        #   [ 2.6146875e+00,-4.6834018e-02, 2.7398547e-04,-5.2419876e-07,],
        #])
        mhz_freqs = 1e3 * freqs # GHz -> MHz
        mhz_freqs = n.array([
           n.ones_like(mhz_freqs), 
           mhz_freqs, 
           mhz_freqs**2, 
           mhz_freqs**3,
        ])
        self.BAm = n.dot(CsAm, mhz_freqs)
        self.BXo = n.dot(CsXo, mhz_freqs)
        self.BSd = n.dot(CsSd, mhz_freqs)
        sim.Beam.__init__(self, freqs, active_chans)
    def select_chans(self, active_chans):
        if active_chans is None: active_chans = n.arange(self.freqs.size)
        sim.Beam.select_chans(self, active_chans)
        self.BAm_sel = self.BAm.take(active_chans, axis=1)
        self.BXo_sel = self.BXo.take(active_chans, axis=1)
        self.BSd_sel = self.BSd.take(active_chans, axis=1)
    def response(self, topo_xyz):
        """Return the beam response across the band for input zenith angle
        (zang) and azimuth (az).  Rotate beam model 90 degrees if pol == 'y'."""
        zang, phi = coord.xyz_to_th_phi(topo_xyz)
        az = -phi
        a = n.array([0,2,4,6,8,10],dtype=n.float)
        a.shape = (1,) + a.shape; az.shape += (1,); zang.shape += (1,)
        a = n.cos(n.dot(az, a))
        a[:,0] = 0.5
        a1 = n.dot(a, self.BAm_sel)
        a2 = n.dot(a, self.BXo_sel)
        a3 = n.dot(a, self.BSd_sel)
        z = (180*zang/n.pi - a2) / a3
        return n.sqrt(a1 * n.exp(-z**2/2))

