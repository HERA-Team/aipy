'''This is a calibration file for data collected at PAPER in Green Bank
on JD 2454574.'''

import aipy as a, numpy as n

class XTalkAntennaArray(a.fit.AntennaArray):
    '''Add xtalk fitting capability to AntennaArray.'''
    def __init__(self, *args, **kwargs):
        a.fit.AntennaArray.__init__(self, *args, **kwargs)
        self.xtalk = {}
        for bl in self.bl_order:
            self.xtalk['xa_%d' % bl] = []
            self.xtalk['xp_%d' % bl] = []
    def get_params(self, ant_prms={'*':'*','aa':'*'}):
        prms = a.fit.AntennaArray.get_params(self, ant_prms=ant_prms)
        if not ant_prms.has_key('aa'): return prms
        for p in ant_prms['aa']:
            if p.startswith('*'): prms.update(self.xtalk)
            else:
                try: 
                    if len(self.xtalk[p]) > 0: prms[p] = self.xtalk[p]
                except(KeyError): pass
        return prms
    def set_params(self, prms):
        a.fit.AntennaArray.set_params(self, prms)
        for key in self.xtalk:
            try: self.xtalk[key] = prms[key]
            except(KeyError): pass
    def sim_cache(self, *args, **kwargs):
        a.fit.AntennaArray.sim_cache(self, *args, **kwargs)
        for bl in self.bl_order:
            apoly = self.xtalk['xa_%d' % bl]
            ppoly = self.xtalk['xp_%d' % bl]
            x = n.polyval(apoly, self.ants[0].beam.afreqs).astype(n.complex)
            x *= n.exp(2*n.pi*1j*n.polyval(ppoly, self.ants[0].beam.afreqs))
            self._cache['x_%d' % bl] = x
    def sim(self, i, j, pol='xx'):
        Vij_f = a.fit.AntennaArray.sim(self, i, j, pol=pol)
        bl = 'x_%d' % a.miriad.ij2bl(i,j)
        return Vij_f + self._cache[bl]

prms = {
    'loc': ('38:25:59.24',  '-79:51:02.1'), # Green Bank, WV
    'antpos':
        [[   0.00,    0.00,    0.00],
         [ 214.11, -137.07, -259.71],
         [ 196.52, -810.36, -253.62],
         [-255.16, -674.34,  308.17],
         [-285.31, -463.01,  350.48],
         [-277.44, -361.92,  342.78],
         [-174.60, -102.56,  218.51],
         [ -76.55,  -20.95,   98.04]],
    'delays': [0.000,-10.722,-0.995, 8.309, 9.467, 2.992, 3.311, 0.000],
    'offsets': [0., 2.335, 3.017, 0., 0., 0., 0., 0.],
    'amps': [ 0.0086] * 8,
    'passbands': n.array([
        [-21.60550147235115, -2096.8027419993487, 820.28013290379204, -106.12782554677528, 5.4544068361089124],
        [-191.91156233339871, 166.41682897705232, -379.92321423032604, 102.73266064450996, -6.1221743594964124],
        [19.814315307496415, -3078.8760681203639, 1179.058479253717, -146.75025988148019, 6.9753639259636344],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ]),
    'beam': a.fit.BeamPolynomial,
    'bm_poly': n.reshape(n.array(
        [3.2706543329899773, -30.060111820705359, 94.699043700476409, 0.25922151678704664, -2.8449917420421791, 10.289468327128318, 0.50420327063121784, -5.6288445551201889, 13.4755446626395, 0.97625132261788883, -12.014627975436618, 35.579699435769172, 0.10747847373840302, -1.8009932273709515, 6.9210556027805019, -0.071016615227851854, 0.30787730434582339, 0.96080735693237029]
    ), (6,3)),
    'aa': {
        'xa_259':
            [0,0,0.18066613425943423, 0.47192565284191368, 0.4180807403246809, -0.081774445102069979],
        'xp_259':
            [0,0,-0.25701327768400983, -0.37646525073375664, -0.50381458302928495, -0.034605410415949445],

    }
}

def get_aa(freqs):
    '''Return the AntennaArray to be used fro simulation.'''
    beam = prms['beam'](freqs)
    try: beam.set_params(prms)
    except(AttributeError): pass
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    assert(len(prms['delays']) == nants and len(prms['offsets']) == nants \
        and len(prms['amps']) == nants and len(prms['passbands']) == nants)
    for pos, dly, off, amp, bp in zip(prms['antpos'], prms['delays'],
            prms['offsets'], prms['amps'], prms['passbands']):
        antennas.append(
            a.fit.Antenna(pos[0],pos[1],pos[2], beam, delay=dly, offset=off,
                amp=amp, bp=bp)
        )
    aa = XTalkAntennaArray(prms['loc'], antennas)
    aa.set_params(prms['aa'])
    return aa

