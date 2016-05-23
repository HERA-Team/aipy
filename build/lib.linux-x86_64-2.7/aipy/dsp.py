import numpy as n
from aipy._cephes import i0
from _dsp import *

WINDOW_FUNC = {
    'blackman': lambda x,L: .42-.5*n.cos(2*n.pi*x/(L-1))+.08*n.cos(4*n.pi*x/(L-1)),
    'blackman-harris': lambda x,L: .35875 - .48829*n.cos(2*n.pi*x/(L-1)) + .14128*n.cos(4*n.pi*x/(L-1)) - .01168*n.cos(6*n.pi*x/(L-1)),
    'gaussian0.4': lambda x,L: n.exp(-0.5 * ((x - (L-1)/2)/(0.4 * (L-1)/2))**2),
    'kaiser2': lambda x,L: i0(n.pi * 2 * n.sqrt(1-(2*x/(L-1) - 1)**2)) / i0(n.pi * 2),
    'kaiser3': lambda x,L: i0(n.pi * 3 * n.sqrt(1-(2*x/(L-1) - 1)**2)) / i0(n.pi * 3),
    'hamming': lambda x,L: .54 - .46 * n.cos(2*n.pi*x/(L-1)),
    'hanning': lambda x,L: .5 - .5 * n.cos(2*n.pi*x/(L-1)),
    'parzen': lambda x,L: 1 - n.abs(L/2. - x) / (L/2.),
    'none': lambda x,L: 1,
}

def gen_window(L, window='hamming'):
    '''Return the specified window (see WINDOW_FUNC) for a length L.'''
    return n.fromfunction(lambda x: WINDOW_FUNC[window](x,L), (L,))
