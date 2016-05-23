import numpy as np
from aipy._cephes import i0
from _dsp import *

WINDOW_FUNC = {
    'blackman': lambda x,L: .42-.5*np.cos(2*np.pi*x/(L-1))+.08*np.cos(4*np.pi*x/(L-1)),
    'blackman-harris': lambda x,L: .35875 - .48829*np.cos(2*np.pi*x/(L-1)) + .14128*np.cos(4*np.pi*x/(L-1)) - .01168*np.cos(6*np.pi*x/(L-1)),
    'gaussian0.4': lambda x,L: np.exp(-0.5 * ((x - (L-1)/2)/(0.4 * (L-1)/2))**2),
    'kaiser2': lambda x,L: i0(np.pi * 2 * np.sqrt(1-(2*x/(L-1) - 1)**2)) / i0(np.pi * 2),
    'kaiser3': lambda x,L: i0(np.pi * 3 * np.sqrt(1-(2*x/(L-1) - 1)**2)) / i0(np.pi * 3),
    'hamming': lambda x,L: .54 - .46 * np.cos(2*np.pi*x/(L-1)),
    'hanning': lambda x,L: .5 - .5 * np.cos(2*np.pi*x/(L-1)),
    'parzen': lambda x,L: 1 - np.abs(L/2. - x) / (L/2.),
    'none': lambda x,L: 1,
}

def gen_window(L, window='hamming'):
    '''Return the specified window (see WINDOW_FUNC) for a length L.'''
    return np.fromfunction(lambda x: WINDOW_FUNC[window](x,L), (L,))
