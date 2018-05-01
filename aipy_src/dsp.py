import numpy as np
from aipy._cephes import i0
from _dsp import *


def tukey(x, L, alpha=0.5):
    """A Tukey, or tapered-cosine, window. See scipy.signal.windows for details. 
    alpha=0 is rectangle, alpha=1 is hann window."""
    if alpha <= 0:
        return np.ones(L, 'd')
    elif alpha >= 1.0:
        return WINDOW_FUNC['hanning'](x, L)

    width = int(np.floor(alpha*(L-1)/2.0))
    x1 = x[0:width+1]
    x2 = x[width+1:L-width-1]
    x3 = x[L-width-1:]

    w1 = 0.5 * (1 + np.cos(np.pi * (-1 + 2.0*x1/alpha/(L-1))))
    w2 = np.ones(x2.shape)
    w3 = 0.5 * (1 + np.cos(np.pi * (-2.0/alpha + 1 + 2.0*x3/alpha/(L-1))))

    w = np.concatenate((w1, w2, w3))

    return w

def barthann(x, L):
    """A Bartlett-Hann window. See scipy.signal.windows for details """
    x = np.arange(0, L)
    fac = np.abs(x / (L - 1.0) - 0.5)
    w = 0.62 - 0.48 * fac + 0.38 * np.cos(2 * np.pi * fac)

    return w

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
    'tukey': tukey,
    'barthann': barthann
    }

def gen_window(L, window='hamming', **kwargs):
    '''Return the specified window (see WINDOW_FUNC) for a length L.'''
    return np.fromfunction(lambda x: WINDOW_FUNC[window](x,L,**kwargs), (L,))
