"""
A package for interfacing to Miriad: the Multichannel Image Reconstruction
Image Analysis and Display package for reducing interferometric data for
radio telescopes.
"""

__version__ = '0.1.1'

import numpy as np, _miriad

def echo(uv, p, d): return p, d

str2pol = {
    'I' :  1,   # Stokes Paremeters
    'Q' :  2,
    'U' :  3,
    'V' :  4,
    'rr': -1,   # Circular Polarizations
    'll': -2,
    'rl': -3,
    'lr': -4,
    'xx': -5,   # Linear Polarizations
    'yy': -6,
    'xy': -7,
    'yx': -8,
}

pol2str = {}
for k in str2pol: pol2str[str2pol[k]] = k

itemtable = {
    'obstype' : 'a',
    'history' : 'a',
    'vartable': 'a',
    #'visdata' : '?',
    #'flags'   : 'i',
    #'wflags'  : 'i',
    #'gains'   : '?',
    'ngains'  : 'i',
    'nfeeds'  : 'i',
    'ntau'    : 'i',
    'nsols'   : 'i',
    'interval': 'd',
    'leakage' : 'c',
    'freq0'   : 'd',
    'freqs'   : '?',
    'bandpass': 'c',
    'nspect0' : 'i',
    'nchan0'  : 'i',
}

data_types = {
    'a' : "ascii (NULL terminated)",
    'r' : "real (32 bin IEEE)",
    'd' : "double (64 bit)",
    'c' : "complex (2 * 32 bit IEEE)",
    'i' : "integer (32 bit two's complement)",
}

#  _   ___     __
# | | | \ \   / /
# | | | |\ \ / / 
# | |_| | \ V /  
#  \___/   \_/   

class UV(_miriad.UV):
    """Top-level interface to a Miriad UV data set."""
    def __init__(self, filename, status='old', corrmode='r'):
        """Open a miriad file.  status can be ('old','new','append').  
        corrmode can be 'r' (float32 data storage) or 'j' (int16 with shared exponent).  Default is 'r'."""
        assert(status in ['old', 'new', 'append'])
        assert(corrmode in ['r', 'j'])
        _miriad.UV.__init__(self, filename, status, corrmode)
        self.status = status
        self.nchan = 4096
        if status == 'old':
            self.vartable = self._gen_vartable()
            self.read(); self.rewind() # Update variables for the user
            try: self.nchan = self['nchan']
            except(KeyError): pass
        else: self.vartable = {'corr':corrmode}
    def _gen_vartable(self):
        """Generate table of variables and types from the vartable header."""
        vartable = {}
        for line in self._rdhd('vartable').split('\n'):
            try:
                type, name = line.split()
                vartable[name] = type
            except(ValueError): pass
        return vartable
    def vars(self):
        """Return a list of available variables."""
        return self.vartable.keys()
    def items(self):
        """Return a list of available header items."""
        items = []
        for i in itemtable:
            try:
                _miriad.hdaccess(self.haccess(i, 'read'))
                items.append(i)
            except(IOError): pass
        return items
    def _rdhd(self, name):
        """Provide read access to header items via low-level calls."""
        itype = itemtable[name]
        if itype == '?': return self._rdhd_special(name)
        h = self.haccess(name, 'read')
        rv = []
        if len(itype) == 1:
            if itype == 'a': offset = 0
            else:
                t, offset = _miriad.hread_init(h)
                assert(itype == t)
            while True:
                try: c, o = _miriad.hread(h, offset, itype)
                except(IOError): break
                if itype == 'a': c = c[:o]
                rv.append(c)
                offset += o
            if itype == 'a': rv = ''.join(rv)
        else:
            t, offset = _miriad.hread_init(h); assert(t == 'b')
            for t in itype:
                v, o = _miread.hread(h, offset, t)
                rv.append(v); offset += o
        _miriad.hdaccess(h)
        if len(rv) == 1: return rv[0]
        elif type(rv) == str: return rv
        else: return np.array(rv)
    def _wrhd(self, name, val):
        """Provide write access to header items via low-level calls."""
        type = itemtable[name]
        if type == '?': return self._wrhd_special(name, val)
        h = self.haccess(name, 'write')
        if len(type) == 1:
            try: len(val)
            except(TypeError): val = [val]
            if type == 'a': offset = 0
            else: offset = _miriad.hwrite_init(h, type)
            for v in val: offset += _miriad.hwrite(h, offset, v, type)
        else:
            offset = _miriad.hwrite_init(h, 'b')
            for v, t in zip(val,type): offset += _miriad.hwrite(h,offset,v,t)
        _miriad.hdaccess(h)
    def _rdhd_special(self, name):
        """Provide read access to special header items of type '?' to _rdhd"""
        if name == 'freqs':
            h = self.haccess(name, 'read')
            c, o = _miriad.hread(h, 0, 'i')
            rv = [c]; offset = 8
            while True:
                try:
                    c, o = _miriad.hread(h, offset, 'i')
                    rv.append(c); offset += 8
                    c, o = _miriad.hread(h, offset, 'd')
                    rv.append(c); offset += 8
                    c, o = _miriad.hread(h, offset, 'd')
                    rv.append(c); offset += 8
                except(IOError): break
            _miriad.hdaccess(h)
            return rv
        else: raise ValueError('Unknown special header: ' + name)
    def _wrhd_special(self, name, val):
        """Provide write access to special header items of type '?' to _rdhd"""
        if name == 'freqs':
            h = self.haccess(name, 'write')
            o = _miriad.hwrite(h, 0, val[0], 'i')
            offset = 8
            for i,v in enumerate(val[1:]):
                if i % 3 == 0: o = _miriad.hwrite(h, offset, v, 'i')
                else: o = _miriad.hwrite(h, offset, v, 'd')
                offset += 8
            _miriad.hdaccess(h)
        else: raise ValueError('Unknown special header: ' + name)
    def __getitem__(self, name):
        """Allow access to variables and header items via uv[name]."""
        try:
            type = self.vartable[name]
            return self._rdvr(name, type)
        except(KeyError):
            type = itemtable[name]
            return self._rdhd(name)
    def __setitem__(self, name, val):
        """Allow setting variables and header items via uv[name] = val."""
        try:
            type = self.vartable[name]
            self._wrvr(name,type,val)
        except(KeyError):
            self._wrhd(name,val)
    def select(self, name, n1, n2, include=1):
        """Choose which data are returned by read().  
            name    This can be: 'decimate','time','antennae','visibility',
                    'uvrange','pointing','amplitude','window','or','dra',
                    'ddec','uvnrange','increment','ra','dec','and', 'clear',
                    'on','polarization','shadow','auto','dazim','delev'
            n1,n2   Generally this is the range of values to select. For
                    'antennae', this is the two antennae pair to select
                    (indexed from 0); a -1 indicates 'all antennae'.
                    For 'decimate', n1 is every Nth integration to use, and
                    n2 is which integration within a block of N to use.
                    For 'shadow', a zero indicates use 'antdiam' variable.
                    For 'on','window','polarization','increment','shadow' only
                    p1 is used.
                    For 'and','or','clear','auto' p1 and p2 are ignored.
            include If true, the data is selected. If false, the data is
                    discarded. Ignored for 'and','or','clear'."""
        if name == 'antennae':
            n1 += 1; n2 += 1
        self._select(name, float(n1), float(n2), int(include))
    def read(self, raw=False):
        """Return the next data record.  Calling this function causes 
        vars to change to reflect the record which this function returns.
        'raw' causes data and flags to be returned seperately."""
        preamble, data, flags, nread = self.raw_read(self.nchan)
        if nread == 0: raise IOError("No data read")
        flags = np.logical_not(flags)
        if raw: return preamble, data, flags
        return preamble, np.ma.array(data, mask=flags)
    def all(self, raw=False):
        """Provide an iterator over preamble, data.  Allows constructs like: 
        for preamble, data in uv.all(): ..."""
        curtime = None
        while True:
            try: yield self.read(raw=raw)
            except(IOError): return
    def write(self, preamble, data, flags=None):
        """Write the next data record.  data must be a complex, masked
        array.  preamble must be (uvw, t, (i,j)), where uvw is an array of 
        u,v,w, t is the Julian date, and (i,j) is an antenna pair."""
        if data is None: return
        if not flags is None: flags = np.logical_not(flags)
        elif len(data.mask.shape) == 0:
            flags = np.ones(data.shape)
            data = data.unmask()
        else:
            flags = np.logical_not(data.mask)
            #data = data.filled(0)
            data = data.data
        self.raw_write(preamble,data.astype(np.complex64),flags.astype(np.int32))
    def init_from_uv(self, uv, override={}, exclude=[]):
        """Initialize header items and variables from another UV.  Those in 
        override will be overwritten by override[k], and tracking will be 
        turned off (meaning they will not be updated in pipe()).  Those in
        exclude are omitted completely."""
        for k in uv.items():
            if k in exclude: continue
            elif k in override: self._wrhd(k, override[k])
            else: self._wrhd(k, uv[k])
        self.vartable = {}
        for k in uv.vars():
            if k in exclude: continue
            # I don't understand why reading 'corr' segfaults miriad,
            # but it does.  This is a cludgy work-around.
            elif k == 'corr': continue
            elif k in override:
                self.vartable[k] = uv.vartable[k]
                self._wrvr(k, uv.vartable[k], override[k])
            else:
                self.vartable[k] = uv.vartable[k]
                self._wrvr(k, uv.vartable[k], uv[k])
                uv.trackvr(k, 'c') # Set to copy when copyvr() called 
    def pipe(self, uv, mfunc=echo, append2hist='', raw=False):
        """Pipe in data from another UV through the function
        mfunc(uv,preamble,data), which should return (preamble,data).  If 
        mfunc is not provided, the dataset will just be cloned, and if the 
        returned data is None, it will be omitted.  The string 'append2hist' 
        will be appended to history."""
        self._wrhd('history', self['history'] + append2hist)
        # Pipe all data through mfunc
        if raw:
            for p,d,f in uv.all(raw=raw):
                np, nd, nf = mfunc(uv, p, d, f)
                self.copyvr(uv)
                self.write(np, nd, nf)
        else:
            for p, d in uv.all():
                np, nd = mfunc(uv, p, d)
                self.copyvr(uv)
                self.write(np, nd)
    def add_var(self, name, type):
        """Add a variable of the specified type to a UV file."""
        self.vartable[name] = type

def bl2ij(bl):
    bl = int(bl)
    if (bl > 65536):
        bl -= 65536
        mant = 2048
    else:
        mant = 256
#AAR    return (bl>>8)-1, (bl&255) - 1
    return bl/mant - 1, bl%mant -1

def ij2bl(i, j):
    if i > j: i,j = j,i
    if j + 1 < 256: return 256*(i+1) + (j+1)
    else: return 2048*(i+1) + (j+1) + 65536

