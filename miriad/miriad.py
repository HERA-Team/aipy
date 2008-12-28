#! /usr/bin/env python

"""
A package for interfacing to Miriad: the Multichannel Image Reconstruction
Image Analysis and Display package for reducing interferometric data for
radio telescopes.

Author: Aaron Parsons
Date: 10/08/06
Revisions:
    arp 10/25/06    Moved miruv.c functions to seperate module 'miruv'.
    arp 10/26/06    Added map_uv to easily clone datasets.  Found some bugs.
    arp 11/02/06    Added support for arbitrary length preambles and 
                    functionality for changing the preamble format.
    arp 11/10/06    Modified map_uv to be able to read chucks of time data,
                    rather than 1 spectrum at a time.
    arp 12/06/06    Added idiom of uv.read_data(); uv.rewind() to init so
                    that vars get initialized.
    arp 01/21/07    Changed buffering of data, flags, preamble in read_data
                    to cure a memory leak, and (possibly) speed reading.
    arp 01/28/07    Added uvupdate, uvselect, and uvtrack to speed up copying
                    and mapping of datasets.
    arp 03/13/07    Changed "map_uv" to "pipe_uv" and removed buffering in
                    it to fix a bug and simplify what is no longer needed
                    now that sim_data only generates single baselines.

Known Issues:
    You have to del(uv) in ipython, or vartable doesn't get written.
    UVVar 'corr' segfaults miriad

Todo:
    Implement uvupdate_c,  uvselect_c, uvscan_c for only reading out the data 
        you want. 
    Implement uvcopyvr_c to speed cloning of datasets.
    Work with bug_c to produce python exceptions instead of sys exits.
"""
__version__ = '0.0.4'

import numpy, miruv

MAX_CHANNELS = 4096
HISTORY_LINE_LEN = 80
MAX_VAR_LEN = 1024
MAX_VARTABLE_LEN = 1024
MAX_PREAMBLE = 9            # This value comes from uvio.c

class TypeCode:
    """Class for interfacing Miriad UV variable types with python."""
    def __init__(self, mircode, bytes, numpytype):
        self.mircode = mircode
        self.bytes = bytes
        self.numpytype = numpytype
    def unpack(self, binstr):
        """Create a numpy array of the binary data read from a UV variable."""
        if self.numpytype is numpy.character:
            return numpy.array(binstr, dtype=self.numpytype)
        return numpy.fromstring(binstr, dtype=self.numpytype)
    def pack(self, *args):
        """Encode a sequence of numbers as a binary string."""
        a = numpy.array(args, dtype=self.numpytype)
        return a.tostring()

class AsciiTypeCode(TypeCode):
    def unpack(self, binstr):
        return binstr
    def pack(self, s):
        return s
        
typecodes = {'a':AsciiTypeCode(  1,  1,  numpy.character),
             'r':TypeCode(  4,  4,  numpy.float32),
             'd':TypeCode(  5,  8,  numpy.float64),
             'c':TypeCode(  6,  8,  numpy.complex64),   # This might be wrong
             'i':TypeCode(  2,  4,  numpy.int32),
             'j':TypeCode(  3,  2,  numpy.int16),       # This segfaults
             None:None,}

item_types = {'obstype' : 'a',
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
         #'freqs'   : '?',
         'bandpass': 'c',
         'nspect0' : 'i',
         'nchan0'  : 'i',}

#  _   ___     _____ _                 
# | | | \ \   / /_ _| |_ ___ _ __ ___  
# | | | |\ \ / / | || __/ _ \ '_ ` _ \ 
# | |_| | \ V /  | || ||  __/ | | | | |
#  \___/   \_/  |___|\__\___|_| |_| |_|

class UVItem:
    """Class for reading/writing UV items."""
    def __init__(self, name, uvhandle, typestr):
        self.typestr = typestr
        self.name = name
        if typestr == 'r':
            self.read = lambda: miruv.rdhdr_c(uvhandle, self.name, 0.)
            self.write = lambda x: miruv.wrhdr_c(uvhandle, self.name, x)
        if typestr == 'd':
            self.read = lambda: miruv.rdhdd_c(uvhandle, self.name, 0.)
            self.write = lambda x: miruv.wrhdd_c(uvhandle, self.name, x)
        if typestr == 'i':
            self.read = lambda: miruv.rdhdi_c(uvhandle, self.name, 0)
            self.write = lambda x: miruv.wrhdi_c(uvhandle, self.name, x)
        if typestr == 'c':
            def read():
                n = miruv.hsize_c_wrap(uvhandle, self.name)
                d = numpy.zeros(n, dtype=numpy.float32)
                n = miruv.rdhdc_c_wrap(uvhandle, self.name, d)
                d = d[:n]
                d.shape = (d.size/2, 2)
                return d[:,0] + d[:,1]*1j
            self.read = read
            def write(x):
                nx = numpy.zeros((x.size, 2), dtype=numpy.float32)
                nx[:,0] = x.real; nx[:,1] = x.imag
                miruv.wrhdc_c_wrap(uvhandle, self.name, nx.flatten())
            self.write = write
        if typestr == 'a':
            def read():
                handle, status = miruv.haccess_c(uvhandle, self.name, 'read')
                if status != 0: return ''
                val = []
                while status == 0:
                    s, status = miruv.hreada_c(handle, MAX_VARTABLE_LEN)
                    if status == 0: val.append(s)
                miruv.hdaccess_c(handle)
                return '\n'.join(val)
            self.read = read
            def write(x):
                handle, status = miruv.haccess_c(uvhandle, self.name, 'write')
                if status != 0: return
                status = miruv.hwritea_c(handle, x, len(x))
                miruv.hdaccess_c(handle)
            self.write = write

#  _   ___     _____ _               _____     _     _      
# | | | \ \   / /_ _| |_ ___ _ __ __|_   _|_ _| |__ | | ___ 
# | | | |\ \ / / | || __/ _ \ '_ ` _ \| |/ _` | '_ \| |/ _ \
# | |_| | \ V /  | || ||  __/ | | | | | | (_| | |_) | |  __/
#  \___/   \_/  |___|\__\___|_| |_| |_|_|\__,_|_.__/|_|\___|

class UVItemTable(dict):
    def __init__(self, uvhandle):
        dict.__init__(self)
        self.uvhandle = uvhandle
        for i in item_types:
            uvi = UVItem(i, uvhandle, item_types[i])
            dict.__setitem__(self, i, uvi)
    def __getitem__(self, k):
        uvi = dict.__getitem__(self, k)
        return uvi.read()
    def __setitem__(self, k, val):
        uvi = self.get(k, UVItem(k, self.uvhandle, item_types[k]))
        uvi.write(val)
        dict.__setitem__(self, k, uvi)
    def gen_uvvartable(self, status):
        vt = UVVarTable(self.uvhandle, status)
        try: vartablestr = self['vartable']
        except(KeyError): return vt
        for v in vartablestr.splitlines():
            typecode = v[0]
            if typecode not in typecodes.keys(): continue
            name = v[2:]
            vt.add_var(name, typecode=typecode)
        return vt

#  _   ___     ____     __         
# | | | \ \   / /\ \   / /_ _ _ __ 
# | | | |\ \ / /  \ \ / / _` | '__|
# | |_| | \ V /    \ V / (_| | |   
#  \___/   \_/      \_/ \__,_|_|   

class UVVar:
    """Class for keeping track of UV variables as we step through Miriad 
    data."""
    def __init__(self, name, uvhandle, typecode=None):
        self.uvhandle = uvhandle
        self.name = name
        self.typecode = typecode
        self.tc = typecodes[typecode]
        self.initialized = False
    def changed(self):
        """Probes to see if this variable changed on the last read.  Meanwhile, 
        it uses the information it gathers to initialize the variable."""
        T, L, U = miruv.uvprobvr_c(self.uvhandle, self.name, ' ')
        if T == ' ':
            raise ValueError('UVVar "%s" does not exist.' % (self.name))
        elif L == 0:
            raise ValueError('UVVar "%s" currently has no value.' % (self.name))
        self.tc = typecodes[T]
        self.len = L
        # Because character variables have a terminating null, we have to
        # ask for one more byte when reading them.
        if T == 'a': self.len += 1
        self.initialized = True
        return U
    def read(self):
        """Return the value of the variable.  Initializes the variable if
        'changed' hasn't been called already."""
        if not self.initialized: self.changed()
        binstr = miruv.uvgetvr_c_wrap(self.uvhandle, self.tc.mircode,
            self.name, MAX_VAR_LEN, self.len, self.tc.bytes)
        val = self.tc.unpack(binstr)
        if type(val) is not str and val.size == 1: val = val[0]
        return val
    def write(self, val):
        """Write a value for a uv variable to a dataset.  I am guessing that
        you better have opened your uv dataset in 'write' mode..."""
        binstr = self.tc.pack(val)
        miruv.uvputvr_c(self.uvhandle, self.tc.mircode, self.name, binstr, 
            len(binstr) / self.tc.bytes)

#  _   ___     ____     __        _____     _     _      
# | | | \ \   / /\ \   / /_ _ _ _|_   _|_ _| |__ | | ___ 
# | | | |\ \ / /  \ \ / / _` | '__|| |/ _` | '_ \| |/ _ \
# | |_| | \ V /    \ V / (_| | |   | | (_| | |_) | |  __/
#  \___/   \_/      \_/ \__,_|_|   |_|\__,_|_.__/|_|\___|

class UVVarTable(dict):
    def __init__(self, uvhandle, status):
        dict.__init__(self)
        self.uvhandle = uvhandle
        self.status = status
    def add_var(self, name, typecode):
        dict.__setitem__(self, name, UVVar(name, self.uvhandle, typecode))
    def __getitem__(self, k):
        v = dict.__getitem__(self, k)
        return v.read()
    def __setitem__(self, k, val):
        if self.status == 'old':
            raise IOError('Cannot write to UV in %s mode.' % (self.status))
        elif not self.has_key(k):
            raise KeyError('Key "%s" does not exist.' % (k) + 
                '  You can add this key to a UVVarTable with vartable.add_var')
        dict.__getitem__(self, k).write(val)
    def set_tracking(self, name, mode):
        """Adjust the tracking status of variable 'name' to 'mode'.  Valid
        switches are 'u' (you care if it is updated), 'c' (it should be
        copied to an output data set), 'uc' (both), or '' (none)."""
        miruv.uvtrack_c(self.uvhandle, name, mode)
    def changed(self):
        """Returns true if any variables in mode 'u' (see set_tracking)
        changed during the last read."""
        return miruv.uvupdate_c(self.uvhandle)

#  _   ___     __
# | | | \ \   / /
# | | | |\ \ / / 
# | |_| | \ V /  
#  \___/   \_/   

class UV:
    """Class for representing a Miriad UV data set."""
    def __init__(self, name, status='old'):
        """Open miriad dataset 'name'.  status can be 'old', 'new', or 
        'append'."""
        if status not in ('old', 'new', 'append'):
            raise ValueError('Unrecognized status: %s' % (status))
        self.status = status
        self.handle = miruv.uvopen_c(name, status)
        # Initialize to use the standard preamble mode: uvw/time/baseline
        self.configure_preamble()
        self.items = UVItemTable(self.handle)
        self.vars = self.items.gen_uvvartable(status)
        # Define buffers for holding data and flags which are reused on reads.
        self._data = numpy.zeros((2*MAX_CHANNELS,), dtype=numpy.float32)
        self._flags = numpy.zeros((MAX_CHANNELS,), dtype=numpy.int32)
        # Read 1 data and then rewind so Miriad will provide values for
        # basic variables right from the start.
        if status == 'old': self.read_data(); self.rewind()
    def configure_preamble(self, preamble_type='uvw/time/baseline'):
        """Configure which data are placed in the preamble returned
        with every read_data() and written with every write_data.
        Some common preamble_types are 'uv/time/baseline', 'uvw/time/baseline',
        'uvw/time/baseline/pol', etc."""
        # Calculate size of preamble:
        self.preamble_size = 0
        for w in preamble_type.split('/'):
            if w == 'uv': self.preamble_size += 2
            elif w == 'uvw': self.preamble_size += 3
            else: self.preamble_size += 1
        # Configure preamble readout mode
        miruv.uvset_c(self.handle, 'preamble', preamble_type, 0, 0., 0., 0.)
        self._preamble = numpy.zeros((self.preamble_size,), dtype=numpy.double)
    def select_data(self, name, p1, p2, include_it=1):
        """Choose which data are returned by 'read_data'.  
            name    This can be one of 'time','antennae','visibility',
                    'uvrange','pointing','amplitude','window','or','dra',
                    'ddec','uvnrange','increment','ra','dec','and', 'clear',
                    'on','polarization','shadow','auto','dazim','delev'
            p1,p2   Generally this is the range of values to select. For
                    'antennae', this is the two antennae pair to select
                    (indexed from 0); a -1 indicates 'all antennae'.
                    For 'shadow', a zero indicates use 'antdiam' variable.
                    For 'on','window','polarization','increment','shadow' only
                    p1 is used.
                    For 'and','or','clear','auto' p1 and p2 are ignored.
            include_it    If true, the data is selected. If false, the data is
                    discarded. Ignored for 'and','or','clear'."""
        if name == 'antennae':
            p1 += 1; p2 += 1
        miruv.uvselect_c(self.handle, str(name), float(p1), float(p2), 
            int(include_it))
    def rewind(self):
        """Reset read-point to beginning of dataset, so that UV.next_data() will
        return the first record."""
        miruv.uvrewind_c(self.handle)
    def read_data(self, use_mask=True):
        """Return the data in the next dataset record.  Calling this function
        causes the values in UVVars to change to reflect the record
        which this function returns."""
        if self.status != 'old':
            raise IOError('Cannot read from UV in %s mode.' % (self.status))
        nread = miruv.uvread_c_wrap(self.handle, self._preamble, 
            self._data, self._flags)
        data = self._data[:2*nread]; data.shape = (data.shape[0]/2, 2)
        data = data[:,0] + data[:,1] * 1j
        if use_mask:
            flags = self._flags[:nread]
            flags = numpy.logical_not(flags)
            data = numpy.ma.array(data, mask=flags)
        return self._preamble.copy(), data
    def write_data(self, preamble, data):
        """Write data as the next UV record.  'data' must be a complex, masked
        array.  'preamble' must be a numpy array of values associated with
        the preamble_type used in config_preamble()."""
        if self.status == 'old':
            raise IOError('Cannot write to UV in %s mode.' % (self.status))
        if preamble.size != self.preamble_size:
            raise ValueError('Expected preamble of length %d.  Got %d.' % \
                (self.preamble_size, preamble.size))
        flags = numpy.logical_not(data.mask).astype(numpy.int32)
        # Changed to use filled (even though it loses the info of what the
        # flagged data was) b/c data.data doesn't work...
        data = data.filled(0)
        fl_data = numpy.zeros((data.size, 2), dtype=numpy.float32)
        #fl_data[:,0] = data.data.real; fl_data[:,1] = data.data.imag
        fl_data[:,0] = data.real; fl_data[:,1] = data.imag
        fl_data.shape = (2*data.size,)
        if len(flags.shape) == 0:
            flags = numpy.ones(data.shape, dtype=numpy.int32)
        miruv.uvwrite_c_wrap(self.handle, preamble, fl_data, flags)
    def update_vars(self):
        """Refresh the list of variables available in this dataset."""
        self.vars.from_vartable(self.items['vartable'].read())
    def keys(self):
        return self.vars.keys() + self.items.keys()
    def __getitem__(self, k):
        if self.vars.has_key(k): return self.vars[k]
        elif self.items.has_key(k): return self.items[k]
        else: raise KeyError('%s not a key in vars or in items.' % (k))
    def __del__(self):
        miruv.uvclose_c(self.handle)

def bl2ij(bl):
    bl = int(bl)
    return (bl>>8)-1, (bl&255) - 1

def ij2bl(i, j):
    return (i+1)<<8 | (j+1)

#        _   _ _ _ _         
#  _   _| |_(_) (_) |_ _   _ 
# | | | | __| | | | __| | | |
# | |_| | |_| | | | |_| |_| |
#  \__,_|\__|_|_|_|\__|\__, |
#                      |___/ 

def init_from_uv(uvi, uvo, append2hist='', override={}, exclude=[]):
    """Initialize 'uvo' from 'uvi'.  Append 'append2hist' to uvo's history.""" 
    for k in uvi.items:
        if k in exclude: continue
        if k in override: uvo.items[k] = override[k]
        else: uvo.items[k] = uvi.items[k]
    uvo.items['history'] += append2hist
    for k in uvi.vars:
        # I don't understand why reading 'corr' segfaults miriad,
        # but it does.  This is a cludgy work-around.
        if k == 'corr': continue
        uvo.vars.add_var(k, dict.__getitem__(uvi.vars, k).typecode)
        if k in override: uvo.vars[k] = override[k]
        else:
            try: uvo.vars[k] = uvi.vars[k]
            except(ValueError): pass

def pipe_uv(uvi, uvo, mfunc=None, append2hist='', init=True, notrack=[],
        override={}, exclude=[]):
    """Pipe one UV dataset (uvi) into another (uvo) through the mapping function
    'mfunc' (if not provided, this will just clone a dataset).
    
    'mfunc' should be a function of 3 args: (uv, preamble, data),
    and should return (preamble, data).  If data is None, it will be
    omitted from 'uvo'.
    'uvi' and 'uvo' should be an opened, 
    but otherwise untouched uv files, unless 'init' is False.
    'append2hist' is a string to add to the history of uvo, if you want.
    'init' will copy the initial values of all items and vars in uvi to
    uvo.  If you don't set this True, then it is your responsibility to
    initialize uvo."""
    # Optionally initialize uvo with appropriate variables, items
    if init:
        init_from_uv(uvi, uvo, append2hist=append2hist, override=override,
            exclude=exclude)
    # Set up uvi to copy all variables when 'uvcopyvr' is called.
    for k in uvi.vars:
        if k == 'corr': continue        # Cludge again --- yuck.
        if k not in notrack: uvi.vars.set_tracking(k, 'c')
    # Pipe all data through mfunc to uvo
    while True:
        p, d = uvi.read_data()
        if d.size == 0: break
        if not mfunc is None:
            np, nd = mfunc(uvi, p, d)
            if nd is None: continue
        else: np, nd = p, d
        miruv.uvcopyvr_c(uvi.handle, uvo.handle)
        uvo.write_data(np, nd)

#  _            _   _                     _     
# | |_ ___  ___| |_| |__   ___ _ __   ___| |__  
# | __/ _ \/ __| __| '_ \ / _ \ '_ \ / __| '_ \ 
# | ||  __/\__ \ |_| |_) |  __/ | | | (__| | | |
#  \__\___||___/\__|_.__/ \___|_| |_|\___|_| |_|

if __name__ == '__main__':
    import os, sys

    if os.path.exists(sys.argv[1]): os.system('rm -r %s' % (sys.argv[1]))

    # Create a test UV file
    time = 12345.
    HIST_STR = 'A test history.'
    BANDPASS = numpy.array([1+1j,2+2j,3+3j], dtype=numpy.complex64)
    uv = UV(sys.argv[1], status='new')
    uv.items['history'] = HIST_STR
    uv.items['ngains'] = 1
    uv.items['bandpass'] = BANDPASS
    a = numpy.ma.array([1+1j, 2+2j, 3+3j], mask=[0, 0, 1], 
        dtype=numpy.complex64)
    uv.vars.add_var('pol', 'i')
    for t in range(100):
      for i in range(8):
        for j in range(8):
            uv.vars['pol'] = i
            bl = float((i+1)<<8 + j+1)
            preamble = numpy.array((1., 1., 1., time, bl), dtype=numpy.double)
            uv.write_data(preamble, a)
      time += 10
    del(uv)

    # Make sure you can read everything from the UV file.
    time = 12345.
    uv = UV(sys.argv[1], status='old')
    assert(uv.items['history'] == HIST_STR)
    assert(uv.items['ngains'] == 1)
    assert(numpy.all(uv.items['bandpass'] == BANDPASS))
    for t in range(100):
      for i in range(8):
        for j in range(8):
            bl = float((i+1)<<8 + j+1)
            preamble = numpy.array((1., 1., 1., time, bl), dtype=numpy.double)
            p, d = uv.read_data()
            if numpy.any(d.filled() != a.filled()) or numpy.any(p != preamble):
                raise ValueError('Data was not written correctly.')
            if uv.vars['pol'] != i: raise ValueError('UVVar read/write failed.')
      time += 10
    del(uv)

    print 'Passed test.'
