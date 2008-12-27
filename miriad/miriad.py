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

Known Issues:
    You have to del(uv) in ipython, or vartable doesn't get written.
    UVVar 'corr' segfaults miriad

Todo:
    Implement uvupdate_c,  uvselect_c, uvscan_c for only reading out the data 
        you want. 
    Implement uvcopyvr_c to speed cloning of datasets.
    Work with bug_c to produce python exceptions instead of sys exits.
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

import numpy, miruv

MAX_CHANNELS = 4096
HISTORY_LINE_LEN = 80
MAX_VAR_LEN = 1024
MAX_VARTABLE_LEN = 1024
MAX_PREAMBLE = 9            # This value comes from uvio.c

#  _____                  ____          _      
# |_   _|   _ _ __   ___ / ___|___   __| | ___ 
#   | || | | | '_ \ / _ \ |   / _ \ / _` |/ _ \
#   | || |_| | |_) |  __/ |__| (_) | (_| |  __/
#   |_| \__, | .__/ \___|\____\___/ \__,_|\___|
#       |___/|_|      

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
    def __init__(self, name, uvhandle, typecode):
        self.uvhandle = uvhandle
        self.name = name
        self.tc = typecodes[typecode]
    def _access(self, mode):
        """Acceptable modes are 'read', 'write', 'append', and 'scratch'."""
        self.handle, status = miruv.haccess_c(self.uvhandle, self.name, mode)
        if status != 0:
            raise ValueError('miruv.haccess_c failed for "%s".' % (self.name))
        self.curmode = mode
    def _daccess(self):
        if self.curmode is not None: miruv.hdaccess_c(self.handle)
        self.curmode = None
    def read(self):
        """Read the (semi-static) value of this item from the UV dataset."""
        self._access('read')
        val = []; status = 0
        while status == 0:
            s, status = miruv.hreada_c(self.handle, MAX_VARTABLE_LEN)
            if status == 0: val.append(self.tc.unpack(s))
        self._daccess()
        if len(val) == 0: return ''
        elif type(val[0]) is str: return '\n'.join(val)
        elif len(val) == 1:
            val = val[0]
            if val.size == 1: val = val[0]
        return val
    def write(self, val):
        """Overwrite the value of this item from the UV dataset."""
        self._access('write')
        binstr = self.tc.pack(val)
        status = miruv.hwritea_c(self.handle, binstr, len(binstr))
        if status != 0:
            raise IOError('Write failed for miruv.hwritea_c.')
        self._daccess()
    def append(self, val):
        """Append to the value of this item in the UV dataset."""
        self._access('append')
        binstr = self.tc.pack(val)
        status = miruv.hwritea_c(self.handle, binstr, len(binstr))
        if status != 0:
            raise IOError('Write failed for miruv.hwritea_c.')
        self._daccess()

#  _   ___     _____ _               _____     _     _      
# | | | \ \   / /_ _| |_ ___ _ __ __|_   _|_ _| |__ | | ___ 
# | | | |\ \ / / | || __/ _ \ '_ ` _ \| |/ _` | '_ \| |/ _ \
# | |_| | \ V /  | || ||  __/ | | | | | | (_| | |_) | |  __/
#  \___/   \_/  |___|\__\___|_| |_| |_|_|\__,_|_.__/|_|\___|

class UVItemTable(dict):
    def __init__(self, uvhandle, status):
        dict.__init__(self)
        self.status = status
        self.uvhandle = uvhandle
        for i in item_types:
            uvi = UVItem(i, uvhandle, item_types[i])
            try: 
                uvi._access('read')
                uvi._daccess()
                dict.__setitem__(self, i, uvi)
            except(ValueError): pass
    def __getitem__(self, k):
        uvi = dict.__getitem__(self, k)
        return uvi.read()
    def __setitem__(self, k, val):
        uvi = self.get(k, UVItem(k, self.uvhandle, item_types[k]))
        uvi.write(val)
        dict.__setitem__(self, k, uvi)
    def gen_uvvartable(self):
        vt = UVVarTable(self.uvhandle, self.status)
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
        self.items = UVItemTable(self.handle, status)
        self.vars = self.items.gen_uvvartable()
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
                    'antennae', this is the two antennae pair to select.
                    For 'antennae', a zero indicates 'all antennae'.
                    For 'shadow', a zero indicates use 'antdiam' variable.
                    For 'on','window','polarization','increment','shadow' only
                    p1 is used.
                    For 'and','or','clear','auto' p1 and p2 are ignored.
            include_it    If true, the data is selected. If false, the data is
                    discarded. Ignored for 'and','or','clear'."""
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
        fl_data = numpy.zeros((data.size, 2), dtype=numpy.float32)
        fl_data[:,0] = data.data.real; fl_data[:,1] = data.data.imag
        fl_data.shape = (2*data.size,)
        flags = numpy.logical_not(data.mask).astype(numpy.int32)
        miruv.uvwrite_c_wrap(self.handle, preamble, fl_data, flags)
    def update_vars(self):
        """Refresh the list of variables available in this dataset."""
        self.vars.from_vartable(self.items['vartable'].read())
    def __getitem__(self, k):
        if self.vars.has_key(k): return self.vars[k]
        elif self.items.has_key(k): return self.items[k]
        else: raise KeyError('%s not a key in vars or in items.' % (k))
    def __del__(self):
        miruv.uvclose_c(self.handle)

#        _   _ _ _ _         
#  _   _| |_(_) (_) |_ _   _ 
# | | | | __| | | | __| | | |
# | |_| | |_| | | | |_| |_| |
#  \__,_|\__|_|_|_|\__|\__, |
#                      |___/ 

def map_uv(uvi, uvo, mfunc=None, append2history='', send_time_blks=False,
        initialize=True, bufvars=[]):
    """Map one UV dataset (uvi) into another (uvo) through the mapping function
    'mfunc' (if not provided, this will just clone a dataset).
    
    'mfunc' should be a function of 4 args: (uv, preamble, data, vars),
    and should return (preamble, data).  If 'send_time_blks' is true, 
    each of these variables will be a list with each item corresponding to
    a read from 'uvi'.  In this case, the return values should also be lists.
    'uvi' and 'uvo' should be an opened, but otherwise untouched uv files.
    'append2history' is a string to add to the history of uvo, if you want.
    'initialize' will copy the initial values of all items and vars in uvi to
    uvo.  If you don't set this true, then it is your responsibility to
    initialize uvo.  'bufvars' is a list of additional variables you want to
    keep track of (this is because in send_time_blks mode, the uv file passed
    to your mfunc may not tell you all you need to know)."""
    if initialize:
        for k in uvi.items: uvo.items[k] = uvi.items[k]
        uvo.items['history'] += append2history
        for k in uvi.vars:
            # I don't understand why reading 'corr' segfaults miriad,
            # but it does.  This is a cludgy work-around.
            if k == 'corr': continue
            uvo.vars.add_var(k, dict.__getitem__(uvi.vars, k).typecode)
            try: uvo.vars[k] = uvi.vars[k]
            except(ValueError): pass
    for k in uvi.vars:
        # Set up uvi to copy all variables when 'uvcopyvr' is called.
        if k == 'corr': continue        # Cludge again --- yuck.
        uvi.vars.set_tracking(k, 'c')
    if send_time_blks: uvi.vars.set_tracking('time', 'uc')
    pbuf, dbuf, vbuf = [], [], {}
    for b in bufvars: vbuf[b] = []
    while True:
        p, d = uvi.read_data()
        if d.size == 0:
            # Send any remaining buffered data
            _process_buf(uvi, uvo, pbuf, dbuf, vbuf, mfunc, 
                blkmode=send_time_blks)
            break
        elif not send_time_blks or uvi.vars.changed():
            _process_buf(uvi, uvo, pbuf, dbuf, vbuf, mfunc, 
                blkmode=send_time_blks)
            pbuf, dbuf = [], []
            vbuf = {}
            for b in bufvars: vbuf[b] = []
        for b in bufvars: vbuf[b].append(uvi.vars[b])
        pbuf.append(p); dbuf.append(d)

def _process_buf(uvi, uvo, pbuf, dbuf, vbuf, mfunc, blkmode=False):
    """Write a buffer of input data to the output uv file.
       'uvi' is the input file.
       'uvo' is the input file.
       'pbuf' is a list of preambles.
       'dbuf' is a list of data.
       'vbuf' is a dictionary of any specified variables & their values.
       'mfunc' is the mapping function.
       'blkmode' specifies if more than 1 spectrum is being sent.""" 
    if len(pbuf) == 0: return
    elif not blkmode:
        if mfunc is not None:
            pbuf, dbuf = mfunc(uvi, pbuf[0], dbuf[0], vbuf)
            pbuf, dbuf = [pbuf], [dbuf]
    else:
        if mfunc is not None: pbuf, dbuf = mfunc(uvi, pbuf, dbuf, vbuf)
    for p, d in zip(pbuf, dbuf):
        miruv.uvcopyvr_c(uvi.handle, uvo.handle)
        uvo.write_data(p, d)

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
    uv = UV(sys.argv[1], status='new')
    uv.items['history'] = HIST_STR
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
    if uv.items['history'] != HIST_STR:
        raise ValueError('UVItem read/write failed.')
    for t in range(100):
      for i in range(8):
        for j in range(8):
            bl = float((i+1)<<8 + j+1)
            preamble = numpy.array((1., 1., 1., time, bl), dtype=numpy.double)
            p, d = uv.read_data()
            if numpy.any(d.filled() != a.filled()) or numpy.any(p != preamble):
                raise ValueError('Data was not written correctly.')
            if uv.vars['pol'] != i:
                raise ValueError('UVVar read/write failed.')
      time += 10
    del(uv)

    print 'Passed test.'
