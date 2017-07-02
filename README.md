# AIPY (Astronomical Interferometry in PYthon)

[![Build Status](https://travis-ci.org/HERA-Team/aipy.svg?branch=master)](https://travis-ci.org/HERA-Team/aipy)
[![Coverage Status](https://coveralls.io/repos/github/HERA-Team/aipy/badge.svg?branch=master)](https://coveralls.io/github/HERA-Team/aipy?branch=master)

# Description
This package collects together tools for radio astronomical interferometry.
In addition to pure-python phasing, calibration, imaging, and
deconvolution code, this package includes interfaces to MIRIAD (a Fortran
interferometry package) and HEALPix (a package for representing spherical 
data sets), and some math/fitting routines from SciPy. 

# On the Web
There are further instructions, documentation, and a FAQ at:
http://github.com/HERA-Team/aipy

# Installation
AIPY depends on the following python packages:
1. numpy >= 1.2
2. pyephem >= 3.7.3
3. astropy, or pyfits >= 1.3
4. matplotlib >= 0.98.3 (optional)
5. matplotlib-basemap >= 0.99 (optional)


## Installing Dependencies
OPTION 1 (safest): Manually install the dependencies.
OPTION 2 (experimental): Open up the AIPY download, and with network connectivity and root access, type:
```
$ install_required.sh
```
and then (if you want matplotlib/basemap):
```
$ install_recommended.sh
```

## Install as User
```
$ python setup.py install
```

## Install as Root
```
$ sudo python setup.py install
```

## Command-line Scripts
Unless you installed as user, these will be in /usr/bin.  For more info
use the "-h" option with any of these commands:

# Documentation
If you want to build html documentation, you'll first need sphinx on 
your system:
```
$ easy_install -U sphinx
```
Then you should be able to cd into the doc directory and run:
```
$ make html
```
The results will appear in doc/build/html.  There are plenty of other 
build options, too.  Many thanks to Jayce Dowell for his work on this.

Enjoy,
Aaron Parsons

-----------------------------------------------------------------------------

# Package Info for Developers
The subpackage "optimize" was copied in from scipy-0.6.0, and then all
code that depended on non-pure-python modules was removed.  If these ever 
need to be updated, download scipy source and copy scipy/scipy/optimize 
into aipy, and then remove any code deemed unnecessary.  Unfortunately, 
then you may need to crawl through the code and replace all "scipy" 
references with "aipy".

The subpackage "_cephes" was copied in from scipy-0.6.0/special, and then
all but the cephes and c_misc code was removed to avoid needing a Fortran
compiler.  _cephesmodule.c needed substantial editing to remove external
dependencies.

The miriad source code (aipy/miriad/mirsrc) was included from 
MIRIAD-4.0.5.  To update, download miriad source and copy $MIR/src/subs/* 
and $MIR/src/inc/* into aipy/miriad/mirsrc.  Not all files are used, but 
include them all anyway.

Healpix source code (aipy/healpix/cxx) was included from Healpix-2.01.
To update, download healpix source and copy src/cxx into aipy/healpix.
