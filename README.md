# AIPY (Astronomical Interferometry in PYthon)

[![Build Status](https://travis-ci.org/HERA-Team/aipy.svg?branch=master)](https://travis-ci.org/HERA-Team/aipy)
[![Coverage Status](https://coveralls.io/repos/github/HERA-Team/aipy/badge.svg?branch=master)](https://coveralls.io/github/HERA-Team/aipy?branch=master)

## Description

This package collects together tools for radio astronomical interferometry.
In addition to pure-python phasing, calibration, imaging, and
deconvolution code, this package includes interfaces to MIRIAD (a Fortran
interferometry package) and HEALPix (a package for representing spherical
data sets), and some math/fitting routines from SciPy.

Instructions, documentation, and a FAQ may be found at
[the aipy GitHub page](http://github.com/HERA-Team/aipy).

## Installation

We strongly recommend installing AIPY with the [conda](https://conda.io/docs/)
packaging tool using the public [conda-forge](https://conda-forge.org/)
package repository:

```
$ conda install -c conda-forge aipy
```

## Python 3 and AIPY 3

The 2.x release series of AIPY supports only Python 2. However, there is now a
**prototype** version of **AIPY 3.0** which adds support for **Python 3**
(while still supporting Python 2). You can install it using the following
[pip](https://pip.pypa.io/en/stable/) command:

```
pip install git+https://github.com/HERA-Team/aipy.git@v3
```

There is not presently a Conda package of this version. Development of the 3.x
series happens on the [v3](https://github.com/HERA-Team/aipy/commits/v3)
branch of this repository.


## Documentation

If you want to build HTML documentation, you'll need to have
[Sphinx](http://www.sphinx-doc.org/) installed. Then, change to the `doc/`
directory and run:

```
$ make html
```

The results will appear in `doc/build/html`.  There are plenty of other
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
