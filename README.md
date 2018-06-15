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

The miriad source code (aipy/miriad/mirsrc) was included from
MIRIAD-4.0.5.  To update, download miriad source and copy $MIR/src/subs/*
and $MIR/src/inc/* into aipy/miriad/mirsrc.  Not all files are used, but
include them all anyway.

Healpix source code (aipy/healpix/cxx) was included from Healpix-2.01.
To update, download healpix source and copy src/cxx into aipy/healpix.
