# AIPY (Astronomical Interferometry in PYthon)

[![Build Status](https://github.com/HERA-Team/aipy/workflows/Run%20Tests/badge.svg)](https://github.com/HERA-Team/aipy/actions)
[![Coverage Status](https://codecov.io/gh/HERA-Team/aipy/branch/master/graph/badge.svg?token=Vrr4XFcE8p)](https://codecov.io/gh/HERA-Team/aipy)

## Description

This package collects together tools for radio astronomical interferometry.
In addition to pure-python phasing, calibration, imaging, and
deconvolution code, this package includes interfaces to MIRIAD (a Fortran
interferometry package) and HEALPix (a package for representing spherical
data sets).

Instructions, documentation, and a FAQ may be found at
[the aipy GitHub page](http://github.com/HERA-Team/aipy).

## Installation

We strongly recommend installing AIPY with the [conda](https://conda.io/docs/)
packaging tool using the public [conda-forge](https://conda-forge.org/)
package repository:

```
$ conda install -c conda-forge aipy
```

As of the 3.0.x version series, AIPY supports both Python 3 and Python 2.

You can also install with `pip` if you wish:

```
$ pip install aipy
```

To install the source code in development mode, use:

```
$ pip install -e .
$ python setup.py build_ext --inplace
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

## Package Info for Developers

The miriad source code (`aipy/miriad/mirsrc`) was included from MIRIAD 4.0.5.
To update, download a MIRIAD distribution and copy `$MIR/src/subs/*` and
`$MIR/src/inc/*` into `aipy/miriad/mirsrc`.

Healpix source code (`aipy/healpix/cxx`) was included from Healpix 2.01. To
update, download a HEALPix distribution and copy `src/cxx` into
`aipy/healpix`.

## Making Releases (for Maintainers)

To make a release of `aipy` (both on Github and PyPI), head to the most current 
[Draft Release](https://github.com/HERA-Team/aipy/releases) and note the *suggested*
release version. Contact the maintainers with your intention to make a release either 
to that version (or, if appropriate, to a different version), and publish the release
via the Github UI. All done!