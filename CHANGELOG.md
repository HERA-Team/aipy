# Next version (not yet released)

- Remove `ez_setup` installation infrastructure; it has been deprecated for a
  long time.
- Allow installing from empty environment
- Use github actions instead of travis
- Use python3 print()
- Use pytest instead of unittest
- Use explicit numpy types

# 3.0.1 (2019 Aug 21)

- Remove Astropy version cap; it was intended to promote simultaneous Python
  2/3 compatibility, but if someone wants to use the newer version, that's
  fine too.
- Fix a bug in the Healpix wrapper
  ([#52](https://github.com/HERA-Team/aipy/pull/52)).
- Fix some lingering Python 2/3 compatibility issues
  ([#51](https://github.com/HERA-Team/aipy/pull/51)).
- Honor `pos_def` in the 1D real CLEAN
  ([#44](https://github.com/HERA-Team/aipy/issues/44))
- Install the Helmboldt catalog data files, which were mistakenly not being
  included in installations
  ([#45](https://github.com/HERA-Team/aipy/issues/45))
- Rename the package directory in the source tree from `aipy_src` to `aipy`.

# 3.0.0rc2 (2018 Aug 27)

- Fix generation of the package description in `setup.py`
  ([#46](https://github.com/HERA-Team/aipy/issues/46),
  [#47](https://github.com/HERA-Team/aipy/issues/47))
- Fix division that should have become integer division
  ([#48](https://github.com/HERA-Team/aipy/issues/48))

# 3.0.0rc1 (2018 Aug 20)

- Make compatible with Python 3; remove scipy and healpix borrows.

# 2.x series

- No changelog was maintained. Consult the Git history.

1.0.0
* Added some performance enhancements
* Miriad.UV files can now store data as both shorts and as float32 (avoiding dynamic range problems)
* Various bugfixes
* Because of stability of 0.X.X series, we've promoted this to a non-beta release!

0.9.0
* Made some major naming changes in preparation for an AIPY publication.  These include
  renaming aipy.ant to aipy.phs, aipy.sim to aipy.amp, aipy.loc to aipy.cal.
* Renamed BeamFlat to just plain Beam
* Changed the -l options (for "loc files") in scripting interfaces to -C (for "cal" files)
* Added support in scripting interface channel selection for constructions of the form
  0_10_2, meaning channels 0 to 10 (including endpoints), counting by 2.
* Refactored code so source spectra are computed internally to RadioBodys, and are
  accessed through the get_jys() method.  Added an update_jys() method for computing
  a source spectrum for specified frequencies without needing an AntennaArray object.
* Added a get_afreqs() method to AntennaArrays to elucidate the opaque incantation;
  "AntennaArray.ants[0].beam.afreqs."
* Added a system support for adding source catalogs through src.py calling catalog
  interfaces stored in the _srcs directory.  Made the Helmboldt catalog the default
  catalog that AIPY comes with (since it's pretty small and geared toward low frequencies).
* Added unit tests (in "test" directory) for the most-used modules, and will continue to
  add them as other modules are revisited.  There is also a script called
  "script_battery.sh" that does basic exercising of most AIPY scripts.
* In Antennas, grouped delay/offset into an array called phsoff which is a phase
  polynomial in frequency.
* Couple of bugfixes on how objects are updated (now made explicit with "update()"
  methods) and fixed an array shape bug in Antenna.refract()
* Fixed a bug in filter_src.py where data were zeroed when a source is below the horizon.
* Transitioned all nomenclature from "Fringe-Rate" to "Delay-Rate" as per
  Parsons & Backer (2009) in AJ.
* Changed how fit parameters are input in scripting module
* Bugfix for when resolution effects are taken into account in
  AntennaArray.gen_phs.
* Bugfix in xrfi.py for flagging by integration
* Improved filtering in filter_src.py to use padding and optionally passband
  information.
* Added script srclist.py
* Added masking by wieght to plot_map.py

0.8.6
* Fixed a data-type error for Healpix maps that aren't doubles
* Added max/drng/cmap parameters to standard command-line options
* Started adding unit tests
* Changed noise option in mdlvis.py to add noise in Jy (not data) units.
* Changed behavior of reducing NSIDE in Healpix maps to sum instead of average.
* Bugfix on extracting source below the horizon in filter_src.py
* Work aound problem with numpy MA arrays not having a mask if entire image
  is valid.
* Added new option (0,1,2)_3 for specifying antennas on the command line.

0.8.5
* Added support to mdlvis.py for simulating data without an input UV file.
* Added "none" mode to xrfi.py to flag only those channels specified manually.
* Added line to alm_wrapper.cpp to support compilation on python < 2.5.

0.8.4
* Changed modmap.py to work with sources of nonzero angular size.  Also added
  support for changing the data type (double,float) of a map.
* Added parameter to mk_img.py to control minimum altitude to allow for
  phase center.  Also added parameter for resuming interrupted faceting.
* Fixed xrfi.py default parameter "intthresh" to avoid excessive flagging
  of integrations.
* Made fitmdl.py skip fitting when there are no data.  Added "snap" mode to
  fit parameters every integration.
* Added text color option to plot_map.py
* Added ionospheric refraction code using 1/freq**2 scaling.  This involved
  some significant refactoring of code, and some changes in fitmdl.py and
  mdlvis.py.
* Changed fitmdl.py to cache UV data, resulting in less time reading files.
  This significantly sped up code, although may become problematic with
  fits involving lots of data.
* Removed some per-baseline computations in ant.py that occured every
  compute, even if unnecessary.  Baseline projections are now computed as
  needed.  Sped up simulations involving few baselines.
* Bugfix to miriad "select" command involving "time" selection: double
  arguments were incorrectly cast as floats.
* Made mk_map.py a little more efficient in calculating map weights.

0.8.3
* Changed how loc-specific source parameters are implemented (now with
  get_catalog() instead of src_prms).  Changed scripts to use this format,
  and scripting.parse_srcs() to return values that can feed get_catalog().

0.8.2
* Built in decimation to UV files in miriad_wrap.cpp.
* Added decimation selection with select() to UV files in miriad.py.
* Bugfix on various scripts to use new clean command.

0.8.1
* Split out cleaning functionality from mk_img.py into cl_img.py.
* Added snapshot capability to mk_img.py
* Bugfix for plot_img.py grid being upside-down.
* Reimplemented clean in C for speed.

0.8.0
* Added resolving support for elliptical sources in resolve_src() in ant.py.
* Reworked how get_baseline() and gen_uvw() functioned to facilitate elliptical
  sources.  This change required updating scripts that use gen_uvw().
* Allowed specifying coordinate epochs for RadioFixedBody's
* Changed how plot_img.py, mk_map.py work.  There is now a general script
  mk_img.py which geberates FITS files (introduced in 0.7.1) for each specified
  source, or for equally-spaced phase centers around the sphere.  This script
  is responsible for all image generation.  mk_map.py now takes these FITS
  images as input and uses them as facet images to make a spherical Healpix
  map.  plot_img.py is now a visualization tool (a la plot_uv.py and
  plot_map.py), and overlays appropriate ra/dec grids on the FITS images.
* Flipped u/v order (and sign of v) in the Img class in img.py to be more
  naturally compatible with FITS image conventions.  Required a change to
  Img.gen_LM() as well.

0.7.2
* Added sorting to plot_uv.py to make plots in sensible baseline order.
* Edited script documentation for plot_map.py
* Changed "-a" in scripting.py so that "-a 0,1,2" includes only baselines that
  involve these antennas without autos (i.e. 0_1,0_2,1_2), rather than all
  baselines that involve any of the antennas. Listing 1 antenna will
  include all baselines involving that antenna.
* Fixed combine_freqs.py to work with new numpy masked arrays.
* Made fitmdl.py have a --sim_auto parameter that rejects auto-correlations
  by default.  This makes fitmdl.py work more nicely with '-a' change above.

0.7.1
* Changed "flat" to "flatten()" in plot_img.py, as per numpy-1.2.
* Made plot_map.py able to function (only for cyl projection) without basemap.
* Added FITS image I/O convenience functions to img.py.  Defined standard
  FITS format based on VLA images.

0.7.0
* Moved functions from xrfi.py script into rfi module
* Bugfix in array shape in BeamFlat response() function
* Ported to numpy>=1.2.  Involved:
  - Update flag_by_int in rfi module for new numpy-1.2.ma behavior
  - Update miriad_wrap.cpp, healpix_wrap.cpp, alm_wrap.cpp for new
    numpy-1.2 API (PyArray_SimpleNew instead of PyArray_FromDims,
    npy_intp instead of int for array dimensions)

0.6.5
* Bugfix to source strength polynomials as a function of LST
* Corrected all 3C source postions in src.py using NED.
* Ported to numpy>=1.1, matplotlib>=0.98, basemap>=0.99.  Involved:
    - use of numpy.histogram with new=True in xrfi.py
    - importing Basemap from mpl_toolkits.basemap in plot_map.py
* Split dependency installation into two scripts: install_required.sh
  which installs the minimum needed to be able to install aipy and run
  "import aipy", and install_recommended.sh which installs matplotlib/basemap
  and is needed for most of the scripts that come with AIPY.

0.6.4
* Attempting to support easy_install/setuptools functionality while maintaining
  the simple building that already works.
* Feature added to plot_map.py to give (ra,dec),flux of location on a map
  when right-clicked.
* Changed how plot_map.py plots to use imshow instead of contourf.

0.6.3
* Added _cephes module (from scipy, cutting out fortran dependencies) which
  implements various scientific functions.  Needed j1 = Bessel function.
* Moved source-resolution estimator (using symmetric disc approx) from phs2src
  into separate AntennaArray.resolve_src().  Changed from incorrect sinc
  function to correct j1 (first order Bessel function).
