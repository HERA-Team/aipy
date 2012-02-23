General Utilities
=================
lst
  :Description: Output LST, Julian date, and Sun position at the provided location

  :Usage: lst jd1 [jd2 [...]]

  :Options: None

compress_uv.py
  :Description: Tarball and compress (using bz2) Miriad UV files.

  :Usage: compress_uv.py [options] file.uv

  :Options: -h, --help    show this help message and exit

            -d, --delete  Delete a uv file after compressing it

            -x, --expand  Inflate tar.bz2 files

combine_freqs.py 
  :Description: A script for reducing the number of channels in a UV data set by coherently adding adjacent channels together.

  :Usage: combine_freqs.py [options] file.uv

  :Options: -h, --help               Show this help message and exit

            -n NCHAN, --nchan=NCHAN  Reduce the number of channels in a spectrum to this number.

            -c, --careful_flag       Flag resultant bin if any component bins are flagged (otherwise, flags only when there are more flagged bins than unflagged bins).

            -d, --dont_flag          Only flag resultant bin if every component bin is flagged (otherwise, uses any data available).

            -u, --unify              Output to a single UV file.

Calibration
===========
apply_bp.py
  :Description: Apply the bandpass function in a UV file to the raw data, and then write to a new UV file which will not have a bandpass file.  Has the (recommended) option of linearizing for quantization gain effects.  For now, these are applied only to the auto-correlations.  Cross-correlation quantization gain is less sensitive to level, so linearization will come later.

  :Usage: apply_bp.py [options] file.uv

  :Options: -h, --help            show this help message and exit

            -l LINEARIZATION, --linearization=LINEARIZATION  Apply the specified quantization linearization function to raw correlator values before applying bandpass.  Options are null, digi, full, and comb.  Default is comb

            -s SCALE, --scale=SCALE  An additional numerical scaling to apply to the data.  Default: 12250000.

flux_cal.py
  :Description: A script for dividing out the passband, primary beam, and/or source spectrum scaling.  When dividing by a primary beam or source spectrum, it is recommended a single source have been isolated in the data set.

  :Usage: flux_cal.py [options] file.uv

  :Options: -h, --help         show this help message and exit

            -C CAL, --cal=CAL  Use specified <cal>.py for calibration information.

            -s SRC, --src=SRC  Phase centers/source catalog entries to use.  Options are "all", "<src_name1>,...", or "<ra XX[:XX:xx]>_<dec XX[:XX:xx]>".

            --cat=CAT          A comma-delimited list of catalogs from which sources are to be drawn.  Default is "helm,misc".  Other available catalogs are listed under aipy._src.  Some catalogs may require a separate data file to be downloaded and installed.

            -b, --beam         Normalize by the primary beam response in the direction of the specified source.

            -p, --passband     Normalize by the passband response. 

            -f, --srcflux      Normalize by the spectrum of the specified source.

Modeling
========
mdlvis.py
  :Description: Models visibilities for various catalog sources and creates a new Miriad UV file containing either the simulated data, or the residual when the model is removed from measured data.

  :Usage: mdlvis.py [options] file.uv

  :Options: -h, --help            show this help message and exit

            -C CAL, --cal=CAL     Use specified <cal>.py for calibration information.

            -s SRC, --src=SRC     Phase centers/source catalog entries to use.  Options are "all", "<src_name1>,...", or "<ra XX[:XX:xx]>_<dec XX[:XX:xx]>".

            --cat=CAT             A comma-delimited list of catalogs from which sources are to be drawn.  Default is "helm,misc".  Other available catalogs are listed under aipy._src.  Some catalogs may require a separate data file to be downloaded and installed.

            -m MODE, --mode=MODE  Operation mode.  Can be "sim" (output simulated data), "sub" (subtract from input data), or "add" (add to input data).  Default is "sim"

            -f, --flag            If outputting a simulated data set, mimic the data flagging of the original dataset.

            -n NOISELEV, --noiselev=NOISELEV  RMS amplitude of noise (Jy) added to each UV sample of simulation.

            --nchan=NCHAN         Number of channels in simulated data if no input data to mimic.  Default is 256

            --sfreq=SFREQ         Start frequency (GHz) in simulated data if no input data to mimic.  Default is 0.075

            --sdf=SDF             Channel spacing (GHz) in simulated data if no input data to mimic.  Default is .150/256

            --inttime=INTTIME     Integration time (s) in simulated data if no input data to mimic.  Default is 10

            --startjd=STARTJD     Julian Date to start observation if no input data to mimic.  Default is 2454600

            --endjd=ENDJD         Julian Date to end observation if no input data to mimic.  Default is 2454601

            --pol=POL             Polarizations to simulate (xx,yy,xy,yx) if starting file from scratch.

filter_src.py
  :Description: A script for filtering using a delay/delay-rate transform.  If a source is specified, will remove/extract that source.  If none is specified, will filter/extract in absolute terms.

  :Usage: filter_src.py [options] file.uv

  :Options: -h, --help         show this help message and exit

            -C CAL, --cal=CAL  Use specified <cal>.py for calibration information.

            -s SRC, --src=SRC  Phase centers/source catalog entries to use.  Options are "all", "<src_name1>,...", or "<ra XX[:XX:xx]>_<dec XX[:XX:xx]>".

            --cat=CAT          A comma-delimited list of catalogs from which sources are to be drawn.  Default is "helm,misc".  Other available catalogs are listed under aipy._src.  Some catalogs may require a separate data file to be downloaded and installed.

            -r DRW, --drw=DRW  The number of delay-rate bins to null.  Default is -1 = no fringe filtering.

            -d DW, --dw=DW     The number of delay bins to null. If -1, uses baseline lengths to generate a sky-pass filter.

            -p, --passband     Divide by the model passband before transforming.

            -e, --extract      Extract the source instead of removing it.

            --clean=CLEAN      Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).

fitmdl.py
  :Description: A script for fitting parameters of a measurement equation given starting parameters in a cal file and a list of sources.  The fitter used here is a steepest-decent filter and does not make use of priors.

  :Usage: fitmdl.py [options] file.uv

  :Options: -h, --help            show this help message and exit

            -a ANT, --ant=ANT     Select ants/baselines to include.  Examples: all (all baselines) auto (of active baselines, only i=j) cross (only i!=j) 0,1,2 (any baseline involving listed ants) 0_2,0_3 (only listed baselines) "(0,1)_(2,3)" (same as 0_2,0_3,1_2,1_3. Quotes help bash deal with parentheses) "(-0,1)_(2,-3)" (exclude 0_2,0_3,1_3 include 1_2).  Default is "cross".

            -p POL, --pol=POL     Choose polarization (xx, yy, xy, yx) to include.

            -c CHAN, --chan=CHAN  Select channels (after any delay/delay-rate transforms) to include.  Examples: all (all channels), 0_10 (channels from 0 to 10, including 0 and 10) 0_10_2 (channels from 0 to 10, counting by 2), 0,10,20_30 (mix of individual channels and ranges).  Default is "all".

            -C CAL, --cal=CAL     Use specified <cal>.py for calibration information.

            -s SRC, --src=SRC     Phase centers/source catalog entries to use.  Options are "all", "<src_name1>,...", or "<ra XX[:XX:xx]>_<dec XX[:XX:xx]>".

            --cat=CAT             A comma-delimited list of catalogs from which sources are to be drawn.  Default is "helm,misc".  Other available catalogs are listed under aipy._src.  Some catalogs may require a separate data file to be downloaded and installed.

            -P PRMS, --prms=PRMS  Parameters (for fitting, usually), can be specified as can be: "obj=prm", "obj=prm/val", "obj=prm/val/sig", "(obj1/obj2)=prm/(val1/val2)/sig", "obj=(prm1/prm2)/val/(sig1/sig2)", comma separated versions of the above, and so on.

            -x DECIMATE, --decimate=DECIMATE  Use only every Nth integration.  Default is 1. 

            --dphs=DECPHS         Offset to use when decimating (i.e. start counting integrations at this number for the purpose of decimation).  Default is 0.

            -S SHPRMS, --shared_prms=SHPRMS  Parameter listing with the same syntax as "-P/--prms" except that all objects listed for a parameter will share an instance of that parameter.

            --snap                Snapshot mode.  Fits parameters separately for each integration.

            -q, --quiet           Be less verbose.

            --maxiter=MAXITER     Maximum # of iterations to run.  Default is infinite.

            --xtol=XTOL           Fractional change sought in it parameters before convergence.  Default 1e-10.

            --ftol=FTOL           Fractional tolerance sought in score before convergence.  Default 1e-10.

            --remem               Remember values from last fit when fitting in snapshot mode.

            --baseport=BASEPORT   Base port # to use for tx/rx.  Each daemon adds it's daemon id to this to determine the actual port used for TCP transactions.

            --daemon=DAEMON       Operate in daemon mode, opening a TCP Server to handle requests on the specified increment to the base port.

            --master=MASTER       Operate in master mode, employing daemon-mode servers to do the work and collecting the results.  Should be a comma delimited list of host:daemonid pairs to contact.  Daemon ID will be added to baseport to determine actual port used for TCP transactions.

            --sim_autos           Use auto-correlations in fitting.  Default is to use only cross-correlations.

Imaging
=======
mk_imag.py
  :Description: This is a general-purpose script for making images from MIRIAD UV files.  Data (optionally selected for baseline, channel) are read from the file, phased to a provided position, normalized for passband/primary beam effects, gridded to a UV matrix, and imaged.

  :Usage: mk_img.py [options] file.uv

  :Options: -h, --help            show this help message and exit

            -a ANT, --ant=ANT     Select ants/baselines to include.  Examples: all (all baselines) auto (of active baselines, only i=j) cross (only i!=j) 0,1,2 (any baseline involving listed ants) 0_2,0_3 (only listed baselines) "(0,1)_(2,3)" (same as 0_2,0_3,1_2,1_3.  Quotes help bash deal with parentheses) "(-0,1)_(2,-3)" (exclude 0_2,0_3,1_3 include 1_2).  Default is "cross".

            -p POL, --pol=POL     Choose polarization (xx, yy, xy, yx) to include.

            -c CHAN, --chan=CHAN  Select channels (after any delay/delay-rate transforms) to include.  Examples: all (all channels), 0_10 (channels from 0 to 10, including 0 and 10) 0_10_2 (channels from 0 to 10, counting by 2), 0,10,20_30 (mix of individual channels and ranges).  Default is "all".

            -C CAL, --cal=CAL     Use specified <cal>.py for calibration information.

            -s SRC, --src=SRC     Phase centers/source catalog entries to use.  Options are "all", "<src_name1>,...", or "<ra XX[:XX:xx]>_<dec XX[:XX:xx]>".

            --cat=CAT             A comma-delimited list of catalogs from which sources are to be drawn.  Default is "helm,misc".  Other available catalogs are listed under aipy._src.  Some catalogs may require a separate data file to be downloaded and installed.

            -x DECIMATE, --decimate=DECIMATE  Use only every Nth integration.  Default is 1. 

            --dphs=DECPHS         Offset to use when decimating (i.e. start counting integrations at this number for the purpose of decimation).  Default is 0.

            -o OUTPUT, --output=OUTPUT  Comma delimited list of data to generate FITS files for.  Can be: dim (dirty image), dbm (dirty beam), uvs (uv sampling), or bms (beam sampling).  Default is dim,dbm.

            --list_facets         List the coordinates of all the pointings that will be used.

            --facets=FACETS       If no src is provided, facet the sphere into this many pointings for making a map.  Default 200.

            --snap=SNAP           Number of integrations to use in "snapshot" images.  Default is to not do snapshoting (i.e. all integrations go into one image).

            --cnt=CNT             Start counting output images from this number.  Default 0.

            --fmt=FMT             A format string for counting successive images written to files.  Default is im%04d (i.e. im0001).

            --skip_phs            Do not phase visibilities before gridding.

            --skip_amp            Do not use amplitude information to normalize visibilities.

            --skip_bm             Do not weight visibilities by the strength of the primary beam.

            --skip=SKIP           Skip this many pointings before starting.  Useful in conjungtion with --cnt for resuming.

            --size=SIZE           Size of maximum UV baseline.

            --res=RES             Resolution of UV matrix.

            --no_w                Don't use W projection.

            --altmin=ALTMIN       Minimum allowed altitude for pointing, in degrees.  When the phase center is lower than this altitude, data is omitted.  Default is 0.

            --minuv=MINUV         Minimum distance from the origin in the UV plane (in wavelengths) for a baseline to be included.  Default is 0.

            --buf_thresh=BUF_THRESH  Maximum amount of data to buffer before gridding.  Excessive gridding takes performance hit, but if buffer exceeds memory available... ouch

cl_img.py
  :Description: This is a general-purpose script for deconvolving dirty images by a corresponding PSF to produce a clean image.

  :Usage: cl_img.py [options] file.dim.fits file.dbm.fits

  :Options: -h, --help            show this help message and exit

            -d DECONV, --deconv=DECONV  Attempt to deconvolve the dirty image by the dirty beam using the specified deconvolver (none,mem,lsq,cln,ann).

            -o OUTPUT, --output=OUTPUT  Comma delimited list of data to generate FITS files for.  Can be: cim (clean image), rim (residual image), or bim (best image = clean + residuals). Default is bim.

            --var=VAR             Starting guess for variance in maximum entropy fit (defaults to variance of dirty image.

            --tol=TOL             Tolerance for successful deconvolution.  For annealing, interpreted as cooling speed.

            --div                 Allow clean to diverge (i.e. allow residual score to increase)

            -r REWGT, --rewgt=REWGT  Reweighting to apply to dim/dbm data before cleaning. Options are: natural, uniform(LEVEL), or radial, where LEVEL is the fractional cutoff for using uniform weighting (recommended range .01 to .1).  Default is natural.

            --maxiter=MAXITER     Number of allowable iterations per deconvolve attempt.
