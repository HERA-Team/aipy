Requirements
============
The requirements for AIPY are:
 * numpy >= 1.2
 * pyephem >= 3.7.2.3
 * astropy >= 1.0, or pyfits >= 1.1
 * matplotlib >= 0.98 [1]_
 * matplotlib-basemap >= 0.99 [1]_

.. [1] Required for some of the included scripts

Installing
==========
Install AIPY by running::

	python setup.py install [--prefix=<prefix>|--user]

If the '--prefix' option is not provided, then the installation tree root 
directory is the same as for the python interpreter used to run setup.py.  
For instance, if the python interpreter is in '/usr/local/bin/python', 
then <prefix> will be set to '/usr/local'.  Otherwise, the explicit <prefix> 
value is taken from the command line option.  The package will install 
files in the following locations:
 * <prefix>/bin
 * <prefix>/lib/python2.6/site-packages
 * <prefix>/share/doc
 * <prefix>/share/install
If an alternate <prefix> value is provided, you should set the PATH 
environment to include directory '<prefix>/bin' and the PYTHONPATH 
environment to include directory '<prefix>/lib/python2.6/site-packages'.

If the '--user' option is provided, then then installation tree root 
directory will be in the current user's home directory. 

Obtaining Catalogs
==================
Most of the catalogs that AIPY provides an interface to, e.g. the Third Cambridge
Catalog, needed to be downloaded from `VizieR <http://vizier.u-strasbg.fr/cgi-bin/VizieR>`_ as
tab-separated files and placed in the AIPY `_src` dirctory before they can be
used.  See the documentation for the :mod:`aipy.src` module for details.


