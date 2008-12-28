c=======================================================================
c - mirconst.h  Include file for various fundamental physical constants.
c
c  History:
c    jm  18dec90  Original code.  Constants taken from the paper
c                 "The Fundamental Physical Constants" by E. Richard
c                 Cohen and Barry N. Taylor (PHYICS TODAY, August 1989).
c ----------------------------------------------------------------------
c  Pi.
      real PI, TWOPI
      double precision DPI, DTWOPI
      parameter (PI = 3.14159265358979323846)
      parameter (DPI = 3.14159265358979323846)
      parameter (TWOPI = 2 * PI)
      parameter (DTWOPI = 2 * DPI)
c ----------------------------------------------------------------------
c  Speed of light (meters/second).
      real CMKS
      double precision DCMKS
      parameter (CMKS = 299792458.0)
      parameter (DCMKS = 299792458.0)
c ----------------------------------------------------------------------
c  Boltzmann constant (Joules/Kelvin).
      real KMKS
      double precision DKMKS
      parameter (KMKS = 1.380658E-23)
      parameter (DKMKS = 1.380658D-23)
c ----------------------------------------------------------------------
c  Planck constant (Joules-second).
      real HMKS
      double precision DHMKS
      parameter (HMKS = 6.6260755E-34)
      parameter (DHMKS = 6.6260755D-34)
c ----------------------------------------------------------------------
c  Planck constant divided by Boltzmann constant (Kelvin/GHz).
      real HOVERK
      double precision DHOVERK
      parameter (HOVERK = 0.04799216)
      parameter (DHOVERK = 0.04799216)
c=======================================================================
