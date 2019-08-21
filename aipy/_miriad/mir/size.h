c=======================================================================
c	[size.h]
c
c  Include file used by calibration software
c
c  UL, RI, ULRI -- parameters which stand for LSB/USB and REAL/IMAG stuff
c
	integer UL, RI, ULRI
	parameter( UL = 2 )
	parameter( RI = 2 )
	parameter( ULRI = 4 )
c
c  SIZE -- lengths of reals/integers in number of bytes (IEEE/2s compl)
c
	integer SIZE
	parameter( SIZE = 4 )
c=======================================================================

