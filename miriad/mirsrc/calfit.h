c=======================================================================
c	[calfit.h]
c*calfit_h -- calibration include file
c&pjt
c:calibration, include
c+
c  lsq()        matrix used for inversion
c  lsqb()       vector, used as right hand side AND solution
c  xsums()      temporary storage use to fill lsq()
c  ysums()      temporary storage use to fill lsq()
c
c
c  MAXLSQ -- the size of arrays and matrices in the least squares fit
c
      INTEGER MAXLSQ
      PARAMETER(MAXLSQ=(MAXANTHC*(MAXORDER+1)))
c
c  lsq, lsqb -- least squares matrix and rhs
c
      DOUBLE PRECISION lsq( MAXLSQ, MAXLSQ ) 
      DOUBLE PRECISION lsqb( MAXLSQ )
c
c  xsums, ysums -- sums-of-squares
c
      DOUBLE PRECISION xsums( 0 : 2 * MAXORDER )
      DOUBLE PRECISION ysums( 0 : MAXORDER )
c
      COMMON /lsq/ xsums, ysums, lsq, lsqb
c--
c=======================================================================

