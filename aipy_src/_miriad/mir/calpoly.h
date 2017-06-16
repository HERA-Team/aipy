c=======================================================================
c	[calpoly.h]
c
c	last update: 26-apr-90	new poly's
c		     18-nov-90  descriptive comments
c		      7-jul-91  new MAXBASHC parameter
c*calpoly_h -- header file for calibration data
c:calibration, include
c&pjt
c+
c
c  order -- current polynomial order for a slot
c
      INTEGER order(MAXSLOT)

c
c  pcount -- number of polynomials in a (slot,baseline) pair
c
      INTEGER pcount(MAXSLOT, MAXBASHC)

c
c  tvalid - validity of time interval
c		(min & max time-range; breakpoints; slotcode; baseline)
c
      REAL tvalid(2,MAXBREAK+1,MAXSLOT,MAXBASHC)
c
c  poly -- actual polynomial coefficients 
c		(order; breakpoints, slotcode, baseline)
c
      REAL poly( 0:MAXORDER, MAXBREAK+1, MAXSLOT, MAXBASHC )

c
c  CALPOLY -- the common block
c
      COMMON /calpoly/ pcount, order, tvalid, poly
c--
c=======================================================================
