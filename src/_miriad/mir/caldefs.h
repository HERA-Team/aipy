c=======================================================================
c [caldefs.h]
c General include file that every calibration routine needs. It does not
c define (space for) data, but defines the constants that should be used
c in defining array-sizes. 
c	18-nov-90 	added PI and TWOPI			PJT
c	 6-jul-91	V6.3 new MAXANTHC and MAXBASHC          pjt
c	26-sep-91       increased MAXUVPNT and MAXBREAK
c	 2-dec-91	V6.4 changed calbflux's dimensionality	PJT
c	19-apr-93	MAXSLOT from 20 to 36 for screwed mode-4  PJT
c			added MAXWIN from 8 (hardcoded) to 16
c	15-apr-94	MAXWIN is now in maxdim.h		PJT/mchw
c       19-dec-95       increase MAXANTHC from 6 to 9 !!!       pjt
c        9-Feb-99       increase MAXANTHC from 9 to 10 !       mchw
c=======================================================================
      INCLUDE 'maxdim.h'
c=======================================================================
      INCLUDE 'size.h'
c=======================================================================
c
c  SVERSION -- current version of cal format (consistency checks)
c              Version number is 10*MAJOR+MINOR
c              Datasets with a MAJOR difference, cannot be read
c              anymore, MINOR differences are normally compatible
c
      INTEGER SVERSION
      PARAMETER ( SVERSION = 64)
c
c MAXANTHC -- special Hat Creek MAXANT, to keep it small for us....
c MAXBASHC -- derived from MAXANTHC = max. number of baselines
c
      INTEGER MAXANTHC, MAXBASHC
      PARAMETER (MAXANTHC=10, MAXBASHC=((MAXANTHC*(MAXANTHC-1))/2))
c
c MAXUVPNT -- maximum number of UV points in time per track
c
      INTEGER   MAXUVPNT
      PARAMETER(MAXUVPNT=2000)
c
c  MAXSRC -- maximum number of sources in multi-source databases
c
      INTEGER MAXSRC
      PARAMETER( MAXSRC = 32 )
c
c  MAXBREAK -- maximum number of breakpoints on an interval
c
      INTEGER MAXBREAK
      PARAMETER( MAXBREAK = 10 )
c
c  MAXPOLY -- maximum number of polynomials on an interval
c             is normally one more than number of breakpoints, right?
c
      INTEGER MAXPOLY
      PARAMETER( MAXPOLY = MAXBREAK+1 )
c
c  MAXSLOT -- maximum number of polynomial slots
c             Allows for simultaneous presence of gcal & pcal format
c             The first 4 for gcal, last 16 for pcal
c             CALIB makes use of slot 5 and 6 for phase diff&mean
      INTEGER MAXSLOT
c      PARAMETER( MAXSLOT = 20 )
      PARAMETER( MAXSLOT = 36 )
c
c  MAXORDER -- highest allowed polynomial order
c             Don't set too high, inversions are already needed in
c             double precision this way.
      INTEGER MAXORDER
      PARAMETER( MAXORDER = 5 )
c
c  MAXWIN -- maximum number of windows - normally equals MAXSPECT = 8
c	     but temporarely increased to 16 to deal with mode=4
c	     problems.
c     INTEGER MAXWIN
c     PARAMETER ( MAXWIN = 16 )
c
c=======================================================================

