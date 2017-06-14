c=======================================================================
c	[calsubs.h]
c
c	Last update:    22-may-90	deleted pmode	PJT
c			18-dec-90	played with  /calmode/ [Arie/Peter]
c			 7-jul-91       new MAXBASHC parameter
c			11-dec-91	documentation formalized in doc format
c
c*calsubs_h -- header file for calibration
c:calibration, include
c&pjt
c+
c
c  time0 -- Julian time offset for all time measurements
	DOUBLE PRECISION time0

c  base(nbl) -- the baseline number and order
	INTEGER nbl, base(MAXBASHC)

c  version -- the read version of cal format
	INTEGER version

c  calsubs -- THE common block
	COMMON /calsubs/ time0, nbl, base, version

c
c  kludge before upgrading to the next CAL-SVERSION
	INTEGER scalmode
c
c	0:	(old system) fitting is done in K
c	1:	fitting is done in data space (K/Jy), i.e. the same is above
c      -1:	fitting is done in gain space (Jy/K)
c
c   For now, scalmode is set in calfit and calflag and such.
c   CAopen should get an extra stamp for this, see calio.h
c
	COMMON /calscal/scalmode
c--
c   Mode of phase/amp calibration in calflag/calib/...
c   [may be used in Arie's calib.for - but is internal for now]
c   nmode=1   amp
c	  2   phase
c         3   diff
c
c	INTEGER nmode
c	COMMON /calmode/nmode
c=======================================================================

