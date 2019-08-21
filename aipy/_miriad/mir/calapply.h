c=======================================================================
c   [calappply.h]
c		dark ages	bs/lgm/pjt, in some form
c		17-dec-91	documentation format
c
c*calapply_h -- header file for calibration data
c:calibration, include
c&pjt
c+
c    Include file for reading and writing uv data files used by
c    subroutines ... READNEXT and WRITNEXT ... ???
c    It also has a common block /calcnt/ keeping track of good and
c    bad scans
c
c    Space used:  appr. 2 * MAXCHAN
c

      DOUBLE PRECISION preamble(4)
      COMPLEX          data(MAXCHAN)
      LOGICAL          flags(MAXCHAN), wflags(MAXWIDE)
      INTEGER          nread,nwcorr
      COMPLEX          wcorr(MAXWIDE)
      COMMON /uvdata/ preamble, data, flags, nread, wcorr, nwcorr,
     *                wflags

      INTEGER calcntg, calcntb
      COMMON /calcnt/calcntg, calcntb
c--
c  History:
c    9feb99 changed wflags and wcorr from 8 to MAXWIDE
c=======================================================================
