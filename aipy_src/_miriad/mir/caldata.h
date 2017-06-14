c=======================================================================
c       [caldata.h]
c       Header file for calibration data
c
c       Last update:    14-mar-90       added Sflag     PJT
c                       20-mar-90       srcindex now srcindex(MAXPNT)   PJT
c                        7-may-90       rdata(ULRI,MAXBASHC,MAXUVPNT)   PJT
c                        4-sep-90       baseline based calib stuff      PJT
c                       18-nov-90       more descriptive comments       PJT
c                       18-dec-90       played with rdata()             PJT
c                        7-jul-91       MAXBASE -> BAXMASHC             PJT
c                        2-dec-91       make fluxes time dependant, not
c                                       source dependant; =V6.4 fmt     PJT
c                       11-dec-91       documentation formalized in doc format
c*caldata_h -- header file for calibration data
c:calibration, include
c&pjt
c+
c   The calibration format on disk is split out in memory over several
c   common blocks, each of them described in a separate include file.
c
c
c <general header> (see also calsubs.h)
c       double time0            time to which all times are offset (JD)
c       int    nbl              number of baselines
c       int    base(nbl)        nbl baselines (notation: A1*256+A2)
c       int    version          version of the calibration format
c
c <data> (this caldata.h include file)
c   rtime                       
c       float  rtime(rcount)             offset times to time0
c   rdata
c       float  rdata(ULRI,nbl,rcount)    correllation data
c   rflag
c       int    rflag(UL,nbl,rcount)      flagged bad? (0=false, 1=true)
c   sindex
c       int    sindex(rcount)            source index number
c   sname
c       char*8 sname(scount)             list of sources
c   

c
c  Rxxxx -- the actual calibration data
c
      INTEGER   rcount
      REAL      rtime(MAXUVPNT)
      REAL      rdata(ULRI, MAXBASHC, MAXUVPNT)
      INTEGER   rflag(UL, MAXBASHC, MAXUVPNT)
c
c  Sxxxx -- index for multiple source calibration format
c           (plstuff was added later and is not quite needed??)
c
      INTEGER   scount
      CHARACTER sname(MAXSRC)*8
      INTEGER   sindex(MAXUVPNT), sflag(MAXUVPNT)
      REAL      plstuff(4,MAXSRC)
c
c  Bxxxx - list of break points
c
      REAL      btime(MAXBREAK, 2, MAXBASHC)
      INTEGER   bcount(2, MAXBASHC)
c
c  Xxxxx - odds and ends, shouldn't really have to exist
c          are merely there for planets
c
      REAL      volts(MAXUVPNT)
c
c  baseline based calibrator fluxes, per time
c
      REAL      calbflux(MAXBASHC,MAXUVPNT)

c
c  the common blocks with all the data
c
      COMMON /caldatr/ rcount, scount, rtime, rdata, rflag,
     -                   sindex, btime, bcount, sflag, plstuff,volts
      COMMON /caldatc/ sname
      COMMON /caldatp/ calbflux
c--
c=======================================================================

