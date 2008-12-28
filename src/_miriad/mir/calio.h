c=======================================================================
c	[calio.h]
c
      INCLUDE 'size.h'	
c*calio -- calibration include file
c&pjt
c:calibration, include
c+
c  calio.h -- Calibration I/O Routines
c
c
c
c  CAnbl   -- number of baselines
c  IRtime  -- Rtime item handle
c  IRdata  -- Rdata item handle
c  IRflag  -- Rflag item handle
c  ISname  -- sname item handle
c  ISindex -- sindex item handle
c  IPspan  -- Polynomial time limit item handle
c  IPdata  -- Polynomial data item handle
c  IPindex -- Polynomial indicies item handle
c  IBdata  -- Breakpoints item handle
c
c
c  Required space:    10*MAXOPEN 
c
      INTEGER MAXOPEN
      PARAMETER( MAXOPEN = 10 )

      INTEGER canbl(MAXOPEN)
      INTEGER irtime(MAXOPEN), irdata(MAXOPEN), irflag(MAXOPEN)
      INTEGER isname(MAXOPEN), isindex(MAXOPEN)
      INTEGER ipspan(MAXOPEN), ipdata(MAXOPEN), ipindex(MAXOPEN)
      INTEGER ibdata(MAXOPEN)
c
c  CA -- the common block
c
      COMMON /ca/ canbl, 
     *              irtime, irdata, irflag, 
     *              isname, isindex,
     *              ipspan, ipdata, ipindex,
     *              ibdata
c--
c=======================================================================

