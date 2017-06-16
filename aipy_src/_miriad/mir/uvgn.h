c
c Variables to handle normal gains and delays associated with a file:
c
c  t1,t2	Pointers into the solution table.
c  nsols	Number of gain solutions in the gains file.
c  ngains	Total number of gains = (ntau + nfeeds) * nants
c  nfeeds	Number of feeds in each gain solution.
c  nants	Number of antennae in each gain solution.
c  ntau		Number of delay terms.
c  gitem	The item of the gains file.
c  solno	Gives the solution number of the solutions in memory.
c  timetab	Gives the time of the gain solutions in memory.
c  gains	The portion of the gain solutions in memory.
c  gflag	Flags of whether the corresponding gain is good.
c------------------------------------------------------------------------
c  MCHW's baseline/channel number based bandpass correction.
c
c  docgains	Apply channel gains to channel data.
c  dowgains	Apply wideband gains to wideband data.
c  nwbase	The number of baselines for complex wideband gains.
c  ncbase	The number of baselines for complex channel gains.
c  nwgains	The number of complex wideband gains.
c  ncgains	The number of complex channel gains.
c  wgains	The complex wideband gains.
c  cgains	The complex channel gains.
c------------------------------------------------------------------------
c  RJS's antenna/frequency based bandpass correction.
c
c  vwide	UV variable handle for wfreq,wwidth.
c  vline	UV variable handle for sfreq,sdf,nschan.
c  dopass
c  aver
c  
c
	include 'maxdim.h'
	integer MAXTAB,MAXFEEDS,MAXGAINS,MAXSPECT
	parameter(MAXTAB=2,MAXFEEDS=2,MAXSPECT=32)
	parameter(MAXGAINS=3*MAXANT)
	integer t1,t2,nsols,nants,nfeeds,ntau,ngains,gitem,solno(MAXTAB)
	double precision timetab(MAXTAB),dtime
	complex gains(MAXGAINS,MAXTAB)
	logical gflag(MAXGAINS,MAXTAB),dogains,dotau
c
	integer ncgains,ncbase,nwgains,nwbase
	logical docgains,dowgains
	integer pCgains,pWgains
c
	logical dopass,aver,first
	integer tno,vwide,vline,nchan,nspect,nschan(MAXSPECT)
	double precision sfreq(MAXSPECT),sdf(MAXSPECT),freq0
	integer pFlags(2),pDat(2),nDat(2),pTab,nTab,pFreq(2),nFreq(2)
c
c
c  The common blocks.
c
	common/UvGnA/timetab,dtime,gains,t1,t2,nsols,nants,nfeeds,
     *	  ngains,ntau,gitem,solno,gflag,dogains,dotau
c
	common/UvGnB/ncgains,nwgains,ncbase,nwbase,pCgains,pWgains,
     *		docgains,dowgains
c
	common/uvGnC/sfreq,sdf,freq0,tno,vwide,vline,nchan,nspect,
     *	  nschan,pFlags,pDat,nDat,pTab,nTab,pFreq,nFreq,dopass,
     *	  aver,first

