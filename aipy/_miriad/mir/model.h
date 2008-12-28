c************************************************
c	model.h - include file for model.for
c------------------------------------------------
	integer maxcals
	parameter(maxcals=32)
	integer ncals,calcur
	logical planet(maxcals)
	character cals(maxcals)*16
	real calflux(maxcals),CalFreq(maxcals)
	common/calcom/ncals,calcur,calflux,calfreq,planet
	common/calcomc/cals
