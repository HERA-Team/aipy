	integer MAXSRC,MAXWTS
	parameter(MAXSRC=10000,MAXWTS=100)
	integer nsrc,nwts
	double precision lmn(3,MAXSRC),freq0
	double precision ucoeff(3),vcoeff(3),wcoeff(3),cosd,sind
	real flux(2,MAXSRC),factor(MAXSRC),inttime,wts(MAXWTS)
	common/modpcom/lmn,ucoeff,vcoeff,wcoeff,cosd,sind,freq0,
     *	  flux,factor,inttime,wts,nwts,nsrc
