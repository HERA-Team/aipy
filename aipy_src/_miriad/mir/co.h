	include 'maxnax.h'
c
c  The various sorts of coordinates.
c
	integer LAT,LON,VELO,FELO,FREQ,LINEAR
	parameter(LAT=1,LON=2,VELO=3,FELO=4,FREQ=5,LINEAR=6)
c
	integer MAXOPEN
	parameter(MAXOPEN=15)
c
	integer Lus(MAXOPEN),nalloc(MAXOPEN)
	integer naxis(MAXOPEN)
	integer ilong(MAXOPEN),ilat(MAXOPEN),ifreq(MAXOPEN)
	double precision crval(MAXNAX,MAXOPEN)
	double precision crpix(MAXNAX,MAXOPEN)
	double precision cdelt(MAXNAX,MAXOPEN)
	double precision restfreq(MAXOPEN),vobs(MAXOPEN)
	double precision epoch(MAXOPEN),obstime(MAXOPEN)
	double precision llcos(MAXOPEN),llsin(MAXOPEN)
	character ctype(MAXNAX,MAXOPEN)*16,coproj(MAXOPEN)*3
	integer cotype(MAXNAX,MAXOPEN)
	logical cellscal(MAXOPEN)
c
	common/cocom/crval,crpix,cdelt,restfreq,vobs,epoch,obstime,
     *		llcos,llsin,
     *		lus,nalloc,naxis,cotype,ilong,ilat,ifreq,cellscal
	common/cocomc/ctype,coproj
