	integer MAXCOUNT
	parameter(MAXCOUNT=101)
	logical doinit,mfs,mosaic
	real lstart,lwidth,lstep,vobs,epoch
	double precision crval1,crval2,crval3,cdelt1,cdelt2,cdelt3
	double precision obsra,obsdec,restfreq,obstime
	character ctype1*12,ctype2*12,ctype3*12,ltype*12,telescop*12,
     *	  source*16,observer*16,pbtype*16
	integer Count,naver,nchan
	logical VChange,VLinear,Rconst,RChange,XChange,YChange,SChange
	common/hdtab1/crval1,crval2,crval3,cdelt1,cdelt2,cdelt3,
     *	  obsra,obsdec,restfreq,obstime,
     *	  lstart,lwidth,lstep,epoch,vobs,
     *	  Count,naver,nchan,doinit,mfs,mosaic,
     *	  VChange,VLinear,Rconst,RChange,XChange,YChange,SChange
	common/hdtab2/ctype1,ctype2,ctype3,ltype,telescop,source,
     *		observer,pbtype
