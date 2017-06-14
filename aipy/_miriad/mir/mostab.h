	integer MAXPNT,MAXHASH
	parameter(MAXPNT=2048,MAXHASH=4*MAXPNT+1)
	integer nxy,pX,pY,coRef
	integer pntno,npnt,vPntUpd,HashSize,Hash(MAXHASH),nx2,ny2
	integer pbObj(MAXPNT)
	logical solar,doinit
	double precision ucoeff(3,MAXPNT),vcoeff(3,MAXPNT)
	real Rms2(MAXPNT),SumWt(MAXPNT),x0(MAXPNT),y0(MAXPNT)
	double precision llmm(2,MAXPNT),radec(2,MAXPNT),radec0(2)
	double precision cdelt1,cdelt2
	character telescop(MAXPNT)*16
	common/mostab1/	llmm,radec,radec0,cdelt1,cdelt2,
     *	  ucoeff,vcoeff,Rms2,SumWt,x0,y0,
     *	  pntno,npnt,nxy,pX,pY,vPntUpd,HashSize,Hash,nx2,ny2,coRef,
     *	  pbObj,
     *	  solar,doinit
	common/mostab2/telescop
