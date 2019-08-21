	integer MAXPNT
	parameter(MAXPNT=2048)
	integer cnvl(MAXPNT),npix,npnt,n1,n2,n1d,n2d,ic,jc
	integer pWrk1,pWrk2,nWrk,pWts1,pWts2,nWts,tno
	logical mosini,dogaus
	real bmaj,bmin,bpa
	character flags*8
c
	common/mccom/bmaj,bmin,bpa,cnvl,tno,npix,npnt,
     *	    n1,n2,n1d,n2d,ic,jc,
     *	    pWrk1,pWrk2,nWrk,pWts1,pWts2,nWts,mosini,dogaus
	common/mccomc/flags
