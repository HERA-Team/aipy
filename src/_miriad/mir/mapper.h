c
c  Gridding related parameters:
c
c    MAXCGF	Maximum size of the convolutional gridding function
c		array.
c    n1,n2	Size of the transform.
c    width	Width of the gridding function in uv-plane cells.
c    ncgf	Size of the cgf array.
c    cgf	The convolutional gridding function.
c    xcorr,ycorr Correction arrays.
c
	include 'maxdim.h'
	integer MAXPNT,MAXCGF,MAXT
	parameter(MAXPNT=2048,MAXCGF=2048,MAXT=5)
c
	real scale(MAXPNT)
	real cgf(MAXCGF),xcorr(MAXDIM),ycorr(MAXDIM),umax,vmax
	integer tscr,nvis
	integer width,ncgf,offcorr,chan1,chan2,npnt,totchan
	integer nchan(MAXT),nx(MAXT),ny(MAXT),nt
	integer n1,n2,nu,nv,u0,v0,nextra,nxc,nyc
	logical ginit
	integer nBuff,pBuff
	character mode*8
	common/mapcom/scale,cgf,xcorr,ycorr,umax,vmax,
     *	  tscr,nvis,width,ncgf,offcorr,chan1,chan2,npnt,totchan,
     *	    nchan,nx,ny,nt,n1,n2,nu,nv,u0,v0,nextra,nxc,nyc,nBuff,pBuff,
     *	  ginit
	common/mapcomc/mode
c
	integer num
	common/mapcom2/num
