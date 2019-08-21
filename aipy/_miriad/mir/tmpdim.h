c=======================================================================
c - MAXBUF tells us how many words of memory we can use for data
c  History:
c    mchw 06mar91  Changed maxbuf on ral
c    mjs  29jun92  maxdim of 2048 at UIUC/astro
c    mjs  18feb93  maxchan of 2048 at UIUC/astro
c    mjs  26feb93  lower maxdim for memory hog pgms
c
	INTEGER   MAXBUF
	PARAMETER(MAXBUF=1000000)
c-----------------------------------------------------------------------
c - MAXDIM is an often used parameter, to indicate maximum size of maps
	INTEGER   MAXDIM
	PARAMETER(MAXDIM=1024)
c-----------------------------------------------------------------------
c		maximum number of antennae (HC=3/6/9/..., WSRT=14, VLA=27)
	INTEGER   MAXANT
	PARAMETER(MAXANT=6)

c		maximum number of baselines
	INTEGER   MAXBASE
	PARAMETER(MAXBASE=((MAXANT*(MAXANT-1))/2))

c		maximum number of channels in spectral data
	INTEGER   MAXCHAN
	PARAMETER(MAXCHAN=513)
c=======================================================================
