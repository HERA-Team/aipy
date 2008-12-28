c
c  uvdat.h - Common blocks and parameters for the uvDat routines.
c
c	History:
c
c	??????? mchw Original version.
c	04aug91 mjs  Replaced MAXANT(S) by including maxdim.h
c       29mar94 nebk Add logical CALMSG
c       10dec95 nebk Extend CALMSG to an array
c
c  General variables:
c
c    sels	The selection specification to pass to SelApply.
c    line,nchan,lstart,lwidth,lstep	The data linetype specification.
c    ref,rstart,rwidth		The reference linetype specification.
c    plmaj,plmin,plangle	Planet parameters.
c    InBuf(k1(i):k2(i))		The name of the i'th input file.
c    doplanet			True if planet processing is to be performed.
c    plinit			True if the planet parameters have been
c				initialised.
c    dosels			True if uvselection is to be applied.
c    dowave			True if u-v coords are to be given in
c				wavelengths.
c    doref			True if there is to be a reference line.
c    dodata			True if we are to perform data linetype
c				processing.
c    docal			Check for a calibration file, and apply
c				the calibration corrections.
c    WillCal			Set true if docal is true and  cal files
c				are present.
c    dopass			Check for a bandpass file, and apply the
c				bandpass corrections.
c    doleak			Check for polarisation leakage file, and
c				apply the polarisation corrections.
c    WillLeak			Set true if doleak is true and leakage files
c				are present.
c    dow			Return "w" in the preamble.
c    nIn			The number of input files given.
c    pnt			The number of the current input file being
c				processed.
c    tno			The handle of the current input file being
c				processed. If 0, then there is no current
c				file.
c    auto			If true, input data must be autocorrelation
c				data.
c    cross			If true, input data must be crosscorrelation
c				data.
c    calmsg                     If true, calibration message already issued
c			        for this file.
c    dogsv			Enable scaling of the variance by the
c				gains.
c    npream			Number of elements in the preamble.
c    idxT			Index of "time" in the preamble.
c    idxBL			Index of "baseline" in the preamble.
c
	integer maxsels,maxNam,maxIn
	parameter(maxsels=1024,maxNam=20000,maxIn=400)
	real sels(maxsels),lstart,lwidth,lstep,lflag,rstart,rwidth,rstep
	real plmaj,plmin,plangle
	logical doplanet,dowave,doref,dodata,docal,dosels,doleak,dopass
	logical PlInit,WillCal,WillLeak,auto,cross,calmsg(maxIn),dow
	logical dogsv
	integer k1(maxIn),k2(maxIn),nchan,npream,idxT,idxBL
	character line*32,ref*32,InBuf*(maxNam)
	integer nIn,pnt,tno
c
c  Variables to handle polarisation conversion.
c
c  nPol	    The number of polarisations desired by the user.
c  nPolF    Number of simultaneous polarisations in the file. This
c	    will be zero if we cannot determine it, or if it varies.
c  SelPol   Logical, being true if the user has used polarisation
c	    selection.
c  SelPol1  Logical, being true if the user selected just 1 polarisation.
c  Pols	    The polarisations desired (in the appropriate codes).
c  WillPol  If true, then polarisation must be determined the hard way.
c  PolCpy   Do polarisation processing, but just copy out the polarisations
c	    found in the file.
c  iPol	    Pols(iPol) is the current polarisation of interest.
c  Snread   The number of channels.
c  SData    The visibility data. The file can contain up to 4 different
c	    sorts of polarisations/Stokes parameters.
c  Sflags   The flags corresponding to SData.
c  Spreambl The saved preamble.
c  ncoeff ) These are used to convert from the polarisation present in the
c  indices) file to the polarisation that the caller desires. See uvGetPol
c  coeffs ) for more info.
c  doaver )
c  GWts     The gain weights to apply to the data.
c  SumWts   Sum of the polarisation weights squared -- for noise calculations.
c  Leaks    Polarisation leakage parameters.
c  nLeaks   Number of polarisation leakage parameters.
c
	include 'maxdim.h'
	integer MAXPOL,MAXPRE
	parameter(MAXPOL=4,MAXPRE=8)
c
	integer PolII,PolI,PolQ,PolU,PolV,PolRR,PolLL,PolRL,PolLR
	integer PolXX,PolYY,PolXY,PolYX,PolQQ,PolUU,PolMin,PolMax
	parameter(PolII=0,PolI=1,PolQ=2,PolU=3,PolV=4,PolRR=-1)
	parameter(PolLL=-2,PolRL=-3,PolLR=-4,PolXX=-5,PolYY=-6)
	parameter(PolXY=-7,PolYX=-8,PolQQ=5,PolUU=6)
	parameter(PolMin=PolYX,PolMax=PolUU)
c
	integer nPol,Pols(MAXPOL),iPol,nPolF
	logical WillPol,SelPol,SelPol1,PolCpy
	double precision Spreambl(MAXPRE)
	integer Snread
	complex SData(MAXCHAN,MAXPOL)
	integer ncoeff(MAXPOL),indices(MAXPOL,MAXPOL)
	logical doaver(MAXPOL),Sflags(MAXCHAN,MAXPOL)
	complex coeffs(MAXPOL,MAXPOL)
	real SumWts(MAXPOL),GWt
	integer nLeaks
	complex Leaks(2,MAXANT)
c
c  The common blocks.
c
	common/UVDatCoA/sels,lstart,lwidth,lstep,lflag,
     *	 rstart,rwidth,rstep,
     *	 plmaj,plmin,plangle,doplanet,dowave,doref,dodata,dosels,dow,
     *	 dogsv,plinit,k1,k2,nchan,nIn,pnt,tno,npream,idxT,idxBL,
     *	 auto,cross,docal,WillCal,doleak,WillLeak,dopass,calmsg
c
	common/UVDatCoB/line,ref,InBuf
c
	common/UvDatCoC/Spreambl,Leaks,coeffs,SumWts,GWt,WillPol,
     *	 PolCpy,
     *	 SelPol,SelPol1,nPol,nPolF,Pols,iPol,Snread,SData,ncoeff,doaver,
     *	 Sflags,indices,nLeaks
