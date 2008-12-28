	include 'maxdim.h'
	include 'maxnax.h'
c
c  Parameters setting internal buffers.
c  NOTE: We should ensure that 13*MAXCHAN is greater than MAXDIM.
c
	integer MAXSIZE,MAXCARDS,MAXOPEN
	parameter(MAXSIZE=13*MAXCHAN,MAXCARDS=16,MAXOPEN=4)
c
c  Parameters associated with the FITS file as a whole.
c  NOTE: DatOff,DatBase and new are used by higher level routines.
c
	integer curlu,curcard
	integer DatSize(MAXOPEN),HdSize(MAXOPEN)
	integer DatOff(MAXOPEN),HdOff(MAXOPEN)
	integer DatBase(MAXOPEN)
	integer item(MAXOPEN)
	integer ncards(MAXOPEN)
	logical new(MAXOPEN),opened(MAXOPEN)
	character carray*(80*MAXCARDS)
c
c  Parameters common to images (XY) and visibility (UV) files.
c
	integer dBypPix
	parameter(dBypPix=-4)
	integer BypPix(MAXOPEN)
	logical float(MAXOPEN)
	real bscale(MAXOPEN),bzero(MAXOPEN)
	integer BlankVal(MAXOPEN)
c
c  For Image files.
c
	integer axes(MAXNAX,MAXOPEN),PixBase(MAXOPEN)
c
c  For Visibility files.
c
	integer MAXRAN
	parameter(MAXRAN=10)
	real scales1(MAXRAN,MAXOPEN),scales2(MAXRAN,MAXOPEN)
	real zeros(MAXRAN,MAXOPEN)
	integer indices1(MAXRAN,MAXOPEN),indices2(MAXRAN,MAXOPEN)
	double precision TimOff(MAXOPEN)
	real WtScal(MAXOPEN)
	integer nRanFile(MAXOPEN),nRanProg(MAXOPEN)
	integer ncomplex(MAXOPEN),pols(MAXOPEN),freqs(MAXOPEN)
	integer visibs(MAXOPEN)
	integer nts1(MAXOPEN),nts2(MAXOPEN),nts3(MAXOPEN)
c
c  Parameters related to arrays returned by the uv routines.
c
	integer uvCrval,uvCdelt,uvCrpix
	integer uvStokes,uvFreq,uvRa,uvDec
	parameter(uvCrval=1,uvCdelt=2,uvCrpix=3)
	parameter(uvStokes=1,uvFreq=2,uvRa=3,uvDec=4)
c
c  Scratch arrays used by all and sundry.
c
	integer array(MAXSIZE),arrayd(MAXSIZE)
	real rarray(MAXSIZE),rarrayd(MAXSIZE)
c
c  The common block which contains this mess.
c
	common/fitscom/TimOff,curlu,curcard,item,bscale,bzero,
     *	  scales1,scales2,zeros,axes,DatOff,HdOff,DatSize,HdSize,
     *	  DatBase,PixBase,BypPix,BlankVal,ncards,nRanFile,nRanProg,
     *	  ncomplex,indices1,indices2,visibs,nts1,nts2,nts3,pols,freqs,
     *	  array,arrayd,WtScal,rarray,rarrayd,
     *	  new,opened,float
	common/fitscomc/carray
c
c  Info to help find an extension table that we are interested in.
c
c  MAXIDX	The size of the index.
c  ExtName	Name of an extension table (e.g. 'AIPS CC ')
c		as given by the EXTNAME keyword.
c  ExtOff	Offset to start of this tables header.
c  nExtOff	Number of entries in the index.
c  ExtNo	The index of the current loaded table. 0 indicates
c		the main header/data.
c
	integer MAXIDX
	parameter(MAXIDX=16)
	character ExtName(MAXIDX,MAXOPEN)*8
	integer ExtOff(MAXIDX,MAXOPEN),nExtOff(MAXOPEN),ExtNo(MAXOPEN)
	common/fitsidx/ExtOff,nExtOff,ExtNo
	common/fitsidxc/ExtName
c
c  Variables for i/o on an "A3DTABLE"
c
c  MAXCOL	Max number of columns a table can have.
c  rows		Number of rows in the table.
c  cols		Number of columns in the table
c  width	The width of a row of a table, in bytes.
c  ColForm	Gives the data type of a column. Its value with be
c		one of the "Form" parameter values.
c  ColCnt	Size of the data in this column, for one row, in bits.
c  ColOff	Offset, in bytes, to the start of the data in each row.
c  ColType	The "TTYPE" value.
c  ColUnits	The "TUNIT" value.
c  
	integer MAXCOL,NFORMS
	parameter(MAXCOL=400,NFORMS=10)
	integer rows(MAXOPEN),cols(MAXOPEN),width(MAXOPEN)
	integer ColForm(MAXCOL,MAXOPEN),ColCnt(MAXCOL,MAXOPEN)
	integer ColOff(MAXCOL,MAXOPEN)
	character ColType(MAXCOL,MAXOPEN)*32
	character ColUnits(MAXCOL,MAXOPEN)*16
c
        integer FormI,FormJ,FormA,FormE,FormD,FormX,FormL,FormC
	integer FormM,FormP
	parameter(FormJ=1,FormI=2,FormA=3,FormE=4,FormD=5,FormX=6,
     *            FormL=7,FormC=8,FormM=9,FormP=10)
	common/fitstab/rows,cols,width,ColForm,ColCnt,ColOff
	common/fitstabc/ColType,ColUnits
