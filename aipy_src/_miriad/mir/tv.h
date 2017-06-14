c ---------------------------------------------------------------------
c  History
c   jm  xxxxxxx Original version.
c  mjs  11apr91 Removed ultra/cray and raster/fx parameters.
c   jm  16apr91 Replaced ultra/cray and raster/fx parameters.
c   jm  20apr91 Added mxas driver.
c   jm  24apr91 Modified TVzscr common block.
c   jm  25apr92 Updated mxas parameters (Port# 5500 -> 5000).
c   jm  02may92 Added mxas parameter Interrogate (XasIntgt=13).
c   jm  15jul92 Added Graph parameter and server xmtv.
c ---------------------------------------------------------------------
c
c  Variables in Common TVcomm:
c    Handle	Some sort of handle or file descriptor used by the TCP software.
c    BufSize	Size of Buffer.
c    LenBuf	The number of integers stored in Buffer.
c    Nack	Number of acknowledgements to expect in the input
c		stream.
c    Protocol	The device type currently being used.
c
c  Variables in Common TVzscr:
c    mxzoom  	The maximum zoom in powers of 2.  If ``mxzoom'' is
c               less than 0, then zoom is linear from 1 to abs(mxzoom).
c    tvscrx  	The current amount of scroll in the X direction.
c    tvscry  	The current amount of scroll in the Y direction.
c
	integer Ivas,Sss,Ivserve,File,VFile,Ultra,Raster,Xas
	parameter(Ivas=1,Sss=2,Ivserve=3,File=4,VFile=5)
	parameter(Ultra=6,Raster=7,Xas=8)
	integer IvPort,SssPort,XasPort
	parameter(IvPort=9,SssPort=5000,XasPort=5000)
c
c  Parameters for the Sun server.
c
	integer SssOpen,SssSplit,SssClear,SssImWrt
	integer SssRCurs,SssRButt,SssWZscr,SssWLut
	integer SssSelpt,SssText,SssWindo,SssRZscr
	integer SssGraph
	parameter(SssOpen=11,SssSplit=46,SssClear=15,SssImWrt=21)
	parameter(SssRCurs=61,SssRButt=62,SssWZscr=83,SssWLut=41)
	parameter(SssSelpt=71,SssText=53,SssWindo=14,SssRZscr=84)
	parameter(SssGraph=45)
c
c  Parameters for the IVAS server.
c
	integer IvClose,IvVPsetU,IvMaImag,IvGPHset,IvInit
	integer IvMoStat,IvVPZScr,IVlocal
	integer PassIn,PassByte,PassInt
	parameter(IvClose=0,IvVPsetU=93,IvMaImag=67,IvGPHset=57)
	parameter(IvMoStat=77,IvVPZScr=94,IVlocal=9)
	parameter(IvInit=62,PassIn=1,PassByte=1,PassInt=3)
c
c  Parameters for the X servers (both the MXAS and XMTV servers
c  honor these requests).
c
	integer XasOpen,XasClose,XasIntgt,XasWindo,XasClear
	integer XasImWrt,XasScale,XasWLut,XasGraph,XasSplit
	integer XasText,XasRCurb
	integer XasSelpt,XasWZscr,XasRZscr
	parameter(XasOpen=11,XasClose=12,XasIntgt=13,XasWindo=14)
	parameter(XasClear=15,XasImWrt=21,XasScale=29)
	parameter(XasWLut=41,XasGraph=45,XasSplit=46)
	parameter(XasText=53,XasRCurb=64)
	parameter(XasSelpt=71,XasWZscr=83,XasRZscr=84)
c
	integer bufsize,BypWrd
#ifdef unicos
	parameter(BufSize=2048,BypWrd=8)
#else
	parameter(BufSize=1024,BypWrd=4)
#endif
	integer buffer(BufSize)
	integer BufLen,handle,protocol,Nack
	integer MxZoom, TvScrx, TvScry
	integer LastMag, LastX, LastY
	common/TVcomm/handle,BufLen,protocol,Nack,buffer
	common /TVzscr/ MxZoom, TvScrx, TvScry, LastMag, LastX, LastY
