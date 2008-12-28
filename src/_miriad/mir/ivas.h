c
c  Graphics processor definitions
c
	integer GPHrFIFOentry,GPHrCmdCtl,GPHrOpMode,GPHrDispCtl
	parameter( GPHrFIFOentry        = 0 )
	parameter( GPHrCmdCtl           = 2 )
	parameter( GPHrOpMode           = 4 )
	parameter( GPHrDispCtl          = 6 )
c
	integer GPHrRasterCount,GPHrHorizSync,GPHrHorizDisp
	integer GPHrVertSync,GPHrVertDisp,GPHrSplitScreen
	integer GPHrBlinkCtl,GPHrHorizWindow,GPHrVertWindow
	integer GPHrGraphicCursor
c
	parameter( GPHrRasterCount       = 128 )
	parameter( GPHrHorizSync         = 130 )
	parameter( GPHrHorizDisp         = 132 )
	parameter( GPHrVertSync          = 134 )
	parameter( GPHrVertDisp          = 136 )
	parameter( GPHrSplitScreen       = 138 )
	parameter( GPHrBlinkCtl          = 144 )
	parameter( GPHrHorizWindow       = 146 )
	parameter( GPHrVertWindow        = 148 )
	parameter( GPHrGraphicCursor     = 152 )
c
	integer GPHrRasterAddr0,GPHrMemWidth0,GPHrStartAddr0
	integer GPHrRasterAddr1,GPHrMemWidth1,GPHrStartAddr1
	integer GPHrRasterAddr2,GPHrMemWidth2,GPHrStartAddr2
	integer GPHrRasterAddr3,GPHrMemWidth3,GPHrStartAddr3
	integer GPHrBlockCur1,GPHrBlockCur2,GPHrCursorDef
	integer GPHrZoomFactor,GPHrLightPenAddr
c
	parameter( GPHrRasterAddr0       = 192 )
	parameter( GPHrMemWidth0         = 194 ) 
	parameter( GPHrStartAddr0        = 196 )
	parameter( GPHrRasterAddr1       = 200 )
	parameter( GPHrMemWidth1         = 202 )
	parameter( GPHrStartAddr1        = 204 )
	parameter( GPHrRasterAddr2       = 208 )
	parameter( GPHrMemWidth2         = 210 )
	parameter( GPHrStartAddr2        = 212 )
	parameter( GPHrRasterAddr3       = 216 )
	parameter( GPHrMemWidth3         = 218 )
	parameter( GPHrStartAddr3        = 220 )
	parameter( GPHrBlockCur1         = 224 )
	parameter( GPHrBlockCur2         = 228 )
	parameter( GPHrCursorDef         = 232 )
	parameter( GPHrZoomFactor        = 234 )
	parameter( GPHrLightPenAddr      = 236 )
c
	integer GPHsWriteEmpty,GPHsWriteReady,GPHsReadReady
	integer GPHsReadFull,GPHsPenStrobe,GPHsCommandEnd
	integer GPHsAreaDetect,GPHsCommandError
c
	parameter( GPHsWriteEmpty        = 1 )
	parameter( GPHsWriteReady        = 2 )
	parameter( GPHsReadReady         = 4 )
	parameter( GPHsReadFull          = 8 )
	parameter( GPHsPenStrobe         = 16 )
	parameter( GPHsCommandEnd        = 32 )
	parameter( GPHsAreaDetect        = 64 )
	parameter( GPHsCommandError      = 128 )
c
	integer GPHarNoChk,GPHarExAbort,GPHarExSupp,GPHarExDetect
	integer GPHarEnAbort,GPHarEnSupp,GPHarEnDetect
c
	parameter( GPHarNoChk            = 0 )
	parameter( GPHarExAbort          = 32 )
	parameter( GPHarExSupp           = 64 )
	parameter( GPHarExDetect         = 96 )
	parameter( GPHarEnAbort          = 160 )
	parameter( GPHarEnSupp           = 192 )
	parameter( GPHarEnDetect         = 224 )
c
	integer GPHclBoth,GPHclOne,GPHclZero,GPHclRAM
c
	parameter( GPHclBoth             = 0 )
	parameter( GPHclOne              = 8 )
	parameter( GPHclZero             = 16 )
	parameter( GPHclRAM              = 24 )
c
	integer GPHopRep,GPHopOR,GPHopAND,GPHopEOR,GPHopEQL,GPHopNEQ
	integer GPHopLT,GPHopGT
c
	parameter( GPHopRep              = 0 )
	parameter( GPHopOR               = 1 )
	parameter( GPHopAND              = 2 )
	parameter( GPHopEOR              = 3 )
	parameter( GPHopEQL              = 4 )
	parameter( GPHopNEQ              = 5 )
	parameter( GPHopLT               = 6 )
	parameter( GPHopGT               = 7 )
c
	integer GPHabs,GPHrel
	parameter( GPHabs              = 0 )
	parameter( GPHrel              = 1024 )
c
c Video pipeline definitions
c
	integer GMdisable,GMcolor,GMroi,GMmSplit4,GMmSplit2,GMmSplit1
	integer GMmSplit3,GMmWindow,GMmBlinkOn,GMmBlinkOff
c
	parameter( GMdisable=-1 )
	parameter( GMcolor  = 1 )
	parameter( GMroi    = 2 )
c
	parameter( GMmSplit4    = 1 )
	parameter( GMmSplit2    = 2 )
	parameter( GMmSplit1    = 4 )
	parameter( GMmSplit3    = 8 )
	parameter( GMmWindow    = 16 )
c
	parameter( GMmBlinkOn      = 32768 )
	parameter( GMmBlinkOff     = 16384 )
c
	integer PassByte,PassRaw,PassInt,PassLong,PassFloat,PassString
	integer PassIn,PassOut,PassNoWait
	parameter( PassIn             = 1 )
	parameter( PassOut            = 2 )
	parameter( PassNoWait         = 4 )
	parameter( PassByte  = 1 )
	parameter( PassRaw   = 2 )
	parameter( PassInt   = 3 )
	parameter( PassLong  = 4 )
	parameter( PassFloat = 5 )
	parameter( PassString= 6 )
c
	integer IVASsHDRmark,IVASsHDRopcode,IVASsHDRshift8
	integer IVASsHDRshiftRd,IVASsHDRbadMap,IVASsHDRbadXaddr
	integer IVASsHDRbadYaddr,IVASsHDRretOver,IVASsHDRbadCount
c
	parameter( IVASsHDRmark       = 2 )
	parameter( IVASsHDRopcode     = 258 )
	parameter( IVASsHDRshift8     = 514 )
	parameter( IVASsHDRshiftRd    = 770 )
	parameter( IVASsHDRbadMap     = 1026 )
	parameter( IVASsHDRbadXaddr   = 1282 )
	parameter( IVASsHDRbadYaddr   = 1538 )
	parameter( IVASsHDRretOver    = 1794 )
	parameter( IVASsHDRbadCount   = 2050 )
c
c  Status Codes:
c    IVAssDMAerror	DMA (1)
c    IVASsINTcomplete	Intrinsics (2)
c    IVASsURTerror	UART (3)
c    IVASsSIGactOver	SIGNAL (4)
c    IVASsGPHerror	Graphics (5)
c    IVASsDAMnoMem	Dyn. Mem. (6)
c
	integer IVASsDMAerror,IVASsINTcomplete,IVASsURTerror
	integer IVASsSIGactOver,IVASsGPHerror,IVASsDAMnoMem
	integer IVASsDAMnoCtxt
c
	parameter( IVASsDMAerror      = 6 )
	parameter( IVASsINTcomplete   = 8 )
	parameter( IVASsURTerror      = 14 )
	parameter( IVASsSIGactOver    = 18 )
	parameter( IVASsGPHerror      = 22 )
	parameter( IVASsDAMnoMem      = 28 )
	parameter( IVASsDAMnoCtxt     = 284 )
c
c  Memory management constants.
c
	integer IVASmaxTransfer,IVASmemorySize,IVASeightBitRm
	integer IVAStwelveBitRM,IVASmaxRM,IVASavailRMmask
c
	parameter( IVASmaxTransfer = 32768 )
	parameter( IVASmemorySize  = 1024 )
	parameter( IVASeightBitRM  = 3 )
	parameter( IVAStwelveBitRM = 2 )
	parameter( IVASmaxRM       = IVASeightBitRM )
	parameter( IVASavailRMmask = 7 )
