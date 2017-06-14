	integer TIME,ANTS,UV,POINT,AMP,VISNO,WINDOW,UVN,OR,DRA,DDEC
	integer INC,POL,ON,DAYTIME,SHADOW,AUTO,FREQ,SOURCE,RA,DEC
        integer BIN,HA,LST,ELEV,DAZIM,DELEV
	integer ITYPE,NSIZE,LOVAL,HIVAL,MDVAL,NTYPES
	parameter(TIME=1,ANTS=2,UV=3,POINT=4,AMP=5,VISNO=6,WINDOW=7)
	parameter(UVN=8,OR=9,DRA=10,DDEC=11,INC=12,POL=13,ON=14)
	parameter(SHADOW=15,AUTO=16,FREQ=17,SOURCE=18)
	parameter(RA=19,DEC=20,BIN=21,HA=22,LST=23,ELEV=24)
        parameter(DAZIM=25,DELEV=26)
	parameter(NTYPES=26,DAYTIME=NTYPES+1)
	parameter(ITYPE=0,NSIZE=1,LOVAL=2,HIVAL=3,MDVAL=4)
	character types(NTYPES)*12
	data types(TIME)   /'time        '/
	data types(ANTS)   /'antennae    '/
	data types(UV)     /'uvrange     '/
	data types(UVN)    /'uvnrange    '/
	data types(POINT)  /'pointing    '/
	data types(AMP)    /'amplitude   '/
	data types(VISNO)  /'visibility  '/
	data types(WINDOW) /'window      '/
	data types(OR)     /'or          '/
	data types(DRA)    /'dra         '/
	data types(DDEC)   /'ddec        '/
	data types(INC)    /'increment   '/
	data types(POL)	   /'polarization'/
	data types(ON)	   /'on          '/
	data types(SHADOW) /'shadow      '/
	data types(AUTO)   /'auto        '/
	data types(FREQ)   /'frequency   '/
	data types(SOURCE) /'source      '/
	data types(RA)     /'ra          '/
	data types(DEC)    /'dec         '/
	data types(BIN)	   /'bin         '/
	data types(HA)     /'ha          '/
	data types(LST)    /'lst         '/
	data types(ELEV)   /'elevation   '/
	data types(DAZIM)  /'dazim       '/
	data types(DELEV)  /'delev       '/

