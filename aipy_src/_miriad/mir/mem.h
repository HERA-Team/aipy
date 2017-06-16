	real MemR(MAXBUF)
	integer MemI(MAXBUF)
	logical MemL(MAXBUF)
	double precision MemD(MAXBUF/2)
	complex MemC(MAXBUF/2)
	equivalence(MemR,MemI,MemL,MemD,MemC)
	common MemR
