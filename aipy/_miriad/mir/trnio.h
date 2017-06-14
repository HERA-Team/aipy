	include 'maxdim.h'
	include 'maxnax.h'
	integer MAXSLOTS
	parameter(MAXSLOTS=6)
	integer size(MAXNAX,MAXSLOTS),lScr(MAXSLOTS),ndim(MAXSLOTS)
	integer memsize(MAXSLOTS),blk(MAXSLOTS),buf(MAXSLOTS)
	integer p(MAXSLOTS)
	logical inuse(MAXSLOTS),inmem(MAXSLOTS),flip(MAXNAX,MAXSLOTS)
	logical pre(MAXSLOTS),major(MAXSLOTS),post(MAXSLOTS)
	common/trancom/lScr,size,ndim,memsize,blk,p,buf,inuse,
     *	  inmem,pre,major,post,flip
c
	real ref(MAXBUF)
	common ref
