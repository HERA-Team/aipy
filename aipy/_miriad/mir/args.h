	integer MAXLEN,MAXARG
	parameter(MAXLEN=1024,MAXARG=128)
	logical init
	integer nargs,range(2,MAXARG)
	character line*(MAXLEN)
	common/argscom/init,nargs,range,line
	data init/.false./
