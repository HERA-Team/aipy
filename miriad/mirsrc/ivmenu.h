	integer maxfield
	parameter(maxfield=32)
	integer nfields,coord(4,maxfield)
	logical select(maxfield)
	character tags(maxfield)*16
	common/menucom/nfields,coord,select
	common/menucomc/tags
