	integer MAXPARMS
	parameter(MAXPARMS=256)
	integer nparms
	character parname(MAXPARMS)*24
	double precision parvalue(MAXPARMS)
	common/obsparc/parname
	common/obsparv/parvalue,nparms
c
c  The Convex compiler seems to (incorrectly) flag as an error save statements
c  for variables that are in common blocks.
c	save parname,parvalue,nparms

