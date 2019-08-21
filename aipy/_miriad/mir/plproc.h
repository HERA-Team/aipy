c  The following are invariant after plSet.
c   planet	Planet number.
c   meandist	Standard distance of the planet from Earth, in km.
c   erad,prad	Equatorial and polar radius of the planet, in radians
c   fake	True if we are faking the ephemeris.
c
c  The following vary with time.
c   tprev	Time (Julian date, in UT1).
c   ldash,mdash The unit vectors, in a planetocentric coordinate system,
c		of the Earth-equatorial l- and m-axes.
c   bpa		Position angle of the planet.
c   bmaj,bmin	Major and minor axes, in radians, at the standard distance.
c   cospa,sinpa	Cosine and sine of the inclination of the axis.
c   fac		Scale factor to convert to the standard distance.
c   nmat	Number of allocated shadow matrices.
c   dist	Distance (in km) of the planet.
c
c  The following are stored for later use:
c   smat	Shadowing matrices.
c  
	integer MAXMAT
	parameter(MAXMAT=1000000)
	double precision ldash(3),mdash(3),n0(3)
	double precision meandist,dist,tprev,smat(3,3,MAXMAT)
	real bmaj,bmin,bpa,fac,cospa,sinpa,erad,prad
	integer planet,nmat
	logical fake
c
	common/plprocom/ldash,mdash,n0,meandist,dist,tprev,smat,
     *	  erad,prad,bmaj,bmin,bpa,fac,cospa,sinpa,planet,nmat,fake
