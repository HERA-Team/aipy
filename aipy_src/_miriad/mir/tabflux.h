      integer NTABLE
      parameter (NTABLE=5000)
c
      character TSOURCE(NTABLE)*40
      character ROOT(NTABLE)*40
      character ALIAS(NTABLE)*40
      integer NTAB, TINDEX(NTABLE)
      integer NALIASES
      real TFREQ(NTABLE), TFLUX(NTABLE), TRMS(NTABLE)
      double precision TDATE(NTABLE)
c
      common / TCOMI /NTAB, TINDEX, NALIASES
      common / TCOMR /TFREQ, TFLUX, TRMS
      common / TCOMD /TDATE
      common / TCOMS /TSOURCE, ROOT, ALIAS
