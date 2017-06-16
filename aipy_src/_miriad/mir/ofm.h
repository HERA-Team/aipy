c-----------------------------------------------------------------------
c Include file for ofm.for     ofm -> "output function memory"
c
c iofm    is the number of the desired ofm
c         1->b&w, 2->rainbow, 3->Tody, 
c         4->colour contours, 5->fixed zero colour conturs, 6->rgb, 
c         7->background, 8->heat
c ci1,ci2 are the min and max colour indices available to PGGRAY
c
c ofmb    is the BASIC 256 level ofm
c na      is the number of colour indices in the current ACTIVE table
c ofma    is the ACTIVE ofm 
c ofms    is  copy of the ACTIVE ofm before any fiddling was done
c savact  is a scratch array for copying the ACTIVE table around
c savsav  is a scratch array for copying the SAVE table around
c fid     is the last transfer function fiddle index array
c
c ofmdun  says that an ofm has been applied to the PGPLOT device
c fidun   says that a fiddle has been done to OFMA
c hedun   says that histogram equalization has been done this
c         call to ofmmod
c nocurs  says that device does not have a cursor
c ci0     is the colour index at the intensity=0 boundary for
c         fixed zero colour contours
c
c tfpic   true if viewport provided & we can draw transfer function plot
c tfvp    is the viewport to use for the transfer function plot
c xt,yt   are the transfer function plot arrays
c tflab   is the title to write on the transfer function plot
c tfcs    is the character size to use for the transfer function plot
c
c-----------------------------------------------------------------------
c
      integer maxlev
      parameter (maxlev=256)
c
      real ofmb(maxlev,3), ofma(maxlev,3), ofms(maxlev,3),
     + savact(maxlev,3), savsav(maxlev,3), xt(maxlev), yt(maxlev),
     + tfvp(4), tfcs
      integer na, iofm, ci1, ci2, fid(maxlev), ci0
      character*4 tflab
      logical fidun, tfpic, hedun, ofmdun, nocurs
c
      common /ofm1/ ofmb, ofma, ofms, savact, savsav
      common /ofm2/ fid, iofm, na, ci1, ci2, ci0
      common /ofm3/ fidun, tfpic, hedun, ofmdun, nocurs
      common /ofm4/ xt, yt, tfvp, tfcs
      common /ofm5/ tflab
