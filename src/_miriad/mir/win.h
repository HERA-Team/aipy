c=======================================================================
c*win_h -- win data commons
c&pjt 
c:graphics, windows, include
c+
c  This set of variables keeps track of a multi-panel of 
c  WinMaxX by WinMaxY windows,
c
c	In the win subroutines, any variable of dimension *4
c	implies that the values are { X-lo, X-hi, Y-lo, Y-hi },
c       in that order. Examples are the arrays win0, win1
c
c

c
c  NXMAX -- the maximum number of horizontal windows
c
      INTEGER NXMAX
      PARAMETER( NXMAX = 40 )
c
c  NYMAX -- the maximum number of vertical windows
c
      INTEGER NYMAX
      PARAMETER( NYMAX = 30 )
c
c  win0 -- the screen coordinates of the windows
c
      REAL win0( NXMAX, NYMAX, 4 )
c
c  win1 -- the user coordinates inside the windows
c
      REAL win1( NXMAX, NYMAX, 4 )
c
c  WinN*** -- the set of windows currently active
c
      INTEGER WinNXlo, WinNXhi, WinNYlo, WinNYhi
c
c  WinM* -- the currently used/available windows
c
      INTEGER WinMaxX, WinMaxY
c
c  /win/ -- common block with all the above data
c
      COMMON /win/ WinNXlo, WinNXhi, WinNYlo, WinNYhi, 
     *             win0, win1,
     *             WinMaxX, WinMaxY
c--
c=======================================================================
