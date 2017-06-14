C  [calpass.h]
C
C	History:
C
C	             Original version
C	02aug91 mjs  Changed to use CALIB-specific MAXBASE/MAXANT
C	19apr93	pjt	assume caldefs.h has been included
C
C   Common and declaration statements for CALPASS program
C
C   Space used: 4*(MAXCHAN+MAXWIN)*MAXBASE ~ 4*MAXCHAN*MAXBASE
C		
        REAL      scalamp(MAXBASHC,MAXCHAN), scalphas(MAXBASHC,MAXCHAN),
     1            freqs(MAXBASHC,MAXCHAN)
        INTEGER   starwin(MAXWIN), chanwin(MAXWIN),baseline,b,nwins
        CHARACTER code(MAXBASHC,MAXWIN)*4
        LOGICAL   passflag(MAXBASHC,MAXCHAN)

        COMMON /passbandr/ scalamp, scalphas, freqs
        COMMON /passbandi/ starwin, chanwin,baseline,b,nwins
        COMMON /passbandc/ code
        COMMON /passbandl/ passflag

