// jhz - may2004 */
// 2004-5    propose the development of SMA code unde miriad */
// 2004-5    start to write the SMAmir2Miriad code (SMALOD) */
// 2004-7-15 read/write out all vis data */
// 2004-7-16 write lst */
// 2004-7-17 correct for spectral order based on s sequence */
// 2004-7-19 source information */
// 2004-7-20 sideband separation */
// 2004-7-22 fix coordinates for j2000 and apparent (observing) */
// 2004-7-23 fix skip beginning and ending scans */
// 2004-7-24 check and correct the flagging */
// 2004-7-25 add a variable sourceid */
// 2004-7-26 Tsys EL AZ information */
// 2004-8-2  Instrument, telescope, observer version */
// 2004-8-4  flip the phase for lsb */
// 2004-11-30 added a function of retrieving daul rx/if.
// 2004-12-10 rename header file for miriad 4.0.4 - pjt 
// 2004-12-16 increased length of pathname from 36 to 64
//            increased length of location from 64 to 81
// 2005-01-11 merged sma_Resample.c to sma_mirRead.c
//            sma_mirRead.c handle both original correlator
//            configuration and resample the spectra to
//            lower and uniform resolution.
// 2005-02-15 decoded veldop from mir file and fix the
//            the reference channel for miriad.
// 2005-02-22 decoded Tsys from mir antenna based file.
// 2005-02-28 read dual receiver data.
// 2005-03-01 added a function to fix the spectral chunk order
//            which is a temporal problem in the mir file
// 2005-03-02 added an option to read engineer file.
// 2005-03-07 changed the chunk order in frequnecy for the first
//            three blocks in the data observed during the
//            690 campaign spring 2005
// 2005-03-07 added the polarization options circular for waveplates
//            and default for linear
// 2005-03-08 changed the options name cirpol to circular
//            to match the same options in miriad program elsewhere.
// 2005-03-08 decoded antenna positions from mir baseline coordinates
//            converted the geocentrical coordinates to Miriad
//            coordinates system. 
// 2005-03-10 made the consistent array lengths for the variables 
//            required by miriad programs such as bee, uvlist etc.
// 2005-03-17 added linear in options for linear polarization;
//            made the default in options for polarization to be
//            nopol.
// 2005-03-18 decoded velocity type from mir.
// 2005-03-18 fixed problems of el and az in uvput
// 2005-03-21 added integration skip in the process of
//            decoding antenna position from baseline vectors.
// 2005-03-23 fixed a bug in decoding antenna positions for
//            the antenna id > reference antenna's id.
// 2005-03-23 fixed a bug in baseline pntr (last bl of
//            the 1st rx overlapping with that of 2nd).
// 2005-03-29 fixed problem in decoding antenna coordinates
//            in the case missing antenna in a random place.
// 2005-03-31 trim the junk tail in source name.
//            the source name truncates to 8 char.
// 2005-04-05 fixed the polarization label conversion for circular case
//            fixed the uv coordinates scaling (by the base frequency fsky[0]);
// 2005-04-28 change the phase conjugate scheme according Taco's log.
// 2005-05-05 add options of computeing radial velocity wrt either barycenter
//            or lsr.
// 2005-05-11 read spectral configuration from the integration
//            that users want to start (nscans).
// 2005-05-23 (PJT) added prototypes from miriad.h (via sma_data.h) 
//             and cleaned up some code because of this, sans indent
// 2005-06-01 (JHZ) fixed a few bugs in incompatible pointer type.
//                  stat in slaCldj is initialized and added  
//                  stat parse after return.
// 2005-06-02 (JHZ) implemented vsource and restfreq from users
//                  inputs.
// 2005-06-08 (JHZ) implemented a feature of handling
//                  arbitary frequency configuration 
//                  for each of the receivers in the case
//                  of dual recievers.
// 2005-06-20 (JHZ) fixed a bug (pointing to a wrong component) in calculate 
//                  site velocity.
// 2005-06-21 (JHZ) fixed a bug  in the status handle of rspokeflshsma_c
// 2005-06-22 (JHZ) fixed all the loose ends (warnings from compilers)
// 2005-06-22 (JHZ) add a feature allowing user' input of restfrequency
// 2005-06-23 (JHZ) add ut var
// 2005-06-27 (JHZ) add initializing blarray
// 2005-07-05 (JHZ) remove a hidden phase flip for the lsb data
// 2005-07-07 (JHZ) updated aperture efficiency and jyperk
// 2005-07-07 (JHZ) add  parsing source name and changing the source
//                  name if the first 8 chars are identical in
//                  any two source name entries from mir data.
// 2005-08-03 (JHZ) fixed a bug in the channel pntr in the
//                  case of resampling the data to a lower channel resolution.
// 2005-08-03 (JHZ) change the apperture efficiency to 0.5 for 340 GHz
// 2005-09-01 (JHZ) add smaveldop in the structure smlodd. 
//                  see the comment around line 1610. the chunk frequency
//                  stored in the header of MIR data appears to be not the
//                  sky frequency that is actually Doppler tracked by the 
//                  on-line system. The chunk frequency is the sky frequency to 
//                  which the part of Doppler velocity (diurnal term and part 
//                  of the annual term) has been corrected.
//                  smaveldop is the residual veldop after taking out the part
//                  that have been corrected to the chunk frequency.
//                  smaveldop will be specially stored in a variable
//                  called velsma. A special patch in the uvredo (miriad task)
//                  is required to compute the "SMA veldop" for other sources
//                  corresponding to the SMA "sky frequency".
//                  So that, the users can use the rest of Miriad task
//                  to properly reduce the SMA data avoiding the problem 
//                  of smearing spectral lines.  
// 2005-09-01 (JHZ) Removed lines for storing velsma;
//                  The residual Doppler velocity can be stored in veldop. 
// 2005-09-25 (JHZ) removed unused variables.
// 2005-09-29 (JHZ) added vsource (from users input) back to veldop.
// 2005-10-10 (JHZ) added online flagging pntr (convert negative wt to
//                  to online-flagging state of 0).
// 2005-10-11 (JHZ) added option of reading antenna positions from
//                  the ASCII file 'antennas'. 
// 2005-10-13 (JHZ) added swapping polarization states (RR<->LL, LR<->RL)
//                  for the circular polarization data taken before
//                  2005-06-10 or julian date 2453531.5
//                  (asked by dan marrone).
// 2005-10-13 (JHZ) add skipping decoding Doppler velocity
//                  for the polarization data old than 2005-06-10
//                  because the velocity entry in the Mir header
//                  was screwed up.
// 2005-10-13 (JHZ) add options of noskip
// 2005-11-30 (JHZ) store J2000 coordinates of the true
//                  pointing position for all the cases other than stored
//                  the true pointing position only
//                  when the true pointing position differs from
//                  the J2000 coordinates of source catalog
//                  position.
// 2005-12-5  (JHZ) obsoleted options=oldpol;
//                  unified the pol-state conversion using
//                  the function ipolmap.
// 2005-12-6  (JHZ) implemented a feature skipping spectral
//                  windows (spskip).
// 2005-12-30 (JHZ) in the new MIR data sets (230/690
//                  observation 051214_05:48:56 for example),
//                  the the number of inregrations (inhid)
//                  in the integration header (in_read) differs from 
//                  that in the baseline header (bl_read).
//                  The baseline header stores additional
//                  integration in the end of the file; this
//                  bug (spotted for the first time)
//                  causes the segmentation fault in smalod.
//                  The problem has been fixed by adding
//                  a parsing sentence for selecting
//                  the smaller  number of total integrations
//                  from the header files (in_read and bl_read).
// 2006-01-9 (JHZ)  The LONG_MAX used in function mfsize
//                  is a system variable depending on the WORDSIZE
//                  of the computers. The function mfsize worked
//                  only for LONG_MAX=2147483647L, 
//                  the value for 32bits computer.
//                  The LONG_MAX is replaced by SMA_LONG_MAX which
//                  is defined to 2147483647L. 
// 2006-01-18 (JHZ) add instruction print out in case
//                  the option rxif is not properly set.
// 2006-01-24 (JHZ) add the options spskip=-1 to take the frequency
//                  configuration of the first integration and
//                  skip the rest of the frequency configuration.
//                  This assumes that the frequency configuration
//                  does not change during the observing run.
//                  This is for loading the old SMA data. 
// 2006-02-03 (JHZ) optimized the memory requirements.
//                  for wts structure and double sideband loading.
// 2006-02-09 (JHZ) added a feature to handle multiple correlator configuration.
// 2006-03-13 (JHZ) added a feature to handle 2003 incompleted-correlator
//                  data.
// 2006-04-07 (JHZ) added an instruction message for properly using nscan.
// 2006-05-19 (JHZ) implemented feature to handle hybrid high spectral
//                  resolution frequency configuration, allowing
//                  presence of empty chunks.
// 2006-06-08 (JHZ) fixed a bug in parsing the source id in the case
//                  no source information given in the mir data
//                  (an on-line bug).
// 2006-06-09 (JHZ) fixed two warning bugs seen from 32bits:
//                  removed zero-length printf and added parathesis
//                  for if-else around line 2534.
//***********************************************************
#include <math.h>
#include <rpc/rpc.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <fcntl.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "sma_data.h"
#define OK         0
// define SMA_LONG_MAX 
# define SMA_LONG_MAX  2147483647L


/* extern variable while read mir data */
char pathname[36];
static FILE* fpin[6];
int nsets[6];
double jday; /* julian day */
struct inh_def   **inh;
struct blh_def   **blh;
/*struct sph_def   **sph;*/
struct codeh_def **cdh;
struct ant_def   **enh;
struct sch_def   **sch;


/* declare the externals struct which is the same to
   the common blocks, atlodd and atlodc */

char sname[64];
smlodd smabuffer;
// initialize 

struct vis { float real;
  float imag;
};

unsigned long mfsize(FILE *);
struct inh_int_def {
  int  a[512];
};


/* prototypes of everything used here */

void rsmirread_c(char *datapath, char *jst[]);
void rsmiriadwrite_c(char *datapath, char *jst[]);
void rssmaflush_c(int scanskip, int scanproc, int sb, int rxif, int dosporder, int doeng, int doflppha);
void rspokeinisma_c(char *kst[], int tno1, int *dosam1, int *doxyp1, int *doop1, int *dohann1, 
int *birdie1, int *dowt1, int *dopmps1, int *dobary1, int *doif1, int *hires1, int *nopol1, 
int *circular1, int *linear1, int *oldpol1, double lat1, double long1, int rsnchan1, 
int refant1, int *dolsr1, double rfreq1, float *vsour1, double *antpos1, int readant1, 
int *noskip1, int *spskip1, int *dsb1, int*mcconfig1, int*nohighspr);
void rspokeflshsma_c(char *kst[]);


int rsgetdata(float smavis[2*7681], int smaflags[7681], int *smanchan, int p, int bl, int sb, int rx);
struct pols *rscntstokes(int npol, int bl, int sb, int rx);
int rsmir_Read(char *datapath, int jstat);
struct inh_def *inh_read(FILE *fpinh);
struct blh_def *blh_read(FILE *fpblh);
unsigned long mfsize(FILE *fp);
struct sph_def *sph_read(FILE *fpsph);
struct codeh_def *cdh_read(FILE *fpcodeh);
struct ant_def *enh_read(FILE *fpeng);
struct sch_def *sch_head_read(FILE *fpsch);
int sch_data_read(FILE *fpsch, long int datalength, short int *data);
char *rar2c(double ra);
char *decr2c(double dec);
int spdecode(struct codeh_def *specCode[]);
float juliandate(struct codeh_def *refdate[]);
double slaCldj(int iy, int im, int id, int sj);
void precess(double jday1, double ra1, double dec1, double jday2, double *ra2pt, double *dec2pt);
void nutate(double jday, double rmean, double dmean, double *rtrueptr, double *dtrueptr);
void nuts(double jday, double *dpsiptr, double *depsptr);
double mobliq(double jday);
void aberrate(double jday, double ra, double dec, double *rappptr, double *dappptr);
void vearth(double jday, double pos[3], double vel[3]);
void elaz(int tno);
void tsysStore(int tno);
double velrad(short dolsr, double time, double raapp, double decapp, double raepo, double decepo, double lst, double lat);
struct lmn *sph2lmn(double ra, double dec);
struct vel *vsite(double phi, double st);
void vsun(double *VEL);
short ipolmap(short input_ipol);

/* interface between fortran and c */
void rsmirread_c(char *datapath, char *jst[])
{ 
  int jstat;
  strcpy(pathname,datapath);
  jstat= (int)*jst;
  jstat = rsmir_Read(pathname,jstat);
  *jst = (char *)jstat; 
}

void rsmiriadwrite_c(char *datapath, char *jst[])
{ 
  int jstat;
  jstat=-1;
  /* open mir files */
  jstat = rsmir_Read(pathname,jstat);
  *jst = (char *)jstat;
  /* then start to read and write data if the files are ok. */
  if (jstat==0) {
    jstat=0;
    jstat = rsmir_Read(pathname,jstat);
  } else {
    printf("file problem\n");
  }
}

void rssmaflush_c(int scanskip,int scanproc,int sb,int rxif,int dosporder,int doeng,int doflppha)
{ /* flush mirdata*/
  char *kst[4];
  int kstat=-1;  
  int tno;
  char telescope[4];
  char instrument[4];
  char observer[16];
  char version[16];

  tno = smabuffer.tno;
  smabuffer.scanskip=scanskip;
  smabuffer.scanproc=scanproc;
  smabuffer.sb = sb;
  smabuffer.rxif= rxif;
  smabuffer.doChunkOrder = dosporder;
  smabuffer.doeng = doeng;
  smabuffer.doConjugate = doflppha;
  *kst = (char *)&kstat;
  /*  read header  */
    rspokeflshsma_c(kst);  
  /*  write ante numbers */
  if(smabuffer.nants!=0) {
    //         uvputvri_c(tno,"nants",&(smabuffer.nants),1);
    /*  write telescope name and other description parameters */
    sprintf(telescope,"SMA");
    sprintf(instrument,"SMA");
    sprintf(observer,"SmaUser");
    sprintf(version, "test");
    uvputvra_c(tno, "telescop", telescope);
    uvputvra_c(tno, "instrume", instrument);
    uvputvra_c(tno, "observer", observer);
    uvputvra_c(tno, "version", version);
  }
}

void rspokeinisma_c(char *kst[], int tno1, int *dosam1, int *doxyp1,
		    int *doop1, int *dohann1, int *birdie1, int *dowt1, int *dopmps1,
		    int *dobary1, int *doif1, int *hires1, int *nopol1, int *circular1,
		    int *linear1, int *oldpol1, double lat1, double long1, int rsnchan1, 
		    int refant1, int *dolsr1, double rfreq1, float *vsour1,
		    double *antpos1, int readant1, int *noskip1, int *spskip1,
	            int *dsb1, int *mcconfig1, int *nohighspr1)
{ 
  /* rspokeflshsma_c == pokeflsh */
  int buffer, i;
  /* initialize the external buffers */   
  strcpy(sname, " ");
  smabuffer.tno    = tno1;
  smabuffer.rsnchan= rsnchan1;
  smabuffer.dosam  = *dosam1;
  smabuffer.doxyp  = *doxyp1;
  smabuffer.opcorr = *doop1;
  smabuffer.dohann = *dohann1;
  smabuffer.doif   = *doif1;
  smabuffer.dobary = *dobary1;
  smabuffer.birdie = *birdie1;
  smabuffer.dowt   = *dowt1;
  smabuffer.dopmps = *dopmps1;
  smabuffer.hires  = *hires1;
  smabuffer.nopol  = *nopol1;
  smabuffer.circular = *circular1;
  smabuffer.linear = *linear1;
  smabuffer.oldpol = *oldpol1;
  smabuffer.lat    = lat1;
  smabuffer.longi  = long1;
  smabuffer.refant = refant1;
  smabuffer.dolsr  = *dolsr1;
  smabuffer.noskip = *noskip1;
  smabuffer.vsource= *vsour1;
  smabuffer.dsb    = *dsb1;
  smabuffer.juldate= -10.00;
  smabuffer.spskip[0] = spskip1[0];
  smabuffer.spskip[1] = spskip1[1];
  smabuffer.mcconfig  = *mcconfig1;
  smabuffer.highrspectra = *nohighspr1;
      printf("User's input vSource = %f km/s\n", smabuffer.vsource);
      if(rfreq1 > 0.00001 || rfreq1 < -0.00001) {
             for (i=0; i<SMIF+1; i++) {
             smabuffer.restfreq[i] = rfreq1;
                                      }
            smabuffer.dorfreq = -1;
      printf("User's input RestFrequency = %f GHz\n", rfreq1);
                    } else {
             smabuffer.dorfreq =  1;
                                      }
             smabuffer.antpos[0] = readant1;
             smabuffer.readant   = readant1;
             for (i=1; i<readant1*3+1; i++) {
             smabuffer.antpos[i]=antpos1[i-1];
                    }
  if(smabuffer.dowt>0) {
    /* call lagwt(wts,2*smcont-2,0.04) */
    /* process weights here. */ 
  }
  smabuffer.newsc = FALSE;
  smabuffer.newfreq = FALSE;
  smabuffer.nants = 0;
  smabuffer.nifs = 0;
  smabuffer.nused = 0;
  smabuffer.tcorr = 0;
  buffer=(int)*kst;
  *kst= OK;
}

void rspokeflshsma_c(char *kst[])
{ /* rspokeflshsma_c== pokeflsh */
  int tno;
  int i1, i2, ifs, p, bl, sb, rx, nchan, nspect;
  int npol,ipnt,ischan[SMIF];
  int tbinhi,ibuff;
  double preamble[5], tdash;
  long int dummy;
  float jyperk, eta, eta_c, eta_a, r_ant=3, pi;
  float vis[2*MAXCHAN];
  int flags[MAXCHAN]; 
  struct pols *polcnt;
  char telescope[4];
  char instrument[4];
  char observer[16];
  char version[16];
// initialization:
  eta_a=0.75;
  tno = smabuffer.tno;
  sb = smabuffer.sb; /* sb=0 for lsb; sb=1 for usb; sb=2 for both */
  rx = smabuffer.rxif;
  if(smabuffer.nused==0) 
    return;
  /* put ants to uvdata */
  if(smabuffer.nants!=0) {
    uvputvri_c(tno,"nants",&(smabuffer.nants),1);
  /*  write telescope name and other description parameters */
    sprintf(telescope,"SMA");
    sprintf(instrument,"SMA");
    sprintf(observer,"SmaUser");
    sprintf(version, "test");
    uvputvra_c(tno, "telescop", telescope);
    uvputvra_c(tno, "instrume", instrument);
    uvputvra_c(tno, "observer", observer);
    uvputvra_c(tno, "version", version);
  }

  if(smabuffer.newfreq>0) {
    if(smabuffer.doif>0) {
      for (ifs=1; ifs < smabuffer.nifs; ifs++) {
	if(smabuffer.nstoke[ifs-1]!=smabuffer.nstoke[0]) 
	  bug_c( 'f', "Number of polarisations differ between IFs. Use options=noif.\n"); 
	for (p=1; p< smabuffer.nstoke[ifs-1]; p++) {
	  if(smabuffer.polcode[ifs-1][p-1][0]!=smabuffer.polcode[0][p-1][0]) 
	    bug_c( 'f', "Polarisation types differ between IFs. Use options=noif.\n");
	}
      }
    }
  } else {
    if(smabuffer.hires > 0) 
      for (ifs=1; ifs<smabuffer.nifs; ifs++){
	if (smabuffer.nbin[ifs]!=smabuffer.nbin[0]) 
	  bug_c( 'f', 
		 "Number of bins in different IFs must agree for options=hires\n");
      }
  } 
  tdash  = smabuffer.time;
  tbinhi = 1;
// convert julian date to ut and store ut 
// ut and julian date on 2000 julian2000=2451544.5
// julain date = 2451544.5 + day of the yr + fraction of day from 0h UT 
// determine julian date
     if(smabuffer.juldate < -1.) { 
     dummy =  (long int) (smabuffer.time);
     smabuffer.juldate = (double) dummy + 0.5;
        }
// convert juldate to ut in radians
     smabuffer.ut = smabuffer.time - smabuffer.juldate;
     smabuffer.ut = smabuffer.ut*DPI*2.0;
// store ut
     uvputvrd_c(tno,"ut",&(smabuffer.ut),1);
// store apparent LST 
     uvputvrd_c(tno,"lst",&(smabuffer.lst),1);
// store elaz data 
     elaz(tno);
// store radial velocity of the observatory w.r.t. the rest frame 
  uvputvrr_c(tno,"veldop",&(smabuffer.veldop),1);
// store the source velocity w.r.t. the rest frame
  uvputvrr_c(tno,"vsource",&(smabuffer.vsource),1);
//
// antenna aperture efficiency from TK
// nearly no elevation dependence at EL > 20 deg
// at 690 GHz there have been no proper measurements yet
// 0.4 would be number assumed.  07-jul-05
//
  switch(smabuffer.rxif) {  
  case 0: eta_a=0.75;     /* jyperk=139. assuming eta_a=0.7 d=6m */
          break;
  case 1: eta_a=0.5;     /* jyperk=194. assuming eta_a=0.5 d=6m */   
          break;
  case 2: eta_a=0.4;     /* jyperk=242. assuming eta_a=0.4 d=6m */  
          break;
             }
   eta_c = 0.88;
   eta   = eta_a*eta_c;
   pi    = (float) DPI;
// calculate jyperk
   jyperk=2.* 1.38e3/pi/(eta*r_ant*r_ant);
   uvputvrr_c(tno,"jyperk",&jyperk,1);
// Handle the case that we are writing the multiple 
// IFs out as multiple records. 
  if(smabuffer.doif!=1&&smabuffer.nifs>1) {
    nspect =ischan[0]= 1;
    for(ifs=0; ifs < smabuffer.nifs; ifs++) {
      printf(" writing headers\n");
      uvputvri_c(tno,"nspect",&nspect,1);
      printf(" search\n");
      uvputvri_c(tno,"npol",  &(smabuffer.nstoke[ifs]),1);
      uvputvri_c(tno,"nschan",&(smabuffer.nfreq[ifs]),1);
      uvputvri_c(tno,"ischan",&ischan[0],1);
      uvputvrd_c(tno,"sfreq", &(smabuffer.sfreq[ifs]),1);
      uvputvrd_c(tno,"sdf",  &(smabuffer.sdf[ifs]),  1);
      uvputvrd_c(tno,"restfreq",&(smabuffer.restfreq[ifs]),1);
      bl=0;
      for(i2=1; i2<smabuffer.nants+1; i2++){
	for(i1=1; i1<i2+1; i2++){
	  preamble[0] = smabuffer.u[bl];
	  preamble[1] = smabuffer.v[bl];
	  preamble[2] = smabuffer.w[bl];
	  preamble[3] = tdash;
	  preamble[4] = 256*i1 + i2;
	  for(p=0; p<smabuffer.nstoke[ifs]; p++){
          ipnt = smabuffer.pnt[ifs][p][bl][0][0];
          printf("ipnt=%d\n", ipnt);
	  if(ipnt>0) 
	  uvputvrr_c(tno,"inttime",&smabuffer.inttime[bl],1);
	  if(smabuffer.opcorr==0) {
/* call opapply(data(ipnt),nfreq(if),fac(if)) */
/* bug_c('f',"This code section had a bug"); */
/* next line had serious bug before may 23 */
// uvwrite_c(tno,preamble,&smabuffer.data[ipnt],flags,smabuffer.nfreq[ifs]); 
	    }
	  }
	  bl++;
	}
      }
    }
  } else {
// Handle the case were we are writing the multiple IFs out as a single record.
// This way is adopted to store SMA data in Miriad format.
    if(smabuffer.newfreq>0) {
      ischan[0] = 1;
// counting total number of spectral channels
      for (ifs = 1; ifs < smabuffer.nifs; ifs++) {
	ischan[ifs] = ischan[ifs-1] + smabuffer.nfreq[ifs];
//        printf("MMMMM %d", smabuffer.nfreq[ifs]);
      }
      uvputvri_c(tno,"nspect",&(smabuffer.nifs),1);
      uvputvri_c(tno,"ischan",&(ischan),smabuffer.nifs);
      uvputvri_c(tno,"nschan",&(smabuffer.nfreq),smabuffer.nifs);
      uvputvrd_c(tno,"sfreq",&(smabuffer.sfreq),smabuffer.nifs);
      uvputvrd_c(tno,"sdf",&(smabuffer.sdf),smabuffer.nifs);
      uvputvrd_c(tno,"restfreq",&(smabuffer.restfreq),smabuffer.nifs);
      uvputvrr_c(tno,"veldop",&(smabuffer.veldop),1);
    }
    uvputvri_c(tno,"tcorr",&(smabuffer.tcorr),1);
// store system temperature 
    tsysStore(tno); 
// store the random parameters and the visibility data
    bl=0;
    for(i2=1; i2<smabuffer.nants+1; i2++) {
      for(i1=1; i1<i2+1; i1++){
	preamble[0] = smabuffer.u[bl];
	preamble[1] = smabuffer.v[bl];
	preamble[2] = smabuffer.w[bl];
	preamble[3] = smabuffer.time; 
	preamble[4] = (double)smabuffer.blcode[bl]; 
	polcnt = rscntstokes(npol, bl, sb, rx);
	npol = polcnt->npol;
	if(npol>0) {
	uvputvri_c(tno,"npol",&npol,1);
	for(p=polcnt->polstart; p<polcnt->polend+1; p++) {
	    nchan = rsgetdata(vis,flags,&nchan, p, bl, sb, rx);    
	    if(nchan>0) {
	    ibuff = smabuffer.polcode[0][p][bl];
	    uvputvri_c(tno,"pol",&ibuff,1);
	    uvputvrr_c(tno,"inttime",&smabuffer.inttime[bl],1);
	    uvwrite_c(tno,preamble, vis,flags,nchan);         
	                }
	                                                  }
                        }               
	bl++;
            }
          }
        }
// reset the frequency processing handle
  smabuffer.newfreq=-1;
// re-initialize the pntr 
  for (ifs=0; ifs<SMIF; ifs++) {
    for (p=0; p<SMPOL; p++) {
      for (bl=0; bl<SMBAS; bl++) {
	for (sb=0; sb<SMSB; sb++) {
	  for (rx=0; rx<SMRX; rx++) {
	    smabuffer.pnt[ifs][p][bl][sb][rx]=0;
            smabuffer.flag[ifs][p][bl][sb][rx]=-1;
	  }
	}
      }
    }
  }
}

int rsgetdata(float smavis[2*MAXCHAN], int smaflags[MAXCHAN], int *smanchan, int p, int bl, int sb, int rx)
{ /* Construct a visibility record constructed from multiple IFs. */
  int nifs=smabuffer.nifs;
  float fac[nifs];
  int n,ipnt,i,nchand, nchan; 
  nchan = 0;
  nchand = 0;
  for (n=0; n<nifs; n++) {
    ipnt = smabuffer.pnt[n][p][bl][sb][rx]; 
//      printf("n ipnt %d %d\n", n, ipnt);
      if(ipnt>0) {
      if(nchan<nchand) {
	for (i=nchan; i<nchand; i++) {
	  smaflags[i] = -1;
	  smavis[2*i] = 0;
	  smavis[2*i+1] = 0;
	}    
	nchan = nchand;
      }
      for (i=nchan; i< nchan+smabuffer.nfreq[n]; i++){
	fac[n]=1000000.;

	smavis[2*i]   =  fac[n]*smabuffer.data[ipnt].real;
	smavis[2*i+1] =  fac[n]*smabuffer.data[ipnt].imag;
//	  (float)pow((double)(-1),(double)(sb+1)); 
        smaflags[i] =  smabuffer.flag[n][p][bl][sb][rx];       
	ipnt++;    
      }
      if(smabuffer.bchan[n]>=1&&smabuffer.bchan[n]<=smabuffer.nfreq[n])
	smaflags[nchan+smabuffer.bchan[n]] = smabuffer.flag[n][p][bl][sb][rx]; 
        nchan = nchan + smabuffer.nfreq[n];
    }
    nchand = nchand + smabuffer.nfreq[n];
  }
  
  if(nchan<nchand&&nchan>0) {
    for(i=nchan+1; i< nchand+1; i++) {
      smaflags[i] = smabuffer.flag[n][p][bl][sb][rx]; 
      smavis[2*i] = 0;
      smavis[2*i+1] = 0;
    }
    nchan = nchand;
  }
  *smanchan=nchan;
  return nchan;
}

struct pols *rscntstokes(int npol, int bl, int sb, int rx)
{ /*Determine the number of valid Stokes records in this record.*/
  int nifs = SMIF;
  int nstoke = SMPOL;
  short valid;
  int p, p1, p2, ifs;
  static struct pols polcnts;

  npol =0;
  p1=-1;
  p2=1;
  for (p=0; p< nstoke; p++) {
    valid = -1;
    for ( ifs =0; ifs<nifs; ifs++) {
      valid = valid;
      if(smabuffer.pnt[ifs][p][bl][sb][rx] > 0) valid = 1;
    }
    if(valid>0&&p1==-1) {p1=p; p2=p1;}
    if(valid>0) {npol++; if(p>p1) p2=p1+1;}
  }
  polcnts.npol=npol;
  polcnts.polstart=p1;
  polcnts.polend=p2;
  return &polcnts;
}          

int rsmir_Read(char *datapath, int jstat)
{
  char location[6][81];
  char pathname[64];
  char filename[6][36];
  char sours[9], smasours[33];
  int set, readSet;
  int file,nfiles = 6;
  int headerbytes[6];
  smEng **smaEngdata;
  int i,j,k,l,m,i0;
  int kk,iinset,lastinset,ireset,reset_id[10];
  long imax,bytepos,nbytes,datalength;
  long *data_start_pos;
  int firstbsl;
  int numberBaselines=0;
  int numberSpectra,numberSidebands,numberRxif;
  int blhid,firstsp,lastsp;
  int inhset,blhset,sphset,inhset1st;
  int spcode[25], frcode[25];
  short int scale;
  short int *shortdata;
  double r,cost,sint,z0,tmp, rar, decr;
  double antpos[3*MAXANT+1];
  int tno, ipnt, max_sourid;
  int kstat;
  char *kst[4];
  char target[6];
  char unknown[7];
  char skipped[8];
  int ntarget;
  time_t  startTime, endTime;
  float trueTime;
  blvector blarray[MAXANT][MAXANT];
  station  antenna[MAXANT];
  source   multisour[MAXSOURCE];
  struct xyz   antxyz[MAXANT];
  int sourceID, phaseSign;
  short oka,okb,okc,okd;
  correlator smaCorr;
  frequency  smaFreq[2];
  uvwPack **uvwbsln;
  visdataBlock  visSMAscan;
  int sphSizeBuffer;
  int nnants, flush;
  int ifpnt, polpnt, blpnt, sbpnt, rxpnt, sblpnt;
  int avenchan,miriadsp[SMIF+1];
  int blid_intchng[MAXINT];  
               /* the baseline id right after integration change */
  int rxlod;   /* the rx to load rxlod=0 -> smabuffer.rx1        
                                 rxlod=1 -> smabuffer.rx2        */
  float avereal, aveimag;
  extern struct inh_def   **inh;
  extern struct blh_def   **blh;
  static struct blh_config **bln;
  static struct sph_def   **sph;
  static struct sph_def   *sph1;
  static struct sph_config **spn;
  extern struct codeh_def **cdh;
  extern struct ant_def   **enh;
  extern struct sch_def   **sch;
  struct bltsys    **tsys;
  struct anttsys   **atsys;
  struct wtt    **wts;
// initialize
  startTime = time(NULL);
  phaseSign = 1;
  flush     = 1;
  polpnt    = 1;
  blpnt     = 1;
  rxlod     = 0;
  lastinset = 0;
  ireset    = 0;
//  
  strcpy(pathname,datapath);
  strcpy(filename[0],"in_read");
  strcpy(filename[1],"bl_read");
  strcpy(filename[2],"sp_read");
  strcpy(filename[3],"codes_read");
  strcpy(filename[4],"eng_read");
  strcpy(filename[5],"sch_read");
  strcpy(target,"target");
  strcpy(unknown,"unknown");
  strcpy(skipped,"skipped!");
  ntarget=0;
  
  /* number of bytes in each type of header. */
  /* sch is variable number of bytes */
  headerbytes[0] = 132;
  headerbytes[1] = 118;
  headerbytes[2] = 100;
  headerbytes[3] = 42;
  headerbytes[4] = 188;
  headerbytes[5] = 0;
  
  switch(jstat) {
  case -3:
    /* Open all the files */
    for (file=0;file<nfiles;file++){
      strcpy(location[file],pathname); 
      strcat(location[file],filename[file]);
      // open the mir files
      fpin[file] = fopen(location[file],"r");
      if (fpin[file] == NULL) {
	printf("Problem opening the file %s\n",location[file]);
	perror("file open problem");
	exit(-1);
      } else {
	printf("Found file %s\n",location[file]);
      }
    }
    // Get the size of each file and compute the number of headers 
    for (file=0;file<nfiles;file++){
      if(file < 5 ){
	imax = mfsize(fpin[file]);
	nsets[file] = imax / headerbytes[file];
      }
    }
    break;
  case 0:   /*read header & vis */
    startTime = time(NULL);
    /* Allocate memory to store all the headers */
    inh = (struct inh_def **) malloc(nsets[0]*sizeof( struct inh_def *));
    for (set=0;set<nsets[0];set++) {
      inh[set] = (struct inh_def *)malloc(sizeof(struct inh_def ));
      if (inh[set] == NULL ){
	printf("ERROR: Memory allocation for inh failed trying to allocate %d bytes\n", nsets[0]*sizeof(struct inh_def));
	exit(-1);
      }
    }
    if (smabuffer.scanproc==0||smabuffer.scanproc<0)
      smabuffer.scanproc = nsets[0] - smabuffer.scanskip-2;
    
    blh = (struct blh_def **) malloc(nsets[1]*sizeof( struct blh_def *));
    for (set=0;set<nsets[1];set++) {
      blh[set] = (struct blh_def *)malloc(sizeof(struct blh_def ));
      if (blh[set] == NULL ){
	printf("ERROR: Memory allocation for blh failed trying to allocate %d bytes\n", nsets[1]*sizeof(struct blh_def));
	exit(-1);
      }
    }
    // new baseline buffer used for configuring spectra.
    // the array length of bln is the same as the total number
    // of integration sets in the mir file.
    bln = (struct blh_config **) malloc(nsets[0]*sizeof( struct blh_config *));
    for (set=0;set<nsets[0];set++) {
      bln[set] = (struct blh_config *)malloc(sizeof(struct blh_config ));
      if (bln[set] == NULL ){
	printf("ERROR: Memory allocation for blh failed trying to allocate %d bytes\n", nsets[0]*sizeof(struct blh_def));
	exit(-1);
      }
    }
    // baseline coordinate buffer.
    uvwbsln = (struct uvwPack **) malloc(nsets[0]*sizeof( struct uvwPack));
    for (set=0;set<nsets[0];set++) {
      uvwbsln[set] = (struct uvwPack *)malloc(sizeof(struct uvwPack));
      if (uvwbsln[set] == NULL ){
	printf("ERROR: Memory allocation for uvwbsln failed for %d bytes\n",nsets[0]*sizeof(struct uvwPack));
	exit(-1);
      }
    }
    
    /* Read the headers */
    for (set=0;set<nsets[0];set++) {
      *inh[set] = *(inh_read(fpin[0])); 
      if (SWAP_ENDIAN) {
	inh[set] =  swap_inh(inh[set]);
	// reverse mapping the inhid wrt integration set number
	//          inid[inh[set]->inhid]=set;
      }
    }
    if (SWAP_ENDIAN) {
      printf("FINISHED READING  IN HEADERS (endian-swapped)\n");
    } else {
      printf("FINISHED READING  IN HEADERS\n");
    }
    for (set=0;set<nsets[1];set++) {
      *blh[set] = *(blh_read(fpin[1]));
      if (SWAP_ENDIAN) {
	blh[set] =  swap_blh(blh[set]);
	// reverse mapping the blhid wrt bl set number
	//          blid[blh[set]->blhid]=set;
      }
    }
    // count sidebands
    numberSidebands = 1;
    for (set=0; set<nsets[1]; set++) {
      if( blh[set]->inhid == inh[smabuffer.scanskip]->inhid
	  &&  blh[set]->isb != blh[set+1]->isb) {
	numberSidebands = 2;
	break; 
      }
    }
    printf("NUMBER OF SIDEBANDS =%d\n",numberSidebands);
    // count receivers
    smaCorr.no_rxif = 1;
    for (set=0; set<nsets[1]; set++) {
      smabuffer.rx1=smabuffer.rx2=blh[set]->irec;
      if( blh[set]->inhid == inh[smabuffer.scanskip]->inhid
	  &&  blh[set]->irec != blh[set+1]->irec) {
        smabuffer.rx2=blh[set+1]->irec;
	smaCorr.no_rxif = 2;
	break;
      }
    }
    printf("NUMBER OF RECEIVERS =%d \n",smaCorr.no_rxif);

   switch(smabuffer.rx1) {
   case 0: printf("rx1->230\n"); break;
   case 1: printf("rx1->340\n"); break;
   case 2: printf("rx1->690\n"); break;
                         }
       
   if(smaCorr.no_rxif == 2)
   switch(smabuffer.rx2) {
   case 0: printf("rx2->230\n"); break;
   case 1: printf("rx2->340\n"); break;
   case 2: printf("rx2->690\n"); break;
                         }    

           if(smabuffer.sb==0) printf("Processing LSB\n");
           if(smabuffer.sb==1) printf("Processing USB\n");
 
    numberRxif = smaCorr.no_rxif;
    // check if the receiver  smabuffer.rxif==0 for all receivers
    if(smabuffer.rxif==-1) goto foundTheRx;
    for (set=0; set<nsets[1]; set++) {
      if( blh[set]->inhid == inh[smabuffer.scanskip]->inhid
	  &&  blh[set]->irec == smabuffer.rxif) 
	goto foundTheRx;       } 
    printf("ERROR: there is no receiver %d in this data set.\n",
	   smabuffer.rxif);
    if(smaCorr.no_rxif == 1)
    printf("       please try it again with rxif=%d\n", smabuffer.rx1);
     if(smaCorr.no_rxif == 2)
     printf("       please try it again with rxif=%d or rxif=%d\n", 
     smabuffer.rx1,  smabuffer.rx2);
    exit(-1);
  foundTheRx:
// assign the rx id to load in the case of dual rx
      if(smabuffer.rxif==-1||smaCorr.no_rxif==1) {rxlod=0;
      } else {
      if(smabuffer.rxif==smabuffer.rx1) rxlod=0;
      if(smabuffer.rxif==smabuffer.rx2) rxlod=1;
            }

    tsys = (struct bltsys **) malloc(nsets[1]*sizeof( struct bltsys *));
    for (set=0;set<nsets[1];set++) {
      tsys[set] = (struct bltsys *)malloc(sizeof(struct bltsys ));
      if (tsys[set] == NULL ){
	printf("ERROR: Memory allocation for tsys failed trying to allocate %d bytes\n",
	       nsets[1]*sizeof(struct bltsys));
	exit(-1);
      }
    }
    atsys = (struct anttsys **) malloc(nsets[0]*sizeof( struct anttsys *));
    for (set=0;set<nsets[0];set++) {
      atsys[set] = (struct anttsys *)malloc(sizeof(struct anttsys ));
      if (atsys[set] == NULL ){
	printf("ERROR: Memory allocation for atsys failed trying to allocate %d bytes\n",
	       nsets[0]*sizeof(struct anttsys));
	exit(-1);
      }
    }
// allocate memory for wts once in the case of do both sb
         if(smabuffer.dsb!=1||(smabuffer.dsb==1&&smabuffer.sb==0)){
    wts = (struct wtt **) malloc(nsets[0]*sizeof( struct wtt *));
    for (set=0;set<nsets[0];set++) {
      wts[set] = (struct wtt *)malloc(sizeof(struct wtt ));
      if (wts[set] == NULL ){
        printf("ERROR: Memory allocation for wts failed trying to allocate %d bytes\n",
               nsets[0]*sizeof(struct wtt));
        exit(-1);
      }
    } }


    /* loading baselines */
    blhset =0;
    { 
      int blnset;
      int blset;
      int inhid_hdr;
      int blhid_hdr;
      blhid_hdr = 0;
      inhid_hdr = inh[0]->inhid;
      blnset = 0;
      if (smabuffer.rxif!=-1) {


   switch(smabuffer.rxif) {
   case 0: printf("to load rx->230 visdata\n"); break;
   case 1: printf("to load rx->340 visdata\n"); break;
   case 2: printf("to load rx->690 visdata\n"); break;
                         }
      } else {
	printf("to load data for all receivers.\n");
      }
      //
      // bln will be used in configuring spectra
      //
      set=0;
      for (blset=0; blset < nsets[1]; blset++) { 
	// loading baseline based structure
	tsys[blset]->blhid=blh[blset]->blhid;
	tsys[blset]->inhid=blh[blset]->inhid;
        tsys[blset]->blsid=blh[blset]->blsid;
	tsys[blset]->isb  =blh[blset]->isb;
	tsys[blset]->irec =blh[blset]->irec;
        tsys[blset]->ipol = ipolmap(blh[blset]->ipol);
	tsys[blset]->itel1=blh[blset]->itel1;
	tsys[blset]->itel2=blh[blset]->itel2;
	
	// assign baseline id handr
	if(blset==0) blhid_hdr = blset;
	if(blh[blset]->inhid!=inh[set]->inhid) {
	  set++;
	  blhid_hdr = blset;  
	}
// check if set exceeds nsets[0]
       if(set==nsets[0]) goto handling_blarray;        
// select side band 
	if(blh[blset]->inhid==inh[set]->inhid) {
// get the baseline id for the first bl in the new integration
        if(blh[blset]->inhid>blh[0]->inhid){
        if(smaCorr.no_rxif==2&&blh[blset-1]->irec==smabuffer.rx1
       &&blh[blset]->irec==smabuffer.rx2) {
        blid_intchng[set] = blh[blset]->blhid;
/*        printf("blid_intchng blh rx sb blh->inhid %d %d %d %d %d\n", 
        blid_intchng[set],
        blh[blset]->blhid, blh[blset]->irec, blh[blset]->isb,
        blh[blset]->inhid);
        printf("blid_intchng blh rx sb blh->inhid %d %d %d %d %d\n",
        blid_intchng[set],
        blh[blset-1]->blhid, blh[blset-1]->irec, blh[blset-1]->isb,
        blh[blset-1]->inhid);
*/
        }}
// choose rx
	  if(blh[blset]->irec==smabuffer.rxif||smabuffer.rxif==-1) {
// for the first set of integration, take the 1st baseline
	    if(set==0&&blh[blset]->isb==smabuffer.sb) {
	      bln[set]->inhid = blh[blset]->inhid;
	      bln[set]->blhid = blh[blset]->blhid;
	      bln[set]->isb   = blh[blset]->isb;
	      bln[set]->irec  = blh[blset]->irec; 
	      inhid_hdr       = blh[blset]->inhid;
/* printf("bln[set]->blhid  bln[set]->irec bln[set]->isb %d %d %d \n",
 bln[set]->blhid,
                 bln[set]->irec, bln[set]->isb); 
*/

	    }
	    // for the successive integration set, take the 1st baseline
	    // right after change of integration.
	    if(blh[blset]->inhid>inhid_hdr&&blh[blset]->isb==smabuffer.sb){
	      bln[set]->inhid = blh[blset]->inhid;
	      bln[set]->blhid = blh[blset]->blhid;
	      bln[set]->isb   = blh[blset]->isb;
	      bln[set]->irec  = blh[blset]->irec;
	      inhid_hdr       = blh[blset]->inhid;
	    }
	  }


	  // loading data to baseline coordinate structure
	  // convert ovro sign convention to miriad convention
	  // by multiplying a negative sign to uvw coordinates
	  uvwbsln[set]->uvwID[blset-blhid_hdr].u = -blh[blset]->u;
	  uvwbsln[set]->uvwID[blset-blhid_hdr].v = -blh[blset]->v;
	  uvwbsln[set]->uvwID[blset-blhid_hdr].w = -blh[blset]->w;
	  uvwbsln[set]->inhid = blh[blset]->inhid;
	  uvwbsln[set]->uvwID[blset-blhid_hdr].blhid = blh[blset]->blhid;
	  uvwbsln[set]->uvwID[blset-blhid_hdr].blsid = blh[blset]->blsid;
	  uvwbsln[set]->uvwID[blset-blhid_hdr].blcode =
	    blh[blset]->itel1*256+blh[blset]->itel2;
	  uvwbsln[set]->uvwID[blset-blhid_hdr].isb = blh[blset]->isb;
	  uvwbsln[set]->uvwID[blset-blhid_hdr].irec = blh[blset]->irec;
//	  uvwbsln[set]->uvwID[blset-blhid_hdr].ipol = blh[blset]->ipol;
	  // polarization
        uvwbsln[set]->uvwID[blset-blhid_hdr].ipol = ipolmap(blh[blset]->ipol);
	  // counting baseline for each integration set
	  uvwbsln[set]->n_bls++;
	}
	//     printf("set numberBaselines %d %d\n", set, numberBaselines);     
	numberBaselines=uvwbsln[set]->n_bls;
      } /* blset */
    }
handling_blarray:
    /* set antennas */
    {
      int bset;
// initialize blarray

            for (i=0; i< MAXANT; i++) {
             for (j=0; j< MAXANT; j++) {
          blarray[i][j].ee = 0.0;
           blarray[i][j].nn = 0.0;
            blarray[i][j].uu = 0.0;
             blarray[i][j].itel1 = 0;
              blarray[i][j].itel2 = 0;
               blarray[i][j].blid  = 0;
              } }
      bset = smabuffer.scanskip*numberBaselines;
      blarray[blh[bset]->itel1][blh[bset]->itel2].ee = blh[bset]->ble ;
      blarray[blh[bset]->itel1][blh[bset]->itel2].nn = blh[bset]->bln ;
      blarray[blh[bset]->itel1][blh[bset]->itel2].uu = blh[bset]->blu ;
      smabuffer.nants=1;
      for (set=bset+1;set<nsets[1];set++) {
	if(blarray[blh[bset]->itel1][blh[bset]->itel2].ee != blh[set]->ble) {
	  blarray[blh[set]->itel1][blh[set]->itel2].ee = blh[set]->ble;
	  blarray[blh[set]->itel1][blh[set]->itel2].nn = blh[set]->bln;
	  blarray[blh[set]->itel1][blh[set]->itel2].uu = blh[set]->blu;
	  blarray[blh[set]->itel1][blh[set]->itel2].itel1 = blh[set]->itel1;
	  blarray[blh[set]->itel1][blh[set]->itel2].itel2 = blh[set]->itel2;
	  blarray[blh[set]->itel1][blh[set]->itel2].blid  = blh[set]->blsid;
//          printf(" ant1 ant2 e n u %d %d %f %f %f\n", blh[set]->itel1,
//          blh[set]->itel2,blh[set]->ble, blh[set]->bln, blh[set]->blu);
	  smabuffer.nants++;       }
	else
	  {smabuffer.nants = (int)((1+sqrt(1.+8.*smabuffer.nants))/2);
	  printf("mirRead: number of antenna =%d are found.\n", smabuffer.nants);
	  goto blload_done;}
      } /* set */
    }
  blload_done:
    free(blh);
    
    if (SWAP_ENDIAN) {
      printf("FINISHED READING  BL HEADERS (endian-swapped)\n");
    } else {
      printf("FINISHED READING  BL HEADERS\n");
    }
    // assign memory to enh
    enh = (struct ant_def **) malloc(nsets[4]*sizeof( struct ant_def *));
    for (set=0;set<nsets[4];set++) {
      enh[set] = (struct ant_def *)malloc(sizeof(struct ant_def ));
      if (enh[set] == NULL ){
	printf("ERROR: Memory allocation for enh failed for %d bytes\n",
	       nsets[4]*sizeof(struct ant_def));
	exit(-1);
      }
    }
    // make sma engineer data buffer
    smaEngdata = (struct smEng **) malloc(nsets[0]*sizeof( struct smEng *));
    for (set=0;set<nsets[0];set++) {
      smaEngdata[set] = (struct smEng *)malloc(sizeof(struct smEng ));
      if (smaEngdata[set] == NULL ){
	printf("ERROR: Memory allocation for smaEngdata failed for %d bytes\n",
	       nsets[0]*sizeof(struct smEng));
	exit(-1);
      }
    }
    // skip engineer data reading because the engineer
    // file was problem for the two receivers case. 05-2-25
    //printf("doeng = %d\n", smabuffer.doeng);
    if (smabuffer.doeng!=1) {
      goto engskip;
    } else {
      // read engineer data
      for (set=0;set<nsets[4];set++) {
	*enh[set] = *(enh_read(fpin[4]));
        if (SWAP_ENDIAN) {
	  enh[set]=swap_enh(enh[set]);
	}         } 
      // store sma engineer data to smaEngdata 
      inhset=0;
      for (set=0;set<nsets[4];set++) {
        if(enh[set]->inhid!=inh[inhset]->inhid) inhset++;
        if(inhset<nsets[0]) {
	  smaEngdata[inhset]->inhid = enh[set]->inhid;
	  smaEngdata[inhset]->ints  = enh[set]->ints;
	  smaEngdata[inhset]->antpad_no[enh[set]->antennaNumber]
	    = enh[set]->padNumber;
	  smaEngdata[inhset]->antenna_no[enh[set]->antennaNumber]
	    = enh[set]->antennaNumber;
	  smaEngdata[inhset]->lst   = enh[set]->lst;
	  smaEngdata[inhset]->dhrs  = enh[set]->dhrs;
	  smaEngdata[inhset]->ha    = enh[set]->ha;
	  smaEngdata[inhset]->el[enh[set]->antennaNumber]
	    = enh[set]->actual_el;
	  smaEngdata[inhset]->az[enh[set]->antennaNumber]
	    = enh[set]->actual_az;
	  smaEngdata[inhset]->tsys[enh[set]->antennaNumber]
	    = enh[set]->tsys;
	  smaEngdata[inhset]->tamb[enh[set]->antennaNumber]
	    = enh[set]->ambient_load_temperature;
	} 
      }
    }
    if (SWAP_ENDIAN) {
      printf("FINISHED READING EN HEADERS (endian-swapped)\n");
    } else {
      printf("FINISHED READING EN HEADERS\n");
    }
  engskip:
    //free(smaEngdata);
    free(enh);
    // initialize the antenna positions
    smabuffer.nants=8;
    for (i=1; i < smabuffer.nants+1; i++) {
      antenna[i].x = 0.;
      antenna[i].y = 0.;
      antenna[i].z = 0.;
      antenna[i].x_phs = 0.;
      antenna[i].y_phs = 0.;
      antenna[i].z_phs = 0.;
      antenna[i].axisoff_x = 0.;
      antenna[i].axisoff_y = 0.;
      antenna[i].axisoff_z = 0.;
    }
    
    // set up the reference antenna
    antenna[smabuffer.refant].x = 0.;
    antenna[smabuffer.refant].y = 0.;
    antenna[smabuffer.refant].z = 0.;
    
    // derive antenna position in local coordinate system
    // mir stores the position in float
    for (i=1; i < smabuffer.nants+1; i++) {
      if(i<smabuffer.refant) {
	antenna[i].x = (double)blarray[i][smabuffer.refant].ee 
	  - antenna[smabuffer.refant].x;
	antenna[i].y = (double)blarray[i][smabuffer.refant].nn 
	  - antenna[smabuffer.refant].y;
	antenna[i].z = (double)blarray[i][smabuffer.refant].uu 
	  - antenna[smabuffer.refant].z;
//        printf(" ant ee nn uu %d %f %f %f\n",i,
//            blarray[i][smabuffer.refant].ee,
//            blarray[i][smabuffer.refant].nn,
//            blarray[i][smabuffer.refant].uu);
                         
      } else {
         if(i==smabuffer.refant) {
         antenna[smabuffer.refant].x = 0.;
         antenna[smabuffer.refant].y = 0.;
         antenna[smabuffer.refant].z = 0.;
//            printf(" ant ee nn uu %d %f %f %f\n",i,
//            blarray[smabuffer.refant][i].ee,
//            blarray[smabuffer.refant][i].nn,
//           blarray[smabuffer.refant][i].uu);
                                 } else {

	
	antenna[i].x = (double)blarray[smabuffer.refant][i].ee
	  - antenna[smabuffer.refant].x;
	antenna[i].x = - antenna[i].x;
	antenna[i].y = (double)blarray[smabuffer.refant][i].nn
	  - antenna[smabuffer.refant].y;
	antenna[i].y = - antenna[i].y;
	antenna[i].z = (double)blarray[smabuffer.refant][i].uu
	  - antenna[smabuffer.refant].z;
	antenna[i].z = - antenna[i].z;
//            printf(" ant ee nn uu %d %f %f %f\n",i,
//            blarray[smabuffer.refant][i].ee,
//            blarray[smabuffer.refant][i].nn,
//            blarray[smabuffer.refant][i].uu);
                           }
      }
    }
    
    
    // calculate the geocentric coordinates from local 
    // coordinates
    {
      struct xyz geocxyz[MAXANT];
      for (i=1; i < smabuffer.nants+1; i++) {
	geocxyz[i].x = (antenna[i].z)*cos(smabuffer.lat)
	  - (antenna[i].y)*sin(smabuffer.lat);
	geocxyz[i].y = (antenna[i].x);
	geocxyz[i].z = (antenna[i].z)*sin(smabuffer.lat)
	  + (antenna[i].y)*cos(smabuffer.lat);
// printf("ant x y z %d %f %f %f\n", i, antenna[i].x,antenna[i].y,antenna[i].z);
      }
      printf("NUMBER OF ANTENNAS =%d\n", smabuffer.nants);
      //
      // maximum antenna number for the array is 8 
      //
      smabuffer.nants = 8;   
      for (i=1; i < smabuffer.nants+1; i++) {
	sprintf(antenna[i].name, "AN%d", i);
      }     
      // the positions of antennas need to check  
      // antpos on 2005 feb 16
      antxyz[1].x = 4.4394950000000000e+00;
      antxyz[2].x = -5.7018977000000000e+00;
      antxyz[3].x = -5.0781509999999996e-01;
      antxyz[4].x = 5.2273398999999996e+00;
      antxyz[5].x = -1.7918341999999999e+01;
      antxyz[6].x = 0.0000000000000000e+00;
      antxyz[7].x = -1.6564844999999998e+01;
      antxyz[8].x = -6.4019909999999998e+00;
      
      antxyz[1].y = -6.3875615000000003e+01; 
      antxyz[2].y = -1.8985175000000002e+01;
      antxyz[3].y = -2.5154143000000001e+01;
      antxyz[4].y = -2.0077679000000000e+01;
      antxyz[5].y = -5.9557980000000001e+01;
      antxyz[6].y =  0.0000000000000000e+00;
      antxyz[7].y = -2.7025637000000000e+01;
      antxyz[8].y = -6.8001356999999999e+01;
      
      antxyz[1].z = -2.1837547000000001e+01; 
      antxyz[2].z = 1.5609299999999999e+01;
      antxyz[3].z =  1.2797959999999999e+00; 
      antxyz[4].z = -1.4844139000000000e+01;
      antxyz[5].z = 3.0068384999999999e+01;
      antxyz[6].z = 0.0000000000000000e+00;
      antxyz[7].z = 3.0769280999999999e+01;
      antxyz[8].z = 3.6370589999999998e+00;

// reading antenna position from ASCII file
      if(smabuffer.readant > 0) {
double xyzpos;
       for (i=1; i < smabuffer.readant+1; i++) {
         geocxyz[i].x = smabuffer.antpos[1+(i-1)*3];
         geocxyz[i].y = smabuffer.antpos[2+(i-1)*3];
         geocxyz[i].z = smabuffer.antpos[3+(i-1)*3];
         xyzpos = geocxyz[i].x+geocxyz[i].y+geocxyz[i].z;
      if(xyzpos==0) smabuffer.refant = i;
                                 }
                  }



    printf("Geocentrical coordinates of antennas (m), reference antenna=%d\n",
	     smabuffer.refant); 
      for (i=1; i < smabuffer.nants+1; i++) {
       
	printf("ANT x y z %s %11.4f %11.4f %11.4f\n",
	       antenna[i].name,
	       geocxyz[i].x,
	       geocxyz[i].y,
	       geocxyz[i].z);
      }
      //
      // convert geocentrical coordinates to equatorial coordinates
      // of miriad system y is local East, z is parallel to pole
      // Units are nanosecs. 
      // write antenna dat to uv file 
      //
      printf("Miriad coordinates of antennas (nanosecs), reference antenna=%d\n",
	     smabuffer.refant);
      r = sqrt(pow(geocxyz[smabuffer.refant].x,2) + 
	       pow(geocxyz[smabuffer.refant].y,2));
      if(r>0) {
	cost = geocxyz[smabuffer.refant].x / r;
	sint = geocxyz[smabuffer.refant].y / r;
	z0   = geocxyz[smabuffer.refant].z; 
      } else {
	cost = 1;
	sint = 0;
	z0 = 0; 
      }
      
      for (i=1; i < smabuffer.nants+1; i++) {      
	tmp  = ( geocxyz[i].x) * cost + (geocxyz[i].y)*sint - r;
	antpos[i-1] 
	  = (1e9/DCMKS) * tmp;
	printf("ANT x y x %s %11.4f ", antenna[i].name, antpos[i-1]);
	tmp  = (-geocxyz[i].x) * sint + (geocxyz[i].y) * cost;
	antpos[i-1+smabuffer.nants] 
	  = (1e9/DCMKS) * tmp;
	printf("%11.4f ", antpos[i-1 + smabuffer.nants]);
	antpos[i-1+2*smabuffer.nants] 
	  = (1e9/DCMKS) * (geocxyz[i].z-z0);
	printf("%11.4f\n", antpos[i-1+2*smabuffer.nants]);
      }
      
    }
    
    tno  = smabuffer.tno;
    nnants = 3*smabuffer.nants;
    uvputvrd_c(tno,"antpos", antpos, nnants);
    // setup source 
    cdh = (struct codeh_def **) malloc(nsets[3]*sizeof( struct codeh_def *));
    for (set=0;set<nsets[3];set++) {
      cdh[set] = (struct codeh_def *)malloc(sizeof(struct codeh_def ));
      if (cdh[set] == NULL ){
	printf("ERROR: Memory allocation for cdh failed for %d bytes.\n",
	       nsets[3]*sizeof(struct codeh_def));
	exit(-1);
      }
    }
    for (set=0;set<nsets[3];set++){
      *cdh[set] = *(cdh_read(fpin[3]));
      if (SWAP_ENDIAN) {
	cdh[set]=swap_cdh(cdh[set]);
      }
    }
    if (SWAP_ENDIAN) {
      printf("FINISHED READING  CD HEADERS (endian-swapped)\n");
    }else {
      printf("FINISHED READING  CD HEADERS\n");
    }
    sourceID = 0;
    for (set=0;set<nsets[3];set++){
      // decode the ids 
      if((cdh[set]->v_name[0]=='b'&&cdh[set]->v_name[1]=='a')&&
        cdh[set]->v_name[2]=='n'){
        if(smabuffer.spskip[0]==0) { 
	spcode[cdh[set]->icode]=spdecode(&cdh[set]);
          } else { 
           if(cdh[set]->icode < smabuffer.spskip[0]) {
            spcode[cdh[set]->icode]=spdecode(&cdh[set]);
              } else {
        spcode[cdh[set]->icode]
        =spdecode(&cdh[set])-smabuffer.spskip[1];
              }
                }
      }
      // decode the julian date for from the observing date 
      if((cdh[set]->v_name[0]=='r'&&cdh[set]->v_name[1]=='e')&&
	 cdh[set]->v_name[2]=='f'){
	jday = juliandate(&cdh[set]);      }
    }
    // decode velocity type 
    // 2 CODE NAME vctype  STRING vlsr  icode 0 ncode 1
    // 3 CODE NAME vctype  STRING cz    icode 1 ncode 1
    // 4 CODE NAME vctype  STRING vhel  icode 2 ncode 1
    // 5 CODE NAME vctype  STRING pla   icode 3 ncode 1
    
    if(inh[1]->ivctype==0)
      strcpy(multisour[sourceID].veltyp, "VELO-LSR");
    if(inh[1]->ivctype==2)
      strcpy(multisour[sourceID].veltyp, "VELO-HEL");
    if(inh[1]->ivctype!=0 && inh[1]->ivctype!=2) {
      printf("ERROR: veltype ivctype=%d is not supported.\n",
	     inh[1]->ivctype);
      exit(-1); 
    }
    
    uvputvra_c(tno,"veltype", multisour[sourceID].veltyp);
// decode the source information 
    for (set=0;set<nsets[3];set++){
    if(cdh[set]->v_name[0]=='s'&&cdh[set]->v_name[1]=='o') {
//	sourceID++;
        sourceID = cdh[set]->icode;
// parsing the source name and trim the junk tail
	for(i=0; i<9; i++) {
	  sours[i]=cdh[set]->code[i];
	  if(cdh[set]->code[i]==32||cdh[set]->code[i]==0||i==8)
	    sours[i]='\0';
	}
           sprintf(multisour[sourceID].name, "%s", sours);
          for(i=2; i< sourceID; i++) {
          if(strcmp(multisour[i].name, sours)==0) {
// copy the original source name to smasours
          for(i=0; i<33; i++) {
          smasours[i]=cdh[set]->code[i];
          if(cdh[set]->code[i]==32||cdh[set]->code[i]==0||i==32)
            smasours[i]='\0';
                               }
         
          oka=okb=okc=okd=0;
          for(i=2; i< sourceID; i++) {
          if(multisour[sourceID].name[7]=='a') oka=-1;
          if(multisour[sourceID].name[7]=='b') okb=-1;
          if(multisour[sourceID].name[7]=='c') okc=-1;
          if(multisour[sourceID].name[7]=='d') okd=-1; 
                      }
          if(oka==0) {sours[7]='a';
          sprintf(multisour[sourceID].name, "%s", sours);
                    } else {           
          if(okb==0){sours[7]='b';
          sprintf(multisour[sourceID].name, "%s", sours);
                     } else {
          if(okc==0) {sours[7]='c';
          sprintf(multisour[sourceID].name, "%s", sours);
                      } else {
          if(okd==0) {sours[7]='d';
          sprintf(multisour[sourceID].name, "%s", sours); }
                              
                             }
                            }
                           }

    printf("Warning: The original name: '%s' is renamed to '%s'\n", 
    smasours, multisour[sourceID].name);
                                                       }
                                                 }
    multisour[sourceID].sour_id = cdh[set]->icode;
	//         printf("cdh[set]->code=%s\n", cdh[set]->code);
    inhset=0;
	while (inh[inhset]->souid!=multisour[sourceID].sour_id)
	  inhset++;
	//         printf("inh[inhset]->rar to c %s inhid=%d\n",
	//         (char *)rar2c(inh[inhset]->rar),inh[inhset]->inhid);
	multisour[sourceID].ra = inh[inhset]->rar;
	multisour[sourceID].dec = inh[inhset]->decr;
	sprintf(multisour[sourceID].equinox, "%f", inh[inhset]->epoch);
	multisour[sourceID].freqid =-1;
	multisour[sourceID].inhid_1st=inh[inhset]->inhid;
	// calculate the apparent position coordinates from j2000 coordinates 
	{ 
	  double obsra,obsdec,r1,d1;
	  double julian2000=2451544.5;
	  precess(julian2000,
		  multisour[sourceID].ra,
		  multisour[sourceID].dec, jday, &obsra, &obsdec);
	  nutate(jday,obsra,obsdec,&r1,&d1);
	  aberrate(jday,r1,d1,&obsra,&obsdec);
	  multisour[sourceID].ra_app   = obsra;
	  multisour[sourceID].dec_app  = obsdec;
	  multisour[sourceID].qual     = 0;
	  multisour[sourceID].pmra     = 0.;
	  multisour[sourceID].pmdec    = 0.;
	  multisour[sourceID].parallax = 0.;
  //         strcpy(multisour[sourceID].veltyp, "lsr");
	  strcpy(multisour[sourceID].veldef, "radio");
	  strcpy(multisour[sourceID].calcode, "c");
	}
      }
    }
// setup correlator  
// sph1 is a single set of spectra, assign memory to it.    
    sph1 = (struct sph_def *) malloc(sizeof( struct sph_def ));
// spn is a buffer for configuring spectra with an array length
// of the total number of integration sets.
    spn  = (struct sph_config **) malloc(nsets[0]*sizeof( struct sph_config *));
    for (set=0; set<nsets[0]; set++) {
      spn[set] = (struct sph_config *)malloc(sizeof(struct sph_config ));
      if (spn[set] == NULL ){
	printf("ERROR: Memory allocation for sph_config failed for %d bytes\n",
	       nsets[0]*sizeof(struct sph_config));
	exit(-1);
      }
    }
//
// determine if perform high resolution mode:
// data obtained earlier than 2006-may-12 is to be done
// with a low resolution mode;
// data obtained later than  2006-may-12 is handled in the high
// resolution mode unless this mode turned off by options=nohspr
//
    if(smabuffer.highrspectra==0&&jday>2453867) {
      smabuffer.highrspectra = 1;
      } else { smabuffer.highrspectra =-1; }
    
    { 
      int inhid_hdr;
      int blhid_hdr;
      int sphid_hdr;
      int blset;
      int inset;
      int nspectra;
      nspectra=0;
      blset     =  0;
      sphid_hdr =  0;
      blhid_hdr = -10;
      inhid_hdr = -10;
      firstsp   = -1;
      lastsp    = -1;
      numberSpectra = 0;
      // define baseline id used in pursing the spectral
      // configuration. 
      // integration set = smabuffer.scanskip
      // blhset=1 , the second baseline of the integration
      blhset=1;
      //      printf("smabuffer.scanskip=%d\n", smabuffer.scanskip);
      blhid = uvwbsln[smabuffer.scanskip]->uvwID[blhset].blhid;
      rewind(fpin[2]);
      // spn starts from 0
      inset = smabuffer.scanskip;
      for (set=sphid_hdr;set<nsets[2]; set++) {
	*sph1 = *(sph_read(fpin[2]));
	if (SWAP_ENDIAN) {
	  sph1 =  swap_sph(sph1);
	}
	if(sph1->blhid==blhid)
	  numberSpectra++;
	if(sph1->blhid==blhid&&numberSpectra==0) {
	  // assume 25 chunks per baseline
	  numberSpectra=25;
	  smabuffer.scanskip++;
	  inset = smabuffer.scanskip;
	}
// load baseline based tsys structure
	if(sph1->blhid==tsys[0]->blhid) nspectra++;
//      printf("nspectra = %d\n", nspectra);
//      printf("sph1->inhid tsys[blset]->inhid sph1->blhid tsys[blset]->blhid blset %d %d %d %d %d\n",
//        sph1->inhid,tsys[blset]->inhid,sph1->blhid,tsys[blset]->blhid,blset);
	if(sph1->blhid==tsys[blset]->blhid&&sph1->inhid==tsys[blset]->inhid) {
// printf("sphset blid iband tsys %d  %d  %d %f \n",set,sph1->blhid, sph1->iband,sph1->tssb);
 tsys[blset]->tssb[sph1->iband] = sph1->tssb;
// loading online flagging information
        if(tsys[blset]->ipol < -4||sph1->iband!=0) {
       wts[inset]->wt[sph1->iband-1][-4-tsys[blset]->ipol][tsys[blset]->blsid][tsys[blset]->isb][tsys[blset]->irec] = sph1->wt; 
                                     }
          else {
         if(sph1->iband!=0)
 wts[inset]->wt[sph1->iband-1][-tsys[blset]->ipol][tsys[blset]->blsid][tsys[blset]->isb][tsys[blset]->irec] = sph1->wt;
           }   
	if(smabuffer.highrspectra !=1){  
        if(sph1->iband==nspectra-1) blset++;
        } else {
        if(sph1->iband==24) blset++;
        
        }
	}
	// purse the spectral configuration
	// check up inhid to work on the same set of integration 
	if(sph1->inhid>inh[inset]->inhid) inset++ ;
	// check up inhid and blhid to work on the same integration set
	//       and the baseline set with the sidebband and rx as
	//       is desired.
	
	if(sph1->inhid==inh[inset]->inhid&&sph1->blhid==bln[inset]->blhid){
	  spn[inset]->sphid                = sph1->sphid;
	  spn[inset]->inhid                = sph1->inhid;
	  spn[inset]->iband[sph1->iband]   = sph1->iband;
	  
	  // velocity with respect to the rest frame
          // given by the SMA on-line system. This is 
          // only meaningful to the line transition at the
          // rest frequency 
	  spn[inset]->vel[sph1->iband]     = sph1->vel;
	  spn[inset]->vres[sph1->iband]    = sph1->vres;
	  spn[inset]->ivtype               = sph1->ivtype;
	  // sky frequency (corrected for a part of the Doppler
          // velocity, the diurnal term and a part of the annual term) 
	  spn[inset]->fsky[sph1->iband]    = sph1->fsky;
	  spn[inset]->fres[sph1->iband]    = sph1->fres;
   if(smaCorr.no_rxif!=2) spn[inset]->nch[sph1->iband][0]  = sph1->nch;
	  spn[inset]->dataoff              = sph1->dataoff;
	  spn[inset]->rfreq[sph1->iband]   = sph1->rfreq;
	  spn[inset]->isb                  = bln[inset]->isb;
	  spn[inset]->irec                 = bln[inset]->irec;
	  spn[inset]->souid                = inh[inset]->souid;
	  if(sph1->iband==0) spn[inset]->basefreq = sph1->fsky;                
	}
      if(smaCorr.no_rxif==2&&sph1->inhid==inh[inset]->inhid) {
      if(sph1->blhid==blid_intchng[inset])
        spn[inset]->nch[sph1->iband][1] = sph1->nch;    
        
      if(sph1->blhid==(blid_intchng[inset]-1))
        spn[inset]->nch[sph1->iband][0] = sph1->nch; 
/*       if(sph1->blhid==blid_intchng[inset])     
       printf("inset blhid nch[%d][0] nch[i][1]  %d %d %d %d \n",
        sph1->iband,
        inset,
        sph1->blhid,
        spn[inset]->nch[sph1->iband][0],  
        spn[inset]->nch[sph1->iband][1]); 
*/                             }
          lastinset=inset;
	if(inset==smabuffer.scanskip+smabuffer.scanproc) { 
	  goto sphend; 
	}
      }
    }
sphend:
     if (smabuffer.scanskip+smabuffer.scanproc>nsets[0]) {
printf("Hits the end of integration %d\n", nsets[0]);
printf("Reset 'nscans=%d, ' and run it again.\n",smabuffer.scanskip);
         exit(-1);
                  }

// handling 2003 incompleted correlator data
if (smabuffer.spskip[0]==-2003) {
double spfreq[25];
double spfreqo[25];
int nchunk=0;
int minspID=25;
         for (i=1;i<SMIF+1;i++) {
         if (spn[smabuffer.scanskip]->fsky[i]>0.0) {
         nchunk++;
         spfreq[nchunk]=spn[smabuffer.scanskip]->fsky[i];
                            }
                }
         for (i=1;i<SMIF+1;i++) {
         if(spfreq[i]<100.0) {
         spfreq[i]=0.0;
         spcode[i]=0;        }
         spfreqo[i] =  spfreq[i];
         if (spcode[i]!=0&&spcode[i]<minspID) minspID=spcode[i];
                           }
         for (i=1;i<SMIF+1;i++) {
          spcode[i]=spcode[i]-minspID+1;
                           }
goto dat2003; }
                
              lastinset=lastinset-1;
              reset_id[0]= smabuffer.scanskip;
              for (iinset=smabuffer.scanskip;
                   iinset < lastinset;
                   iinset++) {
              for (i=0;i<SMIF+1;i++) {
              if (spn[iinset]->nch[i][0]
               !=spn[iinset+1]->nch[i][0]) {
         printf("Warning: The correlator was reconfigured at integrations=%d\n",
                   iinset+1);
                          ireset++;
                    reset_id[ireset] = iinset+1;
         printf("From    -> To\n");
                    for (i=1;i<SMIF+1;i++) {
         printf("s%2d:%3d -> :%3d\n", i, spn[iinset]->nch[i][0],
                 spn[iinset+1]->nch[i][0]);
                                       }
                                           }
                                 }
                              }
           if(smabuffer.mcconfig==0) {
           for (i=1;i<ireset+1;i++) {
           if(smabuffer.scanskip< reset_id[i]) {
    printf("Suggesttion:\n");
    printf("For a single correlator configuration per loading (recommended),\n");
    printf("reset nscans=%d,%d for the set of the first configuration data;\n", smabuffer.scanskip, reset_id[i]-1);
    printf("reset nscans=%d, for the set of the second configuration data.\n",  reset_id[i]);
    printf("Or choose options=mcconfig for multiple correlator configurations per loading (not recommended).\n");
    printf("Try it again.\n");
         exit(-1);
        }
                                    }
                                     }
dat2003:                                                                           
    printf("\n");
    printf("number of non-empty Spectra = %d\n", numberSpectra);
// reset number of spectra
    if(smabuffer.highrspectra==1) numberSpectra=25;
    if (SWAP_ENDIAN) {
      printf("FINISHED READING SP HEADERS (endian-swapped)\n");
    } else {
      printf("FINISHED READING SP HEADERS\n");
    }
// solve for tsys
    {
      int refant;
      int done;
      int pair1;
      int pair2;
      int blset;
// solve for tsys of a reference ante
      set=0;
      refant=0;
      pair1=0;
      pair2=0;
      blset=1;
      for(blset=0; blset<nsets[1]; blset++) {
	if(set==nsets[0]-1) goto next;
	if(tsys[blset]->inhid!=inh[set]->inhid) {set++; refant=0;}
// choose rx
	if(tsys[blset]->irec==smabuffer.rxif||smabuffer.rxif==-1) {
	  if(tsys[blset]->inhid==inh[set]->inhid&&tsys[blset]->isb==0) {
	    if(refant==0) {
	      refant=100;
	      atsys[set]->refant = tsys[blset]->itel1;
	      atsys[set]->tssb[atsys[set]->refant]=tsys[blset]->tssb[0];
	      atsys[set]->refpair1  = tsys[blset]->itel2;
	    }
	  }
	}
      }
    next:
      done=0;
      set=0;
      for(blset=0; blset<nsets[1]; blset++) {
	if(set==nsets[0]-1) goto nextnext;
	if(tsys[blset]->inhid!=inh[set]->inhid) {set++;done=0;}
	if(tsys[blset]->irec==smabuffer.rxif||smabuffer.rxif==-1) {
	  if(tsys[blset]->inhid==inh[set]->inhid&&tsys[blset]->isb==0) {
	    if(done==0) {
	      if(tsys[blset]->itel1==atsys[set]->refant&&atsys[set]->refpair1!=tsys[blset]->itel2) {
		atsys[set]->refpair2=tsys[blset]->itel2;
		atsys[set]->tssb[atsys[set]->refant] *=
		  tsys[blset]->tssb[0];
		done=100;
	      }
	    }
	  }
	}
      }
    nextnext:
      set=0;
      for(blset=0; blset<nsets[1]; blset++) {
	if(set==nsets[0]-1) goto nextnextnext;
	if(tsys[blset]->inhid!=inh[set]->inhid) {set++; refant=0;}
	if(tsys[blset]->irec==smabuffer.rxif||smabuffer.rxif==-1) {
	  if(tsys[blset]->inhid==inh[set]->inhid&&tsys[blset]->isb==0) {
	    if((tsys[blset]->itel1==atsys[set]->refpair2&&tsys[blset]->itel2==atsys[set]->refpair1)||
	       (tsys[blset]->itel2==atsys[set]->refpair2&&tsys[blset]->itel1==atsys[set]->refpair1)) {
	      atsys[set]->tssb[atsys[set]->refant] =
		atsys[set]->tssb[atsys[set]->refant]/tsys[blset]->tssb[0];
	    }
	  }
	}
      }
    nextnextnext:
      set=0;
      for(blset=0; blset<nsets[1]; blset++) {
	if(set==nsets[0]-1) goto nnextnextnext;
	if(tsys[blset]->inhid!=inh[set]->inhid) {set++;}
	if(tsys[blset]->irec==smabuffer.rxif||smabuffer.rxif==-1) {
	  if(tsys[blset]->inhid==inh[set]->inhid&&tsys[blset]->isb==0) {
	    if(tsys[blset]->itel1==atsys[set]->refant) {
	      atsys[set]->tssb[tsys[blset]->itel2] =
		tsys[blset]->tssb[0]*tsys[blset]->tssb[0]
		/ atsys[set]->tssb[atsys[set]->refant];
	    }
	    if(tsys[blset]->itel2==atsys[set]->refant) {
	      atsys[set]->tssb[tsys[blset]->itel1] =
		tsys[blset]->tssb[0]*tsys[blset]->tssb[0]
		/atsys[set]->tssb[atsys[set]->refant];
	    }
	  }
	}
      }
    nnextnextnext:
      //     for (set=0; set < nsets[0]; set++) {
      //        printf("%d tsys %f  %f \n", set, atsys[set]->tssb[3],
      //        atsys[set]->tssb[7]);
      //     printf("%d entsys %f %f --- %f %f\n", set, atsys[set]->tssb[3], 
      //          smaEngdata[set]->tsys[3], atsys[set]->tssb[7],
      //          smaEngdata[set]->tsys[7]); 
      //   printf("%d tsys %f %f \n", atsys[set]->tssb[3],atsys[set]->tssb[7]);
      //}
    printf("Decoded baseline-based Tsys\n");
    }
    // decode the doppler velocity

 if (jday <2453531.5 && smabuffer.circular == 1) {
      printf("Skip decoding the Doppler velocity\n");
     } else {
      double vabsolute;
      double fratio;
      for(set= smabuffer.scanskip; 
	  set < smabuffer.scanskip+smabuffer.scanproc+1; set++){
	// calculate doppler velocity from chunk 1; 
	// the rest chunk give the same value.
        // absolute velocity is based on form (2-229)
        // Astrophysical Formulae by K. R. Lang (Springer-Verlag 1974)
        // p194 
	fratio = (spn[set]->fsky[12]/spn[set]->rfreq[12]);
	vabsolute = (1. - fratio*fratio ) / (1. + fratio*fratio );
	vabsolute = vabsolute*DCMKS/1000.;
        // Derive the Doppler velocity from observer to the LSR
        // by taking out the Radial velocity (at chunk 1) of the source 
        // w.r.t. the LSR or HEL.
        // Based on the information from Taco after the discussion with
        // Dan Marrone and Jim Moran on 05aug26's meeting,
        // There are two operation modes that we have supported:
        // 1) non-planet source and 2) planet source;
        // For mode 1) the chunk frequency recorded in the MIR data header
        // is the sky frequency at the transit of the source being Doppler
        // tracked, i.e., the recorded frequency is the true sky frequency 
        // only corrected for the diurnal term due to the earth rotation.
        // Then, Jun-Hui started to build a patch in UVREDO for handling
        // this SMA specification and he found that additional amount
        // correction for the part of annual term that has been made to 
        // the sky frequency by the SMA online system. Then he tried to
        // decode the reference time corresonding to a zero value that 
        // the online system uses to take out the steady variation in 
        // the sky frequnecy due to annual term. Jun-Hui further consulted 
        // with Taco on 05Aug31. Taco commented that it should be the 
        // time at the source transit but he could be
        // wrong. Based on the observation on 2005Aug01 for the SgrB2 track,
        // Jun-Hui tried to figure out the reference time using veldop calculated
        // from UVREDO and the SMA veldop decoded from the MIR header.
        // Apparently, using the reference time at the source transit,
        // it gives large difference between result from Miriad and
        // the value of the veldop decoded from the SMA data. There are other
        // possible times for the reference: 1) the time at which 
        // the DopplerTrack command issued by the operator (Mon Aug  1 04:06:14 2005 (brownd)
        // for the SgrB2 observation; 2) the time at which the project command
        // is issued by the operator. (in the SgrB2 observations, the project
        // command was issued towice at Mon Aug  1 04:50:35 2005 (brownd)
        // and at Mon Aug  1 05:35:48 2005 (brownd). Among the three possible time,
        // the time 05:35:48 appears to be the closest to the reference time
        // (which gives precision in velocity of 0.0015 km/s). Therefore 
        // two approaches are used to implement the patch in uvredo to calculate
        // the residual Doppler velocity: 
        // a) Giving the reference time from which the offset in the annual 
        // term that has been corrected to the chunk frequency;
        // b) Using veldop decoded from the MIR header (the sky frequency/the 
        // radial velocity at a channel and the corresponding rest frequency). 
        // The value of decoded veldop is to be stored
        // in the variable 'veldop' but it is not the Tracked Doppler velocity;
        // In other document (users guide), it is called as residual Doppler
        // velocity.
        //  
        // For mode 2) the chunk frequency recorded in the MIR data header
        // is the sky frequency corrected for the radial velocity of the planet
        // at the moment when the dopplerTrack command is issued.
        //   
        spn[set]->veldop = vabsolute - spn[set]->vel[12];
//
// add back the radial velocity of the source
//
        spn[set]->veldop = spn[set]->veldop + smabuffer.vsource;
        spn[set]->smaveldop = vabsolute - spn[set]->vel[12]+ 
                 smabuffer.vsource;
// printf("%d veldop=%f %f\n", set, spn[set]->veldop,spn[set]->vel[12]);
      }
    printf("Decoded the Doppler velocity\n");
    }                   
// rewind the data file
    rewind(fpin[5]);
// start from the inhset = smabuffer.scanskip 
    inhset = smabuffer.scanskip;
    numberBaselines=  uvwbsln[inhset]->n_bls,
    printf("here we are!\n");
//
// parsing to handle spectral chunk data assuming
// 24 spectral chunks handled online regardless
// some of the chunks may be skipped;
// or only partial (< 24) chunks handled online
// (early days' data for exmaple).
//
   if(smabuffer.highrspectra==1) {
   smaCorr.n_chunk = SMIF;
    } else { 
// take out 1 for eliminating the continuum channel
   smaCorr.n_chunk = numberSpectra -1;
    }
// end of spectral configuration loading.
    
// initialize system velocity
    for (j=1; j<sourceID+1; j++) {
      for (i=1; i<smaCorr.n_chunk+1; i++) {
	multisour[j].sysvel[i] = 0.0;
      }           
    }
// print out sources observed
    max_sourid=0;
    for (i=1; i< sourceID+1; i++) {
      if(multisour[i].sour_id==0) 
        strcpy(multisour[i].name, "skipped!");
      printf("source: %-21s id=%2d RA=%13s ", 
	     multisour[i].name,
	     multisour[i].sour_id, 
	     (char *)rar2c(multisour[i].ra));
      printf("DEC=%13s\n", (char *)decr2c(multisour[i].dec));
    if(multisour[i].sour_id>max_sourid) max_sourid=multisour[i].sour_id;
    }


    
// now loading the smabuffer, parameters for headers
    
// initialize Tsys
    for (i=1; i<smabuffer.nants+1; i++){
      smabuffer.tsys[i-1]=0;
    } 
    
// initialize the polcode
    for(j=1;j<smaCorr.n_chunk+1;j++) {
           smabuffer.nstoke[j-1]=4;
      for(i=1; i<smabuffer.nstoke[j-1]; i++) {
	for (k=1; k<SMBAS+1; k++) {
	  smabuffer.polcode[j-1][i-1][k-1] = -5;
	}
      }
    }
// initialize the tsys for the polarization components
    for(j=1;j<smabuffer.nants+1;j++) {
      for(i=1;i<smaCorr.n_chunk+1;i++) {
	smabuffer.xtsys[i-1][j-1]=0.;
	smabuffer.ytsys[i-1][j-1]=0.;
	smabuffer.xyphase[i-1][j-1]=0.;
	smabuffer.xyamp[i-1][j-1]=1.;
      }
    }
    
// initialize the sampler for the polarization components

    for(k=0; k<3; k++) {
      for(j=1;j<smabuffer.nants+1;j++) {
	for(i=1;i<smaCorr.n_chunk+1;i++) {
	  smabuffer.xsampler[i-1][j-1][k]=0.;
	  smabuffer.ysampler[i-1][j-1][k]=0.;
	}
      }
    }
// assign actual number of spectral chunks to smabuffer  
//   smabuffer.nifs = smaCorr.n_chunk;
     smabuffer.nifs = numberSpectra-1;
// initialize pnt flags 
    for (i=1; i<SMIF+1; i++)  {
      for (j=1; j<SMPOL+1; j++) {
	for (k=1; k<SMBAS+1; k++) {
	  for (l=1; l<3; l++) { 
	    for (m=1; m<SMRX+1; m++) {
              smabuffer.flag[i-1][j-1][k-1][l-1][m-1]=-1;
	       smabuffer.pnt[i-1][j-1][k-1][l-1][m-1]= 0;
	    }
	  }
	}
      }
    }
// reverse the spectral chunk order for blocks 1 2 3 4 
    if(smabuffer.doChunkOrder==1) {
      frcode[1]=4;
      frcode[2]=3;
      frcode[3]=2;
      frcode[4]=1;
      
      frcode[5]=8;
      frcode[6]=7;
      frcode[7]=6;
      frcode[8]=5;
      
      frcode[9]=12;
      frcode[10]=11;
      frcode[11]=10;
      frcode[12]=9;
      
      frcode[13]=13;
      frcode[14]=14;
      frcode[15]=15;
      frcode[16]=16;
      
      frcode[17]=17;
      frcode[18]=18;
      frcode[19]=19;
      frcode[20]=20;
      
      frcode[21]=21;
      frcode[22]=22;
      frcode[23]=23;
      frcode[24]=24;
    }
    
// print the side band to be processed
    switch(smabuffer.sb) {
    case 0:
      printf("LSB only\n");
      break;
    case 1:
      printf("USB only\n");
    }
// initializing the number vis points to be read    
    smabuffer.nused=0;
    free(cdh);
    free(sph);
    rewind(fpin[3]);
    rewind(fpin[2]);
    rewind(fpin[5]);
    sch = (struct sch_def **) malloc(nsets[0]*sizeof( struct sch_def *));
    for (set=0; set<nsets[0];set++) {
      sch[set] = (struct sch_def *)malloc(sizeof(struct sch_def ));
      if (sch[set] == NULL ){
	printf("ERROR: Memory allocation for sch failed for %d bytes\n",
	       nsets[0]*sizeof(struct sch_def));
	exit(-1);
      }
    }
    
    
// Need an array to hold the starting byte for each integration of data 
// allocate the memory for data_start_pos
    data_start_pos = (long int*)malloc(nsets[0]*sizeof(long int));
    inhset=smabuffer.scanskip;
    inhset=0;
    for (set=inhset;
      set<smabuffer.scanskip+smabuffer.scanproc; set++) {
      data_start_pos[set] = ftell(fpin[5]);
      *sch[set] = *(sch_head_read(fpin[5]));
      if (SWAP_ENDIAN) {
	sch[set]=swap_sch(sch[set]);
      }
      i = fseek(fpin[5],(long int)sch[set]->nbyt,SEEK_CUR);
    }
// rewind vis data file
    rewind(fpin[5]);  
// initilize the handles for baseline, spectral, integration 
    blhset  = -1;
    sphset  =  0;
    readSet =  1;
    numberBaselines = uvwbsln[smabuffer.scanskip]->n_bls;
// print the information on # of baselines,
//                          # of spectral windows,
//                          # of sidebands,
//                          # of receivers to be processed 
 printf("#Baselines=%d #Spectra=%d  #Sidebands=%d #Receivers=%d\n",
	   numberBaselines/2, numberSpectra-1, numberSidebands, 
	   numberRxif);
    
    if(smabuffer.dobary==1) printf("Compute radial velocity wrt barycenter\n");
    if(smabuffer.dolsr==1) printf("Compute radial velocity wrt LSR\n");
    sphSizeBuffer = numberSpectra*numberBaselines; 
//sphSizeBuffer: a number of total spectral sets for usb and lsb together 
//               in each of the integration sets
    firstsp = sphset;
    firstbsl= blhset;
    printf("start to read vis data!!!\n");
// initialize the vis point handle
    ipnt=1;
// assign the start integration set to be processed
  if(smabuffer.scanskip!=0) inhset=smabuffer.scanskip-1;
  if(smabuffer.scanproc!=0&&(nsets[0]>(smabuffer.scanskip+smabuffer.scanproc)))
     nsets[0]=smabuffer.scanskip+smabuffer.scanproc;
// start the processing loop
    while(inhset<(nsets[0]-1)) {
// progressing the integration set
      inhset++;
      visSMAscan.blockID.ints = inh[inhset]->ints;
      visSMAscan.blockID.inhid = inh[inhset]->inhid;
      
      visSMAscan.blockID.sourID = inh[inhset]->souid;
      visSMAscan.blockID.freqID = smaFreq[0].freqid;
      visSMAscan.time.UTCtime = jday+inh[inhset]->dhrs/24.000; /*hrs*/
// loading smabuffer 
      smabuffer.currentscan=inhset-smabuffer.scanskip;
      smabuffer.time = visSMAscan.time.UTCtime;
// handle source information 
      sourceID = visSMAscan.blockID.sourID;
      smabuffer.obsra = multisour[sourceID].ra_app;
      smabuffer.obsdec = multisour[sourceID].dec_app;
      smabuffer.ra = multisour[sourceID].ra;
      smabuffer.dec = multisour[sourceID].dec;
      {
// the source information was used to be obtained 
// from the engineer file.      
//double HA;
//double LST;
//double ra_apparent;
//double delLST;
// HA = LST - ra_apparent
// if(inh[inhset]->inhid==smaEngdata[inhset]->inhid) 
// { smabuffer.lst = smaEngdata[inhset]->lst*DPI/12.0;
// LST= (double) inh[inhset]->ha*DPI/12.0 + smabuffer.obsra;
//  delLST=(smabuffer.lst-LST)*12.0/DPI*3600;
// in unit of hr in mir format; in unit of radian in miriad format 
// printf("%d lst engLst %f %f delLST=%f obsra=%f %f\n", 
// inhset, smabuffer.lst, LST, delLST, smabuffer.obsra, smabuffer.ra);
//  for (i=0; i<smabuffer.nants; i++) {
//       smabuffer.el[i] = smaEngdata[inhset]->el[i+1];
//       smabuffer.az[i] = smaEngdata[inhset]->az[i+1];
//       smabuffer.tsys[i] = smaEngdata[inhset]->tsys[i+1];
// printf("%d %f %f --- %f %f\n",i, smabuffer.el[i], inh[inhset]->el,
//                              smabuffer.az[i], inh[inhset]->az);     
//                                    }
//}
//
// calculate lst
    if (smabuffer.doeng!=1) {
	smabuffer.lst =(double) inh[inhset]->ha*DPI/12.0 + smabuffer.obsra;
      } else {
	smabuffer.lst = smaEngdata[inhset]->lst*DPI/12.0;
      }
//loading el az and tsys to smabuffer
      for (i=0; i<smabuffer.nants; i++) {
// mir inh file gives the mean el and mane az
	smabuffer.el[i] = inh[inhset]->el;
	smabuffer.az[i] = inh[inhset]->az;
	if (smabuffer.doeng!=1) {
	  smabuffer.tsys[i]=atsys[inhset]->tssb[i+1];
	} else {
	  smabuffer.tsys[i] = smaEngdata[inhset]->tsys[i+1];
	}
        }
        }
        
// write source to uvfile 
        if(strncmp(multisour[sourceID].name,skipped,8)!=0)
        if(((strncmp(multisour[sourceID].name,target,6)!=0)&&
        (strncmp(multisour[sourceID].name,unknown,7)!=0))||
          smabuffer.noskip==1) {
	char sour[9];
	strncpy(sour, multisour[sourceID].name, 9);
// uvputvra_c(tno,"source",multisour[sourceID].name);
        uvputvra_c(tno,"source", sour);
	uvputvrd_c(tno,"ra",&(smabuffer.ra),1);
	uvputvrd_c(tno,"dec",&(smabuffer.dec),1);
// store the true pointing position 
	rar = inh[inhset]->rar;
//	if (rar!=smabuffer.ra) 
	  uvputvrd_c(tno,"pntra", &rar, 1);
	decr = inh[inhset]->decr;
//	if (decr!=smabuffer.dec)
        uvputvrd_c(tno,"pntdec",&decr,1);
	uvputvrd_c(tno,"obsra",&(smabuffer.obsra),1);
	uvputvrd_c(tno,"obsdec",&(smabuffer.obsdec),1);
	uvputvra_c(tno,"calcode",multisour[sourceID].calcode);
	uvputvri_c(tno,"sourid", &sourceID, 1);
      }
// configure the frequency for each of the integration set
      smabuffer.newfreq =1;
      smabuffer.veldop = (float) spn[inhset]->veldop;
      smabuffer.smaveldop = (float) spn[inhset]->smaveldop;
// calculate radial velocity to replace the on-line value
      { 
	short dolsr;
	double time   = smabuffer.time;
	double raapp  = smabuffer.obsra;
	double decapp = smabuffer.obsdec;
	double raepo  = smabuffer.ra;
	double decepo = smabuffer.dec;
	double lst    = smabuffer.lst;
	double lat    = smabuffer.lat;
//
// recalculate the radial velocity
//
	if(smabuffer.dolsr==1) {
  dolsr   = 1;
smabuffer.veldop=(float)velrad(dolsr,time, raapp,decapp,raepo,decepo,lst,lat);
	}
	if(smabuffer.dobary==1) {
	  dolsr   = -1;
smabuffer.veldop = (float)velrad(dolsr,time, raapp,decapp,raepo,decepo,lst,lat);
	}
        }
// get vsource which is source radial velocity w.r.t the LSR
// it is difficult to decode vsource from the MIR data
// unless given enough information in the Doppler tracking
// (sideband, chunk, and channel) from the users.  
//         smabuffer.vsource = spn[inhset]->vel[12];
// checkup skipped sp chunks
    for(i=1;i<SMIF+1;i++) {
        miriadsp[i]=0;
                      }
        i0=0;
        for(i=1;i<SMIF+1;i++) {
        if (spn[inhset]->nch[i][rxlod]!=0) {
// miriadsp: an integer array; if the chunk i
// is empty, miriadsp[i]=0, otherwise miriadsp[i]=i. 
        miriadsp[i]= i;
          } }
//
// now handle the frequency configuration for
// each of the integration sets
//
      for(i=1;i<smaCorr.n_chunk+1; i++) {
//         for(i=1;i<20;i++) {
// the reference channel is the first channel in each chunk in miriad
// the reference channel is the center (nch/2+0.5) in each chunk in MIR
// conversion => nch/2+0.5 - 1 = nch-0.5
         if(smabuffer.spskip[0]!=-1) {
         if(smabuffer.doChunkOrder==1) {
//
// spcode[i]:
// for the last three blocks (4,5,6), the chunk order is normal in each block.
// for the first three blocks (1,2,3), the chunk order is reversed 
// in each block.
// 4 3 2 1, 8 7 6 5, 12 11 10 9, 13 14 15 16, 17 18 19 20, 21 22 23 24 
// reverse chunk order for the frequency only using 
// conversion code frcode defined early:
//
	  smabuffer.sfreq[frcode[i]-1] = spn[inhset]->fsky[i]
      - spn[inhset]->fres[i]/1000.0*
	    (spn[inhset]->nch[i][rxlod]/2-0.5);
   	                                } else {
	  smabuffer.sfreq[spcode[i]-1] = spn[inhset]->fsky[i]
      - spn[inhset]->fres[i]/1000.0*
	    (spn[inhset]->nch[i][rxlod]/2-0.5);
	                                       }
//
// handle rest frequency
//
	if(smabuffer.dorfreq==1) 
        smabuffer.restfreq[spcode[i]-1] = spn[inhset]->rfreq[i];
//
// parsing the resampling option
//
        if(smabuffer.rsnchan<0) {
        smabuffer.sdf[spcode[i]-1]   = spn[inhset]->fres[i]/1000.0;
        smabuffer.nfreq[spcode[i]-1] = spn[inhset]->nch[i][rxlod];
	                        } else {
//
// re-sample the channel
//
        smabuffer.sdf[spcode[i]-1]   = spn[inhset]->fres[i]/1000.0*
  	spn[inhset]->nch[i][rxlod]/smabuffer.rsnchan;
	smabuffer.nfreq[spcode[i]-1] = smabuffer.rsnchan;
                                     	}
	smabuffer.basefreq = spn[inhset]->basefreq;
           } else {
//
// take the frequency configuration from the 1st integration set
// assuming no frequency configuration change during the 
// observing track.
//
          inhset1st=smabuffer.scanskip;
// parsing the option for chunk order reversing
          if(smabuffer.doChunkOrder==1) {
// reversing the chunk order
          smabuffer.sfreq[frcode[i]-1] = spn[inhset1st]->fsky[i]
            - spn[inhset1st]->fres[i]/1000.0*
            (spn[inhset1st]->nch[i][rxlod]/2-0.5);
        } else {
          smabuffer.sfreq[spcode[i]-1] = spn[inhset1st]->fsky[i]
            - spn[inhset1st]->fres[i]/1000.0*
          (spn[inhset1st]->nch[i][rxlod]/2-0.5);
        }
// processing the rest frequency
        if(smabuffer.dorfreq==1)
        smabuffer.restfreq[spcode[i]-1] = spn[inhset1st]->rfreq[i];
// parsing the resampling option
        if(smabuffer.rsnchan<0) {
        smabuffer.sdf[spcode[i]-1] = spn[inhset1st]->fres[i]/1000.0;
        smabuffer.nfreq[spcode[i]-1] = spn[inhset1st]->nch[i][rxlod];
        } else {
// re-sample the channel
        smabuffer.sdf[spcode[i]-1] = spn[inhset1st]->fres[i]/1000.0*
        spn[inhset1st]->nch[i][rxlod]/smabuffer.rsnchan;
        smabuffer.nfreq[spcode[i]-1] = smabuffer.rsnchan;
        }
// assign the base frequency
        smabuffer.basefreq = spn[inhset1st]->basefreq;
          }
// assign the rest of the header parameters
	smabuffer.bchan[spcode[i]-1]=1;
	smabuffer.nstoke[spcode[i]-1]=4;
	smabuffer.edge[spcode[i]-1]=0;
	smabuffer.nbin[spcode[i]-1]=1;
// if 2003 data skip the spcode examination
        if(smabuffer.spskip[0]!=-2003)  
// check the spectral window skipping in MIR data
        if(smabuffer.highrspectra!=1) {
        if(smabuffer.spskip[0]!=-1) {
//
// when only one skipping gap occured
//
   if(spcode[i]!=spn[inhset]->iband[i]) {
 printf("\n");
   if(smabuffer.spskip[0]==0) {
 printf("Spotted skipping in spectral chunks starting at spcode=%d iband=%d\n", spcode[i], spn[inhset]->iband[i]);
 printf("Try smalod with keyword spskip=%d,%d again.\n",
    spn[inhset]->iband[i], spcode[i]-spn[inhset]->iband[i]);
       } else {
 printf("The skipping parameter spskip =%d,%d is inconsistent with\n",
           spcode[i], spn[inhset]->iband[i]);
 printf("the spectral chunks skipped in the MIR data!\n");
        }
   bug_c( 'f', "spcode must match with iband!\n");
      }
      }
      }
      }
//
// handling multiple chunk skip gaps
//
// reconfigure 
   if(smabuffer.highrspectra==1) 
   for(i=1;i<SMIF+1; i++) {
   if(spn[inhset]->nch[i][rxlod]==0) {
// for those empty chunks, padding frequency
// header parameters with artifical numbers.
   smabuffer.sdf[spcode[i]-1]=0.104*pow(-1,smabuffer.sb+1);       
   smabuffer.sfreq[spcode[i]-1]=
   spn[inhset]->basefreq+0.084*(spcode[i]-11.5)*pow(-1,smabuffer.sb+1);
   smabuffer.restfreq[spcode[i]-1]=spn[inhset]->basefreq;
   smabuffer.nfreq[spcode[i]-1]=1;
                      }
        }       

   for(j=0; j < numberBaselines; j++) {
	{
	  sblpnt=j;
	  blhset++;
	  visSMAscan.uvblnID = uvwbsln[inhset]->uvwID[j].blcode;
	  visSMAscan.blockID.sbid = uvwbsln[inhset]->uvwID[j].isb; 
	  visSMAscan.blockID.polid = uvwbsln[inhset]->uvwID[j].ipol;
	  sbpnt = visSMAscan.blockID.sbid;
	  rxpnt = uvwbsln[inhset]->uvwID[j].irec;
   if(smabuffer.rxif==uvwbsln[inhset]->uvwID[j].irec||smabuffer.rxif==-1) 
        {
   switch(sbpnt) {
      case 0: blpnt=uvwbsln[inhset]->uvwID[j].blsid;
	      phaseSign=1;
// if required by users
// if data observed earlier than JD 2453488.5 (2005 4 28)
      if(smabuffer.doConjugate==-1||jday<2453488.5) phaseSign=-1;
      if(smabuffer.sb==0) {
smabuffer.u[blpnt] = uvwbsln[inhset]->uvwID[j].u/smabuffer.basefreq*1000.;
smabuffer.v[blpnt] = uvwbsln[inhset]->uvwID[j].v/smabuffer.basefreq*1000.;
smabuffer.w[blpnt] = uvwbsln[inhset]->uvwID[j].w/smabuffer.basefreq*1000.;
	      }
      break;
      case 1: blpnt=uvwbsln[inhset]->uvwID[j].blsid;
      phaseSign= 1;
      if(smabuffer.sb==1) {
   smabuffer.u[blpnt] = uvwbsln[inhset]->uvwID[j].u/smabuffer.basefreq*1000.;
   smabuffer.v[blpnt] = uvwbsln[inhset]->uvwID[j].v/smabuffer.basefreq*1000.;
   smabuffer.w[blpnt] = uvwbsln[inhset]->uvwID[j].w/smabuffer.basefreq*1000.;
	      }
      break;
	    }
	  }
      flush=1;
      if(smabuffer.nopol==1) visSMAscan.blockID.polid=-5;
	  switch(visSMAscan.blockID.polid)  {
	  case  0: polpnt=0; break;
	  case -1: polpnt=1; break;
	  case -2: polpnt=2; break;
	  case -3: polpnt=3; break;
	  case -4: polpnt=4; break;
	  case -5: polpnt=1; break;
	  case -6: polpnt=2; break;
	  case -7: polpnt=3; break;
	  case -8: polpnt=4; break;         }
// loading smabuffer uvw
  smabuffer.blcode[blpnt] = (float) visSMAscan.uvblnID;
// read sph for a complete spectral records assuming that the correlator
// configuration is not changed during the observation 
  if(readSet<= 1&&j==0) { 
    free(sph);
    sph = (struct sph_def **) malloc(sphSizeBuffer*sizeof( struct sph_def *));
    for (set=0; set < sphSizeBuffer; set++) {
    sph[set] = (struct sph_def *)malloc(sizeof(struct sph_def ));
    if (sph[set] == NULL ){
      printf("ERROR: Memory allocation for sph failed for %d bytes\n",
      sphSizeBuffer*sizeof(struct sph_def));
      exit(-1);
	      }
	    }
    for (set=0; set< sphSizeBuffer; set++) {
      *sph[set] = *(sph_read(fpin[2]));
      if (SWAP_ENDIAN) {
       sph[set] =  swap_sph(sph[set]);
	      }
	    }
// set the dataoff for the first spectrum is 0	    
	    sph[0]->dataoff = 0;
               }

// for single rx case
            for (i=0;i<smaCorr.n_chunk+1;i++) 
            {
            sph[i]->nch = spn[inhset]->nch[i][0];
            }

// separate frequency configuration for rx1 and rx2 in dual rx case
          if(smaCorr.no_rxif==2)
          for (i=0;i<smaCorr.n_chunk+1;i++) {
           if(rxpnt==smabuffer.rx1) 
              sph[i]->nch = spn[inhset]->nch[i][0];
           if(rxpnt==smabuffer.rx2)
              sph[i]->nch = spn[inhset]->nch[i][1];
                }

// The data for this spectrum consists of a 5 short int record header 
// and the data  which is a short for each real and a short for 
// each imag vis 
// for the ch0 
	  datalength = 5 + 2*sph[0]->nch;
// Get some memory to hold the data 
	  shortdata = (short int* )malloc(datalength*sizeof(short int));
// Here is the starting byte for the data for this spectrum 
// The 16 is to skip over the data integration header which precedes 
// all the records. We already read this data integration header 
// with sch_head_read. 
// inhset has the integration id 
	  bytepos = 16 + data_start_pos[inhset] + sph[0]->dataoff; 
	  bytepos = bytepos * (sblpnt+1);
// Move to this position in the file and read the data by
// skipping the first a few integration 
	  if(j<1) fseek(fpin[5],bytepos,SEEK_SET);
	  nbytes = sch_data_read(fpin[5],datalength,shortdata);
	  if (SWAP_ENDIAN) {
	    shortdata=swap_sch_data(shortdata, datalength);
	  }
// take a look at the record header 
// intTime in sec, integration time
	  visSMAscan.time.intTime = (double) (shortdata[0]/10.);
// loading smabuffer intTime 
	  smabuffer.inttime[j]= visSMAscan.time.intTime;
	  smabuffer.inttim = visSMAscan.time.intTime;
// There is a different scale factor for each record (spectrum) 
	  scale = shortdata[4];
// Now the data: There is only one channel for the continuum
// which can recalculated and are not stored in Miriad data
	  visSMAscan.bsline.continuum.real = 
          pow(2.,(double)scale)*shortdata[5];
	  visSMAscan.bsline.continuum.imag = 
	  pow(2.,(double)scale)*(-shortdata[6]*phaseSign);
	  free(shortdata);
	  sphset++;  
// The file is positioned at the record header for the next spectrum. 
// There is no 16 byte integration header between records 
	  ifpnt=0;
//	  for (kk=1;kk<numberSpectra; kk++) {
          for (kk=1;kk<smaCorr.n_chunk+1; kk++) {
// assign sp pointr:
	    ifpnt = spcode[kk]-1;
	   datalength = 5 + 2*sph[kk]->nch;
// 5 short in header of each spectral visibility data 
// 2 shorts for real amd imaginary 
// sph[sphset]->nch number of channel in each spectral set
// Get some memory to hold the data 
// memory size = datalength*sizeof(short)) 
//    if(smabuffer.highrspectra==1&&miriadsp[kk]==0) goto chunkskip;
//    shortdata = (short int* )malloc(datalength*sizeof(short int));
//	    nbytes = sch_data_read(fpin[5],datalength,shortdata);
//	    if (SWAP_ENDIAN) {
//	      shortdata=swap_sch_data(shortdata, datalength);
//	    }
     
// There is a different scale factor for each record (spectrum) 
//	    scale = shortdata[4];
// update polcode and pnt to smabuffer  
//
// swapping pol state (RR<->LL, RL<->LR) before 2005-6-10 or 
//                     -1<->-2,-3<->-4
// julian day 2453531.5
// 
            if (jday <2453531.5 && smabuffer.circular == 1) {
            if(visSMAscan.blockID.polid == -1) 
                 smabuffer.polcode[ifpnt][polpnt][blpnt]=-2;
            if(visSMAscan.blockID.polid == -2)
                 smabuffer.polcode[ifpnt][polpnt][blpnt]=-1;
            if(visSMAscan.blockID.polid == -3)
                 smabuffer.polcode[ifpnt][polpnt][blpnt]=-4;
            if(visSMAscan.blockID.polid == -4)
                 smabuffer.polcode[ifpnt][polpnt][blpnt]=-3;
               } else {
	    smabuffer.polcode[ifpnt][polpnt][blpnt]=visSMAscan.blockID.polid;
                   }
	    /* smabuffer.pnt[ifpnt][polpnt][blpnt][0] = ipnt;*/
	    smabuffer.pnt[ifpnt][polpnt][blpnt][sbpnt][rxpnt] = ipnt;
//         printf("kk spcode ifpnt smabuffer.pnt %d %d %d %d\n",kk,spcode[kk],
//            ifpnt,smabuffer.pnt[ifpnt][polpnt][blpnt][sbpnt][rxpnt]);
            if(wts[inhset]->wt[ifpnt][polpnt][blpnt][sbpnt][rxpnt] < 0.)
            smabuffer.flag[ifpnt][polpnt][blpnt][sbpnt][rxpnt] = 0;
            if(smabuffer.highrspectra==1&&miriadsp[kk]==0)
            smabuffer.flag[ifpnt][polpnt][blpnt][sbpnt][rxpnt] = 0;

    if(smabuffer.highrspectra==1&&miriadsp[kk]==0) goto chunkskip;
    shortdata = (short int* )malloc(datalength*sizeof(short int));
          nbytes = sch_data_read(fpin[5],datalength,shortdata);
          if (SWAP_ENDIAN) {
            shortdata=swap_sch_data(shortdata, datalength);
          }

// There is a different scale factor for each record (spectrum)
          scale = shortdata[4];

	    /* Now the channel data.  */
	    /* Make pseudo continuum */
	    avenchan = 0;
	    avereal  = 0.;
	    aveimag  = 0.;
	    for(i=0;i<sph[kk]->nch;i++){
	      if (smabuffer.rsnchan> 0) {
              if(sph[kk]->nch < smabuffer.rsnchan) {
 printf("Error: rsnchan=%d is greater than %d, the number of channels in a chunk,\n",
            smabuffer.rsnchan, sph[kk]->nch);
            printf("       redo it with a smaller rsnchan.\n");
            exit(-1);
                                                    }
// average the channel to the desired resolution 
  avereal = avereal+(float)pow(2.,(double)scale)*shortdata[5+2*i];
// convert ovro sign convention to miriad
  aveimag = aveimag+(float)pow(2.,(double)scale)*(-shortdata[6+2*i]*phaseSign);
	  avenchan++;
            if(avenchan==(sph[kk]->nch/smabuffer.rsnchan)) {
		  avereal = avereal/avenchan;
		  aveimag = aveimag/avenchan;
		  smabuffer.data[ipnt].real=avereal;
		  smabuffer.data[ipnt].imag=aveimag;
		  ipnt++;
		  avenchan = 0;
                  avereal  = 0.;
                  aveimag  = 0.;
		}
	      } else {
		/* loading the original vis with no average */
		// convert ovro sign convention to miriad
	smabuffer.data[ipnt].real=(float)pow(2.,(double)scale)*shortdata[5+2*i];
	smabuffer.data[ipnt].imag=
		  (float)pow(2.,(double)scale)*(-shortdata[6+2*i]*phaseSign);
		ipnt++;    }
	    }
            free(shortdata);
	    chunkskip:
    if(smabuffer.highrspectra==1&&miriadsp[kk]==0) {
            smabuffer.data[ipnt].real=0.0001;
            smabuffer.data[ipnt].imag=0.;
            ipnt++;  
                                                   }
	    sphset++;  /* count for each spectral chunk */
	  }
	    firstsp = sphset;
        }
      }
      readSet++;
      smabuffer.nused=ipnt;
// re-initialize vis point for next integration 
      ipnt=1;
      if(flush==1) {
      if (fmod((readSet-1), 100.)<0.5||(readSet-1)==1)
  printf("set=%4d ints=%4d inhid=%4d time(JulianDay)=%9.5f int=% 4.1f \n",
		 readSet-1,
		 visSMAscan.blockID.ints,
		 visSMAscan.blockID.inhid,
		 visSMAscan.time.UTCtime,
		 visSMAscan.time.intTime);
// call rspokeflshsma_c to store databuffer to uvfile 
	kstat = -1;
	*kst = (char *)&kstat;
        if(smabuffer.noskip!=1) {
        if(strncmp(multisour[sourceID].name,skipped,8)==0)
 printf("Warnning: one scan is skipped at %9.5f due to insufficient source information.\n",visSMAscan.time.UTCtime);
        if(strncmp(multisour[sourceID].name,skipped,8)!=0) {
        if((strncmp(multisour[sourceID].name,target,6)!=0)&&
           (strncmp(multisour[sourceID].name,unknown,7)!=0)) {
	  rspokeflshsma_c(kst);
	                                    } else {
                                            ntarget++;
	                                         }      }

                                } else {
         rspokeflshsma_c(kst);
                      }
      }
    }
    printf("set=%4d ints=%4d inhid=%4d time(JulianDay)=%9.5f int=% 4.1f \n",
	   readSet-1, 
	   visSMAscan.blockID.ints,
	   visSMAscan.blockID.inhid,
	   visSMAscan.time.UTCtime,
	   visSMAscan.time.intTime);
    printf("skipped %d integration scans on `target&unknown'\n",ntarget);

    avenchan=smabuffer.rsnchan;
         for (i=0; i<ireset+1; i++) {
             if(i>0) smabuffer.scanskip = reset_id[i];
    if (smabuffer.rsnchan>0) {
int nnpadding;
      printf("converted vis spectra from the original correlator configuration\n");
      printf("to low and uniform resolution spectra:\n");
      printf("(starting at integrations=%d)\n", smabuffer.scanskip);
      printf("         input     output\n");
      for (kk=1; kk<numberSpectra; kk++) {
 if(smabuffer.highrspectra==1&&spn[smabuffer.scanskip]->nch[kk][rxlod]==0)
{
       nnpadding = 1-avenchan;
printf("warning: each of the empty chunks is padded with one flagged channel.\n");

} else {
       nnpadding = 0;
}

        if(smabuffer.spskip[0]==0) {
	printf("  s%02d     %4d  =>  s%02d %4d\n",kk, 
        spn[smabuffer.scanskip]->nch[kk][rxlod], kk, avenchan+nnpadding);
          } else {
        if(kk < smabuffer.spskip[0]) {
        printf("  s%02d     %4d  =>  s%02d %4d\n",kk,
        spn[smabuffer.scanskip]->nch[kk][rxlod], kk, avenchan+nnpadding);
          } else {
        printf("  s%02d     %4d  =>  s%02d %4d\n",kk+smabuffer.spskip[1],
        spn[smabuffer.scanskip]->nch[kk][rxlod], kk, avenchan+nnpadding);
           }
           }

      }
    } else {
int npadding=0;
      printf("vis spectra from the original correlator configuration: \n");
      printf("(starting at integrations=%d)\n", smabuffer.scanskip);
      printf("         input     output\n");
      for (kk=1; kk<numberSpectra; kk++) {
 if(smabuffer.highrspectra==1&&spn[smabuffer.scanskip]->nch[kk][rxlod]==0)
{ 
       npadding = 1; 
printf("warning: each of the empty chunks is padded with one flagged channel.\n")
;

} else { 
       npadding = 0;
}   
       if(smabuffer.spskip[0]==0) {
        printf("  s%02d     %4d  =>  s%02d %4d\n",kk, 
        spn[smabuffer.scanskip]->nch[kk][rxlod], kk, 
        spn[smabuffer.scanskip]->nch[kk][rxlod]+npadding);
         } else {
        if(kk < smabuffer.spskip[0]) {
        printf("  s%02d     %4d  =>  s%02d %4d\n",kk,
        spn[smabuffer.scanskip]->nch[kk][rxlod], kk,
        spn[smabuffer.scanskip]->nch[kk][rxlod]+npadding); 
        } else {
        printf("  s%02d     %4d  =>  s%02d %4d\n",
        kk+smabuffer.spskip[1],
        spn[smabuffer.scanskip]->nch[kk][rxlod], kk,
        spn[smabuffer.scanskip]->nch[kk][rxlod]+npadding);
               }
         }
      }
    }
     }
    if(smabuffer.spskip[0]!=-2003) {
    if(smabuffer.spskip[0]!=0&&smabuffer.spskip[0]!=-1)
     printf("The MIR s%02d - s%02d contain no data and are skipped!\n",
        smabuffer.spskip[0],
        smabuffer.spskip[0]+smabuffer.spskip[1]-1);}
     else {
     printf("The MIR s%d contains incompleted-correlator data!\n",
       -smabuffer.spskip[0]);
          } 
    printf("Done with data conversion from mir to miriad!\n");
  }
  /* ---------------------------------------------------------------------- */
  endTime = time(NULL);
  trueTime = difftime(endTime, startTime);
  if(jstat==0) fprintf(stderr,
		       "Real time used =%f sec.\n",
		       trueTime);
  jstat=0;
  return jstat;
} /* end of main */



/* This function reads the integration header */
struct inh_def *inh_read(FILE * fpinh)
{
  int nbytes;   /* counts number of bytes written */
  int nobj;     /* the number of objects written by each write */
  
  struct inh_def inh;
  struct inh_def *inhptr;
 
  nbytes = 0;
  nobj = 0;

  /* read first data in set and check for end of file */
  
  nobj += fread(&inh.conid,sizeof(inh.conid),1,fpinh);
  if (nobj == 0) {
    printf("Unexpected end of file in_read\n");
    exit(-1);
  }
  nbytes += sizeof(inh.conid);
  nobj += fread(&inh.icocd,sizeof(inh.icocd),1,fpinh);
  nbytes += sizeof(inh.icocd);
  nobj += fread(&inh.traid,sizeof(inh.traid),1,fpinh);
  nbytes += sizeof(inh.traid);
  nobj += fread(&inh.inhid,sizeof(inh.inhid),1,fpinh);
  nbytes += sizeof(inh.inhid);
  nobj += fread(&inh.ints,sizeof(inh.ints),1,fpinh);
  nbytes += sizeof(inh.ints);
  nobj += fread(&inh.itq,sizeof(inh.itq),1,fpinh);
  nbytes += sizeof(inh.itq);
  nobj += fread(&inh.az,sizeof(inh.az),1,fpinh);
  nbytes += sizeof(inh.az);
  nobj += fread(&inh.el,sizeof(inh.el),1,fpinh);
  nbytes += sizeof(inh.el);
  nobj += fread(&inh.ha,sizeof(inh.ha),1,fpinh);
  nbytes += sizeof(inh.ha);
  nobj += fread(&inh.iut,sizeof(inh.iut),1,fpinh);
  nbytes += sizeof(inh.iut);
  nobj += fread(&inh.iref_time,sizeof(inh.iref_time),1,fpinh);
  nbytes += sizeof(inh.iref_time);
  nobj += fread(&inh.dhrs,sizeof(inh.dhrs),1,fpinh);
  nbytes += sizeof(inh.dhrs);
  nobj += fread(&inh.vc,sizeof(inh.vc),1,fpinh);
  nbytes += sizeof(inh.vc);
  nobj += fread(&inh.ivctype,sizeof(inh.ivctype),1,fpinh);
  nbytes += sizeof(inh.ivctype);
  nobj += fread(&inh.sx,sizeof(inh.sx),1,fpinh);
  nbytes += sizeof(inh.sx);
  nobj += fread(&inh.sy,sizeof(inh.sy),1,fpinh);
  nbytes += sizeof(inh.sy);
  nobj += fread(&inh.sz,sizeof(inh.sz),1,fpinh);
  nbytes += sizeof(inh.sz);
  nobj += fread(&inh.rinteg,sizeof(inh.rinteg),1,fpinh);
  nbytes += sizeof(inh.rinteg);
  nobj += fread(&inh.proid,sizeof(inh.proid),1,fpinh);
  nbytes += sizeof(inh.proid);
  nobj += fread(&inh.souid,sizeof(inh.souid),1,fpinh);
  nbytes += sizeof(inh.souid);
  nobj += fread(&inh.isource,sizeof(inh.isource),1,fpinh);
  nbytes += sizeof(inh.isource);
  nobj += fread(&inh.ipos,sizeof(inh.ipos),1,fpinh);
  nbytes += sizeof(inh.ipos);
  nobj += fread(&inh.offx,sizeof(inh.offx),1,fpinh);
  nbytes += sizeof(inh.offx);
  nobj += fread(&inh.offy,sizeof(inh.offy),1,fpinh);
  nbytes += sizeof(inh.offy);
  nobj += fread(&inh.iofftype,sizeof(inh.iofftype),1,fpinh);
  nbytes += sizeof(inh.iofftype);
  nobj += fread(&inh.ira,sizeof(inh.ira),1,fpinh);
  nbytes += sizeof(inh.ira);
  nobj += fread(&inh.idec,sizeof(inh.idec),1,fpinh);
  nbytes += sizeof(inh.idec);
  nobj += fread(&inh.rar,sizeof(inh.rar),1,fpinh);
  nbytes += sizeof(inh.rar);
  nobj += fread(&inh.decr,sizeof(inh.decr),1,fpinh);
  nbytes += sizeof(inh.decr);
  nobj += fread(&inh.epoch,sizeof(inh.epoch),1,fpinh);
  nbytes += sizeof(inh.epoch);
  nobj += fread(&inh.sflux,sizeof(inh.sflux),1,fpinh);
  nbytes += sizeof(inh.sflux);
  nobj += fread(&inh.size,sizeof(inh.size),1,fpinh);
  nbytes += sizeof(inh.size);
  inhptr = &inh;
  return inhptr;
  
} /* end of function inh_read */


/* This function reads one baseline header */
struct blh_def *blh_read(FILE * fpblh)
{
  int nbytes;   /* counts number of bytes written */
  int nobj;     /* the number of objects written by each write */
  struct blh_def blh;
  struct blh_def *blhptr;
 
  nbytes = 0;
  nobj = 0;
 
 
  nobj += fread(&blh.blhid,sizeof(blh.blhid),1,fpblh);
  if (nobj == 0) {
    printf("Unexpected end of file bl_read\n");
    exit(-1);
  }
  nbytes += sizeof(blh.blhid);
  nobj += fread(&blh.inhid,sizeof(blh.inhid),1,fpblh);
  nbytes += sizeof(blh.inhid);
  nobj += fread(&blh.isb,sizeof(blh.isb),1,fpblh);
  nbytes += sizeof(blh.isb);
  nobj += fread(&blh.ipol,sizeof(blh.ipol),1,fpblh);
  nbytes += sizeof(blh.ipol);
  nobj += fread(&blh.pa,sizeof(blh.pa),1,fpblh);
  nbytes += sizeof(blh.pa);
  nobj += fread(&blh.iaq,sizeof(blh.iaq),1,fpblh);
  nbytes += sizeof(blh.iaq);
  nobj += fread(&blh.ibq,sizeof(blh.ibq),1,fpblh);
  nbytes += sizeof(blh.ibq);
  nobj += fread(&blh.icq,sizeof(blh.icq),1,fpblh);
  nbytes += sizeof(blh.icq);
  nobj += fread(&blh.ioq,sizeof(blh.ioq),1,fpblh);
  nbytes += sizeof(blh.ioq);
  nobj += fread(&blh.irec,sizeof(blh.irec),1,fpblh);
  nbytes += sizeof(blh.irec);
  nobj += fread(&blh.iifc,sizeof(blh.iifc),1,fpblh);
  nbytes += sizeof(blh.iifc);
  nobj += fread(&blh.u,sizeof(blh.u),1,fpblh);
  nbytes += sizeof(blh.u);
  nobj += fread(&blh.v,sizeof(blh.v),1,fpblh);
  nbytes += sizeof(blh.v);
  nobj += fread(&blh.w,sizeof(blh.w),1,fpblh);
  nbytes += sizeof(blh.w);
  nobj += fread(&blh.prbl,sizeof(blh.prbl),1,fpblh);
  nbytes += sizeof(blh.prbl);
  nobj += fread(&blh.angres,sizeof(blh.angres),1,fpblh);
  nbytes += sizeof(blh.angres);
  nobj += fread(&blh.vis,sizeof(blh.vis),1,fpblh);
  nbytes += sizeof(blh.vis);
  nobj += fread(&blh.coh,sizeof(blh.coh),1,fpblh);
  nbytes += sizeof(blh.coh);
  nobj += fread(&blh.sigcoh,sizeof(blh.sigcoh),1,fpblh);
  nbytes += sizeof(blh.sigcoh);
  nobj += fread(&blh.csnr,sizeof(blh.csnr),1,fpblh);
  nbytes += sizeof(blh.csnr);
  nobj += fread(&blh.vflux,sizeof(blh.vflux),1,fpblh);
  nbytes += sizeof(blh.vflux);
  nobj += fread(&blh.cnoise,sizeof(blh.cnoise),1,fpblh);
  nbytes += sizeof(blh.cnoise);
  nobj += fread(&blh.avedhrs,sizeof(blh.avedhrs),1,fpblh);
  nbytes += sizeof(blh.avedhrs);
  nobj += fread(&blh.ampave,sizeof(blh.ampave),1,fpblh);
  nbytes += sizeof(blh.ampave);
  nobj += fread(&blh.phaave,sizeof(blh.phaave),1,fpblh);
  nbytes += sizeof(blh.phaave);
  nobj += fread(&blh.tpvar,sizeof(blh.tpvar),1,fpblh);
  nbytes += sizeof(blh.tpvar);
  nobj += fread(&blh.blsid,sizeof(blh.blsid),1,fpblh);
  nbytes += sizeof(blh.blsid);
  nobj += fread(&blh.itel1,sizeof(blh.itel1),1,fpblh);
  nbytes += sizeof(blh.itel1);
  nobj += fread(&blh.itel2,sizeof(blh.itel2),1,fpblh);
  nbytes += sizeof(blh.itel2);
  nobj += fread(&blh.iblcd,sizeof(blh.iblcd),1,fpblh);
  nbytes += sizeof(blh.iblcd);
  nobj += fread(&blh.ble,sizeof(blh.ble),1,fpblh);
  nbytes += sizeof(blh.ble);
  nobj += fread(&blh.bln,sizeof(blh.bln),1,fpblh);
  nbytes += sizeof(blh.bln);
  nobj += fread(&blh.blu,sizeof(blh.blu),1,fpblh);
  nbytes += sizeof(blh.blu);
  nobj += fread(&blh.soid,sizeof(blh.soid),1,fpblh);
  nbytes += sizeof(blh.soid);
  blhptr = &blh;
  return blhptr;
  
} /* end of blh_read  */


unsigned long mfsize(FILE *fp)
{
  /* Optimization stuff */
  char temp[BUFSIZ];
  static const long DATALENGTH_MAX=SMA_LONG_MAX%2!=0?SMA_LONG_MAX-1:SMA_LONG_MAX;
  long datalength=DATALENGTH_MAX;
  
  unsigned long counter, fsize;
  fsize = 0;
  
  if (fp==NULL) {
    printf("null pointer\n");
    return (unsigned)NULL;
  }

  /* fseek() doesn't signal EOF so i use fread() to detect the end of file */
  for (fseek(fp, datalength-1, SEEK_SET); 
       datalength>0 && fread(temp, 1, 1, fp)==0; 
       fseek(fp, datalength-1, SEEK_SET)) datalength/=128;
  
  fseek(fp, 0, SEEK_SET);
  
  if (datalength==0 && fread(temp, 1, 1, fp)==0) {
    return fsize;
  } else if (datalength==0)
    datalength=BUFSIZ;

  fseek(fp, datalength-1, SEEK_SET);
  /* fseek() doesn't signal EOF so i use fread() to detect the end of file */
  for(counter=0; fread(temp, 1, 1, fp)!=0; ++counter)
    fseek(fp, datalength-1, SEEK_CUR);
  
  fseek(fp, 0, SEEK_SET);
  
  for( ; counter>0; --counter) {
    fseek(fp, datalength, SEEK_CUR);
    fsize += datalength;
  }

  do {
    fsize += datalength=fread(temp, 1, BUFSIZ, fp);
  } while(datalength!=0);
  
  fseek(fp, 0, SEEK_SET);
  
  return fsize;
}


struct sph_def *sph_read(FILE * fpsph)
{
  int nbytes;   /* counts number of bytes written */
  int nobj;     /* the number of objects written by each write */
 
  struct sph_def sph;
  struct sph_def *sphptr; 
  nbytes = 0;
  nobj = 0;
 
  nobj += fread(&sph.sphid,sizeof(sph.sphid),1,fpsph);
  if (nobj == 0) {
    printf("Unexpected end of file sp_read\n");
    exit(-1);
  }
  nbytes += sizeof(sph.sphid);
  nobj += fread(&sph.blhid,sizeof(sph.blhid),1,fpsph);
  nbytes += sizeof(sph.blhid);
  nobj += fread(&sph.inhid,sizeof(sph.inhid),1,fpsph);
  nbytes += sizeof(sph.inhid);
  nobj += fread(&sph.igq,sizeof(sph.igq),1,fpsph);
  nbytes += sizeof(sph.igq);
  nobj += fread(&sph.ipq,sizeof(sph.ipq),1,fpsph);
  nbytes += sizeof(sph.ipq);
  nobj += fread(&sph.iband,sizeof(sph.iband),1,fpsph);
  nbytes += sizeof(sph.iband);
  nobj += fread(&sph.ipstate,sizeof(sph.ipstate),1,fpsph);
  nbytes += sizeof(sph.ipstate);
  nobj += fread(&sph.tau0,sizeof(sph.tau0),1,fpsph);
  nbytes += sizeof(sph.tau0);
  nobj += fread(&sph.vel,sizeof(sph.vel),1,fpsph);
  nbytes += sizeof(sph.vel);
  nobj += fread(&sph.vres,sizeof(sph.vres),1,fpsph);
  nbytes += sizeof(sph.vres);
  nobj += fread(&sph.ivtype,sizeof(sph.ivtype),1,fpsph);
  nbytes += sizeof(sph.ivtype);
  nobj += fread(&sph.fsky,sizeof(sph.fsky),1,fpsph);
  nbytes += sizeof(sph.fsky);
  nobj += fread(&sph.fres,sizeof(sph.fres),1,fpsph);
  nbytes += sizeof(sph.fres);
  nobj += fread(&sph.tssb,sizeof(sph.tssb),1,fpsph);
  nbytes += sizeof(sph.tssb);
  nobj += fread(&sph.integ,sizeof(sph.integ),1,fpsph);
  nbytes += sizeof(sph.integ);
  nobj += fread(&sph.wt,sizeof(sph.wt),1,fpsph);
  nbytes += sizeof(sph.wt); 
  nobj += fread(&sph.itaper,sizeof(sph.itaper),1,fpsph);
  nbytes += sizeof(sph.itaper);
  nobj += fread(&sph.snoise,sizeof(sph.snoise),1,fpsph);
  nbytes += sizeof(sph.snoise);
  nobj += fread(&sph.nch,sizeof(sph.nch),1,fpsph);
  nbytes += sizeof(sph.nch);
  nobj += fread(&sph.nrec,sizeof(sph.nrec),1,fpsph);
  nbytes += sizeof(sph.nrec);
  nobj += fread(&sph.dataoff,sizeof(sph.dataoff),1,fpsph);
  nbytes += sizeof(sph.dataoff);
  nobj += fread(&sph.linid,sizeof(sph.linid),1,fpsph);
  nbytes += sizeof(sph.linid);
  nobj += fread(&sph.itrans,sizeof(sph.itrans),1,fpsph);
  nbytes += sizeof(sph.itrans);
  nobj += fread(&sph.rfreq,sizeof(sph.rfreq),1,fpsph);
  nbytes += sizeof(sph.rfreq);
  nobj += fread(&sph.pasid,sizeof(sph.pasid),1,fpsph);
  nbytes += sizeof(sph.pasid);
  nobj += fread(&sph.gaiidamp,sizeof(sph.gaiidamp),1,fpsph);
  nbytes += sizeof(sph.gaiidamp);
  nobj += fread(&sph.gaiidpha,sizeof(sph.gaiidpha),1,fpsph);
  nbytes += sizeof(sph.gaiidpha);
  nobj += fread(&sph.flcid,sizeof(sph.flcid),1,fpsph);
  nbytes += sizeof(sph.flcid);
  nobj += fread(&sph.atmid,sizeof(sph.atmid),1,fpsph);
  nbytes += sizeof(sph.atmid);
  sphptr = &sph;
  return sphptr;
  
} /* end of sph_write */

/* This function reads  the code or strings header */
struct codeh_def * cdh_read (FILE * fpcodeh)
{
  int nbytes;   /* counts number of bytes written */
  int nobj;     /* the number of objects written by each write */
 
  struct codeh_def codeh;
  struct codeh_def *codehptr;
  nbytes = 0;
  nobj = 0;

  nobj += fread(codeh.v_name,sizeof(codeh.v_name),1,fpcodeh);
  if (nobj == 0) {
    printf("Unexpected end of file cdh_read\n");
    exit(-1);
  }
  nbytes += sizeof(codeh.v_name);
  nobj += fread(&codeh.icode,sizeof(codeh.icode),1,fpcodeh);
  nbytes += sizeof(codeh.icode);
  nobj += fread(codeh.code,sizeof(codeh.code),1,fpcodeh);
  nbytes += sizeof(codeh.code);
  nobj += fread(&codeh.ncode,sizeof(codeh.ncode),1,fpcodeh);
  nbytes += sizeof(codeh.ncode);
  codehptr = &codeh;
  return codehptr;
  
} /* end of codeh_write */ 

/* This function reads the engineering data header */
struct ant_def * enh_read(FILE * fpeng)
{
  int nbytes;   /* counts number of bytes written */
  int nobj;     /* the number of objects written by each write */

  struct ant_def ant;
  struct ant_def *antptr; 
  nbytes = 0;
  nobj = 0;
 
  nobj += fread(&ant.antennaNumber,sizeof(ant.antennaNumber),1,fpeng);
  if (nobj == 0) {
    printf("Unexpected end of file enh_read\n");
    exit(-1);
  }
  nbytes += sizeof(ant.antennaNumber);
  nobj += fread(&ant.padNumber,sizeof(ant.padNumber),1,fpeng);
  nbytes += sizeof(ant.padNumber);
  nobj += fread(&ant.antennaStatus,sizeof(ant.antennaStatus),1,fpeng);
  nbytes += sizeof(ant.antennaStatus);
  nobj += fread(&ant.trackStatus,sizeof(ant.trackStatus),1,fpeng);
  nbytes += sizeof(ant.trackStatus);
  nobj += fread(&ant.commStatus,sizeof(ant.commStatus),1,fpeng);
  nbytes += sizeof(ant.commStatus);
  nobj += fread(&ant.inhid,sizeof(ant.inhid),1,fpeng);
  nbytes += sizeof(ant.inhid);
  nobj += fread(&ant.ints,sizeof(ant.ints),1,fpeng);
  nbytes += sizeof(ant.ints);
  nobj += fread(&ant.dhrs,sizeof(ant.dhrs),1,fpeng);
  nbytes += sizeof(ant.dhrs);
  nobj += fread(&ant.ha,sizeof(ant.ha),1,fpeng);
  nbytes += sizeof(ant.ha);
  nobj += fread(&ant.lst,sizeof(ant.lst),1,fpeng);
  nbytes += sizeof(ant.lst);
  nobj += fread(&ant.pmdaz,sizeof(ant.pmdaz),1,fpeng);
  nbytes += sizeof(ant.pmdaz);
  nobj += fread(&ant.pmdel,sizeof(ant.pmdel),1,fpeng);
  nbytes += sizeof(ant.pmdel);
  nobj += fread(&ant.tiltx,sizeof(ant.tiltx),1,fpeng);
  nbytes += sizeof(ant.tiltx);
  nobj += fread(&ant.tilty,sizeof(ant.tilty),1,fpeng);
  nbytes += sizeof(ant.tilty);
  nobj += fread(&ant.actual_az,sizeof(ant.actual_az),1,fpeng);
  nbytes += sizeof(ant.actual_az);
  nobj += fread(&ant.actual_el,sizeof(ant.actual_el),1,fpeng);
  nbytes += sizeof(ant.actual_el);
  nobj += fread(&ant.azoff,sizeof(ant.azoff),1,fpeng);
  nbytes += sizeof(ant.azoff);
  nobj += fread(&ant.eloff,sizeof(ant.eloff),1,fpeng);
  nbytes += sizeof(ant.eloff);
  nobj += fread(&ant.az_tracking_error,sizeof(ant.az_tracking_error),1,fpeng);
  nbytes += sizeof(ant.az_tracking_error);
  nobj += fread(&ant.el_tracking_error,sizeof(ant.el_tracking_error),1,fpeng);
  nbytes += sizeof(ant.el_tracking_error);
  nobj += fread(&ant.refraction,sizeof(ant.refraction),1,fpeng);
  nbytes += sizeof(ant.refraction);
  nobj += fread(&ant.chopper_x,sizeof(ant.chopper_x),1,fpeng);
  nbytes += sizeof(ant.chopper_x);
  nobj += fread(&ant.chopper_y,sizeof(ant.chopper_y),1,fpeng);
  nbytes += sizeof(ant.chopper_y);
  nobj += fread(&ant.chopper_z,sizeof(ant.chopper_z),1,fpeng);
  nbytes += sizeof(ant.chopper_z);
  nobj += fread(&ant.chopper_angle,sizeof(ant.chopper_angle),1,fpeng);
  nbytes += sizeof(ant.chopper_angle);
  nobj += fread(&ant.tsys,sizeof(ant.tsys),1,fpeng);
  nbytes += sizeof(ant.tsys);
  nobj += fread(&ant.ambient_load_temperature,sizeof(ant.ambient_load_temperature),1,fpeng);        
  nbytes += sizeof(ant.ambient_load_temperature);
  antptr = &ant;
  return antptr;
} 

struct sch_def * sch_head_read(FILE * fpsch)
{
  int nbytes;   /* counts number of bytes written */
  int nobj;     /* the number of objects written by each write */

  struct sch_def sch;
  struct sch_def *schptr; 
  nbytes = 0;
  nobj = 0;
 
  nobj += fread(&sch.inhid,sizeof(sch.inhid),1,fpsch);
  if (nobj == 0) {
    printf("Unexpected end of file sch_head_read\n");
    exit(-1);
  }
  nbytes += sizeof(sch.inhid);
  nobj += fread(sch.form,sizeof(sch.form),1,fpsch);
  nbytes += sizeof(sch.form);
  nobj += fread(&sch.nbyt,sizeof(sch.nbyt),1,fpsch);
  nbytes += sizeof(sch.nbyt);
  nobj += fread(&sch.nbyt_pack,sizeof(sch.nbyt_pack),1,fpsch);
  nbytes += sizeof(sch.nbyt_pack);
  schptr = &sch;
  return schptr;
  
} /* end of sch_write */ 

int sch_data_read(FILE * fpsch, long int datalength, short int * data)
{
  int nbytes;   /* counts number of bytes written */
  int nobj;     /* the number of objects written by each write */

  /*  short buff;*/
  nbytes = 0;
  nobj = 0;
  nobj += fread(data,datalength*sizeof(short int),1,fpsch);
  if (nobj == 0) {
    printf("The current scan (set=% 4d) is being read.\n",
	   smabuffer.currentscan); 
    bug_c('f',"Unexpected end of file shc_data_read. Using nscans\n to select a scan range and run smalod again.\n"); 
    exit(-1);
  }
  return nbytes;
  
} /* end of sch_data_read */

char *rar2c(double ra)
{ 
  static char rac[13];
  int hh, mm;
  float ss;
  hh = (int) (12.0/DPI*ra);
  mm = (int) ((12.0/DPI*ra-hh)*60.0);
  ss = (float) (((12.0/DPI*ra-hh)*60.0-mm)*60.0);
  sprintf(rac,"%02d:%02d:%07.4f", hh,mm,ss);
  rac[13]='\0';
  /*    printf("ra=%s\n", rac);
   */
  return &rac[0];      
}

char *decr2c(double dec)
{
  static char decc[15];
  int dd, am;
  float as;
  dd = (int)(180./DPI*dec);
  am = (int)((180./DPI*dec-dd)*60.0);
  as = (float)(((180./DPI*dec-dd)*60.0-am))*60.0;
  am = (int)fabs(am);
  as = (float)fabs(as);
  sprintf(decc,"% 3d:%02d:%07.4f", dd,am,as);
  decc[15]='\0';
  return &decc[0];
}
 
int spdecode (struct codeh_def *specCode[])
{ 
  int spid;
  char  cspid[13];
  char  prefix[1];
  cspid[13]='\0';
  memcpy(cspid, specCode[0]->code, 12);
  sscanf(cspid, "%1s%d", prefix, &spid);
  return spid;
}

float juliandate (struct codeh_def *refdate[])
{ 
  int i;
  int stat=0;   
  double jdate;
  char  ccaldate[13];
  static char *months[] = {"ill", "Jan","Feb","Mar","Apr","May","Jun","Jul", 
          "Aug","Sep","Oct","Nov","Dec"};
  char yc[4];
  char mc[3];
  int yi,mi,di;

  ccaldate[13]='\0';
  memcpy(ccaldate,refdate[0]->code, 12);
  sscanf(ccaldate, "%s%d%s%d", mc, &di,yc,&yi);
//  printf("Observing Date: %d %s %d\n", yi, mc, di);
  mi=0;
  for (i=1; i<13; i++){
    if (memcmp(mc,months[i], 3)==0) mi=i;
  }
  jdate = slaCldj (yi, mi, di, stat)+2400000.5;
  if(stat==1) {
     printf("bad year   (MJD not computed).");
     exit(-1);
              }
  if(stat==2) {
     printf("bad month   (MJD not computed).");
     exit(-1);
              }
  if(stat==3) {
     printf("bad day   (MJD not computed).");
     exit(-1);
              }
 printf("Observing Date: %d %s %d    Julian Date: %f\n", yi, mc, di, jdate);

  return jdate;
}


double slaCldj ( int iy, int im, int id,  int sj)
     /*
     **  - - - - - - - -
     **   s l a C l d j
     **  - - - - - - - -
     **
     **  Gregorian calendar to Modified Julian Date.
     **
     **  Given:
     **     iy,im,id     int    year, month, day in Gregorian calendar
     **
     **  Returned:
     **     mjd_rtn      double Modified Julian Date (JD-2400000.5) for 0 hrs
     **     sj           int    status:
     **                           0 = OK
     **                           1 = bad year   (MJD not computed)
     **                           2 = bad month  (MJD not computed)
     **                           3 = bad day    (MJD computed)
     **
     **  The year must be -4699 (i.e. 4700BC) or later.
     **
     **  The algorithm is derived from that of Hatcher 1984 (QJRAS 25, 53-55).
     **
     **  Last revision:   29 August 1994
     **
     **  Copyright P.T.Wallace.  All rights reserved.
     */
{
  long iyL, imL, mjd;
  double mjd_rtn;
  /* Month lengths in days */
  static int mtab[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  
  /* Validate year */
  if ( iy < -4699 ) { sj = 1; return 0.0; }
  
  /* Validate month */
  if ( ( im < 1 ) || ( im > 12 ) ) { sj = 2; return 0.0; }
  
  /* Allow for leap year */
  mtab[1] = ( ( ( iy % 4 ) == 0 ) &&
	      ( ( ( iy % 100 ) != 0 ) || ( ( iy % 400 ) == 0 ) ) ) ?
    29 : 28;
  
  /* Validate day */
  sj =  (( id < 1 || id > mtab[im-1] ) ? 3 : 0);
  
  /* Lengthen year and month numbers to avoid overflow */
  iyL = (long) iy;
  imL = (long) im;
  /* Perform the conversion */
  /*djm = (double)*/
  
  mjd =  ( ( 1461L * ( iyL - ( 12L - imL ) / 10L + 4712L ) ) / 4L
	   + ( 306L * ( ( imL + 9L ) % 12L ) + 5L ) / 10L
	   - ( 3L * ( ( iyL - ( 12L - imL ) / 10L + 4900L ) / 100L ) ) / 4L
	   + (long) id - 2399904L );
  mjd_rtn = (double) mjd;
  return mjd_rtn;
}

void precess(double jday1,double ra1,double dec1,double jday2,double *ra2pt,double *dec2pt)
{
  /* jhz 2004-7-23: based on a miriad fortran code, translate into c.
     c  A simple precession routine, to precess from one set of mean
     c  equatorial coordinates (RA,DEC), to another at a different epoch.
     c  This is accurate to order 0.3 arcsec over 50 years.
     c
     c  Reference:
     c    Explanatory Supplement to the Astronomical Almanac, 1993. p 105-106.
     c
     c  NOTE: This does not take account of atmospheric refraction,
     c  nutation, aberration nor gravitational deflection.
     c
     c  Input:
     c    jday1      Julian day of the known epoch.
     c    ra1,dec1   RA,DEC at the jday1 epoch (radians).
     c    jday2      Julian day of the new epoch.
     c  Output:
     c    ra2,dec2   Precessed coordinates (radians) */
  double r0,d0,rm,dm,T,M,N,ra2,dec2;
  T  = (jday1 - 2451545.0)/36525;
  M  = DPI/180 * (1.2812323 + (0.0003879 + 0.0000101*T)*T)*T;
  N  = DPI/180 * (0.5567530 - (0.0001185 + 0.0000116*T)*T)*T;
  rm = ra1 - 0.5*(M + N*sin(ra1)*tan(dec1));
  dm = dec1 - 0.5*N*cos(rm);
  
  /*   J2000 coordinates */
  r0 = ra1 - M - N*sin(rm)*tan(dm);
  d0 = dec1 - N*cos(rm);
  /* Coordinates of the other epoch. */
  T = (jday2 - 2451545.0)/36525.0;
  M = DPI/180.0 * (1.2812323 + (0.0003879 + 0.0000101*T)*T)*T;
  N = DPI/180.0 * (0.5567530 - (0.0001185 + 0.0000116*T)*T)*T;
  rm = r0 + 0.5*(M + N*sin(r0)*tan(d0));
  dm = d0 - 0.5*N*cos(rm);
  ra2 = r0 + M + N*sin(rm)*tan(dm);
  dec2 = d0 + N*cos(rm);
  *ra2pt=ra2;
  *dec2pt=dec2;
}

void nutate(double jday,double rmean,double dmean,double *rtrueptr,double *dtrueptr)
{
  /* jhz 2004-7-23: based on miriad code in f, translate into c
     c  Convert between mean and true equatorial coordinates, by
     c  accounting for nutation.
     c
     c  Input:
     c    jday       Julian day.
     c    rmean,dmean Mean (RA,DEC) at jday.
     c  Output:
     c    rtrue,dtrue True (RA,DEC) at jday.
  */
  double deps,dpsi,eps, rtrue, dtrue;
  double coseps,sineps,sinra,cosra,tandec;
  
  /*  Nutation parameters. */
  nuts(jday,&dpsi,&deps);
  /*  True obliquity. */
  eps = mobliq(jday) + deps;
  /*  Various parameters. */
  sineps = sin(eps);
  coseps = cos(eps);
  sinra  = sin(rmean);
  cosra  = cos(rmean);
  tandec = tan(dmean);   
  
  rtrue = rmean + (coseps + sineps*sinra*tandec)*dpsi
    - cosra*tandec*deps;
  dtrue = dmean + sineps*cosra*dpsi + sinra*deps;
  *rtrueptr = rtrue;
  *dtrueptr = dtrue;
  /*   printf("nutate: r1 d1 %f %f\n", rtrue, dtrue);
   */
}

void nuts(double ljday, double *dpsiptr,double *depsptr)
{
  /* jhz 2004-7-23: based on miriad code in f and translate into c.
     c
     c  Return nutation parameters. The claimed accuracy is 1 arcsec.
     c
     c  Input:
     c    jday       Julian date.
     c  Output:
     c    dpsi,deps  Difference between mean and true ecliptic latitude and
     c               longitude due to nutation, in radians.
     c
     c  Reference:
     c    Explanatory Supplmenet, page 120.
     c--
  */
  double d,t1,t2, dpsi,  deps;
  d = jday - 2451545.0;
  t1 = DPI/180*(125.0 - 0.05295 * d);
  t2 = DPI/180*(200.9 + 1.97129 * d);
  dpsi = DPI/180 * (-0.0048*sin(t1) - 0.0004*sin(t2));
  deps = DPI/180 * ( 0.0026*cos(t1) + 0.0002*cos(t2));
  *dpsiptr=dpsi;
  *depsptr=deps;
}

double mobliq(double jday)
{
  /* jhz 2004-7-23: based on miriad code in f and translate into c.
     c
     c  Return the mean obliquity of the ecliptic.
     c
     c  Input:
     c    jday       Julian day.
     c  Output:
     c    mobliq     Mean obliquity of the ecliptic, in radians.
     c
     c  Reference:
     c    Explanatory Supplement ... page 114.
     c-- */
  double T;
  double vmobliq;
  /* Centuries from J2000*/
  T = (jday - 2451545.0) / 36525.0;
  /* Mean obliquity.*/    
  vmobliq = 84381.448 - (46.8150+(0.00059-0.001813*T)*T)*T;
  vmobliq = DPI/(180.*3600.) * vmobliq;    
  return vmobliq;
}

void aberrate(double jday,double ra,double dec,double *rappptr,double *dappptr)
{
  // jhz 2004-7-23: based on miriad code in f and translate into c.
  //  Account for the effect of annual aberration, to convert
  //  from a true (RA,DEC) to a geocentric apparent (RA,DEC).
  //
  //  Input:
  //    jday       Julian date.
  //    ra,dec     True (RA,DEC).
  //  Output:
  //    rapp,dapp  Geocentric apparent (RA,DEC).
  //
  double  pos[3],vel[3],sinra,sindec,cosra,cosdec, rapp, dapp;
  void vearth();

  vearth(jday,pos,vel);
  sinra = sin(ra);
  cosra = cos(ra);
  sindec = sin(dec);
  cosdec = cos(dec);
  rapp = ra +  (-vel[0]*sinra + vel[1]*cosra)/
    (0.001*DCMKS*cosdec);
  dapp = dec + (-vel[0]*cosra*sindec - vel[1]*sinra*sindec
		+ vel[2]*cosdec)/(0.001*DCMKS);
  *rappptr= rapp;
  *dappptr= dapp;
  
}

void vearth (double jday, double pos[3], double vel[3])
{
  /* jhz 2004-7-23: based on miriad code in f and translate into c.
   *
   *  Approximate heliocentric position and velocity of the Earth
   *  The date and time should really be in the TDB in the Gregorian
   *  calendar, and is interpreted in a manner which is valid between
   *  1900 March 1 and 2100 February 28.
   *
   *  Input:
   *    jday       Time of interest (as Julian day).
   *  Output:
   *    pos        Position, in km.
   *    vel        Velocity, in km/sec.
   *
   *  The Earth heliocentric position/velocity is for mean equator and equinox
   *  of date.
   *
   *  Max/RMS errors 1950-2050:
   *     13/5 E-5 AU = 19200/7600 km in position
   *     47/26 E-10 AU/s = 0.0070/0.0039 km/s in speed
   */
  float twopi, speed, remb, semb;
  double aukm, j1901;
  float YF,T,ELM,GAMMA,EM,ELT,EPS0,DAY;
  float E,ESQ,V,R,ELMM,COSELT,SINEPS,COSEPS,W1,W2,SELMM,CELMM;
  int IY, QUAD;
  twopi = (float) 2*DPI;
  /* Mean orbital speed of Earth, AU/s */
  speed = 1.9913E-7;
  /* Mean Earth:EMB distance and speed, AU and AU/s */
  remb = 3.12E-5;
  semb = 8.31E-11;
  /* AU to km */
  aukm=149.597870e6;
  /* Julian date for 1 January, 1901. */
  j1901=2415385.5;
  /* Whole years & fraction of year, and years since 1900.*/
  QUAD = (int)((jday - j1901) / 1461.);
  IY   = (int)((jday - j1901 - 1461*QUAD) / 365.0);
  DAY  = (float)(jday - j1901 - 1461*QUAD - 365*IY + 1);
  IY   = 4*QUAD + IY + 1;
  YF   = (4*DAY - 4*(1/(fmod(IY,4)+1)) - fmod(IY,4) - 2) / 1461.0;
  T    = IY + YF;
  /* Geometric mean longitude of Sun
   *  (cf 4.881627938+6.283319509911*T MOD 2PI)*/
  ELM  = fmod(4.881628+twopi*YF+0.00013420*T, twopi);
  /*  Mean longitude of perihelion */
  GAMMA=4.908230+3.0005e-4*T;
  /*  Mean anomaly */
  EM=ELM-GAMMA;
  /*  Mean obliquity */
  EPS0=0.40931975-2.27e-6*T;
  /*  Eccentricity  */
  E=0.016751-4.2e-7*T;
  ESQ=E*E;
  /*  True anomaly */
  V=EM+2.0*E*sin(EM)+1.25*ESQ*sin(2.0*EM);
  /*  True ecliptic longitude*/
  ELT=V+GAMMA;
  /*  True distance */
  R=(1.0-ESQ)/(1.0+E*cos(V));
  /*  Moon's mean longitude */
  ELMM=fmod(4.72+83.9971*T, twopi);
  /*  Useful functions */
  COSELT=cos(ELT);
  SINEPS=sin(EPS0);
  COSEPS=cos(EPS0);
  W1=-R*sin(ELT);
  W2=-speed*(COSELT+E*cos(GAMMA));
  SELMM=sin(ELMM);
  CELMM=cos(ELMM);
  /*  Earth position and velocity*/
  pos[0] = aukm * (-R*COSELT-remb*CELMM);
  pos[1] = aukm * (W1-remb*SELMM)*COSEPS;
  pos[2] = aukm * W1*SINEPS;
  vel[0] = aukm * (speed*(sin(ELT)+E*sin(GAMMA))+semb*SELMM);
  vel[1] = aukm * (W2-semb*CELMM)*COSEPS;
  vel[2] = aukm * W2*SINEPS;
}

void elaz(int tno) {
  /* calculate and store mean az and el into uv data */
  int i;
  double mel, maz;
  mel=0;
  maz=0;

  if(smabuffer.nants!=0) {
    for (i=0; i<smabuffer.nants; i++) {
      mel=mel+smabuffer.el[i];
      maz=maz+smabuffer.az[i];
    }
    mel=mel/smabuffer.nants;
    maz=maz/smabuffer.nants;
    /* both az and el in degree */
    //         uvputvrd_c(tno,"antaz",&maz,1);
    //         uvputvrd_c(tno,"antel",&mel,1);
    // bee.for require store antaz and antel for each antenna
    uvputvrd_c(tno,"antaz",&smabuffer.az,smabuffer.nants);
    uvputvrd_c(tno,"antel",&smabuffer.el,smabuffer.nants);
  }
}

void tsysStore(int tno) {
  /* store Tsys to uvdata */
  int cnt, j;
  float tsysbuf[SMANT*SMIF];

  cnt =0;
  /* only one  if */
  for (j=0; j< smabuffer.nants; j++) {
    tsysbuf[cnt]= smabuffer.tsys[j];
    /*      if(smabuffer.tsys[j]<0.) tsysbuf[cnt]=0.;
     */
    cnt++;
  }
  uvputvrr_c(tno,"systemp",&tsysbuf,cnt);
}

double velrad( short dolsr,
               double time,
               double raapp,
               double decapp,
               double raepo,
               double decepo,
               double lst,
               double lat) {
  //  Based on Miriad subroutine velrad in atlod.for
  //  Compute the radial velocity of the observatory, 
  //  in the direction of a source, with respect to either 
  //  LSR or the barycentre. The subroutine is based on
  //  the fortran routine in miriad.
  //
  //  Input:
  //    dolsr      If >0, compute LSR velocity. Otherwise barycentric.
  //    time       Time of interest (Julian date).
  //    raapp,decapp Apparent RA and DEC (radians).
  //    raepo,decepo RA and DEC at the J2000 epoch (radians).
  //    lat        Observatory geodetic latitude (radians).
  //    lst        Local sideral time (radians).
  //  Output:
  //    vel        Radial velocity.
  struct lmn *inlmn;
  double lmn2000[3], lmnapp[3];
  struct vel *velosite;
  double  posearth[3], velsite[3], velearth[3], velsun[3];
  int i;
  double  vel;
  
  // computer barycentric velocity
  //      printf("velrad:\n");
  //      printf("dolsri %d \n", dolsr);
  //      printf("time %f \n", time);
  //      printf("raapp %f \n", raapp);
  //      printf("decapp %f \n", decapp);
  //      printf("raapp %f \n", raepo);
  //      printf("decapp %f \n", decepo);
  //      printf("lst %f \n", lst);
  //      printf("lat %f \n", lat);
  
  inlmn = sph2lmn(raapp,decapp);
  lmnapp[0] = inlmn->lmn1;
  lmnapp[1] = inlmn->lmn2;
  lmnapp[2] = inlmn->lmn3;
  
  velosite = vsite(lat,lst);
  
  velsite[0] = velosite->vx;
  velsite[1] = velosite->vy;
  velsite[2] = velosite->vz;
  vearth(time, posearth, velearth);
  vel =0.;
  for (i=0; i<3; i++) {
    vel = vel - (velsite[i] + velearth[i])*lmnapp[i];
  }
  
  // To compute LSR velocity, we need the source position in J2000 coordinates.
  // Vsun returns the Suns LSR velocity in the J2000 frame. Add this
  // contribution to the velocity we already have.
  
  if(dolsr==1) {
    inlmn= sph2lmn(raepo,decepo);
    lmn2000[0]=inlmn->lmn1;
    lmn2000[1]=inlmn->lmn2;
    lmn2000[2]=inlmn->lmn3;
    vsun(velsun);              
    for (i=0; i<3; i++ ) {
      vel = vel + lmn2000[i]*velsun[i];
    }
  }
  return vel;
}

struct lmn *sph2lmn(double ra, double dec) {
  //Convert from spherical coordinates to direction cosines.
  //Convert spherical coordinates (e.g. ra,dec or long,lat) into
  //direction cosines.
  //Input:
  //    ra,dec     Angles in radians.
  //Output:
  //    lmn        Direction cosines.
  struct lmn *outlmn;
  outlmn = (struct lmn *)malloc(sizeof(struct lmn ));
  outlmn->lmn1 = cos(ra)*cos(dec);
  outlmn->lmn2 = sin(ra)*cos(dec);
  outlmn->lmn3 = sin(dec);
  return outlmn;
}

struct vel *vsite(double phi, double st) {
  // based on miriad fortran subroutine vsite in velocity.for
  //*  Velocity due to Earth rotation 
  //*
  //*  Input:
  //*     PHI       latitude of observing station (geodetic)
  //*     ST        local apparent sidereal time
  //*  Output:
  //*     VEL       velocity in km/s.
  //*
  //*  PHI and ST are all in radians.
  //*  Accuracy:
  //*     The simple algorithm used assumes a spherical Earth and
  //*     an observing station at sea level.  For actual observing
  //*     sites, the error is unlikely to be greater than 0.0005 km/s.
  //*  Sidereal speed of Earth equator, adjusted to compensate for
  //*  the simple algorithm used.  (The true value is 0.4651.)
  //*  in units of km/s
  float espeed=0.4655;
  struct vel *sitevel;
  sitevel = (struct vel *)malloc(sizeof(struct vel ));
  sitevel->vx = -espeed*cos(phi)*sin(st);
  sitevel->vy =  espeed*cos(phi)*cos(st);
  sitevel->vz =  0.0;
  return sitevel;
}


void vsun(double *VEL) {
  // based on Miriad subroutine vsub in velocity.for
  //  Velocity of the Sun with respect to the Local Standard of Rest
  //
  //  Output:
  //     VEL       Velocity of the Sun.
  //------------------------------------------------------------------------
  //  Speed = 20 km/s
  //
  //  Apex = RA 270 deg, Dec +30deg, 1900.0
  //  = 18 07 50.3, +30 00 52, J2000.0
  //
  //  This is expressed in the form of a J2000.0 x,y,z vector:
  //
  //      VA(1) = X = -SPEED*COS(RA)*COS(DEC)
  //     VA(2) = Y = -SPEED*SIN(RA)*COS(DEC)
  //      VA(3) = Z = -SPEED*SIN(DEC)
  //      DATA VA / -0.29000, +17.31726, -10.00141 /
  float VA[3];
  VA[0] = -0.29000;
  VA[1] = +17.31726;
  VA[2] = -10.00141;
  VEL[0] = VA[0];
  VEL[1] = VA[1];
  VEL[2] = VA[2];
}

short ipolmap(short input_ipol) {
// mapping ipol to miriad polarization states id
short iPolmiriad;
int   ipol;
      ipol=input_ipol;
      iPolmiriad = 1;
if(smabuffer.oldpol==1) {
            switch(ipol) {
            case 0: iPolmiriad = 1; break;
            case 1: iPolmiriad =-5; break;
            case 2: iPolmiriad =-7; break;
            case 3: iPolmiriad =-8; break;
            case 4: iPolmiriad =-6; break;
              // convert MIR polarization label used befor sep1,2004 to Miriad
              // used   MIR  actual          Miriad
              //non      0   I                 1
              //RR       1   HH               -5
              //RL       2   HV               -7
              //LR       3   VH               -8
              //LL       4   VV               -6
                      }
                        } else {
if(smabuffer.circular==1)
              switch(ipol) {
              case 1: iPolmiriad =-1; break;
              case 2: iPolmiriad =-3; break;
              case 3: iPolmiriad =-4; break;
              case 4: iPolmiriad =-2; break;
                //MIR     MIRIAD     STATE
                //1       -1         RR
                //2       -3         RL
                //3       -4         LR
                //4       -2         LL
                           }
if(smabuffer.linear==1)
              switch(ipol) {
              case 1: iPolmiriad =-6; break;
              case 2: iPolmiriad =-7; break;
              case 3: iPolmiriad =-8; break;
              case 4: iPolmiriad =-5; break;
// from ram mar7,2005
                // In the case where options=linear then
                // MIR MIRIAD  POL
                // 1   -6      VV (YY)
                // 2   -7      HV (XY)
                // 3   -8      VH (YX)
                // 4   -5      HH (XX)
                          }
                 }
        if(smabuffer.nopol==1) iPolmiriad=-5;

     return iPolmiriad;
}
