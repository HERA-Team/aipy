// jhz 2006-1-31: change the symbolic SMRX 
//                from 5 to 2 (number of sma rx per operation run)
//                so that to reduce the memory requirements
//                for the data structure wt.
// jhz 2006-2-1: the symbolic constant SMRX is also used as the maximum size
//               of the rx-id array 0(230), 1(340), 2(690), 3(??), 4(??)
//               used in the MIR data. It must be 5 (it can be 3 for now but not 2).
//               SMIF can be 24 instead of 48 which was redundant since in the wt array
//               a sideband variable is used. 
// jhz 2006-2-3: add dsb to smlodd.
// jhz 2006-5-17: add highrspectra to smlodd to handle 
//                high spectral resolution mode, allowing 
//                empty chunks.
// jhz 2006-5-18: change MAXCHAN from 7681 to 8217
#include "miriad.h"


struct inh_def {

	int 	conid     ; /* config id #               */
	short	icocd     ; /* config int code           */
	int 	traid     ; /* track id #                */
	int 	inhid     ; /* integration id #          */
	int 	ints      ; /* integration #             */
	short	itq       ; /* tuning int code           */
	float	az        ; /* azimuth                   */
	float	el        ; /* elevation                 */
	float	ha        ; /* hour angle                */
	short	iut       ; /* ut int code               */
	short	iref_time ; /* ref_time int code         */
	double	dhrs      ; /* hrs from ref_time         */
	float	vc        ; /* vcorr for vctype          */
	short	ivctype   ; /* vctype int code           */
	double	sx        ; /* x vec. for bsl.           */
	double	sy        ; /* y vec. for bsl.           */
	double	sz        ; /* z vec. for bsl.           */
	float	rinteg    ; /* actual int time           */
	int 	proid     ; /* project id #              */
	int 	souid     ; /* source id #               */
	short	isource   ; /* source int code           */
	short	ipos      ; /* position int code         */
	float	offx      ; /* offset in x               */
	float	offy      ; /* offset in y               */
	short	iofftype  ; /* offset int code           */
	short	ira       ; /* ra int code               */
	short	idec      ; /* dec int code              */
	double	rar       ; /* ra (radians)              */
	double	decr      ; /* declination (radians)     */
	float	epoch     ; /* epoch for coord.          */
	float	sflux     ; /* flux                      */
	float	size      ; /* source size               */

};

/* inh is 132 bytes if short=2 int=4 float=4 double=8 */

struct blh_def {

	int 	blhid     ; /*	proj. baseline id #       */
	int 	inhid     ; /*	integration id #          */
	short	isb       ; /*	sideband int code         */
	short	ipol      ; /*	polarization int code     */
	float	pa        ; /*	polarization angle        */
	short	iaq       ; /*	amp qual int code         */
	short	ibq       ; /*	baseline qual int code    */
	short	icq       ; /*	coherence qual int code   */
	short	ioq       ; /*	offset qual int code      */
	short	irec      ; /*	receiver int code         */
	short	iifc      ; /*	if channel int code       */
	float	u         ; /*	u coord for bsl (klambda) */
	float	v         ; /*	v coord. for bsl          */
	float	w         ; /*	w coord. for bsl          */
	float	prbl      ; /*	projected baseline        */
	float	angres    ; /*	baseline resolution       */
	float	vis       ; /*	fractional visibility     */
	float	coh       ; /*	coherence estimate        */
	float	sigcoh    ; /*	sigma on coh.             */
	float	csnr      ; /*	cont. SNR                 */
	float	vflux     ; /*	visibility flux           */
	float	cnoise    ; /*	continuum noise           */
	double	avedhrs   ; /*	hrs offset from ref-time  */
	float	ampave    ; /*	ave continuum amp         */
	float	phaave    ; /*	ave continuum phase       */
	float	tpvar     ; /*	total power stability     */
	int 	blsid     ; /*	physical baseline id #    */
	short	itel1     ; /*	tel 1 int code            */
	short	itel2     ; /*	tel 2 int code            */
	short	iblcd     ; /*	baseline int code         */
	float	ble       ; /*	bsl east vector           */
	float	bln       ; /*	bsl north vector          */
	float	blu       ; /*	bsl up vector             */
	int 	soid      ; /*	bsl soln id #             */

};

/* the size of blh_def is 118 bytes */

struct sph_def {

	int 	sphid     ; /*  spectrum id #             */
	int 	blhid     ; /*  proj. baseline id #       */
	int 	inhid     ; /*  integration id #          */
	short	igq       ; /*  gain qual int code        */
	short	ipq       ; /*  passband qual int code    */
	short	iband     ; /*  spectral band int code    */
	short	ipstate   ; /*  pol state int code        */
	float	tau0      ; /*  tau at zenith             */
	double	vel       ; /*  velocity (vctype) (km/s)  */
	float	vres      ; /*  velocity res.             */
	short	ivtype    ; /*  velocity type int code    */
	double	fsky      ; /*  center sky freq.          */
	float	fres      ; /*  frequency res. (MHz)      */
	float	tssb      ; /*  tsys (ssb)                */
	float	integ     ; /*  integration time          */
	float	wt        ; /*  weight (sec/tssb**2)      */
	short	itaper    ; /*  spectrum taper int code   */
	float	snoise    ; /*  th. noise in spectrum     */
	short	nch       ; /*  # channels in spectrum    */
	short	nrec      ; /*  # of records w/i inh#     */
	int 	dataoff   ; /*  byte offset for data      */
	int 	linid     ; /*  spectral line id #        */
	short	itrans    ; /*  transition int code       */
	double	rfreq     ; /*  rest frequency (GHz)      */
	short	pasid     ; /*  passband fn id #          */
	short	gaiidamp  ; /*  amp gain fn id #          */
	short	gaiidpha  ; /*  pha gain fn id #          */
	short	flcid     ; /*  flux cal id #             */
	short	atmid     ; /*  atmospheric coh id #      */

};

/* the size of sph_def is 100 bytes */

struct codeh_def {

	char v_name[12]	; /* label 			*/
	short  icode	; /* index for a code word	*/
	char code[26]	; /* the code word		*/
	short ncode	; /* # chars in code word	*/
};

struct xyz {
        double x;
        double y;
        double z;
};

/* the size of codeh is 42 bytes char=1 short=2 */


/*
;   In the database archive, each record of data consists of 4 header
;   numbers followed by a variable number of complex channels.
;   The beginning of the record data stream is given by dataOFF from SPH.
;   Record format:
;      ITEM     SIZE OFF
;     integ      i*2   0       Integration time in .1secs
;     toff       i*2   2       Time offset from INH ut in secs
;     dnoise     i*4   4       Noise estimate in milliJanskies (will be read 
				by the data processing program as 2 shorts and
				then reconstituted as a long)
;     scale      i*2   8       Scale exponent for data
;     ch[1].real i*2  10       Signed, scaled, short (2-byte integer) data
;     ch[1].imag i*2
;     ch[2].real i*2
;     ch[2].imag i*2
;     ...
;
;   The total size of each record is 10+nch*4 bytes.
;   Each data point is stored as a 2 byte (short) integer. The restored data
;   (in Janskies) are computed as follows:
;                 restored_data = 2**scale * short_data
;
;   The noise estimate (dnoise) is stored as a 4byte integer in millijanskys.
*/

struct sch_def {
	int inhid	;  /* integration id # */
	char form[4]	;  /* indicates the format of the data, always I2-C */
	int nbyt	;  /* the number of bytes in one integration of data */
	int nbyt_pack	;  /* ? */
	short *packdata	;  /* integer array containing the data in the format above	*/
};

/* 
The dimension of the packdata integer array:  
number of points in one record		=	5 + 2* number of channels	
number of records in each spectrum	=	nrec (always 1 for our data for now)
number of spectra in each baseline	=	numberBands (in this code)
number of baselines in each integration	=	numberBaselines (in this code)

the number of sph header structures written is numberBands*numberBaselines, 
and this is called the number of spectra in the idl programs.
*/

struct rh_def {
	int   noiseEst   ;  /* noise estimate for record in (mJy)                  */
	short integTime  ;  /* The integration time of the record in 0.1 secs      */
	short timeAfter  ;  /* Time offset from INH UT time in secs                */                   
};


struct ant_def{
	int antennaNumber;
	int padNumber;
	int antennaStatus;  /* Antenna is ON or OFF LINE */
	int trackStatus;    /* Track is running or not */
	int commStatus;     /* Data for this integration is valid or not */

	int 	inhid     ; /* integration id #          */
	int 	ints      ; /* integration #             */
	double	dhrs      ; /* hrs from ref_time         */
	double	ha        ; /* hour angle                */
	double	lst       ; /* lst                       */

        double pmdaz ;                       /* pointing model correction */
        double pmdel ;
        double tiltx ;
        double tilty ;
 
        double actual_az ;
        double actual_el ;
        double azoff ;
        double eloff ;
        double az_tracking_error ;
        double el_tracking_error ;
        double refraction ;
        double chopper_x ;
        double chopper_y ;
        double chopper_z ;
        double chopper_angle ;
 
        double tsys ;
        double ambient_load_temperature ;

};

struct sph_config {
        int     sphid         ; /*  spectrum id #             */
        int     inhid         ; /*  integration id #          */
        short   iband[25]     ; /*  spectral band int code    */
        double  vel[25]       ; /*  velocity (vctype) (km/s)  */
        float   vres[25]      ; /*  velocity res.             */
        short   ivtype        ; /*  velocity type int code    */        
        double  fsky[25]      ; /*  center sky freq.          */
        float   fres[25]      ; /*  frequency res. (MHz)      */
        short   nch[25][2]       ; /*  # channels in spectrum    */
        int     dataoff       ; /*  byte offset for data      */
        double  rfreq[25]     ; /*  rest frequency (GHz)      */
        double  basefreq      ; /*  determine the basefreq    */
        float   veldop        ; /*  velocity of observatory in the 
                                    direction of the source tracked 
                                    with doppler correction   */
        float   smaveldop     ; /* the residual veldop after taking out
                                   the amount corrected to the 
                                   chunk frequency by the SMA on-line
                                   system. */
        short   isb           ; /*  sideband int code         */
        short   irec          ; /*  receiver int code         */
        int     souid         ; /* source id #               */
};     

struct blh_config {
        int     blhid     ; /*  proj. baseline id #       */
        int     inhid     ; /*  integration id #          */
        short   isb       ; /*  sideband int code         */
        short   ipol      ; /*  polarization int code     */
        short   irec      ; /*  receiver int code         */
};

struct bltsys  {
        int     blhid     ; /*  proj. baseline id #       */
        int     inhid     ; /*  integration id #          */
        int     blsid     ; /*  physical baseline id #    */
        short   isb       ; /*  sideband int code         */
        short   irec      ; /*  receiver int code         */
        short   ipol      ; /*  polarization int code     */
        short   itel1     ; /*  tel 1 int code            */
        short   itel2     ; /*  tel 2 int code            */
        float   tssb[25]  ; /*  tsys (ssb) 25 spectra     */
};

struct anttsys {
        int    inhid      ; /*  integration id #          */
        short   refant    ; /*  reference ante            */
        short   refpair1  ; /*  reference pair 1          */
        short   refpair2  ; /*  reference pair 2          */
        short   isb       ; /*  sideband int code         */
        short   irec      ; /*  receiver int code         */
        float   tssb[30]  ; /*  upto 30 ant               */
};

       
/* the following definition and struct for smalod */

#define DPI 3.14159265358979323846   /* pi */
#define DCMKS          299792458.0   /* Speed of light (meters/second). */

/* although we could use maxdimc.h, it turns out some variables 
 * are tied to the SMA. it's safer to define them here .
 * Of course your real maxdimc.h better have at least as large as
 * these or else i will not be able to process the subsequent data
 */
#define MAXINT 5000  /* maximum number of integration */
#define MAXANT 10
#define MAXCHAN 8217 /* 7681 8193+24=8217*/

#define MAXBAS 90    /* maxant*(maxant-1)/2*2sb */

#define MAXSOURCE 50
#define SMIF 24
#define SMRX  5     /* number of rx per track operattion */ 
#define SMANT 10
#define SMPOL 5 
#define SMBAS 90   /* smant*(smant-1) */
#define SMBIN 1
#define SMSB  2    /* number of size bands */
#define SMCONT 33
#define CONTCH 16 /* number of continuum chan per chaunk */
#define SMADATA 8294400 /* 24*maxchan*smbase */

/* WORDS_BIGENDIAN comes from miriad's sysdep.h */
/* SWAP_ENDIAN is what SMA code originally used  */
#if defined(WORDS_BIGENDIAN)
# define SWAP_ENDIAN     0 /* for big endian computers  (e/g. sparc, ppc) */
#else
# define SWAP_ENDIAN     1 /* for little endian computers (e.g. intel) */
#endif

struct wtt {
        float wt[SMIF][SMPOL][MAXBAS][SMSB][SMRX];
            };
                                                                                              

struct uvw {
        double u;
        double v;
        double w;
        int blhid;
        int blsid;
        int blcode;  /* itel1*256+itel2 */
        int isb;     /* sideband id */
        int ipol;    /* pol code */
        int irec;    /* rx id    */
};
typedef struct uvw uvw;

struct uvwPack {
        uvw uvwID[MAXBAS];
        int inhid;
        int n_bls; /* number of baselines in each integration */
};
typedef struct uvwPack uvwPack;

struct SMAdataHdr {
        double date;
        double time;
        int bsln;
        int arrayId;
        int srcId;
        int freqId;
        float inttim;
        float weight[48];
};
typedef struct SMAdataHdr SMAdataHdr;


struct timeStamp {
        double UTCtime;
        double intTime;
};
typedef struct timeStamp timeStamp;

struct visdata {
         float real;
         float imag;
};
typedef struct visdata visdata;

struct visdataChunk {
        visdata chunk[1];   /*visdata chunk[MAXSMCHAN];*/
};
typedef struct visdataChunk visdataChunk;
 
struct visdataSp {
        float       wt[SMIF];
        visdata continuum;
        visdataChunk sp[SMIF];
};
typedef struct visdataSp visdataSp;

struct packetIdentity {
        int ints;
        int inhid;
        int polid;
        int sbid;
        int sourID;
        int freqID;
};
typedef struct packetIdentity packetIdentity;

struct visdataBlock {
        packetIdentity blockID;
        timeStamp time;        
        int uvblnID; 
        visdataSp bsline; 
};
typedef struct visdataBlock visdataBlock;

struct source {
       char name[24];
       double ra;
       double dec;
       char equinox[8];
       double ra_app;
       double dec_app;
       char calcode[4];
       int qual;
       int freqid;
       double pmra;
       double pmdec;
       float parallax;
       double sysvel[SMIF+1];
       char veltyp[8];
       char veldef[8];
       double restfreq[SMIF+1];
       int sour_id;
       int inhid_1st;
};
typedef struct source source;

struct tsys {
        double time;
        float time_interval;
        int sour_id;
        int antenna_no;
        int array_no;
        int freq_id;
        float tsys1[SMIF+1];
        float tant1[SMIF+1];
        float tsys2[SMIF+1];
        float tant2[SMIF+1];
};
typedef struct tsys tsys;

struct correlator {
        char corr_name[8];
        int no_stkd;
        int stk_1;
        int no_sideband;
        int n_chunk;
        int no_rxif;
};
typedef struct correlator correlator;

struct frequency {
        double chunkBW[SMIF+1];
        double chunkfreq[SMIF+1];
        double chanWidth[SMIF+1];
        int ref_chan[SMIF+1]; /* reference channel for chunk freq */
        int n_chunk_ch[SMIF+1];
        int freqid;
        int sideband_id;
        int polarization[SMIF+1];
        int rxif;    /* rxif id */
};
typedef struct frequency frequency;

struct  blvector {
        float ee;
        float nn;
        float uu;
        short itel1;
        short itel2;
        int blid;
        };
typedef struct blvector blvector;

struct  station {
        char name[8];  /* station name */
        double x;
        double y;
        double z;
        double x_phs;
        double y_phs;
        double z_phs;
        int axistype;
        float axisoff_x;
        float axisoff_y;
        float axisoff_z;
};
typedef struct station station;

struct pols {
        int npol;
        int polstart;
        int polend;
            };

struct lmn {
        double lmn1;
        double lmn2;
        double lmn3;
           };

struct vel {
        double vx;
        double vy;
        double vz;
           };

struct smlodd {
        double sfreq[SMIF+1];
        double sdf[SMIF+1];
        double restfreq[SMIF+1];
        double basefreq;
        double time; /* in Julian day */
        double juldate; /* julilan date */
        double ut;   /* in radians    */
        double lst;  /* in units of radian */
        double obsra;
        double obsdec;
        double lat;
        double longi;
        double ra;
        double dec;
        double el[SMANT+1];
        double az[SMANT+1];
        double antpos[SMANT*3+1];
        visdata data[SMADATA+1];
        float veldop;
        float smaveldop;
        float vsource;
        float xtsys[SMIF+1][SMANT+1];
        float ytsys[SMIF+1][SMANT+1];
        float tsys[SMANT+1];        
        float chi;
        float xyphase[SMIF+1][ SMANT+1]; 
        float xyamp[SMIF+1][SMANT+1]; 
        float xsampler[3][SMIF+1][SMANT+1];
        float ysampler[3][SMIF+1][SMANT+1];
        float u[SMBAS+1];
        float v[SMBAS+1];
        float w[SMBAS+1];
        float blcode[SMBAS+1];    /* ant1*256 + ant2 */
        float inttime[SMBAS+1];
        float inttim;
        float wts[2*SMCONT-2];
        float mdata[5];
        float axisrms[SMANT+1];
        float axismax[SMANT+1];
/*        int pnt[SMIF+1][SMPOL+1][SMBAS+1][SMBIN+1];*/
        int pnt[SMIF+1][SMPOL+1][SMBAS+1][SMSB][SMRX];
        int nbin[SMIF+1];
        int nused;
        int tno;
        int rsnchan;
        int nants;
        int refant;
        int readant;
        int nifs;
        int nfreq[SMIF+1];
        int nstoke[SMIF+1];
        int polcode[SMIF+1][SMPOL+1][SMBAS+1];
        int edge[SMIF+1];
        int bchan[SMIF+1];
        int tcorr;
/*        int flag[SMIF+1][SMPOL+1][SMBAS+1][SMBIN+1];*/
        int flag[SMIF+1][SMPOL+1][SMBAS+1][SMSB][SMRX];   
        int dosw[SMBAS+1];
        int dosam;
        int dohann;
        int birdie;
        int dowt;
        int dopmps;
        int doxyp;
        int opcorr;
        int doif;
        int dobary; 
        int dolsr;
        int noskip;
        int newfreq;
        int dorfreq;
        int hires;
        int nopol;
        int circular;
        int linear;
        int oldpol;
        int doChunkOrder;
        int doConjugate;
        int mflag;
        int newsc;
        int newpnt;
        int sb;   /* side band sb=0 for lsb, 1 for usb, sb for both */
        int rxif; /* corresponding to irec in blh */
        int rx1;  /* the first rx in the data structure in the dual rx case */
        int rx2;  /* the second rx in the data structure in the dual rx case */
        int doeng; /* 1 to read the engine file for Tsys and LST */
        int scanskip;
        int scanproc;
        int currentscan;
        int spskip[2];
        int dsb; 
        int mcconfig;
        int highrspectra; /* highrspectra=1; full spectral chunks (24)
                                           to be processed, allowing
                                           empty chunks.
                             highrspectra=-1; partial spectral chunks
                                           to be processed, the
                                           number of spectral chunks
                                           to be determined by the 
                                           program. */
    };
typedef struct smlodd smlodd;

struct smEng {
        int inhid;
        int ints;
        int sour_id;
        int freq_id;
        int antpad_no[MAXANT];
        int antenna_no[MAXANT];
        double lst;
        double dhrs;
        double ha;
        double el[MAXANT];
        double az[MAXANT];
 /*       float refraction[MAXANT];
        float pmdaz[MAXANT];
        float pmdel[MAXANT];
        float azoff[MAXANT];
        float eloff[MAXANT];
        float chopper_x[MAXANT];
        float chopper_y[MAXANT];
        float chopper_z[MAXANT];
        float chopper_angle[MAXANT];
        float tilt_x[MAXANT];
        float tilt_y[MAXANT];
        float tilt_x_dc[MAXANT];
        float tilt_y_dc[MAXANT];*/
        double tsys[MAXANT];
        double tamb[MAXANT];
};
typedef struct smEng smEng;

/* function declarations */

void reverse2(char *bytes);
void reverse4(char *bytes);
void reverse8(char *bytes);

struct inh_def *swap_inh(struct inh_def *inh_pnr);
struct blh_def *swap_blh(struct blh_def *blh_pnr);
struct sph_def *swap_sph(struct sph_def *sph_pnr);
struct codeh_def *swap_cdh(struct codeh_def *cdh_pnr);
struct ant_def *swap_enh(struct ant_def *enh_pnr);
struct sch_def *swap_sch(struct sch_def *sch_pnr);
short int *swap_sch_data(short int *sch_data_pnr, int datalength);
int *swap_int_data(int *in_data_pnr, int datalength);


