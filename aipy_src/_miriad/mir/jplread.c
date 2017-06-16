// jhz 15dec04 modify jpl catalog handling routines in c;
//             interface with the miriad fortran routine smauvspec.for
//             retrieve jpl line catalog and plot them on xterm device
//             using pgplot routine and compare with the uv spectral
//             line data.

#define CATDATA struct catdata
CATDATA {
  double freq, derr, str, elow;
  int itd, igup, tag, ifmt;
  short iqn[12];
};
char *catfil(int num);

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define OK 1

char path[80];
int  length;

void molselect_c(jplpath, pathlen, mtag, nmline, mname)
char *jplpath;
int pathlen;
int *mtag;
char *nmline[];
char *mname;
{ int iver, moltag, i, nline, j, k, kp;
  char *molname, *catdir(), select[81], moln[16], mbuff[16*500];
  char *molnp;
  float qrot[7];
  extern char path[80];
  extern int length;
       sprintf(path, "%s", jplpath);
       length=pathlen;
//       printf("pathlen=%d jplpath=%s\n", length, path);
       moltag=0; 
       nline=0;  
     while (nxtdir(&moltag)) {
   //    molname = catdir(moltag, &i, qrot, &iver);
     molname = catdir(moltag, &i, qrot, &iver);
     sscanf(molname, "%s\0", moln); 
   //    strncpy(mol[nline].dname,molname,16);
printf("select %6d %16s (<type> y to select it or n to skip it or t to terminate)\n", moltag, molname);
           scanf("%s", select);
         if(select[0]=='n')  printf("skip\n");
          
          if(select[0]=='y') {
                             printf("yes\n");
           mtag[nline] = moltag;
//          printf("%6d %s\n", mtag[nline], moln);
          strncat(mname, (char *)&moln, 16);
          strncat(mname, "\n", 1);
        nline++;
             } 
           if(select[0]=='t') {
                             printf("terminate selection\n");
                             moltag = 147001 ; 
                             }
        }
      *nmline=(char *)nline;
//         kp=0;      
//         for (j=0;j<nline; j++) {
//             k=0;
//             while(mname[kp]!='\n'&&k<16) {
//             moln[k] = mname[kp];
//             k++;
//             kp++;
//             }
//             moln[k]='\0';
//             kp++;
//           printf("\n");
//          printf("%6d moln=%s %d\n", mtag[j], moln, kp);
//           }
//         mname = (char *)&dump;
//        printf("nline=%d\n", nline);
}


void jpllinerd_c(fmx,fmn,strl,fnmol,fmoltag,mxnline,freq,intensity,uqst,lqst,mtag)
 float fmx,fmn,strl;
 int fnmol, fmoltag[fnmol];
 float *freq, *intensity;
 char *mxnline[];
 int *uqst, *lqst, *mtag;
 //char  ***uqst[][][], ***lqst[][][];
{
#define NLIM 10000
#define BIG  99999999.9999
#define CLIGHT 29979.2458
  int cl, i, itag, moltag, nline, iver, cnvt, ll, nmol, first, all;
  char *molname, *entval, *entnam, *catdir(), *fmakeword(), *makeword();
  char fqmin[14], fqmax[14], line[81], tbuf[22], chose[2];
  double fmin, fmax, strlim, f, err, str, stof();
  float outFREQ, outERR, outLGINT, outELO;
  int outDR, outGUP, outTAG, outQNFMT, outQN1[6], outQN2[6];
  char *qbuff, stc, qstr[81];
  int nn, j, ibuff;
  float qrot[7];
  int uqsta[6][NLIM], lqstb[6][NLIM];
/*    strlim = -500.; */
      
      cnvt = nmol = fqmin[0] = fqmax[0] = 0;
//      strlim = -500;
        strlim= strl;
       nline = 0;
//        fmin = 230.537;
//        fmax = 230.539;
        fmin = fmn * 1000;
        fmax = fmx * 1000;
//        printf("fmax=%f fmin=%f strl=%f fnmol=%d\n",
//            fmx, fmn, strl, fnmol);
         for (i=0; i< fnmol; i++){
//         printf("fmoltag=%d\n", fmoltag[i]);
                }
        sprintf(fqmin, "%13.4f", fmin);
        sprintf(fqmax, "%13.4f", fmax);
//        all = 1;
//        if (all) {
//          moltag = 0;
//          nxtdir(&moltag);
//        } else {
//          moltag = atoi(entval);
//               }
      
         for(itag=0; itag< fnmol; itag++) {
// while (moltag) {
          moltag=fmoltag[itag];
          first = 1;
          if ((ll = catfrq(moltag, fqmin, line)) > 0) {
            while (strcmp(fqmax, line) >= 0) {
              if (stof(line + 21, 8) > strlim) {
                if (first) {
                  molname = catdir(moltag, &i, qrot, &iver);
//jhz                  printf("             %6d %s\n", moltag, molname);
                  nmol++;
                           }
                nline++;
                first = 0;
//                puts(line);
//jhz                printf("%s \n size=%d\n" , line, sizeof(line));
                  nn=0;
// freq
                for (i =0; i <13; i++) {
                   qstr[i] = line[i];
                    }
                   qstr[13]='\0';
//                   printf("loop qstr[i]=%s\n",qstr);
                   sscanf(qstr, "%f", &outFREQ);
//                   printf("freq=%f\n", outFREQ);
                   freq[nline-1]= outFREQ;

// err
                  for (i =0; i <8; i++) {
                   qstr[i] = line[i+13];
                    }
                   qstr[8]='\0';
//                   printf("loop qstr[i]=%s\n",qstr);
                   sscanf(qstr, "%f", &outERR);
//                   printf("err=%f\n", outERR);
// LGINT
                   for (i =0; i <8; i++) {
                   qstr[i] = line[i+13+8];
                    }
                   qstr[8]='\0';
//                   printf("loop qstr[i]=%s\n",qstr);
                   sscanf(qstr, "%f\n", &outLGINT);
//                   printf("lgint=%f\n", outLGINT);
                   intensity[nline-1]=outLGINT;
// DR
                  for (i =0; i <2; i++) {
                   qstr[i] = line[i+13+8+8];
                    }
                   qstr[2]='\0';
//                   printf("loop qstr[i]=%s\n",qstr);
                   sscanf(qstr, "%d\n", &outDR);
//                   printf("lgint=%d\n", outDR);
// ELO
                  for (i =0; i <10; i++) {
                   qstr[i] = line[i+13+8+8+2];
                    }
                   qstr[10]='\0';
//                   printf("loop qstr[i]=%s\n",qstr);
                   sscanf(qstr, "%f\n", &outELO);
//                   printf("elo=%f\n", outELO);
// GUP 
                  for (i =0; i <3; i++) {
                   qstr[i] = line[i+13+8+8+2+10];
                    }
                   qstr[3]='\0';
//                   printf("loop qstr[i]=%s\n", qstr);
                   sscanf(qstr, "%d\n", &outGUP);
//                   printf("gup=%d\n", outGUP);
// TAG
                  for (i =0; i <7; i++) {
                   qstr[i] = line[i+13+8+8+2+10+3];
                    }
                   qstr[7]='\0';
//                   printf("loop qstr[i]=%s\n", qstr);
                   sscanf(qstr, "%d\n", &outTAG);
//                   printf("tag=%d\n", outTAG);
                       mtag[nline-1]=outTAG;
// QNFMT

                   for (i =0; i <4; i++) {
                   qstr[i] = line[i+13+8+8+2+10+3+7];
                    }
                   qstr[4]='\0';
//                   printf("loop qstr[i]=%s\n", qstr);
                   sscanf(qstr, "%d\n", &outQNFMT);
//                   printf("qnfmt=%d\n", outQNFMT);
//QN'
                   for  (j=0;j<6;j++) {
                   for  (i =0; i <2; i++){ 
                   qstr[i] = line[i+13+8+8+2+10+3+7+4+j*2];
                    }
                   qstr[2]='\0';
//                   printf("loop qstr[i]=%s\n", qstr);
                   sscanf(qstr, "%d\0", &outQN1[j]);
                   if(qstr[0]!=' '||qstr[1]!=' '){
//                 printf("qn1[%d]=%d\n",j+1, outQN1[j]);
                   uqsta[j][nline-1]=outQN1[j];
                   } else
                     {
//               printf("qn1[%d]=\n",j+1);
                       }
                   }
//QN''
                   for  (j=0; j<6;j++) {
                   for  (i =0; i <2; i++){
                   qstr[i] = line[i+13+8+8+2+10+3+7+4+12+j*2];
                    }
                   qstr[2]='\0';
//                   printf("loop qstr[i]=%s\n", qstr);
                   sscanf(qstr, "%d\0", &outQN2[j]);
                   if(qstr[0]!=' '||qstr[1]!=' ')
                     {
                   lqstb[j][nline-1]=outQN2[j];
//        printf("qn2[%d]=%d %d\n",j+1, lqstb[j][nline-1],outQN2[j]);
                     } else {
                     //printf("qn2[%d]=\n",j+1);
                           }
//		   fputs(line, stdout);
//jhz                  printf("\n");
                if (nline >= NLIM) {
                  printf("%s %d lines.%s\n",
                         "THIS form is currently limited to",
                         nline, " Please limit your search.");
                                   }
                                               }

                      }

              if (catrd(moltag, ++ll, line))
                break;
                                                  }
                                               }
          if (all) {
            if (nxtdir(&moltag) == 0) 
            goto cntinue;
          } else {
            moltag = 0;
          }
           
        }
    if (nmol == 0)
      puts("No lines were found for your search criteria.");
     cntinue:
     *mxnline = (char *)nline;
//     uqst = (int **)&uqsta;           
//     lqst = (int **)&lqstb;

               ibuff=0;
                   for(i=0; i<nline; i++){
                   for(j=0; j<6;     j++){
                  uqst[ibuff]=uqsta[j][i];
                  lqst[ibuff]=lqstb[j][i];
                   ibuff++;
                 //lqst[j][i] = lqstb[j][i];
                                    }   }

     for (i=0; i< nline;i++){
//jhz     printf("Mtag = %d\n", mtag[i]);
//jhz     printf("freq[%d]=%f inte=%f  ", i,freq[i],intensity[i]);
     for (j=0; j<6; j++) {
//jhz       printf("%d ", uqsta[j][i]);
                    }
//jhz       printf(" -> ");
       for (j=0; j<6; j++) {
//jhz      printf("%d ", lqstb[j][i]);
                    }
//jhz      printf("\n"); 
              }
//jhz      printf("end\n");
}

double stof(str, len)
char *str;
int len;
{
  int neg, dp, last, i;
  char c;
  double f;
  static double dpwr[6] = { 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001 };
  neg = 0;
  dp = len;
  last = -1;
  for (i = 0; i < len; i++) {
    c = str[i];
    if (c == '-') {
      neg = 1;
    } else if (c == '.') {
      dp = i;
    } else if (c >= '0' && c <= '9') {
      if (last >= 0) {
        f = f * 10. + (c - '0');
      } else {
        f = c - '0';
      }
      last = i;
    } else if (last >= 0) {
      break;
    }
  }
  if (last < 0)
    return 0.;
  if (neg)
    f = -f;
  dp = last - dp;
  while (dp) {
    i = dp;
    if (i > 6)
      i = 6;
    f = f * dpwr[i - 1];
    dp = dp - i;
  }
  return f;
}
#include <string.h>

#define NSPC 512
static int cat_nmol = 0;
static struct {
  int moltag, nlen, ver;
  float qln[7];
  char cname[16];
} catcom[NSPC + 1], *caterr, *catptr;

int catfrq(molec, cfreq, buff)
int molec;
char *cfreq, *buff;
{
  /*****************************************************************

    CATFRQ CALLS
    CATDIR CALLS CATRD
    CATRD CALLS CATFIL

    FUNCTION TO FIND AND RETURN FIRST LINE FOR SPECIES "MOLEC"
    THAT IS BEYOND FREQ

    CATFRQ <0 IF ERROR FOUND
    CATFRQ =0 IF FREQ GT FREQUENCY OF LAST LINE
    CATFRQ >0 POINTS TO LINE

    BUFF CONTAINS THE LINE INFORMATION

    CFREQ IS CHARACTER STRING IN FORMAT (F13.4) OR EQUIVALENT

  ******************************************************************/
  int line, nline, l, r;
  nline = catlen(molec);
  /*     lines  < L have frequency  < FREQ */
  /*     lines >= R have frequency >= FREQ */
  l = line = 1;
  r = nline + 1;
  buff[0] = 0;
  while (l < r) {
    line = (l + r) >> 1;
    if (catrd(molec,line,buff))
      return (-1);
    if (strncmp(cfreq,buff, 13) >= 0) {
      l = line + 1;
    } else {
      r = line;
    }
  }
  if (l > nline) {
    buff[0] = 0;
    return (0);
  }
  if (r == line)
    return (line);
  line = r;
  if (catrd(molec, line, buff))
    return (-1);
  return (line);
}

int catrd(molec, line, buff)
int molec, line;
char *buff;
/***********************************************************************

     CATRD CALLS
       CATFIL


     SUBROUTINE TO READ CATALOGUE FILE FOR MOLECULE "MOLEC",LINE# "LINE"

       80 CHARACTER BUFFER RETURNED IN BUFF

     ERROR CODE RETURNED (IOSTAT ERROR #)

*************************************************************************/
{
  static int buflen = 80;
  static int omolec = -1;
  static FILE *handle;
  long offset;
  if (molec != omolec) {        /* new molecule */
    if (omolec > 0)
      fclose(handle);
    omolec = molec;
    handle = fopen(catfil(molec), "r");
    if (handle == NULL) {
      buff[0] = 0;
      omolec = -1;
      return (-1);
    }
  }
  if (line <= 0)
    return (-2);
  offset = (line - 1) * (long) (buflen + 1);
  fseek(handle, offset, SEEK_SET);
  if (fread(buff, 1, buflen, handle) == buflen) {
    buff[buflen] = 0;
    return (0);
  } else {
    buff[0] = 0;
    return (1);
  }
}

char *catfil(num)
int num;
/**********************************************************
     FUNCTION WHICH RETURNS FILE NAME CORRESPONDING TO NUM
***********************************************************/
{
/*  static char catdir[] = { "/catalog/catdir.cat" };
  static char catent[] = { "/catalog/c000000.cat" };
                            1234567890123456
/home/miriad/cat/c000000.cat
123456789012345678901234               
../../cat/jpl/c000000.cat
123456789112345678921
*/
//  static char catdir[] = { "../../cat/jpl/catdir.cat" };
//  static char catent[] = { "../../cat/jpl/c000000.cat" };
    static char catdir[80+10];
    static char catent[80+11];
  char *cfield;
  int k, iq;
//          printf("pathlen=%d jplpath=%s\n", length, path);
  if (num == 0) {
      strncpy(catdir,path,length);
      strncpy(catent,path,length);
      strncat(catdir,"catdir.cat",10);
      strncat(catent,"c000000.cat",11);
//      printf("catdir %s catent %s\n", catdir, catent);
    return catdir;}

//  cfield = &catent[16];
//  cfield = &catent[24];
      cfield = &catent[length+7];
  for (k = 0; k < 6; ++k) {
    iq = num;
    num = num / 10;
    --cfield;
    *cfield = (char) ('0' + (iq - num * 10));
  }
  return catent;
}

int catlen(molec)
int molec;
/**************************************************
C
C   SUBROUTINE TO RETURN CATALOGUE ENTRY LENGTH
C
C   MOLEC IS SPECIES TAG
C
****************************************************/
{
  static int fmt[9] = { 6, 7, 7, 7, 7, 7, 7, 7, 2 };
  FILE *cdir;
  char *pbuf;
  double dval[9];
  char buf[82];
  float *qptr;
  int k;
  if (cat_nmol == 0) {          /* initialize */
    cdir = fopen(catfil(0), "r");
    if (cdir == NULL)
      return (0);
    catptr = catcom;
    while (cat_nmol < NSPC) {
      if (fgets(buf, 82, cdir) == NULL)
        break;
      pcard(buf, dval, 1, fmt);
      catptr->moltag = (int) dval[0];
      pbuf = catptr->cname;
      memcpy(pbuf, buf + 6, 14);
      pbuf[14] = 0;
      pcard(buf + 20, dval, 9, fmt);
      catptr->nlen = (int) dval[0];
      if (catptr->moltag == 0 || catptr->nlen == 0)
        continue;
      qptr = catptr->qln;
      for (k = 0; k < 7; k++)
        qptr[k] = (float) dval[k + 1];
      catptr->ver = (int) dval[8];
      ++catptr;
      ++cat_nmol;
    }
    fclose(cdir);
    caterr = catptr;
    strcpy(caterr->cname, " error");
    qptr = caterr->qln;
    for (k = 0; k < 7; k++)
      qptr[k] = 0.;
    caterr->moltag = 0;
    caterr->nlen = 0;
    caterr->ver = 0;
  }
  if (molec > 0) {
    for (k = 0; k < cat_nmol; ++k) {
      if (catptr == caterr)
        catptr = catcom;
      if (molec == catptr->moltag)
        return catptr->nlen;
      ++catptr;
    }
  }
  catptr = caterr;
  return 0;
}

char *catdir(molec, nline, qqln, iver)
int molec, *nline, *iver;
double *qqln;
/*********************************************************
C   SUBROUTINE TO RETURN CATALOGUE DIRECTORY INFORMATION
C
C   MOLEC IS SPECIES TAG
C   CATDIR IS ASCII SPECIES NAME
C   NLINE IS THE NUMBER OF LINES FOR MOLECULE
C   QLN IS THE LOG10 OF THE PARTITION FUNCTION FOR
C        300,225,150,75,37.5,18.75,9.375 K
C   IVER IS THE VERSION NUMBER
*************************************************************/
{
  int k;
  float *qptr;
  *nline = catlen(molec);
  qptr = catptr->qln;
  for (k = 0; k < 7; k++)
    qqln[k] = qptr[k];
  *iver = catptr->ver;
  return catptr->cname;
}

int nxtdir(molec)
int *molec;
/**********************************************************************
C
C     FUNCTION NXTDIR RETURNS THE NUMBER OF REMAINING DIRECTORY ENTRIES
C        AFTER INPUT SPECIES MOLEC
C     ON RETURN MOLEC IS CHANGED TO THE NEXT SPECIES TAG
C
C     IF INPUT SPECIES MOLEC = 0 THEN POSITION IS BEFORE FIRST ENTRY
C
C     IF INPUT SPECIES MOLEC = LAST ENTRY THEN MOLEC=0 ON RETURN
***********************************************************************/
{
  catlen(*molec);
  if (*molec == 0)
    catptr = catcom;
  else if (catptr != caterr)
    ++catptr;
  *molec = catptr->moltag;
  return catptr->nlen;
}

int getcat(buf, pdata)
char *buf;
struct catdata *pdata;
{
  static double dval[8];
  static int fmt[8] = { 13, 8, 8, 2, 10, 3, 7, 4 };
  if (pcard(buf, dval, 8, fmt) < 8)
    return -1;
  pdata->freq = dval[0];
  pdata->derr = dval[1];
  pdata->str = dval[2];
  pdata->itd = (int) dval[3];
  pdata->elow = dval[4];
  pdata->igup = (int) dval[5];
  pdata->tag = (int) dval[6];
  pdata->ifmt = (int) dval[7];
  return (readqn(buf + 55, pdata->iqn, 12));
}
#include <limits.h>
#if INT_MAX > 0x7fff
#define MAXDEC 9
#else
#define MAXDEC 4
#endif
#define TRUE  1
#define FALSE 0
int readqn(qnstr, iqn, n)
const char *qnstr;
short *iqn;
const int n;
{ /* read n quanta from string to iqn */
  const char *pqn;
  int i, ich, ic, ival;

  pqn = qnstr;
  for (i = 0; i < n; ++i) {
    if (*pqn == 0)
      break;
    ich = (*pqn & 0xff);
    ++pqn;
    ic = (*pqn & 0xff);
    ++pqn;
    if (ic == ' ') {
      iqn[i] = 0;
      continue;
    }
    ival = ic - '0';
    if (ival < 0 || ival > 9)
      break;
    if (ich != ' ') {
      if (ich == '-') {
        ival = -ival;
      } else if (ich >= '0' && ich <= '9') {
        ival += (ich - '0') * 10;
      } else if (ich >= 'a' && ich <= 'z') {
        ival = -ival - 10 * (ich - ('a' - 1));
      } else if (ich >= 'A' && ich <= 'Z') {
        ival += 10 * (ich - ('A' - 10));
      } else {
        break;
      }
    }
    iqn[i] = (short) ival;
  }
  return (n - i);
}                               /* readqn */

int pcard(card, val, nval, fmtlen)
const char *card;
double *val;
const int nval;
const int *fmtlen;
{                               /* parses a card for NVAL numbers */
  /*
   * ndec = -2 for comma detected, ndec = -1 indicates looking for 
   * beginning of field, ndec >= 0 indicates in field 
   * pwrflg = 0 indicates mantissa decimal point not detected 
   * pwrflg > 0 indicates mantissa decimal point detected 
   * pwrflg = -1 indicates e, E, d, or D found 
   * pwrflg = -2 indicates in exponent field 
   */

  static double pten[] =
      { 1., 10., 100., 1000., 1e4, 1e5, 1e6, 1e7, 1e8, 1e9 };
  const char *newfield;
  double tmp, fac;
  int ndec, kval, itmp, ipwr, pwrflg, ich, npwr;
  int newnum, neg, negpwr;      /* boolean */

  newfield = (char *) fmtlen;
  if (newfield)
    newfield = card + fmtlen[0];
  newnum = TRUE;
  neg = negpwr = FALSE;
  itmp = ipwr = pwrflg = kval = npwr = 0;
  tmp = fac = 0.;
  ndec = -1;
  while (kval < nval) {
    ich = (*card & 0xff);
    if (ich != 0) {
      if (card == newfield) {
        ich = ',';
      } else {
        card++;
      }
    }
    if (ich >= '0' && ich <= '9') {     /* character is a number */
      if (pwrflg == 0) {
        ++npwr;                 /* count integer part of mantissa */
      } else if (pwrflg > 0) {
        --ipwr;                 /* count fraction part of mantissa */
      } else {                  /* pwrflg < 0 */
        pwrflg = -2;            /* flag indicates digit found in exponent */
      }
      ich -= '0';
      if (ndec <= 0) {
        itmp = ich;
        ndec = 1;               /* first digit */
      } else if (ndec == MAXDEC) {
        if (pwrflg < 0)
          return kval;          /* exponent field is too big */
        /*
         * now is a good time to convert integer to real 
         */
        if (newnum) {
          tmp = (double) itmp;
          newnum = FALSE;
        } else {
          tmp = tmp * pten[ndec] + itmp;
        }
        itmp = ich;
        ndec = 1;
      } else {                  /* accumulate integer */
        ++ndec;
        itmp = itmp * 10 + ich;
      }
    } else if (ich == '.' && pwrflg == 0) {     /* first decimal point */
      if (ndec < 0)
        ndec = 0;
      pwrflg = 1;
    } else if (ich == '-' && pwrflg == -1) {    /* leading - in exponent
                                                 * field */
      negpwr = TRUE;
      ndec = 0;
      pwrflg = -2;
    } else if (ich == '+' && pwrflg == -1) {    /* leading + in exponent
                                                 * field */
      ndec = 0;
      pwrflg = -2;
    } else {                    /* character is not a number or decimal
                                 * point, E+, E- */
      if (ndec >= 0) {          /* save results from number decoding */
        if (pwrflg < 0) {       /* integer follows 'E' */
          pwrflg = 0;
          if (negpwr)
            itmp = -itmp;
          npwr += itmp;
          ipwr += itmp;
        } else {                /* finish up mantissa */
          pwrflg = 0;
          if (newnum) {
            tmp = (double) itmp;
          } else {
            tmp = itmp + tmp * pten[ndec];
          }
          if (ich == 'E' || ich == 'e' || ich == 'D' || ich == 'd') {
            /*
             * look for exponent 
             */
            pwrflg = -1;
            ndec = 0;
            negpwr = FALSE;
          }
        }
        if (pwrflg == 0) {      /* number finished */
          if (npwr < -38) {
            ipwr = 0;
            tmp = 0.;
          } else if (npwr > 37) {
            ipwr = 0;
            tmp = 1.e+37;
          }
          if (ipwr != 0) {      /* scale by powers of 10 */
            if (ipwr < 0) {
              ipwr = -ipwr;
              if (ipwr > 7)
                fac = 1. / pten[8];
              itmp = ipwr & 7;
              if (itmp != 0)
                tmp /= pten[itmp];
            } else {
              fac = pten[8];
              itmp = ipwr & 7;
              if (itmp != 0)
                tmp *= pten[itmp];
            }
            if (ipwr > 7) {
              ipwr = ipwr >> 3;
              if ((ipwr & 1) != 0)
                tmp *= fac;
              while (ipwr > 1) {
                ipwr = ipwr >> 1;
                fac = fac * fac;
                if ((ipwr & 1) != 0)
                  tmp *= fac;
              }
            }
            ipwr = 0;
          }                     /* end scale */
          if (neg)
            tmp = -tmp;
          val[kval] = tmp;
          ++kval;
          if (newfield && kval < nval)
            newfield += fmtlen[kval];
          ndec = -1;
          npwr = 0;
          neg = FALSE;
        }
        itmp = 0;
        newnum = TRUE;
      }                         /* finished save results */
      if (ich == 0)
        return kval;
      /*
       * check for delimiters in new field 
       */
      if (ich == '.') {         /* decimal afer exponent field, assume new 
                                 */
        ndec = 0;
        pwrflg = 1;
      } else if (ich == '-') {  /* leading - in mantissa field */
        ndec = 0;
        neg = TRUE;
      } else if (ich == '+') {  /* leading + in mantissa field */
        ndec = 0;
        neg = FALSE;
      } else if (ich == ',') {
        if (ndec == -2) {       /* second comma */
          ++kval;
          if (newfield && kval < nval)
            newfield += fmtlen[kval];
        }
        ndec = -2;
      } else if (ich == '/') {  /* end of line character */
        kval = nval;
      }
    }
  }
  return kval;
} /* pcard */

#define LF 10
#define CR 13

void getword(word, line, stop)
char *word, *line, stop;
{
  int x = 0, y;

  for (x = 0; ((line[x]) && (line[x] != stop)); x++)
    word[x] = line[x];

  word[x] = '\0';
  if (line[x])
    ++x;
  y = 0;

  while (line[y++] = line[x++]);
}

char *makeword(line, stop)
char *line, stop;
{
  int x = 0, y;
  char *word = (char *) malloc(sizeof(char) * (strlen(line) + 1));

  for (x = 0; ((line[x]) && (line[x] != stop)); x++)
    word[x] = line[x];

  word[x] = '\0';
  if (line[x])
    ++x;
  y = 0;

  while (line[y++] = line[x++]);
  return word;
}

char *fmakeword(f, stop, cl)
FILE *f;
char stop;
int *cl;
{
  int wsize;
  char *word;
  int ll;

  wsize = 102400;
  ll = 0;
  word = (char *) malloc(sizeof(char) * (wsize + 1));

  while (1) {
    word[ll] = (char) fgetc(f);
    if (ll == wsize) {
      word[ll + 1] = '\0';
      wsize += 102400;
      word = (char *) realloc(word, sizeof(char) * (wsize + 1));
    }
    --(*cl);
    if ((word[ll] == stop) || (feof(f)) || (!(*cl))) {
      if (word[ll] != stop)
        ll++;
      word[ll] = '\0';
      return word;
    }
    ++ll;
  }
}

char x2c(what)
char *what;
{
  register char digit;

  digit =
      (what[0] >= 'A' ? ((what[0] & 0xdf) - 'A') + 10 : (what[0] - '0'));
  digit *= 16;
  digit +=
      (what[1] >= 'A' ? ((what[1] & 0xdf) - 'A') + 10 : (what[1] - '0'));
  return (digit);
}

unescape_url(url)
char *url;
{
  register int x, y;

  for (x = 0, y = 0; url[y]; ++x, ++y) {
    if ((url[x] = url[y]) == '%') {
      url[x] = x2c(&url[y + 1]);
      y += 2;
    }
  }
  url[x] = '\0';
}

plustospace(str)
char *str;
{
  register int x;

  for (x = 0; str[x]; x++)
    if (str[x] == '+')
      str[x] = ' ';
}

int rind(s, c)
char *s, c;
{
  register int x;
  for (x = strlen(s) - 1; x != -1; x--)
    if (s[x] == c)
      return x;
  return -1;
}

int getline(s, n, f)
char *s;
int n;
FILE *f;
{
  register int i = 0;

  while (1) {
    s[i] = (char) fgetc(f);

    if (s[i] == CR)
      s[i] = fgetc(f);

    if ((s[i] == 0x4) || (s[i] == LF) || (i == (n - 1))) {
      s[i] = '\0';
      return (feof(f) ? 1 : 0);
    }
    ++i;
  }
}

send_fd(f, fd)
FILE *f, *fd;
{
  int num_chars = 0;
  char c;

  while (1) {
    c = fgetc(f);
    if (feof(f))
      return;
    fputc(c, fd);
  }
}

int ind(s, c)
char *s, c;
{
  register int x;

  for (x = 0; s[x]; x++)
    if (s[x] == c)
      return x;

  return -1;
}

escape_shell_cmd(cmd)
char *cmd;
{
  register int x, y, l;

  l = strlen(cmd);
  for (x = 0; cmd[x]; x++) {
    if (ind("&;`'\"|*?~<>^()[]{}$\\\n", cmd[x]) != -1) {
      for (y = l + 1; y > x; y--)
        cmd[y] = cmd[y - 1];
      l++;                      /* length has been increased */
      cmd[x] = '\\';
      x++;                      /* skip the character */
    }
  }
}
