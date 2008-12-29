/*  jhz 2004-nov-29: collect the c subroutines for libsma2miriad*/
/****************************************************************/
/* .c
   jhz- 2004-jun-15: swapendian subs used in converting
   mir data to miriad format in linux boxes
 
   pjt  2004-dec-10: renamed data_write.h to sma_data.h for MIR 4.0.4
   pjt  2005-may-23: slight re-indendation to align with updated sma code
                     added warning about sizeof(long int) = 4 vs. 8
		     added some assert to make the code fail if this happens

 
   TODO:   <string.h>::bcopy -> memcpy
 */

/* #dfin SWAP_ENDIAN 1  for cfa0 (big endian lf-o-righ incra
    in addr numbr - Moorola 680x0  */
/* #dfin SWAP_ENDIAN 0  for id (inl)  buzz (alpha) and linux (
   lil ndian righ-o-lf Inl 80x86 and Pnium procor */

#include <stdio.h>
#include <assert.h>
#include "sma_data.h"

int check_s2 = sizeof(short);
int check_s4 = sizeof(int);
int check_s8 = sizeof(long);

void
reverse4(char *bytes)
{
   char t;
   t = *bytes;
   *bytes = bytes[3];
   bytes[3] = t;
   t = bytes[1];
   bytes[1] = bytes[2];
   bytes[2] = t;
}

void
reverse8(char *bytes)
{
   char t;

   t = *bytes;
   *bytes = bytes[7];
   bytes[7] = t;
   t = bytes[1];
   bytes[1] = bytes[6];
   bytes[6] = t;
   t = bytes[2];
   bytes[2] = bytes[5];
   bytes[5] = t;
   t = bytes[3];
   bytes[3] = bytes[4];
   bytes[4] = t;

}
void
reverse2(char *bytes)
{
   char t;
   t = *bytes;
   *bytes = bytes[1];
   bytes[1] = t;
}

void
reverse1(char *bytes)
{
   char t;
   t = *bytes;
   *bytes = bytes[0];
}

struct inh_def *swap_inh(struct inh_def *inh_pnr)
{     struct inh_def inh_buff;
      bcopy(inh_pnr, &inh_buff, sizeof(struct inh_def));
      reverse4((char *)(&inh_buff.conid));
      reverse2((char *)(&inh_buff.icocd));
      reverse4((char *)(&inh_buff.traid));
      reverse4((char *)(&inh_buff.inhid));
      reverse4((char *)(&inh_buff.ints));
      reverse2((char *)(&inh_buff.itq));
      reverse4((char *)(&inh_buff.az));
      reverse4((char *)(&inh_buff.el));
      reverse4((char *)(&inh_buff.ha));
      reverse2((char *)(&inh_buff.iut));
      reverse2((char *)(&inh_buff.iref_time));
      reverse8((char *)(&inh_buff.dhrs)); 
      reverse4((char *)(&inh_buff.vc));
      reverse2((char *)(&inh_buff.ivctype));
      reverse8((char *)(&inh_buff.sx));
      reverse8((char *)(&inh_buff.sy));
      reverse8((char *)(&inh_buff.sz));
      reverse4((char *)(&inh_buff.rinteg));
      reverse4((char *)(&inh_buff.proid));
      reverse4((char *)(&inh_buff.souid));
      reverse2((char *)(&inh_buff.isource));
      reverse2((char *)(&inh_buff.ipos));
      reverse4((char *)(&inh_buff.offx));
      reverse4((char *)(&inh_buff.offy));
      reverse2((char *)(&inh_buff.iofftype));
      reverse2((char *)(&inh_buff.ira));
      reverse2((char *)(&inh_buff.idec));
      reverse8((char *)(&inh_buff.rar));
      reverse8((char *)(&inh_buff.decr));
      reverse4((char *)(&inh_buff.epoch));
      reverse4((char *)(&inh_buff.sflux));
      reverse4((char *)(&inh_buff.size));
      bcopy(&inh_buff,inh_pnr, sizeof(struct inh_def));
      return(inh_pnr);
}

struct blh_def *swap_blh(struct blh_def *blh_pnr)
{     struct blh_def blh_buff;
      bcopy(blh_pnr, &blh_buff, sizeof(struct blh_def));
      reverse4((char *)(&blh_buff.blhid));
      reverse4((char *)(&blh_buff.inhid));
      reverse2((char *)(&blh_buff.isb));
      reverse2((char *)(&blh_buff.ipol));
      reverse4((char *)(&blh_buff.pa));
      reverse2((char *)(&blh_buff.iaq));
      reverse2((char *)(&blh_buff.ibq));
      reverse2((char *)(&blh_buff.icq));
      reverse2((char *)(&blh_buff.ioq));
      reverse2((char *)(&blh_buff.irec));
      reverse2((char *)(&blh_buff.iifc));
      reverse4((char *)(&blh_buff.u));
      reverse4((char *)(&blh_buff.v));
      reverse4((char *)(&blh_buff.w));
      reverse4((char *)(&blh_buff.prbl));
      reverse4((char *)(&blh_buff.angres));
      reverse4((char *)(&blh_buff.vis));
      reverse4((char *)(&blh_buff.coh));
      reverse4((char *)(&blh_buff.sigcoh));
      reverse4((char *)(&blh_buff.csnr));
      reverse4((char *)(&blh_buff.vflux));
      reverse4((char *)(&blh_buff.cnoise));
      reverse8((char *)(&blh_buff.avedhrs));
      reverse4((char *)(&blh_buff.ampave));
      reverse4((char *)(&blh_buff.phaave));
      reverse4((char *)(&blh_buff.tpvar));
      reverse4((char *)(&blh_buff.blsid));
      reverse2((char *)(&blh_buff.itel1));
      reverse2((char *)(&blh_buff.itel2));
      reverse2((char *)(&blh_buff.iblcd));
      reverse4((char *)(&blh_buff.ble));
      reverse4((char *)(&blh_buff.bln));
      reverse4((char *)(&blh_buff.blu));
      reverse4((char *)(&blh_buff.soid));
      bcopy(&blh_buff,blh_pnr, sizeof(struct blh_def));
return(blh_pnr);
}

struct sph_def *swap_sph(struct sph_def *sph_pnr)
{    struct sph_def sph_buff;
     bcopy(sph_pnr, &sph_buff, sizeof(struct sph_def));
     reverse4((char *)(&sph_buff.sphid));
     reverse4((char *)(&sph_buff.blhid));
     reverse4((char *)(&sph_buff.inhid));
     reverse2((char *)(&sph_buff.igq));
     reverse2((char *)(&sph_buff.ipq));
     reverse2((char *)(&sph_buff.iband));
     reverse2((char *)(&sph_buff.ipstate));
     reverse4((char *)(&sph_buff.tau0));
     reverse8((char *)(&sph_buff.vel));
     reverse4((char *)(&sph_buff.vres));
     reverse2((char *)(&sph_buff.ivtype));
     reverse8((char *)(&sph_buff.fsky));
     reverse4((char *)(&sph_buff.fres));
     reverse4((char *)(&sph_buff.tssb));
     reverse4((char *)(&sph_buff.integ));
     reverse4((char *)(&sph_buff.wt));
     reverse2((char *)(&sph_buff.itaper));
     reverse4((char *)(&sph_buff.snoise));
     reverse2((char *)(&sph_buff.nch));
     reverse2((char *)(&sph_buff.nrec));
     reverse4((char *)(&sph_buff.dataoff));
     reverse4((char *)(&sph_buff.linid));
     reverse2((char *)(&sph_buff.itrans));
     reverse8((char *)(&sph_buff.rfreq));
     reverse2((char *)(&sph_buff.pasid));
     reverse2((char *)(&sph_buff.gaiidamp));
     reverse2((char *)(&sph_buff.gaiidpha));
     reverse2((char *)(&sph_buff.flcid));
     reverse2((char *)(&sph_buff.atmid));
     bcopy(&sph_buff,sph_pnr, sizeof(struct sph_def));
         return(sph_pnr);
}

struct codeh_def *swap_cdh(struct codeh_def *cdh_pnr)
{    struct codeh_def cdh_buff;
     int i;
     bcopy(cdh_pnr, &cdh_buff, sizeof(struct codeh_def));
     for (i=0; i< 12; i++) {
       reverse1((char *)(&cdh_buff.v_name[i]));
     }
     
     reverse2((char *)(&cdh_buff.icode));
      
     for (i=0; i< 26; i++)     {
       reverse1((char *)(&cdh_buff.code[i]));
     }
     reverse2((char *)(&cdh_buff.ncode));

     bcopy(&cdh_buff,cdh_pnr, sizeof(struct codeh_def));
     return(cdh_pnr);
}
 
struct ant_def *swap_enh(struct ant_def *enh_pnr)
{    struct ant_def enh_buff;
     bcopy(enh_pnr, &enh_buff, sizeof(struct ant_def));
     reverse4((char *)(&enh_buff.antennaNumber));
     reverse4((char *)(&enh_buff.padNumber));
     reverse4((char *)(&enh_buff.antennaStatus));
     reverse4((char *)(&enh_buff.trackStatus));
     reverse4((char *)(&enh_buff.commStatus));
     reverse4((char *)(&enh_buff.inhid));
     reverse4((char *)(&enh_buff.ints));
     reverse8((char *)(&enh_buff.dhrs));
     reverse8((char *)(&enh_buff.ha));
     reverse8((char *)(&enh_buff.lst));
     reverse8((char *)(&enh_buff.pmdaz));
     reverse8((char *)(&enh_buff.pmdel));
     reverse8((char *)(&enh_buff.tiltx));
     reverse8((char *)(&enh_buff.tilty));
     reverse8((char *)(&enh_buff.actual_az));
     reverse8((char *)(&enh_buff.actual_el));
     reverse8((char *)(&enh_buff.azoff));
     reverse8((char *)(&enh_buff.eloff));
     reverse8((char *)(&enh_buff.az_tracking_error));
     reverse8((char *)(&enh_buff.el_tracking_error));
     reverse8((char *)(&enh_buff.refraction));
     reverse8((char *)(&enh_buff.chopper_x));
     reverse8((char *)(&enh_buff.chopper_y));
     reverse8((char *)(&enh_buff.chopper_z));
     reverse8((char *)(&enh_buff.chopper_angle));
     reverse8((char *)(&enh_buff.tsys));
     reverse8((char *)(&enh_buff.ambient_load_temperature));
     bcopy(&enh_buff, enh_pnr, sizeof(struct ant_def));
     return(enh_pnr);
}

struct sch_def *swap_sch(struct sch_def *sch_pnr)
{   struct sch_def sch_buff;
    bcopy(sch_pnr, &sch_buff, sizeof(struct sch_def));
    reverse4((char *)(&sch_buff.inhid));
    reverse1((char *)(&sch_buff.form[0]));
    reverse1((char *)(&sch_buff.form[1]));
    reverse1((char *)(&sch_buff.form[2]));
    reverse1((char *)(&sch_buff.form[3]));
    reverse4((char *)(&sch_buff.nbyt));
    reverse4((char *)(&sch_buff.nbyt_pack));
    bcopy(&sch_buff,sch_pnr, sizeof(struct sch_def));
    return(sch_pnr);
}

short int *swap_sch_data(short int *sch_data_pnr, int datalength)
{    short int sch_data_buff;
     int i;
     assert(check_s2==2);
     for (i=0; i<datalength; i++) {
       sch_data_buff=sch_data_pnr[i];
       reverse2((char *)(&sch_data_buff));
       sch_data_pnr[i]=sch_data_buff;
     }
     return(sch_data_pnr);
}

int *swap_int_data(int *int_data_pnr, int datalength)
{    /* datalength: number of array elements */
     int int_data_buff;
     int i;
     assert(check_s4==4);
     printf("before swap:  i pnr buff :--: before swap: pnr buff\n");
     for (i=0; i<datalength; i++) {
       int_data_buff=int_data_pnr[i];     
       printf("%d %d %d :--:  ", i,
	      int_data_pnr[i], int_data_buff);
       reverse4((char *)(&int_data_buff));
       int_data_pnr[i]=int_data_buff;
       printf("%d %d\n",
	      int_data_pnr[i], int_data_buff);

     }
     return(int_data_pnr);
}

/* WARNING: long int is 4 bytes on IA-32 and 8 bytes on IA-64 */

long int *swap_long_data(long int *long_data_pnr, int datalength)
{    /* datalength: number of array elements */
     long int long_data_buff;
     int i;
     assert(check_s8==8);
     printf("long before swap:  i pnr buff :--: before swap: pnr buff %d\n",
            datalength);
     for (i=0; i<datalength; i++) {
       long_data_buff = long_data_pnr[i];
       printf("%d %d %d :--:  ", i,
	      long_data_pnr[i], long_data_buff);
       reverse8((char *)(&long_data_buff));
       long_data_pnr[i] = long_data_buff;
       printf("%d %d\n",
	      long_data_pnr[i], long_data_buff);
       
     }
     return(long_data_pnr);
}

