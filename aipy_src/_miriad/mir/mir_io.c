/************************************************************************/
/*       MIR -- routines for writing a MIR dataset                      */
/*                                                                      */
/*  History:                                                            */
/*    dnf  16mar06  initial version                                     */
/************************************************************************/

#define CHECK(x,a) if(x) {  bug_c('w',((void)sprintf a,message)); \
			    bugno_c('f',x);			 \
			 }
#define NUMFLS 5        //number of outpur MIR files
#include <stdio.h>
#include <string.h>
#include "miriad.h"

static char message[128];
static FILE* ofls[NUMFLS];

//must convert to shorts here as it is not easy to pass integer*2 from Fortran to C short
short i2s(const int num);

/************************************************************************/
void mirInit_c(const char *f_name)
/**mirInit -- open the output files for a MIR data set                  */
/*&dnf                                                                  */
/*low-level-i/o                                                         */
/*+ FORTRAN call sequence:

         subroutine mirInit(f_name)
         character f_name*(*)

  Create and open the output files for a MIR dataset.

  Input:
    f_name     Name of the parent directory where the files will be 
               written                                                  */
/*--                                                                    */
/*----------------------------------------------------------------------*/
{
    char path[128];
    char location[128];
    char filename[NUMFLS][128];
    int i,iostat;
    strcpy(filename[0],"in_read");
    strcpy(filename[1],"bl_read");
    strcpy(filename[2],"sp_read");
    strcpy(filename[3],"codes_read");
    strcpy(filename[4],"vis_data");

    strcpy(path,f_name);
    dmkdir_c(path,&iostat);
    CHECK(iostat,(message,"Error opening %s, in mirInit",path));
    for(i=0; i<NUMFLS; i++){
	strcpy(location,path);
	strcat(location,filename[i]);
	ofls[i]=fopen(location,"w");
	if(ofls[i] == NULL){
	    bug_c('f',"Cannot open output file in mirInit");
	}
    }
}

/************************************************************************/
void mirClose_c(){
/**mirClose -- Close a MIR dataset                                      */
/*&dnf                                                                  */
/*:low-level-i/o                                                        */
/*+ FORTRAN call sequence

       subroutine mirClose

  This closes a MIR data set                                            */
/*--                                                                    */
/*----------------------------------------------------------------------*/
    int i;
    for(i=0;i<NUMFLS;i++){
	fclose(ofls[i]);
    }
}

/************************************************************************/
void inWrite_c(const int conid, const int icocd, const int traid, const int inhid, const int ints, const int itq, const float az, const float el, const float ha, const int iut, const int iref_time, const double dhrs, const float vc, const int ivctype, const double sx, const double sy, const double sz, const float rinteg, const int proid, const int souid, const int isource, const int ipos, const float offx, const float offy, const int iofftype, const int ira, const int idec, const double rar, const double decr, const float epoch, const float sflux, const float size){
/**inWrite -- write the integration header for a MIR data set           */
/*&dnf                                                                  */
/*:low-level-i/o                                                        */
/*+ FORTRAN call sequence

       subroutine inWrite(conid,icocd,traid,inhid,ints,itq,az,el,ha,
                          iut,ireftime,dhrs,vc,ivctype,sx,sy,sz,rinteg,
                          proid,souid,isource,ipos,offx,offy,iofftype,
                          ira,idec,rar,decr,epoch,sflux,size)
       integer conid,icocd,traid,inhid,ints,itq,iut,ireftime,ivctype
               proid, souid,isource,ipos,iofftype,ira,idec
       real az,el,ha,vc,rinteg,offx,offy,epoch,sflux,size
       double precision dhrs,sx,sy,sz,rar,decr

  This writes an entry into a MIR integration header

  Input:
    conid    configuration id number
    icocd    configuration code number
    traid    track id
    inhid    integration id (globally unique number)
    int      integration number
    itq      tuning code
    az       average antenna azimuth (degrees)
    el       average antenna elevation (degrees)
    ha       hour angle
    iut      ut code
    ireftime reference time code
    dhrs     hours from reference time (hrs)
    vc       velocity offset
    ivctype  type of offset
    sx       directional cosine to source
    sy       directional cosine to source
    sz       directional cosine to source
    rinteg   integration time
    proid    projeci id
    souid    source id
    isource  source code
    ipos     position code
    offx     offset in x (dra)
    offy     offset in y (ddec)
    iofftype offset type 'eq' for MIRIAD data
    ira      ra code
    idec     dec code
    rar      ra (radians)
    decr     dec (radians)
    epoch    epoch of coordinates
    sflux    flux of planet (if applicable)
    size     size of planet (if applicable)                             */
/*--                                                                    */
/*----------------------------------------------------------------------*/
    int nobj=0;
    short svar;

    nobj += fwrite(&conid,sizeof(conid),1,ofls[0]);
    if(nobj == 0){
	bug_c('f',"Unable to write to in_read.");	
    }
    svar=i2s(icocd);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[0]);
    nobj += fwrite(&traid,sizeof(traid),1,ofls[0]);
    nobj += fwrite(&inhid,sizeof(inhid),1,ofls[0]);
    nobj += fwrite(&ints,sizeof(ints),1,ofls[0]);
    svar=i2s(itq);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[0]);
    nobj += fwrite(&az,sizeof(az),1,ofls[0]);
    nobj += fwrite(&el,sizeof(el),1,ofls[0]);
    nobj += fwrite(&ha,sizeof(ha),1,ofls[0]);
    svar=i2s(iut);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[0]);
    svar=i2s(iref_time);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[0]);
    nobj += fwrite(&dhrs,sizeof(dhrs),1,ofls[0]);
    nobj += fwrite(&vc,sizeof(vc),1,ofls[0]);
    svar=i2s(ivctype);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[0]);
    nobj += fwrite(&sx,sizeof(sx),1,ofls[0]);
    nobj += fwrite(&sy,sizeof(sy),1,ofls[0]);
    nobj += fwrite(&sz,sizeof(sz),1,ofls[0]);
    nobj += fwrite(&rinteg,sizeof(rinteg),1,ofls[0]);
    nobj += fwrite(&proid,sizeof(proid),1,ofls[0]);
    nobj += fwrite(&souid,sizeof(souid),1,ofls[0]);
    svar=i2s(isource);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[0]);
    svar=i2s(ipos);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[0]);
    nobj += fwrite(&offx,sizeof(offx),1,ofls[0]);
    nobj += fwrite(&offy,sizeof(offy),1,ofls[0]);
    svar=i2s(iofftype);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[0]);
    svar=i2s(ira);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[0]);
    svar=i2s(idec);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[0]);
    nobj += fwrite(&rar,sizeof(rar),1,ofls[0]);
    nobj += fwrite(&decr,sizeof(decr),1,ofls[0]);
    nobj += fwrite(&epoch,sizeof(epoch),1,ofls[0]);
    nobj += fwrite(&sflux,sizeof(sflux),1,ofls[0]);
    nobj += fwrite(&size,sizeof(size),1,ofls[0]);
}

/************************************************************************/
void blWrite_c(const int blhid, const int inhid, const int isb, const int ipol, const float pa, const int iaq, const int ibq, const int icq, const int ioq, const int irec, const int iifc, const float u, const float v, const float w, const float prbl, const float angres, const float vis, const float coh, const float sigcoh, const float csnr, const float vflux, const float cnoise, const double avedhrs, const float ampav, const float phaave, const float tpvar, const int blsid, const int itel1, const int itel2, const int iblcd, const float ble, const float bln, const float blu, const int soid){
/**blWrite -- write the baseline header for a MIR data set              */
/*&dnf                                                                  */
/*:low-level-i/o                                                        */
/*+ FORTRAN call sequence

       subroutine blWrite(blhid,inhid,isb,ipol,pa,iaq,ibq,icq,ioq,irec,
                         iifc,u,v,w,prbl,angres,vis,coh,sigcoh,csnr,
                         vflux,cnoise,avedhrs,ampav,phaave,tpvar,blsid,
                         itel1,itel2,iblcd,ble,bln,blu,soid)

       integer blhid,inhid,isb,ipol,iaq,ibq,icq,ioq,irec,iifc,blsid,
               itel1,tel2,iblcd
       real pa,u,v,w,prbl,angres,vis,coh,sigcoh,csnr,vflux,cnoise,ampav,
            phaave,tpvar,ble,bln,blu,soid
       double precision avedhrs

  This writes an entry into a MIR baseline header

  Input:
    blhid     baseline header id (globally unique)
    inhid     integration header id (globally unique, same as that for integration header)
    isb       sideband code
    ipol      polarization code
    pa        polarization angle
    iaq       amplitude quality code
    ibq       baseline quality code
    icq       coherence quality code
    ioq       offset quality code
    irec      receiver code
    iifc      if channel code
    u         u in meters
    v         v in meters
    w         w in meters
    prbl      projected baseline in meters
    angres    angular resolution of baseline in arcsec
    vis       fractional visibility of planet (if applicable)
    coh       coherence
    sigcoh    uncertainty of coherence
    csnr      continuum SNR of planet (if applicable)
    vflux     visibility flux of planet (if applicable)
    cnoise    continuum noise
    avedhrs   offset from reference time
    ampave    average wide band amplitude
    phaave    average wide band phase
    tpvar     total power stability
    blsid     baseline id (a1*256 + a2)
    itel1     antenna 1 code
    itel2     antenna 2 code
    iblcd     baseline code
    ble       baseline east vector (meters)
    bln       baseline north vector (meters)
    blu       baseline up vector (meters)
    soid      baseline solution id                                      */
/*--                                                                    */
/*----------------------------------------------------------------------*/
    int nobj=0;
    short svar;

    nobj += fwrite(&blhid,sizeof(blhid),1,ofls[1]);
    if(nobj == 0){
	bug_c('f',"Unable to write to bl_read.");	
    }
    nobj += fwrite(&inhid,sizeof(inhid),1,ofls[1]);
    svar=i2s(isb);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[1]);
    svar=i2s(ipol);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[1]);
    nobj += fwrite(&pa,sizeof(pa),1,ofls[1]);
    svar=i2s(iaq);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[1]);
    svar=i2s(ibq);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[1]);
    svar=i2s(icq);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[1]);
    svar=i2s(ioq);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[1]);
    svar=i2s(irec);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[1]);
    svar=i2s(iifc);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[1]);
    nobj += fwrite(&u,sizeof(u),1,ofls[1]);
    nobj += fwrite(&v,sizeof(v),1,ofls[1]);
    nobj += fwrite(&w,sizeof(w),1,ofls[1]);
    nobj += fwrite(&prbl,sizeof(prbl),1,ofls[1]);
    nobj += fwrite(&angres,sizeof(angres),1,ofls[1]);
    nobj += fwrite(&vis,sizeof(vis),1,ofls[1]);
    nobj += fwrite(&coh,sizeof(coh),1,ofls[1]);
    nobj += fwrite(&sigcoh,sizeof(sigcoh),1,ofls[1]);
    nobj += fwrite(&csnr,sizeof(csnr),1,ofls[1]);
    nobj += fwrite(&vflux,sizeof(vflux),1,ofls[1]);
    nobj += fwrite(&cnoise,sizeof(cnoise),1,ofls[1]);
    nobj += fwrite(&avedhrs,sizeof(avedhrs),1,ofls[1]);
    nobj += fwrite(&ampav,sizeof(ampav),1,ofls[1]);
    nobj += fwrite(&phaave,sizeof(phaave),1,ofls[1]);
    nobj += fwrite(&tpvar,sizeof(tpvar),1,ofls[1]);
    nobj += fwrite(&blsid,sizeof(blsid),1,ofls[1]);
    svar=i2s(itel1);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[1]);
    svar=i2s(itel2);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[1]);
    svar=i2s(iblcd);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[1]);
    nobj += fwrite(&ble,sizeof(ble),1,ofls[1]);
    nobj += fwrite(&bln,sizeof(bln),1,ofls[1]);
    nobj += fwrite(&blu,sizeof(blu),1,ofls[1]);
    nobj += fwrite(&soid,sizeof(soid),1,ofls[1]);   
}

/************************************************************************/
void spWrite_c(const int sphid, const int blhid, const int inhid, const int igq, const int ipq, const int iband, const int ipstate, const float tau0, const double vel, const float vres, const int ivtype, const double fsky, const float fres, const float tssb, const float integ, const float wt, const int itaper, const float snoise, const int nch, const int nrec, const int dataoff, const int linid, const int itrans, const double rfreq, const int pasid, const int gaiidamp, const int gaiidpha, const int flcid, const int atmid){
/**spWrite -- write the spectral window header for a MIR data set       */
/*&dnf                                                                  */
/*:low-level-i/o                                                        */
/*+ FORTRAN call sequence

       subroutine spWrite(sphid,blhid,inhid,igq,ipq,iband,ipstate,tau0,
                          vel,vres,ivtype,fsky,fres,tssb,integ,wt,
                          itaper,snoise,nch,nrec,dataoff,linid,itrans,
                          rfreq,pasid,gaiidamp,gaiidpha,flcid,atmid)

       integer sphid,blhid,inhid,igq,ipq,iband,ipstate,ivtype,itaper,
               nch,nrec,dataoff,linid,itrans,pasid,gaiidamp,gaiidpha,
               flcid,atmid
       real tau0,vres,fres,tssb,integ,wt,snoise
       double precision vel,fsky,rfreq

  This writes an entry into a MIR spectral window header

  Input:
    sphid      spectral header id (globally unique)
    blhid      baseline header id (globally unique, same as that in baseline header)
    inhid      integration header id (globally unique, same as that in integration header)
    igq        gain calibrator code
    ipq        passband calibrator code
    iband      spectral band code
    ipstate    polarization state code
    tau0       opacity
    vel        velocity of center channel
    vres       velocity resolution of channel (km/s)
    ivtype     velocity axis type
    fsky       sky frequency of window (centered in window, GHz)
    fres       frequency resolution (MHz)
    tssb       baseline based single sideband system temperature (K)
    integ      integration time
    wt         weight of data - based on integration time and system temp
    itaper     taper code
    snoise     single channel noise
    nch        number of spectral channels in window
    nrec       number of records
    dataoff    data offset in bytes
    linid      spectral line id
    itrans     transition code
    rfreq      rest frequency (GHz)
    pasid      passband function number
    gaiidamp   amp gain function number
    gaiidpha   phase gain function number
    flcid      flux cal id
    atmid      atm ospheric coherence id                                */
/*--                                                                    */
/*----------------------------------------------------------------------*/
    int nobj=0;
    short svar;

    nobj += fwrite(&sphid,sizeof(sphid),1,ofls[2]);
    if(nobj == 0){
	bug_c('f',"Unable to write to sp_read.");	
    }
    nobj += fwrite(&blhid,sizeof(blhid),1,ofls[2]);
    nobj += fwrite(&inhid,sizeof(inhid),1,ofls[2]);
    svar=i2s(igq);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
    svar=i2s(ipq);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
    svar=i2s(iband);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
    svar=i2s(ipstate);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
    nobj += fwrite(&tau0,sizeof(tau0),1,ofls[2]);
    nobj += fwrite(&vel,sizeof(vel),1,ofls[2]);
    nobj += fwrite(&vres,sizeof(vres),1,ofls[2]);
    svar=i2s(ivtype);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
    nobj += fwrite(&fsky,sizeof(fsky),1,ofls[2]);
    nobj += fwrite(&fres,sizeof(fres),1,ofls[2]);
    nobj += fwrite(&tssb,sizeof(tssb),1,ofls[2]);
    nobj += fwrite(&integ,sizeof(integ),1,ofls[2]);
    nobj += fwrite(&wt,sizeof(wt),1,ofls[2]);
    svar=i2s(itaper);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
    nobj += fwrite(&snoise,sizeof(snoise),1,ofls[2]);
    svar=i2s(nch);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
    svar=i2s(nrec);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
    nobj += fwrite(&dataoff,sizeof(dataoff),1,ofls[2]);
    nobj += fwrite(&linid,sizeof(linid),1,ofls[2]);
    svar=i2s(itrans);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
    nobj += fwrite(&rfreq,sizeof(rfreq),1,ofls[2]);
    svar=i2s(pasid);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
    svar=i2s(gaiidamp);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
    svar=i2s(gaiidpha);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
    svar=i2s(flcid);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
    svar=i2s(atmid);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[2]);
}

/************************************************************************/
void codeWrite_c(const char *v_name, const int icode, const char *code, const int ncode){
/**codeWrite -- write the code index for a MIR data set                 */
/*&dnf                                                                  */
/*:low-level-i/o                                                        */
/*+ FORTRAN call sequence

       subroutine codeWrite(v_name,icode,code,ncode)
       
       character v_name*(*),code*(*)
       integer icode,ncode

  This writes an entry into a MIR code index

  Input:
    v_name   label of code
    icode    index # for code
    code     code value
    ncode    length of code value (# characters)                        */
/*--                                                                    */
/*----------------------------------------------------------------------*/
    int nobj=0;
    int i;
    short svar;
    char temp[27];
    for(i=0;i<27;i++){          //pad the character array
	temp[i]=' ';
    }
    for(i=0;i<strlen(v_name);i++){
	temp[i]=v_name[i];
    }
    nobj += fwrite(temp,12,1,ofls[3]);
    if(nobj == 0){
	bug_c('f',"Unable to write to codes_read.");	
    }
    svar=i2s(icode);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[3]);
    for(i=0;i<27;i++){          //pad the character array
	temp[i]=' ';
    }
    for(i=0;i<strlen(code);i++){
	temp[i]=code[i];
    }
    nobj += fwrite(temp,26,1,ofls[3]);
    svar=i2s(ncode);
    nobj += fwrite(&svar,sizeof(svar),1,ofls[3]);
}

/************************************************************************/
void visWrite_c(const float *re, const float *im, const int numvis, const int startvis, int *nbytes){
/**visWrite -- write the visibility file for a MIR data set             */
/*&dnf                                                                  */
/*:low-level-i/o                                                        */
/*+ FORTRAN call sequence

       subroutine visWrite(re,im,numvis,startvis,nbytes)

       real re(*),im(*)
       integer iscale,time,numvis

  This writes an entry into a MIR visibility file
  Input:
    re       real part of visibility converted to an integer
    im       imaginary part of visibility converted to an integer
    numvis   number of visibilities to write
    startvis starting number of the visibilities to write
    nbytes   number of bytes that the data contain                      */
/*--                                                                    */
/*----------------------------------------------------------------------*/
    int nobj=0;
    int bytes;
    int i;

    *nbytes=0;
    for(i=startvis;i<(numvis+startvis);i++){
	nobj += fwrite(&re[i],(bytes=sizeof(re[i])),1,ofls[4]);
	*nbytes+=bytes;
	nobj += fwrite(&im[i],(bytes=sizeof(im[i])),1,ofls[4]);
	*nbytes+=bytes;
    }
}


//function to convert integer to a short
//this is needed as it is difficult to pass an integer*2 from miriad
//  to the C subroutines. The conversion is easier in C
short i2s(const int num){
    return (short)num;
}
