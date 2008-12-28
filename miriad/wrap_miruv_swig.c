/*  Author: Aaron Parsons
    Date: 10/11/06
    Revisions:
        arp  10/25/06   Renaming files, functions
        arp  01/21/07   Made preamble an array

    Some miriad routines need a little extra wrapping to make nice interfaces
    to python through SWIG.
*/

#include "wrap_miruv_swig.h"

int miriad_error=0;
char miriad_error_str[256];

void uvread_c_wrap(int tno, double *preamble, int n0, float *data, int n1, 
                   int *flags, int n2, int *nread) {
    /* Wrap uvread_c to have a length argument for both data and flags, and
    then check to see that the length of data is twice the length of flags.
    Raises python exception if untrue. */
    if (n1 != 2*n2) {
        miriad_error = 1;
        sprintf(miriad_error_str, 
            "data (size=%d) must be exactly twice as long as flags (size=%d)", 
            n1, n2);
        return;
    }
    uvread_c(tno, preamble, data, flags, n2, nread);
    return;
}

void uvgetvr_c_wrap(int tno, int type, const char *var, char *data,
                    int *maxdata, int nread, int bytesize) {
    /* Wrap uvgetvr_c to have a length argument for data.  Sets this length
    argument to the amount of data expected to get back from uvgetvr_c. */
    *maxdata = nread * bytesize;
    return uvgetvr_c(tno, type, var, data, nread);
}

void uvwrite_c_wrap(int tno, double *preamble, int n0, float *data, int n1,
                    int *flags, int n2) {
    /* Wrap uvwrite_c to have length argument for both data and flags, and
    then check to see that the length of data is twice the length of flags.
    Raises python exception if untrue. */
    if (n1 != 2*n2) {
        miriad_error = 1;
        sprintf(miriad_error_str,
            "data (size=%d) must be exactly twice as long as flags (size=%d)", 
            n1, n2);
        return;
    }
    return uvwrite_c(tno, preamble, data, flags, n2);
}

int miriad_error_exists() {
    /* Rudimentary error-handling system. */
    return miriad_error;
}

char *miriad_error_string() {
    /* Rudimentary error-handling system. */
    return miriad_error_str;
}

void wrhdc_c_wrap(int thandle, const char *keyword, const float *value, int n) {
    int item;
    int iostat,offset;

    if (n % 2 != 0) {
        miriad_error = 1;
        sprintf(miriad_error_str, "value (size=%d) must have even length", n);
        return;
    }

    haccess_c(thandle,&item,keyword,"write",&iostat);     check(iostat);
    hwriteb_c(item,cmplx_item,0,ITEM_HDR_SIZE,&iostat);   check(iostat);
    offset = mroundup(ITEM_HDR_SIZE,H_CMPLX_SIZE);
    hwritec_c(item,value,offset,H_CMPLX_SIZE*n/2,&iostat);  check(iostat);
    hdaccess_c(item,&iostat);                             check(iostat);
}

size_t rdhdc_c_wrap(int thandle, const char *keyword, float *value, int n) {
  int item;
  char s[ITEM_HDR_SIZE];
  int iostat,offset,length;

    if (n % 2 != 0) {
        miriad_error = 1;
        sprintf(miriad_error_str, "value (size=%d) must have even length", n);
        return 0;
    }

/* Firstly assume the variable is missing. Try to get it. If successful
   read it. */

  haccess_c(thandle,&item,keyword,"read",&iostat);  if(iostat)return 0;
  offset = mroundup(ITEM_HDR_SIZE,H_CMPLX_SIZE);
  length = hsize_c(item) - offset;
  if(length > H_CMPLX_SIZE * n / 2){
    length = H_CMPLX_SIZE * n / 2;
  }
  hreadb_c(item,s,0,ITEM_HDR_SIZE,&iostat);       check(iostat);
  iostat = 0;
  if(!memcmp(s,cmplx_item, ITEM_HDR_SIZE)){
    hreadc_c(item,value,offset,length,&iostat);
  }
  check(iostat);
  hdaccess_c(item,&iostat);             check(iostat);
  return length / H_CMPLX_SIZE * 2;
}

size_t hsize_c_wrap(int thandle, const char *keyword) {
  int item,offset,length,iostat;
  haccess_c(thandle,&item,keyword,"read",&iostat);  if(iostat)return 0;
  offset = mroundup(ITEM_HDR_SIZE,H_CMPLX_SIZE);
  length = hsize_c(item) - offset;
  hdaccess_c(item,&iostat);
  return length / H_CMPLX_SIZE * 2;
}

void write_freqs(int thandle, int nspect, int nschan, double sfreq, 
                 double sdf) {
  int item,offset,iostat,i;
  haccess_c(thandle,&item,"freqs","write",&iostat); if(iostat)return;
  offset = 0;
  /* I think I can put anything here, so i'm putting nspect */
  hwritei_c(item,&nspect,offset,4,&iostat); if(iostat)return;
  offset += 8;
  for(i=0; i < nspect; i++) {
      hwritei_c(item,&nschan,offset,4,&iostat); if(iostat)return;
      offset += 8;
      hwrited_c(item,&sfreq,offset,8,&iostat); if(iostat)return;
      offset += 8;
      hwrited_c(item,&sdf,offset,8,&iostat); if(iostat)return;
      offset += 8;
  }
  hdaccess_c(item,&iostat);
}

void read_freqs(int thandle, int *nspect, int *nschan, double *sfreq,
                double *sdf) {
  int item,offset,iostat,i;
  haccess_c(thandle,&item,"freqs","read",&iostat); if(iostat)return;
  offset = 0;
  hreadi_c(item,nspect,offset,4,&iostat); if(iostat)return;
  offset += 8;
  /* for(i=0; i < *nspect; i++) { */
  /* Only supporting that there is 1 set of nschan,sfreq,sdf for all ants */
  for(i=0; i < 1; i++) {
      hreadi_c(item,nschan,offset,4,&iostat); if(iostat)return;
      offset += 8;
      hreadd_c(item,sfreq,offset,8,&iostat); if(iostat)return;
      offset += 8;
      hreadd_c(item,sdf,offset,8,&iostat); if(iostat)return;
      offset += 8;
  }
  hdaccess_c(item,&iostat);
}
