#include "miriad.h"
#include <stdio.h>
#include "Python.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "io.h"

#define check(iostat) if(iostat)bugno_c('f',iostat)
#define MAXSIZE 1024
#define MAXLINE 80

extern void uvread_c_wrap(int tno, double *preamble, int n0, float *data, int n1, int *flags, int n2, int *nread);

extern void uvgetvr_c_wrap(int tno, int type, const char *var, char *data, int *maxdata, int nread, int bytesize);

extern void uvwrite_c_wrap(int tno, double *preamble, int n0, float *data, int n1, int *flags, int n2);

extern int miriad_error_exists();

extern char *miriad_error_string();

extern void wrhdc_c_wrap(int thandle, const char *keyword, const float *value, int n);

extern size_t rdhdc_c_wrap(int thandle, const char *keyword, float *value, int n);

extern size_t hsize_c_wrap(int thandle, const char *keyword);

extern void write_freqs(int thandle, int nspect, int nschan, double sfreq, double sdf);
