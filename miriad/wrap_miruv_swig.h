#include "miriad.h"
#include <stdio.h>
#include "Python.h"

extern void uvread_c_wrap(int tno, double *preamble, int n0, float *data, int n1, int *flags, int n2, int *nread);

extern void uvgetvr_c_wrap(int tno, int type, const char *var, char *data, int *maxdata, int nread, int bytesize);

extern void uvwrite_c_wrap(int tno, double *preamble, int n0, float *data, int n1, int *flags, int n2);

extern int miriad_error_exists();

extern char *miriad_error_string();
