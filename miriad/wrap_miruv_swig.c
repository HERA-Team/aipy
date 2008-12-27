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

