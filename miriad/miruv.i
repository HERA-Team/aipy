%module(package="miruv", docstring="This is a SWIG interface for binding to Miriad I/O libraries") miruv

/*
Author: Aaron Parsons, borrowing from early work by Garrelt Mallema
Date: 10/05/06
Revisions:
    11/02/06 arp Added support for arbitrary length preambles.
    01/21/07 arp Made preamble an array b/c of memory allocation issues
    01/28/07 arp Added uvselect, uvtrack, uvcopyvr
    06/10/07 arp Addio hio_c, but had to change off_t to size_t in hio.c
                 miriad.h, and miruv.i

ToDo:
    Find a way to have bug_c raise python exceptions.

*/

/*
             _             _                
 _   ___   _(_) ___     __| | ___   ___ ___ 
| | | \ \ / / |/ _ \   / _` |/ _ \ / __/ __|
| |_| |\ V /| | (_) | | (_| | (_) | (__\__ \
 \__,_| \_/ |_|\___/   \__,_|\___/ \___|___/
*/

%define UVOPEN_DOCSTRING
"Input:
    name    Name of the directory tree containg the u-v data.
    status  Either 'old', 'new' or 'append'. Old files can be read,
            whereas new and append files can only be written. Append
            files must be pre-existing uv data-sets.
Output:
    tno     Handle of the uv data set. "
%enddef

%define UVCLOSE_DOCSTRING
"This closes a uv data file.
Input:
    tno     Handle of the uv data set. "
%enddef

%define UVREWIND_DOCSTRING
"Rewind a uv file, readying it to be read from the begining again.
Input:
    tno     The uv data file handle.  "
%enddef

%define UVFLUSH_DOCSTRING
"This closes a uv data file.
Input:
    tno     Make sure anything buffered up is flushed to disk. The
            disk file should be readable (up to data written here)
            even if the caller or computer crashes.  "
%enddef

%define UVNEXT_DOCSTRING
"Skip to the next uv data record. On write, this causes an end-of-record
mark to be written. On read, this causes data to be read up until the
next end-of-record mark.
Input:
    tno     The uv data file handle. "
%enddef

%define UVCOPYVR_DOCSTRING
"This copies those variables, in the input uv data set, which have
changed and which are marked as 'copy' ('u' flag of a call to uvtrack).
Inputs:
    tin     File handle of the input uv data set.
    tout    File handle of the output uv data set."
%enddef

%define UVUPDATE_DOCSTRING
"This checks whether any ``important variables'' has been updated in the
last call to uvread or uvscan. Important variables are those flagged
with the 'u' flag in a call to uvtrack.
Input:
    tno     File handle of the uv file to check.  "
%enddef

%define UVRDVR_DOCSTRING
"These routines get the first value of a variable. If the variable is
missing,the default value is returned.
Input:
    tno     The handle of the uv data file.
    varname The name of the variable to return.
    default The default value.
Output:
    data    The returned value of the variable. "
%enddef

%define UVGETVR_DOCSTRING
"These routines return the current value of a uv variable. N gives the size
of elements in the return array. This MUST match with the actual number
of elements in the variable. An exception is for the character string
routine, where the size of the 'adata' string must be strictly greater
than (not equal to!) the size of the string.
Input:
    tno     The handle of the uv data file.
    varname The name of the variable to return.
    n       The number of elements to return. This must agree with
            the size of the variable!
Output:
    data    The returned values of the variable. "
%enddef

%define UVPROBVR_DOCSTRING
"This checks whether a particular variable exists. If it does, this
passes back the type and (current) length of the variable, and whether
it was updated on the last call to uvread or uvscan.
Input:
    tno     The handle of the input uv data file.
    varname The name of the variable to check.
Output:
    type    The type of the variable. If the variable does not
            exist, this is blank. Otherwise it is one of 'a', 'r',
            'i', 'd' or 'c'.
    length  The number of elements in the uv variable. If this is not
            known (which is true if the variable has never been read)
            then this will be zero.
    update  This will be set .true. if this variable was updated
            on the last call to uvread or uvscan. "
%enddef

%define UVPUTVR_DOCSTRING
"These routines write new values for a uv variable. N gives the number
of elements to write.
Input:
    tno     The handle of the uv data file.
    varname The name of the variable to write.
    n       The number of elements to write.
    data    The values of the variable.    "
%enddef

%define UVSELECT_DOCSTRING
"This specifies the portion of the data to be selected by calls to
uvread. Normally data that are not selected, are not returned.
Exceptions are the 'window' and 'amplitude' objects, which cause the
corresponding visibilities to be flagged as bad, but returned anyway.
Input:
    tno     Handle of the uv data file.
    object  This can be one of 'time','antennae','visibility',
            'uvrange','pointing','amplitude','window','or','dra',
            'ddec','uvnrange','increment','ra','dec','and', 'clear',
            'on','polarization','shadow','auto','dazim','delev'
    p1,p2   Generally this is the range of values to select. For
            'antennae', this is the two antennae pair to select.
            For 'antennae', a zero indicates 'all antennae'.
            For 'shadow', a zero indicates use 'antdiam' variable.
            For 'on','window','polarization','increment','shadow' only
            p1 is used.
            For 'and','or','clear','auto' p1 and p2 are ignored.
    flag    If true, the data is selected. If false, the data is
            discarded. Ignored for 'and','or','clear'. "
%enddef

%define UVTRACK_DOCSTRING
"UVTRACK allows the programmer to set switches and flags associated with
a particular uv variable, to allow extra processing steps of that
variable.
Input:
  tno         The handle of the input uv file.
  varname     The name of the variable of interest.
  switches    This is a character string, each character of which
              causes a particular flag or switch to be turned on for
              this particular variable. Valid values are:
               'u'  Remember if this variable gets updated, and  when
                    it gets updated, uvupdate returns .true. the next
                    time it is called.
               'c'  Remember if this variable gets updated, and when
                    it gets updated, cause it to be copied during the
                    next call to uvcopyvr."
%enddef

%define UVWRITE_DOCSTRING
"Write a visibility record to the data file.  Please note uvwrite()
closes the record. Any wideband data should have been written with
uvwwrite() before this call.
Input:
    tno     Handle of the uv data set.
    n       Number of channels to be written.
    preamble    A double array of 4 elements giving u,v, time and
            baseline number (in that order).
    data    A complex array of n elements containing
            the correlation data.
    flags   Logical array of n elements. A true value for
            a channel indicates good data. "
%enddef

%define UVREAD_DOCSTRING
"This reads a single visibility from the data file. This starts by scanning
the uv data stream until a correlation record is found. The correlations
are then converted to complex numbers if necessary, and returned to the
caller. Uvread also performs some massaging (see uvset) and selection
(see uvselect) steps.
Input:
    tno     Handle of the uv data set.
    n       Max number of channels that can be read.
Output:
    preamble    A double array of elements giving things such as
            u,v, time and baseline number. Setable using uvset.
    data    A real array of at least n complex elements (or 2n real
            elements). This returns the correlation data.
    flags   Logical array of at least n elements. A true value for
            a channel indicates good data.
    nread   Number of correlations returned. On end-of-file, zero
            is returned. "
%enddef

%define UVSET_DOCSTRING
"Set up the way uvread behaves. This determines whether uvread returns
correlation channels, wide-band channels or velocity channels. This also
sets up whether u and v are returned in wavelengths or nanosec, and
what planet processing is performed.
Input:
    tno     Handle of the uv data set.
    object  Name of the object that we are setting the type of.
    type    The type of data that the user wants returned.
    n       Some integer parameter.
    p1,p2,p3    Some real parameters.   "
%enddef

/*
 _     _             _                
| |__ (_) ___     __| | ___   ___ ___ 
| '_ \| |/ _ \   / _` |/ _ \ / __/ __|
| | | | | (_) | | (_| | (_) | (__\__ \
|_| |_|_|\___/   \__,_|\___/ \___|___/
*/


%define HACCESS_DOCSTRING
"Miriad data sets consist of a collection of items. Before an item within
a data set can be read/written, etc, it must be 'opened' with the haccess
routine.
Input:
    tno     The handle of the data set.
    keyword The name of the item.
    status  This can be 'read', 'write', 'append' or 'scratch'.
Output:
    ihandle The handle of the opened item. Note that item handles are
            quite distinct from data-set handles.
    iostat  I/O status indicator. 0 indicates success. Other values
            are standard system error numbers.    "
%enddef

%define HDACCESS_DOCSTRING
"This releases an item. It flushes buffers and waits for i/o to complete.
For small items that are entirely in memory, these are saved until
the whole tree is closed before they are written out.
Input:
    itno    The handle of the item to close up.
Output:
    iostat  I/O status indicator. 0 indicates success. Other values
            are standard system error numbers.  "
%enddef

%define HIO_DOCSTRING
"Read or write items of a Miriad data set.
Input:
    ihandle The handle of the item to perform I/O on.
    dowrite Perform a write (as opposed to a read).
    type    Miriad number code for data type.
    offset  The byte offset into the item where I/O is to be
            performed.
    length  The number of bytes to be read.
Output:
    iostat  I/O status indicator. 0 indicates success. -1 indicates
            end-of-file. Other values are standard system
            error numbers. 
Either:
    buf     In write mode, the data to be written.  In read mode, returns
            the data read."
%enddef

%define HREADA_DOCSTRING
""
%enddef

%define HWRITEA_DOCSTRING
""
%enddef

%define RDHD_DOCSTRING
""
%enddef

%define WRHD_DOCSTRING
""
%enddef

/*
           _                                        _       _      
 _ __ ___ (_)_ __ _   ___   __  _ __ ___   ___   __| |_   _| | ___ 
| '_ ` _ \| | '__| | | \ \ / / | '_ ` _ \ / _ \ / _` | | | | |/ _ \
| | | | | | | |  | |_| |\ V /  | | | | | | (_) | (_| | |_| | |  __/
|_| |_| |_|_|_|   \__,_| \_/   |_| |_| |_|\___/ \__,_|\__,_|_|\___|
*/

%{
#define SWIG_FILE_WITH_INIT
#define MAX_PREAMBLE 9
#include "io.h"
#include "miriad.h"
#include "wrap_miruv_swig.h"
%}

%include "typemaps.i"
%include "cstring.i"
%include "numpy.i"

%init %{
  import_array();
%}

%feature("autodoc", "1");

/* Implement a rudimentary error-handling scheme to raise python exceptions.
This will need work to intercept calls to bug_c by Miriad, which write to
stderr and sometimes exit the program. */
%exception {
   $action
   if (miriad_error_exists()) {
      PyErr_SetString(PyExc_RuntimeError, miriad_error_string());
      return NULL;
   }
}

/*
 _   ___   _(_) ___   ___ 
| | | \ \ / / |/ _ \ / __|
| |_| |\ V /| | (_) | (__ 
 \__,_| \_/ |_|\___(_)___|
*/

/* uvopen_c */
%feature("docstring", UVOPEN_DOCSTRING);
%apply int *OUTPUT {int *tno};
extern void uvopen_c (int *tno, const char *name, const char *status);

/* uvrewind_c */
%feature("docstring", UVREWIND_DOCSTRING);
extern void uvrewind_c (int tno);

/* uvclose_c */
%feature("docstring", UVCLOSE_DOCSTRING);
extern void uvclose_c (int tno);

/* uvflush_c */
%feature("docstring", UVFLUSH_DOCSTRING);
extern void uvflush_c (int tno);

/* uvnext_c */
%feature("docstring", UVNEXT_DOCSTRING);
extern void uvnext_c (int tno);

/* uvupdate_c */
%feature("docstring", UVUPDATE_DOCSTRING);
extern int uvupdate_c (int tno);

/* uvset_c */
%feature("docstring", UVSET_DOCSTRING);
extern void uvset_c (int tno, const char *object, const char *type, 
                     int n, double p1, double p2, double p3);

/* uvread_c_wrap */
/* from wrap_miruv_swig.c, which is a thin wrapper over uvread_c */
%feature("docstring", UVREAD_DOCSTRING);
%apply  (double *INPLACE_ARRAY1, int DIM1) {(double *preamble, int n0)};
%apply  (float *INPLACE_ARRAY1, int DIM1) {(float *data, int n1)};
%apply  (int *INPLACE_ARRAY1, int DIM1) {(int *flags, int n2)};
%apply int *OUTPUT {int *nread};
extern void uvread_c_wrap (int tno, double *preamble, int n0, float *data, int n1, int *flags, int n2, int *nread);

/* uvwrite_c_wrap */
/* from wrap_miruv_swig.c, which is a thin wrapper over uvwrite_c */
%feature("docstring", UVWRITE_DOCSTRING);
%apply (double *IN_ARRAY1, int DIM1) {(double *preamble, int n0)};
%apply (float *IN_ARRAY1, int DIM1) {(float *data, int n1)};
%apply (int *IN_ARRAY1, int DIM1) {(int *flags, int n2)};
extern void uvwrite_c_wrap (int tno, double *preamble, int n0, float *data, int n1, int *flags, int n2);

/* uvprobvr_c */
%feature("docstring", UVPROBVR_DOCSTRING);
%apply int *OUTPUT {int *length, int *updated};
%cstring_bounded_mutable(char *type, 64);
extern void uvprobvr_c (int tno, const char *var, char *type, int *length, int *updated);

/* uvgetvr_c */
%feature("docstring", UVGETVR_DOCSTRING);
%cstring_output_withsize(char *data, int *maxdata);
extern void uvgetvr_c_wrap (int tno, int type, const char *var, char *data, int *maxdata, int nread, int bytesize);

/* uvputvr_c */
%feature("docstring", UVPUTVR_DOCSTRING);
extern void uvputvr_c (int tno, int type, const char *var, const char *data, int n);

/* uvselect_c */
%feature("docstring", UVSELECT_DOCSTRING);
extern void uvselect_c (int tno, const char *object, double p1, double p2, int datasel);

/* uvtrack_c */
%feature("docstring", UVTRACK_DOCSTRING);
extern void uvtrack_c (int tno, const char *name, const char *switches);

/* uvcopyvr_c */
%feature("docstring", UVCOPYVR_DOCSTRING);
extern void uvcopyvr_c (int tin, int tout);

/*
 _                    _ _            
| |__   ___  __ _  __| (_) ___   ___ 
| '_ \ / _ \/ _` |/ _` | |/ _ \ / __|
| | | |  __/ (_| | (_| | | (_) | (__ 
|_| |_|\___|\__,_|\__,_|_|\___(_)___|
*/

/* rdhdc_c */
%feature("docstring", RDHD_DOCSTRING);
%apply (float *INPLACE_ARRAY1, int DIM1) {(float *value, int n)};
extern size_t rdhdc_c_wrap (int thandle, const char *keyword, float *value, int n);

/* wrhdr_c */
%feature("docstring", WRHD_DOCSTRING);
extern void wrhdr_c (int tno, const char *keyword, double value);

/* wrhdd_c */
%feature("docstring", WRHD_DOCSTRING);
extern void wrhdd_c (int tno, const char *keyword, double value);

/* wrhdi_c */
%feature("docstring", WRHD_DOCSTRING);
extern void wrhdi_c (int tno, const char *keyword, int value);

/* wrhdc_c */
%feature("docstring", WRHD_DOCSTRING);
%apply (float *IN_ARRAY1, int DIM1) {(const float *value, int n)};
extern void wrhdc_c_wrap (int thandle, const char *keyword, const float *value, int n);

/* rdhdr_c */
%feature("docstring", RDHD_DOCSTRING);
%apply float *OUTPUT {float *value};
extern void rdhdr_c (int tno, const char *keyword, float *value, double defval);

/* rdhdd_c */
%feature("docstring", RDHD_DOCSTRING);
%apply double *OUTPUT {double *value};
extern void rdhdd_c (int tno, const char *keyword, double *value, double defval);

/* rdhdi_c */
%feature("docstring", RDHD_DOCSTRING);
%apply int *OUTPUT {int *value};
extern void rdhdi_c (int tno, const char *keyword, int *value, int defval);

/*%cstring_output_withsize(char *value, int *len);
extern void rdhda_c_wrap (int tno, const char *keyword, char *value, int *len, const char *defval);*/
/*extern void rdhdl_c (int tno, const char *keyword, int8 *value,int8 defval);*/
/*extern void wrhdl_c (int tno, const char *keyword, int8 value);*/
/*extern void wrhda_c (int tno, const char *keyword, const char *value);*/

/*
 _     _
| |__ (_) ___   ___
| '_ \| |/ _ \ / __|
| | | | | (_) | (__
|_| |_|_|\___(_)___|
*/

/* haccess_c */
%feature("docstring", HACCESS_DOCSTRING);
%apply int *OUTPUT {int *ihandle, int *iostat};
extern void haccess_c(int tno, int *ihandle, const char *keyword, const char *status, int *iostat);

/* hdaccess_c */
%feature("docstring", HDACCESS_DOCSTRING);
%apply int *OUTPUT {int *iostat};
extern void hdaccess_c(int ihandle, int *iostat);

/* hreada_c */
%feature("docstring", HREADA_DOCSTRING);
%cstring_output_maxsize(char *line, size_t length);
%apply int *OUTPUT {int *iostat};
extern void hreada_c(int ihandle, char *line, size_t length, int *iostat);

/* hwritea_c */
%feature("docstring", HWRITEA_DOCSTRING);
%apply int *OUTPUT {int *iostat};
extern void hwritea_c(int ihandle, const char *line, size_t length, int *iostat);

extern size_t hsize_c_wrap(int thandle, const char *keyword);

extern void write_freqs(int thandle, int nspect, int nschan, double sfreq, double sdf);

