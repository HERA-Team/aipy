#ifndef _DSP_H_
#define _DSP_H_

#include "grid.h"
#include <Python.h>
#include "numpy/arrayobject.h"

#define QUOTE(s) # s
#define CHK_ARRAY_TYPE(a,type) \
    if (PyArray_TYPE(a) != type) { \
        PyErr_Format(PyExc_ValueError, "type(%s) != %s", \
        QUOTE(a), QUOTE(type)); \
        return NULL; }
#define CHK_ARRAY_DIM(a,i,d) \
    if (PyArray_DIM(a,i) != d) { \
        PyErr_Format(PyExc_ValueError, "dim(%s) != %s", \
        QUOTE(a), QUOTE(d)); \
        return NULL; }
#define RANK(a) a->nd
#define CHK_ARRAY_RANK(a,r) \
    if (RANK(a) != r) { \
        PyErr_Format(PyExc_ValueError, "rank(%s) != %s", \
        QUOTE(a), QUOTE(r)); \
        return NULL; }

#endif
