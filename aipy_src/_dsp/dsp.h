#ifndef _DSP_H_
#define _DSP_H_

#include <Python.h>
#include "grid.h"
#include "numpy/arrayobject.h"

// Python3 compatibility
#if PY_MAJOR_VERSION >= 3
	#define PyCapsule_Type PyCObject_Type
	#define PyInt_AsLong PyLong_AsLong
	#define PyInt_FromLong PyLong_FromLong
	#define PyString_FromString PyUnicode_FromString
char* PyString_AsString(PyObject *ob) {
	PyObject *enc;
	char *cstr;
	enc = PyUnicode_AsEncodedString(ob, "utf-8", "Error");
	if( enc == NULL ) {
		PyErr_Format(PyExc_ValueError, "Cannot encode string");
		return NULL;
	}
	cstr = PyBytes_AsString(enc);
	Py_XDECREF(enc);
	return cstr;
}
#endif

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
#define RANK(a) PyArray_NDIM(a)
#define CHK_ARRAY_RANK(a,r) \
    if (RANK(a) != r) { \
        PyErr_Format(PyExc_ValueError, "rank(%s) != %s", \
        QUOTE(a), QUOTE(r)); \
        return NULL; }

#endif
