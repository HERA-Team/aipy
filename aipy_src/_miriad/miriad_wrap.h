#ifndef _MIRIAD_WRAP_H_
#define _MIRIAD_WRAP_H_

#include <Python.h>
#include "miriad.h"
#include "numpy/arrayobject.h"
#include <string>
#include "hio.h"
#include "io.h"
#include "maxdimc.h"

// Python3 compatibility
#if PY_MAJOR_VERSION >= 3
	#define PyCapsule_Type PyCObject_Type
	#define PyInt_AsLong PyLong_AsLong
	#define PyInt_FromLong PyLong_FromLong
	#define PyInt_Check PyLong_Check
	#define PyString_Check PyUnicode_Check
	#define PyString_Size PyUnicode_GET_LENGTH
	#define PyString_FromString PyUnicode_FromString
	#define PyString_FromStringAndSize PyUnicode_FromStringAndSize
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

// Some miriad macros...
#define PREAMBLE_SIZE 5

//AAR: (use the other Miriad BL convention to accommodate more baselines)
#define GETI(bl) ( (int) (bl) > 65536 ? (((int) bl-65536)/2048 - 1) : (((int) bl >> 8) - 1) )
#define GETJ(bl) ( (int) (bl) > 65536 ? (((int) bl-65536)%2048 - 1 ) : (((int) bl & 255) - 1) )
#define MKBL(i,j) ( i+1 < 256 && j+1 < 256 ? ((float) (((i+1)<<8) | (j+1))) : ((float) (((i+1)*2048) + (j+1+65536))))

#define CHK_IO(i) \
    if (i != 0) { \
        PyErr_Format(PyExc_IOError, "IO failed"); \
        return NULL; }


// Some numpy macros...
#define QUOTE(a) # a
#define IND1(a,i,type) *((type *)((char *)PyArray_DATA(a) + i*PyArray_STRIDES(a)[0]))
#define IND2(a,i,j,type) *((type *)((char *)PyArray_DATA(a)+i*PyArray_STRIDES(a)[0]+j*PyArray_STRIDES(a)[1]))
#define TYPE(a) PyArray_DESCR(a)->type_num
#define CHK_ARRAY_TYPE(a,type) \
    if (TYPE(a) != type) { \
        PyErr_Format(PyExc_ValueError, "type(%s) != %s", \
        QUOTE(a), QUOTE(type)); \
        return NULL; }
#define DIM(a,i) PyArray_DIM(a, i)
#define RANK(a) PyArray_NDIM(a)
#define CHK_ARRAY_RANK(a,r) \
    if (RANK(a) != r) { \
        PyErr_Format(PyExc_ValueError, "rank(%s) != %s", \
        QUOTE(a), QUOTE(r)); \
        return NULL; }
#define CHK_NULL(a) \
    if (a == NULL) { \
        PyErr_Format(PyExc_MemoryError, "Failed to allocate %s", QUOTE(a)); \
        return NULL; }

// Some python macros...
#define CHK_STRING(o) \
    if (!PyString_Check(o)) { \
        PyErr_Format(PyExc_ValueError, "expected a string"); \
        return NULL; }
#define CHK_INT(o) \
    if (!PyInt_Check(o)) { \
        PyErr_Format(PyExc_ValueError, "expected an int"); \
        return NULL; }
#define CHK_LONG(o) \
    if (!PyLong_Check(o)) { \
        PyErr_Format(PyExc_ValueError, "expected a long"); \
        return NULL; }
#define CHK_FLOAT(o) \
    if (!PyFloat_Check(o)) { \
        PyErr_Format(PyExc_ValueError, "expected a float"); \
        return NULL; }
#define CHK_COMPLEX(o) \
    if (!PyComplex_Check(o)) { \
        PyErr_Format(PyExc_ValueError, "expected a complex"); \
        return NULL; }

// Define an error that we can have miriad throw to avoid exiting.
class MiriadError {
  private:
    std::string msg;
  public:
    MiriadError(const std::string &message) : msg (message) {}
    const char* get_message() const { return msg.c_str(); }
};

extern PyTypeObject UVType;

#endif
