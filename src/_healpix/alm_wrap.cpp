/*
 * Author: Aaron Parsons
 * Date: 03/24/08
 * Revisions:
 */

#include <Python.h>
#include "numpy/arrayobject.h"
#include "alm.h"
#include "healpix_base.h"
#include "healpix_map.h"
#include "alm_map_tools.h"
#include "xcomplex.h"

#if PY_VERSION_HEX < 0x02050000
#define lenfunc inquiry
#endif

// Some macros...
#define QUOTE(a) # a

#define IND1(a,i,type) *((type *)(a->data + i*a->strides[0]))
#define IND2(a,i,j,type) *((type *)(a->data+i*a->strides[0]+j*a->strides[1]))

#define TYPE(a) a->descr->type_num
#define CHK_ARRAY_TYPE(a,type) \
    if (TYPE(a) != type) { \
        PyErr_Format(PyExc_ValueError, "type(%s) != %s", \
        QUOTE(a), QUOTE(type)); \
        return NULL; }

#define DIM(a,i) a->dimensions[i]

#define RANK(a) a->nd
#define CHK_ARRAY_RANK(a,r) \
    if (RANK(a) != r) { \
        PyErr_Format(PyExc_ValueError, "rank(%s) != %s", \
        QUOTE(a), QUOTE(r)); \
        return NULL; }

#define CHK_ARRAY_RANK2(a,r1,r2) \
    if (RANK(a) != r1 and RANK(a) != r2) { \
        PyErr_Format(PyExc_ValueError, "rank(%s) != %s or %s", \
        QUOTE(a), QUOTE(r1), QUOTE(r2)); \
        return NULL; }

#define CHK_NULL(a) \
    if (a == NULL) { \
        PyErr_Format(PyExc_MemoryError, "Failed to allocate %s", QUOTE(a)); \
        return NULL; }

#define CHK_COMPLEX(o) \
    if (!PyComplex_Check(o)) { \
        PyErr_Format(PyExc_ValueError, "expected a complex"); \
        return NULL; }

// Some helper functions

void option_err(char *options[]) {
    // Raise a python error if invalid option was provided
    char errstr[256];
    int i = 1;
    strcpy(errstr, "option not in ["); strcat(errstr, options[0]);
    while (options[i] != NULL) {
        strcat(errstr, ",");
        strcat(errstr, options[i]);
        i++;
    }
    strcat(errstr, "]");
    PyErr_Format(PyExc_ValueError, errstr);
}
    
int get_option(char *options[], PyObject *choice) {
    // Return index of choice in options, -1 if failure.  Defaults to 0.
    int i = 0;
    char *c;
    if (choice == NULL) return 0;
    if (!PyString_Check(choice)) {
        option_err(options);
        return -1;
    }
    c = PyString_AsString(choice);
    while (options[i] != NULL) {
        if (strcmp(c, options[i]) == 0) return i;
        i++;
    }
    option_err(options);
    return -1;
}


/*____                           _                    _    
 / ___|_ __ ___  _   _ _ __   __| |_      _____  _ __| | __
| |  _| '__/ _ \| | | | '_ \ / _` \ \ /\ / / _ \| '__| |/ /
| |_| | | | (_) | |_| | | | | (_| |\ V  V / (_) | |  |   < 
 \____|_|  \___/ \__,_|_| |_|\__,_| \_/\_/ \___/|_|  |_|\_\
*/
// Python object that holds instance of Healpix_Base
typedef struct {
    PyObject_HEAD
    Alm<xcomplex<double> > alm;
} AlmObject;

// Deallocate memory when Python object is deleted
static void AlmObject_dealloc(AlmObject* self) {
    self->ob_type->tp_free((PyObject*)self);
}

// Allocate memory for Python object
static PyObject *AlmObject_new(PyTypeObject *type, 
        PyObject *args, PyObject *kwds) {
    AlmObject *self;
    self = (AlmObject *) type->tp_alloc(type, 0);
    return (PyObject *) self;
}

// Initialize object (__init__)
static int AlmObject_init(AlmObject *self, PyObject *args, PyObject *kwds) {
    int lmax=0, mmax=0;
    if (!PyArg_ParseTuple(args, "ii", &lmax, &mmax))
        return -1;
    // Carefully try to create a new healpix_base
    try {
        self->alm = Alm<xcomplex<double> >(lmax, mmax);
        self->alm.SetToZero();
    } catch (Message_error &e) {
        PyErr_Format(PyExc_RuntimeError, e.what());
        return -1;
    }
    return 0;
}

/* ___  _     _           _     __  __      _   _               _     
  / _ \| |__ (_) ___  ___| |_  |  \/  | ___| |_| |__   ___   __| |___ 
 | | | | '_ \| |/ _ \/ __| __| | |\/| |/ _ \ __| '_ \ / _ \ / _` / __|
 | |_| | |_) | |  __/ (__| |_  | |  | |  __/ |_| | | | (_) | (_| \__ \
  \___/|_.__// |\___|\___|\__| |_|  |_|\___|\__|_| |_|\___/ \__,_|___/
           |__/                                                       
*/
// Get a coefficient
static PyObject * AlmObject_get(AlmObject *self, PyObject *args) {
    int l, m, lmax=self->alm.Lmax(), mmax=self->alm.Mmax();
    xcomplex<double> c;
    if (!PyArg_ParseTuple(args, "ii", &l, &m)) return NULL;
    if (l < 0 || l > lmax || m < 0 || m > mmax || m > l) {
        PyErr_Format(PyExc_RuntimeError, "Index out of range");
        return NULL;
    }
    try {
        c = self->alm(l,m);
        return PyComplex_FromDoubles((double)c.real(), (double)c.imag());
    } catch (Message_error &e) {
        PyErr_Format(PyExc_RuntimeError, e.what());
        return NULL;
    }
}

// Set a coefficient
static int AlmObject_set(AlmObject *self, PyObject *ind, PyObject *val) {
    int l, m, lmax=self->alm.Lmax(), mmax=self->alm.Mmax();
    xcomplex<double> c;
    if (!PyArg_ParseTuple(ind, "ii", &l, &m)) return -1;
    if (l < 0 || l > lmax || m < 0 || m > mmax || m > l) {
        PyErr_Format(PyExc_RuntimeError, "Index out of range");
        return -1;
    }
    if (PyComplex_Check(val)) {
        c.Set(PyComplex_RealAsDouble(val), PyComplex_ImagAsDouble(val));
    } else if (PyFloat_Check(val)) {
        c.Set(PyFloat_AsDouble(val), 0);
    } else if (PyInt_Check(val)) {
        c.Set((double) PyInt_AsLong(val), 0);
    } else {
        PyErr_Format(PyExc_ValueError, "Could not convert value to complex");
        return -1;
    }
    try {
        self->alm(l,m) = c;
        return 0;
    } catch (Message_error &e) {
        PyErr_Format(PyExc_RuntimeError, e.what());
        return -1;
    }
}

// Get maximum L
static PyObject * AlmObject_lmax(AlmObject *self) {
    return PyInt_FromLong((long) self->alm.Lmax());
}

// Get maximum M
static PyObject * AlmObject_mmax(AlmObject *self) {
    return PyInt_FromLong((long) self->alm.Mmax());
}

// Clear all coefficients
static PyObject * AlmObject_set_to_zero(AlmObject *self) {
    self->alm.SetToZero();
    Py_INCREF(Py_None);
    return Py_None;
}

// Convert to a Healpix Map
static PyObject * AlmObject_to_map(AlmObject *self, PyObject *args) {
    PyObject *ordering;
    Healpix_Ordering_Scheme scheme = RING;
    int nside;
    npy_intp npix;
    PyArrayObject *rv;
    if (!PyArg_ParseTuple(args, "iO", &nside, &ordering)) return NULL;
    if (strcmp(PyString_AsString(ordering), "NEST") == 0) {
        scheme = NEST;
    } else if (strcmp(PyString_AsString(ordering), "RING") != 0) {
        PyErr_Format(PyExc_ValueError,"ordering must be 'RING' or 'NEST'.");
        return NULL;
    }
    try {
        Healpix_Map<double> map(nside, scheme, SET_NSIDE);
        alm2map<double>(self->alm, map);
        // Transfer map contents into numpy array
        npix = map.Npix();
        rv = (PyArrayObject *) PyArray_SimpleNew(1, &npix, PyArray_DOUBLE);
        CHK_NULL(rv);
        for (int i=0; i < npix; i++) IND1(rv,i,double) = map[i];
        return PyArray_Return(rv);
    } catch (Message_error &e) {
        PyErr_Format(PyExc_RuntimeError, e.what());
        return NULL;
    }
}

// Compute coeffs of a Healpix Map
static PyObject * AlmObject_from_map(AlmObject *self, PyObject *args) {
    int iter, nside, npix;
    PyArrayObject *data;
    Healpix_Base hpb;
    if (!PyArg_ParseTuple(args, "O!i", &PyArray_Type, &data, &iter))
        return NULL;
    if (RANK(data) != 1) {
        PyErr_Format(PyExc_ValueError, "data must have 1 dimension.");
        return NULL;
    }
    nside = hpb.npix2nside(DIM(data,0));
    try {
        Healpix_Map<double> map(nside, RING, SET_NSIDE);
        npix = map.Npix();
        if (TYPE(data) == NPY_BOOL) {
            for (int i=0; i < npix; i++) map[i] = (double) IND1(data,i,bool);
        } else if (TYPE(data) == NPY_SHORT) {
            for (int i=0; i < npix; i++) map[i] = (double) IND1(data,i,short);
        } else if (TYPE(data) == NPY_USHORT) {
            for (int i=0; i < npix; i++) map[i] = (double) IND1(data,i,unsigned short);
        } else if (TYPE(data) == NPY_INT) {
            for (int i=0; i < npix; i++) map[i] = (double) IND1(data,i,int);
        } else if (TYPE(data) == NPY_UINT) {
            for (int i=0; i < npix; i++) map[i] = (double) IND1(data,i,unsigned int);
        } else if (TYPE(data) == NPY_LONG) {
            for (int i=0; i < npix; i++) map[i] = (double) IND1(data,i,long);
        } else if (TYPE(data) == NPY_ULONG) {
            for (int i=0; i < npix; i++) map[i] = (double) IND1(data,i,unsigned long);
        } else if (TYPE(data) == NPY_LONGLONG) {
            for (int i=0; i < npix; i++) map[i] = (double) IND1(data,i,long long);
        } else if (TYPE(data) == NPY_ULONGLONG) {
            for (int i=0; i < npix; i++) map[i] = (double) IND1(data,i,unsigned long long);
        } else if (TYPE(data) == NPY_FLOAT) {
            for (int i=0; i < npix; i++) map[i] = (double) IND1(data,i,float);
        } else if (TYPE(data) == NPY_DOUBLE) {
            for (int i=0; i < npix; i++) map[i] = (double) IND1(data,i,double);
        } else if (TYPE(data) == NPY_LONGDOUBLE) {
            for (int i=0; i < npix; i++) map[i] = (double) IND1(data,i,long double);
        } else {
            PyErr_Format(PyExc_ValueError, "Unsupported data type");
            return NULL;
        }
        map2alm_iter<double>(map, self->alm, iter);
        Py_INCREF(Py_None);
        return Py_None;
    } catch (Message_error &e) {
        PyErr_Format(PyExc_RuntimeError, e.what());
        return NULL;
    }
}

// Return the L,M indices of the data in an Alm
static PyObject * AlmObject_lm_indices(AlmObject *self) {
    PyArrayObject *L, *M;
    int lmax = self->alm.Lmax(), mmax = self->alm.Mmax();
    npy_intp num_alms = ((mmax+1)*(mmax+2))/2 + (mmax+1)*(lmax-mmax);
    int cnt=0;
    L = (PyArrayObject *) PyArray_SimpleNew(1, &num_alms, PyArray_INT);
    M = (PyArrayObject *) PyArray_SimpleNew(1, &num_alms, PyArray_INT);
    CHK_NULL(L); CHK_NULL(M);
    for (int l=0; l<=lmax; ++l) {
        for (int m=0; m<=mmax; ++m) {
            if (m > l) break;
            IND1(L,cnt,int) = l;
            IND1(M,cnt,int) = m;
            cnt++;
        }
    }
    return Py_BuildValue("(OO)", PyArray_Return(L), PyArray_Return(M));
}

// Return the coefficient data in an Alm
static PyObject * AlmObject_get_data(AlmObject *self) {
    PyArrayObject *rv;
    int lmax = self->alm.Lmax(), mmax = self->alm.Mmax();
    npy_intp num_alms = ((mmax+1)*(mmax+2))/2 + (mmax+1)*(lmax-mmax);
    xcomplex<double> c;
    int cnt=0;
    rv = (PyArrayObject *) PyArray_SimpleNew(1, &num_alms, PyArray_CDOUBLE);
    CHK_NULL(rv);
    for (int l=0; l<=lmax; ++l) {
        for (int m=0; m<=mmax; ++m) {
            if (m > l) break;
            c = self->alm(l,m);
            *((double *)(rv->data + cnt*rv->strides[0])) = (double) c.real();
            *((double *)(rv->data + cnt*rv->strides[0] + sizeof(double))) = (double) c.imag();
            cnt++;
        }
    }
    return PyArray_Return(rv);
}

// Set the coefficient data in an Alm
static PyObject * AlmObject_set_data(AlmObject *self, PyObject *args) {
    PyArrayObject *data;
    int lmax = self->alm.Lmax(), mmax = self->alm.Mmax();
    int num_alms = ((mmax+1)*(mmax+2))/2 + (mmax+1)*(lmax-mmax);
    int cnt=0;
    xcomplex<double> c;
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &data)) return NULL;
    if (RANK(data) != 1 || DIM(data,0) != num_alms) {
        PyErr_Format(PyExc_ValueError, "data must have length %d.", num_alms);
        return NULL;
    }
    CHK_ARRAY_TYPE(data, NPY_CDOUBLE);
    for (int l=0; l<=lmax; ++l) {
        for (int m=0; m<=mmax; ++m) {
            if (m > l) break;
            c.Set(*((double *)(data->data + cnt*data->strides[0])),
             *((double *)(data->data + cnt*data->strides[0] + sizeof(double))));
            self->alm(l,m) = c;
            cnt++;
        }
    }
    Py_INCREF(Py_None);
    return Py_None;
}

/*_        __                     _               _   _       
 \ \      / / __ __ _ _ __  _ __ (_)_ __   __ _  | | | |_ __  
  \ \ /\ / / '__/ _` | '_ \| '_ \| | '_ \ / _` | | | | | '_ \ 
   \ V  V /| | | (_| | |_) | |_) | | | | | (_| | | |_| | |_) |
    \_/\_/ |_|  \__,_| .__/| .__/|_|_| |_|\__, |  \___/| .__/ 
                     |_|   |_|            |___/        |_|    
*/
// Bind methods to object
static PyMethodDef AlmObject_methods[] = {
    {"set_to_zero", (PyCFunction)AlmObject_set_to_zero, METH_VARARGS,
        "set_to_zero()\nClear all coefficients."},
    {"lmax", (PyCFunction)AlmObject_lmax, METH_VARARGS,
        "lmax()\nReturn the maximum L."},
    {"mmax", (PyCFunction)AlmObject_mmax, METH_VARARGS,
        "mmax()\nReturn the maximum M."},
    {"to_map", (PyCFunction)AlmObject_to_map, METH_VARARGS,
        "to_map(nside, scheme)\nReturn data for the HealpixMap with the specified nside (power of 2) and ordering scheme ('RING' OR 'NEST') that is generated by these Alm coefficients."},
    {"from_map", (PyCFunction)AlmObject_from_map, METH_VARARGS,
        "from_map(data, iter)\nSet the coefficients of this Alm object (with its specified lmax and mmax) from the data of a HealpixMap in 'RING' mode.  Greater accuracy can be achieved with higher values of iter (iter=1 uses the fastest computation)."},
    {"lm_indices", (PyCFunction)AlmObject_lm_indices, METH_VARARGS,
        "lm_indices()\nReturn the L and M indices of the coefficients contained in Alm, in the order that they are returned by get_data()."},
    {"get_data", (PyCFunction)AlmObject_get_data, METH_VARARGS,
        "get_data()\nReturn all of the coefficients contained in Alm, in the order indexed by lm_indices()."},
    {"set_data", (PyCFunction)AlmObject_set_data, METH_VARARGS,
        "set_data()\nSet all of coefficients contained in Alm, in the order indexed by lm_indices()."},
    {NULL}  /* Sentinel */
};

static PyMappingMethods AlmObject_as_mapping = {
    (lenfunc)NULL,                  /*mp_length*/
    (binaryfunc)AlmObject_get,      /*mp_subscript*/
    (objobjargproc)AlmObject_set,   /*mp_ass_subscript*/
};

static PyTypeObject AlmType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_alm.Alm", /*tp_name*/
    sizeof(AlmObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)AlmObject_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    &AlmObject_as_mapping,     /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "Alm(lmax,mmax)\n Holds coefficients for spherical harmonics up to the specified order, and generates a real-valued HealpixMap from them.  Individual (l,m) coefficients can be accessed by alm[l,m].",       /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    AlmObject_methods,             /* tp_methods */
    0,                     /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)AlmObject_init,  /* tp_init */
    0,                         /* tp_alloc */
    AlmObject_new,       /* tp_new */
};

// Module methods
static PyMethodDef _alm_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

// Module init
PyMODINIT_FUNC init_alm(void) {
    PyObject* m;
    AlmType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&AlmType) < 0) return;
    m = Py_InitModule3("_alm", _alm_methods,
    "This is a hand-written wrapper (by Aaron Parsons) for Healpix_cxx, which was developed at the Max-Planck-Institut fuer Astrophysik and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt (DLR).");
    import_array();
    Py_INCREF(&AlmType);
    PyModule_AddObject(m, "Alm", (PyObject *)&AlmType);
}

