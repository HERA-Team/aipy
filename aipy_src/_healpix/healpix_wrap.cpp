/*
 * This is a hand-written wrapper for Healpix.  Healpix_Base is wrapped into
 * numpy arrays, but instead of interfacing to Healpix_Map, numpy arrays
 * are used to hold map data.  This both makes the Healpix framework more
 * flexible (more data types available, etc.), and keeps data in the hands
 * of the programmer.
 *
 * Author: Aaron Parsons
 * Date: 12/05/07
 * Revisions:
 *      01/23/08    arp     bugfix on get_interpol for memory leak
 *      04/24/08    arp     moved interpol into crd2px functions
 */

#include <Python.h>
#include "numpy/arrayobject.h"
#include "healpix_base.h"
#include "healpix_map.h"
#include "arr.h"
#include "pointing.h"
#include "vec3.h"

#include <cmath>

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

// Some helper functions

/*____                           _                    _    
 / ___|_ __ ___  _   _ _ __   __| |_      _____  _ __| | __
| |  _| '__/ _ \| | | | '_ \ / _` \ \ /\ / / _ \| '__| |/ /
| |_| | | | (_) | |_| | | | | (_| |\ V  V / (_) | |  |   < 
 \____|_|  \___/ \__,_|_| |_|\__,_| \_/\_/ \___/|_|  |_|\_\
*/
// Python object that holds instance of Healpix_Base
typedef struct {
    PyObject_HEAD
    Healpix_Base hpb;
} HPBObject;

// Deallocate memory when Python object is deleted
static void HPBObject_dealloc(HPBObject* self) {
    self->ob_type->tp_free((PyObject*)self);
}

// Allocate memory for Python object and Healpix_Base (__new__)
static PyObject *HPBObject_new(PyTypeObject *type, 
        PyObject *args, PyObject *kwds) {
    HPBObject *self;
    self = (HPBObject *) type->tp_alloc(type, 0);
    return (PyObject *) self;
}

// Initialize object (__init__)
static int HPBObject_init(HPBObject *self, PyObject *args, PyObject *kwds) {
    int nside=-1;
    Healpix_Ordering_Scheme scheme = RING;
    PyObject *scheme_str=NULL;
    static char *kwlist[] = {"nside", "scheme", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds,"|iO", kwlist, \
            &nside, &scheme_str))
        return -1;
    if (scheme_str == NULL) scheme = RING;
    else if (strcmp(PyString_AsString(scheme_str), "NEST") == 0) scheme = NEST;
    else if (strcmp(PyString_AsString(scheme_str), "RING") != 0) {
        PyErr_Format(PyExc_ValueError, "scheme must be 'RING' or 'NEST'.");
        return -1;
    }
    // Carefully try to create a new Healpix_Base
    try {
        if (nside == -1) self->hpb = Healpix_Base();
        else self->hpb = Healpix_Base(nside, scheme, SET_NSIDE);
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
// Thin wrapper over Healpix_Base::npix2nside
static PyObject * HPBObject_npix2nside(HPBObject *self, PyObject *args) {
    int npix;
    if (!PyArg_ParseTuple(args, "i", &npix)) return NULL;
    try {
        return PyInt_FromLong(self->hpb.npix2nside(npix));
    } catch (Message_error &e) {
        PyErr_Format(PyExc_RuntimeError, e.what());
        return NULL;
    }
}

/* Wrapper over Healpix_Base::nest2ring and Healpix_Base::ring2nest
 * to convert an array of pixel indices into output order specified in
 * 'scheme'.  Modifies pixel array in place.
 */
static PyObject * HPBObject_nest_ring_conv(HPBObject *self, PyObject *args) {
    PyArrayObject *px;
    PyObject *scheme;
    // Parse and check input arguments
    if (!PyArg_ParseTuple(args, "O!O", &PyArray_Type, &px, &scheme))
        return NULL;
    CHK_ARRAY_TYPE(px,NPY_LONG);
    CHK_ARRAY_RANK(px,1);
    try {
        if (strcmp(PyString_AsString(scheme), "NEST") == 0) {
            for (int i=0; i < DIM(px,0); i++)
                IND1(px,i,long) = self->hpb.ring2nest(IND1(px,i,long));
        } else if (strcmp(PyString_AsString(scheme), "RING") == 0) {
            for (int i=0; i < DIM(px,0); i++)
                IND1(px,i,long) = self->hpb.nest2ring(IND1(px,i,long));
        } else {
            PyErr_Format(PyExc_ValueError,"scheme must be 'RING' or 'NEST'.");
            return NULL;
        }
    } catch (Message_error &e) {
        PyErr_Format(PyExc_RuntimeError, e.what());
        return NULL;
    }
    Py_INCREF(px);
    return PyArray_Return(px);
}

// Thin wrapper over Healpix_Base::SetNside
static PyObject * HPBObject_SetNside(HPBObject *self, PyObject *args) {
    Healpix_Ordering_Scheme hp_scheme = RING;
    int nside;
    PyObject *scheme = NULL;
    if (!PyArg_ParseTuple(args, "iO", &nside, &scheme)) return NULL;
    if (strcmp(PyString_AsString(scheme), "NEST") == 0) hp_scheme = NEST;
    else if (strcmp(PyString_AsString(scheme), "RING") != 0) {
        PyErr_Format(PyExc_ValueError, "scheme must be 'RING' or 'NEST'.");
        return NULL;
    }
    self->hpb.SetNside(nside, hp_scheme);
    Py_INCREF(Py_None);
    return Py_None;
}
    

/* Wraps ang2pix, and uses arrays to do many at once. */
static PyObject * HPBObject_crd2px(HPBObject *self, PyObject *args,
        PyObject *kwds) {
    int interpolate=0;
    double c1, c2, c3=0;
    fix_arr<int,4> fix_pix;
    fix_arr<double,4> fix_wgt;
    pointing p;
    vec3 v;
    PyArrayObject *crd1, *crd2, *crd3=NULL, *rv, *wgt=NULL;
    PyObject *rv2;
    static char *kwlist[] = {"crd1", "crd2", "crd3", "interpolate", NULL};
    // Parse and check input arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwds,"O!O!|O!i", kwlist, 
            &PyArray_Type, &crd1, &PyArray_Type, &crd2, &PyArray_Type, &crd3,
            &interpolate))
        return NULL;
    CHK_ARRAY_RANK(crd1,1);
    CHK_ARRAY_RANK(crd2,1);
    if (crd3 != NULL) CHK_ARRAY_RANK(crd3,1);
    int sz = DIM(crd1,0);
    if (DIM(crd2,0) != sz || (crd3 != NULL && DIM(crd3,0) != sz)) {
        PyErr_Format(PyExc_RuntimeError, "input crds must have same length.");
        return NULL;
    }
    CHK_ARRAY_TYPE(crd1, NPY_DOUBLE);
    CHK_ARRAY_TYPE(crd2, NPY_DOUBLE);
    if (crd3 != NULL) CHK_ARRAY_TYPE(crd3, NPY_DOUBLE);
    // Make array(s) to hold the results
    if (interpolate == 0) {
        npy_intp dimens[1] = {sz};
        rv = (PyArrayObject *) PyArray_SimpleNew(1, dimens, PyArray_LONG);
        CHK_NULL(rv);
    } else {
        npy_intp dimens[2] = {sz, 4};
        rv = (PyArrayObject *) PyArray_SimpleNew(2, dimens, PyArray_LONG);
        wgt = (PyArrayObject *) PyArray_SimpleNew(2, dimens, PyArray_DOUBLE);
        CHK_NULL(rv);
        CHK_NULL(wgt);
    }     
    // Interpret coordinates
    for (int i=0; i < sz; i++) {
        c1 = IND1(crd1,i,double);
        c2 = IND1(crd2,i,double);
        if (crd3 != NULL) c3 = IND1(crd3,i,double);
        if (!std::isfinite(c1) || !std::isfinite(c2) ||
                (crd3 != NULL && !std::isfinite(c3))) {
                    printf("Warning: encountered NaN/Inf in crd2px\n");
                    c1 = 0; c2 = 0; c3 = 1;
        }
        if (crd3 == NULL) {
            p.theta = c1; p.phi = c2;
        } else {
            v.x = c1; v.y = c2; v.z = c3;
        }
        if (interpolate == 0) {
            if (crd3 == NULL) IND1(rv,i,long) = self->hpb.ang2pix(p);
            else IND1(rv,i,long) = self->hpb.vec2pix(v);
        } else {    // Do interpolation
            if (crd3 != NULL) p = pointing(v);
            self->hpb.get_interpol(p, fix_pix, fix_wgt);
            for (int j=0; j < 4; j++) {
                IND2(rv,i,j,long) = fix_pix[j];
                IND2(wgt,i,j,double) = fix_wgt[j];
            }
        }
    }
    if (interpolate == 0) return PyArray_Return(rv);
    // Otherwise build tuple to return.
    // Make sure to DECREF when using Py_BuildValue() !!
    rv2 = Py_BuildValue("(OO)", PyArray_Return(rv), PyArray_Return(wgt));
    Py_DECREF(rv); Py_DECREF(wgt);
    return rv2;
    
}
    
/* Wraps pix2ang, but adds option of vector output as well.  Similarly
 * uses array I/O to do many at once.
 */
static PyObject * HPBObject_px2crd(HPBObject *self,
        PyObject *args, PyObject *kwds) {
    pointing p;
    vec3 v;
    int ncrd=3;
    PyArrayObject *px, *crd1, *crd2, *crd3;
    static char *kwlist[] = {"px", "ncrd", NULL};
    // Parse and check input arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwds,"O!|i", kwlist, 
            &PyArray_Type, &px, &ncrd))
        return NULL;
    if (ncrd != 2 && ncrd != 3) {
        PyErr_Format(PyExc_ValueError, "ncrd must be 2 or 3.");
        return NULL;
    }
    CHK_ARRAY_RANK(px,1);
    CHK_ARRAY_TYPE(px,NPY_LONG);
    // Make an array to hold the results
    int sz = px->dimensions[0];
    npy_intp dimens[1] = {sz};
    crd1 = (PyArrayObject *) PyArray_SimpleNew(1, dimens, PyArray_DOUBLE);
    crd2 = (PyArrayObject *) PyArray_SimpleNew(1, dimens, PyArray_DOUBLE);
    CHK_NULL(crd1);
    CHK_NULL(crd2);
    if (ncrd == 2) {
        for (int i=0; i < sz; i++) {
            p = self->hpb.pix2ang(IND1(px,i,int));
            IND1(crd1,i,double) = p.theta;
            IND1(crd2,i,double) = p.phi;
        }
        return Py_BuildValue("(OO)",PyArray_Return(crd1),PyArray_Return(crd2));
    } else {
        crd3 = (PyArrayObject *) PyArray_SimpleNew(1, dimens, PyArray_DOUBLE);
        for (int i=0; i < sz; i++) {
            p = self->hpb.pix2ang(IND1(px,i,int));
            v = p.to_vec3();
            IND1(crd1,i,double) = v.x;
            IND1(crd2,i,double) = v.y;
            IND1(crd3,i,double) = v.z;
        }
        return Py_BuildValue("(OOO)", PyArray_Return(crd1),
            PyArray_Return(crd2), PyArray_Return(crd3));
    }
}
        
// Thin wrapper over Healpix_Base::Order
static PyObject * HPBObject_Order(HPBObject *self) {
    return PyInt_FromLong(self->hpb.Order());
}

// Thin wrapper over Healpix_Base::Nside
static PyObject * HPBObject_Nside(HPBObject *self) {
    return PyInt_FromLong(self->hpb.Nside());
}

// Thin wrapper over Healpix_Base::Npix
static PyObject * HPBObject_Npix(HPBObject *self) {
    return PyInt_FromLong(self->hpb.Npix());
}

// Thin wrapper over Healpix_Base::Scheme
static PyObject * HPBObject_Scheme(HPBObject *self) {
    Healpix_Ordering_Scheme scheme = self->hpb.Scheme();
    if (scheme == RING) return PyString_FromString("RING");
    return PyString_FromString("NEST");
}

/*_        __                     _               _   _       
 \ \      / / __ __ _ _ __  _ __ (_)_ __   __ _  | | | |_ __  
  \ \ /\ / / '__/ _` | '_ \| '_ \| | '_ \ / _` | | | | | '_ \ 
   \ V  V /| | | (_| | |_) | |_) | | | | | (_| | | |_| | |_) |
    \_/\_/ |_|  \__,_| .__/| .__/|_|_| |_|\__, |  \___/| .__/ 
                     |_|   |_|            |___/        |_|    
*/
// Bind methods to object
static PyMethodDef HPBObject_methods[] = {
    {"npix2nside", (PyCFunction)HPBObject_npix2nside, METH_VARARGS,
        "npix2nside(npix)\nConvert number of pixels to number of sides."},
    {"nest_ring_conv", (PyCFunction)HPBObject_nest_ring_conv,
        METH_VARARGS,
        "nest_ring_conv(px,scheme)\nTranslate an array of pixel numbers to index data in the scheme specified in 'scheme' ('NEST' or 'RING').  Returns px, which has been modified in place (so beware!)."},
    {"set_nside_scheme", (PyCFunction)HPBObject_SetNside, METH_VARARGS,
        "set_nside_scheme(nside,scheme)\nAdjust Nside and Scheme ('RING' or 'NEST')."},
    {"crd2px", (PyCFunction)HPBObject_crd2px, METH_VARARGS|METH_KEYWORDS,
        "crd2px(c1,c2,c3=None,interpolate=False)\nConvert 1 dimensional arrays of input coordinates to pixel indices. If only c1,c2 provided, then read them as th,phi.  If c1,c2,c3 provided, read them as x,y,z. If interpolate is False, return a single pixel coordinate.  If interpolate is True, return px,wgts where each entry in px contains the 4 pixels adjacent to the specified location, and wgt contains the 4 corresponding weights of those pixels."},
    {"px2crd", (PyCFunction)HPBObject_px2crd,METH_VARARGS|METH_KEYWORDS,
        "px2crd(px,ncrd=3)\nConvert a 1 dimensional input array of pixel numbers to the type of coordinates specified by ncrd.  If ncrd=3 (default), the returned array will have (x,y,z) for each pixel.  Otherwise if ncrd=2, the returned array will have (theta,phi) for each pixel."},
    {"order", (PyCFunction)HPBObject_Order,METH_NOARGS,
        "order()\nReturn the order parameter."},
    {"nside", (PyCFunction)HPBObject_Nside,METH_NOARGS,
        "nside()\nReturn the Nside parameter."},
    {"npix", (PyCFunction)HPBObject_Npix,METH_NOARGS,
        "npix()\nReturn the number of pixels in the map."},
    {"scheme", (PyCFunction)HPBObject_Scheme,METH_NOARGS,
        "scheme()\nReturn the scheme of the map ('NEST' or 'RING')."},
    {NULL}  /* Sentinel */
};

static PyTypeObject HPBType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_healpix.HPM", /*tp_name*/
    sizeof(HPBObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)HPBObject_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "Functionality related to the HEALPix pixelisation.  HealpixBase() or HealpixBase(nside, scheme='RING').",       /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    HPBObject_methods,             /* tp_methods */
    0,                     /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)HPBObject_init,      /* tp_init */
    0,                         /* tp_alloc */
    HPBObject_new,       /* tp_new */
};

// Module methods
static PyMethodDef _healpix_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

// Module init
PyMODINIT_FUNC init_healpix(void) {
    PyObject* m;
    HPBType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&HPBType) < 0) return;
    m = Py_InitModule3("_healpix", _healpix_methods,
    "This is a hand-written wrapper (by Aaron Parsons) for Healpix_cxx, which was developed at the Max-Planck-Institut für Astrophysik and financially supported by the Deutsches Zentrum für Luft- und Raumfahrt (DLR).");
    import_array();
    Py_INCREF(&HPBType);
    PyModule_AddObject(m, "HealpixBase", (PyObject *)&HPBType);
}

