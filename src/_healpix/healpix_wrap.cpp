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
    PyObject *order = NULL;
    static char *options[] = {"RING", "NEST", NULL};
    static char *kwlist[] = {"nside", "ordering", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds,"|iO", kwlist, &nside, &order))
        return -1;
    int use_nest = get_option(options, order);
    if (use_nest == -1) return -1;
    Healpix_Ordering_Scheme scheme = RING;
    if (use_nest == 1) scheme = NEST;
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
 * 'ordering'.  Modifies pixel array in place.
 */
static PyObject * HPBObject_nest_ring_conv(HPBObject *self, PyObject *args) {
    PyArrayObject *px;
    PyObject *ordering;
    // Parse and check input arguments
    if (!PyArg_ParseTuple(args, "O!O", &PyArray_Type, &px, &ordering))
        return NULL;
    CHK_ARRAY_TYPE(px,NPY_LONG);
    CHK_ARRAY_RANK(px,1);
    try {
        if (strcmp(PyString_AsString(ordering), "NEST") == 0) {
            for (int i=0; i < DIM(px,0); i++)
                IND1(px,i,long) = self->hpb.ring2nest(IND1(px,i,long));
        } else if (strcmp(PyString_AsString(ordering), "RING") == 0) {
            for (int i=0; i < DIM(px,0); i++)
                IND1(px,i,long) = self->hpb.nest2ring(IND1(px,i,long));
        } else {
            PyErr_Format(PyExc_ValueError,"ordering must be 'RING' or 'NEST'.");
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
    Healpix_Ordering_Scheme scheme = RING;
    int nside;
    PyObject *ordering = NULL;
    if (!PyArg_ParseTuple(args, "iO", &nside, &ordering)) return NULL;
    if (strcmp(PyString_AsString(ordering), "NEST") == 0) scheme = NEST;
    else if (strcmp(PyString_AsString(ordering), "RING") != 0) {
        PyErr_Format(PyExc_ValueError, "ordering must be 'RING' or 'NEST'.");
        return NULL;
    }
    self->hpb.SetNside(nside, scheme);
    Py_INCREF(Py_None);
    return Py_None;
}
    

/* Wraps together ang2pix and vec2pix for a unified "coordinate input", and 
 * uses input array to do many at once.  If 2nd axis of input array has
 * length 2, it is interpreted as angles for a Pointing.  If it has length
 * 3, it is interpreted as x,y,z for a Vector.
 */
static PyObject * HPBObject_crd2px(HPBObject *self, PyObject *args) {
    PyArrayObject *crd, *rv;
    // Parse and check input arguments
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &crd)) return NULL;
    CHK_ARRAY_RANK(crd,2);
    if (DIM(crd,1) != 2 and DIM(crd,1) != 3) {
        PyErr_Format(PyExc_RuntimeError, "2nd axis must have length 2 or 3.");
        return NULL;
    }
    CHK_ARRAY_TYPE(crd, NPY_DOUBLE);
    // Make an array to hold the results
    int sz = DIM(crd,0);
    int dimens[1] = {sz};
    rv = (PyArrayObject *) PyArray_FromDims(1, dimens, PyArray_LONG);
    CHK_NULL(rv);
    // Interpret coordinates
    if (DIM(crd,1) == 2) {
        pointing p;
        for (int i=0; i < sz; i++) {
            p.theta = IND2(crd,i,0,double);
            p.phi = IND2(crd,i,1,double);
            if (std::isfinite(p.theta) && std::isfinite(p.phi))
                IND1(rv,i,long) = self->hpb.ang2pix(p);
            else {
                printf("Warning: encountered NaN/Inf in crd2px\n");
                IND1(rv,i,long) = 0;
            }
        }
    } else { // DIM(crd,1) == 3
        vec3 v;
        for (int i=0; i < sz; i++) {
            v.x = IND2(crd,i,0,double);
            v.y = IND2(crd,i,1,double);
            v.z = IND2(crd,i,2,double);
            if (std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z))
                IND1(rv,i,long) = self->hpb.vec2pix(v);
            else {
                printf("Warning: encountered NaN/Inf in crd2px\n");
                IND1(rv,i,long) = 0;
            }
        }
    }
    return PyArray_Return(rv);
}

/* Wraps pix2ang, but adds possibility of vector output as well.  Similarly
 * uses array I/O to do many at once.
 */
static PyObject * HPBObject_px2crd(HPBObject *self,
        PyObject *args, PyObject *kwds) {
    pointing p;
    vec3 v;
    PyArrayObject *px, *rv;
    static char *options[] = {"pnt", "vec", NULL};
    static char *kwlist[] = {"px", "crd_type", NULL};
    PyObject *crd_type = NULL;
    // Parse and check input arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwds,"O!|O", kwlist, 
            &PyArray_Type, &px, &crd_type))
        return NULL;
    int return_type = get_option(options, crd_type);
    printf("Return type: %d = %s\n", return_type, options[return_type]);
    if (return_type == -1) return NULL;
    CHK_ARRAY_RANK(px,1);
    CHK_ARRAY_TYPE(px,NPY_LONG);
    // Make an array to hold the results
    int sz = px->dimensions[0];
    if (return_type == 0) {
        int dimens[2] = {sz, 2};
        rv = (PyArrayObject *) PyArray_FromDims(2, dimens, PyArray_DOUBLE);
    } else {
        int dimens[2] = {sz, 3};
        rv = (PyArrayObject *) PyArray_FromDims(2, dimens, PyArray_DOUBLE);
    }
    CHK_NULL(rv);
    // Interpret coordinates
    for (int i=0; i < sz; i++) {
        p = self->hpb.pix2ang(IND1(px,i,int));
        if (return_type == 0) {
            IND2(rv,i,0,double) = p.theta;
            IND2(rv,i,1,double) = p.phi;
        } else {
            v = p.to_vec3();
            IND2(rv,i,0,double) = v.x;
            IND2(rv,i,1,double) = v.y;
            IND2(rv,i,2,double) = v.z;
        }
    }
    return PyArray_Return(rv);
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

/* Wraps get_interpol and get_interpol2 to return nearby pixels and weights
 * for an array of input coordinates, allowing you to choose between them for 
 * RING and NEST schemes.  
 */
static PyObject * HPBObject_get_interpol(HPBObject *self,
        PyObject *args, PyObject *kwds) {
    pointing p;
    vec3 v;
    fix_arr<int,4> fix_pix;
    fix_arr<double,4> fix_wgt;
    PyArrayObject *crd, *px, *wgt;
    PyObject *rv = NULL;
    // Parse and check input arguments
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &crd)) return NULL;
    CHK_ARRAY_RANK2(crd,1,2);
    if (RANK(crd) == 1) {
        CHK_ARRAY_TYPE(crd, NPY_LONG);
    } else if (DIM(crd,1) != 2 and DIM(crd,1) != 3) {
        PyErr_Format(PyExc_RuntimeError, "2nd dim must be 2 or 3");
        return NULL;
    } else {
        CHK_ARRAY_TYPE(crd, NPY_DOUBLE);
    }
    // Make arrays to hold the results
    int sz = DIM(crd,0);
    int dimens[2] = {sz, 4};
    px = (PyArrayObject *) PyArray_FromDims(2, dimens, PyArray_LONG);
    wgt = (PyArrayObject *) PyArray_FromDims(2, dimens, PyArray_DOUBLE);
    CHK_NULL(px);
    CHK_NULL(wgt);
    // Interpret coordinates
    for (int i=0; i < sz; i++) {
        if (RANK(crd) == 1) {
            p = self->hpb.pix2ang(IND1(crd,i,long));
        } else if (DIM(crd,1) == 2) {
            p.theta = IND2(crd,i,0,double);
            p.phi = IND2(crd,i,1,double);
            if (!std::isfinite(p.theta) || !std::isfinite(p.phi)) {
                printf("Warning: encountered NaN/Inf in get_interpol\n");
                p.theta = 0; p.phi = 0;
            }
        } else { // DIM(crd,1) == 3
            v.x = IND2(crd,i,0,double);
            v.y = IND2(crd,i,1,double);
            v.z = IND2(crd,i,2,double);
            if (std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z))
                p = pointing(v);
            else {
                printf("Warning: encountered NaN/Inf in get_interpol\n");
                p.theta = 0; p.phi = 0;
            }
        }
        self->hpb.get_interpol(p, fix_pix, fix_wgt);
        for (int j=0; j < 4; j++) {
            IND2(px,i,j,long) = fix_pix[j];
            IND2(wgt,i,j,double) = fix_wgt[j];
        }
    }
    rv = PyTuple_Pack(2, (PyObject *) px, (PyObject *) wgt);
    // We now have duplicate references, so DECREF to avoid memory leak
    Py_DECREF(px); Py_DECREF(wgt);
    return rv;
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
        "nest_ring_conv(px, ordering)\nTranslate an array of pixel numbers to index data in the scheme specified in 'ordering' ('NEST' or 'RING').  Returns px, which has been modified in place (so beware!)."},
    {"SetNside", (PyCFunction)HPBObject_SetNside, METH_VARARGS,
        "SetNside(nside, ordering)\nAdjust Nside and ordering scheme ('RING' or 'NEST')."},
    {"crd2px", (PyCFunction)HPBObject_crd2px, METH_VARARGS,
        "crd2px(crd)\nConvert a 2 dimensional input array of coordinates to pixel indices.  If 2nd axis has length 2, it is interpreted as (theta,phi) coordinates on a sphere.  If 2nd axis has length 3, it is interpreted as (x,y,z) coordinates, which will point to a pixel on the sphere."},
    {"px2crd", (PyCFunction)HPBObject_px2crd,METH_VARARGS|METH_KEYWORDS,
        "px2crd(px, crd_type='vec')\nConvert a 1 dimensional input array of pixel numbers to the type of coordinates specified by 'crd_type'.  If crd_type='vec', the returned array will have (x,y,z) for each pixel.  Otherwise if crd_type='pnt', the returned array will have (theta,phi) for each pixel."},
    {"Order", (PyCFunction)HPBObject_Order,METH_NOARGS,
        "Order()\nReturn the order parameter."},
    {"Nside", (PyCFunction)HPBObject_Nside,METH_NOARGS,
        "Nside()\nReturn the Nside parameter."},
    {"Npix", (PyCFunction)HPBObject_Npix,METH_NOARGS,
        "Npix()\nReturn the number of pixels in the map."},
    {"Scheme", (PyCFunction)HPBObject_Scheme,METH_NOARGS,
        "Scheme()\nReturn the ordering scheme of the map ('NEST' or 'RING')."},
    {"get_interpol", (PyCFunction)HPBObject_get_interpol,
        METH_VARARGS|METH_KEYWORDS,
        "get_interpol(crd, ordering='RING')\nReturn px,wgts where each entry in px contains the 4 pixels adjacent to the coordinates in 'crd', and wgt contains the 4 corresponding weights of those pixels.  'ordering' should generally represent the current ordering scheme of the map for most efficient computation (although it will still give you the right answer otherwise)."},
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
    "Functionality related to the HEALPix pixelisation.  HPM() or HPM(nside, ordering='RING').",       /* tp_doc */
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
    "This is a hand-written wrapper (by Aaron Parsons) for Healpix_cxx, which was developed at the Max-Planck-Institut fuer Astrophysik and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt (DLR).");
    import_array();
    Py_INCREF(&HPBType);
    PyModule_AddObject(m, "HealpixBase", (PyObject *)&HPBType);
}

