#include "dsp.h"

// Adds data to a at indicies specified in ind.  Checks safety of arrays input.
PyObject *wrap_grid1D_c(PyObject *self, PyObject *args) {
    PyArrayObject *buf, *ind, *dat;
    int rv;
    long footprint=6;
    // Parse arguments and perform sanity check
    if (!PyArg_ParseTuple(args, "O!O!O!|i", &PyArray_Type, &buf, 
            &PyArray_Type, &ind, &PyArray_Type, &dat, &footprint)) 
        return NULL;
    CHK_ARRAY_RANK(buf, 1);
    CHK_ARRAY_RANK(ind, 1);
    CHK_ARRAY_RANK(dat, 1);
    CHK_ARRAY_TYPE(buf, NPY_CFLOAT);
    CHK_ARRAY_TYPE(ind, NPY_FLOAT);
    CHK_ARRAY_TYPE(dat, NPY_CFLOAT);
    if (PyArray_DIM(ind,0) != PyArray_DIM(dat,0)) {
        PyErr_Format(PyExc_ValueError, "Dimensions of ind and dat do not match");
        return NULL;
    }
        
    Py_INCREF(buf);
    Py_INCREF(ind);
    Py_INCREF(dat);
    rv = grid1D_c((float *) PyArray_DATA(buf), (long) PyArray_DIM(buf,0),
                  (float *) PyArray_DATA(ind), 
                  (float *) PyArray_DATA(dat), (long) PyArray_DIM(dat,0), footprint);
    Py_DECREF(buf);
    Py_DECREF(ind);
    Py_DECREF(dat);
    if (rv == 0) {
        Py_INCREF(Py_None);
        return Py_None;
    } else {
        PyErr_Format(PyExc_ValueError, "Invalid indices found.");
        return NULL;
    }
}

PyObject *wrap_grid2D_c(PyObject *self, PyObject *args) {
    PyArrayObject *buf, *ind1, *ind2, *dat;
    int rv;
    long footprint=6;
    // Parse arguments and perform sanity check
    if (!PyArg_ParseTuple(args, "O!O!O!O!|i", &PyArray_Type, &buf, 
            &PyArray_Type, &ind1, &PyArray_Type, &ind2, &PyArray_Type, &dat, &footprint)) 
        return NULL;
    CHK_ARRAY_RANK(buf, 2);
    CHK_ARRAY_RANK(ind1, 1);
    CHK_ARRAY_RANK(ind2, 1);
    CHK_ARRAY_RANK(dat, 1);
    CHK_ARRAY_TYPE(buf, NPY_CFLOAT);
    CHK_ARRAY_TYPE(ind1, NPY_FLOAT);
    CHK_ARRAY_TYPE(ind2, NPY_FLOAT);
    CHK_ARRAY_TYPE(dat, NPY_CFLOAT);
    if (PyArray_DIM(ind1,0) != PyArray_DIM(dat,0) || PyArray_DIM(ind2,0) != PyArray_DIM(dat,0)) {
        PyErr_Format(PyExc_ValueError, "Dimensions of ind and dat do not match");
        return NULL;
    }
        
    Py_INCREF(buf);
    Py_INCREF(ind1);
    Py_INCREF(ind2);
    Py_INCREF(dat);
    rv = grid2D_c((float *) PyArray_DATA(buf), (long) PyArray_DIM(buf,0), (long) PyArray_DIM(buf,1),
                  (float *) PyArray_DATA(ind1), (float *) PyArray_DATA(ind2), 
                  (float *) PyArray_DATA(dat), (long) PyArray_DIM(dat,0), footprint);
    Py_DECREF(buf);
    Py_DECREF(ind1);
    Py_DECREF(ind2);
    Py_DECREF(dat);
    if (rv == 0) {
        Py_INCREF(Py_None);
        return Py_None;
    } else {
        PyErr_Format(PyExc_ValueError, "Invalid indices found.");
        return NULL;
    }
}

PyObject *wrap_degrid2D_c(PyObject *self, PyObject *args) {
    PyArrayObject *buf, *ind1, *ind2, *dat;
    int rv;
    long footprint=6;
    // Parse arguments and perform sanity check
    if (!PyArg_ParseTuple(args, "O!O!O!O!|i", &PyArray_Type, &buf, 
            &PyArray_Type, &ind1, &PyArray_Type, &ind2, &PyArray_Type, &dat, &footprint)) 
        return NULL;
    CHK_ARRAY_RANK(buf, 2);
    CHK_ARRAY_RANK(ind1, 1);
    CHK_ARRAY_RANK(ind2, 1);
    CHK_ARRAY_RANK(dat, 1);
    CHK_ARRAY_TYPE(buf, NPY_CFLOAT);
    CHK_ARRAY_TYPE(ind1, NPY_FLOAT);
    CHK_ARRAY_TYPE(ind2, NPY_FLOAT);
    CHK_ARRAY_TYPE(dat, NPY_CFLOAT);
    if (PyArray_DIM(ind1,0) != PyArray_DIM(dat,0) || PyArray_DIM(ind2,0) != PyArray_DIM(dat,0)) {
        PyErr_Format(PyExc_ValueError, "Dimensions of ind and dat do not match");
        return NULL;
    }
        
    Py_INCREF(buf);
    Py_INCREF(ind1);
    Py_INCREF(ind2);
    Py_INCREF(dat);
    // Being lazy.  should allocate data rather than take it as an argument
    rv = degrid2D_c((float *) PyArray_DATA(buf), (long) PyArray_DIM(buf,0), (long) PyArray_DIM(buf,1),
                  (float *) PyArray_DATA(ind1), (float *) PyArray_DATA(ind2), 
                  (float *) PyArray_DATA(dat), (long) PyArray_DIM(dat,0), footprint);
    Py_DECREF(buf);
    Py_DECREF(ind1);
    Py_DECREF(ind2);
    Py_DECREF(dat);
    if (rv == 0) {
        Py_INCREF(Py_None);
        return Py_None;
    } else {
        PyErr_Format(PyExc_ValueError, "Invalid indices found.");
        return NULL;
    }
}

// Wrap function into module
static PyMethodDef _dsp_methods[] = {
    {"grid1D_c", (PyCFunction)wrap_grid1D_c, METH_VARARGS,
        "grid1D_c(buf,ind,dat,footprint=6)\nTBD."},
    {"grid2D_c", (PyCFunction)wrap_grid2D_c, METH_VARARGS,
        "grid2D_c(buf,ind1,ind2,dat,footprint=6)\nTBD."},
    {"degrid2D_c", (PyCFunction)wrap_degrid2D_c, METH_VARARGS,
        "degrid2D_c(buf,ind1,ind2,dat,footprint=6)\nTBD."},
    {NULL, NULL}
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC init_dsp(void) {
    (void) Py_InitModule("_dsp", _dsp_methods);
    import_array();
};
