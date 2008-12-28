/*
 * Some additional deconvolution functions for AIPY, written in C++.  These are
 * mostly for speed-critical functions. 
 *
 * Author: Aaron Parsons
 * Date: 05/08/08
 * Revisions:
 */

#include <Python.h>
#include "numpy/arrayobject.h"

#define QUOTE(s) # s

#define PNT1(a,i) (a->data + i*a->strides[0])
#define PNT2(a,i,j) (a->data+i*a->strides[0]+j*a->strides[1])
#define IND1(a,i,type) *((type *)PNT1(a,i))
#define IND2(a,i,j,type) *((type *)PNT2(a,i,j))
#define CIND1R(a,i,type) *((type *)PNT1(a,i))
#define CIND1I(a,i,type) *((type *)(PNT1(a,i)+sizeof(type)))
#define CIND2R(a,i,j,type) *((type *)PNT2(a,i,j))
#define CIND2I(a,i,j,type) *((type *)(PNT2(a,i,j)+sizeof(type)))

#define TYPE(a) a->descr->type_num
#define CHK_ARRAY_TYPE(a,type) \
    if (TYPE(a) != type) { \
        PyErr_Format(PyExc_ValueError, "type(%s) != %s", \
        QUOTE(a), QUOTE(type)); \
        return NULL; }

#define DIM(a,i) a->dimensions[i]
#define CHK_ARRAY_DIM(a,i,d) \
    if (DIM(a,i) != d) { \
        PyErr_Format(PyExc_ValueError, "dim(%s) != %s", \
        QUOTE(a), QUOTE(d)); \
        return NULL; }

#define RANK(a) a->nd
#define CHK_ARRAY_RANK(a,r) \
    if (RANK(a) != r) { \
        PyErr_Format(PyExc_ValueError, "rank(%s) != %s", \
        QUOTE(a), QUOTE(r)); \
        return NULL; }


// A template for implementing addition loops for different data types
template<typename T> struct Clean {
    // Does a 1d real-valued clean
    static int clean_1d_r(PyArrayObject *res, PyArrayObject *ker, 
            PyArrayObject *mdl, double gain, int maxiter, double tol) {
        T score=-1, nscore, max=0, mmax, val, mval, step, q=0, mq=0;
        int argmax=0, nargmax, dim=DIM(res,0), wrap_n;
        // Compute gain/phase of kernel
        for (int n=0; n < dim; n++) {
            val = IND1(ker,n,T);
            mval = val * val;
            if (mval > mq) {
                mq = mval;
                q = val;
            }
        }
        q = 1/q;
        // The clean loop
        for (int i=0; i < maxiter; i++) {
            nscore = 0;
            mmax = -1;
            step = (T) gain * max * q;
            IND1(mdl,argmax,T) += step;
            // Take next step and compute score
            for (int n=0; n < dim; n++) {
                wrap_n = (n + argmax) % dim;
                IND1(res,wrap_n,T) -= IND1(ker,n,T) * step;
                val = IND1(res,wrap_n,T);
                mval = val * val;
                nscore += mval;
                if (mval > mmax) {
                    nargmax = wrap_n;
                    max = val;
                    mmax = mval;
                }
            }
            argmax = nargmax;
            nscore = sqrt(nscore / dim);
            if (nscore < tol) {
                // We're done
                return i;
            } else if (score > 0 && nscore > score) {
                // We've diverged: undo last step and give up
                IND1(mdl,argmax,T) -= step;
                for (int n=0; n < dim; n++) {
                    wrap_n = (n + argmax) % dim;
                    IND1(res,wrap_n,T) += IND1(ker,n,T) * step;
                }
                return -i;
            }
            score = nscore;
        }
        return maxiter;
    }

    // Does a 1d complex-valued clean
    static int clean_1d_c(PyArrayObject *res, PyArrayObject *ker, 
            PyArrayObject *mdl, double gain, int maxiter, double tol) {
        T maxr=0, maxi=0, valr, vali, stepr, stepi, qr, qi;
        T score=-1, nscore, mmax, mval, mq=0;
        int argmax=0, nargmax, dim=DIM(res,0), wrap_n;
        // Compute gain/phase of kernel
        for (int n=0; n < dim; n++) {
            valr = CIND1R(ker,n,T);
            vali = CIND1I(ker,n,T);
            mval = valr * valr + vali * vali;
            if (mval > mq) {
                mq = mval;
                qr = valr;
                qi = vali;
            }
        }
        qr /= mq;
        qi = -qi / mq;
        // The clean loop
        for (int i=0; i < maxiter; i++) {
            nscore = 0;
            mmax = -1;
            stepr = (T) gain * (maxr * qr - maxi * qi);
            stepi = (T) gain * (maxr * qi + maxi * qr);
            CIND1R(mdl,argmax,T) += stepr;
            CIND1I(mdl,argmax,T) += stepi;
            // Take next step and compute score
            for (int n=0; n < dim; n++) {
                wrap_n = (n + argmax) % dim;
                CIND1R(res,wrap_n,T) -= CIND1R(ker,n,T) * stepr - \
                                        CIND1I(ker,n,T) * stepi;
                CIND1I(res,wrap_n,T) -= CIND1R(ker,n,T) * stepi + \
                                        CIND1I(ker,n,T) * stepr;
                valr = CIND1R(res,wrap_n,T);
                vali = CIND1I(res,wrap_n,T);
                mval = valr * valr + vali * vali;
                nscore += mval;
                if (mval > mmax) {
                    nargmax = wrap_n;
                    maxr = valr;
                    maxi = vali;
                    mmax = mval;
                }
            }
            argmax = nargmax;
            nscore = sqrt(nscore / dim);
            if (nscore < tol) {
                // We're done
                return i;
            } else if (score > 0 && nscore > score) {
                // We've diverged: undo last step and give up
                CIND1R(mdl,argmax,T) -= stepr;
                CIND1I(mdl,argmax,T) -= stepi;
                for (int n=0; n < dim; n++) {
                    wrap_n = (n + argmax) % dim;
                    CIND1R(res,wrap_n,T) -= CIND1R(ker,n,T) * stepr - \
                                            CIND1I(ker,n,T) * stepi;
                    CIND1I(res,wrap_n,T) -= CIND1R(ker,n,T) * stepi + \
                                            CIND1I(ker,n,T) * stepr;
                }
                return -i;
            }
            score = nscore;
        }
        return maxiter;
    }
};  // END TEMPLATE

// Clean wrapper that handles all different data types
PyObject *clean1d(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyArrayObject *res, *ker, *mdl;
    double gain=.1, tol=.001;
    int maxiter=200, dim, rv;
    static char *kwlist[] = {"res","ker","mdl","gain","maxiter","tol",NULL};
    // Parse arguments and perform sanity check
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!O!|did", kwlist, \
            &PyArray_Type, &res, &PyArray_Type, &ker, &PyArray_Type, &mdl, 
            &gain, &maxiter, &tol)) 
        return NULL;
    CHK_ARRAY_RANK(res, 1); CHK_ARRAY_RANK(ker, 1); CHK_ARRAY_RANK(mdl, 1);
    dim = DIM(res,0);
    CHK_ARRAY_DIM(ker, 0, dim); CHK_ARRAY_DIM(mdl, 0, dim);
    if (TYPE(res) != TYPE(ker) || TYPE(res) != TYPE(mdl)) {
        PyErr_Format(PyExc_ValueError, "array types must match");
        return NULL;
    }
    Py_INCREF(res); Py_INCREF(ker); Py_INCREF(mdl);
    // Use template to implement data loops for all data types
    if (TYPE(res) == NPY_FLOAT) {
        rv = Clean<float>::clean_1d_r(res,ker,mdl,gain,maxiter,tol);
    } else if (TYPE(res) == NPY_DOUBLE) {
        rv = Clean<double>::clean_1d_r(res,ker,mdl,gain,maxiter,tol);
    } else if (TYPE(res) == NPY_LONGDOUBLE) {
        rv = Clean<long double>::clean_1d_r(res,ker,mdl,gain,maxiter,tol);
    } else if (TYPE(res) == NPY_CFLOAT) {
        rv = Clean<float>::clean_1d_c(res,ker,mdl,gain,maxiter,tol);
    } else if (TYPE(res) == NPY_CDOUBLE) {
        rv = Clean<double>::clean_1d_c(res,ker,mdl,gain,maxiter,tol);
    } else if (TYPE(res) == NPY_CLONGDOUBLE) {
        rv = Clean<long double>::clean_1d_c(res,ker,mdl,gain,maxiter,tol);
    } else {
        PyErr_Format(PyExc_ValueError, "Unsupported data type.");
        return NULL;
    }
    Py_DECREF(res); Py_DECREF(ker); Py_DECREF(mdl);
    return Py_BuildValue("(Oi)", PyArray_Return(mdl), rv);
}

// Wrap function into module
static PyMethodDef DeconvMethods[] = {
    {"clean1d", (PyCFunction)clean1d, METH_VARARGS|METH_KEYWORDS,
        "clean1d(res,ker,mdl,gain=.1,maxiter=200,tol=.001)\nPerform a 1 dimensional real deconvolution using the CLEAN algorithm.."},
    {NULL, NULL}
};

PyMODINIT_FUNC init_deconv(void) {
    (void) Py_InitModule("_deconv", DeconvMethods);
    import_array();
};
