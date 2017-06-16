
/* Cephes module version 1.5
 *  This module defines the functions in the cephes and amos libraries as
 *   Numerical python ufunc objects so that they can operate on arbitrary 
 *   NumPy arrays with broadcasting and typecasting rules implemented.
 *  
 *  Copyright 1999  Travis E. Oliphant
 * Revisions 2002 (added functions from cdflib)
 * Revisions 2008 Aaron Parsons: stripped out fortran routines
 */

#include "Python.h"
#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h" 
#include "ufunc_extras.h"
#include "abstract.h"
#include "cephes.h"
#include "c_misc/misc.h"

/* Defined in mtherr in the cephes library */
extern int scipy_special_print_error_messages;
 
#include "cephes_doc.h"


static PyUFuncGenericFunction cephes1_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes1rc_functions[] = { NULL, NULL, NULL, NULL};
static PyUFuncGenericFunction cephes1_2_functions[] = { NULL, NULL, NULL, NULL,};
static PyUFuncGenericFunction cephes1_2c_functions[] = { NULL, NULL,};
static PyUFuncGenericFunction cephes1c_4_functions[] = { NULL, NULL, NULL, NULL };
static PyUFuncGenericFunction cephes1cp_4_functions[] = { NULL, NULL, NULL, NULL};
static PyUFuncGenericFunction cephes1cpb_4_functions[] = { NULL, NULL,};
static PyUFuncGenericFunction cephes2_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes2_2_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes2_4_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes2a_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes2c_functions[] = { NULL, NULL, NULL, NULL };
static PyUFuncGenericFunction cephes2cp_functions[] = { NULL, NULL, NULL, NULL, };
static PyUFuncGenericFunction cephes2cpp_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes3_functions[] = { NULL, NULL, NULL, NULL};
static PyUFuncGenericFunction cephes3a_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes3_2_functions[] = { NULL, NULL,};
static PyUFuncGenericFunction cephes4_functions[] = { NULL, NULL, NULL, NULL,};
static PyUFuncGenericFunction cephes4a_2_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes4_2_functions[] = { NULL, NULL, };
static PyUFuncGenericFunction cephes5_2_functions[] = { NULL, NULL, };

static PyUFuncGenericFunction cephes1c_functions[] = { NULL, NULL, };


static void * ellpj_data[] = { (void *)ellpj, (void *)ellpj,};

static void * expn_data[] = { (void *)expn, (void *)expn, };
static void * kn_data[] = { (void *)kn, (void *)kn, };

static void * pdtrc_data[] = { (void *)pdtrc, (void *)pdtrc, };
static void * pdtr_data[] = { (void *)pdtr, (void *)pdtr, };
static void * pdtri_data[] = { (void *)pdtri, (void *)pdtri, };

static void * shichi_data[] = { (void *)shichi, (void *)shichi, };
static void * sici_data[] = { (void *)sici, (void *)sici, };


static void * yn_data[] = { (void *)yn, (void *)yn, };
static void * smirnov_data[] = { (void *)smirnov, (void *)smirnov, };
static void * smirnovi_data[] = { (void *)smirnovi, (void *)smirnovi, };

static void * bdtrc_data[] = { (void *)bdtrc, (void *)bdtrc, };
static void * bdtr_data[] = { (void *)bdtr, (void *)bdtr, };
static void * bdtri_data[] = { (void *)bdtri, (void *)bdtri, };
static void * btdtr_data[] = { (void *)btdtr, (void *)btdtr, };
static void * btdtri_data[] = { (void *)incbi, (void *)incbi, };

static void * fdtrc_data[] = { (void *)fdtrc, (void *)fdtrc, };
static void * fdtr_data[] = { (void *)fdtr, (void *)fdtr, };
static void * fdtri_data[] = { (void *)fdtri, (void *)fdtri, };

static void * gdtrc_data[] = { (void *)gdtrc, (void *)gdtrc, };
static void * gdtr_data[] = { (void *)gdtr, (void *)gdtr, };
/*
static void * gdtri_data[] = { (void *)gdtri, (void *)gdtri, };
*/
static void * hyp2f0_data[] = { (void *)hyp2f0, (void *)hyp2f0, };
static void * threef0_data[] = { (void *)threef0, (void *)threef0, };
static void * onef2_data[] = { (void *)onef2, (void *)onef2, };

static void * incbet_data[] = { (void *)incbet, (void *)incbet, };
static void * incbi_data[] = { (void *)incbi, (void *)incbi, };

static void * nbdtrc_data[] = { (void *)nbdtrc, (void *)nbdtrc, };
static void * nbdtr_data[] = { (void *)nbdtr, (void *)nbdtr, };
static void * nbdtri_data[] = { (void *)nbdtri, (void *)nbdtri, };

static void * beta_data[] = { (void *)beta, (void *)beta, };
static void * lbeta_data[] = { (void *)lbeta, (void *)lbeta, };
static void * cbrt_data[] = { (void *)cbrt, (void *)cbrt, };
static void * chdtrc_data[] = { (void *)chdtrc, (void *)chdtrc, };
static void * chdtr_data[] = { (void *)chdtr, (void *)chdtr, };
static void * chdtri_data[] = { (void *)chdtri, (void *)chdtri, };
static void * dawsn_data[] = {  (void *)dawsn, (void *)dawsn, };
static void * ellie_data[] = { (void *)ellie, (void *)ellie, };
static void * ellik_data[] = { (void *)ellik, (void *)ellik, };
static void * ellpe_data[] = { (void *)ellpe, (void *)ellpe, };
static void * ellpk_data[] = { (void *)ellpk, (void *)ellpk, };
static void * exp10_data[] = { (void *)exp10, (void *)exp10, };
static void * exp2_data[] = { (void *)exp2, (void *)exp2, };
static void * i0_data[] = { (void *)i0, (void *)i0, };
static void * i0e_data[] = { (void *)i0e, (void *)i0e, };
static void * i1_data[] = { (void *)i1, (void *)i1, };
static void * i1e_data[] = { (void *)i1e, (void *)i1e, };
static void * igamc_data[] = { (void *)igamc, (void *)igamc, };
static void * igam_data[] = { (void *)igam, (void *)igam, };
static void * igami_data[] = { (void *)igami, (void *)igami, };

static void * j0_data[] = { (void *)j0,  (void *)j0,  };
static void * y0_data[] = { (void *)y0, (void *)y0, };
static void * j1_data[] = { (void *)j1,  (void *)j1,  };
static void * y1_data[] = { (void *)y1, (void *)y1, };

static void * k0_data[] = { (void *)k0, (void *)k0, };
static void * k0e_data[] = { (void *)k0e, (void *)k0e, };
static void * k1_data[] = { (void *)k1, (void *)k1, };
static void * k1e_data[] = { (void *)k1e, (void *)k1e, };

static void * ndtr_data[] = { (void *)ndtr, (void *)ndtr, };
static void * erfc_data[] = { (void *)erfc, (void *)erfc, };
static void * ndtri_data[] = { (void *)ndtri, (void *)ndtri, };

static void * round_data[] = { (void *)round, (void *)round, };
static void * sindg_data[] = { (void *)sindg, (void *)sindg, };
static void * cosdg_data[] = { (void *)cosdg, (void *)cosdg, };
static void * radian_data[] = { (void *)radian, (void *)radian, };
static void * tandg_data[] = { (void *)tandg, (void *)tandg, };
static void * cotdg_data[] = { (void *)cotdg, (void *)cotdg, };
static void * log1p_data[] = { (void *)log1p, (void *)log1p, };
static void * expm1_data[] = { (void *)expm1, (void *)expm1, };
static void * cosm1_data[] = { (void *)cosm1, (void *)cosm1, };

static void * spence_data[] = { (void *)spence, (void *)spence, };
/* static void * struve_data[] = { (void *)struve, (void *)struve, };*/


static void * zeta_data[] = { (void *)zeta, (void *)zeta, };
static void * zetac_data[] = { (void *)zetac, (void *)zetac, };

static void * kolmogorov_data[] = { (void *)kolmogorov, (void *)kolmogorov, };
static void * kolmogi_data[] = { (void *)kolmogi, (void *)kolmogi, };


static void * besselpoly_data[] = {(void *)besselpoly, (void *)besselpoly,};


static char cephes_6_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE,};




static char cephes_4_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT, PyArray_FLOAT, PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, PyArray_DOUBLE,};


static char cephes_3_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_FLOAT,   PyArray_DOUBLE,  PyArray_DOUBLE, PyArray_DOUBLE, };
static char cephes_2_types[] = { PyArray_FLOAT,  PyArray_FLOAT,  PyArray_DOUBLE,  PyArray_DOUBLE,  };


/* Some functions needed from ufunc object, so that Py_complex's aren't being returned 
between code possibly compiled with different compilers.
*/

typedef Py_complex ComplexUnaryFunc(Py_complex x);

static void cephes_F_F_As_D_D(char **args, intp *dimensions, intp *steps, void *func) {
    int i; Py_complex x;
    char *ip1=args[0], *op=args[1];
    for(i=0; i<*dimensions; i++, ip1+=steps[0], op+=steps[1]) {
	x.real = ((float *)ip1)[0]; x.imag = ((float *)ip1)[1];
	x = ((ComplexUnaryFunc *)func)(x);
	((float *)op)[0] = (float)x.real;
	((float *)op)[1] = (float)x.imag;
    }
}

static void cephes_D_D(char **args, intp *dimensions, intp *steps, void *func) {
    int i; Py_complex x;
    char *ip1=args[0], *op=args[1];
    for(i=0; i<*dimensions; i++, ip1+=steps[0], op+=steps[1]) {
	x.real = ((double *)ip1)[0]; x.imag = ((double *)ip1)[1];
	x = ((ComplexUnaryFunc *)func)(x);
	((double *)op)[0] = x.real;
	((double *)op)[1] = x.imag;
    }
}

static void Cephes_InitOperators(PyObject *dictionary) {
	PyObject *f;

        cephes1_functions[0] = PyUFunc_f_f_As_d_d;
        cephes1_functions[1] = PyUFunc_d_d;
        cephes1c_functions[0] = cephes_F_F_As_D_D;
	cephes1c_functions[1] = cephes_D_D;
        cephes1rc_functions[0] = PyUFunc_f_f_As_d_d;
        cephes1rc_functions[1] = PyUFunc_d_d;
        cephes1rc_functions[2] = cephes_F_F_As_D_D;
	cephes1rc_functions[3] = cephes_D_D;
        cephes1_2_functions[0] = PyUFunc_f_ff_As_d_dd;
        cephes1_2_functions[1] = PyUFunc_d_dd;
        cephes1_2_functions[2] = PyUFunc_F_FF_As_D_DD;
        cephes1_2_functions[3] = PyUFunc_D_DD;
        cephes1_2c_functions[0] = PyUFunc_f_FF_As_d_DD;
        cephes1_2c_functions[1] = PyUFunc_d_DD;
        cephes1c_4_functions[0] = PyUFunc_f_ffff_As_d_dddd;
        cephes1c_4_functions[1] = PyUFunc_d_dddd;
        cephes1c_4_functions[2] = PyUFunc_F_FFFF_As_D_DDDD;
        cephes1c_4_functions[3] = PyUFunc_D_DDDD;
        cephes1cp_4_functions[0] = PyUFunc_f_ffff_As_D_DDDD;
        cephes1cp_4_functions[1] = PyUFunc_d_dddd_As_D_DDDD;
        cephes1cp_4_functions[2] = PyUFunc_F_FFFF_As_D_DDDD;
        cephes1cp_4_functions[3] = PyUFunc_D_DDDD;
        cephes1cpb_4_functions[0] = PyUFunc_f_FFFF_As_d_DDDD;
        cephes1cpb_4_functions[1] = PyUFunc_d_DDDD;
        cephes2_functions[0] = PyUFunc_ff_f_As_dd_d;
        cephes2_functions[1] = PyUFunc_dd_d;
        cephes2_2_functions[0] = PyUFunc_ff_ff_As_dd_dd;
        cephes2_2_functions[1] = PyUFunc_dd_dd;
        cephes2a_functions[0] = PyUFunc_ff_f_As_id_d;
        cephes2a_functions[1] = PyUFunc_dd_d_As_id_d;
        cephes2c_functions[0] = PyUFunc_ff_f_As_dd_d;
        cephes2c_functions[1] = PyUFunc_dd_d;
        cephes2c_functions[2] = PyUFunc_fF_F_As_dD_D;
        cephes2c_functions[3] = PyUFunc_dD_D;
        cephes2cp_functions[0] = PyUFunc_ff_f_As_dD_D;
        cephes2cp_functions[1] = PyUFunc_dd_d_As_dD_D;
        cephes2cp_functions[2] = PyUFunc_fF_F_As_dD_D;
        cephes2cp_functions[3] = PyUFunc_dD_D;
        cephes2cpp_functions[0] = PyUFunc_fF_F_As_dD_D;
        cephes2cpp_functions[1] = PyUFunc_dD_D;
        cephes2_4_functions[0] = PyUFunc_ff_ffff_As_dd_dddd;
        cephes2_4_functions[1] = PyUFunc_dd_dddd;
        cephes3_functions[0] = PyUFunc_fff_f_As_ddd_d;
        cephes3_functions[1] = PyUFunc_ddd_d;
        cephes3_functions[2] = PyUFunc_ffF_F_As_ddD_D;
        cephes3_functions[3] = PyUFunc_ddD_D;
        cephes3a_functions[0] = PyUFunc_fff_f_As_iid_d;
        cephes3a_functions[1] = PyUFunc_ddd_d_As_iid_d;
        cephes3_2_functions[0] = PyUFunc_fff_ff_As_ddd_dd;
        cephes3_2_functions[1] = PyUFunc_ddd_dd;
        cephes4_functions[0] = PyUFunc_ffff_f_As_dddd_d;
        cephes4_functions[1] = PyUFunc_dddd_d;
        cephes4_functions[2] = PyUFunc_fffF_F_As_dddD_D;
        cephes4_functions[3] = PyUFunc_dddD_D;
        cephes4_2_functions[0] = PyUFunc_ffff_ff_As_dddd_dd;
        cephes4_2_functions[1] = PyUFunc_dddd_dd;
        cephes4a_2_functions[0] = PyUFunc_ffff_ff_As_dddi_dd;
        cephes4a_2_functions[1] = PyUFunc_dddd_dd_As_dddi_dd;
        cephes5_2_functions[0] = PyUFunc_fffff_ff_As_ddddd_dd;
        cephes5_2_functions[1] = PyUFunc_ddddd_dd;
	
	/* Create function objects for each function call and insert
	   them in the dictionary */
	f = PyUFunc_FromFuncAndData(cephes3a_functions, bdtrc_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "bdtrc", bdtrc_doc, 0);
	PyDict_SetItemString(dictionary, "bdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, bdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "bdtr", bdtr_doc, 0);
	PyDict_SetItemString(dictionary, "bdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, bdtri_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "bdtri", bdtri_doc, 0);
	PyDict_SetItemString(dictionary, "bdtri", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, btdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "btdtr", btdtr_doc, 0);
	PyDict_SetItemString(dictionary, "btdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, btdtri_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "btdtri", btdtri_doc, 0);
	PyDict_SetItemString(dictionary, "btdtri", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, fdtrc_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "fdtrc", fdtrc_doc, 0);
        PyDict_SetItemString(dictionary, "fdtrc", f);
        Py_DECREF(f);
        f = PyUFunc_FromFuncAndData(cephes3_functions, fdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "fdtr", fdtr_doc, 0);
        PyDict_SetItemString(dictionary, "fdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, fdtri_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "fdtri", fdtri_doc, 0);
	PyDict_SetItemString(dictionary, "fdtri", f);
 	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, gdtrc_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "gdtrc", gdtrc_doc, 0);
	PyDict_SetItemString(dictionary, "gdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, gdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "gdtr", gdtr_doc, 0);
	PyDict_SetItemString(dictionary, "gdtr", f);
	Py_DECREF(f);



	f = PyUFunc_FromFuncAndData(cephes4a_2_functions, hyp2f0_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "hyp2f0", hyp2f0_doc, 0);
	PyDict_SetItemString(dictionary, "hyp2f0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_2_functions, onef2_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "hyp1f2", hyp1f2_doc, 0);
	PyDict_SetItemString(dictionary, "hyp1f2", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes4_2_functions, threef0_data, cephes_6_types, 2, 4, 2, PyUFunc_None, "hyp3f0", hyp3f0_doc, 0);
	PyDict_SetItemString(dictionary, "hyp3f0", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3_functions, incbet_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "betainc", betainc_doc, 0);
	PyDict_SetItemString(dictionary, "betainc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, incbi_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "betaincinv", betaincinv_doc, 0);
	PyDict_SetItemString(dictionary, "betaincinv", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes3a_functions, nbdtrc_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nbdtrc", nbdtrc_doc, 0);
	PyDict_SetItemString(dictionary, "nbdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, nbdtr_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nbdtr", nbdtr_doc, 0);
	PyDict_SetItemString(dictionary, "nbdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3a_functions, nbdtri_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "nbdtri", nbdtri_doc, 0);
	PyDict_SetItemString(dictionary, "nbdtri", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes2_functions, beta_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "beta", beta_doc, 0);
	PyDict_SetItemString(dictionary, "beta", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, lbeta_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "betaln", betaln_doc, 0);
	PyDict_SetItemString(dictionary, "betaln", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, cbrt_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "cbrt", cbrt_doc, 0);
	PyDict_SetItemString(dictionary, "cbrt", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, chdtrc_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "chdtrc", chdtrc_doc, 0);
	PyDict_SetItemString(dictionary, "chdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, chdtr_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "chdtr", chdtr_doc, 0);
	PyDict_SetItemString(dictionary, "chdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, chdtri_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "chdtri", chdtri_doc, 0);
	PyDict_SetItemString(dictionary, "chdtri", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, dawsn_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "dawsn", dawsn_doc, 0);
	PyDict_SetItemString(dictionary, "dawsn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, ellie_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "ellipeinc", ellipeinc_doc, 0);
	PyDict_SetItemString(dictionary, "ellipeinc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, ellik_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "ellipkinc", ellipkinc_doc, 0);
	PyDict_SetItemString(dictionary, "ellipkinc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, ellpe_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ellipe", ellipe_doc, 0);
	PyDict_SetItemString(dictionary, "ellipe", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, ellpk_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ellipk", ellipk_doc, 0);
	PyDict_SetItemString(dictionary, "ellipk", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, exp10_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "exp10", exp10_doc, 0);
	PyDict_SetItemString(dictionary, "exp10", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, exp2_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "exp2", exp2_doc, 0);
	PyDict_SetItemString(dictionary, "exp2", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, i0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "i0", i0_doc, 0);
	PyDict_SetItemString(dictionary, "i0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, i0e_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "i0e", i0e_doc, 0);
	PyDict_SetItemString(dictionary, "i0e", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, i1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "i1", i1_doc, 0);
	PyDict_SetItemString(dictionary, "i1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, i1e_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "i1e", i1e_doc, 0);
	PyDict_SetItemString(dictionary, "i1e", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2_functions, igamc_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "gammaincc", gammaincc_doc, 0);
	PyDict_SetItemString(dictionary, "gammaincc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, igam_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "gammainc", gammainc_doc, 0);
	PyDict_SetItemString(dictionary, "gammainc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2_functions, igami_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "gammainccinv", gammainccinv_doc, 0);
	PyDict_SetItemString(dictionary, "gammainccinv", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2_4_functions, ellpj_data, cephes_6_types, 2, 2, 4, PyUFunc_None, "ellipj", ellipj_doc, 0);
	PyDict_SetItemString(dictionary, "ellipj", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes2a_functions, expn_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "expn", expn_doc, 0);
	PyDict_SetItemString(dictionary, "expn", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes2a_functions, kn_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "kn", kn_doc, 0);
	PyDict_SetItemString(dictionary, "kn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, pdtrc_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "pdtrc", pdtrc_doc, 0);
	PyDict_SetItemString(dictionary, "pdtrc", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, pdtr_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "pdtr", pdtr_doc, 0);
	PyDict_SetItemString(dictionary, "pdtr", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, pdtri_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "pdtri", pdtri_doc, 0);
	PyDict_SetItemString(dictionary, "pdtri", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, yn_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "yn", yn_doc, 0);
	PyDict_SetItemString(dictionary, "yn", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, smirnov_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "smirnov", smirnov_doc, 0);
	PyDict_SetItemString(dictionary, "smirnov", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes2a_functions, smirnovi_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "smirnovi", smirnovi_doc, 0);
	PyDict_SetItemString(dictionary, "smirnovi", f);
	Py_DECREF(f);



	f = PyUFunc_FromFuncAndData(cephes1_2_functions, shichi_data, cephes_3_types, 2, 1, 2, PyUFunc_None, "shichi", shichi_doc, 0);
	PyDict_SetItemString(dictionary, "shichi", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_2_functions, sici_data, cephes_3_types, 2, 1, 2, PyUFunc_None, "sici", sici_doc, 0);
	PyDict_SetItemString(dictionary, "sici", f);
	Py_DECREF(f);




	f = PyUFunc_FromFuncAndData(cephes1_functions, j0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "j0", j0_doc, 0);
	PyDict_SetItemString(dictionary, "j0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, y0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "y0", y0_doc, 0);
	PyDict_SetItemString(dictionary, "y0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, j1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "j1", j1_doc, 0);
	PyDict_SetItemString(dictionary, "j1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, y1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "y1", y1_doc, 0);
	PyDict_SetItemString(dictionary, "y1", f);
	Py_DECREF(f);



	f = PyUFunc_FromFuncAndData(cephes1_functions, k0_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "k0", k0_doc, 0);
	PyDict_SetItemString(dictionary, "k0", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, k0e_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "k0e", k0e_doc, 0);
	PyDict_SetItemString(dictionary, "k0e", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, k1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "k1", k1_doc, 0);
	PyDict_SetItemString(dictionary, "k1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, k1e_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "k1e", k1e_doc, 0);
	PyDict_SetItemString(dictionary, "k1e", f);
	Py_DECREF(f);



	f = PyUFunc_FromFuncAndData(cephes1_functions, ndtr_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ndtr", ndtr_doc, 0);
	PyDict_SetItemString(dictionary, "ndtr", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, erfc_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "erfc", erfc_doc, 0);
	PyDict_SetItemString(dictionary, "erfc", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, ndtri_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "ndtri", ndtri_doc, 0);
	PyDict_SetItemString(dictionary, "ndtri", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, round_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "round", round_doc, 0);
	PyDict_SetItemString(dictionary, "round", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, sindg_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "sindg", sindg_doc, 0);
	PyDict_SetItemString(dictionary, "sindg", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, cosdg_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "cosdg", cosdg_doc, 0);
	PyDict_SetItemString(dictionary, "cosdg", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes3_functions, radian_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "radian", radian_doc, 0);
	PyDict_SetItemString(dictionary, "radian", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, tandg_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "tandg", tandg_doc, 0);
	PyDict_SetItemString(dictionary, "tandg", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, cotdg_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "cotdg", cotdg_doc, 0);
	PyDict_SetItemString(dictionary, "cotdg", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, log1p_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "log1p", log1p_doc, 0);
	PyDict_SetItemString(dictionary, "log1p", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, expm1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "expm1", expm1_doc, 0);
	PyDict_SetItemString(dictionary, "expm1", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, cosm1_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "cosm1", cosm1_doc, 0);
	PyDict_SetItemString(dictionary, "cosm1", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, spence_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "spence", spence_doc, 0);
	PyDict_SetItemString(dictionary, "spence", f);
	Py_DECREF(f);
	f = PyUFunc_FromFuncAndData(cephes1_functions, zetac_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "zetac", zetac_doc, 0);
	PyDict_SetItemString(dictionary, "zetac", f);
	Py_DECREF(f);




	f = PyUFunc_FromFuncAndData(cephes2_functions, zeta_data, cephes_3_types, 2, 2, 1, PyUFunc_None, "zeta", zeta_doc, 0);
	PyDict_SetItemString(dictionary, "zeta", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, kolmogorov_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "kolmogorov", kolmogorov_doc, 0);
	PyDict_SetItemString(dictionary, "kolmogorov", f);
	Py_DECREF(f);

	f = PyUFunc_FromFuncAndData(cephes1_functions, kolmogi_data, cephes_2_types, 2, 1, 1, PyUFunc_None, "kolmogi", kolmogi_doc, 0);
	PyDict_SetItemString(dictionary, "kolmogi", f);
	Py_DECREF(f);


	f = PyUFunc_FromFuncAndData(cephes3_functions, besselpoly_data, cephes_4_types, 2, 3, 1, PyUFunc_None, "besselpoly", besselpoly_doc, 0);
	PyDict_SetItemString(dictionary, "besselpoly", f);
	Py_DECREF(f);

        










}

static char errprint_doc[] = \
"errprint({flag}) sets the error printing flag for special functions\n" \
"    (from the cephesmodule). The output is the previous state.\n" \
"    With errprint(0) no error messages are shown;\n" \
"    the default is errprint(1).\n" \
"    If no argument is given the current state of\n" \
"    the flag is returned and no change occurs.\n";


static PyObject *errprint_func(PyObject *self, PyObject *args)
{
  int inflag = -37;
  int oldflag = 0;
  if (!PyArg_ParseTuple ( args, "|i;cephes.errprint", &inflag)) return NULL;

  oldflag = scipy_special_print_error_messages;  
  if (inflag != -37) {
    scipy_special_print_error_messages = (inflag != 0);
  }
  return PyInt_FromLong((long) oldflag);
}

  
static struct PyMethodDef methods[] = {
  {"errprint", errprint_func, METH_VARARGS, errprint_doc},
  {NULL,		NULL, 0}		/* sentinel */
};


PyMODINIT_FUNC init_cephes(void) {
  PyObject *m, *d, *s;
  
  /* Create the module and add the functions */
  m = Py_InitModule("_cephes", methods); 

  /* Import the ufunc objects */
  import_array();
  import_ufunc();

  /* Add some symbolic constants to the module */
  d = PyModule_GetDict(m);

  s = PyString_FromString("2.0");
  PyDict_SetItemString(d, "__version__", s);
  Py_DECREF(s);

  /* Add scipy_special_print_error_message global variable */
  /*  No, instead acessible through errprint */

  /* Load the cephes operators into the array module's namespace */
  Cephes_InitOperators(d); 
  
  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module _cephes");
}

