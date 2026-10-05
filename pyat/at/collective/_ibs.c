/*
 * Python interface to the intra-beam scattering growth rates of atibslib.c,
 * also used by the IBSPass pass method.
 */

#include "atcommon.h"
#include "atibslib.c"

#define MODULE_NAME _ibs
#define MODULE_DESCR "Intra-beam scattering growth rates"

static double ibs_lambda[IBS_NQUAD];
static double ibs_weight[IBS_NQUAD];

static PyObject *growth_rates(PyObject *self, PyObject *args)
{
    PyObject *pyoptics, *pyrates;
    PyArrayObject *optics;
    double energy, mass, charge, beta;
    struct ibs_beam b;
    npy_intp dims[1] = {3};

    if (!PyArg_ParseTuple(args, "Oddddddddd", &pyoptics, &energy, &mass,
                          &charge, &beta, &b.npart, &b.emitx, &b.emity,
                          &b.sigma_e, &b.bunch_length)) {
        return NULL;
    }
    optics = (PyArrayObject *)PyArray_FROM_OTF(pyoptics, NPY_DOUBLE, NPY_ARRAY_FARRAY_RO);
    if (optics == NULL) return NULL;
    if (PyArray_NDIM(optics) != 2 || PyArray_DIM(optics, 0) != IBS_NOPTICS) {
        Py_DECREF(optics);
        PyErr_SetString(PyExc_ValueError, "optics must be a (9, N) array");
        return NULL;
    }
    pyrates = PyArray_ZEROS(1, dims, NPY_DOUBLE, 0);
    ibs_growth_rates((int)PyArray_DIM(optics, 1), PyArray_DATA(optics),
                     ibs_lambda, ibs_weight, energy, mass, charge, beta, &b,
                     PyArray_DATA((PyArrayObject *)pyrates));
    Py_DECREF(optics);
    return pyrates;
}

static PyMethodDef IbsMethods[] = {
    {"growth_rates", (PyCFunction)growth_rates, METH_VARARGS,
    PyDoc_STR(
    "growth_rates(optics, energy, mass, charge, beta, npart, emitx, emity,\n"
    "             sigma_e, bunch_length)\n\n"
    "IBS emittance growth rates [1/s], Bjorken-Mtingwa model\n\n"
    "Args:\n"
    "    optics:       (9, N) array: integration weights normalised to 1,\n"
    "                  beta_x, beta_y, alpha_x, alpha_y, D_x, D'_x, D_y, D'_y\n"
    "    energy:       Energy [eV]\n"
    "    mass:         Rest energy [eV]\n"
    "    charge:       Particle charge [e]\n"
    "    beta:         Relativistic beta\n"
    "    npart:        Number of particles in the bunch\n"
    "    emitx:        Horizontal emittance [m]\n"
    "    emity:        Vertical emittance [m]\n"
    "    sigma_e:      Relative momentum spread\n"
    "    bunch_length: RMS bunch length [m]\n\n"
    "Returns:\n"
    "    rates:        (3,) horizontal, vertical and longitudinal emittance\n"
    "                  growth rates [1/s]\n"
    )},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC MOD_INIT(MODULE_NAME)
{
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        STR(MODULE_NAME),           /* m_name */
        PyDoc_STR(MODULE_DESCR),    /* m_doc */
        -1,                         /* m_size */
        IbsMethods,                 /* m_methods */
        NULL,                       /* m_reload */
        NULL,                       /* m_traverse */
        NULL,                       /* m_clear */
        NULL,                       /* m_free */
    };
    PyObject *m = PyModule_Create(&moduledef);
    if (m == NULL) return MOD_ERROR_VAL;
    import_array();
    ibs_quadrature(ibs_lambda, ibs_weight);
    return MOD_SUCCESS_VAL(m);
}
