#ifndef CATPASSS_H
#define CATPASSS_H

#ifndef NPY_MAIN
#define NO_IMPORT_ARRAY
#endif
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL PYAT_ARRAY_API
#include <numpy/ndarrayobject.h>
// #if NPY_ABI_VERSION < 0x02000000
//     #define NPY_RAVEL_AXIS 32
// #endif

#include <stdarg.h>
#include "attypes.h"
#include "atrandom.h"

#define CAPSULE_NAME "at.lattice"

/* define the general signature of a pass function */
typedef struct elem *(*track_function)(
    const PyObject *element,
    struct elem *elemptr,
    double *r_in,
    int num_particles,
    struct parameters *param
);

/* Module state */
typedef struct _at_state {
    struct LibraryListElement *LibraryList;
    double electron_mass;
    double clight;
    /* state buffers for RNGs */
    pcg32_random_t common_state;
    pcg32_random_t thread_state;
    PyObject *particle_type;
    PyObject *element_type;
    int last_turn;
} at_state;

/* Lattice state */
typedef struct _lat_state {
    struct parameters params;
    PyArrayObject *bunch_currents;
    PyArrayObject *bunch_spos;
    PyArrayObject *fillpattern;
    npy_uint32 num_elements;
    struct elem **elemdata_list;
    PyObject **element_list;
    double *elemlength_list;
    track_function *integrator_list;
    PyObject **pyintegrator_list;
/*  PyObject **kwargs_list; */
} lat_state;

PyObject *line_pass(PyObject *self, PyObject *args, PyObject *kwargs);
PyObject *reset_rng(PyObject *self, PyObject *args, PyObject *kwargs);
PyObject *common_rng(PyObject *self);
PyObject *thread_rng(PyObject *self);
int catpass_random_exec(PyObject *self);

static void *error_return(int elem_number, int count, ...)
{
    va_list ap;
    va_start(ap, count);
    for (int i = 0; i < count; i++) {
        Py_XDECREF(va_arg(ap, void *));
    }
    va_end(ap);
    return NULL;
}

#endif /* CATPASS _H */