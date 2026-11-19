#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NO_IMPORT_ARRAY
#include "catpass.h"
#include "atrandom.c"
#ifdef MPI
#include <mpi.h>
#endif /* MPI */

PyObject *reset_rng(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"rank", "seed", NULL};
    uint64_t seed = AT_RNG_STATE;
#ifdef MPI
    int irank;
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    uint64_t rank = irank;
#else
    uint64_t rank = 0;
#endif /* MPI */

    at_state *module_data = (at_state*)PyModule_GetState(self);
    if (module_data == NULL) return NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|$KK", kwlist,
        &rank, &seed)) {
        return NULL;
    }
    pcg32_srandom_r(&(module_data->common_state), seed, AT_RNG_INC);
    pcg32_srandom_r(&(module_data->thread_state), seed, rank);
    Py_RETURN_NONE;
}

PyObject *common_rng(PyObject *self)
{
    at_state *module_data = (at_state*)PyModule_GetState(self);
    if (module_data == NULL) return NULL;

    double drand = atrandd_r(&(module_data->common_state));
    return Py_BuildValue("d", drand);
}

PyObject *thread_rng(PyObject *self)
{
    at_state *module_data = (at_state*)PyModule_GetState(self);
    if (module_data == NULL) return NULL;

    double drand = atrandd_r(&(module_data->thread_state));
    return Py_BuildValue("d", drand);
}

int catpass_random_exec(PyObject *self)
{
#ifdef MPI
    int irank;
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    uint64_t rank = irank;
#else
    uint64_t rank = 0;
#endif /* MPI */

    at_state *module_data = (at_state *)PyModule_GetState(self);
    if (module_data == NULL) return -1;

    pcg32_srandom_r(&(module_data->common_state), AT_RNG_STATE, AT_RNG_INC);
    pcg32_srandom_r(&(module_data->thread_state), AT_RNG_STATE, rank);
    return 0;
}
