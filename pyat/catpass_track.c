#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "catpass.h"
#ifdef _OPENMP
#include <string.h>
#include <omp.h>
#endif /*_OPENMP*/
#include <math.h>
#include <float.h>

#define atPrintf(...) PySys_WriteStdout(__VA_ARGS__)

#define LIMIT_AMPLITUDE		1

static void checkiflost(double *drin, npy_uint32 np, int num_elem, int num_turn,
        int *lossTurn, int *lossElem, bool *lost, double *lossCoords)
{
    for (unsigned int c=0; c<np; c++) {/* Loop over particles */
        if (!lost[c]) {  /* No change if already marked */
            double *r6 = drin+c*6;
            for (unsigned int n=0; n<6; n++) {
                if (!isfinite(r6[n]) || ((fabs(r6[n])>LIMIT_AMPLITUDE)&&n<5)) {
                    if (lossTurn) {
                        lost[c] = 1;
                        lossTurn[c] = num_turn;
                        lossElem[c] = num_elem;
                        memcpy(lossCoords+6*c,r6,6*sizeof(double));
                    }
                    r6[0] = NAN;
                    r6[1] = 0;
                    r6[2] = 0;
                    r6[3] = 0;
                    r6[4] = 0;
                    r6[5] = 0;
                    break;
                }
            }
        }
    }
}

PyObject *line_pass(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {
        "", "", "nturns", "refpts", "turn", "keep_counter",
        "omp_num_threads", "losses", "start_elem", "end_elem", NULL
    };
    PyObject *pyLossTurn = NULL;
    PyObject *pyLossElem = NULL;
    PyObject *pyLost = NULL;
    PyObject *pyLossCoords = NULL;
    int *lossTurn = NULL;
    int *lossElem = NULL;
    bool *lost = NULL;
    double *lossCoords = NULL;
    npy_uint32 *refpts = NULL;
    npy_uint32 nextref;
    unsigned int nextrefindex;
    unsigned int num_refpts = 0;
    npy_intp outdims[4];
    npy_intp pdims[1];
    #ifdef _OPENMP
    int maxthreads;
    #endif /*_OPENMP*/

    /* Positional arguments */
    PyObject *capsule = NULL;
    PyArrayObject *rin = NULL;
    int num_turns=1;
    /* Keyword arguments */
    PyArrayObject *refs=NULL;
    int counter=0;
    int keep_counter=0;
    npy_uint32 omp_num_threads=0;
    int losses=0;
    npy_uint32 start_elem = 0;
    npy_uint32 end_elem = 0;

    at_state *module_data = (at_state *)PyModule_GetState(self);
    if (module_data == NULL) return NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!|iO!$ipIpii", kwlist,
        &PyCapsule_Type, &capsule, &PyArray_Type, &rin, &num_turns,
        &PyArray_Type, &refs, &counter,
        &keep_counter, &omp_num_threads, &losses,
        &start_elem, &end_elem)) return NULL;

    lat_state *lattice_data  = (lat_state *)PyCapsule_GetPointer(capsule, CAPSULE_NAME);
    if (lattice_data == NULL) return NULL;

    struct parameters *param = &(lattice_data->params);
    npy_uint32 num_elements = lattice_data->num_elements;

    if (!keep_counter) param->nturn = counter;

    if (PyArray_DIM(rin,0) != 6) {
        return PyErr_Format(PyExc_ValueError, "rin is not 6D");
    }
    if (PyArray_TYPE(rin) != NPY_DOUBLE) {
        return PyErr_Format(PyExc_ValueError, "rin is not a double array");
    }
    if ((PyArray_FLAGS(rin) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
        return PyErr_Format(PyExc_ValueError, "rin is not Fortran-aligned");
    }
    double *drin = PyArray_DATA(rin);
    npy_uint32 num_particles = (PyArray_SIZE(rin)/6);
    npy_uint32 np6 = num_particles*6;
    if (refs) {
        if (PyArray_TYPE(refs) != NPY_UINT32) {
            return PyErr_Format(PyExc_ValueError, "refpts is not a uint32 array");
        }
        refpts = PyArray_DATA(refs);
        num_refpts = PyArray_SIZE(refs);
    }
    pdims[0] = num_particles;
    outdims[0] = 6;
    outdims[1] = num_particles;
    outdims[2] = num_refpts;
    outdims[3] = num_turns;

    #ifdef _OPENMP
    if ((omp_num_threads > 0) && (num_particles > OMP_PARTICLE_THRESHOLD)) {
        unsigned int nthreads = omp_get_num_procs();
        maxthreads = omp_get_max_threads();
        if (omp_num_threads < nthreads) nthreads = omp_num_threads;
        if (num_particles < nthreads) nthreads = num_particles;
        omp_set_num_threads(nthreads);
    }
    #endif /*_OPENMP*/

    /* Check for partial turn tracking */
    if (end_elem==0) end_elem = num_elements;
    if (num_elements!=0) {
        if (start_elem<0 || start_elem>=num_elements) {
            char errMsg[64];
            sprintf(errMsg,"start_elem out of range, %d not in [0..%d]",start_elem,num_elements-1);
            return PyErr_Format(PyExc_ValueError, errMsg);
        }
        if (end_elem<0 || end_elem>num_elements) {
            char errMsg[64];
            sprintf(errMsg,"end_elem out of range, %d not in [0..%d]",start_elem,num_elements);
            return PyErr_Format(PyExc_ValueError, errMsg);
        }
        if (start_elem>end_elem) {
            return PyErr_Format(PyExc_ValueError, "end_elem must be greater that start_elem");
        }
        if ((end_elem!=num_elements) && (num_turns!=1)) {
            return PyErr_Format(PyExc_ValueError, "partial turn tracking not ending at the end of the ring requires nturns=1");
        }
    } else {
        start_elem = 0;
        end_elem = 0;
    }

    /* Memory allocation for outputs */
    PyObject *rout = PyArray_EMPTY(4, outdims, NPY_DOUBLE, 1);
    double *drout = PyArray_DATA((PyArrayObject *)rout);

    /* Memory allocation for loss information */
    pyLost = PyArray_ZEROS(1, pdims, NPY_BOOL, 1);
    lost = PyArray_DATA((PyArrayObject *)pyLost);
    if (losses) {
        pyLossTurn = PyArray_ZEROS(1, pdims, NPY_UINT32, 1);
        pyLossElem = PyArray_ZEROS(1, pdims, NPY_UINT32, 1);
        pyLossCoords = PyArray_ZEROS(2, outdims, NPY_DOUBLE, 1);
        lossTurn = PyArray_DATA((PyArrayObject *)pyLossTurn);
        lossElem = PyArray_DATA((PyArrayObject *)pyLossElem);
        lossCoords = PyArray_DATA((PyArrayObject *)pyLossCoords);
    }

    for (int turn = 0; turn < num_turns; turn++) {
        PyObject **element = lattice_data->element_list;
        double *elem_length = lattice_data->elemlength_list;
        track_function *integrator = lattice_data->integrator_list;
        PyObject **pyintegrator = lattice_data->pyintegrator_list;
/*      PyObject **kwargs = lattice_data->kwargs_list; */
        struct elem **elemdata = lattice_data->elemdata_list;
        double s_coord = 0.0;

        if (start_elem != 0) {
            /* Update pointers and s_coord when start_elem is not 0 */
            for (npy_uint32 elem_index = 0; elem_index < start_elem; elem_index++ ) {
                s_coord += *elem_length++;
                element++;
                integrator++;
                pyintegrator++;
                elemdata++;
            }
        }

        nextrefindex = 0;
        nextref= (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
        for (npy_uint32 elem_index = start_elem; elem_index < end_elem; elem_index++) {
            param->s_coord = s_coord;
            if (elem_index == nextref) {
                memcpy(drout, drin, np6*sizeof(double));
                drout += np6; /*  shift the location to write to in the output array */
                nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
            }
            /* the actual integrator call */
            if (*pyintegrator) {
                PyObject *res = PyObject_CallFunctionObjArgs(*pyintegrator, rin, *element, NULL);
                if (!res) return error_return(elem_index, 5, rout, pyLost, pyLossTurn, pyLossElem, pyLossCoords);       /* trackFunction failed */
                Py_DECREF(res);
            } else {
                *elemdata = (*integrator)(*element, *elemdata, drin, num_particles, param);
                if (!*elemdata) return error_return(elem_index, 5, rout, pyLost, pyLossTurn, pyLossElem, pyLossCoords); /* trackFunction failed */
            }
            checkiflost(drin, num_particles, elem_index, param->nturn, lossTurn, lossElem, lost, lossCoords);
            s_coord += *elem_length++;
            element++;
            integrator++;
            pyintegrator++;
            elemdata++;
        }
        /* the last element in the ring */
        if (num_elements == nextref) {
            memcpy(drout, drin, np6*sizeof(double));
            drout += np6; /*  shift the location to write to in the output array */
        }
        param->nturn++;
        start_elem = 0;
    }

    #ifdef _OPENMP
    if ((omp_num_threads > 0) && (num_particles > OMP_PARTICLE_THRESHOLD)) {
        omp_set_num_threads(maxthreads);
    }
    #endif /*_OPENMP*/

    if (losses) {
        PyObject *tout = PyTuple_New(2);
        PyObject *dict = PyDict_New();
        PyDict_SetItemString(dict,(char *)"islost",(PyObject *)pyLost);
        PyDict_SetItemString(dict,(char *)"turn",(PyObject *)pyLossTurn);
        PyDict_SetItemString(dict,(char *)"elem",(PyObject *)pyLossElem);
        PyDict_SetItemString(dict,(char *)"coord",(PyObject *)pyLossCoords);
        PyTuple_SetItem(tout, 0, rout);
        PyTuple_SetItem(tout, 1, dict);
        Py_DECREF(pyLost);
        Py_DECREF(pyLossTurn);
        Py_DECREF(pyLossElem);
        Py_DECREF(pyLossCoords);
        return tout;
    } else {
        Py_DECREF(pyLost);
        return rout;
    }
}
