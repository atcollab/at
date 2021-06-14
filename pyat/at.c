/*
 * This file contains the Python interface to AT, compatible with
 * Python 3 only. It provides a module 'atpass' containing one method 'atpass'.
 */
#include <stdarg.h>
#include <Python.h>
#ifdef _OPENMP
#include <string.h>
#include <omp.h>
#endif /*_OPENMP*/
#include "attypes.h"
#include <stdbool.h> 
#include <math.h>
#include <float.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

#define NUMPY_IMPORT_ARRAY_RETVAL NULL
#define NUMPY_IMPORT_ARRAY_TYPE void *

typedef PyObject atElem;

#define ATPY_PASS "trackFunction"

#if defined(PCWIN) || defined(PCWIN64) || defined(_WIN32)
#include <windows.h>
#define LIBRARYHANDLETYPE HINSTANCE
#define FREELIBFCN(libfilename) FreeLibrary((libfilename))
#define LOADLIBFCN(libfilename) LoadLibrary((libfilename))
#define GETTRACKFCN(libfilename) GetProcAddress((libfilename),ATPY_PASS)
#define SEPARATOR "\\"
#define OBJECTEXT ".pyd"
#else
#include <dlfcn.h>
#define LIBRARYHANDLETYPE void *
#define FREELIBFCN(libfilename) dlclose(libfilename)
#define LOADLIBFCN(libfilename) dlopen((libfilename),RTLD_LAZY)
#define GETTRACKFCN(libfilename) dlsym((libfilename),ATPY_PASS)
#define SEPARATOR "/"
#define OBJECTEXT ".so"
#endif

#define MOD_ERROR_VAL NULL
#define MOD_SUCCESS_VAL(val) val
#define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)

#define LIMIT_AMPLITUDE		1

/* define the general signature of a pass function */
typedef struct elem *(*track_function)(const PyObject *element,
                                      struct elem *elemptr,
                                      double *r_in,
                                      int num_particles,
                                      struct parameters *param);

static npy_uint32 num_elements = 0;
static struct elem **elemdata_list = NULL;
static PyObject **element_list = NULL;
static track_function *integrator_list = NULL;
static PyObject **pyintegrator_list = NULL;
static PyObject **kwargs_list = NULL;
static char integrator_path[300];

/* Directly copied from atpass.c */
static struct LibraryListElement {
    const char *MethodName;
    LIBRARYHANDLETYPE LibraryHandle;
    track_function FunctionHandle;
    PyObject *PyFunctionHandle;
    struct LibraryListElement *Next;
} *LibraryList = NULL;

static PyObject *print_error(int elem_number, PyObject *rout)
{
    printf("Error in tracking element %d.\n", elem_number);
    Py_XDECREF(rout);
    return NULL;
}

static PyObject *set_error(PyObject *errtype, const char *fmt, ...)
{
    char buffer[132];
    va_list ap;
    va_start(ap, fmt);
    vsprintf(buffer, fmt, ap);
    PyErr_SetString(errtype, buffer);
    va_end(ap);
    return NULL;
}
/*
 * Recursively search the list to check if the library containing
 * method_name is already loaded. If it is - return the pointer to the
 * list element. If not, return NULL.
 */
static struct LibraryListElement* SearchLibraryList(struct LibraryListElement *head, const char *method_name)
{
    if (head)
        return (strcmp(head->MethodName, method_name)==0) ? head :
            SearchLibraryList(head->Next, method_name);
    else
        return NULL;
}


static void checkiflost(double *drin, npy_uint32 np, int num_elem, int num_turn, 
        double *xnturn, double *xnelem, double *xcoord, double *xlostcoord, 
        bool *xlost, double *histbuf, npy_uint32 ihist, npy_uint32 lhist)
{
    int n, c;
    for (c=0; c<np; c++) {/* Loop over particles */
        if (!xlost[c]) {  /* No change if already marked */
           double *r6 = drin+c*6;
           for (n=0; n<6; n++) {	/* I remove the check on the sixth coordinate N.C. */
                if (!isfinite(r6[n]) || (fabs(r6[n])>LIMIT_AMPLITUDE)) {
                    int h, k=ihist;
                    xlost[c] = 1;
                    xnturn[c] = num_turn;
                    xnelem[c] = num_elem;
                    for (h=0; h<lhist; h++) {
                        if (++k >= lhist) k=0;
                        memcpy(xcoord+6*(np*h+c),histbuf+6*(np*k+c),6*sizeof(double));
                    }
                    memcpy(xlostcoord+6*c,r6,6*sizeof(double));
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

/*
 * Use Python calls to establish the location of the at integrators
 * package.
 */
static PyObject *get_integrators(void) {
    PyObject *at_module, *os_module, *fileobj, *dirname_function, *dirobj;
    at_module = PyImport_ImportModule("at.integrators");
    if (at_module == NULL) return NULL;
    fileobj = PyObject_GetAttrString(at_module, "__file__");
    Py_DECREF(at_module);
    if (fileobj == NULL) return NULL;
    os_module = PyImport_ImportModule("os.path");
    if (os_module == NULL) return NULL;
    dirname_function = PyObject_GetAttrString(os_module, "dirname");
    Py_DECREF(os_module);
    if (dirname_function == NULL) return NULL;
    dirobj = PyObject_CallFunctionObjArgs(dirname_function, fileobj, NULL);
    Py_DECREF(fileobj);
    Py_DECREF(dirname_function);
    return dirobj;
}

/*
 * Query Python for the full extension given to shared objects.
 * This is useful for Python 3, where the extension may not be trivial.
 * If none is defined, return NULL.
 */
static PyObject *get_ext_suffix(void) {
    PyObject *sysconfig_module, *get_config_var_fn, *ext_suffix;
    sysconfig_module = PyImport_ImportModule("distutils.sysconfig");
    if (sysconfig_module == NULL) return NULL;
    get_config_var_fn = PyObject_GetAttrString(sysconfig_module, "get_config_var");
    Py_DECREF(sysconfig_module);
    if (get_config_var_fn == NULL) return NULL;
    ext_suffix = PyObject_CallFunction(get_config_var_fn, "s", "EXT_SUFFIX");
    Py_DECREF(get_config_var_fn);
    return ext_suffix;
}


/*
 * Import the python module for python integrators
 * and return the function object
 */
static PyObject *GetpyFunction(const char *fn_name)
{
  char dest[300];
  strcpy(dest,"at.integrators.");
  strcat(dest,fn_name);
  PyObject *pModule;
  pModule = PyImport_ImportModule(fn_name);
  if (!pModule){
      PyErr_Clear();
      pModule = PyImport_ImportModule(dest);
  }
  if(!pModule){
      return NULL;
  }
  PyObject *pyfunction = PyObject_GetAttrString(pModule, "trackFunction");
  if ((!pyfunction) || !PyCallable_Check(pyfunction)) {
      Py_DECREF(pModule);
      if(pyfunction){
          Py_DECREF(pyfunction);
      }
      return NULL;
  }
  Py_DECREF(pModule);
  return pyfunction;
}

/*
 * Build input positional arguments for python integrators
 */
static PyObject *Buildkwargs(const atElem *ElemData)
{
  PyObject *kwargs;
  kwargs = PyDict_New();
  PyDict_SetItemString(kwargs,(char *)"elem",(PyObject *)ElemData);
  return kwargs;
}

/*
 * Build input keyword arguments for python integrators
 */
static PyObject *Buildargs(double *r_in, int num_particles)
{
  npy_intp outdims[1];
  outdims[0] = 6*num_particles;
  PyObject *rin;
  rin = PyArray_SimpleNewFromData(1, outdims, NPY_DOUBLE, r_in);
  if (!rin){
      printf("PyFuncPass: could not generate pyArray rin");
    }
  return PyTuple_Pack(1,rin);
}

/*
 * Call python integrators
 */
static PyObject *PyIntegratorPass(double *r_in, PyObject *function, PyObject *kwargs, int num_particles)
{
  PyObject *args;
  args = Buildargs(r_in, num_particles);
  PyObject_Call(function, args, kwargs);
  Py_DECREF(args);
  return kwargs;
}

/*
 * Find the correct track function by name.
 */
static struct LibraryListElement* get_track_function(const char *fn_name) {

    struct LibraryListElement *LibraryListPtr = SearchLibraryList(LibraryList, fn_name);

    if (!LibraryListPtr) {
        LIBRARYHANDLETYPE dl_handle=NULL;
        track_function fn_handle = NULL;
        char lib_file[300], buffer[200];
        PyObject *pyfunction = NULL;

        pyfunction = GetpyFunction(fn_name);

        if(!pyfunction){
            snprintf(lib_file, sizeof(lib_file), integrator_path, fn_name);
            dl_handle = LOADLIBFCN(lib_file);
            if (dl_handle) {
                fn_handle = (track_function) GETTRACKFCN(dl_handle);
            }
        }
        
        if((fn_handle==NULL) && (pyfunction==NULL)){
            snprintf(buffer, sizeof(buffer), "PassMethod %s: library, module or trackFunction not found", fn_name);
            if(dl_handle){
                FREELIBFCN(dl_handle);
            }
            if(pyfunction){
                Py_DECREF(pyfunction);
            }
            PyErr_SetString(PyExc_RuntimeError, buffer);
            return NULL;
        }

        LibraryListPtr = (struct LibraryListElement *)malloc(sizeof(struct LibraryListElement));
        LibraryListPtr->MethodName = strcpy(malloc(strlen(fn_name)+1), fn_name);
        LibraryListPtr->LibraryHandle = dl_handle;
        LibraryListPtr->FunctionHandle = fn_handle;
        LibraryListPtr->PyFunctionHandle = pyfunction;
        LibraryListPtr->Next = LibraryList;
        LibraryList = LibraryListPtr;
    }
    return LibraryListPtr;
}

/*
 * Parse the arguments to atpass, set things up, and execute.
 * Arguments:
 *  - line: sequence of elements
 *  - rin: numpy 6-vector of initial conditions
 *  - nturns: int number of turns to simulate
 *  - refpts: numpy uint32 array denoting elements at which to return state
 *  - reuse: whether to reuse the cached state of the ring
 */
static PyObject *at_atpass(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"line","rin","nturns","refpts","reuse","omp_num_threads", NULL};
    static double lattice_length = 0.0;
    int new_lattice;

    PyObject *lattice;
    PyArrayObject *rin;
    PyArrayObject *refs = NULL;
    PyObject *rout;
    double *drin, *drout;
    int num_turns;
    npy_uint32 omp_num_threads=0;
    npy_uint32 num_particles, np6;
    npy_uint32 elem_index;
    npy_uint32 *refpts = NULL;
    npy_uint32 nextref;
    unsigned int nextrefindex;
    unsigned int num_refpts;
    npy_uint32 reuse=0;
    npy_intp outdims[4];
    int turn;
    #ifdef _OPENMP
    int maxthreads;
    #endif /*_OPENMP*/
    struct parameters param;
    struct LibraryListElement *LibraryListPtr;


    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!i|O!II", kwlist, &PyList_Type, &lattice,
        &PyArray_Type, &rin, &num_turns, &PyArray_Type, &refs, &reuse, &omp_num_threads)) {
        return NULL;
    }
    if (PyArray_DIM(rin,0) != 6) {
        PyErr_SetString(PyExc_ValueError, "Numpy array is not 6D");
        return NULL;
    }
    if (PyArray_TYPE(rin) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "rin is not a double array");
        return NULL;
    }
    if ((PyArray_FLAGS(rin) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
        PyErr_SetString(PyExc_ValueError, "rin is not Fortran-aligned");
        return NULL;
    }

    num_particles = (PyArray_SIZE(rin)/6);
    np6 = num_particles*6;
    drin = PyArray_DATA(rin);
    new_lattice = (reuse == 0) ? 1 : 0;

    if (refs) {
        if (PyArray_TYPE(refs) != NPY_UINT32) {
            PyErr_SetString(PyExc_ValueError, "refpts is not a uint32 array");
            return NULL;
        }
        refpts = PyArray_DATA(refs);
        num_refpts = PyArray_SIZE(refs);
    }
    else {
        refpts = NULL;
        num_refpts = 0;
    }
    outdims[0] = 6;
    outdims[1] = num_particles;
    outdims[2] = num_refpts;
    outdims[3] = num_turns;
    rout = PyArray_EMPTY(4, outdims, NPY_DOUBLE, 1);
    drout = PyArray_DATA((PyArrayObject *)rout);


    #ifdef _OPENMP
    if ((omp_num_threads > 0) && (num_particles > OMP_PARTICLE_THRESHOLD)) {
        unsigned int nthreads = omp_get_num_procs();
        maxthreads = omp_get_max_threads();
        if (omp_num_threads < nthreads) nthreads = omp_num_threads;
        if (num_particles < nthreads) nthreads = num_particles;
        omp_set_num_threads(nthreads);
    }
    #endif /*_OPENMP*/

    if (new_lattice) {
        PyObject **element;
        track_function *integrator;
        PyObject **pyintegrator;
        for (elem_index=0; elem_index < num_elements; elem_index++) {
            free(elemdata_list[elem_index]);
            Py_XDECREF(element_list[elem_index]);        /* Release the stored elements, may be NULL if */
        }                                           /* a previous call was interrupted by an error */
        num_elements = PyList_Size(lattice);

        /* Pointer to Element structures used by the tracking function */
        free(elemdata_list);
        elemdata_list = (struct elem **)calloc(num_elements, sizeof(struct elem *));

        /* Pointer to Element list, make sure all pointers are initially NULL */
        free(element_list);
        element_list = (PyObject **)calloc(num_elements, sizeof(PyObject *));

        /* pointer to the list of C integrators */
        integrator_list = (track_function *)realloc(integrator_list, num_elements*sizeof(track_function));

        /* pointer to the list of python integrators, make sure all pointers are initially NULL */
        free(pyintegrator_list);
        pyintegrator_list = (PyObject **)calloc(num_elements, sizeof(PyObject *));

        /* pointer to the list of python integrators kwargs, make sure all pointers are initially NULL */
        free(kwargs_list);
        kwargs_list = (PyObject **)calloc(num_elements, sizeof(PyObject *));

        lattice_length = 0.0;
        element = element_list;
        integrator = integrator_list;
        pyintegrator = pyintegrator_list;
        for (elem_index = 0; elem_index < num_elements; elem_index++) {
            PyObject *el = PyList_GET_ITEM(lattice, elem_index);
            PyObject *PyPassMethod = PyObject_GetAttrString(el, "PassMethod");
            double length;
            if (!PyPassMethod) return print_error(elem_index, rout);     /* No PassMethod */
            LibraryListPtr = get_track_function(PyUnicode_AsUTF8(PyPassMethod));
            if (!LibraryListPtr) return print_error(elem_index, rout);        /* No trackFunction for the given PassMethod */
            length = PyFloat_AsDouble(PyObject_GetAttrString(el, "Length"));
            if (PyErr_Occurred()) PyErr_Clear();
            else lattice_length += length;
            *integrator++ = LibraryListPtr->FunctionHandle;
            *pyintegrator++ = LibraryListPtr->PyFunctionHandle;
            *element++ = el;
            Py_INCREF(el);                          /* Keep a reference to each element in case of reuse */
        }
        new_lattice = 0;
    }

    param.RingLength = lattice_length;
    param.T0 = lattice_length/299792458.0;
    for (turn = 0; turn < num_turns; turn++) {
        PyObject **element = element_list;
        track_function *integrator = integrator_list;
        PyObject **pyintegrator = pyintegrator_list;
        PyObject **kwargs = kwargs_list;
        struct elem **elemdata = elemdata_list;
        param.nturn = turn;
        nextrefindex = 0;
        nextref= (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
        for (elem_index = 0; elem_index < num_elements; elem_index++) {
            if (elem_index == nextref) {
                memcpy(drout, drin, np6*sizeof(double));
                drout += np6; /*  shift the location to write to in the output array */
                nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
            }
            /* the actual integrator call */
            if(*pyintegrator) {
                if(!*kwargs) *kwargs = Buildkwargs(*element);
                *kwargs = PyIntegratorPass(drin, *pyintegrator, *kwargs, num_particles);
                if (!*kwargs) return print_error(elem_index, rout);       /* trackFunction failed */
            }else {
                *elemdata = (*integrator)(*element, *elemdata, drin, num_particles, &param);
                if (!*elemdata) return print_error(elem_index, rout);       /* trackFunction failed */
            }
            element++;
            integrator++;
            pyintegrator++;
            elemdata++;
            kwargs++;
        }
        /* the last element in the ring */
        if (num_elements == nextref) {
            memcpy(drout, drin, np6*sizeof(double));
            drout += np6; /*  shift the location to write to in the output array */
        }
    }

    #ifdef _OPENMP
    if ((omp_num_threads > 0) && (num_particles > OMP_PARTICLE_THRESHOLD)) {
        omp_set_num_threads(maxthreads);
    }
    #endif /*_OPENMP*/
    return rout;
}

static PyObject *at_elempass(PyObject *self, PyObject *args)
{
    PyObject *element;
    PyArrayObject *rin;
    PyObject *PyPassMethod;
    npy_uint32 num_particles;
    track_function fn_handle;
    PyObject *pyfn_handle;
    PyObject *kwargs = NULL;
    double *drin;
    struct parameters param;
    struct elem *elem_data = NULL;
    struct LibraryListElement *LibraryListPtr;

    if (!PyArg_ParseTuple(args, "OO!", &element,  &PyArray_Type, &rin)) {
        return NULL;
    }
    if (PyArray_DIM(rin,0) != 6) {
        return set_error(PyExc_ValueError, "rin is not 6D");
    }
    if (PyArray_TYPE(rin) != NPY_DOUBLE) {
        return set_error(PyExc_ValueError, "rin is not a double array");
    }
    if ((PyArray_FLAGS(rin) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
        return set_error(PyExc_ValueError, "rin is not Fortran-aligned");
    }
    num_particles = (PyArray_SIZE(rin)/6);
    drin = PyArray_DATA(rin);

    param.RingLength = 0.0;
    param.T0 = 0.0;
    param.nturn = 0;

    PyPassMethod = PyObject_GetAttrString(element, "PassMethod");
    if (!PyPassMethod) return NULL;
    LibraryListPtr = get_track_function(PyUnicode_AsUTF8(PyPassMethod));
    fn_handle = LibraryListPtr->FunctionHandle;
    pyfn_handle = LibraryListPtr->PyFunctionHandle;
    if (fn_handle) {
        elem_data = fn_handle(element, elem_data, drin, num_particles, &param);
    }
    else {
        kwargs = Buildkwargs(element);
        kwargs = PyIntegratorPass(drin, pyfn_handle, kwargs, num_particles);
    }
    if ((!elem_data) && (!kwargs)) return NULL;
    free(elem_data);
    if(pyfn_handle) Py_DECREF(pyfn_handle);
    if(kwargs) Py_DECREF(kwargs);
    Py_INCREF(Py_None);
    return Py_None;
}
/* Boilerplate to register methods. */

static PyMethodDef AtMethods[] = {
    {"atpass",  (PyCFunction)at_atpass, METH_VARARGS | METH_KEYWORDS,
    PyDoc_STR("rout = atpass(line, rin, n_turns, refpts=[], reuse=False, omp_num_threads=0)\n\n"
              "Track input particles rin along line for nturns turns.\n"
              "Record 6D phase space at elements corresponding to refpts for each turn.\n\n"
              "line:    list of elements\n"
              "rin:     6 x n_particles Fortran-ordered numpy array.\n"
              "         On return, rin contains the final coordinates of the particles\n"
              "n_turns: number of turns to be tracked\n"
              "refpts:  numpy array of indices of elements where output is desired\n"
              "         0 means entrance of the first element\n"
              "         len(line) means end of the last element\n"
              "reuse:   if True, use previously cached description of the lattice.\n\n"
              "rout:    6 x n_particles x n_refpts x n_turns Fortran-ordered numpy array\n"
              "         of particle coordinates\n"
              )},
    {"elempass",  (PyCFunction)at_elempass, METH_VARARGS,
    PyDoc_STR("rout = elempass(element, rin)\n\n"
              "Track input particles rin through a single element.\n\n"
              "element: AT element\n"
              "rin:     6 x n_particles Fortran-ordered numpy array.\n"
              "         On return, rin contains the final coordinates of the particles\n"
             )},
   {NULL, NULL, 0, NULL}        /* Sentinel */
};

MOD_INIT(atpass)
{
    PyObject *integ_path_obj, *ext_suffix_obj;
    const char *ext_suffix, *integ_path;

    static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "at",         /* m_name */
    PyDoc_STR("Clone of atpass in Accelerator Toolbox"),      /* m_doc */
    -1,           /* m_size */
    AtMethods,    /* m_methods */
    NULL,         /* m_reload */
    NULL,         /* m_traverse */
    NULL,         /* m_clear */
    NULL,         /* m_free */
    };
    PyObject *m = PyModule_Create(&moduledef);

    if (m == NULL) return MOD_ERROR_VAL;
    import_array();

    integ_path_obj = get_integrators();
    if (integ_path_obj == NULL) return MOD_ERROR_VAL;
    ext_suffix_obj = get_ext_suffix();
    if (ext_suffix_obj == NULL) return MOD_ERROR_VAL;
    ext_suffix = (ext_suffix_obj == Py_None) ? OBJECTEXT : PyUnicode_AsUTF8(ext_suffix_obj);
    integ_path = PyUnicode_AsUTF8(integ_path_obj);
    snprintf(integrator_path, sizeof(integrator_path), "%s%s%%s%s", integ_path, SEPARATOR, ext_suffix);
    Py_DECREF(integ_path_obj);
    Py_DECREF(ext_suffix_obj);

    return MOD_SUCCESS_VAL(m);
}
