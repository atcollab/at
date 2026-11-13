/*
 * This file contains the Python interface to AT, compatible with
 * Python 3 only. It provides a module 'atpass' containing 4 python functions:
 * atpass, elempass, isopenmp, ismpi
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#ifdef _OPENMP
#include <string.h>
#include <omp.h>
#endif /*_OPENMP*/
#ifdef MPI
#include <mpi.h>
#endif /* MPI */
#include "attypes.h"
#include <stdbool.h> 
#include <math.h>
#include <float.h>
#include <atrandom.c>

#define atPrintf(...) PySys_WriteStdout(__VA_ARGS__)

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>
#if NPY_ABI_VERSION < 0x02000000
    #define NPY_RAVEL_AXIS 32
#endif

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

#define SYSCONFIG "sysconfig"
#define LIMIT_AMPLITUDE		1
#define C0  	2.99792458e8

/* define the general signature of a pass function */
typedef struct elem *(*track_function)(const PyObject *element,
                                      struct elem *elemptr,
                                      double *r_in,
                                      int num_particles,
                                      struct parameters *param);

/* Lattice state */
typedef struct _lat_state {
    npy_uint32 num_elements;
    struct elem **elemdata_list;
    PyObject **element_list;
    double *elemlength_list;
    track_function *integrator_list;
    PyObject **pyintegrator_list;
    PyObject **kwargs_list;
    double lattice_length;
    int last_turn;
} lat_state;

/* Module state */
typedef struct _at_state {
    /* state buffers for RNGs */
    pcg32_random_t common_state;
    pcg32_random_t thread_state;
} at_state;

/* Global state */
static char integrator_path[300];
static PyObject *particle_type;
static PyObject *element_type;
static struct LibraryListElement {
    const char *MethodName;
    LIBRARYHANDLETYPE LibraryHandle;
    track_function FunctionHandle;
    PyObject *PyFunctionHandle;
    struct LibraryListElement *Next;
} *LibraryList = NULL;


static PyObject *print_error(int elem_number, PyObject *rout)
{
    Py_XDECREF(rout);
    return NULL;
}

static const char *pyprint(PyObject* pyobj) {
    PyObject *pystr = PyObject_Str(pyobj);
    const char* str = PyUnicode_AsUTF8(pystr);
    Py_XDECREF(pystr);
    return str;
}

/*
 * Recursively search the list to check if the library containing
 * method_name is already loaded. If it is - return the pointer to the
 * list element. If not, return NULL.
 */
static struct LibraryListElement *SearchLibraryList(const char *method_name)
{
    struct LibraryListElement *head = LibraryList;
    while (head) {
        if (strcmp(head->MethodName, method_name) == 0) return head;
        head = head->Next;
    }
    return head;
}


static void checkiflost(double *drin, npy_uint32 np, int num_elem, int num_turn, 
        int *xnturn, int *xnelem, bool *xlost, double *xlostcoord)
{
    unsigned int n, c;
    for (c=0; c<np; c++) {/* Loop over particles */
        if (!xlost[c]) {  /* No change if already marked */
           double *r6 = drin+c*6;
           for (n=0; n<6; n++) {
                if (!isfinite(r6[n]) || ((fabs(r6[n])>LIMIT_AMPLITUDE)&&n<5)) {
                    xlost[c] = 1;
                    xnturn[c] = num_turn;
                    xnelem[c] = num_elem;
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


static void setlost(double *drin, npy_uint32 np)
{
    unsigned int n, c;
    for (c=0; c<np; c++) {/* Loop over particles */
        double *r6 = drin+c*6;
        if (isfinite(r6[0])) {  /* No change if already marked */
           for (n=0; n<6; n++) {
                if (!isfinite(r6[n]) || ((fabs(r6[n])>LIMIT_AMPLITUDE)&&n<5)) {
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


/* Get a reference to a python object in a module
   Equivalent to "from module_name import object" */
static PyObject *get_pyobj(const char *module_name, const char *object)
{
    PyObject *pyobj = NULL;
    PyObject *module = PyImport_ImportModule(module_name);
    if (module) {
        pyobj = PyObject_GetAttrString(module, object);
        Py_DECREF(module);
    }
    return pyobj;
}

/* Get the location of the at integrators package. */
static PyObject *get_integrators(void) {
    PyObject *dirobj = NULL;
    PyObject *fileobj = get_pyobj("at.integrators", "__file__");
    if (fileobj) {
        PyObject *dirname_function = get_pyobj("os.path", "dirname");
        if (dirname_function) {
            dirobj = PyObject_CallFunctionObjArgs(dirname_function, fileobj, NULL);
            Py_DECREF(dirname_function);
        }
        Py_DECREF(fileobj);
    }
    return dirobj;
}
/*
 * Query Python for the full extension given to shared objects.
 * This is useful for Python 3, where the extension may not be trivial.
 * If none is defined, return NULL.
 */
static PyObject *get_ext_suffix(void) {
    PyObject *get_config_var_fn, *ext_suffix;
    get_config_var_fn = get_pyobj(SYSCONFIG, "get_config_var");
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
      Py_XDECREF(pyfunction);
      return NULL;
  }
  Py_DECREF(pModule);
  return pyfunction;
}

/*
 * Find the correct track function by name.
 */
static struct LibraryListElement* get_track_function(const char *fn_name) {

    struct LibraryListElement *LibraryListPtr = SearchLibraryList(fn_name);

    if (!LibraryListPtr) {
        LIBRARYHANDLETYPE dl_handle=NULL;
        track_function fn_handle = NULL;

        PyObject *pyfunction = GetpyFunction(fn_name);
        PyErr_Clear();      /* Clear any import error if there is no python integrator */

        if (!pyfunction){
            char lib_file[300];
            snprintf(lib_file, sizeof(lib_file), integrator_path, fn_name);
            dl_handle = LOADLIBFCN(lib_file);
            if (dl_handle) {
                fn_handle = (track_function) GETTRACKFCN(dl_handle);
            }
        }
        
        if (!(fn_handle || pyfunction)) {
            if (dl_handle) FREELIBFCN(dl_handle);
            PyErr_Format(PyExc_RuntimeError,
            "PassMethod %s: library, module or trackFunction not found",
            fn_name);
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

static void set_particle(PyObject *particle, struct parameters *param)
{
    if (particle != NULL) {
        PyObject *prest_energy = PyObject_GetAttrString(particle, "rest_energy");
        PyObject *pcharge = PyObject_GetAttrString(particle, "charge");
        if (prest_energy != NULL) {
            param->rest_energy = PyFloat_AsDouble(prest_energy);
            Py_DECREF(prest_energy);
        }
        if (pcharge != NULL) {
            param->charge = PyFloat_AsDouble(pcharge);
            Py_DECREF(pcharge);
        }
        PyErr_Clear();
    }
}

static void *get_lattice_data(PyObject *lattice)
{
    PyObject *capsule = PyObject_GetAttrString(lattice, "Cdata");
    if (capsule == NULL) return NULL;
    return PyCapsule_GetPointer(capsule, "at.lattice");
}

static void set_current_fillpattern(PyArrayObject *bspos, PyArrayObject *bcurrents,
                             PyArrayObject *fillpattern, struct parameters *param){ 
    if (bcurrents && bspos && fillpattern) {
        PyObject *bcurrentsum = PyArray_Sum(bcurrents, NPY_RAVEL_AXIS,
                                            PyArray_DESCR(bcurrents)->type_num,
                                            NULL);
        param->beam_current = PyFloat_AsDouble(bcurrentsum);
        Py_DECREF(bcurrentsum);
        param->nbunch = PyArray_SIZE(bspos);
        param->bunch_spos = PyArray_DATA(bspos);
        param->bunch_currents = PyArray_DATA(bcurrents);
        param->fillpattern = PyArray_DATA(fillpattern); 
    }
    else {
        param->beam_current=0.0;
        param->nbunch=1;
        param->bunch_spos = (double[1]){0.0};
        param->bunch_currents = (double[1]){0.0};
        param->fillpattern = (double[1]){1.0};
    }
}

static void *free_lattice(lat_state *latdata)
{
    struct elem **elem_data = latdata->elemdata_list;
    PyObject **elem = latdata->element_list;

    atPrintf("Freeing lattice\n");
    for (npy_uint32 elem_index=latdata->num_elements; elem_index > 0; elem_index--) {
        free(*elem_data++);
        Py_XDECREF(*elem);   /* Release the stored elements, may be NULL if */
        elem++;              /* a previous call was interrupted by an error */
    }
    free(latdata->elemdata_list);
    free(latdata->elemlength_list);
    free(latdata->element_list);
    free(latdata->integrator_list);
    free(latdata->pyintegrator_list);
    free(latdata->kwargs_list);
    free(latdata);
    return NULL;
}

static void *allocate_lattice(PyObject *lattice)
{
    PyObject **elem;
    double *elem_length;
    track_function *integrator,integr;
    PyObject **pyintegrator;
    npy_uint32 num_elements = PyList_Size(lattice);
    double length = 0.0;
    struct elem **elem_data;

    atPrintf("Allocating capsule\nnum_elements: %d\n", num_elements);
    lat_state *latdata = (lat_state *)calloc(1, sizeof(lat_state));

     latdata->num_elements = num_elements;

    /* Pointer to Element structures used by the tracking function */
    atPrintf("Allocating elem_data\n");
    elem_data = (struct elem **)calloc(num_elements, sizeof(struct elem *));
    latdata->elemdata_list = elem_data;

    /* Pointer to Element lengths */
    atPrintf("Allocating elem_length\n");
    elem_length = (double *)calloc(num_elements, sizeof(double));
    latdata->elemlength_list = elem_length;

    /* Element list, make sure all pointers are initially NULL */
    atPrintf("Allocating elements pointers\n");
    elem = (PyObject **)calloc(num_elements, sizeof(PyObject *));
    latdata->element_list = elem;

    /* List of C integrators */
    atPrintf("Allocating C integrators\n");
    integrator = (track_function *)calloc(num_elements, sizeof(track_function));
    latdata->integrator_list = integrator;

    /* List of python integrators, make sure all pointers are initially NULL */
    atPrintf("Allocating python integrators\n");
    pyintegrator = (PyObject **)calloc(num_elements, sizeof(PyObject *));
    latdata->pyintegrator_list = pyintegrator;

    /* List of python integrators kwargs, make sure all pointers are initially NULL */
    atPrintf("Allocating integrator data\n");
    latdata->kwargs_list = (PyObject **)calloc(num_elements, sizeof(PyObject *));

    for (npy_uint32 elem_index = 0; elem_index < num_elements; elem_index++) {
        struct LibraryListElement *LibraryListPtr;
        PyObject *pylength;
        PyObject *el = PyList_GET_ITEM(lattice, elem_index);
        PyObject *PyPassMethod = PyObject_GetAttrString(el, "PassMethod");
        double lg;
        if (!PyPassMethod) return free_lattice(latdata);    /* No PassMethod: AttributeError */
        LibraryListPtr = get_track_function(PyUnicode_AsUTF8(PyPassMethod));
        Py_DECREF(PyPassMethod);
        if (!LibraryListPtr)
         return free_lattice(latdata);  /* No trackFunction for the given PassMethod: RuntimeError */
        integr = LibraryListPtr->FunctionHandle;
        pylength = PyObject_GetAttrString(el, "Length");
        if (pylength) {
            lg = PyFloat_AsDouble(pylength);
            Py_XDECREF(pylength);
        }
        else {
            lg = 0.0;
            PyErr_Clear();
        }

        if (integr) {
            struct elem *truc = (*integrator)(el, NULL, NULL, 0, NULL);
            if (truc) {
                atPrintf("Ok\n");
//                *elem_data = truc;
                }
            else {
                atPrintf("Erreur\n");
                return free_lattice(latdata);
                }
        }

        length += lg;
        *integrator++ = integr;
        *pyintegrator++ = LibraryListPtr->PyFunctionHandle;
        *elem++ = el;
        elem_data++;
        *elem_length++ = lg;
        Py_INCREF(el);                   /* Keep a reference to each element in case of reuse */
    }

    return latdata;
    latdata->lattice_length = length;

    return latdata;
}

static void free_capsule(PyObject *capsule)
{
    lat_state *latdata = (lat_state *)PyCapsule_GetPointer(capsule, "at.lattice");
    free_lattice (latdata);
}

static PyObject *at_allocate(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *lattice;
    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &lattice)) return NULL;

    lat_state *latdata = allocate_lattice(lattice);
    if (latdata == NULL) {
        atPrintf("NULL returned");
        return NULL;
        return free_lattice(latdata);
    }

    return PyCapsule_New(latdata, "at.lattice", free_capsule);
}
/*
 * Parse the arguments to atpass, set things up, and execute.
 * Arguments:
 *  - line: sequence of elements
 *  - rin: numpy 6-vector of initial conditions
 *  - nturns: int number of turns to simulate
 *  - refpts: numpy uint32 array denoting elements at which to return state
 *  - reuse: whether to reuse the cached state of the ring
 *  - start_elem: start tracking at element #start_elem
 *  - end_elem: end tracking at element #end_elem (exclusive)
 */
static PyObject *at_atpass(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"line","rin","nturns","refpts","turn",
                             "energy", "particle", "keep_counter",
                             "reuse","omp_num_threads","losses",
                             "bunch_spos", "bunch_currents", "start_elem", "end_elem",
                             "harmonic_number", "fillpattern", NULL};
    PyObject *rout;
    double *drin, *drout;
    PyObject *xnturn = NULL;
    PyObject *xnelem = NULL;
    PyObject *xlost = NULL;
    PyObject *xlostcoord = NULL;
    int *ixnturn = NULL;
    int *ixnelem = NULL;
    bool *bxlost = NULL;
    double *dxlostcoord = NULL;
    npy_uint32 num_particles, np6;
    npy_uint32 *refpts = NULL;
    npy_uint32 nextref;
    unsigned int nextrefindex;
    unsigned int num_refpts;
    npy_intp outdims[4];
    npy_intp pdims[1];
    npy_intp lxdims[2];
    int turn;
    #ifdef _OPENMP
    int maxthreads;
    #endif /*_OPENMP*/
    struct parameters param;
    npy_uint32 num_elements;

    /* Positional arguments */
    PyObject *capsule = NULL;
    PyArrayObject *rin;
    int num_turns;
    /* Keyword arguments */
    PyArrayObject *refs=NULL;
    int counter=0;
    double energy=0.0;
    PyObject *particle=NULL;
    int keep_counter=0;
    npy_uint32 omp_num_threads=0;
    int losses=0;
    PyArrayObject *bspos=NULL;
    PyArrayObject *bcurrents=NULL;
    npy_uint32 start_elem = 0;
    npy_uint32 end_elem = 0;
    int harmonic_number = 1;
    PyArrayObject *fillpattern=NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!i|O!$idO!pIpO!O!iiiO!", kwlist,
        &PyCapsule_Type, &capsule, &PyArray_Type, &rin, &num_turns,
        &PyArray_Type, &refs, &counter,
        &energy, particle_type, &particle,
        &keep_counter, &omp_num_threads, &losses,
        &PyArray_Type, &bspos, &PyArray_Type, &bcurrents,
        &start_elem, &end_elem, &harmonic_number,
        &PyArray_Type, &fillpattern)) {
        return NULL;
    }
    if (PyArray_DIM(rin,0) != 6) {
        return PyErr_Format(PyExc_ValueError, "rin is not 6D");
    }
    if (PyArray_TYPE(rin) != NPY_DOUBLE) {
        return PyErr_Format(PyExc_ValueError, "rin is not a double array");
    }
    if ((PyArray_FLAGS(rin) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
        return PyErr_Format(PyExc_ValueError, "rin is not Fortran-aligned");
    }
    num_particles = (PyArray_SIZE(rin)/6);

    at_state *state = (at_state *)PyModule_GetState(self);
    if (state == NULL) return NULL;
    lat_state *latdata  = (lat_state *)PyCapsule_GetPointer(capsule, "at.lattice");
    if (latdata == NULL) return NULL;

    param.common_rng=&(state->common_state);
    param.thread_rng=&(state->thread_state);
    param.energy=energy;
    param.rest_energy=0.0;
    param.charge=-1.0;
    param.num_turns=num_turns;
    param.bdiff = NULL; 
    param.harmonic_number = harmonic_number;   
    if (keep_counter)
        param.nturn = latdata->last_turn;
    else
        param.nturn = counter;

    /* Particle energy */
    set_particle(particle, &param);
    param.RingLength = latdata->lattice_length;
    if (param.rest_energy == 0.0) {
        param.T0 = param.RingLength/C0;
    } else {
        /*  Only relativistic tracking is implemented: keep everything consistent
        double gamma0 = param.energy/param.rest_energy;
        double betagamma0 = sqrt(gamma0*gamma0 - 1.0);
        double beta0 = betagamma0/gamma0;
        param.T0 = param.RingLength/(beta0*C0);
        */
        param.T0 = param.RingLength/C0;
    }

    set_current_fillpattern(bspos, bcurrents, fillpattern, &param);

    #ifdef _OPENMP
    if ((omp_num_threads > 0) && (num_particles > OMP_PARTICLE_THRESHOLD)) {
        unsigned int nthreads = omp_get_num_procs();
        maxthreads = omp_get_max_threads();
        if (omp_num_threads < nthreads) nthreads = omp_num_threads;
        if (num_particles < nthreads) nthreads = num_particles;
        omp_set_num_threads(nthreads);
    }
    #endif /*_OPENMP*/

    num_elements = latdata->num_elements;

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
    np6 = num_particles*6;
    drin = PyArray_DATA(rin);

    if (refs) {
        if (PyArray_TYPE(refs) != NPY_UINT32) {
            return PyErr_Format(PyExc_ValueError, "refpts is not a uint32 array");
        }
        refpts = PyArray_DATA(refs);
        num_refpts = PyArray_SIZE(refs);
    } else {
        refpts = NULL;
        num_refpts = 0;
    }
    outdims[0] = 6;
    outdims[1] = num_particles;
    outdims[2] = num_refpts;
    outdims[3] = num_turns;
    rout = PyArray_EMPTY(4, outdims, NPY_DOUBLE, 1);
    drout = PyArray_DATA((PyArrayObject *)rout);

    if(losses){
        pdims[0]= num_particles;
        lxdims[0]= 6;
        lxdims[1]= num_particles;
        xnturn = PyArray_EMPTY(1, pdims, NPY_UINT32, 1);
        xnelem = PyArray_EMPTY(1, pdims, NPY_UINT32, 1);
        xlost = PyArray_EMPTY(1, pdims, NPY_BOOL, 1);
        xlostcoord = PyArray_EMPTY(2, lxdims, NPY_DOUBLE, 1);
        ixnturn = PyArray_DATA((PyArrayObject *)xnturn);
        ixnelem = PyArray_DATA((PyArrayObject *)xnelem);
        bxlost = PyArray_DATA((PyArrayObject *)xlost);
        dxlostcoord = PyArray_DATA((PyArrayObject *)xlostcoord);
        unsigned int i;
        static double r0[6];
        for(i=0;i<num_particles;i++){
            bxlost[i]=0;
            ixnturn[i]=0;
            ixnelem[i]=0;
            memcpy(dxlostcoord+6*i,r0,6*sizeof(double));
        }
    }

    for (turn = 0; turn < num_turns; turn++) {

        PyObject **element = latdata->element_list;
        double *elem_length = latdata->elemlength_list;
        track_function *integrator = latdata->integrator_list;
        PyObject **pyintegrator = latdata->pyintegrator_list;
        PyObject **kwargs = latdata->kwargs_list;
        struct elem **elemdata = latdata->elemdata_list;
        double s_coord = 0.0;

        if (start_elem != 0) {
            /* Update pointers and s_coord when start_elem is not 0 */
            for (npy_uint32 elem_index = 0;elem_index<start_elem;elem_index++ ) {
                s_coord += *elem_length++;
                element++;
                integrator++;
                pyintegrator++;
                elemdata++;
            }
        }

        /*PySys_WriteStdout("turn: %i\n", param.nturn);*/
        nextrefindex = 0;
        nextref= (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
        for (npy_uint32 elem_index = start_elem; elem_index < end_elem; elem_index++) {
            param.s_coord = s_coord;
            if (elem_index == nextref) {
                memcpy(drout, drin, np6*sizeof(double));
                drout += np6; /*  shift the location to write to in the output array */
                nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
            }
            /* the actual integrator call */
            if (*pyintegrator) {
                PyObject *res = PyObject_CallFunctionObjArgs(*pyintegrator, rin, *element, NULL);
                if (!res) return print_error(elem_index, rout);       /* trackFunction failed */
                Py_DECREF(res);
            } else {
                *elemdata = (*integrator)(*element, *elemdata, drin, num_particles, &param);
                if (!*elemdata) return print_error(elem_index, rout);       /* trackFunction failed */
            }
            if (losses) {
                checkiflost(drin, num_particles, elem_index, param.nturn, ixnturn, ixnelem, bxlost, dxlostcoord);
            } else {
                setlost(drin, num_particles);
            }
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
        param.nturn++;
        start_elem = 0;
    }
    latdata->last_turn = param.nturn;  /* Store turn number in a static variable */

    #ifdef _OPENMP
    if ((omp_num_threads > 0) && (num_particles > OMP_PARTICLE_THRESHOLD)) {
        omp_set_num_threads(maxthreads);
    }
    #endif /*_OPENMP*/

    if (losses) {
        PyObject *tout = PyTuple_New(2);
        PyObject *dict = PyDict_New();
        PyDict_SetItemString(dict,(char *)"islost",(PyObject *)xlost);
        PyDict_SetItemString(dict,(char *)"turn",(PyObject *)xnturn); 
        PyDict_SetItemString(dict,(char *)"elem",(PyObject *)xnelem);
        PyDict_SetItemString(dict,(char *)"coord",(PyObject *)xlostcoord);
        PyTuple_SetItem(tout, 0, rout);  
        PyTuple_SetItem(tout, 1, dict);
        Py_DECREF(xlost);
        Py_DECREF(xnturn);
        Py_DECREF(xnelem);
        Py_DECREF(xlostcoord);
        return tout;          
    } else {
        return rout;
    }
}

static PyObject *at_elempass(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"element", "rin",
                             "energy", "particle", NULL};
    PyObject *PyPassMethod;
    npy_uint32 num_particles;
    track_function integrator;
    PyObject *pyintegrator;
    double *drin;
    struct parameters param;
    struct LibraryListElement *LibraryListPtr;

    /* Positional arguments */
    PyObject *element;
    PyArrayObject *rin;
    /* Keyword arguments */
    double energy = 0.0;
    PyObject *particle = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!|$dO!", kwlist,
        element_type, &element,  &PyArray_Type, &rin,
        &energy, particle_type, &particle)) {
        return NULL;
    }
    if (PyArray_DIM(rin,0) != 6) {
        return PyErr_Format(PyExc_ValueError, "rin is not 6D");
    }
    if (PyArray_TYPE(rin) != NPY_DOUBLE) {
        return PyErr_Format(PyExc_ValueError, "rin is not a double array");
    }
    if ((PyArray_FLAGS(rin) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
        return PyErr_Format(PyExc_ValueError, "rin is not Fortran-aligned");
    }

    at_state *state = (at_state*)PyModule_GetState(self);
    if (state == NULL) return NULL;

    param.common_rng=&(state->common_state);
    param.thread_rng=&(state->thread_state);
    param.energy=energy;
    param.rest_energy=0.0;
    param.charge=-1.0;
    param.num_turns = 1;
    param.bdiff=NULL;
    param.harmonic_number = 1;
    param.nturn = 0;

    set_current_fillpattern(NULL, NULL, NULL, &param);
    set_particle(particle, &param);

    num_particles = (PyArray_SIZE(rin)/6);
    drin = PyArray_DATA(rin);

    param.RingLength = 0.0;
    param.T0 = 0.0;

    PyPassMethod = PyObject_GetAttrString(element, "PassMethod");
    if (!PyPassMethod) return NULL;
    LibraryListPtr = get_track_function(PyUnicode_AsUTF8(PyPassMethod));
    Py_DECREF(PyPassMethod);
    if (!LibraryListPtr) return NULL;
    integrator = LibraryListPtr->FunctionHandle;
    pyintegrator = LibraryListPtr->PyFunctionHandle;
    if (pyintegrator) {
        PyObject *res = PyObject_CallFunctionObjArgs(pyintegrator, rin, element, NULL);
        if (!res) return NULL;
        Py_DECREF(res);
    } else {
        struct elem *elem_data = integrator(element, NULL, drin, num_particles, &param);
        if (!elem_data) return NULL;
        free(elem_data);
    }
    Py_INCREF(rin);
    return (PyObject *) rin;
}

static PyObject *at_diffmatrix(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"element", "rin", "energy", "particle", NULL};
    PyObject *element;
    double energy = 0.0;
    PyObject *particle = NULL;
    PyArrayObject *rin;
    npy_uint32 num_particles;
    track_function integrator;
    PyObject *pyintegrator;
    PyObject *PyPassMethod;
    const double *drin;
    double orb[6];
    struct parameters param;
    struct LibraryListElement *LibraryListPtr;
    npy_intp dims[2] = {6, 6};
    PyObject *pbdiff = PyArray_ZEROS(2, dims, NPY_DOUBLE ,1);
    double *bdiff = (double *)PyArray_DATA((PyArrayObject *)pbdiff);

    param.nturn = 0;
    param.energy=0.0;
    param.rest_energy=0.0;
    param.charge=-1.0;
    param.bdiff=bdiff;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!|$dO!", kwlist,
        element_type, &element,  &PyArray_Type, &rin,
        &energy, particle_type, &particle)) {
        return NULL;
    }
    if (PyArray_DIM(rin,0) != 6) {
        return PyErr_Format(PyExc_ValueError, "rin is not 6D");
    }
    if (PyArray_TYPE(rin) != NPY_DOUBLE) {
        return PyErr_Format(PyExc_ValueError, "rin is not a double array");
    }
    if ((PyArray_FLAGS(rin) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
        return PyErr_Format(PyExc_ValueError, "rin is not Fortran-aligned");
    }
    num_particles = (PyArray_SIZE(rin)/6);
    if (num_particles != 1) {
        return PyErr_Format(PyExc_ValueError, "Number of particles should be 1");
    }

    param.energy = energy;
    set_particle(particle, &param);

    drin = (const double *)PyArray_DATA(rin);
    for (int i=0; i<6; i++) orb[i] = drin[i];

    param.RingLength = 0.0;
    param.T0 = 0.0;
    param.beam_current=0.0;
    param.nbunch=1;
    param.bunch_spos = (double[1]){0.0};
    param.bunch_currents = (double[1]){0.0};
    for (int i = 0; i < 36; i++) bdiff[i] = 0.0;

    PyPassMethod = PyObject_GetAttrString(element, "PassMethod");
    if (!PyPassMethod) return NULL;
    LibraryListPtr = get_track_function(PyUnicode_AsUTF8(PyPassMethod));
    Py_DECREF(PyPassMethod);

    if (!LibraryListPtr) return NULL;
    integrator = LibraryListPtr->FunctionHandle;
    pyintegrator = LibraryListPtr->PyFunctionHandle;
    if (pyintegrator) {
        PyObject *res = PyObject_CallFunctionObjArgs(pyintegrator, rin, element, NULL);
        if (!res) return NULL;
        Py_DECREF(res);
    } else {
        struct elem *elem_data = integrator(element, NULL, orb, 1, &param);
        if (!elem_data) return NULL;
        free(elem_data);
    }
    return pbdiff;
}

static PyObject *reset_rng(PyObject *self, PyObject *args, PyObject *kwargs)
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

    at_state *state = (at_state*)PyModule_GetState(self);
    if (state == NULL) return NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|$KK", kwlist,
        &rank, &seed)) {
        return NULL;
    }
    pcg32_srandom_r(&(state->common_state), seed, AT_RNG_INC);
    pcg32_srandom_r(&(state->thread_state), seed, rank);
    Py_RETURN_NONE;
}

static PyObject *common_rng(PyObject *self)
{
    at_state *state = (at_state*)PyModule_GetState(self);
    if (state == NULL) return NULL;

    double drand = atrandd_r(&(state->common_state));
    return Py_BuildValue("d", drand);
}

static PyObject *thread_rng(PyObject *self)
{
    at_state *state = (at_state*)PyModule_GetState(self);
    if (state == NULL) return NULL;

    double drand = atrandd_r(&(state->thread_state));
    return Py_BuildValue("d", drand);
}


static int atpass_exec(PyObject *self)
{
#ifdef MPI
    int irank;
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);
    uint64_t rank = irank;
#else
    uint64_t rank = 0;
#endif /* MPI */

    at_state *state = (at_state *)PyModule_GetState(self);
    if (state == NULL) return -1;
    pcg32_srandom_r(&(state->common_state), AT_RNG_STATE, AT_RNG_INC);
    pcg32_srandom_r(&(state->thread_state), AT_RNG_STATE, rank);

    atPrintf("Done exec_atpass\n");
    return 0;
}

/* Method table */

static PyMethodDef AtMethods[] = {
    {"atpass",  (PyCFunction)at_atpass, METH_VARARGS | METH_KEYWORDS,
    PyDoc_STR("atpass(line: Sequence[Element], r_in, n_turns: int, refpts: Uint32_refs = [], "
              "reuse: Optional[bool] = False, omp_num_threads: Optional[int] = 0)\n\n"
              "Track input particles r_in along line for nturns turns.\n"
              "Record 6D phase space at elements corresponding to refpts for each turn.\n\n"
              "Parameters:\n"
              "    line:    list of elements\n"
              "    rin:     6 x n_particles Fortran-ordered numpy array.\n"
              "      On return, rin contains the final coordinates of the particles\n"
              "    n_turns: number of turns to be tracked\n"
              "    refpts:  numpy array of indices of elements where output is desired\n"
              "       0 means entrance of the first element\n"
              "       len(line) means end of the last element\n"
              "    energy:  nominal energy [eV]\n"
              "    particle (Optional[Particle]):  circulating particle\n"
              "    reuse:   if True, use previously cached description of the lattice.\n"
              "    omp_num_threads: number of OpenMP threads (default 0: automatic)\n"
              "    losses:  if True, process losses\n"
              "    start_elem: Start tracking at the specified element. The\n"
              "       behavior is different from :py:func:`..lattice.rotate`, refpts are\n"
              "       still indexed from the begininig of the ring and s_coord starts at\n"
              "       the beginning of the ring. Default 0.\n"
              "    end_elem: End tracking at the specified element (exclusive).\n"
              "       When end_elem is different from the number of lattice elements,\n"
              "       nturns=1 is required. Default number of elements.\n\n"
              "Returns:\n"
              "    rout:    6 x n_particles x n_refpts x n_turns Fortran-ordered numpy array\n"
              "         of particle coordinates\n\n"
              ":meta private:"
              )},
    {"elempass",  (PyCFunction)at_elempass, METH_VARARGS | METH_KEYWORDS,
    PyDoc_STR("elempass(element, r_in)\n\n"
              "Track input particles r_in through a single element.\n\n"
              "Parameters:\n"
              "    element (Element):   AT element\n"
              "    rin:     6 x n_particles Fortran-ordered numpy array.\n"
              "      On return, rin contains the final coordinates of the particles\n"
              "    energy (float):      nominal energy [eV]\n"
              "    particle (Optional[Particle]):  circulating particle\n\n"
              ":meta private:"
            )},
    {"diffusion_matrix",  (PyCFunction)at_diffmatrix, METH_VARARGS | METH_KEYWORDS,
    PyDoc_STR("diffusion_matrix(element, r_in)\n\n"
              "Track a particles r_in through a single element,\n"
              "computes the diffusion matrix\n\n"
              "Parameters:\n"
              "    element (Element):   AT element\n"
              "    rin:     6 x n_particles Fortran-ordered numpy array.\n"
              "      On return, rin contains the final coordinates of the particles\n"
              "    energy (float):      nominal energy [eV]\n"
              "    particle (Optional[Particle]):  circulating particle\n\n"
              ":meta private:"
            )},
    {"reset_rng",  (PyCFunction)reset_rng, METH_VARARGS | METH_KEYWORDS,
    PyDoc_STR("reset_rng(*, rank=0, seed=None)\n\n"
              "Reset the *common* and *thread* random generators.\n\n"
              "The seed is applied unchanged to the \"common\" generator, and modified in a\n"
              "thread-specific way to the \"thread\" generator\n\n"
              "Parameters:\n"
              "    rank (int):    thread identifier (for MPI and python multiprocessing)\n"
              "    seed (int):    single seed for both generators. Default: initial seed\n"
            )},
    {"allocate_lattice",  (PyCFunction)at_allocate, METH_VARARGS,
    PyDoc_STR("allocate_lattice(lattice)\n\n"
              "Allocate memory for a lattice.\n\n"
              "Parameters:\n"
              "    lattice (Lattice):   AT lattice\n\n"
            )},
    {"common_rng",  (PyCFunction)common_rng, METH_NOARGS,
    PyDoc_STR("common_rng()\n\n"
              "Return a double from the *common* generator .\n"
             )},
    {"thread_rng",  (PyCFunction)thread_rng, METH_NOARGS,
    PyDoc_STR("thread_rng()\n\n"
              "Return a double from the *thread* generator .\n"
             )},
   {NULL, NULL, 0, NULL}        /* Sentinel */
};

static PyModuleDef_Slot atpass_slots[] = {
    {Py_mod_exec, atpass_exec},
#ifdef Py_GIL_DISABLED
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "atpass",         /* m_name */
    PyDoc_STR("Clone of atpass in Accelerator Toolbox"),      /* m_doc */
    sizeof(at_state), /* m_size */
    AtMethods,        /* m_methods */
    atpass_slots,     /* m_slots */
    NULL,             /* m_traverse */
    NULL,             /* m_clear */
    NULL,             /* m_free */
};

PyMODINIT_FUNC PyInit_atpass(void)
{
    PyObject *integ_path_obj;
    PyObject *m = PyModuleDef_Init(&moduledef);

    if (m == NULL) return NULL;
    import_array();

    /* Build path for loading Python integrators */
    integ_path_obj = get_integrators();
    if (integ_path_obj) {
        const char *integ_path = PyUnicode_AsUTF8(integ_path_obj);
        PyObject *ext_suffix_obj = get_ext_suffix();
        Py_DECREF(integ_path_obj);
        if (ext_suffix_obj) {
            const char *ext_suffix = (ext_suffix_obj == Py_None) ? OBJECTEXT : PyUnicode_AsUTF8(ext_suffix_obj);
            Py_DECREF(ext_suffix_obj);
            snprintf(integrator_path, sizeof(integrator_path), "%s%s%%s%s", integ_path, SEPARATOR, ext_suffix);
        }
        else {
            return NULL;
        }
    }
    else {
        return NULL;
    }

    /* get Particle type */
    particle_type = get_pyobj("at.lattice", "Particle");
    if (particle_type == NULL) return NULL;

    /* get Element type */
    element_type = get_pyobj("at.lattice", "Element");
    if (element_type == NULL) return NULL;

    atPrintf("Done init_atpass\n");

    return m;
}
