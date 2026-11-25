/*
 * This file contains the Python interface to AT, compatible with
 * Python 3 only. It provides a module 'catpass' containing 4 python functions:
 * line_store, line_pass, elem_pass, diffusion_matrix
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_MAIN
#include "catpass.h"
#include <math.h>
#include <float.h>

#define atPrintf(...) PySys_WriteStdout(__VA_ARGS__)

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

#define GET_ATTRIBUTE(value, convert, object, attrname) \
{\
    PyObject *attribute = PyObject_GetAttrString(object, attrname);\
    if (attribute != NULL) {\
        *value = convert(attribute);\
        Py_DECREF(attribute);\
    }\
    else {\
        PyErr_Clear();\
    }\
}

/* Passmethod descriptor */
struct LibraryListElement {
    const char *MethodName;
    LIBRARYHANDLETYPE LibraryHandle;
    track_function FunctionHandle;  /* C passmethod */
    PyObject *PyFunctionHandle;     /* python passmethod */
    struct LibraryListElement *Next;
};

/* Global state */
static char integrator_path[300];

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
static struct LibraryListElement *SearchLibraryList(struct LibraryListElement *LibraryList, const char *method_name)
{
    struct LibraryListElement *head = LibraryList;
    while (head) {
        if (strcmp(head->MethodName, method_name) == 0) return head;
        head = head->Next;
    }
    return head;
}

/* Get a new reference to a python object in a module
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
    PyObject *pModule;
    pModule = PyImport_ImportModule(fn_name);
    if (!pModule){
        PyErr_Clear();
        char dest[300];
        snprintf(dest, sizeof(dest), "at.integrators.%s", fn_name);
        pModule = PyImport_ImportModule(dest);
    }
    if(!pModule) return error_return(0, 0);
    PyObject *pyfunction = PyObject_GetAttrString(pModule, "trackFunction");
    if ((!pyfunction) || !PyCallable_Check(pyfunction)) return error_return(0, 2, pModule, pyfunction);
     Py_DECREF(pModule);
    return pyfunction;
}

/*
 * Find the correct track function by name.
 */
static struct LibraryListElement* get_track_function(at_state *module_data, PyObject *element)
{
    PyObject *PyPassMethod = PyObject_GetAttrString(element, "PassMethod");
    if (!PyPassMethod) return error_return(0, 0);     /* No PassMethod: AttributeError */
    const char *fn_name = PyUnicode_AsUTF8(PyPassMethod);
    if (fn_name == NULL) return error_return(0, 1, PyPassMethod);
    struct LibraryListElement *LibraryListPtr = SearchLibraryList(module_data->LibraryList, fn_name);

    if (!LibraryListPtr) {
        LIBRARYHANDLETYPE dl_handle=NULL;
        track_function fn_handle = NULL;

        /* Look first for a python integrator */
        PyObject *pyfunction = GetpyFunction(fn_name);
        PyErr_Clear();      /* Clear any import error if there is no python integrator */

        /* Then look for a C integrator */
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
            return error_return(0, 1, PyPassMethod);
        }

        LibraryListPtr = (struct LibraryListElement *)malloc(sizeof(struct LibraryListElement));
        LibraryListPtr->MethodName = strcpy(malloc(strlen(fn_name)+1), fn_name);
        LibraryListPtr->LibraryHandle = dl_handle;
        LibraryListPtr->FunctionHandle = fn_handle;
        LibraryListPtr->PyFunctionHandle = pyfunction;
        LibraryListPtr->Next = module_data->LibraryList;
        module_data->LibraryList = LibraryListPtr;
    }
    Py_DECREF(PyPassMethod);
    return LibraryListPtr;
}

/* Fill the param structure with default values */
static void default_parameters(struct parameters *param, at_state *module_data)
{
    param->nturn = 0;
    param->num_turns = 0;
    param->RingLength = 0.0;
    param->T0 =0.0;
    param->s_coord = 0.0;
    param->energy = 1.0e9;
    param->rest_energy = module_data->electron_mass;
    param->charge = -1.0;
    param->gamma0 = 1.0E-9 * param->energy / __E0;
    param->beta0 = 1.0;
    param->beam_current = 0.0;
    param->nbunch = 1;
    param->bunch_spos = (double[1]){0.0};
    param->bunch_currents = (double[1]){0.0};
    param->common_rng=&(module_data->common_state);
    param->thread_rng=&(module_data->thread_state);
    param->bdiff = NULL;
    param->fillpattern = (double[1]){1.0};
    param->harmonic_number = 1;
}

/* Fill the param structure with data form the Lattice object */
static void lattice_parameters(
    struct parameters *param,
    at_state *module_data,
    PyObject *lattice,
    PyObject *particle,
)
{
     if (particle != NULL)     /* Try argument first */
        Py_INCREF(particle);
    else if (lattice != NULL) /* otherwise take from lattice */
        particle =  PyObject_GetAttrString(lattice, "particle");

    if (particle != NULL) {
        GET_ATTRIBUTE(&(param->rest_energy), PyFloat_AsDouble, particle, "rest_energy");
        GET_ATTRIBUTE(&param->charge, PyFloat_AsDouble, particle, "charge");
        GET_ATTRIBUTE(&param->energy, PyFloat_AsDouble, particle, "energy");
        Py_DECREF(particle);
    }
    else {
        PyErr_Clear();
    }

    if (lattice != NULL) {
       GET_ATTRIBUTE(&param->energy, PyFloat_AsDouble, lattice, "energy");
    }
}

static void energy_parameters(struct parameters *param, double arg_energy)
{
    if (arg_energy > 0.0) {
        param->energy = arg_energy;
    }

    if (param->rest_energy == 0.0) {
        param->T0 = param->RingLength/module_data->clight;
        param->gamma0 = 1.0E-9 * param->energy / __E0;
        param->beta0 = 1.0;
    }
    else {
        double gamma0 = param->energy / param->rest_energy;
        double betagamma0 = sqrt(gamma0*gamma0 - 1.0);
        double beta0 = betagamma0/gamma0;
        param->gamma0 = gamma0;
        param->beta0 = beta0;
        param->T0 = param->RingLength/(beta0*module_data->clight);
    }
}

/* Fill the param structure with filling pattern */
static void set_current_fillpattern(struct parameters *param, lat_state *lattice_data, PyObject *lattice, PyArrayObject *bspos)
{
    PyArrayObject *bcurrents = (PyArrayObject *)PyObject_GetAttrString(lattice, "bunch_currents");
    lattice_data->bunch_currents = bcurrents;
    if (bcurrents != NULL) {
        PyObject *bcurrentsum = PyArray_Sum(bcurrents, NPY_RAVEL_AXIS, NPY_NOTYPE, NULL);
        param->beam_current = PyFloat_AsDouble(bcurrentsum);
        param->bunch_currents = PyArray_DATA(bcurrents);
        param->nbunch = PyArray_SIZE(bcurrents);
        Py_DECREF(bcurrentsum);
    }
    else {
        PyErr_Clear();
    }

    PyArrayObject *fillpattern = (PyArrayObject *)PyObject_GetAttrString(lattice, "fillpattern");
    lattice_data->fillpattern = fillpattern;
    if (fillpattern != NULL) param->fillpattern = PyArray_DATA(fillpattern);
    else PyErr_Clear();

    if (bspos == NULL)
        bspos = (PyArrayObject *)PyObject_GetAttrString(lattice, "bunch_spos");
    else
        Py_INCREF(bspos);

    lattice_data->bunch_spos = bspos;
    if (bspos != NULL) param->bunch_spos = PyArray_DATA(bspos);
    else PyErr_Clear();

    GET_ATTRIBUTE(&(param->harmonic_number), PyLong_AsLong, lattice, "harmonic_number");
}

static void *free_lattice(lat_state *lattice_data)
{
    struct elem **elem_data = lattice_data->elemdata_list;
    PyObject **elem = lattice_data->element_list;

//  atPrintf("Freeing lattice\n");
    for (npy_uint32 elem_index=lattice_data->num_elements; elem_index > 0; elem_index--) {
        free(*elem_data++);
        Py_XDECREF(*elem++);   /* Release the stored elements */
    }
    Py_XDECREF(lattice_data->bunch_currents);
    Py_XDECREF(lattice_data->bunch_spos);
    Py_XDECREF(lattice_data->fillpattern);
    free(lattice_data->elemdata_list);
    free(lattice_data->elemlength_list);
    free(lattice_data->element_list);
    free(lattice_data->integrator_list);
    free(lattice_data->pyintegrator_list);
/*  free(lattice_data->kwargs_list); */
    free(lattice_data);
    return NULL;
}

static void *allocate_lattice(
    at_state *module_data,
    PyObject *lattice,
    double energy,
    PyObject *particle,
    PyArrayObject *bspos
)
{
    PyObject **elem;
    double *elem_length;
    track_function *integrator,integr;
    PyObject **pyintegrator;
    npy_uint32 num_elements = PyList_Size(lattice);
    double length = 0.0;
    struct elem **elem_data;

    lat_state *lattice_data = (lat_state *)calloc(1, sizeof(lat_state));
    struct parameters *param = &(lattice_data->params);
    default_parameters(param, module_data);
    lattice_parameters(param, module_data, lattice, particle);
    set_current_fillpattern(param, lattice_data, lattice, bspos);

    lattice_data->num_elements = num_elements;

    /* Pointer to Element structures used by the tracking function */
    elem_data = (struct elem **)calloc(num_elements, sizeof(struct elem *));
    lattice_data->elemdata_list = elem_data;

    /* Pointer to Element lengths */
    elem_length = (double *)calloc(num_elements, sizeof(double));
    lattice_data->elemlength_list = elem_length;

    /* Element list, make sure all pointers are initially NULL */
    elem = (PyObject **)calloc(num_elements, sizeof(PyObject *));
    lattice_data->element_list = elem;

    /* List of C integrators */
    integrator = (track_function *)calloc(num_elements, sizeof(track_function));
    lattice_data->integrator_list = integrator;

    /* List of python integrators, make sure all pointers are initially NULL */
    pyintegrator = (PyObject **)calloc(num_elements, sizeof(PyObject *));
    lattice_data->pyintegrator_list = pyintegrator;

    /* List of python integrators kwargs, make sure all pointers are initially NULL */
/*  lattice_data->kwargs_list = (PyObject **)calloc(num_elements, sizeof(PyObject *));*/

    for (npy_uint32 elem_index = 0; elem_index < num_elements; elem_index++) {
        double lg = 0.0;
        struct LibraryListElement *track_function;
        PyObject *el = PyList_GET_ITEM(lattice, elem_index);
        track_function = get_track_function(module_data, el);
        if (!track_function) return free_lattice(lattice_data);  /* No trackFunction for the given PassMethod: RuntimeError */
        integr = track_function->FunctionHandle;
        GET_ATTRIBUTE(&lg, PyFloat_AsDouble, el, "Length")

        if (integr) {
            *elem_data = (*integr)(el, NULL, NULL, 0, param);
            if (!*elem_data) {
                return free_lattice(lattice_data);
            }
         }

        length += lg;
        *integrator++ = integr;
        *pyintegrator++ = track_function->PyFunctionHandle;
        *elem++ = el;
        elem_data++;
        *elem_length++ = lg;
        Py_INCREF(el);                   /* Keep a reference to each element in case of reuse */
    }

    param->RingLength = length;
    param->nturn = 0;
    energy_parameters(param, energy);

    return lattice_data;
}

static void free_capsule(PyObject *capsule)
{
    lat_state *lattice_data = (lat_state *)PyCapsule_GetPointer(capsule, "at.lattice");
    free_lattice (lattice_data);
}

static PyObject *line_store(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"","energy", "particle", "bunch_spos", NULL};

    /* Positional arguments */
    PyObject *lattice;

    /* Keyword arguments*/
    double energy=-1.0;
    PyObject *particle=NULL;
    PyArrayObject *bspos=NULL;

    at_state *module_data = (at_state *)PyModule_GetState(self);

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!|$dO!O!", kwlist,
        &PyList_Type, &lattice,
        &energy, module_data->particle_type, &particle,
        &PyArray_Type, &bspos)) return NULL;

    lat_state *lattice_data = allocate_lattice(module_data, lattice, energy, particle, bspos);
    if (lattice_data == NULL) return NULL;
    return PyCapsule_New(lattice_data, CAPSULE_NAME, free_capsule);
}

static PyObject *elem_pass(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"", "", "energy", "particle", NULL};
    track_function integrator;
    PyObject *pyintegrator;
    struct parameters param;
    struct LibraryListElement *track_function;

    /* Positional arguments */
    PyObject *element;
    PyArrayObject *rin;
    /* Keyword arguments */
    double energy = 0.0;
    PyObject *particle = NULL;

    at_state *module_data = (at_state*)PyModule_GetState(self);
    if (module_data == NULL) return NULL;

    default_parameters(&param, module_data);
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!|$dO!", kwlist,
        module_data->element_type, &element,  &PyArray_Type, &rin,
        &energy, module_data->particle_type, &particle)) {
        return NULL;
    }
    lattice_parameters(&param, module_data, NULL, particle);
    energy_parameters(&param, energy);

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

    track_function = get_track_function(module_data, element);
    if (!track_function) return NULL;
    integrator = track_function->FunctionHandle;
    pyintegrator = track_function->PyFunctionHandle;
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
    static char *kwlist[] = {"", "", "energy", "particle", NULL};
    track_function integrator;
    PyObject *pyintegrator;
    double orb[6];
    struct parameters param;
    struct LibraryListElement *track_function;
    npy_intp dims[2] = {6, 6};

    /* Positional parameters */
    PyObject *element;
    PyArrayObject *rin;
    /* Keyword parameters */
    double energy = 0.0;
    PyObject *particle = NULL;

    at_state *module_data = (at_state*)PyModule_GetState(self);
    if (module_data == NULL) return NULL;

    default_parameters(&param, module_data);
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!|$dO!", kwlist,
        module_data->element_type, &element,  &PyArray_Type, &rin,
        &energy, module_data->particle_type, &particle)) {
        return NULL;
    }
    lattice_parameters(&param, module_data, NULL, particle);
    energy_parameters(&param, energy);

    if (PyArray_DIM(rin,0) != 6) {
        return PyErr_Format(PyExc_ValueError, "rin is not 6D");
    }
    if (PyArray_TYPE(rin) != NPY_DOUBLE) {
        return PyErr_Format(PyExc_ValueError, "rin is not a double array");
    }
    if ((PyArray_FLAGS(rin) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
        return PyErr_Format(PyExc_ValueError, "rin is not Fortran-aligned");
    }
    npy_uint32 num_particles = (PyArray_SIZE(rin)/6);
    if (num_particles != 1) {
        return PyErr_Format(PyExc_ValueError, "Number of particles should be 1");
    }
    const double *drin = (const double *)PyArray_DATA(rin);

    for (int i = 0; i < 6; i++) orb[i] = drin[i];

    PyObject *pbdiff = PyArray_ZEROS(2, dims, NPY_DOUBLE ,1);
    double *bdiff = (double *)PyArray_DATA((PyArrayObject *)pbdiff);
    param.bdiff = bdiff;

    track_function = get_track_function(module_data, element);
    if (!track_function) return error_return(0, 1, pbdiff);
    integrator = track_function->FunctionHandle;
    pyintegrator = track_function->PyFunctionHandle;
    if (pyintegrator) {
        PyObject *res = PyObject_CallFunctionObjArgs(pyintegrator, rin, element, NULL);
        if (!res) return error_return(0, 1, pbdiff);
        Py_DECREF(res);
    } else {
        struct elem *elem_data = integrator(element, NULL, orb, 1, &param);
        if (!elem_data) return error_return(0, 1, pbdiff);
        free(elem_data);
    }
    return pbdiff;
}


static int catpass_exec(PyObject *self)
{
    at_state *module_data = (at_state *)PyModule_GetState(self);
    if (module_data == NULL) return -1;

    if (PyArray_ImportNumPyAPI() < 0) return -1;

    module_data -> LibraryList = NULL;

    /* get Particle type */
    PyObject *particle_type = get_pyobj("at.lattice", "Particle");
    if (particle_type == NULL) return -1;
    module_data->particle_type = particle_type;

    /* get Element type */
    PyObject *element_type = get_pyobj("at.lattice", "Element");
    if (element_type == NULL) return -1;
    module_data->element_type = element_type;

    /* get electron mass */
    PyObject *electron_mass = get_pyobj("at.constants", "e_mass");
    if (electron_mass == NULL) return -1;
    module_data->electron_mass = PyFloat_AsDouble(electron_mass);
    Py_DECREF(electron_mass);

    /* get speed of light */
    PyObject *pyclight = get_pyobj("at.constants", "clight");
    if (pyclight == NULL) return -1;
    module_data->clight = PyFloat_AsDouble(pyclight);
    Py_DECREF(pyclight);
    return 0;
}

static void catpass_free(void *self)
{
    at_state *module_data = (at_state *)PyModule_GetState(self);
    Py_XDECREF(module_data->particle_type);
    Py_XDECREF(module_data->element_type);
}

/* Method table */

static PyMethodDef AtMethods[] = {
    {"line_store",  (PyCFunction)line_store, METH_VARARGS | METH_KEYWORDS,
    PyDoc_STR("line_store(lattice, *, energy=None, Particle=None)\n\n"
              "Allocate memory for a lattice.\n\n"
              "Parameters:\n"
              "    lattice (Lattice):          AT lattice\n\n"
              "Keyword Arguments:\n"
              "    energy (float | None):      nominal energy [eV]. Default: from the lattice\n"
              "    particle (Particle | None): circulating particle. Default: from the lattice\n"
            )},
    {"line_pass",  (PyCFunction)line_pass, METH_VARARGS | METH_KEYWORDS,
    PyDoc_STR("(linedata, r_in, n_turns=1, refpts=at.End, *, "
              "turn=0, energy=None, Particle=None, keep_counter=False, losses=False, "
              "start_elem=0, end_elem=nb_elems, omp_num_threads = 0)\n\n"
              "Track input particles *r_in* along *line* for *nturns* turns.\n"
              "Record 6D phase space at elements corresponding to *refpts* for each turn.\n"
              "\n"
              "Parameters:\n"
              "    linedata:           lattice data as returned by *line_store*\n"
              "    rin:                6 x n_particles Fortran-ordered numpy array.\n"
              "      On return, rin contains the final coordinates of the particles\n\n"
              "Keyword Arguments:\n"
              "    n_turns (int):      number of turns to be tracked. Default: 1\n"
              "    refpts (Uint32_refs): numpy array of indices of elements where output is desired\n"
              "       0 means entrance of the first element\n"
              "       len(line) means end of the last element\n"
              "    energy (float | None): nominal energy [eV]. Default: from the lattice\n"
              "    particle (Particle | None): circulating particle. Default: from the lattice\n"
              "    turn (int):         Initial turn number (ignored if *keep_counter* is True)."
              "        Default: 0\n"
              "    keep_counter (int): Start with the current turn number\n"
              "    omp_num_threads (int): number of OpenMP threads (default 0: automatic)\n"
              "    losses (bool):      if True, process losses\n"
              "    start_elem (int):   Start tracking at the specified element. The\n"
              "       behavior is different from :py:meth:`.Lattice.rotate`, *refpts* are\n"
              "       still indexed from the beginning of the ring and *s_coord* starts at\n"
              "       the beginning of the ring. Default 0.\n"
              "    end_elem (int):     End tracking at the specified element (exclusive).\n"
              "       When end_elem is different from the number of lattice elements,\n"
              "       nturns=1 is required. Default number of elements.\n\n"
              "Returns:\n"
              "    rout:    6 x n_particles x n_refpts x n_turns Fortran-ordered numpy array\n"
              "         of particle coordinates\n\n"
              )},
    {"elem_pass",  (PyCFunction)elem_pass, METH_VARARGS | METH_KEYWORDS,
    PyDoc_STR("elem_pass(element, r_in, *, energy=1.0E9, Particle=Particle('relativistic')\n\n"
              "Track input particles *r_in* through a single element.\n\n"
              "Parameters:\n"
              "    element (Element):          AT element\n"
              "    rin:                        6 x n_particles Fortran-ordered numpy array.\n"
              "      On return, rin contains the final coordinates of the particles\n\n"
              "Keyword Arguments:\n"
              "    energy (float | None):      nominal energy [eV]. Default: 1.0E9\n"
              "    particle (Particle | None): circulating particle. Default: relativistic\n\n"
            )},
    {"diffusion_matrix",  (PyCFunction)at_diffmatrix, METH_VARARGS | METH_KEYWORDS,
    PyDoc_STR("diffusion_matrix(element, r_in, *, energy=1.0E9, Particle=Particle('relativistic')\n\n"
              "Track a single particle *r_in* through a single element,\n"
              "computes the diffusion matrix.\n\n"
              "Parameters:\n"
              "    element (Element):          AT element\n"
              "    rin:                        6 x n_particles Fortran-ordered numpy array.\n"
              "      On return, rin contains the final coordinates of the particle\n\n"
              "Keyword Arguments:\n"
              "    energy (float | None):      nominal energy [eV]. Default: 1.0E9\n"
              "    particle (Particle | None): circulating particle. Default: relativistic\n\n"
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

static PyModuleDef_Slot catpass_slots[] = {
    {Py_mod_exec, catpass_exec},
    {Py_mod_exec, catpass_random_exec},
#ifdef Py_GIL_DISABLED
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "catpass",        /* m_name */
    PyDoc_STR("Low-level tracking functions"),      /* m_doc */
    sizeof(at_state), /* m_size */
    AtMethods,        /* m_methods */
    catpass_slots,    /* m_slots */
    NULL,             /* m_traverse */
    NULL,             /* m_clear */
    catpass_free,     /* m_free */
};

PyMODINIT_FUNC PyInit_catpass(void)
{
    PyObject *m = PyModuleDef_Init(&moduledef);

    if (m == NULL) return NULL;

    /* Build path for loading Python integrators */
    PyObject *integ_path_obj = get_integrators();
    if (integ_path_obj) {
        const char *integ_path = PyUnicode_AsUTF8(integ_path_obj);
        PyObject *ext_suffix_obj = get_ext_suffix();
        if (ext_suffix_obj) {
            const char *ext_suffix = (ext_suffix_obj == Py_None) ? OBJECTEXT : PyUnicode_AsUTF8(ext_suffix_obj);
            snprintf(integrator_path, sizeof(integrator_path), "%s%s%%s%s", integ_path, SEPARATOR, ext_suffix);
            Py_DECREF(ext_suffix_obj);
        }
        else {
            return NULL;
        }
        Py_DECREF(integ_path_obj);
    }
    else {
        return NULL;
    }
    return m;
}
