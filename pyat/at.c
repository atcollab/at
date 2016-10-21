/*
 * 1. Register single function atpass
 * 2. Interpret list of objects
 * 3. Go round ring
 * 4. Call functions
 * 5. Return numpy array.
 */

#include <Python.h>
#include "attypes.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

#if PY_MAJOR_VERSION >= 3
#define NUMPY_IMPORT_ARRAY_RETVAL NULL
#define NUMPY_IMPORT_ARRAY_TYPE void *
#else
#define NUMPY_IMPORT_ARRAY_RETVAL
#define NUMPY_IMPORT_ARRAY_TYPE void
#define PyLong_AsLong PyInt_AsLong
#endif

typedef PyObject atElem;

#define ATPY_PASS "trackFunction"

#if defined(PCWIN) || defined(PCWIN64) || defined(_WIN32)
#include <windows.h>
#define LIBRARYHANDLETYPE HINSTANCE
#define FREELIBFCN(libfilename) FreeLibrary((libfilename))
#define LOADLIBFCN(libfilename) LoadLibrary((libfilename))
#define GETTRACKFCN(libfilename) GetProcAddress((libfilename),ATPY_PASS)
#else
#include <dlfcn.h>
#define LIBRARYHANDLETYPE void *
#define FREELIBFCN(libfilename) dlclose(libfilename)
#define LOADLIBFCN(libfilename) dlopen((libfilename),RTLD_LAZY)
#define GETTRACKFCN(libfilename) dlsym((libfilename),ATPY_PASS)
#endif

#if PY_MAJOR_VERSION >= 3
  #define MOD_ERROR_VAL NULL
  #define MOD_SUCCESS_VAL(val) val
  #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#else
  #define MOD_ERROR_VAL
  #define MOD_SUCCESS_VAL(val)
  #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
  #define PyUnicode_AsUTF8 PyString_AsString
#endif

#ifndef INTEGRATOR_PATH
#define INTEGRATOR_PATH "../atintegrators/%s.so"
#endif /*INTEGRATOR_PATH*/

typedef struct elem *(*pass_function)(const PyObject *element, struct elem *elemptr,
        double *r_in, int num_particles, struct parameters *param);

static npy_uint32 num_elements = 0;
static struct elem **elemdata_list = NULL;
static PyObject **element_list = NULL;
static pass_function *integrator_list = NULL;

/* Directly copied from atpass.c */
static struct LibraryListElement {
    const char *MethodName;
    LIBRARYHANDLETYPE LibraryHandle;
    pass_function FunctionHandle;
    struct LibraryListElement *Next;
} *LibraryList = NULL;

static PyObject *print_error(int elem_number, PyObject *rout)
{
    printf("Error in tracking element %d\n", elem_number);
    Py_XDECREF(rout);
    return NULL;
}

static struct LibraryListElement* SearchLibraryList(struct LibraryListElement *head, const char *method_name)
{
    /* recusively search the list to check if the library containing method_name is
     * already loaded. If it is - return the pointer to the list element. If not -
     * return NULL */
    if (head)
        return (strcmp(head->MethodName, method_name)==0) ? head :
            SearchLibraryList(head->Next, method_name);
    else
        return NULL;
}

static pass_function pass_method(char *fn_name) {
    pass_function fn_handle = NULL;
    struct LibraryListElement *LibraryListPtr = SearchLibraryList(LibraryList, fn_name);

    if (LibraryListPtr) {
        fn_handle = LibraryListPtr->FunctionHandle;
    }
    else {
        char lib_file[300], buffer[200];
        LIBRARYHANDLETYPE dl_handle;
        snprintf(lib_file, sizeof(lib_file), INTEGRATOR_PATH, fn_name);
        dl_handle = LOADLIBFCN(lib_file);
        if (dl_handle == NULL) {
            snprintf(buffer, sizeof(buffer), "Cannot load %s", lib_file);
            PyErr_SetString(PyExc_RuntimeError, buffer);
            return NULL;
        }
        fn_handle = GETTRACKFCN(dl_handle);
        if (fn_handle == NULL) {
            snprintf(buffer, sizeof(buffer), "No trackFunction in %s", lib_file);
            FREELIBFCN(dl_handle);
            PyErr_SetString(PyExc_RuntimeError, buffer);
            return NULL;
        }
        LibraryListPtr = (struct LibraryListElement *)malloc(sizeof(struct LibraryListElement));
        LibraryListPtr->MethodName = strcpy(malloc(strlen(fn_name)+1), fn_name);
        LibraryListPtr->LibraryHandle = dl_handle;
        LibraryListPtr->FunctionHandle = fn_handle;
        LibraryListPtr->Next = LibraryList;
        LibraryList = LibraryListPtr;
    }
    return fn_handle;
}

/*
 * Arguments:
 *  - the_ring: sequence of elements
 *  - rin: numpy 6-vector of initial conditions
 *  - num_turns: int number of turns to simulate
 */
static PyObject *at_atpass(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"line","rin","nturns","refpts","reuse", NULL};
    static int new_lattice = 1;
    static double lattice_length = 0.0;

    PyObject *lattice;
    PyArrayObject *rin;
    PyArrayObject *refs = NULL;
    PyObject *rout;
    double *drin, *drout;
    int num_turns;
    npy_uint32 num_particles, np6;
    npy_uint32 numel;
    npy_uint32 *refpts = NULL;
    npy_uint32 nextref;
    unsigned int nextrefindex;
    unsigned int num_refpts;
    unsigned int reuse=0;
    npy_intp outdims[2];
    int turn, nelem;
    struct parameters param;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!i|O!I", kwlist, &PyList_Type, &lattice,
        &PyArray_Type, &rin, &num_turns, &PyArray_Type, &refs, &reuse)) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse arguments to atpass");
        return NULL;
    }
    if (PyArray_DIM(rin,PyArray_NDIM(rin)-1) != 6) {
        PyErr_SetString(PyExc_ValueError, "Numpy array is not 6D");
        return NULL;
    }
    if (PyArray_TYPE(rin) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "rin is not a double array");
        return NULL;
    }
    if ((PyArray_FLAGS(rin) & NPY_ARRAY_CARRAY_RO) != NPY_ARRAY_CARRAY_RO) {
        PyErr_SetString(PyExc_ValueError, "rin is not C aligned");
        return NULL;
    }

    numel = PyList_Size(lattice);
    num_particles = (PyArray_SIZE(rin)/6);
    np6 = num_particles*6;
    drin = PyArray_DATA(rin);

    if (refs) {
        if (PyArray_TYPE(refs) != NPY_UINT32) {
            PyErr_SetString(PyExc_ValueError, "refpts is not a double array");
            return NULL;
        }
        if ((PyArray_FLAGS(refs) & NPY_ARRAY_CARRAY_RO) != NPY_ARRAY_CARRAY_RO) {
            PyErr_SetString(PyExc_ValueError, "refpts is not C aligned");
            return NULL;
        }
        refpts = PyArray_DATA(refs);
        num_refpts = PyArray_SIZE(refs);
        if (num_refpts == 0)
            outdims[0] = num_particles;
        else
            outdims[0] = num_turns*num_refpts*num_particles;
    }
    else {              /* only end of the line */
        refpts = &num_elements;
        num_refpts = 1;
        outdims[0] = num_turns*num_particles;
    }
    outdims[1] = 6;
    rout = PyArray_SimpleNew(2, outdims, NPY_DOUBLE);
    drout = PyArray_DATA((PyArrayObject *)rout);

    if (!reuse) new_lattice = 1;
    if (new_lattice) {
		PyObject **element;
		pass_function *integrator;
        int nelem;
        for (nelem=0; nelem<num_elements; nelem++) {
            free(elemdata_list[nelem]);
            Py_XDECREF(element_list[nelem]);        /* Release the stored elements, may be NULL if */
        }                                           /* a previous call was interrupted by an error */
        num_elements = numel;

        /* Pointer to Element structures used by the tracking function */
        free(elemdata_list);
        elemdata_list = (struct elem **)calloc(num_elements, sizeof(struct elem *));

        /* Pointer to Element list, make sure all pointers are initially NULL */
        free(element_list);
        element_list = (PyObject **)calloc(num_elements, sizeof(PyObject *));

        /* pointer to the list of integrators */
        integrator_list = (pass_function *)realloc(integrator_list, num_elements*sizeof(pass_function));

        lattice_length = 0.0;
        element = element_list;
        integrator = integrator_list;
        for (nelem = 0; nelem < num_elements; nelem++) {
            PyObject *el = PyList_GET_ITEM(lattice, nelem);
            PyObject *PyPassMethod = PyObject_GetAttrString(el, "PassMethod");
			pass_function fn_handle;
			double length;
            if (!PyPassMethod) return print_error(nelem, rout);     /* No PassMethod */
            fn_handle = pass_method(PyUnicode_AsUTF8(PyPassMethod));
            if (!fn_handle) return print_error(nelem, rout);        /* No trackFunction for the given PassMethod */
            length = PyFloat_AsDouble(PyObject_GetAttrString(el, "Length"));
            if (PyErr_Occurred()) PyErr_Clear();
            else lattice_length += length;
            *integrator++ = fn_handle;
            *element++ = el;
            Py_INCREF(el);                          /* Keep a reference to each element in case of reuse */
        }
        new_lattice = 0;
    }

    printf("There are %u elements in the list\n", num_elements);
    printf("There are %u particles\n", num_particles);
    printf("Going for %u turns\n", num_turns);

    param.RingLength = lattice_length;
    param.T0 = lattice_length/299792458.0;
    for (turn = 0; turn < num_turns; turn++) {
        PyObject **element = element_list;
        pass_function *integrator = integrator_list;
        struct elem **elemdata = elemdata_list;
        param.nturn = turn;
        nextrefindex = 0;
        nextref= (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
        for (nelem = 0; nelem < num_elements; nelem++) {
            if (nelem == nextref) {
                memcpy(drout, drin, np6*sizeof(double));
                drout += np6; /*  shift the location to write to in the output array */
                nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
            }
            *elemdata = (*integrator++)(*element++, *elemdata, drin, num_particles, &param);
            if (!*elemdata) return print_error(nelem, rout);       /* trackFunction failed */
            elemdata++;
        }
        if (num_elements == nextref) {
            memcpy(drout, drin, np6*sizeof(double));
            drout += np6; /*  shift the location to write to in the output array */
        }
    }
    if (num_refpts == 0) {
        memcpy(drout, drin, np6*sizeof(double));
        drout += np6; /*  shift the location to write to in the output array */
    }
    return rout;
}

/* Boilerplate to register methods. */

static PyMethodDef AtMethods[] = {
    {"atpass",  (PyCFunction)at_atpass, METH_VARARGS | METH_KEYWORDS,
    PyDoc_STR("atpass(line,rin,nturns)\n\nTrack rin along line for nturns turns")},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

MOD_INIT(atpass)
{
#if PY_MAJOR_VERSION >= 3
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
#else
    PyObject *m = Py_InitModule3("atpass", AtMethods,
        "Clone of atpass in Accelerator Toolbox");
#endif
    if (m == NULL)
       return MOD_ERROR_VAL;
    import_array();
    return MOD_SUCCESS_VAL(m);
}

#if PY_MAJOR_VERSION < 3
int main(int argc, char *argv[]) {
    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(argv[0]);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Add a static module */
    initatpass();
    return 0;
}
#endif
