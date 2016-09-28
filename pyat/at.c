/*
 * 1. Register single function atpass
 * 2. Interpret list of objects
 * 3. Go round ring
 * 4. Call functions
 * 5. Return numpy array.
 */

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>
#include <stdio.h>
#include <stdlib.h>
#include "at.h"

// Linux only
#include <dlfcn.h>

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

#define MAX_ORDER 3
#define MAX_INT_STEPS 5
#ifndef INTEGRATOR_PATH
#define INTEGRATOR_PATH "../atintegrators"
#endif /*INTEGRATOR_PATH*/
#define ATPY_PASS "atpyPass"

typedef int (*pass_function)(double *rin, int num_particles, PyObject *element, struct parameters *param);

/* Directly copied from atpass.c */
static struct LibraryListElement {
    char *LibraryFileName;
    char *MethodName;
    void *FunctionHandle;
    struct LibraryListElement *Next;
} *LibraryList = NULL;

struct LibraryListElement* SearchLibraryList(struct LibraryListElement *head, const char *method_name)
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
        char lib_file[300];
        snprintf(lib_file, sizeof(lib_file), "%s/%s.so", INTEGRATOR_PATH, fn_name);
        void *dl_handle = dlopen(lib_file, RTLD_LAZY);
        if (dl_handle == NULL) {
            PyErr_SetString(PyExc_RuntimeError, dlerror());
            return NULL;
        }
        fn_handle = dlsym(dl_handle, ATPY_PASS);
        if (fn_handle == NULL) {
            PyErr_SetString(PyExc_RuntimeError, dlerror());
            return NULL;
        }
        LibraryListPtr = (struct LibraryListElement *)malloc(sizeof(struct LibraryListElement));
        LibraryListPtr->Next = LibraryList;
        LibraryListPtr->MethodName = fn_name;
        LibraryListPtr->FunctionHandle = fn_handle;
        LibraryList = LibraryListPtr;
    }
    return fn_handle;
}


static int pass_element(double *rin, int num_particles, PyObject *element, struct parameters *param) {
    pass_function fn_handle = NULL;
    PyObject *fn_name_object = PyObject_GetAttrString(element, "PassMethod");
    if (fn_name_object && (fn_handle = pass_method(PyUnicode_AsUTF8(fn_name_object)))) {
        return fn_handle(rin, num_particles, element, param);
    }
    else {
        return -1;
    }
}


/*
 * Arguments:
 *  - the_ring: sequence of elements
 *  - rin: numpy 6-vector of initial conditions
 *  - num_turns: int number of turns to simulate
 */
static PyObject *at_atpass(PyObject *self, PyObject *args) {
    PyObject *element_list;
    PyArrayObject *rin;
    double *drin;
    int num_turns;
    int num_parts;
    int i, j;
    struct parameters param;
    param.nturn = 0;
    param.mode = 0;
    param.T0 = 0;
    param.RingLength = 0;

    if (!PyArg_ParseTuple(args, "O!O!i", &PyList_Type, &element_list, &PyArray_Type, &rin, &num_turns)) {
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
    num_parts = (int)(PyArray_SIZE(rin)/6);
    drin = PyArray_DATA(rin);

    long num_elements = PyList_Size(element_list);
    printf("There are %ld elements in the list\n", num_elements);
    printf("There are %d particles\n", num_parts);
    printf("Going for %d turns\n", num_turns);
    for (i = 0; i < num_turns; i++) {
        param.nturn = i;
        for (j = 0; j < num_elements; j++) {
            PyObject *element = PyList_GET_ITEM(element_list, j);
            if (pass_element(drin, num_parts, element, &param) != 0) {
                char *pass_error_template = "Error occurred during pass method for element %d";
                if (!PyErr_Occurred()) {
                    char pass_error[50];
                    snprintf(pass_error, sizeof(pass_error), pass_error_template, j);
                    PyErr_SetString(PyExc_RuntimeError, pass_error);
                } else {
                    printf(pass_error_template, j);
                    printf(".\n");
                }
                return NULL;
            }
        }
    }
    return Py_BuildValue("i", 1);
}

/* Boilerplate to register methods. */

static PyMethodDef AtMethods[] = {
    {"atpass",  at_atpass, METH_VARARGS,
    PyDoc_STR("atpass(line,rin,nturns)\n\nTrack rin along line for nturns turns")},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

MOD_INIT(at)
{
#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"at",			/* m_name */
	PyDoc_STR("Clone of atpass in Accelerator Toolbox"),      /* m_doc */
	-1,             /* m_size */
	AtMethods,		/* m_methods */
	NULL,			/* m_reload */
	NULL,			/* m_traverse */
	NULL,			/* m_clear */
	NULL,			/* m_free */
    };
    PyObject *m = PyModule_Create(&moduledef);
#else
    PyObject *m = Py_InitModule3("at", AtMethods,
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
    initat();
    return 0;
}
#endif
