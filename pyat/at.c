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

#ifndef INTEGRATOR_PATH
#define INTEGRATOR_PATH "../atintegrators"
#endif /*INTEGRATOR_PATH*/
#define ATPY_PASS "trackFunction"

static int nb_allocated_elements = 0;
static struct elem **ElemStruct_ptr = NULL;

typedef struct elem *(*pass_function)(const PyObject *element, struct elem *elemptr,
        double *r_in, int num_particles, struct parameters *param);

/* Directly copied from atpass.c */
static struct LibraryListElement {
    char *LibraryFileName;
    char *MethodName;
    void *FunctionHandle;
    struct LibraryListElement *Next;
} *LibraryList = NULL;

static PyObject *print_error(int elem_number)
{
    printf("Error in tracking element %d\n", elem_number);
    return NULL;
}

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

    if (1) {
        int n;
        for (n=0; n<nb_allocated_elements; n++) {
            free(ElemStruct_ptr[n]);
        }
        free(ElemStruct_ptr);
        ElemStruct_ptr = (struct elem **)calloc(num_elements, sizeof(struct elem *));
        nb_allocated_elements = num_elements;
    }

    for (i = 0; i < num_turns; i++) {
        struct elem **elemptr = ElemStruct_ptr;
        param.nturn = i;
        for (j = 0; j < num_elements; j++) {
            pass_function fn_handle;
            struct elem *result;
            PyObject *element = PyList_GET_ITEM(element_list, j);
            PyObject *fn_name_object = PyObject_GetAttrString(element, "PassMethod");
            if (!fn_name_object) return print_error(j);     /* No PassMethod */
            fn_handle = pass_method(PyUnicode_AsUTF8(fn_name_object));
            if (!fn_handle) return print_error(j);          /* No trackFunction for the given PassMethod */
            result = fn_handle(element, *elemptr, drin, num_parts, &param);
            if (!result) return print_error(j);             /* trackFunction failed */
            *elemptr = result;
            elemptr++;
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
