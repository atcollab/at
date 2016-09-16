#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

static int array_imported = 0;

#if PY_MAJOR_VERSION >= 3
#define NUMPY_IMPORT_ARRAY_RETVAL NULL
#define NUMPY_IMPORT_ARRAY_TYPE void *
#else
#define NUMPY_IMPORT_ARRAY_RETVAL
#define NUMPY_IMPORT_ARRAY_TYPE void
#define PyLong_AsLong PyInt_AsLong
#endif

static NUMPY_IMPORT_ARRAY_TYPE init_numpy(void) {
    import_array();
    return NUMPY_IMPORT_ARRAY_RETVAL;
}

static long py_get_long(PyObject *element, char *name) {
    return PyLong_AsLong(PyObject_GetAttrString(element, name));
}

static double py_get_double(PyObject *element, char *name) {
    return PyFloat_AsDouble(PyObject_GetAttrString(element, name));
}

static double *numpy_get_double_array(PyObject *element, char *name) {
    if (!array_imported) {
        init_numpy();
        array_imported = 1;
    }
    PyObject *previous_error = PyErr_Occurred();    /* Error occurred in previous arguments ? */
    PyArrayObject *array = (PyArrayObject *) PyObject_GetAttrString(element, name);
    if (array == NULL) {
        if (previous_error == NULL) PyErr_Clear();  /* keep track of previous argument errors */
        return NULL;
    }
    if (!PyArray_Check(array)) {
        printf("%s not an array\n", name);
        return NULL;
    }
    if (PyArray_TYPE(array) != NPY_DOUBLE) {
        printf("%s not double\n", name);
        return NULL;
    }
    if ((PyArray_FLAGS(array) & NPY_ARRAY_CARRAY_RO) != NPY_ARRAY_CARRAY_RO) {
        printf("%s not C aligned\n", name);
        return NULL;
    }
    return PyArray_DATA(array);
}
