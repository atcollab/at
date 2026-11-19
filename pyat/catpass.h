#ifndef CATPASSS_H
#define CATPASSS_H
#include "atrandom.h"

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

PyObject *reset_rng(PyObject *self, PyObject *args, PyObject *kwargs);
PyObject *common_rng(PyObject *self);
PyObject *thread_rng(PyObject *self);
int catpass_random_exec(PyObject *self);

#endif /* CATPASS _H */