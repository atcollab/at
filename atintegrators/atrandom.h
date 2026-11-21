#ifndef ATRANDOM_H
#define ATRANDOM_H
#include <stdbool.h>

struct pcg_state_setseq_64 {    // Internals are *Private*.
    uint64_t state;             // RNG state.  All values are possible.
    uint64_t inc;               // Controls which RNG sequence (stream) is
                                // selected. Must *always* be odd.
    double spare;               // spare value for normal distribution
    bool hasSpare;
};
typedef struct pcg_state_setseq_64 pcg32_random_t;

double atrandd_r(pcg32_random_t* rng);
double atrandn_r(pcg32_random_t* rng, double mean, double stdDev);
int atrandp_r(pcg32_random_t* rng, double lamb);

#ifdef MATLAB_MEX_FILE
#include "atrandom.c"
#endif /* MATLAB_MEX_FILE */

#endif /* ATRANDOM_H */
