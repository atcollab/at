/*
 * PCG Random Number Generation for C.
 *
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *     http://www.pcg-random.org
 */
#include <math.h>
#include <stdint.h>

/* initial RNG definitions */
#define AT_RNG_STATE 0x853c49e6748fea9bULL
#define AT_RNG_INC 0xda3e39cb94b95bdbULL

#define COMMON_PCG32_INITIALIZER   { AT_RNG_STATE, AT_RNG_INC, 0.0, false }
#define THREAD_PCG32_INITIALIZER   { AT_RNG_STATE, 1ULL, 0.0, false }

struct pcg_state_setseq_64 {    // Internals are *Private*.
    uint64_t state;             // RNG state.  All values are possible.
    uint64_t inc;               // Controls which RNG sequence (stream) is
                                // selected. Must *always* be odd.
    double spare;               // spare value for normal distribution
    bool hasSpare;
};
typedef struct pcg_state_setseq_64 pcg32_random_t;

// If you *must* statically initialize it, here's one.

#define PCG32_INITIALIZER   { AT_RNG_STATE, AT_RNG_INC, 0.0, false }

static pcg32_random_t pcg32_global = PCG32_INITIALIZER;

// pcg32_random()
// pcg32_random_r(rng)
//     Generate a uniformly distributed 32-bit random number

static uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    if (!rng) rng = &pcg32_global;
    uint64_t oldstate;
    #pragma omp atomic read
    oldstate = rng->state;
    #pragma omp atomic write
    rng->state = oldstate * 6364136223846793005ULL + rng->inc;
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

static uint32_t pcg32_random()
{
    return pcg32_random_r(&pcg32_global);
}

// pcg32_srandom(initstate, initseq)
// pcg32_srandom_r(rng, initstate, initseq):
//     Seed the rng.  Specified in two parts, state initializer and a
//     sequence selection constant (a.k.a. stream id)

static void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
{
    if (!rng) rng = &pcg32_global;
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
}

static void pcg32_srandom(uint64_t seed, uint64_t seq)
{
    pcg32_srandom_r(&pcg32_global, seed, seq);
}

/* Functions for uniform, normal and Poisson distributions with
   external state variable */

static double atrandd_r(pcg32_random_t* rng)
/* Uniform [0, 1) distribution */
{
    return ldexp(pcg32_random_r(rng), -32);
}

static double atrandn_r(pcg32_random_t* rng, double mean, double stdDev)
/* gaussian distribution */
{
    /* Marsaglia polar method: https://en.wikipedia.org/wiki/Marsaglia_polar_method */

	double u, v, s;

	if (rng->hasSpare) {
		rng->hasSpare = false;
		return mean + stdDev * rng->spare;
	}

	rng->hasSpare = true;
	do {
		u = 2.0 * atrandd_r(rng) - 1.0;
		v = 2.0 * atrandd_r(rng) - 1.0;
		s = u * u + v * v;
	}
	while ((s >= 1.0) || (s == 0.0));
	s = sqrt(-2.0 * log(s) / s);
	rng->spare = v * s;
	return mean + stdDev * u * s;
}

static int atrandp_r(pcg32_random_t* rng, double lamb)
/* poisson distribution */
{
    int pk;

    if (lamb<11) {
        double l = -lamb;
        int k = 0;
        double p = 0;
        do {
            k += 1;
            p += log(atrandd_r(rng));
        } while (p>l);
        pk = k-1;
    }
    else {      /* Gaussian approximation */
        pk = (int)floor(atrandn_r(rng, lamb, sqrt(lamb)));
    }

    return pk;
}

/* Functions for uniform, normal and Poisson distributions with
   internal state variable */

static inline double atrandd(void)
{
    return atrandd_r(&pcg32_global);
}

static inline double atrandn(double mean, double stdDev)
{
    return atrandn_r(&pcg32_global, mean, stdDev);
}

static inline int atrandp(double lamb)
{
    return atrandp_r(&pcg32_global, lamb);
}
