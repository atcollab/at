#ifndef __TRACK_H__
#define __TRACK_H__

/* element interface */

enum { FMAX = 32 };

enum element_type
{
    drift = 0,
    dipole,
    multipole,
    marker
};

struct element
{
    double L;
    double phi;
    double gK;
    double F[FMAX];
    int nF;
    int slices;
    int type;
    int do_multipole_fringe;
};

struct lattice
{
    element * next;
    int N;
};

template<typename T> void track_element(T * x, element * e);

#endif

