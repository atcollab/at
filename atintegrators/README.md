How to write an integrator
==========================

The core functionality of AT is built around 'integrators', which act on
a 6D vector representing a particle in phase space. A number of different
integrators are available, representing different types of elements or
different methods of calculation. Each integrator is known as a 'pass method'.

The components of an integrator are as follows, illustrated by the pass
method DriftPass:


### Element parameters

The names and types required by the element using a pass method are defined
in a struct:

    struct elem 
    {
      double Length;
      /* Optional fields */
      double *R1;
      double *R2;
      double *T1;
      double *T2;
      double *EApertures;
      double *RApertures;
    };


### Physics calculation

Each pass method should have a function that takes the array r_in representing 
the input particles plus the parameters defined in the struct and alters the
particles in place:


    void DriftPass(double *r_in, double le,
               const double *T1, const double *T2,
               const double *R1, const double *R2,
               double *RApertures, double *EApertures,
               int num_particles)
    {
        ...
    }

This function should not use any Matlab or Python functionality.

### Track function

The function trackFunction will be used either by Python or Matlab, and uses
generic functions defined in atelem.c. Those functions have implementations
for both Python and Matlab - the correct ones are chosen based on the
variables MATLAB_MEX_FILE and PYAT:


    #if defined(MATLAB_MEX_FILE) || defined(PYAT)
    ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
                double *r_in, int num_particles, struct parameters *Param)
    {
        ...
    }

For PYAT builds, you need to include

    MODULE_DEF(DriftPass)

This does setup for Python imports.


### Mex function

This is defined only for Matlab builds and allows you to call the integrator
directly from Matlab:

    
    #if defined(MATLAB_MEX_FILE)
    void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        ...
    }

