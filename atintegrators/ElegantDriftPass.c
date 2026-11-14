/*
 * Exact Hamiltonian drift corresponding to elegant's EDRIFT element.
 *
 * elegant stores transverse slopes, while PyAT stores canonical momenta
 * normalized to the reference momentum. For
 *
 *   pz = sqrt((1 + delta)^2 - px^2 - py^2),
 *
 * the equivalent elegant slopes are xp = px/pz and yp = py/pz. Applying
 * elegant's exact drift in PyAT coordinates therefore gives the map below
 * directly, without explicit coordinate transformations.
 *
 * Portions of the algorithm and formulas are based on Elegant source code:
 *   https://github.com/rtsoliday/elegant
 *
 * Elegant is copyright (c) 2002 University of Chicago. See
 * doc/licenseFiles/ElegantLicense.txt for the Elegant license terms.
 *
 * This file adapts the relevant Elegant method to AT/PyAT data structures,
 * canonical coordinates, aperture handling, and pass-method conventions.
 * Not all Elegant element features are implemented.
 */

#include "atelem.c"
#include "atlalib.c"

struct elem {
    double Length;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *EApertures;
    double *RApertures;
};

static int exact_drift(double *r6, double length)
{
    double dp1 = 1.0 + r6[delta_];
    double pz2 = dp1 * dp1
               - r6[px_] * r6[px_]
               - r6[py_] * r6[py_];
    double normalized_length;

    if (pz2 <= 0.0)
        return 0;

    normalized_length = length / sqrt(pz2);
    r6[x_] += r6[px_] * normalized_length;
    r6[y_] += r6[py_] * normalized_length;
    r6[ct_] += dp1 * normalized_length;
    r6[ct_] -= length;
    return 1;
}

static void elegant_drift(double *r_in, const struct elem *Elem,
                          int num_particles)
{
    int c;

    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD * 10) \
        default(none) shared(r_in, Elem, num_particles)
    for (c = 0; c < num_particles; c++) {
        double *r6 = r_in + 6 * c;

        if (!atIsNaN(r6[x_])) {
            if (Elem->T1)
                ATaddvv(r6, Elem->T1);
            if (Elem->R1)
                ATmultmv(r6, Elem->R1);
            if (Elem->RApertures)
                checkiflostRectangularAp(r6, Elem->RApertures);
            if (Elem->EApertures)
                checkiflostEllipticalAp(r6, Elem->EApertures);

            if (!atIsNaN(r6[x_]) && !exact_drift(r6, Elem->Length))
                r6[x_] = atGetNaN();

            if (Elem->RApertures)
                checkiflostRectangularAp(r6, Elem->RApertures);
            if (Elem->EApertures)
                checkiflostEllipticalAp(r6, Elem->EApertures);
            if (Elem->R2)
                ATmultmv(r6, Elem->R2);
            if (Elem->T2)
                ATaddvv(r6, Elem->T2);
        }
    }
}

static struct elem *initialize_elem(const atElem *ElemData, struct elem *Elem)
{
    Elem->Length = atGetDouble(ElemData, "Length"); check_error();
    Elem->R1 = atGetOptionalDoubleArray(ElemData, "R1"); check_error();
    Elem->R2 = atGetOptionalDoubleArray(ElemData, "R2"); check_error();
    Elem->T1 = atGetOptionalDoubleArray(ElemData, "T1"); check_error();
    Elem->T2 = atGetOptionalDoubleArray(ElemData, "T2"); check_error();
    Elem->EApertures =
        atGetOptionalDoubleArray(ElemData, "EApertures"); check_error();
    Elem->RApertures =
        atGetOptionalDoubleArray(ElemData, "RApertures"); check_error();
    return Elem;
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData, struct elem *Elem,
                                      double *r_in, int num_particles,
                                      struct parameters *Param)
{
    (void)Param;

    if (!Elem) {
        Elem = (struct elem *)atMalloc(sizeof(struct elem));
        Elem = initialize_elem(ElemData, Elem); check_error();
    }

    elegant_drift(r_in, Elem, num_particles);
    return Elem;
}

MODULE_DEF(ElegantDriftPass)
#endif

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        struct elem Elem;
        double *r_in;
        int num_particles = mxGetN(prhs[1]);

        if (mxGetM(prhs[1]) != 6)
            mexErrMsgIdAndTxt("AT:WrongArg",
                              "Second argument must be a 6 x N matrix");
        initialize_elem(prhs[0], &Elem);
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        elegant_drift(r_in, &Elem, num_particles);
    } else if (nrhs == 0) {
        plhs[0] = mxCreateCellMatrix(1, 1);
        mxSetCell(plhs[0], 0, mxCreateString("Length"));
        if (nlhs > 1) {
            plhs[1] = mxCreateCellMatrix(6, 1);
            mxSetCell(plhs[1], 0, mxCreateString("R1"));
            mxSetCell(plhs[1], 1, mxCreateString("R2"));
            mxSetCell(plhs[1], 2, mxCreateString("T1"));
            mxSetCell(plhs[1], 3, mxCreateString("T2"));
            mxSetCell(plhs[1], 4, mxCreateString("EApertures"));
            mxSetCell(plhs[1], 5, mxCreateString("RApertures"));
        }
    } else {
        mexErrMsgIdAndTxt("AT:WrongArg", "Needs 0 or 2 arguments");
    }
}
#endif
