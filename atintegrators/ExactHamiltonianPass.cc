/* See ExactHamiltonianPass.m for further notes. */
#include "atelem.c"
#include <string.h>
#include "atlalib.c"

#define AT_MODE
#include "track.cc"

struct elem
{
    double Length;
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    int NumIntSteps;
    /* type is defined in track.h:
        - 0: drift
        - 1: dipole
        - 2: multipole
        - 3: marker
     */
    int Type;
    /* Optional fields */
    double gK;  /* g * K, required for bend */
    double BendingAngle; /* required for bend */
    int MultipoleFringe; /* bool, whether to calculate multipole fringe */
    double *R1;
    double *R2;
    double *T1;
    double *T2;
};

static const char * FirstOptionalName = "gK";

#define NUM_FIELDS_2_REMEMBER (sizeof(FieldNames) / sizeof(FieldNames[0]))

void ExactHamiltonianPass(double *r_in, double le,
                     double *A, double *B,
                     const double *T1, const double *T2,
                     const double *R1, const double *R2,
                     int max_order, int num_int_steps,
                     double phi, int type,
                     double gK, int multipole_fringe,
                     int num_particles)
{
    int c, n;
    double * r6;

    // copy AT structure into pass method structure

    element e = {0};

    for(n = 0; n < max_order+1; n++)
    {
        e.F[2 * n]     = B[n];
        e.F[2 * n + 1] = A[n];
    }

    e.L = le;
    e.phi = phi;
    e.gK = gK;
    e.nF = max_order+1;
    e.slices = num_int_steps;
    e.type = type;
    e.do_multipole_fringe = multipole_fringe;

    for(c = 0; c<num_particles; c++)
    {
        r6 = r_in+c*6;
        if(!atIsNaN(r6[0]))
        {
            /* misalignment at entrance */
            if(T1)
                ATaddvv(r6,T1);
            if(R1)
                ATmultmv(r6,R1);

            track_element(r6, &e);

            /* misalignment at exit */
            if(R2)
                ATmultmv(r6,R2);
            if(T2)
                ATaddvv(r6,T2);
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)

ExportMode struct elem *trackFunction(const atElem *ElemData, struct elem *Elem,
                               double *r_in, int num_particles,
                               struct parameters *Param)

{
    if (!Elem) {
        double le, bending_angle;
        double *polynom_a, *polynom_b;
        long max_order, num_int_steps, type, multipole_fringe;
        double *R1, *R2, *T1, *T2;
        double phi, gK;
        le = atGetDouble(ElemData,"Length"); check_error();
        polynom_a = atGetDoubleArray(ElemData,"PolynomA"); check_error();
        polynom_b = atGetDoubleArray(ElemData,"PolynomB"); check_error();
        max_order = atGetLong(ElemData, "MaxOrder"); check_error();
        num_int_steps = atGetLong(ElemData, "NumIntSteps"); check_error();
        type = atGetLong(ElemData, "Type"); check_error();
        multipole_fringe = atGetOptionalLong(ElemData, "MultipoleFringe", 0); check_error();
        bending_angle = atGetOptionalDouble(ElemData,"BendingAngle", 0.0); check_error();
        gK = atGetOptionalDouble(ElemData,"gK", 0.0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=le;
        Elem->PolynomA=polynom_a;
        Elem->PolynomB=polynom_b;
        Elem->MaxOrder=max_order;
        Elem->NumIntSteps=num_int_steps;
        Elem->Type = type;
        /*optional fields*/
        Elem->MultipoleFringe = multipole_fringe;
        Elem->BendingAngle = bending_angle;
        Elem->gK = gK;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
    }
    ExactHamiltonianPass(r_in, Elem->Length, Elem->PolynomA, Elem->PolynomB,
                               Elem->T1, Elem->T2, Elem->R1, Elem->R2,
                               Elem->MaxOrder, Elem->NumIntSteps,
                               Elem->BendingAngle, Elem->Type, Elem->gK,
                               Elem->MultipoleFringe, num_particles);

    return Elem;
}

MODULE_DEF(ExactHamiltonianPass) /* Dummy module initialisation */

#endif /* defined(PYAT) || defined(MATLAB_MEX_FILE) */


#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double * r_in;
        const mxArray * ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        int m = mxGetM(prhs[1]);
        if (m!=6)
            mexErrMsgTxt("Second argument must be a 6 x N matrix");

        double le, bending_angle;
        double *polynom_a, *polynom_b;
        long max_order, num_int_steps, type, multipole_fringe;
        double *R1, *R2, *T1, *T2;
        double phi, gK;
        le = atGetDouble(ElemData,"Length"); check_error();
        polynom_a = atGetDoubleArray(ElemData,"PolynomA"); check_error();
        polynom_b = atGetDoubleArray(ElemData,"PolynomB"); check_error();
        max_order = atGetLong(ElemData, "MaxOrder"); check_error();
        num_int_steps = atGetLong(ElemData, "NumIntSteps"); check_error();
        type = atGetLong(ElemData, "Type"); check_error();
        /*optional fields*/
        multipole_fringe = atGetOptionalLong(ElemData, "MultipoleFringe", 0); check_error();
        bending_angle = atGetOptionalDouble(ElemData,"BendingAngle", 0.0); check_error();
        gK = atGetOptionalDouble(ElemData,"gK", 0.0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();

        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        ExactHamiltonianPass(r_in, le, polynom_a, polynom_b, T1, T2, R1, R2, max_order,
                     num_int_steps, bending_angle, type, gK, multipole_fringe, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(6,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("PolynomA"));
        mxSetCell(plhs[0],2,mxCreateString("PolynomB"));
        mxSetCell(plhs[0],3,mxCreateString("MaxOrder"));
        mxSetCell(plhs[0],4,mxCreateString("NumIntSteps"));
        mxSetCell(plhs[0],5,mxCreateString("Type"));
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(7,1);
            mxSetCell(plhs[1],0,mxCreateString("MultipoleFringe"));
            mxSetCell(plhs[1],1,mxCreateString("BendingAngle"));
            mxSetCell(plhs[1],2,mxCreateString("gK"));
            mxSetCell(plhs[1],3,mxCreateString("T1"));
            mxSetCell(plhs[1],4,mxCreateString("T2"));
            mxSetCell(plhs[1],5,mxCreateString("R1"));
            mxSetCell(plhs[1],6,mxCreateString("R2"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*MATLAB_MEX_FILE*/
