/* ThinMPolePass
   Accelerator Toolbox 
   A.Terebilo
*/

#include "atelem.c"
#include "atlalib.c"
#include "kick_kn.h"  /* kick */

struct elem
{
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    /* Optional fields */
    double Scaling;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *RApertures;
    double *EApertures;
};

void ThinMPolePass(double *r, double *A, double *B, int max_order,
        double *T1, double *T2,
        double *R1, double *R2,
        double *RApertures, double *EApertures,
        double scaling, int num_particles)
{

    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(none) \
    shared(r,num_particles,A,B,max_order,T1,T2,R1,R2,EApertures,RApertures,scaling)
    for (int c = 0; c<num_particles; c++) { /* Loop over particles */
        double *r6 = r + 6*c;
        if (!atIsNaN(r6[0])) {
            /* Check for change of reference momentum */
            if (scaling != 1.0) ATChangePRef(r6, scaling);
            /*  misalignment at entrance  */
            if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            kick(r6, A, B, max_order, 1.0);
            /*  misalignment at exit */
            if (R2) ATmultmv(r6,R2);
            if (T2) ATaddvv(r6,T2);
            /* Check for change of reference momentum */
            if (scaling != 1.0) ATChangePRef(r6, 1.0/scaling);
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double Scaling;
        int MaxOrder;
        double *PolynomA, *PolynomB, *R1, *R2, *T1, *T2, *EApertures, *RApertures;
        int nl, nc;
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        /*optional fields*/
        Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();

        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->PolynomA=PolynomA;
        Elem->PolynomB=PolynomB;
        Elem->MaxOrder=MaxOrder;
        /*optional fields*/
        Elem->Scaling=Scaling;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
    }
    ThinMPolePass(r_in, Elem->PolynomA, Elem->PolynomB, Elem->MaxOrder,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2,
            Elem->RApertures, Elem->EApertures,
            Elem->Scaling, num_particles);
    return Elem;
}

MODULE_DEF(ThinMPolePass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#ifdef MATLAB_MEX_FILE
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Scaling;
        int MaxOrder;
        double *PolynomA, *PolynomB, *BendingAngle, *R1, *R2, *T1, *T2, *EApertures, *RApertures;
        int nl, nc;
        if (mxGetM(prhs[1]) != 6) mexErrMsgTxt("Second argument must be a 6 x N matrix");

        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        /*optional fields*/
        Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();

        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        ThinMPolePass(r_in, PolynomA, PolynomB, MaxOrder,
            T1, T2, R1, R2, RApertures, EApertures,
            Scaling, num_particles);
    } else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(3,1);
        mxSetCell(plhs[0],0,mxCreateString("PolynomA"));
        mxSetCell(plhs[0],1,mxCreateString("PolynomB"));
        mxSetCell(plhs[0],2,mxCreateString("MaxOrder"));
        if (nlhs>1) {    /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(7,1);
            mxSetCell(plhs[1],0,mxCreateString("T1"));
            mxSetCell(plhs[1],1,mxCreateString("T2"));
            mxSetCell(plhs[1],2,mxCreateString("R1"));
            mxSetCell(plhs[1],3,mxCreateString("R2"));
            mxSetCell(plhs[1],4,mxCreateString("RApertures"));
            mxSetCell(plhs[1],5,mxCreateString("EApertures"));
            mxSetCell(plhs[1],6,mxCreateString("FieldScaling"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
