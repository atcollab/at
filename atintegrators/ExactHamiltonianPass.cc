#include "atelem.c"
#include <string.h>
#include "atlalib.c"

#define AT_MODE
#include "track.cc"

struct elem
{
    double Length;
    double BendingAngle;
    double *PolynomA;
    double *PolynomB;
    int MaxOrder;
    int NumIntSteps;
    int Type;
    /* Optional fields */
    double gK;
    int MultipoleFringe;
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

    for(n = 0; n < max_order; n++)
    {
        e.F[2 * n]     = B[n];
        e.F[2 * n + 1] = A[n];
    }
    
    e.L = le;
    e.phi = phi;
    e.gK = gK;
    e.nF = max_order;
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
    double le, bending_angle;
    double *polynom_a, *polynom_b;
    long max_order, num_int_steps, type, multipole_fringe;
    double *R1, *R2, *T1, *T2;
    double phi, gK;
    le = atGetDouble(ElemData,"Length"); check_error();
    bending_angle = atGetDouble(ElemData,"BendingAngle"); check_error();
    polynom_a = atGetDoubleArray(ElemData,"PolynomA"); check_error();
    polynom_b = atGetDoubleArray(ElemData,"PolynomB"); check_error();
    max_order = atGetLong(ElemData, "MaxOrder"); check_error();
    num_int_steps = atGetLong(ElemData, "NumIntSteps"); check_error();
    type = atGetLong(ElemData, "Type"); check_error();
    multipole_fringe = atGetLong(ElemData, "Type"); check_error();
    phi = atGetDouble(ElemData,"Phi"); check_error();
    gK = atGetDouble(ElemData,"gK"); check_error();
    R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
    R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
    T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
    T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
    ExactHamiltonianPass(r_in, le, polynom_a, polynom_b, T1, T2, R1, R2, max_order,
             num_int_steps, phi, type, gK, multipole_fringe, num_particles);
    Elem = (struct elem*)atMalloc(sizeof(struct elem));
    Elem->Length=le;
    Elem->BendingAngle=bending_angle;
    Elem->PolynomA=polynom_a;
    Elem->PolynomB=polynom_b;
    Elem->MaxOrder=max_order;
    Elem->NumIntSteps=num_int_steps;
    /*optional fields*/
    Elem->R1=R1;
    Elem->R2=R2;
    Elem->T1=T1;
    Elem->T2=T2;

    return Elem;
}

MODULE_DEF(ExactHamiltonianPass) /* Dummy module initialisation */

#endif /* defined(PYAT) || defined(MATLAB_MEX_FILE) */



#if defined(MATLAB_MEX_FILE)

/* generic show required and optional fields */
static void ShowFields(int nlhs, mxArray **plhs)
{
    int n;
    int optional = 0;
    int count_required = 0;
    int count_optional = 0;
    for(n = 0; n < NUM_FIELDS_2_REMEMBER; n++)
    {
        const char * name = FieldNames[n];
        if(strcmp(name, FirstOptionalName) == 0)
        {
            optional = 1;
        }
        if(optional)
        {
            count_optional++;
        }
        else
        {
            count_required++;
        }
    }
    plhs[0] = mxCreateCellMatrix(count_required, 1);
    for(n = 0; n < count_required; n++)
    {
        mxSetCell(plhs[0],n,mxCreateString(FieldNames[n]));
    }
    if(nlhs > 1)
    {   
        plhs[1] = mxCreateCellMatrix(count_optional, 1);
        for(n = 0; n < count_optional; n++)
        {
            mxSetCell(plhs[1],n,mxCreateString(FieldNames[n + count_required]));
        }
    }
}

/* generic mexFunction */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int FieldNumbers[NUM_FIELDS_2_REMEMBER];
    double * r_in;
    if(nrhs >= 2)
    {
	int m = mxGetM(prhs[1]);
	int n = mxGetN(prhs[1]);
        const mxArray * ElemData = prhs[0];
	if(m!=6) 
            mexErrMsgTxt("Second argument must be a 6 x N matrix");
        Lookup_Fields(ElemData, FieldNumbers);
	plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);
        le = atGetDouble(ElemData,"Length"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        max_order = atGetInt(ElemData, "MaxOrder"); check_error();
        num_int_steps = atGetInt(ElemData, "NumIntSteps"); check_error();
        phi = atGetDouble(ElemData,"Phi"); check_error();
        gK = atGetDouble(ElemData,"gK"); check_error();
        multipole_fringe = atGetLong(ElemData, "Type"); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        ExactHamiltonianPass(r_in, le, polynom_a, polynom_b, T1, T2, R1, R2, max_order,
                 num_int_steps, phi, type, gK, multipole_fringe, num_particles)
    }
    else
    {
        ShowFields(nlhs, plhs);
    }
}
#endif
