#include "mex.h"
extern "C"
{
#include "elempass.h"
}
#include <string.h>
#include "atlalib.c"

#define AT_MODE
#include "track.cc"

/* field metadata */

static const char * FieldNames[] =
{
    "PolynomA",
    "PolynomB",
    "MaxOrder", 
    "NumIntSteps",
    "Length",
    "BendingAngle",
    "Type",
    "gK",
    "MultipoleFringe",
    "R1", 
    "R2",
    "T1",
    "T2"
};

static const char * FirstOptionalName = "gK";

#define NUM_FIELDS_2_REMEMBER (sizeof(FieldNames) / sizeof(FieldNames[0]))

static int Lookup_Field(const char * name, const mxArray *ElemData, int optional)
{
    char buffer[1024];
    int fnum = mxGetFieldNumber(ElemData, name);
    if(fnum < 0 && optional == 0)
    {
        sprintf(buffer, 
                "Required field '%s' was not found in the element data structure", name);
        mexErrMsgTxt(buffer);
    }
    return fnum;
}

static double * Get_Pr(const mxArray *ElemData, int *FieldNumbers, int n)
{
    if(FieldNumbers[n] < 0)
        return NULL;
    else 
        return mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[n]));
}

static double Get_Scalar(const mxArray *ElemData, int *FieldNumbers, int n)
{
    if(FieldNumbers[n] < 0)
        return 0;
    else 
        return mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[n]));
}

static void Lookup_Fields(const mxArray *ElemData, int * Index)
{
    int n;
    int optional = 0;
    for(n = 0; n < NUM_FIELDS_2_REMEMBER; n++)
    {
        const char * name = FieldNames[n];
        if(strcmp(name, FirstOptionalName) == 0)
        {
            optional = 1;
        }
        Index[n] = Lookup_Field(name, ElemData, optional);
    }
}

static void BendPass(double * r, int num_particles, const mxArray * ElemData, int * FieldNumbers)
{
    int c, n, p;
    double * r6;

    int f = 0;
    
    /* oh the old rtti problem again */

    double * A           =      Get_Pr     (ElemData, FieldNumbers, f++);
    double * B           =      Get_Pr     (ElemData, FieldNumbers, f++);
    int max_order        = (int)Get_Scalar (ElemData, FieldNumbers, f++);
    int num_int_steps    = (int)Get_Scalar (ElemData, FieldNumbers, f++);
    double le            =      Get_Scalar (ElemData, FieldNumbers, f++);
    double phi           =      Get_Scalar (ElemData, FieldNumbers, f++);
    int type             = (int)Get_Scalar (ElemData, FieldNumbers, f++);

    double gK            =      Get_Scalar (ElemData, FieldNumbers, f++);
    int multipole_fringe = (int)Get_Scalar (ElemData, FieldNumbers, f++);

    double * R1          =      Get_Pr     (ElemData, FieldNumbers, f++);
    double * R2          =      Get_Pr     (ElemData, FieldNumbers, f++);
    double * T1          =      Get_Pr     (ElemData, FieldNumbers, f++);
    double * T2          =      Get_Pr     (ElemData, FieldNumbers, f++);
    
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
        r6 = r+c*6;	
        if(!mxIsNaN(r6[0]))
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

ExportMode int * passFunction(const mxArray *ElemData, int *FieldNumbers,
                              double *r_in, int num_particles, int mode)

{
    int n;
    if(mode == MAKE_LOCAL_COPY)
    {
        /* Allocate memory for integer array of 
           field numbers for faster future reference

           JR pretty uncool, what about the memory leak?
           that's why my AT crashes
        */
        FieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
        Lookup_Fields(ElemData, FieldNumbers);
     }
        
    if(mode == MAKE_LOCAL_COPY || mode == USE_LOCAL_COPY)
    {
        BendPass(r_in, num_particles, ElemData, FieldNumbers);
    }
    
    return FieldNumbers;
}

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
        BendPass(r_in, n, ElemData, FieldNumbers);
    }
    else
    {
        ShowFields(nlhs, plhs);
    }
}
