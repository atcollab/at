/* VariableThinMPolePass
   Accelerator Toolbox 
   S.White simon.white@esrf.fr
*/

#include "atelem.c"
#include "atlalib.c"
#include "driftkick.c"

#define TWOPI  6.28318530717959

struct elem 
{
    double *PolynomA;
    double *PolynomB;
    double *AmplitudeA;
    double *AmplitudeB;
    double Frequency;
    int MaxOrder;
};

void VariableThinMPolePass(double *r, double *A, double *B, int max_order,
                           int num_particles)
{
    int c;
    double *r6;
    for (c = 0;c<num_particles;c++){
        r6 = r+c*6;
        if (!atIsNaN(r6[0])) {
            strthinkick(r6, A, B, 1.0, max_order);
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
			      double *r_in, int num_particles, struct parameters *Param)
{   
    if (!Elem) {
        int MaxOrder;
        double *PolynomA, *PolynomB, *AmplitudeA, *AmplitudeB;
        double FrequencyA, FrequencyB;
        double PhaseA, PhaseB;
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        AmplitudeA=atGetDoubleArray(ElemData,"AmplitudeA"); check_error();
        AmplitudeB=atGetDoubleArray(ElemData,"AmplitudeB"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();        
        FrequencyA=atGetDouble(ElemData,"FrequencyA"); check_error();
        FrequencyB=atGetDouble(ElemData,"FrequencyB"); check_error();
        PhaseA=atGetDouble(ElemData,"PhaseA"); check_error();
        PhaseB=atGetDouble(ElemData,"PhaseB"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->PolynomA=PolynomA;
        Elem->PolynomB=PolynomB;
        Elem->AmplitudeA=AmplitudeA;
        Elem->AmplitudeB=AmplitudeB;        
        Elem->MaxOrder=MaxOrder;
        Elem->FrequencyA=FrequencyA;
        Elem->FrequencyB=FrequencyB;
    }
    int i;
    double t=Param->T0*Param->nturn;
    for(i=0;i<Elem->MaxOrder;i++){
        Elem->PolynomA[i] = Elem->AmplitudeA[i]*sin(TWOPI*Elem->FrequencyA*t+PhaseA);
        Elem->PolynomB[i] = Elem->AmplitudeB[i]*sin(TWOPI*Elem->FrequencyB*t+PhaseB);
    };    
    VariableThinMPolePass(r_in, Elem->PolynomA, Elem->PolynomB, Elem->MaxOrder, num_particles);
    return Elem;
}

MODULE_DEF(VariableThinMPolePass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#ifdef MATLAB_MEX_FILE
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2 ) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        int MaxOrder;
        double *PolynomA, *PolynomB;
        double *AmplitudeA, *AmplitudeB;
        double frequency;
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        AmplitudeA=atGetDoubleArray(ElemData,"AmplitudeA"); check_error();
        AmplitudeB=atGetDoubleArray(ElemData,"AmplitudeB"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();        
        frequency=atGetDouble(ElemData,"Frequency"); check_error();
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        VariableThinMPolePass(r_in, PolynomA, PolynomB, MaxOrder,
                              num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(6,1);
        mxSetCell(plhs[0],0,mxCreateString("MaxOrder"));
        mxSetCell(plhs[0],1,mxCreateString("AmplitudeA"));
        mxSetCell(plhs[0],2,mxCreateString("AmplitudeB"));
        mxSetCell(plhs[0],3,mxCreateString("PolynomA"));
        mxSetCell(plhs[0],4,mxCreateString("PolynomB"));
        mxSetCell(plhs[0],5,mxCreateString("Frequency"));
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
