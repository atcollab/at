/* MatrixTijkPass.c 
   Accelerator Toolbox 
*/

#include "mex.h"
#include "elempass.h"
#include "atlalib.c"

static void ATmultTijk(double *r, const double* T)
/*	multiplies 6-component column vector r by 6x6x6 tensor T:
 * as in r_i=Sum_j(Sum_k(Tijk*r_j*r_k)) 
  The result is stored in the memory area of r !!!
*/

{   int i,j,k;
    double temp[6];
    
    for(i=0;i<6;i++)
    {	temp[i]=0;
        for(j=0;j<6;j++)
        {
            for(k=0;k<6;k++)
                temp[i]+=T[i+j*6+k*36]*r[j]*r[k];
        }
    }
  	for(i=0;i<6;i++)
	r[i]+=temp[i];
} 

void MatrixTijkPass(double *r, const double *M, const double *Tijk,
        const double *T1, const double *T2,
        const double *R1, const double *R2, int num_particles)

{
    double *r6;
    int c;
    
    for (c = 0; c<num_particles; c++) {	/*Loop over particles  */
        r6 = r+c*6;
        if (!mxIsNaN(r6[0])) {
            if (T1 != NULL) ATaddvv(r6, T1);
            if (R1 != NULL) ATmultmv(r6, R1);
            ATmultmv(r6, M);
            ATmultTijk(r6,Tijk);
            if (R2 != NULL) ATmultmv(r6, R2);
            if (T2 != NULL) ATaddvv(r6, T2);
        }
	}
}

#ifndef NOMEX

#include "mxutils.c"

ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
        double *r_in, int num_particles, int mode)

#define NUM_FIELDS_2_REMEMBER 5
{
    double *M, *T;
    double  *pr1, *pr2, *pt1, *pt2;
    
    switch(mode) {
        case MAKE_LOCAL_COPY: 	/* Find field numbers first
         * Save a list of field number in an array
         * and make returnptr point to that array
         */
            FieldNumbers = (int *) mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
            FieldNumbers[0] = GetRequiredFieldNumber(ElemData, "M66");
            FieldNumbers[1] = GetRequiredFieldNumber(ElemData, "Tijk");
            
            
            FieldNumbers[2] = mxGetFieldNumber(ElemData,"R1");
            FieldNumbers[3] = mxGetFieldNumber(ElemData,"R2");
            FieldNumbers[4] = mxGetFieldNumber(ElemData,"T1");
            FieldNumbers[5] = mxGetFieldNumber(ElemData,"T2");
            /* Fall through next section... */
            
        case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
         * The second argument ponter to the array of field
         * numbers is previously created with
         * QuadLinPass( ..., MAKE_LOCAL_COPY)
         */
            M = mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[0]));
            T = mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[1]));
            
            /* Optional fields */
            pr1 = (FieldNumbers[2] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[2])) : NULL;
            pr2 = (FieldNumbers[3] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[3])) : NULL;
            pt1 = (FieldNumbers[4] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[4])) : NULL;
            pt2 = (FieldNumbers[5] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[5])) : NULL;
            break;
            
        default:
            mexErrMsgTxt("No match for calling mode in function MatrixTijkPass\n");
            
    }
    MatrixTijkPass(r_in, M, T, pt1, pt2, pr1, pr2, num_particles);
    return FieldNumbers;
}


void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        double  *pr1, *pr2, *pt1, *pt2;
        mxArray *tmpmxptr;
        
        double *M = mxGetPr(GetRequiredField(prhs[0], "M66"));
        double *T = mxGetPr(GetRequiredField(prhs[0], "Tijk"));
        
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        
        /* Optional arguments */
        tmpmxptr = mxGetField(prhs[0],0,"R1");
        pr1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"R2");
        pr2 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"T1");
        pt1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"T2");
        pt2 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
 		
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        MatrixTijkPass(r_in, M, T, pt1, pt2, pr1, pr2, num_particles);	
	}
    else if (nrhs == 0) {
        /* list of required fields */
	    plhs[0] = mxCreateCellMatrix(1,1);
	    mxSetCell(plhs[0],0,mxCreateString("M66"));
	    mxSetCell(plhs[0],1,mxCreateString("Tijk"));
	    if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(4,1);
            mxSetCell(plhs[1],0,mxCreateString("T1"));
            mxSetCell(plhs[1],1,mxCreateString("T2"));
            mxSetCell(plhs[1],2,mxCreateString("R1"));
            mxSetCell(plhs[1],3,mxCreateString("R2"));
	    }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif
