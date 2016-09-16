#include "at.h"
#include "atlalib.c"

#define TWOPI		6.28318530717959

void DeltaQPass(double *r_in, int num_particles, double alphax, double alphay,
        double betax, double betay, double qpx, double qpy,
        double a1,double a2, double a3,
        const double *T1, const double *T2,
        const double *R1, const double *R2)
{
     /*
     r_in - 6-by-N matrix of initial conditions reshaped into
     1-d array of 6*N elements
     */
    int i;
    double *rtmp;
    double x,xp,y,yp,dpp;
    double jx,jy;
    double gx,gy;
    double dqx,dqy;
    double r11x,r12x,r21x,r22x;
    double r11y,r12y,r21y,r22y;
    bool useT1 = (T1 != NULL);
    bool useT2 = (T2 != NULL);
    bool useR1 = (R1 != NULL);
    bool useR2 = (R2 != NULL);
    double cxy, sxy;
    
    gx = (1+alphax*alphax)/betax;
    gy = (1+alphay*alphay)/betay;
    
    for(i=0; i<num_particles; i++) {
        rtmp = r_in+i*6;
        if(!mxIsNaN(rtmp[0])) {
            /*  misalignment at entrance  */
            if (useT1) ATaddvv(rtmp, T1);
            if (useR1) ATmultmv(rtmp, R1);
            dpp = rtmp[4];
            x = rtmp[0];
            xp = rtmp[1]/(1.0+dpp);
            y = rtmp[2];
            yp = rtmp[3]/(1.0+dpp);
            
            jx = 0.5*(gx*x*x+2.0*alphax*x*xp+betax*xp*xp);
            jy = 0.5*(gy*y*y+2.0*alphay*y*yp+betay*yp*yp);
            
            dqx = qpx*dpp+a1*jx+a2*jy;
            dqy = qpy*dpp+a2*jx+a3*jy;
            
            cxy = cos(TWOPI*dqx);
            sxy = sin(TWOPI*dqx);
            r11x = cxy+alphax*sxy;
            r12x = betax*sxy;
            r21x = -gx*sxy;
            r22x = cxy-alphax*sxy;
            
            cxy = cos(TWOPI*dqy);
            sxy = sin(TWOPI*dqy);
            r11y = cxy+alphay*sxy;
            r12y = betay*sxy;
            r21y = -gy*sxy;
            r22y = cxy-alphay*sxy;
            
            r_in[i*6] = r11x*x+r12x*xp;
            r_in[i*6+1] = (r21x*x+r22x*xp)*(1+dpp);
            r_in[i*6+2] = r11y*y+r12y*yp;
            r_in[i*6+3] = (r21y*y+r22y*yp)*(1+dpp);
            /* Misalignment at exit */
            if (useR2) ATmultmv(rtmp, R2);
            if (useT2) ATaddvv(rtmp, T2);
        }    
    }
}

#ifdef MATLAB_MEX_FILE

#include "elempass.h"
#include "mxutils.c"

ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
                             double *r_in, int num_particles, int mode)
#define NUM_FIELDS_2_REMEMBER 13
{
	double alphax;
    double alphay;
    double betax;
    double betay;
    double qpx;
    double qpy;
    double a1;
    double a2;
    double a3;
    double  *pr1, *pr2, *pt1, *pt2;

    
	switch(mode) {
	    case MAKE_LOCAL_COPY: 	/* Find field numbers first
                                 Save a list of field number in an array
                                 and make returnptr point to that array
                                 */
            FieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
            
            /*  Populate*/
            
            FieldNumbers[0] = GetRequiredFieldNumber(ElemData, "Alphax");
            FieldNumbers[1] = GetRequiredFieldNumber(ElemData, "Alphay");
            FieldNumbers[2] = GetRequiredFieldNumber(ElemData, "Betax");
            FieldNumbers[3] = GetRequiredFieldNumber(ElemData, "Betay");
            FieldNumbers[4] = GetRequiredFieldNumber(ElemData, "Qpx");
            FieldNumbers[5] = GetRequiredFieldNumber(ElemData, "Qpy");
            FieldNumbers[6] = GetRequiredFieldNumber(ElemData, "A1");
            FieldNumbers[7] = GetRequiredFieldNumber(ElemData, "A2");
            FieldNumbers[8] = GetRequiredFieldNumber(ElemData, "A3");
            
            FieldNumbers[9] = mxGetFieldNumber(ElemData,"R1");
            FieldNumbers[10] = mxGetFieldNumber(ElemData,"R2");
            FieldNumbers[11] = mxGetFieldNumber(ElemData,"T1");
            FieldNumbers[12] = mxGetFieldNumber(ElemData,"T2");
            /* Fall through next section... */
            
	    case USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
                                 The second argument ponter to the array of field
                                 numbers is previously created with
                                 QuadLinPass( ..., MAKE_LOCAL_COPY)
                                 */
            alphax = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
            alphay = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
            betax = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
            betay = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[3]));
            qpx = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
            qpy = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[5]));
            a1 = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[6]));
            a2 = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[7]));
            a3 = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[8]));
            
            /* Optional fields */
            pr1 = (FieldNumbers[9] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[9])) : NULL;
            pr2 = (FieldNumbers[10] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[10])) : NULL;
            pt1 = (FieldNumbers[11] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[11])) : NULL;
            pt2 = (FieldNumbers[12] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[12])) : NULL;
            break;
            
	    default:
            mexErrMsgTxt("No match for calling mode in function DeltaQPass\n");
	}
       
	DeltaQPass(r_in, num_particles, alphax, alphay, betax, betay, qpx, qpy,
            a1, a2, a3, pt1, pt2, pr1, pr2);
	return FieldNumbers;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        
        double *r_in;
        double  *pr1, *pr2, *pt1, *pt2;
        mxArray *tmpmxptr;

        double alphax = mxGetScalar(GetRequiredField(prhs[0], "Alphax"));
        double alphay = mxGetScalar(GetRequiredField(prhs[0], "Alphay"));
        double betax = mxGetScalar(GetRequiredField(prhs[0], "Betax"));
        double betay = mxGetScalar(GetRequiredField(prhs[0], "Betay"));
        double qpx = mxGetScalar(GetRequiredField(prhs[0], "Qpx"));
        double qpy = mxGetScalar(GetRequiredField(prhs[0], "Qpy"));
        double a1 = mxGetScalar(GetRequiredField(prhs[0], "A1"));
        double a2 = mxGetScalar(GetRequiredField(prhs[0], "A2"));
        double a3 = mxGetScalar(GetRequiredField(prhs[0], "A3"));
        
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix: particle array");
        
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
        DeltaQPass(r_in, num_particles, alphax, alphay, betax, betay, qpx, qpy,
                a1, a2, a3, pt1, pt2, pr1, pr2);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(9,1);
        mxSetCell(plhs[0],0,mxCreateString("Alphax"));
        mxSetCell(plhs[0],1,mxCreateString("Alphay"));
        mxSetCell(plhs[0],2,mxCreateString("Betax"));
        mxSetCell(plhs[0],3,mxCreateString("Betay"));
        mxSetCell(plhs[0],4,mxCreateString("Qpx"));
        mxSetCell(plhs[0],5,mxCreateString("Qpy"));
        mxSetCell(plhs[0],6,mxCreateString("A1"));
        mxSetCell(plhs[0],7,mxCreateString("A2"));
        mxSetCell(plhs[0],8,mxCreateString("A3"));
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
