/* SolenoidLinearPass.c 
   Accelerator Toolbox 
   Revision 7/17/03
   A.Terebilo terebilo@slac.stanford.edu
*/

#include "atelem.c"
#include "atlalib.c"

struct elem {
    double Length;
    double ks;
    /* Optional fields */
    double *R1;
    double *R2;
    double *T1;
    double *T2;
};

void SolenoidLinearPass(double *r_in, double le, double ks, double *T1, double *T2, double *R1, double *R2, int num_particles)
/* Constant field hard edge model is assumed
   le - physical length
   ks - solenoid strength Bo / (B*rho)
   r_in - 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements 
*/
{	int c;
	double *r6, p_norm, H, S, C, x, xpr, y, ypr, NormL;
	
	if (ks!=0)
	    for (c = 0;c<num_particles;c++) {
		    r6 = r_in+c*6;
			p_norm = 1/(1+r6[4]); 

            /* Misalignment at entrance */
	        if (T1) ATaddvv(r6,T1);
			if (R1) ATmultmv(r6,R1);
			    
			x   = r6[0];
	        xpr = r6[1]*p_norm;
	        y   = r6[2];
	        ypr = r6[3]*p_norm;
			H   = ks*p_norm/2;
			S  = sin(le*H);
			C  = cos(le*H);
			
			r6[0]=     x*C*C         +  xpr*C*S/H       +   y*C*S         +   ypr*S*S/H;
	        r6[1]= ( - x*H*C*S       +  xpr*C*C         -   y*H*S*S       +   ypr*C*S         )/p_norm; 
            r6[2]=   - x*C*S         -  xpr*S*S/H       +   y*C*C         +   ypr*C*S/H;
	        r6[3]= (   x*H*S*S       -  xpr*C*S         -   y*C*S*H       +   ypr*C*C         )/p_norm;
   			r6[5]+= le*(H*H*(x*x+y*y) + 2*H*(xpr*y-ypr*x) +xpr*xpr+ypr*ypr)/2;
			
   			/* Misalignment at exit */	
			if (R2) ATmultmv(r6,R2);
		    if (T2) ATaddvv(r6,T2);
		}
    else /* Drift */
        for (c = 0;c<num_particles;c++) {
		    r6 = r_in+c*6;
			p_norm = 1/(1+r6[4]);  
			NormL  = le*p_norm;
   			r6[0]+= NormL*r6[1];
   			r6[2]+= NormL*r6[3];
   			r6[5]+= NormL*p_norm*(r6[1]*r6[1]+r6[3]*r6[3])/2;
		}
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double Length, ks;
        double *R1, *R2, *T1, *T2;
        Length=atGetDouble(ElemData,"Length"); check_error();
        ks=atGetDouble(ElemData,"K"); check_error();
        /*optional fields*/
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->ks=ks;
        /*optional fields*/
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
    }
	SolenoidLinearPass(r_in, Elem->Length, Elem->ks,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2, num_particles);
    return Elem;
}

MODULE_DEF(SolenoidLinearPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#ifdef MATLAB_MEX_FILE
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2 ) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Length, ks;
        double *R1, *R2, *T1, *T2;
        Length=atGetDouble(ElemData,"Length"); check_error();
        ks=atGetDouble(ElemData,"K"); check_error();
        /*optional fields*/
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
	    SolenoidLinearPass(r_in, Length, ks, T1, T2, R1, R2, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(2,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("K"));
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
#endif /*MATLAB_MEX_FILE*/
