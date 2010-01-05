#include "mex.h"
#include "elempass.h"
#include "atlalib.c"
#include "strdriftkick.c"

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656

/*
static void QuadFringePass(double* r, const double b2)
{
    double  u, ps;

    u = b2/(12.0*(1.0+r[4])); ps = u/(1.0+r[4]);
    r[3] /= 1.0 - 3.0*u*pow(r[2],2); r[2] -= u*pow(r[2],3);
//    if (globval.Cavity_on) x[ct_] -= ps*cube(x[y_])*x[py_];
    r[5] -= ps*pow(r[2],3)*r[3];
	r[1] /= 1.0 + 3.0*u*pow(r[1],2);
//    if (globval.Cavity_on) x[ct_] += ps*cube(x[x_])*x[px_];
    r[5] += ps*pow(r[0],3)*r[1];
	r[0] += u*pow(r[0],3); u = u*3.0; ps = ps*3.0;
    r[2] = exp(-u*pow(r[0],2))*r[2]; r[3] = exp(u*pow(r[0],2))*r[3];
    r[1] += 2.0*u*r[0]*r[2]*r[3];
//    if (globval.Cavity_on) x[ct_] -= ps*sqr(x[x_])*x[y_]*x[py_];
    r[5] -= ps*pow(r[0],2)*r[2]*r[3];
	r[0] = exp(u*pow(r[2],2))*r[0]; r[1] = exp(-u*pow(r[2],2))*r[1];
    r[3] -= 2.0*u*r[2]*r[1]*r[1];
//    if (globval.Cavity_on) x[ct_] += ps*sqr(x[y_])*x[x_]*x[px_];
	r[5] += ps*pow(r[2],2)*r[0]*r[1];
}
*/

static void QuadFringePassP(double* r, const double b2)
{
/*	x=r[0],px=r[1],y=r[2],py=r[3],delta=r[4],ct=r[5] 
	Lee-Whiting's thin lens limit formula as given in p. 390 of "Beam Dynamics..."by E. Forest */
   register double u = b2/(12.0*(1.0+r[4]));
   register double x2 = r[0]*r[0];
   register double z2 = r[2]*r[2];
   register double xz = r[0]*r[2];
   register double gx = u * (x2+3*z2) * r[0];
   register double gz = u * (z2+3*x2) * r[2];
   
   r[0]+=gx;
   r[1]+=3*u*(2*xz*r[3]-(x2+z2)*r[1]);
   r[2]-=gz;
   r[3]-=3*u*(2*xz*r[1]-(x2+z2)*r[3]);
   r[5]-=(gz*r[3] - gx*r[1])/(1+r[4]);
}

static void QuadFringePassN(double* r, const double b2)
{
/*	x=r[0],px=r[1],y=r[2],py=r[3],delta=r[4],ct=r[5] 
	Lee-Whiting's thin lens limit formula as given in p. 390 of "Beam Dynamics..."by E. Forest */
   register double u = b2/(12.0*(1.0+r[4]));
   register double x2 = r[0]*r[0];
   register double z2 = r[2]*r[2];
   register double xz = r[0]*r[2];
   register double gx = u * (x2+3*z2) * r[0];
   register double gz = u * (z2+3*x2) * r[2];
   
   r[0]-=gx;
   r[1]-=3*u*(2*xz*r[3]-(x2+z2)*r[1]);
   r[2]+=gz;
   r[3]+=3*u*(2*xz*r[1]-(x2+z2)*r[3]);
   r[5]+=(gz*r[3] - gx*r[1])/(1+r[4]);
}

void QuadMPoleFringePass(double *r, double le, const double *A, const double *B,
			 int max_order, int num_int_steps,
			 const double *T1, const double *T2,	
			 const double *R1, const double *R2, int num_particles)
{
   double *r6;
   int c, m;
   double norm, NormL1, NormL2;	
   bool useT1 = (T1 != NULL);
   bool useT2 = (T2 != NULL);
   bool useR1 = (R1 != NULL);
   bool useR2 = (R2 != NULL);
   double SL = le/num_int_steps;
   double L1 = SL*DRIFT1;
   double L2 = SL*DRIFT2;
   double K1 = SL*KICK1;
   double K2 = SL*KICK2;
   
   for(c = 0;c<num_particles;c++) {	/*Loop over particles  */
      r6 = r+c*6;
      if (!mxIsNaN(r6[0])) {
	  /*  misalignment at entrance  */
	 if (useT1) ATaddvv(r6, T1);
	 if (useR1) ATmultmv(r6, R1);
	 QuadFringePassP(r6,B[1]);
	 /*  integrator  */
	 for (m=0; m < num_int_steps; m++) { /*  Loop over slices */
	    r6 = r+c*6;	
	    norm = 1/(1+r6[4]);
	    NormL1 = L1*norm;
	    NormL2 = L2*norm;
	    fastdrift(r6, NormL1);
	    strthinkick(r6, A, B,  K1, max_order);
	    fastdrift(r6, NormL2);
	    strthinkick(r6, A, B, K2, max_order);
	    fastdrift(r6, NormL2);
	    strthinkick(r6, A, B,  K1, max_order);
	    fastdrift(r6, NormL1);	
	 }  
	 QuadFringePassN(r6,B[1]);
	 /* Misalignment at exit */
	 if (useR2) ATmultmv(r6, R2);
	 if (useT2) ATaddvv(r6, T2);
      }
   }
}

#ifndef NOMEX
ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
			     double *r_in, int num_particles, int mode)
#define NUM_FIELDS_2_REMEMBER 9
{
   int fnum;
   double *A , *B;
   double  *pr1, *pr2, *pt1, *pt2;   
   int max_order, num_int_steps;
   double le;
   
   switch(mode) {
      case NO_LOCAL_COPY:	/* NOT used in AT1.3 Get fields by names from MATLAB workspace   */
	 break;	
	 
      case MAKE_LOCAL_COPY: 	/* Find field numbers first
       Save a list of field number in an array
       and make returnptr point to that array
       */
	 /* Allocate memory for integer array of 
	  field numbers for faster futurereference
	  */
	 FieldNumbers = (int *) mxCalloc(NUM_FIELDS_2_REMEMBER, sizeof(int));
	 
	 /*  Populate */
	 
	 fnum = mxGetFieldNumber(ElemData,"PolynomA");
	 if(fnum<0) 
	    mexErrMsgTxt("Required field 'PolynomA' was not found in the element data structure"); 
	 FieldNumbers[0] = fnum;
	 
	 fnum = mxGetFieldNumber(ElemData,"PolynomB");
	 if(fnum<0) 
	    mexErrMsgTxt("Required field 'PolynomB' was not found in the element data structure"); 
	 FieldNumbers[1] = fnum;
	 
	 fnum = mxGetFieldNumber(ElemData,"MaxOrder");
	 if (fnum<0) mexErrMsgTxt("Required field 'MaxOrder' was not found in the element data structure"); 
	 FieldNumbers[2] = fnum;
	 
	 fnum = mxGetFieldNumber(ElemData,"NumIntSteps");
	 if (fnum<0) mexErrMsgTxt("Required field 'NumIntSteps' was not found in the element data structure"); 
	 FieldNumbers[3] = fnum;
	 
	 fnum = mxGetFieldNumber(ElemData,"Length");
	 if (fnum<0) mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
	 FieldNumbers[4] = fnum;
	 
	 FieldNumbers[5] = mxGetFieldNumber(ElemData,"R1");
	 FieldNumbers[6] = mxGetFieldNumber(ElemData,"R2");
	 FieldNumbers[7] = mxGetFieldNumber(ElemData,"T1");
	 FieldNumbers[8] = mxGetFieldNumber(ElemData,"T2");
	 /* Fall through next section... */
	 
      case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
       The second argument pointer to the array of field 
       numbers is previously created with 
       QuadMPoleFringePass(..., MAKE_LOCAL_COPY)
       */
	 A = mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[0]));
	 B = mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[1]));
	 max_order = (int)mxGetScalar(mxGetFieldByNumber(ElemData, 0, FieldNumbers[2]));
	 num_int_steps = (int)mxGetScalar(mxGetFieldByNumber(ElemData, 0, FieldNumbers[3]));
	 le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
	 
	 /* Optional fields */
	 pr1 = (FieldNumbers[5] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[5])) : NULL;
	 pr2 = (FieldNumbers[6] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[6])) : NULL;
	 pt1 = (FieldNumbers[7] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[7])) : NULL;
	 pt2 = (FieldNumbers[8] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[8])) : NULL;
	 
	 break;

      default:
	 mexErrMsgTxt("No match for calling mode in function QuadMPoleFringePass\n");

   }
   
   QuadMPoleFringePass(r_in, le, A, B, max_order, num_int_steps,
		       pt1, pt2, pr1, pr2, num_particles);
   return FieldNumbers;
}

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int m,n;
   double *r_in;
   double le, *A, *B, *pr1, *pr2, *pt1, *pt2;  
   int max_order, num_int_steps;
   mxArray *tmpmxptr;
   
   if (nrhs) {
      /* ALLOCATE memory for the output array of the same size as the input  */
      m = mxGetM(prhs[1]);
      n = mxGetN(prhs[1]);
      if(m!=6) 
	 mexErrMsgTxt("Second argument must be a 6 x N matrix");
      
      tmpmxptr =mxGetField(prhs[0],0,"PolynomA");
      if(tmpmxptr)
	 A = mxGetPr(tmpmxptr);
      else
	 mexErrMsgTxt("Required field 'PolynomA' was not found in the element data structure"); 
      
      tmpmxptr =mxGetField(prhs[0],0,"PolynomB");
      if(tmpmxptr)   
	 B = mxGetPr(tmpmxptr);
      else
	 mexErrMsgTxt("Required field 'PolynomB' was not found in the element data structure");
      
      tmpmxptr = mxGetField(prhs[0],0,"MaxOrder");
      if(tmpmxptr)
	 max_order = (int)mxGetScalar(tmpmxptr);
      else
	 mexErrMsgTxt("Required field 'MaxOrder' was not found in the element data structure");
      
      tmpmxptr = mxGetField(prhs[0],0,"NumIntSteps");
      if(tmpmxptr)   
	 num_int_steps = (int)mxGetScalar(tmpmxptr);
      else
	 mexErrMsgTxt("Required field 'NumIntSteps' was not found in the element data structure");    
      
      tmpmxptr = mxGetField(prhs[0],0,"Length");
      if(tmpmxptr)
	 le = mxGetScalar(tmpmxptr);
      else
	 mexErrMsgTxt("Required field 'Length' was not found in the element data structure");    
      
      /* Optional arguments */    
      tmpmxptr = mxGetField(prhs[0],0,"R1");
      pr1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
      
      tmpmxptr = mxGetField(prhs[0],0,"R2");
      pr2 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
      
      tmpmxptr = mxGetField(prhs[0],0,"T1");
      pt1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
      
      tmpmxptr = mxGetField(prhs[0],0,"T2");
      pt2 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
      
      plhs[0] = mxDuplicateArray(prhs[1]);
      r_in = mxGetPr(plhs[0]);
      QuadMPoleFringePass(r_in, le, A, B, max_order, num_int_steps,
			  pt1, pt2, pr1, pr2, n);
   }
   else {
      /* return list of required fields */
      plhs[0] = mxCreateCellMatrix(5,1);
      mxSetCell(plhs[0],0,mxCreateString("Length"));
      mxSetCell(plhs[0],1,mxCreateString("PolynomA"));
      mxSetCell(plhs[0],2,mxCreateString("PolynomB"));
      mxSetCell(plhs[0],3,mxCreateString("MaxOrder"));
      mxSetCell(plhs[0],4,mxCreateString("NumIntSteps"));	    
      
      if (nlhs>1) { /* Required and optional fields */
	 plhs[1] = mxCreateCellMatrix(4,1);
	 mxSetCell(plhs[1],0,mxCreateString("T1"));
	 mxSetCell(plhs[1],1,mxCreateString("T2"));
	 mxSetCell(plhs[1],2,mxCreateString("R1"));
	 mxSetCell(plhs[1],3,mxCreateString("R2"));
      }
   }
}
#endif
