#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "elempass.h"
#include "atlalib.c"


/******************************************************************************/
/* PHYSICS SECTION ************************************************************/

static void matfoc(double L, double Knorm, double *MD, double *M12, double *M21)
{
   double g  = fabs(Knorm);
   double t  = sqrt(g);
   double lt = L*t;
				
   if (Knorm == 0) {
      *MD = 1.0;
      *M12 = L;
      *M21 = 0;
   }
   else if (Knorm > 0) {
	*MD = cos(lt);
	*M12 = sin(lt)/t;
	*M21 = -*M12*g;
   }
   else  {
	*MD = cosh(lt);
	*M12 = sinh(lt)/t;
	*M21 = *M12*g;
   }
}

static void foc6(double *r, double L, double Kx, double Kz)
{
   double M12,M21,M34,M43,MVD,MHD;  /* non-0 elements of transfer matrix */
   double p_norm = 1/(1+r[4]);
   double x   = r[0];
   double xpr = r[1]*p_norm;
   double y   = r[2];
   double ypr = r[3]*p_norm;

   (void) matfoc(L, Kx*p_norm, &MHD, &M12, &M21);
   (void) matfoc(L, Kz*p_norm, &MVD, &M34, &M43);

   r[0]=  MHD*x + M12*xpr;
   r[1]= (M21*x + MHD*xpr)/p_norm;
   r[2]=  MVD*y + M34*ypr;
   r[3]= (M43*y + MVD*ypr)/p_norm;

    /* no change in r[4] (delta) */

   r[5]+= (fabs(Kx)*p_norm*x*x*(L-MHD*M12) - fabs(Kz)*p_norm*y*y*(L-MVD*M34))/4;
   r[5]+= (xpr*xpr*(L+MHD*M12)+ypr*ypr*(L+MVD*M34))/4;
   r[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2;
}

void WiggLinearPass(double *r, double le, double invrho, double kxkz, double *T1, double *T2, double *R1, double *R2, int num_particles)
{
   int c;
   double *r6;
   double kz = 0.5/(1.0+kxkz)*invrho*invrho;
   double kx = kxkz*kz;

   for(c = 0;c<num_particles;c++) {
      r6 = r+c*6;
      if(!mxIsNaN(r6[0]) & mxIsFinite(r6[4])) {
      /* function quad6 internally calculates the square root
         of the energy deviation of the particle 
         To protect against DOMAIN and OVERFLOW error, check if the
         fifth component of the phase spacevector r6[4] is finite
       */
	  /* Misalignment at entrance */
	 if (T1 != NULL) ATaddvv(r6,T1);
	 if (R1 != NULL) ATmultmv(r6,R1);

	 foc6(r6, le, kx, kz);

	 /* Misalignment at exit */	
	 if (T2 != NULL) ATmultmv(r6,R2);
	 if (R2 != NULL) ATaddvv(r6,T2);
      }
   }
}

/********** END PHYSICS SECTION ***********************************************/
/******************************************************************************/
#ifndef NOMEX
/********** WINDOWS DLL GATEWAY SECTION ***************************************/


static int ReqFieldNumber(const mxArray *ElemData, const char *fname)
{
   int fnum = mxGetFieldNumber(ElemData, fname);
   if (fnum < 0) mexErrMsgIdAndTxt("At:MissingField","Required field '%s' was not found in the element data structure", fname);
   return fnum;
}

ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
				double *r_in, int num_particles, int mode)

#define NUM_FIELDS_2_REMEMBER 7

{
   double *pr1, *pr2, *pt1, *pt2 , le, invrho, kxkz;   
   int fnum, inum;

   switch(mode) {
      case NO_LOCAL_COPY:	/* Obsolete in AT1.3 et fields by names from MATLAB workspace  */
         break;	
				
      case MAKE_LOCAL_COPY: 	/* Find field numbers first
				Save a list of field number in an array
				and make returnptr point to that array
				*/
	 FieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
	 inum = 0;
	 
	 FieldNumbers[inum++] = ReqFieldNumber(ElemData, "Length");
	 FieldNumbers[inum++] = ReqFieldNumber(ElemData, "InvRho");

	 FieldNumbers[inum++] = mxGetFieldNumber(ElemData,"KxKz");
	 FieldNumbers[inum++] = mxGetFieldNumber(ElemData,"R1");
	 FieldNumbers[inum++] = mxGetFieldNumber(ElemData,"R2");
	 FieldNumbers[inum++] = mxGetFieldNumber(ElemData,"T1");
	 FieldNumbers[inum++] = mxGetFieldNumber(ElemData,"T2");

      case USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
				The second argument ponter to the array of field 
				numbers is previously created with 
				QuadLinPass( ..., MAKE_LOCAL_COPY)
				*/
	 inum = 0;
	 
	 le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[inum++]));
	 invrho = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[inum++]));

	 /* Optional fields */
	 kxkz = ((fnum=FieldNumbers[inum++]) >= 0) ? mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[fnum])) : 0.0;
	 pr1 = ((fnum=FieldNumbers[inum++]) >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[fnum])) : NULL;
	 pr2 = ((fnum=FieldNumbers[inum++]) >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[fnum])) : NULL;
	 pt1 = ((fnum=FieldNumbers[inum++]) >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[fnum])) : NULL;
	 pt2 = ((fnum=FieldNumbers[inum++]) >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[fnum])) : NULL;
	 break;

   }
   WiggLinearPass(r_in, le, invrho, kxkz, pt1, pt2, pr1, pr2 , num_particles);
   return FieldNumbers;	
}


/********** END WINDOWS DLL GATEWAY SECTION ***************************************/
/********** MATLAB GATEWAY  ***************************************/

static const mxArray *ReqField(const mxArray *ElemData, const char *fname)
{
   const mxArray *tmpmxptr = mxGetField(ElemData, 0, fname);
   if (tmpmxptr == NULL) mexErrMsgIdAndTxt("At:MissingField","Required field '%s was not found in the element data structure", fname);
   return tmpmxptr;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int m, n;
   double *r_in;   
   double  *pr1, *pr2, *pt1, *pt2 , le, invrho, kxkz; 
   mxArray *tmpmxptr;


   if (nrhs) {
	
	/* ALLOCATE memory for the output array of the same size as the input */
	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);

	if (m != 6) mexErrMsgTxt("Second argument must be a 6 x N matrix");
	
	/* Required Fields */
	
	le = mxGetScalar(ReqField(prhs[0], "Length"));
	invrho = mxGetScalar(ReqField(prhs[0], "InvRho"));
	
	/* Optionnal arguments */    
	kxkz = (tmpmxptr = mxGetField(prhs[0],0,"KxKz")) ? mxGetScalar(tmpmxptr) : 0.0;
	pr1 = (tmpmxptr = mxGetField(prhs[0],0,"R1")) ? mxGetPr(tmpmxptr) : NULL;
	pr2 = (tmpmxptr = mxGetField(prhs[0],0,"R2")) ? mxGetPr(tmpmxptr) : NULL;
	pt1 = (tmpmxptr = mxGetField(prhs[0],0,"T1")) ? mxGetPr(tmpmxptr) : NULL;
	pt2 = (tmpmxptr = mxGetField(prhs[0],0,"T2")) ? mxGetPr(tmpmxptr) : NULL;


	plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);
	WiggLinearPass(r_in, le, invrho, kxkz, pt1, pt2, pr1, pr2, n);
   }
   else {   /* return list of required fields */
      plhs[0] = mxCreateCellMatrix(2,1);
      mxSetCell(plhs[0],0,mxCreateString("Length"));
      mxSetCell(plhs[0],1,mxCreateString("InvRho"));

      if (nlhs>1) { /* Required and optional fields */ 
          plhs[1] = mxCreateCellMatrix(5,1);
          mxSetCell(plhs[1],0,mxCreateString("KxKz"));
	  mxSetCell(plhs[1],1,mxCreateString("T1"));
	  mxSetCell(plhs[1],2,mxCreateString("T2"));
	  mxSetCell(plhs[1],3,mxCreateString("R1"));
	  mxSetCell(plhs[1],4,mxCreateString("R2"));
      }
   }
}
#endif
