/* GWigSymplecticPass.c for
   Accelerator Toolbox 
*/

/*
 *---------------------------------------------------------------------------
 * Modification Log:
 * -----------------
 * .02  2003-06-18     J. Li, jing@fel.duke.edu
 *				Cleanup the code
 *
 * .01  2003-04-20     YK Wu, wu@fel.duke.edu
 *				GWiggler interface
 *
 *---------------------------------------------------------------------------
 *  Accelerator Physics Group, Duke FEL Lab, www.fel.duke.edu  
 */

#include "mex.h"
#include "elempass.h"
#include "atlalib.c"
#include <stdlib.h>
#include <math.h>
#include "gwig.h"

#define GWIG

#include "gwig.c"

/******************************************************************************/
/* PHYSICS SECTION ************************************************************/

void GWigInit(struct gwig *Wig, double Ltot, double Lw, double Bmax,
	      int Nstep, int Nmeth, int NHharm, int NVharm,
	      double *pBy, double *pBx, double *T1, double *T2, 
	      double *R1, double *R2)
{
  double *tmppr;
  double design_energy;
  int    i;
  double kw;
  mxArray *E0Field;
  const mxArray *GLOBVALPTR = mexGetVariablePtr("GLOBVAL","global");

  /* Get the design energy for normalization from GLOBVAL.E0 in
   * global workspace 
   */
  if(GLOBVALPTR == NULL)
    mexPrintf("GLOBVALPTR = NULL\n");
  if(GLOBVALPTR != NULL){
    E0Field = mxGetField(GLOBVALPTR,0,"E0");
    if(E0Field != NULL)
      design_energy = mxGetScalar(E0Field)*1e-9; /* convert to GeV */
    else
      mexErrMsgTxt("Field 'E0' is not defined in GLOBVAL structure");
  }
  else
    mexErrMsgTxt("global variable GLOBVAL does not exist");

  Wig->E0 = design_energy;
  Wig->Pmethod = Nmeth;
  Wig->PN = Nstep;
  Wig->Nw = (int)(Ltot / Lw);
  Wig->NHharm = NHharm;
  Wig->NVharm = NVharm;
  Wig->PB0 = Bmax;
  Wig->Lw  = Lw;


  kw = 2.0e0*PI/(Wig->Lw);
  Wig->Zw = 0.0;
  Wig->Aw = 0.0;
  tmppr = pBy;
  for (i = 0; i < NHharm; i++){
    tmppr++;
    Wig->HCw[i] = 0.0;
    Wig->HCw_raw[i] = *tmppr;

    tmppr++;
    Wig->Hkx[i]     = (*tmppr) * kw;

    tmppr++;
    Wig->Hky[i]     = (*tmppr) * kw;

    tmppr++;
    Wig->Hkz[i]     = (*tmppr) * kw;

    tmppr++;
    Wig->Htz[i]     =  *tmppr;

    tmppr++;
  }

  tmppr = pBx;
  for (i = 0; i < NVharm; i++){
    tmppr++;
    Wig->VCw[i] = 0.0;
    Wig->VCw_raw[i] = *tmppr;

    tmppr++;
    Wig->Vkx[i]     = (*tmppr) * kw;

    tmppr++;
    Wig->Vky[i]     = (*tmppr) * kw;

    tmppr++;
    Wig->Vkz[i]     = (*tmppr) * kw;

    tmppr++;
    Wig->Vtz[i]     =  *tmppr;

    tmppr++;
  }
  
  for (i = NHharm ; i< WHmax; i++) {
    Wig->HCw[i] = 0.0;
    Wig->HCw_raw[i] = 0.0;
    Wig->Hkx[i] = 0.0;
    Wig->Hky[i] = 0.0;
    Wig->Hkz[i] = 0.0;
    Wig->Htz[i] = 0.0;
  }
  for (i = NVharm ; i< WHmax; i++) {
    Wig->VCw[i] = 0.0;
    Wig->VCw_raw[i] = 0.0;
    Wig->Vkx[i] = 0.0;
    Wig->Vky[i] = 0.0;
    Wig->Vkz[i] = 0.0;
    Wig->Vtz[i] = 0.0;
  }
}

#define second 2
#define fourth 4
void GWigSymplecticPass(double *r, double Ltot, double Lw, double Bmax,
			int Nstep, int Nmeth, int NHharm, int NVharm,
			double *pBy, double *pBx, double *T1, double *T2, 
			double *R1, double *R2, int num_particles)
{	

  int c, i;
  double *r6;
  double *tmppr;
  struct gwig Wig;

  GWigInit(&Wig, Ltot, Lw, Bmax, Nstep, Nmeth,NHharm,NVharm,pBy,pBx,T1,T2,R1,R2);

  for(c = 0;c<num_particles;c++){	
    r6 = r+c*6;	
    
    if(!mxIsNaN(r6[0]) & mxIsFinite(r6[4]))
    {
    switch (Nmeth) {
      case second :
	GWigPass_2nd(&Wig, r6);
	break;
      case fourth:
	GWigPass_4th(&Wig, r6);
	break;
      default:
	printf("Invalid method ...\n");
	break;
    }
    }  
  }

}

/********** END PHYSICS SECTION ***********************************************/
/******************************************************************************/

ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
			     double *r_in, int num_particles, int mode)

#define NUM_FIELDS_2_REMEMBER 11

{	
  double *pr1, *pr2, *pt1, *pt2;
  int m, n;
  double *pBy, *pBx;
  double Ltot, Lw, Bmax; 
  int Nstep, Nmeth;
  int NHharm, NVharm;
  mxArray *tmpmxptr;
  int *returnptr,fnum;
  int *NewFieldNumbers;

  switch(mode)
    {	
    case NO_LOCAL_COPY:	/* Get fields by names from MATLAB workspace  */
      {	
	tmpmxptr = mxGetField(ElemData,0,"Length");
	if(tmpmxptr)
	  Ltot = mxGetScalar(tmpmxptr);
	else
	  mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
	
	tmpmxptr = mxGetField(ElemData,0,"Lw");
	if(tmpmxptr)
	  Lw = mxGetScalar(tmpmxptr);
	else
	  mexErrMsgTxt("Required field 'Lw' was not found in the element data structure"); 
	
	tmpmxptr = mxGetField(ElemData,0,"Bmax");
	if(tmpmxptr)
	  Bmax = mxGetScalar(tmpmxptr);
	else
	  mexErrMsgTxt("Required field 'Bmax' was not found in the element data structure"); 
	
	tmpmxptr = mxGetField(ElemData,0,"Nstep");
	if(tmpmxptr)
	  Nstep = (int)mxGetScalar(tmpmxptr);
	else
	  mexErrMsgTxt("Required field 'Nstep' was not found in the element data structure"); 
  
	tmpmxptr = mxGetField(ElemData,0,"Nmeth");
	if(tmpmxptr)
	  Nmeth = (int)mxGetScalar(tmpmxptr);
	else
	  mexErrMsgTxt("Required field 'Nmeth' was not found in the element data structure"); 

	tmpmxptr = mxGetField(ElemData,0,"NHharm");
	if(tmpmxptr)
	  NHharm = (int)mxGetScalar(tmpmxptr);
	else
	  mexErrMsgTxt("Required field 'NHharm' was not found in the element data structure"); 
	
	tmpmxptr = mxGetField(ElemData,0,"NVharm");
	if(tmpmxptr)
	  NVharm = (int)mxGetScalar(tmpmxptr);
	else
	  mexErrMsgTxt("Required field 'NVharm' was not found in the element data structure"); 

	tmpmxptr = mxGetField(ElemData,0,"By");
	if(tmpmxptr)
	  pBy = mxGetPr(tmpmxptr);
	else
	  mexErrMsgTxt("Required field 'By' was not found in the element data structure"); 

	tmpmxptr = mxGetField(ElemData,0,"Bx");
	if(tmpmxptr)
	  pBx = mxGetPr(tmpmxptr);
	else
	  mexErrMsgTxt("Required field 'Bx' was not found in the element data structure"); 
  
	tmpmxptr = mxGetField(ElemData,0,"R1");
	if(tmpmxptr)
	  pr1 = mxGetPr(tmpmxptr);
	else
	  mexErrMsgTxt("Required field 'R1' was not found in the element data structure"); 
	
	tmpmxptr = mxGetField(ElemData,0,"R2");
	if(tmpmxptr)
	  pr2 = mxGetPr(tmpmxptr);
	else
	  mexErrMsgTxt("Required field 'R2' was not found in the element data structure"); 
	
	tmpmxptr = mxGetField(ElemData,0,"T1");
	if(tmpmxptr)
	  pt1 = mxGetPr(tmpmxptr);
	else
	  mexErrMsgTxt("Required field 'T1' was not found in the element data structure"); 
	
	tmpmxptr = mxGetField(ElemData,0,"T2");
	if(tmpmxptr)
	  pt2 = mxGetPr(tmpmxptr);
	else
	  mexErrMsgTxt("Required field 'T2' was not found in the element data structure"); 	    

	returnptr = NULL;
	
      }	break;	
    case MAKE_LOCAL_COPY: 	/* Find field numbers first 
				 * Save a list of field number in an array
				 *  and make returnptr point to that array
				 */
      {						
	NewFieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
					
	fnum = mxGetFieldNumber(ElemData,"Length");
	if(fnum<0) 
	  mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
	NewFieldNumbers[0] = fnum;
					
	fnum = mxGetFieldNumber(ElemData,"Lw");
	if(fnum<0) 
	  mexErrMsgTxt("Required field 'Lw' was not found in the element data structure"); 
	NewFieldNumbers[1] = fnum;
					

	fnum = mxGetFieldNumber(ElemData,"Bmax");
	if(fnum<0) 
	  mexErrMsgTxt("Required field 'Bmax' was not found in the element data structure"); 
	NewFieldNumbers[2] = fnum;
					

	fnum = mxGetFieldNumber(ElemData,"Nstep");
	if(fnum<0) 
	  mexErrMsgTxt("Required field 'Nstep' was not found in the element data structure"); 
	NewFieldNumbers[3] = fnum;
					
					
	fnum = mxGetFieldNumber(ElemData,"Nmeth");
	if(fnum<0) 
	  mexErrMsgTxt("Required field 'Nmeth' was not found in the element data structure"); 
	NewFieldNumbers[4] = fnum;
					
	fnum = mxGetFieldNumber(ElemData,"NHharm");
	if(fnum<0) 
	  mexErrMsgTxt("Required field 'NHharm' was not found in the element data structure"); 
	NewFieldNumbers[5] = fnum;

	fnum = mxGetFieldNumber(ElemData,"NVharm");
	if(fnum<0) 
	  mexErrMsgTxt("Required field 'NVharm' was not found in the element data structure"); 
	NewFieldNumbers[6] = fnum;       
         
	fnum = mxGetFieldNumber(ElemData,"By");
	if(fnum<0) 
	  mexErrMsgTxt("Required field 'By' was not found in the element data structure"); 
	NewFieldNumbers[7] = fnum;
	                
	fnum = mxGetFieldNumber(ElemData,"Bx");
	if(fnum<0) 
	  mexErrMsgTxt("Required field 'Bx' was not found in the element data structure"); 
	NewFieldNumbers[8] = fnum;
	                
	fnum = mxGetFieldNumber(ElemData,"R1");
	if(fnum<0) 
	  mexErrMsgTxt("Required field 'R1' was not found in the element data structure"); 
	NewFieldNumbers[9] = fnum;
	
	
	fnum = mxGetFieldNumber(ElemData,"R2");
	if(fnum<0) 
	  mexErrMsgTxt("Required field 'R2' was not found in the element data structure"); 
	NewFieldNumbers[10] = fnum;
	
	
	fnum = mxGetFieldNumber(ElemData,"T1");
	if(fnum<0) 
	  mexErrMsgTxt("Required field 'T1' was not found in the element data structure"); 
	NewFieldNumbers[11] = fnum;
	
	
	fnum = mxGetFieldNumber(ElemData,"T2");
	if(fnum<0) 
	  mexErrMsgTxt("Required field 'T2' was not found in the element data structure"); 
	NewFieldNumbers[12] = fnum;
	
	Ltot   = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[0]));
	Lw     = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[1]));
	Bmax   = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[2]));
	Nstep  = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[3]));
	Nmeth  = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[4]));
	NHharm = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[5]));
	NVharm = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[6]));
	pBy    = mxGetPr(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[7]));
	pBx    = mxGetPr(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[8]));
	pr1    = mxGetPr(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[9]));
	pr2    = mxGetPr(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[10]));
	pt1    = mxGetPr(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[11]));
	pt2    = mxGetPr(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[12]));
	
	returnptr = NewFieldNumbers;
      }	break;
    case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
				   The second argument ponter to the array of
				   field numbers is previously created with 
				   QuadLinPass( ..., MAKE_LOCAL_COPY)
				*/
      {	
	Ltot   = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
	Lw     = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
	Bmax   = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
	Nstep  = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[3]));
	Nmeth  = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
	NHharm = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[5]));
	NVharm = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[6]));
	pBy    = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[7]));
	pBx    = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[8]));
	pr1    = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[9]));
	pr2    = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[10]));
	pt1    = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[11]));
	pt2    = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[12]));
	
	returnptr = FieldNumbers;
      }	break;
    default:
      {	mexErrMsgTxt("No match for calling mode in function QuadLinPass\n");
      }
    }

  GWigSymplecticPass(r_in,Ltot,Lw,Bmax,Nstep,Nmeth,NHharm,NVharm,pBy,pBx,pt1,pt2,pr1,pr2,num_particles);
  return(returnptr);	
}


/********** END WINDOWS DLL GATEWAY SECTION ***************************************/

/********** MATLAB GATEWAY  ***************************************/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
  int m, n;
  double *r_in;
  double *pBy, *pBx;
  double *pr1, *pr2, *pt1, *pt2 ;
  double Ltot, Lw, Bmax; 
  int Nstep, Nmeth;
  int NHharm, NVharm;
  mxArray *tmpmxptr;
	
  /* ALLOCATE memory for the output array of the same size as the input */
  m = mxGetM(prhs[1]);
  n = mxGetN(prhs[1]);
  if(m!=6) 
    {mexErrMsgTxt("Second argument must be a 6 x N matrix");}	
	
  tmpmxptr = mxGetField(prhs[0],0,"Length");
  if(tmpmxptr)
    Ltot = mxGetScalar(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
	
  tmpmxptr = mxGetField(prhs[0],0,"Lw");
  if(tmpmxptr)
    Lw = mxGetScalar(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'Lw' was not found in the element data structure"); 
  
  tmpmxptr = mxGetField(prhs[0],0,"Bmax");
  if(tmpmxptr)
    Bmax = mxGetScalar(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'Bmax' was not found in the element data structure"); 
  
  tmpmxptr = mxGetField(prhs[0],0,"Nstep");
  if(tmpmxptr)
    Nstep = (int)mxGetScalar(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'Nstep' was not found in the element data structure"); 
  
  tmpmxptr = mxGetField(prhs[0],0,"Nmeth");
  if(tmpmxptr)
    Nmeth = (int)mxGetScalar(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'Nmeth' was not found in the element data structure"); 

  tmpmxptr = mxGetField(prhs[0],0,"NHharm");
  if(tmpmxptr)
    NHharm = (int)mxGetScalar(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'NHharm' was not found in the element data structure"); 
  
  tmpmxptr = mxGetField(prhs[0],0,"NVharm");
  if(tmpmxptr)
    NVharm = (int)mxGetScalar(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'NVharm' was not found in the element data structure"); 
  
  tmpmxptr = mxGetField(prhs[0],0,"By");
  if(tmpmxptr)
    pBy = mxGetPr(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'By' was not found in the element data structure"); 

  tmpmxptr = mxGetField(prhs[0],0,"Bx");
  if(tmpmxptr)
    pBx = mxGetPr(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'Bx' was not found in the element data structure"); 
  
  tmpmxptr = mxGetField(prhs[0],0,"R1");
  if(tmpmxptr)
    pr1 = mxGetPr(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'R1' was not found in the element data structure"); 
  
  tmpmxptr = mxGetField(prhs[0],0,"R2");
  if(tmpmxptr)
    pr2 = mxGetPr(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'R2' was not found in the element data structure"); 
  
  tmpmxptr = mxGetField(prhs[0],0,"T1");
  if(tmpmxptr)
    pt1 = mxGetPr(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'T1' was not found in the element data structure"); 
  
  tmpmxptr = mxGetField(prhs[0],0,"T2");
  if(tmpmxptr)
    pt2 = mxGetPr(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'T2' was not found in the element data structure"); 	    

  plhs[0] = mxDuplicateArray(prhs[1]);
  r_in = mxGetPr(plhs[0]);
  
  GWigSymplecticPass(r_in,Ltot,Lw,Bmax,Nstep,Nmeth,NHharm,NVharm,pBy,pBx,pt1,pt2,pr1,pr2,n);
}


