/* CrabCavityPass.c
   for Accelerator Toolbox 
   Xiaobiao Huang, created 11/1/2016
*/

#include "mex.h"
#include "elempass.h"
#include <math.h>
#define TWOPI  6.28318530717959
#define C0  	2.99792458e8 
#include <omp.h>
#include <float.h>
/*#define isnan(x) ((x) != (x))
#define isinf(x) (isnan(x-x))
*/
double grandn()
{
	const double epsilon = DBL_EPSILON;

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
	   return z1 ;

	double u1, u2;
	do
	 {
	   u1 = rand() * (1.0 / RAND_MAX);
	   u2 = rand() * (1.0 / RAND_MAX);
	 }
	while ( u1 <= epsilon );

	z0 = sqrt(-2.0 * log(u1)) * cos(TWOPI * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(TWOPI * u2);
	return z0 ;
	
}
	/*CrabCavityPass(r_in, le, vx/energy,vy/energy, freq, phi,sigphi,sigvv, n);*/
void CrabCavityPass(double *r_in, double le, double nvx, double nvy, double freq, double phi, double sigphi, double sigvv, int num_particles)
/* le - physical length
   nv - peak voltage (V) normalized to the design enegy (eV)
   r is a 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements 
*/
{	int c, c6;
	double halflength , p_norm, NormL;
	double nvxa,nvya, ddpp,omega_t,k;
	
	/* move the random number outside of the for loop, 11/29/2017  */
	double rndph, rndx, rndy;
	rndph = grandn();
	rndx = grandn();
	rndy = grandn();
	
	k = TWOPI*freq/C0;
	if(le == 0)
		{	
				
			for(c = 0;c<num_particles;c++)
			{	c6 = c*6;
			    if(!mxIsNaN(r_in[c6]))
				{
					if(sigphi>0)
						omega_t = k*r_in[c6+5] + phi + sigphi*rndph;
					else
						omega_t = k*r_in[c6+5] + phi;
					if(sigvv>0)
					{	
						nvxa = nvx*(1.0+sigvv*rndx);
						nvya = nvy*(1.0+sigvv*rndy);
					}
					else
					{
						nvxa = nvx;
						nvya = nvy;
					}
					ddpp = -(nvya*r_in[c6+2]+nvxa*r_in[c6+0])*k*cos(omega_t)/(1+r_in[c6+4]);
					
					r_in[c6+1] += nvxa*sin(omega_t);
					r_in[c6+3] += nvya*sin(omega_t);
					r_in[c6+4] += ddpp;
				}
			}
		}
	else
		{	halflength = le/2;
			
			
			for(c = 0;c<num_particles;c++)
			{	c6 = c*6;
				/*if(!mxIsNaN(r_in[c6])) */
				if(!isnan(r_in[c6]))				
				{   p_norm = 1/(1+r_in[c6+4]); 				
				    NormL  = halflength*p_norm;
				    /* Propagate through a drift equal to half cavity length */
				    r_in[c6+0]+= NormL*r_in[c6+1];
   			        r_in[c6+2]+= NormL*r_in[c6+3];
   			        r_in[c6+5]+= NormL*p_norm*(r_in[c6+1]*r_in[c6+1]+r_in[c6+3]*r_in[c6+3])/2;
				
				    /* Apply kicks */
				    if(sigphi>0)
						omega_t = k*r_in[c6+5] + phi + sigphi*rndph;
					else
						omega_t = k*r_in[c6+5] + phi;
					if(sigvv>0)
					{	
						nvxa = nvx*(1.0+sigvv*rndx);
						nvya = nvy*(1.0+sigvv*rndy);
					}
					else
					{
						nvxa = nvx;
						nvya = nvy;
					}
					ddpp = -(nvya*r_in[c6+2]+nvxa*r_in[c6+0])*k*cos(omega_t)/(1+r_in[c6+4]);
					
					r_in[c6+1] += nvxa*sin(omega_t);
					r_in[c6+3] += nvya*sin(omega_t);
					r_in[c6+4] += ddpp;
					
					
				    p_norm = 1/(1+r_in[c6+4]); 				
				    NormL  = halflength*p_norm;
				    /* Propagate through a drift equal to half cavity length */
				    r_in[c6+0]+= NormL*r_in[c6+1];
   			        r_in[c6+2]+= NormL*r_in[c6+3];
   			        r_in[c6+5]+= NormL*p_norm*(r_in[c6+1]*r_in[c6+1]+r_in[c6+3]*r_in[c6+3])/2;
                }
			}
		}

}


ExportMode int* passFunction(const mxArray *ElemData,int *FieldNumbers,
				double *r_in, int num_particles, int mode)


#define NUM_FIELDS_2_REMEMBER 5

{	double le, *volt, freq, energy, phi, vx, vy,sigphi, sigvv;
	int *returnptr;
	int *NewFieldNumbers, fnum;

    
	switch(mode)
		{	case NO_LOCAL_COPY:	/* Obsolete in AT1.3  */
				{	
				   
				}	break;	
			
			case MAKE_LOCAL_COPY: 	/* Find field numbers first 
									    Save a list of field number in an array
										 and make returnptr point to that array
								    */
				{	
					NewFieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
					
					fnum = mxGetFieldNumber(ElemData,"Length");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
					NewFieldNumbers[0] = fnum;
					
					fnum = mxGetFieldNumber(ElemData,"Voltages");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Voltages' was not found in the element data structure"); 
					NewFieldNumbers[1] = fnum;
					
					fnum = mxGetFieldNumber(ElemData,"Energy");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Energy' was not found in the element data structure"); 
					NewFieldNumbers[2] = fnum;
					
					
					fnum = mxGetFieldNumber(ElemData,"Frequency");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Frequency' was not found in the element data structure"); 
					NewFieldNumbers[3] = fnum;
                    
                    /* Optional field TimeLag */
                    fnum = mxGetFieldNumber(ElemData,"Phase");
                    NewFieldNumbers[4] = fnum;
					if(fnum<0) 
					    phi = 0;
                    else
                        phi = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					fnum = mxGetFieldNumber(ElemData,"SigPhi");
                    NewFieldNumbers[5] = fnum;
					if(fnum<0) 
					    sigphi = 0;
                    else
                        sigphi = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					fnum = mxGetFieldNumber(ElemData,"SigVV");
                    NewFieldNumbers[6] = fnum;
					if(fnum<0) 
					    sigvv = 0;
                    else
                        sigvv = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					
					le = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[0]));
					volt = mxGetPr(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[1]));
					vx = volt[0];
					vy = volt[1];
					energy = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[2]));
					freq = mxGetScalar(mxGetFieldByNumber(ElemData,0,NewFieldNumbers[3]));

					returnptr = NewFieldNumbers;
				}	break;

			case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
										 The second argument ponter to the array of field 
										 numbers is previously created with 
										 QuadLinPass( ..., MAKE_LOCAL_COPY)
								    */
											
				{	le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
					volt = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
					vx = volt[0];
					vy = volt[1];
					energy = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
					freq = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[3]));
                    /* Optional field TimeLag */
                    if(FieldNumbers[4]<0) 
					    phi = 0;
					else
                        phi = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
					
					if(FieldNumbers[5]<0) 
					    sigphi = 0;
					else
                        sigphi = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[5]));
					
					if(FieldNumbers[6]<0) 
					    sigvv = 0;
					else
                        sigvv = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[6]));
                    
					returnptr = FieldNumbers;
				}	break;

			default:
				{	mexErrMsgTxt("No match found for calling mode in function CavityPass\n");
				}
	}

	
	CrabCavityPass(r_in, le, vx/energy,vy/energy, freq, phi,sigphi,sigvv, num_particles);
	return(returnptr);
}



void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 	double *volt,freq, energy,vx,vy;
	int m,n;
	double *r_in, le, phi, sigphi, sigvv;   
	mxArray *tmpmxptr;
	
	if(nrhs)
	{

	/* ALLOCATE memory for the output array of the same size as the input */
	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if(m!=6) 
		mexErrMsgTxt("Second argument must be a 6 x N matrix");

	
	tmpmxptr=mxGetField(prhs[0],0,"Length");
	if(tmpmxptr)
		le = mxGetScalar(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 

	tmpmxptr=mxGetField(prhs[0],0,"Phase");
	if(tmpmxptr)
		phi = mxGetScalar(tmpmxptr);
	else
		phi=0;

	tmpmxptr=mxGetField(prhs[0],0,"SigPhi");
	if(tmpmxptr)
		sigphi = mxGetScalar(tmpmxptr);
	else
		sigphi=0;
	
	tmpmxptr=mxGetField(prhs[0],0,"SigVV");
	if(tmpmxptr)
		sigvv = mxGetScalar(tmpmxptr);
	else
		sigvv=0;
					    
	tmpmxptr=mxGetField(prhs[0],0,"Voltages");
	if(tmpmxptr)
	{
		volt = mxGetPr(tmpmxptr);/* two numbers in [Vx, Vy] */
		vx = volt[0];
		vy = volt[1];
	}
	else
		mexErrMsgTxt("Required field 'Voltages' was not found in the element data structure");

    tmpmxptr=mxGetField(prhs[0],0,"Energy");
	if(tmpmxptr)
		energy = mxGetScalar(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'Energy' was not found in the element data structure");
					    
	tmpmxptr=mxGetField(prhs[0],0,"Frequency");
	if(tmpmxptr)
		freq = mxGetScalar(tmpmxptr);
	else
		mexErrMsgTxt("Required field 'Frequency' was not found in the element data structure"); 
    
    plhs[0] = mxDuplicateArray(prhs[1]);
    r_in = mxGetPr(plhs[0]);
	CrabCavityPass(r_in, le, vx/energy,vy/energy, freq, phi,sigphi,sigvv, n);
    }
    else
    {   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(4,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("Voltages"));
	    mxSetCell(plhs[0],2,mxCreateString("Energy"));
	    mxSetCell(plhs[0],3,mxCreateString("Frequency"));
	    if(nlhs>1) /* optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(3,1); 
            mxSetCell(plhs[1],0,mxCreateString("Phase"));
			mxSetCell(plhs[1],1,mxCreateString("SigPhi"));
			mxSetCell(plhs[1],2,mxCreateString("SigVV"));
	    }
	}

}