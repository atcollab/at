#include "mex.h"
#include "elempass.h"
#include "atlalib.c"

#define SQR(X) ((X)*(X))

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656

/*
 * X. Huang, 8/3/2018, modifed to include quadrupole linear fringe field and 
 * use large angle drift propagation
 */

/*
Full linear quadrupole fringe field effect, implemneting formulas in D. Zhou's paper
*/
void quadlinearfringefull(double *r, double *Iminus, double *Iplus, double K1,int flag_edge)
/*Iminus = [I0m,I1m,I2m,I3m,Lambdam], normalized by K1
 *Iplus = [I0p,I1p,I2p,I3p,Lambdap]
 *flag_order = 0 for entrance, 1 for exit
 */
{
    double p_norm,J1x,J1y,J2,J3;
    if ((Iminus!=NULL) && (Iplus!=NULL))
    {
        J1x = Iminus[1]+Iplus[1] - 2./3.*K1*Iminus[3]+0.5*Iplus[0]*(Iminus[2]+Iplus[2]);
        J1y = Iminus[1]+Iplus[1] + 2./3.*K1*Iminus[3]-0.5*Iplus[0]*(Iminus[2]+Iplus[2]);
        J2 = Iminus[2]+Iplus[2];
        J3 = K1*Iminus[2]+(Iminus[4]+Iplus[4])-Iplus[0]*(Iminus[1]+Iplus[1]);
        
        p_norm = 1.0/(1+r[4]);
        /*mexPrintf("I1m=%4.3e\tI1p=%4.3e\tJ1x=%5.4e\tJ2=%4.2e\tJ3=%4.2e\n",Iminus[1],Iplus[1],J1x,J2,J3);*/
        J1x = J1x*p_norm*K1;
        J1y = J1y*p_norm*K1;
        J2 = J2*p_norm*K1;
        J3 = J3*p_norm*K1;
        

        if (flag_edge==0)
        {
            r[1] += J3*r[0];
            r[3] += J3*r[2];
            
            r[0] += J2*r[1];
            r[2] -= J2*r[3];
         
            r[0] /= exp(J1x); 
            r[1] *= exp(J1x);
            r[2] *= exp(J1y);
            r[3] /= exp(J1y);
        }
        
        else
        {
            r[0] *= exp(J1x); 
            r[1] /= exp(J1x);
            r[2] /= exp(J1y);
            r[3] *= exp(J1y);
            
            r[0] += J2*r[1];
            r[2] -= J2*r[3];
            
            r[1] += J3*r[0];
            r[3] += J3*r[2];
        }
        
    }
    
}

/*
Hard edge fringe field effect of the third order terms
See Forest & Milutinovic, NIMA 269 (1988) 474
*/

void quadnonlinearfringe(double *r, double kv)
{
	double p_norm = 1/(1+r[4]);
	
	/* R1 reprensets a rotation of the xOy axes by -pi/4 about the s-axis */
	double R1[] ={1/sqrt(2.0), 0, 1/sqrt(2.0), 0, 0, 0,    0, 1/sqrt(2.0),  0, 1/sqrt(2.0), 0, 0,    
		-1/sqrt(2.0), 0, 1/sqrt(2.0), 0, 0, 0,   0, -1/sqrt(2.0), 0, 1/sqrt(2.0),0,0,
	0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0, 1};
	double R2[] ={1/sqrt(2.0), 0, -1/sqrt(2.0), 0, 0, 0,    0, 1/sqrt(2.0),  0, -1/sqrt(2.0), 0, 0,    
		1/sqrt(2.0), 0, 1/sqrt(2.0), 0, 0, 0,   0, 1/sqrt(2.0), 0, 1/sqrt(2.0),0,0,
	  0, 0, 0, 0, 1, 0,     0, 0, 0, 0, 0, 1};
	double x, px, y ,py;
			
	
	ATmultmv(r,R2);	 
	x   = r[0];
	px  = r[1];
	y   = r[2];
	py  = r[3];
		 
	r[0] = x - kv/6.*y*y*y*p_norm;
	r[3] = py+ kv/2.*y*y*px*p_norm;
	
	x   = r[0];
	py  = r[3];
	
	r[1] = px + kv/2.*x*x*py*p_norm;
	r[2] = y -  kv/6.*x*x*x*p_norm;
	
  ATmultmv(r,R1);	 
}

void ladrift6(double* r, double L)
/* large angle drift, X. Huang, 7/31/2018  
 * Input parameter L is the physical length
     1/(1+delta) normalization is done internally
 * Hamiltonian H = (1+\delta)-sqrt{(1+\delta)^2-p_x^2-p_y^2}, change sign for $\Delta z$ in AT
*/
{	double p_norm = 1./sqrt(SQR(1+r[4])-SQR(r[1])-SQR(r[3])); 
	double NormL  = L*p_norm;   
	r[0]+= NormL*r[1]; 
	r[2]+= NormL*r[3];
	r[5]+= L*(p_norm*(1+r[4])-1.);
}

void strthinkick(double* r, double* A, double* B, double L, int max_order)
/***************************************************************************** 
Calculate and apply a multipole kick to a 6-dimentional
phase space vector in a straight element ( quadrupole)

IMPORTANT !!!
The reference coordinate system is straight but the field expansion may still
contain dipole terms: PolynomA(1), PolynomB(1) - in MATLAB notation,
A[0], B[0] - C,C++ notation


   Note: in the US convention the transverse multipole field is written as:

                         max_order+1
                           ----
                           \                       n-1
	   (B + iB  )/ B rho  =  >   (ia  + b ) (x + iy)
         y    x            /       n    n
	                        ----
                          n=1
	is a polynomial in (x,y) with the highest order = MaxOrder
	

	Using different index notation 
   
                         max_order
                           ----
                           \                       n
	   (B + iB  )/ B rho  =  >   (iA  + B ) (x + iy)
         y    x            /       n    n
	                       ----
                          n=0

	A,B: i=0 ... max_order
   [0] - dipole, [1] - quadrupole, [2] - sextupole ...
   units for A,B[i] = 1/[m]^(i+1)
	Coeficients are stroed in the PolynomA, PolynomB field of the element
	structure in MATLAB

	A[i] (C++,C) =  PolynomA(i+1) (MATLAB) 
	B[i] (C++,C) =  PolynomB(i+1) (MATLAB) 
	i = 0 .. MaxOrder

******************************************************************************/
{  int i;
	double ReSum = B[max_order];
 	double ImSum = A[max_order];
	double ReSumTemp;
    	for(i=max_order-1;i>=0;i--)
        {   ReSumTemp = ReSum*r[0] - ImSum*r[2] + B[i];
            ImSum = ImSum*r[0] +  ReSum*r[2] + A[i];
            ReSum = ReSumTemp;
        }

    r[1] -=  L*ReSum;
    r[3] +=  L*ImSum;
}

void fastdrift(double* r, double NormL)

/*   NormL=(Physical Length)/(1+delta)  is computed externally to speed up calculations
     in the loop if momentum deviation (delta) does not change
     such as in 4-th order symplectic integrator w/o radiation
*/

{   double dx = NormL*r[1];
    double dy = NormL*r[3];
    r[0]+= dx;
    r[2]+= dy;
    r[5]+= NormL*(r[1]*r[1]+r[3]*r[3])/(2*(1+r[4]));
}


void StrMPoleSymplectic4NPass(double *r, double le, double *A, double *B,
					int max_order, int num_int_steps,double *Iminus, double *Iplus, int *flag_quadedge,
					double *T1, double *T2,	
					double *R1, double *R2, int num_particles)
/* 
 * Note Iminus and Iplus are for quadrupole fringe field
 */
{	int c,m;
	
	double *r6;
	bool useT1, useT2, useR1, useR2;
	double SL, L1, L2, K1, K2;
    double kv;
	SL = le/num_int_steps;
	L1 = SL*DRIFT1;
	L2 = SL*DRIFT2;
	K1 = SL*KICK1;
	K2 = SL*KICK2;
	
    if (max_order<2)
    {
        kv = 0;
    }
    else
        kv = B[1];
    
	if(T1==NULL)
	    useT1=false;
	else 
	    useT1=true;  
	    
    if(T2==NULL)
	    useT2=false; 
	else 
	    useT2=true;  
	
	if(R1==NULL)
	    useR1=false; 
	else 
	    useR1=true;  
	    
    if(R2==NULL)
	    useR2=false;
	else 
	    useR2=true;
    
    for(c = 0;c<num_particles;c++)	/*Loop over particles  */
			{	r6 = r+c*6;	
			    if(!mxIsNaN(r6[0]))
			    {   
					/*  misalignment at entrance  */
					if(useT1)
			            ATaddvv(r6,T1);
			        if(useR1)
			            ATmultmv(r6,R1);
                    
                    if(flag_quadedge!=NULL)
                    {
                        if (flag_quadedge[0]==3)
                            quadnonlinearfringe(r6,kv);
                       
                        if ((flag_quadedge[0]==1)|| (flag_quadedge[0]==3))
                            /*quadlinearfringe(r6, -kv*I1a);*/
                            quadlinearfringefull(r6,Iminus,Iplus,kv,0);
                    }
					/*  integrator  */
					for(m=0; m < num_int_steps; m++) /*  Loop over slices */
						{	r6 = r+c*6;	
							
							ladrift6(r6, L1);
           					strthinkick(r6, A, B,  K1, max_order);
           					ladrift6(r6, L2);
           					strthinkick(r6, A, B, K2, max_order);
           					ladrift6(r6, L2);
           					strthinkick(r6, A, B,  K1, max_order);
           					ladrift6(r6, L1);	
						}  
					if(flag_quadedge!=NULL)
                    {
                        if (flag_quadedge[1]==3)
                            quadnonlinearfringe(r6,-kv);
                       
                        if ((flag_quadedge[1]==1)||(flag_quadedge[1]==3))
                            quadlinearfringefull(r6,Iminus,Iplus,kv,1);
                            /*quadlinearfringe(r6, kv*I1b);*/
                    }
                    
                      
					/* Misalignment at exit */	
			        if(useR2)
			            ATmultmv(r6,R2);
		            if(useT2)   
			            ATaddvv(r6,T2);
			    }
			}
}



ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
								double *r_in, int num_particles, int mode)
#define NUM_FIELDS_2_REMEMBER 9
{	int fnum;
	double *A , *B;
	double  *pr1, *pr2, *pt1, *pt2;   
    double *Iminus, *Iplus;
    int flag_quadedge[2]={0,0};
    double *tmpflag;
	int max_order, num_int_steps;
	double le;
	int *returnptr;
	int *NewFieldNumbers;

	switch(mode)
		{	case NO_LOCAL_COPY:	/* NOT used in AT1.3 Get fields by names from MATLAB workspace   */
				{	
				
				}	break;	
	
			case MAKE_LOCAL_COPY: 	/* Find field numbers first
										Save a list of field number in an array
										and make returnptr point to that array
								    */
				{	
					/* Allocate memory for integer array of 
					   field numbers for faster futurereference
		            */
					NewFieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));

					/*  Populate */
					
					fnum = mxGetFieldNumber(ElemData,"PolynomA");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'PolynomA' was not found in the element data structure"); 
					NewFieldNumbers[0] = fnum;
					A = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
					
					fnum = mxGetFieldNumber(ElemData,"PolynomB");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'PolynomB' was not found in the element data structure"); 
					NewFieldNumbers[1] = fnum;
					B = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
					
					
					fnum = mxGetFieldNumber(ElemData,"MaxOrder");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'MaxOrder' was not found in the element data structure"); 
					NewFieldNumbers[2] = fnum;
					max_order = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					
					fnum = mxGetFieldNumber(ElemData,"NumIntSteps");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'NumIntSteps' was not found in the element data structure"); 
					NewFieldNumbers[3] = fnum;
					num_int_steps = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					
					fnum = mxGetFieldNumber(ElemData,"Length");
					if(fnum<0) 
					    mexErrMsgTxt("Required field 'Length' was not found in the element data structure"); 
					NewFieldNumbers[4] = fnum;
					le = mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
					
					fnum = mxGetFieldNumber(ElemData,"R1");
					NewFieldNumbers[5] = fnum;
					if(fnum<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					

					fnum = mxGetFieldNumber(ElemData,"R2");
					NewFieldNumbers[6] = fnum;
					if(fnum<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
					
                    fnum = mxGetFieldNumber(ElemData,"T1");
	                NewFieldNumbers[7] = fnum;
					if(fnum<0)
					    pt1 = NULL;
					else
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
	                
	                fnum = mxGetFieldNumber(ElemData,"T2");
	                NewFieldNumbers[8] = fnum;
					if(fnum<0)
					    pt2 = NULL;
					else
					    pt2 = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
								
					
                    fnum = mxGetFieldNumber(ElemData,"Iminus");
					NewFieldNumbers[9] = fnum;
					if(fnum<0) 
					    Iminus = 0;
					else
					   Iminus = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
					
					fnum = mxGetFieldNumber(ElemData,"Iplus");
					NewFieldNumbers[10] = fnum;
					if(fnum<0) 
					    Iplus = 0;
					else
					   Iplus =  mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
								
				
                    fnum = mxGetFieldNumber(ElemData,"QuadEdgeFlag");
					NewFieldNumbers[11] = fnum;
					if(fnum>=0) /*[a,b], a for entrance edge, b for exit edge. =0 for off, 1 for linear, 3 for linear and nonlinear*/
                    {   
                        tmpflag = mxGetPr(mxGetFieldByNumber(ElemData,0,fnum));
                        /*flag_quadedge[0] = (int) mxGetScalar(mxGetFieldByNumber(ElemData,0,fnum));
                        flag_quadedge[1] = (int) mxGetScalar(mxGetFieldByNumber(ElemData,1,fnum));
                         */
                        flag_quadedge[0] = (int) tmpflag[0];
                        flag_quadedge[1] = (int) tmpflag[1];
                    }
                            
					returnptr = NewFieldNumbers;

				}	break;

			case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
									    The second argument ponter to the array of field 
									    numbers is previously created with 
									    QuadLinPass( ..., MAKE_LOCAL_COPY)
											
									*/
				{	A = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
					B = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
					max_order = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
					num_int_steps = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[3]));
					le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
					
					/* Optional fields */
					if(FieldNumbers[5]<0)
					    pr1 = NULL;
					else
					    pr1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[5]));
					
					if(FieldNumbers[6]<0)
					    pr2 = NULL;
					else
					    pr2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[6]));
					
					    
					if(FieldNumbers[7]<0)
					    pt1 = NULL;
					else    
					    pt1 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[7]));
					    
					if(FieldNumbers[8]<0)
					    pt2 = NULL;
					else 
					    pt2 = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[8]));
                    
                    if(FieldNumbers[9]<0)
					    Iminus = 0;
					else 
					    Iminus = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[9]));
					    
					if(FieldNumbers[10]<0)
					    Iplus = 0;
					else 
					    Iplus = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[10])); 
                    
                    if(FieldNumbers[11]>=0)
                    { 
                      /*[a,b], a for entrance edge, b for exit edge. =0 for off, 1 for linear, 3 for linear and nonlinear*/
                        /*flag_quadedge[0] = (int) mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[11]));
                        flag_quadedge[1] = (int) mxGetScalar(mxGetFieldByNumber(ElemData,1,FieldNumbers[11]));
                        */
                        tmpflag =  mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[11])); 
                        flag_quadedge[0] = (int) tmpflag[0];
                        flag_quadedge[1] = (int) tmpflag[1];
                    }
                    
					returnptr = FieldNumbers;
				}	break;
			default:
				{	mexErrMsgTxt("No match for calling mode in function StrMPoleSymplectic4Pass\n");
				}
		}

    StrMPoleSymplectic4NPass(r_in, le, A, B, max_order, num_int_steps,Iminus,Iplus,flag_quadedge,
									pt1, pt2, pr1, pr2, num_particles);
	return(returnptr);
}

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	int m,n;
	double *r_in;
	double le, *A, *B, *pr1, *pr2, *pt1, *pt2;  
	int max_order, num_int_steps;
    double *Iminus, *Iplus;
    int flag_quadedge[2]={0,0};
    double *tmpflag;
	mxArray *tmpmxptr;

	if(nrhs)
	{
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
					
				    
	/* Optionnal arguments */    
	tmpmxptr = mxGetField(prhs[0],0,"R1");
	if(tmpmxptr)
	    pr1 = mxGetPr(tmpmxptr);
	else
	    pr1=NULL; 
	    
	tmpmxptr = mxGetField(prhs[0],0,"R2");
	if(tmpmxptr)
	    pr2 = mxGetPr(tmpmxptr);
	else
	    pr2=NULL; 

	tmpmxptr = mxGetField(prhs[0],0,"T1");

	
	if(tmpmxptr)
	    pt1=mxGetPr(tmpmxptr);
	else
	    pt1=NULL;

	tmpmxptr = mxGetField(prhs[0],0,"T2");
	if(tmpmxptr)
	    pt2=mxGetPr(tmpmxptr);
	else
	    pt2=NULL;
	
	
    /*Iminus and Iplus are fringe field, Iminus is before hard edge, Iplus is after hard edge */
	tmpmxptr = mxGetField(prhs[0],0,"Iminus"); 
	if(tmpmxptr)
	    Iminus=mxGetPr(tmpmxptr);
	else
	    Iminus=0;
	    
	tmpmxptr = mxGetField(prhs[0],0,"Iplus");  
	if(tmpmxptr)
	    Iplus=mxGetPr(tmpmxptr);
	else
	    Iplus=0;
    
    tmpmxptr = mxGetField(prhs[0],0,"QuadEdgeFlag");  
	if(tmpmxptr)
    {   /* flag_quadedge=mxGetPr(tmpmxptr);*/
        tmpflag = mxGetPr(tmpmxptr);
        flag_quadedge[0] = (int) tmpflag[0];
        flag_quadedge[1] = (int) tmpflag[1];
                        
        /*tmpmxptr = mxGetField(prhs[0],0,"QuadEdgeFlag");
        flag_quadedge[0]=(int)mxGetScalar(tmpmxptr);
        tmpmxptr = mxGetField(prhs[0],1,"QuadEdgeFlag");
        flag_quadedge[1]=(int)mxGetScalar(tmpmxptr);*/
    }   
	plhs[0] = mxDuplicateArray(prhs[1]);
	r_in = mxGetPr(plhs[0]);
	StrMPoleSymplectic4NPass(r_in, le, A, B, max_order, num_int_steps,
									Iminus,Iplus,flag_quadedge,pt1, pt2, pr1, pr2, n);
	}
	else
	{   /* return list of required fields */
	    plhs[0] = mxCreateCellMatrix(5,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("PolynomA"));
	    mxSetCell(plhs[0],2,mxCreateString("PolynomB"));
	    mxSetCell(plhs[0],3,mxCreateString("MaxOrder"));
        mxSetCell(plhs[0],4,mxCreateString("NumIntSteps"));	    
	    
        if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(4,1);
	        mxSetCell(plhs[1],0,mxCreateString("T1"));
	        mxSetCell(plhs[1],1,mxCreateString("T2"));
	        mxSetCell(plhs[1],2,mxCreateString("R1"));
	        mxSetCell(plhs[1],3,mxCreateString("R2"));
            mxSetCell(plhs[1],4,mxCreateString("Iminus"));
            mxSetCell(plhs[1],5,mxCreateString("Iplus"));
            mxSetCell(plhs[1],6,mxCreateString("QuadEdgeFlag"));
	    }
    }
}



