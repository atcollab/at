/* FDW.c

   mex-function to calculate integrated radiation diffusion matrix B defined in [2] 
   for wiggler elements in MATLAB Accelerator Toolbox
   O.Jimenez 3/08/18

   References
   [1] M.Sands 'The Physics of Electron Storage Rings
   [2] Ohmi, Kirata, Oide 'From the beam-envelope matrix to synchrotron
   radiation integrals', Phys.Rev.E  Vol.49 p.751 (1994)
*/

#include "mex.h"
#include "matrix.h"
#include "atlalib.c"
#include <math.h>
#include "gwigR.c"


/* Fourth order-symplectic integrator constants */

#define x1     1.351207191959657328
#define x0    -1.702414383919314656

/* Physical constants used in the calculations */

#define TWOPI	   6.28318530717959
#define CGAMMA 	   8.846056192e-05 			/* [m]/[GeV^3] Ref[1] (4.1)      */
#define M0C2       5.10999060e5				/* Electron rest mass [eV]       */
#define LAMBDABAR  3.86159323e-13			/* Compton wavelength/2pi [m]    */
#define CER   	   2.81794092e-15			/* Classical electron radius [m] */
#define CU         1.323094366892892		/* 55/(24*sqrt(3)) factor        */

#define SQR(X) ((X)*(X))




void wigglerM(struct gwigR *pWig, double* orbit_in, double L, double *M66)
{
    /* Computes the transfer matrix for a wiggler. */
	int i,j,k;
	double H2[36];
	double H[6][6];
	Hessian(pWig, orbit_in ,H2);
	
	for (i=0;i<6;i++){
		for(j=0;j<6;j++){
			k=j+i*6;
			H[i][j]=H2[k];
		}
	}
	/* Initialize M66 to a 1-by-36 zero vector. */
	for(i=0;i<36;i++){
		M66[i]= 0;
	}

	M66[0]   =1+H[1][0];
	M66[1]   =-H[0][0];
	M66[2]   = H[3][0];
	M66[3]   =-H[2][0];
	M66[4]   = H[5][0];
	M66[5]   =-H[4][0];
	
	M66[6]   = L*H[1][0]+H[1][1];
	M66[7]   =1-L*H[0][0]-H[0][1];
	M66[8]   = L*H[3][0]+H[3][1];
	M66[9]   =-L*H[2][0]-H[2][1];
	M66[10]  = L*H[5][0]+H[5][1];
	M66[11]  =-L*H[4][0]-H[4][1];
	
	M66[12]  = H[1][2];
	M66[13]  =-H[0][2];
	M66[14]  =1+H[3][2];
	M66[15]  =-H[2][2];
	M66[16]  = H[5][2];
	M66[17]  =-H[4][2];
	
	M66[18]  = L*H[1][2]+H[1][3];
	M66[19]  =-L*H[0][2]-H[0][3];
	M66[20]  = L*H[3][2]+H[3][3];
	M66[21]  =1-L*H[2][2]-H[2][3];
	M66[22]  = L*H[5][2]+H[5][3];
	M66[23]  =-L*H[4][2]-H[4][3];
	
	M66[24]  = H[1][4];
	M66[25]  =-H[0][4];
	M66[26]  = H[3][4];
	M66[27]  =-H[2][4];
	M66[28]  =1+H[5][4];
	M66[29]  =-H[4][4];
	
	M66[30]  = H[1][5];
	M66[31]  =-H[0][5];
	M66[32]  = H[3][5];
	M66[33]  =-H[2][5];
	M66[34]  = H[5][5];
	M66[35]  =1-H[4][5];

}


void wigglerB(struct gwigR *pWig, double* orbit_in, double L, 
				 double E0, double *B66)

/* Calculate Ohmi's diffusion matrix of a wiggler. 
   Since wigglers have a straight coordinate system: irho=0
   The result is stored in a preallocated 1-dimentional array B66
   (of 36 elements) of matrix B arranged column-by-column
*/

{	
	double ax,ay,kx,ky,axpy,aypx;
	double BB,B3P;
	double H,irho2,B2;
	double Bxy[2];
	double p_norm = 1/(1+orbit_in[4]);
	double p_norm2 = SQR(p_norm);
    double c = 299792458;   /* m s-1 */
    double e = 1;
    double irho = 0;
	double x,px,y,py,D;
	double DLDS;
    int i;

	x =orbit_in[0];
	px=orbit_in[1];
	y =orbit_in[2];
	py=orbit_in[3];
	D =orbit_in[4];
	GWigAx(pWig, orbit_in, &ax, &axpy);
	GWigAy(pWig, orbit_in, &ay, &aypx);
	kx=px-ax;
	ky=py-ay;
          
	DLDS=L*SQR(SQR(1+orbit_in[4]))*(1+((px-ax)*(px-ax)/(2*(1+D)*(1+D)))+((py-ay)*(py-ay)/(2*(1+D)*(1+D))));
	
	//printf("DLDS=%e DLDS2=%e\n",(px)*(px)/(2*(1+D)*(1+D)),(px-ax)*(px-ax)/(2*(1+D)*(1+D)));

  	/* calculate the local  magnetic field in m-1 */  
    GWigB(pWig, orbit_in, Bxy);
    B2 = (Bxy[0]*Bxy[0]) + (Bxy[1]*Bxy[1]);
  
    /* Beam rigidity in T*m */
    H = (pWig->Po)/586.679074042074490;

    /* 1/rho^2 */
    irho2 = B2/(H*H);
	B3P = irho2*sqrt(irho2);

    BB=CU * CER * LAMBDABAR * pow(E0/M0C2,5) * B3P * DLDS;
    /* calculate |B x n|^3 - the third power of the B field component 
	   orthogonal to the normalized velocity vector n
    */
	
   
     
	//printf("BB=%e, BB2=%e \n",(1+((px)*(px)/(2*(1+D)*(1+D)))+((py-ay)*(py-ay)/(2*(1+D)*(1+D)))),(1+((px-ax)*(px-ax)/(2*(1+D)*(1+D)))+((py-ay)*(py-ay)/(2*(1+D)*(1+D)))));

	/* When a 6-by-6 matrix is represented in MATLAB as one-dimentional 
	   array containing 36 elements arranged column-by-column,
	   the relationship between indexes  is
	   [i][j] <---> [i+6*j] 

	*/
	
	/* initialize B66 to 0 */
	for(i=0;i<36;i++)
		B66[i] = 0;
	
	/* Populate B66 */
	B66[7]      = BB*SQR(kx)*p_norm2;
	B66[19]     = BB*kx*ky*p_norm2;
	B66[9]      = B66[19];
	B66[21]     = BB*SQR(ky)*p_norm2;
	B66[10]     = BB*kx*p_norm;
	B66[25]     = B66[10];
	B66[22]     = BB*ky*p_norm;
	B66[27]     = B66[22];
	B66[28]     = BB;
	
	
}


void drift_propagateB(double *orb_in, double L,  double *B)
{	/* Propagate cumulative Ohmi's diffusion matrix B through a drift.
	   B is a (*double) pointer to 1-dimentional array 
	   containing 36 elements of matrix elements arranged column-by-column
	   as in MATLAB representation 

	   The relationship between indexes when a 6-by-6 matrix is 
	   represented in MATLAB as one-dimentional array containing
	   36 elements arranged column-by-column is
	   [i][j] <---> [i+6*j] 
	*/
		
	int m;
		
	double *DRIFTMAT = (double*)mxCalloc(36,sizeof(double));
	for(m=0;m<36;m++)
		DRIFTMAT[m] = 0;
	/* Set diagonal elements to 1	*/
	for(m=0;m<6;m++)
		DRIFTMAT[m*7] = 1;

	DRIFTMAT[6]  =  L/(1+orb_in[4]);
	DRIFTMAT[20] =  DRIFTMAT[6];
	DRIFTMAT[24] = -L*orb_in[1]/SQR(1+orb_in[4]);
	DRIFTMAT[26] = -L*orb_in[3]/SQR(1+orb_in[4]);
	DRIFTMAT[11] =  L*orb_in[1]/SQR(1+orb_in[4]);
	DRIFTMAT[23] =  L*orb_in[3]/SQR(1+orb_in[4]);	
	DRIFTMAT[29] = -L*(SQR(orb_in[1])+SQR(orb_in[3]))/((1+orb_in[4])*SQR(1+orb_in[4]));

	ATsandwichmmt(DRIFTMAT,B);
	mxFree(DRIFTMAT);
	for (m=0;m<36;m++)
		printf("B[%d]=%e \n",m,B[m]);
	
}


void FindElemB(double *orbit_in, double le, double Lw, double Bmax,
               int Nstep, int Nmeth, int NHharm, int NVharm,
               double *pBy, double *pBx, double E0, double *pt1,
               double *pt2, double *PR1, double *PR2, double *BDIFF)

{	/* Find Ohmi's diffusion matrix BDIFF (la integrada) of a wiggler
	   BDIFF - cumulative Ohmi's diffusion is initialized to 0
	   BDIFF is preallocated 1 dimensional array to store 6-by-6 matrix 
	   columnwise.
	*/
	
	int m;	
	double  *MKICK, *BKICK;

	/* 4-th order symplectic integrator constants */
	double dl1,dl0;
    int Niter = Nstep*(le/Lw);
    double SL = Lw/Nstep;
	double zEndPointH[2];
    double zEndPointV[2];
	double ax,ay,axpy,aypx;
	double B[2];
	double E0G;
    struct gwigR pWig;
    int flag;       
    
	dl1 = x1*SL;
    dl0 = x0*SL;
	
	/* Allocate memory for thin kick matrix MKICK
	   and a diffusion matrix BKICK
	*/
 	MKICK = (double*)mxCalloc(36,sizeof(double));
	BKICK = (double*)mxCalloc(36,sizeof(double));
	for(m=0; m < 6; m++)
		{	MKICK[m] = 0;
			BKICK[m] = 0;
		}
	
	/* Transform orbit to a local coordinate system of an element
       BDIFF stays zero	*/
    if(pt1)
        ATaddvv(orbit_in,pt1);	
    if(PR1)
        ATmultmv(orbit_in,PR1);	

   
    /*Generate the wiggler*/ 

  
  zEndPointH[0] = 0;
  zEndPointH[1] = le;
  zEndPointV[0] = 0;
  zEndPointV[1] = le;
  
  E0G = E0 / 1e9;

		  
  GWigInit(&pWig, E0G, le, Lw, Bmax, Nstep, Nmeth,NHharm,NVharm,0,0,zEndPointH,zEndPointV,pBy,pBx,pt1,pt2,PR1,PR2);
    
 /*   GWigInit(&pWig, le, Lw, Bmax, Nstep, Nmeth,NHharm,NVharm,pBy,pBx,
              E0,pt1,pt2,PR1,PR2);
    
 */   
	/* Propagate orbit_in and BDIFF through a 4-th orderintegrator */
        
    flag=1;
    GWigGauge(&pWig, orbit_in, flag);
 
    GWigAx(&pWig, orbit_in, &ax, &axpy);
    GWigAy(&pWig, orbit_in, &ay, &aypx);
    GWigB(&pWig, orbit_in, B);

    orbit_in[1] -= ax;
    orbit_in[3] -= ay;
    GWigRadiationKicks(&pWig, orbit_in, B, SL);
    orbit_in[1] += ax;
    orbit_in[3] += ay;
	
    for (m = 1; m <= Niter; m++ )   /* Loop over slices	*/
		{	   

			   wigglerM(&pWig, orbit_in, dl1, MKICK);
			   wigglerB(&pWig, orbit_in, dl1, E0, BKICK);
			   ATsandwichmmt(MKICK,BKICK);
			   ATaddmm(BKICK,BDIFF);
			   GWigMap_2nd(&pWig, orbit_in, dl1);
                
			   wigglerM(&pWig, orbit_in, dl0, MKICK);
			   wigglerB(&pWig, orbit_in, dl0, E0, BKICK);
			   ATsandwichmmt(MKICK,BKICK);
			   ATaddmm(BKICK,BDIFF);
               GWigMap_2nd(&pWig, orbit_in, dl0);
               
			   wigglerM(&pWig, orbit_in, dl1, MKICK);
			   wigglerB(&pWig, orbit_in, dl1, E0, BKICK);
			   ATsandwichmmt(MKICK,BKICK);
			   ATaddmm(BKICK,BDIFF);
               GWigMap_2nd(&pWig, orbit_in, dl1);
			   
			   GWigAx(&pWig, orbit_in, &ax, &axpy);
			   GWigAy(&pWig, orbit_in, &ay, &aypx);
               GWigB(&pWig, orbit_in, B);
			   orbit_in[1] -= ax;
               orbit_in[3] -= ay;
               GWigRadiationKicks(&pWig, orbit_in, B, SL);
               orbit_in[1] += ax;
               orbit_in[3] += ay;	
		}      
        
        flag=-1;
        GWigGauge(&pWig, orbit_in, flag);
       
		
        if(PR2)
            ATmultmv(orbit_in,PR2);	
        if(pt2)
            ATaddvv(orbit_in,pt2);	
        
        
		mxFree(MKICK);
		mxFree(BKICK);
}


#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* The calling syntax of this mex-function from MATLAB is
   FindWiggRadDiffMatrix(ELEMENT, ORBIT)
   ELEMENT is the element structure with field names consistent with 
           a multipole transverse field model.
   ORBIT is a 6-by-1 vector of the closed orbit at the entrance (calculated elsewhere)
*/
{	int m,n;  
	double E0;
    double Ltot, Lw, Bmax; 
    double *pBy, *pBx;
	double *BDIFF;
    int Nstep, Nmeth;
    int NHharm, NVharm;
    mxArray *tmpmxptr;

	double *orb, *orb0;
	double *pt1, *pt2, *PR1, *PR2;
    


	m = mxGetM(prhs[1]);
	n = mxGetN(prhs[1]);
	if(!(m==6 && n==1))
		mexErrMsgTxt("Second argument must be a 6-by-1 column vector");
    
	/* ALLOCATE memory for the output array */
	plhs[0] = mxCreateDoubleMatrix(6,6,mxREAL);
	BDIFF = mxGetPr(plhs[0]);

    

	orb0 = mxGetPr(prhs[1]);
	/* make local copy of the input closed orbit vector */
	orb = (double*)mxCalloc(6,sizeof(double));
	for(m=0;m<6;m++)
		orb[m] = orb0[m];
    
	/* Retrieve element information */
    
     
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

  tmpmxptr = mxGetField(prhs[0],0,"Energy");
  if(tmpmxptr)
    E0 = mxGetScalar(tmpmxptr);
  else
    mexErrMsgTxt("Required field 'Energy' was not found in the element data structure"); 
  
  tmpmxptr = mxGetField(prhs[0],0,"R1");
  if(tmpmxptr)
    PR1 = mxGetPr(tmpmxptr);
  else
    PR1 = NULL; 
  
  tmpmxptr = mxGetField(prhs[0],0,"R2");
  if(tmpmxptr)
    PR2 = mxGetPr(tmpmxptr);
  else
    PR2 = NULL; 
  
  tmpmxptr = mxGetField(prhs[0],0,"T1");
  if(tmpmxptr)
    pt1 = mxGetPr(tmpmxptr);
  else
    pt1 = NULL; 
  
  tmpmxptr = mxGetField(prhs[0],0,"T2");
  if(tmpmxptr)
    pt2 = mxGetPr(tmpmxptr);
  else
    pt2 = NULL;
    

	FindElemB(orb, Ltot, Lw, Bmax, Nstep, Nmeth, NHharm, NVharm,
                pBy, pBx, E0, pt1, pt2, PR1, PR2, BDIFF);
}
#endif /*MATLAB_MEX_FILE*/

#if defined(PYAT)

#define MODULE_NAME diffmatrix
#define MODULE_DESCR "Computation of the radiation diffusion matrix"

static PyObject *compute_diffmatrix(PyObject *self, PyObject *args) {
    PyObject *pyElem, *pyMatrix;
    PyArrayObject *pyOrbit;
    double *orb0, *bdiff, *retval;
    double orb[6];
    double energy;
    npy_intp outdims[2] = {6, 6};
    int i;

    if (!PyArg_ParseTuple(args, "OO!d", &pyElem, &PyArray_Type, &pyOrbit, &energy)) {
        return NULL;
    }
    if (PyArray_DIM(pyOrbit,0) != 6) {
        PyErr_SetString(PyExc_ValueError, "Orbit is not a (6,) array");
        return NULL;
    }
    if (PyArray_TYPE(pyOrbit) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "Orbit is not a double array");
        return NULL;
    }
    if ((PyArray_FLAGS(pyOrbit) & NPY_ARRAY_FARRAY_RO) != NPY_ARRAY_FARRAY_RO) {
        PyErr_SetString(PyExc_ValueError, "Orbit is not Fortran-aligned");
        return NULL;
    }

    orb0 = PyArray_DATA(pyOrbit);
	/* make local copy of the input closed orbit vector */
	for (i=0;i<6;i++) orb[i] = orb0[i];

	/* ALLOCATE memory for the output array */
    pyMatrix = PyArray_ZEROS(2, outdims, NPY_DOUBLE, 1);
    bdiff = PyArray_DATA((PyArrayObject *)pyMatrix);

    retval = diffmatrix(pyElem, orb, energy, bdiff);
    if (retval == NULL) {
        Py_DECREF(pyMatrix);
        return NULL;
    }
    return pyMatrix;
}

static PyMethodDef AtMethods[] = {
    {"find_mpole_raddiff_matrix",
    (PyCFunction)compute_diffmatrix, METH_VARARGS,
    PyDoc_STR(
    "find_mpole_raddiff_matrix(element, orbit, energy)\n\n"
    "Computes the radiation diffusion matrix B defined in [2]_\n"
    "for multipole elements\n\n"
    "Args:\n"
    "    element:    Lattice element\n"
    "    orbit:      (6,) closed orbit at the entrance of ``element``\n"
    "    energy:     particle energy\n\n"
    "Returns:\n"
    "    diffmatrix: The radiation diffusion matrix\n\n"
    "References:\n"
    "    **[1]** M.Sands, *The Physics of Electron Storage Rings*\n\n"
    "    .. [2] Ohmi, Kirata, Oide, *From the beam-envelope matrix to synchrotron\n"
    "       radiation integrals*, Phys.Rev.E  Vol.49 p.751 (1994)\n"
	)},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC MOD_INIT(MODULE_NAME)
{

#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    STR(MODULE_NAME), /* m_name */
    PyDoc_STR(MODULE_DESCR),      /* m_doc */
    -1,           /* m_size */
    AtMethods,    /* m_methods */
    NULL,         /* m_reload */
    NULL,         /* m_traverse */
    NULL,         /* m_clear */
    NULL,         /* m_free */
    };
    PyObject *m = PyModule_Create(&moduledef);
#else
    PyObject *m = Py_InitModule3(STR(MODULE_NAME), AtMethods,
        MODULE_DESCR);
#endif
    if (m == NULL) return MOD_ERROR_VAL;
    import_array();
    return MOD_SUCCESS_VAL(m);
}

#endif /*PYAT*/