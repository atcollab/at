/*FDW.c for Accelerator Toolbox
 *----------------------------------------------------------------------------
 * Modification Log:
 * -----------------
 * .02  2024-05-09      J. Arenillas, ALBA, jarenillas@axt.email
 *                      Rewriting mex function.
 *						Cleaning up code, using constants in gwig.h.
 *                    	New factors for diffusion matrix computation.
 *						Updated version of radiation kicks in diffusion matrix propagation.
 *
 * .01  2018-06-01      A. Mash'al, ILSF, a-mashal@ipm.ir
 *                      Implementing the diffusion matrix for a wiggler
 *----------------------------------------------------------------------------
 *  mex-function to calculate integrated radiation diffusion matrix B defined in [2]
 *  for wiggler elements in MATLAB Accelerator Toolbox and based on [3].
 *
 * References
 * [1] M.Sands 'The Physics of Electron Storage Rings'
 * [2] Ohmi, Kirata, Oide 'From the beam-envelope matrix to synchrotron
 * radiation integrals', Phys.Rev.E  Vol.49 p.751 (1994)
 * [3] A. Mashal, F. D. Kashani, J. Rahighi, E. Ahmadi, Z. Marti, and O. J. Arranz.
 * Study the effect of insertion devices on electron beam properties by envelope method.
 * Journal of Instrumentation, May 2021.
 */

#include "atelem.c"
#include "atlalib.c"
#include "gwigR.c"
#include <math.h>
#include "gwig.h"

/* Fourth order-symplectic integrator constants */

#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656

/* Physical constants used in the calculations */

#define LAMBDABAR  3.8615926796e-13			/* Compton wavelength/2pi [m]    */
#define CU         1.323094366892892		/* 55/(24*sqrt(3)) factor        */

#define SQR(X) ((X)*(X))


static void wigglerM(struct gwigR *pWig, double* orbit_in, double L, double *M66)
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


static void wigglerB(struct gwigR *pWig, double* orbit_in, double L, double *B66)
{  /* Calculate Ohmi's diffusion matrix of a wiggler.
   Since wigglers have a straight coordinate system: irho=0
   The result is stored in a preallocated 1-dimentional array B66
   (of 36 elements) of matrix B arranged column-by-column
*/

  double ax,ay,kx,ky,axpy,aypx;
	double BB, E;
	double Brho,irho3,B2;
	double Bxyz[3];
	double p_norm = 1/(1+orbit_in[4]);
	double p_norm2 = SQR(p_norm);
	double px,py,D;
	double DLDS;
  int i;

	px= orbit_in[1];
	py= orbit_in[3];
	D = orbit_in[4];
	GWigAx(pWig, orbit_in, &ax, &axpy);
	GWigAy(pWig, orbit_in, &ay, &aypx);
	kx=px-ax;
	ky=py-ay;

	/* Calculate the local  magnetic field in T^2 */
    GWigB(pWig, orbit_in, Bxyz);
    B2 = (Bxyz[0]*Bxyz[0]) + (Bxyz[1]*Bxyz[1]) + (Bxyz[2]*Bxyz[2]);

    /* Beam rigidity in T*m */
	E = (pWig->E0)*(1+orbit_in[4]);
	Brho = E*1E9/clight;

    /* 1/rho^3 */
  irho3 = B2/(Brho*Brho)*sqrt(B2)/Brho;

	DLDS = (1+D)/(sqrt((1+D)*(1+D)-(px-ax)*(px-ax)-(py-ay)*(py-ay)));
	BB = CU*r_e*LAMBDABAR*pow(pWig->Po,5)*irho3*DLDS*L;

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


static void FindElemB(double *orbit_in, double le, double Lw, double Bmax,
               int Nstep, int Nmeth, int NHharm, int NVharm,
               double *By, double *Bx, double E0, double *T1,
               double *T2, double *R1, double *R2, double *BDIFF)
{	/* Find Ohmi's diffusion matrix BDIFF of a wiggler
	   BDIFF - cumulative Ohmi's diffusion is initialized to 0
	   BDIFF is preallocated 1 dimensional array to store 6-by-6 matrix
	   columnwise.
	*/

	int m;
	double  MKICK[36], BKICK[36];

	/* 4-th order symplectic integrator constants */
	int Niter = Nstep*(le/Lw);
  double SL = Lw/Nstep;
	double dl1 = SL*KICK1;
  double dl0 = SL*KICK2;

	/* Allocate memory for thin kick matrix MKICK
	   and a diffusion matrix BKICK
	*/
	for(m=0; m < 6; m++)
		{	MKICK[m] = 0;
			BKICK[m] = 0;
		}

	double zEndPointH[2];
  double zEndPointV[2];
	double ax,ay,axpy,aypx;
	double B[3];
	double E0G;
  struct gwigR pWig;
  int flag;

	/* Transform orbit to a local coordinate system of an element
       BDIFF stays zero	*/
    if(T1) ATaddvv(orbit_in,T1);
    if(R1) ATmultmv(orbit_in,R1);

    /*Generate the wiggler*/
    zEndPointH[0] = 0;
    zEndPointH[1] = le;
    zEndPointV[0] = 0;
    zEndPointV[1] = le;

	/* Energy is defined in the lattice in eV but GeV is used by the gwig code. */
	E0G = E0 / 1e9;
    GWigInit2(&pWig,E0G,le,Lw,Bmax,Nstep,Nmeth,NHharm,NVharm,0,0,zEndPointH,zEndPointV,By,Bx,T1,T2,R1,R2);

    /* Propagate orbit_in and BDIFF through a 4-th order integrator */
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
			wigglerB(&pWig, orbit_in, dl1, BKICK);
			ATsandwichmmt(MKICK,BDIFF);
		  ATaddmm(BKICK,BDIFF);
	   	GWigMap_2nd(&pWig, orbit_in, dl1);

			wigglerM(&pWig, orbit_in, dl0, MKICK);
			wigglerB(&pWig, orbit_in, dl0, BKICK);
			ATsandwichmmt(MKICK,BDIFF);
			ATaddmm(BKICK,BDIFF);
      GWigMap_2nd(&pWig, orbit_in, dl0);

			wigglerM(&pWig, orbit_in, dl1, MKICK);
			wigglerB(&pWig, orbit_in, dl1, BKICK);
			ATsandwichmmt(MKICK,BDIFF);
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

    if(R2) {
	   ATmultmv(orbit_in,R2);
	   ATsandwichmmt(R2,BDIFF);
	}
    if(T2) ATaddvv(orbit_in,T2);
}


static double *wiggdiffmatrix(const atElem *ElemData, double *orb, double energy, double *bdiff)
{
    double Ltot, Lw, Bmax;
	int Nstep, Nmeth;
	int NHharm, NVharm;
    double *R1, *R2, *T1, *T2;
    double *By, *Bx;

    /* Required fields */
	Ltot = atGetOptionalDouble(ElemData, "Length",0.0);
	/* If ELEMENT has a zero length, return zeros matrix end exit */
	if (Ltot == 0.0) return bdiff;

	Nstep = atGetLong(ElemData, "Nstep"); check_error();

	if (atIsNaN(energy)) {
        energy=atGetDouble(ElemData,"Energy"); check_error();
    }
    Lw = atGetDouble(ElemData, "Lw"); check_error();
    Bmax = atGetDouble(ElemData, "Bmax"); check_error();
    Nmeth = atGetLong(ElemData, "Nmeth"); check_error();
    NHharm = atGetLong(ElemData, "NHharm"); check_error();
    NVharm = atGetLong(ElemData, "NVharm"); check_error();
    By = atGetDoubleArray(ElemData, "By"); check_error();
    Bx = atGetDoubleArray(ElemData, "Bx"); check_error();

	/* Optional fields */
    R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
    R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
    T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
    T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();

	FindElemB(orb, Ltot, Lw, Bmax, Nstep, Nmeth, NHharm, NVharm, By, Bx, energy, T1, T2, R1, R2, bdiff);
    return bdiff;
}


#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	 /* The calling syntax of this mex-function from MATLAB is
     FindWiggRadDiffMatrix(ELEMENT, ORBIT)
     ELEMENT is the element structure with field names consistent with
             a multipole transverse field model.
     ORBIT is a 6-by-1 vector of the closed orbit at the entrance (calculated elsewhere)
*/

	mxDouble *orb0;
	const mxArray *mxElem, *mxOrbit;
	double *bdiff;
	double orb[6];
    double energy = mxGetNaN(); /* Default energy value */
    int i, m, n;

	if (nrhs < 2)
        mexErrMsgIdAndTxt("AT:WrongArg", """FDW"" needs at least 2 arguments");
    mxElem = prhs[0];
    mxOrbit = prhs[1];

	m = mxGetM(mxOrbit);
	n = mxGetN(mxOrbit);
	if (!(mxIsDouble(mxOrbit) && (m==6 && n==1)))
		mexErrMsgIdAndTxt("AT:WrongArg", "Orbit must be a double 6x1 vector");

	orb0 = mxGetDoubles(mxOrbit);
	/* make local copy of the input closed orbit vector */
	for (i=0;i<6;i++) orb[i] = orb0[i];

	if (nrhs >= 3) {
	    const mxArray *mxEnergy = prhs[2];
	    if (!(mxIsDouble(mxEnergy) && mxIsScalar(mxEnergy)))
	        mexErrMsgIdAndTxt("AT:WrongArg", "Energy must be a double scalar");
	    energy=mxGetScalar(mxEnergy);
	}

	/* ALLOCATE memory for the output array */
	plhs[0] = mxCreateDoubleMatrix(6,6,mxREAL);
	bdiff = mxGetDoubles(plhs[0]);
    for (i=0; i<36; i++) bdiff[i]=0.0;

    wiggdiffmatrix(mxElem, orb, energy, bdiff);
}
#endif


#if defined(PYAT)

#define MODULE_NAME wiggdiffmatrix
#define MODULE_DESCR "Computation of the radiation diffusion matrix for a wiggler"

static PyObject *compute_wiggdiffmatrix(PyObject *self, PyObject *args) {
    PyObject *pyElem, *pyMatrix;
    PyArrayObject *pyOrbit;
    double *orb0, *bdiff, *retval;
    double orb[6];
    double energy = NAN; /* Default energy value */
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

    retval = wiggdiffmatrix(pyElem, orb, energy, bdiff);
    if (retval == NULL) {
        Py_DECREF(pyMatrix);
        return NULL;
    }
    return pyMatrix;
}

static PyMethodDef AtMethods[] = {
    {"FDW",
    (PyCFunction)compute_wiggdiffmatrix, METH_VARARGS,
    PyDoc_STR(
    "FDW(element, orbit, energy)\n\n"
    "Computes the radiation diffusion matrix B defined in [2]_\n"
    "for wiggler elements\n\n"
    "Args:\n"
    "    element:    Lattice element\n"
    "    orbit:      (6,) closed orbit at the entrance of ``element``\n"
    "    energy:     particle energy\n\n"
    "Returns:\n"
    "    diffmatrix: The radiation diffusion matrix of the wiggler\n\n"
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
