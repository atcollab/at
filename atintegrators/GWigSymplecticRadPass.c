/* GWigSymplecticRadPass.c for
   Accelerator Toolbox
*/

/*
 *---------------------------------------------------------------------------
 * Modification Log:
 * -----------------
 * .04  2024-10-21     L. Farvacque
 *              Merged tracking and diffusion matrices
 * .03  2024-05-06     J. Arenillas, ALBA, jarenillas@axt.email
 *              Adding rotations and translations to wiggler.
 *				Bug fix in wiggler initialisation.
 *				Energy parameter bug fix.
 * .02  2003-06-18     J. Li
 *				Cleanup the code
 *
 * .01  2003-04-20     YK Wu
 *				GWiggler interface
 *
 *---------------------------------------------------------------------------
 *  Accelerator Physics Group, Duke FEL Lab, www.fel.duke.edu
 */

#include "atelem.c"
#include "atlalib.c"
#include "wigrad.c"

struct elem {
    double Energy;
    double Length;
    double Lw;
    double Bmax;
    int Nstep;
    int Nmeth;
    int NHharm;
    int NVharm;
    double *By;
    double *Bx;
    /* Optional fields */
    double *R1;
    double *R2;
    double *T1;
    double *T2;
};

/*****************************************************************************/
/* PHYSICS SECTION ***********************************************************/


static void wigglerM(struct gwig *pWig, double* orbit_in, double L, double *bdiff)
{
    /* Computes the transfer matrix for a wiggler. */
	double M66[36];
	double H2[36];
	Hessian(pWig, orbit_in ,H2);

	M66[0]   = 1.0 + H2[6]; /* H[1, 0] */
	M66[1]   =-H2[0];       /* H[0, 0] */
	M66[2]   = H2[18];      /* H[3, 0] */
	M66[3]   =-H2[12];      /* H[2, 0] */
	M66[4]   = H2[30];      /* H[5, 0] */
	M66[5]   =-H2[24];      /* H[4, 0] */

	M66[6]   = L*H2[6]+H2[7];       /* H[1, 0] + H[1, 1] */
	M66[7]   = 1.0 - L*H2[0]-H2[1]; /* H[0, 0] - H[0, 1] */
	M66[8]   = L*H2[18]+H2[19];     /* H[3, 0] + H[3, 1] */
	M66[9]   =-L*H2[12]-H2[13];     /* H[2, 0] - H[2, 1] */
	M66[10]  = L*H2[30]+H2[31];     /* H[5, 0] + H[5, 1] */
	M66[11]  =-L*H2[24]-H2[25];     /* H[4, 0] - H[4, 1] */

	M66[12]  = H2[8];        /* H[1, 2] */
	M66[13]  =-H2[2];        /* H[0, 2] */
	M66[14]  = 1.0 + H2[20]; /* H[3, 2] */
	M66[15]  =-H2[14];       /* H[2, 2] */
	M66[16]  = H2[32];       /* H[5, 2] */
	M66[17]  =-H2[26];       /* H[4, 2] */

	M66[18]  = L*H2[8]+H2[9];         /* H[1, 2] + H[1, 3] */
	M66[19]  =-L*H2[2]-H2[3];         /* H[0, 2] - H[0, 3] */
	M66[20]  = L*H2[20]+H2[21];       /* H[3, 2] + H[3, 3] */
	M66[21]  = 1.0 - L*H2[14]-H2[15]; /* H[2, 2] - H[2, 3] */
	M66[22]  = L*H2[32]+H2[33];       /* H[5, 2] - H[5, 3] */
	M66[23]  =-L*H2[26]-H2[27];       /* H[4, 2] - H[4, 3] */

	M66[24]  = H2[10];       /* H[1, 4] */
	M66[25]  =-H2[4];        /* H[0, 4] */
	M66[26]  = H2[22];       /* H[3, 4] */
	M66[27]  =-H2[16];       /* H[2, 4] */
	M66[28]  = 1.0 + H2[34]; /* H[5, 4] */
	M66[29]  =-H2[28];       /* H[4, 4] */

	M66[30]  = H2[11];       /* H[1, 5] */
	M66[31]  =-H2[5];        /* H[0, 5] */
	M66[32]  = H2[23];       /* H[3, 5] */
	M66[33]  =-H2[17];       /* H[2, 5] */
	M66[34]  = H2[35];       /* H[5, 5] */
	M66[35]  = 1.0 - H2[29]; /* H[4, 5] */

    ATsandwichmmt(M66, bdiff);
}


static void wigglerB(struct gwig *pWig, double* orbit_in, double L, double *bdiff)
{  /* Calculate Ohmi's diffusion matrix of a wiggler.
   Since wigglers have a straight coordinate system: irho=0
   The result is stored in a preallocated 1-dimentional array B66
   (of 36 elements) of matrix B arranged column-by-column
*/

    double B66[36];
    double ax,ay,kx,ky,axpy,aypx;
	double BB, E;
	double Brho,irho3,B2;
	double Bxyz[3];
    double gamma0 = sqrt(1.0 + pWig->Po*pWig->Po);
	double px= orbit_in[1];
	double py= orbit_in[3];
	double D = orbit_in[4];
	double p_norm = 1.0/(1.0+D);
	double DLDS;

	GWigAx(pWig, orbit_in, &ax, &axpy);
	GWigAy(pWig, orbit_in, &ay, &aypx);
	kx=(px-ax)*p_norm;
	ky=(py-ay)*p_norm;

	/* Calculate the local  magnetic field in T^2 */
    GWigB(pWig, orbit_in, Bxyz);
    B2 = (Bxyz[0]*Bxyz[0]) + (Bxyz[1]*Bxyz[1]) + (Bxyz[2]*Bxyz[2]);

    /* Beam rigidity in T*m */
	Brho = 1.0e9 * pWig->Po * (1.0+D) *__E0 / C0;

    /* 1/rho^3 */
    irho3 = B2/(Brho*Brho)*sqrt(B2)/Brho;

	DLDS = 1.0/sqrt(1.0 - kx*kx - ky*ky);
	BB = DIF_CONST * pow(gamma0, 5) * irho3 * DLDS * L;

	/* When a 6-by-6 matrix is represented in MATLAB as one-dimentional
	   array containing 36 elements arranged column-by-column,
	   the relationship between indexes  is
	   [i][j] <---> [i+6*j]
	*/

	/* initialize B66 to 0 */
	for (int i=0; i<36; i++)
		B66[i] = 0.0;

	/* Populate B66 */
	B66[7]      = BB*kx*kx;  /* B66[1, 1] */
	B66[19]     = BB*kx*ky;  /* B66[1, 3] */
	B66[9]      = B66[19];   /* B66[3, 1] */
	B66[21]     = BB*ky*ky;  /* B66[3, 3] */
	B66[10]     = BB*kx;     /* B66[4, 1] */
	B66[25]     = B66[10];   /* B66[1, 4] */
	B66[22]     = BB*ky;     /* B66[4, 3] */
	B66[27]     = B66[22];   /* B66[3, 4] */
	B66[28]     = BB;        /* B66[4, 4] */

    ATaddmm(B66, bdiff);
}

void GWigSymplecticRadPass(double *orbit_in, double gamma, double le, double Lw,
            double Bmax, int Nstep, int Nmeth, int NHharm, int NVharm,
            double *By, double *Bx, double *T1, double *T2,
            double *R1, double *R2, int num_particles, double *bdiff)
{
	double ax,ay,axpy,aypx;
    struct gwig pWig;
	double B[3];
    double SL = Lw/Nstep;
	double dl1 = SL*KICK1;
    double dl0 = SL*KICK2;
	int Niter = Nstep*(le/Lw);

    GWigInit2(&pWig, gamma,le, Lw, Bmax, Nstep, Nmeth, NHharm, NVharm,0,0,By,Bx,T1,T2,R1,R2);

    for (int c = 0; c<num_particles; c++) { /* Loop over particles */
        double *r6 = orbit_in + c*6;
        pWig.Zw = 0.0;
        if(!atIsNaN(r6[0])) {
			/* Misalignment at entrance */
			if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);
            GWigGauge(&pWig, r6, 1);

            GWigAx(&pWig, r6, &ax, &axpy);
            GWigAy(&pWig, r6, &ay, &aypx);
            GWigB(&pWig, r6, B);
            r6[1] -= ax;
            r6[3] -= ay;
            GWigRadiationKicks(&pWig, r6, B, SL);
            r6[1] += ax;
            r6[3] += ay;

            /* integrator */
            for (int m=0; m < Niter; m++) { /* Loop over slices */
                if (bdiff) {
                    wigglerM(&pWig, r6, dl1, bdiff);
                    wigglerB(&pWig, r6, dl1, bdiff);
                }
                GWigMap_2nd(&pWig, r6, dl1);

                if (bdiff) {
                    wigglerM(&pWig, r6, dl0, bdiff);
                    wigglerB(&pWig, r6, dl0, bdiff);
                }
                GWigMap_2nd(&pWig, r6, dl0);

                if (bdiff) {
                    wigglerM(&pWig, r6, dl1, bdiff);
                    wigglerB(&pWig, r6, dl1, bdiff);
                }
                GWigMap_2nd(&pWig, r6, dl1);

                GWigAx(&pWig, r6, &ax, &axpy);
                GWigAy(&pWig, r6, &ay, &aypx);
                GWigB(&pWig, r6, B);
                r6[1] -= ax;
                r6[3] -= ay;
                GWigRadiationKicks(&pWig, r6, B, SL);
                r6[1] += ax;
                r6[3] += ay;
		    }

            GWigGauge(&pWig, r6, -1);

            if(R2) {
               ATmultmv(r6,R2);
               ATsandwichmmt(R2,bdiff);
            }
            if(T2) ATaddvv(r6,T2);
        }
    }
}

/********** END PHYSICS SECTION **********************************************/
/*****************************************************************************/

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)

{
    double gamma;
    double *bdiff = Param->bdiff;

    if (!Elem) {
        double *R1, *R2, *T1, *T2;
        double *By, *Bx;
        double Ltot, Lw, Bmax, Energy;
        int Nstep, Nmeth;
        int NHharm, NVharm;

        Ltot = atGetDouble(ElemData, "Length"); check_error();
        Lw = atGetDouble(ElemData, "Lw"); check_error();
        Bmax = atGetDouble(ElemData, "Bmax"); check_error();
        Nstep = atGetLong(ElemData, "Nstep"); check_error();
        Nmeth = atGetLong(ElemData, "Nmeth"); check_error();
        NHharm = atGetLong(ElemData, "NHharm"); check_error();
        NVharm = atGetLong(ElemData, "NVharm"); check_error();
        By = atGetDoubleArray(ElemData, "By"); check_error();
        Bx = atGetDoubleArray(ElemData, "Bx"); check_error();
        /* Optional fields */
        Energy = atGetOptionalDouble(ElemData, "Energy",Param->energy); check_error();
        R1 = atGetOptionalDoubleArray(ElemData, "R1"); check_error();
        R2 = atGetOptionalDoubleArray(ElemData, "R2"); check_error();
        T1 = atGetOptionalDoubleArray(ElemData, "T1"); check_error();
        T2 = atGetOptionalDoubleArray(ElemData, "T2"); check_error();

        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Energy=Energy;
        Elem->Length=Ltot;
        Elem->Lw=Lw;
        Elem->Bmax=Bmax;
        Elem->Nstep=Nstep;
        Elem->Nmeth=Nmeth;
        Elem->NHharm=NHharm;
        Elem->NVharm=NVharm;
        Elem->By=By;
        Elem->Bx=Bx;
        /* Optional fields */
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
    }
    gamma = atGamma(Param->energy, Elem->Energy, Param->rest_energy);

    GWigSymplecticRadPass(r_in, gamma, Elem->Length, Elem->Lw, Elem->Bmax,
            Elem->Nstep, Elem->Nmeth, Elem->NHharm, Elem->NVharm,
            Elem->By, Elem->Bx, Elem->T1, Elem->T2, Elem->R1, Elem->R2,
            num_particles, bdiff);
  return Elem;
}

/********** END WINDOWS DLL GATEWAY SECTION **********************************/

/********** MATLAB GATEWAY ***************************************************/

MODULE_DEF(GWigSymplecticRadPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(       int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double rest_energy = 0.0;
        double charge = -1.0;
        double *r_in;
        double Gamma;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double *By, *Bx;
        double *R1, *R2, *T1, *T2;
        double Ltot, Lw, Bmax, Energy;
        int Nstep, Nmeth;
        int NHharm, NVharm;
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");

        Ltot = atGetDouble(ElemData, "Length"); check_error();
        Lw = atGetDouble(ElemData, "Lw"); check_error();
        Bmax = atGetDouble(ElemData, "Bmax"); check_error();
        Nstep = atGetLong(ElemData, "Nstep"); check_error();
        Nmeth = atGetLong(ElemData, "Nmeth"); check_error();
        NHharm = atGetLong(ElemData, "NHharm"); check_error();
        NVharm = atGetLong(ElemData, "NVharm"); check_error();
        By = atGetDoubleArray(ElemData, "By"); check_error();
        Bx = atGetDoubleArray(ElemData, "Bx"); check_error();
        /* Optional fields */
        Energy = atGetOptionalDouble(ElemData, "Energy",0.0); check_error();
        R1 = atGetOptionalDoubleArray(ElemData, "R1"); check_error();
        R2 = atGetOptionalDoubleArray(ElemData, "R2"); check_error();
        T1 = atGetOptionalDoubleArray(ElemData, "T1"); check_error();
        T2 = atGetOptionalDoubleArray(ElemData, "T2"); check_error();
        if (nrhs > 2) atProperties(prhs[2], &Energy, &rest_energy, &charge);

        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        Gamma = atGamma(Energy, Energy, rest_energy);
        r_in = mxGetDoubles(plhs[0]);

        GWigSymplecticRadPass(r_in, Gamma, Ltot, Lw, Bmax, Nstep, Nmeth,
            NHharm, NVharm, By, Bx, T1, T2, R1, R2, num_particles, NULL);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(10,1);
        mxSetCell(plhs[0],0,mxCreateString("Energy"));
        mxSetCell(plhs[0],1,mxCreateString("Length"));
        mxSetCell(plhs[0],2,mxCreateString("Lw"));
        mxSetCell(plhs[0],3,mxCreateString("Bmax"));
        mxSetCell(plhs[0],4,mxCreateString("Nstep"));
        mxSetCell(plhs[0],5,mxCreateString("Nmeth"));
        mxSetCell(plhs[0],6,mxCreateString("NHharm"));
        mxSetCell(plhs[0],7,mxCreateString("NVharm"));
        mxSetCell(plhs[0],8,mxCreateString("By"));
        mxSetCell(plhs[0],9,mxCreateString("Bx"));
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(4,1);
            mxSetCell(plhs[1],0,mxCreateString("R1"));
            mxSetCell(plhs[1],1,mxCreateString("R2"));
            mxSetCell(plhs[1],2,mxCreateString("T1"));
            mxSetCell(plhs[1],3,mxCreateString("T2"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}

#endif /*MATLAB_MEX_FILE*/
