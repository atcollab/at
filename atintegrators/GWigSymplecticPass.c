/* GWigSymplecticPass.c for
   Accelerator Toolbox 
*/

/*
 *---------------------------------------------------------------------------
 * Modification Log:
 * -----------------
 * .03  2024-05-06     J. Arenillas, ALBA, jarenillas@axt.email
 *              Adding rotations and translations to wiggler.
 *				Bug fix in wiggler initialisation.
 *				Energy parameter bug fix.
 *				
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
#include "gwig.c"

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
/*****************************************************************************/

static void GWigPass_2nd(struct gwig *pWig, double *X)
{
  int Nstep = pWig->PN*(pWig->Nw);
  double dl    = pWig->Lw/(pWig->PN);

  for (int i = 1; i <= Nstep; i++) {
    GWigMap_2nd(pWig, X, dl);
  }
}

static void GWigPass_4th(struct gwig *pWig, double *X)
{
  int Nstep = pWig->PN*(pWig->Nw);
  double dl = pWig->Lw/(pWig->PN);
  double dl1 = dl*KICK1;
  double dl0 = dl*KICK2;

  for (int i = 1; i <= Nstep; i++ ) {
    GWigMap_2nd(pWig, X, dl1);
    GWigMap_2nd(pWig, X, dl0);
    GWigMap_2nd(pWig, X, dl1);
  }
}

#define second 2
#define fourth 4
void GWigSymplecticPass(double *r, double gamma, double Ltot, double Lw,
            double Bmax, int Nstep, int Nmeth, int NHharm, int NVharm,
            double *By, double *Bx, double *T1, double *T2,
            double *R1, double *R2, int num_particles)
{
    int c;
    double *r6;
    struct gwig pWig;
    /* Energy is defined in the lattice in eV but GeV is used by the gwig code. */
    GWigInit2(&pWig, gamma,Ltot, Lw, Bmax, Nstep, Nmeth, NHharm, NVharm,0, 0, By,Bx,T1,T2,R1,R2);

    for (c = 0;c<num_particles;c++) {
        r6 = r+c*6;
        pWig.Zw = 0.0;
        if (!atIsNaN(r6[0])) {
			/* Misalignment at entrance */
			if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);
            switch (Nmeth) {
                case second:
                    GWigPass_2nd(&pWig, r6);
                    break;
                case fourth:
                    GWigPass_4th(&pWig, r6);
                    break;
                default:
                    printf("Invalid wiggler integration method %d.\n", Nmeth);
                    break;
            }
			/* Misalignment at exit */
            if (R2) ATmultmv(r6,R2);
            if (T2) ATaddvv(r6,T2);
        }
    }
}

/*****************************************************************************/
/********** END PHYSICS SECTION **********************************************/
/*****************************************************************************/

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)

{
    double gamma;

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
        Energy=atGetOptionalDouble(ElemData,"Energy",Param->energy); check_error();
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

    GWigSymplecticPass(r_in, gamma, Elem->Length, Elem->Lw, Elem->Bmax,
            Elem->Nstep, Elem->Nmeth, Elem->NHharm, Elem->NVharm,
            Elem->By, Elem->Bx, Elem->T1, Elem->T2, Elem->R1, Elem->R2,
            num_particles);
    return Elem;
}


/********** END WINDOWS DLL GATEWAY SECTION **********************************/

/********** MATLAB GATEWAY ***************************************************/

MODULE_DEF(GWigSymplecticPass)        /* Dummy module initialisation */

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
        double *By, *Bx;
        double *R1, *R2, *T1, *T2;
        double Ltot, Lw, Bmax, Energy;
        int Nstep, Nmeth;
        int NHharm, NVharm;
        int num_particles = mxGetN(prhs[1]);
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
        Energy=atGetOptionalDouble(ElemData,"Energy",0.0); check_error();
        R1 = atGetOptionalDoubleArray(ElemData, "R1"); check_error();
        R2 = atGetOptionalDoubleArray(ElemData, "R2"); check_error();
        T1 = atGetOptionalDoubleArray(ElemData, "T1"); check_error();
        T2 = atGetOptionalDoubleArray(ElemData, "T2"); check_error();
        if (nrhs > 2) atProperties(prhs[2], &Energy, &rest_energy, &charge);

        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        Gamma = atGamma(Energy, Energy, rest_energy);
        r_in = mxGetDoubles(plhs[0]);

        GWigSymplecticPass(r_in, Gamma, Ltot, Lw, Bmax, Nstep, Nmeth, NHharm,
            NVharm, By, Bx, T1, T2, R1, R2, num_particles);
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
