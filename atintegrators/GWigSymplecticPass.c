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

void GWigInit(struct gwig *Wig, double design_energy, double Ltot, double Lw,
            double Bmax, int Nstep, int Nmeth, int NHharm, int NVharm,
            double *By, double *Bx, double *T1, double *T2, double *R1,
            double *R2)
{
    double *tmppr;
    int    i;
    double kw;

    Wig->E0 = design_energy;
    Wig->Pmethod = Nmeth;
    Wig->PN = Nstep;
    Wig->Nw = (int)(Ltot / Lw);
    Wig->NHharm = NHharm;
    Wig->NVharm = NVharm;
    Wig->PB0 = Bmax;
    Wig->Lw = Lw;

    kw = 2.0e0*PI/(Wig->Lw);
    Wig->Zw = 0.0;
    Wig->Aw = 0.0;
    tmppr = By;
    for (i = 0; i < NHharm; i++) {
        tmppr++;
        Wig->HCw[i] = 0.0;
        Wig->HCw_raw[i] = *tmppr;
        tmppr++;
        Wig->Hkx[i] = (*tmppr) * kw;
        tmppr++;
        Wig->Hky[i] = (*tmppr) * kw;
        tmppr++;
        Wig->Hkz[i] = (*tmppr) * kw;
        tmppr++;
        Wig->Htz[i] = *tmppr;
        tmppr++;
    }
    tmppr = Bx;
    for (i = 0; i < NVharm; i++) {
        tmppr++;
        Wig->VCw[i] = 0.0;
        Wig->VCw_raw[i] = *tmppr;
        tmppr++;
        Wig->Vkx[i] = (*tmppr) * kw;
        tmppr++;
        Wig->Vky[i] = (*tmppr) * kw;
        tmppr++;
        Wig->Vkz[i] = (*tmppr) * kw;
        tmppr++;
        Wig->Vtz[i] = *tmppr;
        tmppr++;
    }
    for (i = NHharm; i< WHmax; i++) {
        Wig->HCw[i] = 0.0;
        Wig->HCw_raw[i] = 0.0;
        Wig->Hkx[i] = 0.0;
        Wig->Hky[i] = 0.0;
        Wig->Hkz[i] = 0.0;
        Wig->Htz[i] = 0.0;
    }
    for (i = NVharm; i< WHmax; i++) {
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
void GWigSymplecticPass(double *r, double Energy, double Ltot, double Lw,
            double Bmax, int Nstep, int Nmeth, int NHharm, int NVharm,
            double *By, double *Bx, double *T1, double *T2, double *R1,
            double *R2, int num_particles)
{
    int c;
    double *r6;
    struct gwig Wig;
    /* Energy is defined in the lattice in eV but GeV is used by the gwig code. */
    Energy = Energy / 1e9;

    GWigInit(&Wig, Energy, Ltot, Lw, Bmax, Nstep, Nmeth, NHharm, NVharm, By, Bx, T1, T2, R1, R2);

    for (c = 0;c<num_particles;c++) {
        r6 = r+c*6;
        if (!atIsNaN(r6[0])) {
            switch (Nmeth) {
                case second:
                    GWigPass_2nd(&Wig, r6);
                    break;
                case fourth:
                    GWigPass_4th(&Wig, r6);
                    break;
                default:
                    printf("Invalid wiggler integration method %d.\n", Nmeth);
                    break;
            }
        }
    }
}

/********** END PHYSICS SECTION **********************************************/
/*****************************************************************************/

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)

{
    if (!Elem) {
        double *R1, *R2, *T1, *T2;
        double *By, *Bx;
        double Ltot, Lw, Bmax, Energy;
        int Nstep, Nmeth;
        int NHharm, NVharm;

        Energy = atGetDouble(ElemData, "Energy"); check_error();
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
    GWigSymplecticPass(r_in, Elem->Energy, Elem->Length, Elem->Lw, Elem->Bmax,
            Elem->Nstep, Elem->Nmeth, Elem->NHharm, Elem->NVharm, Elem->By,
            Elem->Bx, Elem->T1, Elem->T2, Elem->R1, Elem->R2, num_particles);
    return Elem;
}


/********** END WINDOWS DLL GATEWAY SECTION **********************************/

/********** MATLAB GATEWAY ***************************************************/

MODULE_DEF(GWigSymplecticPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(       int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double *By, *Bx;
        double *R1, *R2, *T1, *T2;
        double Ltot, Lw, Bmax, Energy;
        int Nstep, Nmeth;
        int NHharm, NVharm;

        Energy = atGetDouble(ElemData, "Energy"); check_error();
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
        R1 = atGetOptionalDoubleArray(ElemData, "R1"); check_error();
        R2 = atGetOptionalDoubleArray(ElemData, "R2"); check_error();
        T1 = atGetOptionalDoubleArray(ElemData, "T1"); check_error();
        T2 = atGetOptionalDoubleArray(ElemData, "T2"); check_error();
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        GWigSymplecticPass(r_in, Energy, Ltot, Lw, Bmax, Nstep, Nmeth, NHharm,
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
