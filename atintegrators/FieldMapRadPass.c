/* 
 *  FieldMapRadPass.c
 *  PassMethod for Accelerator Toolbox to track throug a magnet described with a Field Map with radiation
 *
 *  10/07/2018
 *  Nicola Carmignani
 */

#include "atelem.c"
#include "atlalib.c"
#include "driftkickFM.c"		/* fastdrift.c, strthinkick.c */
#include<math.h>

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656


void FieldMapRadPass(double *r,struct elem *Elem,double Brho,double irho,int num_particles)
{
    int c,m;
    double norm, NormL1, NormL2;
    double *r6;
    double SL, L1, L2, K1, K2; /*K1divBrho, K2divBrho;*/
    /*double Brho=Elem->Energy/299792458;*/
    SL = Elem->Length/Elem->NumIntSteps;
    L1 = SL*DRIFT1;
    L2 = SL*DRIFT2;
    K1 = SL*KICK1;
    K2 = SL*KICK2;
    /*K1divBrho = SL*KICK1/Brho;
    K2divBrho = SL*KICK2/Brho;*/
    for (c = 0;c<num_particles;c++)	{   /*Loop over particles  */
        r6 = r+c*6;
        if(!atIsNaN(r6[0])) {
            /*  misalignment at entrance  */
            if (Elem->T1) ATaddvv(r6,Elem->T1);
            if (Elem->R1) ATmultmv(r6,Elem->R1);
            /* Check physical apertures at the entrance of the magnet */
            if (Elem->RApertures) checkiflostRectangularAp(r6,Elem->RApertures);
            if (Elem->EApertures) checkiflostEllipticalAp(r6,Elem->EApertures);
            /*  integrator  */
            if (Elem->BendingAngle)
            {
                for (m=0; m < Elem->NumIntSteps; m++) {  /*  Loop over slices */
                    r6 = r+c*6;
                    ATdrift6(r6,L1);
                    bndthinkickFMrad(r6, Elem, K1, Brho, irho);
                    ATdrift6(r6,L2);
                    bndthinkickFMrad(r6, Elem, K2, Brho, irho);
                    ATdrift6(r6,L2);
                    bndthinkickFMrad(r6, Elem, K1, Brho, irho);
                    ATdrift6(r6,L1);
                }
            }
            else
            {
                for (m=0; m < Elem->NumIntSteps; m++) {  /*  Loop over slices */
                    r6 = r+c*6;
                    ATdrift6(r6,L1);
                    strthinkickFMrad(r6, Elem, K1, Brho);
                    ATdrift6(r6,L2);
                    strthinkickFMrad(r6, Elem, K2, Brho);
                    ATdrift6(r6,L2);
                    strthinkickFMrad(r6, Elem, K1, Brho);
                    ATdrift6(r6,L1);
                }
            }
            /* Check physical apertures at the exit of the magnet */
            if (Elem->RApertures) checkiflostRectangularAp(r6,Elem->RApertures);
            if (Elem->EApertures) checkiflostEllipticalAp(r6,Elem->EApertures);
            /* Misalignment at exit */
            if (Elem->R2) ATmultmv(r6,Elem->R2);
            if (Elem->T2) ATaddvv(r6,Elem->T2);
        }
    }
}



#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    double Brho, irho;
    if (!Elem) {
        double Length;
        int NumIntSteps, Nx, Ny;
        double *LUT_Bx, *LUT_By, *LUT_dBxdx, *LUT_dBydx, *LUT_dBxdy, *LUT_dBydy;
        double *LUT_d2Bxdxdx, *LUT_d2Bxdxdy, *LUT_d2Bxdydy, *LUT_d2Bydxdx, *LUT_d2Bydxdy, *LUT_d2Bydydy;
        double *X, *Y, xmin, ymin, delta_x, delta_y, Energy, BendingAngle;
        double *R1, *R2, *T1, *T2, *EApertures, *RApertures, *KickAngle;
        Length=atGetDouble(ElemData,"Length"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        LUT_Bx=atGetDoubleArray(ElemData,"LUT_Bx"); check_error();
        LUT_By=atGetDoubleArray(ElemData,"LUT_By"); check_error();
        LUT_dBxdx=atGetDoubleArray(ElemData,"LUT_dBxdx"); check_error();
        LUT_dBydx=atGetDoubleArray(ElemData,"LUT_dBydx"); check_error();
        LUT_dBxdy=atGetDoubleArray(ElemData,"LUT_dBxdy"); check_error();
        LUT_dBydy=atGetDoubleArray(ElemData,"LUT_dBydy"); check_error();
        LUT_d2Bxdxdx=atGetDoubleArray(ElemData,"LUT_d2Bxdxdx"); check_error();
        LUT_d2Bxdxdy=atGetDoubleArray(ElemData,"LUT_d2Bxdxdy"); check_error();
        LUT_d2Bxdydy=atGetDoubleArray(ElemData,"LUT_d2Bxdydy"); check_error();
        LUT_d2Bydxdx=atGetDoubleArray(ElemData,"LUT_d2Bydxdx"); check_error();
        LUT_d2Bydxdy=atGetDoubleArray(ElemData,"LUT_d2Bydxdy"); check_error();
        LUT_d2Bydydy=atGetDoubleArray(ElemData,"LUT_d2Bydydy"); check_error();
        X=atGetDoubleArray(ElemData,"X"); check_error();
        Y=atGetDoubleArray(ElemData,"Y"); check_error();
        Nx=atGetLong(ElemData,"Nx"); check_error();
        Ny=atGetLong(ElemData,"Ny"); check_error();
        xmin=atGetDouble(ElemData,"xmin"); check_error();
        ymin=atGetDouble(ElemData,"ymin"); check_error();
        delta_x=atGetDouble(ElemData,"delta_x"); check_error();
        delta_y=atGetDouble(ElemData,"delta_y"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        /*optional fields*/
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();
        BendingAngle=atGetOptionalDouble(ElemData,"BendingAngle",0); check_error();
  
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->NumIntSteps=NumIntSteps;
        Elem->LUT_Bx=LUT_Bx;
        Elem->LUT_By=LUT_By;
        Elem->LUT_dBxdx=LUT_dBxdx;
        Elem->LUT_dBydx=LUT_dBydx;
        Elem->LUT_dBxdy=LUT_dBxdy;
        Elem->LUT_dBydy=LUT_dBydy;
        Elem->LUT_d2Bxdxdx=LUT_d2Bxdxdx;
        Elem->LUT_d2Bxdxdy=LUT_d2Bxdxdy;
        Elem->LUT_d2Bxdydy=LUT_d2Bxdydy;
        Elem->LUT_d2Bydxdx=LUT_d2Bydxdx;
        Elem->LUT_d2Bydxdy=LUT_d2Bydxdy;
        Elem->LUT_d2Bydydy=LUT_d2Bydydy;
        Elem->X=X;
        Elem->Y=Y;
        Elem->Nx=Nx;
        Elem->Ny=Ny;
        Elem->xmin=xmin;
        Elem->ymin=ymin;
        Elem->delta_x=delta_x;
        Elem->delta_y=delta_y;
        Elem->Energy=Energy;
        /*optional fields*/
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
        Elem->KickAngle=KickAngle;
        Elem->BendingAngle=BendingAngle;
    }
    Brho=Elem->Energy/299792458;
    irho = Elem->BendingAngle/Elem->Length;
    FieldMapRadPass(r_in,Elem,Brho,irho,num_particles);
    return Elem;
}

MODULE_DEF(FieldMapRadPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Length, Brho, irho, BendingAngle;
        int NumIntSteps, Nx, Ny;
        double *LUT_Bx, *LUT_By, *LUT_dBxdx, *LUT_dBydx, *LUT_dBxdy, *LUT_dBydy;
        double *LUT_d2Bxdxdx, *LUT_d2Bxdxdy, *LUT_d2Bxdydy, *LUT_d2Bydxdx, *LUT_d2Bydxdy, *LUT_d2Bydydy;
        double *X, *Y, xmin, ymin, delta_x, delta_y, Energy;
        double *R1, *R2, *T1, *T2, *EApertures, *RApertures, *KickAngle;
        Length=atGetDouble(ElemData,"Length"); check_error();
        NumIntSteps=atGetLong(ElemData,"NumIntSteps"); check_error();
        LUT_Bx=atGetDoubleArray(ElemData,"LUT_Bx"); check_error();
        LUT_By=atGetDoubleArray(ElemData,"LUT_By"); check_error();
        LUT_dBxdx=atGetDoubleArray(ElemData,"LUT_dBxdx"); check_error();
        LUT_dBydx=atGetDoubleArray(ElemData,"LUT_dBydx"); check_error();
        LUT_dBxdy=atGetDoubleArray(ElemData,"LUT_dBxdy"); check_error();
        LUT_dBydy=atGetDoubleArray(ElemData,"LUT_dBydy"); check_error();
        LUT_d2Bxdxdx=atGetDoubleArray(ElemData,"LUT_d2Bxdxdx"); check_error();
        LUT_d2Bxdxdy=atGetDoubleArray(ElemData,"LUT_d2Bxdxdy"); check_error();
        LUT_d2Bxdydy=atGetDoubleArray(ElemData,"LUT_d2Bxdydy"); check_error();
        LUT_d2Bydxdx=atGetDoubleArray(ElemData,"LUT_d2Bydxdx"); check_error();
        LUT_d2Bydxdy=atGetDoubleArray(ElemData,"LUT_d2Bydxdy"); check_error();
        LUT_d2Bydydy=atGetDoubleArray(ElemData,"LUT_d2Bydydy"); check_error();
        X=atGetDoubleArray(ElemData,"X"); check_error();
        Y=atGetDoubleArray(ElemData,"Y"); check_error();
        Nx=atGetLong(ElemData,"Nx"); check_error();
        Ny=atGetLong(ElemData,"Ny"); check_error();
        xmin=atGetDouble(ElemData,"xmin"); check_error();
        ymin=atGetDouble(ElemData,"ymin"); check_error();
        delta_x=atGetDouble(ElemData,"delta_x"); check_error();
        delta_y=atGetDouble(ElemData,"delta_y"); check_error();
        Energy=atGetDouble(ElemData,"Energy"); check_error();
        /*optional fields*/
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        KickAngle=atGetOptionalDoubleArray(ElemData,"KickAngle"); check_error();
        BendingAngle=atGetOptionalDouble(ElemData,"BendingAngle",0); check_error();
        
        struct elem *Elem;
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->NumIntSteps=NumIntSteps;
        Elem->LUT_Bx=LUT_Bx;
        Elem->LUT_By=LUT_By;
        Elem->LUT_dBxdx=LUT_dBxdx;
        Elem->LUT_dBydx=LUT_dBydx;
        Elem->LUT_dBxdy=LUT_dBxdy;
        Elem->LUT_dBydy=LUT_dBydy;
        Elem->LUT_d2Bxdxdx=LUT_d2Bxdxdx;
        Elem->LUT_d2Bxdxdy=LUT_d2Bxdxdy;
        Elem->LUT_d2Bxdydy=LUT_d2Bxdydy;
        Elem->LUT_d2Bydxdx=LUT_d2Bydxdx;
        Elem->LUT_d2Bydxdy=LUT_d2Bydxdy;
        Elem->LUT_d2Bydydy=LUT_d2Bydydy;
        Elem->X=X;
        Elem->Y=Y;
        Elem->Nx=Nx;
        Elem->Ny=Ny;
        Elem->xmin=xmin;
        Elem->ymin=ymin;
        Elem->delta_x=delta_x;
        Elem->delta_y=delta_y;
        Elem->Energy=Energy;
        /*optional fields*/
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
        Elem->KickAngle=KickAngle;
        Elem->BendingAngle=BendingAngle;
        
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        Brho=Elem->Energy/299792458;
        irho = Elem->BendingAngle/Elem->Length;
        FieldMapRadPass(r_in,Elem,Brho,irho,num_particles);
        
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(23,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("NumIntSteps"));
        mxSetCell(plhs[0],2,mxCreateString("LUT_Bx"));
        mxSetCell(plhs[0],3,mxCreateString("LUT_By"));
        mxSetCell(plhs[0],4,mxCreateString("LUT_dBxdx"));
        mxSetCell(plhs[0],5,mxCreateString("LUT_dBydx"));
        mxSetCell(plhs[0],6,mxCreateString("LUT_dBxdy"));
        mxSetCell(plhs[0],7,mxCreateString("LUT_dBydy"));
        mxSetCell(plhs[0],8,mxCreateString("LUT_d2Bxdxdx"));
        mxSetCell(plhs[0],9,mxCreateString("LUT_d2Bxdxdy"));
        mxSetCell(plhs[0],10,mxCreateString("LUT_d2Bxdydy"));
        mxSetCell(plhs[0],11,mxCreateString("LUT_d2Bydxdx"));
        mxSetCell(plhs[0],12,mxCreateString("LUT_d2Bydxdy"));
        mxSetCell(plhs[0],13,mxCreateString("LUT_d2Bydydy"));
        mxSetCell(plhs[0],14,mxCreateString("X"));
        mxSetCell(plhs[0],15,mxCreateString("Y"));
        mxSetCell(plhs[0],16,mxCreateString("Nx"));
        mxSetCell(plhs[0],17,mxCreateString("Ny"));
        mxSetCell(plhs[0],18,mxCreateString("xmin"));
        mxSetCell(plhs[0],19,mxCreateString("ymin"));
        mxSetCell(plhs[0],20,mxCreateString("delta_x"));
        mxSetCell(plhs[0],21,mxCreateString("delta_y"));
        mxSetCell(plhs[0],22,mxCreateString("Energy"));
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(6,1);
            mxSetCell(plhs[1], 0,mxCreateString("T1"));
            mxSetCell(plhs[1], 1,mxCreateString("T2"));
            mxSetCell(plhs[1], 2,mxCreateString("R1"));
            mxSetCell(plhs[1], 3,mxCreateString("R2"));
            mxSetCell(plhs[1], 4,mxCreateString("RApertures"));
            mxSetCell(plhs[1], 5,mxCreateString("EApertures"));
            /*mxSetCell(plhs[1], 6,mxCreateString("KickAngle"));*/
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*MATLAB_MEX_FILE*/
