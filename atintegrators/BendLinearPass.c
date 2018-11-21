/* BendLinearPass.c
 * Accelerator Toolbox
 * Revision 6/26/00
 * A.Terebilo terebilo@ssrl.slac.stanford.edu
 */


#include "atelem.c"
#include "atlalib.c"
#include "atphyslib.c"


#define SQR(X) X*X
struct elem
{
    double Length;
    double BendingAngle;
    double EntranceAngle;
    double ExitAngle;
    /* Optional fields */
    double K;
    double ByError;
    double FringeInt1;
    double FringeInt2;
    double FullGap;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
};

void bend6(double* r, double L, double b_angle, double grd, double ByError)
{
    double M12,M21,M34,M43,MVD,MHD;  /*  non-0 elements of transfer matrix */
    double x, xpr,y,ypr,delta;
    double sqrtG1, sqrtG2, arg1, arg2;
    double Kx = b_angle/L; /* Curvature of the design trajectory */
    double p_norm = 1/(1+r[4]);
    double G1 = (Kx*Kx+grd)*p_norm;
    double G2 = -grd*p_norm;
    /* Horizontal Transverse Matrix */
    if(G1==0)		/* Special case Kx^2 + grd = 0 */
    {
        MHD = 1;
        M12 = L;
        M21 = 0;
    }
    else
    {
        if(G1 > 0)
        {
            sqrtG1 = sqrt(G1);
            arg1 = L*sqrtG1;
            MHD = cos(arg1);
            M12 = sin(arg1)/sqrtG1;
            M21 = -sin(arg1)*sqrtG1;
        }
        else
        {
            sqrtG1 = sqrt(-G1);
            arg1 = L*sqrtG1;
            MHD = cosh(arg1);
            M12 = sinh(arg1)/sqrtG1;
            M21 = sinh(arg1)*sqrtG1;
        }
    }
    /*  Vertical Transverse Matrix */
    if(G2==0) /*  No gradient - vertical motion is a drift  */
    {
        MVD = 1;
        M34 = L;
        M43 = 0;
    }
    else
    {
        if(G2 > 0)	/* Vertical focusing */
        {
            sqrtG2 = sqrt(G2);
            arg2 = L*sqrtG2;
            MVD = cos(arg2);;
            M34 = sin(arg2)/sqrtG2;
            M43 = -sin(arg2)*sqrtG2;
        }
        else		/*  Vertical defocusing	*/
        {
            sqrtG2 = sqrt(-G2);
            arg2 = L*sqrtG2;
            MVD = cosh(arg2);
            M34 = sinh(arg2)/sqrtG2;
            M43 = sinh(arg2)*sqrtG2;
        }
    }
    x   = r[0];
    xpr = r[1]*p_norm;
    y   = r[2];
    ypr = r[3]*p_norm;
    delta = r[4];
    r[0]=  MHD*x + M12*xpr ;
    r[1]= (M21*x + MHD*xpr)/p_norm;
    if(G1==0)
    {
        r[0]+= (delta*p_norm-ByError)*L*L*Kx/2;
        r[1]+= (delta*p_norm-ByError)*L*Kx/p_norm ;
    }
    else
    {
        if(G1>0)
        {
            r[0]+= (delta*p_norm-ByError)*(1-cos(arg1))*Kx/G1;
            r[1]+= (delta*p_norm-ByError)*sin(arg1)*Kx/(sqrtG1*p_norm) ;
        }
        else
        {
            r[0]+= (delta*p_norm-ByError)*(1-cosh(arg1))*Kx/G1;
            r[1]+= (delta*p_norm-ByError)*sinh(arg1)*Kx/(sqrtG1*p_norm) ;
        }
    }
    r[2]=  MVD*y + M34*ypr;
    r[3]= (M43*y + MVD*ypr)/p_norm ;
    r[5]+= xpr*xpr*(L+MHD*M12)/4;
    if (G1==0) {
        /* Do nothing */
    }
    else
    {
        r[5]+= (L-MHD*M12)*(x*x*G1+(delta*p_norm-ByError)*(delta*p_norm-ByError)*Kx*Kx/G1
                -2*x*Kx*(delta*p_norm-ByError))/4;
        r[5]+= M12*M21*( x*xpr - xpr*(delta*p_norm-ByError)*Kx/G1)/2;
        r[5]+= Kx*x*M12  +   xpr*(1-MHD)*Kx/G1   +   (delta*p_norm-ByError)*(L-M12)*Kx*Kx/G1;
    }
    r[5]+= ((L-MVD*M34)*y*y*G2 + ypr*ypr*(L+MVD*M34))/4;
    r[5]+= M34*M43*x*xpr/2;
}

void BendLinearPass(double *r, double le, double grd ,double ba, double bye,
        double entrance_angle, double exit_angle,
        double fint1, double fint2, double gap,
        double *T1, double *T2, double *R1, double *R2,
        int num_particles)
        /* Set T1, T2, R1, R2 to NULL pointers to ignore misalignmets
         * Set fint OR gap are 0 to ignore fringe effects
         * Set bye to 0 to ignore ByError*/
{
    double *r6;
    double irho = ba/le;
    int c;
    
    #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(shared) shared(r,num_particles) private(c,r6)
    for (c = 0;c<num_particles;c++) {
        r6 = r+c*6;
        if(!atIsNaN(r6[0]) && atIsFinite(r6[4]))
        /*
         * function bend6 internally calculates the square root
         * of the energy deviation of the particle
         * To protect against DOMAIN and OVERFLOW error, check if the
         * fifth component of the phase spacevector r6[4] is finite
         */
        {
            /* Misalignment at entrance */
            if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);
            /* edge focus */
            edge_fringe_entrance(r6, irho, entrance_angle, fint1, gap, 1);
            bend6(r6, le, ba, grd, bye);
            /* edge focus */
            edge_fringe_exit(r6, irho, exit_angle, fint2, gap, 1);
            /* Misalignment at exit */
            if (R2) ATmultmv(r6,R2);
            if (T2) ATaddvv(r6,T2);
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double Length, BendingAngle, EntranceAngle, ExitAngle, K, ByError, FringeInt1, FringeInt2, FullGap;
        double *R1, *R2, *T1, *T2;
        Length=atGetDouble(ElemData,"Length"); check_error();
        BendingAngle=atGetDouble(ElemData,"BendingAngle"); check_error();
        EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        /*optional fields*/
        K=atGetOptionalDouble(ElemData,"K",0); check_error();
        ByError=atGetOptionalDouble(ElemData,"ByError",0); check_error();
        FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0); check_error();
        FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0); check_error();
        FullGap=atGetOptionalDouble(ElemData,"FullGap",0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->BendingAngle=BendingAngle;
        Elem->EntranceAngle=EntranceAngle;
        Elem->ExitAngle=ExitAngle;
        Elem->K=K;
        Elem->ByError=ByError;
        Elem->FringeInt1=FringeInt1;
        Elem->FringeInt2=FringeInt2;
        Elem->FullGap=FullGap;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
    }
    BendLinearPass(r_in, Elem->Length,Elem->K,Elem->BendingAngle,Elem->ByError,
            Elem->EntranceAngle,Elem->ExitAngle,Elem->FringeInt1,Elem->FringeInt2,
            Elem->FullGap,Elem->T1,Elem->T2,Elem->R1,Elem->R2,num_particles);
    return(Elem);
}

MODULE_DEF(BendLinearPass)        /* Dummy module initialisation */
#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2 ) {
        double Length, BendingAngle, EntranceAngle, ExitAngle, K, ByError, FringeInt1, FringeInt2, FullGap;
        double *R1, *R2, *T1, *T2;
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        Length=atGetDouble(ElemData,"Length"); check_error();
        BendingAngle=atGetDouble(ElemData,"BendingAngle"); check_error();
        EntranceAngle=atGetDouble(ElemData,"EntranceAngle"); check_error();
        ExitAngle=atGetDouble(ElemData,"ExitAngle"); check_error();
        /*optional fields*/
        K=atGetOptionalDouble(ElemData,"K",0); check_error();
        ByError=atGetOptionalDouble(ElemData,"ByError",0); check_error();
        FringeInt1=atGetOptionalDouble(ElemData,"FringeInt1",0); check_error();
        FringeInt2=atGetOptionalDouble(ElemData,"FringeInt2",0); check_error();
        FullGap=atGetOptionalDouble(ElemData,"FullGap",0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        BendLinearPass(r_in,Length,K,BendingAngle,ByError,
                EntranceAngle,ExitAngle,FringeInt1,FringeInt2,
                FullGap,T1,T2,R1,R2,num_particles);
    }
    else
    {   /* return list of required fields */
        plhs[0] = mxCreateCellMatrix(4,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("BendingAngle"));
        mxSetCell(plhs[0],2,mxCreateString("EntranceAngle"));
        mxSetCell(plhs[0],3,mxCreateString("ExitAngle"));
        
        if(nlhs>1) /* Required and optional fields */
        {   plhs[1] = mxCreateCellMatrix(9,1);
            mxSetCell(plhs[1],0,mxCreateString("K"));
            mxSetCell(plhs[1],1,mxCreateString("ByError"));
            mxSetCell(plhs[1],2,mxCreateString("FullGap"));
            mxSetCell(plhs[1],3,mxCreateString("FringeInt1"));
            mxSetCell(plhs[1],4,mxCreateString("FringeInt2"));
            mxSetCell(plhs[1],5,mxCreateString("T1"));
            mxSetCell(plhs[1],6,mxCreateString("T2"));
            mxSetCell(plhs[1],7,mxCreateString("R1"));
            mxSetCell(plhs[1],8,mxCreateString("R2"));
        }
    }
}
#endif /*MATLAB_MEX_FILE*/
