#include "atelem.c"
#include "atlalib.c"

struct elem {
    double Length;
    double InvRho;
    /* Optional fields */
    double KxKz;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
};

/******************************************************************************/
/* PHYSICS SECTION ************************************************************/

static void matfoc(double L, double Knorm, double *MD, double *M12, double *M21)
{
   double g  = fabs(Knorm);
   double t  = sqrt(g);
   double lt = L*t;
				
   if (Knorm == 0) {
      *MD = 1.0;
      *M12 = L;
      *M21 = 0;
   }
   else if (Knorm > 0) {
	*MD = cos(lt);
	*M12 = sin(lt)/t;
	*M21 = -*M12*g;
   }
   else  {
	*MD = cosh(lt);
	*M12 = sinh(lt)/t;
	*M21 = *M12*g;
   }
}

static void foc6(double *r, double L, double Kx, double Kz)
{
   double M12,M21,M34,M43,MVD,MHD;  /* non-0 elements of transfer matrix */
   double p_norm = 1/(1+r[4]);
   double x   = r[0];
   double xpr = r[1]*p_norm;
   double y   = r[2];
   double ypr = r[3]*p_norm;

   (void) matfoc(L, Kx*p_norm, &MHD, &M12, &M21);
   (void) matfoc(L, Kz*p_norm, &MVD, &M34, &M43);

   r[0]=  MHD*x + M12*xpr;
   r[1]= (M21*x + MHD*xpr)/p_norm;
   r[2]=  MVD*y + M34*ypr;
   r[3]= (M43*y + MVD*ypr)/p_norm;

    /* no change in r[4] (delta) */

   r[5]+= (fabs(Kx)*p_norm*x*x*(L-MHD*M12) - fabs(Kz)*p_norm*y*y*(L-MVD*M34))/4;
   r[5]+= (xpr*xpr*(L+MHD*M12)+ypr*ypr*(L+MVD*M34))/4;
   r[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2;
}

void WiggLinearPass(double *r, double le, double invrho, double kxkz, double *T1, double *T2, double *R1, double *R2, int num_particles)
{
    int c;
    double *r6;
    double kz = 0.5/(1.0+kxkz)*invrho*invrho;
    double kx = kxkz*kz;
    
    for (c = 0;c<num_particles;c++) {
        r6 = r+c*6;
        if (!atIsNaN(r6[0])) {
            /* Misalignment at entrance */
            if (T1 != NULL) ATaddvv(r6,T1);
            if (R1 != NULL) ATmultmv(r6,R1);
            
            foc6(r6, le, kx, kz);
            
            /* Misalignment at exit */
            if (T2 != NULL) ATmultmv(r6,R2);
            if (R2 != NULL) ATaddvv(r6,T2);
        }
    }
}

/********** END PHYSICS SECTION ***********************************************/
/******************************************************************************/

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double Length, InvRho, KxKz;
        double*R1, *R2, *T1, *T2;
        Length=atGetDouble(ElemData,"Length"); check_error();
        InvRho=atGetDouble(ElemData,"InvRho"); check_error();
        /*optional fields*/
        KxKz=atGetOptionalDouble(ElemData,"KxKz", 0.0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->InvRho=InvRho;
        /*optional fields*/
        Elem->KxKz=KxKz;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
    }
    WiggLinearPass(r_in, Elem->Length, Elem->InvRho, Elem->KxKz,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2 , num_particles);
    return Elem;
}

MODULE_DEF(WiggLinearPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#ifdef MATLAB_MEX_FILE
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Length, InvRho, KxKz;
        double*R1, *R2, *T1, *T2;
        Length=atGetDouble(ElemData,"Length"); check_error();
        InvRho=atGetDouble(ElemData,"InvRho"); check_error();
        /*optional fields*/
        KxKz=atGetOptionalDouble(ElemData,"KxKz", 0.0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        WiggLinearPass(r_in, Length, InvRho, KxKz, T1, T2, R1, R2 , num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(2,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("InvRho"));
        
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(5,1);
            mxSetCell(plhs[1],0,mxCreateString("KxKz"));
            mxSetCell(plhs[1],1,mxCreateString("T1"));
            mxSetCell(plhs[1],2,mxCreateString("T2"));
            mxSetCell(plhs[1],3,mxCreateString("R1"));
            mxSetCell(plhs[1],4,mxCreateString("R2"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*MATLAB_MEX_FILE*/
