/* CorrectorPass.c 
   Accelerator Toolbox 
   Revision 10/11/01
   A.Terebilo
   CorrectorPass expects an element to have a fields:
   Length, KickAngle
*/

#include "atelem.c"
#include "atlalib.c"

struct elem
{
    double Length;
    double Scaling;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *RApertures;
    double *EApertures;
    double *KickAngle;
};

void CorrectorPass(double *r, double xkick, double ykick, double le,
        double *T1, double *T2,
        double *R1, double *R2,
        double *RApertures, double *EApertures,
        double scaling, int num_particles)
/* xkick, ykick - horizontal and vertical kicks in radians
   r - 6-by-N matrix of initial conditions reshaped into 1-d array of 6*N elements
*/
{
	if (le == 0.0)
        for (int c = 0; c<num_particles; c++) { /* Loop over particles */
            double *r6 = r + 6*c;
		    if (!atIsNaN(r6[0])) {
                /* Check for change of reference momentum */
                if (scaling != 1.0) ATChangePRef(r6, scaling);
                /*  misalignment at entrance  */
                if (T1) ATaddvv(r6,T1);
                if (R1) ATmultmv(r6,R1);
                /* Check physical apertures at the entrance of the magnet */
                if (RApertures) checkiflostRectangularAp(r6,RApertures);
                if (EApertures) checkiflostEllipticalAp(r6,EApertures);
		        r6[1] += xkick;
   		        r6[3] += ykick;
                /* Misalignment at exit */
                if (R2) ATmultmv(r6,R2);
                if (T2) ATaddvv(r6,T2);
                 /* Check for change of reference momentum */
                if (scaling != 1.0) ATChangePRef(r6, 1.0/scaling);
  		    }
		}	
	else
        #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(none) \
        shared(r,num_particles,le,xkick,ykick,T1,T2,R1,R2,EApertures,RApertures,scaling)
        for (int c = 0; c<num_particles; c++) { /* Loop over particles */
            double *r6 = r + 6*c;
		    if(!atIsNaN(r6[0])) {
			    double p_norm, NormL;
                /* Check for change of reference momentum */
                if (scaling != 1.0) ATChangePRef(r6, scaling);
                p_norm = 1.0/(1.0+r6[4]);
                NormL = le*p_norm;
                /*  misalignment at entrance  */
                if (T1) ATaddvv(r6,T1);
                if (R1) ATmultmv(r6,R1);
                /* Check physical apertures at the entrance of the magnet */
                if (RApertures) checkiflostRectangularAp(r6,RApertures);
                if (EApertures) checkiflostEllipticalAp(r6,EApertures);
	            r6[5] += NormL*p_norm*(xkick*xkick/3 + ykick*ykick/3 +
   		            r6[1]*r6[1] + r6[3]*r6[3] +
   		            r6[1]*xkick + r6[3]*ykick)/2;

			    r6[0] += NormL*(r6[1]+xkick/2);
		        r6[1] += xkick;
		        r6[2] += NormL*(r6[3]+ykick/2);
   		        r6[3] += ykick;
                /* Check physical apertures at the exit of the magnet */
                if (RApertures) checkiflostRectangularAp(r6,RApertures);
                if (EApertures) checkiflostEllipticalAp(r6,EApertures);
                /* Misalignment at exit */
                if (R2) ATmultmv(r6,R2);
                if (T2) ATaddvv(r6,T2);
                /* Check for change of reference momentum */
                if (scaling != 1.0) ATChangePRef(r6, 1.0/scaling);
   		    }
		}	
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double Length, Scaling;
        double *R1, *R2, *T1, *T2, *EApertures, *RApertures, *KickAngle;
        Length=atGetDouble(ElemData,"Length"); check_error();
        KickAngle=atGetDoubleArray(ElemData,"KickAngle"); check_error();
        /*optional fields*/
        Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();

        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->KickAngle = KickAngle;
        /*optional fields*/
        Elem->Scaling=Scaling;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
    }
	CorrectorPass(r_in, Elem->KickAngle[0], Elem->KickAngle[1], Elem->Length,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2,
            Elem->RApertures, Elem->EApertures,
	        Elem->Scaling, num_particles);
    return Elem;
}

MODULE_DEF(CorrectorPass)        /* Dummy module initialisation */
#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#ifdef MATLAB_MEX_FILE
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Length, Scaling;
        double *R1, *R2, *T1, *T2, *EApertures, *RApertures, *KickAngle;
        if (mxGetM(prhs[1]) != 6) mexErrMsgTxt("Second argument must be a 6 x N matrix");

        Length=atGetDouble(ElemData,"Length"); check_error();
        KickAngle=atGetDoubleArray(ElemData,"KickAngle"); check_error();
        /*optional fields*/
        Scaling=atGetOptionalDouble(ElemData,"FieldScaling",1.0); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();

        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        CorrectorPass(r_in, KickAngle[0], KickAngle[1], Length,
            T1, T2, R1, R2, RApertures, EApertures,
            Scaling, num_particles);
    } else if (nrhs == 0) {
        /* list of required fields */
	    plhs[0] = mxCreateCellMatrix(2,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("KickAngle"));
        if (nlhs>1) {    /* list of optional fields */
	        plhs[1] = mxCreateCellMatrix(8,0);
            mxSetCell(plhs[1],1,mxCreateString("T1"));
            mxSetCell(plhs[1],2,mxCreateString("T2"));
            mxSetCell(plhs[1],3,mxCreateString("R1"));
            mxSetCell(plhs[1],4,mxCreateString("R2"));
            mxSetCell(plhs[1],5,mxCreateString("RApertures"));
            mxSetCell(plhs[1],6,mxCreateString("EApertures"));
            mxSetCell(plhs[1],7,mxCreateString("FieldScaling"));
	    }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
