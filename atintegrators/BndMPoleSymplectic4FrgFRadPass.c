#include "at.h"
#include "atlalib.c"
#include "atphyslib.c"		/* edge, edge_fringe */
#include "driftkickrad.c"	/* drift6.c, bndthinkickrad.c */
#include "quadfringe.c"		/* QuadFringePassP, QuadFringePassN */

#define DRIFT1    0.6756035959798286638
#define DRIFT2   -0.1756035959798286639
#define KICK1     1.351207191959657328
#define KICK2    -1.702414383919314656

void BndMPoleSymplectic4FrgFRadPass(double *r, double le, double irho, double *A, double *B,
        int max_order, int num_int_steps,
        double entrance_angle, 	double exit_angle,
        double fint1, double fint2, double gap,
        double *T1, double *T2,
        double *R1, double *R2,
        double *RApertures, double *EApertures,
        double E0, int num_particles)
{
    double *r6;
    int c, m;
    bool useT1 = (T1 != NULL);
    bool useT2 = (T2 != NULL);
    bool useR1 = (R1 != NULL);
    bool useR2 = (R2 != NULL);
    double SL = le/num_int_steps;
    double L1 = SL*DRIFT1;
    double L2 = SL*DRIFT2;
    double K1 = SL*KICK1;
    double K2 = SL*KICK2;
    bool useFringe1 = ((fint1!=0) && (gap!=0));
    bool useFringe2 = ((fint2!=0) && (gap!=0));
    
    for (c = 0; c<num_particles; c++) {	/* Loop over particles  */
        r6 = r+c*6;
        /*  misalignment at entrance  */
        if(useT1) ATaddvv(r6,T1);
        if(useR1) ATmultmv(r6,R1);
        /* Check physical apertures at the entrance of the magnet */
        if (RApertures) checkiflostRectangularAp(r6,RApertures);
        if (EApertures) checkiflostEllipticalAp(r6,EApertures);
        /* edge focus */
        if(useFringe1)
            edge_fringe(r6, irho, entrance_angle,fint1,gap);
        else
            edge(r6, irho, entrance_angle);
        /* quadrupole gradient fringe */
        QuadFringePassP(r6,B[1]);
        /* integrator  */
        for(m=0; m < num_int_steps; m++){ /* Loop over slices */
            drift6(r6,L1);
            bndthinkickrad(r6, A, B, K1, irho, E0, max_order);
            drift6(r6,L2);
            bndthinkickrad(r6, A, B, K2, irho, E0, max_order);
            drift6(r6,L2);
            bndthinkickrad(r6, A, B,  K1, irho, E0, max_order);
            drift6(r6,L1);
        }
        /* quadrupole gradient fringe */
        QuadFringePassN(r6,B[1]);
        /* edge focus */
        if(useFringe2)
            edge_fringe(r6, irho, exit_angle,fint2,gap);
        else
            edge(r6, irho, exit_angle);
        /* Check physical apertures at the exit of the magnet */
        if (RApertures) checkiflostRectangularAp(r6,RApertures);
        if (EApertures) checkiflostEllipticalAp(r6,EApertures);
        /* Misalignment at exit */
        if(useR2) ATmultmv(r6,R2);
        if(useT2) ATaddvv(r6,T2);
    }
}


#ifdef MATLAB_MEX_FILE
#include "elempass.h"
#include "mxutils.c"

ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
        double *r_in, int num_particles, int mode)
        
#define NUM_FIELDS_2_REMEMBER 18
        
        
{	double *A , *B;
    double  *pr1, *pr2, *pt1, *pt2, *RApertures, *EApertures, fint1, fint2, gap;
    double entrance_angle, exit_angle;
    double E0;		/* Design energy [eV] */
    int max_order, num_int_steps;
    double le,ba,irho;
    
    switch(mode) {
        case MAKE_LOCAL_COPY: 	/* Find field numbers first
         * Save a list of field number in an array
         * and make returnptr point to that array
         */
            /* Allocate memory for integer array of
             * field numbers for faster future reference
             */
            
            FieldNumbers = (int*)mxCalloc(NUM_FIELDS_2_REMEMBER,sizeof(int));
            
            /* Populate */
            
            FieldNumbers[0] = GetRequiredFieldNumber(ElemData, "PolynomA");
            FieldNumbers[1] = GetRequiredFieldNumber(ElemData, "PolynomB");
            FieldNumbers[2] = GetRequiredFieldNumber(ElemData, "MaxOrder");
            FieldNumbers[3] = GetRequiredFieldNumber(ElemData, "NumIntSteps");
            FieldNumbers[4] = GetRequiredFieldNumber(ElemData, "Length");
            FieldNumbers[5] = GetRequiredFieldNumber(ElemData, "Energy");
            FieldNumbers[6] = GetRequiredFieldNumber(ElemData, "BendingAngle");
            FieldNumbers[7] = GetRequiredFieldNumber(ElemData, "EntranceAngle");
            FieldNumbers[8] = GetRequiredFieldNumber(ElemData, "ExitAngle");
            
            FieldNumbers[9] = mxGetFieldNumber(ElemData,"FringeInt1");
            FieldNumbers[10] = mxGetFieldNumber(ElemData,"FringeInt2");
            FieldNumbers[11] = mxGetFieldNumber(ElemData,"FullGap");
            FieldNumbers[12] = mxGetFieldNumber(ElemData,"R1");
            FieldNumbers[13] = mxGetFieldNumber(ElemData,"R2");
            FieldNumbers[14] = mxGetFieldNumber(ElemData,"T1");
            FieldNumbers[15] = mxGetFieldNumber(ElemData,"T2");
            FieldNumbers[16] = mxGetFieldNumber(ElemData,"RApertures");
            FieldNumbers[17] = mxGetFieldNumber(ElemData,"EApertures");
            /* Fall through next section... */
            
        case	USE_LOCAL_COPY:	/* Get fields from MATLAB using field numbers
         * The second argument ponter to the array of field
         * numbers is previously created with
         * QuadLinPass( ..., MAKE_LOCAL_COPY)
         */
            
            A = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[0]));
            B = mxGetPr(mxGetFieldByNumber(ElemData,0,FieldNumbers[1]));
            max_order = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[2]));
            num_int_steps = (int)mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[3]));
            le = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[4]));
            E0 = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[5]));
            ba = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[6]));
            entrance_angle = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[7]));
            exit_angle = mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[8]));
            
            /* Optional fields */
            
            fint1=(FieldNumbers[9] >= 0) ? mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[9])) : 0;
            fint2=(FieldNumbers[10] >= 0) ? mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[10])) : 0;
            gap = (FieldNumbers[11] >= 0) ? mxGetScalar(mxGetFieldByNumber(ElemData,0,FieldNumbers[11])) : 0;
            pr1 = (FieldNumbers[12] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[12])) : NULL;
            pr2 = (FieldNumbers[13] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[13])) : NULL;
            pt1 = (FieldNumbers[14] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[14])) : NULL;
            pt2 = (FieldNumbers[15] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[15])) : NULL;
            RApertures = (FieldNumbers[16] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[16])) : NULL;
            EApertures = (FieldNumbers[17] >= 0) ? mxGetPr(mxGetFieldByNumber(ElemData, 0, FieldNumbers[17])) : NULL;
            break;
        default:
            mexErrMsgTxt("No match for calling mode in function BndMPoleSymplectic4FrgFRadPass\n");
    }
    
    irho = ba/le;
    
    BndMPoleSymplectic4FrgFRadPass(r_in, le, irho, A, B, max_order, num_int_steps,
            entrance_angle, exit_angle, fint1, fint2, gap, pt1, pt2, pr1, pr2, RApertures, EApertures, E0, num_particles);
    
    return FieldNumbers;
}


void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        double *pr1, *pr2, *pt1, *pt2, *RApertures, *EApertures;
        double fint1, fint2, gap, irho;
        mxArray *tmpmxptr;
        
        double *A = mxGetPr(GetRequiredField(prhs[0], "PolynomA"));
        double *B = mxGetPr(GetRequiredField(prhs[0], "PolynomB"));
        int max_order = (int) mxGetScalar(GetRequiredField(prhs[0], "MaxOrder"));
        int num_int_steps = (int) mxGetScalar(GetRequiredField(prhs[0], "NumIntSteps"));
        double le = mxGetScalar(GetRequiredField(prhs[0], "Length"));
        double E0 = mxGetScalar(GetRequiredField(prhs[0], "Energy"));
        double ba = mxGetScalar(GetRequiredField(prhs[0], "BendingAngle"));
        double entrance_angle = mxGetScalar(GetRequiredField(prhs[0], "EntranceAngle"));
        double exit_angle = mxGetScalar(GetRequiredField(prhs[0], "ExitAngle"));
        int num_particles = mxGetN(prhs[1]);
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        
        /* Optional arguments */
        tmpmxptr = mxGetField(prhs[0],0,"FringeInt1");
        fint1 = tmpmxptr ? mxGetScalar(tmpmxptr) : 0;
        
        tmpmxptr = mxGetField(prhs[0],0,"FringeInt2");
        fint2 = tmpmxptr ? mxGetScalar(tmpmxptr) : 0;
        
        tmpmxptr = mxGetField(prhs[0],0,"FullGap");
        gap = tmpmxptr ? mxGetScalar(tmpmxptr) : 0;
        
        tmpmxptr = mxGetField(prhs[0],0,"R1");
        pr1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"R2");
        pr2 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"T1");
        pt1 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"T2");
        pt2 = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"RApertures");
        RApertures = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        tmpmxptr = mxGetField(prhs[0],0,"EApertures");
        EApertures = tmpmxptr ? mxGetPr(tmpmxptr) : NULL;
        
        irho = ba/le;
        
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        BndMPoleSymplectic4FrgFRadPass(r_in, le, irho, A, B, max_order, num_int_steps,
                entrance_angle, exit_angle, fint1, fint2, gap, pt1, pt2, pr1, pr2, RApertures, EApertures, E0, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(9,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("BendingAngle"));
        mxSetCell(plhs[0],2,mxCreateString("EntranceAngle"));
        mxSetCell(plhs[0],3,mxCreateString("ExitAngle"));
        mxSetCell(plhs[0],4,mxCreateString("PolynomA"));
        mxSetCell(plhs[0],5,mxCreateString("PolynomB"));
        mxSetCell(plhs[0],6,mxCreateString("MaxOrder"));
        mxSetCell(plhs[0],7,mxCreateString("NumIntSteps"));
        mxSetCell(plhs[0],8,mxCreateString("Energy"));
        
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(9,1);
            mxSetCell(plhs[1],0,mxCreateString("FullGap"));
            mxSetCell(plhs[1],1,mxCreateString("FringeInt1"));
            mxSetCell(plhs[1],2,mxCreateString("FringeInt2"));
            mxSetCell(plhs[1],3,mxCreateString("T1"));
            mxSetCell(plhs[1],4,mxCreateString("T2"));
            mxSetCell(plhs[1],5,mxCreateString("R1"));
            mxSetCell(plhs[1],6,mxCreateString("R2"));
            mxSetCell(plhs[1],7,mxCreateString("RApertures"));
            mxSetCell(plhs[1],8,mxCreateString("EApertures"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif
