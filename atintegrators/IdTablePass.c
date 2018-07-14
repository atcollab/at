/* IdTablePass.c
 * Accelerator Toolbox
 * Created: 13/11/08
 * Z.Mart?? zeus@cells.es
 *
 * Based in the matlab routine:
 * WigTablePass.m - The tracking table is described in
 * P. Elleaume, "A new approach to the electron beam dynamics in undulators
 * and wigglers", EPAC92.
 *
 */

#include "atelem.c"
#include "atlalib.c"
#include "interpolate.c"

struct elem {
    double Length;
    double *xkick;
    double *ykick;
    double *xtable;
    double *ytable;
    int n_map;
    int m_map;
    int Nslice;
    /* Optional fields */
    double *xkick1;
    double *ykick1;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
};

double *GLOBAL_x, *GLOBAL_y;
int GLOBAL_m,GLOBAL_n;

static double get_kick(double *r6, double *ktable)
{
    double f;
    /*cubic interpolation*/
    /*splin2(GLOBAL_y,GLOBAL_x,GLOBAL_xkick,GLOBAL_xkick2,GLOBAL_n,GLOBAL_m,y,x,&f);*/
    
    /*biliniar interpolation*/
#ifdef MATLAB_MEX_FILE
    /* Transpose coordinates because kick-table is FORTRAN-ordered */
    linint(GLOBAL_y, GLOBAL_x, ktable, GLOBAL_m, GLOBAL_n, r6[2], r6[0], &f);
#else
    /* Assume kick-table is C-ordered */
    linint(GLOBAL_x, GLOBAL_y, ktable, GLOBAL_m, GLOBAL_n, r6[0], r6[2], &f);
#endif
    return f;
 }

void IdKickMapModelPass(double *r, double le, double *xkick1, double *ykick1,
        double *xkick, double *ykick, double *x, double *y,int n,int m, int Nslice,
        double *T1, double *T2, double *R1, double *R2, int num_particles)
{
    double *r6, deltaxp, deltayp, limitsptr[4];
    int c;
    double L1 = le/(2*Nslice);
    
    /*Act as AperturePass*/
    limitsptr[0]=x[0];
    limitsptr[1]=x[n-1];
    limitsptr[2]=y[0];
    limitsptr[3]=y[m-1];
    
    /*globalize*/
    
    /* For cubic interpolation only*/
    
    GLOBAL_x=x;
    GLOBAL_y=y;
    GLOBAL_m=m; /* y used as colums*/
    GLOBAL_n=n; /* x used as rows*/
    
    for (c=0; c<num_particles; c++) {
        r6 = r+c*6;
        if (!atIsNaN(r6[0])) {
            /* Misalignment at entrance */
            if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
              checkiflostRectangularAp(r6, limitsptr);
            /*Tracking in the main body*/
            for (m=0; m<Nslice; m++) { /* Loop over slices*/
                ATdrift6(r6,L1);
                if (!atIsNaN(r6[0])&&!atIsNaN(r6[2])) {
                    /*The kick from IDs varies quadratically, not linearly, with energy.   */
                    deltaxp = get_kick(r6, xkick)/(1.0+r6[4]);
                    deltayp = get_kick(r6, ykick)/(1.0+r6[4]);
                    if (xkick1)  deltaxp += get_kick(r6, xkick1);
                    if (ykick1)  deltayp += get_kick(r6, ykick1);
                    r6[1] = r6[1] + deltaxp / Nslice;
                    r6[3] = r6[3] + deltayp / Nslice;
                }
                ATdrift6(r6,L1);
            }
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
        int Nslice, n_map,m_map;
        double Length, *xkick, *ykick, *xtable, *ytable;
        double *xkick1, *ykick1;
        double *R1, *R2, *T1, *T2;
        Length=atGetDouble(ElemData,"Length"); check_error();
        xkick=atGetDoubleArray(ElemData,"xkick"); check_error();
        ykick=atGetDoubleArray(ElemData,"ykick"); check_error();
        xtable=atGetDoubleArraySz(ElemData,"xtable", &m_map, &n_map); check_error();
        ytable=atGetDoubleArray(ElemData,"ytable"); check_error();
        Nslice=atGetLong(ElemData,"Nslice"); check_error();
        /*optional fields*/
        xkick1=atGetOptionalDoubleArray(ElemData,"xkick1"); check_error();
        ykick1=atGetOptionalDoubleArray(ElemData,"ykick1"); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->xkick=xkick;
        Elem->ykick=ykick;
        Elem->xtable=xtable;
        Elem->ytable=ytable;
        Elem->n_map=n_map;
        Elem->m_map=m_map;
        Elem->Nslice=Nslice;
        Elem->xkick1=xkick1;
        Elem->ykick1=ykick1;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
    }
    IdKickMapModelPass(r_in, Elem->Length, Elem->xkick1, Elem->ykick1,
            Elem->xkick, Elem->ykick, Elem->xtable, Elem->ytable, Elem->n_map, Elem->m_map,Elem->Nslice,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2, num_particles);
    return Elem;
}

MODULE_DEF(IdTablePass)        /* Dummy module initialisation */
#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#ifdef MATLAB_MEX_FILE
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        int Nslice, n_map,m_map;
        double Length, *xkick, *ykick, *xtable, *ytable;
        double *xkick1, *ykick1;
        double *R1, *R2, *T1, *T2;
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        Length=atGetDouble(ElemData,"Length"); check_error();
        xkick=atGetDoubleArray(ElemData,"xkick"); check_error();
        ykick=atGetDoubleArray(ElemData,"ykick"); check_error();
        xtable=atGetDoubleArraySz(ElemData,"xtable", &m_map, &n_map); check_error();
        ytable=atGetDoubleArray(ElemData,"ytable"); check_error();
        Nslice=atGetLong(ElemData,"Nslice"); check_error();
        /*optional fields*/
        xkick1=atGetOptionalDoubleArray(ElemData,"xkick1"); check_error();
        ykick1=atGetOptionalDoubleArray(ElemData,"ykick1"); check_error();
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetPr(plhs[0]);
        IdKickMapModelPass(r_in, Length, xkick1, ykick1, xkick, ykick, xtable, ytable,
                n_map, m_map, Nslice, T1, T2, R1, R2, num_particles);
    }
    else if (nrhs == 0) {
        /* return list of required fields */
        plhs[0] = mxCreateCellMatrix(6,1);
        mxSetCell(plhs[0],0,mxCreateString("Length"));
        mxSetCell(plhs[0],1,mxCreateString("xkick"));
        mxSetCell(plhs[0],2,mxCreateString("ykick"));
        mxSetCell(plhs[0],3,mxCreateString("xtable"));
        mxSetCell(plhs[0],4,mxCreateString("ytable"));
        mxSetCell(plhs[0],5,mxCreateString("Nslice"));
        
        if (nlhs > 1) {
            /* Required and optional fields */
            plhs[1] = mxCreateCellMatrix(6,1);
            mxSetCell(plhs[1],0,mxCreateString("xkick1"));
            mxSetCell(plhs[1],1,mxCreateString("ykick1"));
            mxSetCell(plhs[1],2,mxCreateString("T1"));
            mxSetCell(plhs[1],3,mxCreateString("T2"));
            mxSetCell(plhs[1],4,mxCreateString("R1"));
            mxSetCell(plhs[1],5,mxCreateString("R2"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif
