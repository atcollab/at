/* IdTablePass.c
 * Accelerator Toolbox
 * Created: 13/11/08
 * Zeus Marti
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
    double *x_map;
    double *y_map;
    int nx_map;
    int ny_map;
    int Nslice;
    /* Optional fields */
    double *xkick1;
    double *ykick1;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
};

double *GLOBAL_x_map, *GLOBAL_y_map;
int GLOBAL_nx_map,GLOBAL_ny_map;

static double get_kick(double *r6, double *ktable)
{
    double f;
    /*cubic interpolation*/
    /*splin2(GLOBAL_y_map,GLOBAL_x_map,GLOBAL_xkick,GLOBAL_xkick2,GLOBAL_n,GLOBAL_nx,y,x,&f);*/
    
    /*biliniar interpolation*/
    /* Transpose coordinates because the kick-table is FORTRAN-ordered */
    linint(GLOBAL_y_map, GLOBAL_x_map, ktable, GLOBAL_ny_map, GLOBAL_nx_map, r6[2], r6[0], &f);
    return f;
 }

void IdKickMapModelPass(double *r, double le, double *xkick1, double *ykick1,
        double *xkick, double *ykick, double *x_map, double *y_map,int nx_map,int ny_map, int Nslice,
        double *T1, double *T2, double *R1, double *R2, int num_particles)
{
    double *r6, deltaxp, deltayp, limitsptr[4];
    int c, ns;
    double L1 = le/(2*Nslice);
    
    /*Act as AperturePass*/
    limitsptr[0]=x_map[0];
    limitsptr[1]=x_map[nx_map-1];
    limitsptr[2]=y_map[0];
    limitsptr[3]=y_map[ny_map-1];
    
    /*globalize*/    
    GLOBAL_x_map=x_map;
    GLOBAL_y_map=y_map;
    GLOBAL_nx_map=nx_map;
    GLOBAL_ny_map=ny_map;
    
    for (c=0; c<num_particles; c++) {
        r6 = r+c*6;
        if (!atIsNaN(r6[0])) {
            /* Misalignment at entrance */
            if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
            checkiflostRectangularAp(r6, limitsptr);
            /*Tracking in the main body*/
            for (ns=0; ns<Nslice; ns++) { /* Loop over slices*/
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
        int Nslice, ny_map,nx_map;
        double Length, *xkick, *ykick, *x_map, *y_map;
        double *xkick1, *ykick1;
        double *R1, *R2, *T1, *T2;
        Length=atGetDouble(ElemData,"Length"); check_error();
        xkick=atGetDoubleArraySz(ElemData,"xkick", &ny_map, &nx_map); check_error();
        /* the third input of atGetDoubleArraySz is a pointer for the 
         * number of rows in the 2-D array, the fourth input is a pointer 
         * for the number of columns
         */
        ykick=atGetDoubleArray(ElemData,"ykick"); check_error();
        x_map=atGetDoubleArray(ElemData,"xtable"); check_error();
        y_map=atGetDoubleArray(ElemData,"ytable"); check_error();
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
        Elem->x_map=x_map;
        Elem->y_map=y_map;
        Elem->nx_map=nx_map;
        Elem->ny_map=ny_map;
        Elem->Nslice=Nslice;
        Elem->xkick1=xkick1;
        Elem->ykick1=ykick1;
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
    }
    IdKickMapModelPass(r_in, Elem->Length, Elem->xkick1, Elem->ykick1,
            Elem->xkick, Elem->ykick, Elem->x_map, Elem->y_map, Elem->nx_map, Elem->ny_map,Elem->Nslice,
            Elem->T1, Elem->T2, Elem->R1, Elem->R2, num_particles);
    return Elem;
}

MODULE_DEF(IdTablePass)        /* Dummy module initialisation */
#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#ifdef MATLAB_MEX_FILE
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs >= 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        int Nslice, ny_map, nx_map;
        double Length, *xkick, *ykick, *x_map, *y_map;
        double *xkick1, *ykick1;
        double *R1, *R2, *T1, *T2;
        Length=atGetDouble(ElemData,"Length"); check_error();
        xkick=atGetDoubleArraySz(ElemData,"xkick", &ny_map, &nx_map); check_error();
        /* the third input of atGetDoubleArraySz is a pointer for the 
         * number of rows in the 2-D array, the fourth input is a pointer 
         * for the number of columns
         */
        ykick=atGetDoubleArray(ElemData,"ykick"); check_error();
        x_map=atGetDoubleArray(ElemData,"xtable"); check_error();
        y_map=atGetDoubleArray(ElemData,"ytable"); check_error();
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
        r_in = mxGetDoubles(plhs[0]);
        IdKickMapModelPass(r_in, Length, xkick1, ykick1, xkick, ykick, x_map, y_map,
                nx_map, ny_map, Nslice, T1, T2, R1, R2, num_particles);
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
