#include "atelem.c"
#include "atlalib.c"

struct elem {
    double dpp;
};

void EnergyLossRadPass(double *r_in, double dpp, int num_particles)
{
    for (int c = 0; c<num_particles; c++) { /*Loop over particles  */
        
        /* Get the pointer to the current particle's state vector. */
        double *r6 = r_in + 6*c;
        
        if (!atIsNaN(r6[0])) {
            double p_norm, ddp;
            double xpr, ypr;

            /* calculate angles from tranverse momenta 	*/
             p_norm = 1.0 + r6[4];
             xpr = r6[1]/p_norm;
             ypr = r6[3]/p_norm;

            /* momentum step */
             ddp = dpp * p_norm *p_norm;    /* dpp: nominal momentum step */
             r6[4] -= ddp;

             /* recalculate momenta from angles after losing energy for radiation 	*/
             p_norm = 1.0 + r6[4];
             r6[1] = xpr*p_norm;
             r6[3] = ypr*p_norm;
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
                double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double EnergyLoss=atGetDouble(ElemData,"EnergyLoss"); check_error();
        double Energy=atGetOptionalDouble(ElemData,"Energy",Param->energy); check_error();
        double dpp = EnergyLoss/Energy;

        if (Param->rest_energy != 0.0) {
            double gamma = Param->energy/Param->rest_energy;
            double beta2 = 1.0 - 1.0/gamma/gamma;
            dpp /= beta2;
        }
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->dpp = dpp;
    }
    EnergyLossRadPass(r_in, Elem->dpp, num_particles);
    return Elem;
}

MODULE_DEF(EnergyLossRadPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     /* Check if the number of input arguments is correct. */
    if (nrhs >= 2) {
        /* Get the input arguments. */
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double EnergyLoss=atGetDouble(ElemData,"EnergyLoss"); check_error();
        double Energy=atGetDouble(ElemData,"Energy"); check_error();
        double dpp = EnergyLoss/Energy;

        /* Check if the second argument is a 6 x N matrix. */
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        
         /* Call the EnergyLossRadPass function to apply the energy loss. */
        EnergyLossRadPass(r_in, dpp , num_particles);

    } else if (nrhs == 0) {
        /* list of required fields */
        /* Create a cell array to store the required fields. */
        plhs[0] = mxCreateCellMatrix(2,1);
        /* Set the elements of the cell array. */
        mxSetCell(plhs[0],0,mxCreateString("EnergyLoss"));
        mxSetCell(plhs[0],1,mxCreateString("Energy"));
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(0,1);
        }
    } else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /*defined(MATLAB_MEX_FILE)*/
