/* CorrectorPass.c 
   Accelerator Toolbox 
   Revision 10/11/01
   A.Terebilo terebilo@ssrl.slac.stanford.edu
   CorrectorPass expects an element to have a fields:
   Length, KickAngle
*/

#include "atelem.c"

struct elem
{
    double Length;
    double *KickAngle;
};

void CorrectorPass(double *r_in, double xkick, double ykick, double len,  int num_particles)
/* xkick, ykick - horizontal and vertical kicks in radiand 
   r - 6-by-N matrix of initial conditions reshaped into 
   1-d array of 6*N elements 
*/
{	int c, c6;
    double NormL, p_norm;
	if (len==0)
	    for(c = 0;c<num_particles;c++) {
	        c6 = c*6;
		    if(!atIsNaN(r_in[c6])) {
		        r_in[c6+1] += xkick;
   		        r_in[c6+3] += ykick; 		    
   		    }
		}	
	else
        #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD*10) default(shared) shared(r_in,num_particles) private(c,c6)
        for(c = 0;c<num_particles;c++) {
            c6 = c*6;
		    if(!atIsNaN(r_in[c6])) {
		        p_norm = 1/(1+r_in[c6+4]);
			    NormL  = len*p_norm;
	            r_in[c6+5] += NormL*p_norm*(xkick*xkick/3 + ykick*ykick/3 +
   		                r_in[c6+1]*r_in[c6+1] + r_in[c6+3]*r_in[c6+3] + 
   		                r_in[c6+1]*xkick + r_in[c6+3]*ykick)/2;

			    r_in[c6]   += NormL*(r_in[c6+1]+xkick/2);
		        r_in[c6+1] += xkick;
		        r_in[c6+2] += NormL*(r_in[c6+3]+ykick/2);
   		        r_in[c6+3] += ykick;
   		    }
		}	
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    if (!Elem) {
        double Length;
        double *KickAngle;
        Length=atGetDouble(ElemData,"Length"); check_error();
        KickAngle=atGetDoubleArray(ElemData, "KickAngle"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Length=Length;
        Elem->KickAngle = KickAngle;
    }
	CorrectorPass(r_in, Elem->KickAngle[0], Elem->KickAngle[1], Elem->Length, num_particles);
    return Elem;
}

MODULE_DEF(CorrectorPass)        /* Dummy module initialisation */
#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#ifdef MATLAB_MEX_FILE
void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        double Length;
        double *KickAngle;
        Length=atGetDouble(ElemData,"Length"); check_error();
        KickAngle=atGetDoubleArray(ElemData, "KickAngle"); check_error();
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        CorrectorPass(r_in, KickAngle[0], KickAngle[1], Length, num_particles);
      }
    else {
	    plhs[0] = mxCreateCellMatrix(2,1);
	    mxSetCell(plhs[0],0,mxCreateString("Length"));
	    mxSetCell(plhs[0],1,mxCreateString("KickAngle"));
	    if(nlhs>1) /* Required and optional fields */ 
	    {   plhs[1] = mxCreateCellMatrix(0,0); /* No optional fields */
	    }
    }
}
#endif /*MATLAB_MEX_FILE*/
