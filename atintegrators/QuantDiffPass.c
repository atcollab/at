#include "atelem.c"
#include <time.h>
#if !(defined PCWIN || defined PCWIN32 || defined PCWIN64 || defined _WIN32)
#include <sys/time.h>
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct elem 
{
  double* Lmatp;
  /* optional */
  int Seed;
};

double drand(void);   /* uniform distribution, (0..1] */
double random_normal(void);  /* normal distribution, centered on 0, std dev 1 */
double generateGaussian(double mean,double stdDev);

void QuantDiffPass(double *r_in, double* Lmatp , int Seed, int nturn, int num_particles)
/* Lmatp 6x6 matrix
 * r_in - 6-by-N matrix of initial conditions reshaped into
 * 1-d array of 6*N elements
 */
{
  int c;
  static int initSeed = 1;	/* 	If this variable is 1, then I initialize the seed
				 * 	to the clock and I change the variable to 0*/
  
  if(Seed)
  {
      if(nturn==0)
      {
          srand(Seed);
      }
  }
  else if(initSeed)
  {
#if !(defined PCWIN || defined PCWIN32 || defined PCWIN64 || defined _WIN32)
      {
      struct timeval time;
      gettimeofday(&time,NULL);
      srand((time.tv_sec * 1000000) + (time.tv_usec));
      }
#endif
      initSeed = 0;
  }

  #pragma omp parallel for if (num_particles > OMP_PARTICLE_THRESHOLD) default(none) \
  shared(r_in,num_particles,Lmatp) \
  private(c)
  for (c = 0; c<num_particles; c++) {
      /*Loop over particles  */
      int i, j;
      double randnorm[6];
      double diffusion[6];
      double *r6 = r_in+c*6;
      for (i=0;i<6;i++) {
          diffusion[i]=0.0;
          /*randnorm[i]=random_normal();*/
          randnorm[i]= generateGaussian(0.0,1.0);
      }
      
      for (i=0;i<6;i++) {
          for (j=0;j<=i;j++) {
              diffusion[i]+=randnorm[j]*Lmatp[i+6*j];
          }
      }
      if (!atIsNaN(r6[0])) {
          r6[0] += diffusion[0];
          r6[1] += diffusion[1];
          r6[2] += diffusion[2];
          r6[3] += diffusion[3];
          r6[4] += diffusion[4];
          r6[5] += diffusion[5];
      }
  }
}


double drand(void)   /* uniform distribution, (0..1] */
{
    return (rand()+1.0)/(RAND_MAX+1.0);
}
double random_normal(void)  /* normal distribution, centered on 0, std dev 1 */
{
    return sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}

double generateGaussian(double mean,double stdDev)
{
	static bool hasSpare = false;
	static double spare;
	static double u, v, s;

	if (hasSpare) {
		hasSpare = false;
		return mean + stdDev * spare;
	}

	hasSpare = true;
	do {
		u = (rand() / ((double) RAND_MAX)) * 2.0 - 1.0;
		v = (rand() / ((double) RAND_MAX)) * 2.0 - 1.0;
		s = u * u + v * v;
	}
	while( (s >= 1.0) || (s == 0.0) );
	s = sqrt(-2.0 * log(s) / s);
	spare = v * s;
	return mean + stdDev * u * s;
}


#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
        double *r_in, int num_particles, struct parameters *Param)
{
    int nturn=Param->nturn;
    if (!Elem) {
        double *Lmatp;
        int Seed;
        Lmatp=atGetDoubleArray(ElemData,"Lmatp"); check_error();
        Seed=atGetOptionalLong(ElemData,"Seed",0); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        Elem->Lmatp=Lmatp;
        Elem->Seed=Seed;   
    }
    QuantDiffPass(r_in, Elem->Lmatp, Elem->Seed, nturn, num_particles);
    return Elem;
}

MODULE_DEF(QuantDiffPass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#if defined(MATLAB_MEX_FILE)

void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        int Seed;
        double *Lmatp;
        Lmatp=atGetDoubleArray(ElemData,"Lmatp"); check_error();
        Seed=atGetOptionalLong(ElemData,"Seed",0); check_error();
        if (mxGetM(prhs[1]) != 6) mexErrMsgIdAndTxt("AT:WrongArg","Second argument must be a 6 x N matrix");
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        QuantDiffPass(r_in, Lmatp, Seed, 0, num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(1,1);
        mxSetCell(plhs[0],0,mxCreateString("Lmatp"));
        if (nlhs>1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(1,1);
            mxSetCell(plhs[1],0,mxCreateString("Seed"));
        }
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
