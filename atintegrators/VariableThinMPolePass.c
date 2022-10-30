/* VariableThinMPolePass
   Accelerator Toolbox 
   S.White simon.white@esrf.fr
*/

#include "atelem.c"
#include "atlalib.c"
#include "driftkick.c"

#define TWOPI  6.28318530717959


struct elemab
{
    double *Amplitude;
    double Frequency;
    double Phase;
    int NSamples;
    double *Func;
};

struct elem 
{
    double *PolynomA;
    double *PolynomB;
    struct elemab *ElemA;
    struct elemab *ElemB;
    int Seed;
    int Mode;
    int MaxOrder;
    double *Ramps;
};


double randn(int seed)
{
  double u, v, w, mult;
  static int initseed=1;

  if(initseed){
      srand(seed);
      initseed=0;
  }

  do{
      u = -1 + ((double) rand () / RAND_MAX) * 2;
      v = -1 + ((double) rand () / RAND_MAX) * 2;
      w = u*u + v*v;
  }
  while (w >= 1 || w == 0);
  mult=sqrt ((-2 * log (w)) / w);
  return u*mult;
}

double get_amp(double amp, double *ramps, double t)
{
    double ampt=amp;
    if(ramps){
        if(t<=ramps[0]){
            ampt=0.0;
        }else if(t<=ramps[1]){
            ampt=amp*(t-ramps[1])/(ramps[2]-ramps[1]);
        }else if(t<=ramps[2]){
            ampt=amp;
        }else if(t<=ramps[3]){
            ampt=amp-amp*(t-ramps[2])/(ramps[3]-ramps[2]);
        }else{
            ampt=0.0;
        }
    }
    return ampt;
}

double get_pol(struct elemab *elem, double *ramps, int mode,
               double t, int turn, int seed, int order)
{
    int idx;
    double ampt, freq, ph, val;
    double *func;
    double *amp = elem->Amplitude;
    if(!amp){
        return 0.0;
    }
    ampt=get_amp(amp[order],ramps,t);
    switch(mode){
    case 0:
        freq = elem->Frequency;
        ph = elem->Phase;
        ampt *= sin(TWOPI*freq*t+ph);
        return ampt;
    case 1:
        val = randn(seed);
        ampt *= val;
        return ampt;
    case 2:
        func = elem->Func;
        idx = turn%elem->NSamples;
        ampt *= func[idx];
        return ampt;
    }
}


void VariableThinMPolePass(double *r, struct elem *Elem, double t0, int turn, int num_particles)
{

    int i, c;
    double *r6;
    double t = t0*turn;

    int maxorder = Elem->MaxOrder;
    double *pola = Elem->PolynomA;
    double *polb = Elem->PolynomB;
    int seed = Elem->Seed;
    int mode = Elem->Mode;
    struct elemab *ElemA = Elem->ElemA;
    struct elemab *ElemB = Elem->ElemB;
    double *ramps = Elem->Ramps;

    for(i=0;i<maxorder;i++){
        pola[i]=get_pol(ElemA, ramps, mode, t, turn, seed, i);
        polb[i]=get_pol(ElemB, ramps, mode, t, turn, seed, i);
    };

    for (c = 0;c<num_particles;c++){
        r6 = r+c*6;
        if (!atIsNaN(r6[0])) {
            strthinkick(r6, pola, polb, 1.0, maxorder);
        }
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem *trackFunction(const atElem *ElemData,struct elem *Elem,
			      double *r_in, int num_particles, struct parameters *Param)
{   
    if (!Elem) {
        int MaxOrder, Mode, Seed, NSamplesA, NSamplesB;
        double *PolynomA, *PolynomB, *AmplitudeA, *AmplitudeB;
        double *Ramps, *FuncA, *FuncB;
        double FrequencyA, FrequencyB;
        double PhaseA, PhaseB;
        struct elemab *ElemA, *ElemB;
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        Mode=atGetLong(ElemData,"Mode"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();
        AmplitudeA=atGetOptionalDoubleArray(ElemData,"AmplitudeA"); check_error();
        AmplitudeB=atGetOptionalDoubleArray(ElemData,"AmplitudeB"); check_error();
        FrequencyA=atGetOptionalDouble(ElemData,"FrequencyA", 0); check_error();
        FrequencyB=atGetOptionalDouble(ElemData,"FrequencyB", 0); check_error();
        PhaseA=atGetOptionalDouble(ElemData,"PhaseA", 0); check_error();
        PhaseB=atGetOptionalDouble(ElemData,"PhaseB", 0); check_error();
        Ramps=atGetOptionalDoubleArray(ElemData, "Ramps"); check_error();
        Seed=atGetOptionalLong(ElemData, "Seed", 0);
        NSamplesA=atGetOptionalLong(ElemData, "NSamplesA", 0);
        NSamplesB=atGetOptionalLong(ElemData, "NSamplesA", 0);
        FuncA=atGetOptionalDoubleArray(ElemData,"FuncA"); check_error();
        FuncB=atGetOptionalDoubleArray(ElemData,"FuncB"); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        ElemA = (struct elemab *)atMalloc(sizeof(struct elemab));
        ElemB = (struct elemab *)atMalloc(sizeof(struct elemab));
        Elem->PolynomA=PolynomA;
        Elem->PolynomB=PolynomB;
        Elem->Ramps=Ramps;
        Elem->Seed=Seed;
        Elem->Mode=Mode;
        Elem->MaxOrder=MaxOrder;
        ElemA->Amplitude=AmplitudeA;
        ElemB->Amplitude=AmplitudeB;
        ElemA->Frequency=FrequencyA;
        ElemB->Frequency=FrequencyB;
        ElemA->Phase=PhaseA;
        ElemB->Phase=PhaseB;
        ElemA->NSamples=NSamplesA;
        ElemB->NSamples=NSamplesB;
        ElemA->Func=FuncA;
        ElemB->Func=FuncB;
        Elem->ElemA = ElemA;
        Elem->ElemB = ElemB;
    }
    double t0=Param->T0;
    int turn=Param->nturn;
    VariableThinMPolePass(r_in, Elem, t0, turn, num_particles);
    return Elem;
}

MODULE_DEF(VariableThinMPolePass)        /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#ifdef MATLAB_MEX_FILE
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 2 ) {
        double *r_in;
        const mxArray *ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        int MaxOrder;
        double *PolynomA, *PolynomB;
        double *AmplitudeA, *AmplitudeB;
        double frequency;
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        AmplitudeA=atGetDoubleArray(ElemData,"AmplitudeA"); check_error();
        AmplitudeB=atGetDoubleArray(ElemData,"AmplitudeB"); check_error();
        PolynomA=atGetDoubleArray(ElemData,"PolynomA"); check_error();
        PolynomB=atGetDoubleArray(ElemData,"PolynomB"); check_error();        
        frequency=atGetDouble(ElemData,"Frequency"); check_error();
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        VariableThinMPolePass(r_in, PolynomA, PolynomB, MaxOrder,
                              num_particles);
    }
    else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(6,1);
        mxSetCell(plhs[0],0,mxCreateString("MaxOrder"));
        mxSetCell(plhs[0],1,mxCreateString("AmplitudeA"));
        mxSetCell(plhs[0],2,mxCreateString("AmplitudeB"));
        mxSetCell(plhs[0],3,mxCreateString("PolynomA"));
        mxSetCell(plhs[0],4,mxCreateString("PolynomB"));
        mxSetCell(plhs[0],5,mxCreateString("Frequency"));
    }
    else {
        mexErrMsgIdAndTxt("AT:WrongArg","Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
