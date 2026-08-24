/* VariableThinMPolePass
   Accelerator Toolbox
   S.White
*/

#include "atconstants.h"
#include "atelem.c"
#include "atlalib.c"
#include "atrandom.c"
#include "driftkick.c"

struct elemab {
    double* Amplitude;
    double Frequency;
    double Phase;
    double Sinmin, Sinmax;
    int NSamples;
    double* Func;
};

struct elem {
    double* PolynomA;
    double* PolynomB;
    double* PolynomAstart;
    double* PolynomBstart;
    struct elemab* ElemA;
    struct elemab* ElemB;
    int Mode;
    int MaxOrder;
    double* Ramps;
    int Periodic;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *EApertures;
    double *RApertures;
    double *Buffer;
    double BufferSize;
    double BufferOffset;
};

double get_amp(double amp, double* ramps, double t)
{
    double ampt = amp;
    if (ramps) {
        if (t <= ramps[0]) {
            ampt = 0.0;
        } else if (t <= ramps[1]) {
            ampt = amp * (t - ramps[0]) / (ramps[1] - ramps[0]);
        } else if (t <= ramps[2]) {
            ampt = amp;
        } else if (t <= ramps[3]) {
            ampt = amp - amp * (t - ramps[2]) / (ramps[3] - ramps[2]);
        } else {
            ampt = 0.0;
        }
    }
    return ampt;
}

double get_pol(
      struct elemab* elem,
      double* ramps, int mode,
      double t,
      int turn,
      int seed,
      int order,
      int periodic,
      pcg32_random_t* rng
      )
{
    int idx;
    double ampt, freq, ph, val;
    double* func;
    double* amp = elem->Amplitude;
    if (!amp) {
        return 0.0;
    }
    ampt = get_amp(amp[order], ramps, turn);
    switch (mode) {
    case 0:
        freq = elem->Frequency;
        ph = elem->Phase;
        val = sin(TWOPI * freq * t + ph);
        if (val < elem->Sinmin) val = elem->Sinmin;
        if (val > elem->Sinmax) val = elem->Sinmax;
        ampt *= val;
        return ampt;
    case 1:
        val = atrandn_r(rng, 0, 1);
        idx = turn - elem->BufferOffset;
        if (idx >= 0 && idx < elem->BufferSize) elem->Buffer[idx] = val;
        ampt *= val;
        return ampt;
    case 2:
        if (periodic || turn < elem->NSamples) {
            func = elem->Func;
            idx = turn % elem->NSamples;
            ampt *= func[idx];
            return ampt;
        } else {
            return 0.0;
        }
    default:
        return 0.0;
    }
}

void VariableThinMPolePass(
      double* r,
      struct elem* Elem,
      double t0,
      int turn, int num_particles,
      pcg32_random_t* rng // pcg32 random stream
    )
{

    int i, c;
    double* r6;
    double t = t0 * turn;

    int maxorder = Elem->MaxOrder;
    int periodic = Elem->Periodic;
    double* pola = Elem->PolynomA;
    double* polb = Elem->PolynomB;
    int mode = Elem->Mode;
    struct elemab* ElemA = Elem->ElemA;
    struct elemab* ElemB = Elem->ElemB;
    double* ramps = Elem->Ramps;


    // offsets at input and output
    double *T1 = Elem->T1;
    double *T2 = Elem->T2;
    // rotations at input and output
    double *R1 = Elem->R1;
    double *R2 = Elem->R2;
    // apertures
    double *RApertures = Elem->RApertures;
    double *EApertures = Elem->EApertures;

    if (mode != 0) {
        for (i = 0; i < maxorder + 1; i++) {
            pola[i] = get_pol(ElemA, ramps, mode, t, turn, seed, i, periodic, rng);
            polb[i] = get_pol(ElemB, ramps, mode, t, turn, seed, i, periodic, rng);
        };
    };

    for (c = 0; c < num_particles; c++) {
        r6 = r + c * 6;
        if (!atIsNaN(r6[0])) {
            if (mode == 0) {
                double tpart = t + r6[5] / C0;
                for (i = 0; i < maxorder + 1; i++) {
                    pola[i] = get_pol(ElemA, ramps, mode, tpart, turn, seed, i, periodic);
                    polb[i] = get_pol(ElemB, ramps, mode, tpart, turn, seed, i, periodic);
                };
            };
            /*  misalignment at entrance  */
            if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            strthinkick(r6, pola, polb, 1.0, maxorder);
            /* Misalignment at exit */
            if (R2) ATmultmv(r6,R2);
            if (T2) ATaddvv(r6,T2);
        }
    }

    /* restore polynoms to initial value. */
    for (i = 0; i < maxorder + 1; i++) {
        pola[i] = Elem->PolynomAstart[i];
        polb[i] = Elem->PolynomBstart[i];
    }
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem* trackFunction(const atElem* ElemData, struct elem* Elem,
    double* r_in, int num_particles, struct parameters* Param)
{
    if (!Elem) {
        int MaxOrder, Mode, NSamplesA, NSamplesB, Periodic;
        double *R1, *R2, *T1, *T2, *EApertures, *RApertures;
        double *PolynomA, *PolynomB, *AmplitudeA, *AmplitudeB;
        double *Ramps, *FuncA, *FuncB;
        double FrequencyA, FrequencyB;
        double PhaseA, PhaseB;
        double Sinmin, Sinmax;
        struct elemab *ElemA, *ElemB;
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
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
        Sinmin=atGetOptionalDouble(ElemData,"Sinmin", -1.1); check_error();
        Sinmax=atGetOptionalDouble(ElemData,"Sinmax", 1.1); check_error();
        Ramps=atGetOptionalDoubleArray(ElemData, "Ramps"); check_error();
        NSamplesA=atGetOptionalLong(ElemData, "NSamplesA", 1); check_error();
        NSamplesB=atGetOptionalLong(ElemData, "NSamplesB", 1); check_error();
        FuncA=atGetOptionalDoubleArray(ElemData,"FuncA"); check_error();
        FuncB=atGetOptionalDoubleArray(ElemData,"FuncB"); check_error();
        BufferA=atGetOptionalDoubleArray(ElemData,"BufferA"); check_error();
        BufferB=atGetOptionalDoubleArray(ElemData,"BufferB"); check_error();
        BufferSizeA=atGetOptionalLong(ElemData, "BufferSizeA", 0); check_error();
        BufferSizeB=atGetOptionalLong(ElemData, "BufferSizeB", 0); check_error();
        Periodic=atGetOptionalLong(ElemData,"Periodic", 1); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        ElemA = (struct elemab*)atMalloc(sizeof(struct elemab));
        ElemB = (struct elemab*)atMalloc(sizeof(struct elemab));
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
        Elem->PolynomA = PolynomA;
        Elem->PolynomB = PolynomB;
        Elem->PolynomAstart = (double *) atMalloc(sizeof(double *));
        Elem->PolynomBstart = (double *) atMalloc(sizeof(double *));
        memcpy(Elem->PolynomAstart, Elem->PolynomA, (MaxOrder+1)*sizeof(double));
        memcpy(Elem->PolynomBstart, Elem->PolynomB, (MaxOrder+1)*sizeof(double));
        Elem->Ramps = Ramps;
        Elem->Mode = Mode;
        Elem->MaxOrder = MaxOrder;
        Elem->Periodic = Periodic;
        ElemA->Amplitude = AmplitudeA;
        ElemB->Amplitude = AmplitudeB;
        ElemA->Frequency = FrequencyA;
        ElemB->Frequency = FrequencyB;
        ElemA->Phase = PhaseA;
        ElemB->Phase = PhaseB;
        ElemA->Sinmin = Sinmin;
        ElemB->Sinmin = Sinmin;
        ElemA->Sinmax = Sinmax;
        ElemB->Sinmax = Sinmax;
        ElemA->NSamples = NSamplesA;
        ElemB->NSamples = NSamplesB;
        ElemA->Func = FuncA;
        ElemB->Func = FuncB;
        Elem->ElemA = ElemA;
        Elem->ElemB = ElemB;
    }
    double t0 = Param->T0;
    int turn = Param->nturn;
    VariableThinMPolePass(r_in, Elem, t0, turn, num_particles, Param->thread_rng);
    return Elem;
}

MODULE_DEF(VariableThinMPolePass) /* Dummy module initialisation */

#endif /*defined(MATLAB_MEX_FILE) || defined(PYAT)*/

#ifdef MATLAB_MEX_FILE
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    if (nrhs >= 2) {
        double* r_in;
        const mxArray* ElemData = prhs[0];
        int num_particles = mxGetN(prhs[1]);
        int MaxOrder, Mode, NSamplesA, NSamplesB, Periodic;
        double *R1, *R2, *T1, *T2, *EApertures, *RApertures;
        double *PolynomA, *PolynomB, *AmplitudeA, *AmplitudeB;
        double *Ramps, *FuncA, *FuncB;
        double FrequencyA, FrequencyB;
        double PhaseA, PhaseB;
        double Sinmin, Sinmax;
        struct elemab ElA, *ElemA = &ElA;
        struct elemab ElB, *ElemB = &ElB;
        struct elem El, *Elem = &El;
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
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
        Sinmin=atGetOptionalDouble(ElemData,"Sinmin", -1.1); check_error();
        Sinmax=atGetOptionalDouble(ElemData,"Sinmax", 1.1); check_error();
        Ramps=atGetOptionalDoubleArray(ElemData, "Ramps"); check_error();
        NSamplesA=atGetOptionalLong(ElemData, "NSamplesA", 0); check_error();
        NSamplesB=atGetOptionalLong(ElemData, "NSamplesB", 0); check_error();
        FuncA=atGetOptionalDoubleArray(ElemData,"FuncA"); check_error();
        FuncB=atGetOptionalDoubleArray(ElemData,"FuncB"); check_error();
        Periodic=atGetOptionalLong(ElemData,"Periodic", 1); check_error();
        Elem->PolynomA = PolynomA;
        Elem->PolynomB = PolynomB;
        Elem->PolynomAstart = (double *) atMalloc(sizeof(double *));
        Elem->PolynomBstart = (double *) atMalloc(sizeof(double *));
        memcpy(Elem->PolynomAstart, Elem->PolynomA, (MaxOrder+1)*sizeof(double));
        memcpy(Elem->PolynomBstart, Elem->PolynomB, (MaxOrder+1)*sizeof(double));
        Elem->Ramps = Ramps;
        Elem->Mode = Mode;
        Elem->MaxOrder = MaxOrder;
        Elem->Periodic = Periodic;
        Elem->R1 = R1;
        Elem->R2 = R2;
        Elem->T1 = T1;
        Elem->T2 = T2;
        Elem->EApertures = EApertures;
        Elem->RApertures = RApertures;
        ElemA->Amplitude = AmplitudeA;
        ElemB->Amplitude = AmplitudeB;
        ElemA->Frequency = FrequencyA;
        ElemB->Frequency = FrequencyB;
        ElemA->Phase = PhaseA;
        ElemB->Phase = PhaseB;
        ElemA->Sinmin = Sinmin;
        ElemB->Sinmin = Sinmin;
        ElemA->Sinmax = Sinmax;
        ElemB->Sinmax = Sinmax;
        ElemA->NSamples = NSamplesA;
        ElemB->NSamples = NSamplesB;
        ElemA->Func = FuncA;
        ElemB->Func = FuncB;
        Elem->ElemA = ElemA;
        Elem->ElemB = ElemB;
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        VariableThinMPolePass(r_in, Elem, 0, 0, num_particles, &pcg32_global);
    } else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(4, 1);
        mxSetCell(plhs[0], 0, mxCreateString("MaxOrder"));
        mxSetCell(plhs[0], 1, mxCreateString("Mode"));
        mxSetCell(plhs[0], 2, mxCreateString("PolynomA"));
        mxSetCell(plhs[0], 3, mxCreateString("PolynomB"));
        if (nlhs > 1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(22, 1);
            mxSetCell(plhs[1], 0, mxCreateString("AmplitudeA"));
            mxSetCell(plhs[1], 1, mxCreateString("AmplitudeB"));
            mxSetCell(plhs[1], 2, mxCreateString("FrequencyA"));
            mxSetCell(plhs[1], 3, mxCreateString("FrequencyB"));
            mxSetCell(plhs[1], 4, mxCreateString("PhaseA"));
            mxSetCell(plhs[1], 5, mxCreateString("PhaseB"));
            mxSetCell(plhs[1], 6, mxCreateString("Ramps"));
            mxSetCell(plhs[1], 7, mxCreateString("FuncA"));
            mxSetCell(plhs[1], 8, mxCreateString("FuncB"));
            mxSetCell(plhs[1], 9, mxCreateString("NSamplesA"));
            mxSetCell(plhs[1], 10, mxCreateString("NSamplesB"));
            mxSetCell(plhs[1], 11, mxCreateString("Periodic"));
            mxSetCell(plhs[1], 12, mxCreateString("T1"));
            mxSetCell(plhs[1], 13, mxCreateString("T2"));
            mxSetCell(plhs[1], 14, mxCreateString("R1"));
            mxSetCell(plhs[1], 15, mxCreateString("R2"));
            mxSetCell(plhs[1], 16, mxCreateString("RApertures"));
            mxSetCell(plhs[1], 17, mxCreateString("EApertures"));
            mxSetCell(plhs[1], 18, mxCreateString("Sinmin"));
            mxSetCell(plhs[1], 19, mxCreateString("Sinmax"));
            mxSetCell(plhs[1], 20, mxCreateString("Buffer"));
            mxSetCell(plhs[1], 21, mxCreateString("Buffersize"));
        }
    } else {
        mexErrMsgIdAndTxt("AT:WrongArg", "Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
