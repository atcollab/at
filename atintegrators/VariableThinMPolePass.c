/* VariableThinMPolePass
   Accelerator Toolbox
   S.White
*/

#include "atconstants.h"
#include "atelem.c"
#include "atlalib.c"
#include "atrandom.c"
#include "driftkick.c"
#include "interpolate.c"
#include <math.h>

struct elemab {
    double* Amplitude;
    double Frequency;
    double Phase;
    double Sinmin, Sinmax;
    int NSamples;
    double* Func;
    double* Finterpolate;
    double* Tinterpolate;
};

struct elem {
    struct elemab ElemA;
    struct elemab ElemB;
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

double get_val(struct elemab* elem, double* ramps, int mode,
    double t, int turn, int order, int periodic, pcg32_random_t* rng)
{
    int idx;
    double ampt, freq, ph, val;
    double* func;
    double* titp;
    double* fitp;
    int nsamples = elem->NSamples;
    ampt = get_amp(1.0, ramps, turn);
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
        ampt *= atrandn_r(rng, 0, 1);
        return ampt;
    case 2:
        if (periodic || turn < nsamples) {
            func = elem->Func;
            idx = turn % nsamples;
            ampt *= func[idx];
            return ampt;
        } else {
            return 0.0;
        }
    case 3:
       titp = elem->Tinterpolate;
       fitp = elem->Finterpolate;
       if (periodic){
         while (t < titp[0]){t = t+titp[nsamples-1];};
         t = fmod(t, titp[nsamples-1]);
       };
       idx = binarySearch(titp , t, nsamples, 0, 0);
       /* checking if t is outside the range of titp */
       if (t < titp[0]){
         val = fitp[0];
       }else if(t > titp[nsamples-1]){
         val = fitp[nsamples-1];
       }else{
         val = interpolTable(fitp, titp, t, idx);
       };
       return ampt *= val;
    default:
        return 0.0;
    }
}

void VariableThinMPolePass(double* r, struct elem* Elem, double t0, int turn, int num_particles,
    pcg32_random_t* rng)
{

    int i, c;
    double* r6;
    double t = t0 * turn;
    double tpart;
    double vala, valb;

    int maxorder = Elem->MaxOrder;
    int periodic = Elem->Periodic;
    int mode = Elem->Mode;
    struct elemab* ElemA = &(Elem->ElemA);
    struct elemab* ElemB = &(Elem->ElemB);
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

    // create thread safe polynoms. Each thread will allocate and free pola and polb
    double *pola = (double *) atCalloc((maxorder+1),sizeof(double));
    double *polb = (double *) atCalloc((maxorder+1),sizeof(double));

    /* mode 0 : sin function */
    /* mode 1 : random value applied to all particles */
    /* mode 2 : custom function */
    /* mode 3 : interpolate */

    if (mode == 1) {
        if (ElemA->Amplitude){
            vala = get_val(ElemA, ramps, mode, 0, turn, i, periodic, rng);
            for (i = 0; i < maxorder + 1; i++) pola[i] = vala * ElemA->Amplitude[i];
        };
        if (ElemB->Amplitude){
            valb = get_val(ElemB, ramps, mode, 0, turn, i, periodic, rng);
            for (i = 0; i < maxorder + 1; i++) polb[i] = valb * ElemB->Amplitude[i];
        };
    };

    for (c = 0; c < num_particles; c++) {
        r6 = r + c * 6;
        if (!atIsNaN(r6[0])) {
            if (mode != 1){
                tpart = t*(mode == 0) + t*(mode == 3) + r6[5] / C0;
                if (ElemA->Amplitude){
                    vala = get_val(ElemA, ramps, mode, tpart, turn, i, periodic, rng);
                    for (i = 0; i < maxorder + 1; i++) pola[i]=vala*ElemA->Amplitude[i];
                };
                if (ElemB->Amplitude){
                    valb = get_val(ElemB, ramps, mode, tpart, turn, i, periodic, rng);
                    for (i = 0; i < maxorder + 1; i++) polb[i]=valb*ElemB->Amplitude[i];
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

    atFree(pola);
    atFree(polb);
}

#if defined(MATLAB_MEX_FILE) || defined(PYAT)
ExportMode struct elem* trackFunction(const atElem* ElemData, struct elem* Elem,
    double* r_in, int num_particles, struct parameters* Param)
{
    if (!Elem) {
        int MaxOrder, Mode, NSamplesA, NSamplesB, Periodic;
        double *R1, *R2, *T1, *T2, *EApertures, *RApertures;
        double *AmplitudeA, *AmplitudeB;
        double *Ramps, *FuncA, *FuncB;
        double *FinterpolateA, *FinterpolateB;
        double *TinterpolateA, *TinterpolateB;
        double FrequencyA, FrequencyB;
        double PhaseA, PhaseB;
        double Sinmin, Sinmax;
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        Mode=atGetLong(ElemData,"Mode"); check_error();
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
        FinterpolateA=atGetOptionalDoubleArray(ElemData,"FinterpolateA"); check_error();
        FinterpolateB=atGetOptionalDoubleArray(ElemData,"FinterpolateB"); check_error();
        TinterpolateA=atGetOptionalDoubleArray(ElemData,"TinterpolateA"); check_error();
        TinterpolateB=atGetOptionalDoubleArray(ElemData,"TinterpolateB"); check_error();
        Periodic=atGetOptionalLong(ElemData,"Periodic", 1); check_error();
        Elem = (struct elem*)atMalloc(sizeof(struct elem));
        struct elemab* ElemA = &(Elem->ElemA);
        struct elemab* ElemB = &(Elem->ElemB);
        Elem->R1=R1;
        Elem->R2=R2;
        Elem->T1=T1;
        Elem->T2=T2;
        Elem->EApertures=EApertures;
        Elem->RApertures=RApertures;
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
        ElemA->Finterpolate = FinterpolateA;
        ElemB->Finterpolate = FinterpolateB;
        ElemA->Tinterpolate = TinterpolateA;
        ElemB->Tinterpolate = TinterpolateB;
    }
    double t0 = Param->T0;
    int turn = Param->nturn;
    VariableThinMPolePass(r_in, Elem, t0, turn, num_particles, Param->common_rng);
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
        double *AmplitudeA, *AmplitudeB;
        double *Ramps, *FuncA, *FuncB;
        double *FinterpolateA, *FinterpolateB;
        double *TinterpolateA, *TinterpolateB;
        double FrequencyA, FrequencyB;
        double PhaseA, PhaseB;
        double Sinmin, Sinmax;
        struct elem El, *Elem = &El;
        struct elemab* ElemA = &(Elem->ElemA);
        struct elemab* ElemB = &(Elem->ElemB);
        R1=atGetOptionalDoubleArray(ElemData,"R1"); check_error();
        R2=atGetOptionalDoubleArray(ElemData,"R2"); check_error();
        T1=atGetOptionalDoubleArray(ElemData,"T1"); check_error();
        T2=atGetOptionalDoubleArray(ElemData,"T2"); check_error();
        EApertures=atGetOptionalDoubleArray(ElemData,"EApertures"); check_error();
        RApertures=atGetOptionalDoubleArray(ElemData,"RApertures"); check_error();
        MaxOrder=atGetLong(ElemData,"MaxOrder"); check_error();
        Mode=atGetLong(ElemData,"Mode"); check_error();
        AmplitudeA=atGetOptionalDoubleArray(ElemData,"AmplitudeA"); check_error();
        AmplitudeB=atGetOptionalDoubleArray(ElemData,"AmplitudeB"); check_error();
        FrequencyA=atGetOptionalDouble(ElemData,"FrequencyA", 0); check_error();
        FrequencyB=atGetOptionalDouble(ElemData,"FrequencyB", 0); check_error();
        PhaseA=atGetOptionalDouble(ElemData,"PhaseA", 0); check_error();
        PhaseB=atGetOptionalDouble(ElemData,"PhaseB", 0); check_error();
        Sinmin=atGetOptionalDouble(ElemData,"Sinmin", -1.1); check_error();
        Sinmax=atGetOptionalDouble(ElemData,"Sinmax", 1.1); check_error();
        Ramps=atGetOptionalDoubleArray(ElemData,"Ramps"); check_error();
        NSamplesA=atGetOptionalLong(ElemData,"NSamplesA", 0); check_error();
        NSamplesB=atGetOptionalLong(ElemData,"NSamplesB", 0); check_error();
        FuncA=atGetOptionalDoubleArray(ElemData,"FuncA"); check_error();
        FuncB=atGetOptionalDoubleArray(ElemData,"FuncB"); check_error();
        FinterpolateA=atGetOptionalDoubleArray(ElemData,"FinterpolateA"); check_error();
        FinterpolateB=atGetOptionalDoubleArray(ElemData,"FinterpolateB"); check_error();
        TinterpolateA=atGetOptionalDoubleArray(ElemData,"TinterpolateA"); check_error();
        TinterpolateB=atGetOptionalDoubleArray(ElemData,"TinterpolateB"); check_error();
        Periodic=atGetOptionalLong(ElemData,"Periodic", 1); check_error();
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
        ElemA->Finterpolate = FinterpolateA;
        ElemB->Finterpolate = FinterpolateB;
        ElemA->Tinterpolate = TinterpolateA;
        ElemB->Tinterpolate = TinterpolateB;
        /* ALLOCATE memory for the output array of the same size as the input  */
        plhs[0] = mxDuplicateArray(prhs[1]);
        r_in = mxGetDoubles(plhs[0]);
        VariableThinMPolePass(r_in, Elem, 0, 0, num_particles, &pcg32_global);
    } else if (nrhs == 0) {
        /* list of required fields */
        plhs[0] = mxCreateCellMatrix(2, 1);
        mxSetCell(plhs[0], 0, mxCreateString("MaxOrder"));
        mxSetCell(plhs[0], 1, mxCreateString("Mode"));
        if (nlhs > 1) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(24, 1);
            mxSetCell(plhs[1], 0, mxCreateString("AmplitudeA"));
            mxSetCell(plhs[1], 1, mxCreateString("AmplitudeB"));
            mxSetCell(plhs[1], 2, mxCreateString("FrequencyA"));
            mxSetCell(plhs[1], 3, mxCreateString("FrequencyB"));
            mxSetCell(plhs[1], 4, mxCreateString("PhaseA"));
            mxSetCell(plhs[1], 5, mxCreateString("PhaseB"));
            mxSetCell(plhs[1], 6, mxCreateString("Ramps"));
            mxSetCell(plhs[1], 7, mxCreateString("FuncA"));
            mxSetCell(plhs[1], 8, mxCreateString("FuncB"));
            mxSetCell(plhs[1], 9, mxCreateString("FinterpolateA"));
            mxSetCell(plhs[1], 10, mxCreateString("FinterpolateB"));
            mxSetCell(plhs[1], 11, mxCreateString("TinterpolateA"));
            mxSetCell(plhs[1], 12, mxCreateString("TinterpolateB"));
            mxSetCell(plhs[1], 13, mxCreateString("NSamplesA"));
            mxSetCell(plhs[1], 14, mxCreateString("NSamplesB"));
            mxSetCell(plhs[1], 15, mxCreateString("Periodic"));
            mxSetCell(plhs[1], 16, mxCreateString("T1"));
            mxSetCell(plhs[1], 17, mxCreateString("T2"));
            mxSetCell(plhs[1], 18, mxCreateString("R1"));
            mxSetCell(plhs[1], 19, mxCreateString("R2"));
            mxSetCell(plhs[1], 20, mxCreateString("RApertures"));
            mxSetCell(plhs[1], 21, mxCreateString("EApertures"));
            mxSetCell(plhs[1], 22, mxCreateString("Sinmin"));
            mxSetCell(plhs[1], 23, mxCreateString("Sinmax"));
        }
    } else {
        mexErrMsgIdAndTxt("AT:WrongArg", "Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
