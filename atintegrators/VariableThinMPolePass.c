/* VariableThinMPolePass
   Accelerator Toolbox
   S.White
*/

/* 2024nov19 oblanco at ALBA CELLS. modified to be compatible with
                                    pyat and AT matlab*/

/*
 This pass method is implements three modes.
 Mode 0 : a sin function with frequency and phase
 Mode 1 : a random kick
 Mode 2 : a custom function defined in every turn

 All modes could be ramped using the flag `ramp:.

 The value of the polynoms A and B are calculated per turn during tracking,
 and also, in modes 0 and 2, per particle to take in to account
 the delay time.

 In mode 0, the sin function could be limited to be between any value above
 `Sinlimit`. For example a half-sin function would be obtained by setting
 `Sinlimit` to zero.

 In mode 1 the stream of pseudo-random values is taken from the
 parameters structure in attypes.h. For more details about the random
 generator in AT:
 https://github.com/atcollab/at/discussions/879

 In mode 2, the values could be periodically applied or not.
*/


#include "atconstants.h"
#include "atelem.c"
#include "atlalib.c"
#include "atrandom.c"
#include "driftkick.c"

/* constants to be used in the Taylor expansion of the custom function */
#define oneoversix 0.166666666666667
#define oneovertwentyfour 0.041666666666667

/* This struct contains the values to set one of the two
   poynoms: A or B */
struct elemab {
    double* Amplitude;
    double Frequency;
    double Sinlimit;
    double Phase;
    int NSamples;
    double* Func;
    double* Funcderiv1;
    double* Funcderiv2;
    double* Funcderiv3;
    double* Funcderiv4;
    double FuncTimeDelay;
};

/* This struct contains the parameters of the element.
 * It uses elemab struct */
struct elem {
    struct elemab* ElemA;
    struct elemab* ElemB;
    double* PolynomA;        // calculated on every turn
    double* PolynomB;        // calculated on every turn
    int MaxOrder;
    int Mode;
    double* Ramps;
    int Periodic;
    double *R1;
    double *R2;
    double *T1;
    double *T2;
    double *EApertures;
    double *RApertures;
};


/* get_amp returns the input value `amp` when ramps is False.
 * If ramps is True, it returns a value linearly interpolated
 * accoding to the ramping turn.*/
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


/* get_pol calculates the PolynomA/B per turn and mode.
   If the mode is 0 or 2, the polynom depends on the particle time
   of delay*/
double get_pol(
    struct elemab* elem,
    double* ramps,
    int mode,
    double t,           // time
    int turn,
    int order,          //max(length of any of the polynoms) - 1
    int periodic,
    pcg32_random_t* rng
    )
{
    int turnidx; // turn index
    double ampt; // amplitude per turn

    // sin mode parameters
    double whole_sin_limit = elem->Sinlimit;
    double freq, ph, sinval;

    // custom mode parameters
    double* func;
    double *funcderiv1, *funcderiv2, *funcderiv3, *funcderiv4;
    double* amp = elem->Amplitude;
    double t2;       // time squared
    double functdelay;

    if (!amp) {
        return 0.0;
    }
    ampt = get_amp(amp[order], ramps, turn);

    switch (mode) {
    case 0:
        freq = elem->Frequency;
        ph = elem->Phase;
        sinval = sin(TWOPI * freq * t + ph);
        // branchless if
        // sinval >= wholelimit  -> sinval
        // else                  -> 0
        ampt = ampt*sinval*(sinval >= whole_sin_limit);
        return ampt;
    case 1:
        ampt *= atrandn_r(rng, 0, 1);
        return ampt;
    case 2:
        if (periodic || turn < elem->NSamples) {
            func = elem->Func;
            funcderiv1 = elem->Funcderiv1;
            funcderiv2 = elem->Funcderiv2;
            funcderiv3 = elem->Funcderiv3;
            funcderiv4 = elem->Funcderiv4;
            functdelay = elem->FuncTimeDelay;
            turnidx = turn % elem->NSamples;

            t = t - functdelay;
            t2 = t*t;
            ampt = ampt*(func[turnidx]
                  + funcderiv1[turnidx]*t
                  + 0.5*funcderiv2[turnidx]*t2
                  + oneoversix*funcderiv3[turnidx]*t2*t
                  + oneovertwentyfour*funcderiv4[turnidx]*t2*t2);
            return ampt;
        } else {
            return 0.0;
        }
    default:
        return 0.0;
    }
}

/* This function tracks a particle through a thin element with
   variable PolynomB and PolynomA per turn.*/
void VariableThinMPolePass(
    double* r, // the particle in 6D
    struct elem* Elem,  //the variable element
    double t0, //time of one revolution
    int turn,  //number of the current turn
    int num_particles, //number of particles
    pcg32_random_t* rng // pcg32 random stream
    )
{

    int i; // order of the polynom, counter
    int c; // particle, counter
    double* r6; // particle 6D coordinates
    double tpart; //time of particle delay

    // particle time offset for the mode
    double time_in_this_mode = 0;

    // setting the element properties
    int mode = Elem->Mode;
    struct elemab* ElemA = Elem->ElemA;
    struct elemab* ElemB = Elem->ElemB;
    double* pola = Elem->PolynomA;
    double* polb = Elem->PolynomB;
    int maxorder = Elem->MaxOrder;
    double* ramps = Elem->Ramps;

    // custom function is periodic
    int periodic = Elem->Periodic;

    // offsets at input and output
    double *T1 = Elem->T1;
    double *T2 = Elem->T2;
    // rotations at input and output
    double *R1 = Elem->R1;
    double *R2 = Elem->R2;
    // apertures
    double *RApertures = Elem->RApertures;
    double *EApertures = Elem->EApertures;

    /* mode 0 : sin function */
    /* mode 1 : random value applied to all particles */
    /* mode 2 : custom function */

    if (mode == 1) {
         /* calculate the polynom to apply on all particles */
         /* the particle delay time is not used */
        for (i = 0; i < maxorder + 1; i++) {
            pola[i] = get_pol(ElemA, ramps, mode, 0, turn, i, periodic, rng);
            polb[i] = get_pol(ElemB, ramps, mode, 0, turn, i, periodic, rng);
        };
    };

    // offset the time when applying the sin function
    // branchless if
    // mode == 0 -> time_in_this_mode = t0 *turn
    // else      -> 0
    time_in_this_mode = t0 * turn * (mode == 0);

    /* cycle over all particles */
    for (c = 0; c < num_particles; c++) {
        r6 = r + c * 6;
        /* check if the particle is alive */
        if (!atIsNaN(r6[0])) {
            /* mode 0 and mode 2 take into account the particle delay time */
            if (mode == 0 || mode == 2) {
                /* modify the time of delay of the particle  */
                tpart = time_in_this_mode + r6[5] / C0;
                /* calculate the polynom A and B components seen by the particle */
                for (i = 0; i < maxorder + 1; i++) {
                    pola[i] = get_pol(ElemA, ramps, mode, tpart, turn, i, periodic, rng);
                    polb[i] = get_pol(ElemB, ramps, mode, tpart, turn, i, periodic, rng);
                };
            };
            /*  misalignment at entrance  */
            if (T1) ATaddvv(r6,T1);
            if (R1) ATmultmv(r6,R1);
            /* Check physical apertures at the entrance of the magnet */
            if (RApertures) checkiflostRectangularAp(r6,RApertures);
            if (EApertures) checkiflostEllipticalAp(r6,EApertures);
            /* track */
            strthinkick(r6, pola, polb, 1.0, maxorder);
            /* Misalignment at exit */
            if (R2) ATmultmv(r6,R2);
            if (T2) ATaddvv(r6,T2);
        }
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
        double *FuncAderiv1, *FuncBderiv1;
        double *FuncAderiv2, *FuncBderiv2;
        double *FuncAderiv3, *FuncBderiv3;
        double *FuncAderiv4, *FuncBderiv4;
        double FuncATimeDelay, FuncBTimeDelay;
        double FrequencyA, FrequencyB;
        double PhaseA, PhaseB;
        double SinlimitA, SinlimitB;
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
        SinlimitA=atGetOptionalDouble(ElemData,"SinlimitA", 0); check_error();
        SinlimitB=atGetOptionalDouble(ElemData,"SinlimitB", 0); check_error();
        Ramps=atGetOptionalDoubleArray(ElemData, "Ramps"); check_error();
        NSamplesA=atGetOptionalLong(ElemData, "NSamplesA", 1); check_error();
        NSamplesB=atGetOptionalLong(ElemData, "NSamplesB", 1); check_error();
        FuncA=atGetOptionalDoubleArray(ElemData,"FuncA"); check_error();
        FuncB=atGetOptionalDoubleArray(ElemData,"FuncB"); check_error();
        FuncAderiv1=atGetOptionalDoubleArray(ElemData,"FuncAderiv1"); check_error();
        FuncBderiv1=atGetOptionalDoubleArray(ElemData,"FuncBderiv1"); check_error();
        FuncAderiv2=atGetOptionalDoubleArray(ElemData,"FuncAderiv2"); check_error();
        FuncBderiv2=atGetOptionalDoubleArray(ElemData,"FuncBderiv2"); check_error();
        FuncAderiv3=atGetOptionalDoubleArray(ElemData,"FuncAderiv3"); check_error();
        FuncBderiv3=atGetOptionalDoubleArray(ElemData,"FuncBderiv3"); check_error();
        FuncAderiv4=atGetOptionalDoubleArray(ElemData,"FuncAderiv4"); check_error();
        FuncBderiv4=atGetOptionalDoubleArray(ElemData,"FuncBderiv4"); check_error();
        FuncATimeDelay=atGetOptionalDouble(ElemData,"FuncATimeDelay", 0); check_error();
        FuncBTimeDelay=atGetOptionalDouble(ElemData,"FuncBTimeDelay", 0); check_error();
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
        ElemA->Sinlimit = SinlimitA;
        ElemB->Sinlimit = SinlimitB;
        ElemA->NSamples = NSamplesA;
        ElemB->NSamples = NSamplesB;
        ElemA->Func = FuncA;
        ElemB->Func = FuncB;
        ElemA->Funcderiv1 = FuncAderiv1;
        ElemB->Funcderiv1 = FuncBderiv1;
        ElemA->Funcderiv2 = FuncAderiv2;
        ElemB->Funcderiv2 = FuncBderiv2;
        ElemA->Funcderiv3 = FuncAderiv3;
        ElemB->Funcderiv3 = FuncBderiv3;
        ElemA->Funcderiv4 = FuncAderiv4;
        ElemB->Funcderiv4 = FuncBderiv4;
        ElemA->FuncTimeDelay = FuncATimeDelay;
        ElemB->FuncTimeDelay = FuncBTimeDelay;
        Elem->ElemA = ElemA;
        Elem->ElemB = ElemB;
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
        double *PolynomA, *PolynomB, *AmplitudeA, *AmplitudeB;
        double *Ramps, *FuncA, *FuncB;
        double *FuncAderiv1, *FuncBderiv1;
        double *FuncAderiv2, *FuncBderiv2;
        double *FuncAderiv3, *FuncBderiv3;
        double *FuncAderiv4, *FuncBderiv4;
        double FuncATimeDelay, FuncBTimeDelay;
        double FrequencyA, FrequencyB;
        double PhaseA, PhaseB;
        double SinlimitA, SinlimitB;
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
        SinlimitA=atGetOptionalDouble(ElemData,"SinlimitA", 0); check_error();
        SinlimitB=atGetOptionalDouble(ElemData,"SinlimitB", 0); check_error();
        Ramps=atGetOptionalDoubleArray(ElemData, "Ramps"); check_error();
        NSamplesA=atGetOptionalLong(ElemData, "NSamplesA", 0); check_error();
        NSamplesB=atGetOptionalLong(ElemData, "NSamplesB", 0); check_error();
        FuncA=atGetOptionalDoubleArray(ElemData,"FuncA"); check_error();
        FuncB=atGetOptionalDoubleArray(ElemData,"FuncB"); check_error();
        FuncAderiv1=atGetOptionalDoubleArray(ElemData,"FuncAderiv1"); check_error();
        FuncBderiv1=atGetOptionalDoubleArray(ElemData,"FuncBderiv1"); check_error();
        FuncAderiv2=atGetOptionalDoubleArray(ElemData,"FuncAderiv2"); check_error();
        FuncBderiv2=atGetOptionalDoubleArray(ElemData,"FuncBderiv2"); check_error();
        FuncAderiv3=atGetOptionalDoubleArray(ElemData,"FuncAderiv3"); check_error();
        FuncBderiv3=atGetOptionalDoubleArray(ElemData,"FuncBderiv3"); check_error();
        FuncAderiv4=atGetOptionalDoubleArray(ElemData,"FuncAderiv4"); check_error();
        FuncBderiv4=atGetOptionalDoubleArray(ElemData,"FuncBderiv4"); check_error();
        FuncATimeDelay=atGetOptionalDouble(ElemData,"FuncATimeDelay", 0); check_error();
        FuncBTimeDelay=atGetOptionalDouble(ElemData,"FuncBTimeDelay", 0); check_error();
        Periodic=atGetOptionalLong(ElemData,"Periodic", 1); check_error();
        Elem->PolynomA = PolynomA;
        Elem->PolynomB = PolynomB;
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
        ElemA->Sinlimit = SinlimitA;
        ElemB->Sinlimit = SinlimitB;
        ElemA->NSamples = NSamplesA;
        ElemB->NSamples = NSamplesB;
        ElemA->Func = FuncA;
        ElemB->Func = FuncB;
        Elem->ElemA = ElemA;
        Elem->ElemB = ElemB;
        ElemA->Funcderiv1 = FuncAderiv1;
        ElemB->Funcderiv1 = FuncBderiv1;
        ElemA->Funcderiv2 = FuncAderiv2;
        ElemB->Funcderiv2 = FuncBderiv2;
        ElemA->Funcderiv3 = FuncAderiv3;
        ElemB->Funcderiv3 = FuncBderiv3;
        ElemA->Funcderiv4 = FuncAderiv4;
        ElemB->Funcderiv4 = FuncBderiv4;
        ElemA->FuncTimeDelay = FuncATimeDelay;
        ElemB->FuncTimeDelay = FuncBTimeDelay;
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
        if (nlhs > 4) {
            /* list of optional fields */
            plhs[1] = mxCreateCellMatrix(23, 1);
            mxSetCell(plhs[1], 0, mxCreateString("AmplitudeA"));
            mxSetCell(plhs[1], 1, mxCreateString("AmplitudeB"));
            mxSetCell(plhs[1], 2, mxCreateString("FrequencyA"));
            mxSetCell(plhs[1], 3, mxCreateString("FrequencyB"));
            mxSetCell(plhs[1], 4, mxCreateString("PhaseA"));
            mxSetCell(plhs[1], 5, mxCreateString("PhaseB"));
            mxSetCell(plhs[1], 6, mxCreateString("SinlimitA"));
            mxSetCell(plhs[1], 7, mxCreateString("SinlimitB"));
            mxSetCell(plhs[1], 8, mxCreateString("Ramps"));
            mxSetCell(plhs[1], 9, mxCreateString("FuncA"));
            mxSetCell(plhs[1], 10, mxCreateString("FuncB"));
            mxSetCell(plhs[1], 11, mxCreateString("FuncAderiv1"));
            mxSetCell(plhs[1], 12, mxCreateString("FuncBderiv1"));
            mxSetCell(plhs[1], 13, mxCreateString("FuncAderiv2"));
            mxSetCell(plhs[1], 14, mxCreateString("FuncBderiv2"));
            mxSetCell(plhs[1], 15, mxCreateString("FuncAderiv3"));
            mxSetCell(plhs[1], 16, mxCreateString("FuncBderiv3"));
            mxSetCell(plhs[1], 17, mxCreateString("FuncAderiv4"));
            mxSetCell(plhs[1], 18, mxCreateString("FuncBderiv4"));
            mxSetCell(plhs[1], 19, mxCreateString("FuncATimeDelay"));
            mxSetCell(plhs[1], 20, mxCreateString("FuncBTimeDelay"));
            mxSetCell(plhs[1], 21, mxCreateString("NSamplesA"));
            mxSetCell(plhs[1], 22, mxCreateString("NSamplesB"));
            mxSetCell(plhs[1], 23, mxCreateString("Periodic"));
            mxSetCell(plhs[1], 24,mxCreateString("T1"));
            mxSetCell(plhs[1], 25,mxCreateString("T2"));
            mxSetCell(plhs[1], 26,mxCreateString("R1"));
            mxSetCell(plhs[1], 27,mxCreateString("R2"));
            mxSetCell(plhs[1], 28,mxCreateString("RApertures"));
            mxSetCell(plhs[1], 29,mxCreateString("EApertures"));
        }
    } else {
        mexErrMsgIdAndTxt("AT:WrongArg", "Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
