/* VariableThinMPolePass
   Accelerator Toolbox
   S.White
*/

/* 2024nov19 oblanco at ALBA CELLS. modified to be compatible with
                                    pyat and AT matlab */

/*
 This pass method implements three modes.
 Mode 0 : a sin function with frequency and phase
 Mode 1 : a random kick
 Mode 2 : a custom function defined in every turn

 All modes could be ramped using the flag `ramp:.

 The value of the polynoms A and B are calculated per turn during tracking,
 and also, in modes 0 and 2, per particle to take in to account
 the individual delay time.

 In mode 0, the sin function could be limited to be between any value above
 `Sinabove`. For example a positive half-sin function would be obtained by
 setting `Sinabove` to zero.

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

/* This struct contains the values to set one of the two
   polynoms: A or B , as a function of the element parameters. */
struct elemab {
    double* Amplitude;
    double Frequency;
    double Sinabove;
    double Phase;
    double* Func;
    double FuncTimeDelay;
    int NSamples;
    int Ktaylor;
    double* Buffer;
    int BufferSize;
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
   If the mode is 0 or 2, the polynom depends on the particle
   time delay*/
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
    int turnidx,i; // index
    double ampt;

    // sin mode parameters
    double whole_sin_above = elem->Sinabove;
    double freq, ph, sinval;

    // random mode parameters
    double random_value;

    // custom mode parameters
    double* func;
    double* amp = elem->Amplitude;
    double functdelay;
    double functot;
    double tpow,thefactorial; // variables in for cycles
    int Ktaylor;

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
        // sinval >= wholeabove  -> sinval
        // else                  -> 0
        ampt = ampt*sinval*(sinval >= whole_sin_above);
        return ampt;
    case 1:
        random_value = atrandn_r(rng, 0, 1);
        ampt *= random_value
        /* save random value into buffer */
        if (turn < elem->BufferSize){
            elem->Buffer[turn] = random_value;
        }
        return ampt;
    case 2:
        if (periodic || turn < elem->NSamples) {
            /* get the function, the delay, and turns */
            func = elem->Func;
            functdelay = elem->FuncTimeDelay;
            turnidx = turn % elem->NSamples;
            Ktaylor = elem->Ktaylor;

            /* calculate the amplitude as a Taylor expansion */
            /* first order is taken directly from the function table */
            t = t - functdelay;
            functot = func[Ktaylor*turnidx];
            tpow = 1;
            thefactorial = 1;
            /* do a taylor expansion if Ktaylor is more than first order */
            for (i=1;i<Ktaylor;i++){
              tpow = tpow * t;
              thefactorial = thefactorial * i;
              /* indexing is fortran-like. We start with columns.
                 cols are taylor components.
                 rows are turn samples.
              */
              functot =  functot + tpow / thefactorial * func[i + Ktaylor*turnidx];
            };
            ampt = ampt * functot;

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
    struct elem* Elem,  // the variable element
    double t0, // time of one revolution
    int turn,  // number of the current turn
    int num_particles, //number of particles
    pcg32_random_t* rng // pcg32 random stream
    )
{

    int i; // order of the polynom, counter
    int c; // particle, counter
    double* r6; // particle 6D coordinates
    double tpart; // particle time delay

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
    /* Initialize element */
    if (!Elem) {
        int MaxOrder, Mode, NSamplesA, NSamplesB, KtaylorA, KtaylorB, Periodic;
        int BufferSizeA, BufferSizeB;
        double *R1, *R2, *T1, *T2, *EApertures, *RApertures;
        double *PolynomA, *PolynomB, *AmplitudeA, *AmplitudeB;
        double *Ramps, *FuncA, *FuncB, *BufferA, *BufferB;
        double FuncATimeDelay, FuncBTimeDelay;
        double FrequencyA, FrequencyB;
        double PhaseA, PhaseB;
        double SinAabove, SinBabove;
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
        SinAabove=atGetOptionalDouble(ElemData,"SinAabove", 0); check_error();
        SinBabove=atGetOptionalDouble(ElemData,"SinBabove", 0); check_error();
        Ramps=atGetOptionalDoubleArray(ElemData, "Ramps"); check_error();
        NSamplesA=atGetOptionalLong(ElemData, "NSamplesA", 1); check_error();
        NSamplesB=atGetOptionalLong(ElemData, "NSamplesB", 1); check_error();
        KtaylorA=atGetOptionalLong(ElemData, "KtaylorA", 1); check_error();
        KtaylorB=atGetOptionalLong(ElemData, "KtaylorB", 1); check_error();
        FuncA=atGetOptionalDoubleArray(ElemData,"FuncA"); check_error();
        FuncB=atGetOptionalDoubleArray(ElemData,"FuncB"); check_error();
        BufferA=atGetOptionalDoubleArray(ElemData,"BufferA"); check_error();
        BufferB=atGetOptionalDoubleArray(ElemData,"BufferB"); check_error();
        BufferSizeA=atGetOptionalLong(ElemData, "BufferSizeA", 0); check_error();
        BufferSizeB=atGetOptionalLong(ElemData, "BufferSizeB", 0); check_error();
        FuncATimeDelay=atGetOptionalDouble(ElemData,"FuncATimeDelay", 0); check_error();
        FuncBTimeDelay=atGetOptionalDouble(ElemData,"FuncBTimeDelay", 0); check_error();
        Periodic=atGetOptionalLong(ElemData,"Periodic", 0); check_error();
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
        ElemA->Sinabove = SinAabove;
        ElemB->Sinabove = SinBabove;
        ElemA->NSamples = NSamplesA;
        ElemB->NSamples = NSamplesB;
        ElemA->Ktaylor = KtaylorA;
        ElemB->Ktaylor = KtaylorB;
        ElemA->Func = FuncA;
        ElemB->Func = FuncB;
        ElemA->Buffer = BufferA;
        ElemB->Buffer = BufferB;
        ElemA->BufferSize = BufferSizeA;
        ElemB->BufferSize = BufferSizeB;
        ElemA->FuncTimeDelay = FuncATimeDelay;
        ElemB->FuncTimeDelay = FuncBTimeDelay;
        Elem->ElemA = ElemA;
        Elem->ElemB = ElemB;
    }
    double t0 = Param->T0;  // revolution time of the nominal ring
    int turn = Param->nturn;  // current turn
    int i;

    /* track */
    VariableThinMPolePass(r_in, Elem, t0, turn, num_particles, Param->thread_rng);
    /* reset the polynom values.
       This is done to save the ring in an unperturbed configuration. */
    for (i = 0; i < Elem->MaxOrder + 1; i++){
        Elem->PolynomA[i] = 0;
        Elem->PolynomB[i] = 0;
    };
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
        int MaxOrder, Mode, NSamplesA, NSamplesB, KtaylorA, KtaylorB, Periodic;
        int BufferSizeA, BufferSizeB;
        double *R1, *R2, *T1, *T2, *EApertures, *RApertures;
        double *PolynomA, *PolynomB, *AmplitudeA, *AmplitudeB;
        double *Ramps, *FuncA, *FuncB, *BufferA, *BufferB;
        double FuncATimeDelay, FuncBTimeDelay;
        double FrequencyA, FrequencyB;
        double PhaseA, PhaseB;
        double SinAabove, SinBabove;
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
        SinAabove=atGetOptionalDouble(ElemData,"SinAabove", 0); check_error();
        SinBabove=atGetOptionalDouble(ElemData,"SinBabove", 0); check_error();
        Ramps=atGetOptionalDoubleArray(ElemData, "Ramps"); check_error();
        NSamplesA=atGetOptionalLong(ElemData, "NSamplesA", 0); check_error();
        NSamplesB=atGetOptionalLong(ElemData, "NSamplesB", 0); check_error();
        KtaylorA=atGetOptionalLong(ElemData, "KtaylorA", 0); check_error();
        KtaylorB=atGetOptionalLong(ElemData, "KtaylorB", 0); check_error();
        FuncA=atGetOptionalDoubleArray(ElemData,"FuncA"); check_error();
        FuncB=atGetOptionalDoubleArray(ElemData,"FuncB"); check_error();
        BufferA=atGetOptionalDoubleArray(ElemData,"BufferA"); check_error();
        BufferB=atGetOptionalDoubleArray(ElemData,"BufferB"); check_error();
        BufferSizeA=atGetOptionalLong(ElemData, "BufferSizeA", 0); check_error();
        BufferSizeB=atGetOptionalLong(ElemData, "BufferSizeB", 0); check_error();
        FuncATimeDelay=atGetOptionalDouble(ElemData,"FuncATimeDelay", 0); check_error();
        FuncBTimeDelay=atGetOptionalDouble(ElemData,"FuncBTimeDelay", 0); check_error();
        Periodic=atGetOptionalLong(ElemData,"Periodic", 0); check_error();
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
        ElemA->Sinabove = SinAabove;
        ElemB->Sinabove = SinBabove;
        ElemA->NSamples = NSamplesA;
        ElemB->NSamples = NSamplesB;
        ElemA->Ktaylor = KtaylorA;
        ElemB->Ktaylor = KtaylorB;
        ElemA->Func = FuncA;
        ElemB->Func = FuncB;
        ElemA->Buffer = BufferA;
        ElemB->Buffer = BufferB;
        ElemA->BufferSize = BufferSizeA;
        ElemB->BufferSize = BufferSizeB;
        Elem->ElemA = ElemA;
        Elem->ElemB = ElemB;
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
            mxSetCell(plhs[1], 6, mxCreateString("SinAabove"));
            mxSetCell(plhs[1], 7, mxCreateString("SinBabove"));
            mxSetCell(plhs[1], 8, mxCreateString("Ramps"));
            mxSetCell(plhs[1], 9, mxCreateString("FuncA"));
            mxSetCell(plhs[1], 10,mxCreateString("FuncB"));
            mxSetCell(plhs[1], 11,mxCreateString("BufferA"));
            mxSetCell(plhs[1], 12,mxCreateString("BufferB"));
            mxSetCell(plhs[1], 13,mxCreateString("BufferSizeA"));
            mxSetCell(plhs[1], 14,mxCreateString("BufferSizeB"));
            mxSetCell(plhs[1], 15,mxCreateString("FuncATimeDelay"));
            mxSetCell(plhs[1], 16,mxCreateString("FuncBTimeDelay"));
            mxSetCell(plhs[1], 17,mxCreateString("NSamplesA"));
            mxSetCell(plhs[1], 18,mxCreateString("NSamplesB"));
            mxSetCell(plhs[1], 19,mxCreateString("KtaylorA"));
            mxSetCell(plhs[1], 20,mxCreateString("KtaylorB"));
            mxSetCell(plhs[1], 21,mxCreateString("Periodic"));
            mxSetCell(plhs[1], 22,mxCreateString("T1"));
            mxSetCell(plhs[1], 23,mxCreateString("T2"));
            mxSetCell(plhs[1], 24,mxCreateString("R1"));
            mxSetCell(plhs[1], 25,mxCreateString("R2"));
            mxSetCell(plhs[1], 26,mxCreateString("RApertures"));
            mxSetCell(plhs[1], 27,mxCreateString("EApertures"));
        }
    } else {
        mexErrMsgIdAndTxt("AT:WrongArg", "Needs 0 or 2 arguments");
    }
}
#endif /* MATLAB_MEX_FILE */
