/*  Tracking engine ATPASS for Accelerator Toolbox 1.3    */

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <mex.h>
#include "attypes.h"
#include "atrandom.c"

/* Get ready for R2018a C matrix API */
#ifndef mxGetDoubles
#define mxGetDoubles mxGetPr
#define mxSetDoubles mxSetPr
typedef double mwDouble;
#endif

#if defined(PCWIN) || defined(PCWIN64)
#include <windows.h>
#define LIBRARYHANDLETYPE HINSTANCE
#define FREELIBFCN(libfilename) FreeLibrary((libfilename))
#define LOADLIBFCN(libfilename) LoadLibrary((libfilename))
#define GETPASSFCN(libfilename) GetProcAddress((libfilename),"passFunction")
#define GETTRACKFCN(libfilename) GetProcAddress((libfilename),"trackFunction")

#else
#include <dlfcn.h>
#define LIBRARYHANDLETYPE void*
#define FREELIBFCN(libfilename) dlclose(libfilename)
#define LOADLIBFCN(libfilename) dlopen((libfilename),RTLD_LAZY)
#define GETPASSFCN(libfilename) dlsym((libfilename),"passFunction")
#define GETTRACKFCN(libfilename) dlsym((libfilename),"trackFunction")

#endif

#define C0  	2.99792458e8

typedef int*(*pass_function)(mxArray*, int*, double*, int, int);
typedef struct elem*(*track_function)(const mxArray*, struct elem*, double*, int, struct parameters*);

/* state buffers for RNGs */
static pcg32_random_t common_state = COMMON_PCG32_INITIALIZER;
static pcg32_random_t thread_state = THREAD_PCG32_INITIALIZER;

static struct LibraryListElement {
    const char *MethodName;
    LIBRARYHANDLETYPE LibraryHandle;
    pass_function PassHandle;
    track_function TrackHandle;
    struct LibraryListElement *Next;
} *LibraryList = NULL;


static struct LibraryListElement* SearchLibraryList(struct LibraryListElement *head, const char *method_name)
{
    /* recusively search the list to check if the library containing method_name is
     * already loaded. If it is - retutn the pointer to the list element. If not -
     * return NULL */
    if (head)
        return (strcmp(head->MethodName, method_name)==0) ? head :
            SearchLibraryList(head->Next, method_name);
    else
        return NULL;
}

static struct LibraryListElement* pass_method(mxArray *mxPassMethod, int nelem)
{
    const char *fn_name = mxArrayToString(mxPassMethod);
    struct LibraryListElement *LibraryListPtr = SearchLibraryList(LibraryList, fn_name);

    if (!LibraryListPtr) {
        int tempint;
        mxArray *mxExist, *mxWhich;
        const char *LibraryFileName;
        mexCallMATLAB(1, &mxExist, 1, &mxPassMethod, "exist");
        tempint = (int)mxGetScalar(mxExist);
        mxDestroyArray(mxExist);
        LibraryListPtr = (struct LibraryListElement*)mxMalloc(sizeof(struct LibraryListElement));
        LibraryListPtr->MethodName = fn_name;
        switch (tempint) {
            case 2: /* m-file on the search path */
                LibraryListPtr->LibraryHandle = NULL;
                LibraryListPtr->PassHandle = NULL;
                LibraryListPtr->TrackHandle = NULL;
                break;
            case 3: /* mex-file not on the search path */
                mexCallMATLAB(1, &mxWhich, 1, &mxPassMethod, "which");
                LibraryFileName=mxArrayToString(mxWhich);
                mxDestroyArray(mxWhich);

                LibraryListPtr->LibraryHandle = LOADLIBFCN(LibraryFileName);
                LibraryListPtr->PassHandle = (pass_function)GETPASSFCN(LibraryListPtr->LibraryHandle);
                LibraryListPtr->TrackHandle = (track_function)GETTRACKFCN(LibraryListPtr->LibraryHandle);
                if ((LibraryListPtr->PassHandle == NULL) && (LibraryListPtr->TrackHandle == NULL)) {
                    FREELIBFCN(LibraryListPtr->LibraryHandle);
                    mexErrMsgIdAndTxt("Atpass:UnknownLibrary",
                            "Element # %d: Library file: %s, No passFunction or trackFunction available",
                            nelem, LibraryFileName);
                }
                mxFree((void *)LibraryFileName);
                break;
            default:
                mexErrMsgIdAndTxt("Atpass:UnknownPassMethod",
                        "Element #%d: PassMethod '%s' is not on MATLAB search path",
                        nelem, LibraryListPtr->MethodName);
        }
        mexMakeMemoryPersistent(LibraryListPtr);
        mexMakeMemoryPersistent((void *)LibraryListPtr->MethodName);
        LibraryListPtr->Next = LibraryList;
        LibraryList = LibraryListPtr;
    }
    else {
        mxFree((void *)fn_name);
    }
    return LibraryListPtr;
}

static double getdouble(const mxArray *input, int index)
{
    if (!(mxIsDouble(input) && mxIsScalar(input)))
        mexErrMsgIdAndTxt("Atpass:WrongParameter","Argument %d must be a scalar double", index);
    return mxGetScalar(input);
}

static double getoptionaldoubleprop(const mxArray *obj, const char *fieldname, double default_value)
{
    const mxArray *field=mxGetProperty(obj, 0, fieldname);
    return (field) ? mxGetScalar(field) : default_value;
}

static void getparticle(const mxArray *part, double *rest_energy, double *charge, int index)
{
    if (mxIsClass(part, "atparticle")) {  /* OK */
        *rest_energy = getoptionaldoubleprop(part, "rest_energy", 0.0);
        *charge = getoptionaldoubleprop(part, "charge", -1.0);
    }
    else {                              /* particle is not a Particle object */
        mexErrMsgIdAndTxt("Atpass:WrongParameter","Argument %d must be an 'atparticle' object", index);
    }
}

/* Input arguments */

#define ELEMENT 0
#define RIN 1
#define ENERGY 2
#define PARTICLE 3
#define CELL_LENGTH 4
#define S_POS 5

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    struct parameters param;
    double orbit[6];

    if (nrhs != 6)  mexErrMsgIdAndTxt("Atpass:ArgNumber", "Needs 6 arguments");
    const mxArray *mxElem = prhs[ELEMENT];
    const mxArray *OrbitData = prhs[1];
    if (mxGetM(OrbitData) * mxGetN(OrbitData) != 6)
        mexErrMsgIdAndTxt("Atpass:WrongParameter", "Argument %n must be a 6 x 1 matrix", RIN);
    const double *drin = mxGetDoubles(OrbitData);
    const mxArray *mxPassMethod = mxGetField(mxElem, 0, "PassMethod");
    if (!mxPassMethod)
        mexErrMsgIdAndTxt("Atpass:MissingPassMethod","Required field 'PassMethod' was not found in the element data structure");
    if (!mxIsChar(mxPassMethod))
        mexErrMsgIdAndTxt("Atpass:WrongPassMethod","'PassMethod' field must be a string");
    struct LibraryListElement *LibraryListPtr = pass_method((mxArray *)mxPassMethod, 0);
    track_function integrator = LibraryListPtr->TrackHandle;
    plhs[0] = mxCreateDoubleMatrix(6, 6, mxREAL);
    double *bdiff = mxGetDoubles(plhs[0]);

    param.common_rng = &common_state;
    param.thread_rng = &thread_state;
    param.energy = getdouble(prhs[ENERGY], ENERGY);
    param.rest_energy = 0.0;
    param.charge = -1.0;
    param.nbunch = 3;
    param.num_turns = 0;
    param.bdiff = bdiff;
    param.nturn = 0;
    getparticle(prhs[PARTICLE], &param.rest_energy, &param.charge, PARTICLE);
    param.RingLength = getdouble(prhs[CELL_LENGTH], CELL_LENGTH);
    if (param.rest_energy == 0.0) {
        param.T0 = param.RingLength/C0;
    }
    else {
        double gamma0 = param.energy/param.rest_energy;
        double betagamma0 = sqrt(gamma0*gamma0 - 1.0);
        double beta0 = betagamma0/gamma0;
        param.T0 = param.RingLength/beta0/C0;
    }
    param.s_coord = getdouble(prhs[S_POS], S_POS);

    for (int i = 0; i < 36; i++)
        bdiff[i] = 0.0;
    for (int i = 0; i < 6; i++)
        orbit[i] = drin[i];

    struct elem *elemdata = (*integrator)(mxElem, NULL, orbit, 1, &param);
    mxFree(elemdata);
}
