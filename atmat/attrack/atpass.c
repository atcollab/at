/*  Tracking engine ATPASS for Accelerator Toolbox 1.3    */

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <mex.h>
#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/
#include "attypes.h"
#include "elempass.h"
#include "atrandom.c"
#include "ringproperties.c"

/* Get ready for R2018a C matrix API */
#ifndef mxGetDoubles
#define mxGetDoubles mxGetPr
#define mxSetDoubles mxSetPr
typedef double mxDouble;
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



#define LATTICE prhs[0]
#define INITCONDITIONS prhs[1]
#define NEWLATTICE prhs[2]
#define NTURNS prhs[3]
#define REFPTS prhs[4]
#define PREHOOK prhs[5]
#define POSTHOOK prhs[6]
#define LHIST prhs[7]
#define NUMTHREADS prhs[8]
#define RINGPROPERTIES prhs[9]
#define TURN prhs[10]
#define KEEPCOUNTER prhs[11]
#define SEED prhs[12]

#define LIMIT_AMPLITUDE		1	/*  if any of the phase space variables (except the sixth N.C.) 
									exceeds this limit it is marked as lost */
#define C0  	2.99792458e8


typedef int*(*pass_function)(mxArray*, int*, double*, int, int);
typedef struct elem*(*track_function)(mxArray*, struct elem*, double*, int, struct parameters*);

static int num_elements = 0;
static struct elem **elemdata_list = NULL;
static mxArray **element_list = NULL;
static double *elemlength_list = NULL;
static track_function *integrator_list = NULL;
static pass_function *pass_list = NULL;
static int **field_numbers_ptr = NULL;

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
        mexCallMATLAB(1,&mxExist,1,&mxPassMethod,"exist");
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
                mexCallMATLAB(1,&mxWhich,1,&mxPassMethod,"which");
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

static void cleanup(void)
{
    struct LibraryListElement *LibraryListPtr=LibraryList;
    int n;
    for (n=0; n<num_elements; n++) {
        mxFree(field_numbers_ptr[n]);
        mxFree(elemdata_list[n]);
    }
    mxFree(field_numbers_ptr);
    mxFree(elemdata_list);
    mxFree(element_list);
    mxFree(pass_list);
    mxFree(integrator_list);
    
    /* Free memory and unload libraries */
    while (LibraryListPtr) {
        FREELIBFCN(LibraryListPtr->LibraryHandle);
        mxFree((void *)LibraryListPtr->MethodName);
        LibraryList = LibraryListPtr->Next;
        mxFree(LibraryListPtr);
        LibraryListPtr = LibraryList;
    }
}

static void checkiflost(mxDouble *drin, int np,
        mxDouble num_elem, mxDouble num_turn, mxDouble *xnturn, mxDouble *xnelem,
        mxDouble *xcoord, mxDouble *xlostcoord, mxLogical *xlost, mxDouble *histbuf, int ihist, int lhist)
{
    int n, c;
    for (c=0; c<np; c++) {/* Loop over particles */
        if (!xlost[c]) {  /* No change if already marked */
           mxDouble *r6 = drin+c*6;
           for (n=0; n<6; n++) {	/* I remove the check on the sixth coordinate N.C. */
                if (!mxIsFinite(r6[n]) || ((fabs(r6[n])>LIMIT_AMPLITUDE)&&n<5)) {
                    int h, k=ihist;
                    xlost[c] = 1;
                    xnturn[c] = num_turn;
                    xnelem[c] = num_elem;
                    for (h=0; h<lhist; h++) {
                        if (++k >= lhist) k=0;
                        memcpy(xcoord+6*(np*h+c),histbuf+6*(np*k+c),6*sizeof(mxDouble));
                    }
                    memcpy(xlostcoord+6*c,r6,6*sizeof(mxDouble));
                    r6[0] = mxGetNaN();
                    r6[1] = 0;
                    r6[2] = 0;
                    r6[3] = 0;
                    r6[4] = 0;
                    r6[5] = 0;
                    break;
                }
            }
        }
    }
}

static mxDouble *passmfile(mxArray *mxPassArg[], mxArray *mxElem)
{
	mxArray *mxPassMethod = mxGetField(mxElem,0,"PassMethod");
    const char *method_name = mxArrayToString(mxPassMethod);
    mxArray *tempmxptr;
    mxDouble *tempdoubleptr;
    mxPassArg[0] = mxElem;
    if (mexCallMATLAB(1,&tempmxptr,2,mxPassArg,method_name) != 0)
        mexErrMsgIdAndTxt("Atpass:PassError","error in evaluating %s",method_name);
    /* Swap data  between two mxArrays */
    tempdoubleptr = mxGetDoubles(tempmxptr);
    mxSetDoubles(tempmxptr,mxGetDoubles(mxPassArg[1]));
    mxSetDoubles(mxPassArg[1],tempdoubleptr);
    mxDestroyArray(tempmxptr);
    mxFree((void *)method_name);
    return tempdoubleptr;
}

static mxDouble *passhook(mxArray *mxPassArg[], mxArray *mxElem, mxArray *func)
{
    mxArray *tempmxptr;
    mxDouble *tempdoubleptr;
    mxPassArg[0] = func;
    mxPassArg[1] = mxElem;
    if (mexCallMATLAB(1,&tempmxptr,5,mxPassArg,"feval") != 0)
        mexErrMsgIdAndTxt("Atpass:HookError","error in evaluating %s","feval");
    /* Swap data  between two mxArrays */
    tempdoubleptr = mxGetDoubles(tempmxptr);
    mxSetDoubles(tempmxptr,mxGetDoubles(mxPassArg[1]));
    mxSetDoubles(mxPassArg[1],tempdoubleptr);
    mxDestroyArray(tempmxptr);
    return tempdoubleptr;
}

/*
@param[in]      [0] LATTICE
@param[in,out]  [1] INITCONDITIONS
@param[in]      [2] NEWLATTICE
@param[in]      [3] NTURNS
@param[in]      [4] REFPTS
@param[in]      [5] PREHOOK
@param[in]      [6] POSTHOOK
@param[in]      [7] LHIST
@param[in]      [8] NUMTHREADS
@param[in]      [9] RINGPROPERTIES
@param[in]     [10] TURN
@param[in]     [11] KEEPCOUNTER
@param[in]     [12] SEED
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    static double lattice_length = 0.0;
    static int last_turn = 0;
    static int valid = 0;
    
    int turn, elem_index, npart, *refpts, num_refpts;
    int nextrefindex, nextref; /* index to the array of refpts */
    mxArray *mxBuffer;
    const char *lossinfo[] = {"lost", "turn", "element", "coordinates_at_loss", "coordinates"};
    mxArray *mxLost, *mxNturn, *mxNelem, *mxCoord, *mxLostCoord, *mxLoss;
    mwSize CoordDims[3] = {6,0,0};
    mwSize LostCoordDims[2] = {6,0};
    mxLogical *xlost;
    mxDouble *xnturn, *xnelem, *xcoord, *xlostcoord;
    mxArray *mxTurn, *mxElmn;
    mxDouble *xturn, *xelmn;
    mxArray *mxPassArg1[5], *mxPre, *mxPost;
    mwSize outsize;
    int pass_mode;

    mxDouble *drout ,*datain, *drin;
    
    struct LibraryListElement *LibraryListPtr;
    int numel  = mxGetNumberOfElements(LATTICE);

    int new_lattice = (mxGetScalar(NEWLATTICE) == 0) ? 0 : 1;
    int num_turns = (int)mxGetScalar(NTURNS);
    int num_particles = mxGetN(INITCONDITIONS);
    int counter = (nrhs >= 11) ? (int)mxGetScalar(TURN) : 0;
    int keep_counter = (nrhs >= 12) ? (int)mxGetScalar(KEEPCOUNTER) : 0;
    int seed = (nrhs >= 13) ? (int)mxGetScalar(SEED) : -1;
    int np6 = num_particles*6;
    int ihist, lhist;
    mxDouble *histbuf = NULL;
    struct parameters param;

    #ifdef _OPENMP
    int maxthreads;
    int omp_num_threads = (nrhs >= 9) ? (int)mxGetScalar(NUMTHREADS) : 0;
    if ((omp_num_threads > 0) && (num_particles > OMP_PARTICLE_THRESHOLD)) {
        int nthreads = omp_get_num_procs();
        maxthreads = omp_get_max_threads();
        if (omp_num_threads < nthreads) nthreads = omp_num_threads;
        if (num_particles < nthreads) nthreads = num_particles;
        omp_set_num_threads(nthreads);
    }
    #endif /*_OPENMP*/

    param.common_rng = &common_state;
    param.thread_rng = &thread_state;
    param.energy = 0.0;
    param.rest_energy = 0.0;
    param.charge = -1.0;
    param.nbunch = 1;
    param.num_turns = num_turns;
    param.bdiff = NULL;
    
    if (seed >= 0) {
        pcg32_srandom_r(&common_state, seed, AT_RNG_INC);
        pcg32_srandom_r(&thread_state, seed, 0);
    }
    if (keep_counter)
        param.nturn = last_turn;
    else
        param.nturn = counter;

    if (nrhs >= 10) {
        atProperties(RINGPROPERTIES, &param.energy, &param.rest_energy, &param.charge);
    }

    if (nlhs >= 2) {
        lhist = (nrhs >= 8) ? (int)mxGetScalar(LHIST) : 1;
        if (lhist < 0) {
            mexErrMsgIdAndTxt("Atpass:WrongParameter","History length must be non-negative");
        }
        else if (lhist > 0) {
            histbuf = mxCalloc(lhist*np6,sizeof(mxDouble));
        }
    }
    else {
        lhist=0;
    }
    if (!mxIsDouble(INITCONDITIONS) || mxIsComplex(INITCONDITIONS)) {
        mexErrMsgIdAndTxt("Atpass:WrongType","RIN must be real");
    }
    if (mxGetM(INITCONDITIONS) != 6) {
        mexErrMsgIdAndTxt("Atpass:WrongSize","RIN must have 6 rows");
    }
    
    mexAtExit(cleanup);

    if (new_lattice || !valid) {
        mxArray **element;
        double *elem_length;
        pass_function *oldintegrator;
        track_function *integrator;
        for (elem_index=0; elem_index<num_elements; elem_index++) { /* free memory from previously used lattice */
            mxFree(field_numbers_ptr[elem_index]);
            mxFree(elemdata_list[elem_index]);
        }        
        num_elements = numel;
        
        /* Pointer to integer maps of Element data fields used by the tracking function */
        mxFree(field_numbers_ptr);	/* Use calloc to ensure uninitialized values are NULL */
        field_numbers_ptr = (int**)mxCalloc(num_elements,sizeof(int*));
        mexMakeMemoryPersistent(field_numbers_ptr);
        
        /* Pointer to Element structures used by the tracking function */
        mxFree(elemdata_list);	/* Use calloc to ensure uninitialized values are NULL */
        elemdata_list = (struct elem**)mxCalloc(num_elements,sizeof(struct elem*));
        mexMakeMemoryPersistent(elemdata_list);

        /* Pointer to Element lengths */
        free(elemlength_list);
        elemlength_list = (double *)calloc(num_elements, sizeof(double));
        
        /* Pointer to Element list */
        element_list = (mxArray **)mxRealloc(element_list, num_elements*sizeof(mxArray *));
        mexMakeMemoryPersistent(element_list);
        
        /* pointer to the list of integrators */
		pass_list = (pass_function*)mxRealloc(pass_list, num_elements*sizeof(pass_function));
		integrator_list = (track_function*)mxRealloc(integrator_list, num_elements*sizeof(track_function));
		mexMakeMemoryPersistent(pass_list);
        mexMakeMemoryPersistent(integrator_list);
        
        lattice_length = 0.0;
        element = element_list;
        elem_length = elemlength_list;
        oldintegrator = pass_list;
        integrator = integrator_list;
        for (elem_index=0; elem_index<num_elements; elem_index++) {
            mxArray *mxElem = mxGetCell(LATTICE,elem_index);
            mxArray *mxPassMethod = mxGetField(mxElem,0,"PassMethod");
            mxArray *mxLength = mxGetField(mxElem, 0, "Length");
            double length;
            if (!mxPassMethod)
                mexErrMsgIdAndTxt("Atpass:MissingPassMethod","Element # %d: Required field 'PassMethod' was not found in the element data structure", elem_index);
            if (!mxIsChar(mxPassMethod))
                mexErrMsgIdAndTxt("Atpass:WrongPassMethod","Element # %d: 'PassMethod' field must be a string", elem_index);            
            if (mxLength) length=mxGetScalar(mxLength); else length = 0.0;
            lattice_length += length;
            LibraryListPtr = pass_method(mxPassMethod, elem_index);
            *oldintegrator++ = LibraryListPtr->PassHandle;
            *integrator++ = LibraryListPtr->TrackHandle;
            *element++ = mxElem;
            *elem_length++ = length;
        }
        pass_mode = MAKE_LOCAL_COPY;
        valid = 0;
    }
    else {
        pass_mode = USE_LOCAL_COPY;
    }

    param.RingLength = lattice_length;
    if (param.rest_energy == 0.0) {
        param.T0 = param.RingLength/C0;
    }
    else {
        double gamma0 = param.energy/param.rest_energy;
        double betagamma0 = sqrt(gamma0*gamma0 - 1.0);
        double beta0 = betagamma0/gamma0;
        param.T0 = param.RingLength/beta0/C0;
    }

    if (nrhs >= 5) {    /* subtract 1 for C numbering: 0 to num_elements-1 */
        int nref;
        mxDouble *dblrefpts = mxGetDoubles(REFPTS);
        num_refpts = (int)mxGetNumberOfElements(REFPTS);
        refpts = mxCalloc(num_refpts,sizeof(int));
        for (nref=0; nref<num_refpts; nref++) refpts[nref] = ((int)dblrefpts[nref])-1;
        if (num_refpts == 0)
            outsize=num_particles;
        else
            outsize=num_particles*num_refpts*num_turns;
    }
    else {              /* only end of the line */
        num_refpts = 1;
        refpts = mxCalloc(num_refpts,sizeof(int));
        refpts[0] = num_elements;
        outsize=num_particles*num_turns;
    }
    
    plhs[0] = mxCreateDoubleMatrix(6,outsize,mxREAL);

    mxLost=mxCreateLogicalMatrix(1,num_particles);          /* lost particle flag, initialized to 0 */
    mxNturn=mxCreateDoubleMatrix(1,num_particles,mxREAL);   /* turn number when lost */
    mxNelem=mxCreateDoubleMatrix(1,num_particles,mxREAL);   /* element number when lost */
    CoordDims[1] = num_particles;
    CoordDims[2] = lhist;
    LostCoordDims[1] = num_particles;
    mxLostCoord=mxCreateNumericArray(2,LostCoordDims,mxDOUBLE_CLASS,mxREAL);   /* Coordinates when particle is lost */
    mxCoord=mxCreateNumericArray(3,CoordDims,mxDOUBLE_CLASS,mxREAL);   /* Coordinates before particle is lost */
    xlost=mxGetLogicals(mxLost);	/* lost particle flag */
    xnturn=mxGetDoubles(mxNturn);        /* turn number when lost */
    xnelem=mxGetDoubles(mxNelem);        /* element number when lost */
    xlostcoord=mxGetDoubles(mxLostCoord);/* Coordinates when lost */
    xcoord=mxGetDoubles(mxCoord);        /* Coordinates before particle is lost */
    if (lhist>=1)
        mxLoss=mxCreateStructMatrix(1,1,5,lossinfo);
    else
        mxLoss=mxCreateStructMatrix(1,1,4,lossinfo);
    mxSetField(mxLoss, 0, lossinfo[0], mxLost);
    mxSetField(mxLoss, 0, lossinfo[1], mxNturn);
    mxSetField(mxLoss, 0, lossinfo[2], mxNelem);
    mxSetField(mxLoss, 0, lossinfo[3], mxLostCoord);
    if (lhist>=1)
        mxSetField(mxLoss, 0, lossinfo[4], mxCoord);
    mxTurn=mxCreateDoubleMatrix(1,1,mxREAL);    /* Current turn number */
    mxElmn=mxCreateDoubleMatrix(1,1,mxREAL);    /* Current element number */
    xturn=mxGetDoubles(mxTurn);                 /* Current turn number */
    xelmn=mxGetDoubles(mxElmn);                 /* Current element number */
    for (npart=0; npart<num_particles; npart++) {
        double val = mxGetNaN();
        xnturn[npart]=val;
        xnelem[npart]=val;
    }
    for (npart=0; npart<lhist*np6; npart++) {
        xcoord[npart] = mxGetNaN();
    }
    for (npart=0; npart<np6; npart++) {
        xlostcoord[npart] = mxGetNaN();
    }
        
    datain  = mxGetDoubles(INITCONDITIONS);
    drout = mxGetDoubles(plhs[0]);
    
    mxBuffer = mxCreateDoubleMatrix(6,num_particles,mxREAL);    /* Current coordinates */
    drin = mxGetDoubles(mxBuffer);              /* Current coordinates */
    memcpy(drin, datain, np6*sizeof(mxDouble));
    
    mxPassArg1[2] = mxBuffer;
    mxPassArg1[3] = mxTurn;
    mxPassArg1[4] = mxElmn;
    mxPre  = ((nrhs >= 6) && !mxIsEmpty(PREHOOK)) ? (mxArray *) PREHOOK : NULL;
    mxPost = ((nrhs >= 7) && !mxIsEmpty(POSTHOOK)) ? (mxArray *) POSTHOOK : NULL;
    
    /* start tracking */
    ihist = 0;
    for (turn=0; turn<num_turns; turn++) {
        mxArray **element = element_list;
        double *elem_length = elemlength_list;
        pass_function *oldintegrator = pass_list;
        track_function *integrator = integrator_list;
        int **field_numbers = field_numbers_ptr;
        struct elem **elemdata= elemdata_list;
        double s_coord = 0.0;

        *xturn = (mxDouble)(param.nturn+1);
        nextrefindex = 0;
        nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
        for (elem_index=0; elem_index<num_elements; elem_index++) {
            *xelmn = (mxDouble)(elem_index+1);
            param.s_coord = s_coord;
            if (elem_index == nextref) {
                memcpy(drout, drin, np6*sizeof(mxDouble));
                drout += np6; /*  shift the location to write to in the output array */
                nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
            }
            if (lhist > 0) {
                memcpy(histbuf+ihist*np6, drin, np6*sizeof(mxDouble));
            }
            if (mxPre) {                            /* Pre-tracking optional function */
                drin=passhook(mxPassArg1, *element, mxPre);
            }
            if (*integrator) {                      /* pointer to a TrackFunction */
				*elemdata = (*integrator)(*element,*elemdata,drin,num_particles,&param);
			}
            else if (*oldintegrator) {                 /* Pointer to a passFunction */
                *field_numbers = (*oldintegrator)(*element,*field_numbers,drin,num_particles,pass_mode);
            }
            else {                                  /* M-File */
                drin=passmfile(mxPassArg1+1, *element);
            }
            if (mxPost) {                           /* Post-tracking optional function */
                drin=passhook(mxPassArg1, *element, mxPost);
            }
            checkiflost(drin,num_particles,*xelmn,*xturn,xnturn,xnelem,xcoord,xlostcoord,xlost,histbuf,ihist,lhist);
            if (++ihist >= lhist) ihist = 0;
            s_coord += *elem_length++;
            element++;
            oldintegrator++;
            integrator++;
            elemdata++;
            field_numbers++;
        }
        if (num_elements == nextref) {
            memcpy(drout, drin, np6*sizeof(mxDouble));
            drout += np6; /*  shift the location to write to in the output array */
            nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
        }
        if (pass_mode == MAKE_LOCAL_COPY) {             /* First turn */
            for (elem_index=0; elem_index<num_elements; elem_index++) {
                if (field_numbers_ptr[elem_index]) mexMakeMemoryPersistent(field_numbers_ptr[elem_index]);
                if (elemdata_list[elem_index]) mexMakeMemoryPersistent(elemdata_list[elem_index]);
            }
        }   
        pass_mode = USE_LOCAL_COPY;
        param.nturn++;
    }
    valid = 1;      /* Tracking successful: the lattice can be reused */
    last_turn = param.nturn;  /* Store turn number in a static variable */

    if (num_refpts == 0) {
        memcpy(drout, drin, np6*sizeof(mxDouble));
        drout += np6; /*  shift the location to write to in the output array */
    }
    
    #ifdef _OPENMP
    if ((omp_num_threads > 0) && (num_particles > OMP_PARTICLE_THRESHOLD)) {
        omp_set_num_threads(maxthreads);
    }
    #endif /*_OPENMP*/


    mxFree(refpts);
    if (lhist > 0) mxFree(histbuf);
    
    if (nlhs >= 2) plhs[1]=mxLoss;
    else mxDestroyArray(mxLoss); /* Destroys also the structure members */
    mxDestroyArray(mxBuffer);
    mxDestroyArray(mxTurn);
    mxDestroyArray(mxElmn);
}
