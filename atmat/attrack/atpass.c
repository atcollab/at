/*  Tracking engine ATPASS for Accelerator Toolbox 1.3    */

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <mex.h>
#include "attypes.h"
#include "elempass.h"

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

#define NO_LOCAL_COPY 		0
#define MAKE_LOCAL_COPY 	1
#define USE_LOCAL_COPY		2

#define LIMIT_AMPLITUDE		1	/*  if any of the phase space variables (except the sixth N.C.) 
									exceedes this limit it is marked as lost */


typedef int*(*MYPROC)(mxArray*, int*, double*, int, int);
typedef struct elem*(*MYPROC2)(mxArray*,struct elem*, double*, int, struct parameters*);

static int num_elements = 0;
static int **field_numbers_ptr = NULL;
static struct elem **elemdata_list = NULL;
static mxArray **element_list = NULL;
static MYPROC *methods_table = NULL;
static MYPROC2 *methods_table2 = NULL;

static struct LibraryListElement {
    const char *MethodName;
    LIBRARYHANDLETYPE LibraryHandle;
    MYPROC FunctionHandle;
    MYPROC2 FunctionHandle2;
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
                LibraryListPtr->FunctionHandle = NULL;
                LibraryListPtr->FunctionHandle2 = NULL;
                break;
            case 3: /* mex-file not on the search path */
                mexCallMATLAB(1,&mxWhich,1,&mxPassMethod,"which");
                LibraryFileName=mxArrayToString(mxWhich);
                mxDestroyArray(mxWhich);

                LibraryListPtr->LibraryHandle = LOADLIBFCN(LibraryFileName);
                LibraryListPtr->FunctionHandle = (MYPROC)GETPASSFCN(LibraryListPtr->LibraryHandle);
                LibraryListPtr->FunctionHandle2 = (MYPROC2)GETTRACKFCN(LibraryListPtr->LibraryHandle);
                if ((LibraryListPtr->FunctionHandle == NULL) && (LibraryListPtr->FunctionHandle2 == NULL)) {
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
    mxFree(methods_table);
    mxFree(methods_table2);
    
    /* Free memory and unload libraries */
    while (LibraryListPtr) {
        FREELIBFCN(LibraryListPtr->LibraryHandle);
        mxFree((void *)LibraryListPtr->MethodName);
        LibraryList = LibraryListPtr->Next;
        mxFree(LibraryListPtr);
        LibraryListPtr = LibraryList;
    }
}

static void checkiflost(double *DblBuffer, int np,
        double num_elem, double num_turn, double *xnturn, double *xnelem,
        double *xcoord, mxLogical *xlost, double *histbuf, int ihist, int lhist)
{
    int n, c;
    for (c=0; c<np; c++) {/* Loop over particles */
        if (!xlost[c]) {  /* No change if already marked */
           double *r6 = DblBuffer+c*6;
           for (n=0; n<5; n++) {	/* I remove the check on the sixth coordinate N.C. */
                if (!mxIsFinite(r6[n]) || (fabs(r6[n])>LIMIT_AMPLITUDE)) {
                    int h, k=ihist;
                    xlost[c] = 1;
                    xnturn[c] = num_turn;
                    xnelem[c] = num_elem;
                    for (h=0; h<lhist; h++) {
                        if (++k >= lhist) k=0;
                        memcpy(xcoord+6*(np*h+c),histbuf+6*(np*k+c),6*sizeof(double));
                    }
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

static double *passmfile(mxArray *mxPassArg[], mxArray *mxElem)
{
	mxArray *mxPassMethod = mxGetField(mxElem,0,"PassMethod");
    const char *method_name = mxArrayToString(mxPassMethod);
    mxArray *tempmxptr;
    double *tempdoubleptr;
    mxPassArg[0] = mxElem;
    if (mexCallMATLAB(1,&tempmxptr,2,mxPassArg,method_name) != 0)
        mexErrMsgIdAndTxt("Atpass:PassError","error in evaluating %s",method_name);
    /* Swap data  between two mxArrays */
    tempdoubleptr = mxGetPr(tempmxptr);
    mxSetPr(tempmxptr,mxGetPr(mxPassArg[1]));
    mxSetPr(mxPassArg[1],tempdoubleptr);
    mxDestroyArray(tempmxptr);
    mxFree((void *)method_name);
    return tempdoubleptr;
}

static double *passhook(mxArray *mxPassArg[], mxArray *mxElem, mxArray *func)
{
    mxArray *tempmxptr;
    double *tempdoubleptr;
    mxPassArg[0] = func;
    mxPassArg[1] = mxElem;
    if (mexCallMATLAB(1,&tempmxptr,5,mxPassArg,"feval") != 0)
        mexErrMsgIdAndTxt("Atpass:HookError","error in evaluating %s","feval");
    /* Swap data  between two mxArrays */
    tempdoubleptr = mxGetPr(tempmxptr);
    mxSetPr(tempmxptr,mxGetPr(mxPassArg[2]));
    mxSetPr(mxPassArg[2],tempdoubleptr);
    mxDestroyArray(tempmxptr);
    return tempdoubleptr;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    static double lattice_length= 0.0;
    
    int turn, nelem, npart, *refpts, num_refpts;
    int nextrefindex, nextref; /* index to the array of refpts */
    mxArray *mxBuffer;
    const char *lossinfo[] = {"lost", "turn", "element", "coordinates"};
    mxArray *mxLost, *mxNturn, *mxNelem, *mxCoord, *mxLoss;
    mwSize CoordDims[3] = {6,0,0};
    mxLogical *xlost;
    double *xnturn, *xnelem, *xcoord;
    mxArray *mxTurn, *mxElmn;
    double *xturn, *xelmn;
    mxArray *mxPassArg1[5], *mxPassArg2[2], *mxPre, *mxPost;
    mwSize outsize;
    int pass_mode;
    
    double *DblPtrDataOut ,*DblPtrDataIn, *DblBuffer;
    
    struct LibraryListElement *LibraryListPtr;
    int numel;
    int *dummy;
    struct elem *dummy2;
    /* reuse set to false causes
     * rebuilding of the persistent data structures
     * and reloading function libraries
     */
    
    /* check if in 'new' or 'reuse' mode */
    static bool new_lattice = true;
    bool reuse = (mxGetScalar(prhs[2]) == 0);
    int num_turns = (int)mxGetScalar(prhs[3]);
    int num_particles = mxGetN(INITCONDITIONS);
    int np6 = num_particles*6;
    int ihist, lhist;
    double *histbuf = NULL;
    struct parameters paramStruct;
	
    if (nlhs >= 2) {
        lhist = (nrhs >= 8) ? (int)mxGetScalar(prhs[7]) : 1;
        if (lhist < 0) {
            mexErrMsgIdAndTxt("Atpass:WrongParameter","History length must be non-negative");
        }
        else if (lhist > 0) {
            histbuf = mxCalloc(lhist*np6,sizeof(double));
        }
    }
    else {
        lhist=0;
    }

    mexAtExit(cleanup);

    /* Pointer to Element data in MATLAB */
    numel  = mxGetNumberOfElements(LATTICE);

    if (!reuse) new_lattice=true;
    if (new_lattice) {
        mxArray **element;
        MYPROC *integrate1;
        MYPROC2 *integrate2;
        for (nelem=0; nelem<num_elements; nelem++) { /* free memory from previously used lattice */
            mxFree(field_numbers_ptr[nelem]);
            mxFree(elemdata_list[nelem]);
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
        
        /* Pointer to Element list */
        element_list = (mxArray **)mxRealloc(element_list, num_elements*sizeof(mxArray *));
        mexMakeMemoryPersistent(element_list);
        
        /* pointer to the list of integrators */
		methods_table = (MYPROC*)mxRealloc(methods_table, num_elements*sizeof(MYPROC));
		methods_table2 = (MYPROC2*)mxRealloc(methods_table2, num_elements*sizeof(MYPROC2));
		mexMakeMemoryPersistent(methods_table);
        mexMakeMemoryPersistent(methods_table2);
        
        lattice_length = 0.0;
        element = element_list;
        integrate1 = methods_table;
        integrate2 = methods_table2;
        for (nelem=0; nelem<num_elements; nelem++) {
            mxArray *mxElem = mxGetCell(LATTICE,nelem);
            mxArray *mxPassMethod = mxGetField(mxElem,0,"PassMethod");
            mxArray *mxLength = mxGetField(mxElem, 0, "Length");
            if (!mxPassMethod)
                mexErrMsgIdAndTxt("Atpass:MissingPassMethod","Element # %d: Required field 'PassMethod' was not found in the element data structure", nelem);
            if (!mxIsChar(mxPassMethod))
                mexErrMsgIdAndTxt("Atpass:WrongPassMethod","Element # %d: 'PassMethod' field must be a string", nelem);            
            if (mxLength) lattice_length+=mxGetScalar(mxLength);
            LibraryListPtr = pass_method(mxPassMethod, nelem);
            *integrate1++ = LibraryListPtr->FunctionHandle;
            *integrate2++ = LibraryListPtr->FunctionHandle2;
            *element++=mxElem;
        }
        pass_mode = MAKE_LOCAL_COPY;
        new_lattice = false;
    }
    else {
        pass_mode = USE_LOCAL_COPY;
    }
	paramStruct.RingLength = lattice_length;
	paramStruct.T0 = lattice_length/299792458;
    if (nrhs >= 5) {    /* subtract 1 for C numbering: 0 to num_elements-1 */
        int nref;
        double *dblrefpts = mxGetPr(prhs[4]);
        num_refpts = (int)mxGetNumberOfElements(prhs[4]);
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
    mxCoord=mxCreateNumericArray(3,CoordDims,mxDOUBLE_CLASS,mxREAL);   /* Coordinates when lost */
    xlost=mxGetLogicals(mxLost);	/* lost particle flag */
    xnturn=mxGetPr(mxNturn);        /* turn number when lost */
    xnelem=mxGetPr(mxNelem);        /* element number when lost */
    xcoord=mxGetPr(mxCoord);        /* Coordinates when lost */
    mxLoss=mxCreateStructMatrix(1,1,4,lossinfo);
    mxSetField(mxLoss, 0, lossinfo[0], mxLost);
    mxSetField(mxLoss, 0, lossinfo[1], mxNturn);
    mxSetField(mxLoss, 0, lossinfo[2], mxNelem);
    mxSetField(mxLoss, 0, lossinfo[3], mxCoord);
    mxTurn=mxCreateDoubleMatrix(1,1,mxREAL);    /* Current turn number */
    mxElmn=mxCreateDoubleMatrix(1,1,mxREAL);    /* Current element number */
    xturn=mxGetPr(mxTurn);                      /* Current turn number */
    xelmn=mxGetPr(mxElmn);                      /* Current element number */
    for (npart=0; npart<num_particles; npart++) {
        double val = mxGetNaN();
        xnturn[npart]=val;
        xnelem[npart]=val;
    }
    for (npart=0; npart<lhist*np6; npart++) {
        xcoord[npart] = mxGetNaN();
    }
    
    DblPtrDataIn  = mxGetPr(INITCONDITIONS);
    DblPtrDataOut = mxGetPr(plhs[0]);
    
    mxBuffer = mxCreateDoubleMatrix(6,num_particles,mxREAL);    /* Current coordinates */
    DblBuffer = mxGetPr(mxBuffer);              /* Current coordinates */
    memcpy(DblBuffer, DblPtrDataIn, np6*sizeof(double));
    
    mxPassArg1[2] = mxBuffer;
    mxPassArg1[3] = mxTurn;
    mxPassArg1[4] = mxElmn;
    mxPassArg2[1] = mxBuffer;
    mxPre  = ((nrhs >= 6) && !mxIsEmpty(prhs[5])) ? (mxArray *) prhs[5] : NULL;
    mxPost = ((nrhs >= 7) && !mxIsEmpty(prhs[6])) ? (mxArray *) prhs[6] : NULL;
    
    /* start tracking */
    ihist = 0;
    for (turn=0; turn<num_turns; turn++) {
        mxArray **element = element_list;
        MYPROC *integrate1 = methods_table;
        MYPROC2 *integrate2 = methods_table2;
        int **field_numbers = field_numbers_ptr;
        struct elem **elemdata= elemdata_list;

        nextrefindex = 0;
        nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
        *xturn = (double)(turn+1);
		paramStruct.nturn = turn;
        for (nelem=0; nelem<num_elements; nelem++) {
            *xelmn = (double)(nelem+1);
            if (nelem == nextref) {
                memcpy(DblPtrDataOut, DblBuffer, np6*sizeof(double));
                DblPtrDataOut += np6; /*  shift the location to write to in the output array */
                nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
            }
            if (lhist > 0) {
                memcpy(histbuf+ihist*np6, DblBuffer, np6*sizeof(double));
            }
            if (mxPre) {                            /* Pre-tracking optional function */
                DblBuffer=passhook(mxPassArg1, *element, mxPre);
            }
            if (*integrate2) {                      /* pointer to a TrackFunction */
				*elemdata = (*integrate2)(*element,*elemdata,DblBuffer,num_particles,&paramStruct);
			}
            else if (*integrate1) {                 /* Pointer to a passFunction */
                *field_numbers = (*integrate1)(*element,*field_numbers,DblBuffer,num_particles,pass_mode);
            }
            else {                                  /* M-File */
                DblBuffer=passmfile(mxPassArg1+1, *element);
            }
            if (mxPost) {                           /* Post-tracking optional function */
                DblBuffer=passhook(mxPassArg1, *element, mxPost);
            }
            checkiflost(DblBuffer,num_particles,*xelmn,*xturn,xnturn,xnelem,xcoord,xlost,histbuf,ihist,lhist);
            if (++ihist >= lhist) ihist = 0;
            integrate1++;
            integrate2++;
            elemdata++;
            field_numbers++;
            element++;
        }
        if (num_elements == nextref) {
            memcpy(DblPtrDataOut, DblBuffer, np6*sizeof(double));
            DblPtrDataOut += np6; /*  shift the location to write to in the output array */
            nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
        }
        if (pass_mode == MAKE_LOCAL_COPY) {             /* First turn */
            for (nelem=0; nelem<num_elements; nelem++) {
                if (field_numbers_ptr[nelem]) mexMakeMemoryPersistent(field_numbers_ptr[nelem]);
                if (elemdata_list[nelem]) mexMakeMemoryPersistent(elemdata_list[nelem]);
            }
        }   
        pass_mode = USE_LOCAL_COPY;
    }
    if (num_refpts == 0) {
        memcpy(DblPtrDataOut, DblBuffer, np6*sizeof(double));
        DblPtrDataOut += np6; /*  shift the location to write to in the output array */
    }
    
    mxFree(refpts);
    if (lhist > 0) mxFree(histbuf);
    
    if (nlhs >= 2) plhs[1]=mxLoss;
    else mxDestroyArray(mxLoss); /* Destroys also the structure members */
    mxDestroyArray(mxBuffer);
    mxDestroyArray(mxTurn);
    mxDestroyArray(mxElmn);
}
