/*  Tracking engine ATPASS for Accelerator Toolbox 1.3    */


#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <mex.h>
#include "at.h"
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


typedef int*(*MYPROC)(mxArray*,int*, double*, int, int);
typedef struct elem*(*MYPROC2)(mxArray*,struct elem*, double*, int, struct parameters*);

static int **field_numbers_ptr = NULL;
static struct elem **ElemStruct_ptr = NULL;
static MYPROC *methods_table = NULL;
static MYPROC2 *methods_table2 = NULL;

static char **method_names = NULL;
static int num_elements = 0;
static struct LibraryListElement {
    char *LibraryFileName;
    char *MethodName;
    LIBRARYHANDLETYPE LibraryHandle;
    MYPROC FunctionHandle;
    MYPROC2 FunctionHandle2;
    struct LibraryListElement *Next;
} *LibraryList = NULL;


struct LibraryListElement* SearchLibraryList(struct LibraryListElement *head, const char *method_name)
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

static void cleanup(void)
{
    struct LibraryListElement *LibraryListPtr=LibraryList;
    int n;
    for (n=0; n<num_elements; n++) {
        mxFree(field_numbers_ptr[n]);
        mxFree(ElemStruct_ptr[n]);
        mxFree(method_names[n]);
    }
    mxFree(field_numbers_ptr);
    mxFree(ElemStruct_ptr);
    mxFree(methods_table);
    mxFree(methods_table2);
    mxFree(method_names);
    
    /* Free memory and unload libraries */
    while (LibraryListPtr) {
        /* mexPrintf("Freeing library: %s\n",LibraryListPtr->LibraryFileName); */
        FREELIBFCN(LibraryListPtr->LibraryHandle);
        mxFree(LibraryListPtr->LibraryFileName);
        mxFree(LibraryListPtr->MethodName);
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

static double *passmfile(mxArray *mxPassArg[], mxArray *elem, const char *method_name)
{
    mxArray *tempmxptr;
    double *tempdoubleptr;
    mxPassArg[0] = elem;
    if (mexCallMATLAB(1,&tempmxptr,2,mxPassArg,method_name) != 0)
        mexErrMsgIdAndTxt("Atpass:PassError","error in evaluating %s",method_name);
    /* Swap data  between two mxArrays */
    tempdoubleptr = mxGetPr(tempmxptr);
    mxSetPr(tempmxptr,mxGetPr(mxPassArg[1]));
    mxSetPr(mxPassArg[1],tempdoubleptr);
    mxDestroyArray(tempmxptr);
    return tempdoubleptr;
}

static double *passhook(mxArray *mxPassArg[], mxArray *elem, mxArray *func)
{
    mxArray *tempmxptr;
    double *tempdoubleptr;
    mxPassArg[0] = func;
    mxPassArg[1] = elem;
    if (mexCallMATLAB(1,&tempmxptr,5,mxPassArg,"feval") != 0)
        mexErrMsgIdAndTxt("Atpass:HookError","error in evaluating %s","feval");
    /* Swap data  between two mxArrays */
    tempdoubleptr = mxGetPr(tempmxptr);
    mxSetPr(tempmxptr,mxGetPr(mxPassArg[2]));
    mxSetPr(mxPassArg[2],tempdoubleptr);
    mxDestroyArray(tempmxptr);
    return tempdoubleptr;
}

static double *passobject(mxArray *mxPassArg[], mxArray *elem)
{
    mxArray *tempmxptr;
    double *tempdoubleptr;
    mxPassArg[0] = mxGetProperty(elem,0,"PassFunction");
    if (mexCallMATLAB(1,&tempmxptr,2,mxPassArg,"feval") != 0)
        mexErrMsgIdAndTxt("Atpass:PassError","error in evaluating %s","feval");
    /* Swap data  between two mxArrays */
    tempdoubleptr = mxGetPr(tempmxptr);
    mxSetPr(tempmxptr,mxGetPr(mxPassArg[1]));
    mxSetPr(mxPassArg[1],tempdoubleptr);
    mxDestroyArray(tempmxptr);
    return tempdoubleptr;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int m,n,p, *refpts, num_refpts;
    int nextrefindex, nextref; /* index to the array of refpts */
    mxArray *tempmxptr1, *tempmxptr2, *mxBuffer;
    const char *lossinfo[] = {"lost", "turn", "element", "coordinates"};
    mxArray *mxLost, *mxNturn, *mxNelem, *mxCoord, *mxLoss;
    mwSize CoordDims[3] = {6,0,0};
    mxLogical *xlost;
    double *xnturn, *xnelem, *xcoord;
    mxArray *mxTurn, *mxElmn;
    mxArray *mxElem, **mxPtrELEM;
    double *xturn, *xelmn;
    mxArray *mxPassArg1[5], *mxPassArg2[2], *mxPre, *mxPost;
    mwSize outsize;
    
    double *DblPtrDataOut ,*DblPtrDataIn, *DblBuffer;
    
    struct LibraryListElement *LibraryListPtr;
    int numel;
    int *dummy;
    struct elem *dummy2;
    /* NewLatticeFlag set to true causes
     * rebuilding of the persistent data structures
     * and reloading function libraries
     */
    
    /* check if in 'new' or 'reuse' mode */
    static bool FirstCallFlag = true;
    bool NewLatticeFlag = (mxGetScalar(prhs[2]) != 0);
    int num_turns = (int)mxGetScalar(prhs[3]);
    int num_particles = mxGetN(INITCONDITIONS);
    int np6 = num_particles*6;
    int ihist, lhist;
    double *histbuf = NULL;
	double RingLength = 0;
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
    mxPtrELEM = (mxArray**)mxMalloc(numel*sizeof(mxArray*));
    for (n=0; n<numel; n++) {
       mxArray *mxLength, *mxElem = mxGetCell(LATTICE,n);
       if (!mxIsStruct(mxElem)) {
           mxArray *mxObj = mxElem;
           mxElem = mxGetProperty(mxObj,0,"AtStruct");
       }
       mxLength = mxGetField(mxElem, 0, "Length");
       if (mxLength != NULL) RingLength+=mxGetScalar(mxLength);
       mxPtrELEM[n] = mxElem;
    }

    if (NewLatticeFlag) FirstCallFlag=true;
    
    if (FirstCallFlag) {
        
        for (n=0; n<num_elements; n++) { /* free memory from previously used lattice */
            mxFree(field_numbers_ptr[n]);
            mxFree(ElemStruct_ptr[n]);
            mxFree(method_names[n]);
        }
        
        /* Allocate persistent memory */
        num_elements = numel;
        
        /* Pointer to method names used by elements */
        mxFree(method_names);       /* Use calloc to ensure uninitialized values are NULL */
        method_names = (char**)mxCalloc(num_elements,sizeof(char*));
        mexMakeMemoryPersistent(method_names);
        
        /* Pointer to integer maps of Element data fields used by the tracking function */
        mxFree(field_numbers_ptr);	/* Use calloc to ensure uninitialized values are NULL */
        field_numbers_ptr = (int**)mxCalloc(num_elements,sizeof(int*));
        mexMakeMemoryPersistent(field_numbers_ptr);
        
        /* Pointer to Element structures used by the tracking function */
        mxFree(ElemStruct_ptr);	/* Use calloc to ensure uninitialized values are NULL */
        ElemStruct_ptr = (struct elem**)mxCalloc(num_elements,sizeof(struct elem*));
        mexMakeMemoryPersistent(ElemStruct_ptr);
        
        /* pointer to library function used for tracking */
		methods_table = (MYPROC*)mxRealloc(methods_table, num_elements*sizeof(MYPROC));
		methods_table2 = (MYPROC2*)mxRealloc(methods_table2, num_elements*sizeof(MYPROC2));
		mexMakeMemoryPersistent(methods_table);
        mexMakeMemoryPersistent(methods_table2);
               
        for (n=0;n<num_elements;n++) {
            mxElem = mxPtrELEM[n];
            if (!(tempmxptr1 = mxGetField(mxElem,0,"PassMethod")))
                mexErrMsgIdAndTxt("Atpass:MissingPassMethod","Element # %d: Required field 'PassMethod' was not found in the element data structure",n);
            
            if (!mxIsChar(tempmxptr1))
                mexErrMsgIdAndTxt("Atpass:WrongPassMethod","Element # %d: 'PassMethod' field must be a string",n);
            
            method_names[n] = mxArrayToString(tempmxptr1);
            mexMakeMemoryPersistent(method_names[n]);
            
            /*  Check if the library function associated with this method is already loaded */
            LibraryListPtr = SearchLibraryList(LibraryList,method_names[n]);
            
            if (LibraryListPtr) {
                methods_table[n] = LibraryListPtr->FunctionHandle;
                methods_table2[n] = LibraryListPtr->FunctionHandle2;
            }
            else {
                int tempint;
                mexCallMATLAB(1,&tempmxptr2,1,&tempmxptr1,"exist");
                tempint = (int)mxGetScalar(tempmxptr2);
                mxDestroyArray(tempmxptr2);
                switch (tempint) {
                    case 2: /* m-file on the search path */
                        LibraryListPtr = (struct LibraryListElement*)mxMalloc(sizeof(struct LibraryListElement));
                        mexMakeMemoryPersistent(LibraryListPtr);
                        LibraryListPtr->Next = LibraryList;
                        LibraryList = LibraryListPtr;
                        
                        LibraryListPtr->MethodName = mxArrayToString(tempmxptr1);
                        LibraryListPtr->LibraryFileName = NULL;
                        LibraryListPtr->FunctionHandle = NULL;
                        LibraryListPtr->FunctionHandle2 = NULL;
                                    
                        methods_table[n] = NULL;
                        methods_table2[n] = NULL;
                        break;
                    case 3: /* mex-file not on the search path */
                        LibraryListPtr = (struct LibraryListElement*)mxMalloc(sizeof(struct LibraryListElement));
                        mexMakeMemoryPersistent(LibraryListPtr);
                        LibraryListPtr->Next = LibraryList;
                        LibraryList = LibraryListPtr;
                        
                        mexCallMATLAB(1,&tempmxptr2,1,&tempmxptr1,"which");
                        LibraryListPtr->MethodName = mxArrayToString(tempmxptr1);
                        LibraryListPtr->LibraryFileName=mxArrayToString(tempmxptr2);
                        mxDestroyArray(tempmxptr2);
                        LibraryListPtr->FunctionHandle = NULL;
                        LibraryListPtr->FunctionHandle2 = NULL;
                        LibraryListPtr->LibraryHandle = LOADLIBFCN(LibraryListPtr->LibraryFileName);
                        LibraryListPtr->FunctionHandle = (MYPROC)GETPASSFCN(LibraryListPtr->LibraryHandle);
                        LibraryListPtr->FunctionHandle2 = (MYPROC2)GETTRACKFCN(LibraryListPtr->LibraryHandle);
                                              
                        if ((LibraryListPtr->FunctionHandle == NULL) && (LibraryListPtr->FunctionHandle2 == NULL)) {
                            FREELIBFCN(LibraryListPtr->LibraryHandle);
                            mexErrMsgIdAndTxt("Atpass:UnknownLibrary",
                                    "Element # %d: Library file: %s, Function: %s: could not be loaded",
                                    n, LibraryListPtr->LibraryFileName, LibraryListPtr->MethodName);
                        }
                        
                        methods_table[n] = LibraryListPtr->FunctionHandle;
                        methods_table2[n] = LibraryListPtr->FunctionHandle2;
                                                
                        break;
                    default:
                        mexErrMsgIdAndTxt("Atpass:UnknownPassMethod",
                                "Element #%d: PassMethod '%s' is not on MATLAB search path",
                                n, method_names[n]);
                }
            }
        }
    }
	paramStruct.RingLength = RingLength;
	paramStruct.T0 = RingLength/299792458;
    if (nrhs >= 5) {    /* subtract 1 for C numbering: 0 to num_elements-1 */
        double *dblrefpts = mxGetPr(prhs[4]);
        num_refpts = (int)mxGetNumberOfElements(prhs[4]);
        refpts = mxCalloc(num_refpts,sizeof(int));
        for (n=0; n<num_refpts; n++) refpts[n] = ((int)dblrefpts[n])-1;
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
  /*CoordDims[1] = lhist;
    CoordDims[2] = num_particles;*/
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
    for (n=0; n<num_particles; n++) {
        double val = mxGetNaN();
        xnturn[n]=val;
        xnelem[n]=val;
    }
    for (n=0; n<lhist*np6; n++) {
        xcoord[n] = mxGetNaN();
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
    nextrefindex = 0;
    nextref= (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
    *xturn = (double)1;
	paramStruct.nturn = 0;
    ihist = 0;

    if (FirstCallFlag) {
        /* If FirstCallFlag go through all elements
         * push particles ans built persistemt data structures (field_numbers_ptr)
         * along the way
         */
        for (n=0; n<num_elements; n++) {
            mxElem = mxPtrELEM[n];
            *xelmn = (double)(n+1);
            if (n == nextref) {
                memcpy(DblPtrDataOut, DblBuffer, np6*sizeof(double));
                DblPtrDataOut += np6; /*  shift the location to write to in the output array */
                nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
            }
            if (lhist > 0) {
                memcpy(histbuf+ihist*np6, DblBuffer, np6*sizeof(double));
            }
            if (mxPre) {                            /* Pre-tracking optional function */
                DblBuffer=passhook(mxPassArg1, mxElem, mxPre);
            }
            if (methods_table[n]) {                 /* Pointer to a loaded library function */
                field_numbers_ptr[n] = (*methods_table[n])(mxElem,field_numbers_ptr[n],DblBuffer,num_particles,MAKE_LOCAL_COPY);
                mexMakeMemoryPersistent(field_numbers_ptr[n]);
            }
            else if (methods_table2[n]) {			/* pointer to a TrackFunction */
				ElemStruct_ptr[n] = (*methods_table2[n])(mxElem,ElemStruct_ptr[n],DblBuffer,num_particles,&paramStruct);
                mexMakeMemoryPersistent(ElemStruct_ptr[n]);
			}else if (mxIsStruct(mxElem)) {    /* M-File */
                DblBuffer=passmfile(mxPassArg1+1, mxElem, method_names[n]);
            }
            else {                                  /* Element class */
                DblBuffer=passobject(mxPassArg2, mxElem);
            }
            if (mxPost) {                           /* Post-tracking optional function */
                DblBuffer=passhook(mxPassArg1, mxElem, mxPost);
            }
            checkiflost(DblBuffer,num_particles,*xelmn,*xturn,xnturn,xnelem,xcoord,xlost,histbuf,ihist,lhist);
            if (++ihist >= lhist) ihist = 0;
        }
        FirstCallFlag=false;
    }
    else {
        for (n=0; n<num_elements; n++) {
            mxElem = mxGetCell(LATTICE,n);
            *xelmn = (double)(n+1);
            if (n == nextref) {
                memcpy(DblPtrDataOut, DblBuffer, np6*sizeof(double));
                DblPtrDataOut += np6; /*  shift the location to write to in the output array */
                nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
            }
            if (lhist > 0) {
                memcpy(histbuf+ihist*np6, DblBuffer, np6*sizeof(double));
            }
            if (mxPre) {                            /* Pre-tracking optional function */
                DblBuffer=passhook(mxPassArg1, mxElem, mxPre);
            }
            if (methods_table[n]) {                 /* Pointer to a loaded library function */
                dummy = (*methods_table[n])(mxElem,field_numbers_ptr[n],DblBuffer,num_particles,USE_LOCAL_COPY);
            }
            else if (methods_table2[n]) {			/* pointer to a trackFunction */
				dummy2 = (*methods_table2[n])(mxElem,ElemStruct_ptr[n],DblBuffer,num_particles,&paramStruct);
			}
            else if (mxIsStruct(mxElem)) {	/* M-File */
                DblBuffer=passmfile(mxPassArg1+1, mxElem, method_names[n]);
            }
            else {                                  /* Element class */
                DblBuffer=passobject(mxPassArg2, mxElem);
            }
            if (mxPost) {                           /* Post-tracking optional function */
                DblBuffer=passhook(mxPassArg1, mxElem, mxPost);
            }
            checkiflost(DblBuffer,num_particles,*xelmn,*xturn,xnturn,xnelem,xcoord,xlost,histbuf,ihist,lhist);
            if (++ihist >= lhist) ihist = 0;
        }
    }
    if (num_elements == nextref) {
        memcpy(DblPtrDataOut, DblBuffer, np6*sizeof(double));
        DblPtrDataOut += np6; /*  shift the location to write to in the output array */
        nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
    }
    
    /* Subsequent turns */
    for (m=1; m<num_turns; m++) {
        nextrefindex = 0;
        nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
        *xturn = (double)(m+1);
		paramStruct.nturn = m;
        for (n=0; n<num_elements; n++) {
            mxElem = mxPtrELEM[n];
            *xelmn = (double)(n+1);
            if (n == nextref) {
                memcpy(DblPtrDataOut, DblBuffer, np6*sizeof(double));
                DblPtrDataOut += np6; /*  shift the location to write to in the output array */
                nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
            }
            if (lhist > 0) {
                memcpy(histbuf+ihist*np6, DblBuffer, np6*sizeof(double));
            }
            if (mxPre) {                            /* Pre-tracking optional function */
                DblBuffer=passhook(mxPassArg1, mxElem, mxPre);
            }
            if (methods_table[n]) {                 /* Pointer to a loaded library function */
                dummy = (*methods_table[n])(mxElem,field_numbers_ptr[n],DblBuffer,num_particles,USE_LOCAL_COPY);
            }
            else if (methods_table2[n]) {			/* pointer to a TrackFunction */
				dummy2 = (*methods_table2[n])(mxElem,ElemStruct_ptr[n],DblBuffer,num_particles,&paramStruct);
			}
            else if (mxIsStruct(mxElem)) {    /* M-File */
                DblBuffer=passmfile(mxPassArg1+1, mxElem, method_names[n]);
            }
            else {                                  /* Element class */
                DblBuffer=passobject(mxPassArg2, mxElem);
            }
            if (mxPost) {                           /* Post-tracking optional function */
                DblBuffer=passhook(mxPassArg1, mxElem, mxPost);
            }
            checkiflost(DblBuffer,num_particles,*xelmn,*xturn,xnturn,xnelem,xcoord,xlost,histbuf,ihist,lhist);
            if (++ihist >= lhist) ihist = 0;
        }
        if (num_elements == nextref) {
            memcpy(DblPtrDataOut, DblBuffer, np6*sizeof(double));
            DblPtrDataOut += np6; /*  shift the location to write to in the output array */
            nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
        }
    }
    if (num_refpts == 0) {
        memcpy(DblPtrDataOut, DblBuffer, np6*sizeof(double));
        DblPtrDataOut += np6; /*  shift the location to write to in the output array */
    }
    
    mxFree(mxPtrELEM);
    mxFree(refpts);
    if (lhist > 0) mxFree(histbuf);
    
    if (nlhs >= 2) plhs[1]=mxLoss;
    else mxDestroyArray(mxLoss); /* Destroys also the structure members */
    mxDestroyArray(mxBuffer);
    mxDestroyArray(mxTurn);
    mxDestroyArray(mxElmn);
}
