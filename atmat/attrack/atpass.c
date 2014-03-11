/*  Tracking engine ATPASS for Accelerator Toolbox 1.3    */

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <mex.h>

#if defined(PCWIN) || defined(PCWIN64)
#include <windows.h>
#define LIBRARYHANDLETYPE HINSTANCE
#define FREELIBFCN(libfilename) FreeLibrary((libfilename))
#define LOADLIBFCN(libfilename) LoadLibrary((libfilename))
#define GETPASSFCN(libfilename) GetProcAddress((libfilename),"passFunction")

#else
#include <dlfcn.h>
#define LIBRARYHANDLETYPE void*
#define FREELIBFCN(libfilename) dlclose(libfilename)
#define LOADLIBFCN(libfilename) dlopen((libfilename),RTLD_LAZY)
#define GETPASSFCN(libfilename) dlsym((libfilename),"passFunction")

#endif

#define LATTICE prhs[0]
#define INITCONDITIONS prhs[1]

#define NO_LOCAL_COPY 		0
#define MAKE_LOCAL_COPY 	1
#define USE_LOCAL_COPY		2

#define LIMIT_AMPLITUDE		1 /* if any of the phase space variables exceedes this limit it is marked as lost */

typedef int*(*MYPROC)(mxArray*,int*, double*, int, int);

static int **field_numbers_ptr = NULL;
static mxArray **mxPtrELEM = NULL;
static MYPROC *methods_table = NULL;
static char **method_names = NULL;
static int num_elements = 0;
static struct LibraryListElement {
    char *LibraryFileName;
    char *MethodName;
    LIBRARYHANDLETYPE LibraryHandle;
    MYPROC FunctionHandle;
    struct LibraryListElement *Next;
} *LibraryList = NULL;


struct LibraryListElement* SearchLibraryList(struct LibraryListElement *head, const char *method_name)
{
    /* recusively search the list to check if the library containing method_name is
     * already loaded. If it is - retutn the pointer to the list element. If not -
     * return NULL */

    if (head) return (strcmp(head->MethodName, method_name)==0) ?
            head : SearchLibraryList(head->Next, method_name);
    else return NULL;
}

static void cleanup(void)
{
    struct LibraryListElement *LibraryListPtr=LibraryList;
    int n;
    for (n=0; n<num_elements; n++) {
        mxFree(field_numbers_ptr[n]);
        mxFree(method_names[n]);
    }
    mxFree(field_numbers_ptr);
    mxFree(methods_table);
    mxFree(method_names);
    mxFree(mxPtrELEM);
    
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

static void checkiflost(double *DblBuffer, int num_particles, double num_elem, double num_turn,
        double *xnturn, double *xnelem, double *xcoord)
{
    int n, c;
    double *r6;
    for (c=0; c<num_particles; c++) {/* Loop over particles */
        r6 = DblBuffer+c*6;
        if (!mxIsNaN(r6[0])) {   /* No change if already marked */
           for (n=0;n<6;n++) {
                if (mxIsInf(r6[n]) || (abs(r6[n])>LIMIT_AMPLITUDE)) {
                    xnturn[c] = num_turn;
                    xnelem[c] = num_elem;
                    memcpy(xcoord+c*6,r6,6*sizeof(double));
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
    mexCallMATLAB(1,&tempmxptr,2,mxPassArg,method_name);
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
    mexCallMATLAB(1,&tempmxptr,5,mxPassArg,"feval");
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
    mexCallMATLAB(1,&tempmxptr,2,mxPassArg,"feval");
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
    mxArray *tempmxptr1, *tempmxptr2, *mxBuffer, *mxLoss;
    mxArray *mxNturn, *mxNelem, *mxCoord, *mxTurn, *mxElem;
    double *xnturn, *xnelem, *xcoord, *xturn, *xelem;
    const char *lossinfo[] = {"turn","element","coordinates"};
    mxArray *mxPassArg1[5], *mxPassArg2[2], *mxPre, *mxPost;
    
    double *DblPtrDataOut ,*DblPtrDataIn, *DblBuffer;
    
    struct LibraryListElement *LibraryListPtr;
    int *dummy;
    
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
    
    mexAtExit(cleanup);
    
    if (NewLatticeFlag || FirstCallFlag) {
        
        for(n=0;n<num_elements;n++) { /* free memory from previously used lattice */
            mxFree(field_numbers_ptr[n]);
            mxFree(method_names[n]);
        }
        
        /* Allocate persistent memory */
        num_elements  = mxGetNumberOfElements(LATTICE);
        
        /* Pointer to method names used by elements */
        method_names = (char**)mxRealloc(method_names, num_elements*sizeof(char*));
        mexMakeMemoryPersistent(method_names);
        
        /* Pointer to integer maps of Elemant data filelds used by the tracking function */
        field_numbers_ptr = (int**)mxRealloc(field_numbers_ptr, num_elements*sizeof(int*));
        mexMakeMemoryPersistent(field_numbers_ptr);
        
        /* Pointer to Element data in MATLAB */
        mxPtrELEM = (mxArray**)mxRealloc(mxPtrELEM, num_elements*sizeof(mxArray*));
        mexMakeMemoryPersistent(mxPtrELEM);
        
        /* pointer to library function used for tracking */
        methods_table = (MYPROC*)mxRealloc(methods_table, num_elements*sizeof(MYPROC));
        mexMakeMemoryPersistent(methods_table);
        
        for (n=0;n<num_elements;n++) {
            mxArray *mxElem = mxGetCell(LATTICE,n);
            
            /*      if (mxIsClass(mxElem,"at.Drift")) {
             * method_names[n] = NULL;
             * methods_table[n] = NULL;
             * }
             * else {*/
            if (!mxIsStruct(mxElem)) {
                mxArray *mxObj = mxElem;
                mxElem = mxGetProperty(mxObj,0,"AtStruct");
                /*mxDestroyArray(mxObj);*/
            }
            if (!(tempmxptr1 = mxGetField(mxElem,0,"PassMethod")))
                mexErrMsgTxt("Required field 'PassMethod' was not found in the element data structure");
            
            if (!mxIsChar(tempmxptr1)) {
                mexPrintf("Element # %d\n",n);
                mexErrMsgTxt("'PassMethod' field must be a string");
            }
            
            method_names[n] = mxArrayToString(tempmxptr1);
            mexMakeMemoryPersistent(method_names[n]);
            
            /*  Check if the library function associated with this method is already loaded */
            LibraryListPtr = SearchLibraryList(LibraryList,method_names[n]);
            
            if (LibraryListPtr)
                methods_table[n] = LibraryListPtr->FunctionHandle;
            else {
                int tempint;
                mexCallMATLAB(1,&tempmxptr2,1,&tempmxptr1,"exist");
                tempint = (int)mxGetScalar(tempmxptr2);
                mxDestroyArray(tempmxptr2);
                switch (tempint) {
                    case 2:
                        /* m-file on the search path */
                        LibraryListPtr = (struct LibraryListElement*)mxMalloc(sizeof(struct LibraryListElement));
                        mexMakeMemoryPersistent(LibraryListPtr);
                        LibraryListPtr->Next = LibraryList;
                        LibraryList = LibraryListPtr;
                        
                        LibraryListPtr->MethodName = mxArrayToString(tempmxptr1);
                        LibraryListPtr->LibraryFileName = NULL;
                        LibraryListPtr->FunctionHandle = NULL;
                        
                        methods_table[n] = NULL;
                        
                        break;
                        
                    case 3:
                        /* mex-file on the search path */
                        LibraryListPtr = (struct LibraryListElement*)mxMalloc(sizeof(struct LibraryListElement));
                        mexMakeMemoryPersistent(LibraryListPtr);
                        LibraryListPtr->Next = LibraryList;
                        LibraryList = LibraryListPtr;
                        
                        mexCallMATLAB(1,&tempmxptr2,1,&tempmxptr1,"which");
                        LibraryListPtr->MethodName = mxArrayToString(tempmxptr1);
                        LibraryListPtr->LibraryFileName=mxArrayToString(tempmxptr2);
                        mxDestroyArray(tempmxptr2);
                        
                        LibraryListPtr->LibraryHandle = LOADLIBFCN(LibraryListPtr->LibraryFileName);
                        LibraryListPtr->FunctionHandle = (MYPROC)GETPASSFCN(LibraryListPtr->LibraryHandle);
                        if(LibraryListPtr->FunctionHandle == NULL) {
                            FREELIBFCN(LibraryListPtr->LibraryHandle);
                            mexPrintf("Element # %d\tLibrary file: %s\tFunction: %s\n",n, LibraryListPtr->LibraryFileName,LibraryListPtr->MethodName);
                            mexErrMsgTxt("Library or function could not be loaded");
                        }
                        
                        methods_table[n] = LibraryListPtr->FunctionHandle;
                        
                        break;
                        
                    default:
                        mexPrintf("Element #%d\tPassMethod: '%s'\n",  n, method_names[n]);
                        mexErrMsgTxt("Specified PassMethod is not on MATLAB search path");
                }
            }
            /*	   }*/
            mxPtrELEM[n]= mxElem;
        }
    }
    
    if (nrhs >= 5) {    /* subtract 1 for C numbering: 0 to num_elements-1 */
        num_refpts = (int)mxGetNumberOfElements(prhs[4]);
        double *dblrefpts = mxGetPr(prhs[4]);
        refpts = mxCalloc(num_refpts,sizeof(int));
        for (n=0; n<num_refpts; n++) refpts[n] = ((int)dblrefpts[n])-1;
    }
    else {              /* only end of the line */
        num_refpts = 1;
        refpts = mxCalloc(num_refpts,sizeof(int));
        refpts[0] = num_elements;
    }
    
    plhs[0] = mxCreateDoubleMatrix(6,num_particles*num_refpts*num_turns,mxREAL);
    
    mxCoord=mxCreateDoubleMatrix(6,num_particles,mxREAL);
    mxNturn=mxCreateDoubleMatrix(1,num_particles,mxREAL);
    mxNelem=mxCreateDoubleMatrix(1,num_particles,mxREAL);
    mxTurn=mxCreateDoubleMatrix(1,1,mxREAL);
    mxElem=mxCreateDoubleMatrix(1,1,mxREAL);
    mxLoss=mxCreateStructMatrix(1,1,3,lossinfo);
    mxSetField(mxLoss, 0, lossinfo[0], mxNturn);
    mxSetField(mxLoss, 0, lossinfo[1], mxNelem);
    mxSetField(mxLoss, 0, lossinfo[2], mxCoord);
    xnturn=mxGetPr(mxNturn);
    xnelem=mxGetPr(mxNelem);
    xcoord=mxGetPr(mxCoord);
    xturn=mxGetPr(mxTurn);
    xelem=mxGetPr(mxElem);
    for (n=0; n<num_particles; n++) {
        double val = mxGetNaN();
        xnturn[n]=val;
        xnelem[n]=val;
        xcoord[6*n+0]=val;
        xcoord[6*n+1]=val;
        xcoord[6*n+2]=val;
        xcoord[6*n+3]=val;
        xcoord[6*n+4]=val;
        xcoord[6*n+5]=val;
    }
    
    DblPtrDataIn  = mxGetPr(INITCONDITIONS);
    DblPtrDataOut = mxGetPr(plhs[0]);
    
    mxBuffer = mxCreateDoubleMatrix(6,num_particles,mxREAL);
    DblBuffer = mxGetPr(mxBuffer);
    memcpy(DblBuffer, DblPtrDataIn, np6*sizeof(double));
    
    mxPassArg1[2] = mxBuffer;
    mxPassArg1[3] = mxTurn;
    mxPassArg1[4] = mxElem;
    mxPassArg2[1] = mxBuffer;
    mxPre  = ((nrhs >= 6) && !mxIsEmpty(prhs[5])) ? (mxArray *) prhs[5] : NULL;
    mxPost = ((nrhs >= 7) && !mxIsEmpty(prhs[6])) ? (mxArray *) prhs[6] : NULL;
    
    /* start tracking */
    nextrefindex = 0;
    nextref= (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
    *xturn = (double)1;
    if (NewLatticeFlag || FirstCallFlag) {
        /* If NewLatticeFlag || FirstCallFlag - go through all elements
         * push particles ans built persistemt data structures (field_numbers_ptr)
         * along the way
         */
        for (n=0; n<num_elements; n++) {
            *xelem = (double)(n+1);
            if (n == nextref) {
                memcpy(DblPtrDataOut, DblBuffer, np6*sizeof(double));
                DblPtrDataOut += np6; /*  shift the location to write to in the output array */
                nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
            }
            if (mxPre) {                            /* Pre-tracking optional function */
                DblBuffer=passhook(mxPassArg1,mxPtrELEM[n], mxPre);
            }
            if (methods_table[n]) {                 /*Pointer to a loaded library function */
                field_numbers_ptr[n] = (*methods_table[n])(mxPtrELEM[n],field_numbers_ptr[n],DblBuffer,num_particles,MAKE_LOCAL_COPY);
                mexMakeMemoryPersistent(field_numbers_ptr[n]);
            }
            else if (mxIsStruct(mxPtrELEM[n])) {    /* M-File */
                DblBuffer=passmfile(mxPassArg1+1, mxPtrELEM[n], method_names[n]);
            }
            else {                                  /* Element class */
                DblBuffer=passobject(mxPassArg2, mxPtrELEM[n]);
            }
            if (mxPost) {                           /* Post-tracking optional function */
                DblBuffer=passhook(mxPassArg1,mxPtrELEM[n], mxPost);
            }
            checkiflost(DblBuffer, num_particles, *xelem, *xturn, xnturn, xnelem, xcoord);
        }
    }
    else {
        for (n=0; n<num_elements; n++) {
            *xelem = (double)(n+1);
            if (n == nextref) {
                memcpy(DblPtrDataOut, DblBuffer, np6*sizeof(double));
                DblPtrDataOut += np6; /*  shift the location to write to in the output array */
                nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
            }
            if (mxPre) {                            /* Pre-tracking optional function */
                DblBuffer=passhook(mxPassArg1,mxPtrELEM[n], mxPre);
            }
            if (methods_table[n]) {                 /*Pointer to a loaded library function */
                dummy = (*methods_table[n])(mxPtrELEM[n],field_numbers_ptr[n],DblBuffer,num_particles,USE_LOCAL_COPY);
            }
            else if (mxIsStruct(mxPtrELEM[n])) {	/* M-File */
                DblBuffer=passmfile(mxPassArg1+1, mxPtrELEM[n], method_names[n]);
            }
            else {                                  /* Element class */
                DblBuffer=passobject(mxPassArg2, mxPtrELEM[n]);
            }
            if (mxPost) {                           /* Post-tracking optional function */
                DblBuffer=passhook(mxPassArg1,mxPtrELEM[n], mxPost);
            }
            checkiflost(DblBuffer, num_particles, *xelem, *xturn, xnturn, xnelem, xcoord);
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
        for (n=0;n<num_elements;n++) {
            *xelem = (double)(n+1);
            if (n == nextref) {
                memcpy(DblPtrDataOut, DblBuffer, np6*sizeof(double));
                DblPtrDataOut += np6; /*  shift the location to write to in the output array */
                nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
            }
            if (mxPre) {                            /* Pre-tracking optional function */
                DblBuffer=passhook(mxPassArg1,mxPtrELEM[n], mxPre);
            }
            if (methods_table[n]) {                 /*Pointer to a loaded library function */
                dummy = (*methods_table[n])(mxPtrELEM[n],field_numbers_ptr[n],DblBuffer,num_particles,USE_LOCAL_COPY);
            }
            else if (mxIsStruct(mxPtrELEM[n])) {    /* M-File */
                DblBuffer=passmfile(mxPassArg1+1, mxPtrELEM[n], method_names[n]);
            }
            else {                                  /* Element class */
                DblBuffer=passobject(mxPassArg2, mxPtrELEM[n]);
            }
            if (mxPost) {                           /* Post-tracking optional function */
                DblBuffer=passhook(mxPassArg1,mxPtrELEM[n], mxPost);
            }
            checkiflost(DblBuffer, num_particles, *xelem, *xturn, xnturn, xnelem, xcoord);
        }
        if (num_elements == nextref) {
            memcpy(DblPtrDataOut, DblBuffer, np6*sizeof(double));
            DblPtrDataOut += np6; /*  shift the location to write to in the output array */
            nextref = (nextrefindex<num_refpts) ? refpts[nextrefindex++] : INT_MAX;
        }
    }
    
    mxFree(refpts);
    
    if (nlhs >= 2) {
        plhs[1]=mxLoss;
    }
    else {
        mxDestroyArray(mxLoss); /* Destroys also the structure members */
        /*mxDestroyArray(mxCoord);
        mxDestroyArray(mxNturn);
        mxDestroyArray(mxNelem);*/
    }
    mxDestroyArray(mxBuffer);
    mxDestroyArray(mxTurn);
    mxDestroyArray(mxElem);
    
    FirstCallFlag=false;
}
