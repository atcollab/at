/*  Tracking engine ATPASS for Accelerator Toolbox 1.3    */
    
#include "mex.h"

#include <string.h>



#if defined(PCWIN)
    #include <windows.h>
    #define MEXEXTENSIONSTRING ".dll"
    #define PASSMETHODSUBDIRSTRING "\\simulator\\element\\"
    #define LIBRARYHANDLETYPE HINSTANCE
    #define FREELIBFCN(libfilename) FreeLibrary((libfilename))
    #define LOADLIBFCN(libfilename) LoadLibrary((libfilename))
    #define GETPASSFCN(libfilename) GetProcAddress((libfilename),"passFunction")
    
#elif defined(ALPHA)
    #include <dlfcn.h>
    #define MEXEXTENSIONSTRING ".mexaxp"
    #define PASSMETHODSUBDIRSTRING "/simulator/element/"
    #define LIBRARYHANDLETYPE void*
    #define FREELIBFCN(libfilename) dlclose(libfilename)
    #define LOADLIBFCN(libfilename) dlopen((libfilename),RTLD_LAZY)    
    #define GETPASSFCN(libfilename) dlsym((libfilename),"passFunction")
    
    
#elif defined(GLNX86)
    #include <dlfcn.h>
    #define MEXEXTENSIONSTRING ".mexglx"
    #define PASSMETHODSUBDIRSTRING "/simulator/element/"
    #define LIBRARYHANDLETYPE void*
    #define FREELIBFCN(libfilename) dlclose(libfilename)
    #define LOADLIBFCN(libfilename) dlopen((libfilename),RTLD_LAZY)    
    #define GETPASSFCN(libfilename) dlsym((libfilename),"passFunction")
        
#elif defined(SOL2)
    #include <dlfcn.h>
    #define MEXEXTENSIONSTRING ".mexsol"
    #define PASSMETHODSUBDIRSTRING "/simulator/element/"
    #define LIBRARYHANDLETYPE void*
    #define FREELIBFCN(libfilename) dlclose(libfilename)
    #define LOADLIBFCN(libfilename) dlopen((libfilename),RTLD_LAZY)    
    #define GETPASSFCN(libfilename) dlsym((libfilename),"passFunction")

/* Add other platforms here macros definition here */    
    
#endif



#define LATTICE prhs[0]
#define INITCONDITIONS prhs[1]



#define NO_LOCAL_COPY 		0
#define MAKE_LOCAL_COPY 	1
#define USE_LOCAL_COPY		2

#define LIMIT_AMPLITUDE		1 /* if any of the phase space variables exceedes this limit it is marked as lost */

int num_elements, **field_numbers_ptr;
bool FirstCallFlag = true;  

mxArray **mxPtrELEM; 

LIBRARYHANDLETYPE *loaded_dlls;

typedef int*(*MYPROC)(mxArray*,int*, double*,int, int);

MYPROC *methods_table;


char **method_names;


struct LibraryListElement
{      
    char *LibraryFileName;
    char *MethodName;
    LIBRARYHANDLETYPE LibraryHandle;
    MYPROC FunctionHandle;
    struct LibraryListElement *Next;
} *LibraryList; 


struct LibraryListElement* SearchLibraryList(struct LibraryListElement *head, char *method_name)
{   /* recusively search the list to check if the library containing method_name is
       already loaded. If it is - retutn the pointer to the list element. If not -
       return NULL 
    */
    
    struct LibraryListElement *lstptr = head;
    if(lstptr==NULL)
        return(NULL);
    else
    {   if(strcmp(lstptr->MethodName,method_name)==0)
             return(lstptr);
        else
            return(SearchLibraryList(lstptr->Next,method_name));
    }   
}

static void cleanup(void)
{   struct LibraryListElement *LibraryListPtr;
    int n;
    for(n=0;n<num_elements;n++) 
    {   mxFree(field_numbers_ptr[n]);
        mxFree(method_names[n]);
    }
    mxFree(field_numbers_ptr);
    mxFree(methods_table);
    mxFree(method_names);
    mxFree(mxPtrELEM);
    
    /* Free memory and unload libraries */
	LibraryListPtr = LibraryList;
	while(LibraryListPtr)
	    {   /* mexPrintf("Freeing library: %s\n",LibraryListPtr->LibraryFileName); */
	        FREELIBFCN(LibraryListPtr->LibraryHandle);
	        mxFree(LibraryListPtr->LibraryFileName);
	        mxFree(LibraryListPtr->MethodName);
	        LibraryList = LibraryListPtr->Next;
            
	        mxFree(LibraryListPtr);
            
	        LibraryListPtr = LibraryList;
        }
}

void markaslost(double *r6)
{	int i;
    
    r6[0] = mxGetNaN();
	
    for(i=1;i<6;i++)
		r6[i] =0;
}



void mexFunction(	int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    int m,n,p,  num_particles, np6, *refpts, num_refpts, tempint, num_turns;
    int nextrefindex=0; /* index to the array of refpts */
    bool  NewLatticeFlag;
    mxArray *tempmxptr, *tempmxptr1, *tempmxptr2, *ptr2mxptr[2], *mxBuffer;
    
    double *DblPtrDataOut ,*DblPtrDataIn, *DblBuffer, *tempdoubleptr;

    struct LibraryListElement *LibraryListPtr;
    int  *dummy;
    
    
    /* NewLatticeFlag set to true causes 
        rebuilding of the persistent data structures
        and reloading function libraries
    */
    if(mxGetScalar(prhs[2]) == 0) /* check if in 'new' or 'reuse' mode */
        NewLatticeFlag = false;
    else
        NewLatticeFlag = true;


    
num_turns = (int)mxGetScalar(prhs[3]);



num_particles = mxGetN(INITCONDITIONS);
np6 = num_particles*6;


mexAtExit(cleanup);

      
if(NewLatticeFlag || FirstCallFlag)
{  
    
    if(!FirstCallFlag) /* free memory from previously used lattice */
    {   
        for(n=0;n<num_elements;n++)
         {  mxFree(field_numbers_ptr[n]);
            mxFree(method_names[n]);
         }
         mxFree(field_numbers_ptr);
         mxFree(methods_table);
         mxFree(method_names);
         mxFree(mxPtrELEM);
    }
    
   
    
    /* Allocate persistent memory */
   
    num_elements  = mxGetNumberOfElements(LATTICE);     
    
 
    
    /* Pointer to method names used by elements */
    method_names = (char**)mxCalloc(num_elements,sizeof(char*));
    mexMakeMemoryPersistent(method_names);
    
    /* Pointer to integer maps of Elemant data filelds used by the tracking function */
    field_numbers_ptr = (int**)mxCalloc(num_elements,sizeof(int*));
    mexMakeMemoryPersistent(field_numbers_ptr);
       
    /* Pointer to Element data in MATLAB */
    mxPtrELEM = (mxArray**)mxCalloc(num_elements,sizeof(mxArray*));
    mexMakeMemoryPersistent(mxPtrELEM);
    
    /* pointer to library function used for tracking */
    methods_table = (MYPROC*)mxCalloc(num_elements, sizeof(MYPROC));
    mexMakeMemoryPersistent(methods_table);

 
    for(n=0;n<num_elements;n++)
    {   mxPtrELEM[n]= mxGetCell(LATTICE,n);
    
        tempmxptr1 = mxGetField(mxPtrELEM[n],0,"PassMethod");
        if(!tempmxptr1)
            mexErrMsgTxt("Required field 'PassMethod' was not found in the element data structure"); 
        
        if(!mxIsChar(tempmxptr1))
        {   mexPrintf("Element # %d\n",n);
            mexErrMsgTxt("'PassMethod' field must be a string"); 
        }
        
        
        method_names[n] = mxArrayToString(tempmxptr1);
                
        
        /*  Check if the library function associated with this method is already loaded */
        LibraryListPtr = SearchLibraryList(LibraryList,method_names[n]);
        
        if(LibraryListPtr)
            methods_table[n] = LibraryListPtr->FunctionHandle;

        else 
        {   mexCallMATLAB(1,&tempmxptr2,1,&tempmxptr1,"exist");

            tempint = (int)mxGetScalar(tempmxptr2);
            mxDestroyArray(tempmxptr2);
            switch (tempint)
            
            {   case 2: 
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
		            if(LibraryListPtr->FunctionHandle == NULL)
		            {   FREELIBFCN(LibraryListPtr->LibraryHandle);
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
    }
} 


if(num_turns > 1)
{   plhs[0] = mxCreateDoubleMatrix(6,num_particles*num_turns,mxREAL);
    num_refpts = 1;
    refpts = mxCalloc(num_refpts,sizeof(int));
    refpts[0] = num_elements+1;
}    
else
{   if(nrhs==5)
    {   num_refpts  =  (int)mxGetNumberOfElements(prhs[4]);
        refpts = mxCalloc(num_refpts,sizeof(int));
        tempdoubleptr = mxGetPr(prhs[4]);
        for(n=0;n<num_refpts;n++)
            refpts[n] = (int)tempdoubleptr[n];
    }
    else
    {   num_refpts = 1;
        refpts = mxCalloc(num_refpts,sizeof(int));
        refpts[0] = num_elements+1;
    }
    plhs[0] = mxCreateDoubleMatrix(6,num_particles*num_refpts,mxREAL);    
}





DblPtrDataIn  = mxGetPr(INITCONDITIONS); 
DblPtrDataOut = mxGetPr(plhs[0]);

mxBuffer = mxCreateDoubleMatrix(6,num_particles,mxREAL);
DblBuffer = mxGetPr(mxBuffer);
for(m=0;m<np6;m++)
    DblBuffer[m]=DblPtrDataIn[m];

ptr2mxptr[1] = mxBuffer;
    

if(refpts[0]==1) /* Include the starting point in the output */
{   for(m=0;m<np6;m++)
        DblPtrDataOut[m]=DblBuffer[m];
    DblPtrDataOut +=  np6;
    nextrefindex++;
}
    






/* start tracking */
if(NewLatticeFlag || FirstCallFlag)
{   /* If NewLatticeFlag || FirstCallFlag - go through all elements 
       push particles ans built persistemt data structures (field_numbers_ptr)
       along the way
    */
       
    for(n=0;n<num_elements;n++)
    {   if(methods_table[n]) /*Pointer to a loaded library function */  
        {   field_numbers_ptr[n] = (*methods_table[n])(mxPtrELEM[n],field_numbers_ptr[n],DblBuffer,num_particles,MAKE_LOCAL_COPY);   
            mexMakeMemoryPersistent(field_numbers_ptr[n]);
        }
        else /* (NULL) m-file */
        {   ptr2mxptr[0] = mxPtrELEM[n];
            mexCallMATLAB(1,&tempmxptr,2,ptr2mxptr,method_names[n]);
            /* Swap data  between two mxArrays */
            tempdoubleptr = mxGetPr(tempmxptr);
            mxSetPr(tempmxptr,mxGetPr(mxBuffer));
            mxSetPr(mxBuffer,tempdoubleptr);
            DblBuffer = tempdoubleptr;

            mxDestroyArray(tempmxptr);

        }
        if(nextrefindex<num_refpts)
            if(n==refpts[nextrefindex]-2)
            {   for(m=0;m<np6;m++)
	                DblPtrDataOut[m]=DblBuffer[m];
		        nextrefindex++;
		        DblPtrDataOut += np6; /*  shift the location to write to in the output array */
	        }

	}
}
else 
{   for(n=0;n<refpts[num_refpts-1]-1;n++)  
    {	     
        if(methods_table[n])
            dummy = (*methods_table[n])(mxPtrELEM[n],field_numbers_ptr[n],DblBuffer,num_particles,USE_LOCAL_COPY);
            
        else
        {   ptr2mxptr[0] = mxPtrELEM[n];
            mexCallMATLAB(1,&tempmxptr,2,ptr2mxptr,method_names[n]);
            /* Swap data  between two mxArrays */
            tempdoubleptr = mxGetPr(tempmxptr);
            mxSetPr(tempmxptr,mxGetPr(mxBuffer));
            mxSetPr(mxBuffer,tempdoubleptr);
            DblBuffer = tempdoubleptr;
            mxDestroyArray(tempmxptr);
        }
            
        if(n==(refpts[nextrefindex]-2))
	    {   for(m=0;m<np6;m++)
		        DblPtrDataOut[m]=DblBuffer[m];
		    nextrefindex++;
		    DblPtrDataOut += np6; /* shift the location to write to in the output array */
	    }
	}
}    


mxFree(refpts);





/* Subsequent turns */

{	for(m=1;m<num_turns;m++)
    {   for(n=0;n<num_elements;n++) 
        {	if(methods_table[n])
                dummy = (*methods_table[n])(mxPtrELEM[n],field_numbers_ptr[n],DblBuffer,num_particles,USE_LOCAL_COPY);
            else
            {   ptr2mxptr[0] = mxPtrELEM[n];
                mexCallMATLAB(1,&tempmxptr,2,ptr2mxptr,method_names[n]);
                /* Swap data  between two mxArrays */
                tempdoubleptr = mxGetPr(tempmxptr);
                mxSetPr(tempmxptr,mxGetPr(mxBuffer));
                mxSetPr(mxBuffer,tempdoubleptr);
                DblBuffer = tempdoubleptr;
                mxDestroyArray(tempmxptr);
            }
        }
        
        for(p=0;p<np6;p+=6) /* check limits */
            if(!mxIsNaN(DblBuffer[p]))
                for(n=0;n<6;n++)
                    if(!mxIsFinite(DblBuffer[p+n]) || DblBuffer[p+n]>LIMIT_AMPLITUDE || DblBuffer[p+n]<-LIMIT_AMPLITUDE)
                    {   markaslost(DblBuffer+p);
                        break;
                    }
		        


		/* copy data to the output */
		for(p=0;p<np6;p++)
		    DblPtrDataOut[p]=DblBuffer[p];
		DblPtrDataOut += np6;
	}
} 
        		            

mxDestroyArray(mxBuffer);

FirstCallFlag=false;

}
