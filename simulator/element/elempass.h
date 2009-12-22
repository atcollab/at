/* Header file for the element pass-functions */

#if defined(PCWIN)
#define ExportMode __declspec(dllexport)
#include <float.h>
#elif defined(GLNX86)
#define ExportMode
#elif defined(ALPHA)
#define ExportMode
#elif defined(SOL2)
#define ExportMode
#else 
/* Default - Windows */
#define ExportMode __declspec(dllexport)
#endif


ExportMode int* passFunction(const mxArray *ElemData, int *FieldNumbers,
				double *r_in, int num_particles, int mode);


#define NO_LOCAL_COPY 		0	/* function retieves element data from MATLAB workspace
								   each time it is called and reterns NULL pointer
								*/

#define MAKE_LOCAL_COPY 	1 	/* function retieves element data from MATLAB workspace
								   allocates memory and makes a persistent local copy
								   of the fields it uses for use on subsequent calls. 
								   Returns a pointer to that structure
                                */

#define USE_LOCAL_COPY		2  /*  Uses the previously stored local copy of the element data */

