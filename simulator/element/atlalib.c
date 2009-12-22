/*   File: atlalib.c
     Matrix and Vector algebra operations for Accelerator Toolbox
     A.Terebilo   8/14/00


1. Use mxCalloc, mxMalloc , mxFree for memory allocation/deallocation
2. Vector and matrix input arguments are (*double) obtained with
   mxGetPr
    Note: MATLAB internally represents matrixes as column-by-column
   1-dimentional arrays. For example  
                  A B C
                  D E F         is represented as    A D G B E H C F I
                  G H I
3. All matrixes are 6-by-6
4. All vectors  are 1-by-6 or 6-by-1
*/
		

#include "mex.h"




void ATmultmv(double *r, const double* A)
/*	multiplies 6-component column vector r by 6x6 matrix R: as in A*r 
  The result is store in the memory area of r !!!
*/

{   int i,j;
    double temp[6];

	for(i=0;i<6;i++)
     	{	temp[i]=0;
			for(j=0;j<6;j++)
        		temp[i]+=A[i+j*6]*r[j];
		}
  	for(i=0;i<6;i++)
	r[i]=temp[i];
} 

void ATmultmv55(double *r, const double* A)
/*	multiplies 5-component column vector r by 5x5 matrix R: as in A*r 
    The result is store in the memory area of r !!!
*/ 
{   int i,j;
    double temp[5];

	for(i=0;i<5;i++)
     	{	temp[i]=0;
			for(j=0;j<5;j++)
        		temp[i]+=A[i+j*5]*r[j];
		}
  	for(i=0;i<5;i++)
	r[i]=temp[i];
} 

void ATaddvv(double *r, double *dr)
{	/*	Add two 6-component vectors vectors.
		The result is store in the memory area of r !!!
    */
	int i;
	for(i=0;i<6;i++)
		r[i]+=dr[i];
}     	

void ATdrift6(double* r, double L)
/*   Input parameter L is the physical length
     1/(1+delta) normalization is done internally
*/
{	double p_norm = 1/(1+r[4]); 
	double NormL  = L*p_norm;   
	r[0]+= NormL*r[1]; 
	r[2]+= NormL*r[3];
	r[5]+= NormL*p_norm*(r[1]*r[1]+r[3]*r[3])/2;
}



void ATtranspm(double *M)
{	/* Transpose matrix M	
	   The result replaces the original matrix 
	*/
	int i,j;
	double temp;
	for(i=1;i<6;i++)
		for(j=0;j<i;j++)
			{	temp = M[i+j*6];
				M[i+j*6] = M[j+i*6];
				M[j+i*6] = temp;
			}
}


void ATmultmm(double *M2 , double *M1)
{	/* Mutrix multiplication M2*M1, the result is stored in the M1 memory area */
	int i,j,k;
	double column_temp[6];


	for (i=0;i<6;i++) 
		{	for (j=0;j<6;j++)
				{	column_temp[j] = 0;
					for (k=0;k<6;k++)
						column_temp[j] += M2[j+6*k]*M1[k+6*i];
				}
			for (j=0;j<6;j++)
					M1[j+6*i] = 	column_temp[j];
		}
}


void ATmultmm55(double *M2 , double *M1)
{	/* Mutrix multiplication M2*M1, the result is stored in the M1 memory area */
	int i,j,k;
	double column_temp[5];


	for (i=0;i<5;i++) 
		{	for (j=0;j<5;j++)
				{	column_temp[j] = 0;
					for (k=0;k<5;k++)
						column_temp[j] += M2[j+5*k]*M1[k+5*i];
				}
			for (j=0;j<5;j++)
					M1[j+5*i] = 	column_temp[j];
		}
}



void ATsandwichmmt(double *M ,double *B)
/* calculates the matrix product M*B*M' (M' = M transposed)
   The result is stored in B memory area
*/
{	int i,j,k;
	double row_temp[6];

	ATmultmm(M,B);
	
	for (i=0;i<6;i++) 
		{	for (j=0;j<6;j++)
				{	row_temp[j] = 0;
					for (k=0;k<6;k++)
						row_temp[j] += B[i+6*k]*M[j+6*k];
				}
			for (j=0;j<6;j++)
					B[i+6*j] = 	row_temp[j];
		}
}

void ATaddmm(double *M2 , double *M1)
/* adds two 6-by-6  matrixes  M1, M2 element-by-element
   The result is stored in M1 memory area
*/
{	int i,j;
	for (i=0;i<6;i++) 
		for (j=0;j<6;j++)
			M1[i+6*j] += M2[i+6*j];

}